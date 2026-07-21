use std::collections::HashSet;
use std::error::Error;
use std::io::{self, BufRead, Write};

use memchr::memchr;

use crate::io_util::{open_buf_reader, open_buf_writer};

/// Counters collected while streaming sequence files.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct ProcessStats {
    pub total_sequences: usize,
    pub written_sequences: usize,
}

/// Process one or more sequence files and extract matching read IDs.
///
/// Multiple inputs are streamed sequentially into a single output to avoid
/// materializing each file in memory. Gzip inputs/outputs are supported when
/// paths end with `.gz` or the input starts with gzip magic bytes.
pub fn process_sequence_files(
    input_files: &[String],
    save_readids: &HashSet<String>,
    output_file: &str,
    exclude: bool,
) -> Result<ProcessStats, Box<dyn Error + Send + Sync>> {
    let mut writer = open_buf_writer(output_file)?;
    let mut stats = ProcessStats::default();

    for input_file in input_files {
        let file_stats = process_sequence_file(input_file, save_readids, exclude, &mut writer)?;
        stats.total_sequences += file_stats.total_sequences;
        stats.written_sequences += file_stats.written_sequences;
    }

    writer.flush()?;
    Ok(stats)
}

/// Extract paired-end reads, writing mates to separate outputs.
pub fn process_paired_sequence_files(
    input_r1: &str,
    input_r2: &str,
    save_readids: &HashSet<String>,
    output_r1: &str,
    output_r2: &str,
    exclude: bool,
) -> Result<ProcessStats, Box<dyn Error + Send + Sync>> {
    let mut writer_r1 = open_buf_writer(output_r1)?;
    let mut writer_r2 = open_buf_writer(output_r2)?;
    let mut reader_r1 = open_buf_reader(input_r1)?;
    let mut reader_r2 = open_buf_reader(input_r2)?;

    let mut pending_r1: Option<Vec<u8>> = None;
    let mut pending_r2: Option<Vec<u8>> = None;
    let mut stats = ProcessStats::default();

    loop {
        let Some(record_r1) = read_next_record(&mut reader_r1, &mut pending_r1, input_r1)? else {
            break;
        };
        let Some(record_r2) = read_next_record(&mut reader_r2, &mut pending_r2, input_r2)? else {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("paired file '{input_r2}' ended before '{input_r1}'"),
            )
            .into());
        };

        stats.total_sequences += 1;
        let id = parse_id(&record_r1.header);
        let should_write = save_readids.contains(id) != exclude;
        if should_write {
            writer_r1.write_all(&record_r1.bytes)?;
            writer_r2.write_all(&record_r2.bytes)?;
            stats.written_sequences += 1;
        }
    }

    if read_next_record(&mut reader_r2, &mut pending_r2, input_r2)?.is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("paired file '{input_r2}' has more records than '{input_r1}'"),
        )
        .into());
    }

    writer_r1.flush()?;
    writer_r2.flush()?;
    Ok(stats)
}

struct SequenceRecord {
    header: Vec<u8>,
    bytes: Vec<u8>,
}

fn process_sequence_file<W: Write>(
    input_file: &str,
    save_readids: &HashSet<String>,
    exclude: bool,
    writer: &mut W,
) -> io::Result<ProcessStats> {
    let mut reader = open_buf_reader(input_file)?;
    let mut pending_header: Option<Vec<u8>> = None;
    let mut stats = ProcessStats::default();

    while let Some(record) = read_next_record(&mut reader, &mut pending_header, input_file)? {
        stats.total_sequences += 1;
        let id = parse_id(&record.header);
        let should_write = save_readids.contains(id) != exclude;
        if should_write {
            writer.write_all(&record.bytes)?;
            stats.written_sequences += 1;
        }
    }

    Ok(stats)
}

fn read_next_record(
    reader: &mut dyn BufRead,
    pending_header: &mut Option<Vec<u8>>,
    input_file: &str,
) -> io::Result<Option<SequenceRecord>> {
    let header = if let Some(header) = pending_header.take() {
        header
    } else {
        let mut line = Vec::new();
        loop {
            line.clear();
            if reader.read_until(b'\n', &mut line)? == 0 {
                return Ok(None);
            }
            if !line.iter().all(u8::is_ascii_whitespace) {
                break;
            }
        }
        line
    };

    let marker = header
        .first()
        .copied()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "empty sequence header"))?;
    if marker != b'>' && marker != b'@' {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("expected FASTA/FASTQ header in '{input_file}'"),
        ));
    }

    let id = parse_id(&header).to_string();
    let mut record = header.clone();

    if marker == b'@' {
        for line_number in 2..=4 {
            let mut line = Vec::new();
            if reader.read_until(b'\n', &mut line)? == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    format!("truncated FASTQ record '{id}' at line {line_number}"),
                ));
            }
            record.extend_from_slice(&line);
        }
    } else {
        loop {
            let mut line = Vec::new();
            if reader.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            if line.first() == Some(&b'>') {
                *pending_header = Some(line);
                break;
            }
            record.extend_from_slice(&line);
        }
    }

    Ok(Some(SequenceRecord {
        header,
        bytes: record,
    }))
}

fn parse_id(line: &[u8]) -> &str {
    if line.len() <= 1 {
        return "";
    }

    let end = line.len()
        - line
            .iter()
            .rev()
            .take_while(|b| **b == b'\n' || **b == b'\r')
            .count();
    let payload = &line[1..end];

    let delimiter_pos = memchr(b' ', payload)
        .or_else(|| memchr(b'\t', payload))
        .or_else(|| memchr(b'/', payload))
        .unwrap_or(payload.len());

    std::str::from_utf8(&payload[..delimiter_pos]).unwrap_or("")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_id_strips_pair_suffix() {
        assert_eq!(parse_id(b"@read1/1\n"), "read1");
        assert_eq!(parse_id(b">read2 description\n"), "read2");
    }
}
