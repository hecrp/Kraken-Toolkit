use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

use memchr::memchr;
use rayon::prelude::*;

// Optimized buffer size constant for efficient I/O operations
const BUFFER_SIZE: usize = 1024 * 1024; // 1MB buffer

/// Process sequence files and extract those matching the specified read IDs
///
/// # Arguments
/// * `input_files` - List of input FASTA/FASTQ files to process
/// * `save_readids` - Set of read IDs to extract or exclude
/// * `output_file` - Path to the output file where matching sequences will be written
/// * `exclude` - If true, excludes the IDs in save_readids; if false, includes them
///
/// # Returns
/// * `Result<(), Box<dyn Error + Send + Sync>>` - Result of the operation
///
/// # Implementation Details
/// Multiple inputs are parsed in parallel and written in input order. A single input is
/// streamed directly to disk to avoid synchronization and unnecessary buffering.
pub fn process_sequence_files(
    input_files: &[String],
    save_readids: &HashSet<String>,
    output_file: &str,
    exclude: bool,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let output = File::create(output_file)?;
    let mut writer = BufWriter::with_capacity(BUFFER_SIZE, output);

    if input_files.len() == 1 {
        process_sequence_file(&input_files[0], save_readids, exclude, &mut writer)?;
    } else {
        let chunks: Result<Vec<Vec<u8>>, io::Error> = input_files
            .par_iter()
            .map(|input_file| {
                let mut output = Vec::new();
                process_sequence_file(input_file, save_readids, exclude, &mut output)?;
                Ok(output)
            })
            .collect();

        for chunk in chunks? {
            writer.write_all(&chunk)?;
        }
    }

    writer.flush()?;
    Ok(())
}

fn process_sequence_file<W: Write>(
    input_file: &str,
    save_readids: &HashSet<String>,
    exclude: bool,
    writer: &mut W,
) -> io::Result<()> {
    let file = File::open(input_file)?;
    let mut reader = BufReader::with_capacity(BUFFER_SIZE, file);
    let mut pending_header: Option<Vec<u8>> = None;

    loop {
        let header = if let Some(header) = pending_header.take() {
            header
        } else {
            let mut line = Vec::new();
            if reader.read_until(b'\n', &mut line)? == 0 {
                break;
            }
            if line.iter().all(u8::is_ascii_whitespace) {
                continue;
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
        let should_write = save_readids.contains(&id) != exclude;
        let mut record = header;

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
                    pending_header = Some(line);
                    break;
                }
                record.extend_from_slice(&line);
            }
        }

        if should_write {
            writer.write_all(&record)?;
        }
    }

    Ok(())
}

/// Extract the sequence ID from a FASTA/FASTQ header line
///
/// # Arguments
/// * `line` - Header line that starts with '>' or '@'
///
/// # Returns
/// * `&str` - Extracted sequence ID
///
/// # Performance Note
/// This function uses the highly optimized memchr library for byte-level
/// searching, avoiding unnecessary UTF-8 validation until the final step.
fn parse_id(line: &[u8]) -> &str {
    if line.len() <= 1 {
        return "";
    }

    // Find the first space or tab after the initial character
    // Using memchr for optimized byte searching instead of iterating character by character
    let delimiter_pos = memchr(b' ', &line[1..])
        .or_else(|| memchr(b'\t', &line[1..]))
        .map(|pos| pos + 1) // Adjust for the offset from &line[1..]
        .unwrap_or(line.len() - 1);

    std::str::from_utf8(&line[1..delimiter_pos]).unwrap_or("")
}
