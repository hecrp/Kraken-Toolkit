use std::collections::HashSet;
use std::error::Error;
use std::io::{BufRead, Write};

use crate::io_util::{open_buf_reader, open_buf_writer};
use crate::tabular::{field_str, parse_f64_field, parse_u32_field, parse_u64_field, split_tabs};

#[derive(Debug, Clone, PartialEq)]
pub struct BrackenRecord {
    pub name: String,
    pub taxonomy_id: u32,
    pub taxonomy_lvl: String,
    pub kraken_assigned_reads: u64,
    pub added_reads: u64,
    pub new_est_reads: u64,
    pub fraction_total_reads: f64,
}

#[derive(Debug)]
pub enum BrackenError {
    Io(std::io::Error),
    Format(String),
}

impl std::fmt::Display for BrackenError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(error) => write!(f, "I/O error: {error}"),
            Self::Format(message) => write!(f, "{message}"),
        }
    }
}

impl Error for BrackenError {}

impl From<std::io::Error> for BrackenError {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error)
    }
}

pub type BrackenResult<T> = Result<T, BrackenError>;

pub fn parse_bracken_file(path: &str) -> BrackenResult<Vec<BrackenRecord>> {
    let mut reader = open_buf_reader(path)?;
    let mut line = Vec::new();
    let mut records = Vec::new();

    // Skip header if present.
    let header_bytes = reader.read_until(b'\n', &mut line)?;
    if header_bytes == 0 {
        return Ok(records);
    }
    let header = field_str(&line).unwrap_or("").to_ascii_lowercase();
    let has_header = header.contains("taxonomy_id") || header.contains("new_est_reads");
    if !has_header {
        if let Some(record) = parse_bracken_line(&line)? {
            records.push(record);
        }
    }

    loop {
        line.clear();
        let bytes = reader.read_until(b'\n', &mut line)?;
        if bytes == 0 {
            break;
        }
        if let Some(record) = parse_bracken_line(&line)? {
            records.push(record);
        }
    }

    Ok(records)
}

pub fn parse_bracken_line(line: &[u8]) -> BrackenResult<Option<BrackenRecord>> {
    if line.iter().all(u8::is_ascii_whitespace) {
        return Ok(None);
    }
    let (fields, count) = split_tabs::<7>(line);
    if count < 7 {
        return Err(BrackenError::Format(format!(
            "expected 7 Bracken columns, found {count}"
        )));
    }

    Ok(Some(BrackenRecord {
        name: field_str(fields[0])
            .ok_or_else(|| BrackenError::Format("invalid name".into()))?
            .to_string(),
        taxonomy_id: parse_u32_field(fields[1])
            .ok_or_else(|| BrackenError::Format("invalid taxonomy_id".into()))?,
        taxonomy_lvl: field_str(fields[2])
            .ok_or_else(|| BrackenError::Format("invalid taxonomy_lvl".into()))?
            .to_string(),
        kraken_assigned_reads: parse_u64_field(fields[3])
            .ok_or_else(|| BrackenError::Format("invalid kraken_assigned_reads".into()))?,
        added_reads: parse_u64_field(fields[4])
            .ok_or_else(|| BrackenError::Format("invalid added_reads".into()))?,
        new_est_reads: parse_u64_field(fields[5])
            .ok_or_else(|| BrackenError::Format("invalid new_est_reads".into()))?,
        fraction_total_reads: parse_f64_field(fields[6])
            .ok_or_else(|| BrackenError::Format("invalid fraction_total_reads".into()))?,
    }))
}

pub fn filter_bracken(
    records: &[BrackenRecord],
    include: Option<&HashSet<u32>>,
    exclude: Option<&HashSet<u32>>,
) -> Vec<BrackenRecord> {
    let mut filtered: Vec<_> = records
        .iter()
        .filter(|record| {
            if let Some(include) = include {
                include.contains(&record.taxonomy_id)
            } else if let Some(exclude) = exclude {
                !exclude.contains(&record.taxonomy_id)
            } else {
                true
            }
        })
        .cloned()
        .collect();

    let total: u64 = filtered.iter().map(|record| record.new_est_reads).sum();
    for record in &mut filtered {
        record.fraction_total_reads = if total > 0 {
            record.new_est_reads as f64 / total as f64
        } else {
            0.0
        };
    }
    filtered.sort_by(|left, right| {
        right
            .new_est_reads
            .cmp(&left.new_est_reads)
            .then_with(|| left.taxonomy_id.cmp(&right.taxonomy_id))
    });
    filtered
}

pub fn write_bracken_file(path: &str, records: &[BrackenRecord]) -> BrackenResult<()> {
    let mut writer = open_buf_writer(path)?;
    writeln!(
        writer,
        "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads"
    )?;
    for record in records {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.10}",
            record.name,
            record.taxonomy_id,
            record.taxonomy_lvl,
            record.kraken_assigned_reads,
            record.added_reads,
            record.new_est_reads,
            record.fraction_total_reads
        )?;
    }
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn filters_and_renormalizes_fractions() {
        let records = vec![
            BrackenRecord {
                name: "Host".into(),
                taxonomy_id: 9606,
                taxonomy_lvl: "S".into(),
                kraken_assigned_reads: 10,
                added_reads: 0,
                new_est_reads: 10,
                fraction_total_reads: 0.1,
            },
            BrackenRecord {
                name: "Bug".into(),
                taxonomy_id: 562,
                taxonomy_lvl: "S".into(),
                kraken_assigned_reads: 90,
                added_reads: 0,
                new_est_reads: 90,
                fraction_total_reads: 0.9,
            },
        ];
        let exclude = HashSet::from([9606_u32]);
        let filtered = filter_bracken(&records, None, Some(&exclude));
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].taxonomy_id, 562);
        assert!((filtered[0].fraction_total_reads - 1.0).abs() < 1e-12);
    }
}
