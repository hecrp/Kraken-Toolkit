use std::fs::File;
use std::io::{BufWriter, Write};

use crate::krk_parser::{KrakenReport, TaxonEntry};

const BUFFER_SIZE: usize = 256 * 1024;

#[derive(Debug, Clone, Copy)]
pub struct LineageOptions {
    pub intermediate_ranks: bool,
    pub use_percentages: bool,
    pub replace_spaces: bool,
    pub display_header: bool,
}

impl Default for LineageOptions {
    fn default() -> Self {
        Self {
            intermediate_ranks: false,
            use_percentages: false,
            replace_spaces: true,
            display_header: false,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MpaRow {
    pub path: String,
    pub value: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct KronaRow {
    pub count: u64,
    pub path: Vec<String>,
}

fn is_standard_rank(rank: &str) -> bool {
    matches!(
        rank.chars().next(),
        Some('U' | 'R' | 'D' | 'K' | 'P' | 'C' | 'O' | 'F' | 'G' | 'S')
    ) && (rank.len() == 1
        || (rank.starts_with('S') && rank.chars().nth(1).is_some_and(|c| c.is_ascii_digit()))
        || (rank.starts_with('G') && rank.chars().nth(1).is_some_and(|c| c.is_ascii_digit())))
}

fn mpa_prefix(rank: &str) -> String {
    match rank.chars().next().unwrap_or('x') {
        'U' => "u".to_string(),
        'R' => "r".to_string(),
        'D' | 'K' => "d".to_string(),
        'P' => "p".to_string(),
        'C' => "c".to_string(),
        'O' => "o".to_string(),
        'F' => "f".to_string(),
        'G' => "g".to_string(),
        'S' => {
            if rank.len() > 1 {
                format!("s{}", &rank[1..])
            } else {
                "s".to_string()
            }
        }
        _ => "x".to_string(),
    }
}

fn format_name(name: &str, replace_spaces: bool) -> String {
    if replace_spaces {
        name.replace(' ', "_")
    } else {
        name.to_string()
    }
}

fn lineage_token(rank: &str, name: &str, replace_spaces: bool) -> String {
    format!(
        "{}__{}",
        mpa_prefix(rank),
        format_name(name, replace_spaces)
    )
}

/// Build KrakenTools-compatible MPA rows for one report.
pub fn mpa_rows(report: &KrakenReport, options: LineageOptions) -> Vec<MpaRow> {
    let mut rows = Vec::new();
    let mut stack: Vec<String> = Vec::new();

    if let Some(ref unclassified) = report.unclassified {
        walk_mpa(unclassified, &mut stack, options, &mut rows);
    }
    walk_mpa(&report.root, &mut stack, options, &mut rows);
    rows
}

fn walk_mpa(
    node: &TaxonEntry,
    stack: &mut Vec<String>,
    options: LineageOptions,
    rows: &mut Vec<MpaRow>,
) {
    // Root is traversed for structure but omitted from MetaPhlAn-style paths.
    let include = node.rank != "R"
        && (options.intermediate_ranks || is_standard_rank(&node.rank) || node.rank == "U");
    let pushed = if include {
        let token = lineage_token(&node.rank, &node.name, options.replace_spaces);
        stack.push(token);
        true
    } else {
        false
    };

    if include {
        let value = if options.use_percentages {
            node.percentage as f64
        } else {
            node.clade_reads as f64
        };
        rows.push(MpaRow {
            path: stack.join("|"),
            value,
        });
    }

    for child in &node.children {
        walk_mpa(child, stack, options, rows);
    }

    if pushed {
        stack.pop();
    }
}

/// Build KrakenTools-compatible Krona rows for one report.
pub fn krona_rows(report: &KrakenReport, options: LineageOptions) -> Vec<KronaRow> {
    let mut rows = Vec::new();
    let mut stack: Vec<String> = Vec::new();
    let mut pending_intermediate = 0_u64;

    if let Some(ref unclassified) = report.unclassified {
        walk_krona(
            unclassified,
            &mut stack,
            options,
            &mut rows,
            &mut pending_intermediate,
        );
    }
    walk_krona(
        &report.root,
        &mut stack,
        options,
        &mut rows,
        &mut pending_intermediate,
    );
    rows
}

fn walk_krona(
    node: &TaxonEntry,
    stack: &mut Vec<String>,
    options: LineageOptions,
    rows: &mut Vec<KronaRow>,
    pending_intermediate: &mut u64,
) {
    let standard = is_standard_rank(&node.rank) || node.rank == "U";
    let include = options.intermediate_ranks || standard;
    let pushed = if include {
        stack.push(format_name(&node.name, false));
        true
    } else {
        *pending_intermediate = pending_intermediate.saturating_add(node.direct_reads);
        false
    };

    if include {
        let count = node.direct_reads.saturating_add(*pending_intermediate);
        *pending_intermediate = 0;
        if count > 0 || node.children.is_empty() {
            rows.push(KronaRow {
                count,
                path: stack.clone(),
            });
        }
    }

    for child in &node.children {
        walk_krona(child, stack, options, rows, pending_intermediate);
    }

    if pushed {
        stack.pop();
    }
}

pub fn write_mpa_file(
    report: &KrakenReport,
    output_file: &str,
    sample_name: Option<&str>,
    options: LineageOptions,
) -> std::io::Result<()> {
    let file = File::create(output_file)?;
    let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);
    if options.display_header {
        let name = sample_name.unwrap_or("Sample");
        writeln!(writer, "#Classification\t{name}")?;
    }
    for row in mpa_rows(report, options) {
        writeln!(writer, "{}\t{:.0}", row.path, row.value)?;
    }
    writer.flush()
}

pub fn write_krona_file(
    report: &KrakenReport,
    output_file: &str,
    options: LineageOptions,
) -> std::io::Result<()> {
    let file = File::create(output_file)?;
    let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);
    for row in krona_rows(report, options) {
        write!(writer, "{}", row.count)?;
        for part in &row.path {
            write!(writer, "\t{part}")?;
        }
        writeln!(writer)?;
    }
    writer.flush()
}

/// Merge multiple single-sample MPA tables by lineage path.
pub fn write_combined_mpa(
    reports: &[(String, KrakenReport)],
    output_file: &str,
    options: LineageOptions,
) -> std::io::Result<()> {
    use std::collections::{BTreeMap, BTreeSet};

    let mut paths = BTreeSet::new();
    let mut values: Vec<BTreeMap<String, f64>> = Vec::with_capacity(reports.len());
    for (_, report) in reports {
        let mut map = BTreeMap::new();
        for row in mpa_rows(report, options) {
            paths.insert(row.path.clone());
            map.insert(row.path, row.value);
        }
        values.push(map);
    }

    let file = File::create(output_file)?;
    let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);
    write!(writer, "#Classification")?;
    for (sample, _) in reports {
        write!(writer, "\t{sample}")?;
    }
    writeln!(writer)?;

    for path in paths {
        write!(writer, "{path}")?;
        for map in &values {
            let value = map.get(&path).copied().unwrap_or(0.0);
            write!(writer, "\t{value:.0}")?;
        }
        writeln!(writer)?;
    }
    writer.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::krk_parser::TaxonEntry;

    fn sample_report() -> KrakenReport {
        let alpha = TaxonEntry::new(
            40.0,
            400,
            400,
            "S".to_string(),
            3,
            "Species alpha".to_string(),
            2,
        );
        let beta = TaxonEntry::new(
            30.0,
            300,
            300,
            "S".to_string(),
            4,
            "Species beta".to_string(),
            2,
        );
        let mut bacteria =
            TaxonEntry::new(90.0, 900, 0, "D".to_string(), 2, "Bacteria".to_string(), 1);
        bacteria.children.push(alpha);
        bacteria.children.push(beta);
        let mut root = TaxonEntry::new(90.0, 900, 0, "R".to_string(), 1, "root".to_string(), 0);
        root.children.push(bacteria);
        let unclassified = TaxonEntry::new(
            10.0,
            100,
            100,
            "U".to_string(),
            0,
            "unclassified".to_string(),
            0,
        );
        KrakenReport {
            index: crate::krk_parser::TaxonIndex::from_tree(&root, Some(&unclassified)),
            root,
            unclassified: Some(unclassified),
        }
    }

    #[test]
    fn mpa_emits_full_lineage_paths() {
        let rows = mpa_rows(&sample_report(), LineageOptions::default());
        assert!(rows
            .iter()
            .any(|row| row.path == "u__unclassified" && row.value == 100.0));
        assert!(rows
            .iter()
            .any(|row| row.path == "d__Bacteria|s__Species_alpha" && row.value == 400.0));
    }

    #[test]
    fn krona_uses_direct_reads_and_path() {
        let rows = krona_rows(&sample_report(), LineageOptions::default());
        assert!(rows.iter().any(|row| {
            row.count == 400
                && row.path
                    == vec![
                        "root".to_string(),
                        "Bacteria".to_string(),
                        "Species alpha".to_string(),
                    ]
        }));
    }
}
