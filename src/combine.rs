use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::krk_parser::{self, TaxonEntry};

const BUFFER_SIZE: usize = 256 * 1024;

/// Combine multiple Kraken2 reports by summing clade/direct reads per taxid.
pub fn combine_kreports(input_files: &[String], output_file: &str) -> Result<(), Box<dyn Error>> {
    if input_files.is_empty() {
        return Err("at least one input report is required".into());
    }

    let mut combined: HashMap<u32, CombinedTaxon> = HashMap::new();
    let mut sample_count = 0_u64;

    for file in input_files {
        let (report, _) = krk_parser::parse_kraken2_report(file)?;
        sample_count += 1;
        merge_report(&report.root, &mut combined, None);
        if let Some(ref unclassified) = report.unclassified {
            merge_report(unclassified, &mut combined, None);
        }
    }

    let mut entries: Vec<_> = combined.into_values().collect();
    entries.sort_by(|left, right| {
        left.depth
            .cmp(&right.depth)
            .then_with(|| left.taxid.cmp(&right.taxid))
    });

    let total_reads = entries
        .iter()
        .filter(|entry| entry.depth == 0)
        .map(|entry| entry.clade_reads)
        .sum::<u64>()
        .max(1);

    let file = File::create(output_file)?;
    let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);
    writeln!(
        writer,
        "# Combined from {sample_count} Kraken2 reports by KrakenClip"
    )?;

    for entry in entries {
        let percentage = (entry.clade_reads as f64 / total_reads as f64) * 100.0;
        let indent = "  ".repeat(entry.depth);
        writeln!(
            writer,
            "{percentage:.2}\t{}\t{}\t{}\t{}\t{}{}",
            entry.clade_reads, entry.direct_reads, entry.rank, entry.taxid, indent, entry.name
        )?;
    }

    writer.flush()?;
    Ok(())
}

#[derive(Clone)]
struct CombinedTaxon {
    taxid: u32,
    name: String,
    rank: String,
    depth: usize,
    clade_reads: u64,
    direct_reads: u64,
}

fn merge_report(
    node: &TaxonEntry,
    combined: &mut HashMap<u32, CombinedTaxon>,
    parent_depth: Option<usize>,
) {
    let depth = parent_depth.map(|d| d + 1).unwrap_or(node.depth);
    combined
        .entry(node.taxid)
        .and_modify(|existing| {
            existing.clade_reads = existing.clade_reads.saturating_add(node.clade_reads);
            existing.direct_reads = existing.direct_reads.saturating_add(node.direct_reads);
        })
        .or_insert(CombinedTaxon {
            taxid: node.taxid,
            name: node.name.clone(),
            rank: node.rank.clone(),
            depth,
            clade_reads: node.clade_reads,
            direct_reads: node.direct_reads,
        });

    for child in &node.children {
        merge_report(child, combined, Some(depth));
    }
}
