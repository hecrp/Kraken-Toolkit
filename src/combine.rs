use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use rayon::prelude::*;

use crate::krk_parser::{self, TaxonEntry};

const BUFFER_SIZE: usize = 256 * 1024;

#[derive(Debug, Clone, Default)]
pub struct CombineOptions {
    pub sample_names: Vec<String>,
    pub display_headers: bool,
    pub no_headers: bool,
    pub only_combined: bool,
}

#[derive(Clone)]
struct CombinedTaxon {
    taxid: u32,
    name: String,
    rank: String,
    depth: usize,
    parent: Option<u32>,
    clade_reads: u64,
    direct_reads: u64,
    per_sample: Vec<(u64, u64)>,
    children: Vec<u32>,
}

/// Combine multiple Kraken2 reports.
pub fn combine_kreports(
    input_files: &[String],
    output_file: &str,
    options: &CombineOptions,
) -> Result<(), Box<dyn Error>> {
    if input_files.is_empty() {
        return Err("at least one input report is required".into());
    }

    let sample_names: Vec<String> = if options.sample_names.is_empty() {
        input_files
            .iter()
            .enumerate()
            .map(|(index, file)| {
                Path::new(file)
                    .file_stem()
                    .and_then(|name| name.to_str())
                    .map(|name| name.to_string())
                    .unwrap_or_else(|| format!("S{}", index + 1))
            })
            .collect()
    } else {
        if options.sample_names.len() != input_files.len() {
            return Err("sample-names count must match input report count".into());
        }
        options.sample_names.clone()
    };

    let parsed: Result<Vec<_>, String> = input_files
        .par_iter()
        .map(|file| {
            let (report, _) = krk_parser::parse_kraken2_report(file)
                .map_err(|error| format!("Error parsing '{file}': {error}"))?;
            Ok(report)
        })
        .collect();
    let reports = parsed?;

    let mut combined: HashMap<u32, CombinedTaxon> = HashMap::new();
    for (sample_idx, report) in reports.iter().enumerate() {
        if let Some(ref unclassified) = report.unclassified {
            merge_report(
                unclassified,
                None,
                sample_idx,
                input_files.len(),
                &mut combined,
            );
        }
        merge_report(
            &report.root,
            None,
            sample_idx,
            input_files.len(),
            &mut combined,
        );
    }

    // Build children lists for preorder emission based on first-seen parent links.
    let taxids: Vec<u32> = combined.keys().copied().collect();
    for taxid in taxids {
        if let Some(parent) = combined.get(&taxid).and_then(|node| node.parent) {
            if let Some(parent_node) = combined.get_mut(&parent) {
                if !parent_node.children.contains(&taxid) {
                    parent_node.children.push(taxid);
                }
            }
        }
    }
    for node in combined.values_mut() {
        node.children.sort_by_key(|child| *child);
    }

    let total_reads = combined
        .values()
        .filter(|entry| entry.depth == 0)
        .map(|entry| entry.clade_reads)
        .sum::<u64>()
        .max(1);

    let file = File::create(output_file)?;
    let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);

    if options.only_combined {
        writeln!(
            writer,
            "# Combined from {} Kraken2 reports by KrakenClip",
            input_files.len()
        )?;
        emit_preorder_combined(&mut writer, &combined, 0, true, total_reads)?;
        if combined.contains_key(&1) {
            emit_preorder_combined(&mut writer, &combined, 1, false, total_reads)?;
        }
    } else {
        if !options.no_headers {
            if options.display_headers {
                write!(writer, "#perc\ttot_all\ttot_lvl")?;
                for name in &sample_names {
                    write!(writer, "\t{name}_all\t{name}_lvl")?;
                }
                writeln!(writer, "\trank\ttaxid\tname")?;
            } else {
                write!(writer, "#perc\ttot_all\ttot_lvl")?;
                for index in 0..sample_names.len() {
                    write!(writer, "\tS{}_all\tS{}_lvl", index + 1, index + 1)?;
                }
                writeln!(writer, "\trank\ttaxid\tname")?;
            }
        }

        emit_preorder_multi(&mut writer, &combined, 0, true, total_reads)?;
        if combined.contains_key(&1) {
            emit_preorder_multi(&mut writer, &combined, 1, false, total_reads)?;
        }
    }

    writer.flush()?;
    Ok(())
}

fn merge_report(
    node: &TaxonEntry,
    parent: Option<u32>,
    sample_idx: usize,
    sample_count: usize,
    combined: &mut HashMap<u32, CombinedTaxon>,
) {
    combined
        .entry(node.taxid)
        .and_modify(|existing| {
            existing.clade_reads = existing.clade_reads.saturating_add(node.clade_reads);
            existing.direct_reads = existing.direct_reads.saturating_add(node.direct_reads);
            if let Some(slot) = existing.per_sample.get_mut(sample_idx) {
                slot.0 = slot.0.saturating_add(node.clade_reads);
                slot.1 = slot.1.saturating_add(node.direct_reads);
            }
        })
        .or_insert_with(|| {
            let mut per_sample = vec![(0, 0); sample_count];
            per_sample[sample_idx] = (node.clade_reads, node.direct_reads);
            CombinedTaxon {
                taxid: node.taxid,
                name: node.name.clone(),
                rank: node.rank.clone(),
                depth: node.depth,
                parent,
                clade_reads: node.clade_reads,
                direct_reads: node.direct_reads,
                per_sample,
                children: Vec::new(),
            }
        });

    for child in &node.children {
        merge_report(child, Some(node.taxid), sample_idx, sample_count, combined);
    }
}

fn emit_preorder_combined(
    writer: &mut BufWriter<File>,
    combined: &HashMap<u32, CombinedTaxon>,
    taxid: u32,
    allow_missing: bool,
    total_reads: u64,
) -> Result<(), Box<dyn Error>> {
    let Some(entry) = combined.get(&taxid) else {
        if allow_missing {
            return Ok(());
        }
        return Ok(());
    };
    let percentage = (entry.clade_reads as f64 / total_reads as f64) * 100.0;
    let indent = "  ".repeat(entry.depth);
    writeln!(
        writer,
        "{percentage:.2}\t{}\t{}\t{}\t{}\t{}{}",
        entry.clade_reads, entry.direct_reads, entry.rank, entry.taxid, indent, entry.name
    )?;
    for child in &entry.children {
        emit_preorder_combined(writer, combined, *child, false, total_reads)?;
    }
    Ok(())
}

fn emit_preorder_multi(
    writer: &mut BufWriter<File>,
    combined: &HashMap<u32, CombinedTaxon>,
    taxid: u32,
    allow_missing: bool,
    total_reads: u64,
) -> Result<(), Box<dyn Error>> {
    let Some(entry) = combined.get(&taxid) else {
        if allow_missing {
            return Ok(());
        }
        return Ok(());
    };
    let percentage = (entry.clade_reads as f64 / total_reads as f64) * 100.0;
    write!(
        writer,
        "{percentage:.2}\t{}\t{}",
        entry.clade_reads, entry.direct_reads
    )?;
    for (clade, direct) in &entry.per_sample {
        write!(writer, "\t{clade}\t{direct}")?;
    }
    let indent = "  ".repeat(entry.depth);
    writeln!(
        writer,
        "\t{}\t{}\t{}{}",
        entry.rank, entry.taxid, indent, entry.name
    )?;
    for child in &entry.children {
        emit_preorder_multi(writer, combined, *child, false, total_reads)?;
    }
    Ok(())
}
