use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::ktaxonomy::KTaxonomy;
use crate::logkrk_parser::aggregate_taxid_counts;

const BUFFER_SIZE: usize = 256 * 1024;

#[derive(Debug, Clone)]
struct ReportNode {
    taxid: u32,
    parent: u32,
    rank: String,
    level: usize,
    name: String,
    direct_reads: u64,
    clade_reads: u64,
    children: Vec<u32>,
}

/// Build a Kraken-style report from a classification log and condensed taxonomy.
pub fn make_kreport(
    log_path: &str,
    taxonomy_path: &str,
    output_path: &str,
    use_read_len: bool,
) -> Result<(), Box<dyn Error>> {
    let taxonomy = KTaxonomy::from_file(taxonomy_path)?;
    let direct_counts = aggregate_taxid_counts(log_path, use_read_len)?;

    let mut nodes: HashMap<u32, ReportNode> = HashMap::new();
    for (taxid, node) in &taxonomy.nodes {
        nodes.insert(
            *taxid,
            ReportNode {
                taxid: *taxid,
                parent: node.parent,
                rank: node.rank.clone(),
                level: node.level,
                name: node.name.clone(),
                direct_reads: *direct_counts.get(taxid).unwrap_or(&0),
                clade_reads: 0,
                children: node.children.clone(),
            },
        );
    }

    // Ensure unclassified exists when present in counts.
    if let Some(&count) = direct_counts.get(&0) {
        nodes.entry(0).or_insert(ReportNode {
            taxid: 0,
            parent: 0,
            rank: "U".to_string(),
            level: 0,
            name: "unclassified".to_string(),
            direct_reads: count,
            clade_reads: 0,
            children: Vec::new(),
        });
    }

    // Include any taxids observed in the log but missing from taxonomy as dangling nodes.
    for (&taxid, &count) in &direct_counts {
        nodes.entry(taxid).or_insert(ReportNode {
            taxid,
            parent: 1,
            rank: "S".to_string(),
            level: 1,
            name: format!("taxid_{taxid}"),
            direct_reads: count,
            clade_reads: 0,
            children: Vec::new(),
        });
    }

    let roots: Vec<u32> = nodes
        .values()
        .filter(|node| node.taxid == node.parent || !nodes.contains_key(&node.parent))
        .map(|node| node.taxid)
        .collect();

    for root in &roots {
        propagate_clade(*root, &mut nodes);
    }

    let total_reads = nodes
        .values()
        .filter(|node| node.level == 0 || node.rank == "U" || node.rank == "R")
        .map(|node| node.clade_reads)
        .sum::<u64>()
        .max(1);

    // Prefer writing unclassified then root preorder when available.
    let file = File::create(output_path)?;
    let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);

    if let Some(unclassified) = nodes.get(&0).cloned() {
        write_node(&mut writer, &unclassified, total_reads)?;
    }

    let root_id = if nodes.contains_key(&1) {
        1
    } else {
        roots.into_iter().find(|id| *id != 0).unwrap_or(1)
    };
    write_preorder(&mut writer, root_id, &nodes, total_reads)?;
    writer.flush()?;
    Ok(())
}

fn propagate_clade(taxid: u32, nodes: &mut HashMap<u32, ReportNode>) -> u64 {
    let children = nodes
        .get(&taxid)
        .map(|node| node.children.clone())
        .unwrap_or_default();
    let mut clade = nodes.get(&taxid).map(|node| node.direct_reads).unwrap_or(0);
    for child in children {
        clade = clade.saturating_add(propagate_clade(child, nodes));
    }
    if let Some(node) = nodes.get_mut(&taxid) {
        node.clade_reads = clade;
    }
    clade
}

fn write_preorder(
    writer: &mut BufWriter<File>,
    taxid: u32,
    nodes: &HashMap<u32, ReportNode>,
    total_reads: u64,
) -> Result<(), Box<dyn Error>> {
    let Some(node) = nodes.get(&taxid).cloned() else {
        return Ok(());
    };
    if taxid != 0 {
        write_node(writer, &node, total_reads)?;
    }
    let mut children = node.children;
    children.sort_by_key(|child| {
        std::cmp::Reverse(nodes.get(child).map(|n| n.clade_reads).unwrap_or(0))
    });
    for child in children {
        write_preorder(writer, child, nodes, total_reads)?;
    }
    Ok(())
}

fn write_node(
    writer: &mut BufWriter<File>,
    node: &ReportNode,
    total_reads: u64,
) -> Result<(), Box<dyn Error>> {
    let percentage = (node.clade_reads as f64 / total_reads as f64) * 100.0;
    let indent = "  ".repeat(node.level);
    writeln!(
        writer,
        "{percentage:.2}\t{}\t{}\t{}\t{}\t{}{}",
        node.clade_reads, node.direct_reads, node.rank, node.taxid, indent, node.name
    )?;
    Ok(())
}
