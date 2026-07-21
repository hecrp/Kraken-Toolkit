use crate::krk_parser::TaxonEntry;
use colored::*;

pub struct TaxonInfo {
    pub taxon: TaxonEntry,
    pub parents: Vec<TaxonEntry>,
    pub children: Vec<TaxonEntry>,
}

pub fn find_taxon_info(root: &TaxonEntry, target_tax_id: u64) -> Option<TaxonInfo> {
    fn find_taxon_recursive(
        entry: &TaxonEntry,
        target_tax_id: u64,
        parents: &mut Vec<TaxonEntry>,
    ) -> Option<TaxonInfo> {
        if entry.taxon_id == target_tax_id {
            let mut all_children = Vec::new();
            collect_children(&entry.children, &mut all_children);

            return Some(TaxonInfo {
                taxon: entry.clone(),
                parents: parents.clone(),
                children: all_children,
            });
        }

        for child in &entry.children {
            parents.push(entry.clone());
            if let Some(info) = find_taxon_recursive(child, target_tax_id, parents) {
                return Some(info);
            }
            parents.pop();
        }

        None
    }

    fn collect_children(entries: &[TaxonEntry], all_children: &mut Vec<TaxonEntry>) {
        for entry in entries {
            all_children.push(entry.clone());
            collect_children(&entry.children, all_children);
        }
    }

    let mut parents = Vec::new();
    find_taxon_recursive(root, target_tax_id, &mut parents)
}

pub fn print_taxon_info(info: &TaxonInfo) {
    println!("{}", "Taxon Information:".bold());
    println!("ID: {}", info.taxon.taxon_id.to_string().cyan());
    println!("Name: {}", info.taxon.name.yellow());
    println!("Level: {}", info.taxon.level);
    println!(
        "Percentage: {}%",
        info.taxon.percentage.to_string().magenta()
    );
    println!(
        "Clade Fragments: {}",
        info.taxon.clade_fragments.to_string().blue()
    );
    println!(
        "Direct Fragments: {}",
        info.taxon.direct_fragments.to_string().blue()
    );
    println!("Rank Code: {}", info.taxon.rank_code.red());

    println!(
        "\n{}",
        "Taxonomic Subtree - C clade fragments - D direct fragments".bold()
    );
    print_taxonomic_tree(&info.parents, &info.taxon, &info.children);
}

fn print_taxonomic_tree(parents: &[TaxonEntry], current: &TaxonEntry, children: &[TaxonEntry]) {
    for (i, parent) in parents.iter().enumerate() {
        println!(
            "{}├── {}: {} (C{}) (D{})",
            "│   ".repeat(i),
            parent.taxon_id,
            parent.name,
            parent.clade_fragments,
            parent.direct_fragments
        );
    }

    let parent_depth = parents.len();
    println!(
        "{}└── {}: {} (C{}) (D{})",
        "│   ".repeat(parent_depth),
        current.taxon_id,
        current.name.green(),
        current.clade_fragments,
        current.direct_fragments
    );

    print_children_tree(children, &"│   ".repeat(parent_depth + 1));
}

fn print_children_tree(children: &[TaxonEntry], prefix: &str) {
    for (i, child) in children.iter().enumerate() {
        let is_last = i == children.len() - 1;
        let current_prefix = if is_last {
            format!("{}└── ", prefix)
        } else {
            format!("{}├── ", prefix)
        };
        println!(
            "{}{}: {} (C{}) (D{})",
            current_prefix,
            child.taxon_id,
            child.name,
            child.clade_fragments,
            child.direct_fragments
        );

        let new_prefix = if is_last {
            format!("{}    ", prefix)
        } else {
            format!("{}│   ", prefix)
        };
        print_children_tree(&child.children, &new_prefix);
    }
}
