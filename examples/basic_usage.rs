//! Example showing how to use KrakenClip as a Rust library.
//!
//! Run with:
//! `cargo run --example basic_usage -- tests/fixtures/report_a.txt`

use std::env;
use std::process;

use krakenclip::abundance_matrix::AbundanceMatrix;
use krakenclip::krk_parser::parse_kraken2_report;
use krakenclip::taxon_query::find_taxon_info;

fn main() {
    let report_path = env::args()
        .nth(1)
        .unwrap_or_else(|| "tests/fixtures/report_a.txt".to_string());

    let (report, parse_seconds) = match parse_kraken2_report(&report_path) {
        Ok(result) => result,
        Err(error) => {
            eprintln!("failed to parse '{report_path}': {error}");
            process::exit(1);
        }
    };

    println!("Parsed {report_path} in {parse_seconds:.6}s");
    println!(
        "Root taxon: {} (taxid {})",
        report.root.name, report.root.taxid
    );
    println!(
        "Indexed taxa with children entries: {}",
        report.index.all_descendants(report.root.taxid).len()
    );

    if let Some(info) = find_taxon_info(&report.root, 3) {
        println!(
            "Found taxon {} with {} parents and {} descendants",
            info.taxon.name,
            info.parents.len(),
            info.children.len()
        );
    }

    let mut matrix = AbundanceMatrix::new("S");
    matrix.add_sample(&report, "example_sample", 0.0, true);
    println!("Species rows available: {}", matrix.rows().len());
}
