use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use flate2::write::GzEncoder;
use flate2::Compression;
use serde_json::Value;
use tempfile::tempdir;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn run(arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_krakenclip"))
        .args(arguments)
        .output()
        .expect("KrakenClip should execute")
}

#[test]
fn help_lists_commands() {
    let output = run(&["--help"]);
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(output.status.success());
    assert!(stdout.contains("analyze"));
    assert!(stdout.contains("extract"));
    assert!(stdout.contains("abundance-matrix"));
    assert!(stdout.contains("combine-kreports"));
}

#[test]
fn analyze_reports_real_metrics_and_rejects_invalid_taxid() {
    let report = fixture("report_a.txt");
    let output = run(&["analyze", report.to_str().unwrap()]);
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(output.status.success());
    assert!(stdout.contains("File parsing time"));
    assert!(!stdout.contains("Hierarchy build time"));

    let invalid = run(&["analyze", report.to_str().unwrap(), "--tax-id", "abc"]);
    assert!(!invalid.status.success());
    assert!(String::from_utf8(invalid.stderr)
        .unwrap()
        .contains("invalid value"));
}

#[test]
fn analyze_propagates_json_write_failures() {
    let directory = tempdir().unwrap();
    let report = fixture("report_a.txt");
    let output = run(&[
        "analyze",
        report.to_str().unwrap(),
        "--json",
        directory.path().to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert!(String::from_utf8(output.stderr)
        .unwrap()
        .contains("Error writing JSON report"));
}

#[test]
fn extracts_complete_fastq_records() {
    let directory = tempdir().unwrap();
    let output_file = directory.path().join("selected.fastq");
    let output = run(&[
        "extract",
        fixture("reads.fastq").to_str().unwrap(),
        fixture("kraken.log").to_str().unwrap(),
        "--output",
        output_file.to_str().unwrap(),
        "--taxids",
        "3",
    ]);
    assert!(output.status.success(), "{:?}", output.stderr);
    assert_eq!(
        fs::read_to_string(output_file).unwrap(),
        "@read1 description\nACGT\n+\n!!!!\n"
    );
}

#[test]
fn fasta_exclusion_preserves_following_headers() {
    let directory = tempdir().unwrap();
    let output_file = directory.path().join("selected.fasta");
    let output = run(&[
        "extract",
        fixture("reads.fasta").to_str().unwrap(),
        fixture("kraken.log").to_str().unwrap(),
        "--output",
        output_file.to_str().unwrap(),
        "--taxids",
        "3",
        "--exclude",
    ]);
    assert!(output.status.success(), "{:?}", output.stderr);
    assert_eq!(
        fs::read_to_string(output_file).unwrap(),
        ">read2\nTGCA\n>read3\nAAAA\n"
    );
}

#[test]
fn rejects_truncated_fastq_records() {
    let directory = tempdir().unwrap();
    let input = directory.path().join("truncated.fastq");
    let output_file = directory.path().join("output.fastq");
    fs::write(&input, "@read1\nACGT\n+\n").unwrap();
    let output = run(&[
        "extract",
        input.to_str().unwrap(),
        fixture("kraken.log").to_str().unwrap(),
        "--output",
        output_file.to_str().unwrap(),
        "--taxids",
        "3",
    ]);
    assert!(!output.status.success());
    assert!(String::from_utf8(output.stderr)
        .unwrap()
        .contains("truncated FASTQ"));
}

#[test]
fn abundance_threshold_matches_output_units() {
    let directory = tempdir().unwrap();
    let proportional = directory.path().join("proportional.tsv");
    let output = run(&[
        "abundance-matrix",
        fixture("report_a.txt").to_str().unwrap(),
        "--output",
        proportional.to_str().unwrap(),
        "--min-abundance",
        "50",
    ]);
    assert!(output.status.success(), "{:?}", output.stderr);
    assert_eq!(
        fs::read_to_string(proportional).unwrap(),
        "Taxon\treport_a\n"
    );

    let absolute = directory.path().join("absolute.tsv");
    let output = run(&[
        "abundance-matrix",
        fixture("report_a.txt").to_str().unwrap(),
        "--output",
        absolute.to_str().unwrap(),
        "--absolute-counts",
        "--min-abundance",
        "350",
    ]);
    assert!(output.status.success(), "{:?}", output.stderr);
    assert!(fs::read_to_string(absolute)
        .unwrap()
        .contains("Species alpha\t400.000000"));
}

#[test]
fn biom_contains_all_samples_with_sparse_data() {
    let directory = tempdir().unwrap();
    let output_file = directory.path().join("matrix.biom");
    let output = run(&[
        "abundance-matrix",
        fixture("report_a.txt").to_str().unwrap(),
        fixture("report_b.txt").to_str().unwrap(),
        "--output",
        output_file.to_str().unwrap(),
        "--format",
        "biom",
    ]);
    assert!(output.status.success(), "{:?}", output.stderr);
    let biom: Value = serde_json::from_slice(&fs::read(output_file).unwrap()).unwrap();
    assert_eq!(biom["shape"], serde_json::json!([3, 2]));
    assert_eq!(biom["matrix_type"], "sparse");
    assert_eq!(biom["columns"].as_array().unwrap().len(), 2);
    assert!(biom["data"]
        .as_array()
        .unwrap()
        .iter()
        .all(|entry| entry.as_array().unwrap().len() == 3));
}

#[test]
fn biom_dense_contains_row_major_values() {
    let directory = tempdir().unwrap();
    let output_file = directory.path().join("matrix.biom");
    let output = run(&[
        "abundance-matrix",
        fixture("report_a.txt").to_str().unwrap(),
        fixture("report_b.txt").to_str().unwrap(),
        "--output",
        output_file.to_str().unwrap(),
        "--format",
        "biom",
        "--biom-matrix-type",
        "dense",
    ]);
    assert!(output.status.success(), "{:?}", output.stderr);
    let biom: Value = serde_json::from_slice(&fs::read(output_file).unwrap()).unwrap();
    assert_eq!(biom["matrix_type"], "dense");
    assert!(biom["data"]
        .as_array()
        .unwrap()
        .iter()
        .all(|row| row.as_array().unwrap().len() == 2));
}

#[test]
fn kingdom_alias_selects_kraken_domain_rank() {
    let directory = tempdir().unwrap();
    let output_file = directory.path().join("domain.tsv");
    let output = run(&[
        "abundance-matrix",
        fixture("report_a.txt").to_str().unwrap(),
        "--output",
        output_file.to_str().unwrap(),
        "--level",
        "K",
    ]);
    assert!(output.status.success(), "{:?}", output.stderr);
    assert!(fs::read_to_string(output_file)
        .unwrap()
        .contains("Bacteria\t90.000000"));
}

#[test]
fn generated_report_has_requested_lines_and_valid_hierarchy() {
    let directory = tempdir().unwrap();
    let report = directory.path().join("generated.txt");
    let output = run(&[
        "generate-test-data",
        "--output",
        report.to_str().unwrap(),
        "--lines",
        "10",
        "--type",
        "dense",
    ]);
    assert!(output.status.success(), "{:?}", output.stderr);
    assert_eq!(fs::read_to_string(&report).unwrap().lines().count(), 10);

    let analysis = run(&["analyze", report.to_str().unwrap(), "--tax-id", "1"]);
    assert!(analysis.status.success(), "{:?}", analysis.stderr);
}

#[test]
fn extract_include_children_uses_report_index() {
    let directory = tempdir().unwrap();
    let output_file = directory.path().join("selected.fastq");
    let stats_file = directory.path().join("stats.csv");
    let output = run(&[
        "extract",
        fixture("reads.fastq").to_str().unwrap(),
        fixture("kraken.log").to_str().unwrap(),
        "--output",
        output_file.to_str().unwrap(),
        "--taxids",
        "2",
        "--include-children",
        "--report",
        fixture("report_a.txt").to_str().unwrap(),
        "--stats-output",
        stats_file.to_str().unwrap(),
    ]);
    assert!(output.status.success(), "{:?}", output.stderr);
    let selected = fs::read_to_string(output_file).unwrap();
    assert!(selected.contains("@read1"));
    assert!(selected.contains("@read2"));
    let stats = fs::read_to_string(stats_file).unwrap();
    assert!(stats.contains("taxid,sequences"));
}

#[test]
fn extract_supports_gzip_and_paired_outputs() {
    let directory = tempdir().unwrap();
    let r1 = directory.path().join("reads_r1.fastq.gz");
    let r2 = directory.path().join("reads_r2.fastq.gz");
    let out1 = directory.path().join("out_r1.fastq.gz");
    let out2 = directory.path().join("out_r2.fastq");

    for (path, suffix) in [(&r1, "/1"), (&r2, "/2")] {
        let file = fs::File::create(path).unwrap();
        let mut encoder = GzEncoder::new(file, Compression::fast());
        write!(
            encoder,
            "@read1{suffix}\nACGT\n+\n!!!!\n@read2{suffix}\nTGCA\n+\n!!!!\n"
        )
        .unwrap();
        encoder.finish().unwrap();
    }

    let output = run(&[
        "extract",
        r1.to_str().unwrap(),
        fixture("kraken.log").to_str().unwrap(),
        "--sequence2",
        r2.to_str().unwrap(),
        "--output",
        out1.to_str().unwrap(),
        "--output2",
        out2.to_str().unwrap(),
        "--taxids",
        "3",
    ]);
    assert!(output.status.success(), "{:?}", output.stderr);
    assert!(out1.exists());
    assert_eq!(
        fs::read_to_string(out2).unwrap(),
        "@read1/2\nACGT\n+\n!!!!\n"
    );
}

#[test]
fn combine_kreports_and_mpa_formats() {
    let directory = tempdir().unwrap();
    let combined = directory.path().join("combined.txt");
    let output = run(&[
        "combine-kreports",
        fixture("report_a.txt").to_str().unwrap(),
        fixture("report_b.txt").to_str().unwrap(),
        "--output",
        combined.to_str().unwrap(),
    ]);
    assert!(output.status.success(), "{:?}", output.stderr);
    let contents = fs::read_to_string(combined).unwrap();
    assert!(contents.contains("Species alpha"));

    let mpa = directory.path().join("matrix.mpa");
    let output = run(&[
        "abundance-matrix",
        fixture("report_a.txt").to_str().unwrap(),
        "--output",
        mpa.to_str().unwrap(),
        "--format",
        "mpa",
    ]);
    assert!(output.status.success(), "{:?}", output.stderr);
    assert!(fs::read_to_string(mpa)
        .unwrap()
        .contains("s__Species_alpha"));
}
