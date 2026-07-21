use std::hint::black_box;
use std::path::PathBuf;

use criterion::{criterion_group, criterion_main, Criterion};
use krakenclip::generate_test_data::generate_data;
use krakenclip::krk_parser::parse_kraken2_report;
use tempfile::tempdir;

fn large_report_path(lines: usize) -> (tempfile::TempDir, PathBuf) {
    let directory = tempdir().expect("tempdir");
    let path = directory.path().join(format!("report_{lines}.txt"));
    generate_data(path.to_str().unwrap(), lines, "dense").expect("generate report");
    (directory, path)
}

fn benchmark_parsing(c: &mut Criterion) {
    c.bench_function("parse kraken2 report fixture", |b| {
        b.iter(|| {
            parse_kraken2_report(black_box("tests/fixtures/report_a.txt"))
                .expect("benchmark fixture should parse")
        })
    });

    let (_dir_100k, path_100k) = large_report_path(100_000);
    let path_100k = path_100k.to_string_lossy().to_string();
    c.bench_function("parse kraken2 report 100k lines", |b| {
        b.iter(|| parse_kraken2_report(black_box(&path_100k)).expect("100k report should parse"))
    });
}

criterion_group!(benches, benchmark_parsing);
criterion_main!(benches);
