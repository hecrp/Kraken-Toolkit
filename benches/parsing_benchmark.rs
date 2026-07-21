use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use krakenclip::krk_parser::parse_kraken2_report;

fn benchmark_parsing(c: &mut Criterion) {
    c.bench_function("parse kraken2 report", |b| {
        b.iter(|| {
            parse_kraken2_report(black_box("tests/fixtures/report_a.txt"))
                .expect("benchmark fixture should parse")
        })
    });
}

criterion_group!(benches, benchmark_parsing);
criterion_main!(benches);
