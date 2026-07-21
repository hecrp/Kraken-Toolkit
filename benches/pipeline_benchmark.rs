use std::collections::HashSet;
use std::fs;
use std::hint::black_box;
use std::io::Write;
use std::path::PathBuf;

use criterion::{criterion_group, criterion_main, Criterion};
use flate2::write::GzEncoder;
use flate2::Compression;
use krakenclip::abundance_matrix::AbundanceMatrix;
use krakenclip::biom::{BiomMatrixType, BiomTable};
use krakenclip::generate_test_data::generate_data;
use krakenclip::krk_parser::parse_kraken2_report;
use krakenclip::logkrk_parser::parse_kraken_output;
use krakenclip::sequence_processor::process_sequence_files;
use tempfile::tempdir;

fn write_synthetic_extract_inputs(dir: &std::path::Path, reads: usize) -> (PathBuf, PathBuf) {
    let log_path = dir.join("reads.kraken");
    let fq_path = dir.join("reads.fastq");
    let mut log = fs::File::create(&log_path).unwrap();
    let mut fq = fs::File::create(&fq_path).unwrap();

    for i in 0..reads {
        let taxid = if i % 10 == 0 { 3 } else { 4 };
        writeln!(log, "C\tread{i}\t{taxid}\t4\t{taxid}:4").unwrap();
        writeln!(fq, "@read{i}").unwrap();
        writeln!(fq, "ACGT").unwrap();
        writeln!(fq, "+").unwrap();
        writeln!(fq, "!!!!").unwrap();
    }

    (log_path, fq_path)
}

fn benchmark_log_and_extract(c: &mut Criterion) {
    let directory = tempdir().unwrap();
    let (log_path, fq_path) = write_synthetic_extract_inputs(directory.path(), 50_000);
    let taxids = HashSet::from([3_u32]);

    c.bench_function("parse kraken log 50k reads", |b| {
        b.iter(|| {
            parse_kraken_output(black_box(log_path.to_str().unwrap()), black_box(&taxids))
                .expect("log parse")
        })
    });

    let readids = parse_kraken_output(log_path.to_str().unwrap(), &taxids).unwrap();
    let output = directory.path().join("out.fastq");
    c.bench_function("extract matching fastq 50k reads", |b| {
        b.iter(|| {
            process_sequence_files(
                black_box(&[fq_path.to_string_lossy().to_string()]),
                black_box(&readids),
                black_box(output.to_str().unwrap()),
                false,
            )
            .expect("extract")
        })
    });

    // Gzip round-trip path for workflow compatibility.
    let gz_out = directory.path().join("out.fastq.gz");
    c.bench_function("extract matching fastq.gz 50k reads", |b| {
        b.iter(|| {
            process_sequence_files(
                black_box(&[fq_path.to_string_lossy().to_string()]),
                black_box(&readids),
                black_box(gz_out.to_str().unwrap()),
                false,
            )
            .expect("extract gz")
        })
    });
}

fn benchmark_abundance_and_biom(c: &mut Criterion) {
    let directory = tempdir().unwrap();
    let mut reports = Vec::new();
    for idx in 0..8 {
        let path = directory.path().join(format!("sample_{idx}.txt"));
        generate_data(path.to_str().unwrap(), 20_000, "dense").unwrap();
        reports.push(path);
    }

    c.bench_function("abundance-matrix 8x20k reports", |b| {
        b.iter(|| {
            let mut matrix = AbundanceMatrix::new("S");
            for (idx, path) in reports.iter().enumerate() {
                let (report, _) = parse_kraken2_report(path.to_str().unwrap()).unwrap();
                matrix.add_sample(&report, &format!("s{idx}"), 0.0, true);
            }
            black_box(matrix.rows().len())
        })
    });

    let mut matrix = AbundanceMatrix::new("S");
    for (idx, path) in reports.iter().enumerate() {
        let (report, _) = parse_kraken2_report(path.to_str().unwrap()).unwrap();
        matrix.add_sample(&report, &format!("s{idx}"), 0.0, true);
    }
    let biom_path = directory.path().join("matrix.biom");
    c.bench_function("biom sparse serialize", |b| {
        b.iter(|| {
            BiomTable::from_abundance_matrix(black_box(&matrix), BiomMatrixType::Sparse)
                .write_json(biom_path.to_str().unwrap())
                .unwrap()
        })
    });
}

fn benchmark_gzip_helper(c: &mut Criterion) {
    let directory = tempdir().unwrap();
    let path = directory.path().join("tiny.fastq.gz");
    let file = fs::File::create(&path).unwrap();
    let mut encoder = GzEncoder::new(file, Compression::fast());
    writeln!(encoder, "@read1").unwrap();
    writeln!(encoder, "ACGT").unwrap();
    writeln!(encoder, "+").unwrap();
    writeln!(encoder, "!!!!").unwrap();
    encoder.finish().unwrap();

    let readids = HashSet::from(["read1".to_string()]);
    let output = directory.path().join("out.fastq");
    c.bench_function("extract from gzip fastq", |b| {
        b.iter(|| {
            process_sequence_files(
                black_box(&[path.to_string_lossy().to_string()]),
                black_box(&readids),
                black_box(output.to_str().unwrap()),
                false,
            )
            .unwrap()
        })
    });
}

criterion_group!(
    benches,
    benchmark_log_and_extract,
    benchmark_abundance_and_biom,
    benchmark_gzip_helper
);
criterion_main!(benches);
