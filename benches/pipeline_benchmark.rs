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
use krakenclip::bracken::{filter_bracken, parse_bracken_file};
use krakenclip::combine::{combine_kreports, CombineOptions};
use krakenclip::diversity::{
    alpha_from_bracken, beta_matrix, load_sample_counts_bracken, AlphaMetric,
};
use krakenclip::generate_test_data::generate_data;
use krakenclip::kreport_builder::make_kreport;
use krakenclip::krk_parser::parse_kraken2_report;
use krakenclip::lineage::{krona_rows, mpa_rows, LineageOptions};
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

fn write_mini_taxonomy(path: &std::path::Path) {
    fs::write(
        path,
        "1\t|\t1\t|\tR\t|\t0\t|\troot\n2\t|\t1\t|\tD\t|\t1\t|\tBacteria\n3\t|\t2\t|\tS\t|\t2\t|\tSpecies alpha\n4\t|\t2\t|\tS\t|\t2\t|\tSpecies beta\n0\t|\t0\t|\tU\t|\t0\t|\tunclassified\n",
    )
    .unwrap();
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

    let large = directory.path().join("large_report.txt");
    generate_data(large.to_str().unwrap(), 100_000, "dense").unwrap();
    let (report, _) = parse_kraken2_report(large.to_str().unwrap()).unwrap();
    c.bench_function("mpa lineage export 100k nodes", |b| {
        b.iter(|| black_box(mpa_rows(black_box(&report), LineageOptions::default()).len()))
    });
    c.bench_function("krona lineage export 100k nodes", |b| {
        b.iter(|| black_box(krona_rows(black_box(&report), LineageOptions::default()).len()))
    });
}

fn benchmark_new_features(c: &mut Criterion) {
    let directory = tempdir().unwrap();
    let bracken = directory.path().join("sample.bracken");
    generate_data(bracken.to_str().unwrap(), 100_000, "bracken").unwrap();
    let records = parse_bracken_file(bracken.to_str().unwrap()).unwrap();
    c.bench_function("bracken filter 100k rows", |b| {
        b.iter(|| {
            let exclude = HashSet::from([10_u32, 11, 12]);
            black_box(filter_bracken(black_box(&records), None, Some(&exclude)).len())
        })
    });
    c.bench_function("alpha shannon 100k bracken rows", |b| {
        b.iter(|| {
            black_box(alpha_from_bracken(
                black_box(&records),
                AlphaMetric::Shannon,
            ))
        })
    });

    let mut bracken_files = Vec::new();
    for idx in 0..8 {
        let path = directory.path().join(format!("b{idx}.bracken"));
        generate_data(path.to_str().unwrap(), 5_000, "bracken").unwrap();
        bracken_files.push(path);
    }
    let samples: Vec<_> = bracken_files
        .iter()
        .enumerate()
        .map(|(idx, path)| {
            load_sample_counts_bracken(path.to_str().unwrap(), &format!("s{idx}")).unwrap()
        })
        .collect();
    c.bench_function("beta bray-curtis 8 samples", |b| {
        b.iter(|| black_box(beta_matrix(black_box(&samples)).len()))
    });

    let mut reports = Vec::new();
    for idx in 0..16 {
        let path = directory.path().join(format!("c{idx}.txt"));
        generate_data(path.to_str().unwrap(), 20_000, "dense").unwrap();
        reports.push(path.to_string_lossy().to_string());
    }
    let combined = directory.path().join("combined.txt");
    c.bench_function("combine-kreports 16x20k", |b| {
        b.iter(|| {
            combine_kreports(
                black_box(&reports),
                black_box(combined.to_str().unwrap()),
                black_box(&CombineOptions::default()),
            )
            .unwrap()
        })
    });

    let log = directory.path().join("reads.log");
    generate_data(log.to_str().unwrap(), 1_000_000, "log").unwrap();
    let taxonomy = directory.path().join("tax.txt");
    write_mini_taxonomy(&taxonomy);
    let made = directory.path().join("made.kreport");
    c.bench_function("make-kreport 1M log lines", |b| {
        b.iter(|| {
            make_kreport(
                black_box(log.to_str().unwrap()),
                black_box(taxonomy.to_str().unwrap()),
                black_box(made.to_str().unwrap()),
                false,
            )
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
    benchmark_new_features,
    benchmark_gzip_helper
);
criterion_main!(benches);
