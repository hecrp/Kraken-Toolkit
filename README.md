# KrakenClip

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Rust](https://img.shields.io/badge/rust-1.86%2B-blue.svg)](https://www.rust-lang.org/)
[![Docker Pulls](https://img.shields.io/docker/pulls/hecrp/krakenclip.svg)](https://hub.docker.com/r/hecrp/krakenclip)
[![Docker Image Size](https://img.shields.io/docker/image-size/hecrp/krakenclip.svg)](https://hub.docker.com/r/hecrp/krakenclip)
[![CI](https://github.com/hecrp/KrakenClip/actions/workflows/ci.yml/badge.svg)](https://github.com/hecrp/KrakenClip/actions/workflows/ci.yml)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/hecrp/krakenclip/graphs/commit-activity)

KrakenClip is a high-performance command-line toolkit written in [Rust](https://www.rust-lang.org/) for processing [Kraken2](https://ccb.jhu.edu/software/kraken2/) reports, classification logs, Bracken tables, and sequence files. It is designed as a fast, dependency-free alternative for common post-processing workflows used in metagenomic pipelines.

KrakenClip implements functionality inspired by [KrakenTools](https://github.com/jenniferlu717/KrakenTools). Credit for the original ideas and workflows belongs to KrakenTools; if you use KrakenClip in your research, please also cite:

> Lu J, Rincon N, Wood DE, Breitwieser FP, Pockrandt C, Langmead B, Salzberg SL, Steinegger M. Metagenome analysis using the Kraken software suite. Nature Protocols, doi: 10.1038/s41596-022-00738-y (2022)

## Features

- **Standalone binary** with no runtime Python/Biopython/NumPy dependencies
- **Fast report parsing** with indexed taxonomy and flat streaming iterators
- **Sequence extraction** from FASTA/FASTQ, including gzip, paired-end, and `--max`
- **Hierarchical MPA/Krona exports** with lineage paths compatible with KrakenTools semantics
- **Multi-sample abundance matrices** in TSV and BIOM (sparse/dense)
- **Parallel report combination** (`combine-kreports`) with per-sample columns or aggregated output
- **Bracken filtering** and **alpha/beta diversity** (Shannon, Berger–Parker, Simpson, Fisher, Bray–Curtis)
- **`make-kreport`** from classification logs + condensed taxonomy
- **Criterion benches** and synthetic generators for reports, Bracken tables, and logs
- **Docker multi-arch image** for x86_64 and ARM64

## Getting Started

### Prerequisites

- Rust 1.86 or newer
- Cargo

### Install from source

```bash
git clone https://github.com/hecrp/krakenclip.git
cd krakenclip
cargo build --release --locked
./target/release/krakenclip --help
```

### Docker

```bash
docker build -t krakenclip .
docker run --rm -v /path/to/local/data:/data krakenclip analyze /data/report.txt
docker pull hecrp/krakenclip:latest
```

## Quick examples

```bash
# Inspect a Kraken2 report
krakenclip analyze sample.kreport --tax-id 562 --json sample.json

# Extract reads (optionally hierarchical / paired / gzip)
krakenclip extract reads.fastq.gz kraken.log \
  --output selected.fastq.gz \
  --taxids 2 \
  --include-children \
  --report sample.kreport

# Hierarchical MPA / Krona
krakenclip abundance-matrix sample.kreport --output sample.mpa --format mpa
krakenclip abundance-matrix sample.kreport --output sample.krona --format krona

# Combine reports with per-sample columns
krakenclip combine-kreports a.kreport b.kreport \
  --output combined.txt \
  --display-headers

# Bracken filter + diversity
krakenclip filter-bracken --input sample.bracken --output no_host.bracken --exclude 9606
krakenclip alpha-diversity --input no_host.bracken --type shannon
krakenclip beta-diversity --input s1.bracken s2.bracken --output beta.tsv --type bracken

# Rebuild a kreport from a classification log
krakenclip make-kreport \
  --log sample.kraken \
  --taxonomy taxonomy.txt \
  --output rebuilt.kreport
```

## Usage

```text
Commands:
  analyze             Analyzes a Kraken2 report
  extract             Extracts sequences based on Kraken2 results
  abundance-matrix    Generates taxonomic abundance matrices / MPA / Krona
  combine-kreports    Combines multiple Kraken2 reports into one
  filter-bracken      Filters Bracken abundance tables
  alpha-diversity     Computes alpha diversity from Bracken output
  beta-diversity      Computes Bray–Curtis beta diversity
  make-kreport        Builds a Kraken report from a log + taxonomy
  generate-test-data  Generates synthetic reports/logs/Bracken tables
```

### `analyze`

```bash
krakenclip analyze sample.kreport --tax-id 562 --json sample.json
```

### `extract`

Supports gzip inputs/outputs, paired-end (`--sequence2/--output2`), hierarchy (`--include-children/--include-parents` with `--report`), exclusion, statistics, and `--max`.

### `abundance-matrix`

| Format | Behavior |
|--------|----------|
| `tsv` | Taxon × sample matrix at a selected rank (`--level`) |
| `biom` | BIOM 1.0.0 JSON; sparse by default (`--biom-matrix-type`) |
| `mpa` | Full lineage paths (`d__Bacteria\|s__Escherichia_coli`), multi-sample merge by path |
| `krona` | `direct_reads` + taxonomy path fragments (single report) |

MPA/Krona options: `--intermediate-ranks`, `--percentages`, `--keep-spaces`, `--display-header`.

### `combine-kreports`

- Default: multi-sample table with `tot_all`/`tot_lvl` and per-sample columns
- `--only-combined`: classic aggregated kreport
- `--sample-names`, `--display-headers`, `--no-headers`
- Reports are parsed in parallel and emitted in taxonomic preorder

### `filter-bracken`

```bash
krakenclip filter-bracken --input sample.bracken --output filtered.bracken --include 562,1280
krakenclip filter-bracken --input sample.bracken --output filtered.bracken --exclude 9606
```

Recalculates `fraction_total_reads` over retained `new_est_reads` and sorts by abundance.

### `alpha-diversity`

Metrics: `shannon`, `berger-parker`, `simpson`, `inverse-simpson`, `fisher`.

```bash
krakenclip alpha-diversity --input sample.bracken --type shannon
```

### `beta-diversity`

```bash
krakenclip beta-diversity --input a.bracken b.bracken --output beta.tsv --type bracken
krakenclip beta-diversity --input a.kreport b.kreport --output beta.tsv --type kreport --level S
krakenclip beta-diversity --input table.tsv --output beta.tsv --type tsv --cols 0,1
```

Samples are loaded in parallel; Bray–Curtis pairs are computed on the upper triangle.

### `make-kreport`

Requires a condensed taxonomy file in `make_ktaxonomy.py` format (`taxid\t|\tparent\t|\trank\t|\tlevel\t|\tname`).

```bash
krakenclip make-kreport --log sample.kraken --taxonomy taxonomy.txt -o sample.kreport
krakenclip make-kreport --log sample.kraken --taxonomy taxonomy.txt -o sample.kreport --use-read-len
```

Counts are aggregated from the log without storing read IDs.

### `generate-test-data`

```bash
krakenclip generate-test-data -o report.txt -l 100000 -t dense
krakenclip generate-test-data -o sample.bracken -l 100000 -t bracken
krakenclip generate-test-data -o sample.log -l 1000000 -t log
```

## Performance

Design invariants:

- streaming I/O with large buffers and optional gzip
- numeric taxids and early filtering
- flat kreport iteration when a full tree is unnecessary
- Rayon only for multi-file workloads (abundance, combine, beta)
- no heavy numeric dependencies

```bash
cargo bench --bench parsing_benchmark -- --quick
cargo bench --bench pipeline_benchmark -- --quick
```

Benchmarks cover report parsing, extract, abundance/BIOM, MPA/Krona lineage export, Bracken filter/alpha, beta diversity, combine-kreports, and make-kreport.

## Library usage

```bash
cargo run --example basic_usage -- tests/fixtures/report_a.txt
```

Public modules include `krk_parser`, `lineage`, `bracken`, `diversity`, `combine`, `ktaxonomy`, `kreport_builder`, `logkrk_parser`, `sequence_processor`, `abundance_matrix`, and `biom`.

## Compatibility with KrakenTools

Implemented from public behavior of KrakenTools v1.2.1 (not by copying GPL code).

| KrakenTools | KrakenClip | Notes |
|-------------|------------|-------|
| `extract_kraken_reads.py` | `extract` | gzip, paired-end, `--max`; format preserved (not forced to FASTA) |
| `kreport2mpa.py` / `combine_mpa.py` | `abundance-matrix --format mpa` | Full lineages; multi-sample merge by path |
| `kreport2krona.py` | `abundance-matrix --format krona` | Direct reads + path; single report |
| `combine_kreports.py` | `combine-kreports` | Multi-sample columns + `--only-combined` |
| `filter_bracken_out.py` | `filter-bracken` | include/exclude + renormalized fractions |
| `alpha_diversity.py` | `alpha-diversity` | Shannon/BP/Simpson/InvSimpson/Fisher |
| `beta_diversity.py` | `beta-diversity` | Bray–Curtis for bracken/kreport/tsv |
| `make_kreport.py` | `make-kreport` | Needs condensed taxonomy input |
| `make_ktaxonomy.py` | — | Not implemented (offline DB tooling) |
| `fix_unmapped.py` | — | Not implemented |

## Testing

```bash
cargo test --locked
cargo clippy --locked --all-targets --all-features -- -D warnings
cargo bench --locked --no-run
```

Fixtures under `tests/fixtures/` cover reports, logs, FASTA/FASTQ, Bracken, and condensed taxonomy.

## License

MIT
