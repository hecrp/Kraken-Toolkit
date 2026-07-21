# KrakenClip

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Rust](https://img.shields.io/badge/rust-1.86%2B-blue.svg)](https://www.rust-lang.org/)
[![Docker Pulls](https://img.shields.io/docker/pulls/hecrp/krakenclip.svg)](https://hub.docker.com/r/hecrp/krakenclip)
[![Docker Image Size](https://img.shields.io/docker/image-size/hecrp/krakenclip.svg)](https://hub.docker.com/r/hecrp/krakenclip)
[![CI](https://github.com/hecrp/KrakenClip/actions/workflows/ci.yml/badge.svg)](https://github.com/hecrp/KrakenClip/actions/workflows/ci.yml)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/hecrp/krakenclip/graphs/commit-activity)

KrakenClip is a high-performance command-line toolkit written in [Rust](https://www.rust-lang.org/) for processing [Kraken2](https://ccb.jhu.edu/software/kraken2/) reports, classification logs, and sequence files. It is designed as a fast, dependency-free alternative for the most common post-processing workflows used in metagenomic pipelines.

KrakenClip implements functionality inspired by [KrakenTools](https://github.com/jenniferlu717/KrakenTools). Credit for the original ideas and workflows belongs to KrakenTools; if you use KrakenClip in your research, please also cite:

> Lu J, Rincon N, Wood DE, Breitwieser FP, Pockrandt C, Langmead B, Salzberg SL, Steinegger M. Metagenome analysis using the Kraken software suite. Nature Protocols, doi: 10.1038/s41596-022-00738-y (2022)

## Features

- **Standalone binary** with no runtime Python or Biopython dependencies
- **Fast report parsing** with indexed taxonomy for parent/child expansion
- **Sequence extraction** from FASTA/FASTQ, including gzip and paired-end inputs
- **Multi-sample abundance matrices** in TSV, BIOM, MPA, and Krona formats
- **Report combination** across samples (`combine-kreports`)
- **Built-in performance tooling** via Criterion benches and synthetic test-data generation
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

Build locally:

```bash
docker build -t krakenclip .
docker run --rm krakenclip --help
```

Process local files by mounting a data directory:

```bash
docker run --rm -v /path/to/local/data:/data krakenclip analyze /data/report.txt
```

Or pull the multi-architecture image from Docker Hub:

```bash
docker pull hecrp/krakenclip:latest
```

The published image runs natively on both x86_64 and ARM architectures, including Apple Silicon.

## Quick examples

```bash
# Inspect a Kraken2 report and optionally emit JSON
krakenclip analyze sample.kreport --tax-id 562 --json sample.json

# Extract reads assigned to one or more taxids
krakenclip extract reads.fastq kraken.log \
  --output selected.fastq \
  --taxids 562,1280

# Extract a taxon plus all descendants
krakenclip extract reads.fastq kraken.log \
  --output bacteria.fastq \
  --taxids 2 \
  --include-children \
  --report sample.kreport

# Paired-end + gzip extraction
krakenclip extract reads_R1.fastq.gz kraken.log \
  --sequence2 reads_R2.fastq.gz \
  --output selected_R1.fastq.gz \
  --output2 selected_R2.fastq.gz \
  --taxids 562

# Build a multi-sample abundance matrix
krakenclip abundance-matrix sample_a.kreport sample_b.kreport \
  --output matrix.tsv \
  --level S

# Export BIOM (sparse by default), MPA, or Krona
krakenclip abundance-matrix *.kreport --output matrix.biom --format biom
krakenclip abundance-matrix *.kreport --output matrix.mpa --format mpa
krakenclip abundance-matrix sample.kreport --output krona.txt --format krona

# Combine multiple Kraken2 reports
krakenclip combine-kreports sample_*.kreport --output combined.kreport
```

## Usage

```text
Usage: krakenclip <COMMAND>

Commands:
  analyze             Analyzes a Kraken2 report
  extract             Extracts sequences based on Kraken2 results
  abundance-matrix    Generates taxonomic abundance matrices from multiple reports
  combine-kreports    Combines multiple Kraken2 reports into one
  generate-test-data  Generates test data for performance testing
  help                Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

### `analyze`

Inspect a Kraken2 report, query a taxon, and optionally export JSON.

```text
Usage: krakenclip analyze [OPTIONS] <REPORT>

Arguments:
  <REPORT>  Kraken2 report file

Options:
      --json <JSON>        Generate JSON output
      --tax-id <TAXON_ID>  Look for a specific taxon by ID
  -h, --help               Print help
```

Example:

```bash
krakenclip analyze sample.kreport --tax-id 562 --json sample.json
```

The command prints parse time and an approximate memory delta after processing.

### `extract`

Extract or exclude reads from FASTA/FASTQ files using a Kraken2 classification log.

```text
Usage: krakenclip extract [OPTIONS] --output <OUTPUT> --taxids <TAXIDS> <SEQUENCE> <LOG>

Arguments:
  <SEQUENCE>  Input FASTA/FASTQ file (optionally .gz)
  <LOG>       Kraken2 log file (optionally .gz)

Options:
  -o, --output <OUTPUT>              Output file for extracted sequences (optionally .gz)
      --sequence2 <SEQUENCE2>        Mate/pair FASTA/FASTQ for paired-end extraction [alias: --s2]
      --output2 <OUTPUT2>            Mate/pair output file (required with --sequence2) [alias: --o2]
      --report <REPORT>              Kraken2 report file (required for hierarchy options)
      --taxids <TAXIDS>              Comma-separated list of taxids to extract
      --include-children             Include sequences from all descendant taxa
      --include-parents              Include sequences from all ancestor taxa
      --exclude                      Exclude sequences matching the specified taxids
      --stats-output <STATS_OUTPUT>  Generate a statistics file with detailed information
  -h, --help                         Print help
```

#### Hierarchical extraction

- `--include-children` keeps the requested taxids and all descendants
- `--include-parents` keeps the requested taxids and all ancestors
- Both options require `--report`
- Expansion uses a taxonomy index built while parsing the report

#### Paired-end and gzip

- Provide `--sequence2/--s2` together with `--output2/--o2`
- Mate IDs may include `/1` and `/2` suffixes; matching is done on the shared read ID
- `.gz` inputs/outputs are detected from the file extension or gzip magic bytes

#### Statistics

`--stats-output` writes a CSV file with:

- total sequences seen during the extraction pass
- sequences extracted per taxid
- percentages relative to extracted and total input
- whether each taxid was requested directly or added by hierarchy expansion

Sequence totals are counted while streaming the input, so extraction does not require a second full-file scan.

### `abundance-matrix`

Build taxonomic abundance tables from one or more Kraken2 reports. Reports are parsed in parallel.

```text
Usage: krakenclip abundance-matrix [OPTIONS] --output <OUTPUT> <INPUT>...

Arguments:
  <INPUT>...  Input Kraken2 report files (can be multiple)

Options:
  -o, --output <OUTPUT>
          Output file for the abundance matrix
      --format <FORMAT>
          Output format (tsv, biom, mpa, or krona) [default: tsv]
      --biom-matrix-type <BIOM_MATRIX_TYPE>
          BIOM matrix encoding when --format biom (dense or sparse) [default: sparse]
      --level <LEVEL>
          Taxonomic level for aggregating abundances [default: S]
      --min-abundance <MIN_ABUNDANCE>
          Minimum abundance threshold (0.0-100.0) [default: 0.0]
      --normalize
          Normalize abundances to percentages during processing
      --include-unclassified
          Include unclassified sequences in the matrix
      --proportions
          Transform counts to proportions
      --absolute-counts
          Use absolute read counts without converting to proportions
  -h, --help
          Print help
```

#### Output formats

| Format | Description |
|--------|-------------|
| `tsv` | Tab-separated taxon × sample matrix (default) |
| `biom` | BIOM 1.0.0 JSON table; sparse by default, dense optional |
| `mpa` | MetaPhlAn-style abundance table |
| `krona` | Simple Krona text input (`magnitude` + taxon name) |

#### Abundance units

- Proportions/percentages are the default
- Use `--absolute-counts` for raw clade read counts
- `--min-abundance` is interpreted in the same units as the output
- Supported levels: `S`, `G`, `F`, `O`, `C`, `P`, `D` (`K` is accepted as an alias for `D`)

Examples:

```bash
krakenclip abundance-matrix *.kreport -o matrix.tsv --level G
krakenclip abundance-matrix *.kreport -o matrix.biom --format biom --biom-matrix-type sparse
krakenclip abundance-matrix sample.kreport -o matrix.tsv --absolute-counts --min-abundance 100
```

### `combine-kreports`

Combine multiple Kraken2 reports by summing clade and direct reads for shared taxids.

```text
Usage: krakenclip combine-kreports --output <OUTPUT> <INPUT>...

Arguments:
  <INPUT>...  Input Kraken2 report files

Options:
  -o, --output <OUTPUT>  Combined output report path
  -h, --help             Print help
```

Example:

```bash
krakenclip combine-kreports sample_a.kreport sample_b.kreport -o combined.kreport
```

### `generate-test-data`

Generate synthetic Kraken2-style reports for benchmarking and stress testing.

```text
Usage: krakenclip generate-test-data --output <OUTPUT> --lines <LINES> --type <TYPE>

Options:
  -o, --output <OUTPUT>  Output file path
  -l, --lines <LINES>    Number of lines to generate
  -t, --type <TYPE>      Type of data to generate (wide, deep, fragments, dense, etc.)
  -h, --help             Print help
```

Useful `--type` values include `wide`, `deep`, `dense`, `fragments`, `extreme`, `unbalanced`, and `mixed`.

Example:

```bash
krakenclip generate-test-data -o big.kreport -l 100000 -t dense
```

## Performance

KrakenClip focuses on the hot paths that dominate Kraken2 post-processing:

- numeric taxid filtering in log parsing
- single-pass sequence extraction with optional statistics
- indexed taxonomy for hierarchical expansion
- parallel multi-sample report parsing
- compact JSON serialization for report/BIOM exports

Run the Criterion suites locally:

```bash
cargo bench --bench parsing_benchmark -- --quick
cargo bench --bench pipeline_benchmark -- --quick
```

Example `--quick` results from one development host (illustrative only; re-measure on your hardware):

| Benchmark | Time |
|-----------|------|
| parse kraken2 report fixture | ~12 µs |
| parse kraken2 report ~100k lines | ~11.5 ms |
| parse kraken log 50k reads | ~2.1 ms |
| extract matching FASTQ 50k reads | ~10.4 ms |
| abundance-matrix 8×20k reports | ~240 ms |
| BIOM sparse serialize | ~445 ms |

Do not treat fixture-only timings as representative of production metagenomic workloads, and do not claim speedups versus KrakenTools without a side-by-side comparison on the same machine and inputs.

## Library usage

KrakenClip can also be used as a Rust library:

```bash
cargo run --example basic_usage -- tests/fixtures/report_a.txt
```

The example shows how to:

1. parse a Kraken2 report
2. query a taxon
3. inspect the taxonomy index
4. build an in-memory abundance matrix

Public modules include `krk_parser`, `logkrk_parser`, `sequence_processor`, `abundance_matrix`, `biom`, `combine`, and `taxon_query`.

## Compatibility with KrakenTools

KrakenClip is a focused, high-performance subset of KrakenTools workflows rather than a full 1:1 replacement.

| KrakenTools script | KrakenClip equivalent | Notes |
|--------------------|-----------------------|-------|
| `extract_kraken_reads.py` | `extract` | Supports paired-end and gzip |
| `kreport2mpa.py` / `combine_mpa.py` | `abundance-matrix --format mpa` | Multi-sample MPA-style table |
| `kreport2krona.py` | `abundance-matrix --format krona` | Simple Krona text output |
| `combine_kreports.py` | `combine-kreports` | Sums shared taxids across reports |
| Bracken filters | — | Not implemented |
| Alpha/beta diversity scripts | — | Not implemented |
| Taxonomy DB builders (`make_ktaxonomy.py`, `make_kreport.py`) | — | Not implemented |

## Testing

```bash
cargo test --locked
cargo clippy --locked --all-targets --all-features -- -D warnings
cargo bench --locked --no-run
```

Integration fixtures live in `tests/fixtures/` and cover report analysis, extraction, exclusion, hierarchy expansion, gzip/paired-end extraction, abundance thresholds, BIOM sparse/dense export, MPA output, and report combination.

## License

MIT
