# KrakenClip

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Rust](https://img.shields.io/badge/rust-1.86%2B-blue.svg)](https://www.rust-lang.org/)
[![Docker Pulls](https://img.shields.io/docker/pulls/hecrp/krakenclip.svg)](https://hub.docker.com/r/hecrp/krakenclip)
[![Docker Image Size](https://img.shields.io/docker/image-size/hecrp/krakenclip.svg)](https://hub.docker.com/r/hecrp/krakenclip)
[![CI](https://github.com/hecrp/KrakenClip/actions/workflows/ci.yml/badge.svg)](https://github.com/hecrp/KrakenClip/actions/workflows/ci.yml)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/hecrp/krakenclip/graphs/commit-activity)

KrakenClip is a high-performance command-line utility written in [Rust](https://www.rust-lang.org/) for processing and analyzing [Kraken2](https://ccb.jhu.edu/software/kraken2/) bioinformatics software reports and log files. This toolkit focuses on fast and efficient processing of classifier outputs, making it an ideal choice for large datasets or bioinformatics pipelines as it works as a standalone binary with no external dependencies.

KrakenClip implements several functionalities inspired by [KrakenTools](https://github.com/jenniferlu717/KrakenTools), which was the pioneering suite of scripts for handling Kraken2 outputs. The merit for these functionalities belongs to KrakenTools, and if you use KrakenClip in your research, please consider citing KrakenTools as well:

Lu J, Rincon N, Wood DE, Breitwieser FP, Pockrandt C, Langmead B, Salzberg SL, Steinegger M. Metagenome analysis using the Kraken software suite. Nature Protocols, doi: 10.1038/s41596-022-00738-y (2022)

The inspiration to develop KrakenClip came from creating a high-performance version of KrakenTools while improving Rust programming skills. KrakenTools was extensively used during my PhD and integrated into metagenomic data analysis pipelines, making it an essential tool in my research workflow. This experience highlighted both the utility and the potential performance limitations of the original Python-based tools, motivating the development of this Rust implementation.

## Getting Started

### Prerequisites

- Rust programming language 1.86 or newer
- Cargo (Rust's package manager)

### Installation

1. Clone the repository:
   ```
   git clone https://github.com/hecrp/krakenclip.git
   cd krakenclip
   ```

2. Build the project:
   ```
   cargo build --release --locked
   ```

3. The executable will be available in `target/release/krakenclip`
    ```
    ./target/release/krakenclip --help
    ```

### Docker

You can also use the provided Dockerfile to build and run KrakenClip in a container:

1. Build the Docker image:
   ```
   docker build -t krakenclip .
   ```

2. Run the Docker container:
   ```
   docker run --rm krakenclip
   ```

3. To process files, mount volumes and run the toolkit:
   ```
   docker run --rm -v /path/to/local/data:/data krakenclip analyze /data/report.txt
   ```

4. Alternatively, pull the pre-built multi-architecture image from Docker Hub:
   ```
   docker pull hecrp/krakenclip:latest
   ```
   
   The image on Docker Hub is built with multi-architecture support, allowing it to run transparently on both x86 and ARM architectures, including Apple Silicon where it has been developed and tested. This means the same image will work natively on any system without manual configuration.

## Usage

The basic syntax for using KrakenClip is shown below:

```
USAGE:
    krakenclip [SUBCOMMAND]

OPTIONS:
    -h, --help       Print help information
    -V, --version    Print version information

SUBCOMMANDS:
    analyze               Analyze Kraken2 report
    extract               Extract sequences based on Kraken2 results
    abundance-matrix      Generate taxonomic abundance matrices from multiple reports
    generate-test-data    Generate test data for performance testing
    help                  Print this message or the help of the given subcommand(s)
```

### Analyze Module

Used to analyze and extract information from Kraken2 reports:

```
USAGE:
    krakenclip analyze [OPTIONS] <REPORT>

ARGS:
    <REPORT>    Kraken2 report file

OPTIONS:
    -h, --help                 Print help information
        --json <JSON>          Generate JSON output
        --tax-id <TAXON_ID>    Search for a specific taxon by ID
```

### Extract Module

Used to extract sequences based on Kraken2 results:

```
USAGE:
    krakenclip extract [OPTIONS] --output <OUTPUT> --taxids <TAXIDS> <SEQUENCE> <LOG>

ARGS:
    <SEQUENCE>                Input FASTA/FASTQ file
    <LOG>                     Kraken2 log file

OPTIONS:
    -h, --help                Print help information
    -o, --output <OUTPUT>     Output file for extracted sequences
        --report <REPORT>     Kraken2 report file (required for hierarchy (--include) options)
        --taxids <TAXIDS>     Comma-separated list of taxids to extract
        --include-children    Include sequences from all descendant taxa
        --include-parents     Include sequences from all ancestor taxa
        --exclude             Exclude sequences matching the specified taxids (inverse operation)
        --stats-output <FILE> Generate a statistics file with detailed extraction information
```

#### Hierarchical Extraction

The Extract module now supports hierarchical taxonomic extraction with two key options:

- **`--include-children`**: Extracts sequences from the specified taxids AND all their descendant taxa.

- **`--include-parents`**: Extracts sequences from the specified taxids AND all their ancestor taxa.

Both options require providing a Kraken2 report file with the `--report` option, as the taxonomic hierarchy information is extracted from there.

#### Statistics Report

The new `--stats-output` option generates a comprehensive markdown-formatted statistics file that includes:

- Total sequence counts (extracted vs. input)
- Breakdown of extracted sequences by taxid
- Percentage of sequences per taxid relative to total extracted and total input
- Distinction between original taxids and those added through hierarchical expansion (expanded)
- Summary statistics for original vs. expanded taxa

### Abundance Matrix Module

Used to generate taxonomic abundance matrices from Kraken2 reports:

```
USAGE:
    krakenclip abundance-matrix [OPTIONS] --output <o> <INPUT>...

ARGS:
    <INPUT>...               Input Kraken2 report files (can be multiple)

OPTIONS:
    -h, --help               Print help information
    -o, --output <o>         Output file for the abundance matrix
        --format <FORMAT>    Output format: tsv (default) or biom [default: tsv]
        --level <LEVEL>      Taxonomic level to aggregate abundances (S=species, G=genus, F=family,
                             O=order, C=class, P=phylum, D=domain; K is a D alias) [default: S]
        --min-abundance <MIN> Minimum abundance threshold [default: 0.0]
        --normalize          Normalize abundances to percentages during processing
        --include-unclassified Include unclassified sequences in the matrix
        --proportions        Transform counts to proportions (default behavior)
        --absolute-counts    Use absolute read counts without converting to proportions
```

#### Features
- Generates a matrix of taxonomic abundances across multiple samples
- Supports two output formats:
  - **TSV (default)**: Standard tab-separated values format
  - **BIOM**: Biological Observation Matrix format (v1.0.0), including all input samples in one row-major table
- Supports all standard Kraken2 taxonomic levels from species (`S`) to domain (`D`); `K` is retained as a compatibility alias for `D`
- Optional abundance threshold filtering. The threshold is a percentage for proportional output and a read count with `--absolute-counts`
- **Uses proportions (percentages) by default** for better comparability between samples
- Two options for handling abundance values:
  - **Proportions (default)**: Shows relative abundance as percentages
  - **Absolute counts**: Shows raw read counts (use `--absolute-counts` to enable)
- Complete handling of unclassified reads with `--include-unclassified`