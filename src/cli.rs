use crate::abundance_matrix::{validate_taxonomic_level, AbundanceMatrix};
use crate::biom::{BiomMatrixType, BiomTable};
use crate::combine;
use crate::generate_test_data;
use crate::krk_parser;
use crate::logkrk_parser;
use crate::sequence_processor;
use crate::taxon_query::{find_taxon_info, print_taxon_info};
use chrono;
use clap::{Args, Parser, Subcommand};
use memory_stats::memory_stats;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

const BUFFER_SIZE: usize = 512 * 1024;

/// KrakenClip - High-performance Kraken2 data processing toolkit
#[derive(Parser)]
#[command(version = env!("CARGO_PKG_VERSION"), about = "A high-performance toolkit for processing Kraken2 reports, logs, and sequence files")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Analyzes a Kraken2 report
    Analyze(AnalyzeArgs),

    /// Extracts sequences based on Kraken2 results
    Extract(ExtractArgs),

    /// Generates taxonomic abundance matrices from multiple reports
    #[command(name = "abundance-matrix")]
    AbundanceMatrix(AbundanceMatrixArgs),

    /// Combines multiple Kraken2 reports into one
    #[command(name = "combine-kreports")]
    CombineKreports(CombineKreportsArgs),

    /// Generates test data for performance testing
    #[command(name = "generate-test-data")]
    GenerateTestData(GenerateTestDataArgs),
}

#[derive(Args)]
struct AnalyzeArgs {
    /// Kraken2 report file
    report: String,

    /// Generate JSON output
    #[arg(long)]
    json: Option<String>,

    /// Look for a specific taxon by ID
    #[arg(long = "tax-id")]
    taxon_id: Option<u32>,
}

#[derive(Args)]
struct ExtractArgs {
    /// Input FASTA/FASTQ file (optionally .gz)
    sequence: String,

    /// Kraken2 log file (optionally .gz)
    log: String,

    /// Output file for extracted sequences (optionally .gz)
    #[arg(short, long)]
    output: String,

    /// Optional mate/pair FASTA/FASTQ file for paired-end extraction
    #[arg(long = "sequence2", visible_alias = "s2")]
    sequence2: Option<String>,

    /// Output file for mate/pair sequences (required with --sequence2)
    #[arg(long = "output2", visible_alias = "o2")]
    output2: Option<String>,

    /// Kraken2 report file (required for hierarchy options)
    #[arg(long)]
    report: Option<String>,

    /// Comma-separated list of taxids to extract
    #[arg(long)]
    taxids: String,

    /// Include sequences from all descendant taxa
    #[arg(long = "include-children")]
    include_children: bool,

    /// Include sequences from all ancestor taxa
    #[arg(long = "include-parents")]
    include_parents: bool,

    /// Exclude sequences matching the specified taxids
    #[arg(long)]
    exclude: bool,

    /// Generate a statistics file with detailed information
    #[arg(long = "stats-output")]
    stats_output: Option<String>,
}

#[derive(Args)]
struct AbundanceMatrixArgs {
    /// Input Kraken2 report files (can be multiple)
    #[arg(required = true)]
    input: Vec<String>,

    /// Output file for the abundance matrix
    #[arg(short, long)]
    output: String,

    /// Output format (tsv, biom, mpa, or krona)
    #[arg(long, default_value = "tsv")]
    format: String,

    /// BIOM matrix encoding when --format biom (dense or sparse)
    #[arg(long = "biom-matrix-type", default_value = "sparse")]
    biom_matrix_type: String,

    /// Taxonomic level for aggregating abundances
    #[arg(long, default_value = "S")]
    level: String,

    /// Minimum abundance threshold (0.0-100.0)
    #[arg(long = "min-abundance", default_value = "0.0")]
    min_abundance: f64,

    /// Normalize abundances to percentages during processing
    #[arg(long, conflicts_with = "absolute_counts")]
    normalize: bool,

    /// Include unclassified sequences in the matrix
    #[arg(long = "include-unclassified")]
    include_unclassified: bool,

    /// Transform counts to proportions
    #[arg(long, conflicts_with = "absolute_counts")]
    proportions: bool,

    /// Use absolute read counts without converting to proportions
    #[arg(long = "absolute-counts")]
    absolute_counts: bool,
}

#[derive(Args)]
struct CombineKreportsArgs {
    /// Input Kraken2 report files
    #[arg(required = true)]
    input: Vec<String>,

    /// Combined output report path
    #[arg(short, long)]
    output: String,
}

#[derive(Args)]
struct GenerateTestDataArgs {
    /// Output file path
    #[arg(short, long)]
    output: String,

    /// Number of lines to generate
    #[arg(short, long)]
    lines: usize,

    /// Type of data to generate (wide, deep, fragments, dense, etc.)
    #[arg(short, long)]
    r#type: String,
}

pub fn run_cli() {
    let cli = Cli::parse();

    let result = match cli.command {
        Commands::Analyze(args) => run_analyze(args),
        Commands::Extract(args) => run_extract(args),
        Commands::AbundanceMatrix(args) => run_abundance_matrix(args),
        Commands::CombineKreports(args) => run_combine_kreports(args),
        Commands::GenerateTestData(args) => run_generate_test_data(args),
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

fn run_analyze(args: AnalyzeArgs) -> Result<(), Box<dyn Error>> {
    let before_memory = memory_stats().map(|s| s.physical_mem).unwrap_or(0);
    let start_time = Instant::now();

    let (report, parse_time) = krk_parser::parse_kraken2_report(&args.report)
        .map_err(|e| format!("Error parsing Kraken2 report file '{}': {}", args.report, e))?;

    if let Some(taxon_id) = args.taxon_id {
        match find_taxon_info(&report.root, taxon_id as u64) {
            Some(taxon_info) => print_taxon_info(&taxon_info),
            None => println!("Taxon with ID {} not found", taxon_id),
        }
    }

    if let Some(json_output) = args.json {
        krk_parser::write_json_report(&report, &json_output)
            .map_err(|e| format!("Error writing JSON report '{}': {}", json_output, e))?;
        println!("JSON report written to {}", json_output);
    }

    let after_memory = memory_stats().map(|s| s.physical_mem).unwrap_or(0);
    let memory_used = after_memory.saturating_sub(before_memory);
    let total_time = start_time.elapsed().as_secs_f64();

    println!("Total time: {:.6} seconds", total_time);
    println!("File parsing time: {:.6} seconds", parse_time);
    println!("Memory usage: {} bytes", memory_used);

    Ok(())
}

fn run_extract(args: ExtractArgs) -> Result<(), Box<dyn Error>> {
    if args.sequence2.is_some() != args.output2.is_some() {
        return Err(
            "Error: --sequence2/--s2 and --output2/--o2 must be provided together for paired-end extraction."
                .into(),
        );
    }

    let taxids: HashSet<u32> = args
        .taxids
        .split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| {
            s.parse::<u32>()
                .map_err(|_| format!("Invalid taxid '{s}': expected an unsigned integer"))
        })
        .collect::<Result<_, _>>()?;

    if taxids.is_empty() {
        return Err("Error: at least one taxid is required".into());
    }

    let original_taxids = taxids.clone();
    let mut expanded_taxids = taxids.clone();

    if (args.include_children || args.include_parents) && args.report.is_some() {
        let report_file = args.report.clone().unwrap();
        println!("Reading taxonomy from report file: {}", report_file);

        let (report, _) = krk_parser::parse_kraken2_report(&report_file)
            .map_err(|e| format!("Error parsing Kraken2 report file '{}': {}", report_file, e))?;

        let mut expanded_count = 0;
        let mut new_taxids = HashSet::with_capacity(taxids.len() * 10);
        new_taxids.extend(taxids.iter().copied());

        for &taxid in &taxids {
            if args.include_children {
                println!("Including children for taxid: {}", taxid);
                let children = report.index.all_descendants(taxid);
                expanded_count += children.len();
                new_taxids.extend(children);
            }

            if args.include_parents {
                println!("Including parents for taxid: {}", taxid);
                let parents = report.index.all_ancestors(taxid);
                expanded_count += parents.len();
                new_taxids.extend(parents);
            }
        }

        expanded_taxids = new_taxids;
        println!(
            "Expanded to {} taxids (added {} through hierarchy)",
            expanded_taxids.len(),
            expanded_count
        );
    } else if (args.include_children || args.include_parents) && args.report.is_none() {
        return Err("Error: A report file (--report) is required when using --include-children or --include-parents options.".into());
    }

    let want_stats = args.stats_output.is_some();
    let mut taxid_readid_map: HashMap<u32, HashSet<String>> = HashMap::new();

    let readids = if want_stats {
        logkrk_parser::parse_kraken_output_with_taxids(
            &args.log,
            &expanded_taxids,
            &mut taxid_readid_map,
        )
        .map_err(|e| format!("Error parsing Kraken2 log file: {}", e))?
    } else {
        logkrk_parser::parse_kraken_output(&args.log, &expanded_taxids)
            .map_err(|e| format!("Error parsing Kraken2 log file: {}", e))?
    };

    let stats = if let (Some(sequence2), Some(output2)) = (&args.sequence2, &args.output2) {
        sequence_processor::process_paired_sequence_files(
            &args.sequence,
            sequence2,
            &readids,
            &args.output,
            output2,
            args.exclude,
        )
        .map_err(|e| std::io::Error::other(e.to_string()))?
    } else {
        sequence_processor::process_sequence_files(
            std::slice::from_ref(&args.sequence),
            &readids,
            &args.output,
            args.exclude,
        )
        .map_err(|e| std::io::Error::other(e.to_string()))?
    };

    println!("Sequences extracted successfully to {}", args.output);
    if let Some(output2) = &args.output2 {
        println!("Paired sequences written to {}", output2);
    }
    println!(
        "{} sequences matching {} taxids",
        stats.written_sequences,
        expanded_taxids.len()
    );
    println!(
        "Total {} sequences: {}",
        if args.exclude {
            "excluded"
        } else {
            "extracted"
        },
        stats.written_sequences
    );

    if let Some(ref stats_file) = args.stats_output {
        generate_statistics_file(
            stats_file,
            &taxid_readid_map,
            &original_taxids,
            stats.total_sequences,
            &args,
        )?;
        println!("Statistics written to {}", stats_file);
    }

    Ok(())
}

fn generate_statistics_file(
    stats_file: &str,
    taxid_readid_map: &HashMap<u32, HashSet<String>>,
    original_taxids: &HashSet<u32>,
    total_sequences: usize,
    args: &ExtractArgs,
) -> Result<(), Box<dyn Error>> {
    let file = File::create(stats_file)?;
    let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);

    let mut total_extracted = 0;
    let mut total_original_taxids = 0;
    let mut total_expanded_taxids = 0;
    let mut orig_sequences = 0;
    let mut expanded_sequences = 0;
    let mut stats: Vec<(u32, usize, bool)> = Vec::with_capacity(taxid_readid_map.len());

    for (taxid, readids) in taxid_readid_map {
        let count = readids.len();
        total_extracted += count;
        let is_original = original_taxids.contains(taxid);
        if is_original {
            total_original_taxids += 1;
            orig_sequences += count;
        } else {
            total_expanded_taxids += 1;
            expanded_sequences += count;
        }
        stats.push((*taxid, count, is_original));
    }

    stats.sort_by_key(|entry| std::cmp::Reverse(entry.1));

    let percent_extracted = if total_sequences > 0 {
        (total_extracted as f64 / total_sequences as f64) * 100.0
    } else {
        0.0
    };

    writeln!(writer, "# KrakenClip Extraction Statistics")?;
    writeln!(
        writer,
        "# Date: {}",
        chrono::Local::now().format("%Y-%m-%d %H:%M:%S")
    )?;
    writeln!(writer, "# Input file: {}", args.sequence)?;
    if let Some(ref sequence2) = args.sequence2 {
        writeln!(writer, "# Paired input file: {}", sequence2)?;
    }
    writeln!(writer, "# Kraken output: {}", args.log)?;
    if let Some(ref report) = args.report {
        writeln!(writer, "# Kraken report: {}", report)?;
    }
    writeln!(writer, "# Include children: {}", args.include_children)?;
    writeln!(writer, "# Include parents: {}", args.include_parents)?;
    writeln!(writer, "# Exclude mode: {}", args.exclude)?;
    writeln!(writer, "# Total sequences in input: {}", total_sequences)?;
    writeln!(
        writer,
        "# Total sequences extracted: {} ({:.2}%)",
        total_extracted, percent_extracted
    )?;
    writeln!(
        writer,
        "# Original taxids: {} (sequences: {})",
        total_original_taxids, orig_sequences
    )?;
    writeln!(
        writer,
        "# Expanded taxids: {} (sequences: {})",
        total_expanded_taxids, expanded_sequences
    )?;
    writeln!(writer)?;
    writeln!(
        writer,
        "taxid,sequences,percent_of_extracted,percent_of_total,is_original"
    )?;

    for (taxid, count, is_original) in stats {
        let percent_of_extracted = if total_extracted > 0 {
            (count as f64 / total_extracted as f64) * 100.0
        } else {
            0.0
        };
        let percent_of_total = if total_sequences > 0 {
            (count as f64 / total_sequences as f64) * 100.0
        } else {
            0.0
        };
        writeln!(
            writer,
            "{},{},{:.4},{:.4},{}",
            taxid, count, percent_of_extracted, percent_of_total, is_original
        )?;
    }

    writer.flush()?;
    Ok(())
}

fn run_abundance_matrix(args: AbundanceMatrixArgs) -> Result<(), Box<dyn Error>> {
    if !validate_taxonomic_level(&args.level) {
        return Err(format!(
            "The taxonomic level '{}' is not valid. Use D/K, P, C, O, F, G or S.",
            args.level
        )
        .into());
    }
    if args.min_abundance < 0.0 {
        return Err("Minimum abundance cannot be negative".into());
    }

    let supported = ["tsv", "biom", "mpa", "krona"];
    if !supported.contains(&args.format.as_str()) {
        return Err(format!(
            "Unsupported output format '{}'. Use 'tsv', 'biom', 'mpa', or 'krona'.",
            args.format
        )
        .into());
    }

    let biom_matrix_type = match args.biom_matrix_type.as_str() {
        "dense" => BiomMatrixType::Dense,
        "sparse" => BiomMatrixType::Sparse,
        other => {
            return Err(
                format!("Unsupported BIOM matrix type '{other}'. Use 'dense' or 'sparse'.").into(),
            )
        }
    };

    let mut matrix = AbundanceMatrix::new(&args.level);
    matrix.set_force_include_unclassified(args.include_unclassified);
    let proportional = args.normalize || args.proportions || !args.absolute_counts;

    let mut used_sample_names = HashSet::new();
    let jobs: Vec<(String, String)> = args
        .input
        .iter()
        .enumerate()
        .map(|(index, file)| {
            let default_name = format!("sample_{}", index + 1);
            let base_name = Path::new(file)
                .file_stem()
                .and_then(|name| name.to_str())
                .unwrap_or(&default_name);
            let sample_name = if used_sample_names.insert(base_name.to_string()) {
                base_name.to_string()
            } else {
                format!("{base_name}_{}", index + 1)
            };
            (file.clone(), sample_name)
        })
        .collect();

    let parsed: Result<Vec<_>, String> = jobs
        .par_iter()
        .map(|(file, sample_name)| {
            println!("Processing sample: {sample_name}");
            let (report, _) = krk_parser::parse_kraken2_report(file)
                .map_err(|error| format!("Error parsing Kraken2 report file '{file}': {error}"))?;
            Ok((sample_name.clone(), report))
        })
        .collect();

    for (sample_name, report) in parsed? {
        matrix.add_sample(&report, &sample_name, args.min_abundance, proportional);
    }

    match args.format.as_str() {
        "tsv" => {
            matrix
                .write_matrix(&args.output)
                .map_err(|error| format!("Error generating abundance matrix: {error}"))?;
            println!(
                "Abundance matrix successfully generated in: {}",
                args.output
            );
        }
        "biom" => {
            BiomTable::from_abundance_matrix(&matrix, biom_matrix_type)
                .write_json(&args.output)
                .map_err(|error| format!("Error generating BIOM output: {error}"))?;
            println!(
                "BIOM format output successfully generated in: {}",
                args.output
            );
        }
        "mpa" => {
            matrix
                .write_mpa(&args.output)
                .map_err(|error| format!("Error generating MPA output: {error}"))?;
            println!(
                "MPA format output successfully generated in: {}",
                args.output
            );
        }
        "krona" => {
            matrix
                .write_krona(&args.output)
                .map_err(|error| format!("Error generating Krona output: {error}"))?;
            println!(
                "Krona format output successfully generated in: {}",
                args.output
            );
        }
        _ => unreachable!("format was validated above"),
    }

    Ok(())
}

fn run_combine_kreports(args: CombineKreportsArgs) -> Result<(), Box<dyn Error>> {
    combine::combine_kreports(&args.input, &args.output)?;
    println!("Combined report written to {}", args.output);
    Ok(())
}

fn run_generate_test_data(args: GenerateTestDataArgs) -> Result<(), Box<dyn Error>> {
    match generate_test_data::generate_data(&args.output, args.lines, &args.r#type) {
        Ok(_) => println!("Test data generated successfully"),
        Err(e) => return Err(format!("Error generating test data: {}", e).into()),
    }
    Ok(())
}
