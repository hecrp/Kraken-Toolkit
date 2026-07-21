use crate::abundance_matrix::{validate_taxonomic_level, AbundanceMatrix};
use crate::biom::{BiomMatrixType, BiomTable};
use crate::bracken::{filter_bracken, parse_bracken_file, write_bracken_file};
use crate::combine::{self, CombineOptions};
use crate::diversity::{
    alpha_from_bracken, beta_matrix, default_sample_name, load_sample_counts_bracken,
    load_sample_counts_kreport, load_sample_counts_tsv, write_beta_matrix, AlphaMetric,
};
use crate::generate_test_data;
use crate::kreport_builder;
use crate::krk_parser;
use crate::lineage::{self, LineageOptions};
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

#[derive(Parser)]
#[command(version = env!("CARGO_PKG_VERSION"), about = "A high-performance toolkit for processing Kraken2 reports, logs, and sequence files")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Analyze(AnalyzeArgs),
    Extract(ExtractArgs),
    #[command(name = "abundance-matrix")]
    AbundanceMatrix(AbundanceMatrixArgs),
    #[command(name = "combine-kreports")]
    CombineKreports(CombineKreportsArgs),
    #[command(name = "filter-bracken")]
    FilterBracken(FilterBrackenArgs),
    #[command(name = "alpha-diversity")]
    AlphaDiversity(AlphaDiversityArgs),
    #[command(name = "beta-diversity")]
    BetaDiversity(BetaDiversityArgs),
    #[command(name = "make-kreport")]
    MakeKreport(MakeKreportArgs),
    #[command(name = "generate-test-data")]
    GenerateTestData(GenerateTestDataArgs),
}

#[derive(Args)]
struct AnalyzeArgs {
    report: String,
    #[arg(long)]
    json: Option<String>,
    #[arg(long = "tax-id")]
    taxon_id: Option<u32>,
}

#[derive(Args)]
struct ExtractArgs {
    sequence: String,
    log: String,
    #[arg(short, long)]
    output: String,
    #[arg(long = "sequence2", visible_alias = "s2")]
    sequence2: Option<String>,
    #[arg(long = "output2", visible_alias = "o2")]
    output2: Option<String>,
    #[arg(long)]
    report: Option<String>,
    #[arg(long)]
    taxids: String,
    #[arg(long = "include-children")]
    include_children: bool,
    #[arg(long = "include-parents")]
    include_parents: bool,
    #[arg(long)]
    exclude: bool,
    #[arg(long = "stats-output")]
    stats_output: Option<String>,
    /// Stop after writing this many matching records
    #[arg(long = "max")]
    max_reads: Option<usize>,
}

#[derive(Args)]
struct AbundanceMatrixArgs {
    #[arg(required = true)]
    input: Vec<String>,
    #[arg(short, long)]
    output: String,
    #[arg(long, default_value = "tsv")]
    format: String,
    #[arg(long = "biom-matrix-type", default_value = "sparse")]
    biom_matrix_type: String,
    #[arg(long, default_value = "S")]
    level: String,
    #[arg(long = "min-abundance", default_value = "0.0")]
    min_abundance: f64,
    #[arg(long, conflicts_with = "absolute_counts")]
    normalize: bool,
    #[arg(long = "include-unclassified")]
    include_unclassified: bool,
    #[arg(long, conflicts_with = "absolute_counts")]
    proportions: bool,
    #[arg(long = "absolute-counts")]
    absolute_counts: bool,
    /// Include intermediate ranks in MPA/Krona lineages
    #[arg(long = "intermediate-ranks")]
    intermediate_ranks: bool,
    /// Use report percentages instead of clade reads for MPA
    #[arg(long = "percentages")]
    percentages: bool,
    /// Keep spaces in MPA taxon names
    #[arg(long = "keep-spaces")]
    keep_spaces: bool,
    /// Write an MPA header line
    #[arg(long = "display-header")]
    display_header: bool,
}

#[derive(Args)]
struct CombineKreportsArgs {
    #[arg(required = true)]
    input: Vec<String>,
    #[arg(short, long)]
    output: String,
    #[arg(long = "sample-names", value_delimiter = ',')]
    sample_names: Vec<String>,
    #[arg(long = "display-headers")]
    display_headers: bool,
    #[arg(long = "no-headers")]
    no_headers: bool,
    #[arg(long = "only-combined")]
    only_combined: bool,
}

#[derive(Args)]
struct FilterBrackenArgs {
    #[arg(short, long)]
    input: String,
    #[arg(short, long)]
    output: String,
    #[arg(long = "include", value_delimiter = ',')]
    include: Vec<u32>,
    #[arg(long = "exclude", value_delimiter = ',')]
    exclude: Vec<u32>,
}

#[derive(Args)]
struct AlphaDiversityArgs {
    #[arg(short, long)]
    input: String,
    /// Alpha metric: shannon, berger-parker, simpson, inverse-simpson, fisher
    #[arg(long = "type", default_value = "shannon")]
    metric: String,
}

#[derive(Args)]
struct BetaDiversityArgs {
    #[arg(short = 'i', long = "input", required = true, num_args = 1..)]
    input: Vec<String>,
    #[arg(short, long)]
    output: String,
    /// Input type: bracken, kreport, or tsv
    #[arg(long = "type", default_value = "bracken")]
    input_type: String,
    /// Taxonomic level filter for kreport inputs
    #[arg(long)]
    level: Option<String>,
    /// Category,count columns for generic TSV inputs
    #[arg(long = "cols", default_value = "0,1")]
    cols: String,
}

#[derive(Args)]
struct MakeKreportArgs {
    #[arg(short = 'k', long = "log")]
    log: String,
    #[arg(short = 't', long = "taxonomy")]
    taxonomy: String,
    #[arg(short, long)]
    output: String,
    #[arg(long = "use-read-len")]
    use_read_len: bool,
}

#[derive(Args)]
struct GenerateTestDataArgs {
    #[arg(short, long)]
    output: String,
    #[arg(short, long)]
    lines: usize,
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
        Commands::FilterBracken(args) => run_filter_bracken(args),
        Commands::AlphaDiversity(args) => run_alpha_diversity(args),
        Commands::BetaDiversity(args) => run_beta_diversity(args),
        Commands::MakeKreport(args) => run_make_kreport(args),
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
    println!(
        "Total time: {:.6} seconds",
        start_time.elapsed().as_secs_f64()
    );
    println!("File parsing time: {:.6} seconds", parse_time);
    println!(
        "Memory usage: {} bytes",
        after_memory.saturating_sub(before_memory)
    );
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
                let children = report.index.all_descendants(taxid);
                expanded_count += children.len();
                new_taxids.extend(children);
            }
            if args.include_parents {
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
    let mut readids = if want_stats {
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

    if let Some(max_reads) = args.max_reads {
        if readids.len() > max_reads {
            let mut limited = HashSet::with_capacity(max_reads);
            for id in readids.into_iter().take(max_reads) {
                limited.insert(id);
            }
            readids = limited;
        }
    }

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
    writeln!(writer, "# Kraken output: {}", args.log)?;
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

    let lineage_options = LineageOptions {
        intermediate_ranks: args.intermediate_ranks,
        use_percentages: args.percentages,
        replace_spaces: !args.keep_spaces,
        display_header: args.display_header,
    };

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
    let reports = parsed?;

    match args.format.as_str() {
        "mpa" => {
            lineage::write_combined_mpa(&reports, &args.output, lineage_options)?;
            println!(
                "MPA format output successfully generated in: {}",
                args.output
            );
        }
        "krona" => {
            if reports.len() != 1 {
                return Err(
                    "Krona export currently supports a single input report (KrakenTools kreport2krona semantics)"
                        .into(),
                );
            }
            lineage::write_krona_file(&reports[0].1, &args.output, lineage_options)?;
            println!(
                "Krona format output successfully generated in: {}",
                args.output
            );
        }
        "tsv" | "biom" => {
            let biom_matrix_type = match args.biom_matrix_type.as_str() {
                "dense" => BiomMatrixType::Dense,
                "sparse" => BiomMatrixType::Sparse,
                other => {
                    return Err(format!(
                        "Unsupported BIOM matrix type '{other}'. Use 'dense' or 'sparse'."
                    )
                    .into())
                }
            };
            let mut matrix = AbundanceMatrix::new(&args.level);
            matrix.set_force_include_unclassified(args.include_unclassified);
            let proportional = args.normalize || args.proportions || !args.absolute_counts;
            for (sample_name, report) in &reports {
                matrix.add_sample(report, sample_name, args.min_abundance, proportional);
            }
            if args.format == "tsv" {
                matrix.write_matrix(&args.output)?;
                println!(
                    "Abundance matrix successfully generated in: {}",
                    args.output
                );
            } else {
                BiomTable::from_abundance_matrix(&matrix, biom_matrix_type)
                    .write_json(&args.output)?;
                println!(
                    "BIOM format output successfully generated in: {}",
                    args.output
                );
            }
        }
        _ => unreachable!(),
    }
    Ok(())
}

fn run_combine_kreports(args: CombineKreportsArgs) -> Result<(), Box<dyn Error>> {
    let options = CombineOptions {
        sample_names: args.sample_names,
        display_headers: args.display_headers,
        no_headers: args.no_headers,
        only_combined: args.only_combined,
    };
    combine::combine_kreports(&args.input, &args.output, &options)?;
    println!("Combined report written to {}", args.output);
    Ok(())
}

fn run_filter_bracken(args: FilterBrackenArgs) -> Result<(), Box<dyn Error>> {
    if !args.include.is_empty() && !args.exclude.is_empty() {
        return Err("--include and --exclude are mutually exclusive".into());
    }
    let records = parse_bracken_file(&args.input)?;
    let include = if args.include.is_empty() {
        None
    } else {
        Some(args.include.into_iter().collect::<HashSet<_>>())
    };
    let exclude = if args.exclude.is_empty() {
        None
    } else {
        Some(args.exclude.into_iter().collect::<HashSet<_>>())
    };
    let filtered = filter_bracken(&records, include.as_ref(), exclude.as_ref());
    write_bracken_file(&args.output, &filtered)?;
    println!(
        "Wrote {} Bracken records to {}",
        filtered.len(),
        args.output
    );
    Ok(())
}

fn run_alpha_diversity(args: AlphaDiversityArgs) -> Result<(), Box<dyn Error>> {
    let metric = AlphaMetric::parse(&args.metric)
        .ok_or_else(|| format!("unsupported alpha metric '{}'", args.metric))?;
    let records = parse_bracken_file(&args.input)?;
    let value = alpha_from_bracken(&records, metric);
    println!("{value:.6}");
    Ok(())
}

fn run_beta_diversity(args: BetaDiversityArgs) -> Result<(), Box<dyn Error>> {
    let samples: Result<Vec<_>, Box<dyn Error + Send + Sync>> = args
        .input
        .par_iter()
        .enumerate()
        .map(|(index, path)| {
            let name = default_sample_name(path, index);
            match args.input_type.as_str() {
                "bracken" => load_sample_counts_bracken(path, &name)
                    .map_err(|e| -> Box<dyn Error + Send + Sync> { e.to_string().into() }),
                "kreport" => load_sample_counts_kreport(path, &name, args.level.as_deref())
                    .map_err(|e| -> Box<dyn Error + Send + Sync> { e.to_string().into() }),
                "tsv" => {
                    let mut parts = args.cols.split(',');
                    let category = parts
                        .next()
                        .ok_or("cols must be category,count")?
                        .parse::<usize>()
                        .map_err(|_| "invalid category column")?;
                    let count = parts
                        .next()
                        .ok_or("cols must be category,count")?
                        .parse::<usize>()
                        .map_err(|_| "invalid count column")?;
                    load_sample_counts_tsv(path, &name, category, count)
                        .map_err(|e| -> Box<dyn Error + Send + Sync> { e.to_string().into() })
                }
                other => Err(format!("unsupported beta input type '{other}'").into()),
            }
        })
        .collect();
    let samples = samples.map_err(|e| e.to_string())?;
    let matrix = beta_matrix(&samples);
    write_beta_matrix(&samples, &matrix, &args.output)?;
    println!("Beta diversity matrix written to {}", args.output);
    Ok(())
}

fn run_make_kreport(args: MakeKreportArgs) -> Result<(), Box<dyn Error>> {
    kreport_builder::make_kreport(&args.log, &args.taxonomy, &args.output, args.use_read_len)?;
    println!("Kraken report written to {}", args.output);
    Ok(())
}

fn run_generate_test_data(args: GenerateTestDataArgs) -> Result<(), Box<dyn Error>> {
    generate_test_data::generate_data(&args.output, args.lines, &args.r#type)?;
    println!("Test data generated successfully");
    Ok(())
}
