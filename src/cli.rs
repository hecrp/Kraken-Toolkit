use crate::abundance_matrix::{validate_taxonomic_level, AbundanceMatrix};
use crate::biom::BiomTable;
use crate::generate_test_data;
use crate::krk_parser;
use crate::logkrk_parser;
use crate::sequence_processor;
use crate::taxon_query::{find_taxon_info, print_taxon_info};
use chrono;
use clap::{Args, Parser, Subcommand};
use memory_stats::memory_stats;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::time::Instant;

const BUFFER_SIZE: usize = 512 * 1024; // 512KB buffer for I/O

/// KrakenClip - High-performance Kraken2 data processing toolkit
#[derive(Parser)]
#[command(version = env!("CARGO_PKG_VERSION"), about = "A high-performance toolkit for processing Kraken2 reports, logs, and sequence files")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

/// Available subcommands
#[derive(Subcommand)]
enum Commands {
    /// Analyzes a Kraken2 report
    Analyze(AnalyzeArgs),

    /// Extracts sequences based on Kraken2 results
    Extract(ExtractArgs),

    /// Generates taxonomic abundance matrices from multiple reports
    #[command(name = "abundance-matrix")]
    AbundanceMatrix(AbundanceMatrixArgs),

    /// Generates test data for performance testing
    #[command(name = "generate-test-data")]
    GenerateTestData(GenerateTestDataArgs),
}

/// Arguments for the 'analyze' command
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

/// Arguments for the 'extract' command
#[derive(Args)]
struct ExtractArgs {
    /// Input FASTA/FASTQ file
    sequence: String,

    /// Kraken2 log file
    log: String,

    /// Output file for extracted sequences
    #[arg(short, long)]
    output: String,

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

/// Arguments for the 'abundance-matrix' command
#[derive(Args)]
struct AbundanceMatrixArgs {
    /// Input Kraken2 report files (can be multiple)
    #[arg(required = true)]
    input: Vec<String>,

    /// Output file for the abundance matrix
    #[arg(short, long)]
    output: String,

    /// Output format (tsv or biom)
    #[arg(long, default_value = "tsv")]
    format: String,

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

/// Arguments for the 'generate-test-data' command
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

/// Main function to run the command-line interface
pub fn run_cli() {
    let cli = Cli::parse();

    let result = match cli.command {
        Commands::Analyze(args) => run_analyze(args),
        Commands::Extract(args) => run_extract(args),
        Commands::AbundanceMatrix(args) => run_abundance_matrix(args),
        Commands::GenerateTestData(args) => run_generate_test_data(args),
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}

/// Implements the 'analyze' command
fn run_analyze(args: AnalyzeArgs) -> Result<(), Box<dyn Error>> {
    let before_memory = memory_stats().map(|s| s.physical_mem).unwrap_or(0);
    let start_time = Instant::now();

    // Parse the Kraken2 report with proper error handling
    let (report, parse_time) = krk_parser::parse_kraken2_report(&args.report)
        .map_err(|e| format!("Error parsing Kraken2 report file '{}': {}", args.report, e))?;

    // If a specific taxon ID search was specified
    if let Some(taxon_id) = args.taxon_id {
        match find_taxon_info(&report.root, taxon_id as u64) {
            Some(taxon_info) => print_taxon_info(&taxon_info),
            None => println!("Taxon with ID {} not found", taxon_id),
        }
    }

    // Generate JSON output if requested
    if let Some(json_output) = args.json {
        krk_parser::write_json_report(&report, &json_output)
            .map_err(|e| format!("Error writing JSON report '{}': {}", json_output, e))?;
        println!("JSON report written to {}", json_output);
    }

    // Show performance metrics
    let after_memory = memory_stats().map(|s| s.physical_mem).unwrap_or(0);
    let memory_used = after_memory.saturating_sub(before_memory);
    let total_time = start_time.elapsed().as_secs_f64();

    println!("Total time: {:.6} seconds", total_time);
    println!("File parsing time: {:.6} seconds", parse_time);
    println!("Memory usage: {} bytes", memory_used);

    Ok(())
}

/// Implements the 'extract' command
fn run_extract(args: ExtractArgs) -> Result<(), Box<dyn Error>> {
    // Parsear los taxids de la línea de comandos
    let taxids: HashSet<String> = args
        .taxids
        .split(',')
        .map(|s| s.trim().to_string())
        .collect();

    // Usamos el conjunto original para las estadísticas
    let original_taxids = taxids.clone();

    // Almacena los readids después de la consulta
    let readids;

    // Almacena los mapeos de taxid a readids
    let mut taxid_readid_map: HashMap<String, HashSet<String>> = HashMap::new();

    // Crear un conjunto expandido para almacenar todos los taxids
    // (originales + padres/hijos)
    let mut expanded_taxids = taxids.clone();

    // If we need to include children or parents, we need the report file
    if (args.include_children || args.include_parents) && args.report.is_some() {
        let report_file = args.report.clone().unwrap();

        println!("Reading taxonomy from report file: {}", report_file);

        // Parse the Kraken2 report to get the taxonomic tree
        let (report, _) = match krk_parser::parse_kraken2_report(&report_file) {
            Ok(result) => result,
            Err(e) => {
                return Err(
                    format!("Error parsing Kraken2 report file '{}': {}", report_file, e).into(),
                )
            }
        };

        // Expand the set of taxids according to the options
        let mut expanded_count = 0;

        // Estimate initial size for the expanded set
        // For large taxonomies, could be 10x the original size
        let mut new_taxids = HashSet::with_capacity(if taxids.len() < 10 {
            // For small sets of taxids, we might have many children
            taxids.len() * 100
        } else {
            // For larger sets, the relative expansion is usually smaller
            taxids.len() * 10
        });

        for taxid in &taxids {
            let taxid_num = taxid.parse::<u32>().unwrap_or(0);

            // Add the original taxid
            new_taxids.insert(taxid.clone());

            if args.include_children {
                // If include_children is true, add all descendant taxids
                println!("Including children for taxid: {}", taxid);
                let mut child_taxids = HashSet::new();
                find_all_child_taxids(&report.root, taxid_num, &mut child_taxids);
                expanded_count += child_taxids.len();

                for child_taxid in child_taxids {
                    new_taxids.insert(child_taxid.to_string());
                }
            }

            if args.include_parents {
                // If include_parents is true, add all ancestor taxids
                println!("Including parents for taxid: {}", taxid);
                let mut parent_taxids = HashSet::new();
                find_all_parent_taxids(&report.root, taxid_num, &mut parent_taxids);
                expanded_count += parent_taxids.len();

                for parent_taxid in parent_taxids {
                    new_taxids.insert(parent_taxid.to_string());
                }
            }
        }

        // Replace the original taxids with the expanded set
        expanded_taxids = new_taxids;
        println!(
            "Expanded to {} taxids (added {} through hierarchy)",
            expanded_taxids.len(),
            expanded_count
        );
    } else if (args.include_children || args.include_parents) && args.report.is_none() {
        // Return an error instead of just a warning
        return Err("Error: A report file (--report) is required when using --include-children or --include-parents options.".into());
    }

    // Extract the sequences
    match logkrk_parser::parse_kraken_output_with_taxids(
        &args.log,
        &expanded_taxids,
        &mut taxid_readid_map,
    ) {
        Ok(ids) => {
            readids = ids;
            // Cambiar la siguiente línea si total_sequences no se usa después
            let total_sequences = count_sequences_in_file(&args.sequence)?;

            match sequence_processor::process_sequence_files(
                std::slice::from_ref(&args.sequence),
                &readids,
                &args.output,
                args.exclude,
            ) {
                Ok(_) => {
                    println!("Sequences extracted successfully to {}", args.output);
                    println!(
                        "{} sequences matching {} taxids",
                        readids.len(),
                        expanded_taxids.len()
                    );
                    println!(
                        "Total {} sequences: {}",
                        if args.exclude {
                            "excluded"
                        } else {
                            "extracted"
                        },
                        readids.len()
                    );
                }
                Err(e) => return Err(Box::new(std::io::Error::other(e.to_string()))),
            }

            // Generate statistics file if requested
            if let Some(ref stats_file) = args.stats_output {
                match generate_statistics_file(
                    stats_file,
                    &taxid_readid_map,
                    &original_taxids,
                    total_sequences,
                    &args,
                ) {
                    Ok(_) => println!("Statistics written to {}", stats_file),
                    Err(e) => eprintln!("Error writing statistics: {}", e),
                }
            }
        }
        Err(e) => return Err(format!("Error parsing Kraken2 log file: {}", e).into()),
    }

    Ok(())
}

/// Counts the total sequences in a FASTA/FASTQ file
fn count_sequences_in_file(filename: &str) -> Result<usize, Box<dyn Error>> {
    let file = File::open(filename)?;

    // Use a large buffer (512 KB) for efficient reading
    let mut reader = BufReader::with_capacity(BUFFER_SIZE, file);
    let mut buffer = Vec::with_capacity(1024);

    // Detect the file format by reading the first character
    reader.read_until(b'\n', &mut buffer)?;
    match buffer.first() {
        Some(b'>') => count_fasta_sequences(reader),
        Some(b'@') => count_fastq_sequences(reader),
        _ => Err("Unknown format or empty file".into()), // Empty file or error
    }
}

/// Counts sequences in a FASTA file
fn count_fasta_sequences(mut reader: BufReader<File>) -> Result<usize, Box<dyn Error>> {
    // Counters for statistics
    let mut sequence_count = 1; // We already read the first line
    let mut buffer = Vec::with_capacity(1024);

    // Identify bytes for sequence identifiers
    let header_char = b'>';

    // FASTA file
    // Count each line starting with '>'
    loop {
        buffer.clear();
        let bytes_read = reader.read_until(b'\n', &mut buffer)?;

        if bytes_read == 0 {
            break;
        }

        if !buffer.is_empty() && buffer[0] == header_char {
            sequence_count += 1;
        }
    }

    Ok(sequence_count)
}

/// Counts sequences in a FASTQ file
fn count_fastq_sequences(mut reader: BufReader<File>) -> Result<usize, Box<dyn Error>> {
    // Counters for statistics
    let mut sequence_count = 1; // We already read the first line
    let mut line_count = 1;
    let mut buffer = Vec::with_capacity(1024);

    // FASTQ file
    // In FASTQ each record has 4 lines, only count the first one
    loop {
        buffer.clear();
        let bytes_read = reader.read_until(b'\n', &mut buffer)?;

        if bytes_read == 0 {
            break;
        }

        line_count += 1;

        if line_count % 4 == 1 {
            if !buffer.is_empty() && buffer[0] == b'@' {
                sequence_count += 1;
            } else {
                eprintln!("Warning: Unrecognized sequence file format, assuming FASTA");
                return count_fasta_sequences(reader);
            }
        }
    }

    Ok(sequence_count)
}

/// Generates a detailed statistics file
fn generate_statistics_file(
    stats_file: &str,
    taxid_readid_map: &HashMap<String, HashSet<String>>,
    original_taxids: &HashSet<String>,
    total_sequences: usize,
    args: &ExtractArgs,
) -> Result<(), Box<dyn Error>> {
    let file = File::create(stats_file)?;

    // Create a large buffer for writing
    let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);

    // Calculate statistics once to avoid repeated calculations
    let mut total_extracted = 0;
    let mut total_original_taxids = 0;
    let mut total_expanded_taxids = 0;

    // Counters for summary statistics
    let mut orig_sequences = 0;
    let mut expanded_sequences = 0;

    // Pre-calculate taxid counts and statistics all at once
    // This avoids traversing the HashMap multiple times
    let mut stats: Vec<(String, usize, bool)> = Vec::with_capacity(taxid_readid_map.len());

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

        stats.push((taxid.clone(), count, is_original));
    }

    // Sort by number of sequences (descending)
    stats.sort_by_key(|entry| std::cmp::Reverse(entry.1));

    // Calculate percentages for the summary
    let percent_extracted = if total_sequences > 0 {
        (total_extracted as f64 / total_sequences as f64) * 100.0
    } else {
        0.0
    };

    // Write metadata as comments (lines starting with #)
    // These will be recognized as comments by pandas, R, and other tools
    writeln!(writer, "# KrakenClip Extraction Statistics")?;
    writeln!(
        writer,
        "# Date: {}",
        chrono::Local::now().format("%Y-%m-%d %H:%M:%S")
    )?;
    writeln!(writer, "# Input file: {}", args.sequence)?;
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

    // Write the CSV header - common CSV format with clear field names
    // Use descriptive column names compatible with data analysis software
    writeln!(
        writer,
        "taxid,sequences,percent_of_extracted,percent_of_total,is_original"
    )?;

    // Write detailed statistics as CSV rows
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

    // Ensure all data is written
    writer.flush()?;

    Ok(())
}

/// Recursively finds all child taxids
fn find_all_child_taxids(
    node: &krk_parser::TaxonEntry,
    target_taxid: u32,
    taxids: &mut HashSet<String>,
) {
    // Pre-allocate space to avoid reallocations
    // We estimate that each node might have approximately 10 children with a tree depth of 10
    if taxids.capacity() < 100 {
        taxids.reserve(100);
    }

    // Check if this node is the target
    if node.taxid == target_taxid {
        // If found, add all children recursively
        add_all_children_to_taxids(node, taxids);
        return;
    }

    // Otherwise, search in the children
    for child in &node.children {
        find_all_child_taxids(child, target_taxid, taxids);
    }
}

/// Adds all children of a node to the taxid set
fn add_all_children_to_taxids(node: &krk_parser::TaxonEntry, taxids: &mut HashSet<String>) {
    // Reserve space to avoid frequent reallocations
    if taxids.capacity() < taxids.len() + node.children.len() {
        taxids.reserve(node.children.len() * 2);
    }

    for child in &node.children {
        taxids.insert(child.taxid.to_string());
        // Process the children recursively
        add_all_children_to_taxids(child, taxids);
    }
}

/// Recursively finds all parent taxids
fn find_all_parent_taxids(
    root: &krk_parser::TaxonEntry,
    target_taxid: u32,
    taxids: &mut HashSet<String>,
) {
    // Use a vector with pre-allocated capacity for the path
    // Most taxonomies don't exceed 30-40 levels of depth
    let mut path = Vec::with_capacity(50);

    // Find parents using DFS with path tracking
    find_parents_dfs(root, target_taxid, &mut path, taxids);
}

fn find_parents_dfs(
    node: &krk_parser::TaxonEntry,
    target_taxid: u32,
    path: &mut Vec<u32>,
    taxids: &mut HashSet<String>,
) -> bool {
    // Add current node to the path
    path.push(node.taxid);

    // If this is the target node, add all parents from the path
    if node.taxid == target_taxid {
        // Reserve space for the number of parents to add
        if taxids.capacity() < taxids.len() + path.len() {
            taxids.reserve(path.len());
        }

        // Add all nodes from the path except the last (current node)
        for &parent_id in path.iter().take(path.len() - 1) {
            taxids.insert(parent_id.to_string());
        }

        path.pop(); // Restore path for backtracking
        return true;
    }

    // Search in the children
    for child in &node.children {
        if find_parents_dfs(child, target_taxid, path, taxids) {
            path.pop(); // Restore path for backtracking
            return true;
        }
    }

    // Target node not found in this subtree
    path.pop(); // Restore path for backtracking
    false
}

/// Implements the 'abundance-matrix' command
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
    if args.format != "tsv" && args.format != "biom" {
        return Err(format!(
            "Unsupported output format '{}'. Use 'tsv' or 'biom'.",
            args.format
        )
        .into());
    }

    let mut matrix = AbundanceMatrix::new(&args.level);
    matrix.set_force_include_unclassified(args.include_unclassified);
    let proportional = args.normalize || args.proportions || !args.absolute_counts;
    let mut used_sample_names = HashSet::new();

    for (index, file) in args.input.iter().enumerate() {
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

        println!("Processing sample: {sample_name}");
        let (report, _) = krk_parser::parse_kraken2_report(file)
            .map_err(|error| format!("Error parsing Kraken2 report file '{file}': {error}"))?;
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
            BiomTable::from_abundance_matrix(&matrix)
                .write_json(&args.output)
                .map_err(|error| format!("Error generating BIOM output: {error}"))?;
            println!(
                "BIOM format output successfully generated in: {}",
                args.output
            );
        }
        _ => unreachable!("format was validated above"),
    }

    Ok(())
}

/// Implements the 'generate-test-data' command
fn run_generate_test_data(args: GenerateTestDataArgs) -> Result<(), Box<dyn Error>> {
    // Add aggregated information as needed
    match generate_test_data::generate_data(&args.output, args.lines, &args.r#type) {
        Ok(_) => println!("Test data generated successfully"),
        Err(e) => return Err(format!("Error generating test data: {}", e).into()),
    }

    Ok(())
}
