use rand::RngExt;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;

/// Structure to define generation parameters
pub struct GeneratorParams {
    pub output_file: String,
    pub num_lines: usize,
    pub max_depth: usize,
    pub max_children: usize,
    pub max_fragments: u64,
}

impl Default for GeneratorParams {
    fn default() -> Self {
        Self {
            output_file: "test_data_large.txt".to_string(),
            num_lines: 100_000,
            max_depth: 30,
            max_children: 50,
            max_fragments: 10_000_000,
        }
    }
}

pub fn generate_data(
    output_file: &str,
    num_lines: usize,
    data_type: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    match data_type {
        "bracken" => return generate_bracken_data(output_file, num_lines),
        "log" => return generate_log_data(output_file, num_lines),
        _ => {}
    }

    let mut params = GeneratorParams {
        output_file: output_file.to_string(),
        num_lines,
        ..Default::default()
    };

    // Configure parameters according to the required data type
    match data_type {
        "wide" => {
            // Wide but shallow taxonomy
            params.max_depth = 5;
            params.max_children = 100;
        }
        "deep" => {
            // Narrow but very deep taxonomy
            params.max_depth = 50;
            params.max_children = 3;
        }
        "fragments" => {
            // Many fragments per node
            params.max_fragments = 100_000_000;
        }
        "dense" => {
            // Dense taxonomy with few nodes
            params.max_depth = 10;
            params.max_children = 10;
            params.max_fragments = 5_000_000;
        }
        "lines" => {
            // Many lines but few fragments
            params.num_lines = num_lines * 3;
            params.max_fragments = 1_000_000;
        }
        "random" => {
            // Random distribution
            let mut rng = rand::rng();
            params.max_depth = rng.random_range(5..40);
            params.max_children = rng.random_range(3..80);
            params.max_fragments = rng.random_range(1_000_000..50_000_000);
        }
        "complex" => {
            // Complex taxonomy
            params.max_depth = 35;
            params.max_children = 70;
            params.max_fragments = 30_000_000;
        }
        "extreme" => {
            // Extremely large
            params.max_depth = 50;
            params.max_children = 100;
            params.max_fragments = 100_000_000;
        }
        "unbalanced" => {
            // Unbalanced taxonomy (varying depth and children)
            params.max_depth = 60;
            params.max_children = 40;
        }
        "mixed" => {
            // Mixed characteristics
            params.max_depth = 25;
            params.max_children = 60;
            params.max_fragments = 20_000_000;
        }
        _ => {
            // Default - balanced
            params.max_depth = 15;
            params.max_children = 20;
            params.max_fragments = 5_000_000;
        }
    }

    // Generate the test file with the configured parameters
    generate_test_data(&params)
}

/// Generate test data file for performance testing
///
/// Creates a customizable test data file with random taxonomic structure.
pub fn generate_test_data(params: &GeneratorParams) -> Result<(), Box<dyn std::error::Error>> {
    if params.num_lines < 2 {
        return Err("at least two lines are required for unclassified and root entries".into());
    }

    let start = Instant::now();

    let file = File::create(&params.output_file)?;
    let mut writer = BufWriter::new(file);

    let max_fragments = params.max_fragments;
    let total_fragments = max_fragments / 10; // Use 10% for unclassified

    // Write unclassified line
    writeln!(
        writer,
        "0.00\t{}\t{}\tU\t0\tunclassified",
        total_fragments, total_fragments
    )?;

    // Write root node
    writeln!(
        writer,
        "100.00\t{}\t0\tR\t1\troot",
        max_fragments - total_fragments
    )?;

    // Generate and distribute fragments among children
    let mut rng = rand::rng();
    let child_fragments = max_fragments - total_fragments;
    let mut remaining_lines = params.num_lines - 2;

    let mut total_lines = 2; // Already wrote 2 lines
    let fragments_generated = distribute_fragments(
        &mut writer,
        child_fragments,
        1,
        params.max_depth,
        params.max_children,
        &mut remaining_lines,
        &mut rng,
    )?;

    total_lines += fragments_generated.0;

    let elapsed = start.elapsed();
    println!("Generated {} lines in {:.2?}", total_lines, elapsed);

    Ok(())
}

// Simple weighted distribution struct
struct WeightedDistribution {
    weights: Vec<f64>,
    total: f64,
}

impl WeightedDistribution {
    fn new(count: usize) -> Self {
        let mut rng = rand::rng();
        let mut weights = Vec::with_capacity(count);
        let mut total = 0.0;

        for _ in 0..count {
            let weight = rng.random_range(1.0..10.0);
            weights.push(weight);
            total += weight;
        }

        Self { weights, total }
    }
}

// Distribute fragments recursively among taxonomy nodes
//
// Returns a tuple of (lines_generated, fragments_assigned)
fn distribute_fragments<R: RngExt>(
    writer: &mut BufWriter<File>,
    fragments: u64,
    depth: usize,
    max_depth: usize,
    max_children: usize,
    remaining_lines: &mut usize,
    rng: &mut R,
) -> Result<(usize, u64), Box<dyn std::error::Error>> {
    if fragments == 0 || depth > max_depth || *remaining_lines == 0 {
        return Ok((0, 0));
    }

    // Generate a random number of children based on depth
    let base_children = if depth < 3 {
        let half = (max_children / 2).max(1);
        half + rng.random_range(0..half)
    } else {
        rng.random_range(0..max_children.max(1))
    };

    // Reduce children as depth increases
    let depth_factor = 1.0 - (depth as f64 / max_depth as f64) * 0.8;
    let num_children =
        ((base_children as f64 * depth_factor).round() as usize).min(*remaining_lines);

    if num_children == 0 {
        return Ok((0, 0));
    }

    // Create weighted distribution for fragment allocation
    let distribution = WeightedDistribution::new(num_children);

    // Assign fragments proportionally in O(children), not O(fragments).
    let mut child_fragments = vec![0_u64; num_children];
    let direct_fragments = fragments / 5; // 20% direct, 80% in children
    let remaining = fragments - direct_fragments;

    let mut assigned = 0_u64;
    for (index, weight) in distribution.weights.iter().enumerate() {
        let allocation = ((remaining as f64) * (*weight / distribution.total)) as u64;
        child_fragments[index] = allocation;
        assigned += allocation;
    }
    for index in 0..(remaining - assigned) as usize {
        child_fragments[index % num_children] += 1;
    }

    // Generate random taxon IDs for children
    let mut child_ids = Vec::with_capacity(num_children);
    for _ in 0..num_children {
        let id = rng.random_range(10..1_000_000);
        child_ids.push(id);
    }

    // Generate taxonomy ranks based on depth
    let rank = match depth {
        1 => "D", // Domain
        2 => "P", // Phylum
        3 => "C", // Class
        4 => "O", // Order
        5 => "F", // Family
        6 => "G", // Genus
        _ => "S", // Species
    };

    // Write lines for all children
    let mut total_lines = 0;
    let mut total_fragments = 0;

    for i in 0..num_children {
        let taxid = child_ids[i];
        let fragments_this_child = child_fragments[i];
        if fragments_this_child == 0 {
            continue;
        }
        if *remaining_lines == 0 {
            break;
        }

        let percentage = (fragments_this_child as f64 / fragments as f64) * 100.0;
        let name = format!("taxon_{}", taxid);
        let indent = "  ".repeat(depth);

        // Write the child's line
        writeln!(
            writer,
            "{:.2}\t{}\t{}\t{}\t{}\t{}{}",
            percentage,
            fragments_this_child,
            rng.random_range(0..fragments_this_child + 1), // Direct fragments
            rank,
            taxid,
            indent,
            name
        )?;
        total_lines += 1;
        *remaining_lines -= 1;

        // Recursively distribute fragments to this child's descendants
        let (child_lines, _child_frags) = distribute_fragments(
            writer,
            fragments_this_child,
            depth + 1,
            max_depth,
            max_children,
            remaining_lines,
            rng,
        )?;

        total_lines += child_lines;
        total_fragments += fragments_this_child;
    }

    Ok((total_lines, total_fragments + direct_fragments))
}

fn generate_bracken_data(
    output_file: &str,
    num_lines: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::create(output_file)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads"
    )?;
    let mut rng = rand::rng();
    let mut total = 0_u64;
    let mut rows = Vec::with_capacity(num_lines);
    for i in 0..num_lines {
        let reads = rng.random_range(1..1_000_u64);
        total += reads;
        rows.push((i as u32 + 10, reads));
    }
    for (taxid, reads) in rows {
        let fraction = reads as f64 / total as f64;
        writeln!(
            writer,
            "taxon_{taxid}\t{taxid}\tS\t{reads}\t0\t{reads}\t{fraction:.10}"
        )?;
    }
    writer.flush()?;
    println!("Generated {num_lines} Bracken rows");
    Ok(())
}

fn generate_log_data(
    output_file: &str,
    num_lines: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::create(output_file)?;
    let mut writer = BufWriter::new(file);
    let mut rng = rand::rng();
    for i in 0..num_lines {
        let taxid = if i % 17 == 0 {
            0
        } else {
            rng.random_range(2..500)
        };
        let status = if taxid == 0 { "U" } else { "C" };
        writeln!(writer, "{status}\tread{i}\t{taxid}\t100\t{taxid}:100")?;
    }
    writer.flush()?;
    println!("Generated {num_lines} Kraken log lines");
    Ok(())
}

/// Run the generator with the specified parameters
#[allow(dead_code)]
pub fn run_generator(
    output_file: &str,
    num_lines: usize,
    max_depth: usize,
    max_children: usize,
    max_fragments: u64,
) -> Result<(), Box<dyn std::error::Error>> {
    let params = GeneratorParams {
        output_file: output_file.to_string(),
        num_lines,
        max_depth,
        max_children,
        max_fragments,
    };

    generate_test_data(&params)
}
