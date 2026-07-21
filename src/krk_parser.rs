use memchr::memchr_iter;
use std::fs::File;
use std::io::Read;
use std::io::Write;
use std::path::Path;
use std::time::Instant;

// Optimized constants for performance-critical operations
// Buffer size is carefully chosen for optimal memory usage vs. throughput
const BUFFER_SIZE: usize = 512 * 1024; // 512KB balances memory usage and read performance
                                       // Raw byte literals are used for faster comparison operations (avoids UTF-8 decoding)
const TAB_CHAR: u8 = b'\t';
const SPACE_CHAR: u8 = b' ';
const NEWLINE_CHAR: u8 = b'\n';

/// Structure representing a node in the taxonomic tree
/// Each node contains comprehensive information about a specific taxon
#[derive(Debug, Clone)]
pub struct TaxonEntry {
    pub percentage: f32,   // Percentage of reads in the sample assigned to this clade
    pub clade_reads: u64,  // Total reads assigned to this clade and its descendants
    pub direct_reads: u64, // Reads assigned directly to this taxon (not descendants)
    pub rank: String,      // Taxonomic rank (e.g., "P" for phylum, "G" for genus)
    pub taxid: u32,        // NCBI Taxonomy identifier as u32 for memory efficiency
    pub name: String,      // Scientific name of the taxon
    #[allow(dead_code)]
    pub depth: usize, // Depth in the taxonomy tree (number of ancestors)
    pub children: Vec<TaxonEntry>, // Child nodes (descendant taxa)

    // Additional fields for API compatibility with other formats
    pub level: usize, // Hierarchy level (same as depth, maintained for compatibility)
    pub taxon_id: u64, // Taxon ID as u64 (for compatibility with external systems)
    pub clade_fragments: u64, // Alias for clade_reads (terminology varies by tool)
    pub direct_fragments: u64, // Alias for direct_reads (terminology varies by tool)
    pub rank_code: String, // Alias for rank (maintained for compatibility)
}

/// Main structure for representing a complete Kraken report
/// Contains the full taxonomic hierarchy.
#[derive(Debug)]
pub struct KrakenReport {
    pub root: TaxonEntry,                 // Root node of the taxonomic tree
    pub unclassified: Option<TaxonEntry>, // Special node for unclassified sequences (if present)
}

impl TaxonEntry {
    // Create a new empty taxonomy node
    fn new(
        percentage: f32,
        clade_reads: u64,
        direct_reads: u64,
        rank: String,
        taxid: u32,
        name: String,
        depth: usize,
    ) -> Self {
        Self {
            percentage,
            clade_reads,
            direct_reads,
            rank: rank.clone(), // Clone is needed as we use the value twice
            taxid,
            name,
            depth,
            children: Vec::new(),
            // Initialize additional fields for compatibility
            level: depth,
            taxon_id: taxid as u64,
            clade_fragments: clade_reads,
            direct_fragments: direct_reads,
            rank_code: rank,
        }
    }

    /// Add a child node to this taxon
    /// This establishes the parent-child relationship in the taxonomy tree
    #[allow(dead_code)]
    fn add_child(&mut self, child: TaxonEntry) {
        self.children.push(child);
    }
}

impl Default for KrakenReport {
    /// Create a default empty Kraken report with just a root node
    /// This implements the Default trait for KrakenReport
    fn default() -> Self {
        Self {
            root: TaxonEntry::new(0.0, 0, 0, "R".to_string(), 1, "root".to_string(), 0),
            unclassified: None,
        }
    }
}

/// String cache to avoid memory duplication for common rank codes
/// This significantly reduces memory usage for large reports
pub struct StringCache {
    rank_codes: Vec<String>,
}

impl StringCache {
    /// Create a new string cache pre-populated with common rank codes
    ///
    /// # Performance optimization
    /// Pre-populating with common values avoids allocations for frequent codes
    fn new() -> Self {
        // Pre-populate with common rank codes used in taxonomy
        // D=Domain, P=Phylum, C=Class, O=Order, F=Family, G=Genus, S=Species, U=Unclassified
        let mut rank_codes = Vec::with_capacity(10);
        for code in ["D", "P", "C", "O", "F", "G", "S", "U"].iter() {
            rank_codes.push((*code).to_string());
        }
        Self { rank_codes }
    }

    /// Get a rank code string from cache or add if not present
    ///
    /// # Arguments
    /// * `code` - The rank code to retrieve or add
    ///
    /// # Returns
    /// A clone of the cached string (reduces total allocations)
    ///
    /// # Performance characteristics
    /// Uses linear search which is efficient for small collections like this
    /// For the small set of rank codes, this outperforms a HashMap due to lower overhead
    fn get_rank_code(&mut self, code: &str) -> String {
        // First check if code exists in cache to avoid allocation
        for existing in &self.rank_codes {
            if existing == code {
                return existing.clone();
            }
        }
        // If not in cache, add it
        let code_string = code.to_string();
        self.rank_codes.push(code_string.clone());
        code_string
    }
}

// Highly optimized line parser using pointers and avoiding unnecessary allocations
#[inline(always)]
pub fn parse_line(
    line: &[u8],
    string_cache: &mut StringCache,
    line_number: Option<usize>,
) -> Option<TaxonEntry> {
    // Use fast tab search with memchr for vectorized byte searching
    // This is significantly faster than manual iteration
    let tab_positions: Vec<usize> = memchr_iter(TAB_CHAR, line).take(5).collect();
    if tab_positions.len() < 5 {
        if let Some(line_num) = line_number {
            eprintln!("Warning: Line {} does not have enough tab-separated fields (expected at least 5, found {})", 
                     line_num, tab_positions.len());
        }
        return None; // Not enough fields
    }

    // Optimization: calculate field positions directly from tab positions
    // Avoids string splits which would create new allocations
    let field_starts = [
        0,
        tab_positions[0] + 1,
        tab_positions[1] + 1,
        tab_positions[2] + 1,
        tab_positions[3] + 1,
        tab_positions[4] + 1,
    ];
    let field_ends = [
        tab_positions[0],
        tab_positions[1],
        tab_positions[2],
        tab_positions[3],
        tab_positions[4],
        line.len(),
    ];

    // Fast calculation of indentation level (two spaces per level)
    // This determines the taxon's position in the hierarchy
    let name_start = field_starts[5];
    let mut level = 0;
    let mut i = name_start;

    // Count spaces quickly in pairs (Kraken uses two spaces per level)
    // This optimization is faster than counting individual spaces
    while i + 1 < line.len() && line[i] == SPACE_CHAR && line[i + 1] == SPACE_CHAR {
        level += 1;
        i += 2;
    }

    // Extract name efficiently by slicing the original buffer
    // Avoids allocating a new string until necessary
    let name_bytes = &line[name_start + level * 2..];
    let name = match std::str::from_utf8(name_bytes) {
        Ok(s) => s.trim(),
        Err(e) => {
            if let Some(line_num) = line_number {
                eprintln!(
                    "Warning: UTF-8 encoding error in taxon name at line {}: {}",
                    line_num, e
                );
            }
            ""
        }
    };

    let percentage_bytes = &line[field_starts[0]..field_ends[0]];
    let percentage = match std::str::from_utf8(percentage_bytes) {
        Ok(s) => s.parse::<f32>().unwrap_or_else(|_| {
            if let Some(line_num) = line_number {
                eprintln!(
                    "Warning: Failed to parse percentage value '{}' at line {}",
                    s, line_num
                );
            }
            0.0
        }),
        Err(_) => {
            if let Some(line_num) = line_number {
                eprintln!(
                    "Warning: Invalid UTF-8 in percentage field at line {}",
                    line_num
                );
            }
            0.0
        }
    };

    // Parse integer counts directly from byte slices
    let clade_bytes = &line[field_starts[1]..field_ends[1]];
    let clade_fragments = match std::str::from_utf8(clade_bytes) {
        Ok(s) => s.parse::<u64>().unwrap_or_else(|_| {
            if let Some(line_num) = line_number {
                eprintln!(
                    "Warning: Failed to parse clade reads value '{}' at line {}",
                    s, line_num
                );
            }
            0
        }),
        Err(_) => {
            if let Some(line_num) = line_number {
                eprintln!(
                    "Warning: Invalid UTF-8 in clade reads field at line {}",
                    line_num
                );
            }
            0
        }
    };

    let direct_bytes = &line[field_starts[2]..field_ends[2]];
    let direct_fragments = match std::str::from_utf8(direct_bytes) {
        Ok(s) => s.parse::<u64>().unwrap_or_else(|_| {
            if let Some(line_num) = line_number {
                eprintln!(
                    "Warning: Failed to parse direct reads value '{}' at line {}",
                    s, line_num
                );
            }
            0
        }),
        Err(_) => {
            if let Some(line_num) = line_number {
                eprintln!(
                    "Warning: Invalid UTF-8 in direct reads field at line {}",
                    line_num
                );
            }
            0
        }
    };

    // Optimization for rank codes using string cache to avoid duplication
    // This significantly reduces memory usage for large reports with many identical ranks
    let rank_bytes = &line[field_starts[3]..field_ends[3]];
    let rank_str = match std::str::from_utf8(rank_bytes) {
        Ok(s) => s.trim(),
        Err(_) => {
            if let Some(line_num) = line_number {
                eprintln!("Warning: Invalid UTF-8 in rank field at line {}", line_num);
            }
            ""
        }
    };
    let rank_code = string_cache.get_rank_code(rank_str);

    // Parse taxon ID
    let taxon_bytes = &line[field_starts[4]..field_ends[4]];
    let taxon_id = match std::str::from_utf8(taxon_bytes) {
        Ok(s) => s.parse::<u64>().unwrap_or_else(|_| {
            if let Some(line_num) = line_number {
                eprintln!(
                    "Warning: Failed to parse taxon ID value '{}' at line {}",
                    s, line_num
                );
            }
            0
        }),
        Err(_) => {
            if let Some(line_num) = line_number {
                eprintln!(
                    "Warning: Invalid UTF-8 in taxon ID field at line {}",
                    line_num
                );
            }
            0
        }
    };

    // Construct and return the taxon entry
    Some(TaxonEntry {
        level,
        percentage,
        clade_fragments,
        direct_fragments,
        rank_code: rank_code.clone(),
        taxid: taxon_id as u32,
        name: name.to_string(),
        depth: level,
        children: Vec::new(),
        taxon_id,
        clade_reads: clade_fragments,
        direct_reads: direct_fragments,
        rank: rank_code,
    })
}

// Building hierarchy with optimal memory reservation
pub fn build_hierarchy_optimized(entries: Vec<TaxonEntry>) -> Vec<TaxonEntry> {
    if entries.is_empty() {
        return Vec::new();
    }

    // Pre-allocate space for the final hierarchy
    // Most reports have few top-level nodes compared to total entries
    let mut hierarchy = Vec::with_capacity(entries.len() / 8 + 1);

    // Use a stack with pre-allocated capacity
    // Taxonomic hierarchies rarely exceed 20 levels in depth
    let mut stack: Vec<TaxonEntry> = Vec::with_capacity(20);

    for entry in entries {
        // Process entries based on their level in the hierarchy
        // This implements a non-recursive tree construction algorithm

        // Optimization to avoid chain reactions of push/pop
        // Pop entries from stack until we find the parent for current entry
        while !stack.is_empty() && stack.last().unwrap().level >= entry.level {
            let popped = stack.pop().unwrap();
            if let Some(parent) = stack.last_mut() {
                // Add the popped entry as a child of its parent

                // Reserve capacity if needed to avoid reallocations during push
                // Most taxonomic nodes have relatively few direct children
                if parent.children.is_empty() {
                    parent.children.reserve(4); // Most nodes have fewer than 4 children
                }
                parent.children.push(popped);
            } else {
                // If stack is empty, this is a root-level node
                hierarchy.push(popped);
            }
        }
        stack.push(entry);
    }

    // Process remaining elements in the stack (cleanup phase)
    // This handles nodes that remain after all entries are processed
    while let Some(entry) = stack.pop() {
        if let Some(parent) = stack.last_mut() {
            parent.children.push(entry);
        } else {
            hierarchy.push(entry);
        }
    }

    hierarchy
}

/// Custom buffer implementation for optimized file reading
/// Provides block-based reading with minimal allocations
struct OptimizedBuffer {
    buffer: Box<[u8]>, // Fixed-size buffer using Box for heap allocation
    pos: usize,        // Current position in the buffer
    cap: usize,        // Current capacity (filled bytes) in the buffer
    file: File,        // Source file
}

impl OptimizedBuffer {
    fn new(file: File) -> Self {
        Self {
            // Use boxed slice for more efficient memory layout
            buffer: vec![0; BUFFER_SIZE].into_boxed_slice(),
            pos: 0,
            cap: 0,
            file,
        }
    }

    fn fill_buffer(&mut self) -> std::io::Result<usize> {
        self.pos = 0;
        self.cap = self.file.read(&mut self.buffer)?;
        Ok(self.cap)
    }

    fn read_line(&mut self, line_buffer: &mut Vec<u8>) -> std::io::Result<bool> {
        line_buffer.clear();

        // Check if buffer needs to be filled
        if self.pos >= self.cap {
            match self.fill_buffer() {
                Ok(0) => return Ok(false), // EOF
                Ok(_) => {}                // Successfully filled buffer
                Err(e) => {
                    return Err(std::io::Error::new(
                        e.kind(),
                        format!("Failed to read from file: {}", e),
                    ))
                }
            }
        }

        loop {
            // Fast scan for newline character in current buffer
            let mut i = self.pos;
            while i < self.cap {
                if self.buffer[i] == NEWLINE_CHAR {
                    // Found line ending - extract the line
                    line_buffer.extend_from_slice(&self.buffer[self.pos..i]);
                    self.pos = i + 1;
                    return Ok(true);
                }
                i += 1;
            }

            // No line break found, copy rest of buffer
            line_buffer.extend_from_slice(&self.buffer[self.pos..self.cap]);

            // Fill the buffer again
            match self.fill_buffer() {
                Ok(0) => {
                    // EOF - return true if we read any data
                    return Ok(!line_buffer.is_empty());
                }
                Ok(_) => {} // Successfully filled buffer
                Err(e) => {
                    return Err(std::io::Error::new(
                        e.kind(),
                        format!("Failed to read next buffer block: {}", e),
                    ))
                }
            }
        }
    }
}

// Main analysis function with block processing
pub fn parse_kraken2_report(file_path: &str) -> Result<(KrakenReport, f64), std::io::Error> {
    // Replace expect with proper error handling
    let file = File::open(file_path)?;

    // Estimate size for pre-allocation based on file size
    // Average line length in Kraken reports is ~50 bytes
    let file_size = file.metadata().map(|m| m.len() as usize).unwrap_or(0);

    let mut buffer = OptimizedBuffer::new(file);

    // Start timing the parsing process (excluding file opening)
    let start_time = Instant::now();
    let mut string_cache = StringCache::new();

    // Pre-allocate line buffer with reasonable size
    let mut line_buffer = Vec::with_capacity(1024);

    // Pre-allocate entries vector based on estimated line count
    let mut entries = Vec::with_capacity(file_size / 50); // Approximate estimation

    // Parse file line by line with line number tracking
    let mut line_number = 1;
    while buffer.read_line(&mut line_buffer)? {
        if let Some(entry) = parse_line(&line_buffer, &mut string_cache, Some(line_number)) {
            entries.push(entry);
        }
        line_number += 1;
    }

    // Build hierarchical representation from flat entries
    let hierarchy = build_hierarchy_optimized(entries);

    // Determine if there's an "unclassified" node
    // This is typically the first node in Kraken reports
    let unclassified = if !hierarchy.is_empty()
        && hierarchy[0].level == 0
        && hierarchy[0].name == "unclassified"
    {
        Some(hierarchy[0].clone())
    } else {
        None
    };

    // Get root node (typically the second node if unclassified is present)
    let root = if hierarchy.len() > 1 || (hierarchy.len() == 1 && unclassified.is_some()) {
        let root_index = if unclassified.is_some() { 1 } else { 0 };
        if root_index < hierarchy.len() {
            hierarchy[root_index].clone()
        } else {
            // Fallback to default root if hierarchy structure is unexpected
            TaxonEntry::new(0.0, 0, 0, "R".to_string(), 1, "root".to_string(), 0)
        }
    } else {
        // Create default root if no suitable node found
        TaxonEntry::new(0.0, 0, 0, "R".to_string(), 1, "root".to_string(), 0)
    };

    // Calculate total parsing time
    let duration = start_time.elapsed().as_secs_f64();

    // Return the constructed report and parsing duration
    Ok((KrakenReport { unclassified, root }, duration))
}

// Function to write the report in JSON format
pub fn write_json_report(report: &KrakenReport, output_path: &str) -> std::io::Result<()> {
    use serde_json;

    // Recursive function to convert a node to JSON format
    fn node_to_json(node: &TaxonEntry) -> serde_json::Value {
        // Create a JSON object with all relevant fields
        let mut json = serde_json::json!({
            "name": node.name,
            "taxid": node.taxid,
            "rank": node.rank,
            "percentage": node.percentage,
            "clade_reads": node.clade_reads,
            "direct_reads": node.direct_reads,
            "level": node.level,
            "children": []
        });

        // Add children recursively
        let children = node.children.iter().map(node_to_json).collect::<Vec<_>>();
        json["children"] = serde_json::Value::Array(children);

        json
    }

    // Convert the entire report to JSON
    let mut json = serde_json::json!({
        "root": node_to_json(&report.root)
    });

    // Add unclassified node if it exists
    if let Some(ref unclassified) = report.unclassified {
        json["unclassified"] = node_to_json(unclassified);
    }

    // Write to file using buffered I/O for performance
    let file = std::fs::File::create(Path::new(output_path))?;
    let mut writer = std::io::BufWriter::new(file);

    // Use pretty printing for human-readable output
    serde_json::to_writer_pretty(&mut writer, &json)?;
    writer.flush()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_line() {
        let test_line = b"50.00\t1000\t500\tP\t123\t  Bacteria";
        let mut cache = StringCache::new();
        let result = parse_line(test_line, &mut cache, None);
        assert!(result.is_some());
        let entry = result.unwrap();
        assert_eq!(entry.level, 1);
        assert_eq!(entry.percentage, 50.0);
        assert_eq!(entry.clade_fragments, 1000);
        assert_eq!(entry.direct_fragments, 500);
        assert_eq!(entry.rank_code, "P");
        assert_eq!(entry.taxid, 123);
        assert_eq!(entry.name, "Bacteria");
    }

    #[test]
    fn test_build_hierarchy() {
        let entries = vec![
            TaxonEntry {
                level: 0,
                percentage: 100.0,
                clade_fragments: 1000,
                direct_fragments: 0,
                rank_code: "D".to_string(),
                taxid: 1,
                name: "Root".to_string(),
                depth: 0,
                children: Vec::new(),
                // Add missing fields for the test (these would normally be set by new())
                taxon_id: 1,
                clade_reads: 1000,
                direct_reads: 0,
                rank: "D".to_string(),
            },
            TaxonEntry {
                level: 1,
                percentage: 80.0,
                clade_fragments: 800,
                direct_fragments: 200,
                rank_code: "P".to_string(),
                taxid: 2,
                name: "Bacteria".to_string(),
                depth: 1,
                children: Vec::new(),
                // Add missing fields for the test
                taxon_id: 2,
                clade_reads: 800,
                direct_reads: 200,
                rank: "P".to_string(),
            },
            TaxonEntry {
                level: 2,
                percentage: 60.0,
                clade_fragments: 600,
                direct_fragments: 100,
                rank_code: "C".to_string(),
                taxid: 3,
                name: "Proteobacteria".to_string(),
                depth: 2,
                children: Vec::new(),
                // Add missing fields for the test
                taxon_id: 3,
                clade_reads: 600,
                direct_reads: 100,
                rank: "C".to_string(),
            },
        ];

        let hierarchy = build_hierarchy_optimized(entries);
        assert_eq!(hierarchy.len(), 1);
        assert_eq!(hierarchy[0].name, "Root");
        assert_eq!(hierarchy[0].children.len(), 1);
        assert_eq!(hierarchy[0].children[0].name, "Bacteria");
        assert_eq!(hierarchy[0].children[0].children.len(), 1);
        assert_eq!(hierarchy[0].children[0].children[0].name, "Proteobacteria");
    }
}
