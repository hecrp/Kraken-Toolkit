use memchr::memchr_iter;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Read;
use std::io::Write;
use std::path::Path;
use std::time::Instant;

const BUFFER_SIZE: usize = 512 * 1024;
const TAB_CHAR: u8 = b'\t';
const SPACE_CHAR: u8 = b' ';
const NEWLINE_CHAR: u8 = b'\n';

/// Node in the taxonomic tree produced from a Kraken2 report.
#[derive(Debug, Clone)]
pub struct TaxonEntry {
    pub percentage: f32,
    pub clade_reads: u64,
    pub direct_reads: u64,
    pub rank: String,
    pub taxid: u32,
    pub name: String,
    pub depth: usize,
    pub children: Vec<TaxonEntry>,
}

impl TaxonEntry {
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
            rank,
            taxid,
            name,
            depth,
            children: Vec::new(),
        }
    }

    #[inline]
    pub fn level(&self) -> usize {
        self.depth
    }

    #[inline]
    pub fn taxon_id(&self) -> u64 {
        self.taxid as u64
    }

    #[inline]
    pub fn clade_fragments(&self) -> u64 {
        self.clade_reads
    }

    #[inline]
    pub fn direct_fragments(&self) -> u64 {
        self.direct_reads
    }

    #[inline]
    pub fn rank_code(&self) -> &str {
        &self.rank
    }
}

/// Indexed taxonomic hierarchy for O(1) parent/child lookups.
#[derive(Debug, Default, Clone)]
pub struct TaxonIndex {
    parents: HashMap<u32, u32>,
    children: HashMap<u32, Vec<u32>>,
}

impl TaxonIndex {
    pub fn from_tree(root: &TaxonEntry, unclassified: Option<&TaxonEntry>) -> Self {
        let mut index = Self::default();
        if let Some(node) = unclassified {
            Self::index_node(node, None, &mut index);
        }
        Self::index_node(root, None, &mut index);
        index
    }

    fn index_node(node: &TaxonEntry, parent: Option<u32>, index: &mut Self) {
        if let Some(parent_id) = parent {
            index.parents.insert(node.taxid, parent_id);
        }
        if !node.children.is_empty() {
            index.children.insert(
                node.taxid,
                node.children.iter().map(|child| child.taxid).collect(),
            );
        }
        for child in &node.children {
            Self::index_node(child, Some(node.taxid), index);
        }
    }

    pub fn all_descendants(&self, taxid: u32) -> HashSet<u32> {
        let mut result = HashSet::new();
        let mut stack = Vec::new();
        if let Some(children) = self.children.get(&taxid) {
            stack.extend(children.iter().copied());
        }
        while let Some(current) = stack.pop() {
            if result.insert(current) {
                if let Some(children) = self.children.get(&current) {
                    stack.extend(children.iter().copied());
                }
            }
        }
        result
    }

    pub fn all_ancestors(&self, taxid: u32) -> HashSet<u32> {
        let mut result = HashSet::new();
        let mut current = taxid;
        while let Some(&parent) = self.parents.get(&current) {
            result.insert(parent);
            current = parent;
        }
        result
    }
}

/// Complete Kraken report with hierarchy and lookup index.
#[derive(Debug)]
pub struct KrakenReport {
    pub root: TaxonEntry,
    pub unclassified: Option<TaxonEntry>,
    pub index: TaxonIndex,
}

impl Default for KrakenReport {
    fn default() -> Self {
        let root = TaxonEntry::new(0.0, 0, 0, "R".to_string(), 1, "root".to_string(), 0);
        let index = TaxonIndex::from_tree(&root, None);
        Self {
            root,
            unclassified: None,
            index,
        }
    }
}

/// Cache of common taxonomic rank codes.
pub struct StringCache {
    rank_codes: Vec<String>,
}

impl StringCache {
    fn new() -> Self {
        let mut rank_codes = Vec::with_capacity(16);
        for code in ["R", "D", "P", "C", "O", "F", "G", "S", "U"].iter() {
            rank_codes.push((*code).to_string());
        }
        Self { rank_codes }
    }

    fn get_rank_code(&mut self, code: &str) -> String {
        for existing in &self.rank_codes {
            if existing == code {
                return existing.clone();
            }
        }
        let code_string = code.to_string();
        self.rank_codes.push(code_string.clone());
        code_string
    }
}

#[inline(always)]
pub fn parse_line(
    line: &[u8],
    string_cache: &mut StringCache,
    line_number: Option<usize>,
) -> Option<TaxonEntry> {
    let mut tab_positions = [0_usize; 5];
    let mut tab_count = 0_usize;
    for pos in memchr_iter(TAB_CHAR, line) {
        if tab_count == 5 {
            break;
        }
        tab_positions[tab_count] = pos;
        tab_count += 1;
    }
    if tab_count < 5 {
        if let Some(line_num) = line_number {
            eprintln!(
                "Warning: Line {} does not have enough tab-separated fields (expected at least 5, found {})",
                line_num, tab_count
            );
        }
        return None;
    }

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

    let name_start = field_starts[5];
    let mut level = 0;
    let mut i = name_start;
    while i + 1 < line.len() && line[i] == SPACE_CHAR && line[i + 1] == SPACE_CHAR {
        level += 1;
        i += 2;
    }

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

    let clade_bytes = &line[field_starts[1]..field_ends[1]];
    let clade_reads = match std::str::from_utf8(clade_bytes) {
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
    let direct_reads = match std::str::from_utf8(direct_bytes) {
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
    let rank = string_cache.get_rank_code(rank_str);

    let taxon_bytes = &line[field_starts[4]..field_ends[4]];
    let taxid = match std::str::from_utf8(taxon_bytes) {
        Ok(s) => s.parse::<u32>().unwrap_or_else(|_| {
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

    Some(TaxonEntry {
        percentage,
        clade_reads,
        direct_reads,
        rank,
        taxid,
        name: name.to_string(),
        depth: level,
        children: Vec::new(),
    })
}

pub fn build_hierarchy_optimized(entries: Vec<TaxonEntry>) -> Vec<TaxonEntry> {
    if entries.is_empty() {
        return Vec::new();
    }

    let mut hierarchy = Vec::with_capacity(entries.len() / 8 + 1);
    let mut stack: Vec<TaxonEntry> = Vec::with_capacity(20);

    for entry in entries {
        while !stack.is_empty() && stack.last().unwrap().depth >= entry.depth {
            let popped = stack.pop().unwrap();
            if let Some(parent) = stack.last_mut() {
                if parent.children.is_empty() {
                    parent.children.reserve(4);
                }
                parent.children.push(popped);
            } else {
                hierarchy.push(popped);
            }
        }
        stack.push(entry);
    }

    while let Some(entry) = stack.pop() {
        if let Some(parent) = stack.last_mut() {
            parent.children.push(entry);
        } else {
            hierarchy.push(entry);
        }
    }

    hierarchy
}

struct OptimizedBuffer {
    buffer: Box<[u8]>,
    pos: usize,
    cap: usize,
    file: File,
}

impl OptimizedBuffer {
    fn new(file: File) -> Self {
        Self {
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

        if self.pos >= self.cap {
            match self.fill_buffer() {
                Ok(0) => return Ok(false),
                Ok(_) => {}
                Err(e) => {
                    return Err(std::io::Error::new(
                        e.kind(),
                        format!("Failed to read from file: {}", e),
                    ))
                }
            }
        }

        loop {
            let mut i = self.pos;
            while i < self.cap {
                if self.buffer[i] == NEWLINE_CHAR {
                    line_buffer.extend_from_slice(&self.buffer[self.pos..i]);
                    self.pos = i + 1;
                    return Ok(true);
                }
                i += 1;
            }

            line_buffer.extend_from_slice(&self.buffer[self.pos..self.cap]);

            match self.fill_buffer() {
                Ok(0) => return Ok(!line_buffer.is_empty()),
                Ok(_) => {}
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

pub fn parse_kraken2_report(file_path: &str) -> Result<(KrakenReport, f64), std::io::Error> {
    let file = File::open(file_path)?;
    let file_size = file.metadata().map(|m| m.len() as usize).unwrap_or(0);
    let mut buffer = OptimizedBuffer::new(file);
    let start_time = Instant::now();
    let mut string_cache = StringCache::new();
    let mut line_buffer = Vec::with_capacity(1024);
    let mut entries = Vec::with_capacity(file_size / 50);

    let mut line_number = 1;
    while buffer.read_line(&mut line_buffer)? {
        if let Some(entry) = parse_line(&line_buffer, &mut string_cache, Some(line_number)) {
            entries.push(entry);
        }
        line_number += 1;
    }

    let mut hierarchy = build_hierarchy_optimized(entries);

    let unclassified = if !hierarchy.is_empty()
        && hierarchy[0].depth == 0
        && hierarchy[0].name == "unclassified"
    {
        Some(hierarchy.remove(0))
    } else {
        None
    };

    let root = if !hierarchy.is_empty() {
        hierarchy.remove(0)
    } else {
        TaxonEntry::new(0.0, 0, 0, "R".to_string(), 1, "root".to_string(), 0)
    };

    let index = TaxonIndex::from_tree(&root, unclassified.as_ref());
    let duration = start_time.elapsed().as_secs_f64();
    Ok((
        KrakenReport {
            unclassified,
            root,
            index,
        },
        duration,
    ))
}

pub fn write_json_report(report: &KrakenReport, output_path: &str) -> std::io::Result<()> {
    fn node_to_json(node: &TaxonEntry) -> serde_json::Value {
        let mut json = serde_json::json!({
            "name": node.name,
            "taxid": node.taxid,
            "rank": node.rank,
            "percentage": node.percentage,
            "clade_reads": node.clade_reads,
            "direct_reads": node.direct_reads,
            "level": node.depth,
            "children": []
        });
        let children = node.children.iter().map(node_to_json).collect::<Vec<_>>();
        json["children"] = serde_json::Value::Array(children);
        json
    }

    let mut json = serde_json::json!({
        "root": node_to_json(&report.root)
    });

    if let Some(ref unclassified) = report.unclassified {
        json["unclassified"] = node_to_json(unclassified);
    }

    let file = std::fs::File::create(Path::new(output_path))?;
    let mut writer = std::io::BufWriter::new(file);
    serde_json::to_writer(&mut writer, &json)?;
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
        assert_eq!(entry.depth, 1);
        assert_eq!(entry.percentage, 50.0);
        assert_eq!(entry.clade_reads, 1000);
        assert_eq!(entry.direct_reads, 500);
        assert_eq!(entry.rank, "P");
        assert_eq!(entry.taxid, 123);
        assert_eq!(entry.name, "Bacteria");
    }

    #[test]
    fn test_build_hierarchy() {
        let entries = vec![
            TaxonEntry::new(100.0, 1000, 0, "D".to_string(), 1, "Root".to_string(), 0),
            TaxonEntry::new(
                80.0,
                800,
                200,
                "P".to_string(),
                2,
                "Bacteria".to_string(),
                1,
            ),
            TaxonEntry::new(
                60.0,
                600,
                100,
                "C".to_string(),
                3,
                "Proteobacteria".to_string(),
                2,
            ),
        ];

        let hierarchy = build_hierarchy_optimized(entries);
        assert_eq!(hierarchy.len(), 1);
        assert_eq!(hierarchy[0].name, "Root");
        assert_eq!(hierarchy[0].children.len(), 1);
        assert_eq!(hierarchy[0].children[0].name, "Bacteria");
        assert_eq!(hierarchy[0].children[0].children.len(), 1);
        assert_eq!(hierarchy[0].children[0].children[0].name, "Proteobacteria");
    }

    #[test]
    fn taxon_index_resolves_ancestors_and_descendants() {
        let child = TaxonEntry::new(50.0, 50, 50, "S".to_string(), 3, "alpha".to_string(), 2);
        let mut bacteria =
            TaxonEntry::new(90.0, 90, 0, "D".to_string(), 2, "Bacteria".to_string(), 1);
        bacteria.children.push(child);
        let mut root = TaxonEntry::new(100.0, 100, 0, "R".to_string(), 1, "root".to_string(), 0);
        root.children.push(bacteria);
        let index = TaxonIndex::from_tree(&root, None);
        assert_eq!(index.all_descendants(2), HashSet::from([3]));
        assert_eq!(index.all_ancestors(3), HashSet::from([2, 1]));
    }
}
