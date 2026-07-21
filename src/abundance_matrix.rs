use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::OnceLock;

use crate::krk_parser::{KrakenReport, TaxonEntry};

const BUFFER_SIZE: usize = 256 * 1024;
const TAXON_LEVELS: &[(&str, &str)] = &[
    ("D", "domain"),
    ("K", "domain (compatibility alias)"),
    ("P", "phylum"),
    ("C", "class"),
    ("O", "order"),
    ("F", "family"),
    ("G", "genus"),
    ("S", "species"),
];
const UNCLASSIFIED_NAME: &str = "Unclassified";

#[derive(Debug)]
pub enum AbundanceMatrixError {
    IoError(std::io::Error),
    InvalidLevel(String),
}

impl std::fmt::Display for AbundanceMatrixError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::IoError(error) => write!(f, "I/O error: {error}"),
            Self::InvalidLevel(level) => write!(f, "invalid taxonomic level: {level}"),
        }
    }
}

impl Error for AbundanceMatrixError {}

impl From<std::io::Error> for AbundanceMatrixError {
    fn from(error: std::io::Error) -> Self {
        Self::IoError(error)
    }
}

pub type AbundanceResult<T> = Result<T, AbundanceMatrixError>;

#[derive(Clone, Debug)]
struct TaxonMeta {
    name: String,
    rank: String,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MatrixRow {
    pub taxid: u32,
    pub name: String,
    pub rank: String,
    pub values: Vec<f64>,
}

#[derive(Default)]
pub struct AbundanceMatrix {
    taxon_abundances: HashMap<u32, HashMap<String, f64>>,
    taxon_meta: HashMap<u32, TaxonMeta>,
    samples: HashSet<String>,
    level: String,
    sample_totals: HashMap<String, f64>,
    force_include_unclassified: bool,
    cached_sample_names: OnceLock<Vec<String>>,
    cached_rows: OnceLock<Vec<MatrixRow>>,
}

impl AbundanceMatrix {
    pub fn new(level: &str) -> Self {
        Self {
            level: canonical_level(level).to_string(),
            ..Self::default()
        }
    }

    pub fn set_force_include_unclassified(&mut self, include: bool) {
        self.force_include_unclassified = include;
        self.invalidate_cache();
    }

    fn invalidate_cache(&mut self) {
        self.cached_sample_names = OnceLock::new();
        self.cached_rows = OnceLock::new();
    }

    /// Add a sample. When `proportional` is true, values and the threshold are
    /// percentages; otherwise both are absolute read counts.
    pub fn add_sample(
        &mut self,
        report: &KrakenReport,
        sample_name: &str,
        min_abundance: f64,
        proportional: bool,
    ) {
        self.invalidate_cache();
        self.samples.insert(sample_name.to_string());
        let total_reads = report.root.clade_reads as f64
            + report
                .unclassified
                .as_ref()
                .map_or(0.0, |entry| entry.clade_reads as f64);
        self.sample_totals
            .insert(sample_name.to_string(), total_reads);

        if self.force_include_unclassified {
            if let Some(unclassified) = &report.unclassified {
                let abundance =
                    abundance_value(unclassified.clade_reads, total_reads, proportional);
                if abundance >= min_abundance {
                    self.insert(
                        unclassified.taxid,
                        UNCLASSIFIED_NAME,
                        &unclassified.rank,
                        sample_name,
                        abundance,
                    );
                }
            }
        }

        self.process_node(
            &report.root,
            sample_name,
            min_abundance,
            proportional,
            total_reads,
        );
    }

    fn process_node(
        &mut self,
        node: &TaxonEntry,
        sample_name: &str,
        min_abundance: f64,
        proportional: bool,
        total_reads: f64,
    ) {
        if node.rank == self.level {
            let abundance = abundance_value(node.clade_reads, total_reads, proportional);
            if abundance >= min_abundance {
                self.insert(node.taxid, &node.name, &node.rank, sample_name, abundance);
            }
        }

        for child in &node.children {
            self.process_node(child, sample_name, min_abundance, proportional, total_reads);
        }
    }

    fn insert(&mut self, taxid: u32, name: &str, rank: &str, sample_name: &str, abundance: f64) {
        self.taxon_meta.entry(taxid).or_insert_with(|| TaxonMeta {
            name: name.to_string(),
            rank: rank.to_string(),
        });
        self.taxon_abundances
            .entry(taxid)
            .or_default()
            .insert(sample_name.to_string(), abundance);
    }

    pub fn transform_to_proportions(&mut self) {
        self.invalidate_cache();
        for sample_abundances in self.taxon_abundances.values_mut() {
            for (sample, abundance) in sample_abundances {
                if let Some(total) = self.sample_totals.get(sample).filter(|total| **total > 0.0) {
                    *abundance = (*abundance / total) * 100.0;
                }
            }
        }
    }

    pub fn sample_names(&self) -> Vec<String> {
        self.cached_sample_names
            .get_or_init(|| {
                let mut samples: Vec<_> = self.samples.iter().cloned().collect();
                samples.sort();
                samples
            })
            .clone()
    }

    pub fn rows(&self) -> Vec<MatrixRow> {
        self.cached_rows
            .get_or_init(|| {
                let samples = self.sample_names();
                let mut rows: Vec<_> = self
                    .taxon_abundances
                    .iter()
                    .filter_map(|(taxid, abundances)| {
                        let meta = self.taxon_meta.get(taxid)?;
                        Some(MatrixRow {
                            taxid: *taxid,
                            name: meta.name.clone(),
                            rank: meta.rank.clone(),
                            values: samples
                                .iter()
                                .map(|sample| *abundances.get(sample).unwrap_or(&0.0))
                                .collect(),
                        })
                    })
                    .collect();
                rows.sort_by(|left, right| {
                    let left_unclassified = left.name == UNCLASSIFIED_NAME;
                    let right_unclassified = right.name == UNCLASSIFIED_NAME;
                    right_unclassified
                        .cmp(&left_unclassified)
                        .then_with(|| left.name.cmp(&right.name))
                        .then_with(|| left.taxid.cmp(&right.taxid))
                });
                rows
            })
            .clone()
    }

    pub fn write_matrix(&self, output_file: &str) -> AbundanceResult<()> {
        let file = File::create(output_file)?;
        let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);
        let samples = self.sample_names();

        write!(writer, "Taxon")?;
        for sample in &samples {
            write!(writer, "\t{sample}")?;
        }
        writeln!(writer)?;

        for row in self.rows() {
            write!(writer, "{}", row.name)?;
            for abundance in row.values {
                write!(writer, "\t{abundance:.6}")?;
            }
            writeln!(writer)?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Write a MetaPhlAn-style abundance table using taxonomic names as paths.
    pub fn write_mpa(&self, output_file: &str) -> AbundanceResult<()> {
        let file = File::create(output_file)?;
        let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);
        let samples = self.sample_names();

        write!(writer, "#SampleID")?;
        for sample in &samples {
            write!(writer, "\t{sample}")?;
        }
        writeln!(writer)?;

        for row in self.rows() {
            let prefix = mpa_rank_prefix(&row.rank);
            write!(writer, "{prefix}__{}", row.name.replace(' ', "_"))?;
            for abundance in row.values {
                write!(writer, "\t{abundance:.6}")?;
            }
            writeln!(writer)?;
        }

        writer.flush()?;
        Ok(())
    }

    /// Write a simple Krona-compatible TSV (count, taxonomy path fragments).
    pub fn write_krona(&self, output_file: &str) -> AbundanceResult<()> {
        let file = File::create(output_file)?;
        let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);

        // Krona text input expects magnitude first; use the first sample column.
        for row in self.rows() {
            let magnitude = row.values.first().copied().unwrap_or(0.0);
            writeln!(writer, "{magnitude:.6}\t{}", row.name)?;
        }

        writer.flush()?;
        Ok(())
    }
}

fn mpa_rank_prefix(rank: &str) -> &'static str {
    match rank.chars().next().unwrap_or('x') {
        'D' | 'K' => "k",
        'P' => "p",
        'C' => "c",
        'O' => "o",
        'F' => "f",
        'G' => "g",
        'S' => "s",
        _ => "x",
    }
}

fn abundance_value(reads: u64, total_reads: f64, proportional: bool) -> f64 {
    if proportional {
        if total_reads > 0.0 {
            (reads as f64 / total_reads) * 100.0
        } else {
            0.0
        }
    } else {
        reads as f64
    }
}

fn canonical_level(level: &str) -> &str {
    if level == "K" {
        "D"
    } else {
        level
    }
}

pub fn validate_taxonomic_level(level: &str) -> bool {
    TAXON_LEVELS.iter().any(|(code, _)| *code == level)
}

pub fn get_taxonomic_level_name(level: &str) -> Option<&str> {
    TAXON_LEVELS
        .iter()
        .find(|(code, _)| *code == level)
        .map(|(_, name)| *name)
}
