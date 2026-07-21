use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use rayon::prelude::*;

use crate::bracken::{parse_bracken_file, BrackenRecord};
use crate::io_util::open_buf_reader;
use crate::krk_parser::{collect_flat_taxons, FlatTaxon};
use crate::tabular::split_tabs_all;

const BUFFER_SIZE: usize = 256 * 1024;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlphaMetric {
    Shannon,
    BergerParker,
    Simpson,
    InverseSimpson,
    Fisher,
}

impl AlphaMetric {
    pub fn parse(name: &str) -> Option<Self> {
        match name.to_ascii_lowercase().as_str() {
            "sh" | "shannon" => Some(Self::Shannon),
            "bp" | "berger-parker" | "berger_parker" => Some(Self::BergerParker),
            "si" | "simpson" => Some(Self::Simpson),
            "isi" | "inverse-simpson" | "inverse_simpson" => Some(Self::InverseSimpson),
            "fisher" | "fisher-alpha" | "fisher_alpha" => Some(Self::Fisher),
            _ => None,
        }
    }
}

#[derive(Debug, Default, Clone)]
struct AlphaAccumulators {
    total: f64,
    richness: u64,
    shannon_term: f64,
    max_count: f64,
    simpson_num: f64,
}

impl AlphaAccumulators {
    fn observe(&mut self, count: f64) {
        if count <= 0.0 {
            return;
        }
        self.total += count;
        self.richness += 1;
        self.shannon_term += count * count.ln();
        self.max_count = self.max_count.max(count);
        self.simpson_num += count * (count - 1.0);
    }

    fn finish(self, metric: AlphaMetric) -> f64 {
        if self.total <= 0.0 {
            return 0.0;
        }
        match metric {
            AlphaMetric::Shannon => {
                let h = self.shannon_term / self.total - self.total.ln();
                -h
            }
            AlphaMetric::BergerParker => self.max_count / self.total,
            AlphaMetric::Simpson => {
                let denom = self.total * (self.total - 1.0);
                if denom <= 0.0 {
                    0.0
                } else {
                    1.0 - (self.simpson_num / denom)
                }
            }
            AlphaMetric::InverseSimpson => {
                let denom = self.total * (self.total - 1.0);
                if denom <= 0.0 {
                    0.0
                } else {
                    let d = self.simpson_num / denom;
                    if d <= 0.0 {
                        0.0
                    } else {
                        1.0 / d
                    }
                }
            }
            AlphaMetric::Fisher => fisher_alpha(self.total, self.richness as f64),
        }
    }
}

/// Fisher's alpha via bounded Newton / bisection on S = alpha * ln(1 + N/alpha).
fn fisher_alpha(n: f64, s: f64) -> f64 {
    if n <= 0.0 || s <= 0.0 {
        return 0.0;
    }
    if s >= n {
        return s;
    }

    let f = |alpha: f64| alpha * (1.0 + n / alpha).ln() - s;
    let mut lo = 1e-12;
    let mut hi = n.max(s);
    while f(hi) < 0.0 {
        hi *= 2.0;
        if hi > 1e18 {
            return hi;
        }
    }

    let mut alpha = ((lo + hi) / 2.0).max(1e-12);
    for _ in 0..64 {
        let value = f(alpha);
        if value.abs() < 1e-10 {
            break;
        }
        let derivative = (1.0 + n / alpha).ln() - (n / (alpha + n));
        if derivative.abs() < 1e-18 {
            break;
        }
        let next = alpha - value / derivative;
        if next <= lo || next >= hi || !next.is_finite() {
            if value > 0.0 {
                hi = alpha;
            } else {
                lo = alpha;
            }
            alpha = (lo + hi) / 2.0;
        } else {
            if value > 0.0 {
                hi = alpha;
            } else {
                lo = alpha;
            }
            alpha = next;
        }
    }
    alpha
}

pub fn alpha_from_counts<I>(counts: I, metric: AlphaMetric) -> f64
where
    I: IntoIterator<Item = f64>,
{
    let mut acc = AlphaAccumulators::default();
    for count in counts {
        acc.observe(count);
    }
    acc.finish(metric)
}

pub fn alpha_from_bracken(records: &[BrackenRecord], metric: AlphaMetric) -> f64 {
    alpha_from_counts(
        records.iter().map(|record| record.new_est_reads as f64),
        metric,
    )
}

#[derive(Debug, Clone)]
pub struct SampleCounts {
    pub name: String,
    pub counts: HashMap<u32, u64>,
    pub total: u64,
}

pub fn load_sample_counts_bracken(
    path: &str,
    sample_name: &str,
) -> Result<SampleCounts, Box<dyn Error>> {
    let records = parse_bracken_file(path)?;
    let mut counts = HashMap::new();
    let mut total = 0_u64;
    for record in records {
        counts.insert(record.taxonomy_id, record.new_est_reads);
        total = total.saturating_add(record.new_est_reads);
    }
    Ok(SampleCounts {
        name: sample_name.to_string(),
        counts,
        total,
    })
}

pub fn load_sample_counts_kreport(
    path: &str,
    sample_name: &str,
    level: Option<&str>,
) -> Result<SampleCounts, Box<dyn Error>> {
    let rows = collect_flat_taxons(path, level)?;
    let mut counts = HashMap::new();
    let mut total = 0_u64;
    for row in rows {
        // Prefer direct reads for diversity comparisons when available.
        let value = if row.direct_reads > 0 {
            row.direct_reads
        } else {
            row.clade_reads
        };
        counts.insert(row.taxid, value);
        total = total.saturating_add(value);
    }
    Ok(SampleCounts {
        name: sample_name.to_string(),
        counts,
        total,
    })
}

pub fn load_sample_counts_tsv(
    path: &str,
    sample_name: &str,
    category_col: usize,
    count_col: usize,
) -> Result<SampleCounts, Box<dyn Error>> {
    let mut reader = open_buf_reader(path)?;
    let mut line = Vec::new();
    let mut counts = HashMap::new();
    let mut total = 0_u64;
    let mut next_id = 1_u32;
    let mut labels: HashMap<String, u32> = HashMap::new();

    loop {
        line.clear();
        let bytes = std::io::BufRead::read_until(&mut reader, b'\n', &mut line)?;
        if bytes == 0 {
            break;
        }
        if line.iter().all(u8::is_ascii_whitespace) {
            continue;
        }
        let mut field_bufs = Vec::new();
        split_tabs_all(&line, &mut field_bufs);
        let owned_fields: Vec<String> = field_bufs
            .iter()
            .map(|field| String::from_utf8_lossy(field).into_owned())
            .collect();
        if owned_fields.len() <= category_col.max(count_col) {
            continue;
        }
        let category = owned_fields[category_col].clone();
        if category.is_empty() || category.starts_with('#') || category.eq_ignore_ascii_case("name")
        {
            continue;
        }
        let count = owned_fields[count_col]
            .parse::<u64>()
            .ok()
            .or_else(|| {
                owned_fields[count_col]
                    .parse::<f64>()
                    .ok()
                    .map(|v| v.round() as u64)
            })
            .unwrap_or(0);
        if count == 0 {
            continue;
        }
        let taxid = if let Ok(id) = category.parse::<u32>() {
            id
        } else {
            *labels.entry(category).or_insert_with(|| {
                let id = next_id;
                next_id += 1;
                id
            })
        };
        *counts.entry(taxid).or_insert(0) += count;
        total = total.saturating_add(count);
    }

    Ok(SampleCounts {
        name: sample_name.to_string(),
        counts,
        total,
    })
}

pub fn bray_curtis(left: &SampleCounts, right: &SampleCounts) -> f64 {
    if left.total == 0 && right.total == 0 {
        return 0.0;
    }
    let mut shared = 0_u64;
    for (taxid, left_count) in &left.counts {
        if let Some(right_count) = right.counts.get(taxid) {
            shared += (*left_count).min(*right_count);
        }
    }
    let denom = (left.total + right.total) as f64;
    if denom <= 0.0 {
        0.0
    } else {
        1.0 - ((2.0 * shared as f64) / denom)
    }
}

pub fn beta_matrix(samples: &[SampleCounts]) -> Vec<Vec<f64>> {
    let n = samples.len();
    let mut matrix = vec![vec![0.0; n]; n];
    let pairs: Vec<(usize, usize)> = (0..n)
        .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
        .collect();
    let values: Vec<(usize, usize, f64)> = pairs
        .par_iter()
        .map(|&(i, j)| (i, j, bray_curtis(&samples[i], &samples[j])))
        .collect();
    for (i, j, value) in values {
        matrix[i][j] = value;
        matrix[j][i] = value;
    }
    matrix
}

pub fn write_beta_matrix(
    samples: &[SampleCounts],
    matrix: &[Vec<f64>],
    output_file: &str,
) -> Result<(), Box<dyn Error>> {
    let file = File::create(output_file)?;
    let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);
    write!(writer, "sample")?;
    for sample in samples {
        write!(writer, "\t{}", sample.name)?;
    }
    writeln!(writer)?;
    for (i, sample) in samples.iter().enumerate() {
        write!(writer, "{}", sample.name)?;
        for value in &matrix[i] {
            write!(writer, "\t{value:.3}")?;
        }
        writeln!(writer)?;
    }
    writer.flush()?;
    Ok(())
}

pub fn default_sample_name(path: &str, index: usize) -> String {
    Path::new(path)
        .file_stem()
        .and_then(|name| name.to_str())
        .map(|name| name.to_string())
        .unwrap_or_else(|| format!("sample_{}", index + 1))
}

#[allow(dead_code)]
pub fn flat_taxon_count(row: &FlatTaxon) -> u64 {
    if row.direct_reads > 0 {
        row.direct_reads
    } else {
        row.clade_reads
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shannon_of_even_two_species() {
        let value = alpha_from_counts([50.0, 50.0], AlphaMetric::Shannon);
        assert!((value - 2.0_f64.ln()).abs() < 1e-9);
    }

    #[test]
    fn bray_curtis_identical_and_disjoint() {
        let a = SampleCounts {
            name: "a".into(),
            counts: HashMap::from([(1, 10), (2, 5)]),
            total: 15,
        };
        let identical = a.clone();
        let b = SampleCounts {
            name: "b".into(),
            counts: HashMap::from([(3, 15)]),
            total: 15,
        };
        assert!((bray_curtis(&a, &identical) - 0.0).abs() < 1e-12);
        assert!((bray_curtis(&a, &b) - 1.0).abs() < 1e-12);
    }
}
