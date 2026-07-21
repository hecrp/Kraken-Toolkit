use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};

use chrono::Local;
use serde_json::{json, Value};

use crate::abundance_matrix::AbundanceMatrix;

const BUFFER_SIZE: usize = 256 * 1024;

#[derive(Debug)]
pub enum BiomError {
    IoError(std::io::Error),
    JsonError(serde_json::Error),
}

impl std::fmt::Display for BiomError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::IoError(error) => write!(f, "I/O error: {error}"),
            Self::JsonError(error) => write!(f, "JSON error: {error}"),
        }
    }
}

impl Error for BiomError {}

impl From<std::io::Error> for BiomError {
    fn from(error: std::io::Error) -> Self {
        Self::IoError(error)
    }
}

impl From<serde_json::Error> for BiomError {
    fn from(error: serde_json::Error) -> Self {
        Self::JsonError(error)
    }
}

pub type BiomResult<T> = Result<T, BiomError>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BiomMatrixType {
    Dense,
    Sparse,
}

#[derive(Debug, PartialEq)]
pub struct BiomRow {
    pub taxid: u32,
    pub name: String,
    pub rank: String,
}

pub struct BiomTable {
    pub dense_data: Option<Vec<Vec<f64>>>,
    pub sparse_data: Option<Vec<[Value; 3]>>,
    pub rows: Vec<BiomRow>,
    pub column_ids: Vec<String>,
    pub date: String,
    pub matrix_type: BiomMatrixType,
}

impl BiomTable {
    pub fn from_abundance_matrix(matrix: &AbundanceMatrix, matrix_type: BiomMatrixType) -> Self {
        let matrix_rows = matrix.rows();
        let column_ids = matrix.sample_names();
        let (dense_data, sparse_data) = match matrix_type {
            BiomMatrixType::Dense => (
                Some(
                    matrix_rows
                        .iter()
                        .map(|row| row.values.clone())
                        .collect::<Vec<_>>(),
                ),
                None,
            ),
            BiomMatrixType::Sparse => {
                let mut sparse = Vec::new();
                for (row_idx, row) in matrix_rows.iter().enumerate() {
                    for (col_idx, value) in row.values.iter().enumerate() {
                        if *value != 0.0 {
                            sparse.push([
                                Value::from(row_idx),
                                Value::from(col_idx),
                                Value::from(*value),
                            ]);
                        }
                    }
                }
                (None, Some(sparse))
            }
        };

        Self {
            dense_data,
            sparse_data,
            rows: matrix_rows
                .into_iter()
                .map(|row| BiomRow {
                    taxid: row.taxid,
                    name: row.name,
                    rank: row.rank,
                })
                .collect(),
            column_ids,
            date: Local::now().to_rfc3339(),
            matrix_type,
        }
    }

    pub fn write_json(&self, output_file: &str) -> BiomResult<()> {
        let file = File::create(output_file)?;
        let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);
        let matrix_type = match self.matrix_type {
            BiomMatrixType::Dense => "dense",
            BiomMatrixType::Sparse => "sparse",
        };
        let data = match self.matrix_type {
            BiomMatrixType::Dense => json!(self.dense_data.clone().unwrap_or_default()),
            BiomMatrixType::Sparse => json!(self.sparse_data.clone().unwrap_or_default()),
        };

        let biom_json = json!({
            "id": "krakenclip",
            "format": "Biological Observation Matrix 1.0.0",
            "format_url": "http://biom-format.org/documentation/format_versions/biom-1.0.html",
            "type": "OTU table",
            "generated_by": format!("KrakenClip {}", env!("CARGO_PKG_VERSION")),
            "date": self.date,
            "matrix_type": matrix_type,
            "matrix_element_type": "float",
            "shape": [self.rows.len(), self.column_ids.len()],
            "data": data,
            "rows": self.rows.iter().map(|row| {
                json!({
                    "id": row.taxid.to_string(),
                    "metadata": {
                        "name": row.name,
                        "rank": row.rank,
                    }
                })
            }).collect::<Vec<_>>(),
            "columns": self.column_ids.iter().map(|id| {
                json!({
                    "id": id,
                    "metadata": {}
                })
            }).collect::<Vec<_>>()
        });

        serde_json::to_writer(&mut writer, &biom_json)?;
        writer.flush()?;
        Ok(())
    }
}
