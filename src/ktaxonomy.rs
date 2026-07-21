use std::collections::HashMap;
use std::error::Error;
use std::io::BufRead;

use crate::io_util::open_buf_reader;
use crate::tabular::parse_u32_field;

#[derive(Debug, Clone)]
pub struct KTaxonomyNode {
    pub taxid: u32,
    pub parent: u32,
    pub rank: String,
    pub level: usize,
    pub name: String,
    pub children: Vec<u32>,
}

#[derive(Debug, Default)]
pub struct KTaxonomy {
    pub nodes: HashMap<u32, KTaxonomyNode>,
}

impl KTaxonomy {
    /// Load condensed taxonomy produced by make_ktaxonomy.py:
    /// `taxid\t|\tparent\t|\trank\t|\tlevel\t|\tname`
    pub fn from_file(path: &str) -> Result<Self, Box<dyn Error>> {
        let mut reader = open_buf_reader(path)?;
        let mut line = Vec::new();
        let mut taxonomy = Self::default();

        loop {
            line.clear();
            let bytes = reader.read_until(b'\n', &mut line)?;
            if bytes == 0 {
                break;
            }
            if line.iter().all(u8::is_ascii_whitespace) || line.starts_with(b"#") {
                continue;
            }

            let text = String::from_utf8_lossy(&line).trim().to_string();
            let parts: Vec<String> = if text.contains("\t|\t") {
                text.split("\t|\t").map(|part| part.to_string()).collect()
            } else {
                text.split('\t').map(|part| part.to_string()).collect()
            };
            if parts.len() < 5 {
                continue;
            }

            let taxid = parts[0].parse::<u32>()?;
            let parent = parts[1].parse::<u32>().unwrap_or(taxid);
            let rank = parts[2].clone();
            let level = parts[3].parse::<usize>().unwrap_or(0);
            let name = parts[4].clone();
            taxonomy.nodes.insert(
                taxid,
                KTaxonomyNode {
                    taxid,
                    parent,
                    rank,
                    level,
                    name,
                    children: Vec::new(),
                },
            );
        }

        let parents: Vec<(u32, u32)> = taxonomy
            .nodes
            .values()
            .map(|node| (node.taxid, node.parent))
            .collect();
        for (taxid, parent) in parents {
            if taxid != parent {
                if let Some(parent_node) = taxonomy.nodes.get_mut(&parent) {
                    parent_node.children.push(taxid);
                }
            }
        }

        Ok(taxonomy)
    }

    pub fn get(&self, taxid: u32) -> Option<&KTaxonomyNode> {
        self.nodes.get(&taxid)
    }
}

#[inline]
pub fn parse_taxid_field(bytes: &[u8]) -> Option<u32> {
    parse_u32_field(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn loads_pipe_separated_taxonomy() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(
            file,
            "1\t|\t1\t|\tR\t|\t0\t|\troot\n2\t|\t1\t|\tD\t|\t1\t|\tBacteria\n3\t|\t2\t|\tS\t|\t2\t|\tSpecies alpha"
        )
        .unwrap();
        let tax = KTaxonomy::from_file(file.path().to_str().unwrap()).unwrap();
        assert_eq!(tax.get(3).unwrap().parent, 2);
        assert_eq!(tax.get(1).unwrap().children, vec![2]);
    }
}
