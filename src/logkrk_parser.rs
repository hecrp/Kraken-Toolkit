use memchr::memchr;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::io::BufRead;

use crate::io_util::open_buf_reader;

const BUFFER_SIZE: usize = 512 * 1024;
const TAB_CHAR: u8 = b'\t';
const LF_CHAR: u8 = b'\n';

#[derive(Debug)]
pub enum KrakenParseError {
    IoError(std::io::Error),
    Utf8Error(std::str::Utf8Error),
    #[allow(dead_code)]
    MalformedLine(String),
}

impl std::fmt::Display for KrakenParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::IoError(e) => write!(f, "I/O error: {}", e),
            Self::Utf8Error(e) => write!(f, "UTF-8 decoding error: {}", e),
            Self::MalformedLine(s) => write!(f, "Malformed Kraken log line: {}", s),
        }
    }
}

impl Error for KrakenParseError {}

impl From<std::io::Error> for KrakenParseError {
    fn from(e: std::io::Error) -> Self {
        Self::IoError(e)
    }
}

impl From<std::str::Utf8Error> for KrakenParseError {
    fn from(e: std::str::Utf8Error) -> Self {
        Self::Utf8Error(e)
    }
}

type KrakenResult<T> = Result<T, KrakenParseError>;

/// Lightweight parse that only collects matching read IDs.
pub fn parse_kraken_output(
    kraken_output: &str,
    save_taxids: &HashSet<u32>,
) -> Result<HashSet<String>, Box<dyn Error>> {
    Ok(parse_kraken_log(kraken_output, save_taxids, None)?)
}

/// Parse a Kraken log and optionally populate taxid -> read ID statistics.
pub fn parse_kraken_output_with_taxids(
    kraken_output: &str,
    save_taxids: &HashSet<u32>,
    taxid_readid_map: &mut HashMap<u32, HashSet<String>>,
) -> KrakenResult<HashSet<String>> {
    parse_kraken_log(kraken_output, save_taxids, Some(taxid_readid_map))
}

fn parse_kraken_log(
    kraken_output: &str,
    save_taxids: &HashSet<u32>,
    mut taxid_readid_map: Option<&mut HashMap<u32, HashSet<String>>>,
) -> KrakenResult<HashSet<String>> {
    let estimated_results = save_taxids.len().saturating_mul(1000).max(1024);
    let mut save_readids = HashSet::with_capacity(estimated_results);
    let mut reader = open_buf_reader(kraken_output)?;
    let mut buffer = Vec::with_capacity(BUFFER_SIZE);
    let mut tab_positions = Vec::with_capacity(4);

    loop {
        buffer.clear();
        tab_positions.clear();

        let bytes_read = reader.read_until(LF_CHAR, &mut buffer)?;
        if bytes_read == 0 {
            break;
        }

        if buffer.last() == Some(&LF_CHAR) {
            buffer.pop();
        }

        let mut pos = 0;
        while let Some(offset) = memchr(TAB_CHAR, &buffer[pos..]) {
            let absolute_pos = pos + offset;
            tab_positions.push(absolute_pos);
            pos = absolute_pos + 1;
            if tab_positions.len() >= 3 {
                break;
            }
        }

        if tab_positions.len() < 2 {
            continue;
        }

        let taxid_start = tab_positions[1] + 1;
        let taxid_end = if tab_positions.len() > 2 {
            tab_positions[2]
        } else {
            buffer.len()
        };

        let taxid_str = std::str::from_utf8(&buffer[taxid_start..taxid_end])?;
        let Ok(taxid) = taxid_str.parse::<u32>() else {
            continue;
        };

        if !save_taxids.contains(&taxid) {
            continue;
        }

        let readid_start = tab_positions[0] + 1;
        let readid_end = tab_positions[1];
        let readid = std::str::from_utf8(&buffer[readid_start..readid_end])?.to_string();

        if let Some(map) = taxid_readid_map.as_mut() {
            map.entry(taxid)
                .or_insert_with(|| {
                    let mut set = HashSet::new();
                    set.reserve(256);
                    set
                })
                .insert(readid.clone());
        }

        save_readids.insert(readid);
    }

    Ok(save_readids)
}
