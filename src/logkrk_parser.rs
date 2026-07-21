use memchr::memchr;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Optimized buffer size for efficient file reading
const BUFFER_SIZE: usize = 512 * 1024; // 512KB

/// Constant characters for fast byte-level comparisons
const TAB_CHAR: u8 = b'\t';
const LF_CHAR: u8 = b'\n';

/// Specific error type for Kraken file processing
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

/// Specialized result type for Kraken parsing functions
type KrakenResult<T> = Result<T, KrakenParseError>;

/// Parses a Kraken output file and returns a set of read IDs
/// that match the specified taxids.
///
/// This version is maintained for backward compatibility.
///
/// # Arguments
/// * `kraken_output` - Path to the Kraken output file
/// * `save_taxids` - Set of taxonomic IDs to look for
///
/// # Returns
/// * `Result<HashSet<String>, Box<dyn Error>>` - Set of matching read IDs
#[allow(dead_code)]
pub fn parse_kraken_output(
    kraken_output: &str,
    save_taxids: &HashSet<String>,
) -> Result<HashSet<String>, Box<dyn Error>> {
    let mut taxid_readid_map = HashMap::new();
    let result =
        parse_kraken_output_with_taxids(kraken_output, save_taxids, &mut taxid_readid_map)?;
    Ok(result)
}

/// Parses a Kraken output file and returns a set of read IDs
/// that match the specified taxids. Additionally, it records which
/// taxid corresponds to each read ID in the taxid_readid_map.
///
/// This version is optimized for performance using zero-copy operations when possible.
///
/// # Arguments
/// * `kraken_output` - Path to the Kraken output file
/// * `save_taxids` - Set of taxonomic IDs to look for
/// * `taxid_readid_map` - Map to record the relationship between taxids and read IDs
///
/// # Returns
/// * `Result<HashSet<String>, Box<dyn Error>>` - Set of matching read IDs
///
/// # Implementation Details
/// This function uses several optimization techniques:
/// - Preallocates memory based on expected result sizes
/// - Uses memchr for efficient byte searching
/// - Employs zero-copy string extraction where possible
/// - Reuses buffers to minimize memory allocations
pub fn parse_kraken_output_with_taxids(
    kraken_output: &str,
    save_taxids: &HashSet<String>,
    taxid_readid_map: &mut HashMap<String, HashSet<String>>,
) -> KrakenResult<HashSet<String>> {
    // Estimate expected result size to avoid reallocations
    let estimated_results = save_taxids.len() * 1000;
    let mut save_readids = HashSet::with_capacity(estimated_results);

    // Open the file and create an optimized buffered reader
    let file = File::open(kraken_output)?;
    let mut reader = BufReader::with_capacity(BUFFER_SIZE, file);

    // Reusable buffer to minimize memory allocations
    let mut buffer = Vec::with_capacity(BUFFER_SIZE);

    // Pre-reserve space for tab positions
    let mut tab_positions = Vec::with_capacity(4);

    loop {
        // Clear buffer before next read
        buffer.clear();
        tab_positions.clear();

        // Read until newline character
        let bytes_read = reader.read_until(LF_CHAR, &mut buffer)?;
        if bytes_read == 0 {
            // End of file reached
            break;
        }

        // Remove the newline character for cleaner processing
        if buffer.last() == Some(&LF_CHAR) {
            buffer.pop();
        }

        // Find tab positions using optimized memchr
        let mut pos = 0;
        while let Some(offset) = memchr(TAB_CHAR, &buffer[pos..]) {
            let absolute_pos = pos + offset;
            tab_positions.push(absolute_pos);
            pos = absolute_pos + 1;

            // Exit early once we have enough tabs
            if tab_positions.len() >= 3 {
                break;
            }
        }

        // We need at least 2 tabs to get readid and taxid
        if tab_positions.len() < 2 {
            continue;
        }

        // Extract taxid (field 3)
        let taxid_start = tab_positions[1] + 1;
        let taxid_end = if tab_positions.len() > 2 {
            tab_positions[2]
        } else {
            buffer.len()
        };

        // Convert byte slice to UTF-8 without allocating a new String
        // This is a zero-copy operation that improves performance
        let taxid = std::str::from_utf8(&buffer[taxid_start..taxid_end])?;

        // Check if this taxid is one we're interested in
        if save_taxids.contains(taxid) {
            // Extract readid (field 2)
            let readid_start = tab_positions[0] + 1;
            let readid_end = tab_positions[1];

            // We need to allocate a String because we'll store it in the HashSet
            let readid = std::str::from_utf8(&buffer[readid_start..readid_end])?.to_string();

            // Store the readid in our result set
            save_readids.insert(readid.clone());

            // Record which taxid this readid belongs to
            taxid_readid_map
                .entry(taxid.to_string())
                .or_insert_with(|| {
                    let mut set = HashSet::new();
                    set.reserve(1000); // Pre-allocate space for efficiency
                    set
                })
                .insert(readid);
        }
    }

    Ok(save_readids)
}
