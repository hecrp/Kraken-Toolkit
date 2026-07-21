use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

const BUFFER_SIZE: usize = 1024 * 1024;

/// Returns true when the path or file contents indicate gzip compression.
pub fn is_gzip_path(path: &str) -> bool {
    Path::new(path)
        .extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}

fn looks_like_gzip(file: &mut File) -> io::Result<bool> {
    let mut magic = [0_u8; 2];
    let read = file.read(&mut magic)?;
    use std::io::Seek;
    file.rewind()?;
    Ok(read == 2 && magic == [0x1f, 0x8b])
}

/// Open a plain or gzip-compressed file for buffered reading.
pub fn open_buf_reader(path: &str) -> io::Result<Box<dyn BufRead + Send>> {
    let mut file = File::open(path)?;
    let gzip = is_gzip_path(path) || looks_like_gzip(&mut file)?;
    if gzip {
        let decoder = MultiGzDecoder::new(file);
        Ok(Box::new(BufReader::with_capacity(BUFFER_SIZE, decoder)))
    } else {
        Ok(Box::new(BufReader::with_capacity(BUFFER_SIZE, file)))
    }
}

/// Open a plain or gzip-compressed file for buffered writing.
pub fn open_buf_writer(path: &str) -> io::Result<Box<dyn Write + Send>> {
    let file = File::create(path)?;
    if is_gzip_path(path) {
        let encoder = GzEncoder::new(file, Compression::fast());
        Ok(Box::new(BufWriter::with_capacity(BUFFER_SIZE, encoder)))
    } else {
        Ok(Box::new(BufWriter::with_capacity(BUFFER_SIZE, file)))
    }
}
