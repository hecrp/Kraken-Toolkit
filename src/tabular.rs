use memchr::memchr_iter;

const TAB: u8 = b'\t';

/// Split a byte line into tab-separated field ranges without allocating.
///
/// Returns up to `MAX` field start/end pairs. Trailing newline/CR are ignored.
pub fn split_tabs<const MAX: usize>(line: &[u8]) -> ([&[u8]; MAX], usize) {
    let trimmed = trim_line_ending(line);
    let mut fields = [b"" as &[u8]; MAX];
    let mut count = 0_usize;
    let mut start = 0_usize;

    for pos in memchr_iter(TAB, trimmed) {
        if count >= MAX {
            break;
        }
        fields[count] = &trimmed[start..pos];
        count += 1;
        start = pos + 1;
    }

    if count < MAX {
        fields[count] = &trimmed[start..];
        count += 1;
    }

    (fields, count)
}

/// Collect all tab-separated fields into a reusable vector of slices.
pub fn split_tabs_all<'a>(line: &'a [u8], out: &mut Vec<&'a [u8]>) {
    out.clear();
    let trimmed = trim_line_ending(line);
    let mut start = 0_usize;
    for pos in memchr_iter(TAB, trimmed) {
        out.push(&trimmed[start..pos]);
        start = pos + 1;
    }
    out.push(&trimmed[start..]);
}

#[inline]
pub fn trim_line_ending(line: &[u8]) -> &[u8] {
    let mut end = line.len();
    while end > 0 && (line[end - 1] == b'\n' || line[end - 1] == b'\r') {
        end -= 1;
    }
    &line[..end]
}

#[inline]
pub fn parse_u64_field(bytes: &[u8]) -> Option<u64> {
    std::str::from_utf8(bytes).ok()?.parse().ok()
}

#[inline]
pub fn parse_u32_field(bytes: &[u8]) -> Option<u32> {
    std::str::from_utf8(bytes).ok()?.parse().ok()
}

#[inline]
pub fn parse_f64_field(bytes: &[u8]) -> Option<f64> {
    std::str::from_utf8(bytes).ok()?.parse().ok()
}

#[inline]
pub fn field_str(bytes: &[u8]) -> Option<&str> {
    std::str::from_utf8(bytes).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splits_fixed_fields() {
        let (fields, count) = split_tabs::<4>(b"a\tb\tc\td\n");
        assert_eq!(count, 4);
        assert_eq!(fields[0], b"a");
        assert_eq!(fields[3], b"d");
    }

    #[test]
    fn splits_all_fields() {
        let mut out = Vec::new();
        split_tabs_all(b"one\ttwo\tthree\r\n", &mut out);
        assert_eq!(
            out,
            vec![b"one".as_ref(), b"two".as_ref(), b"three".as_ref()]
        );
    }
}
