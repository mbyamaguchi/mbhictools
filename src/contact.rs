//! Parallel reading of a sparse contact file (`bin1<TAB>bin2<TAB>score`).
//!
//! Assumed input:
//!   - TSV whose first line is a header (`bin1  bin2  score`)
//!   - upper triangle only (`bin1 <= bin2`)
//!   - `score` is an integer >= 1
//!   - `bin2 - bin1` is within some limit (e.g. 1 Mbp worth of bins)
//!
//! Integer scores let callers accumulate in u64 rather than f64: exact, and cheap
//! enough to add atomically into one shared grid, so no per-thread copies.
//!
//! The file is mmapped, split on line boundaries and parsed with rayon. Records are
//! never retained, so memory stays proportional to the caller's grid.

use std::fs::File;
use std::path::Path;

use memmap2::Mmap;
use rayon::prelude::*;

/// Rough bytes per parallel task, rounded up to a line boundary.
const CHUNK_BYTES: usize = 8 << 20;

/// Result of one pass.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ScanStats {
    /// Lines parsed successfully.
    pub rows: u64,
    /// Lines skipped because they could not be parsed.
    pub malformed: u64,
}

#[derive(Debug)]
pub enum Error {
    Io(std::io::Error),
    /// Too many bad lines, suggesting the file is not in the expected format.
    TooManyMalformed {
        path: String,
        malformed: u64,
        rows: u64,
    },
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Io(e) => write!(f, "I/O error: {e}"),
            Error::TooManyMalformed {
                path,
                malformed,
                rows,
            } => write!(
                f,
                "{path}: too many unparsable lines ({malformed} of {}). \
                 Expected a TSV of `bin1<TAB>bin2<TAB>score`",
                malformed + rows
            ),
        }
    }
}

impl std::error::Error for Error {}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e)
    }
}

/// Scan the file once, calling `visit(bin1, bin2, score)` for every record.
///
/// `visit` runs on many threads at once and not in file order, so it must only do
/// atomic or commutative work (add, min, max).
pub fn visit<F>(path: &Path, visit: F) -> Result<ScanStats, Error>
where
    F: Fn(u32, u32, u32) + Sync,
{
    let file = File::open(path)?;
    // SAFETY: mapped read-only. Truncating the file during the scan would be
    // undefined behaviour; input files are treated as immutable.
    let mmap = unsafe { Mmap::map(&file)? };

    let body = skip_header(&mmap);
    let stats = chunk_ranges(body)
        .into_par_iter()
        .map(|(start, end)| visit_chunk(&body[start..end], &visit))
        .reduce(ScanStats::default, |a, b| ScanStats {
            rows: a.rows + b.rows,
            malformed: a.malformed + b.malformed,
        });

    // A few bad lines are tolerable; mostly bad means the wrong format.
    if stats.malformed > 0 && stats.malformed * 2 > stats.rows {
        return Err(Error::TooManyMalformed {
            path: path.display().to_string(),
            malformed: stats.malformed,
            rows: stats.rows,
        });
    }
    Ok(stats)
}

/// Treat a first line that does not start with a digit as a header.
fn skip_header(buf: &[u8]) -> &[u8] {
    match buf.first() {
        Some(b) if b.is_ascii_digit() => buf,
        _ => match memchr_newline(buf, 0) {
            Some(nl) => &buf[nl + 1..],
            None => &buf[buf.len()..],
        },
    }
}

/// Split `buf` into line-aligned ranges that can be parsed independently.
fn chunk_ranges(buf: &[u8]) -> Vec<(usize, usize)> {
    let mut ranges = Vec::with_capacity(buf.len() / CHUNK_BYTES + 1);
    let mut start = 0;
    while start < buf.len() {
        // Extend to the next newline so no chunk straddles a line.
        let target = (start + CHUNK_BYTES).min(buf.len());
        let end = match memchr_newline(buf, target) {
            Some(nl) => nl + 1,
            None => buf.len(),
        };
        ranges.push((start, end));
        start = end;
    }
    ranges
}

fn visit_chunk<F>(chunk: &[u8], visit: &F) -> ScanStats
where
    F: Fn(u32, u32, u32) + Sync,
{
    let mut stats = ScanStats::default();
    for line in chunk.split(|&b| b == b'\n') {
        let line = trim_cr(line);
        if line.is_empty() {
            continue;
        }
        match parse_line(line) {
            Some((b1, b2, score)) => {
                stats.rows += 1;
                visit(b1, b2, score);
            }
            None => stats.malformed += 1,
        }
    }
    stats
}

/// Parse one `bin1<TAB>bin2<TAB>score` line, ignoring anything after column 3.
fn parse_line(line: &[u8]) -> Option<(u32, u32, u32)> {
    let mut fields = line.split(|&b| b == b'\t');
    let b1 = parse_u32(fields.next()?)?;
    let b2 = parse_u32(fields.next()?)?;
    let score = parse_u32(fields.next()?)?;
    Some((b1, b2, score))
}

/// Parse a decimal unsigned integer. Empty, non-digit or overflowing input is `None`.
fn parse_u32(field: &[u8]) -> Option<u32> {
    if field.is_empty() {
        return None;
    }
    let mut n: u32 = 0;
    for &b in field {
        let d = b.wrapping_sub(b'0');
        if d > 9 {
            return None;
        }
        n = n.checked_mul(10)?.checked_add(d as u32)?;
    }
    Some(n)
}

fn trim_cr(line: &[u8]) -> &[u8] {
    match line.last() {
        Some(b'\r') => &line[..line.len() - 1],
        _ => line,
    }
}

/// Offset of the first `\n` at or after `from`.
fn memchr_newline(buf: &[u8], from: usize) -> Option<usize> {
    buf[from..]
        .iter()
        .position(|&b| b == b'\n')
        .map(|i| from + i)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    fn write_temp(name: &str, content: &[u8]) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(name);
        std::fs::write(&path, content).unwrap();
        path
    }

    #[test]
    fn parses_a_line() {
        assert_eq!(parse_line(b"27\t79\t1"), Some((27, 79, 1)));
        assert_eq!(parse_line(b"62842\t62842\t6"), Some((62842, 62842, 6)));
    }

    #[test]
    fn rejects_malformed_lines() {
        assert_eq!(parse_line(b"27\t79"), None, "missing column");
        assert_eq!(parse_line(b"27\t79\t1.5"), None, "score is not an integer");
        assert_eq!(parse_line(b"27\t79\t"), None, "empty score");
        assert_eq!(parse_line(b"a\tb\tc"), None, "not numeric");
        assert_eq!(parse_line(b"27\t79\t-1"), None, "negative score");
    }

    #[test]
    fn parse_u32_rejects_overflow() {
        assert_eq!(parse_u32(b"4294967295"), Some(u32::MAX));
        assert_eq!(parse_u32(b"4294967296"), None);
    }

    #[test]
    fn visits_every_record_exactly_once() {
        // The score sum shows every record was seen, exactly once.
        let path = write_temp(
            "mbhictools_visit.txt",
            b"bin1\tbin2\tscore\n1\t2\t3\n4\t5\t6\n7\t8\t9\n",
        );
        let sum = AtomicU64::new(0);
        let rows = AtomicU64::new(0);
        let stats = visit(&path, |_, _, s| {
            sum.fetch_add(s as u64, Ordering::Relaxed);
            rows.fetch_add(1, Ordering::Relaxed);
        })
        .unwrap();

        assert_eq!(stats.rows, 3);
        assert_eq!(stats.malformed, 0);
        assert_eq!(rows.load(Ordering::Relaxed), 3);
        assert_eq!(sum.load(Ordering::Relaxed), 18);
    }

    #[test]
    fn handles_crlf_and_missing_trailing_newline() {
        let path = write_temp(
            "mbhictools_crlf.txt",
            b"bin1\tbin2\tscore\r\n1\t2\t3\r\n4\t5\t6",
        );
        let count = AtomicU64::new(0);
        let stats = visit(&path, |_, _, _| {
            count.fetch_add(1, Ordering::Relaxed);
        })
        .unwrap();
        assert_eq!(stats.rows, 2);
        assert_eq!(count.load(Ordering::Relaxed), 2);
    }

    #[test]
    fn accepts_file_without_header() {
        let path = write_temp("mbhictools_noheader.txt", b"1\t2\t3\n");
        let stats = visit(&path, |_, _, _| {}).unwrap();
        assert_eq!(stats.rows, 1, "a first line starting with a digit is data");
    }

    #[test]
    fn rejects_wrong_format() {
        let path = write_temp("mbhictools_wrong.txt", b"a,b,c\nd,e,f\ng,h,i\n");
        assert!(matches!(
            visit(&path, |_, _, _| {}),
            Err(Error::TooManyMalformed { .. })
        ));
    }
}
