//! A chromosome length table (`chr<TAB>bp`) and its mapping to global bins.
//!
//! Bins in a contact file are 1-based indices over the concatenated genome, not bp.
//! Each chromosome spans `ceil(bp / bin_size)` of them.
//!
//! Knowing the mapping allows two things: labelling axes in bp, and dropping contacts
//! that cross a chromosome boundary, whose "distance" is meaningless once the genome
//! has been laid out in a line.

use std::path::Path;

/// One chromosome and the global bins it spans, [start_bin, end_bin] inclusive.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chrom {
    pub name: String,
    pub length_bp: u64,
    pub bins: u32,
    pub start_bin: u32,
    pub end_bin: u32,
}

/// A chromosome table that can map a bin back to its chromosome.
#[derive(Debug, Clone)]
pub struct ChromIndex {
    chroms: Vec<Chrom>,
    /// First bin of each chromosome, ascending, for binary search.
    starts: Vec<u32>,
    bin_size: u32,
}

#[derive(Debug)]
pub enum Error {
    Io(std::io::Error),
    Parse {
        path: String,
        line: usize,
        reason: String,
    },
    Empty(String),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Io(e) => write!(f, "I/O error: {e}"),
            Error::Parse { path, line, reason } => {
                write!(f, "{path}:{line}: {reason} (expected `name<TAB>length_bp`)")
            }
            Error::Empty(path) => write!(f, "{path}: no chromosomes found"),
        }
    }
}

impl std::error::Error for Error {}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e)
    }
}

impl ChromIndex {
    /// Read a length table and assign bin ranges at `bin_size` bp per bin.
    pub fn load(path: &Path, bin_size: u32) -> Result<Self, Error> {
        let text = std::fs::read_to_string(path)?;
        Self::parse(&text, bin_size, &path.display().to_string())
    }

    fn parse(text: &str, bin_size: u32, path: &str) -> Result<Self, Error> {
        assert!(bin_size > 0, "bin_size must be positive");
        let mut chroms = Vec::new();
        let mut next_start: u32 = 1; // global bins are 1-based

        for (i, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let (name, bp) = line.split_once('\t').ok_or_else(|| Error::Parse {
                path: path.to_string(),
                line: i + 1,
                reason: "no tab-separated second column".into(),
            })?;
            let length_bp: u64 = bp.trim().parse().map_err(|_| Error::Parse {
                path: path.to_string(),
                line: i + 1,
                reason: format!("length `{}` is not an integer", bp.trim()),
            })?;
            if length_bp == 0 {
                return Err(Error::Parse {
                    path: path.to_string(),
                    line: i + 1,
                    reason: "length is zero".into(),
                });
            }

            let bins = length_bp.div_ceil(bin_size as u64) as u32;
            chroms.push(Chrom {
                name: name.trim().to_string(),
                length_bp,
                bins,
                start_bin: next_start,
                end_bin: next_start + bins - 1,
            });
            next_start += bins;
        }

        if chroms.is_empty() {
            return Err(Error::Empty(path.to_string()));
        }
        let starts = chroms.iter().map(|c| c.start_bin).collect();
        Ok(ChromIndex {
            chroms,
            starts,
            bin_size,
        })
    }

    /// Total length in bp.
    pub fn total_bp(path: &Path) -> Result<u64, Error> {
        // Independent of bin_size, so read at a nominal 1 bp/bin.
        let idx = Self::load(path, 1)?;
        Ok(idx.chroms.iter().map(|c| c.length_bp).sum())
    }

    pub fn chroms(&self) -> &[Chrom] {
        &self.chroms
    }

    pub fn bin_size(&self) -> u32 {
        self.bin_size
    }

    /// Bins in the whole genome (the last chromosome's `end_bin`).
    pub fn total_bins(&self) -> u32 {
        self.chroms.last().map_or(0, |c| c.end_bin)
    }

    /// Index of the chromosome containing `bin`, or `None` if out of range.
    pub fn chrom_of(&self, bin: u32) -> Option<usize> {
        if bin < self.starts[0] || bin > self.total_bins() {
            return None;
        }
        // starts ascends, so the answer is the last chromosome starting at or below bin.
        Some(match self.starts.binary_search(&bin) {
            Ok(i) => i,
            Err(i) => i - 1,
        })
    }

    /// Are both bins on the same chromosome? False if either is out of range.
    pub fn same_chrom(&self, bin1: u32, bin2: u32) -> bool {
        match (self.chrom_of(bin1), self.chrom_of(bin2)) {
            (Some(a), Some(b)) => a == b,
            _ => false,
        }
    }

    /// Global bins where one chromosome meets the next, for drawing dividers.
    pub fn boundaries(&self) -> Vec<u32> {
        self.chroms.iter().skip(1).map(|c| c.start_bin).collect()
    }
}

/// Estimate `bin_size` from total bp / max bin.
///
/// Bins are a 1-based running index, so the bins covering the genome are about the
/// max bin. Rounding absorbs the error from per-chromosome ceilings (200.05 -> 200).
pub fn estimate_bin_size(total_bp: u64, max_bin: u32) -> u32 {
    assert!(max_bin > 0, "max_bin must be positive");
    let raw = total_bp as f64 / max_bin as f64;
    // Rounding to a step one decade below `raw` lands on the round sizes actually
    // used in practice (100, 200, 500, 1000, ...).
    let step = 10f64.powi((raw.log10().floor() as i32 - 1).max(0));
    ((raw / step).round() * step) as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    const POMBE: &str = "I\t5579133\nII\t4539804\nIII\t2452883\n";

    fn pombe(bin_size: u32) -> ChromIndex {
        ChromIndex::parse(POMBE, bin_size, "test").unwrap()
    }

    #[test]
    fn assigns_contiguous_bin_ranges() {
        let idx = pombe(200);
        let c = idx.chroms();
        assert_eq!(c.len(), 3);
        // ceil(5579133/200) = 27896, ceil(4539804/200) = 22700, ceil(2452883/200) = 12265
        assert_eq!((c[0].bins, c[0].start_bin, c[0].end_bin), (27896, 1, 27896));
        assert_eq!(
            (c[1].bins, c[1].start_bin, c[1].end_bin),
            (22700, 27897, 50596)
        );
        assert_eq!(
            (c[2].bins, c[2].start_bin, c[2].end_bin),
            (12265, 50597, 62861)
        );
        assert_eq!(idx.total_bins(), 62861);
    }

    #[test]
    fn maps_bins_back_to_chromosomes() {
        let idx = pombe(200);
        assert_eq!(idx.chrom_of(1), Some(0));
        assert_eq!(idx.chrom_of(27896), Some(0), "end of I");
        assert_eq!(idx.chrom_of(27897), Some(1), "start of II");
        assert_eq!(idx.chrom_of(50596), Some(1));
        assert_eq!(idx.chrom_of(50597), Some(2));
        assert_eq!(idx.chrom_of(62861), Some(2), "end of the genome");
        assert_eq!(idx.chrom_of(0), None, "bins are 1-based");
        assert_eq!(idx.chrom_of(62862), None, "past the end");
    }

    #[test]
    fn detects_boundary_crossing_pairs() {
        let idx = pombe(200);
        assert!(idx.same_chrom(27890, 27896), "both on I");
        assert!(!idx.same_chrom(27896, 27897), "I to II");
        assert!(!idx.same_chrom(100, 62861), "I to III");
        assert!(!idx.same_chrom(1, 62862), "one is out of range");
    }

    #[test]
    fn reports_boundaries_between_chromosomes() {
        assert_eq!(pombe(200).boundaries(), vec![27897, 50597]);
    }

    #[test]
    fn estimates_pombe_bin_size() {
        // Real data: 12,571,820 bp over a max bin of 62,843 = 200.05 bp/bin.
        assert_eq!(estimate_bin_size(12_571_820, 62_843), 200);
    }

    #[test]
    fn estimates_round_bin_sizes() {
        assert_eq!(estimate_bin_size(12_571_820, 12_571), 1000);
        assert_eq!(estimate_bin_size(12_571_820, 2_514), 5000);
    }

    #[test]
    fn rejects_malformed_tables() {
        assert!(matches!(
            ChromIndex::parse("I 5579133\n", 200, "t"),
            Err(Error::Parse { .. })
        ));
        assert!(matches!(
            ChromIndex::parse("I\tabc\n", 200, "t"),
            Err(Error::Parse { .. })
        ));
        assert!(matches!(
            ChromIndex::parse("", 200, "t"),
            Err(Error::Empty(_))
        ));
    }

    #[test]
    fn skips_blank_and_comment_lines() {
        let idx = ChromIndex::parse("# comment\n\nI\t5579133\n", 200, "t").unwrap();
        assert_eq!(idx.chroms().len(), 1);
    }
}
