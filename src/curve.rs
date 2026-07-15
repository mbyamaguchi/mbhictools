//! The contact-frequency versus genomic-distance curve, P(s).
//!
//! For intra-chromosomal contacts at a separation of `s` bins,
//!
//! ```text
//!   P(s) = contacts(s) / pairs(s)
//! ```
//!
//! where `contacts(s)` sums the scores of every observed pair at that separation and
//! `pairs(s)` counts the pairs that could have been observed. The denominator is the
//! whole difficulty: the input is sparse, so a pair absent from the file contributes
//! nothing to the numerator and must still be counted below the line.
//!
//! # Zeros are not stored, they are counted
//!
//! An absent pair is not missing data, and it does not need to be materialised as a
//! zero row. `pairs(s)` is a property of the genome, not of the file, so it is
//! computed from the chromosome table: a chromosome of `N` bins holds `N - s` pairs
//! at separation `s`.
//!
//! # Except where a bin cannot be measured at all
//!
//! That closed form assumes every bin is measurable, and roughly 3.6% of them are
//! not: centromeres, rDNA and other repeats never receive a read, at any depth. A
//! pair touching one of those is undefined, not zero. Counting it in the denominator
//! while it can never reach the numerator biases P(s) downwards.
//!
//! So [`mask`] marks bins that carry no coverage, and [`pairs_per_distance`] counts
//! only pairs of measurable bins. That costs an O(bins * smax) sweep instead of a
//! closed form, which is cheap next to reading the file.

use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

use rayon::prelude::*;

use crate::chrom::ChromIndex;
use crate::contact;

/// Per-distance totals, indexed by separation in bins.
#[derive(Debug, Clone)]
pub struct Totals {
    /// Summed score at each separation.
    pub contacts: Vec<u64>,
    /// Pairs that could have been observed at each separation.
    pub pairs: Vec<u64>,
    /// Largest separation covered, inclusive.
    pub smax: u32,
    /// Rows dropped for crossing a chromosome boundary.
    pub inter_chrom: u64,
    /// Rows dropped for touching an unmeasurable bin.
    pub masked: u64,
}

/// One plotted point: a separation, or a range of them once log binned.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    /// Separations covered, inclusive.
    pub s_lo: u32,
    pub s_hi: u32,
    /// Representative separation in bins: the pair-weighted geometric mean, which is
    /// the centre of mass of the denominator on a log axis.
    pub s: f64,
    pub contacts: u64,
    pub pairs: u64,
    /// `contacts / pairs`. Zero when the range was observed but held no contacts.
    pub prob: f64,
}

/// Total score touching each bin, indexed by global bin (element 0 is unused).
pub fn coverage(path: &Path, total_bins: u32) -> Result<Vec<u64>, contact::Error> {
    let cov: Vec<AtomicU64> = (0..=total_bins as usize)
        .map(|_| AtomicU64::new(0))
        .collect();
    let add = |b: u32, score: u32| {
        if let Some(c) = cov.get(b as usize) {
            c.fetch_add(score as u64, Ordering::Relaxed);
        }
    };
    contact::visit(path, |b1, b2, score| {
        add(b1, score);
        // The diagonal is one pair, so do not count it against the bin twice.
        if b2 != b1 {
            add(b2, score);
        }
    })?;
    Ok(cov.into_iter().map(AtomicU64::into_inner).collect())
}

/// Mark bins carrying at least `min_coverage` score as measurable.
///
/// Bin 0 does not exist (bins are 1-based) and is always false.
pub fn mask(cov: &[u64], min_coverage: u64) -> Vec<bool> {
    let threshold = min_coverage.max(1);
    cov.iter().map(|&c| c >= threshold).collect()
}

/// Sum scores per separation, keeping only intra-chromosomal pairs of valid bins.
pub fn contacts_per_distance(
    path: &Path,
    idx: &ChromIndex,
    valid: Option<&[bool]>,
    smax: u32,
) -> Result<(Vec<u64>, u64, u64), contact::Error> {
    let num: Vec<AtomicU64> = (0..=smax as usize).map(|_| AtomicU64::new(0)).collect();
    let inter = AtomicU64::new(0);
    let masked = AtomicU64::new(0);

    let measurable = |b: u32| valid.is_none_or(|v| v.get(b as usize).copied().unwrap_or(false));

    contact::visit(path, |b1, b2, score| {
        if !idx.same_chrom(b1, b2) {
            inter.fetch_add(1, Ordering::Relaxed);
            return;
        }
        if !measurable(b1) || !measurable(b2) {
            masked.fetch_add(1, Ordering::Relaxed);
            return;
        }
        let s = b2.abs_diff(b1);
        if s <= smax {
            num[s as usize].fetch_add(score as u64, Ordering::Relaxed);
        }
    })?;

    Ok((
        num.into_iter().map(AtomicU64::into_inner).collect(),
        inter.load(Ordering::Relaxed),
        masked.load(Ordering::Relaxed),
    ))
}

/// Count the pairs that could have been observed at each separation.
///
/// Without a mask this is the closed form `N - s` per chromosome. With one, only
/// pairs whose both bins are measurable count, which needs an actual sweep.
pub fn pairs_per_distance(idx: &ChromIndex, valid: Option<&[bool]>, smax: u32) -> Vec<u64> {
    (0..=smax)
        .into_par_iter()
        .map(|s| {
            let mut n = 0u64;
            for c in idx.chroms() {
                if c.bins <= s {
                    continue; // too short to hold this separation
                }
                match valid {
                    None => n += (c.bins - s) as u64,
                    Some(v) => {
                        for i in c.start_bin..=(c.end_bin - s) {
                            let ok = v.get(i as usize).copied().unwrap_or(false)
                                && v.get((i + s) as usize).copied().unwrap_or(false);
                            n += ok as u64;
                        }
                    }
                }
            }
            n
        })
        .collect()
}

/// Read the file once and build the per-distance totals.
pub fn totals(
    path: &Path,
    idx: &ChromIndex,
    valid: Option<&[bool]>,
    smax: u32,
) -> Result<Totals, contact::Error> {
    let (contacts, inter_chrom, masked) = contacts_per_distance(path, idx, valid, smax)?;
    let pairs = pairs_per_distance(idx, valid, smax);
    Ok(Totals {
        contacts,
        pairs,
        smax,
        inter_chrom,
        masked,
    })
}

/// Reduce the totals to plottable points.
///
/// `bins_per_decade` log bins the separations; `None` keeps every separation.
///
/// Log binning is what makes P(s) readable. Per-separation points crowd the right of
/// a log axis and grow noisy exactly where pairs are scarcest, so the convention is
/// to pool them into geometrically spaced ranges.
///
/// Within a range the numerator and denominator are summed separately, and divided
/// once. Averaging the per-separation P(s) instead would weight a separation holding
/// a handful of pairs the same as one holding thousands.
pub fn points(t: &Totals, min_s: u32, bins_per_decade: Option<f64>) -> Vec<Point> {
    let min_s = min_s.max(1); // s = 0 is the diagonal and has no place on a log axis
    if min_s > t.smax {
        return Vec::new();
    }
    let edges = match bins_per_decade {
        Some(per_decade) => log_edges(min_s, t.smax, per_decade),
        None => (min_s..=t.smax + 1).collect(),
    };

    let mut out = Vec::with_capacity(edges.len());
    for w in edges.windows(2) {
        let (lo, hi) = (w[0], w[1]); // [lo, hi)
        let mut contacts = 0u64;
        let mut pairs = 0u64;
        let mut weight = 0f64;
        let mut ln_sum = 0f64;

        for s in lo..hi {
            contacts += t.contacts[s as usize];
            let p = t.pairs[s as usize];
            pairs += p;
            if p > 0 {
                weight += p as f64;
                ln_sum += p as f64 * (s as f64).ln();
            }
        }
        if pairs == 0 {
            continue; // nothing could have been observed here
        }
        out.push(Point {
            s_lo: lo,
            s_hi: hi - 1,
            s: (ln_sum / weight).exp(),
            contacts,
            pairs,
            prob: contacts as f64 / pairs as f64,
        });
    }
    out
}

/// Geometrically spaced, integral, strictly increasing edges over [min_s, max_s + 1].
///
/// Rounding to integers collapses the closest edges near the left, where a decade
/// spans only a few separations; duplicates are dropped rather than left empty.
fn log_edges(min_s: u32, max_s: u32, per_decade: f64) -> Vec<u32> {
    assert!(min_s >= 1, "min_s must be at least 1");
    assert!(max_s >= min_s, "max_s must be at least min_s");
    assert!(per_decade > 0.0, "bins per decade must be positive");

    let (lo, hi) = (min_s as f64, (max_s + 1) as f64);
    let n = ((hi / lo).log10() * per_decade).ceil().max(1.0) as usize;

    let mut edges = Vec::with_capacity(n + 1);
    for i in 0..=n {
        let v = ((lo * 10f64.powf(i as f64 / per_decade)).round() as u32).min(max_s + 1);
        if edges.last().is_none_or(|&last| v > last) {
            edges.push(v);
        }
    }
    if edges.last().is_some_and(|&last| last < max_s + 1) {
        edges.push(max_s + 1);
    }
    edges
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pombe(bin_size: u32) -> ChromIndex {
        ChromIndex::parse("I\t5579133\nII\t4539804\nIII\t2452883\n", bin_size, "test").unwrap()
    }

    /// A tiny genome: one chromosome of 10 bins.
    fn tiny() -> ChromIndex {
        ChromIndex::parse("c\t10\n", 1, "test").unwrap()
    }

    #[test]
    fn counts_pairs_with_the_closed_form() {
        let idx = tiny(); // 10 bins
        let p = pairs_per_distance(&idx, None, 4);
        assert_eq!(p, vec![10, 9, 8, 7, 6], "N - s pairs at separation s");
    }

    #[test]
    fn closed_form_matches_the_real_genome() {
        let idx = pombe(200);
        let p = pairs_per_distance(&idx, None, 5000);
        // 27896 + 22700 + 12265 bins at s = 0, three fewer at s = 1, and so on.
        assert_eq!(p[0], 62861);
        assert_eq!(p[1], 62861 - 3);
        assert_eq!(p[5000], (27896 - 5000) + (22700 - 5000) + (12265 - 5000));
        // Matches the total measured on the real data.
        assert_eq!(p.iter().sum::<u64>(), 276_860_361);
    }

    #[test]
    fn skips_chromosomes_shorter_than_the_separation() {
        let idx = pombe(200); // III holds 12265 bins
        let p = pairs_per_distance(&idx, None, 12_300);
        assert_eq!(
            p[12_300],
            (27896 - 12_300) + (22700 - 12_300),
            "III cannot reach this separation"
        );
    }

    #[test]
    fn mask_excludes_pairs_touching_unmeasurable_bins() {
        let idx = tiny(); // bins 1..=10
        // Bin 5 carries no coverage. valid is indexed by bin, so element 0 is unused.
        let cov = vec![0, 7, 7, 7, 7, 0, 7, 7, 7, 7, 7];
        let valid = mask(&cov, 1);
        assert!(!valid[0], "bin 0 does not exist");
        assert!(!valid[5], "bin 5 is unmeasurable");

        let p = pairs_per_distance(&idx, Some(&valid), 2);
        // s = 0: the 9 measurable bins pair with themselves.
        assert_eq!(p[0], 9);
        // s = 1: 9 pairs exist, but (4,5) and (5,6) are lost.
        assert_eq!(p[1], 7);
        // s = 2: 8 pairs exist, but (3,5) and (5,7) are lost.
        assert_eq!(p[2], 6);
    }

    #[test]
    fn masking_lowers_the_denominator() {
        let idx = tiny();
        let all_valid = vec![
            false, true, true, true, true, true, true, true, true, true, true,
        ];
        assert_eq!(
            pairs_per_distance(&idx, Some(&all_valid), 4),
            pairs_per_distance(&idx, None, 4),
            "an all-valid mask must agree with the closed form"
        );
    }

    #[test]
    fn coverage_counts_the_diagonal_once() {
        let path = std::env::temp_dir().join("mbhictools_curve_cov.txt");
        std::fs::write(&path, "bin1\tbin2\tscore\n1\t1\t5\n1\t2\t3\n").unwrap();
        let cov = coverage(&path, 3).unwrap();
        assert_eq!(cov[1], 8, "5 from the diagonal, 3 from the pair with bin 2");
        assert_eq!(cov[2], 3);
        assert_eq!(cov[3], 0);
    }

    #[test]
    fn computes_p_of_s_end_to_end() {
        let idx = tiny(); // one chromosome, bins 1..=10
        let path = std::env::temp_dir().join("mbhictools_curve_ps.txt");
        // Two pairs at s = 1 with scores 4 and 6; nine pairs exist at s = 1.
        std::fs::write(&path, "bin1\tbin2\tscore\n1\t2\t4\n5\t6\t6\n").unwrap();

        let t = totals(&path, &idx, None, 3).unwrap();
        assert_eq!(t.contacts[1], 10);
        assert_eq!(t.pairs[1], 9);
        assert_eq!(t.contacts[2], 0, "nothing observed at s = 2");
        assert_eq!(t.pairs[2], 8, "but eight pairs could have been");

        let pts = points(&t, 1, None);
        assert_eq!(pts[0].s_lo, 1);
        assert!((pts[0].prob - 10.0 / 9.0).abs() < 1e-12);
        assert_eq!(
            pts[1].prob, 0.0,
            "an observed range with no contacts is P = 0"
        );
    }

    #[test]
    fn drops_inter_chromosomal_rows() {
        let idx = pombe(200);
        let path = std::env::temp_dir().join("mbhictools_curve_inter.txt");
        // 27896 ends chromosome I and 27897 starts II, so this pair spans them.
        std::fs::write(&path, "bin1\tbin2\tscore\n27896\t27897\t9\n100\t101\t4\n").unwrap();

        let t = totals(&path, &idx, None, 5).unwrap();
        assert_eq!(t.inter_chrom, 1);
        assert_eq!(t.contacts[1], 4, "only the intra-chromosomal pair counts");
    }

    #[test]
    fn log_edges_increase_and_span_the_range() {
        let e = log_edges(1, 5000, 10.0);
        assert_eq!(e[0], 1);
        assert_eq!(*e.last().unwrap(), 5001, "the end is exclusive");
        assert!(e.windows(2).all(|w| w[1] > w[0]), "strictly increasing");
    }

    #[test]
    fn log_edges_collapse_duplicates_near_one() {
        // At 10 per decade the first few edges all round to 1, 1, 2, 2, ...
        let e = log_edges(1, 100, 10.0);
        assert!(e.windows(2).all(|w| w[1] > w[0]), "no empty ranges");
        assert_eq!(e[0], 1);
    }

    /// The point of log binning: pool the numerator and denominator, divide once.
    #[test]
    fn log_binning_pools_totals_rather_than_averaging_ratios() {
        let t = Totals {
            //          s = 0  1     2
            contacts: vec![0, 1, 100],
            pairs: vec![0, 1, 1000],
            smax: 2,
            inter_chrom: 0,
            masked: 0,
        };
        // One bin covering s = 1 and 2.
        let pts = points(&t, 1, Some(0.5));
        assert_eq!(pts.len(), 1);
        let p = pts[0];
        assert_eq!((p.contacts, p.pairs), (101, 1001));
        assert!(
            (p.prob - 101.0 / 1001.0).abs() < 1e-12,
            "must pool totals, not average 1/1 and 100/1000"
        );
        // Averaging the ratios would have given (1.0 + 0.1) / 2 = 0.55.
        assert!(p.prob < 0.2);
    }

    #[test]
    fn representative_s_leans_towards_the_denominator() {
        let t = Totals {
            contacts: vec![0, 0, 0, 0],
            pairs: vec![0, 1, 1, 1000],
            smax: 3,
            inter_chrom: 0,
            masked: 0,
        };
        let pts = points(&t, 1, Some(0.1)); // one wide bin over s = 1..=3
        assert_eq!(pts.len(), 1);
        assert!(
            pts[0].s > 2.9,
            "s = 3 holds nearly all the pairs, so it should dominate (got {})",
            pts[0].s
        );
    }

    #[test]
    fn skips_ranges_where_nothing_could_be_observed() {
        let t = Totals {
            contacts: vec![0, 0],
            pairs: vec![0, 0],
            smax: 1,
            inter_chrom: 0,
            masked: 0,
        };
        assert!(points(&t, 1, None).is_empty());
    }

    #[test]
    fn handles_a_min_s_past_the_end() {
        let t = Totals {
            contacts: vec![0, 5],
            pairs: vec![0, 5],
            smax: 1,
            inter_chrom: 0,
            masked: 0,
        };
        assert!(points(&t, 99, None).is_empty());
    }
}
