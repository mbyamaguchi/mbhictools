//! Aggregating contacts into a display grid in rotated coordinates.
//!
//! # Why rotate
//!
//! The input holds only the upper triangle, and only within a distance limit
//! (1 Mbp by default). A square contact map would therefore be mostly empty by
//! construction: only a thin band along the diagonal carries data.
//!
//! Rotating each contact by 45 degrees,
//!
//! ```text
//!   xr = (bin1 + bin2) / 2   genomic position (midpoint of the two bins)
//!   yr = (bin2 - bin1) / 2   half the interaction distance
//! ```
//!
//! turns that band into a triangle whose base is the genome and whose height is half
//! the distance limit, so no pixels are wasted.
//!
//! # Choosing the pixel count
//!
//! Let `dx` be the bins per pixel. Pixels are square (`dy == dx`), which keeps the
//! figure isotropic in data units so its aspect ratio reports distance honestly.
//!
//! Raising `nx` too far breaks this. Rotated contacts do not land anywhere
//! continuous; they land on a lattice, and that lattice is a checkerboard. `bin1` and
//! `bin2` are integers, so `xr` and `yr` are multiples of 0.5 — but not every
//! combination occurs, because `xr + yr = bin2` is always an integer. Only points
//! with integral `xr + yr` exist.
//!
//! Against that lattice, a `dx`-square pixel contains:
//!
//! ```text
//!   dx = 1     exactly 2 points in every pixel (one all-integer, one all-half)
//!   dx = 0.5   1 point where (i+j) is even, 0 where odd: a full checkerboard
//!   dx >= 1    2*dx^2 on average, and more uniform the larger dx is
//! ```
//!
//! So `dx = 1 bin/pixel` is a hard floor. Below it, empty pixels alternate into
//! moiré: the figure gets coarser even though no data changed, and it happens while
//! "increasing the resolution". [`GridSpec::aliasing`] detects it.

use std::path::Path;
use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};

use crate::chrom::ChromIndex;
use crate::contact;

/// Grid geometry, all in bins (never bp).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GridSpec {
    /// Genomic range shown, [x0, x1), in `xr`.
    pub x0: f64,
    pub x1: f64,
    /// Bins per pixel, the same both ways (square pixels).
    pub dx: f64,
    /// Pixels across.
    pub nx: usize,
    /// Pixels up; follows from `ymax / dx`.
    pub ny: usize,
    /// Triangle height, `max_distance / 2`.
    pub ymax: f64,
    /// Largest `|bin2 - bin1|` shown.
    pub max_distance: f64,
}

/// Diagnosis of pixels too fine for the lattice, i.e. `nx` set too high.
#[derive(Debug, Clone, PartialEq)]
pub struct Aliasing {
    /// Current bins per pixel (< 1).
    pub dx: f64,
    /// Largest `nx` that avoids moiré.
    pub max_nx: usize,
}

impl std::fmt::Display for Aliasing {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "a pixel spans only {:.3} bins (< 1 bin/pixel).\n\
             Rotated contacts exist only where `xr + yr` is an integer, so at this\n\
             resolution empty pixels alternate into checkerboard moire.\n\
             Use --nx {} or less for this range",
            self.dx, self.max_nx
        )
    }
}

/// Extent of the data, used to resolve an omitted range.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Extent {
    pub bin_min: u32,
    pub bin_max: u32,
    pub dist_max: u32,
}

/// An aggregated grid. `cells[py * nx + px]` holds the summed score.
#[derive(Debug, Clone)]
pub struct Grid {
    pub spec: GridSpec,
    cells: Vec<u64>,
    /// Rows added to the grid.
    pub counted: u64,
    /// Rows outside the view.
    pub out_of_range: u64,
    /// Rows dropped for crossing a chromosome boundary (0 when not filtering).
    pub inter_chrom: u64,
}

impl GridSpec {
    /// Derive the geometry from a range and a pixel width.
    ///
    /// `ny` follows from `ymax / dx` to keep pixels square. Height is deliberately
    /// not settable on its own: that would let the aspect ratio misreport distance.
    pub fn new(x0: f64, x1: f64, max_distance: f64, nx: usize) -> Self {
        assert!(nx >= 1, "nx must be at least 1");
        assert!(x1 > x0, "empty range: [{x0}, {x1})");
        assert!(max_distance > 0.0, "max_distance must be positive");

        let dx = (x1 - x0) / nx as f64;
        let ymax = max_distance / 2.0;
        let ny = (ymax / dx).ceil().max(1.0) as usize;
        GridSpec {
            x0,
            x1,
            dx,
            nx,
            ny,
            ymax,
            max_distance,
        }
    }

    /// Report moire if the pixels are finer than the lattice supports.
    pub fn aliasing(&self) -> Option<Aliasing> {
        if self.dx >= 1.0 {
            return None;
        }
        // The largest nx with dx = (x1 - x0) / nx >= 1.
        let max_nx = (self.x1 - self.x0).floor().max(1.0) as usize;
        Some(Aliasing {
            dx: self.dx,
            max_nx,
        })
    }

    pub fn pixels(&self) -> usize {
        self.nx * self.ny
    }

    /// The pixel a contact lands in, or `None` if outside the view.
    ///
    /// Position is half-open [x0, x1); distance is closed [0, max_distance]. Closing
    /// distance mirrors the input spec ("within 1 Mbp"): half-open would discard
    /// every contact sitting exactly on the limit (`bin2 - bin1 == 5000` in the data).
    #[inline]
    fn pixel_of(&self, bin1: u32, bin2: u32) -> Option<(usize, usize)> {
        let (b1, b2) = (bin1 as f64, bin2 as f64);
        let xr = (b1 + b2) * 0.5;
        let yr = (b2 - b1) * 0.5;
        if xr < self.x0 || xr >= self.x1 || yr < 0.0 || yr > self.ymax {
            return None;
        }
        let px = ((xr - self.x0) / self.dx) as usize;
        // Rounding, or yr sitting exactly on ymax, can overshoot the top row by one.
        let py = ((yr / self.dx) as usize).min(self.ny - 1);
        if px >= self.nx {
            return None;
        }
        Some((px, py))
    }

    /// Genomic coordinate of a pixel centre, in bins.
    pub fn pixel_center(&self, px: usize, py: usize) -> (f64, f64) {
        (
            self.x0 + (px as f64 + 0.5) * self.dx,
            (py as f64 + 0.5) * self.dx,
        )
    }
}

impl Grid {
    /// Build a grid from pre-aggregated cells. Scan counters are left at zero.
    pub fn from_parts(spec: GridSpec, cells: Vec<u64>) -> Self {
        assert_eq!(
            cells.len(),
            spec.pixels(),
            "cells do not match the pixel count"
        );
        Grid {
            spec,
            cells,
            counted: 0,
            out_of_range: 0,
            inter_chrom: 0,
        }
    }

    /// `cells[py * nx + px]`.
    pub fn get(&self, px: usize, py: usize) -> u64 {
        self.cells[py * self.spec.nx + px]
    }

    pub fn cells(&self) -> &[u64] {
        &self.cells
    }

    /// Values of the non-empty pixels, for computing quantiles.
    pub fn nonzero(&self) -> Vec<u64> {
        self.cells.iter().copied().filter(|&v| v > 0).collect()
    }
}

/// Scan the file once to measure the data.
///
/// Used to resolve an omitted range before allocating the grid. Scores are not read,
/// so this is lighter than the aggregating pass.
pub fn scan_extent(path: &Path) -> Result<Extent, contact::Error> {
    let bin_min = AtomicU32::new(u32::MAX);
    let bin_max = AtomicU32::new(0);
    let dist_max = AtomicU32::new(0);

    let stats = contact::visit(path, |b1, b2, _| {
        bin_min.fetch_min(b1, Ordering::Relaxed);
        bin_max.fetch_max(b2, Ordering::Relaxed);
        dist_max.fetch_max(b2.saturating_sub(b1), Ordering::Relaxed);
    })?;

    if stats.rows == 0 {
        return Err(contact::Error::TooManyMalformed {
            path: path.display().to_string(),
            malformed: stats.malformed,
            rows: 0,
        });
    }
    Ok(Extent {
        bin_min: bin_min.load(Ordering::Relaxed),
        bin_max: bin_max.load(Ordering::Relaxed),
        dist_max: dist_max.load(Ordering::Relaxed),
    })
}

/// Scan the file once and aggregate it into `spec`.
///
/// Memory is the grid alone (nx * ny * 8 bytes), independent of the row count.
/// Integer scores add exactly, and atomically, into one shared grid.
///
/// `chrom_filter` drops contacts crossing a chromosome boundary.
pub fn build(
    path: &Path,
    spec: GridSpec,
    chrom_filter: Option<&ChromIndex>,
) -> Result<Grid, contact::Error> {
    let cells: Vec<AtomicU64> = (0..spec.pixels()).map(|_| AtomicU64::new(0)).collect();
    let counted = AtomicU64::new(0);
    let out_of_range = AtomicU64::new(0);
    let inter_chrom = AtomicU64::new(0);

    contact::visit(path, |b1, b2, score| {
        if let Some(idx) = chrom_filter
            && !idx.same_chrom(b1, b2)
        {
            inter_chrom.fetch_add(1, Ordering::Relaxed);
            return;
        }
        match spec.pixel_of(b1, b2) {
            Some((px, py)) => {
                cells[py * spec.nx + px].fetch_add(score as u64, Ordering::Relaxed);
                counted.fetch_add(1, Ordering::Relaxed);
            }
            None => {
                out_of_range.fetch_add(1, Ordering::Relaxed);
            }
        }
    })?;

    Ok(Grid {
        spec,
        cells: cells.into_iter().map(AtomicU64::into_inner).collect(),
        counted: counted.load(Ordering::Relaxed),
        out_of_range: out_of_range.load(Ordering::Relaxed),
        inter_chrom: inter_chrom.load(Ordering::Relaxed),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn derives_square_pixels() {
        // 1000 bins over 100 pixels = 10 bins/pixel; a 200 bin limit is 100 bins high.
        let s = GridSpec::new(0.0, 1000.0, 200.0, 100);
        assert_eq!(s.dx, 10.0);
        assert_eq!(s.ymax, 100.0);
        assert_eq!(s.ny, 10);
        assert_eq!(s.pixels(), 1000);
    }

    #[test]
    fn rounds_ny_up_to_cover_the_top() {
        // ymax / dx = 105 / 10 = 10.5, so 11 pixels rather than clipping the top.
        let s = GridSpec::new(0.0, 1000.0, 210.0, 100);
        assert_eq!(s.ny, 11);
    }

    #[test]
    fn maps_contacts_to_rotated_pixels() {
        let s = GridSpec::new(0.0, 1000.0, 200.0, 100);
        // bin1 = bin2 = 5: xr = 5, yr = 0, on the diagonal.
        assert_eq!(s.pixel_of(5, 5), Some((0, 0)));
        // bin1 = 100, bin2 = 120: xr = 110, yr = 10.
        assert_eq!(s.pixel_of(100, 120), Some((11, 1)));
        // xr at or past x1 is outside.
        assert_eq!(s.pixel_of(1000, 1000), None);
        // Distance exactly on the limit (yr == ymax) belongs to the top row.
        assert_eq!(s.pixel_of(0, 200), Some((10, 9)));
        assert_eq!(s.pixel_of(0, 202), None, "past the limit");
    }

    #[test]
    fn pixel_center_matches_mapping() {
        let s = GridSpec::new(0.0, 1000.0, 200.0, 100);
        assert_eq!(s.pixel_center(0, 0), (5.0, 5.0));
        assert_eq!(s.pixel_center(11, 1), (115.0, 15.0));
    }

    /// The module doc's claim: at dx = 1 every pixel holds exactly 2 lattice points.
    #[test]
    fn unit_pixels_receive_exactly_two_lattice_points_each() {
        let s = GridSpec::new(0.0, 40.0, 20.0, 40);
        assert_eq!(s.dx, 1.0);
        assert_eq!(s.aliasing(), None, "dx = 1 is the floor, not below it");

        let mut count = vec![0usize; s.pixels()];
        for b1 in 0..=60u32 {
            for b2 in b1..=60u32 {
                if let Some((px, py)) = s.pixel_of(b1, b2) {
                    count[py * s.nx + px] += 1;
                }
            }
        }
        // Inside the triangle, away from the clipped edges, always exactly 2.
        for py in 0..5 {
            for px in 10..30 {
                assert_eq!(count[py * s.nx + px], 2, "pixel ({px}, {py})");
            }
        }
    }

    /// Below dx = 1, pixels start missing the lattice entirely: checkerboard moire.
    #[test]
    fn subunit_pixels_leave_empty_gaps() {
        let s = GridSpec::new(0.0, 40.0, 20.0, 80);
        assert_eq!(s.dx, 0.5);

        let mut count = vec![0usize; s.pixels()];
        for b1 in 0..=60u32 {
            for b2 in b1..=60u32 {
                if let Some((px, py)) = s.pixel_of(b1, b2) {
                    count[py * s.nx + px] += 1;
                }
            }
        }
        let empty = count.iter().filter(|&&c| c == 0).count();
        assert!(
            empty > s.pixels() / 3,
            "dx = 0.5 should empty about half the pixels (got {empty} of {})",
            s.pixels()
        );
    }

    #[test]
    fn reports_aliasing_with_a_usable_nx_limit() {
        let s = GridSpec::new(0.0, 1000.0, 200.0, 4000);
        let a = s.aliasing().expect("dx = 0.25 aliases");
        assert_eq!(a.dx, 0.25);
        assert_eq!(a.max_nx, 1000, "1000 bins allow at most nx = 1000");
    }

    #[test]
    fn no_aliasing_at_coarse_resolution() {
        // A 62861 bin genome over 4000 pixels is about 15.7 bins/pixel.
        let s = GridSpec::new(0.0, 62861.0, 5000.0, 4000);
        assert!(s.aliasing().is_none());
        assert_eq!(s.ny, 160);
    }

    #[test]
    fn aggregates_scores_per_pixel() {
        let path = std::env::temp_dir().join("mbhictools_grid.txt");
        // Two rows landing in one pixel, plus one elsewhere.
        std::fs::write(
            &path,
            "bin1\tbin2\tscore\n100\t120\t3\n101\t121\t4\n5\t5\t7\n",
        )
        .unwrap();

        let s = GridSpec::new(0.0, 1000.0, 200.0, 100);
        let g = build(&path, s, None).unwrap();

        assert_eq!(g.get(11, 1), 7, "scores in one pixel are summed");
        assert_eq!(g.get(0, 0), 7);
        assert_eq!(g.counted, 3);
        assert_eq!(g.out_of_range, 0);
    }

    #[test]
    fn excludes_contacts_outside_the_view() {
        let path = std::env::temp_dir().join("mbhictools_grid_oor.txt");
        std::fs::write(&path, "bin1\tbin2\tscore\n100\t120\t3\n5000\t5000\t9\n").unwrap();

        let g = build(&path, GridSpec::new(0.0, 1000.0, 200.0, 100), None).unwrap();
        assert_eq!(g.counted, 1);
        assert_eq!(g.out_of_range, 1);
        assert_eq!(g.nonzero(), vec![3]);
    }

    #[test]
    fn scans_extent_in_one_pass() {
        let path = std::env::temp_dir().join("mbhictools_extent.txt");
        std::fs::write(
            &path,
            "bin1\tbin2\tscore\n27\t79\t1\n28\t5028\t2\n60\t100\t1\n",
        )
        .unwrap();

        let e = scan_extent(&path).unwrap();
        assert_eq!(
            e,
            Extent {
                bin_min: 27,
                bin_max: 5028,
                dist_max: 5000
            }
        );
    }
}
