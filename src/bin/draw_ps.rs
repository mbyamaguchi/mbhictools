//! `draw_ps`: draws the Hi-C distance curve P(s).
//!
//! The computation lives in `mbhictools_rs::curve`; this wires it up and plots it.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Parser, ValueEnum};
use plotters::coord::combinators::IntoLogRange;
use plotters::prelude::*;

use mbhictools_rs::chrom::{self, ChromIndex};
use mbhictools_rs::curve::{self, Point};
use mbhictools_rs::font;
use mbhictools_rs::grid;

/// Draw the contact frequency versus genomic distance curve, P(s).
///
/// For intra-chromosomal pairs at separation s, P(s) is the summed score divided by
/// the number of pairs that could have been observed. Absent pairs are counted in
/// that denominator without being stored: it follows from the chromosome lengths.
/// Bins that never receive coverage are excluded from it instead, since a pair
/// touching one is undefined rather than zero.
///
/// Separations are pooled into geometrically spaced bins, summing numerator and
/// denominator separately, which is the convention for reading P(s) on a log axis.
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Contact data (TSV: bin1 bin2 score)
    input: PathBuf,

    /// Output PNG
    #[arg(short, long, default_value = "distance_curve.png")]
    output: PathBuf,

    /// Chromosome lengths (chr<TAB>bp). Required: the denominator comes from them
    #[arg(long)]
    lengths: PathBuf,

    /// bp per bin. Estimated from --lengths and the max bin when omitted
    #[arg(long)]
    bin_size: Option<u32>,

    /// Log bins per decade of distance
    #[arg(long, default_value_t = 10.0)]
    bins_per_decade: f64,

    /// Report every separation instead of log binning them
    #[arg(long)]
    no_logbin: bool,

    /// Bins below this total score are treated as unmeasurable and excluded
    #[arg(long, default_value_t = 1)]
    min_coverage: u64,

    /// Count every bin, even ones that never receive coverage
    #[arg(long)]
    no_mask: bool,

    /// Smallest separation plotted (bins). s = 0 is the diagonal and never plotted
    #[arg(long, default_value_t = 1)]
    min_s: u32,

    /// Largest separation (bins). Defaults to the data's largest
    #[arg(long)]
    max_s: Option<u32>,

    /// Axis scaling
    #[arg(long, value_enum, default_value_t = Scale::LogLog)]
    scale: Scale,

    /// Write the curve to TSV (s_lo s_hi s_bin s_bp contacts pairs prob)
    #[arg(long)]
    dump: Option<PathBuf>,

    /// Figure title. Defaults to the input file name
    #[arg(long)]
    title: Option<String>,

    /// Font family for labels. Auto-selected when omitted
    #[arg(long)]
    font: Option<String>,

    /// Figure width in pixels
    #[arg(long, default_value_t = 900)]
    width: u32,

    /// Figure height in pixels
    #[arg(long, default_value_t = 700)]
    height: u32,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum Scale {
    /// Both axes logarithmic, the convention for P(s)
    LogLog,
    /// Logarithmic distance, linear P(s)
    SemiLogX,
    /// Both axes linear
    Linear,
}

fn main() -> ExitCode {
    match run(Args::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    if args.bins_per_decade <= 0.0 {
        return Err("--bins-per-decade must be positive".into());
    }

    // --- Pass 1: extent, to estimate bin_size and bound the separations ---
    eprintln!("[1/4] scanning extent: {}", args.input.display());
    let extent = grid::scan_extent(&args.input)?;
    eprintln!(
        "      bins {}..{}, max distance {} bins",
        extent.bin_min, extent.bin_max, extent.dist_max
    );

    let bin_size = match args.bin_size {
        Some(bs) if bs > 0 => bs,
        Some(_) => return Err("--bin-size must be at least 1".into()),
        None => {
            let total_bp = ChromIndex::total_bp(&args.lengths)?;
            let bs = chrom::estimate_bin_size(total_bp, extent.bin_max);
            eprintln!(
                "      bin_size estimated at {bs} bp/bin ({total_bp} bp over max bin {})",
                extent.bin_max
            );
            bs
        }
    };
    let idx = ChromIndex::load(&args.lengths, bin_size)?;
    eprintln!(
        "      {} chromosomes over {} bins ({} bp/bin)",
        idx.chroms().len(),
        idx.total_bins(),
        bin_size
    );

    let longest = idx.chroms().iter().map(|c| c.bins).max().unwrap_or(1) - 1;
    let smax = args
        .max_s
        .unwrap_or(extent.dist_max)
        .min(longest)
        .max(args.min_s.max(1));

    // --- Pass 2: coverage, to find bins that cannot be measured ---
    let valid = if args.no_mask {
        eprintln!("[2/4] skipping the coverage mask (--no-mask)");
        None
    } else {
        eprintln!("[2/4] measuring per-bin coverage");
        let cov = curve::coverage(&args.input, idx.total_bins())?;
        let v = curve::mask(&cov, args.min_coverage);
        // Element 0 is the non-existent bin 0, so exclude it from the tally.
        let usable = v[1..].iter().filter(|&&b| b).count();
        let total = idx.total_bins() as usize;
        eprintln!(
            "      {} of {} bins measurable, {} excluded ({:.2}%)",
            usable,
            total,
            total - usable,
            (total - usable) as f64 / total as f64 * 100.0
        );
        Some(v)
    };

    // --- Pass 3: the curve ---
    eprintln!("[3/4] summing contacts and pairs per distance (s <= {smax})");
    let totals = curve::totals(&args.input, &idx, valid.as_deref(), smax)?;
    eprintln!(
        "      {} rows crossed chromosomes, {} touched an excluded bin",
        totals.inter_chrom, totals.masked
    );

    let per_decade = if args.no_logbin {
        None
    } else {
        Some(args.bins_per_decade)
    };
    let points = curve::points(&totals, args.min_s, per_decade);
    if points.is_empty() {
        return Err("no separation holds any pair; check --min-s / --max-s".into());
    }
    eprintln!(
        "      {} points over s = {}..{} bins",
        points.len(),
        points[0].s_lo,
        points[points.len() - 1].s_hi
    );

    if let Some(path) = &args.dump {
        dump(path, &points, bin_size)?;
        eprintln!("      wrote the curve to {}", path.display());
    }

    // --- Pass 4: plot ---
    eprintln!("[4/4] drawing: {}", args.output.display());
    let title = args
        .title
        .clone()
        .unwrap_or_else(|| default_title(&args.input));
    let family = font::pick(args.font.as_deref());
    eprintln!("      font: {family}");
    plot(&args, &points, bin_size, &title, &family)?;
    eprintln!("done: {}", args.output.display());
    Ok(())
}

/// Build the chart, draw the curve, and label the axes.
///
/// The three scalings produce different coordinate types, so the body is shared by a
/// macro rather than by a generic: the trait bounds buy nothing here.
macro_rules! draw_chart {
    ($root:expr, $args:expr, $xy:expr, $x_range:expr, $y_range:expr, $title:expr, $family:expr) => {{
        let mut chart = ChartBuilder::on($root)
            .margin(18)
            .margin_right(30)
            .caption($title, ($family, 17))
            .x_label_area_size(56)
            .y_label_area_size(78)
            .build_cartesian_2d($x_range, $y_range)?;

        chart
            .configure_mesh()
            .x_desc("Genomic distance s (bp)")
            .y_desc("Contact frequency P(s)")
            .axis_desc_style(($family, 15))
            .label_style(($family, 12))
            .light_line_style(RGBColor(0xE8, 0xE8, 0xE8))
            .draw()?;

        chart.draw_series(LineSeries::new(
            $xy.iter().copied(),
            RGBColor(0xCB, 0x18, 0x1D).stroke_width(2),
        ))?;
    }};
}

fn plot(
    args: &Args,
    points: &[Point],
    bin_size: u32,
    title: &str,
    family: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let log_y = args.scale == Scale::LogLog;

    // A log axis cannot show zero. P(s) = 0 means the range was measurable but held
    // no contacts, which is a real observation, so say how many points it costs.
    let xy: Vec<(f64, f64)> = points
        .iter()
        .filter(|p| !log_y || p.prob > 0.0)
        .map(|p| (p.s * bin_size as f64, p.prob))
        .collect();
    let dropped = points.len() - xy.len();
    if dropped > 0 {
        eprintln!("      {dropped} points with P(s) = 0 omitted (a log axis cannot show zero)");
    }
    if xy.is_empty() {
        return Err("every point has P(s) = 0; try --scale linear".into());
    }

    let (x_lo, x_hi) = bounds(xy.iter().map(|p| p.0));
    let (y_lo, y_hi) = bounds(xy.iter().map(|p| p.1));

    let root = BitMapBackend::new(&args.output, (args.width, args.height)).into_drawing_area();
    root.fill(&WHITE)?;

    match args.scale {
        Scale::LogLog => {
            // Pad by a decade fraction so the extreme points are not on the frame.
            let x = (x_lo * 0.8..x_hi * 1.25).log_scale();
            let y = (y_lo * 0.7..y_hi * 1.4).log_scale();
            draw_chart!(&root, args, xy, x, y, title, family);
        }
        Scale::SemiLogX => {
            let x = (x_lo * 0.8..x_hi * 1.25).log_scale();
            let y = 0.0..y_hi * 1.05;
            draw_chart!(&root, args, xy, x, y, title, family);
        }
        Scale::Linear => {
            let x = 0.0..x_hi * 1.02;
            let y = 0.0..y_hi * 1.05;
            draw_chart!(&root, args, xy, x, y, title, family);
        }
    }

    root.present()?;
    Ok(())
}

/// Smallest and largest of a non-empty series.
fn bounds(vs: impl Iterator<Item = f64>) -> (f64, f64) {
    vs.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
        (lo.min(v), hi.max(v))
    })
}

fn default_title(input: &Path) -> String {
    let name = input.file_stem().map_or_else(
        || input.display().to_string(),
        |s| s.to_string_lossy().into_owned(),
    );
    format!("Hi-C distance curve P(s): {name}")
}

/// Write the curve as TSV. `s_bin` is the representative separation the point plots
/// at; `s_lo`/`s_hi` bound the separations it pools.
fn dump(path: &Path, points: &[Point], bin_size: u32) -> std::io::Result<()> {
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);
    writeln!(w, "s_lo\ts_hi\ts_bin\ts_bp\tcontacts\tpairs\tprob")?;
    for p in points {
        writeln!(
            w,
            "{}\t{}\t{:.4}\t{:.1}\t{}\t{}\t{:.8e}",
            p.s_lo,
            p.s_hi,
            p.s,
            p.s * bin_size as f64,
            p.contacts,
            p.pairs,
            p.prob
        )?;
    }
    w.flush()
}
