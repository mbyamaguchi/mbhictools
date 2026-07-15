//! `draw_hic`: draws a Hi-C contact map to PNG.
//!
//! The aggregation and rendering live in the `mbhictools_rs` library; this only
//! wires them together.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;

use mbhictools_rs::chrom::{self, ChromIndex};
use mbhictools_rs::font;
use mbhictools_rs::grid::{self, GridSpec};
use mbhictools_rs::render::{self, Labels, Palette, Scale};

/// Draw a Hi-C contact map from sparse contacts (bin1<TAB>bin2<TAB>score).
///
/// Assumes the upper triangle only, integer scores >= 1, and distances within a
/// limit. Because distance is capped, the map is drawn not as a square but as a
/// triangle rotated 45 degrees: genomic position across, interaction distance up.
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Contact data (TSV: bin1 bin2 score)
    input: PathBuf,

    /// Output PNG
    #[arg(short, long, default_value = "hicmap.png")]
    output: PathBuf,

    /// Pixels across (the resolution). Below 1 bin per pixel is rejected
    #[arg(long, default_value_t = 4000)]
    nx: usize,

    /// Start of the genomic range shown (bins). Defaults to all data
    #[arg(long)]
    x_start: Option<f64>,

    /// End of the genomic range shown (bins, exclusive). Defaults to all data
    #[arg(long)]
    x_end: Option<f64>,

    /// Largest bin distance shown (bins, inclusive). Defaults to the data's largest
    #[arg(long)]
    max_distance: Option<f64>,

    /// Chromosome lengths (chr<TAB>bp). Labels axes in bp, drops contacts crossing
    /// a chromosome boundary, and marks the boundaries
    #[arg(long)]
    lengths: Option<PathBuf>,

    /// bp per bin. Estimated from --lengths and the max bin when omitted
    #[arg(long)]
    bin_size: Option<u32>,

    /// Keep contacts that cross a chromosome boundary (only affects --lengths)
    #[arg(long)]
    keep_inter: bool,

    /// Colour by the raw summed score instead of log10
    #[arg(long)]
    no_log: bool,

    /// Quantile at which to clip the top (0..1)
    #[arg(long, default_value_t = 0.99)]
    trim_quantile: f64,

    /// Do not clip the top at all
    #[arg(long)]
    no_trim: bool,

    /// Figure title. Defaults to the input file name
    #[arg(long)]
    title: Option<String>,

    /// Font family for labels. Auto-selected when omitted
    #[arg(long)]
    font: Option<String>,

    /// Write the aggregated grid to TSV (px py x y raw value)
    #[arg(long)]
    dump_grid: Option<PathBuf>,

    /// Cap on total pixels. The grid costs 8 bytes each
    #[arg(long, default_value_t = 40_000_000)]
    max_pixels: usize,

    /// Render even below 1 bin per pixel, moire and all
    #[arg(long)]
    allow_aliasing: bool,
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
    if !(0.0..=1.0).contains(&args.trim_quantile) {
        return Err(format!(
            "--trim-quantile must be within 0..1 (got {})",
            args.trim_quantile
        )
        .into());
    }
    if args.nx == 0 {
        return Err("--nx must be at least 1".into());
    }

    // --- Pass 1: measure the data, to resolve the range and estimate bin_size ---
    eprintln!("[1/3] scanning extent: {}", args.input.display());
    let extent = grid::scan_extent(&args.input)?;
    eprintln!(
        "      bins {}..{}, max distance {} bins",
        extent.bin_min, extent.bin_max, extent.dist_max
    );

    // --- bin_size and the chromosome table ---
    let bin_size = resolve_bin_size(&args, &extent)?;
    let chroms = match (&args.lengths, bin_size) {
        (Some(path), Some(bs)) => {
            let idx = ChromIndex::load(path, bs)?;
            eprintln!(
                "      {} chromosomes over {} bins ({} bp/bin)",
                idx.chroms().len(),
                idx.total_bins(),
                bs
            );
            if extent.bin_max > idx.total_bins() {
                eprintln!(
                    "warning: max bin {} exceeds the {} bins of the chromosome table; \
                     --bin-size may be wrong",
                    extent.bin_max,
                    idx.total_bins()
                );
            }
            Some(idx)
        }
        _ => None,
    };

    // --- Resolve the view ---
    // xr = (bin1 + bin2) / 2 peaks at bin_max, so the half-open end is one past it.
    let x0 = args.x_start.unwrap_or(extent.bin_min as f64);
    let x1 = args.x_end.unwrap_or(extent.bin_max as f64 + 1.0);
    if x1 <= x0 {
        return Err(format!("empty range: [{x0}, {x1})").into());
    }
    let max_distance = args.max_distance.unwrap_or(extent.dist_max as f64).max(1.0);

    let spec = GridSpec::new(x0, x1, max_distance, args.nx);
    report_geometry(&spec, bin_size);

    // --- Check the pixel count ---
    if let Some(a) = spec.aliasing() {
        if args.allow_aliasing {
            eprintln!("warning: {a}");
        } else {
            return Err(format!("{a}\n(pass --allow-aliasing to render anyway)").into());
        }
    }
    if spec.pixels() > args.max_pixels {
        return Err(format!(
            "{} pixels exceeds the cap of {} (the grid would need {:.1} GB).\n\
             Lower --nx or raise --max-pixels",
            spec.pixels(),
            args.max_pixels,
            spec.pixels() as f64 * 8.0 / 1e9
        )
        .into());
    }

    // --- Pass 2: aggregate ---
    eprintln!("[2/3] aggregating ({} x {} pixels)", spec.nx, spec.ny);
    let filter = if args.keep_inter {
        None
    } else {
        chroms.as_ref()
    };
    let g = grid::build(&args.input, spec, filter)?;
    eprintln!(
        "      {} rows counted, {} out of range{}",
        g.counted,
        g.out_of_range,
        if g.inter_chrom > 0 {
            format!(", {} dropped crossing chromosomes", g.inter_chrom)
        } else {
            String::new()
        }
    );
    if g.counted == 0 {
        return Err(
            "no contacts fall in the view; check --x-start / --x-end / --max-distance".into(),
        );
    }

    // --- Transform and draw ---
    let scale = Scale {
        log: !args.no_log,
        trim_quantile: if args.no_trim {
            None
        } else {
            Some(args.trim_quantile)
        },
    };
    let values = render::transform(&g, &scale);
    eprintln!(
        "      {} of {} pixels filled ({:.1}%), colour range {:.3}..{:.3}",
        values.filled,
        spec.pixels(),
        values.filled as f64 / spec.pixels() as f64 * 100.0,
        values.vmin,
        values.vmax
    );

    if let Some(path) = &args.dump_grid {
        dump_grid(path, &g, &values)?;
        eprintln!("      wrote the grid to {}", path.display());
    }

    eprintln!("[3/3] drawing: {}", args.output.display());
    let font = font::pick(args.font.as_deref());
    eprintln!("      font: {font}");
    let labels = Labels {
        title: args
            .title
            .clone()
            .unwrap_or_else(|| default_title(&args.input)),
        bin_size,
        legend: if scale.log {
            "log10(score)".into()
        } else {
            "score".into()
        },
        font,
    };
    render::render(
        &args.output,
        &g,
        &values,
        &Palette::default(),
        &labels,
        chroms.as_ref(),
    )?;
    eprintln!("done: {}", args.output.display());
    Ok(())
}

/// Resolve bp/bin: explicit, else estimated from the lengths, else unknown.
fn resolve_bin_size(
    args: &Args,
    extent: &grid::Extent,
) -> Result<Option<u32>, Box<dyn std::error::Error>> {
    if let Some(bs) = args.bin_size {
        if bs == 0 {
            return Err("--bin-size must be at least 1".into());
        }
        return Ok(Some(bs));
    }
    let Some(path) = &args.lengths else {
        return Ok(None); // without lengths there is no way to convert to bp
    };
    let total_bp = ChromIndex::total_bp(path)?;
    let bs = chrom::estimate_bin_size(total_bp, extent.bin_max);
    eprintln!(
        "      bin_size estimated at {bs} bp/bin ({total_bp} bp over max bin {})",
        extent.bin_max
    );
    Ok(Some(bs))
}

fn report_geometry(spec: &GridSpec, bin_size: Option<u32>) {
    let bp = |v: f64| match bin_size {
        Some(bs) => format!(" ({:.2} Mb)", v * bs as f64 / 1e6),
        None => String::new(),
    };
    eprintln!(
        "      showing bins {:.0}..{:.0}{}, max distance {:.0} bins{}",
        spec.x0,
        spec.x1,
        bp(spec.x1 - spec.x0),
        spec.max_distance,
        bp(spec.max_distance)
    );
    let per_px = match bin_size {
        Some(bs) => format!(" = {:.0} bp/pixel", spec.dx * bs as f64),
        None => String::new(),
    };
    eprintln!("      resolution {:.3} bins/pixel{}", spec.dx, per_px);
}

fn default_title(input: &Path) -> String {
    let name = input.file_stem().map_or_else(
        || input.display().to_string(),
        |s| s.to_string_lossy().into_owned(),
    );
    format!("Hi-C contact map: {name}")
}

/// Write the aggregate as TSV, skipping empty pixels. x and y are the pixel centre
/// in bins; y is the rotated yr, i.e. half the distance.
fn dump_grid(path: &Path, g: &grid::Grid, values: &render::Values) -> std::io::Result<()> {
    let file = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(file);
    writeln!(w, "px\tpy\tx\ty\traw\tvalue")?;
    for py in 0..g.spec.ny {
        for px in 0..g.spec.nx {
            let Some(v) = values.v[py * g.spec.nx + px] else {
                continue;
            };
            let (x, y) = g.spec.pixel_center(px, py);
            writeln!(w, "{px}\t{py}\t{x}\t{y}\t{}\t{v:.6}", g.get(px, py))?;
        }
    }
    w.flush()
}
