//! `draw_hic`: Hi-C contact map を PNG へ描く CLI。
//!
//! 集計と描画の実体は `mbhictools_rs` ライブラリ側にある。ここはその組み立てだけ。

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Parser;

use mbhictools_rs::chrom::{self, ChromIndex};
use mbhictools_rs::font;
use mbhictools_rs::grid::{self, GridSpec};
use mbhictools_rs::render::{self, Labels, Palette, Scale};

/// スパース接触データ (bin1<TAB>bin2<TAB>score) から Hi-C contact map を描く。
///
/// 入力は上三角成分のみ・score >= 1 の整数・bin 間距離が上限以内、を前提とする。
/// 距離が制限されているため、正方形ではなく 45 度回した三角形として描く
/// (横 = ゲノム位置、縦 = 相互作用距離)。
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// 接触データ (TSV: bin1 bin2 score)
    input: PathBuf,

    /// 出力 PNG
    #[arg(short, long, default_value = "hicmap.png")]
    output: PathBuf,

    /// 横方向の画素数 (= 解像度)。1 画素が 1 bin を下回る指定は拒否する
    #[arg(long, default_value_t = 4000)]
    nx: usize,

    /// 表示するゲノム範囲の開始 (bin)。省略でデータ全域
    #[arg(long)]
    x_start: Option<f64>,

    /// 表示するゲノム範囲の終了 (bin, この値は含まない)。省略でデータ全域
    #[arg(long)]
    x_end: Option<f64>,

    /// 表示する最大の bin 間距離 (bin, この値を含む)。省略でデータ内の最大距離
    #[arg(long)]
    max_distance: Option<f64>,

    /// 染色体長ファイル (chr<TAB>bp)。指定すると軸を bp で目盛り、
    /// 染色体境界を跨ぐ接触を除外して境界線を引く
    #[arg(long)]
    lengths: Option<PathBuf>,

    /// bp/bin。省略時は --lengths と最大 bin から推定する
    #[arg(long)]
    bin_size: Option<u32>,

    /// 染色体境界を跨ぐ接触も残す (--lengths 指定時のみ意味を持つ)
    #[arg(long)]
    keep_inter: bool,

    /// log10 をとらず生の score 総和で色をつける
    #[arg(long)]
    no_log: bool,

    /// 上側をクリップする分位点 (0..1)
    #[arg(long, default_value_t = 0.99)]
    trim_quantile: f64,

    /// 上側のクリップを行わない
    #[arg(long)]
    no_trim: bool,

    /// 図の見出し。省略で入力ファイル名
    #[arg(long)]
    title: Option<String>,

    /// ラベルに使うフォントファミリ。省略で使えるものを自動選択
    #[arg(long)]
    font: Option<String>,

    /// 集計したグリッドを TSV へ書き出す (px py x y raw value)
    #[arg(long)]
    dump_grid: Option<PathBuf>,

    /// 総画素数の上限。グリッドは 8 バイト/画素を使う
    #[arg(long, default_value_t = 40_000_000)]
    max_pixels: usize,

    /// 1 画素が 1 bin 未満でも描く (モアレを承知のうえで)
    #[arg(long)]
    allow_aliasing: bool,
}

fn main() -> ExitCode {
    match run(Args::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("エラー: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    if !(0.0..=1.0).contains(&args.trim_quantile) {
        return Err(format!(
            "--trim-quantile は 0..1 で指定してください (指定値: {})",
            args.trim_quantile
        )
        .into());
    }
    if args.nx == 0 {
        return Err("--nx は 1 以上で指定してください".into());
    }

    // --- 1 パス目: データの広がりを調べる (範囲の自動決定と bin_size 推定のため) ---
    eprintln!("[1/3] データの広がりを走査中: {}", args.input.display());
    let extent = grid::scan_extent(&args.input)?;
    eprintln!(
        "      bin {}..{}, 最大距離 {} bin",
        extent.bin_min, extent.bin_max, extent.dist_max
    );

    // --- bin_size と染色体テーブル ---
    let bin_size = resolve_bin_size(&args, &extent)?;
    let chroms = match (&args.lengths, bin_size) {
        (Some(path), Some(bs)) => {
            let idx = ChromIndex::load(path, bs)?;
            eprintln!(
                "      染色体 {} 本 / 全 {} bin ({} bp/bin)",
                idx.chroms().len(),
                idx.total_bins(),
                bs
            );
            if extent.bin_max > idx.total_bins() {
                eprintln!(
                    "警告: データの最大 bin {} が染色体テーブルの全 bin 数 {} を超えています。\n\
                     --bin-size が実際と食い違っている可能性があります",
                    extent.bin_max,
                    idx.total_bins()
                );
            }
            Some(idx)
        }
        _ => None,
    };

    // --- 表示範囲を決める ---
    // xr = (bin1 + bin2) / 2 の最大は bin_max なので、半開区間の上端は +1 して含める。
    let x0 = args.x_start.unwrap_or(extent.bin_min as f64);
    let x1 = args.x_end.unwrap_or(extent.bin_max as f64 + 1.0);
    if x1 <= x0 {
        return Err(format!("表示範囲が空です: [{x0}, {x1})").into());
    }
    let max_distance = args.max_distance.unwrap_or(extent.dist_max as f64).max(1.0);

    let spec = GridSpec::new(x0, x1, max_distance, args.nx);
    report_geometry(&spec, bin_size);

    // --- 画素数の検査 ---
    if let Some(a) = spec.aliasing() {
        if args.allow_aliasing {
            eprintln!("警告: {a}");
        } else {
            return Err(format!("{a}\n(承知のうえで描くなら --allow-aliasing)").into());
        }
    }
    if spec.pixels() > args.max_pixels {
        return Err(format!(
            "画素数 {} が上限 {} を超えます (グリッドに {:.1} GB 必要)。\n\
             --nx を下げるか --max-pixels を上げてください",
            spec.pixels(),
            args.max_pixels,
            spec.pixels() as f64 * 8.0 / 1e9
        )
        .into());
    }

    // --- 2 パス目: グリッドへ集計 ---
    eprintln!("[2/3] グリッドへ集計中 ({} x {} 画素)", spec.nx, spec.ny);
    let filter = if args.keep_inter {
        None
    } else {
        chroms.as_ref()
    };
    let g = grid::build(&args.input, spec, filter)?;
    eprintln!(
        "      集計 {} 行 / 範囲外 {} 行{}",
        g.counted,
        g.out_of_range,
        if g.inter_chrom > 0 {
            format!(" / 染色体を跨ぎ除外 {} 行", g.inter_chrom)
        } else {
            String::new()
        }
    );
    if g.counted == 0 {
        return Err("表示範囲に接触が 1 つもありません。--x-start / --x-end / --max-distance を確認してください".into());
    }

    // --- 値変換と描画 ---
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
        "      値のある画素 {} / {} ({:.1}%), 色域 {:.3}..{:.3}",
        values.filled,
        spec.pixels(),
        values.filled as f64 / spec.pixels() as f64 * 100.0,
        values.vmin,
        values.vmax
    );

    if let Some(path) = &args.dump_grid {
        dump_grid(path, &g, &values)?;
        eprintln!("      グリッドを書き出しました: {}", path.display());
    }

    eprintln!("[3/3] 描画中: {}", args.output.display());
    let font = font::pick(args.font.as_deref());
    eprintln!(
        "      フォント: {}{}",
        font.family,
        if font.cjk {
            ""
        } else {
            " (日本語非対応のため英語ラベル)"
        }
    );
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
    eprintln!("完了: {}", args.output.display());
    Ok(())
}

/// bp/bin を決める。明示指定 > 染色体長からの推定 > 不明 (軸は bin 単位) の順。
fn resolve_bin_size(
    args: &Args,
    extent: &grid::Extent,
) -> Result<Option<u32>, Box<dyn std::error::Error>> {
    if let Some(bs) = args.bin_size {
        if bs == 0 {
            return Err("--bin-size は 1 以上で指定してください".into());
        }
        return Ok(Some(bs));
    }
    let Some(path) = &args.lengths else {
        return Ok(None); // 染色体長が無ければ bp へ換算できない
    };
    let total_bp = ChromIndex::total_bp(path)?;
    let bs = chrom::estimate_bin_size(total_bp, extent.bin_max);
    eprintln!(
        "      bin_size 推定: {bs} bp/bin (総 {total_bp} bp / 最大 bin {})",
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
        "      表示範囲 bin {:.0}..{:.0}{} / 最大距離 {:.0} bin{}",
        spec.x0,
        spec.x1,
        bp(spec.x1 - spec.x0),
        spec.max_distance,
        bp(spec.max_distance)
    );
    let per_px = match bin_size {
        Some(bs) => format!(" = {:.0} bp/画素", spec.dx * bs as f64),
        None => String::new(),
    };
    eprintln!("      解像度 {:.3} bin/画素{}", spec.dx, per_px);
}

fn default_title(input: &Path) -> String {
    let name = input.file_stem().map_or_else(
        || input.display().to_string(),
        |s| s.to_string_lossy().into_owned(),
    );
    format!("Hi-C contact map: {name}")
}

/// 集計結果を TSV で書き出す。値のない画素は行を出さない。
/// x, y は画素中心の座標 (bin 単位)。y は回転座標 yr (= 距離 / 2)。
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
