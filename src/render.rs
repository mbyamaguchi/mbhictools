//! グリッドの値変換・配色と、PNG への描画。

use std::path::Path;

use plotters::coord::Shift;
use plotters::prelude::*;
use plotters::style::text_anchor::{HPos, Pos, VPos};

use crate::chrom::ChromIndex;
use crate::font::Font;
use crate::grid::Grid;

/// 図の余白 (画素)。作図領域がちょうど nx * ny 画素になるよう、
/// 図全体の大きさをこれらから逆算する。
const MARGIN_TOP: u32 = 46;
const LABEL_LEFT: u32 = 86;
const LABEL_BOTTOM: u32 = 54;
const CBAR_GAP: u32 = 24;
const CBAR_WIDTH: u32 = 18;
const CBAR_LABEL: u32 = 66;
const MARGIN_RIGHT: u32 = 16;

/// 作図領域の右に確保する幅 (カラーバーぶん)。
const RIGHT_TOTAL: u32 = CBAR_GAP + CBAR_WIDTH + CBAR_LABEL + MARGIN_RIGHT;

/// 値変換の指定。0 の画素は常に「データなし」として背景色のままにする。
#[derive(Debug, Clone, Copy)]
pub struct Scale {
    /// log10 をとる。Hi-C の score は数桁にわたるので既定で有効。
    pub log: bool,
    /// 上側をこの分位点でクリップする。少数の極端に強い画素に
    /// 色域を占有されるのを防ぐ。`None` でクリップなし。
    pub trim_quantile: Option<f64>,
}

impl Default for Scale {
    fn default() -> Self {
        Scale {
            log: true,
            trim_quantile: Some(0.99),
        }
    }
}

/// 変換後の画素値。`v[i]` は `grid.cells()[i]` に対応し、`None` はデータなし。
#[derive(Debug, Clone)]
pub struct Values {
    pub v: Vec<Option<f64>>,
    pub vmin: f64,
    pub vmax: f64,
    /// 値のある画素数。
    pub filled: usize,
}

/// 0 を「データなし」にし、log10 をとり、上側分位点でクリップする。
pub fn transform(grid: &Grid, scale: &Scale) -> Values {
    let convert = |c: u64| -> Option<f64> {
        if c == 0 {
            return None;
        }
        let x = c as f64;
        Some(if scale.log { x.log10() } else { x })
    };

    let mut v: Vec<Option<f64>> = grid.cells().iter().map(|&c| convert(c)).collect();

    let mut present: Vec<f64> = v.iter().flatten().copied().collect();
    if present.is_empty() {
        return Values {
            v,
            vmin: 0.0,
            vmax: 1.0,
            filled: 0,
        };
    }
    present.sort_by(f64::total_cmp);

    let vmin = present[0];
    let mut vmax = present[present.len() - 1];
    if let Some(q) = scale.trim_quantile {
        let cap = quantile_type7(&present, q);
        for x in v.iter_mut().flatten() {
            *x = x.min(cap);
        }
        vmax = cap;
    }
    // 全画素が同値でも色域が潰れないようにする。
    if vmax <= vmin {
        vmax = vmin + 1.0;
    }
    let filled = present.len();
    Values {
        v,
        vmin,
        vmax,
        filled,
    }
}

/// R の `quantile(type = 7)` (既定) と同じ定義。`sorted` は昇順で非空であること。
fn quantile_type7(sorted: &[f64], p: f64) -> f64 {
    assert!(!sorted.is_empty());
    assert!(
        (0.0..=1.0).contains(&p),
        "分位点は 0..1 の範囲で指定すること"
    );
    let h = (sorted.len() - 1) as f64 * p;
    let lo = h.floor() as usize;
    let hi = (lo + 1).min(sorted.len() - 1);
    sorted[lo] + (h - lo as f64) * (sorted[hi] - sorted[lo])
}

/// 配色。R 版の既定と同じ白 -> 濃赤 (ColorBrewer Reds 系)。
#[derive(Debug, Clone)]
pub struct Palette {
    ramp: Vec<RGBColor>,
}

const DEFAULT_STOPS: [(u8, u8, u8); 6] = [
    (0xFF, 0xFF, 0xFF),
    (0xFF, 0xF5, 0xF0),
    (0xFC, 0xBB, 0xA1),
    (0xFB, 0x6A, 0x4A),
    (0xCB, 0x18, 0x1D),
    (0x67, 0x00, 0x0D),
];

impl Default for Palette {
    fn default() -> Self {
        Palette::from_stops(&DEFAULT_STOPS)
    }
}

impl Palette {
    /// 色見本を線形補間して 256 段のランプにする。
    pub fn from_stops(stops: &[(u8, u8, u8)]) -> Self {
        assert!(stops.len() >= 2, "色見本は 2 色以上必要");
        const N: usize = 256;
        let ramp = (0..N)
            .map(|i| {
                let t = i as f64 / (N - 1) as f64 * (stops.len() - 1) as f64;
                let lo = (t.floor() as usize).min(stops.len() - 2);
                let f = t - lo as f64;
                let (a, b) = (stops[lo], stops[lo + 1]);
                let mix = |x: u8, y: u8| (x as f64 + (y as f64 - x as f64) * f).round() as u8;
                RGBColor(mix(a.0, b.0), mix(a.1, b.1), mix(a.2, b.2))
            })
            .collect();
        Palette { ramp }
    }

    /// [0, 1] に正規化された値を色へ。範囲外は両端に丸める。
    pub fn at(&self, t: f64) -> RGBColor {
        let i = (t.clamp(0.0, 1.0) * (self.ramp.len() - 1) as f64).round() as usize;
        self.ramp[i]
    }
}

/// 図の見出しや軸の単位。
#[derive(Debug, Clone)]
pub struct Labels {
    pub title: String,
    /// bp 単位で目盛る場合の bp/bin。`None` なら bin 単位で目盛る。
    pub bin_size: Option<u32>,
    /// 凡例の見出し (例: `log10(score)`)。
    pub legend: String,
    /// 描画に使うフォント。日本語を持たないフォントなら軸名も英語にする。
    pub font: Font,
}

impl Labels {
    /// 横軸の名前。
    fn x_desc(&self, unit: &str) -> String {
        if self.font.cjk {
            format!("ゲノム位置 ({unit})")
        } else {
            format!("Genomic position ({unit})")
        }
    }

    /// 縦軸の名前。
    fn y_desc(&self, unit: &str) -> String {
        if self.font.cjk {
            format!("相互作用距離 ({unit})")
        } else {
            format!("Interaction distance ({unit})")
        }
    }
}

/// 背景色 (データのない画素の色)。
const BACKGROUND: RGBColor = WHITE;

/// グリッドを PNG へ描く。作図領域はちょうど `nx * ny` 画素になる。
pub fn render(
    out: &Path,
    grid: &Grid,
    values: &Values,
    palette: &Palette,
    labels: &Labels,
    chroms: Option<&ChromIndex>,
) -> Result<(), Box<dyn std::error::Error>> {
    let spec = grid.spec;
    let (nx, ny) = (spec.nx as u32, spec.ny as u32);
    let width = LABEL_LEFT + nx + RIGHT_TOTAL;
    let height = MARGIN_TOP + ny + LABEL_BOTTOM;

    let root = BitMapBackend::new(out, (width, height)).into_drawing_area();
    root.fill(&BACKGROUND)?;

    // 軸の単位。bin_size があれば bp、なければ bin のまま。
    let unit = labels.bin_size.map_or(1.0, |b| b as f64);
    let x_range = (spec.x0 * unit)..(spec.x1 * unit);
    // 縦軸は回転座標 yr ではなく実際の距離 (= 2 * yr) で目盛る。
    let y_range = 0.0..(spec.max_distance * unit);

    let mut chart = ChartBuilder::on(&root)
        .margin_top(MARGIN_TOP)
        .margin_bottom(0)
        .margin_left(0)
        .margin_right(RIGHT_TOTAL)
        .x_label_area_size(LABEL_BOTTOM)
        .y_label_area_size(LABEL_LEFT)
        .build_cartesian_2d(x_range.clone(), y_range.clone())?;

    // 逆算した図の大きさどおりなら、作図領域は画素と 1 対 1 に対応するはず。
    let area = chart.plotting_area().dim_in_pixel();
    debug_assert_eq!(area, (nx, ny), "作図領域がグリッドと一致していません");

    let axis_unit = if labels.bin_size.is_some() {
        "bp"
    } else {
        "bin"
    };
    let family = labels.font.family.as_str();
    chart
        .configure_mesh()
        .disable_mesh()
        .x_desc(labels.x_desc(axis_unit))
        .y_desc(labels.y_desc(axis_unit))
        .x_label_formatter(&|v| format_pos(*v, labels.bin_size.is_some()))
        .y_label_formatter(&|v| format_pos(*v, labels.bin_size.is_some()))
        .axis_desc_style((family, 15))
        .label_style((family, 12))
        .draw()?;

    // グリッドを RGB バッファへ焼いてから一括で貼る。画素ごとに描くより速く、
    // かつ作図領域の画素と 1 対 1 で対応することが保証できる。
    let buf = rasterize(grid, values, palette);
    let image = BitMapElement::with_owned_buffer((x_range.start, y_range.end), (nx, ny), buf)
        .ok_or("ラスタバッファの大きさが作図領域と一致しません")?;
    chart.plotting_area().draw(&image)?;

    if let Some(idx) = chroms {
        draw_chrom_boundaries(&mut chart, grid, idx, unit)?;
    }

    root.draw(&Rectangle::new(
        [
            (LABEL_LEFT as i32, MARGIN_TOP as i32),
            ((LABEL_LEFT + nx) as i32, (MARGIN_TOP + ny) as i32),
        ],
        BLACK.stroke_width(1),
    ))?;

    draw_title(&root, &labels.title, family, width)?;
    draw_colorbar(&root, values, palette, &labels.legend, family, nx, ny)?;

    root.present()?;
    Ok(())
}

/// グリッドを RGB (3 バイト / 画素) のバッファへ。行 0 が図の上端になるよう
/// py を上下反転する (py = 0 は対角線 = 図の下端)。
fn rasterize(grid: &Grid, values: &Values, palette: &Palette) -> Vec<u8> {
    let spec = grid.spec;
    let span = values.vmax - values.vmin;
    let mut buf = vec![0u8; spec.nx * spec.ny * 3];

    for py in 0..spec.ny {
        let row = spec.ny - 1 - py; // 上下反転
        for px in 0..spec.nx {
            let color = match values.v[py * spec.nx + px] {
                Some(x) => palette.at((x - values.vmin) / span),
                None => BACKGROUND,
            };
            let o = (row * spec.nx + px) * 3;
            buf[o] = color.0;
            buf[o + 1] = color.1;
            buf[o + 2] = color.2;
        }
    }
    buf
}

/// 染色体境界を縦線で示す。回転座標では境界を跨ぐ領域は三角形になるが、
/// 境界の位置を示すには縦線で足りる。
fn draw_chrom_boundaries<DB: DrawingBackend>(
    chart: &mut ChartContext<
        '_,
        DB,
        Cartesian2d<plotters::coord::types::RangedCoordf64, plotters::coord::types::RangedCoordf64>,
    >,
    grid: &Grid,
    idx: &ChromIndex,
    unit: f64,
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    let spec = grid.spec;
    let style = RGBColor(0x40, 0x40, 0x40).mix(0.55).stroke_width(1);
    for b in idx.boundaries() {
        let x = b as f64;
        if x <= spec.x0 || x >= spec.x1 {
            continue;
        }
        chart.draw_series(LineSeries::new(
            [(x * unit, 0.0), (x * unit, spec.max_distance * unit)],
            style,
        ))?;
    }
    Ok(())
}

fn draw_title<DB: DrawingBackend>(
    root: &DrawingArea<DB, Shift>,
    title: &str,
    family: &str,
    width: u32,
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    let style = TextStyle::from((family, 17).into_font()).color(&BLACK);
    root.draw_text(
        title,
        &style.pos(Pos::new(HPos::Center, VPos::Center)),
        (width as i32 / 2, 22),
    )?;
    Ok(())
}

/// 作図領域の右に縦向きのカラーバーを描く。
fn draw_colorbar<DB: DrawingBackend>(
    root: &DrawingArea<DB, Shift>,
    values: &Values,
    palette: &Palette,
    legend: &str,
    family: &str,
    nx: u32,
    ny: u32,
) -> Result<(), Box<dyn std::error::Error>>
where
    DB::ErrorType: 'static,
{
    // カラーバーが潰れないよう、作図領域が低くても最低限の高さを確保する。
    let bar_h = ny.max(120);
    let top = MARGIN_TOP as i32 + (ny as i32 - bar_h as i32) / 2;
    let left = (LABEL_LEFT + nx + CBAR_GAP) as i32;
    let right = left + CBAR_WIDTH as i32;

    // 1 画素ずつの横帯で階調をつくる (上が vmax)。
    for i in 0..bar_h {
        let t = 1.0 - i as f64 / (bar_h - 1) as f64;
        let y = top + i as i32;
        root.draw(&Rectangle::new(
            [(left, y), (right, y + 1)],
            palette.at(t).filled(),
        ))?;
    }
    root.draw(&Rectangle::new(
        [(left, top), (right, top + bar_h as i32)],
        BLACK.stroke_width(1),
    ))?;

    let text = TextStyle::from((family, 12).into_font()).color(&BLACK);
    let ticks = [
        (top, values.vmax),
        (top + bar_h as i32 / 2, (values.vmin + values.vmax) / 2.0),
        (top + bar_h as i32, values.vmin),
    ];
    for (y, v) in ticks {
        root.draw(&PathElement::new(
            [(right, y), (right + 4, y)],
            BLACK.stroke_width(1),
        ))?;
        root.draw_text(
            &format!("{v:.2}"),
            &text.pos(Pos::new(HPos::Left, VPos::Center)),
            (right + 7, y),
        )?;
    }

    let label = TextStyle::from((family, 13).into_font()).color(&BLACK);
    root.draw_text(
        legend,
        &label.pos(Pos::new(HPos::Center, VPos::Bottom)),
        (left + CBAR_WIDTH as i32 / 2, top - 8),
    )?;
    Ok(())
}

/// 座標を短く読める形に。bp なら kb / Mb へ丸める。
fn format_pos(v: f64, as_bp: bool) -> String {
    if !as_bp {
        return format!("{v:.0}");
    }
    let a = v.abs();
    if a >= 1e6 {
        format!("{:.1} Mb", v / 1e6)
    } else if a >= 1e3 {
        format!("{:.0} kb", v / 1e3)
    } else {
        format!("{v:.0}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::{Grid, GridSpec};

    fn grid_from(cells: Vec<u64>, nx: usize) -> Grid {
        let ny = cells.len() / nx;
        // 画素数だけを合わせた最小の spec (値変換のテストは幾何に依らない)。
        let spec = GridSpec {
            x0: 0.0,
            x1: nx as f64,
            dx: 1.0,
            nx,
            ny,
            ymax: ny as f64,
            max_distance: ny as f64 * 2.0,
        };
        Grid::from_parts(spec, cells)
    }

    #[test]
    fn treats_zero_as_missing() {
        let g = grid_from(vec![0, 10, 0, 100], 2);
        let v = transform(
            &g,
            &Scale {
                log: false,
                trim_quantile: None,
            },
        );
        assert_eq!(v.v, vec![None, Some(10.0), None, Some(100.0)]);
        assert_eq!(v.filled, 2);
    }

    #[test]
    fn takes_log10() {
        let g = grid_from(vec![1, 10, 100, 1000], 2);
        let v = transform(
            &g,
            &Scale {
                log: true,
                trim_quantile: None,
            },
        );
        assert_eq!(v.v, vec![Some(0.0), Some(1.0), Some(2.0), Some(3.0)]);
        assert_eq!(v.vmin, 0.0);
        assert_eq!(v.vmax, 3.0);
    }

    #[test]
    fn clips_at_the_upper_quantile() {
        // 1..=10 の 0.5 分位 (type 7) は 5.5。それを超える値は 5.5 に丸められる。
        let g = grid_from((1..=10).collect(), 10);
        let v = transform(
            &g,
            &Scale {
                log: false,
                trim_quantile: Some(0.5),
            },
        );
        assert_eq!(v.vmax, 5.5);
        assert_eq!(v.v[9], Some(5.5), "最大値 10 はクリップされる");
        assert_eq!(v.v[0], Some(1.0), "下側は変わらない");
    }

    /// R の `quantile(1:10, p)` (type 7) と一致すること。
    #[test]
    fn quantile_matches_r_type7() {
        let x: Vec<f64> = (1..=10).map(|i| i as f64).collect();
        assert_eq!(quantile_type7(&x, 0.0), 1.0);
        assert_eq!(quantile_type7(&x, 1.0), 10.0);
        assert_eq!(quantile_type7(&x, 0.5), 5.5);
        assert!((quantile_type7(&x, 0.99) - 9.91).abs() < 1e-9);
        assert!((quantile_type7(&x, 0.25) - 3.25).abs() < 1e-9);
    }

    #[test]
    fn handles_an_empty_grid() {
        let g = grid_from(vec![0, 0, 0, 0], 2);
        let v = transform(&g, &Scale::default());
        assert_eq!(v.filled, 0);
        assert!(v.vmax > v.vmin, "色域が潰れないこと");
    }

    #[test]
    fn keeps_a_usable_range_when_all_pixels_are_equal() {
        let g = grid_from(vec![5, 5, 5, 5], 2);
        let v = transform(
            &g,
            &Scale {
                log: false,
                trim_quantile: Some(0.99),
            },
        );
        assert!(v.vmax > v.vmin);
    }

    #[test]
    fn palette_spans_from_first_to_last_stop() {
        let p = Palette::default();
        assert_eq!(p.at(0.0), RGBColor(0xFF, 0xFF, 0xFF), "下端は白");
        assert_eq!(p.at(1.0), RGBColor(0x67, 0x00, 0x0D), "上端は濃赤");
        assert_eq!(p.at(-5.0), p.at(0.0), "範囲外は丸める");
        assert_eq!(p.at(5.0), p.at(1.0));
    }

    #[test]
    fn rasterizes_with_the_diagonal_at_the_bottom() {
        // 2x2。py=0 (対角) が図の下段、py=1 が上段に来ること。
        let g = grid_from(vec![1, 1, 1000, 1000], 2);
        let v = transform(
            &g,
            &Scale {
                log: true,
                trim_quantile: None,
            },
        );
        let p = Palette::default();
        let buf = rasterize(&g, &v, &p);

        assert_eq!(buf.len(), 2 * 2 * 3);
        let top_left = (buf[0], buf[1], buf[2]);
        let bottom_left = (buf[6], buf[7], buf[8]);
        assert_eq!(
            top_left,
            (0x67, 0x00, 0x0D),
            "上段 (py=1) は値 1000 -> 濃赤"
        );
        assert_eq!(bottom_left, (0xFF, 0xFF, 0xFF), "下段 (py=0) は値 1 -> 白");
    }

    #[test]
    fn missing_pixels_get_the_background() {
        let g = grid_from(vec![0, 5], 2);
        let v = transform(
            &g,
            &Scale {
                log: false,
                trim_quantile: None,
            },
        );
        let buf = rasterize(&g, &v, &Palette::default());
        assert_eq!(
            (buf[0], buf[1], buf[2]),
            (0xFF, 0xFF, 0xFF),
            "0 の画素は背景色"
        );
    }

    #[test]
    fn formats_positions() {
        assert_eq!(format_pos(2_500_000.0, true), "2.5 Mb");
        assert_eq!(format_pos(12_000.0, true), "12 kb");
        assert_eq!(format_pos(500.0, true), "500");
        assert_eq!(format_pos(4000.0, false), "4000");
    }
}
