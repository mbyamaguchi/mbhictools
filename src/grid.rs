//! 接触データを回転座標の表示グリッド (= 画素) へ集計する。
//!
//! # なぜ回転座標か
//!
//! 入力は上三角成分のみ、かつ bin 間距離が上限 (既定では 1 Mbp) 以内に限られる。
//! そのため素直な正方形の contact map を描くと、面積の大半が定義上まったく
//! データのない領域になる。対角線に沿う細い帯だけが意味を持つ。
//!
//! そこで各接触を 45 度回した座標
//!
//! ```text
//!   xr = (bin1 + bin2) / 2   ゲノム上の位置 (2 つの bin の中点)
//!   yr = (bin2 - bin1) / 2   相互作用距離の半分
//! ```
//!
//! へ写す。すると帯は「底辺 = ゲノム、高さ = 最大距離 / 2」の三角形になり、
//! 画素が無駄にならない。
//!
//! # 画素数の決め方 (このモジュールの肝)
//!
//! `dx` を 1 画素あたりの bin 数とする。画素は正方形 (`dy == dx`) にとる。
//! これはデータ単位で等方、つまり図の縦横比が実際の距離を正しく表すということ。
//!
//! ここで `nx` を上げ過ぎると破綻する。回転後の接触点は連続ではなく格子上にあり、
//! しかもその格子は市松模様だからである。`bin1, bin2` が整数なので
//! `xr, yr` は 0.5 の倍数だが、任意の組が現れるわけではない:
//! `xr + yr = bin2` は必ず整数になるため、`xr + yr` が整数の点しか存在しない。
//!
//! この格子に対して 1 画素が `dx` 四方のとき、画素に入る格子点の数は:
//!
//! ```text
//!   dx = 1     : どの画素もちょうど 2 点 (xr,yr が共に整数の点と、共に半整数の点)
//!   dx = 0.5   : (i+j) が偶数の画素だけ 1 点、奇数の画素は 0 点 → 完全な市松模様
//!   dx >= 1    : 平均 2*dx^2 点。dx が大きいほど画素間のばらつきは相対的に小さい
//! ```
//!
//! つまり `dx = 1 bin/画素` が厳密な下限で、これを下回ると空画素が交互に並ぶ
//! モアレが出る。データが増えたわけでもないのに図が粗く見えるという最悪の失敗で、
//! しかも「解像度を上げた」つもりで起きる。[`GridSpec::aliasing`] がこれを検出する。

use std::path::Path;
use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};

use crate::chrom::ChromIndex;
use crate::contact;

/// 表示グリッドの幾何。すべて bin 単位 (bp ではない)。
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GridSpec {
    /// 表示するゲノム範囲 [x0, x1) (xr の範囲)。
    pub x0: f64,
    pub x1: f64,
    /// 1 画素あたりの bin 数。縦横とも同じ (正方形画素)。
    pub dx: f64,
    /// 横方向の画素数。
    pub nx: usize,
    /// 縦方向の画素数。`ymax / dx` から決まる。
    pub ny: usize,
    /// 三角形の高さ = max_distance / 2 (回転座標での縦の上限)。
    pub ymax: f64,
    /// 表示する最大の bin 間距離 |bin2 - bin1|。
    pub max_distance: f64,
}

/// 画素が細か過ぎるときの診断。`nx` を上げ過ぎた結果として起きる。
#[derive(Debug, Clone, PartialEq)]
pub struct Aliasing {
    /// 現在の 1 画素あたりの bin 数 (< 1)。
    pub dx: f64,
    /// モアレを起こさない `nx` の上限。
    pub max_nx: usize,
}

impl std::fmt::Display for Aliasing {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "1 画素が {:.3} bin しかありません (< 1 bin/画素)。\n\
             回転後の接触点は `xr + yr` が整数の格子上にしかないため、\n\
             この解像度では空の画素が交互に並ぶ市松模様のモアレが出ます。\n\
             この範囲では --nx {} 以下を指定してください",
            self.dx, self.max_nx
        )
    }
}

/// データの広がり。`x_range` / `max_distance` を省略したときの自動決定に使う。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Extent {
    pub bin_min: u32,
    pub bin_max: u32,
    pub dist_max: u32,
}

/// 集計済みグリッド。`cells[py * nx + px]` に score の総和が入る。
#[derive(Debug, Clone)]
pub struct Grid {
    pub spec: GridSpec,
    cells: Vec<u64>,
    /// グリッドに実際に加算された接触の行数。
    pub counted: u64,
    /// 表示範囲外だった行数。
    pub out_of_range: u64,
    /// 染色体境界を跨ぐために除外した行数 (フィルタ無効時は 0)。
    pub inter_chrom: u64,
}

impl GridSpec {
    /// 表示範囲と横画素数から幾何を決める。
    ///
    /// `ny` は正方形画素 (`dy == dx`) になるよう `ymax / dx` から導く。縦の画素数を
    /// 独立に指定できないのは意図的で、そうしないと図の縦横比が距離を偽ってしまう。
    pub fn new(x0: f64, x1: f64, max_distance: f64, nx: usize) -> Self {
        assert!(nx >= 1, "nx は 1 以上であること");
        assert!(x1 > x0, "表示範囲が空です: [{x0}, {x1})");
        assert!(max_distance > 0.0, "max_distance は正であること");

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

    /// 画素が細か過ぎて回転格子のモアレが出るなら、その診断を返す。
    pub fn aliasing(&self) -> Option<Aliasing> {
        if self.dx >= 1.0 {
            return None;
        }
        // dx >= 1 すなわち (x1 - x0) / nx >= 1 となる最大の nx。
        let max_nx = (self.x1 - self.x0).floor().max(1.0) as usize;
        Some(Aliasing {
            dx: self.dx,
            max_nx,
        })
    }

    /// 総画素数。
    pub fn pixels(&self) -> usize {
        self.nx * self.ny
    }

    /// 接触 (bin1, bin2) が落ちる画素。表示範囲外なら `None`。
    ///
    /// ゲノム位置は半開区間 [x0, x1)、距離は閉区間 [0, max_distance] で扱う。
    /// 距離だけ閉じているのは「距離が 1 Mbp 以下・以内」という入力の仕様を
    /// そのまま写すため。半開にすると上限ちょうどの接触 (実データでは
    /// `bin2 - bin1 == 5000` の行) を丸ごと落としてしまう。
    #[inline]
    fn pixel_of(&self, bin1: u32, bin2: u32) -> Option<(usize, usize)> {
        let (b1, b2) = (bin1 as f64, bin2 as f64);
        let xr = (b1 + b2) * 0.5;
        let yr = (b2 - b1) * 0.5;
        if xr < self.x0 || xr >= self.x1 || yr < 0.0 || yr > self.ymax {
            return None;
        }
        let px = ((xr - self.x0) / self.dx) as usize;
        // 除算の丸めや、yr == ymax ちょうどの点が最上段を 1 つ踏み越えうる。
        let py = ((yr / self.dx) as usize).min(self.ny - 1);
        if px >= self.nx {
            return None;
        }
        Some((px, py))
    }

    /// 画素の中心のゲノム座標 (bin 単位)。
    pub fn pixel_center(&self, px: usize, py: usize) -> (f64, f64) {
        (
            self.x0 + (px as f64 + 0.5) * self.dx,
            (py as f64 + 0.5) * self.dx,
        )
    }
}

impl Grid {
    /// 集計済みの画素値からグリッドを組み立てる。走査の統計は 0 になる。
    pub fn from_parts(spec: GridSpec, cells: Vec<u64>) -> Self {
        assert_eq!(
            cells.len(),
            spec.pixels(),
            "cells の長さが画素数と一致しません"
        );
        Grid {
            spec,
            cells,
            counted: 0,
            out_of_range: 0,
            inter_chrom: 0,
        }
    }

    /// `cells[py * nx + px]`。
    pub fn get(&self, px: usize, py: usize) -> u64 {
        self.cells[py * self.spec.nx + px]
    }

    pub fn cells(&self) -> &[u64] {
        &self.cells
    }

    /// 0 でない画素の値だけを集めたもの (分位点の計算用)。
    pub fn nonzero(&self) -> Vec<u64> {
        self.cells.iter().copied().filter(|&v| v > 0).collect()
    }
}

/// ファイルを 1 パス走査してデータの広がりを調べる。
///
/// `x_range` や `max_distance` を省略したときに、グリッドを確保する前の
/// 範囲決定に使う (score は見ないので集計本体より軽い)。
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

/// ファイルを 1 パス走査して `spec` のグリッドへ集計する。
///
/// メモリはグリッドぶん (nx * ny * 8 バイト) のみで、行数には依存しない。
/// score は整数なので u64 のアトミック加算で厳密かつ並列に足し込める。
///
/// `chrom_filter` を渡すと、染色体境界を跨ぐ接触を除外する。
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
        // ゲノム 1000 bin を 100 画素 -> 10 bin/画素。最大距離 200 bin -> 高さ 100 bin -> 10 画素。
        let s = GridSpec::new(0.0, 1000.0, 200.0, 100);
        assert_eq!(s.dx, 10.0);
        assert_eq!(s.ymax, 100.0);
        assert_eq!(s.ny, 10);
        assert_eq!(s.pixels(), 1000);
    }

    #[test]
    fn rounds_ny_up_to_cover_the_top() {
        // ymax / dx = 105 / 10 = 10.5 -> 端を切り落とさないよう 11 画素。
        let s = GridSpec::new(0.0, 1000.0, 210.0, 100);
        assert_eq!(s.ny, 11);
    }

    #[test]
    fn maps_contacts_to_rotated_pixels() {
        let s = GridSpec::new(0.0, 1000.0, 200.0, 100);
        // bin1=bin2=5 -> xr=5, yr=0 -> 対角 (px=0, py=0)
        assert_eq!(s.pixel_of(5, 5), Some((0, 0)));
        // bin1=100, bin2=120 -> xr=110, yr=10 -> px=11, py=1
        assert_eq!(s.pixel_of(100, 120), Some((11, 1)));
        // xr が右端 (x1) 以上は範囲外
        assert_eq!(s.pixel_of(1000, 1000), None);
        // 距離が上限ちょうど (yr == ymax) の接触は最上段に含める (閉区間)
        assert_eq!(s.pixel_of(0, 200), Some((10, 9)));
        // 上限を超えたものは範囲外
        assert_eq!(s.pixel_of(0, 202), None);
    }

    #[test]
    fn pixel_center_matches_mapping() {
        let s = GridSpec::new(0.0, 1000.0, 200.0, 100);
        assert_eq!(s.pixel_center(0, 0), (5.0, 5.0));
        assert_eq!(s.pixel_center(11, 1), (115.0, 15.0));
    }

    /// dx = 1 のとき、回転格子はどの画素にもちょうど 2 点を与える (モジュール doc の主張)。
    #[test]
    fn unit_pixels_receive_exactly_two_lattice_points_each() {
        let s = GridSpec::new(0.0, 40.0, 20.0, 40);
        assert_eq!(s.dx, 1.0);
        assert_eq!(
            s.aliasing(),
            None,
            "dx = 1 は下限ちょうどで、モアレは出ない"
        );

        let mut count = vec![0usize; s.pixels()];
        for b1 in 0..=60u32 {
            for b2 in b1..=60u32 {
                if let Some((px, py)) = s.pixel_of(b1, b2) {
                    count[py * s.nx + px] += 1;
                }
            }
        }
        // 三角形の内側 (端の切れる画素を避ける) はすべて 2 点ちょうど。
        for py in 0..5 {
            for px in 10..30 {
                assert_eq!(count[py * s.nx + px], 2, "画素 ({px}, {py})");
            }
        }
    }

    /// dx < 1 では格子点の入らない画素が生じる (= 市松模様のモアレ)。
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
            "dx = 0.5 では約半数の画素が空になるはず (実際 {empty} / {})",
            s.pixels()
        );
    }

    #[test]
    fn reports_aliasing_with_a_usable_nx_limit() {
        let s = GridSpec::new(0.0, 1000.0, 200.0, 4000);
        let a = s.aliasing().expect("dx = 0.25 なのでモアレが出る");
        assert_eq!(a.dx, 0.25);
        assert_eq!(a.max_nx, 1000, "1000 bin の範囲なら nx <= 1000");
    }

    #[test]
    fn no_aliasing_at_coarse_resolution() {
        // 全ゲノム 62861 bin を 4000 画素 -> 約 15.7 bin/画素。
        let s = GridSpec::new(0.0, 62861.0, 5000.0, 4000);
        assert!(s.aliasing().is_none());
        assert_eq!(s.ny, 160);
    }

    #[test]
    fn aggregates_scores_per_pixel() {
        let path = std::env::temp_dir().join("mbhictools_grid.txt");
        // 同じ画素 (px=11, py=1) に落ちる 2 行と、別の画素に落ちる 1 行。
        std::fs::write(
            &path,
            "bin1\tbin2\tscore\n100\t120\t3\n101\t121\t4\n5\t5\t7\n",
        )
        .unwrap();

        let s = GridSpec::new(0.0, 1000.0, 200.0, 100);
        let g = build(&path, s, None).unwrap();

        assert_eq!(g.get(11, 1), 7, "同じ画素の score は足し合わされる");
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
