//! 染色体長テーブル (`chr<TAB>bp`) と、global bin 番号との対応。
//!
//! 接触ファイルの bin は染色体を連結したゲノム全体の通し番号 (1 始まり) であり、
//! bp ではない。各染色体の bin 数は `ceil(bp / bin_size)` で決まる。
//!
//! この対応が分かると 2 つのことができる:
//!   - 軸を bin ではなく bp で目盛る
//!   - 染色体境界を跨ぐ (= inter-chromosomal な) 接触を除外する。
//!     ゲノムを一直線に並べた図では、境界を跨ぐペアの「距離」に意味はない。

use std::path::Path;

/// 1 本の染色体と、それが占める global bin の範囲 [start_bin, end_bin] (1 始まり・両端含む)。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Chrom {
    pub name: String,
    pub length_bp: u64,
    pub bins: u32,
    pub start_bin: u32,
    pub end_bin: u32,
}

/// 染色体テーブル。bin -> 染色体の逆引きができる。
#[derive(Debug, Clone)]
pub struct ChromIndex {
    chroms: Vec<Chrom>,
    /// 各染色体の先頭 bin (昇順)。bin からの逆引きに二分探索で使う。
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
            Error::Io(e) => write!(f, "入出力エラー: {e}"),
            Error::Parse { path, line, reason } => {
                write!(
                    f,
                    "{path}:{line}: {reason} (`染色体名<TAB>長さ(bp)` を期待)"
                )
            }
            Error::Empty(path) => write!(f, "{path}: 染色体が 1 本も読めませんでした"),
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
    /// 染色体長ファイルを読み、`bin_size` bp/bin として bin 範囲を割り当てる。
    pub fn load(path: &Path, bin_size: u32) -> Result<Self, Error> {
        let text = std::fs::read_to_string(path)?;
        Self::parse(&text, bin_size, &path.display().to_string())
    }

    fn parse(text: &str, bin_size: u32, path: &str) -> Result<Self, Error> {
        assert!(bin_size > 0, "bin_size は正であること");
        let mut chroms = Vec::new();
        let mut next_start: u32 = 1; // global bin は 1 始まり

        for (i, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let (name, bp) = line.split_once('\t').ok_or_else(|| Error::Parse {
                path: path.to_string(),
                line: i + 1,
                reason: "TAB 区切りの 2 列がありません".into(),
            })?;
            let length_bp: u64 = bp.trim().parse().map_err(|_| Error::Parse {
                path: path.to_string(),
                line: i + 1,
                reason: format!("長さ `{}` を整数として読めません", bp.trim()),
            })?;
            if length_bp == 0 {
                return Err(Error::Parse {
                    path: path.to_string(),
                    line: i + 1,
                    reason: "長さが 0 です".into(),
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

    /// 染色体長の合計 (bp)。
    pub fn total_bp(path: &Path) -> Result<u64, Error> {
        // bin_size に依存しないので、仮の 1 bp/bin で読んで合計だけ取る。
        let idx = Self::load(path, 1)?;
        Ok(idx.chroms.iter().map(|c| c.length_bp).sum())
    }

    pub fn chroms(&self) -> &[Chrom] {
        &self.chroms
    }

    pub fn bin_size(&self) -> u32 {
        self.bin_size
    }

    /// ゲノム全体の bin 数 (= 最終染色体の end_bin)。
    pub fn total_bins(&self) -> u32 {
        self.chroms.last().map_or(0, |c| c.end_bin)
    }

    /// `bin` を含む染色体の添字。テーブルの範囲外なら `None`。
    pub fn chrom_of(&self, bin: u32) -> Option<usize> {
        if bin < self.starts[0] || bin > self.total_bins() {
            return None;
        }
        // starts は昇順。bin 以下で最大の start を持つ染色体が答え。
        Some(match self.starts.binary_search(&bin) {
            Ok(i) => i,
            Err(i) => i - 1,
        })
    }

    /// 2 つの bin が同一染色体に属するか (どちらかが範囲外なら false)。
    pub fn same_chrom(&self, bin1: u32, bin2: u32) -> bool {
        match (self.chrom_of(bin1), self.chrom_of(bin2)) {
            (Some(a), Some(b)) => a == b,
            _ => false,
        }
    }

    /// 染色体の境目にあたる global bin 座標 (最初の染色体の始点と最後の終点は含まない)。
    /// 図に区切り線を引くのに使う。
    pub fn boundaries(&self) -> Vec<u32> {
        self.chroms.iter().skip(1).map(|c| c.start_bin).collect()
    }
}

/// `bin_size` を `総 bp / 最大 bin` から推定する。
///
/// bin は 1 始まりの通し番号なので、ゲノムを覆う bin 数 ≒ 最大 bin。
/// 端数の丸めぶん誤差が出るため、10 の位で丸めた値を返す (例: 200.05 -> 200)。
pub fn estimate_bin_size(total_bp: u64, max_bin: u32) -> u32 {
    assert!(max_bin > 0, "最大 bin は正であること");
    let raw = total_bp as f64 / max_bin as f64;
    // 1, 10, 100, ... のうち raw の桁に合う刻みで丸めると、実際に使われる
    // キリの良い bin_size (100, 200, 500, 1000, ...) に一致しやすい。
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
        assert_eq!(idx.chrom_of(27896), Some(0), "染色体 I の末端");
        assert_eq!(idx.chrom_of(27897), Some(1), "染色体 II の先頭");
        assert_eq!(idx.chrom_of(50596), Some(1));
        assert_eq!(idx.chrom_of(50597), Some(2));
        assert_eq!(idx.chrom_of(62861), Some(2), "ゲノム末端");
        assert_eq!(idx.chrom_of(0), None, "bin は 1 始まり");
        assert_eq!(idx.chrom_of(62862), None, "範囲外");
    }

    #[test]
    fn detects_boundary_crossing_pairs() {
        let idx = pombe(200);
        assert!(idx.same_chrom(27890, 27896), "どちらも染色体 I");
        assert!(!idx.same_chrom(27896, 27897), "I と II を跨ぐ");
        assert!(!idx.same_chrom(100, 62861), "I と III を跨ぐ");
        assert!(!idx.same_chrom(1, 62862), "範囲外を含む");
    }

    #[test]
    fn reports_boundaries_between_chromosomes() {
        assert_eq!(pombe(200).boundaries(), vec![27897, 50597]);
    }

    #[test]
    fn estimates_pombe_bin_size() {
        // 実データ: 総 12,571,820 bp / 最大 bin 62,843 = 200.05 bp/bin
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
