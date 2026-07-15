//! スパース接触ファイル (`bin1<TAB>bin2<TAB>score`) の並列読み込み。
//!
//! 入力の前提 (本ツールが対象とするファイルの仕様):
//!   - TSV, 1 行目はヘッダ (`bin1  bin2  score`)
//!   - 上三角成分のみ (`bin1 <= bin2`)
//!   - `score` は 1 以上の整数
//!   - bin 間距離 `bin2 - bin1` は一定値以下 (例: 1 Mbp 相当)
//!
//! score が整数であることから、集計は f64 ではなく u64 で厳密に行える。
//! これは丸め誤差がないだけでなく、アトミック加算による共有グリッドへの
//! 並列集計 (=スレッドごとのグリッド複製が不要) を可能にする。
//!
//! ファイルは mmap して行境界で分割し、rayon で並列にパースする。
//! 全レコードを保持することはなく、メモリはグリッドぶんのみ。

use std::fs::File;
use std::path::Path;

use memmap2::Mmap;
use rayon::prelude::*;

/// 並列パース時に 1 タスクへ渡すおおよそのバイト数。行境界へ切り上げられる。
const CHUNK_BYTES: usize = 8 << 20;

/// 1 パスの走査結果。`visit` へ渡した closure の呼ばれた回数と、壊れた行の数。
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ScanStats {
    /// 正常にパースできた行数。
    pub rows: u64,
    /// パースできず読み飛ばした行数。
    pub malformed: u64,
}

/// 読み込み時に起こりうるエラー。
#[derive(Debug)]
pub enum Error {
    Io(std::io::Error),
    /// 壊れた行が多すぎる (ファイル形式の取り違えが疑われる)。
    TooManyMalformed {
        path: String,
        malformed: u64,
        rows: u64,
    },
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Io(e) => write!(f, "入出力エラー: {e}"),
            Error::TooManyMalformed {
                path,
                malformed,
                rows,
            } => write!(
                f,
                "{path}: パースできない行が多すぎます ({malformed} / {} 行)。\
                 `bin1<TAB>bin2<TAB>score` 形式の TSV か確認してください",
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

/// ファイル全体を 1 パス走査し、各レコードで `visit(bin1, bin2, score)` を呼ぶ。
///
/// `visit` は複数スレッドから同時に呼ばれる。呼び出し順は入力順とは限らない。
/// 集計はアトミックか、可換な操作 (加算・min・max) に限ること。
pub fn visit<F>(path: &Path, visit: F) -> Result<ScanStats, Error>
where
    F: Fn(u32, u32, u32) + Sync,
{
    let file = File::open(path)?;
    // SAFETY: 読み込み専用に mmap する。走査中に外部からファイルが切り詰められると
    // 未定義動作になりうるが、解析対象の入力ファイルは不変とみなす。
    let mmap = unsafe { Mmap::map(&file)? };

    let body = skip_header(&mmap);
    let stats = chunk_ranges(body)
        .into_par_iter()
        .map(|(start, end)| visit_chunk(&body[start..end], &visit))
        .reduce(ScanStats::default, |a, b| ScanStats {
            rows: a.rows + b.rows,
            malformed: a.malformed + b.malformed,
        });

    // 数行の破損は許容するが、大半が読めないなら形式の取り違えとして落とす。
    if stats.malformed > 0 && stats.malformed * 2 > stats.rows {
        return Err(Error::TooManyMalformed {
            path: path.display().to_string(),
            malformed: stats.malformed,
            rows: stats.rows,
        });
    }
    Ok(stats)
}

/// 先頭行が数字で始まらなければヘッダとみなして読み飛ばす。
fn skip_header(buf: &[u8]) -> &[u8] {
    match buf.first() {
        Some(b) if b.is_ascii_digit() => buf,
        _ => match memchr_newline(buf, 0) {
            Some(nl) => &buf[nl + 1..],
            None => &buf[buf.len()..],
        },
    }
}

/// `buf` を行境界に揃った範囲へ分割する。各範囲は独立にパースできる。
fn chunk_ranges(buf: &[u8]) -> Vec<(usize, usize)> {
    let mut ranges = Vec::with_capacity(buf.len() / CHUNK_BYTES + 1);
    let mut start = 0;
    while start < buf.len() {
        // 目標位置から次の改行まで進め、行を跨がないようにする。
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

/// `bin1<TAB>bin2<TAB>score` を 1 行パースする。3 列目より後ろは無視する。
fn parse_line(line: &[u8]) -> Option<(u32, u32, u32)> {
    let mut fields = line.split(|&b| b == b'\t');
    let b1 = parse_u32(fields.next()?)?;
    let b2 = parse_u32(fields.next()?)?;
    let score = parse_u32(fields.next()?)?;
    Some((b1, b2, score))
}

/// 10 進の符号なし整数をパースする。空文字列・非数字・桁溢れは `None`。
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

/// `from` 以降で最初の `\n` の位置。
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
        assert_eq!(parse_line(b"27\t79"), None, "列が足りない");
        assert_eq!(parse_line(b"27\t79\t1.5"), None, "score が整数でない");
        assert_eq!(parse_line(b"27\t79\t"), None, "score が空");
        assert_eq!(parse_line(b"a\tb\tc"), None, "数字でない");
        assert_eq!(parse_line(b"27\t79\t-1"), None, "負の score");
    }

    #[test]
    fn parse_u32_rejects_overflow() {
        assert_eq!(parse_u32(b"4294967295"), Some(u32::MAX));
        assert_eq!(parse_u32(b"4294967296"), None);
    }

    #[test]
    fn visits_every_record_exactly_once() {
        // ヘッダ + 末尾改行あり。score の総和で全件訪問を確認する。
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
        assert_eq!(stats.rows, 1, "数字で始まる先頭行はデータとして読む");
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
