//! Hi-C contact map の可視化。
//!
//! スパースな接触データ (`bin1<TAB>bin2<TAB>score`) を 1 パスで表示グリッドへ
//! 集計し、回転座標の三角形 contact map として PNG に描く。
//!
//! - [`contact`] 入力ファイルの並列パース
//! - [`chrom`]   染色体長テーブルと global bin の対応
//! - [`grid`]    回転座標への集計と画素数の決定
//! - [`render`]  値変換・配色と PNG 出力
//! - [`font`]    ラベル描画に使えるフォントの選択

pub mod chrom;
pub mod contact;
pub mod font;
pub mod grid;
pub mod render;
