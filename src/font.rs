//! 描画に使うフォントの選択。
//!
//! # なぜ必要か
//!
//! plotters に `sans-serif` を指定すると fontconfig の代替解決に委ねられる。
//! ところが環境によっては、これが Latin の字送り幅を 0 として返すフォント
//! (例: Droid Sans Fallback) に解決されてしまう。するとフォントの読み込み自体は
//! 成功したまま、目盛りの数値だけが幅 0 に潰れて消える。エラーは出ない。
//!
//! そこで候補を順に当たり、実際に文字幅を測って正気なものを選ぶ。
//! 併せて、選ばれたフォントが日本語を持つかどうかも持ち回る。持たない環境で
//! 日本語ラベルを出すと豆腐になるので、その場合は英語ラベルへ切り替える。

use plotters::prelude::*;

/// 選ばれたフォント。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Font {
    /// plotters へ渡すフォントファミリ名。
    pub family: String,
    /// 日本語 (CJK) の字形を持つか。ラベルの言語を決めるのに使う。
    pub cjk: bool,
}

/// 候補と、それが CJK を持つか。上から順に試す。
///
/// CJK の有無は名前で決め打ちする。字形を実際に持つかは文字幅からは判定できず
/// (字形の無いフォントは豆腐を返し、豆腐にも幅があるため)、名前で見るしかない。
const CANDIDATES: &[(&str, bool)] = &[
    // 日本語が出せるもの
    ("Noto Sans CJK JP", true),
    ("Source Han Sans JP", true),
    ("IPAPGothic", true),
    ("IPAGothic", true),
    ("TakaoPGothic", true),
    ("VL PGothic", true),
    ("Yu Gothic", true),
    ("Meiryo", true),
    ("Hiragino Sans", true),
    ("MS Gothic", true),
    // Latin のみ。日本語は出せないが、数値ラベルは正しく出る
    ("DejaVu Sans", false),
    ("Liberation Sans", false),
    ("Arial", false),
    ("Helvetica", false),
    ("sans-serif", false),
];

/// 幅の検査に使う文字列と、その想定される最小幅の目安。
const PROBE: &str = "0123456789";
const PROBE_SIZE: f64 = 20.0;
/// 1 文字あたりこれを下回る字送りは、フォントの計量が壊れているとみなす。
const MIN_ADVANCE_PER_CHAR: f64 = PROBE_SIZE / 4.0;

/// 使えるフォントを選ぶ。
///
/// `preferred` を指定した場合はそれだけを検査し、駄目なら候補へ落ちる。
/// どれも通らなければ最後の手段として `sans-serif` を返す (数値が潰れうるが、
/// 描画そのものは続行できる)。
pub fn pick(preferred: Option<&str>) -> Font {
    if let Some(family) = preferred {
        if metrics_ok(family) {
            // 指定されたフォントの CJK 有無は、候補表に載っていれば従い、
            // 未知なら「持つ」とみなす (利用者が意図して選んだものを尊重する)。
            let cjk = CANDIDATES
                .iter()
                .find(|(name, _)| name.eq_ignore_ascii_case(family))
                .is_none_or(|(_, cjk)| *cjk);
            return Font {
                family: family.to_string(),
                cjk,
            };
        }
        eprintln!(
            "警告: フォント `{family}` は文字幅を正しく返しません。既定の候補から選び直します"
        );
    }

    for &(family, cjk) in CANDIDATES {
        if metrics_ok(family) {
            return Font {
                family: family.to_string(),
                cjk,
            };
        }
    }

    eprintln!(
        "警告: 文字幅を正しく返すフォントが見つかりませんでした。ラベルが潰れる可能性があります"
    );
    Font {
        family: "sans-serif".into(),
        cjk: false,
    }
}

/// そのフォントが Latin の字送りをまともに返すか。
///
/// 壊れたフォントは幅 0 付近を返す。1 文字あたりの字送りが極端に狭くないことを見る。
fn metrics_ok(family: &str) -> bool {
    let font: FontDesc = (family, PROBE_SIZE).into_font();
    match font.layout_box(PROBE) {
        Ok(((x0, _), (x1, _))) => {
            let width = (x1 - x0) as f64;
            width >= MIN_ADVANCE_PER_CHAR * PROBE.chars().count() as f64
        }
        Err(_) => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// このモジュールの存在理由: 壊れたフォントを弾けること。
    /// Droid Sans Fallback は Latin の字送りを 0 付近で返すので落ちるはず。
    #[test]
    fn rejects_fonts_with_broken_latin_metrics() {
        if !metrics_ok("DejaVu Sans") {
            // フォントの入っていない環境ではこの検査自体が成立しない。
            eprintln!("DejaVu Sans が無いためスキップ");
            return;
        }
        assert!(
            !metrics_ok("Droid Sans Fallback"),
            "Latin の字送りが潰れるフォントは採用してはいけない"
        );
    }

    #[test]
    fn picks_a_font_with_usable_metrics() {
        let f = pick(None);
        assert!(metrics_ok(&f.family) || f.family == "sans-serif");
    }

    #[test]
    fn falls_back_when_the_preferred_font_is_broken() {
        if !metrics_ok("DejaVu Sans") {
            eprintln!("フォントが無いためスキップ");
            return;
        }
        let f = pick(Some("Droid Sans Fallback"));
        assert_ne!(f.family, "Droid Sans Fallback", "壊れた指定は採用しない");
    }

    #[test]
    fn honours_a_working_preferred_font() {
        if !metrics_ok("DejaVu Sans") {
            eprintln!("フォントが無いためスキップ");
            return;
        }
        let f = pick(Some("DejaVu Sans"));
        assert_eq!(
            f,
            Font {
                family: "DejaVu Sans".into(),
                cjk: false
            }
        );
    }
}
