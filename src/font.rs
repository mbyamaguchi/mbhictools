//! Picking a font for labels.
//!
//! Passing `sans-serif` to plotters defers to fontconfig, which on some systems
//! resolves to a font reporting zero advance widths (e.g. Droid Sans Fallback).
//! The font then loads without error, but every tick label collapses to zero width
//! and silently disappears. So instead of trusting a family name, measure it.

use plotters::prelude::*;

/// Families to try, in order.
const CANDIDATES: &[&str] = &[
    "DejaVu Sans",
    "Liberation Sans",
    "Arial",
    "Helvetica",
    "Noto Sans",
    "sans-serif",
];

/// Text used to measure a font, and the size it is measured at.
const PROBE: &str = "0123456789";
const PROBE_SIZE: f64 = 20.0;
/// An advance narrower than this per character means the metrics are broken.
const MIN_ADVANCE_PER_CHAR: f64 = PROBE_SIZE / 4.0;

/// Pick a font family whose metrics work.
///
/// `preferred` is measured first and falls back to the candidates if it fails.
/// If nothing passes, return `sans-serif` as a last resort: labels may collapse,
/// but rendering still completes.
pub fn pick(preferred: Option<&str>) -> String {
    if let Some(family) = preferred {
        if metrics_ok(family) {
            return family.to_string();
        }
        eprintln!("warning: font `{family}` reports broken metrics; falling back");
    }

    for &family in CANDIDATES {
        if metrics_ok(family) {
            return family.to_string();
        }
    }

    eprintln!("warning: no font with usable metrics found; labels may be unreadable");
    "sans-serif".into()
}

/// Does this family report sane advance widths?
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

    /// Skip when no font is installed at all: the checks below cannot hold there.
    fn fonts_available() -> bool {
        if metrics_ok("DejaVu Sans") {
            return true;
        }
        eprintln!("skipped: DejaVu Sans is not installed");
        false
    }

    /// Why this module exists: Droid Sans Fallback reports ~zero Latin advances.
    #[test]
    fn rejects_fonts_with_broken_metrics() {
        if !fonts_available() {
            return;
        }
        assert!(!metrics_ok("Droid Sans Fallback"));
    }

    #[test]
    fn picks_a_font_with_usable_metrics() {
        let f = pick(None);
        assert!(metrics_ok(&f) || f == "sans-serif");
    }

    #[test]
    fn falls_back_when_the_preferred_font_is_broken() {
        if !fonts_available() {
            return;
        }
        assert_ne!(pick(Some("Droid Sans Fallback")), "Droid Sans Fallback");
    }

    #[test]
    fn honours_a_working_preferred_font() {
        if !fonts_available() {
            return;
        }
        assert_eq!(pick(Some("DejaVu Sans")), "DejaVu Sans");
    }
}
