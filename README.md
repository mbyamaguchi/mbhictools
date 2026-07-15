# mbhictools_rs

Rust tools for Hi-C contact maps. Aggregation and rendering live in the
`mbhictools_rs` library; commands live in `src/bin/`.

| Command | Purpose |
| --- | --- |
| `draw_hic` | Draw a triangular contact map as PNG |
| `draw_ps` | Draw the distance curve P(s) as PNG |

## draw_hic

Aggregates sparse contacts (`bin1<TAB>bin2<TAB>score`) into a display grid in one
pass over the file and renders them in rotated coordinates. Memory scales with the
grid, not the row count, so 100M+ rows are fine. Ported from an earlier R
implementation.

```sh
cargo build --release

# Whole genome
./target/release/draw_hic contacts.txt \
    --lengths chrom_lengths.txt --nx 4000 -o hicmap.png

# Zoom into the first chromosome (ranges are in bins, not bp)
./target/release/draw_hic contacts.txt \
    --lengths chrom_lengths.txt \
    --x-start 1 --x-end 27897 --nx 3000 -o chr1.png
```

See `draw_hic --help` for all options.

### Input assumptions

| Assumption | Consequence |
| --- | --- |
| Upper triangle only (`bin1 <= bin2`) | The lower triangle is symmetric and not stored |
| `score` is an integer >= 1 | Scores sum exactly in u64, no rounding error |
| Distance within a limit (e.g. 1 Mbp) | Only a band along the diagonal holds data |

Bins are 1-based indices over the concatenated genome in `bin_size` bp steps, not bp.
Given `--lengths` (`chr<TAB>bp`), `bin_size` is estimated as total bp / max bin; axes
are then labelled in bp, contacts crossing chromosome boundaries are dropped, and the
boundaries are marked. For the example data above, that is 200 bp/bin and a
5000 bin = 1 Mbp limit.

### Why a triangle

With distance capped, a square map is mostly empty by construction. Rotating each
contact by 45 degrees,

```
xr = (bin1 + bin2) / 2    genomic position (midpoint of the two bins)
yr = (bin2 - bin1) / 2    half the interaction distance
```

turns the band into a triangle whose base is the genome and whose height is half the
distance limit, wasting no pixels. The y axis is labelled with the real distance
(`2 * yr`).

### Pixel count

Raising `--nx` too far degrades the figure rather than resolving it. Rotated contacts
land on a checkerboard lattice: `bin1` and `bin2` are integers, so `xr` and `yr` are
multiples of 0.5, but `xr + yr = bin2` is always integral, so only points with
integral `xr + yr` exist. With pixels `dx` bins square:

| `dx` | Lattice points per pixel |
| --- | --- |
| 1 | Exactly 2 in every pixel (uniform) |
| 0.5 | 1 where `i+j` is even, 0 where odd — a full checkerboard |
| >= 1 | `2 * dx^2` on average |

So `dx = 1 bin/pixel` is a hard floor; below it, empty pixels alternate into moiré.
Measured: at `dx = 0.25` the pixel count is 16x higher but filled pixels only rise
from 133k to 190k, so occupancy drops from 66.4% to 5.9%. `draw_hic` detects this,
refuses to render by default, and reports the largest `--nx` usable for the range.
`--allow-aliasing` overrides.

Pixels are always square (`dy == dx`), so `ny` follows from `--nx` and the range.
Height is deliberately not settable on its own, as that would let the aspect ratio
misrepresent distance.

## draw_ps

Computes the contact frequency versus genomic distance curve,
`P(s) = contacts(s) / pairs(s)`, over intra-chromosomal pairs, and plots it log-log.

```sh
./target/release/draw_ps contacts.txt \
    --lengths chrom_lengths.txt -o distance_curve.png --dump curve.tsv

# Only the range unaffected by the pipeline's short-range filter (see below)
./target/release/draw_ps contacts.txt \
    --lengths chrom_lengths.txt --min-s 51 -o clean.png
```

`--lengths` is required here: the denominator comes from it.

### Log binning

Per-separation points crowd the right of a log axis and are noisiest where pairs are
scarcest, so separations are pooled into geometrically spaced ranges
(`--bins-per-decade`, default 10; `--no-logbin` to disable).

Within a range the numerator and denominator are summed separately and divided once.
Averaging the per-separation P(s) would weight a separation holding a handful of
pairs the same as one holding thousands. Each point plots at the pair-weighted
geometric mean of its range, the centre of mass of its denominator on a log axis.

Gaps are not interpolated. A range with pairs but no contacts is a real observation
of zero, so it stays 0 in the TSV and is dropped from a log axis with a count on
stderr, rather than being filled in with a fabricated value.

### A known artifact in this data

Both input files step by a factor of 1.92 at exactly s = 10 kb: contacts jump while
pairs stay flat. Two independent datasets sharing the same threshold and the same
near-2x factor point at the upstream pipeline, which appears to discard two of the
four read orientations below 10 kb. P(s) below that is depressed accordingly. Pass
`--min-s 51` to plot only the unaffected range, where the slope is -1.22, within the
usual -1 to -1.5.

## Layout

| Module | Role |
| --- | --- |
| `contact` | Parallel parsing (mmap + rayon, split on line boundaries) |
| `chrom` | Chromosome lengths and their global bin ranges |
| `grid` | Aggregation in rotated coordinates; pixel geometry |
| `curve` | Contacts and pairs per distance; log binning |
| `render` | Value transform (log10, quantile clip), palette, PNG |
| `font` | Picking a usable font for labels |

Aggregation uses u64 atomic adds into one shared grid. Integer scores mean no
per-thread grid copies are needed and the result is exact: the summed `raw` of a
dumped grid matches the summed `score` of the input.

Passing `sans-serif` to plotters defers to fontconfig, which on some systems resolves
to a font reporting zero Latin advance widths (e.g. Droid Sans Fallback); tick labels
then vanish with no error. `font` measures text to pick a sane candidate, and falls
back to English axis labels where no CJK font exists. `--font` overrides.

