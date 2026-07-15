//! Hi-C contact map visualisation.
//!
//! Aggregates sparse contacts (`bin1<TAB>bin2<TAB>score`) into a display grid in one
//! pass and draws them as a triangular map in rotated coordinates.
//!
//! - [`contact`] parallel parsing of the input file
//! - [`chrom`]   chromosome lengths and their global bin ranges
//! - [`grid`]    aggregation in rotated coordinates; pixel geometry
//! - [`curve`]   the contact-frequency versus distance curve, P(s)
//! - [`render`]  value transform, palette and PNG output
//! - [`font`]    picking a font whose metrics actually work

pub mod chrom;
pub mod contact;
pub mod curve;
pub mod font;
pub mod grid;
pub mod render;
