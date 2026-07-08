//! `zstd-c` backend: a direct re-export of the `zstd`/`zstd-safe` crates.

pub(crate) use zstd::stream;
pub(crate) use zstd::zstd_safe;
pub(crate) use zstd::{Decoder, Encoder};
