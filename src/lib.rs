#![cfg_attr(docsrs, feature(doc_cfg))]
#![cfg_attr(docsrs, doc(auto_cfg))]
#![doc = include_str!("../README.md")]
//!
//! # Library
//!
//! The `binseq` library provides efficient APIs for working with the [BINSEQ](https://www.biorxiv.org/content/10.1101/2025.04.08.647863v2) file format family:
//!
//! - **[`cbq`]** — the recommended variant: columnar, block-compressed, lossless by default
//! - **[`bq`]** and **[`vbq`]** — earlier variants, still fully supported
//! - [`BinseqRecord`], [`BinseqReader`], and [`BinseqWriter`] — variant-agnostic record, reader, and writer abstractions
//! - [`ParallelProcessor`] — parallel processing of records with arbitrary per-record logic
//! - [`Policy`] — configurable handling of invalid nucleotides when encoding BQ/VBQ (CBQ stores `N` natively)
//!
//! # Example: Parallel Processing
//!
//! ```
//! use binseq::Result;
//! use binseq::prelude::*;
//!
//! #[derive(Clone, Default)]
//! pub struct Processor {
//!     // per-thread state
//! }
//!
//! impl ParallelProcessor for Processor {
//!     fn process_record<B: BinseqRecord>(&mut self, record: B) -> Result<()> {
//!         // per-record logic
//!         Ok(())
//!     }
//!
//!     fn on_batch_complete(&mut self) -> Result<()> {
//!         // per-batch logic (e.g. flush to a shared writer)
//!         Ok(())
//!     }
//! }
//!
//! fn main() -> Result<()> {
//!     // any BINSEQ variant (*.cbq, *.bq, *.vbq); format is sniffed from content
//!     let path = "./data/subset.cbq";
//!
//!     let reader = BinseqReader::new(path)?;
//!     let processor = Processor::default();
//!
//!     // process records in parallel with 8 threads
//!     reader.process_parallel(processor, 8)?;
//!     Ok(())
//! }
//! ```

#![allow(clippy::module_inception)]

/// CBQ - columnar variable-length records, optional quality scores and headers
pub mod cbq;

/// BQ - fixed-length records, no quality scores
pub mod bq;

/// VBQ - variable-length records, optional quality scores, compressed blocks
pub mod vbq;

/// Shared nucleotide encoder for BQ/VBQ writers
mod encoder;

/// Error definitions
pub mod error;

/// Parallel processing
mod parallel;

/// Invalid nucleotide policy
mod policy;

/// Record types and traits shared between BINSEQ variants
mod record;

/// Commonly used types and traits
pub mod prelude;

/// Write operations generic over the BINSEQ variant
pub mod write;

/// Utilities for working with BINSEQ files
pub mod utils;

pub use error::{Error, IntoBinseqError, Result};
pub use parallel::{BinseqReader, ParallelProcessor, ParallelReader};
pub use policy::{Policy, RNG_SEED};
pub use record::{BinseqRecord, SequencingRecord, SequencingRecordBuilder};
pub use write::{BinseqWriter, BinseqWriterBuilder};

/// Re-export `bitnuc::BitSize`
pub use bitnuc_deprec::BitSize;

/// Default quality score for BINSEQ readers without quality scores
pub(crate) const DEFAULT_QUALITY_SCORE: u8 = b'?';
