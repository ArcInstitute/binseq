#![doc = include_str!("../README.md")]
//!
//! # BINSEQ
//!
//! The `binseq` library provides efficient APIs for working with the [BINSEQ](https://www.biorxiv.org/content/10.1101/2025.04.08.647863v2) file format family.
//!
//! It offers methods to read and write BINSEQ files, providing:
//!
//! - Compact multi-bit encoding and decoding of nucleotide sequences through [`bitnuc`](https://docs.rs/bitnuc/latest/bitnuc/)
//! - Support for both single and paired-end sequences
//! - Abstract [`BinseqRecord`] trait for representing records from all variants
//! - Abstract [`BinseqReader`] enum for processing records from all variants
//! - Abstract [`BinseqWriter`] enum for writing records to all variants
//! - Parallel processing capabilities for arbitrary tasks through the [`ParallelProcessor`] trait.
//! - Configurable [`Policy`] for handling invalid nucleotides (BQ/VBQ, CBQ natively supports `N` nucleotides)
//!
//! ## Recent additions (v0.9.0):
//!
//! ### New variant: CBQ
//! **[`cbq`]** is a new variant of BINSEQ that solves many of the pain points around VBQ.
//! The CBQ format is a columnar-block-based format that offers improved compression and faster processing speeds compared to VBQ.
//! It natively supports `N` nucleotides and avoids the need for additional 4-bit encoding.
//!
//! ### Improved interface for writing records
//! **[`BinseqWriter`]** provides a unified interface for writing records generically to BINSEQ files.
//! This makes use of the new [`SequencingRecord`] which provides a cleaner builder API for writing records to BINSEQ files.
//!
//! ## Recent VBQ Format Changes (v0.7.0+)
//!
//! The VBQ format has undergone significant improvements:
//!
//! - **Embedded Index**: VBQ files now contain their index data embedded at the end of the file,
//!   eliminating separate `.vqi` index files and improving portability.
//! - **Headers Support**: Optional sequence identifiers/headers can be stored with each record.
//! - **Extended Capacity**: u64 indexing supports files with more than 4 billion records.
//! - **Multi-bit Encoding**: Support for both 2-bit and 4-bit nucleotide encodings.
//!
//! Legacy VBQ files are automatically migrated to the new format when accessed.
//!
//! # Example: Memory-mapped Access
//!
//! ```
//! use binseq::Result;
//! use binseq::prelude::*;
//!
//! #[derive(Clone, Default)]
//! pub struct Processor {
//!     // Define fields here
//! }
//!
//! impl ParallelProcessor for Processor {
//!     fn process_record<B: BinseqRecord>(&mut self, record: B) -> Result<()> {
//!         // Implement per-record logic here
//!         Ok(())
//!     }
//!
//!     fn on_batch_complete(&mut self) -> Result<()> {
//!         // Implement per-batch logic here
//!         Ok(())
//!     }
//! }
//!
//! fn main() -> Result<()> {
//!     // provide an input path (*.bq or *.vbq)
//!     let path = "./data/subset.bq";
//!
//!     // open a reader
//!     let reader = BinseqReader::new(path)?;
//!
//!     // initialize a processor
//!     let processor = Processor::default();
//!
//!     // process the records in parallel with 8 threads
//!     reader.process_parallel(processor, 8)?;
//!     Ok(())
//! }
//! ```

#![allow(clippy::module_inception)]

/// BQ - fixed length records, no quality scores
pub mod bq;

/// Error definitions
pub mod error;

/// Parallel processing
mod parallel;

/// Invalid nucleotide policy
mod policy;

/// Record types and traits shared between BINSEQ variants
mod record;

/// VBQ - Variable length records, optional quality scores, compressed blocks
pub mod vbq;

/// CBQ - Columnar variable length records, optional quality scores and headers
pub mod cbq;

/// Prelude - Commonly used types and traits
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
pub use bitnuc::BitSize;

/// Default quality score for BINSEQ readers without quality scores
pub(crate) const DEFAULT_QUALITY_SCORE: u8 = b'?';
