//! # CBQ Format
//!
//! CBQ is a high-performance binary format built around blocked columnar storage.
//! It optimizes for storage efficiency and parallel processing of records.
//!
//! ## Overview
//!
//! CBQ was built to solve the rough edges of VBQ.
//! It keeps the blocked structure of VBQ, but instead of interleaving the internal data of all records in the block, it stores each attribute in a separate column.
//! Each of these columns are then ZSTD compressed and optionally decoded when reading.
//!
//! This has a few benefits and advantages over VBQ:
//!
//! 1. Better compression ratios for each individual attribute.
//! 2. Significantly faster throughput for reading (easier decompression + pay-per-use decompression).
//! 3. Simple record parsing and manipulation.
//!
//! Notably this format *only* performs two-bit encoding of sequences.
//! However, it tracks the positions of all ambiguous nucleotides (`N`) within the sequence.
//! When it is decoded and the two-bit encoded sequence is decoded back to nucleotides, the `N` positions are backfilled with `N`.
//!
//! To make use of the sparse-but-clustered nature of the `N`-positions, we make use of an Elias-Fano encoding of the `N`-positions.
//! This encoding is then used to efficiently store and retrieve the positions of `N`s within the sequence.
//!
//! ## File Structure
//!
//! A CBQ file consists of a [`FileHeader`](cbq::FileHeader), followed by record blocks and an embedded [`Index`](cbq::Index).
//! Each record block is composed of a [`BlockHeader`](cbq::BlockHeader) which provides metadata about the block, and a [`ColumnarBlock`](cbq::ColumnarBlock) containing the actual data.
//!
//! The [`IndexHeader`](cbq::IndexHeader) and [`IndexFooter`](cbq::IndexFooter) are used to locate and access the data within the file when reading as memory mapped.
//!
//! ```text
//! ┌───────────────────┐
//! │    File Header    │ 64 bytes
//! ├───────────────────┤
//! │   Block Header    │ 96 bytes
//! ├───────────────────┤
//! │                   │
//! │   Block Records   │ Variable size
//! │                   │
//! ├───────────────────┤
//! │       ...         │ More blocks
//! ├───────────────────┤
//! │    Index Header   │ 24 bytes
//! ├───────────────────┤
//! │ Compressed Index  │ Variable size
//! ├───────────────────┤
//! │    Index Footer   │ 16 bytes
//! └───────────────────┘
//! ```
//!
//! ## Block Format
//!
//! The blocks on-disk are stored as ZSTD compressed data.
//! Each column is ZSTD compressed and stored contiguously next to each other.
//!
//! The [BlockHeader](cbq::BlockHeader) contains the compressed sizes of each of the columns as well as the relevant information for their uncompressed sizes.
//!
//! ```text
//! [BlockHeader][col1][col2][col3]...[BlockHeader][col1][col2][col3]...
//! ```
//!
//! The order of columns in the block is as follows:
//!
//! 1. `z_seq_len` - sequence lengths
//! 2. `z_header_len` - header lengths (optional)
//! 3. `z_npos` - Elias-Fano encoded positions of N's (optional)
//! 4. `z_seq` - sequence data (2-bit encoded)
//! 5. `z_flags` - flags (optional)
//! 6. `z_headers` - sequence headers (optional)
//! 7. `z_qual` - sequence quality scores (optional)

mod core;
mod read;
mod write;

pub use core::{
    BlockHeader, BlockRange, ColumnarBlock, FileHeader, FileHeaderBuilder, Index, IndexFooter,
    IndexHeader, RefRecord, RefRecordIter, SequencingRecord, SequencingRecordBuilder,
};
pub use read::{MmapReader, Reader};
pub use write::ColumnarBlockWriter;

/// The magic number for CBQ files.
pub const FILE_MAGIC: &[u8; 7] = b"CBQFILE";

/// The magic number for CBQ blocks.
pub const BLOCK_MAGIC: &[u8; 3] = b"BLK";

/// The magic number for CBQ index files.
pub const INDEX_MAGIC: &[u8; 8] = b"CBQINDEX";

/// The current file version.
pub const FILE_VERSION: u8 = 1;

/// The default block size.
pub const DEFAULT_BLOCK_SIZE: u64 = 1024 * 1024;

/// The default compression level.
pub const DEFAULT_COMPRESSION_LEVEL: u64 = 0;
