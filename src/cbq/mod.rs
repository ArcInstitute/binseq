//! # CBQ Format
//!
//! CBQ is the recommended BINSEQ variant: a binary format built around blocked
//! columnar storage, optimized for storage efficiency and parallel processing.
//!
//! ## Overview
//!
//! Records are grouped into blocks; within each block every record attribute
//! (lengths, sequences, quality scores, headers, flags) is stored as a separate
//! column, and each column is ZSTD-compressed independently. Compared to the
//! row-based VBQ this gives better compression ratios, faster reads
//! (pay-per-use decompression), and simpler record parsing.
//!
//! Sequences are two-bit encoded and lossless by default: the positions of
//! ambiguous nucleotides (`N`) are tracked in an Elias-Fano encoded column and
//! backfilled on decode.
//!
//! ## Usage
//!
//! Write CBQ files through [`BinseqWriterBuilder`](crate::BinseqWriterBuilder)
//! (or [`ColumnarBlockWriter`](crate::cbq::ColumnarBlockWriter) directly), and read
//! them with [`MmapReader`](crate::cbq::MmapReader) via the
//! [`ParallelProcessor`](crate::ParallelProcessor) trait — see the crate-level example.
//!
//! ## File Structure
//!
//! A CBQ file consists of a [`FileHeader`](crate::cbq::FileHeader), followed by record
//! blocks and an embedded [`Index`](crate::cbq::Index). Each record block is a
//! [`BlockHeader`](crate::cbq::BlockHeader) with block metadata followed by a
//! [`ColumnarBlock`](crate::cbq::ColumnarBlock) of data. The
//! [`IndexHeader`](crate::cbq::IndexHeader) and [`IndexFooter`](crate::cbq::IndexFooter)
//! locate the index for memory-mapped reading.
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
//! Each column is ZSTD-compressed and stored contiguously; the
//! [`BlockHeader`](crate::cbq::BlockHeader) records the compressed and
//! uncompressed sizes of every column.
//!
//! ```text
//! [BlockHeader][col1][col2][col3]...[BlockHeader][col1][col2][col3]...
//! ```
//!
//! Column order:
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
    IndexHeader, RefRecord, RefRecordIter,
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
