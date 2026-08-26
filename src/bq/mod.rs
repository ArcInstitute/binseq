//! # BQ Format
//!
//! BQ (`*.bq`) is the BINSEQ variant for **fixed-length** records without
//! quality scores. Its uniform record size gives constant-time random access
//! with minimal overhead. For variable-length records and optional quality
//! scores use [`cbq`](crate::cbq).
//!
//! For detailed information on the file format, see our [paper](https://www.biorxiv.org/content/10.1101/2025.04.08.647863v1).
//!
//! ## Usage
//!
//! ### Reading
//!
//! ```rust
//! use binseq::{bq, BinseqRecord};
//!
//! let reader = bq::MmapReader::new("./data/subset.bq").unwrap();
//! let num_records = reader.num_records();
//!
//! // Random access to any record in the file
//! let record = reader.get(num_records / 2).unwrap();
//!
//! // Decode the 2-bit encoded sequence back to bytes
//! let mut sbuf = Vec::new();
//! record.decode_s(&mut sbuf).unwrap();
//! ```
//!
//! ### Writing
//!
//! ```rust
//! use binseq::{bq, SequencingRecordBuilder};
//! use std::io::Cursor;
//!
//! // BQ is fixed-length: sequence lengths are set in the header
//! let header = bq::FileHeaderBuilder::new().slen(64).build().unwrap();
//! let mut writer = bq::WriterBuilder::default()
//!     .header(header)
//!     .build(Cursor::new(Vec::new()))
//!     .unwrap();
//!
//! let seq = [b'A'; 64];
//! let record = SequencingRecordBuilder::default()
//!     .s_seq(&seq)
//!     .build()
//!     .unwrap();
//! writer.push(record).unwrap();
//! writer.flush().unwrap();
//! ```
//!
//! Paired records work the same way: set `xlen` on the header and `x_seq` on
//! the record. For streaming over arbitrary readers/writers (e.g. network
//! sockets) see [`StreamReader`](crate::bq::StreamReader) and the
//! `network_streaming` example.
//!
//! ## File Format
//!
//! A BQ file is a fixed-size header followed by densely packed records.
//!
//! ### Header (32 bytes)
//!
//! | Offset | Size (bytes) | Name     | Description                  | Type   |
//! | ------ | ------------ | -------- | ---------------------------- | ------ |
//! | 0      | 4            | magic    | Magic number (0x42534551)    | uint32 |
//! | 4      | 1            | format   | Format version (currently 2) | uint8  |
//! | 5      | 4            | slen     | Sequence length (primary)    | uint32 |
//! | 9      | 4            | xlen     | Sequence length (secondary)  | uint32 |
//! | 13     | 19           | reserved | Reserved for future use      | bytes  |
//!
//! ### Records
//!
//! Each record is a flag field (8 bytes, uint64, implementation-defined) followed
//! by the encoded sequence data (`ceil(N/32) * 8` bytes for sequence length `N`).
//! The leading flag enables filtering without touching sequence data, and the
//! uniform record size makes random access a constant-time offset calculation.
//!
//! ### Encoding
//!
//! Nucleotides are 2-bit encoded (A=00, C=01, G=10, T=11) into little-endian
//! u64 words of 32 bases each, zero-padded in the final word. Non-ACGT characters
//! are unsupported; see [`Policy`](crate::Policy) for handling options and
//! [`bitnuc`] for implementation details.

mod header;
mod reader;
mod writer;

pub use header::{FILE_MAGIC, FileHeader, FileHeaderBuilder, SIZE_HEADER};
pub use reader::{MmapReader, RefRecord, StreamReader};
pub use writer::{Writer, WriterBuilder};
