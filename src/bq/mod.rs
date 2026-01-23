//! # bq
//!
//! *.bq files are BINSEQ variants for **fixed-length** records and **does not support quality scores**.
//!
//! For variable-length records and optional quality scores use the [`vbq`](crate::vbq) module.
//!
//! This module contains the utilities for reading, writing, and interacting with BINSEQ files.
//!
//! For detailed information on the file format, see our [paper](https://www.biorxiv.org/content/10.1101/2025.04.08.647863v1).
//!
//! ## Usage
//!
//! ### Reading
//! ```rust
//! use binseq::{bq, BinseqRecord};
//! use rand::{thread_rng, Rng};
//!
//! let path = "./data/subset.bq";
//! let reader = bq::MmapReader::new(path).unwrap();
//!
//! // We can easily determine the number of records in the file
//! let num_records = reader.num_records();
//!
//! // We have random access to any record within the range
//! let random_index = thread_rng().gen_range(0..num_records);
//! let record = reader.get(random_index).unwrap();
//!
//! // We can easily decode the (2bit)encoded sequence back to a sequence of bytes
//! let mut sbuf = Vec::new();
//! let mut xbuf = Vec::new();
//!
//! record.decode_s(&mut sbuf);
//! if record.is_paired() {
//!     record.decode_x(&mut xbuf);
//! }
//! ```
//!
//! ### Writing
//!
//! #### Writing unpaired sequences
//!
//! ```rust
//! use binseq::{bq, SequencingRecordBuilder};
//! use std::io::Cursor;
//!
//! // Create an in-memory buffer for output
//! let output_handle = Cursor::new(Vec::new());
//!
//! // Initialize our BINSEQ header (64 bp, only primary)
//! let header = bq::FileHeaderBuilder::new().slen(64).build().unwrap();
//!
//! // Initialize our BINSEQ writer
//! let mut writer = bq::WriterBuilder::default()
//!     .header(header)
//!     .build(output_handle)
//!     .unwrap();
//!
//! // Generate a random sequence
//! let seq = [b'A'; 64];
//!
//! // Build a record and write it to the file
//! let record = SequencingRecordBuilder::default()
//!     .s_seq(&seq)
//!     .flag(0)
//!     .build()
//!     .unwrap();
//! writer.push(record).unwrap();
//!
//! // Flush the writer
//! writer.flush().unwrap();
//! ```
//!
//! #### Writing paired sequences
//!
//! ```rust
//! use binseq::{bq, SequencingRecordBuilder};
//! use std::io::Cursor;
//!
//! // Create an in-memory buffer for output
//! let output_handle = Cursor::new(Vec::new());
//!
//! // Initialize our BINSEQ header (64 bp and 128bp)
//! let header = bq::FileHeaderBuilder::new().slen(64).xlen(128).build().unwrap();
//!
//! // Initialize our BINSEQ writer
//! let mut writer = bq::WriterBuilder::default()
//!     .header(header)
//!     .build(output_handle)
//!     .unwrap();
//!
//! // Generate paired sequences
//! let primary = [b'A'; 64];
//! let secondary = [b'C'; 128];
//!
//! // Build a paired record and write it to the file
//! let record = SequencingRecordBuilder::default()
//!     .s_seq(&primary)
//!     .x_seq(&secondary)
//!     .flag(0)
//!     .build()
//!     .unwrap();
//! writer.push(record).unwrap();
//!
//! // Flush the writer
//! writer.flush().unwrap();
//! ```
//!
//! # Example: Streaming Access
//!
//! ```
//! use binseq::{Policy, Result, BinseqRecord, SequencingRecordBuilder};
//! use binseq::bq::{FileHeaderBuilder, StreamReader, StreamWriterBuilder};
//! use std::io::{BufReader, Cursor};
//!
//! fn main() -> Result<()> {
//!     // Create a header for sequences of length 100
//!     let header = FileHeaderBuilder::new().slen(100).build()?;
//!
//!     // Create a stream writer
//!     let mut writer = StreamWriterBuilder::default()
//!         .header(header)
//!         .buffer_capacity(8192)
//!         .build(Cursor::new(Vec::new()))?;
//!
//!     // Write sequences
//!     let sequence = b"ACGT".repeat(25); // 100 nucleotides
//!     let record = SequencingRecordBuilder::default()
//!         .s_seq(&sequence)
//!         .flag(0)
//!         .build()?;
//!     writer.push(record)?;
//!
//!     // Get the inner buffer
//!     let buffer = writer.into_inner()?;
//!     let data = buffer.into_inner();
//!
//!     // Create a stream reader
//!     let mut reader = StreamReader::new(BufReader::new(Cursor::new(data)));
//!
//!     // Process records as they arrive
//!     while let Some(record) = reader.next_record() {
//!         // Process each record
//!         let record = record?;
//!         let flag = record.flag();
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! ## BQ file format
//!
//! A BINSEQ file consists of two sections:
//!
//! 1. Fixed-size header (32 bytes)
//! 2. Record data section
//!
//! ### Header Format (32 bytes total)
//!
//! | Offset | Size (bytes) | Name     | Description                  | Type   |
//! | ------ | ------------ | -------- | ---------------------------- | ------ |
//! | 0      | 4            | magic    | Magic number (0x42534551)    | uint32 |
//! | 4      | 1            | format   | Format version (currently 2) | uint8  |
//! | 5      | 4            | slen     | Sequence length (primary)    | uint32 |
//! | 9      | 4            | xlen     | Sequence length (secondary)  | uint32 |
//! | 13     | 19           | reserved | Reserved for future use      | bytes  |
//!
//! ### Record Format
//!
//! Each record consists of a:
//!
//! 1. Flag field (8 bytes, uint64)
//! 2. Sequence data (ceil(N/32) \* 8 bytes, where N is sequence length)
//!
//! The flag field is implementation-defined and can be used for filtering, metadata, or other purposes. The placement of the flag field at the start of each record enables efficient filtering without reading sequence data.
//!
//! Total record size = 8 + (ceil(N/32) \* 8) bytes, where N is sequence length
//!
//! ## Encoding
//!
//! - Each nucleotide is encoded using 2 bits:
//!   - A = 00
//!   - C = 01
//!   - G = 10
//!   - T = 11
//! - Non-ATCG characters are **unsupported**.
//! - Sequences are stored in Little-Endian order
//! - The final u64 of sequence data is padded with zeros if the sequence length is not divisible by 32
//!
//! See [`bitnuc`] for 2bit implementation details.
//!
//! ## bq implementation Notes
//!
//! - Sequences are stored in u64 chunks, each holding up to 32 bases
//! - Random access to any record can be calculated as:
//!   - record_size = 8 + (ceil(sequence_length/32) \* 8)
//!   - record_start = 16 + (record_index \* record_size)
//! - Total number of records can be calculated as: (file_size - 16) / record_size
//! - Flag field placement allows for efficient filtering strategies:
//!   - Records can be skipped based on flag values without reading sequence data
//!   - Flag checks can be vectorized for parallel processing
//!   - Memory access patterns are predictable for better cache utilization
//!
//! ## Example Storage Requirements
//!
//! Common sequence lengths:
//!
//! - 32bp reads:
//!   - Sequence: 1 \* 8 = 8 bytes (fits in one u64)
//!   - Flag: 8 bytes
//!   - Total per record: 16 bytes
//! - 100bp reads:
//!   - Sequence: 4 \* 8 = 32 bytes (requires four u64s)
//!   - Flag: 8 bytes
//!   - Total per record: 40 bytes
//! - 150bp reads:
//!   - Sequence: 5 \* 8 = 40 bytes (requires five u64s)
//!   - Flag: 8 bytes
//!   - Total per record: 48 bytes
//!
//! ## Validation
//!
//! Implementations should verify:
//!
//! 1. Correct magic number
//! 2. Compatible version number
//! 3. Sequence length is greater than 0
//! 4. File size minus header (32 bytes) is divisible by the record size
//!
//! ## Future Considerations
//!
//! - The 19 reserved bytes in the header allow for future format extensions
//! - The 64-bit flag field provides space for implementation-specific features such as:
//!   - Quality score summaries
//!   - Filtering flags
//!   - Read group identifiers
//!   - Processing state
//!   - Count data

mod header;
mod reader;
mod writer;

pub use header::{FileHeader, FileHeaderBuilder, SIZE_HEADER};
pub use reader::{MmapReader, RefRecord, StreamReader};
pub use writer::{Encoder, StreamWriter, StreamWriterBuilder, Writer, WriterBuilder};
