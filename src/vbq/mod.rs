//! # VBQ Format
//!
//! VBQ is a binary format for variable-length nucleotide sequences with optional
//! quality scores, headers, and per-block ZSTD compression.
//!
//! **VBQ is superseded by [`cbq`](crate::cbq)**, which improves on its compression
//! and read throughput. VBQ remains fully supported for existing files, but new
//! projects should use CBQ.
//!
//! Records are row-based: each record stores its fields contiguously
//! (2-bit or 4-bit encoded sequences), organized into independently
//! compressed blocks with an embedded index for random access.
//! For more information on the format, please refer to our [preprint](https://www.biorxiv.org/content/10.1101/2025.04.08.647863v1).
//!
//! ## File Structure
//!
//! A VBQ file consists of a 32-byte header followed by record blocks and an embedded index:
//!
//! ```text
//! ┌───────────────────┐
//! │    File Header    │ 32 bytes
//! ├───────────────────┤
//! │   Block Header    │ 32 bytes
//! ├───────────────────┤
//! │                   │
//! │   Block Records   │ Variable size
//! │                   │
//! ├───────────────────┤
//! │       ...         │ More blocks
//! ├───────────────────┤
//! │ Compressed Index  │ Variable size
//! ├───────────────────┤
//! │   Index Size      │ 8 bytes (u64)
//! ├───────────────────┤
//! │ Index End Magic   │ 8 bytes
//! └───────────────────┘
//! ```
//!
//! ## Record Format
//!
//! Each record contains the following fields in order:
//!
//! * Flag field (8 bytes)
//! * Primary sequence length (8 bytes)
//! * Extended sequence length (8 bytes, 0 if not paired)
//! * Primary sequence data (2-bit or 4-bit encoded)
//! * Extended sequence data (optional, for paired-end)
//! * Primary quality scores (optional, if `qual` flag set)
//! * Extended quality scores (optional, if paired and `qual` flag set)
//! * Primary header length + data (8 bytes + UTF-8, if `headers` flag set)
//! * Extended header length + data (8 bytes + UTF-8, if paired and `headers` flag set)
//!
//! ## Usage Example
//!
//! ```
//! use std::fs::File;
//! use std::io::BufWriter;
//! use binseq::vbq::{FileHeaderBuilder, WriterBuilder, MmapReader};
//! use binseq::{BinseqRecord, SequencingRecordBuilder};
//!
//! /*
//!    WRITING
//! */
//!
//! // Create a header for sequences with quality scores and headers
//! let header = FileHeaderBuilder::new()
//!     .qual(true)
//!     .compressed(true)
//!     .headers(true)
//!     .build();
//!
//! // Create a writer
//! let file = File::create("example.vbq").unwrap();
//! let mut writer = WriterBuilder::default()
//!     .header(header)
//!     .build(BufWriter::new(file))
//!     .unwrap();
//!
//! // Write a sequence with quality scores and header
//! let record = SequencingRecordBuilder::default()
//!     .s_seq(b"ACGTACGT")
//!     .s_qual(b"IIIIFFFF")
//!     .s_header(b"sequence_001")
//!     .build()
//!     .unwrap();
//! writer.push(record).unwrap();
//! writer.finish().unwrap();
//!
//! /*
//!    READING
//! */
//!
//! // Read the sequences back
//! let mut reader = MmapReader::new("example.vbq").unwrap();
//! let mut block = reader.new_block();
//!
//! // Process blocks one at a time
//! let mut seq_buffer = Vec::new();
//! while reader.read_block_into(&mut block).unwrap() {
//!     for record in block.iter() {
//!         record.decode_s(&mut seq_buffer).unwrap();
//!         let header = record.sheader();
//!         println!("Header: {}", std::str::from_utf8(header).unwrap());
//!         println!("Sequence: {}", std::str::from_utf8(&seq_buffer).unwrap());
//!         println!("Quality: {}", std::str::from_utf8(record.squal()).unwrap());
//!         seq_buffer.clear();
//!     }
//! }
//! # std::fs::remove_file("example.vbq").unwrap_or(());
//! ```

mod header;
mod index;
mod reader;
mod writer;

pub use header::{BlockHeader, FILE_MAGIC, FileHeader, FileHeaderBuilder};
pub use index::{BlockIndex, BlockRange};
pub use reader::{MmapReader, RecordBlock, RecordBlockIter, RefRecord};
pub use writer::{Writer, WriterBuilder};
