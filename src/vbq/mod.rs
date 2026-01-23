//! # VBQ Format
//!
//! VBQ is a high-performance binary format for variable-length nucleotide sequences
//! that optimizes both storage efficiency and parallel processing capabilities.
//!
//! For more information on the format, please refer to our [preprint](https://www.biorxiv.org/content/10.1101/2025.04.08.647863v1).
//!
//! ## Overview
//!
//! VBQ extends the core principles of BINSEQ to accommodate:
//!
//! * **Variable-length sequences**: Unlike BINSEQ which requires fixed-length reads, VBQ can store
//!   sequences of any length, making it suitable for technologies like PacBio and Oxford Nanopore.
//!
//! * **Quality scores**: Optional storage of quality scores alongside nucleotide data when needed.
//!
//! * **Sequence headers**: Optional storage of sequence identifiers/headers with each record.
//!
//! * **Block-based organization**: Data is organized into fixed-size independent record blocks
//!   for efficient parallel processing.
//!
//! * **Compression**: Optional ZSTD compression of individual blocks balances storage
//!   efficiency with processing speed.
//!
//! * **Paired-end support**: Native support for paired sequences without needing multiple files.
//!
//! * **Multi-bit encoding**: Support for 2-bit and 4-bit nucleotide encodings.
//!
//! * **Embedded index**: Self-contained files with embedded index data for efficient random access.
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
//! * Primary header length (8 bytes, if `headers` flag set)
//! * Primary header data (UTF-8 string, if `headers` flag set)
//! * Extended header length (8 bytes, if paired and `headers` flag set)
//! * Extended header data (UTF-8 string, if paired and `headers` flag set)
//!
//! ## Recent Format Changes (v0.7.0+)
//!
//! * **Embedded Index**: Index data is now stored within the VBQ file itself, eliminating
//!   improving portability.
//! * **Headers Support**: Optional sequence identifiers can be stored with each record.
//! * **Extended Capacity**: u64 indexing supports files with more than 4 billion records.
//! * **Multi-bit Encoding**: Support for both 2-bit and 4-bit nucleotide encodings.
//!
//! ## Performance Characteristics
//!
//! VBQ is designed for high-throughput parallel processing:
//!
//! * Independent blocks enable true parallel processing without synchronization
//! * Memory-mapped access provides efficient I/O
//! * Embedded index enables fast random access without auxiliary files
//! * Multi-bit encoding (2-bit/4-bit) optimizes storage for different use cases
//! * Optional ZSTD compression reduces file size with minimal performance impact
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

pub use header::{BlockHeader, FileHeader, FileHeaderBuilder};
pub use index::{BlockIndex, BlockRange};
pub use reader::{MmapReader, RecordBlock, RecordBlockIter, RefRecord};
pub use writer::{Writer, WriterBuilder};
