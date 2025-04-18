#![doc = include_str!("../README.md")]
//!
//! # Overview
//!
//! The `binseq` library provides efficient tools for working with binary-encoded
//! nucleotide sequences. It offers:
//!
//! - Compact 2-bit encoding of nucleotide sequences
//! - Memory-mapped file access for efficient reading
//! - Parallel processing capabilities
//! - Configurable policies for handling invalid nucleotides
//! - Support for both single and paired-end sequences
//! - SIMD acceleration for faster encoding/decoding (on supported CPUs)
//!
//! # Core Components
//!
//! - [`BinseqWriter`]: Writes sequences to binary format
//! - [`MmapReader`]: Reads sequences using memory mapping
//! - [`BinseqHeader`]: Defines file format and sequence lengths
//! - [`Policy`]: Configures invalid nucleotide handling
//! - [`ParallelProcessor`]: Enables parallel sequence processing
//!
//! # Example
//!
//! ```
//! use binseq::{BinseqHeader, BinseqWriterBuilder, MmapReader, Policy, Result};
//! use std::io::Cursor;
//!
//! fn main() -> Result<()> {
//!     // Create a writer for sequences of length 100
//!     let header = BinseqHeader::new(100);
//!
//!     let mut writer = BinseqWriterBuilder::default()
//!         .header(header)
//!         .build(Cursor::new(Vec::new()))?;
//!
//!     // Write a sequence
//!     let sequence = b"ACGT".repeat(25); // 100 nucleotides
//!     writer.write_nucleotides(0, &sequence)?;
//!
//!     Ok(())
//! }
//! ```

#![allow(clippy::module_inception)]

pub mod error;
mod header;
pub mod nuc;
mod parallel;
mod policy;
mod reader;
mod simd;
mod utils;
pub mod writer;

pub use error::{Error, HeaderError, ReadError, Result, WriteError};
pub use header::{BinseqHeader, SIZE_HEADER};
pub use parallel::ParallelProcessor;
pub use policy::{Policy, RNG_SEED};
pub use reader::{MmapReader, RefRecord};
pub use simd::is_simd_supported;
pub use utils::expected_file_size;
pub use writer::{BinseqWriter, BinseqWriterBuilder, Encoder};

// #[cfg(test)]
// mod testing {

//     use super::*;
//     use anyhow::Result;
//     use nucgen::Sequence;
//     use std::io::Cursor;

//     #[test]
//     fn test_binseq_short() -> Result<()> {
//         let header = BinseqHeader::new(16);
//         let mut writer = BinseqWriter::new(Cursor::new(Vec::new()), header)?;

//         let sequence = b"ACGTACGTACGTACGT";
//         writer.write_nucleotides(0, sequence)?;

//         let cursor = writer.into_inner().into_inner();
//         let mut reader = SingleReader::new(cursor.as_slice())?;
//         let record = reader.next().unwrap()?;
//         assert_eq!(record.flag(), 0);
//         let bitseq = record.sequence()[0];
//         let readout = bitnuc::from_2bit_alloc(bitseq, 16)?;
//         assert_eq!(&readout, sequence);

//         Ok(())
//     }

//     #[test]
//     fn test_binseq_short_multiple() -> Result<()> {
//         let header = BinseqHeader::new(16);
//         let mut writer = BinseqWriter::new(Cursor::new(Vec::new()), header)?;

//         let sequence = b"ACGTACGTACGTACGT";
//         writer.write_nucleotides(0, sequence)?;
//         writer.write_nucleotides(0, sequence)?;
//         writer.write_nucleotides(0, sequence)?; // write 3 times

//         let cursor = writer.into_inner().into_inner();
//         let mut reader = SingleReader::new(cursor.as_slice())?;

//         for _ in 0..3 {
//             let record = reader.next().unwrap()?;
//             assert_eq!(record.flag(), 0);
//             let bitseq = record.sequence()[0];
//             let dbuf = bitnuc::from_2bit_alloc(bitseq, 16)?;
//             assert_eq!(&dbuf, sequence);
//         }

//         Ok(())
//     }

//     #[test]
//     fn test_binseq_long() -> Result<()> {
//         let header = BinseqHeader::new(40);
//         let mut writer = BinseqWriter::new(Cursor::new(Vec::new()), header)?;

//         let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
//         writer.write_nucleotides(0, sequence)?;

//         let cursor = writer.into_inner().into_inner();
//         let mut reader = SingleReader::new(cursor.as_slice())?;
//         let record = reader.next().unwrap()?;
//         assert_eq!(record.flag(), 0);

//         let dbuf = record.decode_alloc()?;
//         assert_eq!(&dbuf, sequence);

//         Ok(())
//     }

//     #[test]
//     fn test_binseq_long_multiple() -> Result<()> {
//         let header = BinseqHeader::new(40);
//         let mut writer = BinseqWriter::new(Cursor::new(Vec::new()), header)?;

//         let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
//         writer.write_nucleotides(0, sequence)?;
//         writer.write_nucleotides(0, sequence)?;
//         writer.write_nucleotides(0, sequence)?; // write 3 times

//         let cursor = writer.into_inner().into_inner();
//         let mut reader = SingleReader::new(cursor.as_slice())?;

//         for _ in 0..3 {
//             let record = reader.next().unwrap()?;
//             assert_eq!(record.flag(), 0);

//             let dbuf = record.decode_alloc()?;
//             assert_eq!(&dbuf, sequence);
//         }

//         Ok(())
//     }

//     #[test]
//     fn test_n_in_sequence() -> Result<()> {
//         let header = BinseqHeader::new(40);
//         let mut writer = BinseqWriter::new(Cursor::new(Vec::new()), header)?;

//         let sequence = b"ACGTACGTACGTACNTACGTACGTACGTACGTACGTACGT";
//         writer.write_nucleotides(0, sequence)?;

//         let cursor = writer.into_inner().into_inner();
//         let mut reader = SingleReader::new(cursor.as_slice())?;
//         let record = reader.next();
//         dbg!(&record);
//         assert!(record.is_none());

//         Ok(())
//     }

//     #[test]
//     fn test_n_in_sequence_policy_ignore() -> Result<()> {
//         let header = BinseqHeader::new(10);
//         let sequence = b"NNNNNNNNNN";

//         let mut writer =
//             BinseqWriter::new_with_policy(Cursor::new(Vec::new()), header, Policy::IgnoreSequence)?;
//         writer.write_nucleotides(0, sequence)?;
//         let cursor = writer.into_inner().into_inner();
//         let mut reader = SingleReader::new(cursor.as_slice())?;
//         let record = reader.next();
//         assert!(record.is_none());
//         Ok(())
//     }

//     #[test]
//     fn test_n_in_sequence_policy_error() -> Result<()> {
//         let header = BinseqHeader::new(10);
//         let sequence = b"NNNNNNNNNN";

//         let mut writer =
//             BinseqWriter::new_with_policy(Cursor::new(Vec::new()), header, Policy::BreakOnInvalid)?;
//         let result = writer.write_nucleotides(0, sequence);
//         assert!(result.is_err());
//         Ok(())
//     }

//     #[test]
//     fn test_n_in_sequence_policy_random() -> Result<()> {
//         let header = BinseqHeader::new(10);
//         let sequence = b"NNNNNNNNNN";

//         let mut writer =
//             BinseqWriter::new_with_policy(Cursor::new(Vec::new()), header, Policy::RandomDraw)?;
//         writer.write_nucleotides(0, sequence)?;

//         let cursor = writer.into_inner().into_inner();
//         let mut reader = SingleReader::new(cursor.as_slice())?;
//         let record = reader.next().unwrap()?;
//         assert_eq!(record.flag(), 0);
//         let dbuf = record.decode_alloc()?;
//         assert_eq!(dbuf.len(), 10);
//         for c in dbuf.iter() {
//             assert_ne!(*c, b'N');
//         }
//         Ok(())
//     }

//     #[test]
//     fn test_n_in_sequence_policy_set_to_a() -> Result<()> {
//         let header = BinseqHeader::new(10);
//         let sequence = b"NNNNNNNNNN";

//         let mut writer =
//             BinseqWriter::new_with_policy(Cursor::new(Vec::new()), header, Policy::SetToA)?;
//         writer.write_nucleotides(0, sequence)?;

//         let cursor = writer.into_inner().into_inner();
//         let mut reader = SingleReader::new(cursor.as_slice())?;
//         let record = reader.next().unwrap()?;
//         assert_eq!(record.flag(), 0);
//         let dbuf = record.decode_alloc()?;
//         assert_eq!(dbuf.len(), 10);
//         for c in dbuf.iter() {
//             assert_eq!(*c, b'A');
//         }
//         Ok(())
//     }

//     #[test]
//     fn test_n_in_sequence_policy_set_to_c() -> Result<()> {
//         let header = BinseqHeader::new(10);
//         let sequence = b"NNNNNNNNNN";

//         let mut writer =
//             BinseqWriter::new_with_policy(Cursor::new(Vec::new()), header, Policy::SetToC)?;
//         writer.write_nucleotides(0, sequence)?;

//         let cursor = writer.into_inner().into_inner();
//         let mut reader = SingleReader::new(cursor.as_slice())?;
//         let record = reader.next().unwrap()?;
//         assert_eq!(record.flag(), 0);
//         let dbuf = record.decode_alloc()?;
//         assert_eq!(dbuf.len(), 10);
//         for c in dbuf.iter() {
//             assert_eq!(*c, b'C');
//         }
//         Ok(())
//     }

//     #[test]
//     fn test_n_in_sequence_policy_set_to_g() -> Result<()> {
//         let header = BinseqHeader::new(10);
//         let sequence = b"NNNNNNNNNN";

//         let mut writer =
//             BinseqWriter::new_with_policy(Cursor::new(Vec::new()), header, Policy::SetToG)?;
//         writer.write_nucleotides(0, sequence)?;

//         let cursor = writer.into_inner().into_inner();
//         let mut reader = SingleReader::new(cursor.as_slice())?;
//         let record = reader.next().unwrap()?;
//         assert_eq!(record.flag(), 0);
//         let dbuf = record.decode_alloc()?;
//         assert_eq!(dbuf.len(), 10);
//         for c in dbuf.iter() {
//             assert_eq!(*c, b'G');
//         }
//         Ok(())
//     }

//     #[test]
//     fn test_n_in_sequence_policy_set_to_t() -> Result<()> {
//         let header = BinseqHeader::new(10);
//         let sequence = b"NNNNNNNNNN";

//         let mut writer =
//             BinseqWriter::new_with_policy(Cursor::new(Vec::new()), header, Policy::SetToT)?;
//         writer.write_nucleotides(0, sequence)?;

//         let cursor = writer.into_inner().into_inner();
//         let mut reader = SingleReader::new(cursor.as_slice())?;
//         let record = reader.next().unwrap()?;
//         assert_eq!(record.flag(), 0);
//         let dbuf = record.decode_alloc()?;
//         assert_eq!(dbuf.len(), 10);
//         for c in dbuf.iter() {
//             assert_eq!(*c, b'T');
//         }
//         Ok(())
//     }

//     fn valid_reconstruction(seq_len: usize, num_records: usize) -> Result<()> {
//         let mut rng = rand::thread_rng();
//         let mut sequence = Sequence::new();

//         // stores the original sequences
//         let mut seq_vec = Vec::new();

//         // write the sequences to a binseq file
//         // and store the original sequences
//         let header = BinseqHeader::new(seq_len as u32);
//         let mut writer = BinseqWriter::new(Cursor::new(Vec::new()), header)?;
//         for _ in 0..num_records {
//             sequence.fill_buffer(&mut rng, seq_len);
//             seq_vec.push(sequence.bytes().to_vec());
//             writer.write_nucleotides(0, sequence.bytes())?;
//         }

//         // Verify that the file size is as expected
//         let cursor = writer.into_inner().into_inner();
//         let file_size = cursor.len();
//         let expected_size = expected_file_size(num_records, seq_len);
//         assert_eq!(file_size, expected_size);

//         // read the sequences back from the binseq file
//         // and compare them to the original sequences
//         // stored in seq_vec
//         let mut reader = SingleReader::new(cursor.as_slice())?;
//         for seq in seq_vec.iter() {
//             let record = reader.next().unwrap()?;
//             assert_eq!(record.flag(), 0);
//             let dbuf = record.decode_alloc()?;
//             assert_eq!(&dbuf, seq);
//         }
//         assert!(reader.next().is_none());

//         Ok(())
//     }

//     #[test]
//     fn test_reconstruction() -> Result<()> {
//         // test various sequence lengths
//         for n_bases in [10, 32, 50, 64, 100, 1000] {
//             for n_records in [1, 10, 32, 100, 256, 1000] {
//                 valid_reconstruction(n_bases, n_records)?;
//             }
//         }
//         Ok(())
//     }
// }
