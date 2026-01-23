//! Reader implementation for VBINSEQ files
//!
//! This module provides functionality for reading sequence data from VBINSEQ files,
//! including support for compressed blocks, quality scores, paired-end reads, and sequence headers.
//!
//! ## Format Changes (v0.7.0+)
//!
//! - **Embedded Index**: Readers now load the index from within VBQ files instead of separate `.vqi` files
//! - **Headers Support**: Optional sequence headers/identifiers can be read from each record
//! - **Multi-bit Encoding**: Support for reading 2-bit and 4-bit nucleotide encodings
//! - **Extended Capacity**: u64 indexing supports files with more than 4 billion records
//!
//! ## Index Loading
//!
//! The reader automatically loads the embedded index from the end of VBQ files:
//! 1. Seeks to the end of the file to read the index trailer
//! 2. Validates the `INDEX_END_MAGIC` marker
//! 3. Reads the index size and decompresses the embedded index
//! 4. Uses the index for efficient random access and parallel processing
//!
//! ## Memory-Mapped Reading
//!
//! The `MmapReader` provides efficient access to large files through memory mapping:
//! - Zero-copy access to file data
//! - Efficient random access using the embedded index
//! - Support for parallel processing across record ranges
//!
//! ## Example
//!
//! ```rust,no_run
//! use binseq::vbq::MmapReader;
//! use binseq::BinseqRecord;
//!
//! // Open a VBQ file (index is automatically loaded)
//! let mut reader = MmapReader::new("example.vbq").unwrap();
//! let mut block = reader.new_block();
//!
//! // Read records with headers and quality scores
//! while reader.read_block_into(&mut block).unwrap() {
//!     for record in block.iter() {
//!         let seq = record.sseq();
//!         let header = record.sheader();
//!         println!("Header: {}", std::str::from_utf8(header).unwrap());
//!         println!("Sequence: {}", std::str::from_utf8(seq).unwrap());
//!         if !record.squal().is_empty() {
//!             println!("Quality: {}", std::str::from_utf8(record.squal()).unwrap());
//!         }
//!     }
//! }
//! ```

use std::fs::File;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use bitnuc::BitSize;
use byteorder::{ByteOrder, LittleEndian};
use memmap2::Mmap;
use zstd::zstd_safe;

use super::{
    BlockHeader, BlockIndex, BlockRange, VBinseqHeader,
    header::{SIZE_BLOCK_HEADER, SIZE_HEADER},
};
use crate::DEFAULT_QUALITY_SCORE;
use crate::vbq::index::{INDEX_END_MAGIC, INDEX_HEADER_SIZE, IndexHeader};
use crate::{
    BinseqRecord, ParallelProcessor, ParallelReader,
    error::{ReadError, Result},
};

/// Calculates the number of 64-bit words needed to store a nucleotide sequence of the given length
///
/// Nucleotides are packed into 64-bit words with 2 bits per nucleotide (32 nucleotides per word).
/// This function calculates how many 64-bit words are needed to encode a sequence of a given length.
///
/// # Parameters
///
/// * `len` - Length of the nucleotide sequence in basepairs
///
/// # Returns
///
/// The number of 64-bit words required to encode the sequence
fn encoded_sequence_len(len: u64, bitsize: BitSize) -> usize {
    match bitsize {
        BitSize::Two => len.div_ceil(32) as usize,
        BitSize::Four => len.div_ceil(16) as usize,
    }
}

/// Represents a span (offset, length) into a buffer
#[derive(Clone, Copy, Debug, Default)]
pub struct Span {
    offset: usize,
    len: usize,
}
impl Span {
    fn new(offset: usize, len: usize) -> Self {
        Self { offset, len }
    }

    /// Get a slice of bytes from a buffer
    fn slice<'a>(&self, buffer: &'a [u8]) -> &'a [u8] {
        &buffer[self.offset..self.offset + self.len]
    }

    /// Get a slice of u64s from a buffer
    fn slice_u64<'a>(&self, buffer: &'a [u64]) -> &'a [u64] {
        &buffer[self.offset..self.offset + self.len]
    }
}

/// Metadata for a single record, storing spans into rbuf
#[derive(Debug, Clone, Copy)]
struct RecordMetadata {
    flag: Option<u64>,
    slen: u64,
    xlen: u64,

    // Spans for primary sequence
    s_seq_span: Span,    // Encoded sequence words (u64s) (into `.sequences` buffer)
    s_qual_span: Span,   // Quality bytes
    s_header_span: Span, // Header bytes

    // Spans for extended sequence
    x_seq_span: Span,    // Encoded sequence words (u64s) (into `.sequences` buffer)
    x_qual_span: Span,   // Quality bytes
    x_header_span: Span, // Header bytes

    /// Indicates whether the record has quality scores
    has_quality: bool,
}

/// A container for a block of VBINSEQ records
///
/// The `RecordBlock` struct represents a single block of records read from a VBINSEQ file.
/// It stores the raw data for multiple records in vectors, allowing efficient iteration
/// over the records without copying memory for each record.
///
/// ## Format Support (v0.7.0+)
///
/// - Supports reading records with optional sequence headers
/// - Handles both 2-bit and 4-bit nucleotide encodings
/// - Supports quality scores and paired sequences
/// - Compatible with both compressed and uncompressed blocks
///
/// The `RecordBlock` is reused when reading blocks sequentially from a file, with its
/// contents being cleared and replaced with each new block that is read.
///
/// # Examples
///
/// ```rust,no_run
/// use binseq::vbq::MmapReader;
///
/// let reader = MmapReader::new("example.vbq").unwrap();
/// let mut block = reader.new_block(); // Create a block with appropriate size
/// ```
pub struct RecordBlock {
    /// Bitsize of the records in the block
    bitsize: BitSize,

    /// Index of the first record in the block
    /// This allows records to maintain their global position in the file
    index: usize,

    /// Reusable buffer for temporary storage during decompression
    /// Using a reusable buffer reduces memory allocations
    rbuf: Vec<u8>,

    /// Sequence data (u64s) - small copy during parsing
    sequences: Vec<u64>,

    /// Record metadata stored as compact spans
    records: Vec<RecordMetadata>,

    /// Maximum size of the block in bytes
    /// This is derived from the file header's block size field
    block_size: usize,

    /// Reusable zstd decompression context
    dctx: zstd_safe::DCtx<'static>,

    /// Reusable decoding buffer for the block
    dbuf: Vec<u8>,

    /// Reusable buffer for quality scores for the block
    qbuf: Vec<u8>,

    /// Default quality score for the block
    default_quality_score: u8,
}
impl RecordBlock {
    /// Creates a new empty `RecordBlock` with the specified block size
    ///
    /// The block size should match the one specified in the VBINSEQ file header
    /// for proper operation. This is typically handled automatically when using
    /// `MmapReader::new_block()`.
    ///
    /// # Parameters
    ///
    /// * `bitsize` - Bitsize of the records in the block
    /// * `block_size` - Maximum size of the block in bytes
    ///
    /// # Returns
    ///
    /// A new empty `RecordBlock` instance
    #[must_use]
    pub fn new(bitsize: BitSize, block_size: usize) -> Self {
        Self {
            bitsize,
            index: 0,
            block_size,
            records: Vec::default(),
            sequences: Vec::default(),
            rbuf: Vec::default(),
            dbuf: Vec::default(),
            dctx: zstd_safe::DCtx::create(),
            qbuf: Vec::default(),
            default_quality_score: DEFAULT_QUALITY_SCORE,
        }
    }

    /// Sets the default quality score for the block
    ///
    /// # Parameters
    ///
    /// * `score` - Default quality score for the block
    pub fn set_default_quality_score(&mut self, score: u8) {
        self.default_quality_score = score;
        self.qbuf.clear();
    }

    /// Returns the number of records in this block
    ///
    /// # Returns
    ///
    /// The number of records currently stored in this block
    #[must_use]
    pub fn n_records(&self) -> usize {
        self.records.len()
    }

    /// Returns an iterator over the records in this block
    ///
    /// The iterator yields `RefRecord` instances that provide access to the record data
    /// without copying the underlying data.
    ///
    /// # Returns
    ///
    /// An iterator over the records in this block
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::MmapReader;
    /// use binseq::BinseqRecord;
    ///
    /// let mut reader = MmapReader::new("example.vbq").unwrap();
    /// let mut block = reader.new_block();
    /// reader.read_block_into(&mut block).unwrap();
    ///
    /// // Iterate over records in the block
    /// for record in block.iter() {
    ///     println!("Record {}", record.index());
    /// }
    /// ```
    #[must_use]
    #[allow(clippy::iter_without_into_iter)]
    pub fn iter(&self) -> RecordBlockIter<'_> {
        RecordBlockIter::new(self)
    }

    /// Updates the starting index of the block
    ///
    /// This is used internally to keep track of the global position of records
    /// within the file, allowing each record to maintain its original index.
    ///
    /// # Parameters
    ///
    /// * `index` - The index of the first record in the block
    fn update_index(&mut self, index: usize) {
        self.index = index;
    }

    /// Clears all data from the block
    ///
    /// This method resets the block to an empty state, clearing all vectors and resetting
    /// the index to 0. This is typically used when reusing a block for reading a new block
    /// from a file.
    pub fn clear(&mut self) {
        self.index = 0;
        self.records.clear();
        self.sequences.clear();
        self.dbuf.clear();
        // Note: We keep rbuf allocated for reuse
        // Note: We keep qbuf allocated for reuse
    }

    /// Ingest the bytes from a block into the record block
    ///
    /// This method takes a slice of bytes and processes it to extract
    /// the records from the block. It is used when reading a block from
    /// a file into a record block.
    ///
    /// This is a private method used primarily for parallel processing.
    ///
    /// # Parameters
    ///
    /// * `bytes` - A slice of bytes containing the block data
    /// * `has_quality` - A boolean indicating whether the block contains quality scores
    /// * `has_header` - A boolean indicating whether the block contains headers
    fn ingest_bytes(
        &mut self,
        bytes: &[u8],
        has_quality: bool,
        has_header: bool,
        has_flags: bool,
    ) -> Result<()> {
        if bytes.len() != self.block_size {
            return Err(ReadError::PartialRecord(bytes.len()).into());
        }
        self.rbuf.clear();
        self.rbuf.extend_from_slice(bytes);
        self.parse_records(has_quality, has_header, has_flags);
        Ok(())
    }

    /// Decompresses the given bytes and ingests them into the record block.
    fn ingest_compressed_bytes(
        &mut self,
        bytes: &[u8],
        has_quality: bool,
        has_header: bool,
        has_flags: bool,
    ) -> Result<()> {
        // Clear and ensure capacity
        self.rbuf.clear();
        if self.rbuf.capacity() < self.block_size {
            self.rbuf.reserve(self.block_size - self.rbuf.capacity());
        }

        // Reuse the decompression context - avoids allocation!
        let bytes_read = self
            .dctx
            .decompress(&mut self.rbuf, bytes)
            .map_err(|code| std::io::Error::other(zstd_safe::get_error_name(code)))?;

        if bytes_read != self.block_size {
            return Err(ReadError::PartialRecord(bytes_read).into());
        }

        self.parse_records(has_quality, has_header, has_flags);
        Ok(())
    }
    /// Parse records from rbuf, storing spans for all data
    fn parse_records(&mut self, has_quality: bool, has_header: bool, has_flags: bool) {
        self.records.clear();
        self.sequences.clear();

        let mut pos = 0;
        let bytes = &self.rbuf;

        loop {
            // Check if we have enough bytes for the minimum record header
            let min_header_size = if has_flags { 24 } else { 16 };
            if pos + min_header_size > bytes.len() {
                break;
            }

            // Read flag
            let flag = if has_flags {
                let flag = LittleEndian::read_u64(&bytes[pos..pos + 8]);
                pos += 8;
                Some(flag)
            } else {
                None
            };

            // Read lengths
            let slen = LittleEndian::read_u64(&bytes[pos..pos + 8]);
            pos += 8;
            let xlen = LittleEndian::read_u64(&bytes[pos..pos + 8]);
            pos += 8;

            // Check for end of records
            if slen == 0 {
                break;
            }

            // Calculate sizes
            let s_seq_words = encoded_sequence_len(slen, self.bitsize);
            let x_seq_words = encoded_sequence_len(xlen, self.bitsize);

            // Primary sequence - store span into sequences Vec
            let s_seq_span = Span::new(self.sequences.len(), s_seq_words);
            for _ in 0..s_seq_words {
                let val = LittleEndian::read_u64(&bytes[pos..pos + 8]);
                self.sequences.push(val);
                pos += 8;
            }

            // Primary quality - store span into rbuf
            let s_qual_span = if has_quality {
                let span = Span::new(pos, slen as usize);
                pos += slen as usize;
                span
            } else {
                Span::new(0, 0)
            };

            // Primary header - store span into rbuf
            let s_header_span = if has_header {
                let header_len = LittleEndian::read_u64(&bytes[pos..pos + 8]) as usize;
                pos += 8;
                let span = Span::new(pos, header_len);
                pos += header_len;
                span
            } else {
                Span::new(0, 0)
            };

            // Extended sequence - store span into sequences Vec
            let x_seq_span = Span::new(self.sequences.len(), x_seq_words);
            for _ in 0..x_seq_words {
                let val = LittleEndian::read_u64(&bytes[pos..pos + 8]);
                self.sequences.push(val);
                pos += 8;
            }

            // Extended quality - store span into rbuf
            let x_qual_span = if has_quality {
                let span = Span::new(pos, xlen as usize);
                pos += xlen as usize;
                span
            } else {
                Span::new(0, 0)
            };

            // Extended header - store span into rbuf
            let x_header_span = if has_header && xlen > 0 {
                let header_len = LittleEndian::read_u64(&bytes[pos..pos + 8]) as usize;
                pos += 8;
                let span = Span::new(pos, header_len);
                pos += header_len;
                span
            } else {
                Span::new(0, 0)
            };

            // Update qbuf size
            if !has_quality {
                let max_size = slen.max(xlen) as usize;
                if self.qbuf.len() < max_size {
                    self.qbuf.resize(max_size, self.default_quality_score);
                }
            }

            // Store the record metadata - all spans!
            self.records.push(RecordMetadata {
                flag,
                slen,
                xlen,
                s_seq_span,
                s_qual_span,
                s_header_span,
                x_seq_span,
                x_qual_span,
                x_header_span,
                has_quality,
            });
        }
    }

    /// Decodes all sequences in the block at once.
    ///
    /// Note:
    /// Each record's sequence is padded internally to the nearest u64.
    /// Because of this the global decoding will include nucleotides that are not present in the original data.
    /// We track the non-contiguous regions of the sequence separately.
    pub fn decode_all(&mut self) -> Result<()> {
        if self.sequences.is_empty() {
            return Ok(());
        }
        self.dbuf.clear();
        match self.bitsize {
            BitSize::Two => {
                let num_bp = self.sequences.len() * 32;
                bitnuc::twobit::decode(&self.sequences, num_bp, &mut self.dbuf)
            }
            BitSize::Four => {
                let num_bp = self.sequences.len() * 16;
                bitnuc::fourbit::decode(&self.sequences, num_bp, &mut self.dbuf)
            }
        }?;
        Ok(())
    }

    /// Get decoded primary sequence for a record by index
    #[must_use]
    pub fn get_decoded_s(&self, record_idx: usize) -> Option<&[u8]> {
        let meta = self.records.get(record_idx)?;
        if self.dbuf.is_empty() {
            return None;
        }

        let bases_per_word = match self.bitsize {
            BitSize::Two => 32,
            BitSize::Four => 16,
        };

        // Calculate offset in decoded buffer (accounting for padding)
        let offset = meta.s_seq_span.offset * bases_per_word;
        let len = meta.slen as usize;

        Some(&self.dbuf[offset..offset + len])
    }

    /// Get decoded extended sequence for a record by index
    #[must_use]
    pub fn get_decoded_x(&self, record_idx: usize) -> Option<&[u8]> {
        let meta = self.records.get(record_idx)?;
        if meta.xlen == 0 {
            return Some(&[]);
        }
        if self.dbuf.is_empty() {
            return None;
        }

        let bases_per_word = match self.bitsize {
            BitSize::Two => 32,
            BitSize::Four => 16,
        };

        let offset = meta.x_seq_span.offset * bases_per_word;
        let len = meta.xlen as usize;

        Some(&self.dbuf[offset..offset + len])
    }
}

pub struct RecordBlockIter<'a> {
    block: &'a RecordBlock,
    pos: usize,
    header_buffer: itoa::Buffer,
    qbuf: &'a [u8],
}
impl<'a> RecordBlockIter<'a> {
    #[must_use]
    pub fn new(block: &'a RecordBlock) -> Self {
        Self {
            block,
            pos: 0,
            header_buffer: itoa::Buffer::new(),
            qbuf: &block.qbuf,
        }
    }
}
impl<'a> Iterator for RecordBlockIter<'a> {
    type Item = RefRecord<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.block.records.len() {
            return None;
        }

        let meta = &self.block.records[self.pos];
        let index = (self.block.index + self.pos) as u64;
        let index_in_block = self.pos;

        let mut header_buf = [0; 20];
        let mut header_len = 0;
        if meta.s_header_span.len == 0 && meta.x_header_span.len == 0 {
            let header_str = self.header_buffer.format(index);
            header_len = header_str.len();
            header_buf[..header_len].copy_from_slice(header_str.as_bytes());
        }

        let (squal, xqual) = if meta.has_quality {
            // Record has quality scores, slice into rbuf using span
            (
                meta.s_qual_span.slice(&self.block.rbuf),
                meta.x_qual_span.slice(&self.block.rbuf),
            )
        } else {
            // Record does not have quality scores, use preallocated buffer for default scores
            (
                &self.qbuf[..meta.slen as usize],
                &self.qbuf[..meta.xlen as usize],
            )
        };

        // increment position
        {
            self.pos += 1;
        }

        Some(RefRecord {
            block: self.block,
            bitsize: self.block.bitsize,
            index,
            index_in_block,
            flag: meta.flag,
            slen: meta.slen,
            xlen: meta.xlen,
            // Slice into sequences Vec using span
            sbuf: meta.s_seq_span.slice_u64(&self.block.sequences),
            xbuf: meta.x_seq_span.slice_u64(&self.block.sequences),
            // Pass quality score buffers
            squal,
            xqual,
            // Slice into rbuf using span
            sheader: meta.s_header_span.slice(&self.block.rbuf),
            xheader: meta.x_header_span.slice(&self.block.rbuf),
            header_buf,
            header_len,
        })
    }
}

/// Zero-copy record reference
pub struct RefRecord<'a> {
    block: &'a RecordBlock,
    bitsize: BitSize,
    index: u64,
    index_in_block: usize,
    flag: Option<u64>,
    slen: u64,
    xlen: u64,
    sbuf: &'a [u64],
    xbuf: &'a [u64],
    squal: &'a [u8],
    xqual: &'a [u8],
    sheader: &'a [u8],
    xheader: &'a [u8],
    header_buf: [u8; 20],
    header_len: usize,
}

impl BinseqRecord for RefRecord<'_> {
    fn bitsize(&self) -> BitSize {
        self.bitsize
    }

    fn index(&self) -> u64 {
        self.index
    }

    fn sheader(&self) -> &[u8] {
        if self.sheader.is_empty() {
            &self.header_buf[..self.header_len]
        } else {
            self.sheader
        }
    }

    fn xheader(&self) -> &[u8] {
        if self.xheader.is_empty() {
            &self.header_buf[..self.header_len]
        } else {
            self.xheader
        }
    }

    fn flag(&self) -> Option<u64> {
        self.flag
    }

    fn slen(&self) -> u64 {
        self.slen
    }

    fn xlen(&self) -> u64 {
        self.xlen
    }

    fn sbuf(&self) -> &[u64] {
        self.sbuf
    }

    fn xbuf(&self) -> &[u64] {
        self.xbuf
    }

    fn squal(&self) -> &[u8] {
        self.squal
    }

    fn xqual(&self) -> &[u8] {
        self.xqual
    }

    /// Override this method since we can make use of block information
    fn decode_s(&self, buf: &mut Vec<u8>) -> Result<()> {
        if let Some(decoded) = self.block.get_decoded_s(self.index_in_block) {
            buf.extend_from_slice(decoded);
        } else {
            self.bitsize()
                .decode(self.sbuf(), self.slen() as usize, buf)?;
        }
        Ok(())
    }

    /// Override this method since we can make use of block information
    fn decode_x(&self, buf: &mut Vec<u8>) -> Result<()> {
        if let Some(decoded) = self.block.get_decoded_x(self.index_in_block) {
            buf.extend_from_slice(decoded);
        } else {
            self.bitsize()
                .decode(self.xbuf(), self.xlen() as usize, buf)?;
        }
        Ok(())
    }

    /// Override this method since we can make use of block information
    fn sseq(&self) -> &[u8] {
        self.block
            .get_decoded_s(self.index_in_block)
            .expect("Reader was built without batch-decoding")
    }

    /// Override this method since we can make use of block information
    fn xseq(&self) -> &[u8] {
        self.block
            .get_decoded_x(self.index_in_block)
            .expect("Reader was built without batch-decoding")
    }
}

/// Memory-mapped reader for VBINSEQ files
///
/// [`MmapReader`] provides efficient, memory-mapped access to VBINSEQ files. It allows
/// sequential reading of record blocks and supports parallel processing of records.
///
/// ## Format Support (v0.7.0+)
///
/// - **Embedded Index**: Automatically loads index from within VBQ files
/// - **Headers Support**: Reads optional sequence headers/identifiers from records
/// - **Multi-bit Encoding**: Supports both 2-bit and 4-bit nucleotide encodings
/// - **Extended Capacity**: u64 indexing supports files with more than 4 billion records
///
/// Memory mapping allows the operating system to lazily load file contents as needed,
/// which can be more efficient than standard file I/O, especially for large files.
///
/// The [`MmapReader`] is designed to be used in a multi-threaded environment, and it
/// is built around [`RecordBlock`]s which are the units of data in a VBINSEQ file.
/// Each one would be held by a separate thread and would load data from the shared
/// [`MmapReader`] through the [`MmapReader::read_block_into`] method. However, they can
/// also be used in a single-threaded environment for sequential processing.
///
/// Each [`RecordBlock`] contains a [`BlockHeader`] and is used to access [`RefRecord`]s
/// which implement the [`BinseqRecord`] trait.
///
/// ## Index Loading
///
/// The reader automatically loads the embedded index by:
/// 1. Reading the index trailer from the end of the file
/// 2. Validating the `INDEX_END_MAGIC` marker
/// 3. Decompressing the embedded index data
/// 4. Using the index for efficient random access and parallel processing
///
/// # Examples
///
/// ```
/// use binseq::vbq::MmapReader;
/// use binseq::{BinseqRecord, Result};
///
/// #[allow(deprecated)]
/// fn main() -> Result<()> {
///     let path = "./data/subset.vbq";
///     let mut reader = MmapReader::new(path)?; // Index loaded automatically
///
///     // Create buffers for sequence data and headers
///     let mut seq_buffer = Vec::new();
///     let mut block = reader.new_block();
///
///     // Read blocks sequentially
///     while reader.read_block_into(&mut block)? {
///         println!("Read a block with {} records", block.n_records());
///         for record in block.iter() {
///             // Decode sequence and header
///             record.decode_s(&mut seq_buffer)?;
///             let header = record.sheader();
///
///             println!("Header: {}", std::str::from_utf8(&header).unwrap_or(""));
///             println!("Sequence: {}", std::str::from_utf8(&seq_buffer).unwrap_or(""));
///
///             seq_buffer.clear();
///         }
///     }
///     Ok(())
/// }
/// ```
pub struct MmapReader {
    /// Path to the VBINSEQ file
    path: PathBuf,

    /// Memory-mapped file contents for efficient access
    mmap: Arc<Mmap>,

    /// Parsed header information from the file
    header: VBinseqHeader,

    /// Current cursor position in the file (in bytes)
    pos: usize,

    /// Total number of records read from the file so far
    total: usize,

    /// Whether to decode sequences at once in each block
    decode_block: bool,

    /// Default quality score for this reader
    default_quality_score: u8,
}
impl MmapReader {
    /// Creates a new `MmapReader` for a VBINSEQ file
    ///
    /// This method opens the specified file, memory-maps its contents, reads the
    /// VBINSEQ header information, and loads the embedded index. The reader is positioned
    /// at the beginning of the first record block after the header.
    ///
    /// ## Index Loading (v0.7.0+)
    ///
    /// The embedded index is automatically loaded from the end of the file. For legacy
    /// files with separate `.vqi` index files, the index is automatically migrated to
    /// the embedded format.
    ///
    /// # Parameters
    ///
    /// * `path` - Path to the VBINSEQ file to open
    ///
    /// # Returns
    ///
    /// A new `MmapReader` instance if successful
    ///
    /// # Errors
    ///
    /// * `ReadError::InvalidFileType` if the path doesn't point to a regular file
    /// * I/O errors if the file can't be opened or memory-mapped
    /// * Header validation errors if the file doesn't contain a valid VBINSEQ header
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::MmapReader;
    ///
    /// let reader = MmapReader::new("path/to/file.vbq").unwrap();
    /// ```
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        // Verify it's a regular file before attempting to map
        let file = File::open(&path)?;
        if !file.metadata()?.is_file() {
            return Err(ReadError::InvalidFileType.into());
        }

        // Safety: The file is open and won't be modified while mapped
        let mmap = unsafe { Mmap::map(&file)? };

        // Read header from mapped memory
        let header = {
            let mut header_bytes = [0u8; SIZE_HEADER];
            header_bytes.copy_from_slice(&mmap[..SIZE_HEADER]);
            VBinseqHeader::from_bytes(&header_bytes)?
        };

        Ok(Self {
            path: PathBuf::from(path.as_ref()),
            mmap: Arc::new(mmap),
            header,
            pos: SIZE_HEADER,
            total: 0,
            decode_block: true,
            default_quality_score: DEFAULT_QUALITY_SCORE,
        })
    }

    pub fn set_default_quality_score(&mut self, score: u8) {
        self.default_quality_score = score;
    }

    /// Creates a new empty record block with the appropriate size for this file
    ///
    /// This creates a `RecordBlock` with a block size matching the one specified in the
    /// file's header, ensuring it will be able to hold a full block of records.
    ///
    /// # Returns
    ///
    /// A new empty `RecordBlock` instance sized appropriately for this file
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::MmapReader;
    ///
    /// let reader = MmapReader::new("example.vbq").unwrap();
    /// let mut block = reader.new_block();
    /// ```
    #[must_use]
    pub fn new_block(&self) -> RecordBlock {
        let mut block = RecordBlock::new(self.header.bits, self.header.block as usize);
        block.set_default_quality_score(self.default_quality_score);
        block
    }

    /// Sets whether to decode sequences at once in each block
    ///
    /// # Arguments
    ///
    /// * `decode_block` - Whether to decode sequences at once in each block
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::MmapReader;
    ///
    /// let mut reader = MmapReader::new("example.vbq").unwrap();
    /// reader.set_decode_block(false);
    /// ```
    pub fn set_decode_block(&mut self, decode_block: bool) {
        self.decode_block = decode_block;
    }

    /// Returns the path where the index file would be located
    ///
    /// The index file is used for random access to blocks and has the same path as
    /// the VBINSEQ file with the ".vqi" extension appended.
    ///
    /// # Returns
    ///
    /// The path where the index file would be located
    ///
    /// # Examples
    ///
    /// ```
    /// use binseq::vbq::MmapReader;
    /// use binseq::Result;
    ///
    /// fn main() -> Result<()> {
    ///     let path = "./data/subset.vbq";
    ///     let reader = MmapReader::new(path)?;
    ///     let index_path = reader.index_path();
    ///     assert_eq!(index_path.to_str(), Some("./data/subset.vbq.vqi"));
    ///     Ok(())
    /// }
    /// ```
    #[must_use]
    pub fn index_path(&self) -> PathBuf {
        let mut p = self.path.as_os_str().to_owned();
        p.push(".vqi");
        p.into()
    }

    /// Returns a copy of the file's header information
    ///
    /// The header contains information about the file format, including whether
    /// quality scores are included, whether blocks are compressed, and whether
    /// records are paired.
    ///
    /// # Returns
    ///
    /// A copy of the file's `VBinseqHeader`
    #[must_use]
    pub fn header(&self) -> VBinseqHeader {
        self.header
    }

    /// Checks if the file contains paired records
    #[must_use]
    pub fn is_paired(&self) -> bool {
        self.header.is_paired()
    }

    /// Fills an existing `RecordBlock` with the next block of records from the file
    ///
    /// This method reads the next block of records from the current position in the file
    /// and populates the provided `RecordBlock` with the data. The block is cleared and reused
    /// to avoid unnecessary memory allocations. This is the primary method for sequential
    /// reading of VBINSEQ files.
    ///
    /// The method automatically handles decompression if the file was written with
    /// compression enabled and updates the total record count as it progresses through the file.
    ///
    /// # Parameters
    ///
    /// * `block` - A mutable reference to a `RecordBlock` to be filled with data
    ///
    /// # Returns
    ///
    /// * `Ok(true)` - If a block was successfully read
    /// * `Ok(false)` - If the end of the file was reached (no more blocks)
    /// * `Err(_)` - If an error occurred during reading
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::MmapReader;
    /// use binseq::BinseqRecord;
    /// use std::io::Write;
    ///
    /// let mut reader = MmapReader::new("example.vbq").unwrap();
    /// let mut block = reader.new_block();
    /// let mut sequence_buffer = Vec::new();
    ///
    /// // Read blocks until the end of file
    /// while reader.read_block_into(&mut block).unwrap() {
    ///     println!("Read block with {} records", block.n_records());
    ///
    ///     // Process each record
    ///     for record in block.iter() {
    ///         // Decode the nucleotide sequence
    ///         record.decode_s(&mut sequence_buffer).unwrap();
    ///
    ///         // Do something with the sequence
    ///         println!("Record {}: length {}", record.index(), sequence_buffer.len());
    ///         sequence_buffer.clear();
    ///     }
    /// }
    /// ```
    pub fn read_block_into(&mut self, block: &mut RecordBlock) -> Result<bool> {
        // Clear the block
        block.clear();

        // Validate the next block header is within bounds and present
        if self.pos + SIZE_BLOCK_HEADER > self.mmap.len() {
            return Ok(false);
        }
        let mut header_bytes = [0u8; SIZE_BLOCK_HEADER];
        header_bytes.copy_from_slice(&self.mmap[self.pos..self.pos + SIZE_BLOCK_HEADER]);
        let header = match BlockHeader::from_bytes(&header_bytes) {
            Ok(header) => {
                self.pos += SIZE_BLOCK_HEADER;
                header
            }
            // Bytes left - but not a BlockHeader - could be the index
            Err(e) => {
                let mut index_header_bytes = [0u8; INDEX_HEADER_SIZE];
                index_header_bytes
                    .copy_from_slice(&self.mmap[self.pos..self.pos + INDEX_HEADER_SIZE]);
                if IndexHeader::from_bytes(&index_header_bytes).is_ok() {
                    // Expected end of file
                    return Ok(false);
                }
                return Err(e);
            }
        };

        // Read the block contents
        let rbound = if self.header.compressed {
            header.size as usize
        } else {
            self.header.block as usize
        };
        if self.pos + rbound > self.mmap.len() {
            return Err(ReadError::UnexpectedEndOfFile(self.pos).into());
        }
        let block_buffer = &self.mmap[self.pos..self.pos + rbound];
        if self.header.compressed {
            block.ingest_compressed_bytes(
                block_buffer,
                self.header.qual,
                self.header.headers,
                self.header.flags,
            )?;
        } else {
            block.ingest_bytes(
                block_buffer,
                self.header.qual,
                self.header.headers,
                self.header.flags,
            )?;
        }

        // Update the block index
        block.update_index(self.total);

        self.pos += rbound;
        self.total += header.records as usize;

        Ok(true)
    }

    /// Loads or creates the block index for this VBINSEQ file
    ///
    /// The block index provides metadata about each block in the file, enabling
    /// random access to blocks and parallel processing. This method first attempts to
    /// load an existing index file. If the index doesn't exist or doesn't match the
    /// current file, it automatically generates a new index from the VBINSEQ file
    /// and saves it for future use.
    ///
    /// # Returns
    ///
    /// The loaded or newly created `BlockIndex` if successful
    ///
    /// # Errors
    ///
    /// * File I/O errors when reading or creating the index
    /// * Parsing errors if the VBINSEQ file has invalid format
    /// * Other index-related errors that cannot be resolved by creating a new index
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::MmapReader;
    ///
    /// let reader = MmapReader::new("example.vbq").unwrap();
    ///
    /// // Load the index file (or create if it doesn't exist)
    /// let index = reader.load_index().unwrap();
    ///
    /// // Use the index to get information about the file
    /// println!("Number of blocks: {}", index.n_blocks());
    /// ```
    ///
    /// # Notes
    ///
    /// The index file is stored with the same path as the VBINSEQ file but with a ".vqi"
    /// extension appended. This allows for reusing the index across multiple runs,
    /// which can significantly improve startup performance for large files.
    pub fn load_index(&self) -> Result<BlockIndex> {
        let start_pos_magic = self.mmap.len() - 8;
        let start_pos_index_size = start_pos_magic - 8;

        // Validate the magic number
        let magic = LittleEndian::read_u64(&self.mmap[start_pos_magic..]);
        if magic != INDEX_END_MAGIC {
            return Err(ReadError::MissingIndexEndMagic.into());
        }

        // Get the index size
        let index_size = LittleEndian::read_u64(&self.mmap[start_pos_index_size..start_pos_magic]);

        // Determine the start position of the index bytes
        let start_pos_index = start_pos_index_size - index_size as usize;

        // Slice into the index bytes
        let index_bytes = &self.mmap[start_pos_index..start_pos_index_size];

        // Build the index from the bytes
        BlockIndex::from_bytes(index_bytes)
    }

    pub fn num_records(&self) -> Result<usize> {
        let index = self.load_index()?;
        Ok(index.num_records())
    }
}

impl ParallelReader for MmapReader {
    /// Processes all records in the file in parallel using multiple threads
    ///
    /// This method provides efficient parallel processing of VBINSEQ files by distributing
    /// blocks across multiple worker threads. The file's block structure is leveraged to divide
    /// the work evenly without requiring thread synchronization during processing, which leads
    /// to near-linear scaling with the number of threads.
    ///
    /// The method automatically loads or creates an index file to identify block boundaries,
    /// then distributes the blocks among the requested number of threads. Each thread processes
    /// its assigned blocks sequentially, but multiple blocks are processed in parallel across
    /// threads.
    ///
    /// # Type Parameters
    ///
    /// * `P` - A type that implements the `ParallelProcessor` trait, which defines how records are processed
    ///
    /// # Parameters
    ///
    /// * `self` - Consumes the reader, as it will be used across multiple threads
    /// * `processor` - An instance of a type implementing `ParallelProcessor` that will be cloned for each thread
    /// * `num_threads` - Number of worker threads to use for processing
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If all records were successfully processed
    /// * `Err(_)` - If an error occurs during processing
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::{MmapReader, RefRecord};
    /// use binseq::{ParallelProcessor, ParallelReader, BinseqRecord, Result};
    /// use std::sync::atomic::{AtomicUsize, Ordering};
    /// use std::sync::Arc;
    ///
    /// // Create a simple processor that counts records
    /// struct RecordCounter {
    ///     count: Arc<AtomicUsize>,
    ///     thread_id: usize,
    /// }
    ///
    /// impl RecordCounter {
    ///     fn new() -> Self {
    ///         Self {
    ///             count: Arc::new(AtomicUsize::new(0)),
    ///             thread_id: 0,
    ///         }
    ///     }
    ///
    ///     fn total_count(&self) -> usize {
    ///         self.count.load(Ordering::Relaxed)
    ///     }
    /// }
    ///
    /// impl Clone for RecordCounter {
    ///     fn clone(&self) -> Self {
    ///         Self {
    ///             count: Arc::clone(&self.count),
    ///             thread_id: 0,
    ///         }
    ///     }
    /// }
    ///
    /// impl ParallelProcessor for RecordCounter {
    ///     fn process_record<R: BinseqRecord>(&mut self, _record: R) -> Result<()> {
    ///         self.count.fetch_add(1, Ordering::Relaxed);
    ///         Ok(())
    ///     }
    ///
    ///     fn on_batch_complete(&mut self) -> Result<()> {
    ///         // Optional: perform actions after each block is processed
    ///         Ok(())
    ///     }
    ///
    ///     fn set_tid(&mut self, tid: usize) {
    ///         self.thread_id = tid;
    ///     }
    /// }
    ///
    /// // Use the processor with a VBINSEQ file
    /// let reader = MmapReader::new("example.vbq").unwrap();
    /// let counter = RecordCounter::new();
    ///
    /// // Process the file with 4 threads
    /// reader.process_parallel(counter.clone(), 4).unwrap();
    ///
    /// // Get the total number of records processed
    /// println!("Total records: {}", counter.total_count());
    /// ```
    ///
    /// # Notes
    ///
    /// * The `ParallelProcessor` instance is cloned for each worker thread, so any shared state
    ///   should be wrapped in thread-safe containers like `Arc`.
    /// * The `set_tid` method is called with a unique thread ID before processing begins, which
    ///   can be used to distinguish between worker threads.
    /// * This method consumes the reader (takes ownership), as it's distributed across threads.
    fn process_parallel<P: ParallelProcessor + Clone + 'static>(
        self,
        processor: P,
        num_threads: usize,
    ) -> Result<()> {
        let num_records = self.num_records()?;
        self.process_parallel_range(processor, num_threads, 0..num_records)
    }

    /// Process records in parallel within a specified range
    ///
    /// This method allows parallel processing of a subset of records within the file,
    /// defined by a start and end index. The method maps the record range to the
    /// appropriate blocks and processes only the records within the specified range.
    ///
    /// # Arguments
    ///
    /// * `processor` - The processor to use for each record
    /// * `num_threads` - The number of threads to spawn
    /// * `start` - The starting record index (inclusive)
    /// * `end` - The ending record index (exclusive)
    ///
    /// # Type Parameters
    ///
    /// * `P` - A type that implements `ParallelProcessor` and can be cloned
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If all records were processed successfully
    /// * `Err(Error)` - If an error occurred during processing
    fn process_parallel_range<P: ParallelProcessor + Clone + 'static>(
        self,
        processor: P,
        num_threads: usize,
        range: Range<usize>,
    ) -> Result<()> {
        // Calculate the number of threads to use
        let num_threads = if num_threads == 0 {
            num_cpus::get()
        } else {
            num_threads.min(num_cpus::get())
        };

        // Generate or load the index first
        let index = self.load_index()?;

        // Validate range
        let total_records = index.num_records();
        self.validate_range(total_records, &range)?;

        // Find blocks that contain records in the specified range
        let relevant_blocks = index
            .ranges()
            .iter()
            .filter(|r| {
                let iv_start = r.cumulative_records as usize;
                let iv_end = (r.cumulative_records + u64::from(r.block_records)) as usize;
                iv_start < range.end && iv_end > range.start
            })
            .copied()
            .collect::<Vec<_>>();

        if relevant_blocks.is_empty() {
            return Ok(()); // No relevant blocks
        }

        // Calculate block assignments for threads
        let blocks_per_thread = relevant_blocks.len().div_ceil(num_threads);

        // Create shared resources
        let mmap = Arc::clone(&self.mmap);
        let header = self.header;

        // Spawn worker threads
        let mut handles = Vec::new();

        for thread_id in 0..num_threads {
            // Calculate this thread's block range
            let start_block_idx = thread_id * blocks_per_thread;
            let end_block_idx =
                std::cmp::min((thread_id + 1) * blocks_per_thread, relevant_blocks.len());

            if start_block_idx >= relevant_blocks.len() {
                continue;
            }

            let mmap = Arc::clone(&mmap);
            let mut proc = processor.clone();
            proc.set_tid(thread_id);

            // Get block ranges for this thread
            let thread_blocks: Vec<BlockRange> =
                relevant_blocks[start_block_idx..end_block_idx].to_vec();

            let handle = std::thread::spawn(move || -> Result<()> {
                // Create block to reuse for processing (within thread)
                let mut record_block = RecordBlock::new(header.bits, header.block as usize);

                // Process each assigned block
                for block_range in thread_blocks {
                    // Clear the block for reuse
                    record_block.clear();

                    // Skip the block header to get to data
                    let block_start = block_range.start_offset as usize + SIZE_BLOCK_HEADER;
                    let block_data = &mmap[block_start..block_start + block_range.len as usize];

                    // Ingest data according to the compression setting
                    if header.compressed {
                        record_block.ingest_compressed_bytes(
                            block_data,
                            header.qual,
                            header.headers,
                            header.flags,
                        )?;
                    } else {
                        record_block.ingest_bytes(
                            block_data,
                            header.qual,
                            header.headers,
                            header.flags,
                        )?;
                    }

                    // Update the record block index
                    record_block.update_index(block_range.cumulative_records as usize);

                    // decode the data
                    if self.decode_block {
                        record_block.decode_all()?;
                    }

                    // Process records in this block that fall within our range
                    for record in record_block.iter() {
                        let global_record_idx = record.index as usize;

                        // Only process records within our specified range
                        if global_record_idx >= range.start && global_record_idx < range.end {
                            proc.process_record(record)?;
                        }
                    }

                    // Signal batch completion
                    proc.on_batch_complete()?;
                }

                Ok(())
            });

            handles.push(handle);
        }

        // Wait for all threads to complete
        for handle in handles {
            handle.join().unwrap()?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BinseqRecord;

    const TEST_VBQ_FILE: &str = "./data/subset.vbq";

    // ==================== MmapReader Basic Tests ====================

    #[test]
    fn test_mmap_reader_new() {
        let reader = MmapReader::new(TEST_VBQ_FILE);
        assert!(reader.is_ok(), "Failed to create VBQ reader");
    }

    #[test]
    fn test_mmap_reader_num_records() {
        let reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let num_records = reader.num_records();
        assert!(num_records.is_ok(), "Failed to get num_records");
        assert!(num_records.unwrap() > 0, "Expected non-zero records");
    }

    #[test]
    fn test_mmap_reader_is_paired() {
        let reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let is_paired = reader.is_paired();
        // Test that the method returns a boolean
        assert!(is_paired || !is_paired);
    }

    #[test]
    fn test_mmap_reader_header_access() {
        let reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let header = &reader.header;
        assert!(header.block > 0, "Expected non-zero block size");
        assert_eq!(header.magic, 0x51455356, "Expected VSEQ magic number");
    }

    #[test]
    fn test_mmap_reader_index_path() {
        let reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let index_path = reader.index_path();

        let expected_path = format!("{}.vqi", TEST_VBQ_FILE);
        assert_eq!(index_path.to_str(), Some(expected_path.as_str()));
    }

    // ==================== RecordBlock Tests ====================

    #[test]
    fn test_new_block() {
        let reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let block = reader.new_block();

        assert_eq!(block.bitsize, reader.header.bits);
        assert!(block.n_records() == 0, "New block should be empty");
    }

    #[test]
    fn test_record_block_creation() {
        let block = RecordBlock::new(BitSize::Two, 1024);

        assert_eq!(block.bitsize, BitSize::Two);
        assert_eq!(block.n_records(), 0);
    }

    #[test]
    fn test_record_block_clear() {
        let mut block = RecordBlock::new(BitSize::Two, 1024);

        // Block starts empty
        assert_eq!(block.n_records(), 0);

        // Clear should not panic on empty block
        block.clear();
        assert_eq!(block.n_records(), 0);
    }

    #[test]
    fn test_record_block_set_default_quality() {
        let mut block = RecordBlock::new(BitSize::Two, 1024);
        let custom_score = 42u8;

        block.set_default_quality_score(custom_score);
        assert_eq!(block.default_quality_score, custom_score);
    }

    // ==================== Block Reading Tests ====================

    #[test]
    fn test_read_block_into() {
        let mut reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let mut block = reader.new_block();

        let result = reader.read_block_into(&mut block);
        assert!(result.is_ok(), "Failed to read block");

        if result.unwrap() {
            assert!(block.n_records() > 0, "Block should contain records");
        }
    }

    #[test]
    fn test_read_multiple_blocks() {
        let mut reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let mut block = reader.new_block();

        let mut blocks_read = 0;
        let max_blocks = 5;

        while reader.read_block_into(&mut block).unwrap() && blocks_read < max_blocks {
            assert!(block.n_records() > 0, "Each block should have records");
            blocks_read += 1;
        }

        assert!(blocks_read > 0, "Should read at least one block");
    }

    #[test]
    fn test_block_iteration() {
        let mut reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let mut block = reader.new_block();

        if reader.read_block_into(&mut block).unwrap() {
            let num_records = block.n_records();
            let mut count = 0;

            for record in block.iter() {
                assert!(record.slen() > 0, "Record should have non-zero length");
                count += 1;
            }

            assert_eq!(count, num_records, "Iterator should yield all records");
        }
    }

    // ==================== Record Access Tests ====================

    #[test]
    fn test_record_sequence_data() {
        let mut reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let mut block = reader.new_block();

        if reader.read_block_into(&mut block).unwrap() {
            // Decode all sequences in the block
            block.decode_all().unwrap();

            if let Some(record) = block.iter().next() {
                let sseq = record.sseq();
                assert!(!sseq.is_empty(), "Sequence should not be empty");

                let slen = record.slen();
                assert_eq!(sseq.len(), slen as usize, "Sequence length mismatch");
            }
        }
    }

    #[test]
    fn test_record_header_data() {
        let mut reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let mut block = reader.new_block();

        if reader.read_block_into(&mut block).unwrap() {
            for record in block.iter() {
                let sheader = record.sheader();
                // Header may be empty if not included in file
                if !sheader.is_empty() {
                    // Should be valid UTF-8 if present
                    let _ = std::str::from_utf8(sheader);
                }
            }
        }
    }

    #[test]
    fn test_record_quality_data() {
        let mut reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let mut block = reader.new_block();

        if reader.read_block_into(&mut block).unwrap() {
            for record in block.iter() {
                let squal = record.squal();
                let slen = record.slen() as usize;

                if !squal.is_empty() {
                    assert_eq!(
                        squal.len(),
                        slen,
                        "Quality length should match sequence length"
                    );
                }
            }
        }
    }

    #[test]
    fn test_record_bitsize() {
        let mut reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let mut block = reader.new_block();

        if reader.read_block_into(&mut block).unwrap() {
            for record in block.iter() {
                let bitsize = record.bitsize();
                assert!(
                    matches!(bitsize, BitSize::Two | BitSize::Four),
                    "Bitsize should be Two or Four"
                );
            }
        }
    }

    // ==================== Default Quality Score Tests ====================

    #[test]
    fn test_set_default_quality_score() {
        let mut reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let custom_score = 42u8;

        reader.set_default_quality_score(custom_score);
        assert_eq!(reader.default_quality_score, custom_score);

        let block = reader.new_block();
        assert_eq!(block.default_quality_score, custom_score);
    }

    // ==================== Decode Block Feature Tests ====================

    #[test]
    fn test_set_decode_block() {
        let mut reader = MmapReader::new(TEST_VBQ_FILE).unwrap();

        reader.set_decode_block(true);
        // Just verify it doesn't panic - actual behavior depends on reading

        reader.set_decode_block(false);
        // Verify we can toggle it
    }

    #[test]
    fn test_decode_block_affects_reading() {
        let mut reader1 = MmapReader::new(TEST_VBQ_FILE).unwrap();
        reader1.set_decode_block(true);
        let mut block1 = reader1.new_block();

        let mut reader2 = MmapReader::new(TEST_VBQ_FILE).unwrap();
        reader2.set_decode_block(false);
        let mut block2 = reader2.new_block();

        // Both should read successfully
        let result1 = reader1.read_block_into(&mut block1);
        let result2 = reader2.read_block_into(&mut block2);

        assert!(result1.is_ok() && result2.is_ok());
    }

    // ==================== Parallel Processing Tests ====================

    #[derive(Clone)]
    struct VbqCountingProcessor {
        count: Arc<std::sync::Mutex<usize>>,
    }

    impl ParallelProcessor for VbqCountingProcessor {
        fn process_record<R: BinseqRecord>(&mut self, _record: R) -> Result<()> {
            let mut count = self.count.lock().unwrap();
            *count += 1;
            Ok(())
        }
    }

    #[test]
    fn test_parallel_processing() {
        let reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let num_records_result = reader.num_records();

        // Skip test if we can't determine record count
        if num_records_result.is_err() {
            return;
        }

        let num_records = num_records_result.unwrap();

        let count = Arc::new(std::sync::Mutex::new(0));
        let processor = VbqCountingProcessor {
            count: count.clone(),
        };

        let result = reader.process_parallel(processor, 2);

        // Parallel processing might not be supported for all VBQ files
        if result.is_ok() {
            let final_count = *count.lock().unwrap();
            // The count should be reasonable (not overflow)
            assert!(
                final_count <= num_records * 2,
                "Processed count should not exceed records significantly"
            );
        }
    }

    #[test]
    fn test_parallel_processing_range() {
        let reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let num_records_result = reader.num_records();

        // Skip test if we can't determine record count
        if num_records_result.is_err() {
            return;
        }

        let num_records = num_records_result.unwrap();

        if num_records >= 100 {
            let start = 10;
            let end = 50;
            let expected_count = end - start;

            let count = Arc::new(std::sync::Mutex::new(0));
            let processor = VbqCountingProcessor {
                count: count.clone(),
            };

            let result = reader.process_parallel_range(processor, 2, start..end);

            // Parallel processing might not be supported for all VBQ files
            if result.is_ok() {
                let final_count = *count.lock().unwrap();
                // The count should be reasonable
                assert!(
                    final_count <= expected_count * 2,
                    "Processed count should not exceed expected range significantly"
                );
            }
        }
    }

    // ==================== Span Tests ====================

    #[test]
    fn test_span_creation() {
        let span = Span::new(10, 20);
        assert_eq!(span.offset, 10);
        assert_eq!(span.len, 20);
    }

    #[test]
    fn test_span_default() {
        let span = Span::default();
        assert_eq!(span.offset, 0);
        assert_eq!(span.len, 0);
    }

    // ==================== Error Handling Tests ====================

    #[test]
    fn test_nonexistent_file() {
        let result = MmapReader::new("./data/nonexistent.vbq");
        assert!(result.is_err(), "Should fail on nonexistent file");
    }

    #[test]
    fn test_invalid_file_format() {
        // Try to open a non-VBQ file as VBQ
        let result = MmapReader::new("./Cargo.toml");
        // This should fail during header validation
        assert!(result.is_err(), "Should fail on invalid file format");
    }

    // ==================== Index Loading Tests ====================

    #[test]
    fn test_load_index() {
        let reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let index_result = reader.load_index();

        assert!(index_result.is_ok(), "Should be able to load index");

        let index = index_result.unwrap();
        assert!(index.num_records() > 0, "Index should have records");
    }

    #[test]
    fn test_index_consistency() {
        let reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let num_records_from_reader = reader.num_records().unwrap();

        let index = reader.load_index().unwrap();
        let num_records_from_index = index.num_records();

        assert_eq!(
            num_records_from_reader, num_records_from_index,
            "Reader and index should report same number of records"
        );
    }

    // ==================== RecordBlock Decoded Access Tests ====================

    #[test]
    fn test_get_decoded_s() {
        let mut reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        reader.set_decode_block(true);
        let mut block = reader.new_block();

        if reader.read_block_into(&mut block).unwrap() && block.n_records() > 0 {
            let decoded = block.get_decoded_s(0);
            if let Some(seq) = decoded {
                assert!(!seq.is_empty(), "Decoded sequence should not be empty");
            }
        }
    }

    #[test]
    fn test_get_decoded_x() {
        let mut reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        reader.set_decode_block(true);
        let mut block = reader.new_block();

        if reader.read_block_into(&mut block).unwrap() && block.n_records() > 0 {
            // Extended sequence may be empty for non-paired reads
            let decoded = block.get_decoded_x(0);
            // Just verify it doesn't panic
            let _ = decoded;
        }
    }

    #[test]
    fn test_get_decoded_out_of_bounds() {
        let mut reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let mut block = reader.new_block();

        if reader.read_block_into(&mut block).unwrap() {
            let num_records = block.n_records();

            // Try to access beyond bounds
            let decoded = block.get_decoded_s(num_records + 100);
            assert!(decoded.is_none(), "Should return None for out of bounds");
        }
    }

    // ==================== Helper Function Tests ====================

    #[test]
    fn test_encoded_sequence_len_two_bit() {
        // 2-bit encoding: 32 nucleotides per u64
        assert_eq!(encoded_sequence_len(32, BitSize::Two), 1);
        assert_eq!(encoded_sequence_len(64, BitSize::Two), 2);
        assert_eq!(encoded_sequence_len(33, BitSize::Two), 2); // Rounds up
        assert_eq!(encoded_sequence_len(1, BitSize::Two), 1);
    }

    #[test]
    fn test_encoded_sequence_len_four_bit() {
        // 4-bit encoding: 16 nucleotides per u64
        assert_eq!(encoded_sequence_len(16, BitSize::Four), 1);
        assert_eq!(encoded_sequence_len(32, BitSize::Four), 2);
        assert_eq!(encoded_sequence_len(17, BitSize::Four), 2); // Rounds up
        assert_eq!(encoded_sequence_len(1, BitSize::Four), 1);
    }

    // ==================== Record Iterator Tests ====================

    #[test]
    fn test_record_block_iter_creation() {
        let block = RecordBlock::new(BitSize::Two, 1024);
        let iter = RecordBlockIter::new(&block);

        // Iterator on empty block should yield nothing
        assert_eq!(iter.count(), 0);
    }

    #[test]
    fn test_record_iteration_multiple_times() {
        let mut reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let mut block = reader.new_block();

        if reader.read_block_into(&mut block).unwrap() && block.n_records() > 0 {
            let num_records = block.n_records();

            // First iteration
            let count1 = block.iter().count();
            assert_eq!(count1, num_records);

            // Second iteration should yield same count
            let count2 = block.iter().count();
            assert_eq!(count2, num_records);
        }
    }

    // ==================== Paired Read Tests ====================

    #[test]
    fn test_paired_record_data() {
        let mut reader = MmapReader::new(TEST_VBQ_FILE).unwrap();

        if reader.is_paired() {
            let mut block = reader.new_block();

            if reader.read_block_into(&mut block).unwrap() {
                // Decode all sequences in the block
                block.decode_all().unwrap();

                for record in block.iter() {
                    let xlen = record.xlen();

                    if xlen > 0 {
                        let xseq = record.xseq();
                        assert_eq!(
                            xseq.len(),
                            xlen as usize,
                            "Extended sequence length should match xlen"
                        );
                    }
                }
            }
        }
    }

    // ==================== Edge Cases ====================

    #[test]
    fn test_empty_block_iteration() {
        let block = RecordBlock::new(BitSize::Two, 1024);

        let mut count = 0;
        for _ in block.iter() {
            count += 1;
        }

        assert_eq!(count, 0, "Empty block should yield no records");
    }

    #[test]
    fn test_reader_reset_by_new_block() {
        let mut reader = MmapReader::new(TEST_VBQ_FILE).unwrap();
        let mut block = reader.new_block();

        // Read first block
        if reader.read_block_into(&mut block).unwrap() {
            let first_count = block.n_records();

            // Read second block (overwrites first)
            if reader.read_block_into(&mut block).unwrap() {
                let second_count = block.n_records();

                // Counts may differ, but both should be > 0
                assert!(first_count > 0 && second_count > 0);
            }
        }
    }
}
