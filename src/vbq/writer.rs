//! Writer implementation for VBINSEQ files
//!
//! This module provides functionality for writing sequence data to VBINSEQ files,
//! including support for compression, quality scores, paired-end reads, and sequence headers.
//!
//! The VBINSEQ writer implements a block-based approach where records are packed
//! into fixed-size blocks. Each block has a header containing metadata about the
//! records it contains. Blocks may be optionally compressed using zstd compression.
//!
//! ## Format Changes (v0.7.0+)
//!
//! - **Embedded Index**: Writers now automatically embed an index at the end of the file
//! - **Headers Support**: Optional sequence headers/identifiers can be written with each record
//! - **Multi-bit Encoding**: Support for 2-bit and 4-bit nucleotide encodings
//! - **Extended Capacity**: u64 indexing supports more than 4 billion records
//!
//! ## File Structure Written
//!
//! ```text
//! [File Header][Data Blocks][Compressed Index][Index Size][Index End Magic]
//! ```
//!
//! The writer automatically:
//! 1. Writes the file header
//! 2. Writes data blocks as records are added
//! 3. Builds an index during writing
//! 4. On `finish()`, compresses and embeds the index at the end of the file
//!
//! # Example
//!
//! ```rust,no_run
//! use binseq::vbq::{VBinseqWriterBuilder, VBinseqHeaderBuilder};
//! use std::fs::File;
//!
//! // Create a VBINSEQ file writer with headers and compression
//! let file = File::create("example.vbq").unwrap();
//! let header = VBinseqHeaderBuilder::new()
//!     .block(128 * 1024)
//!     .qual(true)
//!     .compressed(true)
//!     .headers(true)
//!     .flags(true)
//!     .build();
//!
//! let mut writer = VBinseqWriterBuilder::default()
//!     .header(header)
//!     .build(file)
//!     .unwrap();
//!
//! // Write a nucleotide sequence with quality scores and header
//! let sequence = b"ACGTACGTACGT";
//! let quality = b"IIIIIIIIIIII";
//! let header_str = b"sequence_001";
//! writer.write_record(Some(0), Some(header_str), sequence, Some(quality)).unwrap();
//!
//! // Must call finish() to write the embedded index
//! writer.finish().unwrap();
//! ```

use std::io::Write;

use bitnuc::BitSize;
use byteorder::{LittleEndian, WriteBytesExt};
use rand::SeedableRng;
use rand::rngs::SmallRng;
use zstd::stream::copy_encode;

use super::header::{BlockHeader, VBinseqHeader};
use crate::error::{Result, WriteError};
use crate::policy::{Policy, RNG_SEED};
use crate::vbq::header::{SIZE_BLOCK_HEADER, SIZE_HEADER};
use crate::vbq::index::{INDEX_END_MAGIC, IndexHeader};
use crate::vbq::{BlockIndex, BlockRange};

/// Calculates the storage size in bytes required for a record without quality scores
///
/// This function calculates the total size needed to store a record in the VBINSEQ format,
/// including the flag, sequence lengths, and the encoded sequence data. The formula
/// used is: `S = w(Cs + Cx + 3)` where:
///
/// - `w`: Word size (8 bytes)
/// - `Cs`: Chunk size of the primary sequence in 64-bit words
/// - `Cx`: Chunk size of the extended sequence in 64-bit words (for paired-end reads)
/// - `3`: Additional words for flag, primary length, and extended length
///
/// # Parameters
///
/// * `schunk` - Number of 64-bit words needed for the primary sequence
/// * `xchunk` - Number of 64-bit words needed for the extended sequence (0 for single-end)
///
/// # Returns
///
/// The total size in bytes needed to store the record
pub fn record_byte_size(schunk: usize, xchunk: usize, has_flags: bool) -> usize {
    8 * (schunk + xchunk + if has_flags { 3 } else { 2 })
}

fn record_byte_size_quality_header(
    schunk: usize,
    xchunk: usize,
    squal: usize,
    xqual: usize,
    sheader: usize,
    xheader: usize,
    has_flags: bool,
) -> usize {
    // counting the header length bytes (u64)
    let bytes_sheader = if sheader > 0 { sheader + 8 } else { 0 };
    let bytes_xheader = if xheader > 0 { xheader + 8 } else { 0 };
    record_byte_size(schunk, xchunk, has_flags) + squal + xqual + bytes_sheader + bytes_xheader
}

/// A builder for creating configured `VBinseqWriter` instances
///
/// This builder provides a fluent interface for configuring and creating a
/// `VBinseqWriter` with customized settings. It allows specifying the file header,
/// encoding policy, and whether to operate in headless mode.
///
/// # Examples
///
/// ```rust,no_run
/// use binseq::vbq::{VBinseqWriterBuilder, VBinseqHeaderBuilder};
/// use binseq::Policy;
/// use std::fs::File;
///
/// // Create a writer with custom settings
/// let file = File::create("example.vbq").unwrap();
/// let mut writer = VBinseqWriterBuilder::default()
///     .header(VBinseqHeaderBuilder::new()
///         .block(65536)
///         .qual(true)
///         .compressed(true)
///         .build())
///     .policy(Policy::IgnoreSequence)
///     .build(file)
///     .unwrap();
///
/// // Use the writer...
/// ```
#[derive(Default)]
pub struct VBinseqWriterBuilder {
    /// Header of the file
    header: Option<VBinseqHeader>,
    /// Optional policy for encoding
    policy: Option<Policy>,
    /// Optional headless mode (used in parallel writing)
    headless: Option<bool>,
}
impl VBinseqWriterBuilder {
    /// Sets the header for the VBINSEQ file
    ///
    /// The header defines the file format parameters such as block size, whether
    /// the file contains quality scores, paired-end reads, and compression settings.
    ///
    /// # Parameters
    ///
    /// * `header` - The `VBinseqHeader` to use for the file
    ///
    /// # Returns
    ///
    /// The builder with the header configured
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::{VBinseqWriterBuilder, VBinseqHeaderBuilder};
    ///
    /// // Create a header with 64KB blocks and quality scores
    /// let header = VBinseqHeaderBuilder::new()
    ///     .block(65536)
    ///     .qual(true)
    ///     .paired(true)
    ///     .compressed(true)
    ///     .build();
    ///
    /// let builder = VBinseqWriterBuilder::default().header(header);
    /// ```
    #[must_use]
    pub fn header(mut self, header: VBinseqHeader) -> Self {
        self.header = Some(header);
        self
    }

    /// Sets the encoding policy for nucleotide sequences
    ///
    /// The policy determines how sequences are encoded into the binary format.
    /// Different policies offer trade-offs between compression ratio and compatibility
    /// with different types of sequence data.
    ///
    /// # Parameters
    ///
    /// * `policy` - The encoding policy to use
    ///
    /// # Returns
    ///
    /// The builder with the encoding policy configured
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::{VBinseqWriterBuilder};
    /// use binseq::Policy;
    ///
    /// let builder = VBinseqWriterBuilder::default().policy(Policy::IgnoreSequence);
    /// ```
    #[must_use]
    pub fn policy(mut self, policy: Policy) -> Self {
        self.policy = Some(policy);
        self
    }

    /// Sets whether to operate in headless mode
    ///
    /// In headless mode, the writer does not write a file header. This is useful
    /// when creating part of a file that will be merged with other parts later,
    /// such as in parallel writing scenarios.
    ///
    /// # Parameters
    ///
    /// * `headless` - Whether to operate in headless mode
    ///
    /// # Returns
    ///
    /// The builder with the headless mode configured
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::VBinseqWriterBuilder;
    ///
    /// // Create a headless writer for parallel writing
    /// let builder = VBinseqWriterBuilder::default().headless(true);
    /// ```
    #[must_use]
    pub fn headless(mut self, headless: bool) -> Self {
        self.headless = Some(headless);
        self
    }

    /// Builds a `VBinseqWriter` with the configured settings
    ///
    /// This finalizes the builder and creates a new `VBinseqWriter` instance using
    /// the provided writer and the configured settings. If any settings were not
    /// explicitly set, default values will be used.
    ///
    /// # Parameters
    ///
    /// * `inner` - The underlying writer where data will be written
    ///
    /// # Returns
    ///
    /// * `Ok(VBinseqWriter)` - A configured `VBinseqWriter` ready for use
    /// * `Err(_)` - If an error occurred while initializing the writer
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::VBinseqWriterBuilder;
    /// use std::fs::File;
    ///
    /// let file = File::create("example.vbq").unwrap();
    /// let mut writer = VBinseqWriterBuilder::default()
    ///     .build(file)
    ///     .unwrap();
    /// ```
    pub fn build<W: Write>(self, inner: W) -> Result<VBinseqWriter<W>> {
        VBinseqWriter::new(
            inner,
            self.header.unwrap_or_default(),
            self.policy.unwrap_or_default(),
            self.headless.unwrap_or(false),
        )
    }
}

/// Writer for VBINSEQ format files
///
/// The `VBinseqWriter` handles writing nucleotide sequence data to VBINSEQ files in a
/// block-based format. It manages the file structure, compression settings, and ensures
/// data is properly encoded and organized.
///
/// ## File Structure
///
/// A VBINSEQ file consists of:
/// 1. A file header that defines parameters like block size and compression settings
/// 2. A series of blocks, each with:
///    - A block header with metadata (e.g., record count)
///    - A collection of encoded records
///
/// Each block is filled with records until either the block is full or no more complete
/// records can fit. The writer automatically handles block boundaries and creates new
/// blocks as needed.
///
/// ## Usage
///
/// The writer supports multiple formats:
/// - Single-end sequences with or without quality scores
/// - Paired-end sequences with or without quality scores
///
/// It's recommended to use the `VBinseqWriterBuilder` to create and configure a writer
/// instance with the appropriate settings.
///
/// ```rust,no_run
/// use binseq::vbq::{VBinseqWriterBuilder, VBinseqHeader};
/// use std::fs::File;
///
/// // Create a writer for single-end reads
/// let file = File::create("example.vbq").unwrap();
/// let mut writer = VBinseqWriterBuilder::default()
///     .header(VBinseqHeader::default())
///     .build(file)
///     .unwrap();
///
/// // Write a sequence
/// let sequence = b"ACGTACGTACGT";
/// writer.write_record(None, None, sequence, None).unwrap();
///
/// // Writer automatically flushes when dropped
/// ```
#[derive(Clone)]
pub struct VBinseqWriter<W: Write> {
    /// Inner Writer
    inner: W,

    /// Header of the file
    header: VBinseqHeader,

    /// Encoder for nucleotide sequences
    encoder: Encoder,

    /// Pre-initialized writer for compressed blocks
    cblock: BlockWriter,

    /// Growable buffer for the block ranges found
    ranges: Vec<BlockRange>,

    /// Total bytes written to this writer
    bytes_written: usize,

    /// Total records written to this writer
    records_written: usize,

    /// Determines if index is already written
    index_written: bool,
}
impl<W: Write> VBinseqWriter<W> {
    pub fn new(inner: W, header: VBinseqHeader, policy: Policy, headless: bool) -> Result<Self> {
        let mut wtr = Self {
            inner,
            header,
            encoder: Encoder::with_policy(header.bits, policy),
            cblock: BlockWriter::new(header.block as usize, header.compressed, header.flags),
            ranges: Vec::new(),
            bytes_written: 0,
            records_written: 0,
            index_written: false,
        };
        if !headless {
            wtr.init()?;
        }
        Ok(wtr)
    }

    /// Initializes the writer by writing the file header
    ///
    /// This method is called automatically during creation unless headless mode is enabled.
    /// It writes the `VBinseqHeader` to the underlying writer.
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If the header was successfully written
    /// * `Err(_)` - If an error occurred during writing
    fn init(&mut self) -> Result<()> {
        self.header.write_bytes(&mut self.inner)?;
        self.bytes_written += SIZE_HEADER;
        Ok(())
    }

    /// Checks if the writer is configured for paired-end reads
    ///
    /// This method returns whether the writer expects paired-end reads based on the
    /// header settings. If true, you should use `write_paired_nucleotides` instead of
    /// `write_nucleotides` to write sequences.
    ///
    /// # Returns
    ///
    /// `true` if the writer is configured for paired-end reads, `false` otherwise
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::{VBinseqWriterBuilder, VBinseqHeader};
    /// use std::fs::File;
    ///
    /// // Create a header for paired-end reads
    /// let mut header = VBinseqHeader::default();
    /// header.paired = true;
    ///
    /// let file = File::create("paired_reads.vbq").unwrap();
    /// let writer = VBinseqWriterBuilder::default()
    ///     .header(header)
    ///     .build(file)
    ///     .unwrap();
    ///
    /// assert!(writer.is_paired());
    /// ```
    pub fn is_paired(&self) -> bool {
        self.header.paired
    }

    /// Checks if the writer is configured for quality scores
    ///
    /// This method returns whether the writer expects quality scores based on the
    /// header settings. If true, you should use methods that include quality scores
    /// (`write_nucleotides_with_quality` or `write_paired_nucleotides_with_quality`)
    /// to write sequences.
    ///
    /// # Returns
    ///
    /// `true` if the writer is configured for quality scores, `false` otherwise
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::{VBinseqWriterBuilder, VBinseqHeader};
    /// use std::fs::File;
    ///
    /// // Create a header for sequences with quality scores
    /// let mut header = VBinseqHeader::default();
    /// header.qual = true;
    ///
    /// let file = File::create("reads_with_quality.vbq").unwrap();
    /// let writer = VBinseqWriterBuilder::default()
    ///     .header(header)
    ///     .build(file)
    ///     .unwrap();
    ///
    /// assert!(writer.has_quality());
    /// ```
    pub fn has_quality(&self) -> bool {
        self.header.qual
    }

    pub fn has_headers(&self) -> bool {
        self.header.headers
    }

    pub fn write_record(
        &mut self,
        flag: Option<u64>,
        header: Option<&[u8]>,
        sequence: &[u8],
        quality: Option<&[u8]>,
    ) -> Result<bool> {
        if self.is_paired() {
            return Err(WriteError::PairedFlagSet.into());
        }

        // ignore the header if not set
        let header = if header.is_none() && self.header.headers {
            return Err(WriteError::HeaderFlagSet.into());
        } else if header.is_some() && !self.header.headers {
            None
        } else {
            header
        };

        // ignore the quality if not set
        let quality = if quality.is_none() && self.header.qual {
            return Err(WriteError::QualityFlagSet.into());
        } else if quality.is_some() && !self.header.qual {
            None
        } else {
            quality
        };

        // encode the sequence
        if let Some(sbuffer) = self.encoder.encode_single(sequence)? {
            let record_size = record_byte_size_quality_header(
                sbuffer.len(),
                0,
                quality.map_or(0, <[u8]>::len),
                0,
                header.map_or(0, <[u8]>::len),
                0,
                self.header.flags,
            );
            if self.cblock.exceeds_block_size(record_size)? {
                impl_flush_block(
                    &mut self.inner,
                    &mut self.cblock,
                    &mut self.ranges,
                    &mut self.bytes_written,
                    &mut self.records_written,
                )?;
            }

            // Write the flag, length, and sequence to the block
            self.cblock.write_record(
                flag,
                sequence.len() as u64,
                0,
                sbuffer,
                quality,
                header,
                None,
                None,
                None,
            )?;

            // Return true if the sequence was successfully written
            Ok(true)
        } else {
            // Silently ignore sequences that fail encoding
            Ok(false)
        }
    }

    pub fn write_paired_record(
        &mut self,
        flag: Option<u64>,
        s_header: Option<&[u8]>,
        s_sequence: &[u8],
        s_qual: Option<&[u8]>,
        x_header: Option<&[u8]>,
        x_sequence: &[u8],
        x_qual: Option<&[u8]>,
    ) -> Result<bool> {
        if !self.is_paired() {
            return Err(WriteError::PairedFlagNotSet.into());
        }

        let s_header = if s_header.is_none() && self.header.headers {
            return Err(WriteError::HeaderFlagSet.into());
        } else if s_header.is_some() && !self.header.headers {
            None
        } else {
            s_header
        };
        let x_header = if x_header.is_none() && self.header.headers {
            return Err(WriteError::HeaderFlagSet.into());
        } else if x_header.is_some() && !self.header.headers {
            None
        } else {
            x_header
        };

        let s_qual = if s_qual.is_none() && self.header.qual {
            return Err(WriteError::QualityFlagSet.into());
        } else if s_qual.is_some() && !self.header.qual {
            None
        } else {
            s_qual
        };

        let x_qual = if x_qual.is_none() && self.header.qual {
            return Err(WriteError::QualityFlagSet.into());
        } else if x_qual.is_some() && !self.header.qual {
            None
        } else {
            x_qual
        };

        // encode the sequences
        if let Some((sbuffer, xbuffer)) = self.encoder.encode_paired(s_sequence, x_sequence)? {
            // Check if the current block can handle the next record
            let record_size = record_byte_size_quality_header(
                sbuffer.len(),
                xbuffer.len(),
                s_qual.map_or(0, <[u8]>::len),
                x_qual.map_or(0, <[u8]>::len),
                s_header.map_or(0, <[u8]>::len),
                x_header.map_or(0, <[u8]>::len),
                self.header.flags,
            );
            if self.cblock.exceeds_block_size(record_size)? {
                impl_flush_block(
                    &mut self.inner,
                    &mut self.cblock,
                    &mut self.ranges,
                    &mut self.bytes_written,
                    &mut self.records_written,
                )?;
            }

            // Write the flag, length, sequence, and quality scores to the block
            self.cblock.write_record(
                flag,
                s_sequence.len() as u64,
                x_sequence.len() as u64,
                sbuffer,
                s_qual,
                s_header,
                Some(xbuffer),
                x_qual,
                x_header,
            )?;

            // Return true if the record was successfully written
            Ok(true)
        } else {
            // Return false if the record was not successfully written
            Ok(false)
        }
    }

    /// Finishes writing and flushes all data to the underlying writer
    ///
    /// This method should be called when you're done writing to ensure all data
    /// is properly flushed to the underlying writer. It's automatically called
    /// when the writer is dropped, but calling it explicitly allows you to handle
    /// any errors that might occur during flushing.
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If all data was successfully flushed
    /// * `Err(_)` - If an error occurred during flushing
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::{VBinseqWriterBuilder, VBinseqHeader};
    /// use std::fs::File;
    ///
    /// let file = File::create("example.vbq").unwrap();
    /// let mut writer = VBinseqWriterBuilder::default()
    ///     .build(file)
    ///     .unwrap();
    ///
    /// // Write some sequences...
    /// let sequence = b"ACGTACGTACGT";
    /// writer.write_record(None, None, sequence, None).unwrap();
    ///
    /// // Manually finish and check for errors
    /// if let Err(e) = writer.finish() {
    ///     eprintln!("Error flushing data: {}", e);
    /// }
    /// ```
    pub fn finish(&mut self) -> Result<()> {
        impl_flush_block(
            &mut self.inner,
            &mut self.cblock,
            &mut self.ranges,
            &mut self.bytes_written,
            &mut self.records_written,
        )?;
        self.inner.flush()?;

        if !self.index_written {
            self.write_index()?;
            self.index_written = true;
        }
        Ok(())
    }

    /// Provides a mutable reference to the inner writer
    fn by_ref(&mut self) -> &mut W {
        self.inner.by_ref()
    }

    /// Provides a mutable reference to the `BlockWriter`
    fn cblock_mut(&mut self) -> &mut BlockWriter {
        &mut self.cblock
    }

    /// Ingests data from another `VBinseqWriter` that uses a `Vec<u8>` as its inner writer
    ///
    /// This method is particularly useful for parallel processing, where multiple writers
    /// might be writing to memory buffers and need to be combined into a single file. It
    /// transfers all complete blocks and any partial blocks from the other writer into this one.
    ///
    /// The method clears the other writer's buffer after ingestion, allowing it to be reused.
    ///
    /// # Parameters
    ///
    /// * `other` - Another `VBinseqWriter` whose inner writer is a `Vec<u8>`
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If ingestion was successful
    /// * `Err(_)` - If an error occurred during ingestion or if the headers are incompatible
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The headers of the two writers are not compatible (`WriteError::IncompatibleHeaders`)
    /// - An I/O error occurred during data transfer
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::{VBinseqWriterBuilder, VBinseqHeader};
    /// use std::fs::File;
    ///
    /// // Create a file writer
    /// let file = File::create("combined.vbq").unwrap();
    /// let mut file_writer = VBinseqWriterBuilder::default()
    ///     .build(file)
    ///     .unwrap();
    ///
    /// // Create a memory writer
    /// let mut mem_writer = VBinseqWriterBuilder::default()
    ///     .build(Vec::new())
    ///     .unwrap();
    ///
    /// // Write some data to the memory writer
    /// mem_writer.write_record(None, None, b"ACGTACGT", None).unwrap();
    ///
    /// // Ingest data from memory writer into file writer
    /// file_writer.ingest(&mut mem_writer).unwrap();
    /// ```
    pub fn ingest(&mut self, other: &mut VBinseqWriter<Vec<u8>>) -> Result<()> {
        if self.header != other.header {
            return Err(WriteError::IncompatibleHeaders(self.header, other.header).into());
        }

        // Write complete blocks from other directly
        // and clear the other (mimics reading)
        {
            self.inner.write_all(other.by_ref())?;
            other.by_ref().clear();
        }

        // Pull the ranges from the other writer and update their statistics
        {
            for range in other.ranges.drain(..) {
                // Build the updated range with main-file specific information
                let updated_range = BlockRange::new(
                    self.bytes_written as u64, // Current position in main file
                    range.len,
                    range.block_records,
                    self.records_written as u64, // Current number of records written in main file
                );

                self.ranges.push(updated_range);

                // Update counters incrementally for each range
                self.bytes_written += (range.len + SIZE_BLOCK_HEADER as u64) as usize;
                self.records_written += range.block_records as usize;
            }

            // reset the other writer
            other.bytes_written = 0;
            other.records_written = 0;
        }

        // Ingest incomplete block from other
        {
            let header = self.cblock.ingest(other.cblock_mut(), &mut self.inner)?;
            if !header.is_empty() {
                let range = BlockRange::new(
                    self.bytes_written as u64,
                    header.size,
                    header.records,
                    self.records_written as u64,
                );
                self.ranges.push(range);
                self.bytes_written += header.size_with_header();
                self.records_written += header.records as usize;
            }
        }
        Ok(())
    }

    pub fn write_index(&mut self) -> Result<()> {
        // Build the index
        let index_header = IndexHeader::new(self.bytes_written as u64);
        let block_index = BlockIndex {
            header: index_header,
            ranges: self.ranges.clone(),
        };

        // Write the index to a temporary buffer
        let mut buffer = Vec::new();
        block_index.write_bytes(&mut buffer)?;

        // Determine the number of bytes written to the buffer
        let n_bytes = buffer.len() as u64;

        // Write the index to the underlying writer
        self.inner.write_all(&buffer)?;

        // Write the number of bytes written to the index
        self.inner.write_u64::<LittleEndian>(n_bytes)?;

        // Write the index footer magic
        self.inner.write_u64::<LittleEndian>(INDEX_END_MAGIC)?;

        Ok(())
    }
}

fn impl_flush_block<W: Write>(
    writer: &mut W,
    cblock: &mut BlockWriter,
    ranges: &mut Vec<BlockRange>,
    bytes_written: &mut usize,
    records_written: &mut usize,
) -> Result<()> {
    let block_header = cblock.flush(writer)?;
    let range = BlockRange::new(
        *bytes_written as u64,
        block_header.size,
        block_header.records,
        *records_written as u64,
    );
    ranges.push(range);
    *bytes_written += block_header.size_with_header();
    *records_written += block_header.records as usize;
    Ok(())
}

impl<W: Write> Drop for VBinseqWriter<W> {
    fn drop(&mut self) {
        self.finish()
            .expect("VBinseqWriter: Failed to finish writing");
    }
}

#[derive(Clone)]
struct BlockWriter {
    /// Current position in the block
    pos: usize,
    /// Tracks all record start positions in the block
    starts: Vec<usize>,
    /// Virtual block size
    block_size: usize,
    /// Compression level
    level: i32,
    /// Uncompressed buffer
    ubuf: Vec<u8>,
    /// Compressed buffer
    zbuf: Vec<u8>,
    /// Reusable padding buffer
    padding: Vec<u8>,
    /// Compression flag
    /// If false, the block is written uncompressed
    compress: bool,
    /// Has flags
    has_flags: bool,
}
impl BlockWriter {
    fn new(block_size: usize, compress: bool, has_flags: bool) -> Self {
        Self {
            pos: 0,
            starts: Vec::default(),
            block_size,
            level: 3,
            ubuf: Vec::with_capacity(block_size),
            zbuf: Vec::with_capacity(block_size),
            padding: vec![0; block_size],
            compress,
            has_flags,
        }
    }

    fn exceeds_block_size(&self, record_size: usize) -> Result<bool> {
        if record_size > self.block_size {
            return Err(WriteError::RecordSizeExceedsMaximumBlockSize(
                record_size,
                self.block_size,
            )
            .into());
        }
        Ok(self.pos + record_size > self.block_size)
    }

    #[allow(clippy::too_many_arguments)]
    fn write_record(
        &mut self,
        flag: Option<u64>,
        slen: u64,
        xlen: u64,
        sbuf: &[u64],
        squal: Option<&[u8]>,
        sheader: Option<&[u8]>,
        xbuf: Option<&[u64]>,
        xqual: Option<&[u8]>,
        xheader: Option<&[u8]>,
    ) -> Result<()> {
        // Tracks the record start position
        self.starts.push(self.pos);

        // Write the flag
        if self.has_flags {
            self.write_flag(flag.unwrap_or(0))?;
        }

        // Write the lengths
        self.write_length(slen)?;
        self.write_length(xlen)?;

        // Write the primary sequence and optional quality
        self.write_buffer(sbuf)?;
        if let Some(qual) = squal {
            self.write_u8buf(qual)?;
        }
        if let Some(sheader) = sheader {
            self.write_length(sheader.len() as u64)?;
            self.write_u8buf(sheader)?;
        }

        // Write the optional extended sequence and optional quality
        if let Some(xbuf) = xbuf {
            self.write_buffer(xbuf)?;
        }
        if let Some(qual) = xqual {
            self.write_u8buf(qual)?;
        }
        if let Some(xheader) = xheader {
            self.write_length(xheader.len() as u64)?;
            self.write_u8buf(xheader)?;
        }

        Ok(())
    }

    fn write_flag(&mut self, flag: u64) -> Result<()> {
        self.ubuf.write_u64::<LittleEndian>(flag)?;
        self.pos += 8;
        Ok(())
    }

    fn write_length(&mut self, length: u64) -> Result<()> {
        self.ubuf.write_u64::<LittleEndian>(length)?;
        self.pos += 8;
        Ok(())
    }

    fn write_buffer(&mut self, ebuf: &[u64]) -> Result<()> {
        ebuf.iter()
            .try_for_each(|&x| self.ubuf.write_u64::<LittleEndian>(x))?;
        self.pos += 8 * ebuf.len();
        Ok(())
    }

    fn write_u8buf(&mut self, buf: &[u8]) -> Result<()> {
        self.ubuf.write_all(buf)?;
        self.pos += buf.len();
        Ok(())
    }

    fn flush_compressed<W: Write>(&mut self, inner: &mut W) -> Result<BlockHeader> {
        // Encode the block
        copy_encode(self.ubuf.as_slice(), &mut self.zbuf, self.level)?;

        // Build a block header (this is variably sized in the compressed case)
        let header = BlockHeader::new(self.zbuf.len() as u64, self.starts.len() as u32);

        // Write the block header and compressed block
        header.write_bytes(inner)?;
        inner.write_all(&self.zbuf)?;

        Ok(header)
    }

    fn flush_uncompressed<W: Write>(&mut self, inner: &mut W) -> Result<BlockHeader> {
        // Build a block header (this is static in size in the uncompressed case)
        let header = BlockHeader::new(self.block_size as u64, self.starts.len() as u32);

        // Write the block header and uncompressed block
        header.write_bytes(inner)?;
        inner.write_all(&self.ubuf)?;

        Ok(header)
    }

    fn flush<W: Write>(&mut self, inner: &mut W) -> Result<BlockHeader> {
        // Skip if the block is empty
        if self.pos == 0 {
            return Ok(BlockHeader::empty());
        }

        // Finish out the block with padding
        let bytes_to_next_start = self.block_size - self.pos;
        self.ubuf.write_all(&self.padding[..bytes_to_next_start])?;

        // Flush the block (implemented differently based on compression)
        let header = if self.compress {
            self.flush_compressed(inner)
        } else {
            self.flush_uncompressed(inner)
        }?;

        // Reset the position and buffers
        self.clear();

        Ok(header)
    }

    fn clear(&mut self) {
        self.pos = 0;
        self.starts.clear();
        self.ubuf.clear();
        self.zbuf.clear();
    }

    /// Ingests *all* bytes from another `BlockWriter`.
    ///
    /// Because both block sizes should be equivalent the process should take
    /// at most two steps.
    ///
    /// I.e. the bytes can either all fit directly into self.ubuf or an intermediate
    /// flush step is required.
    fn ingest<W: Write>(&mut self, other: &mut Self, inner: &mut W) -> Result<BlockHeader> {
        if self.block_size != other.block_size {
            return Err(
                WriteError::IncompatibleBlockSizes(self.block_size, other.block_size).into(),
            );
        }
        // Number of available bytes in buffer (self)
        let remaining = self.block_size - self.pos;

        // Quick ingestion (take all without flush)
        if other.pos <= remaining {
            self.ingest_all(other)?;
            Ok(BlockHeader::empty())
        } else {
            self.ingest_subset(other)?;
            let header = self.flush(inner)?;
            self.ingest_all(other)?;
            Ok(header)
        }
    }

    /// Takes all bytes from the other into self
    ///
    /// Do not call this directly - always go through `ingest`
    fn ingest_all(&mut self, other: &mut Self) -> Result<()> {
        let n_bytes = other.pos;

        // Drain bounded bytes from other (clearing them in the process)
        self.ubuf.write_all(other.ubuf.drain(..).as_slice())?;

        // Take starts from other (shifting them in the process)
        other
            .starts
            .drain(..)
            .for_each(|start| self.starts.push(start + self.pos));

        // Left shift all remaining starts in other
        other.starts.iter_mut().for_each(|x| {
            *x -= n_bytes;
        });

        // Shift position cursors
        self.pos += n_bytes;

        // Clear the other for good measure
        other.clear();

        Ok(())
    }

    /// Takes as many bytes as possible from the other into self
    ///
    /// Do not call this directly - always go through `ingest`
    fn ingest_subset(&mut self, other: &mut Self) -> Result<()> {
        let remaining = self.block_size - self.pos;
        let (start_index, end_byte) = other
            .starts
            .iter()
            .enumerate()
            .take_while(|(_idx, x)| **x <= remaining)
            .last()
            .map(|(idx, x)| (idx, *x))
            .unwrap();

        // Drain bounded bytes from other (clearing them in the process)
        self.ubuf
            .write_all(other.ubuf.drain(0..end_byte).as_slice())?;

        // Take starts from other (shifting them in the process)
        other
            .starts
            .drain(0..start_index)
            .for_each(|start| self.starts.push(start + self.pos));

        // Left shift all remaining starts in other
        other.starts.iter_mut().for_each(|x| {
            *x -= end_byte;
        });

        // Shift position cursors
        self.pos += end_byte;
        other.pos -= end_byte;

        Ok(())
    }
}

/// Encapsulates the logic for encoding sequences into a binary format.
#[derive(Clone)]
pub struct Encoder {
    /// Bitsize of the nucleotides
    bitsize: BitSize,

    /// Reusable buffers for all nucleotides (written as 2-bit after conversion)
    sbuffer: Vec<u64>,
    xbuffer: Vec<u64>,

    /// Reusable buffers for invalid nucleotide sequences
    s_ibuf: Vec<u8>,
    x_ibuf: Vec<u8>,

    /// Invalid Nucleotide Policy
    policy: Policy,

    /// Random Number Generator
    rng: SmallRng,
}

impl Encoder {
    /// Initialize a new encoder with the given policy.
    pub fn with_policy(bitsize: BitSize, policy: Policy) -> Self {
        Self {
            bitsize,
            policy,
            sbuffer: Vec::default(),
            xbuffer: Vec::default(),
            s_ibuf: Vec::default(),
            x_ibuf: Vec::default(),
            rng: SmallRng::seed_from_u64(RNG_SEED),
        }
    }

    /// Encodes a single sequence as 2-bit.
    ///
    /// Will return `None` if the sequence is invalid and the policy does not allow correction.
    pub fn encode_single(&mut self, primary: &[u8]) -> Result<Option<&[u64]>> {
        // Fill the buffer with the bit representation of the nucleotides
        self.clear();
        if self.bitsize.encode(primary, &mut self.sbuffer).is_err() {
            self.clear();
            if self
                .policy
                .handle(primary, &mut self.s_ibuf, &mut self.rng)?
            {
                self.bitsize.encode(&self.s_ibuf, &mut self.sbuffer)?;
            } else {
                return Ok(None);
            }
        }
        Ok(Some(&self.sbuffer))
    }

    /// Encodes a pair of sequences as 2-bit.
    ///
    /// Will return `None` if either sequence is invalid and the policy does not allow correction.
    pub fn encode_paired(
        &mut self,
        primary: &[u8],
        extended: &[u8],
    ) -> Result<Option<(&[u64], &[u64])>> {
        self.clear();
        if self.bitsize.encode(primary, &mut self.sbuffer).is_err()
            || self.bitsize.encode(extended, &mut self.xbuffer).is_err()
        {
            self.clear();
            if self
                .policy
                .handle(primary, &mut self.s_ibuf, &mut self.rng)?
                && self
                    .policy
                    .handle(extended, &mut self.x_ibuf, &mut self.rng)?
            {
                self.bitsize.encode(&self.s_ibuf, &mut self.sbuffer)?;
                self.bitsize.encode(&self.x_ibuf, &mut self.xbuffer)?;
            } else {
                return Ok(None);
            }
        }
        Ok(Some((&self.sbuffer, &self.xbuffer)))
    }

    /// Clear all buffers and reset the encoder.
    pub fn clear(&mut self) {
        self.sbuffer.clear();
        self.xbuffer.clear();
        self.s_ibuf.clear();
        self.x_ibuf.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vbq::{VBinseqHeaderBuilder, header::SIZE_HEADER};

    #[test]
    fn test_headless_writer() -> super::Result<()> {
        let writer = VBinseqWriterBuilder::default()
            .headless(true)
            .build(Vec::new())?;
        assert_eq!(writer.inner.len(), 0);

        let writer = VBinseqWriterBuilder::default()
            .headless(false)
            .build(Vec::new())?;
        assert_eq!(writer.inner.len(), SIZE_HEADER);

        Ok(())
    }

    #[test]
    fn test_ingest_empty_writer() -> super::Result<()> {
        // Test ingesting from an empty writer
        let header = VBinseqHeaderBuilder::new().build();

        // Create a source writer that's empty
        let mut source = VBinseqWriterBuilder::default()
            .header(header)
            .headless(true)
            .build(Vec::new())?;

        // Create a destination writer
        let mut dest = VBinseqWriterBuilder::default()
            .header(header)
            .headless(true)
            .build(Vec::new())?;

        // Ingest from source to dest
        dest.ingest(&mut source)?;

        // Both writers should be empty
        let source_vec = source.by_ref();
        let dest_vec = dest.by_ref();

        assert_eq!(source_vec.len(), 0);
        assert_eq!(dest_vec.len(), 0);

        Ok(())
    }

    #[test]
    fn test_ingest_single_record() -> super::Result<()> {
        // Test ingesting a single record
        let header = VBinseqHeaderBuilder::new().build();

        // Create a source writer with a single record
        let mut source = VBinseqWriterBuilder::default()
            .header(header)
            .headless(true)
            .build(Vec::new())?;

        // Write a single sequence
        let seq = b"ACGTACGTACGT";
        source.write_record(
            Some(1), // flag
            None,    // header
            seq,     // sequence
            None,    // quality
        )?;

        // We have not crossed a boundary
        assert!(source.by_ref().is_empty());

        // Create a destination writer
        let mut dest = VBinseqWriterBuilder::default()
            .header(header)
            .headless(true)
            .build(Vec::new())?;

        // Ingest from source to dest
        dest.ingest(&mut source)?;

        // Source should be empty, dest should have content
        let source_vec = source.by_ref();
        assert_eq!(source_vec.len(), 0);

        // Source ubuffer should be empty as well
        let source_ubuf = &source.cblock.ubuf;
        assert!(source_ubuf.is_empty());

        // The destination vec will be empty because we haven't hit a buffer limit
        let dest_vec = dest.by_ref();
        assert!(dest_vec.is_empty());

        // The destination ubuffer should have some data however
        let dest_ubuf = &dest.cblock.ubuf;
        assert!(!dest_ubuf.is_empty());

        Ok(())
    }

    #[test]
    fn test_ingest_multi_record() -> super::Result<()> {
        // Test ingesting a single record
        let header = VBinseqHeaderBuilder::new().build();

        // Create a source writer with a single record
        let mut source = VBinseqWriterBuilder::default()
            .header(header)
            .headless(true)
            .build(Vec::new())?;

        // Write multiple sequences
        for _ in 0..30 {
            let seq = b"ACGTACGTACGT";
            source.write_record(Some(1), None, seq, None)?;
        }
        // We have not crossed a boundary
        assert!(source.by_ref().is_empty());

        // Create a destination writer
        let mut dest = VBinseqWriterBuilder::default()
            .header(header)
            .headless(true)
            .build(Vec::new())?;

        // Ingest from source to dest
        dest.ingest(&mut source)?;

        // Source should be empty, dest should have content
        let source_vec = source.by_ref();
        assert_eq!(source_vec.len(), 0);

        // Source ubuffer should be empty as well
        let source_ubuf = &source.cblock.ubuf;
        assert!(source_ubuf.is_empty());

        // The destination vec will be empty because we haven't hit a buffer limit
        let dest_vec = dest.by_ref();
        assert!(dest_vec.is_empty());

        // The destination ubuffer should have some data however
        let dest_ubuf = &dest.cblock.ubuf;
        assert!(!dest_ubuf.is_empty());

        Ok(())
    }

    #[test]
    fn test_ingest_block_boundary() -> super::Result<()> {
        // Test ingesting a single record
        let header = VBinseqHeaderBuilder::new().build();

        // Create a source writer with a single record
        let mut source = VBinseqWriterBuilder::default()
            .header(header)
            .headless(true)
            .build(Vec::new())?;

        // Write multiple sequences (will cross boundary)
        for _ in 0..30000 {
            let seq = b"ACGTACGTACGT";
            source.write_record(Some(1), None, seq, None)?;
        }

        // We have crossed a boundary
        assert!(!source.by_ref().is_empty());

        // Create a destination writer
        let mut dest = VBinseqWriterBuilder::default()
            .header(header)
            .headless(true)
            .build(Vec::new())?;

        // Ingest from source to dest
        dest.ingest(&mut source)?;

        // Source should be empty, dest should have content
        let source_vec = source.by_ref();
        assert_eq!(source_vec.len(), 0);

        // Source ubuffer should be empty as well
        let source_ubuf = &source.cblock.ubuf;
        assert!(source_ubuf.is_empty());

        // The destination vec will not be empty because we hit a buffer limit
        let dest_vec = dest.by_ref();
        assert!(!dest_vec.is_empty());

        // The destination ubuffer should have some data however
        let dest_ubuf = &dest.cblock.ubuf;
        assert!(!dest_ubuf.is_empty());

        Ok(())
    }

    #[test]
    fn test_ingest_with_quality_scores() -> super::Result<()> {
        // Test ingesting records with quality scores
        let source_header = VBinseqHeaderBuilder::new().qual(true).build();
        let dest_header = VBinseqHeaderBuilder::new().qual(true).build();

        // Create a source writer with quality scores
        let mut source = VBinseqWriterBuilder::default()
            .header(source_header)
            .headless(true)
            .build(Vec::new())?;

        // Write sequences with quality scores
        for i in 0..5 {
            let seq = b"ACGTACGTACGT";
            // Simple quality scores (all the same for this test)
            let qual = vec![40; seq.len()];
            source.write_record(Some(i), None, seq, Some(&qual))?;
        }

        // Create a destination writer
        let mut dest = VBinseqWriterBuilder::default()
            .header(dest_header)
            .headless(true)
            .build(Vec::new())?;

        // Ingest from source to dest
        dest.ingest(&mut source)?;

        // Verify source is cleared
        let source_vec = source.by_ref();
        assert_eq!(source_vec.len(), 0);

        // Verify destination has content in ubuf
        let dest_ubuf = &dest.cblock.ubuf;
        assert!(!dest_ubuf.is_empty());

        Ok(())
    }

    #[test]
    fn test_ingest_with_compression() -> super::Result<()> {
        // Test ingesting a single record
        let header = VBinseqHeaderBuilder::new().compressed(true).build();

        // Create a source writer with a single record
        let mut source = VBinseqWriterBuilder::default()
            .header(header)
            .headless(true)
            .build(Vec::new())?;

        // Write multiple sequences (will cross boundary)
        for _ in 0..30000 {
            let seq = b"ACGTACGTACGT";
            source.write_record(Some(1), None, seq, None)?;
        }

        // Create a destination writer
        let mut dest = VBinseqWriterBuilder::default()
            .header(header)
            .headless(true)
            .build(Vec::new())?;

        // Ingest from source to dest
        dest.ingest(&mut source)?;

        // Source should be empty, dest should have content
        let source_vec = source.by_ref();
        assert_eq!(source_vec.len(), 0);

        // Source ubuffer should be empty as well
        let source_ubuf = &source.cblock.ubuf;
        assert!(source_ubuf.is_empty());

        // The destination vec will not be empty because we hit a buffer limit
        let dest_vec = dest.by_ref();
        assert!(!dest_vec.is_empty());

        // The destination ubuffer should have some data however
        let dest_ubuf = &dest.cblock.ubuf;
        assert!(!dest_ubuf.is_empty());

        Ok(())
    }

    #[test]
    fn test_ingest_incompatible_headers() -> super::Result<()> {
        let source_header = VBinseqHeaderBuilder::new().build();
        let dest_header = VBinseqHeaderBuilder::new().qual(true).build();

        // Create a source writer with quality scores
        let mut source = VBinseqWriterBuilder::default()
            .header(source_header)
            .headless(true)
            .build(Vec::new())?;

        // Create a destination writer
        let mut dest = VBinseqWriterBuilder::default()
            .header(dest_header)
            .headless(true)
            .build(Vec::new())?;

        // Ingest from source to dest (will error)
        assert!(dest.ingest(&mut source).is_err());

        Ok(())
    }

    #[test]
    #[allow(clippy::identity_op)]
    fn test_record_byte_size() {
        let size = record_byte_size(2, 0, true);
        assert_eq!(size, 8 * (2 + 0 + 3)); // 40 bytes

        let size = record_byte_size(4, 8, true);
        assert_eq!(size, 8 * (4 + 8 + 3)); // 128 bytes
    }
}
