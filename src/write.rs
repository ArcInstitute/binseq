//! Unified writer interface for BINSEQ formats
//!
//! This module provides a unified [`BinseqWriter`] enum that abstracts over the three
//! BINSEQ format writers (BQ, VBQ, CBQ), allowing format-agnostic writing of sequence data.
//!
//! # Example
//!
//! ```rust
//! use binseq::{write::{BinseqWriter, BinseqWriterBuilder, Format}, SequencingRecordBuilder};
//! use std::io::Cursor;
//!
//! // Create a VBQ writer with quality scores and headers
//! let mut writer = BinseqWriterBuilder::new(Format::Vbq)
//!     .paired(false)
//!     .quality(true)
//!     .headers(true)
//!     .build(Cursor::new(Vec::new()))
//!     .unwrap();
//!
//! // Write a record
//! let record = SequencingRecordBuilder::default()
//!     .s_seq(b"ACGTACGT")
//!     .s_qual(b"IIIIIIII")
//!     .s_header(b"seq1")
//!     .build()
//!     .unwrap();
//!
//! writer.push(record).unwrap();
//! writer.finish().unwrap();
//! ```
//!
//! # Parallel Writing
//!
//! For parallel writing scenarios, use `headless(true)` for thread-local writers
//! and `ingest()` to merge them into a global writer:
//!
//! ```rust,no_run
//! use binseq::{write::{BinseqWriter, BinseqWriterBuilder, Format}, SequencingRecordBuilder};
//! use std::fs::File;
//!
//! // Global writer (writes header)
//! let mut global = BinseqWriterBuilder::new(Format::Vbq)
//!     .paired(false)
//!     .build(File::create("output.vbq").unwrap())
//!     .unwrap();
//!
//! // Thread-local writer (headless, Vec<u8> buffer)
//! let mut local = global.new_headless_buffer().unwrap();
//!
//! // Write to local buffer
//! let record = SequencingRecordBuilder::default()
//!     .s_seq(b"ACGTACGT")
//!     .build()
//!     .unwrap();
//! local.push(record).unwrap();
//!
//! // Merge into global writer
//! global.ingest(&mut local).unwrap();
//! global.finish().unwrap();
//! ```

use std::{io::Write, str::FromStr};

use crate::{BitSize, Policy, Result, SequencingRecord, bq, cbq, error::WriteError, vbq};

/// Output format for BINSEQ files
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Format {
    /// BQ format - fixed length records, no quality scores
    Bq,
    /// VBQ format - variable length records, optional quality scores
    Vbq,
    /// CBQ format - columnar variable length records, optional quality scores
    #[default]
    Cbq,
}
impl FromStr for Format {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "bq" | "BQ" | "b" => Ok(Self::Bq),
            "vbq" | "VBQ" | "v" => Ok(Self::Vbq),
            "cbq" | "CBQ" | "c" => Ok(Self::Cbq),
            _ => Err(format!("Unknown format: {}", s)),
        }
    }
}

impl Format {
    /// Returns the file extension for this format (including the dot)
    #[must_use]
    pub fn extension(&self) -> &'static str {
        match self {
            Self::Bq => ".bq",
            Self::Vbq => ".vbq",
            Self::Cbq => ".cbq",
        }
    }
}

/// Builder for creating [`BinseqWriter`] instances
///
/// This builder provides a unified interface for configuring writers across all
/// BINSEQ formats. Settings that don't apply to a particular format are silently
/// ignored.
///
/// # Format-specific behavior
///
/// | Setting | BQ | VBQ | CBQ |
/// |---------|:--:|:---:|:---:|
/// | `quality(true)` | ignored | applied | applied |
/// | `headers(true)` | ignored | applied | applied |
/// | `compression(true)` | ignored | applied | applied |
/// | `compression_level(n)` | ignored | ignored | applied |
/// | `block_size(n)` | ignored | applied | applied |
/// | `bitsize(b)` | applied | applied | ignored |
/// | `slen(n)` | **required** | ignored | ignored |
/// | `xlen(n)` | required if paired | ignored | ignored |
/// | `policy(p)` | applied | applied | ignored |
/// | `headless(true)` | applied | applied | applied |
#[derive(Debug, Clone)]
pub struct BinseqWriterBuilder {
    format: Format,
    paired: bool,
    quality: bool,
    headers: bool,
    flags: bool,
    compression: bool,
    compression_level: Option<i32>,
    block_size: Option<usize>,
    policy: Option<Policy>,
    headless: bool,
    bitsize: Option<BitSize>,
    slen: Option<u32>,
    xlen: Option<u32>,
}

impl BinseqWriterBuilder {
    /// Create a new builder for the specified format
    #[must_use]
    pub fn new(format: Format) -> Self {
        Self {
            format,
            paired: false,
            quality: false,
            headers: false,
            flags: false,
            compression: true,
            compression_level: None,
            block_size: None,
            policy: None,
            headless: false,
            bitsize: None,
            slen: None,
            xlen: None,
        }
    }

    /// Set whether records are paired-end
    #[must_use]
    pub fn paired(mut self, paired: bool) -> Self {
        self.paired = paired;
        self
    }

    /// Set whether to store quality scores (ignored for BQ)
    #[must_use]
    pub fn quality(mut self, quality: bool) -> Self {
        self.quality = quality;
        self
    }

    /// Set whether to store sequence headers (ignored for BQ)
    #[must_use]
    pub fn headers(mut self, headers: bool) -> Self {
        self.headers = headers;
        self
    }

    /// Set whether to store flags
    #[must_use]
    pub fn flags(mut self, flags: bool) -> Self {
        self.flags = flags;
        self
    }

    /// Set whether to compress data (ignored for BQ)
    #[must_use]
    pub fn compression(mut self, compression: bool) -> Self {
        self.compression = compression;
        self
    }

    /// Set the compression level (only applies to CBQ)
    #[must_use]
    pub fn compression_level(mut self, level: i32) -> Self {
        self.compression_level = Some(level);
        self
    }

    /// Set the block size in bytes (ignored for BQ)
    #[must_use]
    pub fn block_size(mut self, size: usize) -> Self {
        self.block_size = Some(size);
        self
    }

    /// Set the policy for handling invalid nucleotides (ignored for CBQ)
    #[must_use]
    pub fn policy(mut self, policy: Policy) -> Self {
        self.policy = Some(policy);
        self
    }

    /// Set whether to operate in headless mode (for parallel writing)
    #[must_use]
    pub fn headless(mut self, headless: bool) -> Self {
        self.headless = headless;
        self
    }

    /// Set the bit size for nucleotide encoding (ignored for CBQ)
    #[must_use]
    pub fn bitsize(mut self, bitsize: BitSize) -> Self {
        self.bitsize = Some(bitsize);
        self
    }

    /// Set the primary sequence length (required for BQ, ignored for VBQ/CBQ)
    #[must_use]
    pub fn slen(mut self, len: u32) -> Self {
        self.slen = Some(len);
        self
    }

    /// Set the extended sequence length (required for paired BQ, ignored for VBQ/CBQ)
    #[must_use]
    pub fn xlen(mut self, len: u32) -> Self {
        self.xlen = Some(len);
        self
    }

    /// Build the writer
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - Format is BQ and `slen` is not set
    /// - Format is BQ, `paired` is true, but `xlen` is not set
    pub fn build<W: Write>(self, writer: W) -> Result<BinseqWriter<W>> {
        match self.format {
            Format::Bq => self.build_bq(writer),
            Format::Vbq => self.build_vbq(writer),
            Format::Cbq => self.build_cbq(writer),
        }
    }

    fn build_bq<W: Write>(self, writer: W) -> Result<BinseqWriter<W>> {
        let slen = self.slen.ok_or(WriteError::MissingSequenceLength {
            exp_primary: true,
            exp_extended: self.paired,
            obs_primary: self.slen.is_some(),
            obs_extended: self.xlen.is_some(),
        })?;
        let xlen = if self.paired || self.xlen.is_some_and(|x| x > 0) {
            self.xlen.ok_or(WriteError::MissingSequenceLength {
                exp_primary: true,
                exp_extended: true,
                obs_primary: self.slen.is_some(),
                obs_extended: self.xlen.is_some(),
            })?
        } else {
            0
        };

        let mut header_builder = bq::BinseqHeaderBuilder::new().slen(slen).xlen(xlen);

        if let Some(bitsize) = self.bitsize {
            header_builder = header_builder.bitsize(bitsize);
        }

        header_builder = header_builder.flags(self.flags);

        let header = header_builder.build()?;

        let inner = bq::BinseqWriterBuilder::default()
            .header(header)
            .policy(self.policy.unwrap_or_default())
            .headless(self.headless)
            .build(writer)?;

        Ok(BinseqWriter::Bq(inner))
    }

    fn build_vbq<W: Write>(self, writer: W) -> Result<BinseqWriter<W>> {
        let mut header_builder = vbq::VBinseqHeaderBuilder::new()
            .paired(self.paired)
            .qual(self.quality)
            .headers(self.headers)
            .flags(self.flags)
            .compressed(self.compression);

        if let Some(block_size) = self.block_size {
            header_builder = header_builder.block(block_size as u64);
        }

        if let Some(bitsize) = self.bitsize {
            header_builder = header_builder.bitsize(bitsize);
        }

        let header = header_builder.build();

        let inner = vbq::VBinseqWriterBuilder::default()
            .header(header)
            .policy(self.policy.unwrap_or_default())
            .headless(self.headless)
            .build(writer)?;

        Ok(BinseqWriter::Vbq(inner))
    }

    fn build_cbq<W: Write>(self, writer: W) -> Result<BinseqWriter<W>> {
        let header = cbq::FileHeaderBuilder::default()
            .is_paired(self.paired)
            .with_qualities(self.quality)
            .with_headers(self.headers)
            .with_flags(self.flags)
            .with_optional_block_size(self.block_size)
            .with_optional_compression_level(self.compression_level.map(|level| level as usize))
            .build();

        let inner = if self.headless {
            cbq::ColumnarBlockWriter::new_headless(writer, header)?
        } else {
            cbq::ColumnarBlockWriter::new(writer, header)?
        };

        Ok(BinseqWriter::Cbq(inner))
    }
}

/// Unified writer for BINSEQ formats
///
/// This enum wraps the three format-specific writers (BQ, VBQ, CBQ) and provides
/// a unified interface for writing sequence data.
pub enum BinseqWriter<W: Write> {
    /// BQ format writer
    Bq(bq::BinseqWriter<W>),
    /// VBQ format writer
    Vbq(vbq::VBinseqWriter<W>),
    /// CBQ format writer
    Cbq(cbq::ColumnarBlockWriter<W>),
}

impl<W: Write> BinseqWriter<W> {
    /// Push a record to the writer
    ///
    /// Returns `Ok(true)` if the record was written successfully, or `Ok(false)`
    /// if the record was skipped due to invalid nucleotides (based on the configured
    /// policy). CBQ always returns `Ok(true)` as it handles N's explicitly.
    ///
    /// # Errors
    ///
    /// Returns an error if there's an I/O error or if the record doesn't match
    /// the writer's configuration (e.g., paired record to unpaired writer).
    pub fn push(&mut self, record: SequencingRecord) -> Result<bool> {
        match self {
            Self::Bq(w) => w.push(record),
            Self::Vbq(w) => w.push(record),
            Self::Cbq(w) => w.push(record),
        }
    }

    /// Finish writing and flush any remaining data
    ///
    /// For VBQ and CBQ formats, this writes the embedded index. For BQ, this
    /// is equivalent to `flush()`.
    ///
    /// # Errors
    ///
    /// Returns an error if there's an I/O error writing the final data.
    pub fn finish(&mut self) -> Result<()> {
        match self {
            Self::Bq(w) => w.flush().map_err(Into::into),
            Self::Vbq(w) => w.finish(),
            Self::Cbq(w) => w.finish(),
        }
    }

    /// Returns the format of this writer
    #[must_use]
    pub fn format(&self) -> Format {
        match self {
            Self::Bq(_) => Format::Bq,
            Self::Vbq(_) => Format::Vbq,
            Self::Cbq(_) => Format::Cbq,
        }
    }

    /// Returns whether this writer is configured for paired-end records
    #[must_use]
    pub fn is_paired(&self) -> bool {
        match self {
            Self::Bq(w) => w.is_paired(),
            Self::Vbq(w) => w.is_paired(),
            Self::Cbq(w) => w.header().is_paired(),
        }
    }

    /// Returns whether this writer stores quality scores
    ///
    /// Always returns `false` for BQ format.
    #[must_use]
    pub fn has_quality(&self) -> bool {
        match self {
            Self::Bq(_) => false,
            Self::Vbq(w) => w.has_quality(),
            Self::Cbq(w) => w.header().has_qualities(),
        }
    }

    /// Returns whether this writer stores sequence headers
    ///
    /// Always returns `false` for BQ format.
    #[must_use]
    pub fn has_headers(&self) -> bool {
        match self {
            Self::Bq(_) => false,
            Self::Vbq(w) => w.has_headers(),
            Self::Cbq(w) => w.header().has_headers(),
        }
    }
}

impl<W: Write + Clone> Clone for BinseqWriter<W> {
    fn clone(&self) -> Self {
        match self {
            Self::Bq(w) => Self::Bq(w.clone()),
            Self::Vbq(w) => Self::Vbq(w.clone()),
            Self::Cbq(w) => Self::Cbq(w.clone()),
        }
    }
}

impl<W: Write> BinseqWriter<W> {
    /// Ingest records from a headless `Vec<u8>` writer into this writer
    ///
    /// This is used in parallel writing scenarios where thread-local writers
    /// buffer to `Vec<u8>` and then get merged into a global writer.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The source and destination writers have different formats
    /// - The source and destination writers have incompatible headers
    /// - There's an I/O error during ingestion
    pub fn ingest(&mut self, other: &mut BinseqWriter<Vec<u8>>) -> Result<()> {
        match (self, other) {
            (Self::Bq(dst), BinseqWriter::Bq(src)) => dst.ingest(src),
            (Self::Vbq(dst), BinseqWriter::Vbq(src)) => dst.ingest(src),
            (Self::Cbq(dst), BinseqWriter::Cbq(src)) => dst.ingest(src),
            _ => Err(WriteError::FormatMismatch.into()),
        }
    }
}

impl<W: Write> BinseqWriter<W> {
    /// Create a new headless writer with the same configuration, using a `Vec<u8>` buffer
    ///
    /// This is useful for parallel writing scenarios where each thread has its own
    /// buffer that gets merged into a global writer via `ingest()`.
    ///
    /// # Errors
    ///
    /// Returns an error if the writer cannot be created.
    pub fn new_headless_buffer(&self) -> Result<BinseqWriter<Vec<u8>>> {
        match self {
            Self::Bq(w) => {
                let inner = bq::BinseqWriterBuilder::default()
                    .header(w.header())
                    .policy(w.policy())
                    .headless(true)
                    .build(Vec::new())?;
                Ok(BinseqWriter::Bq(inner))
            }
            Self::Vbq(w) => {
                let inner = vbq::VBinseqWriterBuilder::default()
                    .header(w.header())
                    .policy(w.policy())
                    .headless(true)
                    .build(Vec::new())?;
                Ok(BinseqWriter::Vbq(inner))
            }
            Self::Cbq(w) => {
                let inner = cbq::ColumnarBlockWriter::new_headless(Vec::new(), w.header())?;
                Ok(BinseqWriter::Cbq(inner))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SequencingRecordBuilder;
    use std::io::Cursor;

    #[test]
    fn test_format_extension() {
        assert_eq!(Format::Bq.extension(), ".bq");
        assert_eq!(Format::Vbq.extension(), ".vbq");
        assert_eq!(Format::Cbq.extension(), ".cbq");
    }

    #[test]
    fn test_build_bq_writer() -> Result<()> {
        let writer = BinseqWriterBuilder::new(Format::Bq)
            .slen(100)
            .paired(false)
            .build(Cursor::new(Vec::new()))?;

        assert_eq!(writer.format(), Format::Bq);
        assert!(!writer.is_paired());
        assert!(!writer.has_quality());
        assert!(!writer.has_headers());
        Ok(())
    }

    #[test]
    fn test_build_bq_writer_paired() -> Result<()> {
        let writer = BinseqWriterBuilder::new(Format::Bq)
            .slen(100)
            .xlen(150)
            .paired(true)
            .build(Cursor::new(Vec::new()))?;

        assert_eq!(writer.format(), Format::Bq);
        assert!(writer.is_paired());
        Ok(())
    }

    #[test]
    fn test_build_bq_missing_slen() {
        let result = BinseqWriterBuilder::new(Format::Bq)
            .paired(false)
            .build(Cursor::new(Vec::new()));

        assert!(result.is_err());
    }

    #[test]
    fn test_build_bq_paired_missing_xlen() {
        let result = BinseqWriterBuilder::new(Format::Bq)
            .slen(100)
            .paired(true)
            .build(Cursor::new(Vec::new()));

        assert!(result.is_err());
    }

    #[test]
    fn test_build_vbq_writer() -> Result<()> {
        let writer = BinseqWriterBuilder::new(Format::Vbq)
            .paired(true)
            .quality(true)
            .headers(true)
            .build(Cursor::new(Vec::new()))?;

        assert_eq!(writer.format(), Format::Vbq);
        assert!(writer.is_paired());
        assert!(writer.has_quality());
        assert!(writer.has_headers());
        Ok(())
    }

    #[test]
    fn test_build_cbq_writer() -> Result<()> {
        let writer = BinseqWriterBuilder::new(Format::Cbq)
            .paired(false)
            .quality(true)
            .headers(true)
            .compression_level(3)
            .build(Cursor::new(Vec::new()))?;

        assert_eq!(writer.format(), Format::Cbq);
        assert!(!writer.is_paired());
        assert!(writer.has_quality());
        assert!(writer.has_headers());
        Ok(())
    }

    #[test]
    fn test_push_and_finish_vbq() -> Result<()> {
        let mut writer = BinseqWriterBuilder::new(Format::Vbq)
            .paired(false)
            .quality(false)
            .headers(false)
            .build(Cursor::new(Vec::new()))?;

        let record = SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGTACGT")
            .build()?;

        let written = writer.push(record)?;
        assert!(written);

        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_push_and_finish_cbq() -> Result<()> {
        let mut writer = BinseqWriterBuilder::new(Format::Cbq)
            .paired(false)
            .quality(false)
            .headers(false)
            .build(Cursor::new(Vec::new()))?;

        let record = SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGTACGT")
            .build()?;

        let written = writer.push(record)?;
        assert!(written);

        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_push_and_finish_bq() -> Result<()> {
        let mut writer = BinseqWriterBuilder::new(Format::Bq)
            .slen(12)
            .paired(false)
            .build(Cursor::new(Vec::new()))?;

        let record = SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGTACGT")
            .build()?;

        let written = writer.push(record)?;
        assert!(written);

        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_new_headless_buffer_vbq() -> Result<()> {
        let global = BinseqWriterBuilder::new(Format::Vbq)
            .paired(true)
            .quality(true)
            .headers(true)
            .build(Cursor::new(Vec::new()))?;

        let local = global.new_headless_buffer()?;

        assert_eq!(local.format(), Format::Vbq);
        assert!(local.is_paired());
        assert!(local.has_quality());
        assert!(local.has_headers());
        Ok(())
    }

    #[test]
    fn test_new_headless_buffer_cbq() -> Result<()> {
        let global = BinseqWriterBuilder::new(Format::Cbq)
            .paired(false)
            .quality(true)
            .build(Cursor::new(Vec::new()))?;

        let local = global.new_headless_buffer()?;

        assert_eq!(local.format(), Format::Cbq);
        assert!(!local.is_paired());
        assert!(local.has_quality());
        Ok(())
    }

    #[test]
    fn test_new_headless_buffer_bq() -> Result<()> {
        let global = BinseqWriterBuilder::new(Format::Bq)
            .slen(100)
            .xlen(150)
            .paired(true)
            .build(Cursor::new(Vec::new()))?;

        let local = global.new_headless_buffer()?;

        assert_eq!(local.format(), Format::Bq);
        assert!(local.is_paired());
        Ok(())
    }

    #[test]
    fn test_ingest_vbq() -> Result<()> {
        let mut global = BinseqWriterBuilder::new(Format::Vbq)
            .paired(false)
            .quality(false)
            .headers(false)
            .build(Cursor::new(Vec::new()))?;

        let mut local = global.new_headless_buffer()?;

        // Write to local
        let record = SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGTACGT")
            .build()?;
        local.push(record)?;

        // Ingest into global
        global.ingest(&mut local)?;
        global.finish()?;

        Ok(())
    }

    #[test]
    fn test_ingest_cbq() -> Result<()> {
        let mut global = BinseqWriterBuilder::new(Format::Cbq)
            .paired(false)
            .quality(false)
            .headers(false)
            .build(Cursor::new(Vec::new()))?;

        let mut local = global.new_headless_buffer()?;

        // Write to local
        let record = SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGTACGT")
            .build()?;
        local.push(record)?;

        // Ingest into global
        global.ingest(&mut local)?;
        global.finish()?;

        Ok(())
    }

    #[test]
    fn test_ingest_bq() -> Result<()> {
        let mut global = BinseqWriterBuilder::new(Format::Bq)
            .slen(12)
            .paired(false)
            .build(Cursor::new(Vec::new()))?;

        let mut local = global.new_headless_buffer()?;

        // Write to local
        let record = SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGTACGT")
            .build()?;
        local.push(record)?;

        // Ingest into global
        global.ingest(&mut local)?;
        global.finish()?;

        Ok(())
    }

    #[test]
    fn test_ingest_format_mismatch() -> Result<()> {
        let mut global = BinseqWriterBuilder::new(Format::Vbq)
            .paired(false)
            .build(Cursor::new(Vec::new()))?;

        let mut local = BinseqWriterBuilder::new(Format::Cbq)
            .paired(false)
            .headless(true)
            .build(Vec::new())?;

        let result = global.ingest(&mut local);
        assert!(result.is_err());

        Ok(())
    }
}
