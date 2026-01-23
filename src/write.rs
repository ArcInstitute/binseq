//! Unified writer interface for BINSEQ formats
//!
//! This module provides a unified `BinseqWriter` enum that abstracts over the three
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
            _ => Err(format!("Unknown format: {s}")),
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
    pub(crate) format: Format,
    pub(crate) paired: bool,
    quality: bool,
    headers: bool,
    flags: bool,
    compression: bool,
    compression_level: Option<i32>,
    block_size: Option<usize>,
    policy: Option<Policy>,
    headless: bool,
    bitsize: Option<BitSize>,
    pub(crate) slen: Option<u32>,
    pub(crate) xlen: Option<u32>,
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

    /// Sets the corresponding values for this builder given an existing BQ header
    #[must_use]
    pub fn from_bq_header(header: bq::FileHeader) -> Self {
        Self {
            format: Format::Bq,
            slen: Some(header.slen),
            xlen: (header.xlen > 0).then_some(header.xlen),
            bitsize: Some(header.bits),
            paired: header.is_paired(),
            flags: header.flags,
            compression: false,
            headers: false,
            quality: false,
            compression_level: None,
            block_size: None,
            headless: false,
            policy: None,
        }
    }

    /// Sets the corresponding values for this builder given an existing VBQ header
    #[must_use]
    pub fn from_vbq_header(header: vbq::FileHeader) -> Self {
        Self {
            format: Format::Vbq,
            slen: None,
            xlen: None,
            flags: header.flags,
            quality: header.qual,
            paired: header.paired,
            bitsize: Some(header.bits),
            headers: header.headers,
            compression: header.compressed,
            block_size: Some(header.block as usize),
            policy: None,
            compression_level: None,
            headless: false,
        }
    }

    /// Sets the corresponding values for this builder given an existing CBQ header
    #[must_use]
    pub fn from_cbq_header(header: cbq::FileHeader) -> Self {
        Self {
            format: Format::Cbq,
            flags: header.has_flags(),
            quality: header.has_qualities(),
            headers: header.has_headers(),
            paired: header.is_paired(),
            block_size: Some(header.block_size as usize),
            compression_level: Some(header.compression_level as i32),
            compression: false,
            slen: None,
            xlen: None,
            bitsize: None,
            policy: None,
            headless: false,
        }
    }

    /// Encode FASTX file(s) to BINSEQ format
    ///
    /// This method returns a [`FastxEncoderBuilder`] that allows you to configure
    /// the input source and threading options before executing the encoding.
    ///
    /// This is an alternative to [`build`](Self::build) that directly processes
    /// FASTX files using parallel processing.
    ///
    /// # Availability
    ///
    /// This method is only available when the `paraseq` feature is enabled.
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use binseq::write::{BinseqWriterBuilder, Format};
    /// use std::fs::File;
    ///
    /// // Encode from stdin to VBQ
    /// let writer = BinseqWriterBuilder::new(Format::Vbq)
    ///     .quality(true)
    ///     .headers(true)
    ///     .encode_fastx(File::create("output.vbq")?)
    ///     .input_stdin()
    ///     .threads(8)
    ///     .run()?;
    ///
    /// // Encode paired-end reads
    /// let writer = BinseqWriterBuilder::new(Format::Vbq)
    ///     .quality(true)
    ///     .encode_fastx(File::create("output.vbq")?)
    ///     .input_paired("R1.fastq", "R2.fastq")
    ///     .run()?;
    /// # Ok::<(), binseq::Error>(())
    /// ```
    #[cfg(feature = "paraseq")]
    #[must_use]
    pub fn encode_fastx<W: Write + Send + 'static>(
        self,
        output: W,
    ) -> crate::utils::FastxEncoderBuilder {
        crate::utils::FastxEncoderBuilder::new(self, Box::new(output))
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

        let mut header_builder = bq::FileHeaderBuilder::new().slen(slen).xlen(xlen);

        if let Some(bitsize) = self.bitsize {
            header_builder = header_builder.bitsize(bitsize);
        }

        header_builder = header_builder.flags(self.flags);

        let header = header_builder.build()?;

        let inner = bq::WriterBuilder::default()
            .header(header)
            .policy(self.policy.unwrap_or_default())
            .headless(self.headless)
            .build(writer)?;

        Ok(BinseqWriter::Bq(inner))
    }

    fn build_vbq<W: Write>(self, writer: W) -> Result<BinseqWriter<W>> {
        let mut header_builder = vbq::FileHeaderBuilder::new()
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
    Bq(bq::Writer<W>),
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
            Self::Bq(w) => w.flush(),
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
                let inner = bq::WriterBuilder::default()
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

    // ==================== Record Specification Tests ====================
    //
    // These tests verify that writers correctly handle records with different
    // levels of specification relative to the writer's configuration:
    // - Under-specified: record is missing data the writer needs (should error)
    // - Over-specified: record has extra data the writer ignores (should succeed)
    // - Correctly-specified: record matches writer config exactly (should succeed)

    /// Helper to create a minimal single-end record (sequence only)
    fn minimal_single_record() -> SequencingRecord<'static> {
        SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGTACGTACGTACGTACGTACGTACGT")
            .build()
            .unwrap()
    }

    /// Helper to create a minimal paired record (sequences only)
    fn minimal_paired_record() -> SequencingRecord<'static> {
        SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGTACGTACGTACGTACGTACGTACGT")
            .x_seq(b"TGCATGCATGCATGCATGCATGCATGCATGCA")
            .build()
            .unwrap()
    }

    /// Helper to create a fully-specified single-end record
    fn full_single_record() -> SequencingRecord<'static> {
        SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGTACGTACGTACGTACGTACGTACGT")
            .s_qual(b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII")
            .s_header(b"read1")
            .flag(42u64)
            .build()
            .unwrap()
    }

    /// Helper to create a fully-specified paired record
    fn full_paired_record() -> SequencingRecord<'static> {
        SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGTACGTACGTACGTACGTACGTACGT")
            .s_qual(b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII")
            .s_header(b"read1")
            .x_seq(b"TGCATGCATGCATGCATGCATGCATGCATGCA")
            .x_qual(b"JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ")
            .x_header(b"read2")
            .flag(42u64)
            .build()
            .unwrap()
    }

    // ==================== VBQ Tests ====================

    #[test]
    fn test_vbq_single_minimal_writer_minimal_record() -> Result<()> {
        // Writer: single-end, no quality, no headers, no flags
        // Record: single-end, no quality, no headers, no flags
        // Expected: success (correctly specified)
        let mut writer = BinseqWriterBuilder::new(Format::Vbq)
            .paired(false)
            .quality(false)
            .headers(false)
            .flags(false)
            .build(Cursor::new(Vec::new()))?;

        let record = minimal_single_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_vbq_single_minimal_writer_full_record() -> Result<()> {
        // Writer: single-end, no quality, no headers, no flags
        // Record: single-end, with quality, headers, flags
        // Expected: success (over-specified - extra data ignored)
        let mut writer = BinseqWriterBuilder::new(Format::Vbq)
            .paired(false)
            .quality(false)
            .headers(false)
            .flags(false)
            .build(Cursor::new(Vec::new()))?;

        let record = full_single_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_vbq_single_full_writer_minimal_record() -> Result<()> {
        // Writer: single-end, with quality, headers, flags
        // Record: single-end, no quality, no headers, no flags
        // Expected: error (under-specified)
        let mut writer = BinseqWriterBuilder::new(Format::Vbq)
            .paired(false)
            .quality(true)
            .headers(true)
            .flags(true)
            .build(Cursor::new(Vec::new()))?;

        let record = minimal_single_record();
        let result = writer.push(record);
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn test_vbq_single_full_writer_full_record() -> Result<()> {
        // Writer: single-end, with quality, headers, flags
        // Record: single-end, with quality, headers, flags
        // Expected: success (correctly specified)
        let mut writer = BinseqWriterBuilder::new(Format::Vbq)
            .paired(false)
            .quality(true)
            .headers(true)
            .flags(true)
            .build(Cursor::new(Vec::new()))?;

        let record = full_single_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_vbq_paired_writer_single_record() -> Result<()> {
        // Writer: paired
        // Record: single-end
        // Expected: error (under-specified - missing R2)
        let mut writer = BinseqWriterBuilder::new(Format::Vbq)
            .paired(true)
            .quality(false)
            .headers(false)
            .flags(false)
            .build(Cursor::new(Vec::new()))?;

        let record = minimal_single_record();
        let result = writer.push(record);
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn test_vbq_single_writer_paired_record() -> Result<()> {
        // Writer: single-end
        // Record: paired
        // Expected: success (over-specified - R2 ignored)
        let mut writer = BinseqWriterBuilder::new(Format::Vbq)
            .paired(false)
            .quality(false)
            .headers(false)
            .flags(false)
            .build(Cursor::new(Vec::new()))?;

        let record = minimal_paired_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_vbq_paired_minimal_writer_paired_full_record() -> Result<()> {
        // Writer: paired, no quality, no headers, no flags
        // Record: paired, with quality, headers, flags
        // Expected: success (over-specified)
        let mut writer = BinseqWriterBuilder::new(Format::Vbq)
            .paired(true)
            .quality(false)
            .headers(false)
            .flags(false)
            .build(Cursor::new(Vec::new()))?;

        let record = full_paired_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_vbq_paired_full_writer_paired_full_record() -> Result<()> {
        // Writer: paired, with quality, headers, flags
        // Record: paired, with quality, headers, flags
        // Expected: success (correctly specified)
        let mut writer = BinseqWriterBuilder::new(Format::Vbq)
            .paired(true)
            .quality(true)
            .headers(true)
            .flags(true)
            .build(Cursor::new(Vec::new()))?;

        let record = full_paired_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    // ==================== CBQ Tests ====================

    #[test]
    fn test_cbq_single_minimal_writer_minimal_record() -> Result<()> {
        // Writer: single-end, no quality, no headers, no flags
        // Record: single-end, no quality, no headers, no flags
        // Expected: success (correctly specified)
        let mut writer = BinseqWriterBuilder::new(Format::Cbq)
            .paired(false)
            .quality(false)
            .headers(false)
            .flags(false)
            .build(Cursor::new(Vec::new()))?;

        let record = minimal_single_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_cbq_single_minimal_writer_full_record() -> Result<()> {
        // Writer: single-end, no quality, no headers, no flags
        // Record: single-end, with quality, headers, flags
        // Expected: success (over-specified - extra data ignored)
        let mut writer = BinseqWriterBuilder::new(Format::Cbq)
            .paired(false)
            .quality(false)
            .headers(false)
            .flags(false)
            .build(Cursor::new(Vec::new()))?;

        let record = full_single_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_cbq_single_full_writer_minimal_record() -> Result<()> {
        // Writer: single-end, with quality, headers, flags
        // Record: single-end, no quality, no headers, no flags
        // Expected: error (under-specified)
        let mut writer = BinseqWriterBuilder::new(Format::Cbq)
            .paired(false)
            .quality(true)
            .headers(true)
            .flags(true)
            .build(Cursor::new(Vec::new()))?;

        let record = minimal_single_record();
        let result = writer.push(record);
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn test_cbq_single_full_writer_full_record() -> Result<()> {
        // Writer: single-end, with quality, headers, flags
        // Record: single-end, with quality, headers, flags
        // Expected: success (correctly specified)
        let mut writer = BinseqWriterBuilder::new(Format::Cbq)
            .paired(false)
            .quality(true)
            .headers(true)
            .flags(true)
            .build(Cursor::new(Vec::new()))?;

        let record = full_single_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_cbq_paired_writer_single_record() -> Result<()> {
        // Writer: paired
        // Record: single-end
        // Expected: error (under-specified - missing R2)
        let mut writer = BinseqWriterBuilder::new(Format::Cbq)
            .paired(true)
            .quality(false)
            .headers(false)
            .flags(false)
            .build(Cursor::new(Vec::new()))?;

        let record = minimal_single_record();
        let result = writer.push(record);
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn test_cbq_single_writer_paired_record() -> Result<()> {
        // Writer: single-end
        // Record: paired
        // Expected: success (over-specified - R2 ignored)
        let mut writer = BinseqWriterBuilder::new(Format::Cbq)
            .paired(false)
            .quality(false)
            .headers(false)
            .flags(false)
            .build(Cursor::new(Vec::new()))?;

        let record = minimal_paired_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_cbq_paired_minimal_writer_paired_full_record() -> Result<()> {
        // Writer: paired, no quality, no headers, no flags
        // Record: paired, with quality, headers, flags
        // Expected: success (over-specified)
        let mut writer = BinseqWriterBuilder::new(Format::Cbq)
            .paired(true)
            .quality(false)
            .headers(false)
            .flags(false)
            .build(Cursor::new(Vec::new()))?;

        let record = full_paired_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_cbq_paired_full_writer_paired_full_record() -> Result<()> {
        // Writer: paired, with quality, headers, flags
        // Record: paired, with quality, headers, flags
        // Expected: success (correctly specified)
        let mut writer = BinseqWriterBuilder::new(Format::Cbq)
            .paired(true)
            .quality(true)
            .headers(true)
            .flags(true)
            .build(Cursor::new(Vec::new()))?;

        let record = full_paired_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    // ==================== BQ Tests ====================
    // Note: BQ format has fixed-length sequences and doesn't support headers

    #[test]
    fn test_bq_single_minimal_writer_minimal_record() -> Result<()> {
        // Writer: single-end, no quality, no flags
        // Record: single-end, no quality, no flags
        // Expected: success (correctly specified)
        let mut writer = BinseqWriterBuilder::new(Format::Bq)
            .slen(32)
            .paired(false)
            .quality(false)
            .flags(false)
            .build(Cursor::new(Vec::new()))?;

        let record = minimal_single_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_bq_single_minimal_writer_full_record() -> Result<()> {
        // Writer: single-end, no quality, no flags
        // Record: single-end, with quality, headers, flags
        // Expected: success (over-specified - extra data ignored)
        let mut writer = BinseqWriterBuilder::new(Format::Bq)
            .slen(32)
            .paired(false)
            .quality(false)
            .flags(false)
            .build(Cursor::new(Vec::new()))?;

        let record = full_single_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_bq_single_with_quality_writer_minimal_record() -> Result<()> {
        // Writer: single-end, with quality (note: BQ ignores quality setting)
        // Record: single-end, no quality
        // Expected: success (BQ format doesn't support quality scores, setting is ignored)
        let mut writer = BinseqWriterBuilder::new(Format::Bq)
            .slen(32)
            .paired(false)
            .quality(true) // This is ignored for BQ format
            .build(Cursor::new(Vec::new()))?;

        // BQ always reports has_quality as false
        assert!(!writer.has_quality());

        let record = minimal_single_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_bq_single_with_quality_writer_full_record() -> Result<()> {
        // Writer: single-end, with quality
        // Record: single-end, with quality
        // Expected: success (correctly specified)
        let mut writer = BinseqWriterBuilder::new(Format::Bq)
            .slen(32)
            .paired(false)
            .quality(true)
            .build(Cursor::new(Vec::new()))?;

        let record = full_single_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_bq_paired_writer_single_record() -> Result<()> {
        // Writer: paired
        // Record: single-end
        // Expected: error (under-specified - missing R2)
        let mut writer = BinseqWriterBuilder::new(Format::Bq)
            .slen(32)
            .xlen(32)
            .paired(true)
            .quality(false)
            .build(Cursor::new(Vec::new()))?;

        let record = minimal_single_record();
        let result = writer.push(record);
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn test_bq_single_writer_paired_record() -> Result<()> {
        // Writer: single-end
        // Record: paired
        // Expected: success (over-specified - R2 ignored)
        let mut writer = BinseqWriterBuilder::new(Format::Bq)
            .slen(32)
            .paired(false)
            .quality(false)
            .build(Cursor::new(Vec::new()))?;

        let record = minimal_paired_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_bq_paired_minimal_writer_paired_full_record() -> Result<()> {
        // Writer: paired, no quality, no flags
        // Record: paired, with quality, headers, flags
        // Expected: success (over-specified)
        let mut writer = BinseqWriterBuilder::new(Format::Bq)
            .slen(32)
            .xlen(32)
            .paired(true)
            .quality(false)
            .flags(false)
            .build(Cursor::new(Vec::new()))?;

        let record = full_paired_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_bq_paired_full_writer_paired_full_record() -> Result<()> {
        // Writer: paired, with quality, flags
        // Record: paired, with quality, headers, flags
        // Expected: success (correctly specified, headers ignored for BQ)
        let mut writer = BinseqWriterBuilder::new(Format::Bq)
            .slen(32)
            .xlen(32)
            .paired(true)
            .quality(true)
            .flags(true)
            .build(Cursor::new(Vec::new()))?;

        let record = full_paired_record();
        assert!(writer.push(record)?);
        writer.finish()?;
        Ok(())
    }

    // ==================== Configured Size Calculation Tests ====================

    #[test]
    fn test_configured_size_cbq_single_minimal() {
        let record = minimal_single_record();
        // 32 nucleotides = 1 u64 word = 8 bytes
        let size = record.configured_size_cbq(false, false, false, false);
        assert_eq!(size, 8);
    }

    #[test]
    fn test_configured_size_cbq_single_with_flags() {
        let record = full_single_record();
        // 32 nucleotides = 8 bytes + 8 bytes flag
        let size = record.configured_size_cbq(false, true, false, false);
        assert_eq!(size, 16);
    }

    #[test]
    fn test_configured_size_cbq_single_with_all() {
        let record = full_single_record();
        // 32 nucleotides = 8 bytes
        // + 8 bytes flag
        // + 5 bytes header ("read1")
        // + 32 bytes quality
        let size = record.configured_size_cbq(false, true, true, true);
        assert_eq!(size, 8 + 8 + 5 + 32);
    }

    #[test]
    fn test_configured_size_cbq_paired_minimal() {
        let record = full_paired_record();
        // s_seq: 32 nucleotides = 8 bytes
        // x_seq: 32 nucleotides = 8 bytes
        let size = record.configured_size_cbq(true, false, false, false);
        assert_eq!(size, 16);
    }

    #[test]
    fn test_configured_size_cbq_paired_with_all() {
        let record = full_paired_record();
        // s_seq: 32 nucleotides = 8 bytes
        // x_seq: 32 nucleotides = 8 bytes
        // flag: 8 bytes
        // s_header: 5 bytes ("read1")
        // x_header: 5 bytes ("read2")
        // s_qual: 32 bytes
        // x_qual: 32 bytes
        let size = record.configured_size_cbq(true, true, true, true);
        assert_eq!(size, 8 + 8 + 8 + 5 + 5 + 32 + 32);
    }

    #[test]
    fn test_configured_size_cbq_paired_record_single_writer() {
        // A paired record being written to a single-end writer
        // should only count R1 data
        let record = full_paired_record();
        let size = record.configured_size_cbq(false, true, true, true);
        // Only s_seq (8) + flag (8) + s_header (5) + s_qual (32)
        assert_eq!(size, 8 + 8 + 5 + 32);
    }

    #[test]
    fn test_configured_size_vbq_single_minimal() {
        use bitnuc::BitSize;
        let record = minimal_single_record();
        // s_len (8) + x_len (8) + s_seq (32 nucs = 1 word = 8 bytes)
        let size = record.configured_size_vbq(false, false, false, false, BitSize::Two);
        assert_eq!(size, 16 + 8);
    }

    #[test]
    fn test_configured_size_vbq_single_with_flags() {
        use bitnuc::BitSize;
        let record = full_single_record();
        // s_len (8) + x_len (8) + flag (8) + s_seq (8)
        let size = record.configured_size_vbq(false, true, false, false, BitSize::Two);
        assert_eq!(size, 16 + 8 + 8);
    }

    #[test]
    fn test_configured_size_vbq_single_with_all() {
        use bitnuc::BitSize;
        let record = full_single_record();
        // s_len (8) + x_len (8) + flag (8) + s_seq (8) + s_qual (32) + s_header_len (8) + s_header (5)
        let size = record.configured_size_vbq(false, true, true, true, BitSize::Two);
        assert_eq!(size, 16 + 8 + 8 + 32 + 8 + 5);
    }

    #[test]
    fn test_configured_size_vbq_paired_minimal() {
        use bitnuc::BitSize;
        let record = full_paired_record();
        // s_len (8) + x_len (8) + s_seq (8) + x_seq (8)
        let size = record.configured_size_vbq(true, false, false, false, BitSize::Two);
        assert_eq!(size, 16 + 8 + 8);
    }

    #[test]
    fn test_configured_size_vbq_paired_with_all() {
        use bitnuc::BitSize;
        let record = full_paired_record();
        // s_len (8) + x_len (8) + flag (8) + s_seq (8) + x_seq (8)
        // + s_qual (32) + x_qual (32)
        // + s_header_len (8) + s_header (5) + x_header_len (8) + x_header (5)
        let size = record.configured_size_vbq(true, true, true, true, BitSize::Two);
        assert_eq!(size, 16 + 8 + 8 + 8 + 32 + 32 + 8 + 5 + 8 + 5);
    }

    #[test]
    fn test_configured_size_vbq_paired_record_single_writer() {
        use bitnuc::BitSize;
        // A paired record being written to a single-end writer
        // should only count R1 data
        let record = full_paired_record();
        let size = record.configured_size_vbq(false, true, true, true, BitSize::Two);
        // s_len (8) + x_len (8) + flag (8) + s_seq (8) + s_qual (32) + s_header_len (8) + s_header (5)
        assert_eq!(size, 16 + 8 + 8 + 32 + 8 + 5);
    }

    #[test]
    fn test_configured_size_vbq_four_bit_encoding() {
        use bitnuc::BitSize;
        let record = minimal_single_record();
        // With 4-bit encoding: 2 nucleotides per byte, 16 per word
        // 32 nucleotides = 2 words = 16 bytes
        // s_len (8) + x_len (8) + s_seq (16)
        let size = record.configured_size_vbq(false, false, false, false, BitSize::Four);
        assert_eq!(size, 16 + 16);
    }

    // ==================== Multiple Records Tests ====================

    #[test]
    fn test_vbq_multiple_records_mixed_specification() -> Result<()> {
        // Writer configured minimally, records over-specified
        let mut writer = BinseqWriterBuilder::new(Format::Vbq)
            .paired(false)
            .quality(false)
            .headers(false)
            .flags(false)
            .build(Cursor::new(Vec::new()))?;

        // Push minimal record
        assert!(writer.push(minimal_single_record())?);
        // Push full record (over-specified, should work)
        assert!(writer.push(full_single_record())?);
        // Push paired record (over-specified, R2 ignored)
        assert!(writer.push(full_paired_record())?);

        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_cbq_multiple_records_mixed_specification() -> Result<()> {
        // Writer configured minimally, records over-specified
        let mut writer = BinseqWriterBuilder::new(Format::Cbq)
            .paired(false)
            .quality(false)
            .headers(false)
            .flags(false)
            .build(Cursor::new(Vec::new()))?;

        // Push minimal record
        assert!(writer.push(minimal_single_record())?);
        // Push full record (over-specified, should work)
        assert!(writer.push(full_single_record())?);
        // Push paired record (over-specified, R2 ignored)
        assert!(writer.push(full_paired_record())?);

        writer.finish()?;
        Ok(())
    }

    #[test]
    fn test_bq_multiple_records_mixed_specification() -> Result<()> {
        // Writer configured minimally, records over-specified
        let mut writer = BinseqWriterBuilder::new(Format::Bq)
            .slen(32)
            .paired(false)
            .quality(false)
            .flags(false)
            .build(Cursor::new(Vec::new()))?;

        // Push minimal record
        assert!(writer.push(minimal_single_record())?);
        // Push full record (over-specified, should work)
        assert!(writer.push(full_single_record())?);
        // Push paired record (over-specified, R2 ignored)
        assert!(writer.push(full_paired_record())?);

        writer.finish()?;
        Ok(())
    }
}
