//! Binary sequence writer module
//!
//! This module provides functionality for writing nucleotide sequences to binary files
//! in a compact 2-bit format. It includes support for:
//! - Single and paired sequence writing
//! - Invalid nucleotide handling with configurable policies
//! - Efficient buffering and encoding
//! - Headless mode for parallel writing

use std::io::Write;

use super::FileHeader;
use crate::{
    Policy, SequencingRecord,
    encoder::Encoder,
    error::{Result, WriteError},
};

/// Writes a buffer of u64 values to a writer in little-endian format
fn write_buffer<W: Write>(writer: &mut W, ebuf: &[u64]) -> Result<()> {
    ebuf.iter()
        .try_for_each(|&x| writer.write_all(&x.to_le_bytes()))?;
    Ok(())
}

/// Builder for creating configured `Writer` instances
///
/// This builder provides a flexible way to create writers with various
/// configurations. It follows the builder pattern, allowing for optional
/// settings to be specified in any order.
///
/// # Examples
///
/// ```
/// # use binseq::{Policy, Result};
/// # use binseq::bq::{FileHeaderBuilder, WriterBuilder};
/// # fn main() -> Result<()> {
/// let header = FileHeaderBuilder::new().slen(100).build()?;
/// let writer = WriterBuilder::default()
///     .header(header)
///     .policy(Policy::SetToA)
///     .headless(false)
///     .build(Vec::new())?;
/// # Ok(())
/// # }
/// ```
#[derive(Default)]
pub struct WriterBuilder {
    /// Required header defining sequence lengths and format
    header: Option<FileHeader>,
    /// Optional policy for handling invalid nucleotides
    policy: Option<Policy>,
    /// Optional headless mode for parallel writing scenarios
    headless: Option<bool>,
}
impl WriterBuilder {
    #[must_use]
    pub fn header(mut self, header: FileHeader) -> Self {
        self.header = Some(header);
        self
    }

    #[must_use]
    pub fn policy(mut self, policy: Policy) -> Self {
        self.policy = Some(policy);
        self
    }

    #[must_use]
    pub fn headless(mut self, headless: bool) -> Self {
        self.headless = Some(headless);
        self
    }

    pub fn build<W: Write>(self, inner: W) -> Result<Writer<W>> {
        let Some(header) = self.header else {
            return Err(WriteError::MissingHeader.into());
        };
        Writer::new(
            inner,
            header,
            self.policy.unwrap_or_default(),
            self.headless.unwrap_or(false),
        )
    }
}

/// High-level writer for binary sequence files
///
/// This writer provides a convenient interface for writing nucleotide sequences
/// to binary files in a compact format. It handles sequence encoding, invalid
/// nucleotide processing, and file format compliance.
///
/// The writer can operate in two modes:
/// - Normal mode: Writes the header followed by records
/// - Headless mode: Writes only records (useful for parallel writing)
///
/// # Type Parameters
///
/// * `W` - The underlying writer type that implements `Write`
#[derive(Clone)]
pub struct Writer<W: Write> {
    /// The underlying writer for output
    inner: W,

    /// Header defining sequence lengths and format
    header: FileHeader,

    /// Encoder for converting sequences to binary format
    encoder: Encoder,

    /// Whether this writer is in headless mode
    /// When true, the header is not written to the output
    headless: bool,
}
impl<W: Write> Writer<W> {
    /// Creates a new `Writer` instance with specified configuration
    ///
    /// This is a low-level constructor. For a more convenient way to create a
    /// `Writer`, use the `WriterBuilder` struct.
    ///
    /// # Arguments
    ///
    /// * `inner` - The underlying writer to write to
    /// * `header` - The header defining sequence lengths and format
    /// * `policy` - The policy for handling invalid nucleotides
    /// * `headless` - Whether to skip writing the header (for parallel writing)
    ///
    /// # Returns
    ///
    /// * `Ok(Writer)` - A new writer instance
    /// * `Err(Error)` - If writing the header fails
    ///
    /// # Examples
    ///
    /// ```
    /// # use binseq::bq::{FileHeaderBuilder, Writer};
    /// # use binseq::{Result, Policy};
    /// # fn main() -> Result<()> {
    /// let header = FileHeaderBuilder::new().slen(100).build()?;
    /// let writer = Writer::new(
    ///     Vec::new(),
    ///     header,
    ///     Policy::default(),
    ///     false
    /// )?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(mut inner: W, header: FileHeader, policy: Policy, headless: bool) -> Result<Self> {
        if !headless {
            header.write_bytes(&mut inner)?;
        }
        Ok(Self {
            inner,
            header,
            encoder: Encoder::with_policy(header.bits, policy),
            headless,
        })
    }

    /// Returns whether the header is paired-end.
    pub fn is_paired(&self) -> bool {
        self.header.is_paired()
    }

    /// Returns the header of the writer
    pub fn header(&self) -> FileHeader {
        self.header
    }

    /// Returns the N-policy of the writer
    pub fn policy(&self) -> Policy {
        self.encoder.policy()
    }

    /// Writes a record using the unified [`SequencingRecord`] API
    ///
    /// This method provides a consistent interface with VBQ and CBQ writers.
    /// Note that BQ format does not support quality scores or headers - these
    /// fields from the record will be ignored.
    ///
    /// # Arguments
    ///
    /// * `record` - A [`SequencingRecord`] containing the sequence data to write
    ///
    /// # Returns
    ///
    /// * `Ok(true)` if the record was written successfully
    /// * `Ok(false)` if the record was skipped due to invalid nucleotides
    /// * `Err(_)` if writing failed
    ///
    /// # Examples
    ///
    /// ```
    /// # use binseq::bq::{FileHeaderBuilder, WriterBuilder};
    /// # use binseq::{Result, SequencingRecordBuilder};
    /// # fn main() -> Result<()> {
    /// let header = FileHeaderBuilder::new().slen(8).build()?;
    /// let mut writer = WriterBuilder::default()
    ///     .header(header)
    ///     .build(Vec::new())?;
    ///
    /// let record = SequencingRecordBuilder::default()
    ///     .s_seq(b"ACGTACGT")
    ///     .flag(42)
    ///     .build()?;
    ///
    /// writer.push(record)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn push(&mut self, record: SequencingRecord) -> Result<bool> {
        let has_flag = self.header.flags;
        if has_flag {
            self.inner
                .write_all(&record.flag.unwrap_or(0).to_le_bytes())?;
        }

        // Check paired status - writer can require paired (record must have R2),
        // but if writer is single-end, we simply ignore any R2 data in the record.
        if self.header.is_paired() && !record.is_paired() {
            return Err(WriteError::ConfigurationMismatch {
                attribute: "paired",
                expected: self.header.is_paired(),
                actual: record.is_paired(),
            }
            .into());
        }

        // BQ records are fixed-length: validate against the header
        if record.s_seq.len() != self.header.slen as usize {
            return Err(WriteError::UnexpectedSequenceLength {
                expected: self.header.slen,
                got: record.s_seq.len(),
            }
            .into());
        }
        if self.header.is_paired() {
            let xlen = record.x_seq.unwrap_or_default().len();
            if xlen != self.header.xlen as usize {
                return Err(WriteError::UnexpectedSequenceLength {
                    expected: self.header.xlen,
                    got: xlen,
                }
                .into());
            }
        }

        if self.header.is_paired() {
            if let Some((sbuffer, xbuffer)) = self
                .encoder
                .encode_paired(record.s_seq, record.x_seq.unwrap_or_default())?
            {
                write_buffer(&mut self.inner, sbuffer)?;
                write_buffer(&mut self.inner, xbuffer)?;
                Ok(true)
            } else {
                Ok(false)
            }
        } else if let Some(buffer) = self.encoder.encode_single(record.s_seq)? {
            write_buffer(&mut self.inner, buffer)?;
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Consumes the writer and returns the underlying writer
    ///
    /// This is useful when you need to access the underlying writer after
    /// writing is complete, for example to get the contents of a `Vec<u8>`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use binseq::bq::{FileHeaderBuilder, WriterBuilder};
    /// # use binseq::Result;
    /// # fn main() -> Result<()> {
    /// let header = FileHeaderBuilder::new().slen(100).build()?;
    /// let writer = WriterBuilder::default()
    ///     .header(header)
    ///     .build(Vec::new())?;
    ///
    /// // After writing sequences...
    /// let bytes = writer.into_inner();
    /// # Ok(())
    /// # }
    /// ```
    pub fn into_inner(self) -> W {
        self.inner
    }

    /// Gets a mutable reference to the underlying writer
    ///
    /// This allows direct access to the underlying writer while retaining
    /// ownership of the `Writer`.
    pub fn by_ref(&mut self) -> &mut W {
        &mut self.inner
    }

    /// Flushes any buffered data to the underlying writer
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If the flush was successful
    /// * `Err(Error)` - If flushing failed
    pub fn flush(&mut self) -> Result<()> {
        self.inner.flush()?;
        Ok(())
    }

    /// Checks if this writer is in headless mode
    ///
    /// In headless mode, the writer does not write the header to the output.
    /// This is useful for parallel writing scenarios where only one writer
    /// should write the header.
    ///
    /// # Returns
    ///
    /// `true` if the writer is in headless mode, `false` otherwise
    pub fn is_headless(&self) -> bool {
        self.headless
    }

    /// Ingests the contents of another writer's buffer
    ///
    /// This method is used in parallel writing scenarios to combine the output
    /// of multiple writers. It takes the contents of another writer's buffer
    /// and writes them to this writer's output.
    ///
    /// # Arguments
    ///
    /// * `other` - Another writer whose underlying writer is a `Vec<u8>`
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If the contents were successfully ingested
    /// * `Err(Error)` - If writing the contents failed
    pub fn ingest(&mut self, other: &mut Writer<Vec<u8>>) -> Result<()> {
        let other_inner = other.by_ref();
        self.inner.write_all(other_inner)?;
        other_inner.clear();
        Ok(())
    }
}

#[cfg(test)]
mod testing {

    use std::{fs::File, io::BufWriter};

    use super::*;
    use crate::SequencingRecordBuilder;
    use crate::bq::{FileHeaderBuilder, SIZE_HEADER};

    #[test]
    fn test_headless() -> Result<()> {
        let inner = Vec::new();
        let mut writer = WriterBuilder::default()
            .header(FileHeaderBuilder::new().slen(32).build()?)
            .headless(true)
            .build(inner)?;
        assert!(writer.is_headless());
        let inner = writer.by_ref();
        assert!(inner.is_empty());
        Ok(())
    }

    #[test]
    fn test_not_headless() -> Result<()> {
        let inner = Vec::new();
        let mut writer = WriterBuilder::default()
            .header(FileHeaderBuilder::new().slen(32).build()?)
            .build(inner)?;
        assert!(!writer.is_headless());
        let inner = writer.by_ref();
        assert_eq!(inner.len(), SIZE_HEADER);
        Ok(())
    }

    #[test]
    fn test_stdout() -> Result<()> {
        let writer = WriterBuilder::default()
            .header(FileHeaderBuilder::new().slen(32).build()?)
            .build(std::io::stdout())?;
        assert!(!writer.is_headless());
        Ok(())
    }

    #[test]
    fn test_to_path() -> Result<()> {
        let path = "test_to_path.file";
        let inner = File::create(path).map(BufWriter::new)?;
        let mut writer = WriterBuilder::default()
            .header(FileHeaderBuilder::new().slen(32).build()?)
            .build(inner)?;
        assert!(!writer.is_headless());
        let inner = writer.by_ref();
        inner.flush()?;

        // delete file
        std::fs::remove_file(path)?;

        Ok(())
    }

    // ==================== Length Validation Tests ====================

    #[test]
    fn test_push_wrong_length() -> Result<()> {
        let mut writer = WriterBuilder::default()
            .header(FileHeaderBuilder::new().slen(8).build()?)
            .build(Vec::new())?;
        let record = SequencingRecordBuilder::default().s_seq(b"ACGT").build()?;
        assert!(writer.push(record).is_err());
        Ok(())
    }

    #[test]
    fn test_push_wrong_extended_length() -> Result<()> {
        let mut writer = WriterBuilder::default()
            .header(FileHeaderBuilder::new().slen(8).xlen(8).build()?)
            .build(Vec::new())?;
        let record = SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGT")
            .x_seq(b"ACGT")
            .build()?;
        assert!(writer.push(record).is_err());
        Ok(())
    }

    // ==================== WriterBuilder Tests ====================

    #[test]
    fn test_writer_builder_missing_header() {
        let result = WriterBuilder::default().build(Vec::new());
        assert!(result.is_err());
    }

    // ==================== push() Tests ====================

    #[test]
    fn test_push_single() -> Result<()> {
        let mut writer = WriterBuilder::default()
            .header(FileHeaderBuilder::new().slen(8).build()?)
            .build(Vec::new())?;
        let record = SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGT")
            .build()?;
        assert!(writer.push(record)?);
        Ok(())
    }

    #[test]
    fn test_push_single_invalid_skipped() -> Result<()> {
        let mut writer = WriterBuilder::default()
            .header(FileHeaderBuilder::new().slen(8).build()?)
            .build(Vec::new())?;
        let record = SequencingRecordBuilder::default()
            .s_seq(b"NNNNNNNN")
            .build()?;
        assert!(!writer.push(record)?);
        Ok(())
    }

    #[test]
    fn test_push_paired() -> Result<()> {
        let mut writer = WriterBuilder::default()
            .header(FileHeaderBuilder::new().slen(8).xlen(8).build()?)
            .build(Vec::new())?;
        let record = SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGT")
            .x_seq(b"TTGGCCAA")
            .build()?;
        assert!(writer.push(record)?);
        Ok(())
    }

    #[test]
    fn test_push_paired_invalid_skipped() -> Result<()> {
        let mut writer = WriterBuilder::default()
            .header(FileHeaderBuilder::new().slen(8).xlen(8).build()?)
            .build(Vec::new())?;
        let record = SequencingRecordBuilder::default()
            .s_seq(b"NNNNNNNN")
            .x_seq(b"TTGGCCAA")
            .build()?;
        assert!(!writer.push(record)?);
        Ok(())
    }

    #[test]
    fn test_push_paired_mismatch() -> Result<()> {
        let mut writer = WriterBuilder::default()
            .header(FileHeaderBuilder::new().slen(8).xlen(8).build()?)
            .build(Vec::new())?;
        // Writer expects paired records but record has no x_seq
        let record = SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGT")
            .build()?;
        let result = writer.push(record);
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn test_push_with_flag() -> Result<()> {
        let mut writer = WriterBuilder::default()
            .header(FileHeaderBuilder::new().slen(8).flags(true).build()?)
            .build(Vec::new())?;
        let record = SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGT")
            .flag(99)
            .build()?;
        assert!(writer.push(record)?);
        Ok(())
    }

    // ==================== Writer Misc Tests ====================

    #[test]
    fn test_writer_flush() -> Result<()> {
        let mut writer = WriterBuilder::default()
            .header(FileHeaderBuilder::new().slen(8).build()?)
            .build(Vec::new())?;
        writer.flush()?;
        Ok(())
    }

    #[test]
    fn test_writer_ingest() -> Result<()> {
        let mut main_writer = WriterBuilder::default()
            .header(FileHeaderBuilder::new().slen(8).build()?)
            .build(Vec::new())?;

        let header = main_writer.header();
        let mut other_writer = WriterBuilder::default()
            .header(header)
            .headless(true)
            .build(Vec::new())?;

        let record = SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGT")
            .build()?;
        other_writer.push(record)?;

        main_writer.ingest(&mut other_writer)?;
        assert!(other_writer.by_ref().is_empty());
        assert_eq!(main_writer.by_ref().len(), SIZE_HEADER + 8);
        Ok(())
    }

    #[test]
    fn test_writer_policy_accessor() -> Result<()> {
        let writer = WriterBuilder::default()
            .header(FileHeaderBuilder::new().slen(8).build()?)
            .policy(Policy::SetToA)
            .build(Vec::new())?;
        assert!(matches!(writer.policy(), Policy::SetToA));
        Ok(())
    }
}
