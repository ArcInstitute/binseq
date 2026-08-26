//! Binary sequence writer module
//!
//! This module provides functionality for writing nucleotide sequences to binary files
//! in a compact 2-bit format. It includes support for:
//! - Single and paired sequence writing
//! - Invalid nucleotide handling with configurable policies
//! - Efficient buffering and encoding
//! - Headless mode for parallel writing

use std::io::Write;

use rand::{SeedableRng, rngs::SmallRng};

use super::FileHeader;
use crate::{
    Policy, RNG_SEED, SequencingRecord,
    error::{Result, WriteError},
};

/// Writes a buffer of u64 values to a writer in little-endian format
fn write_buffer<W: Write>(writer: &mut W, ebuf: &[u64]) -> Result<()> {
    ebuf.iter()
        .try_for_each(|&x| writer.write_all(&x.to_le_bytes()))?;
    Ok(())
}

/// Encodes nucleotide sequences into a compact 2-bit binary format
///
/// The `Encoder` handles the conversion of nucleotide sequences (A, C, G, T)
/// into a compact binary representation where each nucleotide is stored using
/// 2 bits. It also handles invalid nucleotides according to a configurable policy.
///
/// The encoder maintains internal buffers to avoid repeated allocations during
/// encoding operations. These buffers are reused across multiple encode calls
/// and are cleared automatically when needed.
#[derive(Clone)]
pub struct Encoder {
    /// Header containing sequence length and format information
    header: FileHeader,

    /// Buffers for storing encoded nucleotides in 2-bit format
    /// Each u64 can store 32 nucleotides (64 bits / 2 bits per nucleotide)
    sbuffer: Vec<u64>, // Primary sequence buffer
    xbuffer: Vec<u64>, // Extended sequence buffer

    /// Temporary buffers for handling invalid nucleotides
    /// These store the processed sequences after policy application
    s_ibuf: Vec<u8>, // Primary sequence invalid buffer
    x_ibuf: Vec<u8>, // Extended sequence invalid buffer

    /// Policy for handling invalid nucleotides during encoding
    policy: Policy,

    /// Random number generator for the `RandomDraw` policy
    /// Seeded with `RNG_SEED` for reproducibility
    rng: SmallRng,
}
impl Encoder {
    /// Creates a new encoder with default invalid nucleotide policy
    ///
    /// # Arguments
    ///
    /// * `header` - The header defining sequence lengths and format
    ///
    /// # Examples
    ///
    /// ```
    /// # use binseq::bq::{FileHeaderBuilder, Encoder};
    /// let header = FileHeaderBuilder::new().slen(100).build().unwrap();
    /// let encoder = Encoder::new(header);
    /// ```
    #[must_use]
    pub fn new(header: FileHeader) -> Self {
        Self::with_policy(header, Policy::default())
    }

    /// Creates a new encoder with a specific invalid nucleotide policy
    ///
    /// # Arguments
    ///
    /// * `header` - The header defining sequence lengths and format
    /// * `policy` - The policy for handling invalid nucleotides
    ///
    /// # Examples
    ///
    /// ```
    /// # use binseq::bq::{FileHeaderBuilder, Encoder};
    /// # use binseq::Policy;
    /// let header = FileHeaderBuilder::new().slen(100).build().unwrap();
    /// let encoder = Encoder::with_policy(header, Policy::SetToA);
    /// ```
    #[must_use]
    pub fn with_policy(header: FileHeader, policy: Policy) -> Self {
        Self {
            header,
            policy,
            sbuffer: Vec::default(),
            xbuffer: Vec::default(),
            s_ibuf: Vec::default(),
            x_ibuf: Vec::default(),
            rng: SmallRng::seed_from_u64(RNG_SEED),
        }
    }

    /// Returns whether the header is paired-end.
    #[must_use]
    pub fn is_paired(&self) -> bool {
        self.header.is_paired()
    }

    /// Encodes a single sequence as 2-bit.
    ///
    /// Will return `None` if the sequence is invalid and the policy does not allow correction.
    pub fn encode_single(&mut self, primary: &[u8]) -> Result<Option<&[u64]>> {
        if primary.len() != self.header.slen as usize {
            return Err(WriteError::UnexpectedSequenceLength {
                expected: self.header.slen,
                got: primary.len(),
            }
            .into());
        }

        // Fill the buffer with the 2-bit representation of the nucleotides
        self.clear();
        if self.header.bits.encode(primary, &mut self.sbuffer).is_err() {
            self.clear();
            if self
                .policy
                .handle(primary, &mut self.s_ibuf, &mut self.rng)?
            {
                self.header.bits.encode(&self.s_ibuf, &mut self.sbuffer)?;
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
        if primary.len() != self.header.slen as usize {
            return Err(WriteError::UnexpectedSequenceLength {
                expected: self.header.slen,
                got: primary.len(),
            }
            .into());
        }
        if extended.len() != self.header.xlen as usize {
            return Err(WriteError::UnexpectedSequenceLength {
                expected: self.header.xlen,
                got: extended.len(),
            }
            .into());
        }

        self.clear();
        if self.header.bits.encode(primary, &mut self.sbuffer).is_err()
            || self
                .header
                .bits
                .encode(extended, &mut self.xbuffer)
                .is_err()
        {
            self.clear();
            if self
                .policy
                .handle(primary, &mut self.s_ibuf, &mut self.rng)?
                && self
                    .policy
                    .handle(extended, &mut self.x_ibuf, &mut self.rng)?
            {
                self.header.bits.encode(&self.s_ibuf, &mut self.sbuffer)?;
                self.header.bits.encode(&self.x_ibuf, &mut self.xbuffer)?;
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
            encoder: Encoder::with_policy(header, policy),
            headless,
        })
    }

    /// Returns whether the header is paired-end.
    pub fn is_paired(&self) -> bool {
        self.encoder.is_paired()
    }

    /// Returns the header of the writer
    pub fn header(&self) -> FileHeader {
        self.encoder.header
    }

    /// Returns the N-policy of the writer
    pub fn policy(&self) -> Policy {
        self.encoder.policy
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
        let has_flag = self.encoder.header.flags;
        if has_flag {
            self.inner
                .write_all(&record.flag.unwrap_or(0).to_le_bytes())?;
        }

        // Check paired status - writer can require paired (record must have R2),
        // but if writer is single-end, we simply ignore any R2 data in the record.
        if self.encoder.header.is_paired() && !record.is_paired() {
            return Err(WriteError::ConfigurationMismatch {
                attribute: "paired",
                expected: self.encoder.header.is_paired(),
                actual: record.is_paired(),
            }
            .into());
        }

        if self.encoder.header.is_paired() {
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

    // ==================== Encoder Tests ====================

    #[test]
    fn test_encoder_new() {
        let header = FileHeaderBuilder::new().slen(8).build().unwrap();
        let encoder = Encoder::new(header);
        assert!(matches!(encoder.policy, Policy::IgnoreSequence));
    }

    #[test]
    fn test_encoder_encode_single_wrong_length() {
        let header = FileHeaderBuilder::new().slen(8).build().unwrap();
        let mut encoder = Encoder::new(header);
        let result = encoder.encode_single(b"ACGT");
        assert!(result.is_err());
    }

    #[test]
    fn test_encoder_encode_single_invalid_ignored() {
        let header = FileHeaderBuilder::new().slen(8).build().unwrap();
        let mut encoder = Encoder::with_policy(header, Policy::IgnoreSequence);
        let result = encoder.encode_single(b"ACGTNNNN").unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_encoder_encode_single_invalid_corrected() {
        let header = FileHeaderBuilder::new().slen(8).build().unwrap();
        let mut encoder = Encoder::with_policy(header, Policy::SetToA);
        let result = encoder.encode_single(b"ACGTNNNN").unwrap();
        assert!(result.is_some());
    }

    #[test]
    fn test_encoder_encode_paired_wrong_primary_length() {
        let header = FileHeaderBuilder::new().slen(8).xlen(8).build().unwrap();
        let mut encoder = Encoder::new(header);
        let result = encoder.encode_paired(b"ACGT", b"ACGTACGT");
        assert!(result.is_err());
    }

    #[test]
    fn test_encoder_encode_paired_wrong_extended_length() {
        let header = FileHeaderBuilder::new().slen(8).xlen(8).build().unwrap();
        let mut encoder = Encoder::new(header);
        let result = encoder.encode_paired(b"ACGTACGT", b"ACGT");
        assert!(result.is_err());
    }

    #[test]
    fn test_encoder_encode_paired_invalid_ignored() {
        let header = FileHeaderBuilder::new().slen(8).xlen(8).build().unwrap();
        let mut encoder = Encoder::with_policy(header, Policy::IgnoreSequence);
        let result = encoder.encode_paired(b"ACGTNNNN", b"ACGTACGT").unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_encoder_encode_paired_invalid_corrected() {
        let header = FileHeaderBuilder::new().slen(8).xlen(8).build().unwrap();
        let mut encoder = Encoder::with_policy(header, Policy::SetToA);
        let result = encoder.encode_paired(b"ACGTNNNN", b"NNNNACGT").unwrap();
        assert!(result.is_some());
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
