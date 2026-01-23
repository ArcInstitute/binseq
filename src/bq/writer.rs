//! Binary sequence writer module
//!
//! This module provides functionality for writing nucleotide sequences to binary files
//! in a compact 2-bit format. It includes support for:
//! - Single and paired sequence writing
//! - Invalid nucleotide handling with configurable policies
//! - Efficient buffering and encoding
//! - Headless mode for parallel writing

use std::io::{BufWriter, Write};

use byteorder::{LittleEndian, WriteBytesExt};
use rand::{SeedableRng, rngs::SmallRng};

use super::FileHeader;
use crate::{
    Policy, RNG_SEED, SequencingRecord,
    error::{Result, WriteError},
};

/// Writes a single flag value to a writer in little-endian format
///
/// # Arguments
///
/// * `writer` - Any type that implements the `Write` trait
/// * `flag` - The 64-bit flag value to write
///
/// # Returns
///
/// * `Ok(())` - If the flag was successfully written
/// * `Err(Error)` - If writing to the writer failed
pub fn write_flag<W: Write>(writer: &mut W, flag: u64) -> Result<()> {
    writer.write_u64::<LittleEndian>(flag)?;
    Ok(())
}

/// Writes a buffer of u64 values to a writer in little-endian format
///
/// This function is used to write encoded sequence data to the output.
/// Each u64 in the buffer contains up to 32 nucleotides in 2-bit format.
///
/// # Arguments
///
/// * `writer` - Any type that implements the `Write` trait
/// * `ebuf` - The buffer of u64 values to write
///
/// # Returns
///
/// * `Ok(())` - If the buffer was successfully written
/// * `Err(Error)` - If writing to the writer failed
pub fn write_buffer<W: Write>(writer: &mut W, ebuf: &[u64]) -> Result<()> {
    ebuf.iter()
        .try_for_each(|&x| writer.write_u64::<LittleEndian>(x))?;
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

    /// Writes a single record to the output
    ///
    /// This method encodes and writes a primary sequence along with an associated flag.
    ///
    /// # Arguments
    ///
    /// * `flag` - A 64-bit flag value associated with the sequence
    /// * `primary` - The nucleotide sequence to write
    ///
    /// # Returns
    ///
    /// * `Ok(true)` if the record was written successfully
    /// * `Ok(false)` if the record was not written because it was empty
    /// * `Err(WriteError::FlagSet)` if the flag is set but no flag value is provided
    #[deprecated]
    pub fn write_record(&mut self, flag: Option<u64>, primary: &[u8]) -> Result<bool> {
        let has_flag = self.encoder.header.flags;
        if let Some(sbuffer) = self.encoder.encode_single(primary)? {
            if has_flag {
                write_flag(&mut self.inner, flag.unwrap_or(0))?;
            }
            write_buffer(&mut self.inner, sbuffer)?;
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Writes a paired record to the output
    ///
    /// This method writes a paired record to the output. It takes a flag, primary sequence, and extended sequence as input.
    /// If the flag is set but no flag value is provided, it returns an error.
    /// Otherwise, it writes the encoded single and extended sequences to the output and returns true.
    ///
    /// # Arguments
    /// * `flag` - The flag value to write to the output
    /// * `primary` - The primary sequence to encode and write to the output
    /// * `extended` - The extended sequence to encode and write to the output
    ///
    /// # Returns
    /// * `Result<bool>` - A result indicating whether the write was successful or not
    #[deprecated]
    pub fn write_paired_record(
        &mut self,
        flag: Option<u64>,
        primary: &[u8],
        extended: &[u8],
    ) -> Result<bool> {
        let has_flag = self.encoder.header.flags;
        if let Some((sbuffer, xbuffer)) = self.encoder.encode_paired(primary, extended)? {
            if has_flag {
                write_flag(&mut self.inner, flag.unwrap_or(0))?;
            }
            write_buffer(&mut self.inner, sbuffer)?;
            write_buffer(&mut self.inner, xbuffer)?;
            Ok(true)
        } else {
            Ok(false)
        }
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
            write_flag(&mut self.inner, record.flag().unwrap_or(0))?;
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

    /// Creates a new encoder with the same configuration as this writer
    ///
    /// This is useful when you need a separate encoder instance for parallel
    /// processing or other scenarios where you need independent encoding.
    /// The new encoder is initialized with a cleared state.
    ///
    /// # Returns
    ///
    /// A new `Encoder` instance with the same configuration but cleared buffers
    pub fn new_encoder(&self) -> Encoder {
        let mut encoder = self.encoder.clone();
        encoder.clear();
        encoder
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

/// A streaming writer for binary sequence data
///
/// This writer buffers data before writing it to the underlying writer,
/// providing efficient streaming capabilities suitable for:
/// - Writing to network connections
/// - Processing very large datasets
/// - Pipeline processing
///
/// The `StreamWriter` is a specialized version of `Writer` that
/// adds internal buffering and is optimized for streaming scenarios.
pub struct StreamWriter<W: Write> {
    /// The underlying writer for processing sequences
    writer: Writer<BufWriter<W>>,
}

impl<W: Write> StreamWriter<W> {
    /// Creates a new `StreamWriter` with the default buffer size
    ///
    /// This constructor initializes a `StreamWriter` with an 8K buffer
    /// for efficient writing to the underlying writer.
    ///
    /// # Arguments
    ///
    /// * `inner` - The writer to write binary sequence data to
    /// * `header` - The header defining sequence lengths and format
    /// * `policy` - The policy for handling invalid nucleotides
    /// * `headless` - Whether to skip writing the header
    ///
    /// # Returns
    ///
    /// * `Ok(StreamWriter)` - A new streaming writer
    /// * `Err(Error)` - If initialization fails
    pub fn new(inner: W, header: FileHeader, policy: Policy, headless: bool) -> Result<Self> {
        Self::with_capacity(inner, 8192, header, policy, headless)
    }

    /// Creates a new `StreamWriter` with a specified buffer capacity
    ///
    /// This constructor allows customizing the buffer size based on
    /// expected usage patterns and performance requirements.
    ///
    /// # Arguments
    ///
    /// * `inner` - The writer to write binary sequence data to
    /// * `capacity` - The size of the internal buffer in bytes
    /// * `header` - The header defining sequence lengths and format
    /// * `policy` - The policy for handling invalid nucleotides
    /// * `headless` - Whether to skip writing the header
    ///
    /// # Returns
    ///
    /// * `Ok(StreamWriter)` - A new streaming writer with the specified buffer capacity
    /// * `Err(Error)` - If initialization fails
    pub fn with_capacity(
        inner: W,
        capacity: usize,
        header: FileHeader,
        policy: Policy,
        headless: bool,
    ) -> Result<Self> {
        let buffered = BufWriter::with_capacity(capacity, inner);
        let writer = Writer::new(buffered, header, policy, headless)?;

        Ok(Self { writer })
    }

    #[deprecated(note = "use `push` method with SequencingRecord instead")]
    pub fn write_record(&mut self, flag: Option<u64>, primary: &[u8]) -> Result<bool> {
        #[allow(deprecated)]
        self.writer.write_record(flag, primary)
    }

    #[deprecated(note = "use `push` method with SequencingRecord instead")]
    pub fn write_paired_record(
        &mut self,
        flag: Option<u64>,
        primary: &[u8],
        extended: &[u8],
    ) -> Result<bool> {
        #[allow(deprecated)]
        self.writer.write_paired_record(flag, primary, extended)
    }

    /// Writes a record using the unified [`SequencingRecord`] API
    pub fn push(&mut self, record: SequencingRecord) -> Result<bool> {
        self.writer.push(record)
    }

    /// Flushes any buffered data to the underlying writer
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If the flush was successful
    /// * `Err(Error)` - If flushing failed
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()
    }

    /// Consumes the streaming writer and returns the inner writer after flushing
    ///
    /// This method is useful when you need access to the underlying writer
    /// after all writing is complete.
    ///
    /// # Returns
    ///
    /// * `Ok(W)` - The inner writer after flushing all data
    /// * `Err(Error)` - If flushing failed
    pub fn into_inner(self) -> Result<W> {
        // First unwrap the writer inner (BufWriter<W>)
        let bufw = self.writer.into_inner();
        // Now unwrap the BufWriter to get W
        match bufw.into_inner() {
            Ok(inner) => Ok(inner),
            Err(e) => Err(std::io::Error::from(e).into()),
        }
    }
}

/// Builder for `StreamWriter` instances
///
/// This builder provides a convenient way to create and configure `StreamWriter`
/// instances with custom buffer sizes and other settings.
#[derive(Default)]
pub struct StreamWriterBuilder {
    /// Required header defining sequence lengths and format
    header: Option<FileHeader>,
    /// Optional policy for handling invalid nucleotides
    policy: Option<Policy>,
    /// Optional headless mode for parallel writing scenarios
    headless: Option<bool>,
    /// Optional buffer capacity setting
    buffer_capacity: Option<usize>,
}

impl StreamWriterBuilder {
    /// Sets the header for the writer
    #[must_use]
    pub fn header(mut self, header: FileHeader) -> Self {
        self.header = Some(header);
        self
    }

    /// Sets the policy for handling invalid nucleotides
    #[must_use]
    pub fn policy(mut self, policy: Policy) -> Self {
        self.policy = Some(policy);
        self
    }

    /// Sets headless mode (whether to skip writing the header)
    #[must_use]
    pub fn headless(mut self, headless: bool) -> Self {
        self.headless = Some(headless);
        self
    }

    /// Sets the buffer capacity for the writer
    #[must_use]
    pub fn buffer_capacity(mut self, capacity: usize) -> Self {
        self.buffer_capacity = Some(capacity);
        self
    }

    /// Builds a `StreamWriter` with the configured settings
    ///
    /// # Arguments
    ///
    /// * `inner` - The writer to write binary sequence data to
    ///
    /// # Returns
    ///
    /// * `Ok(StreamWriter)` - A new streaming writer with the specified configuration
    /// * `Err(Error)` - If building the writer fails
    pub fn build<W: Write>(self, inner: W) -> Result<StreamWriter<W>> {
        let Some(header) = self.header else {
            return Err(WriteError::MissingHeader.into());
        };

        let capacity = self.buffer_capacity.unwrap_or(8192);
        StreamWriter::with_capacity(
            inner,
            capacity,
            header,
            self.policy.unwrap_or_default(),
            self.headless.unwrap_or(false),
        )
    }
}

#[cfg(test)]
mod testing {

    use std::{fs::File, io::BufWriter};

    use super::*;
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

    #[test]
    fn test_stream_writer() -> Result<()> {
        let inner = Vec::new();
        let writer = StreamWriterBuilder::default()
            .header(FileHeaderBuilder::new().slen(32).build()?)
            .buffer_capacity(16384)
            .build(inner)?;

        // Convert back to Vec to verify it works
        let inner = writer.into_inner()?;
        assert_eq!(inner.len(), SIZE_HEADER);
        Ok(())
    }
}
