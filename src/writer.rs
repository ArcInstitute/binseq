//! Binary sequence writer module
//!
//! This module provides functionality for writing nucleotide sequences to binary files
//! in a compact 2-bit format. It includes support for:
//! - Single and paired sequence writing
//! - Invalid nucleotide handling with configurable policies
//! - Efficient buffering and encoding
//! - Headless mode for parallel writing

use std::io::Write;

use byteorder::{LittleEndian, WriteBytesExt};
use rand::{rngs::SmallRng, SeedableRng};

use crate::{error::WriteError, BinseqHeader, Policy, Result, RNG_SEED};

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
    header: BinseqHeader,

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

    /// Random number generator for the RandomDraw policy
    /// Seeded with RNG_SEED for reproducibility
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
    /// # use binseq::{BinseqHeader, Encoder};
    /// let header = BinseqHeader::new(100); // For sequences of length 100
    /// let encoder = Encoder::new(header);
    /// ```
    pub fn new(header: BinseqHeader) -> Self {
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
    /// # use binseq::{BinseqHeader, Encoder, Policy};
    /// let header = BinseqHeader::new(100);
    /// let encoder = Encoder::with_policy(header, Policy::SetToA);
    /// ```
    pub fn with_policy(header: BinseqHeader, policy: Policy) -> Self {
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
        if crate::nuc::encode(primary, &mut self.sbuffer).is_err() {
            self.clear();
            if self
                .policy
                .handle(primary, &mut self.s_ibuf, &mut self.rng)?
            {
                crate::nuc::encode(&self.s_ibuf, &mut self.sbuffer)?;
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
        if crate::nuc::encode(primary, &mut self.sbuffer).is_err()
            || crate::nuc::encode(extended, &mut self.xbuffer).is_err()
        {
            self.clear();
            if self
                .policy
                .handle(primary, &mut self.s_ibuf, &mut self.rng)?
                && self
                    .policy
                    .handle(extended, &mut self.x_ibuf, &mut self.rng)?
            {
                crate::nuc::encode(&self.s_ibuf, &mut self.sbuffer)?;
                crate::nuc::encode(&self.x_ibuf, &mut self.xbuffer)?;
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

/// Builder for creating configured `BinseqWriter` instances
///
/// This builder provides a flexible way to create writers with various
/// configurations. It follows the builder pattern, allowing for optional
/// settings to be specified in any order.
///
/// # Examples
///
/// ```
/// # use binseq::{BinseqHeader, BinseqWriterBuilder, Policy, Result};
/// # fn main() -> Result<()> {
/// let header = BinseqHeader::new(100);
/// let writer = BinseqWriterBuilder::default()
///     .header(header)
///     .policy(Policy::SetToA)
///     .headless(false)
///     .build(Vec::new())?;
/// # Ok(())
/// # }
/// ```
#[derive(Default)]
pub struct BinseqWriterBuilder {
    /// Required header defining sequence lengths and format
    header: Option<BinseqHeader>,
    /// Optional policy for handling invalid nucleotides
    policy: Option<Policy>,
    /// Optional headless mode for parallel writing scenarios
    headless: Option<bool>,
}
impl BinseqWriterBuilder {
    pub fn header(mut self, header: BinseqHeader) -> Self {
        self.header = Some(header);
        self
    }

    pub fn policy(mut self, policy: Policy) -> Self {
        self.policy = Some(policy);
        self
    }

    pub fn headless(mut self, headless: bool) -> Self {
        self.headless = Some(headless);
        self
    }

    pub fn build<W: Write>(self, inner: W) -> Result<BinseqWriter<W>> {
        let Some(header) = self.header else {
            return Err(WriteError::MissingHeader.into());
        };
        BinseqWriter::new(
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
pub struct BinseqWriter<W: Write> {
    /// The underlying writer for output
    inner: W,

    /// Encoder for converting sequences to binary format
    encoder: Encoder,

    /// Whether this writer is in headless mode
    /// When true, the header is not written to the output
    headless: bool,
}
impl<W: Write> BinseqWriter<W> {
    /// Creates a new `BinseqWriter` instance with specified configuration
    ///
    /// This is a low-level constructor. For a more convenient way to create a
    /// `BinseqWriter`, use the `BinseqWriterBuilder` struct.
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
    /// * `Ok(BinseqWriter)` - A new writer instance
    /// * `Err(Error)` - If writing the header fails
    ///
    /// # Examples
    ///
    /// ```
    /// # use binseq::{BinseqHeader, BinseqWriter, Policy, Result};
    /// # fn main() -> Result<()> {
    /// let header = BinseqHeader::new(100);
    /// let writer = BinseqWriter::new(
    ///     Vec::new(),
    ///     header,
    ///     Policy::default(),
    ///     false
    /// )?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn new(mut inner: W, header: BinseqHeader, policy: Policy, headless: bool) -> Result<Self> {
        if !headless {
            header.write_bytes(&mut inner)?;
        }
        Ok(Self {
            inner,
            encoder: Encoder::with_policy(header, policy),
            headless,
        })
    }

    /// Writes a single nucleotide sequence to the output
    ///
    /// This method encodes and writes a primary sequence along with an associated flag.
    /// The sequence must match the length specified in the header.
    ///
    /// # Arguments
    ///
    /// * `flag` - A 64-bit flag value associated with the sequence
    /// * `primary` - The nucleotide sequence to write
    ///
    /// # Returns
    ///
    /// * `Ok(true)` - If the sequence was successfully written
    /// * `Ok(false)` - If the sequence was skipped due to invalid nucleotides
    /// * `Err(Error)` - If an error occurred during writing
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// * The sequence length doesn't match the header
    /// * Writing to the underlying writer fails
    pub fn write_nucleotides(&mut self, flag: u64, primary: &[u8]) -> Result<bool> {
        if let Some(sbuffer) = self.encoder.encode_single(primary)? {
            write_flag(&mut self.inner, flag)?;
            write_buffer(&mut self.inner, sbuffer)?;
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Writes a pair of nucleotide sequences to the output
    ///
    /// This method encodes and writes both a primary sequence and an extended sequence
    /// (e.g., quality scores) along with an associated flag. Both sequences must
    /// match their respective lengths specified in the header.
    ///
    /// # Arguments
    ///
    /// * `flag` - A 64-bit flag value associated with the sequences
    /// * `primary` - The primary nucleotide sequence
    /// * `extended` - The extended sequence (e.g., quality scores)
    ///
    /// # Returns
    ///
    /// * `Ok(true)` - If both sequences were successfully written
    /// * `Ok(false)` - If the sequences were skipped due to invalid nucleotides
    /// * `Err(Error)` - If an error occurred during writing
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// * Either sequence length doesn't match the header
    /// * Writing to the underlying writer fails
    pub fn write_paired(&mut self, flag: u64, primary: &[u8], extended: &[u8]) -> Result<bool> {
        if let Some((sbuffer, xbuffer)) = self.encoder.encode_paired(primary, extended)? {
            write_flag(&mut self.inner, flag)?;
            write_buffer(&mut self.inner, sbuffer)?;
            write_buffer(&mut self.inner, xbuffer)?;
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
    /// # use binseq::{BinseqHeader, BinseqWriterBuilder, Result};
    /// # fn main() -> Result<()> {
    /// let header = BinseqHeader::new(100);
    /// let writer = BinseqWriterBuilder::default()
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
    /// ownership of the `BinseqWriter`.
    pub fn by_ref(&mut self) -> &mut W {
        self.inner.by_ref()
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
    pub fn ingest(&mut self, other: &mut BinseqWriter<Vec<u8>>) -> Result<()> {
        let other_inner = other.by_ref();
        self.inner.write_all(other_inner.by_ref())?;
        other_inner.clear();
        Ok(())
    }
}

#[cfg(test)]
mod testing {

    use std::{fs::File, io::BufWriter};

    use super::*;
    use crate::SIZE_HEADER;

    #[test]
    fn test_headless() -> Result<()> {
        let inner = Vec::new();
        let mut writer = BinseqWriterBuilder::default()
            .header(BinseqHeader::new(32))
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
        let mut writer = BinseqWriterBuilder::default()
            .header(BinseqHeader::new(32))
            .build(inner)?;
        assert!(!writer.is_headless());
        let inner = writer.by_ref();
        assert_eq!(inner.len(), SIZE_HEADER);
        Ok(())
    }

    #[test]
    fn test_stdout() -> Result<()> {
        let writer = BinseqWriterBuilder::default()
            .header(BinseqHeader::new(32))
            .build(std::io::stdout())?;
        assert!(!writer.is_headless());
        Ok(())
    }

    #[test]
    fn test_to_path() -> Result<()> {
        let path = "test_to_path.file";
        let inner = File::create(path).map(BufWriter::new)?;
        let mut writer = BinseqWriterBuilder::default()
            .header(BinseqHeader::new(32))
            .build(inner)?;
        assert!(!writer.is_headless());
        let inner = writer.by_ref();
        inner.flush()?;

        // delete file
        std::fs::remove_file(path)?;

        Ok(())
    }
}
