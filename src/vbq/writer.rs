//! Writer implementation for VBQ files.
//!
//! Records are packed into fixed-size blocks, optionally zstd-compressed,
//! producing:
//!
//! ```text
//! [File Header][Data Blocks][Compressed Index][Index Size][Index End Magic]
//! ```
//!
//! Call `finish()` when done so the embedded index is written.

use std::io::Write;

use zstd::stream::copy_encode;

use super::header::{BlockHeader, FileHeader};
use crate::SequencingRecord;
use crate::encoder::Encoder;
use crate::error::{Result, WriteError};
use crate::policy::Policy;
use crate::vbq::header::{SIZE_BLOCK_HEADER, SIZE_HEADER};
use crate::vbq::index::{INDEX_END_MAGIC, IndexHeader};
use crate::vbq::{BlockIndex, BlockRange};

/// Builder for configured [`Writer`] instances.
#[derive(Default)]
pub struct WriterBuilder {
    /// Header of the file
    header: Option<FileHeader>,
    /// Optional policy for encoding
    policy: Option<Policy>,
    /// Optional headless mode (used in parallel writing)
    headless: Option<bool>,
}
impl WriterBuilder {
    /// Sets the file header (block size, quality, pairing, compression, etc.).
    #[must_use]
    pub fn header(mut self, header: FileHeader) -> Self {
        self.header = Some(header);
        self
    }

    /// Sets the encoding policy for invalid nucleotides.
    #[must_use]
    pub fn policy(mut self, policy: Policy) -> Self {
        self.policy = Some(policy);
        self
    }

    /// Sets headless mode: skips writing the file header, for parts that
    /// will be merged into another file later (e.g. parallel writing).
    #[must_use]
    pub fn headless(mut self, headless: bool) -> Self {
        self.headless = Some(headless);
        self
    }

    /// Builds a `Writer` around `inner`, using defaults for unset options.
    pub fn build<W: Write>(self, inner: W) -> Result<Writer<W>> {
        Writer::new(
            inner,
            self.header.unwrap_or_default(),
            self.policy.unwrap_or_default(),
            self.headless.unwrap_or(false),
        )
    }
}

/// Block-based writer for VBQ files.
///
/// Fills fixed-size blocks with encoded records, starting a new block when a
/// record would overflow the current one. Create via [`WriterBuilder`]; call
/// [`Writer::finish`] (also invoked on drop) to flush and embed the index.
#[derive(Clone)]
pub struct Writer<W: Write> {
    /// Inner Writer
    inner: W,

    /// Header of the file
    header: FileHeader,

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
impl<W: Write> Writer<W> {
    pub fn new(inner: W, header: FileHeader, policy: Policy, headless: bool) -> Result<Self> {
        let mut wtr = Self {
            inner,
            header,
            encoder: Encoder::with_policy(header.bits, policy),
            cblock: BlockWriter::new(
                header.block as usize,
                header.compressed,
                header.flags,
                header.qual,
                header.headers,
            ),
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

    /// Writes the file header; called on creation unless headless.
    fn init(&mut self) -> Result<()> {
        self.header.write_bytes(&mut self.inner)?;
        self.bytes_written += SIZE_HEADER;
        Ok(())
    }

    /// Returns whether the writer is configured for paired-end reads.
    pub fn is_paired(&self) -> bool {
        self.header.paired
    }

    /// Returns the header of the writer
    pub fn header(&self) -> FileHeader {
        self.header
    }

    /// Returns the N-policy of the writer
    pub fn policy(&self) -> Policy {
        self.encoder.policy()
    }

    /// Returns whether the writer is configured for quality scores.
    pub fn has_quality(&self) -> bool {
        self.header.qual
    }

    pub fn has_headers(&self) -> bool {
        self.header.headers
    }

    /// Flush the current block and record its range in the index.
    ///
    /// Free-standing over the individual fields so callers holding a borrow of
    /// `self.encoder` (the encoded record buffers) can still flush.
    fn flush_block_split(
        inner: &mut W,
        cblock: &mut BlockWriter,
        ranges: &mut Vec<BlockRange>,
        bytes_written: &mut usize,
        records_written: &mut usize,
    ) -> Result<()> {
        let block_header = cblock.flush(inner)?;
        ranges.push(BlockRange::new(
            *bytes_written as u64,
            block_header.size,
            block_header.records,
            *records_written as u64,
        ));
        *bytes_written += block_header.size_with_header();
        *records_written += block_header.records as usize;
        Ok(())
    }

    /// Flush the current block and record its range in the index
    fn flush_block(&mut self) -> Result<()> {
        Self::flush_block_split(
            &mut self.inner,
            &mut self.cblock,
            &mut self.ranges,
            &mut self.bytes_written,
            &mut self.records_written,
        )
    }

    /// Writes a [`SequencingRecord`], the unified API shared with the BQ and CBQ writers.
    ///
    /// Returns `Ok(false)` if the record was skipped by the encoding policy
    /// (e.g. invalid nucleotides).
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::{WriterBuilder, FileHeaderBuilder};
    /// use binseq::SequencingRecordBuilder;
    /// use std::fs::File;
    ///
    /// let header = FileHeaderBuilder::new().qual(true).headers(true).build();
    /// let mut writer = WriterBuilder::default()
    ///     .header(header)
    ///     .build(File::create("example.vbq").unwrap())
    ///     .unwrap();
    ///
    /// let record = SequencingRecordBuilder::default()
    ///     .s_seq(b"ACGTACGT")
    ///     .s_qual(b"IIIIFFFF")
    ///     .s_header(b"seq_001")
    ///     .flag(42)
    ///     .build()
    ///     .unwrap();
    ///
    /// writer.push(record).unwrap();
    /// writer.finish().unwrap();
    /// ```
    pub fn push(&mut self, record: SequencingRecord) -> Result<bool> {
        // Check paired status - writer can require paired (record must have R2),
        // but if writer is single-end, we simply ignore any R2 data in the record.
        if self.header.paired && !record.is_paired() {
            return Err(WriteError::ConfigurationMismatch {
                attribute: "paired",
                expected: self.header.paired,
                actual: record.is_paired(),
            }
            .into());
        }

        // For qualities and headers: the writer can require them (record must have them),
        // but if the writer doesn't need them, we simply ignore any extra data in the record.
        if self.header.qual && !record.has_qualities() {
            return Err(WriteError::ConfigurationMismatch {
                attribute: "qual",
                expected: self.header.qual,
                actual: record.has_qualities(),
            }
            .into());
        }
        if self.header.headers && !record.has_headers() {
            return Err(WriteError::ConfigurationMismatch {
                attribute: "headers",
                expected: self.header.headers,
                actual: record.has_headers(),
            }
            .into());
        }

        let record_size = record.configured_size_vbq(
            self.header.paired,
            self.header.flags,
            self.header.headers,
            self.header.qual,
            self.header.bits,
        );

        // encode the sequence(s); a `None` means the record was skipped by policy
        let encoded = if self.header.is_paired() {
            self.encoder
                .encode_paired(record.s_seq, record.x_seq.unwrap_or_default())?
                .map(|(sbuffer, xbuffer)| (sbuffer, Some(xbuffer)))
        } else {
            self.encoder
                .encode_single(record.s_seq)?
                .map(|sbuffer| (sbuffer, None))
        };
        let Some((sbuffer, xbuffer)) = encoded else {
            return Ok(false);
        };

        // flush the current block if this record would overflow it; done only
        // once we know the record will actually be written, so policy-skipped
        // records never fragment blocks (split-borrow: `encoded` holds the encoder)
        if self.cblock.exceeds_block_size(record_size)? {
            Self::flush_block_split(
                &mut self.inner,
                &mut self.cblock,
                &mut self.ranges,
                &mut self.bytes_written,
                &mut self.records_written,
            )?;
        }

        self.cblock.write_record(&record, sbuffer, xbuffer)?;
        Ok(true)
    }

    /// Flushes remaining data and writes the embedded index.
    ///
    /// Called automatically on drop, but calling it explicitly lets you handle
    /// errors. Idempotent: the index is only written once.
    pub fn finish(&mut self) -> Result<()> {
        // Flush any remaining data in the current block
        self.flush_block()?;
        self.inner.flush()?;

        // Always write the index - this is critical for VBQ file validity
        // The index_written flag prevents double-writing on subsequent finish() calls
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

    /// Transfers all blocks (complete and partial) from an in-memory writer
    /// into this one, clearing the other for reuse. Used to combine per-thread
    /// buffers in parallel writing.
    ///
    /// Errors if the two writers' headers are incompatible.
    pub fn ingest(&mut self, other: &mut Writer<Vec<u8>>) -> Result<()> {
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

    fn write_index(&mut self) -> Result<()> {
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
        self.inner.write_all(&n_bytes.to_le_bytes())?;

        // Write the index footer magic
        self.inner.write_all(&INDEX_END_MAGIC.to_le_bytes())?;

        Ok(())
    }
}

impl<W: Write> Drop for Writer<W> {
    fn drop(&mut self) {
        self.finish().expect("Writer: Failed to finish writing");
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
    /// Has quality scores
    has_qualities: bool,
    /// Has headers
    has_headers: bool,
}
impl BlockWriter {
    fn new(
        block_size: usize,
        compress: bool,
        has_flags: bool,
        has_qualities: bool,
        has_headers: bool,
    ) -> Self {
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
            has_qualities,
            has_headers,
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

    fn write_record(
        &mut self,
        record: &SequencingRecord,
        sbuf: &[u64],
        xbuf: Option<&[u64]>,
    ) -> Result<()> {
        // Tracks the record start position
        self.starts.push(self.pos);

        // Write the flag (only if configured)
        if self.has_flags {
            self.write_flag(record.flag.unwrap_or(0))?;
        }

        // Write the lengths
        self.write_length(record.s_seq.len() as u64)?;
        self.write_length(record.x_seq.map_or(0, <[u8]>::len) as u64)?;

        // Write the primary sequence
        self.write_buffer(sbuf)?;

        // Write primary quality (only if configured)
        if self.has_qualities
            && let Some(qual) = record.s_qual
        {
            self.write_u8buf(qual)?;
        }

        // Write primary header (only if configured)
        if self.has_headers
            && let Some(sheader) = record.s_header
        {
            self.write_length(sheader.len() as u64)?;
            self.write_u8buf(sheader)?;
        }

        // Write the optional extended sequence
        if let Some(xbuf) = xbuf {
            self.write_buffer(xbuf)?;
        }

        // Write extended quality (only if configured)
        if self.has_qualities
            && let Some(qual) = record.x_qual
        {
            self.write_u8buf(qual)?;
        }

        // Write extended header (only if configured)
        if self.has_headers
            && let Some(xheader) = record.x_header
        {
            self.write_length(xheader.len() as u64)?;
            self.write_u8buf(xheader)?;
        }

        Ok(())
    }

    fn write_flag(&mut self, flag: u64) -> Result<()> {
        self.ubuf.write_all(&flag.to_le_bytes())?;
        self.pos += 8;
        Ok(())
    }

    fn write_length(&mut self, length: u64) -> Result<()> {
        self.ubuf.write_all(&length.to_le_bytes())?;
        self.pos += 8;
        Ok(())
    }

    fn write_buffer(&mut self, ebuf: &[u64]) -> Result<()> {
        ebuf.iter()
            .try_for_each(|&x| self.ubuf.write_all(&x.to_le_bytes()))?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SequencingRecordBuilder;
    use crate::utils::read_u64_le;
    use crate::vbq::{FileHeaderBuilder, header::SIZE_HEADER};

    #[test]
    fn test_policy_skipped_records_do_not_fragment_blocks() -> super::Result<()> {
        // Two writers fed the same valid records must produce identical output,
        // even when one is also fed invalid (policy-skipped) records that would
        // overflow the current block: a skipped record must never trigger a flush.
        let header = FileHeaderBuilder::new().block(2048).build();
        let mut plain = WriterBuilder::default().header(header).build(Vec::new())?;
        let mut skipped = WriterBuilder::default().header(header).build(Vec::new())?;

        let valid = b"ACGT".repeat(25); // 100 bp
        let invalid = b"N".repeat(300); // always skipped by IgnoreSequence
        for _ in 0..200 {
            let invalid_record = SequencingRecordBuilder::default().s_seq(&invalid).build()?;
            assert!(!skipped.push(invalid_record)?);

            let valid_record = SequencingRecordBuilder::default().s_seq(&valid).build()?;
            assert!(plain.push(valid_record)?);
            let valid_record = SequencingRecordBuilder::default().s_seq(&valid).build()?;
            assert!(skipped.push(valid_record)?);
        }
        plain.finish()?;
        skipped.finish()?;

        assert!(plain.ranges.len() > 1, "test should span multiple blocks");
        assert_eq!(plain.ranges.len(), skipped.ranges.len());
        assert_eq!(plain.inner, skipped.inner);
        Ok(())
    }

    #[test]
    fn test_headless_writer() -> super::Result<()> {
        let writer = WriterBuilder::default().headless(true).build(Vec::new())?;
        assert_eq!(writer.inner.len(), 0);

        let writer = WriterBuilder::default().headless(false).build(Vec::new())?;
        assert_eq!(writer.inner.len(), SIZE_HEADER);

        Ok(())
    }

    #[test]
    fn test_ingest_empty_writer() -> super::Result<()> {
        // Test ingesting from an empty writer
        let header = FileHeaderBuilder::new().build();

        // Create a source writer that's empty
        let mut source = WriterBuilder::default()
            .header(header)
            .headless(true)
            .build(Vec::new())?;

        // Create a destination writer
        let mut dest = WriterBuilder::default()
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
        let header = FileHeaderBuilder::new().build();

        // Create a source writer with a single record
        let mut source = WriterBuilder::default()
            .header(header)
            .headless(true)
            .build(Vec::new())?;

        // Write a single sequence
        let record = SequencingRecordBuilder::default()
            .s_seq(b"ACGTACGTACGT")
            .flag(1)
            .build()?;
        source.push(record)?;

        // We have not crossed a boundary
        assert!(source.by_ref().is_empty());

        // Create a destination writer
        let mut dest = WriterBuilder::default()
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
        let header = FileHeaderBuilder::new().build();

        // Create a source writer with a single record
        let mut source = WriterBuilder::default()
            .header(header)
            .headless(true)
            .build(Vec::new())?;

        // Write multiple sequences
        for _ in 0..30 {
            let record = SequencingRecordBuilder::default()
                .s_seq(b"ACGTACGTACGT")
                .flag(1)
                .build()?;
            source.push(record)?;
        }
        // We have not crossed a boundary
        assert!(source.by_ref().is_empty());

        // Create a destination writer
        let mut dest = WriterBuilder::default()
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
        let header = FileHeaderBuilder::new().build();

        // Create a source writer with a single record
        let mut source = WriterBuilder::default()
            .header(header)
            .headless(true)
            .build(Vec::new())?;

        // Write multiple sequences (will cross boundary)
        for _ in 0..30000 {
            let record = SequencingRecordBuilder::default()
                .s_seq(b"ACGTACGTACGT")
                .flag(1)
                .build()?;
            source.push(record)?;
        }

        // We have crossed a boundary
        assert!(!source.by_ref().is_empty());

        // Create a destination writer
        let mut dest = WriterBuilder::default()
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
        let source_header = FileHeaderBuilder::new().qual(true).build();
        let dest_header = FileHeaderBuilder::new().qual(true).build();

        // Create a source writer with quality scores
        let mut source = WriterBuilder::default()
            .header(source_header)
            .headless(true)
            .build(Vec::new())?;

        // Write sequences with quality scores
        let seq = b"ACGTACGTACGT";
        let qual = vec![40u8; seq.len()];
        for i in 0..5 {
            let record = SequencingRecordBuilder::default()
                .s_seq(seq)
                .s_qual(&qual)
                .flag(i)
                .build()?;
            source.push(record)?;
        }

        // Create a destination writer
        let mut dest = WriterBuilder::default()
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
        let header = FileHeaderBuilder::new().compressed(true).build();

        // Create a source writer with a single record
        let mut source = WriterBuilder::default()
            .header(header)
            .headless(true)
            .build(Vec::new())?;

        // Write multiple sequences (will cross boundary)
        for _ in 0..30000 {
            let record = SequencingRecordBuilder::default()
                .s_seq(b"ACGTACGTACGT")
                .flag(1)
                .build()?;
            source.push(record)?;
        }

        // Create a destination writer
        let mut dest = WriterBuilder::default()
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
        let source_header = FileHeaderBuilder::new().build();
        let dest_header = FileHeaderBuilder::new().qual(true).build();

        // Create a source writer with quality scores
        let mut source = WriterBuilder::default()
            .header(source_header)
            .headless(true)
            .build(Vec::new())?;

        // Create a destination writer
        let mut dest = WriterBuilder::default()
            .header(dest_header)
            .headless(true)
            .build(Vec::new())?;

        // Ingest from source to dest (will error)
        assert!(dest.ingest(&mut source).is_err());

        Ok(())
    }

    #[test]
    fn test_index_always_written_on_finish() -> super::Result<()> {
        use crate::vbq::index::INDEX_END_MAGIC;

        // Create a writer with some records
        let header = FileHeaderBuilder::new().build();
        let mut writer = WriterBuilder::default().header(header).build(Vec::new())?;

        // Write some records
        for i in 0..10 {
            let record = SequencingRecordBuilder::default()
                .s_seq(b"ACGTACGTACGT")
                .flag(i)
                .build()?;
            writer.push(record)?;
        }

        // Finish the writer
        writer.finish()?;

        // Get the written bytes
        let bytes = &writer.inner;

        // Verify the file ends with the index magic number
        assert!(bytes.len() >= 8, "File is too short to contain index");
        let magic_offset = bytes.len() - 8;
        let magic = read_u64_le(&bytes[magic_offset..]);
        assert_eq!(
            magic, INDEX_END_MAGIC,
            "Index magic number not found at end of file"
        );

        // Verify we can read the index size
        assert!(bytes.len() >= 16, "File is too short to contain index size");
        let size_offset = bytes.len() - 16;
        let index_size = read_u64_le(&bytes[size_offset..size_offset + 8]);
        assert!(index_size > 0, "Index size should be greater than 0");

        // Verify the index size makes sense (should be less than total file size)
        assert!(
            index_size < bytes.len() as u64,
            "Index size is larger than file"
        );

        Ok(())
    }

    #[test]
    fn test_finish_idempotent() -> super::Result<()> {
        use crate::vbq::index::INDEX_END_MAGIC;

        // Create a writer
        let header = FileHeaderBuilder::new().build();
        let mut writer = WriterBuilder::default().header(header).build(Vec::new())?;

        // Write some records
        for i in 0..10 {
            let record = SequencingRecordBuilder::default()
                .s_seq(b"ACGTACGTACGT")
                .flag(i)
                .build()?;
            writer.push(record)?;
        }

        // Call finish() multiple times
        writer.finish()?;
        let size_after_first_finish = writer.inner.len();

        writer.finish()?;
        let size_after_second_finish = writer.inner.len();

        writer.finish()?;
        let size_after_third_finish = writer.inner.len();

        // All sizes should be the same - index should only be written once
        assert_eq!(size_after_first_finish, size_after_second_finish);
        assert_eq!(size_after_second_finish, size_after_third_finish);

        // Verify only one index magic number at the end
        let bytes = &writer.inner;
        let magic_offset = bytes.len() - 8;
        let magic = read_u64_le(&bytes[magic_offset..]);
        assert_eq!(magic, INDEX_END_MAGIC);

        Ok(())
    }
}
