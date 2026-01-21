use std::io;

use bitnuc::BitSize;
use bytemuck::{cast_slice, cast_slice_mut};
use sucds::Serializable;
use sucds::mii_sequences::{EliasFano, EliasFanoBuilder};
use zstd::stream::copy_decode;
use zstd::zstd_safe;

use crate::cbq::core::utils::sized_compress;
use crate::error::{CbqError, WriteError};
use crate::{BinseqRecord, DEFAULT_QUALITY_SCORE, Result};

use super::utils::{Span, calculate_offsets, extension_read, resize_uninit, slice_and_increment};
use super::{BlockHeader, BlockRange, FileHeader};
use crate::SequencingRecord;

/// A block of records where all data is stored in separate columns.
#[derive(Clone, Default)]
pub struct ColumnarBlock {
    /// Separate columns for each data type
    seq: Vec<u8>,
    flags: Vec<u64>,
    headers: Vec<u8>,
    qual: Vec<u8>,

    /// Length of sequences for each record
    pub(crate) l_seq: Vec<u64>,
    /// Length of headers for each record
    pub(crate) l_headers: Vec<u64>,
    /// Position of all N's in the sequence
    pub(crate) npos: Vec<u64>,

    /// Reusable buffer for encoding sequences
    ebuf: Vec<u64>,

    /// An Elias-Fano encoding for the N-positions
    pub(crate) ef: Option<EliasFano>,
    /// Reusable buffer for encoding Elias-Fano struct
    pub(crate) ef_bytes: Vec<u8>,
    /// Length of serialized Elias-Fano encoding in bytes
    pub(crate) len_nef: usize,

    // Reusable zstd compression buffer for columnar data
    pub(crate) z_seq_len: Vec<u8>,
    pub(crate) z_header_len: Vec<u8>,
    pub(crate) z_npos: Vec<u8>,
    pub(crate) z_seq: Vec<u8>,
    pub(crate) z_flags: Vec<u8>,
    pub(crate) z_headers: Vec<u8>,
    pub(crate) z_qual: Vec<u8>,

    // reusable offset buffers
    l_seq_offsets: Vec<u64>,
    l_header_offsets: Vec<u64>,

    /// Number of records in the block
    ///
    /// A record is a logical unit of data.
    /// If the records are paired sequences this is the number of pairs.
    pub(crate) num_records: usize,

    /// Number of sequences in the block
    ///
    /// This is the same as the number of records for unpaired sequences.
    /// For paired sequences it will be twice the number of records.
    pub(crate) num_sequences: usize,

    /// Total nucleotides in this block
    pub(crate) nuclen: usize,
    /// Number of npos positions
    pub(crate) num_npos: usize,
    /// Current size of this block (virtual)
    current_size: usize,

    /// Reusable buffer for missing quality scores
    qbuf: Vec<u8>,
    default_quality_score: u8,

    /// The file header (used for block configuration)
    ///
    /// Not to be confused with the `BlockHeader`
    pub(crate) header: FileHeader,
}
impl ColumnarBlock {
    /// Create a new columnar block with the given block size
    #[must_use]
    pub fn new(header: FileHeader) -> Self {
        Self {
            header,
            default_quality_score: DEFAULT_QUALITY_SCORE,
            ..Default::default()
        }
    }

    /// Update the default quality score for this block
    pub fn set_default_quality_score(&mut self, score: u8) {
        self.default_quality_score = score;
        self.qbuf.clear();
    }

    fn is_empty(&self) -> bool {
        self.current_size == 0
    }

    /// Clears the internal data structures
    pub(crate) fn clear(&mut self) {
        // clear index counters
        {
            self.nuclen = 0;
            self.num_sequences = 0;
            self.num_records = 0;
            self.current_size = 0;
            self.num_npos = 0;
        }

        // clear spans
        {
            self.l_seq.clear();
            self.l_headers.clear();
            self.l_seq_offsets.clear();
            self.l_header_offsets.clear();
        }

        // clear vectors
        {
            self.seq.clear();
            self.flags.clear();
            self.headers.clear();
            self.qual.clear();
            self.npos.clear();
            self.ef = None;
        }

        // clear encodings
        {
            self.ebuf.clear();
            self.z_seq_len.clear();
            self.z_header_len.clear();
            self.z_npos.clear();
            self.z_seq.clear();
            self.z_flags.clear();
            self.z_headers.clear();
            self.z_qual.clear();
            self.ef_bytes.clear();
        }
    }

    fn add_sequence(&mut self, record: &SequencingRecord) {
        self.l_seq.push(record.s_seq.len() as u64);
        self.seq.extend_from_slice(record.s_seq);
        self.num_sequences += 1;
        if let Some(x_seq) = record.x_seq {
            self.l_seq.push(x_seq.len() as u64);
            self.seq.extend_from_slice(x_seq);
            self.num_sequences += 1;
        }

        // keep the sequence size up to date
        self.nuclen = self.seq.len();
    }

    fn add_flag(&mut self, record: &SequencingRecord) {
        if let Some(flag) = record.flag {
            self.flags.push(flag);
        }
    }

    fn add_headers(&mut self, record: &SequencingRecord) {
        if let Some(header) = record.s_header {
            self.l_headers.push(header.len() as u64);
            self.headers.extend_from_slice(header);
        }
        if let Some(header) = record.x_header {
            self.l_headers.push(header.len() as u64);
            self.headers.extend_from_slice(header);
        }
    }

    /// Note: this does not check if quality scores are different lengths from sequence
    fn add_quality(&mut self, record: &SequencingRecord) {
        if let Some(qual) = record.s_qual {
            self.qual.extend_from_slice(qual);
        }
        if let Some(qual) = record.x_qual {
            self.qual.extend_from_slice(qual);
        }
    }

    /// Calculate the usage of the block as a percentage
    #[must_use]
    pub fn usage(&self) -> f64 {
        self.current_size as f64 / self.header.block_size as f64
    }

    pub(crate) fn can_fit(&self, record: &SequencingRecord<'_>) -> bool {
        self.current_size + record.size() <= self.header.block_size as usize
    }

    pub(crate) fn can_ingest(&self, other: &Self) -> bool {
        self.current_size + other.current_size <= self.header.block_size as usize
    }

    /// Ensure that the record can be pushed into the block
    fn validate_record(&self, record: &SequencingRecord) -> Result<()> {
        if !self.can_fit(record) {
            if record.size() > self.header.block_size as usize {
                return Err(WriteError::RecordSizeExceedsMaximumBlockSize(
                    record.size(),
                    self.header.block_size as usize,
                )
                .into());
            }
            return Err(CbqError::BlockFull {
                current_size: self.current_size,
                record_size: record.size(),
                block_size: self.header.block_size as usize,
            }
            .into());
        }

        if record.is_paired() != self.header.is_paired() {
            return Err(WriteError::ConfigurationMismatch {
                attribute: "paired",
                expected: self.header.is_paired(),
                actual: record.is_paired(),
            }
            .into());
        }

        if record.has_flags() != self.header.has_flags() {
            return Err(WriteError::ConfigurationMismatch {
                attribute: "flags",
                expected: self.header.has_flags(),
                actual: record.has_flags(),
            }
            .into());
        }

        if record.has_headers() != self.header.has_headers() {
            return Err(WriteError::ConfigurationMismatch {
                attribute: "headers",
                expected: self.header.has_headers(),
                actual: record.has_headers(),
            }
            .into());
        }

        if record.has_qualities() != self.header.has_qualities() {
            return Err(WriteError::ConfigurationMismatch {
                attribute: "qualities",
                expected: self.header.has_qualities(),
                actual: record.has_qualities(),
            }
            .into());
        }
        Ok(())
    }

    pub fn push(&mut self, record: SequencingRecord) -> Result<()> {
        self.validate_record(&record)?;

        self.add_sequence(&record);
        self.add_flag(&record);
        self.add_headers(&record);
        self.add_quality(&record);
        self.current_size += record.size();
        self.num_records += 1;

        Ok(())
    }

    /// Returns the expected length of the encoded sequence buffer
    ///
    /// This is deterministically calculated based on the sequence length and the encoding scheme.
    fn ebuf_len(&self) -> usize {
        self.nuclen.div_ceil(32)
    }

    /// Encode the sequence into a compressed representation
    fn encode_sequence(&mut self) -> Result<()> {
        bitnuc::twobit::encode_with_invalid(&self.seq, &mut self.ebuf)?;
        Ok(())
    }

    /// Find all positions of 'N' in the sequence
    fn fill_npos(&mut self) -> Result<()> {
        self.npos
            .extend(memchr::memchr_iter(b'N', &self.seq).map(|i| i as u64));
        self.num_npos = self.npos.len();

        // build Elias-Fano encoding for N positions
        if self.npos.is_empty() {
            self.ef = None;
            Ok(())
        } else {
            let mut ef_builder = EliasFanoBuilder::new(self.seq.len(), self.npos.len())?;
            ef_builder.extend(self.npos.iter().map(|idx| *idx as usize))?;
            let ef = ef_builder.build();

            self.ef = Some(ef);
            Ok(())
        }
    }

    /// Convert all ambiguous bases back to N
    fn backfill_npos(&mut self) {
        if let Some(ef) = self.ef.as_ref() {
            ef.iter(0).for_each(|idx| {
                if let Some(base) = self.seq.get_mut(idx) {
                    *base = b'N';
                }
            });
        }
    }

    /// Compress all native columns into compressed representation
    fn compress_columns(&mut self, cctx: &mut zstd_safe::CCtx) -> Result<()> {
        // compress sequence lengths

        sized_compress(&mut self.z_seq_len, cast_slice(&self.l_seq), cctx)?;

        if !self.headers.is_empty() {
            sized_compress(&mut self.z_header_len, cast_slice(&self.l_headers), cctx)?;
        }

        // compress N-positions (Elias-Fano encoded)
        if let Some(ef) = self.ef.as_ref() {
            ef.serialize_into(&mut self.ef_bytes)?;
            self.len_nef = self.ef_bytes.len();
            sized_compress(&mut self.z_npos, &self.ef_bytes, cctx)?;
        }

        // compress sequence
        sized_compress(&mut self.z_seq, cast_slice(&self.ebuf), cctx)?;

        // compress flags
        if !self.flags.is_empty() {
            sized_compress(&mut self.z_flags, cast_slice(&self.flags), cctx)?;
        }

        // compress headers
        if !self.headers.is_empty() {
            sized_compress(&mut self.z_headers, cast_slice(&self.headers), cctx)?;
        }

        // compress quality
        if !self.qual.is_empty() {
            sized_compress(&mut self.z_qual, cast_slice(&self.qual), cctx)?;
        }

        Ok(())
    }

    /// Decompress all columns back to native representation
    pub fn decompress_columns(&mut self) -> Result<()> {
        // decompress sequence lengths
        {
            self.l_seq.resize(self.num_sequences, 0);
            copy_decode(self.z_seq_len.as_slice(), cast_slice_mut(&mut self.l_seq))?;
        }

        // decompress header lengths
        if !self.z_header_len.is_empty() {
            self.l_headers.resize(self.num_sequences, 0);
            copy_decode(
                self.z_header_len.as_slice(),
                cast_slice_mut(&mut self.l_headers),
            )?;
        }

        // decompress npos
        if !self.z_npos.is_empty() {
            self.ef_bytes.resize(self.len_nef, 0);
            copy_decode(self.z_npos.as_slice(), &mut self.ef_bytes)?;

            let ef = EliasFano::deserialize_from(self.ef_bytes.as_slice())?;
            self.num_npos = ef.len();
            self.ef = Some(ef);
        }

        // decompress sequence
        {
            self.ebuf.resize(self.ebuf_len(), 0);
            copy_decode(self.z_seq.as_slice(), cast_slice_mut(&mut self.ebuf))?;

            bitnuc::twobit::decode(&self.ebuf, self.nuclen, &mut self.seq)?;
            self.backfill_npos();
        }

        // decompress flags
        if !self.z_flags.is_empty() {
            self.flags.resize(self.num_records, 0);
            copy_decode(self.z_flags.as_slice(), cast_slice_mut(&mut self.flags))?;
        }

        // decompress headers
        if !self.z_headers.is_empty() {
            copy_decode(self.z_headers.as_slice(), &mut self.headers)?;
        }

        // decompress quality scores
        if !self.z_qual.is_empty() {
            copy_decode(self.z_qual.as_slice(), &mut self.qual)?;
        }

        // calculate offsets
        {
            calculate_offsets(&self.l_seq, &mut self.l_seq_offsets);
            calculate_offsets(&self.l_headers, &mut self.l_header_offsets);
        }

        Ok(())
    }

    fn write<W: io::Write>(&mut self, writer: &mut W) -> Result<()> {
        writer.write_all(&self.z_seq_len)?;
        writer.write_all(&self.z_header_len)?;
        writer.write_all(&self.z_npos)?;
        writer.write_all(&self.z_seq)?;
        writer.write_all(&self.z_flags)?;
        writer.write_all(&self.z_headers)?;
        writer.write_all(&self.z_qual)?;
        Ok(())
    }

    pub fn flush_to<W: io::Write>(
        &mut self,
        writer: &mut W,
        cctx: &mut zstd_safe::CCtx,
    ) -> Result<Option<BlockHeader>> {
        if self.is_empty() {
            return Ok(None);
        }

        // encode all sequences at once
        self.encode_sequence()?;

        // fill npos
        self.fill_npos()?;

        // compress each column
        self.compress_columns(cctx)?;

        // build the block header
        let header = BlockHeader::from_block(self);
        // eprintln!("{header:?}");

        // write the block header
        header.write(writer)?;

        // write the internal state to the inner writer
        self.write(writer)?;

        // clear the internal state
        self.clear();

        Ok(Some(header))
    }

    pub fn read_from<R: io::Read>(&mut self, reader: &mut R, header: BlockHeader) -> Result<()> {
        // clears the internal state
        self.clear();

        // reload the internal state from the reader
        self.nuclen = header.nuclen as usize;
        self.num_records = header.num_records as usize;
        self.num_sequences = header.num_sequences as usize;
        self.len_nef = header.len_nef as usize;

        extension_read(reader, &mut self.z_seq_len, header.len_z_seq_len as usize)?;
        extension_read(
            reader,
            &mut self.z_header_len,
            header.len_z_header_len as usize,
        )?;
        extension_read(reader, &mut self.z_npos, header.len_z_npos as usize)?;
        extension_read(reader, &mut self.z_seq, header.len_z_seq as usize)?;
        extension_read(reader, &mut self.z_flags, header.len_z_flags as usize)?;
        extension_read(reader, &mut self.z_headers, header.len_z_headers as usize)?;
        extension_read(reader, &mut self.z_qual, header.len_z_qual as usize)?;
        Ok(())
    }

    pub fn decompress_from_bytes(
        &mut self,
        bytes: &[u8],
        header: BlockHeader,
        dctx: &mut zstd_safe::DCtx,
    ) -> Result<()> {
        // clears the internal state
        self.clear();

        // reload the internal state from the header
        self.nuclen = header.nuclen as usize;
        self.num_records = header.num_records as usize;
        self.num_sequences = header.num_sequences as usize;
        self.len_nef = header.len_nef as usize;

        let mut byte_offset = 0;

        // decompress sequence lengths
        {
            resize_uninit(&mut self.l_seq, self.num_sequences);
            dctx.decompress(
                cast_slice_mut(&mut self.l_seq),
                slice_and_increment(&mut byte_offset, header.len_z_seq_len, bytes),
            )
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;

            // update default quality score buffer size
            self.l_seq.iter().for_each(|len| {
                if *len as usize > self.qbuf.len() {
                    self.qbuf.resize(*len as usize, self.default_quality_score);
                }
            })
        }

        // decompress header lengths
        if header.len_z_header_len > 0 {
            resize_uninit(&mut self.l_headers, self.num_sequences);
            dctx.decompress(
                cast_slice_mut(&mut self.l_headers),
                slice_and_increment(&mut byte_offset, header.len_z_header_len, bytes),
            )
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;
        }

        // calculate offsets
        {
            calculate_offsets(&self.l_seq, &mut self.l_seq_offsets);
            calculate_offsets(&self.l_headers, &mut self.l_header_offsets);
        }

        // decompress npos
        if header.len_z_npos > 0 {
            resize_uninit(&mut self.ef_bytes, self.len_nef);
            dctx.decompress(
                &mut self.ef_bytes,
                slice_and_increment(&mut byte_offset, header.len_z_npos, bytes),
            )
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;

            // reinitialize the EliasFano encoding
            let ef = EliasFano::deserialize_from(self.ef_bytes.as_slice())?;
            self.num_npos = ef.len();
            self.ef = Some(ef);
        }

        // decompress sequence
        {
            let ebuf_len = self.ebuf_len();
            resize_uninit(&mut self.ebuf, ebuf_len);
            dctx.decompress(
                cast_slice_mut(&mut self.ebuf),
                slice_and_increment(&mut byte_offset, header.len_z_seq, bytes),
            )
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;

            bitnuc::twobit::decode(&self.ebuf, self.nuclen, &mut self.seq)?;
            self.backfill_npos();
        }

        // decompress flags
        if header.len_z_flags > 0 {
            resize_uninit(&mut self.flags, self.num_records);
            dctx.decompress(
                cast_slice_mut(&mut self.flags),
                slice_and_increment(&mut byte_offset, header.len_z_flags, bytes),
            )
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;
        }

        // decompress headers
        if header.len_z_headers > 0 {
            let headers_len = (self.l_header_offsets.last().copied().unwrap_or(0)
                + self.l_headers.last().copied().unwrap_or(0))
                as usize;
            resize_uninit(&mut self.headers, headers_len);
            dctx.decompress(
                &mut self.headers,
                slice_and_increment(&mut byte_offset, header.len_z_headers, bytes),
            )
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;
        }

        // decompress quality scores
        if header.len_z_qual > 0 {
            resize_uninit(&mut self.qual, self.nuclen);
            dctx.decompress(
                &mut self.qual,
                slice_and_increment(&mut byte_offset, header.len_z_qual, bytes),
            )
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;
        }

        Ok(())
    }

    pub(crate) fn take_incomplete(&mut self, other: &Self) -> Result<()> {
        if !self.can_ingest(other) {
            return Err(CbqError::CannotIngestBlock {
                self_block_size: self.header.block_size as usize,
                other_block_size: other.header.block_size as usize,
            }
            .into());
        }

        // increment attributes
        {
            self.nuclen += other.nuclen;
            self.num_records += other.num_records;
            self.num_sequences += other.num_sequences;
            self.current_size += other.current_size;
        }

        // extend data
        {
            self.seq.extend_from_slice(&other.seq);
            self.flags.extend_from_slice(&other.flags);
            self.headers.extend_from_slice(&other.headers);
            self.qual.extend_from_slice(&other.qual);
            self.l_seq.extend_from_slice(&other.l_seq);
            self.l_headers.extend_from_slice(&other.l_headers);
        }

        {
            // Note:
            //
            // Remaining buffers and attributes are left untouched.
            // These are not modified because they aren't used mid-writing
            // and are populated during the flush step.
        }

        Ok(())
    }

    #[must_use]
    pub fn iter_records(&self, range: BlockRange) -> RefRecordIter<'_> {
        RefRecordIter {
            block: self,
            range,
            qbuf: &self.qbuf,
            index: 0,
            is_paired: self.header.is_paired(),
            has_headers: self.header.has_headers(),
            header_buffer: itoa::Buffer::new(),
        }
    }
}

/// A zero-copy iterator over [`RefRecord`](crate::cbq::RefRecord)s in a [`ColumnarBlock`](crate::cbq::ColumnarBlock)
pub struct RefRecordIter<'a> {
    /// The block containing the records
    block: &'a ColumnarBlock,

    /// The record range of this block
    range: BlockRange,

    /// Record index within the block
    index: usize,

    /// Convenience attribute if block is paired
    is_paired: bool,

    /// Convenience attribute if block has headers
    has_headers: bool,

    /// Preallocated buffer for quality scores
    qbuf: &'a [u8],

    /// Preallocated itoa buffer for converting global record index to string
    header_buffer: itoa::Buffer,
}
impl<'a> Iterator for RefRecordIter<'a> {
    type Item = RefRecord<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.block.num_records {
            None
        } else {
            // Calculate the actual array index
            let seq_idx = if self.is_paired {
                self.index * 2
            } else {
                self.index
            };

            let sseq_span =
                Span::new_u64(self.block.l_seq_offsets[seq_idx], self.block.l_seq[seq_idx]);
            let sheader_span = if self.has_headers {
                Some(Span::new_u64(
                    self.block.l_header_offsets[seq_idx],
                    self.block.l_headers[seq_idx],
                ))
            } else {
                None
            };
            let xseq_span = if self.is_paired {
                Some(Span::new_u64(
                    self.block.l_seq_offsets[seq_idx + 1],
                    self.block.l_seq[seq_idx + 1],
                ))
            } else {
                None
            };
            let xheader_span = if self.is_paired && self.has_headers {
                Some(Span::new_u64(
                    self.block.l_header_offsets[seq_idx + 1],
                    self.block.l_headers[seq_idx + 1],
                ))
            } else {
                None
            };

            let global_index =
                self.range.cumulative_records as usize - self.block.num_records + self.index;

            let rr_index = RefRecordIndex::new(global_index, &mut self.header_buffer);

            let record = RefRecord {
                block: self.block,
                index: self.index,
                qbuf: &self.qbuf,
                global_index,
                sseq_span,
                sheader_span,
                xseq_span,
                xheader_span,
                rr_index,
            };

            self.index += 1;
            Some(record)
        }
    }
}

/// A convenience struct for creating global indices as `&[u8]` buffers
#[derive(Clone, Copy)]
struct RefRecordIndex {
    index_buf: [u8; 20],
    index_len: usize,
}
impl RefRecordIndex {
    fn new(index: usize, itoa_buf: &mut itoa::Buffer) -> Self {
        let mut index_buf = [0u8; 20];
        let header_str = itoa_buf.format(index);
        let index_len = header_str.len();
        index_buf[..index_len].copy_from_slice(header_str.as_bytes());
        Self {
            index_buf,
            index_len,
        }
    }

    fn as_bytes(&self) -> &[u8] {
        &self.index_buf[..self.index_len]
    }
}

/// A reference to a record in a [`ColumnarBlock`](crate::cbq::ColumnarBlock) that implements the [`BinseqRecord`](crate::BinseqRecord) trait
#[derive(Clone, Copy)]
pub struct RefRecord<'a> {
    /// A reference to the block containing this record
    block: &'a ColumnarBlock,

    /// Preallocated buffer for quality scores
    qbuf: &'a [u8],

    /// Local index of this record within the block
    index: usize,

    /// Global index of this record in the file
    global_index: usize,

    /// Span of the primary sequence within the block
    sseq_span: Span,

    /// Span of the extended sequence within the block
    xseq_span: Option<Span>,

    /// Span of the primary header within the block
    sheader_span: Option<Span>,

    /// Span of the extended header within the block
    xheader_span: Option<Span>,

    /// A buffer to the name of this record when not storing headers
    rr_index: RefRecordIndex,
}
impl BinseqRecord for RefRecord<'_> {
    fn bitsize(&self) -> BitSize {
        BitSize::Two
    }

    fn index(&self) -> u64 {
        self.global_index as u64
    }

    fn flag(&self) -> Option<u64> {
        self.block.flags.get(self.index).copied()
    }

    fn is_paired(&self) -> bool {
        self.xseq_span.is_some()
    }

    fn sheader(&self) -> &[u8] {
        if let Some(span) = self.sheader_span {
            &self.block.headers[span.range()]
        } else {
            self.rr_index.as_bytes()
        }
    }

    fn xheader(&self) -> &[u8] {
        if let Some(span) = self.xheader_span {
            &self.block.headers[span.range()]
        } else {
            self.rr_index.as_bytes()
        }
    }

    fn sbuf(&self) -> &[u64] {
        unimplemented!("sbuf is not implemented for cbq")
    }

    fn xbuf(&self) -> &[u64] {
        unimplemented!("xbuf is not implemented for cbq")
    }

    fn slen(&self) -> u64 {
        self.sseq_span.len() as u64
    }

    fn xlen(&self) -> u64 {
        self.xseq_span.map_or(0, |span| span.len() as u64)
    }

    fn decode_s(&self, buf: &mut Vec<u8>) -> crate::Result<()> {
        buf.extend_from_slice(self.sseq());
        Ok(())
    }

    fn decode_x(&self, buf: &mut Vec<u8>) -> crate::Result<()> {
        buf.extend_from_slice(self.xseq());
        Ok(())
    }

    fn sseq(&self) -> &[u8] {
        &self.block.seq[self.sseq_span.range()]
    }

    fn xseq(&self) -> &[u8] {
        self.xseq_span
            .map_or(&[], |span| &self.block.seq[span.range()])
    }

    fn has_quality(&self) -> bool {
        self.block.header.has_qualities()
    }

    fn squal(&self) -> &[u8] {
        if self.has_quality() {
            &self.block.qual[self.sseq_span.range()]
        } else {
            &self.qbuf[..self.slen() as usize]
        }
    }

    fn xqual(&self) -> &[u8] {
        if self.has_quality()
            && let Some(span) = self.xseq_span
        {
            &self.block.qual[span.range()]
        } else {
            &self.qbuf[..self.xlen() as usize]
        }
    }
}
