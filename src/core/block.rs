use std::io;

use anyhow::{Result, bail};
use bytemuck::{cast_slice, cast_slice_mut};
use zstd::stream::{copy_decode, copy_encode};
use zstd::zstd_safe;

use super::utils::{Span, calculate_offsets, extension_read, resize_uninit, slice_and_increment};
use super::{BlockHeader, BlockRange, FileHeader, SequencingRecord};

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
    pub(crate) num_records: usize,
    /// Total nucleotides in this block
    pub(crate) nuclen: usize,
    /// Number of npos positions
    pub(crate) num_npos: usize,
    /// Current size of this block (virtual)
    current_size: usize,

    /// The file header (used for block configuration)
    ///
    /// Not to be confused with the `BlockHeader`
    pub(crate) header: FileHeader,
}
impl ColumnarBlock {
    /// Create a new columnar block with the given block size
    pub fn new(header: FileHeader) -> Self {
        Self {
            header,
            ..Default::default()
        }
    }

    fn is_empty(&self) -> bool {
        self.current_size == 0
    }

    /// Clears the internal data structures
    fn clear(&mut self) {
        // clear index counters
        {
            self.nuclen = 0;
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
        }
    }

    fn add_sequence(&mut self, record: &SequencingRecord) {
        self.l_seq.push(record.s_seq.len() as u64);
        self.seq.extend_from_slice(&record.s_seq);
        if let Some(x_seq) = record.x_seq {
            self.l_seq.push(x_seq.len() as u64);
            self.seq.extend_from_slice(x_seq);
        }

        // keep the sequence size up to date
        self.nuclen = self.seq.len();
    }

    fn add_flag(&mut self, record: &SequencingRecord) {
        record.flag.map(|flag| self.flags.push(flag));
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

    pub(crate) fn can_fit(&self, record: &SequencingRecord<'_>) -> bool {
        self.current_size + record.size() <= self.header.block_size as usize
    }

    pub fn push(&mut self, record: SequencingRecord) -> Result<()> {
        if !self.can_fit(&record) {
            bail!("Block is full")
        }

        self.add_sequence(&record);
        self.add_flag(&record);
        self.add_headers(&record);
        self.add_quality(&record);
        if record.is_paired() {
            self.num_records += 2;
        } else {
            self.num_records += 1;
        }
        self.current_size += record.size();

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
    fn fill_npos(&mut self) {
        self.npos
            .extend(memchr::memchr_iter(b'N', &self.seq).map(|i| i as u64));
        self.num_npos = self.npos.len() as usize;
    }

    /// Convert all ambiguous bases back to N
    fn backfill_npos(&mut self) {
        self.npos.iter().for_each(|idx| {
            self.seq.get_mut(*idx as usize).map(|base| *base = b'N');
        });
    }

    /// Compress all native columns into compressed representation
    fn compress_columns(&mut self) -> Result<()> {
        // compress sequence lengths
        copy_encode(
            cast_slice(&self.l_seq),
            &mut self.z_seq_len,
            self.header.compression_level as i32,
        )?;

        if self.headers.len() > 0 {
            copy_encode(
                cast_slice(&self.l_headers),
                &mut self.z_header_len,
                self.header.compression_level as i32,
            )?;
        }

        // compress npos
        if self.npos.len() > 0 {
            copy_encode(
                cast_slice(&self.npos),
                &mut self.z_npos,
                self.header.compression_level as i32,
            )?;
        }

        // compress sequence
        copy_encode(
            cast_slice(&self.ebuf),
            &mut self.z_seq,
            self.header.compression_level as i32,
        )?;

        // compress flags
        if self.flags.len() > 0 {
            copy_encode(
                cast_slice(&self.flags),
                &mut self.z_flags,
                self.header.compression_level as i32,
            )?;
        }

        // compress headers
        if self.headers.len() > 0 {
            copy_encode(
                cast_slice(&self.headers),
                &mut self.z_headers,
                self.header.compression_level as i32,
            )?;
        }

        // compress quality
        if self.qual.len() > 0 {
            copy_encode(
                cast_slice(&self.qual),
                &mut self.z_qual,
                self.header.compression_level as i32,
            )?;
        }

        Ok(())
    }

    /// Decompress all columns back to native representation
    pub fn decompress_columns(&mut self) -> Result<()> {
        // decompress sequence lengths
        {
            self.l_seq.resize(self.num_records, 0);
            copy_decode(self.z_seq_len.as_slice(), cast_slice_mut(&mut self.l_seq))?;
        }

        // decompress header lengths
        if self.z_header_len.len() > 0 {
            self.l_headers.resize(self.num_records, 0);
            copy_decode(
                self.z_header_len.as_slice(),
                cast_slice_mut(&mut self.l_headers),
            )?;
        }

        // decompress npos
        if self.z_npos.len() > 0 {
            self.npos.resize(self.num_npos, 0);
            copy_decode(self.z_npos.as_slice(), cast_slice_mut(&mut self.npos))?;
        }

        // decompress sequence
        {
            self.ebuf.resize(self.ebuf_len(), 0);
            copy_decode(self.z_seq.as_slice(), cast_slice_mut(&mut self.ebuf))?;

            bitnuc::twobit::decode(&self.ebuf, self.nuclen, &mut self.seq)?;
            self.backfill_npos();
        }

        // decompress flags
        if self.z_flags.len() > 0 {
            self.flags.resize(self.num_records, 0);
            copy_decode(self.z_flags.as_slice(), cast_slice_mut(&mut self.flags))?;
        }

        // decompress headers
        if self.z_headers.len() > 0 {
            copy_decode(self.z_headers.as_slice(), &mut self.headers)?;
        }

        // decompress quality scores
        if self.z_qual.len() > 0 {
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

    pub fn flush_to<W: io::Write>(&mut self, writer: &mut W) -> Result<Option<BlockHeader>> {
        if self.is_empty() {
            return Ok(None);
        }

        // encode all sequences at once
        self.encode_sequence()?;

        // fill npos
        self.fill_npos();

        // compress each column
        self.compress_columns()?;

        // build the block header
        let header = BlockHeader::from_block(&self);
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
        self.num_npos = header.num_npos as usize;

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
        self.num_npos = header.num_npos as usize;

        let mut byte_offset = 0;

        // decompress sequence lengths
        {
            resize_uninit(&mut self.l_seq, self.num_records);
            dctx.decompress(
                cast_slice_mut(&mut self.l_seq),
                slice_and_increment(&mut byte_offset, header.len_z_seq_len, bytes),
            )
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;
        }

        // decompress header lengths
        if header.len_z_header_len > 0 {
            resize_uninit(&mut self.l_headers, self.num_records);
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
            resize_uninit(&mut self.npos, self.num_npos);
            dctx.decompress(
                cast_slice_mut(&mut self.npos),
                slice_and_increment(&mut byte_offset, header.len_z_npos, bytes),
            )
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;
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

    pub fn iter_records(&self, range: BlockRange) -> RefRecordIter<'_> {
        RefRecordIter {
            block: self,
            range,
            index: 0,
            is_paired: self.header.is_paired(),
            has_headers: self.header.has_headers(),
        }
    }
}

pub struct RefRecordIter<'a> {
    block: &'a ColumnarBlock,
    range: BlockRange,
    index: usize,
    is_paired: bool,
    has_headers: bool,
}
impl<'a> Iterator for RefRecordIter<'a> {
    type Item = RefRecord<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.block.num_records {
            None
        } else {
            let sseq_span = Span::new_u64(
                self.block.l_seq_offsets[self.index],
                self.block.l_seq[self.index],
            );
            let sheader_span = if self.has_headers {
                Some(Span::new_u64(
                    self.block.l_header_offsets[self.index],
                    self.block.l_headers[self.index],
                ))
            } else {
                None
            };
            let xseq_span = if self.is_paired {
                Some(Span::new_u64(
                    self.block.l_seq_offsets[self.index + 1],
                    self.block.l_seq[self.index + 1],
                ))
            } else {
                None
            };
            let xheader_span = if self.is_paired && self.has_headers {
                Some(Span::new_u64(
                    self.block.l_header_offsets[self.index + 1],
                    self.block.l_headers[self.index + 1],
                ))
            } else {
                None
            };

            let record = RefRecord {
                block: self.block,
                range: self.range,
                index: self.index,
                sseq_span,
                sheader_span,
                xseq_span,
                xheader_span,
            };

            self.index += 1 + self.is_paired as usize;
            Some(record)
        }
    }
}

#[derive(Clone, Copy)]
pub struct RefRecord<'a> {
    /// A reference to the block containing this record
    block: &'a ColumnarBlock,

    /// The block range
    range: BlockRange,

    /// Local index of this record within the block
    index: usize,

    /// Span of the primary sequence within the block
    sseq_span: Span,

    /// Span of the extended sequence within the block
    xseq_span: Option<Span>,

    /// Span of the primary header within the block
    sheader_span: Option<Span>,

    /// Span of the extended header within the block
    xheader_span: Option<Span>,
}
impl<'a> binseq::BinseqRecord for RefRecord<'a> {
    fn bitsize(&self) -> binseq::BitSize {
        binseq::BitSize::Two
    }

    fn index(&self) -> u64 {
        self.range.cumulative_records - (self.block.num_records + self.index) as u64
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
            &[]
        }
    }

    fn xheader(&self) -> &[u8] {
        if let Some(span) = self.xheader_span {
            &self.block.headers[span.range()]
        } else {
            &[]
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

    fn decode_s(&self, buf: &mut Vec<u8>) -> binseq::Result<()> {
        buf.extend_from_slice(self.sseq());
        Ok(())
    }

    fn decode_x(&self, buf: &mut Vec<u8>) -> binseq::Result<()> {
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
            &[]
        }
    }

    fn xqual(&self) -> &[u8] {
        if self.has_quality()
            && let Some(span) = self.xseq_span
        {
            &self.block.qual[span.range()]
        } else {
            &[]
        }
    }
}
