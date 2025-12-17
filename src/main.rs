use std::{fs, io, path::Path, sync::Arc, thread};

use anyhow::{Result, bail};
use binseq::{BinseqRecord, ParallelProcessor, ParallelReader};
use bytemuck::{Pod, Zeroable, cast_slice, checked::cast_slice_mut};
use memmap2::Mmap;
use paraseq::{Record, fastx};
use parking_lot::Mutex;
use zstd::{
    stream::{copy_decode, copy_encode},
    zstd_safe,
};

pub const FILE_MAGIC: &[u8; 7] = b"CBQFILE";
pub const BLOCK_MAGIC: &[u8; 3] = b"BLK";
pub const INDEX_MAGIC: &[u8; 8] = b"CBQINDEX";

pub const FILE_VERSION: u8 = 1;
pub const DEFAULT_BLOCK_SIZE: u64 = 1024 * 1024;
pub const DEFAULT_COMPRESSION_LEVEL: u64 = 0;

/// Records are paired
pub const PRESENCE_PAIRED: u64 = 1 << 0;

/// Records have quality scores
pub const PRESENCE_QUALITIES: u64 = 1 << 1;

/// Records have headers
pub const PRESENCE_HEADERS: u64 = 1 << 2;

/// Records have flags
pub const PRESENCE_FLAGS: u64 = 1 << 3;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Zeroable, Pod)]
#[repr(C)]
pub struct FileHeader {
    // File Type Metadata (8 bytes)
    /// File magic number
    magic: [u8; 7],
    /// File version number
    version: u8,

    // Data presence flags (8 bytes)
    /// A bitfield indicating which data fields are present in the file
    presence_flags: u64,

    // Configuration (16 bytes)
    /// compression level
    compression_level: u64,
    /// block size in bytes
    block_size: u64,

    /// Reserved for future use
    reserved: [u8; 32],
}
impl Default for FileHeader {
    fn default() -> Self {
        let mut header = Self {
            magic: *FILE_MAGIC,
            version: FILE_VERSION,
            presence_flags: 0,
            compression_level: DEFAULT_COMPRESSION_LEVEL,
            block_size: DEFAULT_BLOCK_SIZE,
            reserved: [0; 32],
        };
        header.set_headers();
        header.set_qualities();
        header
    }
}

/// Flag getters and setters
impl FileHeader {
    pub fn set_paired(&mut self) {
        self.presence_flags |= PRESENCE_PAIRED;
    }
    pub fn set_qualities(&mut self) {
        self.presence_flags |= PRESENCE_QUALITIES;
    }
    pub fn set_headers(&mut self) {
        self.presence_flags |= PRESENCE_HEADERS;
    }
    pub fn set_flags(&mut self) {
        self.presence_flags |= PRESENCE_FLAGS;
    }
    pub fn is_paired(&self) -> bool {
        self.presence_flags & PRESENCE_PAIRED != 0
    }
    pub fn has_qualities(&self) -> bool {
        self.presence_flags & PRESENCE_QUALITIES != 0
    }
    pub fn has_headers(&self) -> bool {
        self.presence_flags & PRESENCE_HEADERS != 0
    }
    pub fn has_flags(&self) -> bool {
        self.presence_flags & PRESENCE_FLAGS != 0
    }
}

impl FileHeader {
    pub fn as_bytes(&self) -> &[u8] {
        bytemuck::bytes_of(self)
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let header: Self = *bytemuck::from_bytes(bytes);
        if header.magic != *FILE_MAGIC {
            bail!("Invalid file magic")
        }
        Ok(header)
    }
}

#[derive(Clone, Default)]
pub struct ColumnarBlock {
    /// Separate columns for each data type
    seq: Vec<u8>,
    flags: Vec<u64>,
    headers: Vec<u8>,
    qual: Vec<u8>,

    /// Length of sequences for each record
    l_seq: Vec<u64>,
    /// Length of headers for each record
    l_headers: Vec<u64>,
    /// Position of all N's in the sequence
    npos: Vec<u64>,

    /// Reusable buffer for encoding sequences
    ebuf: Vec<u64>,

    // Reusable zstd compression buffer for columnar data
    z_seq_len: Vec<u8>,
    z_header_len: Vec<u8>,
    z_npos: Vec<u8>,
    z_seq: Vec<u8>,
    z_flags: Vec<u8>,
    z_headers: Vec<u8>,
    z_qual: Vec<u8>,

    // reusable offset buffers
    l_seq_offsets: Vec<u64>,
    l_header_offsets: Vec<u64>,

    /// Number of records in the block
    num_records: usize,
    /// Total nucleotides in this block
    nuclen: usize,
    /// Number of npos positions
    num_npos: usize,
    /// Current size of this block (virtual)
    current_size: usize,

    /// The file header (used for block configuration)
    ///
    /// Not to be confused with the `BlockHeader`
    header: FileHeader,
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

    fn can_fit(&self, record: &SequencingRecord<'_>) -> bool {
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
    fn decompress_columns(&mut self) -> Result<()> {
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
            self.l_seq.resize(self.num_records, 0);
            dctx.decompress(
                cast_slice_mut(&mut self.l_seq),
                slice_and_increment(&mut byte_offset, header.len_z_seq_len, bytes),
            )
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;
        }

        // decompress header lengths
        if header.len_z_header_len > 0 {
            self.l_headers.resize(self.num_records, 0);
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
            self.npos.resize(self.num_npos, 0);
            dctx.decompress(
                cast_slice_mut(&mut self.npos),
                slice_and_increment(&mut byte_offset, header.len_z_npos, bytes),
            )
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;
        }

        // decompress sequence
        {
            self.ebuf.resize(self.ebuf_len(), 0);
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
            self.flags.resize(self.num_records, 0);
            dctx.decompress(
                cast_slice_mut(&mut self.flags),
                slice_and_increment(&mut byte_offset, header.len_z_flags, bytes),
            )
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;
        }

        // decompress headers
        if header.len_z_headers > 0 {
            self.headers.resize(
                (self.l_header_offsets.last().copied().unwrap_or(0)
                    + self.l_headers.last().copied().unwrap_or(0)) as usize,
                0,
            );
            dctx.decompress(
                &mut self.headers,
                slice_and_increment(&mut byte_offset, header.len_z_headers, bytes),
            )
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;
        }

        // decompress quality scores
        if header.len_z_qual > 0 {
            self.qual.resize(self.nuclen, 0);
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
        }
    }
}

fn extension_read<R: io::Read>(reader: &mut R, dst: &mut Vec<u8>, size: usize) -> Result<()> {
    dst.resize(size, 0);
    reader.read_exact(dst)?;
    Ok(())
}

fn slice_and_increment<'a>(offset: &mut usize, len: u64, bytes: &'a [u8]) -> &'a [u8] {
    let slice = &bytes[*offset..*offset + len as usize];
    *offset += len as usize;
    slice
}

fn calculate_offsets(values: &[u64], offsets: &mut Vec<u64>) {
    offsets.clear();
    offsets.push(0);
    for i in 1..values.len() {
        offsets.push(offsets[i - 1] + values[i - 1]);
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Span {
    offset: usize,
    length: usize,
}
impl Span {
    pub fn new(offset: usize, length: usize) -> Self {
        Span { offset, length }
    }

    pub fn new_u64(offset: u64, length: u64) -> Self {
        Span {
            offset: offset as usize,
            length: length as usize,
        }
    }

    pub fn range(&self) -> std::ops::Range<usize> {
        self.offset..self.offset + self.length
    }

    pub fn len(&self) -> usize {
        self.length
    }
}

pub struct RefRecordIter<'a> {
    block: &'a ColumnarBlock,
    range: BlockRange,
    index: usize,
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
            let sheader_span = if self.block.header.has_headers() {
                Some(Span::new_u64(
                    self.block.l_header_offsets[self.index],
                    self.block.l_headers[self.index],
                ))
            } else {
                None
            };
            let xseq_span = if self.block.header.is_paired() {
                Some(Span::new_u64(
                    self.block.l_seq_offsets[self.index + 1],
                    self.block.l_seq[self.index + 1],
                ))
            } else {
                None
            };
            let xheader_span = if self.block.header.is_paired() && self.block.header.has_headers() {
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

            if self.block.header.is_paired() {
                self.index += 2;
            } else {
                self.index += 1;
            }
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
impl<'a> RefRecord<'a> {
    /// Returns the paired index of this record within the block
    ///
    /// Note: Paired records are stored sequentially but the `RefRecord` struct's index is
    /// is the index of the Pair not the index of the individual records themselves.
    pub fn p_idx(&self) -> usize {
        if self.is_paired() {
            self.index * 2
        } else {
            self.index
        }
    }
}
impl<'a> BinseqRecord for RefRecord<'a> {
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
        self.block.header.is_paired()
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

#[derive(Clone)]
struct ColumnarBlockWriter<W: io::Write> {
    /// Internal writer for the block
    inner: W,

    /// The CBQ file header
    #[allow(dead_code)]
    header: FileHeader,

    /// A reusable block for this writer
    block: ColumnarBlock,

    /// Offsets of the blocks written by this writer
    offsets: Vec<BlockHeader>,
}
impl<W: io::Write> ColumnarBlockWriter<W> {
    pub fn new(inner: W, header: FileHeader) -> Result<Self> {
        // Build the writer
        let mut writer = Self {
            inner,
            header,
            block: ColumnarBlock::new(header),
            offsets: Vec::default(),
        };

        // Ensure the header is written to the file
        writer.inner.write_all(header.as_bytes())?;

        Ok(writer)
    }

    pub fn push(&mut self, record: SequencingRecord) -> Result<()> {
        if !self.block.can_fit(&record) {
            self.flush()?;
        }
        self.block.push(record)?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.block.flush_to(&mut self.inner)?.map(|header| {
            self.offsets.push(header);
        });
        Ok(())
    }

    pub fn finish(&mut self) -> Result<()> {
        self.flush()?;
        self.write_index()?;
        Ok(())
    }

    fn write_index(&mut self) -> Result<()> {
        let index = Index::from_block_headers(&self.offsets);
        let z_index = index.encoded()?;
        let header = IndexHeader::new(index.size(), z_index.len() as u64);
        let footer = IndexFooter::new(z_index.len() as u64);

        // Write the index to the inner writer
        {
            self.inner.write_all(header.as_bytes())?;
            self.inner.write_all(&z_index)?;
            self.inner.write_all(footer.as_bytes())?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, Zeroable, Pod)]
#[repr(C)]
struct IndexHeader {
    /// Magic number identifying the index format
    magic: [u8; 8],

    /// Number of bytes in the uncompressed index
    u_bytes: u64,

    /// Number of bytes in the compressed index
    z_bytes: u64,
}
impl IndexHeader {
    /// Creates a new index header
    pub fn new(u_bytes: u64, z_bytes: u64) -> Self {
        Self {
            magic: *INDEX_MAGIC,
            u_bytes,
            z_bytes,
        }
    }

    pub fn as_bytes(&self) -> &[u8] {
        bytemuck::bytes_of(self)
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let header: Self = *bytemuck::from_bytes(bytes);
        if header.magic != *INDEX_MAGIC {
            bail!("Invalid index header magic");
        }
        Ok(header)
    }
}

#[derive(Debug, Clone, Copy, Zeroable, Pod)]
#[repr(C)]
struct IndexFooter {
    /// Number of bytes in the compressed index
    bytes: u64,

    /// Magic number identifying the index format
    magic: [u8; 8],
}

impl IndexFooter {
    /// Creates a new index footer
    pub fn new(bytes: u64) -> Self {
        Self {
            bytes,
            magic: *INDEX_MAGIC,
        }
    }
    pub fn as_bytes(&self) -> &[u8] {
        bytemuck::bytes_of(self)
    }
    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let footer: Self = *bytemuck::from_bytes(bytes);
        if footer.magic != *INDEX_MAGIC {
            bail!("Invalid index footer magic");
        }
        Ok(footer)
    }
}

/// An index of block ranges for quick lookups
#[derive(Clone)]
pub struct Index {
    ranges: Vec<BlockRange>,
}
impl Index {
    /// Builds the index from a list of block headers
    pub fn from_block_headers(block_headers: &[BlockHeader]) -> Self {
        let mut offset = size_of::<FileHeader>() as u64;
        let mut cumulative_records = 0;
        let mut ranges = Vec::default();
        for block_header in block_headers {
            let range = BlockRange::new(offset, cumulative_records + block_header.num_records);
            offset += (size_of::<BlockHeader>() + block_header.block_len()) as u64;
            cumulative_records += block_header.num_records;
            ranges.push(range);
        }
        Self { ranges }
    }

    /// Returns the byte representation of the index
    pub fn as_bytes(&self) -> &[u8] {
        bytemuck::cast_slice(&self.ranges)
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let ranges = match bytemuck::try_cast_slice(bytes) {
            Ok(ranges) => ranges.to_vec(),
            Err(_) => bail!("Failed to cast bytes to Index"),
        };
        Ok(Self { ranges })
    }

    /// Returns the size of the index in bytes
    pub fn size(&self) -> u64 {
        self.as_bytes().len() as u64
    }

    /// Encodes the index into a ZSTD-compressed byte array
    pub fn encoded(&self) -> Result<Vec<u8>> {
        let mut encoded = Vec::default();
        copy_encode(self.as_bytes(), &mut encoded, 0)?;
        Ok(encoded)
    }

    /// Returns the number of records in the index
    pub fn num_records(&self) -> usize {
        self.ranges
            .last()
            .map_or(0, |range| range.cumulative_records as usize)
    }

    /// Returns the number of blocks in the index
    pub fn num_blocks(&self) -> usize {
        self.ranges.len()
    }

    pub fn iter_blocks(&self) -> BlockIter<'_> {
        BlockIter {
            index: self,
            pos: 0,
        }
    }
}

pub struct BlockIter<'a> {
    index: &'a Index,
    pos: usize,
}
impl Iterator for BlockIter<'_> {
    type Item = BlockRange;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.index.num_blocks() {
            None
        } else {
            let block = self.index.ranges[self.pos];
            self.pos += 1;
            Some(block)
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Zeroable, Pod, Default)]
#[repr(C)]
pub struct BlockRange {
    /// Byte offset of this block
    offset: u64,

    /// Number of records up to and including this block
    cumulative_records: u64,
}
impl BlockRange {
    pub fn new(offset: u64, cumulative_records: u64) -> Self {
        Self {
            offset,
            cumulative_records,
        }
    }
}

#[derive(Copy, Clone, Pod, Zeroable, Debug, PartialEq, Eq, Hash)]
#[repr(C)]
pub struct BlockHeader {
    magic: [u8; 3],
    version: u8,
    padding: [u8; 4],

    // length of compressed columns
    len_z_seq_len: u64,
    len_z_header_len: u64,
    len_z_npos: u64,
    len_z_seq: u64,
    len_z_flags: u64,
    len_z_headers: u64,
    len_z_qual: u64,

    // full decoded length of the sequence block
    nuclen: u64,

    // number of npos positions
    num_npos: u64,

    // number of records in the block
    num_records: u64,
}
impl BlockHeader {
    pub fn from_block(block: &ColumnarBlock) -> Self {
        Self {
            magic: *BLOCK_MAGIC,
            version: 1,
            padding: [42; 4],
            len_z_seq_len: block.z_seq_len.len() as u64,
            len_z_header_len: block.z_header_len.len() as u64,
            len_z_npos: block.z_npos.len() as u64,
            len_z_seq: block.z_seq.len() as u64,
            len_z_flags: block.z_flags.len() as u64,
            len_z_headers: block.z_headers.len() as u64,
            len_z_qual: block.z_qual.len() as u64,
            nuclen: block.nuclen as u64,
            num_npos: block.num_npos as u64,
            num_records: block.num_records as u64,
        }
    }

    /// Calculate the length of the block in bytes.
    #[allow(dead_code)]
    pub fn block_len(&self) -> usize {
        (self.len_z_seq_len
            + self.len_z_header_len
            + self.len_z_npos
            + self.len_z_seq
            + self.len_z_flags
            + self.len_z_headers
            + self.len_z_qual) as usize
    }

    pub fn as_bytes(&self) -> &[u8] {
        bytemuck::bytes_of(self)
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let header: Self = *bytemuck::from_bytes(bytes);
        if header.magic != *BLOCK_MAGIC {
            bail!("Invalid Block Header found")
        }
        Ok(header)
    }

    pub fn write<W: io::Write>(&self, writer: &mut W) -> Result<(), std::io::Error> {
        writer.write_all(self.as_bytes())
    }
}

#[derive(Clone, Default)]
pub struct SequencingRecord<'a> {
    s_seq: &'a [u8],
    s_qual: Option<&'a [u8]>,
    s_header: Option<&'a [u8]>,
    x_seq: Option<&'a [u8]>,
    x_qual: Option<&'a [u8]>,
    x_header: Option<&'a [u8]>,
    flag: Option<u64>,
}
impl<'a> SequencingRecord<'a> {
    pub fn new(
        s_seq: &'a [u8],
        s_qual: Option<&'a [u8]>,
        s_header: Option<&'a [u8]>,
        x_seq: Option<&'a [u8]>,
        x_qual: Option<&'a [u8]>,
        x_header: Option<&'a [u8]>,
        flag: Option<u64>,
    ) -> Self {
        Self {
            s_seq,
            s_qual,
            s_header,
            x_seq,
            x_qual,
            x_header,
            flag,
        }
    }

    /// Returns the size of the record in bytes
    pub fn size(&self) -> usize {
        self.s_seq.len()
            + self.s_qual.map_or(0, |q| q.len())
            + self.s_header.map_or(0, |h| h.len())
            + self.x_seq.map_or(0, |q| q.len())
            + self.x_qual.map_or(0, |q| q.len())
            + self.x_header.map_or(0, |h| h.len())
            + self.flag.map_or(0, |f| f.to_le_bytes().len())
    }

    pub fn is_paired(&self) -> bool {
        self.x_seq.is_some()
    }
}

pub struct Reader<R: io::Read> {
    inner: R,
    block: ColumnarBlock,
    iheader: Option<IndexHeader>,
}
impl<R: io::Read> Reader<R> {
    pub fn new(mut inner: R) -> Result<Self> {
        let mut header_buf = [0u8; size_of::<FileHeader>()];
        inner.read_exact(&mut header_buf)?;
        let header = FileHeader::from_bytes(&header_buf)?;

        Ok(Self {
            inner,
            block: ColumnarBlock::new(header),
            iheader: None,
        })
    }

    pub fn read_block(&mut self) -> Result<bool> {
        let mut iheader_buf = [0u8; size_of::<IndexHeader>()];
        let mut diff_buf = [0u8; size_of::<BlockHeader>() - size_of::<IndexHeader>()];
        let mut header_buf = [0u8; size_of::<BlockHeader>()];

        // Attempt to read the index header
        match self.inner.read_exact(&mut iheader_buf) {
            Ok(_) => {}
            Err(e) => {
                if e.kind() == io::ErrorKind::UnexpectedEof {
                    return Ok(false);
                } else {
                    return Err(e.into());
                }
            }
        }

        // The stream is exhausted, no more blocks to read
        if let Ok(iheader) = IndexHeader::from_bytes(&iheader_buf) {
            self.iheader = Some(iheader);
            return Ok(false);
        } else {
            // attempt to read the rest of the block header
            match self.inner.read_exact(&mut diff_buf) {
                Ok(_) => {}
                Err(e) => {
                    if e.kind() == io::ErrorKind::UnexpectedEof {
                        return Ok(false);
                    } else {
                        return Err(e.into());
                    }
                }
            }
            header_buf[..iheader_buf.len()].copy_from_slice(&iheader_buf);
            header_buf[iheader_buf.len()..].copy_from_slice(&diff_buf);
        }

        let header = BlockHeader::from_bytes(&header_buf)?;
        self.block.read_from(&mut self.inner, header)?;

        Ok(true)
    }

    pub fn read_index(&mut self) -> Result<Option<Index>> {
        let Some(header) = self.iheader else {
            return Ok(None);
        };
        let mut z_index_buf = Vec::new();
        let mut index_buf = Vec::new();
        let mut footer_buf = [0u8; size_of::<IndexFooter>()];

        // Read the index data from the reader
        z_index_buf.resize(header.z_bytes as usize, 0);

        // Reads the compressed index data
        self.inner.read_exact(&mut z_index_buf)?;
        copy_decode(z_index_buf.as_slice(), &mut index_buf)?;
        let index = Index::from_bytes(&index_buf)?;

        // Read the footer data from the reader
        self.inner.read_exact(&mut footer_buf)?;
        let _footer = IndexFooter::from_bytes(&footer_buf)?;

        Ok(Some(index))
    }
}

pub struct MmapReader {
    inner: Arc<Mmap>,
    index: Arc<Index>,

    /// Reusable record block
    block: ColumnarBlock,

    /// Reusable decompression context
    dctx: zstd_safe::DCtx<'static>,
}
impl Clone for MmapReader {
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
            index: self.index.clone(),
            block: self.block.clone(),
            dctx: zstd_safe::DCtx::create(),
        }
    }
}
impl MmapReader {
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = fs::File::open(path)?;

        // Load the mmap
        let inner = unsafe { Mmap::map(&file) }?;

        // Build the header
        let header = FileHeader::from_bytes(&inner[..size_of::<FileHeader>()])?;

        // build the index
        let index = {
            // Load the index footer
            let footer_start = inner.len() - size_of::<IndexFooter>();
            let index_footer = IndexFooter::from_bytes(&inner[footer_start..])?;

            // Find the coordinates of the compressed index
            let z_index_start = footer_start - index_footer.bytes as usize;
            let z_index_slice = &inner[z_index_start..footer_start];

            // Decompress the index
            let mut index_buf = Vec::default();
            copy_decode(z_index_slice, &mut index_buf)?;

            // Load the index
            Index::from_bytes(&index_buf)
        }?;

        Ok(Self {
            inner: Arc::new(inner),
            index: Arc::new(index),
            block: ColumnarBlock::new(header),
            dctx: zstd_safe::DCtx::create(),
        })
    }

    pub fn num_records(&self) -> usize {
        self.index.num_records()
    }

    fn load_block(&mut self, range: BlockRange) -> Result<()> {
        let header_start = range.offset as usize;
        let header_end = size_of::<BlockHeader>() + header_start;
        let block_header = {
            let mut block_header_buf = [0u8; size_of::<BlockHeader>()];
            block_header_buf.copy_from_slice(&self.inner[header_start..header_end]);
            BlockHeader::from_bytes(&block_header_buf)
        }?;

        let data_end = header_end + block_header.block_len();
        let block_data_slice = &self.inner[header_end..data_end];
        self.block
            .decompress_from_bytes(&block_data_slice, block_header, &mut self.dctx)?;
        Ok(())
    }
}
impl ParallelReader for MmapReader {
    fn process_parallel<P: ParallelProcessor + Clone + 'static>(
        self,
        processor: P,
        num_threads: usize,
    ) -> binseq::Result<()> {
        let num_threads = if num_threads == 0 {
            num_cpus::get()
        } else {
            num_threads.min(num_cpus::get())
        };

        let blocks_per_thread = self.index.num_blocks().div_ceil(num_threads);

        let mut handles = Vec::new();
        for thread_id in 0..num_threads {
            let start_block_idx = thread_id * blocks_per_thread;
            let end_block_idx = ((thread_id + 1) * blocks_per_thread).min(self.index.num_blocks());
            let mut t_reader = self.clone();
            let mut t_proc = processor.clone();

            let thread_handle = thread::spawn(move || -> binseq::Result<()> {
                // Pull all block ranges for this thread
                let t_block_ranges: Vec<_> = t_reader
                    .index
                    .iter_blocks()
                    .skip(start_block_idx)
                    .take(end_block_idx - start_block_idx)
                    .collect();

                for range in t_block_ranges {
                    t_reader.load_block(range)?;
                    for record in t_reader.block.iter_records(range) {
                        t_proc.process_record(record)?;
                    }
                    t_proc.on_batch_complete()?;
                }

                Ok(())
            });
            handles.push(thread_handle);
        }

        for handle in handles {
            handle.join().unwrap()?;
        }
        Ok(())
    }

    fn process_parallel_range<P: ParallelProcessor + Clone + 'static>(
        self,
        _processor: P,
        _num_threads: usize,
        _range: std::ops::Range<usize>,
    ) -> binseq::Result<()> {
        unimplemented!()
    }
}

fn write_file(ipath: &str, opath: &str) -> Result<()> {
    let handle = io::BufWriter::new(fs::File::create(opath)?);
    let header = FileHeader::default();
    let mut writer = ColumnarBlockWriter::new(handle, header)?;

    let mut reader = fastx::Reader::from_path(ipath)?;
    let mut rset = reader.new_record_set();
    while rset.fill(&mut reader)? {
        for res in rset.iter() {
            let record = res?;
            let seq = record.seq();
            let ref_record = SequencingRecord::new(
                &seq,
                record.qual(),
                Some(record.id()),
                None,
                None,
                None,
                None,
            );
            writer.push(ref_record)?;
        }
    }
    writer.finish()?;

    Ok(())
}

fn read_file(ipath: &str) -> Result<()> {
    let rhandle = fs::File::open(ipath).map(io::BufReader::new)?;
    let mut reader = Reader::new(rhandle)?;

    while reader.read_block()? {
        reader.block.decompress_columns()?;
    }
    reader.read_index()?;
    Ok(())
}

type BoxedWriter = Box<dyn io::Write + Send>;

#[derive(Clone)]
pub struct Processor {
    l_records: usize,
    l_buf: Vec<u8>,

    records: Arc<Mutex<usize>>,
    writer: Arc<Mutex<BoxedWriter>>,
}
impl Processor {
    pub fn new(writer: BoxedWriter) -> Self {
        Self {
            l_records: 0,
            l_buf: Vec::new(),
            records: Arc::new(Mutex::new(0)),
            writer: Arc::new(Mutex::new(writer)),
        }
    }

    pub fn n_records(&self) -> usize {
        *self.records.lock()
    }
}
impl ParallelProcessor for Processor {
    fn process_record<R: BinseqRecord>(&mut self, record: R) -> binseq::Result<()> {
        write_fastq(
            &mut self.l_buf,
            record.sheader(),
            record.sseq(),
            record.squal(),
        )?;
        if record.is_paired() {
            write_fastq(
                &mut self.l_buf,
                record.xheader(),
                record.xseq(),
                record.xqual(),
            )?;
        }
        self.l_records += 1;
        Ok(())
    }
    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        {
            let mut writer = self.writer.lock();
            writer.write_all(&self.l_buf)?;
            writer.flush()?;
        }
        self.l_buf.clear();

        *self.records.lock() += self.l_records;
        self.l_records = 0;
        Ok(())
    }
}

fn write_fastq<W: io::Write>(writer: &mut W, header: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(header)?;
    writer.write_all(b"\n")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n")?;
    writer.write_all(b"+\n")?;
    writer.write_all(qual)?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn read_mmap(ipath: &str) -> Result<()> {
    let reader = MmapReader::new(ipath)?;
    let proc = Processor::new(Box::new(io::stdout()));
    reader.process_parallel(proc.clone(), 0)?;
    println!("Number of records: {}", proc.n_records());
    Ok(())
}

fn main() -> Result<()> {
    let ipath = "./data/some.fq.gz";
    let opath = "./data/some.cbq";

    eprintln!("Writing file {} - reading from {}", opath, ipath);
    write_file(ipath, opath)?;

    eprintln!("Reading file {}", opath);
    read_file(opath)?;

    eprintln!("Reading file {} using memory mapping", opath);
    read_mmap(opath)?;

    Ok(())
}
