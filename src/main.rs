use std::{fs, io};

use anyhow::Result;
use bytemuck::{Pod, Zeroable, cast_slice};
use paraseq::{Record, fastx};
use zstd::stream::copy_encode;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Pod, Zeroable)]
#[repr(C)]
struct Span {
    offset: u64,
    length: u64,
}
impl Span {
    pub fn new(offset: usize, length: usize) -> Self {
        Self {
            offset: offset as u64,
            length: length as u64,
        }
    }

    pub fn from_extension<T: Clone>(buffer: &mut Vec<T>, addition: &[T]) -> Self {
        let offset = buffer.len();
        let length = addition.len();
        buffer.extend_from_slice(addition);
        Self::new(offset, length)
    }

    pub fn from_push<T: Clone>(buffer: &mut Vec<T>, item: T) -> Self {
        let offset = buffer.len();
        buffer.push(item);
        Self::new(offset, 1)
    }
}

#[derive(Clone)]
struct ColumnarBlock<W: io::Write> {
    /// Internal writer for the block
    inner: W,

    /// Separate columns for each data type
    seq: Vec<u8>,
    flags: Vec<u64>,
    headers: Vec<u8>,
    qual: Vec<u8>,

    /// Record spans for columns
    seq_spans: Vec<Span>,
    flag_spans: Vec<Span>,
    header_spans: Vec<Span>,
    qual_spans: Vec<Span>,

    /// Reusable buffer for encoding sequences
    ebuf: Vec<u64>,

    /// Reusable zstd compression buffer for columnar data
    z_seq: Vec<u8>,
    z_flags: Vec<u8>,
    z_headers: Vec<u8>,
    z_qual: Vec<u8>,

    /// Total nucleotides in this block
    nuclen: usize,
    /// Current size of this block (virtual)
    current_size: usize,
    /// Maximum size of this block (virtual)
    block_size: usize,
}
impl<W: io::Write> ColumnarBlock<W> {
    pub fn new(inner: W, block_size: usize) -> Self {
        Self {
            inner,
            block_size,
            current_size: 0,
            nuclen: 0,
            seq: Vec::default(),
            flags: Vec::default(),
            headers: Vec::default(),
            qual: Vec::default(),
            flag_spans: Vec::default(),
            seq_spans: Vec::default(),
            header_spans: Vec::default(),
            qual_spans: Vec::default(),
            ebuf: Vec::default(),
            z_seq: Vec::default(),
            z_flags: Vec::default(),
            z_headers: Vec::default(),
            z_qual: Vec::default(),
        }
    }

    fn clear(&mut self) {
        self.nuclen = 0;
        self.current_size = 0;

        // clear spans
        {
            self.seq_spans.clear();
            self.flag_spans.clear();
            self.header_spans.clear();
            self.qual_spans.clear();
        }

        // clear vectors
        {
            self.seq.clear();
            self.flags.clear();
            self.headers.clear();
            self.qual.clear();
        }

        // clear encodings
        {
            self.ebuf.clear();
            self.z_seq.clear();
            self.z_flags.clear();
            self.z_headers.clear();
            self.z_qual.clear();
        }
    }

    fn add_sequence(&mut self, record: &SequencingRecord) {
        self.seq_spans
            .push(Span::from_extension(&mut self.seq, record.s_seq));
        if let Some(x_seq) = record.x_seq {
            self.seq_spans
                .push(Span::from_extension(&mut self.seq, x_seq));
        }

        // keep the sequence size up to date
        self.nuclen = self.seq.len();
    }

    fn add_flag(&mut self, record: &SequencingRecord) {
        if let Some(flag) = record.flag {
            self.flag_spans.push(Span::from_push(&mut self.flags, flag));
        }
    }

    fn add_headers(&mut self, record: &SequencingRecord) {
        if let Some(header) = record.s_header {
            self.header_spans
                .push(Span::from_extension(&mut self.headers, header));
        }
        if let Some(header) = record.x_header {
            self.header_spans
                .push(Span::from_extension(&mut self.headers, header));
        }
    }

    fn add_quality(&mut self, record: &SequencingRecord) {
        if let Some(qual) = record.s_qual {
            self.qual_spans
                .push(Span::from_extension(&mut self.qual, qual));
        }
        if let Some(qual) = record.x_qual {
            self.qual_spans
                .push(Span::from_extension(&mut self.qual, qual));
        }
    }

    pub fn push(&mut self, record: SequencingRecord) -> Result<()> {
        if self.current_size + record.size() > self.block_size {
            self.flush()?;
        }

        self.add_sequence(&record);
        self.add_flag(&record);
        self.add_headers(&record);
        self.add_quality(&record);
        self.current_size += record.size();

        Ok(())
    }

    fn encode_sequence(&mut self) -> Result<()> {
        bitnuc::twobit::encode(&self.seq, &mut self.ebuf)?;
        Ok(())
    }

    fn compress_columns(&mut self) -> Result<()> {
        // compress sequence
        copy_encode(cast_slice(&self.ebuf), &mut self.z_seq, 3)?;

        // compress flags
        if self.flags.len() > 0 {
            copy_encode(cast_slice(&self.flags), &mut self.z_flags, 3)?;
        }

        // compress headers
        if self.headers.len() > 0 {
            copy_encode(cast_slice(&self.headers), &mut self.z_headers, 3)?;
        }

        // compress quality
        if self.qual.len() > 0 {
            copy_encode(cast_slice(&self.qual), &mut self.z_qual, 3)?;
        }

        Ok(())
    }

    fn write_to_inner(&mut self) -> Result<()> {
        // write all spans
        {
            self.inner
                .write_all(bytemuck::cast_slice(&self.seq_spans))?;
            self.inner
                .write_all(bytemuck::cast_slice(&self.flag_spans))?;
            self.inner
                .write_all(bytemuck::cast_slice(&self.header_spans))?;
            self.inner
                .write_all(bytemuck::cast_slice(&self.qual_spans))?;
        }

        // write all compressed buffers
        {
            self.inner.write_all(&self.z_seq)?;
            self.inner.write_all(&self.z_flags)?;
            self.inner.write_all(&self.z_headers)?;
            self.inner.write_all(&self.z_qual)?;
        }

        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        eprintln!("Flushing block!");

        // encode all sequences at once
        self.encode_sequence()?;

        // compress each column
        self.compress_columns()?;

        // build the block header
        let header = BlockHeader::from_writer_state(&self);

        // write the block header
        header.write(&mut self.inner)?;

        // write the internal state to the inner writer
        self.write_to_inner()?;

        // clear the internal state
        self.clear();

        Ok(())
    }
}

#[derive(Copy, Clone, Pod, Zeroable, Debug, PartialEq, Eq, Hash)]
#[repr(C)]
struct BlockHeader {
    magic: [u8; 3],
    version: u8,
    padding: [u8; 4],

    // length of spans
    len_span_seq: u32,
    len_span_flags: u32,
    len_span_headers: u32,
    len_span_qual: u32,

    // length of compressed columns
    len_z_seq: u32,
    len_z_flags: u32,
    len_z_headers: u32,
    len_z_qual: u32,

    // full decoded length of the sequence block
    nuclen: u64,
}
impl BlockHeader {
    pub fn from_writer_state<W: io::Write>(writer: &ColumnarBlock<W>) -> Self {
        Self {
            magic: *b"CBQ",
            version: 1,
            padding: [42; 4],
            len_span_seq: writer.seq_spans.len() as u32,
            len_span_flags: writer.flag_spans.len() as u32,
            len_span_headers: writer.header_spans.len() as u32,
            len_span_qual: writer.qual_spans.len() as u32,
            len_z_seq: writer.z_seq.len() as u32,
            len_z_flags: writer.z_flags.len() as u32,
            len_z_headers: writer.z_headers.len() as u32,
            len_z_qual: writer.z_qual.len() as u32,
            nuclen: writer.nuclen as u64,
        }
    }

    pub fn as_bytes(&self) -> &[u8] {
        bytemuck::bytes_of(self)
    }

    pub fn write<W: io::Write>(&self, writer: &mut W) -> Result<(), std::io::Error> {
        writer.write_all(self.as_bytes())
    }
}

#[derive(Clone, Default)]
struct SequencingRecord<'a> {
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
}

fn main() -> Result<()> {
    let path = "./data/some.fq.gz";
    let opath = "./data/some.cbq";
    let handle = io::BufWriter::new(fs::File::create(opath)?);
    let mut writer = ColumnarBlock::new(handle, 1024 * 1024);

    let mut reader = fastx::Reader::from_path(path)?;
    let mut rset = reader.new_record_set();
    while rset.fill(&mut reader)? {
        for res in rset.iter() {
            let record = res?;
            let seq = record.seq();
            let ref_record = SequencingRecord::new(&seq, None, None, None, None, None, None);
            writer.push(ref_record)?;
        }
    }
    writer.flush()?;
    Ok(())
}
