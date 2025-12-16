use std::{io, sync::Arc};

use anyhow::Result;
use bytemuck::cast_slice;
use paraseq::fastx;
use parking_lot::Mutex;
use zstd::stream::copy_encode;

struct Encoder<W: io::Write> {
    writer: Arc<Mutex<W>>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
struct Span {
    offset: usize,
    length: usize,
}
impl Span {
    pub fn new(offset: usize, length: usize) -> Self {
        Self { offset, length }
    }

    pub fn from_extension<T: Clone>(buffer: &mut Vec<T>, addition: &[T]) -> Self {
        let offset = buffer.len();
        let length = addition.len();
        buffer.extend_from_slice(addition);
        Self { offset, length }
    }

    pub fn from_push<T: Clone>(buffer: &mut Vec<T>, item: T) -> Self {
        let offset = buffer.len();
        buffer.push(item);
        Self { offset, length: 1 }
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

    /// Maximum size of this block (virtual)
    seq_size: usize,
    current_size: usize,
    block_size: usize,
}
impl<W: io::Write> ColumnarBlock<W> {
    pub fn new(inner: W, block_size: usize) -> Self {
        Self {
            inner,
            block_size,
            current_size: 0,
            seq_size: 0,
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

    fn add_sequence(&mut self, record: &SequencingRecord) {
        self.seq_spans
            .push(Span::from_extension(&mut self.seq, record.s_seq));
        if let Some(x_seq) = record.x_seq {
            self.seq_spans
                .push(Span::from_extension(&mut self.seq, x_seq));
        }

        // keep the sequence size up to date
        self.seq_size = self.seq.len();
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

    pub fn push(&mut self, record: &SequencingRecord) -> Result<()> {
        if self.current_size + record.size() > self.block_size {
            self.flush()?;
        }

        self.add_sequence(record);
        self.add_flag(record);
        self.add_headers(record);
        self.add_quality(record);

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

    pub fn flush(&mut self) -> Result<()> {
        // encode all sequences at once
        self.encode_sequence()?;

        // compress each column
        self.compress_columns()?;

        Ok(())
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
    let path = "./data/some.fq";
    let reader = fastx::Reader::from_path(path)?;
    Ok(())
}
