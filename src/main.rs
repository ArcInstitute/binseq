use std::{fs, io};

use anyhow::Result;
use bytemuck::{Pod, Zeroable, cast_slice};
use paraseq::{Record, fastx};
use zstd::stream::copy_encode;

#[derive(Clone)]
struct ColumnarBlock<W: io::Write> {
    /// Internal writer for the block
    inner: W,

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

    /// Reusable zstd compression buffer for columnar data
    z_seq_len: Vec<u8>,
    z_header_len: Vec<u8>,
    z_npos: Vec<u8>,
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

            // data buffers
            seq: Vec::default(),
            flags: Vec::default(),
            headers: Vec::default(),
            qual: Vec::default(),

            // span buffer
            l_seq: Vec::default(),
            l_headers: Vec::default(),

            // Position of all N's in the sequence
            npos: Vec::default(),

            // encoding buffer
            ebuf: Vec::default(),

            // compression buffers
            z_seq_len: Vec::default(),
            z_header_len: Vec::default(),
            z_npos: Vec::default(),
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
            self.l_seq.clear();
            self.l_headers.clear();
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
        bitnuc::twobit::encode_with_invalid(&self.seq, &mut self.ebuf)?;
        Ok(())
    }

    fn fill_npos(&mut self) {
        self.npos
            .extend(memchr::memchr_iter(b'N', &self.seq).map(|i| i as u64))
    }

    fn compress_columns(&mut self) -> Result<()> {
        // compress sequence lengths
        copy_encode(cast_slice(&self.l_seq), &mut self.z_seq_len, 0)?;

        if self.headers.len() > 0 {
            copy_encode(cast_slice(&self.l_headers), &mut self.z_header_len, 0)?;
        }

        // compress npos
        if self.npos.len() > 0 {
            copy_encode(cast_slice(&self.npos), &mut self.z_npos, 0)?;
        }

        // compress sequence
        copy_encode(cast_slice(&self.ebuf), &mut self.z_seq, 0)?;

        // compress flags
        if self.flags.len() > 0 {
            copy_encode(cast_slice(&self.flags), &mut self.z_flags, 0)?;
        }

        // compress headers
        if self.headers.len() > 0 {
            copy_encode(cast_slice(&self.headers), &mut self.z_headers, 0)?;
        }

        // compress quality
        if self.qual.len() > 0 {
            copy_encode(cast_slice(&self.qual), &mut self.z_qual, 0)?;
        }

        Ok(())
    }

    fn write_to_inner(&mut self) -> Result<()> {
        self.inner.write_all(&self.z_seq_len)?;
        self.inner.write_all(&self.z_header_len)?;
        self.inner.write_all(&self.z_npos)?;
        self.inner.write_all(&self.z_seq)?;
        self.inner.write_all(&self.z_flags)?;
        self.inner.write_all(&self.z_headers)?;
        self.inner.write_all(&self.z_qual)?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        // encode all sequences at once
        self.encode_sequence()?;

        // fill npos
        self.fill_npos();

        // compress each column
        self.compress_columns()?;

        // build the block header
        let header = BlockHeader::from_writer_state(&self);
        eprintln!("{header:?}");

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
}
impl BlockHeader {
    pub fn from_writer_state<W: io::Write>(writer: &ColumnarBlock<W>) -> Self {
        Self {
            magic: *b"CBQ",
            version: 1,
            padding: [42; 4],
            len_z_seq_len: writer.z_seq_len.len() as u64,
            len_z_header_len: writer.z_header_len.len() as u64,
            len_z_npos: writer.z_npos.len() as u64,
            len_z_seq: writer.z_seq.len() as u64,
            len_z_flags: writer.z_flags.len() as u64,
            len_z_headers: writer.z_headers.len() as u64,
            len_z_qual: writer.z_qual.len() as u64,
            nuclen: writer.nuclen as u64,
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
            + self.len_z_qual
            + self.nuclen) as usize
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
    writer.flush()?;
    Ok(())
}
