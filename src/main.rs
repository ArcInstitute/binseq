use std::{fs, io};

use anyhow::{Result, bail};
use bytemuck::{Pod, Zeroable, cast_slice, checked::cast_slice_mut};
use paraseq::{Record, fastx};
use zstd::stream::{copy_decode, copy_encode};

#[derive(Clone, Default)]
struct ColumnarBlock {
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

    /// Length of the encoded sequence
    ebuf_len: usize,
    /// Number of records in the block
    num_records: usize,
    /// Total nucleotides in this block
    nuclen: usize,
    /// Number of npos positions
    num_npos: usize,
    /// Current size of this block (virtual)
    current_size: usize,
    /// Maximum size of this block (virtual)
    block_size: usize,
}
impl ColumnarBlock {
    /// Create a new columnar block with the given block size
    pub fn new(block_size: usize) -> Self {
        Self {
            block_size,
            ..Default::default()
        }
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
        self.current_size + record.size() <= self.block_size
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

    /// Encode the sequence into a compressed representation
    fn encode_sequence(&mut self) -> Result<()> {
        bitnuc::twobit::encode_with_invalid(&self.seq, &mut self.ebuf)?;
        self.ebuf_len = self.ebuf.len();
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
            self.ebuf.resize(self.ebuf_len, 0);
            copy_decode(self.z_seq.as_slice(), cast_slice_mut(&mut self.ebuf))?;

            self.seq.resize(self.nuclen, 0);
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

    pub fn flush_to<W: io::Write>(&mut self, writer: &mut W) -> Result<()> {
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

        Ok(())
    }

    pub fn read_from<R: io::Read>(&mut self, reader: &mut R, header: BlockHeader) -> Result<()> {
        // clears the internal state
        self.clear();

        // reload the internal state from the reader
        self.nuclen = header.nuclen as usize;
        self.num_records = header.num_records as usize;
        self.num_npos = header.num_npos as usize;
        self.ebuf_len = header.ebuf_len as usize;

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
}

fn extension_read<R: io::Read>(reader: &mut R, dst: &mut Vec<u8>, size: usize) -> Result<()> {
    dst.resize(size, 0);
    reader.read_exact(dst)?;
    Ok(())
}

#[derive(Clone)]
struct ColumnarBlockWriter<W: io::Write> {
    /// Internal writer for the block
    inner: W,

    /// The block for this writer
    block: ColumnarBlock,
}
impl<W: io::Write> ColumnarBlockWriter<W> {
    pub fn new(inner: W, block_size: usize) -> Self {
        Self {
            inner,
            block: ColumnarBlock::new(block_size),
        }
    }

    pub fn push(&mut self, record: SequencingRecord) -> Result<()> {
        if !self.block.can_fit(&record) {
            self.block.flush_to(&mut self.inner)?;
        }
        self.block.push(record)?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.block.flush_to(&mut self.inner)?;
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

    // length of the encoded sequence
    ebuf_len: u64,

    // number of npos positions
    num_npos: u64,

    // number of records in the block
    num_records: u64,
}
impl BlockHeader {
    pub fn from_block(block: &ColumnarBlock) -> Self {
        Self {
            magic: *b"CBQ",
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
            ebuf_len: block.ebuf.len() as u64,
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

    pub fn from_bytes(bytes: &[u8]) -> Self {
        *bytemuck::from_bytes(bytes)
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

    pub fn is_paired(&self) -> bool {
        self.x_seq.is_some()
    }
}

pub struct Reader<R: io::Read> {
    inner: R,
    block: ColumnarBlock,
}
impl<R: io::Read> Reader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            block: ColumnarBlock::new(0),
        }
    }

    pub fn read_block(&mut self) -> Result<bool> {
        let mut header_buf = [0u8; size_of::<BlockHeader>()];

        // Read the block header from the reader
        match self.inner.read_exact(&mut header_buf) {
            Ok(_) => {}
            Err(e) => {
                if e.kind() == io::ErrorKind::UnexpectedEof {
                    return Ok(false);
                } else {
                    return Err(e.into());
                }
            }
        }
        let header = BlockHeader::from_bytes(&header_buf);
        self.block.read_from(&mut self.inner, header)?;
        // eprintln!("{:?}", header);

        Ok(true)
    }
}

fn write_file(ipath: &str, opath: &str) -> Result<()> {
    let handle = io::BufWriter::new(fs::File::create(opath)?);
    let mut writer = ColumnarBlockWriter::new(handle, 1024 * 1024);

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
    writer.flush()?;

    Ok(())
}

fn read_file(ipath: &str) -> Result<()> {
    let rhandle = fs::File::open(ipath).map(io::BufReader::new)?;
    let mut reader = Reader::new(rhandle);

    while reader.read_block()? {
        reader.block.decompress_columns()?;
        // break;
    }
    Ok(())
}

fn main() -> Result<()> {
    let ipath = "./data/some.fq.gz";
    let opath = "./data/some.cbq";

    eprintln!("Writing file {} - reading from {}", opath, ipath);
    write_file(ipath, opath)?;

    eprintln!("Reading file {}", opath);
    read_file(opath)?;

    Ok(())
}
