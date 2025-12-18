use std::{fs, io, sync::Arc};

use anyhow::Result;
use binseq::{BinseqRecord, ParallelProcessor, ParallelReader};
use cbq::{BlockRange, ColumnarBlockWriter, FileHeader, MmapReader, SequencingRecord};
use paraseq::{Record, fastx};
use parking_lot::Mutex;

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
    let mut reader = cbq::Reader::new(rhandle)?;
    let mut writer = io::BufWriter::new(io::stdout());

    let mut total_records = 0;
    while let Some(header) = reader.read_block()? {
        reader.block.decompress_columns()?;
        let range = BlockRange::new(0, total_records + header.num_records);
        for record in reader.block.iter_records(range) {
            write_fastq(&mut writer, record.sheader(), record.sseq(), record.squal())?;
        }
        total_records += header.num_records;
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
