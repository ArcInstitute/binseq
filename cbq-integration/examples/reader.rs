use std::{io, sync::Arc};

use anyhow::Result;
use binseq::prelude::*;
use cbq::MmapReader;
use clap::Parser;
use parking_lot::Mutex;

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

#[derive(Parser)]
struct Args {
    #[clap(required = true)]
    input: String,

    #[clap(short, long)]
    output: Option<String>,

    #[clap(short = 'T', long, default_value_t = 0)]
    threads: usize,
}

fn match_output(path: Option<&str>) -> Result<BoxedWriter> {
    match path {
        Some(path) => {
            let handle = std::fs::File::create(path).map(io::BufWriter::new)?;
            Ok(Box::new(handle))
        }
        None => Ok(Box::new(io::stdout())),
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    let reader = MmapReader::new(&args.input)?;
    let handle = match_output(args.output.as_deref())?;
    let proc = Processor::new(handle);
    reader.process_parallel(proc.clone(), 0)?;
    eprintln!("Number of records: {}", proc.n_records());
    Ok(())
}
