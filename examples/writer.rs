use std::{fs, io, sync::Arc};

use anyhow::Result;
use cbq::{ColumnarBlockWriter, FileHeader, SequencingRecordBuilder};
use clap::Parser;
use paraseq::{
    Record, fastx,
    prelude::{PairedParallelProcessor, ParallelProcessor, ParallelReader},
};
use parking_lot::Mutex;

type BoxedWriter = Box<dyn io::Write + Send>;

#[derive(Clone)]
pub struct ParallelWriter {
    local_writer: cbq::ColumnarBlockWriter<Vec<u8>>,
    local_num_records: usize,
    writer: Arc<Mutex<cbq::ColumnarBlockWriter<BoxedWriter>>>,
    num_records: Arc<Mutex<usize>>,
}
impl ParallelWriter {
    pub fn new(handle: BoxedWriter, header: FileHeader) -> Result<Self> {
        let local_writer = ColumnarBlockWriter::new_headless(Vec::default(), header);
        let writer = ColumnarBlockWriter::new(handle, header)?;
        Ok(Self {
            local_writer,
            local_num_records: 0,
            writer: Arc::new(Mutex::new(writer)),
            num_records: Arc::new(Mutex::new(0)),
        })
    }
    pub fn finish(&mut self) -> Result<()> {
        let mut writer = self.writer.lock();
        writer.finish()?;
        Ok(())
    }
}
impl<R: Record> ParallelProcessor<R> for ParallelWriter {
    fn process_record(&mut self, record: R) -> paraseq::Result<()> {
        let seq = &record.seq();
        let seq_record = SequencingRecordBuilder::default()
            .s_seq(seq)
            .opt_s_qual(record.qual())
            .s_header(record.id())
            .build()?;
        self.local_writer.push(seq_record)?;
        self.local_num_records += 1;
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        {
            let mut writer = self.writer.lock();
            writer.ingest(&mut self.local_writer)?;
        }

        {
            *self.num_records.lock() += self.local_num_records;
            self.local_num_records = 0;
        }

        Ok(())
    }
}
impl<R: Record> PairedParallelProcessor<R> for ParallelWriter {
    fn process_record_pair(&mut self, r1: R, r2: R) -> paraseq::Result<()> {
        let s_seq = &r1.seq();
        let x_seq = &r2.seq();
        let seq_record = SequencingRecordBuilder::default()
            .s_seq(s_seq)
            .opt_s_qual(r1.qual())
            .s_header(r1.id())
            .x_seq(x_seq)
            .opt_x_qual(r2.qual())
            .x_header(r2.id())
            .build()?;
        self.local_writer.push(seq_record)?;
        self.local_num_records += 1;
        Ok(())
    }
    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        {
            let mut writer = self.writer.lock();
            writer.ingest(&mut self.local_writer)?;
        }

        {
            *self.num_records.lock() += self.local_num_records;
            self.local_num_records = 0;
        }

        Ok(())
    }
}

#[derive(Parser)]
struct Args {
    /// Input FASTX (single or paired)
    #[clap(required = true, num_args=1..2)]
    input: Vec<String>,

    /// Output cbq file
    #[clap(short, long)]
    output: String,

    #[clap(short = 'T', long, default_value_t = 0)]
    threads: usize,
}

pub fn main() -> Result<()> {
    let args = Args::parse();

    let handle = Box::new(fs::File::create(&args.output).map(io::BufWriter::new)?);
    let mut header = FileHeader::default();
    if args.input.len() == 2 {
        header.set_paired();
    }
    let mut proc = ParallelWriter::new(handle, header)?;

    if args.input.len() == 2 {
        eprintln!("Processing paired-end FASTX files");
        let r1 = fastx::Reader::from_path(&args.input[0])?;
        let r2 = fastx::Reader::from_path(&args.input[1])?;
        r1.process_parallel_paired(r2, &mut proc, args.threads)?;
    } else {
        eprintln!("Processing single-end FASTX file");
        let reader = fastx::Reader::from_path(&args.input[0])?;
        reader.process_parallel(&mut proc, args.threads)?;
    }
    proc.finish()?;
    Ok(())
}
