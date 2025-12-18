use std::{
    io::{self, BufWriter},
    sync::Arc,
};

use anyhow::Result;
use cbq::{ColumnarBlockWriter, FileHeader, SequencingRecordBuilder};
use paraseq::{
    Record, fastx,
    prelude::{ParallelProcessor, ParallelReader},
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

pub fn main() -> Result<()> {
    let path = std::env::args()
        .nth(1)
        .expect("Usage: writer <input_path>.fastx");
    let reader = fastx::Reader::from_path(path)?;
    let handle = Box::new(BufWriter::new(io::stdout()));
    let header = FileHeader::default();
    let mut proc = ParallelWriter::new(handle, header)?;
    reader.process_parallel(&mut proc, 0)?;
    proc.finish()?;
    Ok(())
}
