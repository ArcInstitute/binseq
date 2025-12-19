use std::fs::File;
use std::io::{BufWriter, Write, stdout};
use std::sync::Arc;

use anyhow::Result;
use binseq::prelude::*;

use parking_lot::Mutex;

/// A struct for decoding BINSEQ data back to FASTQ format.
#[derive(Clone)]
pub struct Decoder {
    /// Reusable context
    ctx: Ctx,

    /// local output buffer
    local_writer: Vec<u8>,

    /// global output buffer
    global_writer: Arc<Mutex<Box<dyn Write + Send>>>,

    /// Local count of records
    local_count: usize,

    /// global count of records
    global_count: Arc<Mutex<usize>>,
}

impl Decoder {
    #[must_use]
    pub fn new(writer: Box<dyn Write + Send>) -> Self {
        let global_writer = Arc::new(Mutex::new(writer));
        Decoder {
            local_writer: Vec::new(),
            ctx: Ctx::default(),
            local_count: 0,
            global_writer,
            global_count: Arc::new(Mutex::new(0)),
        }
    }

    #[must_use]
    pub fn num_records(&self) -> usize {
        *self.global_count.lock()
    }
}
impl ParallelProcessor for Decoder {
    fn process_record<R: BinseqRecord>(&mut self, record: R) -> binseq::Result<()> {
        self.ctx.fill(&record)?;
        write_fastq_parts(
            &mut self.local_writer,
            self.ctx.sheader(),
            self.ctx.sbuf(),
            self.ctx.squal(),
        )?;

        // write extended fastq to local buffer
        if record.is_paired() {
            write_fastq_parts(
                &mut self.local_writer,
                self.ctx.xheader(),
                &self.ctx.xbuf(),
                self.ctx.xqual(),
            )?;
        }

        self.local_count += 1;
        Ok(())
    }

    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        // Lock the mutex to write to the global buffer
        {
            let mut lock = self.global_writer.lock();
            lock.write_all(&self.local_writer)?;
            lock.flush()?;
        }
        // Lock the mutex to update the number of records
        {
            let mut global_count = self.global_count.lock();
            *global_count += self.local_count;
        }

        // Clear the local buffer and reset the local record count
        self.local_writer.clear();
        self.local_count = 0;
        Ok(())
    }
}

#[allow(clippy::missing_errors_doc)]
pub fn write_fastq_parts<W: Write>(
    writer: &mut W,
    index: &[u8],
    sequence: &[u8],
    quality: &[u8],
) -> Result<(), std::io::Error> {
    writer.write_all(b"@seq.")?;
    writer.write_all(index)?;
    writer.write_all(b"\n")?;
    writer.write_all(sequence)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(quality)?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn match_output(path: Option<&str>) -> Result<Box<dyn Write + Send>> {
    if let Some(path) = path {
        let writer = File::create(path).map(BufWriter::new)?;
        Ok(Box::new(writer))
    } else {
        let stdout = stdout();
        Ok(Box::new(BufWriter::new(stdout)))
    }
}

fn main() -> Result<()> {
    let file = std::env::args()
        .nth(1)
        .unwrap_or("./data/subset.bq".to_string());
    let n_threads = std::env::args().nth(2).unwrap_or("1".to_string()).parse()?;

    let reader = BinseqReader::new(&file)?;
    let writer = match_output(None)?;
    let proc = Decoder::new(writer);

    reader.process_parallel(proc.clone(), n_threads)?;
    eprintln!("Read {} records", proc.num_records());

    Ok(())
}
