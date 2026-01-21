use std::{
    fs::File,
    io::BufWriter,
    sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    },
};

use anyhow::{Result, bail};
use binseq::{
    bq::{self, BinseqHeaderBuilder},
    prelude::*,
};
use nucgen::Sequence;

#[derive(Clone, Default)]
pub struct MyProcessor {
    local_counter: usize,
    counter: Arc<AtomicUsize>,
}
impl MyProcessor {
    #[must_use]
    pub fn counter(&self) -> usize {
        self.counter.load(Ordering::Relaxed)
    }
}
impl ParallelProcessor for MyProcessor {
    fn process_record<R: BinseqRecord>(&mut self, _record: R) -> binseq::Result<()> {
        self.local_counter += 1;
        Ok(())
    }
    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        self.counter
            .fetch_add(self.local_counter, Ordering::Relaxed);
        self.local_counter = 0;
        Ok(())
    }
}

fn mmap_processing(binseq_path: &str, n_threads: usize) -> Result<()> {
    let reader = BinseqReader::new(binseq_path)?;
    let proc = MyProcessor::default();
    reader.process_parallel(proc.clone(), n_threads)?;
    Ok(())
}

pub fn main() -> Result<()> {
    let binseq_path_single = "./data/test.bq";
    let binseq_path_paired = "./data/test_paired.bq";
    let r1_size = 150;
    let r2_size = 300;
    let num_seq = 1_000_000;

    time_it(
        || {
            write_single(binseq_path_single, num_seq, r1_size)?;
            Ok(())
        },
        "write_single",
    );

    time_it(
        || {
            write_paired(binseq_path_paired, num_seq, r1_size, r2_size)?;
            Ok(())
        },
        "write_paired",
    );

    for n_threads in 1..=16 {
        if n_threads > 1 && n_threads % 2 != 0 {
            continue;
        }
        time_it(
            || {
                mmap_processing(binseq_path_single, n_threads)?;
                Ok(())
            },
            &format!("single - mmap_parallel_processing ({n_threads})"),
        );
    }
    for n_threads in 1..=16 {
        if n_threads > 1 && n_threads % 2 != 0 {
            continue;
        }
        time_it(
            || {
                mmap_processing(binseq_path_paired, n_threads)?;
                Ok(())
            },
            &format!("paired - mmap_parallel_processing ({n_threads})"),
        );
    }

    Ok(())
}

fn time_it<F>(f: F, name: &str)
where
    F: Fn() -> Result<()>,
{
    let now = std::time::Instant::now();
    f().unwrap();
    let elapsed = now.elapsed();
    eprintln!("Elapsed time ({name}): {elapsed:?}");
}

fn write_single(binseq_path: &str, num_seq: usize, seq_size: usize) -> Result<()> {
    // Open the output file
    let header = BinseqHeaderBuilder::new().slen(seq_size as u32).build()?;
    let out_handle = File::create(binseq_path).map(BufWriter::new)?;
    let mut writer = bq::BinseqWriterBuilder::default()
        .header(header)
        .build(out_handle)?;

    // Write the binary sequence
    let mut sequence = Sequence::new();
    let mut rng = rand::rng();
    for _ in 0..num_seq {
        sequence.fill_buffer(&mut rng, seq_size);
        if !writer.write_record(Some(0), sequence.bytes())? {
            bail!("Error writing nucleotides")
        }
    }
    writer.flush()?;
    eprintln!("Finished writing {num_seq} records to path: {binseq_path}");
    Ok(())
}

fn write_paired(binseq_path: &str, num_seq: usize, r1_size: usize, r2_size: usize) -> Result<()> {
    // Open the output file
    let header = bq::BinseqHeaderBuilder::new()
        .slen(r1_size as u32)
        .xlen(r2_size as u32)
        .build()?;
    let out_handle = File::create(binseq_path).map(BufWriter::new)?;
    let mut writer = bq::BinseqWriterBuilder::default()
        .header(header)
        .build(out_handle)?;

    // Write the binary sequence
    let mut r1 = Sequence::new();
    let mut r2 = Sequence::new();
    let mut rng = rand::rng();
    for _ in 0..num_seq {
        r1.fill_buffer(&mut rng, r1_size);
        r2.fill_buffer(&mut rng, r2_size);

        if !writer.write_paired_record(Some(0), r1.bytes(), r2.bytes())? {
            bail!("Error writing nucleotides")
        }
    }
    writer.flush()?;
    eprintln!("Finished writing {num_seq} records to path: {binseq_path}");
    Ok(())
}
