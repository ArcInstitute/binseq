use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use clap::Parser;

use binseq::{BinseqReader, BinseqRecord, ParallelProcessor, ParallelReader, Result};

#[derive(Clone)]
struct RangeProcessor {
    counter: Arc<AtomicUsize>,
    tid: Option<usize>,
    range_start: usize,
    range_end: usize,
}

impl RangeProcessor {
    fn new(range_start: usize, range_end: usize) -> Self {
        Self {
            counter: Arc::new(AtomicUsize::new(0)),
            tid: None,
            range_start,
            range_end,
        }
    }

    fn count(&self) -> usize {
        self.counter.load(Ordering::Relaxed)
    }
}

impl ParallelProcessor for RangeProcessor {
    fn process_record<R: BinseqRecord>(&mut self, record: R) -> Result<()> {
        let count = self.counter.fetch_add(1, Ordering::Relaxed);

        // Print progress every 10,000 records
        if count.is_multiple_of(10_000)
            && let Some(tid) = self.tid
        {
            println!(
                "Thread {}: Processed {} records (Range: {}-{}, Index: {}, Len: {})",
                tid,
                count + 1,
                self.range_start,
                self.range_end,
                record.index(),
                record.sseq().len(),
            );
        }

        Ok(())
    }

    fn set_tid(&mut self, tid: usize) {
        self.tid = Some(tid);
    }

    fn get_tid(&self) -> Option<usize> {
        self.tid
    }

    fn on_batch_complete(&mut self) -> Result<()> {
        if let Some(tid) = self.tid {
            println!("Thread {tid} completed batch processing");
        }
        Ok(())
    }
}

#[derive(Parser)]
struct Cli {
    path: PathBuf,

    #[clap(short = 'T', long, default_value_t = 0)]
    threads: usize,

    /// Optional start index for processing range
    #[clap(long, default_value_t = 0)]
    start: usize,

    /// Optional end index for processing range
    #[clap(long, default_value_t = 10_000)]
    end: usize,
}

fn main() -> Result<()> {
    let args = Cli::parse();

    // Create reader to get total record count
    let reader = BinseqReader::new(&args.path)?;
    let total_records = reader.num_records()?;

    println!("File: {}", args.path.display());
    println!("Total records in file: {total_records}");

    // Validate range
    if args.start >= total_records {
        eprintln!(
            "Error: Start index {} is >= total records {total_records}",
            args.start
        );
        std::process::exit(1);
    }
    if args.end > total_records {
        eprintln!(
            "Warning: End index {} is > total records {total_records}, clamping to {total_records}",
            args.end
        );
    }
    let end = args.end.min(total_records);

    if args.start >= end {
        eprintln!(
            "Error: Start index {} must be < end index {end}",
            args.start
        );
        std::process::exit(1);
    }

    println!(
        "Processing range: {} to {} ({} records)",
        args.start,
        end,
        end - args.start
    );
    println!("Using {} threads", args.threads);
    println!();

    // Demonstrate processing the full file
    println!("=== Processing full file ===");
    let reader_full = BinseqReader::new(&args.path)?;
    let processor_full = RangeProcessor::new(0, total_records);
    let start_time = std::time::Instant::now();

    reader_full.process_parallel(processor_full.clone(), args.threads)?;

    let elapsed_full = start_time.elapsed();
    println!("Full file processing completed!");
    println!("Records processed: {}", processor_full.count());
    println!("Time taken: {elapsed_full:.2?}");
    println!();

    // Demonstrate processing a specific range
    println!("=== Processing specific range ===");
    let reader_range = BinseqReader::new(args.path)?;
    let processor_range = RangeProcessor::new(args.start, end);
    let start_time = std::time::Instant::now();

    reader_range.process_parallel_range(processor_range.clone(), args.threads, args.start..end)?;

    let elapsed_range = start_time.elapsed();
    println!("Range processing completed!");
    println!("Records processed: {}", processor_range.count());
    println!("Expected records: {}", end - args.start);
    println!("Time taken: {elapsed_range:.2?}");

    // Compare performance
    if processor_range.count() > 0 && processor_full.count() > 0 {
        let full_rate = processor_full.count() as f64 / elapsed_full.as_secs_f64();
        let range_rate = processor_range.count() as f64 / elapsed_range.as_secs_f64();
        println!();
        println!("=== Performance Comparison ===");
        println!("Full file rate: {full_rate:.0} records/sec");
        println!("Range rate: {range_rate:.0} records/sec");

        if range_rate > full_rate {
            println!(
                "Range processing was {:.1}x faster per record",
                range_rate / full_rate
            );
        } else {
            println!(
                "Full file processing was {:.1}x faster per record",
                full_rate / range_rate
            );
        }
    }

    Ok(())
}
