//! Conditional decoding of CBQ columns.
//!
//! CBQ stores each record attribute in its own ZSTD-compressed column
//! (sequence, flags, headers, quality scores). A processor that only touches
//! *some* of those attributes pays for decompressing all of them by default.
//!
//! [`DecompressionOptions`] lets you tell the reader up-front which columns it
//! can leave compressed. Skipped columns are stepped over in the block payload,
//! so the columns you *do* want still decode correctly.
//!
//! IMPORTANT: skipping a column means the buffer backing it is empty for the
//! whole block. Calling the accessor for a skipped column
//! (`sseq`/`xseq`, `sheader`/`xheader`, `squal`/`xqual`) will panic on an
//! out-of-bounds slice; `flag()` returns `None`. Only skip what your processor
//! genuinely never reads.
//!
//! Note: `slen()`/`xlen()` come from the sequence-length column, which is
//! always decoded, so record lengths remain available even with every optional
//! column skipped.
//!
//! Usage:
//!
//! ```text
//! cargo run --release --example conditional_decoding -- [path.cbq] [num_threads]
//! ```

use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use binseq::Result;
use binseq::cbq::{DecompressionOptions, MmapReader};
use binseq::prelude::*;

use clap::Parser;

/// A processor that only ever looks at nucleotides.
///
/// It never calls `sheader()`, `squal()`, or `flag()`, so the header, quality,
/// and flag columns can all stay compressed.
#[derive(Clone, Default)]
struct GcCounter {
    gc: Arc<AtomicUsize>,
    total: Arc<AtomicUsize>,

    // thread-local accumulators (flushed on batch completion)
    local_gc: usize,
    local_total: usize,
}
impl GcCounter {
    fn gc_fraction(&self) -> f64 {
        let total = self.total.load(Ordering::Relaxed);
        if total == 0 {
            0.0
        } else {
            self.gc.load(Ordering::Relaxed) as f64 / total as f64
        }
    }
}
impl ParallelProcessor for GcCounter {
    fn process_record<R: BinseqRecord>(&mut self, record: R) -> Result<()> {
        for seq in [record.sseq(), record.xseq()] {
            self.local_gc += seq
                .iter()
                .filter(|b| matches!(b, b'G' | b'C' | b'g' | b'c'))
                .count();
            self.local_total += seq.len();
        }
        Ok(())
    }

    fn on_batch_complete(&mut self) -> Result<()> {
        self.gc.fetch_add(self.local_gc, Ordering::Relaxed);
        self.total.fetch_add(self.local_total, Ordering::Relaxed);
        self.local_gc = 0;
        self.local_total = 0;
        Ok(())
    }
}

/// A processor that only needs record *lengths* - no columns at all.
///
/// Sequence lengths live in a column that is always decoded, so this one can
/// skip every optional column in the block.
#[derive(Clone, Default)]
struct LengthHistogram {
    nucleotides: Arc<AtomicUsize>,
    records: Arc<AtomicUsize>,
    local_nuc: usize,
    local_records: usize,
}
impl ParallelProcessor for LengthHistogram {
    fn process_record<R: BinseqRecord>(&mut self, record: R) -> Result<()> {
        self.local_nuc += (record.slen() + record.xlen()) as usize;
        self.local_records += 1;
        Ok(())
    }

    fn on_batch_complete(&mut self) -> Result<()> {
        self.nucleotides
            .fetch_add(self.local_nuc, Ordering::Relaxed);
        self.records
            .fetch_add(self.local_records, Ordering::Relaxed);
        self.local_nuc = 0;
        self.local_records = 0;
        Ok(())
    }
}

#[derive(Parser)]
struct Cli {
    path: PathBuf,

    #[clap(short = 'T', long, default_value_t = 0)]
    threads: usize,

    #[clap(short = 's', long)]
    only_seq: bool,

    #[clap(short = 'd', long, conflicts_with = "only_seq")]
    no_dctx: bool,
}

fn main() -> Result<()> {
    let args = Cli::parse();

    let mut reader = MmapReader::new(&args.path)?;

    let start = Instant::now();
    if args.no_dctx {
        let proc = LengthHistogram::default();

        let mut opt = DecompressionOptions::default();
        opt.set_skip_flags(true);
        opt.set_skip_header(true);
        opt.set_skip_sequence(true);
        opt.set_skip_quality(true);

        reader.process_parallel(proc.clone(), args.threads)?;
        let length_elapsed = start.elapsed();
        println!(
            "{} records / {} nucleotides (no columns decoded) in {length_elapsed:?}",
            proc.records.load(Ordering::Relaxed),
            proc.nucleotides.load(Ordering::Relaxed),
        );
    } else {
        let proc = GcCounter::default();

        if args.only_seq {
            let mut opt = DecompressionOptions::default();
            opt.set_skip_flags(true);
            opt.set_skip_header(true);
            opt.set_skip_quality(true);
            reader.set_decompression_options(opt);
        }

        reader.process_parallel(proc.clone(), args.threads)?;
        let elapsed = start.elapsed();

        println!("GC fraction: {}", proc.gc_fraction());
        println!("Elapsed: {elapsed:?}");
    }
    Ok(())
}
