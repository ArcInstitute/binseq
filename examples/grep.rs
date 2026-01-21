use std::sync::Arc;

use anyhow::Result;
use binseq::prelude::*;
use clap::Parser;
use memchr::memmem::Finder;
use parking_lot::Mutex;

#[derive(Clone)]
pub struct GrepCounter {
    // (thread) local variables
    local_count: usize,

    // search pattern (using memchr::memmem::Finder for fast searching)
    pattern: Finder<'static>,

    // global variables
    count: Arc<Mutex<usize>>,
}
impl GrepCounter {
    #[must_use]
    pub fn new(pattern: &[u8]) -> Self {
        Self {
            pattern: Finder::new(pattern).into_owned(),
            local_count: 0,
            count: Arc::new(Mutex::new(0)),
        }
    }

    fn match_sequence(&self, seq: &[u8]) -> bool {
        self.pattern.find(seq).is_some()
    }

    fn pprint(&self) {
        println!("Matching records: {}", self.count.lock());
    }
}
impl ParallelProcessor for GrepCounter {
    fn process_record<R: binseq::BinseqRecord>(&mut self, record: R) -> binseq::Result<()> {
        if self.match_sequence(&record.sseq()) || self.match_sequence(&record.xseq()) {
            self.local_count += 1;
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        *self.count.lock() += self.local_count;
        self.local_count = 0;
        Ok(())
    }
}

#[derive(Parser)]
struct Args {
    /// Input BINSEQ path to grep
    #[clap(required = true)]
    input: String,

    /// Pattern to search for (either sseq or xseq)
    #[clap(required = true)]
    pattern: String,

    /// Threads to use [0: auto]
    #[clap(short = 'T', long, default_value_t = 0)]
    threads: usize,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let reader = BinseqReader::new(&args.input)?;
    let counter = GrepCounter::new(args.pattern.as_bytes());
    reader.process_parallel(counter.clone(), args.threads)?;
    counter.pprint();
    Ok(())
}
