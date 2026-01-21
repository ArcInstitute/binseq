use std::ops::Range;
use std::path::Path;

use crate::{BinseqRecord, Result, bq, cbq, error::ExtensionError, vbq};

/// An enum abstraction for BINSEQ readers that can process records in parallel
///
/// This is a convenience enum that can be used for general workflows where the
/// distinction between BQ and VBQ readers is not important.
///
/// For more specialized workflows see [`bq::MmapReader`] and [`vbq::MmapReader`].
pub enum BinseqReader {
    Bq(bq::MmapReader),
    Vbq(vbq::MmapReader),
    Cbq(cbq::MmapReader),
}
impl BinseqReader {
    pub fn new(path: &str) -> Result<Self> {
        let pathbuf = Path::new(path);
        match pathbuf.extension() {
            Some(ext) => match ext.to_str() {
                Some("bq") => Ok(Self::Bq(bq::MmapReader::new(path)?)),
                Some("vbq") => Ok(Self::Vbq(vbq::MmapReader::new(path)?)),
                Some("cbq") => Ok(Self::Cbq(cbq::MmapReader::new(path)?)),
                _ => Err(ExtensionError::UnsupportedExtension(path.to_string()).into()),
            },
            None => Err(ExtensionError::UnsupportedExtension(path.to_string()).into()),
        }
    }

    /// Set whether to decode sequences at once in each block
    ///
    /// Note: This setting applies to VBQ readers only.
    pub fn set_decode_block(&mut self, decode_block: bool) {
        match self {
            Self::Bq(_) | Self::Cbq(_) => {
                // no-op
            }
            Self::Vbq(reader) => reader.set_decode_block(decode_block),
        }
    }

    pub fn set_default_quality_score(&mut self, score: u8) {
        match self {
            Self::Bq(reader) => reader.set_default_quality_score(score),
            Self::Vbq(reader) => reader.set_default_quality_score(score),
            Self::Cbq(reader) => reader.set_default_quality_score(score),
        }
    }

    #[must_use]
    pub fn is_paired(&self) -> bool {
        match self {
            Self::Bq(reader) => reader.is_paired(),
            Self::Vbq(reader) => reader.is_paired(),
            Self::Cbq(reader) => reader.is_paired(),
        }
    }

    pub fn num_records(&self) -> Result<usize> {
        match self {
            Self::Bq(reader) => Ok(reader.num_records()),
            Self::Vbq(reader) => reader.num_records(),
            Self::Cbq(reader) => Ok(reader.num_records()),
        }
    }

    /// Process records in parallel within a specified range
    ///
    /// This method allows parallel processing of a subset of records within the file,
    /// defined by a start and end index. The range is distributed across the specified
    /// number of threads.
    ///
    /// # Arguments
    ///
    /// * `processor` - The processor to use for each record
    /// * `num_threads` - The number of threads to spawn
    /// * `start` - The starting record index (inclusive)
    /// * `end` - The ending record index (exclusive)
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If all records were processed successfully
    /// * `Err(Error)` - If an error occurred during processing
    pub fn process_parallel_range<P: ParallelProcessor + Clone + 'static>(
        self,
        processor: P,
        num_threads: usize,
        range: Range<usize>,
    ) -> Result<()> {
        match self {
            Self::Bq(reader) => reader.process_parallel_range(processor, num_threads, range),
            Self::Vbq(reader) => reader.process_parallel_range(processor, num_threads, range),
            Self::Cbq(reader) => reader.process_parallel_range(processor, num_threads, range),
        }
    }
}
impl ParallelReader for BinseqReader {
    fn process_parallel<P: ParallelProcessor + Clone + 'static>(
        self,
        processor: P,
        num_threads: usize,
    ) -> Result<()> {
        let num_records = self.num_records()?;
        self.process_parallel_range(processor, num_threads, 0..num_records)
    }

    fn process_parallel_range<P: ParallelProcessor + Clone + 'static>(
        self,
        processor: P,
        num_threads: usize,
        range: Range<usize>,
    ) -> Result<()> {
        match self {
            Self::Bq(reader) => reader.process_parallel_range(processor, num_threads, range),
            Self::Vbq(reader) => reader.process_parallel_range(processor, num_threads, range),
            Self::Cbq(reader) => reader.process_parallel_range(processor, num_threads, range),
        }
    }
}

/// Trait for BINSEQ readers that can process records in parallel
///
/// This is implemented by the **reader** not by the **processor**.
/// For the **processor**, see the [`ParallelProcessor`] trait.
pub trait ParallelReader {
    fn process_parallel<P: ParallelProcessor + Clone + 'static>(
        self,
        processor: P,
        num_threads: usize,
    ) -> Result<()>;

    /// Process records in parallel within a specified range
    ///
    /// This method allows parallel processing of a subset of records within the file,
    /// defined by a start and end index. The range is distributed across the specified
    /// number of threads.
    ///
    /// # Arguments
    ///
    /// * `processor` - The processor to use for each record
    /// * `num_threads` - The number of threads to spawn
    /// * `range` - The range of record indices to process
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If all records were processed successfully
    /// * `Err(Error)` - If an error occurred during processing
    fn process_parallel_range<P: ParallelProcessor + Clone + 'static>(
        self,
        processor: P,
        num_threads: usize,
        range: Range<usize>,
    ) -> Result<()>;
}

/// Trait for types that can process records in parallel.
///
/// This is implemented by the **processor** not by the **reader**.
/// For the **reader**, see the [`ParallelReader`] trait.
pub trait ParallelProcessor: Send + Clone {
    /// Process a single record
    fn process_record<R: BinseqRecord>(&mut self, record: R) -> Result<()>;

    /// Called when a thread finishes processing its batch
    /// Default implementation does nothing
    #[allow(unused_variables)]
    fn on_batch_complete(&mut self) -> Result<()> {
        Ok(())
    }

    /// Set the thread ID for this processor
    ///
    /// Each thread should call this method with its own unique ID.
    #[allow(unused_variables)]
    fn set_tid(&mut self, _tid: usize) {
        // Default implementation does nothing
    }

    /// Get the thread ID for this processor
    fn get_tid(&self) -> Option<usize> {
        None
    }
}
