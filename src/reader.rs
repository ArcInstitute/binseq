//! Binary sequence reader module
//!
//! This module provides functionality for reading binary sequence files using memory mapping
//! for efficient access. It supports both sequential and parallel processing of records,
//! with configurable record layouts for different sequence types.

use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use bytemuck::cast_slice;
use memmap2::Mmap;

use crate::error::{ReadError, Result};
use crate::header::{BinseqHeader, SIZE_HEADER};
use crate::ParallelProcessor;

/// A reference to a binary sequence record in a memory-mapped file
///
/// This struct provides a view into a single record within a binary sequence file,
/// allowing access to the record's components (sequence data, flags, etc.) without
/// copying the data from the memory-mapped file.
///
/// The record's data is stored in a compact binary format where:
/// - The first u64 contains flags
/// - Subsequent u64s contain the primary sequence data
/// - If present, final u64s contain the extended sequence data
#[derive(Clone, Copy)]
pub struct RefRecord<'a> {
    /// The position (index) of this record in the file (0-based record index, not byte offset)
    id: usize,
    /// The underlying u64 buffer representing the record's binary data
    buffer: &'a [u64],
    /// The configuration that defines the layout and size of record components
    config: RecordConfig,
}
impl<'a> RefRecord<'a> {
    /// Creates a new record reference
    ///
    /// # Arguments
    ///
    /// * `id` - The record's position in the file (0-based record index, not byte offset)
    /// * `buffer` - The u64 slice containing the record's binary data
    /// * `config` - Configuration defining the record's layout
    ///
    /// # Panics
    ///
    /// Panics if the buffer length doesn't match the expected size from the config
    pub fn new(id: usize, buffer: &'a [u64], config: RecordConfig) -> Self {
        assert_eq!(buffer.len(), config.record_size_u64());
        Self { id, buffer, config }
    }
    /// Returns the record's position (index) in the file
    pub fn id(&self) -> usize {
        self.id
    }

    /// Returns the record's flag value
    ///
    /// The flag is stored in the first u64 of the record and can contain
    /// various metadata about the sequence.
    pub fn flag(&self) -> u64 {
        self.buffer[0]
    }

    /// Returns a reference to the primary sequence data buffer
    ///
    /// The returned slice contains the compressed nucleotide sequence,
    /// where each u64 stores up to 32 nucleotides.
    pub fn sbuf(&self) -> &[u64] {
        &self.buffer[1..1 + self.config.schunk]
    }

    /// Returns a reference to the extended sequence data buffer
    ///
    /// The returned slice contains the compressed extended data (e.g., quality scores),
    /// where each u64 stores 32 values.
    pub fn xbuf(&self) -> &[u64] {
        &self.buffer[1 + self.config.schunk..]
    }
    /// Decodes the primary sequence into a byte buffer
    ///
    /// This method decompresses the binary sequence data into standard nucleotide
    /// representation (A, C, G, T).
    ///
    /// # Arguments
    ///
    /// * `dbuf` - The buffer to store the decoded sequence
    ///
    /// # Errors
    ///
    /// Returns an error if the decoding process fails
    pub fn decode_s(&self, dbuf: &mut Vec<u8>) -> Result<()> {
        crate::nuc::decode(self.sbuf(), self.config.slen, dbuf)?;
        Ok(())
    }

    /// Decodes the extended sequence data into a byte buffer
    ///
    /// This method decompresses the binary extended data (e.g., quality scores)
    /// into its original representation.
    ///
    /// # Arguments
    ///
    /// * `dbuf` - The buffer to store the decoded data
    ///
    /// # Errors
    ///
    /// Returns an error if the decoding process fails
    pub fn decode_x(&self, dbuf: &mut Vec<u8>) -> Result<()> {
        crate::nuc::decode(self.xbuf(), self.config.xlen, dbuf)?;
        Ok(())
    }
    /// Returns whether this record contains extended sequence data
    ///
    /// A record is considered paired if it has a non-zero extended sequence length.
    pub fn paired(&self) -> bool {
        self.config.paired()
    }

    /// Returns the record's configuration
    ///
    /// The configuration defines the layout and size of the record's components.
    pub fn config(&self) -> RecordConfig {
        self.config
    }
}

/// Configuration for binary sequence record layout
///
/// This struct defines the size and layout of binary sequence records,
/// including both primary sequence data and optional extended data.
/// It handles the translation between sequence lengths in base pairs
/// and the number of u64 chunks needed to store the compressed data.
#[derive(Clone, Copy)]
pub struct RecordConfig {
    /// The primary sequence length in base pairs
    slen: usize,
    /// The extended sequence length in base pairs
    xlen: usize,
    /// The number of u64 chunks needed to store the primary sequence
    /// (each u64 stores 32 nucleotides)
    schunk: usize,
    /// The number of u64 chunks needed to store the extended sequence
    /// (each u64 stores 32 values)
    xchunk: usize,
}
impl RecordConfig {
    /// Creates a new record configuration
    ///
    /// This constructor initializes a configuration for a binary sequence record
    /// with specified primary and extended sequence lengths.
    ///
    /// # Arguments
    ///
    /// * `slen` - The length of primary sequences in the file
    /// * `xlen` - The length of secondary/extended sequences in the file
    ///
    /// # Returns
    ///
    /// A new `RecordConfig` instance with the specified sequence lengths
    pub fn new(slen: usize, xlen: usize) -> Self {
        Self {
            slen,
            xlen,
            schunk: slen.div_ceil(32),
            xchunk: xlen.div_ceil(32),
        }
    }

    /// Creates a new record configuration from a header
    ///
    /// This constructor initializes a configuration based on a header that contains
    /// the sequence lengths for primary and extended sequences.
    ///
    /// # Arguments
    ///
    /// * `header` - A reference to a `BinseqHeader` containing sequence lengths
    ///
    /// # Returns
    ///
    /// A new `RecordConfig` instance with the sequence lengths from the header
    pub fn from_header(header: &BinseqHeader) -> Self {
        Self::new(header.slen as usize, header.xlen as usize)
    }

    /// Returns whether this record contains extended sequence data
    ///
    /// A record is considered paired if it has a non-zero extended sequence length.
    pub fn paired(&self) -> bool {
        self.xlen > 0
    }

    /// Returns the primary sequence length in base pairs
    ///
    /// This method returns the length of the primary sequence in base pairs.
    pub fn slen(&self) -> usize {
        self.slen
    }

    /// Returns the extended sequence length in base pairs
    ///
    /// This method returns the length of the extended sequence in base pairs.
    pub fn xlen(&self) -> usize {
        self.xlen
    }

    /// Returns the number of u64 chunks needed to store the primary sequence
    ///
    /// This method returns the number of u64 chunks required to store the primary
    /// sequence, where each u64 stores 32 nucleotides.
    pub fn schunk(&self) -> usize {
        self.schunk
    }

    /// Returns the number of u64 chunks needed to store the extended sequence
    ///
    /// This method returns the number of u64 chunks required to store the extended
    /// sequence, where each u64 stores 32 values.
    pub fn xchunk(&self) -> usize {
        self.xchunk
    }

    /// Returns the full record size in bytes (u8):
    /// 8 * (schunk + xchunk + 1 (flag))
    pub fn record_size_bytes(&self) -> usize {
        8 * self.record_size_u64()
    }

    /// Returns the full record size in u64
    /// schunk + xchunk + 1 (flag)
    pub fn record_size_u64(&self) -> usize {
        self.schunk + self.xchunk + 1
    }
}

/// A memory-mapped reader for binary sequence files
///
/// This reader provides efficient access to binary sequence files by memory-mapping
/// them instead of performing traditional I/O operations. It supports both
/// sequential access to individual records and parallel processing of records
/// across multiple threads.
///
/// The reader ensures thread-safety through the use of `Arc` for sharing the
/// memory-mapped data between threads.
pub struct MmapReader {
    /// Memory mapped file contents, wrapped in Arc for thread-safe sharing
    mmap: Arc<Mmap>,

    /// Binary sequence file header containing format information
    header: BinseqHeader,

    /// Configuration defining the layout of records in the file
    config: RecordConfig,
}

impl MmapReader {
    /// Creates a new memory-mapped reader for a binary sequence file
    ///
    /// This method opens the file, memory-maps its contents, and validates
    /// the file structure to ensure it contains valid binary sequence data.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the binary sequence file
    ///
    /// # Returns
    ///
    /// * `Ok(MmapReader)` - A new reader if the file is valid
    /// * `Err(Error)` - If the file is invalid or cannot be opened
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// * The file cannot be opened
    /// * The file is not a regular file
    /// * The file header is invalid
    /// * The file size doesn't match the expected size based on the header
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        // Verify input file is a file before attempting to map
        let file = File::open(path)?;
        if !file.metadata()?.is_file() {
            return Err(ReadError::IncompatibleFile.into());
        }

        // Safety: the file is open and won't be modified while mapped
        let mmap = unsafe { Mmap::map(&file)? };

        // Read header from mapped memory
        let header = BinseqHeader::from_buffer(&mmap)?;

        // Record configuraration
        let config = RecordConfig::from_header(&header);

        // Immediately validate the size of the file against the expected byte size of records
        if (mmap.len() - SIZE_HEADER) % config.record_size_bytes() != 0 {
            return Err(ReadError::FileTruncation(mmap.len()).into());
        }

        Ok(Self {
            mmap: Arc::new(mmap),
            header,
            config,
        })
    }

    /// Returns the total number of records in the file
    ///
    /// This is calculated by subtracting the header size from the total file size
    /// and dividing by the size of each record.
    pub fn num_records(&self) -> usize {
        (self.mmap.len() - SIZE_HEADER) / self.config.record_size_bytes()
    }

    /// Returns a copy of the binary sequence file header
    ///
    /// The header contains format information and sequence length specifications.
    pub fn header(&self) -> BinseqHeader {
        self.header
    }

    /// Returns a reference to a specific record
    ///
    /// # Arguments
    ///
    /// * `idx` - The index of the record to retrieve (0-based)
    ///
    /// # Returns
    ///
    /// * `Ok(RefRecord)` - A reference to the requested record
    /// * `Err(Error)` - If the index is out of bounds
    ///
    /// # Errors
    ///
    /// Returns an error if the requested index is beyond the number of records in the file
    pub fn get(&self, idx: usize) -> Result<RefRecord> {
        if idx > self.num_records() {
            return Err(ReadError::OutOfRange(idx, self.num_records()).into());
        }
        let lbound = SIZE_HEADER + (idx * self.config.record_size_bytes());
        let rbound = lbound + self.config.record_size_bytes();
        let bytes = &self.mmap[lbound..rbound];
        let buffer = cast_slice(bytes);
        Ok(RefRecord::new(idx, buffer, self.config))
    }
}

/// Default batch size for parallel processing
///
/// This constant defines how many records each thread processes at a time
/// during parallel processing operations.
pub const BATCH_SIZE: usize = 1024;

/// Parallel processing implementation for memory-mapped readers
impl MmapReader {
    /// Processes all records in parallel using multiple threads
    ///
    /// This method distributes the records across the specified number of threads
    /// and processes them using the provided processor. Each thread receives its
    /// own clone of the processor and processes a contiguous chunk of records.
    ///
    /// # Arguments
    ///
    /// * `processor` - The processor to use for handling records
    /// * `num_threads` - The number of threads to use for processing
    ///
    /// # Type Parameters
    ///
    /// * `P` - A type that implements `ParallelProcessor` and can be cloned
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If all records were processed successfully
    /// * `Err(Error)` - If an error occurred during processing
    pub fn process_parallel<P: ParallelProcessor + Clone + 'static>(
        self,
        processor: P,
        num_threads: usize,
    ) -> Result<()> {
        // Calculate number of records for each thread
        let num_records = self.num_records();
        let records_per_thread = num_records.div_ceil(num_threads);

        // Arc self
        let reader = Arc::new(self);

        // Build thread handles
        let mut handles = Vec::new();
        for tid in 0..num_threads {
            let mut processor = processor.clone();
            let reader = reader.clone();
            processor.set_tid(tid);

            let handle = std::thread::spawn(move || -> Result<()> {
                let start_idx = tid * records_per_thread;
                let end_idx = (start_idx + records_per_thread).min(num_records);

                for (batch_idx, idx) in (start_idx..end_idx).enumerate() {
                    let record = reader.get(idx)?;
                    processor.process_record(record)?;

                    if batch_idx % BATCH_SIZE == 0 {
                        processor.on_batch_complete()?;
                    }
                }
                processor.on_batch_complete()?;

                Ok(())
            });

            handles.push(handle);
        }

        for handle in handles {
            handle
                .join()
                .expect("Error joining handle (1)")
                .expect("Error joining handle (2)");
        }

        Ok(())
    }
}
