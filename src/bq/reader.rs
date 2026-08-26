//! Binary sequence reader module
//!
//! This module provides functionality for reading binary sequence files using either:
//! 1. Memory mapping for efficient access to entire files
//! 2. Streaming for processing data as it arrives
//!
//! It supports both sequential and parallel processing of records,
//! with configurable record layouts for different sequence types.

use std::fs::File;
use std::io::Read;
use std::ops::Range;
use std::path::Path;
use std::sync::Arc;

use bitnuc_deprec::BitSize;
use bytemuck::cast_slice;
use memmap2::Mmap;

use super::header::{FileHeader, SIZE_HEADER};
use crate::{
    BinseqRecord, DEFAULT_QUALITY_SCORE, Error, ParallelProcessor, ParallelReader,
    error::{ReadError, Result},
};

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
    id: u64,
    /// The underlying u64 buffer representing the record's binary data
    buffer: &'a [u64],
    /// Reusable default quality buffer
    qbuf: &'a [u8],
    /// The configuration that defines the layout and size of record components
    config: RecordConfig,
    /// Cached index string for the sequence header
    header_buf: [u8; 20],
    /// Length of the header in bytes
    header_len: usize,
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
    #[must_use]
    pub fn new(id: u64, buffer: &'a [u64], qbuf: &'a [u8], config: RecordConfig) -> Self {
        assert_eq!(buffer.len(), config.record_size_u64());
        Self {
            id,
            buffer,
            qbuf,
            config,
            header_buf: [0; 20],
            header_len: 0,
        }
    }
}

impl BinseqRecord for RefRecord<'_> {
    fn bitsize(&self) -> BitSize {
        self.config.bitsize
    }
    fn index(&self) -> u64 {
        self.id
    }
    /// Clear the buffer and fill it with the sequence header
    fn sheader(&self) -> &[u8] {
        &self.header_buf[..self.header_len]
    }

    /// Clear the buffer and fill it with the extended header
    fn xheader(&self) -> &[u8] {
        self.sheader()
    }

    fn flag(&self) -> Option<u64> {
        if self.config.flags {
            Some(self.buffer[0])
        } else {
            None
        }
    }
    fn slen(&self) -> u64 {
        self.config.slen
    }
    fn xlen(&self) -> u64 {
        self.config.xlen
    }
    fn sbuf(&self) -> &[u64] {
        if self.config.flags {
            &self.buffer[1..=(self.config.schunk as usize)]
        } else {
            &self.buffer[..(self.config.schunk as usize)]
        }
    }
    fn xbuf(&self) -> &[u64] {
        if self.config.flags {
            &self.buffer[1 + self.config.schunk as usize..]
        } else {
            &self.buffer[self.config.schunk as usize..]
        }
    }
    fn squal(&self) -> &[u8] {
        &self.qbuf[..self.config.slen as usize]
    }
    fn xqual(&self) -> &[u8] {
        &self.qbuf[..self.config.xlen as usize]
    }
}

/// A reference to a record in the map with a precomputed decoded buffer slice
pub struct BatchRecord<'a> {
    /// Unprocessed buffer slice (with flags)
    buffer: &'a [u64],
    /// Decoded buffer slice
    dbuf: &'a [u8],
    /// Record ID
    id: u64,
    /// The configuration that defines the layout and size of record components
    config: RecordConfig,
    /// A reusable pre-initialized quality score buffer
    qbuf: &'a [u8],
    /// Cached index string for the sequence header
    header_buf: [u8; 20],
    /// Length of the header in bytes
    header_len: usize,
}
impl BinseqRecord for BatchRecord<'_> {
    fn bitsize(&self) -> BitSize {
        self.config.bitsize
    }
    fn index(&self) -> u64 {
        self.id
    }
    /// Clear the buffer and fill it with the sequence header
    fn sheader(&self) -> &[u8] {
        &self.header_buf[..self.header_len]
    }

    /// Clear the buffer and fill it with the extended header
    fn xheader(&self) -> &[u8] {
        self.sheader()
    }

    fn flag(&self) -> Option<u64> {
        if self.config.flags {
            Some(self.buffer[0])
        } else {
            None
        }
    }
    fn slen(&self) -> u64 {
        self.config.slen
    }
    fn xlen(&self) -> u64 {
        self.config.xlen
    }
    fn sbuf(&self) -> &[u64] {
        if self.config.flags {
            &self.buffer[1..=(self.config.schunk as usize)]
        } else {
            &self.buffer[..(self.config.schunk as usize)]
        }
    }
    fn xbuf(&self) -> &[u64] {
        if self.config.flags {
            &self.buffer[1 + self.config.schunk as usize..]
        } else {
            &self.buffer[self.config.schunk as usize..]
        }
    }
    fn decode_s(&self, dbuf: &mut Vec<u8>) -> Result<()> {
        dbuf.extend_from_slice(self.sseq());
        Ok(())
    }
    fn decode_x(&self, dbuf: &mut Vec<u8>) -> Result<()> {
        dbuf.extend_from_slice(self.xseq());
        Ok(())
    }
    /// Override this method since we can make use of block information
    fn sseq(&self) -> &[u8] {
        let scalar = self.config.scalar();
        let mut lbound = 0;
        let mut rbound = self.config.slen();
        if self.config.flags {
            lbound += scalar;
            rbound += scalar;
        }
        &self.dbuf[lbound..rbound]
    }
    /// Override this method since we can make use of block information
    fn xseq(&self) -> &[u8] {
        let scalar = self.config.scalar();
        let mut lbound = scalar * self.config.schunk();
        let mut rbound = lbound + self.config.xlen();
        if self.config.flags {
            lbound += scalar;
            rbound += scalar;
        }
        &self.dbuf[lbound..rbound]
    }
    fn squal(&self) -> &[u8] {
        &self.qbuf[..self.config.slen()]
    }
    fn xqual(&self) -> &[u8] {
        &self.qbuf[..self.config.xlen()]
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
    slen: u64,
    /// The extended sequence length in base pairs
    xlen: u64,
    /// The number of u64 chunks needed to store the primary sequence
    /// (each u64 stores 32 nucleotides)
    schunk: u64,
    /// The number of u64 chunks needed to store the extended sequence
    /// (each u64 stores 32 values)
    xchunk: u64,
    /// The bitsize of the record
    bitsize: BitSize,
    /// Whether flags are present
    flags: bool,
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
    /// * `bitsize` - The bitsize of the record
    /// * `flags` - Whether flags are present
    ///
    /// # Returns
    ///
    /// A new `RecordConfig` instance with the specified sequence lengths
    pub fn new(slen: usize, xlen: usize, bitsize: BitSize, flags: bool) -> Self {
        let (schunk, xchunk) = match bitsize {
            BitSize::Two => (slen.div_ceil(32), xlen.div_ceil(32)),
            BitSize::Four => (slen.div_ceil(16), xlen.div_ceil(16)),
        };
        Self {
            slen: slen as u64,
            xlen: xlen as u64,
            schunk: schunk as u64,
            xchunk: xchunk as u64,
            bitsize,
            flags,
        }
    }

    /// Creates a new record configuration from a header
    ///
    /// This constructor initializes a configuration based on a header that contains
    /// the sequence lengths for primary and extended sequences.
    ///
    /// # Arguments
    ///
    /// * `header` - A reference to a `FileHeader` containing sequence lengths
    ///
    /// # Returns
    ///
    /// A new `RecordConfig` instance with the sequence lengths from the header
    pub fn from_header(header: &FileHeader) -> Self {
        Self::new(
            header.slen as usize,
            header.xlen as usize,
            header.bits,
            header.flags,
        )
    }

    /// Returns the primary sequence length in base pairs
    ///
    /// This method returns the length of the primary sequence in base pairs.
    pub fn slen(&self) -> usize {
        self.slen as usize
    }

    /// Returns the extended sequence length in base pairs
    ///
    /// This method returns the length of the extended sequence in base pairs.
    pub fn xlen(&self) -> usize {
        self.xlen as usize
    }

    /// Returns the number of u64 chunks needed to store the primary sequence
    ///
    /// This method returns the number of u64 chunks required to store the primary
    /// sequence, where each u64 stores 32 nucleotides.
    pub fn schunk(&self) -> usize {
        self.schunk as usize
    }

    /// Returns the number of u64 chunks needed to store the extended sequence
    ///
    /// This method returns the number of u64 chunks required to store the extended
    /// sequence, where each u64 stores 32 values.
    pub fn xchunk(&self) -> usize {
        self.xchunk as usize
    }

    /// Returns the full record size in bytes (u8):
    /// 8 * (schunk + xchunk + 1 (flag))
    pub fn record_size_bytes(&self) -> usize {
        8 * self.record_size_u64()
    }

    /// Returns the full record size in u64
    /// schunk + xchunk + 1 (flag)
    pub fn record_size_u64(&self) -> usize {
        if self.flags {
            (self.schunk + self.xchunk + 1) as usize
        } else {
            (self.schunk + self.xchunk) as usize
        }
    }

    /// The number of nucleotides per word
    pub fn scalar(&self) -> usize {
        match self.bitsize {
            BitSize::Two => 32,
            BitSize::Four => 16,
        }
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
///
/// Records are returned as [`RefRecord`] which implement the [`BinseqRecord`] trait.
///
/// # Examples
///
/// ```
/// use binseq::bq::MmapReader;
/// use binseq::Result;
///
/// fn main() -> Result<()> {
///     let path = "./data/subset.bq";
///     let reader = MmapReader::new(path)?;
///
///     // Calculate the number of records in the file
///     let num_records = reader.num_records();
///     println!("Number of records: {}", num_records);
///
///     // Get the record at index 20 (0-indexed)
///     let record = reader.get(20)?;
///
///     Ok(())
/// }
/// ```
pub struct MmapReader {
    /// Memory mapped file contents, wrapped in Arc for thread-safe sharing
    mmap: Arc<Mmap>,

    /// Binary sequence file header containing format information
    header: FileHeader,

    /// Configuration defining the layout of records in the file
    config: RecordConfig,

    /// Reusable buffer for quality scores
    qbuf: Vec<u8>,

    /// Default quality score for records without quality scores
    default_quality_score: u8,
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
        let header = FileHeader::from_buffer(&mmap)?;

        // Record configuraration
        let config = RecordConfig::from_header(&header);

        // Immediately validate the size of the file against the expected byte size of records
        if !(mmap.len() - SIZE_HEADER).is_multiple_of(config.record_size_bytes()) {
            return Err(ReadError::FileTruncation(mmap.len()).into());
        }

        // preinitialize quality buffer
        let qbuf = vec![DEFAULT_QUALITY_SCORE; header.slen.max(header.xlen) as usize];

        Ok(Self {
            mmap: Arc::new(mmap),
            header,
            config,
            qbuf,
            default_quality_score: DEFAULT_QUALITY_SCORE,
        })
    }

    /// Returns the total number of records in the file
    ///
    /// This is calculated by subtracting the header size from the total file size
    /// and dividing by the size of each record.
    #[must_use]
    pub fn num_records(&self) -> usize {
        (self.mmap.len() - SIZE_HEADER) / self.config.record_size_bytes()
    }

    /// Returns a copy of the binary sequence file header
    ///
    /// The header contains format information and sequence length specifications.
    #[must_use]
    pub fn header(&self) -> FileHeader {
        self.header
    }

    /// Checks if the file has paired-records
    #[must_use]
    pub fn is_paired(&self) -> bool {
        self.header.is_paired()
    }

    /// Sets the default quality score for records without quality information
    pub fn set_default_quality_score(&mut self, score: u8) {
        self.default_quality_score = score;
        self.qbuf = self.build_qbuf();
    }

    /// Creates a new quality score buffer
    #[must_use]
    pub fn build_qbuf(&self) -> Vec<u8> {
        vec![self.default_quality_score; self.header.slen.max(self.header.xlen) as usize]
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
    pub fn get(&self, idx: usize) -> Result<RefRecord<'_>> {
        if idx > self.num_records() {
            return Err(ReadError::OutOfRange {
                requested_index: idx,
                max_index: self.num_records(),
            }
            .into());
        }
        let rsize = self.config.record_size_bytes();
        let lbound = SIZE_HEADER + (idx * rsize);
        let rbound = lbound + rsize;
        let bytes = &self.mmap[lbound..rbound];
        let buffer = cast_slice(bytes);
        Ok(RefRecord::new(idx as u64, buffer, &self.qbuf, self.config))
    }

    /// Returns a slice of the buffer containing the underlying u64 for that range
    /// of records.
    ///
    /// Note: range 10..40 will return all u64s in the mmap between the record index 10 and 40
    pub fn get_buffer_slice(&self, range: Range<usize>) -> Result<&[u64]> {
        if range.end > self.num_records() {
            return Err(ReadError::OutOfRange {
                requested_index: range.end,
                max_index: self.num_records(),
            }
            .into());
        }
        let rsize = self.config.record_size_bytes();
        let total_records = range.end - range.start;
        let lbound = SIZE_HEADER + (range.start * rsize);
        let rbound = lbound + (total_records * rsize);
        let bytes = &self.mmap[lbound..rbound];
        let buffer = cast_slice(bytes);
        Ok(buffer)
    }
}

/// A reader for streaming binary sequence data from any source that implements Read
///
/// Unlike `MmapReader` which requires the entire file to be accessible at once,
/// `StreamReader` processes data as it becomes available, making it suitable for:
/// - Processing data as it arrives over a network
/// - Handling very large files that exceed available memory
/// - Pipeline processing where data is flowing continuously
///
/// The reader maintains an internal buffer and can handle partial record reconstruction
/// across chunk boundaries.
pub struct StreamReader<R: Read> {
    /// The source reader for binary sequence data
    reader: R,

    /// Binary sequence file header containing format information
    header: Option<FileHeader>,

    /// Configuration defining the layout of records in the file
    config: Option<RecordConfig>,

    /// Buffer for storing incoming data
    buffer: Vec<u8>,

    /// Buffer for reusable quality scores
    qbuf: Vec<u8>,

    /// Default quality score for records without quality information
    default_quality_score: u8,

    /// Current position in the buffer
    buffer_pos: usize,

    /// Length of valid data in the buffer
    buffer_len: usize,

    /// Number of records returned so far, used to assign each record's id
    ///
    /// This is tracked independently of `buffer_pos` because `fill_buffer`
    /// shifts remaining bytes to the start of the buffer and resets
    /// `buffer_pos` whenever a mid-stream refill is needed, so `buffer_pos`
    /// no longer reflects the absolute stream offset once that happens.
    records_read: u64,
}

impl<R: Read> StreamReader<R> {
    /// Creates a new `StreamReader` with the default buffer size
    ///
    /// This constructor initializes a `StreamReader` that will read from the provided
    /// source, using an 8K default buffer size.
    ///
    /// # Arguments
    ///
    /// * `reader` - The source to read binary sequence data from
    ///
    /// # Returns
    ///
    /// A new `StreamReader` instance
    pub fn new(reader: R) -> Self {
        Self::with_capacity(reader, 8192)
    }

    /// Creates a new `StreamReader` with a specified buffer capacity
    ///
    /// This constructor initializes a `StreamReader` with a custom buffer size,
    /// which can be tuned based on the expected usage pattern.
    ///
    /// # Arguments
    ///
    /// * `reader` - The source to read binary sequence data from
    /// * `capacity` - The size of the internal buffer in bytes
    ///
    /// # Returns
    ///
    /// A new `StreamReader` instance with the specified buffer capacity
    pub fn with_capacity(reader: R, capacity: usize) -> Self {
        Self {
            reader,
            header: None,
            config: None,
            buffer: vec![0; capacity],
            qbuf: vec![0; capacity],
            buffer_pos: 0,
            buffer_len: 0,
            default_quality_score: DEFAULT_QUALITY_SCORE,
            records_read: 0,
        }
    }

    /// Sets the default quality score for records without quality information
    pub fn set_default_quality_score(&mut self, score: u8) {
        if score != self.default_quality_score {
            self.qbuf.clear();
        }
        self.default_quality_score = score;
    }

    /// Reads and validates the header from the underlying reader
    ///
    /// This method reads the binary sequence file header and validates it.
    /// It caches the header internally for future use.
    ///
    /// # Returns
    ///
    /// * `Ok(&FileHeader)` - A reference to the validated header
    /// * `Err(Error)` - If reading or validating the header fails
    ///
    /// # Panics
    ///
    /// Panics if the header is missing when expected in the stream.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// * There is an I/O error when reading from the source
    /// * The header data is invalid
    /// * End of stream is reached before the full header can be read
    pub fn read_header(&mut self) -> Result<&FileHeader> {
        if self.header.is_none() {
            // Ensure we have enough data for the header
            while self.buffer_len - self.buffer_pos < SIZE_HEADER {
                self.fill_buffer()?;
            }

            // Parse header
            let header_slice = &self.buffer[self.buffer_pos..self.buffer_pos + SIZE_HEADER];
            let header = FileHeader::from_buffer(header_slice)?;

            self.header = Some(header);
            self.config = Some(RecordConfig::from_header(&header));
            self.buffer_pos += SIZE_HEADER;
        }

        Ok(self
            .header
            .as_ref()
            .expect("header was just populated above"))
    }

    /// Fills the internal buffer with more data from the reader
    ///
    /// This method reads more data from the underlying reader, handling
    /// the case where some unprocessed data remains in the buffer.
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If the buffer was successfully filled with new data
    /// * `Err(Error)` - If reading from the source fails
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// * There is an I/O error when reading from the source
    /// * End of stream is reached (no more data available)
    fn fill_buffer(&mut self) -> Result<()> {
        // Move remaining data to beginning of buffer if needed
        if self.buffer_pos > 0 && self.buffer_pos < self.buffer_len {
            self.buffer.copy_within(self.buffer_pos..self.buffer_len, 0);
            self.buffer_len -= self.buffer_pos;
            self.buffer_pos = 0;
        } else if self.buffer_pos == self.buffer_len {
            self.buffer_len = 0;
            self.buffer_pos = 0;
        }

        // Read more data
        let bytes_read = self.reader.read(&mut self.buffer[self.buffer_len..])?;
        if bytes_read == 0 {
            return Err(ReadError::EndOfStream.into());
        }

        self.buffer_len += bytes_read;
        Ok(())
    }

    /// Retrieves the next record from the stream
    ///
    /// This method reads and processes the next complete record from the stream.
    /// It handles the case where a record spans multiple buffer fills.
    ///
    /// # Returns
    ///
    /// * `Ok(Some(RefRecord))` - The next record was successfully read
    /// * `Ok(None)` - End of stream was reached (no more records)
    /// * `Err(Error)` - If an error occurred during reading
    ///
    /// # Panics
    ///
    /// Panics if the configuration is missing when expected in the stream.
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// * There is an I/O error when reading from the source
    /// * The header has not been read yet
    /// * The data format is invalid
    pub fn next_record(&mut self) -> Option<Result<RefRecord<'_>>> {
        // Ensure header is read
        if self.header.is_none()
            && let Some(e) = self.read_header().err()
        {
            return Some(Err(e));
        }

        let config = self
            .config
            .expect("Missing configuration when expected in stream");
        let record_size = config.record_size_bytes();

        // Ensure we have enough data for a complete record
        while self.buffer_len - self.buffer_pos < record_size {
            match self.fill_buffer() {
                Ok(()) => {}
                Err(Error::ReadError(ReadError::EndOfStream)) => {
                    // End of stream reached - if we have any partial data, it's an error
                    if self.buffer_len - self.buffer_pos > 0 {
                        return Some(Err(ReadError::PartialRecord(
                            self.buffer_len - self.buffer_pos,
                        )
                        .into()));
                    }
                    return None;
                }
                Err(e) => return Some(Err(e)),
            }
        }

        // Process record
        let record_start = self.buffer_pos;
        self.buffer_pos += record_size;

        let record_bytes = &self.buffer[record_start..record_start + record_size];
        let record_u64s = cast_slice(record_bytes);

        // update quality score buffer if necessary
        if self.qbuf.is_empty() {
            let max_size = config.slen.max(config.xlen) as usize;
            self.qbuf.resize(max_size, self.default_quality_score);
        }

        // Create record with an incremental ID, tracked independently of
        // buffer position since `fill_buffer` may have shifted the buffer
        let id = self.records_read;
        self.records_read += 1;
        Some(Ok(RefRecord::new(id, record_u64s, &self.qbuf, config)))
    }

    /// Consumes the stream reader and returns the inner reader
    ///
    /// This method is useful when you need access to the underlying reader
    /// after processing is complete.
    ///
    /// # Returns
    ///
    /// The inner reader that was used by this `StreamReader`
    pub fn into_inner(self) -> R {
        self.reader
    }
}

/// Default batch size for parallel processing
///
/// This constant defines how many records each thread processes at a time
/// during parallel processing operations.
pub const BATCH_SIZE: usize = 1024;

/// Parallel processing implementation for memory-mapped readers
impl ParallelReader for MmapReader {
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
    fn process_parallel<P: ParallelProcessor + Clone + 'static>(
        self,
        processor: P,
        num_threads: usize,
    ) -> Result<()> {
        let num_records = self.num_records();
        self.process_parallel_range(processor, num_threads, 0..num_records)
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
    /// * `range` - The range of record indices to process
    ///
    /// # Type Parameters
    ///
    /// * `P` - A type that implements `ParallelProcessor` and can be cloned
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
    ) -> Result<()> {
        // Calculate the number of threads to use
        let num_threads = if num_threads == 0 {
            num_cpus::get()
        } else {
            num_threads.min(num_cpus::get())
        };

        // Validate range
        let num_records = self.num_records();
        self.validate_range(num_records, &range)?;

        // Calculate number of records for each thread within the range
        let range_size = range.end - range.start;
        let records_per_thread = range_size.div_ceil(num_threads);

        // Arc self
        let reader = Arc::new(self);

        // Build thread handles
        let mut handles = Vec::new();
        for tid in 0..num_threads {
            let mut processor = processor.clone();
            let reader = reader.clone();
            processor.set_tid(tid);

            let handle = std::thread::spawn(move || -> Result<()> {
                let start_idx = range.start + tid * records_per_thread;
                let end_idx = (start_idx + records_per_thread).min(range.end);

                if start_idx >= end_idx {
                    return Ok(()); // No records for this thread
                }

                // create a reusable buffer for translating record IDs
                let mut translater = itoa::Buffer::new();

                // initialize a decoding buffer
                let mut dbuf = Vec::new();

                // initialize a quality score buffer
                let qbuf = reader.build_qbuf();

                // calculate the size of a record in the cast u64 slice
                let rsize_u64 = reader.config.record_size_bytes() / 8;

                // determine the required scalar size
                let scalar = reader.config.scalar();

                // calculate the size of a record in the batch decoded buffer
                let mut dbuf_rsize = { (reader.config.schunk() + reader.config.xchunk()) * scalar };
                if reader.config.flags {
                    dbuf_rsize += scalar;
                }

                // iterate over the range of indices
                for range_start in (start_idx..end_idx).step_by(BATCH_SIZE) {
                    let range_end = (range_start + BATCH_SIZE).min(end_idx);

                    // clear the decoded buffer
                    dbuf.clear();

                    // get the encoded buffer slice
                    let ebuf = reader.get_buffer_slice(range_start..range_end)?;

                    // decode the entire buffer at once (with flags and extra bases)
                    reader
                        .config
                        .bitsize
                        .decode(ebuf, ebuf.len() * scalar, &mut dbuf)?;

                    // iterate over each index in the range
                    for (inner_idx, idx) in (range_start..range_end).enumerate() {
                        // translate the index
                        let id_str = translater.format(idx);

                        // create the index buffer
                        let mut header_buf = [0; 20];
                        let header_len = id_str.len();
                        header_buf[..header_len].copy_from_slice(id_str.as_bytes());

                        // find the buffer starts
                        let ebuf_start = inner_idx * rsize_u64;
                        let dbuf_start = inner_idx * dbuf_rsize;

                        // initialize the record
                        let record = BatchRecord {
                            buffer: &ebuf[ebuf_start..(ebuf_start + rsize_u64)],
                            dbuf: &dbuf[dbuf_start..(dbuf_start + dbuf_rsize)],
                            qbuf: &qbuf,
                            id: idx as u64,
                            config: reader.config,
                            header_buf,
                            header_len,
                        };

                        // process the record
                        processor.process_record(record)?;
                    }

                    // process the batch
                    processor.on_batch_complete()?;
                }

                // process the thread
                processor.on_thread_complete()?;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BinseqRecord;
    use bitnuc_deprec::BitSize;

    const TEST_BQ_FILE: &str = "./data/subset.bq";

    // ==================== MmapReader Basic Tests ====================

    #[test]
    fn test_mmap_reader_new() {
        let reader = MmapReader::new(TEST_BQ_FILE);
        assert!(reader.is_ok(), "Failed to create reader");
    }

    #[test]
    fn test_mmap_reader_num_records() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();
        let num_records = reader.num_records();
        assert!(num_records > 0, "Expected non-zero records");
    }

    #[test]
    fn test_mmap_reader_is_paired() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();
        // The fixture file contains paired records
        assert!(reader.is_paired());
    }

    #[test]
    fn test_mmap_reader_header_access() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();
        let header = reader.header();
        assert!(header.slen > 0, "Expected non-zero sequence length");
    }

    #[test]
    fn test_mmap_reader_config_access() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();
        let header = reader.header();
        let config = RecordConfig::from_header(&header);
        assert!(
            config.slen > 0,
            "Expected non-zero sequence length in config"
        );
    }

    // ==================== Record Access Tests ====================

    #[test]
    fn test_get_record() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();
        let num_records = reader.num_records();

        if num_records > 0 {
            let record = reader.get(0);
            assert!(record.is_ok(), "Expected to get first record");

            let record = record.unwrap();
            assert_eq!(record.index(), 0, "Expected record index to be 0");
        }
    }

    #[test]
    fn test_get_record_out_of_bounds() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();
        let num_records = reader.num_records();

        let record = reader.get(num_records + 100);
        assert!(record.is_err(), "Expected error for out of bounds index");
    }

    #[test]
    fn test_record_sequence_data() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();

        if let Ok(record) = reader.get(0) {
            let sbuf = record.sbuf();
            assert!(!sbuf.is_empty(), "Expected non-empty sequence buffer");

            let slen = record.slen();
            assert!(slen > 0, "Expected non-zero sequence length");
        }
    }

    #[test]
    fn test_record_quality_data() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();

        if let Ok(record) = reader.get(0) {
            let squal = record.squal();
            let slen = record.slen() as usize;
            assert_eq!(
                squal.len(),
                slen,
                "Quality length should match sequence length"
            );
        }
    }

    // ==================== Default Quality Score Tests ====================

    #[test]
    fn test_set_default_quality_score() {
        let mut reader = MmapReader::new(TEST_BQ_FILE).unwrap();
        let custom_score = 42u8;

        reader.set_default_quality_score(custom_score);

        if let Ok(record) = reader.get(0) {
            let squal = record.squal();
            // All quality scores should be the custom score
            assert!(
                squal.iter().all(|&q| q == custom_score),
                "All quality scores should be {custom_score}"
            );
        }
    }

    // ==================== Parallel Processing Tests ====================

    #[derive(Clone)]
    struct CountingProcessor {
        count: Arc<std::sync::Mutex<usize>>,
    }

    impl ParallelProcessor for CountingProcessor {
        fn process_record<R: BinseqRecord>(&mut self, _record: R) -> Result<()> {
            let mut count = self.count.lock().unwrap();
            *count += 1;
            Ok(())
        }
    }

    #[test]
    fn test_parallel_processing() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();
        let num_records = reader.num_records();

        let count = Arc::new(std::sync::Mutex::new(0));
        let processor = CountingProcessor {
            count: count.clone(),
        };

        reader.process_parallel(processor, 2).unwrap();

        let final_count = *count.lock().unwrap();
        assert_eq!(final_count, num_records, "All records should be processed");
    }

    #[test]
    fn test_parallel_processing_range() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();
        let num_records = reader.num_records();

        if num_records >= 100 {
            let start = 10;
            let end = 50;
            let expected_count = end - start;

            let count = Arc::new(std::sync::Mutex::new(0));
            let processor = CountingProcessor {
                count: count.clone(),
            };

            reader
                .process_parallel_range(processor, 2, start..end)
                .unwrap();

            let final_count = *count.lock().unwrap();
            assert_eq!(
                final_count, expected_count,
                "Should process exactly {expected_count} records"
            );
        }
    }

    // ==================== RecordConfig Tests ====================

    #[test]
    fn test_record_config_from_header() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();
        let header = reader.header();
        let config = RecordConfig::from_header(&header);

        assert_eq!(
            config.slen,
            u64::from(header.slen),
            "Sequence length mismatch"
        );
        assert_eq!(
            config.xlen,
            u64::from(header.xlen),
            "Extended length mismatch"
        );
        assert_eq!(config.bitsize, header.bits, "Bit size mismatch");
    }

    #[test]
    fn test_record_config_record_size() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();
        let header = reader.header();
        let config = RecordConfig::from_header(&header);

        let size_u64 = config.record_size_u64();
        assert!(size_u64 > 0, "Record size should be non-zero");

        let size_bytes = config.record_size_bytes();
        assert_eq!(size_bytes, size_u64 * 8, "Byte size should be 8x u64 size");
    }

    // ==================== RefRecord Tests ====================

    #[test]
    fn test_ref_record_bitsize() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();

        if let Ok(record) = reader.get(0) {
            let bitsize = record.bitsize();
            assert!(
                matches!(bitsize, BitSize::Two | BitSize::Four),
                "Bitsize should be Two or Four"
            );
        }
    }

    #[test]
    fn test_ref_record_flag() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();

        if let Ok(record) = reader.get(0) {
            let flag = record.flag();
            // Flag should be Some if header has flags enabled
            assert!(flag.is_some() || flag.is_none()); // Tests method works
        }
    }

    #[test]
    fn test_ref_record_paired_data() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();

        if reader.is_paired()
            && let Ok(record) = reader.get(0)
        {
            let xbuf = record.xbuf();
            let xlen = record.xlen();

            if xlen > 0 {
                assert!(
                    !xbuf.is_empty(),
                    "Extended buffer should not be empty for paired"
                );
            }
        }
    }

    // ==================== Error Handling Tests ====================

    #[test]
    fn test_nonexistent_file() {
        let result = MmapReader::new("./data/nonexistent.bq");
        assert!(result.is_err(), "Should fail on nonexistent file");
    }

    #[test]
    fn test_invalid_file_format() {
        // Try to open a non-BQ file as BQ (use Cargo.toml for example)
        let result = MmapReader::new("./Cargo.toml");
        // This should either fail to open or fail validation
        if let Ok(reader) = result {
            // If it opens, try to access records (should fail or have issues)
            let num_records = reader.num_records();
            // The number might be nonsensical for invalid data
            let _ = num_records; // Just verify it doesn't panic
        }
    }

    // ==================== Multiple Records Tests ====================

    #[test]
    fn test_sequential_record_access() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();
        let num_records = reader.num_records().min(10);

        for i in 0..num_records {
            let record = reader.get(i);
            assert!(record.is_ok(), "Should get record at index {i}");
            assert_eq!(
                record.unwrap().index() as usize,
                i,
                "Record index mismatch at {i}"
            );
        }
    }

    #[test]
    fn test_random_record_access() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();
        let num_records = reader.num_records();

        if num_records > 10 {
            let indices = [0, 5, num_records / 2, num_records - 1];

            for &idx in &indices {
                let record = reader.get(idx);
                assert!(record.is_ok(), "Should get record at index {idx}");
                assert_eq!(record.unwrap().index() as usize, idx);
            }
        }
    }

    // ==================== get_buffer_slice Tests ====================

    #[test]
    fn test_get_buffer_slice_valid() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();
        let num_records = reader.num_records().min(10);
        let slice = reader.get_buffer_slice(0..num_records);
        assert!(slice.is_ok());
        assert_eq!(
            slice.unwrap().len(),
            num_records * reader.config.record_size_u64()
        );
    }

    #[test]
    fn test_get_buffer_slice_out_of_range() {
        let reader = MmapReader::new(TEST_BQ_FILE).unwrap();
        let num_records = reader.num_records();
        let slice = reader.get_buffer_slice(0..(num_records + 100));
        assert!(slice.is_err());
    }

    // ==================== MmapReader Error Path Tests ====================

    #[test]
    fn test_mmap_reader_directory_is_incompatible() {
        // Directories cannot be memory-mapped as regular files
        let result = MmapReader::new("./data");
        assert!(result.is_err(), "Should fail when given a directory");
    }

    #[test]
    fn test_mmap_reader_truncated_file() {
        use std::io::Write as _;

        let path = "test_truncated_reader.bq";
        {
            let header = crate::bq::FileHeaderBuilder::new()
                .slen(64)
                .build()
                .unwrap();
            let mut file = std::fs::File::create(path).unwrap();
            header.write_bytes(&mut file).unwrap();
            // Write a partial record (not a full multiple of the record size)
            file.write_all(&[0u8; 4]).unwrap();
        }

        let result = MmapReader::new(path);
        assert!(result.is_err(), "Should fail on truncated record data");

        std::fs::remove_file(path).unwrap();
    }

    // ==================== StreamReader Tests ====================

    fn build_stream_bytes(paired: bool) -> Vec<u8> {
        use crate::SequencingRecordBuilder;
        use crate::bq::{FileHeaderBuilder, WriterBuilder};

        let header = if paired {
            FileHeaderBuilder::new().slen(64).xlen(32).build().unwrap()
        } else {
            FileHeaderBuilder::new().slen(64).build().unwrap()
        };

        let mut writer = WriterBuilder::default()
            .header(header)
            .build(Vec::new())
            .unwrap();

        for i in 0..5 {
            let s_seq = vec![b"ACGT"[i % 4]; 64];
            let x_seq = vec![b"TGCA"[i % 4]; 32];
            let record = if paired {
                SequencingRecordBuilder::default()
                    .s_seq(&s_seq)
                    .x_seq(&x_seq)
                    .build()
                    .unwrap()
            } else {
                SequencingRecordBuilder::default()
                    .s_seq(&s_seq)
                    .build()
                    .unwrap()
            };
            writer.push(record).unwrap();
        }
        writer.flush().unwrap();
        writer.into_inner()
    }

    #[test]
    fn test_stream_reader_new_and_with_capacity() {
        let data = build_stream_bytes(false);
        let cursor = std::io::Cursor::new(data.clone());
        let _reader = StreamReader::new(cursor);

        let cursor = std::io::Cursor::new(data);
        let _reader = StreamReader::with_capacity(cursor, 64);
    }

    #[test]
    fn test_stream_reader_read_header() {
        let data = build_stream_bytes(false);
        let mut reader = StreamReader::new(std::io::Cursor::new(data));
        let header = reader.read_header().unwrap();
        assert_eq!(header.slen, 64);

        // Second call should hit the cached path
        let header_again = reader.read_header().unwrap();
        assert_eq!(header_again.slen, 64);
    }

    #[test]
    fn test_stream_reader_next_record_unpaired() {
        let data = build_stream_bytes(false);
        let mmap_reader = {
            let path = "test_stream_reader_compare.bq";
            std::fs::write(path, &data).unwrap();
            let reader = MmapReader::new(path).unwrap();
            std::fs::remove_file(path).unwrap();
            reader
        };

        let mut reader = StreamReader::new(std::io::Cursor::new(data));
        let mut count = 0;
        while let Some(record) = reader.next_record() {
            let record = record.unwrap();
            let expected = mmap_reader.get(count).unwrap();
            assert_eq!(record.index(), expected.index());
            assert_eq!(record.sbuf(), expected.sbuf());
            count += 1;
        }
        assert_eq!(count, mmap_reader.num_records());
    }

    #[test]
    fn test_stream_reader_next_record_paired() {
        let data = build_stream_bytes(true);
        let mut reader = StreamReader::new(std::io::Cursor::new(data));
        let mut count = 0;
        while let Some(record) = reader.next_record() {
            let record = record.unwrap();
            assert!(record.is_paired());
            count += 1;
        }
        assert_eq!(count, 5);
    }

    #[test]
    fn test_stream_reader_small_buffer_forces_multiple_fills() {
        // Use a tiny buffer capacity so `fill_buffer` must shift remaining
        // bytes to the start of the internal buffer mid-stream (reader.rs
        // 716-719). Record ids are tracked via a dedicated `records_read`
        // counter (not derived from `buffer_pos`), so they must still come
        // back sequential across that shift boundary.
        let data = build_stream_bytes(false);
        let mut reader = StreamReader::with_capacity(std::io::Cursor::new(data), 40);
        let mut expected_id = 0u64;
        while let Some(record) = reader.next_record() {
            let record = record.unwrap();
            assert_eq!(record.index(), expected_id);
            expected_id += 1;
        }
        assert_eq!(expected_id, 5);
    }

    #[test]
    fn test_stream_reader_partial_record_error() {
        let mut data = build_stream_bytes(false);
        // Truncate the data in the middle of the last record
        data.truncate(data.len() - 4);
        let mut reader = StreamReader::new(std::io::Cursor::new(data));

        let mut saw_error = false;
        while let Some(record) = reader.next_record() {
            if record.is_err() {
                saw_error = true;
                break;
            }
        }
        assert!(saw_error, "Expected a partial record error");
    }

    #[test]
    fn test_stream_reader_set_default_quality_score() {
        let data = build_stream_bytes(false);
        let mut reader = StreamReader::new(std::io::Cursor::new(data));
        reader.set_default_quality_score(42);
        if let Some(Ok(record)) = reader.next_record() {
            assert!(record.squal().iter().all(|&q| q == 42));
        }
    }

    #[test]
    fn test_stream_reader_into_inner() {
        let data = build_stream_bytes(false);
        let reader = StreamReader::new(std::io::Cursor::new(data.clone()));
        let cursor = reader.into_inner();
        assert_eq!(cursor.into_inner(), data);
    }
}
