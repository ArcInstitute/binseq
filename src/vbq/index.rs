//! # VBQ Index Format
//!
//! This module implements the embedded index format for VBQ files.
//!
//! ## Format Changes (v0.7.0+)
//!
//! **BREAKING CHANGE**: The VBQ index is now embedded at the end of VBQ files instead of
//! being stored in separate `.vqi` files. This improves portability and eliminates the
//! need to manage auxiliary files.
//!
//! ## Embedded Index Structure
//!
//! The index is located at the end of the VBQ file with this layout:
//!
//! ```text
//! [VBQ Data Blocks][Compressed Index][Index Size (u64)][INDEX_END_MAGIC (u64)]
//! ```
//!
//! Where:
//! - **Compressed Index**: ZSTD-compressed index data (`IndexHeader` + `BlockRanges`)
//! - **Index Size**: 8 bytes indicating size of compressed index data
//! - **`INDEX_END_MAGIC`**: 8 bytes (`0x444E455845444E49` = "INDEXEND")
//!
//! ## Index Contents
//!
//! The compressed index contains:
//! 1. **`IndexHeader`** (32 bytes): Metadata about the indexed file
//! 2. **`BlockRange` entries** (32 bytes each): One per data block
//!
//! ## Key Changes from v0.6.x
//!
//! - Index moved from separate `.vqi` files into VBQ files
//! - Cumulative record counts changed from `u32` to `u64`
//! - Support for files with more than 4 billion records

use std::{
    fs::File,
    io::{BufReader, BufWriter, Cursor, Read, Write},
    path::Path,
};

use byteorder::{ByteOrder, LittleEndian};
use zstd::{Decoder, Encoder};

use super::{
    BlockHeader, FileHeader,
    header::{SIZE_BLOCK_HEADER, SIZE_HEADER},
};
use crate::error::{IndexError, Result};

/// Size of `BlockRange` in bytes
pub const SIZE_BLOCK_RANGE: usize = 32;
/// Size of `IndexHeader` in bytes
pub const INDEX_HEADER_SIZE: usize = 32;
/// Magic number to designate index (VBQINDEX)
#[allow(clippy::unreadable_literal)]
pub const INDEX_MAGIC: u64 = 0x5845444e49514256;
/// Magic number to designate end of index (INDEXEND)
#[allow(clippy::unreadable_literal)]
pub const INDEX_END_MAGIC: u64 = 0x444E455845444E49;
/// Index Block Reservation
pub const INDEX_RESERVATION: [u8; 4] = [42; 4];

/// Descriptor of the dimensions of a block in a VBQ file
///
/// A `BlockRange` contains metadata about a single block within a VBQ file,
/// including its position, size, and record count. This information enables
/// efficient random access to blocks without scanning the entire file.
///
/// Block ranges are stored in a `BlockIndex` to form a complete index of a VBQ file.
/// Each range is serialized to a fixed-size 32-byte structure when stored in the embedded index.
///
/// ## Format Changes (v0.7.0+)
///
/// - `cumulative_records` field changed from `u32` to `u64`
/// - Supports files with more than 4 billion records
/// - Reserved bytes reduced from 8 to 4 bytes
///
/// # Examples
///
/// ```rust
/// use binseq::vbq::BlockRange;
///
/// // Create a new block range
/// let range = BlockRange::new(
///     1024,                  // Starting offset in the file (bytes)
///     8192,                  // Length of the block (bytes)
///     1000,                  // Number of records in this block
///     5000                   // Cumulative number of records up to this block (now u64)
/// );
///
/// // Use the range information
/// println!("Block starts at byte {}", range.start_offset);
/// println!("Block contains {} records", range.block_records);
/// ```
#[derive(Debug, Clone, Copy)]
pub struct BlockRange {
    /// File offset where the block starts (in bytes, including headers)
    ///
    /// This is the absolute byte position in the file where this block begins,
    /// including the file header and block header.
    ///
    /// (8 bytes in serialized form)
    pub start_offset: u64,

    /// Length of the block data in bytes
    ///
    /// This is the size of the block data, not including the block header.
    /// For compressed blocks, this is the compressed size.
    ///
    /// (8 bytes in serialized form)
    pub len: u64,

    /// Number of records contained in this block
    ///
    /// (4 bytes in serialized form)
    pub block_records: u32,

    /// Cumulative number of records up to this block
    ///
    /// This allows efficient determination of which block contains a specific record
    /// by index without scanning through all previous blocks.
    ///
    /// **BREAKING CHANGE (v0.7.0+)**: Changed from u32 to u64 to support files
    /// with more than 4 billion records.
    ///
    /// (8 bytes in serialized form)
    pub cumulative_records: u64,

    /// Reserved bytes for future extensions
    pub reservation: [u8; 4],
}
impl BlockRange {
    /// Creates a new `BlockRange` with the specified parameters
    ///
    /// # Parameters
    ///
    /// * `start_offset` - The byte offset in the file where this block starts
    /// * `len` - The length of the block data in bytes
    /// * `block_records` - The number of records contained in this block
    /// * `cumulative_records` - The total number of records up to and including this block
    ///
    /// # Returns
    ///
    /// A new `BlockRange` instance with the specified parameters
    ///
    /// # Examples
    ///
    /// ```rust
    /// use binseq::vbq::BlockRange;
    ///
    /// // Create a new block range for a block starting at byte 1024
    /// let range = BlockRange::new(1024, 8192, 1000, 5000);
    /// ```
    #[must_use]
    pub fn new(start_offset: u64, len: u64, block_records: u32, cumulative_records: u64) -> Self {
        Self {
            start_offset,
            len,
            block_records,
            cumulative_records,
            reservation: INDEX_RESERVATION,
        }
    }

    /// Serializes the block range to a binary format and writes it to the provided writer
    ///
    /// This method serializes the `BlockRange` to a fixed-size 32-byte structure and
    /// writes it to the provided writer. The serialized format is:
    /// - Bytes 0-7: `start_offset` (u64, little endian)
    /// - Bytes 8-15: len (u64, little endian)
    /// - Bytes 16-19: `block_records` (u32, little endian)
    /// - Bytes 20-23: `cumulative_records` (u32, little endian)
    /// - Bytes 24-31: reservation (8 bytes)
    ///
    /// # Parameters
    ///
    /// * `writer` - The destination to write the serialized block range to
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If the block range was successfully written
    /// * `Err(_)` - If an error occurred during writing
    pub fn write_bytes<W: Write>(&self, writer: &mut W) -> Result<()> {
        let mut buf = [0; SIZE_BLOCK_RANGE];
        LittleEndian::write_u64(&mut buf[0..8], self.start_offset);
        LittleEndian::write_u64(&mut buf[8..16], self.len);
        LittleEndian::write_u32(&mut buf[16..20], self.block_records);
        LittleEndian::write_u64(&mut buf[20..28], self.cumulative_records);
        buf[28..].copy_from_slice(&self.reservation);
        writer.write_all(&buf)?;
        Ok(())
    }

    /// Deserializes a `BlockRange` from a fixed-size buffer
    ///
    /// This method deserializes a `BlockRange` from a 32-byte buffer in the format
    /// used by `write_bytes`. It's typically used when reading an index file.
    ///
    /// # Parameters
    ///
    /// * `buffer` - A fixed-size buffer containing a serialized `BlockRange`
    ///
    /// # Returns
    ///
    /// A new `BlockRange` with the values read from the buffer
    ///
    /// # Format
    ///
    /// The buffer is expected to contain:
    /// - Bytes 0-7: `start_offset` (u64, little endian)
    /// - Bytes 8-15: len (u64, little endian)
    /// - Bytes 16-19: `block_records` (u32, little endian)
    /// - Bytes 20-27: `cumulative_records` (u64, little endian)
    /// - Bytes 28-31: reservation (ignored, default value used)
    #[must_use]
    pub fn from_exact(buffer: &[u8; SIZE_BLOCK_RANGE]) -> Self {
        Self {
            start_offset: LittleEndian::read_u64(&buffer[0..8]),
            len: LittleEndian::read_u64(&buffer[8..16]),
            block_records: LittleEndian::read_u32(&buffer[16..20]),
            cumulative_records: LittleEndian::read_u64(&buffer[20..28]),
            reservation: INDEX_RESERVATION,
        }
    }

    /// Deserializes a `BlockRange` from a slice of bytes
    ///
    /// This is a convenience method that copies the first 32 bytes from the provided slice
    /// into a fixed-size buffer and then calls `from_exact`. It's useful when reading from
    /// a larger buffer that contains multiple serialized `BlockRange` instances.
    ///
    /// # Parameters
    ///
    /// * `buffer` - A slice containing at least 32 bytes with a serialized `BlockRange`
    ///
    /// # Returns
    ///
    /// A new `BlockRange` with the values read from the buffer
    ///
    /// # Panics
    ///
    /// This method will panic if the buffer is less than 32 bytes long.
    #[must_use]
    pub fn from_bytes(buffer: &[u8]) -> Self {
        let mut buf = [0; SIZE_BLOCK_RANGE];
        buf.copy_from_slice(buffer);
        Self::from_exact(&buf)
    }
}

/// Header for a VBQ index file
///
/// The `IndexHeader` contains metadata about an index file, including a magic number
/// for validation and the size of the indexed file. This allows verifying that an index
/// file matches its corresponding VBQ file.
///
/// The header has a fixed size of 32 bytes to ensure compatibility across versions.
#[derive(Debug, Clone, Copy)]
pub struct IndexHeader {
    /// Magic number to designate the index file ("VBQINDEX" in ASCII)
    ///
    /// This is used to verify that a file is indeed a VBQ index file.
    /// (8 bytes in serialized form)
    magic: u64,

    /// Total size of the indexed VBQ file in bytes
    ///
    /// This is used to verify that the index matches the file it references.
    /// (8 bytes in serialized form)
    bytes: u64,

    /// Reserved bytes for future extensions
    ///
    /// (16 bytes in serialized form)
    reserved: [u8; INDEX_HEADER_SIZE - 16],
}
impl IndexHeader {
    /// Creates a new index header for a VBQ file of the specified size
    ///
    /// # Parameters
    ///
    /// * `bytes` - The total size of the VBQ file being indexed, in bytes
    ///
    /// # Returns
    ///
    /// A new `IndexHeader` instance with the appropriate magic number and size
    pub fn new(bytes: u64) -> Self {
        Self {
            magic: INDEX_MAGIC,
            bytes,
            reserved: [42; INDEX_HEADER_SIZE - 16],
        }
    }
    /// Reads an index header from the provided reader
    ///
    /// This method reads 32 bytes from the provided reader and deserializes them
    /// into an `IndexHeader`. It validates the magic number to ensure that the file
    /// is indeed a VBQ index file.
    ///
    /// # Parameters
    ///
    /// * `reader` - The source from which to read the header
    ///
    /// # Returns
    ///
    /// * `Ok(Self)` - If the header was successfully read and has a valid magic number
    /// * `Err(_)` - If an error occurred during reading or the magic number is invalid
    ///
    /// # Format
    ///
    /// The header is expected to be 32 bytes with the following structure:
    /// - Bytes 0-7: magic number (u64, little endian, must be `INDEX_MAGIC`)
    /// - Bytes 8-15: file size in bytes (u64, little endian)
    /// - Bytes 16-31: reserved for future extensions
    pub fn from_reader<R: Read>(reader: &mut R) -> Result<Self> {
        let mut buffer = [0; INDEX_HEADER_SIZE];
        reader.read_exact(&mut buffer)?;
        let magic = LittleEndian::read_u64(&buffer[0..8]);
        let bytes = LittleEndian::read_u64(&buffer[8..16]);
        let Ok(reserved) = buffer[16..INDEX_HEADER_SIZE].try_into() else {
            return Err(IndexError::InvalidReservedBytes.into());
        };
        if magic != INDEX_MAGIC {
            return Err(IndexError::InvalidMagicNumber(magic).into());
        }
        Ok(Self {
            magic,
            bytes,
            reserved,
        })
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let mut buffer = [0; INDEX_HEADER_SIZE];
        buffer.copy_from_slice(&bytes[..INDEX_HEADER_SIZE]);
        Self::from_reader(&mut Cursor::new(buffer))
    }

    /// Serializes the index header to a binary format and writes it to the provided writer
    ///
    /// This method serializes the `IndexHeader` to a fixed-size 32-byte structure and
    /// writes it to the provided writer. This is typically used when saving an index to a file.
    ///
    /// # Parameters
    ///
    /// * `writer` - The destination to write the serialized header to
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If the header was successfully written
    /// * `Err(_)` - If an error occurred during writing
    ///
    /// # Format
    ///
    /// The header is serialized as:
    /// - Bytes 0-7: magic number (u64, little endian)
    /// - Bytes 8-15: file size in bytes (u64, little endian)
    /// - Bytes 16-31: reserved for future extensions
    pub fn write_bytes<W: Write>(&self, writer: &mut W) -> Result<()> {
        let mut buffer = [0; INDEX_HEADER_SIZE];
        LittleEndian::write_u64(&mut buffer[0..8], self.magic);
        LittleEndian::write_u64(&mut buffer[8..16], self.bytes);
        buffer[16..].copy_from_slice(&self.reserved);
        writer.write_all(&buffer)?;
        Ok(())
    }
}

/// Complete index for a VBQ file
///
/// A `BlockIndex` contains metadata about a VBQ file and all of its blocks,
/// enabling efficient random access and parallel processing. It consists of an
/// `IndexHeader` and a collection of `BlockRange` entries, one for each block in
/// the file.
///
/// The index can be created by scanning a VBQ file or loaded from a previously
/// created index file. Once loaded, it provides information about block locations,
/// sizes, and record counts.
///
/// # Examples
///
/// ```rust,no_run
/// use binseq::vbq::{BlockIndex, MmapReader};
/// use std::path::Path;
///
/// // Create an index from a VBQ file
/// let vbq_path = Path::new("example.vbq");
/// let index = BlockIndex::from_vbq(vbq_path).unwrap();
///
/// // Save the index for future use
/// let index_path = Path::new("example.vbq.vqi");
/// index.save_to_path(index_path).unwrap();
///
/// // Use the index with a reader for parallel processing
/// let reader = MmapReader::new(vbq_path).unwrap();
/// println!("File contains {} blocks", index.n_blocks());
/// ```
#[derive(Debug, Clone)]
pub struct BlockIndex {
    /// Header containing metadata about the indexed file
    pub(crate) header: IndexHeader,

    /// Collection of block ranges, one for each block in the file
    pub(crate) ranges: Vec<BlockRange>,
}
impl BlockIndex {
    /// Creates a new empty block index with the specified header
    ///
    /// # Parameters
    ///
    /// * `header` - The index header containing metadata about the indexed file
    ///
    /// # Returns
    ///
    /// A new empty `BlockIndex` instance
    #[must_use]
    pub fn new(header: IndexHeader) -> Self {
        Self {
            header,
            ranges: Vec::default(),
        }
    }
    /// Returns the number of blocks in the indexed file
    ///
    /// # Returns
    ///
    /// The number of blocks in the VBQ file described by this index
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::BlockIndex;
    /// use std::path::Path;
    ///
    /// let index = BlockIndex::from_path(Path::new("example.vbq.vqi")).unwrap();
    /// println!("The file contains {} blocks", index.n_blocks());
    /// ```
    #[must_use]
    pub fn n_blocks(&self) -> usize {
        self.ranges.len()
    }

    /// Writes the collection of `BlockRange` to a file
    /// Saves the index to a file
    ///
    /// This writes the index header and all block ranges to a file, which can be loaded
    /// later to avoid rescanning the VBQ file. The index is compressed to reduce
    /// storage space.
    ///
    /// # Parameters
    ///
    /// * `path` - The path where the index file should be saved
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If the index was successfully saved
    /// * `Err(_)` - If an error occurred during saving
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::BlockIndex;
    /// use std::path::Path;
    ///
    /// // Create an index from a VBQ file
    /// let index = BlockIndex::from_vbq(Path::new("example.vbq")).unwrap();
    ///
    /// // Save it for future use
    /// index.save_to_path(Path::new("example.vbq.vqi")).unwrap();
    /// ```
    pub fn save_to_path<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let mut writer = File::create(path).map(BufWriter::new)?;
        self.header.write_bytes(&mut writer)?;
        let mut writer = Encoder::new(writer, 3)?.auto_finish();
        self.write_range(&mut writer)?;
        writer.flush()?;
        Ok(())
    }

    /// Write the index to an output buffer
    pub fn write_bytes<W: Write>(&self, writer: &mut W) -> Result<()> {
        self.header.write_bytes(writer)?;
        let mut writer = Encoder::new(writer, 3)?.auto_finish();
        self.write_range(&mut writer)?;
        writer.flush()?;
        Ok(())
    }

    /// Write the collection of `BlockRange` to an output handle
    /// Writes all block ranges to the provided writer
    ///
    /// This method is used internally by `save_to_path` to write the block ranges
    /// to an index file. It can also be used to serialize an index to any destination
    /// that implements `Write`.
    ///
    /// # Parameters
    ///
    /// * `writer` - The destination to write the block ranges to
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If all block ranges were successfully written
    /// * `Err(_)` - If an error occurred during writing
    pub fn write_range<W: Write>(&self, writer: &mut W) -> Result<()> {
        self.ranges
            .iter()
            .filter(|range| range.block_records > 0)
            .try_for_each(|range| -> Result<()> { range.write_bytes(writer) })
    }

    /// Adds a block range to the index
    ///
    /// This method is used internally during index creation to add information
    /// about each block in the file. Blocks are typically added in order.
    ///
    /// # Parameters
    ///
    /// * `range` - The block range to add to the index
    fn add_range(&mut self, range: BlockRange) {
        self.ranges.push(range);
    }

    /// Creates a new index by scanning a VBQ file
    ///
    /// This method memory-maps the specified VBQ file and scans it block by block
    /// to create an index. The index can then be saved to a file for future use, enabling
    /// efficient random access without rescanning the file.
    ///
    /// # Parameters
    ///
    /// * `path` - Path to the VBQ file to index
    ///
    /// # Returns
    ///
    /// * `Ok(Self)` - A new `BlockIndex` containing information about all blocks in the file
    /// * `Err(_)` - If an error occurred during file opening, validation, or scanning
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::BlockIndex;
    /// use std::path::Path;
    ///
    /// // Create an index from a VBQ file
    /// let index = BlockIndex::from_vbq(Path::new("example.vbq")).unwrap();
    ///
    /// // Save the index for future use
    /// index.save_to_path(Path::new("example.vbq.vqi")).unwrap();
    ///
    /// // Get statistics about the file
    /// println!("File contains {} blocks", index.n_blocks());
    ///
    /// // Analyze the record distribution
    /// if let Some(last_range) = index.ranges().last() {
    ///     println!("Total records: {}", last_range.cumulative_records);
    ///     println!("Average records per block: {}",
    ///              last_range.cumulative_records as f64 / index.n_blocks() as f64);
    /// }
    /// ```
    ///
    /// # Notes
    ///
    /// This method uses memory mapping for efficiency, which allows the operating system
    /// to load only the needed portions of the file into memory as they are accessed.
    pub fn from_vbq<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { memmap2::Mmap::map(&file)? };
        let file_size = mmap.len();

        // Read header from mapped memory (unused but checks for validity)
        let _header = {
            let mut header_bytes = [0u8; SIZE_HEADER];
            header_bytes.copy_from_slice(&mmap[..SIZE_HEADER]);
            FileHeader::from_bytes(&header_bytes)?
        };

        // Initialize position after the header
        let mut pos = SIZE_HEADER;

        // Initialize the collection
        let index_header = IndexHeader::new(file_size as u64);
        let mut index = BlockIndex::new(index_header);

        // Find all block headers
        let mut record_total = 0;
        while pos < mmap.len() {
            let block_header = {
                let mut header_bytes = [0u8; SIZE_BLOCK_HEADER];
                header_bytes.copy_from_slice(&mmap[pos..pos + SIZE_BLOCK_HEADER]);
                BlockHeader::from_bytes(&header_bytes)?
            };
            index.add_range(BlockRange::new(
                pos as u64,
                block_header.size,
                block_header.records,
                record_total,
            ));
            pos += SIZE_BLOCK_HEADER + block_header.size as usize;
            record_total += u64::from(block_header.records);
        }

        Ok(index)
    }

    /// Reads an index from a path
    ///
    /// # Panics
    /// Panics if the path is not a valid UTF-8 string.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let Some(upstream_file) = path.as_ref().to_str().unwrap().strip_suffix(".vqi") else {
            return Err(IndexError::MissingUpstreamFile(
                path.as_ref().to_string_lossy().to_string(),
            )
            .into());
        };
        let upstream_handle = File::open(upstream_file)?;
        let mmap = unsafe { memmap2::Mmap::map(&upstream_handle)? };
        let file_size = mmap.len() as u64;

        let mut file_handle = File::open(path).map(BufReader::new)?;
        let index_header = IndexHeader::from_reader(&mut file_handle)?;
        if index_header.bytes != file_size {
            return Err(IndexError::ByteSizeMismatch(file_size, index_header.bytes).into());
        }
        let buffer = {
            let mut buffer = Vec::new();
            let mut decoder = Decoder::new(file_handle)?;
            decoder.read_to_end(&mut buffer)?;
            buffer
        };

        let mut ranges = Self::new(index_header);
        let mut pos = 0;
        while pos < buffer.len() {
            let bound = pos + SIZE_BLOCK_RANGE;
            let range = BlockRange::from_bytes(&buffer[pos..bound]);
            ranges.add_range(range);
            pos += SIZE_BLOCK_RANGE;
        }

        Ok(ranges)
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let index_header = IndexHeader::from_bytes(bytes)?;
        let buffer = {
            let mut buffer = Vec::new();
            let mut decoder = Decoder::new(Cursor::new(&bytes[INDEX_HEADER_SIZE..]))?;
            decoder.read_to_end(&mut buffer)?;
            buffer
        };

        let mut ranges = Self::new(index_header);
        let mut pos = 0;
        while pos < buffer.len() {
            let bound = pos + SIZE_BLOCK_RANGE;
            let range = BlockRange::from_bytes(&buffer[pos..bound]);
            ranges.add_range(range);
            pos += SIZE_BLOCK_RANGE;
        }

        Ok(ranges)
    }

    /// Get a reference to the internal ranges
    /// Returns a reference to the collection of block ranges
    ///
    /// This provides access to the metadata for all blocks in the indexed file,
    /// which can be used for operations like parallel processing or random access.
    ///
    /// # Returns
    ///
    /// A slice containing all `BlockRange` entries in this index
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// use binseq::vbq::BlockIndex;
    /// use std::path::Path;
    ///
    /// let index = BlockIndex::from_path(Path::new("example.vbq.vqi")).unwrap();
    ///
    /// // Examine the ranges to determine which blocks to process
    /// for (i, range) in index.ranges().iter().enumerate() {
    ///     println!("Block {}: {} records at offset {}",
    ///              i, range.block_records, range.start_offset);
    /// }
    /// ```
    #[must_use]
    pub fn ranges(&self) -> &[BlockRange] {
        &self.ranges
    }

    pub fn pprint(&self) {
        self.ranges.iter().for_each(|range| {
            println!(
                "{}\t{}\t{}\t{}",
                range.start_offset, range.len, range.block_records, range.cumulative_records
            );
        });
    }

    /// Returns the total number of records in the dataset
    #[must_use]
    pub fn num_records(&self) -> usize {
        self.ranges
            .iter()
            .next_back()
            .map(|r| (r.cumulative_records + u64::from(r.block_records)) as usize)
            .unwrap_or_default()
    }
}
