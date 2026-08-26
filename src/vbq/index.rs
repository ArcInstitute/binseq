//! # VBQ Index Format
//!
//! This module implements the embedded index format for VBQ files.
//!
//! ## Format Changes (v0.7.0+)
//!
//! **BREAKING CHANGE**: The VBQ index is now embedded at the end of VBQ files,
//! improving portability and eliminating the need to manage auxiliary files.
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
//! - Index is now embedded in VBQ files
//! - Cumulative record counts changed from `u32` to `u64`
//! - Support for files with more than 4 billion records

use std::io::{Cursor, Read, Write};

use zstd::{Decoder, Encoder};

use crate::utils::{read_u32_le, read_u64_le};

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
        buf[0..8].copy_from_slice(&self.start_offset.to_le_bytes());
        buf[8..16].copy_from_slice(&self.len.to_le_bytes());
        buf[16..20].copy_from_slice(&self.block_records.to_le_bytes());
        buf[20..28].copy_from_slice(&self.cumulative_records.to_le_bytes());
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
            start_offset: read_u64_le(&buffer[0..8]),
            len: read_u64_le(&buffer[8..16]),
            block_records: read_u32_le(&buffer[16..20]),
            cumulative_records: read_u64_le(&buffer[20..28]),
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
        let magic = read_u64_le(&buffer[0..8]);
        let bytes = read_u64_le(&buffer[8..16]);
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
        buffer[0..8].copy_from_slice(&self.magic.to_le_bytes());
        buffer[8..16].copy_from_slice(&self.bytes.to_le_bytes());
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
/// The index is embedded at the end of VBQ files and can be loaded using
/// `MmapReader::load_index()`. Once loaded, it provides information about block
/// locations, sizes, and record counts.
///
/// # Examples
///
/// ```rust,no_run
/// use binseq::vbq::MmapReader;
/// use std::path::Path;
///
/// // Load the embedded index from a VBQ file
/// let reader = MmapReader::new(Path::new("example.vbq")).unwrap();
/// let index = reader.load_index().unwrap();
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
    /// use binseq::vbq::{BlockIndex, MmapReader};
    /// use std::path::Path;
    ///
    /// let reader = MmapReader::new(Path::new("example.vbq")).unwrap();
    /// let index = reader.load_index().unwrap();
    /// println!("The file contains {} blocks", index.n_blocks());
    /// ```
    #[must_use]
    pub fn n_blocks(&self) -> usize {
        self.ranges.len()
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
    /// This method is used internally to write the block ranges to the embedded index.
    /// It can also be used to serialize an index to any destination that implements `Write`.
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
    /// use binseq::vbq::MmapReader;
    /// use std::path::Path;
    ///
    /// let reader = MmapReader::new(Path::new("example.vbq")).unwrap();
    /// let index = reader.load_index().unwrap();
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_num_records_empty_index() {
        let index = BlockIndex::new(IndexHeader::new(0));
        assert_eq!(index.num_records(), 0);
        assert_eq!(index.n_blocks(), 0);
    }

    // ==================== IndexHeader Tests ====================

    #[test]
    fn test_index_header_from_reader_invalid_magic() {
        let mut buffer = [0u8; INDEX_HEADER_SIZE];
        buffer[0..8].copy_from_slice(&0xDEAD_BEEF_u64.to_le_bytes());
        let result = IndexHeader::from_reader(&mut Cursor::new(buffer));
        assert!(result.is_err());
    }

    #[test]
    fn test_index_header_roundtrip() {
        let header = IndexHeader::new(12345);
        let mut buffer = Vec::new();
        header.write_bytes(&mut buffer).unwrap();
        let parsed = IndexHeader::from_reader(&mut Cursor::new(buffer)).unwrap();
        assert_eq!(parsed.bytes, 12345);
        assert_eq!(parsed.magic, INDEX_MAGIC);
    }

    // ==================== BlockIndex round-trip through write_bytes/from_bytes ====================

    #[test]
    fn test_block_index_write_and_from_bytes_roundtrip() {
        let mut index = BlockIndex::new(IndexHeader::new(999));
        index.add_range(BlockRange::new(32, 100, 5, 0));
        index.add_range(BlockRange::new(164, 200, 7, 5));

        let mut buffer = Vec::new();
        index.write_bytes(&mut buffer).unwrap();

        let parsed = BlockIndex::from_bytes(&buffer).unwrap();
        assert_eq!(parsed.n_blocks(), 2);
        assert_eq!(parsed.num_records(), 12);
        assert_eq!(parsed.ranges()[0].start_offset, 32);
        assert_eq!(parsed.ranges()[1].cumulative_records, 5);
    }

    #[test]
    fn test_block_index_write_range_skips_empty_blocks() {
        let mut index = BlockIndex::new(IndexHeader::new(0));
        index.add_range(BlockRange::new(32, 0, 0, 0));
        index.add_range(BlockRange::new(64, 50, 3, 0));

        let mut buffer = Vec::new();
        index.write_range(&mut buffer).unwrap();
        assert_eq!(buffer.len(), SIZE_BLOCK_RANGE);
    }
}
