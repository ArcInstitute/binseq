//! Embedded index for VBQ files (v0.7.0+).
//!
//! The index sits at the end of the file:
//!
//! ```text
//! [VBQ Data Blocks][Compressed Index][Index Size (u64)][INDEX_END_MAGIC (u64)]
//! ```
//!
//! The compressed section is ZSTD-compressed and holds an `IndexHeader`
//! (32 bytes) followed by one 32-byte `BlockRange` per data block.

use std::io::{Cursor, Read, Write};

use crate::utils::{read_u32_le, read_u64_le};
use zstd::{Decoder, Encoder};

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

/// Position, size, and record counts of a single block; 32 bytes serialized.
///
/// Stored in a `BlockIndex` to enable random access without scanning the file.
#[derive(Debug, Clone, Copy)]
pub struct BlockRange {
    /// Absolute file offset where the block (header included) starts (8 bytes)
    pub start_offset: u64,

    /// Block data length in bytes, excluding the block header; compressed size if compressed (8 bytes)
    pub len: u64,

    /// Number of records in this block (4 bytes)
    pub block_records: u32,

    /// Cumulative number of records before this block (8 bytes)
    pub cumulative_records: u64,

    /// Reserved bytes for future extensions (4 bytes)
    pub reservation: [u8; 4],
}
impl BlockRange {
    /// Creates a new `BlockRange`.
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

    /// Serializes the block range as 32 bytes (little-endian fields) to a writer.
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

    /// Deserializes a `BlockRange` from bytes (layout mirrors `write_bytes`).
    ///
    /// # Panics
    ///
    /// Panics if the buffer is less than 28 bytes long.
    #[must_use]
    pub fn from_bytes(buffer: &[u8]) -> Self {
        Self {
            start_offset: read_u64_le(&buffer[0..8]),
            len: read_u64_le(&buffer[8..16]),
            block_records: read_u32_le(&buffer[16..20]),
            cumulative_records: read_u64_le(&buffer[20..28]),
            reservation: INDEX_RESERVATION,
        }
    }
}

/// 32-byte header of the embedded index: `INDEX_MAGIC` (8 bytes),
/// indexed file size (8 bytes), and 16 reserved bytes.
#[derive(Debug, Clone, Copy)]
pub struct IndexHeader {
    /// Total size of the indexed VBQ file in bytes
    bytes: u64,
}
impl IndexHeader {
    /// Creates an index header for a VBQ file of the given size in bytes.
    pub fn new(bytes: u64) -> Self {
        Self { bytes }
    }
    /// Parses an index header from bytes, validating the magic number.
    pub fn from_bytes(buffer: &[u8]) -> Result<Self> {
        let magic = read_u64_le(&buffer[0..8]);
        if magic != INDEX_MAGIC {
            return Err(IndexError::InvalidMagicNumber(magic).into());
        }
        Ok(Self {
            bytes: read_u64_le(&buffer[8..16]),
        })
    }

    /// Serializes the index header as 32 bytes to a writer.
    pub fn write_bytes<W: Write>(self, writer: &mut W) -> Result<()> {
        let mut buffer = [42; INDEX_HEADER_SIZE];
        buffer[0..8].copy_from_slice(&INDEX_MAGIC.to_le_bytes());
        buffer[8..16].copy_from_slice(&self.bytes.to_le_bytes());
        writer.write_all(&buffer)?;
        Ok(())
    }
}

/// Complete embedded index of a VBQ file: an `IndexHeader` plus one
/// `BlockRange` per block.
///
/// Loaded from the end of a VBQ file with `MmapReader::load_index()`.
///
/// # Examples
///
/// ```rust,no_run
/// use binseq::vbq::MmapReader;
/// use std::path::Path;
///
/// let reader = MmapReader::new(Path::new("example.vbq")).unwrap();
/// let index = reader.load_index().unwrap();
/// println!("File contains {} blocks", index.n_blocks());
/// ```
#[derive(Debug, Clone)]
pub struct BlockIndex {
    /// Header containing metadata about the indexed file
    pub(crate) header: IndexHeader,

    /// Block ranges, one per block in the file
    pub(crate) ranges: Vec<BlockRange>,
}
impl BlockIndex {
    /// Creates a new empty block index with the given header.
    #[must_use]
    pub fn new(header: IndexHeader) -> Self {
        Self {
            header,
            ranges: Vec::default(),
        }
    }
    /// Returns the number of blocks in the indexed file.
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

    /// Writes all non-empty block ranges to a writer.
    pub fn write_range<W: Write>(&self, writer: &mut W) -> Result<()> {
        self.ranges
            .iter()
            .filter(|range| range.block_records > 0)
            .try_for_each(|range| -> Result<()> { range.write_bytes(writer) })
    }

    /// Adds a block range to the index.
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

    /// Returns the block ranges in this index.
    #[must_use]
    pub fn ranges(&self) -> &[BlockRange] {
        &self.ranges
    }

    /// Prints a tab-separated summary of each block range to stdout
    pub fn pprint(&self) {
        for range in &self.ranges {
            println!(
                "{}\t{}\t{}\t{}",
                range.start_offset, range.len, range.block_records, range.cumulative_records
            );
        }
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
    fn test_index_header_from_bytes_invalid_magic() {
        let mut buffer = [0u8; INDEX_HEADER_SIZE];
        buffer[0..8].copy_from_slice(&0xDEAD_BEEF_u64.to_le_bytes());
        let result = IndexHeader::from_bytes(&buffer);
        assert!(result.is_err());
    }

    #[test]
    fn test_index_header_roundtrip() {
        let header = IndexHeader::new(12345);
        let mut buffer = Vec::new();
        header.write_bytes(&mut buffer).unwrap();
        let parsed = IndexHeader::from_bytes(&buffer).unwrap();
        assert_eq!(parsed.bytes, 12345);
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
