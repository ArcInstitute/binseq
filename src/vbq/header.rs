//! VBQ file and block header definitions.
//!
//! `FileHeader` opens the file; a `BlockHeader` precedes each record block.
//! Both are fixed-size (32 bytes) and carry magic numbers for validation.

use std::io::{Read, Write};

use bitnuc_deprec::BitSize;

use crate::{
    error::{HeaderError, ReadError, Result},
    utils::{read_u32_le, read_u64_le},
};

/// Magic bytes at the start of a VBQ file on disk
pub const FILE_MAGIC: [u8; 4] = *b"VSEQ";
const MAGIC: u32 = u32::from_le_bytes(FILE_MAGIC);

/// Block magic number: "BLOCKSEQ" in ASCII (0x5145534B434F4C42)
#[allow(clippy::unreadable_literal)]
const BLOCK_MAGIC: u64 = 0x5145534B434F4C42;

/// Current format version number
const FORMAT: u8 = 1;

/// Size of the file header in bytes
pub const SIZE_HEADER: usize = 32;

/// Size of each block header in bytes
pub const SIZE_BLOCK_HEADER: usize = 32;

/// Default virtual block size in bytes (128KB)
pub const BLOCK_SIZE: u64 = 128 * 1024;

/// Reserved bytes in the file header (placeholder value 42)
pub const RESERVED_BYTES: [u8; 13] = [42; 13];

/// Reserved bytes in block headers (placeholder value 42)
pub const RESERVED_BYTES_BLOCK: [u8; 12] = [42; 12];

#[derive(Default, Debug, Clone, Copy)]
pub struct FileHeaderBuilder {
    qual: Option<bool>,
    block: Option<u64>,
    compressed: Option<bool>,
    paired: Option<bool>,
    bitsize: Option<BitSize>,
    headers: Option<bool>,
    flags: Option<bool>,
}
impl FileHeaderBuilder {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
    #[must_use]
    pub fn qual(mut self, qual: bool) -> Self {
        self.qual = Some(qual);
        self
    }
    #[must_use]
    pub fn block(mut self, block: u64) -> Self {
        self.block = Some(block);
        self
    }
    #[must_use]
    pub fn compressed(mut self, compressed: bool) -> Self {
        self.compressed = Some(compressed);
        self
    }
    #[must_use]
    pub fn paired(mut self, paired: bool) -> Self {
        self.paired = Some(paired);
        self
    }
    #[must_use]
    pub fn bitsize(mut self, bitsize: BitSize) -> Self {
        self.bitsize = Some(bitsize);
        self
    }
    #[must_use]
    pub fn headers(mut self, headers: bool) -> Self {
        self.headers = Some(headers);
        self
    }
    #[must_use]
    pub fn flags(mut self, flags: bool) -> Self {
        self.flags = Some(flags);
        self
    }
    #[must_use]
    pub fn build(self) -> FileHeader {
        FileHeader {
            block: self.block.unwrap_or(BLOCK_SIZE),
            qual: self.qual.unwrap_or(false),
            compressed: self.compressed.unwrap_or(false),
            paired: self.paired.unwrap_or(false),
            bits: self.bitsize.unwrap_or_default(),
            headers: self.headers.unwrap_or(false),
            flags: self.flags.unwrap_or(false),
            ..FileHeader::default()
        }
    }
}

/// 32-byte header at the start of every VBQ file.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FileHeader {
    /// Magic number "VSEQ" (4 bytes)
    pub magic: u32,

    /// Format version (1 byte)
    pub format: u8,

    /// Virtual (uncompressed) block size in bytes (8 bytes)
    pub block: u64,

    /// Whether quality scores are included (1 byte)
    pub qual: bool,

    /// Whether blocks are ZSTD compressed (1 byte)
    pub compressed: bool,

    /// Whether records contain paired sequences (1 byte)
    pub paired: bool,

    /// Bits per nucleotide: 2-bit standard or 4-bit ambiguous (1 byte)
    pub bits: BitSize,

    /// Whether length-prefixed sequence headers are included (1 byte)
    pub headers: bool,

    /// Whether per-record flags are included (1 byte)
    pub flags: bool,

    /// Reserved bytes for future extensions (13 bytes)
    pub reserved: [u8; 13],
}
impl Default for FileHeader {
    /// Default block size (128KB), 2-bit encoding, all features disabled.
    fn default() -> Self {
        Self {
            magic: MAGIC,
            format: FORMAT,
            block: BLOCK_SIZE,
            qual: false,
            compressed: false,
            paired: false,
            headers: false,
            flags: false,
            bits: BitSize::default(),
            reserved: RESERVED_BYTES,
        }
    }
}
impl FileHeader {
    /// Parses a header from a 32-byte buffer, validating magic and version.
    pub fn from_bytes(buffer: &[u8; SIZE_HEADER]) -> Result<Self> {
        let magic = read_u32_le(&buffer[0..4]);
        if magic != MAGIC {
            return Err(HeaderError::InvalidMagicNumber(magic).into());
        }
        let format = buffer[4];
        if format != FORMAT {
            return Err(HeaderError::InvalidFormatVersion(format).into());
        }
        let block = read_u64_le(&buffer[5..13]);
        let qual = buffer[13] != 0;
        let compressed = buffer[14] != 0;
        let paired = buffer[15] != 0;
        let bits = match buffer[16] {
            0 | 2 | 42 => BitSize::Two,
            4 => BitSize::Four,
            x => return Err(HeaderError::InvalidBitSize(x).into()),
        };
        let headers = match buffer[17] {
            0 | 42 => false, // backwards compatibility
            _ => true,
        };
        let flags = buffer[18] != 0;
        let Ok(reserved) = buffer[19..32].try_into() else {
            return Err(HeaderError::InvalidReservedBytes.into());
        };
        Ok(Self {
            magic,
            format,
            block,
            qual,
            compressed,
            paired,
            bits,
            headers,
            flags,
            reserved,
        })
    }

    /// Serializes the header as 32 bytes to a writer.
    pub fn write_bytes<W: Write>(&self, writer: &mut W) -> Result<()> {
        let mut buffer = [0u8; SIZE_HEADER];
        buffer[0..4].copy_from_slice(&self.magic.to_le_bytes());
        buffer[4] = self.format;
        buffer[5..13].copy_from_slice(&self.block.to_le_bytes());
        buffer[13] = self.qual.into();
        buffer[14] = self.compressed.into();
        buffer[15] = self.paired.into();
        buffer[16] = self.bits.into();
        buffer[17] = self.headers.into();
        buffer[18] = self.flags.into();
        buffer[19..32].copy_from_slice(&self.reserved);
        writer.write_all(&buffer)?;
        Ok(())
    }

    /// Reads and parses a 32-byte header from a reader.
    pub fn from_reader<R: Read>(reader: &mut R) -> Result<Self> {
        let mut buffer = [0u8; SIZE_HEADER];
        reader.read_exact(&mut buffer)?;
        Self::from_bytes(&buffer)
    }

    #[must_use]
    pub fn is_paired(&self) -> bool {
        self.paired
    }
}

/// 32-byte header preceding each block in a VBQ file.
#[derive(Clone, Copy, Debug)]
pub struct BlockHeader {
    /// Magic number "BLOCKSEQ" (8 bytes)
    pub magic: u64,

    /// Actual on-disk block size in bytes; differs from the virtual size when compressed (8 bytes)
    pub size: u64,

    /// Number of records in this block (4 bytes)
    pub records: u32,

    /// Reserved bytes for future extensions (12 bytes)
    pub reserved: [u8; 12],
}
impl BlockHeader {
    /// Creates a new block header with the given on-disk size and record count.
    #[must_use]
    pub fn new(size: u64, records: u32) -> Self {
        Self {
            magic: BLOCK_MAGIC,
            size,
            records,
            reserved: RESERVED_BYTES_BLOCK,
        }
    }

    #[must_use]
    pub fn empty() -> Self {
        Self {
            magic: BLOCK_MAGIC,
            size: 0,
            records: 0,
            reserved: RESERVED_BYTES_BLOCK,
        }
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.size == 0 && self.records == 0
    }

    /// Serializes the block header as 32 bytes to a writer.
    pub fn write_bytes<W: Write>(&self, writer: &mut W) -> Result<()> {
        let mut buffer = [0u8; SIZE_BLOCK_HEADER];
        buffer[0..8].copy_from_slice(&self.magic.to_le_bytes());
        buffer[8..16].copy_from_slice(&self.size.to_le_bytes());
        buffer[16..20].copy_from_slice(&self.records.to_le_bytes());
        buffer[20..].copy_from_slice(&self.reserved);
        writer.write_all(&buffer)?;
        Ok(())
    }

    /// Parses a block header from a 32-byte buffer, validating the magic number.
    pub fn from_bytes(buffer: &[u8; SIZE_BLOCK_HEADER]) -> Result<Self> {
        let magic = read_u64_le(&buffer[0..8]);
        if magic != BLOCK_MAGIC {
            return Err(ReadError::InvalidBlockMagicNumber(magic, 0).into());
        }
        let size = read_u64_le(&buffer[8..16]);
        let records = read_u32_le(&buffer[16..20]);
        Ok(Self::new(size, records))
    }

    #[must_use]
    pub fn size_with_header(&self) -> usize {
        self.size as usize + SIZE_BLOCK_HEADER
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ==================== FileHeaderBuilder Tests ====================

    #[test]
    fn test_builder_block_and_bitsize() {
        let header = FileHeaderBuilder::new()
            .block(4096)
            .bitsize(BitSize::Four)
            .build();
        assert_eq!(header.block, 4096);
        assert_eq!(header.bits, BitSize::Four);
    }

    #[test]
    fn test_builder_defaults() {
        let header = FileHeaderBuilder::new().build();
        assert_eq!(header.block, BLOCK_SIZE);
        assert!(!header.qual);
        assert!(!header.compressed);
        assert!(!header.paired);
        assert!(!header.headers);
        assert!(!header.flags);
    }

    // ==================== FileHeader Constructor Tests ====================

    #[test]
    fn test_file_header_default() {
        let header = FileHeader::default();
        assert_eq!(header.block, BLOCK_SIZE);
        assert!(!header.is_paired());
    }

    // ==================== FileHeader from_bytes/from_reader Tests ====================

    #[test]
    fn test_file_header_roundtrip() {
        let header = FileHeaderBuilder::new()
            .qual(true)
            .paired(true)
            .headers(true)
            .flags(true)
            .build();
        let mut buffer = Vec::new();
        header.write_bytes(&mut buffer).unwrap();
        let mut cursor = std::io::Cursor::new(buffer);
        let parsed = FileHeader::from_reader(&mut cursor).unwrap();
        assert_eq!(parsed, header);
    }

    #[test]
    fn test_file_header_from_bytes_four_bit() {
        let header = FileHeaderBuilder::new().bitsize(BitSize::Four).build();
        let mut buffer = [0u8; SIZE_HEADER];
        {
            let mut cursor = std::io::Cursor::new(&mut buffer[..]);
            header.write_bytes(&mut cursor).unwrap();
        }
        let parsed = FileHeader::from_bytes(&buffer).unwrap();
        assert_eq!(parsed.bits, BitSize::Four);
    }

    #[test]
    fn test_file_header_from_bytes_invalid_magic() {
        let buffer = [0u8; SIZE_HEADER];
        let result = FileHeader::from_bytes(&buffer);
        assert!(result.is_err());
    }

    #[test]
    fn test_file_header_from_bytes_invalid_format() {
        let header = FileHeader::default();
        let mut buffer = [0u8; SIZE_HEADER];
        {
            let mut cursor = std::io::Cursor::new(&mut buffer[..]);
            header.write_bytes(&mut cursor).unwrap();
        }
        buffer[4] = 99;
        let result = FileHeader::from_bytes(&buffer);
        assert!(result.is_err());
    }

    #[test]
    fn test_file_header_from_bytes_invalid_bitsize() {
        let header = FileHeader::default();
        let mut buffer = [0u8; SIZE_HEADER];
        {
            let mut cursor = std::io::Cursor::new(&mut buffer[..]);
            header.write_bytes(&mut cursor).unwrap();
        }
        buffer[16] = 99;
        let result = FileHeader::from_bytes(&buffer);
        assert!(result.is_err());
    }

    #[test]
    fn test_file_header_from_reader_truncated() {
        let mut cursor = std::io::Cursor::new(vec![0u8; 5]);
        let result = FileHeader::from_reader(&mut cursor);
        assert!(result.is_err());
    }

    // ==================== BlockHeader Tests ====================

    #[test]
    fn test_block_header_from_bytes_invalid_magic() {
        let buffer = [0u8; SIZE_BLOCK_HEADER];
        let result = BlockHeader::from_bytes(&buffer);
        assert!(result.is_err());
    }

    #[test]
    fn test_block_header_roundtrip() {
        let header = BlockHeader::new(2048, 42);
        let mut buffer = Vec::new();
        header.write_bytes(&mut buffer).unwrap();
        let mut fixed = [0u8; SIZE_BLOCK_HEADER];
        fixed.copy_from_slice(&buffer);
        let parsed = BlockHeader::from_bytes(&fixed).unwrap();
        assert_eq!(parsed.size, 2048);
        assert_eq!(parsed.records, 42);
        assert!(!parsed.is_empty());
    }
}
