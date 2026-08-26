//! Fixed-size file header for BQ files.

use bitnuc_deprec::BitSize;
use std::io::{Read, Write};

use crate::{
    error::{HeaderError, Result},
    utils::read_u32_le,
};

/// The magic bytes at the start of a BQ file on disk.
pub const FILE_MAGIC: [u8; 4] = *b"BSEQ";
const MAGIC: u32 = u32::from_le_bytes(FILE_MAGIC);

/// Current format version
const FORMAT: u8 = 1;

/// Size of the header in bytes
pub const SIZE_HEADER: usize = 32;

/// Reserved bytes in the header (future use)
pub const RESERVED: [u8; 17] = [42; 17];

#[derive(Debug, Clone, Copy)]
pub struct FileHeaderBuilder {
    slen: Option<u32>,
    xlen: Option<u32>,
    bitsize: Option<BitSize>,
    flags: Option<bool>,
}
impl Default for FileHeaderBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl FileHeaderBuilder {
    #[must_use]
    pub fn new() -> Self {
        FileHeaderBuilder {
            slen: None,
            xlen: None,
            bitsize: None,
            flags: None,
        }
    }
    #[must_use]
    pub fn slen(mut self, slen: u32) -> Self {
        self.slen = Some(slen);
        self
    }
    #[must_use]
    pub fn xlen(mut self, xlen: u32) -> Self {
        self.xlen = Some(xlen);
        self
    }
    #[must_use]
    pub fn bitsize(mut self, bitsize: BitSize) -> Self {
        self.bitsize = Some(bitsize);
        self
    }
    #[must_use]
    pub fn flags(mut self, flags: bool) -> Self {
        self.flags = Some(flags);
        self
    }
    pub fn build(self) -> Result<FileHeader> {
        Ok(FileHeader {
            magic: MAGIC,
            format: FORMAT,
            slen: if let Some(slen) = self.slen {
                slen
            } else {
                return Err(HeaderError::MissingSequenceLength.into());
            },
            xlen: self.xlen.unwrap_or(0),
            bits: self.bitsize.unwrap_or_default(),
            flags: self.flags.unwrap_or(false),
            reserved: RESERVED,
        })
    }
}

/// Fixed 32-byte header for BQ files.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FileHeader {
    /// Magic number identifying the file format (4 bytes)
    pub magic: u32,

    /// Format version (1 byte)
    pub format: u8,

    /// Primary sequence length (4 bytes)
    pub slen: u32,

    /// Secondary sequence length (4 bytes)
    pub xlen: u32,

    /// Bits per nucleotide, 2 or 4 (1 byte)
    pub bits: BitSize,

    /// Whether all records carry a flag attribute (1 byte)
    pub flags: bool,

    /// Reserved for future use (17 bytes)
    pub reserved: [u8; 17],
}
impl FileHeader {
    /// Checks if the file is paired
    #[must_use]
    pub fn is_paired(&self) -> bool {
        self.xlen > 0
    }

    /// Parses and validates a header from a fixed-size byte array
    pub fn from_bytes(buffer: &[u8; SIZE_HEADER]) -> Result<Self> {
        let magic = read_u32_le(&buffer[0..4]);
        if magic != MAGIC {
            return Err(HeaderError::InvalidMagicNumber(magic).into());
        }
        let format = buffer[4];
        if format != FORMAT {
            return Err(HeaderError::InvalidFormatVersion(format).into());
        }
        let slen = read_u32_le(&buffer[5..9]);
        let xlen = read_u32_le(&buffer[9..13]);
        let bits = match buffer[13] {
            0 | 2 | 42 => BitSize::Two,
            4 => BitSize::Four,
            x => return Err(HeaderError::InvalidBitSize(x).into()),
        };
        let flags = buffer[14] != 0;
        let Ok(reserved) = buffer[15..32].try_into() else {
            return Err(HeaderError::InvalidReservedBytes.into());
        };
        Ok(Self {
            magic,
            format,
            slen,
            xlen,
            bits,
            flags,
            reserved,
        })
    }

    /// Parses a header from the first `SIZE_HEADER` bytes of a buffer
    pub fn from_buffer(buffer: &[u8]) -> Result<Self> {
        let mut bytes = [0u8; SIZE_HEADER];
        if buffer.len() < SIZE_HEADER {
            return Err(HeaderError::InvalidSize(buffer.len(), SIZE_HEADER).into());
        }
        bytes.copy_from_slice(&buffer[..SIZE_HEADER]);
        Self::from_bytes(&bytes)
    }

    /// Serializes the header and writes it to a writer
    pub fn write_bytes<W: Write>(&self, writer: &mut W) -> Result<()> {
        let mut buffer = [0u8; SIZE_HEADER];
        buffer[0..4].copy_from_slice(&self.magic.to_le_bytes());
        buffer[4] = self.format;
        buffer[5..9].copy_from_slice(&self.slen.to_le_bytes());
        buffer[9..13].copy_from_slice(&self.xlen.to_le_bytes());
        buffer[13] = self.bits.into();
        buffer[14] = self.flags.into();
        buffer[15..32].copy_from_slice(&self.reserved);
        writer.write_all(&buffer)?;
        Ok(())
    }

    /// Reads and parses a header from a reader
    pub fn from_reader<R: Read>(reader: &mut R) -> Result<Self> {
        let mut buffer = [0u8; SIZE_HEADER];
        reader.read_exact(&mut buffer)?;
        Self::from_bytes(&buffer)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ==================== FileHeaderBuilder Tests ====================

    #[test]
    fn test_builder_default() {
        let builder = FileHeaderBuilder::default();
        // Missing slen should fail to build
        assert!(builder.build().is_err());
    }

    #[test]
    fn test_builder_bitsize() {
        let header = FileHeaderBuilder::new()
            .slen(64)
            .bitsize(BitSize::Four)
            .build()
            .unwrap();
        assert_eq!(header.bits, BitSize::Four);
    }

    #[test]
    fn test_builder_missing_slen() {
        let result = FileHeaderBuilder::new().xlen(10).build();
        assert!(result.is_err());
    }

    // ==================== from_bytes Tests ====================

    #[test]
    fn test_from_bytes_invalid_format_version() {
        let header = FileHeaderBuilder::new().slen(32).build().unwrap();
        let mut buffer = [0u8; SIZE_HEADER];
        let mut cursor = std::io::Cursor::new(&mut buffer[..]);
        header.write_bytes(&mut cursor).unwrap();
        buffer[4] = 99; // corrupt format version
        let result = FileHeader::from_bytes(&buffer);
        assert!(result.is_err());
    }

    #[test]
    fn test_from_bytes_four_bit_size() {
        let header = FileHeaderBuilder::new()
            .bitsize(BitSize::Four)
            .slen(32)
            .build()
            .unwrap();
        let mut buffer = [0u8; SIZE_HEADER];
        let mut cursor = std::io::Cursor::new(&mut buffer[..]);
        header.write_bytes(&mut cursor).unwrap();
        let parsed = FileHeader::from_bytes(&buffer).unwrap();
        assert_eq!(parsed.bits, BitSize::Four);
    }

    #[test]
    fn test_from_bytes_invalid_bitsize() {
        let header = FileHeaderBuilder::new().slen(32).build().unwrap();
        let mut buffer = [0u8; SIZE_HEADER];
        let mut cursor = std::io::Cursor::new(&mut buffer[..]);
        header.write_bytes(&mut cursor).unwrap();
        buffer[13] = 99; // corrupt bitsize byte
        let result = FileHeader::from_bytes(&buffer);
        assert!(result.is_err());
    }

    #[test]
    fn test_from_bytes_invalid_magic() {
        let buffer = [0u8; SIZE_HEADER];
        let result = FileHeader::from_bytes(&buffer);
        assert!(result.is_err());
    }

    // ==================== from_buffer Tests ====================

    #[test]
    fn test_from_buffer_too_small() {
        let buffer = [0u8; 10];
        let result = FileHeader::from_buffer(&buffer);
        assert!(result.is_err());
    }

    #[test]
    fn test_from_buffer_valid() {
        let header = FileHeaderBuilder::new().slen(32).build().unwrap();
        let mut buffer = Vec::new();
        header.write_bytes(&mut buffer).unwrap();
        buffer.extend_from_slice(&[0u8; 16]); // trailing data beyond header
        let parsed = FileHeader::from_buffer(&buffer).unwrap();
        assert_eq!(parsed.slen, 32);
    }

    // ==================== from_reader Tests ====================

    #[test]
    fn test_from_reader_valid() {
        let header = FileHeaderBuilder::new()
            .slen(32)
            .xlen(16)
            .flags(true)
            .build()
            .unwrap();
        let mut buffer = Vec::new();
        header.write_bytes(&mut buffer).unwrap();
        let mut cursor = std::io::Cursor::new(buffer);
        let parsed = FileHeader::from_reader(&mut cursor).unwrap();
        assert_eq!(parsed, header);
    }

    #[test]
    fn test_from_reader_truncated() {
        let mut cursor = std::io::Cursor::new(vec![0u8; 5]);
        let result = FileHeader::from_reader(&mut cursor);
        assert!(result.is_err());
    }
}
