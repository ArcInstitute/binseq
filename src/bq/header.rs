//! Header module for the binseq library
//!
//! This module provides the header structure and functionality for binary sequence files.
//! The header contains metadata about the binary sequence data, including format version,
//! sequence length, and other information necessary for proper interpretation of the data.

use bitnuc_deprec::BitSize;
use byteorder::{ByteOrder, LittleEndian};
use std::io::{Read, Write};

use crate::error::{HeaderError, Result};

/// Current magic number: "BSEQ" in ASCII (in little-endian byte order)
///
/// This is used to identify binary sequence files and verify file integrity.
#[allow(clippy::unreadable_literal)]
const MAGIC: u32 = 0x51455342;

/// The magic bytes as they appear at the start of a BQ file on disk.
///
/// Used to identify BQ files by content rather than by file extension.
pub const FILE_MAGIC: [u8; 4] = MAGIC.to_le_bytes();

/// Current format version of the binary sequence file format
///
/// This version number allows for future format changes while maintaining backward compatibility.
const FORMAT: u8 = 1;

/// Size of the header in bytes
///
/// The header has a fixed size to ensure consistent reading and writing of binary sequence files.
pub const SIZE_HEADER: usize = 32;

/// Reserved bytes in the header
///
/// These bytes are reserved for future use and should be set to a consistent value.
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

/// Header structure for binary sequence files
///
/// The `FileHeader` contains metadata about the binary sequence data stored in a file,
/// including format information, sequence lengths, and space for future extensions.
///
/// The total size of this structure is 32 bytes, with a fixed layout to ensure
/// consistent reading and writing across different platforms.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FileHeader {
    /// Magic number to identify the file format
    ///
    /// 4 bytes
    pub magic: u32,

    /// Version of the file format
    ///
    /// 1 byte
    pub format: u8,

    /// Length of all sequences in the file
    ///
    /// 4 bytes
    pub slen: u32,

    /// Length of secondary sequences in the file
    ///
    /// 4 bytes
    pub xlen: u32,

    /// Number of bits per nucleotide (currently 2 or 4)
    ///
    /// 1 byte
    pub bits: BitSize,

    /// All records have a flag attribute
    ///
    /// 1 byte
    pub flags: bool,

    /// Reserve remaining bytes for future use
    ///
    /// 17 bytes
    pub reserved: [u8; 17],
}
impl FileHeader {
    /// Creates a new header with the specified sequence length
    ///
    /// This constructor initializes a standard header with the given sequence length,
    /// setting the magic number and format version to their default values.
    /// The extended sequence length (xlen) is set to 0.
    ///
    /// # Arguments
    ///
    /// * `bits` - The number of bits per nucleotide (currently 2 or 4)
    /// * `slen` - The length of sequences in the file
    /// * `flags` - The flags for the header
    ///
    /// # Returns
    ///
    /// A new `FileHeader` instance
    #[must_use]
    pub fn new(bits: BitSize, slen: u32, flags: bool) -> Self {
        Self {
            magic: MAGIC,
            format: FORMAT,
            slen,
            xlen: 0,
            bits,
            flags,
            reserved: RESERVED,
        }
    }

    /// Creates a new header with both primary and extended sequence lengths
    ///
    /// This constructor initializes a header for files that contain both primary
    /// and secondary sequence data, such as quality scores or annotations.
    ///
    /// # Arguments
    ///
    /// * `bits` - The number of bits per nucleotide (currently 2 or 4)
    /// * `slen` - The length of primary sequences in the file
    /// * `xlen` - The length of secondary/extended sequences in the file
    /// * `flags` - The flags for the header
    ///
    /// # Returns
    ///
    /// A new `FileHeader` instance with extended sequence information
    #[must_use]
    pub fn new_extended(bits: BitSize, slen: u32, xlen: u32, flags: bool) -> Self {
        Self {
            magic: MAGIC,
            format: FORMAT,
            slen,
            xlen,
            bits,
            flags,
            reserved: RESERVED,
        }
    }

    /// Sets the bitsize of the header
    pub fn set_bitsize(&mut self, bits: BitSize) {
        self.bits = bits;
    }

    /// Checks if the file is paired
    #[must_use]
    pub fn is_paired(&self) -> bool {
        self.xlen > 0
    }

    /// Parses a header from a fixed-size byte array
    ///
    /// This method validates the magic number and format version before constructing
    /// a header instance. If validation fails, appropriate errors are returned.
    ///
    /// # Arguments
    ///
    /// * `buffer` - A byte array of exactly `SIZE_HEADER` bytes containing the header data
    ///
    /// # Returns
    ///
    /// * `Ok(FileHeader)` - A valid header parsed from the buffer
    /// * `Err(Error)` - If the buffer contains invalid header data
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// * The magic number is incorrect
    /// * The format version is unsupported
    /// * The reserved bytes are invalid
    pub fn from_bytes(buffer: &[u8; SIZE_HEADER]) -> Result<Self> {
        let magic = LittleEndian::read_u32(&buffer[0..4]);
        if magic != MAGIC {
            return Err(HeaderError::InvalidMagicNumber(magic).into());
        }
        let format = buffer[4];
        if format != FORMAT {
            return Err(HeaderError::InvalidFormatVersion(format).into());
        }
        let slen = LittleEndian::read_u32(&buffer[5..9]);
        let xlen = LittleEndian::read_u32(&buffer[9..13]);
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

    /// Parses a header from an arbitrarily sized buffer
    ///
    /// This method extracts the header from the beginning of a buffer that may be larger
    /// than the header size. It checks that the buffer is at least as large as the header
    /// before attempting to parse it.
    ///
    /// # Arguments
    ///
    /// * `buffer` - A byte slice containing at least `SIZE_HEADER` bytes
    ///
    /// # Returns
    ///
    /// * `Ok(FileHeader)` - A valid header parsed from the buffer
    /// * `Err(Error)` - If the buffer is too small or contains invalid header data
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// * The buffer is smaller than `SIZE_HEADER`
    /// * The header data is invalid (see `from_bytes` for validation details)
    pub fn from_buffer(buffer: &[u8]) -> Result<Self> {
        let mut bytes = [0u8; SIZE_HEADER];
        if buffer.len() < SIZE_HEADER {
            return Err(HeaderError::InvalidSize(buffer.len(), SIZE_HEADER).into());
        }
        bytes.copy_from_slice(&buffer[..SIZE_HEADER]);
        Self::from_bytes(&bytes)
    }

    /// Writes the header to a writer
    ///
    /// This method serializes the header to its binary representation and writes it
    /// to the provided writer.
    ///
    /// # Arguments
    ///
    /// * `writer` - Any type that implements the `Write` trait
    ///
    /// # Returns
    ///
    /// * `Ok(())` - If the header was successfully written
    /// * `Err(Error)` - If writing to the writer failed
    ///
    /// # Errors
    ///
    /// Returns an error if writing to the writer fails (typically an I/O error).
    pub fn write_bytes<W: Write>(&self, writer: &mut W) -> Result<()> {
        let mut buffer = [0u8; SIZE_HEADER];
        LittleEndian::write_u32(&mut buffer[0..4], self.magic);
        buffer[4] = self.format;
        LittleEndian::write_u32(&mut buffer[5..9], self.slen);
        LittleEndian::write_u32(&mut buffer[9..13], self.xlen);
        buffer[13] = self.bits.into();
        buffer[14] = self.flags.into();
        buffer[15..32].copy_from_slice(&self.reserved);
        writer.write_all(&buffer)?;
        Ok(())
    }

    /// Reads a header from a reader
    ///
    /// This method reads exactly `SIZE_HEADER` bytes from the provided reader and
    /// parses them into a header structure.
    ///
    /// # Arguments
    ///
    /// * `reader` - Any type that implements the `Read` trait
    ///
    /// # Returns
    ///
    /// * `Ok(FileHeader)` - A valid header read from the reader
    /// * `Err(Error)` - If reading from the reader failed or the header data is invalid
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// * Reading from the reader fails (typically an I/O error)
    /// * The header data is invalid (see `from_bytes` for validation details)
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

    // ==================== FileHeader Constructor Tests ====================

    #[test]
    fn test_header_new() {
        let header = FileHeader::new(BitSize::Two, 100, true);
        assert_eq!(header.slen, 100);
        assert_eq!(header.xlen, 0);
        assert!(header.flags);
        assert!(!header.is_paired());
    }

    #[test]
    fn test_header_new_extended() {
        let header = FileHeader::new_extended(BitSize::Four, 100, 50, false);
        assert_eq!(header.slen, 100);
        assert_eq!(header.xlen, 50);
        assert_eq!(header.bits, BitSize::Four);
        assert!(header.is_paired());
    }

    #[test]
    fn test_set_bitsize() {
        let mut header = FileHeader::new(BitSize::Two, 100, false);
        header.set_bitsize(BitSize::Four);
        assert_eq!(header.bits, BitSize::Four);
    }

    // ==================== from_bytes Tests ====================

    #[test]
    fn test_from_bytes_invalid_format_version() {
        let header = FileHeader::new(BitSize::Two, 32, false);
        let mut buffer = [0u8; SIZE_HEADER];
        let mut cursor = std::io::Cursor::new(&mut buffer[..]);
        header.write_bytes(&mut cursor).unwrap();
        buffer[4] = 99; // corrupt format version
        let result = FileHeader::from_bytes(&buffer);
        assert!(result.is_err());
    }

    #[test]
    fn test_from_bytes_four_bit_size() {
        let header = FileHeader::new(BitSize::Four, 32, false);
        let mut buffer = [0u8; SIZE_HEADER];
        let mut cursor = std::io::Cursor::new(&mut buffer[..]);
        header.write_bytes(&mut cursor).unwrap();
        let parsed = FileHeader::from_bytes(&buffer).unwrap();
        assert_eq!(parsed.bits, BitSize::Four);
    }

    #[test]
    fn test_from_bytes_invalid_bitsize() {
        let header = FileHeader::new(BitSize::Two, 32, false);
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
        let header = FileHeader::new(BitSize::Two, 32, false);
        let mut buffer = Vec::new();
        header.write_bytes(&mut buffer).unwrap();
        buffer.extend_from_slice(&[0u8; 16]); // trailing data beyond header
        let parsed = FileHeader::from_buffer(&buffer).unwrap();
        assert_eq!(parsed.slen, 32);
    }

    // ==================== from_reader Tests ====================

    #[test]
    fn test_from_reader_valid() {
        let header = FileHeader::new_extended(BitSize::Two, 32, 16, true);
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
