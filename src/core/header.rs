use anyhow::{Result, bail};
use bytemuck::{Pod, Zeroable};

use super::{DEFAULT_BLOCK_SIZE, DEFAULT_COMPRESSION_LEVEL, FILE_MAGIC, FILE_VERSION};

/// Records are paired
pub const PRESENCE_PAIRED: u64 = 1 << 0;

/// Records have quality scores
pub const PRESENCE_QUALITIES: u64 = 1 << 1;

/// Records have headers
pub const PRESENCE_HEADERS: u64 = 1 << 2;

/// Records have flags
pub const PRESENCE_FLAGS: u64 = 1 << 3;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Zeroable, Pod)]
#[repr(C)]
pub struct FileHeader {
    // File Type Metadata (8 bytes)
    /// File magic number
    magic: [u8; 7],
    /// File version number
    pub version: u8,

    // Data presence flags (8 bytes)
    /// A bitfield indicating which data fields are present in the file
    pub presence_flags: u64,

    // Configuration (16 bytes)
    /// compression level
    pub compression_level: u64,
    /// block size in bytes
    pub block_size: u64,

    /// Reserved for future use
    reserved: [u8; 32],
}
impl Default for FileHeader {
    fn default() -> Self {
        let mut header = Self {
            magic: *FILE_MAGIC,
            version: FILE_VERSION,
            presence_flags: 0,
            compression_level: DEFAULT_COMPRESSION_LEVEL,
            block_size: DEFAULT_BLOCK_SIZE,
            reserved: [0; 32],
        };
        header.set_headers();
        header.set_qualities();
        header
    }
}

/// Flag getters and setters
impl FileHeader {
    pub fn set_paired(&mut self) {
        self.presence_flags |= PRESENCE_PAIRED;
    }
    pub fn set_qualities(&mut self) {
        self.presence_flags |= PRESENCE_QUALITIES;
    }
    pub fn set_headers(&mut self) {
        self.presence_flags |= PRESENCE_HEADERS;
    }
    pub fn set_flags(&mut self) {
        self.presence_flags |= PRESENCE_FLAGS;
    }
    pub fn is_paired(&self) -> bool {
        self.presence_flags & PRESENCE_PAIRED != 0
    }
    pub fn has_qualities(&self) -> bool {
        self.presence_flags & PRESENCE_QUALITIES != 0
    }
    pub fn has_headers(&self) -> bool {
        self.presence_flags & PRESENCE_HEADERS != 0
    }
    pub fn has_flags(&self) -> bool {
        self.presence_flags & PRESENCE_FLAGS != 0
    }
}

impl FileHeader {
    pub fn as_bytes(&self) -> &[u8] {
        bytemuck::bytes_of(self)
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let header: Self = *bytemuck::from_bytes(bytes);
        if header.magic != *FILE_MAGIC {
            bail!("Invalid file magic")
        }
        Ok(header)
    }
}
