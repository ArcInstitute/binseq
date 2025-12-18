use std::fmt::Display;

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

    #[inline]
    pub fn is_paired(&self) -> bool {
        self.presence_flags & PRESENCE_PAIRED != 0
    }
    #[inline]
    pub fn has_qualities(&self) -> bool {
        self.presence_flags & PRESENCE_QUALITIES != 0
    }
    #[inline]
    pub fn has_headers(&self) -> bool {
        self.presence_flags & PRESENCE_HEADERS != 0
    }
    #[inline]
    pub fn has_flags(&self) -> bool {
        self.presence_flags & PRESENCE_FLAGS != 0
    }
}

impl Display for FileHeader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CBQ {{ version: {}, paired: {}, qualities: {}, headers: {}, flags: {}, block_size: {}, compression: {} }}",
            self.version,
            self.is_paired(),
            self.has_qualities(),
            self.has_headers(),
            self.has_flags(),
            self.block_size,
            self.compression_level,
        )
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

#[derive(Default)]
pub struct FileHeaderBuilder {
    compression_level: Option<usize>,
    block_size: Option<usize>,
    is_paired: Option<bool>,
    with_headers: Option<bool>,
    with_flags: Option<bool>,
    with_qualities: Option<bool>,
}

impl FileHeaderBuilder {
    pub fn with_compression_level(&mut self, compression_level: usize) -> &mut Self {
        self.compression_level = Some(compression_level);
        self
    }

    pub fn with_block_size(&mut self, block_size: usize) -> &mut Self {
        self.block_size = Some(block_size);
        self
    }

    pub fn is_paired(&mut self, is_paired: bool) -> &mut Self {
        self.is_paired = Some(is_paired);
        self
    }

    pub fn with_flags(&mut self, with_flags: bool) -> &mut Self {
        self.with_flags = Some(with_flags);
        self
    }

    pub fn with_headers(&mut self, with_headers: bool) -> &mut Self {
        self.with_headers = Some(with_headers);
        self
    }

    pub fn with_qualities(&mut self, with_qualities: bool) -> &mut Self {
        self.with_qualities = Some(with_qualities);
        self
    }

    pub fn build(&self) -> FileHeader {
        let mut header = FileHeader {
            magic: *FILE_MAGIC,
            version: FILE_VERSION,
            compression_level: self
                .compression_level
                .map_or(DEFAULT_COMPRESSION_LEVEL, |level| level as u64),
            block_size: self
                .block_size
                .map_or(DEFAULT_BLOCK_SIZE, |size| size as u64),
            presence_flags: 0,
            reserved: [0; 32],
        };

        // default to unpaired
        match self.is_paired {
            Some(true) => header.set_paired(),
            _ => {}
        }

        // default to using headers
        match self.with_headers {
            Some(false) => {}
            _ => header.set_headers(),
        };

        // default to not using flags
        match self.with_flags {
            Some(true) => header.set_flags(),
            _ => {}
        };

        // default to using qualities
        match self.with_qualities {
            Some(false) => {}
            _ => header.set_qualities(),
        };

        header
    }
}
