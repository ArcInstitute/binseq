use bytemuck::{Pod, Zeroable};
use zstd::stream::copy_encode;

use crate::{Result, error::CbqError};

use super::{BlockHeader, FileHeader, INDEX_MAGIC};

/// The header for a compressed index.
///
/// This is stored identically in memory and on disk.
#[derive(Debug, Clone, Copy, Zeroable, Pod)]
#[repr(C)]
pub struct IndexHeader {
    /// Magic number identifying the index format
    magic: [u8; 8],

    /// Number of bytes in the uncompressed index
    pub(crate) u_bytes: u64,

    /// Number of bytes in the compressed index
    pub(crate) z_bytes: u64,
}
impl IndexHeader {
    /// Creates a new index header
    #[must_use]
    pub fn new(u_bytes: u64, z_bytes: u64) -> Self {
        Self {
            magic: *INDEX_MAGIC,
            u_bytes,
            z_bytes,
        }
    }

    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        bytemuck::bytes_of(self)
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let header: Self = *bytemuck::from_bytes(bytes);
        if header.magic != *INDEX_MAGIC {
            return Err(CbqError::InvalidIndexHeaderMagic.into());
        }
        Ok(header)
    }
}

/// The footer for a compressed index.
///
/// This is stored identically in memory and on disk.
#[derive(Debug, Clone, Copy, Zeroable, Pod)]
#[repr(C)]
pub struct IndexFooter {
    /// Number of bytes in the compressed index
    pub(crate) bytes: u64,

    /// Magic number identifying the index format
    magic: [u8; 8],
}

impl IndexFooter {
    /// Creates a new index footer
    #[must_use]
    pub fn new(bytes: u64) -> Self {
        Self {
            bytes,
            magic: *INDEX_MAGIC,
        }
    }
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        bytemuck::bytes_of(self)
    }
    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let footer: Self = *bytemuck::from_bytes(bytes);
        if footer.magic != *INDEX_MAGIC {
            return Err(CbqError::InvalidIndexFooterMagic.into());
        }
        Ok(footer)
    }
}

/// An index of block ranges for quick lookups
#[derive(Clone)]
pub struct Index {
    ranges: Vec<BlockRange>,
}
impl Index {
    /// Builds the index from a list of block headers
    #[must_use]
    pub fn from_block_headers(block_headers: &[BlockHeader]) -> Self {
        let mut offset = size_of::<FileHeader>() as u64;
        let mut cumulative_records = 0;
        let mut ranges = Vec::default();
        for block_header in block_headers {
            let range = BlockRange::new(offset, cumulative_records + block_header.num_records);
            offset += (size_of::<BlockHeader>() + block_header.block_len()) as u64;
            cumulative_records += block_header.num_records;
            ranges.push(range);
        }
        Self { ranges }
    }

    /// Returns the byte representation of the index
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        bytemuck::cast_slice(&self.ranges)
    }

    /// Builds the index from a byte slice
    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let ranges = match bytemuck::try_cast_slice(bytes) {
            Ok(ranges) => ranges.to_vec(),
            Err(_) => return Err(CbqError::IndexCastingError.into()),
        };
        Ok(Self { ranges })
    }

    /// Returns the size of the index in bytes
    #[must_use]
    pub fn size(&self) -> u64 {
        self.as_bytes().len() as u64
    }

    /// Encodes the index into a ZSTD-compressed byte array
    pub fn encoded(&self) -> Result<Vec<u8>> {
        let mut encoded = Vec::default();
        copy_encode(self.as_bytes(), &mut encoded, 0)?;
        Ok(encoded)
    }

    /// Returns the number of records in the index
    #[must_use]
    pub fn num_records(&self) -> usize {
        self.ranges
            .last()
            .map_or(0, |range| range.cumulative_records as usize)
    }

    /// Returns the number of blocks in the index
    #[must_use]
    pub fn num_blocks(&self) -> usize {
        self.ranges.len()
    }

    #[must_use]
    pub fn iter_blocks(&self) -> BlockIter<'_> {
        BlockIter {
            index: self,
            pos: 0,
        }
    }

    #[must_use]
    pub fn average_block_size(&self) -> f64 {
        let mut block_iter = self.iter_blocks();
        let mut last_block = match block_iter.next() {
            Some(block) => block,
            None => return 0.0,
        };
        let mut total_size = 0.0;
        let mut count = 0;
        for block in block_iter {
            let last_block_size = block.offset - last_block.offset;
            total_size += last_block_size as f64;
            count += 1;
            last_block = block;
        }
        total_size / f64::from(count)
    }

    pub fn pprint(&self) {
        for block in self.iter_blocks() {
            println!("{block:?}");
        }
    }
}

pub struct BlockIter<'a> {
    index: &'a Index,
    pos: usize,
}
impl Iterator for BlockIter<'_> {
    type Item = BlockRange;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.index.num_blocks() {
            None
        } else {
            let block = self.index.ranges[self.pos];
            self.pos += 1;
            Some(block)
        }
    }
}

/// A struct representing a block range in a CBQ file and stored in the [`Index`](crate::cbq::Index)
///
/// This is stored identically in memory and on disk.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Zeroable, Pod, Default)]
#[repr(C)]
pub struct BlockRange {
    /// Byte offset of this block
    pub(crate) offset: u64,

    /// Number of records up to and including this block
    pub(crate) cumulative_records: u64,
}
impl BlockRange {
    #[must_use]
    pub fn new(offset: u64, cumulative_records: u64) -> Self {
        Self {
            offset,
            cumulative_records,
        }
    }
}
