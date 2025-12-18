mod core;
mod read;
mod write;

pub use core::{BlockHeader, BlockRange, FileHeader, SequencingRecord};
pub use read::{MmapReader, Reader};
pub use write::ColumnarBlockWriter;

pub const FILE_MAGIC: &[u8; 7] = b"CBQFILE";
pub const BLOCK_MAGIC: &[u8; 3] = b"BLK";
pub const INDEX_MAGIC: &[u8; 8] = b"CBQINDEX";

pub const FILE_VERSION: u8 = 1;
pub const DEFAULT_BLOCK_SIZE: u64 = 1024 * 1024;
pub const DEFAULT_COMPRESSION_LEVEL: u64 = 0;
