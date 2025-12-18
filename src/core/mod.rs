mod block;
mod block_header;
mod header;
mod index;
mod sequencing_record;
pub(crate) mod utils;

pub use block::ColumnarBlock;
pub use block_header::BlockHeader;
pub use header::{FileHeader, FileHeaderBuilder};
pub use index::{BlockRange, Index, IndexFooter, IndexHeader};
pub use sequencing_record::{SequencingRecord, SequencingRecordBuilder};

use super::{
    BLOCK_MAGIC, DEFAULT_BLOCK_SIZE, DEFAULT_COMPRESSION_LEVEL, FILE_MAGIC, FILE_VERSION,
    INDEX_MAGIC,
};
