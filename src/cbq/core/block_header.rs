use std::io;

use bytemuck::{Pod, Zeroable};

use crate::{IntoBinseqError, Result, error::CbqError};

use super::{BLOCK_MAGIC, ColumnarBlock};

#[derive(Copy, Clone, Pod, Zeroable, Debug, PartialEq, Eq, Hash)]
#[repr(C)]
pub struct BlockHeader {
    magic: [u8; 3],
    version: u8,
    padding: [u8; 4],

    // length of compressed columns
    pub(crate) len_z_seq_len: u64,
    pub(crate) len_z_header_len: u64,
    pub(crate) len_z_npos: u64,
    pub(crate) len_z_seq: u64,
    pub(crate) len_z_flags: u64,
    pub(crate) len_z_headers: u64,
    pub(crate) len_z_qual: u64,

    // full decoded length of the sequence block
    pub(crate) nuclen: u64,

    // length of uncompressed N-positions (Elias-Fano encoded)
    pub(crate) len_nef: u64,

    // number of records in the block
    pub num_records: u64,
}
impl BlockHeader {
    pub fn from_block(block: &ColumnarBlock) -> Self {
        Self {
            magic: *BLOCK_MAGIC,
            version: 1,
            padding: [42; 4],
            len_z_seq_len: block.z_seq_len.len() as u64,
            len_z_header_len: block.z_header_len.len() as u64,
            len_z_npos: block.z_npos.len() as u64,
            len_z_seq: block.z_seq.len() as u64,
            len_z_flags: block.z_flags.len() as u64,
            len_z_headers: block.z_headers.len() as u64,
            len_z_qual: block.z_qual.len() as u64,
            nuclen: block.nuclen as u64,
            len_nef: block.len_nef as u64,
            num_records: block.num_records as u64,
        }
    }

    /// Calculate the length of the block in bytes.
    #[allow(dead_code)]
    pub fn block_len(&self) -> usize {
        (self.len_z_seq_len
            + self.len_z_header_len
            + self.len_z_npos
            + self.len_z_seq
            + self.len_z_flags
            + self.len_z_headers
            + self.len_z_qual) as usize
    }

    pub fn as_bytes(&self) -> &[u8] {
        bytemuck::bytes_of(self)
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        let header: Self = *bytemuck::from_bytes(bytes);
        if header.magic != *BLOCK_MAGIC {
            return Err(CbqError::InvalidBlockHeaderMagic.into());
        }
        Ok(header)
    }

    pub fn write<W: io::Write>(&self, writer: &mut W) -> Result<()> {
        writer
            .write_all(self.as_bytes())
            .map_err(IntoBinseqError::into_binseq_error)
    }
}
