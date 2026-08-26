//! Utility modules for working with BINSEQ files

#[cfg(feature = "paraseq")]
pub mod fastx;

#[cfg(feature = "paraseq")]
pub use fastx::FastxEncoderBuilder;

/// Read a little-endian u64 from the start of a byte slice
pub(crate) fn read_u64_le(b: &[u8]) -> u64 {
    u64::from_le_bytes(b[..8].try_into().unwrap())
}

/// Read a little-endian u32 from the start of a byte slice
pub(crate) fn read_u32_le(b: &[u8]) -> u32 {
    u32::from_le_bytes(b[..4].try_into().unwrap())
}
