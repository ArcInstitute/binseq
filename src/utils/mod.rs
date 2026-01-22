//! Utility modules for working with BINSEQ files

#[cfg(feature = "paraseq")]
pub mod fastx;

#[cfg(feature = "paraseq")]
pub use fastx::FastxEncoderBuilder;
