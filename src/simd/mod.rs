//! SIMD implementations for faster nucleotide processing
//!
//! Contains vectorized implementations for nucleotide encoding and decoding
//! using SIMD instructions on supported platforms.

#[cfg(target_arch = "aarch64")]
pub mod aarch64;

/// Check if SIMD is supported and not disabled by environment
pub fn is_simd_supported() -> bool {
    // Check for environment override
    if std::env::var("DISABLE_SIMD").is_ok() {
        return false;
    }

    #[cfg(target_arch = "aarch64")]
    {
        aarch64::is_neon_supported()
    }

    #[cfg(not(target_arch = "aarch64"))]
    {
        false
    }
}
