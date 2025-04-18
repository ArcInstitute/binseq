//! Nucleotide encoding and decoding module with SIMD acceleration
//!
//! This module provides accelerated functions for encoding and decoding nucleotides
//! using SIMD instructions when available. It wraps the underlying `bitnuc` crate
//! and transparently uses SIMD operations on compatible hardware, falling back to
//! scalar operations when necessary.

use crate::error::Error;
use crate::simd::is_simd_supported;
use crate::Result;

/// Encodes nucleotides to 2-bit representation with potential SIMD acceleration
///
/// This function is a direct replacement for `bitnuc::encode`, but it will use
/// SIMD instructions when available for better performance.
///
/// # Arguments
///
/// * `input` - A slice of ASCII nucleotides (A, C, G, T)
/// * `output` - A mutable slice to store the encoded nucleotides
///
/// # Returns
///
/// * `Ok(())` - If encoding was successful
/// * `Err(Error)` - If invalid nucleotides were found in the input
///
/// # Example
///
/// ```
/// use binseq::nuc;
///
/// let sequence = b"ACGTACGT";
/// let mut output = vec![0u64; 2]; // Allocate enough space for the encoded sequence
/// nuc::encode(sequence, &mut output).unwrap();
/// ```
pub fn encode(input: &[u8], output: &mut Vec<u64>) -> Result<()> {
    if input.len() > 1_000
        && is_simd_supported()
        && crate::simd::aarch64::encode(input, output).is_ok()
    {
        return Ok(());
    }

    // Use scalar implementation
    bitnuc::encode(input, output).map_err(Error::from)
}

/// Decodes 2-bit nucleotides to ASCII with potential SIMD acceleration
///
/// This function is a direct replacement for `bitnuc::decode`, but it will use
/// SIMD instructions when available for better performance.
///
/// # Arguments
///
/// * `input` - A slice of encoded 2-bit nucleotides
/// * `len` - The number of nucleotides to decode
/// * `output` - A mutable vector to store the decoded ASCII nucleotides
///
/// # Returns
///
/// * `Ok(())` - If decoding was successful
/// * `Err(Error)` - If an error occurred during decoding
///
/// # Example
///
/// ```
/// use binseq::nuc;
///
/// let encoded = [0xe4]; // Encoded "ACGT" (binary: 11 10 01 00)
/// let mut output = Vec::new();
/// nuc::decode(&encoded, 4, &mut output).unwrap();
/// assert_eq!(output, b"ACGT");
/// ```
pub fn decode(input: &[u64], len: usize, output: &mut Vec<u8>) -> Result<()> {
    if input.len() > 1_000
        && is_simd_supported()
        && crate::simd::aarch64::decode(input, len, output).is_ok()
    {
        return Ok(());
    }

    bitnuc::decode(input, len, output).map_err(Error::from)
}

#[cfg(test)]
mod tests {
    use std::vec;

    use super::*;

    #[test]
    fn test_encode_decode_simple() {
        let sequence = b"ACGTACGT";
        let mut encoded = vec![0u64; 2]; // Allocate enough space for the encoded sequence
        let mut decoded = Vec::new();
        encode(sequence, &mut encoded).unwrap();
        decode(&encoded, sequence.len(), &mut decoded).unwrap();
        assert_eq!(sequence, decoded.as_slice());
    }

    #[test]
    fn test_encode_invalid() {
        let sequence = b"ACGTNACGT"; // Contains 'N'
        let encoded = [0u64; 1];
        let mut vec_encoded = vec![0; encoded.len()];
        vec_encoded.resize(encoded.len(), 0);
        let result = bitnuc::encode(sequence, &mut vec_encoded);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Invalid nucleotide base: 78"
        );
    }

    #[cfg(target_arch = "aarch64")]
    #[test]
    fn test_longer_sequence() {
        let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut encoded = [0u64; 2];
        unsafe {
            crate::simd::aarch64::encode_nucleotides_simd(sequence, &mut encoded).unwrap();
        }
        let mut decoded = Vec::new();
        decode(&encoded, sequence.len(), &mut decoded).unwrap();
        assert_eq!(sequence, decoded.as_slice());
    }

    #[cfg(target_arch = "aarch64")]
    #[test]
    fn test_differnt() {
        let sequence = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut encoded = Vec::new();
        let mut decoded = Vec::new();
        crate::simd::aarch64::encode(sequence, &mut encoded).unwrap();
        decode(&encoded, sequence.len(), &mut decoded).unwrap();
        assert_eq!(sequence, decoded.as_slice());
    }
}
