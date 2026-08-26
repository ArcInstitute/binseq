//! Shared nucleotide encoder used by the BQ and VBQ writers.

use bitnuc_deprec::BitSize;
use rand::{SeedableRng, rngs::SmallRng};

use crate::{Policy, RNG_SEED, error::Result};

/// Encodes nucleotide sequences into a compact multi-bit binary format.
///
/// The encoder maintains internal buffers to avoid repeated allocations and
/// handles invalid nucleotides according to a configurable [`Policy`].
#[derive(Clone)]
pub(crate) struct Encoder {
    /// Bitsize of the nucleotides
    bitsize: BitSize,

    /// Reusable buffers for encoded nucleotides
    sbuffer: Vec<u64>,
    xbuffer: Vec<u64>,

    /// Reusable buffers for invalid nucleotide sequences
    s_ibuf: Vec<u8>,
    x_ibuf: Vec<u8>,

    /// Invalid Nucleotide Policy
    policy: Policy,

    /// Random Number Generator (seeded with `RNG_SEED` for reproducibility)
    rng: SmallRng,
}

impl Encoder {
    /// Initialize a new encoder with the given policy.
    pub fn with_policy(bitsize: BitSize, policy: Policy) -> Self {
        Self {
            bitsize,
            policy,
            sbuffer: Vec::default(),
            xbuffer: Vec::default(),
            s_ibuf: Vec::default(),
            x_ibuf: Vec::default(),
            rng: SmallRng::seed_from_u64(RNG_SEED),
        }
    }

    /// Returns the invalid-nucleotide policy of this encoder.
    pub fn policy(&self) -> Policy {
        self.policy
    }

    /// Encodes a single sequence.
    ///
    /// Will return `None` if the sequence is invalid and the policy does not allow correction.
    pub fn encode_single(&mut self, primary: &[u8]) -> Result<Option<&[u64]>> {
        self.clear();
        if self.bitsize.encode(primary, &mut self.sbuffer).is_err() {
            self.clear();
            if self
                .policy
                .handle(primary, &mut self.s_ibuf, &mut self.rng)?
            {
                self.bitsize.encode(&self.s_ibuf, &mut self.sbuffer)?;
            } else {
                return Ok(None);
            }
        }
        Ok(Some(&self.sbuffer))
    }

    /// Encodes a pair of sequences.
    ///
    /// Will return `None` if either sequence is invalid and the policy does not allow correction.
    pub fn encode_paired(
        &mut self,
        primary: &[u8],
        extended: &[u8],
    ) -> Result<Option<(&[u64], &[u64])>> {
        self.clear();
        if self.bitsize.encode(primary, &mut self.sbuffer).is_err()
            || self.bitsize.encode(extended, &mut self.xbuffer).is_err()
        {
            self.clear();
            if self
                .policy
                .handle(primary, &mut self.s_ibuf, &mut self.rng)?
                && self
                    .policy
                    .handle(extended, &mut self.x_ibuf, &mut self.rng)?
            {
                self.bitsize.encode(&self.s_ibuf, &mut self.sbuffer)?;
                self.bitsize.encode(&self.x_ibuf, &mut self.xbuffer)?;
            } else {
                return Ok(None);
            }
        }
        Ok(Some((&self.sbuffer, &self.xbuffer)))
    }

    /// Clear all buffers and reset the encoder.
    pub fn clear(&mut self) {
        self.sbuffer.clear();
        self.xbuffer.clear();
        self.s_ibuf.clear();
        self.x_ibuf.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_single_invalid_ignored() {
        let mut encoder = Encoder::with_policy(BitSize::Two, Policy::IgnoreSequence);
        assert!(encoder.encode_single(b"ACGTNNNN").unwrap().is_none());
    }

    #[test]
    fn test_encode_single_invalid_corrected() {
        let mut encoder = Encoder::with_policy(BitSize::Two, Policy::SetToA);
        assert!(encoder.encode_single(b"ACGTNNNN").unwrap().is_some());
    }

    #[test]
    fn test_encode_paired_invalid_ignored() {
        let mut encoder = Encoder::with_policy(BitSize::Two, Policy::IgnoreSequence);
        assert!(
            encoder
                .encode_paired(b"ACGTNNNN", b"ACGTACGT")
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn test_encode_paired_invalid_corrected() {
        let mut encoder = Encoder::with_policy(BitSize::Two, Policy::SetToA);
        assert!(
            encoder
                .encode_paired(b"ACGTNNNN", b"NNNNACGT")
                .unwrap()
                .is_some()
        );
    }
}
