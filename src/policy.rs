//! Policies for handling invalid nucleotides during BQ/VBQ encoding.

use rand::Rng;

use crate::error::{Result, WriteError};

/// Seed for the RNG used by `RandomDraw`, for reproducibility across runs
pub const RNG_SEED: u64 = 42;

/// Policy for handling invalid nucleotides (anything other than A, C, G, T) during encoding
///
/// Defaults to `IgnoreSequence`, which skips offending sequences.
/// Only applies to BQ/VBQ; CBQ stores `N`s natively.
#[derive(Debug, Clone, Copy, Default)]
pub enum Policy {
    /// Skip sequences containing invalid nucleotides (default policy)
    #[default]
    IgnoreSequence,

    /// Fail with an error when invalid nucleotides are encountered
    BreakOnInvalid,

    /// Replace invalid nucleotides with randomly chosen valid nucleotides (A, C, G, or T)
    RandomDraw,

    /// Replace all invalid nucleotides with 'A'
    SetToA,

    /// Replace all invalid nucleotides with 'C'
    SetToC,

    /// Replace all invalid nucleotides with 'G'
    SetToG,

    /// Replace all invalid nucleotides with 'T'
    SetToT,
}
impl Policy {
    /// Replace invalid nucleotides with a specific nucleotide
    fn fill_with_known(sequence: &[u8], val: u8, ibuf: &mut Vec<u8>) {
        for &n in sequence {
            ibuf.push(match n {
                b'A' | b'C' | b'G' | b'T' => n,
                _ => val,
            });
        }
    }

    /// Replace invalid nucleotides with randomly chosen valid nucleotides
    fn fill_with_random<R: Rng>(sequence: &[u8], rng: &mut R, ibuf: &mut Vec<u8>) {
        for &n in sequence {
            ibuf.push(match n {
                b'A' | b'C' | b'G' | b'T' => n,
                _ => b"ACGT"[rng.random_range(0..4)],
            });
        }
    }

    /// Apply the policy to a sequence, writing the corrected result into `ibuf` (cleared first)
    ///
    /// Returns `Ok(true)` if the sequence should be encoded, `Ok(false)` if it should
    /// be skipped (`IgnoreSequence`), or an error (`BreakOnInvalid`).
    ///
    /// # Examples
    ///
    /// ```
    /// # use binseq::{Policy, Result};
    /// # fn main() -> Result<()> {
    /// let mut output = Vec::new();
    /// let mut rng = rand::rng();
    ///
    /// let should_process = Policy::SetToA.handle(b"ACGTNX", &mut output, &mut rng)?;
    ///
    /// assert!(should_process);
    /// assert_eq!(output, b"ACGTAA");
    /// # Ok(())
    /// # }
    /// ```
    pub fn handle<R: Rng>(&self, sequence: &[u8], ibuf: &mut Vec<u8>, rng: &mut R) -> Result<bool> {
        ibuf.clear();
        match self {
            Self::IgnoreSequence => Ok(false),
            Self::BreakOnInvalid => {
                let seq_str = std::str::from_utf8(sequence)?.to_string();
                Err(WriteError::InvalidNucleotideSequence(seq_str).into())
            }
            Self::RandomDraw => {
                Self::fill_with_random(sequence, rng, ibuf);
                Ok(true)
            }
            Self::SetToA => {
                Self::fill_with_known(sequence, b'A', ibuf);
                Ok(true)
            }
            Self::SetToC => {
                Self::fill_with_known(sequence, b'C', ibuf);
                Ok(true)
            }
            Self::SetToG => {
                Self::fill_with_known(sequence, b'G', ibuf);
                Ok(true)
            }
            Self::SetToT => {
                Self::fill_with_known(sequence, b'T', ibuf);
                Ok(true)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn run(policy: Policy, seq: &[u8]) -> (Result<bool>, Vec<u8>) {
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);
        let result = policy.handle(seq, &mut output, &mut rng);
        (result, output)
    }

    #[test]
    fn test_default_policy() {
        assert!(matches!(Policy::default(), Policy::IgnoreSequence));
    }

    #[test]
    fn test_ignore_sequence_policy() {
        let (result, output) = run(Policy::IgnoreSequence, b"ACGTNX");
        assert!(!result.unwrap()); // Should return false to skip this sequence
        assert!(output.is_empty());
    }

    #[test]
    fn test_break_on_invalid_policy() {
        let (result, _) = run(Policy::BreakOnInvalid, b"ACGTNX");
        assert!(matches!(
            result.unwrap_err(),
            crate::error::Error::WriteError(WriteError::InvalidNucleotideSequence(_))
        ));
    }

    #[test]
    fn test_set_to_policies() {
        // (policy, input, expected output)
        let cases: &[(Policy, &[u8], &[u8])] = &[
            (Policy::SetToA, b"ACGTNX", b"ACGTAA"),
            (Policy::SetToC, b"ACGTNX", b"ACGTCC"),
            (Policy::SetToG, b"ACGTNX", b"ACGTGG"),
            (Policy::SetToT, b"ACGTNX", b"ACGTTT"),
            // valid nucleotides remain unchanged
            (Policy::SetToA, b"ACGTACGT", b"ACGTACGT"),
            // empty sequence
            (Policy::SetToA, b"", b""),
            // all invalid
            (Policy::SetToG, b"NNNXXX", b"GGGGGG"),
            // lowercase and ambiguous codes are invalid
            (Policy::SetToA, b"acgt", b"AAAA"),
            (Policy::SetToC, b"AcGt", b"ACGC"),
            (Policy::SetToT, b"RYWSMK", b"TTTTTT"),
        ];
        for &(policy, seq, expected) in cases {
            let (result, output) = run(policy, seq);
            assert!(result.unwrap(), "{policy:?} on {seq:?}");
            assert_eq!(output, expected, "{policy:?} on {seq:?}");
        }
    }

    #[test]
    fn test_random_draw_policy() {
        let (result, output) = run(Policy::RandomDraw, b"ACGTNX");
        assert!(result.unwrap());
        assert_eq!(&output[0..4], b"ACGT"); // valid prefix unchanged
        assert!(output[4..].iter().all(|n| b"ACGT".contains(n)));
    }

    #[test]
    fn test_random_draw_deterministic_with_seed() {
        let (_, output1) = run(Policy::RandomDraw, b"NNNN");
        let (_, output2) = run(Policy::RandomDraw, b"NNNN");
        assert_eq!(output1, output2); // same seed, same output
    }

    #[test]
    fn test_buffer_cleared_before_processing() {
        let policy = Policy::SetToA;
        let mut output = vec![b'X', b'Y', b'Z']; // pre-filled buffer
        let mut rng = StdRng::seed_from_u64(RNG_SEED);
        policy.handle(b"ACGT", &mut output, &mut rng).unwrap();
        assert_eq!(output, b"ACGT");

        policy.handle(b"TT", &mut output, &mut rng).unwrap();
        assert_eq!(output, b"TT"); // cleared again on the next call
    }
}
