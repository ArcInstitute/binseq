//! Nucleotide sequence validation and correction policies
//!
//! This module provides policies for handling invalid nucleotides in sequences
//! during encoding operations. Different policies allow for ignoring, rejecting,
//! or correcting sequences with invalid nucleotides.

use rand::Rng;

use crate::error::{Result, WriteError};

/// A global seed for the random number generator used in randomized policies
///
/// This seed ensures reproducible behavior when using the `RandomDraw` policy
/// across different runs of the program.
pub const RNG_SEED: u64 = 42;

/// Policy for handling invalid nucleotide sequences during encoding
///
/// When encoding sequences into binary format, non-standard nucleotides (anything
/// other than A, C, G, or T) may be encountered. This enum defines different
/// strategies for handling such invalid nucleotides.
///
/// The default policy is `IgnoreSequence`, which skips sequences containing
/// invalid nucleotides.
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
    /// Helper method to replace invalid nucleotides with a specific nucleotide
    ///
    /// This internal method processes a sequence and replaces any non-standard
    /// nucleotides (anything other than A, C, G, or T) with the specified value.
    ///
    /// # Arguments
    ///
    /// * `sequence` - The input sequence to process
    /// * `val` - The replacement nucleotide (should be one of A, C, G, or T)
    /// * `ibuf` - The output buffer to store the processed sequence
    fn fill_with_known(sequence: &[u8], val: u8, ibuf: &mut Vec<u8>) {
        for &n in sequence {
            ibuf.push(match n {
                b'A' | b'C' | b'G' | b'T' => n,
                _ => val,
            });
        }
    }

    /// Helper method to replace invalid nucleotides with random valid nucleotides
    ///
    /// This internal method processes a sequence and replaces any non-standard
    /// nucleotides with randomly chosen valid nucleotides (A, C, G, or T).
    ///
    /// # Arguments
    ///
    /// * `sequence` - The input sequence to process
    /// * `rng` - The random number generator to use for selecting replacement nucleotides
    /// * `ibuf` - The output buffer to store the processed sequence
    ///
    /// # Type Parameters
    ///
    /// * `R` - A type that implements the `Rng` trait from the `rand` crate
    fn fill_with_random<R: Rng>(sequence: &[u8], rng: &mut R, ibuf: &mut Vec<u8>) {
        for &n in sequence {
            ibuf.push(match n {
                b'A' | b'C' | b'G' | b'T' => n,
                _ => match rng.random_range(0..4) {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    3 => b'T',
                    _ => unreachable!(),
                },
            });
        }
    }

    /// Process a sequence according to the selected policy for handling invalid nucleotides
    ///
    /// This method applies the policy to the given sequence, handling any invalid nucleotides
    /// according to the policy's rules. It first clears the input buffer to ensure that it is empty,
    /// then processes the sequence accordingly.
    ///
    /// # Arguments
    ///
    /// * `sequence` - The nucleotide sequence to be processed
    /// * `ibuf` - The buffer to store the processed sequence (will be cleared first)
    /// * `rng` - The random number generator (used only with the `RandomDraw` policy)
    ///
    /// # Returns
    ///
    /// * `Ok(true)` - If the sequence was processed and should be encoded
    /// * `Ok(false)` - If the sequence should be skipped (for `IgnoreSequence` policy)
    /// * `Err(Error)` - If an error occurred (for `BreakOnInvalid` policy when invalid nucleotides are found)
    ///
    /// # Type Parameters
    ///
    /// * `R` - A type that implements the `Rng` trait from the `rand` crate
    ///
    /// # Examples
    ///
    /// ```
    /// # use binseq::{Policy, Result};
    /// # use rand::thread_rng;
    /// # fn main() -> Result<()> {
    /// let policy = Policy::SetToA;
    /// let sequence = b"ACGTNX";
    /// let mut output = Vec::new();
    /// let mut rng = thread_rng();
    ///
    /// let should_process = policy.handle(sequence, &mut output, &mut rng)?;
    ///
    /// assert!(should_process);
    /// assert_eq!(output, b"ACGTAA");
    /// # Ok(())
    /// # }
    /// ```
    pub fn handle<R: Rng>(&self, sequence: &[u8], ibuf: &mut Vec<u8>, rng: &mut R) -> Result<bool> {
        // First clears the input buffer to ensure that it is empty.
        ibuf.clear();

        // Returns a boolean indicating whether the sequence should be processed further.
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

    // ==================== Basic Policy Tests ====================

    #[test]
    fn test_default_policy() {
        let policy = Policy::default();
        assert!(matches!(policy, Policy::IgnoreSequence));
    }

    #[test]
    fn test_ignore_sequence_policy() {
        let policy = Policy::IgnoreSequence;
        let sequence = b"ACGTNX";
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        let should_process = policy.handle(sequence, &mut output, &mut rng).unwrap();

        assert!(!should_process); // Should return false to skip this sequence
        assert!(output.is_empty()); // Output buffer should be empty
    }

    #[test]
    fn test_break_on_invalid_policy() {
        let policy = Policy::BreakOnInvalid;
        let sequence = b"ACGTNX";
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        let result = policy.handle(sequence, &mut output, &mut rng);

        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            crate::error::Error::WriteError(WriteError::InvalidNucleotideSequence(_))
        ));
    }

    #[test]
    fn test_break_on_invalid_with_valid_sequence() {
        let policy = Policy::BreakOnInvalid;
        let sequence = b"ACGT";
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        let result = policy.handle(sequence, &mut output, &mut rng);

        // Valid sequences should error because handle() doesn't validate for BreakOnInvalid
        // It only returns an error immediately
        assert!(result.is_err());
    }

    // ==================== Set-to-Specific-Nucleotide Tests ====================

    #[test]
    fn test_set_to_a_policy() {
        let policy = Policy::SetToA;
        let sequence = b"ACGTNX";
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        let should_process = policy.handle(sequence, &mut output, &mut rng).unwrap();

        assert!(should_process); // Should return true to process this sequence
        assert_eq!(output, b"ACGTAA"); // N and X should be replaced with A
    }

    #[test]
    fn test_set_to_c_policy() {
        let policy = Policy::SetToC;
        let sequence = b"ACGTNX";
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        let should_process = policy.handle(sequence, &mut output, &mut rng).unwrap();

        assert!(should_process);
        assert_eq!(output, b"ACGTCC"); // N and X should be replaced with C
    }

    #[test]
    fn test_set_to_g_policy() {
        let policy = Policy::SetToG;
        let sequence = b"ACGTNX";
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        let should_process = policy.handle(sequence, &mut output, &mut rng).unwrap();

        assert!(should_process);
        assert_eq!(output, b"ACGTGG"); // N and X should be replaced with G
    }

    #[test]
    fn test_set_to_t_policy() {
        let policy = Policy::SetToT;
        let sequence = b"ACGTNX";
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        let should_process = policy.handle(sequence, &mut output, &mut rng).unwrap();

        assert!(should_process);
        assert_eq!(output, b"ACGTTT"); // N and X should be replaced with T
    }

    #[test]
    fn test_all_valid_nucleotides_unchanged() {
        let policy = Policy::SetToA;
        let sequence = b"ACGTACGT";
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        let should_process = policy.handle(sequence, &mut output, &mut rng).unwrap();

        assert!(should_process);
        assert_eq!(output, b"ACGTACGT"); // All valid, should remain unchanged
    }

    // ==================== Random Draw Tests ====================

    #[test]
    fn test_random_draw_policy() {
        let policy = Policy::RandomDraw;
        let sequence = b"ACGTNX";
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        let should_process = policy.handle(sequence, &mut output, &mut rng).unwrap();

        assert!(should_process);
        assert_eq!(output.len(), 6); // Same length as input
        // First 4 nucleotides should be unchanged
        assert_eq!(&output[0..4], b"ACGT");
        // Last 2 should be valid nucleotides (A, C, G, or T)
        assert!(matches!(output[4], b'A' | b'C' | b'G' | b'T'));
        assert!(matches!(output[5], b'A' | b'C' | b'G' | b'T'));
    }

    #[test]
    fn test_random_draw_deterministic_with_seed() {
        let policy = Policy::RandomDraw;
        let sequence = b"NNNN";
        let mut output1 = Vec::new();
        let mut output2 = Vec::new();
        let mut rng1 = StdRng::seed_from_u64(RNG_SEED);
        let mut rng2 = StdRng::seed_from_u64(RNG_SEED);

        policy.handle(sequence, &mut output1, &mut rng1).unwrap();
        policy.handle(sequence, &mut output2, &mut rng2).unwrap();

        // Same seed should produce same output
        assert_eq!(output1, output2);
    }

    // ==================== Buffer Clearing Tests ====================

    #[test]
    fn test_buffer_cleared_before_processing() {
        let policy = Policy::SetToA;
        let sequence = b"ACGT";
        let mut output = vec![b'X', b'Y', b'Z']; // Pre-fill buffer
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        policy.handle(sequence, &mut output, &mut rng).unwrap();

        // Buffer should be cleared and only contain new data
        assert_eq!(output, b"ACGT");
    }

    #[test]
    fn test_multiple_calls_clear_buffer() {
        let policy = Policy::SetToA;
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        policy.handle(b"ACGT", &mut output, &mut rng).unwrap();
        assert_eq!(output, b"ACGT");

        policy.handle(b"TT", &mut output, &mut rng).unwrap();
        assert_eq!(output, b"TT"); // Should only contain second sequence
    }

    // ==================== Edge Case Tests ====================

    #[test]
    fn test_empty_sequence() {
        let policy = Policy::SetToA;
        let sequence = b"";
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        let should_process = policy.handle(sequence, &mut output, &mut rng).unwrap();

        assert!(should_process);
        assert!(output.is_empty());
    }

    #[test]
    fn test_all_invalid_nucleotides() {
        let policy = Policy::SetToG;
        let sequence = b"NNNXXX";
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        let should_process = policy.handle(sequence, &mut output, &mut rng).unwrap();

        assert!(should_process);
        assert_eq!(output, b"GGGGGG"); // All should be replaced with G
    }

    #[test]
    fn test_policy_clone() {
        let policy1 = Policy::SetToA;
        let policy2 = policy1;

        // Should be able to use both (tests Copy trait)
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        policy1.handle(b"NT", &mut output, &mut rng).unwrap();
        assert_eq!(output, b"AT");

        policy2.handle(b"NT", &mut output, &mut rng).unwrap();
        assert_eq!(output, b"AT");
    }

    #[test]
    fn test_policy_debug() {
        let policy = Policy::SetToA;
        let debug_str = format!("{:?}", policy);
        assert!(debug_str.contains("SetToA"));
    }

    // ==================== Various Invalid Character Tests ====================

    #[test]
    fn test_lowercase_nucleotides_treated_as_invalid() {
        let policy = Policy::SetToA;
        let sequence = b"acgt"; // lowercase
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        policy.handle(sequence, &mut output, &mut rng).unwrap();

        // Lowercase nucleotides should be treated as invalid
        assert_eq!(output, b"AAAA");
    }

    #[test]
    fn test_mixed_case_nucleotides() {
        let policy = Policy::SetToC;
        let sequence = b"AcGt";
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        policy.handle(sequence, &mut output, &mut rng).unwrap();

        assert_eq!(output, b"ACGC"); // Only uppercase are valid
    }

    #[test]
    fn test_ambiguous_nucleotide_codes() {
        let policy = Policy::SetToT;
        let sequence = b"RYWSMK"; // R, Y, W, S, M, K are ambiguous codes
        let mut output = Vec::new();
        let mut rng = StdRng::seed_from_u64(RNG_SEED);

        policy.handle(sequence, &mut output, &mut rng).unwrap();

        assert_eq!(output, b"TTTTTT"); // All ambiguous codes replaced with T
    }
}
