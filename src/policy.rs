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
