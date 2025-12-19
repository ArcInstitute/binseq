use std::error::Error as StdError;

/// Custom Result type for binseq operations, wrapping the custom [`Error`] type
pub type Result<T> = std::result::Result<T, Error>;

/// The main error type for the binseq library, encompassing all possible error cases
/// that can occur during binary sequence operations.
#[derive(thiserror::Error, Debug)]
#[error(transparent)]
pub enum Error {
    /// Errors related to file and block headers
    #[error("Error processing header: {0}")]
    HeaderError(#[from] HeaderError),

    /// Errors related to the CBQ format
    #[error("Error processing CBQ: {0}")]
    CbqError(#[from] CbqError),

    /// Errors that occur during write operations
    #[error("Error writing file: {0}")]
    WriteError(#[from] WriteError),

    /// Errors that occur during read operations
    #[error("Error reading file: {0}")]
    ReadError(#[from] ReadError),

    /// Errors that occur during build operations
    #[error("Error building file: {0}")]
    BuilderError(#[from] BuilderError),

    /// Errors related to file indexing
    #[error("Error processing Index: {0}")]
    IndexError(#[from] IndexError),

    /// Standard I/O errors
    #[error("Error with IO: {0}")]
    IoError(#[from] std::io::Error),

    /// UTF-8 conversion errors
    #[error("Error with UTF8: {0}")]
    Utf8Error(#[from] std::str::Utf8Error),

    /// Errors related to missing extensions
    ExtensionError(#[from] ExtensionError),

    /// Errors from the bitnuc dependency for nucleotide encoding/decoding
    #[error("Bitnuc error: {0}")]
    BitnucError(#[from] bitnuc::Error),

    /// Conversion errors from anyhow errors
    #[error("Generic error: {0}")]
    AnyhowError(#[from] anyhow::Error),

    /// Generic errors for other unexpected situations
    #[error("Generic error: {0}")]
    GenericError(#[from] Box<dyn StdError + Send + Sync>),
}
impl Error {
    /// Checks if the error is an index mismatch error
    ///
    /// This is useful for determining if a file's index is out of sync with its content,
    /// which might require rebuilding the index.
    ///
    /// # Returns
    ///
    /// * `true` if the error is an `IndexError::ByteSizeMismatch`
    /// * `false` for all other error types
    #[must_use]
    pub fn is_index_mismatch(&self) -> bool {
        match self {
            Self::IndexError(err) => err.is_mismatch(),
            _ => false,
        }
    }
}

/// Errors specific to processing and validating binary sequence headers
#[derive(thiserror::Error, Debug)]
pub enum HeaderError {
    /// The magic number in the header does not match the expected value
    ///
    /// # Arguments
    /// * `u32` - The invalid magic number that was found
    #[error("Invalid magic number: {0}")]
    InvalidMagicNumber(u32),

    /// The format version in the header is not supported
    ///
    /// # Arguments
    /// * `u8` - The unsupported version number that was found
    #[error("Invalid format version: {0}")]
    InvalidFormatVersion(u8),

    /// The reserved bytes in the header contain unexpected values
    #[error("Invalid reserved bytes")]
    InvalidReservedBytes,

    /// The bits in the header contain unexpected values
    #[error("Invalid bit size found in header: {0} - expecting [2,4]")]
    InvalidBitSize(u8),

    /// The size of the data does not match what was specified in the header
    ///
    /// # Arguments
    /// * First `usize` - The actual number of bytes provided
    /// * Second `usize` - The expected number of bytes according to the header
    #[error("Invalid number of bytes provided: {0}. Expected: {1}")]
    InvalidSize(usize, usize),
}

/// Errors that can occur while reading binary sequence data
#[derive(thiserror::Error, Debug)]
pub enum ReadError {
    /// The file being read is not a regular file (e.g., it might be a directory or special file)
    #[error("File is not regular")]
    IncompatibleFile,

    /// The file appears to be truncated or corrupted
    ///
    /// # Arguments
    /// * `usize` - The byte position where the truncation was detected
    #[error(
        "Number of bytes in file does not match expectation - possibly truncated at byte pos {0}"
    )]
    FileTruncation(usize),

    /// Attempted to access a record index that is beyond the available range
    ///
    /// # Arguments
    /// * First `usize` - The requested record index
    /// * Second `usize` - The maximum available record index
    #[error("Requested record index ({0}) is out of record range ({1})")]
    OutOfRange(usize, usize),

    /// End of stream was reached while reading
    #[error("End of stream reached")]
    EndOfStream,

    /// A partial record was encountered at the end of a stream
    ///
    /// # Arguments
    /// * `usize` - The number of bytes read in the partial record
    #[error("Partial record at end of stream ({0} bytes)")]
    PartialRecord(usize),

    /// When a block header contains an invalid magic number
    ///
    /// The first parameter is the invalid magic number, the second is the position in the file
    #[error("Unexpected Block Magic Number found: {0} at position {1}")]
    InvalidBlockMagicNumber(u64, usize),

    /// When trying to read a block but reaching the end of the file unexpectedly
    ///
    /// The parameter is the position in the file where the read was attempted
    #[error("Unable to find an expected full block at position {0}")]
    UnexpectedEndOfFile(usize),

    /// When the file metadata doesn't match the expected VBINSEQ format
    #[error("Unexpected file metadata")]
    InvalidFileType,

    /// Missing the index end magic number
    #[error("Missing index end magic number")]
    MissingIndexEndMagic,
}

#[derive(thiserror::Error, Debug)]
pub enum BuilderError {
    #[error("Missing sequence length")]
    MissingSlen,
}

/// Errors that can occur while writing binary sequence data
#[derive(thiserror::Error, Debug)]
pub enum WriteError {
    /// The length of the sequence being written does not match what was specified in the header
    ///
    /// # Fields
    /// * `expected` - The sequence length specified in the header
    /// * `got` - The actual length of the sequence being written
    #[error("Sequence length ({got}) does not match the header ({expected})")]
    UnexpectedSequenceLength { expected: u32, got: usize },

    /// The sequence contains invalid nucleotide characters
    ///
    /// # Arguments
    /// * `String` - Description of the invalid nucleotides found
    #[error("Invalid nucleotides found in sequence: {0}")]
    InvalidNucleotideSequence(String),

    /// Attempted to write data without first setting up the header
    #[error("Missing header in writer builder")]
    MissingHeader,

    /// When trying to write data without quality scores but the header specifies they should be present
    #[error("Quality flag is set in header but trying to write without quality scores.")]
    QualityFlagSet,

    /// When trying to write data without a pair but the header specifies paired records
    #[error("Paired flag is set in header but trying to write without record pair.")]
    PairedFlagSet,

    /// When trying to write quality scores but the header specifies they are not present
    #[error("Quality flag not set in header but trying to write quality scores.")]
    QualityFlagNotSet,

    /// When trying to write paired data but the header doesn't specify paired records
    #[error("Paired flag not set in header but trying to write with record pair.")]
    PairedFlagNotSet,

    /// When trying to write data without headers but the header specifies they should be present
    #[error("Header flag is set in header but trying to write without headers.")]
    HeaderFlagSet,

    /// When a record is too large to fit in a block of the configured size
    ///
    /// The first parameter is the record size, the second is the maximum block size
    #[error(
        "Encountered a record with embedded size {0} but the maximum block size is {1}. Rerun with increased block size."
    )]
    RecordSizeExceedsMaximumBlockSize(usize, usize),
    /// When trying to ingest blocks with different sizes than expected
    ///
    /// The first parameter is the expected size, the second is the found size
    #[error(
        "Incompatible block sizes encountered in BlockWriter Ingest. Found ({1}) Expected ({0})"
    )]
    IncompatibleBlockSizes(usize, usize),

    /// When trying to ingest data with an incompatible header
    ///
    /// The first parameter is the expected header, the second is the found header
    #[error("Incompatible headers found in VBinseqWriter::ingest. Found ({1:?}) Expected ({0:?})")]
    IncompatibleHeaders(crate::vbq::VBinseqHeader, crate::vbq::VBinseqHeader),
}

/// Errors related to VBINSEQ file indexing
///
/// These errors occur when there are issues with the index of a VBINSEQ file,
/// such as corruption or mismatches with the underlying file.
#[derive(thiserror::Error, Debug)]
pub enum IndexError {
    /// When the magic number in the index doesn't match the expected value
    ///
    /// The parameter is the invalid magic number that was found
    #[error("Invalid magic number: {0}")]
    InvalidMagicNumber(u64),

    /// When the index references a file that doesn't exist
    ///
    /// The parameter is the missing file path
    #[error("Index missing upstream file path: {0}")]
    MissingUpstreamFile(String),

    /// When the size of the file doesn't match what the index expects
    ///
    /// The first parameter is the actual file size, the second is the expected size
    #[error("Mismatch in size between upstream size: {0} and expected index size {1}")]
    ByteSizeMismatch(u64, u64),

    /// Invalid reserved bytes in the index header
    #[error("Invalid reserved bytes in index header")]
    InvalidReservedBytes,
}
impl IndexError {
    /// Checks if this error indicates a mismatch between the index and file
    ///
    /// This is useful to determine if the index needs to be rebuilt.
    ///
    /// # Returns
    ///
    /// * `true` for `ByteSizeMismatch` errors
    /// * `true` for any other error type (this behavior is likely a bug and should be fixed)
    #[must_use]
    pub fn is_mismatch(&self) -> bool {
        matches!(self, Self::ByteSizeMismatch(_, _) | _) // Note: this appears to always return true regardless of error type
    }
}

#[derive(thiserror::Error, Debug)]
pub enum CbqError {
    #[error(
        "Record size ({record_size}) exceeds maximum block size ({max_block_size}) - Try increasing block size."
    )]
    ExceedsMaximumBlockSize {
        max_block_size: usize,
        record_size: usize,
    },

    #[error("Cannot ingest block of size {other_block_size} into block of size {self_block_size}")]
    CannotIngestBlock {
        self_block_size: usize,
        other_block_size: usize,
    },

    /// Attempting to write a record into a full block
    #[error(
        "Block(size: {block_size}) will be exceeded by record size {record_size}. Current size: {current_size}"
    )]
    BlockFull {
        current_size: usize,
        record_size: usize,
        block_size: usize,
    },

    #[error(
        "Cannot push record ({attribute}: {actual}) with block configuration ({attribute}: {expected})"
    )]
    ConfigurationMismatch {
        attribute: &'static str,
        expected: bool,
        actual: bool,
    },

    #[error("Invalid block header MAGIC found")]
    InvalidBlockHeaderMagic,

    #[error("Invalid file header MAGIC found")]
    InvalidFileHeaderMagic,

    #[error("Invalid index header MAGIC found")]
    InvalidIndexHeaderMagic,

    #[error("Invalid index footer MAGIC found")]
    InvalidIndexFooterMagic,

    #[error("Unable to cast bytes to Index - likely an alignment error")]
    IndexCastingError,

    #[error("SequenceRecordBuilder failed on build due to missing primary sequence (`s_seq`)")]
    MissingSequenceOnSequencingRecord,
}

#[derive(thiserror::Error, Debug)]
pub enum ExtensionError {
    /// When the extension is not supported
    #[error("Unsupported extension in path: {0}")]
    UnsupportedExtension(String),
}

/// Trait for converting arbitrary errors into `Error`
pub trait IntoBinseqError {
    fn into_binseq_error(self) -> Error;
}

// Implement conversion for Box<dyn Error>
impl<E> IntoBinseqError for E
where
    E: StdError + Send + Sync + 'static,
{
    fn into_binseq_error(self) -> Error {
        Error::GenericError(Box::new(self))
    }
}

mod testing {
    #[allow(unused)]
    use super::*;
    use thiserror::Error;

    #[allow(unused)]
    #[derive(Error, Debug)]
    pub enum MyError {
        #[error("Custom error: {0}")]
        CustomError(String),
    }

    #[test]
    fn test_into_binseq_error() {
        let my_error = MyError::CustomError(String::from("some error"));
        let binseq_error = my_error.into_binseq_error();
        assert!(matches!(binseq_error, Error::GenericError(_)));
    }
}
