use std::error::Error as StdError;

/// Custom Result type for binseq operations, wrapping the custom [`Error`] type
pub type Result<T> = std::result::Result<T, Error>;

/// The main error type for the binseq library, encompassing all possible error cases
/// that can occur during binary sequence operations.
#[derive(thiserror::Error, Debug)]
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

    /// Errors related to file indexing
    #[error("Error processing Index: {0}")]
    IndexError(#[from] IndexError),

    /// Standard I/O errors
    #[error("Error with IO: {0}")]
    IoError(#[from] std::io::Error),

    /// UTF-8 conversion errors
    #[error("Error with UTF8: {0}")]
    Utf8Error(#[from] std::str::Utf8Error),

    /// Errors related to determining the BINSEQ format of a file
    #[error("Error determining BINSEQ format: {0}")]
    FormatError(#[from] FormatError),

    /// Errors from the deprecated bitnuc dependency for nucleotide encoding/decoding (bq/vbq)
    #[error("Bitnuc error: {0}")]
    BitnucError(#[from] bitnuc_deprec::Error),

    /// Errors from the bitnuc dependency for nucleotide encoding/decoding (cbq)
    #[error("Bitnuc error: {0}")]
    BitnucEncodingError(#[from] bitnuc::BitnucError),

    /// Conversion errors from anyhow errors
    #[cfg(feature = "anyhow")]
    #[error("Generic error: {0}")]
    AnyhowError(#[from] anyhow::Error),

    /// Generic errors for other unexpected situations
    #[error("Generic error: {0}")]
    GenericError(#[from] Box<dyn StdError + Send + Sync>),
}

/// Errors specific to processing and validating binary sequence headers
#[derive(thiserror::Error, Debug)]
pub enum HeaderError {
    /// The magic number in the header does not match the expected value
    #[error("Invalid magic number: {0}")]
    InvalidMagicNumber(u32),

    /// The format version in the header is not supported
    #[error("Invalid format version: {0}")]
    InvalidFormatVersion(u8),

    /// The reserved bytes in the header contain unexpected values
    #[error("Invalid reserved bytes")]
    InvalidReservedBytes,

    /// The bits in the header contain unexpected values
    #[error("Invalid bit size found in header: {0} - expecting [2,4]")]
    InvalidBitSize(u8),

    /// The size of the data does not match what was specified in the header
    #[error("Invalid number of bytes provided: {0}. Expected: {1}")]
    InvalidSize(usize, usize),

    /// A required sequence length was not provided when building the header
    #[error("Missing sequence length")]
    MissingSequenceLength,
}

/// Errors that can occur while reading binary sequence data
#[derive(thiserror::Error, Debug)]
pub enum ReadError {
    /// The file being read is not a regular file (e.g., it might be a directory or special file)
    #[error("File is not regular")]
    IncompatibleFile,

    /// The file appears to be truncated or corrupted
    #[error(
        "Number of bytes in file does not match expectation - possibly truncated at byte pos {0}"
    )]
    FileTruncation(usize),

    /// Attempted to access a record index that is beyond the available range
    #[error("Requested record index ({requested_index}) is out of record range ({max_index})")]
    OutOfRange {
        requested_index: usize,
        max_index: usize,
    },

    #[error("Invalid range specified: start ({start}) is greater than end ({end})")]
    InvalidRange { start: usize, end: usize },

    /// End of stream was reached while reading
    #[error("End of stream reached")]
    EndOfStream,

    /// A partial record was encountered at the end of a stream
    #[error("Partial record at end of stream ({0} bytes)")]
    PartialRecord(usize),

    /// A block header contains an invalid magic number
    #[error("Unexpected Block Magic Number found: {0} at position {1}")]
    InvalidBlockMagicNumber(u64, usize),

    /// Reached the end of the file unexpectedly while reading a block
    #[error("Unable to find an expected full block at position {0}")]
    UnexpectedEndOfFile(usize),

    /// The file metadata doesn't match the expected VBQ format
    #[error("Unexpected file metadata")]
    InvalidFileType,

    /// Missing the index end magic number
    #[error("Missing index end magic number")]
    MissingIndexEndMagic,
}

/// Errors that can occur while writing binary sequence data
#[derive(thiserror::Error, Debug)]
pub enum WriteError {
    /// Error between configuration of writer and incoming sequencing record
    #[error(
        "Cannot push record ({attribute}: {actual}) with writer configuration ({attribute}: {expected})"
    )]
    ConfigurationMismatch {
        attribute: &'static str,
        expected: bool,
        actual: bool,
    },

    #[error("Cannot ingest writer with incompatible formats")]
    FormatMismatch,

    #[error(
        "Missing required sequence length, expected (primary: {exp_primary}, extended: {exp_extended}), got (primary: {obs_primary}, extended: {obs_extended})"
    )]
    MissingSequenceLength {
        exp_primary: bool,
        exp_extended: bool,
        obs_primary: bool,
        obs_extended: bool,
    },

    /// The length of the sequence being written does not match the header
    #[error("Sequence length ({got}) does not match the header ({expected})")]
    UnexpectedSequenceLength { expected: u32, got: usize },

    /// The sequence contains invalid nucleotide characters
    #[error("Invalid nucleotides found in sequence: {0}")]
    InvalidNucleotideSequence(String),

    /// Attempted to write data without first setting up the header
    #[error("Missing header in writer builder")]
    MissingHeader,

    /// A record is too large to fit in a block of the configured size
    #[error(
        "Encountered a record with embedded size {0} but the maximum block size is {1}. Rerun with increased block size."
    )]
    RecordSizeExceedsMaximumBlockSize(usize, usize),

    /// Attempted to ingest blocks with a different size than expected
    #[error(
        "Incompatible block sizes encountered in BlockWriter Ingest. Found ({1}) Expected ({0})"
    )]
    IncompatibleBlockSizes(usize, usize),

    /// Attempted to ingest data with an incompatible header
    #[error("Incompatible headers found in vbq::Writer::ingest. Found ({1:?}) Expected ({0:?})")]
    IncompatibleHeaders(crate::vbq::FileHeader, crate::vbq::FileHeader),

    /// A `SequencingRecord` was built without a primary sequence
    #[error("SequencingRecordBuilder requires a primary sequence (s_seq)")]
    MissingSequence,

    /// FASTX encoding was attempted on an empty input file
    #[cfg(feature = "paraseq")]
    #[error("Empty FASTX file")]
    EmptyFastxFile,

    /// FASTX encoding was attempted without any input
    #[cfg(feature = "paraseq")]
    #[error("Builder not provided with any input")]
    MissingInput,
}

/// Errors related to VBQ file indexing
#[derive(thiserror::Error, Debug)]
pub enum IndexError {
    /// The magic number in the index doesn't match the expected value
    #[error("Invalid magic number: {0}")]
    InvalidMagicNumber(u64),

    /// Invalid reserved bytes in the index header
    #[error("Invalid reserved bytes in index header")]
    InvalidReservedBytes,
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
}

#[derive(thiserror::Error, Debug)]
pub enum FormatError {
    /// The BINSEQ format could not be determined from a file's magic bytes
    #[error("Unable to determine BINSEQ format from magic bytes in file: {0}")]
    UnrecognizedMagicBytes(String),
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

#[cfg(test)]
mod testing {
    use super::*;
    use thiserror::Error;

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

    #[test]
    fn test_error_from_sub_error() {
        let error: Error = HeaderError::InvalidMagicNumber(0x1234).into();
        assert!(matches!(error, Error::HeaderError(_)));
    }
}
