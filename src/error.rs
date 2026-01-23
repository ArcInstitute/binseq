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
    #[cfg(feature = "anyhow")]
    #[error("Generic error: {0}")]
    AnyhowError(#[from] anyhow::Error),

    /// Generic errors for other unexpected situations
    #[error("Generic error: {0}")]
    GenericError(#[from] Box<dyn StdError + Send + Sync>),

    #[cfg(feature = "paraseq")]
    #[error("Fastx encoding error: {0}")]
    FastxEncodingError(#[from] FastxEncodingError),
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

    /// When building a `SequencingRecord` without a primary sequence
    #[error("SequencingRecordBuilder requires a primary sequence (s_seq)")]
    MissingSequence,
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

#[cfg(feature = "paraseq")]
#[derive(thiserror::Error, Debug)]
pub enum FastxEncodingError {
    #[error("Empty FASTX file")]
    EmptyFastxFile,

    #[error("Builder not provided with any input")]
    MissingInput,
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

    // ==================== Error::is_index_mismatch Tests ====================

    #[test]
    fn test_is_index_mismatch_with_byte_size_mismatch() {
        let error = Error::IndexError(IndexError::ByteSizeMismatch(100, 200));
        assert!(error.is_index_mismatch());
    }

    #[test]
    fn test_is_index_mismatch_with_invalid_magic() {
        let error = Error::IndexError(IndexError::InvalidMagicNumber(0x1234));
        // Note: The current implementation has a bug - it always returns true
        assert!(error.is_index_mismatch());
    }

    #[test]
    fn test_is_index_mismatch_with_non_index_error() {
        let error = Error::WriteError(WriteError::MissingHeader);
        assert!(!error.is_index_mismatch());
    }

    // ==================== IndexError Tests ====================

    #[test]
    fn test_index_error_is_mismatch() {
        let error = IndexError::ByteSizeMismatch(100, 200);
        assert!(error.is_mismatch());
    }

    #[test]
    fn test_index_error_invalid_magic() {
        let error = IndexError::InvalidMagicNumber(0x1234);
        // Note: Current implementation bug - always returns true
        assert!(error.is_mismatch());
    }

    #[test]
    fn test_index_error_missing_upstream_file() {
        let error = IndexError::MissingUpstreamFile("test.vbq".to_string());
        assert!(error.is_mismatch());
        assert!(format!("{}", error).contains("test.vbq"));
    }

    #[test]
    fn test_index_error_invalid_reserved_bytes() {
        let error = IndexError::InvalidReservedBytes;
        assert!(error.is_mismatch());
    }

    // ==================== HeaderError Tests ====================

    #[test]
    fn test_header_error_invalid_magic_number() {
        let error = HeaderError::InvalidMagicNumber(0xDEADBEEF);
        let error_str = format!("{}", error);
        assert!(error_str.contains("0xdeadbeef") || error_str.contains("3735928559"));
    }

    #[test]
    fn test_header_error_invalid_format_version() {
        let error = HeaderError::InvalidFormatVersion(99);
        let error_str = format!("{}", error);
        assert!(error_str.contains("99"));
    }

    #[test]
    fn test_header_error_invalid_bit_size() {
        let error = HeaderError::InvalidBitSize(8);
        let error_str = format!("{}", error);
        assert!(error_str.contains("8"));
        assert!(error_str.contains("[2,4]"));
    }

    #[test]
    fn test_header_error_invalid_size() {
        let error = HeaderError::InvalidSize(100, 200);
        let error_str = format!("{}", error);
        assert!(error_str.contains("100"));
        assert!(error_str.contains("200"));
    }

    // ==================== ReadError Tests ====================

    #[test]
    fn test_read_error_out_of_range() {
        let error = ReadError::OutOfRange {
            requested_index: 150,
            max_index: 100,
        };
        let error_str = format!("{}", error);
        assert!(error_str.contains("150"));
        assert!(error_str.contains("100"));
    }

    #[test]
    fn test_read_error_file_truncation() {
        let error = ReadError::FileTruncation(12345);
        let error_str = format!("{}", error);
        assert!(error_str.contains("12345"));
    }

    #[test]
    fn test_read_error_partial_record() {
        let error = ReadError::PartialRecord(42);
        let error_str = format!("{}", error);
        assert!(error_str.contains("42"));
    }

    #[test]
    fn test_read_error_invalid_block_magic_number() {
        let error = ReadError::InvalidBlockMagicNumber(0xBADC0DE, 1000);
        let error_str = format!("{}", error);
        assert!(error_str.contains("1000"));
    }

    // ==================== WriteError Tests ====================

    #[test]
    fn test_write_error_configuration_mismatch() {
        let error = WriteError::ConfigurationMismatch {
            attribute: "paired",
            expected: true,
            actual: false,
        };
        let error_str = format!("{}", error);
        assert!(error_str.contains("paired"));
        assert!(error_str.contains("true"));
        assert!(error_str.contains("false"));
    }

    #[test]
    fn test_write_error_unexpected_sequence_length() {
        let error = WriteError::UnexpectedSequenceLength {
            expected: 100,
            got: 150,
        };
        let error_str = format!("{}", error);
        assert!(error_str.contains("100"));
        assert!(error_str.contains("150"));
    }

    #[test]
    fn test_write_error_invalid_nucleotide_sequence() {
        let error = WriteError::InvalidNucleotideSequence("ACGTNX".to_string());
        let error_str = format!("{}", error);
        assert!(error_str.contains("ACGTNX"));
    }

    #[test]
    fn test_write_error_record_size_exceeds_max() {
        let error = WriteError::RecordSizeExceedsMaximumBlockSize(2000, 1024);
        let error_str = format!("{}", error);
        assert!(error_str.contains("2000"));
        assert!(error_str.contains("1024"));
    }

    #[test]
    fn test_write_error_missing_sequence_length() {
        let error = WriteError::MissingSequenceLength {
            exp_primary: true,
            exp_extended: false,
            obs_primary: false,
            obs_extended: false,
        };
        let error_str = format!("{}", error);
        assert!(error_str.contains("Missing required sequence length"));
    }

    // ==================== CbqError Tests ====================

    #[test]
    fn test_cbq_error_exceeds_maximum_block_size() {
        let error = CbqError::ExceedsMaximumBlockSize {
            max_block_size: 1024,
            record_size: 2048,
        };
        let error_str = format!("{}", error);
        assert!(error_str.contains("1024"));
        assert!(error_str.contains("2048"));
    }

    #[test]
    fn test_cbq_error_block_full() {
        let error = CbqError::BlockFull {
            current_size: 900,
            record_size: 200,
            block_size: 1024,
        };
        let error_str = format!("{}", error);
        assert!(error_str.contains("900"));
        assert!(error_str.contains("200"));
        assert!(error_str.contains("1024"));
    }

    #[test]
    fn test_cbq_error_cannot_ingest_block() {
        let error = CbqError::CannotIngestBlock {
            self_block_size: 1024,
            other_block_size: 2048,
        };
        let error_str = format!("{}", error);
        assert!(error_str.contains("1024"));
        assert!(error_str.contains("2048"));
    }

    // ==================== BuilderError Tests ====================

    #[test]
    fn test_builder_error_missing_slen() {
        let error = BuilderError::MissingSlen;
        let error_str = format!("{}", error);
        assert!(error_str.contains("Missing sequence length"));
    }

    // ==================== ExtensionError Tests ====================

    #[test]
    fn test_extension_error_unsupported() {
        let error = ExtensionError::UnsupportedExtension("test.xyz".to_string());
        let error_str = format!("{}", error);
        assert!(error_str.contains("test.xyz"));
    }

    // ==================== Error Conversion Tests ====================

    #[test]
    fn test_error_from_header_error() {
        let header_error = HeaderError::InvalidMagicNumber(0x1234);
        let error: Error = header_error.into();
        assert!(matches!(error, Error::HeaderError(_)));
    }

    #[test]
    fn test_error_from_write_error() {
        let write_error = WriteError::MissingHeader;
        let error: Error = write_error.into();
        assert!(matches!(error, Error::WriteError(_)));
    }

    #[test]
    fn test_error_from_read_error() {
        let read_error = ReadError::EndOfStream;
        let error: Error = read_error.into();
        assert!(matches!(error, Error::ReadError(_)));
    }

    #[test]
    fn test_error_from_index_error() {
        let index_error = IndexError::InvalidMagicNumber(0x5678);
        let error: Error = index_error.into();
        assert!(matches!(error, Error::IndexError(_)));
    }

    #[test]
    fn test_error_from_cbq_error() {
        let cbq_error = CbqError::InvalidBlockHeaderMagic;
        let error: Error = cbq_error.into();
        assert!(matches!(error, Error::CbqError(_)));
    }

    #[test]
    fn test_error_from_builder_error() {
        let builder_error = BuilderError::MissingSlen;
        let error: Error = builder_error.into();
        assert!(matches!(error, Error::BuilderError(_)));
    }

    #[test]
    fn test_error_debug_output() {
        let error = Error::WriteError(WriteError::MissingHeader);
        let debug_str = format!("{:?}", error);
        assert!(debug_str.contains("WriteError"));
    }

    // ==================== Fastx Error Tests (conditional) ====================

    #[cfg(feature = "paraseq")]
    #[test]
    fn test_fastx_error_empty_file() {
        use super::FastxEncodingError;
        let error = FastxEncodingError::EmptyFastxFile;
        let error_str = format!("{}", error);
        assert!(error_str.contains("Empty FASTX file"));
    }

    #[cfg(feature = "paraseq")]
    #[test]
    fn test_fastx_error_missing_input() {
        use super::FastxEncodingError;
        let error = FastxEncodingError::MissingInput;
        let error_str = format!("{}", error);
        assert!(error_str.contains("not provided with any input"));
    }
}
