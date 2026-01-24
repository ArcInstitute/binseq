use crate::{BitSize, Result, error::WriteError};

/// A zero-copy record used to write sequences to binary sequence files.
///
/// This struct provides a unified API for writing records to all binseq formats
/// (BQ, VBQ, and CBQ). It uses borrowed references for zero-copy efficiency.
///
/// # Example
///
/// ```
/// use binseq::SequencingRecordBuilder;
///
/// let record = SequencingRecordBuilder::default()
///     .s_seq(b"ACGTACGT")
///     .s_qual(b"IIIIFFFF")
///     .s_header(b"seq_001")
///     .flag(42)
///     .build()
///     .unwrap();
/// ```
#[derive(Clone, Copy, Default)]
pub struct SequencingRecord<'a> {
    pub(crate) s_seq: &'a [u8],
    pub(crate) s_qual: Option<&'a [u8]>,
    pub(crate) s_header: Option<&'a [u8]>,
    pub(crate) x_seq: Option<&'a [u8]>,
    pub(crate) x_qual: Option<&'a [u8]>,
    pub(crate) x_header: Option<&'a [u8]>,
    pub(crate) flag: Option<u64>,
}

impl<'a> SequencingRecord<'a> {
    #[inline]
    #[must_use]
    pub fn new(
        s_seq: &'a [u8],
        s_qual: Option<&'a [u8]>,
        s_header: Option<&'a [u8]>,
        x_seq: Option<&'a [u8]>,
        x_qual: Option<&'a [u8]>,
        x_header: Option<&'a [u8]>,
        flag: Option<u64>,
    ) -> Self {
        Self {
            s_seq,
            s_qual,
            s_header,
            x_seq,
            x_qual,
            x_header,
            flag,
        }
    }

    /// Returns the primary sequence
    #[inline]
    #[must_use]
    pub fn s_seq(&self) -> &'a [u8] {
        self.s_seq
    }

    /// Returns the primary quality scores if present
    #[inline]
    #[must_use]
    pub fn s_qual(&self) -> Option<&'a [u8]> {
        self.s_qual
    }

    /// Returns the primary header if present
    #[inline]
    #[must_use]
    pub fn s_header(&self) -> Option<&'a [u8]> {
        self.s_header
    }

    /// Returns the extended/paired sequence if present
    #[inline]
    #[must_use]
    pub fn x_seq(&self) -> Option<&'a [u8]> {
        self.x_seq
    }

    /// Returns the extended quality scores if present
    #[inline]
    #[must_use]
    pub fn x_qual(&self) -> Option<&'a [u8]> {
        self.x_qual
    }

    /// Returns the extended header if present
    #[inline]
    #[must_use]
    pub fn x_header(&self) -> Option<&'a [u8]> {
        self.x_header
    }

    /// Returns the flag if present
    #[inline]
    #[must_use]
    pub fn flag(&self) -> Option<u64> {
        self.flag
    }

    /// Returns the configured size of this record for CBQ format.
    ///
    /// CBQ uses columnar storage so there are no per-record length prefixes.
    /// This calculates the size based on writer configuration, ignoring any
    /// extra data in the record that the writer won't use.
    #[inline]
    #[must_use]
    pub fn configured_size_cbq(
        &self,
        is_paired: bool,
        has_flags: bool,
        has_headers: bool,
        has_qualities: bool,
    ) -> usize {
        // CBQ uses 2-bit encoding: 4 nucleotides per byte, 32 per u64 word
        const NUCS_PER_WORD: usize = 32;

        let mut size = 0;

        // Sequence size (encoded into u64 words)
        let s_chunks = self.s_seq.len().div_ceil(NUCS_PER_WORD);
        size += s_chunks * 8;

        // Extended sequence (only if writer is configured for paired)
        if is_paired {
            let x_chunks = self.x_seq.map_or(0, |x| x.len().div_ceil(NUCS_PER_WORD));
            size += x_chunks * 8;
        }

        // Flag size (only if writer is configured for flags)
        if has_flags {
            size += 8; // u64
        }

        // Header size (only if writer is configured for headers)
        if has_headers {
            size += self.s_header.map_or(0, <[u8]>::len);
            if is_paired {
                size += self.x_header.map_or(0, <[u8]>::len);
            }
        }

        // Quality size (only if writer is configured for qualities)
        if has_qualities {
            size += self.s_qual.map_or(0, <[u8]>::len);
            if is_paired {
                size += self.x_qual.map_or(0, <[u8]>::len);
            }
        }

        size
    }

    /// Returns the configured size of this record for VBQ format.
    ///
    /// VBQ uses a row-based format with length prefixes for each field.
    /// This calculates the size based on writer configuration, ignoring any
    /// extra data in the record that the writer won't use.
    ///
    /// The VBQ record layout is:
    /// - Flag (8 bytes, if `has_flags`)
    /// - `s_len` (8 bytes)
    /// - `x_len` (8 bytes)
    /// - `s_seq` (encoded, rounded up to 8-byte words)
    /// - `s_qual` (raw bytes, if `has_qualities`)
    /// - `s_header_len` + `s_header` (8 + len bytes, if `has_headers` and `s_header` present)
    /// - `x_seq` (encoded, rounded up to 8-byte words, if paired)
    /// - `x_qual` (raw bytes, if `has_qualities` and paired)
    /// - `x_header_len` + `x_header` (8 + len bytes, if `has_headers` and `x_header` present)
    #[inline]
    #[must_use]
    pub fn configured_size_vbq(
        &self,
        is_paired: bool,
        has_flags: bool,
        has_headers: bool,
        has_qualities: bool,
        bitsize: BitSize,
    ) -> usize {
        // Calculate how many nucleotides fit per byte for the given bitsize
        let nucs_per_byte = if matches!(bitsize, BitSize::Two) {
            4
        } else {
            2
        };
        // VBQ packs sequences into u64 words
        let nucs_per_word = nucs_per_byte * 8;

        let mut size = 0;

        // Length prefixes: s_len and x_len (always present)
        size += 16; // 2 * u64

        // Flag (8 bytes, if has_flags)
        if has_flags {
            size += 8;
        }

        // Primary sequence (encoded into u64 words)
        let s_chunks = self.s_seq.len().div_ceil(nucs_per_word);
        size += s_chunks * 8;

        // Extended sequence (only if writer is configured for paired)
        if is_paired {
            let x_chunks = self.x_seq.map_or(0, |x| x.len().div_ceil(nucs_per_word));
            size += x_chunks * 8;
        }

        // Quality scores (raw bytes, only if writer configured for qualities)
        if has_qualities {
            size += self.s_qual.map_or(0, <[u8]>::len);
            if is_paired {
                size += self.x_qual.map_or(0, <[u8]>::len);
            }
        }

        // Headers (length prefix + raw bytes, only if writer configured for headers)
        if has_headers {
            if let Some(h) = self.s_header {
                size += 8 + h.len(); // length prefix + header bytes
            }
            if is_paired && let Some(h) = self.x_header {
                size += 8 + h.len(); // length prefix + header bytes
            }
        }

        size
    }

    #[inline]
    #[must_use]
    pub fn is_paired(&self) -> bool {
        self.x_seq.is_some()
    }

    #[inline]
    #[must_use]
    pub fn has_flags(&self) -> bool {
        self.flag.is_some()
    }

    #[inline]
    #[must_use]
    pub fn has_headers(&self) -> bool {
        self.s_header.is_some() || self.x_header.is_some()
    }

    #[inline]
    #[must_use]
    pub fn has_qualities(&self) -> bool {
        self.s_qual.is_some() || self.x_qual.is_some()
    }
}

/// A convenience builder struct for creating a [`SequencingRecord`]
///
/// # Example
///
/// ```
/// use binseq::SequencingRecordBuilder;
///
/// // Build a simple unpaired record
/// let record = SequencingRecordBuilder::default()
///     .s_seq(b"ACGTACGT")
///     .build()
///     .unwrap();
///
/// // Build a paired record with quality scores
/// let paired = SequencingRecordBuilder::default()
///     .s_seq(b"ACGTACGT")
///     .s_qual(b"IIIIFFFF")
///     .x_seq(b"TGCATGCA")
///     .x_qual(b"FFFFHHHH")
///     .flag(1)
///     .build()
///     .unwrap();
/// ```
#[derive(Default)]
pub struct SequencingRecordBuilder<'a> {
    s_seq: Option<&'a [u8]>,
    s_qual: Option<&'a [u8]>,
    s_header: Option<&'a [u8]>,
    x_seq: Option<&'a [u8]>,
    x_qual: Option<&'a [u8]>,
    x_header: Option<&'a [u8]>,
    flag: Option<u64>,
}

impl<'a> SequencingRecordBuilder<'a> {
    /// Sets the primary sequence (required)
    #[must_use]
    pub fn s_seq(mut self, s_seq: &'a [u8]) -> Self {
        self.s_seq = Some(s_seq);
        self
    }

    /// Sets the primary quality scores
    #[must_use]
    pub fn s_qual(mut self, s_qual: &'a [u8]) -> Self {
        self.s_qual = Some(s_qual);
        self
    }

    /// Sets the primary quality scores from an Option
    #[must_use]
    pub fn opt_s_qual(mut self, s_qual: Option<&'a [u8]>) -> Self {
        self.s_qual = s_qual;
        self
    }

    /// Sets the primary header
    #[must_use]
    pub fn s_header(mut self, s_header: &'a [u8]) -> Self {
        self.s_header = Some(s_header);
        self
    }

    /// Sets the primary header from an Option
    #[must_use]
    pub fn opt_s_header(mut self, s_header: Option<&'a [u8]>) -> Self {
        self.s_header = s_header;
        self
    }

    /// Sets the extended/paired sequence
    #[must_use]
    pub fn x_seq(mut self, x_seq: &'a [u8]) -> Self {
        self.x_seq = Some(x_seq);
        self
    }

    /// Sets the extended/paired sequence from an Option
    #[must_use]
    pub fn opt_x_seq(mut self, x_seq: Option<&'a [u8]>) -> Self {
        self.x_seq = x_seq;
        self
    }

    /// Sets the extended quality scores
    #[must_use]
    pub fn x_qual(mut self, x_qual: &'a [u8]) -> Self {
        self.x_qual = Some(x_qual);
        self
    }

    /// Sets the extended quality scores from an Option
    #[must_use]
    pub fn opt_x_qual(mut self, x_qual: Option<&'a [u8]>) -> Self {
        self.x_qual = x_qual;
        self
    }

    /// Sets the extended header
    #[must_use]
    pub fn x_header(mut self, x_header: &'a [u8]) -> Self {
        self.x_header = Some(x_header);
        self
    }

    /// Sets the extended header from an Option
    #[must_use]
    pub fn opt_x_header(mut self, x_header: Option<&'a [u8]>) -> Self {
        self.x_header = x_header;
        self
    }

    /// Sets the flag value
    #[must_use]
    pub fn flag(mut self, flag: u64) -> Self {
        self.flag = Some(flag);
        self
    }

    /// Sets the flag value from an Option
    #[must_use]
    pub fn opt_flag(mut self, flag: Option<u64>) -> Self {
        self.flag = flag;
        self
    }

    /// Builds the `SequencingRecord`
    ///
    /// # Errors
    ///
    /// Returns an error if the primary sequence (`s_seq`) is not set.
    pub fn build(self) -> Result<SequencingRecord<'a>> {
        let Some(s_seq) = self.s_seq else {
            return Err(WriteError::MissingSequence.into());
        };
        Ok(SequencingRecord {
            s_seq,
            s_qual: self.s_qual,
            s_header: self.s_header,
            x_seq: self.x_seq,
            x_qual: self.x_qual,
            x_header: self.x_header,
            flag: self.flag,
        })
    }
}
