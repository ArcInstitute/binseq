use bitnuc_deprec::BitSize;

use crate::Result;

/// Record trait shared between BINSEQ variants.
///
/// Exposes public methods for accessing internal data.
/// Interfaces with the [`bitnuc`] crate for decoding sequences.
///
/// Implemented by [`bq::RefRecord`](crate::bq::RefRecord) and [`vbq::RefRecord`](crate::vbq::RefRecord).
///
/// Used to interact with [`ParallelProcessor`](crate::ParallelProcessor) for easy parallel processing.
pub trait BinseqRecord {
    /// Returns the bitsize of the record (number of bits per nucleotide)
    fn bitsize(&self) -> BitSize;

    /// Returns the global index of the record.
    fn index(&self) -> u64;

    /// Returns the flag value of this record
    fn flag(&self) -> Option<u64>;

    /// Returns the header of this record
    fn sheader(&self) -> &[u8];

    /// Returns the header of the extended/paired sequence (empty if not paired)
    fn xheader(&self) -> &[u8];

    /// Returns the length of the primary sequence of this record
    fn slen(&self) -> u64;

    /// Returns the length of the extended sequence of this record
    fn xlen(&self) -> u64;

    /// Returns a reference to the **encoded** primary sequence of this record
    fn sbuf(&self) -> &[u64];

    /// Returns a reference to the **encoded** extended sequence of this record.
    ///
    /// Empty if no extended sequence is present.
    fn xbuf(&self) -> &[u64];

    /// Returns a reference to the quality scores of the primary sequence of this record.
    ///
    /// Empty if no quality scores are present.
    fn squal(&self) -> &[u8] {
        &[]
    }

    /// Returns a reference to the quality scores of the extended sequence of this record.
    ///
    /// Empty if no quality scores are present.
    fn xqual(&self) -> &[u8] {
        &[]
    }

    /// Decodes the primary sequence of this record into the provided buffer.
    fn decode_s(&self, buf: &mut Vec<u8>) -> Result<()> {
        self.bitsize()
            .decode(self.sbuf(), self.slen() as usize, buf)?;
        Ok(())
    }

    /// Decodes the extended sequence of this record into the provided buffer.
    fn decode_x(&self, buf: &mut Vec<u8>) -> Result<()> {
        self.bitsize()
            .decode(self.xbuf(), self.xlen() as usize, buf)?;
        Ok(())
    }

    /// Returns a reference to the primary decoded sequence of this record.
    ///
    /// This is not available on all types that implement the `Record` trait.
    /// It should be available on types that implement it in this library however.
    fn sseq(&self) -> &[u8] {
        unimplemented!("This record does not implement direct sequence access");
    }

    /// Returns a reference to the extended decoded sequence of this record.
    ///
    /// This may not be available on all types that implement the `Record` trait.
    /// It should be available on types that implement it in this library however.
    fn xseq(&self) -> &[u8] {
        unimplemented!("This record does not implement direct sequence access");
    }

    /// Decodes the primary sequence of this record into a newly allocated buffer.
    ///
    /// Not advised to use this function as it allocates a new buffer every time.
    fn decode_s_alloc(&self) -> Result<Vec<u8>> {
        let mut buf = Vec::with_capacity(self.slen() as usize);
        self.decode_s(&mut buf)?;
        Ok(buf)
    }

    /// Decodes the extended sequence of this record into a newly allocated buffer.
    ///
    /// Not advised to use this function as it allocates a new buffer every time.
    fn decode_x_alloc(&self) -> Result<Vec<u8>> {
        let mut buf = Vec::with_capacity(self.xlen() as usize);
        self.decode_x(&mut buf)?;
        Ok(buf)
    }

    /// A convenience function to check if the record is paired.
    fn is_paired(&self) -> bool {
        self.xlen() > 0
    }

    /// A convenience function to check if record has associated quality scores
    fn has_quality(&self) -> bool {
        !self.squal().is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Minimal implementation exercising only the trait's required methods,
    /// so that the default-method implementations get covered.
    struct MockRecord {
        bitsize: BitSize,
        index: u64,
        flag: Option<u64>,
        sbuf: Vec<u64>,
        xbuf: Vec<u64>,
        slen: u64,
        xlen: u64,
        squal: Vec<u8>,
    }

    impl BinseqRecord for MockRecord {
        fn bitsize(&self) -> BitSize {
            self.bitsize
        }
        fn index(&self) -> u64 {
            self.index
        }
        fn flag(&self) -> Option<u64> {
            self.flag
        }
        fn sheader(&self) -> &[u8] {
            b"seq_header"
        }
        fn xheader(&self) -> &[u8] {
            b""
        }
        fn slen(&self) -> u64 {
            self.slen
        }
        fn xlen(&self) -> u64 {
            self.xlen
        }
        fn sbuf(&self) -> &[u64] {
            &self.sbuf
        }
        fn xbuf(&self) -> &[u64] {
            &self.xbuf
        }
        fn squal(&self) -> &[u8] {
            &self.squal
        }
    }

    fn unpaired_record() -> MockRecord {
        let seq = b"ACGTACGTAC";
        let mut sbuf = Vec::new();
        BitSize::Two.encode(seq, &mut sbuf).unwrap();
        MockRecord {
            bitsize: BitSize::Two,
            index: 7,
            flag: Some(3),
            sbuf,
            xbuf: Vec::new(),
            slen: seq.len() as u64,
            xlen: 0,
            squal: Vec::new(),
        }
    }

    fn paired_record() -> MockRecord {
        let sseq = b"ACGTACGTAC";
        let xseq = b"TTGGCCAATT";
        let mut sbuf = Vec::new();
        let mut xbuf = Vec::new();
        BitSize::Two.encode(sseq, &mut sbuf).unwrap();
        BitSize::Two.encode(xseq, &mut xbuf).unwrap();
        MockRecord {
            bitsize: BitSize::Two,
            index: 1,
            flag: None,
            sbuf,
            xbuf,
            slen: sseq.len() as u64,
            xlen: xseq.len() as u64,
            squal: vec![b'I'; sseq.len()],
        }
    }

    #[test]
    fn test_default_index_and_flag() {
        let record = unpaired_record();
        assert_eq!(record.index(), 7);
        assert_eq!(record.flag(), Some(3));
    }

    #[test]
    fn test_decode_s() {
        let record = unpaired_record();
        let mut buf = Vec::new();
        record.decode_s(&mut buf).unwrap();
        assert_eq!(buf, b"ACGTACGTAC");
    }

    #[test]
    fn test_decode_x() {
        let record = paired_record();
        let mut buf = Vec::new();
        record.decode_x(&mut buf).unwrap();
        assert_eq!(buf, b"TTGGCCAATT");
    }

    #[test]
    fn test_decode_s_alloc() {
        let record = unpaired_record();
        let buf = record.decode_s_alloc().unwrap();
        assert_eq!(buf, b"ACGTACGTAC");
    }

    #[test]
    fn test_decode_x_alloc() {
        let record = paired_record();
        let buf = record.decode_x_alloc().unwrap();
        assert_eq!(buf, b"TTGGCCAATT");
    }

    #[test]
    #[should_panic(expected = "does not implement direct sequence access")]
    fn test_sseq_default_panics() {
        let record = unpaired_record();
        let _ = record.sseq();
    }

    #[test]
    #[should_panic(expected = "does not implement direct sequence access")]
    fn test_xseq_default_panics() {
        let record = unpaired_record();
        let _ = record.xseq();
    }

    #[test]
    fn test_is_paired() {
        assert!(!unpaired_record().is_paired());
        assert!(paired_record().is_paired());
    }

    #[test]
    fn test_has_quality() {
        assert!(!unpaired_record().has_quality());
        assert!(paired_record().has_quality());
    }

    #[test]
    fn test_default_xqual_is_empty() {
        let record = unpaired_record();
        assert!(record.xqual().is_empty());
    }
}
