use auto_impl::auto_impl;
use bitnuc::BitSize;

use crate::Result;

/// Record trait shared between BINSEQ variants.
///
/// Exposes public methods for accessing internal data.
/// Interfaces with the [`bitnuc`] crate for decoding sequences.
///
/// Implemented by [`bq::RefRecord`](crate::bq::RefRecord) and [`vbq::RefRecord`](crate::vbq::RefRecord).
///
/// Used to interact with [`ParallelProcessor`](crate::ParallelProcessor) for easy parallel processing.
#[auto_impl(&, &mut)]
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
