use std::io;

use anyhow::Result;
use zstd::zstd_safe;

use crate::cbq::core::{
    BlockHeader, ColumnarBlock, FileHeader, Index, IndexFooter, IndexHeader, SequencingRecord,
};

pub struct ColumnarBlockWriter<W: io::Write> {
    /// Internal writer for the block
    inner: W,

    /// A reusable block for this writer
    block: ColumnarBlock,

    /// All block headers written by this writer
    headers: Vec<BlockHeader>,

    /// Compression context for the thread
    cctx: zstd_safe::CCtx<'static>,
}
impl<W: io::Write + Clone> Clone for ColumnarBlockWriter<W> {
    fn clone(&self) -> Self {
        let mut writer = Self {
            inner: self.inner.clone(),
            block: self.block.clone(),
            headers: self.headers.clone(),
            cctx: zstd_safe::CCtx::create(),
        };
        writer
            .init_compressor()
            .expect("Failed to set compression level in writer clone");
        writer
    }
}
impl<W: io::Write> ColumnarBlockWriter<W> {
    /// Creates a new writer with the header written to the inner writer
    pub fn new(inner: W, header: FileHeader) -> Result<Self> {
        // Build the writer
        let mut writer = Self::new_headless(inner, header)?;

        // Ensure the header is written to the file
        writer.inner.write_all(header.as_bytes())?;

        Ok(writer)
    }

    /// Creates a new writer without writing the header to the inner writer
    pub fn new_headless(inner: W, header: FileHeader) -> Result<Self> {
        let mut writer = Self {
            inner,
            block: ColumnarBlock::new(header),
            headers: Vec::default(),
            cctx: zstd_safe::CCtx::create(),
        };

        // Set the compression level for this writer
        writer.init_compressor()?;

        Ok(writer)
    }

    /// Sets the compression level for Writer
    ///
    /// Note: only used on init, shouldn't be set by the user
    fn init_compressor(&mut self) -> Result<()> {
        // Initialize the compressor with the compression level
        self.cctx
            .set_parameter(zstd_safe::CParameter::CompressionLevel(
                self.block.header.compression_level as i32,
            ))
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;

        // Set long distance matching
        self.cctx
            .set_parameter(zstd_safe::CParameter::EnableLongDistanceMatching(true))
            .map_err(|e| io::Error::other(zstd_safe::get_error_name(e)))?;
        Ok(())
    }

    pub fn header(&self) -> FileHeader {
        self.block.header
    }

    /// Calculate the usage of the block as a percentage
    pub fn usage(&self) -> f64 {
        self.block.usage()
    }

    pub fn push(&mut self, record: SequencingRecord) -> Result<()> {
        if !self.block.can_fit(&record) {
            self.flush()?;
        }
        self.block.push(record)?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        if let Some(header) = self.block.flush_to(&mut self.inner, &mut self.cctx)? {
            self.headers.push(header);
        }
        Ok(())
    }

    pub fn finish(&mut self) -> Result<()> {
        self.flush()?;
        self.write_index()?;
        Ok(())
    }

    fn write_index(&mut self) -> Result<()> {
        let index = Index::from_block_headers(&self.headers);
        let z_index = index.encoded()?;
        let header = IndexHeader::new(index.size(), z_index.len() as u64);
        let footer = IndexFooter::new(z_index.len() as u64);

        // Write the index to the inner writer
        {
            self.inner.write_all(header.as_bytes())?;
            self.inner.write_all(&z_index)?;
            self.inner.write_all(footer.as_bytes())?;
        }
        Ok(())
    }

    pub fn ingest(&mut self, other: &mut ColumnarBlockWriter<Vec<u8>>) -> Result<()> {
        // Write all completed blocks from the other
        self.inner.write_all(other.inner_data())?;
        // eprintln!(
        //     "Wrote {} bytes from completed blocks",
        //     other.inner_data().len()
        // );

        // Take all headers from the other
        self.headers.extend_from_slice(&other.headers);

        // Attempt to ingest the incomplete block from the other
        if self.block.can_ingest(&other.block) {
            // eprintln!("Can ingest incomplete block");
            self.block.take_incomplete(&other.block)?;

        // Make space by flushing the current block
        // Then ingest the incomplete block from the other
        } else {
            // eprintln!("Cannot ingest incomplete block");
            self.flush()?;
            self.block.take_incomplete(&other.block)?;
        }

        // Clear the other's inner data and offsets
        other.clear_inner_data();

        Ok(())
    }
}

/// Specialized implementation when using a local `Vec<u8>` as the inner data structure
impl ColumnarBlockWriter<Vec<u8>> {
    pub fn inner_data(&self) -> &[u8] {
        &self.inner
    }

    pub fn clear_inner_data(&mut self) {
        self.inner.clear();
        self.headers.clear();
        self.block.clear();
    }

    /// Returns the number of bytes written to the inner data structure
    pub fn bytes_written(&self) -> usize {
        self.inner.len()
    }
}
