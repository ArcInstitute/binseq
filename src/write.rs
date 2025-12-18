use std::io;

use anyhow::Result;

use crate::core::{
    BlockHeader, ColumnarBlock, FileHeader, Index, IndexFooter, IndexHeader, SequencingRecord,
};

#[derive(Clone)]
pub struct ColumnarBlockWriter<W: io::Write> {
    /// Internal writer for the block
    inner: W,

    /// A reusable block for this writer
    block: ColumnarBlock,

    /// Offsets of the blocks written by this writer
    offsets: Vec<BlockHeader>,
}
impl<W: io::Write> ColumnarBlockWriter<W> {
    /// Creates a new writer with the header written to the inner writer
    pub fn new(inner: W, header: FileHeader) -> Result<Self> {
        // Build the writer
        let mut writer = Self {
            inner,
            block: ColumnarBlock::new(header),
            offsets: Vec::default(),
        };

        // Ensure the header is written to the file
        writer.inner.write_all(header.as_bytes())?;

        Ok(writer)
    }

    /// Creates a new writer without writing the header to the inner writer
    pub fn new_headless(inner: W, header: FileHeader) -> Self {
        Self {
            inner,
            block: ColumnarBlock::new(header),
            offsets: Vec::default(),
        }
    }

    pub fn push(&mut self, record: SequencingRecord) -> Result<()> {
        if !self.block.can_fit(&record) {
            self.flush()?;
        }
        self.block.push(record)?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.block.flush_to(&mut self.inner)?.map(|header| {
            self.offsets.push(header);
        });
        Ok(())
    }

    pub fn finish(&mut self) -> Result<()> {
        self.flush()?;
        self.write_index()?;
        Ok(())
    }

    fn write_index(&mut self) -> Result<()> {
        let index = Index::from_block_headers(&self.offsets);
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
        self.inner.write_all(&other.inner_data())?;

        // Take all offsets from the other
        self.offsets.extend_from_slice(&other.offsets);

        // Attempt to ingest the incomplete block from the other
        if self.block.can_ingest(&other.block) {
            self.block.take_incomplete(&other.block)?;

        // Make space by flushing the current block
        // Then ingest the incomplete block from the other
        } else {
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
        self.offsets.clear();
        self.block.clear();
    }
}
