use std::io;

use anyhow::Result;

use crate::core::{
    BlockHeader, ColumnarBlock, FileHeader, Index, IndexFooter, IndexHeader, SequencingRecord,
};

#[derive(Clone)]
pub struct ColumnarBlockWriter<W: io::Write> {
    /// Internal writer for the block
    inner: W,

    /// The CBQ file header
    #[allow(dead_code)]
    header: FileHeader,

    /// A reusable block for this writer
    block: ColumnarBlock,

    /// Offsets of the blocks written by this writer
    offsets: Vec<BlockHeader>,
}
impl<W: io::Write> ColumnarBlockWriter<W> {
    pub fn new(inner: W, header: FileHeader) -> Result<Self> {
        // Build the writer
        let mut writer = Self {
            inner,
            header,
            block: ColumnarBlock::new(header),
            offsets: Vec::default(),
        };

        // Ensure the header is written to the file
        writer.inner.write_all(header.as_bytes())?;

        Ok(writer)
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
}
