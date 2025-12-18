use std::{fs, io, path::Path, sync::Arc, thread};

use anyhow::Result;
use binseq::{ParallelProcessor, ParallelReader};
use memmap2::Mmap;
use zstd::{stream::copy_decode, zstd_safe};

use crate::core::{
    BlockHeader, BlockRange, ColumnarBlock, FileHeader, Index, IndexFooter, IndexHeader,
};

pub struct Reader<R: io::Read> {
    inner: R,
    pub block: ColumnarBlock,
    iheader: Option<IndexHeader>,
}
impl<R: io::Read> Reader<R> {
    pub fn new(mut inner: R) -> Result<Self> {
        let mut header_buf = [0u8; size_of::<FileHeader>()];
        inner.read_exact(&mut header_buf)?;
        let header = FileHeader::from_bytes(&header_buf)?;

        Ok(Self {
            inner,
            block: ColumnarBlock::new(header),
            iheader: None,
        })
    }

    pub fn read_block(&mut self) -> Result<Option<BlockHeader>> {
        let mut iheader_buf = [0u8; size_of::<IndexHeader>()];
        let mut diff_buf = [0u8; size_of::<BlockHeader>() - size_of::<IndexHeader>()];
        let mut header_buf = [0u8; size_of::<BlockHeader>()];

        // Attempt to read the index header
        match self.inner.read_exact(&mut iheader_buf) {
            Ok(_) => {}
            Err(e) => {
                if e.kind() == io::ErrorKind::UnexpectedEof {
                    // no more bytes, the stream is exhausted
                    return Ok(None);
                } else {
                    return Err(e.into());
                }
            }
        }

        // The stream is exhausted, no more blocks to read
        if let Ok(iheader) = IndexHeader::from_bytes(&iheader_buf) {
            self.iheader = Some(iheader);
            return Ok(None);
        } else {
            // attempt to read the rest of the block header
            match self.inner.read_exact(&mut diff_buf) {
                Ok(_) => {}
                Err(e) => {
                    return Err(e.into());
                }
            }
            header_buf[..iheader_buf.len()].copy_from_slice(&iheader_buf);
            header_buf[iheader_buf.len()..].copy_from_slice(&diff_buf);
        }

        let header = BlockHeader::from_bytes(&header_buf)?;
        self.block.read_from(&mut self.inner, header)?;

        Ok(Some(header))
    }

    pub fn read_index(&mut self) -> Result<Option<Index>> {
        let Some(header) = self.iheader else {
            return Ok(None);
        };
        let mut z_index_buf = Vec::new();
        let mut index_buf = Vec::new();
        let mut footer_buf = [0u8; size_of::<IndexFooter>()];

        // Read the index data from the reader
        z_index_buf.resize(header.z_bytes as usize, 0);

        // Reads the compressed index data
        self.inner.read_exact(&mut z_index_buf)?;
        copy_decode(z_index_buf.as_slice(), &mut index_buf)?;
        let index = Index::from_bytes(&index_buf)?;

        // Read the footer data from the reader
        self.inner.read_exact(&mut footer_buf)?;
        let _footer = IndexFooter::from_bytes(&footer_buf)?;

        Ok(Some(index))
    }
}

pub struct MmapReader {
    inner: Arc<Mmap>,
    index: Arc<Index>,

    /// Reusable record block
    block: ColumnarBlock,

    /// Reusable decompression context
    dctx: zstd_safe::DCtx<'static>,
}
impl Clone for MmapReader {
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
            index: self.index.clone(),
            block: self.block.clone(),
            dctx: zstd_safe::DCtx::create(),
        }
    }
}
impl MmapReader {
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = fs::File::open(path)?;

        // Load the mmap
        let inner = unsafe { Mmap::map(&file) }?;

        // Build the header
        let header = FileHeader::from_bytes(&inner[..size_of::<FileHeader>()])?;

        // build the index
        let index = {
            // Load the index footer
            let footer_start = inner.len() - size_of::<IndexFooter>();
            let mut footer_buf = [0u8; size_of::<IndexFooter>()];
            footer_buf.copy_from_slice(&inner[footer_start..]);
            let index_footer = IndexFooter::from_bytes(&footer_buf)?;

            // Find the coordinates of the compressed index
            let z_index_start = footer_start - index_footer.bytes as usize;
            let z_index_slice = &inner[z_index_start..footer_start];

            // Decompress the index
            let mut index_buf = Vec::default();
            copy_decode(z_index_slice, &mut index_buf)?;

            // Load the index
            Index::from_bytes(&index_buf)
        }?;

        Ok(Self {
            inner: Arc::new(inner),
            index: Arc::new(index),
            block: ColumnarBlock::new(header),
            dctx: zstd_safe::DCtx::create(),
        })
    }

    pub fn num_records(&self) -> usize {
        self.index.num_records()
    }

    fn load_block(&mut self, range: BlockRange) -> Result<()> {
        let header_start = range.offset as usize;
        let header_end = size_of::<BlockHeader>() + header_start;
        let block_header = {
            let mut block_header_buf = [0u8; size_of::<BlockHeader>()];
            block_header_buf.copy_from_slice(&self.inner[header_start..header_end]);
            BlockHeader::from_bytes(&block_header_buf)
        }?;

        let data_end = header_end + block_header.block_len();
        let block_data_slice = &self.inner[header_end..data_end];
        self.block
            .decompress_from_bytes(&block_data_slice, block_header, &mut self.dctx)?;
        Ok(())
    }
}
impl ParallelReader for MmapReader {
    fn process_parallel<P: ParallelProcessor + Clone + 'static>(
        self,
        processor: P,
        num_threads: usize,
    ) -> binseq::Result<()> {
        let num_records = self.num_records();
        self.process_parallel_range(processor, num_threads, 0..num_records)
    }

    fn process_parallel_range<P: ParallelProcessor + Clone + 'static>(
        self,
        processor: P,
        num_threads: usize,
        range: std::ops::Range<usize>,
    ) -> binseq::Result<()> {
        let num_threads = if num_threads == 0 {
            num_cpus::get()
        } else {
            num_threads.min(num_cpus::get())
        };

        // validate range
        let total_records = self.num_records();
        if range.start >= total_records || range.end > total_records || range.start > range.end {
            return Ok(()); // nothing to do
        }

        let mut iv_start = 0;
        let relevant_blocks = self
            .index
            .iter_blocks()
            .filter(|block| {
                let iv_end = block.cumulative_records as usize;
                let relevant = iv_start <= range.end && iv_end > range.start;
                iv_start = iv_end;
                relevant
            })
            .collect::<Vec<_>>();
        let num_blocks = relevant_blocks.len();

        if relevant_blocks.is_empty() {
            return Ok(()); // nothing to do
        }

        let blocks_per_thread = num_blocks.div_ceil(num_threads);

        let mut handles = Vec::new();
        for thread_id in 0..num_threads {
            let start_block_idx = thread_id * blocks_per_thread;
            let end_block_idx = ((thread_id + 1) * blocks_per_thread).min(num_blocks);

            let mut t_reader = self.clone();
            let mut t_proc = processor.clone();

            // pull all block ranges for this thread
            let t_block_ranges = relevant_blocks
                .iter()
                .skip(start_block_idx)
                .take(end_block_idx - start_block_idx)
                .copied()
                .collect::<Vec<_>>();

            let thread_handle = thread::spawn(move || -> binseq::Result<()> {
                for range in t_block_ranges {
                    t_reader.load_block(range)?;
                    for record in t_reader.block.iter_records(range) {
                        t_proc.process_record(record)?;
                    }
                    t_proc.on_batch_complete()?;
                }
                Ok(())
            });
            handles.push(thread_handle);
        }

        for handle in handles {
            handle.join().unwrap()?;
        }
        Ok(())
    }
}
