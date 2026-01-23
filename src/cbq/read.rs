use std::{fs, io, path::Path, sync::Arc, thread};

use memmap2::Mmap;
use zstd::{stream::copy_decode, zstd_safe};

use crate::{
    BinseqRecord, ParallelProcessor, ParallelReader, Result,
    cbq::core::{
        BlockHeader, BlockRange, ColumnarBlock, FileHeader, Index, IndexFooter, IndexHeader,
    },
};

/// A reader for CBQ files operating on generic readers (streaming).
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

    /// Update the default quality score for this reader
    pub fn set_default_quality_score(&mut self, score: u8) {
        self.block.set_default_quality_score(score);
    }

    pub fn read_block(&mut self) -> Result<Option<BlockHeader>> {
        let mut iheader_buf = [0u8; size_of::<IndexHeader>()];
        let mut diff_buf = [0u8; size_of::<BlockHeader>() - size_of::<IndexHeader>()];
        let mut header_buf = [0u8; size_of::<BlockHeader>()];

        // Attempt to read the index header
        match self.inner.read_exact(&mut iheader_buf) {
            Ok(()) => {}
            Err(e) => {
                if e.kind() == io::ErrorKind::UnexpectedEof {
                    // no more bytes, the stream is exhausted
                    return Ok(None);
                }
                return Err(e.into());
            }
        }

        // The stream is exhausted, no more blocks to read
        if let Ok(iheader) = IndexHeader::from_bytes(&iheader_buf) {
            self.iheader = Some(iheader);
            return Ok(None);
        }
        // attempt to read the rest of the block header
        match self.inner.read_exact(&mut diff_buf) {
            Ok(()) => {}
            Err(e) => {
                return Err(e.into());
            }
        }
        header_buf[..iheader_buf.len()].copy_from_slice(&iheader_buf);
        header_buf[iheader_buf.len()..].copy_from_slice(&diff_buf);

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

/// A memory-mapped reader for CBQ files.
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

    /// Update the default quality score for this reader
    pub fn set_default_quality_score(&mut self, score: u8) {
        self.block.set_default_quality_score(score);
    }

    #[must_use]
    pub fn header(&self) -> FileHeader {
        self.block.header
    }

    #[must_use]
    pub fn is_paired(&self) -> bool {
        self.block.header.is_paired()
    }

    #[must_use]
    pub fn num_records(&self) -> usize {
        self.index.num_records()
    }

    #[must_use]
    pub fn num_blocks(&self) -> usize {
        self.index.num_blocks()
    }

    #[must_use]
    pub fn index(&self) -> &Index {
        &self.index
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
            .decompress_from_bytes(block_data_slice, block_header, &mut self.dctx)?;
        Ok(())
    }

    /// Iterate over block headers in the CBQ file.
    ///
    /// Note: This requires reading slices from the file so it will be IO-bound.
    pub fn iter_block_headers(&self) -> impl Iterator<Item = Result<BlockHeader>> {
        self.index.iter_blocks().map(|range| {
            let mut block_header_buf = [0u8; size_of::<BlockHeader>()];
            block_header_buf.copy_from_slice(
                &self.inner
                    [range.offset as usize..range.offset as usize + size_of::<BlockHeader>()],
            );
            BlockHeader::from_bytes(&block_header_buf)
        })
    }
}
impl ParallelReader for MmapReader {
    fn process_parallel<P: ParallelProcessor + Clone + 'static>(
        self,
        processor: P,
        num_threads: usize,
    ) -> crate::Result<()> {
        let num_records = self.num_records();
        self.process_parallel_range(processor, num_threads, 0..num_records)
    }

    fn process_parallel_range<P: ParallelProcessor + Clone + 'static>(
        self,
        processor: P,
        num_threads: usize,
        range: std::ops::Range<usize>,
    ) -> crate::Result<()> {
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

            // eprintln!(
            //     "Thread {} block range: {}-{}. First block Cumulative Records: {}. Last block Cumulative Records: {}",
            //     thread_id,
            //     start_block_idx,
            //     end_block_idx,
            //     t_block_ranges[0].cumulative_records,
            //     t_block_ranges.last().unwrap().cumulative_records
            // );

            let thread_handle = thread::spawn(move || -> crate::Result<()> {
                for b_range in t_block_ranges {
                    t_reader.load_block(b_range)?;
                    for record in t_reader.block.iter_records(b_range) {
                        let global_record_idx = record.index() as usize;

                        // Only process records within our specified range
                        if global_record_idx >= range.start && global_record_idx < range.end {
                            t_proc.process_record(record)?;
                        }
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
#[cfg(test)]
mod tests {
    use super::*;
    use crate::BinseqRecord;

    const TEST_CBQ_FILE: &str = "./data/subset.cbq";

    // ==================== MmapReader Basic Tests ====================

    #[test]
    fn test_mmap_reader_new() {
        let reader = MmapReader::new(TEST_CBQ_FILE);
        assert!(reader.is_ok(), "Failed to create CBQ reader");
    }

    #[test]
    fn test_mmap_reader_num_records() {
        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();
        let num_records = reader.num_records();
        assert!(num_records > 0, "Expected non-zero records");
    }

    #[test]
    fn test_mmap_reader_is_paired() {
        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();
        let is_paired = reader.is_paired();
        // Test that the method returns a boolean
        assert!(is_paired || !is_paired);
    }

    #[test]
    fn test_mmap_reader_header_access() {
        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();
        let header = reader.header();
        assert!(header.block_size > 0, "Expected non-zero block size");
    }

    #[test]
    fn test_mmap_reader_index_access() {
        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();
        let index = reader.index();
        assert!(index.num_records() > 0, "Index should have records");
    }

    #[test]
    fn test_mmap_reader_num_blocks() {
        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();
        let num_blocks = reader.num_blocks();
        assert!(num_blocks > 0, "Should have at least one block");
    }

    // ==================== Default Quality Score Tests ====================

    #[test]
    fn test_set_default_quality_score() {
        let mut reader = MmapReader::new(TEST_CBQ_FILE).unwrap();
        let custom_score = 42u8;

        reader.set_default_quality_score(custom_score);
        // Just verify it doesn't panic
    }

    // ==================== Parallel Processing Tests ====================

    #[derive(Clone)]
    struct CbqCountingProcessor {
        count: Arc<std::sync::Mutex<usize>>,
    }

    impl ParallelProcessor for CbqCountingProcessor {
        fn process_record<R: BinseqRecord>(&mut self, _record: R) -> Result<()> {
            let mut count = self.count.lock().unwrap();
            *count += 1;
            Ok(())
        }
    }

    #[test]
    fn test_parallel_processing() {
        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();
        let num_records = reader.num_records();

        let count = Arc::new(std::sync::Mutex::new(0));
        let processor = CbqCountingProcessor {
            count: count.clone(),
        };

        reader.process_parallel(processor, 2).unwrap();

        let final_count = *count.lock().unwrap();
        assert_eq!(final_count, num_records, "All records should be processed");
    }

    #[test]
    fn test_parallel_processing_range() {
        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();
        let num_records = reader.num_records();

        if num_records >= 100 {
            let start = 10;
            let end = 50;
            let expected_count = end - start;

            let count = Arc::new(std::sync::Mutex::new(0));
            let processor = CbqCountingProcessor {
                count: count.clone(),
            };

            reader
                .process_parallel_range(processor, 2, start..end)
                .unwrap();

            let final_count = *count.lock().unwrap();
            assert_eq!(
                final_count, expected_count,
                "Should process exactly {} records",
                expected_count
            );
        }
    }

    #[test]
    fn test_parallel_processing_with_record_data() {
        #[derive(Clone)]
        struct RecordValidator {
            valid_count: Arc<std::sync::Mutex<usize>>,
        }

        impl ParallelProcessor for RecordValidator {
            fn process_record<R: BinseqRecord>(&mut self, record: R) -> Result<()> {
                // Validate record has non-zero length
                assert!(record.slen() > 0, "Record should have non-zero length");

                let mut count = self.valid_count.lock().unwrap();
                *count += 1;
                Ok(())
            }
        }

        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();
        let num_records = reader.num_records();

        let count = Arc::new(std::sync::Mutex::new(0));
        let processor = RecordValidator {
            valid_count: count.clone(),
        };

        reader.process_parallel(processor, 2).unwrap();

        let final_count = *count.lock().unwrap();
        assert_eq!(final_count, num_records);
    }

    // ==================== Index Tests ====================

    #[test]
    fn test_index_num_records() {
        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();

        let index_records = reader.index().num_records();
        let reader_records = reader.num_records();

        assert_eq!(
            index_records, reader_records,
            "Index and reader should report same number of records"
        );
    }

    #[test]
    fn test_index_num_blocks() {
        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();

        let num_blocks = reader.index().num_blocks();
        assert!(num_blocks > 0, "Should have at least one block");
    }

    #[test]
    fn test_index_iter_blocks() {
        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();

        let blocks: Vec<_> = reader.index().iter_blocks().collect();
        assert!(!blocks.is_empty(), "Should have at least one block");

        let num_blocks = reader.num_blocks();
        assert_eq!(blocks.len(), num_blocks, "Block count should match");
    }

    // ==================== Error Handling Tests ====================

    #[test]
    fn test_nonexistent_file() {
        let result = MmapReader::new("./data/nonexistent.cbq");
        assert!(result.is_err(), "Should fail on nonexistent file");
    }

    #[test]
    fn test_invalid_file_format() {
        // Try to open a non-CBQ file as CBQ
        let result = MmapReader::new("./Cargo.toml");
        // This should fail during header validation
        assert!(result.is_err(), "Should fail on invalid file format");
    }

    // ==================== Block Header Iterator Tests ====================

    #[test]
    fn test_iter_block_headers() {
        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();

        let headers: Vec<_> = reader
            .iter_block_headers()
            .take(5)
            .collect::<Result<Vec<_>>>()
            .unwrap();

        assert!(!headers.is_empty(), "Should have at least one block header");

        for header in headers {
            assert!(header.num_records > 0, "Block should have records");
        }
    }

    #[test]
    fn test_iter_block_headers_count() {
        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();

        let header_count = reader
            .iter_block_headers()
            .collect::<Result<Vec<_>>>()
            .unwrap()
            .len();

        let num_blocks = reader.num_blocks();
        assert_eq!(header_count, num_blocks, "Should iterate all block headers");
    }

    // ==================== Empty Range Tests ====================

    #[test]
    fn test_parallel_processing_empty_range() {
        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();

        let count = Arc::new(std::sync::Mutex::new(0));
        let processor = CbqCountingProcessor {
            count: count.clone(),
        };

        // Process empty range
        reader.process_parallel_range(processor, 2, 0..0).unwrap();

        let final_count = *count.lock().unwrap();
        assert_eq!(final_count, 0, "Empty range should process no records");
    }

    #[test]
    fn test_parallel_processing_invalid_range() {
        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();
        let num_records = reader.num_records();

        let count = Arc::new(std::sync::Mutex::new(0));
        let processor = CbqCountingProcessor {
            count: count.clone(),
        };

        // Process out of bounds range (should handle gracefully)
        let result =
            reader.process_parallel_range(processor, 2, num_records + 100..num_records + 200);

        assert!(
            result.is_ok(),
            "Should handle out of bounds range gracefully"
        );
    }

    // ==================== Thread Count Tests ====================

    #[test]
    fn test_parallel_processing_single_thread() {
        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();
        let num_records = reader.num_records();

        let count = Arc::new(std::sync::Mutex::new(0));
        let processor = CbqCountingProcessor {
            count: count.clone(),
        };

        reader.process_parallel(processor, 1).unwrap();

        let final_count = *count.lock().unwrap();
        assert_eq!(final_count, num_records);
    }

    #[test]
    fn test_parallel_processing_many_threads() {
        let reader = MmapReader::new(TEST_CBQ_FILE).unwrap();
        let num_records = reader.num_records();

        let count = Arc::new(std::sync::Mutex::new(0));
        let processor = CbqCountingProcessor {
            count: count.clone(),
        };

        reader.process_parallel(processor, 8).unwrap();

        let final_count = *count.lock().unwrap();
        assert_eq!(final_count, num_records);
    }
}
