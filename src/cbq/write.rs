use std::io;

use zstd::zstd_safe;

use crate::{
    Result, SequencingRecord,
    cbq::core::{BlockHeader, ColumnarBlock, FileHeader, Index, IndexFooter, IndexHeader},
};

/// Writer for CBQ files operating on generic writers (streaming).
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

    /// Push a record to the writer
    ///
    /// Returns `Ok(true)` if the record was written successfully.
    /// CBQ handles N's explicitly in its encoding, so records are never skipped.
    pub fn push(&mut self, record: SequencingRecord) -> Result<bool> {
        if !self.block.can_fit(&record) {
            self.flush()?;
        }
        self.block.push(record)?;
        Ok(true)
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

    /// Ingest only the *completed* (already-compressed) blocks from `other`.
    ///
    /// Unlike [`ingest`](Self::ingest), this never touches either writer's
    /// incomplete block, so it performs no zstd compression. The work done
    /// under a global lock is reduced to a `write_all` of pre-compressed bytes
    /// plus a header copy — compression has already been paid for on the worker
    /// thread when `other`'s blocks were flushed in `push`.
    ///
    /// `other` keeps building its incomplete block across calls; only its
    /// completed-block buffer and headers are drained.
    pub fn ingest_completed(&mut self, other: &mut ColumnarBlockWriter<Vec<u8>>) -> Result<()> {
        if other.headers.is_empty() {
            return Ok(()); // short-circuit
        }

        // Write all completed blocks from the other
        self.inner.write_all(other.inner_data())?;

        // Take all headers from the other
        self.headers.extend_from_slice(&other.headers);

        // Clear only the drained completed-block state, leaving the incomplete
        // block intact so the worker thread can keep accumulating into it.
        other.clear_completed_data();

        Ok(())
    }

    /// Ingests only the *incomplete* (non-compressed) blocks from the `other`.
    ///
    /// This should not be used in isolation and should be handled from the [`ingest`](Self::ingest) API only
    /// to avoid any mistakes.
    ///
    /// [`ingest_completed`](Self::ingest_completed) should always be called first.
    fn ingest_incompleted(&mut self, other: &mut ColumnarBlockWriter<Vec<u8>>) -> Result<()> {
        if other.block.num_records == 0 {
            return Ok(()); // short-circuit
        }

        // Attempt to ingest the incomplete block from the other
        if !self.block.can_ingest(&other.block) {
            // Make space by flushing the current block
            // Then ingest the incomplete block from the other
            self.flush()?;
        }
        self.block.take_incomplete(&other.block)?;

        // clear the drained incomplete-block state
        other.clear_incomplete_data();

        Ok(())
    }

    pub fn ingest(&mut self, other: &mut ColumnarBlockWriter<Vec<u8>>) -> Result<()> {
        self.ingest_completed(other)?;
        self.ingest_incompleted(other)?;
        Ok(())
    }
}

/// Specialized implementation when using a local `Vec<u8>` as the inner data structure
impl ColumnarBlockWriter<Vec<u8>> {
    #[must_use]
    pub fn inner_data(&self) -> &[u8] {
        &self.inner
    }

    /// Clears only the completed-block state (compressed bytes + headers),
    /// leaving the incomplete block intact.
    ///
    /// Used by [`ingest_completed`](ColumnarBlockWriter::ingest_completed) so a
    /// worker thread can keep accumulating records into its in-progress block
    /// across batches.
    pub fn clear_completed_data(&mut self) {
        self.inner.clear();
        self.headers.clear();
    }

    /// Clears the incomplete-block state
    pub fn clear_incomplete_data(&mut self) {
        self.block.clear();
    }

    /// Returns the number of bytes written to the inner data structure
    #[must_use]
    pub fn bytes_written(&self) -> usize {
        self.inner.len()
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;
    use crate::{
        BinseqRecord, SequencingRecordBuilder,
        cbq::{BlockRange, Reader, core::FileHeaderBuilder},
    };

    /// Build a `FileHeader` for sequence-only records with the given block size.
    fn header(block_size: usize) -> FileHeader {
        FileHeaderBuilder::default()
            .is_paired(false)
            .with_headers(false)
            .with_qualities(false)
            .with_flags(false)
            .with_block_size(block_size)
            .build()
    }

    /// Build a sequence-only record from a sequence slice.
    fn record(seq: &[u8]) -> SequencingRecord<'_> {
        SequencingRecordBuilder::default()
            .s_seq(seq)
            .build()
            .expect("failed to build record")
    }

    /// Read every sequence back from a finished CBQ byte buffer.
    fn read_all_sequences(bytes: Vec<u8>) -> Vec<Vec<u8>> {
        let mut reader = Reader::new(Cursor::new(bytes)).expect("failed to open reader");
        let mut out = Vec::new();
        let mut cumulative = 0u64;
        while let Some(block_header) = reader.read_block().expect("failed to read block") {
            cumulative += block_header.num_records;
            // `read_block` only loads the compressed columns; decode them before iterating.
            reader
                .block
                .decompress_columns()
                .expect("failed to decompress block");
            let range = BlockRange::new(0, cumulative);
            for rec in reader.block.iter_records(range) {
                out.push(rec.sseq().to_vec());
            }
        }
        out
    }

    /// 64 distinct fixed-length sequences (each unique by index).
    fn sample_sequences() -> Vec<Vec<u8>> {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        (0..64u32)
            .map(|i| {
                (0..40)
                    .map(|j| BASES[(i as usize + j) % 4])
                    .collect::<Vec<u8>>()
            })
            .collect()
    }

    /// `ingest_completed` must leave the source's incomplete block untouched
    /// while draining its completed blocks and headers.
    #[test]
    fn test_ingest_completed_preserves_incomplete_block() -> Result<()> {
        // Small block size so a handful of records fills multiple blocks.
        let block_size = 64;
        let mut global = ColumnarBlockWriter::new(Vec::new(), header(block_size))?;
        let mut local = ColumnarBlockWriter::new_headless(Vec::new(), header(block_size))?;

        // Push enough records to flush several completed blocks plus a partial tail.
        let seqs = sample_sequences();
        for seq in &seqs {
            local.push(record(seq))?;
        }

        // There must be both completed blocks (compressed bytes) and an
        // in-progress incomplete block at this point.
        assert!(
            !local.inner_data().is_empty(),
            "expected completed blocks before ingest_completed"
        );
        assert!(
            local.block.num_records > 0,
            "expected a non-empty incomplete block before ingest_completed"
        );

        let incomplete_records_before = local.block.num_records;
        let completed_headers = local.headers.len();

        global.ingest_completed(&mut local)?;

        // The completed-block buffer and headers are drained from the source...
        assert!(
            local.inner_data().is_empty(),
            "completed bytes should be drained"
        );
        assert!(local.headers.is_empty(), "headers should be drained");

        // ...but the incomplete block is left intact for further accumulation.
        assert_eq!(
            local.block.num_records, incomplete_records_before,
            "incomplete block must be preserved across ingest_completed"
        );

        // The global writer received exactly the completed headers (no flush of
        // its own empty incomplete block).
        assert_eq!(global.headers.len(), completed_headers);
        assert_eq!(
            global.block.num_records, 0,
            "ingest_completed must not touch the global incomplete block"
        );

        Ok(())
    }

    /// The parallel pattern (`ingest_completed` per batch, then a final full
    /// `ingest`) must round-trip every record in order.
    #[test]
    fn test_batched_ingest_completed_then_finish_roundtrips() -> Result<()> {
        let block_size = 64;
        let mut global = ColumnarBlockWriter::new(Vec::new(), header(block_size))?;
        let mut local = ColumnarBlockWriter::new_headless(Vec::new(), header(block_size))?;

        let seqs = sample_sequences();

        // Simulate batches: push a chunk, drain completed blocks, repeat.
        for chunk in seqs.chunks(7) {
            for seq in chunk {
                local.push(record(seq))?;
            }
            global.ingest_completed(&mut local)?;
        }

        // Final thread completion: drain the residual incomplete block.
        global.ingest(&mut local)?;
        global.finish()?;

        // The source has been fully drained.
        assert!(local.inner_data().is_empty());
        assert_eq!(local.block.num_records, 0);

        let read_back = read_all_sequences(global.inner);
        assert_eq!(read_back, seqs, "round-trip mismatch after batched ingest");

        Ok(())
    }

    /// A single full `ingest` (no prior `ingest_completed`) must round-trip.
    #[test]
    fn test_full_ingest_roundtrips() -> Result<()> {
        let block_size = 64;
        let mut global = ColumnarBlockWriter::new(Vec::new(), header(block_size))?;
        let mut local = ColumnarBlockWriter::new_headless(Vec::new(), header(block_size))?;

        let seqs = sample_sequences();
        for seq in &seqs {
            local.push(record(seq))?;
        }

        global.ingest(&mut local)?;
        global.finish()?;

        let read_back = read_all_sequences(global.inner);
        assert_eq!(read_back, seqs);

        Ok(())
    }

    /// Multiple local writers draining into one global writer (as distinct
    /// threads would) must preserve every record without loss or duplication.
    ///
    /// Note: input *order* is not preserved across independent locals. When one
    /// local's full `ingest` leaves a tail in the global incomplete block, a
    /// later local's `ingest_completed` writes its completed blocks ahead of
    /// that buffered tail. The guarantee is multiset equality, not sequence
    /// equality — which matches real parallel writing, where threads merge into
    /// the global in lock-acquisition order rather than input order.
    #[test]
    fn test_multiple_locals_ingest_into_global() -> Result<()> {
        let block_size = 64;
        let mut global = ColumnarBlockWriter::new(Vec::new(), header(block_size))?;

        let seqs = sample_sequences();
        let mut expected = Vec::new();

        // Three "threads", each owning a third of the records.
        for group in seqs.chunks(seqs.len().div_ceil(3)) {
            let mut local = ColumnarBlockWriter::new_headless(Vec::new(), header(block_size))?;
            for seq in group {
                local.push(record(seq))?;
                expected.push(seq.clone());
            }
            // per-batch drain followed by a thread-final full ingest
            global.ingest_completed(&mut local)?;
            global.ingest(&mut local)?;
        }

        global.finish()?;

        let mut read_back = read_all_sequences(global.inner);
        read_back.sort();
        expected.sort();
        assert_eq!(read_back, expected, "records lost or duplicated on merge");

        Ok(())
    }

    /// `ingest_completed` on a source with no completed blocks is a no-op for
    /// the global writer and preserves the source's incomplete block.
    #[test]
    fn test_ingest_completed_no_completed_blocks() -> Result<()> {
        // Large block size so a few small records never fill a block.
        let block_size = 1 << 20;
        let mut global = ColumnarBlockWriter::new(Vec::new(), header(block_size))?;
        let mut local = ColumnarBlockWriter::new_headless(Vec::new(), header(block_size))?;

        local.push(record(b"ACGTACGTACGT"))?;
        assert!(
            local.inner_data().is_empty(),
            "no block should have flushed"
        );

        let records_before = local.block.num_records;
        global.ingest_completed(&mut local)?;

        assert_eq!(local.block.num_records, records_before);
        assert!(global.headers.is_empty());
        assert_eq!(global.block.num_records, 0);

        Ok(())
    }

    // ==================== N-nucleotide round-trip (issue #94) ====================
    //
    // CBQ stores `N` as a 2-bit placeholder `A` plus an Elias-Fano index of the
    // N-positions that is used to backfill `N` on decode. These tests guard the
    // serial `Reader`/`decompress_columns` path, which previously dropped every N
    // (decoding it back as `A`). See ArcInstitute/binseq#94.

    /// Build a `FileHeader` with quality scores and headers enabled.
    fn header_with_quality_headers(block_size: usize) -> FileHeader {
        FileHeaderBuilder::default()
            .is_paired(false)
            .with_headers(true)
            .with_qualities(true)
            .with_flags(false)
            .with_block_size(block_size)
            .build()
    }

    /// Read back every record's `(sequence, quality, header)` from a finished CBQ buffer.
    fn read_all_records(bytes: Vec<u8>) -> Vec<(Vec<u8>, Vec<u8>, Vec<u8>)> {
        let mut reader = Reader::new(Cursor::new(bytes)).expect("failed to open reader");
        let mut out = Vec::new();
        let mut cumulative = 0u64;
        while let Some(block_header) = reader.read_block().expect("failed to read block") {
            cumulative += block_header.num_records;
            reader
                .block
                .decompress_columns()
                .expect("failed to decompress block");
            let range = BlockRange::new(0, cumulative);
            for rec in reader.block.iter_records(range) {
                out.push((
                    rec.sseq().to_vec(),
                    rec.squal().to_vec(),
                    rec.sheader().to_vec(),
                ));
            }
        }
        out
    }

    /// Sequences with `N` at every notable position: leading, interior, trailing,
    /// all-N, and a no-N control. Direct reproduction of ArcInstitute/binseq#94.
    #[test]
    fn test_n_roundtrip_serial_edge_positions() -> Result<()> {
        let seqs: Vec<&[u8]> = vec![
            b"NACGTACGT", // leading
            b"ACGTNACGT", // interior
            b"ACGTACGTN", // trailing
            b"NNNNNNNN",  // all-N
            b"ACGTACGT",  // no-N control
        ];

        // Large block size: all records land in a single block.
        let mut writer = ColumnarBlockWriter::new(Vec::new(), header(1 << 20))?;
        for &seq in &seqs {
            writer.push(record(seq))?;
        }
        writer.finish()?;

        let read_back = read_all_sequences(writer.inner);
        let expected: Vec<Vec<u8>> = seqs.iter().map(|s| s.to_vec()).collect();
        assert_eq!(read_back, expected, "N positions must round-trip verbatim");
        Ok(())
    }

    /// Many N-containing sequences with a small block size, forcing multiple
    /// blocks. Guards per-block Elias-Fano correctness and `ef_bytes`/`len_nef`
    /// reuse across `clear()`.
    #[test]
    fn test_n_roundtrip_serial_multi_block() -> Result<()> {
        let seqs: Vec<Vec<u8>> = (0..64u32)
            .map(|i| {
                if i % 5 == 0 {
                    // occasionally all-N
                    vec![b'N'; 8]
                } else {
                    // vary N placement so per-block position sets differ
                    let mut s = vec![b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T'];
                    let pos = i as usize % s.len();
                    s[pos] = b'N';
                    s
                }
            })
            .collect();

        let mut writer = ColumnarBlockWriter::new(Vec::new(), header(64))?;
        for seq in &seqs {
            writer.push(record(seq))?;
        }
        writer.finish()?;

        // Sanity: the small block size must actually span more than one block.
        assert!(
            writer.headers.len() > 1,
            "expected multiple blocks, got {}",
            writer.headers.len()
        );

        let read_back = read_all_sequences(writer.inner);
        assert_eq!(
            read_back, seqs,
            "N round-trip must hold across multiple blocks"
        );
        Ok(())
    }

    /// Mirrors the exact reproduction in ArcInstitute/binseq#94: quality + headers
    /// enabled, with N at leading / all / interior positions.
    #[test]
    fn test_n_roundtrip_with_quality_and_headers() -> Result<()> {
        let records: [(&[u8], &[u8], &[u8]); 3] = [
            (b"r1", b"NACGTACGT", b"IIIIIIIII"),
            (b"r2", b"NNNNNNNN", b"IIIIIIII"),
            (b"r3", b"ACGTNACGT", b"IIIIIIIII"),
        ];

        let mut writer =
            ColumnarBlockWriter::new(Vec::new(), header_with_quality_headers(1 << 20))?;
        for (h, s, q) in records {
            let rec = SequencingRecordBuilder::default()
                .s_header(h)
                .s_seq(s)
                .s_qual(q)
                .build()?;
            writer.push(rec)?;
        }
        writer.finish()?;

        let read_back = read_all_records(writer.inner);
        let expected: Vec<(Vec<u8>, Vec<u8>, Vec<u8>)> = records
            .iter()
            .map(|(h, s, q)| (s.to_vec(), q.to_vec(), h.to_vec()))
            .collect();
        assert_eq!(
            read_back, expected,
            "seq/qual/header must round-trip with N present"
        );
        Ok(())
    }
}
