//! Compares the two zstd backends (`zstd-c` vs `zstd-rs`, see `src/zstdlib.rs`)
//! end-to-end through the public CBQ writer/reader, since the backend is an
//! internal implementation detail rather than part of the public API.
//!
//! Run under each backend to compare:
//! ```sh
//! cargo bench --no-default-features --features paraseq,anyhow,zstd-c
//! cargo bench --no-default-features --features paraseq,anyhow,zstd-rs
//! ```

use std::hint::black_box;
use std::io::Cursor;

use binseq::SequencingRecordBuilder;
use binseq::cbq::{BlockRange, ColumnarBlockWriter, FileHeaderBuilder, Reader};
use criterion::{Criterion, criterion_group, criterion_main};

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
const NUM_RECORDS: usize = 5_000;
const SEQ_LEN: usize = 150;

fn make_sequences(n: usize, len: usize) -> Vec<Vec<u8>> {
    (0..n)
        .map(|i| (0..len).map(|j| BASES[(i + j) % 4]).collect())
        .collect()
}

fn write_block(seqs: &[Vec<u8>], quals: &[u8], compression_level: usize) -> Vec<u8> {
    let header = FileHeaderBuilder::default()
        .is_paired(false)
        .with_headers(false)
        .with_qualities(true)
        .with_flags(false)
        // large enough that every record lands in a single block, so a run
        // measures one compress/decompress pass rather than many small ones.
        .with_block_size(seqs.len() * (SEQ_LEN + 64))
        .with_compression_level(compression_level)
        .build();
    let mut writer = ColumnarBlockWriter::new(Vec::new(), header).expect("failed to create writer");
    for seq in seqs {
        let record = SequencingRecordBuilder::default()
            .s_seq(seq)
            .s_qual(quals)
            .build()
            .expect("failed to build record");
        writer.push(record).expect("failed to push record");
    }
    writer.finish().expect("failed to finish writer");
    writer.inner_data().to_vec()
}

fn read_all_records(bytes: &[u8]) -> usize {
    let mut reader = Reader::new(Cursor::new(bytes)).expect("failed to open reader");
    let mut cumulative = 0u64;
    let mut count = 0;
    while let Some(block_header) = reader.read_block().expect("failed to read block") {
        cumulative += block_header.num_records;
        reader
            .block
            .decompress_columns()
            .expect("failed to decompress block");
        let range = BlockRange::new(0, cumulative);
        count += reader.block.iter_records(range).count();
    }
    count
}

fn bench_compress(c: &mut Criterion) {
    let seqs = make_sequences(NUM_RECORDS, SEQ_LEN);
    let quals = vec![b'I'; SEQ_LEN];

    let mut group = c.benchmark_group("cbq_compress");
    for level in [1usize, 3, 9] {
        group.bench_function(format!("level_{level}"), |b| {
            b.iter(|| black_box(write_block(black_box(&seqs), black_box(&quals), level)));
        });
    }
    group.finish();
}

fn bench_decompress(c: &mut Criterion) {
    let seqs = make_sequences(NUM_RECORDS, SEQ_LEN);
    let quals = vec![b'I'; SEQ_LEN];
    let bytes = write_block(&seqs, &quals, 3);

    let mut group = c.benchmark_group("cbq_decompress");
    group.bench_function("read_all_records", |b| {
        b.iter(|| black_box(read_all_records(black_box(&bytes))));
    });
    group.finish();
}

criterion_group!(benches, bench_compress, bench_decompress);
criterion_main!(benches);
