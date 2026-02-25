# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```sh
cargo build                          # debug build
cargo test                           # run all tests
cargo test --release                 # run all tests in release mode
cargo fmt --check                    # check formatting
cargo clippy                         # lint
cargo test <test_name>               # run a single test by name
cargo run --release --example <name> -- <args>  # run an example
```

Available examples: read, write, auto-write, grep, parallel_range, streaming, network_streaming. Test data lives in `data/` (subset.bq, subset.vbq, subset.cbq, subset_R1.fastq.gz, subset_R2.fastq.gz).

## Architecture

binseq is a library (no binary targets) for reading and writing binary DNA sequence file formats. The CLI tool is in a separate repo (`bqtools`).

### Three format variants

Each lives in its own module with a reader and writer:

- **BQ** (`src/bq/`) — Fixed-length records, 2-bit nucleotide encoding, no quality scores. Simplest and most compact.
- **VBQ** (`src/vbq/`) — Variable-length records, row-based blocks, optional quality scores and headers, zstd compression, embedded index for random access.
- **CBQ** (`src/cbq/`) — Variable-length records, columnar block storage, zstd compression, tracks N bases natively with Elias-Fano encoding. Recommended format.

### Unified API

- `BinseqReader` enum (`src/parallel.rs`) — dispatches over Bq/Vbq/Cbq MmapReaders for reading via memory-mapped I/O.
- `BinseqWriter` enum + `BinseqWriterBuilder` (`src/write.rs`) — dispatches over Bq/Vbq/Cbq writers with a builder for configuration.

### Key traits and types

- `BinseqRecord` (`src/record/binseq_record.rs`) — trait for reading records; implemented by each format's RefRecord.
- `SequencingRecord` + `SequencingRecordBuilder` (`src/record/sequencing_record.rs`) — zero-copy record type for writing, uses borrowed references.
- `ParallelReader` / `ParallelProcessor` (`src/parallel.rs`) — traits for parallel range-based processing across threads.
- `Policy` (`src/policy.rs`) — how to handle invalid nucleotides (ignore, break, random draw, set to specific base).
- Error hierarchy (`src/error.rs`) — `thiserror`-based enums: HeaderError, ReadError, WriteError, CbqError, BuilderError, IndexError, ExtensionError.

### Dependencies of note

- `bitnuc` — nucleotide 2-bit and 4-bit encoding/decoding.
- `paraseq` — parallel FASTX file parsing (optional, enabled by default).
- `zstd` — block-level compression for VBQ and CBQ.
- `memmap2` — memory-mapped file reading.

## Conventions

- Rust edition 2024.
- Clippy pedantic enabled (`cast_possible_truncation` and `missing_errors_doc` allowed).
- Native CPU target set in `.cargo/config.toml` (`-C target-cpu=native`).
- All tests are inline `#[cfg(test)]` modules — no separate `tests/` directory.
- Default features: `paraseq` and `anyhow`. FASTX encoding utilities (`src/utils/fastx.rs`) require the `paraseq` feature.
