# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.9.6] - 2026-08-12

### Changed

- `cbq` uses the `bitnuc` implementation with `u8`-based encoding with better SIMD perf
  - Does **not** change on-disk representation - block reading/writing is still done with u64 byte-boundaries and the format is unchanged.
- `cbq` blocks now use `bitnuc::ambiguous_bases` to find noncanonical bases instead of just `N` and adds them to the `npos` index.
- dependency both on `bitnuc-0.4` and `bitnuc-0.5` to keep 4-bit support for `bq` and `vbq` which is unsupported in `bitnuc-0.5.x`.

## [0.9.5] - 2026-08-10

### Fixed

- Support empty CBQ files in the `cbq` reader which previously failed due to index alignment error

## [0.9.4] - 2026-07-15

### Fixed

- Bug in the `bq` stream reader that used an incorrect record id for small buffers.

### Added

- `BinseqReader::new` now determines a file's BINSEQ format (BQ, VBQ, or CBQ) by sniffing its
  magic bytes instead of relying on the file extension, so it works regardless of how the file
  is named.
- `Format::sniff` for identifying a BINSEQ format from a byte buffer, plus public `FILE_MAGIC`
  constants on `bq`, `vbq`, and `cbq`.

### Changed

- Added clippy checks and additional lint allowances to CI, plus general style fixes.
- Improved test coverage throughout the library.
- Renamed `ExtensionError` to `FormatError` (`UnrecognizedMagicBytes` variant) to reflect
  detection now being based on file content rather than the file extension.

## [0.9.3] - 2026-07-01

### Changed

- `Binseq` input is now accepted as a path rather than a pre-opened handle.

### Fixed

- Use `reserve` instead of `resize` to avoid unnecessary zeroing/copies (issue #94).
- Ensure `len_nef` is cleared correctly.

### Tests

- Added coverage for the sequential reader and N-containing sequences.

## [0.9.1] - 2026-06-25

### Added

- Support for writing only completed blocks, enabling on-thread compression and
  reducing lock contention for parallel writers.
- `on_thread_complete` callback support.
- `ingest_completed` is now passed through, with behavior changes scoped to `cbq`.

### Tests

- Added coverage for `cbq` ingestion.

## [0.9.0] - 2026-01-23

### Added

- Integrated `cbq` (compressed BINSEQ) into the core library as a first-class format.
- Builder-based construction of a file from an existing header.
- Auto-encode a FASTX file (or a pair of FASTX files) directly into a BINSEQ format via
  a builder method.
- A generic writer abstraction over BINSEQ files.
- Default quality score handling across all readers.
- Elias-Fano encoding for N-positions instead of raw `u64` indices, reducing index size.
- An iterator over block headers, and convenience accessors for index metadata.

### Changed

- **Breaking:** renamed `BINSEQ`/`VBINSEQ` naming throughout to `BQ`/`VBQ`
  (e.g. `BinseqHeader` → `bq::FileHeader`), and renamed several `vbq`/`bq`-specific
  headers and writer structs for clarity.
- Reworked the write API around a sequence record struct.
- Reworked the internal index to track record count and sequence count separately.
- Replaced `anyhow`-based errors throughout the crate with crate-native error types.
- Removed the external `vbq` index file in favor of writing the index into the file itself.

### Fixed

- Several small hotfixes to the write path.
- Incorrect binary file used in a test fixture.

## [0.8.1] - 2025-12-10

### Fixed

- Incorrect decoding of the extended buffer on batch decoding in `bq` when a full
  chunk of primary data was present (issue #75).

## [0.8.0] - 2025-12-10

### Added

- `sseq`/`xseq` functions on the record trait for use during batch decoding, implemented
  for both `bq` and `vbq`.
- Additional usage examples.

### Changed

- **Breaking:** improved `bq` decoding via batch decoding, keeping data in SIMD form longer.
- Improved header access logic for zero-copy reads.

### Fixed

- Correct index returned from decoding.
- Correct slicing when per-record flags are set in `bq`.

## [0.7.8] - 2025-12-09

### Changed

- `vbq` reader is now zero-copy or single-copy depending on compression status.
- Simplified the encoding API.
- Reusable decompression context for improved throughput.
- Updated all dependencies.

### Fixed

- Validate the size of bytes provided to the non-compressed reader.

## [0.7.7] - 2025-12-08

### Changed

- Split the sequencing context into multiple traits with distinct variants
  (renamed `FullCtx` to `Ctx`).
- Removed unsafe code from the context implementation in favor of a resize-copy pattern.

### Fixed

- Logical error in quality-score padding (`ensure_quality`).

### Docs

- Documented the different context variants.

## [0.7.6] - 2025-12-08

### Added

- A basic sequencing context struct, with mutable references to its internals.
- Ability to fill in missing quality scores when required.
- A generic error-conversion trait, made public.

## [0.7.5] - 2025-10-02

### Added

- BINSEQ v2 format.
- Optional headers in `vbq` files, with the index moved into the file itself
  (no more auxiliary index file).
- Support for 4-bit encodings, with a `bitsize` re-export so downstream users don't
  need to depend on `bitnuc` directly.
- Header builder pattern for constructing `bq` and `vbq` files.

### Changed

- Per-record flags are now optional across `bq`, `vbq`, and the record API.
- `vbq` index now uses `u64` for cumulative record counts.

## [0.6.5] - 2025-07-01

Baseline release for this changelog. Earlier history predates version tagging.

[Unreleased]: https://github.com/arcinstitute/binseq/compare/v0.9.3...HEAD
[0.9.3]: https://github.com/arcinstitute/binseq/compare/v0.9.1...v0.9.3
[0.9.1]: https://github.com/arcinstitute/binseq/compare/v0.9.0...v0.9.1
[0.9.0]: https://github.com/arcinstitute/binseq/compare/v0.8.1...v0.9.0
[0.8.1]: https://github.com/arcinstitute/binseq/compare/v0.8.0...v0.8.1
[0.8.0]: https://github.com/arcinstitute/binseq/compare/v0.7.8...v0.8.0
[0.7.8]: https://github.com/arcinstitute/binseq/compare/v0.7.7...v0.7.8
[0.7.7]: https://github.com/arcinstitute/binseq/compare/v0.7.6...v0.7.7
[0.7.6]: https://github.com/arcinstitute/binseq/compare/v0.7.5...v0.7.6
[0.7.5]: https://github.com/arcinstitute/binseq/compare/v0.6.5...v0.7.5
[0.6.5]: https://github.com/arcinstitute/binseq/releases/tag/v0.6.5
