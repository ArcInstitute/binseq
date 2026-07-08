//! Internal zstd backend selection.
//!
//! `binseq` supports two interchangeable zstd backends, selected at compile time by
//! Cargo feature:
//!
//! - `zstd-c` (default): the [`zstd`] crate, which links the C reference
//!   implementation via `zstd-sys`/`cc`.
//! - `zstd-rs`: a hand-written wrapper around Trifecta Tech Foundation's
//!   [`libzstd-rs-sys`](https://crates.io/crates/libzstd-rs-sys), a c2rust
//!   transpilation of the same library with no C compiler dependency.
//!
//! Both backends expose an identical surface (`zstd_safe::{CCtx, DCtx, CParameter,
//! compress_bound, get_error_name}`, `stream::{copy_encode, copy_decode}`, and
//! `{Encoder, Decoder}`), so call sites elsewhere in the crate are backend-agnostic:
//! they just `use crate::zstdlib::...` instead of `use zstd::...`.
//!
//! Enable exactly one of the two features - see the individual backend modules for
//! implementation notes.

#[cfg(all(feature = "zstd-c", feature = "zstd-rs"))]
compile_error!("features `zstd-c` and `zstd-rs` are mutually exclusive - pick one zstd backend");
#[cfg(not(any(feature = "zstd-c", feature = "zstd-rs")))]
compile_error!("one of `zstd-c` or `zstd-rs` must be enabled to select a zstd backend");

#[cfg(feature = "zstd-c")]
mod backend_c;
#[cfg(feature = "zstd-c")]
pub(crate) use backend_c::*;

#[cfg(feature = "zstd-rs")]
mod backend_rs;
#[cfg(feature = "zstd-rs")]
pub(crate) use backend_rs::*;

/// Cross-backend interop guarantee: `zstd-rs` (Trifecta's `libzstd-rs-sys`) is a
/// pre-1.0 spike backend. The whole point of offering it as a drop-in alternative to
/// `zstd-c` is that files aren't tied to whichever backend wrote them - both produce
/// standard zstd frames, so a file written by one build must be readable by the
/// other. That guarantee is what this module pins down, so a future change to either
/// backend that breaks the bitstream (not just self-consistency) gets caught.
///
/// `zstd-c` and `zstd-rs` are mutually exclusive Cargo features, so a single test
/// binary can never link both at once to compress-then-decompress in-process.
/// Instead, this decodes two *fixture files* - one written by each backend ahead of
/// time - with whichever backend the current build uses, and asserts they decode to
/// identical records. Run under both feature sets (this repo's CI does: the default
/// `build` job covers `zstd-c`, `build_zstd_rs` covers `zstd-rs`), that covers all
/// four cells of the write x read backend matrix without ever needing both backends
/// in one binary.
///
/// Regenerate the fixtures (only needed if the on-disk format changes) with:
/// ```sh
/// cargo run --release --example write -- data/subset_R1.fastq.gz -o data/interop_zstd_c.cbq -T 1
/// cargo run --release --example write -- data/subset_R1.fastq.gz -o data/interop_zstd_c.vbq -T 1
/// cargo run --release --no-default-features --features paraseq,anyhow,zstd-rs \
///     --example write -- data/subset_R1.fastq.gz -o data/interop_zstd_rs.cbq -T 1
/// cargo run --release --no-default-features --features paraseq,anyhow,zstd-rs \
///     --example write -- data/subset_R1.fastq.gz -o data/interop_zstd_rs.vbq -T 1
/// ```
#[cfg(test)]
mod tests {
    use std::sync::{Arc, Mutex};

    use crate::{BinseqReader, BinseqRecord, ParallelProcessor, ParallelReader, Result};

    type RecordTuple = (Vec<u8>, Vec<u8>, Vec<u8>);

    #[derive(Clone, Default)]
    struct Collector {
        records: Arc<Mutex<Vec<RecordTuple>>>,
    }

    impl ParallelProcessor for Collector {
        fn process_record<R: BinseqRecord>(&mut self, record: R) -> Result<()> {
            self.records.lock().unwrap().push((
                record.sheader().to_vec(),
                record.sseq().to_vec(),
                record.squal().to_vec(),
            ));
            Ok(())
        }
    }

    /// Decodes every record in `path`, sorted for order-independent comparison
    /// (decoding is parallelized, so record order isn't guaranteed to match between
    /// two reads).
    fn decode_all(path: &str) -> Vec<RecordTuple> {
        let reader = BinseqReader::new(path).expect("failed to open fixture");
        let collector = Collector::default();
        reader
            .process_parallel(collector.clone(), 0)
            .expect("failed to decode fixture");
        let mut records = Arc::try_unwrap(collector.records)
            .expect("no other references")
            .into_inner()
            .unwrap();
        records.sort();
        records
    }

    #[test]
    fn cbq_fixtures_decode_identically_across_backends() {
        let from_c = decode_all("data/interop_zstd_c.cbq");
        let from_rs = decode_all("data/interop_zstd_rs.cbq");
        assert!(!from_c.is_empty(), "fixture must contain records");
        assert_eq!(
            from_c, from_rs,
            "a file written by one zstd backend must decode identically under the other"
        );
    }

    #[test]
    fn vbq_fixtures_decode_identically_across_backends() {
        let from_c = decode_all("data/interop_zstd_c.vbq");
        let from_rs = decode_all("data/interop_zstd_rs.vbq");
        assert!(!from_c.is_empty(), "fixture must contain records");
        assert_eq!(
            from_c, from_rs,
            "a file written by one zstd backend must decode identically under the other"
        );
    }
}
