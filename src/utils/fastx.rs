//! Parallel encoding of FASTX (FASTA/FASTQ) files into BINSEQ formats via `paraseq`.

use std::{
    io::{Read, Write},
    path::{Path, PathBuf},
    sync::Arc,
};

use paraseq::{
    Record, fastx,
    prelude::{IntoProcessError, PairedParallelProcessor, ParallelProcessor, ParallelReader},
};
use parking_lot::Mutex;

use crate::{
    BinseqWriter, BinseqWriterBuilder, IntoBinseqError, Result, SequencingRecordBuilder,
    error::WriteError,
};

type BoxedRead = Box<dyn Read + Send>;
type BoxedWrite = Box<dyn Write + Send>;

/// Input source for FASTX encoding
#[derive(Debug, Clone)]
enum FastxInput {
    /// Read from stdin
    Stdin,
    /// Read from a single file
    Single(PathBuf),
    /// Read from paired files (R1, R2)
    Paired(PathBuf, PathBuf),
}

/// Builder for encoding FASTX files to BINSEQ format
///
/// Created by [`BinseqWriterBuilder::encode_fastx`]; configures the input source
/// and threading before running the encoding.
///
/// Can be ordered or unordered;
/// though unordered takes better advantage of parallelism, ordering preserves
/// the input order in the fastx.
///
/// # Example
///
/// ```rust,no_run
/// use binseq::write::{BinseqWriterBuilder, Format};
/// use std::fs::File;
///
/// // Encode paired-end FASTQ to CBQ
/// BinseqWriterBuilder::new(Format::Cbq)
///     .quality(true)
///     .headers(true)
///     .encode_fastx(Box::new(File::create("output.cbq")?))
///     .input_paired("R1.fastq", "R2.fastq")
///     .threads(8)
///     .run()?;
/// # Ok::<(), binseq::Error>(())
/// ```
pub struct FastxEncoderBuilder {
    builder: BinseqWriterBuilder,
    output: BoxedWrite,
    input: Option<FastxInput>,
    ordered: bool,
    threads: usize,
}

impl FastxEncoderBuilder {
    /// Create a new encoder builder
    pub(crate) fn new(builder: BinseqWriterBuilder, output: BoxedWrite) -> Self {
        Self {
            builder,
            output,
            input: None,
            threads: 0,     // 0 means use all available cores
            ordered: false, // default to unordered for speed
        }
    }

    /// Read from a single FASTX file
    #[must_use]
    pub fn input<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.input = Some(FastxInput::Single(path.as_ref().to_path_buf()));
        self
    }

    /// Read from stdin
    #[must_use]
    pub fn input_stdin(mut self) -> Self {
        self.input = Some(FastxInput::Stdin);
        self
    }

    /// Read from paired FASTX files (R1, R2); sets the writer to paired mode
    #[must_use]
    pub fn input_paired<P: AsRef<Path>>(mut self, r1: P, r2: P) -> Self {
        self.input = Some(FastxInput::Paired(
            r1.as_ref().to_path_buf(),
            r2.as_ref().to_path_buf(),
        ));
        // Automatically set paired mode
        self.builder = self.builder.paired(true);
        self
    }

    /// Set the number of threads for parallel processing (0: all available cores)
    #[must_use]
    pub fn threads(mut self, n: usize) -> Self {
        self.threads = n;
        self
    }

    /// Set whether the output should be ordered as input (small perf cost, switch off for speed)
    #[must_use]
    pub fn ordered(mut self, ordered: bool) -> Self {
        self.ordered = ordered;
        self
    }

    /// Execute the FASTX encoding, consuming the builder
    pub fn run(mut self) -> Result<()> {
        let (r1, r2) = match self.input {
            Some(FastxInput::Single(path)) => {
                // build interleaved reader
                let mut reader =
                    fastx::Reader::from_path(path).map_err(IntoBinseqError::into_binseq_error)?;
                // Only probe for an extended length when the writer is configured
                // as paired (i.e. the single input is interleaved); otherwise a
                // second read's length would mark fixed-length formats as paired.
                let (slen, xlen) = detect_seq_len(&mut reader, self.builder.paired)?;
                self.builder = self.builder.slen(slen as u32).xlen(xlen as u32);
                (reader, None)
            }
            Some(FastxInput::Stdin) => {
                let mut reader =
                    fastx::Reader::from_stdin().map_err(IntoBinseqError::into_binseq_error)?;
                // Only probe for an extended length when the writer is configured
                // as paired (i.e. the single input is interleaved); otherwise a
                // second read's length would mark fixed-length formats as paired.
                let (slen, xlen) = detect_seq_len(&mut reader, self.builder.paired)?;
                self.builder = self.builder.slen(slen as u32).xlen(xlen as u32);
                (reader, None)
            }
            Some(FastxInput::Paired(path1, path2)) => {
                // build interleaved reader
                let mut reader1 =
                    fastx::Reader::from_path(path1).map_err(IntoBinseqError::into_binseq_error)?;
                let mut reader2 =
                    fastx::Reader::from_path(path2).map_err(IntoBinseqError::into_binseq_error)?;
                let (slen, _) = detect_seq_len(&mut reader1, false)?;
                let (xlen, _) = detect_seq_len(&mut reader2, false)?;
                self.builder = self.builder.slen(slen as u32).xlen(xlen as u32);
                (reader1, Some(reader2))
            }
            None => return Err(WriteError::MissingInput.into()),
        };

        let writer = self.builder.build(self.output)?;
        let paired = writer.is_paired();
        let mut encoder = Encoder::new(writer, self.ordered)?;
        match (paired, r2) {
            (true, Some(r2)) => r1.process_parallel_paired(r2, &mut encoder, self.threads),
            (true, None) => r1.process_parallel_interleaved(&mut encoder, self.threads),
            (false, _) => r1.process_parallel(&mut encoder, self.threads),
        }
        .map_err(IntoBinseqError::into_binseq_error)?;
        encoder.finish()?;

        Ok(())
    }
}

fn detect_seq_len(
    reader: &mut fastx::Reader<BoxedRead>,
    interleaved: bool,
) -> Result<(usize, usize)> {
    // Initialize the record set
    let mut rset = reader.new_record_set();
    rset.fill(reader)
        .map_err(IntoBinseqError::into_binseq_error)?;

    let (slen, xlen) = {
        let mut rset_iter = rset.iter();
        let mut next_len = || -> Result<usize> {
            let rec = rset_iter
                .next()
                .ok_or(WriteError::EmptyFastxFile)?
                .map_err(IntoBinseqError::into_binseq_error)?;
            Ok(rec.seq().len())
        };

        let slen = next_len()?;
        let xlen = if interleaved { next_len()? } else { 0 };
        (slen, xlen)
    };

    reader
        .reload(&mut rset)
        .map_err(IntoBinseqError::into_binseq_error)?;
    Ok((slen, xlen))
}

/// Parallel encoder implementing `paraseq`'s processor traits
#[derive(Clone)]
struct Encoder {
    /// Global writer (shared across threads)
    writer: Arc<Mutex<BinseqWriter<Box<dyn Write + Send>>>>,
    /// Thread-local writer buffer
    thread_writer: BinseqWriter<Vec<u8>>,
    /// Whether the output should follow same order as input
    ordered: bool,
}

impl Encoder {
    /// Create a new encoder with a global writer
    pub fn new(writer: BinseqWriter<Box<dyn Write + Send>>, ordered: bool) -> Result<Self> {
        let thread_writer = writer.new_headless_buffer()?;
        Ok(Self {
            writer: Arc::new(Mutex::new(writer)),
            thread_writer,
            ordered,
        })
    }
    /// Finish the stream on the global writer
    pub fn finish(&mut self) -> Result<()> {
        self.writer.lock().finish()?;
        Ok(())
    }
    fn flush_batch(&mut self) -> paraseq::Result<()> {
        let mut writer = self.writer.lock();
        if self.ordered {
            // can't take full advantage of parallelism if we have to order the output
            writer.ingest(&mut self.thread_writer)
        } else {
            writer.ingest_completed(&mut self.thread_writer)
        }
        .map_err(IntoProcessError::into_process_error)
    }
}

impl<Rf: Record> ParallelProcessor<Rf> for Encoder {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        let seq = record.seq();
        let seq_record = SequencingRecordBuilder::default()
            .s_header(record.id())
            .s_seq(&seq)
            .opt_s_qual(record.qual())
            .build()
            .map_err(IntoProcessError::into_process_error)?;
        self.thread_writer
            .push(seq_record)
            .map_err(IntoProcessError::into_process_error)?;
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        self.flush_batch()
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        self.writer
            .lock()
            .ingest(&mut self.thread_writer)
            .map_err(IntoProcessError::into_process_error)?;
        Ok(())
    }

    fn requires_ordering(&self) -> bool {
        self.ordered
    }
}

impl<Rf: Record> PairedParallelProcessor<Rf> for Encoder {
    fn process_record_pair(&mut self, record1: Rf, record2: Rf) -> paraseq::Result<()> {
        let sseq = record1.seq();
        let xseq = record2.seq();
        let seq_record = SequencingRecordBuilder::default()
            .s_header(record1.id())
            .s_seq(&sseq)
            .opt_s_qual(record1.qual())
            .x_header(record2.id())
            .x_seq(&xseq)
            .opt_x_qual(record2.qual())
            .build()
            .map_err(IntoProcessError::into_process_error)?;

        self.thread_writer
            .push(seq_record)
            .map_err(IntoProcessError::into_process_error)?;
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        self.flush_batch()
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        self.writer
            .lock()
            .ingest(&mut self.thread_writer)
            .map_err(IntoProcessError::into_process_error)?;
        Ok(())
    }

    fn requires_ordering(&self) -> bool {
        self.ordered
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::write::Format;
    use std::io::Cursor;

    const FASTQ_R1_PATH: &str = "./data/subset_R1.fastq.gz";
    const FASTQ_R2_PATH: &str = "./data/subset_R2.fastq.gz";

    #[test]
    fn test_encoder_builder_construction() {
        let builder = BinseqWriterBuilder::new(Format::Cbq);
        let handle = Box::new(Cursor::new(Vec::new()));
        let encoder_builder = FastxEncoderBuilder::new(builder, handle);

        assert!(encoder_builder.input.is_none());
        assert_eq!(encoder_builder.threads, 0);
        assert!(!encoder_builder.ordered);
    }

    #[test]
    fn test_encoder_builder_ordered_setter() {
        let builder = BinseqWriterBuilder::new(Format::Cbq);
        let handle = Box::new(Cursor::new(Vec::new()));
        let encoder_builder = FastxEncoderBuilder::new(builder, handle).ordered(false);

        assert!(!encoder_builder.ordered);
    }

    #[test]
    fn test_encoder_builder_input_methods() {
        let builder = BinseqWriterBuilder::new(Format::Cbq);
        let handle = Box::new(Cursor::new(Vec::new()));
        let encoder_builder = FastxEncoderBuilder::new(builder, handle)
            .input("test.fastq")
            .threads(4);

        assert!(matches!(encoder_builder.input, Some(FastxInput::Single(_))));
        assert_eq!(encoder_builder.threads, 4);
    }

    #[test]
    fn test_encoder_builder_stdin() {
        let builder = BinseqWriterBuilder::new(Format::Cbq);
        let handle = Box::new(Cursor::new(Vec::new()));
        let encoder_builder = FastxEncoderBuilder::new(builder, handle).input_stdin();

        assert!(matches!(encoder_builder.input, Some(FastxInput::Stdin)));
    }

    #[test]
    fn test_encoder_builder_single() {
        let builder = BinseqWriterBuilder::new(Format::Cbq);
        let handle = Box::new(Cursor::new(Vec::new()));
        let encoder_builder = FastxEncoderBuilder::new(builder, handle).input(FASTQ_R1_PATH);

        assert!(matches!(encoder_builder.input, Some(FastxInput::Single(_))));

        // Run the encoder builder and assert that it is successful
        assert!(encoder_builder.run().is_ok());
    }

    #[test]
    fn test_encoder_builder_paired() {
        let builder = BinseqWriterBuilder::new(Format::Cbq);
        let handle = Box::new(Cursor::new(Vec::new()));
        let encoder_builder =
            FastxEncoderBuilder::new(builder, handle).input_paired(FASTQ_R1_PATH, FASTQ_R2_PATH);

        assert!(matches!(
            encoder_builder.input,
            Some(FastxInput::Paired(_, _))
        ));
        // Should automatically set paired mode
        assert!(encoder_builder.builder.paired);

        // Run the encoder builder and assert that it is successful
        assert!(encoder_builder.run().is_ok());
    }

    #[derive(Clone, Default)]
    struct HeaderCollector {
        headers: Arc<parking_lot::Mutex<Vec<(u64, String)>>>,
    }
    impl crate::ParallelProcessor for HeaderCollector {
        fn process_record<R: crate::BinseqRecord>(&mut self, record: R) -> crate::Result<()> {
            let header = String::from_utf8_lossy(record.sheader()).into_owned();
            self.headers.lock().push((record.index(), header));
            Ok(())
        }
    }

    /// Encodes a synthetic multi-threaded FASTQ input with `ordered(true)` and confirms
    /// the written BINSEQ records come back out in the same order as the input, even
    /// though multiple threads raced to produce them.
    #[test]
    fn test_encoder_builder_ordered_preserves_record_order() {
        use crate::{BinseqReader, ParallelReader as DecodeReader};
        use std::fmt::Write as _;
        use std::sync::Arc;

        const N_RECORDS: usize = 4_000;
        const SEQ_LEN: usize = 32;

        let mut fastq = String::new();
        for i in 0..N_RECORDS {
            let base = b"ACGT"[i % 4] as char;
            let seq: String = std::iter::repeat_n(base, SEQ_LEN).collect();
            let qual: String = std::iter::repeat_n('F', SEQ_LEN).collect();
            let _ = writeln!(fastq, "@read_{i:06}\n{seq}\n+\n{qual}");
        }

        let pid = std::process::id();
        let input_path = std::env::temp_dir().join(format!("binseq_ordered_input_{pid}.fastq"));
        let output_path = std::env::temp_dir().join(format!("binseq_ordered_output_{pid}.cbq"));
        std::fs::write(&input_path, &fastq).unwrap();

        let builder = BinseqWriterBuilder::new(Format::Cbq).headers(true);
        let handle = Box::new(std::fs::File::create(&output_path).unwrap());
        let result = FastxEncoderBuilder::new(builder, handle)
            .input(&input_path)
            .threads(4)
            .ordered(true)
            .run();
        std::fs::remove_file(&input_path).unwrap();
        assert!(result.is_ok());

        let reader = BinseqReader::new(&output_path).unwrap();
        let processor = HeaderCollector::default();
        let headers = processor.headers.clone();
        reader.process_parallel(processor, 4).unwrap();
        std::fs::remove_file(&output_path).unwrap();

        let mut results = Arc::try_unwrap(headers).unwrap().into_inner();
        results.sort_by_key(|(idx, _)| *idx);

        assert_eq!(results.len(), N_RECORDS);
        for (i, (idx, header)) in results.iter().enumerate() {
            assert_eq!(*idx, i as u64);
            assert_eq!(header, &format!("read_{i:06}"));
        }
    }
}
