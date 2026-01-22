//! FASTX encoding utilities for converting FASTX files to BINSEQ formats
//!
//! This module provides utilities for encoding FASTX (FASTA/FASTQ) files into
//! BINSEQ formats using parallel processing via the `paraseq` crate.

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
    error::FastxEncodingError,
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
/// This builder is created by calling [`BinseqWriterBuilder::encode_fastx`] and
/// provides a fluent interface for configuring the input source and threading options.
///
/// # Example
///
/// ```rust,no_run
/// use binseq::write::{BinseqWriterBuilder, Format};
/// use std::fs::File;
///
/// // Encode from stdin to VBQ
/// let writer = BinseqWriterBuilder::new(Format::Vbq)
///     .quality(true)
///     .headers(true)
///     .encode_fastx(Box::new(File::create("output.vbq")?))
///     .input_stdin()
///     .threads(8)
///     .run()?;
/// # Ok::<(), binseq::Error>(())
/// ```
pub struct FastxEncoderBuilder {
    builder: BinseqWriterBuilder,
    output: BoxedWrite,
    input: Option<FastxInput>,
    threads: usize,
}

impl FastxEncoderBuilder {
    /// Create a new encoder builder
    pub(crate) fn new(builder: BinseqWriterBuilder, output: BoxedWrite) -> Self {
        Self {
            builder,
            output,
            input: None,
            threads: 0, // 0 means use all available cores
        }
    }

    /// Read from a single FASTX file
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// # use binseq::write::{BinseqWriterBuilder, Format};
    /// # use std::fs::File;
    /// BinseqWriterBuilder::new(Format::Vbq)
    ///     .encode_fastx(Box::new(File::create("output.vbq")?))
    ///     .input("input.fastq")
    ///     .run()?;
    /// # Ok::<(), binseq::Error>(())
    /// ```
    pub fn input<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.input = Some(FastxInput::Single(path.as_ref().to_path_buf()));
        self
    }

    /// Read from stdin
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// # use binseq::write::{BinseqWriterBuilder, Format};
    /// # use std::fs::File;
    /// BinseqWriterBuilder::new(Format::Vbq)
    ///     .encode_fastx(Box::new(File::create("output.vbq")?))
    ///     .input_stdin()
    ///     .run()?;
    /// # Ok::<(), binseq::Error>(())
    /// ```
    #[must_use]
    pub fn input_stdin(mut self) -> Self {
        self.input = Some(FastxInput::Stdin);
        self
    }

    /// Read from paired FASTX files (R1, R2)
    ///
    /// This automatically sets the writer to paired mode.
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// # use binseq::write::{BinseqWriterBuilder, Format};
    /// # use std::fs::File;
    /// BinseqWriterBuilder::new(Format::Vbq)
    ///     .encode_fastx(Box::new(File::create("output.vbq")?))
    ///     .input_paired("R1.fastq", "R2.fastq")
    ///     .run()?;
    /// # Ok::<(), binseq::Error>(())
    /// ```
    pub fn input_paired<P: AsRef<Path>>(mut self, r1: P, r2: P) -> Self {
        self.input = Some(FastxInput::Paired(
            r1.as_ref().to_path_buf(),
            r2.as_ref().to_path_buf(),
        ));
        // Automatically set paired mode
        self.builder = self.builder.paired(true);
        self
    }

    /// Set the number of threads for parallel processing
    ///
    /// If not set or set to 0, uses all available CPU cores.
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// # use binseq::write::{BinseqWriterBuilder, Format};
    /// # use std::fs::File;
    /// BinseqWriterBuilder::new(Format::Vbq)
    ///     .encode_fastx(Box::new(File::create("output.vbq")?))
    ///     .input("input.fastq")
    ///     .threads(8)
    ///     .run()?;
    /// # Ok::<(), binseq::Error>(())
    /// ```
    #[must_use]
    pub fn threads(mut self, n: usize) -> Self {
        self.threads = n;
        self
    }

    /// Execute the FASTX encoding
    ///
    /// This consumes the builder and returns a `BinseqWriter` that has been
    /// populated with all records from the input FASTX file(s).
    ///
    /// # Errors
    ///
    /// Returns an error if:
    /// - The input files cannot be read
    /// - The FASTX format is invalid
    /// - The writer configuration is incompatible with the input
    /// - For BQ format with stdin input (cannot detect sequence length)
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// # use binseq::write::{BinseqWriterBuilder, Format};
    /// # use std::fs::File;
    /// let writer = BinseqWriterBuilder::new(Format::Vbq)
    ///     .encode_fastx(Box::new(File::create("output.vbq")?))
    ///     .input("input.fastq")
    ///     .run()?;
    /// # Ok::<(), binseq::Error>(())
    /// ```
    pub fn run(mut self) -> Result<()> {
        let (r1, r2) = match self.input {
            Some(FastxInput::Single(path)) => {
                // build interleaved reader
                let mut reader =
                    fastx::Reader::from_path(path).map_err(IntoBinseqError::into_binseq_error)?;
                let (slen, xlen) = detect_seq_len(&mut reader, true)?;
                self.builder = self.builder.slen(slen as u32).xlen(xlen as u32);
                (reader, None)
            }
            Some(FastxInput::Stdin) => {
                let mut reader =
                    fastx::Reader::from_stdin().map_err(IntoBinseqError::into_binseq_error)?;
                let (slen, xlen) = detect_seq_len(&mut reader, true)?;
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
            None => return Err(FastxEncodingError::MissingInput.into()),
        };

        let writer = self.builder.build(self.output)?;
        if writer.is_paired() {
            if let Some(r2) = r2 {
                encode_paired(writer, r1, r2, self.threads)?;
            } else {
                encode_interleaved(writer, r1, self.threads)?;
            }
        } else {
            encode_single_file(writer, r1, self.threads)?;
        }

        Ok(())
    }
}

/// Encode single-end reads from a file
fn encode_single_file(
    writer: BinseqWriter<BoxedWrite>,
    reader: fastx::Reader<BoxedRead>,
    threads: usize,
) -> Result<()> {
    let mut encoder = Encoder::new(writer)?;
    reader
        .process_parallel(&mut encoder, threads)
        .map_err(IntoBinseqError::into_binseq_error)?;
    encoder.finish()?;
    Ok(())
}

/// Encode paired-end reads from interleaved file
fn encode_interleaved(
    writer: BinseqWriter<BoxedWrite>,
    reader: fastx::Reader<BoxedRead>,
    threads: usize,
) -> Result<()> {
    let mut encoder = Encoder::new(writer)?;
    reader
        .process_parallel_interleaved(&mut encoder, threads)
        .map_err(IntoBinseqError::into_binseq_error)?;
    encoder.finish()?;
    Ok(())
}

/// Encode paired-end reads from files
fn encode_paired(
    writer: BinseqWriter<BoxedWrite>,
    r1: fastx::Reader<BoxedRead>,
    r2: fastx::Reader<BoxedRead>,
    threads: usize,
) -> Result<()> {
    let mut encoder = Encoder::new(writer)?;
    r1.process_parallel_paired(r2, &mut encoder, threads)
        .map_err(IntoBinseqError::into_binseq_error)?;
    encoder.finish()?;
    Ok(())
}

fn detect_seq_len(
    reader: &mut fastx::Reader<BoxedRead>,
    interleaved: bool,
) -> Result<(usize, usize)> {
    // Initialze the record set
    let mut rset = reader.new_record_set();
    rset.fill(reader)
        .map_err(IntoBinseqError::into_binseq_error)?;

    let (slen, xlen) = if interleaved {
        let mut rset_iter = rset.iter();
        let Some(Ok(slen)) = rset_iter.next().map(|r| -> Result<usize> {
            let rec = r.map_err(IntoBinseqError::into_binseq_error)?;
            Ok(rec.seq().len())
        }) else {
            return Err(FastxEncodingError::EmptyFastxFile.into());
        };
        let Some(Ok(xlen)) = rset_iter.next().map(|r| -> Result<usize> {
            let rec = r.map_err(IntoBinseqError::into_binseq_error)?;
            Ok(rec.seq().len())
        }) else {
            return Err(FastxEncodingError::EmptyFastxFile.into());
        };
        (slen, xlen)
    } else {
        let mut rset_iter = rset.iter();
        let Some(Ok(slen)) = rset_iter.next().map(|r| -> Result<usize> {
            let rec = r.map_err(IntoBinseqError::into_binseq_error)?;
            Ok(rec.seq().len())
        }) else {
            return Err(FastxEncodingError::EmptyFastxFile.into());
        };
        (slen, 0)
    };
    reader
        .reload(&mut rset)
        .map_err(IntoBinseqError::into_binseq_error)?;
    Ok((slen, xlen))
}

/// Parallel encoder for FASTX records to BINSEQ format
///
/// This struct implements the `ParallelProcessor` and `PairedParallelProcessor`
/// traits from `paraseq` to enable efficient parallel encoding of FASTX files.
#[derive(Clone)]
struct Encoder {
    /// Global writer (shared across threads)
    writer: Arc<Mutex<BinseqWriter<Box<dyn Write + Send>>>>,
    /// Thread-local writer buffer
    thread_writer: BinseqWriter<Vec<u8>>,
}

impl Encoder {
    /// Create a new encoder with a global writer
    pub fn new(writer: BinseqWriter<Box<dyn Write + Send>>) -> Result<Self> {
        let thread_writer = writer.new_headless_buffer()?;
        Ok(Self {
            writer: Arc::new(Mutex::new(writer)),
            thread_writer,
        })
    }
    /// Finish the stream on the global writer
    pub fn finish(&mut self) -> Result<()> {
        self.writer.lock().finish()?;
        Ok(())
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
        self.writer
            .lock()
            .ingest(&mut self.thread_writer)
            .map_err(IntoProcessError::into_process_error)?;
        Ok(())
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
        self.writer
            .lock()
            .ingest(&mut self.thread_writer)
            .map_err(IntoProcessError::into_process_error)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::write::Format;
    use std::io::Cursor;

    #[test]
    fn test_encoder_builder_construction() {
        let builder = BinseqWriterBuilder::new(Format::Vbq);
        let handle = Box::new(Cursor::new(Vec::new()));
        let encoder_builder = FastxEncoderBuilder::new(builder, handle);

        assert!(encoder_builder.input.is_none());
        assert_eq!(encoder_builder.threads, 0);
    }

    #[test]
    fn test_encoder_builder_input_methods() {
        let builder = BinseqWriterBuilder::new(Format::Vbq);
        let handle = Box::new(Cursor::new(Vec::new()));
        let encoder_builder = FastxEncoderBuilder::new(builder, handle)
            .input("test.fastq")
            .threads(4);

        assert!(matches!(encoder_builder.input, Some(FastxInput::Single(_))));
        assert_eq!(encoder_builder.threads, 4);
    }

    #[test]
    fn test_encoder_builder_stdin() {
        let builder = BinseqWriterBuilder::new(Format::Vbq);
        let handle = Box::new(Cursor::new(Vec::new()));
        let encoder_builder = FastxEncoderBuilder::new(builder, handle).input_stdin();

        assert!(matches!(encoder_builder.input, Some(FastxInput::Stdin)));
    }

    #[test]
    fn test_encoder_builder_paired() {
        let builder = BinseqWriterBuilder::new(Format::Vbq);
        let handle = Box::new(Cursor::new(Vec::new()));
        let encoder_builder =
            FastxEncoderBuilder::new(builder, handle).input_paired("r1.fastq", "r2.fastq");

        assert!(matches!(
            encoder_builder.input,
            Some(FastxInput::Paired(_, _))
        ));
        // Should automatically set paired mode
        assert!(encoder_builder.builder.paired);
    }
}
