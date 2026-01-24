use std::{
    io::{BufWriter, Read},
    sync::Arc,
};

use anyhow::{Result, bail};
use binseq::{
    SequencingRecordBuilder,
    write::{BinseqWriter, BinseqWriterBuilder, Format},
};
use bitnuc::BitSize;
use clap::Parser;
use paraseq::{
    Record, fastx,
    prelude::{IntoProcessError, PairedParallelProcessor, ParallelProcessor, ParallelReader},
};
use parking_lot::Mutex;

type BoxedWriter = Box<dyn std::io::Write + Send>;

#[derive(Parser)]
struct Args {
    /// Input FASTX to encode into BINSEQ format
    #[clap(required = true)]
    input: String,

    /// Input FASTX to encode into BINSEQ format (R2)
    #[clap(required = false)]
    input2: Option<String>,

    /// Output file path for BINSEQ format
    #[clap(short = 'o', long)]
    output: Option<String>,

    /// Default prefix for writing BINSEQ: `<prefix>.<ext>`
    #[clap(short = 'p', long, default_value = "output")]
    prefix: String,

    /// Format of the output BINSEQ file
    ///
    /// [bq: bq|BQ|b, vbq: vbq|VBQ|v, cbq: cbq|CBQ|c]
    #[clap(short = 'f', long)]
    format: Option<Format>,

    /// Exclude quality information in BINSEQ output
    ///
    /// (bq ignores quality always)
    #[clap(short = 'Q', long)]
    exclude_quality: bool,

    /// Exclude sequence headers in BINSEQ output
    ///
    /// (bq ignores headers always)
    #[clap(short = 'H', long)]
    exclude_headers: bool,

    /// Compression level for BINSEQ output (0: auto)
    #[clap(long, default_value_t = 0)]
    compression_level: i32,

    /// Default BITSIZE for BINSEQ output (2: 2bit, 4: 4bit)
    #[clap(long, default_value_t = 2)]
    bitsize: u8,

    /// Default BLOCKSIZE in KB for BINSEQ output (vbq,cbq)
    #[clap(long, default_value_t = 128)]
    blocksize: usize,

    /// Number of threads to use for parallel processing, 0: all available
    #[clap(short = 'T', long, default_value = "0")]
    threads: usize,
}
impl Args {
    /// Determines the output format based on the file extension or the provided format
    fn format(&self) -> Format {
        if let Some(format) = self.format {
            format
        } else {
            if let Some(output) = &self.output {
                match output.split(".").last() {
                    Some("bq") => Format::Bq,
                    Some("vbq") => Format::Vbq,
                    Some("cbq") => Format::Cbq,
                    _ => Format::default(),
                }
            } else {
                Format::default()
            }
        }
    }
    fn bitsize(&self) -> BitSize {
        match self.bitsize {
            4 => BitSize::Four,
            _ => BitSize::Two,
        }
    }

    /// Creates an output file handle
    fn ohandle(&self) -> Result<BoxedWriter> {
        let path = if let Some(output) = &self.output {
            output.to_string()
        } else {
            format!("{}{}", &self.prefix, self.format().extension())
        };
        let ofile = std::fs::File::create(path).map(BufWriter::new)?;
        Ok(Box::new(ofile))
    }

    fn is_paired(&self) -> bool {
        self.input2.is_some()
    }
}

/// Calculates the sequence length of the first record in the reader
fn get_seq_len<R: Read>(reader: &mut fastx::Reader<R>) -> Result<usize> {
    let mut rset = reader.new_record_set();
    rset.fill(reader)?;

    let slen = if let Some(record) = rset.iter().next() {
        let record = record?;
        record.seq().len()
    } else {
        bail!("No records found in reader");
    };

    reader.reload(&mut rset)?;

    Ok(slen)
}

#[derive(Clone)]
struct Encoder {
    /// global writer
    writer: Arc<Mutex<BinseqWriter<BoxedWriter>>>,
    thread_writer: BinseqWriter<Vec<u8>>,
}
impl Encoder {
    pub fn new(writer: BinseqWriter<BoxedWriter>) -> Result<Self> {
        let thread_writer = writer.new_headless_buffer()?;
        Ok(Self {
            writer: Arc::new(Mutex::new(writer)),
            thread_writer,
        })
    }
    pub fn finish(&mut self) -> Result<()> {
        self.writer.lock().finish()?;
        Ok(())
    }
}
impl<Rf: paraseq::Record> ParallelProcessor<Rf> for Encoder {
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

impl<Rf: paraseq::Record> PairedParallelProcessor<Rf> for Encoder {
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

fn encode_paired(args: &Args) -> Result<()> {
    let mut r1 = fastx::Reader::from_path(&args.input)?;
    let mut r2 = fastx::Reader::from_path(&args.input2.as_ref().expect("Missing input2"))?;
    let ohandle = args.ohandle()?;

    // prepare writer
    let writer = {
        let format = args.format();
        let mut builder = BinseqWriterBuilder::new(format)
            .headers(!args.exclude_headers)
            .quality(!args.exclude_quality)
            .compression_level(args.compression_level)
            .bitsize(args.bitsize())
            .paired(true)
            .block_size(args.blocksize * 1024);

        // BQ requires a fixed sequence length from init time
        if matches!(format, Format::Bq) {
            builder = builder.slen(get_seq_len(&mut r1)? as u32);
            builder = builder.xlen(get_seq_len(&mut r2)? as u32);
        }

        builder.build(ohandle)?
    };

    let mut encoder = Encoder::new(writer)?;
    r1.process_parallel_paired(r2, &mut encoder, args.threads)?;
    encoder.finish()?;

    Ok(())
}

fn encode_single(args: &Args) -> Result<()> {
    let mut reader = fastx::Reader::from_path(&args.input)?;
    let ohandle = args.ohandle()?;

    // prepare writer
    let writer = {
        let format = args.format();
        let mut builder = BinseqWriterBuilder::new(format)
            .headers(!args.exclude_headers)
            .quality(!args.exclude_quality)
            .compression_level(args.compression_level)
            .bitsize(args.bitsize())
            .block_size(args.blocksize * 1024);

        // BQ requires a fixed sequence length from init time
        if matches!(format, Format::Bq) {
            builder = builder.slen(get_seq_len(&mut reader)? as u32);
        }

        builder.build(ohandle)?
    };

    let mut encoder = Encoder::new(writer)?;
    reader.process_parallel(&mut encoder, args.threads)?;
    encoder.finish()?;

    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();
    if args.is_paired() {
        encode_paired(&args)
    } else {
        encode_single(&args)
    }
}
