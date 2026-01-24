use std::{fs::File, io::BufWriter};

use anyhow::Result;
use binseq::{BinseqWriterBuilder, write::Format};
use bitnuc::BitSize;
use clap::Parser;

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
        let ofile = File::create(path).map(BufWriter::new)?;
        Ok(Box::new(ofile))
    }

    fn is_paired(&self) -> bool {
        self.input2.is_some()
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    let handle = args.ohandle()?;
    let builder = BinseqWriterBuilder::new(args.format())
        .bitsize(args.bitsize())
        .block_size(args.blocksize * 1024)
        .headers(!args.exclude_headers)
        .quality(!args.exclude_quality)
        .compression_level(args.compression_level)
        .encode_fastx(handle);
    if args.is_paired() {
        builder.input_paired(&args.input, args.input2.as_ref().unwrap())
    } else {
        builder.input(&args.input)
    }
    .threads(args.threads)
    .run()?;

    Ok(())
}
