use anyhow::{bail, Result};
use seq_io::fastq::{Reader, Record};
use std::{
    fs::File,
    io::{BufReader, BufWriter},
};

use binseq::{
    BinseqHeader, BinseqRead, BinseqRecord, BinseqWriter, PairedRead, PairedReader, SingleReader,
};

fn read_write_single(fastq_path: &str, binseq_path: &str, seq_size: usize) -> Result<()> {
    // Open the input FASTQ file
    let (in_handle, _comp) = niffler::from_path(fastq_path)?;

    // Open the output file
    let header = BinseqHeader::new(seq_size as u32);
    let out_handle = File::create(binseq_path).map(BufWriter::new)?;
    let mut writer = BinseqWriter::new(out_handle, header, false)?;

    let mut all_sequences = Vec::new();

    // Write the binary sequence
    let mut reader = Reader::new(in_handle);
    let mut num_records_write = 0;
    let mut skipped_records = 0;
    while let Some(record) = reader.next() {
        let record = record?;
        let seq = record.seq();
        if writer.write_nucleotides(0, seq)? {
            num_records_write += 1;
            all_sequences.push(seq.to_vec());
        } else {
            skipped_records += 1;
        }
    }
    writer.flush()?;
    eprintln!(
        "Finished writing {} records to path: {}",
        num_records_write, binseq_path
    );
    eprintln!("Skipped {} records", skipped_records);

    // Read the binary sequence
    let bufreader = File::open(binseq_path).map(BufReader::new)?;
    let mut reader = SingleReader::new(bufreader)?;
    let mut num_records_read = 0;
    let mut record_buffer = Vec::new();
    while let Some(record) = reader.next() {
        let record = record?;

        record.decode(&mut record_buffer)?;

        // Check if the decoded sequence matches the original
        let buf_str = std::str::from_utf8(&record_buffer)?;
        let seq_str = std::str::from_utf8(&all_sequences[num_records_read])?;
        assert_eq!(buf_str, seq_str);

        num_records_read += 1;
        record_buffer.clear();
    }
    eprintln!("Finished reading {} records", num_records_read);
    eprintln!(
        "Difference in total records: {}",
        num_records_write - num_records_read
    );
    eprintln!("Number of records in vec: {}", all_sequences.len());

    Ok(())
}

fn read_write_paired(
    fastq_path_r1: &str,
    fastq_path_r2: &str,
    binseq_path: &str,
    seq_size_r1: usize,
    seq_size_r2: usize,
) -> Result<()> {
    // Open the input FASTQ files

    let in_buf_r1 = File::open(fastq_path_r1).map(BufReader::new)?;
    let in_buf_r2 = File::open(fastq_path_r2).map(BufReader::new)?;

    let (in_handle_r1, _comp) = niffler::get_reader(Box::new(in_buf_r1))?;
    let (in_handle_r2, _comp) = niffler::get_reader(Box::new(in_buf_r2))?;

    // Create the header
    let header = BinseqHeader::new_extended(seq_size_r1 as u32, seq_size_r2 as u32);

    // Open the output handle
    let out_handle = File::create(binseq_path).map(BufWriter::new)?;

    // Create the writer
    let mut writer = BinseqWriter::new(out_handle, header, false)?;

    // Open the FASTQ readers
    let mut reader_r1 = Reader::new(in_handle_r1);
    let mut reader_r2 = Reader::new(in_handle_r2);

    // Write the binary sequence
    let mut num_records = 0;
    let mut num_skipped = 0;

    let mut r1_storage = Vec::new();
    let mut r2_storage = Vec::new();

    loop {
        let (record_r1, record_r2) = match (reader_r1.next(), reader_r2.next()) {
            (Some(r1), Some(r2)) => (r1?, r2?),
            (None, None) => break,
            _ => bail!("Mismatched number of records in R1 and R2"),
        };

        let seq_r1 = record_r1.seq();
        let seq_r2 = record_r2.seq();

        if writer.write_paired(0, seq_r1, seq_r2)? {
            num_records += 1;
            r1_storage.push(seq_r1.to_vec());
            r2_storage.push(seq_r2.to_vec());
        } else {
            num_skipped += 1;
        }
    }
    writer.flush()?;
    eprintln!("Finished writing {} records", num_records);
    eprintln!("Skipped {} records", num_skipped);

    // Read the binary sequence

    let bufreader = File::open(binseq_path).map(BufReader::new)?;
    let mut reader = PairedReader::new(bufreader)?;
    let mut sbuf = Vec::new();
    let mut xbuf = Vec::new();

    let mut n_processed = 0;
    while let Some(pair) = reader.next_paired() {
        let pair = pair?;
        pair.decode_s(&mut sbuf)?;
        pair.decode_x(&mut xbuf)?;

        // Check if the decoded sequence matches the original
        let s_str = std::str::from_utf8(&sbuf)?;
        let x_str = std::str::from_utf8(&xbuf)?;

        let s_exp = std::str::from_utf8(&r1_storage[n_processed])?;
        let x_exp = std::str::from_utf8(&r2_storage[n_processed])?;

        assert_eq!(s_str, s_exp);
        assert_eq!(x_str, x_exp);

        sbuf.clear();
        xbuf.clear();

        n_processed += 1;
    }
    eprintln!("Finished reading {} records", n_processed);

    Ok(())
}

fn main() -> Result<()> {
    // INPUT ARGUMENTS
    let fastq_path_r1 = "./data/subset_R1.fastq.gz"; // exists
    let fastq_path_r2 = "./data/subset_R2.fastq.gz"; // exists
    let binseq_path_r1 = "./data/subset_R1.bq"; // created
    let binseq_path_r2 = "./data/subset_R2.bq"; // created
    let binseq_path = "./data/subset.bq"; // created
    let seq_size_r1 = 28; // a priori known
    let seq_size_r2 = 90; // a priori known

    read_write_single(fastq_path_r1, binseq_path_r1, seq_size_r1)?;
    read_write_single(fastq_path_r2, binseq_path_r2, seq_size_r2)?;
    read_write_paired(
        fastq_path_r1,
        fastq_path_r2,
        binseq_path,
        seq_size_r1,
        seq_size_r2,
    )?;

    Ok(())
}
