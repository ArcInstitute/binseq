use std::io::{BufReader, Cursor};

use binseq::bq::{FileHeaderBuilder, StreamReader, StreamWriterBuilder};
use binseq::{BinseqRecord, Policy, Result};

fn main() -> Result<()> {
    // Create a header for sequences of length 100
    let header = FileHeaderBuilder::new().slen(100).build()?;

    // Create some example sequence data
    let sequence = b"ACGT".repeat(25); // 100 nucleotides

    // Create a stream writer with a memory buffer as destination
    let mut writer = StreamWriterBuilder::default()
        .header(header)
        .policy(Policy::RandomDraw) // Use random nucleotides for invalid bases
        .buffer_capacity(4096) // Use 4K buffer
        .build(Cursor::new(Vec::new()))?;

    // Write the sequence with flag 0
    #[allow(deprecated)]
    writer.write_record(Some(0), &sequence)?;

    // Write the sequence with flag 1
    #[allow(deprecated)]
    writer.write_record(Some(1), &sequence)?;

    // Flush and get the buffer
    let buffer = writer.into_inner()?;
    let buffer_inner = buffer.into_inner();

    println!("Wrote {} bytes to buffer", buffer_inner.len());

    // Now read from the buffer using the streaming reader
    let cursor = Cursor::new(buffer_inner);
    let buf_reader = BufReader::new(cursor);

    // Create a stream reader
    let mut reader = StreamReader::new(buf_reader);

    // Read and display the header
    let header = reader.read_header()?;
    println!("Read header: sequence length = {}", header.slen);

    // Read records one by one
    let mut count = 0;
    while let Some(record) = reader.next_record() {
        let record = record?;
        println!("Record {}: flag = {:?}", count, record.flag());
        count += 1;
    }

    println!("Read {count} records");

    Ok(())
}
