use std::io::{BufReader, BufWriter};
use std::net::{TcpListener, TcpStream};
use std::thread;

use binseq::bq::{BinseqHeader, BinseqHeaderBuilder, StreamReader, StreamWriterBuilder};
use binseq::{BinseqRecord, Policy, Result};

fn server(header: BinseqHeader, sequence: &[u8]) -> Result<()> {
    // Create a listener on localhost:3000
    let listener = TcpListener::bind("127.0.0.1:3000").expect("Failed to bind to address");
    println!("Server listening on 127.0.0.1:3000");

    // Accept one connection
    let (stream, _) = listener.accept().expect("Failed to accept connection");
    println!("Client connected");

    let stream = BufWriter::new(stream);

    // Create a stream writer with the network stream as destination
    let mut writer = StreamWriterBuilder::default()
        .header(header)
        .policy(Policy::RandomDraw)
        .buffer_capacity(16384) // Larger buffer for network I/O
        .build(stream)?;

    // Write sequences in a loop
    for i in 0..10 {
        #[allow(deprecated)]
        writer.write_record(Some(i), sequence)?;
        println!("Server: Sent record {i}");

        // Simulate delay between records
        thread::sleep(std::time::Duration::from_millis(100));
    }

    // Ensure flush on drop
    writer.flush()?;
    println!("Server: All records sent");

    Ok(())
}

fn client() -> Result<()> {
    // Wait a moment for the server to start
    thread::sleep(std::time::Duration::from_millis(500));

    // Connect to the server
    let stream = TcpStream::connect("127.0.0.1:3000").expect("Failed to connect to server");
    println!("Connected to server");

    // Create a buffered reader for the stream
    let reader = BufReader::new(stream);

    // Create a streaming reader
    let mut reader = StreamReader::new(reader);

    // Read the header
    let header = reader.read_header()?;
    println!(
        "Client: Received header with sequence length = {}",
        header.slen
    );

    // Read records as they arrive
    let mut count = 0;
    while let Some(record) = reader.next_record() {
        let record = record?;
        println!(
            "Client: Received record {} with flag = {:?}",
            count,
            record.flag()
        );
        count += 1;
    }

    println!("Client: Received {count} records total");

    Ok(())
}

fn main() -> Result<()> {
    // Create a header for sequences of length 100
    let header = BinseqHeaderBuilder::new().slen(100).build()?;

    // Create some example sequence data
    let sequence = b"ACGT".repeat(25); // 100 nucleotides

    // Spawn the server in a separate thread
    let server_thread = thread::spawn(move || {
        if let Err(e) = server(header, &sequence) {
            eprintln!("Server error: {e:?}");
        }
    });

    // Run the client in the main thread
    if let Err(e) = client() {
        eprintln!("Client error: {e:?}");
    }

    // Wait for the server to finish
    server_thread.join().unwrap();

    Ok(())
}
