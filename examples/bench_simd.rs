use binseq::{is_simd_supported, nuc};
use std::time::{Duration, Instant};
use rand::Rng;

fn generate_random_sequence(len: usize) -> Vec<u8> {
    let mut sequence = Vec::with_capacity(len);
    let nucleotides = [b'A', b'C', b'G', b'T'];
    let mut rng = rand::thread_rng();
    
    for _ in 0..len {
        let idx = rng.gen_range(0..4);
        sequence.push(nucleotides[idx]);
    }
    
    sequence
}

struct BenchmarkResult {
    sequence_length: usize,
    iterations: usize,
    bitnuc_duration: Duration,
    simd_duration: Duration,
    bitnuc_speed: f64,
    simd_speed: f64,
    speedup: f64,
}

fn benchmark_encode(sequence: &[u8], iterations: usize) -> BenchmarkResult {
    // Verify the sequence only contains valid nucleotides
    for &b in sequence {
        match b {
            b'A' | b'C' | b'G' | b'T' => {},
            other => panic!("Invalid nucleotide in sequence: {}", other),
        }
    }
    let chunks = (sequence.len() + 31) / 32;
    
    let mut bitnuc_output = Vec::with_capacity(chunks);
    bitnuc_output.resize(chunks, 0);
    
    let start = Instant::now();
    for _ in 0..iterations {
        bitnuc_output.clear();
        bitnuc_output.resize(chunks, 0);
        bitnuc::encode(sequence, &mut bitnuc_output).unwrap();
    }
    let bitnuc_duration = start.elapsed();
    
    let mut simd_output = vec![0u64; chunks];
    let start = Instant::now();
    for _ in 0..iterations {
        simd_output.clear();
        simd_output.resize(chunks, 0);
        nuc::encode(sequence, &mut simd_output).unwrap();
    }
    let simd_duration = start.elapsed();
    
    // Calculate processing speeds in millions of nucleotides per second
    let bitnuc_speed = (sequence.len() * iterations) as f64 / bitnuc_duration.as_secs_f64() / 1_000_000.0;
    let simd_speed = (sequence.len() * iterations) as f64 / simd_duration.as_secs_f64() / 1_000_000.0;
    let speedup = bitnuc_duration.as_secs_f64() / simd_duration.as_secs_f64();
    
    BenchmarkResult {
        sequence_length: sequence.len(),
        iterations,
        bitnuc_duration,
        simd_duration,
        bitnuc_speed,
        simd_speed,
        speedup,
    }
}

fn benchmark_decode(sequence: &[u8], iterations: usize) -> BenchmarkResult {
    // Verify the sequence only contains valid nucleotides
    for &b in sequence {
        match b {
            b'A' | b'C' | b'G' | b'T' => {},
            other => panic!("Invalid nucleotide in sequence: {}", other),
        }
    }
    let chunks = (sequence.len() + 31) / 32;
    let mut encoded = vec![0u64; chunks];
    bitnuc::encode(sequence, &mut encoded).unwrap();
    
    let mut decoded = Vec::new();
    let start = Instant::now();
    for _ in 0..iterations {
        decoded.clear();
        bitnuc::decode(&encoded, sequence.len(), &mut decoded).unwrap();
    }
    let bitnuc_duration = start.elapsed();
    
    let mut decoded = Vec::new();
    let start = Instant::now();
    for _ in 0..iterations {
        decoded.clear();
        nuc::decode(&encoded, sequence.len(), &mut decoded).unwrap();
    }
    let simd_duration = start.elapsed();
    
    // Calculate processing speeds in millions of nucleotides per second
    let bitnuc_speed = (sequence.len() * iterations) as f64 / bitnuc_duration.as_secs_f64() / 1_000_000.0;
    let simd_speed = (sequence.len() * iterations) as f64 / simd_duration.as_secs_f64() / 1_000_000.0;
    let speedup = bitnuc_duration.as_secs_f64() / simd_duration.as_secs_f64();
    
    BenchmarkResult {
        sequence_length: sequence.len(),
        iterations,
        bitnuc_duration,
        simd_duration,
        bitnuc_speed,
        simd_speed,
        speedup,
    }
}

fn print_table_header(operation: &str) {
    println!("\n{} Benchmark Results", operation);
    println!("{:-<106}", "");
    println!("| {:<10} | {:<10} | {:<15} | {:<15} | {:<12} | {:<12} | {:<8} |", 
        "Length", "Iterations", "Original (ms)", "SIMD (ms)", "Orig (Mnuc/s)", "SIMD (Mnuc/s)", "Speedup");
    println!("{:-<106}", "");
}

fn print_table_row(result: &BenchmarkResult) {
    println!("| {:<10} | {:<10} | {:<15} | {:<15} | {:<13.2} | {:<13.2} | {:<8.2} |",
        result.sequence_length,
        result.iterations,
        result.bitnuc_duration.as_millis(),
        result.simd_duration.as_millis(),
        result.bitnuc_speed,
        result.simd_speed,
        result.speedup);
}

fn print_table_footer() {
    println!("{:-<106}", "");
}

fn main() {
    // Print architecture information (ARM only)
    println!("Architecture: aarch64 (ARM64)");
    println!("SIMD support detected: {}", is_simd_supported());
    
    // Test with different sequence lengths
    let sequence_lengths = [
        32,      // Exactly one SIMD vector on x86_64 or 2 vectors on aarch64
        64,      // Two SIMD vectors on x86_64 or 4 vectors on aarch64
        100,     // Not a multiple of vector size
        1000,    // Longer sequence
        10000,   // Much longer sequence
        100000,  // Very long sequence
        1000000, // Extremely long sequence
        5000000, // Extremely Extremely long sequence
    ];
    
    let iterations = 1000;
    
    // Collect encode benchmark results
    let mut encode_results = Vec::new();
    for &len in &sequence_lengths {
        let sequence = generate_random_sequence(len);
        encode_results.push(benchmark_encode(&sequence, iterations));
    }
    
    // Collect decode benchmark results
    let mut decode_results = Vec::new();
    for &len in &sequence_lengths {
        let sequence = generate_random_sequence(len);
        decode_results.push(benchmark_decode(&sequence, iterations));
    }
    
    // Print encode results table
    print_table_header("Encode");
    for result in &encode_results {
        print_table_row(result);
    }
    print_table_footer();
    
    // Print decode results table
    print_table_header("Decode");
    for result in &decode_results {
        print_table_row(result);
    }
    print_table_footer();
}