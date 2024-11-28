use rand::Rng;
use std::io::Write;
use xz2::write::{XzDecoder, XzEncoder};
use zstd;

fn generate_biased_bits(length: usize, prob_one: f64) -> Vec<u8> {
    let mut rng = rand::thread_rng();
    let mut bits = Vec::with_capacity(length / 8 + 1);
    let mut current_byte = 0u8;
    let mut bit_count = 0;

    for _ in 0..length {
        if rng.gen_bool(prob_one) {
            current_byte |= 1 << (7 - (bit_count % 8));
        }
        bit_count += 1;

        if bit_count % 8 == 0 {
            bits.push(current_byte);
            current_byte = 0;
        }
    }

    // Handle remaining bits if length is not divisible by 8
    if bit_count % 8 != 0 {
        bits.push(current_byte);
    }

    bits
}

fn count_ones(data: &[u8], length: usize) -> usize {
    let mut count = 0;
    for (i, &byte) in data.iter().enumerate() {
        let bits_to_check = if (i + 1) * 8 > length { length % 8 } else { 8 };

        for j in 0..bits_to_check {
            if (byte >> (7 - j)) & 1 == 1 {
                count += 1;
            }
        }
    }
    count
}

fn calculate_entropy(p: f64) -> f64 {
    if p == 0.0 || p == 1.0 {
        return 0.0;
    }
    -(p * p.log2() + (1.0 - p) * (1.0 - p).log2())
}

fn compress_zstd(data: &[u8]) -> Vec<u8> {
    zstd::encode_all(data, 3).unwrap()
}

fn compress_lzma(data: &[u8]) -> Vec<u8> {
    let mut encoder = XzEncoder::new(Vec::new(), 6);
    encoder.write_all(data).unwrap();
    encoder.finish().unwrap()
}

#[derive(Debug)]
struct CompressionStats {
    entropy: f64,
    theoretical_bits: f64,
    zstd_bits: usize,
    lzma_bits: usize,
    information_equivalent_length: f64,
}

fn analyze_compression(length: usize, prob_one: f64, iterations: usize) -> CompressionStats {
    let mut total_zstd_bits = 0;
    let mut total_lzma_bits = 0;

    for _ in 0..iterations {
        let data = generate_biased_bits(length, prob_one);

        let compressed_zstd = compress_zstd(&data);
        let compressed_lzma = compress_lzma(&data);

        total_zstd_bits += compressed_zstd.len() * 8;
        total_lzma_bits += compressed_lzma.len() * 8;
    }

    let entropy = calculate_entropy(prob_one);
    let theoretical_bits = entropy * length as f64;

    let information_equivalent_lenght = 1.44 * length as f64 * -(1. - prob_one).ln();

    CompressionStats {
        entropy,
        theoretical_bits,
        zstd_bits: total_zstd_bits / iterations,
        lzma_bits: total_lzma_bits / iterations,
        information_equivalent_length: information_equivalent_lenght,
    }
}

fn main() {
    let lengths = [1000, 1_0000, 100_000, 1_000_000];
    let probabilities = [0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.5];
    let iterations = 10;

    // output storage overhead vs p = 1/2 representation

    println!("Compression Analysis with Entropy");
    println!("=====================================================================================================");
    println!("Length\tP(=1)\tEntropy\tTheor.\tZstd\tLzma\tZstd/Theor.\tOverhead\tFree Space");
    println!("-----------------------------------------------------------------------------------------------------");

    for &length in &lengths {
        for &prob in &probabilities {
            let stats = analyze_compression(length, prob, iterations);

            println!(
                "{}\t{:.2}\t{:.3}\t{:.0}b\t{}b\t{}b\t{:.2}\t\t{:.2}x\t\t{:.2}x",
                length,
                prob,
                stats.entropy,
                stats.theoretical_bits,
                stats.zstd_bits,
                stats.lzma_bits,
                // stats.lzma_bits,
                stats.zstd_bits as f64 / stats.theoretical_bits,
                // length as f64 / stats.theoretical_bits,
                stats.zstd_bits.min(stats.lzma_bits) as f64 / stats.information_equivalent_length,
                length as f64 / stats.information_equivalent_length,
                // stats.lzma_bits as f64 / stats.theoretical_bits,
            );
        }
        println!("-----------------------------------------------------------------------------------------------------");
    }
}
