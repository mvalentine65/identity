use std::collections::HashMap;
use std::env;
use std::fs;
use std::fs::File;
use std::io::Write;

fn read_fasta_file(filename: &str) -> HashMap<String, String> {
    let mut sequences = HashMap::new();
    let mut sequence_id = None;
    let mut sequence = String::new();
    let contents = fs::read_to_string(filename).expect("Failed to read file");

    for line in contents.lines() {
        if line.starts_with('>') {
            if let Some(id) = sequence_id.take() {
                sequences.insert(id, sequence);
                sequence = String::new();
            }
            sequence_id = Some(line[1..].to_string());
        } else {
            sequence.push_str(line.trim());
        }
    }

    if let Some(id) = sequence_id {
        sequences.insert(id, sequence);
    }

    sequences
}

fn pairwise_identity(seq1: &[u8], seq2: &[u8]) -> f64 {
    if seq1.len() != seq2.len() {
        panic!("Sequences must be of equal length");
    }

    let mut identical = 0;
    let mut aligned_length = 0;

    for (a, b) in seq1.zip(seq2) {
        if a == '-' || b == '-' {
            continue;
        }
        aligned_length += 1;
        if a == b {
            identical += 1;
        }
    }

    if aligned_length == 0 {
        return 0.0;
    }

    (identical as f64 / aligned_length as f64) * 100.0
}

fn average_pairwise_identity(sequences: &HashMap<String, String>) -> f64 {
    // let sequence_ids = sequences.keys().cloned().collect::<Vec<_>>();
    let seqs: Vec<&[u8]> = sequences.values().map(|x| x.as_bytes()).collect();
    let mut percents = Vec::<f64>::with_capacity(seqs.len());
    for column in 0..seqs[0].len() {
        let mut num_pairs = 0;
        let mut total_identity = 0;
        for row1 in 0..(seqs.len() - 1) {
            let a: u8 = seqs[row1][column];
            if a == b'-' { continue; }
            for row2 in (row1 + 1)..seqs.len() {
                let b: u8 = seqs[row2][column];
                if b == b'-' { continue; }
                total_identity += 1;
                if a == b {
                    num_pairs += 1
                }
            }
        }
        if total_identity == 0 {
            percents.push(0.0)
        } else {
            percents.push(num_pairs as f64 / total_identity as f64)
        }
    }
    match percents.len() == 0 {
        false => (100.0 * percents.iter().sum::<f64>()) / percents.len() as f64,
        true => 0.0_f64
    }
}

fn overlap_identity_filter(sequences: &HashMap<String, String>, threshold: f64) -> Vec<(&String, &String)>{
    let mut results = vec![];

    for (id1, seq1) in sequences {
        let mut identities = vec![];
        // Find the non-gap interval of seq1
        let (start1, end1) = find_non_gap_positions(seq1.as_bytes());

        for (id2, seq2) in sequences {
            if id1 == id2 {
                continue;
            }

        // Find the non-gap interval of seq2
            let (start2, end2) = find_non_gap_positions(seq2.as_bytes());

            let overlap_start = start1.max(start2);
            let overlap_end = end1.min(end2);
            if overlap_start <= overlap_end {
                let identity = pairwise_identity(&seq1.as_bytes()[overlap_start..=overlap_end], &seq2.as_bytes()[overlap_start..=overlap_end]);
                if identity >= threshold {
                    identities.push(identity);
                }
            }
        }
        match identities.iter().sum::<f64>()/identities.len() as f64 > threshold {
            true => results.push((id1, seq1)),
            false => {},
        }
    }
    results
}



// Find the positions of the first and last non-gap characters in a byte slice
fn find_non_gap_positions(seq: &[u8]) -> (usize, usize) {
    let start = seq.iter().position(|&c| c != b'-').unwrap_or(seq.len());
    let end = seq.iter().rposition(|&c| c != b'-').unwrap_or(0);
    (start, end)
}

// Calculate pairwise identities between a query sequence and a vector of subject sequences
fn pairwise_identities(query: &str, subjects: &[&str]) -> Vec<f64> {
    let query_bytes = query.as_bytes();
    let query_len = query_bytes.len();
    let mut pairwise_identities = Vec::with_capacity(subjects.len());
    for subject in subjects {
        let subject_bytes = subject.as_bytes();
        let mut identity_count = 0;
        for (q, s) in query_bytes.iter().zip(subject_bytes) {
            if *q != b'-' && *s != b'-' && q == s {
                identity_count += 1;
            }
        }
        pairwise_identities.push(identity_count as f64 / query_len as f64 * 100.0);
    }
    pairwise_identities
}


fn write_filtered_fasta_file(filename: &str, sequences: Vec<(&String, &String)>) -> std::io::Result<()> {
    let mut file = File::create(filename)?;
    for (id, seq) in sequences {
        file.write_all(b">")?;
        file.write_all(id.as_bytes())?;
        file.write_all(b"\n")?;
        file.write_all(seq.as_bytes())?;
        file.write_all(b"\n")?;
    }
    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        eprintln!("Usage: {} <fasta_file> <avg_identity_threshold> <overlap_identity_threshold>", args[0]);
        std::process::exit(1);
    }

    let filename = &args[1];

    let avg_identity_threshold = args[2].parse::<f64>().unwrap_or_else(|_| {
        eprintln!("Error: Invalid average identiy threshold argument");
        std::process::exit(1);
    });
    let ov_identity_threshold = args[3].parse::<f64>().unwrap_or_else(|_| {
        eprintln!("Error: Invalid overlap threshold value");
        std::process::exit(1);
    });
    let sequences = read_fasta_file(filename);
    let avg_identity = average_pairwise_identity(&sequences);

    println!("Average pairwise identity: {:.2}%", avg_identity);
    if avg_identity < avg_identity_threshold {
        println!("{} below identity threshold {}",filename, avg_identity_threshold);
        std::process::exit(0);
    }

    println!("Filtering sequences with identity above {}", avg_identity_threshold);
    let filtered_sequences = overlap_identity_filter(&sequences, ov_identity_threshold);
    match write_filtered_fasta_file(filename, filtered_sequences) {
        Ok(_) => {},
        Err(_) => {
            eprintln!("Error writing filtered sequences");
            std::process::exit(1);
        },
    }
}
