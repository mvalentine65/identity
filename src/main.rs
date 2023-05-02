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

fn filter_sequences(sequences: &HashMap<String, String>, max_disagreements: usize) -> Vec<(&String,&String)> {
    let seqs: Vec<&[u8]> = sequences.values().map(|x| x.as_bytes()).collect();
    let consensus = make_consensus(&seqs);
    let mut filtered_seqs = Vec::new();
    for (id, seq) in sequences.iter() {
        let seq_bytes = seq.as_bytes();
        let mut disagreements = 0;
        for (i, &base) in seq_bytes.iter().enumerate() {
            if base != b'-' && base != consensus[i] {
                    disagreements += 1;
            }
        }
        if disagreements <= max_disagreements {
            filtered_seqs.push((id, seq));
        }
    }
    filtered_seqs
}

fn make_consensus(seqs: &[&[u8]]) -> Vec<u8> {
    let mut consensus = Vec::with_capacity(seqs[0].len());
    for column in 0..seqs[0].len() {
        let mut counts = vec![0; 256];
        for seq in seqs.iter() {
            let c = seq[column];
            if c != b'-' {
                counts[c as usize] += 1;
            }
        }
        let (majority, _) = counts.iter().enumerate().max_by_key(|&(_, count)| count).unwrap();
        consensus.push(majority as u8);
    }
    consensus
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
        eprintln!("Usage: {} <fasta_file> <identity_threshold> <max_disagrees>", args[0]);
        std::process::exit(1);
    }

    let filename = &args[1];
    let max_disagrees = args[3].parse::<usize>().unwrap_or_else(|_| {
        eprintln!("Error: Invalid maximum disagree value");
        std::process::exit(1);
    });
    let identity_threshold = args[2].parse::<f64>().unwrap_or_else(|_| {
        eprintln!("Error: Invalid identiy threshold argument");
        std::process::exit(1);
    });

    let sequences = read_fasta_file(filename);
    let avg_identity = average_pairwise_identity(&sequences);

    println!("Average pairwise identity: {:.2}%", avg_identity);
    if avg_identity < identity_threshold {
        println!("{} below identity threshold {}",filename, identity_threshold);
        std::process::exit(0);
    }

    println!("Filtering sequences with identity above {}", identity_threshold);
    let filtered_sequences = filter_sequences(&sequences, max_disagrees);
    match write_filtered_fasta_file(filename, filtered_sequences) {
        Ok(_) => {},
        Err(_) => {
            eprintln!("Error writing filtered sequences");
            std::process::exit(1);
        },
    }
}
