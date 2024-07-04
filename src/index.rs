// pub async fn new(
//     mut reader: noodles::fastq::Reader<BufReader<GzDecoder<File>>>,
//     min_freq: Option<f64>,
//     max_len: Option<usize>,
// ) -> Self {
//     let reads: Vec<FastqRecord> = reader.records().filter_map(|record| record.ok()).collect();

//     let (seq_counts, total_count) =
//         reads
//             .iter()
//             .fold((HashMap::new(), 0), |(mut counts, read_count), read| {
//                 *counts.entry(read.sequence().to_owned()).or_insert(0) += 1;
//                 (counts, read_count + 1)
//             });

//     let unique_seqs: HashMap<Vec<u8>, f64> = seq_counts
//         .into_iter()
//         .map(|(seq, count)| (seq, (count as f64) / (total_count as f64)))
//         .collect();

//     let min_freq = min_freq.unwrap_or(0.0);
// }
