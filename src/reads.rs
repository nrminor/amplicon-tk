// #![warn(missing_docs)]

//!

use crate::primers::AmpliconScheme;
use color_eyre::eyre::Result;
use noodles::fastq::Record as FastqRecord;

pub trait Trimming {
    fn run_trimming(&mut self, primer_pairs: AmpliconScheme) -> Result<Self>
    where
        Self: std::marker::Sized;
}

// impl Trimming for Reads {
//     fn run_trimming(&mut self, amplicon_scheme: AmpliconScheme) -> Result<Self> {
//         let primer_pairs = amplicon_scheme.scheme;
//         let record_with_primers: Vec<(&FastqRecord, &PrimerPair)> = self
//             .reads
//             .iter()
//             .filter_map(|record| {
//                 record
//                     .amplicon(&primer_pairs)
//                     .map(|primers_found| (record, primers_found))
//             })
//             .collect();

//         let trimmed_reads: Vec<FastqRecord> = record_with_primers
//             .iter()
//             .filter_map(|(record, primers)| record.trim_to_primers(primers).ok())
//             .flatten()
//             .collect();

//         let new_dataset = Reads {
//             reads: trimmed_reads,
//             unique_seqs: self.unique_seqs.clone(),
//             min_freq: self.min_freq,
//             max_len: self.max_len,
//         };

//         Ok(new_dataset)
//     }
// }

pub trait Filtering {
    fn run_filters(&mut self) -> Result<Self>
    where
        Self: std::marker::Sized;
}

pub trait Sorting {
    fn sort_reads(&mut self) -> Result<Vec<Vec<FastqRecord>>>;
}
