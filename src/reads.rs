use crate::{
    primers::{AmpliconScheme, PrimerPair},
    record::{FindAmplicons, Trim},
};
use color_eyre::eyre::Result;
use noodles::fastq::Record as FastqRecord;
use std::{collections::HashMap, fs::File, io::BufReader};

// #![warn(missing_docs)]

#[derive(Debug)]
pub struct Reads {
    pub reads: Vec<FastqRecord>,
    pub unique_seqs: HashMap<Vec<u8>, f64>,
    pub min_freq: f64,
    pub max_len: Option<usize>,
}

impl Reads {
    /// Creates a new [`Reads`].
    pub fn new(
        mut reader: noodles::fastq::Reader<BufReader<File>>,
        min_freq: Option<f64>,
        max_len: Option<usize>,
    ) -> Self {
        let reads: Vec<FastqRecord> = reader.records().filter_map(|record| record.ok()).collect();

        let (seq_counts, total_count) =
            reads
                .iter()
                .fold((HashMap::new(), 0), |(mut counts, read_count), read| {
                    *counts.entry(read.sequence().to_owned()).or_insert(0) += 1;
                    (counts, read_count + 1)
                });

        let unique_seqs: HashMap<Vec<u8>, f64> = seq_counts
            .into_iter()
            .map(|(seq, count)| (seq, (count as f64) / (total_count as f64)))
            .collect();

        let min_freq = min_freq.unwrap_or(0.0);

        Self {
            reads,
            unique_seqs,
            min_freq,
            max_len,
        }
    }
}

pub trait Trimming {
    /// .
    ///
    /// # Errors
    ///
    /// This function will return an error if .
    fn run_trimming(&mut self, primer_pairs: AmpliconScheme) -> Result<Self>
    where
        Self: std::marker::Sized;
}

impl Trimming for Reads {
    fn run_trimming(&mut self, amplicon_scheme: AmpliconScheme) -> Result<Self> {
        let primer_pairs = amplicon_scheme.scheme;
        let record_with_primers: Vec<(&FastqRecord, &PrimerPair)> = self
            .reads
            .iter()
            .filter_map(|record| {
                record
                    .amplicon(&primer_pairs)
                    .map(|primers_found| (record, primers_found))
            })
            .collect();

        let trimmed_reads: Vec<FastqRecord> = record_with_primers
            .iter()
            .filter_map(|(record, primers)| record.trim_to_primers(primers).ok())
            .flatten()
            .collect();

        let new_dataset = Reads {
            reads: trimmed_reads,
            unique_seqs: self.unique_seqs.clone(),
            min_freq: self.min_freq,
            max_len: self.max_len,
        };

        Ok(new_dataset)
    }
}

pub trait Filtering {
    /// .
    ///
    /// # Errors
    ///
    /// This function will return an error if .
    fn run_filters(&mut self) -> Result<Self>
    where
        Self: std::marker::Sized;
}

impl Filtering for Reads {
    fn run_filters(&mut self) -> Result<Self> {
        let filtered_reads = self
            .reads
            .iter()
            .filter(|read| {
                let seq = &read.sequence().to_owned();
                let freq = self.unique_seqs.get(seq);
                matches!(freq, Some(freq) if *freq >= self.min_freq)
            })
            .filter(|read| match self.max_len {
                Some(len) => len < read.sequence().len(),
                None => false,
            })
            .map(|record| record.to_owned())
            .collect();

        Ok(Reads {
            reads: filtered_reads,
            unique_seqs: self.unique_seqs.clone(),
            min_freq: self.min_freq,
            max_len: self.max_len,
        })
    }
}

pub trait Sorting {
    /// .
    ///
    /// # Errors
    ///
    /// This function will return an error if .
    fn sort_reads(&mut self) -> Result<Vec<Vec<FastqRecord>>>;
}

impl Sorting for Reads {
    fn sort_reads(&mut self) -> Result<Vec<Vec<FastqRecord>>> {
        let mut sorted_reads: HashMap<Vec<u8>, Vec<FastqRecord>> = HashMap::new();

        for key in self.unique_seqs.keys() {
            sorted_reads.insert(key.clone(), Vec::new());
        }

        assert_eq!(&self.unique_seqs.keys().len(), &sorted_reads.keys().len());

        self.reads.iter().for_each(|read| {
            let seq = read.sequence().to_owned();
            sorted_reads.entry(seq).or_default().push(read.clone());
        });

        let read_stacks: Vec<Vec<FastqRecord>> = sorted_reads.into_values().collect();

        Ok(read_stacks)
    }
}
