use std::collections::HashMap;

use noodles::fastq::Record as FastqRecord;

use futures::io;

use crate::prelude::{Parseable, RecordStream};

pub struct FilterSettings<'a, 'b> {
    pub min_freq: &'a f64,
    pub max_len: &'a usize,
    pub unique_seqs: &'b HashMap<Vec<u8>, f64>,
}

impl<'a, 'b> FilterSettings<'a, 'b> {
    pub fn new(
        min_freq: &'a Option<f64>,
        max_len: &'a Option<usize>,
        unique_seqs: &'b Option<HashMap<Vec<u8>, f64>>,
    ) -> Option<FilterSettings<'a, 'b>> {
        match (min_freq, max_len, unique_seqs) {
            (Some(min_freq), Some(max_len), Some(unique_seqs)) => Some(FilterSettings {
                min_freq,
                max_len,
                unique_seqs,
            }),
            (Some(min_freq), None, Some(unique_seqs)) => Some(FilterSettings {
                min_freq,
                max_len: &123456789,
                unique_seqs,
            }),
            (None, Some(max_len), Some(unique_seqs)) => Some(FilterSettings {
                min_freq: &0.0,
                max_len,
                unique_seqs,
            }),
            (_, _, _) => None,
        }
    }
}

/// Trait `Filtering` handles special filtering requests that require context from the
/// entire sequencing read dataset, e.g., filtering by frequency of each unique sequence.
/// Unlike trait `Trimming`, `Filtering` is not intrinsically amplicon-aware. Instead,
/// it pulls the information it needs from a required index file that can be computed
/// ahead of time with `amplicon-tk index`.
pub trait Filtering: Sized {
    type RecordType: Parseable;

    fn run_filters(self) -> impl std::future::Future<Output = io::Result<Self>>;
}

impl<'a> Filtering for RecordStream<'a, FastqRecord> {
    type RecordType = FastqRecord;
    async fn run_filters(self) -> io::Result<Self> {
        Ok(self)
    }
}
