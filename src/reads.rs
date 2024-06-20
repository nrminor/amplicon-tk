use derive_new::new;
use serde::{Deserialize, Serialize};

use crate::Record;

// #![warn(missing_docs)]

#[derive(Debug, new)]
pub struct Reads<'a> {
    pub reads: Vec<Record<'a>>,
    pub unique_seqs: Vec<SeqFreq>,
    #[new(value = "0.0")]
    pub min_freq: f32,
}

#[derive(Debug, new, Serialize, Deserialize)]
pub struct SeqFreq {
    _seq: String,
    _freq: f32,
}

// pub trait Freq {
//     fn reduce_to_unique_reads() {}

//     fn count_uniques() {}
// }

pub trait Filter {
    // filter_by_freq() {}

    // filter_by_qual() {}

    // filter_by_len() {}

    // filter_invalid_reads() {}
}

pub trait Sort {}

// impl FreqCalculator for Reads {}
