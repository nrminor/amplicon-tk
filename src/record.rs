// #![warn(missing_docs)]

//!

use color_eyre::eyre::Result;
use noodles::fastq::Record as FastqRecord;

use crate::primers::PrimerPair;

///
pub trait FindAmplicons<'a, 'b> {
    ///
    fn matches_forward(&'a self, pair: &'b PrimerPair) -> bool;

    ///
    fn matches_reverse(&'a self, pair: &'b PrimerPair) -> bool;

    /// .
    fn amplicon(&'a self, primerpairs: &'b [PrimerPair]) -> Option<&'b PrimerPair>;
}

impl<'a, 'b> FindAmplicons<'a, 'b> for FastqRecord {
    fn matches_forward(&'a self, pair: &'b PrimerPair) -> bool {
        self.sequence()
            .windows(pair.fwd.len())
            .any(|window| window.eq(pair.fwd.as_bytes()))
            || self
                .sequence()
                .windows(pair.fwd.len())
                .any(|window| window.eq(pair.fwd_rc.as_bytes()))
    }

    fn matches_reverse(&'a self, pair: &'b PrimerPair) -> bool {
        self.sequence()
            .windows(pair.fwd.len())
            .any(|window| window.eq(pair.rev.as_bytes()))
            || self
                .sequence()
                .windows(pair.fwd.len())
                .any(|window| window.eq(pair.rev_rc.as_bytes()))
    }

    fn amplicon(&'a self, primerpairs: &'b [PrimerPair]) -> Option<&'b PrimerPair> {
        let amplicon_match: Vec<&PrimerPair> = primerpairs
            .iter()
            .filter(|pair| self.matches_forward(pair) && self.matches_reverse(pair))
            .collect();

        match amplicon_match.len() {
            1 => Some(amplicon_match[0]),
            _ => None,
        }
    }
}

pub async fn trim_records<'a, 'b>(
    record: &'a mut FastqRecord,
    primers: &'b PrimerPair,
) -> Result<Option<&'a FastqRecord>> {
    let seq_str = std::str::from_utf8(record.sequence())?;
    match (&seq_str.find(&primers.fwd), &seq_str.find(&primers.rev)) {
        (Some(fwd_idx), Some(rev_idx)) => {
            let new_start = fwd_idx + primers.fwd.len();
            let new_end = rev_idx;

            *record.sequence_mut() = record.sequence()[new_start..*new_end].to_vec();
            *record.quality_scores_mut() = record.quality_scores()[new_start..*new_end].to_vec();

            Ok(Some(record))
        }
        _ => Ok(None),
    }
}
