// #![warn(missing_docs)]

//!

use color_eyre::eyre::Result;
use itertools::Itertools;
use noodles::fastq::Record as FastqRecord;

use crate::{
    primers::{PossiblePrimers, PrimerPair},
    reads::FilterSettings,
};

///
pub trait FindAmplicons<'a, 'b> {
    ///
    fn forward_match(&'a self, pair: &'b PossiblePrimers) -> Option<&'b str>;

    ///
    fn reverse_match(&'a self, pair: &'b PossiblePrimers) -> Option<&'b str>;

    /// .
    fn find_amplicon(
        &'a self,
        primerpairs: &'b [PossiblePrimers],
    ) -> impl futures::Future<Output = Option<PrimerPair>>;

    ///
    fn trim_to_amplicon(
        self,
        primers: PrimerPair,
    ) -> impl futures::Future<Output = Result<Option<Self>>>
    where
        Self: Sized;

    ///
    fn whether_to_write(
        &'a self,
        filters: &'b Option<FilterSettings>,
    ) -> impl futures::Future<Output = bool>;
}

impl<'a, 'b> FindAmplicons<'a, 'b> for FastqRecord {
    fn forward_match(&'a self, pair: &'b PossiblePrimers) -> Option<&'b str> {
        if self
            .sequence()
            .windows(pair.fwd.len())
            .any(|window| window.eq(pair.fwd.as_bytes()))
        {
            Some(&pair.fwd)
        } else if self
            .sequence()
            .windows(pair.fwd.len())
            .any(|window| window.eq(pair.fwd_rc.as_bytes()))
        {
            Some(&pair.fwd_rc)
        } else {
            None
        }
    }

    fn reverse_match(&'a self, pair: &'b PossiblePrimers) -> Option<&'b str> {
        if self
            .sequence()
            .windows(pair.rev.len())
            .any(|window| window.eq(pair.rev.as_bytes()))
        {
            Some(&pair.rev)
        } else if self
            .sequence()
            .windows(pair.rev.len())
            .any(|window| window.eq(pair.rev_rc.as_bytes()))
        {
            Some(&pair.rev_rc)
        } else {
            None
        }
    }

    async fn find_amplicon(&'a self, primerpairs: &'b [PossiblePrimers]) -> Option<PrimerPair> {
        let mut amplicon_match: Vec<PrimerPair> = primerpairs
            .iter()
            .filter_map(|pair| {
                let maybe_fwd = self.forward_match(pair);
                let maybe_rev = self.reverse_match(pair);

                match (maybe_fwd, maybe_rev) {
                    (Some(fwd), Some(rev)) => Some(PrimerPair {
                        fwd: fwd.to_string(),
                        rev: rev.to_string(),
                    }),
                    _ => None,
                }
            })
            .unique()
            .collect();

        match (amplicon_match.len(), amplicon_match.pop()) {
            (1, Some(success)) => Some(success),
            _ => None,
        }
    }

    async fn trim_to_amplicon(mut self, primers: PrimerPair) -> Result<Option<Self>> {
        let seq_str = std::str::from_utf8(self.sequence())?;
        match (&seq_str.find(&primers.fwd), &seq_str.find(&primers.rev)) {
            (Some(fwd_idx), Some(rev_idx)) => {
                let new_start = fwd_idx + primers.fwd.len();
                let new_end = rev_idx;

                if &new_start >= new_end {
                    return Ok(None);
                }

                *self.sequence_mut() = self.sequence()[new_start..*new_end].to_vec();
                *self.quality_scores_mut() = self.quality_scores()[new_start..*new_end].to_vec();

                Ok(Some(self))
            }
            _ => Ok(None),
        }
    }

    async fn whether_to_write(&'a self, filters: &'b Option<FilterSettings<'_, '_>>) -> bool {
        if let Some(filters) = filters {
            let seq = self.sequence().to_vec();
            let seq_len = seq.len();
            if let Some(freq) = filters.unique_seqs.get(&seq) {
                freq >= filters.min_freq && &seq_len <= filters.max_len
            } else {
                false
            }
        } else {
            true
        }
    }
}
