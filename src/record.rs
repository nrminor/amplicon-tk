// #![warn(missing_docs)]

//!

use itertools::Itertools;
use noodles::fastq::Record as FastqRecord;
use pretty_assertions::assert_eq;

use crate::{
    primers::{AmpliconBounds, PossiblePrimers},
    reads::FilterSettings,
};

///
pub trait FindAmplicons<'a, 'b> {
    ///
    fn find_primer_match(&'a self, primer: &'b str, rc_primer: &'b str) -> Option<usize>;

    /// .
    fn find_amplicon(
        &'a self,
        primerpairs: &'b [PossiblePrimers],
    ) -> impl futures::Future<Output = Option<AmpliconBounds>>;

    ///
    fn to_bounds(self, bounds: AmpliconBounds) -> impl futures::Future<Output = Self>
    where
        Self: Sized;

    ///
    fn whether_to_write(
        &'a self,
        filters: &'b Option<FilterSettings>,
    ) -> impl futures::Future<Output = bool>;
}

impl<'a, 'b> FindAmplicons<'a, 'b> for FastqRecord {
    fn find_primer_match(&'a self, primer: &'b str, rc_primer: &'b str) -> Option<usize> {
        assert_eq!(
            primer.len(),
            rc_primer.len(),
            "Primer length mismatch:\n{}\n{}",
            primer,
            rc_primer
        );
        let primer_hit = self
            .sequence()
            .windows(primer.len())
            .position(|window| window.eq(primer.as_bytes()));
        let rc_primer_hit = self
            .sequence()
            .windows(rc_primer.len())
            .position(|window| window.eq(rc_primer.as_bytes()));
        match (primer_hit, rc_primer_hit) {
            (Some(_), Some(_)) => None, // ambiguous case where both a primer and its reverse complement are found, which should be rare
            (Some(hit), None) => Some(hit),
            (None, Some(hit)) => Some(hit),
            (None, None) => None,
        }
    }

    async fn find_amplicon(&'a self, primerpairs: &'b [PossiblePrimers]) -> Option<AmpliconBounds> {
        let mut amplicon_match: Vec<AmpliconBounds> = primerpairs
            .iter()
            .filter_map(|pair| {
                let maybe_fwd = self.find_primer_match(&pair.fwd, &pair.fwd_rc);
                let maybe_rev = self.find_primer_match(&pair.rev, &pair.rev_rc);
                match (maybe_fwd, maybe_rev) {
                    (Some(fwd), Some(rev)) => {
                        let (amplicon_start, amplicon_stop) = match fwd < rev {
                            true => (fwd + pair.fwd.len() - 1, rev),
                            false => (rev + pair.rev.len() - 1, fwd),
                        };
                        let amplicon_len = amplicon_stop - amplicon_start;
                        if amplicon_len > pair.fwd.len()
                            && amplicon_len > pair.rev.len()
                            && amplicon_stop != amplicon_start
                        {
                            Some(AmpliconBounds {
                                start: amplicon_start,
                                stop: amplicon_stop,
                            })
                        } else {
                            None
                        }
                    }
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

    async fn to_bounds(mut self, bounds: AmpliconBounds) -> Self {
        *self.sequence_mut() = self.sequence()[bounds.start..bounds.stop].to_vec();
        *self.quality_scores_mut() = self.quality_scores()[bounds.start..bounds.stop].to_vec();
        assert_eq!(
            self.sequence().len(),
            self.quality_scores().len(),
            "Trimming mistake encountered where the number of bases does not equal the number of quality scores:\n{}\n{}",
            String::from_utf8(self.sequence().to_vec()).unwrap(),
            String::from_utf8(self.quality_scores().to_vec()).unwrap(),
        );
        self
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
