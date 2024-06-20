use color_eyre::eyre::Result;
use derive_new::new;

use crate::{primers::PrimerPair, Trim};

// #![warn(missing_docs)]

#[derive(Debug, new, Clone)]
pub struct Record<'a> {
    pub id: &'a [u8],
    pub seq: &'a [u8],
    pub qual: Option<&'a [u8]>,
}

impl<'a> Trim<'a> for Record<'a> {
    fn trim_fastq_to_primers(self, primers: &'a PrimerPair) -> Result<Option<Record<'a>>> {
        let seq_str = std::str::from_utf8(self.seq)?;
        let qual_str = std::str::from_utf8(self.qual.unwrap())?;

        match (&seq_str.find(primers.fwd), &seq_str.find(primers.rev)) {
            (Some(fwd_idx), Some(rev_idx)) => {
                let new_start = fwd_idx + primers.fwd.len();
                let new_end = rev_idx - primers.rev.len();
                let new_seq = seq_str[new_start..new_end].as_bytes();
                let new_qual = qual_str[new_start..new_end].as_bytes();

                let trimmed_record = Record::new(self.id, new_seq, Some(new_qual));

                Ok(Some(trimmed_record))
            }
            _ => Ok(None),
        }
    }

    fn trim_fasta_to_primers(self, primers: &'a PrimerPair) -> Result<Option<Record<'a>>> {
        let seq_str = std::str::from_utf8(self.seq)?;

        match (&seq_str.find(primers.fwd), &seq_str.find(primers.rev)) {
            (Some(fwd_idx), Some(rev_idx)) => {
                let new_start = fwd_idx + primers.fwd.len();
                let new_end = rev_idx - primers.rev.len();
                let new_seq = seq_str[new_start..new_end].as_bytes();

                let trimmed_record = Record::new(self.id, new_seq, None);

                Ok(Some(trimmed_record))
            }
            _ => Ok(None),
        }
    }
}
