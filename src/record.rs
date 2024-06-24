use color_eyre::eyre::Result;
use noodles::fastq::Record as FastqRecord;

use crate::primers::PrimerPair;

// #![warn(missing_docs)]

/// .
///
/// # Panics
///
/// Panics if .
pub fn find_amplicon<'a>(
    record: &'a FastqRecord,
    primerpairs: &'a [PrimerPair<'a>],
) -> Option<&'a PrimerPair<'a>> {
    let amplicon_match: Vec<&PrimerPair> = primerpairs
        .iter()
        .filter(|pair| {
            std::str::from_utf8(record.sequence())
                .unwrap()
                .contains(pair.fwd)
                && std::str::from_utf8(record.sequence())
                    .unwrap()
                    .contains(pair.rev)
        })
        .collect();

    match amplicon_match.len() {
        1 => Some(amplicon_match[0]),
        _ => None,
    }
}

pub trait Trim<'a> {
    /// .
    ///
    /// # Errors
    ///
    /// This function will return an error if .
    fn trim_to_primers(&self, primers: &'a PrimerPair) -> Result<Option<Self>>
    where
        Self: Sized;
}

impl<'a> Trim<'a> for FastqRecord {
    fn trim_to_primers(&self, primers: &'a PrimerPair) -> Result<Option<Self>> {
        let seq_str = std::str::from_utf8(self.sequence())?;
        match (&seq_str.find(primers.fwd), &seq_str.find(primers.rev)) {
            (Some(fwd_idx), Some(rev_idx)) => {
                let new_start = fwd_idx + primers.fwd.len();
                let new_end = rev_idx - primers.rev.len();

                let trimmed_record = FastqRecord::new(
                    *self.definition(),
                    self.sequence()[new_start..new_end].to_vec(),
                    self.quality_scores()[new_start..new_end].to_vec();
                );

                Ok(Some(trimmed_record))
            }
            _ => Ok(None),
        }
    }
}

// impl<'a> Trim<'a> for FastaRecord {
//     fn trim_to_primers(&mut self, primers: &'a PrimerPair) -> Result<Option<&Self>> {
//         let seq_str = std::str::from_utf8(self.sequence().as_ref())?;
//         match (&seq_str.find(primers.fwd), &seq_str.find(primers.rev)) {
//             (Some(fwd_idx), Some(rev_idx)) => {
//                 let new_start = fwd_idx + primers.fwd.len();
//                 let new_end = rev_idx - primers.rev.len();

//                 let new_seq = Sequence::from(self.sequence().as_ref()[new_start..new_end].to_vec());

//                 let trimmed_record =
//                     noodles::fasta::Record::new(self.definition().clone(), new_seq);

//                 Ok(Some(trimmed_record))
//             }
//             _ => Ok(None),
//         }
//     }
// }
