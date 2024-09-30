// #![warn(missing_docs)]

//!

use std::io::BufReader;
use std::{collections::HashMap, fs::File};

use color_eyre::eyre::{eyre, Result};
use derive_new::new;
use noodles::bed::Reader as BedReader;
use noodles::fasta::io::Reader as FastaReader;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

struct PrimerSeq<'a> {
    primer_name: String,
    primer_seq: &'a str,
}

///
#[derive(Debug, new, Hash, Serialize, Deserialize, Eq, PartialEq)]
pub struct PossiblePrimers {
    /// The name or label of the amplicon
    pub amplicon: String,

    /// The forward primer sequence in 5' to 3' orientation
    pub fwd: String,

    /// The reverse complement of the forward primer sequence
    pub fwd_rc: String,

    /// The reverse primer sequence in 5' to 3' orientation
    pub rev: String,

    /// The reverse complement of the reverse primer sequence
    pub rev_rc: String,
}

///
#[derive(Debug, Hash, Serialize, Deserialize, Eq, PartialEq, Clone)]
pub struct AmpliconBounds {
    pub start: usize,
    pub stop: usize,
}

///
#[derive(Debug, Hash, Serialize, Deserialize, PartialEq)]
pub struct AmpliconScheme {
    ///
    pub scheme: Vec<PossiblePrimers>,
}

impl AmpliconScheme {
    pub fn hash_amplicon_scheme(&self) -> Result<String> {
        let encoded_scheme: Vec<u8> = bincode::serialize(self)?;
        let mut hasher = Sha256::new();
        hasher.update(&encoded_scheme);
        let hash = format!("{:?}", hasher.finalize());

        Ok(hash)
    }
}

/// .
///
/// # Errors
///
/// This function will return an error if .
pub async fn ref_to_dict(
    ref_file: &mut FastaReader<BufReader<File>>,
) -> Result<HashMap<Vec<u8>, Vec<u8>>> {
    let ref_dict = ref_file
        .records()
        .filter_map(|record| record.ok())
        .map(|record| {
            let name = record.name().to_owned();
            let sequence = record.sequence().as_ref().to_owned();
            (name, sequence)
        })
        .collect();
    Ok(ref_dict)
}

///
fn get_reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .flat_map(|base| match base {
            'A' => Some('T'),
            'T' => Some('A'),
            'G' => Some('C'),
            'C' => Some('G'),
            'U' => Some('A'),
            _ => None,
        })
        .rev()
        .collect::<String>()
}

///
async fn collect_primer_seqs(
    mut bed: BedReader<BufReader<File>>,
    ref_dict: &HashMap<Vec<u8>, Vec<u8>>,
) -> Result<Vec<PrimerSeq>> {
    let all_primer_seqs: Vec<PrimerSeq> = bed
        .records()
        .filter_map(|record| record.ok())
        .map(|record: noodles::bed::Record<4>| -> Result<PrimerSeq> {
            // define the primer name and amplicon name
            let primer_name = record.name().unwrap().to_string();

            // define the ref name and start and stop positions
            let ref_name = record.reference_sequence_name().as_bytes().to_owned();
            let start_pos = record.start_position().get();
            let stop_pos = record.end_position().get();

            // pull in the sequence from the ref hashmap
            let seq = ref_dict.get(&ref_name).unwrap();
            match stop_pos <= seq.len() {
                true => {
                    let primer_seq_bytes = &seq[start_pos..stop_pos];
                    let primer_seq = std::str::from_utf8(primer_seq_bytes)?;

                    Ok(PrimerSeq {
                        primer_name,
                        primer_seq,
                    })
                }
                false => {
                    eprintln!("{}", seq.len());
                    let message = format!(
                        "Positions {} and {} for {} are not present in the reference sequence, {}. The reference sequence is:\n\n{}\n",
                        &start_pos,
                        &stop_pos,
                        &primer_name,
                        String::from_utf8(ref_name)?,
                        String::from_utf8(seq.clone())?
                    );
                    eprintln!("{}", &message);
                    Err(eyre!(message))
                }
            }
        })
        .filter_map(|primer_seq| primer_seq.ok())
        .collect();
    Ok(all_primer_seqs)
}

pub trait DefineAmplicons {
    fn define_amplicons<'a>(
        self,
        ref_dict: &'a HashMap<Vec<u8>, Vec<u8>>,
        fwd_suffix: &'a str,
        rev_suffix: &'a str,
    ) -> impl futures::Future<Output = Result<AmpliconScheme>>;
}

impl DefineAmplicons for BedReader<BufReader<File>> {
    async fn define_amplicons<'a>(
        self,
        ref_dict: &'a HashMap<Vec<u8>, Vec<u8>>,
        fwd_suffix: &'a str,
        rev_suffix: &'a str,
    ) -> Result<AmpliconScheme> {
        let all_primer_seqs = collect_primer_seqs(self, ref_dict).await?;

        let amplicons = all_primer_seqs
            .iter()
            .map(|primer_seq| {
                primer_seq
                    .primer_name
                    .replace(fwd_suffix, "")
                    .replace(rev_suffix, "")
            })
            .collect::<Vec<String>>();

        let scheme = amplicons
            .into_iter()
            .filter_map(|amplicon| {
                let primers = all_primer_seqs
                    .iter()
                    .filter(|primer| primer.primer_name.contains(&amplicon))
                    .collect::<Vec<&PrimerSeq>>();

                if primers.len() != 2 {
                    return None;
                }

                let fwd_hits = primers
                    .iter()
                    .filter(|primer| primer.primer_name.contains(fwd_suffix))
                    .collect::<Vec<&&PrimerSeq>>();
                let fwd = fwd_hits.first();

                let rev_hits = primers
                    .iter()
                    .filter(|primer| primer.primer_name.contains(rev_suffix))
                    .collect::<Vec<&&PrimerSeq>>();
                let rev = rev_hits.first();

                if let (Some(fwd), Some(rev)) = (fwd, rev) {
                    let fwd_rc = get_reverse_complement(fwd.primer_seq);
                    let rev_rc = get_reverse_complement(rev.primer_seq);
                    let pair = PossiblePrimers {
                        amplicon,
                        fwd: fwd.primer_seq.to_owned(),
                        fwd_rc,
                        rev: rev.primer_seq.to_owned(),
                        rev_rc,
                    };
                    Some(pair)
                } else {
                    None
                }
            })
            .collect::<Vec<PossiblePrimers>>();

        Ok(AmpliconScheme { scheme })
    }
}

/// .
///
/// # Panics
///
/// Panics if .
///
/// # Errors
///
/// This function will return an error if .
pub async fn define_amplicons<'a>(
    bed: BedReader<BufReader<File>>,
    ref_dict: &'a HashMap<Vec<u8>, Vec<u8>>,
    fwd_suffix: &'a str,
    rev_suffix: &'a str,
) -> Result<AmpliconScheme> {
    let all_primer_seqs = collect_primer_seqs(bed, ref_dict).await?;

    let amplicons = all_primer_seqs
        .iter()
        .map(|primer_seq| {
            primer_seq
                .primer_name
                .replace(fwd_suffix, "")
                .replace(rev_suffix, "")
        })
        .collect::<Vec<String>>();

    let scheme = amplicons
        .into_iter()
        .filter_map(|amplicon| {
            let primers = all_primer_seqs
                .iter()
                .filter(|primer| primer.primer_name.contains(&amplicon))
                .collect::<Vec<&PrimerSeq>>();

            if primers.len() != 2 {
                return None;
            }

            let fwd_hits = primers
                .iter()
                .filter(|primer| primer.primer_name.contains(fwd_suffix))
                .collect::<Vec<&&PrimerSeq>>();
            let fwd = fwd_hits.first();

            let rev_hits = primers
                .iter()
                .filter(|primer| primer.primer_name.contains(rev_suffix))
                .collect::<Vec<&&PrimerSeq>>();
            let rev = rev_hits.first();

            if let (Some(fwd), Some(rev)) = (fwd, rev) {
                let fwd_rc = get_reverse_complement(fwd.primer_seq);
                let rev_rc = get_reverse_complement(rev.primer_seq);
                let pair = PossiblePrimers {
                    amplicon,
                    fwd: fwd.primer_seq.to_owned(),
                    fwd_rc,
                    rev: rev.primer_seq.to_owned(),
                    rev_rc,
                };
                Some(pair)
            } else {
                None
            }
        })
        .collect::<Vec<PossiblePrimers>>();

    Ok(AmpliconScheme { scheme })
}
