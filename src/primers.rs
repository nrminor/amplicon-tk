use std::io::BufReader;
use std::{collections::HashMap, fs::File};

use color_eyre::eyre::Result;

use derive_new::new;
use noodles::bed::Reader as BedReader;
use noodles::fasta::io::Reader as FastaReader;

// #![warn(missing_docs)]

struct PrimerSeq<'a> {
    primer_name: String,
    primer_seq: &'a str,
}

#[derive(Debug, new)]
pub struct PrimerPair<'a> {
    pub amplicon: String,
    pub fwd: &'a str,
    pub rev: &'a str,
}

#[derive(Debug)]
pub struct AmpliconScheme<'a> {
    pub scheme: Vec<PrimerPair<'a>>,
}

/// .
///
/// # Errors
///
/// This function will return an error if .
pub fn ref_to_dict(
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

/// .
///
/// # Panics
///
/// Panics if .
///
/// # Errors
///
/// This function will return an error if .
pub fn define_amplicons<'a>(
    mut bed: BedReader<BufReader<File>>,
    ref_dict: &'a HashMap<Vec<u8>, Vec<u8>>,
    fwd_suffix: &'a str,
    rev_suffix: &'a str,
) -> Result<AmpliconScheme<'a>> {
    let all_primer_seqs: Vec<PrimerSeq> = bed
        .records()
        .filter_map(|record| record.ok())
        .map(|record: noodles::bed::Record<4>| -> Result<PrimerSeq> {
            // define the primer name and amplicon name
            let primer_name = record.name().unwrap().to_string();

            // define the ref name and start and stop positions
            let ref_name = record.reference_sequence_name().as_bytes().to_owned();
            let start_pos = record.start_position().get() - 1;
            let stop_pos = record.end_position().get() - 1;

            // pull in the sequence from the ref hashmap
            let primer_seq_bytes = &ref_dict[&ref_name][start_pos..stop_pos];
            let primer_seq = std::str::from_utf8(primer_seq_bytes)?;

            Ok(PrimerSeq {
                primer_name,
                primer_seq,
            })
        })
        .filter_map(|primer_seq| primer_seq.ok())
        .collect();

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

            let fwd = primers
                .iter()
                .filter(|primer| primer.primer_name.contains(fwd_suffix))
                .collect::<Vec<&&PrimerSeq>>()[0];

            let rev = primers
                .iter()
                .filter(|primer| primer.primer_name.contains(rev_suffix))
                .collect::<Vec<&&PrimerSeq>>()[0];

            let pair = PrimerPair {
                amplicon,
                fwd: fwd.primer_seq,
                rev: rev.primer_seq,
            };

            Some(pair)
        })
        .collect::<Vec<PrimerPair>>();

    Ok(AmpliconScheme { scheme })
}
