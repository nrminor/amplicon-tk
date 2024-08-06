use flate2::bufread::GzDecoder;
use noodles::fastq::Reader as FastqReader;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::io::{Read, Write};
use std::path::Path;
use std::{collections::HashMap, fs::File, io::BufReader};

use color_eyre::eyre::Result;

use crate::amplicons::AmpliconScheme;
use crate::io::FastqGz;
use crate::io::{Fastq, SupportedFormat};
use crate::reads::sync_trimming;

#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct IndexFormat {
    hash: String,
    pub unique_seqs: HashMap<Vec<u8>, f64>,
}

pub trait Index: SupportedFormat {
    type Reader: Unpin + Send;
    fn index(
        self,
        reader: Self::Reader,
        scheme: AmpliconScheme,
        input_file: &Path,
    ) -> impl futures::Future<Output = Result<()>>;
    fn load_index(
        &self,
        input_file: &Path,
        current_hash: &str,
    ) -> Result<Option<HashMap<Vec<u8>, f64>>> {
        let index_filename = format!("{}.ampidx", input_file.to_string_lossy());
        let index_file = File::open(&index_filename);
        let potential_index = match index_file {
            Err(_) => None,
            Ok(mut file) => {
                let mut buffer = Vec::new();
                file.read_to_end(&mut buffer)?;
                let index: IndexFormat = serde_cbor::from_slice(&buffer)?;
                match index.hash.eq(current_hash) {
                    true => Some(index),
                    false => {
                        eprintln!(
                            "An index for the current sample, {}, was found, but it was built with a different primer scheme. As such, filtering cannot be performed. Please rerun indexing before attempting to filter.",
                            &index_filename
                        );
                        None
                    }
                }
            }
        };
        let unique_seqs = if let Some(index) = potential_index {
            Some(index.unique_seqs)
        } else {
            None
        };

        Ok(unique_seqs)
    }
}

impl Index for Fastq {
    type Reader = FastqReader<BufReader<File>>;
    async fn index(
        self,
        mut reader: Self::Reader,
        scheme: AmpliconScheme,
        input_file: &Path,
    ) -> Result<()> {
        // hash the amplicon scheme
        let hash = scheme.hash_amplicon_scheme()?;

        // collect the reads into an eager vector
        let reads = reader.records().filter_map(|record| record.ok());

        // trim them down based on the amplicon scheme
        let reads = sync_trimming(reads, &scheme).await?;

        // use the trimmed sequences to find and count unique amplicon sequences
        let (seq_counts, total_count) =
            reads
                .iter()
                .fold((HashMap::new(), 0), |(mut counts, read_count), read| {
                    *counts.entry(read.sequence().to_owned()).or_insert(0) += 1;
                    (counts, read_count + 1)
                });

        // compute the prevalence for each sequence
        let unique_seqs: HashMap<Vec<u8>, f64> = seq_counts
            .into_iter()
            .map(|(seq, count)| (seq, (count as f64) / (total_count as f64)))
            .collect();
        let format = IndexFormat { hash, unique_seqs };

        let serialized_index = serde_cbor::to_vec(&format)?;

        let index_filename = format!("{}.ampidx", input_file.to_string_lossy());
        let mut file = File::create(index_filename)?;
        file.write_all(&serialized_index)?;

        Ok(())
    }
}

impl Index for FastqGz {
    type Reader = FastqReader<BufReader<GzDecoder<BufReader<File>>>>;
    async fn index(
        self,
        mut reader: Self::Reader,
        scheme: AmpliconScheme,
        input_file: &Path,
    ) -> Result<()> {
        // hash the amplicon scheme
        let encoded_scheme: Vec<u8> = bincode::serialize(&scheme)?;
        let mut hasher = Sha256::new();
        hasher.update(&encoded_scheme);
        let hash = format!("{:?}", hasher.finalize());

        // collect the reads into a lazy iterator
        let reads = reader.records().filter_map(|record| record.ok());

        // trim them down based on the amplicon scheme
        let reads = sync_trimming(reads, &scheme).await?;

        // use the trimmed sequences to find and count unique amplicon sequences
        let (seq_counts, total_count) =
            reads
                .iter()
                .fold((HashMap::new(), 0), |(mut counts, read_count), read| {
                    *counts.entry(read.sequence().to_owned()).or_insert(0) += 1;
                    (counts, read_count + 1)
                });

        // compute the prevalence for each sequence
        let unique_seqs: HashMap<Vec<u8>, f64> = seq_counts
            .into_iter()
            .map(|(seq, count)| (seq, (count as f64) / (total_count as f64)))
            .collect();
        let format = IndexFormat { hash, unique_seqs };

        let serialized = serde_cbor::to_vec(&format)?;

        let index_filename = format!("{}.ampidx", input_file.to_string_lossy());
        let mut file = File::create(index_filename)?;
        file.write_all(&serialized)?;

        Ok(())
    }
}
