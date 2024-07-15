use flate2::bufread::GzDecoder;
use noodles::fastq::Reader as FastqReader;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::io::Write;
use std::path::Path;
use std::{collections::HashMap, fs::File, io::BufReader};

use color_eyre::eyre::Result;

use crate::io::FastqGz;
use crate::io::{Fastq, SupportedFormat};
use crate::primers::AmpliconScheme;
use crate::reads::sync_trimming;

#[derive(Debug, Serialize, Deserialize)]
struct IndexFormat {
    hash: String,
    unique_seqs: HashMap<Vec<u8>, f64>,
}

pub trait Index: SupportedFormat {
    type Reader: Unpin + Send;
    fn index(self, reader: Self::Reader, scheme: AmpliconScheme, input_file: &Path) -> Result<()>;
    // fn check_for_index(self);
    // fn load_index(self);
}

impl Index for Fastq {
    type Reader = FastqReader<BufReader<File>>;
    fn index(
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

        // collect the reads into an eager vector
        let reads = reader.records().filter_map(|record| record.ok());

        // trim them down based on the amplicon scheme
        let handle = tokio::runtime::Handle::current();
        let reads = handle.block_on(sync_trimming(reads, &scheme))?;

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
    fn index(
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

        // trim them down based on the amplicon scheme, which will require temporarily blocking
        // the tokio asynchronous runtime
        let handle = tokio::runtime::Handle::current();
        let reads = handle.block_on(sync_trimming(reads, &scheme))?;

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
