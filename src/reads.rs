// #![warn(missing_docs)]

//!

use futures::TryStreamExt;
use futures::{future::join_all, Future};
use noodles::fastq::Record as FastqRecord;
use std::{collections::HashMap, path::Path};

use crate::{
    io::{Fastq, FastqGz, Init, SeqWriter, SupportedFormat},
    primers::AmpliconScheme,
    record::FindAmplicons,
};
use color_eyre::eyre::Result;

pub struct FilterSettings<'a, 'b> {
    pub min_freq: &'a f64,
    pub max_len: &'a usize,
    pub unique_seqs: &'b HashMap<Vec<u8>, f64>,
}

impl<'a, 'b> FilterSettings<'a, 'b> {
    pub fn new(
        min_freq: &'a Option<f64>,
        max_len: &'a Option<usize>,
        unique_seqs: &'b Option<HashMap<Vec<u8>, f64>>,
    ) -> Option<FilterSettings<'a, 'b>> {
        match (min_freq, max_len, unique_seqs) {
            (Some(min_freq), Some(max_len), Some(unique_seqs)) => Some(FilterSettings {
                min_freq,
                max_len,
                unique_seqs,
            }),
            (Some(min_freq), None, Some(unique_seqs)) => Some(FilterSettings {
                min_freq,
                max_len: &123456789,
                unique_seqs,
            }),
            (None, Some(max_len), Some(unique_seqs)) => Some(FilterSettings {
                min_freq: &0.0,
                max_len,
                unique_seqs,
            }),
            (_, _, _) => None,
        }
    }
}

pub trait Trimming: SupportedFormat {
    type Record;
    fn trim(
        self,
        input_path: &Path,
        output_path: &Path,
        scheme: AmpliconScheme,
        _filters: Option<FilterSettings>,
    ) -> impl Future<Output = Result<()>>;
}

impl Trimming for Fastq {
    type Record = FastqRecord;
    async fn trim(
        self,
        input_path: &Path,
        output_path: &Path,
        scheme: AmpliconScheme,
        filters: Option<FilterSettings<'_, '_>>,
    ) -> Result<()> {
        let (mut reader, format) = self.init(input_path).await?;
        let mut records = reader.records();
        let mut writer = format.read_writer(output_path).await?;

        // iterate through records asynchronously, find amplicon hits, and trim them down to
        // exclude primers and anything that extends beyond them
        while let Some(record) = records.try_next().await? {
            let amplicon_hit = record.find_amplicon(&scheme.scheme).await;
            if let Some(hit) = amplicon_hit {
                let trimmed = record.trim_to_amplicon(hit).await?;
                match trimmed {
                    Some(trimmed_record) => match trimmed_record.whether_to_write(&filters).await {
                        true => writer.write_record(&trimmed_record).await?,
                        false => continue,
                    },
                    _ => continue,
                }
            } else {
                continue;
            }
        }

        // Finalize the written contents to make sure the file is not corrupted
        format.finalize_write(writer).await?;

        Ok(())
    }
}

impl Trimming for FastqGz {
    type Record = FastqRecord;
    async fn trim(
        self,
        input_path: &Path,
        output_path: &Path,
        scheme: AmpliconScheme,
        filters: Option<FilterSettings<'_, '_>>,
    ) -> Result<()> {
        let (mut reader, format) = self.init(input_path).await?;
        let mut records = reader.records();
        let mut writer = format.read_writer(output_path).await?;

        // iterate through records asynchronously, find amplicon hits, and trim them down to
        // exclude primers and anything that extends beyond them
        while let Some(record) = records.try_next().await? {
            let amplicon_hit = record.find_amplicon(&scheme.scheme).await;
            if let Some(hit) = amplicon_hit {
                let trimmed = record.trim_to_amplicon(hit).await?;
                match trimmed {
                    Some(record) => match record.whether_to_write(&filters).await {
                        true => writer.write_record(&record).await?,
                        false => continue,
                    },
                    _ => continue,
                }
            } else {
                continue;
            }
        }

        // Finalize the written contents to make sure the file is not corrupted
        format.finalize_write(writer).await?;

        Ok(())
    }
}

pub trait Sorting: SupportedFormat {
    fn sort_reads(self) -> impl Future<Output = Result<Self>>
    where
        Self: std::marker::Sized;
}

pub async fn sync_trimming<I>(reads: I, scheme: &AmpliconScheme) -> Result<Vec<FastqRecord>>
where
    I: IntoIterator<Item = FastqRecord>,
{
    // trim them down based on the amplicon scheme
    let reads = reads.into_iter().map(|record| async move {
        if let Some(hit) = record.find_amplicon(&scheme.scheme).await {
            let trimmed_record = record.trim_to_amplicon(hit).await?;
            Ok(trimmed_record)
        } else {
            Ok(None)
        }
    });

    let collected = join_all(reads)
        .await
        .into_iter()
        .filter_map(|result: Result<Option<FastqRecord>>| result.ok().flatten())
        .collect();

    Ok(collected)
}
