// #![warn(missing_docs)]

//!

use futures::{future::join_all, Future};
use noodles::fastq::Record as FastqRecord;
use std::path::Path;

use crate::{
    io::{Fastq, FastqGz, Init, SeqWriter, SupportedFormat},
    primers::AmpliconScheme,
    record::{trim_records, FindAmplicons},
};
use color_eyre::eyre::Result;
use futures::TryStreamExt;

pub trait Trimming: SupportedFormat {
    fn run_trimming(
        self,
        input_path: &Path,
        output_path: &Path,
        scheme: AmpliconScheme,
    ) -> impl Future<Output = Result<()>>;
}

impl Trimming for Fastq {
    async fn run_trimming(
        self,
        input_path: &Path,
        output_path: &Path,
        scheme: AmpliconScheme,
    ) -> Result<()> {
        let (mut reader, format) = self.init(input_path).await?;
        let mut records = reader.records();

        // iterate through records asynchronously, find amplicon hits, and trim them down to exclude
        // primers and anything that extends beyond them
        while let Some(mut record) = records.try_next().await? {
            let amplicon_hit = record.amplicon(&scheme.scheme);
            if let Some(hit) = amplicon_hit {
                trim_records(&mut record, hit).await?;
            } else {
                continue;
            }
        }

        // inititialize a writer and finalize the contents going into this (this currently does
        // nothing but demonstrates the syntax necessary to finish a write job)
        let writer = format.read_writer(output_path).await?;
        format.finalize_write(writer).await?;

        Ok(())
    }
}

impl Trimming for FastqGz {
    async fn run_trimming(
        self,
        input_path: &Path,
        output_path: &Path,
        scheme: AmpliconScheme,
    ) -> Result<()> {
        let (mut reader, format) = self.init(input_path).await?;
        let mut records = reader.records();

        // iterate through records asynchronously, find amplicon hits, and trim them down to exclude
        // primers and anything that extends beyond them
        while let Some(mut record) = records.try_next().await? {
            let amplicon_hit = record.amplicon(&scheme.scheme);
            if let Some(hit) = amplicon_hit {
                trim_records(&mut record, hit).await?;
            } else {
                continue;
            }
        }

        // inititialize a writer and finalize the contents going into this (this currently does
        // nothing but demonstrates the syntax necessary to finish a write job)
        let writer = format.read_writer(output_path).await?;
        format.finalize_write(writer).await?;

        Ok(())
    }
}

pub trait Filtering: SupportedFormat {
    fn run_filters(&mut self) -> impl Future<Output = Result<Self>>
    where
        Self: std::marker::Sized;
}

pub trait Sorting: SupportedFormat {
    fn sort_reads(&mut self) -> impl Future<Output = Result<Self>>
    where
        Self: std::marker::Sized;
}

pub async fn sync_trimming<I>(reads: I, scheme: &AmpliconScheme) -> Result<Vec<FastqRecord>>
where
    I: IntoIterator<Item = FastqRecord>,
{
    // trim them down based on the amplicon scheme
    let reads = reads.into_iter().map(|mut record| async move {
        if let Some(hit) = record.amplicon(&scheme.scheme) {
            let _ = trim_records(&mut record, hit).await?;
        }
        Ok(record)
    });

    join_all(reads).await.into_iter().collect()
}

// pub async fn trim_bam_records(
//     input_path: &Path,
//     output_path: &Path,
//     scheme: AmpliconScheme<'_>,
// ) -> Result<()> {
//     let (mut reader, format) = init_bam(input_path).await?;
//     let mut records = reader.records();

//     // iterate through records asynchronously, find amplicon hits, and trim them down to exclude
//     // primers and anything that extends beyond them
//     while let Some(mut record) = records.try_next().await? {
//         let amplicon_hit = record.amplicon(&scheme.scheme);
//         if let Some(hit) = amplicon_hit {
//             trim_records(&mut record, hit).await?;
//         } else {
//             continue;
//         }
//     }

//     // inititialize a writer and finalize the contents going into this (this currently does
//     // nothing but demonstrates the syntax necessary to finish a write job)
//     let writer = format.read_writer(output_path).await?;
//     format.finalize_write(writer).await?;

//     Ok(())
// }
