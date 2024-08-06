// #![warn(missing_docs)]

//!

use async_compression::tokio::write::GzipEncoder;
use futures::TryStreamExt;
use futures::{future::join_all, Future};
use noodles::fastq::AsyncWriter as FastqWriter;
use noodles::fastq::Record as FastqRecord;
use std::mem;
use std::path::Path;
use std::sync::Arc;
use tokio::{fs::File, io::BufWriter, sync::Mutex};
use tracing::info;

use crate::filtering::FilterSettings;
use crate::io::RecordParser;
use crate::{
    amplicons::AmpliconScheme,
    io::{Fastq, FastqGz, SeqWriter, SupportedFormat},
    record::FindAmplicons,
};
use color_eyre::eyre::Result;

pub trait Trimming<R, W>: SupportedFormat {
    fn trim(
        self,
        reader: R,
        output_path: &Path,
        scheme: Arc<AmpliconScheme>,
        filters: Arc<Option<FilterSettings>>,
    ) -> impl Future<Output = Result<()>>
    where
        R: RecordParser,
        for<'a, 'b> R::Record: FindAmplicons<'a, 'b> + Unpin;
}

impl<R, W> Trimming<R, W> for Fastq
where
    R: RecordParser,
    for<'a, 'b> R::Record: FindAmplicons<'a, 'b> + Unpin,
{
    async fn trim(
        self,
        mut reader: R,
        output_path: &Path,
        scheme: Arc<AmpliconScheme>,
        filters: Arc<Option<FilterSettings<'_, '_>>>,
    ) -> Result<()> {
        let records = reader.parse_records();
        let writer = self.read_writer(output_path).await?;
        let safe_writer = Arc::from(Mutex::from(writer));

        let handle = tokio::runtime::Handle::current();
        let workers = handle.metrics().num_workers();
        info!("{workers} worker threads allocated for processing records.");

        records
            .try_for_each_concurrent(workers, |record| {
                let scheme = Arc::clone(&scheme);
                let filters = Arc::clone(&filters);
                let _writer_instance = Arc::clone(&safe_writer);

                async move {
                    let amplicon_hit = record.find_amplicon(&scheme.scheme).await;
                    if let Some(hit) = amplicon_hit {
                        let trimmed = record.to_bounds(hit).await;
                        match trimmed.whether_to_write(&filters).await {
                            true => {
                                // let mut writer = writer_instance.lock().await;
                                // writer.write_trimmed(&trimmed).await?;
                                todo!();
                                // Ok(())
                            }
                            false => Ok(()),
                        }
                    } else {
                        Ok(())
                    }
                }
            })
            .await?;

        // Finalize the written contents to make sure the file is not corrupted
        let mut final_writer = safe_writer.lock().await;
        let extracted_writer = mem::replace(
            &mut *final_writer,
            FastqWriter::new(BufWriter::new(File::open(output_path).await?)),
        );
        drop(final_writer);
        let final_contents = extracted_writer.into_inner();
        self.finalize_write(final_contents).await?;

        Ok(())
    }
}

impl<R, W> Trimming<R, W> for FastqGz {
    // type Record = FastqRecord;
    async fn trim(
        self,
        mut reader: R,
        output_path: &Path,
        scheme: Arc<AmpliconScheme>,
        filters: Arc<Option<FilterSettings<'_, '_>>>,
    ) -> Result<()>
    where
        R: RecordParser,
        for<'a, 'b> R::Record: FindAmplicons<'a, 'b> + Unpin,
    {
        let records = reader.parse_records();
        let writer = self.read_writer(output_path).await?;
        let safe_writer = Arc::from(Mutex::from(writer));

        let handle = tokio::runtime::Handle::current();
        let workers = handle.metrics().num_workers();
        info!("{workers} worker threads allocated for processing records.");

        records
            .try_for_each_concurrent(workers, |record| {
                let scheme = Arc::clone(&scheme);
                let filters = Arc::clone(&filters);
                let _writer_instance = Arc::clone(&safe_writer);

                async move {
                    let amplicon_hit = record.find_amplicon(&scheme.scheme).await;
                    if let Some(hit) = amplicon_hit {
                        let trimmed = record.to_bounds(hit).await;
                        match trimmed.whether_to_write(&filters).await {
                            true => {
                                // let mut writer = writer_instance.lock().await;
                                // writer.write_trimmed(&trimmed).await?;
                                // Ok(())
                                todo!();
                            }
                            false => Ok(()),
                        }
                    } else {
                        Ok(())
                    }
                }
            })
            .await?;

        // Finalize the written contents to make sure the file is not corrupted
        let mut final_writer = safe_writer.lock().await;
        let extracted_writer = mem::replace(
            &mut *final_writer,
            FastqWriter::new(GzipEncoder::new(BufWriter::new(
                File::open(output_path).await?,
            ))),
        );
        drop(final_writer);
        let final_contents = extracted_writer.into_inner();
        self.finalize_write(final_contents).await?;

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
            let trimmed_record = record.to_bounds(hit).await;
            Ok(Some(trimmed_record))
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
