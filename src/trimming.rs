use std::{pin::Pin, sync::Arc};

use futures::{io, TryStreamExt};
use noodles::fastq::Record as FastqRecord;
use tokio::runtime::Handle;
use tracing::info;

use crate::{amplicons::AmpliconScheme, prelude::RecordStream, record::FindAmplicons};

/// Trait `Trimming` is implemented for the various supported Record streams, and brings
/// with it the ability to find and trim sequences and quality scores down to only
/// the positions between forward and reverse primers. By default, it only retains reads
/// where a single primer-pair was identically matched. As such, `Trimming` also does
/// a form of filtering that is not handled by the `Filtering` trait.
pub trait Trimming: Sized {
    fn trim_to_amplicons(
        self,
        scheme: Arc<AmpliconScheme>,
    ) -> impl std::future::Future<Output = io::Result<Self>>;
}

impl<'a> Trimming for RecordStream<'a, FastqRecord> {
    async fn trim_to_amplicons(mut self, scheme: Arc<AmpliconScheme>) -> io::Result<Self> {
        let workers = Handle::current().metrics().num_workers();
        info!("{workers} worker threads allocated for processing records.");

        // use in-place mutation to trim the original record memory allocations
        let pinned_stream = Pin::new(&mut self);
        pinned_stream
            .project()
            .inner
            .as_mut()
            .try_for_each_concurrent(None, |record| {
                let scheme = Arc::clone(&scheme);
                async move {
                    let amplicon_hit = record.find_amplicon(&scheme.scheme).await;
                    if let Some(hit) = amplicon_hit {
                        // side-effect! in-place mutation
                        record.to_bounds(hit).await;
                    }
                    Ok(())
                }
            })
            .await?;

        Ok(self)
    }
}
