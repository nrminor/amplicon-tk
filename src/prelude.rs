//! `prelude` contains the core traits, implementations, and structs for file I/O and
//! record streaming.

use std::mem;
use std::sync::Arc;
use std::{fmt, io, path::Path};
use tokio::{io::AsyncWriteExt, sync::Mutex};

use async_compression::tokio::{bufread::GzipDecoder, write::GzipEncoder};
use clap::ValueEnum;
use color_eyre::eyre::Result;
use futures::{Stream, TryStreamExt};
use noodles::{bam::Record as BamRecord, fastq::Record as FastqRecord};

use noodles::fastq::{AsyncReader as FastqReader, AsyncWriter as FastqWriter};
use pin_project_lite::pin_project;
use tokio::{
    fs::File,
    io::{BufReader, BufWriter},
};

#[derive(ValueEnum, Debug, Clone, PartialEq)]
///
pub enum SupportedTypes {
    /// Read from Gzip- or BGzip-compressed FASTQ files.
    FASTQGZ,

    /// Read from uncompressed FASTQ files.
    FASTQ,

    /// Read from BAM files.
    BAM,
}

///
impl fmt::Display for SupportedTypes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                SupportedTypes::FASTQGZ => ".fastq.gz",
                SupportedTypes::FASTQ => ".fastq",
                SupportedTypes::BAM => ".bam",
            }
        )
    }
}

impl SupportedTypes {
    ///
    pub fn from_file_name(file_name: &Path) -> Option<Self> {
        if let Some(extension) = file_name.extension() {
            match extension.to_str().unwrap_or("") {
                "gz" => {
                    if file_name.to_str().unwrap_or("").ends_with(".fastq.gz") {
                        return Some(SupportedTypes::FASTQGZ);
                    }
                }
                "fastq" => return Some(SupportedTypes::FASTQ),
                "bam" => return Some(SupportedTypes::BAM),
                _ => return None,
            }
        }
        None
    }
}

/// `SeqReader` is the first of a few container types used in `amplicon-tk` to generically
/// support multiple input data formats. With a few trait bounds, `SeqReader` uses
/// parametric polymorphism to support containing readers for Gzipped FASTQs,
/// uncompressed FASTQs, and BAM files, with support for more formats coming in the
/// future.
pub struct SeqReader<R>
where
    R: Unpin + Send + 'static,
{
    pub inner: R,
}

impl SeqReader<FastqReader<BufReader<GzipDecoder<BufReader<File>>>>> {
    ///
    pub async fn new_fastq_gz(input_path: &Path) -> Result<Self> {
        let input_file = File::open(input_path).await?;
        let reader = BufReader::new(input_file);
        let decoder = GzipDecoder::new(reader);
        let decode_reader = BufReader::new(decoder);
        let full_reader = SeqReader {
            inner: FastqReader::new(decode_reader),
        };

        Ok(full_reader)
    }

    ///
    pub fn records(&mut self) -> impl Stream<Item = io::Result<FastqRecord>> + '_ {
        self.inner.records()
    }
}

impl SeqReader<FastqReader<BufReader<File>>> {
    ///
    pub async fn new_fastq(input_path: &Path) -> Result<Self> {
        let input_file = File::open(input_path).await?;
        let reader = BufReader::new(input_file);
        let full_reader = SeqReader {
            inner: FastqReader::new(reader),
        };

        Ok(full_reader)
    }

    ///
    pub fn records(&mut self) -> impl Stream<Item = io::Result<FastqRecord>> + '_ {
        self.inner.records()
    }
}

/// `parseable` constrains which bioinformatic data formats can be processed and offers a
/// run_filters getter. Future versions will bring in convenience sequence name and sequence
/// bases methods.
pub trait Parseable: Sized {
    fn get(self) -> Self;
}

impl Parseable for FastqRecord {
    fn get(self) -> Self {
        self
    }
}

impl Parseable for BamRecord {
    fn get(self) -> Self {
        self
    }
}

pin_project! {
    /// RecordStream is the core container type used to make fluent interfaces between the
    /// several steps in `amplicon-tk`. Ultimately, it contains a lazy, asynchronous stream
    /// of bioinformatic data entries or "records". `RecordStream` is generic to the input
    /// data type and can be used with FASTQ entries or BAM entries, with potential for
    /// more types in the future.
    pub struct RecordStream<'a, P>
    where
        P: Parseable,
    {
        #[pin]
        pub inner: Box<dyn Stream<Item = io::Result<P>> + Send + Unpin + 'a>,
    }
}

impl<'a, P: Parseable> RecordStream<'a, P> {
    ///
    pub fn new<S>(stream: S) -> Self
    where
        S: Stream<Item = io::Result<P>> + Send + Unpin + 'a,
    {
        RecordStream {
            inner: Box::new(stream),
        }
    }
}

impl<'a> RecordStream<'a, FastqRecord> {
    ///
    pub async fn from_fastq(
        records: impl Stream<Item = io::Result<FastqRecord>> + Send + Unpin + 'a,
    ) -> io::Result<Self> {
        Ok(Self::new(records))
    }
}

///
pub trait SeqWriter {
    fn write_records(
        self,
        output_type: SupportedTypes,
        output_path: &Path,
    ) -> impl std::future::Future<Output = Result<()>>;
}

impl<'a> SeqWriter for RecordStream<'a, FastqRecord> {
    async fn write_records(self, output_type: SupportedTypes, output_path: &Path) -> Result<()> {
        match output_type {
            SupportedTypes::BAM => todo!(),
            SupportedTypes::FASTQ => self.write_fastq(output_path).await,
            SupportedTypes::FASTQGZ => self.write_fastq_gz(output_path).await,
        }
    }
}

impl<'a> RecordStream<'a, FastqRecord> {
    async fn write_fastq_gz(self, output_path: &Path) -> Result<()> {
        let output_file = File::create(output_path).await?;
        let writer = BufWriter::new(output_file);
        let encoder = GzipEncoder::new(writer);
        let fastq_writer = FastqWriter::new(encoder);
        let safe_writer = Arc::from(Mutex::from(fastq_writer));

        self.inner
            .try_for_each(|record| {
                let writer_instance = Arc::clone(&safe_writer);
                async move {
                    let mut writer = writer_instance.lock().await;
                    writer.write_record(&record).await?;
                    Ok(())
                }
            })
            .await?;

        //
        let mut final_writer = safe_writer.lock().await;
        let extracted_writer = mem::replace(
            &mut *final_writer,
            FastqWriter::new(GzipEncoder::new(BufWriter::new(
                File::open(output_path).await?,
            ))),
        );
        drop(final_writer);
        let mut final_contents = extracted_writer.into_inner();
        final_contents.flush().await?;
        final_contents.shutdown().await?;

        Ok(())
    }
    async fn write_fastq(self, output_path: &Path) -> Result<()> {
        let output_file = File::create(output_path).await?;
        let writer = BufWriter::new(output_file);
        let fastq_writer = FastqWriter::new(writer);
        let safe_writer = Arc::from(Mutex::from(fastq_writer));

        self.inner
            .try_for_each(|record| {
                let writer_instance = Arc::clone(&safe_writer);
                async move {
                    let mut writer = writer_instance.lock().await;
                    writer.write_record(&record).await?;
                    Ok(())
                }
            })
            .await?;

        //
        let mut final_writer = safe_writer.lock().await;
        let extracted_writer = mem::replace(
            &mut *final_writer,
            FastqWriter::new(BufWriter::new(File::open(output_path).await?)),
        );
        drop(final_writer);
        let mut final_contents = extracted_writer.into_inner();
        final_contents.flush().await?;
        final_contents.shutdown().await?;

        Ok(())
    }
}
