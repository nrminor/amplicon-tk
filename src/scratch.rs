use std::fmt;
use std::{
    io,
    path::{Path, PathBuf},
};

use async_compression::tokio::bufread::GzipDecoder;
use clap::{Parser, ValueEnum};
use color_eyre::eyre::{eyre, Result};
use futures::Stream;
use noodles::bam::Record as BamRecord;
use noodles::fastq::Record as FastqRecord;

use noodles::fastq::AsyncReader as FastqReader;
use tokio::{fs::File, io::BufReader};

use crate::cli::Commands;

#[derive(ValueEnum, Debug, Clone, PartialEq)]
pub enum SupportedTypes {
    /// Read from Gzip- or BGzip-compressed FASTQ files.
    FASTQGZ,

    /// Read from uncompressed FASTQ files.
    FASTQ,

    /// Read from BAM files.
    BAM,
}

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
    async fn new_fastq_gz(input_path: &Path) -> Result<Self> {
        let input_file = File::open(input_path).await?;
        let reader = BufReader::new(input_file);
        let decoder = GzipDecoder::new(reader);
        let decode_reader = BufReader::new(decoder);
        let full_reader = SeqReader {
            inner: FastqReader::new(decode_reader),
        };

        Ok(full_reader)
    }

    fn records(&mut self) -> impl Stream<Item = io::Result<FastqRecord>> + '_ {
        self.inner.records()
    }
}

impl SeqReader<FastqReader<BufReader<File>>> {
    async fn new_fastq(input_path: &Path) -> Result<Self> {
        let input_file = File::open(input_path).await?;
        let reader = BufReader::new(input_file);
        let full_reader = SeqReader {
            inner: FastqReader::new(reader),
        };

        Ok(full_reader)
    }

    fn records(&mut self) -> impl Stream<Item = io::Result<FastqRecord>> + '_ {
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

/// RecordStream is the core container type used to make fluent interfaces between the
/// several steps in `amplicon-tk`. Ultimately, it contains a lazy, asynchronous stream
/// of bioinformatic data entries or "records". `RecordStream` is generic to the input
/// data type and can be used with FASTQ entries or BAM entries, with potential for
/// more types in the future.
pub struct RecordStream<'a, P>
where
    P: Parseable,
{
    #[allow(dead_code)]
    inner: Box<dyn Stream<Item = io::Result<P>> + Send + Unpin + 'a>,
}

impl<'a, P: Parseable> RecordStream<'a, P> {
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
    pub async fn from_fastq(
        records: impl Stream<Item = io::Result<FastqRecord>> + Send + Unpin + 'a,
    ) -> io::Result<Self> {
        Ok(Self::new(records))
    }
}

/// Trait `Trimming` is implemented for the various supported Record streams, and brings
/// with it the ability to find and trim sequences and quality scores down to only
/// the positions between forward and reverse primers. By default, it only retains reads
/// where a single primer-pair was identically matched. As such, `Trimming` also does
/// a form of filtering that is not handled by the `Filtering` trait.
pub trait Trimming: Sized {
    type RecordType: Parseable;

    fn trim_to_amplicons(self) -> impl std::future::Future<Output = io::Result<Self>>;
}

impl<'a> Trimming for RecordStream<'a, FastqRecord> {
    type RecordType = FastqRecord;
    async fn trim_to_amplicons(self) -> io::Result<Self> {
        Ok(self)
    }
}

/// Trait `Filtering` handles special filtering requests that require context from the
/// entire sequencing read dataset, e.g., filtering by frequency of each unique sequence.
/// Unlike trait `Trimming`, `Filtering` is not intrinsically amplicon-aware. Instead,
/// it pulls the information it needs from a required index file that can be computed
/// ahead of time with `amplicon-tk index`.
pub trait Filtering: Sized {
    type RecordType: Parseable;

    fn run_filters(self) -> impl std::future::Future<Output = io::Result<Self>>;
}

impl<'a> Filtering for RecordStream<'a, FastqRecord> {
    type RecordType = FastqRecord;
    async fn run_filters(self) -> io::Result<Self> {
        Ok(self)
    }
}

#[allow(dead_code)]
async fn _main() -> Result<()> {
    let cli = crate::cli::Cli::parse();

    #[allow(unused_variables)]
    match &cli.command {
        Some(Commands::Trim {
            input_file,
            bed_file,
            fasta_ref,
            keep_multi,
            left_suffix,
            right_suffix,
            min_freq,
            expected_len,
            output,
        }) => {
            todo!()
        }
        _ => eprintln!("ignored"),
    }

    // this is starting to come together as a very nice looking
    // fluent interface
    let test_fastq = PathBuf::from("test.fastq");

    // parse the supplied filetype
    let file_type = SupportedTypes::from_file_name(&test_fastq);

    // early return if it's unsupported
    if file_type.is_none() {
        return Err(eyre!("Unsupported file type supplied in {:?}", file_type));
    }

    // use pattern matching to process the reads based on their record type
    match file_type.unwrap_or(SupportedTypes::FASTQ) {
        SupportedTypes::FASTQ => {
            let mut reader = SeqReader::new_fastq(&test_fastq).await?;
            let records = reader.records();
            let _results = RecordStream::from_fastq(records)
                .await?
                .trim_to_amplicons()
                .await?
                .run_filters()
                .await?;
        }
        SupportedTypes::FASTQGZ => {
            let mut reader = SeqReader::new_fastq_gz(&test_fastq).await?;
            let records = reader.records();
            let _results = RecordStream::from_fastq(records)
                .await?
                .trim_to_amplicons()
                .await?
                .run_filters()
                .await?;
        }
        _ => eprintln!("BAM and other files are not yet supported."),
    }

    Ok(())
}
