//! Module `io` contains traits and trait implementations to handle the reading and writing
//! of various bioinformatic file formats, including FASTQ files, FASTA files, and BED files.
//! Where possible, care has been taken to make these readers and writers generic over possible
//! readers and writers, e.g., over compressed or uncompressed inputs and outputs.
//!
//! Traits exported here rely heavily on the parsers implemented in the excellend `noodles` crate.

// #![warn(missing_docs)]

use std::path::Path;

use async_compression::tokio::bufread::GzipDecoder;
use async_compression::tokio::write::GzipEncoder;
use color_eyre::eyre::eyre;
use color_eyre::eyre::Result;
use noodles::bam::AsyncReader as BamReader;
use noodles::bam::AsyncWriter as BamWriter;
use noodles::bed::io::Reader as BedReader;
use noodles::bgzf::AsyncReader as BgzfReader;
use noodles::bgzf::AsyncWriter as BgzfWriter;
use noodles::fasta::io::Reader as FastaReader;
use noodles::fastq::AsyncReader as FastqReader;
use noodles::fastq::AsyncWriter as FastqWriter;
use tokio::io::AsyncWriteExt;
use tokio::io::BufWriter;
use tokio::{fs::File, io::BufReader};

// supported sequencing read formats
pub struct FastqGz;
pub struct Fastq;
pub struct Bam;

pub enum InputType {
    FASTQGZ(FastqGz),
    FASTQ(Fastq),
    BAM(Bam),
}

impl InputType {
    pub fn extension(&self) -> String {
        match self {
            InputType::FASTQGZ(_) => String::from(".fastq.gz"),
            InputType::FASTQ(_) => String::from(".fastq"),
            InputType::BAM(_) => String::from(".bam"),
        }
    }
}

pub enum OutputType {
    FASTQGZ(FastqGz),
    FASTQ(Fastq),
    BAM(Bam),
}

// supported input primer and reference formats
pub struct Bed;
pub struct Fasta;
pub struct Genbank;

pub enum PrimerType {
    BED,
}

// implementing marker traits to constrain which formats are representable
pub trait SupportedFormat {}
impl SupportedFormat for FastqGz {}
impl SupportedFormat for Fastq {}
impl SupportedFormat for Bam {}

pub trait PrimerFormat {}
impl PrimerFormat for Bed {}

pub trait RefFormat {}
impl RefFormat for Fasta {}

pub trait SeqReader {
    type Format: SupportedFormat;
    type Reader: Unpin + Send;
    fn read_reads(&self, input_path: &Path) -> impl futures::Future<Output = Result<Self::Reader>>;
}

impl SeqReader for FastqGz {
    type Format = FastqGz;
    type Reader = FastqReader<BufReader<GzipDecoder<BufReader<File>>>>;
    async fn read_reads(&self, input_path: &Path) -> Result<Self::Reader> {
        let input_file = File::open(input_path).await?;
        let reader = BufReader::new(input_file);
        let decoder = GzipDecoder::new(reader);
        let decode_reader = BufReader::new(decoder);
        let fastq = FastqReader::new(decode_reader);

        Ok(fastq)
    }
}

impl SeqReader for Fastq {
    type Format = Fastq;
    type Reader = FastqReader<BufReader<File>>;
    async fn read_reads(&self, input_path: &Path) -> Result<Self::Reader> {
        let input_file = File::open(input_path).await?;
        let reader = BufReader::new(input_file);
        let fastq = FastqReader::new(reader);

        Ok(fastq)
    }
}

impl SeqReader for Bam {
    type Format = Bam;
    type Reader = BamReader<BgzfReader<File>>;
    async fn read_reads(&self, input_path: &Path) -> Result<Self::Reader> {
        let input_file = File::open(input_path).await?;
        let bam = BamReader::new(input_file);

        Ok(bam)
    }
}

pub trait PrimerReader {
    type Format: PrimerFormat;
    type Reader;
    fn read_primers(&self, input_path: &Path) -> Result<Self::Reader>;
}

impl PrimerReader for Bed {
    type Format = Bed;
    type Reader = BedReader<std::io::BufReader<std::fs::File>>;
    fn read_primers(&self, input_path: &Path) -> Result<Self::Reader> {
        let reader = std::fs::File::open(input_path)
            .map(std::io::BufReader::new)
            .map(noodles::bed::io::Reader::new)?;

        Ok(reader)
    }
}

pub trait RefReader: RefFormat {
    type Reader;
    fn read_ref(&self, input_path: &Path) -> Result<Self::Reader>;
}

impl RefReader for Fasta {
    type Reader = FastaReader<std::io::BufReader<std::fs::File>>;
    fn read_ref(&self, input_path: &Path) -> Result<Self::Reader> {
        let reader = std::fs::File::open(input_path)
            .map(std::io::BufReader::new)
            .map(noodles::fasta::io::Reader::new)?;

        Ok(reader)
    }
}

pub trait SeqWriter: SupportedFormat {
    type Writer: Unpin + Send;
    fn read_writer(
        &self,
        output_file_path: &Path,
    ) -> impl futures::Future<Output = Result<Self::Writer>>;
    fn finalize_write(&self, writer: Self::Writer) -> impl futures::Future<Output = Result<()>>;
}

impl SeqWriter for FastqGz {
    type Writer = FastqWriter<GzipEncoder<BufWriter<File>>>;
    async fn read_writer(&self, output_file_path: &Path) -> Result<Self::Writer> {
        let output_file = File::create(output_file_path).await?;
        let writer = BufWriter::new(output_file);
        let encoder = GzipEncoder::new(writer);
        let fastq_writer = FastqWriter::new(encoder);

        Ok(fastq_writer)
    }
    async fn finalize_write(&self, writer: Self::Writer) -> Result<()> {
        let mut final_contents = writer.into_inner();
        final_contents.shutdown().await?;
        Ok(())
    }
}

impl SeqWriter for Fastq {
    type Writer = FastqWriter<BufWriter<File>>;
    async fn read_writer(&self, output_path: &Path) -> Result<Self::Writer> {
        let output_file = File::create(output_path).await?;
        let writer = BufWriter::new(output_file);
        let fastq_writer = FastqWriter::new(writer);

        Ok(fastq_writer)
    }
    async fn finalize_write(&self, writer: Self::Writer) -> Result<()> {
        let mut final_contents = writer.into_inner();
        final_contents.flush().await?;
        Ok(())
    }
}

impl SeqWriter for Bam {
    type Writer = BamWriter<BgzfWriter<File>>;
    async fn read_writer(&self, output_path: &Path) -> Result<Self::Writer> {
        let output_file = File::create(output_path).await?;
        let bam_writer = BamWriter::new(output_file);

        Ok(bam_writer)
    }
    async fn finalize_write(&self, writer: Self::Writer) -> Result<()> {
        let mut final_contents = writer.into_inner();
        final_contents.shutdown().await?;
        Ok(())
    }
}

pub async fn io_selector(input_path: &Path) -> Result<InputType> {
    match input_path.try_exists() {
        Ok(_) => (),
        Err(_) => return Err(eyre!("The provided file {:?} does not exist.", input_path)),
    }

    let extension = input_path.extension();
    if let Some(ext) = extension {
        match ext.to_str().unwrap_or("") {
            "gz" => Ok(InputType::FASTQGZ(FastqGz)),
            "fastq" => Ok(InputType::FASTQ(Fastq)),
            "bam" => Ok(InputType::BAM(Bam)),
            _ => Err(eyre!("Unsupported file type provided: {:?}", input_path)),
        }
    } else {
        Err(eyre!(
            "Could not determine an extension from the provided file name: {:?}.",
            input_path
        ))
    }
}

pub trait Init: SupportedFormat {
    type Reader;
    fn init(self, input_path: &Path) -> impl futures::Future<Output = Result<(Self::Reader, Self)>>
    where
        Self: std::marker::Sized;
}

impl Init for FastqGz {
    type Reader = FastqReader<BufReader<GzipDecoder<BufReader<File>>>>;
    async fn init(self, input_path: &Path) -> Result<(Self::Reader, Self)>
    where
        Self: std::marker::Sized,
    {
        let reader = self.read_reads(input_path).await?;
        Ok((reader, self))
    }
}

impl Init for Fastq {
    type Reader = FastqReader<BufReader<File>>;
    async fn init(self, input_path: &Path) -> Result<(Self::Reader, Self)>
    where
        Self: std::marker::Sized,
    {
        let reader = self.read_reads(input_path).await?;
        Ok((reader, self))
    }
}
