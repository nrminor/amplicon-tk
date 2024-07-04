//! Module `io` contains traits and trait implementations to handle the reading and writing
//! of various bioinformatic file formats, including FASTQ files, FASTA files, and BED files.
//! Where possible, care has been taken to make these readers and writers generic over possible
//! readers and writers, e.g., over compressed or uncompressed inputs and outputs.
//!
//! Traits exported here rely heavily on the parsers implemented in the excellend `noodles` crate.

// #![warn(missing_docs)]

use std::path::Path;
use std::path::PathBuf;

use async_compression::tokio::bufread::GzipDecoder;
use async_compression::tokio::write::GzipEncoder;
use color_eyre::eyre::eyre;
use color_eyre::eyre::Result;
use futures::TryStreamExt;
use noodles::bam::AsyncReader as BamReader;
use noodles::bed::io::Reader as BedReader;
use noodles::bgzf::AsyncReader as BgzfReader;
use noodles::cram::AsyncReader as CramReader;
use noodles::fasta::io::Reader as FastaReader;
use noodles::fastq::AsyncReader as FastqReader;
use noodles::fastq::AsyncWriter as FastqWriter;
use noodles::sam::AsyncReader as SamReader;
use tokio::io::AsyncWriteExt;
use tokio::io::BufWriter;
use tokio::{fs::File, io::BufReader};

use crate::primers::define_amplicons;
use crate::primers::ref_to_dict;
use crate::record::trim_fq_records;
use crate::record::FindAmplicons;

pub enum InputType {
    FASTQGZ,
    FASTQ,
    BAM,
    SAM,
    CRAM,
}

pub enum PrimerType {
    BED,
}

pub enum OutputType {
    FASTQGZ,
    FASTQ,
    BAM,
    SAM,
    CRAM,
}

// supported sequencing read formats
pub struct FastqGz;
pub struct Fastq;
pub struct Bam;
pub struct Sam;
pub struct Cram;

// supported input primer and reference formats
pub struct Bed;
pub struct Fasta;
pub struct Genbank;

// implementing a marker trait to constrain which formats are representable
pub trait SupportedFormat {}
impl SupportedFormat for FastqGz {}
impl SupportedFormat for Fastq {}
impl SupportedFormat for Bam {}
impl SupportedFormat for Sam {}
impl SupportedFormat for Cram {}

pub trait PrimerFormat {}
impl PrimerFormat for Bed {}

pub trait RefFormat {}
impl RefFormat for Fasta {}

pub trait SeqReader {
    type Format: SupportedFormat;
    type Reader: Unpin + Send;
    async fn read_reads(&self, input_path: &Path) -> Result<Self::Reader>;
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

impl SeqReader for Sam {
    type Format = Sam;
    type Reader = SamReader<BufReader<File>>;
    async fn read_reads(&self, input_path: &Path) -> Result<Self::Reader> {
        let input_file = File::open(input_path).await?;
        let reader = BufReader::new(input_file);
        let sam = SamReader::new(reader);

        Ok(sam)
    }
}

impl SeqReader for Cram {
    type Format = Cram;
    type Reader = CramReader<File>;
    async fn read_reads(&self, input_path: &Path) -> Result<Self::Reader> {
        let input_file = File::open(input_path).await?;
        let cram = CramReader::new(input_file);

        Ok(cram)
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

pub trait RefReader {
    type Format: RefFormat;
    type Reader;
    fn read_ref(&self, input_path: &Path) -> Result<Self::Reader>;
}

impl RefReader for Fasta {
    type Format = Fasta;
    type Reader = FastaReader<std::io::BufReader<std::fs::File>>;
    fn read_ref(&self, input_path: &Path) -> Result<Self::Reader> {
        let reader = std::fs::File::open(input_path)
            .map(std::io::BufReader::new)
            .map(noodles::fasta::io::Reader::new)?;

        Ok(reader)
    }
}

pub trait SeqWriter {
    type Format: SupportedFormat;
    type Writer: Unpin + Send;
    async fn read_writer(&self, output_file_path: &Path) -> Result<Self::Writer>;
    async fn finalize_write(&self, writer: Self::Writer) -> Result<()>;
}

impl SeqWriter for FastqGz {
    type Format = FastqGz;
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
    type Format = Fastq;
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

async fn _test(
    input_path: &Path,
    output_path: &Path,
    input_type: InputType,
    output_type: OutputType,
) -> Result<()> {
    // defining input and output types for the reads
    let input_type = match input_type {
        InputType::FASTQGZ => {
            let file_type = FastqGz;
            Ok(file_type)
        }
        _ => Err(eyre!("Other formats than .fastq.gz not yet supported!")),
    }?;
    let output_type = match output_type {
        OutputType::FASTQGZ => {
            let file_type = FastqGz;
            Ok(file_type)
        }
        _ => Err(eyre!("Other formats than .fastq.gz not yet supported!")),
    }?;

    // pulling in the primers
    let primer_type = Bed;
    let test_bed_path = PathBuf::from("test.bed");
    let bed = primer_type.read_primers(&test_bed_path)?;
    let left_suffix = "_LEFT";
    let right_suffix = "_RIGHT";

    // pulling in the reference
    let ref_type = Fasta;
    let test_ref_path = PathBuf::from("test.fasta");
    let mut fasta = ref_type.read_ref(&test_ref_path)?;

    // convert the reference to a hashmap and use it to pull in the primer pairs for each
    // amplicon
    let ref_dict = ref_to_dict(&mut fasta).await?;
    let scheme = define_amplicons(bed, &ref_dict, left_suffix, right_suffix).await?;

    // initializing a stream that will hold the reads to be processed lazily and asynchronously
    let mut reader = input_type.read_reads(input_path).await?;
    let mut records = reader.records();

    // iterate through records asynchronously, find amplicon hits, and trim them down to exclude
    // primers and anything that extends beyond them
    while let Some(mut record) = records.try_next().await? {
        let amplicon_hit = record.amplicon(&scheme.scheme);
        if let Some(hit) = amplicon_hit {
            trim_fq_records(&mut record, hit).await?;
        } else {
            continue;
        }
    }

    // inititialize a writer and finalize the contents going into this (this currently does
    // nothing but demonstrates the syntax necessary to finish a write job)
    let writer = output_type.read_writer(output_path).await?;
    output_type.finalize_write(writer).await?;

    Ok(())
}
