//! Module `io` contains traits and trait implementations to handle the reading and writing
//! of various bioinformatic file formats, including FASTQ files, FASTA files, and BED files.
//! Where possible, care has been taken to make these readers and writers generic over possible
//! readers and writers, e.g., over compressed or uncompressed inputs and outputs.
//!
//! Traits exported here rely heavily on the parses implemented in the excellend `noodles` crate.

#![warn(missing_docs)]

use std::io::Write;
use std::{fs::File, io::BufReader, path::Path};

use color_eyre::eyre::Result;

use noodles::bed::Reader as BedReader;
use noodles::fasta::Reader as FaReader;
use noodles::fastq::io::Writer as FqWriter;
use noodles::fastq::Reader as FqReader;

use crate::reads::Reads;

///
pub trait SeqReader {
    ///
    fn read_fq(input_path: &Path) -> Result<Self>
    where
        Self: std::marker::Sized;
}

///
pub trait PrimerReader {
    ///
    fn read_bed(input_path: &Path) -> Result<Self>
    where
        Self: std::marker::Sized;
}

///
pub trait RefReader {
    ///
    fn read_ref(input_path: &Path) -> Result<Self>
    where
        Self: std::marker::Sized;
}

///
pub trait SeqWriter {
    ///
    fn write_fq(self, fq_dataset: Reads) -> Result<()>;
}

impl SeqReader for FqReader<BufReader<File>> {
    fn read_fq(input_path: &Path) -> Result<Self> {
        let reader = File::open(input_path)
            .map(BufReader::new)
            .map(noodles::fastq::io::Reader::new)?;

        Ok(reader)
    }
}

impl PrimerReader for BedReader<BufReader<File>> {
    fn read_bed(input_path: &Path) -> Result<Self> {
        let reader = File::open(input_path)
            .map(BufReader::new)
            .map(noodles::bed::io::Reader::new)?;

        Ok(reader)
    }
}

impl RefReader for FaReader<BufReader<File>> {
    fn read_ref(input_path: &Path) -> Result<Self> {
        let reader = File::open(input_path)
            .map(BufReader::new)
            .map(noodles::fasta::io::Reader::new)?;

        Ok(reader)
    }
}

impl<W: Write> SeqWriter for FqWriter<W> {
    fn write_fq(mut self, fq_dataset: Reads) -> Result<()> {
        for record in fq_dataset.reads {
            self.write_record(&record)?;
        }
        Ok(())
    }
}
