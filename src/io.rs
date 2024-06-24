use std::{fs::File, io::BufReader, path::Path};

use color_eyre::eyre::Result;

pub trait SeqReader {
    fn read_fq(input_path: &Path) -> Result<Self>
    where
        Self: std::marker::Sized;
    // fn read_fa(input_path: &Path) -> Result<Self>
    // where
    //     Self: std::marker::Sized;
}

pub trait SeqWriter {
    fn write_fq(self) -> Result<()>;
    fn write_fa(self) -> Result<()>;
}

pub trait PrimerReader {
    fn read_bed(input_path: &Path) -> Result<Self>
    where
        Self: std::marker::Sized;
}

pub trait RefReader {
    fn read_ref(input_path: &Path) -> Result<Self>
    where
        Self: std::marker::Sized;
}

impl SeqReader for noodles::fastq::Reader<BufReader<File>> {
    fn read_fq(input_path: &Path) -> Result<Self> {
        let reader = File::open(input_path)
            .map(BufReader::new)
            .map(noodles::fastq::io::Reader::new)?;

        Ok(reader)
    }
}

// impl SeqWriter for

impl PrimerReader for noodles::bed::Reader<BufReader<File>> {
    fn read_bed(input_path: &Path) -> Result<Self> {
        let reader = File::open(input_path)
            .map(BufReader::new)
            .map(noodles::bed::io::Reader::new)?;

        Ok(reader)
    }
}

impl RefReader for noodles::fasta::Reader<BufReader<File>> {
    fn read_ref(input_path: &Path) -> Result<Self> {
        let reader = File::open(input_path)
            .map(BufReader::new)
            .map(noodles::fasta::io::Reader::new)?;

        Ok(reader)
    }
}
