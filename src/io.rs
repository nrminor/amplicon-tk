use std::{fs::File, io::BufReader, path::Path};

use crate::reads::Reads;
use color_eyre::eyre::Result;

pub trait SeqReader<'a, R> {
    fn read_fq(self) -> Result<Reads>;
    fn read_fa(self) -> Result<Reads>;
}

/// .
///
/// # Errors
///
/// This function will return an error if .
pub fn read_fq(input_path: &Path) -> Result<noodles::fastq::Reader<BufReader<File>>> {
    let reader = File::open(input_path)
        .map(BufReader::new)
        .map(noodles::fastq::io::Reader::new)?;

    Ok(reader)
}

pub trait SeqWriter {
    fn write_fq(self) -> Result<()>;
    fn write_fa(self) -> Result<()>;
}
