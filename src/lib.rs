use color_eyre::eyre::Result;
use per_record::Record;
use primers::PrimerPair;

pub mod cli;
pub mod per_record;
pub mod primers;
pub mod py_api;
pub mod reads;

// #![warn(missing_docs)]

pub trait Trim<'a> {
    fn trim_fastq_to_primers(self, primers: &'a PrimerPair) -> Result<Option<Record<'a>>>;
    fn trim_fasta_to_primers(self, primers: &'a PrimerPair) -> Result<Option<Record<'a>>>;
}
