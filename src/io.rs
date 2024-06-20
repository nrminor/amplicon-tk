use crate::reads::Reads;
use color_eyre::eyre::Result;

pub trait SeqReader<'a> {
    fn read_fq(self) -> Result<Reads<'a>>;
    fn read_fa(self) -> Result<Reads<'a>>;
}

pub trait SeqWriter {
    fn write_fq(self) -> Result<()>;
    fn write_fa(self) -> Result<()>;
}
