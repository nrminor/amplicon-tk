use crate::reads::Reads;
use color_eyre::eyre::Result;

pub trait SeqReader<'a, R> {
    fn read_fq(self) -> Result<Reads<R>>;
    fn read_fa(self) -> Result<Reads<R>>;
}

pub trait SeqWriter {
    fn write_fq(self) -> Result<()>;
    fn write_fa(self) -> Result<()>;
}
