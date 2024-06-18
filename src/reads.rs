use color_eyre::eyre::Result;
use needletail::parser::SequenceRecord;

pub struct SeqFreq {
    _seq: String,
    _freq: f32,
}

pub struct Reads<'a> {
    pub reads: Vec<SequenceRecord<'a>>,
    pub unique_seqs: Vec<SeqFreq>,
    pub min_freq: f32,
}

pub trait FreqCalculator {}

pub trait Trim {
    fn trim_fwd_primer<'a>(
        record: &'a mut SequenceRecord<'a>,
        _primer: &'a str,
    ) -> Result<&'a SequenceRecord<'a>> {
        Ok(record)
    }

    fn trim_rev_primer<'a>(
        record: &'a mut SequenceRecord<'a>,
        _primer: &'a str,
    ) -> Result<&'a SequenceRecord<'a>> {
        Ok(record)
    }
}

// impl FreqCalculator for Reads {}

// impl Trim for Reads {}
