use derive_new::new;

// #![warn(missing_docs)]

#[derive(Debug, new)]
pub struct PrimerPair<'a> {
    // pub amplicon: &'a str,
    pub fwd: &'a str,
    pub rev: &'a str,
}

#[derive(Debug)]
pub struct AmpliconScheme<'a> {
    scheme: Vec<&'a PrimerPair<'a>>,
}

pub trait PrimerGrouping {
    fn define_amplicons() {}

    fn get_primer_seqs() {}
}

// impl PrimerGrouping for AmpliconScheme {}
