use derive_new::new;

// #![warn(missing_docs)]

#[derive(Debug, new)]
pub struct PrimerPair<'a> {
    pub fwd: &'a str,
    pub rev: &'a str,
}
