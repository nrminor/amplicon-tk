//! `amplicon-tk` implements three subcommands—trim, sort, and consensus—that, unlike many
//! other tools, are amplicon-aware. By that, we mean that amplicon-tk enforces a simple rule:
//! any given read must contain both primers in at least one amplicon. This ensures that all
//! reads in the resulting dataset will correspond to one, complete amplicon, meaning that PCR
//! chimeras and other artifacts will be removed.
//!
//! Crate `amplicon-tk` contains modules for file I/O (`io`), primer-handling (`primers`), read-
//! handling (`reads`), individual record-handling `record`, consensus sequence-calling
//! (`consensus`), the command-line interface (`cli`), and a work-in-progress Python interface.

use std::{fs::File, path::PathBuf, sync::Arc};

#[allow(unused_imports)]
use amplicon_tk::{
    amplicons::{define_amplicons, ref_to_dict},
    cli::{self, Commands},
    index::Index,
    io::{io_selector, Bed, Fasta, InputType, PrimerReader, RefReader, SeqReader},
    reads::Trimming,
};
use amplicon_tk::{
    amplicons::{AmpliconScheme, DefineAmplicons},
    filtering::FilterSettings,
    io::Tsv,
};
use clap::Parser;
use color_eyre::{
    eyre::{eyre, Result},
    owo_colors::OwoColorize,
};
use flate2::bufread::GzDecoder;
use pyo3::exceptions::PyArithmeticError;
use tracing_subscriber::EnvFilter;

#[tokio::main]
async fn main() -> Result<()> {
    // set up the color-eyre display and tracer
    setup()?;

    // parse command line arguments and use a match statement to determine behavior based on the
    // provided subcommand
    let cli = cli::Cli::parse();
    match &cli.command {
        Some(Commands::Index {
            input_file: _,
            bed_file: _,
            fasta_ref: _,
            left_suffix: _,
            right_suffix: _,
        }) => {
            todo!()
        }
        Some(Commands::Extract {
            input_file,
            bed_file,
            fasta_ref,
            primer_fasta,
            primer_table,
            left_suffix,
            right_suffix,
            output,
            demux,
        }) => {
            // initialize the primer scheme given the provided input files
            let _scheme = if let (Some(bed), Some(fasta)) = (bed_file, fasta_ref) {
                let primers = Bed::read_primers(bed)?;
                let mut parsed_ref = Fasta::read_ref(fasta)?;
                let ref_dict = ref_to_dict(&mut parsed_ref).await?;
                primers
                    .define_amplicons(&ref_dict, left_suffix, right_suffix)
                    .await
            } else if let (Some(fasta), Some(tsv)) = (primer_fasta, primer_table) {
                let mut primer_seqs = Fasta::read_primers(fasta)?;
                let _seq_dict = ref_to_dict(&mut primer_seqs).await?;
                let _parsed_tsv = Tsv::read_primers(tsv)?;
                todo!()
            } else {
                Err(eyre!("Either `--bed_file` and `--fasta_ref` or `--primer_fasta` and `--primer_table` must be provided. Please double check that one of those pairs of arguments were specified before trying again."))
            }?;

            // define input and output types for the reads
            let input_type = io_selector(input_file).await?;
            let output_name = format!("{}{}", output, input_type.extension());
            let _output_path = PathBuf::from(output_name);

            // Use pattern-matching to handle the input based on what type it is
            match input_type {
                InputType::FASTQGZ(_) => todo!(),
                InputType::FASTQ(_) => todo!(),
                InputType::BAM(_) => todo!(),
            }
        }
        Some(Commands::Trim {
            input_file,
            bed_file,
            fasta_ref,
            keep_multi: _,
            left_suffix,
            right_suffix,
            min_freq,
            expected_len,
            output,
        }) => {
            // pull in the primers
            let bed = Bed::read_primers(bed_file)?;

            // pull in the reference
            let mut fasta = Fasta::read_ref(fasta_ref)?;

            // convert the reference to a hashmap and use it to pull in the primer pairs for each
            // amplicon
            let ref_dict = ref_to_dict(&mut fasta).await?;
            let scheme = define_amplicons(bed, &ref_dict, left_suffix, right_suffix).await?;

            // hash the current primer scheme to compare with a potential index
            let current_hash = scheme.hash_amplicon_scheme()?;
            let _safe_scheme = Arc::from(scheme);

            // define input and output types for the reads
            let input_type = io_selector(input_file).await?;
            let output_name = format!("{}{}", output, input_type.extension());
            let _output_path = PathBuf::from(output_name);
            // still need to work out how to select different input and output types

            // based on the file type, run lazy, asynchronous trimming with the appropriate record type
            match input_type {
                InputType::FASTQGZ(supported_type) => {
                    // attempt to retrieve a set of unique sequences from an index to use with filtering
                    let unique_seqs = supported_type.load_index(input_file, &current_hash)?;

                    // bundle the requested filter settings. These settings will be None if no unique sequences
                    // could be retrieved from the index
                    let filters = FilterSettings::new(min_freq, expected_len, &unique_seqs);
                    let _safe_filters = Arc::from(filters);

                    // load an appropriate reader
                    let mut _reader = supported_type.read_seq_reads(input_file).await?;

                    // perform trimming based on the supported type
                    todo!()
                    // supported_type
                    //     .trim(input_file, &output_path, safe_scheme, safe_filters)
                    //     .await?
                }
                InputType::FASTQ(supported_type) => {
                    let unique_seqs = supported_type.load_index(input_file, &current_hash)?;
                    let filters = FilterSettings::new(min_freq, expected_len, &unique_seqs);
                    let _safe_filters = Arc::from(filters);
                    todo!();
                    // supported_type
                    //     .trim(input_file, &output_path, safe_scheme, safe_filters)
                    //     .await?
                }
                InputType::BAM(_supported_type) => {
                    eprintln!("Unaligned BAM inputs are not yet supported but will be soon!")
                }
            };
        }
        Some(Commands::Sort {
            input_file: _,
            bed_file: _,
            primer_file: _,
            ref_file: _,
            min_freq: _,
            keep_multi: _,
        }) => {
            eprintln!("{}\n", cli::INFO);
            eprintln!("\nSorting is not yet ready for use, but it will be available soon!")
        }
        Some(Commands::Consensus {
            input_file: _,
            bed_file: _,
            primer_file: _,
            ref_file: _,
            min_freq: _,
            keep_multi: _,
            output: _,
        }) => {
            eprintln!("{}\n", cli::INFO);
            eprintln!("\nAmplicon consensus calling is not yet ready for use, but it will be available soon!")
        }
        None => {
            eprintln!("{}\n", cli::INFO);
        }
    }

    Ok(())
}

fn setup() -> Result<()> {
    // if std::env::var("RUST_LIB_BACKTRACE").is_err() {
    //     std::env::set_var("RUST_LIB_BACKTRACE", "1")
    // }
    // color_eyre::install()?;

    // if std::env::var("RUST_LOG").is_err() {
    //     std::env::set_var("RUST_LOG", "info")
    // }
    tracing_subscriber::fmt::fmt()
        .with_env_filter(EnvFilter::from_default_env())
        .init();

    Ok(())
}
