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
    cli::{self, Commands},
    index::Index,
    io::{io_selector, Bed, Fasta, InputType, PrimerReader, RefReader, SeqReader},
    primers::{define_amplicons, ref_to_dict},
    reads::{FilterSettings, Trimming},
};
use clap::Parser;
use color_eyre::eyre::Result;
use flate2::bufread::GzDecoder;
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
            input_file,
            bed_file,
            fasta_ref,
            left_suffix,
            right_suffix,
        }) => {
            // defining input and output types for the reads
            let input_type = io_selector(input_file).await?;

            // pulling in the primers
            let primer_type = Bed;
            let bed = primer_type.read_primers(bed_file)?;

            // pulling in the reference
            let ref_type = Fasta;
            let mut fasta = ref_type.read_ref(fasta_ref)?;

            // convert the reference to a hashmap and use it to pull in the primer pairs for each
            // amplicon
            let ref_dict = ref_to_dict(&mut fasta).await?;
            let scheme = define_amplicons(bed, &ref_dict, left_suffix, right_suffix).await?;

            // based on the input filetype, open, decode, and parse the sequence read records
            // lazily and use them to create an index
            match input_type {
                InputType::FASTQGZ(supported_type) => {
                    let opened_file = File::open(input_file)?;
                    let buffer_raw = std::io::BufReader::new(opened_file);
                    let decoded = GzDecoder::new(buffer_raw);
                    let decoded_buffer = std::io::BufReader::new(decoded);
                    let reader = noodles::fastq::Reader::new(decoded_buffer);
                    supported_type.index(reader, scheme, input_file).await?;
                }
                InputType::FASTQ(supported_type) => {
                    let opened_file = File::open(input_file)?;
                    let buffer = std::io::BufReader::new(opened_file);
                    let reader = noodles::fastq::Reader::new(buffer);
                    supported_type.index(reader, scheme, input_file).await?;
                }
                InputType::BAM(_supported_type) => {
                    eprintln!("Unaligned BAM inputs are not yet supported but will be soon!")
                }
            };
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
            let primer_type = Bed;
            let bed = primer_type.read_primers(bed_file)?;

            // pull in the reference
            let ref_type = Fasta;
            let mut fasta = ref_type.read_ref(fasta_ref)?;

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
                    let mut _reader = supported_type.read_reads(input_file).await?;

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
    if std::env::var("RUST_LIB_BACKTRACE").is_err() {
        std::env::set_var("RUST_LIB_BACKTRACE", "1")
    }
    color_eyre::install()?;

    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info")
    }
    tracing_subscriber::fmt::fmt()
        .with_env_filter(EnvFilter::from_default_env())
        .init();

    Ok(())
}
