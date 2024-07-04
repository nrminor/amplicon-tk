//! `amplicon-tk` implements three subcommands—trim, sort, and consensus—that, unlike many
//! other tools, are amplicon-aware. By that, we mean that amplicon-tk enforces a simple rule:
//! any given read must contain both primers in at least one amplicon. This ensures that all
//! reads in the resulting dataset will correspond to one, complete amplicon, meaning that PCR
//! chimeras and other artifacts will be removed.
//!
//! Crate `amplicon-tk` contains modules for file I/O (`io`), primer-handling (`primers`), read-
//! handling (`reads`), individual record-handling `record`, consensus sequence-calling (`consensus`),
//! the command-line interface (`cli`), and a work-in-progress Python interface.

use amplicon_tk::cli::{self, Commands};
use clap::Parser;
use color_eyre::eyre::Result;
use tracing_subscriber::EnvFilter;

/// .
///
/// # Errors
///
/// This function will return an error if .
#[tokio::main]
async fn main() -> Result<()> {
    setup()?;

    let cli = cli::Cli::parse();
    match &cli.command {
        Some(Commands::Trim {
            input_file: _,
            bed_file: _,
            fasta_ref: _,
            keep_multi: _,
            left_suffix: _,
            right_suffix: _,
            min_freq: _,
            expected_len: _,
            output: _,
        }) => {
            eprintln!("{}\n", cli::INFO);
            eprintln!("Trimming is not yet ready for use, but it will be available soon!")
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

/// .
///
/// # Errors
///
/// This function will return an error if .
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
