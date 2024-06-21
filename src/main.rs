use amplicon_tk::cli::{self, Commands};
use clap::Parser;
use color_eyre::eyre::Result;
use tracing_subscriber::EnvFilter;

// #![warn(missing_docs)]

#[tokio::main]
async fn main() -> Result<()> {
    setup()?;

    let cli = cli::Cli::parse();
    match &cli.command {
        Some(Commands::Trim {
            input_file,
            bed_file,
            primer_file,
            ref_file,
            keep_multi,
            output,
        }) => {
            eprintln!("{}\n", cli::INFO);
        }
        Some(Commands::Sort {
            input_file,
            bed_file,
            primer_file,
            ref_file,
            min_freq,
            keep_multi,
        }) => {
            eprintln!("{}\n", cli::INFO);
        }
        Some(Commands::Consensus {
            input_file,
            bed_file,
            primer_file,
            ref_file,
            min_freq,
            keep_multi,
            output,
        }) => {
            eprintln!("{}\n", cli::INFO);
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
