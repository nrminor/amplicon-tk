use std::fs::File;

use amplicon_tk::cli::{self, Commands};
use amplicon_tk::io::{PrimerReader, RefReader, SeqReader};
use amplicon_tk::primers::{define_amplicons, ref_to_dict};
use amplicon_tk::reads::{Reads, Trimming};
use clap::Parser;
use color_eyre::eyre::Result;
use flate2::Compression;
use noodles::bed::Reader as BedReader;
use noodles::fasta::Reader as FaReader;
use noodles::fastq::io::Writer as FqWriter;
use noodles::fastq::Reader as FqReader;
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
            ref_file,
            keep_multi: _,
            fwd_suffix,
            rev_suffix,
            freq_min,
            len,
            output,
        }) => {
            let fastq = FqReader::read_fq(input_file)?;
            let bed = BedReader::read_bed(bed_file)?;
            let mut fasta = FaReader::read_ref(ref_file)?;

            let ref_dict = ref_to_dict(&mut fasta)?;
            let scheme = define_amplicons(bed, &ref_dict, fwd_suffix, rev_suffix)?;

            let trimmed_reads = Reads::new(fastq, *freq_min, *len).run_trimming(scheme)?;

            let file = File::create(output)?;
            let encoder = flate2::write::GzEncoder::new(file, Compression::default());
            let mut writer = FqWriter::new(encoder);

            for record in trimmed_reads.reads {
                writer.write_record(&record)?;
            }
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
