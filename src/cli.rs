use std::path::PathBuf;

use clap::{Parser, Subcommand};

pub const INFO: &str = r"

  ____  ___ ___  ____  _      ____   __   ___   ____          ______  __  _
 /    ||   |   ||    \| |    |    | /  ] /   \ |    \        |      ||  |/ ]
|  o  || _   _ ||  o  ) |     |  | /  / |     ||  _  | _____ |      ||  ' /
|     ||  \_/  ||   _/| |___  |  |/  /  |  O  ||  |  ||     ||_|  |_||    \
|  _  ||   |   ||  |  |     | |  /   \_ |     ||  |  ||_____|  |  |  |     \
|  |  ||   |   ||  |  |     | |  \     ||     ||  |  |         |  |  |  .  |
|__|__||___|___||__|  |_____||____\____| \___/ |__|__|         |__|  |__|\_|

ampliton-tk: A Command Line and Python Interface for Amplicon-Aware FASTQ Operations
====================================================================================

amplicon-tk implements three subcommands—trim, sort, and consensus—that, unlike many
other tools, are amplicon-aware. By that, we mean that amplicon-tk enforces a simple rule:
any given read must contain both primers in at least one amplicon. This ensures that all
reads in the resulting dataset will correspond to one, complete amplicon, meaning that PCR
chimeras and other artifacts will be removed.

";

#[derive(Parser)]
#[clap(name = "amplicon-tk")]
#[clap(about = INFO)]
#[clap(version = "v0.1.0")]
pub struct Cli {
    #[command(flatten)]
    pub verbose: clap_verbosity_flag::Verbosity,

    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Subcommand)]
pub enum Commands {
    #[clap(
        about = "Trim a set of reads down to only those reads that contain a complete amplicon.",
        aliases = &["tr", "tirm", "trm", "tri", "tm"])]
    Trim {
        /// Input FASTQ file (optionally compressed with gzip or bgzip)
        #[arg(short, long, required = true)]
        input_file: PathBuf,

        /// Input BED file of primer coordinates
        #[arg(short, long, required = false)]
        bed_file: PathBuf,

        /// Primer sequences in FASTA format
        #[arg(short, long, required = true)]
        primer_file: PathBuf,

        /// Reference sequence in FASTA format
        #[arg(short, long, required = false)]
        ref_file: PathBuf,

        /// Whether to keep reads that contain multiple pairs of primers
        #[arg(short, long, required = false, default_value_t = false)]
        keep_multi: bool,

        /// Output file name
        #[arg(short, long, required = false, default_value = "trimmed.fastq.gz")]
        output: String,
    },

    #[clap(
            about = "Sort reads representing each amplicon into their own FASTQs, one per amplicon.",
            aliases = &["so", "srt", "st", "srot"])]
    Sort {
        /// Input FASTQ file (optionally compressed with gzip or bgzip)
        #[arg(short, long, required = true)]
        input_file: PathBuf,

        /// Input BED file of primer coordinates
        #[arg(short, long, required = false)]
        bed_file: PathBuf,

        /// Primer sequences in FASTA format
        #[arg(short, long, required = true)]
        primer_file: PathBuf,

        /// Reference sequence in FASTA format
        #[arg(short, long, required = false)]
        ref_file: PathBuf,

        /// Minimum frequency for variations of the same amplicon
        #[arg(short, long, required = false, default_value_t = 0.0)]
        min_freq: f32,

        /// Whether to keep reads that contain multiple pairs of primers
        #[arg(short, long, required = false, default_value_t = false)]
        keep_multi: bool,
    },

    #[clap(
            about = "Sort reads representing each amplicon into their own sets, call a consensus sequence for each set, and save it into an output FASTA file.",
            aliases = &["cons", "co", "cd", "consseq", "cseq", "cnsns"])]
    Consensus {
        /// Input FASTQ file (optionally compressed with gzip or bgzip)
        #[arg(short, long, required = true)]
        input_file: PathBuf,

        /// Input BED file of primer coordinates
        #[arg(short, long, required = false)]
        bed_file: PathBuf,

        /// Primer sequences in FASTA format
        #[arg(short, long, required = true)]
        primer_file: PathBuf,

        /// Reference sequence in FASTA format
        #[arg(short, long, required = false)]
        ref_file: PathBuf,

        /// Minimum frequency for variations of the same amplicon
        #[arg(short, long, required = false, default_value_t = 0.0)]
        min_freq: f32,

        /// Whether to keep reads that contain multiple pairs of primers
        #[arg(short, long, required = false, default_value_t = false)]
        keep_multi: bool,

        /// Output file name
        #[arg(short, long, required = false, default_value = "amplicons.fasta")]
        output: String,
    },
}
