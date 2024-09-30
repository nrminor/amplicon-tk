use std::path::PathBuf;

use bincode::Options;
use clap::{Parser, Subcommand};

pub const INFO: &str = r"

  ____  ___ ___  ____  _      ____   __   ___   ____          ______  __  _
 /    ||   |   ||    \| |    |    | /  ] /   \ |    \        |      ||  |/ ]
|  o  || _   _ ||  o  ) |     |  | /  / |     ||  _  | _____ |      ||  ' /
|     ||  \_/  ||   _/| |___  |  |/  /  |  O  ||  |  ||     ||_|  |_||    \
|  _  ||   |   ||  |  |     | |  /   \_ |     ||  |  ||_____|  |  |  |     \
|  |  ||   |   ||  |  |     | |  \     ||     ||  |  |         |  |  |  .  |
|__|__||___|___||__|  |_____||____\____| \___/ |__|__|         |__|  |__|\_|

amplicon-tk: A Command Line and Python Interface for Amplicon-Aware FASTQ Operations
====================================================================================

amplicon-tk implements a series of helpful subcommands that, unlike many other tools, are
amplicon-aware. By that, we mean that amplicon-tk enforces a simple rule: any given read
must contain both primers in at least one amplicon. This ensures that all reads in the
resulting dataset will correspond to one, complete amplicon, meaning that PCR chimeras
and other artifacts will be removed.

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
        about = "Index FASTQ records and record unique amplicon statistics. Indexing implicitly \
        finds and trims primers before identifying unique amplicon sequences.",
        aliases = &["id", "ind", "idx", "ix"])]
    Index {
        /// Input FASTQ file (optionally compressed with gzip or bgzip)
        #[arg(short, long, required = true)]
        input_file: PathBuf,

        /// Input BED file of primer coordinates
        #[arg(short, long, required = true)]
        bed_file: PathBuf,

        /// Reference sequence in FASTA format
        #[arg(short, long, required = true)]
        fasta_ref: PathBuf,

        /// The suffix used to identify forward primers in the provided BED file
        #[arg(short, long, required = false, default_value = "_LEFT")]
        left_suffix: String,

        /// The suffix used to identify reverse primers in the provided BED file
        #[arg(short, long, required = false, default_value = "_RIGHT")]
        right_suffix: String,
    },

    #[clap(
        about = "Extract only those reads that represent complete amplicons, optionally demultiplexing \
        into separate files. Use this subcommand when you want to filter or sort to amplicons without \
        trimming off ends containing primers/barcodes.",
        aliases = &["x", "xtr", "extra", "demux", "demultiplex", "d"]
    )]
    Extract {
        /// Input FASTQ, gzipped FASTQ, or BAM file to query for amplicons.
        #[arg(short, long, required = true)]
        input_file: PathBuf,

        /// Input BED file of primer coordinates
        #[arg(short, long, required = false)]
        bed_file: Option<PathBuf>,

        /// Reference sequence in FASTA format. Required if a primer bed was provided.
        #[arg(short, long, required = false)]
        fasta_ref: Option<PathBuf>,

        /// FASTA file of primer sequences and labels.
        #[arg(short = 'f', long, required = false)]
        primer_fasta: Option<PathBuf>,

        /// Tab-delimited text file pairing up primers. Required if a primer FASTA was provided instead of a primer BED.
        #[arg(short = 't', long, required = false)]
        primer_table: Option<PathBuf>,

        /// The suffix used to identify forward primers in the provided BED file
        #[arg(short, long, required = false, default_value = "_LEFT")]
        left_suffix: String,

        /// The suffix used to identify reverse primers in the provided BED file
        #[arg(short, long, required = false, default_value = "_RIGHT")]
        right_suffix: String,

        /// Output file name
        #[arg(short, long, required = false, default_value = "extracted_amplicons")]
        output: String,

        /// Whether to demultiplex each amplicon/barcode pair into its own output file.
        #[arg(short, long, required = false, default_value_t = false)]
        demux: bool,
    },

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

        /// Reference sequence in FASTA format. Required if a primer bed was provided.
        #[arg(short, long, required = false)]
        fasta_ref: PathBuf,

        // /// FASTA file of primer sequences and labels.
        // #[arg(short = "f", long, required = false)]
        // primer_fasta: PathBuf,

        // /// Tab-delimited text file pairing up primers. Required if a primer FASTA was provided instead of a primer BED.
        // #[arg(short = "t", long, required = false)]
        // primer_table: PathBuf,
        /// Whether to keep reads that contain multiple pairs of primers
        #[arg(short, long, required = false, default_value_t = false)]
        keep_multi: bool,

        /// The suffix used to identify forward primers in the provided BED file
        #[arg(short, long, required = false, default_value = "_LEFT")]
        left_suffix: String,

        /// The suffix used to identify reverse primers in the provided BED file
        #[arg(short, long, required = false, default_value = "_RIGHT")]
        right_suffix: String,

        /// The minimum allowed frequency for amplicon variants
        #[arg(short, long, required = false)]
        min_freq: Option<f64>,

        /// Whether to filter by an expected maximum length for amplicons in this scheme
        #[arg(short, long, required = false)]
        expected_len: Option<usize>,

        /// Output file name
        #[arg(short, long, required = false, default_value = "trimmed")]
        output: String,
    },

    #[clap(
            about = "Trim and sort reads representing each amplicon into their own FASTQs, one per amplicon. \
            Indexing with `amplicon-tk index` must be performed before sorting.",
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
            about = "Trim and sort reads representing each amplicon into their own sets, call a consensus \
            sequence for each set, and save it into an output FASTA file. Indexing with `amplicon-tk index` \
            must be performed before calling consensus amplicons.",
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
