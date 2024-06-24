# `amplicon-tk`: A Command Line and Python Interface for Amplicon-aware FASTQ Operations

`amplicon-tk` implements three subcommands—trim, sort, and consensus—that, unlike many
other tools, are amplicon-aware. By that, we mean that amplicon-tk enforces a simple rule:
any given read must contain both primers in at least one amplicon. This ensures that all
reads in the resulting dataset will correspond to one, complete amplicon, meaning that PCR
chimeras and other artifacts will be removed.

```
Usage: amplicon-tk [OPTIONS] [COMMAND]

Commands:
  trim       Trim a set of reads down to only those reads that contain a complete amplicon.
  sort       Sort reads representing each amplicon into their own FASTQs, one per amplicon.
  consensus  Sort reads representing each amplicon into their own sets, call a consensus sequence for each set, and save it into an output FASTA file.
  help       Print this message or the help of the given subcommand(s)

Options:
  -v, --verbose...  Increase logging verbosity
  -q, --quiet...    Decrease logging verbosity
  -h, --help        Print help
  -V, --version     Print version
```

Feature roadmap:

-   method for pulling primer sequences based on coordinates in a BED file and sequence in a reference FASTA
-   doc-tested and "black-box" unit-tested read sequence and quality score trimming
-   trimming and other read operations that are generic to FASTQs and FASTAs
-   method for finding unique sequences alongside frequencies of their instances in the dataset ✅
-   method for sorting out reads into individual "stacks" that can be written out to individual FASTQs ✅
-   readers and writers for all of the above file formats
-   async throughout
-   Python library interface via Pyo3
