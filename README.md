# `amplicon-tk`: A Command Line and Python Interface for Amplicon-aware FASTQ Operations

[![Rust CI](https://github.com/nrminor/amplicon-tk/actions/workflows/ci.yml/badge.svg)](https://github.com/nrminor/amplicon-tk/actions/workflows/ci.yml)

`amplicon-tk` is currently mostly a concept (the subcommand trim is almost finished) but aims to implement three subcommands—trim, sort, and consensus—that, unlike many
other tools, are amplicon-aware. By that, we mean that `amplicon-tk` enforces a simple rule:
any given read must contain both primers in at least one amplicon. This ensures that all
reads in the resulting dataset will correspond to one, complete amplicon, meaning that PCR
chimeras and other artifacts will be removed.

> [!WARNING] > `amplicon-tk` is not ready to be used outside its developer team. An official release of multi-platform binaries alongside publication on [PyPI](https://pypi.org/) and [crates.io](https://crates.io/) will mark when the toolkit is ready for widespread usage. We may also choose to make it available on Bioconda, so stay tuned! See the feature roadmap below for our progress implementing the toolkit.

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

### Feature roadmap:

-   [ ] method for pulling primer sequences based on coordinates in a BED file and sequence in a reference FASTA
-   [x] doc-tested and "black-box" unit-tested read sequence and quality score trimming
-   [ ] trimming and other read operations that are generic to FASTQs and FASTAs
-   [x] method for finding unique sequences alongside frequencies of their instances in the dataset
-   [x] method for sorting out reads into individual "stacks" that can be written out to individual FASTQs
-   [x] readers and writers for all of the above file formats
-   [ ] a method for using primer-pair hits to write one FASTQ per amplicon
-   [ ] a method for calling frequency-aware consensus sequences for each variation of an amplicon and writing those "amplitypes" to their own sequences in a FASTA
-   [ ] async throughout
-   [ ] Python library interface via Pyo3
