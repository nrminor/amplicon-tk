// use clap::Parser;
// use color_eyre::eyre::{eyre, Result};

// use crate::{cli::Commands, prelude::*};
// use std::path::PathBuf;

// #[allow(dead_code)]
// async fn _main() -> Result<()> {
//     let cli = crate::cli::Cli::parse();

//     #[allow(unused_variables)]
//     match &cli.command {
//         Some(Commands::Trim {
//             input_file,
//             bed_file,
//             fasta_ref,
//             keep_multi,
//             left_suffix,
//             right_suffix,
//             min_freq,
//             expected_len,
//             output,
//         }) => {
//             todo!()
//         }
//         _ => eprintln!("ignored"),
//     }

//     // this is starting to come together as a very nice looking
//     // fluent interface
//     let test_fastq = PathBuf::from("test.fastq");

//     // parse the supplied filetype
//     let file_type = SupportedTypes::from_file_name(&test_fastq);

//     // early return if it's unsupported
//     if file_type.is_none() {
//         return Err(eyre!("Unsupported file type supplied in {:?}", file_type));
//     }

//     // use pattern matching to process the reads based on their record type
//     match file_type.unwrap_or(SupportedTypes::FASTQ) {
//         SupportedTypes::FASTQ => {
//             let mut reader = SeqReader::new_fastq(&test_fastq).await?;
//             let records = reader.records();
//             let _results = RecordStream::from_fastq(records)
//                 .await?
//                 .trim_to_amplicons()
//                 .await?
//                 .run_filters()
//                 .await?;
//         }
//         SupportedTypes::FASTQGZ => {
//             let mut reader = SeqReader::new_fastq_gz(&test_fastq).await?;
//             let records = reader.records();
//             let _results = RecordStream::from_fastq(records)
//                 .await?
//                 .trim_to_amplicons()
//                 .await?
//                 .run_filters()
//                 .await?;
//         }
//         _ => eprintln!("BAM and other files are not yet supported."),
//     }

//     Ok(())
// }
