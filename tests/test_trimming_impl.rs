// use color_eyre::eyre::Result;

// #[test]
// fn test_with_primer_match() -> Result<()> {
//     let seq_str: &str;
//     let qual_str: &str;
//     let expected_seq: &str;
//     let expected_qual: &str;
//     let fwd: &str;
//     let rev: &str;
//     let _ = match (&seq_str.find(&fwd), &seq_str.find(&rev)) {
//         (Some(fwd_idx), Some(rev_idx)) => {
//             let new_start = fwd_idx + fwd.len();
//             let new_end = rev_idx - rev.len();
//             let trimmed_seq = &seq_str[new_start..new_end];
//             let trimmed_qual = &qual_str[new_start..new_end];

//             assert_eq!(trimmed_seq, expected_seq);
//             assert_eq!(trimmed_qual, expected_qual);
//             ()
//         }
//         _ => (),
//     };
//     Ok(())
// }

// #[test]
// fn test_with_primer_mismatch() -> Result<()> {
//     let seq_str: &str;
//     let qual_str: &str;
//     let expected_seq: &str;
//     let expected_qual: &str;
//     let fwd: &str;
//     let rev: &str;
//     let _ = match (&seq_str.find(&fwd), &seq_str.find(&rev)) {
//         (Some(fwd_idx), Some(rev_idx)) => {
//             let new_start = fwd_idx + fwd.len();
//             let new_end = rev_idx - rev.len();
//             let trimmed_seq = &seq_str[new_start..new_end];
//             let trimmed_qual = &qual_str[new_start..new_end];

//             assert_eq!(trimmed_seq, expected_seq);
//             assert_eq!(trimmed_qual, expected_qual);
//             ()
//         }
//         _ => (),
//     };
//     Ok(())
// }
