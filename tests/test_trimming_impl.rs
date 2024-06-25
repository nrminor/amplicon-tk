use color_eyre::eyre::Result;

#[test]
fn test_with_primer_match() -> Result<()> {
    let seq_str: &str =
        "TGTTTCCACTGGAGGATACTCACCCCTCTTGCACTCAAGTTAAACAGTTTCCAAAGCGTACTATGGTTAAGCCACAGCCT";
    let qual_str: &str =
        "445656:11DHHGJPSHFDCDDOMIBD@?@DDD><<<<FFLDFGIJCIKJIKFGSOSCC=;98782-,-..112299:B=";
    assert_eq!(seq_str.len(), qual_str.len());

    let expected_seq: &str = "ACTCACCCCTCTTGCACTCAAGTTAAACAGTTTCCAAAGCG";
    let expected_qual: &str = "FDCDDOMIBD@?@DDD><<<<FFLDFGIJCIKJIKFGSOSC";
    let fwd: &str = "TGGAGGAT";
    let rev: &str = "TACTATGG";
    if let (Some(fwd_idx), Some(rev_idx)) = (&seq_str.find(fwd), &seq_str.find(rev)) {
        let new_start = fwd_idx + fwd.len();
        let new_end = rev_idx;
        eprintln!("Forward end index: {new_start}");
        eprintln!("Reverse start index: {new_end}");
        let trimmed_seq = &seq_str[new_start..*new_end];
        eprintln!("Trimmed sequence: {trimmed_seq}");
        let trimmed_qual = &qual_str[new_start..*new_end];
        eprintln!("Trimmed quality scores: {trimmed_qual}");

        assert_eq!(
            trimmed_seq.len(),
            trimmed_qual.len(),
            "Somehow the sequence and quality score are not the same length."
        );
        assert_eq!(
            trimmed_seq,
            expected_seq,
            "The trimmed sequence looks different than the expected sequence:\nTrimmed: {trimmed_seq}\nExpected: {expected_seq}"
        );
        assert_eq!(
            trimmed_qual,
            expected_qual,
            "The trimmed sequence looks different than the expected sequence:\nTrimmed: {trimmed_qual}\nExpected: {expected_qual}"
        );
        assert!(
            !trimmed_seq.contains(fwd),
            "The trimmed sequence still contains the forward primer."
        );
        assert!(
            !trimmed_seq.contains(fwd),
            "The trimmed sequence still contains the forward primer."
        );
    }
    Ok(())
}

#[test]
fn test_with_no_primer_match() -> Result<()> {
    let seq_str: &str =
        "TGTTTCCACTGGAGGATACTCACCCCTCTTGCACTCAAGTTAAACAGTTTCCAAAGCGTACTATGGTTAAGCCACAGCCT";
    let qual_str: &str =
        "445656:11DHHGJPSHFDCDDOMIBD@?@DDD><<<<FFLDFGIJCIKJIKFGSOSCC=;98782-,-..112299:B=";
    assert_eq!(seq_str.len(), qual_str.len());

    let expected_seq: &str =
        "TGTTTCCACTGGAGGATACTCACCCCTCTTGCACTCAAGTTAAACAGTTTCCAAAGCGTACTATGGTTAAGCCACAGCCT";
    let expected_qual: &str =
        "445656:11DHHGJPSHFDCDDOMIBD@?@DDD><<<<FFLDFGIJCIKJIKFGSOSCC=;98782-,-..112299:B=";
    let fwd: &str = "TGGAGGAT";
    let rev: &str = "TACTATGG";
    let test_result =
        if let (Some(fwd_idx), Some(rev_idx)) = (&seq_str.find(fwd), &seq_str.find(rev)) {
            let new_start = fwd_idx + fwd.len();
            let new_end = rev_idx;
            eprintln!("Forward end index: {new_start}");
            eprintln!("Reverse start index: {new_end}");
            let trimmed_seq = &seq_str[new_start..*new_end];
            eprintln!("Trimmed sequence: {trimmed_seq}");
            let trimmed_qual = &qual_str[new_start..*new_end];
            eprintln!("Trimmed quality scores: {trimmed_qual}");

            Some((trimmed_seq, trimmed_qual))
        } else {
            None
        };

    match test_result {
        Some((trimmed_seq, trimmed_qual)) => {
            assert_eq!(
                trimmed_seq.len(),
                trimmed_qual.len(),
                "Somehow the sequence and quality score are not the same length."
            );
            assert_eq!(
                    trimmed_seq,
                    expected_seq,
                    "The trimmed sequence looks different than the expected sequence:\nTrimmed: {trimmed_seq}\nExpected: {expected_seq}"
                );
            assert_eq!(
                    trimmed_qual,
                    expected_qual,
                    "The trimmed sequence looks different than the expected sequence:\nTrimmed: {trimmed_qual}\nExpected: {expected_qual}"
                );
            assert!(
                !trimmed_seq.contains(fwd),
                "The trimmed sequence still contains the forward primer."
            );
            assert!(
                !trimmed_seq.contains(fwd),
                "The trimmed sequence still contains the forward primer."
            );
        }
        None => eprintln!("No matches found"),
    }

    Ok(())
}
