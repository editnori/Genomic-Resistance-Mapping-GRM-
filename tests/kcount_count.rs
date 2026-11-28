//! Integration tests for `grm::kcount::superkmer` and `grm kcount count`:
//! - slow reference counter parity across `k = 5..=64` on plain/gzip/multi-member
//!   gzip FASTA and FASTQ, CRLF, lowercase, N/IUPAC runs, palindromes, homopolymer
//!   runs > 255 bp, short/empty contigs, and empty files;
//! - exact `--min-qual 20` boundary testing (`Q19`/`Q20`/`Q21`);
//! - exact QC metrics (`gc_fraction`, `total_bases`, `valid_bases`,
//!   `usable_read_len_median`, `singleton_fraction`, `n50`) and `peak_coverage`
//!   histogram parity against a slow reference;
//! - `u32::MAX` saturation refusal and `--allow-saturation` clamping;
//! - spill vs in-memory SHA-256 byte identity and 1-thread vs 8-thread SHA-256
//!   byte identity;
//! - `.kset` v1 parity against `grm count` at `k = 5, 21, 31, 32`;
//! - `.prof.bin` per-position count parity, anonymized `r000001` IDs vs
//!   `--keep-record-ids`, and `inspect --verify`;
//! - end-to-end CLI checks (`sample` default name, `InputModality`, `--input-type`
//!   mismatch errors, `--minimizer-len` validation, extension checks, `SIGINT`/
//!   `SIGTERM` spill cleanup, and dead-PID stale spill folder cleanup).

use flate2::Compression;
use flate2::write::GzEncoder;
use grm::dna::KmerSpace;
use grm::kcount::format::{
    DEFAULT_SAMPLE_NAME, GrKsetReader, InputModality, KtabReader, convert_to_v1, hex_bytes,
    read_profiles_verified,
};
use grm::kcount::key::{KmerSpace64, KmerSpace128, MAX_K_U64};
use grm::kcount::superkmer::{
    CountConfig, InputTypeChoice, count_inputs, finalize_u32_count, merge_counted_runs,
};
use grm::kset::{self, KsetReader};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use ring::digest::{SHA256, digest};
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

fn grm_bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_grm"))
}

fn sha256_file(path: &Path) -> String {
    let bytes = fs::read(path).unwrap();
    hex_bytes(digest(&SHA256, &bytes).as_ref())
}

fn gzip_bytes(data: &[u8]) -> Vec<u8> {
    let mut enc = GzEncoder::new(Vec::new(), Compression::fast());
    enc.write_all(data).unwrap();
    enc.finish().unwrap()
}

fn rc_dna(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => unreachable!(),
        })
        .collect()
}

fn canonical_kmer_str(s: &str) -> String {
    let rc = rc_dna(s);
    if s <= rc.as_str() { s.to_string() } else { rc }
}

struct SlowRefResult {
    counts: BTreeMap<String, u32>,
    total_reads: u64,
    total_bases: u64,
    valid_bases: u64,
    gc_bases: u64,
    masked_bases: u64,
    raw_lens: Vec<u64>,
    usable_lens: Vec<u64>,
    per_seq_profiles: Vec<(String, u64, Vec<u32>)>,
}

/// Pure slow reference parser and counter directly over ASCII strings.
fn slow_count_records(
    records: &[(String, String, Option<String>)],
    k: usize,
    min_qual: u8,
    min_count: u32,
) -> SlowRefResult {
    let mut counts: BTreeMap<String, u32> = BTreeMap::new();
    let qual_floor = 33u8.saturating_add(min_qual);
    let mut total_bases = 0u64;
    let mut valid_bases = 0u64;
    let mut gc_bases = 0u64;
    let mut masked_bases = 0u64;
    let mut raw_lens = Vec::new();
    let mut usable_lens = Vec::new();

    for (_, seq, qual) in records {
        let sb = seq.as_bytes();
        let qb = qual.as_deref().map(str::as_bytes);
        total_bases += sb.len() as u64;
        raw_lens.push(sb.len() as u64);

        let mut clean: Vec<Option<u8>> = Vec::with_capacity(sb.len());
        let mut rec_valid = 0u64;
        for (i, &b) in sb.iter().enumerate() {
            let up = b.to_ascii_uppercase();
            let low_q = min_qual > 0 && qb.is_some_and(|q| q[i] < qual_floor);
            if !low_q && matches!(up, b'A' | b'C' | b'G' | b'T') {
                rec_valid += 1;
                valid_bases += 1;
                if matches!(up, b'C' | b'G') {
                    gc_bases += 1;
                }
                clean.push(Some(up));
            } else {
                masked_bases += 1;
                clean.push(None);
            }
        }
        usable_lens.push(rec_valid);

        if clean.len() >= k {
            for w in clean.windows(k) {
                if w.iter().all(Option::is_some) {
                    let s: String = w.iter().map(|c| c.unwrap() as char).collect();
                    let can = canonical_kmer_str(&s);
                    *counts.entry(can).or_insert(0) += 1;
                }
            }
        }
    }

    let kept: BTreeMap<String, u32> = counts
        .iter()
        .filter(|&(_, &c)| c >= min_count)
        .map(|(s, &c)| (s.clone(), c))
        .collect();

    let mut per_seq_profiles = Vec::new();
    for (name, seq, qual) in records {
        let sb = seq.as_bytes();
        let qb = qual.as_deref().map(str::as_bytes);
        let mut clean: Vec<Option<u8>> = Vec::with_capacity(sb.len());
        for (i, &b) in sb.iter().enumerate() {
            let up = b.to_ascii_uppercase();
            let low_q = min_qual > 0 && qb.is_some_and(|q| q[i] < qual_floor);
            if !low_q && matches!(up, b'A' | b'C' | b'G' | b'T') {
                clean.push(Some(up));
            } else {
                clean.push(None);
            }
        }
        let mut pos_counts = Vec::new();
        if clean.len() >= k {
            for w in clean.windows(k) {
                if w.iter().all(Option::is_some) {
                    let s: String = w.iter().map(|c| c.unwrap() as char).collect();
                    let can = canonical_kmer_str(&s);
                    pos_counts.push(*kept.get(&can).unwrap_or(&0));
                } else {
                    pos_counts.push(0);
                }
            }
        }
        per_seq_profiles.push((name.clone(), sb.len() as u64, pos_counts));
    }

    SlowRefResult {
        counts,
        total_reads: records.len() as u64,
        total_bases,
        valid_bases,
        gc_bases,
        masked_bases,
        raw_lens,
        usable_lens,
        per_seq_profiles,
    }
}

fn decode_ktab_to_map(path: &Path, k: usize) -> BTreeMap<String, u32> {
    let reader = KtabReader::open_verified(path).unwrap();
    assert_eq!(reader.k, k);
    let mut map = BTreeMap::new();
    if k <= MAX_K_U64 {
        let space = KmerSpace64::new(k);
        for (mixed, c) in reader.read_all::<u64>().unwrap() {
            let dec = space.decode(space.unmix(mixed));
            let prev = map.insert(dec, c);
            assert!(prev.is_none());
        }
    } else {
        let space = KmerSpace128::new(k);
        for (mixed, c) in reader.read_all::<u128>().unwrap() {
            let dec = space.decode(space.unmix(mixed));
            let prev = map.insert(dec, c);
            assert!(prev.is_none());
        }
    }
    map
}

fn random_dna(rng: &mut ChaCha8Rng, len: usize) -> String {
    (0..len)
        .map(|_| b"ACGT"[rng.random_range(0..4)] as char)
        .collect()
}

// ---------------------------------------------------------------------------
// 1. Slow reference counter parity across k = 5..=64 and tricky formats
// ---------------------------------------------------------------------------

#[test]
fn slow_reference_parity_across_k_5_to_64_and_edge_inputs() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(202_601);

    let half32 = random_dna(&mut rng, 32);
    let pal64 = format!("{half32}{}", rc_dna(&half32));
    let half16 = random_dna(&mut rng, 16);
    let pal32 = format!("{half16}{}", rc_dna(&half16));
    let repeat_unit = random_dna(&mut rng, 80);

    let fasta_records: Vec<(String, String, Option<String>)> = vec![
        (
            "c_mix".into(),
            format!(
                "acgtACGTnnNN{}RYSWKM{}",
                random_dna(&mut rng, 140),
                random_dna(&mut rng, 90)
            ),
            None,
        ),
        ("c_pal64".into(), format!("{pal64}{pal64}{pal32}"), None),
        ("c_hom255".into(), "A".repeat(255), None),
        ("c_hom256".into(), "C".repeat(256), None),
        ("c_hom600".into(), "G".repeat(600), None),
        ("c_hom1000".into(), "T".repeat(1000), None),
        ("c_rep".into(), repeat_unit.repeat(7), None),
        ("c_short".into(), "ACG".into(), None),
        ("c_empty".into(), String::new(), None),
    ];

    let mut fasta_text = String::new();
    for (idx, (name, seq, _)) in fasta_records.iter().enumerate() {
        fasta_text.push('>');
        fasta_text.push_str(name);
        fasta_text.push_str(" extra desc\r\n");
        for chunk in seq.as_bytes().chunks(60) {
            fasta_text.push_str(std::str::from_utf8(chunk).unwrap());
            if idx + 1 < fasta_records.len() {
                fasta_text.push_str("\r\n");
            } else {
                fasta_text.push('\n');
            }
        }
    }
    if fasta_text.ends_with('\n') {
        fasta_text.pop();
    }
    let fasta_plain = dir.path().join("edge.fa");
    fs::write(&fasta_plain, fasta_text.as_bytes()).unwrap();

    let mid = fasta_text.len() / 2;
    let fasta_mgz = dir.path().join("edge_multi.fa.gz");
    let multi_gz = [
        gzip_bytes(&fasta_text.as_bytes()[..mid]),
        gzip_bytes(&fasta_text.as_bytes()[mid..]),
    ]
    .concat();
    fs::write(&fasta_mgz, &multi_gz).unwrap();

    let mut fq_records: Vec<(String, String, Option<String>)> = Vec::new();
    let mut fq_text = String::new();
    for i in 0..80 {
        let seq = if i % 10 == 0 {
            format!(
                "{}NNN{}",
                random_dna(&mut rng, 40),
                random_dna(&mut rng, 40)
            )
        } else if i == 1 {
            "ACG".to_string()
        } else {
            random_dna(&mut rng, 110)
        };
        let mut qual: String = (0..seq.len())
            .map(|j| {
                if j == 0 && i % 2 == 0 {
                    '@'
                } else if j == 0 {
                    '+'
                } else {
                    (33u8 + rng.random_range(10u8..=40u8)) as char
                }
            })
            .collect();
        if seq.is_empty() {
            qual.clear();
        }
        fq_text.push_str(&format!("@read_{i} lane=1\r\n{seq}\r\n+\r\n{qual}\r\n"));
        fq_records.push((format!("read_{i}"), seq, Some(qual)));
    }
    let fq_gz = dir.path().join("reads.fq.gz");
    fs::write(&fq_gz, gzip_bytes(fq_text.as_bytes())).unwrap();

    let empty_path = dir.path().join("empty.fa");
    fs::write(&empty_path, b"").unwrap();
    let empty_ktab = dir.path().join("empty.ktab");
    let empty_cfg = CountConfig {
        k: 31,
        out_path: empty_ktab.clone(),
        ..Default::default()
    };
    let empty_sum = count_inputs(&[empty_path], &empty_cfg).unwrap();
    assert_eq!(empty_sum.distinct_retained, 0);
    assert_eq!(empty_sum.modality, InputModality::Unknown);
    assert_eq!(empty_sum.histogram.peak_coverage, 0);
    let empty_reader = KtabReader::open_verified(&empty_ktab).unwrap();
    assert_eq!(empty_reader.n_kmers, 0);
    assert_eq!(empty_reader.modality, InputModality::Unknown);

    for &k in &[5usize, 15, 21, 31, 32, 33, 45, 61, 64] {
        let slow_fa = slow_count_records(&fasta_records, k, 0, 1);
        let ktab_plain = dir.path().join(format!("fa_{k}_plain.ktab"));
        let ktab_mgz = dir.path().join(format!("fa_{k}_mgz.ktab"));
        let grk_mgz = dir.path().join(format!("fa_{k}_mgz.grk"));

        let s_plain = count_inputs(
            std::slice::from_ref(&fasta_plain),
            &CountConfig {
                k,
                out_path: ktab_plain.clone(),
                ..Default::default()
            },
        )
        .unwrap();
        let s_mgz = count_inputs(
            std::slice::from_ref(&fasta_mgz),
            &CountConfig {
                k,
                out_path: ktab_mgz.clone(),
                gr_kset_out: Some(grk_mgz.clone()),
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(s_plain.modality, InputModality::Assembly);
        assert_eq!(s_mgz.modality, InputModality::Assembly);
        assert_eq!(sha256_file(&ktab_plain), sha256_file(&ktab_mgz));
        assert_eq!(
            decode_ktab_to_map(&ktab_plain, k),
            slow_fa.counts,
            "FASTA k={k}"
        );

        let gr = GrKsetReader::open_verified(&grk_mgz).unwrap();
        assert_eq!(gr.modality, InputModality::Assembly);
        assert_eq!(gr.n_kmers, slow_fa.counts.len() as u64);

        let slow_fq = slow_count_records(&fq_records, k, 20, 1);
        let ktab_fq = dir.path().join(format!("fq_{k}.ktab"));
        let s_fq = count_inputs(
            std::slice::from_ref(&fq_gz),
            &CountConfig {
                k,
                min_qual: 20,
                out_path: ktab_fq.clone(),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(s_fq.modality, InputModality::Reads);
        assert_eq!(
            decode_ktab_to_map(&ktab_fq, k),
            slow_fq.counts,
            "FASTQ k={k}"
        );
    }
}

// ---------------------------------------------------------------------------
// 2. Exact --min-qual 20 boundary test (Q19 masked, Q20 & Q21 kept!)
// ---------------------------------------------------------------------------

#[test]
fn min_qual_20_exact_boundary_q19_masked_q20_q21_kept() {
    let dir = tempfile::tempdir().unwrap();
    let k = 11usize;
    let s1 = "ACGTACGTACGTA";
    let s2 = "TGCATGCATGCAT";
    let s3 = "GGGGCCCCAAAAT";

    let q19 = (33u8 + 19) as char;
    let q20 = (33u8 + 20) as char;
    let q21 = (33u8 + 21) as char;

    let mut qual_s1 = vec![q21; s1.len()];
    qual_s1[6] = q19;
    let qual_s1_str: String = qual_s1.into_iter().collect();
    let qual_s2_str: String = std::iter::repeat_n(q20, s2.len()).collect();
    let qual_s3_str: String = std::iter::repeat_n(q21, s3.len()).collect();

    let records = vec![
        (
            "q19_mid".to_string(),
            s1.to_string(),
            Some(qual_s1_str.clone()),
        ),
        (
            "q20_all".to_string(),
            s2.to_string(),
            Some(qual_s2_str.clone()),
        ),
        (
            "q21_all".to_string(),
            s3.to_string(),
            Some(qual_s3_str.clone()),
        ),
    ];
    let fq_path = dir.path().join("qual_edge.fq");
    let mut text = String::new();
    for (id, s, q) in &records {
        text.push_str(&format!("@{id}\n{s}\n+\n{}\n", q.as_ref().unwrap()));
    }
    fs::write(&fq_path, text).unwrap();

    let slow = slow_count_records(&records, k, 20, 1);
    let ktab_path = dir.path().join("qual_edge.ktab");
    let prof_path = dir.path().join("qual_edge.prof.bin");
    let sum = count_inputs(
        &[fq_path],
        &CountConfig {
            k,
            min_qual: 20,
            out_path: ktab_path.clone(),
            profile_out: Some(prof_path.clone()),
            ..Default::default()
        },
    )
    .unwrap();

    assert_eq!(sum.qc.masked_bases, 1);
    assert_eq!(sum.qc.valid_bases, 38);
    assert_eq!(sum.qc.total_bases, 39);
    assert_eq!(decode_ktab_to_map(&ktab_path, k), slow.counts);

    let (_, profs) = read_profiles_verified(&prof_path).unwrap();
    assert_eq!(profs.len(), 3);
    assert_eq!(profs[0].counts, vec![0, 0, 0]);
    assert!(profs[1].counts.iter().all(|&c| c > 0));
    assert!(profs[2].counts.iter().all(|&c| c > 0));
    for (got, (_, want_len, want_counts)) in profs.iter().zip(&slow.per_seq_profiles) {
        assert_eq!(got.seq_len, *want_len);
        assert_eq!(got.counts, *want_counts);
    }
}

// ---------------------------------------------------------------------------
// 3. Exact QC metrics, peak_coverage histogram, and split long FASTA contig
// ---------------------------------------------------------------------------

#[test]
fn qc_metrics_and_peak_coverage_match_slow_reference_exactly() {
    let dir = tempfile::tempdir().unwrap();
    let k = 15usize;

    let mut rng = ChaCha8Rng::seed_from_u64(4242);
    let mut records: Vec<(String, String, Option<String>)> = Vec::new();
    let peak_seq = "GGGGCCCCGGCCGGCCGGCC".to_string();
    for i in 0..5 {
        records.push((format!("peak_{i}"), peak_seq.clone(), Some("I".repeat(20))));
    }
    let minor_seq = "AAAATTTTAAAACCCC".to_string();
    for i in 0..3 {
        records.push((
            format!("minor_{i}"),
            minor_seq.clone(),
            Some("I".repeat(16)),
        ));
    }
    let sing_a = format!("NN{}NN", random_dna(&mut rng, 20));
    let sing_b = format!("NN{}NN", random_dna(&mut rng, 20));
    let sing_c = random_dna(&mut rng, 22);
    records.push(("s_a".into(), sing_a.clone(), Some("I".repeat(sing_a.len()))));
    records.push(("s_b".into(), sing_b.clone(), Some("I".repeat(sing_b.len()))));
    records.push(("s_c".into(), sing_c.clone(), Some("I".repeat(sing_c.len()))));

    let fq_path = dir.path().join("qc_exact.fq");
    let mut text = String::new();
    for (id, s, q) in &records {
        text.push_str(&format!("@{id}\n{s}\n+\n{}\n", q.as_ref().unwrap()));
    }
    fs::write(&fq_path, text).unwrap();

    let slow = slow_count_records(&records, k, 0, 1);
    let ktab_path = dir.path().join("qc_exact.ktab");
    let hist_path = dir.path().join("qc_exact.hist.json");
    let sum = count_inputs(
        &[fq_path],
        &CountConfig {
            k,
            out_path: ktab_path,
            hist_out: Some(hist_path.clone()),
            ..Default::default()
        },
    )
    .unwrap();

    assert_eq!(sum.qc.total_reads, slow.total_reads);
    assert_eq!(sum.qc.total_bases, slow.total_bases);
    assert_eq!(sum.qc.valid_bases, slow.valid_bases);
    assert_eq!(sum.qc.masked_bases, slow.masked_bases);
    let want_gc = slow.gc_bases as f64 / slow.valid_bases as f64;
    assert!(
        (sum.qc.gc_fraction - want_gc).abs() < 1e-12,
        "gc got {} want {want_gc}",
        sum.qc.gc_fraction
    );
    let mut usable_sorted = slow.usable_lens.clone();
    usable_sorted.sort_unstable();
    let n = usable_sorted.len();
    let want_median = if n % 2 == 1 {
        usable_sorted[n / 2] as f64
    } else {
        (usable_sorted[n / 2 - 1] + usable_sorted[n / 2]) as f64 / 2.0
    };
    assert_eq!(sum.qc.usable_read_len_median, want_median);
    assert!(
        (sum.qc.usable_read_len_mean - (slow.valid_bases as f64 / slow.total_reads as f64)).abs()
            < 1e-12
    );
    assert_eq!(
        sum.qc.usable_read_len_p10,
        usable_sorted[((n - 1) as f64 * 0.1).round() as usize] as f64
    );
    assert_eq!(
        sum.qc.usable_read_len_p90,
        usable_sorted[((n - 1) as f64 * 0.9).round() as usize] as f64
    );

    let mut raw_desc = slow.raw_lens.clone();
    raw_desc.sort_unstable_by(|a, b| b.cmp(a));
    let mut run = 0u64;
    let mut want_n50 = 0u64;
    for &l in &raw_desc {
        run += l;
        if run * 2 >= slow.total_bases {
            want_n50 = l;
            break;
        }
    }
    assert_eq!(sum.qc.n50, want_n50);

    let want_f1 = slow.counts.values().filter(|&&c| c == 1).count() as u64;
    assert_eq!(sum.histogram.f1, want_f1);
    assert_eq!(sum.histogram.distinct_total, slow.counts.len() as u64);
    assert!(
        (sum.qc.singleton_fraction - (want_f1 as f64 / slow.counts.len() as f64)).abs() < 1e-12
    );
    assert!(sum.histogram.f1 > 6);
    assert_eq!(sum.histogram.peak_coverage, 5);
    assert!(hist_path.exists());

    // Peak coverage tie-break test: equal frequency (2 distinct k-mers each)
    // at c=4 and c=9 -> smaller count (4) wins.
    let seq_c4 = "ACGTTTGCAGTCGATG".to_string();
    let seq_c9 = "TGCAAACCGGTAAGTC".to_string();
    let mut tie_fa = String::new();
    for i in 0..4 {
        tie_fa.push_str(&format!(">a{i}\n{seq_c4}\n"));
    }
    for i in 0..9 {
        tie_fa.push_str(&format!(">b{i}\n{seq_c9}\n"));
    }
    let tie_path = dir.path().join("tie.fa");
    fs::write(&tie_path, tie_fa).unwrap();
    let tie_sum = count_inputs(
        &[tie_path],
        &CountConfig {
            k,
            out_path: dir.path().join("tie.ktab"),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(tie_sum.histogram.peak_coverage, 4);

    // Singletons-only peak_coverage == 1.
    let s_only_path = dir.path().join("singletons.fa");
    fs::write(&s_only_path, format!(">s\n{}\n", random_dna(&mut rng, 40))).unwrap();
    let s_only_sum = count_inputs(
        &[s_only_path],
        &CountConfig {
            k,
            out_path: dir.path().join("singletons.ktab"),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(s_only_sum.histogram.peak_coverage, 1);

    // Long single FASTA contig exceeding the 256 KiB chunk split threshold at
    // `--max-memory-mb 1` verifies overlap continuity and exact QC totals.
    let long_seq = format!(
        "{}NNNN{}",
        random_dna(&mut rng, 350_000),
        random_dna(&mut rng, 350_000)
    );
    let long_fa = dir.path().join("long_contig.fa");
    fs::write(&long_fa, format!(">long1\n{long_seq}\n")).unwrap();
    let long_ktab_spill = dir.path().join("long_spill.ktab");
    let long_ktab_mem = dir.path().join("long_mem.ktab");
    let s_long_spill = count_inputs(
        std::slice::from_ref(&long_fa),
        &CountConfig {
            k: 31,
            input_type: InputTypeChoice::Assembly,
            max_memory_mb: 1,
            out_path: long_ktab_spill.clone(),
            ..Default::default()
        },
    )
    .unwrap();
    let s_long_mem = count_inputs(
        std::slice::from_ref(&long_fa),
        &CountConfig {
            k: 31,
            input_type: InputTypeChoice::Assembly,
            max_memory_mb: 512,
            out_path: long_ktab_mem.clone(),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(s_long_spill.qc.total_reads, 1);
    assert_eq!(s_long_spill.qc.total_bases, 700_004);
    assert_eq!(s_long_spill.qc.valid_bases, 700_000);
    assert_eq!(s_long_spill.qc.masked_bases, 4);
    assert_eq!(s_long_spill.qc.n50, 700_004);
    assert_eq!(s_long_spill.qc, s_long_mem.qc);
    assert_eq!(s_long_spill.histogram, s_long_mem.histogram);
    assert_eq!(sha256_file(&long_ktab_spill), sha256_file(&long_ktab_mem));
}

// ---------------------------------------------------------------------------
// 4. u32::MAX saturation refusal & --allow-saturation clamping (Finding 5)
// ---------------------------------------------------------------------------

#[test]
fn u32_max_saturation_refuses_by_default_and_clamps_when_enabled() {
    assert_eq!(
        finalize_u32_count(4_294_967_295, false).unwrap(),
        (u32::MAX, false)
    );
    let err = finalize_u32_count(4_294_967_296, false).unwrap_err();
    assert!(format!("{err:#}").contains("u32::MAX"));
    assert!(format!("{err:#}").contains("--allow-saturation"));
    assert_eq!(
        finalize_u32_count(5_000_000_000, true).unwrap(),
        (u32::MAX, true)
    );

    let runs_u64: Vec<Vec<(u64, u64)>> = vec![
        vec![(11, 3_000_000_000), (22, 100)],
        vec![(11, 2_000_000_000), (22, 50)],
    ];
    assert!(merge_counted_runs(&runs_u64, false).is_err());
    let (merged64, sat64) = merge_counted_runs(&runs_u64, true).unwrap();
    assert_eq!(sat64, 1);
    assert_eq!(merged64, vec![(11, u32::MAX), (22, 150)]);

    let runs_u128: Vec<Vec<(u128, u64)>> = vec![
        vec![(1u128 << 90, 3_000_000_000)],
        vec![(1u128 << 90, 2_000_000_000)],
    ];
    assert!(merge_counted_runs(&runs_u128, false).is_err());
    let (merged128, sat128) = merge_counted_runs(&runs_u128, true).unwrap();
    assert_eq!(sat128, 1);
    assert_eq!(merged128, vec![(1u128 << 90, u32::MAX)]);
}

// ---------------------------------------------------------------------------
// 5. Spill vs in-memory AND 1-thread vs 8-thread SHA-256 byte identity
// ---------------------------------------------------------------------------

#[test]
fn spill_vs_in_memory_and_1_vs_8_threads_are_byte_for_byte_identical() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(9900);

    let base_chr = random_dna(&mut rng, 150_000);
    let fq_path = dir.path().join("noisy.fq");
    {
        let mut f = File::create(&fq_path).unwrap();
        let bytes = base_chr.as_bytes();
        let qual = "I".repeat(150);
        for i in 0..15_000 {
            let pos = rng.random_range(0..(bytes.len() - 150));
            let mut read = bytes[pos..pos + 150].to_vec();
            if i % 3 == 0 {
                let err_pos = rng.random_range(0..150);
                read[err_pos] = b"ACGT"[(read[err_pos] as usize + 1) % 4];
            }
            writeln!(f, "@r{i}").unwrap();
            f.write_all(&read).unwrap();
            writeln!(f, "\n+\n{qual}").unwrap();
        }
    }

    let spill_tmp = dir.path().join("spill_tmp");
    let mem_ktab = dir.path().join("mem.ktab");
    let mem_grk = dir.path().join("mem.grk");
    let mem_prof = dir.path().join("mem.prof.bin");
    let mem_hist = dir.path().join("mem.hist.json");

    let spill_ktab = dir.path().join("spill.ktab");
    let spill_grk = dir.path().join("spill.grk");
    let spill_prof = dir.path().join("spill.prof.bin");
    let spill_hist = dir.path().join("spill.hist.json");

    let s_mem = count_inputs(
        std::slice::from_ref(&fq_path),
        &CountConfig {
            k: 31,
            max_memory_mb: 512,
            out_path: mem_ktab.clone(),
            gr_kset_out: Some(mem_grk.clone()),
            profile_out: Some(mem_prof.clone()),
            hist_out: Some(mem_hist.clone()),
            ..Default::default()
        },
    )
    .unwrap();

    let s_spill = count_inputs(
        std::slice::from_ref(&fq_path),
        &CountConfig {
            k: 31,
            max_memory_mb: 1,
            tmp_dir: Some(spill_tmp.clone()),
            out_path: spill_ktab.clone(),
            gr_kset_out: Some(spill_grk.clone()),
            profile_out: Some(spill_prof.clone()),
            hist_out: Some(spill_hist.clone()),
            ..Default::default()
        },
    )
    .unwrap();

    assert!(!s_mem.spilled_to_disk);
    assert!(s_spill.spilled_to_disk);
    assert_eq!(s_mem.fnv64, s_spill.fnv64);
    assert_eq!(s_mem.keys_sha256, s_spill.keys_sha256);
    assert_eq!(s_mem.counts_sha256, s_spill.counts_sha256);
    assert_eq!(s_mem.sketch_hex, s_spill.sketch_hex);
    assert_eq!(s_mem.histogram, s_spill.histogram);
    assert_eq!(s_mem.qc, s_spill.qc);
    assert_eq!(sha256_file(&mem_ktab), sha256_file(&spill_ktab));
    assert_eq!(sha256_file(&mem_grk), sha256_file(&spill_grk));
    assert_eq!(sha256_file(&mem_prof), sha256_file(&spill_prof));
    assert_eq!(sha256_file(&mem_hist), sha256_file(&spill_hist));
    let leftovers: Vec<_> = fs::read_dir(&spill_tmp).unwrap().collect();
    assert!(
        leftovers.is_empty(),
        "leftover spill entries: {leftovers:?}"
    );

    for (label, max_mem) in [("inmem", "512"), ("spill", "1")] {
        let mut hashes = Vec::new();
        for threads in ["1", "8"] {
            let ktab = dir.path().join(format!("t{threads}_{label}.ktab"));
            let grk = dir.path().join(format!("t{threads}_{label}.grk"));
            let prof = dir.path().join(format!("t{threads}_{label}.prof.bin"));
            let hist = dir.path().join(format!("t{threads}_{label}.hist.json"));
            let meta = dir.path().join(format!("t{threads}_{label}.meta.json"));
            let out = grm_bin()
                .env("RAYON_NUM_THREADS", threads)
                .args([
                    "kcount",
                    "count",
                    "--input",
                    fq_path.to_str().unwrap(),
                    "--k",
                    "31",
                    "--max-memory-mb",
                    max_mem,
                    "--out",
                    ktab.to_str().unwrap(),
                    "--gr-kset-out",
                    grk.to_str().unwrap(),
                    "--profile-out",
                    prof.to_str().unwrap(),
                    "--hist-out",
                    hist.to_str().unwrap(),
                    "--meta-out",
                    meta.to_str().unwrap(),
                ])
                .output()
                .unwrap();
            assert!(
                out.status.success(),
                "threads={threads} {label} failed: {}",
                String::from_utf8_lossy(&out.stderr)
            );
            hashes.push((
                sha256_file(&ktab),
                sha256_file(&grk),
                sha256_file(&prof),
                sha256_file(&hist),
                sha256_file(&meta),
            ));
        }
        assert_eq!(
            hashes[0], hashes[1],
            "1-thread vs 8-thread differed ({label})"
        );
    }
}

// ---------------------------------------------------------------------------
// 6. .kset v1 parity against grm count at k = 5, 21, 31, 32
// ---------------------------------------------------------------------------

#[test]
fn kset_v1_parity_against_grm_count_at_k_5_21_31_32() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(5213132);

    let fa_text = format!(
        ">ctg1\n{}NNN{}rys{}\n>ctg2\n{}\n>short\nACGT\n",
        random_dna(&mut rng, 15_000),
        random_dna(&mut rng, 5_000).to_ascii_lowercase(),
        random_dna(&mut rng, 10_000),
        random_dna(&mut rng, 12_000),
    );
    let fa_path = dir.path().join("asm.fna");
    fs::write(&fa_path, &fa_text).unwrap();

    for &k in &[5usize, 21, 31, 32] {
        let space = KmerSpace::new(k);
        let (v1_keys, _) = kset::count_file(&fa_path, &space, 1).unwrap();
        let v1_path = dir.path().join(format!("ref_k{k}.kset"));
        kset::write(&v1_path, "sample", &space, &v1_keys).unwrap();
        let v1_reader = KsetReader::open(&v1_path).unwrap();

        let ktab_path = dir.path().join(format!("kc_k{k}.ktab"));
        let grk_path = dir.path().join(format!("kc_k{k}.grk"));
        let sum = count_inputs(
            std::slice::from_ref(&fa_path),
            &CountConfig {
                k,
                sample_name: "sample".to_string(),
                out_path: ktab_path.clone(),
                gr_kset_out: Some(grk_path.clone()),
                ..Default::default()
            },
        )
        .unwrap();

        assert_eq!(sum.distinct_retained, v1_reader.n_kmers, "n_kmers k={k}");
        assert_eq!(sum.fnv64, kset::digest(&v1_keys), "fnv64 k={k}");

        let want_sketch: Vec<String> = v1_reader
            .sketch(1000)
            .unwrap()
            .into_iter()
            .map(|v| format!("{v:016x}"))
            .collect();
        assert_eq!(sum.sketch_hex, want_sketch, "sketch k={k}");

        let back_v1 = dir.path().join(format!("back_k{k}.kset"));
        convert_to_v1(&ktab_path, &back_v1, "sample").unwrap();
        assert_eq!(
            sha256_file(&v1_path),
            sha256_file(&back_v1),
            "convert-v1 bytes differed at k={k}"
        );
    }

    let root = dir.path().join("grm_root");
    fs::create_dir_all(root.join("fasta")).unwrap();
    fs::create_dir_all(root.join("species/syn")).unwrap();
    fs::write(root.join("fasta/g1.fna.gz"), gzip_bytes(fa_text.as_bytes())).unwrap();
    fs::write(root.join("species/syn/genomes.txt"), "g1\n").unwrap();

    let count_out = grm_bin()
        .args([
            "--data",
            root.to_str().unwrap(),
            "count",
            "--species",
            "syn",
            "--k",
            "31",
            "--allow-missing-catalog",
        ])
        .output()
        .unwrap();
    assert!(
        count_out.status.success(),
        "grm count failed: {}",
        String::from_utf8_lossy(&count_out.stderr)
    );
    let grm_kset = root.join("species/syn/ksets/g1.kset");
    assert!(grm_kset.exists(), "grm count did not create g1.kset");
    let kc_ktab = dir.path().join("cli_g1.ktab");
    let kc_out = grm_bin()
        .args([
            "kcount",
            "count",
            "--input",
            root.join("fasta/g1.fna.gz").to_str().unwrap(),
            "--k",
            "31",
            "--name",
            "g1",
            "--out",
            kc_ktab.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(kc_out.status.success());

    let i1: serde_json::Value = serde_json::from_slice(
        &grm_bin()
            .args([
                "kcount",
                "inspect",
                "--input",
                grm_kset.to_str().unwrap(),
                "--verify",
            ])
            .output()
            .unwrap()
            .stdout,
    )
    .unwrap();
    let i2: serde_json::Value = serde_json::from_slice(
        &grm_bin()
            .args([
                "kcount",
                "inspect",
                "--input",
                kc_ktab.to_str().unwrap(),
                "--verify",
            ])
            .output()
            .unwrap()
            .stdout,
    )
    .unwrap();
    assert_eq!(i1["n_kmers"], i2["n_kmers"]);
    assert_eq!(i1["fnv64"], i2["fnv64"]);
    assert_eq!(i1["keys_sha256"], i2["keys_sha256"]);
    assert_eq!(i2["verified"], true);
}

// ---------------------------------------------------------------------------
// 7. .prof.bin per-position parity, anonymized IDs vs --keep-record-ids,
//    and inspect --verify
// ---------------------------------------------------------------------------

#[test]
fn profile_bin_matches_slow_reference_and_anonymizes_ids_by_default() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(777);

    let shared = random_dna(&mut rng, 80);
    let records: Vec<(String, String, Option<String>)> = vec![
        (
            "SECRET_READ_ALPHA".into(),
            format!("{shared}NNNN{}", random_dna(&mut rng, 70)),
            Some("I".repeat(154)),
        ),
        (
            "SECRET_READ_BETA".into(),
            shared.clone(),
            Some("I".repeat(80)),
        ),
        ("SECRET_SHORT".into(), "ACGT".into(), Some("IIII".into())),
        ("SECRET_EMPTY".into(), String::new(), Some(String::new())),
    ];
    let fq_path = dir.path().join("prof_in.fq");
    let mut fq_text = String::new();
    for (id, s, q) in &records {
        fq_text.push_str(&format!("@{id}\n{s}\n+\n{}\n", q.as_ref().unwrap()));
    }
    fs::write(&fq_path, fq_text).unwrap();

    for &k in &[21usize, 31, 45, 64] {
        let slow = slow_count_records(&records, k, 0, 2);
        let ktab = dir.path().join(format!("p_{k}.ktab"));
        let prof_anon = dir.path().join(format!("p_{k}_anon.prof.bin"));
        let prof_keep = dir.path().join(format!("p_{k}_keep.prof.bin"));

        count_inputs(
            std::slice::from_ref(&fq_path),
            &CountConfig {
                k,
                min_count: 2,
                keep_record_ids: false,
                out_path: ktab.clone(),
                profile_out: Some(prof_anon.clone()),
                ..Default::default()
            },
        )
        .unwrap();

        let (hdr, got_anon) = read_profiles_verified(&prof_anon).unwrap();
        assert_eq!(hdr.k, k);
        assert_eq!(got_anon.len(), records.len());
        for (i, (got, (_, want_len, want_counts))) in
            got_anon.iter().zip(&slow.per_seq_profiles).enumerate()
        {
            assert_eq!(got.name, format!("r{:06}", i + 1));
            assert_eq!(got.seq_len, *want_len);
            assert_eq!(got.counts, *want_counts, "k={k} record {i}");
        }
        let raw_anon = fs::read(&prof_anon).unwrap();
        assert!(!raw_anon.windows(11).any(|w| w == b"SECRET_READ"));

        count_inputs(
            std::slice::from_ref(&fq_path),
            &CountConfig {
                k,
                min_count: 2,
                keep_record_ids: true,
                out_path: ktab.clone(),
                profile_out: Some(prof_keep.clone()),
                ..Default::default()
            },
        )
        .unwrap();
        let (_, got_keep) = read_profiles_verified(&prof_keep).unwrap();
        for (got, (want_name, want_len, want_counts)) in got_keep.iter().zip(&slow.per_seq_profiles)
        {
            assert_eq!(got.name, *want_name);
            assert_eq!(got.seq_len, *want_len);
            assert_eq!(got.counts, *want_counts);
        }

        let insp_out = grm_bin()
            .args([
                "kcount",
                "inspect",
                "--input",
                prof_anon.to_str().unwrap(),
                "--verify",
            ])
            .output()
            .unwrap();
        assert!(insp_out.status.success());
        let j: serde_json::Value = serde_json::from_slice(&insp_out.stdout).unwrap();
        assert_eq!(j["verified"], true);
        assert_eq!(j["k"], k);
        assert_eq!(j["n_sequences"], records.len());
        assert_eq!(j["total_bases"], slow.total_bases);
    }

    let grk_only = dir.path().join("pres.grk");
    let grk_prof = dir.path().join("pres.prof.bin");
    count_inputs(
        &[fq_path],
        &CountConfig {
            k: 31,
            min_count: 2,
            presence_only_out: true,
            out_path: grk_only,
            profile_out: Some(grk_prof.clone()),
            ..Default::default()
        },
    )
    .unwrap();
    let (_, gp) = read_profiles_verified(&grk_prof).unwrap();
    assert!(gp[1].counts.iter().all(|&c| c == 1));
}

// ---------------------------------------------------------------------------
// 8. End-to-end CLI contracts (Findings 7, 9, 13, 24, 25)
// ---------------------------------------------------------------------------

#[test]
fn cli_count_contracts_no_filename_leak_modality_minimizer_and_extensions() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(8080);

    let secret_fq = dir.path().join("DO_NOT_LEAK_STEM_999.fq");
    let seq = random_dna(&mut rng, 200);
    fs::write(
        &secret_fq,
        format!("@SECRET_SEQ_HEADER\n{seq}\n+\n{}\n", "I".repeat(200)),
    )
    .unwrap();

    let secret_fa = dir.path().join("ANOTHER_PRIVATE_STEM.fna");
    fs::write(&secret_fa, format!(">CONTIG_PRIVATE\n{seq}\n")).unwrap();

    let ktab = dir.path().join("out.ktab");
    let grk = dir.path().join("out.grk");
    let hist = dir.path().join("out.hist.json");
    let meta = dir.path().join("out.meta.json");
    let prof = dir.path().join("out.prof.bin");

    let run = grm_bin()
        .args([
            "kcount",
            "count",
            "--input",
            secret_fq.to_str().unwrap(),
            "--k",
            "31",
            "--out",
            ktab.to_str().unwrap(),
            "--gr-kset-out",
            grk.to_str().unwrap(),
            "--hist-out",
            hist.to_str().unwrap(),
            "--meta-out",
            meta.to_str().unwrap(),
            "--profile-out",
            prof.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "cli count failed: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let stdout_str = String::from_utf8(run.stdout).unwrap();
    let j: serde_json::Value = serde_json::from_str(&stdout_str).unwrap();
    assert_eq!(j["sample"], DEFAULT_SAMPLE_NAME);
    assert_eq!(j["k"], 31);
    assert_eq!(j["modality"], "reads");
    assert_eq!(j["saturated_counts"], 0);
    assert_eq!(j["distinct_retained"], 170);
    assert!(!stdout_str.contains("DO_NOT_LEAK_STEM_999"));
    assert!(!stdout_str.contains("SECRET_SEQ_HEADER"));

    for p in [&ktab, &grk, &hist, &meta, &prof] {
        let raw = fs::read(p).unwrap();
        assert!(
            !raw.windows(16).any(|w| w == b"DO_NOT_LEAK_STEM"),
            "filename leaked into {}",
            p.display()
        );
        assert!(
            !raw.windows(17).any(|w| w == b"SECRET_SEQ_HEADER"),
            "record id leaked into {}",
            p.display()
        );
    }

    for p in [&ktab, &grk] {
        let insp = grm_bin()
            .args([
                "kcount",
                "inspect",
                "--input",
                p.to_str().unwrap(),
                "--verify",
            ])
            .output()
            .unwrap();
        assert!(insp.status.success());
        let ij: serde_json::Value = serde_json::from_slice(&insp.stdout).unwrap();
        assert_eq!(ij["verified"], true);
        assert_eq!(ij["sample"], "sample");
        assert_eq!(ij["modality"], "reads");
    }

    let fa_ktab = dir.path().join("fa.ktab");
    let fa_run = grm_bin()
        .args([
            "kcount",
            "count",
            "--input",
            secret_fa.to_str().unwrap(),
            "--k",
            "31",
            "--out",
            fa_ktab.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(fa_run.status.success());
    let fa_j: serde_json::Value = serde_json::from_slice(&fa_run.stdout).unwrap();
    assert_eq!(fa_j["modality"], "assembly");
    assert_eq!(
        KtabReader::open_verified(&fa_ktab).unwrap().modality,
        InputModality::Assembly
    );

    let bad_reads_on_fa = grm_bin()
        .args([
            "kcount",
            "count",
            "--input",
            secret_fa.to_str().unwrap(),
            "--input-type",
            "reads",
            "--out",
            dir.path().join("bad1.ktab").to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!bad_reads_on_fa.status.success());
    assert!(
        String::from_utf8_lossy(&bad_reads_on_fa.stderr)
            .contains("conflicts with --input-type reads")
    );

    let bad_asm_on_fq = grm_bin()
        .args([
            "kcount",
            "count",
            "--input",
            secret_fq.to_str().unwrap(),
            "--input-type",
            "assembly",
            "--out",
            dir.path().join("bad2.ktab").to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!bad_asm_on_fq.status.success());
    assert!(
        String::from_utf8_lossy(&bad_asm_on_fq.stderr)
            .contains("conflicts with --input-type assembly")
    );

    let help = grm_bin()
        .args(["kcount", "count", "--help"])
        .output()
        .unwrap();
    assert!(help.status.success());
    let help_str = String::from_utf8(help.stdout).unwrap();
    assert!(
        help_str.contains("min(k, 31)"),
        "help missing min(k, 31): {help_str}"
    );
    assert!(!help_str.to_ascii_lowercase().contains("uncapped"));
    for k_str in ["32", "33", "64"] {
        let bad_m = grm_bin()
            .args([
                "kcount",
                "count",
                "--input",
                secret_fa.to_str().unwrap(),
                "--k",
                k_str,
                "--minimizer-len",
                "32",
                "--out",
                dir.path().join("badm.ktab").to_str().unwrap(),
            ])
            .output()
            .unwrap();
        assert!(!bad_m.status.success());
        assert!(String::from_utf8_lossy(&bad_m.stderr).contains("min(k, 31)"));
    }

    let bad_ext = grm_bin()
        .args([
            "kcount",
            "count",
            "--input",
            secret_fa.to_str().unwrap(),
            "--presence-only",
            "--out",
            dir.path().join("wrong.ktab").to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!bad_ext.status.success());
    assert!(String::from_utf8_lossy(&bad_ext.stderr).contains("expected .grk"));

    let ok_pres = dir.path().join("ok_pres.grk");
    let ok_pres_run = grm_bin()
        .args([
            "kcount",
            "count",
            "--input",
            secret_fa.to_str().unwrap(),
            "--presence-only",
            "--out",
            ok_pres.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(ok_pres_run.status.success());
    assert_eq!(GrKsetReader::open_verified(&ok_pres).unwrap().n_kmers, 170);
}

// ---------------------------------------------------------------------------
// 9. Stale dead-PID spill cleanup & SIGTERM/SIGINT cleanup (Finding 12)
// ---------------------------------------------------------------------------

#[cfg(unix)]
#[test]
fn stale_dead_pid_spill_cleaned_and_signals_remove_active_spill_dirs() {
    let dir = tempfile::tempdir().unwrap();
    let spill_root = dir.path().join("spills");
    fs::create_dir_all(&spill_root).unwrap();

    let mut short = Command::new("true").spawn().unwrap();
    let dead_pid = short.id();
    let _ = short.wait().unwrap();

    let dead_dir = spill_root.join(format!("grm-kcount-{dead_pid}-0"));
    fs::create_dir_all(&dead_dir).unwrap();
    fs::write(dead_dir.join("p0000.skm"), b"stale").unwrap();

    let live_pid = std::process::id();
    let live_dir = spill_root.join(format!("grm-kcount-{live_pid}-999999"));
    fs::create_dir_all(&live_dir).unwrap();
    fs::write(live_dir.join("keep.skm"), b"live").unwrap();

    let tiny_fa = dir.path().join("tiny.fa");
    fs::write(&tiny_fa, ">c1\nACGTACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
    let tiny_out = dir.path().join("tiny.ktab");

    let run = grm_bin()
        .args([
            "kcount",
            "count",
            "--input",
            tiny_fa.to_str().unwrap(),
            "--k",
            "21",
            "--tmp-dir",
            spill_root.to_str().unwrap(),
            "--out",
            tiny_out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(run.status.success());
    assert!(
        !dead_dir.exists(),
        "dead-PID spill directory should have been removed"
    );
    assert!(
        live_dir.exists(),
        "live-PID spill directory must NOT be removed"
    );
    let _ = fs::remove_dir_all(&live_dir);

    let mut rng = ChaCha8Rng::seed_from_u64(1212);
    let big_fq = dir.path().join("sig_input.fq");
    {
        let mut f = File::create(&big_fq).unwrap();
        let qual = "I".repeat(150);
        let pool: Vec<String> = (0..64).map(|_| random_dna(&mut rng, 150)).collect();
        for i in 0..120_000 {
            let s = &pool[i & 63];
            writeln!(f, "@r{i}\n{s}\n+\n{qual}").unwrap();
        }
    }

    for (sig_name, sig_num) in [("SIGTERM", libc::SIGTERM), ("SIGINT", libc::SIGINT)] {
        let sig_spills = dir.path().join(format!("sig_spills_{sig_name}"));
        fs::create_dir_all(&sig_spills).unwrap();
        let sig_out = dir.path().join(format!("sig_{sig_name}.ktab"));

        let mut child = grm_bin()
            .env("RAYON_NUM_THREADS", "1")
            .args([
                "kcount",
                "count",
                "--input",
                big_fq.to_str().unwrap(),
                "--k",
                "31",
                "--max-memory-mb",
                "1",
                "--tmp-dir",
                sig_spills.to_str().unwrap(),
                "--out",
                sig_out.to_str().unwrap(),
            ])
            .stdout(Stdio::null())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        let cpid = child.id();
        let prefix = format!("grm-kcount-{cpid}-");

        let deadline = Instant::now() + Duration::from_secs(15);
        let mut saw_spill = false;
        while Instant::now() < deadline {
            let found = fs::read_dir(&sig_spills)
                .unwrap()
                .flatten()
                .any(|e| e.file_name().to_string_lossy().starts_with(&prefix));
            if found {
                saw_spill = true;
                // SAFETY: Deliver test termination signal to our spawned child process.
                unsafe {
                    libc::kill(cpid as libc::pid_t, sig_num);
                }
                break;
            }
            if child.try_wait().unwrap().is_some() {
                break;
            }
            std::thread::sleep(Duration::from_millis(5));
        }

        let output = child.wait_with_output().unwrap();
        assert!(
            saw_spill,
            "{sig_name}: spill directory never appeared before child exited"
        );
        assert!(
            !output.status.success(),
            "{sig_name}: expected non-zero exit after signal"
        );
        let leftovers: Vec<PathBuf> = fs::read_dir(&sig_spills)
            .unwrap()
            .flatten()
            .map(|e| e.path())
            .collect();
        assert!(
            leftovers.is_empty(),
            "{sig_name}: spill directory leaked after signal: {leftovers:?}"
        );
        assert!(
            !sig_out.exists(),
            "{sig_name}: partial output table left behind"
        );
        assert!(
            !PathBuf::from(format!("{}.part", sig_out.display())).exists(),
            "{sig_name}: .part file left behind"
        );
    }
}
