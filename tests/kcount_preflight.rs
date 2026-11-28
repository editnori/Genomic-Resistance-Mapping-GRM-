//! Exact value, `--min-qual` cutoff, HyperLogLog small/large cardinality,
//! modality-calibrated disk footprint, SHA-256 seed provenance, and end-to-end
//! CLI tests for `grm::kcount::preflight` and `grm kcount preflight`.

use flate2::Compression;
use flate2::write::GzEncoder;
use grm::dna::KmerSpace;
use grm::kcount::format::{DEFAULT_SAMPLE_NAME, InputModality, write_gr_kset};
use grm::kcount::key::KmerSpace128;
use grm::kcount::preflight::{
    BASE_SEED, PreflightOptions, preflight_inputs, preflight_inputs_with_options,
};
use grm::kcount::spss::archive_from_fasta;
use grm::kset;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use ring::digest::{SHA256, digest};
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::Write;
use std::process::Command;

fn rand_dna(rng: &mut ChaCha8Rng, len: usize) -> Vec<u8> {
    (0..len).map(|_| b"ACGT"[rng.random_range(0..4)]).collect()
}

fn grm_bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_grm"))
}

// ---------------------------------------------------------------------------
// 1. Verified SHA-256 derivation of all 4 odd ntHash base seeds (Finding 20)
// ---------------------------------------------------------------------------

#[test]
fn preflight_nthash_seeds_are_odd_and_match_sha256_v1() {
    for (idx, &seed) in BASE_SEED.iter().enumerate() {
        let label = format!("GRM-KCOUNT-NTHASH-SEED-V1:{idx}");
        let h = digest(&SHA256, label.as_bytes());
        let expected = u64::from_le_bytes(h.as_ref()[..8].try_into().unwrap()) | 1;
        assert_eq!(seed, expected, "seed {idx} diverged from {label}");
        assert_eq!(seed & 1, 1, "seed {idx} must be odd");
    }
}

// ---------------------------------------------------------------------------
// 2. Exact `--min-qual` boundary at Q19 / Q20 / Q21 + exact GC fraction
//    (kills preflight.rs:121 `<` vs `<=` and preflight.rs:306 `gc_fraction`!)
// ---------------------------------------------------------------------------

#[test]
fn preflight_min_qual_exact_boundary_q19_q20_q21_and_gc_fraction() {
    let dir = tempfile::tempdir().unwrap();
    let k = 11usize;

    // Build 3 deterministic non-overlapping 20-bp reads with asymmetric base
    // composition (r0 has 4 A, 8 C, 5 G, 3 T -> 13 GC; r1 has 4 A, 6 C, 6 G,
    // 4 T -> 12 GC; r2 has 4 A, 6 C, 6 G, 4 T -> 12 GC with 'G' at index 10).
    // Read 0: all Q20 (Phred+33 = 53 = b'5') -> ALL 10 k-mers kept!
    // Read 1: all Q21 (Phred+33 = 54 = b'6') -> ALL 10 k-mers kept!
    // Read 2: base 10 ('G') set to Q19 (Phred+33 = 52 = b'4'), rest Q20 ->
    //         base 10 is MASKED, splitting into [0..10] (len 10 < 11) and
    //         [11..20] (len 9 < 11) -> 0 k-mers emitted, and 1 G removed from GC!
    let r0 = b"ACGTCCGGACCCGGTTCCAA"; // 20 bp: 4 A, 8 C, 5 G, 3 T -> 13 GC
    let r1 = b"CCGGTTAACCGGCCAATTGG"; // 20 bp: 4 A, 6 C, 6 G, 4 T -> 12 GC
    let r2 = b"TTGGCCAACCGGAATTCCGG"; // 20 bp: 4 A, 6 C, 6 G, 4 T -> idx 10 is 'G'
    assert_eq!(r0.len(), 20);
    assert_eq!(r1.len(), 20);
    assert_eq!(r2.len(), 20);
    assert_eq!(r2[10], b'G');
    assert_eq!(r0.iter().filter(|&&b| b == b'C' || b == b'G').count(), 13);
    assert_eq!(r1.iter().filter(|&&b| b == b'C' || b == b'G').count(), 12);
    assert_eq!(r2.iter().filter(|&&b| b == b'C' || b == b'G').count(), 12);

    let q20 = "5".repeat(20); // Phred 20 (ASCII 53)
    let q21 = "6".repeat(20); // Phred 21 (ASCII 54)
    let mut q_mask = q20.clone().into_bytes();
    q_mask[10] = b'4'; // Phred 19 (ASCII 52)
    let q_mask_str = String::from_utf8(q_mask).unwrap();

    // Repeat r0 3 times (counts = 3 -> f2+, peak = 3) and r1 once (counts = 1 -> f1)
    // and r2 once (masked at idx 10 -> 0 k-mers at min_qual=20!).
    let mut fq = String::new();
    for rep in 0..3 {
        fq.push_str(&format!(
            "@r0_{rep}\n{}\n+\n{q20}\n",
            std::str::from_utf8(r0).unwrap()
        ));
    }
    fq.push_str(&format!(
        "@r1\n{}\n+\n{q21}\n",
        std::str::from_utf8(r1).unwrap()
    ));
    fq.push_str(&format!(
        "@r2_q19\n{}\n+\n{q_mask_str}\n",
        std::str::from_utf8(r2).unwrap()
    ));

    let fq_plain = dir.path().join("qual_edge.fq");
    fs::write(&fq_plain, &fq).unwrap();
    let fq_gz = dir.path().join("qual_edge.fq.gz");
    {
        let mut enc = GzEncoder::new(File::create(&fq_gz).unwrap(), Compression::fast());
        enc.write_all(fq.as_bytes()).unwrap();
        enc.finish().unwrap();
    }

    // Compute exact reference canonical k-mer counts at min_qual = 20.
    let space = KmerSpace::new(k);
    let mut exact: HashMap<u64, u32> = HashMap::new();
    for _ in 0..3 {
        space.for_each_canonical(r0, |x| *exact.entry(x).or_insert(0) += 1);
    }
    space.for_each_canonical(r1, |x| *exact.entry(x).or_insert(0) += 1);
    let exact_f0 = exact.len() as u64;
    let exact_f1 = exact.values().filter(|&&c| c == 1).count() as u64;
    let exact_f2 = exact.values().filter(|&&c| c >= 2).count() as u64;
    let exact_occ: u64 = exact.values().map(|&c| c as u64).sum();
    assert_eq!(exact_occ, 3 * 10 + 10);

    let pf_plain = preflight_inputs(std::slice::from_ref(&fq_plain), k, 20).unwrap();
    let pf_gz = preflight_inputs(&[fq_gz], k, 20).unwrap();
    assert_eq!(pf_plain, pf_gz, "plain and .gz FASTQ preflight must match");

    assert_eq!(pf_plain.sample, DEFAULT_SAMPLE_NAME);
    assert_eq!(pf_plain.total_records, 5);
    assert_eq!(pf_plain.total_bases, 100);
    // Only the single Q19 base at r2[10] is masked; all Q20 and Q21 bases stay valid!
    assert_eq!(pf_plain.valid_bases, 99);
    // Unmasked GC bases: 3*13 (r0) + 12 (r1) + 11 (r2 without masked G) = 62 of 99.
    let want_gc = 62.0 / 99.0;
    assert!(
        (pf_plain.gc_fraction - want_gc).abs() < 1e-12,
        "gc_fraction {} != exact {want_gc}",
        pf_plain.gc_fraction
    );
    assert!(pf_plain.sample_exact);
    assert_eq!(pf_plain.subsample_shift, 0);
    assert_eq!(pf_plain.total_kmer_occurrences, exact_occ);
    assert_eq!(pf_plain.f0_est, exact_f0);
    assert_eq!(pf_plain.f1_est, exact_f1);
    assert_eq!(pf_plain.f2_plus_est, exact_f2);
    assert_eq!(pf_plain.peak_coverage_est, 3);

    // At min_qual = 21, ALL Q20 bases (r0 and r2) are masked; only r1 (Q21) remains!
    let mut r1_exact: HashMap<u64, u32> = HashMap::new();
    space.for_each_canonical(r1, |x| *r1_exact.entry(x).or_insert(0) += 1);
    let r1_f0 = r1_exact.len() as u64;
    let r1_f1 = r1_exact.values().filter(|&&c| c == 1).count() as u64;
    let r1_f2 = r1_exact.values().filter(|&&c| c >= 2).count() as u64;
    let pf_q21 = preflight_inputs(&[fq_plain], k, 21).unwrap();
    assert_eq!(pf_q21.valid_bases, 20);
    assert_eq!(pf_q21.total_kmer_occurrences, 10);
    assert_eq!(pf_q21.f0_est, r1_f0);
    assert_eq!(pf_q21.f1_est, r1_f1);
    assert_eq!(pf_q21.f2_plus_est, r1_f2);
    assert!((pf_q21.gc_fraction - 12.0 / 20.0).abs() < 1e-12);
}

// ---------------------------------------------------------------------------
// 3. Small (LinearCounting) & Large (raw HLL + subsampling) cardinality checks
//    (kills preflight.rs:152 HLL rho, :169 LinearCounting, and :322 f0_est!)
// ---------------------------------------------------------------------------

#[test]
fn preflight_small_and_large_hll_estimates_verify_linear_counting_and_subsampling() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(141_421);
    let k = 31usize;
    let space = KmerSpace::new(k);

    // Case A (Small: F0 = 2,500 << 65,536, well inside HLL LinearCounting regime
    // raw <= 2.5 * 16,384 = 40,960):
    // Without LinearCounting (preflight.rs:169), raw HLL returns ~11,800 (+370% error!).
    let small_seq = rand_dna(&mut rng, 2_530);
    let small_fa = dir.path().join("small.fna");
    fs::write(&small_fa, [b">s\n".as_slice(), &small_seq, b"\n"].concat()).unwrap();
    let (small_vals, _) = kset::count_file(&small_fa, &space, 1).unwrap();
    let small_f0 = small_vals.len() as u64;
    assert_eq!(small_f0, 2_500);

    let pf_small = preflight_inputs(&[small_fa], k, 0).unwrap();
    assert!(pf_small.sample_exact);
    assert_eq!(pf_small.subsample_shift, 0);
    assert_eq!(pf_small.f0_est, small_f0);
    assert_eq!(pf_small.f1_est, small_f0);
    assert_eq!(pf_small.f2_plus_est, 0);
    let small_hll_rel = (pf_small.f0_hll_est as f64 - small_f0 as f64).abs() / small_f0 as f64;
    assert!(
        small_hll_rel < 0.05,
        "small HLL LinearCounting estimate {} too far from {small_f0} ({:.2}%)",
        pf_small.f0_hll_est,
        small_hll_rel * 100.0
    );

    // Case B (Large: F0 = 120,000 > 65,536, above LinearCounting cutoff 40,960):
    // Forces BOTH raw HLL estimation (preflight.rs:152) AND adaptive prefix
    // subsampling (`sample_exact == false`, `subsample_shift >= 1`, `f0_est == f0_hll_est`,
    // killing preflight.rs:322!).
    // Bias bases to 35% C, 25% G, 20% A, 20% T (60% GC) to also guard GC on large runs.
    let weighted_bases = b"AAAACCCCCCCGGGGGTTTT"; // 4 A, 7 C, 5 G, 4 T -> 60% GC
    let large_seq: Vec<u8> = (0..120_030)
        .map(|_| weighted_bases[rng.random_range(0..weighted_bases.len())])
        .collect();
    let large_fa = dir.path().join("large.fna");
    fs::write(&large_fa, [b">l\n".as_slice(), &large_seq, b"\n"].concat()).unwrap();
    let (large_vals, _) = kset::count_file(&large_fa, &space, 1).unwrap();
    let large_f0 = large_vals.len() as u64;

    let pf_large = preflight_inputs(&[large_fa], k, 0).unwrap();
    assert!(!pf_large.sample_exact);
    assert!(pf_large.subsample_shift >= 1);
    assert_eq!(pf_large.f0_est, pf_large.f0_hll_est);

    let large_hll_rel = (pf_large.f0_hll_est as f64 - large_f0 as f64).abs() / large_f0 as f64;
    assert!(
        large_hll_rel < 0.04,
        "large HLL estimate {} too far from exact {large_f0} ({:.2}%)",
        pf_large.f0_hll_est,
        large_hll_rel * 100.0
    );
    let large_f1_rel = (pf_large.f1_est as f64 - large_f0 as f64).abs() / large_f0 as f64;
    assert!(
        large_f1_rel < 0.05,
        "subsampled f1_est {} too far from {large_f0} ({:.2}%)",
        pf_large.f1_est,
        large_f1_rel * 100.0
    );
    assert!((pf_large.gc_fraction - 0.60).abs() < 0.01);
}

// ---------------------------------------------------------------------------
// 4. Disk footprint estimates vs ACTUAL written .kset, .grk, and .spss files
//    across BOTH clean multi-contig assembly AND 1%-error read set!
//    (kills preflight.rs:336, :345 and verifies Finding 16 modality scaling!)
// ---------------------------------------------------------------------------

#[test]
fn preflight_disk_size_estimates_match_actual_files_for_assembly_and_noisy_reads() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(314_159);
    let k = 31usize;
    let space = KmerSpace::new(k);

    // --- A) Clean multi-contig assembly (100 contigs x 1,000 bp) ---
    let asm_fa = dir.path().join("synth_asm.fna");
    {
        let mut w = std::io::BufWriter::new(File::create(&asm_fa).unwrap());
        for c in 0..100usize {
            writeln!(w, ">contig_{c}").unwrap();
            w.write_all(&rand_dna(&mut rng, 1_000)).unwrap();
            writeln!(w).unwrap();
        }
        w.flush().unwrap();
    }
    let (asm_vals, _) = kset::count_file(&asm_fa, &space, 1).unwrap();
    let asm_v1 = dir.path().join("asm.kset");
    kset::write(&asm_v1, DEFAULT_SAMPLE_NAME, &space, &asm_vals).unwrap();
    let asm_v1_bytes = fs::metadata(&asm_v1).unwrap().len();

    let asm_grk = dir.path().join("asm.grk");
    write_gr_kset(
        &asm_grk,
        DEFAULT_SAMPLE_NAME,
        k,
        InputModality::Assembly,
        &asm_vals,
    )
    .unwrap();
    let asm_grk_bytes = fs::metadata(&asm_grk).unwrap().len();

    let asm_spss = dir.path().join("asm.spss");
    let asm_arch = archive_from_fasta(
        &asm_fa,
        &asm_spss,
        k,
        DEFAULT_SAMPLE_NAME,
        InputModality::Assembly,
    )
    .unwrap();

    let pf_asm = preflight_inputs(std::slice::from_ref(&asm_fa), k, 0).unwrap();
    assert_eq!(pf_asm.spss_archive_profile, "assembly");
    let v1_err = (pf_asm.est_kset_v1_bytes.unwrap() as f64 - asm_v1_bytes as f64).abs()
        / asm_v1_bytes as f64;
    let grk_err =
        (pf_asm.est_gr_kset_bytes as f64 - asm_grk_bytes as f64).abs() / asm_grk_bytes as f64;
    let spss_asm_err = (pf_asm.est_spss_archive_bytes as f64 - asm_arch.archive_bytes as f64).abs()
        / asm_arch.archive_bytes as f64;
    assert!(
        v1_err < 0.12,
        "assembly .kset v1 estimate {} vs actual {asm_v1_bytes} ({:.2}%)",
        pf_asm.est_kset_v1_bytes.unwrap(),
        v1_err * 100.0
    );
    assert!(
        grk_err < 0.05,
        "assembly .grk estimate {} vs actual {asm_grk_bytes} ({:.2}%)",
        pf_asm.est_gr_kset_bytes,
        grk_err * 100.0
    );
    assert!(
        spss_asm_err < 0.10,
        "assembly .spss estimate {} vs actual {} ({:.2}%)",
        pf_asm.est_spss_archive_bytes,
        asm_arch.archive_bytes,
        spss_asm_err * 100.0
    );

    // At k = 61 (> 32), est_kset_v1_bytes MUST be None and est_gr_kset_bytes
    // must match actual k = 61 .grk bytes within 5%.
    let s61 = KmerSpace128::new(61);
    let mut keys61 = Vec::new();
    grm::fasta::for_each_record(&asm_fa, |_, seq| {
        s61.for_each_canonical(seq, |x| keys61.push(s61.mix(x)));
    })
    .unwrap();
    keys61.sort_unstable();
    keys61.dedup();
    let asm_grk61 = dir.path().join("asm61.grk");
    write_gr_kset(
        &asm_grk61,
        DEFAULT_SAMPLE_NAME,
        61,
        InputModality::Assembly,
        &keys61,
    )
    .unwrap();
    let asm_grk61_bytes = fs::metadata(&asm_grk61).unwrap().len();
    let pf_asm61 = preflight_inputs(&[asm_fa], 61, 0).unwrap();
    assert!(pf_asm61.est_kset_v1_bytes.is_none());
    let grk61_err =
        (pf_asm61.est_gr_kset_bytes as f64 - asm_grk61_bytes as f64).abs() / asm_grk61_bytes as f64;
    assert!(
        grk61_err < 0.05,
        "k=61 .grk estimate {} vs actual {asm_grk61_bytes} ({:.2}%)",
        pf_asm61.est_gr_kset_bytes,
        grk61_err * 100.0
    );

    // --- B) Synthetic 1%-error FASTQ read set (25 kb genome, ~20x depth, 150 bp reads) ---
    let genome = rand_dna(&mut rng, 25_000);
    let reads_fq = dir.path().join("noisy_reads.fq");
    let read_len = 150usize;
    let n_reads = 3_400usize; // ~20.4x depth
    let qual_line = "I".repeat(read_len);
    let mut exact_reads: HashMap<u64, u32> = HashMap::new();
    {
        let mut w = std::io::BufWriter::new(File::create(&reads_fq).unwrap());
        let mut buf = vec![0u8; read_len];
        for r in 0..n_reads {
            let start = rng.random_range(0..=(genome.len() - read_len));
            buf.copy_from_slice(&genome[start..start + read_len]);
            for base in &mut buf {
                // 1% substitution error rate.
                if rng.random_range(0..100u32) == 0 {
                    let alt = match *base {
                        b'A' => b"CGT"[rng.random_range(0..3)],
                        b'C' => b"AGT"[rng.random_range(0..3)],
                        b'G' => b"ACT"[rng.random_range(0..3)],
                        _ => b"ACG"[rng.random_range(0..3)],
                    };
                    *base = alt;
                }
            }
            space.for_each_canonical(&buf, |x| *exact_reads.entry(x).or_insert(0) += 1);
            writeln!(
                w,
                "@read_{r}\n{}\n+\n{qual_line}",
                std::str::from_utf8(&buf).unwrap()
            )
            .unwrap();
        }
        w.flush().unwrap();
    }

    let true_f0 = exact_reads.len() as u64;
    let true_f1 = exact_reads.values().filter(|&&c| c == 1).count() as u64;
    let true_f2 = exact_reads.values().filter(|&&c| c >= 2).count() as u64;

    let mut reads_keys: Vec<u64> = exact_reads.keys().map(|&x| space.mix(x)).collect();
    reads_keys.sort_unstable();
    let reads_grk = dir.path().join("reads.grk");
    write_gr_kset(
        &reads_grk,
        DEFAULT_SAMPLE_NAME,
        k,
        InputModality::Reads,
        &reads_keys,
    )
    .unwrap();
    let reads_grk_bytes = fs::metadata(&reads_grk).unwrap().len();

    let reads_spss = dir.path().join("reads.spss");
    let reads_arch = archive_from_fasta(
        &reads_fq,
        &reads_spss,
        k,
        DEFAULT_SAMPLE_NAME,
        InputModality::Reads,
    )
    .unwrap();

    let pf_reads = preflight_inputs(&[reads_fq], k, 0).unwrap();
    assert_eq!(pf_reads.spss_archive_profile, "reads");
    let f0_rel = (pf_reads.f0_est as f64 - true_f0 as f64).abs() / true_f0 as f64;
    let f1_rel = (pf_reads.f1_est as f64 - true_f1 as f64).abs() / true_f1 as f64;
    let f2_rel = (pf_reads.f2_plus_est as f64 - true_f2 as f64).abs() / true_f2 as f64;
    assert!(f0_rel < 0.04, "reads f0 error {:.2}%", f0_rel * 100.0);
    assert!(f1_rel < 0.05, "reads f1 error {:.2}%", f1_rel * 100.0);
    assert!(f2_rel < 0.08, "reads f2+ error {:.2}%", f2_rel * 100.0);
    assert!(
        (8..=22).contains(&pf_reads.peak_coverage_est),
        "expected peak coverage around 12..18 for 20x 1%-error reads, got {}",
        pf_reads.peak_coverage_est
    );

    let reads_grk_err =
        (pf_reads.est_gr_kset_bytes as f64 - reads_grk_bytes as f64).abs() / reads_grk_bytes as f64;
    let reads_spss_err = (pf_reads.est_spss_archive_bytes as f64 - reads_arch.archive_bytes as f64)
        .abs()
        / reads_arch.archive_bytes as f64;
    assert!(
        reads_grk_err < 0.05,
        "reads .grk estimate {} vs actual {reads_grk_bytes} ({:.2}%)",
        pf_reads.est_gr_kset_bytes,
        reads_grk_err * 100.0
    );
    assert!(
        reads_spss_err < 0.15,
        "reads .spss estimate {} vs actual {} ({:.2}%, actual {:.4} B/k-mer)",
        pf_reads.est_spss_archive_bytes,
        reads_arch.archive_bytes,
        reads_spss_err * 100.0,
        reads_arch.bytes_per_kmer
    );
    // Modality scaling check (Finding 16): reads per-k-mer SPSS estimate MUST
    // exceed the clean assembly per-k-mer estimate (not a fixed 0.30 B/k-mer!).
    let asm_est_bpk = pf_asm.est_spss_archive_bytes as f64 / pf_asm.f0_est as f64;
    let reads_est_bpk = pf_reads.est_spss_archive_bytes as f64 / pf_reads.f0_est as f64;
    assert!(
        reads_est_bpk > asm_est_bpk + 0.12,
        "reads est {reads_est_bpk:.4} B/k-mer must clearly exceed assembly est {asm_est_bpk:.4} B/k-mer"
    );
    assert!(
        reads_arch.bytes_per_kmer > 0.40 && asm_arch.bytes_per_kmer < 0.30,
        "actual archive B/k-mer: reads={:.4}, asm={:.4}",
        reads_arch.bytes_per_kmer,
        asm_arch.bytes_per_kmer
    );

    // Explicit modality override switches spss_archive_profile cleanly.
    let forced_reads_on_asm = preflight_inputs_with_options(
        &[dir.path().join("synth_asm.fna")],
        &PreflightOptions {
            sample: "custom_reads".to_string(),
            k,
            min_qual: 0,
            modality: InputModality::Reads,
        },
    )
    .unwrap();
    assert_eq!(forced_reads_on_asm.sample, "custom_reads");
    assert_eq!(forced_reads_on_asm.spss_archive_profile, "reads");
    assert!(forced_reads_on_asm.est_spss_archive_bytes > pf_asm.est_spss_archive_bytes);

    let forced_asm_on_reads = preflight_inputs_with_options(
        &[dir.path().join("noisy_reads.fq")],
        &PreflightOptions {
            sample: "custom_asm".to_string(),
            k,
            min_qual: 0,
            modality: InputModality::Assembly,
        },
    )
    .unwrap();
    assert_eq!(forced_asm_on_reads.sample, "custom_asm");
    assert_eq!(forced_asm_on_reads.spss_archive_profile, "assembly");
    assert_ne!(
        forced_asm_on_reads.est_spss_archive_bytes,
        pf_reads.est_spss_archive_bytes
    );
}

// ---------------------------------------------------------------------------
// 5. End-to-end CLI test for `grm kcount preflight`
//    - verifies all JSON fields
//    - defaults sample name to "sample" with zero filename stem leak
//    - rejects invalid k (4 and 65) and malformed FASTQ
// ---------------------------------------------------------------------------

#[test]
fn cli_preflight_json_output_default_sample_name_no_leak_and_error_cases() {
    let dir = tempfile::tempdir().unwrap();
    let secret_stem = "arbitrary_preflight_stem_88888";
    let fq_gz = dir.path().join(format!("{secret_stem}.fq.gz"));
    let seq = "ACGTACGTCCGGTTAACCGGTTAACCGGAATTCCGG"; // 36 bp (20 GC -> 20/36)
    assert_eq!(seq.len(), 36);
    assert_eq!(seq.bytes().filter(|&b| b == b'C' || b == b'G').count(), 20);
    let qual = "I".repeat(seq.len());
    let fq_text = format!("@r1\n{seq}\n+\n{qual}\n@r2\n{seq}\n+\n{qual}\n");
    {
        let mut enc = GzEncoder::new(File::create(&fq_gz).unwrap(), Compression::fast());
        enc.write_all(fq_text.as_bytes()).unwrap();
        enc.finish().unwrap();
    }

    // 1) Default CLI call (stdout JSON): sample == "sample", NO filename leak!
    let out = grm_bin()
        .args([
            "kcount",
            "preflight",
            "--input",
            fq_gz.to_str().unwrap(),
            "--k",
            "21",
            "--min-qual",
            "20",
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "preflight CLI failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    assert!(
        !stdout.contains(secret_stem),
        "preflight stdout leaked input filename stem"
    );
    let j: serde_json::Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(j["sample"], DEFAULT_SAMPLE_NAME);
    assert_eq!(j["k"], 21);
    assert_eq!(j["min_qual"], 20);
    assert_eq!(j["total_records"], 2);
    assert_eq!(j["total_bases"], 72);
    assert_eq!(j["valid_bases"], 72);
    assert!((j["gc_fraction"].as_f64().unwrap() - (20.0 / 36.0)).abs() < 1e-12);
    assert_eq!(j["total_kmer_occurrences"], 32);
    assert_eq!(j["sample_exact"], true);
    assert_eq!(j["subsample_shift"], 0);
    let space21 = KmerSpace::new(21);
    let mut seq_exact: HashMap<u64, u32> = HashMap::new();
    for _ in 0..2 {
        space21.for_each_canonical(seq.as_bytes(), |x| *seq_exact.entry(x).or_insert(0) += 1);
    }
    let want_f0 = seq_exact.len() as u64;
    let want_f1 = seq_exact.values().filter(|&&c| c == 1).count() as u64;
    let want_f2 = seq_exact.values().filter(|&&c| c >= 2).count() as u64;

    assert_eq!(j["f0_est"], want_f0);
    assert_eq!(j["f1_est"], want_f1);
    assert_eq!(j["f2_plus_est"], want_f2);
    assert_eq!(j["peak_coverage_est"], 2);
    assert_eq!(j["spss_archive_profile"], "reads");
    assert!(j["est_kset_v1_bytes"].as_u64().unwrap() > 0);
    assert!(j["est_gr_kset_bytes"].as_u64().unwrap() > 0);
    assert!(j["est_spss_archive_bytes"].as_u64().unwrap() > 0);
    assert!(j.get("path").is_none());
    assert!(j.get("input").is_none());

    // 2) Explicit --name, --modality, and --out file.
    let pf_out = dir.path().join("preflight.json");
    let out = grm_bin()
        .args([
            "kcount",
            "preflight",
            "--input",
            fq_gz.to_str().unwrap(),
            "--k",
            "31",
            "--name",
            "explicit_tag",
            "--modality",
            "assembly",
            "--out",
            pf_out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    let file_text = fs::read_to_string(&pf_out).unwrap();
    assert!(!file_text.contains(secret_stem));
    let j2: serde_json::Value = serde_json::from_str(&file_text).unwrap();
    assert_eq!(j2["sample"], "explicit_tag");
    assert_eq!(j2["k"], 31);
    assert_eq!(j2["spss_archive_profile"], "assembly");

    // 3) Invalid k (4 and 65) MUST fail cleanly.
    for bad_k in ["4", "65"] {
        let res = grm_bin()
            .args([
                "kcount",
                "preflight",
                "--input",
                fq_gz.to_str().unwrap(),
                "--k",
                bad_k,
            ])
            .output()
            .unwrap();
        assert!(!res.status.success(), "preflight must reject --k {bad_k}");
    }

    // 4) Malformed FASTQ (mismatched quality length) MUST fail cleanly.
    let bad_fq = dir.path().join("bad.fq");
    fs::write(&bad_fq, "@r1\nACGTACGT\n+\nIII\n").unwrap();
    let res = grm_bin()
        .args([
            "kcount",
            "preflight",
            "--input",
            bad_fq.to_str().unwrap(),
            "--k",
            "5",
        ])
        .output()
        .unwrap();
    assert!(
        !res.status.success(),
        "preflight must reject malformed FASTQ"
    );
}
