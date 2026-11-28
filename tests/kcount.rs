//! Comprehensive correctness, parity, 128-bit discrimination, uncapped counting,
//! dosage missingness, SPSS archive compression, spill parity, preflight, and
//! streaming set-algebra tests for `grm::kcount`.

use flate2::Compression;
use flate2::write::GzEncoder;
use grm::dna::KmerSpace;
use grm::kcount::dosage::{
    CorePanel, DosageConfig, DosageStatus, InputType, TargetPanel, extract_sample_dosage,
};
use grm::kcount::format::{
    GrKsetReader, KtabReader, convert_to_v1, convert_v1_to_gr, read_profiles,
};
use grm::kcount::key::{KmerSpace128, MAX_K};
use grm::kcount::ops::{
    MatrixExportConfig, MatrixFormat, MatrixMode, SetOp, export_matrix, run_set_op,
};
use grm::kcount::preflight::preflight_inputs;
use grm::kcount::spss::{archive_from_table, expand_to_kset_v1, expand_to_ktab};
use grm::kcount::superkmer::{CountConfig, count_inputs};
use grm::kcount::{KcountArgs, run_cli};
use grm::kset;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::collections::{BTreeMap, HashMap};
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn rc_dna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            other => other,
        })
        .collect()
}

fn rand_dna(rng: &mut ChaCha8Rng, len: usize) -> Vec<u8> {
    (0..len).map(|_| b"ACGT"[rng.random_range(0..4)]).collect()
}

// Exact bit-for-bit parity vs `crate::kset::count_into` and `kset::digest`
// for k = 5, 21, 31, 32 on multi-contig FASTA with Ns, repeats, and reverse
// complements, at both `min_count = 1` and `min_count = 3`.
#[test]
fn parity_with_kset_v1_across_u64_widths() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(42);
    let c1 = rand_dna(&mut rng, 12_000);
    let repeat = c1[500..650].to_vec();
    let repeat_rc = rc_dna(&repeat);
    let mut c2 = rand_dna(&mut rng, 5_000);
    c2.extend_from_slice(b"NNNNRY");
    c2.extend_from_slice(&repeat);
    c2.extend_from_slice(b"nnn");
    c2.extend_from_slice(&repeat_rc);
    c2.extend_from_slice(b"acgtACGT");
    c2.extend_from_slice(&repeat);

    let fna = dir.path().join("multi.fna");
    let mut f = File::create(&fna).unwrap();
    writeln!(f, ">contig1 desc").unwrap();
    for chunk in c1.chunks(80) {
        f.write_all(chunk).unwrap();
        writeln!(f).unwrap();
    }
    writeln!(f, ">contig2").unwrap();
    f.write_all(&c2).unwrap();
    writeln!(f).unwrap();
    drop(f);

    for k in [5usize, 21, 31, 32] {
        for min_count in [1u32, 3] {
            let space = KmerSpace::new(k);
            let mut want = Vec::new();
            let mut scratch = Vec::new();
            kset::count_into(&fna, &space, min_count, &mut want, &mut scratch).unwrap();
            let want_digest = kset::digest(&want);

            let out = dir.path().join(format!("k{k}_m{min_count}.ktab"));
            let cfg = CountConfig {
                k,
                sample_name: format!("s_{k}"),
                min_count,
                out_path: out.clone(),
                ..Default::default()
            };
            let summary = count_inputs(std::slice::from_ref(&fna), &cfg).unwrap();
            assert_eq!(
                summary.fnv64, want_digest,
                "digest mismatch k={k} m={min_count}"
            );
            assert_eq!(summary.distinct_retained as usize, want.len());

            let reader = KtabReader::open(&out).unwrap();
            let entries = reader.read_all::<u64>().unwrap();
            let got_keys: Vec<u64> = entries.iter().map(|&(key, _)| key).collect();
            assert_eq!(got_keys, want, "keys mismatch k={k} m={min_count}");
        }
    }
}

// Exhaustive 128-bit round-trip mixing and canonical decode test covering
// EVERY k from 33 to 64.
#[test]
fn u128_mix_unmix_and_canonical_decode_all_k_33_to_64() {
    let mut rng = ChaCha8Rng::seed_from_u64(3364);
    for k in 33..=MAX_K {
        let space = KmerSpace128::new(k);
        let mask = space.mask();
        assert_eq!(space.bits() as usize, 2 * k);

        // Edge values plus random draws across full 2k-bit space.
        let mut samples = vec![
            0u128,
            1,
            2,
            3,
            mask,
            mask - 1,
            mask >> 1,
            (mask >> 1) + 1,
            1u128 << 64,
            (1u128 << 64) - 1,
            (1u128 << 64) + 1,
        ];
        for _ in 0..500 {
            let v = (((rng.random::<u64>() as u128) << 64) | (rng.random::<u64>() as u128)) & mask;
            samples.push(v);
        }
        for x in samples {
            assert!(x <= mask);
            let y = space.mix(x);
            assert!(y <= mask, "k={k}: mix escaped mask");
            assert_eq!(space.unmix(y), x, "k={k}: unmix(mix(x)) != x");
            assert_eq!(space.mix(space.unmix(x)), x, "k={k}: mix(unmix(x)) != x");

            let decoded = space.decode(x);
            assert_eq!(decoded.len(), k);
            let rc_str = String::from_utf8(rc_dna(decoded.as_bytes())).unwrap();
            let rc_val = space.reverse_complement(x);
            assert_eq!(space.decode(rc_val), rc_str, "k={k}: SWAR rc mismatch");
            let canon = space.canonical(&decoded).unwrap();
            let canon_rc = space.canonical(&rc_str).unwrap();
            assert_eq!(canon, canon_rc, "k={k}: strand canonical mismatch");
            assert_eq!(canon, x.min(rc_val), "k={k}: canonical not min(fwd, rc)");
        }
    }
}

// 61-mer high-bit discrimination test: pairs identical in their first 32
// bases (differing only in bases 33..=61) AND pairs identical in their last
// 32 bases (identical low 64 bits, differing only above bit 64), plus reverse
// complements, verified through counting, `.ktab`, `spss-archive`, and `op`.
#[test]
fn k61_high_bit_and_low_bit_discrimination_zero_collisions() {
    let dir = tempfile::tempdir().unwrap();
    let space = KmerSpace128::new(61);

    // Pair A: identical first 32 bases, differing ONLY in bases 33..=61.
    let prefix32 = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let suffix29_a = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    let suffix29_b = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAC";
    let kmer_a1 = format!("{prefix32}{suffix29_a}");
    let kmer_a2 = format!("{prefix32}{suffix29_b}");
    assert_eq!(kmer_a1.len(), 61);
    assert_eq!(kmer_a2.len(), 61);
    assert_eq!(&kmer_a1[..32], &kmer_a2[..32]);

    // Pair B: identical last 32 bases (identical low 64 encoded bits!), differing ONLY in bases 1..=29.
    let suffix32 = "TGCATGCATGCATGCATGCATGCATGCATGCA";
    let prefix29_a = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
    let prefix29_b = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCG";
    let kmer_b1 = format!("{prefix29_a}{suffix32}");
    let kmer_b2 = format!("{prefix29_b}{suffix32}");
    assert_eq!(kmer_b1.len(), 61);
    assert_eq!(kmer_b2.len(), 61);
    assert_eq!(&kmer_b1[29..], &kmer_b2[29..]);

    // Verify all 4 canonical keys and their mixes are distinct, and confirm Pair B shares low 64 raw bits.
    let raw_pack = |s: &str| -> u128 {
        s.bytes().fold(0u128, |acc, b| {
            let c = match b {
                b'A' => 0u128,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => panic!("bad base"),
            };
            (acc << 2) | c
        })
    };
    assert_eq!(raw_pack(&kmer_b1) as u64, raw_pack(&kmer_b2) as u64);
    assert_ne!(raw_pack(&kmer_b1), raw_pack(&kmer_b2));

    let c_a1 = space.canonical(&kmer_a1).unwrap();
    let c_a2 = space.canonical(&kmer_a2).unwrap();
    let c_b1 = space.canonical(&kmer_b1).unwrap();
    let c_b2 = space.canonical(&kmer_b2).unwrap();
    let all_c = [c_a1, c_a2, c_b1, c_b2];
    for i in 0..4 {
        for j in (i + 1)..4 {
            assert_ne!(all_c[i], all_c[j]);
            assert_ne!(space.mix(all_c[i]), space.mix(all_c[j]));
        }
        assert_eq!(space.unmix(space.mix(all_c[i])), all_c[i]);
    }

    // Build two FASTAs (including reverse complements) and verify exact counting & ops.
    let rc_a1 = String::from_utf8(rc_dna(kmer_a1.as_bytes())).unwrap();
    let rc_b2 = String::from_utf8(rc_dna(kmer_b2.as_bytes())).unwrap();

    let fa1 = dir.path().join("set1.fna");
    fs::write(
        &fa1,
        format!(">a1\n{kmer_a1}\n>a1_rc\n{rc_a1}\n>a2\n{kmer_a2}\n>b1\n{kmer_b1}\n"),
    )
    .unwrap();
    let fa2 = dir.path().join("set2.fna");
    fs::write(
        &fa2,
        format!(">a2\n{kmer_a2}\n>b2\n{kmer_b2}\n>b2_rc\n{rc_b2}\n"),
    )
    .unwrap();

    let tab1 = dir.path().join("set1.ktab");
    let tab2 = dir.path().join("set2.ktab");
    let s1 = count_inputs(
        std::slice::from_ref(&fa1),
        &CountConfig {
            k: 61,
            sample_name: "set1".to_string(),
            out_path: tab1.clone(),
            ..Default::default()
        },
    )
    .unwrap();
    let s2 = count_inputs(
        std::slice::from_ref(&fa2),
        &CountConfig {
            k: 61,
            sample_name: "set2".to_string(),
            out_path: tab2.clone(),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(s1.distinct_retained, 3);
    assert_eq!(s2.distinct_retained, 2);

    let r1: HashMap<u128, u32> = KtabReader::open(&tab1)
        .unwrap()
        .read_all::<u128>()
        .unwrap()
        .into_iter()
        .collect();
    assert_eq!(r1[&space.mix(c_a1)], 2); // fwd + rc merged
    assert_eq!(r1[&space.mix(c_a2)], 1);
    assert_eq!(r1[&space.mix(c_b1)], 1);

    let r2: HashMap<u128, u32> = KtabReader::open(&tab2)
        .unwrap()
        .read_all::<u128>()
        .unwrap()
        .into_iter()
        .collect();
    assert_eq!(r2[&space.mix(c_a2)], 1);
    assert_eq!(r2[&space.mix(c_b2)], 2);

    // SPSS archive round-trip at k = 61.
    let spss1 = dir.path().join("set1.spss");
    let arch = archive_from_table(&tab1, &spss1).unwrap();
    assert_eq!(arch.n_kmers, 3);
    let exp_tab = dir.path().join("set1_exp.ktab");
    let exp = expand_to_ktab(&spss1, &exp_tab).unwrap();
    assert_eq!(exp.fnv64, arch.fnv64);
    let exp_keys: Vec<u128> = KtabReader::open(&exp_tab)
        .unwrap()
        .read_all::<u128>()
        .unwrap()
        .into_iter()
        .map(|(k, _)| k)
        .collect();
    let mut orig_keys: Vec<u128> = r1.keys().copied().collect();
    orig_keys.sort_unstable();
    assert_eq!(exp_keys, orig_keys);

    // Streaming set ops at k = 61 distinguish all 4 keys cleanly.
    let union_out = dir.path().join("u61.ktab");
    let inter_out = dir.path().join("i61.ktab");
    let diff_out = dir.path().join("d61.ktab");
    let sym_out = dir.path().join("s61.ktab");

    let su = run_set_op(
        SetOp::Union,
        &[tab1.clone(), tab2.clone()],
        &union_out,
        "u",
        false,
    )
    .unwrap();
    let si = run_set_op(
        SetOp::Intersect,
        &[tab1.clone(), tab2.clone()],
        &inter_out,
        "i",
        false,
    )
    .unwrap();
    let sd = run_set_op(
        SetOp::Diff,
        &[tab1.clone(), tab2.clone()],
        &diff_out,
        "d",
        false,
    )
    .unwrap();
    let ss = run_set_op(
        SetOp::SymDiff,
        &[tab1.clone(), tab2.clone()],
        &sym_out,
        "s",
        false,
    )
    .unwrap();

    assert_eq!(su.n_kmers_out, 4);
    assert_eq!(si.n_kmers_out, 1); // only a2
    assert_eq!(sd.n_kmers_out, 2); // a1, b1
    assert_eq!(ss.n_kmers_out, 3); // a1, b1, b2
    let inter_entries = KtabReader::open(&inter_out)
        .unwrap()
        .read_all::<u128>()
        .unwrap();
    assert_eq!(inter_entries, vec![(space.mix(c_a2), 1)]);
}

// Uncapped u32 counting (1, 2, 3, 255, 256, 300, 40000, 70000), exact
// 10.0x plasmid `count / D_core` recovery, explicit `None` on low core
// coverage (never silent zero), and core panel anti-leakage guardrails.
#[test]
fn uncapped_u32_counts_and_fold_safe_dosage_extraction() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(2026);
    let k = 21usize;
    let space = KmerSpace::new(k);

    // Generate 18 distinct non-overlapping 21-mers (10 core + 8 exact-count targets).
    let mut distinct_kmers: Vec<String> = Vec::new();
    let mut seen_canon = std::collections::HashSet::new();
    while distinct_kmers.len() < 18 {
        let s = String::from_utf8(rand_dna(&mut rng, k)).unwrap();
        let c = space.canonical(&s).unwrap();
        if seen_canon.insert(c) {
            distinct_kmers.push(space.decode(c));
        }
    }
    let core_kmers = distinct_kmers[0..10].to_vec();
    let test_counts: [u32; 8] = [1, 2, 3, 255, 256, 300, 40_000, 70_000];
    let target_kmers = distinct_kmers[10..18].to_vec();

    // Sample A (high core coverage: D_core = 30, plasmid target at 300 -> 10.0x dosage!).
    let fa_high = dir.path().join("high.fna");
    {
        let mut w = std::io::BufWriter::new(File::create(&fa_high).unwrap());
        for (i, ck) in core_kmers.iter().enumerate() {
            for r in 0..30 {
                writeln!(w, ">core_{i}_{r}\n{ck}").unwrap();
            }
        }
        for (i, (tk, &cnt)) in target_kmers.iter().zip(&test_counts).enumerate() {
            for r in 0..cnt {
                writeln!(w, ">t_{i}_{r}\n{tk}").unwrap();
            }
        }
        w.flush().unwrap();
    }

    let tab_high = dir.path().join("high.ktab");
    let meta_high = dir.path().join("high.meta.json");
    let sum_high = count_inputs(
        std::slice::from_ref(&fa_high),
        &CountConfig {
            k,
            sample_name: "high_cov".to_string(),
            out_path: tab_high.clone(),
            meta_out: Some(meta_high.clone()),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(sum_high.max_count_retained, 70_000);
    assert_eq!(sum_high.histogram.max_count, 70_000);

    // Verify exact uncapped u32 recovery from `.ktab` on disk.
    let map_high: HashMap<u64, u32> = KtabReader::open(&tab_high)
        .unwrap()
        .read_all::<u64>()
        .unwrap()
        .into_iter()
        .collect();
    for (tk, &want_c) in target_kmers.iter().zip(&test_counts) {
        let key = space.mix(space.canonical(tk).unwrap());
        assert_eq!(
            map_high[&key], want_c,
            "count saturation or mismatch for {want_c}"
        );
    }

    let core_panel = CorePanel::new(
        k,
        InputType::Reads,
        Some("train-fold-0".to_string()),
        core_kmers.clone(),
    )
    .unwrap();
    let mut targets_map = BTreeMap::new();
    targets_map.insert("plasmid_10x".to_string(), vec![target_kmers[5].clone()]); // count 300 / 30 = 10.0x
    targets_map.insert("high_copy".to_string(), vec![target_kmers[7].clone()]); // count 70000 / 30
    let target_panel = TargetPanel {
        k,
        targets: targets_map,
    };

    let report_high = extract_sample_dosage(
        &tab_high,
        &core_panel,
        Some(&target_panel),
        Some(&sum_high.qc),
        &DosageConfig::default(),
    )
    .unwrap();
    assert_eq!(report_high.dosage_status, DosageStatus::Ok);
    assert!((report_high.d_core - 30.0).abs() < 1e-12);
    assert!((report_high.core_breadth - 1.0).abs() < 1e-12);
    let p10 = report_high
        .targets
        .iter()
        .find(|t| t.target == "plasmid_10x")
        .unwrap();
    assert!((p10.normalized_dosage.unwrap() - 10.0).abs() < 1e-12);
    assert!((p10.kmers[0].dosage.unwrap() - 10.0).abs() < 1e-12);

    // Sample B (low core depth D_core = 4.0 < 10.0, even though plasmid is present at 300!):
    // MUST emit `insufficient_core_coverage` and `None` (JSON null), NEVER 0.0!
    let fa_low = dir.path().join("low.fna");
    {
        let mut w = std::io::BufWriter::new(File::create(&fa_low).unwrap());
        for (i, ck) in core_kmers.iter().enumerate() {
            for r in 0..4 {
                writeln!(w, ">core_{i}_{r}\n{ck}").unwrap();
            }
        }
        for r in 0..300 {
            writeln!(w, ">p_{r}\n{}", target_kmers[5]).unwrap();
        }
        w.flush().unwrap();
    }
    let tab_low = dir.path().join("low.ktab");
    let sum_low = count_inputs(
        std::slice::from_ref(&fa_low),
        &CountConfig {
            k,
            sample_name: "low_cov".to_string(),
            out_path: tab_low.clone(),
            ..Default::default()
        },
    )
    .unwrap();
    let report_low = extract_sample_dosage(
        &tab_low,
        &core_panel,
        Some(&target_panel),
        Some(&sum_low.qc),
        &DosageConfig::default(),
    )
    .unwrap();
    assert_eq!(
        report_low.dosage_status,
        DosageStatus::InsufficientCoreCoverage
    );
    assert!((report_low.d_core - 4.0).abs() < 1e-12);
    for t in &report_low.targets {
        assert!(
            t.normalized_dosage.is_none(),
            "low coverage must emit None, never 0.0"
        );
        assert!(t.detected_normalized_dosage.is_none());
        for ke in &t.kmers {
            assert!(ke.dosage.is_none());
        }
    }
    let json_val = serde_json::to_value(&report_low).unwrap();
    assert_eq!(json_val["dosage_status"], "insufficient_core_coverage");
    assert!(json_val["targets"][0]["normalized_dosage"].is_null());

    // Anti-leakage guardrails: mismatched input_type and mismatched k must fail.
    let bad_type_cfg = DosageConfig {
        input_type: InputType::Assembly,
        ..Default::default()
    };
    assert!(extract_sample_dosage(&tab_high, &core_panel, None, None, &bad_type_cfg).is_err());
    let wrong_k_panel = CorePanel::new(
        31,
        InputType::Reads,
        None,
        vec!["ACGTACGTACGTACGTACGTACGTACGTACG".to_string()],
    )
    .unwrap();
    assert!(
        extract_sample_dosage(
            &tab_high,
            &wrong_k_panel,
            None,
            None,
            &DosageConfig::default()
        )
        .is_err()
    );

    // Verify matrix-export --mode dosage preserves explicit missingness (NA + missing=1).
    let mat_dir = dir.path().join("mat_dosage_dense");
    let manifest = export_matrix(
        &[tab_high.clone(), tab_low.clone()],
        &mat_dir,
        &MatrixExportConfig {
            mode: MatrixMode::Dosage,
            format: MatrixFormat::Dense,
            min_support: 1,
            max_support: None,
            core_panel: Some(core_panel.clone()),
            dosage: DosageConfig::default(),
        },
    )
    .unwrap();
    assert_eq!(manifest.n_samples, 2);
    assert_eq!(manifest.n_missing_samples, 1);
    assert!(!manifest.samples[0].missing);
    assert!(manifest.samples[1].missing);
    let mat_text = fs::read_to_string(mat_dir.join("matrix.tsv")).unwrap();
    assert!(mat_text.contains("\tNA\n"));

    let mat_sparse_dir = dir.path().join("mat_dosage_sparse");
    let manifest_sp = export_matrix(
        &[tab_high, tab_low],
        &mat_sparse_dir,
        &MatrixExportConfig {
            mode: MatrixMode::Dosage,
            format: MatrixFormat::Sparse,
            min_support: 1,
            max_support: None,
            core_panel: Some(core_panel),
            dosage: DosageConfig::default(),
        },
    )
    .unwrap();
    assert_eq!(manifest_sp.n_missing_samples, 1);
    let coo_text = fs::read_to_string(mat_sparse_dir.join("matrix.coo.tsv")).unwrap();
    assert!(coo_text.contains("\t1\tNA\t1\n"));
}

// Compression & speed benchmark test comparing `.kset` v1 (`GRMKSET1`),
// Golomb-Rice sorted keys (`GRMKGR1`), and experimental 2-bit SPSS archive
// (`GRMSPS1`) on a synthetic multi-contig bacterial-density genome.
#[test]
fn spss_archive_beats_sorted_integer_bound_and_round_trips_kset_v1() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(314159);
    let k = 31usize;
    let space = KmerSpace::new(k);

    // 400,000 bp synthetic multi-contig bacterial-style genome with repeats and Ns.
    let fna = dir.path().join("synth_bact.fna");
    {
        let mut w = std::io::BufWriter::new(File::create(&fna).unwrap());
        let repeat = rand_dna(&mut rng, 1_500);
        for c in 0..8 {
            let mut contig = rand_dna(&mut rng, 48_000);
            contig.extend_from_slice(&repeat);
            contig.extend_from_slice(b"NNNNN");
            contig.extend_from_slice(&rand_dna(&mut rng, 500));
            writeln!(w, ">contig_{c}").unwrap();
            w.write_all(&contig).unwrap();
            writeln!(w).unwrap();
        }
        w.flush().unwrap();
    }

    let (values, _) = kset::count_file(&fna, &space, 1).unwrap();
    let orig_digest = kset::digest(&values);
    let v1_path = dir.path().join("genome.kset");
    kset::write(&v1_path, "synth", &space, &values).unwrap();
    let v1_bytes = fs::metadata(&v1_path).unwrap().len();
    let v1_bpk = v1_bytes as f64 / values.len() as f64;

    // Golomb-Rice sorted mixed integer baseline (GRMKGR1).
    let gr_path = dir.path().join("genome.grk");
    let n_gr = convert_v1_to_gr(&v1_path, &gr_path).unwrap();
    assert_eq!(n_gr as usize, values.len());
    let gr_bytes = fs::metadata(&gr_path).unwrap().len();
    let gr_bpk = gr_bytes as f64 / values.len() as f64;
    assert_eq!(
        GrKsetReader::open(&gr_path)
            .unwrap()
            .read_all::<u64>()
            .unwrap(),
        values
    );
    assert!(
        gr_bytes < v1_bytes,
        "Golomb-Rice ({gr_bpk:.3} B/k-mer) must beat kset v1 ({v1_bpk:.3} B/k-mer)"
    );

    // Convert GRMKGR1 back to GRMKSET1 losslessly.
    let gr_back_v1 = dir.path().join("from_gr.kset");
    convert_to_v1(&gr_path, &gr_back_v1).unwrap();
    assert_eq!(fs::read(&gr_back_v1).unwrap(), fs::read(&v1_path).unwrap());

    // Experimental SPSS greedy simplitig archive (GRMSPS1).
    let spss_path = dir.path().join("genome.spss");
    let arch = archive_from_table(&v1_path, &spss_path).unwrap();
    assert_eq!(arch.fnv64, orig_digest);
    assert!(
        arch.bytes_per_kmer < 0.40,
        "SPSS archive expected < 0.40 B/k-mer, got {:.4} B/k-mer",
        arch.bytes_per_kmer
    );
    assert!(arch.archive_bytes * 10 < v1_bytes);

    let rebuilt_v1 = dir.path().join("rebuilt.kset");
    let exp = expand_to_kset_v1(&spss_path, &rebuilt_v1).unwrap();
    assert_eq!(exp.fnv64, orig_digest);
    assert_eq!(fs::read(&rebuilt_v1).unwrap(), fs::read(&v1_path).unwrap());

    eprintln!(
        "compression benchmark (N={} k={k}): kset_v1={:.3} B/k-mer, gr_sorted={:.3} B/k-mer, spss_archive={:.4} B/k-mer (simplitigs={}, archive_ms={:.1}, expand_ms={:.1})",
        values.len(),
        v1_bpk,
        gr_bpk,
        arch.bytes_per_kmer,
        arch.n_simplitigs,
        arch.elapsed_ms,
        exp.elapsed_ms
    );
}

// Spill-to-disk parity: forcing `--max-memory-mb 1` triggers partition
// spilling and yields bit-for-bit identical digests, counts, and histogram.
#[test]
fn bounded_memory_spill_matches_in_memory_bit_for_bit() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(808);
    let fna = dir.path().join("big.fna");
    {
        let mut w = std::io::BufWriter::new(File::create(&fna).unwrap());
        // ~2.2 MB of sequence (~650 KB of 2-bit super-k-mers) exceeds 512 KB budget at 1 MiB.
        for i in 0..44 {
            writeln!(w, ">c{i}").unwrap();
            w.write_all(&rand_dna(&mut rng, 50_000)).unwrap();
            writeln!(w).unwrap();
        }
        w.flush().unwrap();
    }

    let mem_out = dir.path().join("mem.ktab");
    let spill_out = dir.path().join("spill.ktab");
    let s_mem = count_inputs(
        std::slice::from_ref(&fna),
        &CountConfig {
            k: 31,
            sample_name: "spilltest".to_string(),
            max_memory_mb: 256,
            out_path: mem_out.clone(),
            ..Default::default()
        },
    )
    .unwrap();
    let s_spill = count_inputs(
        std::slice::from_ref(&fna),
        &CountConfig {
            k: 31,
            sample_name: "spilltest".to_string(),
            max_memory_mb: 1,
            out_path: spill_out.clone(),
            tmp_dir: Some(dir.path().join("tmp")),
            ..Default::default()
        },
    )
    .unwrap();

    assert!(!s_mem.spilled_to_disk);
    assert!(s_spill.spilled_to_disk);
    assert_eq!(s_mem.fnv64, s_spill.fnv64);
    assert_eq!(s_mem.keys_sha256, s_spill.keys_sha256);
    assert_eq!(s_mem.counts_sha256, s_spill.counts_sha256);
    assert_eq!(s_mem.histogram, s_spill.histogram);
    assert_eq!(s_mem.sketch_hex, s_spill.sketch_hex);
    assert_eq!(fs::read(&mem_out).unwrap(), fs::read(&spill_out).unwrap());
}

// FASTQ quality masking (`--min-qual`), gzip `.fq.gz` parity, read QC
// metrics, FastK-style `.prof.bin` per-sequence profile, and preflight.
#[test]
fn fastq_quality_masking_gzip_profile_and_preflight() {
    let dir = tempfile::tempdir().unwrap();
    let k = 15usize;
    let base_seq = "AACCTTAGGCCTGAACTGGTTCCAGGCTAATC"; // 32 bp, non-periodic
    let high_q = "I".repeat(base_seq.len()); // Phred 40
    // Low-quality base at index 15 ('!' = Phred 0) splits into two runs: [0..15] (len 15) and [16..32] (len 16).
    let mut masked_q = high_q.clone().into_bytes();
    masked_q[15] = b'!';
    let masked_q_str = String::from_utf8(masked_q).unwrap();

    let mut fq_text = String::new();
    for i in 0..5 {
        fq_text.push_str(&format!("@solid_{i}\n{base_seq}\n+\n{high_q}\n"));
    }
    fq_text.push_str(&format!("@masked_read\n{base_seq}\n+\n{masked_q_str}\n"));
    // Add one singleton read with unique k-mers.
    let uniq_seq = "ATATATATATATATAGCGCTCGATCGATCGAT";
    fq_text.push_str(&format!(
        "@uniq\n{uniq_seq}\n+\n{}\n",
        "I".repeat(uniq_seq.len())
    ));

    let fq_plain = dir.path().join("reads.fq");
    fs::write(&fq_plain, &fq_text).unwrap();
    let fq_gz = dir.path().join("reads.fq.gz");
    {
        let mut enc = GzEncoder::new(File::create(&fq_gz).unwrap(), Compression::fast());
        enc.write_all(fq_text.as_bytes()).unwrap();
        enc.finish().unwrap();
    }

    let out_plain = dir.path().join("plain.ktab");
    let out_gz = dir.path().join("gz.ktab");
    let prof_path = dir.path().join("reads.prof.bin");
    let hist_path = dir.path().join("reads.hist.json");

    let s_plain = count_inputs(
        std::slice::from_ref(&fq_plain),
        &CountConfig {
            k,
            sample_name: "fq".to_string(),
            min_count: 1,
            min_qual: 20,
            out_path: out_plain.clone(),
            hist_out: Some(hist_path.clone()),
            profile_out: Some(prof_path.clone()),
            ..Default::default()
        },
    )
    .unwrap();
    let s_gz = count_inputs(
        std::slice::from_ref(&fq_gz),
        &CountConfig {
            k,
            sample_name: "fq".to_string(),
            min_count: 1,
            min_qual: 20,
            out_path: out_gz.clone(),
            ..Default::default()
        },
    )
    .unwrap();

    assert_eq!(s_plain.fnv64, s_gz.fnv64);
    assert_eq!(s_plain.counts_sha256, s_gz.counts_sha256);
    assert!(s_plain.qc.is_fastq);
    assert_eq!(s_plain.qc.total_reads, 7);
    assert_eq!(s_plain.qc.masked_bases, 1);
    assert!(s_plain.qc.singleton_fraction > 0.0);

    // Check per-sequence profile: masked_read must have count 0 on positions crossing index 15!
    let (pk, profiles) = read_profiles(&prof_path).unwrap();
    assert_eq!(pk, k);
    assert_eq!(profiles.len(), 7);
    let masked_prof = &profiles[5];
    assert_eq!(masked_prof.name, "masked_read");
    assert_eq!(masked_prof.counts.len(), base_seq.len() - k + 1);
    assert_eq!(masked_prof.counts[0], 6); // 5 solid + 1 valid prefix [0..15]
    for pos in 1..16 {
        assert_eq!(masked_prof.counts[pos], 0, "pos {pos} spans masked base 15");
    }
    assert_eq!(masked_prof.counts[16], 6); // valid suffix [16..31]
    assert_eq!(masked_prof.counts[17], 6); // valid suffix [17..32]

    // Verify preflight agrees closely with exact histogram.
    let pf = preflight_inputs(&[fq_gz], k, 20).unwrap();
    assert!(pf.sample_exact);
    assert_eq!(pf.f0_est, s_plain.histogram.distinct_total);
    assert_eq!(pf.f1_est, s_plain.histogram.f1);
    assert_eq!(pf.f2_plus_est, s_plain.histogram.f2_plus);
    assert_eq!(
        pf.total_kmer_occurrences,
        s_plain.histogram.total_occurrences
    );
}

// End-to-end CLI smoke test exercising every `grm kcount` subcommand via
// `clap::Parser` + `grm::kcount::run_cli`.
#[test]
fn cli_subcommands_end_to_end() {
    use clap::Parser;

    #[derive(Parser)]
    struct Wrapper {
        #[command(flatten)]
        args: KcountArgs,
    }

    let dir = tempfile::tempdir().unwrap();
    let fna = dir.path().join("mini.fna");
    let mut rng = ChaCha8Rng::seed_from_u64(999);
    let seq = rand_dna(&mut rng, 4_000);
    fs::write(&fna, [b">c1\n".as_slice(), &seq, b"\n"].concat()).unwrap();

    let tab = dir.path().join("mini.ktab");
    let spss = dir.path().join("mini.spss");
    let kset_v1 = dir.path().join("mini.kset");
    let gr = dir.path().join("mini.grk");
    let pf_out = dir.path().join("pf.json");
    let insp_out = dir.path().join("inspect.json");

    let run = |argv: &[&str]| {
        let w = Wrapper::parse_from(argv);
        run_cli(w.args).unwrap();
    };

    run(&[
        "kcount",
        "count",
        "--input",
        fna.to_str().unwrap(),
        "--k",
        "31",
        "--out",
        tab.to_str().unwrap(),
    ]);
    run(&[
        "kcount",
        "archive",
        "--input",
        tab.to_str().unwrap(),
        "--out",
        spss.to_str().unwrap(),
    ]);
    run(&[
        "kcount",
        "expand",
        "--input",
        spss.to_str().unwrap(),
        "--out",
        kset_v1.to_str().unwrap(),
        "--format",
        "kset-v1",
    ]);
    run(&[
        "kcount",
        "convert-v1",
        "--input",
        kset_v1.to_str().unwrap(),
        "--out",
        gr.to_str().unwrap(),
    ]);
    run(&[
        "kcount",
        "preflight",
        "--input",
        fna.to_str().unwrap(),
        "--k",
        "31",
        "--out",
        pf_out.to_str().unwrap(),
    ]);
    run(&[
        "kcount",
        "inspect",
        "--input",
        spss.to_str().unwrap(),
        "--verify",
        "--out",
        insp_out.to_str().unwrap(),
    ]);

    assert!(tab.exists() && spss.exists() && kset_v1.exists() && gr.exists());
    let insp: serde_json::Value = serde_json::from_slice(&fs::read(&insp_out).unwrap()).unwrap();
    assert_eq!(insp["verified"], true);
    assert_eq!(insp["k"], 31);
    let _ = PathBuf::from(&insp_out);
}

// Streaming set algebra across mixed `.kset` v1, `.ktab`, and `.grk` inputs,
// plus deterministic `pa` and `counts` matrix exports.
#[test]
fn mixed_format_set_ops_and_pa_counts_matrix_export() {
    let dir = tempfile::tempdir().unwrap();
    let k = 15usize;
    let space = KmerSpace::new(k);
    let k1 = "AAAAAAAAAAAAAAA";
    let k2 = "CCCCCCCCCCCCCCC";
    let k3 = "ACGTACGTACGTACG";
    let k4 = "TGCATGCATGCATGC";

    let fa_a = dir.path().join("a.fna");
    let fa_b = dir.path().join("b.fna");
    let fa_c = dir.path().join("c.fna");
    fs::write(&fa_a, format!(">1\n{k1}\n>2\n{k2}\n>2b\n{k2}\n>3\n{k3}\n")).unwrap();
    fs::write(&fa_b, format!(">2\n{k2}\n>3\n{k3}\n>4\n{k4}\n")).unwrap();
    fs::write(&fa_c, format!(">3\n{k3}\n>4\n{k4}\n")).unwrap();

    // Sample A as `.ktab`, Sample B as `.kset` v1, Sample C as `.grk`.
    let tab_a = dir.path().join("a.ktab");
    count_inputs(
        &[fa_a],
        &CountConfig {
            k,
            sample_name: "a".to_string(),
            out_path: tab_a.clone(),
            ..Default::default()
        },
    )
    .unwrap();

    let (vals_b, _) = kset::count_file(&fa_b, &space, 1).unwrap();
    let v1_b = dir.path().join("b.kset");
    kset::write(&v1_b, "b", &space, &vals_b).unwrap();

    let gr_c = dir.path().join("c.grk");
    count_inputs(
        &[fa_c],
        &CountConfig {
            k,
            sample_name: "c".to_string(),
            presence_only_out: true,
            out_path: gr_c.clone(),
            ..Default::default()
        },
    )
    .unwrap();

    let inputs = vec![tab_a.clone(), v1_b.clone(), gr_c.clone()];
    let u_out = dir.path().join("u.ktab");
    let i_out = dir.path().join("i.ktab");
    let d_out = dir.path().join("d.ktab");
    let s_out = dir.path().join("s.ktab");

    let su = run_set_op(SetOp::Union, &inputs, &u_out, "u", false).unwrap();
    let si = run_set_op(SetOp::Intersect, &inputs, &i_out, "i", false).unwrap();
    let sd = run_set_op(SetOp::Diff, &inputs, &d_out, "d", false).unwrap();
    let ss = run_set_op(SetOp::SymDiff, &inputs, &s_out, "s", false).unwrap();

    assert_eq!(su.n_kmers_out, 4); // k1, k2, k3, k4
    assert_eq!(si.n_kmers_out, 1); // only k3 in all three
    assert_eq!(sd.n_kmers_out, 1); // only k1 in A \ (B U C)
    assert_eq!(ss.n_kmers_out, 2); // k1 (in 1), k3 (in 3) -> odd count!

    let m_pa = export_matrix(
        &inputs,
        &dir.path().join("mat_pa"),
        &MatrixExportConfig {
            mode: MatrixMode::Pa,
            format: MatrixFormat::Dense,
            min_support: 1,
            max_support: None,
            core_panel: None,
            dosage: DosageConfig::default(),
        },
    )
    .unwrap();
    assert_eq!(m_pa.n_features, 4);
    assert_eq!(m_pa.n_missing_samples, 0);

    let m_cnt = export_matrix(
        &inputs,
        &dir.path().join("mat_cnt"),
        &MatrixExportConfig {
            mode: MatrixMode::Counts,
            format: MatrixFormat::Sparse,
            min_support: 2,
            max_support: None,
            core_panel: None,
            dosage: DosageConfig::default(),
        },
    )
    .unwrap();
    assert_eq!(m_cnt.n_features, 3); // k2, k3, k4 have support >= 2
}

// Corrupt SPSS archives (flipped base / broken digest) and truncated tables are rejected.
#[test]
fn corrupt_archive_and_tables_are_rejected() {
    let dir = tempfile::tempdir().unwrap();
    let fna = dir.path().join("g.fna");
    let mut rng = ChaCha8Rng::seed_from_u64(555);
    fs::write(
        &fna,
        [b">c\n".as_slice(), &rand_dna(&mut rng, 2_000), b"\n"].concat(),
    )
    .unwrap();

    let tab = dir.path().join("g.ktab");
    count_inputs(
        &[fna],
        &CountConfig {
            k: 31,
            sample_name: "g".to_string(),
            out_path: tab.clone(),
            ..Default::default()
        },
    )
    .unwrap();

    let spss = dir.path().join("g.spss");
    archive_from_table(&tab, &spss).unwrap();
    let mut corrupted = fs::read(&spss).unwrap();
    let last = corrupted.len() - 1;
    corrupted[last] ^= 0x03;
    let bad_spss = dir.path().join("bad.spss");
    fs::write(&bad_spss, &corrupted).unwrap();
    assert!(expand_to_kset_v1(&bad_spss, &dir.path().join("out.kset")).is_err());

    let tab_bytes = fs::read(&tab).unwrap();
    let trunc_tab = dir.path().join("trunc.ktab");
    fs::write(&trunc_tab, &tab_bytes[..tab_bytes.len() - 3]).unwrap();
    assert!(KtabReader::open(&trunc_tab).is_err());
}
