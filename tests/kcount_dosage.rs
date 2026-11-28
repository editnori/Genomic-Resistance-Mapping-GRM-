//! Integration and CLI tests for modality-verified core-panel dosage
//! extraction (`grm kcount dosage`), training-fold core-panel construction
//! (`grm kcount build-core-panel`), and dosage matrix export
//! (`grm kcount matrix-export --mode dosage`).

use grm::dna::KmerSpace;
use grm::kcount::dosage::{
    CorePanel, CorePanelBuildConfig, DosageConfig, DosageQcVector, DosageStatus, InputType,
    TargetPanel, build_training_fold_core_panel, extract_sample_dosage, load_qc_metrics,
};
use grm::kcount::format::{DEFAULT_SAMPLE_NAME, InputModality, write_gr_kset, write_ktab};
use grm::kcount::key::{KmerSpace64, KmerSpace128};
use grm::kcount::ops::{MatrixExportConfig, MatrixLayout, MatrixMode, export_matrix};
use grm::kset;
use grm::matrix::REVIEWED_COLUMNS;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn grm_bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_grm"))
}

fn gen_distinct_kmers_u64(seed: u64, k: usize, n: usize) -> Vec<String> {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let sp = KmerSpace64::new(k);
    let mut seen = HashSet::new();
    let mut out = Vec::with_capacity(n);
    while out.len() < n {
        let s: String = (0..k)
            .map(|_| b"ACGT"[rng.random_range(0..4)] as char)
            .collect();
        let c = sp.canonical(&s).unwrap();
        if seen.insert(c) {
            out.push(sp.decode(c));
        }
    }
    out
}

fn gen_distinct_kmers_u128(seed: u64, k: usize, n: usize) -> Vec<String> {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let sp = KmerSpace128::new(k);
    let mut seen = HashSet::new();
    let mut out = Vec::with_capacity(n);
    while out.len() < n {
        let s: String = (0..k)
            .map(|_| b"ACGT"[rng.random_range(0..4)] as char)
            .collect();
        let c = sp.canonical(&s).unwrap();
        if seen.insert(c) {
            out.push(sp.decode(c));
        }
    }
    out
}

fn write_ktab_u64(
    path: &Path,
    name: &str,
    k: usize,
    modality: InputModality,
    items: &[(String, u32)],
) {
    let sp = KmerSpace64::new(k);
    let mut entries: Vec<(u64, u32)> = items
        .iter()
        .map(|(s, c)| (sp.mix(sp.canonical(s).unwrap()), *c))
        .collect();
    entries.sort_unstable_by_key(|&(key, _)| key);
    write_ktab(path, name, k, modality, 1, 0, &entries).unwrap();
}

fn write_ktab_u128(
    path: &Path,
    name: &str,
    k: usize,
    modality: InputModality,
    items: &[(String, u32)],
) {
    let sp = KmerSpace128::new(k);
    let mut entries: Vec<(u128, u32)> = items
        .iter()
        .map(|(s, c)| (sp.mix(sp.canonical(s).unwrap()), *c))
        .collect();
    entries.sort_unstable_by_key(|&(key, _)| key);
    write_ktab(path, name, k, modality, 1, 0, &entries).unwrap();
}

#[test]
fn modality_sha256_fold_id_and_k_mismatch_guardrails() {
    let dir = tempfile::tempdir().unwrap();
    let k = 21;
    let kmers = gen_distinct_kmers_u64(7001, k, 12);
    let core_kmers = kmers[0..8].to_vec();
    let target_kmers = kmers[8..12].to_vec();

    let panel_reads = CorePanel::new(
        k,
        InputType::Reads,
        Some("fold-3".to_string()),
        core_kmers.clone(),
    )
    .unwrap();
    let panel_asm = CorePanel::new(
        k,
        InputType::Assembly,
        Some("fold-3".to_string()),
        core_kmers.clone(),
    )
    .unwrap();

    let items: Vec<(String, u32)> = core_kmers.iter().map(|s| (s.clone(), 20)).collect();
    let reads_tab = dir.path().join("sample_reads.ktab");
    let asm_tab = dir.path().join("sample_asm.ktab");
    let unk_tab = dir.path().join("sample_unk.ktab");
    write_ktab_u64(&reads_tab, "r", k, InputModality::Reads, &items);
    write_ktab_u64(&asm_tab, "a", k, InputModality::Assembly, &items);
    write_ktab_u64(&unk_tab, "u", k, InputModality::Unknown, &items);

    // 1) Reads .ktab MUST be rejected against Assembly panel EVEN WHEN
    //    expected_input_type == Some(InputType::Assembly) agrees with panel!
    let err = extract_sample_dosage(
        &reads_tab,
        &panel_asm,
        None,
        None,
        &DosageConfig {
            expected_input_type: Some(InputType::Assembly),
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(
        format!("{err:#}").contains("modality mismatch"),
        "expected stored table modality mismatch error, got: {err:#}"
    );

    // 2) Assembly .ktab rejected against Reads panel even when flag agrees with panel.
    let err = extract_sample_dosage(
        &asm_tab,
        &panel_reads,
        None,
        None,
        &DosageConfig {
            expected_input_type: Some(InputType::Reads),
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(format!("{err:#}").contains("modality mismatch"));

    // 3) Unknown modality .ktab rejected against both Reads and Assembly panels.
    for panel in [&panel_reads, &panel_asm] {
        let err = extract_sample_dosage(&unk_tab, panel, None, None, &DosageConfig::default())
            .unwrap_err();
        assert!(
            format!("{err:#}").contains("unknown input modality"),
            "expected unknown modality rejection, got: {err:#}"
        );
    }

    // 4) Matching table modality with mismatched CLI flag rejected.
    let err = extract_sample_dosage(
        &reads_tab,
        &panel_reads,
        None,
        None,
        &DosageConfig {
            expected_input_type: Some(InputType::Assembly),
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(format!("{err:#}").contains("input_type flag mismatch"));

    // 5) expected_fold_id match succeeds; mismatch or missing fold_id fails.
    let ok = extract_sample_dosage(
        &reads_tab,
        &panel_reads,
        None,
        None,
        &DosageConfig {
            expected_fold_id: Some("fold-3".to_string()),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(ok.dosage_status, DosageStatus::Ok);
    assert_eq!(ok.fold_id.as_deref(), Some("fold-3"));
    assert_eq!(ok.table_modality, InputModality::Reads);

    let err_fold = extract_sample_dosage(
        &reads_tab,
        &panel_reads,
        None,
        None,
        &DosageConfig {
            expected_fold_id: Some("fold-9".to_string()),
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(format!("{err_fold:#}").contains("fold_id mismatch"));

    let panel_no_fold = CorePanel::new(k, InputType::Reads, None, core_kmers.clone()).unwrap();
    assert!(
        extract_sample_dosage(
            &reads_tab,
            &panel_no_fold,
            None,
            None,
            &DosageConfig {
                expected_fold_id: Some("fold-3".to_string()),
                ..Default::default()
            },
        )
        .is_err()
    );

    // 6) Missing/empty/mismatched panel_sha256 rejected unconditionally.
    let mut bad_sha_panel = panel_reads.clone();
    bad_sha_panel.panel_sha256 = "f".repeat(64);
    let err_sha = extract_sample_dosage(
        &reads_tab,
        &bad_sha_panel,
        None,
        None,
        &DosageConfig::default(),
    )
    .unwrap_err();
    assert!(format!("{err_sha:#}").contains("sha256 mismatch"));

    // 7) Table k mismatch rejected with explicit "k mismatch".
    let wrong_k_panel = CorePanel::new(
        31,
        InputType::Reads,
        None,
        gen_distinct_kmers_u64(7002, 31, 4),
    )
    .unwrap();
    let err_k = extract_sample_dosage(
        &reads_tab,
        &wrong_k_panel,
        None,
        None,
        &DosageConfig::default(),
    )
    .unwrap_err();
    assert!(
        format!("{err_k:#}").contains("k mismatch"),
        "expected k mismatch error, got: {err_k:#}"
    );

    // 8) Target panel k mismatch rejected EVEN WITH EMPTY targets map (kills tp.k check!).
    let empty_wrong_k_target = TargetPanel {
        k: 31,
        targets: BTreeMap::new(),
    };
    let err_tk = extract_sample_dosage(
        &reads_tab,
        &panel_reads,
        Some(&empty_wrong_k_target),
        None,
        &DosageConfig::default(),
    )
    .unwrap_err();
    assert!(
        format!("{err_tk:#}").contains("k mismatch"),
        "expected target panel k mismatch error, got: {err_tk:#}"
    );

    // Non-empty target panel with matching k succeeds.
    let mut tm = BTreeMap::new();
    tm.insert("t0".to_string(), target_kmers);
    let ok_tp = TargetPanel { k, targets: tm };
    assert!(
        extract_sample_dosage(
            &reads_tab,
            &panel_reads,
            Some(&ok_tp),
            None,
            &DosageConfig::default(),
        )
        .is_ok()
    );
}

#[test]
fn cutoffs_even_median_partial_core_target_math_and_null_qc() {
    let dir = tempfile::tempdir().unwrap();
    let k = 21;
    let kmers = gen_distinct_kmers_u64(8008, k, 15);
    let core_kmers = kmers[0..10].to_vec();
    let target_kmers = kmers[10..15].to_vec();
    let panel = CorePanel::new(
        k,
        InputType::Reads,
        Some("fold-0".to_string()),
        core_kmers.clone(),
    )
    .unwrap();

    // Exact even-length median + partial core presence + min_count cutoff test:
    // Core k-mers 0..4 have counts [10, 20, 30, 40] -> detected median = 25.0!
    // Core k-mer 4 has count 4 (below min_count=5, so NOT detected!).
    // Core k-mers 5..10 are ABSENT (missing -> count 0, exercises unwrap_or(&0)!).
    // Let's also build an 8-of-10 breadth boundary sample.
    let mut even_items: Vec<(String, u32)> = vec![
        (core_kmers[0].clone(), 10),
        (core_kmers[1].clone(), 20),
        (core_kmers[2].clone(), 30),
        (core_kmers[3].clone(), 40),
        (core_kmers[4].clone(), 4), // below min_count=5!
        // Target k-mers: [0, 0, 50, 100, 200] ->
        //   median_raw_count = 50.0, mean_raw_count = 70.0,
        //   detected_median = 100.0,
        //   normalized_dosage = 50.0 / 25.0 = 2.0,
        //   detected_normalized_dosage = 100.0 / 25.0 = 4.0!
        (target_kmers[2].clone(), 50),
        (target_kmers[3].clone(), 100),
        (target_kmers[4].clone(), 200),
    ];
    let tab_even = dir.path().join("even_median.ktab");
    write_ktab_u64(
        &tab_even,
        "unused_name",
        k,
        InputModality::Reads,
        &even_items,
    );

    let mut targets_map = BTreeMap::new();
    targets_map.insert("plasmid_mixed".to_string(), target_kmers.clone());
    let tp = TargetPanel {
        k,
        targets: targets_map,
    };

    // With min_count = 5 and min_core_breadth = 0.40:
    // 4 of 10 core k-mers have count >= 5 -> breadth = 0.40 (exact boundary!), d_core = 25.0!
    let rep = extract_sample_dosage(
        &tab_even,
        &panel,
        Some(&tp),
        None,
        &DosageConfig {
            min_count: 5,
            min_core_depth: 25.0,
            min_core_breadth: 0.40,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(rep.sample, DEFAULT_SAMPLE_NAME);
    assert_eq!(rep.dosage_status, DosageStatus::Ok);
    assert_eq!(rep.core_detected, 4);
    assert!((rep.core_breadth - 0.40).abs() < 1e-12);
    assert!(
        (rep.d_core - 25.0).abs() < 1e-12,
        "even-length median [10,20,30,40] must be 25.0, got {}",
        rep.d_core
    );
    // All 10 core counts sorted: [0, 0, 0, 0, 0, 4, 10, 20, 30, 40] -> median = (0 + 4) / 2 = 2.0.
    assert!((rep.d_core_all_median - 2.0).abs() < 1e-12);

    // Verify target math distinguishes normalized_dosage (2.0), detected_normalized_dosage (4.0), and mean (70.0).
    let t = &rep.targets[0];
    assert_eq!(t.n_kmers, 5);
    assert_eq!(t.detected_kmers, 3);
    assert!((t.breadth - 0.60).abs() < 1e-12);
    assert!((t.median_raw_count - 50.0).abs() < 1e-12);
    assert!((t.mean_raw_count - 70.0).abs() < 1e-12);
    assert!((t.normalized_dosage.unwrap() - 2.0).abs() < 1e-12);
    assert!((t.detected_normalized_dosage.unwrap() - 4.0).abs() < 1e-12);
    assert!((t.kmers[0].dosage.unwrap() - 0.0).abs() < 1e-12);
    assert!((t.kmers[2].dosage.unwrap() - 2.0).abs() < 1e-12);
    assert!((t.kmers[4].dosage.unwrap() - 8.0).abs() < 1e-12);

    // Omitted QC: EVERY field in JSON `qc` MUST be null (ZERO 0 / 0.0 fields!).
    let rep_json = serde_json::to_value(&rep).unwrap();
    let qc_obj = rep_json["qc"].as_object().unwrap();
    assert!(!qc_obj.is_empty());
    for (field, val) in qc_obj {
        assert!(
            val.is_null(),
            "omitted --qc-json produced non-null {field}={val}"
        );
    }

    // Just above depth threshold (min_core_depth = 25.01 > 25.0) -> InsufficientCoreCoverage, all dosages null!
    let rep_low_depth = extract_sample_dosage(
        &tab_even,
        &panel,
        Some(&tp),
        None,
        &DosageConfig {
            min_count: 5,
            min_core_depth: 25.01,
            min_core_breadth: 0.40,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        rep_low_depth.dosage_status,
        DosageStatus::InsufficientCoreCoverage
    );
    assert!(rep_low_depth.targets[0].normalized_dosage.is_none());
    assert!(
        rep_low_depth.targets[0]
            .detected_normalized_dosage
            .is_none()
    );
    assert!(
        rep_low_depth.targets[0]
            .kmers
            .iter()
            .all(|e| e.dosage.is_none())
    );
    let low_json = serde_json::to_value(&rep_low_depth).unwrap();
    assert!(low_json["targets"][0]["normalized_dosage"].is_null());
    assert!(low_json["targets"][0]["detected_normalized_dosage"].is_null());
    for entry in low_json["targets"][0]["kmers"].as_array().unwrap() {
        assert!(entry["dosage"].is_null());
    }

    // Just above breadth threshold (min_core_breadth = 0.41 > 0.40) -> InsufficientCoreCoverage!
    // Also prove min_count=4 vs min_count=5 changes core_detected from 4 to 5 (breadth 0.50 -> Ok)!
    let rep_low_breadth = extract_sample_dosage(
        &tab_even,
        &panel,
        Some(&tp),
        None,
        &DosageConfig {
            min_count: 5,
            min_core_depth: 10.0,
            min_core_breadth: 0.41,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        rep_low_breadth.dosage_status,
        DosageStatus::InsufficientCoreCoverage
    );

    let rep_mc4 = extract_sample_dosage(
        &tab_even,
        &panel,
        Some(&tp),
        None,
        &DosageConfig {
            min_count: 4,
            min_core_depth: 10.0,
            min_core_breadth: 0.50,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(rep_mc4.core_detected, 5);
    assert_eq!(rep_mc4.dosage_status, DosageStatus::Ok);
    assert!((rep_mc4.d_core - 20.0).abs() < 1e-12); // median of [4, 10, 20, 30, 40]

    // Populated QC JSON: present values emitted, absent values stay null.
    let qc_file = dir.path().join("qc.json");
    fs::write(
        &qc_file,
        r#"{"qc": {"total_reads": 12345, "gc_fraction": 0.51, "singleton_fraction": 0.08}}"#,
    )
    .unwrap();
    let qc: DosageQcVector = load_qc_metrics(&qc_file).unwrap();
    let rep_with_qc = extract_sample_dosage(
        &tab_even,
        &panel,
        None,
        Some(&qc),
        &DosageConfig {
            min_count: 5,
            min_core_depth: 10.0,
            min_core_breadth: 0.40,
            ..Default::default()
        },
    )
    .unwrap();
    let jqc = serde_json::to_value(&rep_with_qc).unwrap();
    assert_eq!(jqc["qc"]["total_reads"], 12_345);
    assert_eq!(jqc["qc"]["gc_fraction"], 0.51);
    assert_eq!(jqc["qc"]["singleton_fraction"], 0.08);
    assert!(jqc["qc"]["usable_read_len_median"].is_null());
    assert!(jqc["qc"]["n50"].is_null());

    // Keep even_items referenced to test 8-of-10 default breadth cutoff cleanly.
    even_items.push((core_kmers[5].clone(), 25));
    even_items.push((core_kmers[6].clone(), 25));
    even_items.push((core_kmers[7].clone(), 25));
    let tab_8of10 = dir.path().join("b80.ktab");
    write_ktab_u64(&tab_8of10, "b80", k, InputModality::Reads, &even_items);
    let rep_80 =
        extract_sample_dosage(&tab_8of10, &panel, None, None, &DosageConfig::default()).unwrap();
    assert_eq!(rep_80.core_detected, 8);
    assert_eq!(rep_80.dosage_status, DosageStatus::Ok);
}

#[test]
fn build_core_panel_excludes_multicopy_plasmids_and_supports_u64_and_u128() {
    let dir = tempfile::tempdir().unwrap();

    // --- k = 21 (u64): 3 training tables ---
    // 10 single-copy core k-mers at count ~30 (all 3 tables)
    // 2 high-copy plasmid k-mers at count 180 in all 3 tables (6x median -> excluded!)
    // 1 plasmid k-mer high-copy (90 = 3x) in ONLY training table 2 (still excluded!)
    // 2 rare k-mers present in only 1 of 3 tables (prevalence 0.33 < 0.90 -> excluded!)
    let k = 21;
    let pool = gen_distinct_kmers_u64(9001, k, 15);
    let core = &pool[0..10];
    let plasmid_all = &pool[10..12];
    let plasmid_one = &pool[12];
    let rare = &pool[13..15];

    let mut train_paths = Vec::new();
    for t_idx in 0..3 {
        let mut items: Vec<(String, u32)> = core
            .iter()
            .enumerate()
            .map(|(j, s)| (s.clone(), 28 + ((j + t_idx) % 5) as u32))
            .collect();
        for p in plasmid_all {
            items.push((p.clone(), 180));
        }
        items.push((plasmid_one.clone(), if t_idx == 2 { 90 } else { 30 }));
        if t_idx == 0 {
            items.push((rare[0].clone(), 30));
        } else if t_idx == 1 {
            items.push((rare[1].clone(), 30));
        }
        let path = dir.path().join(format!("train_{t_idx}.ktab"));
        write_ktab_u64(&path, "sample", k, InputModality::Reads, &items);
        train_paths.push(path);
    }

    let panel = build_training_fold_core_panel(
        &train_paths,
        &CorePanelBuildConfig {
            k,
            input_type: InputType::Reads,
            min_prevalence: 0.90,
            max_copy_factor: 1.5,
            max_kmers: 100,
            fold_id: Some("train-fold-1".to_string()),
        },
    )
    .unwrap();
    panel.validate().unwrap();
    assert_eq!(panel.k, k);
    assert_eq!(panel.input_type, InputType::Reads);
    assert_eq!(panel.fold_id.as_deref(), Some("train-fold-1"));
    assert_eq!(panel.panel_sha256.len(), 64);

    // Exact 10 single-copy core k-mers retained; all 3 plasmid k-mers and both rare k-mers excluded!
    let mut want_core = core.to_vec();
    want_core.sort();
    assert_eq!(panel.kmers, want_core);
    for p in plasmid_all {
        assert!(!panel.kmers.contains(p));
    }
    assert!(!panel.kmers.contains(plasmid_one));
    for r in rare {
        assert!(!panel.kmers.contains(r));
    }

    // max_kmers cap respected deterministically.
    let capped = build_training_fold_core_panel(
        &train_paths,
        &CorePanelBuildConfig {
            k,
            input_type: InputType::Reads,
            min_prevalence: 0.90,
            max_copy_factor: 1.5,
            max_kmers: 4,
            fold_id: None,
        },
    )
    .unwrap();
    assert_eq!(capped.kmers.len(), 4);
    capped.validate().unwrap();

    // Cross-modality rejection in build_training_fold_core_panel.
    assert!(
        build_training_fold_core_panel(
            &train_paths,
            &CorePanelBuildConfig {
                k,
                input_type: InputType::Assembly,
                min_prevalence: 0.90,
                max_copy_factor: 1.5,
                max_kmers: 100,
                fold_id: None,
            },
        )
        .is_err()
    );

    // REVIEWED_COLUMNS + 1 rejected before opening files.
    let too_many: Vec<PathBuf> = (0..=REVIEWED_COLUMNS)
        .map(|i| dir.path().join(format!("nonexistent_{i}.ktab")))
        .collect();
    let err_cap = build_training_fold_core_panel(
        &too_many,
        &CorePanelBuildConfig {
            k,
            input_type: InputType::Reads,
            min_prevalence: 0.90,
            max_copy_factor: 1.5,
            max_kmers: 100,
            fold_id: None,
        },
    )
    .unwrap_err();
    assert!(format!("{err_cap:#}").contains("4096"));

    // --- k = 45 (u128) round trip through build_training_fold_core_panel + extract_sample_dosage ---
    let k45 = 45;
    let pool45 = gen_distinct_kmers_u128(9045, k45, 8);
    let core45 = &pool45[0..6];
    let plasmid45 = &pool45[6..8];
    let w1 = dir.path().join("w1.ktab");
    let w2 = dir.path().join("w2.ktab");
    let mut items_w1: Vec<(String, u32)> = core45.iter().map(|s| (s.clone(), 40)).collect();
    items_w1.push((plasmid45[0].clone(), 200)); // 5x -> excluded
    items_w1.push((plasmid45[1].clone(), 400)); // 10x -> excluded
    write_ktab_u128(&w1, "sample", k45, InputModality::Reads, &items_w1);
    write_ktab_u128(&w2, "sample", k45, InputModality::Reads, &items_w1);

    let panel45 = build_training_fold_core_panel(
        &[w1.clone(), w2],
        &CorePanelBuildConfig {
            k: k45,
            input_type: InputType::Reads,
            min_prevalence: 1.0,
            max_copy_factor: 1.5,
            max_kmers: 50,
            fold_id: Some("k45-fold".to_string()),
        },
    )
    .unwrap();
    assert_eq!(panel45.kmers.len(), 6);
    let mut tp45_map = BTreeMap::new();
    tp45_map.insert("p45".to_string(), vec![plasmid45[1].clone()]);
    let tp45 = TargetPanel {
        k: k45,
        targets: tp45_map,
    };
    let rep45 =
        extract_sample_dosage(&w1, &panel45, Some(&tp45), None, &DosageConfig::default()).unwrap();
    assert_eq!(rep45.dosage_status, DosageStatus::Ok);
    assert!((rep45.d_core - 40.0).abs() < 1e-12);
    assert!((rep45.targets[0].normalized_dosage.unwrap() - 10.0).abs() < 1e-12);
}

#[test]
fn matrix_export_dosage_dense_and_sparse_without_na_bloat() {
    let dir = tempfile::tempdir().unwrap();
    let k = 21;
    let pool = gen_distinct_kmers_u64(4444, k, 12);
    let core = pool[0..10].to_vec();
    let plasmid = pool[10].clone();
    let accessory = pool[11].clone();

    let panel = CorePanel::new(
        k,
        InputType::Reads,
        Some("fold-2".to_string()),
        core.clone(),
    )
    .unwrap();

    // Sample 0: high core depth (d_core = 20.0, breadth = 1.0), plasmid = 100 (5.0x), accessory = 40 (2.0x)
    let mut s0_items: Vec<(String, u32)> = core.iter().map(|s| (s.clone(), 20)).collect();
    s0_items.push((plasmid.clone(), 100));
    s0_items.push((accessory.clone(), 40));
    let s0 = dir.path().join("sensitive_name_0.ktab");
    write_ktab_u64(&s0, "sensitive_hdr_0", k, InputModality::Reads, &s0_items);

    // Sample 1: low core depth (d_core = 3.0 < 10.0), plasmid = 150 -> InsufficientCoreCoverage!
    let mut s1_items: Vec<(String, u32)> = core.iter().map(|s| (s.clone(), 3)).collect();
    s1_items.push((plasmid.clone(), 150));
    let s1 = dir.path().join("sensitive_name_1.ktab");
    write_ktab_u64(&s1, "sensitive_hdr_1", k, InputModality::Reads, &s1_items);

    // 1) Dense dosage layout: column 1 has explicit "NA" in every feature row!
    let dense_dir = dir.path().join("dosage_dense");
    let man_d = export_matrix(
        &[s0.clone(), s1.clone()],
        &dense_dir,
        &MatrixExportConfig {
            mode: MatrixMode::Dosage,
            layout: MatrixLayout::Dense,
            core_panel: Some(panel.clone()),
            dosage: DosageConfig {
                expected_fold_id: Some("fold-2".to_string()),
                ..Default::default()
            },
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(man_d.n_samples, 2);
    assert_eq!(man_d.n_features, 12);
    assert_eq!(man_d.n_missing_samples, 1);
    assert!(!man_d.samples[0].missing_core_coverage);
    assert!(man_d.samples[1].missing_core_coverage);
    assert_eq!(man_d.samples[0].dosage_status, Some(DosageStatus::Ok));
    assert_eq!(
        man_d.samples[1].dosage_status,
        Some(DosageStatus::InsufficientCoreCoverage)
    );

    let dense_mat = fs::read_to_string(dense_dir.join("matrix.tsv")).unwrap();
    let dense_lines: Vec<&str> = dense_mat.lines().collect();
    assert_eq!(dense_lines[0], "feature_idx\ts000001\ts000002");
    assert_eq!(dense_lines.len(), 13);
    for row in &dense_lines[1..] {
        let cols: Vec<&str> = row.split('\t').collect();
        assert_eq!(cols.len(), 3);
        assert_eq!(cols[2], "NA", "dense missing column must write explicit NA");
    }
    // Check plasmid dosage 100 / 20 = 5.000000 in sample 0!
    let feat_text = fs::read_to_string(dense_dir.join("features.tsv")).unwrap();
    let plasmid_row_idx: usize = feat_text
        .lines()
        .skip(1)
        .find_map(|l| {
            let cols: Vec<&str> = l.split('\t').collect();
            (cols[2] == plasmid).then(|| cols[0].parse().unwrap())
        })
        .unwrap();
    let plasmid_cols: Vec<&str> = dense_lines[plasmid_row_idx + 1].split('\t').collect();
    assert_eq!(plasmid_cols[1], "5.000000");
    assert_eq!(plasmid_cols[2], "NA");

    // 2) Sparse dosage layout: ZERO "NA" rows in matrix.tsv, ZERO rows for sample_idx=1,
    //    and missing sample cleanly recorded in missing_mask.tsv + manifest.json!
    let sparse_dir = dir.path().join("dosage_sparse");
    let man_s = export_matrix(
        &[s0.clone(), s1.clone()],
        &sparse_dir,
        &MatrixExportConfig {
            mode: MatrixMode::Dosage,
            layout: MatrixLayout::Sparse,
            core_panel: Some(panel.clone()),
            dosage: DosageConfig::default(),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(man_s.n_features, 12);
    assert_eq!(man_s.n_missing_samples, 1);

    let sparse_mat = fs::read_to_string(sparse_dir.join("matrix.tsv")).unwrap();
    assert!(
        !sparse_mat.contains("NA"),
        "sparse dosage matrix.tsv must not emit NA bloat rows"
    );
    let mut sparse_rows = sparse_mat.lines();
    assert_eq!(
        sparse_rows.next().unwrap(),
        "feature_idx\tsample_idx\tvalue"
    );
    let mut count_rows = 0usize;
    let mut saw_plasmid_5x = false;
    for line in sparse_rows {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 3);
        assert_eq!(
            cols[1], "0",
            "sparse dosage matrix.tsv must omit missing sample_idx=1"
        );
        if cols[0].parse::<usize>().unwrap() == plasmid_row_idx {
            assert_eq!(cols[2], "5.000000");
            saw_plasmid_5x = true;
        }
        count_rows += 1;
    }
    assert_eq!(count_rows, 12); // Only the 12 valid entries of sample 0!
    assert!(saw_plasmid_5x);

    let mask_text = fs::read_to_string(sparse_dir.join("missing_mask.tsv")).unwrap();
    let mask_lines: Vec<&str> = mask_text.lines().collect();
    assert_eq!(
        mask_lines[0],
        "sample_idx\tsample\tmissing_core_coverage\tdosage_status\td_core\tcore_breadth"
    );
    assert!(mask_lines[1].starts_with("0\ts000001\t0\tok\t20.000000\t1.000000"));
    assert!(
        mask_lines[2].starts_with("1\ts000002\t1\tinsufficient_core_coverage\t3.000000\t1.000000")
    );

    // 3) Rejects .kset v1, .grk, wrong modality .ktab, and REVIEWED_COLUMNS + 1 in dosage mode.
    let sp = KmerSpace::new(k);
    let mut v1_keys: Vec<u64> = core
        .iter()
        .map(|s| sp.mix(sp.canonical(s).unwrap()))
        .collect();
    v1_keys.sort_unstable();
    let v1_path = dir.path().join("bad.kset");
    kset::write(&v1_path, "v1", &sp, &v1_keys).unwrap();
    let gr_path = dir.path().join("bad.grk");
    write_gr_kset(&gr_path, "gr", k, InputModality::Reads, &v1_keys).unwrap();
    let asm_path = dir.path().join("bad_asm.ktab");
    write_ktab_u64(&asm_path, "asm", k, InputModality::Assembly, &s0_items);

    for bad_in in [v1_path, gr_path, asm_path] {
        assert!(
            export_matrix(
                &[bad_in],
                &dir.path().join("bad_dosage_out"),
                &MatrixExportConfig {
                    mode: MatrixMode::Dosage,
                    layout: MatrixLayout::Sparse,
                    core_panel: Some(panel.clone()),
                    dosage: DosageConfig::default(),
                    ..Default::default()
                },
            )
            .is_err()
        );
    }

    let too_many: Vec<PathBuf> = (0..=REVIEWED_COLUMNS)
        .map(|i| dir.path().join(format!("no_file_{i}.ktab")))
        .collect();
    let err_cap = export_matrix(
        &too_many,
        &dir.path().join("never_created"),
        &MatrixExportConfig {
            mode: MatrixMode::Dosage,
            layout: MatrixLayout::Sparse,
            core_panel: Some(panel),
            dosage: DosageConfig::default(),
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(format!("{err_cap:#}").contains("4096"));
    assert!(!dir.path().join("never_created").exists());
}

#[test]
fn cli_build_core_panel_dosage_and_matrix_export_dosage_end_to_end() {
    let dir = tempfile::tempdir().unwrap();
    let k = 21;
    let pool = gen_distinct_kmers_u64(5555, k, 12);
    let core = &pool[0..10];
    let plasmid = &pool[10];
    let rare = &pool[11];

    let secret_train = "hidden_train_genome_777";
    let secret_eval = "hidden_eval_isolate_888";
    let secret_low = "hidden_lowcov_isolate_999";

    let mut tr1_items: Vec<(String, u32)> = core.iter().map(|s| (s.clone(), 30)).collect();
    tr1_items.push((plasmid.clone(), 150)); // 5x -> excluded by max_copy_factor=1.5
    tr1_items.push((rare.clone(), 30)); // only in 1 of 2 -> excluded by min_prevalence=0.95
    let tr1 = dir.path().join(format!("{secret_train}_1.ktab"));
    write_ktab_u64(&tr1, secret_train, k, InputModality::Reads, &tr1_items);

    let mut tr2_items: Vec<(String, u32)> = core.iter().map(|s| (s.clone(), 30)).collect();
    tr2_items.push((plasmid.clone(), 150));
    let tr2 = dir.path().join(format!("{secret_train}_2.ktab"));
    write_ktab_u64(&tr2, secret_train, k, InputModality::Reads, &tr2_items);

    let eval_hi = dir.path().join(format!("{secret_eval}.ktab"));
    let mut hi_items: Vec<(String, u32)> = core.iter().map(|s| (s.clone(), 20)).collect();
    hi_items.push((plasmid.clone(), 200)); // 10.0x
    write_ktab_u64(&eval_hi, secret_eval, k, InputModality::Reads, &hi_items);

    let eval_lo = dir.path().join(format!("{secret_low}.ktab"));
    let mut lo_items: Vec<(String, u32)> = core.iter().map(|s| (s.clone(), 2)).collect();
    lo_items.push((plasmid.clone(), 200));
    write_ktab_u64(&eval_lo, secret_low, k, InputModality::Reads, &lo_items);

    let asm_wrong = dir.path().join("wrong_mod.ktab");
    write_ktab_u64(&asm_wrong, "asm", k, InputModality::Assembly, &hi_items);

    // 1) CLI build-core-panel
    let panel_path = dir.path().join("core_panel.json");
    let res = grm_bin()
        .args([
            "kcount",
            "build-core-panel",
            "--input",
            tr1.to_str().unwrap(),
            "--input",
            tr2.to_str().unwrap(),
            "--k",
            "21",
            "--input-type",
            "reads",
            "--min-prevalence",
            "0.95",
            "--max-copy-factor",
            "1.5",
            "--fold-id",
            "cli-fold-0",
            "--out",
            panel_path.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        res.status.success(),
        "build-core-panel failed: {}",
        String::from_utf8_lossy(&res.stderr)
    );
    let panel_text = fs::read_to_string(&panel_path).unwrap();
    assert!(!panel_text.contains(secret_train));
    let loaded_panel = CorePanel::load(&panel_path).unwrap();
    assert_eq!(loaded_panel.kmers.len(), 10);
    assert!(!loaded_panel.kmers.contains(plasmid));
    assert!(!loaded_panel.kmers.contains(rare));

    // Write targets JSON.
    let targets_path = dir.path().join("targets.json");
    let mut tmap = BTreeMap::new();
    tmap.insert("plasmid_10x".to_string(), vec![plasmid.clone()]);
    fs::write(
        &targets_path,
        serde_json::to_vec_pretty(&TargetPanel { k, targets: tmap }).unwrap(),
    )
    .unwrap();

    // 2) CLI dosage without --qc-json and without --name -> default "sample", ALL qc fields null, 0 leaks.
    let dosage_out = dir.path().join("dosage.json");
    let res_d = grm_bin()
        .args([
            "kcount",
            "dosage",
            "--table",
            eval_hi.to_str().unwrap(),
            "--core-panel",
            panel_path.to_str().unwrap(),
            "--targets",
            targets_path.to_str().unwrap(),
            "--expected-fold-id",
            "cli-fold-0",
            "--out",
            dosage_out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        res_d.status.success(),
        "cli dosage failed: {}",
        String::from_utf8_lossy(&res_d.stderr)
    );
    let dosage_text = fs::read_to_string(&dosage_out).unwrap();
    assert!(!dosage_text.contains(secret_eval));
    let dj: serde_json::Value = serde_json::from_str(&dosage_text).unwrap();
    assert_eq!(dj["sample"], DEFAULT_SAMPLE_NAME);
    assert_eq!(dj["dosage_status"], "ok");
    assert_eq!(dj["table_modality"], "reads");
    assert!((dj["d_core"].as_f64().unwrap() - 20.0).abs() < 1e-12);
    assert!((dj["targets"][0]["normalized_dosage"].as_f64().unwrap() - 10.0).abs() < 1e-12);
    for (field, val) in dj["qc"].as_object().unwrap() {
        assert!(val.is_null(), "CLI omitted --qc-json set {field}={val}");
    }

    // 3) CLI dosage rejects Assembly table against Reads panel EVEN WHEN `--input-type reads` is passed!
    let res_bad_mod = grm_bin()
        .args([
            "kcount",
            "dosage",
            "--table",
            asm_wrong.to_str().unwrap(),
            "--core-panel",
            panel_path.to_str().unwrap(),
            "--input-type",
            "reads",
        ])
        .output()
        .unwrap();
    assert!(!res_bad_mod.status.success());
    assert!(String::from_utf8_lossy(&res_bad_mod.stderr).contains("modality mismatch"));

    // 4) CLI matrix-export --mode dosage --layout sparse (0 NA bloat rows, 0 leaks, default s000001/s000002).
    let mat_dir = dir.path().join("cli_dosage_sparse");
    let res_m = grm_bin()
        .args([
            "kcount",
            "matrix-export",
            "--input",
            eval_hi.to_str().unwrap(),
            "--input",
            eval_lo.to_str().unwrap(),
            "--out-dir",
            mat_dir.to_str().unwrap(),
            "--mode",
            "dosage",
            "--layout",
            "sparse",
            "--core-panel",
            panel_path.to_str().unwrap(),
            "--expected-fold-id",
            "cli-fold-0",
        ])
        .output()
        .unwrap();
    assert!(
        res_m.status.success(),
        "cli matrix-export dosage failed: {}",
        String::from_utf8_lossy(&res_m.stderr)
    );
    let m_stdout = String::from_utf8(res_m.stdout).unwrap();
    assert!(!m_stdout.contains(secret_eval));
    assert!(!m_stdout.contains(secret_low));
    let mj: serde_json::Value = serde_json::from_str(&m_stdout).unwrap();
    assert_eq!(mj["n_missing_samples"], 1);
    assert_eq!(mj["samples"][0]["sample"], "s000001");
    assert_eq!(mj["samples"][1]["sample"], "s000002");
    assert_eq!(mj["samples"][1]["missing_core_coverage"], true);

    for fname in [
        "features.tsv",
        "matrix.tsv",
        "missing_mask.tsv",
        "manifest.json",
    ] {
        let text = fs::read_to_string(mat_dir.join(fname)).unwrap();
        assert!(!text.contains(secret_eval), "{fname} leaked {secret_eval}");
        assert!(!text.contains(secret_low), "{fname} leaked {secret_low}");
    }
    let sp_mat = fs::read_to_string(mat_dir.join("matrix.tsv")).unwrap();
    assert!(!sp_mat.contains("NA"));
    assert!(sp_mat.contains("10.000000"));
}
