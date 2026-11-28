//! Integration and CLI tests for streaming k-mer set algebra (`grm kcount op`)
//! and bounded matrix export (`grm kcount matrix-export`).

use grm::dna::KmerSpace;
use grm::kcount::format::{
    DEFAULT_SAMPLE_NAME, GrKsetReader, InputModality, KtabReader, MAX_COUNT, hex_bytes,
    write_gr_kset, write_ktab,
};
use grm::kcount::key::{KmerSpace64, KmerSpace128};
use grm::kcount::ops::{
    MatrixExportConfig, MatrixLayout, MatrixMode, SetOp, SetOpConfig, export_matrix, run_set_op,
};
use grm::kset;
use grm::matrix::REVIEWED_COLUMNS;
use ring::digest::{Context as Sha256Context, SHA256};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn grm_bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_grm"))
}

fn sha256_file(path: &Path) -> String {
    let bytes = fs::read(path).unwrap();
    let mut ctx = Sha256Context::new(&SHA256);
    ctx.update(&bytes);
    hex_bytes(ctx.finish().as_ref())
}

fn build_ktab_u64(
    path: &Path,
    name: &str,
    k: usize,
    modality: InputModality,
    items: &[(&str, u32)],
) {
    let sp = KmerSpace64::new(k);
    let mut entries: Vec<(u64, u32)> = items
        .iter()
        .map(|&(s, c)| (sp.mix(sp.canonical(s).unwrap()), c))
        .collect();
    entries.sort_unstable_by_key(|&(key, _)| key);
    let sat = entries.iter().filter(|&&(_, c)| c == MAX_COUNT).count() as u64;
    write_ktab(path, name, k, modality, 1, sat, &entries).unwrap();
}

fn build_grk_u64(path: &Path, name: &str, k: usize, modality: InputModality, kmers: &[&str]) {
    let sp = KmerSpace64::new(k);
    let mut keys: Vec<u64> = kmers
        .iter()
        .map(|&s| sp.mix(sp.canonical(s).unwrap()))
        .collect();
    keys.sort_unstable();
    write_gr_kset(path, name, k, modality, &keys).unwrap();
}

fn build_v1_u64(path: &Path, name: &str, k: usize, kmers: &[&str]) {
    let sp = KmerSpace::new(k);
    let mut keys: Vec<u64> = kmers
        .iter()
        .map(|&s| sp.mix(sp.canonical(s).unwrap()))
        .collect();
    keys.sort_unstable();
    kset::write(path, name, &sp, &keys).unwrap();
}

fn build_ktab_u128(
    path: &Path,
    name: &str,
    k: usize,
    modality: InputModality,
    items: &[(&str, u32)],
) {
    let sp = KmerSpace128::new(k);
    let mut entries: Vec<(u128, u32)> = items
        .iter()
        .map(|&(s, c)| (sp.mix(sp.canonical(s).unwrap()), c))
        .collect();
    entries.sort_unstable_by_key(|&(key, _)| key);
    let sat = entries.iter().filter(|&&(_, c)| c == MAX_COUNT).count() as u64;
    write_ktab(path, name, k, modality, 1, sat, &entries).unwrap();
}

fn build_grk_u128(path: &Path, name: &str, k: usize, modality: InputModality, kmers: &[&str]) {
    let sp = KmerSpace128::new(k);
    let mut keys: Vec<u128> = kmers
        .iter()
        .map(|&s| sp.mix(sp.canonical(s).unwrap()))
        .collect();
    keys.sort_unstable();
    write_gr_kset(path, name, k, modality, &keys).unwrap();
}

#[test]
fn set_ops_exact_counts_and_saturation_for_u64_and_u128() {
    let dir = tempfile::tempdir().unwrap();

    // --- k = 21 (u64 path), mixing .ktab, .kset v1, and .grk across 2 and 3 tables ---
    let k21 = 21;
    let sp21 = KmerSpace64::new(k21);
    let a = "AAAAAAAAAAAAAAAAAAAAA";
    let b = "CCCCCCCCCCCCCCCCCCCCC";
    let c = "ACGTACGTACGTACGTACGTA";
    let d = "TGCATGCATGCATGCATGCAT";
    let e = "GCTAGCTAGCTAGCTAGCTAG";

    let t1 = dir.path().join("t1.ktab");
    let t2 = dir.path().join("t2.ktab");
    let t3_v1 = dir.path().join("t3.kset");
    let t3_gr = dir.path().join("t3.grk");

    build_ktab_u64(
        &t1,
        "s1",
        k21,
        InputModality::Reads,
        &[(a, 5), (b, 10), (c, 4), (e, 9)],
    );
    build_ktab_u64(
        &t2,
        "s2",
        k21,
        InputModality::Reads,
        &[(a, 7), (b, 3), (d, 11)],
    );
    build_v1_u64(&t3_v1, "s3", k21, &[b, c, d]); // presence-only (count = 1 each)
    build_grk_u64(&t3_gr, "s3", k21, InputModality::Reads, &[b, c, d]);

    let key21 = |s: &str| sp21.mix(sp21.canonical(s).unwrap());

    // 2-table union: exact sums 5 + 7 = 12, 10 + 3 = 13, 4, 11, 9; modality Reads.
    let u2_path = dir.path().join("u2.ktab");
    let u2 = run_set_op(
        SetOp::Union,
        &[t1.clone(), t2.clone()],
        &u2_path,
        &SetOpConfig::default(),
    )
    .unwrap();
    assert_eq!(u2.n_kmers_out, 5);
    assert_eq!(u2.max_count_out, 13);
    assert_eq!(u2.saturated_counts, 0);
    assert!(!u2.counts_saturated);
    assert_eq!(u2.modality, InputModality::Reads);
    let u2_map: HashMap<u64, u32> = KtabReader::open_verified(&u2_path)
        .unwrap()
        .read_all::<u64>()
        .unwrap()
        .into_iter()
        .collect();
    assert_eq!(u2_map[&key21(a)], 12);
    assert_eq!(u2_map[&key21(b)], 13);
    assert_eq!(u2_map[&key21(c)], 4);
    assert_eq!(u2_map[&key21(d)], 11);
    assert_eq!(u2_map[&key21(e)], 9);

    // 2-table intersect: min counts (a -> min(5,7)=5, b -> min(10,3)=3).
    let i2_path = dir.path().join("i2.ktab");
    let i2 = run_set_op(
        SetOp::Intersect,
        &[t1.clone(), t2.clone()],
        &i2_path,
        &SetOpConfig::default(),
    )
    .unwrap();
    assert_eq!(i2.n_kmers_out, 2);
    assert_eq!(i2.max_count_out, 5);
    let i2_map: HashMap<u64, u32> = KtabReader::open_verified(&i2_path)
        .unwrap()
        .read_all::<u64>()
        .unwrap()
        .into_iter()
        .collect();
    assert_eq!(i2_map[&key21(a)], 5);
    assert_eq!(i2_map[&key21(b)], 3);

    // 2-table diff (t1 \ t2): c (4), e (9).
    let d2_path = dir.path().join("d2.ktab");
    let d2 = run_set_op(
        SetOp::Diff,
        &[t1.clone(), t2.clone()],
        &d2_path,
        &SetOpConfig::default(),
    )
    .unwrap();
    assert_eq!(d2.n_kmers_out, 2);
    let d2_map: HashMap<u64, u32> = KtabReader::open_verified(&d2_path)
        .unwrap()
        .read_all::<u64>()
        .unwrap()
        .into_iter()
        .collect();
    assert_eq!(d2_map[&key21(c)], 4);
    assert_eq!(d2_map[&key21(e)], 9);

    // 2-table symdiff: c (4), d (11), e (9).
    let s2_path = dir.path().join("s2.ktab");
    let s2 = run_set_op(
        SetOp::SymDiff,
        &[t1.clone(), t2.clone()],
        &s2_path,
        &SetOpConfig::default(),
    )
    .unwrap();
    assert_eq!(s2.n_kmers_out, 3);
    let s2_map: HashMap<u64, u32> = KtabReader::open_verified(&s2_path)
        .unwrap()
        .read_all::<u64>()
        .unwrap()
        .into_iter()
        .collect();
    assert_eq!(s2_map[&key21(c)], 4);
    assert_eq!(s2_map[&key21(d)], 11);
    assert_eq!(s2_map[&key21(e)], 9);

    // 3-table mixed formats (.ktab + .kset v1 + .grk):
    // t1={a:5,b:10,c:4,e:9}, t3_v1={b:1,c:1,d:1}, t2={a:7,b:3,d:11}
    // Wait: mix t1(.ktab), t3_v1(.kset v1), t3_gr(.grk) AND t1, t2, t3_v1!
    let inputs3 = vec![t1.clone(), t2.clone(), t3_v1.clone()];
    // Union: a=12, b=14, c=5, d=12, e=9; modality Unknown because t3_v1 has no modality.
    let u3_path = dir.path().join("u3.ktab");
    let u3 = run_set_op(SetOp::Union, &inputs3, &u3_path, &SetOpConfig::default()).unwrap();
    assert_eq!(u3.modality, InputModality::Unknown);
    let u3_map: HashMap<u64, u32> = KtabReader::open_verified(&u3_path)
        .unwrap()
        .read_all::<u64>()
        .unwrap()
        .into_iter()
        .collect();
    assert_eq!(u3_map[&key21(a)], 12);
    assert_eq!(u3_map[&key21(b)], 14);
    assert_eq!(u3_map[&key21(c)], 5);
    assert_eq!(u3_map[&key21(d)], 12);
    assert_eq!(u3_map[&key21(e)], 9);

    // Intersect of 3: only b (min(10, 3, 1) = 1).
    let i3_path = dir.path().join("i3.ktab");
    let i3 = run_set_op(
        SetOp::Intersect,
        &inputs3,
        &i3_path,
        &SetOpConfig::default(),
    )
    .unwrap();
    assert_eq!(i3.n_kmers_out, 1);
    assert_eq!(i3.max_count_out, 1);
    assert!(i3.presence_only);
    let i3_entries = KtabReader::open_verified(&i3_path)
        .unwrap()
        .read_all::<u64>()
        .unwrap();
    assert_eq!(i3_entries, vec![(key21(b), 1)]);

    // Diff of 3 (t1 \ (t2 U t3_v1)): only e (9).
    let d3_path = dir.path().join("d3.ktab");
    let d3 = run_set_op(SetOp::Diff, &inputs3, &d3_path, &SetOpConfig::default()).unwrap();
    assert_eq!(d3.n_kmers_out, 1);
    let d3_entries = KtabReader::open_verified(&d3_path)
        .unwrap()
        .read_all::<u64>()
        .unwrap();
    assert_eq!(d3_entries, vec![(key21(e), 9)]);

    // SymDiff of 3 (odd occurrence counts): b (in 3 -> sum 14), e (in 1 -> sum 9).
    let s3_path = dir.path().join("s3.ktab");
    let s3 = run_set_op(SetOp::SymDiff, &inputs3, &s3_path, &SetOpConfig::default()).unwrap();
    assert_eq!(s3.n_kmers_out, 2);
    let s3_map: HashMap<u64, u32> = KtabReader::open_verified(&s3_path)
        .unwrap()
        .read_all::<u64>()
        .unwrap()
        .into_iter()
        .collect();
    assert_eq!(s3_map[&key21(b)], 14);
    assert_eq!(s3_map[&key21(e)], 9);

    // Mixing t1 (.ktab, Reads), t2 (.ktab, Reads), and t3_gr (.grk, Reads) preserves Reads.
    let u3_reads_path = dir.path().join("u3_reads.grk");
    let u3_reads = run_set_op(
        SetOp::Union,
        &[t1, t2, t3_gr],
        &u3_reads_path,
        &SetOpConfig {
            presence_only_out: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(u3_reads.modality, InputModality::Reads);
    assert_eq!(u3_reads.format, "GRMKGR1");
    assert_eq!(
        GrKsetReader::open_verified(&u3_reads_path)
            .unwrap()
            .modality,
        InputModality::Reads
    );

    // --- k = 45 (u128 path), mixing .ktab and .grk across 2 and 3 tables ---
    let k45 = 45;
    let sp45 = KmerSpace128::new(k45);
    let w1 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA";
    let w2 = "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT";
    let w3 = "GGATCCGGATCCGGATCCGGATCCGGATCCGGATCCGGATCCGGA";
    let w4 = "AATTCCAATTCCAATTCCAATTCCAATTCCAATTCCAATTCCAAT";
    let key45 = |s: &str| sp45.mix(sp45.canonical(s).unwrap());

    let wt1 = dir.path().join("w1.ktab");
    let wt2 = dir.path().join("w2.ktab");
    let wt3 = dir.path().join("w3.grk");
    build_ktab_u128(
        &wt1,
        "w1",
        k45,
        InputModality::Assembly,
        &[(w1, 6), (w2, 2_500_000_000), (w3, 8)],
    );
    build_ktab_u128(
        &wt2,
        "w2",
        k45,
        InputModality::Assembly,
        &[(w1, 9), (w2, 2_500_000_000), (w4, 3)],
    );
    build_grk_u128(&wt3, "w3", k45, InputModality::Assembly, &[w1, w2, w4]);

    // 3-table symdiff overflows u32 on w2 (2.5e9 + 2.5e9 + 1 > u32::MAX): refused by default.
    let w_sym_fail = dir.path().join("w_sym_fail.ktab");
    assert!(
        run_set_op(
            SetOp::SymDiff,
            &[wt1.clone(), wt2.clone(), wt3.clone()],
            &w_sym_fail,
            &SetOpConfig::default(),
        )
        .is_err()
    );
    // Allowed with allow_saturation: w1 (6+9+1=16), w2 (clamped u32::MAX), w3 (8); w4 is in 2 (even -> excluded).
    let w_sym_ok = dir.path().join("w_sym_ok.ktab");
    let ws = run_set_op(
        SetOp::SymDiff,
        &[wt1.clone(), wt2.clone(), wt3.clone()],
        &w_sym_ok,
        &SetOpConfig {
            allow_saturation: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(ws.n_kmers_out, 3);
    assert_eq!(ws.saturated_counts, 1);
    assert!(ws.counts_saturated);
    assert_eq!(ws.modality, InputModality::Assembly);
    let ws_map: HashMap<u128, u32> = KtabReader::open_verified(&w_sym_ok)
        .unwrap()
        .read_all::<u128>()
        .unwrap()
        .into_iter()
        .collect();
    assert_eq!(ws_map[&key45(w1)], 16);
    assert_eq!(ws_map[&key45(w2)], u32::MAX);
    assert_eq!(ws_map[&key45(w3)], 8);

    // Intersect at k = 45 across 3 tables: w1 (min(6,9,1)=1), w2 (min(2.5e9,2.5e9,1)=1).
    let wi_path = dir.path().join("w_inter.ktab");
    let wi = run_set_op(
        SetOp::Intersect,
        &[wt1.clone(), wt2.clone(), wt3.clone()],
        &wi_path,
        &SetOpConfig::default(),
    )
    .unwrap();
    assert_eq!(wi.n_kmers_out, 2);
    assert_eq!(wi.max_count_out, 1);

    // Diff at k = 45: only w3 (8).
    let wd_path = dir.path().join("w_diff.ktab");
    let wd = run_set_op(
        SetOp::Diff,
        &[wt1, wt2, wt3],
        &wd_path,
        &SetOpConfig::default(),
    )
    .unwrap();
    assert_eq!(wd.n_kmers_out, 1);
    let wd_entries = KtabReader::open_verified(&wd_path)
        .unwrap()
        .read_all::<u128>()
        .unwrap();
    assert_eq!(wd_entries, vec![(key45(w3), 8)]);
}

#[test]
fn unequal_k_rejected_across_all_ops_and_matrix_export() {
    let dir = tempfile::tempdir().unwrap();
    let k21_path = dir.path().join("k21.ktab");
    let k31_path = dir.path().join("k31.ktab");
    build_ktab_u64(
        &k21_path,
        "sample",
        21,
        InputModality::Assembly,
        &[("ACGTACGTACGTACGTACGTA", 2)],
    );
    build_ktab_u64(
        &k31_path,
        "sample",
        31,
        InputModality::Assembly,
        &[("ACGTACGTACGTACGTACGTACGTACGTACG", 2)],
    );

    let pair = vec![k21_path, k31_path];
    for op in [SetOp::Union, SetOp::Intersect, SetOp::Diff, SetOp::SymDiff] {
        let err = run_set_op(
            op,
            &pair,
            &dir.path().join(format!("{op:?}.ktab")),
            &SetOpConfig::default(),
        )
        .unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("k mismatch"),
            "op {op:?} did not report k mismatch: {msg}"
        );
    }

    let err = export_matrix(
        &pair,
        &dir.path().join("mat"),
        &MatrixExportConfig::default(),
    )
    .unwrap_err();
    assert!(
        format!("{err:#}").contains("k mismatch"),
        "matrix-export did not report k mismatch: {err:#}"
    );
}

#[test]
fn matrix_export_exact_cells_sparse_and_dense_prevalence_and_cap() {
    let dir = tempfile::tempdir().unwrap();
    let k = 21;
    let sp = KmerSpace64::new(k);
    let k1 = "AAAAAAAAAAAAAAAAAAAAA"; // support 1: [5, 0, 0]
    let k2 = "CCCCCCCCCCCCCCCCCCCCC"; // support 2: [2, 8, 0]
    let k3 = "ACGTACGTACGTACGTACGTA"; // support 3: [4, 1, 9]
    let k4 = "TGCATGCATGCATGCATGCAT"; // support 2: [0, 3, 6]

    let t0 = dir.path().join("secret_genome_alpha.ktab");
    let t1 = dir.path().join("secret_genome_beta.ktab");
    let t2 = dir.path().join("secret_genome_gamma.ktab");
    build_ktab_u64(
        &t0,
        "secret_header_a",
        k,
        InputModality::Reads,
        &[(k1, 5), (k2, 2), (k3, 4)],
    );
    build_ktab_u64(
        &t1,
        "secret_header_b",
        k,
        InputModality::Reads,
        &[(k2, 8), (k3, 1), (k4, 3)],
    );
    build_ktab_u64(
        &t2,
        "secret_header_c",
        k,
        InputModality::Reads,
        &[(k3, 9), (k4, 6)],
    );

    let inputs = vec![t0, t1, t2];

    // Expected counts indexed by canonical decoded k-mer sequence.
    let expected_counts: HashMap<String, [u32; 3]> = [
        (sp.decode(sp.canonical(k1).unwrap()), [5, 0, 0]),
        (sp.decode(sp.canonical(k2).unwrap()), [2, 8, 0]),
        (sp.decode(sp.canonical(k3).unwrap()), [4, 1, 9]),
        (sp.decode(sp.canonical(k4).unwrap()), [0, 3, 6]),
    ]
    .into_iter()
    .collect();

    // 1) Dense PA (all 4 features), default sample names s000001..s000003.
    let dense_pa_dir = dir.path().join("dense_pa");
    let man = export_matrix(
        &inputs,
        &dense_pa_dir,
        &MatrixExportConfig {
            mode: MatrixMode::Pa,
            layout: MatrixLayout::Dense,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(man.n_samples, 3);
    assert_eq!(man.n_features, 4);
    assert_eq!(man.n_missing_samples, 0);
    assert_eq!(
        man.samples
            .iter()
            .map(|s| s.sample.as_str())
            .collect::<Vec<_>>(),
        vec!["s000001", "s000002", "s000003"]
    );
    assert_eq!(
        man.sha256_features,
        sha256_file(&dense_pa_dir.join("features.tsv"))
    );
    assert_eq!(
        man.sha256_matrix,
        sha256_file(&dense_pa_dir.join("matrix.tsv"))
    );
    assert_eq!(
        man.sha256_missing_mask,
        sha256_file(&dense_pa_dir.join("missing_mask.tsv"))
    );

    // Check zero leaks of filenames or embedded .ktab header names across ALL outputs.
    for fname in [
        "features.tsv",
        "matrix.tsv",
        "missing_mask.tsv",
        "manifest.json",
    ] {
        let content = fs::read_to_string(dense_pa_dir.join(fname)).unwrap();
        for forbidden in [
            "secret_genome_alpha",
            "secret_genome_beta",
            "secret_genome_gamma",
            "secret_header_a",
            "secret_header_b",
            "secret_header_c",
        ] {
            assert!(
                !content.contains(forbidden),
                "{fname} leaked forbidden token {forbidden}"
            );
        }
    }

    let feat_text = fs::read_to_string(dense_pa_dir.join("features.tsv")).unwrap();
    let feat_rows: Vec<Vec<&str>> = feat_text
        .lines()
        .skip(1)
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(feat_rows.len(), 4);

    let mat_text = fs::read_to_string(dense_pa_dir.join("matrix.tsv")).unwrap();
    let mut mat_lines = mat_text.lines();
    assert_eq!(
        mat_lines.next().unwrap(),
        "feature_idx\ts000001\ts000002\ts000003"
    );
    for (row_idx, line) in mat_lines.enumerate() {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 4);
        assert_eq!(cols[0], row_idx.to_string());
        let seq = feat_rows[row_idx][2];
        let want_cnt = expected_counts[seq];
        for s in 0..3 {
            let want_pa = if want_cnt[s] > 0 { "1" } else { "0" };
            assert_eq!(
                cols[s + 1],
                want_pa,
                "dense PA mismatch on {seq} sample {s}"
            );
        }
    }

    // 2) Dense Counts with custom sample names.
    let dense_cnt_dir = dir.path().join("dense_cnt");
    let man_cnt = export_matrix(
        &inputs,
        &dense_cnt_dir,
        &MatrixExportConfig {
            mode: MatrixMode::Counts,
            layout: MatrixLayout::Dense,
            sample_names: Some(vec!["gA".into(), "gB".into(), "gC".into()]),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(man_cnt.n_features, 4);
    let cnt_text = fs::read_to_string(dense_cnt_dir.join("matrix.tsv")).unwrap();
    let mut cnt_lines = cnt_text.lines();
    assert_eq!(cnt_lines.next().unwrap(), "feature_idx\tgA\tgB\tgC");
    for (row_idx, line) in cnt_lines.enumerate() {
        let cols: Vec<&str> = line.split('\t').collect();
        let seq = feat_rows[row_idx][2];
        let want_cnt = expected_counts[seq];
        for s in 0..3 {
            assert_eq!(
                cols[s + 1],
                want_cnt[s].to_string(),
                "dense counts mismatch on {seq} sample {s}"
            );
        }
    }

    // 3) Sparse PA: exact non-zero entries ONLY (zero absent entries emitted!).
    let sparse_pa_dir = dir.path().join("sparse_pa");
    export_matrix(
        &inputs,
        &sparse_pa_dir,
        &MatrixExportConfig {
            mode: MatrixMode::Pa,
            layout: MatrixLayout::Sparse,
            ..Default::default()
        },
    )
    .unwrap();
    let spa_text = fs::read_to_string(sparse_pa_dir.join("matrix.tsv")).unwrap();
    let mut spa_lines = spa_text.lines();
    assert_eq!(spa_lines.next().unwrap(), "feature_idx\tsample_idx\tvalue");
    let mut observed_pa: HashMap<(usize, usize), &str> = HashMap::new();
    for line in spa_lines {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 3);
        let fi: usize = cols[0].parse().unwrap();
        let si: usize = cols[1].parse().unwrap();
        assert_eq!(cols[2], "1");
        assert!(observed_pa.insert((fi, si), cols[2]).is_none());
    }
    // Total present pairs across k1(1), k2(2), k3(3), k4(2) = 8 exact entries!
    assert_eq!(observed_pa.len(), 8);
    for (fi, frow) in feat_rows.iter().enumerate() {
        let want_cnt = expected_counts[frow[2]];
        for (si, &c) in want_cnt.iter().enumerate() {
            assert_eq!(
                observed_pa.contains_key(&(fi, si)),
                c > 0,
                "sparse PA presence mismatch at ({fi}, {si})"
            );
        }
    }

    // 4) Sparse Counts with --min-prevalence 0.65 (keeps support >= 2 of 3 -> k2, k3, k4).
    let sparse_cnt_dir = dir.path().join("sparse_cnt_prev");
    let man_sp = export_matrix(
        &inputs,
        &sparse_cnt_dir,
        &MatrixExportConfig {
            mode: MatrixMode::Counts,
            layout: MatrixLayout::Sparse,
            min_prevalence: Some(0.65),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(man_sp.n_features, 3);
    let sp_feats_text = fs::read_to_string(sparse_cnt_dir.join("features.tsv")).unwrap();
    let sp_feats: Vec<Vec<&str>> = sp_feats_text
        .lines()
        .skip(1)
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(sp_feats.len(), 3);
    for frow in &sp_feats {
        assert_ne!(frow[2], sp.decode(sp.canonical(k1).unwrap()));
    }
    let sc_text = fs::read_to_string(sparse_cnt_dir.join("matrix.tsv")).unwrap();
    let mut sc_observed: HashMap<(usize, usize), u32> = HashMap::new();
    for line in sc_text.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        let fi: usize = cols[0].parse().unwrap();
        let si: usize = cols[1].parse().unwrap();
        let val: u32 = cols[2].parse().unwrap();
        assert!(val > 0, "sparse counts must never emit 0");
        sc_observed.insert((fi, si), val);
    }
    // k2(2) + k3(3) + k4(2) = 7 exact non-zero entries.
    assert_eq!(sc_observed.len(), 7);
    for (fi, frow) in sp_feats.iter().enumerate() {
        let want_cnt = expected_counts[frow[2]];
        for (si, &c) in want_cnt.iter().enumerate() {
            if c > 0 {
                assert_eq!(sc_observed[&(fi, si)], c);
            } else {
                assert!(!sc_observed.contains_key(&(fi, si)));
            }
        }
    }

    // 5) REVIEWED_COLUMNS + 1 inputs rejected without opening files.
    let too_many: Vec<PathBuf> = (0..=REVIEWED_COLUMNS)
        .map(|i| dir.path().join(format!("missing_path_{i}.ktab")))
        .collect();
    assert!(
        export_matrix(
            &too_many,
            &dir.path().join("never"),
            &MatrixExportConfig::default()
        )
        .is_err()
    );
    assert!(!dir.path().join("never").exists());
}

#[test]
fn cli_op_and_matrix_export_end_to_end() {
    let dir = tempfile::tempdir().unwrap();
    let k = 21;
    let k1 = "AAAAAAAAAAAAAAAAAAAAA";
    let k2 = "CCCCCCCCCCCCCCCCCCCCC";
    let k3 = "ACGTACGTACGTACGTACGTA";

    let secret_a = "private_isolate_path_11111";
    let secret_b = "private_isolate_path_22222";
    let in_a = dir.path().join(format!("{secret_a}.ktab"));
    let in_b = dir.path().join(format!("{secret_b}.ktab"));
    build_ktab_u64(
        &in_a,
        secret_a,
        k,
        InputModality::Reads,
        &[(k1, 5), (k2, 3_000_000_000)],
    );
    build_ktab_u64(
        &in_b,
        secret_b,
        k,
        InputModality::Reads,
        &[(k1, 7), (k2, 2_000_000_000), (k3, 4)],
    );

    // 1) CLI union without --allow-saturation fails non-zero.
    let u_fail = dir.path().join("u_fail.ktab");
    let res = grm_bin()
        .args([
            "kcount",
            "op",
            "union",
            "--input",
            in_a.to_str().unwrap(),
            "--input",
            in_b.to_str().unwrap(),
            "--out",
            u_fail.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!res.status.success());
    assert!(!u_fail.exists());

    // 2) CLI union with --allow-saturation succeeds, defaults sample name to "sample",
    //    records saturation, and leaks 0 input filename/path strings.
    let u_ok = dir.path().join("u_ok.ktab");
    let res = grm_bin()
        .args([
            "kcount",
            "op",
            "union",
            "--input",
            in_a.to_str().unwrap(),
            "--input",
            in_b.to_str().unwrap(),
            "--out",
            u_ok.to_str().unwrap(),
            "--allow-saturation",
        ])
        .output()
        .unwrap();
    assert!(
        res.status.success(),
        "cli union failed: {}",
        String::from_utf8_lossy(&res.stderr)
    );
    let stdout = String::from_utf8(res.stdout).unwrap();
    assert!(!stdout.contains(secret_a));
    assert!(!stdout.contains(secret_b));
    let j: serde_json::Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(j["sample"], DEFAULT_SAMPLE_NAME);
    assert_eq!(j["modality"], "reads");
    assert_eq!(j["n_kmers_out"], 3);
    assert_eq!(j["max_count_out"], u32::MAX);
    assert_eq!(j["saturated_counts"], 1);
    assert_eq!(j["counts_saturated"], true);
    let r_u = KtabReader::open_verified(&u_ok).unwrap();
    assert_eq!(r_u.name, DEFAULT_SAMPLE_NAME);
    assert_eq!(r_u.saturated_counts, 1);

    // 3) Extension mismatch refused by CLI.
    let res_ext = grm_bin()
        .args([
            "kcount",
            "op",
            "intersect",
            "--input",
            in_a.to_str().unwrap(),
            "--input",
            in_b.to_str().unwrap(),
            "--out",
            dir.path().join("wrong.grk").to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(!res_ext.status.success());
    assert!(String::from_utf8_lossy(&res_ext.stderr).contains("disagrees"));

    // 4) CLI intersect --presence-only, diff, and symdiff write verified outputs with default sample "sample".
    let i_gr = dir.path().join("inter.grk");
    let res_i = grm_bin()
        .args([
            "kcount",
            "op",
            "intersect",
            "--input",
            in_a.to_str().unwrap(),
            "--input",
            in_b.to_str().unwrap(),
            "--out",
            i_gr.to_str().unwrap(),
            "--presence-only",
        ])
        .output()
        .unwrap();
    assert!(res_i.status.success());
    assert_eq!(
        GrKsetReader::open_verified(&i_gr).unwrap().name,
        DEFAULT_SAMPLE_NAME
    );

    let d_tab = dir.path().join("diff.ktab");
    let res_d = grm_bin()
        .args([
            "kcount",
            "op",
            "diff",
            "--input",
            in_b.to_str().unwrap(),
            "--input",
            in_a.to_str().unwrap(),
            "--out",
            d_tab.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(res_d.status.success());
    let dj: serde_json::Value = serde_json::from_slice(&res_d.stdout).unwrap();
    assert_eq!(dj["n_kmers_out"], 1);
    assert_eq!(dj["max_count_out"], 4);
    assert_eq!(dj["sample"], DEFAULT_SAMPLE_NAME);

    let sd_tab = dir.path().join("symdiff.ktab");
    let res_sd = grm_bin()
        .args([
            "kcount",
            "op",
            "symdiff",
            "--input",
            in_a.to_str().unwrap(),
            "--input",
            in_b.to_str().unwrap(),
            "--out",
            sd_tab.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(res_sd.status.success());
    let sdj: serde_json::Value = serde_json::from_slice(&res_sd.stdout).unwrap();
    assert_eq!(sdj["n_kmers_out"], 1);
    assert_eq!(sdj["max_count_out"], 4);
    assert_eq!(sdj["sample"], DEFAULT_SAMPLE_NAME);

    // 5) CLI matrix-export default naming ("s000001", "s000002"), --min-prevalence, 0 leaks.
    let mat_dir = dir.path().join("cli_mat");
    let res_m = grm_bin()
        .args([
            "kcount",
            "matrix-export",
            "--input",
            in_a.to_str().unwrap(),
            "--input",
            in_b.to_str().unwrap(),
            "--out-dir",
            mat_dir.to_str().unwrap(),
            "--mode",
            "counts",
            "--layout",
            "dense",
            "--min-prevalence",
            "1.0",
        ])
        .output()
        .unwrap();
    assert!(
        res_m.status.success(),
        "cli matrix-export failed: {}",
        String::from_utf8_lossy(&res_m.stderr)
    );
    let m_stdout = String::from_utf8(res_m.stdout).unwrap();
    assert!(!m_stdout.contains(secret_a));
    assert!(!m_stdout.contains(secret_b));
    let mj: serde_json::Value = serde_json::from_str(&m_stdout).unwrap();
    assert_eq!(mj["n_features"], 2); // k1 and k2 have prevalence 1.0; k3 is 0.5
    assert_eq!(mj["samples"][0]["sample"], "s000001");
    assert_eq!(mj["samples"][1]["sample"], "s000002");

    for fname in [
        "features.tsv",
        "matrix.tsv",
        "missing_mask.tsv",
        "manifest.json",
    ] {
        let text = fs::read_to_string(mat_dir.join(fname)).unwrap();
        assert!(!text.contains(secret_a), "{fname} leaked {secret_a}");
        assert!(!text.contains(secret_b), "{fname} leaked {secret_b}");
    }
    let mat_text = fs::read_to_string(mat_dir.join("matrix.tsv")).unwrap();
    assert!(mat_text.starts_with("feature_idx\ts000001\ts000002\n"));
}
