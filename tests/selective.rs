//! Planted-fault leakage and invariant tests for nested selective calling
//! (`grm::learn::selective`).

use grm::bits::{BitMatrix, Bits};
use grm::intervals::{self, ALPHA};
use grm::layout::Layout;
use grm::learn::eval::{EvalConfig, Role, evaluate, evaluate_assigned};
use grm::learn::fold_masks;
use grm::learn::gbt::GbtParams;
use grm::learn::selective::{
    AbstainReason, COVERAGE_LADDER, CallTally, ConformalThresholds, CoverageThresholds,
    DistanceMatrix, FORMAT, FoldCalibration, JointModelParams, LoadedDistances, POLICIES,
    PRIMARY_COVERAGE, SPARSE_DENOMINATOR_MIN, SelectiveResult, evaluate_split, retention_count,
    select_matched_coverage, write_sketches_bin,
};
use grm::learn::tree::ForestParams;
use grm::species_check::{self, GenomeSketch};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

/// SYNTHETIC cohort: 150 genomes in 30 lineages of 5, a third resistant,
/// informative + noisy features, and a symmetric pairwise Mash distance matrix
/// where within-lineage distances are small (~0.004) and between-lineage
/// distances vary smoothly across [0.02, 0.14].
fn synthetic_cohort() -> (BitMatrix, Bits, Vec<u32>, DistanceMatrix) {
    let n = 150;
    let y = Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0));
    let groups: Vec<u32> = (0..n as u32).map(|i| i / 5).collect();
    let mut rng = ChaCha8Rng::seed_from_u64(20261127);
    let mut x = BitMatrix::with_cols(n);
    for f in 0..10 {
        let present: Vec<usize> = (0..n)
            .filter(|&i| {
                if f < 3 {
                    let flip = if groups[i] >= 24 { 0.30 } else { 0.12 };
                    y.contains(i) != rng.random_bool(flip)
                } else {
                    rng.random_bool(0.4)
                }
            })
            .collect();
        x.push_row(Bits::from_indices(n, present).words());
    }
    let mut values = vec![0.0f32; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = if groups[i] == groups[j] {
                0.002 + 0.0005 * ((i + j) % 5) as f32
            } else {
                let gi = groups[i] as f32;
                let gj = groups[j] as f32;
                let gap = (gi - gj).abs() / 30.0;
                0.015 + 0.12 * gap + 0.001 * ((i * 7 + j * 13) % 11) as f32
            };
            values[i * n + j] = d;
            values[j * n + i] = d;
        }
    }
    let dist = DistanceMatrix::from_dense(n, species_check::K, species_check::SKETCH, values)
        .expect("dense matrix");
    (x, y, groups, dist)
}

fn fast_cfg() -> EvalConfig {
    EvalConfig {
        folds: 3,
        inner_folds: 2,
        forest: ForestParams {
            n_trees: 8,
            ..Default::default()
        },
        gbt: GbtParams {
            rounds: 8,
            ..Default::default()
        },
        screened_features: 8,
        lasso_lambdas: 5,
        ..Default::default()
    }
}

fn assert_opt_f64_bits_eq(a: Option<f64>, b: Option<f64>, ctx: &str) {
    match (a, b) {
        (None, None) => {}
        (Some(x), Some(y)) => assert_eq!(
            x.to_bits(),
            y.to_bits(),
            "{ctx}: float bits differ ({x} vs {y})"
        ),
        _ => panic!("{ctx}: option mismatch ({a:?} vs {b:?})"),
    }
}

fn assert_calibration_bit_identical(a: &FoldCalibration, b: &FoldCalibration, ctx: &str) {
    assert_eq!(a.fold, b.fold, "{ctx}: fold");
    assert_eq!(
        (a.train_n, a.train_r, a.train_s, a.inner_null_fits),
        (b.train_n, b.train_r, b.train_s, b.inner_null_fits),
        "{ctx}: train counts"
    );
    assert_eq!(
        a.inner_unavailable_reason, b.inner_unavailable_reason,
        "{ctx}: inner reason"
    );
    assert_eq!(
        a.thresholds.len(),
        b.thresholds.len(),
        "{ctx}: thresholds len"
    );
    for (ta, tb) in a.thresholds.iter().zip(&b.thresholds) {
        let CoverageThresholds {
            target_coverage: ca,
            tau_blind: ba,
            t_lo: la,
            t_hi: ha,
            tau_margin: ma,
            tau_distance: da,
            tau_joint: ja,
        } = ta;
        let CoverageThresholds {
            target_coverage: cb,
            tau_blind: bb,
            t_lo: lb,
            t_hi: hb,
            tau_margin: mb,
            tau_distance: db,
            tau_joint: jb,
        } = tb;
        assert_eq!(ca.to_bits(), cb.to_bits(), "{ctx}: target_coverage");
        assert_eq!(ba.to_bits(), bb.to_bits(), "{ctx} c={ca}: tau_blind");
        assert_eq!(la.to_bits(), lb.to_bits(), "{ctx} c={ca}: t_lo");
        assert_eq!(ha.to_bits(), hb.to_bits(), "{ctx} c={ca}: t_hi");
        assert_eq!(ma.to_bits(), mb.to_bits(), "{ctx} c={ca}: tau_margin");
        assert_opt_f64_bits_eq(*da, *db, &format!("{ctx} c={ca}: tau_distance"));
        assert_opt_f64_bits_eq(*ja, *jb, &format!("{ctx} c={ca}: tau_joint"));
    }
    assert_eq!(a.conformal.len(), b.conformal.len(), "{ctx}: conformal len");
    for (ca, cb) in a.conformal.iter().zip(&b.conformal) {
        let ConformalThresholds {
            alpha: aa,
            q_r: qra,
            q_s: qsa,
            t_lo: loa,
            t_hi: hia,
        } = ca;
        let ConformalThresholds {
            alpha: ab,
            q_r: qrb,
            q_s: qsb,
            t_lo: lob,
            t_hi: hib,
        } = cb;
        assert_eq!(aa.to_bits(), ab.to_bits(), "{ctx}: alpha");
        assert_opt_f64_bits_eq(*qra, *qrb, &format!("{ctx} alpha={aa}: q_r"));
        assert_opt_f64_bits_eq(*qsa, *qsb, &format!("{ctx} alpha={aa}: q_s"));
        assert_opt_f64_bits_eq(*loa, *lob, &format!("{ctx} alpha={aa}: conformal t_lo"));
        assert_opt_f64_bits_eq(*hia, *hib, &format!("{ctx} alpha={aa}: conformal t_hi"));
    }
    match (&a.joint_model, &b.joint_model) {
        (None, None) => {}
        (Some(ja), Some(jb)) => {
            let JointModelParams {
                intercept: i1,
                weight_margin: wm1,
                weight_distance: wd1,
                mean_margin: mm1,
                sd_margin: sm1,
                mean_distance: md1,
                sd_distance: sd1,
                l2_lambda: l1,
                fallback_equal_weights: fb1,
            } = ja;
            let JointModelParams {
                intercept: i2,
                weight_margin: wm2,
                weight_distance: wd2,
                mean_margin: mm2,
                sd_margin: sm2,
                mean_distance: md2,
                sd_distance: sd2,
                l2_lambda: l2,
                fallback_equal_weights: fb2,
            } = jb;
            assert_eq!(i1.to_bits(), i2.to_bits(), "{ctx}: joint intercept");
            assert_eq!(wm1.to_bits(), wm2.to_bits(), "{ctx}: joint weight_margin");
            assert_eq!(wd1.to_bits(), wd2.to_bits(), "{ctx}: joint weight_distance");
            assert_eq!(mm1.to_bits(), mm2.to_bits(), "{ctx}: joint mean_margin");
            assert_eq!(sm1.to_bits(), sm2.to_bits(), "{ctx}: joint sd_margin");
            assert_eq!(md1.to_bits(), md2.to_bits(), "{ctx}: joint mean_distance");
            assert_eq!(sd1.to_bits(), sd2.to_bits(), "{ctx}: joint sd_distance");
            assert_eq!(l1.to_bits(), l2.to_bits(), "{ctx}: joint l2_lambda");
            assert_eq!(fb1, fb2, "{ctx}: joint fallback_equal_weights");
        }
        _ => panic!("{ctx}: joint_model presence differs"),
    }
}

/// CRITICAL TEST 1:
/// Mutate outer test labels `y[E_f]` and outer test distances `D[E_f, *]` (and
/// `D[*, E_f]`) and assert every inner-calibrated threshold (`t_lo`, `t_hi`,
/// `tau_d`, `tau_joint`, `tau_blind`, Mondrian conformal quantiles) and joint
/// risk-model parameter is bit-for-bit identical across all learners on drawn
/// (`random`, `lineage`) and supplied (`study-whole`) splits.
#[test]
fn critical_leakage_mutating_outer_test_labels_and_distances_preserves_calibration_bit_for_bit() {
    let (x, y, groups, dist) = synthetic_cohort();
    let cfg = fast_cfg();
    let out = evaluate("synth__drug", &x, &y, &groups, &[], 3, None, &cfg).unwrap();
    let loaded = LoadedDistances::from_matrix(dist.clone(), "synthetic");
    let n = x.cols;
    let all: Vec<usize> = (0..n).collect();

    // 1A: Drawn splits (`random` and `lineage`).
    for split in &out.splits {
        assert!(split.unavailable_reason.is_none());
        let trains: Vec<Bits> = fold_masks(&split.test_fold, split.folds)
            .into_iter()
            .map(|(tr, _)| tr)
            .collect();
        let baseline = evaluate_split(
            "synth__drug",
            "synth",
            "drug",
            &split.scheme,
            &x,
            &y,
            &groups,
            &[],
            &trains,
            &all,
            &split.test_fold,
            &split.learners,
            None,
            None,
            3,
            &cfg,
            &loaded,
        )
        .unwrap();

        // For each outer fold f, mutate ONLY outer test fold E_f (labels AND distances AND features).
        for f in 0..split.folds {
            let mut y_mut = y.clone();
            let mut d_mut = dist.clone();
            let mut x_mut = BitMatrix::with_cols(n);
            let in_test = |i: usize| split.test_fold[i] == f;
            for r in 0..x.rows {
                let row = (0..n).filter(|&c| x.get(r, c) != in_test(c));
                x_mut.push_row(Bits::from_indices(n, row).words());
            }
            for u in 0..n {
                if in_test(u) {
                    if y.contains(u) {
                        y_mut = y_mut.and_not(&Bits::from_indices(n, [u]));
                    } else {
                        y_mut.insert(u);
                    }
                    // Overwrite every distance involving u with extreme adversarial values.
                    for j in 0..n {
                        if j != u {
                            d_mut.set(u, j, 0.999);
                            d_mut.set(j, u, 0.999);
                        }
                    }
                }
            }
            let loaded_mut = LoadedDistances::from_matrix(d_mut, "synthetic");
            let mutated = evaluate_split(
                "synth__drug",
                "synth",
                "drug",
                &split.scheme,
                &x_mut,
                &y_mut,
                &groups,
                &[],
                &trains,
                &all,
                &split.test_fold,
                &split.learners,
                None,
                None,
                3,
                &cfg,
                &loaded_mut,
            )
            .unwrap();

            for (la, lb) in baseline.learners.iter().zip(&mutated.learners) {
                assert_eq!(la.learner, lb.learner);
                assert_calibration_bit_identical(
                    &la.calibrations[f],
                    &lb.calibrations[f],
                    &format!("scheme={} learner={} fold={f}", split.scheme, la.learner),
                );
            }
        }
    }

    // 1B: Supplied fold-file split (`study-whole`) with `Role::Unused` samples too.
    let roles: Vec<Role> = groups
        .iter()
        .map(|&g| {
            if g >= 27 {
                Role::Unused
            } else {
                Role::Test((g % 3) as usize)
            }
        })
        .collect();
    let assigned = evaluate_assigned(
        "synth__drug",
        "study-whole",
        &x,
        &y,
        &groups,
        &[],
        &roles,
        3,
        &cfg,
    )
    .unwrap();
    let trains: Vec<Bits> = (0..3)
        .map(|f| {
            Bits::from_indices(
                n,
                (0..n).filter(|&i| {
                    matches!(roles[i], Role::Test(g) if g != f) || roles[i] == Role::Train
                }),
            )
        })
        .collect();
    let base_study = evaluate_split(
        "synth__drug",
        "synth",
        "drug",
        "study-whole",
        &x,
        &y,
        &groups,
        &[],
        &trains,
        &assigned.tested,
        &assigned.result.test_fold,
        &assigned.result.learners,
        None,
        Some("0".repeat(64)),
        assigned.min_support,
        &cfg,
        &loaded,
    )
    .unwrap();

    // Mutate fold 0's test samples AND all Unused samples in labels, features, and distances.
    let mut y_mut = y.clone();
    let mut d_mut = dist.clone();
    let mut x_mut = BitMatrix::with_cols(n);
    let outside_t0 = |i: usize| !trains[0].contains(i);
    for r in 0..x.rows {
        let row = (0..n).filter(|&c| x.get(r, c) != outside_t0(c));
        x_mut.push_row(Bits::from_indices(n, row).words());
    }
    for u in 0..n {
        if outside_t0(u) {
            if y.contains(u) {
                y_mut = y_mut.and_not(&Bits::from_indices(n, [u]));
            } else {
                y_mut.insert(u);
            }
            for j in 0..n {
                if j != u {
                    d_mut.set(u, j, 0.888);
                    d_mut.set(j, u, 0.888);
                }
            }
        }
    }
    let loaded_mut = LoadedDistances::from_matrix(d_mut, "synthetic");
    let mut_study = evaluate_split(
        "synth__drug",
        "synth",
        "drug",
        "study-whole",
        &x_mut,
        &y_mut,
        &groups,
        &[],
        &trains,
        &assigned.tested,
        &assigned.result.test_fold,
        &assigned.result.learners,
        None,
        Some("0".repeat(64)),
        assigned.min_support,
        &cfg,
        &loaded_mut,
    )
    .unwrap();
    for (la, lb) in base_study.learners.iter().zip(&mut_study.learners) {
        assert_calibration_bit_identical(
            &la.calibrations[0],
            &lb.calibrations[0],
            &format!("study-whole learner={} fold=0", la.learner),
        );
    }

    // Sanity check (negative control): mutating training distances inside T_0 MUST move tau_distance.
    let mut d_train_mut = dist;
    for i in trains[0].iter_ones() {
        for j in trains[0].iter_ones() {
            if i != j {
                d_train_mut.set(i, j, (d_train_mut.get(i, j) * 2.0).min(0.95));
            }
        }
    }
    let loaded_train_mut = LoadedDistances::from_matrix(d_train_mut, "synthetic");
    let train_mut_study = evaluate_split(
        "synth__drug",
        "synth",
        "drug",
        "study-whole",
        &x,
        &y,
        &groups,
        &[],
        &trains,
        &assigned.tested,
        &assigned.result.test_fold,
        &assigned.result.learners,
        None,
        Some("0".repeat(64)),
        assigned.min_support,
        &cfg,
        &loaded_train_mut,
    )
    .unwrap();
    let d0_before = base_study.learners[0].calibrations[0].thresholds[3]
        .tau_distance
        .unwrap();
    let d0_after = train_mut_study.learners[0].calibrations[0].thresholds[3]
        .tau_distance
        .unwrap();
    assert_ne!(
        d0_before.to_bits(),
        d0_after.to_bits(),
        "mutating training distances should change calibrated tau_distance"
    );
}

/// CRITICAL TEST 2:
/// Verify exact `floor(c * n)` retention count on every outer fold in
/// `matched_coverage` mode across the full `COVERAGE_LADDER`, under tied scores
/// (including binary `scm` scores and completely constant tied scores), with
/// zero dependence on outer test labels.
#[test]
fn critical_matched_coverage_retains_exact_floor_c_n_under_ties_and_label_flips() {
    let (x, y, groups, dist) = synthetic_cohort();
    let cfg = fast_cfg();
    let out = evaluate("synth__drug", &x, &y, &groups, &[], 3, None, &cfg).unwrap();
    let loaded = LoadedDistances::from_matrix(dist, "synthetic");
    let n = x.cols;
    let all: Vec<usize> = (0..n).collect();

    for split in &out.splits {
        let trains: Vec<Bits> = fold_masks(&split.test_fold, split.folds)
            .into_iter()
            .map(|(tr, _)| tr)
            .collect();
        let res = evaluate_split(
            "synth__drug",
            "synth",
            "drug",
            &split.scheme,
            &x,
            &y,
            &groups,
            &[],
            &trains,
            &all,
            &split.test_fold,
            &split.learners,
            None,
            None,
            3,
            &cfg,
            &loaded,
        )
        .unwrap();

        // Check every learner (including binary-output `scm`), every policy, every coverage rung.
        for learner in &res.learners {
            for policy in &learner.policies {
                assert!(
                    policy.available,
                    "{} {} {} unavailable: {:?}",
                    split.scheme, learner.learner, policy.policy, policy.unavailable_reason
                );
                assert_eq!(policy.matched_coverage.len(), COVERAGE_LADDER.len());
                for (ci, pt) in policy.matched_coverage.iter().enumerate() {
                    let c = COVERAGE_LADDER[ci];
                    assert_eq!(pt.target_coverage, c);
                    assert_eq!(pt.folds.len(), split.folds);
                    let mut sum_called = 0usize;
                    let mut sum_n = 0usize;
                    for fold in &pt.folds {
                        let want = retention_count(c, fold.n);
                        assert_eq!(
                            fold.called,
                            want,
                            "scheme={} learner={} policy={} c={c} fold={}: got {} != floor({c} * {}) = {want}",
                            split.scheme,
                            learner.learner,
                            policy.policy,
                            fold.fold,
                            fold.called,
                            fold.n
                        );
                        assert_eq!(fold.called, fold.called_r + fold.called_s);
                        assert_eq!(fold.called_r, fold.tp + fold.fn_);
                        assert_eq!(fold.called_s, fold.tn + fold.fp);
                        assert_eq!(
                            fold.abstain_r + fold.abstain_s,
                            fold.abstain_reasons.total()
                        );
                        sum_called += fold.called;
                        sum_n += fold.n;
                    }
                    assert_eq!(pt.accounting.called, sum_called);
                    assert_eq!(pt.accounting.n, sum_n);
                }
            }
        }
    }

    // Direct tie & label-independence stress test on awkward fold sizes (including n=140 at c=0.95).
    for fold_n in [17usize, 23, 37, 140] {
        // Every score completely tied at 0.5 (binary SCM-style tie).
        let tied: Vec<(usize, f64)> = (0..fold_n).map(|pos| (pos, 0.5)).collect();
        for &c in &COVERAGE_LADDER {
            let kept_a = select_matched_coverage(&tied, c, 42, "lineage", 1);
            let kept_b = select_matched_coverage(&tied, c, 42, "lineage", 1);
            assert_eq!(kept_a, kept_b, "tie-breaking must be deterministic");
            let count = kept_a.iter().filter(|&&k| k).count();
            assert_eq!(
                count,
                retention_count(c, fold_n),
                "n={fold_n} c={c}: exact floor count failed"
            );
        }
    }

    // Label-flip invariance on outer fold 0 in evaluate_split:
    // Flipping every truth label in E_0 keeps the exact same samples selected in
    // matched_coverage, so called_r and called_s swap identically!
    let lineage = out.splits.iter().find(|s| s.scheme == "lineage").unwrap();
    let trains: Vec<Bits> = fold_masks(&lineage.test_fold, lineage.folds)
        .into_iter()
        .map(|(tr, _)| tr)
        .collect();
    let mut y_flip = y.clone();
    for i in 0..n {
        if lineage.test_fold[i] == 0 {
            if y.contains(i) {
                y_flip = y_flip.and_not(&Bits::from_indices(n, [i]));
            } else {
                y_flip.insert(i);
            }
        }
    }
    let res_orig = evaluate_split(
        "synth__drug",
        "synth",
        "drug",
        "lineage",
        &x,
        &y,
        &groups,
        &[],
        &trains,
        &all,
        &lineage.test_fold,
        &lineage.learners,
        None,
        None,
        3,
        &cfg,
        &loaded,
    )
    .unwrap();
    let res_flip = evaluate_split(
        "synth__drug",
        "synth",
        "drug",
        "lineage",
        &x,
        &y_flip,
        &groups,
        &[],
        &trains,
        &all,
        &lineage.test_fold,
        &lineage.learners,
        None,
        None,
        3,
        &cfg,
        &loaded,
    )
    .unwrap();
    for (la, lb) in res_orig.learners.iter().zip(&res_flip.learners) {
        for (pa, pb) in la.policies.iter().zip(&lb.policies) {
            for (ca, cb) in pa.matched_coverage.iter().zip(&pb.matched_coverage) {
                let fa = &ca.folds[0];
                let fb = &cb.folds[0];
                assert_eq!(fa.called, fb.called);
                assert_eq!(
                    fa.called_r, fb.called_s,
                    "{} {} c={}: flipping y[E_0] changed retained sample set",
                    la.learner, pa.policy, ca.target_coverage
                );
                assert_eq!(
                    fa.called_s, fb.called_r,
                    "{} {} c={}: flipping y[E_0] changed retained sample set",
                    la.learner, pa.policy, ca.target_coverage
                );
            }
        }
    }
}

/// CRITICAL TEST 3:
/// Verify `coverage_r`, `coverage_s`, `vme_called`, `vme_all`, `me_called`,
/// `me_all`, `npv`, `for_s`, `ppv`, `balanced_accuracy_called`, per-100
/// clinical accounting, sparse-denominator threshold (25), Clopper-Pearson
/// 95% CIs, and `conformal_both` vs `conformal_empty` math against exact
/// hand-computed reference cases.
#[test]
fn critical_accounting_and_conformal_math_matches_hand_computed_reference() {
    // Reference case 1: non-sparse counts (n_r = 40, n_s = 50, called_r = 32, called_s = 40).
    let mut tally = CallTally::default();
    for _ in 0..28 {
        tally.record_call(true, true); // TP = 28
    }
    for _ in 0..4 {
        tally.record_call(true, false); // FN = 4
    }
    for _ in 0..8 {
        tally.record_abstain(true, AbstainReason::UncertainScore); // abstain_r = 8
    }
    for _ in 0..35 {
        tally.record_call(false, false); // TN = 35
    }
    for _ in 0..5 {
        tally.record_call(false, true); // FP = 5
    }
    for _ in 0..6 {
        tally.record_abstain(false, AbstainReason::UncertainScore);
    }
    for _ in 0..4 {
        tally.record_abstain(false, AbstainReason::OodMashDistance);
    }
    let acc = tally.finish(false);

    assert_eq!((acc.n, acc.n_r, acc.n_s), (90, 40, 50));
    assert_eq!((acc.called, acc.called_r, acc.called_s), (72, 32, 40));
    assert_eq!((acc.abstain_r, acc.abstain_s), (8, 10));
    assert_eq!((acc.tp, acc.fp, acc.tn, acc.fn_), (28, 5, 35, 4));
    assert_eq!(acc.abstain_reasons.uncertain_score, 14);
    assert_eq!(acc.abstain_reasons.ood_mash_distance, 4);
    assert_eq!(acc.abstain_reasons.total(), 18);
    assert!(!acc.sparse_denominator);
    assert!(!acc.sparse_denominator_r);
    assert!(!acc.sparse_denominator_s);
    assert!(acc.set_true_label_coverage.is_none());

    assert!((acc.coverage.unwrap() - 72.0 / 90.0).abs() < 1e-12);
    assert!((acc.coverage_r.unwrap() - 32.0 / 40.0).abs() < 1e-12);
    assert!((acc.coverage_s.unwrap() - 40.0 / 50.0).abs() < 1e-12);
    assert!((acc.vme_called.rate.unwrap() - 4.0 / 32.0).abs() < 1e-12);
    assert!((acc.vme_all.rate.unwrap() - 4.0 / 40.0).abs() < 1e-12);
    assert!((acc.me_called.rate.unwrap() - 5.0 / 40.0).abs() < 1e-12);
    assert!((acc.me_all.rate.unwrap() - 5.0 / 50.0).abs() < 1e-12);
    assert!((acc.npv.rate.unwrap() - 35.0 / 39.0).abs() < 1e-12);
    assert!((acc.for_s.rate.unwrap() - 4.0 / 39.0).abs() < 1e-12);
    assert!((acc.npv.rate.unwrap() + acc.for_s.rate.unwrap() - 1.0).abs() < 1e-12);
    assert!((acc.ppv.rate.unwrap() - 28.0 / 33.0).abs() < 1e-12);
    assert!((acc.balanced_accuracy_called.unwrap() - 0.875).abs() < 1e-12);

    // Exact 95% Clopper-Pearson interval checks against crate::intervals.
    let (v_lo, v_hi) = intervals::clopper_pearson(4, 32, ALPHA).unwrap();
    assert!((acc.vme_called.ci_lo.unwrap() - v_lo).abs() < 1e-12);
    assert!((acc.vme_called.ci_hi.unwrap() - v_hi).abs() < 1e-12);
    let (n_lo, n_hi) = intervals::clopper_pearson(35, 39, ALPHA).unwrap();
    assert!((acc.npv.ci_lo.unwrap() - n_lo).abs() < 1e-12);
    assert!((acc.npv.ci_hi.unwrap() - n_hi).abs() < 1e-12);
    assert!((acc.for_s.ci_lo.unwrap() - (1.0 - n_hi)).abs() < 1e-12);
    assert!((acc.for_s.ci_hi.unwrap() - (1.0 - n_lo)).abs() < 1e-12);

    // Per-100 clinical accounting sums to 100.0 for each class.
    assert!((acc.per_100_r.caught.unwrap() - 70.0).abs() < 1e-12);
    assert!((acc.per_100_r.missed.unwrap() - 10.0).abs() < 1e-12);
    assert!((acc.per_100_r.deferred.unwrap() - 20.0).abs() < 1e-12);
    assert!((acc.per_100_s.cleared.unwrap() - 70.0).abs() < 1e-12);
    assert!((acc.per_100_s.overcalled.unwrap() - 10.0).abs() < 1e-12);
    assert!((acc.per_100_s.deferred.unwrap() - 20.0).abs() < 1e-12);

    // Reference case 2: exact sparse-denominator boundary (24 vs 25).
    assert_eq!(SPARSE_DENOMINATOR_MIN, 25);
    let mut sparse = CallTally::default();
    for _ in 0..24 {
        sparse.record_call(true, true);
    }
    for _ in 0..25 {
        sparse.record_call(false, false);
    }
    let acc_sparse = sparse.finish(false);
    assert!(acc_sparse.sparse_denominator_r);
    assert!(!acc_sparse.sparse_denominator_s);
    assert!(acc_sparse.sparse_denominator);
    assert!(acc_sparse.vme_called.sparse_denominator);
    assert!(!acc_sparse.me_called.sparse_denominator);

    sparse.record_call(true, true); // now called_r = 25, called_s = 25
    let acc_exact25 = sparse.finish(false);
    assert!(!acc_exact25.sparse_denominator_r);
    assert!(!acc_exact25.sparse_denominator_s);
    assert!(!acc_exact25.sparse_denominator);
    assert!(!acc_exact25.vme_called.sparse_denominator);

    // Reference case 3: Conformal mixed batch containing {R}, {S}, {R, S} (both) and {} (empty).
    let mut conf = CallTally::default();
    conf.record_conformal(true, true, false); // {R}, true=R -> TP, covered
    conf.record_conformal(true, false, true); // {S}, true=R -> FN, NOT covered
    conf.record_conformal(true, true, true); // {R,S}, true=R -> conformal_both, covered
    conf.record_conformal(true, false, false); // {}, true=R -> conformal_empty, NOT covered
    conf.record_conformal(false, false, true); // {S}, true=S -> TN, covered
    conf.record_conformal(false, true, false); // {R}, true=S -> FP, NOT covered
    conf.record_conformal(false, true, true); // {R,S}, true=S -> conformal_both, covered
    conf.record_conformal(false, false, false); // {}, true=S -> conformal_empty, NOT covered
    let acc_c = conf.finish(true);
    assert_eq!((acc_c.tp, acc_c.fp, acc_c.tn, acc_c.fn_), (1, 1, 1, 1));
    assert_eq!(acc_c.abstain_reasons.conformal_both, 2);
    assert_eq!(acc_c.abstain_reasons.conformal_empty, 2);
    assert_eq!(acc_c.abstain_reasons.total(), 4);
    assert_eq!((acc_c.abstain_r, acc_c.abstain_s), (2, 2));
    // Covered = TP(1) + TN(1) + conformal_both(2) = 4 out of 8 (conformal_empty does NOT cover).
    let set_cov = acc_c.set_true_label_coverage.as_ref().unwrap();
    assert_eq!((set_cov.x, set_cov.n), (4, 8));
    assert!((set_cov.rate.unwrap() - 0.5).abs() < 1e-12);
    assert!((acc_c.coverage_r.unwrap() - 0.5).abs() < 1e-12);
    assert!((acc_c.coverage_s.unwrap() - 0.5).abs() < 1e-12);
    assert!((acc_c.vme_called.rate.unwrap() - 0.5).abs() < 1e-12);
    assert!((acc_c.vme_all.rate.unwrap() - 0.25).abs() < 1e-12);
    assert!((acc_c.npv.rate.unwrap() - 0.5).abs() < 1e-12);
    assert!((acc_c.for_s.rate.unwrap() - 0.5).abs() < 1e-12);
}

/// Verify missing sketches gracefully fall back to score-only policies while
/// marking distance policies unavailable, and emitted `grm-selective-v1` JSON
/// contains zero genome IDs under the `scripts/c7-outside/scan.py` token rule.
#[test]
fn missing_sketches_fallback_and_zero_id_hygiene() {
    let (x, y, groups, _) = synthetic_cohort();
    let cfg = fast_cfg();
    let out = evaluate(
        "escherichia_coli__meropenem",
        &x,
        &y,
        &groups,
        &[],
        3,
        None,
        &cfg,
    )
    .unwrap();
    let lineage = out.splits.iter().find(|s| s.scheme == "lineage").unwrap();
    let trains: Vec<Bits> = fold_masks(&lineage.test_fold, lineage.folds)
        .into_iter()
        .map(|(tr, _)| tr)
        .collect();
    let all: Vec<usize> = (0..x.cols).collect();

    // Realistic BV-BRC-shaped synthetic genome IDs that must NEVER appear in JSON.
    let fake_ids: Vec<String> = (0..x.cols).map(|i| format!("562.{}", 90000 + i)).collect();
    let dir = tempfile::tempdir().unwrap();
    let layout = Layout::new(dir.path());
    let loaded_missing = grm::learn::selective::load_cohort_distances(
        &layout,
        "escherichia_coli",
        &fake_ids,
        &Bits::ones(x.cols),
    )
    .unwrap();
    assert!(!loaded_missing.status.available);

    let res = evaluate_split(
        "escherichia_coli__meropenem",
        "escherichia_coli",
        "meropenem",
        "lineage",
        &x,
        &y,
        &groups,
        &[],
        &trains,
        &all,
        &lineage.test_fold,
        &lineage.learners,
        None,
        None,
        3,
        &cfg,
        &loaded_missing,
    )
    .unwrap();

    assert_eq!(res.format, FORMAT);
    assert_eq!(res.evaluation_revision, grm::analysis::EVALUATION_REVISION);
    assert_eq!(res.constants.primary_coverage, PRIMARY_COVERAGE);
    assert_eq!(res.constants.sparse_denominator_min, SPARSE_DENOMINATOR_MIN);
    assert!(!res.distances.available);

    for lr in &res.learners {
        let names: Vec<&str> = lr.policies.iter().map(|p| p.policy.as_str()).collect();
        assert_eq!(names, POLICIES);
        assert!(lr.policies[0].available); // outcome_blind
        assert!(lr.policies[1].available); // score_margin
        assert!(!lr.policies[2].available); // mash_distance
        assert!(!lr.policies[3].available); // margin_plus_distance
        assert!(lr.policies[2].unavailable_reason.is_some());
        assert!(lr.policies[3].unavailable_reason.is_some());
        assert!(lr.mondrian_conformal_available);
    }

    let json = serde_json::to_string(&res).unwrap();
    // Round-trip check
    let parsed: SelectiveResult = serde_json::from_str(&json).unwrap();
    assert_eq!(parsed.format, FORMAT);

    // Zero-ID check matching scripts/c7-outside/scan.py
    for id in &fake_ids {
        assert!(
            !json.contains(id),
            "genome ID {id} leaked into selective JSON"
        );
    }

    // Also verify cached species-check JSON sketches AND sketches.bin sidecar work cleanly.
    let bin_entries: Vec<(String, Vec<u64>)> = fake_ids
        .iter()
        .enumerate()
        .map(|(i, id)| {
            let base = (groups[i] as u64) * 500 + (i as u64 % 5) * 5;
            let hashes: Vec<u64> = (0..20).map(|k| base + k * 3 + 1).collect();
            (id.clone(), hashes)
        })
        .collect();
    // Put half in sketches.bin and half in species-check/sketches/<id>.json.
    let bin_path = layout.species_dir("escherichia_coli").join("sketches.bin");
    write_sketches_bin(
        &bin_path,
        species_check::K,
        species_check::SKETCH,
        &bin_entries[..75],
    )
    .unwrap();
    for (id, hashes) in &bin_entries[75..] {
        let p = species_check::sketch_path(&layout, id);
        std::fs::create_dir_all(p.parent().unwrap()).unwrap();
        let gs = GenomeSketch {
            genome_id: id.clone(),
            k: species_check::K,
            sketch_size: species_check::SKETCH,
            length: 500_000,
            gc: 0.51,
            sketch: hashes.clone(),
        };
        std::fs::write(&p, serde_json::to_vec(&gs).unwrap()).unwrap();
    }
    let loaded_ok = grm::learn::selective::load_cohort_distances(
        &layout,
        "escherichia_coli",
        &fake_ids,
        &Bits::ones(x.cols),
    )
    .unwrap();
    assert!(loaded_ok.status.available);
    let res_ok = evaluate_split(
        "escherichia_coli__meropenem",
        "escherichia_coli",
        "meropenem",
        "lineage",
        &x,
        &y,
        &groups,
        &[],
        &trains,
        &all,
        &lineage.test_fold,
        &lineage.learners,
        None,
        None,
        3,
        &cfg,
        &loaded_ok,
    )
    .unwrap();
    let json_ok = serde_json::to_string(&res_ok).unwrap();
    for lr in &res_ok.learners {
        assert!(lr.policies.iter().all(|p| p.available));
    }
    for id in &fake_ids {
        assert!(
            !json_ok.contains(id),
            "genome ID {id} leaked into distance-enabled selective JSON"
        );
    }
}
