//! Learners on a supplied fold plan.
//!
//! `evaluate` draws its own outer folds, inner folds and tuning holdouts.
//! `evaluate_plan` takes all three from a plan instead (G's `FoldPlan`, or a
//! `CohortPlan` for an arm whose cohort differs from the reference arm's) and
//! draws none. It scores one split, fitting every learner once per outer
//! fold, on a matrix the caller builds for that fold's training samples, so
//! an ablation arm can assemble its feature blocks per fold (contract §16,
//! §18). Each learner fits exactly as in `evaluate`: the forest, boosting and
//! the lasso screen `screened_features` on their own fitting samples.
//!
//! With G's repeat-0 plan of `EvalConfig::default()` and the matrix
//! `training_candidates` gives, a split's result equals `evaluate`'s (test
//! `the_default_plan_reproduces_evaluate`). Final models, the MIC regression
//! and conformal diagnostics are not part of an arm's score and are not fitted.

use super::eval::{self, EvalConfig, LEARNERS, LassoContext, LearnerResult, Scored, SplitResult};
use super::gbt::{self, GbtParams, Objective};
use super::{lasso, screen};
use crate::ablate::plan::{self, CohortPlan, FoldPlan, Scheme, TuningHoldouts};
use crate::bits::{BitMatrix, Bits};
use anyhow::{Context, Result, ensure};
use std::time::Instant;

/// The lasso path's stopping size, `evaluate`'s literal.
const LASSO_MAX_NONZERO: usize = 400;

/// The partitions one split is scored on.
#[derive(Clone, Debug)]
pub struct PlanView<'a> {
    pub scheme: Scheme,
    pub k: usize,
    pub inner_k: usize,
    /// Seeds every outer fold's random draws inside the learners, through
    /// `plan::fold_seed`.
    pub repeat_seed: u64,
    /// Each sample's outer test fold, or None for a sample that is never
    /// scored.
    pub outer: Vec<Option<usize>>,
    /// Per outer fold, samples that neither train nor are scored there (a
    /// cohort plan's `held_out`); empty for a `FoldPlan`.
    pub held_out: &'a [Vec<usize>],
    pub inner: &'a [Vec<Option<usize>>],
    pub holdouts: &'a [TuningHoldouts],
    /// Why the plan cannot be scored, or None.
    pub unavailable_reason: Option<String>,
    /// Seeds the learners' own draws (the forest's and boosting's) with
    /// `plan::fold_seed(learner_seed, f)` instead of the plan's repeat seed,
    /// leaving every partition as it is; None in every plan. Only the
    /// calibration twin of the input comparison sets it.
    pub learner_seed: Option<u64>,
}

impl<'a> From<&'a FoldPlan> for PlanView<'a> {
    fn from(p: &'a FoldPlan) -> Self {
        PlanView {
            scheme: p.scheme,
            k: p.k,
            inner_k: p.inner_k,
            repeat_seed: p.repeat_seed,
            outer: p.outer.iter().map(|&f| Some(f)).collect(),
            held_out: &[],
            inner: &p.inner,
            holdouts: &p.holdouts,
            unavailable_reason: p.unavailable_reason.clone(),
            learner_seed: None,
        }
    }
}

impl<'a> From<&'a CohortPlan> for PlanView<'a> {
    fn from(p: &'a CohortPlan) -> Self {
        PlanView {
            scheme: p.scheme,
            k: p.k,
            inner_k: p.inner_k,
            repeat_seed: p.repeat_seed,
            outer: p.outer.clone(),
            held_out: &p.held_out,
            inner: &p.inner,
            holdouts: &p.holdouts,
            unavailable_reason: p.unavailable_reason.clone(),
            learner_seed: None,
        }
    }
}

impl PlanView<'_> {
    /// The samples scored, in sample order: those with an outer test fold.
    pub fn tested(&self) -> Vec<usize> {
        (0..self.outer.len())
            .filter(|&i| self.outer[i].is_some())
            .collect()
    }

    /// Outer fold `f`'s training samples: every sample outside its test fold
    /// and not held out of it.
    pub fn train(&self, f: usize) -> Bits {
        let n = self.outer.len();
        let mut train = Bits::from_indices(n, (0..n).filter(|&i| self.outer[i] != Some(f)));
        for &i in self.held_out.get(f).map(Vec::as_slice).unwrap_or(&[]) {
            train.remove(i);
        }
        train
    }

    fn check(&self, n: usize) -> Result<()> {
        ensure!(
            self.outer.len() == n,
            "the plan covers {} samples, not {n}",
            self.outer.len()
        );
        ensure!(
            self.inner.len() == self.k && self.holdouts.len() == self.k,
            "the plan has {} inner plans and {} holdout pairs for {} outer folds",
            self.inner.len(),
            self.holdouts.len(),
            self.k
        );
        for (f, inner) in self.inner.iter().enumerate() {
            ensure!(
                inner.len() == n,
                "outer fold {f}'s inner plan covers {} of {n} samples",
                inner.len()
            );
        }
        Ok(())
    }
}

/// One split scored on a plan: the result over the tested samples, and
/// which samples those are.
pub struct PlannedSplit {
    /// `test_fold` and each learner's `scores` are over `tested`, in order.
    pub result: SplitResult,
    /// The scored samples, in sample order (`PlanView::tested`).
    pub tested: Vec<usize>,
    /// Per outer fold, the resistant share of its training samples: the
    /// score a no-call sample of that fold carries in a prediction row.
    pub train_prevalence: Vec<f64>,
}

/// Scores every learner once per outer fold of `plan`. `fold_x(f, train)`
/// returns the feature matrix outer fold `f` fits on, over every sample,
/// built from `train` alone; it must already be the fold's candidate set,
/// since no support filter is applied here. An unavailable plan returns a
/// split with its reason and no learners, as `evaluate` does for a lineage
/// split it cannot score.
pub fn evaluate_plan(
    dataset: &str,
    y: &Bits,
    view: &PlanView,
    cfg: &EvalConfig,
    fold_x: &mut dyn FnMut(usize, &Bits) -> Result<BitMatrix>,
) -> Result<PlannedSplit> {
    let scheme = view.scheme.as_str();
    let tested = view.tested();
    if let Some(reason) = &view.unavailable_reason {
        return Ok(PlannedSplit {
            result: SplitResult {
                scheme: scheme.to_string(),
                folds: view.k,
                n_groups: 0,
                test_fold: vec![],
                learners: vec![],
                unavailable_reason: Some(reason.clone()),
            },
            tested,
            train_prevalence: vec![],
        });
    }
    let n = y.len();
    view.check(n)?;
    let mut fits: Vec<Vec<Scored>> = LEARNERS.iter().map(|_| Vec::new()).collect();
    let mut seconds = vec![0.0; LEARNERS.len()];
    let mut train_prevalence = Vec::with_capacity(view.k);
    for f in 0..view.k {
        let train = view.train(f);
        ensure!(train.count() > 0, "outer fold {f} has no training samples");
        train_prevalence.push(train.and(y).count() as f64 / train.count() as f64);
        let x = fold_x(f, &train).with_context(|| format!("{scheme} outer fold {f}"))?;
        ensure!(
            x.cols == n,
            "outer fold {f}'s matrix covers {} samples, not {n}",
            x.cols
        );
        let inner: Vec<usize> = view.inner[f]
            .iter()
            .map(|g| g.unwrap_or(usize::MAX))
            .collect();
        let held = |val: &[usize]| Bits::from_indices(n, val.iter().copied());
        let (gbt_val, lasso_val) = (held(&view.holdouts[f].gbt), held(&view.holdouts[f].lasso));
        let seed = plan::fold_seed(view.learner_seed.unwrap_or(view.repeat_seed), f);
        let context = LassoContext {
            dataset,
            scheme,
            outer_fold: Some(f),
        };
        for (li, learner) in LEARNERS.iter().enumerate() {
            let t0 = Instant::now();
            let fit = match *learner {
                "scm" => eval::fit_scm(&x, y, &train, &inner, cfg),
                "cart" => eval::fit_cart(&x, y, &train, &inner, cfg),
                "rf" => eval::fit_rf(&x, y, &train, cfg, seed).0,
                "gbt" => fit_gbt(&x, y, &train, &gbt_val, cfg, seed),
                "lasso" => fit_lasso(&x, y, &train, &lasso_val, cfg, context)?,
                _ => unreachable!(),
            };
            fits[li].push(fit);
            seconds[li] += t0.elapsed().as_secs_f64();
        }
    }
    let truth: Vec<bool> = tested.iter().map(|&i| y.contains(i)).collect();
    let test_fold: Vec<usize> = tested
        .iter()
        .map(|&i| view.outer[i].expect("a tested sample has an outer fold"))
        .collect();
    let learners = LEARNERS
        .iter()
        .zip(fits)
        .zip(seconds)
        .map(|((name, fits), seconds)| {
            let fits: Vec<Scored> = fits
                .into_iter()
                .map(|s| Scored {
                    scores: s.scores.map(|all| tested.iter().map(|&i| all[i]).collect()),
                    ..s
                })
                .collect();
            eval::pool_folds(name, &truth, &test_fold, &fits, seconds)
                .with_context(|| format!("{scheme} split"))
        })
        .collect::<Result<Vec<LearnerResult>>>()?;
    Ok(PlannedSplit {
        result: SplitResult {
            scheme: scheme.to_string(),
            folds: view.k,
            // A plan carries no group ids; the plan file and the arm's
            // inputs record them.
            n_groups: 0,
            test_fold,
            learners,
            unavailable_reason: None,
        },
        tested,
        train_prevalence,
    })
}

/// `eval::fit_gbt` with its early-stopping holdout taken from the plan.
fn fit_gbt(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    val: &Bits,
    cfg: &EvalConfig,
    seed: u64,
) -> Scored {
    let fit_m = train.and_not(val);
    let feats = screen(x, y, &fit_m, cfg.screened_features);
    let m = gbt::fit(
        x,
        &Objective::Logistic { y },
        &fit_m,
        val,
        &feats,
        &GbtParams { seed, ..cfg.gbt },
    );
    Scored {
        scores: eval::gbt_scores(&m, x),
        complexity: m.trees.len() as f64,
        lasso_fit: None,
    }
}

/// `eval::fit_lasso` with its penalty-selection holdout taken from the plan.
fn fit_lasso(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    val: &Bits,
    cfg: &EvalConfig,
    context: LassoContext<'_>,
) -> Result<Scored> {
    let fit_m = train.and_not(val);
    let feats = screen(x, y, &fit_m, cfg.screened_features);
    let path = lasso::fit_path_with_status(
        x,
        y,
        &fit_m,
        &feats,
        cfg.lasso_lambdas,
        LASSO_MAX_NONZERO,
        &mut |_| {},
    );
    let (scored, _) = eval::fit_lasso_from_path(x, y, val, path, context, |lambdas| {
        lasso::refit_path_with_status(x, y, train, &feats, lambdas, &mut |_| {})
    })?;
    Ok(scored)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ablate::plan::{GroupIds, fold_plan};
    use crate::learn::training_candidates;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    /// 90 samples in 15 lineages; feature 0 tracks y, the rest are noise.
    fn data() -> (BitMatrix, Bits, Vec<u32>, Vec<Vec<usize>>) {
        let n = 90;
        let mut rng = ChaCha8Rng::seed_from_u64(5);
        let lineage: Vec<u32> = (0..n as u32).map(|i| i / 6).collect();
        let y = Bits::from_indices(
            n,
            (0..n).filter(|&i| lineage[i].is_multiple_of(2) || i.is_multiple_of(7)),
        );
        let mut x = BitMatrix::with_cols(n);
        for f in 0..60 {
            let row = Bits::from_indices(
                n,
                (0..n).filter(|&i| {
                    if f == 0 {
                        y.contains(i) != rng.random_bool(0.1)
                    } else {
                        rng.random_bool(0.3)
                    }
                }),
            );
            x.push_row(row.words());
        }
        // Samples 0 and 1 are k-mer-identical, in one lineage.
        let duplicates = vec![vec![0, 1]];
        let first = x.clone();
        let mut x2 = BitMatrix::with_cols(n);
        for r in 0..first.rows {
            let mut row = Bits::from_words(n, first.row(r).to_vec());
            if row.contains(0) {
                row.insert(1);
            } else {
                row.remove(1);
            }
            x2.push_row(row.words());
        }
        let y = {
            let mut y = y;
            if y.contains(0) {
                y.insert(1);
            } else {
                y.remove(1);
            }
            y
        };
        (x2, y, lineage, duplicates)
    }

    #[test]
    fn the_default_plan_reproduces_evaluate() {
        let (x, y, lineage, duplicates) = data();
        let cfg = EvalConfig {
            forest: tree_params(),
            ..EvalConfig::default()
        };
        let min_support = 3;
        let out =
            eval::evaluate("t", &x, &y, &lineage, &duplicates, min_support, None, &cfg).unwrap();
        let groups = GroupIds::new(lineage.clone(), duplicates.clone(), None).unwrap();
        for split in &out.splits {
            let scheme = if split.scheme == "random" {
                Scheme::Random
            } else {
                Scheme::Lineage
            };
            let p =
                fold_plan(&y, &groups, scheme, cfg.folds, cfg.inner_folds, cfg.seed, 0).unwrap();
            let view = PlanView::from(&p);
            let planned = evaluate_plan("t", &y, &view, &cfg, &mut |_, train| {
                Ok(training_candidates(&x, train, min_support).0)
            })
            .unwrap();
            assert_eq!(planned.tested, (0..y.len()).collect::<Vec<_>>());
            assert_eq!(planned.result.unavailable_reason, split.unavailable_reason);
            assert_eq!(
                planned.result.test_fold, split.test_fold,
                "{}",
                split.scheme
            );
            assert_eq!(planned.result.learners.len(), split.learners.len());
            for (a, b) in planned.result.learners.iter().zip(&split.learners) {
                assert_eq!(a.learner, b.learner);
                let bits = |s: &[f64]| s.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
                assert_eq!(
                    bits(&a.scores),
                    bits(&b.scores),
                    "{} {}",
                    split.scheme,
                    a.learner
                );
                assert_eq!(a.fold_complexity, b.fold_complexity);
                assert_eq!(
                    serde_json::to_value(&a.metrics).unwrap(),
                    serde_json::to_value(&b.metrics).unwrap()
                );
            }
        }
    }

    fn tree_params() -> crate::learn::tree::ForestParams {
        crate::learn::tree::ForestParams {
            n_trees: 20,
            ..Default::default()
        }
    }

    /// An arm over the reference's samples `keep` (reference indices) plus
    /// `extra` genomes of its own: extra j has lineage `extra_lineage(j)`,
    /// resistant when j is even, and reference sample j's features.
    struct Arm {
        ids: Vec<String>,
        y: Bits,
        groups: GroupIds,
        x: BitMatrix,
    }

    fn arm(
        x: &BitMatrix,
        y: &Bits,
        lineage: &[u32],
        keep: &[usize],
        extra: usize,
        extra_lineage: impl Fn(usize) -> u32,
    ) -> Arm {
        let n = keep.len() + extra;
        let source: Vec<usize> = keep.iter().copied().chain(0..extra).collect();
        let resistant = |a: usize| {
            if a < keep.len() {
                y.contains(keep[a])
            } else {
                (a - keep.len()).is_multiple_of(2)
            }
        };
        let lineages = (0..n)
            .map(|a| {
                if a < keep.len() {
                    lineage[keep[a]]
                } else {
                    extra_lineage(a - keep.len())
                }
            })
            .collect();
        let duplicates = if keep.starts_with(&[0, 1]) {
            vec![vec![0, 1]]
        } else {
            vec![]
        };
        let mut arm_x = BitMatrix::with_cols(n);
        for r in 0..x.rows {
            let row = x.row(r);
            arm_x.push_row(
                Bits::from_indices(n, (0..n).filter(|&a| crate::bits::get(row, source[a]))).words(),
            );
        }
        Arm {
            ids: (0..n)
                .map(|a| {
                    if a < keep.len() {
                        format!("g{}", keep[a])
                    } else {
                        format!("extra{}", a - keep.len())
                    }
                })
                .collect(),
            y: Bits::from_indices(n, (0..n).filter(|&a| resistant(a))),
            groups: GroupIds::new(lineages, duplicates, None).unwrap(),
            x: arm_x,
        }
    }

    impl Arm {
        /// The species grouping of this arm and the reference cohort `ids`.
        fn species(
            &self,
            ids: &[String],
            lineage: &[u32],
            duplicates: &[Vec<usize>],
        ) -> plan::SpeciesGroups {
            let r = GroupIds::new(lineage.to_vec(), duplicates.to_vec(), None).unwrap();
            plan::SpeciesGroups::of_cohorts(&[(ids, &r), (&self.ids, &self.groups)])
        }
    }

    fn reference(y: &Bits, lineage: &[u32], duplicates: &[Vec<usize>], scheme: Scheme) -> FoldPlan {
        let groups = GroupIds::new(lineage.to_vec(), duplicates.to_vec(), None).unwrap();
        let cfg = EvalConfig::default();
        fold_plan(y, &groups, scheme, cfg.folds, cfg.inner_folds, cfg.seed, 0).unwrap()
    }

    #[test]
    fn a_cohort_plan_scores_a_class_balanced_common_set_and_trains_on_the_rest() {
        let (x, y, lineage, duplicates) = data();
        let cfg = EvalConfig {
            forest: tree_params(),
            ..EvalConfig::default()
        };
        let ids: Vec<String> = (0..y.len()).map(|i| format!("g{i}")).collect();
        let r = reference(&y, &lineage, &duplicates, Scheme::Random);
        // The arm keeps every resistant reference sample but only every third
        // susceptible one, so its common genomes are mostly resistant, and
        // adds 30 genomes of its own.
        let keep: Vec<usize> = (0..y.len())
            .filter(|&i| y.contains(i) || i.is_multiple_of(3))
            .collect();
        let a = arm(&x, &y, &lineage, &keep, 30, |j| 100 + j as u32);
        let cohort = plan::cohort_plan(
            &r,
            &ids,
            &y,
            &a.ids,
            &a.y,
            &a.species(&ids, &lineage, &duplicates),
            None,
        )
        .unwrap();
        cohort.verify().unwrap();
        assert_eq!(
            (cohort.reference_sha256.as_str(), cohort.balance_seed),
            (r.sha256.as_str(), r.repeat_seed)
        );
        assert_eq!(
            plan::cohort_plan(
                &r,
                &ids,
                &y,
                &a.ids,
                &a.y,
                &a.species(&ids, &lineage, &duplicates),
                None
            )
            .unwrap(),
            cohort
        );
        let common = keep.len();
        for f in 0..cfg.folds {
            let [wr, ws] = cohort.reference_test_classes[f];
            let scored: Vec<usize> = (0..common)
                .filter(|&i| cohort.outer[i] == Some(f))
                .collect();
            let sr = scored.iter().filter(|&&i| a.y.contains(i)).count();
            let ss = scored.len() - sr;
            // Susceptible genomes are the scarce class here: all are kept, and
            // the resistant ones are cut to the reference fold's ratio.
            let fold_s = (0..common)
                .filter(|&i| r.outer[keep[i]] == f && !a.y.contains(i))
                .count();
            assert_eq!(ss, fold_s, "outer fold {f}");
            assert_eq!(
                sr,
                ((ss * wr) as f64 / ws as f64).round() as usize,
                "outer fold {f}"
            );
            assert!(scored.iter().all(|&i| r.outer[keep[i]] == f));
            let dropped: Vec<usize> = (0..common)
                .filter(|&i| r.outer[keep[i]] == f && cohort.outer[i].is_none())
                .collect();
            assert_eq!(cohort.held_out[f], dropped);
            assert_eq!(
                (cohort.balance_dropped[f], cohort.lineage_dropped[f]),
                (dropped.len(), 0)
            );
            assert!(cohort.balance_dropped[f] > 0, "outer fold {f}");
        }
        assert!(cohort.outer[common..].iter().all(Option::is_none));
        let view = PlanView::from(&cohort);
        let mut trained = vec![];
        let planned = evaluate_plan("t", &a.y, &view, &cfg, &mut |f, train| {
            trained.push((f, train.clone()));
            Ok(training_candidates(&a.x, train, 3).0)
        })
        .unwrap();
        assert_eq!(planned.tested, view.tested());
        for (f, train) in trained {
            assert!(
                (common..common + 30).all(|i| train.contains(i)),
                "outer fold {f}"
            );
            // Every common genome the reference tests in f is out of training.
            assert!((0..common).all(|i| train.contains(i) == (r.outer[keep[i]] != f)));
        }
    }

    /// On the lineage split, an arm-only genome of a scored genome's lineage
    /// never trains in that fold; this test fails without that drop.
    #[test]
    fn a_cohort_arm_never_trains_on_a_test_lineage_on_the_lineage_split() {
        let (x, y, lineage, duplicates) = data();
        let ids: Vec<String> = (0..y.len()).map(|i| format!("g{i}")).collect();
        let keep: Vec<usize> = (0..60).collect();
        // Extra genomes in the lineages of the kept ones (lineages 0 to 9).
        let a = arm(&x, &y, &lineage, &keep, 30, |j| (j % 10) as u32);
        let r = reference(&y, &lineage, &duplicates, Scheme::Lineage);
        let cohort = plan::cohort_plan(
            &r,
            &ids,
            &y,
            &a.ids,
            &a.y,
            &a.species(&ids, &lineage, &duplicates),
            None,
        )
        .unwrap();
        cohort.verify().unwrap();
        let lineage_of = |i: usize| a.groups.lineage()[i];
        let mut dropped = 0;
        for f in 0..r.k {
            let scored: std::collections::BTreeSet<u32> = (0..a.y.len())
                .filter(|&i| cohort.outer[i] == Some(f))
                .map(lineage_of)
                .collect();
            let train = cohort.train(f);
            let leaking: Vec<usize> = (60..90)
                .filter(|&i| train.contains(i) && scored.contains(&lineage_of(i)))
                .collect();
            assert!(leaking.is_empty(), "outer fold {f} trains on {leaking:?}");
            let expect = (60..90)
                .filter(|&i| scored.contains(&lineage_of(i)))
                .count();
            assert_eq!(cohort.lineage_dropped[f], expect, "outer fold {f}");
            dropped += expect;
        }
        assert!(dropped > 0);
        // The random split keeps every arm-only genome in training.
        let random = reference(&y, &lineage, &duplicates, Scheme::Random);
        let cohort = plan::cohort_plan(
            &random,
            &ids,
            &y,
            &a.ids,
            &a.y,
            &a.species(&ids, &lineage, &duplicates),
            None,
        )
        .unwrap();
        assert!(cohort.lineage_dropped.iter().all(|&d| d == 0));
    }

    /// The calibration arms of the input comparison: the default's cohort sent
    /// through the cohort path (with its own subsample seed) is the default,
    /// bit for bit; the learner-seed twin keeps every partition and changes
    /// only the forest's and boosting's draws.
    #[test]
    fn the_calibration_arms_are_exchangeable_with_the_default() {
        let (x, y, lineage, duplicates) = data();
        let cfg = EvalConfig {
            forest: tree_params(),
            ..EvalConfig::default()
        };
        let ids: Vec<String> = (0..y.len()).map(|i| format!("g{i}")).collect();
        let groups = GroupIds::new(lineage.clone(), duplicates.clone(), None).unwrap();
        let fold_x = |_: usize, train: &Bits| Ok(training_candidates(&x, train, 3).0);
        for scheme in [Scheme::Random, Scheme::Lineage] {
            // The comparison's reference plan is keyed on genome id.
            let r = plan::keyed_fold_plan(&ids, &y, &groups, scheme, 5, 3, 42, 0).unwrap();
            if !r.is_available() {
                assert_eq!(scheme, Scheme::Lineage, "{:?}", r.unavailable_reason);
                continue;
            }
            let base =
                evaluate_plan("t", &y, &PlanView::from(&r), &cfg, &mut fold_x.clone()).unwrap();
            let species = plan::SpeciesGroups::of_cohorts(&[(&ids, &groups)]);
            let cohort = plan::cohort_plan(&r, &ids, &y, &ids, &y, &species, Some(7)).unwrap();
            assert_eq!(cohort.balance_seed, 7);
            assert!(cohort.held_out.iter().all(Vec::is_empty), "{scheme}");
            assert_eq!(cohort.inner, r.inner);
            assert_eq!(cohort.holdouts, r.holdouts);
            let path = evaluate_plan("t", &y, &PlanView::from(&cohort), &cfg, &mut fold_x.clone())
                .unwrap();
            let mut view = PlanView::from(&r);
            view.learner_seed = Some(999);
            let twin = evaluate_plan("t", &y, &view, &cfg, &mut fold_x.clone()).unwrap();
            let bits = |s: &[f64]| s.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
            for ((a, b), c) in base
                .result
                .learners
                .iter()
                .zip(&path.result.learners)
                .zip(&twin.result.learners)
            {
                assert_eq!(bits(&a.scores), bits(&b.scores), "{scheme} {}", a.learner);
                // scm, cart and the lasso draw nothing; the forest always
                // does (boosting may keep no tree either way).
                let same = bits(&a.scores) == bits(&c.scores);
                match a.learner.as_str() {
                    "rf" => assert!(!same, "{scheme} rf"),
                    "gbt" => {}
                    other => assert!(same, "{scheme} {other}"),
                }
            }
            assert_eq!(base.result.test_fold, twin.result.test_fold);
        }
    }

    /// A genome of a scored genome's k-mer duplicate group never trains in
    /// that fold, on either split; this test fails without the duplicate drop.
    #[test]
    fn a_cohort_arm_never_trains_on_a_copy_of_a_scored_genome() {
        let (x, y, lineage, duplicates) = data();
        let ids: Vec<String> = (0..y.len()).map(|i| format!("g{i}")).collect();
        let keep: Vec<usize> = (0..60).collect();
        // Extra genome 0 (arm sample 60) is a copy of reference sample 2:
        // same lineage, same label, one duplicate group with it.
        let mut a = arm(&x, &y, &lineage, &keep, 10, |j| 100 + j as u32);
        let mut lin = a.groups.lineage().to_vec();
        lin[60] = lin[2];
        if a.y.contains(2) {
            a.y.insert(60);
        } else {
            a.y.remove(60);
        }
        a.groups = GroupIds::new(lin, vec![vec![0, 1], vec![2, 60]], None).unwrap();
        for scheme in [Scheme::Random, Scheme::Lineage] {
            let r = reference(&y, &lineage, &duplicates, scheme);
            let cohort = plan::cohort_plan(
                &r,
                &ids,
                &y,
                &a.ids,
                &a.y,
                &a.species(&ids, &lineage, &duplicates),
                None,
            )
            .unwrap();
            cohort.verify().unwrap();
            let Some(f) = cohort.outer[2] else {
                continue;
            };
            assert!(
                !cohort.train(f).contains(60),
                "{scheme}: the copy trains in its original's fold"
            );
            assert!(cohort.held_out[f].contains(&60), "{scheme}");
            if scheme == Scheme::Random {
                assert!(cohort.duplicate_dropped[f] >= 1, "{scheme}");
            }
            for g in (0..r.k).filter(|&g| g != f) {
                assert!(cohort.train(g).contains(60) || cohort.held_out[g].contains(&60));
            }
        }
    }

    /// The comparison's fold plans are keyed on genome identity (main's ruling
    /// on REVIEW-INPUT2 N1): adding genomes moves no shared genome's outer
    /// fold, inner fold or holdout membership.
    #[test]
    fn keyed_folds_do_not_depend_on_the_other_genomes() {
        let (_, y, lineage, duplicates) = data();
        let n = y.len();
        let ids: Vec<String> = (0..n).map(|i| format!("g{i:03}")).collect();
        let full = GroupIds::new(lineage.clone(), duplicates.clone(), None).unwrap();
        let keep: Vec<usize> = (0..n).filter(|i| i % 17 != 5).collect();
        let sub_ids: Vec<String> = keep.iter().map(|&i| ids[i].clone()).collect();
        let sub_y =
            Bits::from_indices(keep.len(), (0..keep.len()).filter(|&j| y.contains(keep[j])));
        let sub = GroupIds::new(
            keep.iter().map(|&i| lineage[i]).collect(),
            vec![vec![0, 1]],
            None,
        )
        .unwrap();
        for scheme in [Scheme::Random, Scheme::Lineage] {
            let a = plan::keyed_fold_plan(&ids, &y, &full, scheme, 5, 3, 42, 0).unwrap();
            let b = plan::keyed_fold_plan(&sub_ids, &sub_y, &sub, scheme, 5, 3, 42, 0).unwrap();
            assert_eq!(a.format, plan::KEYED_PLAN_FORMAT);
            a.verify().unwrap();
            for (j, &i) in keep.iter().enumerate() {
                assert_eq!(a.outer[i], b.outer[j], "{scheme} outer fold of {}", ids[i]);
            }
            if !(a.is_available() && b.is_available()) {
                assert_eq!(scheme, Scheme::Lineage);
                continue;
            }
            for f in 0..5 {
                for (j, &i) in keep.iter().enumerate() {
                    assert_eq!(
                        a.inner[f][i], b.inner[f][j],
                        "{scheme} inner fold of {} in {f}",
                        ids[i]
                    );
                    for (ha, hb) in [
                        (&a.holdouts[f].gbt, &b.holdouts[f].gbt),
                        (&a.holdouts[f].lasso, &b.holdouts[f].lasso),
                    ] {
                        assert_eq!(
                            ha.contains(&i),
                            hb.contains(&j),
                            "{scheme} holdout of {} in {f}",
                            ids[i]
                        );
                    }
                }
            }
        }
        // The draw is keyed on ids, not on positions: reversing the sample order
        // moves every sample with its own fold.
        let rev: Vec<usize> = (0..n).rev().collect();
        let rev_ids: Vec<String> = rev.iter().map(|&i| ids[i].clone()).collect();
        let rev_y = Bits::from_indices(n, (0..n).filter(|&j| y.contains(rev[j])));
        let rev_dup = vec![vec![n - 2, n - 1]];
        let rev_groups =
            GroupIds::new(rev.iter().map(|&i| lineage[i]).collect(), rev_dup, None).unwrap();
        let a = plan::keyed_fold_plan(&ids, &y, &full, Scheme::Random, 5, 3, 42, 0).unwrap();
        let r = plan::keyed_fold_plan(&rev_ids, &rev_y, &rev_groups, Scheme::Random, 5, 3, 42, 0)
            .unwrap();
        for (j, &i) in rev.iter().enumerate() {
            assert_eq!(a.outer[i], r.outer[j]);
        }
    }

    /// A keyed plan whose folds can't all hold both classes is refused, with the
    /// reason, rather than scored.
    #[test]
    fn a_keyed_plan_is_refused_when_a_fold_lacks_a_class() {
        let n = 60;
        let ids: Vec<String> = (0..n).map(|i| format!("g{i}")).collect();
        let y = Bits::from_indices(n, [0, 1]);
        let groups = GroupIds::new((0..n as u32).collect(), vec![], None).unwrap();
        let p = plan::keyed_fold_plan(&ids, &y, &groups, Scheme::Random, 5, 3, 42, 0).unwrap();
        let why = p
            .unavailable_reason
            .clone()
            .expect("a fold without a resistant sample");
        assert!(why.contains("invalid random outer evaluation"), "{why}");
        assert!(!p.is_available());
        p.verify().unwrap();
    }
}
