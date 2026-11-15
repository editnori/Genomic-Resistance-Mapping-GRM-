//! Experiment, off by default: refit the outer folds whose boosted model
//! kept no tree or whose lasso kept no weight (null fits, see
//! docs/equal-scoring.md), under fallback settings, and write what each
//! fallback scores on the fold's test samples to a side file.
//!
//! Nothing here changes a report. `grm learn --exp-null-refit FILE` turns it
//! on for that process; without the flag `enabled` is false and evaluation
//! runs exactly as before. With it, evaluation still scores every null fold
//! at its training prevalence, as revision 14 does, and only appends one
//! JSON line per null fold to FILE. docs/null-fold-refit.md describes the
//! fallbacks and the plan they serve.
//!
//! The fallbacks redo the production fit's own steps on the fold's training
//! samples, with the same seed, holdout and screened features:
//! - lasso: the tuning path, then the refit on the whole training fold down
//!   through every path point, so each penalty's refit is recorded. The
//!   production rule selected the point of lowest holdout log loss; on a
//!   null fold that is point 0, the largest penalty, where every weight is
//!   zero. The record keeps every point, so any fixed rule can be read off.
//! - boosting: the same fit with early stopping off, all `rounds` trees on
//!   the fitting samples, with the validation loss after each round.

use super::eval::{EvalConfig, GBT_HOLDOUT, LASSO_HOLDOUT, LassoContext, tuning_unit};
use super::gbt::{self, GbtParams, Objective};
use super::{lasso, screen, tuning_holdout};
use crate::bits::{BitMatrix, Bits};
use anyhow::{Context, Result};
use serde::Serialize;
use std::io::Write;
use std::path::Path;
use std::sync::{Mutex, OnceLock};

static SINK: OnceLock<Mutex<std::fs::File>> = OnceLock::new();

#[cfg(test)]
thread_local! {
    /// Tests only: while `capture` runs on a thread, that thread's records
    /// go here and not to FILE, so a unit test reads its own records
    /// without setting FILE for every other test of the process.
    static CAPTURED: std::cell::RefCell<Option<Vec<String>>> =
        const { std::cell::RefCell::new(None) };
}

/// Turns the experiment on for this process: every later null fold appends
/// one line to `path`, which is created if missing. A second call fails.
pub fn enable(path: &Path) -> Result<()> {
    let file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)
        .with_context(|| format!("opening {}", path.display()))?;
    SINK.set(Mutex::new(file))
        .map_err(|_| anyhow::anyhow!("the null-fold refit file is already set"))
}

/// Whether `enable` has run in this process (or, in a unit test, whether
/// `capture` runs on this thread).
pub fn enabled() -> bool {
    SINK.get().is_some() || capturing()
}

#[cfg(not(test))]
fn capturing() -> bool {
    false
}

#[cfg(test)]
fn capturing() -> bool {
    CAPTURED.with(|c| c.borrow().is_some())
}

/// Tests only: runs `f` with the experiment on for this thread alone, and
/// returns `f`'s result with the lines it would have appended to FILE.
#[cfg(test)]
pub(crate) fn capture<R>(f: impl FnOnce() -> R) -> (R, Vec<String>) {
    CAPTURED.with(|c| *c.borrow_mut() = Some(Vec::new()));
    let result = f();
    let lines = CAPTURED.with(|c| c.borrow_mut().take()).unwrap_or_default();
    (result, lines)
}

/// One penalty point of a null fold's lasso path and its refit.
#[derive(Clone, Debug, Serialize)]
pub struct LassoPoint {
    pub lambda: f64,
    /// Distinct non-zero weights of the path model, fitted without the holdout.
    pub path_weights: usize,
    pub path_converged: bool,
    /// Class-balanced log loss of the path model on the holdout, the
    /// production selection criterion (`lasso::log_loss`).
    pub holdout_log_loss: f64,
    /// Balanced accuracy of the path model's calls on the holdout at 0.5.
    pub holdout_balanced_accuracy: f64,
    /// Non-zero weights of the refit on the whole training fold.
    pub refit_weights: usize,
    pub refit_converged: bool,
    /// The refit's score for each sample of `NullRefit::test`, in order.
    pub test_scores: Vec<f64>,
}

#[derive(Clone, Debug, Serialize)]
pub struct LassoArms {
    pub holdout_n: usize,
    pub holdout_resistant: usize,
    /// The point the production rule selects: lowest holdout log loss among
    /// converged points, ties to the first. On a null fold it is 0.
    pub production_point: Option<usize>,
    /// The last point converged on both the path and the refit: the
    /// smallest penalty reached.
    pub smallest_point: Option<usize>,
    /// The production rule over the converged points other than 0 that keep
    /// a weight on the path.
    pub next_best_point: Option<usize>,
    pub points: Vec<LassoPoint>,
}

#[derive(Clone, Debug, Serialize)]
pub struct GbtArms {
    pub holdout_n: usize,
    pub holdout_resistant: usize,
    /// Class-weighted validation loss of the 0.5 prior, then after each of
    /// the `rounds` trees.
    pub holdout_loss: Vec<f64>,
    /// Trees early stopping keeps with no patience limit (`gbt::best_round`).
    pub best_round: usize,
    /// Scores on `NullRefit::test` of the model with every round's tree.
    pub full_test_scores: Vec<f64>,
    /// Scores on `NullRefit::test` with the first `best_round` trees, or
    /// None when `best_round` is 0.
    pub best_round_test_scores: Option<Vec<f64>>,
}

/// One null fold, as written to the side file.
#[derive(Clone, Debug, Serialize)]
pub struct NullRefit {
    pub dataset: String,
    pub scheme: String,
    pub fold: usize,
    pub learner: String,
    pub train_n: usize,
    pub train_resistant: usize,
    /// The outer fold's test samples: every sample outside `train`.
    pub test: Vec<usize>,
    pub lasso: Option<LassoArms>,
    pub gbt: Option<GbtArms>,
}

fn balanced_accuracy(scores: &[f64], y: &Bits, mask: &Bits) -> f64 {
    1.0 - super::tree::balanced_error(scores, y, mask)
}

fn pick(scores: &[f64], idx: &[usize]) -> Vec<f64> {
    idx.iter().map(|&i| scores[i]).collect()
}

/// The lasso fallbacks for one outer fold, from `fit_lasso`'s own steps.
#[allow(clippy::too_many_arguments)]
pub(crate) fn lasso_arms(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    groups: Option<&[u32]>,
    cfg: &EvalConfig,
    seed: u64,
    context: LassoContext<'_>,
    test: &[usize],
) -> Result<LassoArms> {
    let unit = tuning_unit(context.scheme);
    let (fit_m, val) = tuning_holdout(train, y, groups, unit, LASSO_HOLDOUT, seed)?;
    let feats = screen(x, y, &fit_m, cfg.screened_features);
    let path =
        lasso::fit_path_with_status(x, y, &fit_m, &feats, cfg.lasso_lambdas, 400, &mut |_| {});
    let lambdas: Vec<f64> = path.iter().map(|f| f.model.lambda).collect();
    let refit = if lambdas.is_empty() {
        vec![]
    } else {
        lasso::refit_path_with_status(x, y, train, &feats, &lambdas, &mut |_| {})
    };
    let points: Vec<LassoPoint> = path
        .iter()
        .enumerate()
        .map(|(k, p)| {
            let on_x = p.model.predict_proba(x);
            let r = refit.get(k);
            LassoPoint {
                lambda: p.model.lambda,
                path_weights: p.model.coef.len(),
                path_converged: p.convergence.converged,
                holdout_log_loss: lasso::log_loss(&on_x, y, &val),
                holdout_balanced_accuracy: balanced_accuracy(&on_x, y, &val),
                refit_weights: r.map_or(0, |r| r.model.coef.len()),
                refit_converged: r.is_some_and(|r| r.convergence.converged),
                test_scores: r.map_or_else(Vec::new, |r| pick(&r.model.predict_proba(x), test)),
            }
        })
        .collect();
    let lowest = |keep: &dyn Fn(usize, &LassoPoint) -> bool| {
        points
            .iter()
            .enumerate()
            .filter(|(k, p)| keep(*k, p))
            .min_by(|a, b| a.1.holdout_log_loss.total_cmp(&b.1.holdout_log_loss))
            .map(|(k, _)| k)
    };
    Ok(LassoArms {
        holdout_n: val.count() as usize,
        holdout_resistant: val.and(y).count() as usize,
        production_point: lowest(&|_, p| p.path_converged),
        smallest_point: points
            .iter()
            .rposition(|p| p.path_converged && p.refit_converged),
        next_best_point: lowest(&|k, p| k > 0 && p.path_converged && p.path_weights > 0),
        points,
    })
}

/// The boosting fallbacks for one outer fold, from `fit_gbt`'s own steps:
/// `seed` draws the early-stopping holdout and `draws` the column samples.
#[allow(clippy::too_many_arguments)]
pub(crate) fn gbt_arms(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    groups: Option<&[u32]>,
    unit: &str,
    cfg: &EvalConfig,
    seed: u64,
    draws: u64,
    test: &[usize],
) -> Result<GbtArms> {
    let (fit_m, val) = tuning_holdout(train, y, groups, unit, GBT_HOLDOUT, seed)?;
    let feats = screen(x, y, &fit_m, cfg.screened_features);
    let (full, trace) = gbt::fit_traced(
        x,
        &Objective::Logistic { y },
        &fit_m,
        &val,
        &feats,
        &GbtParams {
            seed: draws,
            ..cfg.gbt
        },
        false,
    );
    let best = gbt::best_round(&trace);
    let best_scores = (best > 0).then(|| {
        let kept = gbt::Gbt {
            trees: full.trees[..best].to_vec(),
            ..full.clone()
        };
        pick(&kept.predict_proba(x), test)
    });
    Ok(GbtArms {
        holdout_n: val.count() as usize,
        holdout_resistant: val.and(y).count() as usize,
        holdout_loss: trace,
        best_round: best,
        full_test_scores: pick(&full.predict_proba(x), test),
        best_round_test_scores: best_scores,
    })
}

/// Fits `learner`'s fallbacks on one null outer fold and appends its line.
/// Called by evaluation only when `enabled`, with the fold's seeds as
/// `fit_fold` derives them: `seed` for every partition, `draws` for
/// boosting's column samples (equal unless `EvalConfig::learner_seed` is set).
#[allow(clippy::too_many_arguments)]
pub(crate) fn record_fold(
    learner: &str,
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    groups: Option<&[u32]>,
    cfg: &EvalConfig,
    seed: u64,
    draws: u64,
    context: LassoContext<'_>,
) -> Result<()> {
    let test: Vec<usize> = (0..x.cols).filter(|&i| !train.contains(i)).collect();
    let fold = context.outer_fold.expect("a null fold is an outer fold");
    let rec = NullRefit {
        dataset: context.dataset.to_string(),
        scheme: context.scheme.to_string(),
        fold,
        learner: learner.to_string(),
        train_n: train.count() as usize,
        train_resistant: train.and(y).count() as usize,
        lasso: match learner {
            "lasso" => Some(lasso_arms(x, y, train, groups, cfg, seed, context, &test)?),
            _ => None,
        },
        gbt: match learner {
            "gbt" => Some(gbt_arms(
                x,
                y,
                train,
                groups,
                tuning_unit(context.scheme),
                cfg,
                seed,
                draws,
                &test,
            )?),
            _ => None,
        },
        test,
    };
    write(&rec)
}

fn write(rec: &NullRefit) -> Result<()> {
    let mut line = serde_json::to_string(rec)?;
    line.push('\n');
    #[cfg(test)]
    if keep_captured(&line) {
        return Ok(());
    }
    let Some(sink) = SINK.get() else {
        return Ok(());
    };
    let mut file = sink.lock().expect("the null-fold refit file lock");
    file.write_all(line.as_bytes())?;
    file.flush()?;
    Ok(())
}

/// Tests only: keeps `line` when `capture` runs on this thread, and says
/// whether it did.
#[cfg(test)]
fn keep_captured(line: &str) -> bool {
    CAPTURED.with(|c| match c.borrow_mut().as_mut() {
        Some(lines) => {
            lines.push(line.to_string());
            true
        }
        None => false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::learn::eval::{fit_gbt_for_tests, fit_lasso_for_tests};

    fn config() -> EvalConfig {
        EvalConfig {
            folds: 2,
            gbt: GbtParams {
                rounds: 40,
                patience: 5,
                max_depth: 2,
                ..Default::default()
            },
            lasso_lambdas: 12,
            ..Default::default()
        }
    }

    /// 230 training samples and 30 test samples, one in nine resistant. The
    /// one feature marks the resistant fitting samples and the susceptible
    /// holdout samples, so every weight raises the holdout's log loss and
    /// the production rule keeps none: the released lasso cause.
    fn null_lasso_fixture() -> (BitMatrix, Bits, Bits, Vec<usize>) {
        let n = 260;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 9 == 0));
        let train = Bits::from_indices(n, 0..230);
        let (fit, val) = tuning_holdout(&train, &y, None, "group", LASSO_HOLDOUT, 42).unwrap();
        let test_resistant = Bits::from_indices(n, (230..n).filter(|i| i % 9 == 0));
        let mut x = BitMatrix::with_cols(n);
        x.push_row(fit.and(&y).or(&val.and_not(&y)).or(&test_resistant).words());
        (x, y, train, (230..n).collect())
    }

    #[test]
    fn a_null_lasso_fold_records_every_points_refit_and_both_arms() {
        let (x, y, train, test) = null_lasso_fixture();
        let cfg = config();
        let context = LassoContext {
            dataset: "fixture__fixturemycin",
            scheme: "random",
            outer_fold: Some(0),
        };
        // Production: a null fit, the holdout chose point 0.
        let (scored, audit) = fit_lasso_for_tests(&x, &y, &train, &cfg, 42, context);
        assert!(scored.is_none());
        assert_eq!(audit, Some(0));
        let arms = lasso_arms(&x, &y, &train, None, &cfg, 42, context, &test).unwrap();
        assert_eq!(arms.production_point, Some(0));
        assert!(arms.points.len() > 2);
        assert_eq!(arms.points[0].path_weights, 0);
        assert_eq!(arms.points[0].refit_weights, 0);
        // Every point's record is the production refit down to that point,
        // on the whole training fold.
        let (fit, val) = tuning_holdout(&train, &y, None, "group", LASSO_HOLDOUT, 42).unwrap();
        let feats = screen(&x, &y, &fit, cfg.screened_features);
        let lambdas: Vec<f64> = arms.points.iter().map(|p| p.lambda).collect();
        for k in 0..lambdas.len() {
            let refit = lasso::refit_path(&x, &y, &train, &feats, &lambdas[..=k]);
            let m = refit.last().unwrap();
            let p = &arms.points[k];
            assert_eq!(p.refit_weights, m.coef.len(), "point {k}");
            assert_eq!(
                p.test_scores,
                pick(&m.predict_proba(&x), &test),
                "point {k}"
            );
            assert_eq!(p.test_scores.len(), test.len());
        }
        let path = lasso::fit_path(&x, &y, &fit, &feats, cfg.lasso_lambdas, 400);
        for (k, m) in path.iter().enumerate() {
            let loss = lasso::log_loss(&m.predict_proba(&x), &y, &val);
            assert_eq!(arms.points[k].holdout_log_loss, loss, "point {k}");
        }
        // The smallest penalty is the last point; the next best is the
        // lowest holdout loss after point 0, and a later point fits.
        assert_eq!(arms.smallest_point, Some(arms.points.len() - 1));
        let next = arms.next_best_point.unwrap();
        assert!(next > 0);
        for (k, p) in arms.points.iter().enumerate().skip(1) {
            assert!(
                p.holdout_log_loss >= arms.points[next].holdout_log_loss,
                "point {k}"
            );
            if k < next {
                assert!(p.holdout_log_loss > arms.points[next].holdout_log_loss);
            }
        }
        assert!(arms.points[arms.smallest_point.unwrap()].refit_weights > 0);
        assert_eq!(arms.holdout_n, val.count() as usize);
        assert_eq!(arms.holdout_resistant, val.and(&y).count() as usize);
    }

    #[test]
    fn the_gbt_fit_with_early_stopping_is_the_production_fit_and_off_keeps_every_round() {
        // Half resistant. The first matrix's feature marks resistant fitting
        // samples and susceptible holdout samples, so early stopping keeps no
        // tree; the second's marks most resistant samples, so it keeps some.
        let n = 80;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let all = Bits::ones(n);
        let (fit, val) = tuning_holdout(&all, &y, None, "group", GBT_HOLDOUT, 42).unwrap();
        let mut null_x = BitMatrix::with_cols(n);
        null_x.push_row(fit.and(&y).or(&val.and_not(&y)).words());
        let mut fitted_x = BitMatrix::with_cols(n);
        fitted_x
            .push_row(Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0 && i % 10 != 0)).words());
        fitted_x.push_row(Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0)).words());
        let cfg = config();
        let params = GbtParams {
            seed: 42,
            ..cfg.gbt
        };
        let obj = Objective::Logistic { y: &y };
        let json = |m: &gbt::Gbt| serde_json::to_string(m).unwrap();
        for (x, null) in [(&null_x, true), (&fitted_x, false)] {
            let feats = screen(x, &y, &fit, cfg.screened_features);
            let production = gbt::fit(x, &obj, &fit, &val, &feats, &params);
            assert_eq!(production.trees.is_empty(), null);
            assert_eq!(fit_gbt_for_tests(x, &y, &all, &cfg, 42).is_none(), null);
            let (stopped, trace) = gbt::fit_traced(x, &obj, &fit, &val, &feats, &params, true);
            assert_eq!(json(&stopped), json(&production));
            assert!(trace.len() < cfg.gbt.rounds + 1, "null {null}");
            assert_eq!(gbt::best_round(&trace), production.trees.len());
            let (full, full_trace) = gbt::fit_traced(x, &obj, &fit, &val, &feats, &params, false);
            assert_eq!(full.trees.len(), cfg.gbt.rounds);
            assert_eq!(full_trace.len(), cfg.gbt.rounds + 1);
            // The stopped run is the start of the full one.
            assert_eq!(full_trace[..trace.len()], trace[..]);
            let kept = gbt::Gbt {
                trees: full.trees[..production.trees.len()].to_vec(),
                ..full.clone()
            };
            assert_eq!(json(&kept), json(&production));
            let arms = gbt_arms(x, &y, &all, None, "group", &cfg, 42, 42, &[0, 1, 2]).unwrap();
            assert_eq!(arms.holdout_loss, full_trace);
            let best = gbt::best_round(&full_trace);
            assert_eq!(arms.best_round, best);
            assert_eq!(
                arms.full_test_scores,
                pick(&full.predict_proba(x), &[0, 1, 2])
            );
            assert_eq!(arms.best_round_test_scores.is_some(), best > 0);
            assert_eq!(arms.holdout_n, val.count() as usize);
        }
    }

    #[test]
    fn best_round_is_the_last_round_that_lowers_the_loss_by_more_than_1e_9() {
        assert_eq!(best_round_of(&[0.7, 0.7, 0.71]), 0);
        assert_eq!(best_round_of(&[0.7, 0.7 - 1e-10, 0.7 - 2e-10]), 0);
        assert_eq!(best_round_of(&[0.7, 0.8, 0.6, 0.65, 0.5, 0.5]), 4);
        assert_eq!(best_round_of(&[0.7, 0.69]), 1);
    }

    fn best_round_of(trace: &[f64]) -> usize {
        gbt::best_round(trace)
    }

    fn fixture_context(scheme: &str) -> LassoContext<'_> {
        LassoContext {
            dataset: "fixture__fixturemycin",
            scheme,
            outer_fold: Some(1),
        }
    }

    /// Balanced accuracy of the calls at 0.5 on `mask`, counted here rather
    /// than through `tree::balanced_error`.
    fn balanced_accuracy_by_hand(scores: &[f64], y: &Bits, mask: &Bits) -> f64 {
        let (mut tp, mut p, mut tn, mut n) = (0.0, 0.0, 0.0, 0.0);
        for i in mask.iter_ones() {
            if y.contains(i) {
                p += 1.0;
                tp += f64::from(u8::from(scores[i] >= 0.5));
            } else {
                n += 1.0;
                tn += f64::from(u8::from(scores[i] < 0.5));
            }
        }
        (tp / p + tn / n) / 2.0
    }

    /// Checks every point of `arms` against the production steps redone
    /// here: the tuning holdout of `seed`, features screened on its fitting
    /// part, the path on that part, the refit on the whole training fold.
    /// Returns the screened features.
    fn assert_points_redo_production(
        x: &BitMatrix,
        y: &Bits,
        train: &Bits,
        cfg: &EvalConfig,
        test: &[usize],
        arms: &LassoArms,
    ) -> Vec<u32> {
        let (fit, val) = tuning_holdout(train, y, None, "group", LASSO_HOLDOUT, 42).unwrap();
        let feats = screen(x, y, &fit, cfg.screened_features);
        let path = lasso::fit_path(x, y, &fit, &feats, cfg.lasso_lambdas, 400);
        assert_eq!(arms.points.len(), path.len());
        let lambdas: Vec<f64> = path.iter().map(|m| m.lambda).collect();
        for (k, (m, p)) in path.iter().zip(&arms.points).enumerate() {
            let on_x = m.predict_proba(x);
            assert_eq!(p.lambda, m.lambda, "point {k}");
            assert_eq!(p.path_weights, m.coef.len(), "point {k}");
            assert_eq!(
                p.holdout_log_loss,
                lasso::log_loss(&on_x, y, &val),
                "point {k}"
            );
            let by_hand = balanced_accuracy_by_hand(&on_x, y, &val);
            assert!(
                (p.holdout_balanced_accuracy - by_hand).abs() < 1e-12,
                "point {k}"
            );
            let refit = lasso::refit_path(x, y, train, &feats, &lambdas[..=k]);
            let r = refit.last().unwrap();
            assert_eq!(p.refit_weights, r.coef.len(), "point {k}");
            assert_eq!(p.test_scores, pick(&r.predict_proba(x), test), "point {k}");
        }
        assert_eq!(arms.holdout_n, val.count() as usize);
        assert_eq!(arms.holdout_resistant, val.and(y).count() as usize);
        feats
    }

    /// A record counts the fold's training genomes and its resistant ones,
    /// names its test genomes, and holds the arms of its learner only. The
    /// counts are worked by hand: genomes 0 to 229 train, every ninth is
    /// resistant (0, 9, ..., 225), and 230 to 259 are the test genomes.
    #[test]
    fn a_record_counts_the_training_genomes_and_names_the_test_genomes() {
        let (x, y, train, test) = null_lasso_fixture();
        let cfg = config();
        let context = fixture_context("random");
        let (done, lines) =
            capture(|| record_fold("lasso", &x, &y, &train, None, &cfg, 42, 42, context));
        done.unwrap();
        assert_eq!(lines.len(), 1);
        assert!(lines[0].ends_with('\n'));
        let r: serde_json::Value = serde_json::from_str(&lines[0]).unwrap();
        assert_eq!(r["dataset"], "fixture__fixturemycin");
        assert_eq!(r["scheme"], "random");
        assert_eq!(r["fold"], 1);
        assert_eq!(r["learner"], "lasso");
        assert_eq!(r["train_n"], 230);
        assert_eq!(r["train_resistant"], 26);
        assert_eq!(r["test"], serde_json::json!((230..260).collect::<Vec<_>>()));
        assert!(r["gbt"].is_null());
        let arms = lasso_arms(&x, &y, &train, None, &cfg, 42, context, &test).unwrap();
        assert_eq!(r["lasso"], serde_json::to_value(&arms).unwrap());
    }

    /// Each point's holdout balanced accuracy is the path model's, read on
    /// the holdout and not on the fitting part. On this fixture the two
    /// differ at some point, so reading the wrong part fails the test.
    #[test]
    fn each_points_holdout_balanced_accuracy_is_read_on_the_holdout() {
        let (x, y, train, test) = null_lasso_fixture();
        let cfg = config();
        let context = fixture_context("random");
        let arms = lasso_arms(&x, &y, &train, None, &cfg, 42, context, &test).unwrap();
        let feats = assert_points_redo_production(&x, &y, &train, &cfg, &test, &arms);
        let (fit, val) = tuning_holdout(&train, &y, None, "group", LASSO_HOLDOUT, 42).unwrap();
        let path = lasso::fit_path(&x, &y, &fit, &feats, cfg.lasso_lambdas, 400);
        let differ = path.iter().any(|m| {
            let on_x = m.predict_proba(&x);
            let gap = balanced_accuracy_by_hand(&on_x, &y, &fit)
                - balanced_accuracy_by_hand(&on_x, &y, &val);
            gap.abs() > 1e-9
        });
        assert!(differ);
    }

    /// The lasso arms screen features on the tuning holdout's fitting part
    /// alone, as production does. Feature 0 marks the resistant fitting
    /// genomes; feature 1 marks every resistant training genome but two
    /// fitting ones. Screening one feature picks 0 on the fitting part and
    /// 1 on the whole training fold, and the two paths differ.
    #[test]
    fn the_lasso_arms_screen_features_on_the_fitting_part_alone() {
        let (_, y, train, test) = null_lasso_fixture();
        let n = y.len();
        let (fit, val) = tuning_holdout(&train, &y, None, "group", LASSO_HOLDOUT, 42).unwrap();
        let mut x = BitMatrix::with_cols(n);
        x.push_row(fit.and(&y).words());
        let missed: Vec<usize> = fit.and(&y).iter_ones().take(2).collect();
        let most = train.and(&y).and_not(&Bits::from_indices(n, missed));
        x.push_row(most.words());
        let cfg = EvalConfig {
            screened_features: 1,
            ..config()
        };
        assert_eq!(screen(&x, &y, &fit, 1), [0]);
        assert_eq!(screen(&x, &y, &train, 1), [1]);
        let context = fixture_context("random");
        let arms = lasso_arms(&x, &y, &train, None, &cfg, 42, context, &test).unwrap();
        assert_eq!(
            assert_points_redo_production(&x, &y, &train, &cfg, &test, &arms),
            [0]
        );
        let on_train = lasso::fit_path(&x, &y, &fit, &[1], cfg.lasso_lambdas, 400);
        let losses = |path: &[lasso::LassoModel]| -> Vec<f64> {
            let loss = |m: &lasso::LassoModel| lasso::log_loss(&m.predict_proba(&x), &y, &val);
            path.iter().map(loss).collect()
        };
        let on_fit = lasso::fit_path(&x, &y, &fit, &[0], cfg.lasso_lambdas, 400);
        assert_ne!(losses(&on_train), losses(&on_fit));
    }

    /// The best-round arm keeps the first `best_round` trees: its scores are
    /// those of a production fit with no patience limit, which keeps that
    /// many, and one tree more would score otherwise. Feature 0 marks most
    /// resistant genomes and four rows are noise (seed 4, chosen so that
    /// the best round comes before the last), so later rounds fit noise.
    #[test]
    fn the_best_round_arm_keeps_best_round_trees_and_scores_with_them() {
        use rand::{Rng, SeedableRng};
        let n = 80;
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let all = Bits::ones(n);
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(4);
        let mut x = BitMatrix::with_cols(n);
        x.push_row(Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0 && i % 10 != 0)).words());
        for _ in 0..4 {
            let present: Vec<usize> = (0..n).filter(|_| rng.random_bool(0.5)).collect();
            x.push_row(Bits::from_indices(n, present).words());
        }
        let cfg = config();
        let test: Vec<usize> = (0..n).step_by(7).collect();
        let arms = gbt_arms(&x, &y, &all, None, "group", &cfg, 42, 42, &test).unwrap();
        let best = arms.best_round;
        assert!(0 < best && best < cfg.gbt.rounds, "best round {best}");
        let (fit, val) = tuning_holdout(&all, &y, None, "group", GBT_HOLDOUT, 42).unwrap();
        let feats = screen(&x, &y, &fit, cfg.screened_features);
        let unlimited = GbtParams {
            seed: 42,
            patience: cfg.gbt.rounds,
            ..cfg.gbt
        };
        let obj = Objective::Logistic { y: &y };
        let kept = gbt::fit(&x, &obj, &fit, &val, &feats, &unlimited);
        assert_eq!(kept.trees.len(), best);
        let scores = pick(&kept.predict_proba(&x), &test);
        assert_eq!(arms.best_round_test_scores.as_deref(), Some(&scores[..]));
        let (full, _) = gbt::fit_traced(&x, &obj, &fit, &val, &feats, &unlimited, false);
        let one_more = gbt::Gbt {
            trees: full.trees[..=best].to_vec(),
            ..full.clone()
        };
        assert_ne!(pick(&one_more.predict_proba(&x), &test), scores);
    }
}
