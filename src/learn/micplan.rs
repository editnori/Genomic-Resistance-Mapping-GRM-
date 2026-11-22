//! Experiment, off by default (`grm learn --exp-mic`): MIC models scored on a
//! fold plan the caller supplies, in place of the MIC regression's own folds
//! (docs/mic-plan.md).
//!
//! The released model's MIC regression (`eval::evaluate_mic`) draws its own
//! stratified folds over the genomes with an MIC, so a copy or a relative of
//! a test genome may train its fold. Here the folds come from a plan: the
//! lineage split `grm learn --seed` draws, or a fold file. Two arms can be
//! fitted on it, and a third is always computed:
//!
//! - `released`: `evaluate_mic`'s computation, fold for fold, on the plan's
//!   folds: candidates, holdout, screen and censored boosting as there. On
//!   the folds `evaluate_mic` draws it gives its predictions bit for bit.
//! - `ordered`: one logistic boosted fit per dilution threshold the training
//!   genomes hold, a censored value entering only the thresholds it decides,
//!   each probability put back on its training prior, made non-increasing by
//!   pool-adjacent-violators and read off at one half.
//! - `median`: the lower median of the training genomes' exact MICs, the
//!   constant a model has to beat.
//!
//! Each fold also gets the breakpoint `infer_breakpoint` fits on its training
//! genomes. Every fit and every breakpoint reads only the fold's training
//! genomes with an MIC; a test genome is only scored. No arm is compared here
//! and nothing is written to a dataset report.

use super::eval::{EvalConfig, Role, draw_seed};
use super::gbt::{self, GbtParams, Objective};
use super::{holdout, screen, screen_regression, training_candidates};
use crate::bits::{BitMatrix, Bits};
use crate::curate::{Breakpoint, Call, Censor, Mic, infer_breakpoint};
use serde::Serialize;

/// The `format` of every result `grm learn --exp-mic` writes.
pub const FORMAT: &str = "grm-mic-plan-result-v1";

/// The release's gate for an MIC model (`evaluate_mic`): at least this many
/// genomes with an MIC, and `MIN_EXACT` of them exact.
pub const MIN_MIC: usize = 100;
pub const MIN_EXACT: usize = 50;

/// Seed offset between the ordered model's thresholds: threshold j's
/// holdout and boosting draws are seeded `seed + fold + 1000 * (j + 1)`.
pub const THRESHOLD_SEED_STEP: u64 = 1000;

/// A model `grm learn --exp-mic` can fit.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Arm {
    Released,
    Ordered,
}

/// One fold of a plan as the MIC arms see it.
#[derive(Clone, Debug, Default, PartialEq, Serialize)]
pub struct FoldFit {
    /// Training genomes with an MIC, and the exact ones among them.
    pub train_mic: usize,
    pub train_exact: usize,
    /// Test genomes with an MIC, the exact ones, and the resistant ones.
    pub test_mic: usize,
    pub test_exact: usize,
    pub test_resistant: usize,
    /// `infer_breakpoint` on the training genomes' MICs and reported calls.
    pub breakpoint: Option<Breakpoint>,
    /// The training median, when the training genomes hold an exact MIC.
    pub median: Option<i32>,
    /// Trees the released arm's fit kept.
    pub released_trees: Option<usize>,
    pub ordered: Option<OrderedFold>,
}

/// The ordered model in one fold.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize)]
pub struct OrderedFold {
    /// Least and greatest step of the training genomes.
    pub lo: i32,
    pub hi: i32,
    /// hi - lo thresholds, each fitted, held constant at its prior (one
    /// class, or a holdout part without one), or decided by no training
    /// genome.
    pub thresholds: usize,
    pub fitted: usize,
    pub constant: usize,
    pub undecided: usize,
    /// Fitted thresholds whose model kept no tree, so q is the prior.
    pub null_fits: usize,
    /// Test genomes whose probabilities pool-adjacent-violators changed.
    pub pav_changed: usize,
}

/// One test genome with an MIC, in sample order, with no genome id. It
/// stays on the machine that ran the fit; only counts leave it.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct TestedMic {
    pub fold: usize,
    pub log2: i32,
    pub censor: Censor,
    pub resistant: bool,
    /// The released arm's margin, a predicted log2 MIC.
    pub released: Option<f64>,
    /// The ordered arm's predicted step.
    pub ordered: Option<i32>,
    /// The fold's training median.
    pub median: Option<i32>,
}

/// The arms on one plan. Per-sample vectors are in sample order; None for a
/// sample no fold tests, without an MIC, or of an arm that is unavailable.
#[derive(Clone, Debug, PartialEq)]
pub struct PlanFit {
    pub arms: Vec<Arm>,
    /// Why no arm was fitted (the MIC gate, or the plan's own fault).
    pub unavailable: Option<String>,
    pub released_unavailable: Option<String>,
    pub ordered_unavailable: Option<String>,
    pub median_unavailable: Option<String>,
    /// Genomes with an MIC that the plan tests or trains, and exact ones.
    pub n_mic: usize,
    pub n_exact: usize,
    pub folds: Vec<FoldFit>,
    pub fold_of: Vec<Option<usize>>,
    pub released: Vec<Option<f64>>,
    pub ordered: Vec<Option<i32>>,
    pub median: Vec<Option<i32>>,
}

impl PlanFit {
    /// Each tested genome with an MIC, in sample order.
    pub fn tested(&self, mic: &[Option<Mic>], y: &Bits) -> Vec<TestedMic> {
        (0..mic.len())
            .filter_map(|i| {
                let (fold, m) = (self.fold_of[i]?, mic[i]?);
                Some(TestedMic {
                    fold,
                    log2: m.log2,
                    censor: m.censor,
                    resistant: y.contains(i),
                    released: self.released[i],
                    ordered: self.ordered[i],
                    median: self.median[i],
                })
            })
            .collect()
    }
}

/// The lineage plan at `cfg.seed`: the outer folds of the lineage split
/// `evaluate` draws (`grouped_folds` over the cohort's lineages), every
/// sample tested in one fold.
pub fn lineage_roles(y: &Bits, groups: &[u32], cfg: &EvalConfig) -> Vec<Role> {
    super::grouped_folds(y, groups, cfg.folds, cfg.seed)
        .into_iter()
        .map(Role::Test)
        .collect()
}

/// None when every fold's test genomes and training genomes hold both
/// classes, else the first fold that does not. A fold plan is not refused
/// for it: the MIC arms read the calls only through the breakpoint.
pub fn class_check(y: &Bits, roles: &[Role], folds: usize) -> Option<String> {
    (0..folds).find_map(|f| {
        let (mut test, mut train) = ([0usize; 2], [0usize; 2]);
        for (i, r) in roles.iter().enumerate() {
            let side = match r {
                Role::Test(g) if *g == f => &mut test,
                Role::Test(_) | Role::Train => &mut train,
                Role::Unused => continue,
            };
            side[y.contains(i) as usize] += 1;
        }
        [("test", test), ("training", train)]
            .into_iter()
            .find(|(_, c)| c[0] == 0 || c[1] == 0)
            .map(|(name, _)| format!("fold {f}'s {name} genomes lack one or both classes"))
    })
}

/// A plan no arm is fitted on, for `why`.
pub fn refused(why: String, arms: &[Arm], mic: &[Option<Mic>], roles: &[Role]) -> PlanFit {
    let n = mic.len();
    let used = |i: &usize| mic[*i].is_some() && roles[*i] != Role::Unused;
    let mut arms = arms.to_vec();
    arms.sort();
    arms.dedup();
    PlanFit {
        arms,
        unavailable: Some(why),
        released_unavailable: None,
        ordered_unavailable: None,
        median_unavailable: None,
        n_mic: (0..n).filter(used).count(),
        n_exact: (0..n)
            .filter(used)
            .filter(|&i| mic[i].unwrap().censor == Censor::Exact)
            .count(),
        folds: Vec::new(),
        fold_of: vec![None; n],
        released: vec![None; n],
        ordered: vec![None; n],
        median: vec![None; n],
    }
}

/// A genome's dilution step: v for an exact MIC or one at most v, v + 1 for
/// one above v.
pub fn step(m: Mic) -> i32 {
    match m.censor {
        Censor::Exact | Censor::AtMost => m.log2,
        Censor::Above => m.log2 + 1,
    }
}

/// Whether `m` decides "MIC above t", and how: Some(true) when it is above,
/// Some(false) when not, None when its censoring leaves it open.
pub fn decides(m: Mic, t: i32) -> Option<bool> {
    match m.censor {
        Censor::Exact => Some(m.log2 > t),
        Censor::AtMost => (t >= m.log2).then_some(false),
        Censor::Above => (t <= m.log2).then_some(true),
    }
}

/// Pool-adjacent-violators with equal weights: the non-increasing sequence
/// nearest to `q` in squared error, and whether it differs from `q`.
pub fn non_increasing(q: &[f64]) -> (Vec<f64>, bool) {
    let mut blocks: Vec<(f64, usize)> = Vec::with_capacity(q.len());
    let mut merged = false;
    for &v in q {
        blocks.push((v, 1));
        while blocks.len() > 1 {
            let (s1, c1) = blocks[blocks.len() - 1];
            let (s0, c0) = blocks[blocks.len() - 2];
            if s0 / c0 as f64 >= s1 / c1 as f64 {
                break;
            }
            blocks.pop();
            *blocks.last_mut().unwrap() = (s0 + s1, c0 + c1);
            merged = true;
        }
    }
    let out = blocks
        .into_iter()
        .flat_map(|(s, c)| std::iter::repeat_n(s / c as f64, c))
        .collect();
    (out, merged)
}

/// The ordered model's predicted step for one genome: lo plus the number of
/// thresholds whose monotone probability of "above" is one half or more,
/// and whether making it monotone changed `q`.
pub fn read_off(lo: i32, q: &[f64]) -> (i32, bool) {
    let (m, changed) = non_increasing(q);
    (lo + m.iter().filter(|&&p| p >= 0.5).count() as i32, changed)
}

/// The lower median of `values`: with m values sorted, the one at position
/// ceil(m / 2), counting from 1 (lane micbase's B2).
pub fn lower_median(values: &mut [i32]) -> Option<i32> {
    values.sort_unstable();
    (!values.is_empty()).then(|| values[values.len().div_ceil(2) - 1])
}

fn logit(p: f64) -> f64 {
    (p / (1.0 - p)).ln()
}

/// What a threshold of the ordered model gives a fold's test genomes.
enum Threshold {
    /// No training genome decides it.
    Undecided,
    /// Held at the deciding genomes' share above it.
    Constant(f64),
    /// Fitted: each sample's probability, prior corrected.
    Fitted(Vec<f64>),
}

/// Fits the arms on the plan `roles` gives over `folds` outer folds. `x` is
/// the dataset's projection, `y` the reported resistant calls, `mic` each
/// sample's MIC, and `min_support` the floor candidates are chosen at. Only
/// samples with an MIC reach a fit: in fold f, those tested in another fold
/// or training in all of them train, and those tested in f are scored. A
/// fold with no MIC test genome is fitted by no arm.
#[allow(clippy::too_many_arguments)]
pub fn fit_plan(
    x: &BitMatrix,
    y: &Bits,
    mic: &[Option<Mic>],
    roles: &[Role],
    folds: usize,
    min_support: u32,
    cfg: &EvalConfig,
    arms: &[Arm],
) -> PlanFit {
    let n = x.cols;
    assert_eq!(y.len(), n);
    assert_eq!(mic.len(), n);
    assert_eq!(roles.len(), n);
    let mut arms = arms.to_vec();
    arms.sort();
    arms.dedup();
    let used = |i: usize| mic[i].is_some() && roles[i] != Role::Unused;
    let n_mic = (0..n).filter(|&i| used(i)).count();
    let n_exact = (0..n)
        .filter(|&i| used(i) && mic[i].unwrap().censor == Censor::Exact)
        .count();
    let mut out = PlanFit {
        arms: arms.clone(),
        unavailable: None,
        released_unavailable: None,
        ordered_unavailable: None,
        median_unavailable: None,
        n_mic,
        n_exact,
        folds: Vec::new(),
        fold_of: vec![None; n],
        released: vec![None; n],
        ordered: vec![None; n],
        median: vec![None; n],
    };
    if n_mic < MIN_MIC || n_exact < MIN_EXACT {
        out.unavailable = Some(format!(
            "{n_mic} genomes with an MIC and {n_exact} exact, under the gate of {MIN_MIC} and {MIN_EXACT}"
        ));
        return out;
    }
    let value: Vec<f64> = (0..n)
        .map(|i| mic[i].map_or(0.0, |v| v.log2 as f64))
        .collect();
    let censor: Vec<Censor> = (0..n)
        .map(|i| mic[i].map_or(Censor::Exact, |v| v.censor))
        .collect();
    for f in 0..folds {
        let test = Bits::from_indices(
            n,
            (0..n).filter(|&i| mic[i].is_some() && roles[i] == Role::Test(f)),
        );
        let train = Bits::from_indices(
            n,
            (0..n).filter(|&i| {
                mic[i].is_some()
                    && match roles[i] {
                        Role::Test(g) => g != f,
                        Role::Train => true,
                        Role::Unused => false,
                    }
            }),
        );
        let exact = |m: &Bits| {
            m.iter_ones()
                .filter(|&i| censor[i] == Censor::Exact)
                .count()
        };
        let mut fold = FoldFit {
            train_mic: train.count() as usize,
            train_exact: exact(&train),
            test_mic: test.count() as usize,
            test_exact: exact(&test),
            test_resistant: test.and(y).count() as usize,
            ..Default::default()
        };
        if test.count() == 0 {
            out.folds.push(fold);
            continue;
        }
        for i in test.iter_ones() {
            out.fold_of[i] = Some(f);
        }
        let pairs: Vec<(Mic, Call)> = train
            .iter_ones()
            .map(|i| {
                (
                    mic[i].unwrap(),
                    if y.contains(i) { Call::R } else { Call::S },
                )
            })
            .collect();
        fold.breakpoint = infer_breakpoint(&pairs);
        let mut exact_train: Vec<i32> = train
            .iter_ones()
            .filter(|&i| censor[i] == Censor::Exact)
            .map(|i| mic[i].unwrap().log2)
            .collect();
        fold.median = lower_median(&mut exact_train);
        match fold.median {
            Some(c) => test.iter_ones().for_each(|i| out.median[i] = Some(c)),
            None => {
                out.median_unavailable.get_or_insert_with(|| {
                    format!("fold {f}'s training genomes hold no exact MIC")
                });
            }
        }
        if train.count() == 0 {
            let why = format!("fold {f} tests genomes with an MIC and trains on none");
            out.released_unavailable.get_or_insert_with(|| why.clone());
            out.ordered_unavailable.get_or_insert(why);
            out.folds.push(fold);
            continue;
        }
        let (cx, _) = training_candidates(x, &train, min_support);
        if arms.contains(&Arm::Released) && out.released_unavailable.is_none() {
            let (fit_m, val) = holdout(&train, y, 0.15, cfg.seed.wrapping_add(f as u64));
            if fit_m.count() == 0 || val.count() == 0 {
                out.released_unavailable = Some(format!(
                    "fold {f}'s training genomes give an empty fitting or holdout part"
                ));
            } else {
                let feats = screen_regression(&cx, &value, &fit_m, cfg.screened_features);
                let model = gbt::fit(
                    &cx,
                    &Objective::Censored {
                        value: &value,
                        censor: &censor,
                    },
                    &fit_m,
                    &val,
                    &feats,
                    &GbtParams {
                        seed: draw_seed(cfg).wrapping_add(f as u64),
                        tie: cfg.tie_rule,
                        ..cfg.gbt
                    },
                );
                fold.released_trees = Some(model.trees.len());
                let pred = model.margin(&cx);
                for i in test.iter_ones() {
                    out.released[i] = Some(pred[i]);
                }
            }
        }
        if arms.contains(&Arm::Ordered) && out.ordered_unavailable.is_none() {
            match fit_ordered_fold(&cx, mic, &train, &test, f, cfg) {
                Ok((summary, steps)) => {
                    for (i, s) in test.iter_ones().zip(steps) {
                        out.ordered[i] = Some(s);
                    }
                    fold.ordered = Some(summary);
                }
                Err(why) => out.ordered_unavailable = Some(why),
            }
        }
        out.folds.push(fold);
    }
    let tested: Vec<usize> = (0..n).filter(|&i| out.fold_of[i].is_some()).collect();
    if arms.contains(&Arm::Released)
        && out.released_unavailable.is_none()
        && tested
            .iter()
            .any(|&i| !out.released[i].is_some_and(f64::is_finite))
    {
        out.released_unavailable = Some("a prediction is not a finite number".into());
    }
    if tested.is_empty() {
        out.unavailable = Some("the plan tests no genome with an MIC".into());
    }
    if out.released_unavailable.is_some() {
        out.released.iter_mut().for_each(|p| *p = None);
    }
    if out.median_unavailable.is_some() {
        out.median.iter_mut().for_each(|p| *p = None);
    }
    if out.ordered_unavailable.is_some() {
        out.ordered.iter_mut().for_each(|p| *p = None);
    }
    out
}

/// The ordered model in fold `f`: its summary and each test genome's
/// predicted step, in sample order.
fn fit_ordered_fold(
    cx: &BitMatrix,
    mic: &[Option<Mic>],
    train: &Bits,
    test: &Bits,
    f: usize,
    cfg: &EvalConfig,
) -> Result<(OrderedFold, Vec<i32>), String> {
    let n = cx.cols;
    let steps: Vec<i32> = train.iter_ones().map(|i| step(mic[i].unwrap())).collect();
    let (lo, hi) = (*steps.iter().min().unwrap(), *steps.iter().max().unwrap());
    let mut summary = OrderedFold {
        lo,
        hi,
        thresholds: (hi - lo) as usize,
        ..Default::default()
    };
    let mut fits = Vec::with_capacity(summary.thresholds);
    for (j, t) in (lo..hi).enumerate() {
        let calls: Vec<(usize, bool)> = train
            .iter_ones()
            .filter_map(|i| decides(mic[i].unwrap(), t).map(|a| (i, a)))
            .collect();
        if calls.is_empty() {
            summary.undecided += 1;
            fits.push(Threshold::Undecided);
            continue;
        }
        let deciding = Bits::from_indices(n, calls.iter().map(|&(i, _)| i));
        let target = Bits::from_indices(n, calls.iter().filter(|c| c.1).map(|&(i, _)| i));
        let above = target.count() as usize;
        let share = above as f64 / calls.len() as f64;
        if above == 0 || above == calls.len() {
            summary.constant += 1;
            fits.push(Threshold::Constant(share));
            continue;
        }
        let seed = cfg
            .seed
            .wrapping_add(f as u64)
            .wrapping_add(THRESHOLD_SEED_STEP.wrapping_mul(j as u64 + 1));
        let (fit_m, val) = holdout(&deciding, &target, 0.15, seed);
        let both = |m: &Bits| {
            let a = m.and(&target).count();
            a > 0 && a < m.count()
        };
        if !both(&fit_m) || !both(&val) {
            summary.constant += 1;
            fits.push(Threshold::Constant(share));
            continue;
        }
        let prior = fit_m.and(&target).count() as f64 / fit_m.count() as f64;
        let feats = screen(cx, &target, &fit_m, cfg.screened_features);
        let model = gbt::fit(
            cx,
            &Objective::Logistic { y: &target },
            &fit_m,
            &val,
            &feats,
            &GbtParams {
                seed: draw_seed(cfg)
                    .wrapping_add(f as u64)
                    .wrapping_add(THRESHOLD_SEED_STEP.wrapping_mul(j as u64 + 1)),
                tie: cfg.tie_rule,
                ..cfg.gbt
            },
        );
        summary.fitted += 1;
        let null = model.trees.is_empty();
        summary.null_fits += null as usize;
        let q = if null {
            vec![prior; n]
        } else {
            let shift = logit(prior);
            model
                .margin(cx)
                .into_iter()
                .map(|m| gbt::sigmoid(m + shift))
                .collect()
        };
        fits.push(Threshold::Fitted(q));
    }
    if !fits.is_empty() && fits.iter().all(|t| matches!(t, Threshold::Undecided)) {
        return Err(format!(
            "no threshold of fold {f} is decided by a training genome"
        ));
    }
    // An undecided threshold takes the nearest lower decided one, failing
    // that the nearest higher.
    let source: Vec<usize> = (0..fits.len())
        .map(|j| {
            let decided = |k: &usize| !matches!(fits[*k], Threshold::Undecided);
            (0..=j)
                .rev()
                .find(decided)
                .or_else(|| (j..fits.len()).find(decided))
                .unwrap()
        })
        .collect();
    let mut predicted = Vec::with_capacity(test.count() as usize);
    for i in test.iter_ones() {
        let q: Vec<f64> = source
            .iter()
            .map(|&k| match &fits[k] {
                Threshold::Constant(p) => *p,
                Threshold::Fitted(q) => q[i],
                Threshold::Undecided => unreachable!(),
            })
            .collect();
        let (s, changed) = read_off(lo, &q);
        summary.pav_changed += changed as usize;
        predicted.push(s);
    }
    Ok((summary, predicted))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::learn::grouped_folds;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    fn mic(log2: i32, censor: Censor) -> Mic {
        Mic { log2, censor }
    }

    /// A synthetic dataset of `n` samples: 24 random features, an MIC that
    /// rises two steps with each of the first two, censored at both ends of
    /// the range, the reported call resistant above 1, and 40 lineages.
    fn fixture(n: usize) -> (BitMatrix, Bits, Vec<Option<Mic>>, Vec<u32>) {
        let mut rng = ChaCha8Rng::seed_from_u64(7);
        let rows: Vec<Vec<bool>> = (0..24)
            .map(|_| (0..n).map(|_| rng.random_bool(0.4)).collect())
            .collect();
        let mut x = BitMatrix::with_cols(n);
        for r in &rows {
            x.push_row(Bits::from_indices(n, (0..n).filter(|&i| r[i])).words());
        }
        let mut mics = Vec::with_capacity(n);
        for (i, (&a, &b)) in rows[0].iter().zip(&rows[1]).enumerate() {
            let v = -2 + 2 * a as i32 + 3 * b as i32 + rng.random_range(-1..=1);
            let m = if v <= -2 {
                mic(-2, Censor::AtMost)
            } else if v >= 4 {
                mic(3, Censor::Above)
            } else {
                mic(v, Censor::Exact)
            };
            mics.push((i % 9 != 0).then_some(m));
        }
        let y = Bits::from_indices(
            n,
            (0..n).filter(|&i| {
                rows[0][i] && rows[1][i]
                    || (i % 9 == 0 && rows[1][i])
                    || mics[i].is_some_and(|m| step(m) > 1)
            }),
        );
        let groups = (0..n).map(|i| (i % 40) as u32).collect();
        (x, y, mics, groups)
    }

    fn config(seed: u64) -> EvalConfig {
        EvalConfig {
            seed,
            gbt: GbtParams {
                rounds: 30,
                max_depth: 2,
                ..Default::default()
            },
            screened_features: 12,
            ..Default::default()
        }
    }

    fn both() -> Vec<Arm> {
        vec![Arm::Released, Arm::Ordered]
    }

    #[test]
    fn a_censored_value_enters_only_the_thresholds_it_decides() {
        // Exact 1 decides every threshold.
        assert_eq!(decides(mic(1, Censor::Exact), 0), Some(true));
        assert_eq!(decides(mic(1, Censor::Exact), 1), Some(false));
        assert_eq!(decides(mic(1, Censor::Exact), 5), Some(false));
        // At most 1: not above 1 or higher, open below.
        assert_eq!(decides(mic(1, Censor::AtMost), 0), None);
        assert_eq!(decides(mic(1, Censor::AtMost), 1), Some(false));
        assert_eq!(decides(mic(1, Censor::AtMost), 3), Some(false));
        // Above 1: above 1 and lower, open from 2.
        assert_eq!(decides(mic(1, Censor::Above), -3), Some(true));
        assert_eq!(decides(mic(1, Censor::Above), 1), Some(true));
        assert_eq!(decides(mic(1, Censor::Above), 2), None);
        assert_eq!(step(mic(1, Censor::Above)), 2);
        assert_eq!(step(mic(1, Censor::AtMost)), 1);
    }

    #[test]
    fn the_read_off_is_the_median_step_of_the_monotone_sequence() {
        let (m, changed) = non_increasing(&[0.9, 0.25, 0.75, 0.1]);
        assert_eq!(m, vec![0.9, 0.5, 0.5, 0.1]);
        assert!(changed);
        // Monotone after: each value is at most the one before it.
        assert!(m.windows(2).all(|w| w[0] >= w[1]));
        assert_eq!(read_off(-1, &[0.9, 0.25, 0.75, 0.1]), (2, true));
        assert_eq!(read_off(-1, &[0.9, 0.6, 0.4, 0.1]), (1, false));
        assert_eq!(read_off(3, &[]), (3, false));
        assert_eq!(read_off(0, &[0.2, 0.8]), (2, true));
    }

    #[test]
    fn the_training_median_is_the_lower_median() {
        assert_eq!(lower_median(&mut [3, 1, 2, 4]), Some(2));
        assert_eq!(lower_median(&mut [5, 1, 3]), Some(3));
        assert_eq!(lower_median(&mut []), None);
    }

    #[test]
    fn the_released_arm_on_the_released_folds_is_evaluate_mic_bit_for_bit() {
        let (x, y, mics, _) = fixture(220);
        let cfg = config(42);
        let input = super::super::eval::MicInput {
            mic: &mics,
            breakpoint: None,
        };
        let released = super::super::eval::evaluate_mic(&x, &y, 3, &input, &cfg).unwrap();
        let with_mic = Bits::from_indices(220, (0..220).filter(|&i| mics[i].is_some()));
        let folds = super::super::training_folds(&y, &with_mic, None, cfg.folds, cfg.seed);
        let roles: Vec<Role> = folds
            .iter()
            .map(|&f| {
                if f == usize::MAX {
                    Role::Unused
                } else {
                    Role::Test(f)
                }
            })
            .collect();
        let fit = fit_plan(&x, &y, &mics, &roles, cfg.folds, 3, &cfg, &[Arm::Released]);
        assert_eq!(fit.released_unavailable, None);
        for i in 0..220 {
            match fit.released[i] {
                Some(p) => assert_eq!(p.to_bits(), released.predicted[i].to_bits()),
                None => assert!(released.predicted[i].is_nan()),
            }
        }
        let bps: Vec<_> = fit.folds.iter().map(|f| f.breakpoint).collect();
        assert_eq!(bps, released.fold_breakpoints);
    }

    /// The plan fitted twice: as drawn, and with fold 0's test genomes
    /// given another MIC (exact 6) and a resistant call. Fold 0's fits,
    /// breakpoint and median read only the other folds, so none of them may
    /// change.
    fn changed_test_fold() -> (PlanFit, PlanFit, Vec<Role>) {
        let (x, y, mics, groups) = fixture(240);
        let cfg = config(42);
        let roles = lineage_roles(&y, &groups, &cfg);
        let a = fit_plan(&x, &y, &mics, &roles, 5, 3, &cfg, &both());
        let mut mics2 = mics.clone();
        let mut y2 = y.clone();
        for i in 0..240 {
            if roles[i] == Role::Test(0) && mics[i].is_some() {
                mics2[i] = Some(mic(6, Censor::Exact));
                y2.insert(i);
            }
        }
        let b = fit_plan(&x, &y2, &mics2, &roles, 5, 3, &cfg, &both());
        (a, b, roles)
    }

    #[test]
    fn fold_breakpoints_and_medians_read_no_test_genome() {
        let (a, b, _) = changed_test_fold();
        assert_eq!(a.folds[0].breakpoint, b.folds[0].breakpoint);
        assert_eq!(a.folds[0].median, b.folds[0].median);
        assert!(a.folds[0].breakpoint.is_some());
    }

    #[test]
    fn predictions_read_no_test_genome() {
        let (a, b, roles) = changed_test_fold();
        assert!(a.released_unavailable.is_none() && a.ordered_unavailable.is_none());
        for (i, r) in roles.iter().enumerate() {
            if *r == Role::Test(0) {
                assert_eq!(a.released[i], b.released[i], "sample {i}");
                assert_eq!(a.ordered[i], b.ordered[i], "sample {i}");
            }
        }
    }

    #[test]
    fn the_lineage_plan_keeps_each_lineage_in_one_fold_and_follows_the_seed() {
        let (_, y, _, groups) = fixture(240);
        let at = |s| lineage_roles(&y, &groups, &config(s));
        assert_ne!(at(42), at(43));
        for s in [42, 43] {
            let roles = at(s);
            let want: Vec<Role> = grouped_folds(&y, &groups, 5, s)
                .into_iter()
                .map(Role::Test)
                .collect();
            assert_eq!(roles, want);
            for g in 0..40u32 {
                let f: Vec<Role> = (0..240)
                    .filter(|&i| groups[i] == g)
                    .map(|i| roles[i])
                    .collect();
                assert!(f.windows(2).all(|w| w[0] == w[1]), "lineage {g}");
            }
        }
    }

    #[test]
    fn the_fits_follow_the_seed() {
        let (x, y, mics, groups) = fixture(240);
        let roles = lineage_roles(&y, &groups, &config(42));
        let fit = |s| fit_plan(&x, &y, &mics, &roles, 5, 3, &config(s), &both());
        let (a, b) = (fit(42), fit(43));
        assert_ne!(a.released, b.released);
        assert_ne!(a.ordered, b.ordered);
        assert_eq!(a.median, b.median);
        assert_eq!(fit(42), a);
    }

    #[test]
    fn ordered_predictions_are_steps_inside_the_training_range() {
        let (x, y, mics, groups) = fixture(240);
        let cfg = config(42);
        let roles = lineage_roles(&y, &groups, &cfg);
        let fit = fit_plan(&x, &y, &mics, &roles, 5, 3, &cfg, &[Arm::Ordered]);
        assert!(fit.ordered_unavailable.is_none());
        assert!(fit.released.iter().all(Option::is_none));
        for (f, fold) in fit.folds.iter().enumerate() {
            let o = fold.ordered.as_ref().unwrap();
            assert_eq!(o.thresholds, o.fitted + o.constant + o.undecided);
            assert_eq!((o.lo, o.hi), (-2, 4));
            for i in 0..240 {
                if fit.fold_of[i] == Some(f) {
                    let s = fit.ordered[i].unwrap();
                    assert!((o.lo..=o.hi).contains(&s));
                }
            }
        }
        // Most exact MICs are read within one step.
        let t = fit.tested(&mics, &y);
        let exact: Vec<_> = t.iter().filter(|r| r.censor == Censor::Exact).collect();
        let hit = exact
            .iter()
            .filter(|r| (r.ordered.unwrap() - r.log2).abs() <= 1)
            .count();
        assert!(hit * 10 >= exact.len() * 8, "{hit} of {}", exact.len());
    }

    #[test]
    fn one_dilution_step_gives_its_step_and_no_threshold() {
        let (x, y, _, groups) = fixture(240);
        let mics: Vec<Option<Mic>> = (0..240).map(|_| Some(mic(2, Censor::Exact))).collect();
        let cfg = config(42);
        let roles = lineage_roles(&y, &groups, &cfg);
        let fit = fit_plan(&x, &y, &mics, &roles, 5, 3, &cfg, &[Arm::Ordered]);
        assert!(fit.ordered.iter().all(|&p| p == Some(2)));
        assert!(
            fit.folds
                .iter()
                .all(|f| f.ordered.as_ref().unwrap().thresholds == 0)
        );
        assert!(fit.median.iter().all(|&p| p == Some(2)));
    }

    #[test]
    fn a_fold_without_an_mic_test_genome_is_not_fitted_and_the_gate_holds() {
        let (x, y, mut mics, groups) = fixture(240);
        let cfg = config(42);
        let roles = lineage_roles(&y, &groups, &cfg);
        for i in 0..240 {
            if roles[i] == Role::Test(4) {
                mics[i] = None;
            }
        }
        let fit = fit_plan(&x, &y, &mics, &roles, 5, 3, &cfg, &both());
        assert_eq!(fit.folds[4].test_mic, 0);
        assert!(fit.folds[4].breakpoint.is_none() && fit.folds[4].ordered.is_none());
        assert!(fit.released_unavailable.is_none());
        let few: Vec<Option<Mic>> = (0..240)
            .map(|i| {
                if i < 99 {
                    mics[i].or(Some(mic(1, Censor::Exact)))
                } else {
                    None
                }
            })
            .collect();
        let gated = fit_plan(&x, &y, &few, &roles, 5, 3, &cfg, &both());
        assert!(gated.unavailable.is_some());
        assert!(gated.released.iter().all(Option::is_none));
    }
}
