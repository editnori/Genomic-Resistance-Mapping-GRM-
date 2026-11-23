//! Experiment, off by default (plan block C4, docs/screen-arms.md): the
//! screen and the floor, as arms of `grm learn`.
//!
//! `grm learn --exp-screen cmh` ranks each fit's candidates for the forest,
//! boosting and the lasso by the Cochran-Mantel-Haenszel statistic over
//! lineage strata in place of the chi-square; `--screened-features` sets its
//! budget as it sets the default's. `grm learn --exp-floor N` fixes the
//! per-fit support floor at N in place of max(3, n/200). Without either flag
//! `active` is None, no context is installed, `screen_fit` is
//! `learn::screen`, and `grm learn` writes what it wrote before, byte for
//! byte.
//!
//! The strata are the release's lineage clusters, the `groups` `evaluate`
//! and `evaluate_assigned` take, restricted to each fit's own fitting
//! samples: the screen reads no other sample. A stratum among them with one
//! sample or one class adds nothing to the statistic, since its
//! hypergeometric variance is 0; a feature with no informative stratum, or
//! constant inside each one, gets 0. Candidates rank by the statistic, high
//! first, ties by the lower row, as `learn::screen` ranks them.
//!
//! With either flag on, `evaluate` and `evaluate_assigned` install a context
//! on their own thread (every fit there runs on it) and write
//! `<dir>/<dataset>.<tag>-s<seed>.inputs.json` at the end: counts and fit
//! times only, no genome id, feature, k-mer or position.

use crate::bits::{BitMatrix, Bits};
use anyhow::{Context, Result, ensure};
use rayon::prelude::*;
use serde::Serialize;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::time::Instant;

/// The result format. A change to what a result holds changes it.
pub const FORMAT: &str = "grm-c4-screen-1";

/// How the forest, boosting and the lasso rank a fit's candidates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Rule {
    /// The release's 2x2 Pearson chi-square (`learn::screen`).
    ChiSquare,
    /// The Cochran-Mantel-Haenszel statistic over lineage strata.
    Cmh,
}

impl Rule {
    pub fn name(self) -> &'static str {
        match self {
            Rule::ChiSquare => "chi-square",
            Rule::Cmh => "cmh",
        }
    }
}

/// The experiment's settings, as `enable` set them.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Settings {
    pub rule: Rule,
    /// The per-fit support floor in place of max(3, n/200); None keeps it.
    pub floor: Option<u32>,
    /// Where each inputs file goes: the `--exp-pick` folder.
    pub out: PathBuf,
}

impl Settings {
    /// The arm's name in its inputs file: `cmh` or `floor`.
    pub fn arm(&self) -> &'static str {
        if self.rule == Rule::Cmh { "cmh" } else { "floor" }
    }
}

static SINK: OnceLock<Settings> = OnceLock::new();

#[cfg(test)]
thread_local! {
    static TEST_SETTINGS: RefCell<Option<Settings>> = const { RefCell::new(None) };
}

/// Turns the experiment on for this process. A second call fails, and so
/// does a setting that changes nothing or both things.
pub fn enable(settings: Settings) -> Result<()> {
    ensure!(
        (settings.rule == Rule::Cmh) != settings.floor.is_some(),
        "an arm of C4 changes the screen or the floor, one of them"
    );
    std::fs::create_dir_all(&settings.out)
        .with_context(|| format!("create {}", settings.out.display()))?;
    SINK.set(settings)
        .map_err(|_| anyhow::anyhow!("the C4 screen arm is already set"))
}

/// The arm this process runs, if any.
pub fn active() -> Option<Settings> {
    #[cfg(test)]
    if let Some(s) = TEST_SETTINGS.with(|t| t.borrow().clone()) {
        return Some(s);
    }
    SINK.get().cloned()
}

/// Runs `f` with `settings` as this thread's arm: tests only.
#[cfg(test)]
pub(crate) fn with_settings<T>(settings: Settings, f: impl FnOnce() -> T) -> T {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            TEST_SETTINGS.with(|t| *t.borrow_mut() = None);
        }
    }
    TEST_SETTINGS.with(|t| *t.borrow_mut() = Some(settings));
    let _reset = Reset;
    f()
}

/// The floor the active arm fixes, if any (`eval::support_floor`).
pub fn floor() -> Option<u32> {
    active().and_then(|s| s.floor)
}

/// One screen's ranking and counts.
#[derive(Clone, Debug, PartialEq)]
pub struct Screened {
    /// The kept candidate rows, best first.
    pub feats: Vec<u32>,
    /// Kept rows whose statistic is 0.
    pub kept_at_zero: usize,
    /// Strata among the fitting samples.
    pub strata: usize,
    /// Of those, the ones with two or more samples and both classes.
    pub informative: usize,
}

/// The strata of the fitting samples `fit`: for each sample, the index of
/// its informative stratum (two or more fitting samples, both classes) in
/// ascending order of the lineage id, or None; and the (n, r) of each
/// informative stratum, with the count of every stratum among `fit`.
/// `lineage` is read for fitting samples only; `y` too.
pub fn fit_strata(
    y: &Bits,
    fit: &Bits,
    lineage: &[u32],
) -> (Vec<Option<u32>>, Vec<(u64, u64)>, usize) {
    assert_eq!(lineage.len(), fit.len());
    let mut by: BTreeMap<u32, (u64, u64)> = BTreeMap::new();
    for i in fit.iter_ones() {
        let e = by.entry(lineage[i]).or_default();
        e.0 += 1;
        e.1 += y.contains(i) as u64;
    }
    let all = by.len();
    let mut index: BTreeMap<u32, u32> = BTreeMap::new();
    let mut nr = Vec::new();
    for (&id, &(n, r)) in &by {
        if n >= 2 && r > 0 && r < n {
            index.insert(id, nr.len() as u32);
            nr.push((n, r));
        }
    }
    let mut slot = vec![None; fit.len()];
    for i in fit.iter_ones() {
        slot[i] = index.get(&lineage[i]).copied();
    }
    (slot, nr, all)
}

/// The Cochran-Mantel-Haenszel statistic of every row of `x` on the
/// fitting samples `fit`, over the informative strata of `lineage`, with no
/// continuity correction: D² / V, D = Σ (a − m r / n) and
/// V = Σ m (n − m) r (n − r) / (n² (n − 1)) in ascending stratum order, and
/// 0 when V is 0. With one stratum it is (n − 1)/n times `learn::screen`'s
/// chi-square.
pub fn cmh_statistics(x: &BitMatrix, y: &Bits, fit: &Bits, lineage: &[u32]) -> (Vec<f64>, usize, usize) {
    assert_eq!(x.cols, fit.len());
    let (slot, nr, all) = fit_strata(y, fit, lineage);
    let informative = nr.len();
    // Fitting samples in an informative stratum: the only ones that count.
    let inside = Bits::from_indices(fit.len(), (0..fit.len()).filter(|&i| slot[i].is_some()));
    let inside_y = inside.and(y);
    let stats = (0..x.rows)
        .into_par_iter()
        .map_init(
            || (vec![0u64; informative], vec![0u64; informative]),
            |(m, a), f| {
                m.iter_mut().for_each(|v| *v = 0);
                a.iter_mut().for_each(|v| *v = 0);
                let row = x.row(f);
                for (wi, (&r, (&w, &wy))) in row
                    .iter()
                    .zip(inside.words().iter().zip(inside_y.words()))
                    .enumerate()
                {
                    let mut present = r & w;
                    while present != 0 {
                        let i = wi * 64 + present.trailing_zeros() as usize;
                        m[slot[i].expect("an inside sample has a stratum") as usize] += 1;
                        present &= present - 1;
                    }
                    let mut resistant = r & wy;
                    while resistant != 0 {
                        let i = wi * 64 + resistant.trailing_zeros() as usize;
                        a[slot[i].expect("an inside sample has a stratum") as usize] += 1;
                        resistant &= resistant - 1;
                    }
                }
                cmh_from_counts(&nr, m, a)
            },
        )
        .collect();
    (stats, all, informative)
}

/// The statistic from each informative stratum's (n, r) and the feature's
/// present (m) and present-and-resistant (a) counts, in stratum order.
pub fn cmh_from_counts(nr: &[(u64, u64)], m: &[u64], a: &[u64]) -> f64 {
    let (mut d, mut v) = (0.0f64, 0.0f64);
    for ((&(n, r), &m), &a) in nr.iter().zip(m).zip(a) {
        d += a as f64 - (m * r) as f64 / n as f64;
        v += (m * (n - m)) as f64 * (r * (n - r)) as f64 / (n * n * (n - 1)) as f64;
    }
    if v > 0.0 { d * d / v } else { 0.0 }
}

/// The top `k` rows of `x` by the Cochran-Mantel-Haenszel statistic on the
/// fitting samples `fit` over the lineage strata, ties by the lower row, as
/// `learn::screen` ranks by the chi-square.
pub fn cmh_screen(x: &BitMatrix, y: &Bits, fit: &Bits, lineage: &[u32], k: usize) -> Screened {
    let (stats, strata, informative) = cmh_statistics(x, y, fit, lineage);
    let mut scored: Vec<(f64, u32)> = stats.into_iter().zip(0u32..).collect();
    let k = k.min(scored.len());
    if k < scored.len() {
        scored.select_nth_unstable_by(k, |a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        scored.truncate(k);
    }
    scored.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
    Screened {
        kept_at_zero: scored.iter().filter(|s| s.0 == 0.0).count(),
        feats: scored.into_iter().map(|(_, f)| f).collect(),
        strata,
        informative,
    }
}

// ---------- the context of one dataset's fits ----------

/// One outer fold's candidates.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct FoldNote {
    pub split: String,
    pub fold: usize,
    pub train: usize,
    /// Candidates past the floor (`training_candidates`).
    pub candidates: usize,
}

/// One screen of one fit.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct ScreenNote {
    pub split: String,
    /// None for the final models.
    pub fold: Option<usize>,
    pub learner: String,
    pub rule: &'static str,
    /// The fit's fitting samples.
    pub fitting: usize,
    /// The candidates the screen ranked.
    pub candidates: usize,
    /// The columns it kept.
    pub kept: usize,
    /// Of those, the ones at a statistic of 0 (the adjusted screen only).
    pub kept_at_zero: Option<usize>,
    /// Strata among the fitting samples, and the informative ones (the
    /// adjusted screen only).
    pub strata: Option<usize>,
    pub informative_strata: Option<usize>,
}

/// One outer fit's time: wall and CPU seconds, and the process's peak
/// resident size after it.
#[derive(Clone, Debug, Serialize, PartialEq)]
pub struct FitTime {
    pub split: String,
    pub fold: usize,
    pub learner: String,
    pub seconds: f64,
    pub cpu_seconds: f64,
    pub peak_kb: u64,
}

struct Context {
    rule: Rule,
    lineage: Vec<u32>,
    label: (String, Option<usize>, String),
    folds: Vec<FoldNote>,
    screens: Vec<ScreenNote>,
    times: Vec<FitTime>,
}

thread_local! {
    static CONTEXT: RefCell<Option<Context>> = const { RefCell::new(None) };
}

/// Removes this thread's context when dropped.
pub struct Installed(());

impl Drop for Installed {
    fn drop(&mut self) {
        CONTEXT.with(|c| *c.borrow_mut() = None);
    }
}

/// Installs a context for one dataset's fits on this thread when an arm is
/// active, with each sample's lineage; None otherwise, and then nothing
/// here is called.
pub fn install(lineage: &[u32]) -> Option<Installed> {
    let settings = active()?;
    CONTEXT.with(|c| {
        *c.borrow_mut() = Some(Context {
            rule: settings.rule,
            lineage: lineage.to_vec(),
            label: (String::new(), None, String::new()),
            folds: Vec::new(),
            screens: Vec::new(),
            times: Vec::new(),
        })
    });
    Some(Installed(()))
}

fn installed() -> bool {
    CONTEXT.with(|c| c.borrow().is_some())
}

/// Names the fit the next screens belong to.
pub fn label(split: &str, fold: Option<usize>, learner: &str) {
    CONTEXT.with(|c| {
        if let Some(ctx) = c.borrow_mut().as_mut() {
            ctx.label = (split.to_string(), fold, learner.to_string());
        }
    });
}

/// Records an outer fold's training samples and candidates.
pub fn note_fold(split: &str, fold: usize, train: usize, candidates: usize) {
    CONTEXT.with(|c| {
        if let Some(ctx) = c.borrow_mut().as_mut() {
            ctx.folds.push(FoldNote {
                split: split.to_string(),
                fold,
                train,
                candidates,
            });
        }
    });
}

/// The process's CPU seconds (user and system, every thread) and peak
/// resident size in KB.
fn usage() -> (f64, u64) {
    // SAFETY: getrusage writes the struct it is given and reads nothing else.
    let mut u: libc::rusage = unsafe { std::mem::zeroed() };
    let ok = unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut u) } == 0;
    if !ok {
        return (f64::NAN, 0);
    }
    let secs = |t: libc::timeval| t.tv_sec as f64 + t.tv_usec as f64 / 1e6;
    (secs(u.ru_utime) + secs(u.ru_stime), u.ru_maxrss as u64)
}

/// A started fit's clock, when a context is installed.
pub struct FitClock {
    wall: Instant,
    cpu: f64,
}

pub fn start_fit() -> Option<FitClock> {
    installed().then(|| FitClock {
        wall: Instant::now(),
        cpu: usage().0,
    })
}

/// Records an outer fit's time from its clock.
pub fn end_fit(clock: Option<FitClock>, split: &str, fold: usize, learner: &str) {
    let Some(clock) = clock else { return };
    let seconds = clock.wall.elapsed().as_secs_f64();
    let (cpu, peak) = usage();
    CONTEXT.with(|c| {
        if let Some(ctx) = c.borrow_mut().as_mut() {
            ctx.times.push(FitTime {
                split: split.to_string(),
                fold,
                learner: learner.to_string(),
                seconds,
                cpu_seconds: cpu - clock.cpu,
                peak_kb: peak,
            });
        }
    });
}

/// The screen of one fit of the forest, boosting or the lasso: with no
/// context, `learn::screen`, unchanged; with one, the arm's rule on the
/// fitting samples, its counts recorded under the current label.
pub fn screen_fit(x: &BitMatrix, y: &Bits, fit: &Bits, k: usize) -> Vec<u32> {
    let rule = CONTEXT.with(|c| c.borrow().as_ref().map(|ctx| ctx.rule));
    let Some(rule) = rule else {
        return super::screen(x, y, fit, k);
    };
    let (feats, note) = match rule {
        Rule::ChiSquare => {
            let feats = super::screen(x, y, fit, k);
            let kept = feats.len();
            (feats, (None, None, None, kept))
        }
        Rule::Cmh => {
            let s = CONTEXT.with(|c| {
                let ctx = c.borrow();
                let ctx = ctx.as_ref().expect("installed");
                assert_eq!(
                    ctx.lineage.len(),
                    x.cols,
                    "the strata cover every sample of the matrix"
                );
                cmh_screen(x, y, fit, &ctx.lineage, k)
            });
            let kept = s.feats.len();
            (
                s.feats,
                (
                    Some(s.kept_at_zero),
                    Some(s.strata),
                    Some(s.informative),
                    kept,
                ),
            )
        }
    };
    CONTEXT.with(|c| {
        if let Some(ctx) = c.borrow_mut().as_mut() {
            let (split, fold, learner) = ctx.label.clone();
            ctx.screens.push(ScreenNote {
                split,
                fold,
                learner,
                rule: rule.name(),
                fitting: fit.count() as usize,
                candidates: x.rows,
                kept: note.3,
                kept_at_zero: note.0,
                strata: note.1,
                informative_strata: note.2,
            });
        }
    });
    feats
}

/// One dataset's inputs file: counts and fit times only.
#[derive(Debug, Serialize)]
pub struct InputsResult {
    pub format: &'static str,
    pub evaluation_revision: u32,
    pub dataset: String,
    pub tag: String,
    pub seed: u64,
    pub arm: &'static str,
    pub rule: &'static str,
    pub budget: usize,
    pub min_support: u32,
    pub samples: usize,
    /// The dataset matrix's rows: the candidates in.
    pub dataset_rows: usize,
    pub folds: Vec<FoldNote>,
    pub screens: Vec<ScreenNote>,
    pub fit_times: Vec<FitTime>,
}

pub fn file_name(dataset: &str, tag: &str, seed: u64) -> String {
    format!("{dataset}.{tag}-s{seed}.inputs.json")
}

/// Writes this thread's context as `<out>/<dataset>.<tag>-s<seed>.inputs.json`.
#[allow(clippy::too_many_arguments)]
pub fn finish(
    dataset: &str,
    tag: &str,
    seed: u64,
    budget: usize,
    min_support: u32,
    samples: usize,
    dataset_rows: usize,
) -> Result<()> {
    let Some(settings) = active() else {
        return Ok(());
    };
    let ctx = CONTEXT
        .with(|c| c.borrow_mut().take())
        .context("no C4 context is installed")?;
    let result = InputsResult {
        format: FORMAT,
        evaluation_revision: crate::analysis::EVALUATION_REVISION,
        dataset: dataset.to_string(),
        tag: tag.to_string(),
        seed,
        arm: settings.arm(),
        rule: ctx.rule.name(),
        budget,
        min_support,
        samples,
        dataset_rows,
        folds: ctx.folds,
        screens: ctx.screens,
        fit_times: ctx.times,
    };
    write_json(&settings.out.join(file_name(dataset, tag, seed)), &result)?;
    eprintln!(
        "    {dataset}: C4 arm {}: rule {}, budget {budget}, floor {min_support}, {} screens recorded",
        settings.arm(),
        result.rule,
        result.screens.len()
    );
    Ok(())
}

fn write_json(path: &Path, value: &impl Serialize) -> Result<()> {
    let bytes = serde_json::to_vec_pretty(value)?;
    let part = path.with_extension("json.part");
    std::fs::write(&part, &bytes).with_context(|| format!("write {}", part.display()))?;
    std::fs::rename(&part, path).with_context(|| format!("rename to {}", path.display()))
}

#[cfg(test)]
mod tests {
    //! SYNTHETIC tests: made-up matrices, labels and strata, no real genome.
    use super::*;
    use crate::bits;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    fn matrix(rows: &[Vec<bool>]) -> BitMatrix {
        let n = rows[0].len();
        let mut x = BitMatrix::with_cols(n);
        for r in rows {
            x.push_row(Bits::from_indices(n, (0..n).filter(|&i| r[i])).words());
        }
        x
    }

    fn labels(v: &[bool]) -> Bits {
        Bits::from_indices(v.len(), (0..v.len()).filter(|&i| v[i]))
    }

    fn random_cohort(seed: u64, n: usize, rows: usize, strata: u32) -> (BitMatrix, Bits, Vec<u32>) {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let lineage: Vec<u32> = (0..n).map(|_| rng.random_range(0..strata) * 7 + 3).collect();
        let y: Vec<bool> = (0..n).map(|_| rng.random_bool(0.45)).collect();
        let x: Vec<Vec<bool>> = (0..rows)
            .map(|_| {
                let p = rng.random_range(0.05..0.95);
                (0..n).map(|_| rng.random_bool(p)).collect()
            })
            .collect();
        (matrix(&x), labels(&y), lineage)
    }

    /// The chi-square `learn::screen` ranks by, for the one-stratum check.
    fn pearson(x: &BitMatrix, y: &Bits, fit: &Bits, f: usize) -> f64 {
        let pos = fit.and(y);
        let neg = fit.and_not(y);
        let (np, nn) = (pos.count() as f64, neg.count() as f64);
        let a = bits::and_count(x.row(f), pos.words()) as f64;
        let b = bits::and_count(x.row(f), neg.words()) as f64;
        let (c, d) = (np - a, nn - b);
        let denom = (a + b) * (c + d) * np * nn;
        if denom == 0.0 { 0.0 } else { (np + nn) * (a * d - b * c).powi(2) / denom }
    }

    #[test]
    fn with_one_stratum_the_statistic_is_the_chi_square_times_n_less_one_over_n() {
        let (x, y, _) = random_cohort(7, 90, 60, 1);
        let fit = Bits::from_indices(90, (0..90).filter(|i| i % 5 != 0));
        let one = vec![4u32; 90];
        let (stats, strata, informative) = cmh_statistics(&x, &y, &fit, &one);
        assert_eq!((strata, informative), (1, 1));
        let n = fit.count() as f64;
        for (f, s) in stats.iter().enumerate() {
            let want = pearson(&x, &y, &fit, f) * (n - 1.0) / n;
            assert!((s - want).abs() <= 1e-9 * want.max(1.0), "row {f}: {s} against {want}");
        }
        // So the two screens keep the same rows, in the same order.
        assert_eq!(cmh_screen(&x, &y, &fit, &one, 12).feats, super::super::screen(&x, &y, &fit, 12));
    }

    #[test]
    fn a_stratum_of_one_class_or_one_sample_changes_no_statistic() {
        let (x, y, lineage) = random_cohort(11, 120, 40, 6);
        let fit = Bits::ones(120);
        let base = cmh_statistics(&x, &y, &fit, &lineage).0;
        // Move every resistant sample of one stratum to a stratum of its own
        // (one class), and one susceptible sample to another (one sample).
        let mut moved = lineage.clone();
        for i in 0..120 {
            if lineage[i] == 3 && y.contains(i) {
                moved[i] = 1_000;
            }
        }
        let lone = (0..120).find(|&i| !y.contains(i) && lineage[i] == 10).unwrap();
        moved[lone] = 2_000;
        let (after, strata, informative) = cmh_statistics(&x, &y, &fit, &moved);
        assert!(strata >= informative + 2);
        // The strata left lose those samples, so compare with them left out of the fit.
        let mut keep = fit.clone();
        for i in 0..120 {
            if moved[i] == 1_000 || moved[i] == 2_000 {
                keep.remove(i);
            }
        }
        let without = cmh_statistics(&x, &y, &keep, &moved).0;
        assert_eq!(after, without);
        assert_ne!(after, base);
    }

    #[test]
    fn a_pure_lineage_marker_scores_zero_and_the_chi_square_tops_it() {
        // Stratum 0 is mostly resistant and stratum 1 mostly susceptible; row 0
        // is present in stratum 0 alone, row 1 tracks resistance inside both.
        let n = 40;
        let lineage: Vec<u32> = (0..n).map(|i| (i >= 20) as u32).collect();
        let y: Vec<bool> = (0..n).map(|i| if i < 20 { i % 5 != 0 } else { i % 5 == 0 }).collect();
        let marker: Vec<bool> = (0..n).map(|i| i < 20).collect();
        let inside: Vec<bool> = (0..n).map(|i| y[i] != (i % 7 == 1)).collect();
        let x = matrix(&[marker, inside]);
        let y = labels(&y);
        let fit = Bits::ones(n);
        let (stats, _, informative) = cmh_statistics(&x, &y, &fit, &lineage);
        assert_eq!(informative, 2);
        assert_eq!(stats[0], 0.0);
        assert!(stats[1] > 0.0);
        assert_eq!(super::super::screen(&x, &y, &fit, 1), vec![0]);
        assert_eq!(cmh_screen(&x, &y, &fit, &lineage, 1).feats, vec![1]);
    }

    #[test]
    fn no_sample_outside_the_fit_changes_the_screen() {
        let (x, y, lineage) = random_cohort(19, 150, 80, 9);
        let fit = Bits::from_indices(150, (0..150).filter(|i| i % 4 != 1));
        let want = cmh_screen(&x, &y, &fit, &lineage, 20);
        // Flip every label and every presence bit of the samples outside it,
        // and move them to other strata.
        let mut y2 = y.clone();
        let mut lineage2 = lineage.clone();
        let mut x2 = BitMatrix::with_cols(150);
        for f in 0..x.rows {
            // Each outside bit flipped, each fitting bit kept.
            let row = Bits::from_indices(150, (0..150).filter(|&i| x.get(f, i) == fit.contains(i)));
            x2.push_row(row.words());
        }
        for i in (0..150).filter(|&i| !fit.contains(i)) {
            if y2.contains(i) { y2.remove(i) } else { y2.insert(i) }
            lineage2[i] = 99;
        }
        assert_eq!(cmh_screen(&x2, &y2, &fit, &lineage2, 20), want);
    }

    #[test]
    fn the_strata_do_not_change_with_the_labels() {
        let (_, y, lineage) = random_cohort(23, 100, 1, 5);
        let fit = Bits::from_indices(100, (0..100).filter(|i| i % 3 != 0));
        let (slot, _, all) = fit_strata(&y, &fit, &lineage);
        let flipped = Bits::from_indices(100, (0..100).filter(|&i| !y.contains(i)));
        let (slot2, _, all2) = fit_strata(&flipped, &fit, &lineage);
        // Flipping every label keeps each stratum's two classes, so the
        // strata and their order are the same.
        assert_eq!((slot, all), (slot2, all2));
    }

    #[test]
    fn ties_break_by_the_lower_row_and_the_budget_fills_with_zeros() {
        let n = 30;
        let lineage: Vec<u32> = (0..n).map(|i| (i % 3) as u32).collect();
        let y: Vec<bool> = (0..n).map(|i| i % 2 == 0).collect();
        let a: Vec<bool> = (0..n).map(|i| i % 2 == 0 && i % 4 != 2).collect();
        let constant_in_strata: Vec<bool> = (0..n).map(|i| i % 3 == 0).collect();
        let x = matrix(&[constant_in_strata.clone(), a.clone(), constant_in_strata, a]);
        let s = cmh_screen(&x, &labels(&y), &Bits::ones(n), &lineage, 4);
        assert_eq!(s.feats, vec![1, 3, 0, 2]);
        assert_eq!(s.kept_at_zero, 2);
    }

    #[test]
    fn the_statistic_matches_the_plain_python_implementation_on_synthetic_tables() {
        // tests/c4_screen/cmh_cases.json: written by scripts/c4-screen/cmh_plain.py
        // in exact rationals; tests/c4_screen/test_cmh_plain.py holds the file to it.
        let cases: serde_json::Value =
            serde_json::from_str(include_str!("../../tests/c4_screen/cmh_cases.json")).unwrap();
        let cases = cases["cases"].as_array().unwrap();
        assert!(cases.len() >= 5);
        for case in cases {
            let parse = |key: &str| -> Vec<bool> {
                case[key].as_str().unwrap().bytes().map(|b| b == b'1').collect()
            };
            let y = labels(&parse("labels"));
            let fit = labels(&parse("fit"));
            let lineage: Vec<u32> = case["strata"].as_array().unwrap().iter().map(|v| v.as_u64().unwrap() as u32).collect();
            let rows: Vec<Vec<bool>> = case["rows"]
                .as_array()
                .unwrap()
                .iter()
                .map(|r| r.as_str().unwrap().bytes().map(|b| b == b'1').collect())
                .collect();
            let x = matrix(&rows);
            let (stats, strata, informative) = cmh_statistics(&x, &y, &fit, &lineage);
            assert_eq!(strata as u64, case["strata_in_fit"].as_u64().unwrap());
            assert_eq!(informative as u64, case["informative"].as_u64().unwrap());
            for (f, want) in case["statistics"].as_array().unwrap().iter().enumerate() {
                let want = want.as_f64().unwrap();
                assert!((stats[f] - want).abs() <= 1e-9 * want.max(1.0), "case {} row {f}: {} against {want}", case["name"], stats[f]);
            }
            let k = case["k"].as_u64().unwrap() as usize;
            let order: Vec<u32> = case["top"].as_array().unwrap().iter().map(|v| v.as_u64().unwrap() as u32).collect();
            let got = cmh_screen(&x, &y, &fit, &lineage, k);
            assert_eq!(got.feats, order, "case {}", case["name"]);
            assert_eq!(got.kept_at_zero as u64, case["top_at_zero"].as_u64().unwrap());
        }
    }

    /// The negative control (docs/screen-arms.md), on a SYNTHETIC cohort:
    /// labels permuted within strata, 200 times; a permutation rejects when
    /// the top 100 hold more lineage markers than a one-sided hypergeometric
    /// test at 0.05 allows. The adjusted screen passes at 15 or fewer of 200.
    /// The chi-square's count is printed beside it and decides nothing.
    #[test]
    fn negative_control_permuted_labels_within_strata_carry_no_lineage_signal() {
        let (adjusted, chi) = negative_control(200);
        eprintln!("negative control: adjusted screen rejects in {adjusted} of 200 permutations (cut 15); chi-square in {chi} of 200");
        assert!(adjusted <= 15, "the adjusted screen rejects in {adjusted} of 200");
    }

    /// (adjusted rejections, chi-square rejections) over `perms` permutations.
    pub(crate) fn negative_control(perms: u64) -> (usize, usize) {
        let (strata, per, rows, top) = (20usize, 20usize, 2000usize, 100usize);
        let n = strata * per;
        let mut rng = ChaCha8Rng::seed_from_u64(2026_11_23);
        let lineage: Vec<u32> = (0..n).map(|i| (i / per) as u32).collect();
        // Stratum k holds round(per * share) resistant samples, shares 0.1 to 0.9.
        let resistant: Vec<usize> = (0..strata)
            .map(|k| ((0.1 + 0.8 * k as f64 / (strata - 1) as f64) * per as f64).round() as usize)
            .collect();
        // Half the rows are lineage markers: present at 0.8 in a random half of
        // the strata and 0.2 in the rest; the others at 0.5 everywhere. Their
        // order is shuffled.
        let mut kinds: Vec<bool> = (0..rows).map(|r| r < rows / 2).collect();
        use rand::seq::SliceRandom;
        kinds.shuffle(&mut rng);
        let mut x = BitMatrix::with_cols(n);
        for &marker in &kinds {
            let high: Vec<bool> = if marker {
                let mut h: Vec<bool> = (0..strata).map(|k| k < strata / 2).collect();
                h.shuffle(&mut rng);
                h
            } else {
                vec![false; strata]
            };
            let row = Bits::from_indices(
                n,
                (0..n).filter(|&i| {
                    let p = if !marker { 0.5 } else if high[i / per] { 0.8 } else { 0.2 };
                    rng.random_bool(p)
                }),
            );
            x.push_row(row.words());
        }
        let markers = kinds.iter().filter(|&&k| k).count();
        let critical = hypergeometric_critical(rows, markers, top, 0.05);
        let fit = Bits::ones(n);
        let (mut adjusted, mut chi) = (0, 0);
        for p in 0..perms {
            let mut prng = ChaCha8Rng::seed_from_u64(1_000 + p);
            let mut y = Bits::zeros(n);
            for k in 0..strata {
                let mut members: Vec<usize> = (k * per..(k + 1) * per).collect();
                members.shuffle(&mut prng);
                for &i in &members[..resistant[k]] {
                    y.insert(i);
                }
            }
            let count = |feats: &[u32]| feats.iter().filter(|&&f| kinds[f as usize]).count();
            adjusted += (count(&cmh_screen(&x, &y, &fit, &lineage, top).feats) >= critical) as usize;
            chi += (count(&super::super::screen(&x, &y, &fit, top)) >= critical) as usize;
        }
        (adjusted, chi)
    }

    /// The smallest c with P(X >= c) <= alpha, X hypergeometric: `draws`
    /// from `total` holding `good`.
    fn hypergeometric_critical(total: usize, good: usize, draws: usize, alpha: f64) -> usize {
        let ln_choose = |n: usize, k: usize| -> f64 {
            (1..=k).map(|i| ((n - k + i) as f64).ln() - (i as f64).ln()).sum()
        };
        let pmf = |x: usize| -> f64 {
            if x > good || draws - x > total - good {
                return 0.0;
            }
            (ln_choose(good, x) + ln_choose(total - good, draws - x) - ln_choose(total, draws)).exp()
        };
        let mut tail = 0.0;
        for c in (0..=draws).rev() {
            tail += pmf(c);
            if tail > alpha {
                return c + 1;
            }
        }
        0
    }
}
