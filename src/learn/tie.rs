//! Experiment: how the tree learners break ties between equally good
//! features. Off by default; production is `TieRule::Float`.
//!
//! `Float`, the rule of every published report: the decision tree, the
//! forest and the boosted trees compare each candidate's 64-bit gain, take
//! the larger, and give bit-equal gains to the lowest-numbered feature. One
//! side of a split is summed and the other found by subtraction, so two
//! features whose gains are equal in exact arithmetic can get 64-bit gains
//! that differ in the last digits, and then the larger rounded gain wins
//! whatever its number.
//!
//! `Exact`, the experiment (`grm learn --exp-tie-rule exact`), makes the
//! stated tie rule true for cart and rf; for gbt, only between splits with
//! the same two groups.
//! - Decision tree and forest: a candidate's gain is a fraction of integer
//!   counts (`SplitImpurity`), compared exactly; equal fractions go to the
//!   lowest-numbered feature. Which candidates may split (gain above 1e-12)
//!   is still decided on the 64-bit gain, as in production. A fit checks
//!   once, before its first node, that its counts are inside the bound this
//!   arithmetic needs (`check_counts`).
//! - Boosted trees: the rule moves a split only to a lower feature with the
//!   same two groups. The node takes the production pick, then the
//!   lowest-numbered of the tree's candidates that divides the node's
//!   fitting samples into the same two groups as that pick (one's present
//!   side is the other's present or absent side). Candidates with equal
//!   gains and different groups still follow rounding: in the first round
//!   of a logistic fit the gradients are equal within a class, so gbt
//!   candidates with equal counts have equal gains; this rule does not treat
//!   them as equal (`gbt::tests::equal_counts_on_other_samples_follow_the_rounded_gain_under_both_rules`).
//!
//! With `--exp-tie-rule exact` a report adds `eval_config.tie_rule:
//! "exact"` and the input fingerprint hashes it, so `grm learn` refits it
//! and the metric audit flags it. The site reads neither field, so such a
//! report is kept out of `results/`: `grm report` stops at it, and an exact
//! run at the release's seeds needs `--exp-tie-out` to name another folder.
//!
//! A feature is a column number in the learner's own candidate matrix.
//!
//! `start_trace` (`grm learn --exp-tie-trace FILE`) writes one JSON line per
//! fit: its nodes searched, the nodes where the two rules pick different
//! features, and a digest of each tree. Tracing computes both rules' picks
//! at every searched node but splits on the configured rule's, so a traced
//! fit is the same model as an untraced one.

use crate::bits::Bits;
use crate::fnv::Fnv;
use serde::Serialize;
use std::cmp::Ordering;
use std::io::Write;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, Ordering as AtomicOrdering};

/// The tie rule of the tree learners. `Float` is production's.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum TieRule {
    #[default]
    Float,
    Exact,
}

impl TieRule {
    pub fn is_float(&self) -> bool {
        *self == TieRule::Float
    }
}

/// The class-weighted Gini impurity of a split's two sides, as a fraction
/// `num / den` of integer counts, up to a positive factor that is the same
/// for every candidate at a node. A lower impurity is a higher gain.
///
/// With class weights T/(2P) and T/(2N), T = P + N the fit's samples and P,
/// N its resistant and susceptible ones, a side holding p resistant and q
/// susceptible samples has weighted impurity T * p * q / (p * N + q * P).
/// Dividing by T, the two sides sum to
/// (a b D2 + c d D1) / (D1 D2), D1 = a N + b P, D2 = c N + d P,
/// with (a, b) the present side's counts and (c, d) the absent side's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SplitImpurity {
    num: u128,
    den: u128,
}

impl SplitImpurity {
    /// `present` and `absent` are (resistant, susceptible) counts of the two
    /// sides; `fit` is the fit's (resistant, susceptible) counts, which the
    /// fit has passed through `check_counts` once before its first node.
    pub fn new(present: (u32, u32), absent: (u32, u32), fit: (u32, u32)) -> Self {
        let (p_all, n_all) = (u64::from(fit.0), u64::from(fit.1));
        // Every product below stays in range for counts below `COUNT_BOUND`.
        let side = |(p, q): (u32, u32)| -> (u128, u128) {
            let (p, q) = (u64::from(p), u64::from(q));
            ((p * q) as u128, (p * n_all + q * p_all) as u128)
        };
        let (pq1, d1) = side(present);
        let (pq2, d2) = side(absent);
        // An empty side has no impurity: its term is 0, not 0 / 0.
        match (d1, d2) {
            (0, 0) => SplitImpurity { num: 0, den: 1 },
            (0, _) => SplitImpurity { num: pq2, den: d2 },
            (_, 0) => SplitImpurity { num: pq1, den: d1 },
            _ => SplitImpurity {
                num: pq1 * d2 + pq2 * d1,
                den: d1 * d2,
            },
        }
    }

    /// Exact order of the two fractions.
    pub fn cmp_exact(&self, other: &Self) -> Ordering {
        mul_wide(self.num, other.den).cmp(&mul_wide(other.num, self.den))
    }
}

/// Each class count of a fit must be below this for `SplitImpurity`'s
/// products to stay in range: a side's denominator is then below 2^63 and
/// each fraction's terms below 2^126.
pub const COUNT_BOUND: u32 = 1 << 31;

/// Stops a fit whose (resistant, susceptible) counts are not both below
/// `COUNT_BOUND`. A fit under the exact rule calls it once, before its first
/// node, in release builds too; no candidate is checked.
pub fn check_counts(fit: (u32, u32)) {
    assert!(
        fit.0 < COUNT_BOUND && fit.1 < COUNT_BOUND,
        "the exact tie rule needs fewer than {COUNT_BOUND} resistant and \
         {COUNT_BOUND} susceptible fitting samples; this fit has {} and {}",
        fit.0,
        fit.1
    );
}

/// The 256-bit product of two 128-bit numbers, as (high, low) halves, which
/// compare as the product does.
fn mul_wide(a: u128, b: u128) -> (u128, u128) {
    const LOW: u128 = u64::MAX as u128;
    let (a1, a0) = (a >> 64, a & LOW);
    let (b1, b0) = (b >> 64, b & LOW);
    let (p00, p01, p10, p11) = (a0 * b0, a0 * b1, a1 * b0, a1 * b1);
    let mid = (p00 >> 64) + (p01 & LOW) + (p10 & LOW);
    let low = (p00 & LOW) | (mid << 64);
    let high = p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64);
    (high, low)
}

/// A decision tree or forest candidate: 64-bit gain, feature, exact impurity.
pub(crate) type Candidate = (f64, u32, SplitImpurity);

/// Production's choice between two candidates: the larger 64-bit gain, then
/// the lower feature.
pub(crate) fn better_float<T: Copy>(a: (f64, u32, T), b: (f64, u32, T)) -> (f64, u32, T) {
    if a.0 > b.0 || (a.0 == b.0 && a.1 < b.1) {
        a
    } else {
        b
    }
}

/// `better_float` on (gain, feature) pairs, which carry no counts.
pub(crate) fn better_pair(a: (f64, u32), b: (f64, u32)) -> (f64, u32) {
    let (gain, feature, ()) = better_float((a.0, a.1, ()), (b.0, b.1, ()));
    (gain, feature)
}

/// The exact rule's choice: the lower exact impurity, then the lower feature.
pub(crate) fn better_exact(a: Candidate, b: Candidate) -> Candidate {
    match a.2.cmp_exact(&b.2) {
        Ordering::Less => a,
        Ordering::Greater => b,
        Ordering::Equal if a.1 < b.1 => a,
        Ordering::Equal => b,
    }
}

/// Whether `row` divides the samples of `mask` into the same two groups as
/// `present`, the node samples carrying another feature: its present side is
/// `present` or the rest of `mask`.
pub(crate) fn same_partition(mask: &[u64], row: &[u64], present: &[u64]) -> bool {
    let (mut same, mut swapped) = (true, true);
    for ((&m, &r), &s) in mask.iter().zip(row).zip(present) {
        let r = r & m;
        same &= r == s;
        swapped &= r == m & !s;
        if !(same || swapped) {
            return false;
        }
    }
    true
}

/// The lowest-numbered of `candidates` that divides `mask`'s samples into the
/// same two groups as `pick` (`pick` itself when none is lower).
pub(crate) fn lowest_same_partition<'a>(
    mask: &Bits,
    row_of: impl Fn(u32) -> &'a [u64],
    candidates: &[u32],
    pick: u32,
) -> u32 {
    let present = mask.and_slice(row_of(pick));
    candidates
        .iter()
        .copied()
        .filter(|&f| f < pick && same_partition(mask.words(), row_of(f), present.words()))
        .min()
        .unwrap_or(pick)
}

/// What a traced fit counts at its searched nodes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct Tally {
    /// Nodes whose split search ran.
    pub searched: u64,
    /// Decision tree and forest: nodes where two or more candidates share the
    /// best exact gain.
    pub exact_tied: u64,
    /// Decision tree and forest: nodes where the two rules pick different
    /// features of equal exact gain.
    pub differ_tie: u64,
    /// Decision tree and forest: nodes where 64-bit rounding put a feature of
    /// lower exact gain first.
    pub differ_order: u64,
    /// Boosted trees: nodes where a lower-numbered candidate divides the
    /// samples as the production pick does.
    pub differ_partition: u64,
}

impl std::ops::AddAssign for Tally {
    fn add_assign(&mut self, o: Self) {
        self.searched += o.searched;
        self.exact_tied += o.exact_tied;
        self.differ_tie += o.differ_tie;
        self.differ_order += o.differ_order;
        self.differ_partition += o.differ_partition;
    }
}

impl Tally {
    /// Counts one searched decision tree or forest node from both rules'
    /// picks and every admissible candidate.
    pub(crate) fn tree_node(
        &mut self,
        candidates: &[Candidate],
        float: Option<Candidate>,
        exact: Option<Candidate>,
    ) {
        self.searched += 1;
        let (Some(float), Some(exact)) = (float, exact) else {
            return;
        };
        let best = candidates
            .iter()
            .filter(|c| c.2.cmp_exact(&exact.2) == Ordering::Equal)
            .count();
        self.exact_tied += (best > 1) as u64;
        if float.1 != exact.1 {
            match float.2.cmp_exact(&exact.2) {
                Ordering::Equal => self.differ_tie += 1,
                _ => self.differ_order += 1,
            }
        }
    }
}

static TRACING: AtomicBool = AtomicBool::new(false);
static TRACE: Mutex<Option<std::io::BufWriter<std::fs::File>>> = Mutex::new(None);

/// Starts writing one JSON line per tree-learner fit to `path`.
pub fn start_trace(path: &std::path::Path) -> std::io::Result<()> {
    let file = std::fs::File::create(path)?;
    *TRACE.lock().unwrap() = Some(std::io::BufWriter::new(file));
    TRACING.store(true, AtomicOrdering::SeqCst);
    Ok(())
}

pub(crate) fn tracing() -> bool {
    TRACING.load(AtomicOrdering::Relaxed)
}

/// One traced fit.
#[derive(Serialize)]
pub(crate) struct FitRecord {
    /// "cart" (one tree), "rf" (a forest) or "gbt" (a boosted model).
    pub learner: &'static str,
    /// "logistic" or "censored" for boosted trees.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub objective: Option<&'static str>,
    pub rule: TieRule,
    pub seed: u64,
    pub fit_samples: usize,
    pub max_depth: u32,
    /// Boosted trees: rounds grown before early stopping.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rounds: Option<usize>,
    /// FNV-1a digest of each tree kept, in order.
    pub trees: Vec<String>,
    pub tally: Tally,
}

pub(crate) fn record(rec: &FitRecord) {
    let mut guard = TRACE.lock().unwrap();
    if let Some(w) = guard.as_mut() {
        serde_json::to_writer(&mut *w, rec).expect("trace record");
        w.write_all(b"\n").expect("trace write");
        w.flush().expect("trace flush");
    }
}

/// A digest of a tree, from the values `fields` lists for its nodes.
pub(crate) fn digest(fields: impl Iterator<Item = u64>) -> String {
    let mut h = Fnv::new();
    for v in fields {
        h.u64(v);
    }
    h.hex()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The weighted impurity of `tree.rs`, in 64-bit arithmetic.
    fn float_impurity(present: (u32, u32), absent: (u32, u32), fit: (u32, u32)) -> f64 {
        let total = f64::from(fit.0 + fit.1);
        let (wp, wn) = (
            total / (2.0 * f64::from(fit.0)),
            total / (2.0 * f64::from(fit.1)),
        );
        let side = |(p, q): (u32, u32)| {
            let (a, b) = (f64::from(p) * wp, f64::from(q) * wn);
            if a + b == 0.0 {
                0.0
            } else {
                2.0 * a * b / (a + b)
            }
        };
        side(present) + side(absent)
    }

    #[test]
    #[should_panic(expected = "the exact tie rule needs fewer than 2147483648 resistant")]
    fn the_count_bound_stops_a_fit_at_two_to_the_31() {
        check_counts((COUNT_BOUND - 1, COUNT_BOUND));
    }

    #[test]
    fn counts_just_under_the_bound_compare_exactly() {
        let b = COUNT_BOUND - 1;
        check_counts((b, b));
        let fit = (b, b);
        let pure = SplitImpurity::new((b, 0), (0, b), fit);
        let one_off = SplitImpurity::new((b, 1), (0, b - 1), fit);
        let mirrored = SplitImpurity::new((0, b), (b, 0), fit);
        assert_eq!(pure.cmp_exact(&one_off), Ordering::Less);
        assert_eq!(one_off.cmp_exact(&pure), Ordering::Greater);
        assert_eq!(pure.cmp_exact(&mirrored), Ordering::Equal);
        // Both are one half: b^2 / (2 b^2), and h / (2 b) + (b - h) / (2 b).
        let mixed = SplitImpurity::new((b, b), (0, 0), fit);
        let half = SplitImpurity::new((b / 2, b / 2), (b - b / 2, b - b / 2), fit);
        assert_eq!(mixed.cmp_exact(&half), Ordering::Equal);
    }

    #[test]
    fn exact_impurity_orders_splits_as_the_weighted_gini_does() {
        // Fit of 7 resistant and 5 susceptible samples; worked by hand:
        // present (3, 1) absent (2, 2): 3/(15+7) + 4/(10+14) = 3/22 + 1/6,
        // present (5, 0) absent (0, 3): 0, a pure split.
        let fit = (7, 5);
        let a = SplitImpurity::new((3, 1), (2, 2), fit);
        assert_eq!(
            a.cmp_exact(&SplitImpurity { num: 20, den: 66 }),
            Ordering::Equal
        );
        let pure = SplitImpurity::new((5, 0), (0, 3), fit);
        assert_eq!(
            pure.cmp_exact(&SplitImpurity { num: 0, den: 1 }),
            Ordering::Equal
        );
        assert_eq!(pure.cmp_exact(&a), Ordering::Less);
        // Over every split of a node with 4 resistant and 3 susceptible
        // samples, the exact order agrees with the 64-bit one wherever the
        // 64-bit values are more than rounding apart.
        let node = (4u32, 3u32);
        let mut splits = Vec::new();
        for a in 0..=node.0 {
            for b in 0..=node.1 {
                splits.push(((a, b), (node.0 - a, node.1 - b)));
            }
        }
        for &(p1, q1) in &splits {
            for &(p2, q2) in &splits {
                let (e1, e2) = (
                    SplitImpurity::new(p1, q1, fit),
                    SplitImpurity::new(p2, q2, fit),
                );
                let (f1, f2) = (float_impurity(p1, q1, fit), float_impurity(p2, q2, fit));
                if (f1 - f2).abs() > 1e-9 {
                    assert_eq!(
                        e1.cmp_exact(&e2),
                        f1.total_cmp(&f2),
                        "{p1:?} {q1:?} {p2:?} {q2:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn a_split_and_its_mirror_have_equal_exact_impurity() {
        // A feature and its complement on the node's samples are one split
        // with its sides swapped.
        let fit = (40, 23);
        for (present, absent) in [((9, 4), (31, 19)), ((1, 0), (39, 23)), ((17, 11), (23, 12))] {
            let a = SplitImpurity::new(present, absent, fit);
            let b = SplitImpurity::new(absent, present, fit);
            assert_eq!(a.cmp_exact(&b), Ordering::Equal);
        }
    }

    #[test]
    fn mul_wide_is_the_full_product() {
        let max = u128::MAX;
        assert_eq!(mul_wide(max, max), (max - 1, 1));
        assert_eq!(mul_wide(1 << 64, 1 << 64), (1, 0));
        assert_eq!(mul_wide(3, 5), (0, 15));
        assert_eq!(mul_wide(1 << 127, 2), (1, 0));
    }

    #[test]
    fn large_counts_compare_without_overflow() {
        // Counts near 2^30, where the cross products need 256 bits.
        let fit = (1 << 30, (1 << 30) - 7);
        let a = SplitImpurity::new(
            (1 << 29, 1 << 28),
            ((1 << 29), (1 << 30) - 7 - (1 << 28)),
            fit,
        );
        let b = SplitImpurity::new(
            ((1 << 29) + 1, 1 << 28),
            ((1 << 29) - 1, (1 << 30) - 7 - (1 << 28)),
            fit,
        );
        assert_eq!(a.cmp_exact(&a), Ordering::Equal);
        assert_ne!(a.cmp_exact(&b), Ordering::Equal);
        assert_eq!(a.cmp_exact(&b), b.cmp_exact(&a).reverse());
    }

    #[test]
    fn exact_choice_takes_the_lower_feature_of_equal_fractions() {
        let fit = (6, 6);
        let x = SplitImpurity::new((4, 1), (2, 5), fit);
        let mirror = SplitImpurity::new((2, 5), (4, 1), fit);
        // The mirror arrives with a 64-bit gain one unit higher, as the
        // subtraction can make it; the exact rule still takes feature 3.
        let a = (0.5, 3, x);
        let b = (0.5 + f64::EPSILON, 9, mirror);
        assert_eq!(better_float(a, b).1, 9);
        assert_eq!(better_exact(a, b).1, 3);
        assert_eq!(better_exact(b, a).1, 3);
        let worse = (0.9, 1, SplitImpurity::new((3, 3), (3, 3), fit));
        assert_eq!(better_exact(worse, a).1, 3);
    }

    #[test]
    fn same_partition_sees_equal_and_mirrored_sides_only() {
        let mask = Bits::from_indices(8, [0, 1, 2, 3, 5]);
        let f = Bits::from_indices(8, [0, 1, 6]);
        let present = mask.and(&f);
        let same = Bits::from_indices(8, [0, 1, 4, 7]);
        let mirror = Bits::from_indices(8, [2, 3, 5, 6]);
        let other = Bits::from_indices(8, [0, 2]);
        let w = |b: &Bits| b.words().to_vec();
        assert!(same_partition(mask.words(), &w(&same), present.words()));
        assert!(same_partition(mask.words(), &w(&mirror), present.words()));
        assert!(!same_partition(mask.words(), &w(&other), present.words()));
    }

    /// The tied count is of candidates at the best exact impurity, which is
    /// the exact pick's, not the float pick's, when rounding put a worse
    /// split first.
    #[test]
    fn tally_counts_ties_at_the_best_exact_impurity() {
        let fit = (6, 6);
        let x = SplitImpurity::new((4, 1), (2, 5), fit);
        let mirror = SplitImpurity::new((2, 5), (4, 1), fit);
        let worse = SplitImpurity::new((3, 2), (3, 4), fit);
        let worse_mirror = SplitImpurity::new((3, 4), (3, 2), fit);
        assert_eq!(x.cmp_exact(&worse), Ordering::Less);
        // The best exact impurity is shared, the float pick's is not.
        let mut t = Tally::default();
        let c = [
            (0.4 + f64::EPSILON, 2, worse),
            (0.4, 5, x),
            (0.4, 7, mirror),
        ];
        t.tree_node(&c, Some(c[0]), Some(c[1]));
        assert_eq!((t.exact_tied, t.differ_order), (1, 1));
        // The float pick's impurity is shared, the best exact one is not.
        let mut t = Tally::default();
        let d = [
            (0.4 + f64::EPSILON, 2, worse),
            (0.4, 3, worse_mirror),
            (0.39, 5, x),
        ];
        t.tree_node(&d, Some(d[0]), Some(d[2]));
        assert_eq!((t.exact_tied, t.differ_order), (0, 1));
    }

    #[test]
    fn tally_sorts_differences_into_ties_and_rounding_order() {
        let fit = (6, 6);
        let x = SplitImpurity::new((4, 1), (2, 5), fit);
        let mirror = SplitImpurity::new((2, 5), (4, 1), fit);
        let worse = SplitImpurity::new((3, 2), (3, 4), fit);
        let mut t = Tally::default();
        let c = [(0.5, 3, x), (0.5 + f64::EPSILON, 9, mirror)];
        t.tree_node(&c, Some(c[1]), Some(c[0]));
        let d = [(0.4, 2, worse), (0.4 - f64::EPSILON, 5, x)];
        t.tree_node(&d, Some(d[0]), Some(d[1]));
        t.tree_node(&[], None, None);
        assert_eq!(
            t,
            Tally {
                searched: 3,
                exact_tied: 1,
                differ_tie: 1,
                differ_order: 1,
                differ_partition: 0
            }
        );
    }
}
