"""scripts/mechanism/combiner.py, the script of docs/combiner-result.md.

Every report, genome id, label, score, species and antibiotic here is SYNTHETIC and made up in the test. Each count
below is worked by hand from the fixture's design: a report is a few groups of genomes that share their five
scores, copied five times with one copy per stored outer fold, so every fold's fit sees the same shares.
"""

import contextlib
import fractions
import io
import json
import math
import os
import random
import tempfile
import unittest
from unittest import mock

from tests.mechanism import common as c
from tests.mechanism.test_learner_overlap import Report, lay_out

lo = c.load("learner_overlap")
os2 = c.load("overlap_second")
common = c.load("common")
cb = c.load("combiner")
FIVE = ("scm", "cart", "rf", "gbt", "lasso")
SPECIES = ["Synthia alba", "Synthia beta", "Synthia gamma", "Synthia delta", "Synthia epsilon", "Synthia zeta",
           "Synthia eta"]
LOW = (0.2,) * 5
HIGH = (0.8,) * 5


def one(j, high=0.8, low=0.2):
    """Learner j at `high`, the other four at `low`."""
    return tuple(high if k == j else low for k in range(5))


# rf (learner 2) alone calls one group of 9 resistant and 1 susceptible genomes; one other learner calls each of
# four groups of 1 resistant and 2 susceptible; none calls 1 resistant and 4 susceptible. Per copy 14 resistant and
# 13 susceptible, 135 genomes in all. The mean misses every resistant genome (sums 1.0 and 1.6). The weighted mean's
# optimum is all weight on rf: at that vertex the slope of rf's group, -0.6*9/0.8 + 0.6*1/0.2 = -3.75, is under each
# other group's slope at 0, -0.6*1/0.2 + 0.6*2/0.8 = -1.5. The logistic combiner fits each group's own share (six
# groups, six parameters): 0.2, 0.9 and 1/3. So both arms call rf's group alone: 45 resistant and 5 susceptible.
LED = [(LOW, 1, 4), (one(2), 9, 1)] + [(one(j), 1, 2) for j in (0, 1, 3, 4)]

# Symmetric in the five and in score against 1 - score with the labels swapped: every form calls a genome resistant
# where the five sum to more than 2.5. Per copy 30 resistant and 30 susceptible; called: 15 + 9 resistant, 5 + 1
# susceptible.
BLOCK = ([(LOW, 1, 9)] + [(one(j), 1, 3) for j in range(5)] + [(one(j, 0.2, 0.8), 3, 1) for j in range(5)]
         + [(HIGH, 9, 1)])

# Six groups, six parameters: the logistic combiner fits each group's share, 0.1, 0.75, 1/3, 2/3, 7/15 and 0.6.
# 17 resistant of 40 per copy, a prevalence of 0.425, so the group at 7/15 lies between the prevalence and 0.5.
SATURATED = [(LOW, 1, 9), (one(0), 3, 1), (one(1), 1, 2), (one(2), 2, 1), (one(3), 7, 8), (one(4), 3, 2)]

# scm scores every resistant genome exactly 1 and every susceptible one exactly 0, and the other four learners call
# one susceptible group each: the logistic combiner has no finite answer. Per copy 3 resistant and 7 susceptible.
# The mean misses every resistant genome (sum 1.8) and calls no susceptible one (sums 0.8 and 1.4).
SCM0 = (0.0, 0.2, 0.2, 0.2, 0.2)
SEPARABLE = ([(SCM0, 0, 3), ((1.0, 0.2, 0.2, 0.2, 0.2), 3, 0)]
             + [(tuple(0.0 if k == 0 else (0.8 if k == j else 0.2) for k in range(5)), 0, 1) for j in (1, 2, 3, 4)])

# scm scores every susceptible genome exactly 1 and every resistant one exactly 0; the other four agree on each
# genome, at 0.2 and 0.8. At the equal weights the score is 0.36 and 0.64 and the loss -log 0.64; at scm's vertex
# every score is clipped, the loss is -log 0.0001 and its slope in every weight is 0.
ANTI = [((1.0, 0.2, 0.2, 0.2, 0.2), 0, 2), ((0.0, 0.8, 0.8, 0.8, 0.8), 2, 0)]

# rf against the other four, which agree: only rf's weight moves the loss.
LEAN_A, LEAN_B = (0.1, 0.1, 0.9, 0.1, 0.1), (0.9, 0.9, 0.1, 0.9, 0.9)
LEAN = [(LEAN_A, 5, 1), (LEAN_B, 1, 3)]

# gbt is 1 - rf on every genome, in values a float holds exactly, and no column has one value or is another's copy:
# the intercept, rf and gbt are tied, so the first Newton step has no finite solution. Per copy 13 resistant of 24.
Q, T3 = 0.25, 0.75
AFFINE = [((Q, Q, Q, T3, Q), 1, 4), ((Q, Q, T3, Q, Q), 9, 1), ((T3, Q, Q, T3, Q), 1, 2), ((Q, T3, Q, T3, Q), 1, 2),
          ((Q, Q, Q, T3, T3), 1, 2)]

# 2 resistant of 15 per copy: 10 resistant genomes in five copies, the least a report may hold.
TEN = [(LOW, 1, 4), (one(2), 1, 1)] + [(one(j), 0, 2) for j in (0, 1, 3, 4)]

# Eight groups, made up from seeded draws, whose fit from zero takes a fifth Newton step that raises the loss at 1,
# 1/2 and 1/4 of its length: the search halves it to 1/32, and the fit converges in 11 steps at a loss of 0.125.
STEEP = [((0.01, 0.5, 0.01, 0.5, 0.05), False, 1), ((0.1, 0.25, 0.95, 0.1, 0.75), True, 50),
         ((0.1, 0.5, 0.99, 0.5, 0.99), False, 10), ((0.25, 0.25, 0.75, 0.99, 0.9), False, 10),
         ((0.25, 0.25, 0.9, 0.9, 0.95), True, 1), ((0.75, 0.9, 0.9, 0.1, 0.5), True, 10),
         ((0.95, 0.5, 0.5, 0.99, 0.25), True, 3), ((0.99, 0.1, 0.95, 0.75, 0.5), False, 50)]

# What fit_logistic returned at 7897a215, before R30, on fitting genomes with neither a one-score column nor two
# columns the same: float.hex of the intercept and five coefficients, the steps and the end. Folds 0 to 4 of
# made_up(400, 7); SEPARABLE's ten genomes; LED's 27 genomes with a cap of 2 Newton steps.
BEFORE_R30 = {
    0: (["-0x1.fc72552b76d5ep+4", "0x1.ae93e30627613p+3", "0x1.cb0e987efffcep+2", "0x1.5054a4ff892e3p+4",
         "0x1.5e3f7fd60e8a0p+4", "0x1.a6e4b4929fad8p+2"], 12, "converged"),
    1: (["-0x1.3cb759b041aaap+5", "0x1.1d5b4024f9119p+4", "0x1.0d44a5c4779c8p+3", "0x1.4e27ce6b5e71cp+4",
         "0x1.e7c6c143c23edp+4", "0x1.299b49dc60a06p+3"], 13, "converged"),
    2: (["-0x1.1788a4ef7704ep+6", "0x1.a05b70031b0a6p+4", "0x1.39780d9590260p+4", "0x1.4fda92fab290fp+5",
         "0x1.9cc6c548204f3p+5", "0x1.1a38d2c02bab8p+4"], 14, "converged"),
    3: (["-0x1.38b97d931e628p+5", "0x1.19b38f83d3d91p+4", "0x1.80fbc4892bc2cp+2", "0x1.d261a93ed1c1ap+4",
         "0x1.cffb8d5e552b6p+4", "0x1.4a138eba13d1cp+2"], 13, "converged"),
    4: (["-0x1.3d69eeb00f6f7p+5", "0x1.387fc0a7b0983p+4", "0x1.550b391618af2p+3", "0x1.c4565c70a205ep+4",
         "0x1.6b264365b675fp+4", "0x1.4f715ce695672p+2"], 13, "converged"),
    "separable": (["-0x1.b33f0e9691511p+4", "0x1.b33f0e9691514p+5", "-0x1.9658efef6617ap-48",
                   "-0x1.e39e9037c4ccep-48", "-0x1.23cf9d58102c0p-47", "-0x1.ba76b5eb86585p-47"], 26,
                  "no_finite_step"),
    "cap": (["-0x1.b8bc2720f717fp+1", "0x1.23cbf67f99b05p+0", "0x1.23cbf67f99b07p+0", "0x1.716c26a96aeeep+2",
             "0x1.23cbf67f99b04p+0", "0x1.23cbf67f99b01p+0"], 2, "cap"),
}


def made_up(n, seed, strength=(220, 160, 280, 250, 120)):
    """A SYNTHETIC report from seeded integer draws: about 35 in 100 genomes resistant, five scores in thousandths
    that lean with the label and share a draw, folds by genome order. No two columns are the same and none has one
    value."""
    rng = random.Random(seed)
    S, y, folds = [], [], []
    for i in range(n):
        t = rng.random() < 0.35
        shared = rng.randrange(-150, 151)
        S.append(tuple(min(999, max(1, 460 + (a if t else -a) + shared + rng.randrange(-300, 301))) / 1000
                       for a in strength))
        y.append(t)
        folds.append(i % 5)
    return S, y, folds


def reviewers_report(n, seed, prev=0.35, strength=(1.6, 1.2, 2.0, 1.8, 1.0), shared=0.8):
    """The review's SYNTHETIC maker of finding 1 (synth.continuous, review of 25 Nov), as it was written: five
    learners, each a noisy sigmoid of the label with a shared noise term."""
    def sig(z):
        return 1.0 / (1.0 + math.exp(-z)) if z >= 0 else math.exp(z) / (1.0 + math.exp(z))

    rng = random.Random(seed)
    S, y, folds = [], [], []
    for i in range(n):
        t = rng.random() < prev
        common_ = rng.gauss(0, shared)
        S.append(tuple(sig(a * (1 if t else -1) + common_ + rng.gauss(0, 1.0) - 0.3) for a in strength))
        y.append(t)
        folds.append(i % 5)
    order = list(range(n))
    rng.shuffle(order)
    return [S[i] for i in order], [y[i] for i in order], [folds[k] for k in range(n)]


def with_column(S, j, value):
    """Learner j's score replaced on every genome by value(s), the genome's five scores."""
    return [tuple(value(s) if k == j else x for k, x in enumerate(s)) for s in S]


def fitting(S, y, folds, f):
    return cb.groups_of(S, y, cb.training_genomes(folds, f))


def gradient(beta, groups):
    """The mean gradient of the logistic loss in all six numbers, at beta."""
    n = sum(k for _, _, k in groups)
    g = [0.0] * 6
    for s, t, k in groups:
        r = k * (cb.logistic_score(beta, s) - (1.0 if t else 0.0))
        for i, x in enumerate((1.0,) + s):
            g[i] += r * x
    return [v / n for v in g]


def designed(species, antibiotic, groups, tag, copies=5):
    labels, scores, folds = [], [], []
    for copy in range(copies):
        for s, resistant, susceptible in groups:
            for t in [True] * resistant + [False] * susceptible:
                labels.append(t)
                scores.append(s)
                folds.append(copy)
    r = Report(species, antibiotic, sum(labels), len(labels) - sum(labels), tag=tag)
    r.labels, r.catalog, r.folds = labels, list(labels), folds
    r.scores = {n: [s[k] for s in scores] for k, n in enumerate(FIVE)}
    r.random = dict(r.scores)
    return r


def from_arrays(species, antibiotic, S, y, folds, tag):
    """A SYNTHETIC report from a score tuple, label and fold per genome."""
    r = Report(species, antibiotic, sum(y), len(y) - sum(y), tag=tag)
    r.labels, r.catalog, r.folds = list(y), list(y), list(folds)
    r.scores = {n: [s[k] for s in S] for k, n in enumerate(FIVE)}
    r.random = dict(r.scores)
    return r


def as_arrays(groups, copies=5):
    S, y, folds = [], [], []
    for copy in range(copies):
        for s, resistant, susceptible in groups:
            for t in [True] * resistant + [False] * susceptible:
                S.append(s)
                y.append(t)
                folds.append(copy)
    return S, y, folds


def kept_set():
    return [designed(SPECIES[i], "drug a", LED, tag=i + 1) for i in range(6)]


def mixed_set():
    return [designed(SPECIES[i], "drug b", BLOCK, tag=i + 1) for i in range(6)] + [
        designed(SPECIES[6], "drug a", LED, tag=7)]


def run(module, *argv, **kw):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = module.main([str(x) for x in argv], **kw)
    return rc, buf.getvalue()


def gate_line(out, n):
    lines = [line for line in out.split("\n") if line.startswith(f"gate {n} ")]
    return lines[0] if lines else ""


def dataset_like(species, mean, arm):
    """A row of datasets.tsv as the record writes it: `mean` and `arm` are (tp, fn, fp, tn)."""
    row = {"species": species}
    for form, counts in (("mean", mean), ("weighted", arm), ("logistic", arm)):
        for k, v in zip(("tp", "fn", "fp", "tn"), counts):
            row[f"{form}_{k}"] = v
        row[f"{form}_called"] = sum(counts)
        row[f"{form}_log_loss"] = "0.5"
    return row


class Case(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.top = self.tmp.name
        self.first = os.path.join(self.top, "first")
        self.second = os.path.join(self.top, "second")
        self.reading = os.path.join(self.top, "reading")
        self.record = os.path.join(self.top, "record")

    def tearDown(self):
        self.tmp.cleanup()

    def records(self, reports):
        """The learner overlap's record and its second pass's, built from the same SYNTHETIC reports."""
        d = lay_out(self.top, reports)
        n = len(reports)
        for args in (("read", d, "--out", os.path.join(self.top, "r1"), "--expect", n),
                     ("write", os.path.join(self.top, "r1"), "--out", self.first)):
            rc, out = run(lo, *args)
            self.assertEqual(rc, 0, out)
        for args in (("read", d, "--first", self.first, "--out", os.path.join(self.top, "r2"), "--expect", n),
                     ("write", os.path.join(self.top, "r2"), "--first", self.first, "--out", self.second)):
            rc, out = run(os2, *args)
            self.assertEqual(rc, 0, out)
        return d

    def read(self, reports=None, before=None, expect=None, count=None):
        """The read on SYNTHETIC reports; gate 1 is handed the test's own count, the number of reports unless
        `count` says otherwise, since the command line holds a read to the freeze's 268."""
        reports = kept_set() if reports is None else reports
        d = self.records(reports)
        if before:
            before(d)
        return d, run(cb, "read", d, "--first", self.first, "--second", self.second, "--out", self.reading,
                      "--expect", len(reports) if expect is None else expect,
                      reports=len(reports) if count is None else count)

    def read_write(self, reports=None):
        d, (rc, out) = self.read(reports)
        self.assertEqual(rc, 0, out)
        rc, out2 = run(cb, "write", self.reading, "--reports", d, "--out", self.record)
        self.assertEqual(rc, 0, out2)
        return d, out + out2

    def rows(self):
        return {(r["species"], r["antibiotic"]): r for r in lo.read_tsv(os.path.join(self.record, "datasets.tsv"))}

    def tamper(self, change):
        path = os.path.join(self.reading, cb.READING)
        with open(path, encoding="utf-8") as f:
            reading = json.load(f)
        change(reading)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(reading, f)

    def edit_reports(self, d, change):
        """Change the lineage split of every SYNTHETIC report file after the records are built."""
        for name in os.listdir(d):
            path = os.path.join(d, name)
            with open(path, encoding="utf-8") as f:
                doc = json.load(f)
            for split in doc["splits"]:
                if split["scheme"] == "lineage":
                    change(doc, split)
            with open(path, "w", encoding="utf-8") as f:
                json.dump(doc, f)

    def read_safely(self, reports, before=None):
        """The read, with any crash turned into a failure: a stop must be the read's own named line."""
        try:
            return self.read(reports, before=before)[1]
        except Exception as e:      # noqa: BLE001 - a crash is the fault this test is meant to catch
            self.fail(f"the read crashed with {type(e).__name__} in place of its own stop line")


class Forms(unittest.TestCase):
    def test_a_score_of_exactly_the_cutoff_is_a_call(self):
        self.assertEqual(cb.CUTOFF, 0.5)
        self.assertTrue(cb.is_call(0.5))
        self.assertFalse(cb.is_call(0.4999999))
        self.assertTrue(cb.mean_call((0.5,) * 5))

    def test_the_incumbent_is_the_mean_not_the_majority(self):
        self.assertTrue(cb.mean_call((0.9, 0.9, 0.3, 0.3, 0.3)))       # sum 2.7, two of five call
        self.assertFalse(cb.mean_call((0.6, 0.6, 0.6, 0.1, 0.1)))      # sum 2.0, three of five call

    def test_the_weighted_mean_stays_inside_the_five_scores(self):
        # weights a hair under one in sum: the plain weighted sum of five scores of 0.5 is 0.4999999999999995
        w = [0.2 * (1 - 1e-15)] * 5
        self.assertTrue(cb.is_call(cb.weighted_score(w, (0.5,) * 5)))
        self.assertFalse(cb.is_call(cb.weighted_score([0.2] * 5, (0.49, 0.1, 0.2, 0.3, 0.4))))

    def test_log_loss_is_clipped(self):
        # a score of 1.0 on a susceptible genome costs -log(1e-4), not an infinite loss
        self.assertEqual((cb.clipped(1.0), cb.clipped(0.0)), (1.0 - 1e-4, 1e-4))
        got = cb.log_loss([True, False], [1.0, 1.0])
        self.assertAlmostEqual(got, (-math.log(1 - 1e-4) - math.log(1e-4)) / 2, places=12)
        self.assertAlmostEqual(cb.log_loss([True], [0.0]), -math.log(1e-4), places=12)

    def test_measures_are_worked_from_the_counts(self):
        counts = {"tp": 120, "fn": 30, "fp": 20, "tn": 130, "called": 300, "log_loss": 0.25}
        self.assertAlmostEqual(cb.measure(counts, "recall"), 0.8, places=15)                 # 120 / 150
        self.assertAlmostEqual(cb.measure(counts, "specificity"), 130 / 150, places=15)
        self.assertAlmostEqual(cb.measure(counts, "vme"), 0.2, places=15)                    # 30 / 150, not 20 / 150
        self.assertAlmostEqual(cb.measure(counts, "balanced_accuracy"), (0.8 + 130 / 150) / 2, places=15)
        self.assertEqual(cb.measure(counts, "log_loss"), 0.25)

    def test_balanced_accuracy_is_not_coverage(self):
        # every genome called, so coverage is 1; the balanced accuracy is 0.8
        counts = {"tp": 8, "fn": 2, "fp": 2, "tn": 8, "called": 20, "log_loss": 0.1}
        self.assertEqual(cb.coverage(counts, 20), 1)
        self.assertAlmostEqual(cb.measure(counts, "balanced_accuracy"), 0.8, places=15)

    def test_coverage_is_the_share_given_a_call(self):
        # four genomes, one given no call: coverage 3/4, while 2 of 4 are right
        counts = cb.confusion([True, True, False, False], [True, False, False, None])
        self.assertEqual(counts, {"tp": 1, "fn": 1, "fp": 0, "tn": 1, "called": 3})
        self.assertEqual(cb.coverage(counts, 4), fractions.Fraction(3, 4))

    def test_plan_gives_the_frozen_forms(self):
        plan = cb.plan()
        self.assertEqual(plan["forms"], ["mean", "weighted", "logistic"])
        self.assertEqual((plan["cutoff"], plan["clip"], plan["margin"]), (0.5, 1e-4, 0.02))
        self.assertEqual((plan["resamples"], plan["interval_values"], plan["min_species"]), (4000, [101, 3900], 6))
        self.assertEqual(plan["incumbent"], "math.fsum(five scores) >= 2.5")
        # R19's stop rule: a full backtracking search, down to a step of 1e-20, from the equal weights, capped at 5,000
        # steps; Newton's method with no penalty, capped at 100 steps
        w, lg = plan["weighted_fit"], plan["logistic_fit"]
        self.assertEqual((w["cap"], w["start"], w["step_floor"]), (5000, "equal weights", 1e-20))
        self.assertEqual((lg["cap"], lg["penalty"]), (100, None))
        # R30: the logistic fit's stop values, the end that stops the read, the 268 and the flag's fold size
        self.assertEqual((lg["gradient"], lg["step"], lg["pivot"], lg["start"], lg["stops_the_read"]),
                         (1e-10, 1e-8, 1e-13, "zero", "no_finite_first_step"))
        self.assertEqual((plan["reports"], plan["flag_fold"]), (268, 2))


class Fits(unittest.TestCase):
    def test_the_weighted_mean_fits_plain_log_loss(self):
        # the loss is 8 * -log(p) + 2 * -log(1 - p) with p = 0.1 + 0.8 * w_rf on group A, so p = 0.8 and w_rf = 0.875;
        # with each class weighted to one half it would be p = 19/24 and w_rf = 0.8646
        w, steps, end = cb.fit_weighted(cb.groups_of(*as_arrays(LEAN, copies=1)[:2], range(10)))
        self.assertEqual(end, "converged")
        self.assertAlmostEqual(w[2], 0.875, places=6)
        for j in (0, 1, 3, 4):
            self.assertAlmostEqual(w[j], 0.125 / 4, places=6)
        self.assertAlmostEqual(math.fsum(w), 1.0, places=12)

    def test_the_weights_stay_on_the_simplex(self):
        got = cb.project_simplex([0.5, 0.5, 0.5, 0.0, 0.0])
        for a, b in zip(got, [1 / 3, 1 / 3, 1 / 3, 0.0, 0.0]):
            self.assertAlmostEqual(a, b, places=15)
        self.assertEqual(cb.project_simplex([2.0, 0.0, 0.0, 0.0, 0.0]), [1.0, 0.0, 0.0, 0.0, 0.0])
        w, _, end = cb.fit_weighted(cb.groups_of(*as_arrays(LED, copies=1)[:2], range(27)))
        self.assertEqual(end, "converged")
        self.assertEqual(w, [0.0, 0.0, 1.0, 0.0, 0.0])

    def test_the_logistic_fits_each_groups_share(self):
        S, y, folds = as_arrays(SATURATED)
        probs, recs = cb.crossfit("logistic", S, y, folds)
        self.assertEqual([r["end"] for r in recs], ["converged"] * 5)
        shares = {s: p for s, p in zip(S, probs)}
        for s, share in zip([g[0] for g in SATURATED], (0.1, 0.75, 1 / 3, 2 / 3, 7 / 15, 0.6)):
            self.assertAlmostEqual(shares[s], share, places=6)
        # called at 0.5: the groups at 0.75, 2/3 and 0.6; the group at 7/15 is not called
        self.assertEqual(cb.confusion(y, [cb.is_call(p) for p in probs]),
                         {"tp": 40, "fn": 45, "fp": 20, "tn": 95, "called": 200})

    def test_a_separable_fold_ends_with_no_finite_step(self):
        # no finite maximum: the Newton system loses its last finite solution before 100 steps, and the last finite
        # iterate stands; it calls the 3 resistant genomes of each copy and none of the 7 susceptible
        S, y, folds = as_arrays(SEPARABLE)
        beta, steps, end = cb.fit_logistic(cb.groups_of(S, y, range(10)))
        self.assertEqual(end, "no_finite_step")
        self.assertLess(steps, cb.L_CAP)
        self.assertTrue(cb.finite(*beta))
        probs, recs = cb.crossfit("logistic", S, y, folds)
        self.assertEqual([(r["end"], r["not_finite"]) for r in recs], [("no_finite_step", 0)] * 5)
        self.assertEqual(cb.confusion(y, [cb.is_call(p) for p in probs]),
                         {"tp": 15, "fn": 0, "fp": 0, "tn": 35, "called": 50})

    def test_a_fit_at_its_cap_ends_there(self):
        # one weighted step and two Newton steps: neither fit is done, and each stands where its cap stops it
        groups = cb.groups_of(*as_arrays(LED, copies=1)[:2], range(27))
        with mock.patch.object(cb, "W_CAP", 1), mock.patch.object(cb, "L_CAP", 2):
            w, w_steps, w_end = cb.fit_weighted(groups)
            beta, l_steps, l_end = cb.fit_logistic(groups)
        self.assertEqual((w_steps, w_end, l_steps, l_end), (1, "cap", 2, "cap"))
        self.assertLess(cb.weighted_loss_grad(w, groups)[0], cb.weighted_loss_grad([0.2] * 5, groups)[0])

    def test_the_weighted_fit_is_never_above_the_equal_mean(self):
        # it starts from the equal weights and takes only a lower loss; on ANTI it leaves scm, whose vertex is the
        # worst point of the simplex at -log 0.0001
        for design in (LED, BLOCK, SATURATED, SEPARABLE, LEAN, ANTI):
            S, y, _ = as_arrays(design, copies=1)
            groups = cb.groups_of(S, y, range(len(y)))
            w, _, end = cb.fit_weighted(groups)
            self.assertEqual(end, "converged")
            self.assertLessEqual(cb.weighted_loss_grad(w, groups)[0], cb.weighted_loss_grad([0.2] * 5, groups)[0])
        S, y, _ = as_arrays(ANTI, copies=1)
        groups = cb.groups_of(S, y, range(4))
        self.assertAlmostEqual(cb.weighted_loss_grad([0.2] * 5, groups)[0], -math.log(0.64), places=12)
        self.assertAlmostEqual(cb.weighted_loss_grad([1.0, 0, 0, 0, 0], groups)[0], -math.log(cb.CLIP), places=12)
        self.assertLess(cb.weighted_loss_grad(cb.fit_weighted(groups)[0], groups)[0], -math.log(0.64))
        # only a strictly lower loss is taken: a point of the same loss ends the search
        self.assertFalse(cb.lower(0.3, 0.3))
        self.assertTrue(cb.lower(0.29, 0.3))

    def test_a_fit_on_neither_case_is_the_fit_before_r30(self):
        # no column with one value and no two columns the same: bit for bit what 7897a215 returned, in all three ends
        S, y, folds = made_up(400, 7)
        cases = {f: fitting(S, y, folds, f) for f in range(5)}
        cases["separable"] = cb.groups_of(*as_arrays(SEPARABLE, copies=1)[:2], range(10))
        for name, groups in cases.items():
            self.assertEqual(cb.column_sets(groups), ([], [[j] for j in range(5)]), name)
            beta, steps, end = cb.fit_logistic(groups)
            self.assertEqual(([b.hex() for b in beta], steps, end), BEFORE_R30[name], name)
        with mock.patch.object(cb, "L_CAP", 2):
            beta, steps, end = cb.fit_logistic(cb.groups_of(*as_arrays(LED, copies=1)[:2], range(27)))
        self.assertEqual(([b.hex() for b in beta], steps, end), BEFORE_R30["cap"])

    def test_a_one_score_column_is_set_aside(self):
        # lasso at 0.3 on every genome (the review's finding 1): its coefficient is exactly 0.0, and the fit is a
        # maximum-likelihood fit of all six numbers. With rf's copy as gbt too, both rules apply in one fold
        S, y, folds = made_up(400, 7)
        one_value = with_column(S, 4, lambda s: 0.3)
        for S2, sets in ((one_value, [[0], [1], [2], [3]]), (with_column(one_value, 3, lambda s: s[2]),
                                                             [[0], [1], [2, 3]])):
            for f in range(5):
                groups = fitting(S2, y, folds, f)
                self.assertEqual(cb.column_sets(groups), ([4], sets))
                beta, steps, end = cb.fit_logistic(groups)
                self.assertEqual(end, "converged")
                self.assertEqual(beta[5], 0.0)
                self.assertLess(max(abs(v) for v in gradient(beta, groups)), 1e-9)
                self.assertTrue(all(b != 0.0 for b in beta[:5]))

    def test_columns_the_same_are_fitted_as_one(self):
        # gbt a copy of rf on every genome: one column fitted, its coefficient in two even parts, and a
        # maximum-likelihood fit of all six numbers
        S, y, folds = made_up(400, 7)
        S2 = with_column(S, 3, lambda s: s[2])
        for f in range(5):
            groups = fitting(S2, y, folds, f)
            self.assertEqual(cb.column_sets(groups), ([], [[0], [1], [2, 3], [4]]))
            beta, steps, end = cb.fit_logistic(groups)
            self.assertEqual(end, "converged")
            self.assertEqual(beta[3], beta[4])
            self.assertNotEqual(beta[3], 0.0)
            self.assertLess(max(abs(v) for v in gradient(beta, groups)), 1e-9)

    def test_a_held_out_genome_is_scored_from_both_columns(self):
        # gbt equals rf on every genome but fold 0's: fold 0's fit sees them the same, and a held-out genome whose two
        # scores differ is scored from both, each with half the coefficient the one column was fitted with. That
        # coefficient is the fit's on the four other columns: the fit with gbt at one value gives it
        S, y, folds = made_up(400, 7)
        S2 = [tuple((x if g == 0 else s[2]) if j == 3 else x for j, x in enumerate(s)) for s, g in zip(S, folds)]
        groups = fitting(S2, y, folds, 0)
        beta, _, end = cb.fit_logistic(groups)
        whole, _, _ = cb.fit_logistic([(tuple(0.5 if j == 3 else x for j, x in enumerate(s)), t, k)
                                       for s, t, k in groups])
        self.assertEqual(end, "converged")
        self.assertEqual((beta[3], beta[4]), (whole[3] / 2, whole[3] / 2))
        self.assertEqual(whole[4], 0.0)
        probs, recs = cb.crossfit("logistic", S2, y, folds)
        self.assertEqual((recs[0]["columns_aside"], recs[0]["columns_shared"]), (0, 2))
        held = [i for i, g in enumerate(folds) if g == 0 and S2[i][2] != S2[i][3]]
        self.assertGreater(len(held), 70)
        for i in held:
            s = S2[i]
            self.assertEqual(probs[i], cb.logistic_score(beta, s))
            self.assertAlmostEqual(probs[i], cb.logistic_score(whole[:3] + [whole[3], 0.0, whole[5]],
                                                               s[:2] + ((s[2] + s[3]) / 2, 0.0, s[4])), places=12)
            self.assertNotEqual(probs[i], cb.logistic_score(whole, s))      # not from rf alone

    def test_swapping_two_columns_the_same_moves_no_held_out_score(self):
        # rf and gbt the same on fold 0's fitting genomes and apart on fold 0's: swapping the two learners leaves
        # every held-out score of fold 0 the same to the bit. Ordering all five another way solves the Newton system
        # in another order, and its rounding moves no score by 1e-11 or more
        S, y, folds = made_up(400, 7)
        S2 = [tuple((x if g == 0 else s[2]) if j == 3 else x for j, x in enumerate(s)) for s, g in zip(S, folds)]
        before, _ = cb.crossfit("logistic", S2, y, folds)
        swapped, _ = cb.crossfit("logistic", [(s[0], s[1], s[3], s[2], s[4]) for s in S2], y, folds)
        self.assertEqual([p for p, g in zip(swapped, folds) if g == 0], [p for p, g in zip(before, folds) if g == 0])
        order = (4, 3, 0, 2, 1)
        moved, recs = cb.crossfit("logistic", [tuple(s[j] for j in order) for s in S2], y, folds)
        self.assertEqual((recs[0]["columns_shared"], recs[0]["end"]), (2, "converged"))
        for a, b, g in zip(moved, before, folds):
            if g == 0:
                self.assertLess(abs(a - b), 1e-11)

    def test_the_same_and_one_score_are_read_without_tolerance(self):
        # a copy of rf one float step apart on one fitting genome is not the same column, and lasso at 0.3 but one
        # float step apart on one fitting genome does not have one value: neither is set aside or shared, and the
        # fit on five columns, close to tied, has no finite first step
        S, y, folds = made_up(400, 7)
        i0 = folds.index(1)
        for j, S2 in ((3, with_column(S, 3, lambda s: s[2])), (4, with_column(S, 4, lambda s: 0.3))):
            S2[i0] = tuple(math.nextafter(x, 1.0) if k == j else x for k, x in enumerate(S2[i0]))
            groups = fitting(S2, y, folds, 0)
            self.assertEqual(cb.column_sets(groups), ([], [[k] for k in range(5)]))
            self.assertEqual(cb.fit_logistic(groups)[1:], (0, cb.NO_FIRST_STEP))
            _, recs = cb.crossfit("logistic", S2, y, folds)
            self.assertEqual((recs[0]["columns_aside"], recs[0]["columns_shared"], recs[0]["end"]),
                             (0, 0, cb.NO_FIRST_STEP))

    def test_a_first_step_with_no_finite_solution_is_not_scored(self):
        # gbt is 1 - rf: no column is set aside or shared, the first Newton step has no finite solution, and no
        # held-out genome is scored, so none scores 0.5 from coefficients left at 0
        S, y, folds = as_arrays(AFFINE)
        beta, steps, end = cb.fit_logistic(cb.groups_of(S, y, range(24)))
        self.assertEqual((beta, steps, end), ([0.0] * 6, 0, cb.NO_FIRST_STEP))
        probs, recs = cb.crossfit("logistic", S, y, folds)
        self.assertEqual([(r["end"], r["columns_aside"], r["columns_shared"]) for r in recs],
                         [(cb.NO_FIRST_STEP, 0, 0)] * 5)
        self.assertEqual(probs, [None] * 120)
        # a later step with no finite solution still lets the last finite iterate stand (R19 item 4)
        S, y, folds = as_arrays(SEPARABLE)
        probs, recs = cb.crossfit("logistic", S, y, folds)
        self.assertEqual({r["end"] for r in recs}, {"no_finite_step"})
        self.assertNotIn(None, probs)

    def test_the_step_search_halves_under_a_quarter(self):
        # STEEP's fifth Newton step raises the loss at 1, 1/2 and 1/4 of its length: the search goes on halving to a
        # lower loss, and the fit converges at a loss of 0.125, from log 2 at zero. A search that stops at 1/2 takes
        # that step at 1/4 whatever the loss, and its fit ends with no finite step at a loss over 30,000, calling 20
        # susceptible genomes resistant that this fit does not
        groups = sorted(STEEP)
        beta, steps, end = cb.fit_logistic(groups)
        self.assertEqual((steps, end), (11, "converged"))
        self.assertAlmostEqual(cb.logistic_nll(beta, groups), 0.125073, places=6)
        self.assertEqual([cb.is_call(cb.logistic_score(beta, s)) for s, _, _ in groups],
                         [False, True, False, False, False, True, True, False])

    def test_the_newton_system_is_pivoted(self):
        # scm's scores are tiny but not one value, 1e-8 to 4e-8: the system's pivot on scm's own row is under 1e-13
        # while its column holds larger entries. With pivoting the fit converges in 9 steps; without, its first step
        # has no finite solution
        S, y = [], []
        for i in range(60):
            t = i % 3 == 0
            S.append(((1 + i % 4) * 1e-8,) + tuple(min(95, max(5, (65 if t else 35) + (i * m) % 61 - 30)) / 100
                                                     for m in (7, 11, 13, 17)))
            y.append(t)
        groups = cb.groups_of(S, y, range(60))
        self.assertEqual(cb.column_sets(groups), ([], [[k] for k in range(5)]))
        self.assertEqual(cb.fit_logistic(groups)[1:], (9, "converged"))


class Crossfit(Case):
    def test_each_fold_is_fitted_on_the_other_folds(self):
        # the training count is the genomes inside the groups the fit is handed (12 groups of 108 genomes), beside
        # the fold's training list
        S, y, folds = as_arrays(LED)
        for arm in cb.ARMS:
            _, recs = cb.crossfit(arm, S, y, folds)
            self.assertEqual([(r["fold"], r["listed"], r["trained"], r["tested"], r["shared"]) for r in recs],
                             [(f, 108, 108, 27, 0) for f in range(5)])
        self.assertEqual(len(fitting(S, y, folds, 0)), 12)
        # what the fit receives, not only the fold record: flip every label of one fold, and that fold's held-out
        # scores do not move by one bit, while the other folds' do
        for arm in cb.ARMS:
            before, _ = cb.crossfit(arm, S, y, folds)
            for f in range(5):
                flipped = [(not t) if g == f else t for t, g in zip(y, folds)]
                after, _ = cb.crossfit(arm, S, flipped, folds)
                self.assertEqual([p for p, g in zip(after, folds) if g == f],
                                 [p for p, g in zip(before, folds) if g == f], (arm, f))
                self.assertNotEqual([p for p, g in zip(after, folds) if g != f],
                                    [p for p, g in zip(before, folds) if g != f], (arm, f))

    def test_gate_2_counts_no_shared_genome_on_the_standard_read(self):
        _, (rc, out) = self.read()
        self.assertEqual(rc, 0, out)
        self.assertEqual(gate_line(out, 2), "gate 2 cross-fitting: 60 fold fits, 0 whose fit was handed another "
                                            "count of genomes than its training list, or whose training and held-out "
                                            "genomes overlap or do not add up to the report: pass")


class Counts(Case):
    def test_the_rf_led_report(self):
        # each genome's catalogue call is the opposite of its label: no count below moves, since no catalogue call
        # enters a fit or a count
        reports = kept_set()
        for rep in reports:
            rep.catalog = [not t for t in rep.labels]
        self.read_write(reports)
        r = self.rows()[(SPECIES[0], "drug a")]
        self.assertEqual({k: r[k] for k in ("genomes", "resistant", "susceptible", "mean_miss_by_1to4",
                                            "mean_miss_by_all5")},
                         {"genomes": "135", "resistant": "70", "susceptible": "65", "mean_miss_by_1to4": "65",
                          "mean_miss_by_all5": "5"})
        self.assertEqual([r[f"mean_{k}"] for k in cb.COUNTS], ["0", "70", "0", "65", "135"])
        for arm in cb.ARMS:
            self.assertEqual([r[f"{arm}_{k}"] for k in cb.COUNTS], ["45", "25", "5", "60", "135"])
            self.assertEqual([r[f"{arm}_{t}"] for t in cb.TRANSITIONS], ["45", "0", "0", "0", "5"])
            self.assertEqual([r[f"{arm}_folds_{e}"] for e in cb.ENDS[arm]], ["5"] + ["0"] * (len(cb.ENDS[arm]) - 1))
            self.assertEqual(r[f"{arm}_marked"], "no")
        # weighted: 6 genomes at -log 0.2 and 21 at -log 0.8 per copy of 27
        self.assertAlmostEqual(float(r["weighted_log_loss"]), (6 * math.log(5) + 21 * math.log(1.25)) / 27, places=10)
        # transitions.tsv over all six reports: 420 resistant, all missed by the incumbent, 390 by 1 to 4 of the five
        # and 30 by all five; 390 susceptible, none called by the incumbent; each arm recovers 270 and newly calls 30
        cols = ("datasets", "resistant", "mean_miss", "mean_miss_by_1to4", "mean_miss_by_all5", "susceptible",
                "mean_wrong") + cb.TRANSITIONS
        got = {t["arm"]: [t[k] for k in cols] for t in lo.read_tsv(os.path.join(self.record, "transitions.tsv"))
               if t["scope"] == "all"}
        want = ["6", "420", "420", "390", "30", "390", "0", "270", "0", "0", "0", "30"]
        self.assertEqual(got, {arm: want for arm in cb.ARMS})

    def test_the_symmetric_report(self):
        self.read_write(mixed_set())
        r = self.rows()[(SPECIES[0], "drug b")]
        self.assertEqual((r["genomes"], r["mean_miss_by_1to4"], r["mean_miss_by_all5"]), ("300", "25", "5"))
        for form in cb.FORMS:
            self.assertEqual([r[f"{form}_{k}"] for k in cb.COUNTS], ["120", "30", "30", "120", "300"])
        for arm in cb.ARMS:
            self.assertEqual([r[f"{arm}_{t}"] for t in cb.TRANSITIONS], ["0"] * 5)

    def test_the_transitions_split_by_how_many_of_the_five_call(self):
        labels = [True, True, True, True, False, False]
        mean = [False, False, False, True, True, False]
        arm = [True, True, False, False, False, True]
        called_by = [0, 1, 4, 3, 3, 1]
        self.assertEqual(cb.transitions(labels, mean, arm, called_by),
                         {"recovered_1to4": 1, "recovered_all5": 1, "lost": 1, "cleared": 1, "new_wrong": 1})

    def test_print_gives_the_rule_table(self):
        self.read_write()
        rc, out = run(cb, "print", self.record)
        self.assertEqual(rc, 0, out)
        for arm in cb.ARMS:
            self.assertIn(f"| {arm} | all | 6 | 6 | 0.2830 | 0.2830 to 0.2830 | -0.6429 | -0.6429 to -0.6429 | 6 "
                          "| kept | kept |", out)
            self.assertIn(f"| {arm} | without marked | 6 | 6 | 0.2830 | 0.2830 to 0.2830 | -0.6429 | -0.6429 to "
                          "-0.6429 | 6 | kept | decides nothing |", out)
            self.assertIn(f"| {arm} | without one score fold | 6 | 6 | 0.2830 | 0.2830 to 0.2830 | -0.6429 | -0.6429 "
                          "to -0.6429 | 6 | kept | kept |", out)

    def test_a_report_is_flagged_by_one_fold(self):
        # flagged when some learner gives every genome of some fold of 2 or more genomes one score, read with ==
        S, y, folds = made_up(40, 3)
        self.assertFalse(cb.one_score_fold(S, folds))
        in_fold_0 = [tuple(0.3 if (j == 4 and g == 0) else x for j, x in enumerate(s)) for s, g in zip(S, folds)]
        self.assertTrue(cb.one_score_fold(in_fold_0, folds))
        in_fold_0[0] = tuple(math.nextafter(x, 1.0) if j == 4 else x for j, x in enumerate(in_fold_0[0]))
        self.assertFalse(cb.one_score_fold(in_fold_0, folds))
        # a fold of one genome: one score by itself, not a flag
        self.assertFalse(cb.one_score_fold(S, [9 if i == 0 else g for i, g in enumerate(folds)]))
        self.assertTrue(cb.one_score_fold(with_column(S, 1, lambda s: 0.6), folds))


class Rule(Case):
    def test_the_gain_must_reach_the_margin(self):
        self.assertTrue(cb.keep(0.02, 0.001, -0.01))
        self.assertFalse(cb.keep(0.0199, 0.001, -0.01))

    def test_the_gain_interval_must_lie_above_0(self):
        self.assertFalse(cb.keep(0.03, 0.0, -0.01))
        self.assertFalse(cb.keep(0.03, -0.01, -0.02))
        self.assertFalse(cb.keep(0.03, None, -0.01))

    def test_a_higher_vme_blocks_the_arm(self):
        self.assertFalse(cb.keep(0.03, 0.01, 0.001))
        self.assertTrue(cb.keep(0.03, 0.01, 0.0))
        # the rule reads the lower end of the interval: in 4 species the arm gains 0.3 and its very major error rate
        # falls by 0.1, in 3 it gains 0.2 and the rate rises by 0.1. The interval of the change runs from -0.1 to 0.1,
        # so it does not lie wholly above 0, and the arm is kept
        rows = ([dataset_like(SPECIES[i], (5, 5, 5, 5), (6, 4, 0, 10)) for i in range(4)]
                + [dataset_like(SPECIES[i], (5, 5, 5, 5), (4, 6, 0, 10)) for i in range(4, 7)])
        got = cb.rule(rows, "weighted")
        self.assertEqual((got["median_gain"], got["vme_lo"], got["vme_hi"], got["reading"]), (0.3, -0.1, 0.1, "kept"))

    def test_kept_and_no_difference_shown(self):
        # six reports gaining 103/364 in balanced accuracy each, the very major error rate down by 9/14 each
        self.read_write()
        rules = {r["arm"]: r for r in lo.read_tsv(os.path.join(self.record, "rule.tsv")) if r["count"] == "all"}
        for arm in cb.ARMS:
            self.assertEqual({k: rules[arm][k] for k in ("datasets", "species", "median_gain", "gain_lo", "gain_hi",
                                                         "vme_median", "vme_lo", "gain_at_least_margin", "reading")},
                             {"datasets": "6", "species": "6", "median_gain": "0.282967032967",
                              "gain_lo": "0.282967032967", "gain_hi": "0.282967032967", "vme_median": "-0.642857142857",
                              "vme_lo": "-0.642857142857", "gain_at_least_margin": "6", "reading": "kept"})

    def test_a_median_gain_of_0_is_no_difference_shown(self):
        # six symmetric reports at a gain of 0 and one at 103/364: the median gain is 0
        self.read_write(mixed_set())
        rules = {r["arm"]: r for r in lo.read_tsv(os.path.join(self.record, "rule.tsv")) if r["count"] == "all"}
        for arm in cb.ARMS:
            self.assertEqual((rules[arm]["median_gain"], rules[arm]["gain_at_least_margin"], rules[arm]["reading"]),
                             ("0.000000000000", "1", "no difference shown"))
        self.assertEqual(len(lo.read_tsv(os.path.join(self.record, "pooled.tsv"))), (1 + 7) * 2 * 5)

    def test_the_second_count_leaves_out_the_marked_reports(self):
        # the six rf-led reports and one SEPARABLE report, whose 5 logistic folds end with no finite step. Over all 7
        # the median gain is 103/364 (six at 103/364, one at 1/2, the arms calling SEPARABLE without error) and it
        # decides; the logistic arm's second count leaves the marked report out, 6 reports, and decides nothing
        self.read_write(kept_set() + [designed(SPECIES[6], "drug c", SEPARABLE, tag=7)])
        rules = {(r["arm"], r["count"]): r for r in lo.read_tsv(os.path.join(self.record, "rule.tsv"))}
        self.assertEqual(list(rules), [(a, k) for a in cb.ARMS for k in cb.COUNTS_OF_RULE])
        got = {key: (r["decides"], r["datasets"], r["species"], r["median_gain"], r["outcome"], r["reading"])
               for key, r in rules.items()}
        self.assertEqual(got, {
            ("weighted", "all"): ("yes", "7", "7", "0.282967032967", "kept", "kept"),
            ("weighted", "without_marked"): ("no", "7", "7", "0.282967032967", "kept", "decides nothing"),
            ("weighted", "without_one_score_fold"): ("yes", "7", "7", "0.282967032967", "kept", "kept"),
            ("logistic", "all"): ("yes", "7", "7", "0.282967032967", "kept", "kept"),
            ("logistic", "without_marked"): ("no", "6", "6", "0.282967032967", "kept", "decides nothing"),
            ("logistic", "without_one_score_fold"): ("yes", "7", "7", "0.282967032967", "kept", "kept")})
        r = self.rows()[(SPECIES[6], "drug c")]
        self.assertEqual([r[k] for k in ("weighted_marked", "logistic_marked", "logistic_folds_no_finite_step")],
                         ["no", "yes", "5"])
        for arm in cb.ARMS:
            self.assertEqual([r[f"{arm}_{k}"] for k in cb.COUNTS], ["15", "0", "0", "35", "50"])
        with open(os.path.join(self.record, "gates.txt"), encoding="utf-8") as f:
            self.assertEqual([x for x in f.read().split("\n") if x.startswith("reading ")],
                             ["reading weighted: kept", "reading logistic: kept"])

    def test_the_third_count_decides(self):
        # the six rf-led reports, one with lasso at 0.2 on every genome of fold 0: its calls do not move, and it is
        # flagged. Over all 6 the rule keeps each arm; without the flagged report 5 species are left, too few for an
        # interval, so the third count's outcome is no difference shown, and so is each arm's reading
        reports = kept_set()
        flat = reports[0]
        flat.scores["lasso"] = [0.2 if g == 0 else x for x, g in zip(flat.scores["lasso"], flat.folds)]
        flat.random = dict(flat.scores)
        self.read_write(reports)
        rows = self.rows()
        self.assertEqual([rows[(sp, "drug a")]["one_score_fold"] for sp in SPECIES[:6]], ["yes"] + ["no"] * 5)
        self.assertEqual([rows[(SPECIES[0], "drug a")][f"{a}_{k}"] for a in cb.ARMS for k in cb.COUNTS],
                         ["45", "25", "5", "60", "135"] * 2)
        rules = {(r["arm"], r["count"]): r for r in lo.read_tsv(os.path.join(self.record, "rule.tsv"))}
        for arm in cb.ARMS:
            got = [(rules[arm, k]["datasets"], rules[arm, k]["species"], rules[arm, k]["outcome"],
                    rules[arm, k]["reading"]) for k in cb.COUNTS_OF_RULE]
            self.assertEqual(got, [("6", "6", "kept", "no difference shown"), ("6", "6", "kept", "decides nothing"),
                                   ("5", "5", "no difference shown", "no difference shown")])
        with open(os.path.join(self.record, "gates.txt"), encoding="utf-8") as f:
            self.assertEqual([x for x in f.read().split("\n") if x.startswith("reading ")],
                             ["reading weighted: no difference shown", "reading logistic: no difference shown"])

    def test_the_record_says_the_vme_rule_in_r19s_words(self):
        self.read_write()
        with open(os.path.join(self.record, "gates.txt"), encoding="utf-8") as f:
            self.assertIn("rule: the very major error rate is not higher on that interval: the species-clustered "
                          "interval of the median change in the very major error rate does not lie wholly above 0\n",
                          f.read())


class Pooled(unittest.TestCase):
    def species_rows(self):
        # species A: 20 reports gaining 1; B to F: one report each gaining 0; G: one report gaining 0.5
        gain1 = ((0, 2, 2, 0), (2, 0, 0, 2))
        gain0 = ((1, 1, 1, 1), (1, 1, 1, 1))
        half = ((0, 2, 2, 0), (2, 0, 2, 0))
        return ([dataset_like("A", *gain1) for _ in range(20)] + [dataset_like(sp, *gain0) for sp in "BCDEF"]
                + [dataset_like("G", *half)])

    def test_the_interval_resamples_species(self):
        # resampling species leaves A out of a third of the draws, so the 101st median is 0; resampling the 26
        # reports almost never draws 13 that are not A's, so it would be 1
        p = cb.pooled_row(self.species_rows(), "weighted", "balanced_accuracy", "all")
        self.assertEqual((p["median"], p["lo"], p["hi"]), (1.0, 0.0, 1.0))

    def test_the_interval_is_seeded_by_the_rows_name(self):
        names = []
        orig = common.species_interval

        def record(by_species, name, resamples=common.RESAMPLES):
            names.append(name)
            return orig(by_species, name, resamples)

        with mock.patch.object(common, "species_interval", record):
            cb.pooled_rows(self.species_rows())
        self.assertEqual(names, [f"combiner|{a}|{m}" for a in cb.ARMS for m in cb.MEASURES])

    def test_better_worse_and_same_at_the_margin(self):
        self.assertEqual([cb.side(d, "balanced_accuracy") for d in (0.02, 0.0199, -0.02)], ["better", "same", "worse"])
        self.assertEqual([cb.side(d, "vme") for d in (0.02, -0.0199, -0.02)], ["worse", "same", "better"])
        self.assertEqual(cb.side(-0.02, "log_loss"), "better")

    def test_a_change_is_rounded_to_12_decimals_before_the_margin(self):
        # 4 resistant and 25 susceptible genomes; the arm clears 3 susceptible genomes, the incumbent 2. The balanced
        # accuracies are 0.06 and 0.04, a gain of exactly 0.02, which a float subtraction gives as
        # 0.019999999999999997: rounded to 12 decimals each of the 7 reports gains 0.02, at the margin
        rows = [dataset_like(SPECIES[i], (0, 4, 23, 2), (0, 4, 22, 3)) for i in range(7)]
        self.assertLess(cb.measure(cb.form_counts(rows[0], "weighted"), "balanced_accuracy")
                        - cb.measure(cb.form_counts(rows[0], "mean"), "balanced_accuracy"), 0.02)
        p = cb.pooled_row(rows, "weighted", "balanced_accuracy", "all")
        self.assertEqual((p["median"], p["better"], p["same"]), (0.02, 7, 0))
        self.assertEqual(cb.rule(rows, "weighted")["gain_at_least_margin"], 7)


class Gates(Case):
    def test_the_standard_read_passes_every_gate(self):
        self.read_write()
        with open(os.path.join(self.record, "gates.txt"), encoding="utf-8") as f:
            text = f.read()
        lines = [line for line in text.split("\n") if line.startswith("gate ")]
        self.assertEqual([line.split(" ")[1] for line in lines], [str(n) for n in range(1, 9)])
        self.assertTrue(all(line.endswith(": pass") for line in lines), text)
        self.assertEqual(sorted(os.listdir(self.record)), sorted(list(cb.RECORD_FILES) + ["SHA256SUMS"]))

    def test_gate_1_stops_on_a_report_set_that_differs(self):
        def drop_one(d):
            os.remove(os.path.join(d, sorted(os.listdir(d))[0]))

        # 5 reports read and 5 expected: only the set differs from the first record's
        _, (rc, out) = self.read(before=drop_one, expect=5, count=5)
        self.assertEqual(rc, 1, out)
        self.assertIn("5 reports read, 5 expected, 6 pooled in the first record, 1 in one and not the other: STOP",
                      gate_line(out, 1))
        self.assertFalse(os.path.exists(self.reading))

    def test_gate_1_holds_the_read_to_its_count(self):
        # the six reports are the first record's pooled rows, one for one, but the read expects 7
        _, (rc, out) = self.read(count=7)
        self.assertEqual(rc, 1, out)
        self.assertIn("6 reports read, 7 expected, 6 pooled in the first record, 0 in one and not the other: STOP",
                      gate_line(out, 1))
        self.assertFalse(os.path.exists(self.reading))

    def test_the_command_line_holds_the_read_to_268(self):
        # the command line has no option for the count: the six SYNTHETIC reports stop at gate 1, where 268 are
        # expected, and an option that tries to set it is refused
        d = self.records(kept_set())
        argv = ("read", d, "--first", self.first, "--second", self.second, "--out", self.reading, "--expect", 6)
        rc, out = run(cb, *argv)
        self.assertEqual(rc, 1, out)
        self.assertIn("6 reports read, 268 expected,", gate_line(out, 1))
        self.assertTrue(gate_line(out, 1).endswith(": STOP"), out)
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            run(cb, *argv, "--reports", 6)
        self.assertFalse(os.path.exists(self.reading))

    def test_gate_1_reads_a_report_with_exactly_10_resistant(self):
        # a report of 10 resistant genomes is pooled by the learner overlap, so the read takes it
        reports = kept_set() + [designed(SPECIES[6], "drug e", TEN, tag=7)]
        self.read_write(reports)
        self.assertEqual(self.rows()[(SPECIES[6], "drug e")]["resistant"], "10")

    def test_gate_2_stops_a_fit_on_its_own_genomes(self):
        with mock.patch.object(cb, "training_genomes", lambda folds, f: list(range(len(folds)))):
            _, (rc, out) = self.read()
        self.assertEqual(rc, 1, out)
        self.assertTrue(gate_line(out, 2).endswith(": STOP"), out)

    def test_gate_2_stops_a_fit_handed_every_genome(self):
        # the review's fault R1, in memory: the fold's training list is the other folds, but its fit is handed every
        # genome of the report, 135 in place of 108
        orig = cb.groups_of
        with mock.patch.object(cb, "groups_of", lambda S, y, keep: orig(S, y, range(len(y)))):
            _, (rc, out) = self.read()
        self.assertEqual(rc, 1, out)
        self.assertIn("gate 2 cross-fitting: 60 fold fits, 60 whose fit was handed another count of genomes than its "
                      "training list", gate_line(out, 2))
        self.assertTrue(gate_line(out, 2).endswith(": STOP"), out)

    def test_a_genome_with_no_fold_number_stops_the_read(self):
        def no_fold(d):
            self.edit_reports(d, lambda doc, split: split["test_fold"].__setitem__(3, None)
                              if doc["species"] == SPECIES[2] else None)

        rc, out = self.read_safely(kept_set(), before=no_fold)
        self.assertEqual(rc, 1, out)
        self.assertIn(f"stop folds: {SPECIES[2]} with drug a: a genome with no fold number: STOP", out)
        self.assertFalse(os.path.exists(self.reading))

    def test_a_report_with_one_fold_stops_the_read(self):
        def one_fold(d):
            self.edit_reports(d, lambda doc, split: split.update(test_fold=[0] * len(split["test_fold"]))
                              if doc["species"] == SPECIES[4] else None)

        rc, out = self.read_safely(kept_set(), before=one_fold)
        self.assertEqual(rc, 1, out)
        self.assertIn(f"stop folds: {SPECIES[4]} with drug a: one fold, so no fold has other folds' genomes to fit on: "
                      "STOP", out)
        self.assertFalse(os.path.exists(self.reading))

    def test_gate_3_stops_on_an_incumbent_that_differs_from_vote_tsv(self):
        # each of the four counts gate 3 compares, one more in vote.tsv, on a fresh set of records each time
        for column in ("mean_miss", "mean_wrong", "mean_miss_by_1to4", "mean_miss_by_all5"):
            def one_more(d, column=column):
                path = os.path.join(self.second, "vote.tsv")
                with open(path, encoding="utf-8") as f:
                    head = f.readline().rstrip("\n").split("\t")
                rows = lo.read_tsv(path)
                rows[0][column] = str(int(rows[0][column]) + 1)
                lo.write_tsv(path, head, rows)

            _, (rc, out) = self.read(before=one_more)
            self.assertEqual(rc, 1, (column, out))
            self.assertTrue(gate_line(out, 3).endswith(": STOP"), (column, out))
            self.tearDown()
            self.setUp()

    def test_gate_4_stops_on_weights_off_the_simplex(self):
        with mock.patch.object(cb, "fit_weighted", lambda groups: ([0.3] * 5, 1, "converged")):
            _, (rc, out) = self.read()
        self.assertEqual(rc, 1, out)
        self.assertTrue(gate_line(out, 4).endswith(": STOP"), out)

    def test_gate_5_counts_a_fold_with_no_finite_step(self):
        # a fold with no finite step is counted and its report marked, and the read goes on
        reports = kept_set() + [designed(SPECIES[6], "drug c", SEPARABLE, tag=7)]
        _, (rc, out) = self.read(reports)
        self.assertEqual(rc, 0, out)
        self.assertEqual(gate_line(out, 5), "gate 5 fit ends: 35 weighted fits, 35 converged, 0 at the cap; 35 logistic "
                                            "fits, 30 converged, 0 at the cap, 5 at a later step with no finite "
                                            "solution, 0 whose first Newton step has no finite solution, 0 with a "
                                            "score column set aside or shared; 0 reports marked for the weighted arm, "
                                            "1 for the logistic arm; 0 held-out scores not finite: pass")

    def test_gate_5_counts_and_marks_a_fold_at_the_cap(self):
        # one weighted step and two Newton steps: all 30 folds of each arm end at the cap, all 6 reports are marked
        # for both arms, and the second count holds no report
        with mock.patch.object(cb, "W_CAP", 1), mock.patch.object(cb, "L_CAP", 2):
            _, out = self.read_write()
        self.assertEqual(gate_line(out, 5), "gate 5 fit ends: 30 weighted fits, 0 converged, 30 at the cap; 30 logistic "
                                            "fits, 0 converged, 30 at the cap, 0 at a later step with no finite "
                                            "solution, 0 whose first Newton step has no finite solution, 0 with a "
                                            "score column set aside or shared; 6 reports marked for the weighted arm, "
                                            "6 for the logistic arm; 0 held-out scores not finite: pass")
        for r in self.rows().values():
            self.assertEqual([r[k] for k in ("weighted_folds_cap", "logistic_folds_cap", "weighted_marked",
                                             "logistic_marked")], ["5", "5", "yes", "yes"])
        rules = lo.read_tsv(os.path.join(self.record, "rule.tsv"))
        self.assertEqual([(r["count"], r["datasets"]) for r in rules],
                         [("all", "6"), ("without_marked", "0"), ("without_one_score_fold", "6")] * 2)

    def test_gate_5_stops_on_a_first_step_with_no_finite_solution(self):
        # AFFINE's five logistic folds have no finite first step: gate 5 counts them and stops, and the genomes of
        # those folds have no call, so gate 6 stops too
        reports = kept_set() + [designed(SPECIES[6], "drug d", AFFINE, tag=7)]
        _, (rc, out) = self.read(reports)
        self.assertEqual(rc, 1, out)
        self.assertEqual(gate_line(out, 5), "gate 5 fit ends: 35 weighted fits, 35 converged, 0 at the cap; 35 logistic "
                                            "fits, 30 converged, 0 at the cap, 0 at a later step with no finite "
                                            "solution, 5 whose first Newton step has no finite solution, 0 with a "
                                            "score column set aside or shared; 0 reports marked for the weighted arm, "
                                            "1 for the logistic arm; 0 held-out scores not finite: STOP")
        self.assertTrue(gate_line(out, 6).endswith("1 form and report pairs with a genome given no call: STOP"), out)
        self.assertFalse(os.path.exists(self.reading))

    def test_the_reviewers_two_reports_pass_gate_5(self):
        # the review's finding 1: gbt a copy of rf, and lasso at 0.3, on 600 genomes, 206 resistant. Every logistic
        # fold converges with a column shared or set aside; none is marked; no held-out score is 0.5
        S, y, folds = reviewers_report(600, 31)
        self.assertEqual((len(y), sum(y)), (600, 206))
        copy, one_value = with_column(S, 3, lambda s: s[2]), with_column(S, 4, lambda s: 0.3)
        for S2 in (copy, one_value):
            probs, recs = cb.crossfit("logistic", S2, y, folds)
            self.assertEqual({r["end"] for r in recs}, {"converged"})
            self.assertEqual(sum(1 for p in probs if p == 0.5), 0)
        reports = [from_arrays(SPECIES[0], "drug f", copy, y, folds, tag=1),
                   from_arrays(SPECIES[1], "drug f", one_value, y, folds, tag=2)]
        _, out = self.read_write(reports)
        self.assertEqual(gate_line(out, 5), "gate 5 fit ends: 10 weighted fits, 10 converged, 0 at the cap; 10 logistic "
                                            "fits, 10 converged, 0 at the cap, 0 at a later step with no finite "
                                            "solution, 0 whose first Newton step has no finite solution, 10 with a "
                                            "score column set aside or shared; 0 reports marked for the weighted arm, "
                                            "0 for the logistic arm; 0 held-out scores not finite: pass")
        for key in ((SPECIES[0], "drug f"), (SPECIES[1], "drug f")):
            r = self.rows()[key]
            self.assertEqual([r[k] for k in ("logistic_folds_converged", "logistic_folds_aside_or_shared",
                                             "logistic_marked")], ["5", "5", "no"])
        with open(os.path.join(self.reading, cb.READING), encoding="utf-8") as f:
            folds_read = {d["species"]: d["logistic_folds"] for d in json.load(f)["datasets"]}
        self.assertEqual({(r["columns_aside"], r["columns_shared"]) for r in folds_read[SPECIES[0]]}, {(0, 2)})
        self.assertEqual({(r["columns_aside"], r["columns_shared"]) for r in folds_read[SPECIES[1]]}, {(1, 0)})

    def test_write_stops_on_a_reading_made_under_other_rules(self):
        # a reading whose cutoff, whose counts of the rule, or whose report flags are not this script's: write stops
        # and writes nothing
        for change in (lambda reading: reading.update(cutoff=0.4),
                       lambda reading: reading.update(counts_of_rule=["all", "without_marked"]),
                       lambda reading: reading["datasets"][0].pop("one_score_fold")):
            d, (rc, out) = self.read()
            self.assertEqual(rc, 0, out)
            self.tamper(change)
            rc, out = run(cb, "write", self.reading, "--reports", d, "--out", self.record)
            self.assertEqual((rc, out), (1, "the reading was made under other rules\n"))
            self.assertFalse(os.path.exists(self.record))
            self.tearDown()
            self.setUp()

    def test_gate_5_stops_on_a_held_out_score_that_is_not_finite(self):
        with mock.patch.object(cb, "logistic_score", lambda beta, s: math.nan):
            _, (rc, out) = self.read()
        self.assertEqual(rc, 1, out)
        self.assertTrue(gate_line(out, 5).endswith("; 810 held-out scores not finite: STOP"), out)
        self.assertFalse(os.path.exists(self.reading))

    def test_gate_6_stops_on_a_genome_given_no_call(self):
        orig = cb.crossfit

        def holed(arm, S, y, folds):
            probs, recs = orig(arm, S, y, folds)
            return [None] + list(probs[1:]), recs

        with mock.patch.object(cb, "crossfit", holed):
            _, (rc, out) = self.read()
        self.assertEqual(rc, 1, out)
        self.assertTrue(gate_line(out, 6).endswith(": STOP"), out)

    def test_gate_7_stops_on_arm_counts_that_do_not_add_up(self):
        d, (rc, out) = self.read()
        self.assertEqual(rc, 0, out)
        self.tamper(lambda reading: reading["datasets"][0]["weighted_transitions"].update(lost=1))
        rc, out = run(cb, "write", self.reading, "--reports", d, "--out", self.record)
        self.assertEqual(rc, 1, out)
        self.assertTrue(out.startswith("gate 7 rows:") and ": STOP" in out, out)
        self.assertFalse(os.path.exists(self.record))

    def test_gate_7_stops_on_a_row_per_genome(self):
        orig = cb.report_counts

        def per_genome(labels, S, folds):
            row, gate = orig(labels, S, folds)
            return dict(row, calls=[int(t) for t in labels]), gate

        with mock.patch.object(cb, "report_counts", per_genome):
            _, (rc, out) = self.read()
        self.assertEqual(rc, 1, out)
        self.assertIn("gate 7 rows: 6 lists in the reading longer than 5: STOP", out)
        self.assertFalse(os.path.exists(self.reading))

    def test_gate_8_stops_on_a_genome_id_in_the_reading(self):
        # ids with no id shape, and one of them inside a species name: only the search of the id list finds it
        reports = kept_set()
        for r in reports:
            r.ids = [f"q{r.tag}z{i}" for i in range(len(r.labels))]
        reports[0].species = "Synthia q1z0"
        _, (rc, out) = self.read(reports)
        self.assertEqual(rc, 1, out)
        self.assertIn("gate 8 ids: 1 files holding a genome id, 0 id-shaped tokens or barred keys in the reading: STOP",
                      out)

    def test_gate_8_stops_on_an_id_shaped_token(self):
        d, (rc, out) = self.read()
        self.assertEqual(rc, 0, out)
        self.tamper(lambda reading: reading.update(note="12" + "3.45"))      # an id-shaped token, put together here
        rc, out = run(cb, "write", self.reading, "--reports", d, "--out", self.record)
        self.assertEqual(rc, 1, out)
        self.assertTrue(out.startswith("gate 8 ids:") and "0 holding one as a substring, 1 id-shaped" in out, out)

    def test_a_decimal_cell_stays_under_10(self):
        self.assertEqual(cb.cell(9.5), "9.500000000000")
        self.assertEqual(cb.cell(-0.25), "-0.250000000000")
        for x in (float(10), -21 / 2):
            with self.assertRaises(cb.Stop):
                cb.cell(x)


class Ids(Case):
    def test_no_id_is_held_after_parsing(self):
        d = lay_out(self.top, kept_set()[:1])
        report = cb.load_report(os.path.join(d, os.listdir(d)[0]))
        keys = []

        def walk(v):
            if isinstance(v, dict):
                keys.extend(v)
                for x in v.values():
                    walk(x)
            elif isinstance(v, list):
                for x in v:
                    walk(x)

        walk(report)
        self.assertNotIn("id", keys)
        self.assertIn("resistant", keys)

    def test_the_reading_holds_no_row_per_genome(self):
        _, (rc, out) = self.read()
        self.assertEqual(rc, 0, out)
        with open(os.path.join(self.reading, cb.READING), encoding="utf-8") as f:
            reading = json.load(f)
        self.assertEqual(cb.long_lists(reading), 0)
        self.assertEqual(len(reading["datasets"]), 6)
        for r in reading["datasets"]:
            self.assertEqual(set(r) - {"species", "antibiotic", "genomes", "resistant", "susceptible",
                                       "mean_miss_by_1to4", "mean_miss_by_all5", "mean", "weighted", "logistic",
                                       "weighted_folds", "logistic_folds", "weighted_transitions",
                                       "logistic_transitions", "weighted_marked", "logistic_marked",
                                       "one_score_fold"}, set())


if __name__ == "__main__":
    unittest.main()
