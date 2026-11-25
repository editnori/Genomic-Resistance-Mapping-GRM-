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

# scm alone calls every resistant genome and no susceptible one: the logistic combiner has no finite answer.
SEPARABLE = [(LOW, 0, 3), (one(0), 3, 0)] + [(one(j), 0, 1) for j in (1, 2, 3, 4)]

# rf against the other four, which agree: only rf's weight moves the loss.
LEAN_A, LEAN_B = (0.1, 0.1, 0.9, 0.1, 0.1), (0.9, 0.9, 0.1, 0.9, 0.9)
LEAN = [(LEAN_A, 5, 1), (LEAN_B, 1, 3)]


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


def run(module, *argv):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = module.main([str(x) for x in argv])
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

    def read(self, reports=None, before=None, expect=None):
        reports = kept_set() if reports is None else reports
        d = self.records(reports)
        if before:
            before(d)
        return d, run(cb, "read", d, "--first", self.first, "--second", self.second, "--out", self.reading,
                      "--expect", len(reports) if expect is None else expect)

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


class Fits(unittest.TestCase):
    def test_the_weighted_mean_fits_plain_log_loss(self):
        # the loss is 8 * -log(p) + 2 * -log(1 - p) with p = 0.1 + 0.8 * w_rf on group A, so p = 0.8 and w_rf = 0.875;
        # with each class weighted to one half it would be p = 19/24 and w_rf = 0.8646
        w, steps, ok = cb.fit_weighted(cb.groups_of(*as_arrays(LEAN, copies=1)[:2], range(10)))
        self.assertTrue(ok)
        self.assertAlmostEqual(w[2], 0.875, places=6)
        for j in (0, 1, 3, 4):
            self.assertAlmostEqual(w[j], 0.125 / 4, places=6)
        self.assertAlmostEqual(math.fsum(w), 1.0, places=12)

    def test_the_weights_stay_on_the_simplex(self):
        got = cb.project_simplex([0.5, 0.5, 0.5, 0.0, 0.0])
        for a, b in zip(got, [1 / 3, 1 / 3, 1 / 3, 0.0, 0.0]):
            self.assertAlmostEqual(a, b, places=15)
        self.assertEqual(cb.project_simplex([2.0, 0.0, 0.0, 0.0, 0.0]), [1.0, 0.0, 0.0, 0.0, 0.0])
        w, _, ok = cb.fit_weighted(cb.groups_of(*as_arrays(LED, copies=1)[:2], range(27)))
        self.assertTrue(ok)
        self.assertEqual(w, [0.0, 0.0, 1.0, 0.0, 0.0])

    def test_the_logistic_fits_each_groups_share(self):
        S, y, folds = as_arrays(SATURATED)
        probs, recs = cb.crossfit("logistic", S, y, folds)
        self.assertTrue(all(r["converged"] for r in recs))
        shares = {s: p for s, p in zip(S, probs)}
        for s, share in zip([g[0] for g in SATURATED], (0.1, 0.75, 1 / 3, 2 / 3, 7 / 15, 0.6)):
            self.assertAlmostEqual(shares[s], share, places=6)
        # called at 0.5: the groups at 0.75, 2/3 and 0.6; the group at 7/15 is not called
        self.assertEqual(cb.confusion(y, [cb.is_call(p) for p in probs]),
                         {"tp": 40, "fn": 45, "fp": 20, "tn": 95, "called": 200})

    def test_a_separable_fit_does_not_converge(self):
        _, _, ok = cb.fit_logistic(cb.groups_of(*as_arrays(SEPARABLE, copies=1)[:2], range(10)))
        self.assertFalse(ok)


class Crossfit(Case):
    def test_each_fold_is_fitted_on_the_other_folds(self):
        S, y, folds = as_arrays(LED)
        for arm in cb.ARMS:
            _, recs = cb.crossfit(arm, S, y, folds)
            self.assertEqual([(r["fold"], r["trained"], r["tested"], r["shared"]) for r in recs],
                             [(f, 108, 27, 0) for f in range(5)])

    def test_gate_2_counts_no_shared_genome_on_the_standard_read(self):
        _, (rc, out) = self.read()
        self.assertEqual(rc, 0, out)
        self.assertEqual(gate_line(out, 2), "gate 2 cross-fitting: 60 fold fits, 0 whose training and held-out genomes "
                                            "overlap or do not add up to the report: pass")


class Counts(Case):
    def test_the_rf_led_report(self):
        self.read_write()
        r = self.rows()[(SPECIES[0], "drug a")]
        self.assertEqual({k: r[k] for k in ("genomes", "resistant", "susceptible", "mean_miss_by_1to4",
                                            "mean_miss_by_all5")},
                         {"genomes": "135", "resistant": "70", "susceptible": "65", "mean_miss_by_1to4": "65",
                          "mean_miss_by_all5": "5"})
        self.assertEqual([r[f"mean_{k}"] for k in cb.COUNTS], ["0", "70", "0", "65", "135"])
        for arm in cb.ARMS:
            self.assertEqual([r[f"{arm}_{k}"] for k in cb.COUNTS], ["45", "25", "5", "60", "135"])
            self.assertEqual([r[f"{arm}_{t}"] for t in cb.TRANSITIONS], ["45", "0", "0", "0", "5"])
            self.assertEqual(r[f"{arm}_folds_converged"], "5")
        # weighted: 6 genomes at -log 0.2 and 21 at -log 0.8 per copy of 27
        self.assertAlmostEqual(float(r["weighted_log_loss"]), (6 * math.log(5) + 21 * math.log(1.25)) / 27, places=10)

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
            self.assertIn(f"| {arm} | 6 | 6 | 0.2830 | 0.2830 to 0.2830 | -0.6429 | -0.6429 to -0.6429 | 6 | kept |",
                          out)


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

    def test_kept_and_no_difference_shown(self):
        # six reports gaining 103/364 in balanced accuracy each, the very major error rate down by 9/14 each
        self.read_write()
        rules = {r["arm"]: r for r in lo.read_tsv(os.path.join(self.record, "rule.tsv"))}
        for arm in cb.ARMS:
            self.assertEqual({k: rules[arm][k] for k in ("datasets", "species", "median_gain", "gain_lo", "gain_hi",
                                                         "vme_median", "vme_lo", "gain_at_least_margin", "reading")},
                             {"datasets": "6", "species": "6", "median_gain": "0.282967032967",
                              "gain_lo": "0.282967032967", "gain_hi": "0.282967032967", "vme_median": "-0.642857142857",
                              "vme_lo": "-0.642857142857", "gain_at_least_margin": "6", "reading": "kept"})

    def test_a_median_gain_of_0_is_no_difference_shown(self):
        # six symmetric reports at a gain of 0 and one at 103/364: the median gain is 0
        self.read_write(mixed_set())
        rules = {r["arm"]: r for r in lo.read_tsv(os.path.join(self.record, "rule.tsv"))}
        for arm in cb.ARMS:
            self.assertEqual((rules[arm]["median_gain"], rules[arm]["gain_at_least_margin"], rules[arm]["reading"]),
                             ("0.000000000000", "1", "no difference shown"))
        self.assertEqual(len(lo.read_tsv(os.path.join(self.record, "pooled.tsv"))), (1 + 7) * 2 * 5)


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

        _, (rc, out) = self.read(before=drop_one, expect=5)
        self.assertEqual(rc, 1, out)
        self.assertTrue(gate_line(out, 1).endswith(": STOP"), out)
        self.assertFalse(os.path.exists(self.reading))

    def test_gate_2_stops_a_fit_on_its_own_genomes(self):
        with mock.patch.object(cb, "training_genomes", lambda folds, f: list(range(len(folds)))):
            _, (rc, out) = self.read()
        self.assertEqual(rc, 1, out)
        self.assertTrue(gate_line(out, 2).endswith(": STOP"), out)

    def test_gate_3_stops_on_an_incumbent_that_differs_from_vote_tsv(self):
        def one_more_miss(d):
            path = os.path.join(self.second, "vote.tsv")
            with open(path, encoding="utf-8") as f:
                head = f.readline().rstrip("\n").split("\t")
            rows = lo.read_tsv(path)
            rows[0]["mean_miss"] = str(int(rows[0]["mean_miss"]) + 1)
            lo.write_tsv(path, head, rows)

        _, (rc, out) = self.read(before=one_more_miss)
        self.assertEqual(rc, 1, out)
        self.assertTrue(gate_line(out, 3).endswith(": STOP"), out)

    def test_gate_4_stops_on_weights_off_the_simplex(self):
        with mock.patch.object(cb, "fit_weighted", lambda groups: ([0.3] * 5, 1, True)):
            _, (rc, out) = self.read()
        self.assertEqual(rc, 1, out)
        self.assertTrue(gate_line(out, 4).endswith(": STOP"), out)

    def test_gate_5_stops_a_fit_that_does_not_converge(self):
        reports = kept_set() + [designed(SPECIES[6], "drug c", SEPARABLE, tag=7)]
        _, (rc, out) = self.read(reports)
        self.assertEqual(rc, 1, out)
        self.assertTrue(gate_line(out, 5).endswith(": STOP"), out)
        self.assertIn("35 logistic fits, 5 not converged", gate_line(out, 5))

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
                                       "logistic_transitions"}, set())


if __name__ == "__main__":
    unittest.main()
