"""Tests for scripts/input-comparison/second_reading.py: the calibration count, the panel-wide BH family held to
stats.finalize_sets, the firm rule, the sign tests behind a prediction, the arm table, and what the command refuses.
Every report here is J's own finalized records under a manifest cut down to the keys the reading takes.

usage: python3 -B -m unittest discover -s tests/input_comparison -t .
"""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests.input_comparison import test_run as base

second_reading = base.load("second_reading")
stats = second_reading.stats
# The retuned twin judges scm and the lasso, the seed twin judges rf, and nothing judges gbt.
TWINS = {"twin-t": ["scm", "lasso"], "twin-s": ["rf"]}


def contrast(dataset, arm, learner="rf", *, p=0.5, estimate=0.05, half=0.02, auroc=0.0, family="C",
             split="lineage"):
    """J's pending record for default -> arm: a 95% interval of estimate +/- half and a 90% interval of
    +/- 0.8 half. estimate None is an undefined contrast. Family D is descriptive and carries no p."""
    defined = estimate is not None
    r = dict.fromkeys(stats.RECORD_KEYS)
    r.update(id=f"{dataset}/{family}:{learner}:{split}:default->{arm}:log_loss", family=family, baseline="default",
             treatment=arm, learner=learner, split=split, endpoint="log_loss", dataset=dataset, margin=0.02,
             estimate=estimate, ci=[estimate - half, estimate + half] if defined else None,
             ci_tost=[estimate - 0.8 * half, estimate + 0.8 * half] if defined else None,
             ci_method="repeat_fold_t" if defined else "none", p=p if defined and family == "C" else None,
             undefined_reason=None if defined else "one unit: a t interval needs at least 2",
             auroc_delta=auroc if defined else None, adjust="bh" if family == "C" else "none", unit="repeat_fold",
             n_units=15, _pending=True)
    return r


def absent(dataset, arm, learner="rf"):
    """The id of a planned contrast J didn't produce."""
    return contrast(dataset, arm, learner)["id"]


def twin_ids(planned):
    """The planned ids judged on their own p, as run.pending_contrasts names them: the twins' family C ids."""
    return [pid for pid in planned if pid.split("->")[1].split(":")[0] in TWINS and "/C:" in pid]


def report(species, records, missing=(), *, down=(), unreliable=(), own_p=True, **manifest):
    """A report as K writes it, cut down to what the reading takes. `records` are finalized inside the species
    as run.py does it, the twins on their own p; own_p False puts the twins inside BH, as a manifest without
    `sentinel_adjust` says. `missing` are planned ids with no record, and `down` the (dataset, learner) pairs
    where the default has no lineage result."""
    planned = [r["id"] for r in records] + list(missing)
    final = stats.finalize(records, planned, dict.fromkeys(missing, "not produced"),
                           unadjusted=twin_ids(planned) if own_p else ())
    if own_p:
        manifest = {"sentinel_adjust": "none", **manifest}
    datasets = sorted({pid.split("/")[0] for pid in planned})
    runs = [{"arm": "default", "dataset": d, "learner": learner, "split": "lineage",
             "status": "unavailable" if (d, learner) in down else "ok"}
            for d in datasets for learner in ("scm", "lasso", "rf", "gbt")]
    return {"contrasts": final, "planned_contrasts": planned, "undefined_contrasts": [], "unreliable": list(unreliable),
            "not_produced": [{"id": pid, "reason": "not produced"} for pid in missing],
            "manifest": {"datasets": [{"id": d, "species": species} for d in datasets], "runs": runs,
                         "sentinel_arms": sorted(TWINS), "sentinel_learners": TWINS,
                         "source_commit": "c" * 40, "binary_sha256": "b" * 64, **manifest}}


def quiet_twins(datasets, learners=("scm", "rf")):
    """Twin contrasts that don't reject, for the learners each twin judges."""
    return [contrast(d, twin, learner, p=0.6, estimate=0.001) for d in datasets
            for twin, judged in TWINS.items() for learner in judged if learner in learners]


class CalibrationCases(unittest.TestCase):
    def test_a_twin_contrast_counts_when_it_has_a_decision_to_protect(self):
        records = [
            # scm: one rejection in three made, and a fourth the twin has no result for (the default has one).
            contrast("a1", "twin-t", "scm", p=0.01), contrast("a2", "twin-t", "scm"), contrast("a3", "twin-t", "scm"),
            # The seed twin doesn't judge scm, so its rejection there is no evidence about scm.
            contrast("a1", "twin-s", "scm", p=0.001),
            # rf: one rejection in three made; on a4 the default itself has no lineage result.
            contrast("a1", "twin-s", "rf", p=0.01), contrast("a2", "twin-s", "rf"), contrast("a3", "twin-s", "rf"),
            # lasso: one made, one undefined. gbt: no twin judges it.
            contrast("a1", "twin-t", "lasso"), contrast("a2", "twin-t", "lasso", estimate=None),
            contrast("a1", "cap", "gbt"),
            # Rule 10 reads a twin's lineage contrast on a dataset that decides nothing too. It is family D and
            # has no test, so it is no rejection, though the default ran there.
            contrast("a5", "twin-s", "rf", family="D"),
        ]
        missing = [absent("a4", "twin-t", "scm"), absent("a4", "twin-s", "rf")]
        marked = [{"species": "A", "learner": "scm", "sentinel": missing[0], "reason": "not produced",
                   "contrasts": [r["id"] for r in records[:3]]}]
        got = second_reading.read_species(report("A", records, missing, down=[("a4", "rf")], unreliable=marked))
        tails = {learner: c.pop("tail") for learner, c in got["calibration"].items()}
        self.assertEqual(got["calibration"], {
            "scm": {"twins": ["twin-t"], "rejections": 2, "of": 4, "flagged": True, "left_out": 0, "no_test": 1,
                    "rule10_contrasts": 3},
            "rf": {"twins": ["twin-s"], "rejections": 1, "of": 3, "flagged": False, "left_out": 1, "no_test": 0,
                   "rule10_contrasts": 0},
            "lasso": {"twins": ["twin-t"], "rejections": 1, "of": 2, "flagged": False, "left_out": 0, "no_test": 1,
                      "rule10_contrasts": 0},
            "gbt": {"twins": [], "rejections": 0, "of": 0, "flagged": None, "left_out": 0, "no_test": 0,
                    "rule10_contrasts": 0},
        })
        self.assertAlmostEqual(tails["scm"], 1 - 0.95 ** 4 - 4 * 0.05 * 0.95 ** 3, places=12)  # 0.0140
        self.assertAlmostEqual(tails["rf"], 1 - 0.95 ** 3, places=12)  # 0.1426
        self.assertAlmostEqual(tails["lasso"], 1 - 0.95 ** 2, places=12)
        self.assertIsNone(tails["gbt"])

    def test_one_rejection_of_one_contrast_is_not_flagged(self):
        # P(1 or more | 1, 0.05) is 0.05 itself, and the flag needs a tail under 0.05.
        got = second_reading.read_species(report("A", [contrast("a1", "twin-s", "rf", p=0.01)]))["calibration"]["rf"]
        self.assertEqual((got["rejections"], got["of"], got["flagged"]), (1, 1, False))
        self.assertAlmostEqual(got["tail"], 0.05, places=12)
        two = [contrast(d, "twin-s", "rf", p=0.01) for d in ("a1", "a2")]
        self.assertIs(second_reading.read_species(report("A", two))["calibration"]["rf"]["flagged"], True)

    def test_a_twin_inside_bh_is_still_counted_on_its_own_p(self):
        # Three contrasts in one BH family: the twin's 0.03 adjusts to 0.09, so the report leaves it undecided.
        records = [contrast("a1", "twin-s", "rf", p=0.03), contrast("a1", "cap", p=0.5), contrast("a2", "cap", p=0.6)]
        rep = report("A", records, own_p=False)
        twin = rep["contrasts"][0]
        self.assertEqual((twin["adjust"], twin["decision"]), ("bh", "inconclusive"))
        self.assertAlmostEqual(twin["p_adj"], 0.09, places=12)
        got = second_reading.read([rep])
        count = got["species"][0]["calibration"]["rf"]
        self.assertEqual((count["rejections"], count["of"]), (1, 1))
        self.assertEqual(got["panel"]["family"], {"tested": 3, "entered_at_p_1": 0, "size": 3})

    def test_without_sentinel_learners_every_twin_judges_every_learner(self):
        rep = report("A", [contrast("a1", "twin-s", "scm", p=0.001), contrast("a1", "twin-t", "scm")])
        del rep["manifest"]["sentinel_learners"]
        self.assertEqual(second_reading.read_species(rep)["calibration"]["scm"]["twins"], ["twin-s", "twin-t"])
        self.assertEqual(second_reading.read_species(rep)["calibration"]["scm"]["rejections"], 1)

    def test_a_twin_the_manifest_gives_no_learner_judges_none(self):
        rep = report("A", [contrast("a1", "twin-s", "rf", p=0.001)], sentinel_learners={"twin-t": ["scm"]})
        got = second_reading.read_species(rep)["calibration"]["rf"]
        self.assertEqual((got["twins"], got["of"], got["flagged"]), ([], 0, None))

    def test_the_counts_add_over_the_species(self):
        a = report("A", [contrast("a1", "twin-s", "rf", p=0.01), contrast("a2", "twin-s", "rf")])
        b = report("B", [contrast(f"b{i}", "twin-s", "rf") for i in range(8)] + [contrast("b1", "twin-t", "scm")])
        got = second_reading.read([a, b])["panel"]["calibration"]
        self.assertEqual({k: got["rf"][k] for k in ("rejections", "of", "flagged")},
                         {"rejections": 1, "of": 10, "flagged": False})
        self.assertAlmostEqual(got["rf"]["tail"], 1 - 0.95 ** 10, places=12)
        low, high = got["rf"]["ci95"]  # the exact interval for 1 of 10
        self.assertAlmostEqual(low, 0.0025286, places=6)
        self.assertAlmostEqual(high, 0.4450161, places=6)
        self.assertEqual((got["scm"]["rejections"], got["scm"]["of"]), (0, 1))


class PanelCases(unittest.TestCase):
    def sets(self):
        """Two species. A's 0.04 is decided inside A's family of 6 (0.04 * 6 / 5 = 0.048) and lost in the panel's
        family of 9 (0.04 * 9 / 6 = 0.06)."""
        a = [contrast("a1", "cap", p=0.001), contrast("a2", "cap", p=0.04), contrast("a3", "cap", p=0.2),
             contrast("a1", "cap", "scm", p=0.001), contrast("a1", "cap", "gbt", p=0.001),
             contrast("a3", "bud", p=0.002, estimate=-0.05)]
        a += [contrast(d, "twin-t", "scm", p=0.01) for d in ("a1", "a2")] + quiet_twins(["a1", "a2", "a3"], ("rf",))
        b = [contrast("b1", "cap", p=0.004), contrast("b2", "cap", p=0.5)] + quiet_twins(["b1"])
        b += [contrast("b1", "cap", split="random", family="D")]  # descriptive: in no BH family
        return [(a, []), (b, [absent("b3", "cap")])]

    def test_the_panel_family_is_finalize_sets_over_the_species(self):
        sets = self.sets()
        reports = [report(species, records, missing) for species, (records, missing) in zip("AB", sets)]
        panel = second_reading.read(reports)["panel"]
        planned = [[r["id"] for r in records] + missing for records, missing in sets]
        joint = stats.finalize_sets(
            [(records, ids, dict.fromkeys(missing, "not produced")) for (records, missing), ids in zip(sets, planned)],
            unadjusted=twin_ids(planned[0] + planned[1]))
        expected = {r["id"]: r for one in joint for r in one if r["adjust"] == "bh"}
        inside = {x["id"]: x["decision"] for rep in reports for x in rep["contrasts"]}
        decided = [pid for pid, r in expected.items() if {inside[pid], r["decision"]} != {"inconclusive"}]
        self.assertEqual([d["id"] for d in panel["decisions"]], decided)
        self.assertEqual({d["id"]: (d["p_adj_panel"], d["in_panel"]) for d in panel["decisions"]},
                         {pid: (expected[pid]["p_adj"], expected[pid]["decision"]) for pid in decided})
        self.assertEqual(panel["family"], {"tested": 8, "entered_at_p_1": 1, "size": 9})
        # The fixture separates the two families: a2's contrast is improves inside A and undecided in the panel.
        lost = absent("a2", "cap")
        self.assertEqual((inside[lost], expected[lost]["decision"]), ("improves", "inconclusive"))

    def test_a_firm_decision_holds_in_both_families_on_an_unflagged_learner(self):
        reports = [report(species, records, missing) for species, (records, missing) in zip("AB", self.sets())]
        panel = second_reading.read(reports)["panel"]
        firm = {d["id"]: (d["learner_flagged"], d["firm"]) for d in panel["decisions"]}
        self.assertEqual(firm, {
            absent("a1", "cap"): (False, True),             # rf, improves in both
            absent("a2", "cap"): (False, False),            # rf, decided inside A only
            absent("a1", "cap", "scm"): (True, False),      # scm's twin rejected 2 of 2
            absent("a1", "cap", "gbt"): (None, False),      # no twin judges gbt
            absent("a3", "bud"): (False, True),             # rf, worse in both
            absent("b1", "cap"): (False, True),
        })
        self.assertEqual(panel["firm"], {"improves": 2, "worse": 1})
        self.assertEqual((panel["one_family_only"], panel["flagged_or_unjudged"]), (1, 2))


def two_species():
    """X: cap is better on five datasets by the mean over three learners (0.01), though two of the three are below
    zero; bud is worse on four and exactly zero on the fifth. Y: cap is worse on its one dataset."""
    x = [contrast(f"x{i}", "cap", learner, estimate=e) for i in range(1, 6)
         for learner, e in (("rf", 0.05), ("scm", -0.01), ("lasso", -0.01))]
    x += [contrast(f"x{i}", "bud", estimate=-0.02 if i < 5 else 0.0) for i in range(1, 6)]
    y = [contrast("y1", "cap", learner, estimate=-0.05) for learner in ("rf", "scm")]
    return [report("X", x, [absent("x1", "none")]), report("Y", y)]


class SignTestCases(unittest.TestCase):
    def test_a_prediction_is_a_one_sided_sign_test_over_datasets_of_the_named_species(self):
        reports = two_species()
        got = second_reading.read(reports, [("cap", "better")], ["X"])["panel"]["predictions"]
        self.assertEqual(got, {"species": ["X"], "pass": True, "tests": [
            {"arm": "cap", "predicted": "better", "better": 5, "worse": 0, "unchanged": 0, "p_one_sided": 0.5 ** 5}]})
        # With Y's dataset the same prediction fails: 5 of 6, p = 7 / 64.
        everywhere = second_reading.read(reports, [("cap", "better")])["panel"]["predictions"]
        self.assertEqual((everywhere["species"], everywhere["pass"]), (["X", "Y"], False))
        self.assertAlmostEqual(everywhere["tests"][0]["p_one_sided"], 7 / 64, places=12)
        # A tie is dropped, so 4 of 4 worse gives 1 / 16, and one failed prediction fails the set.
        both = second_reading.read(reports, [("cap", "better"), ("bud", "worse")], ["X"])["panel"]["predictions"]
        self.assertEqual(both["tests"][1], {"arm": "bud", "predicted": "worse", "better": 0, "worse": 4,
                                            "unchanged": 1, "p_one_sided": 0.5 ** 4})
        self.assertFalse(both["pass"])

    def test_the_by_dataset_sign_test_is_two_sided(self):
        result = second_reading.read(two_species())
        x = result["species"][0]["by_dataset"]
        self.assertEqual({k: x["cap"][k] for k in ("datasets", "better", "worse", "unchanged", "p")},
                         {"datasets": 5, "better": 5, "worse": 0, "unchanged": 0, "p": 2 * 0.5 ** 5})
        self.assertAlmostEqual(x["cap"]["median"], 0.01, places=12)
        # The test takes the larger side, whichever it is: bud is worse on 4 of 4, with one tie dropped.
        self.assertEqual(x["bud"], {"datasets": 5, "better": 0, "worse": 4, "unchanged": 1, "median": -0.02,
                                    "p": 2 * 0.5 ** 4})
        self.assertEqual(x["none"], {"datasets": 0, "better": 0, "worse": 0, "unchanged": 0, "median": None,
                                     "p": None})
        whole = result["panel"]["by_dataset"]["cap"]
        self.assertEqual((whole["datasets"], whole["better"], whole["worse"]), (6, 5, 1))
        self.assertAlmostEqual(whole["p"], 2 * 7 / 64, places=12)
        # One each way is no evidence: twice the tail is over 1, and p is 1.
        self.assertEqual(second_reading.sign_test([0.1, -0.2, 0.0]),
                         {"datasets": 3, "better": 1, "worse": 1, "unchanged": 1, "median": 0.0, "p": 1.0})


class ArmTableCases(unittest.TestCase):
    def test_the_arm_row_counts_the_size_of_the_change(self):
        records = [contrast("a1", "cap", estimate=0.0, half=0.0125, auroc=0.01),   # zero, 90% interval +/- 0.01
                   contrast("a2", "cap", estimate=0.03, half=0.02, auroc=0.03, p=0.001),    # beyond, improves
                   contrast("a3", "cap", estimate=-0.02, half=0.001, auroc=0.02, p=0.002),  # at the margin, worse
                   # 90% interval 0.002 to 0.018, inside; the 95% one ends at the margin itself
                   contrast("a4", "cap", estimate=0.01, half=0.01, auroc=0.02),
                   contrast("a5", "cap", estimate=None),
                   contrast("a1", "cap", estimate=0.5, split="random", family="D"),  # descriptive: not in the table
                   contrast("a1", "twin-s", estimate=0.07, auroc=0.04, p=0.001)]     # an arm with one contrast
        got = second_reading.read_species(report("A", records, [absent("a6", "cap"), absent("a1", "bud")]))
        row = got["arms"]["cap"]
        for key, value in (("q1", -0.005), ("median", 0.005), ("q3", 0.015), ("auroc_median", 0.02)):
            self.assertAlmostEqual(row.pop(key), value, places=12)
        self.assertEqual(row, {"twin": False, "planned": 6, "made": 4, "improves": 1, "worse": 1, "inside": 2,
                               "beyond": 2, "zero": 1})
        self.assertEqual(got["arms"]["twin-s"], {"twin": True, "planned": 1, "made": 1, "improves": 1, "worse": 0,
                                                 "inside": 0, "beyond": 1, "zero": 0, "q1": 0.07, "median": 0.07,
                                                 "q3": 0.07, "auroc_median": 0.04})
        self.assertEqual(got["arms"]["bud"], {"twin": False, "planned": 1, "made": 0, "improves": 0, "worse": 0,
                                              "inside": 0, "beyond": 0, "zero": 0, "q1": None, "median": None,
                                              "q3": None, "auroc_median": None})
        self.assertEqual((got["planned"], got["made"], got["datasets"], got["unanswered"]), (8, 5, 6, ["a5", "a6"]))


class TextCases(unittest.TestCase):
    def test_the_text_puts_each_number_beside_its_label(self):
        # No two numbers here are equal, so one printed in another's place changes the text.
        species = {
            "species": "X", "datasets": 31, "unanswered": ["x1", "x2"], "planned": 33, "made": 32,
            "arms": {"cap": {"twin": False, "planned": 11, "made": 12, "improves": 13, "worse": 14, "inside": 15,
                             "beyond": 16, "zero": 17, "q1": 0.18, "median": 0.19, "q3": 0.2, "auroc_median": 0.21}},
            "by_dataset": {"cap": {"datasets": 22, "better": 23, "worse": 24, "unchanged": 3, "median": 0.26,
                                   "p": 0.25}},
            "calibration": {"rf": {"twins": ["twin-s"], "rejections": 4, "of": 27, "tail": 0.28, "flagged": False,
                                   "left_out": 29, "no_test": 30, "rule10_contrasts": 34},
                            "gbt": {"twins": [], "rejections": 0, "of": 0, "tail": None, "flagged": None,
                                    "left_out": 0, "no_test": 0, "rule10_contrasts": 35}}}
        panel = {
            "family": {"tested": 41, "entered_at_p_1": 42, "size": 83}, "firm": {"improves": 43, "worse": 44},
            "one_family_only": 45, "flagged_or_unjudged": 46,
            "decisions": [{"id": "x1/C:rf:lineage:default->cap:log_loss", "species": "X", "in_species": "improves",
                           "in_panel": "inconclusive", "learner_flagged": False, "firm": False, "estimate": 0.47,
                           "p_adj_species": 0.048, "p_adj_panel": 0.06, "auroc_delta": 0.49},
                          {"id": "x2/C:rf:lineage:default->cap:log_loss", "species": "X", "in_species": "worse",
                           "in_panel": "worse", "learner_flagged": False, "firm": True, "estimate": -0.5,
                           "p_adj_species": 0.00012345, "p_adj_panel": 0.0006789, "auroc_delta": None}],
            "calibration": {"rf": {"rejections": 5, "of": 51, "tail": 0.52, "flagged": False, "ci95": [0.053, 0.54]},
                            "gbt": {"rejections": 0, "of": 0, "tail": None, "flagged": None, "ci95": None}},
            "by_dataset": {"cap": {"datasets": 55, "better": 56, "worse": 57, "unchanged": 6, "median": 0.58,
                                   "p": 0.59}},
            "predictions": {"species": ["X", "Y"], "pass": False, "tests": [
                {"arm": "cap", "predicted": "better", "better": 61, "worse": 62, "unchanged": 63,
                 "p_one_sided": 0.64}]}}
        self.assertEqual(second_reading.render({"species": [species], "panel": panel}).split("\n"), [
            "X: 31 datasets, 2 the test can't answer; 32 of 33 family C contrasts made",
            "arm                    plan made impr worse inside beyond zero       q1   median       q3    auroc |"
            " datasets better worse  sign p",
            "cap                      11   12   13    14     15     16   17   0.1800   0.1900   0.2000   0.2100 |"
            "       22     23    24  0.2500",
            "calibration rf     twin-s: 4 of 27 (tail 0.2800, flagged False); 29 left out, 30 counted with no test;"
            " rule 10 marks 34 contrasts",
            "calibration gbt    no twin: 0 of 0 (tail -, flagged None); 0 left out, 0 counted with no test;"
            " rule 10 marks 35 contrasts",
            "",
            "panel: one BH family of 83 (41 tested, 42 at p = 1)",
            "firm decisions: 43 improves, 44 worse; 45 decided under one family only; 46 on a flagged or unjudged"
            " learner",
            "       species improves     panel inconclusive x1/C:rf:lineage:default->cap:log_loss estimate 0.4700"
            " p_adj 0.048 / 0.06 auroc 0.4900 learner flagged False",
            "  firm species worse        panel worse        x2/C:rf:lineage:default->cap:log_loss estimate -0.5000"
            " p_adj 0.0001234 / 0.0006789 auroc - learner flagged False",
            "calibration rf     over the species: 5 of 51 (tail 0.5200, flagged False), 95% interval 0.0530 to 0.5400",
            "calibration gbt    over the species: 0 of 0 (tail -, flagged None)",
            "by dataset cap                     56 better  57 worse of  55, median 0.5800, sign p 0.5900",
            "predictions on X, Y: pass False",
            "  cap predicted better: 61 better, 62 worse, 63 unchanged, one-sided p 0.6400"])

    def test_without_a_prediction_the_text_ends_at_the_sign_tests(self):
        text = second_reading.render(second_reading.read(two_species()))
        self.assertTrue(text.endswith("by dataset none                     0 better   0 worse of   0, median -, "
                                      "sign p -"), text[-200:])


class RefusalCases(unittest.TestCase):
    def test_reports_that_are_not_one_comparison_are_refused(self):
        a, b = report("A", [contrast("a1", "cap")]), report("B", [contrast("b1", "cap")])
        mixed = report("A", [contrast("a1", "cap"), contrast("b1", "cap")])
        mixed["manifest"]["datasets"][1]["species"] = "B"
        cases = [([mixed], {}, "must hold one species, and this one holds A, B"),
                 ([a, a], {}, "these have more: A"),
                 ([a, report("B", [contrast("b1", "cap")], source_commit="d" * 40)], {}, "differ in source_commit"),
                 ([a, report("B", [contrast("b1", "cap")], binary_sha256="e" * 64)], {}, "differ in binary_sha256"),
                 ([a, report("B", [contrast("b1", "cap")], synthetic=True)], {}, "differ in synthetic"),
                 ([a, b], {"predictions": [("cap", "better")], "predict_on": ["C"]}, "no report: C"),
                 ([a, b], {"predictions": [("cup", "better")], "predict_on": ["A"]}, "names arm cup")]
        for reports, more, message in cases:
            with self.subTest(message=message), self.assertRaisesRegex(SystemExit, message):
                second_reading.read(reports, **more)


class CommandCases(unittest.TestCase):
    def test_the_command_reads_validated_reports_and_writes_the_reading(self):
        x, y = two_species()
        loaded = {"x.json": x, "y.json": y}
        with tempfile.TemporaryDirectory() as tmp, \
                mock.patch.object(second_reading.k_report, "load_report", side_effect=loaded.__getitem__) as load, \
                contextlib.redirect_stdout(io.StringIO()) as out:
            dest = Path(tmp) / "reading.json"
            second_reading.main(["x.json", "y.json", "--predict", "cap=better", "--predict-on", "X",
                                 "--json", str(dest)])
            written = dest.read_text()
        self.assertEqual([c.args for c in load.call_args_list], [("x.json",), ("y.json",)])
        expected = second_reading.read([x, y], [("cap", "better")], ["X"])
        self.assertEqual(written, json.dumps(expected, indent=1, sort_keys=True))
        self.assertEqual(out.getvalue(), second_reading.render(expected) + "\n")

    def test_the_command_refuses_a_prediction_it_cannot_read(self):
        for argv, message in ((["x.json", "--predict", "cap=up"], "is not ARM=better or ARM=worse"),
                              (["x.json", "--predict-on", "X"], "--predict-on needs a --predict")):
            with self.subTest(argv=argv), contextlib.redirect_stderr(io.StringIO()) as err, \
                    self.assertRaises(SystemExit) as stop:
                second_reading.main(argv)
            self.assertEqual(stop.exception.code, 2)
            self.assertIn(message, err.getvalue())


if __name__ == "__main__":
    unittest.main()
