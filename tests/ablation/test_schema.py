"""Validators of the ablation-report schema against the ablation data contract
v1.5: good records pass, each malformed field fails, null shapes, the decision
rule and claim levels hold.

usage: python3 -I -B -m unittest discover -s tests/ablation -t .
K's mirrored constants are checked against J's stats.py: the tree's own
scripts/ablation/stats.py, or the copy GRM_TASK_J_STATS names.
The good_* builders are shared with test_report.py.
"""
import copy
import datetime
import importlib.util
import os
import sys
import unittest
from pathlib import Path

from scripts.ablation import schema
from scripts.ablation.schema import SchemaError

ROOT = Path(__file__).resolve().parents[2]
H1, H2, H3, H4 = ("1" * 64, "2" * 64, "3" * 64, "4" * 64)
LINEAGE_PLANS, RANDOM_PLANS = ["a" * 64, "c" * 64], ["b" * 64, "d" * 64]
METHOD = "cluster_jackknife_t"


def good_block(**over):
    return {"id": "K", "kind": "Anonymous", "version": "k=31", "source_sha256": H1, "priority": 5,
            "screen_cap": 4000, "catalog_date": None, "circularity_risk": False, "unavailable_reason": None,
            **over}


def good_dataset(**over):
    return {"id": "ecoli__amp", "species": "Escherichia coli", "partition": "confirmatory",
            "n_samples": 400, "group_counts": {"lineage": 12}, **over}


def good_run(**over):
    return {"id": "r1", "dataset": "ecoli__amp", "arm": "K", "blocks": ["K"],
            "learner": "lasso", "split": "lineage", "repeats": 2, "fold_plan_sha256": list(LINEAGE_PLANS),
            "status": "ok", "reason": None, "metrics": {"balanced_accuracy": 0.81, "auroc": None},
            "fit_count": 10, "seconds": 12.5, "predictions": "predictions/r1.jsonl", **over}


def good_manifest():
    """Four datasets over three species (one development), H-shaped blocks K,
    C and an unavailable P, and runs on a lineage and a random split, two of
    them unavailable."""
    datasets = [
        good_dataset(),
        good_dataset(id="ecoli__cip", n_samples=380),
        good_dataset(id="kpneu__amp", species="Klebsiella pneumoniae", n_samples=150,
                     group_counts={"lineage": 5}),
        good_dataset(id="mtb__inh", species="Mycobacterium tuberculosis", partition="development",
                     n_samples=256, group_counts={"lineage": 4, "sublineage": 9}),
    ]
    blocks = [
        good_block(),
        good_block(id="C", kind="CuratedExact", version="2026-09-15 catalog", source_sha256=H2, priority=0,
                   screen_cap=None, catalog_date="2026-09-15", circularity_risk=True),
        good_block(id="P", kind="Presence", version="", source_sha256="", priority=2, screen_cap=None,
                   unavailable_reason="replicon database not licensed for this run"),
    ]
    runs, n = [], 0
    for d in ("ecoli__amp", "ecoli__cip", "kpneu__amp", "mtb__inh"):
        for split, plans in (("lineage", LINEAGE_PLANS), ("random", RANDOM_PLANS)):
            for arm, arm_blocks in (("K", ["K"]), ("K+C", ["K", "C"])):
                n += 1
                runs.append(good_run(id=f"r{n}", dataset=d, arm=arm, blocks=arm_blocks, split=split,
                                     fold_plan_sha256=list(plans), predictions=f"predictions/r{n}.jsonl"))
    for u in runs:
        if u["dataset"] == "kpneu__amp" and u["arm"] == "K+C" and u["split"] == "lineage":
            u.update(status="unavailable", reason="one-class test part in fold 2", metrics={},
                     fit_count=0, seconds=0.4, predictions=None)
    runs.append(good_run(id="r-plasmid", arm="K+P", blocks=["K", "P"], status="unavailable",
                         reason="block P unavailable", metrics={}, fit_count=0, seconds=0, predictions=None,
                         fold_plan_sha256=None))
    return {
        "schema_version": 1,
        "protocol_sha256": H3,
        "source_commit": "c90ebaf48998dc244806ac435c3706bfd87c9a9b",
        "binary_sha256": H4,
        "blocks": blocks,
        "datasets": datasets,
        "missingness": [
            {"dataset": "ecoli__amp", "block": "C", "rate_R": 0.02, "rate_S": 0.11,
             "rate_diff_p": 0.02, "indicator_auroc": 0.36, "indicator_auroc_oriented": 0.64, "flag": True},
            {"dataset": "mtb__inh", "block": "K", "rate_R": 0.0, "rate_S": 0.0, "rate_diff_p": None,
             "indicator_auroc": None, "indicator_auroc_oriented": None, "flag": False},
        ],
        "group_proxy_audit": [
            {"dataset": "ecoli__amp", "block": "K", "metric": "lineage_auroc", "value": 0.97, "flag": True},
        ],
        "runs": runs,
        "compute": {"planned_fit_units": 270, "actual_seconds": 400.0, "tiers_completed": ["smoke", "pilot"]},
    }


def key(family="C", baseline="K", treatment="K+C", endpoint="log_loss", learner="lasso", split="lineage"):
    return f"{family}:{learner}:{split}:{baseline}->{treatment}:{endpoint}"


def good_contrast(**over):
    """Family C on the primary endpoint, decision improves under the contract rule.
    The id is built from the record's dataset and key fields unless given; a given
    id sets the dataset and key fields it implies, unless they are also given."""
    rec = {
        "dataset": "ecoli__amp", "family": "C", "baseline": "K", "treatment": "K+C",
        "learner": "lasso", "split": "lineage", "endpoint": "log_loss",
        "estimate": 0.031, "ci": [0.012, 0.050], "ci_tost": [0.015, 0.047], "ci_method": METHOD,
        "unit": "cluster", "n_units": 64, "few_clusters": False, "p": 0.004, "p_adj": 0.012,
        "adjust": "bh", "margin": 0.02, "decision": "improves", "claim_level": "predictive_increment",
        "undefined_resamples": 0, "undefined_reason": None, "interval_label": "95% lineage-cluster jackknife-t CI",
        "pairing": "same_plan", "auroc_delta": 0.012, "auroc_ci": [-0.004, 0.028],
    }
    if "id" in over:
        dataset, _, rest = over["id"].partition("/")
        parts = schema.parse_key(rest) or {}
        rec.update({"dataset": dataset, **parts})
        if parts.get("family") in schema.FAMILY_ADJUST:
            rec["adjust"] = schema.FAMILY_ADJUST[parts["family"]]
    if isinstance(over.get("family"), str) and "adjust" not in over:
        rec["adjust"] = schema.FAMILY_ADJUST.get(over["family"], rec["adjust"])
    rec.update(over)
    if rec["endpoint"] != "log_loss" and "auroc_delta" not in over:  # J carries it beside log_loss only
        rec.update(auroc_delta=None, auroc_ci=None)
    rec.setdefault("id", f"{rec['dataset']}/" + key(*(rec[k] for k in ("family", "baseline", "treatment",
                                                                       "endpoint", "learner", "split"))))
    return {k: rec[k] for k in ["id", *[k for k in rec if k != "id"]]}


# Interval and test shapes per decision, valid in family C.
SHAPES = {
    "improves": {"estimate": 0.03, "ci": [0.01, 0.05], "ci_tost": [0.015, 0.045], "p": 0.001, "p_adj": 0.002},
    "improves_small": {"estimate": 0.01, "ci": [0.004, 0.016], "ci_tost": [0.005, 0.015], "p": 0.001,
                       "p_adj": 0.002},
    "worse": {"estimate": -0.03, "ci": [-0.05, -0.01], "ci_tost": [-0.045, -0.015], "p": 0.001, "p_adj": 0.002},
    "inside_margin": {"estimate": 0.0, "ci": [-0.01, 0.01], "ci_tost": [-0.008, 0.008], "p": 0.5, "p_adj": 0.6},
    "inconclusive": {"estimate": 0.0, "ci": [-0.03, 0.03], "ci_tost": [-0.025, 0.025], "p": 0.5, "p_adj": 0.6},
    "inconclusive_off_zero": {"estimate": 0.03, "ci": [0.01, 0.05], "ci_tost": [0.015, 0.045],
                              "p": 0.04, "p_adj": 0.08},
}
UNDEFINED = {"estimate": None, "ci": None, "ci_tost": None, "p": None, "p_adj": None, "decision": "inconclusive",
             "undefined_reason": "vme undefined: no resistant (y=1) samples"}


def decision_of(shape, family):
    s = SHAPES[shape]
    return schema.expected_decision(s["ci"], None if family == "D" else s["p_adj"])


def shaped(shape, family="C", **over):
    values = dict(SHAPES[shape])
    if family == "D":
        values.update(p=None, p_adj=None)
    return good_contrast(**{"family": family, "adjust": schema.FAMILY_ADJUST[family],
                            "decision": decision_of(shape, family), **values, **over})


def good_pooled(**over):
    return {
        "contrast": key("A", endpoint="r2_ll"), "family": "A", "unit": "species", "n_units": 3, "n_datasets": 5,
        "n_excluded": 1,
        "excluded": [{"dataset": "kpneu__amp", "reason": "K+C run unavailable: one-class test part"}],
        "estimate": 0.02, "ci": [-0.004, 0.044], "ci_tost": [0.0, 0.04], "p": 0.25, "p_adj": 0.5,
        "adjust": "holm", "margin": 0.02, "decision": "inconclusive",
        "interval_label": "95% species-level jackknife-t CI", **over,
    }


def all_excluded_pooled(**over):
    return good_pooled(**{"n_units": 0, "n_datasets": 1, "estimate": None, "ci": None, "ci_tost": None,
                          "p": None, "p_adj": None, **over})


def good_feature(**over):
    return {"dataset": "ecoli__amp", "column": "blaTEM-1", "blocks": ["C"],
            "selection_freq": 0.93, "sign_stable": True, "claim_level": "association_stable",
            "evidence_refs": [], **over}


def good_prediction(**over):
    return {"dataset": "ecoli__amp", "arm": "K", "learner": "lasso", "split": "lineage", "repeat": 0,
            "sample_id": "SAMN0001", "y": 1, "fold": 2, "score": 0.83, "called": True, "cluster": "L1",
            "fold_plan_sha256": LINEAGE_PLANS[0], **over}


def note(**over):
    return {"by": "Fixture Annotator", "date": "2026-10-30", "level": "causal_supported",
            "previous_level": "association_stable", **over}


MISSING = object()
NAN, INF = float("nan"), float("inf")
HUGE = 10 ** 400
# Wrong values common to every field of a given type.
NOT_TEXT = [None, "", "  ", 3, ["x"]]
NOT_ID = NOT_TEXT + ["has space", "a/b", "-lead"]
NOT_BLOCK_ID = NOT_TEXT + ["has space", "a/b", "-lead", "1K", "K.v2", "K+C", "_x"]
NOT_INT0 = [None, -1, 1.5, True, "3", 2 ** 60]
NOT_NUM = [NAN, INF, -INF, True, "0.1", HUGE]
NOT_PROB = NOT_NUM + [-0.01, 1.01]
NOT_BOOL = [None, 0, 1, "true"]
NOT_HEX = [None, "", "A" * 64, "1" * 63, "1" * 65, "g" * 64]

PREDICTION_BAD = {
    "dataset": NOT_ID, "arm": NOT_TEXT, "learner": NOT_TEXT, "split": NOT_TEXT,
    "repeat": NOT_INT0, "sample_id": NOT_TEXT, "y": [None, 2, -1, True, False, 1.0, "1"],
    "fold": NOT_INT0, "score": NOT_PROB + [None], "called": NOT_BOOL,
    "cluster": [None, "", True, -1, 1.5, ["L1"]], "fold_plan_sha256": NOT_HEX,
}
CONTRAST_BAD = {
    "id": NOT_TEXT + ["ecoli__amp", "ecoli__amp/", "/" + key(), "bad id/" + key(), "ecoli__amp/A-C",
                      "ecoli__amp/" + key(endpoint="log loss"), "ecoli__amp/" + key(endpoint="log/loss"),
                      "ecoli__amp/" + key(endpoint="log`loss"), "ecoli__amp/" + key(learner="las:so"),
                      "ecoli__amp/E:lasso:lineage:K->K+C:log_loss", "ecoli__amp/C:lasso:lineage:K+C:log_loss",
                      "ecoli__amp/C:lasso:lineage:K->K->C:log_loss", "ecoli__amp/C:lasso:lineage:->K+C:log_loss",
                      "ecoli__amp/C:lasso:lineage:K->:log_loss", "ecoli__amp/C::lineage:K->K+C:log_loss",
                      "ecoli__amp/C:lasso:lineage:K->K+C:"],
    "dataset": NOT_ID + ["ecoli__cip"],
    "family": [None, "E", "a", "", 1, "A", "B", "D"], "baseline": NOT_TEXT + ["K+C"],
    "treatment": NOT_TEXT + ["K"], "learner": NOT_TEXT + ["rf"], "split": NOT_TEXT + ["random"],
    "endpoint": NOT_TEXT + ["auroc"], "estimate": NOT_NUM + [None],
    "ci": [None, [], [0.1], [0.1, 0.2, 0.3], [0.2, 0.1], [NAN, 0.1], [0.0, INF], [True, 1.0], "0,1", [0, HUGE],
           ],
    "ci_tost": [None, [0.2, 0.1], [0.1], "x", [-0.5, 0.5], [0.011, 0.047], [0.015, 0.051]],
    "ci_method": NOT_TEXT + ["cluster bootstrap, percentile, 95%, B=2000", "Cluster_jackknife_t",
                             "cluster_jackknife_t ", "none", "cluster_bootstrap_percentile_B300"],
    "unit": NOT_TEXT, "n_units": NOT_INT0 + [0], "few_clusters": NOT_BOOL + [True],
    "p": NOT_PROB + [None], "p_adj": NOT_PROB + [None, 0.001, 0.004 - 1e-13],
    "adjust": [None, "bonferroni", "Holm", "", "holm", "none"],
    "margin": [None, NAN, -0.01, 0, True, "0.02"],
    "decision": [None, "better", "no_value", "", "inconclusive", "worse", "equivalent"],
    "claim_level": [None, "causal", "", "Predictive_increment", "association_stable", "causal_supported"],
    "undefined_resamples": NOT_INT0 + [1], "undefined_reason": ["", 3, []], "interval_label": NOT_TEXT,
    "pairing": NOT_TEXT + ["nearest", "Same_plan", "same plan"],
    "auroc_delta": NOT_NUM, "auroc_ci": [[], [0.1], [0.2, 0.1], "x", [0.02, 0.03]],
}
POOLED_BAD = {
    "contrast": NOT_TEXT + ["A C", "A/C", "A-C", key("A", endpoint="log_loss"), key("A", endpoint="r2`ll"),
                           key("A", endpoint="r2_ll")[:-6]],
    "family": [None, "E", "a", "B", "C", "D"], "unit": [None, "dataset", "Species"],
    "n_units": NOT_INT0, "n_datasets": NOT_INT0, "n_excluded": NOT_INT0,
    "excluded": [None, "kpneu__amp", [{"dataset": "kpneu__amp"}], [{"dataset": "kpneu__amp", "reason": ""}],
                 [{"dataset": "x y", "reason": "r"}],
                 [{"dataset": "kpneu__amp", "reason": "r", "extra": 1}]],
    "estimate": [NAN, True, "0.1", None, HUGE], "ci": [[0.2, 0.1], [0.1], "x", None],
    "ci_tost": [[0.1], None, [-0.5, 0.5], [0.0, 0.045]],
    "p": [NAN, 1.5, -0.1, None], "p_adj": [NAN, 1.5, -0.1, None, 0.01], "adjust": [None, "bh", "none"],
    "margin": [None, 0, -1, NAN], "decision": [None, "unknown", "improves", "equivalent"],
    "interval_label": NOT_TEXT,
}
FEATURE_BAD = {
    "dataset": NOT_ID, "column": NOT_TEXT + ["bla`TEM", "a``b", "`"],
    "blocks": [None, [], "C", ["C", "C"], ["bad block"], [""], ["1K"], ["K+C"]],
    "selection_freq": NOT_PROB, "sign_stable": NOT_BOOL,
    "claim_level": [None, "mechanism", "", "mechanism_candidate"],
    "evidence_refs": [None, "doi:x", [""], ["a", "a"], [1]],
}


def errors_of(validate, rec):
    try:
        validate(rec)
    except SchemaError as e:
        return e.errors
    return []


class FieldCases(unittest.TestCase):
    """Good records validate; every listed bad value, a missing field and an
    unknown field are each refused with the field named in the error path."""

    def assert_fields(self, validate, good, bad):
        validate(good())
        for field, values in bad.items():
            for value in values + [MISSING]:
                rec = good()
                if value is MISSING:
                    del rec[field]
                else:
                    rec[field] = value
                with self.subTest(field=field, value=value):
                    errors = errors_of(validate, rec)
                    self.assertTrue(errors, "accepted")
                    self.assertTrue(any(f".{field}:" in e or f".{field}[" in e or f".{field}." in e
                                        for e in errors), errors)
        rec = good()
        rec["unexpected"] = 1
        with self.assertRaisesRegex(SchemaError, "unexpected: is not a field"):
            validate(rec)
        for not_record in (None, [], "x"):
            with self.assertRaises(SchemaError):
                validate(not_record)

    def test_prediction(self):
        self.assertEqual(set(PREDICTION_BAD), set(good_prediction()))
        self.assert_fields(schema.validate_prediction, good_prediction, PREDICTION_BAD)

    def test_contrast(self):
        self.assertEqual(set(CONTRAST_BAD), set(good_contrast()))
        self.assert_fields(schema.validate_contrast, good_contrast, CONTRAST_BAD)

    def test_pooled(self):
        self.assertEqual(set(POOLED_BAD), set(good_pooled()))
        self.assert_fields(schema.validate_pooled, good_pooled, POOLED_BAD)

    def test_feature(self):
        self.assertEqual(set(FEATURE_BAD), set(good_feature()))
        self.assert_fields(schema.validate_feature, good_feature, FEATURE_BAD)

    def test_contrast_has_exactly_the_contract_fields(self):
        self.assertEqual(set(good_contrast()), {
            "id", "family", "baseline", "treatment", "learner", "split", "endpoint", "estimate", "ci",
            "ci_tost", "ci_method", "unit", "n_units", "few_clusters", "p", "p_adj", "adjust", "margin",
            "decision", "claim_level", "undefined_resamples", "undefined_reason", "interval_label", "dataset",
            "pairing", "auroc_delta", "auroc_ci"})
        self.assertEqual(len(good_contrast()), 27)

    def test_pooled_has_exactly_the_contract_fields(self):
        self.assertEqual(set(good_pooled()), {
            "contrast", "family", "unit", "n_units", "n_datasets", "n_excluded", "excluded", "estimate", "ci",
            "ci_tost", "p", "p_adj", "adjust", "margin", "decision", "interval_label"})

    def test_contrast_id_parts(self):
        self.assertEqual(schema.contrast_dataset("ecoli__amp/C:lr:lineage:K->K+C:log_loss"), "ecoli__amp")
        self.assertEqual(schema.contrast_key("ecoli__amp/C:lr:lineage:K->K+C:log_loss"),
                         "C:lr:lineage:K->K+C:log_loss")
        self.assertEqual(schema.parse_key("B:lr:random:K+C->K:r2_ll"),
                         {"family": "B", "learner": "lr", "split": "random", "baseline": "K+C", "treatment": "K",
                          "endpoint": "r2_ll"})
        self.assertEqual(schema.parse_key("C:lr:lineage:a-b->c:log_loss")["baseline"], "a-b")
        self.assertEqual({k: schema.parse_key("C:lr:lineage:a-->b:log_loss")[k] for k in ("baseline", "treatment")},
                         {"baseline": "a-", "treatment": "b"})
        for bad in ("A-C", "C:lr:lineage:K->K+C", "E:lr:s:K->C:x", "C:lr:s:K`->C:x", "C:lr:s:K->C:x/y", None,
                    "C:lr:s:K->C->D:x", "C:lr:s:->C:x", "C:lr:s:K->:x", "C::s:K->C:x", "C:lr:s:K->C:"):
            self.assertIsNone(schema.parse_key(bad), bad)

    def test_pooled_cross_fields(self):
        cases = [
            ({"n_excluded": 0}, "must equal the number of excluded"),
            ({"n_datasets": 0}, "must count every excluded dataset"),
            ({"n_units": 5}, "species must come from the included"),
            ({"n_units": 0}, "species must come from the included"),
            ({"excluded": [{"dataset": "a", "reason": "r"}, {"dataset": "a", "reason": "s"}],
              "n_excluded": 2, "n_datasets": 3}, "repeats dataset"),
            ({"n_datasets": 1, "n_units": 0}, "estimate: must be null when there is no estimate"),
        ]
        for over, message in cases:
            with self.subTest(over=over), self.assertRaisesRegex(SchemaError, message):
                schema.validate_pooled(good_pooled(**over))
        schema.validate_pooled(all_excluded_pooled())
        with self.assertRaisesRegex(SchemaError, "decision: is improves but the contract rule gives inconclusive"):
            schema.validate_pooled(all_excluded_pooled(decision="improves"))


class NullShapeCases(unittest.TestCase):
    """H1: the contract's undefined-endpoint and family D shapes."""

    def test_undefined_endpoint(self):
        for family in schema.CONTRAST_FAMILIES:
            with self.subTest(family=family):
                schema.validate_contrast(good_contrast(family=family, **UNDEFINED))
        for field in ("estimate", "ci", "ci_tost", "p", "p_adj"):
            value = good_contrast()[field]
            with self.subTest(field=field), self.assertRaisesRegex(SchemaError, f"{field}: must be null"):
                schema.validate_contrast(good_contrast(**{**UNDEFINED, field: value}))
        with self.assertRaisesRegex(SchemaError, "decision: must be one of"):
            schema.validate_contrast(good_contrast(**{**UNDEFINED, "decision": "equivalent"}))

    def test_ci_method_is_an_exact_token(self):
        """Sections 10 and 12 (closes L3): equality against exactly cluster_jackknife_t, the input comparison's
        repeat_fold_corrected_t (main's rulings after VERIFY-H1STATS and the failed calibration) and none."""
        self.assertEqual(schema.CI_METHODS, ("cluster_jackknife_t", "repeat_fold_t", "repeat_fold_corrected_t", "none"))
        for method in schema.CI_METHODS:
            schema.validate_contrast(good_contrast(**{**UNDEFINED, "ci_method": method}))
        schema.validate_contrast(good_contrast(ci_method="cluster_jackknife_t"))
        schema.validate_contrast(good_contrast(ci_method="repeat_fold_corrected_t", unit="repeat_fold", few_clusters=False))
        with self.assertRaisesRegex(SchemaError, "few_clusters: must be n_units < 20, and false on unit repeat_fold"):
            schema.validate_contrast(good_contrast(ci_method="repeat_fold_corrected_t", unit="repeat_fold", few_clusters=True))
        for method in ("cluster_bootstrap_percentile_B2000", "cluster_bootstrap_percentile_B300", "jackknife_t"):
            with self.subTest(method=method), self.assertRaisesRegex(SchemaError, "ci_method: must be one of"):
                schema.validate_contrast(good_contrast(ci_method=method))
        with self.assertRaisesRegex(SchemaError, "a defined record names its interval method"):
            schema.validate_contrast(good_contrast(ci_method="none"))
        # The re-check's false positive: a reason that is a substring of the method is fine.
        schema.validate_contrast(good_contrast(**{**UNDEFINED, "undefined_reason": "cluster"}))
        # The re-check's paraphrase: free text in ci_method is refused.
        for method in ("no interval: vme cannot be computed, no resistant samples",
                       "no interval: " + UNDEFINED["undefined_reason"], "cluster bootstrap, percentile, 95%, B=300"):
            with self.subTest(method=method), self.assertRaisesRegex(SchemaError, "ci_method: must be one of"):
                schema.validate_contrast(good_contrast(**{**UNDEFINED, "ci_method": method}))

    def test_intervals_nest(self):
        """Section 10 (closes L2): ci[0] <= ci_tost[0] <= ci_tost[1] <= ci[1], on contrasts and pooled records."""
        schema.validate_contrast(good_contrast(ci=[0.015, 0.047]))
        for ci, tost in (([0.016, 0.05], [0.015, 0.047]), ([0.012, 0.046], [0.015, 0.047]),
                         ([0.012, 0.05], [-0.5, 0.5]), ([0.012, 0.05], [0.011, 0.051])):
            with self.subTest(ci=ci, tost=tost), self.assertRaisesRegex(SchemaError, "ci_tost: must nest inside ci"):
                schema.validate_contrast(good_contrast(ci=ci, ci_tost=tost))
        with self.assertRaisesRegex(SchemaError, "ci_tost: must nest inside ci"):
            schema.validate_pooled(good_pooled(ci_tost=[-0.01, 0.04]))
        # The re-check's forged interval: narrowing only ci_tost cannot escape ci.
        forged = good_contrast(**SHAPES["inconclusive"], decision="inconclusive")
        forged["ci_tost"] = [-0.008, 0.008]
        schema.validate_contrast(forged)
        forged["ci"] = [-0.007, 0.03]
        with self.assertRaisesRegex(SchemaError, "must nest"):
            schema.validate_contrast(forged)

    def test_defined_endpoint_needs_its_values(self):
        for field in ("estimate", "ci", "ci_tost"):
            with self.subTest(field=field), self.assertRaisesRegex(SchemaError, f"{field}: is required"):
                schema.validate_contrast(good_contrast(**{field: None}))

    def test_family_d_has_no_p_values(self):
        schema.validate_contrast(shaped("inside_margin", "D"))
        schema.validate_contrast(shaped("inconclusive_off_zero", "D"))
        for field in ("p", "p_adj"):
            with self.subTest(field=field), self.assertRaisesRegex(SchemaError, f"{field}: must be null in family D"):
                schema.validate_contrast(shaped("inside_margin", "D", **{field: 0.5}))
        with self.assertRaisesRegex(SchemaError, "adjust: family D uses none"):
            schema.validate_contrast(shaped("inside_margin", "D", adjust="holm"))
        for field in ("p", "p_adj"):
            with self.subTest(field=field), self.assertRaisesRegex(SchemaError, f"{field}: is required in family C"):
                schema.validate_contrast(shaped("inside_margin", "C", **{field: None}))

    def test_adjust_follows_family(self):
        for family in schema.CONTRAST_FAMILIES:
            adjust = schema.FAMILY_ADJUST[family]
            schema.validate_contrast(shaped("inside_margin", family))
            for other in set(schema.ADJUSTMENTS) - {adjust}:
                with self.subTest(family=family, adjust=other), self.assertRaisesRegex(SchemaError, "adjust"):
                    schema.validate_contrast(shaped("inside_margin", family, adjust=other))
        for family in schema.POOLED_FAMILIES:
            with self.subTest(family=family), self.assertRaisesRegex(SchemaError, "adjust"):
                schema.validate_pooled(good_pooled(contrast=key(family, endpoint="r2_ll"), family=family,
                                                   adjust="bh"))
        self.assertEqual(schema.FAMILY_ADJUST, {"A": "holm", "B": "holm", "C": "bh", "D": "none"})

    def test_p_adj_may_equal_p_but_never_fall_below(self):
        """Holm and BH both leave the largest p unchanged, so equality is allowed; no slack below p (I5)."""
        schema.validate_contrast(good_contrast(p=0.004, p_adj=0.004))
        for p_adj in (0.003, 0.004 - 1e-13, 0.0039999999999999):
            with self.subTest(p_adj=p_adj), self.assertRaisesRegex(SchemaError, "cannot be below p"):
                schema.validate_contrast(good_contrast(p=0.004, p_adj=p_adj))

    def test_dataset_repeats_the_id_prefix(self):
        with self.assertRaisesRegex(SchemaError, "dataset: must equal the id prefix ecoli__amp"):
            schema.validate_contrast(good_contrast(id="ecoli__amp/" + key(), dataset="ecoli__cip"))


class FamilyAndKeyCases(unittest.TestCase):
    """Sections 9 and 10: family membership, primary endpoints, and keys that agree with their records."""

    def test_contrasts_are_family_c_or_d(self):
        for family in ("A", "B"):
            with self.subTest(family=family), self.assertRaisesRegex(
                    SchemaError, "families A and B hold pooled records only"):
                schema.validate_contrast(good_contrast(family=family))
        schema.validate_contrast(shaped("inside_margin", "D"))

    def test_family_c_uses_the_primary_endpoint(self):
        for endpoint in ("auroc", "auprc", "balanced_accuracy", "vme", "me", "r2_ll"):
            with self.subTest(endpoint=endpoint), self.assertRaisesRegex(
                    SchemaError, "family C uses the primary endpoint log_loss only"):
                schema.validate_contrast(good_contrast(endpoint=endpoint))
            schema.validate_contrast(shaped("inside_margin", "D", endpoint=endpoint))

    def test_key_agrees_with_the_record(self):
        for field, value in (("baseline", "K+C"), ("treatment", "K"), ("learner", "rf"), ("split", "random"),
                             ("endpoint", "auroc"), ("family", "D")):
            rec = good_contrast()
            rec[field] = value
            with self.subTest(field=field), self.assertRaisesRegex(
                    SchemaError, rf"contrast.{field}: disagrees with the id's key"):
                schema.validate_contrast(rec)

    def test_backticks_are_refused(self):
        """Section 10 (closes L1): no backticks in contrast keys or feature columns."""
        for part in ("learner", "split", "baseline", "treatment", "endpoint"):
            rec = good_contrast(**{part: "a``b"})
            with self.subTest(part=part), self.assertRaisesRegex(SchemaError, "contrast.id: must match"):
                schema.validate_contrast(rec)
        with self.assertRaisesRegex(SchemaError, "pooled.contrast: must match"):
            schema.validate_pooled(good_pooled(contrast="A:lasso:lineage:K->K``[x](y)``:r2_ll"))
        for column in ("a``b", "`", "x`[click](x.invalid)`"):
            with self.subTest(column=column), self.assertRaisesRegex(SchemaError, "feature.column: must match"):
                schema.validate_feature(good_feature(column=column))

    def test_pooled_family_follows_its_key(self):
        """I3: the family its key implies, even for an all-excluded pool."""
        for family in schema.POOLED_FAMILIES:
            endpoint = "auroc" if family == "D" else "r2_ll"
            values = {"p": None, "p_adj": None} if family == "D" else {}
            schema.validate_pooled(good_pooled(contrast=key(family, endpoint=endpoint), family=family,
                                               adjust=schema.FAMILY_ADJUST[family], **values))
        with self.assertRaisesRegex(SchemaError, "pooled.family: must equal the family its key implies, A"):
            schema.validate_pooled(all_excluded_pooled(family="D", adjust="none"))
        with self.assertRaisesRegex(SchemaError, "family C holds per-dataset contrasts only"):
            schema.validate_pooled(good_pooled(contrast=key("C", endpoint="r2_ll"), family="C", adjust="bh"))
        for family in ("A", "B"):
            with self.subTest(family=family), self.assertRaisesRegex(
                    SchemaError, f"family {family} pools the primary endpoint r2_ll only"):
                schema.validate_pooled(good_pooled(contrast=key(family, endpoint="log_loss"), family=family))

    def test_no_adjusted_family_is_shared(self):
        """Section 10: a pooled record never shares an adjusted family with per-dataset records."""
        shared = (set(schema.CONTRAST_FAMILIES) & set(schema.POOLED_FAMILIES))
        self.assertEqual(shared, {"D"})
        self.assertEqual(schema.FAMILY_ADJUST["D"], "none")


class EquivalenceWithdrawnCases(unittest.TestCase):
    """Section 15: no record is ever equivalent; the decision set is improves, worse, inconclusive."""

    def test_decision_set(self):
        self.assertEqual(schema.DECISIONS, ("improves", "worse", "inconclusive"))
        self.assertFalse(hasattr(schema, "MIN_EQUIVALENCE_UNITS"))
        self.assertFalse(hasattr(schema, "equivalence_assessable"))

    def test_equivalent_is_refused_everywhere(self):
        for family in schema.CONTRAST_FAMILIES:
            for shape in SHAPES:
                for n in (1, 59, 60, 64, 1000):
                    rec = shaped(shape, family, n_units=n, few_clusters=n < 20, decision="equivalent")
                    with self.subTest(family=family, shape=shape, n=n), self.assertRaisesRegex(
                            SchemaError, "decision: must be one of"):
                        schema.validate_contrast(rec)
        with self.assertRaisesRegex(SchemaError, "decision: must be one of"):
            schema.validate_contrast(good_contrast(**{**UNDEFINED, "decision": "equivalent"}))
        values = {**SHAPES["inside_margin"], "n_datasets": 400, "n_excluded": 0, "excluded": []}
        for n in (1, 10, 60, 320):
            with self.subTest(n=n), self.assertRaisesRegex(SchemaError, "decision: must be one of"):
                schema.validate_pooled(good_pooled(**values, n_units=n, decision="equivalent"))

    def test_inside_the_margin_is_inconclusive(self):
        """The former equivalent shape is inconclusive at every unit count, contrast or pooled."""
        for family in schema.CONTRAST_FAMILIES:
            for n in (10, 59, 60, 61, 1000):
                for unit in ("cluster", "lineage"):
                    rec = shaped("inside_margin", family, n_units=n, few_clusters=n < 20, unit=unit)
                    self.assertEqual(rec["decision"], "inconclusive")
                    with self.subTest(family=family, n=n, unit=unit):
                        schema.validate_contrast(rec)
        values = {**SHAPES["inside_margin"], "n_datasets": 400, "n_excluded": 0, "excluded": []}
        for n in (2, 10, 60, 320):
            with self.subTest(n=n):
                schema.validate_pooled(good_pooled(**values, n_units=n, decision="inconclusive"))

    def test_constants_match_j(self):
        """J's stats.py is the source; K mirrors its values. Section 15 removed MIN_EQUIVALENCE_UNITS.
        J's file is the tree's own unless GRM_TASK_J_STATS names another copy (R's L6: a plain run checks)."""
        path = Path(os.environ.get("GRM_TASK_J_STATS") or ROOT / "scripts" / "ablation" / "stats.py")
        spec = importlib.util.spec_from_file_location("grm_task_j_stats", path)
        st = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = st  # J's dataclasses look their module up while it loads
        self.addCleanup(sys.modules.pop, spec.name, None)
        spec.loader.exec_module(st)
        self.assertFalse(hasattr(st, "MIN_EQUIVALENCE_UNITS"))
        for ps in ([0.01, 0.02, 0.03], [0.03, 0.01, 0.02], [0.2, 0.2, 1.0], [0.004, 0.03, 1.0, 0.3, 0.5, 1.0], [0.7]):
            with self.subTest(ps=ps):
                self.assertEqual(schema.holm(ps), st.holm(ps))
                self.assertEqual(schema.bh(ps), st.bh(ps))
        self.assertEqual(schema.DECISIONS, st.DECISIONS)
        self.assertEqual((schema.FEW_CLUSTERS, schema.MAX_UNDEFINED_SHARE, schema.UNDEFINED_SHARE_REASON),
                         (st.FEW_CLUSTERS, st.MAX_UNDEFINED_SHARE, st.UNDEFINED_SHARE_REASON))
        self.assertEqual(set(schema.CI_METHODS), set(st.CI_METHODS))
        self.assertEqual(schema.FAMILY_ADJUST, st.FAMILIES)
        self.assertEqual(set(schema.CONTRAST_FAMILIES), set(st.CONTRAST_FAMILIES))
        self.assertEqual(set(schema.POOLED_FAMILIES), set(st.POOLED_FAMILIES))
        self.assertEqual(schema.PRIMARY_CONTRAST_ENDPOINT, st.PRIMARY)
        self.assertEqual(schema.PRIMARY_POOLED_ENDPOINT, st.POOLED_ENDPOINT)
        self.assertEqual(schema.ALPHA, st.ALPHA)
        self.assertEqual(set(good_contrast()), set(st.RECORD_KEYS))
        self.assertEqual(schema.PAIRINGS, st.PAIRINGS)
        self.assertEqual(set(good_pooled()), set(st.POOLED_KEYS))
        self.assertEqual(set(good_prediction()), set(st.FIELDS))
        key = st.contrast_key("C", "lr", "lineage", "base", "base+C", "log_loss")
        self.assertEqual(schema.parse_key(key), {"family": "C", "learner": "lr", "split": "lineage",
                                                 "baseline": "base", "treatment": "base+C", "endpoint": "log_loss"})
        for ci in (None, [-0.01, 0.01], [0.0, 0.05], [0.001, 0.05], [-0.05, 0.0], [-0.05, -0.001]):
            for p_adj in (None, 0.0, 0.049, 0.05, 0.5, 1.0):
                with self.subTest(ci=ci, p_adj=p_adj):
                    self.assertEqual(schema.expected_decision(ci, p_adj), st.decide(ci, p_adj))


class FewSpeciesCases(unittest.TestCase):
    """Section 17: a defined pooled record needs 2 contributing species and p >= 2/2^S - 1e-12."""

    def pooled(self, n, **over):
        return good_pooled(**{"n_units": n, "n_datasets": max(n, 1) + 1, "n_excluded": 1, **over})

    def test_constants(self):
        self.assertEqual((schema.MIN_POOLED_SPECIES, schema.FLOOR_TOLERANCE, schema.ADJUST_TOLERANCE),
                         (2, 1e-12, 1e-12))
        self.assertEqual([schema.p_floor(s) for s in (2, 3, 5, 6)], [0.5, 0.25, 0.0625, 0.03125])

    def test_fewer_than_two_species_is_refused_when_defined(self):
        for family, extra in (("A", {}), ("B", {"contrast": key("B", "K+C", "K", "r2_ll")}),
                              ("D", {"contrast": key("D", endpoint="r2_ll"), "adjust": "none", "p": None,
                                     "p_adj": None})):
            rec = self.pooled(1, family=family, p=1.0, p_adj=1.0, **extra) if family != "D" else \
                self.pooled(1, family=family, **extra)
            with self.subTest(family=family), self.assertRaisesRegex(
                    SchemaError, "n_units: a defined pooled record needs at least 2 contributing species"):
                schema.validate_pooled(rec)
        schema.validate_pooled(all_excluded_pooled())  # undefined, with its exclusions listed
        schema.validate_pooled(self.pooled(2, p=0.5, p_adj=0.5))

    def test_p_floor_in_family_b(self):
        """TK-L2: the floor binds family B as well as A."""
        b = {"contrast": key("B", "K+C", "K", "r2_ll"), "family": "B"}
        schema.validate_pooled(self.pooled(3, p=0.25, p_adj=0.5, **b))
        with self.assertRaisesRegex(SchemaError, r"pooled.p: is below the exact sign-flip floor 2/2\^S = 0.25 for S = 3"):
            schema.validate_pooled(self.pooled(3, estimate=0.03, ci=[0.01, 0.05], ci_tost=[0.015, 0.045], p=0.01,
                                               p_adj=0.02, decision="improves", **b))

    def test_p_floor(self):
        for n, ok, low in ((2, 0.5, 0.49), (3, 0.25, 0.2), (5, 0.0625, 0.01), (6, 0.03125, 0.03)):
            schema.validate_pooled(self.pooled(n, p=ok, p_adj=min(1.0, 2 * ok)))
            with self.subTest(n=n), self.assertRaisesRegex(SchemaError, rf"pooled.p: is below the exact sign-flip "
                                                                        rf"floor 2/2\^S = .* for S = {n} species"):
                schema.validate_pooled(self.pooled(n, p=low, p_adj=min(1.0, 2 * ok)))
        schema.validate_pooled(self.pooled(3, p=0.25 - 1e-13, p_adj=0.5))
        with self.assertRaisesRegex(SchemaError, "below the exact sign-flip floor"):
            schema.validate_pooled(self.pooled(3, p=0.25 - 1e-11, p_adj=0.5))

    def test_pooled_improves_or_worse_below_six_species_is_refused(self):
        for n in (2, 3, 4, 5):
            for shape in ("improves", "worse"):
                values = {**SHAPES[shape], "decision": shape}
                with self.subTest(n=n, shape=shape), self.assertRaisesRegex(SchemaError, "below the exact sign-flip"):
                    schema.validate_pooled(self.pooled(n, **values))
        floor6 = schema.p_floor(6)
        schema.validate_pooled(self.pooled(6, estimate=0.03, ci=[0.01, 0.05], ci_tost=[0.015, 0.045], p=floor6,
                                           p_adj=floor6, decision="improves"))


class AdjustmentCases(unittest.TestCase):
    """Section 17: K's Holm and BH (J's arithmetic in plain Python)."""

    def test_holm_and_bh_differ(self):
        self.assertEqual(schema.holm([0.01, 0.02, 0.03]), [0.03, 0.04, 0.04])
        self.assertEqual(schema.bh([0.01, 0.02, 0.03]), [0.03, 0.03, 0.03])

    def test_input_order_ties_and_cap(self):
        self.assertEqual(schema.holm([0.03, 0.01, 0.02]), [0.04, 0.03, 0.04])
        self.assertEqual(schema.bh([0.03, 0.01, 0.02]), [0.03, 0.03, 0.03])
        self.assertEqual(schema.holm([0.2, 0.2, 1.0]), [0.6000000000000001, 0.6000000000000001, 1.0])
        self.assertEqual(schema.bh([0.2, 0.2, 1.0]), [0.30000000000000004, 0.30000000000000004, 1.0])
        self.assertEqual(schema.holm([0.6, 0.7]), [1.0, 1.0])
        self.assertEqual(schema.bh([0.004, 0.03, 1.0, 0.30, 0.5, 1.0]), [0.024, 0.09, 1.0, 0.6, 0.75, 1.0])
        self.assertEqual((schema.holm([]), schema.bh([]), schema.holm([0.3]), schema.bh([0.3])), ([], [], [0.3], [0.3]))
        self.assertEqual(schema.ADJUSTERS, {"holm": schema.holm, "bh": schema.bh})


class JackknifeCases(unittest.TestCase):
    """Section 13: undefined replicates under the cluster jackknife, and few_clusters."""

    def undefined(self, reason, u, k):
        return good_contrast(**{**UNDEFINED, "undefined_reason": reason, "undefined_resamples": u, "n_units": k,
                                "few_clusters": k < 20, "ci_method": "none"})

    def test_js_three_reasons_are_accepted(self):
        schema.validate_contrast(self.undefined(schema.ONE_CLUSTER_REASON, 1, 1))
        schema.validate_contrast(self.undefined("more than 10% of resamples undefined", 7, 64))
        schema.validate_contrast(self.undefined("6 of 64 leave-one-cluster-out estimates undefined; the jackknife "
                                                "needs all", 6, 64))
        schema.validate_contrast(self.undefined("1 of 10 leave-one-cluster-out estimates undefined; the jackknife "
                                                "needs all", 1, 10))
        schema.validate_contrast(self.undefined("vme undefined: no resistant (y=1) samples", 64, 64))

    def test_reasons_agree_with_the_counts(self):
        cases = [
            (schema.ONE_CLUSTER_REASON, 1, 2, "says one cluster"),
            ("more than 10% of resamples undefined", 6, 64, "says more than 10%"),
            ("more than 10% of resamples undefined", 1, 10, "says more than 10%"),
            ("5 of 64 leave-one-cluster-out estimates undefined; the jackknife needs all", 6, 64,
             "must name undefined_resamples of n_units"),
            ("6 of 60 leave-one-cluster-out estimates undefined; the jackknife needs all", 6, 64,
             "must name undefined_resamples of n_units"),
            ("7 of 64 leave-one-cluster-out estimates undefined; the jackknife needs all", 7, 64,
             "is for 1 to 10% undefined replicates"),
            ("0 of 64 leave-one-cluster-out estimates undefined; the jackknife needs all", 0, 64,
             "is for 1 to 10% undefined replicates"),
        ]
        for reason, u, k, message in cases:
            with self.subTest(reason=reason, u=u, k=k), self.assertRaisesRegex(SchemaError, message):
                schema.validate_contrast(self.undefined(reason, u, k))

    def test_defined_record_has_no_undefined_replicate(self):
        with self.assertRaisesRegex(SchemaError, "must be 0 on a defined record"):
            schema.validate_contrast(good_contrast(undefined_resamples=1))

    def test_few_clusters_is_below_20(self):
        for n, few in ((19, True), (20, False), (1, True), (64, False)):
            schema.validate_contrast(good_contrast(n_units=n, few_clusters=few, decision="improves"))
            with self.subTest(n=n), self.assertRaisesRegex(SchemaError, "few_clusters: must be n_units < 20"):
                schema.validate_contrast(good_contrast(n_units=n, few_clusters=not few))
        self.assertEqual((schema.FEW_CLUSTERS, schema.MIN_TESTABLE_SPECIES), (20, 6))


class DecisionRuleCases(unittest.TestCase):
    """M4: every decision is re-derived with the contract rule."""

    def test_rule(self):
        rule = schema.expected_decision
        self.assertEqual(rule([0.01, 0.05], 0.049), "improves")
        self.assertEqual(rule([0.01, 0.05], 0.05), "inconclusive")
        self.assertEqual(rule([0.0, 0.05], 0.001), "inconclusive")
        self.assertEqual(rule([-0.05, -0.01], 0.001), "worse")
        self.assertEqual(rule([-0.05, 0.0], 0.001), "inconclusive")
        self.assertEqual(rule([-0.05, 0.05], 0.5), "inconclusive")
        self.assertEqual(rule([-0.05, 0.05], 0.001), "inconclusive")
        self.assertEqual(rule([0.001, 0.005], 0.001), "improves")
        # Family D has no p: always inconclusive, even off zero (section 15).
        self.assertEqual(rule([0.001, 0.005], None), "inconclusive")
        self.assertEqual(rule([-0.05, -0.03], None), "inconclusive")
        self.assertEqual(rule(None, None), "inconclusive")
        self.assertEqual(rule(None, 0.001), "inconclusive")

    def test_every_wrong_decision_is_refused(self):
        for family in schema.CONTRAST_FAMILIES:
            for shape in SHAPES:
                good = shaped(shape, family)
                schema.validate_contrast(good)
                for wrong in set(schema.DECISIONS) - {good["decision"]}:
                    with self.subTest(family=family, shape=shape, wrong=wrong):
                        with self.assertRaisesRegex(SchemaError, "the contract rule gives"):
                            schema.validate_contrast({**good, "decision": wrong})

    def test_reviewer_probes(self):
        with self.assertRaisesRegex(SchemaError, "is improves but the contract rule gives inconclusive"):
            schema.validate_contrast(good_contrast(ci=[-0.010, 0.026], ci_tost=[-0.007, 0.023], p=0.4, p_adj=0.8))
        with self.assertRaisesRegex(SchemaError, "decision: must be one of"):
            schema.validate_contrast(good_contrast(decision="equivalent", ci=[-0.5, 0.5], ci_tost=[-0.4, 0.4],
                                                   p=0.9, p_adj=0.9))

    def test_family_d_decides_no_direction(self):
        for shape in ("improves", "worse"):
            rec = shaped(shape, "D")
            self.assertEqual(rec["decision"], "inconclusive")
            schema.validate_contrast(rec)
            with self.subTest(shape=shape), self.assertRaisesRegex(SchemaError, "the contract rule gives"):
                schema.validate_contrast({**rec, "decision": shape})

    def test_pooled_rule(self):
        for shape in SHAPES:
            values = {**SHAPES[shape], "n_units": 12, "n_datasets": 12, "n_excluded": 0, "excluded": []}
            decision = decision_of(shape, "A")
            schema.validate_pooled(good_pooled(**values, decision=decision))
            with self.subTest(shape=shape), self.assertRaisesRegex(SchemaError, "the contract rule gives"):
                schema.validate_pooled(good_pooled(**values, decision="improves" if decision != "improves"
                                                   else "worse"))


class StrictValueCases(unittest.TestCase):
    """L4: huge numbers and NaN raise SchemaError, never another exception."""

    def test_huge_and_non_finite(self):
        for value in (HUGE, -HUGE, 2 ** 54, NAN, INF):
            with self.subTest(value=value), self.assertRaises(SchemaError):
                schema.validate_contrast(good_contrast(estimate=value))
            with self.subTest(value=value), self.assertRaises(SchemaError):
                schema.validate_contrast(good_contrast(ci=[0, value]))
        with self.assertRaises(SchemaError):
            schema.validate_contrast(good_contrast(n_units=HUGE))

    def test_canonical_json(self):
        self.assertEqual(schema.canonical({"b": [1, 2.5], "a": "é"}), '{"a":"\\u00e9","b":[1,2.5]}')
        for bad, message in (({1: "x"}, "key 1 is not a string"), ({"a": NAN}, "not canonical JSON"),
                             ({"a": {2: 1}}, "a: key 2")):
            with self.subTest(bad=bad), self.assertRaisesRegex(SchemaError, message):
                schema.canonical(bad)
        self.assertEqual(len(schema.sha256({"a": 1})), 64)


class ManifestCases(unittest.TestCase):
    def setUp(self):
        self.m = good_manifest()

    def refuse(self, message):
        with self.assertRaisesRegex(SchemaError, message):
            schema.validate_manifest(self.m)

    def test_good(self):
        self.assertIs(schema.validate_manifest(self.m), self.m)

    def test_top_level_fields(self):
        bad = {
            "schema_version": [None, 2, "1", True], "protocol_sha256": NOT_HEX,
            "source_commit": [None, "c90ebaf", "C" * 40, "1" * 41], "binary_sha256": NOT_HEX,
            "blocks": [None, {}], "datasets": [None, []], "missingness": [None, {}],
            "group_proxy_audit": [None, "x"], "runs": [None, []], "compute": [None, []],
        }
        self.assertEqual(set(bad), set(self.m))
        for field, values in bad.items():
            for value in values + [MISSING]:
                m = good_manifest()
                if value is MISSING:
                    del m[field]
                else:
                    m[field] = value
                with self.subTest(field=field, value=value):
                    with self.assertRaisesRegex(SchemaError, f"manifest.{field}"):
                        schema.validate_manifest(m)
        schema.validate_manifest({**good_manifest(), "source_commit": "f" * 64})

    def assert_nested(self, section, index, bad):
        for field, values in bad.items():
            for value in values + [MISSING]:
                m = good_manifest()
                rec = m[section][index] if index is not None else m[section]
                if value is MISSING:
                    del rec[field]
                else:
                    rec[field] = value
                with self.subTest(section=section, field=field, value=value):
                    with self.assertRaisesRegex(SchemaError, rf"{field}"):
                        schema.validate_manifest(m)
        m = good_manifest()
        (m[section][index] if index is not None else m[section])["extra"] = 1
        with self.assertRaisesRegex(SchemaError, "extra: is not a field"):
            schema.validate_manifest(m)

    def test_block_fields(self):
        self.assert_nested("blocks", 1, {
            "id": NOT_BLOCK_ID, "kind": NOT_TEXT + ["curated", "anonymous", "SyntheticCount"],
            "version": [None, "", "  ", 3],
            "source_sha256": [None, "", "A" * 64, "1" * 63, 3],
            "priority": [None, -1, 6, True, 1.5, 3], "screen_cap": [0, -1, True, "4000", 1.5],
            "catalog_date": ["2026-13-01", "2026-9-15", "15/09/2026", "2026-09-15T00:00", 20260915],
            "circularity_risk": [None, "high", "none", 0, 1], "unavailable_reason": ["", 1],
        })

    def test_block_is_h_shaped(self):
        """M5: H's serde JSON validates as is, for available and unavailable blocks."""
        h_available = {"catalog_date": None, "circularity_risk": False, "id": "K", "kind": "Anonymous",
                       "priority": 5, "screen_cap": 4000, "source_sha256": "0123456789abcdef" * 4,
                       "unavailable_reason": None, "version": "k=31 min_support=1 patterns=3 genomes=5"}
        h_unavailable = {"id": "C", "kind": "CuratedExact", "version": "", "source_sha256": "", "priority": 0,
                         "screen_cap": None, "catalog_date": None, "circularity_risk": True,
                         "unavailable_reason": "no source"}
        for block in (h_available, h_unavailable, {**h_unavailable, "source_sha256": H1, "version": "v1"}):
            self.m["blocks"] = [block]
            self.m["missingness"], self.m["group_proxy_audit"] = [], []
            self.m["runs"] = [good_run(blocks=[], arm="none")]
            schema.validate_manifest(self.m)
        self.m["blocks"] = [{**h_unavailable, "source_sha256": "abc"}]
        self.refuse("must be empty or 64 lowercase hex")
        self.m["blocks"] = [{**h_available, "version": ""}]
        self.refuse("must be non-empty for an available block")
        self.m["blocks"] = [{**h_available, "source_sha256": ""}]
        self.refuse("must be 64 lowercase hex characters for an available block")
        self.m["blocks"] = [{**h_available, "priority": 4}]
        self.refuse("kind Anonymous has priority 5")
        self.m["blocks"] = [{**h_available, "kind": "SyntheticCount", "priority": 3}]
        self.refuse("kind: must be one of CuratedExact, Hierarchical, Presence, Count, Context, Anonymous")
        for kind, priority in schema.H_KIND_PRIORITY.items():
            self.m["blocks"] = [{**h_available, "kind": kind, "priority": priority}]
            schema.validate_manifest(self.m)

    def test_block_ids(self):
        """Contract v1.1 section 8: block ids match [A-Za-z][A-Za-z0-9_-]*, and L's noise block validates."""
        self.assertEqual(schema.BLOCK_ID.pattern, "[A-Za-z][A-Za-z0-9_-]*")
        noise = good_block(id="noise", kind="Anonymous", priority=5, screen_cap=None, version="synthetic-v1",
                           source_sha256=schema.sha256([{"name": "noise.0"}]))
        self.m["blocks"].append(noise)
        self.m["runs"].append(good_run(id="r-noise", arm="K+noise", blocks=["K", "noise"]))
        self.m["missingness"].append({"dataset": "ecoli__amp", "block": "noise", "rate_R": 0.0, "rate_S": 0.0,
                                      "rate_diff_p": None, "indicator_auroc": None,
                                      "indicator_auroc_oriented": None, "flag": False})
        schema.validate_manifest(self.m)
        schema.validate_feature(good_feature(blocks=["noise", "C_x-1"], claim_level="predictive_increment"))
        for bad in NOT_BLOCK_ID[5:]:
            with self.subTest(bad=bad), self.assertRaisesRegex(SchemaError, r"blocks\[0\]: must match"):
                schema.validate_feature(good_feature(blocks=[bad]))

    def test_dataset_fields(self):
        self.assert_nested("datasets", 0, {
            "id": NOT_ID, "species": NOT_TEXT, "partition": [None, "test", "Confirmatory"],
            "n_samples": NOT_INT0, "group_counts": [None, [], {"lineage": -1}, {"lineage": 1.5}, {"": 3}],
        })

    def test_missingness_fields(self):
        self.assert_nested("missingness", 0, {
            "dataset": NOT_ID, "block": NOT_BLOCK_ID, "rate_R": NOT_PROB, "rate_S": NOT_PROB,
            "rate_diff_p": [NAN, INF, -0.1, 1.5, True, "0.01", HUGE],
            "indicator_auroc": [NAN, 1.2, True, None, 0.5],
            "indicator_auroc_oriented": [NAN, 0.4, 1.2, True, None, 0.7], "flag": NOT_BOOL + [False],
        })

    def test_missingness_flag_rule(self):
        """Contract v1.1 section 8: flag == (oriented > 0.60) or (rate_diff_p not null and < 0.01)."""
        row = self.m["missingness"][0]
        cases = [  # raw, oriented, rate_diff_p, flag, valid
            (0.36, 0.64, 0.02, True, True), (0.64, 0.64, None, True, True), (0.36, 0.64, 0.02, False, False),
            (0.45, 0.55, 0.02, False, True), (0.45, 0.55, 0.02, True, False), (0.45, 0.55, 0.009, True, True),
            (0.45, 0.55, 0.009, False, False), (0.45, 0.55, 0.01, False, True), (0.45, 0.55, 0.01, True, False),
            (0.6, 0.6, None, False, True), (0.6, 0.6, None, True, False), (None, None, None, False, True),
            (None, None, None, True, False), (None, None, 0.001, True, True), (None, None, 0.001, False, False),
            (0.4, 0.6, 0.5, False, True), (0.39, 0.61, 0.5, True, True), (0.39, 0.61, 0.5, False, False),
        ]
        for raw, oriented, p, flag, valid in cases:
            row.update(indicator_auroc=raw, indicator_auroc_oriented=oriented, rate_diff_p=p, flag=flag)
            with self.subTest(raw=raw, p=p, flag=flag):
                if valid:
                    schema.validate_manifest(self.m)
                else:
                    self.refuse(rf"missingness\[0\].flag: must be {str(not flag).lower()}")
        row.update(indicator_auroc=0.3, indicator_auroc_oriented=0.6, rate_diff_p=None, flag=False)
        self.refuse("must equal max")

    def test_group_proxy_fields(self):
        self.assert_nested("group_proxy_audit", 0, {
            "dataset": NOT_ID, "block": NOT_BLOCK_ID, "metric": NOT_TEXT, "value": [NAN, True, "0.9"],
            "flag": NOT_BOOL,
        })

    def test_run_fields(self):
        self.assert_nested("runs", 0, {
            "id": NOT_ID, "dataset": NOT_ID, "arm": NOT_TEXT, "blocks": [None, "K", ["K", "K"], ["bad id"], ["K.v2"]],
            "learner": NOT_TEXT, "split": NOT_TEXT, "repeats": NOT_INT0 + [0],
            "fold_plan_sha256": [None, [], ["x", "x"], "a" * 64, [LINEAGE_PLANS[0]]],
            "status": [None, "failed", "OK"], "reason": ["", 3],
            "metrics": [None, [], {"auroc": NAN}, {"auroc": "0.8"}], "fit_count": NOT_INT0,
            "seconds": NOT_NUM + [-1], "predictions": ["", 3],
        })

    def test_compute_fields(self):
        self.assert_nested("compute", None, {
            "planned_fit_units": NOT_INT0, "actual_seconds": NOT_NUM + [-0.1],
            "tiers_completed": [None, "smoke", ["smoke", "smoke"], [""]],
        })

    def test_run_status_rules(self):
        ok, unavailable = self.m["runs"][0], self.m["runs"][-1]
        cases = [
            (ok, {"reason": "x"}, "must be null for an ok run"),
            (ok, {"predictions": None}, "predictions: is required for an ok run"),
            (ok, {"fold_plan_sha256": None}, "fold_plan_sha256: is required for an ok run"),
            (ok, {"fit_count": 0}, "at least one fit"),
            (unavailable, {"reason": None}, "must give its reason"),
        ]
        for run, over, message in cases:
            saved = dict(run)
            run.update(over)
            with self.subTest(over=over):
                self.refuse(message)
            run.clear()
            run.update(saved)

    def test_ok_run_cannot_use_unavailable_block(self):
        self.m["runs"][-1].update(status="ok", reason=None, predictions="p.jsonl",
                                  fold_plan_sha256=list(LINEAGE_PLANS), fit_count=1)
        self.refuse("ok run uses unavailable block P")

    def test_cross_references(self):
        cases = [
            (lambda m: m["runs"][0].update(dataset="nope"), "unknown dataset nope"),
            (lambda m: m["runs"][0].update(blocks=["K", "nope"]), "unknown block nope"),
            (lambda m: m["missingness"][0].update(block="nope"), "unknown block nope"),
            (lambda m: m["group_proxy_audit"][0].update(dataset="nope"), "unknown dataset nope"),
            (lambda m: m["runs"][1].update(id="r1"), "repeats run id"),
            (lambda m: m["runs"][1].update(arm="K", blocks=["K"]), r"repeats \(dataset, arm, learner, split\)"),
            (lambda m: m["datasets"].append(good_dataset()), "repeats dataset id"),
            (lambda m: m["blocks"].append(good_block()), "repeats block id"),
            (lambda m: m["missingness"].append(dict(m["missingness"][0])), r"repeats \(dataset, block\)"),
            (lambda m: m["group_proxy_audit"].append(dict(m["group_proxy_audit"][0])), "repeats"),
            (lambda m: m["runs"][1].update(blocks=["C"]), "arm K\\+C has different blocks"),
            (lambda m: m["runs"][1].update(fold_plan_sha256=["e" * 64, LINEAGE_PLANS[1]]),
             r"fold_plan_sha256\[0\]: differs from another run's plan"),
        ]
        for mutate, message in cases:
            m = good_manifest()
            mutate(m)
            with self.subTest(message=message), self.assertRaisesRegex(SchemaError, message):
                schema.validate_manifest(m)

    def test_errors_are_all_reported(self):
        self.m["protocol_sha256"] = "x"
        self.m["runs"][0]["seconds"] = -1
        with self.assertRaises(SchemaError) as caught:
            schema.validate_manifest(self.m)
        self.assertEqual(len(caught.exception.errors), 2, caught.exception.errors)


class PredictionSetCases(unittest.TestCase):
    """M6: J's row fields, score in [0, 1], one fold plan per (dataset, split, repeat)."""

    def rows(self):
        out = []
        for arm in ("K", "K+C"):
            for repeat in range(2):
                for i, y in enumerate((0, 1, 1)):
                    out.append(good_prediction(arm=arm, repeat=repeat, sample_id=f"S{i}", y=y, fold=(i + repeat) % 3,
                                               score=0.8 if y else 0.2, called=bool(y), cluster=f"L{i}",
                                               fold_plan_sha256=LINEAGE_PLANS[repeat]))
        return out

    def test_good(self):
        rows = self.rows()
        self.assertEqual(schema.validate_predictions(rows, good_manifest()), rows)
        self.assertEqual(schema.validate_predictions(iter(rows)), rows)

    def test_called(self):
        schema.validate_prediction(good_prediction(called=False, score=0.83))  # a no-call row
        schema.validate_prediction(good_prediction(called=True, score=0.5))
        with self.assertRaisesRegex(SchemaError, "a called row needs score >= 0.5"):
            schema.validate_prediction(good_prediction(called=True, score=0.49))
        schema.validate_prediction(good_prediction(cluster=7))

    def test_refusals(self):
        cases = [
            (lambda r: r.append(dict(r[0])), "repeats key"),
            (lambda r: r[3].update(y=1 - r[3]["y"]), "label"),
            (lambda r: r[3].update(cluster="L9"), "cluster"),
            (lambda r: r[6].update(fold=(r[6]["fold"] + 1) % 3), "fold"),
            (lambda r: r[6].update(fold_plan_sha256="f" * 64), "fold plan"),
            (lambda r: r[0].update(score=1.5), "score: must be in"),
        ]
        for mutate, message in cases:
            rows = self.rows()
            mutate(rows)
            with self.subTest(message=message), self.assertRaisesRegex(SchemaError, message):
                schema.validate_predictions(rows)

    def test_manifest_links(self):
        m = good_manifest()
        cases = [
            ({"arm": "K+P"}, "belongs to no ok run"),
            ({"dataset": "kpneu__amp", "arm": "K+C"}, "belongs to no ok run"),
            ({"repeat": 2, "fold_plan_sha256": LINEAGE_PLANS[0]}, "has 2 repeats"),
            ({"repeat": 1}, "differs from run r1's plan for this repeat"),
        ]
        for over, message in cases:
            with self.subTest(over=over), self.assertRaisesRegex(SchemaError, message):
                schema.validate_predictions([good_prediction(**over)], m)


class ClaimLevelCases(unittest.TestCase):
    def test_enum_order(self):
        self.assertEqual(schema.CLAIM_LEVELS, ("predictive_increment", "association_stable",
                                               "mechanism_candidate", "causal_supported"))

    def test_contrasts_always_predictive_increment(self):
        schema.validate_contrast(good_contrast())
        for level in ("association_stable", "mechanism_candidate", "causal_supported"):
            with self.subTest(level=level), self.assertRaisesRegex(SchemaError, "contrasts always carry"):
                schema.validate_contrast(good_contrast(claim_level=level))
        for extra in ({"evidence_refs": ["x"]}, {"annotations": [note()]}):
            with self.subTest(extra=extra), self.assertRaisesRegex(SchemaError, "is not a field"):
                schema.validate_contrast(good_contrast(**extra))
        with self.assertRaisesRegex(SchemaError, "only features take manual levels"):
            schema.annotate_claim(good_contrast(), "contrast", "causal_supported", ["x"], "A", "2026-10-30")

    def test_pipeline_feature_allows_at_most_association_stable(self):
        schema.validate_feature(good_feature(claim_level="predictive_increment"))
        schema.validate_feature(good_feature(claim_level="association_stable"))
        for level in ("mechanism_candidate", "causal_supported"):
            with self.subTest(level=level), self.assertRaisesRegex(SchemaError, "pipeline output may set at most"):
                schema.validate_feature(good_feature(claim_level=level, evidence_refs=["doi:10.1/x"]))

    def test_stable_sign_needed_from_association_stable_up(self):
        with self.assertRaisesRegex(SchemaError, "association_stable requires sign_stable"):
            schema.validate_feature(good_feature(sign_stable=False))
        schema.validate_feature(good_feature(sign_stable=False, claim_level="predictive_increment"))
        unstable = good_feature(sign_stable=False, claim_level="predictive_increment")
        for level in schema.MANUAL_LEVELS:
            with self.subTest(level=level), self.assertRaisesRegex(SchemaError, f"{level} requires sign_stable"):
                schema.annotate_claim(unstable, "feature", level, ["x"], "A", "2026-10-30")

    def test_causal_supported_needs_evidence_refs(self):
        with self.assertRaisesRegex(SchemaError, "causal_supported requires non-empty evidence_refs"):
            schema.validate_feature(good_feature(claim_level="causal_supported", annotations=[note()]))

    def test_manual_path(self):
        out = schema.annotate_claim(good_feature(), "feature", "causal_supported",
                                    ["PMID:12345678"], "Layth", "2026-10-30")
        self.assertEqual(out["claim_level"], "causal_supported")
        self.assertEqual(out["evidence_refs"], ["PMID:12345678"])
        self.assertEqual(out["annotations"], [{"by": "Layth", "date": "2026-10-30", "level": "causal_supported",
                                               "previous_level": "association_stable"}])
        self.assertEqual(good_feature()["claim_level"], "association_stable")
        again = schema.annotate_claim(out, "feature", "mechanism_candidate", ["PMID:12345678", "doi:10.1/y"],
                                      "Reviewer B", "2026-10-30")
        self.assertEqual(again["evidence_refs"], ["PMID:12345678", "doi:10.1/y"])
        self.assertEqual([n["previous_level"] for n in again["annotations"]],
                         ["association_stable", "causal_supported"])
        schema.validate_feature(again)

    def test_manual_path_refusals(self):
        future = (datetime.date.today() + datetime.timedelta(days=1)).isoformat()
        cases = [
            (("feature", "causal_supported", [], "Layth", "2026-10-30"), "requires non-empty evidence_refs"),
            (("feature", "mechanism_candidate", [], "Layth", "2026-10-30"), "requires non-empty evidence_refs"),
            (("feature", "association_stable", ["x"], "Layth", "2026-10-30"), "level must be one of"),
            (("feature", "causal_supported", ["x"], "", "2026-10-30"), r"annotations\[0\].by"),
            (("feature", "causal_supported", ["x"], "Layth", "30/10/2026"), r"annotations\[0\].date"),
            (("feature", "causal_supported", ["x"], "Layth", future), "must not be in the future"),
            (("feature", "causal_supported", ["x"], "Layth", "2999-12-31"), "must not be in the future"),
            (("feature", "causal_supported", [""], "Layth", "2026-10-30"), "evidence_refs"),
            (("feature", "causal_supported", "PMID:123", "Layth", "2026-10-30"), "must be a list of strings"),
            (("feature", "causal_supported", ["x", 1], "Layth", "2026-10-30"), "must be a list of strings"),
            (("pooled", "causal_supported", ["x"], "Layth", "2026-10-30"), "only features"),
        ]
        for args, message in cases:
            with self.subTest(args=args), self.assertRaisesRegex(SchemaError, message):
                schema.annotate_claim(good_feature(), *args)

    def test_annotation_chain(self):
        """L2: previous_level must be the level actually replaced."""
        base = good_feature(claim_level="causal_supported", evidence_refs=["x"])
        schema.validate_feature({**base, "annotations": [note()]})
        cases = [
            ([note(previous_level="causal_supported")], r"annotations\[0\].previous_level: the first annotation"),
            ([note(previous_level="mechanism_candidate")], r"annotations\[0\].previous_level"),
            ([note(level="mechanism_candidate"), note(previous_level="association_stable")],
             r"annotations\[1\].previous_level: must equal the level"),
            ([note(level="mechanism_candidate")], "claim_level: must equal the last annotation's level"),
            ([note(date="2026-10-29"), note(date="2026-10-28", previous_level="causal_supported")],
             r"annotations\[1\].date: must not precede"),
            ([], "must not be empty"),
            ([note(level="association_stable")], r"annotations\[0\].level"),
            ([{**note(), "approved": True}], "approved: is not a field"),
        ]
        for notes, message in cases:
            with self.subTest(message=message), self.assertRaisesRegex(SchemaError, message):
                schema.validate_feature({**base, "annotations": notes})
        with self.assertRaisesRegex(SchemaError, "only record a manual"):
            schema.validate_feature(good_feature(annotations=[note()], evidence_refs=["x"]))

    def test_inputs_are_not_mutated(self):
        rec = good_feature(evidence_refs=["a"])
        before = copy.deepcopy(rec)
        out = schema.annotate_claim(rec, "feature", "causal_supported", ["x"], "Layth", "2026-10-30")
        self.assertEqual(rec, before)
        out["evidence_refs"].append("y")
        self.assertEqual(rec, before)


def paired_manifest(rule="common_outer_test_set_balanced", own=(H1, H2), reference="K"):
    """good_manifest() plus an arm `cap` on ecoli__amp's lasso lineage run, paired to `reference`."""
    m = good_manifest()
    m["runs"].append(good_run(id="r-cap", arm="cap", predictions="predictions/r-cap.jsonl", pairing={
        "rule": rule, "reference_arm": reference, "arm_plan_sha256": None if own is None else list(own)}))
    return m


class PairingCases(unittest.TestCase):
    """DECISIONS V-1's pairing rules on runs, contrasts and prediction rows."""

    def test_a_paired_run_names_its_rule_reference_and_own_plans(self):
        schema.validate_manifest(paired_manifest())
        schema.validate_manifest(paired_manifest(rule="k31_groups", own=None))
        schema.validate_manifest(paired_manifest(rule="common_outer_test_set"))  # older records
        for m, field in ((paired_manifest(own=None), "arm_plan_sha256"),
                         (paired_manifest(rule="common_outer_test_set", own=None), "arm_plan_sha256"),
                         (paired_manifest(rule="k31_groups"), "arm_plan_sha256"),
                         (paired_manifest(own=(H1,)), "arm_plan_sha256"),
                         (paired_manifest(reference="cap"), "reference_arm"),
                         (paired_manifest(reference="nobody"), "reference_arm"),
                         (paired_manifest(rule="nearest"), "rule")):
            with self.assertRaises(SchemaError, msg=field) as e:
                schema.validate_manifest(m)
            self.assertIn(field, str(e.exception))

    def test_an_auroc_delta_sits_only_beside_a_log_loss_decision(self):
        schema.validate_contrast(good_contrast(auroc_delta=None, auroc_ci=None))
        for over in ({"auroc_ci": None}, {"auroc_delta": None}, {"auroc_delta": 0.05},
                     {"id": "ecoli__amp/" + key(family="D", endpoint="auroc")}):
            with self.assertRaises(SchemaError, msg=over):
                schema.validate_contrast(good_contrast(**over))

    def test_a_contrast_names_one_of_the_rules(self):
        schema.validate_contrast(good_contrast(pairing="k31_groups"))
        for bad in (None, "", "nearest", 1):
            with self.assertRaises(SchemaError):
                schema.validate_contrast(good_contrast(pairing=bad))
        rec = good_contrast()
        del rec["pairing"]
        with self.assertRaises(SchemaError):
            schema.validate_contrast(rec)

    def test_paired_rows_cover_the_reference_samples_or_a_subset(self):
        base = [good_prediction(sample_id=s, repeat=r, fold_plan_sha256=LINEAGE_PLANS[r])
                for s in ("a", "b") for r in (0, 1)]
        cap = [dict(row, arm="cap") for row in base if row["sample_id"] == "a"]
        schema.validate_predictions(base + cap, paired_manifest())
        with self.assertRaises(SchemaError) as e:
            schema.validate_predictions(base + cap, paired_manifest(rule="k31_groups", own=None))
        self.assertIn("needs the same samples", str(e.exception))
        extra = [dict(row, sample_id="c", cluster="L2") for row in cap]
        with self.assertRaises(SchemaError) as e:
            schema.validate_predictions(base + cap + extra, paired_manifest())
        self.assertIn("a subset of them", str(e.exception))


if __name__ == "__main__":
    unittest.main()
