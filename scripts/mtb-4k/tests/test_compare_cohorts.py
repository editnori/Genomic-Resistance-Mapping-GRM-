"""compare_cohorts.py on synthetic schema-2 / revision-8 reports.

usage: python3 -I scripts/mtb-4k/tests/test_compare_cohorts.py
"""
import hashlib
import json
import random
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import compare_cohorts as cc  # noqa: E402

SLUG = "mycobacterium_tuberculosis"
N = 200
GENOMES = [f"G{i:04d}" for i in range(N)]
PREFIX, ADDED = GENOMES[:100], GENOMES[100:]
_PICK = set(sorted(GENOMES, key=lambda g: hashlib.sha256(g.encode()).hexdigest())[:100])
MATCHED = [g for g in GENOMES if g in _PICK]
CONFIG = {"folds": 5, "seed": 42}


def labels(drug):
    """alpha: balanced everywhere; beta: no resistant isolate after G0099;
    gamma: balanced, but its matched2k report is invalid."""
    out = {}
    for i, g in enumerate(GENOMES):
        if drug == "beta" and i >= 100:
            out[g] = False
        else:
            out[g] = i % 3 == 0
    return out


def scores_for(drug, genomes, cohort):
    rng = random.Random(f"{drug}-{cohort}")
    truth = labels(drug)
    out = []
    for g in genomes:
        p = 0.8 if truth[g] else 0.2
        s = min(1.0, max(0.0, p + rng.uniform(-0.45, 0.45)))
        out.append(None if rng.random() < 0.03 else s)
    return out


def metrics_of(truth, scores):
    m = cc.recompute(truth, scores, [0] * len(truth))
    if m["n_called"] == 0:
        return None
    return {k: m[k] for k in ("tp", "fp", "tn", "fn")} | {"n": m["n_called"], "balanced_accuracy": m["balanced_accuracy"]}


def report(drug, genomes, cohort, *, valid=True, break_metrics=False):
    truth = labels(drug)
    ids = [g for g in genomes if g in truth]
    t = [truth[g] for g in ids]
    splits = []
    for scheme in cc.SCHEMES:
        folds = [i % 5 for i in range(len(ids))]
        learners = []
        for name in cc.LEARNERS:
            sc = scores_for(f"{drug}-{scheme}-{name}", ids, cohort)
            m = metrics_of(t, sc)
            if break_metrics and name == "rf" and m:
                m["tp"] += 1
            nc_r = sum(1 for x, s in zip(t, sc) if s is None and x)
            nc_s = sum(1 for x, s in zip(t, sc) if s is None and not x)
            row = {"learner": name, "metrics": m, "scores": sc, "no_call": {"resistant": nc_r, "susceptible": nc_s}}
            if name == "lasso":
                row["lasso_fits"] = [{"path_points": 30, "path_unconverged": 0,
                                      "selected_convergence": {"converged": True},
                                      "refit_convergence": {"converged": True}} for _ in range(5)]
            learners.append(row)
        splits.append({"scheme": scheme, "folds": 5, "n_groups": 5, "test_fold": folds, "learners": learners})
    return {
        "schema": 2,
        "evaluation": {"revision": 8, "valid": valid, "reasons": [] if valid else ["final model has an unconverged two-class lasso fit"]},
        "input_fingerprint": f"{cohort}-{drug}", "eval_config": CONFIG, "id": f"{SLUG}__{drug}", "antibiotic": drug,
        "counts": {"samples": len(ids), "resistant": sum(t), "susceptible": len(t) - sum(t)},
        "samples": [{"id": g, "resistant": truth[g]} for g in ids],
        "splits": splits,
    }


def write_root(base, cohort, genomes, reports):
    sp = Path(base) / cohort / "species" / SLUG
    (sp / "results").mkdir(parents=True)
    (sp / "genomes.txt").write_text("\n".join(genomes) + "\n")
    for rep in reports:
        (sp / "results" / f"{rep['id']}.json").write_text(json.dumps(rep))
    return str(Path(base) / cohort)


class Compare(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.saved = (cc.MIN_CLASS, cc.EXPECTED_ADDED_ELIGIBLE, cc.EXPECTED_SPARSE_ADDED)
        cc.MIN_CLASS = 10
        cc.EXPECTED_ADDED_ELIGIBLE = {"alpha", "gamma"}
        cc.EXPECTED_SPARSE_ADDED = {"beta": "test"}

    def tearDown(self):
        cc.MIN_CLASS, cc.EXPECTED_ADDED_ELIGIBLE, cc.EXPECTED_SPARSE_ADDED = self.saved
        self.tmp.cleanup()

    def args(self, gamma_matched_valid=False, break_metrics=False, cross_exact=0, accept=None):
        b = self.tmp.name
        roots = {
            "prefix": write_root(b, "prefix", PREFIX, [report(d, PREFIX, "prefix") for d in ("alpha", "beta", "gamma")]),
            "added": write_root(b, "added", ADDED, [report(d, ADDED, "added", break_metrics=break_metrics) for d in ("alpha", "gamma")]),
            "full": write_root(b, "full", GENOMES, [report(d, GENOMES, "full") for d in ("alpha", "beta", "gamma")]),
            "matched2k": write_root(b, "matched2k", MATCHED, [report("alpha", MATCHED, "matched2k"), report("beta", MATCHED, "matched2k"),
                                                              report("gamma", MATCHED, "matched2k", valid=gamma_matched_valid)]),
        }
        receipts = []
        for name in roots:
            p = Path(b) / f"{name}-checker.json"
            p.write_text("{}")
            receipts.append(f"{name}={p}")
        screen = Path(b) / "screen.json"
        screen.write_text(json.dumps({"cross_cohort_exact": cross_exact, "cross_cohort_near": 0}))
        return ["--duplicates", str(screen), "--out", str(Path(b) / "out.json")] + \
            [f"--cohort={k}={v}" for k, v in roots.items()] + [f"--checker-receipt={r}" for r in receipts] + \
            (["--accept-reviewed-duplicates", accept] if accept else [])

    def parse(self, argv):
        sys.argv = ["compare_cohorts.py"] + argv
        return cc.main()

    def test_statistics(self):
        self.assertAlmostEqual(cc.clopper_pearson(0, 10)[1], 0.30850, places=4)
        lo, hi = cc.clopper_pearson(5, 10)
        self.assertAlmostEqual(lo, 0.18709, places=4)
        self.assertAlmostEqual(hi, 0.81291, places=4)
        self.assertAlmostEqual(cc.clopper_pearson(10, 10)[0], 0.69150, places=4)
        self.assertIsNone(cc.clopper_pearson(0, 0))
        self.assertAlmostEqual(cc.mcnemar_exact(0, 5), 0.0625)
        self.assertEqual(cc.mcnemar_exact(0, 0), 1.0)
        self.assertEqual(cc.holm({"a": 0.01, "b": 0.04, "c": 0.03}), {"a": 0.03, "b": 0.06, "c": 0.06})

    def test_invalid_in_one_cohort_is_excluded_there_and_stops_for_review(self):
        rc = self.parse(self.args())
        self.assertEqual(rc, 3)
        out = json.loads((Path(self.tmp.name) / "out.json").read_text())
        self.assertEqual(list(out["excluded"]), ["gamma"])
        self.assertEqual(out["excluded"]["gamma"][0]["cohort"], "matched2k")
        drugs = lambda name: [r["drug"] for r in out["comparisons"][name]["rows"]]
        # gamma stays in comparisons without matched2k, and leaves the two with it.
        self.assertIn("gamma", drugs("full_vs_prefix"))
        self.assertIn("gamma", drugs("added_vs_prefix"))
        self.assertNotIn("gamma", drugs("matched2k_vs_prefix"))
        self.assertNotIn("gamma", drugs("full_vs_matched2k"))
        self.assertNotIn("gamma", [r["drug"] for r in out["size_matched_decomposition"]])
        # beta is not interpretable before/after and not eligible in added.
        self.assertEqual(list(out["not_interpretable_for_before_after"]), ["beta"])
        self.assertIn("beta", out["not_eligible"]["added"])
        row = next(r for r in out["comparisons"]["full_vs_prefix"]["rows"] if r["drug"] == "beta")
        self.assertFalse(row["before_after_interpretable"])
        alpha = next(r for r in out["comparisons"]["full_vs_prefix"]["rows"] if r["drug"] == "alpha")
        paired = alpha["lineage"]["paired"]["rf"]
        self.assertEqual(paired["isolates"], len([g for g in PREFIX if g in labels("alpha")]))
        self.assertEqual(sum(paired["transitions"].values()), paired["isolates"])
        self.assertIn("mcnemar_p_holm", paired)
        self.assertEqual(out["primary_split"], "lineage")

    def test_all_valid_completes(self):
        self.assertEqual(self.parse(self.args(gamma_matched_valid=True)), 0)
        out = json.loads((Path(self.tmp.name) / "out.json").read_text())
        self.assertFalse(out["excluded"])
        self.assertEqual(len(out["size_matched_decomposition"]), 3)

    def test_a_report_whose_metrics_do_not_match_its_scores_is_invalid(self):
        self.assertEqual(self.parse(self.args(gamma_matched_valid=True, break_metrics=True)), 3)
        out = json.loads((Path(self.tmp.name) / "out.json").read_text())
        reasons = out["excluded"]["alpha"][0]["reasons"]
        self.assertTrue(any("stored confusion matrix differs" in r for r in reasons), reasons)

    def test_exact_cross_cohort_duplicates_stop_unless_accepted(self):
        self.assertEqual(self.parse(self.args(gamma_matched_valid=True, cross_exact=1)), 2)
        screen = Path(self.tmp.name) / "screen.json"
        digest = hashlib.sha256(screen.read_bytes()).hexdigest()
        (Path(self.tmp.name) / "out.json").unlink(missing_ok=True)
        self.tmp.cleanup()
        self.tmp = tempfile.TemporaryDirectory()
        args = self.args(gamma_matched_valid=True, cross_exact=1)  # same screen bytes
        self.assertEqual(self.parse(args + ["--accept-reviewed-duplicates", digest]), 0)

    def test_missing_checker_receipt_stops(self):
        args = [a for a in self.args(gamma_matched_valid=True) if not a.startswith("--checker-receipt=matched2k")]
        self.assertEqual(self.parse(args), 2)


if __name__ == "__main__":
    unittest.main()
