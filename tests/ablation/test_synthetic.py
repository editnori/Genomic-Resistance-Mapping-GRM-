"""The synthetic ablation cohort generator: format, config checks,
determinism across workers, dict order and PYTHONHASHSEED, and the planted
truths, checked with a numpy-only ridge logistic regression.

usage: python3 -I -B -m unittest discover -s tests/ablation -t .
The module is loaded by path, so the scripts.ablation package need not exist.
"""
import collections
import copy
import fractions
import functools
import hashlib
import importlib.util
import json
import math
import os
import pathlib
import subprocess
import sys
import types
import unittest

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
MODULE = ROOT / "scripts" / "ablation" / "synthetic.py"
_spec = importlib.util.spec_from_file_location("grm_ablation_synthetic", MODULE)
synthetic = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(synthetic)

SEED = 20261030
# Reference checks average over these seeds so one draw cannot carry them.
REFERENCE_SEEDS = (SEED, 1, 2, 3)


@functools.lru_cache(maxsize=None)
def _cached(config_json, seed):
    return synthetic.to_json(synthetic.generate(json.loads(config_json), seed))


def cohort(config=None, seed=SEED):
    """A fresh copy of a cached cohort, so tests may mutate it."""
    return json.loads(_cached(json.dumps(config or {}, sort_keys=True), seed))


def labels(c):
    return np.array([s["y"] for s in c["samples"]], dtype=float)


def lineage_index(c):
    ids = sorted({s["lineage"] for s in c["samples"]})
    return np.array([ids.index(s["lineage"]) for s in c["samples"]])


def column(c, name):
    return next(col for col in c["columns"] if col["name"] == name)


def missing_of(raw):
    (kind, body), = raw.items()
    if kind == "Binary":
        return [not a for a in body["assessable"]]
    return [v is None for v in body["values" if kind == "Count" else "codes"]]


def matrix(c, block, missing=False):
    """Features of one block: present flags, counts (missing as 0) and one-hot
    categorical codes; or, with missing=True, not-assessable indicators."""
    feats = []
    for col in c["columns"]:
        if col["block"] != block:
            continue
        (kind, body), = col["raw"].items()
        if missing:
            feats.append(missing_of(col["raw"]))
        elif kind == "Binary":
            feats.append(body["present"])
        elif kind == "Count":
            feats.append([0 if v is None else v for v in body["values"]])
        else:
            feats.extend([v == k for v in body["codes"]] for k in range(len(body["levels"])))
    return np.array(feats, dtype=float).T


def auroc(score, y):
    score, y = np.asarray(score, dtype=float), np.asarray(y, dtype=bool)
    n1 = y.sum()
    n0 = len(y) - n1
    _, inv, cnt = np.unique(score, return_inverse=True, return_counts=True)
    ranks = (np.cumsum(cnt) - (cnt - 1) / 2.0)[inv]
    return float((ranks[y].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def fit_logistic(X, y, l2=1.0):
    """Ridge logistic regression by Newton's method; the intercept is unpenalized."""
    A = np.hstack([np.ones((len(X), 1)), X])
    pen = np.full(A.shape[1], l2)
    pen[0] = 0.0
    w = np.zeros(A.shape[1])
    for _ in range(50):
        p = 1.0 / (1.0 + np.exp(-A @ w))
        step = np.linalg.solve((A * (p * (1 - p))[:, None]).T @ A + np.diag(pen), A.T @ (p - y) + pen * w)
        w -= step
        if np.abs(step).max() < 1e-8:
            break
    return w


def cv_auroc(X, y, folds):
    """Mean test-fold AUROC of the logistic model over the folds."""
    out = []
    for f in np.unique(folds):
        train, test = folds != f, folds == f
        w = fit_logistic(X[train], y[train])
        out.append(auroc(w[0] + X[test] @ w[1:], y[test]))
    return float(np.mean(out))


def within_lineage_auroc(X, y, folds, lineage):
    """Mean AUROC within each lineage of the out-of-fold logistic scores:
    the estimand cv_auroc gives under leave-one-lineage-out, so the two splits
    can be compared like for like."""
    score = np.empty(len(y))
    for f in np.unique(folds):
        train, test = folds != f, folds == f
        w = fit_logistic(X[train], y[train])
        score[test] = w[0] + X[test] @ w[1:]
    return float(np.mean([auroc(score[lineage == k], y[lineage == k]) for k in np.unique(lineage)]))


def random_folds(y, k=4):
    """Stratified random folds from a fixed test-side generator."""
    rng = np.random.default_rng(0)
    folds = np.empty(len(y), dtype=int)
    for cls in (0, 1):
        idx = np.flatnonzero(y == cls)
        folds[rng.permutation(idx)] = np.arange(len(idx)) % k
    return folds


@functools.lru_cache(maxsize=None)
def _fisher_weights(n1, n0, k):
    """Exact integer hypergeometric weights C(n1, x) C(n0, k - x) over the
    support; each binomial is stepped exactly from the one before."""
    lo, hi = max(0, k - n0), min(k, n1)
    left, right = math.comb(n1, lo), math.comb(n0, k - lo)
    weights = [left * right]
    for x in range(lo, hi):
        left = left * (n1 - x) // (x + 1)
        right = right * (k - x) // (n0 - k + x + 1)
        weights.append(left * right)
    return lo, tuple(weights)


def brute_fisher(a, b, c, d, tol=fractions.Fraction(1, 10 ** 7)):
    """Two-sided Fisher exact p from exact integer hypergeometric weights,
    keeping tables with weight <= (1 + tol) times the observed one; the
    exact Fraction, as a float."""
    lo, weights = _fisher_weights(a + b, c + d, a + c)
    cut = weights[a - lo] * (1 + tol)
    return float(fractions.Fraction(sum(w for w in weights if w <= cut), sum(weights)))


def unit_table(c, miss_any):
    """Section 11 units from the cohort JSON: a duplicate group per class, or
    a lone sample; the 2x2 table (R missing, R complete, S missing, S complete)."""
    units = {}
    for i, s in enumerate(c["samples"]):
        key = (s["y"], s["duplicate"] if s["duplicate"] is not None else f"#{i}")
        units[key] = units.get(key, False) or bool(miss_any[i])
    count = collections.Counter((cls, m) for (cls, _), m in units.items())
    return count[(1, True)], count[(1, False)], count[(0, True)], count[(0, False)]


def effect(c, scenario):
    return next(e for e in c["truth"]["effects"] if e["scenario"] == scenario)


class Format(unittest.TestCase):
    """The cohort layout, column format and truth manifest."""

    def setUp(self):
        self.c = cohort()

    def test_validates_and_has_the_documented_shape(self):
        c = self.c
        synthetic.validate_cohort(c)
        self.assertEqual((c["format"], c["format_version"], c["seed"]), ("grm-ablation-synthetic", 1, SEED))
        self.assertEqual(len(c["samples"]), 400)
        self.assertEqual([b["id"] for b in c["blocks"]], ["C", "H", "P", "N", "X", "K", "noise"])
        self.assertEqual({b["id"]: b["kind"] for b in c["blocks"]}, synthetic.BLOCK_KINDS)
        self.assertEqual(c["config"], synthetic.resolve_config(None))
        counts = {}
        for s in c["samples"]:
            counts[s["lineage"]] = counts.get(s["lineage"], 0) + 1
        self.assertEqual(counts, {f"L{k}": 50 for k in range(1, 9)})

    def test_x_is_unavailable_with_no_columns(self):
        x = next(b for b in self.c["blocks"] if b["id"] == "X")
        self.assertEqual((x["version"], x["source_sha256"]), ("", ""))
        self.assertEqual(x["unavailable_reason"], synthetic.X_UNAVAILABLE)
        self.assertNotIn("X", {col["block"] for col in self.c["columns"]})
        self.assertTrue(all(b["unavailable_reason"] is None for b in self.c["blocks"] if b["id"] != "X"))

    def test_raw_columns_use_the_one_key_interchange_form(self):
        kinds = {}
        for col in self.c["columns"]:
            self.assertEqual(set(col), {"name", "block", "raw"})
            (kind, body), = col["raw"].items()
            kinds.setdefault(col["block"], set()).add(kind)
            self.assertEqual(set(body), {"Binary": {"present", "assessable"}, "Count": {"values"},
                                         "Categorical": {"codes", "levels"}}[kind])
        self.assertEqual(kinds, {"C": {"Binary"}, "H": {"Binary", "Categorical"}, "P": {"Binary"},
                                 "N": {"Count"}, "K": {"Binary"}, "noise": {"Binary", "Count"}})
        # Each column round-trips through canonical JSON unchanged.
        for col in self.c["columns"]:
            self.assertEqual(json.loads(json.dumps(col["raw"])), col["raw"])

    def test_column_names_follow_h_check_block(self):
        """Section 14: every name starts with "<block>:", as H's check_block
        (blocks.rs:830-836) requires, uses / between levels, and has no
        backtick or whitespace."""
        for config in ({}, {"missingness_shortcut": {"block": "H"}, "noise": {"enabled": False}}):
            names = [col["name"] for col in cohort(config)["columns"]]
            self.assertEqual(len(set(names)), len(names))
            for col in cohort(config)["columns"]:
                self.assertTrue(col["name"].startswith(col["block"] + ":"), col["name"])
                self.assertRegex(col["name"], r"^(C|H|P|N|K|noise):[A-Za-z0-9]+(/[A-Za-z0-9]+)?$")
        first = {col["block"]: col["name"] for col in reversed(self.c["columns"])}
        self.assertEqual(first, {"C": "C:causal", "H": "H:fam001", "P": "P:gene001", "N": "N:copy001",
                                 "K": "K:pat0001", "noise": "noise:bin001"})
        self.assertIn("H:fam001/allele", [col["name"] for col in self.c["columns"]])
        self.assertIn("noise:cnt001", [col["name"] for col in self.c["columns"]])
        self.assertEqual(synthetic.display_name("H.fam001.allele"), "H:fam001/allele")

    def test_hierarchy_binary_and_allele_agree(self):
        for col in self.c["columns"]:
            if col["block"] == "H" and "Binary" in col["raw"]:
                parent = col["raw"]["Binary"]
                child = column(self.c, col["name"] + "/allele")["raw"]["Categorical"]
                self.assertEqual(child["levels"][0], "absent")
                for p, a, code in zip(parent["present"], parent["assessable"], child["codes"]):
                    self.assertEqual(a, code is not None)
                    self.assertEqual(p, code is not None and code > 0)

    def test_k_and_noise_are_always_assessable_and_the_rest_lose_a_few(self):
        for col in self.c["columns"]:
            miss = np.mean(missing_of(col["raw"]))
            if col["block"] in ("K", "noise"):
                self.assertEqual(miss, 0.0, col["name"])
        for block in ("C", "H", "N"):
            rate = matrix(self.c, block, missing=True).mean()
            self.assertTrue(0.005 < rate < 0.05, (block, rate))

    def test_truth_names_every_effect_with_canaries_for_j_k_and_integration(self):
        t = self.c["truth"]
        self.assertEqual([e["scenario"] for e in t["effects"]], list(synthetic.SCENARIOS))
        self.assertEqual(t["disabled"], [])
        names = {col["name"] for col in self.c["columns"]}
        for e in t["effects"]:
            self.assertTrue(set(e["columns"]) <= names, e["scenario"])
            self.assertEqual([x.split(":")[0] for x in e["canaries"]], ["J", "K", "integration"])
            self.assertEqual(e["params"], self.c["config"][e["scenario"]])
        self.assertEqual(effect(self.c, "portable_signal")["columns"], ["C:causal"])
        proxy = effect(self.c, "group_proxy")
        self.assertEqual(len(proxy["columns"]), 8 * 4)
        self.assertEqual(sorted(sum(proxy["markers_by_lineage"].values(), [])), sorted(proxy["columns"]))
        self.assertEqual(effect(self.c, "missingness_shortcut")["columns"],
                         [col["name"] for col in self.c["columns"] if col["block"] == "P"])
        planted = {n for e in t["effects"] for n in e["columns"]}
        self.assertEqual(set(t["null_columns"]), names - planted)
        self.assertTrue(all(n.startswith(("C:call", "H:", "N:", "K:")) for n in t["null_columns"]))

    def test_markers_tag_their_lineage(self):
        lin = lineage_index(self.c)
        for k, (lid, cols) in enumerate(sorted(effect(self.c, "group_proxy")["markers_by_lineage"].items())):
            for name in cols:
                present = np.array(column(self.c, name)["raw"]["Binary"]["present"])
                agree = np.mean(present == (lin == k))
                self.assertGreater(agree, 0.95, (lid, name))

    def test_validate_cohort_refuses_malformed_cohorts(self):
        def broken(mutate):
            c = cohort()
            mutate(c)
            return c
        cases = {
            "present implies assessable": lambda c: c["columns"][0]["raw"]["Binary"].update(
                present=[True] * 400, assessable=[False] * 400),
            "values": lambda c: column(c, "N:copy001")["raw"]["Count"]["values"].__setitem__(0, -1),
            "codes": lambda c: column(c, "H:fam001/allele")["raw"]["Categorical"]["codes"].__setitem__(0, 9),
            "entries": lambda c: column(c, "P:gene001")["raw"]["Binary"]["present"].pop(),
            "y": lambda c: c["samples"][0].update(y=2),
            "unique": lambda c: c["samples"][1].update(sample_id=c["samples"][0]["sample_id"]),
            "exactly one": lambda c: c["columns"][0]["raw"].update(Count={"values": []}),
            "duplicate group": lambda c: c["samples"][0].update(duplicate="D999"),
            "source_sha256 must hash": lambda c: c["blocks"][2].update(source_sha256="0" * 64),
            "priority": lambda c: c["blocks"][0].update(priority=5),
            "keys must be id": lambda c: c["blocks"][0].update(columns=[]),
            "grouped by block": lambda c: c["columns"].append(c["columns"].pop(0)),
            "available block": lambda c: c["columns"][0].update(block="X"),
            "flag": lambda c: c["missingness"][0].update(flag=not c["missingness"][0]["flag"]),
            "oriented": lambda c: c["missingness"][0].update(indicator_auroc_oriented=0.5),
            "u32": lambda c: column(c, "N:copy001")["raw"]["Count"]["values"].__setitem__(0, 2 ** 32),
            "index levels": lambda c: column(c, "H:fam001/allele")["raw"]["Categorical"]["codes"].__setitem__(0, -1),
            "single-class": lambda c: [s.update(y=0) for s in c["samples"]],
            "section 14": lambda c: column(c, "C:causal").update(name="C.causal"),
            "name: must": lambda c: column(c, "P:gene001").update(name="P:gene 001"),
            "no backtick": lambda c: column(c, "P:gene001").update(name="P:`gene001`"),
            "<block>:<name>": lambda c: column(c, "P:gene001").update(name="P:"),
        }
        for needle, mutate in cases.items():
            with self.subTest(needle):
                with self.assertRaisesRegex(ValueError, needle):
                    synthetic.validate_cohort(broken(mutate))
        with self.assertRaises(ValueError):
            synthetic.validate_cohort({"format": "x"})


class Duplicates(unittest.TestCase):

    def test_pairs_are_same_lineage_and_identical_in_every_column(self):
        c = cohort()
        groups = effect(c, "duplicates")["groups"]
        self.assertEqual(len(groups), 20)  # floor(0.1 * 400 / 2)
        self.assertEqual(groups, sorted(groups))
        flat = [i for g in groups for i in g]
        self.assertEqual(len(flat), len(set(flat)))
        by_id = {}
        for i, s in enumerate(c["samples"]):
            if s["duplicate"] is not None:
                by_id.setdefault(s["duplicate"], []).append(i)
        self.assertEqual(sorted(by_id.values()), groups)
        for i, j in groups:
            self.assertLess(i, j)
            a, b = c["samples"][i], c["samples"][j]
            self.assertEqual((a["lineage"], a["y"]), (b["lineage"], b["y"]))
            for col in c["columns"]:
                (_, body), = col["raw"].items()
                for key, vec in body.items():
                    if key != "levels":
                        self.assertEqual(vec[i], vec[j], (col["name"], i, j))
        self.assertEqual(effect(c, "duplicates")["realized"]["fraction"], 0.1)

    def test_fraction_is_configurable_and_off_means_no_groups(self):
        c = cohort({"duplicates": {"fraction": 0.3}})
        self.assertEqual(len(effect(c, "duplicates")["groups"]), 60)
        off = cohort({"duplicates": {"enabled": False}})
        self.assertTrue(all(s["duplicate"] is None for s in off["samples"]))
        self.assertEqual(off["truth"]["disabled"], ["duplicates"])


class Config(unittest.TestCase):

    def test_bad_configs_are_refused_with_the_field_named(self):
        cases = {
            "config.bogus": {"bogus": 1},
            "config.noise.bogus": {"noise": {"bogus": 1}},
            "config.n_samples": {"n_samples": 10.0},
            "config.n_lineages": {"n_lineages": 1},
            "config.portable_signal.frequency": {"portable_signal": {"frequency": 1.0}},
            "config.portable_signal.enabled": {"portable_signal": {"enabled": 1}},
            "config.group_proxy.marker_flip": {"group_proxy": {"marker_flip": float("nan")}},
            "config.missingness_shortcut.block": {"missingness_shortcut": {"block": "K"}},
            "config.missingness_shortcut.rate_R": {"missingness_shortcut": {"rate_R": True}},
            "config.duplicates.fraction": {"duplicates": {"fraction": -0.1}},
            "config.noise:": {"noise": {"n_binary": 0, "n_count": 0}},
            "config.blocks": {"blocks": 3},
        }
        for needle, cfg in cases.items():
            with self.subTest(needle):
                with self.assertRaises(ValueError) as caught:
                    synthetic.generate(cfg, 0)
                self.assertIn(needle, str(caught.exception))
        with self.assertRaisesRegex(ValueError, "block P has no columns"):
            synthetic.resolve_config({"blocks": {"P": {"n_columns": 0}}})
        for seed in (-1, True, 1.0, 2 ** 63):
            with self.assertRaisesRegex(ValueError, "seed"):
                synthetic.generate(None, seed)
        with self.assertRaisesRegex(ValueError, "workers"):
            synthetic.generate(None, 0, workers=0)

    def test_numbers_are_normalized_so_int_and_float_agree(self):
        a = synthetic.generate({"portable_signal": {"logit": 2}, "duplicates": {"fraction": 0}}, 5)
        b = synthetic.generate({"portable_signal": {"logit": 2.0}, "duplicates": {"fraction": 0.0}}, 5)
        self.assertEqual(synthetic.to_json(a), synthetic.to_json(b))

    def test_minimal_cohort_without_scenarios(self):
        cfg = {s: {"enabled": False} for s in synthetic.SCENARIOS}
        c = synthetic.generate(cfg, 3)
        synthetic.validate_cohort(c)
        self.assertEqual(c["truth"]["effects"], [])
        self.assertEqual(c["truth"]["disabled"], list(synthetic.SCENARIOS))
        self.assertNotIn("noise", [b["id"] for b in c["blocks"]])
        self.assertNotIn("C:causal", [col["name"] for col in c["columns"]])


def _reverse_keys(obj):
    if isinstance(obj, dict):
        return {k: _reverse_keys(obj[k]) for k in reversed(list(obj))}
    return obj


class Determinism(unittest.TestCase):

    def test_same_config_and_seed_give_identical_canonical_bytes(self):
        a = synthetic.to_json(synthetic.generate(None, SEED))
        b = synthetic.to_json(synthetic.generate(None, SEED))
        self.assertEqual(a, b)
        self.assertEqual(a, json.dumps(json.loads(a), sort_keys=True, separators=(",", ":")))
        self.assertNotEqual(a, synthetic.to_json(synthetic.generate(None, SEED + 1)))
        self.assertEqual(synthetic.sha256(synthetic.generate(None, SEED)), hashlib.sha256(a.encode()).hexdigest())

    def test_worker_count_does_not_change_the_bytes(self):
        ref = synthetic.sha256(synthetic.generate(None, SEED, workers=1))
        for workers in (2, 4, 7):
            self.assertEqual(synthetic.sha256(synthetic.generate(None, SEED, workers=workers)), ref, workers)

    def test_config_dict_order_does_not_change_the_bytes(self):
        cfg = {"missingness_shortcut": {"block": "H", "rate_S": 0.4}, "n_samples": 240,
               "blocks": {"K": {"n_background": 10}, "C": {"n_null": 3}}}
        self.assertEqual(synthetic.sha256(synthetic.generate(cfg, 9)),
                         synthetic.sha256(synthetic.generate(_reverse_keys(cfg), 9)))
        self.assertEqual(synthetic.sha256(synthetic.generate(None, 9)),
                         synthetic.sha256(synthetic.generate(_reverse_keys(synthetic.DEFAULT_CONFIG), 9)))

    def test_pythonhashseed_does_not_change_the_bytes(self):
        """Subprocesses without -I/-E so PYTHONHASHSEED applies; the first
        command proves the two hash seeds really differ."""
        ref = synthetic.sha256(synthetic.generate(None, SEED))
        probes, digests = [], []
        for hash_seed, workers in (("0", "1"), ("4242", "3")):
            env = {**os.environ, "PYTHONHASHSEED": hash_seed}
            env.pop("PYTHONPATH", None)
            probes.append(subprocess.run([sys.executable, "-B", "-s", "-c", "print(hash('grm-task-l'))"],
                                         env=env, capture_output=True, text=True, check=True).stdout)
            digests.append(subprocess.run(
                [sys.executable, "-B", "-s", str(MODULE), "--seed", str(SEED), "--workers", workers, "--sha256"],
                env=env, cwd=ROOT, capture_output=True, text=True, check=True, timeout=60).stdout.strip())
        self.assertNotEqual(probes[0], probes[1])
        self.assertEqual(digests, [ref, ref])

    def test_turning_a_scenario_off_leaves_unrelated_columns_identical(self):
        base = {col["name"]: col for col in cohort()["columns"]}
        no_noise = cohort({"noise": {"enabled": False}})
        self.assertEqual(no_noise["samples"], cohort()["samples"])
        self.assertEqual([col for col in no_noise["columns"]],
                         [col for col in cohort()["columns"] if col["block"] != "noise"])
        no_proxy = cohort({"group_proxy": {"enabled": False}})
        proxy_cols = set(effect(cohort(), "group_proxy")["columns"])
        background = [col["raw"] for col in base.values() if col["block"] == "K" and col["name"] not in proxy_cols]
        self.assertEqual([col["raw"] for col in no_proxy["columns"] if col["block"] == "K"], background)
        # Without the causal call y changes, and with it P's y-dependent
        # assessable masks; every other column keeps its bytes.
        no_portable = cohort({"portable_signal": {"enabled": False}})
        self.assertNotIn("C:causal", [col["name"] for col in no_portable["columns"]])
        self.assertNotEqual(labels(no_portable).tolist(), labels(cohort()).tolist())
        for col in no_portable["columns"]:
            if col["block"] != "P":
                self.assertEqual(col["raw"], base[col["name"]]["raw"], col["name"])


def load_k_schema():
    """Task K's schema.py, read-only: this repo's copy once K is merged, else
    the file named by GRM_TASK_K_SCHEMA. Compiled from source so no bytecode
    is written beside it."""
    for candidate in (ROOT / "scripts" / "ablation" / "schema.py", os.environ.get("GRM_TASK_K_SCHEMA")):
        if candidate and pathlib.Path(candidate).is_file():
            path = pathlib.Path(candidate)
            module = types.ModuleType("grm_ablation_schema_k")
            module.__file__ = str(path)
            exec(compile(path.read_text(encoding="utf-8"), str(path), "exec"), module.__dict__)
            return module
    return None


def k_manifest(c, schema):
    """A minimal K run manifest around the cohort's blocks and missingness
    audit: one ok run per available block and an unavailable run for X."""
    ds, plan = "synthetic", ["a" * 64]
    runs = []
    for b in c["blocks"]:
        ok = b["unavailable_reason"] is None
        runs.append({"id": f"r-{b['id']}", "dataset": ds, "arm": b["id"], "blocks": [b["id"]],
                     "learner": "ridge", "split": "lineage", "repeats": 1,
                     "fold_plan_sha256": plan if ok else None, "status": "ok" if ok else "unavailable",
                     "reason": None if ok else f"block {b['id']} unavailable",
                     "metrics": {"auroc": 0.5} if ok else {}, "fit_count": 8 if ok else 0,
                     "seconds": 1.0 if ok else 0, "predictions": f"predictions/{b['id']}.jsonl" if ok else None})
    proxy = [e for e in c["truth"]["effects"] if e["scenario"] == "group_proxy"]
    return {
        "schema_version": schema.SCHEMA_VERSION, "protocol_sha256": "3" * 64,
        "source_commit": "c90ebaf48998dc244806ac435c3706bfd87c9a9b", "binary_sha256": "4" * 64,
        "blocks": c["blocks"],
        "datasets": [{"id": ds, "species": "Synthetic cohort", "partition": "development",
                      "n_samples": len(c["samples"]), "group_counts": {"lineage": c["config"]["n_lineages"]}}],
        "missingness": [{"dataset": ds, **r} for r in c["missingness"]],
        "group_proxy_audit": [{"dataset": ds, "block": "K", "metric": "lineage_auroc",
                               "value": e["realized"]["lineage_auroc"], "flag": True} for e in proxy],
        "runs": runs,
        "compute": {"planned_fit_units": 8 * len(runs), "actual_seconds": 1.0, "tiers_completed": ["smoke"]},
    }


# Canonical sha256 of generate(config, seed), re-pinned for contract v1.7. A
# change here means the cohort bytes changed: say why in the change.
PINNED = {
    ("{}", SEED): "07fe1a44dda348da8a4a46a07acf463853105f565bda9b65a4215da9a6cae5cf",
    ("{}", 7): "8fdbbc60618f53b06978e83c9dff789e5c419131876ebd5e33bf182e7ae87d55",
    ('{"missingness_shortcut": {"block": "H"}, "noise": {"enabled": false}}', 11): "1d0ed8b1f6e256fa8c1d8784008f86e58ec485fe7bc630b49f93d5e68c6f59a9",
}


class ContractV11(unittest.TestCase):
    """Block specs (section 5), the missingness audit (section 6), canonical
    JSON (section 7) and the section 8 rulings, plus K's own validators."""

    def test_block_specs(self):
        c = cohort()
        expected = {"C": ("CuratedExact", 0), "H": ("Hierarchical", 1), "P": ("Presence", 2), "N": ("Count", 3),
                    "X": ("Context", 4), "K": ("Anonymous", 5), "noise": ("Anonymous", 5)}
        for b in c["blocks"]:
            self.assertEqual(set(b), set(synthetic.BLOCK_SPEC_KEYS))
            self.assertEqual((b["kind"], b["priority"]), expected[b["id"]])
            self.assertEqual(b["screen_cap"], 4000 if b["id"] == "K" else None)
            self.assertIsNone(b["catalog_date"])
            self.assertIs(b["circularity_risk"], False)
            if b["id"] == "X":
                continue
            own = [col for col in c["columns"] if col["block"] == b["id"]]
            text = json.dumps(own, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)
            self.assertEqual(b["version"], "synthetic-v1")
            self.assertEqual(b["source_sha256"], hashlib.sha256(text.encode("utf-8")).hexdigest())
            self.assertRegex(b["source_sha256"], r"^[0-9a-f]{64}$")
            self.assertIsNone(b["unavailable_reason"])

    def test_missingness_audit_records(self):
        for seed in REFERENCE_SEEDS:
            c = cohort(seed=seed)
            y = labels(c).astype(bool)
            audit = c["missingness"]
            self.assertEqual([r["block"] for r in audit], ["C", "H", "P", "N"])
            for r in audit:
                self.assertEqual(set(r), set(synthetic.MISSINGNESS_KEYS))
                # Units: one per column, except an H family counts once.
                units = np.array([missing_of(col["raw"]) for col in c["columns"] if col["block"] == r["block"]
                                  and not (r["block"] == "H" and "Categorical" in col["raw"])]).T
                k1, k0, n1, n0 = units[y].sum(), units[~y].sum(), units[y].size, units[~y].size
                self.assertAlmostEqual(r["rate_R"], k1 / n1, places=6)
                self.assertAlmostEqual(r["rate_S"], k0 / n0, places=6)
                table = unit_table(c, units.any(axis=1))
                exact = brute_fisher(*table)  # 6 significant digits: relative error <= 5e-6
                self.assertLessEqual(abs(r["rate_diff_p"] - exact), 5e-6 * exact, (r, table))
                self.assertTrue(0 <= r["rate_diff_p"] <= 1)
                self.assertAlmostEqual(r["indicator_auroc"], auroc(units.sum(axis=1), y), places=6)
                self.assertAlmostEqual(r["indicator_auroc_oriented"],
                                       max(r["indicator_auroc"], 1 - r["indicator_auroc"]), places=12)
                rule = r["indicator_auroc_oriented"] > 0.60 or r["rate_diff_p"] < 0.01
                self.assertIs(r["flag"], rule)
                # Planted in P only; the base-rate blocks are nulls on these seeds.
                self.assertIs(r["flag"], r["block"] == "P", (seed, r))

    def test_missingness_audit_follows_the_config(self):
        none = cohort({"base_missing_rate": 0, "missingness_shortcut": {"enabled": False}})
        self.assertEqual(none["missingness"], [])
        only = cohort({"base_missing_rate": 0, "missingness_shortcut": {"block": "N"}})
        self.assertEqual([r["block"] for r in only["missingness"]], ["N"])
        self.assertTrue(only["missingness"][0]["flag"])
        # A small rate difference over many units: p sets the flag, AUROC does not.
        by_p = synthetic.generate({"n_samples": 2000, "base_missing_rate": 0.0,
                                   "blocks": {"P": {"n_columns": 1}},
                                   "missingness_shortcut": {"rate_R": 0.08, "rate_S": 0.04}}, 2)
        (r,) = by_p["missingness"]
        self.assertLessEqual(r["indicator_auroc_oriented"], 0.60, r)
        self.assertLess(r["rate_diff_p"], 0.01, r)
        self.assertIs(r["flag"], True)
        quiet = synthetic.generate({"n_samples": 16, "n_lineages": 2, "base_missing_rate": 0.0,
                                    "missingness_shortcut": {"rate_R": 0.0, "rate_S": 0.0},
                                    "duplicates": {"enabled": False}}, 1)
        (r,) = quiet["missingness"]
        self.assertEqual((r["rate_R"], r["rate_S"], r["rate_diff_p"], r["flag"]), (0.0, 0.0, None, False))
        self.assertEqual((r["indicator_auroc"], r["indicator_auroc_oriented"]), (0.5, 0.5))

    def test_noise_block_spec_follows_section_8(self):
        noise = next(b for b in cohort()["blocks"] if b["id"] == "noise")
        self.assertEqual({k: noise[k] for k in ("id", "kind", "priority", "screen_cap", "version")},
                         {"id": "noise", "kind": "Anonymous", "priority": 5, "screen_cap": None,
                          "version": "synthetic-v1"})

    def test_canonical_json_refuses_non_string_keys(self):
        for bad in ({1: "a"}, {"a": [{"b": {2: 0}}]}, [{"x": 1}, {(1, 2): 0}]):
            with self.assertRaisesRegex(ValueError, "is not a string"):
                synthetic.to_json(bad)
        self.assertEqual(synthetic.to_json({"b": 1, "a": [True, None, 0.5]}), '{"a":[true,null,0.5],"b":1}')

    def test_pinned_hashes(self):
        for (config_json, seed), digest in PINNED.items():
            with self.subTest(config=config_json, seed=seed):
                self.assertEqual(synthetic.sha256(synthetic.generate(json.loads(config_json), seed)), digest)

    def test_cohorts_validate_against_task_k_schema(self):
        schema = load_k_schema()
        if schema is None:
            self.skipTest("Task K's schema.py is not in this repo; set GRM_TASK_K_SCHEMA to its path")
        configs = [{}, {"missingness_shortcut": {"block": "H"}, "noise": {"enabled": False}},
                   {s: {"enabled": False} for s in synthetic.SCENARIOS}]
        for config in configs:
            for seed in REFERENCE_SEEDS[:2]:
                with self.subTest(config=config, seed=seed):
                    schema.validate_manifest(k_manifest(cohort(config, seed), schema))
        good = k_manifest(cohort(), schema)
        bad_flag = copy.deepcopy(good)
        bad_flag["missingness"][0]["flag"] = not bad_flag["missingness"][0]["flag"]
        bad_priority = copy.deepcopy(good)
        bad_priority["blocks"][0]["priority"] = 3
        for bad, needle in ((bad_flag, "flag"), (bad_priority, "priority")):
            with self.assertRaises(schema.SchemaError) as caught:
                schema.validate_manifest(bad)
            self.assertTrue(any(needle in e for e in caught.exception.errors), caught.exception.errors)


class ContractV14(unittest.TestCase):
    """Section 11 rulings and the Task L review's fix round 1."""

    def test_categorical_codes_above_u16_are_refused(self):
        def col(codes):
            return synthetic._Column("H:fam001/allele", "H", "Categorical", np.array(codes),
                                     np.zeros(len(codes), dtype=bool), ["absent", "a1"])
        self.assertEqual(col([0, 65535]).raw()["Categorical"]["codes"], [0, 65535])
        with self.assertRaisesRegex(ValueError, "65535"):
            col([0, 65536]).raw()
        c = cohort()
        allele = column(c, "H:fam001/allele")["raw"]["Categorical"]
        allele["levels"] += [f"x{k}" for k in range(65537 - len(allele["levels"]))]
        allele["codes"][0] = 65536
        with self.assertRaisesRegex(ValueError, r"\[0, 65535\]"):
            synthetic.validate_cohort(c)

    def test_single_class_cohort_is_refused(self):
        base = {"n_samples": 8, "n_lineages": 2, "portable_signal": {"frequency": 1e-9, "logit": 20},
                "group_proxy": {"enabled": False}}
        for noise in (True, False):
            with self.subTest(noise=noise):
                with self.assertRaisesRegex(ValueError, "y = 0; a single-class"):
                    synthetic.generate({**base, "noise": {"enabled": noise}}, 0)
        # All R.
        all_r = {**base, "portable_signal": {"frequency": 1 - 1e-9, "logit": 20}}
        with self.assertRaisesRegex(ValueError, "y = 1; a single-class"):
            synthetic.generate(all_r, 0)
        # Two classes before the duplicate copy, one after: the refusal reads
        # y after the copy.
        dup = {"n_samples": 8, "n_lineages": 2, "group_proxy": {"enabled": False},
               "portable_signal": {"frequency": 1e-9, "logit": 3.0}, "duplicates": {"fraction": 1.0}}
        _, _, before = synthetic._phenotype(synthetic.resolve_config(dup), 7)
        self.assertTrue(before.any() and not before.all(), before)
        with self.assertRaisesRegex(ValueError, "y = 0; a single-class"):
            synthetic.generate(dup, 7)

    def test_fisher_exact_matches_exact_binomials(self):
        self.assertAlmostEqual(synthetic.fisher_exact_p(1, 9, 0, 390), 0.025, places=12)
        self.assertAlmostEqual(synthetic.fisher_exact_p(1, 19, 0, 380), 0.05, places=12)
        self.assertEqual(synthetic.fisher_exact_p(0, 5, 0, 7), 1.0)
        rng = np.random.default_rng(11)
        for _ in range(300):
            a, b, c, d = (int(v) for v in rng.integers(0, 25, size=4))
            if a + b == 0 or c + d == 0:
                continue
            exact = brute_fisher(a, b, c, d)
            got = synthetic.fisher_exact_p(a, b, c, d)
            self.assertLessEqual(abs(got - exact), 1e-12 * max(exact, 1e-300), (a, b, c, d))
            self.assertLessEqual(abs(got - synthetic.fisher_exact_p(c, d, a, b)), 1e-12 * max(got, 1e-300))
        # Large tables stay finite and agree with the exact value.
        self.assertLessEqual(abs(synthetic.fisher_exact_p(30, 970, 60, 940) - brute_fisher(30, 970, 60, 940)),
                             1e-9 * brute_fisher(30, 970, 60, 940))

    def test_fisher_tolerance_is_1e7_relative(self):
        """The re-check's near-tie tables: inside, a neighbour more likely by
        just under 1e-7 is kept; outside, one just over 1e-7 is not, though a
        looser tolerance would keep it."""
        one_e3 = fractions.Fraction(1, 10 ** 3)
        inside, outside = (3483, 3409, 2973, 2908), (3745, 3909, 2505, 2613)
        self.assertNotEqual(brute_fisher(*inside), brute_fisher(*inside, tol=0))
        self.assertNotEqual(brute_fisher(*outside), brute_fisher(*outside, tol=one_e3))
        for table in (inside, outside):
            exact = brute_fisher(*table)
            self.assertLessEqual(abs(synthetic.fisher_exact_p(*table) - exact), 1e-12 * exact, table)

    def test_fisher_balanced_large_tables_do_not_overflow(self):
        """Balanced tables with N >= 2000, where probabilities built from the
        low end of the support overflow."""
        for table in ((400, 600, 600, 400), (1100, 900, 900, 1100), (2300, 2700, 2700, 2300)):
            exact = brute_fisher(*table)
            got = synthetic.fisher_exact_p(*table)
            self.assertLess(exact, 1e-9, table)
            self.assertLessEqual(abs(got - exact), 1e-12 * exact, (table, got, exact))
        self.assertAlmostEqual(synthetic.fisher_exact_p(400, 600, 600, 400) / 4.3e-19, 1, delta=0.05)

    def test_rate_diff_p_counts_a_duplicate_group_once_per_class(self):
        p = synthetic._rate_diff_p
        # Six R samples, the first four in one group, all missing; ten S samples complete.
        y = [True] * 6 + [False] * 10
        miss = [True] * 6 + [False] * 10
        grouped = ["D001"] * 4 + [None] * 12
        self.assertAlmostEqual(p(miss, y, grouped), float(f"{brute_fisher(3, 0, 0, 10):.6g}"), places=12)
        self.assertAlmostEqual(p(miss, y, [None] * 16), float(f"{brute_fisher(6, 0, 0, 10):.6g}"), places=12)
        # A group spanning both classes is one unit in each class.
        mixed = ["D001", None, None, None, None, None, "D001"] + [None] * 9
        self.assertAlmostEqual(p(miss, y, mixed), float(f"{brute_fisher(6, 0, 0, 10):.6g}"), places=12)
        # A unit is missing when any member is: one missing and one complete R
        # member make one missing R unit. Two R units, five complete S units.
        part = p([True] + [False] * 7, [True] * 3 + [False] * 5, ["D001", "D001"] + [None] * 6)
        self.assertEqual(part, float(f"{brute_fisher(1, 1, 0, 5):.6g}"))
        self.assertEqual(part, 0.285714)
        # Null cases: no missing unit, every unit missing, an empty class.
        self.assertIsNone(p([False] * 16, y, [None] * 16))
        self.assertIsNone(p([True] * 16, y, [None] * 16))
        self.assertIsNone(p(miss, [True] * 16, [None] * 16))

    def test_flag_rule_boundaries(self):
        flag = synthetic.missingness_flag
        self.assertIs(flag(0.60, None), False)
        self.assertIs(flag(0.600001, None), True)
        self.assertIs(flag(0.5, 0.01), False)
        self.assertIs(flag(0.5, 0.009999), True)
        self.assertIs(flag(None, None), False)
        self.assertIs(flag(0.60, 0.01), False)

    def test_canonical_json_refuses_tuples(self):
        for bad in ((1, 2), {"a": ({1: "x"},)}, [[({"k": 1},)]]):
            with self.assertRaisesRegex(ValueError, "tuple"):
                synthetic.to_json(bad)

    def test_group_proxy_truth_reports_marker_fidelity_and_reference_scores(self):
        for seed in REFERENCE_SEEDS:
            r = effect(cohort(seed=seed), "group_proxy")["realized"]
            self.assertGreater(r["marker_lineage_auroc_mean"], 0.95, r)
            self.assertGreater(r["marker_lineage_auroc_min"], 0.9, r)
            self.assertGreater(r["marker_agreement"], 0.95, r)
            # Scoring by lineage prevalence is the best lineage-constant score.
            self.assertGreaterEqual(r["lineage_prevalence_auroc"], r["lineage_auroc"], r)
        c = cohort({"group_proxy": {"marker_flip": 0.5}})
        e = effect(c, "group_proxy")
        r = e["realized"]
        # Recompute each marker's AUROC for its own lineage from the JSON.
        lin = lineage_index(c)
        per_marker = [auroc(column(c, name)["raw"]["Binary"]["present"], lin == k)
                      for k, (_, names) in enumerate(sorted(e["markers_by_lineage"].items())) for name in names]
        self.assertAlmostEqual(r["marker_lineage_auroc_min"], min(per_marker), places=6)
        self.assertAlmostEqual(r["marker_lineage_auroc_mean"], float(np.mean(per_marker)), places=6)
        self.assertLess(r["marker_lineage_auroc_min"], r["marker_lineage_auroc_mean"] - 0.03, r)
        self.assertLess(abs(r["marker_lineage_auroc_mean"] - 0.5), 0.05, r)
        self.assertLess(abs(r["marker_agreement"] - 0.5), 0.05, r)
        y = labels(c)
        prev = {lid: v for lid, v in r["prevalence_by_lineage"].items()}
        score = np.array([prev[s["lineage"]] for s in c["samples"]])
        self.assertAlmostEqual(r["lineage_prevalence_auroc"], auroc(score, y), places=5)


class PlantedTruths(unittest.TestCase):
    """Reference checks with ridge logistic regression; each metric is the mean
    over REFERENCE_SEEDS. The lineage split is leave-one-lineage-out, scored by
    AUROC within the held-out lineage; the random split is stratified 4-fold."""

    @classmethod
    def setUpClass(cls):
        rows = []
        for seed in REFERENCE_SEEDS:
            c = cohort(seed=seed)
            y = labels(c)
            lin, rnd = lineage_index(c), random_folds(y)
            XC, XK, XZ = matrix(c, "C"), matrix(c, "K"), matrix(c, "noise")
            XM = matrix(c, "P", missing=True)
            rows.append({
                "C_lineage": cv_auroc(XC, y, lin), "C_random": cv_auroc(XC, y, rnd),
                "C_random_within": within_lineage_auroc(XC, y, rnd, lin),
                "K_lineage": cv_auroc(XK, y, lin), "K_random": cv_auroc(XK, y, rnd),
                "noise_lineage": cv_auroc(XZ, y, lin), "noise_random": cv_auroc(XZ, y, rnd),
                "missing_P_random": cv_auroc(XM, y, rnd),
                "missing_P_count": auroc(XM.sum(axis=1), y),
                **{f"missing_{b}_count": auroc(matrix(c, b, missing=True).sum(axis=1), y) for b in "CHN"},
            })
        cls.rows = rows
        cls.m = {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}

    def test_portable_signal_holds_under_lineage_holdout(self):
        m = self.m
        self.assertGreater(m["C_lineage"], 0.72, m)
        # Two-sided, like for like: within-lineage AUROC under lineage holdout
        # stays within 0.05 of the same estimand on the random split, in
        # every seed (seeds 0-199: max |difference| 0.0399 at seed 71).
        for r in self.rows:
            self.assertLess(abs(r["C_lineage"] - r["C_random_within"]), 0.05, r)
        self.assertGreater(m["C_lineage"] - m["K_lineage"], 0.15, m)
        # Equal across lineages: C:causal alone predicts y in every lineage.
        for seed in REFERENCE_SEEDS:
            c = cohort(seed=seed)
            by_lineage = effect(c, "portable_signal")["realized"]["auroc_by_lineage"]
            self.assertGreater(min(by_lineage.values()), 0.6, (seed, by_lineage))
            col = column(c, "C:causal")["raw"]["Binary"]
            y, lin = labels(c), lineage_index(c)
            ok = np.array(col["assessable"])
            present = np.array(col["present"], dtype=float)
            for k, lid in enumerate(sorted(by_lineage)):
                sel = ok & (lin == k)
                self.assertAlmostEqual(auroc(present[sel], y[sel]), by_lineage[lid], places=6)

    def test_group_proxy_predicts_on_random_split_and_falls_away_under_lineage(self):
        m = self.m
        self.assertGreater(m["K_random"], 0.68, m)
        self.assertLess(abs(m["K_lineage"] - 0.5), 0.07, m)
        self.assertGreater(m["K_random"] - m["K_lineage"], 0.15, m)
        self.assertTrue(all(r["K_random"] - r["K_lineage"] > 0.08 for r in self.rows), self.rows)
        c = cohort()
        realized = effect(c, "group_proxy")["realized"]
        self.assertGreater(realized["lineage_auroc"], 0.65)
        y = labels(c)
        off = cohort({"group_proxy": {"enabled": False}})
        self.assertLess(abs(cv_auroc(matrix(off, "K"), labels(off), random_folds(labels(off))) - 0.5), 0.08)
        self.assertGreater(cv_auroc(matrix(c, "K"), y, random_folds(y)), 0.65)

    def test_noise_block_is_near_the_null(self):
        m = self.m
        self.assertLess(abs(m["noise_random"] - 0.5), 0.06, m)
        self.assertLess(abs(m["noise_lineage"] - 0.5), 0.06, m)
        for seed in REFERENCE_SEEDS:
            c = cohort(seed=seed)
            self.assertLess(effect(c, "noise")["realized"]["max_abs_auroc_minus_half"], 0.12, seed)

    def test_missingness_indicators_alone_predict_y(self):
        for r in self.rows:
            self.assertGreater(r["missing_P_random"], 0.6, r)
            self.assertGreater(max(r["missing_P_count"], 1 - r["missing_P_count"]), 0.6, r)
            for b in "CHN":
                self.assertLess(abs(r[f"missing_{b}_count"] - 0.5), 0.08, (b, r))
        c = cohort()
        realized = effect(c, "missingness_shortcut")["realized"]
        self.assertLess(abs(realized["rate_R"] - 0.05), 0.03, realized)
        self.assertLess(abs(realized["rate_S"] - 0.35), 0.05, realized)
        self.assertGreater(realized["indicator_auroc_oriented"], 0.6)
        self.assertTrue(realized["flag"])
        self.assertAlmostEqual(realized["indicator_auroc"],
                               auroc(matrix(c, "P", missing=True).sum(axis=1), labels(c)), places=6)
        self.assertEqual(realized, {k: v for k, v in c["missingness"][2].items() if k != "block"})
        off = cohort({"missingness_shortcut": {"enabled": False}})
        self.assertLess(abs(auroc(matrix(off, "P", missing=True).sum(axis=1), labels(off)) - 0.5), 0.08)

    def test_missingness_can_target_another_block(self):
        c = cohort({"missingness_shortcut": {"block": "N"}})
        y = labels(c)
        self.assertGreater(cv_auroc(matrix(c, "N", missing=True), y, random_folds(y)), 0.6)
        self.assertLess(abs(auroc(matrix(c, "P", missing=True).sum(axis=1), y) - 0.5), 0.08)
        self.assertEqual(effect(c, "missingness_shortcut")["block"], "N")


if __name__ == "__main__":
    unittest.main()
