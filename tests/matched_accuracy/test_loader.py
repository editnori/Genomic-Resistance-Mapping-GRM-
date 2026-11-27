"""The report loader of scripts/matched-accuracy/stage1.py: the named keys alone, and no score reachable.

Each report here is SYNTHETIC (tests/matched_accuracy/made_up.py): its scores are made-up numbers, and its roc,
metrics, conformal and no_call values carry made-up marker texts.
"""

import ast
import json
import tempfile
import unittest
from pathlib import Path

from tests.matched_accuracy import common, made_up

s1 = common.stage1


def reachable(obj, seen=None):
    """Every object reachable from obj through dict keys and values, list, tuple and set items and attributes."""
    seen = set() if seen is None else seen
    if id(obj) in seen:
        return
    seen.add(id(obj))
    yield obj
    items = []
    if isinstance(obj, dict):
        items = list(obj.keys()) + list(obj.values())
    elif isinstance(obj, (list, tuple, set, frozenset)):
        items = list(obj)
    elif hasattr(obj, "__dict__"):
        items = list(vars(obj).values())
    for x in items:
        yield from reachable(x, seen)


def loaded(doc):
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "report.json"
        p.write_text(json.dumps(doc), encoding="utf-8")
        return s1.load_report(p)


class Loader(unittest.TestCase):
    def setUp(self):
        rel = common.release()
        self.doc = common.report_of(rel, "synthia_alpha__drug_a")
        self.rep = loaded(self.doc)

    def test_the_loader_returns_the_named_keys_alone(self):
        self.assertEqual(set(self.rep), set(s1.REPORT_KEYS) | {"samples", "lineage_split"})
        for s in self.rep["samples"]:
            self.assertEqual(set(s), set(s1.ROW_KEYS))
        self.assertEqual(set(self.rep["lineage_split"]), set(s1.SPLIT_KEYS) | {"learners"})
        for e in self.rep["lineage_split"]["learners"]:
            self.assertEqual(set(e), {"learner", "scores_length"})

    def test_no_score_or_value_of_roc_metrics_conformal_or_no_call_is_reachable(self):
        scores = {x for sp in self.doc["splits"] for e in sp["learners"] for x in e["scores"]}
        self.assertTrue(scores)
        for obj in reachable(self.rep):
            self.assertNotIsInstance(obj, float)
            if isinstance(obj, str):
                for word in made_up.SENTINELS:
                    self.assertNotIn(word, obj)
                self.assertNotIn(obj, ("scores", "roc", "metrics", "conformal", "no_call", "random"))
            if isinstance(obj, (int, float)) and not isinstance(obj, bool):
                self.assertNotIn(obj, scores)

    def test_the_loader_keeps_each_scores_length_alone(self):
        n = len(self.doc["samples"])
        got = [(e["learner"], e["scores_length"]) for e in self.rep["lineage_split"]["learners"]]
        self.assertEqual(got, [(x, n) for x in made_up.FIVE + ("catalog",)])

    def test_the_loader_keeps_the_lineage_split_and_drops_the_random_one(self):
        lineage = next(s for s in self.doc["splits"] if s["scheme"] == "lineage")
        self.assertEqual(self.rep["lineage_split"]["test_fold"], lineage["test_fold"])
        self.assertEqual(self.rep["lineage_split"]["folds"], 5)
        self.assertIsNone(self.rep["lineage_split"]["unavailable_reason"])

    def test_a_report_without_learners_keeps_its_unavailable_reason(self):
        doc = common.report_of(common.release(), "synthia_delta__drug_c")
        rep = loaded(doc)
        self.assertEqual(rep["lineage_split"]["learners"], [])
        self.assertIsNone(rep["lineage_split"]["test_fold"])
        self.assertTrue(rep["lineage_split"]["unavailable_reason"])
        self.assertFalse(s1.enters(rep))


class Source(unittest.TestCase):
    """The script's text: no code path reads, keeps, sums or compares a score."""

    @classmethod
    def setUpClass(cls):
        cls.source = common.SCRIPT.read_text(encoding="utf-8")
        cls.tree = ast.parse(cls.source)

    def test_the_word_scores_stands_in_code_once_inside_a_length(self):
        strings = [n for n in ast.walk(self.tree) if isinstance(n, ast.Constant) and n.value == "scores"]
        self.assertEqual(len(strings), 2)   # e.get("scores") and e["scores"], both in the one length
        line = self.source.split("\n")[strings[0].lineno - 1]
        self.assertEqual(strings[0].lineno, strings[1].lineno)
        self.assertIn('len(e["scores"]) if isinstance(e.get("scores"), list)', line)

    def test_the_script_names_no_value_of_a_learner_but_its_name_and_scores(self):
        for word in ("roc", "metrics", "conformal", "no_call", "fold_balanced_accuracy", "fold_complexity",
                     "catalog_resistant", "random"):
            consts = [n for n in ast.walk(self.tree) if isinstance(n, ast.Constant) and n.value == word]
            self.assertEqual(consts, [], word)

    def test_the_script_uses_the_standard_library_alone(self):
        imported = {a.name.split(".")[0] for n in ast.walk(self.tree) if isinstance(n, ast.Import) for a in n.names}
        imported |= {n.module.split(".")[0] for n in ast.walk(self.tree) if isinstance(n, ast.ImportFrom)}
        self.assertLessEqual(imported, {"argparse", "collections", "datetime", "fractions", "gzip", "hashlib", "json",
                                        "os", "re", "sys", "pathlib"})

    def test_the_script_writes_no_stage_2(self):
        for word in ("numpy", "PCG64", "percentile", "bootstrap", "resample", "statistics", "mean("):
            self.assertNotIn(word, self.source)


if __name__ == "__main__":
    unittest.main()
