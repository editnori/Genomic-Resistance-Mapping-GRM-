"""scripts/matched-accuracy/planted_faults.py: each fault's text occurs once in its file, plants code that compiles,
and names one test of this suite that runs alone. These tests plant nothing and run no fault; the runner does."""

import ast
import importlib.util
import tempfile
import unittest
from pathlib import Path

from tests.matched_accuracy import common

ROOT = common.ROOT
# The rules of stage 1, each with at least one fault. The checks that stop the run are the first six.
RULES = ("seal", "line count", "folds", "lineage cluster", "copy group", "lab call", "loader", "learners", "reasons",
         "BioProject", "drug names", "tuple", "intermediate", "groups", "pairs", "reading B", "no group", "weights",
         "shares", "typing", "row list", "scan", "export", "dropped figure", "sha256", "input shape", "write")


def load_faults():
    spec = importlib.util.spec_from_file_location("matched_accuracy_faults", common.FAULTS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


faults = load_faults()


def test_ids():
    """Every test id of this suite, read from the sources: tests.matched_accuracy.<module>.<class>.<test>."""
    ids = []
    for path in sorted((ROOT / common.SUITE).glob("test_*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for c in tree.body:
            if isinstance(c, ast.ClassDef):
                ids += [f"tests.matched_accuracy.{path.stem}.{c.name}.{f.name}" for f in c.body
                        if isinstance(f, ast.FunctionDef) and f.name.startswith("test_")]
    return ids


class Faults(unittest.TestCase):
    def test_each_faults_text_occurs_once_in_its_file(self):
        for n, _, _, path, old, new, _ in faults.FAULTS:
            text = (ROOT / path).read_text(encoding="utf-8")
            self.assertEqual(text.count(old), 1, f"fault {n}: {path}")
            self.assertNotEqual(old, new, f"fault {n}")

    def test_each_planted_file_still_compiles(self):
        # A fault that broke the file's syntax would fail its test for that alone.
        for n, _, _, path, old, new, _ in faults.FAULTS:
            compile((ROOT / path).read_text(encoding="utf-8").replace(old, new), path, "exec")

    def test_each_faults_test_is_one_test_of_this_suite(self):
        ids = test_ids()
        for n, *_, test in faults.FAULTS:
            self.assertIn(test, ids, f"fault {n}")
            # unittest -k matches a substring of the test's id: the runner must select this one test alone.
            self.assertEqual([i for i in ids if test in i], [test], f"fault {n}")

    def test_the_faults_are_named_once_and_cover_every_rule(self):
        names = [f[0] for f in faults.FAULTS]
        self.assertEqual(names, [f"P{i:02d}" for i in range(1, len(names) + 1)])
        self.assertEqual({f[1] for f in faults.FAULTS}, set(RULES))
        self.assertEqual({f[3] for f in faults.FAULTS}, {"scripts/matched-accuracy/stage1.py"})

    def test_a_fault_counts_only_when_its_test_fails_by_an_assertion(self):
        cases = (("FAILED (failures=1)", (1, 0)), ("FAILED (errors=1)", (0, 1)),
                 ("FAILED (failures=2, errors=1)", (2, 1)), ("OK", (0, 0)), ("Ran 1 test", None))
        for last, want in cases:
            with tempfile.TemporaryDirectory() as d:
                log = Path(d) / "x.log"
                log.write_text(f"test_x ... ok\n\n{'-' * 70}\nRan 1 test in 0.1s\n\n{last}\n")
                self.assertEqual(faults.verdict(log), want, last)

    def test_the_runner_copies_the_tree_and_runs_each_test_as_ci_does(self):
        source = common.FAULTS.read_text(encoding="utf-8")
        self.assertIn('for top in ("scripts", "tests", "reference", "docs"):', source)
        self.assertIn("run_test = _c9().run_test", source)
        self.assertIn('result = "failed by an assertion"', source)


if __name__ == "__main__":
    unittest.main()
