"""scripts/c3-inputs/m_faults.py: each fault's text occurs once in its file, plants code that runs, and has its test.

A runner stops at the first fault whose text it cannot find once, so these
tests hold every fault to the tree as it stands. They plant nothing and run
no fault; m_faults.py does that.
"""

import ast
import subprocess
import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import common

faults = common.load("m_faults")
ROOT = common.ROOT
# The new rules of the run on the datasets of M, each with a fault.
RULES = ("list", "floor", "check species", "Same build", "Same results", "Same counts", "second reading",
         "first run", "Holm", "counts", "tables", "pack", "cost", "run script")


def tests_of(module):
    """The class.test names of one test module under tests/c3_inputs, read from its source."""
    tree = ast.parse((ROOT / "tests" / "c3_inputs" / f"{module}.py").read_text())
    return {f"{c.name}.{f.name}" for c in tree.body if isinstance(c, ast.ClassDef)
            for f in c.body if isinstance(f, ast.FunctionDef) and f.name.startswith("test_")}


class MFaults(unittest.TestCase):
    def test_each_faults_text_occurs_once_in_its_file(self):
        for n, _, _, path, old, new, _ in faults.FAULTS:
            text = (ROOT / path).read_text()
            self.assertEqual(text.count(old), 1, f"fault {n}: {path}")
            self.assertNotEqual(old, new, f"fault {n}")

    def test_each_planted_file_still_parses(self):
        # A fault that broke the file's syntax would fail its test for that alone.
        for n, _, _, path, old, new, _ in faults.FAULTS:
            planted = (ROOT / path).read_text().replace(old, new)
            if path.endswith(".py"):
                compile(planted, path, "exec")
            else:
                with tempfile.TemporaryDirectory() as d:
                    p = Path(d) / "planted.sh"
                    p.write_text(planted)
                    self.assertEqual(subprocess.run(["bash", "-n", str(p)], capture_output=True).returncode, 0, n)

    def test_each_faults_test_exists_in_a_module_of_the_run_on_m(self):
        for n, *_, test in faults.FAULTS:
            parts = test.split(".")
            self.assertEqual(parts[:2], ["tests", "c3_inputs"], f"fault {n}")
            self.assertTrue(parts[2].startswith("test_m_"), f"fault {n}")
            self.assertIn(".".join(parts[3:]), tests_of(parts[2]), f"fault {n}")

    def test_the_faults_are_named_once_and_cover_every_rule_and_file(self):
        names = [f[0] for f in faults.FAULTS]
        self.assertEqual(len(names), len(set(names)))
        self.assertEqual({f[1] for f in faults.FAULTS}, set(RULES))
        self.assertEqual({f[3] for f in faults.FAULTS},
                         {faults.LIST, faults.CHECKS, faults.SECOND, faults.READING, faults.PACK, faults.COST, faults.RUN})

    def test_no_test_name_of_the_run_on_m_holds_another(self):
        # unittest -k matches a substring of the test's id: a runner selecting one test must run that one alone.
        ids = [f"tests.c3_inputs.{p.stem}.{t}" for p in sorted((ROOT / "tests" / "c3_inputs").glob("test_*.py"))
               for t in sorted(tests_of(p.stem))]
        for n, *_, test in faults.FAULTS:
            self.assertEqual([i for i in ids if test in i], [test], f"fault {n}")


if __name__ == "__main__":
    unittest.main()
