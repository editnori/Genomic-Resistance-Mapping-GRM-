"""scripts/c3-inputs/page_faults.py: each fault's text occurs once in its file, and its test exists.

A runner stops at the first fault whose text it cannot find once, so these
tests hold every fault to the tree as it stands. They plant nothing and run
no fault; page_faults.py does that.
"""

import ast
import unittest

from tests.c3_inputs import common

faults = common.load("page_faults")
ROOT = common.ROOT
# The gates the page's brief names, each with a fault.
GATES = ("species", "M", "own datasets", "twin", "study", "constant", "draws", "conditions", "binary")
RUN = ("run: SHA256SUMS", "run: DONE and binary.txt", "run: expect.py", "run: refusals", "run: genome ids")


def tests_of(module):
    """The class.test names of one test module under tests/c3_inputs, read from its source."""
    tree = ast.parse((ROOT / "tests" / "c3_inputs" / f"{module}.py").read_text())
    return {f"{c.name}.{f.name}" for c in tree.body if isinstance(c, ast.ClassDef)
            for f in c.body if isinstance(f, ast.FunctionDef) and f.name.startswith("test_")}


class Faults(unittest.TestCase):
    def test_each_faults_text_occurs_once_in_its_file(self):
        for n, _, _, path, old, new, _ in faults.FAULTS:
            text = (ROOT / path).read_text()
            self.assertEqual(text.count(old), 1, f"fault {n}: {path}")
            self.assertNotEqual(old, new, f"fault {n}")

    def test_each_faults_test_exists(self):
        for n, *_, test in faults.FAULTS:
            parts = test.split(".")
            self.assertEqual(parts[:2], ["tests", "c3_inputs"], f"fault {n}")
            self.assertIn(".".join(parts[3:]), tests_of(parts[2]), f"fault {n}")

    def test_the_faults_are_named_once_and_cover_every_gate(self):
        names = [f[0] for f in faults.FAULTS]
        self.assertEqual(len(names), len(set(names)))
        gates = {f[1] for f in faults.FAULTS}
        self.assertLessEqual(set(GATES) | set(RUN), gates)
        self.assertGreaterEqual(len(faults.FAULTS), 8)

    def test_each_fault_plants_in_the_table_script_or_the_run_checks(self):
        self.assertEqual({f[3] for f in faults.FAULTS}, {faults.PAGE, faults.CHECKS})


if __name__ == "__main__":
    unittest.main()
