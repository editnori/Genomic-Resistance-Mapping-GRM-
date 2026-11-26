"""scripts/c7-outside/page_faults.py: each fault's text occurs once in its file, its test exists, and the faults
cover every gate; the block's own planted_faults.py still finds each of its texts once; and no test name holds
another, since a runner's -k is a substring.

A runner stops at the first fault whose text it cannot find once, so these tests hold every fault to the tree as it
stands. They plant nothing and run no fault; page_faults.py and planted_faults.py do that.
"""

import ast
import unittest

from tests.c7_outside import common

faults = common.load("page_faults")
planted = common.load("planted_faults")
page = common.load("page_tables")
ROOT = common.ROOT
RUN = ("run: species", "run: SHA256SUMS", "run: DONE", "run: expect.py", "run: genome ids", "run: features",
       "run: score.log")


def tests_of(module):
    """The class.test names of one test module under tests/c7_outside, read from its source."""
    tree = ast.parse((ROOT / "tests" / "c7_outside" / f"{module}.py").read_text())
    return {f"{c.name}.{f.name}" for c in tree.body if isinstance(c, ast.ClassDef)
            for f in c.body if isinstance(f, ast.FunctionDef) and f.name.startswith("test_")}


class PageFaults(unittest.TestCase):
    def test_each_page_faults_text_occurs_once_in_its_file(self):
        for n, _, _, path, old, new, _ in faults.FAULTS:
            text = (ROOT / path).read_text()
            self.assertEqual(text.count(old), 1, f"fault {n}: {path}")
            self.assertNotEqual(old, new, f"fault {n}")

    def test_each_page_faults_test_exists(self):
        for n, *_, test in faults.FAULTS:
            parts = test.split(".")
            self.assertEqual(parts[:2], ["tests", "c7_outside"], f"fault {n}")
            self.assertIn(".".join(parts[3:]), tests_of(parts[2]), f"fault {n}")

    def test_the_page_faults_are_named_once_and_cover_every_gate(self):
        names = [f[0] for f in faults.FAULTS]
        self.assertEqual(len(names), len(set(names)))
        gates = {f[1] for f in faults.FAULTS}
        self.assertLessEqual({name for _, name, _ in page.GATES} | set(RUN), gates)
        self.assertGreaterEqual(len(faults.FAULTS), 8)

    def test_each_page_fault_plants_in_the_table_script_or_the_run_checks(self):
        self.assertEqual({f[3] for f in faults.FAULTS}, {faults.PAGE, faults.CHECKS})


class BlockFaults(unittest.TestCase):
    def test_planted_faults_finds_each_of_its_texts_once(self):
        self.assertEqual(len(planted.FAULTS), 13)
        for n, _, path, old, new, _, _ in planted.FAULTS:
            self.assertEqual((ROOT / path).read_text().count(old), 1, f"fault {n}: {path}")
            self.assertNotEqual(old, new, f"fault {n}")


class TestNames(unittest.TestCase):
    def test_no_test_id_in_the_suite_holds_another(self):
        ids = []
        for path in sorted((ROOT / "tests" / "c7_outside").glob("test_*.py")):
            ids += [f"tests.c7_outside.{path.stem}.{name}" for name in tests_of(path.stem)]
        for a in ids:
            held = [b for b in ids if a != b and a in b]
            self.assertEqual(held, [], a)


if __name__ == "__main__":
    unittest.main()
