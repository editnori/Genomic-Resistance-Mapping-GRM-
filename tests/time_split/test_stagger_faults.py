"""The planted faults of scripts/time-split/stagger_faults.py: their count, each fault's text once in its file, and a
test of this suite for each that `-k` picks out alone. Reads the runner and the files it names; plants nothing."""

import importlib
import re
import unittest

from tests.time_split import stagger_fixture as fx

runner = fx.load(fx.FAULTS, "stagger_faults_under_test")
COUNT = 66


class Faults(unittest.TestCase):
    def test_the_count(self):
        self.assertEqual(len(runner.FAULTS), COUNT)

    def test_each_name_once(self):
        names = [f[0] for f in runner.FAULTS]
        self.assertEqual(len(names), len(set(names)))
        for n in names:
            self.assertRegex(n, r"^[a-z0-9-]+$")

    def test_each_faults_text_occurs_once_in_its_file(self):
        for n, _, path, old, new, _ in runner.FAULTS:
            with self.subTest(fault=n):
                self.assertNotEqual(old, new)
                self.assertIn(path, (runner.DESIGN, runner.ANALYZE))
                self.assertEqual((fx.ROOT / path).read_text(encoding="utf-8").count(old), 1, path)

    def test_each_fault_names_one_test_of_this_suite(self):
        ids = []
        for module in ("tests.time_split.test_stagger_design", "tests.time_split.test_stagger_analyze"):
            suite = unittest.TestLoader().loadTestsFromModule(importlib.import_module(module))
            stack = [suite]
            while stack:
                for t in stack.pop():
                    if isinstance(t, unittest.TestSuite):
                        stack.append(t)
                    else:
                        ids.append(t.id())
        for n, _, _, _, _, test in runner.FAULTS:
            with self.subTest(fault=n):
                self.assertEqual([i for i in ids if test in i], [test])

    def test_the_ci_page_gives_the_count(self):
        text = (fx.ROOT / "docs" / "ci.md").read_text(encoding="utf-8")
        self.assertEqual(re.findall(r"`scripts/time-split/stagger_faults\.py` lists (\d+) planted faults", text),
                         [str(COUNT)])


if __name__ == "__main__":
    unittest.main()
