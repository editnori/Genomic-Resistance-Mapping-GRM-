"""The planted faults of scripts/time-split/stagger_faults.py: their count, each fault's text once in its file, a test
of this suite for each that `-k` picks out alone, and the stop on a barred out dir before any directory is made.
Reads the runner and the files it names; plants nothing."""

import importlib
import os
import re
import shutil
import tempfile
import unittest
import unittest.mock
from pathlib import Path

from tests.time_split import stagger_fixture as fx

runner = fx.load(fx.FAULTS, "stagger_faults_under_test")
COUNT = 103


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
                self.assertIn(path, (runner.DESIGN, runner.ANALYZE, runner.RUNNER))
                self.assertEqual((fx.ROOT / path).read_text(encoding="utf-8").count(old), 1, path)

    def test_each_fault_names_one_test_of_this_suite(self):
        ids = []
        for module in ("tests.time_split.test_stagger_design", "tests.time_split.test_stagger_analyze",
                       "tests.time_split.test_stagger_faults"):
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

    def test_the_page_gives_the_count(self):
        self.assertEqual(re.findall(r"(\d+) planted faults", fx.PAGE.read_text(encoding="utf-8")), [str(COUNT)])


class OutDir(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="stagger-faults-"))

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_a_barred_out_dir_stops_before_any_directory_is_made(self):
        sealed = self.tmp / "sealed-place"
        sealed.mkdir()
        os.symlink(sealed, self.tmp / "plain-link")
        cases = [(self.tmp / "Sealed-set" / "out", "sealed"), (self.tmp / "held-out" / "out", "held-out"),
                 (self.tmp / "LATETEST" / "out", "latetest"), (self.tmp / "plain-link" / "out", "sealed")]
        before = sorted(p.name for p in self.tmp.rglob("*"))
        # No git archive and no test run: a runner that reaches them has made its directories already.
        with unittest.mock.patch.object(runner.subprocess, "run", side_effect=AssertionError("a command ran")):
            for out, word in cases:
                with self.subTest(out=out.parent.name):
                    code, text, err = fx.run(runner, out, "cut-floor")
                    self.assertEqual((code, text), (1, ""))
                    self.assertEqual(err, f"stagger_faults: stopped: the out dir path holds '{word}'; nothing under "
                                          "it is opened\n")
                    self.assertFalse(out.exists())
        self.assertEqual(sorted(p.name for p in self.tmp.rglob("*")), before)
        self.assertEqual(list(sealed.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
