"""Each planted fault of scripts/c2-cap/planted_faults.py and scripts/c2-cap/page_faults.py still applies to this
tree. A runner stops at the first fault whose text to change does not occur exactly once in its file, so a line
another block adds or moves stops it; and a runner's -k is a substring, so each fault's test name must select its
one test. Reads the committed runners and the files they name; reads no data and plants nothing.
"""

import unittest

from tests.c2_cap import common


class Faults(unittest.TestCase):
    def test_each_planted_faults_text_occurs_once_in_its_file(self):
        faults = common.load("planted_faults").FAULTS
        self.assertEqual(len(faults), 18)
        for fault in faults:
            n, path, old = fault[0], fault[2], fault[3]
            with self.subTest(fault=n):
                self.assertEqual((common.ROOT / path).read_text(encoding="utf-8").count(old), 1, path)

    def test_each_page_faults_text_occurs_once_in_its_file(self):
        faults = common.load("page_faults").FAULTS
        self.assertEqual(len(faults), 19)
        self.assertEqual(len({f[0] for f in faults}), len(faults))
        for fault in faults:
            n, path, old = fault[0], fault[2], fault[3]
            with self.subTest(fault=n):
                self.assertEqual((common.ROOT / path).read_text(encoding="utf-8").count(old), 1, path)

    def test_each_fault_runners_test_name_selects_one_test(self):
        loader = unittest.TestLoader()
        ids = []
        for suite in loader.discover(str(common.ROOT / "tests" / "c2_cap"), top_level_dir=str(common.ROOT)):
            stack = [suite]
            while stack:
                s = stack.pop()
                if isinstance(s, unittest.TestSuite):
                    stack.extend(s)
                else:
                    ids.append(s.id())
        for runner in ("planted_faults", "page_faults"):
            for fault in common.load(runner).FAULTS:
                with self.subTest(runner=runner, fault=fault[0]):
                    self.assertEqual(sum(fault[6] in i for i in ids), 1, fault[6])


if __name__ == "__main__":
    unittest.main()
