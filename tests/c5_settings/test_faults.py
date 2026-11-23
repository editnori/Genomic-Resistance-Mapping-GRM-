"""Each planted fault of scripts/c5-settings/planted_faults.py still applies to this tree. The runner stops at the
first fault whose text to change does not occur exactly once in its file, so a line another block moves (C4 routes
the forest's screen through src/learn/screen_arm.rs) stops it. Reads the committed runner and the files it names;
reads no data and plants nothing.
"""

import unittest

from tests.c5_settings import common


class Faults(unittest.TestCase):
    def test_each_faults_text_occurs_once_in_its_file(self):
        faults = common.load("planted_faults").FAULTS
        self.assertTrue(faults)
        for fault in faults:
            n, path, old = fault[0], fault[2], fault[3]
            with self.subTest(fault=n):
                self.assertEqual((common.ROOT / path).read_text(encoding="utf-8").count(old), 1, path)


if __name__ == "__main__":
    unittest.main()
