"""scripts/mechanism/planted_faults.py: each planted fault is caught by its test (docs/mechanism-split.md).

The faults are planted in the scripts in memory only; every input is SYNTHETIC.
"""

import unittest

from tests.mechanism import common as c

faults = c.load("planted_faults")


class Planted(unittest.TestCase):
    def test_every_fault_is_caught_by_its_test(self):
        rows = faults.run_all()
        self.assertEqual([r["fault"] for r in rows], list(range(1, 12)))
        for r in rows:
            with self.subTest(fault=r["fault"]):
                self.assertEqual((r["as_written"], r["with_fault"], r["caught"]), ("pass", "fail", "yes"), r["planted"])


if __name__ == "__main__":
    unittest.main()
