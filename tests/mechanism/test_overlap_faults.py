"""scripts/mechanism/overlap_faults.py: each fault planted in the learner overlap's table script or its second pass
is caught by its test.

The faults are planted in memory only; every input is SYNTHETIC.
"""

import unittest

from tests.mechanism import common as c

faults = c.load("overlap_faults")


class Planted(unittest.TestCase):
    def test_every_fault_is_caught_by_its_test(self):
        rows = faults.run_all()
        self.assertEqual([r["fault"] for r in rows], [f"O{i}" for i in range(1, 29)])
        for r in rows:
            with self.subTest(fault=r["fault"]):
                self.assertEqual((r["as_written"], r["with_fault"], r["caught"]), ("pass", "fail", "yes"), r["planted"])


if __name__ == "__main__":
    unittest.main()
