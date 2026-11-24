"""scripts/mechanism/year_faults.py: each fault planted in the mechanism year read's script is caught by its test.

The faults are planted in memory only; every input is SYNTHETIC.
"""

import unittest

from tests.mechanism import common as c

faults = c.load("year_faults")


class Planted(unittest.TestCase):
    def test_every_fault_is_caught_by_its_test(self):
        rows = faults.run_all()
        self.assertEqual([r["fault"] for r in rows], [f"Y{i}" for i in range(1, 17)])
        for r in rows:
            with self.subTest(fault=r["fault"]):
                self.assertEqual((r["as_written"], r["with_fault"], r["caught"]), ("pass", "fail", "yes"), r["planted"])


if __name__ == "__main__":
    unittest.main()
