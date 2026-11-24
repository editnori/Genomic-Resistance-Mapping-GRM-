"""scripts/c5-settings/page_faults.py: each fault planted in the table script and the run's check script is caught
by its test, and each fault's text occurs once in its file. tests/c5_settings/test_faults.py checks that the block's
own planted_faults.py still finds each of its texts once.

The faults are planted in memory only; every input is SYNTHETIC.
"""

import unittest

from tests.c5_settings import common

faults = common.load("page_faults")
page = common.load("page_tables")


class Planted(unittest.TestCase):
    def test_each_fault_text_occurs_once_in_its_file(self):
        for n, _, path, old, new, _, _ in faults.FAULTS:
            with self.subTest(fault=n):
                text = (common.ROOT / path).read_text(encoding="utf-8")
                self.assertEqual(text.count(old), 1, path)
                self.assertNotEqual(old, new)

    def test_every_gate_of_the_table_script_has_a_fault(self):
        gates = {f[5] for f in faults.FAULTS if f[2] == faults.PAGE and f[5] != "table"}
        self.assertEqual(gates, {f"{n} {name}" for n, name, _ in page.GATES})
        self.assertEqual(len({f[0] for f in faults.FAULTS}), len(faults.FAULTS))

    def test_every_fault_is_caught_by_its_test(self):
        rows = faults.run_all()
        self.assertEqual(len(rows), 48)
        for r in rows:
            with self.subTest(fault=r["fault"]):
                self.assertEqual((r["as_written"], r["with_fault"], r["caught"]), ("pass", "fail", "yes"), r["planted"])


if __name__ == "__main__":
    unittest.main()
