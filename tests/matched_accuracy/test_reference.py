"""The committed texts and pins of stage 1: the protocol, the amendment and the values the rulings and brief give."""

import hashlib
import unittest

from tests.matched_accuracy import common

s1 = common.stage1
# The sha256 sums the lead's rulings R49 (the protocol) and R51 (the amendment) name.
PROTOCOL_SHA256 = "e4c124466f7b85c649ea4684004831beb80002174e51317ad0ee7476a84401b0"
AMENDMENT_SHA256 = "634ffb8db8587a34eba879465a414a518d4d1afdaf730bca69b87b5c782fe27e"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Reference(unittest.TestCase):
    def test_the_protocol_and_the_amendment_are_the_bytes_the_rulings_name(self):
        self.assertEqual(sha(common.REFERENCE / "PROTOCOL-1.0.md"), PROTOCOL_SHA256)
        self.assertEqual(sha(common.REFERENCE / "AMENDMENT-1.md"), AMENDMENT_SHA256)
        self.assertEqual(sorted(p.name for p in common.REFERENCE.iterdir()), ["AMENDMENT-1.md", "PROTOCOL-1.0.md"])

    def test_the_pins_are_the_values_of_the_brief_and_rulings(self):
        self.assertEqual(s1.RECORDS_SHA256, "d231fa768a31fc4ce4c8a6bb7024d05d632621846070688adb892fbfe340fd78")
        self.assertEqual(s1.RECORDS_LINES, 1285111)
        self.assertEqual(s1.DROPPED_RECORDED, 252148)
        self.assertEqual(s1.GENOMES_SHA256, "84ae57a777cb234aa0a39102ab7f9ca9448e52d1474654179f26165644c25f64")
        self.assertEqual(s1.SEAL, "2026-11-05T19:10:00")
        self.assertEqual(s1.FIVE, ("scm", "cart", "rf", "gbt", "lasso"))
        self.assertEqual(s1.FOLDS, 5)

    def test_the_sums_list_of_the_release_is_the_committed_one(self):
        self.assertEqual(s1.REPORT_SUMS, common.ROOT / "reference" / "lineage-control" / "rev14-274-20261112.reports.sha256")
        self.assertEqual(sha(s1.REPORT_SUMS), s1.REPORT_SUMS_SHA256)
        lines = s1.REPORT_SUMS.read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(lines), 274)
        self.assertTrue(all(line.split("  ")[1].startswith("./datasets/") for line in lines))


if __name__ == "__main__":
    unittest.main()
