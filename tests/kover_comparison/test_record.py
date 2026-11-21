"""reference/kover-comparison: its hashes, columns and size, the derived files, and the file it is not."""

import hashlib
import json
import re
import unittest

from tests.kover_comparison.common import LEGACY, RECORD, load_script, read_tsv

record = load_script()
# reference/kover-summary.json as it stands on main; the record leaves it alone.
LEGACY_SHA256 = "4197fa5cca20bb21a68a891e02832fa53f3908b559e1a76f815a77842637e902"


def record_files():
    return sorted(p.name for p in RECORD.iterdir() if p.is_file())


class Record(unittest.TestCase):
    def test_sha256sums_whole(self):
        lines = (RECORD / "SHA256SUMS").read_text(encoding="utf-8").rstrip("\n").split("\n")
        listed = {}
        for line in lines:
            digest, name = line.split("  ", 1)
            listed[name] = digest
        self.assertEqual(sorted(listed), [n for n in record_files() if n != "SHA256SUMS"])
        for name, digest in listed.items():
            self.assertEqual(hashlib.sha256((RECORD / name).read_bytes()).hexdigest(), digest, name)

    def test_the_files_and_nothing_else(self):
        self.assertEqual(record_files(), sorted(record.ROW_FILES + record.DERIVED + ["README.md", "SHA256SUMS"]))
        self.assertEqual([p for p in RECORD.iterdir() if not p.is_file()], [])

    def test_under_3_mb(self):
        self.assertLess(sum((RECORD / n).stat().st_size for n in record_files()), 3 * 1024 * 1024)

    def test_columns(self):
        for name, cols in record.COLUMNS.items():
            self.assertEqual(read_tsv(RECORD / name, header_only=True), cols, name)
            for row in read_tsv(RECORD / name):
                self.assertEqual(len(row), len(cols), name)
                self.assertNotIn(None, row.values(), name)

    def test_record_is_what_the_script_writes(self):
        try:
            files = record.assemble(record.from_record(RECORD))
        except record.Stop as exc:
            self.fail("the script stops on the record: %s" % exc)
        self.assertEqual(record.check(files, RECORD), [])

    def test_datasets_file_is_the_one_the_plan_hashed(self):
        plan = [r for r in read_tsv(RECORD / "order.tsv") if r["event"] == "plan_hashed"]
        self.assertEqual([r["file"] for r in plan], record.PLAN_FILES)
        digest = hashlib.sha256((RECORD / "datasets.tsv").read_bytes()).hexdigest()
        self.assertEqual(plan[1], {"utc": plan[0]["utc"], "event": "plan_hashed", "file": "results/datasets.tsv",
                                   "sha256": digest})
        fit = [r for r in read_tsv(RECORD / "order.tsv") if r["event"] == "first_kover_fit"]
        self.assertLess(plan[0]["utc"], fit[0]["utc"])

    def test_sizes(self):
        datasets = read_tsv(RECORD / "datasets.tsv")
        totals = {r["key"]: r["value"] for r in read_tsv(RECORD / "totals.tsv")}
        self.assertEqual(len(datasets), 12)
        self.assertEqual(len(read_tsv(RECORD / "folds.tsv")), 60)
        self.assertEqual(len(read_tsv(RECORD / "per_dataset.tsv")), 12)
        self.assertEqual(len(read_tsv(RECORD / "ties.tsv")), int(totals["kover_one_rule_folds"]))
        self.assertEqual(sum(int(d["genomes"]) for d in datasets), int(totals["genomes"]))
        self.assertEqual(len({d["species"] for d in datasets}), 12)

    def test_no_row_per_genome(self):
        totals = {r["key"]: r["value"] for r in read_tsv(RECORD / "totals.tsv")}
        for name in record.COLUMNS:
            rows = read_tsv(RECORD / name)
            self.assertLessEqual(len(rows), int(totals["folds"]), name)
            for col in read_tsv(RECORD / name, header_only=True):
                self.assertIsNone(re.search(r"(?i)(genome_id|sample|biosample|assembly|accession|kmer|sequence|"
                                            r"^rule$|rule_text|term|weight|importance|^genome$)", col), (name, col))

    def test_readings_fixed_and_exploratory(self):
        summary = {r["measure"]: r for r in read_tsv(RECORD / "summary.tsv")}
        self.assertEqual({m: r["reading"] for m, r in summary.items()}, record.READING)
        self.assertEqual({m: int(r["seed"]) for m, r in summary.items()}, record.SEEDS)

    def test_the_review_s_figures(self):
        # The figures the review's narrower reading gives, at its 3 decimals.
        summary = {r["measure"]: r for r in read_tsv(RECORD / "summary.tsv")}
        got = {m: tuple("%.3f" % float(summary[m][c]) for c in ("value", "low", "high"))
               for m in ("specificity", "recall", "balanced_accuracy", "ties_equal")}
        self.assertEqual(got, {"specificity": ("-0.104", "-0.215", "-0.022"), "recall": ("0.001", "-0.047", "0.046"),
                               "balanced_accuracy": ("-0.052", "-0.094", "-0.018"),
                               "ties_equal": ("-0.035", "-0.060", "-0.013")})

    def test_the_historical_file_is_left_alone(self):
        self.assertEqual(hashlib.sha256(LEGACY.read_bytes()).hexdigest(), LEGACY_SHA256)
        self.assertEqual(len(json.loads(LEGACY.read_text(encoding="utf-8"))), 73)
        self.assertNotEqual(LEGACY.parent / LEGACY.stem, RECORD)


if __name__ == "__main__":
    unittest.main()
