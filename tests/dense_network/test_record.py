"""reference/dense-network: its hashes, columns and size, the derived files, and the figures the experiment and
the review printed."""

import hashlib
import re
import unittest

from tests.dense_network.common import RECORD, load_script, read_tsv

record = load_script()


def record_files():
    return sorted(p.name for p in RECORD.iterdir() if p.is_file())


def keyed(name):
    return {r["key"]: r["value"] for r in read_tsv(RECORD / name)}


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
        self.assertEqual(sorted(record.COLUMNS), sorted(n for n in record_files() if n.endswith(".tsv")))

    def test_under_3_mb(self):
        self.assertLess(sum((RECORD / n).stat().st_size for n in record_files()), 3 * 1024 * 1024)

    def test_columns(self):
        for name, cols in record.COLUMNS.items():
            self.assertEqual(read_tsv(RECORD / name, header_only=True), cols, name)
            for row in read_tsv(RECORD / name):
                self.assertEqual(len(row), len(cols), name)
                self.assertNotIn(None, row.values(), name)

    def test_readme_names_each_file_and_its_columns(self):
        readme = (RECORD / "README.md").read_text(encoding="utf-8")
        for name in record_files():
            self.assertIn("| `%s` |" % name, readme)
            if name.endswith(".tsv"):
                cols = read_tsv(RECORD / name, header_only=True)
                self.assertIn("\n- `%s`: %s.\n" % (name, ", ".join("`%s`" % c for c in cols)), readme)

    def test_record_is_what_the_script_writes(self):
        try:
            files, _, _ = record.assemble(record.from_record(RECORD))
        except record.Stop as exc:
            self.fail("the script stops on the record: %s" % exc)
        self.assertEqual(record.check(files, RECORD), [])

    def test_one_row_per_dataset_split_and_learner(self):
        datasets = read_tsv(RECORD / "datasets.tsv")
        totals = keyed("totals.tsv")
        names = [d["dataset"] for d in datasets]
        keys = [(r["dataset"], r["split"], r["learner"]) for r in read_tsv(RECORD / "measures.tsv")]
        self.assertEqual(len(keys), len(set(keys)))
        want = {(d, s, w) for d in names for s in record.SPLITS for w in record.LEARNERS if w != "catalog"}
        self.assertEqual({k for k in keys if k[2] != "catalog"}, want)
        catalog = {k[0] for k in keys if k[2] == "catalog"}
        self.assertEqual(len(catalog), int(totals["catalog_datasets"]))
        self.assertEqual(len(keys), len(record.SPLITS) * (len(names) * (len(record.LEARNERS) - 1) + len(catalog)))

    def test_no_row_per_genome(self):
        datasets = read_tsv(RECORD / "datasets.tsv")
        fewest = min(int(d["genomes"]) for d in datasets)
        for name in record.COLUMNS:
            rows = read_tsv(RECORD / name)
            if name != "measures.tsv":
                self.assertLess(len(rows), fewest, name)
            for col in read_tsv(RECORD / name, header_only=True):
                self.assertIsNone(re.search(r"(?i)(genome_id|sample|biosample|assembly|accession|kmer|sequence|"
                                            r"feature|weight|importance|^genome$|^row$)", col), (name, col))

    def test_sizes(self):
        datasets = read_tsv(RECORD / "datasets.tsv")
        totals = keyed("totals.tsv")
        self.assertEqual(len(datasets), 24)
        self.assertEqual([d["group"] for d in datasets], ["largest_gap"] * 12 + ["drawn"] * 12)
        self.assertEqual(len({d["species"] for d in datasets}), int(totals["species"]))
        self.assertEqual(int(totals["species"]), 12)
        self.assertEqual(sum(1 for d in datasets if d["gate_a"] == d["gate_b"] == "yes"), 24)
        self.assertEqual(int(keyed("cost.tsv")["fits"]), 24 * 2 * 5 * 3)

    def test_order(self):
        order = read_tsv(RECORD / "order.tsv")
        self.assertEqual([r["event"] for r in order], record.ORDER_EVENTS)
        utc = [r["utc"] for r in order]
        self.assertEqual(utc, sorted(utc))
        self.assertRegex(order[0]["sha256"], r"^[0-9a-f]{64}$")
        self.assertEqual([r["sha256"] for r in order[1:]], [""] * (len(order) - 1))

    def test_settings_are_the_script_s(self):
        self.assertEqual([(r["key"], r["value"]) for r in read_tsv(RECORD / "settings.tsv")], record.SETTINGS)
        self.assertEqual(set(record.SETTINGS_IN_PLAN) - {k for k, _ in record.SETTINGS}, set())

    def test_the_experiment_s_printed_answer(self):
        # The experiment's answer and the review's recompute, as printed at 4 decimals.
        summary = {(r["split"], r["comparator"], r["measure"]): r for r in read_tsv(RECORD / "summary.tsv")}

        def at4(split, who, k):
            r = summary[(split, who, k)]
            return tuple("%+.4f" % float(r[c]) for c in ("difference_median", "low", "high"))

        self.assertEqual(at4("lineage", "best5", "balanced_accuracy"), ("-0.0593", "-0.1323", "-0.0445"))
        self.assertEqual(at4("lineage", "best5", "recall"), ("-0.1020", "-0.3561", "+0.0000"))
        self.assertEqual(at4("lineage", "best5", "log_loss"), ("+0.0210", "-0.2621", "+0.2091"))
        self.assertEqual(at4("lineage", "cart", "balanced_accuracy"), ("-0.0123", "-0.0522", "-0.0007"))
        self.assertEqual(at4("random", "best5", "balanced_accuracy"), ("-0.0426", "-0.0640", "-0.0256"))
        b = summary[("lineage", "best5", "balanced_accuracy")]
        self.assertEqual((b["positive"], b["negative"], b["datasets"], b["species"]), ("1", "23", "24", "12"))
        rule = read_tsv(RECORD / "rule.tsv")
        self.assertEqual(["%.4f" % float(rule[i]["value"]) for i in (2, 4)], ["0.1154", "1.8384"])
        self.assertEqual(rule[2]["needs"], "at least 0.203049, the comparator's")
        self.assertEqual(rule[5]["value"], "299972")
        self.assertEqual([r["holds"] for r in rule], ["no", "no", "no", "no", "yes", "yes", "no"])
        self.assertEqual(rule[-1]["value"], "no gain established")

    def test_the_review_s_figures(self):
        totals = keyed("totals.tsv")
        below = {w: "%.4f" % float(totals["below_best5_lineage_" + w]) for w in record.FIVE}
        self.assertEqual(below, {"rf": "0.0242", "cart": "0.0475", "scm": "0.0517", "gbt": "0.0847", "lasso": "0.0876"})
        self.assertEqual((totals["site_best_is_best5_random"], totals["site_best_is_best5_lineage"]), ("24", "14"))
        # The record's own specificity figure, which the second review counted at +0.025 (-0.106 to +0.181).
        summary = {(r["split"], r["comparator"], r["measure"]): r for r in read_tsv(RECORD / "summary.tsv")}
        spec = summary[("lineage", "best5", "specificity")]
        self.assertEqual(tuple("%+.3f" % float(spec[c]) for c in ("difference_median", "low", "high")),
                         ("+0.025", "-0.106", "+0.181"))
        self.assertEqual((spec["positive"], spec["datasets"]), ("13", "24"))
        order = {r["event"]: r["utc"] for r in read_tsv(RECORD / "order.tsv")}
        self.assertEqual(order["plan_file_time"], "2026-11-15T12:19:40Z")
        self.assertEqual(order["first_network_result"], "2026-11-15T12:21:31Z")


if __name__ == "__main__":
    unittest.main()
