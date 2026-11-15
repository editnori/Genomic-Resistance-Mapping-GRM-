"""reference/other-database: checksums whole, derived files equal derive(), totals agree, sizes and forms."""

import hashlib
import pathlib
import re
import unittest

from tests.other_database.common import RECORD, load_script, read_tsv

record = load_script()
TABLES = sorted(record.COLUMNS)
LIMIT = 4 * 1024 * 1024


def on_disk():
    return {p.relative_to(RECORD).as_posix(): p.read_bytes() for p in sorted(RECORD.rglob("*")) if p.is_file()}


class Record(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.files = on_disk()
        cls.text = {k: v.decode("utf-8") for k, v in cls.files.items()}
        cls.datasets = record.parse(cls.text["datasets.tsv"])
        cls.readings = record.parse(cls.text["readings.tsv"])

    def test_checksums_whole(self):
        sums = dict(reversed(line.split("  ")) for line in self.text["SHA256SUMS"].splitlines())
        self.assertEqual(set(sums), set(self.files) - {"SHA256SUMS"})
        for name, digest in sums.items():
            self.assertEqual(hashlib.sha256(self.files[name]).hexdigest(), digest, name)

    def test_under_four_megabytes(self):
        self.assertLess(sum(len(v) for v in self.files.values()), LIMIT)

    def test_columns_and_order(self):
        for name in TABLES:
            with self.subTest(name=name):
                self.assertEqual(read_tsv(RECORD / name, header_only=True), record.COLUMNS[name])
                rows = record.parse(self.text[name])
                self.assertTrue(rows)
                keys = [tuple(r[k] for k in record.KEYS[name]) for r in rows]
                self.assertEqual(keys, sorted(keys))
                self.assertEqual(len(keys), len(set(keys)), "a key is written twice")
                self.assertEqual(self.text[name], record.tsv(name, rows))

    def test_derived_files_equal_derive(self):
        learners, totals, differences = record.derive(self.datasets, self.readings)
        self.assertEqual(self.text["learners.tsv"], record.tsv("learners.tsv", learners))
        self.assertEqual(self.text["totals.tsv"], record.tsv("totals.tsv", totals))
        self.assertEqual(self.text["differences.tsv"], record.tsv("differences.tsv", differences))

    def test_totals_agree_with_learners(self):
        learners = [r for r in record.parse(self.text["learners.tsv"]) if r["part"] == "eight"]
        for total in record.parse(self.text["totals.tsv"]):
            rows = [r for r in learners if r["learner"] == total["learner"]]
            self.assertEqual(int(total["species"]), len(rows))
            for k in ("estimable", "holds", "recall_below", "log_loss_worse", "specificity_under_half"):
                self.assertEqual(int(total[k]), sum(int(r[k]) for r in rows), (total["learner"], k))
            self.assertEqual(int(total["plan_meets"]), sum(r["rule"] == "meets" for r in rows))
            for v in ("plan", "L", "C", "D"):
                self.assertEqual(sum(int(total[f"{v}_{x}"]) for x in ("meets", "fails", "not_assessable")), len(rows))

    def test_learners_count_datasets(self):
        for r in record.parse(self.text["learners.tsv"]):
            est = [d for d in self.datasets if (d["species"], d["learner"], d["estimable"]) == (r["species"], r["learner"], "yes")]
            self.assertEqual(int(r["estimable"]), len(est))
            self.assertEqual(int(r["holds"]), sum(d["holds"] == "yes" for d in est))
            self.assertEqual(r["rule"], record.rule(int(r["holds"]), int(r["estimable"])))

    def test_dataset_cells_agree(self):
        for d in self.datasets:
            with self.subTest(dataset=d["dataset"], learner=d["learner"]):
                R, S, tp, fn, tn, fp = (int(d[k]) for k in ("R", "S", "TP", "FN", "TN", "FP"))
                self.assertEqual((tp + fn, tn + fp), (R, S))
                est = R >= record.MIN_CLASS and S >= record.MIN_CLASS and d["lineage_figure"] == "yes"
                self.assertEqual(d["estimable"], "yes" if est else "no")
                if not est:
                    self.assertEqual(d["holds"], "")
                    continue
                self.assertEqual(d["recall"], f"{tp / R:.4f}")
                self.assertEqual(d["specificity"], f"{tn / S:.4f}")
                self.assertEqual((d["recall_lo"], d["recall_hi"]), tuple(f"{v:.4f}" for v in record.clopper_pearson(tp, R)))
                held = d["recall_reading"] != "below" and d["log_loss_reading"] != "worse"
                self.assertEqual(d["holds"], "yes" if held else "no")
                self.assertLessEqual(float(d["log_loss_lo"]), float(d["log_loss"]))
                self.assertLessEqual(float(d["log_loss"]), float(d["log_loss_hi"]))
                below = float(d["recall_hi"]) < float(d["lineage_recall"])
                self.assertEqual(d["recall_reading"] == "below", below)
                worse = float(d["log_loss_lo"]) > float(d["lineage_log_loss"])
                self.assertEqual(d["log_loss_reading"] == "worse", worse)

    def test_readings_cover_estimable(self):
        est = {(d["dataset"], d["learner"]) for d in self.datasets if d["estimable"] == "yes"}
        self.assertEqual({(r["dataset"], r["learner"]) for r in self.readings}, est)
        for r in self.readings:
            self.assertEqual(r["status"], "exploratory" if r["part"] == "pilot" else "fixed in advance")
            self.assertEqual(r["T_holds"] != "", r["part"] == "pilot")

    def test_genomes_hold_accessions_only(self):
        self.assertEqual(read_tsv(RECORD / "genomes.tsv", header_only=True), ["species", "assembly", "biosample", "bioprojects"])
        species = set(record.NAMES)
        for r in record.parse(self.text["genomes.tsv"]):
            self.assertIn(r["species"], species)
            self.assertRegex(r["assembly"], r"^GC[AF]_\d{9}\.\d+$")
            self.assertRegex(r["biosample"], r"^SAM[NED][A-Z]?\d+$")
            self.assertRegex(r["bioprojects"], r"^(PRJ[NED][A-Z]\d+(;PRJ[NED][A-Z]\d+)*)?$")
        scored = {r["species"]: int(r["scored"]) for r in record.parse(self.text["species.tsv"])}
        rows = record.parse(self.text["genomes.tsv"])
        self.assertEqual({s: sum(r["species"] == s for r in rows) for s in scored}, scored)

    def test_inputs_named_by_root(self):
        for line in self.text["inputs.sha256"].splitlines():
            digest, path = line.split("  ")
            self.assertRegex(digest, r"^[0-9a-f]{64}$")
            self.assertIn(path.split("/")[0], ("lanes", "release", "data-root"))
            self.assertNotIn("..", path)
        self.assertEqual(self.text["inputs.sha256"], "".join(sorted(self.text["inputs.sha256"].splitlines(keepends=True), key=lambda x: x.split("  ")[1])))

    def test_plans_match_their_copies(self):
        for r in record.parse(self.text["plans.tsv"]):
            self.assertEqual(hashlib.sha256(self.files["plan/" + r["file"]]).hexdigest(), r["sha256"])
            inputs = dict(reversed(line.split("  ")) for line in self.text["inputs.sha256"].splitlines())
            self.assertEqual(inputs["lanes/" + dict((p[0], p[2]) for p in record.PLANS)[r["file"]]], r["sha256"])
            self.assertRegex(r["fixed"], r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ$")

    def test_settings(self):
        s = {r["key"]: r["value"] for r in record.parse(self.text["settings.tsv"])}
        self.assertEqual(self.text["settings.tsv"], record.tsv("settings.tsv", record.settings_rows(s["release"], s["evaluation_revision"])))
        self.assertEqual((s["release"], s["evaluation_revision"]), ("rev14-274-20261112", "14"))

    def test_readme_rebuilt(self):
        rows_of = {k: str(len(record.parse(self.text[k]))) for k in TABLES}
        rows_of["inputs.sha256"] = str(len(self.text["inputs.sha256"].splitlines()))
        settings = record.parse(self.text["settings.tsv"])
        self.assertEqual(self.text["README.md"], record.readme(sorted(self.files), rows_of, settings))

    def test_no_absolute_path(self):
        # plan/ holds the lanes' plan files byte for byte, as their sha256 fixed them.
        for name, text in self.text.items():
            if name.startswith("plan/"):
                continue
            self.assertNotRegex(text, r"(?<![\w.])/(home|tmp|Users|mnt)/", name)
            self.assertIsNone(re.search(r"[A-Z]:\\\\", text), name)


if __name__ == "__main__":
    unittest.main()
