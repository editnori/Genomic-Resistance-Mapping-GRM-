"""reference/other-database: checksums whole, derived files equal derive(), each file's cells agree with each other, sizes and forms."""

import hashlib
import math
import pathlib
import re
import unittest

from tests.other_database.common import RECORD, load_script, read_tsv

record = load_script()
TABLES = sorted(record.COLUMNS)
LIMIT = 4 * 1024 * 1024


def on_disk():
    return {p.relative_to(RECORD).as_posix(): p.read_bytes() for p in sorted(RECORD.rglob("*")) if p.is_file()}


def compared(lo, hi, ref_lo, ref_hi, low, high):
    """The readings a 4-decimal interval allows against a 4-decimal reference interval.

    Rounding keeps order, so a cell above another means the value is above it; equal cells allow both."""
    lo, hi, ref_lo, ref_hi = (float(x) for x in (lo, hi, ref_lo, ref_hi))
    out = set()
    if hi <= ref_lo:
        out.add(low)
    if lo >= ref_hi:
        out.add(high)
    if hi >= ref_lo and lo <= ref_hi:
        out.add("consistent")
    return out


def held(*readings):
    return "yes" if "below" not in readings and "worse" not in readings else "no"


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

    def test_reading_cells_agree(self):
        """Each reading of readings.tsv from its own bounds and its reference, and each *_holds from its two readings."""
        ds = {(d["dataset"], d["learner"]): d for d in self.datasets}
        for r in self.readings:
            d = ds[(r["dataset"], r["learner"])]
            with self.subTest(dataset=r["dataset"], learner=r["learner"]):
                # L: the test's interval against the lineage figure's own interval.
                self.assertIn(r["L_recall"], compared(d["recall_lo"], d["recall_hi"], r["L_lineage_recall_lo"], r["L_lineage_recall_hi"], "below", "above"))
                self.assertIn(r["L_log_loss"], compared(d["log_loss_lo"], d["log_loss_hi"], r["L_lineage_log_loss_lo"], r["L_lineage_log_loss_hi"], "better", "worse"))
                self.assertEqual(r["L_holds"], held(r["L_recall"], r["L_log_loss"]))
                # C: the class-balanced interval against the lineage split's balanced log loss; recall as the plan reads it.
                self.assertIn(r["C_log_loss"], compared(r["C_log_loss_lo"], r["C_log_loss_hi"], r["C_lineage_log_loss"], r["C_lineage_log_loss"], "better", "worse"))
                self.assertLessEqual(float(r["C_log_loss_lo"]), float(r["C_log_loss_hi"]))
                self.assertEqual(r["C_holds"], held(d["recall_reading"], r["C_log_loss"]))
                # D: the BioProject intervals against the lineage figures.
                self.assertIn(r["D_recall"], compared(r["D_recall_lo"], r["D_recall_hi"], d["lineage_recall"], d["lineage_recall"], "below", "above"))
                self.assertIn(r["D_log_loss"], compared(r["D_log_loss_lo"], r["D_log_loss_hi"], d["lineage_log_loss"], d["lineage_log_loss"], "better", "worse"))
                self.assertEqual(r["D_holds"], held(r["D_recall"], r["D_log_loss"]))
                # T, the pilot's two-sample reading.
                if r["part"] == "pilot":
                    self.assertEqual(r["T_holds"], held(r["T_recall"], r["T_log_loss"]))
                self.assertEqual(int(r["errors"]), int(d["FN"]) + int(d["FP"]))

    def test_plan_reading_from_its_bounds(self):
        """recall_reading and log_loss_reading of datasets.tsv allow ties of the 4-decimal cells, as the readings do."""
        for d in self.datasets:
            if d["estimable"] == "yes":
                self.assertIn(d["recall_reading"], compared(d["recall_lo"], d["recall_hi"], d["lineage_recall"], d["lineage_recall"], "below", "above"))
                self.assertIn(d["log_loss_reading"], compared(d["log_loss_lo"], d["log_loss_hi"], d["lineage_log_loss"], d["lineage_log_loss"], "better", "worse"))

    def test_order_cells_agree(self):
        """order.tsv: the yes or no cells from the counts and times beside them, and the plan file's time from plans.tsv."""
        fixed = {r["file"]: r["fixed"] for r in record.parse(self.text["plans.tsv"])}
        for r in record.parse(self.text["order.tsv"]):
            with self.subTest(species=r["species"]):
                if r["assigned"] == "by BioSample":
                    self.assertEqual(r["asked_equals_scored"], "yes" if r["biosamples_asked"] == r["scored_biosamples"] else "no")
                    self.assertEqual(r["plan_file"], "ADDENDUM2.md")
                else:
                    self.assertEqual((r["biosamples_asked"], r["scored_biosamples"], r["asked_equals_scored"], r["plan_file"]), ("", "", "", "PLAN.md"))
                self.assertEqual(r["first_at_or_after_predictions"], "yes" if r["first_label_request"] >= r["predictions_written"] else "no")
                self.assertEqual(r["first_after_plan_file"], "yes" if r["first_label_request"] > fixed[r["plan_file"]] else "no")
                self.assertLessEqual(r["first_label_request"], r["last_label_request"])
                self.assertLessEqual(int(r["answered_200"]), int(r["label_requests"]))

    def test_each_genome_once(self):
        rows = record.parse(self.text["genomes.tsv"])
        for key in ("biosample", "assembly"):
            seen = [r[key] for r in rows]
            self.assertEqual(len(seen), len(set(seen)), f"a {key} is given twice")

    def test_seeds_file_agrees(self):
        """seeds.tsv: seed 1 is the plan's result, and the species rows are exactly the species results its dataset rows change."""
        seeds = record.parse(self.text["seeds.tsv"])
        names = [c for c in record.COLUMNS["seeds.tsv"] if c.startswith("seed_")]
        self.assertEqual(names, [f"seed_{s}" for s in (record.SEED,) + record.OTHER_SEEDS])
        ds = {(d["dataset"], d["learner"]): d for d in self.datasets}
        learners = {(r["species"], r["learner"], r["part"]): r for r in record.parse(self.text["learners.tsv"])}
        moved = [r for r in seeds if r["level"] == "dataset"]
        for r in moved:
            d = ds[(r["dataset"], r["learner"])]
            self.assertEqual((d["estimable"], d["species"], d["part"], r["seed_1"]), ("yes", r["species"], r["part"], d["holds"]))
            self.assertGreater(len({r[c] for c in names}), 1)
        # Every other estimable dataset holds as at seed 1; the species results follow from the counts.
        expect = []
        for (slug, ln, part), row in sorted(learners.items()):
            est = int(row["estimable"])
            results = []
            for c in names:
                n = int(row["holds"]) + sum((r[c] == "yes") - (r["seed_1"] == "yes") for r in moved if (r["species"], r["learner"]) == (slug, ln))
                results.append(record.rule(n, est))
            if len(set(results)) > 1:
                expect.append((slug, ln, part, results))
        got = [(r["species"], r["learner"], r["part"], [r[c] for c in names]) for r in seeds if r["level"] == "species"]
        self.assertEqual(sorted(got), expect)

    def test_lineage_at_share_by_hand(self):
        """The lineage split's error rate and SCM's log loss at the test's resistant share, from the 4-decimal cells."""
        hi, lo = -math.log(record.EPS), -math.log(1 - record.EPS)
        for d in self.datasets:
            if d["estimable"] != "yes":
                self.assertEqual((d["lineage_error_rate_at_share"], d["lineage_log_loss_at_share"]), ("", ""))
                continue
            share = int(d["R"]) / (int(d["R"]) + int(d["S"]))
            error = share * (1 - float(d["lineage_recall"])) + (1 - share) * (1 - float(d["lineage_specificity"]))
            self.assertAlmostEqual(float(d["lineage_error_rate_at_share"]), error, delta=1.5e-4)
            if d["learner"] == "scm":
                self.assertAlmostEqual(float(d["lineage_log_loss_at_share"]), error * hi + (1 - error) * lo, delta=1.5e-3)
        # N. gonorrhoeae azithromycin, SCM: 95 R and 84 S, lineage recall 0.4240 and specificity 0.9732.
        azi = next(d for d in self.datasets if (d["dataset"], d["learner"]) == ("neisseria_gonorrhoeae__azithromycin", "scm"))
        self.assertEqual((azi["R"], azi["S"], azi["lineage_recall"], azi["lineage_specificity"]), ("95", "84", "0.4240", "0.9732"))
        error = 95 / 179 * (1 - 0.4240) + 84 / 179 * (1 - 0.9732)
        self.assertAlmostEqual(error, 0.318275, delta=1e-6)
        self.assertEqual(azi["lineage_error_rate_at_share"], "0.3183")
        self.assertAlmostEqual(float(azi["lineage_log_loss_at_share"]), error * hi + (1 - error) * lo, delta=1e-3)

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

    def test_absolute_path_only_in_the_plan_copy(self):
        # plan/ has the plan files byte for byte, as their sha256 fixed them; PLAN.md names the data root's BV-BRC directory once.
        for name, text in self.text.items():
            found = re.findall(r"(?<![\w.])/(?:home|tmp|Users|mnt)/", text)
            self.assertEqual(len(found), 1 if name == "plan/PLAN.md" else 0, name)
            self.assertIsNone(re.search(r"[A-Z]:\\\\", text), name)


if __name__ == "__main__":
    unittest.main()
