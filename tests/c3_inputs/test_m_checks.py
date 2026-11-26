"""scripts/c3-inputs/m_checks.py: the three checks before any dataset of M runs.

Same counts runs on SYNTHETIC fits, score files and reports (m_common.py).
The check species and Same build read the committed lists.tsv and
page-binaries.tsv, and Same results reads the committed page-datasets.tsv:
its rows are rebuilt from the record itself, so the check must give them
back, and one cell off must stop it. No fit is read.
"""

import copy
import shutil
import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import common, m_common

checks = common.load("m_checks")
reading = common.load("m_reading")
expect = common.load("expect")
REFERENCE = common.ROOT / "reference" / "input-arms"
GONORRHOEAE = "neisseria_gonorrhoeae"


class CheckSpecies(unittest.TestCase):
    def test_the_check_species_on_the_committed_list_is_neisseria_gonorrhoeae(self):
        self.assertEqual(checks.check_species(checks.read_tsv(common.KIT / "lists.tsv")), (GONORRHOEAE, 6))

    def test_a_tie_on_first_run_datasets_goes_to_the_first_by_name(self):
        rows = [{"species_slug": s, "in_M": m, "first_run": f} for s, m, f in
                (("b_sp", "1", "0"), ("b_sp", "0", "1"), ("a_sp", "1", "0"), ("a_sp", "0", "1"), ("c_sp", "0", "1"))]
        self.assertEqual(checks.check_species(rows), ("a_sp", 1))


class SameBuild(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Read here, not at import: page_faults.py runs this suite in a copy of scripts/ and tests/ alone.
        cls.rows = checks.read_tsv(REFERENCE / "page-binaries.tsv")
        cls.code = next(r["code"] for r in cls.rows if r["species_slug"] == GONORRHOEAE)

    def test_the_recorded_code_sha256_at_the_first_run_commit_passes(self):
        self.assertIn("equal to page-binaries.tsv", checks.same_build(self.rows, GONORRHOEAE, checks.FIRST_RUN_COMMIT, self.code))

    def test_a_build_whose_code_sha256_differs_stops(self):
        other = ("0" if self.code[0] != "0" else "1") + self.code[1:]
        with self.assertRaises(checks.Stop):
            checks.same_build(self.rows, GONORRHOEAE, checks.FIRST_RUN_COMMIT, other)
        with self.assertRaises(checks.Stop):
            checks.same_build(self.rows, GONORRHOEAE, "0" * 40, self.code)


class SameResults(unittest.TestCase):
    """The check species' rows of page-datasets.tsv, rebuilt as datasets.tsv rows from the record itself."""

    def rows(self):
        return [r for r in reading.record_rows(REFERENCE) if r["species_slug"] == GONORRHOEAE]

    def test_the_records_own_rows_are_given_back_at_12_decimals(self):
        got = checks.same_results(REFERENCE, GONORRHOEAE, self.rows())
        self.assertIn("rows of neisseria_gonorrhoeae in page-datasets.tsv given back at 12 decimals", got)

    def test_one_cell_off_in_the_twelfth_decimal_stops(self):
        rows = self.rows()
        r = next(r for r in rows if r["arm"] == "unitigs" and r["balanced_accuracy"] not in ("", None))
        r["balanced_accuracy"] = repr(round(float(r["balanced_accuracy"]) + 2e-12, 12))
        with self.assertRaises(checks.Stop):
            checks.same_results(REFERENCE, GONORRHOEAE, rows)

    def test_a_missing_draw_or_a_missing_row_stops(self):
        rows = self.rows()
        r = next(r for r in rows if r["arm"] == "k31-aa")
        r["scored_draws"] = str(int(r["draws"]) - 1)
        with self.assertRaises(checks.Stop):
            checks.same_results(REFERENCE, GONORRHOEAE, rows)
        rows = [r for r in self.rows() if not (r["arm"] == "program" and r["split"] == "lineage")]
        with self.assertRaises(checks.Stop):
            checks.same_results(REFERENCE, GONORRHOEAE, rows)

    def test_a_record_that_differs_from_its_sha256sums_stops(self):
        with tempfile.TemporaryDirectory() as d:
            ref = Path(d) / "ref"
            shutil.copytree(REFERENCE, ref)
            text = (ref / "page-datasets.tsv").read_text()
            (ref / "page-datasets.tsv").write_text(text.replace("\t0.0\n", "\t0.0 \n", 1))
            with self.assertRaises(checks.Stop):
                checks.same_results(ref, GONORRHOEAE, self.rows())


class SameCounts(unittest.TestCase):
    def run_of(self, top, change=None, refused=()):
        rows = expect.rows(m_common.SPECIES, [m_common.M_SET], [m_common.M_SET], [m_common.M_SET], [42, 43], [45],
                           ["program", "k31-counts"])
        gs = m_common.genomes()
        out = m_common.lay_out(top, rows, lambda e: gs, refused=refused)
        if change:
            change(top)
        return out

    def test_every_fold_comes_back_from_the_scores_and_the_kept_reports(self):
        with tempfile.TemporaryDirectory() as d:
            exp, *rest = self.run_of(d)
            fits, kept, refused = checks.same_counts(exp, *rest)
            self.assertEqual((fits, kept, refused), (len(exp), 3, 0))

    def test_one_fold_count_off_by_one_genome_stops(self):
        def change(top):
            path = Path(top) / "picks" / "k31-counts" / f"{m_common.M_SET}.lineage-s43.pick.json"
            doc = checks.load(path)
            doc["folds"][2]["outer"][3]["tp"] += 1
            doc["folds"][2]["outer"][3]["fn"] -= 1
            m_common.write_json(path, doc)
        with tempfile.TemporaryDirectory() as d:
            exp, *rest = self.run_of(d, change)
            with self.assertRaisesRegex(checks.Stop, "fold 2, gbt, tp"):
                checks.same_counts(exp, *rest)

    def test_a_log_loss_sum_off_stops(self):
        def change(top):
            path = Path(top) / "picks" / "k31" / f"{m_common.M_SET}.study-whole-s42.pick.json"
            doc = checks.load(path)
            doc["folds"][0]["outer"][4]["loss_sum"] += 1e-6
            m_common.write_json(path, doc)
        with tempfile.TemporaryDirectory() as d:
            exp, *rest = self.run_of(d, change)
            with self.assertRaisesRegex(checks.Stop, "loss_sum"):
                checks.same_counts(exp, *rest)

    def test_a_kept_report_unlike_its_score_file_stops(self):
        def change(top):
            path = Path(top) / "reports" / "twin" / f"{m_common.M_SET}.lineage-s45.json"
            doc = checks.load(path)
            doc["splits"][1]["test_fold"][0] = (doc["splits"][1]["test_fold"][0] + 1) % m_common.FOLDS
            m_common.write_json(path, doc)
        with tempfile.TemporaryDirectory() as d:
            exp, *rest = self.run_of(d, change)
            with self.assertRaisesRegex(checks.Stop, "kept report"):
                checks.same_counts(exp, *rest)

    def test_a_fit_with_no_score_file_or_no_kept_report_stops(self):
        for gone in (Path("scores") / "program" / f"{m_common.M_SET}.lineage-s42.scores.json",
                     Path("reports") / "k31" / f"{m_common.M_SET}.lineage-s43.json"):
            with tempfile.TemporaryDirectory() as d:
                exp, *rest = self.run_of(d, lambda top: (Path(top) / gone).unlink())
                with self.assertRaises(checks.Stop):
                    checks.same_counts(exp, *rest)

    def test_an_arms_study_refusal_holds_no_scores_and_is_passed_over(self):
        name = f"program/{m_common.M_SET}.study-whole-s42.pick.json"
        with tempfile.TemporaryDirectory() as d:
            exp, *rest = self.run_of(d, refused={name})
            self.assertEqual(checks.same_counts(exp, *rest)[2], 1)

    def test_the_folds_counts_from_scores_are_grms_cutoff_rule(self):
        # Called R at 0.5 and above; the log loss clips each score to [1e-4, 1 - 1e-4].
        c = checks.calls([True, True, False, False], [0.5, 0.0, 0.4999, 1.0])
        self.assertEqual((c["tp"], c["fn"], c["tn"], c["fp"]), (1, 1, 1, 1))
        import math
        want = -(math.log(0.5) + math.log(1e-4) + math.log(1 - 0.4999) + math.log(1e-4))
        self.assertAlmostEqual(c["loss_sum"], want, places=12)


if __name__ == "__main__":
    unittest.main()
