"""The checks of scripts/matched-accuracy/stage1.py that stop the run: one line, exit 1 and nothing written.

Each test changes one input of the SYNTHETIC release of tests/matched_accuracy/made_up.py.
"""

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests.matched_accuracy import common, made_up
from tests.matched_accuracy.common import report_of

AL_A = "synthia_alpha__drug_a"


def append(name, data=b" "):
    def change(top):
        with open(top / name, "ab") as f:
            f.write(data)
    return change


class Sums(common.Case, unittest.TestCase):
    def test_a_lab_records_file_other_than_its_sha256_stops(self):
        self.assert_stops(self.run_stage1(after=append("genome_amr.jsonl")), "the lab records file differs from its sha256")

    def test_a_second_genome_copy_other_than_its_sha256_stops(self):
        run = self.run_stage1(after=append("master-genome-meta.jsonl.gz"))
        self.assert_stops(run, "the genome table, second copy file differs from its sha256")

    def test_a_drug_names_file_other_than_its_sha256_stops(self):
        self.assert_stops(self.run_stage1(after=append("drug_names.json")), "the drug names file differs from its sha256")

    def test_a_label_table_other_than_its_sha256_stops(self):
        self.assert_stops(self.run_stage1(after=append("phenotypes.jsonl")), "the label table file differs from its sha256")

    def test_a_5_nov_genome_table_other_than_its_sha256_stops(self):
        self.assert_stops(self.run_stage1(after=append("genomes.jsonl")),
                          "the genome table, 5 Nov file differs from its sha256")

    def test_a_report_other_than_its_line_of_the_sums_list_stops(self):
        run = self.run_stage1(after=append(f"reports/{AL_A}.json"))
        self.assert_stops(run, f"{AL_A}.json differs from its line of the release's sums list")

    def test_a_reports_directory_with_one_report_more_stops(self):
        def more(top):
            (top / "reports" / "synthia_alpha__drug_z.json").write_text("{}")
        self.assert_stops(self.run_stage1(after=more), "the reports directory holds 13 reports")

    def test_a_sums_list_other_than_its_sha256_stops(self):
        run = self.run_stage1(pins={"REPORT_SUMS_SHA256": "0" * 64})
        self.assert_stops(run, "the release's sums list of reports differs from its sha256")

    def test_every_sha256_is_checked_before_any_other_read(self):
        # Records that are not JSON and a label table other than its sha256: the sha256 line comes first.
        rel = common.release()
        rel["records_bytes"] = b"not json\n"

        def change(top):
            append("phenotypes.jsonl")(top)
        run = self.run_stage1(rel, after=change)
        self.assert_stops(run, "the label table file differs from its sha256")


class Records(common.Case, unittest.TestCase):
    def test_a_record_inserted_at_the_seal_stops(self):
        rel = common.release()
        rel["records"][5]["date_inserted"] = "2026-11-05T19:10:00Z"
        self.assert_stops(self.run_stage1(rel), "a lab record was inserted at or after 2026-11-05T19:10:00")

    def test_a_record_inserted_after_the_seal_stops(self):
        rel = common.release()
        rel["records"][-1]["date_inserted"] = "2026-11-21T00:00:00.000Z"
        self.assert_stops(self.run_stage1(rel), "a lab record was inserted at or after 2026-11-05T19:10:00")

    def test_a_record_inserted_in_the_last_second_before_the_seal_passes(self):
        rel = common.release()
        rel["records"][5]["date_inserted"] = "2026-11-05T19:09:59.999999Z"
        self.assert_passes(self.run_stage1(rel))

    def test_a_record_with_a_blank_insertion_date_stops(self):
        rel = common.release()
        rel["records"][5]["date_inserted"] = " "
        self.assert_stops(self.run_stage1(rel), "a lab record's date_inserted is blank or not a time to the second")

    def test_a_line_count_other_than_the_pinned_one_stops(self):
        for change in (1, -1):
            with self.subTest(change=change):
                probe = self.run_stage1()
                lines = probe.pins["RECORDS_LINES"]
                run = self.run_stage1(pins={"RECORDS_LINES": lines + change})
                self.assert_stops(run, f"the lab records hold {lines} lines, not {lines + change}")

    def test_a_record_name_missing_from_the_drug_names_file_stops(self):
        rel = common.release()
        rel["records"][3]["antibiotic"] = "Drug Y"
        self.assert_stops(self.run_stage1(rel), "a lab record's antibiotic name is not in the drug names file")


class Reports(common.Case, unittest.TestCase):
    def split(self, rel, dataset=AL_A):
        return next(s for s in report_of(rel, dataset)["splits"] if s["scheme"] == "lineage")

    def test_a_test_fold_shorter_than_the_genome_rows_stops(self):
        rel = common.release()
        self.split(rel)["test_fold"].pop()
        self.assert_stops(self.run_stage1(rel), f"{AL_A}: test_fold is not a list as long as the report's 10 genome rows")

    def test_a_scores_list_longer_than_the_genome_rows_stops(self):
        rel = common.release()
        self.split(rel)["learners"][2]["scores"].append(0.5)
        self.assert_stops(self.run_stage1(rel), f"{AL_A}: the scores of rf are not a list as long as the report's 10")

    def test_a_lineage_split_with_other_than_five_folds_stops(self):
        def four_folds(rel):
            self.split(rel)["folds"] = 4

        def fold_four_empty(rel):
            fold = self.split(rel)["test_fold"]
            fold[:] = [3 if f == 4 else f for f in fold]
        for change in (four_folds, fold_four_empty):
            with self.subTest(change=change.__name__):
                rel = common.release()
                change(rel)
                self.assert_stops(self.run_stage1(rel), f"{AL_A}: the lineage split has other than 5 folds")

    def test_a_lineage_cluster_in_two_folds_stops(self):
        rel = common.release()
        self.split(rel)["test_fold"][1] = 1     # position 2, lineage 10, moves to fold 1; position 1 stays in fold 0
        self.assert_stops(self.run_stage1(rel), f"{AL_A}: a lineage cluster sits in two folds")

    def test_a_copy_group_in_two_folds_stops(self):
        rel = common.release()
        report_of(rel, AL_A)["kset_duplicates"] = [[0, 2]]   # positions 1 and 3: folds 0 and 1
        self.assert_stops(self.run_stage1(rel), f"{AL_A}: a copy group sits in two folds")

    def test_a_report_with_four_of_the_five_learners_stops(self):
        rel = common.release()
        learners = self.split(rel)["learners"]
        learners[:] = [e for e in learners if e["learner"] != "gbt"]
        self.assert_stops(self.run_stage1(rel), f"{AL_A}: the lineage split holds 4 of the five learners")


class Calls(common.Case, unittest.TestCase):
    def test_a_row_whose_call_differs_from_the_label_table_stops(self):
        rel = common.release()
        p = next(p for p in rel["labels"] if (p["genome_id"], p["antibiotic"]) == ("gz-al-01", "drug_a"))
        p["reported"] = "S"
        self.assert_stops(self.run_stage1(rel), f"{AL_A}: a row's lab call differs from the label table's call")

    def test_a_row_whose_pair_the_label_table_lacks_stops(self):
        rel = common.release()
        rel["labels"] = [p for p in rel["labels"] if (p["genome_id"], p["antibiotic"]) != ("gz-al-02", "drug_a")]
        self.assert_stops(self.run_stage1(rel), f"{AL_A}: a row's lab call differs from the label table's call")

    def test_a_row_whose_table_pair_is_labelled_from_an_mic_stops(self):
        rel = common.release()
        p = next(p for p in rel["labels"] if (p["genome_id"], p["antibiotic"]) == ("gz-al-03", "drug_a"))
        p["label_from_mic"] = True
        self.assert_stops(self.run_stage1(rel), f"{AL_A}: a row's lab call differs from the label table's call")


class Tables(common.Case, unittest.TestCase):
    def test_a_genome_table_row_inserted_at_the_seal_stops(self):
        rel = common.release()
        rel["genomes"][3]["date_inserted"] = "2026-11-05T19:10:00Z"
        self.assert_stops(self.run_stage1(rel), "a row of the genome table, second copy, is not inserted before")


    def test_a_second_genome_copy_that_does_not_read_to_its_end_stops(self):
        rel = common.release()
        rel["genomes_bytes"] = made_up.gz_bytes(made_up.jsonl_bytes(rel["genomes"]))[:-8]
        self.assert_stops(self.run_stage1(rel), "the genome table, second copy is not a whole gzip file")

    def test_a_record_that_is_not_json_stops(self):
        rel = common.release()
        rel["records_bytes"] = made_up.jsonl_bytes(rel["records"]) + b"not json\n"
        self.assert_stops(self.run_stage1(rel), "line 128 of the lab records is not JSON")


class Groups(common.Case, unittest.TestCase):
    def test_no_group_in_reading_a_stops(self):
        rel = common.release()
        rel["records"] = []
        self.assert_stops(self.run_stage1(rel), "no group qualifies in reading A")

    def test_no_group_in_reading_b_stops(self):
        rel = common.release()
        rel["records"] += [made_up.record(r["genome_id"], r["antibiotic"], call="Intermediate")
                           for r in list(rel["records"]) if r["antibiotic"].strip() in made_up.RAW.values()]
        self.assert_stops(self.run_stage1(rel), "no group qualifies in reading B")


class Output(common.Case, unittest.TestCase):
    def test_an_output_directory_that_holds_an_export_stops(self):
        def before(top):
            (top.parent / "out" / "export").mkdir(parents=True)
        run = self.run_stage1(after=before)
        self.assertEqual(run.code, 1)
        self.assertEqual(run.lines, ["stage1: stopped: the output directory already holds local or export"])
        self.assertEqual(sorted(p.name for p in run.out.iterdir()), ["export"])
        self.assertEqual(list((run.out / "export").iterdir()), [])


class Command(unittest.TestCase):
    """The script as the real run calls it, under -I, with its real pins."""

    def test_the_command_with_its_real_pins_stops_on_made_up_reports(self):
        with tempfile.TemporaryDirectory() as d:
            argv, _ = made_up.write_release(Path(d) / "in", common.release())
            out = Path(d) / "out"
            p = subprocess.run([sys.executable, "-I", "-B", str(common.SCRIPT)] + argv + ["--out", str(out)],
                               capture_output=True, text=True, env={"PATH": os.environ.get("PATH", "/usr/bin:/bin")})
            self.assertEqual(p.returncode, 1, p.stderr)
            self.assertEqual(p.stderr, "")
            self.assertEqual(p.stdout.splitlines(),
                             ["stage1: stopped: the reports directory holds 12 reports; the release's sums list names "
                              "274 and 286 differ"])
            self.assertFalse(out.exists())

    def test_the_command_without_its_arguments_prints_its_usage(self):
        p = subprocess.run([sys.executable, "-I", "-B", str(common.SCRIPT)], capture_output=True, text=True)
        self.assertEqual(p.returncode, 2)
        self.assertEqual(p.stdout, "")
        self.assertIn("usage: stage1.py", p.stderr)


if __name__ == "__main__":
    unittest.main()
