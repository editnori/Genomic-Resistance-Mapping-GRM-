"""scripts/c7-outside/page_tables.py on a SYNTHETIC reading: each gate refuses what it names, and the tables hold
what docs/outside-learners-result.md prints. Every species, dataset, call and count is made up (page_fixture.py)."""

import contextlib
import csv
import hashlib
import io
import json
import re
import shutil
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

from tests.c7_outside import page_fixture as F

page = F.page


def rows_of(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def edit(path, change):
    """Rewrite one tsv after change(rows) changes its rows in place; the columns stay."""
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns, rows = reader.fieldnames, list(reader)
    change(rows)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for r in rows:
            handle.write("\t".join(r[c] for c in columns) + "\n")


class Base(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        F.build(cls.root / "run")
        F.build(cls.root / "withdrawn", withdrawn=True)
        cls.n = 0

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def reading(self, which="run"):
        """A fresh copy of one SYNTHETIC reading, to change."""
        type(self).n += 1
        dst = self.root / f"copy-{self.n}"
        shutil.copytree(self.root / which / "reading", dst / "reading")
        return dst / "reading"

    def lists(self, which="run"):
        return self.root / which / "lists.tsv"

    def write(self, reading, which="run"):
        out = reading.parent / "out"
        page.write(reading, self.lists(which), out)
        return out

    def refused(self, reading, gate, which="run"):
        with self.assertRaisesRegex(page.Stop, f"^gate {gate}: "):
            self.write(reading, which)

    def out(self, which="run"):
        return self.write(self.reading(which), which)


class Gates(Base):
    def test_the_synthetic_readings_pass_every_gate(self):
        for which in ("run", "withdrawn"):
            out = self.out(which)
            self.assertTrue((out / "SHA256SUMS").is_file())
            self.assertEqual(len((out / "gates.txt").read_text().splitlines()), len(page.GATES))

    def test_a_folder_outside_the_reading_is_refused(self):
        reading = self.reading()
        (reading / "notes").mkdir()
        self.refused(reading, "1 reading")

    def test_a_reading_whose_checks_failed_is_refused(self):
        reading = self.reading()
        text = (reading / "checks.txt").read_text().replace("PASS features", "FAIL features")
        (reading / "checks.txt").write_text(text)
        self.refused(reading, "1 reading")

    def test_an_arm_count_tsv_marks_withdrawn_with_no_reason_is_refused(self):
        reading = self.reading("withdrawn")
        (reading / "withdrawn.tsv").write_text("arm\treason\ttime\n")
        self.refused(reading, "1 reading", "withdrawn")

    def test_a_second_pass_count_with_no_arm_withdrawn_is_refused(self):
        reading = self.reading()
        text = (reading / "checks.txt").read_text().replace(
            "beside no newer DONE", "beside no newer DONE; the elastic net withdrawn: 1 species ran a second pass, 1 "
            "second-pass copies set aside", 1)
        (reading / "checks.txt").write_text(text)
        self.refused(reading, "1 reading")

    def test_a_withdrawn_net_with_no_second_pass_count_is_refused(self):
        reading = self.reading("withdrawn")
        lines = (reading / "checks.txt").read_text().splitlines(keepends=True)
        lines = [x.split("; the elastic net withdrawn:")[0] + "\n" if x.startswith("PASS species:") else x
                 for x in lines]
        (reading / "checks.txt").write_text("".join(lines))
        self.refused(reading, "1 reading", "withdrawn")

    def test_the_second_pass_count_goes_to_the_page(self):
        reading = self.reading("withdrawn")
        text = (reading / "checks.txt").read_text().replace("0 species ran a second pass, 0 second-pass",
                                                            "2 species ran a second pass, 1 second-pass", 1)
        (reading / "checks.txt").write_text(text)
        out = self.write(reading, "withdrawn")
        self.assertEqual([r["second_passes"] for r in rows_of(out / "page-withdrawn.tsv")], ["2"])
        self.assertIn("after 2 species had run a second pass with it", (out / "README.md").read_text())
        self.assertEqual([r["second_passes"] for r in rows_of(self.out("withdrawn") / "page-withdrawn.tsv")], ["0"])

    def test_a_first_copy_set_aside_is_never_read_as_a_species(self):
        reading = self.reading()
        shutil.copytree(reading / "score" / "synthia_coli", reading / "score" / "synthia_coli.20261123T080000Z")
        with self.assertRaisesRegex(page.Stop, "^gate 2 species: .*a first copy of synthia_coli set aside"):
            self.write(reading)

    def test_a_withdrawal_time_that_is_no_utc_time_is_refused(self):
        reading = self.reading("withdrawn")
        text = (reading / "withdrawn.tsv").read_text()
        (reading / "withdrawn.tsv").write_text(text.replace(F.WITHDRAWN["elasticnet"][1], "23 Nov"))
        self.refused(reading, "1 reading", "withdrawn")

    def test_a_row_of_another_species_in_a_species_table_is_refused(self):
        reading = self.reading()
        other = rows_of(reading / "score" / "fictus_aureus" / "datasets.tsv")[0]
        edit(reading / "score" / "synthia_coli" / "datasets.tsv", lambda rows: rows.append(dict(other)))
        self.refused(reading, "2 species")

    def test_an_antibiotic_named_otherwise_than_in_lists_tsv_is_refused(self):
        reading = self.reading()

        def rename(rows):
            rows[0]["antibiotic"] = "drug z"
        edit(reading / "score" / "synthia_coli" / "datasets.tsv", rename)
        self.refused(reading, "3 M")

    def test_a_done_of_another_species_is_refused(self):
        reading = self.reading()
        done = reading / "run" / "synthia_coli" / "DONE"
        done.write_text(done.read_text().replace("DONE synthia_coli ", "DONE fictus_aureus ", 1))
        self.refused(reading, "9 passes")

    def test_a_constant_mark_the_report_does_not_hold_is_refused(self):
        reading = self.reading()
        folder = reading / "score" / "fictus_aureus"
        arm, dataset, split, seed = F.CONSTANT

        def unmark(rows):
            for r in rows:
                if (r["arm"], r["dataset"], r["split"], r["seed"]) == (arm, dataset, split, str(seed)):
                    r["constant"] = "0"
        edit(folder / "draws.tsv", unmark)
        draws = [page.typed_draw(r) for r in rows_of(folder / "draws.tsv")]
        F.C9.write_tsv(folder / "datasets.tsv", F.score.dataset_rows(draws), F.score.DATASET_COLUMNS)
        with self.assertRaisesRegex(page.Stop, "^gate 13 reports: .* otherwise than its report"):
            self.write(reading)

    def test_a_plans_log_seed_count_other_than_the_draws_is_refused(self):
        reading = self.reading()
        path = reading / "run" / "synthia_coli" / "plans.log"
        path.write_text(re.sub(r"(\d+) dataset seeds", lambda m: f"{int(m[1]) + 1} dataset seeds", path.read_text()))
        self.refused(reading, "14 checks")

    def test_a_libraries_txt_short_of_a_library_is_refused(self):
        reading = self.reading()
        path = reading / "run" / "synthia_coli" / "libraries.txt"
        path.write_text("".join(path.read_text().splitlines(keepends=True)[:-1]))
        self.refused(reading, "15 binary")

    def test_a_species_outside_the_first_run_is_refused(self):
        reading = self.reading()
        for top in ("score", "run"):
            shutil.copytree(reading / top / "synthia_coli", reading / top / "nullia_species")
        self.refused(reading, "2 species")

    def test_a_dataset_of_m_is_never_counted(self):
        reading = self.reading()

        def add_m(rows):
            r = dict(rows[0], dataset="synthia_coli__drug_m", antibiotic="drug m")
            rows.append(r)
        edit(reading / "score" / "synthia_coli" / "datasets.tsv", add_m)
        self.refused(reading, "3 M")

    def test_an_elastic_net_read_against_gbt_is_refused(self):
        reading = self.reading()

        def to_gbt(rows):
            for r in rows:
                if r["arm"] == "elasticnet":
                    r["learner"] = "gbt"
        edit(reading / "score" / "fictus_aureus" / "datasets.tsv", to_gbt)
        self.refused(reading, "4 pairs")

    def test_an_arm_with_no_default_at_its_seeds_is_refused(self):
        reading = self.reading()
        edit(reading / "score" / "fictus_aureus" / "datasets.tsv",
             lambda rows: rows.__setitem__(slice(None), [r for r in rows if r["arm"] != page.ARM_BASE]))
        self.refused(reading, "4 pairs")

    def test_a_twin_row_at_a_draw_seed_is_refused(self):
        reading = self.reading()

        def move(rows):
            for r in rows:
                if r["arm"] == "twin" and r["seed"] == "45":
                    r["seed"] = "44"
        edit(reading / "score" / "synthia_coli" / "draws.tsv", move)
        self.refused(reading, "5 twin")

    def test_conditions_that_read_the_twin_as_an_arm_are_refused(self):
        reading = self.reading()
        edit(reading / "count" / "conditions.tsv", lambda rows: rows.append(dict(rows[0], arm="twin")))
        self.refused(reading, "5 twin")

    def test_a_study_run_moved_onto_the_lineage_split_is_refused(self):
        reading = self.reading()

        def move(rows):
            for r in rows:
                if r["split"] == "study" and r["arm"] == "xgboost":
                    r["split"] = "lineage"
        edit(reading / "score" / "synthia_coli" / "draws.tsv", move)
        self.refused(reading, "6 study")

    def test_a_family_short_of_a_pair_is_refused(self):
        reading = self.reading()
        edit(reading / "count" / "count.tsv", lambda rows: rows.__setitem__(slice(None), [
            r for r in rows if not (r["split"] == "lineage" and r["pair"] == "twin-k31" and r["learner"] == "lasso"
                                    and r["measure"] == "balanced_accuracy")]))
        self.refused(reading, "7 family")

    def test_a_withdrawn_arms_report_is_refused(self):
        reading = self.reading("withdrawn")
        folder = reading / "run" / "synthia_coli" / "outside" / "elasticnet"
        folder.mkdir(parents=True)
        (folder / "synthia_coli__drug_a.lineage-s42.outside.json").write_text("{}")
        self.refused(reading, "8 withdrawn", "withdrawn")

    def test_a_withdrawn_pair_off_p_one_is_refused(self):
        reading = self.reading("withdrawn")

        def off(rows):
            for r in rows:
                if r["pair"] == "elasticnet-k31":
                    r["sign_p"] = "0.5"
        edit(reading / "count" / "count.tsv", off)
        self.refused(reading, "8 withdrawn", "withdrawn")

    def test_a_first_pass_alone_never_holds_an_elastic_net_result(self):
        reading = self.reading()
        done = reading / "run" / "madeup_pneumoniae" / "DONE"
        done.write_text(done.read_text().replace("elasticnet run", "elasticnet leave-out"))
        self.refused(reading, "9 passes")

    def test_a_first_pass_alone_is_refused_unless_the_elastic_net_is_withdrawn(self):
        reading = self.reading()
        slug = "madeup_pneumoniae"
        done = reading / "run" / slug / "DONE"
        done.write_text(done.read_text().replace("elasticnet run", "elasticnet leave-out"))
        shutil.rmtree(reading / "run" / slug / "outside" / "elasticnet")
        for name in ("draws.tsv", "datasets.tsv"):
            edit(reading / "score" / slug / name, lambda rows: rows.__setitem__(slice(None), [
                r for r in rows if r["arm"] != "elasticnet"]))
        self.refused(reading, "9 passes")

    def test_a_dataset_left_out_of_datasets_tsv_is_refused(self):
        reading = self.reading()
        edit(reading / "score" / "synthia_coli" / "datasets.tsv", lambda rows: rows.__setitem__(slice(None), [
            r for r in rows if r["dataset"] != "synthia_coli__drug_b"]))
        self.refused(reading, "10 rows")

    def test_a_missing_draw_with_figures_is_refused(self):
        reading = self.reading()

        def miss(rows):
            for r in rows:
                if r["arm"] == "lightgbm" and r["seed"] == "43" and r["dataset"] == "synthia_coli__drug_b":
                    r["scored"], r["why"] = "0", "draw missing"
        edit(reading / "score" / "synthia_coli" / "draws.tsv", miss)
        self.refused(reading, "11 missing")

    def test_a_dataset_figure_that_is_not_its_draws_mean_is_refused(self):
        reading = self.reading()

        def move(rows):
            for r in rows:
                if r["arm"] == "xgboost" and r["split"] == "lineage" and r["dataset"] == "synthia_coli__drug_b":
                    r["balanced_accuracy"] = "0.5"
        edit(reading / "score" / "synthia_coli" / "datasets.tsv", move)
        self.refused(reading, "12 draws")

    def test_a_report_whose_fold_count_is_not_its_folds_sum_is_refused(self):
        reading = self.reading()
        path = reading / "run" / "synthia_coli" / "outside" / "xgboost" / "synthia_coli__drug_b.lineage-s42.outside.json"
        doc = json.loads(path.read_text())
        doc["one_tree_folds"] = 3
        path.write_text(json.dumps(doc))
        self.refused(reading, "13 reports")

    def test_a_report_with_no_draw_is_refused(self):
        reading = self.reading()
        folder = reading / "run" / "synthia_coli" / "outside" / "xgboost"
        shutil.copyfile(folder / "synthia_coli__drug_b.lineage-s42.outside.json",
                        folder / "synthia_coli__drug_b.lineage-s45.outside.json")
        self.refused(reading, "13 reports")

    def test_a_report_with_an_unlisted_key_is_refused(self):
        reading = self.reading()
        path = reading / "run" / "fictus_aureus" / "outside" / "lightgbm" / "fictus_aureus__drug_c.lineage-s43.outside.json"
        doc = json.loads(path.read_text())
        doc["feature_names"] = []
        path.write_text(json.dumps(doc))
        self.refused(reading, "13 reports")

    def test_a_failed_score_check_is_refused(self):
        reading = self.reading()

        def fail(rows):
            rows[0]["passed"] = "0"
        edit(reading / "score" / "fictus_aureus" / "checks.tsv", fail)
        self.refused(reading, "14 checks")

    def test_a_plans_count_short_of_its_seeds_is_refused(self):
        reading = self.reading()
        path = reading / "run" / "fictus_aureus" / "plans.log"
        path.write_text("plans: 19 of 20 checks pass, 6 dataset seeds\n")
        self.refused(reading, "14 checks")

    def test_a_done_with_another_result_count_is_refused(self):
        reading = self.reading()
        done = reading / "run" / "synthia_coli" / "DONE"
        words = done.read_text().split()
        words[words.index("results") + 1] = "7"
        done.write_text(" ".join(words) + "\n")
        self.refused(reading, "15 binary")

    def test_a_count_figure_score_py_did_not_write_is_refused(self):
        reading = self.reading()

        def move(rows):
            for r in rows:
                if r["split"] == "lineage" and r["pair"] == "xgboost-k31" and r["measure"] == "recall":
                    r["better"], r["same"] = "4", "1"
        edit(reading / "count" / "count.tsv", move)
        self.refused(reading, "17 count")

    def test_a_beside_count_written_as_deciding_is_refused(self):
        data = page.check_reading(self.root / "run" / "reading", page.load(self.lists()))
        tables = page.tables_of(data)
        for r in tables["page-count.tsv"]:
            if r["file"] == "count-no-unconverged.tsv":
                r["decides"] = 1
        with self.assertRaisesRegex(page.Stop, "^gate 16 beside"):
            page.check_tables(data, tables)

    def test_an_out_dir_that_is_not_empty_is_refused(self):
        reading = self.reading()
        (reading.parent / "out").mkdir()
        (reading.parent / "out" / "x").write_text("x")
        with self.assertRaisesRegex(page.Stop, "not empty"):
            self.write(reading)


class Tables(Base):
    def setUp(self):
        # Written for each test, so that a table that fails its gate fails the test that reads it.
        self.run_out = self.out("run")
        self.w_out = self.out("withdrawn")

    def table(self, name, which="run"):
        return rows_of((self.run_out if which == "run" else self.w_out) / name)

    def test_the_count_files_are_score_pys_own_copied_and_summed(self):
        sums = dict(reversed(line.split("  ", 1)) for line in (self.run_out / "SHA256SUMS").read_text().splitlines())
        for name in page.COUNT_FILES:
            data = (self.root / "run" / "reading" / "count" / name).read_bytes()
            self.assertEqual((self.run_out / "count" / name).read_bytes(), data)
            self.assertEqual(sums[f"count/{name}"], hashlib.sha256(data).hexdigest())
        files = sorted(p.relative_to(self.run_out).as_posix() for p in self.run_out.rglob("*") if p.is_file())
        self.assertEqual(sorted(list(sums) + ["SHA256SUMS"]), sorted(files))
        self.assertIn("score.py's own tables", (self.run_out / "README.md").read_text())

    def test_the_species_table_counts_on_the_primary_measures_orientation(self):
        data = types.SimpleNamespace(species=["synthia_coli"], species_name={"synthia_coli": "Synthia coli"},
                                     withdrawn={})
        rows = [{"dataset": f"synthia_coli__drug_{x}", "split": "lineage", "pair": "xgboost-gbt", "difference": d}
                for x, d in (("a", 0.05), ("b", 0.04), ("c", -0.03), ("d", 0.01))]
        got = [r for r in page.page_species(data, rows) if r["split"] == "lineage" and r["pair"] == "xgboost-gbt"]
        self.assertEqual([(r["better"], r["worse"], r["same"]) for r in got], [(2, 1, 1)])

    def test_score_pys_conditions_are_read_off_the_deciding_count(self):
        paths = sorted((self.root / "run" / "reading" / "score").glob("*/datasets.tsv"))
        rows = [r for p in paths for r in F.C9.read_tsv(p)]
        self.assertNotEqual(F.score.count(rows, []), F.score.count(rows, [], keep_constant=False))
        with tempfile.TemporaryDirectory() as temp, \
                mock.patch.object(F.score, "conditions", wraps=F.score.conditions) as seen, \
                contextlib.redirect_stdout(io.StringIO()):
            F.score.main(["--count", str(Path(temp) / "count")] + [str(p) for p in paths])
        self.assertEqual(seen.call_count, 1)
        self.assertEqual(seen.call_args[0][0], F.score.count(rows, []))

    def test_the_count_rows_decide_for_count_tsv_alone(self):
        rows = self.table("page-count.tsv")
        self.assertEqual({r["file"] for r in rows if r["decides"] == "1"}, {"count.tsv"})
        self.assertEqual({r["file"] for r in rows if r["decides"] == "0"},
                         {"count-no-constant.tsv", "count-no-unconverged.tsv"})

    def test_the_pairs_are_named_with_their_learner_and_what_they_are_read_against(self):
        rows = [r for r in self.table("page-count.tsv") if r["file"] == "count.tsv" and r["split"] == "lineage"
                and r["measure"] == "balanced_accuracy"]
        self.assertEqual([(r["pair"], r["against"]) for r in rows],
                         [("xgboost-gbt", "k31-at-arm-seeds"), ("lightgbm-gbt", "k31-at-arm-seeds"),
                          ("elasticnet-lasso", "k31-at-arm-seeds"), ("twin-gbt", "k31"), ("twin-lasso", "k31")])
        self.assertEqual([r["meets_1"] for r in rows], ["1", "0", "0", "0", "0"])

    def test_conditions_one_and_two_stand_on_the_lineage_primary_rows_only(self):
        for r in self.table("page-count.tsv"):
            on = r["file"] == "count.tsv" and r["split"] == "lineage" and r["measure"] == "balanced_accuracy"
            self.assertEqual(r["meets_1"] != "", on)

    def test_the_beside_condition_decides_nothing(self):
        rows = self.table("page-conditions.tsv")
        self.assertEqual(len(rows), 3 * 7)
        beside = [r for r in rows if r["decides"] == "0"]
        self.assertEqual([(r["arm"], r["condition"]) for r in beside],
                         [("xgboost", "3"), ("lightgbm", "3"), ("elasticnet", "3")])
        self.assertTrue(all(r["reads"].startswith("beside") for r in beside))

    def test_every_species_keeps_a_row_per_split_and_pair(self):
        rows = self.table("page-species.tsv")
        self.assertEqual(len(rows), 3 * (5 + 3))
        self.assertTrue(all(r["decides"] == "0" for r in rows))
        self.assertEqual({r["species"] for r in rows}, {"Synthia coli", "Fictus aureus", "Madeup pneumoniae"})
        empty = [r for r in rows if r["species_slug"] == "madeup_pneumoniae" and r["split"] == "study"]
        self.assertEqual([r["datasets"] for r in empty], ["0", "0", "0"])

    def test_a_study_row_stays_on_the_study_split(self):
        rows = {(r["species_slug"], r["split"], r["pair"]): r for r in self.table("page-species.tsv")}
        self.assertEqual(rows[("synthia_coli", "lineage", "xgboost-gbt")]["datasets"], "2")
        self.assertEqual(rows[("synthia_coli", "study", "xgboost-gbt")]["datasets"], "1")
        self.assertEqual(rows[("fictus_aureus", "study", "xgboost-gbt")]["not_scored"], "1")

    def test_a_dataset_row_names_its_species_antibiotic_and_genomes(self):
        rows = self.table("page-datasets.tsv")
        r = next(r for r in rows if r["dataset"] == "synthia_coli__drug_b" and r["pair"] == "xgboost-gbt")
        self.assertEqual((r["species"], r["antibiotic"], r["genomes"], r["arm_draws_scored"], r["draws"]),
                         ("Synthia coli", "drug b", "120", "3", "3"))
        self.assertAlmostEqual(float(r["difference"]), float(r["arm_value"]) - float(r["default_value"]), places=12)
        self.assertNotIn("synthia_coli__drug_m", {r["dataset"] for r in rows})

    def test_a_dataset_not_scored_holds_no_difference(self):
        rows = [r for r in self.table("page-datasets.tsv") if r["dataset"] == F.UNAVAILABLE and r["split"] == "study"]
        self.assertEqual(len(rows), 3)
        for r in rows:
            self.assertEqual((r["arm_value"], r["difference"], r["recall_difference"]), ("", "", ""))

    def test_an_arm_is_read_against_the_default_at_its_own_seeds(self):
        data = page.check_reading(self.root / "run" / "reading", page.load(self.lists()))
        base = {}
        for r in data.score["synthia_coli"]["datasets.tsv"]:
            if r["arm"] == "k31":
                r["balanced_accuracy"] = "0.123"
            if r["arm"] == page.ARM_BASE:
                base[(r["dataset"], r["split"], r["learner"])] = r["balanced_accuracy"]
        for r in page.page_datasets(data):
            if r["dataset"].startswith("synthia_coli") and r["pair"] != "twin-gbt" and r["pair"] != "twin-lasso" \
                    and r["default_value"] is not None:
                self.assertEqual(r["default_value"], float(base[(r["dataset"], r["split"], r["learner"])]))

    def test_the_fold_sums_carry_their_denominators(self):
        rows = {(r["arm"], r["split"]): r for r in self.table("page-folds.tsv")}
        self.assertEqual((rows[("xgboost", "lineage")]["one_tree_folds"], rows[("xgboost", "lineage")]["folds"]),
                         ("1", "75"))
        self.assertEqual(rows[("lightgbm", "lineage")]["constant_folds"], "1")
        self.assertEqual((rows[("elasticnet", "lineage")]["unconverged_refits"],
                          rows[("elasticnet", "lineage")]["reports_unconverged"]), ("1", "1"))
        self.assertEqual(rows[("xgboost", "study")]["reports_unavailable"], "1")

    def test_one_library_row_when_every_species_agrees(self):
        rows = self.table("page-libraries.tsv")
        self.assertEqual(len(rows), 1)
        self.assertEqual((rows[0]["xgboost"], rows[0]["lightgbm"], rows[0]["sklearn"], rows[0]["species_count"]),
                         ("3.0.5", "4.6.0", "1.7.2", "3"))
        self.assertIn("One binary", (self.run_out / "README.md").read_text())

    def test_versions_that_differ_are_named_with_each_species(self):
        reading = self.reading()
        path = reading / "run" / "fictus_aureus" / "libraries.txt"
        path.write_text(path.read_text().replace("lightgbm 4.6.0", "lightgbm 4.5.0"))
        out = self.write(reading)
        rows = rows_of(out / "page-libraries.tsv")
        self.assertEqual([(r["species"], r["lightgbm"]) for r in rows],
                         [("Fictus aureus", "4.5.0"), ("Madeup pneumoniae", "4.6.0"), ("Synthia coli", "4.6.0")])

    def test_two_binaries_are_named_with_each_species(self):
        reading = self.reading()
        path = reading / "run" / "fictus_aureus" / "DONE"
        path.write_text(path.read_text().replace(F.CODE, hashlib.sha256(b"SYNTHETIC other").hexdigest()))
        out = self.write(reading)
        self.assertEqual(len({r["code"] for r in rows_of(out / "page-binaries.tsv")}), 2)
        self.assertIn("2 binaries", (out / "README.md").read_text())
        self.assertIn("Fictus aureus", (out / "README.md").read_text())

    def test_the_checks_table_holds_score_py_and_plans_py_checks(self):
        rows = {(r["source"], r["check"]): r for r in self.table("page-checks.tsv")}
        self.assertEqual(rows[("plans.py", "the plan's folds by genome are the default's")]["made"], "15")
        self.assertEqual(rows[("plans.py", "the input fingerprint is the release report's")]["made"], "5")
        self.assertEqual(rows[("plans.py", "the export at seed 43 is read")]["passed"], "5")
        self.assertTrue(all(r["passed"] == r["made"] for r in rows.values()))

    def test_the_fit_times_are_sums_and_medians_per_arm(self):
        rows = {r["arm"]: r for r in self.table("page-times.tsv")}
        self.assertEqual((rows["xgboost"]["fits"], rows["xgboost"]["reports"]), ("17", "17"))
        self.assertEqual((rows["elasticnet"]["wall_s_sum"], rows["elasticnet"]["cpu_s_median"]), ("680.0", "39.0"))
        self.assertEqual(set(self.table("page-times.tsv")[0]), set(page.TIME_COLUMNS))

    def test_a_withdrawn_arm_keeps_its_pair_at_p_one_with_its_reason_and_time(self):
        rows = [r for r in self.table("page-count.tsv", "withdrawn") if r["pair"] == "elasticnet-lasso"]
        self.assertTrue(rows)
        self.assertTrue(all((r["withdrawn"], r["entered"], r["sign_p"]) == ("1", "0", "1.0") for r in rows))
        w = self.table("page-withdrawn.tsv", "withdrawn")
        self.assertEqual([(r["arm"], r["reason"], r["time"]) for r in w],
                         [("elasticnet", F.WITHDRAWN["elasticnet"][0], F.WITHDRAWN["elasticnet"][1])])
        self.assertEqual(self.table("page-withdrawn.tsv"), [])
        folds = [r for r in self.table("page-folds.tsv", "withdrawn") if r["arm"] == "elasticnet"]
        self.assertEqual([(r["withdrawn"], r["reports"]) for r in folds], [("1", "0"), ("1", "0")])

    def test_a_withdrawn_arm_keeps_every_table_its_shape(self):
        for name, columns in page.TABLE_COLUMNS.items():
            if name in ("page-withdrawn.tsv", "page-datasets.tsv", "page-binaries.tsv"):
                continue
            a, w = self.table(name), self.table(name, "withdrawn")
            self.assertEqual(len(a), len(w), name)
        printed = page.tables(self.w_out)
        self.assertIn("withdrawn", "\n".join(printed["count"]))
        self.assertIn("withdrawn", "\n".join(printed["folds"]))

    def test_no_genome_id_or_feature_name_is_written(self):
        ids = {g for slug, _, _ in F.SPECIES for g in F.genome_ids(slug)}
        for out in (self.run_out, self.w_out):
            for p in out.rglob("*"):
                if p.is_file():
                    text = p.read_text()
                    self.assertFalse({m.group(0) for m in F.checks.ID_TOKEN.finditer(text)} & ids, p.name)
                    self.assertFalse(any(pattern.search(text) for pattern in F.checks.FEATURE_TOKENS), p.name)

    def test_the_page_tables_print_from_the_out_dir_alone(self):
        printed = page.tables(self.run_out)
        self.assertEqual(tuple(printed), page.PRINTS)
        self.assertIn("*Synthia coli*", "\n".join(printed["species"]))
        self.assertIn("| `xgboost-gbt` | `k31-at-arm-seeds` | 5 of 5 |", "\n".join(printed["count"]))


if __name__ == "__main__":
    unittest.main()
