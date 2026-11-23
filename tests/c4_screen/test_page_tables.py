"""scripts/c4-screen/page_tables.py: the gates and tables of docs/screen-arms-result.md, on SYNTHETIC readings.

Every species, antibiotic, dataset, figure and binary here is made up in the test (tests/c4_screen/page_fixture.py),
and the count tables are score.py --count's own output on them. One test reads the committed list
scripts/c4-screen/lists.tsv, which holds no score.
"""

import contextlib
import hashlib
import io
import shutil
import tempfile
import unittest
from pathlib import Path

from tests.c4_screen import common
from tests.c4_screen import page_fixture as fx

page = common.load("page_tables")
FILES = ["README.md", "SHA256SUMS", "checks.txt", "count/conditions.tsv", "count/count-floor-all.tsv",
         "count/count-no-constant.tsv", "count/count.tsv", "count/withdrawn.txt", "gates.txt", "page-beside.tsv",
         "page-binaries.tsv", "page-checks.tsv", "page-conditions.tsv", "page-count.tsv", "page-datasets.tsv",
         "page-fit-times.tsv", "page-species.tsv", "reading.txt"]


class Built(unittest.TestCase):
    """One SYNTHETIC reading per class, written once; each test copies what it changes."""

    build = {}

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.root = fx.build(Path(cls.tmp.name) / "base", **cls.build)
        cls.n = 0

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def fresh(self):
        type(self).n += 1
        return fx.copy(self.root, Path(self.tmp.name) / f"copy{self.n}")

    def write(self, root):
        out = root / f"out{self.n}"
        return page.write(root / "reading", root / "lists.tsv", root / "checks.txt", out), out

    def stops(self, root, n):
        name = dict((g[0], g[1]) for g in page.GATES)[n]
        with self.assertRaisesRegex(page.Stop, f"^gate {n} {name}: "):
            page.check_reading(root / "reading", fx.read(root / "lists.tsv"), root / "checks.txt")


class Tables(Built):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.out = cls.root / "out"
        cls.files = page.write(cls.root / "reading", cls.root / "lists.tsv", cls.root / "checks.txt", cls.out)

    def table(self, name):
        return fx.read(self.out / name)

    def test_writes_each_file_once_with_its_sha256(self):
        self.assertEqual(sorted(self.files + ["SHA256SUMS"]), FILES)
        lines = (self.out / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
        self.assertEqual([line.split("  ", 1)[1] for line in lines], [f for f in FILES if f != "SHA256SUMS"])
        for line in lines:
            sha, name = line.split("  ", 1)
            self.assertEqual(hashlib.sha256((self.out / name).read_bytes()).hexdigest(), sha, name)
        for name in page.COUNT_FILES + ("withdrawn.txt",):
            self.assertEqual((self.out / "count" / name).read_bytes(),
                             (self.root / "reading" / "count" / name).read_bytes(), name)

    def test_gates_txt_names_the_sixteen_gates_in_order(self):
        lines = (self.out / "gates.txt").read_text(encoding="utf-8").splitlines()
        self.assertEqual([int(line.split(" ", 1)[0]) for line in lines], list(range(1, 17)))
        self.assertEqual(len({line.split(" ")[1] for line in lines}), 16)

    def test_refuses_an_out_dir_that_is_not_empty(self):
        root = self.fresh()
        (root / "taken").mkdir()
        (root / "taken" / "x").write_text("x", encoding="utf-8")
        with self.assertRaisesRegex(page.Stop, "is not empty"):
            page.write(root / "reading", root / "lists.tsv", root / "checks.txt", root / "taken")

    def test_the_count_decides_and_the_family_is_the_21_lineage_pairs(self):
        rows = self.table("page-count.tsv")
        self.assertTrue(rows and all(r["decides"] == "1" for r in rows))
        family = [r for r in rows if r["family"] == "1"]
        self.assertEqual(len(family), 21)
        self.assertEqual(sorted({(r["pair"], r["set"]) for r in family}),
                         sorted([("cmh-16000-k31", "all"), ("cmh-4000-k31", "all"), ("floor-3-k31", "800-or-more"),
                                 ("twin-k31", "all"), ("twin-k31", "800-or-more")]))
        for r in family:
            want = ("1", "1") if r["pair"] == "cmh-4000-k31" else ("0", "0")
            self.assertEqual((r["meets_1"], r["meets_2"]), want, (r["pair"], r["learner"]))
        self.assertTrue(all(r["meets_1"] == r["meets_2"] == "" for r in rows if r["family"] == "0"))

    def test_a_screen_arm_is_counted_on_rf_gbt_and_lasso_only(self):
        for name in ("page-count.tsv", "page-beside.tsv"):
            for r in self.table(name):
                if r["pair"].startswith("cmh-"):
                    self.assertIn(r["learner"], ("rf", "gbt", "lasso"), name)
        cmh = {r["learner"] for r in self.table("page-datasets.tsv") if r["pair"].startswith("cmh-")}
        self.assertEqual(cmh, {"rf", "gbt", "lasso"})

    def test_the_floor_arm_decides_over_800_genomes_or_more_and_floor_all_stands_beside(self):
        floor = [r for r in self.table("page-count.tsv") if r["pair"] == "floor-3-k31"]
        self.assertTrue(floor and all(r["set"] == "800-or-more" for r in floor))
        self.assertEqual({(r["split"], r["datasets"]) for r in floor}, {("lineage", "6"), ("study", "3")})
        beside = self.table("page-beside.tsv")
        self.assertTrue(beside and all(r["decides"] == "0" for r in beside))
        whole = [r for r in beside if r["beside"] == "floor-all"]
        self.assertTrue(whole and all(r["set"] == "all" and r["holm_p"] == "" for r in whole))
        self.assertEqual({(r["split"], r["datasets"]) for r in whole}, {("lineage", "12"), ("study", "6")})

    def test_the_stricter_reading_decides_nothing_and_no_twin_arm_is_read(self):
        rows = self.table("page-conditions.tsv")
        self.assertEqual({r["arm"] for r in rows}, {"cmh-4000", "cmh-16000", "floor-3"})
        for r in rows:
            self.assertEqual(r["decides"], "0" if r["reading"].startswith("stricter") else "1", r["reading"])
        self.assertEqual({r["reading"] for r in rows if r["decides"] == "1"},
                         {"the learners it touches, every dataset", "ruling 15, 800 genomes or more"})
        met = {(r["arm"], r["condition"]): r["met"] for r in rows if r["decides"] == "1"}
        self.assertEqual([met[("cmh-4000", c)] for c in "12345"], ["1"] * 5)
        self.assertEqual(met[("floor-3", "1 and 2")], "0")

    def test_species_are_named_in_full_and_their_counts_decide_nothing(self):
        rows = self.table("page-species.tsv")
        self.assertEqual({r["species"] for r in rows}, {"Synthia alpha", "Synthia beta", "Synthia gamma"})
        self.assertTrue(all(r["decides"] == "0" for r in rows))
        # Per species: 21 lineage pairs and learners, 11 on the study split.
        self.assertEqual(len(rows), 3 * (21 + 11))
        keys = [(r["species_slug"], r["split"], r["pair"], r["set"], r["learner"]) for r in rows]
        self.assertEqual(len(set(keys)), len(keys))
        beta = {(r["split"], r["pair"], r["learner"]): r for r in rows if r["species"] == "Synthia beta"
                and r["set"] == "all"}
        r = beta[("lineage", "cmh-4000-k31", "rf")]
        self.assertEqual((r["datasets"], r["entered"], r["not_scored"], r["better"]), ("4", "3", "1", "3"))

    def test_a_missing_draw_is_never_read_as_zero(self):
        rows = [r for r in self.table("page-datasets.tsv") if r["dataset"] == fx.UNSCORED and r["split"] == "lineage"]
        self.assertTrue(rows)
        for r in rows:
            self.assertEqual((r["value"], r["default_value"], r["difference"], r["recall_difference"]),
                             ("", "", "", ""), r["pair"])
            self.assertEqual(r["draws_scored"], "0")

    def test_a_refused_study_run_is_kept_and_marked(self):
        d, arm = fx.REFUSED
        rows = {(r["dataset"], r["split"], r["pair"], r["learner"]): r for r in self.table("page-datasets.tsv")}
        for learner in ("rf", "gbt", "lasso"):
            r = rows[(d, "study", f"{arm}-k31", learner)]
            self.assertEqual(r["refused"], "1")
            self.assertNotEqual(r["value"], "")
        self.assertEqual({r["refused"] for k, r in rows.items() if k[:3] != (d, "study", f"{arm}-k31")}, {"0"})
        count = {(r["split"], r["pair"], r["learner"], r["measure"]): r for r in self.table("page-count.tsv")}
        self.assertEqual(count[("study", f"{arm}-k31", "rf", "balanced_accuracy")]["constant"], "1")

    def test_each_dataset_names_its_species_and_antibiotic_and_says_what_is_counted(self):
        rows = self.table("page-datasets.tsv")
        self.assertNotIn(fx.M_DATASET, {r["dataset"] for r in rows})
        for r in rows:
            lrow = {x["dataset"]: x for x in fx.list_rows()}[r["dataset"]]
            self.assertEqual((r["species"], r["antibiotic"], r["genomes"]),
                             (lrow["species"], lrow["antibiotic"], str(lrow["samples"])))
            moves = int(lrow["samples"]) >= 800
            self.assertEqual(r["floor_moves"], str(int(moves)))
            self.assertEqual(r["counted"], "0" if r["pair"] == "floor-3-k31" and not moves else "1")
        twin = {r["split"] for r in rows if r["pair"] == "twin-k31"}
        self.assertEqual(twin, {"lineage"})

    def test_the_equality_checks_are_counted_passed_of_made(self):
        rows = {r["check"]: r for r in self.table("page-checks.tsv")}
        screen = rows["a screen arm's scm and cart figures equal the default's"]
        floor = rows["below 800 genomes the floor arm's figures equal the default's"]
        # scm and cart of two screen arms on 12 lineage and 6 study datasets, less the refused run's two.
        self.assertEqual((screen["passed"], screen["made"]), ("70", "70"))
        # Five learners on 6 lineage datasets under 800 genomes and 3 study ones.
        self.assertEqual((floor["passed"], floor["made"]), ("45", "45"))
        self.assertEqual(rows["exists"]["source"], "score.py checks.tsv")

    def test_fit_times_are_sums_and_medians_per_arm_and_learner(self):
        rows = {(r["arm"], r["learner"]): r for r in self.table("page-fit-times.tsv")}
        self.assertEqual(len(rows), 15)
        r = rows[("cmh-4000", "rf")]
        self.assertEqual((r["fits"], r["seconds_sum"], r["seconds_median"], r["cpu_seconds_sum"], r["peak_kb_median"]),
                         ("12", "18.0", "1.5", "36.0", "2048.0"))
        self.assertEqual(fx.header(self.out / "page-fit-times.tsv"), page.FIT_COLUMNS)

    def test_one_binary_is_named_once(self):
        readme = (self.out / "README.md").read_text(encoding="utf-8")
        self.assertIn(f"All 3 species ran one binary by its code sha256, `{fx.CODE}`", readme)
        self.assertNotIn("{", readme)

    def test_print_names_species_and_antibiotics_in_full(self):
        printed = page.tables(self.out)
        self.assertEqual(list(printed), ["conditions", "family", "count", "measures", "beside", "checks", "fit-times",
                                         "binaries", "species", "datasets"])
        text = "\n".join(line for lines in printed.values() for line in lines)
        for name in ("*Synthia alpha*", "*Synthia beta*", "*Synthia gamma*", "synthetic drug 1"):
            self.assertIn(name, text)
        self.assertNotIn("synthia_", text)
        with contextlib.redirect_stdout(io.StringIO()) as out:
            self.assertEqual(page.main(["print", str(self.out)]), 0)
        self.assertIn("[datasets]", out.getvalue())


class Gates(Built):
    def test_the_clean_reading_passes_every_gate(self):
        root = self.fresh()
        with contextlib.redirect_stdout(io.StringIO()) as out:
            code = page.main(["write", str(root / "reading"), str(root / "lists.tsv"), str(root / "checks.txt"),
                              str(root / "out")])
        self.assertEqual(code, 0)
        self.assertEqual(out.getvalue(), "gates: 16 of 16 pass; wrote 18 files\n")

    def test_a_score_log_in_the_reading_is_refused(self):
        root = self.fresh()
        (root / "reading" / "score" / "synthia_alpha" / "score.log").write_text("checks\n", encoding="utf-8")
        self.stops(root, 1)

    def test_a_run_check_that_failed_is_refused(self):
        root = self.fresh()
        (root / "checks.txt").write_text("FAIL every species' SHA256SUMS\nCHECKS 0 of 1\n", encoding="utf-8")
        self.stops(root, 1)

    def test_a_species_missing_from_the_reading_is_refused(self):
        root = self.fresh()
        for f in (root / "reading" / "score" / "synthia_gamma").iterdir():
            f.unlink()
        (root / "reading" / "score" / "synthia_gamma").rmdir()
        self.stops(root, 2)

    def test_a_dataset_of_m_in_the_reading_is_refused(self):
        root = self.fresh()
        for name in ("datasets.tsv", "draws.tsv"):
            path = root / "reading" / "score" / "synthia_alpha" / name
            rows = fx.read(path)
            more = [dict(r, dataset=fx.M_DATASET, antibiotic="synthetic drug m", genomes="500") for r in rows
                    if r["dataset"] == "synthia_alpha__synthetic_drug_1"]
            fx.rewrite(path, rows + more, fx.header(path))
        self.stops(root, 3)

    def test_a_study_row_on_the_lineage_split_is_refused(self):
        root = self.fresh()
        path = root / "reading" / "score" / "synthia_alpha" / "draws.tsv"
        rows = fx.read(path)
        k = next(i for i, r in enumerate(rows) if r["role"] == "study")
        rows[k]["split"] = "lineage"
        fx.rewrite(path, rows, fx.header(path))
        self.stops(root, 4)

    def test_the_twin_read_as_an_arm_is_refused(self):
        root = self.fresh()
        path = root / "reading" / "count" / "conditions.tsv"
        rows = fx.read(path)
        fx.rewrite(path, rows + [dict(rows[0], arm="twin")], fx.header(path))
        self.stops(root, 5)

    def test_a_draw_not_scored_with_a_figure_is_refused(self):
        root = self.fresh()
        path = root / "reading" / "score" / "synthia_beta" / "draws.tsv"
        rows = fx.read(path)
        k = next(i for i, r in enumerate(rows) if r["scored"] == "0")
        rows[k]["balanced_accuracy"] = "0.0"
        fx.rewrite(path, rows, fx.header(path))
        self.stops(root, 6)

    def test_a_constant_model_draw_not_marked_is_refused(self):
        root = self.fresh()
        path = root / "reading" / "score" / "synthia_gamma" / "draws.tsv"
        rows = fx.read(path)
        for r in rows:
            if r["constant"] == "1":
                r["constant"] = "0"
        fx.rewrite(path, rows, fx.header(path))
        self.stops(root, 7)

    def test_a_dataset_figure_that_is_not_the_mean_of_its_draws_is_refused(self):
        root = self.fresh()
        path = root / "reading" / "score" / "synthia_alpha" / "datasets.tsv"
        rows = fx.read(path)
        k = next(i for i, r in enumerate(rows) if (r["arm"], r["learner"], r["split"]) == ("cmh-4000", "rf", "lineage"))
        rows[k]["balanced_accuracy"] = "0.875"
        fx.rewrite(path, rows, fx.header(path))
        self.stops(root, 8)

    def test_a_check_of_the_run_that_failed_is_refused(self):
        root = self.fresh()
        path = root / "reading" / "score" / "synthia_beta" / "checks.tsv"
        rows = fx.read(path)
        rows[0]["passed"] = "0"
        fx.rewrite(path, rows, fx.header(path))
        self.stops(root, 9)

    def test_a_count_row_of_a_screen_arm_on_scm_is_refused(self):
        root = self.fresh()
        path = root / "reading" / "count" / "count.tsv"
        rows = fx.read(path)
        extra = next(r for r in rows if (r["pair"], r["learner"], r["measure"]) == ("cmh-4000-k31", "rf", "recall"))
        fx.rewrite(path, rows + [dict(extra, learner="scm")], fx.header(path))
        self.stops(root, 10)

    def test_the_floor_arm_counted_over_every_dataset_is_refused(self):
        root = self.fresh()
        path = root / "reading" / "count" / "count.tsv"
        rows = fx.read(path)
        for r in rows:
            if r["pair"] == "floor-3-k31":
                r["set"] = "all"
                r["datasets"] = "12" if r["split"] == "lineage" else "6"
        fx.rewrite(path, rows, fx.header(path))
        self.stops(root, 11)

    def test_a_binary_with_no_code_sha256_is_refused(self):
        root = self.fresh()
        (root / "reading" / "score" / "synthia_beta" / "binary.txt").write_text(
            f"commit {fx.COMMIT}\ngrm {fx.GRM}\n", encoding="utf-8")
        self.stops(root, 13)

    def test_a_condition_that_the_count_does_not_bear_is_refused(self):
        root = self.fresh()
        path = root / "reading" / "count" / "conditions.tsv"
        rows = fx.read(path)
        k = next(i for i, r in enumerate(rows) if (r["arm"], r["condition"]) == ("cmh-4000", "1"))
        rows[k]["met"] = "0"
        fx.rewrite(path, rows, fx.header(path))
        self.stops(root, 14)

    def test_a_count_figure_that_score_py_does_not_write_is_refused(self):
        root = self.fresh()
        path = root / "reading" / "count" / "count.tsv"
        rows = fx.read(path)
        k = next(i for i, r in enumerate(rows) if (r["pair"], r["measure"]) == ("twin-k31", "specificity"))
        rows[k]["species"] = "7"
        fx.rewrite(path, rows, fx.header(path))
        self.stops(root, 15)

    def test_fit_times_of_the_twin_are_refused(self):
        root = self.fresh()
        path = root / "reading" / "score" / "synthia_alpha" / "fits.tsv"
        rows = fx.read(path)
        fx.rewrite(path, rows + [dict(rows[0], arm="twin")], fx.header(path))
        self.stops(root, 16)


class GatesHeld(Built):
    """Gates 2, 3, 4, 5 and 16 on a reading a fault in each would let pass (SYNTHETIC)."""

    def test_a_species_folder_the_first_run_does_not_hold_is_refused(self):
        root = self.fresh()
        shutil.copytree(root / "reading" / "score" / "synthia_gamma", root / "reading" / "score" / "synthia_delta")
        self.stops(root, 2)

    def test_a_genome_count_other_than_the_lists_is_refused(self):
        root = self.fresh()
        path = root / "reading" / "score" / "synthia_alpha" / "datasets.tsv"
        rows = fx.read(path)
        for r in rows:
            if r["dataset"] == "synthia_alpha__synthetic_drug_1":
                r["genomes"] = str(int(r["genomes"]) + 1)
        fx.rewrite(path, rows, fx.header(path))
        self.stops(root, 3)

    def test_a_study_run_at_a_draw_seed_is_refused(self):
        root = self.fresh()
        path = root / "reading" / "score" / "synthia_alpha" / "draws.tsv"
        rows = fx.read(path)
        for r in rows:
            if r["role"] == "study":
                r["seed"] = "43"
        fx.rewrite(path, rows, fx.header(path))
        self.stops(root, 4)

    def test_the_twin_at_a_draw_seed_is_refused(self):
        root = self.fresh()
        path = root / "reading" / "score" / "synthia_alpha" / "draws.tsv"
        rows = fx.read(path)
        k = next(i for i, r in enumerate(rows) if r["role"] == "twin")
        rows[k]["seed"] = "44"
        fx.rewrite(path, rows, fx.header(path))
        self.stops(root, 5)

    def test_a_fit_with_no_cpu_time_or_no_peak_is_refused(self):
        for column in ("cpu_seconds", "peak_kb"):
            root = self.fresh()
            path = root / "reading" / "score" / "synthia_alpha" / "fits.tsv"
            rows = fx.read(path)
            rows[0][column] = ""
            fx.rewrite(path, rows, fx.header(path))
            with self.subTest(column=column):
                self.stops(root, 16)


class TieAndFloorSet(unittest.TestCase):
    """The species table's tally and floor set, on made-up differences and lists rows."""

    def test_a_difference_of_0_02_is_the_same_and_more_is_better_or_worse(self):
        self.assertEqual(page.tally([0.02, -0.02, 0.021, -0.021, 0.0]), (1, 1, 3))
        self.assertEqual(page.tally([0.02, -0.021], "log_loss"), (1, 0, 1))

    def test_a_dataset_of_exactly_800_genomes_is_in_the_floor_set(self):
        rows = [{"dataset": "x800", "samples": "800", "floor_moves": "1"},
                {"dataset": "x799", "samples": "799", "floor_moves": "0"}]
        self.assertEqual(page.set_datasets(rows, page.MOVES), ["x800"])
        self.assertEqual(page.set_datasets(rows, page.ALL), ["x800", "x799"])
        self.assertEqual((bool(page.score.floor_moves("800")), bool(page.score.floor_moves("799"))), (True, False))


class TwinSets(Built):
    """Condition 3 reads the twin over the arm's own set: every dataset for a screen arm, 800 or more for the floor."""

    def test_each_arm_reads_the_twin_over_its_own_set(self):
        rows = fx.read(self.root / "reading" / "count" / "count.tsv")
        for r in rows:
            if (r["split"], r["pair"], r["measure"]) == ("lineage", "twin-k31", page.PRIMARY):
                r["interval_lo"] = "0.05" if r["set"] == page.MOVES else "-0.05"
                r["holm_p"] = "1.0"
        read = page.read_conditions(rows)
        self.assertFalse(read[(page.FLOOR_ARM, page.FLOOR_READING, "3")])
        self.assertTrue(read[("cmh-4000", page.SCREEN_READING, "3")])


class CartOff(Built):
    build = {"cart_off": True}

    def test_a_screen_arm_whose_cart_differs_from_the_default_is_refused(self):
        self.stops(self.fresh(), 9)


class Withdrawn(Built):
    build = {"withdrawn": ("cmh-16000",)}

    def test_a_withdrawn_arm_keeps_its_pairs_in_the_family_at_p_1(self):
        root = self.fresh()
        files, out = self.write(root)
        family = [r for r in fx.read(out / "page-count.tsv") if r["family"] == "1"]
        self.assertEqual(len(family), 21)
        gone = [r for r in family if r["pair"] == "cmh-16000-k31"]
        self.assertEqual([(r["withdrawn"], r["datasets"], r["sign_p"]) for r in gone], [("1", "0", "1.0")] * 3)
        self.assertNotIn("cmh-16000", {r["pair"][:9] for r in fx.read(out / "page-species.tsv")})
        self.assertNotIn("cmh-16000", {r["arm"] for r in fx.read(out / "page-fit-times.tsv")})
        self.assertEqual((out / "count" / "withdrawn.txt").read_text(encoding="utf-8"), "cmh-16000\n")

    def test_a_withdrawn_arm_dropped_from_the_family_is_refused(self):
        root = self.fresh()
        path = root / "reading" / "count" / "count.tsv"
        rows = [r for r in fx.read(path) if not (r["pair"] == "cmh-16000-k31" and r["holm_p"] != "")]
        fx.rewrite(path, rows, fx.header(path))
        self.stops(root, 12)


class Binaries(Built):
    build = {"codes": {"synthia_beta": "d" * 64}}

    def test_two_binaries_are_named_each_with_its_species(self):
        root = self.fresh()
        files, out = self.write(root)
        readme = (out / "README.md").read_text(encoding="utf-8")
        self.assertIn("The 3 species ran 2 binaries by their code sha256; `page-binaries.tsv` names each binary with "
                      "its species.", readme)
        rows = {r["species"]: r["code"] for r in fx.read(out / "page-binaries.tsv")}
        self.assertEqual(rows, {"Synthia alpha": fx.CODE, "Synthia beta": "d" * 64, "Synthia gamma": fx.CODE})
        printed = "\n".join(page.tables(out)["binaries"])
        self.assertIn("| *Synthia beta* |", printed)


class CommittedList(unittest.TestCase):
    def test_the_first_run_is_233_datasets_in_23_species_none_in_m(self):
        rows = fx.read(common.KIT / "lists.tsv")
        run = page.first_run(rows)
        self.assertEqual((len(run), len({r["species_slug"] for r in run})), (233, 23))
        self.assertEqual(run, [r for r in rows if r["first_run"] == "1"])
        self.assertFalse([r for r in run if r["in_M"] == "1"])
        self.assertEqual(sum(r["floor_moves"] == "1" for r in run), 92)
        self.assertEqual(sum(r["in_S"] == "1" for r in run), 55)


if __name__ == "__main__":
    unittest.main()
