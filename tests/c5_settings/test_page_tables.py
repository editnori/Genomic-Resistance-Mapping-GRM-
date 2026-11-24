"""scripts/c5-settings/page_tables.py: the gates of plan block C5's tables, each able to fail, and the tables.

Every reading here is SYNTHETIC (page_fixture.py); the one real file read is scripts/c5-settings/lists.tsv.
"""

import csv
import tempfile
import unittest
from pathlib import Path

from tests.c5_settings import common
from tests.c5_settings import page_fixture as fx

page = common.load("page_tables")


def tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def put(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


class Case(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.top = Path(self.tmp.name)
        self.out = self.top / "out"

    def tearDown(self):
        self.tmp.cleanup()

    def build(self, **kw):
        self.reading, self.lists, self.checks = fx.build(self.top, **kw)
        return self.reading

    def write(self):
        return page.write(self.reading, self.lists, self.checks, self.out)

    def refused(self, pattern):
        with self.assertRaisesRegex(page.Stop, pattern):
            self.write()


class RealList(unittest.TestCase):
    def test_the_first_run_is_55_datasets_in_8_species_none_in_m(self):
        run = page.first_run(tsv(common.KIT / "lists.tsv"))
        by = {}
        for r in run:
            by[r["species"]] = by.get(r["species"], 0) + 1
        self.assertEqual(by, {"Acinetobacter baumannii": 2, "Enterobacter cloacae": 8, "Escherichia coli": 15,
                              "Klebsiella pneumoniae": 18, "Mycobacterium tuberculosis": 1,
                              "Pseudomonas aeruginosa": 6, "Salmonella enterica": 1, "Staphylococcus aureus": 4})
        self.assertEqual(len(run), 55)
        self.assertFalse([r for r in run if r["in_M"] != "0"])
        self.assertEqual(len([r for r in tsv(common.KIT / "lists.tsv") if r["in_M"] == "1"]), 14)

    def test_every_first_run_dataset_reads_as_a_choices_file_name(self):
        for r in page.first_run(tsv(common.KIT / "lists.tsv")):
            name = f"{r['dataset']}.study-whole-s42.f0.choices.json"
            self.assertTrue(page.CHOICE_NAME.fullmatch(name), name)


class Gates(Case):
    def test_a_clean_reading_passes_every_gate_and_writes_its_files(self):
        self.build()
        files = self.write()
        self.assertIn("count/count.tsv", files)
        for name in ("page-conditions.tsv", "page-count.tsv", "page-beside.tsv", "page-scm.tsv", "page-choices.tsv",
                     "page-trials.tsv", "page-cutoffs.tsv", "page-species.tsv", "page-datasets.tsv",
                     "page-checks.tsv", "page-binaries.tsv", "gates.txt", "README.md", "reading.txt", "checks.txt"):
            self.assertIn(name, files)
        sums = (self.out / "SHA256SUMS").read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(sums), len(files))
        self.assertEqual(len((self.out / "gates.txt").read_text(encoding="utf-8").splitlines()), 16)
        for name in page.COUNT_FILES:
            self.assertEqual((self.out / "count" / name).read_bytes(), (self.reading / "count" / name).read_bytes())

    def test_a_score_log_in_the_reading_is_refused(self):
        self.build()
        (self.reading / "score" / "synthia_alpha" / "score.log").write_text("SYNTHETIC\n", encoding="utf-8")
        self.refused(r"gate 1 files: score.log stays out")

    def test_checks_txt_with_a_failed_check_is_refused(self):
        self.build()
        self.checks.write_text("FAIL SYNTHETIC\nCHECKS 0 of 1\n", encoding="utf-8")
        self.refused(r"gate 1 files")

    def test_a_species_missing_from_the_reading_is_refused(self):
        self.build(species=["synthia_alpha", "synthia_beta"])
        self.refused(r"gate 2 species: score/ holds 2 species of the first run's 3")

    def test_a_dataset_of_m_on_the_list_is_never_counted(self):
        self.build(m_in_run=True)
        self.refused(r"gate 2 species: the species' datasets.tsv files hold 5 datasets, the first run 4")

    def test_a_study_row_on_the_lineage_split_is_refused(self):
        def change(slug, draws):
            for r in draws:
                if (slug, r["split"], r["arm"], r["seed"], r["learner"]) == ("synthia_beta", "lineage", "default", 42,
                                                                               "scm"):
                    r["role"] = "study"
        self.build(change_draws=change)
        self.refused(r"gate 4 study")

    def test_the_twin_read_as_an_arm_is_refused(self):
        self.build()
        path = self.reading / "count" / "conditions.tsv"
        rows = tsv(path)
        put(path, rows + [dict(r, arm="twin") for r in rows if r["arm"] == "search"])
        self.refused(r"gate 5 twin: conditions.tsv reads cutoff, search, twin as the arms")

    def test_the_plain_best_in_a_holm_family_is_refused(self):
        self.build()
        path = self.reading / "count" / "count-plain-best.tsv"
        put(path, [dict(r, holm_p="0.5") for r in tsv(path)])
        self.refused(r"gate 6 plain")

    def test_a_missing_draw_read_as_zero_is_refused(self):
        def change(slug, draws):
            for r in draws:
                if r["why"] == "missing draw":
                    r["balanced_accuracy"] = 0.0
        self.build(missing=True, change_draws=change)
        self.refused(r"gate 7 missing: score/synthia_alpha/draws.tsv: .* a draw not scored has a figure")

    def test_a_dataset_scored_with_a_draw_missing_is_refused(self):
        self.build(missing=True)
        path = self.reading / "score" / "synthia_alpha" / "datasets.tsv"
        put(path, [dict(r, scored="1") if r["scored_draws"] != r["draws"] else r for r in tsv(path)])
        self.refused(r"gate 7 missing")

    def test_an_unmarked_constant_model_row_is_refused(self):
        def change(slug, draws):
            for r in draws:
                if r["constant"] == 1 and r["learner"] == "rf":
                    r["constant"] = 0
        self.build(constant=True, change_draws=change)
        self.refused(r"gate 8 constant")

    def test_a_draw_twice_is_refused(self):
        def change(slug, draws):
            if slug == "synthia_gamma":
                draws.append(dict(draws[0]))
        self.build(change_draws=change)
        self.refused(r"gate 9 draws: score/synthia_gamma/draws.tsv holds other draws .* 1 twice")

    def test_a_failed_check_is_refused(self):
        self.build()
        path = self.reading / "score" / "synthia_beta" / "checks.tsv"
        rows = tsv(path)
        rows[0]["passed"] = "0"
        put(path, rows)
        self.refused(r"gate 10 checks")

    def test_the_cutoff_arms_log_loss_apart_from_the_defaults_is_refused(self):
        def change(slug, draws):
            for r in draws:
                if r["arm"] == "cutoff" and r["learner"] == "cart" and r.get("log_loss") is not None:
                    r["log_loss"] += 1 / 64
        self.build(change_draws=change)
        self.refused(r"gate 11 cutoff: the cutoff arm's log loss equals the default's")

    def test_scms_cutoff_moved_from_half_is_refused(self):
        def change(key, doc):
            doc["learners"][0]["cutoff"] = 0.45
            return doc
        self.build(change_choice=change)
        self.refused(r"gate 11 cutoff: scm's cutoff is 0.5")

    def test_a_choices_file_missing_inside_a_scored_draw_is_refused(self):
        self.build()
        (self.reading / "score" / "synthia_gamma" / "choices"
         / "synthia_gamma__synthetic_drug_4.lineage-s44.f2.choices.json").unlink()
        self.refused(r"gate 12 choices: score/synthia_gamma: synthia_gamma__synthetic_drug_4 lineage seed 44 has "
                     r"choices files for folds \[0, 1, 3, 4\]")

    def test_a_scored_draw_short_of_its_defaults_folds_is_refused(self):
        self.build()
        (self.reading / "score" / "synthia_gamma" / "choices"
         / "synthia_gamma__synthetic_drug_4.lineage-s44.f4.choices.json").unlink()
        self.refused(r"gate 12 choices: score/synthia_gamma: synthia_gamma__synthetic_drug_4 lineage seed 44 has "
                     r"choices files for 4 folds, its default's result 5")

    def test_a_scored_draw_with_no_choices_file_is_refused(self):
        def change(key, doc):
            return None if key[:3] == ("synthia_beta__synthetic_drug_3", "study-whole", 42) else doc
        self.build(change_choice=change)
        self.refused(r"gate 12 choices: .* has no choices file and its arms are not a missing draw")

    def test_a_pick_against_the_default_anchored_rule_is_refused(self):
        def change(key, doc):
            doc["learners"][1]["picked_trial"] = 3
            return doc
        self.build(change_choice=change)
        self.refused(r"gate 12 choices: .* cart: the pick 3 does not follow the default-anchored rule")

    def test_a_setting_in_the_counts_of_what_the_search_chose_is_refused(self):
        row = {"split": "lineage", "learner": "rf", "trial": 7, "picked": {"max_depth": 4321}, "plain_best": 0}
        with self.assertRaisesRegex(page.Stop, r"gate 13 settings: page-trials.tsv: picked holds"):
            page.check_counts_only("page-trials.tsv", [row], page.TRIAL_COLUMNS)

    def test_a_column_beside_the_counts_is_refused(self):
        with self.assertRaisesRegex(page.Stop, r"gate 13 settings"):
            page.check_counts_only("page-trials.tsv", [], page.TRIAL_COLUMNS + ["picked_setting"])

    def test_a_binary_txt_without_its_code_is_refused(self):
        self.build()
        (self.reading / "score" / "synthia_beta" / "binary.txt").write_text(
            f"commit {fx.COMMIT}\ngrm {fx.GRM}\n", encoding="utf-8")
        self.refused(r"gate 14 binary: score/synthia_beta/binary.txt does not name")

    def test_conditions_tsv_apart_from_the_count_is_refused(self):
        self.build()
        path = self.reading / "count" / "conditions.tsv"
        put(path, [dict(r, met=str(1 - int(r["met"])))
                   if (r["arm"], r["reading"], r["condition"]) == ("search", "ruling 15", "3") else r
                   for r in tsv(path)])
        self.refused(r"gate 15 conditions: search, ruling 15, condition 3")

    def test_a_count_that_is_not_score_pys_is_refused(self):
        self.build()
        path = self.reading / "count" / "count.tsv"
        put(path, [dict(r, species="9") if (r["split"], r["measure"]) == ("study", "specificity") else r
                   for r in tsv(path)])
        self.refused(r"gate 16 count: count/count.tsv is not what score.py --count writes")

    def test_an_antibiotic_named_otherwise_than_in_lists_tsv_is_refused(self):
        self.build()
        path = self.reading / "score" / "synthia_beta" / "datasets.tsv"
        put(path, [dict(r, antibiotic="synthetic drug 33") for r in tsv(path)])
        self.refused(r"gate 3 datasets: synthia_beta__synthetic_drug_3 is named otherwise than in lists.tsv")

    def test_a_study_run_at_a_draws_seed_is_refused(self):
        def change(slug, draws):
            for r in draws:
                if slug == "synthia_beta" and r["role"] == "study":
                    r["seed"] = 43
        self.build(change_draws=change)
        self.refused(r"gate 4 study: score/synthia_beta/draws.tsv: synthia_beta__synthetic_drug_3 \S+: a study run at "
                     r"seed 43")

    def test_a_constant_model_marked_on_the_default_is_refused(self):
        def change(slug, draws):
            for r in draws:
                if (slug, r["split"], r["arm"]) == ("synthia_beta", "study", "default"):
                    r["constant"], r["why"] = 1, "constant model: refused for no candidate"
        self.build(constant=True, change_draws=change)
        self.refused(r"gate 8 constant: score/synthia_beta/draws.tsv: synthia_beta__synthetic_drug_3 study default "
                     r"\w+: a constant model marked 1")

    def test_the_cutoff_arms_log_loss_counted_worse_is_refused(self):
        self.build()
        path = self.reading / "count" / "count.tsv"
        put(path, [dict(r, worse="1") if (r["split"], r["pair"], r["measure"]) == ("lineage", "cutoff-default",
                                                                                  "log_loss") else r
                   for r in tsv(path)])
        self.refused(r"gate 11 cutoff: count.tsv counts the cutoff arm's log loss on the lineage split for \w+ as 0 "
                     r"better and 1 worse")

    def test_a_choices_file_of_another_number_of_random_trials_is_refused(self):
        def change(key, doc):
            if key == ("synthia_gamma__synthetic_drug_4", "lineage", 44, 2):
                doc["random_trials"] = 20
            return doc
        self.build(change_choice=change)
        self.refused(r"gate 12 choices: score/synthia_gamma: synthia_gamma__synthetic_drug_4 lineage seed 44 fold 2: "
                     r"the sampler seed, gamma, grid or random trials are not the design's")

    def test_condition_3_is_not_met_when_the_twin_meets_condition_1_in_three_learners(self):
        def rows(twin_meets):
            out = []
            for arm in ("search", "cutoff", "twin"):
                for i, learner in enumerate(fx.LEARNERS):
                    lo = "0.01" if arm == "twin" and i < twin_meets else "-0.01"
                    for measure in ("balanced_accuracy", "recall"):
                        out.append({"split": "lineage", "pair": f"{arm}-default", "learner": learner,
                                    "measure": measure, "interval_lo": lo, "interval_hi": "0.03", "holm_p": "",
                                    "better": "0", "worse": "0"})
            return out
        for twin_meets, met in ((2, True), (3, False)):
            read = page.read_conditions(rows(twin_meets))
            self.assertEqual((read[("search", "ruling 15", "3")], read[("cutoff", "ruling 15", "3")]), (met, met),
                             twin_meets)


class Tables(Case):
    def setUp(self):
        super().setUp()

    def written(self, **kw):
        self.build(**kw)
        self.write()
        return {p.name: tsv(p) for p in self.out.glob("page-*.tsv")}

    def test_ruling_15_decides_and_the_stricter_reading_stands_beside(self):
        rows = self.written()["page-conditions.tsv"]
        self.assertEqual({r["arm"] for r in rows}, {"search", "cutoff"})
        for r in rows:
            self.assertEqual(r["decides"], "1" if r["reading"] == "ruling 15" else "0", r)

    def test_the_cutoff_arms_log_loss_rows_decide_nothing(self):
        rows = self.written()["page-count.tsv"]
        marked = [r for r in rows if r["pair"] == "cutoff-default" and r["measure"] == "log_loss"]
        self.assertEqual(len(marked), 10)
        for r in marked:
            self.assertEqual((r["decides"], r["better"], r["worse"], r["note"]),
                             ("0", "0", "0", page.EQUAL_BY_CONSTRUCTION))
        self.assertTrue(all(r["decides"] == "1" for r in rows if r not in marked))

    def test_holm_runs_over_the_fifteen_lineage_primary_rows(self):
        rows = self.written()["page-count.tsv"]
        family = [r for r in rows if r["family"] == "1"]
        self.assertEqual(len(family), 15)
        self.assertEqual({(r["split"], r["measure"]) for r in family}, {("lineage", "balanced_accuracy")})
        self.assertEqual({r["pair"] for r in family}, {"search-default", "cutoff-default", "twin-default"})
        self.assertFalse([r for r in rows if r["split"] == "study" and r["pair"] == "twin-default"])
        twin = [r for r in rows if r["pair"] == "twin-default"]
        self.assertEqual(len(twin), 30)
        self.assertTrue(all(r["note"] == page.TWIN_NOTE for r in twin))

    def test_the_plain_best_and_the_count_without_constant_models_decide_nothing(self):
        tables = self.written(constant=True)
        rows = tables["page-beside.tsv"]
        self.assertEqual({r["beside"] for r in rows}, {"plain-best", "no-constant"})
        self.assertTrue(all((r["decides"], r["family"], r["holm_p"]) == ("0", "0", "") for r in rows
                            if r["beside"] == "plain-best"))
        self.assertTrue(all(r["decides"] == "0" for r in rows))
        self.assertFalse([r for r in tables["page-count.tsv"] if r["pair"] == "plain-best-default"])

    def test_scm_stands_beside_the_synthetic_loss(self):
        rows = self.written()["page-scm.tsv"]
        self.assertEqual({(r["split"], r["pair"]) for r in rows},
                         {(s, p) for s in ("lineage", "study") for p in ("search-default", "plain-best-default")})
        lin = [r for r in rows if r["split"] == "lineage"]
        self.assertTrue(all((r["synthetic_loss"], r["synthetic_loss_lo"], r["synthetic_loss_hi"])
                            == ("0.0126", "0.0088", "0.0163") for r in lin))
        self.assertTrue(all(r["decides"] == ("1" if r["pair"] == "search-default" else "0") for r in rows))

    def test_what_the_search_chose_is_counted_per_learner(self):
        tables = self.written()
        rows = {(r["split"], r["learner"]): r for r in tables["page-choices.tsv"]}
        # 4 datasets, 3 lineage draws each, 5 folds; 4 study runs, 5 folds.
        self.assertEqual(rows[("lineage", "rf")]["folds"], "60")
        self.assertEqual(rows[("study", "rf")]["folds"], "20")
        self.assertEqual(rows[("lineage", "scm")]["kept_trial_0"], "60")
        self.assertEqual(rows[("lineage", "rf")]["picked_other"], "60")
        self.assertEqual(rows[("lineage", "cart")]["plain_best_differs"], "60")
        self.assertEqual(rows[("lineage", "gbt")]["trial_0_failed"], "60")
        self.assertEqual(rows[("lineage", "lasso")]["every_trial_failed"], "60")
        self.assertEqual(rows[("lineage", "scm")]["cutoff_moved"], "0")
        trials = {(r["split"], r["learner"], r["trial"]): r for r in tables["page-trials.tsv"]}
        self.assertEqual(len(tables["page-trials.tsv"]), 2 * 5 * 31)
        self.assertEqual(trials[("lineage", "rf", "7")]["picked"], "60")
        self.assertEqual(trials[("lineage", "cart", "3")]["plain_best"], "60")
        cut = {(r["split"], r["learner"], r["cutoff"]): r["folds"] for r in tables["page-cutoffs.tsv"]}
        self.assertEqual(len(tables["page-cutoffs.tsv"]), 2 * 5 * 19)
        self.assertEqual(cut[("lineage", "scm", "0.5")], "60")
        self.assertEqual(cut[("study", "cart", "0.45")], "20")

    def test_the_counts_of_what_the_search_chose_hold_no_setting(self):
        self.build()
        self.write()
        for p in self.out.rglob("*"):
            if p.is_file():
                text = p.read_text(encoding="utf-8")
                self.assertNotIn("4321", text, p.name)
                self.assertNotIn("max_depth", text, p.name)

    def test_each_species_is_named_in_full_and_decides_nothing(self):
        rows = self.written()["page-species.tsv"]
        self.assertEqual({r["species"] for r in rows}, {"Synthia alpha", "Synthia beta", "Synthia gamma"})
        self.assertTrue(all(r["decides"] == "0" for r in rows))
        alpha = [r for r in rows if r["species"] == "Synthia alpha" and r["pair"] == "search-default"
                 and r["split"] == "lineage" and r["learner"] == "rf"]
        self.assertEqual([(r["datasets"], r["entered"], r["better"]) for r in alpha], [("2", "2", "2")])

    def test_each_dataset_with_its_species_antibiotic_genomes_and_values(self):
        rows = self.written(constant=True)["page-datasets.tsv"]
        one = [r for r in rows if (r["dataset"], r["split"], r["pair"], r["learner"])
               == ("synthia_alpha__synthetic_drug_2", "lineage", "search-default", "rf")]
        self.assertEqual(len(one), 1)
        r = one[0]
        self.assertEqual((r["species"], r["antibiotic"], r["genomes"]), ("Synthia alpha", "synthetic drug 2", "130"))
        self.assertEqual((r["default_value"], r["value"], r["difference"]), ("0.75", "0.8125", "0.0625"))
        self.assertEqual((r["draws"], r["draws_scored"], r["default_draws_scored"]), ("3", "3", "3"))
        refused = [r for r in rows if r["study_refused"] == "1"]
        self.assertEqual({(r["dataset"], r["split"], r["pair"]) for r in refused},
                         {("synthia_beta__synthetic_drug_3", "study", "search-default")})
        cut = [r for r in rows if r["pair"] == "cutoff-default"]
        self.assertTrue(all(r["log_loss_difference"] == "0.0" and r["note"] == page.EQUAL_BY_CONSTRUCTION
                            for r in cut))

    def test_a_missing_draw_is_never_read_as_zero_in_the_tables(self):
        tables = self.written(missing=True)
        rows = [r for r in tables["page-datasets.tsv"] if (r["dataset"], r["split"], r["pair"])
                == ("synthia_alpha__synthetic_drug_1", "lineage", "search-default")]
        self.assertEqual(len(rows), 5)
        for r in rows:
            self.assertEqual((r["value"], r["difference"], r["draws_scored"]), ("", "", "2"))
        printed = page.tables(self.out)["datasets"]
        line = [x for x in printed if "synthetic drug 1" in x and "| lineage |" in x][0]
        self.assertIn("none / none / none / none / none", line)
        count = [r for r in tables["page-count.tsv"] if (r["split"], r["pair"], r["learner"], r["measure"])
                 == ("lineage", "search-default", "rf", "balanced_accuracy")][0]
        self.assertEqual((count["datasets"], count["entered"], count["not_scored"]), ("4", "3", "1"))

    def test_one_binary_is_named_once_and_two_with_their_species(self):
        self.build()
        self.write()
        self.assertIn("All 3 species ran one binary", (self.out / "README.md").read_text(encoding="utf-8"))
        self.tmp.cleanup()
        self.tmp = tempfile.TemporaryDirectory()
        self.top, self.out = Path(self.tmp.name), Path(self.tmp.name) / "out"
        self.build()
        (self.reading / "score" / "synthia_gamma" / "binary.txt").write_text(
            f"commit {fx.COMMIT}\ngrm {fx.GRM}\ncode {'d' * 64}\n", encoding="utf-8")
        self.write()
        readme = (self.out / "README.md").read_text(encoding="utf-8")
        self.assertIn("The 3 species ran 2 binaries by their code sha256", readme)
        printed = "\n".join(page.tables(self.out)["binaries"])
        for name in ("*Synthia alpha*", "*Synthia beta*", "*Synthia gamma*"):
            self.assertIn(name, printed)

    def test_the_printed_tables_name_species_in_full_and_recall_beside_log_loss(self):
        self.build()
        self.write()
        printed = page.tables(self.out)
        self.assertIn("*Synthia beta*", "\n".join(printed["species"]))
        self.assertIn("synthetic drug 4", "\n".join(printed["datasets"]))
        head = printed["measures"][0]
        self.assertLess(head.index("Recall"), head.index("Log loss"))
        self.assertIn("0.0126 (0.0088 to 0.0163)", "\n".join(printed["scm"]))
        self.assertIn("no, equal by construction", "\n".join(printed["count"]))
        self.assertIn("no, beside", "\n".join(printed["conditions"]))

    def test_a_species_with_a_missing_draw_counts_its_datasets_beside_those_entered(self):
        rows = self.written(missing=True)["page-species.tsv"]
        alpha = [r for r in rows if (r["species"], r["split"], r["pair"], r["learner"])
                 == ("Synthia alpha", "lineage", "search-default", "rf")]
        self.assertEqual([(r["datasets"], r["entered"], r["not_scored"]) for r in alpha], [("2", "1", "1")])

    def test_a_species_difference_of_minus_0_02_is_counted_the_same(self):
        def change(slug, draws):
            for r in draws:
                if (slug, r["split"], r["arm"]) == ("synthia_gamma", "lineage", "search"):
                    r["balanced_accuracy"] = 0.73
        rows = self.written(change_draws=change)["page-species.tsv"]
        gamma = [r for r in rows if (r["species"], r["split"], r["pair"], r["learner"])
                 == ("Synthia gamma", "lineage", "search-default", "rf")]
        self.assertEqual([(r["entered"], r["better"], r["worse"], r["same"]) for r in gamma], [("1", "0", "0", "1")])

    def test_the_cutoff_arms_log_loss_rows_beside_are_marked_equal_by_construction(self):
        rows = self.written()["page-beside.tsv"]
        marked = [r for r in rows if (r["beside"], r["pair"], r["measure"]) == ("no-constant", "cutoff-default",
                                                                                "log_loss")]
        self.assertEqual(len(marked), 10)
        self.assertTrue(all(r["note"] == page.EQUAL_BY_CONSTRUCTION for r in marked))
        self.assertTrue(all(r["note"] == "" for r in rows if r not in marked))

    def test_a_plain_best_at_trial_0_kept_is_counted_as_the_pick(self):
        rows = {(r["split"], r["learner"]): r for r in self.written()["page-choices.tsv"]}
        self.assertEqual({l: rows[("lineage", l)]["plain_best_is_pick"] for l in fx.LEARNERS},
                         {"scm": "60", "cart": "0", "rf": "60", "gbt": "0", "lasso": "0"})

    def test_an_out_dir_that_is_not_empty_is_refused(self):
        self.build()
        self.out.mkdir()
        (self.out / "x").write_text("x", encoding="utf-8")
        with self.assertRaisesRegex(page.Stop, "is not empty"):
            self.write()


class UnavailableDraw(Case):
    """A default draw that came back unavailable, as score.py marks each row there: split unavailable."""

    DATASET = fx.dataset_name("synthia_beta", "synthetic drug 3")

    def test_an_unavailable_default_draw_is_not_scored_and_leaves_the_count(self):
        self.build(unavailable=True)
        self.write()
        rows = [r for r in tsv(self.out / "page-datasets.tsv") if (r["dataset"], r["split"], r["pair"])
                == (self.DATASET, "lineage", "search-default")]
        self.assertEqual(len(rows), 5)
        for r in rows:
            self.assertEqual((r["value"], r["difference"], r["draws_scored"]), ("", "", "2"))
        for pair in ("search-default", "cutoff-default", "twin-default"):
            count = [r for r in tsv(self.out / "page-count.tsv") if (r["split"], r["pair"], r["learner"], r["measure"])
                     == ("lineage", pair, "rf", "balanced_accuracy")][0]
            self.assertEqual((count["datasets"], count["entered"], count["not_scored"]), ("4", "3", "1"), pair)

    def test_an_unavailable_default_that_owes_choices_files_is_refused(self):
        self.build(unavailable=True)
        folds = self.reading / "score" / "synthia_beta" / "folds.tsv"
        put(folds, [dict(r, folds=str(fx.FOLDS)) for r in tsv(folds)])
        self.refused(r"gate 12 choices: score/synthia_beta: .* lineage seed 43 has no choices file")

    def test_a_column_of_0_in_every_row_is_not_printed(self):
        self.build()
        self.write()
        printed = page.tables(self.out)
        self.assertNotIn("Constant model", printed["count"][0])
        self.assertNotIn("Study run refused", printed["datasets"][0])
        self.assertEqual(len(printed["count"][0].split(" | ")), len(printed["count"][2].split(" | ")))
        self.tmp.cleanup()
        self.tmp = tempfile.TemporaryDirectory()
        self.top, self.out = Path(self.tmp.name), Path(self.tmp.name) / "out"
        self.build(constant=True)
        self.write()
        printed = page.tables(self.out)
        self.assertIn("Constant model", printed["count"][0])
        self.assertIn("Study run refused", printed["datasets"][0])
        self.assertIn("| yes |", "\n".join(printed["datasets"]))


if __name__ == "__main__":
    unittest.main()
