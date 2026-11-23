"""SYNTHETIC tests of the result page's table script, scripts/c2-cap/page_tables.py (docs/cap-arms-result.md).

Every reading here is made up by tests/c2_cap/page_common.py: the 8 species of the run with made-up datasets and
figures. Each gate has a test whose reading breaks it, and scripts/c2-cap/page_faults.py plants a fault in each
gate and table rule for the test named for it to catch.
"""

import json
import tempfile
import unittest
from pathlib import Path

from tests.c2_cap import common, page_common as pc

page = common.load("page_tables")


class Base(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.r = pc.Reading(self.tmp).write()
        self.out = self.tmp / "out"

    def tearDown(self):
        self._tmp.cleanup()

    def write(self):
        return page.write(self.r.count, self.r.runs, self.r.lists, self.out)

    def stops_at(self, gate):
        with self.assertRaises(page.Stop) as cm:
            self.write()
        self.assertTrue(str(cm.exception).startswith(f"gate {gate} "), str(cm.exception))
        self.assertFalse((self.out / "SHA256SUMS").exists())

    def rows(self, name):
        return pc.read(self.out / name)

    def draws(self, slug):
        return self.r.runs / slug / "score" / "draws.tsv"

    def rewrite_draws(self, slug, fn):
        """draws.tsv through fn, datasets.tsv as the scorer makes it from them, and the count again."""
        rows = fn(pc.read(self.draws(slug)))
        typed = [page.typed_draw(r) for r in rows]
        pc.C9.write_tsv(self.draws(slug), typed, pc.score.DRAW_COLUMNS)
        pc.C9.write_tsv(self.r.runs / slug / "score" / "datasets.tsv", pc.C3.dataset_rows(typed), pc.C3.DATASET_COLUMNS)
        self.r.recount()


class Gates(Base):
    def test_a_clean_reading_passes_every_gate_and_writes_its_tables(self):
        files = self.write()
        for name in ("count/count.tsv", "count/conditions.tsv", "page-reading.tsv", "page-conditions.tsv",
                     "page-count.tsv", "page-count-without-few.tsv", "page-few-added.tsv", "page-species.tsv",
                     "page-datasets.tsv", "page-not-fitted.tsv", "page-study.tsv", "page-checks.tsv",
                     "page-binaries.tsv", "page-times.tsv", "gates.txt", "README.md"):
            self.assertIn(name, files)
        self.assertEqual((self.out / "count" / "count.tsv").read_bytes(), (self.r.count / "count.tsv").read_bytes())
        self.assertEqual(len((self.out / "gates.txt").read_text().splitlines()), len(page.GATES))
        listed = [line.split("  ", 1)[1] for line in (self.out / "SHA256SUMS").read_text().splitlines()]
        self.assertEqual(listed, files)

    def test_a_species_folder_missing_from_the_runs_is_refused(self):
        import shutil
        shutil.rmtree(self.r.runs / self.r.slugs[5])
        self.stops_at(2)

    def test_a_dataset_of_m_in_a_species_scored_tables_is_refused(self):
        slug = self.r.slugs[0]
        m = f"{slug}__drug_m"

        def add_m(rows):
            extra = [dict(r, dataset=m, antibiotic="drug m") for r in rows if r["dataset"] == f"{slug}__drug_b"]
            return rows + extra
        self.rewrite_draws(slug, add_m)
        self.stops_at(3)

    def test_a_cohort_json_past_4096_genomes_is_refused(self):
        slug = self.r.slugs[4]
        p = self.r.runs / slug / "results" / "cohorts.json"
        doc = json.loads(p.read_text())
        doc["arms"]["cap-all"]["genomes"] = 4097
        p.write_text(json.dumps(doc))
        self.stops_at(4)

    def test_a_binary_txt_with_no_code_sha256_is_refused(self):
        (self.r.runs / self.r.slugs[2] / "binary.txt").write_text(f"commit {pc.COMMIT}\ngrm {pc.GRM}\n")
        self.stops_at(5)

    def test_a_none_added_dataset_scored_under_cap_all_is_refused(self):
        slug, d = next(iter(self.r.none_added))

        def add(rows):
            return rows + [dict(r, arm="cap-all", fit="cap-all") for r in rows
                           if r["dataset"] == d and r["arm"] == "cap-2000" and r["split"] == "lineage"]
        self.rewrite_draws(slug, add)
        self.stops_at(6)

    def test_a_withdrawn_arms_condition_read_as_an_arm_that_ran_is_refused(self):
        p = self.r.count / "conditions.tsv"
        self.r.edit(p, lambda rows: [dict(r, withdrawn="0", met="0", learners="0") if r["arm"] == "cap-4000" else r
                                     for r in rows])
        self.stops_at(7)

    def test_a_twin_row_from_a_wider_arms_fit_is_refused(self):
        slug = self.r.slugs[6]
        self.rewrite_draws(slug, lambda rows: [dict(r, fit="cap-all") if r["arm"] == "twin" else r for r in rows])
        self.stops_at(8)

    def test_a_study_draw_put_on_the_lineage_split_is_refused(self):
        slug = self.r.slugs[0]
        self.r.edit(self.draws(slug), lambda rows: [dict(r, split="lineage") if r["split"] == "study" else r
                                                    for r in rows])
        self.stops_at(9)

    def test_an_unavailable_lineage_split_scored_on_one_arm_is_refused(self):
        slug, d = self.r.slugs[3], self.r.refused

        def score_cap_all(rows):
            out = []
            for r in rows:
                if r["dataset"] == d and r["split"] == "lineage" and r["arm"] == "cap-all":
                    r = dict(r, scored="1", why="", **{k: pc.C9.fmt(v) for k, v in pc.calls("cap-all", r["learner"]).items()})
                out.append(r)
            return out
        self.rewrite_draws(slug, score_cap_all)
        self.stops_at(10)

    def test_a_draw_missing_from_the_scorers_tables_is_refused(self):
        slug = self.r.slugs[7]
        d = f"{slug}__drug_a"
        self.rewrite_draws(slug, lambda rows: [r for r in rows if not (r["dataset"] == d and r["arm"] == "cap-all"
                                                                       and r["slot"] == "2" and r["learner"] == "rf")])
        self.stops_at(11)

    def test_a_cap_all_row_with_no_check_of_its_scored_genomes_is_refused(self):
        slug = self.r.slugs[5]
        rows = [c for c in self.r.check_rows(slug) if not c["row"].endswith("__drug_b lineage slot 1")]
        pc.tsv(self.r.runs / slug / "score" / "checks.tsv", rows, ["row", "check", "passed", "detail"])
        self.stops_at(12)

    def test_a_count_taken_without_the_few_added_datasets_is_refused(self):
        rows = [r for r in self.r.all_rows() if not (r["arm"] == "cap-all" and r["dataset"] == self.r.few)]
        self.r.recount(rows)
        self.stops_at(13)

    def test_a_count_cell_that_is_not_the_scorers_is_refused(self):
        def change(rows):
            for r in rows:
                if r["pair"] == "cap-all-cap-2000" and r["measure"] == "recall" and r["learner"] == "gbt":
                    r["better"] = str(int(r["better"]) + 1)
            return rows
        self.r.edit(self.r.count / "count.tsv", change)
        self.stops_at(14)

    def test_an_out_dir_that_is_not_empty_is_refused(self):
        self.out.mkdir()
        (self.out / "x").write_text("x")
        with self.assertRaises(page.Stop):
            self.write()


class Tables(Base):
    def test_species_and_datasets_are_named_in_full(self):
        self.write()
        species = self.rows("page-species.tsv")
        self.assertEqual([r["species"] for r in species], [self.r.name[s] for s in sorted(self.r.slugs)])
        for r in self.rows("page-datasets.tsv"):
            self.assertEqual(r["species"], self.r.name[r["dataset"].split("__")[0]])
            self.assertEqual(r["antibiotic"], r["dataset"].split("__")[1].replace("_", " "))

    def test_the_study_split_decides_nothing_and_condition_5_is_not_read(self):
        self.write()
        conds = self.rows("page-conditions.tsv")
        five = [r for r in conds if r["arm"] == "cap-all" and r["reading"] == "ruling 15" and r["condition"] == "5"]
        self.assertEqual(five[0]["decides"], "decides nothing: it reads the study split")
        self.assertEqual([r["condition"] for r in self.rows("page-reading.tsv")], ["1 and 2", "3", "4", "all deciding"])
        study = [r for r in self.rows("page-count.tsv") if r["split"] == "study" and r["withdrawn"] == "0"]
        self.assertTrue(study)
        for r in study:
            self.assertEqual(r["reading"], "the study split: decides nothing")
        for r in self.rows("page-study.tsv"):
            self.assertEqual(r["reading"], "the study split: decides nothing")

    def test_cap_4000_is_named_withdrawn_and_holds_no_figure(self):
        self.write()
        for r in self.rows("page-conditions.tsv"):
            if r["arm"] == "cap-4000":
                self.assertEqual((r["met"], r["decides"]), ("", "withdrawn before any fit: not read"))
        mid = [r for r in self.rows("page-count.tsv") if r["pair"] == "cap-4000-cap-2000"]
        self.assertTrue(mid)
        for r in mid:
            self.assertEqual((r["reading"], r["entered"], r["meets_1"], r["median_difference"]),
                             ("withdrawn before any fit", "0", "", ""))

    def test_the_few_added_datasets_stay_in_the_count_and_the_count_beside_decides_nothing(self):
        self.write()
        few = self.rows("page-few-added.tsv")
        self.assertEqual([(r["species"], r["antibiotic"], r["added_labelled"], r["added_resistant"],
                           r["added_susceptible"]) for r in few], [(self.r.name[self.r.slugs[2]], "drug a", "5", "2", "3")])
        count = {(r["split"], r["pair"], r["learner"], r["measure"]): r for r in self.rows("page-count.tsv")}
        beside = {(r["split"], r["pair"], r["learner"], r["measure"]): r for r in self.rows("page-count-without-few.tsv")}
        key = ("lineage", "cap-all-cap-2000", "lasso", "balanced_accuracy")
        self.assertEqual(count[key]["reading"], "the count")
        self.assertEqual(int(count[key]["datasets"]), int(beside[key]["datasets"]) + 1)
        scorer = {(r["split"], r["pair"], r["learner"], r["measure"]): r for r in pc.read(self.r.count / "count.tsv")}
        self.assertEqual(count[key]["datasets"], scorer[key]["datasets"])
        for r in beside.values():
            self.assertIn(r["reading"], ("beside: decides nothing", "withdrawn before any fit"))

    def test_an_unavailable_lineage_split_is_named_with_its_reason_and_holds_no_figure(self):
        self.write()
        gone = [r for r in self.rows("page-not-fitted.tsv") if r["dataset"] == self.r.refused]
        self.assertEqual(len(gone), 1)
        self.assertEqual((gone[0]["arm"], gone[0]["split"], gone[0]["refused_slots"]),
                         ("every arm and the twin", "lineage", "3"))
        self.assertTrue(gone[0]["reason"].startswith("lineage split unavailable"))
        rows = [r for r in self.rows("page-datasets.tsv") if r["dataset"] == self.r.refused and r["split"] == "lineage"]
        self.assertEqual(len(rows), 5)
        for r in rows:
            self.assertEqual((r["role"], r["default"], r["cap_all"], r["difference"], r["twin"]),
                             ("lineage split unavailable", "", "", "", ""))

    def test_a_none_added_dataset_is_named_and_holds_no_cap_all_figure(self):
        self.write()
        slug, d = next(iter(self.r.none_added))
        named = [r for r in self.rows("page-not-fitted.tsv") if r["dataset"] == d]
        self.assertEqual([(r["arm"], r["reason"]) for r in named], [("cap-all", "none-added: no added genome labelled")])
        for r in self.rows("page-datasets.tsv"):
            if r["dataset"] == d:
                self.assertEqual((r["role"], r["cap_all"], r["difference"]), ("none-added", "", ""))

    def test_species_on_two_binaries_are_each_named_with_their_binary(self):
        self.r.binary[self.r.slugs[7]] = (pc.COMMIT, pc.GRM, "e" * 64)
        self.r.write_species(self.r.slugs[7])
        self.write()
        rows = self.rows("page-binaries.tsv")
        self.assertEqual([(r["code"], r["species"]) for r in rows], [(pc.CODE, "7"), ("e" * 64, "1")])
        self.assertEqual(rows[1]["names"], self.r.name[self.r.slugs[7]])

    def test_the_per_dataset_figures_are_the_scorers_means_and_differences(self):
        self.write()
        row = [r for r in self.rows("page-datasets.tsv")
               if r["dataset"] == f"{self.r.slugs[0]}__drug_a" and r["split"] == "lineage" and r["learner"] == "scm"][0]
        self.assertEqual((row["default"], row["cap_all"], row["difference"], row["twin"], row["twin_difference"]),
                         ("0.8", "0.84", "0.04", "0.8", "0.0"))
        self.assertEqual((row["draws"], row["default_draws_scored"], row["cap_all_draws_scored"]), ("3", "3", "3"))
        self.assertEqual((row["added_labelled"], row["added_resistant"], row["added_susceptible"]), ("400", "200", "200"))

    def test_the_species_count_takes_lineage_pairs_only(self):
        self.write()
        sp = {r["species_slug"]: r for r in self.rows("page-species.tsv")}
        # 2 datasets by 5 learners; the first species' drug_a also has a study row, which stays out.
        self.assertEqual((sp[self.r.slugs[0]]["cap_all_pairs"], sp[self.r.slugs[0]]["cap_all_better"]), ("10", "10"))
        self.assertEqual((sp[self.r.slugs[1]]["cap_all_pairs"], sp[self.r.slugs[1]]["none_added"]), ("5", "1"))
        self.assertEqual((sp[self.r.slugs[3]]["cap_all_pairs"], sp[self.r.slugs[3]]["cap_all_not_scored"]), ("5", "5"))

    def test_the_checks_and_times_are_summed_per_name_and_arm(self):
        self.write()
        checks = {r["check"]: r for r in self.rows("page-checks.tsv")}
        same = checks["the scored genomes of each outer fold are the default's"]
        # 8 species by 2 datasets by 3 slots, less the none-added dataset's 3 and the unavailable one's 3.
        self.assertEqual((same["passed"], same["made"], same["species"]), ("42", "42", "8"))
        times = {(r["arm"], r["split"]): r for r in self.rows("page-times.tsv")}
        self.assertEqual((times[("cap-all", "lineage")]["fits"], times[("cap-all", "lineage")]["wall_s_sum"]),
                         ("15", "300.0"))
        self.assertEqual(times[("cap-2000", "lineage")]["wall_s_median"], "10.0")

    def test_the_print_holds_each_table(self):
        self.write()
        printed = page.tables(self.out)
        self.assertEqual(set(printed), {"reading", "conditions", "family", "count", "without few", "few added",
                                        "species", "not fitted", "datasets", "checks", "binaries", "times"})
        self.assertIn(f"| *{self.r.name[self.r.slugs[2]]}* | drug a | 5 | 2 | 3 |", printed["few added"])


if __name__ == "__main__":
    unittest.main()
