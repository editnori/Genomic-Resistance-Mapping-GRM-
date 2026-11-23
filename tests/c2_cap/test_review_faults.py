"""SYNTHETIC tests for faults the review planted in scripts/c2-cap/page_tables.py and run_checks.py that no test of
the lane caught.

Each test here fails when the one rule it names is broken: a gate or check that lets a wrong copy through, or a table
rule that reads a figure another way. The reading is tests/c2_cap/page_common.py's made-up one, or a made-up row
built in the test. Every figure, name and count here is made up; none comes from a run.
"""

import json
import unittest

from tests.c2_cap import common, page_common as pc, test_page_tables as tpt, test_review_checks as trc

page = common.load("page_tables")
checks = trc.checks

# Made-up summed calls 0.01 above the default's balanced accuracy: positive, and not better by the 0.02 margin.
NUDGED = pc.C9.measures({"n": 100, "tp": 40, "fn": 10, "tn": 41, "fp": 9, "loss_sum": 39.5})


class ReviewGates(tpt.Base):
    def test_review_gate_4_refuses_a_wider_cohort_short_of_the_defaults_genomes(self):
        p = self.r.runs / self.r.slugs[5] / "results" / "cohorts.json"
        doc = json.loads(p.read_text(encoding="utf-8"))
        doc["arms"]["cap-all"]["holds_default"] = 1999
        p.write_text(json.dumps(doc), encoding="utf-8")
        self.stops_at(4)

    def test_review_gate_12_refuses_a_failed_check_of_the_scorer(self):
        slug = self.r.slugs[6]
        rows = self.r.check_rows(slug)
        rows[0]["passed"] = 0
        pc.tsv(self.r.runs / slug / "score" / "checks.tsv", rows, ["row", "check", "passed", "detail"])
        self.stops_at(12)

    def test_review_gate_3_refuses_a_dataset_named_otherwise_than_in_lists_tsv(self):
        slug = self.r.slugs[7]
        d = f"{slug}__drug_b"
        self.r.edit(self.draws(slug), lambda rows: [dict(r, antibiotic="drug z") if r["dataset"] == d else r
                                                    for r in rows])
        self.stops_at(3)

    def test_review_gate_11_refuses_a_dataset_row_that_is_not_the_mean_of_its_draws(self):
        slug = self.r.slugs[0]
        p = self.r.runs / slug / "score" / "datasets.tsv"

        def change(rows):
            for r in rows:
                if (r["dataset"], r["split"], r["arm"], r["learner"]) == (f"{slug}__drug_b", "lineage", "cap-all", "scm"):
                    r["balanced_accuracy"] = str(float(r["balanced_accuracy"]) + 0.01)
            return rows
        self.r.edit(p, change)
        self.r.recount()
        self.stops_at(11)


class ReviewReading(tpt.Base):
    def test_review_the_reading_is_the_deciding_counts_where_the_few_added_dataset_swings_it(self):
        # cap-all 0.01 above the default on 5 of its 14 lineage datasets: same by the margin, so 9 are better, the
        # few-added one among them, and the sign test's p after Holm sits either side of 0.05 with it and without it.
        for slug, drugs in ((self.r.slugs[4], ("drug_a", "drug_b")), (self.r.slugs[5], ("drug_a", "drug_b")),
                            (self.r.slugs[6], ("drug_a",))):
            ds = {f"{slug}__{d}" for d in drugs}
            self.rewrite_draws(slug, lambda rows, ds=ds: [
                dict(r, **{k: pc.C9.fmt(v) for k, v in NUDGED.items()})
                if r["dataset"] in ds and r["arm"] == "cap-all" and r["split"] == "lineage" else r for r in rows])
        self.write()
        runs = {s: page.read_species(self.r.runs / s) for s in self.r.slugs}
        lrows = page.load(self.r.lists)
        names = {r["dataset"]: (r["species"], r["antibiotic"]) for r in page.first_run(lrows)}
        few = page.few_added(runs, names)
        self.assertEqual([r["dataset"] for r in few], [self.r.few])
        beside = {r["condition"]: str(r["met"])
                  for r in page.page_reading(page.counted(page.without_few(page.all_datasets(runs), few))[1])}
        self.assertEqual(beside["1 and 2"], "0")
        got = {r["condition"]: r["met"] for r in self.rows("page-reading.tsv")}
        self.assertEqual(got["1 and 2"], "1")
        want = {r["condition"]: str(r["met"]) for r in page.page_reading(page.load(self.r.count / "conditions.tsv"))}
        self.assertEqual(got, want)


class ReviewTables(unittest.TestCase):
    def test_review_a_difference_of_exactly_the_margin_is_neither_better_nor_worse(self):
        self.assertEqual(page.tally([0.02, -0.02, 0.020000000001, -0.020000000001, 0.0]), (1, 1, 3))
        self.assertEqual(page.tally([0.02, -0.02, -0.03], "log_loss"), (1, 0, 2))

    def test_review_the_few_added_rule_takes_fewer_than_20_labelled_genomes(self):
        names = {f"x__{d}": ("Synthetic species", d) for d in ("drug_19", "drug_20", "drug_0")}
        added = [{"dataset": "x__drug_19", "arm": "cap-all", "none_added": "0", "added_labelled": "19",
                  "added_resistant": "9", "added_susceptible": "10"},
                 {"dataset": "x__drug_20", "arm": "cap-all", "none_added": "0", "added_labelled": "20",
                  "added_resistant": "10", "added_susceptible": "10"},
                 {"dataset": "x__drug_0", "arm": "cap-all", "none_added": "1", "added_labelled": "0",
                  "added_resistant": "0", "added_susceptible": "0"},
                 {"dataset": "x__drug_20", "arm": "cap-4000", "none_added": "0", "added_labelled": "5",
                  "added_resistant": "2", "added_susceptible": "3"}]
        got = page.few_added({"x": {"added": added}}, names)
        self.assertEqual([(r["dataset"], r["added_labelled"]) for r in got], [("x__drug_19", 19)])

    def test_review_the_cohort_as_fitted_takes_out_the_added_genomes_that_failed_the_fetch(self):
        arms = {"cap-2000": {"genomes": 2000}, "cap-4000": {"genomes": 4000}, "cap-all": {"genomes": 4096}}
        run = {"cohorts": {"arms": arms}, "fetch": {"left_wider_cohorts": {"cap-4000": 7, "cap-all": 7}}}
        self.assertEqual(page.cohort_sizes(run), {"cap-2000": 2000, "cap-4000": 3993, "cap-all": 4089})
        taken = {"cohorts": {"arms": {"cap-2000": {"genomes": 2000}, "cap-all": {"genomes": 4089, "fetch_failed": 7}}},
                 "fetch": {"left_wider_cohorts": {"cap-all": 7}}}
        self.assertEqual(page.cohort_sizes(taken), {"cap-2000": 2000, "cap-all": 4089})

    def test_review_a_study_dataset_with_no_added_genome_that_trains_is_named_not_fitted(self):
        d = "x__drug_s"
        expect = [{"dataset": d, "species_slug": "x", "arm": "cap-2000", "split": "study", "slot": "-",
                   "role": "study", "fit": "cap-2000"},
                  {"dataset": d, "species_slug": "x", "arm": "cap-all", "split": "study", "slot": "-",
                   "role": "none-added", "fit": "cap-all"}]
        run = {"expect": expect, "added": [{"dataset": d, "arm": "cap-all", "none_added": "0"}], "draws.tsv": [],
               "refused": {}}
        got = page.page_not_fitted({"x": run}, [{"dataset": d, "species": "Synthetic species", "antibiotic": "drug s"}])
        self.assertEqual([(r["dataset"], r["species"], r["antibiotic"], r["arm"], r["split"], r["reason"]) for r in got],
                         [(d, "Synthetic species", "drug s", "cap-all", "study",
                           "none-added: no added genome trains on the study split")])


class ReviewRunChecks(trc.Copies):
    def done(self, slug, named=None, sums=None):
        f = self.folder(slug)
        b = checks.binary_of(f)
        (f / "DONE").write_text(f"DONE {named or slug} 2026-11-23T14:00:00Z binary {b['grm']} code {b['code']} fits 4 "
                                f"SHA256SUMS {sums or checks.sha256(f / 'SHA256SUMS')}\n", encoding="utf-8")

    def test_review_a_done_that_names_another_species_fails(self):
        self.done(self.r.slugs[2], named=self.r.slugs[3])
        got = self.fails("DONE")
        self.assertIn("7 of 8 name their species", got["DONE"])

    def test_review_a_done_that_names_another_sha256sums_fails(self):
        self.done(self.r.slugs[4], sums="0" * 64)
        got = self.fails("DONE")
        self.assertIn("7 of 8 name their species", got["DONE"])

    def test_review_a_wider_cohort_short_of_the_defaults_genomes_or_its_lines_fails(self):
        slug = self.r.slugs[6]
        p = self.folder(slug) / "results" / "cohorts.json"
        clean = p.read_text(encoding="utf-8")
        for key, value in (("holds_default", 1999), ("lines", 4095)):
            with self.subTest(key=key):
                doc = json.loads(clean)
                doc["arms"]["cap-all"][key] = value
                p.write_text(json.dumps(doc), encoding="utf-8")
                self.seal(slug)
                got = self.fails("cohorts")
                self.assertIn("cap-all does not hold the default's 2,000 or its list's lines differ", got["cohorts"])
        p.write_text(clean, encoding="utf-8")
        self.seal(slug)
        got, ok = self.run_checks()
        self.assertTrue(ok, got["cohorts"])

    def test_review_a_cohorts_json_whose_ranked_count_is_not_the_rules_fails(self):
        slug = self.r.slugs[1]
        self.selection[slug] = (4999, 2000, 2096)
        got = self.fails("cohorts")
        self.assertIn("ranked 5000, not 4999", got["cohorts"])


if __name__ == "__main__":
    unittest.main()
