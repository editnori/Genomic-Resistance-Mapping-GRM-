"""scripts/c2-cap/lists_check.py, expect.py and the committed lists (docs/cap-arms.md, Datasets).

The species, datasets and marks here are SYNTHETIC, except that one test
reads the committed scripts/c2-cap/lists.tsv and species.tsv.
"""

import csv
import unittest

from tests.c2_cap import common

check = common.load("lists_check")
expect = common.load("expect")

SP = common.SPECIES
LISTS = [{"dataset": f"{SP}__drug_a", "species_slug": SP, "first_run": "1"},
         {"dataset": f"{SP}__drug_b", "species_slug": SP, "first_run": "1"},
         {"dataset": f"{SP}__drug_c", "species_slug": SP, "first_run": "0"},
         {"dataset": "other__drug_a", "species_slug": "other", "first_run": "1"}]


def marks(m):
    return [{"dataset": r["dataset"], "pairs_S+I": "2" if r["dataset"] in m else "0", "pairs_I+R": "0"} for r in LISTS]


class ListsCheck(unittest.TestCase):
    def test_first_run_datasets_of_the_species_pass(self):
        self.assertEqual(check.strays(LISTS, marks(set()), SP, [f"{SP}__drug_a", f"{SP}__drug_b"]), {})

    def test_a_dataset_of_m_is_not_counted_in_the_first_run(self):
        # drug_b is first run by lists.tsv, but the label count marks it: its labels move with the I rule.
        self.assertEqual(check.strays(LISTS, marks({f"{SP}__drug_b"}), SP, [f"{SP}__drug_a", f"{SP}__drug_b"]),
                         {"M": 1})

    def test_a_dataset_the_marks_do_not_list_is_counted_as_m(self):
        self.assertEqual(check.strays(LISTS, marks(set())[1:], SP, [f"{SP}__drug_a"]), {"M": 1})

    def test_other_strays_are_counted_by_kind(self):
        got = check.strays(LISTS, marks(set()), SP, [f"{SP}__drug_a", f"{SP}__drug_a", f"{SP}__drug_c", "other__drug_a"])
        self.assertEqual(got, {"twice": 1, "not first run": 2, "other species": 1})


class CommittedLists(unittest.TestCase):
    def read(self, name):
        with open(common.KIT / name, encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    def test_the_committed_lists_hold_166_cap_datasets_and_143_in_the_first_run(self):
        rows = self.read("lists.tsv")
        self.assertEqual(len(rows), 166)
        self.assertEqual(sum(r["first_run"] == "1" for r in rows), 143)
        self.assertTrue(all((r["first_run"] == "1") != (r["in_M"] == "1") for r in rows))
        species = self.read("species.tsv")
        self.assertEqual(len(species), 8)
        self.assertTrue(all(s["cohort"] == "2000" for s in species))
        self.assertEqual(sum(int(s["first_run"]) for s in species), 143)


ROOTS = {"fits": {"cap-all": "cap-all"}, "same_as": {"cap-4000": "cap-all"}}


class Expect(unittest.TestCase):
    def rows(self, **kw):
        args = dict(slug=SP, datasets=["d1", "d2"], study=["d2"], draws=3, twin=3, arms=["cap-4000", "cap-all"],
                    splits=["lineage", "study"], roots=ROOTS, have_file=lambda a, d: True)
        args.update(kw)
        return expect.rows(**args)

    def test_each_dataset_lists_draws_the_twin_and_each_arm(self):
        got = self.rows()
        lineage = [r for r in got if r[0] == "d1"]
        self.assertEqual(len(lineage), 3 + 3 + 3 + 3)
        twin = [r for r in lineage if r[5] == "twin"]
        self.assertEqual([(r[2], r[4], r[6]) for r in twin], [("twin", "3", "cap-2000"), ("twin", "4", "cap-2000"),
                                                               ("twin", "5", "cap-2000")])
        self.assertEqual({r[6] for r in lineage if r[2] == "cap-4000"}, {"cap-all"})
        study = [r for r in got if r[3] == "study"]
        self.assertEqual([(r[0], r[2], r[5]) for r in study], [("d2", "cap-2000", "study"),
                                                                ("d2", "cap-4000", "same-cohort"),
                                                                ("d2", "cap-all", "study")])

    def test_an_arm_sharing_another_arms_cohort_is_counted_under_that_arm_only(self):
        got = self.rows()
        self.assertEqual({r[5] for r in got if r[2] == "cap-4000"}, {"same-cohort"})
        self.assertEqual({r[5] for r in got if r[2] == "cap-all"}, {"draw", "study"})
        apart = self.rows(roots={"fits": {"cap-4000": "cap-4000", "cap-all": "cap-all"}, "same_as": {}})
        self.assertEqual({r[5] for r in apart if r[2] == "cap-4000"}, {"draw", "study"})

    def test_an_arm_with_no_fold_file_adds_no_genome_and_is_not_fitted(self):
        got = self.rows(have_file=lambda a, d: False)
        self.assertEqual([r[5] for r in got if r[3] == "study"], ["study", "none-added", "none-added"])

    def test_an_arm_with_no_added_genome_labelled_is_none_added_on_both_splits(self):
        got = self.rows(arms=["cap-all"], roots={"fits": {"cap-all": "cap-all"}, "same_as": {}},
                        empty={("cap-all", "d2")})
        self.assertEqual({(r[3], r[5]) for r in got if r[0] == "d2" and r[2] == "cap-all"},
                         {("lineage", "none-added"), ("study", "none-added")})
        self.assertEqual({r[5] for r in got if r[0] == "d1" and r[2] == "cap-all"}, {"draw"})
        self.assertEqual({r[5] for r in got if r[2] == "cap-2000"}, {"draw", "study"})

    def test_the_none_added_pairs_are_read_from_the_cohort_count(self):
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "added-labelled.tsv"
            p.write_text("dataset\tspecies\tantibiotic\tarm\tlabelled\tadded\tadded_labelled\tadded_resistant\t"
                         "added_susceptible\tnone_added\n"
                         "d1\tSynthia coli\tdrug a\tcap-all\t50\t30\t0\t0\t0\t1\n"
                         "d2\tSynthia coli\tdrug b\tcap-all\t60\t30\t2\t1\t1\t0\n", encoding="utf-8")
            self.assertEqual(expect.none_added(p), {("cap-all", "d1")})

    def test_bad_counts_and_arms_are_refused(self):
        for kw in ({"draws": 0}, {"twin": -1}, {"arms": ["cap-8000"]}, {"arms": ["cap-all", "cap-all"]},
                   {"study": ["d3"]}, {"roots": {"fits": {}, "same_as": {}}}):
            with self.subTest(kw=kw), self.assertRaises(ValueError):
                self.rows(**kw)


if __name__ == "__main__":
    unittest.main()
