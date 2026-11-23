"""scripts/c2-cap/cohorts.py and roots.py: the wider cohorts by the release's rule, inside the gate (docs/cap-arms.md).

Every genome, call, dataset and cohort here is SYNTHETIC.
"""

import tempfile
import unittest
from pathlib import Path

from tests.c2_cap import common

cohorts = common.load("cohorts")
roots = common.load("roots")


def calls(n, species=common.SPECIES):
    """n SYNTHETIC genomes with calls on two drugs; genome i has a rarer call the larger i is not."""
    out = []
    for i in range(n):
        g = f"7{i:05d}.1"
        out.append({"genome_id": g, "species": species, "antibiotic": "drug a", "reported": "R" if i % 7 == 0 else "S"})
        if i % 3 == 0:
            out.append({"genome_id": g, "species": species, "antibiotic": "drug b", "reported": "R" if i % 2 else "S"})
    return out


DATASETS = [{"species": common.SPECIES, "antibiotic": "drug a"}, {"species": common.SPECIES, "antibiotic": "drug b"}]


class Rank(unittest.TestCase):
    def test_the_rank_scores_each_call_by_its_class_size_and_breaks_ties_by_id(self):
        ph = [{"genome_id": "1.1", "species": common.SPECIES, "antibiotic": "drug a", "reported": "R"},
              {"genome_id": "2.1", "species": common.SPECIES, "antibiotic": "drug a", "reported": "S"},
              {"genome_id": "3.1", "species": common.SPECIES, "antibiotic": "drug a", "reported": "S"},
              {"genome_id": "3.1", "species": common.SPECIES, "antibiotic": "drug b", "reported": "S"},
              {"genome_id": "4.1", "species": common.SPECIES, "antibiotic": "drug a", "reported": "S", "conflict": True},
              {"genome_id": "5.1", "species": common.SPECIES, "antibiotic": "drug a", "reported": "R",
               "label_from_mic": True},
              {"genome_id": "6.1", "species": "other", "antibiotic": "drug a", "reported": "R"}]
        # 1.1: 1/1; 2.1: 1/2; 3.1: 1/2 + 1/1; the conflict, the MIC call and the other species rank nowhere.
        self.assertEqual(cohorts.rank(ph, DATASETS, common.SPECIES), ["3.1", "1.1", "2.1"])

    def test_the_cohorts_nest_the_default_inside_each_wider_arm(self):
        order = cohorts.rank(calls(5000), DATASETS, common.SPECIES)
        release = order[:2000]
        arms = cohorts.cohorts(order, release, [4000, 4096])
        self.assertEqual(sorted(arms), ["cap-2000", "cap-4000", "cap-all"])
        self.assertEqual((len(arms["cap-4000"]), len(arms["cap-all"])), (4000, 4096))
        self.assertTrue(set(arms["cap-2000"]) <= set(arms["cap-4000"]) <= set(arms["cap-all"]))
        cohorts.check(arms, release, order)


class Gate(unittest.TestCase):
    def setUp(self):
        self.order = cohorts.rank(calls(5000), DATASETS, common.SPECIES)
        self.release = self.order[:2000]

    def test_a_genome_the_rebuilt_table_does_not_hold_is_refused(self):
        arms = cohorts.cohorts(self.order, self.release, [4000, 4096])
        arms["cap-all"] = sorted(arms["cap-all"][:-1] + ["8888888.8"])
        with self.assertRaisesRegex(ValueError, "does not rank"):
            cohorts.check(arms, self.release, self.order)

    def test_a_cohort_past_4096_is_refused(self):
        with self.assertRaisesRegex(ValueError, "past 4,096"):
            cohorts.cohorts(self.order, self.release, [4000, 4097])
        arms = cohorts.cohorts(self.order, self.release, [4000, 4096])
        arms["cap-all"] = self.order[:4097]
        with self.assertRaisesRegex(ValueError, "past 4,096"):
            cohorts.check(arms, self.release, self.order)

    def test_a_wider_cohort_that_drops_a_default_genome_is_refused(self):
        arms = cohorts.cohorts(self.order, self.release, [4000, 4096])
        arms["cap-4000"] = sorted(set(arms["cap-4000"]) - {self.release[5]} | {self.order[4000]})
        with self.assertRaisesRegex(ValueError, "drops 1 genomes of the default"):
            cohorts.check(arms, self.release, self.order)

    def test_a_default_genome_the_rule_does_not_rank_is_refused(self):
        with self.assertRaisesRegex(ValueError, "not ranked"):
            cohorts.cohorts(self.order, self.release[:-1] + ["8888888.8"], [4000])

    def test_a_cap_no_wider_than_the_default_is_refused(self):
        with self.assertRaisesRegex(ValueError, "not wider"):
            cohorts.cohorts(self.order, self.release, [2000])

    def test_a_species_with_fewer_ranked_genomes_gives_equal_wider_cohorts(self):
        order = cohorts.rank(calls(2900), DATASETS, common.SPECIES)
        arms = cohorts.cohorts(order, order[:2000], [4000, 4096])
        self.assertEqual(arms["cap-4000"], arms["cap-all"])
        self.assertEqual(len(arms["cap-all"]), 2900)


class Roots(unittest.TestCase):
    def plan(self, n_mid, n_all, arms):
        with tempfile.TemporaryDirectory() as tmp:
            ids = [f"7{i:05d}.1" for i in range(n_all)]
            common.write_ids(Path(tmp) / "cap-2000.txt", ids[:20])
            common.write_ids(Path(tmp) / "cap-4000.txt", ids[:n_mid])
            common.write_ids(Path(tmp) / "cap-all.txt", ids)
            return roots.plan(tmp, arms)

    def test_equal_wider_cohorts_share_one_root_and_one_fit(self):
        p = self.plan(30, 30, "cap-4000,cap-all")
        self.assertEqual((p["roots"], p["fits"], p["same_as"]),
                         (["cap-all", "cap-2000"], {"cap-all": "cap-all"}, {"cap-4000": "cap-all"}))
        p = self.plan(30, 30, "cap-4000")
        self.assertEqual((p["roots"], p["fits"], p["same_as"]), (["cap-all", "cap-2000"], {"cap-4000": "cap-all"}, {}))

    def test_distinct_wider_cohorts_get_a_root_and_a_fit_each(self):
        p = self.plan(25, 30, "cap-4000,cap-all")
        self.assertEqual((p["roots"], p["fits"], p["same_as"]),
                         (["cap-all", "cap-2000", "cap-4000"], {"cap-4000": "cap-4000", "cap-all": "cap-all"}, {}))
        self.assertEqual(p["genomes"], {"cap-2000": 20, "cap-4000": 25, "cap-all": 30})


if __name__ == "__main__":
    unittest.main()
