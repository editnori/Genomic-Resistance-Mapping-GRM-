"""Item 1 to 3: the population, the label rules and the exclusions."""

import unittest

from tests.outside_study import fixtures as F
from tests.outside_study.fixtures import S, A, B


class Population(unittest.TestCase):
    def test_primary_and_same_study_are_apart(self):
        w = F.world()
        primary = S.population(w["candidates"], F.SPECIES, "primary")
        same = S.population(w["candidates"], F.SPECIES, "same-study")
        self.assertEqual(len(primary), 126)
        self.assertEqual(sorted(same), F.SAME_STUDY)
        self.assertFalse(set(primary) & set(same))
        self.assertNotIn(F.DROPPED, primary)

    def test_only_genomes_on_an_eligible_dataset(self):
        w = F.world()
        w["candidates"][0]["on_eligible_dataset"] = "no"
        self.assertNotIn(w["candidates"][0]["genome_id"], S.population(w["candidates"], F.SPECIES, "primary"))

    def test_eligible_datasets_per_population(self):
        w = F.world()
        self.assertEqual(S.eligible(w["datasets"], F.SPECIES, "primary"), [A, B])
        self.assertEqual(S.eligible(w["datasets"], F.SPECIES, "same-study"), [A])

    def test_a_genome_listed_twice_is_refused(self):
        w = F.world()
        w["candidates"].append(dict(w["candidates"][0]))
        with self.assertRaises(S.Refusal):
            S.population(w["candidates"], F.SPECIES, "primary")


class Labels(unittest.TestCase):
    def setUp(self):
        self.w = F.world()
        self.pop = S.population(self.w["candidates"], F.SPECIES, "primary")
        self.lab = S.labels(self.w["pairs"], self.pop, "primary", [A, B])

    def test_only_r_and_s_are_labels(self):
        for g in ("1001.124", "1001.125", "1001.126"):
            self.assertNotIn((A, g), self.lab, g)
        self.assertIs(self.lab[(A, "1001.121")], True)
        self.assertIs(self.lab[(A, "1001.122")], False)
        self.assertEqual(sum(self.lab[(A, g["id"])] for g in F.new_genomes()), 60)

    def test_counts_per_dataset(self):
        a = [v for (d, _), v in self.lab.items() if d == A]
        b = [v for (d, _), v in self.lab.items() if d == B]
        self.assertEqual((sum(a), len(a) - sum(a)), (62, 61))
        self.assertEqual((sum(b), len(b) - sum(b)), (25, 25))

    def test_a_non_eligible_dataset_gives_no_label(self):
        self.assertFalse([k for k in self.lab if k[0] == F.C])

    def test_same_study_labels_stay_in_their_own_population(self):
        same = S.population(self.w["candidates"], F.SPECIES, "same-study")
        lab = S.labels(self.w["pairs"], same, "same-study", [A])
        self.assertEqual(sorted(g for _, g in lab), F.SAME_STUDY)
        self.assertFalse(set(g for _, g in lab) & set(g for _, g in self.lab))

    def test_intermediate_flagged_as_frozen_is_refused(self):
        p = next(p for p in self.w["pairs"] if p["genome_id"] == "1001.124")
        p["in_frozen_test"] = "yes"
        with self.assertRaises(S.Refusal):
            S.labels(self.w["pairs"], self.pop, "primary", [A, B])

    def test_r_not_flagged_as_frozen_is_refused(self):
        p = next(p for p in self.w["pairs"] if p["genome_id"] == "1001.2" and p["dataset_id"] == A)
        p["in_frozen_test"] = "no"
        with self.assertRaises(S.Refusal):
            S.labels(self.w["pairs"], self.pop, "primary", [A, B])

    def test_frozen_pair_of_a_genome_outside_the_population_is_refused(self):
        p = next(p for p in self.w["pairs"] if p["genome_id"] == F.SAME_STUDY[0])
        p["in_frozen_test"] = "yes"
        with self.assertRaises(S.Refusal):
            S.labels(self.w["pairs"], self.pop, "primary", [A, B])

    def test_a_pair_listed_twice_is_refused(self):
        self.w["pairs"].append(dict(self.w["pairs"][0]))
        with self.assertRaises(S.Refusal):
            S.labels(self.w["pairs"], self.pop, "primary", [A, B])

    def test_scoring_the_same_study_set_alone(self):
        # The same-study set is scored by itself, from its own run. A work
        # directory holding the primary genomes is refused for it.
        with self.assertRaisesRegex(S.Refusal, "outside the frozen list"):
            S.score(F.inputs(F.world(), population="same-study"))


class Exclusions(unittest.TestCase):
    def ex(self, availability, verdicts, counted, ids=("g",)):
        return S.exclusions(list(ids), availability, verdicts, counted)

    def test_in_order_the_first_reason_wins(self):
        a_no = {"g": {"fasta": "no", "detail": "absent"}}
        a_yes = {"g": {"fasta": "yes", "detail": ""}}
        fail = {"g": {"verdict": "fail", "reason": "far"}}
        ok = {"g": {"verdict": "pass", "reason": ""}}
        not_counted = {"g": {"counted": "no", "detail": "x"}}
        counted = {"g": {"counted": "yes", "detail": ""}}
        self.assertEqual(self.ex(a_no, fail, not_counted)[0]["g"][0], "no sequence")
        self.assertEqual(self.ex(a_yes, fail, not_counted)[0]["g"][0], "species check failed")
        self.assertEqual(self.ex(a_yes, ok, not_counted)[0]["g"][0], "not counted")
        self.assertEqual(self.ex(a_yes, ok, counted)[0], {})
        self.assertEqual(S.EXCLUSION_ORDER, ("no sequence", "species check failed", "not counted"))

    def test_borderline_stays_in_and_is_flagged(self):
        out, borderline = self.ex({"g": {"fasta": "yes", "detail": ""}}, {"g": {"verdict": "borderline", "reason": ""}},
                                  {"g": {"counted": "yes", "detail": ""}})
        self.assertEqual((out, borderline), ({}, {"g"}))

    def test_missing_records_are_refused(self):
        a = {"g": {"fasta": "yes", "detail": ""}}
        with self.assertRaises(S.Refusal):
            self.ex({}, {}, {})
        with self.assertRaises(S.Refusal):
            self.ex(a, {}, {})
        with self.assertRaises(S.Refusal):
            self.ex(a, {"g": {"verdict": "pass", "reason": ""}}, {})
        with self.assertRaises(S.Refusal):
            self.ex(a, {"g": {"verdict": "maybe", "reason": ""}}, {"g": {"counted": "yes", "detail": ""}})

    def test_exclusions_never_read_labels(self):
        # Flipping every label leaves the exclusions as they were.
        w = F.world()
        before = S.score(F.inputs(w))["genomes"]
        for p in w["pairs"]:
            p["call"] = {"R": "S", "S": "R"}.get(p["call"], p["call"])
        after = S.score(F.inputs(w))["genomes"]
        self.assertEqual([(g["genome_id"], g["exclusion"]) for g in before],
                         [(g["genome_id"], g["exclusion"]) for g in after])


if __name__ == "__main__":
    unittest.main()
