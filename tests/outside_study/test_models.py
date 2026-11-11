"""Item 4 and 9: the final models replayed from the reports' published bits,
and the known-gene rule replayed from gene calls."""

import copy
import unittest

from tests.outside_study import fixtures as F
from tests.outside_study.fixtures import S, A, B


def cohort_rows(d, sha, learners=None):
    found = {g: {f: g in F.COHORT_F[f] for f in (1, 2)} for g in F.COHORT}
    return [F.make_row(d, sha, g, lr, found[g]) for g in F.COHORT for lr in (learners or F.LEARNERS[d])]


class FinalModels(unittest.TestCase):
    def setUp(self):
        self.ra, self.rb = F.make_report(A), F.make_report(B)
        self.sa, self.sb = F.sha(F.report_bytes(self.ra)), F.sha(F.report_bytes(self.rb))

    def test_complete_learners(self):
        self.assertEqual(S.complete_learners(self.ra), ["scm", "cart", "lasso"])
        self.assertEqual(S.complete_learners(self.rb), ["scm", "cart"])

    def test_bits_are_read_little_endian_per_sample(self):
        b = F.bits({"999.1", "999.5"})
        self.assertEqual([S.sample_bit(b, i) for i in range(8)],
                         [True, False, False, False, True, False, False, False])

    def test_final_calls(self):
        for found, want in (({1: True, 2: False}, ("resistant", None)), ({1: False, 2: True}, ("susceptible", None))):
            self.assertEqual(S.final_call(self.ra, "scm", found), want)
        self.assertEqual(S.final_call(self.ra, "cart", {1: True}), ("resistant", 0.9))
        self.assertEqual(S.final_call(self.ra, "cart", {1: False}), ("susceptible", 0.1))
        self.assertEqual(S.final_call(self.rb, "cart", {}), ("susceptible", 0.3))
        call, p = S.final_call(self.ra, "lasso", {1: False, 2: True})
        self.assertEqual(call, "susceptible")
        self.assertAlmostEqual(p, 1 / (1 + 2.718281828459045 ** 0.7))
        r = copy.deepcopy(self.ra)
        r["models"]["scm"] = {"kind": "disjunction", "rules": [{"feature": 1, "presence": False},
                                                               {"feature": 2, "presence": True}]}
        self.assertEqual(S.final_call(r, "scm", {1: True, 2: False})[0], "susceptible")
        self.assertEqual(S.final_call(r, "scm", {1: True, 2: True})[0], "resistant")
        r["models"]["lasso"]["coef"] = []
        self.assertEqual(S.final_call(r, "lasso", {}), (None, None))

    def test_cohort_replay_passes_on_the_published_bits(self):
        problems, n = S.replay_problems(self.ra, self.sa, cohort_rows(A, self.sa))
        self.assertEqual((problems, n), ([], 8))
        problems, n = S.replay_problems(self.rb, self.sb, cohort_rows(B, self.sb))
        self.assertEqual((problems, n), ([], 8))

    def test_cohort_replay_finds_a_wrong_bit(self):
        rows = cohort_rows(A, self.sa)
        r = next(r for r in rows if r["genome"] == "999.5" and r["learner"] == "scm")
        r["features"][0]["found"] = False
        r["call"] = "susceptible"
        problems, _ = S.replay_problems(self.ra, self.sa, rows)
        self.assertTrue(any("published bit" in p for p in problems), problems)

    def test_cohort_replay_finds_a_missing_learner_and_a_stranger(self):
        rows = [r for r in cohort_rows(A, self.sa) if not (r["genome"] == "999.3" and r["learner"] == "lasso")]
        stranger = copy.deepcopy(rows[0])
        stranger["genome"] = "1001.1"
        problems, _ = S.replay_problems(self.ra, self.sa, rows + [stranger])
        self.assertTrue(any("999.3: rows for" in p for p in problems), problems)
        self.assertTrue(any("not a sample" in p for p in problems), problems)


class KnownGene(unittest.TestCase):
    def setUp(self):
        self.report = F.make_report(A)

    def test_token_match(self):
        self.assertTrue(S.token_match("C:BETA-LACTAM", "BETA-LACTAM", "PENICILLIN"))
        self.assertTrue(S.token_match("S:METHICILLIN", "BETA-LACTAM", "CEPHALOSPORIN/METHICILLIN"))
        self.assertTrue(S.token_match("C:QUINOLONE", "AMINOGLYCOSIDE/ QUINOLONE", "X"))
        self.assertFalse(S.token_match("S:BETA-LACTAM", "BETA-LACTAM", "PENICILLIN"))
        self.assertFalse(S.token_match("C:BETA", "BETA-LACTAM", "PENICILLIN"))

    def test_explain_needs_every_rule_group(self):
        rule = [["C:BETA-LACTAM"], ["S:METHICILLIN", "S:OXACILLIN"]]
        self.assertEqual(S.explain(rule, [F.BLAZ, F.MECA]), (True, ["blaZ", "mecA"]))
        self.assertEqual(S.explain(rule, [F.BLAZ]), (False, ["blaZ"]))
        self.assertEqual(S.explain(rule, []), (False, []))

    def test_replay_on_the_cohort_passes(self):
        self.assertEqual(S.known_gene_replay(self.report, F.COHORT_HITS), [])

    def test_replay_fails_on_a_changed_call(self):
        hits = dict(F.COHORT_HITS)
        hits["999.7"] = [F.BLAZ]
        problems = S.known_gene_replay(self.report, hits)
        self.assertTrue(any("999.7: groups" in p for p in problems), problems)
        self.assertTrue(any("999.7: catalog_resistant True" in p for p in problems), problems)

    def test_replay_fails_on_the_table_and_counts(self):
        r = copy.deepcopy(self.report)
        r["catalog"]["groups"][1]["counts_for_rule"] = False
        self.assertIn("the catalog.groups table does not replay", S.known_gene_replay(r, F.COHORT_HITS))
        r = copy.deepcopy(self.report)
        r["catalog"]["metrics"]["tp"] = 3
        self.assertIn("catalog tp 2, published 3", S.known_gene_replay(r, F.COHORT_HITS))

    def test_replay_needs_every_cohort_genome(self):
        hits = dict(F.COHORT_HITS)
        del hits["999.8"]
        self.assertIn("999.8: no gene calls counted", S.known_gene_replay(self.report, hits))
        self.assertEqual(S.known_gene_replay(self.report, hits, every_sample=False), [])

    def test_no_rule_published(self):
        self.assertEqual(S.known_gene_replay(F.make_report(B), F.COHORT_HITS),
                         ["the report publishes no catalog rule"])


if __name__ == "__main__":
    unittest.main()
