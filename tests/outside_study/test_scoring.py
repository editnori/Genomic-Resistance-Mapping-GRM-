"""The tables and the reading on the hand-built world (fixtures.py), each
number worked out by hand from the fixture's design."""

import unittest

from tests.outside_study import fixtures as F
from tests.outside_study.fixtures import S, A, B


def by(rows, **want):
    out = [r for r in rows if all(r[k] == v for k, v in want.items())]
    if len(out) != 1:
        raise AssertionError(f"{len(out)} rows match {want}")
    return out[0]


class Scoring(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.res = S.score(F.inputs(F.world()))

    def metric(self, d, learner):
        return by(self.res["metrics"], dataset=d, learner=learner)

    def test_confusion_counts_per_dataset_and_learner(self):
        for learner in ("scm", "cart", "lasso"):
            r = self.metric(A, learner)
            self.assertEqual((r["tp"], r["fn"], r["fp"], r["tn"]), (48, 12, 7, 53), learner)
        r = self.metric(A, "known_gene")
        self.assertEqual((r["tp"], r["fn"], r["fp"], r["tn"]), (50, 10, 0, 60))
        r = self.metric(B, "scm")
        self.assertEqual((r["tp"], r["fn"], r["fp"], r["tn"]), (22, 3, 3, 22))
        r = self.metric(B, "cart")
        self.assertEqual((r["tp"], r["fn"], r["fp"], r["tn"]), (0, 25, 0, 25))

    def test_learners_scored(self):
        # drugB's lasso is incomplete and it has no catalog rule; forests and
        # boosted trees are never scored.
        got = sorted((r["dataset"], r["learner"]) for r in self.res["metrics"])
        self.assertEqual(got, sorted([(A, "scm"), (A, "cart"), (A, "lasso"), (A, "known_gene"),
                                      (B, "scm"), (B, "cart")]))

    def test_estimates(self):
        r = self.metric(A, "scm")
        self.assertAlmostEqual(r["sensitivity"], 48 / 60)
        self.assertAlmostEqual(r["specificity"], 53 / 60)
        self.assertAlmostEqual(r["balanced_accuracy"], (48 / 60 + 53 / 60) / 2)
        self.assertAlmostEqual(r["precision"], 48 / 55)
        self.assertAlmostEqual(r["vme"], 12 / 60)
        self.assertAlmostEqual(r["me"], 7 / 60)

    def test_intervals_are_clopper_pearson_on_the_right_counts(self):
        r = self.metric(A, "scm")
        lo, hi = S.clopper_pearson(48, 60)
        self.assertEqual((r["sensitivity_lo"], r["sensitivity_hi"]), (lo, hi))
        lo, hi = S.clopper_pearson(53, 60)
        self.assertEqual((r["specificity_lo"], r["specificity_hi"]), (lo, hi))
        for m in ("sensitivity", "specificity"):
            self.assertLess(r[f"{m}_lo"], r[m])
            self.assertLess(r[m], r[f"{m}_hi"])
        # 48 of 60 to six places (scipy.stats.beta.ppf gives 0.676700, 0.892159).
        self.assertEqual(S.fmt(r["sensitivity_lo"]), "0.676700")
        self.assertEqual(S.fmt(r["sensitivity_hi"]), "0.892159")
        k = self.metric(A, "known_gene")
        self.assertEqual(k["specificity_hi"], 1.0)
        self.assertAlmostEqual(k["specificity_lo"], 0.025 ** (1 / 60), delta=1e-12)

    def test_references_carried_from_both_splits(self):
        r = self.metric(A, "scm")
        self.assertEqual((r["lineage_sensitivity"], r["lineage_specificity"]), (0.95, 0.90))
        self.assertEqual((r["random_sensitivity"], r["random_specificity"]), (0.85, 0.80))
        self.assertEqual((r["lineage_n"], r["lineage_tp"], r["lineage_fn"], r["lineage_fp"], r["lineage_tn"]),
                         (200, 95, 5, 10, 90))
        self.assertEqual((r["random_tp"], r["random_tn"]), (85, 80))
        k = self.metric(A, "known_gene")
        self.assertEqual(k["lineage_sensitivity"], 0.50, "the known-gene rule is set beside the report's catalog row")
        b = self.metric(B, "scm")
        self.assertIsNone(b["lineage_sensitivity"])
        self.assertIn("unavailable", b["lineage_note"])
        self.assertEqual(b["random_sensitivity"], 0.70)

    def test_readings(self):
        want = {
            ("scm", "sensitivity"): "below", ("scm", "specificity"): "consistent",
            ("cart", "sensitivity"): "consistent", ("cart", "specificity"): "above",
            ("lasso", "sensitivity"): "above", ("lasso", "specificity"): "below",
            ("known_gene", "sensitivity"): "above", ("known_gene", "specificity"): "consistent",
        }
        for (learner, m), v in want.items():
            self.assertEqual(self.metric(A, learner)[f"reading_{m}"], v, (learner, m))
        for learner in ("scm", "cart"):
            for m in ("sensitivity", "specificity"):
                self.assertEqual(self.metric(B, learner)[f"reading_{m}"], "not read: no lineage-split estimate")

    def test_reading_rule_at_its_edges(self):
        self.assertEqual(S.reading(0.9, 0.7, 0.9), "consistent")
        self.assertEqual(S.reading(0.7, 0.7, 0.9), "consistent")
        self.assertEqual(S.reading(0.900001, 0.7, 0.9), "below")
        self.assertEqual(S.reading(0.699999, 0.7, 0.9), "above")

    def test_reading_tallies(self):
        got = {t["learner"]: (t["below"], t["consistent"], t["above"], t["not_read"])
               for t in self.res["reading_tallies"]}
        self.assertEqual(got, {"scm": (1, 1, 0, 2), "cart": (0, 1, 1, 2), "lasso": (1, 0, 1, 0),
                               "known_gene": (0, 1, 1, 0)})

    def test_dataset_counts(self):
        d = by(self.res["dataset_counts"], dataset=A)
        self.assertEqual((d["labelled_resistant"], d["labelled_susceptible"]), (62, 61))
        self.assertEqual((d["scored_resistant"], d["scored_susceptible"]), (60, 60))
        self.assertEqual((d["excluded_no_sequence"], d["excluded_species_check_failed"], d["excluded_not_counted"]),
                         (1, 1, 1))
        self.assertEqual(d["scored_borderline"], 1)
        self.assertTrue(d["reportable"])
        self.assertEqual(d["known_gene"], "scored")
        b = by(self.res["dataset_counts"], dataset=B)
        self.assertEqual((b["scored_resistant"], b["scored_susceptible"]), (25, 25))
        self.assertTrue(b["reportable"], "25 of each class is enough")
        self.assertFalse(b["lasso_complete"])
        self.assertEqual(b["known_gene"], "no rule published")

    def test_counts_only_under_25_of_a_class(self):
        res = S.score(F.inputs(F.world(drug_b_resistant=24)))
        d = by(res["dataset_counts"], dataset=B)
        self.assertEqual((d["scored_resistant"], d["scored_susceptible"]), (24, 25))
        self.assertFalse(d["reportable"])
        r = by(res["metrics"], dataset=B, learner="scm")
        self.assertEqual((r["tp"], r["fn"], r["fp"], r["tn"]), (21, 3, 3, 22))
        self.assertNotIn("sensitivity", r)
        self.assertEqual(r["reading_sensitivity"], "not read: counts only")
        self.assertEqual(r["reading_specificity"], "not read: counts only")

    def test_genome_level_record(self):
        res = self.res
        self.assertEqual((res["genomes_frozen"], res["genomes_scored"], res["genomes_borderline_scored"]),
                         (126, 123, 1))
        self.assertEqual(res["exclusions"], {"no sequence": 1, "species check failed": 1, "not counted": 1})
        g = {r["genome_id"]: r for r in res["genomes"]}
        self.assertEqual((g["1001.121"]["exclusion"], g["1001.122"]["exclusion"], g["1001.123"]["exclusion"]),
                         ("no sequence", "species check failed", "not counted"))
        self.assertTrue(g[F.BORDERLINE]["borderline"] and g[F.BORDERLINE]["scored"])
        self.assertEqual([r["genome_id"] for r in res["genomes"]][:3], ["1001.1", "1001.2", "1001.3"])
        self.assertEqual(res["genomes"][-1]["genome_id"], "1001.126", "numeric id order, not text order")
        self.assertNotIn("2001.1", g, "the same-study set is not in the primary population")

    def test_by_study_country_and_st(self):
        rows = [r for r in self.res["groups"] if r["dataset"] == A and r["learner"] == "scm"]
        bp = {r["group"]: (r["scored_resistant"], r["scored_susceptible"], r["tp"], r["fn"], r["fp"], r["tn"])
              for r in rows if r["grouping"] == "bioproject"}
        self.assertEqual(bp, {"P1": (30, 30, 27, 3, 3, 27), "P2": (20, 20, 16, 4, 4, 16),
                              "P3": (10, 10, 5, 5, 0, 10)})
        country = {r["group"]: (r["tp"], r["fn"], r["fp"], r["tn"]) for r in rows if r["grouping"] == "country"}
        self.assertEqual(country, {"Atlantis": (27, 3, 3, 27), "Lemuria": (21, 9, 4, 26)})
        # MLST: the candidates' mlst column against the cohort's STs
        # (MLST.testus.1 and .2); the fixture cycles .1, .7 and none by id.
        st = {r["group"]: r["scored_resistant"] + r["scored_susceptible"]
              for r in rows if r["grouping"] == "mlst_in_cohort"}
        want = {"yes": 0, "no": 0, "no ST": 0}
        for g in F.new_genomes():
            want[{"MLST.testus.1": "yes", "MLST.testus.7": "no", "": "no ST"}[g["mlst"]]] += 1
        self.assertEqual(st, want)

    def test_leave_one_study_out(self):
        rows = {r["left_out"]: r for r in self.res["leave_one_study_out"]
                if r["dataset"] == A and r["learner"] == "scm"}
        self.assertEqual(sorted(rows), ["P1", "P2", "P3"], "P4 has no labelled genome on drugA")
        self.assertAlmostEqual(rows["P1"]["sensitivity"], 21 / 30)
        self.assertAlmostEqual(rows["P2"]["sensitivity"], 32 / 40)
        self.assertAlmostEqual(rows["P3"]["sensitivity"], 43 / 50)
        self.assertAlmostEqual(rows["P1"]["specificity"], 26 / 30)
        self.assertAlmostEqual(rows["P2"]["specificity"], 37 / 40)
        self.assertAlmostEqual(rows["P3"]["specificity"], 43 / 50)
        self.assertEqual((rows["P1"]["scored_resistant"], rows["P1"]["scored_susceptible"]), (30, 30))
        r = self.metric(A, "scm")
        self.assertEqual(r["studies"], 3)
        self.assertEqual((r["loso_sensitivity_n"], r["loso_sensitivity_min"], r["loso_sensitivity_max"]),
                         (3, 21 / 30, 43 / 50))
        self.assertEqual((r["loso_specificity_min"], r["loso_specificity_max"]), (43 / 50, 37 / 40))
        self.assertNotIn("loso_sensitivity_lo", r, "no study-level interval")
        k = self.metric(A, "known_gene")
        self.assertEqual((k["loso_sensitivity_min"], k["loso_sensitivity_max"]), (20 / 30, 1.0))
        # drugB: the one study left out leaves nothing; no estimate, no range.
        b = self.metric(B, "scm")
        lb = by(self.res["leave_one_study_out"], dataset=B, learner="scm")
        self.assertEqual((lb["left_out"], lb["scored_resistant"], lb["reportable"]), ("P1", 0, False))
        self.assertEqual((b["loso_sensitivity_n"], b["loso_sensitivity_min"]), (0, None))

    def test_leave_one_study_out_needs_25_of_each_class_in_the_rest(self):
        w = F.world()
        # Put P2's genomes in P1: P1 then has 50 R and 50 S, P3 10 and 10.
        for c in w["candidates"]:
            if c["bioproject"] == "P2":
                c["bioproject"] = "P1"
        res = S.score(F.inputs(w))
        rows = {r["left_out"]: r for r in res["leave_one_study_out"] if r["dataset"] == A and r["learner"] == "scm"}
        self.assertEqual(sorted(rows), ["P1", "P3"])
        self.assertEqual((rows["P1"]["scored_resistant"], rows["P1"]["reportable"]), (10, False))
        self.assertNotIn("sensitivity", rows["P1"])
        self.assertEqual((rows["P3"]["scored_resistant"], rows["P3"]["reportable"]), (50, True))
        self.assertAlmostEqual(rows["P3"]["sensitivity"], 43 / 50)
        r = by(res["metrics"], dataset=A, learner="scm")
        self.assertEqual((r["studies"], r["loso_sensitivity_n"]), (2, 1))
        self.assertEqual((r["loso_sensitivity_min"], r["loso_sensitivity_max"]), (43 / 50, 43 / 50))

    def test_feature_carriage(self):
        f1 = by(self.res["features"], dataset=A, learner="scm", feature=1)
        self.assertEqual((f1["new_resistant_carriers"], f1["new_resistant"]), (48, 60))
        self.assertEqual((f1["new_susceptible_carriers"], f1["new_susceptible"]), (7, 60))
        self.assertAlmostEqual(f1["new_resistant_share"], 0.8)
        self.assertEqual((f1["cohort_resistant_carriers"], f1["cohort_resistant"]), (2, 4))
        self.assertEqual((f1["cohort_susceptible_carriers"], f1["cohort_susceptible"]), (1, 4))
        self.assertEqual(f1["kmer"], F.KMER[1])
        f2 = by(self.res["features"], dataset=A, learner="lasso", feature=2)
        self.assertEqual((f2["new_resistant_carriers"], f2["new_susceptible_carriers"]), (30, 30))
        self.assertEqual((f2["cohort_resistant_carriers"], f2["cohort_susceptible_carriers"]), (2, 0))
        self.assertFalse([r for r in self.res["features"] if r["learner"] == "known_gene"])
        self.assertFalse([r for r in self.res["features"] if r["dataset"] == B and r["learner"] == "cart"])

    def test_reading_text(self):
        text = S.reading_text(self.res)
        self.assertIn(S.NAME, text)
        for word in S.FORBIDDEN_WORDS:
            self.assertNotIn(word, text.lower())
        self.assertIn("scm: below 1, consistent 1, above 0, not read 2.", text)
        self.assertIn("Sensitivity 80.0% [67.7%, 89.2%], lineage split 95.0%, random split 85.0%: below.", text)
        self.assertIn("Genomes excluded, species check failed: 1.", text)
        self.assertNotIn("REHEARSAL", text)

    def test_reading_text_refuses_forbidden_words(self):
        res = dict(self.res, species="clin" + "ical_testus")
        with self.assertRaises(S.Refusal):
            S.reading_text(res)


if __name__ == "__main__":
    unittest.main()
