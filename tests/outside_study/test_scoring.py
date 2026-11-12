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
        for learner in ("scm", "cart"):
            r = self.metric(A, learner)
            self.assertEqual((r["tp"], r["fn"], r["fp"], r["tn"]), (48, 12, 7, 53), learner)
        # The lasso differs from the SCM on 1001.4 (R) and 1001.31 (S).
        r = self.metric(A, "lasso")
        self.assertEqual((r["tp"], r["fn"], r["fp"], r["tn"]), (47, 13, 6, 54))
        r = self.metric(A, "known_gene")
        self.assertEqual((r["tp"], r["fn"], r["fp"], r["tn"]), (50, 10, 0, 60))
        # drugB labels 1001.24 S and 1001.54 R, the other way round from drugA.
        r = self.metric(B, "scm")
        self.assertEqual((r["tp"], r["fn"], r["fp"], r["tn"]), (21, 4, 4, 21))
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
        self.assertEqual((r["tp"], r["fn"], r["fp"], r["tn"]), (20, 4, 4, 21))
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
        # 1001.4's "MLST.testus.-" names no ST.
        self.assertEqual(st, {"yes": 40, "no": 39, "no ST": 41})
        g = {r["genome_id"]: r for r in self.res["genomes"]}
        self.assertEqual((g[F.NO_ST_DASH]["mlst"], g[F.NO_ST_DASH]["mlst_in_cohort"]), ("MLST.testus.-", "no ST"))

    def test_no_st(self):
        for v in ("", "  ", "MLST.testus.-", "MLST.testus.", "MLST.testus. - "):
            self.assertTrue(S.no_st(v), repr(v))
        for v in ("MLST.testus.1", "MLST.testus.0", "MLST.testus.123"):
            self.assertFalse(S.no_st(v), v)

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
        # P1 holds the most genomes of both classes and can't be left out:
        # the one reportable fold (P3 left out) still includes it.
        a = by(res["dataset_counts"], dataset=A)
        self.assertEqual((a["top_study_resistant"], a["top_study_susceptible"], a["loso_folds_reportable"]),
                         ("P1", "P1", 1))
        self.assertEqual((a["top_study_resistant_fold_reportable"], a["top_study_susceptible_fold_reportable"]),
                         (False, False))
        text = S.reading_text(res)
        line = ("  With P1 left out, fewer than 25 resistant or 25 susceptible genomes remain, so no estimate can be "
                "shown without that study, and every leave-one-study-out range for this dataset includes it.")
        self.assertEqual(text.count(line), 1, "once for the dataset, not once per class")
        self.assertIn("Of 2 studies, 1 can be left out with 25 of each class remaining.\n" + line, text)

    def test_largest_study_of_one_class_only(self):
        w = F.world()
        # P2's susceptible genomes go to P3: P1 30 R and 30 S, P2 20 R, P3
        # 10 R and 30 S. P1 holds the most R, and P1 and P3 tie on S, so P1
        # is named for both. With P1 left out 30 R and 30 S remain.
        moved = {g["id"] for g in F.new_genomes() if g["project"] == "P2" and not g["resistant"]}
        self.assertEqual(len(moved), 20)
        for c in w["candidates"]:
            if c["genome_id"] in moved:
                c["bioproject"] = "P3"
        res = S.score(F.inputs(w))
        a = by(res["dataset_counts"], dataset=A)
        self.assertEqual((a["top_study_resistant"], a["top_study_susceptible"]), ("P1", "P1"))
        self.assertEqual((a["top_study_resistant_fold_reportable"], a["top_study_susceptible_fold_reportable"]),
                         (True, True))
        self.assertNotIn("no estimate can be shown without that study", S.reading_text(res))

    def test_feature_carriage(self):
        f1 = by(self.res["features"], dataset=A, learner="scm", feature=1)
        self.assertEqual((f1["new_resistant_carriers"], f1["new_resistant"]), (48, 60))
        self.assertEqual((f1["new_susceptible_carriers"], f1["new_susceptible"]), (7, 60))
        self.assertAlmostEqual(f1["new_resistant_share"], 0.8)
        self.assertEqual((f1["cohort_resistant_carriers"], f1["cohort_resistant"]), (2, 4))
        self.assertEqual((f1["cohort_susceptible_carriers"], f1["cohort_susceptible"]), (1, 4))
        self.assertEqual(f1["kmer"], F.KMER[1])
        f2 = by(self.res["features"], dataset=A, learner="lasso", feature=2)
        self.assertEqual((f2["new_resistant_carriers"], f2["new_susceptible_carriers"]), (1, 1))
        self.assertEqual((f2["cohort_resistant_carriers"], f2["cohort_susceptible_carriers"]), (2, 0))
        self.assertFalse([r for r in self.res["features"] if r["learner"] == "known_gene"])
        self.assertFalse([r for r in self.res["features"] if r["dataset"] == B and r["learner"] == "cart"])

    def test_studies_per_dataset(self):
        a = by(self.res["dataset_counts"], dataset=A)
        self.assertEqual((a["top_study_resistant"], a["top_study_resistant_count"], a["top_study_resistant_share"]),
                         ("P1", 30, 0.5))
        self.assertEqual((a["top_study_susceptible"], a["top_study_susceptible_count"]), ("P1", 30))
        self.assertEqual((a["studies"], a["loso_folds_reportable"]), (3, 3))
        self.assertEqual((a["top_study_resistant_fold_reportable"], a["top_study_susceptible_fold_reportable"]),
                         (True, True))
        b = by(self.res["dataset_counts"], dataset=B)
        self.assertEqual((b["top_study_resistant"], b["top_study_resistant_share"],
                          b["top_study_susceptible_share"]), ("P1", 1.0, 1.0))
        self.assertEqual((b["studies"], b["loso_folds_reportable"]), (1, 0))
        self.assertEqual((b["top_study_resistant_fold_reportable"], b["top_study_susceptible_fold_reportable"]),
                         (False, False))
        self.assertEqual((a["studies_not_left_out"], b["studies_not_left_out"]), ("", "P1"))
        self.assertEqual((a["scored_no_bioproject"], b["scored_no_bioproject"]), (0, 0))
        text = S.reading_text(self.res)
        self.assertNotIn("no estimate can be shown without that study", text,
                         "a dataset with one study gets the one-study sentence only")
        self.assertNotIn("no BioProject recorded", text)
        self.assertNotIn("measured mostly on different studies", text)
        self.assertIn("  By study: P1 holds 30 of the 60 resistant genomes (50.0%); P1 holds 30 of the 60 susceptible "
                      "genomes (50.0%). Of 3 studies, 3 can be left out with 25 of each class remaining.", text)
        self.assertIn("  By study: P1 holds 25 of the 25 resistant genomes (100.0%); P1 holds 25 of the 25 susceptible "
                      "genomes (100.0%). Of 1 study, 0 can be left out with 25 of each class remaining.\n"
                      "  Every scored genome is from one study, so the estimate is that study's "
                      "and no study-to-study range can be shown.", text)

    def restudied(self, study, **world_args):
        """The world with each labelled genome's BioProject given by
        study(genome, its place among the genomes of its class)."""
        w = F.world(**world_args)
        place, new = {True: 0, False: 0}, {}
        for g in F.new_genomes():
            new[g["id"]] = study(g, place[g["resistant"]])
            place[g["resistant"]] += 1
        for c in w["candidates"]:
            c["bioproject"] = new.get(c["genome_id"], c["bioproject"])
        return w

    def test_two_studies_and_no_fold_is_not_one_study(self):
        # X holds 36 R and 24 S, Y 24 R and 36 S. Neither can be left out.
        w = self.restudied(lambda g, i: "X" if i < (36 if g["resistant"] else 24) else "Y")
        res = S.score(F.inputs(w))
        a = by(res["dataset_counts"], dataset=A)
        self.assertEqual((a["top_study_resistant"], a["top_study_resistant_count"]), ("X", 36))
        self.assertEqual((a["top_study_susceptible"], a["top_study_susceptible_count"]), ("Y", 36))
        self.assertEqual((a["studies"], a["loso_folds_reportable"], a["studies_not_left_out"]), (2, 0, "X;Y"))
        text = S.reading_text(res)
        self.assertNotIn("is from one study", text, "two studies are never called one")
        self.assertIn("Of 2 studies, 0 can be left out with 25 of each class remaining.\n"
                      "  More than half of the resistant genomes are from X and more than half of the susceptible "
                      "genomes from Y, so sensitivity and specificity are measured mostly on different studies.\n"
                      "  No leave-one-study-out fold is reportable, so no study-to-study range can be shown.\n", text)

    def test_exactly_half_is_not_more_than_half(self):
        # On drugA, X holds 30 of the 60 R genomes, exactly half, and Y 40
        # of the 60 S genomes. "More than half" is true of Y only, so
        # nothing is said.
        def study(r_in_x):
            return lambda g, i: (("X" if i < r_in_x else "Y" if i < 50 else "Z") if g["resistant"]
                                 else ("Y" if i < 40 else "X"))

        def drug_a(res):
            return S.reading_text(res).split(f"{A}:")[1].split(f"{B}:")[0]
        res = S.score(F.inputs(self.restudied(study(30))))
        a = by(res["dataset_counts"], dataset=A)
        self.assertEqual((a["top_study_resistant"], a["top_study_resistant_count"], a["top_study_resistant_share"]),
                         ("X", 30, 0.5))
        self.assertEqual((a["top_study_susceptible"], a["top_study_susceptible_count"]), ("Y", 40))
        self.assertNotIn("measured mostly on different studies", drug_a(res))
        # With one more R genome in X, 31 of 60, it is said.
        res = S.score(F.inputs(self.restudied(study(31))))
        a = by(res["dataset_counts"], dataset=A)
        self.assertEqual((a["top_study_resistant"], a["top_study_resistant_count"]), ("X", 31))
        self.assertIn("  More than half of the resistant genomes are from X and more than half of the susceptible "
                      "genomes from Y, so sensitivity and specificity are measured mostly on different studies.\n",
                      drug_a(res))

    def test_every_study_that_cannot_be_left_out_is_named(self):
        # drugA keeps 60 R and 40 S labels (P2's 20 S genomes become I).
        # R: 20 each in X, Y and Z. S: 17 in Y, 16 in Z, 7 in X. Without Y
        # 23 S remain and without Z 24, so neither can be left out; without
        # X 33 remain. X is named for the most R (a tie, first name), Y for
        # the most S. Z is the largest study of neither class.
        def study(g, i):
            if g["resistant"]:
                return "XYZ"[i // 20]
            if g["project"] == "P2":
                return "P2"
            j = i if g["project"] == "P1" else i - 20
            return "Y" if j < 17 else "Z" if j < 33 else "X"
        w = self.restudied(study)
        unlabelled = {g["id"] for g in F.new_genomes() if g["project"] == "P2" and not g["resistant"]}
        for p in w["pairs"]:
            if p["genome_id"] in unlabelled and p["dataset_id"] == A:
                p.update(call="I", records_S=0, records_I=1, in_frozen_test="no")
        res = S.score(F.inputs(w))
        a = by(res["dataset_counts"], dataset=A)
        self.assertEqual((a["scored_resistant"], a["scored_susceptible"]), (60, 40))
        self.assertEqual((a["top_study_resistant"], a["top_study_resistant_count"]), ("X", 20))
        self.assertEqual((a["top_study_susceptible"], a["top_study_susceptible_count"]), ("Y", 17))
        self.assertEqual((a["top_study_resistant_fold_reportable"], a["top_study_susceptible_fold_reportable"]),
                         (True, False))
        self.assertEqual((a["studies"], a["loso_folds_reportable"], a["studies_not_left_out"]), (3, 1, "Y;Z"))
        self.assertEqual((a["top_study_resistant_tied"], a["top_study_susceptible_tied"]), ("Y;Z", ""))
        text = S.reading_text(res)
        self.assertIn("  By study: X holds 20 of the 60 resistant genomes (33.3%), as many as Y and Z; "
                      "Y holds 17 of the 40 susceptible genomes (42.5%). Of 3 studies, 1 can be left out", text)
        for b in "YZ":
            self.assertEqual(text.count(f"  With {b} left out, fewer than 25 resistant or 25 susceptible genomes "
                                        "remain, so no estimate can be shown without that study"), 1, b)
        self.assertNotIn("With X left out", text)

    def test_a_tie_for_the_most_genomes_names_every_study(self):
        # Y and Z each hold 30 R and 30 S: Y is named, and Z beside it.
        w = self.restudied(lambda g, i: "Y" if i < 30 else "Z")
        res = S.score(F.inputs(w))
        a = by(res["dataset_counts"], dataset=A)
        self.assertEqual((a["top_study_resistant"], a["top_study_resistant_tied"]), ("Y", "Z"))
        self.assertEqual((a["top_study_susceptible"], a["top_study_susceptible_tied"]), ("Y", "Z"))
        self.assertIn("  By study: Y holds 30 of the 60 resistant genomes (50.0%), as many as Z; Y holds 30 of the "
                      "60 susceptible genomes (50.0%), as many as Z. Of 2 studies, 2 can be left out",
                      S.reading_text(res))
        # No tie in the plain world.
        a = by(self.res["dataset_counts"], dataset=A)
        self.assertEqual((a["top_study_resistant_tied"], a["top_study_susceptible_tied"]), ("", ""))
        self.assertNotIn("as many as", S.reading_text(self.res))

    def test_genomes_with_no_bioproject_are_counted_and_said(self):
        # P2's and P3's genomes lose their BioProject: one study, P1, and
        # 60 genomes in none. Leaving P1 out leaves those 60 (30 R, 30 S).
        w = self.restudied(lambda g, i: "P1" if g["project"] == "P1" else "")
        res = S.score(F.inputs(w))
        a = by(res["dataset_counts"], dataset=A)
        self.assertEqual((a["studies"], a["loso_folds_reportable"], a["scored_no_bioproject"]), (1, 1, 60))
        self.assertEqual(a["studies_not_left_out"], "")
        text = S.reading_text(res)
        self.assertIn("Of 1 study, 1 can be left out with 25 of each class remaining.\n"
                      "  Scored genomes with no BioProject recorded, counted in no study: 60.\n", text)
        block = text.split(f"{A}:")[1].split(f"{B}:")[0]
        self.assertNotIn("is from one study", block, "60 of the 120 genomes are in no study")

    def test_an_excluded_genome_in_no_study_is_not_a_scored_one(self):
        # The three labelled genomes that are excluded (no sequence, the
        # species check, not counted) lose their BioProject. None of them is
        # scored, so no scored genome is in no study.
        w = F.world()
        gone = {g for g, e in F.EXTRA.items() if e["label"] in ("R", "S")}
        self.assertEqual(len(gone), 3)
        for c in w["candidates"]:
            if c["genome_id"] in gone:
                c["bioproject"] = ""
        res = S.score(F.inputs(w))
        a = by(res["dataset_counts"], dataset=A)
        self.assertEqual((a["labelled_resistant"], a["labelled_susceptible"]), (62, 61))
        self.assertEqual((a["scored_resistant"], a["scored_susceptible"], a["scored_no_bioproject"]), (60, 60, 0))
        self.assertNotIn("counted in no study", S.reading_text(res))

    def test_a_counts_only_dataset_gets_no_sentence_about_its_estimate(self):
        res = S.score(F.inputs(F.world(drug_b_resistant=24)))
        self.assertFalse(by(res["dataset_counts"], dataset=B)["reportable"])
        block = S.reading_text(res).split(f"{B}:")[1]
        self.assertIn("Of 1 study, 0 can be left out with 25 of each class remaining.", block)
        self.assertNotIn("the estimate", block)
        self.assertNotIn("study-to-study range", block)

    def test_reference_kind(self):
        for learner in ("scm", "cart", "lasso"):
            r = self.metric(A, learner)
            self.assertEqual(r["reference_kind"], "cross-validation, lineage split")
            self.assertFalse(r["reference_at_0_or_1"], learner)
        # drugB's report has no lineage split, and its reference says so.
        for learner in ("scm", "cart"):
            r = self.metric(B, learner)
            self.assertIn("unavailable", r["lineage_note"])
            self.assertEqual(r["reference_kind"], "cross-validation, random split only (lineage split unavailable)")
        self.assertIn("    Reference figure: cross-validation, random split only (lineage split unavailable)", S.reading_text(self.res))
        k = self.metric(A, "known_gene")
        self.assertEqual(k["reference_kind"], "fixed rule, whole cohort")
        self.assertTrue(k["reference_at_0_or_1"], "the catalog's lineage specificity is 1.0")
        text = S.reading_text(self.res)
        self.assertIn("    Reference figure: fixed rule, whole cohort; a reference sensitivity or specificity "
                      "is exactly 0 or 1.", text)
        self.assertIn("    Reference figure: cross-validation, lineage split.", text)

    def test_paired_with_the_known_gene_rule(self):
        got = {r["learner"]: r for r in self.res["paired"]}
        self.assertEqual(sorted(got), ["cart", "lasso", "scm"], "drugB has no known-gene rule")
        cols = ("genomes", "both", "learner_only", "known_gene_only", "neither")
        for learner, r_want, s_want in (("scm", (60, 43, 5, 7, 5), (60, 0, 7, 0, 53)),
                                        ("cart", (60, 43, 5, 7, 5), (60, 0, 7, 0, 53)),
                                        ("lasso", (60, 42, 5, 8, 5), (60, 0, 6, 0, 54))):
            r = got[learner]
            self.assertEqual(tuple(r[f"resistant_{c}"] for c in cols), r_want, learner)
            self.assertEqual(tuple(r[f"susceptible_{c}"] for c in cols), s_want, learner)

    def test_duplicates_by_kmer_set(self):
        got = [(r["genome_id"], r["other_genome_id"], r["other_is"]) for r in self.res["duplicates"]]
        self.assertEqual(got, [("1001.5", "999.2", "cohort"), ("1001.7", "1001.8", "scored new genome"),
                               ("1001.8", "1001.7", "scored new genome")])
        self.assertEqual(self.res["duplicate_counts"],
                         {"scored genomes with the k-mer set of a cohort genome": 1,
                          "scored genomes with the k-mer set of another scored genome": 2})
        self.assertEqual(self.res["genomes_scored"], 123, "a duplicate is listed, never excluded")

    def test_duplicate_counts_are_genomes_not_rows(self):
        # 1001.9 gets the k-mer set of 1001.7 and 1001.8: three genomes, six rows.
        inputs = F.inputs(F.world())
        inputs["digests"]["new"]["1001.9"] = inputs["digests"]["new"]["1001.7"]
        res = S.score(inputs)
        self.assertEqual(sum(1 for r in res["duplicates"] if r["other_is"] == "scored new genome"), 6)
        self.assertEqual(res["duplicate_counts"]["scored genomes with the k-mer set of another scored genome"], 3)

    def test_the_0_or_1_flag_reads_both_splits(self):
        # Only the random split's figure is exactly 1 here.
        w = F.world()
        split = by(w["reports"][A]["splits"], scheme="random")
        by(split["learners"], learner="scm")["metrics"]["specificity"] = 1.0
        w["shas"][A] = F.sha(F.report_bytes(w["reports"][A]))
        for r in w["rows"][A] + w["replay_rows"][A]:
            r["report_sha256"] = w["shas"][A]
        res = S.score(F.inputs(w))
        self.assertTrue(by(res["metrics"], dataset=A, learner="scm")["reference_at_0_or_1"])
        self.assertFalse(by(res["metrics"], dataset=A, learner="cart")["reference_at_0_or_1"])

    def test_reading_text(self):
        text = S.reading_text(self.res)
        self.assertIn(S.NAME, text)
        for word in S.FORBIDDEN_WORDS:
            self.assertNotIn(word, text.lower())
        self.assertIn("scm: below 1, consistent 1, above 0, not read 2.", text)
        self.assertIn("Sensitivity 80.0% [67.7%, 89.2%], lineage split 95.0%, random split 85.0%: below.", text)
        self.assertIn("Genomes excluded, species check failed: 1.", text)
        self.assertNotIn("REHEARSAL", text)
        for sentence in (
            "Sensitivity and specificity are over the genomes a learner called; genomes with no call are counted "
            "in their own column.",
            "The datasets share genomes and studies, so they are not independent of each other.",
            "Precision depends on the share of resistant genomes, which differs between the cohort and the new genomes.",
            "Scored genomes with the same k-mer set as a cohort genome: 1; as another scored genome: 2",
            "Stage attempts: setup 1, fetch 1, species-check 1, count 1, known-gene 1, predict 1.",
        ):
            self.assertIn(sentence, text)

    def test_forbidden_words(self):
        for word in ("validated", "generalizes", "clinical", "proves", "proven"):
            self.assertIn(word, S.FORBIDDEN_WORDS)
            res = dict(self.res, species=word + "_testus")
            with self.assertRaises(S.Refusal):
                S.reading_text(res)

    def test_reading_text_refuses_forbidden_words(self):
        res = dict(self.res, species="clin" + "ical_testus")
        with self.assertRaises(S.Refusal):
            S.reading_text(res)



class Variants(unittest.TestCase):
    """Worlds changed in one way each, so that a fault in one place of the
    scoring changes a number these tests know."""

    def test_known_gene_calls_are_read_per_dataset(self):
        # drugB's catalog rule calls every genome resistant; drugA's does not.
        res = S.score(F.inputs(F.world(b_catalog=True)))
        r = by(res["metrics"], dataset=B, learner="known_gene")
        self.assertEqual((r["tp"], r["fn"], r["fp"], r["tn"]), (25, 0, 25, 0))
        r = by(res["metrics"], dataset=A, learner="known_gene")
        self.assertEqual((r["tp"], r["fn"], r["fp"], r["tn"]), (50, 10, 0, 60))

    def test_counts_only_is_decided_after_the_exclusions(self):
        # 25 R labelled on drugB, one of them with no sequence: 24 scored.
        res = S.score(F.inputs(F.world(drug_b_resistant=24, b_extra_label="1001.121")))
        d = by(res["dataset_counts"], dataset=B)
        self.assertEqual((d["labelled_resistant"], d["scored_resistant"], d["excluded_no_sequence"]), (25, 24, 1))
        self.assertFalse(d["reportable"])
        self.assertEqual(by(res["metrics"], dataset=B, learner="scm")["reading_sensitivity"], "not read: counts only")

    def test_a_fold_needs_25_of_each_class_left(self):
        # Six of P2's S genomes moved to P1: leaving P1 out leaves 30 R, 24 S.
        w = F.world()
        moved = {f"1001.{i}" for i in range(81, 87)}
        for c in w["candidates"]:
            if c["genome_id"] in moved:
                c["bioproject"] = "P1"
        res = S.score(F.inputs(w))
        rows = {r["left_out"]: r for r in res["leave_one_study_out"] if r["dataset"] == A and r["learner"] == "scm"}
        self.assertEqual((rows["P1"]["scored_resistant"], rows["P1"]["scored_susceptible"]), (30, 24))
        self.assertFalse(rows["P1"]["reportable"])
        r = by(res["metrics"], dataset=A, learner="scm")
        self.assertEqual(r["loso_sensitivity_n"], 2)
        self.assertEqual(by(res["dataset_counts"], dataset=A)["loso_folds_reportable"], 2)

    def test_no_calls_are_counted_apart(self):
        res = S.score(F.inputs(F.world(lasso_empty=True)))
        r = by(res["metrics"], dataset=A, learner="lasso")
        self.assertEqual((r["tp"], r["fn"], r["fp"], r["tn"], r["no_call_resistant"], r["no_call_susceptible"]),
                         (0, 0, 0, 0, 60, 60))
        self.assertEqual(r["reading_sensitivity"], "not read: no interval")
        self.assertIsNone(r["sensitivity"])
        self.assertFalse([p for p in res["paired"] if p["learner"] == "lasso" and p["resistant_genomes"]])

    def test_an_empty_bioproject_or_country_is_never_a_study(self):
        w = F.world()
        blank = {f"1001.{i}" for i in range(91, 121)}  # all of P3
        for c in w["candidates"]:
            if c["genome_id"] in blank:
                c["bioproject"], c["isolation_country"] = "", " "
        res = S.score(F.inputs(w))
        rows = [r for r in res["leave_one_study_out"] if r["dataset"] == A and r["learner"] == "scm"]
        self.assertEqual(sorted(r["left_out"] for r in rows), ["P1", "P2"])
        self.assertEqual(by(res["metrics"], dataset=A, learner="scm")["studies"], 2)
        groups = {(r["grouping"], r["group"]) for r in res["groups"] if r["dataset"] == A and r["learner"] == "scm"}
        self.assertIn(("bioproject", "no BioProject recorded"), groups)
        self.assertIn(("country", "no country recorded"), groups)
        g = {r["genome_id"]: r for r in res["genomes"]}
        self.assertEqual((g["1001.100"]["bioproject"], g["1001.100"]["country"]),
                         ("no BioProject recorded", "no country recorded"))
        self.assertEqual(by(res["dataset_counts"], dataset=A)["studies"], 2)


if __name__ == "__main__":
    unittest.main()
