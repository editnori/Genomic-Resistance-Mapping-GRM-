"""Every refusal of the scoring: rows from another contract, revision or
report; genomes outside the frozen list or excluded; missing, duplicate or
inconsistent rows; rows for learners or datasets that are not scored; and
the known-gene record."""

import copy
import json
import unittest

from tests.outside_study import fixtures as F
from tests.outside_study.fixtures import S, A, B, C


class Refusals(unittest.TestCase):
    def setUp(self):
        self.w = F.world()

    def refused(self, pattern, rehearsal=False):
        with self.assertRaisesRegex(S.Refusal, pattern):
            S.score(F.inputs(self.w, rehearsal=rehearsal))

    def row(self, d=A, genome="1001.2", learner="scm"):
        return next(r for r in self.w["rows"][d] if r["genome"] == genome and r["learner"] == learner)

    def test_the_world_itself_scores(self):
        S.score(F.inputs(self.w))

    def test_another_contract(self):
        self.row()["contract"] = "grm-predict-2"
        self.refused("contract")

    def test_another_revision(self):
        self.row()["evaluation_revision"] = 11
        self.refused("revision 11")

    def test_a_row_naming_another_report(self):
        self.row()["report_sha256"] = "0" * 64
        self.refused("names report sha256")

    def test_a_row_for_another_dataset(self):
        self.row()["dataset"] = B
        self.refused("a row for dataset")

    def test_a_row_at_another_k(self):
        self.row()["k"] = F.K + 1
        self.refused("another k")

    def test_predicted_genomes_that_are_not_the_survivors(self):
        inputs = F.inputs(self.w)
        inputs["predict_survivors"] = inputs["predict_survivors"][:-1]
        with self.assertRaisesRegex(S.Refusal, "survivors"):
            S.score(inputs)

    def test_a_row_with_another_fingerprint(self):
        self.row()["input_fingerprint"] = "beef"
        self.refused("fingerprint")

    def test_a_genome_outside_the_frozen_list(self):
        r = copy.deepcopy(self.row())
        r["genome"] = "1001.999"
        self.w["rows"][A].append(r)
        self.refused("outside the frozen list")

    def test_a_requested_genome_outside_the_frozen_list(self):
        self.w["availability"]["1001.999"] = {"genome_id": "1001.999", "fasta": "no", "detail": ""}
        self.refused("outside the frozen list")

    def test_a_frozen_genome_never_requested(self):
        del self.w["availability"]["1001.121"]
        self.refused("never requested")

    def test_an_excluded_genome_with_a_prediction(self):
        r = copy.deepcopy(self.row())
        r["genome"] = "1001.122"
        self.w["rows"][A].append(r)
        self.refused("excluded genome 1001.122")

    def test_a_scored_genome_without_a_row(self):
        self.w["rows"][A].remove(self.row(learner="cart"))
        self.refused("no cart row for scored genome 1001.2")

    def test_two_rows_for_one_genome_and_learner(self):
        self.w["rows"][A].append(copy.deepcopy(self.row()))
        self.refused("two rows")

    def test_a_call_that_does_not_follow_from_the_features(self):
        r = self.row()
        r["call"] = "susceptible" if r["call"] == "resistant" else "resistant"
        self.refused("the model gives")

    def test_a_score_that_does_not_follow_from_the_features(self):
        self.row(learner="lasso")["score"] += 1e-6
        self.refused("score")

    def test_a_feature_list_that_is_not_the_model(self):
        self.row(learner="lasso")["features"].pop()
        self.refused("the model reads")

    def test_a_kmer_that_is_not_the_published_one(self):
        self.row()["features"][0]["kmer"] = "CCCCC"
        self.refused("not the published one")

    def test_rows_for_a_dataset_not_eligible(self):
        rep = self.w["reports"][C]
        self.w["rows"][C] = [F.make_row(C, self.w["shas"][C], "1001.2", lr, {1: True, 2: False})
                             for lr in ("scm", "cart", "lasso")]
        self.assertTrue(rep)
        self.refused("not eligible")

    def test_a_row_for_an_incomplete_lasso(self):
        r = copy.deepcopy(self.row(d=B))
        r["learner"] = "lasso"
        self.w["rows"][B].append(r)
        self.refused("learner 'lasso'")

    def test_a_row_for_a_forest(self):
        r = copy.deepcopy(self.row())
        r["learner"] = "rf"
        self.w["rows"][A].append(r)
        self.refused("learner 'rf'")

    def test_known_gene_rule_not_the_reports(self):
        self.w["known_gene"]["datasets"][A]["rule"] = json.dumps([["C:QUINOLONE"]])
        self.refused("not the report's catalog.rule")

    def test_known_gene_call_missing(self):
        del self.w["known_gene"]["calls"][(A, "1001.2")]
        self.refused("no known-gene call")

    def test_known_gene_call_for_an_excluded_genome(self):
        self.w["known_gene"]["calls"][(A, "1001.122")] = {"catalog_resistant": "no"}
        self.refused("excluded genome 1001.122")

    def test_known_gene_call_outside_the_list(self):
        self.w["known_gene"]["calls"][(A, "1001.999")] = {"catalog_resistant": "no"}
        self.refused("outside the frozen list")

    def test_known_gene_record_missing_a_dataset(self):
        del self.w["known_gene"]["datasets"][B]
        self.refused("no row for")

    def test_an_eligible_dataset_without_a_report(self):
        inputs = F.inputs(self.w)
        del inputs["reports"][B]
        with self.assertRaisesRegex(S.Refusal, "no report for eligible dataset"):
            S.score(inputs)

    def test_rehearsal_rows_must_be_the_published_bits(self):
        # In a rehearsal the genomes are cohort samples: their feature values
        # must be their published bits, and their calls the final models'.
        self.w = F.rehearsal_world()
        res = S.score(F.inputs(self.w, rehearsal=True))
        self.assertTrue(res["rehearsal"])
        self.assertIn("REHEARSAL ON COHORT GENOMES: MEASURES NOTHING.", S.reading_text(res))
        r = next(r for r in self.w["rows"][A] if r["genome"] == "999.5" and r["learner"] == "scm")
        r["features"][0]["found"] = False
        r["call"] = "susceptible"
        self.refused("rehearsal replay", rehearsal=True)

    def test_rehearsal_label_on_a_report_the_genome_is_not_in(self):
        # Outside a rehearsal no such check runs; inside one, a labelled
        # stand-in genome must be a sample of the report.
        with self.assertRaisesRegex(S.Refusal, "not one of its samples"):
            S.score(F.inputs(self.w, rehearsal=True))


class SameStudy(unittest.TestCase):
    def test_scored_alone_with_its_own_datasets(self):
        w = F.world()
        ids = F.SAME_STUDY
        w["availability"] = {g: {"genome_id": g, "fasta": "yes", "detail": ""} for g in ids}
        w["verdicts"] = {g: {"genome_id": g, "verdict": "pass", "reason": ""} for g in ids}
        w["counted"] = {g: {"genome_id": g, "counted": "yes", "detail": "", "kset_digest": F.digest_of(g)}
                        for g in ids}
        w["rows"] = {A: [F.make_row(A, w["shas"][A], g, lr, {1: i % 2 == 0, 2: False})
                         for i, g in enumerate(ids) for lr in ("scm", "cart", "lasso")]}
        w["known_gene"]["calls"] = {(A, g): {"catalog_resistant": "no"} for g in ids}
        w["survivors"] = list(ids)
        res = S.score(F.inputs(w, population="same-study"))
        self.assertEqual(res["datasets"], [A])
        self.assertEqual(res["population"], "same-study")
        r = next(r for r in res["metrics"] if r["learner"] == "scm")
        self.assertEqual((r["tp"], r["fn"], r["fp"], r["tn"], r["reportable"]), (2, 0, 0, 2, False))
        # Rows for drugB, eligible only in the primary population, are refused.
        w["rows"][B] = []
        with self.assertRaisesRegex(S.Refusal, "not eligible in the same-study population"):
            S.score(F.inputs(w, population="same-study"))


class Reports(unittest.TestCase):
    def setUp(self):
        import tempfile
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name
        self.report = F.make_report(A)
        self.data = F.report_bytes(self.report)

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, report=None):
        data = self.data if report is None else F.report_bytes(report)
        with open(f"{self.dir}/{A}.json", "wb") as f:
            f.write(data)
        return F.sha(data)

    def test_the_pinned_report_loads(self):
        sha = self.write()
        r, got = S.load_report(self.dir, A, {f"{A}.json": sha})
        self.assertEqual((r["id"], got), (A, sha))

    def test_a_report_that_is_not_the_pinned_one(self):
        self.write()
        with self.assertRaisesRegex(S.Refusal, "is not the pinned"):
            S.load_report(self.dir, A, {f"{A}.json": "0" * 64})

    def test_a_report_without_a_pin(self):
        self.write()
        with self.assertRaisesRegex(S.Refusal, "no pinned report"):
            S.load_report(self.dir, A, {})

    def test_another_revision_or_an_invalid_evaluation(self):
        for ev in ({"revision": 11, "valid": True}, {"revision": 12, "valid": False}):
            r = copy.deepcopy(self.report)
            r["evaluation"] = ev
            sha = self.write(r)
            with self.assertRaisesRegex(S.Refusal, "evaluation revision"):
                S.load_report(self.dir, A, {f"{A}.json": sha})

    def test_another_schema(self):
        r = copy.deepcopy(self.report)
        r["schema"] = 1
        sha = self.write(r)
        with self.assertRaisesRegex(S.Refusal, "schema"):
            S.load_report(self.dir, A, {f"{A}.json": sha})


if __name__ == "__main__":
    unittest.main()
