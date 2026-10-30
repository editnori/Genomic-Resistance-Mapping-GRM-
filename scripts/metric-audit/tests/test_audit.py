"""audit.py on synthetic runs: each fault must raise its own issue, and
each consistent run must audit clean.

Each case writes a /root/grm-shaped tree (report, counted.json, manifest
and published copy, backup, pipeline log) and audits it with a launch
time before the files were written.

usage: python3 -I tests/test_audit.py
"""
import argparse
import copy
import json
import os
import subprocess
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixtures as fx  # noqa: E402
from curation_cases import CurationCases

sys.path.insert(0, fx.LANE)
import audit  # noqa: E402

LAUNCH = "2026-01-01T00:00:00Z"
DEFAULT = object()


def ba(report, scheme, learner):
    return fx.row(report, scheme, learner)["metrics"]["balanced_accuracy"]


def run_log(report, winner="rf"):
    """The log grm writes for a recomputed dataset: the recompute reason, then the done line."""
    return [fx.recompute_line(), fx.done_line(winner, None if winner is None else ba(report, "random", winner))]


def empty_gbt(report, scheme="random"):
    gbt = fx.row(report, scheme, "gbt")
    gbt["scores"] = [None] * fx.N
    gbt["fold_complexity"] = [0.0] * 5
    fx.rescore(report)
    return gbt


def catalog_panel(report):
    calls = [1.0 if s["groups"] else 0.0 for s in report["samples"]]
    truth = [s["resistant"] for s in report["samples"]]
    m = fx.rates(truth, calls, False)
    return {"id": "fixture_panel", "name": "Fixture panel", "mirrors": "fixture", "targets": ["blaFIX"], "metrics": m,
            "detected": [["blaFIX", m["tp"], m["fp"]]]}


def panel_row(report, panel, best):
    return {"panel": panel["id"], "name": panel["name"], "dataset": report["id"], "species": report["species"],
            "antibiotic": report["antibiotic"], "metrics": fx.lite({"metrics": panel["metrics"]}),
            "best_model": fx.lite(fx.row(report, "random", best)), "best_learner": best}


class Case(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fg = audit.FoldGen()

    @classmethod
    def tearDownClass(cls):
        cls.fg.p.stdin.close()
        cls.fg.p.wait()
        cls.fg.p.stdout.close()

    def run_audit(self, report, index=None, backup=DEFAULT, log_lines=None, revision=7, launch=LAUNCH,
                  reshaped=False, published=True, digests=None):
        if backup is DEFAULT:
            backup = fx.make_report(revision=6, forest_variant=1)
        if log_lines is None:
            log_lines = run_log(report)
        with tempfile.TemporaryDirectory() as tmp:
            paths = fx.audit_root(tmp, report, index, backup, log_lines, published, digests)
            args = argparse.Namespace(root=paths["root"], revision=revision, reshaped=reshaped, foldgen=audit.FOLDGEN,
                                      backup=None if reshaped else paths["backup"],
                                      log=None if reshaped else paths["log"], launch=None if reshaped else launch)
            audit.species_finished.clear()
            audit.configure(args)
            done, recomputed = audit.done_datasets()
            return audit.audit_file(paths["report"], self.fg, {}, done, recomputed, audit.load_manifest(paths["root"]))

    def assertClean(self, rec):
        self.assertEqual(rec["issues"], [])

    def assertIssue(self, rec, *texts):
        dumped = [json.dumps(i) for i in rec["issues"]]
        for text in texts:
            self.assertTrue(any(text in d for d in dumped), "no issue with %r in %s" % (text, rec["issues"]))

    def assertNoIssue(self, rec, text):
        self.assertFalse(any(text in json.dumps(i) for i in rec["issues"]), rec["issues"])


class ConsistentRuns(Case):
    def test_a_revision7_run_audits_clean(self):
        rep = fx.make_report()
        rec = self.run_audit(rep)
        self.assertClean(rec)
        self.assertEqual(sorted((f["scheme"], f["learner"], f["fold"], f["fold_ba"], f["complexity"]) for f in rec["null_folds"]),
                         [("lineage", "gbt", 2, None, 0.0), ("random", "lasso", 4, None, 0.0)])
        self.assertEqual((rec["best"]["learner"], rec["manifest_compared"], rec["rf_cards_checked"]), ("rf", True, 3))
        self.assertEqual(rec["log_done"]["form"], "full coverage")
        self.assertIn("random/rf", rec["vs_backup"]["learners_with_changed_scores"])
        self.assertIn("lineage/rf", rec["vs_backup"]["learners_with_changed_scores"])

    def test_a_row_without_any_call_audits_clean_with_null_manifest_rates(self):
        rep = fx.make_report()
        empty_gbt(rep)
        index = fx.manifest(rep, "rf")
        self.assertIsNone(index["datasets"][0]["learners"][3]["random"]["balanced_accuracy"])
        rec = self.run_audit(rep, index)
        self.assertClean(rec)
        gbt = next(x for x in rec["learners"] if (x["scheme"], x["learner"]) == ("random", "gbt"))
        self.assertEqual((gbt["called"], gbt["balanced_accuracy"], gbt["no_call"]), (0, None, {"resistant": 40, "susceptible": 40}))

    def test_a_dataset_new_to_the_run_has_no_recompute_line(self):
        rep = fx.make_report()
        self.assertClean(self.run_audit(rep, backup=None, log_lines=run_log(rep)[1:]))

    def test_a_reshaped_revision6_report_audits_without_a_run(self):
        rep = fx.make_report()
        rep["evaluation"]["revision"] = 6
        del rep["eval_config"]
        index = fx.manifest(rep, "rf", evaluation_revision=6)
        self.assertClean(self.run_audit(rep, index, revision=6, reshaped=True))

    def test_the_all_half_null_fold_test_still_holds_for_schema1(self):
        rep = fx.make_report(revision=6)
        index = fx.manifest(rep, "gbt", schema=1, evaluation_revision=6)
        log = [fx.old_done_line("gbt", ba(rep, "random", "gbt"))]
        rec = self.run_audit(rep, index, backup=None, log_lines=log, revision=6)
        self.assertClean(rec)
        self.assertEqual(sorted((f["scheme"], f["learner"], f["fold"], f["fold_ba"], f["complexity"]) for f in rec["null_folds"]),
                         [("lineage", "gbt", 2, 0.5, 0.0), ("random", "lasso", 4, 0.5, 0.0)])

    def test_a_fitted_model_scoring_a_whole_fold_half_is_a_call(self):
        rep = fx.make_report()
        cart = fx.row(rep, "random", "cart")
        cart["scores"] = [0.5 if fx.FOLDS["random"][i] == 1 else s for i, s in enumerate(cart["scores"])]
        fx.rescore(rep)
        rec = self.run_audit(rep, fx.manifest(rep, "rf"))
        self.assertClean(rec)
        self.assertEqual([(f["learner"], f["fold"]) for f in rec["half_called_folds"]], [("cart", 1)])
        self.assertEqual(len(rec["null_folds"]), 2)


class NullScores(Case):
    def test_no_call_one_short(self):
        rep = fx.make_report()
        fx.row(rep, "lineage", "gbt")["no_call"]["resistant"] -= 1
        self.assertIssue(self.run_audit(rep), "lineage/gbt/no_call", "no_call differs from the null scores")

    def test_no_call_missing_beside_null_scores(self):
        rep = fx.make_report()
        del fx.row(rep, "random", "lasso")["no_call"]
        self.assertIssue(self.run_audit(rep), "random/lasso/no_call")

    def test_a_null_score_in_a_schema1_report(self):
        rep = fx.make_report(revision=6)
        fx.row(rep, "random", "rf")["scores"][0] = None
        rec = self.run_audit(rep, fx.manifest(rep, "gbt", schema=1, evaluation_revision=6), backup=None,
                             log_lines=[fx.old_done_line("gbt", ba(rep, "random", "gbt"))], revision=6)
        self.assertIssue(rec, "null score in a row that calls every isolate")

    def test_a_null_catalog_score(self):
        rep = fx.make_report()
        fx.row(rep, "random", "catalog")["scores"][0] = None
        self.assertIssue(self.run_audit(rep), "null score in a row that calls every isolate")

    def test_a_null_fold_with_a_nonzero_complexity(self):
        rep = fx.make_report()
        fx.row(rep, "lineage", "gbt")["fold_complexity"][2] = 5.0
        self.assertIssue(self.run_audit(rep), "lineage/gbt/fold 2", "null fold without a zero gbt or lasso complexity")

    def test_a_zero_complexity_fold_that_was_scored(self):
        rep = fx.make_report()
        fx.row(rep, "random", "gbt")["fold_complexity"][0] = 0.0
        self.assertIssue(self.run_audit(rep), "random/gbt/fold 0", "null fit scored")

    def test_an_rf_null_fold(self):
        rep = fx.make_report()
        rf = fx.row(rep, "random", "rf")
        rf["scores"] = [None if fx.FOLDS["random"][i] == 0 else s for i, s in enumerate(rf["scores"])]
        fx.rescore(rep)
        self.assertIssue(self.run_audit(rep), "random/rf/fold 0", "null fold without a zero gbt or lasso complexity")

    def test_a_fold_called_in_part(self):
        rep = fx.make_report()
        gbt = fx.row(rep, "lineage", "gbt")
        gbt["scores"][fx.FOLDS["lineage"].index(1)] = None
        fx.rescore(rep)
        self.assertIssue(self.run_audit(rep), "fold called in part")


class CalledOnlyMetrics(Case):
    def scored_half(self, rep, scheme, learner):
        row = fx.row(rep, scheme, learner)
        truth = [s["resistant"] for s in rep["samples"]]
        return fx.rates(truth, [0.5 if s is None else s for s in row["scores"]], True), row

    def test_n_over_every_isolate(self):
        rep = fx.make_report()
        fx.row(rep, "lineage", "gbt")["metrics"]["n"] = fx.N
        self.assertIssue(self.run_audit(rep), '"where": "lineage/gbt/pooled", "field": "n"')

    def test_auprc_with_uncalled_isolates_scored_half(self):
        rep = fx.make_report()
        half, row = self.scored_half(rep, "lineage", "gbt")
        row["metrics"]["auprc"] = half["auprc"]
        self.assertIssue(self.run_audit(rep), '"where": "lineage/gbt/pooled", "field": "auprc"')

    def test_auroc_with_uncalled_isolates_scored_half(self):
        rep = fx.make_report()
        half, row = self.scored_half(rep, "random", "lasso")
        row["metrics"]["auroc"] = half["auroc"]
        self.assertIssue(self.run_audit(rep), '"where": "random/lasso/pooled", "field": "auroc"')

    def test_roc_with_uncalled_isolates_scored_half(self):
        rep = fx.make_report()
        row = fx.row(rep, "lineage", "gbt")
        row["roc"] = fx.roc([s["resistant"] for s in rep["samples"]], [0.5 if s is None else s for s in row["scores"]])
        self.assertIssue(self.run_audit(rep), "lineage/gbt/roc")

    def test_zero_metrics_on_a_row_without_calls(self):
        rep = fx.make_report()
        gbt = empty_gbt(rep)
        gbt["metrics"] = dict(fx.row(rep, "random", "rf")["metrics"], n=0, tp=0, fp=0, tn=0, fn=0)
        self.assertIssue(self.run_audit(rep, fx.manifest(rep, "rf")), '"where": "random/gbt/pooled", "field": "metrics"')

    def test_half_fold_balanced_accuracy_on_a_null_fold(self):
        rep = fx.make_report()
        fx.row(rep, "lineage", "gbt")["fold_balanced_accuracy"][2] = 0.5
        self.assertIssue(self.run_audit(rep), "lineage/gbt/fold_ba")


class SplitSets(Case):
    def test_invalid_split_sets_record_issues_in_both_schemas(self):
        cases = [
            ("empty", []),
            ("duplicate random", [0, 0]),
            ("missing random", [1]),
            ("duplicate lineage", [1, 1]),
            ("missing lineage", [0]),
            ("extra random", [0, 0, 1]),
            ("extra lineage", [0, 1, 1]),
            ("unknown scheme", [0, 1]),
        ]
        for revision in (6, 7):
            for case, positions in cases:
                with self.subTest(revision=revision, case=case):
                    rep = fx.make_report(revision=revision)
                    winner = "rf" if revision == 7 else "gbt"
                    index = fx.manifest(rep, winner, schema=rep["schema"], evaluation_revision=revision)
                    log = run_log(rep) if revision == 7 else [fx.old_done_line(winner, ba(rep, "random", winner))]
                    rep["splits"] = [copy.deepcopy(rep["splits"][i]) for i in positions]
                    if case == "unknown scheme":
                        rep["splits"][1]["scheme"] = "other"
                    rec = self.run_audit(rep, index, backup=None, log_lines=log[-1:], revision=revision)
                    self.assertIssue(rec, "split schemes")

    def test_random_must_be_available_in_both_schemas(self):
        for revision in (6, 7):
            with self.subTest(revision=revision):
                rep = fx.make_report(revision=revision)
                index = fx.manifest(rep, None, schema=rep["schema"], evaluation_revision=revision)
                rep["splits"][0] = {"scheme": "random", "folds": 5, "n_groups": 16,
                                    "learners": [], "unavailable_reason": "fixture: unavailable"}
                rec = self.run_audit(rep, index, backup=None, log_lines=[], revision=revision)
                self.assertIssue(rec, "split schemes")

    def test_lineage_unavailable_remains_valid_in_both_schemas(self):
        for revision in (6, 7):
            with self.subTest(revision=revision):
                rep = fx.make_report(revision=revision, lineage_available=False)
                index = fx.manifest(rep, "gbt", schema=rep["schema"], evaluation_revision=revision)
                line = fx.done_line if revision == 7 else fx.old_done_line
                rec = self.run_audit(rep, index, backup=None,
                                     log_lines=[line("gbt", ba(rep, "random", "gbt"))], revision=revision)
                self.assertClean(rec)
                self.assertEqual(rec["unavailable"], [{"scheme": "lineage", "reason": "fixture: lineage split unavailable"}])

    def test_split_order_does_not_matter_in_both_schemas(self):
        for revision in (6, 7):
            with self.subTest(revision=revision):
                rep = fx.make_report(revision=revision)
                rep["splits"].reverse()
                winner = "rf" if revision == 7 else "gbt"
                index = fx.manifest(rep, winner, schema=rep["schema"], evaluation_revision=revision)
                line = fx.done_line if revision == 7 else fx.old_done_line
                rec = self.run_audit(rep, index, backup=None,
                                     log_lines=[line(winner, ba(rep, "random", winner))], revision=revision)
                self.assertClean(rec)


class TestFold(Case):
    def test_two_entries_swapped(self):
        rep = fx.make_report()
        folds = rep["splits"][0]["test_fold"]
        a = folds.index(0)
        b = folds.index(1)
        folds[a], folds[b] = folds[b], folds[a]
        rec = self.run_audit(rep)
        self.assertIssue(rec, "random/test_fold", '"positions": 2')

    def test_missing(self):
        rep = fx.make_report()
        del rep["splits"][1]["test_fold"]
        self.assertIssue(self.run_audit(rep), "lineage/test_fold", '"saved_len": null')

    def test_saved_in_a_schema1_report(self):
        rep = fx.make_report(revision=6)
        rep["splits"][0]["test_fold"] = list(fx.FOLDS["random"])
        rec = self.run_audit(rep, fx.manifest(rep, "gbt", schema=1, evaluation_revision=6), backup=None,
                             log_lines=[fx.old_done_line("gbt", ba(rep, "random", "gbt"))], revision=6)
        self.assertIssue(rec, "schema-1 split saves test_fold")

    def test_saved_on_an_unavailable_split(self):
        rep = fx.make_report(lineage_available=False)
        rep["splits"][1]["test_fold"] = list(fx.FOLDS["lineage"])
        rec = self.run_audit(rep, fx.manifest(rep, "gbt"), log_lines=run_log(rep, "gbt"))
        self.assertIssue(rec, "unavailable split saves test_fold")

    def test_folds_are_rebuilt_from_the_recorded_seed(self):
        rep = fx.make_report()
        rep["eval_config"]["seed"] = 43
        self.assertIssue(self.run_audit(rep), "random/test_fold", "lineage/test_fold")


class KsetDuplicates(Case):
    def setUp(self):
        # Audit the fixture's revision-7 reports as an isolating revision.
        saved = audit.DUPLICATE_ISOLATION_REVISION
        audit.DUPLICATE_ISOLATION_REVISION = 7
        self.addCleanup(setattr, audit, "DUPLICATE_ISOLATION_REVISION", saved)

    def test_an_isolated_pair_audits_clean(self):
        rep, digests = fx.isolated_report()
        self.assertClean(self.run_audit(rep, digests=digests))

    def test_a_listed_pair_left_across_two_folds(self):
        rep, digests = fx.isolated_report()
        b = rep["kset_duplicates"][0][1]
        rep["splits"][0]["test_fold"][b] = fx.FOLDS["random"][b]
        fx.rescore(rep)
        self.assertIssue(self.run_audit(rep, digests=digests), "random/test_fold", '"positions": 1')

    def test_a_tie_sent_to_the_higher_index_members_fold(self):
        rep, digests = fx.isolated_report()
        a, b = rep["kset_duplicates"][0]
        folds = rep["splits"][0]["test_fold"]
        folds[a] = folds[b] = fx.FOLDS["random"][b]
        fx.rescore(rep)
        self.assertIssue(self.run_audit(rep, digests=digests), "random/test_fold", '"positions": 2')

    def test_a_listed_pair_whose_digests_differ(self):
        rep, _ = fx.isolated_report()
        self.assertIssue(self.run_audit(rep), "kset_duplicates differ")

    def test_an_unlisted_pair_with_one_digest(self):
        _, digests = fx.isolated_report()
        self.assertIssue(self.run_audit(fx.make_report(), digests=digests), "kset_duplicates differ")

    def test_reports_before_the_isolating_revision_list_no_pair(self):
        audit.DUPLICATE_ISOLATION_REVISION = 8
        rep, digests = fx.isolated_report()
        self.assertClean(self.run_audit(fx.make_report(), digests=digests))
        self.assertIssue(self.run_audit(rep, digests=digests), "kset_duplicates differ")


class SchemaAndConfig(Case):
    def test_revision7_in_schema1(self):
        rep = fx.make_report()
        rep["schema"] = 1
        self.assertIssue(self.run_audit(rep), "schema does not match the revision")

    def test_revision6_in_schema2_outside_a_reshape(self):
        rep = fx.make_report()
        rep["evaluation"]["revision"] = 6
        self.assertIssue(self.run_audit(rep, revision=6), "schema does not match the revision")

    def test_a_recorded_forest_seed(self):
        rep = fx.make_report()
        rep["eval_config"]["forest"]["seed"] = 17
        self.assertIssue(self.run_audit(rep), "eval_config.forest fields")

    def test_a_recorded_boosting_seed(self):
        rep = fx.make_report()
        rep["eval_config"]["gbt"]["seed"] = 23
        self.assertIssue(self.run_audit(rep), "eval_config.gbt fields")

    def test_an_unknown_eval_config_field(self):
        rep = fx.make_report()
        rep["eval_config"]["rf_seed"] = 59
        self.assertIssue(self.run_audit(rep), '"what": "eval_config fields"')

    def test_a_revision7_report_without_eval_config(self):
        rep = fx.make_report()
        del rep["eval_config"]
        self.assertIssue(self.run_audit(rep), "without eval_config")


class LogLines(Case):
    def test_the_partial_coverage_top_scorer_logged_as_best(self):
        rep = fx.make_report()
        self.assertIssue(self.run_audit(rep, log_lines=run_log(rep, "gbt")), "log best learner differs")

    def test_no_winner_logged_beside_a_full_coverage_learner(self):
        rep = fx.make_report()
        self.assertIssue(self.run_audit(rep, log_lines=run_log(rep, None)), "log best learner differs")

    def test_the_old_form_for_a_schema2_report(self):
        rep = fx.make_report()
        log = [fx.recompute_line(), fx.old_done_line("rf", ba(rep, "random", "rf"))]
        self.assertIssue(self.run_audit(rep, log_lines=log), "done line form does not match the schema")

    def test_no_done_line(self):
        rep = fx.make_report()
        self.assertIssue(self.run_audit(rep, log_lines=[fx.recompute_line()]), "file without a done log line")

    def test_no_recompute_line_over_a_backup_report(self):
        rep = fx.make_report()
        self.assertIssue(self.run_audit(rep, log_lines=run_log(rep)[1:]), "revision-7 file without a recompute log line")

    def test_a_recompute_line_without_a_backup_report(self):
        rep = fx.make_report()
        self.assertIssue(self.run_audit(rep, backup=None), "recompute logged for a dataset without a backup report")

    def test_the_log_parser_reads_each_form_reason_and_species(self):
        with tempfile.TemporaryDirectory() as tmp:
            log = os.path.join(tmp, "log")
            with open(log, "w") as f:
                f.write("\n".join([
                    fx.recompute_line("evaluation audit revision changed", "a__x"),
                    fx.recompute_line("report has no input fingerprint", "b__y"),
                    fx.recompute_line("inputs changed since the report was written", "c__z"),
                    fx.old_done_line("catalog", 0.948, "a__x"),
                    fx.done_line("rf", 0.8124, "b__y"),
                    fx.done_line(None, None, "c__z"),
                    "fixture_species: 3 dataset reports written",
                ]) + "\n")
            audit.LOG = log
            audit.species_finished.clear()
            done, recomputed = audit.done_datasets()
        self.assertEqual({k: (v["form"], v["learner"], v["ba"]) for k, v in done.items()},
                         {"a__x": ("any row", "catalog", "0.948"), "b__y": ("full coverage", "rf", "0.812"),
                          "c__z": ("full coverage", None, None)})
        self.assertEqual(recomputed, {"a__x": "evaluation audit revision changed", "b__y": "report has no input fingerprint",
                                      "c__z": "inputs changed since the report was written"})
        self.assertEqual(audit.species_finished, {"fixture_species"})


class Provenance(Case):
    def test_a_file_older_than_launch(self):
        self.assertIssue(self.run_audit(fx.make_report(), launch="2999-01-01T00:00:00Z"), "revision-7 file older than launch")

    def test_a_file_identical_to_its_backup(self):
        rep = fx.make_report()
        self.assertIssue(self.run_audit(rep, backup=copy.deepcopy(rep)), "revision-7 file identical to its backup")

    def test_an_earlier_revision_written_after_launch_and_logged(self):
        rep = fx.make_report(revision=6)
        log = [fx.old_done_line("gbt", ba(rep, "random", "gbt")), "fixture_species: 1 dataset reports written"]
        rec = self.run_audit(rep, fx.manifest(rep, "gbt", schema=1, evaluation_revision=6), backup=None, log_lines=log)
        self.assertIssue(rec, "non-revision-7 file written after launch", "log reports a recompute for a non-revision-7 file",
                         "species finished in this run but file is not revision 7")


class Forest(Case):
    def test_forest_scores_unchanged_from_the_backup(self):
        rep = fx.make_report()
        rec = self.run_audit(rep, backup=fx.make_report(revision=6, forest_variant=0))
        self.assertIssue(rec, '"where": "random/rf", "what": "forest scores unchanged from the backup"',
                         '"where": "lineage/rf", "what": "forest scores unchanged from the backup"')

    def test_a_feature_outside_the_dataset_rows(self):
        rep = fx.make_report()
        rep["models"]["rf_importance"][0]["feature"] = fx.FEATURES
        self.assertIssue(self.run_audit(rep), "rf_importance feature is not a dataset row")

    def test_a_top_feature_without_a_card(self):
        rep = fx.make_report()
        rep["features"].pop(1)
        self.assertIssue(self.run_audit(rep), "rf_importance feature without a feature card")

    def test_a_card_whose_counts_differ_from_its_bits(self):
        rep = fx.make_report()
        rep["features"][0]["n_resistant"] += 1
        self.assertIssue(self.run_audit(rep), "rf feature card counts differ from its bits")

    def test_importance_out_of_order(self):
        rep = fx.make_report()
        imp = rep["models"]["rf_importance"]
        imp[1], imp[2] = imp[2], imp[1]
        self.assertIssue(self.run_audit(rep), "rf_importance not in weight then feature order")

    def test_a_zero_weight(self):
        rep = fx.make_report()
        rep["models"]["rf_importance"][2]["weight"] = 0.0
        self.assertIssue(self.run_audit(rep), "rf_importance weight not positive")


class Manifest(Case):
    def test_the_partial_coverage_top_scorer_named_best(self):
        rep = fx.make_report()
        self.assertIssue(self.run_audit(rep, fx.manifest(rep, "gbt")), "manifest best learner differs")

    def test_no_best_learner_beside_a_full_coverage_learner(self):
        rep = fx.make_report()
        self.assertIssue(self.run_audit(rep, fx.manifest(rep, None)), "manifest best learner differs")

    def test_zero_rates_for_a_row_without_calls(self):
        rep = fx.make_report()
        empty_gbt(rep)
        index = fx.manifest(rep, "rf")
        index["datasets"][0]["learners"][3]["random"].update(balanced_accuracy=0.0, sensitivity=0.0, specificity=0.0, vme=0.0, me=0.0)
        self.assertIssue(self.run_audit(rep, index), '"where": "manifest/random/gbt", "what": "manifest rates differ"')

    def test_null_rates_for_a_row_with_calls(self):
        rep = fx.make_report()
        index = fx.manifest(rep, "rf")
        index["datasets"][0]["learners"][2]["lineage"].update(balanced_accuracy=None)
        self.assertIssue(self.run_audit(rep, index), "manifest/lineage/rf")

    def test_rates_without_no_call(self):
        rep = fx.make_report()
        index = fx.manifest(rep, "rf")
        del index["datasets"][0]["learners"][4]["random"]["no_call"]
        self.assertIssue(self.run_audit(rep, index), "manifest/random/lasso")

    def test_a_schema2_report_in_a_schema1_manifest(self):
        rep = fx.make_report()
        self.assertIssue(self.run_audit(rep, fx.manifest(rep, "gbt", schema=1)), "schema-2 report in a schema-1 manifest")

    def test_a_manifest_built_from_another_copy_is_not_compared(self):
        rep = fx.make_report()
        rec = self.run_audit(rep, fx.manifest(rep, "gbt"), published=False)
        self.assertFalse(rec["manifest_compared"])
        self.assertNoIssue(rec, "manifest")

    def test_panel_rows_name_the_best_learner(self):
        rep = fx.make_report()
        panel = catalog_panel(rep)
        rep["panels"] = [panel]
        index = fx.manifest(rep, "rf")
        index["panels"] = [panel_row(rep, panel, "rf")]
        self.assertClean(self.run_audit(rep, index))
        index["panels"] = [panel_row(rep, panel, "gbt")]
        self.assertIssue(self.run_audit(rep, index), "panel rows do not follow the best learner")

    def test_a_dataset_without_a_best_learner_has_no_panel_rows(self):
        rep = fx.make_report()
        panel = catalog_panel(rep)
        rep["panels"] = [panel]
        for name in ("scm", "cart", "rf"):
            fx.row(rep, "random", name)["no_call"] = {"resistant": 1, "susceptible": 0}
        index = fx.manifest(rep, None)
        rec = self.run_audit(rep, index, log_lines=run_log(rep, None))
        self.assertNoIssue(rec, "manifest best learner differs")
        self.assertNoIssue(rec, "panel rows")
        self.assertNoIssue(rec, "log best learner differs")
        index["panels"] = [panel_row(rep, panel, "rf")]
        self.assertIssue(self.run_audit(rep, index, log_lines=run_log(rep, None)), "panel rows do not follow the best learner")


class Probes(Case):
    def probe(self, report, no_call, best, schema=2, evaluation_revision=7):
        """The report with a probe panel (40 R = 33 TP + 4 FN + 3 uncalled) and its manifest."""
        panel = {"id": "gc_gyra", "name": "fixture", "mirrors": "fixture", "targets": [],
                 "metrics": dict(audit.from_counts(33, 3, 35, 4), auroc=None, auprc=None),
                 "detected": [["mutation detected", 33, 3]]}
        if no_call:
            panel["no_call"] = {"resistant": 3, "susceptible": 2}
        report["panels"].append(panel)
        index = fx.manifest(report, best, schema, evaluation_revision)
        index["panels"] = [panel_row(report, panel, best)]
        return report, index

    def test_a_revision7_probe_without_no_call(self):
        rep, index = self.probe(fx.make_report(), False, "rf")
        self.assertIssue(self.run_audit(rep, index), "revision-7 probe panel without no_call")

    def test_a_revision7_probe_with_no_call(self):
        self.assertClean(self.run_audit(*self.probe(fx.make_report(), True, "rf")))

    def test_a_revision5_probe_without_no_call(self):
        rep, index = self.probe(fx.make_report(revision=5), False, "gbt", schema=1, evaluation_revision=5)
        rec = self.run_audit(rep, index, backup=None, log_lines=[fx.old_done_line("gbt", ba(rep, "random", "gbt"))], revision=5)
        self.assertNoIssue(rec, "probe panel without no_call")
        self.assertIssue(rec, "probe R partition")


class Rules(unittest.TestCase):
    """Winner rules, rates, and null handling, called directly."""

    def test_no_learner_that_calls_every_isolate_means_no_winner(self):
        rep = fx.make_report()
        for name in ("scm", "cart", "rf"):
            fx.row(rep, "lineage", name)["no_call"] = {"resistant": 0, "susceptible": 1}
        self.assertIsNone(audit.full_coverage_best(rep))

    def test_an_unavailable_lineage_split_does_not_require_a_learner_row(self):
        rep = fx.make_report()
        rep["splits"][1]["learners"] = []
        rep["splits"][1]["unavailable_reason"] = "Too few lineages"
        self.assertEqual(audit.full_coverage_best(rep)[0], "gbt")

    def test_a_missing_lineage_row_excludes_its_learner(self):
        rep = fx.make_report()
        rep["splits"][1]["learners"] = [r for r in rep["splits"][1]["learners"] if r["learner"] != "gbt"]
        self.assertEqual(audit.full_coverage_best(rep)[0], "rf")

    def test_a_tie_goes_to_the_later_learner(self):
        rep = fx.make_report()
        fx.row(rep, "random", "cart")["metrics"]["balanced_accuracy"] = ba(rep, "random", "rf")
        self.assertEqual(audit.full_coverage_best(rep)[0], "rf")

    def test_the_schema2_rule_leaves_out_the_catalog(self):
        rep = fx.make_report()
        fx.row(rep, "random", "catalog")["metrics"]["balanced_accuracy"] = 0.999
        self.assertEqual(audit.full_coverage_best(rep)[0], "rf")

    def test_schema1_rules_differ_on_the_catalog(self):
        rep = fx.make_report(revision=6)
        fx.row(rep, "random", "catalog")["metrics"]["balanced_accuracy"] = 0.999
        self.assertEqual(audit.any_row_best(rep)[0], "catalog")
        self.assertEqual(audit.max_random_best(rep)[0], "gbt")

    def test_the_schema1_manifest_rule_reads_null_metrics_as_zero(self):
        rep = fx.make_report()
        empty_gbt(rep)
        self.assertEqual(audit.max_random_best(rep), ("rf", ba(rep, "random", "rf")))

    def test_metrics_leave_out_null_scores_and_are_null_without_calls(self):
        self.assertIsNone(audit.metrics([True, False], [None, None], True))
        m = audit.metrics([True, False, True], [0.9, 0.1, None], True)
        self.assertEqual((m["n"], m["tp"], m["tn"], m["fn"], m["auroc"], m["auprc"]), (2, 1, 1, 0, 1.0, 1.0))

    def test_null_fold_tests_follow_the_schema(self):
        self.assertTrue(audit.is_null_fold(1, [0.5, 0.5]))
        self.assertFalse(audit.is_null_fold(2, [0.5, 0.5]))
        self.assertTrue(audit.is_null_fold(2, [None, None]))
        self.assertFalse(audit.is_null_fold(2, [None, 0.2]))
        self.assertFalse(audit.is_null_fold(2, []))

    def test_manifest_rates_keep_null(self):
        row = {"metrics": None, "no_call": {"resistant": 2, "susceptible": 1}}
        self.assertEqual(audit.lite(row, 2), {"balanced_accuracy": None, "sensitivity": None, "specificity": None,
                                              "vme": None, "me": None, "auroc": None, "no_call": row["no_call"]})
        self.assertFalse(audit.lite_differs(audit.lite(row, 2), fx.lite(row)))
        self.assertTrue(audit.lite_differs(dict(audit.lite(row, 2), vme=0.0), audit.lite(row, 2)))



class Curation(CurationCases, unittest.TestCase):
    count_key = "files"

    def run_scope(self, pattern="*"):
        out = Path(self.paths["root"], "scope-audit.json")
        result = subprocess.run(
            [sys.executable, "-I", os.path.join(fx.LANE, "audit.py"), str(out), pattern,
             "--revision", "7", "--root", self.paths["root"], "--backup", self.paths["backup"],
             "--log", str(self.log), "--launch", LAUNCH, "--foldgen", fx.FOLDGEN],
            capture_output=True, text=True)
        payload = json.loads(out.read_text()) if out.exists() else {}
        return result, {**payload.get("summary", {}), "records": payload.get("records", [])}

    def test_explicit_species_pattern_requires_results_without_run_inputs(self):
        rid = self.add_curated("unrun_species")
        self.assert_missing(rid, pattern="unrun_*")

    def test_species_pattern_limits_curation_and_file_checks_together(self):
        self.add_curated("another_species")
        other = self.data / "species/another_species"
        other.mkdir()
        (other / "genomes.txt").write_text("fixture-genome\n")
        (other / "counted.json").write_text("{}")
        _, summary = self.assert_passed(pattern=fx.SLUG)
        self.assertEqual(summary["not_run"], [])
        self.assertEqual(summary["missing_curated_results"], [])


if __name__ == "__main__":
    unittest.main(verbosity=2)
