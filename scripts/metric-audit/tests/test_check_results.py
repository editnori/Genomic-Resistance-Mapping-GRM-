"""check-results.py on synthetic reports: each invalid report must fail on
its own check, named in the assertion message, and each valid one pass.

usage: python3 -I tests/test_check_results.py
"""
import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixtures as fx  # noqa: E402
from curation_cases import CurationCases

spec = importlib.util.spec_from_file_location("check_results", fx.CHECKER)
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


def ba(report, scheme, learner):
    return fx.row(report, scheme, learner)["metrics"]["balanced_accuracy"]


class Case(unittest.TestCase):
    def check(self, report, log_lines=None, args=()):
        with tempfile.TemporaryDirectory() as tmp:
            root = fx.checker_root(tmp, report)
            extra = list(args)
            if log_lines is not None:
                log = os.path.join(tmp, "pipeline.log")
                Path(log).write_text("".join(line + "\n" for line in log_lines))
                extra += ["--log", log]
            result = fx.run_checker(root, *extra)
            receipt = None
            if result.returncode == 0:
                receipt = json.loads(Path(root, "logs", "result-validation-receipt.json").read_text())
            return result, receipt

    def passes(self, report, log_lines=None, args=()):
        result, receipt = self.check(report, log_lines, args)
        self.assertEqual(result.returncode, 0, result.stderr[-1500:])
        return receipt["reports"][0]

    def fails(self, report, reason, log_lines=None, args=()):
        result, _ = self.check(report, log_lines, args)
        self.assertNotEqual(result.returncode, 0, "accepted a report that should fail on: " + reason)
        last = result.stderr.strip().splitlines()[-1]
        self.assertTrue(last.startswith("AssertionError"), result.stderr[-1500:])
        self.assertIn(reason, last)


def rf_line(report):
    return fx.done_line("rf", ba(report, "random", "rf"))


def revision8_report():
    """A schema-2 report with the mandatory clean revision-8 lasso evidence."""
    rep = fx.make_report(revision=8)

    def status():
        return {"irls_steps": 3, "sweep_cap_hits": 0, "converged": True}

    def fit():
        return {"path_points": 30, "path_unconverged": 0, "selected_point": 2,
                "selected_convergence": status(), "refit_convergence": status()}

    for split in rep["splits"]:
        if split.get("unavailable_reason") is None:
            fx.row(rep, split["scheme"], "lasso")["lasso_fits"] = [fit() for _ in range(split["folds"])]
    rep["models"]["lasso"] = {"intercept": -1.0, "coef": [[4, 0.5]], "lambda": 0.01}
    rep["models"]["lasso_fit"] = fit()
    return rep


def empty_lasso_fit():
    """The record for a deliberate empty path: no penalty admits a feature."""
    return {"path_points": 0, "path_unconverged": 0, "selected_point": None,
            "selected_convergence": None, "refit_convergence": None}


class ValidReports(Case):
    def test_revision7_report_with_null_folds_passes_and_names_rf(self):
        rep = fx.make_report()
        entry = self.passes(rep, [rf_line(rep)])
        self.assertEqual((entry["schema"], entry["rows_with_no_call"], entry["best_learner"], entry["log_checked"]), (2, 2, "rf", True))

    def test_revision6_schema1_report_passes_with_its_old_done_line(self):
        rep = fx.make_report(revision=6)
        # The schema-1 log rule takes the best random-split row, catalog included.
        entry = self.passes(rep, [fx.old_done_line("gbt", ba(rep, "random", "gbt"))])
        self.assertEqual((entry["schema"], entry["rows_with_no_call"], entry["best_learner"]), (1, 0, "gbt"))

    def test_a_row_without_any_call_passes_with_null_metrics(self):
        rep = fx.make_report()
        gbt = fx.row(rep, "random", "gbt")
        gbt["scores"] = [None] * fx.N
        gbt["fold_complexity"] = [0.0] * 5
        fx.rescore(rep)
        self.assertIsNone(gbt["metrics"])
        entry = self.passes(rep, [rf_line(rep)])
        self.assertEqual(entry["rows_with_no_call"], 3)

    def test_an_unavailable_lineage_split_leaves_random_rows_to_decide(self):
        rep = fx.make_report(lineage_available=False)
        # gbt's null fold was on the lineage split, so gbt now calls every isolate that was evaluated.
        entry = self.passes(rep, [fx.done_line("gbt", ba(rep, "random", "gbt"))])
        self.assertEqual(entry["best_learner"], "gbt")

    def test_a_revision6_report_needs_no_done_line(self):
        self.passes(fx.make_report(revision=6), [fx.done_line(None, None, "other__drug")])

    def test_the_last_done_line_for_a_dataset_counts(self):
        rep = fx.make_report()
        self.passes(rep, [fx.done_line("gbt", ba(rep, "random", "gbt")), rf_line(rep)])

    def test_reshaped_accepts_schema2_reports_of_revision6(self):
        rep = fx.make_report()
        rep["evaluation"]["revision"] = 6
        self.passes(rep, args=["--reshaped"])
        self.fails(rep, "schema does not match the revision")


class NullScores(Case):
    def test_a_null_score_in_a_schema1_report_fails(self):
        rep = fx.make_report(revision=6)
        fx.row(rep, "random", "gbt")["scores"][0] = None
        self.fails(rep, "null score or no_call on a row that calls every isolate")

    def test_a_null_catalog_score_fails(self):
        rep = fx.make_report()
        fx.row(rep, "random", "catalog")["scores"][3] = None
        self.fails(rep, "null score or no_call on a row that calls every isolate")

    def test_a_null_rf_fold_fails_even_with_matching_no_call(self):
        rep = fx.make_report()
        rf = fx.row(rep, "random", "rf")
        rf["scores"] = [None if fx.FOLDS["random"][i] == 0 else s for i, s in enumerate(rf["scores"])]
        fx.rescore(rep)
        self.fails(rep, "only gbt and lasso fits can be null")

    def test_a_fold_called_in_part_fails(self):
        rep = fx.make_report()
        gbt = fx.row(rep, "lineage", "gbt")
        gbt["scores"][fx.FOLDS["lineage"].index(1)] = None
        fx.rescore(rep)
        self.fails(rep, "a fold is called in part")


class NoCallCounts(Case):
    def test_no_call_one_short_fails(self):
        rep = fx.make_report()
        fx.row(rep, "lineage", "gbt")["no_call"]["resistant"] -= 1
        self.fails(rep, "no_call differs from the null scores")

    def test_no_call_missing_beside_null_scores_fails(self):
        rep = fx.make_report()
        del fx.row(rep, "lineage", "gbt")["no_call"]
        self.fails(rep, "no_call differs from the null scores")

    def test_no_call_missing_on_a_fully_called_schema2_row_fails(self):
        rep = fx.make_report()
        del fx.row(rep, "random", "rf")["no_call"]
        self.fails(rep, "no_call differs from the null scores")

    def test_no_call_on_the_catalog_row_fails(self):
        rep = fx.make_report()
        fx.row(rep, "random", "catalog")["no_call"] = {"resistant": 0, "susceptible": 0}
        self.fails(rep, "null score or no_call on a row that calls every isolate")

    def test_no_call_on_a_schema1_row_fails(self):
        rep = fx.make_report(revision=6)
        fx.row(rep, "random", "rf")["no_call"] = {"resistant": 0, "susceptible": 0}
        self.fails(rep, "null score or no_call on a row that calls every isolate")


class CalledOnlyMetrics(Case):
    def test_n_counting_uncalled_isolates_fails(self):
        rep = fx.make_report()
        fx.row(rep, "lineage", "gbt")["metrics"]["n"] = fx.N
        self.fails(rep, "confusion of the called isolates")

    def test_uncalled_resistant_isolates_counted_as_misses_fails(self):
        rep = fx.make_report()
        gbt = fx.row(rep, "lineage", "gbt")
        gbt["metrics"]["fn"] += gbt["no_call"]["resistant"]
        gbt["metrics"]["n"] += gbt["no_call"]["resistant"]
        self.fails(rep, "confusion of the called isolates")

    def test_balanced_accuracy_over_every_isolate_fails(self):
        rep = fx.make_report()
        gbt = fx.row(rep, "lineage", "gbt")
        half = [0.5 if s is None else s for s in gbt["scores"]]
        gbt["metrics"]["balanced_accuracy"] = checker.balanced_accuracy([s["resistant"] for s in rep["samples"]], half)
        self.fails(rep, "balanced accuracy")


class EmptyCallSets(Case):
    def empty_gbt(self):
        rep = fx.make_report()
        gbt = fx.row(rep, "random", "gbt")
        gbt["scores"] = [None] * fx.N
        gbt["fold_complexity"] = [0.0] * 5
        fx.rescore(rep)
        return rep, gbt

    def test_zero_metrics_on_a_row_without_calls_fail(self):
        rep, gbt = self.empty_gbt()
        gbt["metrics"] = fx.row(rep, "random", "rf")["metrics"] | {"n": 0, "tp": 0, "fp": 0, "tn": 0, "fn": 0}
        self.fails(rep, "a row without calls needs null metrics")

    def test_roc_points_on_a_row_without_calls_fail(self):
        rep, gbt = self.empty_gbt()
        gbt["roc"] = [[0.0, 0.0], [1.0, 1.0]]
        self.fails(rep, "a row without calls needs null metrics")

    def test_fold_balanced_accuracy_on_a_row_without_calls_fails(self):
        rep, gbt = self.empty_gbt()
        gbt["fold_balanced_accuracy"] = [0.5] * 5
        self.fails(rep, "a row without calls needs null metrics")

    def test_null_metrics_on_a_row_with_calls_fail(self):
        rep = fx.make_report()
        fx.row(rep, "lineage", "gbt")["metrics"] = None
        self.fails(rep, "null metrics on a row with calls")

    def test_a_null_fold_with_a_balanced_accuracy_fails(self):
        rep = fx.make_report()
        fx.row(rep, "lineage", "gbt")["fold_balanced_accuracy"][2] = 0.5
        self.fails(rep, "a fold without calls needs null balanced accuracy")


class SplitSets(Case):
    def test_invalid_split_sets_fail_in_both_schemas(self):
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
                    rep["splits"] = [rep["splits"][i] for i in positions]
                    if case == "unknown scheme":
                        rep["splits"][1]["scheme"] = "other"
                    self.fails(rep, "split schemes")

    def test_random_must_be_available_in_both_schemas(self):
        for revision in (6, 7):
            with self.subTest(revision=revision):
                rep = fx.make_report(revision=revision)
                rep["splits"][0] = {"scheme": "random", "folds": 5, "n_groups": 16,
                                    "learners": [], "unavailable_reason": "fixture: unavailable"}
                self.fails(rep, "unavailable split")

    def test_lineage_unavailable_remains_valid_in_both_schemas(self):
        for revision in (6, 7):
            with self.subTest(revision=revision):
                rep = fx.make_report(revision=revision, lineage_available=False)
                entry = self.passes(rep)
                self.assertEqual(entry["best_learner"], "gbt")

    def test_split_order_does_not_matter_in_both_schemas(self):
        for revision in (6, 7):
            with self.subTest(revision=revision):
                rep = fx.make_report(revision=revision)
                rep["splits"].reverse()
                self.passes(rep)


class TestFold(Case):
    def test_a_missing_test_fold_fails(self):
        rep = fx.make_report()
        del rep["splits"][0]["test_fold"]
        self.fails(rep, "'test_fold'")

    def test_a_short_test_fold_fails(self):
        rep = fx.make_report()
        rep["splits"][0]["test_fold"].pop()
        self.fails(rep, "'test_fold'")

    def test_a_fold_id_past_the_last_fold_fails(self):
        rep = fx.make_report()
        rep["splits"][0]["test_fold"][0] = 5
        self.fails(rep, "'test_fold'")

    def test_two_isolates_swapped_between_folds_fail_on_fold_balanced_accuracy(self):
        rep = fx.make_report()
        folds = rep["splits"][0]["test_fold"]
        # Two resistant isolates in different folds, the first called wrongly by rf and the second rightly.
        a = next(i for i in range(fx.N) if fx.TRUTH[i] and fx.wrong("random", "rf", i))
        b = next(i for i in range(fx.N) if folds[i] != folds[a] and fx.TRUTH[i] and not fx.wrong("random", "rf", i))
        folds[a], folds[b] = folds[b], folds[a]
        self.fails(rep, "fold balanced accuracy")

    def test_a_lineage_split_across_two_folds_fails(self):
        rep = fx.make_report()
        folds = rep["splits"][1]["test_fold"]
        folds[0] = (folds[0] + 1) % 5
        self.fails(rep, "a lineage spans test folds")

    def test_a_test_fold_holding_one_class_fails(self):
        rep = fx.make_report()
        folds = rep["splits"][0]["test_fold"]
        for i in range(fx.N):
            if folds[i] == 4 and not fx.TRUTH[i]:
                folds[i] = 3
        self.fails(rep, "test fold lacks a class")

    def test_test_fold_in_a_schema1_report_fails(self):
        rep = fx.make_report(revision=6)
        rep["splits"][0]["test_fold"] = list(fx.FOLDS["random"])
        self.fails(rep, "test_fold in a schema-1 report")

    def test_test_fold_on_an_unavailable_split_fails(self):
        rep = fx.make_report(lineage_available=False)
        rep["splits"][1]["test_fold"] = list(fx.FOLDS["lineage"])
        self.fails(rep, "unavailable split")


class KsetDuplicates(Case):
    def test_an_isolated_pair_before_revision_10_fails_on_the_revision_gate(self):
        # A well-formed isolated pair, so only the revision gate refuses it.
        rep, _ = fx.isolated_report()
        self.fails(rep, "kset_duplicates in a report before revision 10")

    def test_a_pair_across_two_test_folds_fails(self):
        rep, _ = fx.isolated_report()
        b = rep["kset_duplicates"][0][1]
        rep["splits"][0]["test_fold"][b] = fx.FOLDS["random"][b]
        fx.rescore(rep)
        self.fails(rep, "a kset duplicate group spans test folds")

    def test_a_pair_across_two_lineages_fails(self):
        rep = fx.make_report()
        rep["kset_duplicates"] = [[0, 5]]
        self.fails(rep, "a kset duplicate group spans lineages")

    def test_malformed_groups_fail(self):
        for groups, reason in [([], "malformed kset_duplicates"), ([[3]], "malformed kset_duplicates"),
                               ([[2, 1]], "malformed kset_duplicates"), ([[0, 1], [1, 2]], "overlapping"),
                               ([[5, 6], [0, 1]], "out of order"), ([[78, 80]], "out of range")]:
            with self.subTest(groups=groups):
                rep = fx.make_report()
                rep["kset_duplicates"] = groups
                self.fails(rep, reason)


class SchemaAndRevision(Case):
    def test_revision7_in_schema1_fails(self):
        rep = fx.make_report()
        rep["schema"] = 1
        self.fails(rep, "schema does not match the revision")

    def test_revision6_in_schema2_fails_without_reshaped(self):
        rep = fx.make_report()
        rep["evaluation"]["revision"] = 6
        self.fails(rep, "schema does not match the revision")

    def test_reshaped_refuses_schema1(self):
        self.fails(fx.make_report(revision=6), "schema does not match the revision", args=["--reshaped"])

    def test_a_report_without_schema_fails(self):
        rep = fx.make_report()
        del rep["schema"]
        self.fails(rep, "schema does not match the revision")

    def test_revision8_report_with_clean_solver_status_passes(self):
        rep = revision8_report()
        entry = self.passes(rep, [rf_line(rep)])
        self.assertEqual((entry["revision"], entry["schema"], entry["best_learner"]), (8, 2, "rf"))

    def test_future_revision9_fails(self):
        rep = revision8_report()
        rep["evaluation"]["revision"] = 9
        self.fails(rep, "'revision'")


class LassoSolverStatus(Case):
    def test_missing_per_fold_lasso_status_fails(self):
        rep = revision8_report()
        del fx.row(rep, "random", "lasso")["lasso_fits"]
        self.fails(rep, "revision-8 lasso fit status is required")

    def test_missing_final_lasso_status_fails(self):
        rep = revision8_report()
        del rep["models"]["lasso_fit"]
        self.fails(rep, "revision-8 lasso fit status is required")

    def test_malformed_convergence_boolean_fails(self):
        rep = revision8_report()
        fx.row(rep, "lineage", "lasso")["lasso_fits"][1]["selected_convergence"]["converged"] = "true"
        self.fails(rep, "malformed lasso convergence status")

    def test_malformed_path_count_fails(self):
        rep = revision8_report()
        fx.row(rep, "random", "lasso")["lasso_fits"][0]["path_unconverged"] = "0"
        self.fails(rep, "malformed lasso fit status")

    def test_unconverged_path_point_fails(self):
        rep = revision8_report()
        fx.row(rep, "random", "lasso")["lasso_fits"][2]["path_unconverged"] = 1
        self.fails(rep, "unconverged two-class lasso fit")

    def test_unconverged_selected_point_fails(self):
        rep = revision8_report()
        fx.row(rep, "lineage", "lasso")["lasso_fits"][0]["selected_convergence"]["converged"] = False
        self.fails(rep, "unconverged two-class lasso fit")

    def test_unconverged_final_refit_fails(self):
        rep = revision8_report()
        rep["models"]["lasso_fit"]["refit_convergence"]["converged"] = False
        self.fails(rep, "unconverged two-class lasso fit")

    def test_selected_point_without_refit_status_fails(self):
        rep = revision8_report()
        rep["models"]["lasso_fit"]["refit_convergence"] = None
        self.fails(rep, "malformed lasso fit status")

    def test_empty_final_path_with_a_null_model_passes(self):
        rep = revision8_report()
        rep["models"]["lasso"] = None
        rep["models"]["lasso_fit"] = empty_lasso_fit()
        self.passes(rep, [rf_line(rep)])

    def test_saved_lasso_model_without_a_refit_record_fails(self):
        rep = revision8_report()
        rep["models"]["lasso_fit"] = empty_lasso_fit()
        self.fails(rep, "lasso model does not match its solver status")

    def test_a_wrong_number_of_fold_records_fails(self):
        for edit in ("drop", "add"):
            with self.subTest(edit=edit):
                rep = revision8_report()
                fits = fx.row(rep, "random", "lasso")["lasso_fits"]
                if edit == "drop":
                    fits.pop()
                else:
                    fits.append(dict(fits[0]))
                self.fails(rep, "malformed lasso fit status")

    def test_a_nonempty_path_without_a_selected_point_fails(self):
        # Random fold 4 is a null fit, so its calls cannot expose the record.
        rep = revision8_report()
        fx.row(rep, "random", "lasso")["lasso_fits"][4] = empty_lasso_fit() | {"path_points": 30}
        self.fails(rep, "malformed lasso fit status")


class LassoRecordsAgreeWithScores(Case):
    """Each revision-8 fold record against that fold's calls and weight count."""

    def test_an_empty_path_on_an_uncalled_fold_passes(self):
        rep = revision8_report()
        fx.row(rep, "random", "lasso")["lasso_fits"][4] = empty_lasso_fit()
        self.passes(rep, [rf_line(rep)])

    def test_an_empty_path_on_a_called_fold_fails(self):
        rep = revision8_report()
        fx.row(rep, "random", "lasso")["lasso_fits"][0] = empty_lasso_fit()
        self.fails(rep, "lasso fold without a selected point has calls")

    def test_a_weight_count_that_disagrees_with_the_calls_fails(self):
        for fold, weights in ((0, 0.0), (4, 12.0)):
            with self.subTest(fold=fold):
                rep = revision8_report()
                fx.row(rep, "random", "lasso")["fold_complexity"][fold] = weights
                self.fails(rep, "lasso fold calls do not match its weight count")

    def test_a_short_weight_count_list_fails(self):
        rep = revision8_report()
        fx.row(rep, "lineage", "lasso")["fold_complexity"].pop()
        self.fails(rep, "fold complexity length")

    def test_revision7_rows_need_no_agreement(self):
        rep = fx.make_report()
        fx.row(rep, "random", "lasso")["fold_complexity"][0] = 0.0
        self.passes(rep, [rf_line(rep)])


class LearnerRows(Case):
    def test_a_lineage_split_missing_the_gbt_row_fails(self):
        rep = fx.make_report()
        rep["splits"][1]["learners"] = [r for r in rep["splits"][1]["learners"] if r["learner"] != "gbt"]
        self.fails(rep, "learner rows")

    def test_a_split_missing_the_catalog_row_fails(self):
        rep = fx.make_report()
        rep["splits"][0]["learners"].pop()
        self.fails(rep, "learner rows")


class DoneLines(Case):
    def test_the_partial_coverage_top_scorer_fails_as_the_logged_winner(self):
        rep = fx.make_report()
        self.fails(rep, "done line names another best learner", [fx.done_line("gbt", ba(rep, "random", "gbt"))])
        self.fails(rep, "done line names another best learner", [fx.done_line("lasso", ba(rep, "random", "lasso"))])

    def test_no_winner_logged_beside_a_full_coverage_learner_fails(self):
        self.fails(fx.make_report(), "done line names another best learner", [fx.done_line(None, None)])

    def test_another_balanced_accuracy_fails(self):
        self.fails(fx.make_report(), "done line names another best learner", [fx.done_line("rf", 0.949)])

    def test_the_old_form_for_a_revision7_report_fails(self):
        rep = fx.make_report()
        self.fails(rep, "done line form does not match the report schema", [fx.old_done_line("rf", ba(rep, "random", "rf"))])

    def test_the_new_form_for_a_schema1_report_fails(self):
        rep = fx.make_report(revision=6)
        self.fails(rep, "done line form does not match the report schema", [fx.done_line("gbt", ba(rep, "random", "gbt"))])

    def test_an_unrecognized_done_line_fails(self):
        self.fails(fx.make_report(), "done line form does not match the report schema", ["    %s done in 3s (best: rf)" % fx.ID])

    def test_a_revision7_report_without_a_done_line_fails(self):
        self.fails(fx.make_report(), "revision-7 report without a done line", [fx.recompute_line()])


class Rules(unittest.TestCase):
    """The best-learner rules and the log parser, called directly."""

    def test_no_learner_that_calls_every_isolate_means_no_winner(self):
        rep = fx.make_report()
        for name in fx.LEARNERS:
            fx.row(rep, "random", name)["no_call"] = {"resistant": 1, "susceptible": 0}
        self.assertIsNone(checker.full_coverage_best(rep))

    def test_an_unavailable_lineage_split_does_not_require_a_learner_row(self):
        rep = fx.make_report()
        rep["splits"][1]["learners"] = []
        rep["splits"][1]["unavailable_reason"] = "Too few lineages"
        self.assertEqual(checker.full_coverage_best(rep)[0], "gbt")

    def test_a_missing_lineage_row_excludes_its_learner(self):
        rep = fx.make_report()
        rep["splits"][1]["learners"] = [r for r in rep["splits"][1]["learners"] if r["learner"] != "gbt"]
        self.assertEqual(checker.full_coverage_best(rep)[0], "rf")

    def test_a_missing_random_row_excludes_its_learner(self):
        rep = fx.make_report()
        rep["splits"][0]["learners"] = [r for r in rep["splits"][0]["learners"] if r["learner"] != "rf"]
        self.assertEqual(checker.full_coverage_best(rep)[0], "cart")

    def test_a_tie_goes_to_the_later_learner(self):
        rep = fx.make_report()
        fx.row(rep, "random", "scm")["metrics"]["balanced_accuracy"] = ba(rep, "random", "rf")
        self.assertEqual(checker.full_coverage_best(rep)[0], "rf")
        fx.row(rep, "random", "cart")["metrics"]["balanced_accuracy"] = 0.99
        fx.row(rep, "random", "rf")["metrics"]["balanced_accuracy"] = 0.99
        self.assertEqual(checker.full_coverage_best(rep)[0], "rf")

    def test_null_metrics_exclude_a_row_without_no_call(self):
        rep = fx.make_report()
        rf = fx.row(rep, "random", "rf")
        del rf["no_call"]
        rf["metrics"] = None
        self.assertEqual(checker.full_coverage_best(rep)[0], "cart")

    def test_the_schema2_rule_leaves_out_the_catalog(self):
        rep = fx.make_report()
        fx.row(rep, "random", "catalog")["metrics"]["balanced_accuracy"] = 0.999
        self.assertEqual(checker.full_coverage_best(rep)[0], "rf")

    def test_the_schema1_rule_counts_the_catalog(self):
        rep = fx.make_report(revision=6)
        fx.row(rep, "random", "catalog")["metrics"]["balanced_accuracy"] = 0.999
        self.assertEqual(checker.any_row_best(rep), ("catalog", 0.999))

    def test_done_lines_reads_each_form_and_keeps_the_last(self):
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp, "log")
            log.write_text("\n".join([
                "klebsiella_pneumoniae: 38 dataset reports written",
                fx.old_done_line("catalog", 0.948, "a__x"),
                fx.done_line("rf", 0.8124, "b__y"),
                fx.done_line(None, None, "c__z"),
                "    d__w done in 9s (something else)",
                fx.recompute_line(dataset="b__y"),
                fx.done_line("gbt", 0.7, "a__x"),
            ]) + "\n")
            self.assertEqual(checker.done_lines(log), {
                "a__x": ("full coverage", "gbt", "0.700"),
                "b__y": ("full coverage", "rf", "0.812"),
                "c__z": ("full coverage", None, None),
                "d__w": ("unrecognized", "something else", None),
            })



class Curation(CurationCases, unittest.TestCase):
    count_key = "reports"

    def test_explicit_receipt_path_keeps_the_data_root_read_only(self):
        out = Path(self.paths["root"], "external-receipt.json")
        result = fx.run_checker(str(self.data), "--log", str(self.log), "--out", str(out))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(json.loads(out.read_text())["reports"]), 1)
        self.assertFalse((self.data / "logs").exists())

    def reused_report_with_skip(self):
        rid = self.add_curated(name="class_minimum")
        self.log.write_text("    " + rid + ": too few labeled genomes after download, skipped\n")
        return rid

    def test_skip_log_alone_accepts_reused_reports_without_done_lines(self):
        rid = self.reused_report_with_skip()
        result = fx.run_checker(str(self.data), "--skip-log", str(self.log))
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = json.loads(result.stdout)
        self.assertEqual((summary["reports"], summary["log_checked"]), (1, 0))
        self.assertEqual(summary["class_minimum_skipped"], [{
            "file": str(self.result_path.with_name(rid + ".json")),
            "reason": "too few labeled genomes after download, skipped",
        }])
        receipt = json.loads((self.data / "logs/result-validation-receipt.json").read_text())
        self.assertIsNone(receipt["log"])
        self.assertEqual(receipt["skip_log"], str(self.log))
        self.assertFalse(receipt["reports"][0]["log_checked"])

    def test_reused_reports_with_a_missing_result_fail_without_skip_evidence(self):
        rid = self.reused_report_with_skip()
        result = fx.run_checker(str(self.data))
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("missing results for curated dataset", result.stderr)
        self.assertIn(rid + ".json", result.stderr)

    def test_skip_log_skip_is_revoked_by_a_later_done_or_recompute_line(self):
        rid = self.reused_report_with_skip()
        skipped = self.log.read_text()
        for later in [fx.done_line(None, None, rid), fx.recompute_line(dataset=rid)]:
            with self.subTest(later=later):
                self.log.write_text(skipped + later + "\n")
                result = fx.run_checker(str(self.data), "--skip-log", str(self.log))
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertIn("missing results for curated dataset", result.stderr)
                self.assertIn(rid + ".json", result.stderr)

    def test_skip_log_does_not_disable_explicit_done_line_validation(self):
        self.reused_report_with_skip()
        result = fx.run_checker(str(self.data), "--skip-log", str(self.log), "--log", str(self.log))
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("revision-7 report without a done line", result.stderr)

    def test_separate_skip_and_done_logs_record_and_validate_both_sources(self):
        rid = self.add_curated(name="class_minimum")
        skip_log = self.log.with_name("skip-only.log")
        skip_log.write_text("    " + rid + ": too few labeled genomes after download, skipped\n")
        result = fx.run_checker(str(self.data), "--log", str(self.log), "--skip-log", str(skip_log))
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = json.loads(result.stdout)
        self.assertEqual((summary["reports"], summary["log_checked"]), (1, 1))
        self.assertEqual(len(summary["class_minimum_skipped"]), 1)
        receipt = json.loads((self.data / "logs/result-validation-receipt.json").read_text())
        self.assertEqual(receipt["log"], str(self.log))
        self.assertEqual(receipt["skip_log"], str(skip_log))

    def run_scope(self):
        result = fx.run_checker(str(self.data), "--log", str(self.log))
        summary = json.loads(result.stdout) if result.returncode == 0 else {}
        if result.returncode == 0:
            receipt = json.loads((self.data / "logs/result-validation-receipt.json").read_text())
            for key in ("not_in_current_curation", "not_run", "class_minimum_skipped", "missing_curated_results"):
                self.assertEqual(receipt.get(key), summary.get(key))
        return result, summary


if __name__ == "__main__":
    unittest.main(verbosity=2)
