"""Run the same curation and run-scope regressions through both checker CLIs."""
import json
import tempfile
from pathlib import Path

import fixtures as fx


class CurationCases:
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.report = fx.make_report()
        lines = [fx.recompute_line(), fx.done_line(
            "rf", fx.row(self.report, "random", "rf")["metrics"]["balanced_accuracy"])]
        self.paths = fx.audit_root(tmp.name, self.report,
                                   backup=fx.make_report(revision=6, forest_variant=1),
                                   log_lines=lines)
        self.data = Path(self.paths["root"], "data")
        self.curated = self.data / "curated/datasets.json"
        self.result_path = Path(self.paths["report"])
        self.log = Path(self.paths["log"])

    def add_curated(self, slug=fx.SLUG, name="missing"):
        rid = slug + "__" + name
        rows = json.loads(self.curated.read_text())
        rows.append({"id": rid, "species_slug": slug})
        self.curated.write_text(json.dumps(rows))
        return rid

    def append_log(self, *lines):
        with self.log.open("a") as stream:
            stream.write("\n".join(lines) + "\n")

    def assert_passed(self, count=1, **kwargs):
        result, summary = self.run_scope(**kwargs)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(summary.get("errors"), summary)
        self.assertFalse(summary.get("issue_count"), summary)
        self.assertEqual(summary[self.count_key], count)
        return result, summary

    def assert_missing(self, rid, **kwargs):
        result, summary = self.run_scope(**kwargs)
        self.assertNotEqual(result.returncode, 0, result.stdout)
        evidence = result.stderr + json.dumps(summary)
        self.assertIn("missing results for curated dataset", evidence)
        self.assertIn(rid + ".json", evidence)

    def test_dropped_results_are_listed_without_reading_or_counting_them(self):
        invalid = self.result_path.with_name(fx.SLUG + "__dropped.json")
        invalid.write_text("not even JSON")
        stale = self.data / "species/retired_species/results/retired_species__old.json"
        fx.write_json(str(stale), self.report)
        result, summary = self.assert_passed()
        self.assertEqual(summary["not_in_current_curation"], [
            {"file": str(invalid), "reason": "not in current curation"},
            {"file": str(stale), "reason": "not in current curation"},
        ])
        self.assertIn(str(invalid) + ": not in current curation", result.stderr)
        self.assertIn(str(stale) + ": not in current curation", result.stderr)

    def test_missing_curated_result_fails_in_a_species_with_results(self):
        rid = self.add_curated()
        self.assert_missing(rid)

    def test_all_results_missing_fails_when_both_run_inputs_exist(self):
        self.result_path.unlink()
        (self.result_path.parent.parent / "genomes.txt").write_text("fixture-genome\n")
        self.assert_missing(fx.ID)

    def test_unrun_curated_species_are_listed_without_being_counted(self):
        for slug, marker in [("genomes_only", "genomes.txt"), ("counted_only", "counted.json")]:
            self.add_curated(slug)
            folder = self.data / "species" / slug
            folder.mkdir()
            (folder / marker).write_text("[]")
        self.add_curated("absent_species")
        result, summary = self.assert_passed()
        expected = ["absent_species", "counted_only", "genomes_only"]
        self.assertEqual(summary["not_run"], [{"species": s, "reason": "not run"} for s in expected])
        for slug in expected:
            self.assertIn(slug + ": not run", result.stderr)

    def test_exact_class_minimum_skip_excuses_only_the_missing_dataset(self):
        rid = self.add_curated()
        self.append_log("    " + rid + ": too few labeled genomes after download, skipped")
        _, summary = self.assert_passed()
        self.assertEqual(summary["class_minimum_skipped"], [{
            "file": str(self.result_path.with_name(rid + ".json")),
            "reason": "too few labeled genomes after download, skipped",
        }])

    def test_class_minimum_skip_does_not_hide_another_missing_dataset(self):
        skipped = self.add_curated(name="skipped")
        missing = self.add_curated()
        self.append_log("    " + skipped + ": too few labeled genomes after download, skipped")
        self.assert_missing(missing)

    def test_near_match_skip_lines_do_not_excuse_missing_results(self):
        rid = self.add_curated()
        original = self.log.read_text()
        for line in [
            "    " + rid + ": too few labeled genomes after download, skipped (maybe)",
            "warning: " + rid + ": too few labeled genomes after download, skipped",
            "    other__missing: too few labeled genomes after download, skipped",
        ]:
            with self.subTest(line=line):
                self.log.write_text(original + line + "\n")
                self.assert_missing(rid)

    def test_a_later_done_or_recompute_invalidates_a_class_minimum_skip(self):
        rid = self.add_curated()
        self.append_log("    " + rid + ": too few labeled genomes after download, skipped")
        skipped_log = self.log.read_text()
        for later in [fx.done_line(None, None, rid), fx.recompute_line(dataset=rid)]:
            with self.subTest(later=later):
                self.log.write_text(skipped_log + later + "\n")
                self.assert_missing(rid)

    def test_a_later_class_minimum_skip_excuses_a_removed_result(self):
        rid = self.add_curated()
        self.append_log(fx.done_line(None, None, rid),
                        "    " + rid + ": too few labeled genomes after download, skipped")
        _, summary = self.assert_passed()
        self.assertEqual(len(summary["class_minimum_skipped"]), 1)

    def test_a_skip_cannot_exempt_an_existing_invalid_curated_report(self):
        self.report["counts"]["samples"] += 1
        fx.write_json(str(self.result_path), self.report)
        self.append_log("    " + fx.ID + ": too few labeled genomes after download, skipped")
        result, summary = self.run_scope()
        self.assertTrue(result.returncode != 0 or summary.get("issue_count"), result.stdout)
        self.assertNotIn(": not in current curation", result.stderr)

    def test_a_missing_curation_file_fails_instead_of_auditing_everything(self):
        self.curated.unlink()
        result, _ = self.run_scope()
        self.assertNotEqual(result.returncode, 0, result.stdout)
        self.assertIn("datasets.json", result.stderr)

    def test_a_result_in_another_species_cannot_replace_the_curated_result(self):
        (self.result_path.parent.parent / "genomes.txt").write_text("fixture-genome\n")
        wrong = self.data / "species/wrong_species/results" / self.result_path.name
        wrong.parent.mkdir(parents=True)
        self.result_path.rename(wrong)
        self.assert_missing(fx.ID)

    def test_an_existing_report_must_match_its_curated_path(self):
        self.report["id"] = "other__dataset"
        fx.write_json(str(self.result_path), self.report)
        result, summary = self.run_scope()
        self.assertTrue(result.returncode != 0 or summary.get("issue_count"), result.stdout)
        self.assertIn("report identity differs from its curated path", result.stderr + json.dumps(summary))

    def test_empty_curation_lists_all_old_results_without_counting_them(self):
        self.curated.write_text("[]")
        _, summary = self.assert_passed(count=0)
        self.assertEqual(summary["not_in_current_curation"],
                         [{"file": str(self.result_path), "reason": "not in current curation"}])
