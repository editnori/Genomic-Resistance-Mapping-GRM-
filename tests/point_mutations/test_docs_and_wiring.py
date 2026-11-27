"""Documentation, dependency hygiene, and CI wiring checks for Track B."""

from __future__ import annotations

import pathlib
import re
import unittest

from tests.point_mutations.common import (
    CI_DOC,
    FOUR_TRACK_DOC,
    INDEX_DOC,
    PM_PLAN_DOC,
    ROOT,
    SCRIPTS_DIR,
    SUITE,
    WORKFLOW,
)


def suite_count() -> int:
    loader = unittest.TestLoader()
    return loader.discover(start_dir=str(ROOT / SUITE), top_level_dir=str(ROOT)).countTestCases()


class DocsAndHygieneTests(unittest.TestCase):
    def test_four_track_plan_contains_all_locked_rules(self):
        text = FOUR_TRACK_DOC.read_text(encoding="utf-8")
        self.assertNotIn("\u2014", text)
        for token in (
            "EVALUATION_REVISION",
            "REPORT_SCHEMA",
            "LearnerResult.conformal == []",
            "genomes.jsonl",
            "check-results.py",
            "amrfinderplus",
            "resfinder_pointfinder",
            "card_rgi",
            "staramr",
            "kleborate",
            "tbprofiler",
            "mykrobe",
            "kover",
            "xgboost",
            "lightgbm",
            "elasticnet",
            "AlphaGenome 2",
            "ESM-2",
            "Evo 2",
            "ProkBERT",
            "Bacformer",
            "fnv64_k31",
            "sha256_k31",
            "inner_pick.rs",
            "d_inner",
            "d_outer",
            "c = 0.85",
            "{1.00, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50}",
            "floor(c * n)",
            "coverage_r",
            "coverage_s",
            "FN / n_r",
        ):
            self.assertIn(token, text, f"missing locked rule token in four-track-plan.md: {token}")

    def test_point_mutation_plan_contains_required_protocol_sections(self):
        text = PM_PLAN_DOC.read_text(encoding="utf-8")
        self.assertNotIn("\u2014", text)
        for token in (
            "Not run yet.",
            "mechanism-split.md",
            "mechanism-split-result.md",
            "outside-learners.md",
            "dataset-export.md",
            "fold-file.md",
            "Mycobacterium tuberculosis",
            "grm-rep-exposure-v1",
            "grm-rep-compare-v1",
            "AlphaGenome 2",
            "retrospective_only",
            "certified_unseen_sequence",
            "unseen_variant_count",
            "ProteinVariantEncoder",
            "BurdenEncoder",
            "FoundationScoreAdapter",
        ):
            self.assertIn(token, text, f"missing protocol token in point-mutation-plan.md: {token}")

    def test_docs_index_links_both_new_docs(self):
        idx = INDEX_DOC.read_text(encoding="utf-8")
        self.assertIn("[four-track-plan.md](four-track-plan.md)", idx)
        self.assertIn("[point-mutation-plan.md](point-mutation-plan.md)", idx)
        self.assertNotIn("\u2014", idx)

    def test_scripts_and_tests_use_stdlib_and_numpy_only_with_no_skips(self):
        allowed = {
            "__future__",
            "argparse",
            "copy",
            "hashlib",
            "importlib",
            "json",
            "math",
            "numpy",
            "os",
            "pathlib",
            "re",
            "subprocess",
            "sys",
            "tempfile",
            "typing",
            "unittest",
            "tests",
        }
        paths = sorted(SCRIPTS_DIR.glob("*.py")) + sorted(pathlib.Path(ROOT / SUITE).glob("*.py"))
        for path in paths:
            source = path.read_text(encoding="utf-8")
            with self.subTest(file=path.name):
                imports = set(re.findall(r"^(?:import|from) ([a-zA-Z_]+)", source, re.M))
                self.assertLessEqual(imports, allowed, f"{path.name} has disallowed imports: {imports - allowed}")
                self.assertNotIn("url" + "lib", source)
                self.assertNotIn("sock" + "et", source)
                self.assertIsNone(re.search(r"[\"'](/(home|Users|mnt|tmp|root|opt|srv)/|[A-Za-z]:\\)", source))
                if path.parent.name == "point_mutations":
                    self.assertNotIn("skip" + "Test", source)
                    self.assertNotIn("unittest." + "skip", source)

    def test_ci_and_ci_doc_wire_suite(self):
        wf = WORKFLOW.read_text(encoding="utf-8")
        job = wf.split("\n  python:\n", 1)[1].split("\n    steps:", 1)[0]
        self.assertIn(f"          - suite: {SUITE}\n            top: .\n", job)
        self.assertEqual(job.count(f"suite: {SUITE}\n"), 1)

        ci_text = CI_DOC.read_text(encoding="utf-8")
        rows = [line for line in ci_text.split("\n") if line.startswith("| `python`")]
        self.assertEqual(len(rows), 1)
        m = re.search(rf"`{re.escape(SUITE)}` (\d+)(,| \|)", rows[0])
        self.assertIsNotNone(m, "docs/ci.md python row does not list tests/point_mutations")
        assert m is not None
        self.assertEqual(int(m.group(1)), suite_count())
        self.assertEqual(ci_text.count(f"The `{SUITE}` suite"), 1)


if __name__ == "__main__":
    unittest.main()
