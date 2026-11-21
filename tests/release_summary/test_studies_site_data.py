"""Tests for scripts/release-summary/studies_site_data.py and web/src/data/additional-studies.json."""
from __future__ import annotations

import json
import pathlib
import re
import subprocess
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "release-summary" / "studies_site_data.py"
BUNDLED = ROOT / "web" / "src" / "data" / "additional-studies.json"


class StudiesSiteDataTest(unittest.TestCase):
    def test_bundled_json_matches_generator(self) -> None:
        res = subprocess.run(
            [sys.executable, "-I", str(SCRIPT), "--check", str(BUNDLED)],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(res.returncode, 0, res.stderr)

    def test_headline_counts_and_structure(self) -> None:
        data = json.loads(BUNDLED.read_text(encoding="utf-8"))
        self.assertEqual(data["evaluation_revision"], 14)
        self.assertEqual(data["time_split"]["datasets"], 24)
        self.assertEqual(data["time_split"]["species"], 4)
        self.assertEqual(data["time_split"]["total_log_loss_worse"], 5)
        self.assertEqual(data["time_split"]["total_balanced_accuracy_worse"], 5)
        self.assertEqual(data["time_split"]["total_recall_worse"], 2)

        self.assertEqual(data["kover"]["datasets"], 12)
        self.assertEqual(data["kover"]["genomes"], 1118)
        self.assertEqual(data["kover"]["balanced_accuracy"]["lower"], 8)

        self.assertEqual(data["dense_network"]["datasets"], 24)
        self.assertEqual(data["dense_network"]["species"], 12)
        self.assertEqual(data["dense_network"]["outcome"], "no gain established")
        self.assertEqual(data["dense_network"]["vs_best5"]["balanced_accuracy"]["positive"], 1)

        self.assertEqual(data["amrfinderplus"]["datasets"], 235)
        self.assertEqual(data["amrfinderplus"]["species"], 21)
        self.assertEqual(data["amrfinderplus"]["program_vs_known_gene_rule"]["ahead"], 94)
        self.assertEqual(data["amrfinderplus"]["program_vs_known_gene_rule"]["ahead_by_10pts"], 48)
        self.assertEqual(data["amrfinderplus"]["program_vs_known_gene_rule"]["ahead_by_10pts_quinolone"], 20)

        self.assertEqual(data["tbprofiler"]["by_name_endpoints"], 16)
        self.assertEqual(data["tbprofiler"]["graded_endpoints"], 14)
        self.assertEqual(data["tbprofiler"]["by_name_rows_with_call"], 20918)
        self.assertEqual(data["tbprofiler"]["by_name_rows_without_call"], 819)
        self.assertEqual(data["tbprofiler"]["graded_vs_learners_lineage"]["scm"]["tool_higher"], 8)
        self.assertEqual(data["tbprofiler"]["graded_vs_learners_lineage"]["rf"]["tool_higher"], 13)

    def test_no_genome_ids_sequences_or_private_paths(self) -> None:
        raw = BUNDLED.read_text(encoding="utf-8")
        self.assertTrue(raw.endswith("\n"))
        self.assertNotIn("\r", raw)
        self.assertIsNone(re.search(r"[ACGT]{20,}", raw))
        self.assertIsNone(re.search(r"\b(SAMN|SAMEA|SAMD|GCA_|GCF_)\d+", raw))
        self.assertIsNone(re.search("/" + "home/" + "|/" + "Users/" + "|C:\\\\" + "Users", raw))


if __name__ == "__main__":
    unittest.main()
