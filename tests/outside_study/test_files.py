"""The scoring script end to end on files: the pins, the frozen gate, the
outputs and byte-identical reruns."""

import contextlib
import io
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests.outside_study import fixtures as F
from tests.outside_study.fixtures import S

OUTPUTS = ("genomes.tsv", "datasets.tsv", "metrics.tsv", "reading_counts.tsv", "groups.tsv",
           "leave_one_study_out.tsv", "features.tsv", "reading.txt", "results.json")


class Files(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.kit = F.frozen_kit(self.root / "kit")
        self.paths = F.write_world(F.world(), self.root / "world")

    def tearDown(self):
        self.tmp.cleanup()

    def overrides(self):
        return {"pins": self.paths["pins"], "table_pins": self.paths["table_pins"], "kit_root": self.kit}

    def argv(self, out, *extra):
        p = self.paths
        return ["--work", str(p["work"]), "--tables", str(p["tables"]), "--reports", str(p["reports"]),
                "--species", F.SPECIES, "--population", "primary", "--out", str(out), *extra]

    def main(self, out, *extra, overrides=None):
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            code = S.main(self.argv(out, *extra), overrides=overrides or self.overrides())
        return code, err.getvalue()

    def test_scores_and_writes_every_output(self):
        out = self.root / "out"
        code, err = self.main(out)
        self.assertEqual(code, 0, err)
        self.assertEqual(sorted(p.name for p in out.iterdir()), sorted(OUTPUTS))
        self.assertIn("123 of 126 genomes scored on 2 datasets", err)
        res = (out / "results.json").read_text()
        self.assertIn(F.git(self.kit, "rev-parse", "HEAD").stdout.strip(), res)
        metrics = S.read_tsv(out / "metrics.tsv")
        scm = next(r for r in metrics if r["dataset"] == F.A and r["learner"] == "scm")
        self.assertEqual((scm["tp"], scm["fn"], scm["fp"], scm["tn"]), ("48", "12", "7", "53"))
        self.assertEqual((scm["sensitivity"], scm["sensitivity_lo"], scm["reading_sensitivity"]),
                         ("0.800000", "0.676700", "below"))
        self.assertEqual(S.read_tsv(out / "reading_counts.tsv")[0],
                         {"learner": "scm", "below": "1", "consistent": "1", "above": "0", "not_read": "2"})

    def test_reruns_are_byte_identical(self):
        a, b = self.root / "a", self.root / "b"
        self.assertEqual(self.main(a)[0], 0)
        self.assertEqual(self.main(b)[0], 0)
        for name in OUTPUTS:
            self.assertEqual((a / name).read_bytes(), (b / name).read_bytes(), name)

    def test_reruns_are_byte_identical_across_hash_seeds(self):
        a = self.root / "a"
        self.assertEqual(self.main(a)[0], 0)
        code = (
            "import sys; sys.path.insert(0, sys.argv[1]);"
            "from tests.outside_study import fixtures as F;"
            "import json; o = json.loads(sys.argv[2]);"
            "sys.exit(F.S.main(o['argv'], overrides=o['overrides']))"
        )
        import json
        for seed in ("0", "1", "4242"):
            out = self.root / f"seed{seed}"
            payload = json.dumps({"argv": self.argv(out),
                                  "overrides": {k: str(v) for k, v in self.overrides().items()}})
            env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"), "PYTHONHASHSEED": seed}
            p = subprocess.run([sys.executable, "-I", "-c", code, str(F.ROOT), payload],
                               env=env, capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            for name in OUTPUTS:
                self.assertEqual((a / name).read_bytes(), (out / name).read_bytes(), f"{name} seed {seed}")

    def refused(self, pattern, *extra, overrides=None):
        out = self.root / "refused"
        code, err = self.main(out, *extra, overrides=overrides)
        self.assertEqual(code, 2, err)
        self.assertRegex(err, pattern)
        self.assertFalse(out.exists(), "a refusal writes nothing")

    def test_the_frozen_tables_must_match_their_pins(self):
        with open(self.paths["tables"] / "pairs.tsv", "a") as f:
            f.write("\t".join(["1001.1", F.C, "S"] + [""] * 9) + "\n")
        self.refused("the frozen table pairs.tsv has sha256")

    def test_reports_must_match_their_pins(self):
        p = self.paths["reports"] / f"{F.A}.json"
        p.write_bytes(p.read_bytes().replace(b'"schema": 2', b'"schema":  2'))
        self.refused("is not the pinned")

    def test_a_cohort_genome_in_the_frozen_list_is_refused(self):
        rows = S.read_tsv(self.paths["tables"] / "candidates.tsv")
        rows[0]["genome_id"] = "999.3"
        S.write_tsv(self.paths["tables"] / "candidates.tsv", F.CAND_COLS, rows)
        o = self.overrides()
        o["table_pins"] = self.root / "tables-repinned.sha256"
        o["table_pins"].write_text(
            "".join(f"{S.sha256_file(self.paths['tables'] / n)}  {n}\n" for n in S.TABLE_FILES))
        self.refused("999.3 is a cohort genome", overrides=o)

    def test_known_gene_unavailable_only_in_a_rehearsal(self):
        work = self.paths["work"]
        (work / "known-gene" / "unavailable.txt").write_text("no catalog\n")
        self.refused("only a rehearsal may go on without them")

    def test_a_rehearsal_takes_cohort_genomes_only(self):
        self.refused("is not a cohort genome", "--rehearsal")

    def test_rows_naming_another_report(self):
        p = self.paths["work"] / "predict" / f"{F.A}.jsonl"
        p.write_text(p.read_text().replace(F.sha(F.report_bytes(F.make_report(F.A))), "0" * 64, 1))
        self.refused("names report sha256")

    def test_no_option_moves_the_pins(self):
        for flag in ("--pins", "--kit-root"):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                S.main(self.argv(self.root / "x") + [flag, "elsewhere"])


class CommittedPins(unittest.TestCase):
    """The repository's own pin file against the repository's own files, so
    a change to a frozen file that leaves its pin behind fails here and not
    first on the box."""

    def test_the_pin_file_names_the_frozen_files_and_matches_them(self):
        pins = S.read_pins(F.ROOT / "scripts/outside-study/FROZEN.sha256")
        self.assertEqual(sorted(pins), sorted(S.FROZEN_FILES))
        for rel in S.FROZEN_FILES:
            self.assertEqual(pins[rel], S.sha256_file(F.ROOT / rel), rel)


class FrozenGate(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.kit = F.frozen_kit(Path(self.tmp.name) / "kit")

    def tearDown(self):
        self.tmp.cleanup()

    def test_passes_on_the_committed_files(self):
        gate = S.frozen_gate(self.kit)
        self.assertEqual(gate["commit"], F.git(self.kit, "rev-parse", "HEAD").stdout.strip())
        self.assertEqual(sorted(gate["files"]), sorted(S.FROZEN_FILES))

    def test_a_changed_file_is_refused(self):
        (self.kit / "scripts/outside-study/score.py").write_text("score, changed\n")
        with self.assertRaisesRegex(S.Refusal, "score.py has sha256"):
            S.frozen_gate(self.kit)

    def test_a_changed_plan_is_refused(self):
        (self.kit / "docs/outside-study-test.md").write_text("another plan\n")
        with self.assertRaisesRegex(S.Refusal, "outside-study-test.md has sha256"):
            S.frozen_gate(self.kit)

    def test_pins_changed_to_match_but_not_committed(self):
        (self.kit / "scripts/outside-study/run.py").write_text("run, changed\n")
        (self.kit / "scripts/outside-study/FROZEN.sha256").write_text("".join(
            f"{S.sha256_file(self.kit / rel)}  {rel}\n" for rel in S.FROZEN_FILES))
        with self.assertRaisesRegex(S.Refusal, "differ from the commit"):
            S.frozen_gate(self.kit)

    def test_files_never_committed(self):
        F.git(self.kit, "rm", "-q", "--cached", "docs/outside-study-test.md")
        F.git(self.kit, "commit", "-q", "-m", "drop")
        with self.assertRaisesRegex(S.Refusal, "not all committed"):
            S.frozen_gate(self.kit)

    def test_no_frozen_file(self):
        (self.kit / "scripts/outside-study/FROZEN.sha256").unlink()
        with self.assertRaisesRegex(S.Refusal, "not frozen"):
            S.frozen_gate(self.kit)


if __name__ == "__main__":
    unittest.main()
