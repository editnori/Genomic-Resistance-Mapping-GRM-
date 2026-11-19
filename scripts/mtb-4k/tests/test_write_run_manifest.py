"""write_run_manifest.py on the canonical fixtures.

usage: FOLDGEN=<foldgen binary> python3 -I scripts/mtb-4k/tests/test_write_run_manifest.py

The fixture builder writes its own manifest; here each root's is removed and rewritten by the
script from the root alone, which must give the same manifest, pass report_contract, and refuse
a stale report, a missing report or skip line, and an existing manifest.
"""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import profile_boot as pb  # noqa: E402  (MTB4K_TEST_PROFILE, before report_contract loads)
import fixtures as fx  # noqa: E402
import report_contract as vc  # noqa: E402
import write_run_manifest as wm  # noqa: E402


class Write(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.rules = fx.make_rules()

    def tearDown(self):
        self.tmp.cleanup()

    def root(self, name, mutate=None):
        """A fixture root, optionally with one report mutated, and without its manifest."""
        root = fx.build_cohort(Path(self.tmp.name) / name, name, self.rules)
        if mutate:
            fx.mutate_report(root, vc.DRUG_IDS["amikacin"], mutate)
        built = json.loads((root / "validation-manifest.json").read_text())
        (root / "validation-manifest.json").unlink()
        return root, built

    def write(self, root, name, log=True):
        argv = [str(root), "--cohort", name, "--source-commit", fx.SOURCE, "--binary", str(root / "synthetic-binary")]
        if log:
            argv += ["--log", str(root / "pipeline.log")]
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = wm.main(pb.argv(argv), rules=self.rules)
        return code, err.getvalue()

    def test_each_cohorts_manifest_is_rebuilt_from_the_root_alone(self):
        for name in fx.COHORTS:
            with self.subTest(name):
                root, built = self.root(name)
                self.assertEqual(self.write(root, name), (0, ""))
                written = json.loads((root / "validation-manifest.json").read_text())
                self.assertEqual(written["reports"], built["reports"])
                self.assertEqual(sorted(written["class_gate_unavailable"], key=lambda r: r["id"]),
                                 sorted(built["class_gate_unavailable"], key=lambda r: r["id"]))
                self.assertEqual({e["path"] for e in written["inputs"]}, {e["path"] for e in built["inputs"]})
                vc.validate_manifest(root / "validation-manifest.json", self.rules)
                self.assertFalse((root / "validation-manifest.json.part").exists())
                code, err = self.write(root, name)
                self.assertEqual(code, 2)
                self.assertIn("written once", err)

    def test_a_stale_report_is_refused(self):
        root, _ = self.root("prefix", lambda r: r.update(input_fingerprint="0" * 16))
        code, err = self.write(root, "prefix")
        self.assertEqual(code, 2)
        self.assertIn("stale report", err)
        self.assertFalse((root / "validation-manifest.json").exists())

    def test_a_missing_report_or_skip_line_is_refused(self):
        root, _ = self.root("added")
        log = root / "pipeline.log"
        log.write_text("\n".join(line for line in log.read_text().splitlines() if "amoxicillin" not in line) + "\n")
        code, err = self.write(root, "added")
        self.assertEqual(code, 2)
        self.assertIn("no class-gate skip line", err)
        (root / "data/species" / vc.SPECIES / "results" / f"{vc.DRUG_IDS['amikacin']}.json").unlink()
        code, err = self.write(root, "added", log=False)
        self.assertEqual(code, 2)
        self.assertIn("there is no report", err)

    def test_a_manifest_that_fails_the_contract_is_not_left_behind(self):
        root, _ = self.root("prefix", lambda r: r["splits"][0]["learners"][1]["metrics"].update(balanced_accuracy=0.5))
        code, err = self.write(root, "prefix")
        self.assertEqual(code, 2)
        self.assertIn("balanced_accuracy", err)
        self.assertFalse((root / "validation-manifest.json").exists())
        self.assertFalse((root / "validation-manifest.json.part").exists())


if __name__ == "__main__":
    unittest.main()
