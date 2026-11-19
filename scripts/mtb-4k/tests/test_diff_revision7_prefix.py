"""diff_revision7_prefix.py on synthetic revision-8 reports and revision-7 counterparts.

usage: FOLDGEN=<foldgen binary> python3 -I scripts/mtb-4k/tests/test_diff_revision7_prefix.py
"""
import contextlib
import copy
import io
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import diff_revision7_prefix as diff  # noqa: E402
import fixtures as fx  # noqa: E402

SLUG = "mycobacterium_tuberculosis"


def as_revision7(report):
    """What the same evaluation wrote at revision 7: no lasso solver evidence, another lasso fit."""
    r = copy.deepcopy(report)
    r["evaluation"]["revision"] = 7
    r["models"].pop("lasso_fit")
    r["models"]["lasso"] = {"intercept": -0.5, "coef": [[3, 0.25]], "lambda": 0.02}
    for split in r["splits"]:
        for row in split["learners"]:
            if row["learner"] == "lasso":
                row.pop("lasso_fits", None)
                row["seconds"] = 99.0
                row["scores"] = [None if v is None else min(1.0, v + 0.01) for v in row["scores"]]
    for split in r["splits"]:
        for row in split["learners"]:
            row["seconds"] = row["seconds"] + 1.0
    r["seconds"] = 999.0
    r["features"] = [{"feature": 1}]
    return r


class Diff(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        rules = fx.make_rules({"prefix": ("amikacin", "isoniazid")})
        cls.r8 = fx.build_cohort(Path(cls.tmp.name) / "r8", "prefix", rules) / "data"
        # The fixtures stamp the suite's revision (MTB4K_TEST_PROFILE); the script under test compares revision 8 with 7.
        for path in (cls.r8 / "species" / SLUG / "results").glob("*.json"):
            report = json.loads(path.read_text())
            report["evaluation"]["revision"] = 8
            path.write_text(json.dumps(report))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def revision7(self, label, edit=None):
        root = Path(self.tmp.name) / label / "data"
        (root / "species" / SLUG / "results").mkdir(parents=True)
        for path in sorted((self.r8 / "species" / SLUG / "results").glob("*.json")):
            r7 = as_revision7(json.loads(path.read_text()))
            if edit and path.stem.endswith("amikacin"):
                edit(r7)
            (root / "species" / SLUG / "results" / path.name).write_text(json.dumps(r7))
        return root

    def run_diff(self, r7, new=None, profile=None):
        out = r7.parent / "diff.json"
        argv = ["--revision8", str(new or self.r8), "--revision7", str(r7), "--out", str(out)] + (["--profile", profile] if profile else [])
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as err:
            code = diff.main(argv)
        self.last_stderr = err.getvalue()
        return code, json.loads(out.read_text()) if out.exists() else None

    def test_only_the_lasso_may_differ(self):
        code, out = self.run_diff(self.revision7("same"))
        self.assertEqual(code, 0, out)
        self.assertEqual((out["datasets"], out["identical"]), (2, 2))
        self.assertTrue(all(r["lasso_rows_allowed_to_differ"] == 2 for r in out["rows"]))

    def test_a_non_lasso_change_or_another_fingerprint_is_a_difference(self):
        for label, edit, where in [
            ("cart score", lambda r: r["splits"][0]["learners"][1]["scores"].__setitem__(0, 0.5), "splits.random.learners.cart"),
            ("fingerprint", lambda r: r.update(input_fingerprint="0" * 16), "input_fingerprint"),
            ("catalog", lambda r: r["catalog"].update(version="other"), "catalog"),
            ("a sample label", lambda r: r["samples"][0].update(resistant=not r["samples"][0]["resistant"]), "samples"),
            ("a test fold", lambda r: r["splits"][1]["test_fold"].__setitem__(0, (r["splits"][1]["test_fold"][0] + 1) % 5), "splits.lineage.test_fold"),
            ("final rf importance", lambda r: r["models"].update(rf_importance=[{"feature": 1, "weight": 1.0}]), "models.rf_importance")]:
            with self.subTest(label):
                code, out = self.run_diff(self.revision7("d-" + label.replace(" ", "-"), edit))
                self.assertEqual(code, 1, label)
                row = next(r for r in out["rows"] if r["id"].endswith("amikacin"))
                self.assertTrue(any(d.startswith(where) for d in row["differences"]), row["differences"])

    def revision_n_side(self, revision):
        """The revision-8 side restamped to `revision`, as revision N's prefix would be: only lasso-derived parts may move."""
        root = Path(self.tmp.name) / f"rN-{revision}" / "data"
        if not root.exists():
            shutil.copytree(self.r8, root, symlinks=True)
            for path in (root / "species" / SLUG / "results").glob("*.json"):
                report = json.loads(path.read_text())
                report["evaluation"]["revision"] = revision
                report["features"] = [{"feature": 7, "note": "lasso top-30 annotation"}]
                for split in report["splits"]:
                    for row in split["learners"]:
                        if row["learner"] == "lasso":
                            row["scores"] = [None if v is None else max(0.0, v - 0.02) for v in row["scores"]]
                path.write_text(json.dumps(report))
        return root

    @contextlib.contextmanager
    def revision_n(self, value):
        saved = dict(diff.rev_n.VALUES)
        diff.rev_n.VALUES["REVISION_N"] = value
        try:
            yield
        finally:
            diff.rev_n.VALUES.clear()
            diff.rev_n.VALUES.update(saved)

    def test_profile_rn_takes_revision_n_and_keeps_the_field_rules(self):
        new = self.revision_n_side(9)
        with self.revision_n(9):
            code, out = self.run_diff(self.revision7("rN-same"), new, "rN")
            self.assertEqual(code, 0, self.last_stderr)
            self.assertEqual((out["profile"], out["evaluation_revision"], out["identical"]), ("rN", 9, 2))
            code, out = self.run_diff(self.revision7("rN-cart", lambda r: r["splits"][0]["learners"][1]["scores"].__setitem__(0, 0.5)), new, "rN")
            self.assertEqual(code, 1)
            row = next(r for r in out["rows"] if r["id"].endswith("amikacin"))
            self.assertTrue(any(d.startswith("splits.random.learners.cart") for d in row["differences"]), row["differences"])
            # The revision-8 side is not revision N, and revision N is not the r8 profile's.
            self.assertEqual(self.run_diff(self.revision7("rN-r8side"), self.r8, "rN")[0], 2)
            self.assertIn("expected revision 9 against revision 7", self.last_stderr)
        self.assertEqual(self.run_diff(self.revision7("r8-on-rN"), new, "r8")[0], 2)
        self.assertIn("expected revision 8 against revision 7", self.last_stderr)

    def kset_side(self, label, revision, field):
        """A revision-`revision` side whose counted.json gives two samples one kset digest; field(i, j) is what each
        report lists (None: leave the field out)."""
        root = Path(self.tmp.name) / label / "data"
        shutil.copytree(self.r8, root, symlinks=True)
        sp = root / "species" / SLUG
        first = json.loads(next((sp / "results").glob("*amikacin.json")).read_text())
        a, b = first["samples"][0]["id"], first["samples"][1]["id"]
        counted = json.loads((sp / "counted.json").read_text())
        by = {g["genome_id"]: g for g in counted["genomes"]}
        by[b]["kset_digest"] = by[a]["kset_digest"]
        (sp / "counted.json").write_text(json.dumps(counted))
        for path in (sp / "results").glob("*.json"):
            report = json.loads(path.read_text())
            report["evaluation"]["revision"] = revision
            ids = [s["id"] for s in report["samples"]]
            value = field(ids.index(a), ids.index(b))
            if value is not None:
                report["kset_duplicates"] = value
            path.write_text(json.dumps(report))
        return root

    def test_kset_duplicates_only_from_the_threshold_and_only_as_the_digests_give(self):
        ten = diff.DUPLICATE_ISOLATION_REVISION
        with self.revision_n(ten):
            code, out = self.run_diff(self.revision7("k-ok"), self.kset_side("k-ok-new", ten, lambda i, j: [[i, j]]), "rN")
            self.assertEqual(code, 0, self.last_stderr)
            for label, field in (("k-missing", lambda i, j: None), ("k-wrong", lambda i, j: [[i, j + 1]])):
                code, out = self.run_diff(self.revision7(label), self.kset_side(label + "-new", ten, field), "rN")
                self.assertEqual(code, 1, label)
                self.assertIn("kset_duplicates (differs from the root's counted.json kset digests)", out["rows"][0]["differences"])
        with self.revision_n(ten - 1):
            code, out = self.run_diff(self.revision7("k-below"), self.kset_side("k-below-new", ten - 1, lambda i, j: [[i, j]]), "rN")
            self.assertEqual(code, 1)
            self.assertIn(f"kset_duplicates (not before revision {ten})", out["rows"][0]["differences"])
        r7 = self.revision7("k-on-7", lambda r: r.update(kset_duplicates=[[0, 1]]))
        code, out = self.run_diff(r7)
        self.assertEqual(code, 1)
        self.assertIn("kset_duplicates (on the revision-7 side)", next(r for r in out["rows"] if r["id"].endswith("amikacin"))["differences"])

    def test_an_unfilled_rn_refuses(self):
        new = self.revision_n_side(9)
        with self.revision_n("@@REVISION_N@@"):
            code, out = self.run_diff(self.revision7("rN-unfilled"), new, "rN")
        self.assertEqual((code, out), (2, None))
        self.assertIn("diff_revision7_prefix: STOP: unfilled placeholder(s): @@REVISION_N@@", self.last_stderr)

    def test_wrong_revisions_or_dataset_sets_stop(self):
        r7 = self.revision7("wrong-revision", lambda r: r["evaluation"].update(revision=8))
        self.assertEqual(self.run_diff(r7)[0], 2)
        r7 = self.revision7("missing-dataset")
        next((r7 / "species" / SLUG / "results").glob("*isoniazid.json")).unlink()
        self.assertEqual(self.run_diff(r7)[0], 2)


if __name__ == "__main__":
    unittest.main()
