"""artifacts.py on synthetic releases: as-of resolution at the release-day
boundaries, refusal of missing hashes and releases, exact manifests, and the
unseen-gene encoder.

usage: python3 -I -B -m unittest discover -s tests/temporal -t .
"""
import datetime as dt
import hashlib
import importlib.util
import json
import random
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def load(name):
    """Load scripts/temporal/<name>.py by path under the key canaries.py uses."""
    key = f"grm_temporal_{name}"
    module = sys.modules.get(key)
    if module is None:
        spec = importlib.util.spec_from_file_location(key, ROOT / "scripts" / "temporal" / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return module


A = load("artifacts")


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


ROWS = [
    {"kind": "catalog", "release_id": "cat-v1", "released_at": "2026-01-01", "sha256": sha("cat-v1"),
     "vocabulary": ["katG_S315T", "rpoB_S450L"]},
    {"kind": "catalog", "release_id": "cat-v2", "released_at": "2026-03-01", "sha256": sha("cat-v2"),
     "vocabulary": ["katG_S315T", "rpoB_S450L", "inhA_c-15t"]},
    {"kind": "catalog", "release_id": "cat-v3", "released_at": "2026-06-01", "sha256": sha("cat-v3"),
     "vocabulary": ["katG_S315T", "inhA_c-15t", "embB_M306V"]},
    {"kind": "annotation", "release_id": "ann-1", "released_at": "2025-12-01", "sha256": sha("ann-1"),
     "vocabulary": ["geneA", "geneB"]},
    {"kind": "annotation", "release_id": "ann-2", "released_at": "2026-03-02", "sha256": sha("ann-2"),
     "vocabulary": ["geneA", "geneB", "geneC"]},
]


def registry(rows=ROWS):
    return A.ArtifactRegistry.from_records(rows)


def ids(resolution):
    return {release.kind: release.release_id for release in resolution.releases}


class Resolve(unittest.TestCase):
    def test_release_day_is_included_and_day_before_is_not(self):
        reg = registry()
        self.assertEqual(ids(reg.resolve("2026-02-28")), {"catalog": "cat-v1", "annotation": "ann-1"})
        self.assertEqual(ids(reg.resolve("2026-03-01")), {"catalog": "cat-v2", "annotation": "ann-1"})
        self.assertEqual(ids(reg.resolve("2026-03-02")), {"catalog": "cat-v2", "annotation": "ann-2"})
        self.assertEqual(ids(reg.resolve("2026-05-31")), {"catalog": "cat-v2", "annotation": "ann-2"})
        self.assertEqual(ids(reg.resolve("2026-06-01")), {"catalog": "cat-v3", "annotation": "ann-2"})
        self.assertEqual(ids(reg.resolve(dt.date(2026, 6, 2))), {"catalog": "cat-v3", "annotation": "ann-2"})

    def test_resolution_is_sorted_by_kind_and_carries_cutoff(self):
        resolution = registry().resolve("2026-03-01")
        self.assertEqual([r.kind for r in resolution.releases], ["annotation", "catalog"])
        self.assertEqual(resolution.cutoff, dt.date(2026, 3, 1))
        self.assertEqual(resolution.release("catalog").ref, "catalog:cat-v2@2026-03-01")

    def test_missing_release_is_refused(self):
        with self.assertRaisesRegex(A.ArtifactError, "catalog: no release on or before 2025-12-31"):
            registry().resolve("2025-12-31")
        with self.assertRaisesRegex(A.ArtifactError, "annotation: no release on or before 2025-11-30"):
            registry().resolve("2025-11-30")

    def test_missing_hash_is_refused_not_skipped(self):
        for missing in (None, ""):
            rows = ROWS + [{"kind": "catalog", "release_id": "cat-v4", "released_at": "2026-07-01",
                            "sha256": missing, "vocabulary": ["katG_S315T"]}]
            reg = registry(rows)
            with self.assertRaisesRegex(A.ArtifactError, r"catalog:cat-v4@2026-07-01: sha256 is missing"):
                reg.resolve("2026-07-01")
            with self.assertRaises(A.ArtifactError):
                reg.vocabulary("2026-08-01")
            with self.assertRaises(A.ArtifactError):
                reg.manifest("2026-07-01")
            # Before the hashless release, resolution is unaffected.
            self.assertEqual(ids(reg.resolve("2026-06-30"))["catalog"], "cat-v3")

    def test_registration_refuses_bad_rows(self):
        bad = [
            dict(ROWS[0], sha256="ABC"),
            dict(ROWS[0], sha256=sha("x").upper()),
            dict(ROWS[0], released_at="20260101"),
            dict(ROWS[0], released_at="2026-02-30"),
            dict(ROWS[0], released_at="2026-01-01T00:00:00"),
            dict(ROWS[0], kind="panel"),
            dict(ROWS[0], vocabulary="katG"),
            dict(ROWS[0], vocabulary=["katG", ""]),
            dict(ROWS[0], vocabulary=[A.UNKNOWN]),
            dict(ROWS[0], extra=1),
            {k: v for k, v in ROWS[0].items() if k != "sha256"},
        ]
        for row in bad:
            with self.subTest(row=row), self.assertRaises(A.ArtifactError):
                registry([row])

    def test_duplicate_ids_and_same_day_releases_are_refused(self):
        with self.assertRaisesRegex(A.ArtifactError, "duplicate catalog release_id"):
            registry(ROWS + [dict(ROWS[0], released_at="2026-09-01")])
        with self.assertRaisesRegex(A.ArtifactError, "would be ambiguous"):
            registry(ROWS + [dict(ROWS[0], release_id="cat-v1b", sha256=sha("b"))])

    def test_cutoff_must_be_a_date(self):
        reg = registry()
        for cutoff in ("2026-3-1", "2026-03-01T00:00", dt.datetime(2026, 3, 1), 20260301, None):
            with self.subTest(cutoff=cutoff), self.assertRaises(A.ArtifactError):
                reg.resolve(cutoff)

    def test_verify_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            good = Path(tmp, "catalog.tsv")
            good.write_bytes(b"katG\n")
            rows = [dict(ROWS[0], sha256=hashlib.sha256(b"katG\n").hexdigest(), path=str(good)), ROWS[3]]
            reg = registry(rows)
            with self.assertRaisesRegex(A.ArtifactError, "annotation:ann-1@2025-12-01: no path to verify"):
                reg.resolve("2026-01-01", verify_files=True)
            rows[1] = dict(ROWS[3], sha256=hashlib.sha256(b"geneA\n").hexdigest(), path=str(Path(tmp, "ann")))
            Path(tmp, "ann").write_bytes(b"geneA\n")
            self.assertEqual(ids(registry(rows).resolve("2026-01-01", verify_files=True)),
                             {"catalog": "cat-v1", "annotation": "ann-1"})
            good.write_bytes(b"katG tampered\n")
            with self.assertRaisesRegex(A.ArtifactError, "does not match its sha256"):
                registry(rows).resolve("2026-01-01", verify_files=True)
            good.unlink()
            with self.assertRaisesRegex(A.ArtifactError, "is not a file"):
                registry(rows).resolve("2026-01-01", verify_files=True)

    def test_load_from_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp, "releases.json")
            path.write_text(json.dumps(ROWS))
            self.assertEqual(A.ArtifactRegistry.load(path).manifest("2026-03-02"),
                             registry().manifest("2026-03-02"))
            path.write_text(json.dumps({"releases": ROWS}))
            with self.assertRaises(A.ArtifactError):
                A.ArtifactRegistry.load(path)


class Vocabulary(unittest.TestCase):
    def test_vocabulary_as_of(self):
        reg = registry()
        self.assertEqual(reg.vocabulary("2026-03-01", "catalog"),
                         {"katG_S315T", "rpoB_S450L", "inhA_c-15t"})
        self.assertEqual(reg.vocabulary("2026-03-01", "annotation"), {"geneA", "geneB"})
        self.assertEqual(reg.vocabulary("2026-03-01"),
                         {"katG_S315T", "rpoB_S450L", "inhA_c-15t", "geneA", "geneB"})
        self.assertNotIn("embB_M306V", reg.vocabulary("2026-05-31"))
        self.assertIn("embB_M306V", reg.vocabulary("2026-06-01"))
        # The newest release replaces, not extends: v3 retired rpoB_S450L.
        self.assertNotIn("rpoB_S450L", reg.vocabulary("2026-06-01"))
        self.assertIsInstance(reg.vocabulary("2026-03-01"), frozenset)


class Manifest(unittest.TestCase):
    def test_manifest_is_exactly_the_resolution(self):
        self.assertEqual(registry().manifest("2026-03-01"), {
            "cutoff": "2026-03-01",
            "artifacts": {
                "annotation": {"release_id": "ann-1", "released_at": "2025-12-01", "sha256": sha("ann-1")},
                "catalog": {"release_id": "cat-v2", "released_at": "2026-03-01", "sha256": sha("cat-v2")},
            }})

    def test_declared_must_match_exactly(self):
        reg = registry()
        declared = reg.manifest("2026-03-01")
        reg.check_declared("2026-03-01", json.loads(json.dumps(declared)))
        later = reg.manifest("2026-06-01")
        cases = {
            "later catalog": (dict(declared, artifacts=dict(declared["artifacts"],
                                                            catalog=later["artifacts"]["catalog"])), "catalog"),
            "wrong hash": (dict(declared, artifacts=dict(declared["artifacts"], annotation=dict(
                declared["artifacts"]["annotation"], sha256=sha("other")))), "annotation"),
            "missing kind": (dict(declared, artifacts={"catalog": declared["artifacts"]["catalog"]}), "annotation"),
            "extra kind": (dict(declared, artifacts=dict(declared["artifacts"], panel={})), "panel"),
            "wrong cutoff": (dict(declared, cutoff="2026-03-02"), "cutoff"),
            "extra key": (dict(declared, note="x"), "manifest"),
            "extra entry key": (dict(declared, artifacts=dict(declared["artifacts"], catalog=dict(
                declared["artifacts"]["catalog"], path="x"))), "catalog"),
        }
        for name, (bad, where) in cases.items():
            with self.subTest(name):
                mismatches = A.manifest_mismatches(declared, bad)
                self.assertIn(where, [w for w, _ in mismatches])
                with self.assertRaises(A.ArtifactError):
                    reg.check_declared("2026-03-01", bad)
        self.assertEqual(A.manifest_mismatches(declared, None)[0][0], "manifest")


class Encode(unittest.TestCase):
    VOCAB = frozenset({"katG", "rpoB"})

    def test_policies(self):
        self.assertEqual(A.encode_genes(["rpoB", "novel", "katG", "other"], self.VOCAB), ("katG", "rpoB", A.UNKNOWN))
        self.assertEqual(A.encode_genes(["rpoB", "novel"], self.VOCAB, "drop"), ("rpoB",))
        self.assertEqual(A.encode_genes(["novel"], self.VOCAB, "drop"), ())
        self.assertEqual(A.encode_genes(["katG"], self.VOCAB), ("katG",))
        with self.assertRaises(A.ArtifactError):
            A.encode_genes(["katG"], self.VOCAB, "keep")
        with self.assertRaises(A.ArtifactError):
            A.encode_genes("katG", self.VOCAB)

    def test_vocabulary_is_never_extended(self):
        vocab = set(self.VOCAB)
        A.encode_genes(["novel", "katG"], vocab)
        self.assertEqual(vocab, self.VOCAB)


class Determinism(unittest.TestCase):
    def test_registration_order_does_not_matter(self):
        expected = json.dumps([registry().manifest(day) for day in ("2026-03-01", "2026-06-01")])
        for seed in range(5):
            rows = list(ROWS)
            random.Random(seed).shuffle(rows)
            reg = registry(rows)
            self.assertEqual(json.dumps([reg.manifest(day) for day in ("2026-03-01", "2026-06-01")]), expected)
            self.assertEqual(reg.resolve("2026-03-01"), registry().resolve("2026-03-01"))
            self.assertEqual(reg.releases(), registry().releases())


if __name__ == "__main__":
    unittest.main()
