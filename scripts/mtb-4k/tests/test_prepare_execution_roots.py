"""prepare_execution_roots.py on a synthetic candidate root and cohort.

usage: python3 -I scripts/mtb-4k/tests/test_prepare_execution_roots.py
"""
import contextlib
import csv
import gzip
import hashlib
import io
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import prepare_execution_roots as prep  # noqa: E402


def candidate(base: Path) -> Path:
    c = base / "candidate"
    for sub in ["catalog", "curated", "bvbrc/features"]:
        (c / sub).mkdir(parents=True)
    (c / "catalog" / "version.txt").write_text("v1")
    (c / "curated" / "datasets.json").write_text("[]")
    for name in prep.BVBRC_RECORDS:
        (c / "bvbrc" / name).write_text(f"{name}\n")
    for g in ["g1", "g2"]:
        (c / "bvbrc" / "features" / f"{g}.json").write_text(f'[{{"gene": "{g}"}}]')
    return c


def cohort(base: Path):
    """4,000 one-line FASTA in two stores, the cohort evidence files, and a
    candidate whose annotation cache holds 40 cohort genomes and 30 others."""
    ids = [f"1773.{100000 + i}" for i in range(4000)]
    stores = [base / "store-a", base / "store-b"]
    for s in stores:
        s.mkdir()
    rows = []
    for i, gid in enumerate(ids):
        path = stores[i % 2] / f"{gid}.fna.gz"
        path.write_bytes(gzip.compress(f">c\nACGT{i}\n".encode(), mtime=0))
        b = path.read_bytes()
        rows.append({"genome_id": gid, "size_bytes": len(b), "sha256": hashlib.sha256(b).hexdigest(),
                     "canonical_source_path": str(path), "read_hash_reproducible": "True",
                     "unchanged_since_resolution": "True"})
    ev = base / "evidence"
    (ev / "results").mkdir(parents=True)
    (ev / "proofs").mkdir()
    (ev / "results" / "genomes-4000.txt").write_text("".join(f"{g}\n" for g in ids))
    (ev / "proofs" / "frozen-candidate-2000.txt").write_text("".join(f"{g}\n" for g in ids[:2000]))
    (ev / "results" / "added-2000.txt").write_text("".join(f"{g}\n" for g in ids[2000:]))
    with open(ev / "results" / "fasta-manifest.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    c = candidate(base)
    for gid in ids[::100] + [f"other.{i}" for i in range(30)]:
        (c / "bvbrc" / "features" / f"{gid}.json").write_text("x" * 1000)
    return ids, stores, ev, c


class Prepare(unittest.TestCase):
    def setUp(self):
        self.saved = (prep.EVIDENCE, prep.CANDIDATE, prep.ALLOWED_FASTA_PARENTS, prep.EXPECTED, prep.free_bytes)

    def tearDown(self):
        prep.EVIDENCE, prep.CANDIDATE, prep.ALLOWED_FASTA_PARENTS, prep.EXPECTED, prep.free_bytes = self.saved

    def redirect(self, base: Path):
        ids, stores, ev, c = cohort(base)
        h = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()  # noqa: E731
        prep.EVIDENCE, prep.CANDIDATE = ev, c
        prep.ALLOWED_FASTA_PARENTS = set(stores)
        prep.EXPECTED = {"full": (h(ev / "results/genomes-4000.txt"), 4000),
                         "prefix": (h(ev / "proofs/frozen-candidate-2000.txt"), 2000),
                         "added": (h(ev / "results/added-2000.txt"), 2000)}
        return ids, stores, ev, c

    def main(self, *argv):
        saved = sys.argv
        sys.argv = ["prepare_execution_roots.py", *map(str, argv)]
        out = io.StringIO()
        try:
            with contextlib.redirect_stdout(out):
                return prep.main(), out.getvalue()
        finally:
            sys.argv = saved

    def test_matched_selection_is_half_fixed_and_in_frozen_order(self):
        ids = [f"1773.{i}" for i in range(4000)]
        chosen = prep.matched_selection(ids)
        self.assertEqual(len(chosen), 2000)
        self.assertEqual(chosen, [g for g in ids if g in set(chosen)])
        self.assertEqual(chosen, prep.matched_selection(list(ids)))
        first, second = set(ids[:2000]), set(ids[2000:])
        # A hash-chosen half draws from both halves of the frozen order.
        self.assertTrue(800 < len(first & set(chosen)) < 1200, len(first & set(chosen)))
        self.assertTrue(800 < len(second & set(chosen)) < 1200)

    def test_bvbrc_is_lane_local_and_frozen_inputs_verify(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            c = candidate(base)
            manifest = prep.frozen_manifest(c)
            self.assertEqual(len(manifest), 2 + 2 + 2)
            data = base / "run" / "cohorts" / "full" / "data"
            data.mkdir(parents=True)
            self.assertEqual(prep.prepare_bvbrc(data, c, manifest, ["g1", "g2", "g3"])["files"], 2)
            bvbrc = data / "bvbrc"
            self.assertFalse(bvbrc.is_symlink())
            self.assertFalse((bvbrc / "features").is_symlink())
            self.assertTrue((bvbrc / "genomes.jsonl").is_symlink())
            # A fetched annotation lands in the lane's copy only.
            (bvbrc / "features" / "g3.json").write_text("[]")
            run = base / "run"
            (run / "frozen-inputs.json").write_text(json.dumps(
                {"candidate": str(c), "sha256": manifest, "fasta": {}, "cohort_genomes_sha256": {}}))
            self.assertEqual(prep.verify_frozen(run), 0)
            # A write into the frozen cache is caught.
            (c / "bvbrc" / "features" / "g4.json").write_text("[]")
            self.assertEqual(prep.verify_frozen(run), 1)

    def test_only_a_training_cohorts_own_annotations_are_copied(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            c = candidate(base)
            manifest = prep.frozen_manifest(c)
            for name, ids, want in [("one", ["g2", "absent"], ["g2.json"]), ("pilot", [], [])]:
                data = base / name
                data.mkdir()
                got = prep.prepare_bvbrc(data, c, manifest, ids)
                self.assertEqual(sorted(os.listdir(data / "bvbrc" / "features")), want)
                self.assertEqual(got["bytes"], sum((c / "bvbrc/features" / f).stat().st_size for f in want))

    def test_roots_copy_each_cohorts_annotations_and_verify_fasta_and_links(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            ids, stores, ev, c = self.redirect(base)
            run = base / "run"
            rc, _ = self.main(run)
            self.assertEqual(rc, 0)
            receipt = json.loads((run / "preparation.json").read_text())
            seeded = receipt["features_seeded_per_cohort"]
            cached = set(ids[::100])
            prefix = set(ids[:2000])
            matched = set(prep.matched_selection(ids))
            self.assertEqual(seeded["full"]["files"], 40)
            self.assertEqual(seeded["prefix"]["files"], len(cached & prefix))
            self.assertEqual(seeded["added"]["files"], len(cached - prefix))
            self.assertEqual(seeded["matched2k"]["files"], len(cached & matched))
            self.assertEqual(seeded["pilot256"], {"files": 0, "bytes": 0})
            self.assertEqual(seeded["pilot256-stock"], {"files": 0, "bytes": 0})
            self.assertEqual(receipt["disk"]["annotation_bytes_copied"], sum(s["bytes"] for s in seeded.values()))
            self.assertEqual(receipt["fasta_sha256_checked"], 4000)
            for name in receipt["cohorts"]:
                rc, out = self.main("--check-lane", run / "cohorts" / name / "data")
                self.assertEqual(rc, 0, out)
            self.assertEqual(self.main("--verify-frozen", run)[0], 0)
            # A FASTA edited in the linked store is a stop.
            victim = Path(stores[1]) / f"{ids[1]}.fna.gz"
            original = victim.read_bytes()
            victim.write_bytes(original + b"\0")
            rc, out = self.main("--verify-frozen", run)
            self.assertEqual(rc, 1)
            self.assertEqual(json.loads(out)["fasta_changed"], [str(victim)])
            victim.write_bytes(original)
            # So is a FASTA link pointed elsewhere, even at identical bytes.
            link = run / "cohorts" / "added" / "data" / "fasta" / f"{ids[2001]}.fna.gz"
            twin = base / "twin.fna.gz"
            twin.write_bytes(Path(os.readlink(link)).read_bytes())
            link.unlink()
            os.symlink(twin, link)
            rc, out = self.main("--verify-frozen", run)
            self.assertEqual(rc, 1)
            self.assertEqual(json.loads(out)["links_changed"], [str(link)])
            rc, out = self.main("--check-lane", run / "cohorts" / "added" / "data")
            self.assertEqual(rc, 1, out)
            # A cohort list edited after preparation is a stop.
            genomes = run / "cohorts" / "prefix" / "data" / "species" / prep.SLUG / "genomes.txt"
            genomes.write_text(genomes.read_text() + "extra\n")
            rc, out = self.main("--verify-frozen", run)
            self.assertIn("cohort list changed", " ".join(json.loads(out)["links_changed"]))
            # An existing run root is refused.
            with self.assertRaises(SystemExit):
                self.main(run)

    def test_the_disk_floor_stops_preparation_before_anything_is_written(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            self.redirect(base)
            run = base / "run"
            # 40 cohort files of 1,000 bytes go to full, about 20 each to
            # prefix, added and matched2k: one byte short of the floor stops.
            copy_bytes = sum(p.stat().st_size for name in prep.TRAINING
                             for p in prep.cached_features(prep.CANDIDATE, self.cohort_ids(name)))
            prep.free_bytes = lambda _: prep.MIN_FREE_AFTER + copy_bytes - 1
            with self.assertRaises(SystemExit) as stop:
                self.main(run)
            self.assertIn("below the", str(stop.exception))
            self.assertFalse(os.path.lexists(run))
            prep.free_bytes = lambda _: prep.MIN_FREE_AFTER + copy_bytes
            self.assertEqual(self.main(run)[0], 0)

    def cohort_ids(self, name):
        ids = (prep.EVIDENCE / "results/genomes-4000.txt").read_text().splitlines()
        return {"prefix": ids[:2000], "added": ids[2000:], "full": ids, "matched2k": prep.matched_selection(ids)}[name]

    def test_check_lane_refuses_links_where_grm_writes(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            ids, _, _, _ = self.redirect(base)
            run = base / "run"
            self.assertEqual(self.main(run)[0], 0)
            data = run / "cohorts" / "full" / "data"
            outside = base / "outside"
            outside.mkdir()
            for planted in [data / "report", data / "species" / prep.SLUG / "results",
                            data / "bvbrc" / "features" / f"{ids[0]}.json"]:
                if planted.exists():
                    planted.unlink()
                os.symlink(outside, planted)
                rc, out = self.main("--check-lane", data)
                self.assertEqual(rc, 1, planted)
                self.assertIn(str(planted), out)
                planted.unlink()
            fifo = data / "species" / prep.SLUG / "pipe"
            os.mkfifo(fifo)
            self.assertEqual(self.main("--check-lane", data)[0], 1)


if __name__ == "__main__":
    unittest.main()
