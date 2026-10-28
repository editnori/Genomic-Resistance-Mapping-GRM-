"""prepare_execution_roots.py helpers on a synthetic candidate root.

usage: python3 -I scripts/mtb-4k/tests/test_prepare_execution_roots.py
"""
import json
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


class Prepare(unittest.TestCase):
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
            self.assertEqual(prep.prepare_bvbrc(data, c, manifest), 2)
            bvbrc = data / "bvbrc"
            self.assertFalse(bvbrc.is_symlink())
            self.assertFalse((bvbrc / "features").is_symlink())
            self.assertTrue((bvbrc / "genomes.jsonl").is_symlink())
            # A fetched annotation lands in the lane's copy only.
            (bvbrc / "features" / "g3.json").write_text("[]")
            run = base / "run"
            (run / "frozen-inputs.json").write_text(json.dumps({"candidate": str(c), "sha256": manifest}))
            self.assertEqual(prep.verify_frozen(run), 0)
            # A write into the frozen cache is caught.
            (c / "bvbrc" / "features" / "g4.json").write_text("[]")
            self.assertEqual(prep.verify_frozen(run), 1)


if __name__ == "__main__":
    unittest.main()
