"""pilot_invalidation.py against the real GRM binary (Linux), on copies of a real interrupted root.

usage: GRM=target/release/grm python3 -I scripts/mtb-4k/tests/test_pilot_invalidation.py

The interrupted root is made as the plan makes it: a disk-sparse run at width 64 over 128 genomes,
held after block 0 folds (genome 64's FASTA is a pipe), interrupted by pilot_interrupt.py under a
guard stand-in, then copied with `cp -a` into the four copies the plan takes (L5). In each copy the
pipe is replaced by the genome's real bytes, so a resume that is not refused would run to the end
instead of blocking. The helper must report all four checks refused, change nothing in any copy,
and fail when a check is not refused (the width check at the run's own width) or a copy is not
the interrupted state.
"""
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import pilot_invalidation as piv  # noqa: E402
import pilot_interrupt as pi  # noqa: E402
import test_pilot_interrupt as tpi  # noqa: E402

SLUG = tpi.SLUG
HELPER = HERE.parent / "pilot_invalidation.py"


def invalidate(copies, receipt, width=64):
    return subprocess.run([sys.executable, "-I", "-B", str(HELPER), "--copies", str(copies), "--grm", tpi.grm(), "--width", str(width),
                           "--receipt", str(receipt), "--timeout", "300", "--grm-arg=--allow-missing-catalog"], capture_output=True, text=True)


class Invalidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        base = Path(cls.tmp.name)
        cohort = base / "pilot"
        root, ids, texts = tpi.build_root(cohort)
        fifo = root / "fasta" / f"{ids[64]}.fna.gz"
        fifo.unlink()
        os.mkfifo(fifo)
        resources = base / "guard.json"
        wrapper = subprocess.Popen([sys.executable, "-I", "-B", "-c", tpi.WRAPPER, str(resources), "", *tpi.sparse(root)],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        p = tpi.helper(root, base / "interrupt.json", wrapper.pid, resources)
        wrapper.wait(timeout=60)
        assert p.returncode == 0, p.stderr
        cls.copies = base / "copies"
        cls.copies.mkdir()
        for name in piv.COPIES + ("control",):
            subprocess.run(["cp", "-a", str(cohort), str(cls.copies / name)], check=True)
            pipe = cls.copies / name / "data" / "fasta" / fifo.name
            pipe.unlink()
            pipe.write_bytes(texts[ids[64]])
        cls.original = {n: piv.tree(cls.copies / n / "data" / "species" / SLUG) for n in piv.COPIES}
        cls.partial = (piv.run_dir(cls.copies / "forge" / "data" / "species" / SLUG) / "block-00000.partial").read_bytes()
        cls.receipt = cls.copies.parent / "invalidation.json"
        cls.outcome = invalidate(cls.copies, cls.receipt)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_all_four_checks_are_refused_and_change_nothing(self):
        self.assertEqual(self.outcome.returncode, 0, self.outcome.stderr)
        got = json.loads(self.receipt.read_text())
        self.assertTrue(got["refused_as_stated"])
        self.assertEqual(sorted(got["checks"]), sorted(piv.COPIES))
        for name, row in got["checks"].items():
            with self.subTest(name):
                self.assertNotIn(row["run"]["rc"], (0, None))
                self.assertEqual(row["changed"], [])
        cmd = got["checks"]["width"]["run"]["cmd"]
        self.assertEqual(cmd[cmd.index("--block-width") + 1], "128")  # not the run's 64
        self.assertEqual(got["checks"]["forge"]["run"]["cmd"][-1], "--allow-missing-catalog")
        self.assertIn("unfinished disk-sparse run", got["checks"]["width"]["run"]["stderr"])
        self.assertIn("does not hold that genome's kset", got["checks"]["forge"]["run"]["stderr"])
        self.assertTrue(got["checks"]["link-matrix"]["target_kept"] and got["checks"]["link-kset"]["target_kept"])

    def test_the_forged_partial_differs_by_one_bit_and_its_receipt_agrees(self):
        folder = piv.run_dir(self.copies / "forge" / "data" / "species" / SLUG)
        forged = (folder / "block-00000.partial").read_bytes()
        self.assertEqual(len(forged), len(self.partial))
        diff = [(a ^ b) for a, b in zip(forged, self.partial) if a != b]
        self.assertEqual(len(diff), 1)
        self.assertEqual(bin(diff[0]).count("1"), 1)
        receipt = json.loads((folder / "block-00000.partial.receipt.json").read_text())
        self.assertEqual(receipt["sha256"], hashlib.sha256(forged).hexdigest())
        self.assertEqual(receipt["bytes"], len(forged))

    def test_only_the_mutation_differs_from_the_copy_as_taken(self):
        got = json.loads(self.receipt.read_text())
        for name in piv.COPIES:
            with self.subTest(name):
                now = piv.tree(self.copies / name / "data" / "species" / SLUG)
                changed = piv.changes(self.original[name], now)
                want = {"width": [], "forge": None, "link-matrix": ["matrix"], "link-kset": None}[name]
                if name == "forge":
                    self.assertTrue(changed and all(c.endswith((".partial", ".partial.receipt.json")) for c in changed), changed)
                elif name == "link-kset":
                    self.assertEqual(changed, [str(Path(got["checks"][name]["mutation"]["link"]).relative_to(self.copies / name / "data" / "species" / SLUG))])
                else:
                    self.assertEqual(changed, want)

    def test_a_resume_that_is_not_refused_fails_the_check(self):
        # The width check at the run's own width: GRM resumes and finishes, so the check must fail.
        args = type("A", (), {"slug": SLUG, "width": 64, "grm": tpi.grm(), "timeout": 300, "grm_arg": ["--allow-missing-catalog"]})()
        row = piv.check(self.copies / "control", args, None, 64, "unfinished disk-sparse run")
        self.assertEqual(row["run"]["rc"], 0)
        self.assertIn("was not refused", row["error"])
        self.assertTrue((self.copies / "control" / "data" / "species" / SLUG / "matrix" / "matrix.json").exists())
        # And a copy that is no longer the interrupted state is refused before any run.
        row = piv.check(self.copies / "control", args, None, 128, "unfinished disk-sparse run")
        self.assertIn("not the interrupted state", row["error"])
        self.assertNotIn("run", row)

    def test_missing_copies_or_an_existing_receipt_stop(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = invalidate(Path(tmp), Path(tmp) / "r.json")
            self.assertEqual(p.returncode, 1)
            self.assertIn("missing", json.loads((Path(tmp) / "r.json").read_text())["error"])
        before = self.receipt.read_bytes()
        with self.assertRaises(FileExistsError):
            piv.main(["--copies", str(self.copies), "--grm", tpi.grm(), "--width", "64", "--receipt", str(self.receipt),
                      "--grm-arg=--allow-missing-catalog"])
        self.assertEqual(self.receipt.read_bytes(), before)

    def test_drop_one_kmer_reads_the_partial_format(self):
        data = bytearray(self.partial)
        column, offset = piv.drop_one_kmer(data)
        self.assertGreaterEqual(column, 0)
        self.assertLess(column, 64)
        self.assertEqual(sum(a != b for a, b in zip(data, self.partial)), 1)
        with self.assertRaises(ValueError):
            piv.drop_one_kmer(bytearray(b"NOTAPART" + bytes(64)))


if __name__ == "__main__":
    unittest.main()
