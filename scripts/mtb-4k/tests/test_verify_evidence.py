#!/usr/bin/env python3
"""verify_evidence.py on a flat copy of the cohort-r2 evidence files, as push-data carries them into bindings/.

usage: python3 -I scripts/mtb-4k/tests/test_verify_evidence.py
"""
import contextlib, csv, hashlib, io, sys, tempfile, unittest
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import prepare_execution_roots as prep  # noqa: E402
import verify_evidence as ve  # noqa: E402


def sha(b):
    return hashlib.sha256(b).hexdigest()


class VerifyEvidence(unittest.TestCase):
    def setUp(self):
        self.saved = (dict(prep.EXPECTED), ve.EVIDENCE_SUMS_SHA256)
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        ids = [f"1773.{100000 + i}" for i in range(4000)]
        files = {"full": "".join(f"{g}\n" for g in ids).encode(), "prefix": "".join(f"{g}\n" for g in ids[:2000]).encode(),
                 "added": "".join(f"{g}\n" for g in ids[2000:]).encode()}
        buf = io.StringIO()
        w = csv.writer(buf, delimiter="\t", lineterminator="\n")
        w.writerow(["genome_id", "size_bytes", "sha256", "second_pass_size_bytes", "second_pass_sha256", "read_hash_reproducible"])
        for g in ids:
            w.writerow([g, 10, "a" * 64, 10, "a" * 64, "True"])
        files["fasta-manifest"] = buf.getvalue().encode()
        self.files = files
        self.write()

    def write(self):
        for key, data in self.files.items():
            (self.root / prep.EVIDENCE_FILES[key]).write_bytes(data)
        sums = "".join(f"{sha(self.files[k])}  ./{p}\n" for k, p in ve.EVIDENCE_PATHS.items()) + f"{'0' * 64}  ./logs/x.txt\n"
        (self.root / ve.EVIDENCE_SUMS).write_text(sums)
        ve.EVIDENCE_SUMS_SHA256 = sha(sums.encode())
        prep.EXPECTED.update({k: (sha(d), 4000 if k in ("full", "fasta-manifest") else 2000) for k, d in self.files.items()})

    def tearDown(self):
        expected, ve.EVIDENCE_SUMS_SHA256 = self.saved
        prep.EXPECTED.clear(); prep.EXPECTED.update(expected)
        self.tmp.cleanup()

    def run_main(self):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            return ve.main(["--root", str(self.root)]), out.getvalue(), err.getvalue()

    def test_the_flat_files_verify(self):
        rc, out, _ = self.run_main()
        self.assertEqual(rc, 0, out)

    def test_the_pins_are_the_evidences(self):
        self.assertEqual(self.saved[1], "d58008c92ccdc07fde911c6b0cf055d5e473c7416dc4887185c8500b82fef833")
        self.assertEqual(set(ve.EVIDENCE_PATHS), set(prep.EVIDENCE_FILES))

    def test_each_change_stops(self):
        cases = {
            "sha256sums_pinned": lambda: (self.root / ve.EVIDENCE_SUMS).write_text((self.root / ve.EVIDENCE_SUMS).read_text() + "\n"),
            "added_as_sha256sums_lists": lambda: (self.root / prep.EVIDENCE_FILES["added"]).write_bytes(self.files["added"] + b"x\n"),
            "full_as_preparation_pins": lambda: prep.EXPECTED.update({"full": ("0" * 64, 4000)}),
        }
        for check, change in cases.items():
            with self.subTest(check):
                self.write()
                change()
                rc, out, err = self.run_main()
                self.assertEqual(rc, 2)
                self.assertIn(check, out)
                self.assertIn("verify_evidence: STOP", err)

    def test_the_cohort_checks_stop(self):
        self.files["prefix"] = self.files["added"]
        self.write()
        rc, out, _ = self.run_main()
        self.assertEqual(rc, 2)
        self.assertIn("frozen_2k_exact_prefix", out)

    def test_a_missing_file_is_a_stop_line(self):
        (self.root / prep.EVIDENCE_FILES["fasta-manifest"]).unlink()
        rc, out, err = self.run_main()
        self.assertEqual(rc, 2)
        self.assertIn("cannot read the evidence", err)
        self.assertNotIn("Traceback", err)


if __name__ == "__main__":
    unittest.main(verbosity=2)
