"""scripts/c4-screen/run_checks.py: the checks of the run's copies, on SYNTHETIC copies.

Every species, dataset, result file, genome id and binary here is made up in the test (tests/c4_screen/page_fixture.py).
"""

import builtins
import contextlib
import io
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests.c4_screen import common
from tests.c4_screen import page_fixture as fx

checks = common.load("run_checks")
NAMES = ["the first run's species, each once",
         "each species ended DONE, any FAILED older than its DONE",
         "each species' SHA256SUMS holds every file of the run's rule, each at its sha256",
         "score/score.log stays out of SHA256SUMS",
         "each DONE names its species, its binary.txt's binary and code, its results and its SHA256SUMS",
         "every species ran one commit and one code sha256",
         "each lists.txt names its species' first-run datasets, seeds 42, 43 and 44, twin 45, 46 and 47 and the three "
         "arms",
         "each expect.tsv is expect.py's list from lists.tsv",
         "each species holds every listed result or refusal record and inputs file, and no other",
         "each refusal record is a study run scored as the constant model on five learners, and each constant-model "
         "row a refusal record's",
         "no genome id in any file that will be read or published"]


def reseal(out, slug):
    """SHA256SUMS and DONE made again after a change, as run_species.sh would have made them."""
    fx.sums(out)
    fx.done(out, slug, len(list((out / "results" / "picks").rglob("*.pick.json"))))


class Copies(unittest.TestCase):
    build = {}

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.root = fx.build_copies(Path(cls.tmp.name) / "base", **cls.build)
        cls.n = 0

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def fresh(self):
        type(self).n += 1
        return fx.copy(self.root, Path(self.tmp.name) / f"copy{self.n}")

    def run_checks(self, root):
        lines, aside = checks.check(root / "copies", root / "lists.tsv", root / "release")
        return {name: (passed, detail) for passed, name, detail in lines}, aside

    def failed(self, root):
        got, _ = self.run_checks(root)
        return [name for name, (passed, _) in got.items() if not passed]


class Check(Copies):
    def test_clean_copies_pass_every_check(self):
        root = self.fresh()
        got, aside = self.run_checks(root)
        self.assertEqual(list(got), NAMES)
        self.assertEqual([n for n, (p, _) in got.items() if not p], [])
        self.assertEqual(aside, [])
        self.assertEqual(got[NAMES[-2]][1], "1 refusal records, 5 constant-model rows")
        self.assertTrue(got[NAMES[-1]][1].startswith("0 genome ids in "))
        with contextlib.redirect_stdout(io.StringIO()) as out:
            code = checks.main(["check", str(root / "copies"), str(root / "lists.tsv"), str(root / "release"),
                                str(root / "checks.txt")])
        self.assertEqual(code, 0)
        self.assertEqual(out.getvalue().splitlines()[-1], "CHECKS 11 of 11")
        self.assertEqual((root / "checks.txt").read_text(encoding="utf-8"), out.getvalue())

    def test_a_species_folder_more_is_refused(self):
        root = self.fresh()
        shutil.copytree(root / "copies" / "synthia_alpha", root / "copies" / "synthia_delta")
        self.assertEqual(self.failed(root), [NAMES[0]])

    def test_a_species_with_no_done_is_refused(self):
        root = self.fresh()
        (root / "copies" / "synthia_beta" / "DONE").unlink()
        self.assertIn(NAMES[1], self.failed(root))

    def test_a_failed_older_than_done_is_named_and_set_aside(self):
        root = self.fresh()
        (root / "copies" / "synthia_beta" / "FAILED").write_text(
            "FAILED synthia_beta stopped fit 2026-11-23T12:00:00Z\n", encoding="utf-8")
        got, aside = self.run_checks(root)
        self.assertEqual([n for n, (p, _) in got.items() if not p], [])
        self.assertEqual(aside, ["synthia_beta/FAILED at 2026-11-23T12:00:00Z, older than its DONE at "
                                 "2026-11-23T13:00:00Z"])
        lines = checks.report(*checks.check(root / "copies", root / "lists.tsv", root / "release"))
        self.assertIn("SET ASIDE synthia_beta/FAILED at 2026-11-23T12:00:00Z", "\n".join(lines))
        self.assertEqual(lines[-1], "CHECKS 11 of 11")

    def test_a_failed_newer_than_done_is_refused(self):
        root = self.fresh()
        (root / "copies" / "synthia_beta" / "FAILED").write_text(
            "FAILED synthia_beta stopped fit 2026-11-23T13:30:00Z\n", encoding="utf-8")
        self.assertEqual(self.failed(root), [NAMES[1]])

    def test_a_file_that_differs_from_its_sha256_is_refused(self):
        root = self.fresh()
        with open(root / "copies" / "synthia_gamma" / "check" / "refit.tsv", "a", encoding="utf-8") as f:
            f.write("changed\t0\n")
        self.assertEqual(self.failed(root), [NAMES[2]])

    def test_score_log_in_sha256sums_is_refused(self):
        root = self.fresh()
        out = root / "copies" / "synthia_alpha"
        with open(out / "SHA256SUMS", "a", encoding="utf-8") as f:
            f.write(f"{checks.sha256(out / 'score' / 'score.log')}  score/score.log\n")
        fx.done(out, "synthia_alpha", len(list((out / "results" / "picks").rglob("*.pick.json"))))
        self.assertIn(NAMES[3], self.failed(root))

    def test_a_done_naming_other_sha256sums_is_refused(self):
        root = self.fresh()
        path = root / "copies" / "synthia_alpha" / "DONE"
        words = path.read_text(encoding="utf-8").split()
        path.write_text(" ".join(words[:-1] + ["0" * 64]) + "\n", encoding="utf-8")
        self.assertEqual(self.failed(root), [NAMES[4]])

    def test_an_expect_tsv_that_is_not_expect_py_s_is_refused(self):
        root = self.fresh()
        out = root / "copies" / "synthia_alpha"
        lines = (out / "expect.tsv").read_text(encoding="utf-8").splitlines()
        (out / "expect.tsv").write_text("\n".join(lines[:-1]) + "\n", encoding="utf-8")
        reseal(out, "synthia_alpha")
        self.assertEqual(self.failed(root), [NAMES[7]])

    def test_a_missing_result_is_refused(self):
        root = self.fresh()
        out = root / "copies" / "synthia_beta"
        next((out / "results" / "picks" / "twin").glob("*.pick.json")).unlink()
        reseal(out, "synthia_beta")
        self.assertEqual(self.failed(root), [NAMES[8]])

    def test_a_result_more_than_listed_is_refused(self):
        root = self.fresh()
        out = root / "copies" / "synthia_beta"
        (out / "results" / "picks" / "twin" / "synthia_beta__synthetic_drug_1.lineage-s48.pick.json").write_text(
            "{}\n", encoding="utf-8")
        reseal(out, "synthia_beta")
        self.assertEqual(self.failed(root), [NAMES[8]])

    def test_a_refusal_record_with_no_constant_model_rows_is_refused(self):
        root = self.fresh()
        out = root / "copies" / "synthia_gamma"
        path = out / "score" / "draws.tsv"
        rows = fx.read(path)
        for r in rows:
            if r["constant"] == "1":
                r.update(constant="0", why="")
        fx.rewrite(path, rows, fx.header(path))
        reseal(out, "synthia_gamma")
        self.assertEqual(self.failed(root), [NAMES[9]])

    def test_a_genome_id_in_a_table_is_refused(self):
        root = self.fresh()
        out = root / "copies" / "synthia_alpha"
        with open(out / "check" / "refit.tsv", "a", encoding="utf-8") as f:
            f.write(f"{fx.GENOME_IDS[0]}\t1\n")
        reseal(out, "synthia_alpha")
        got, _ = self.run_checks(root)
        self.assertEqual([n for n, (p, _) in got.items() if not p], [NAMES[10]])
        self.assertTrue(got[NAMES[10]][1].startswith("1 genome ids in "))

    def test_score_log_is_never_opened(self):
        root = self.fresh()
        (root / "copies" / "synthia_alpha" / "score" / "score.log").write_text(f"{fx.GENOME_IDS[0]}\n",
                                                                                encoding="utf-8")
        self.assertEqual(self.failed(root), [])


class Unopened(Copies):
    """score/score.log is never opened, even when a SHA256SUMS names it (SYNTHETIC copies)."""

    def test_score_log_named_in_sha256sums_is_refused_and_never_opened(self):
        root = self.fresh()
        out = root / "copies" / "synthia_alpha"
        with open(out / "SHA256SUMS", "a", encoding="utf-8") as f:
            f.write(f"{'0' * 64}  score/score.log\n")
        fx.done(out, "synthia_alpha", len(list((out / "results" / "picks").rglob("*.pick.json"))))
        opened = []
        real_open, real_io_open = builtins.open, io.open

        def guard(real):
            def wrapped(file, *args, **kwargs):
                if str(file).endswith("score/score.log"):
                    opened.append(str(file))
                return real(file, *args, **kwargs)
            return wrapped

        with mock.patch("builtins.open", guard(real_open)), mock.patch("io.open", guard(real_io_open)):
            failed = self.failed(root)
        self.assertIn(NAMES[3], failed)
        self.assertEqual(opened, [])


class Timeless(Copies):
    def test_a_failed_with_no_time_is_refused_and_not_set_aside(self):
        for text in ("", "FAILED synthia_beta stopped fit\n", "FAILED synthia_beta stopped fit 12:00\n"):
            root = self.fresh()
            (root / "copies" / "synthia_beta" / "FAILED").write_text(text, encoding="utf-8")
            got, aside = self.run_checks(root)
            self.assertEqual([n for n, (p, _) in got.items() if not p], [NAMES[1]], text)
            self.assertEqual(aside, [], text)


class Prefixed(Copies):
    def test_a_genome_id_after_a_letter_or_underscore_is_refused(self):
        for text in (f"genome_{fx.GENOME_IDS[0]}", f"id{fx.GENOME_IDS[1]}"):
            root = self.fresh()
            out = root / "copies" / "synthia_alpha"
            with open(out / "check" / "refit.tsv", "a", encoding="utf-8") as f:
                f.write(f"{text}\t1\n")
            reseal(out, "synthia_alpha")
            got, _ = self.run_checks(root)
            self.assertEqual([n for n, (p, _) in got.items() if not p], [NAMES[10]], text)
            self.assertTrue(got[NAMES[10]][1].startswith("1 genome ids in "), text)

    def test_a_longer_number_holding_a_genome_id_is_not_one(self):
        root = self.fresh()
        out = root / "copies" / "synthia_alpha"
        with open(out / "check" / "refit.tsv", "a", encoding="utf-8") as f:
            f.write(f"1{fx.GENOME_IDS[0]}\t0.{fx.GENOME_IDS[0].replace('.', '')}\n")
        reseal(out, "synthia_alpha")
        self.assertEqual(self.failed(root), [])


class TwoBinaries(Copies):
    build = {"codes": {"synthia_beta": "d" * 64}}

    def test_two_code_sha256_are_refused(self):
        root = self.fresh()
        got, _ = self.run_checks(root)
        self.assertEqual([n for n, (p, _) in got.items() if not p], [NAMES[5]])
        self.assertEqual(got[NAMES[5]][1], "1 commits, 2 code sha256, 1 binaries")


class Gather(Copies):
    def test_gather_copies_five_files_a_species_and_never_score_log(self):
        root = self.fresh()
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(checks.main(["check", str(root / "copies"), str(root / "lists.tsv"),
                                          str(root / "release"), str(root / "checks.txt")]), 0)
        n = checks.gather(root / "copies", root / "lists.tsv", root / "checks.txt", root / "reading")
        self.assertEqual(n, 15)
        got = sorted(p.relative_to(root / "reading").as_posix() for p in (root / "reading").rglob("*") if p.is_file())
        self.assertEqual(got, sorted(f"score/{s}/{f}" for s, _ in fx.SPECIES
                                     for f in ("checks.tsv", "draws.tsv", "datasets.tsv", "fits.tsv", "binary.txt")))
        self.assertEqual((root / "reading" / "score" / "synthia_alpha" / "fits.tsv").read_bytes(),
                         (root / "copies" / "synthia_alpha" / "logs" / "fits.tsv").read_bytes())

    def test_gather_refuses_checks_that_did_not_pass(self):
        root = self.fresh()
        (root / "checks.txt").write_text("FAIL x\nCHECKS 10 of 11\n", encoding="utf-8")
        with self.assertRaises(SystemExit):
            checks.gather(root / "copies", root / "lists.tsv", root / "checks.txt", root / "reading")
        self.assertFalse((root / "reading").exists())


if __name__ == "__main__":
    unittest.main()
