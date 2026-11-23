"""scripts/c5-settings/run_checks.py: the checks of the run's copies before any score is read, on SYNTHETIC copies.

Every copy, genome id and file here is made up (page_fixture.py).
"""

import builtins
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests.c5_settings import common
from tests.c5_settings import page_fixture as fx

checks = common.load("run_checks")


class Case(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.top = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def run_checks(self, change=None, after=None):
        copies, lists, release = fx.copies(self.top, change)
        if after:
            after(copies)
        lines, aside = checks.check(copies, lists, release)
        return checks.report(lines, aside)

    def line(self, report, start):
        found = [x for x in report if x.split(" ", 1)[1].startswith(start)]
        self.assertEqual(len(found), 1, start)
        return found[0]


class Check(Case):
    def test_clean_copies_pass_every_check(self):
        report = self.run_checks()
        self.assertEqual(report[-1], "CHECKS 11 of 11")
        self.assertFalse([x for x in report if x.startswith("FAIL") or x.startswith("SET ASIDE")])
        self.assertIn("0 refusal records, 0 constant-model rows", self.line(report, "the refusal records"))

    def test_an_archive_beside_the_species_folders_is_not_read(self):
        report = self.run_checks(after=lambda copies: (copies / "back.tar").write_bytes(b"SYNTHETIC"))
        self.assertEqual(report[-1], "CHECKS 11 of 11")

    def test_a_species_more_than_the_first_runs_is_refused(self):
        report = self.run_checks(after=lambda copies: (copies / "synthia_delta").mkdir())
        self.assertTrue(self.line(report, "the first run's species").startswith("FAIL"))

    def test_a_failed_file_newer_than_done_is_refused(self):
        def after(copies):
            (copies / "synthia_beta" / "FAILED").write_text("FAILED synthia_beta 1 score 2026-11-23T13:00:00Z\n",
                                                            encoding="utf-8")
        report = self.run_checks(after=after)
        self.assertTrue(self.line(report, "each species ended DONE").startswith("FAIL"))

    def test_a_failed_file_older_than_done_is_named_and_set_aside(self):
        def after(copies):
            (copies / "synthia_beta" / "FAILED").write_text("FAILED synthia_beta 1 fit 2026-11-23T11:00:00Z\n",
                                                            encoding="utf-8")
        report = self.run_checks(after=after)
        self.assertIn("SET ASIDE synthia_beta/FAILED at 2026-11-23T11:00:00Z, older than its DONE at "
                      "2026-11-23T12:00:00Z", report)
        self.assertEqual(report[-1], "CHECKS 11 of 11")

    def test_a_failed_file_with_no_time_is_refused_and_not_set_aside(self):
        for i, text in enumerate(("", "FAILED synthia_beta 1 score\n", "FAILED synthia_beta 1 score 12:00\n")):
            self.top = Path(self.tmp.name) / f"timeless{i}"
            self.top.mkdir()

            def after(copies, text=text):
                (copies / "synthia_beta" / "FAILED").write_text(text, encoding="utf-8")
            report = self.run_checks(after=after)
            self.assertTrue(self.line(report, "each species ended DONE").startswith("FAIL"), text)
            self.assertFalse([x for x in report if x.startswith("SET ASIDE")], text)

    def test_a_file_changed_after_its_sha256sums_is_refused(self):
        def after(copies):
            (copies / "synthia_gamma" / "check" / "refit.tsv").write_text("check\tpassed\nSYNTHETIC\t0\n",
                                                                         encoding="utf-8")
        report = self.run_checks(after=after)
        self.assertIn("1 files differ from SHA256SUMS", self.line(report, "each species' SHA256SUMS"))

    def test_score_log_in_sha256sums_is_refused(self):
        def after(copies):
            path = copies / "synthia_alpha" / "SHA256SUMS"
            path.write_text(path.read_text(encoding="utf-8") + f"{'0' * 64}  score/score.log\n", encoding="utf-8")
        report = self.run_checks(after=after)
        self.assertTrue(self.line(report, "score/score.log stays out").startswith("FAIL"))

    def test_two_code_sha256_are_refused(self):
        def change(slug, out):
            if slug == "synthia_gamma":
                (out / "binary.txt").write_text(f"commit {fx.COMMIT}\ngrm {fx.GRM}\ncode {'d' * 64}\n",
                                                encoding="utf-8")
        report = self.run_checks(change)
        self.assertTrue(self.line(report, "every species ran one commit").startswith("FAIL"))
        self.assertTrue(self.line(report, "each DONE names").startswith("PASS"))

    def test_an_expect_tsv_that_is_not_expect_pys_is_refused(self):
        def change(slug, out):
            if slug == "synthia_beta":
                lines = (out / "expect.tsv").read_text(encoding="utf-8").splitlines()
                (out / "expect.tsv").write_text("\n".join(lines[:-1]) + "\n", encoding="utf-8")
        report = self.run_checks(change)
        self.assertTrue(self.line(report, "each expect.tsv").startswith("FAIL"))

    def test_a_missing_choices_file_is_refused(self):
        def change(slug, out):
            if slug == "synthia_beta":
                (out / "results" / "picks" / "choices"
                 / "synthia_beta__synthetic_drug_3.study-whole-s42.f3.choices.json").unlink()
        report = self.run_checks(change)
        self.assertIn("1 missing, 0 more", self.line(report, "each species holds every default"))

    def test_a_genome_id_in_a_choices_file_is_refused(self):
        def change(slug, out):
            if slug == "synthia_alpha":
                p = out / "results" / "picks" / "choices" / "synthia_alpha__synthetic_drug_1.lineage-s42.f0.choices.json"
                doc = json.loads(p.read_text(encoding="utf-8"))
                doc["tuning"] = fx.genome("synthia_alpha", 2)
                p.write_text(json.dumps(doc), encoding="utf-8")
        report = self.run_checks(change)
        line = self.line(report, "no genome id")
        self.assertTrue(line.startswith("FAIL"))
        self.assertIn(": 1 genome ids", line)

    def test_a_genome_id_after_a_letter_or_underscore_is_refused(self):
        for i, text in enumerate((f"genome_{fx.genome('synthia_alpha', 2)}", f"id{fx.genome('synthia_alpha', 3)}")):
            self.top = Path(self.tmp.name) / f"prefixed{i}"
            self.top.mkdir()

            def change(slug, out, text=text):
                if slug == "synthia_alpha":
                    with open(out / "check" / "refit.tsv", "a", encoding="utf-8") as f:
                        f.write(f"{text}\t1\n")
            line = self.line(self.run_checks(change), "no genome id")
            self.assertTrue(line.startswith("FAIL"), text)
            self.assertIn(": 1 genome ids", line)

    def test_a_longer_number_holding_a_genome_id_is_not_one(self):
        gid = fx.genome("synthia_alpha", 2)

        def change(slug, out):
            if slug == "synthia_alpha":
                with open(out / "check" / "refit.tsv", "a", encoding="utf-8") as f:
                    f.write(f"1{gid}\t0.{gid.replace('.', '')}\n")
        self.assertTrue(self.line(self.run_checks(change), "no genome id").startswith("PASS"))

    def test_a_refusal_record_with_no_constant_model_rows_is_refused(self):
        def change(slug, out):
            if slug == "synthia_gamma":
                (out / "results" / "picks" / "default"
                 / "synthia_gamma__synthetic_drug_4.study-whole-s42.refused.json").write_text("{}", encoding="utf-8")
        report = self.run_checks(change)
        line = self.line(report, "the refusal records")
        self.assertTrue(line.startswith("FAIL"))
        self.assertIn("1 refusal records, 0 constant-model rows", line)


class Unopened(Case):
    """score/score.log is never opened, even when a SHA256SUMS names it (SYNTHETIC copies)."""

    def test_score_log_named_in_sha256sums_is_refused_and_never_opened(self):
        def after(copies):
            path = copies / "synthia_alpha" / "SHA256SUMS"
            path.write_text(path.read_text(encoding="utf-8") + f"{'0' * 64}  score/score.log\n", encoding="utf-8")

        opened = []
        real_open, real_io_open = builtins.open, io.open

        def guard(real):
            def wrapped(file, *args, **kwargs):
                if str(file).endswith("score/score.log"):
                    opened.append(str(file))
                return real(file, *args, **kwargs)
            return wrapped

        copies, lists, release = fx.copies(self.top)
        after(copies)
        with mock.patch("builtins.open", guard(real_open)), mock.patch("io.open", guard(real_io_open)):
            lines, aside = checks.check(copies, lists, release)
        report = checks.report(lines, aside)
        self.assertTrue(self.line(report, "each species' SHA256SUMS").startswith("FAIL"))
        self.assertTrue(self.line(report, "score/score.log stays out").startswith("FAIL"))
        self.assertEqual(opened, [])


class Gather(Case):
    def test_gather_copies_four_files_a_species_and_the_choices_and_never_score_log(self):
        copies, lists, release = fx.copies(self.top)
        lines, aside = checks.check(copies, lists, release)
        path = self.top / "checks.txt"
        path.write_text("\n".join(checks.report(lines, aside)) + "\n", encoding="utf-8")
        n = checks.gather(copies, lists, path, self.top / "reading")
        # 4 datasets: 3 lineage draws and 1 study run each, 5 folds; five files a species with folds.tsv.
        self.assertEqual(n, 3 * 5 + 4 * 4 * 5)
        folds = [r for p in sorted((self.top / "reading" / "score").glob("*/folds.tsv")) for r in fx.C3.read_tsv(p)]
        self.assertEqual(len(folds), 4 * 4)
        self.assertEqual({r["folds"] for r in folds}, {str(fx.FOLDS)})
        got = sorted(p.relative_to(self.top / "reading").as_posix() for p in (self.top / "reading").rglob("*")
                     if p.is_file())
        self.assertFalse([g for g in got if g.endswith("score.log")])
        self.assertIn("score/synthia_beta/binary.txt", got)

    def test_gather_refuses_checks_that_did_not_pass(self):
        copies, lists, _ = fx.copies(self.top)
        path = self.top / "checks.txt"
        path.write_text("FAIL x\nCHECKS 0 of 1\n", encoding="utf-8")
        with self.assertRaises(SystemExit):
            checks.gather(copies, lists, path, self.top / "reading")


if __name__ == "__main__":
    unittest.main()
