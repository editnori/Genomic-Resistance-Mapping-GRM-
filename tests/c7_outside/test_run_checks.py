"""scripts/c7-outside/run_checks.py on SYNTHETIC copies: each check fails on what it names, a first copy set aside
and a FAILED older than DONE are named and set aside by time, and gather writes the reading only when every check
passes. Every species, dataset, genome id and file here is made up (page_fixture.py)."""

import builtins
import io
import shutil
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests.c7_outside import page_fixture as F

checks = F.checks


class RunChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        F.build(cls.root / "run")
        F.build(cls.root / "withdrawn", withdrawn=True)
        cls.n = 0

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def back(self, which="run"):
        type(self).n += 1
        dst = self.root / f"copy-{self.n}"
        shutil.copytree(self.root / which / "back", dst / "back")
        return dst / "back"

    def run_on(self, back, which="run", withdrawn=None):
        return checks.run(back, self.root / which / "lists.tsv", self.root / which / "release", withdrawn)

    def lines(self, rep):
        return "\n".join(rep.lines)

    def test_the_synthetic_copies_pass_every_check(self):
        rep = self.run_on(self.back())
        self.assertEqual(rep.failed, 0, self.lines(rep))
        self.assertEqual(rep.lines[-1], "CHECKS 6 of 6")
        rep = self.run_on(self.back("withdrawn"), "withdrawn", F.WITHDRAWN)
        self.assertEqual(rep.failed, 0, self.lines(rep))

    def test_a_first_copy_older_than_its_done_is_named_and_set_aside(self):
        back = self.back()
        shutil.copytree(back / "synthia_coli", back / "synthia_coli.20261123T080000Z")
        rep = self.run_on(back)
        self.assertEqual(rep.failed, 0, self.lines(rep))
        self.assertIn("synthia_coli.20261123T080000Z: a first copy of synthia_coli, older than its DONE",
                      self.lines(rep))

    def test_a_first_copy_newer_than_its_done_fails(self):
        back = self.back()
        shutil.copytree(back / "synthia_coli", back / "synthia_coli.20261123T120000Z")
        self.assertIn("FAIL species", self.lines(self.run_on(back)))

    def test_a_copy_of_a_species_outside_the_run_fails(self):
        back = self.back()
        shutil.copytree(back / "synthia_coli", back / "nullia_species")
        self.assertIn("FAIL species", self.lines(self.run_on(back)))

    def test_a_sha256sums_other_than_the_one_done_names_fails(self):
        back = self.back()
        path = back / "synthia_coli" / "SHA256SUMS"
        path.write_text("".join(reversed(path.read_text().splitlines(keepends=True))))
        text = self.lines(self.run_on(back))
        self.assertIn("FAIL SHA256SUMS", text)
        self.assertIn("lines hold", text)

    def test_a_lists_txt_whose_arms_are_not_its_done_pass_fails(self):
        back = self.back("withdrawn")
        done = back / "fictus_aureus" / "DONE"
        done.write_text(done.read_text().replace("elasticnet leave-out", "elasticnet run"))
        self.assertIn("FAIL expect.py", self.lines(self.run_on(back, "withdrawn", F.WITHDRAWN)))

    def test_a_species_with_no_release_list_fails(self):
        back = self.back()
        release = back.parent / "release"
        shutil.copytree(self.root / "run" / "release", release)
        (release / "species" / "synthia_coli.json").unlink()
        rep = checks.run(back, self.root / "run" / "lists.tsv", release)
        self.assertIn("FAIL genome ids", self.lines(rep))

    def test_a_failed_older_than_its_done_is_set_aside(self):
        back = self.back()
        (back / "fictus_aureus" / "FAILED").write_text("FAILED fictus_aureus at fits 2026-11-23T09:00:00Z\n")
        rep = self.run_on(back)
        self.assertEqual(rep.failed, 0, self.lines(rep))
        self.assertIn("fictus_aureus: FAILED at fits at 2026-11-23T09:00:00Z, before DONE", self.lines(rep))

    def test_a_failed_newer_than_its_done_fails(self):
        back = self.back()
        (back / "fictus_aureus" / "FAILED").write_text("FAILED fictus_aureus stopped at fits 2026-11-23T11:00:00Z\n")
        self.assertIn("FAIL DONE", self.lines(self.run_on(back)))

    def test_a_failed_whose_time_is_no_utc_second_fails(self):
        back = self.back()
        (back / "fictus_aureus" / "FAILED").write_text("FAILED fictus_aureus at fits 08:00\n")
        self.assertIn("FAIL DONE", self.lines(self.run_on(back)))

    def test_a_done_whose_time_is_no_utc_second_fails(self):
        back = self.back()
        shutil.copytree(back / "synthia_coli", back / "synthia_coli.20261123T080000Z")
        done = back / "synthia_coli" / "DONE"
        done.write_text(done.read_text().replace(F.DONE_AT, "9999", 1))
        text = self.lines(self.run_on(back))
        self.assertIn("FAIL species", text)
        self.assertIn("FAIL DONE", text)

    def test_a_changed_result_fails_its_sha256sums(self):
        back = self.back()
        path = back / "synthia_coli" / "results" / "outside" / "xgboost" / "synthia_coli__drug_b.lineage-s43.outside.json"
        path.write_text(path.read_text().replace('"rounds_kept": 35', '"rounds_kept": 34', 1))
        self.assertIn("FAIL SHA256SUMS", self.lines(self.run_on(back)))

    def test_a_first_pass_beside_second_passes_fails_done(self):
        back = self.back()
        done = back / "madeup_pneumoniae" / "DONE"
        done.write_text(done.read_text().replace("elasticnet run", "elasticnet leave-out"))
        text = self.lines(self.run_on(back))
        self.assertIn("FAIL DONE", text)
        self.assertIn("elastic net leave-out: madeup_pneumoniae", text)

    def test_the_first_pass_everywhere_fails_unless_the_elastic_net_is_withdrawn(self):
        text = self.lines(self.run_on(self.back("withdrawn"), "withdrawn"))
        self.assertIn("FAIL DONE", text)
        self.assertIn("PASS expect.py", text)

    def test_two_code_sha256_fail_done(self):
        back = self.back()
        done = back / "fictus_aureus" / "DONE"
        done.write_text(done.read_text().replace(F.CODE, F.sha(b"SYNTHETIC other")))
        self.assertIn("FAIL DONE", self.lines(self.run_on(back)))

    def test_an_expect_tsv_that_is_not_expect_pys_fails(self):
        back = self.back()
        path = back / "synthia_coli" / "expect.tsv"
        path.write_text("".join(path.read_text().splitlines(keepends=True)[:-1]))
        F.resum(back / "synthia_coli")
        rep = self.run_on(back)
        self.assertEqual([line for line in rep.lines if line.startswith("FAIL")][0][:15], "FAIL expect.py:")
        self.assertEqual(rep.failed, 1, self.lines(rep))

    def test_a_report_expect_tsv_does_not_list_fails(self):
        back = self.back()
        folder = back / "synthia_coli" / "results" / "outside" / "lightgbm"
        shutil.copyfile(folder / "synthia_coli__drug_b.lineage-s42.outside.json",
                        folder / "synthia_coli__drug_b.lineage-s45.outside.json")
        self.assertIn("FAIL expect.py", self.lines(self.run_on(back)))

    def test_a_genome_id_in_a_log_fails(self):
        back = self.back()
        (back / "synthia_coli" / "logs" / "run.log").write_text(f"SYNTHETIC: {F.genome_ids('synthia_coli')[1]}\n")
        self.assertIn("FAIL genome ids", self.lines(self.run_on(back)))

    def test_a_genome_id_joined_to_a_word_fails(self):
        for text in (f"genome_{F.genome_ids('synthia_coli')[1]}", f"id{F.genome_ids('synthia_coli')[2]}"):
            back = self.back()
            (back / "synthia_coli" / "logs" / "run.log").write_text(f"SYNTHETIC: {text}\n")
            self.assertIn("FAIL genome ids", self.lines(self.run_on(back)), text)

    def test_a_column_name_in_a_log_fails(self):
        back = self.back()
        (back / "fictus_aureus" / "logs" / "run.log").write_text("SYNTHETIC: Column_12\n")
        self.assertIn("FAIL features", self.lines(self.run_on(back)))

    def test_a_unitig_in_a_log_fails(self):
        back = self.back()
        (back / "fictus_aureus" / "logs" / "run.log").write_text("SYNTHETIC: " + "ACGT" * 8 + "\n")
        self.assertIn("FAIL features", self.lines(self.run_on(back)))

    def test_gather_writes_nothing_when_a_check_fails(self):
        back = self.back()
        (back / "fictus_aureus" / "logs" / "run.log").write_text("SYNTHETIC: Column_12\n")
        reading = back.parent / "reading"
        rep = checks.gather(back, self.root / "run" / "lists.tsv", self.root / "run" / "release", reading)
        self.assertGreater(rep.failed, 0)
        self.assertFalse(reading.exists())

    def test_gather_writes_the_reading_page_tables_reads(self):
        reading = self.root / "run" / "reading"
        self.assertEqual(sorted(p.name for p in reading.iterdir()),
                         ["checks.txt", "count", "reading.txt", "run", "score", "withdrawn.tsv"])
        self.assertEqual(sorted(p.name for p in (reading / "run" / "synthia_coli").iterdir()),
                         ["DONE", "libraries.txt", "outside", "plans.log", "timing.tsv"])
        self.assertEqual((reading / "checks.txt").read_text().splitlines()[-1], "CHECKS 6 of 6")
        self.assertEqual((self.root / "withdrawn" / "reading" / "withdrawn.tsv").read_text(),
                         "arm\treason\ttime\nelasticnet\t" + "\t".join(F.WITHDRAWN["elasticnet"]) + "\n")

    def test_gather_refuses_a_withdrawal_time_later_than_its_clock(self):
        back = self.back("withdrawn")
        reading = back.parent / "reading"
        with self.assertRaisesRegex(SystemExit, "later than this clock"):
            checks.gather(back, self.root / "withdrawn" / "lists.tsv", self.root / "withdrawn" / "release", reading,
                          {"elasticnet": ("SYNTHETIC: cost", "2099-01-01T00:00Z")})
        self.assertFalse(reading.exists())
        self.assertTrue(checks.later_than_now("2026-11-23T09:00:01Z", checks.datetime(2026, 11, 23, 9, 0,
                                                                                       tzinfo=checks.timezone.utc)))
        self.assertFalse(checks.later_than_now("2026-11-23T09:00Z", checks.datetime(2026, 11, 23, 9, 0,
                                                                                     tzinfo=checks.timezone.utc)))

    def test_a_withdrawn_arm_needs_a_reason_and_a_utc_time(self):
        for argv in (["--withdrawn", "elasticnet", "", "2026-11-23T09:00Z"],
                     ["--withdrawn", "elasticnet", "cost", "23 Nov"], ["--withdrawn", "twin", "cost", "2026-11-23T09:00Z"],
                     ["--withdrawn", "elasticnet", "cost"]):
            with self.assertRaises(SystemExit):
                checks.parse(argv)
        self.assertEqual(checks.parse(["--withdrawn", "elasticnet", "cost", "2026-11-23T09:00Z"])[1],
                         {"elasticnet": ("cost", "2026-11-23T09:00Z")})

    def second_passes(self, back):
        """The elastic net withdrawn after two species began a second pass: synthia_coli's second pass came back with
        its own DONE, and madeup_pneumoniae's is still running, so it has no folder. Each second-pass copy holds a
        genome id and a matrix column name, so that reading one would fail a check."""
        for slug in ("synthia_coli", "madeup_pneumoniae"):
            (back / slug).rename(back / f"{slug}.20261123T120000Z")
        shutil.copytree(back / "synthia_coli.20261123T120000Z", back / "synthia_coli")
        done = back / "synthia_coli" / "DONE"
        done.write_text(done.read_text().replace(F.DONE_AT, "2026-11-23T15:00:00Z")
                        .replace("elasticnet leave-out", "elasticnet run"))
        (back / "synthia_coli" / "logs" / "run.log").write_text(
            f"SYNTHETIC: {F.genome_ids('synthia_coli')[2]} Column_7\n")
        return back / "synthia_coli"

    def test_late_withdrawal_reads_each_first_copy_and_never_opens_a_second_pass(self):
        back = self.back("withdrawn")
        second = self.second_passes(back)
        opened, _open, _io_open, _tar_open = [], builtins.open, io.open, tarfile.open

        def watch(real):
            def opener(name, *a, **k):
                opened.append(str(name))
                return real(name, *a, **k)
            return opener

        with mock.patch.object(builtins, "open", watch(_open)), mock.patch.object(io, "open", watch(_io_open)), \
                mock.patch.object(tarfile, "open", watch(_tar_open)):
            rep = self.run_on(back, "withdrawn", F.WITHDRAWN)
        text = self.lines(rep)
        self.assertEqual(rep.failed, 0, text)
        self.assertIn("the elastic net withdrawn: 2 species ran a second pass, 1 second-pass copies set aside", text)
        self.assertIn("synthia_coli: the second pass of synthia_coli, after the elastic net was withdrawn", text)
        self.assertIn("madeup_pneumoniae.20261123T120000Z: the copy of madeup_pneumoniae", text)
        self.assertEqual([x for x in opened if x == str(second) or x.startswith(f"{second}/")], [])
        self.assertTrue([x for x in opened if "synthia_coli.20261123T120000Z" in x])
        reading = back.parent / "reading"
        rep = checks.gather(back, self.root / "withdrawn" / "lists.tsv", self.root / "withdrawn" / "release", reading,
                            F.WITHDRAWN)
        self.assertEqual(rep.failed, 0, self.lines(rep))
        self.assertIn("elasticnet leave-out", (reading / "run" / "synthia_coli" / "DONE").read_text())
        self.assertEqual(sorted(p.name for p in (reading / "run").iterdir()), sorted(s for s, _, _ in F.SPECIES))

    def test_late_withdrawal_takes_no_copy_whose_done_says_the_net_ran(self):
        back = self.back()
        (back / "synthia_coli").rename(back / "synthia_coli.20261123T080000Z")
        shutil.copytree(back / "synthia_coli.20261123T080000Z", back / "synthia_coli")
        text = self.lines(self.run_on(back, "run", F.WITHDRAWN))
        self.assertIn("FAIL DONE", text)
        self.assertIn("the elastic net withdrawn: 0 species ran a second pass", text)


class ScoreLog(unittest.TestCase):
    """Its own copies, built in the test: a fault that reads score.log stops the build, and so this test."""

    def test_score_log_is_neither_read_nor_copied(self):
        with tempfile.TemporaryDirectory() as temp:
            root = F.build(Path(temp) / "run")
            back = root / "back"
            copy = back / "synthia_coli"
            tables = checks.score_files(copy)
            tables["score/score.log"] = f"SYNTHETIC {F.genome_ids('synthia_coli')[0]} Column_3\n".encode()
            F.pack(copy / "score.tar", tables)
            reading = root / "reading-again"
            rep = checks.gather(back, root / "lists.tsv", root / "release", reading)
            self.assertEqual(rep.failed, 0, "\n".join(rep.lines))
            self.assertEqual(sorted(p.name for p in (reading / "score" / "synthia_coli").iterdir()),
                             ["checks.tsv", "datasets.tsv", "draws.tsv"])
            self.assertFalse([p for p in reading.rglob("*.log") if p.name != "plans.log"])


if __name__ == "__main__":
    unittest.main()
