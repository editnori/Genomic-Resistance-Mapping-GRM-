"""scripts/c7-outside/run_checks.py on SYNTHETIC copies: each check fails on what it names, a first copy set aside
and a FAILED older than DONE are named and set aside by time, and gather writes the reading only when every check
passes. Every species, dataset, genome id and file here is made up (page_fixture.py)."""

import shutil
import tempfile
import unittest
from pathlib import Path

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

    def test_a_withdrawn_arm_needs_a_reason_and_a_utc_time(self):
        for argv in (["--withdrawn", "elasticnet", "", "2026-11-23T09:00Z"],
                     ["--withdrawn", "elasticnet", "cost", "23 Nov"], ["--withdrawn", "twin", "cost", "2026-11-23T09:00Z"],
                     ["--withdrawn", "elasticnet", "cost"]):
            with self.assertRaises(SystemExit):
                checks.parse(argv)
        self.assertEqual(checks.parse(["--withdrawn", "elasticnet", "cost", "2026-11-23T09:00Z"])[1],
                         {"elasticnet": ("cost", "2026-11-23T09:00Z")})


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
