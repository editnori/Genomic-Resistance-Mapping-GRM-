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
rescore = common.load("rescore")
STOPPED = "synthia_gamma"
STOPPED_AT = "2026-11-23T12:30:00Z"


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

    def test_a_genome_id_in_failed_is_refused(self):
        def after(copies):
            (copies / "synthia_beta" / "FAILED").write_text(
                f"FAILED synthia_beta 1 fit 2026-11-23T11:00:00Z {fx.genome('synthia_beta', 1)}\n", encoding="utf-8")
        line = self.line(self.run_checks(after=after), "no genome id")
        self.assertTrue(line.startswith("FAIL"), line)
        self.assertIn(": 1 genome ids", line)

    def test_a_default_of_four_outer_folds_owes_four_choices_files(self):
        def change(slug, out):
            if slug != "synthia_beta":
                return
            rows = fx.expect_rows(slug, ["synthia_beta__synthetic_drug_3"])
            picks = out / "results" / "picks"
            for e in rows:
                if (e["split"], int(e["seed"])) != ("lineage", 43):
                    continue
                if e["arm"] == "default":
                    (picks / e["file"]).write_text(json.dumps({"format": "grm-c9-pick-1",
                                                               "folds": [{"fold": k} for k in range(4)]}),
                                                   encoding="utf-8")
                elif e["arm"] != "twin" and (picks / e["file"].replace("<k>", "4")).exists():
                    (picks / e["file"].replace("<k>", "4")).unlink()
        report = self.run_checks(change)
        self.assertTrue(self.line(report, "each species holds every default").startswith("PASS"))
        self.assertEqual(report[-1], "CHECKS 11 of 11")

    def test_a_done_naming_another_choices_count_is_refused(self):
        def after(copies):
            p = copies / "synthia_gamma" / "DONE"
            text = p.read_text(encoding="utf-8")
            n = int(text.split(" choices ", 1)[1].split(" ", 1)[0])
            p.write_text(text.replace(f" choices {n} ", f" choices {n + 1} "), encoding="utf-8")
        line = self.line(self.run_checks(after=after), "each DONE names")
        self.assertTrue(line.startswith("FAIL"), line)
        self.assertIn("synthia_gamma: DONE does not agree", line)

    def test_a_lists_txt_naming_other_twin_seeds_is_refused(self):
        def change(slug, out):
            if slug == "synthia_alpha":
                p = out / "lists.txt"
                p.write_text(p.read_text(encoding="utf-8").replace("twin\t45,46,47", "twin\t45,46"), encoding="utf-8")
        line = self.line(self.run_checks(change), "each lists.txt names")
        self.assertTrue(line.startswith("FAIL"), line)
        self.assertTrue(line.endswith(": synthia_alpha"), line)

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


def stop_at_score(copies, slug=STOPPED, drives=20, ended=None, lead_extra=""):
    """Make a species' copy one whose run stopped at its score step, as Klebsiella pneumoniae's did: FAILED at score,
    no DONE and no SHA256SUMS, score/ holding checks.tsv and score.log, logs/run.log with its drives, and
    SHA256SUMS.lead over every file but score/score.log. Returns the lead's sha256 and the folder of the rescored
    tables (each species' draws.tsv as the copies hold it, the stopped one's made again)."""
    out = copies / slug
    draws = (out / "score" / "draws.tsv").read_bytes()
    for name in ("DONE", "SHA256SUMS", "score/draws.tsv", "score/datasets.tsv"):
        (out / name).unlink()
    (out / "FAILED").write_text(f"FAILED {slug} at score {STOPPED_AT}\n", encoding="utf-8")
    lines = [f"2026-11-23T10:00:00Z start {slug}, 4 jobs, 2 drives at once",
             f"2026-11-23T10:01:00Z drives: {drives} outer folds, {drives} to run, 2 at once, 2 threads each"]
    lines += ended if ended is not None else [f"2026-11-23T11:00:00Z drive SYNTHETIC.lineage-s42.f{k} ok"
                                              for k in range(drives)]
    lines += ["2026-11-23T12:29:00Z score: start", "2026-11-23T12:30:00Z score: stopped; its output is in score/score.log",
              "2026-11-23T12:30:00Z FAILED at score"]
    (out / "logs" / "run.log").write_text("\n".join(lines) + "\n", encoding="utf-8")
    names = sorted(q.relative_to(out).as_posix() for q in out.rglob("*") if q.is_file())
    (out / "SHA256SUMS.lead").write_text("".join(f"{fx.sha256(out / n)}  {n}\n" for n in names
                                                 if n != "score/score.log") + lead_extra, encoding="utf-8")
    tables = copies.parent / "rescored"
    for q in sorted(copies.iterdir()):
        (tables / q.name).mkdir(parents=True)
        (tables / q.name / "draws.tsv").write_bytes(draws if q.name == slug else (q / "score" / "draws.tsv").read_bytes())
    return fx.sha256(out / "SHA256SUMS.lead"), tables


class Allowance(Case):
    """--failed-at-score: one species' run stopped at its score step (SYNTHETIC copies)."""

    def allowed(self, after=None, flag=True, **kw):
        copies, lists, release = fx.copies(self.top)
        lead, tables = stop_at_score(copies, **kw)
        if after:
            after(copies)
        lines, aside = checks.check(copies, lists, release, STOPPED if flag else None, lead if flag else None,
                                    tables)
        return checks.report(lines, aside)

    def test_a_species_stopped_at_score_passes_under_the_named_allowance(self):
        report = self.allowed()
        self.assertEqual(report[-1], "CHECKS 15 of 15", "\n".join(report))
        self.assertIn(f"FAILED at score {STOPPED_AT}; DONE 0, SHA256SUMS 0", self.line(report, f"{STOPPED}: FAILED"))
        self.assertIn("20 to run over 1 drives lines, 20 ended ok (20 named once), 0 ended otherwise",
                      self.line(report, f"{STOPPED}: logs/run.log"))
        self.assertIn("score/score.log named 0; its sha256 is the one given", self.line(report, f"{STOPPED}: every line"))
        self.assertIn(f"commit {fx.COMMIT}, code {fx.CODE}", self.line(report, f"{STOPPED}: binary.txt"))
        self.assertIn(f"2 of 3; {STOPPED} stopped at its score step", self.line(report, "each species ended DONE"))

    def test_a_species_stopped_at_score_fails_without_the_allowance(self):
        report = self.allowed(flag=False)
        self.assertTrue(self.line(report, "each species ended DONE").startswith("FAIL"))
        self.assertFalse([x for x in report if x.split(" ", 1)[1].startswith(f"{STOPPED}: ")])

    def test_the_allowance_refuses_a_lead_sums_of_another_sha256(self):
        copies, lists, release = fx.copies(self.top)
        _, tables = stop_at_score(copies)
        lines, aside = checks.check(copies, lists, release, STOPPED, "0" * 64, tables)
        line = self.line(checks.report(lines, aside), f"{STOPPED}: every line")
        self.assertTrue(line.startswith("FAIL") and line.endswith("is not the one given"), line)

    def test_the_allowance_refuses_a_drive_ended_otherwise(self):
        ended = [f"2026-11-23T11:00:00Z drive SYNTHETIC.lineage-s42.f{k} ok" for k in range(20)]
        ended.append("2026-11-23T11:00:00Z drive SYNTHETIC.lineage-s43.f0 FAILED, rc 1")
        line = self.line(self.allowed(ended=ended), f"{STOPPED}: logs/run.log")
        self.assertTrue(line.startswith("FAIL"), line)
        self.assertIn("20 ended ok (20 named once), 1 ended otherwise", line)

    def test_the_allowance_refuses_fewer_drives_ended_ok_than_to_run(self):
        ended = [f"2026-11-23T11:00:00Z drive SYNTHETIC.lineage-s42.f{k} ok" for k in range(19)]
        line = self.line(self.allowed(ended=ended), f"{STOPPED}: logs/run.log")
        self.assertTrue(line.startswith("FAIL"), line)
        self.assertIn("20 to run over 1 drives lines, 19 ended ok", line)

    def test_the_allowance_refuses_score_log_in_the_lead_sums_and_never_opens_it(self):
        opened = []
        real_open = builtins.open

        def wrapped(file, *args, **kwargs):
            if str(file).endswith("score/score.log"):
                opened.append(str(file))
            return real_open(file, *args, **kwargs)
        copies, lists, release = fx.copies(self.top)
        lead, tables = stop_at_score(copies, lead_extra=f"{'0' * 64}  score/score.log\n")
        with mock.patch("builtins.open", wrapped):
            lines, aside = checks.check(copies, lists, release, STOPPED, lead, tables)
        line = self.line(checks.report(lines, aside), f"{STOPPED}: every line")
        self.assertTrue(line.startswith("FAIL"), line)
        self.assertIn("score/score.log named 1", line)
        self.assertEqual(opened, [])

    def test_the_allowance_refuses_a_failed_line_at_another_step(self):
        def after(copies):
            (copies / STOPPED / "FAILED").write_text(f"FAILED {STOPPED} at drives {STOPPED_AT}\n", encoding="utf-8")
        line = self.line(self.allowed(after=after), f"{STOPPED}: FAILED")
        self.assertTrue(line.startswith("FAIL"), line)

    def test_the_allowance_refuses_a_done_beside_its_failed(self):
        def after(copies):
            (copies / STOPPED / "DONE").write_text("DONE SYNTHETIC\n", encoding="utf-8")
        line = self.line(self.allowed(after=after), f"{STOPPED}: FAILED")
        self.assertTrue(line.startswith("FAIL") and line.endswith("DONE 1, SHA256SUMS 0"), line)

    def test_the_allowance_refuses_another_commit(self):
        copies, lists, release = fx.copies(self.top)
        path = copies / STOPPED / "binary.txt"
        path.write_text(f"commit {'d' * 40}\ngrm {fx.GRM}\ncode {fx.CODE}\n", encoding="utf-8")
        lead, tables = stop_at_score(copies)
        report = checks.report(*checks.check(copies, lists, release, STOPPED, lead, tables))
        self.assertTrue(self.line(report, f"{STOPPED}: binary.txt").startswith("FAIL"))
        self.assertTrue(self.line(report, "every species ran one commit").startswith("FAIL"))

    def test_the_allowance_refuses_a_file_changed_under_the_lead_sums(self):
        def after(copies):
            with open(copies / STOPPED / "check" / "refit.tsv", "a", encoding="utf-8") as f:
                f.write("SYNTHETIC\t0\n")
        line = self.line(self.allowed(after=after), f"{STOPPED}: every line")
        self.assertTrue(line.startswith("FAIL"), line)
        self.assertIn("1 differ or absent", line)


class Times(Case):
    """A time that equals a genome id as text is counted apart (SYNTHETIC copies)."""

    def test_a_seconds_inner_value_equal_to_an_id_is_counted_apart(self):
        gid = fx.genome("synthia_alpha", 2)

        def change(slug, out):
            if slug == "synthia_alpha":
                p = out / "results" / "picks" / "choices" / "synthia_alpha__synthetic_drug_1.lineage-s42.f0.choices.json"
                text = p.read_text(encoding="utf-8")
                self.assertIn('"seconds_inner": 1.0', text)
                p.write_text(text.replace('"seconds_inner": 1.0', f'"seconds_inner": {gid}', 1), encoding="utf-8")
        line = self.line(self.run_checks(change), "no genome id")
        self.assertTrue(line.startswith("PASS"), line)
        self.assertIn(": 0 genome ids", line)
        self.assertIn("counted apart: choices seconds_inner 1", line)

    def test_a_user_time_line_equal_to_an_id_is_counted_apart(self):
        gid = fx.genome("synthia_beta", 1)

        def change(slug, out):
            if slug == "synthia_beta":
                (out / "logs" / "drive").mkdir()
                (out / "logs" / "drive" / "SYNTHETIC.time").write_text(
                    f"\tCommand being timed: SYNTHETIC\n\tUser time (seconds): {gid}\n"
                    f"\tElapsed (wall clock) time (h:mm:ss or m:ss): 1:{gid}\n", encoding="utf-8")
        line = self.line(self.run_checks(change), "no genome id")
        self.assertTrue(line.startswith("PASS"), line)
        self.assertIn("time Elapsed (wall clock) time 1, time User time (seconds) 1", line)

    def test_an_id_on_another_line_of_a_time_file_still_fails(self):
        gid = fx.genome("synthia_beta", 1)

        def change(slug, out):
            if slug == "synthia_beta":
                (out / "logs" / "drive").mkdir()
                (out / "logs" / "drive" / "SYNTHETIC.time").write_text(
                    f"\tMaximum resident set size (kbytes): {gid}\n", encoding="utf-8")
        line = self.line(self.run_checks(change), "no genome id")
        self.assertTrue(line.startswith("FAIL"), line)
        self.assertIn(": 1 genome ids", line)

    def test_an_id_in_another_log_still_fails(self):
        gid = fx.genome("synthia_gamma", 3)

        def change(slug, out):
            if slug == "synthia_gamma":
                (out / "logs" / "run.log").write_text(f"2026-11-23T10:00:00Z SYNTHETIC {gid}\n", encoding="utf-8")
        line = self.line(self.run_checks(change), "no genome id")
        self.assertTrue(line.startswith("FAIL"), line)


class Tables(Case):
    """gather --tables and rescore.py, for a reading taken from tables scored again (SYNTHETIC copies)."""

    def test_gather_with_tables_takes_the_rescored_tables(self):
        copies, lists, release = fx.copies(self.top)
        lines, aside = checks.check(copies, lists, release)
        path = self.top / "checks.txt"
        path.write_text("\n".join(checks.report(lines, aside)) + "\n", encoding="utf-8")
        tables = self.top / "rescored"
        for slug, _, _ in fx.SPECIES:
            (tables / slug).mkdir(parents=True)
            for t in checks.TABLES:
                (tables / slug / t).write_text(f"SYNTHETIC {slug} {t}\n", encoding="utf-8")
        (tables / "rescore.txt").write_text("PASS SYNTHETIC\nSAME 6 of 6\n", encoding="utf-8")
        checks.gather(copies, lists, path, self.top / "reading", tables)
        for slug, _, _ in fx.SPECIES:
            for t in checks.TABLES:
                self.assertEqual((self.top / "reading" / "score" / slug / t).read_text(encoding="utf-8"),
                                 f"SYNTHETIC {slug} {t}\n")
        (tables / "rescore.txt").write_text("PASS SYNTHETIC\nSAME 5 of 6\n", encoding="utf-8")
        with self.assertRaises(SystemExit):
            checks.gather(copies, lists, path, self.top / "again", tables)

    def rescore_copies(self, change_sums=None):
        """Two SYNTHETIC copies of a made-up species from test_score.write_run: one finished, its three tables as
        score.py writes them with their SHA256SUMS lines, and one stopped at its score step."""
        from tests.c5_settings import test_score as ts
        lists = self.top / "lists.tsv"
        lists.write_text("dataset\tspecies_slug\tfirst_run\tin_M\nx\tsynthia_done\t1\t0\ny\tsynthia_stopped\t1\t0\n",
                         encoding="utf-8")
        for slug in ("synthia_done", "synthia_stopped"):
            made = self.top / "made" / slug
            made.mkdir(parents=True)
            args = ts.write_run(made)
            out = self.top / "copies" / slug
            (out / "results").mkdir(parents=True)
            (out / "expect.tsv").write_bytes((made / "expect.tsv").read_bytes())
            (made / "picks").rename(out / "results" / "picks")
            if slug == "synthia_done":
                rc, _ = ts.score([str(out / "expect.tsv"), str(out / "results" / "picks"), args[2], str(out / "score")])
                self.assertEqual(rc, 0)
                lines = [f"{fx.sha256(out / 'score' / t)}  score/{t}" for t in checks.TABLES]
                if change_sums:
                    lines = change_sums(lines)
                (out / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="utf-8")
        return self.top / "copies", lists, self.top / "made" / "synthia_done" / "release"

    def test_rescore_holds_each_finished_species_tables_to_its_sha256sums(self):
        copies, lists, release = self.rescore_copies()
        lines = rescore.rescore(copies, lists, release, self.top / "rescored")
        self.assertEqual(lines[-1], "SAME 3 of 3", "\n".join(lines))
        self.assertIn("PASS synthia_stopped: score.py ended 0; 3 tables, not held: the copy has no SHA256SUMS", lines)
        self.assertTrue(checks.rescored(self.top / "rescored"))

    def test_rescore_refuses_a_table_unlike_its_sha256sums_line(self):
        copies, lists, release = self.rescore_copies(lambda ls: [("0" * 64) + ls[1][64:] if i == 1 else x
                                                                 for i, x in enumerate(ls)])
        lines = rescore.rescore(copies, lists, release, self.top / "rescored")
        self.assertEqual(lines[-1], "SAME 2 of 3", "\n".join(lines))
        self.assertIn("FAIL synthia_done: score.py ended 0; 2 of 3 tables the same as its SHA256SUMS lines; differs: "
                      "draws.tsv", lines)
        self.assertFalse(checks.rescored(self.top / "rescored"))


if __name__ == "__main__":
    unittest.main()
