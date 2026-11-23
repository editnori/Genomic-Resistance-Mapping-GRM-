"""SYNTHETIC tests of the review's fixes to scripts/c2-cap/run_checks.py, before any C2 score is read.

Each copy is tests/c2_cap/test_run_checks.py's made-up run, and where a test says so it is laid out as c2-drive.sh
brings a species back: score/ packed as score.tar beside the rest, not unpacked. Every id, figure and sha256 here is
made up; no figure, id or count comes from a run.
"""

import pathlib
import shutil
import tarfile
import unittest
from unittest import mock

from tests.c2_cap import test_run_checks as base

checks = base.checks


class Copies(unittest.TestCase):
    setUp = base.Run.setUp
    tearDown = base.Run.tearDown
    folder = base.Run.folder
    complete = base.Run.complete
    seal = base.Run.seal
    run_checks = base.Run.run_checks
    fails = base.Run.fails

    def as_driver(self, slugs=None):
        """score/ packed as score.tar, the folder removed: the copy as c2-drive.sh leaves it."""
        for slug in self.r.slugs if slugs is None else slugs:
            f = self.folder(slug)
            with tarfile.open(f / "score.tar", "w") as archive:
                archive.add(f / "score", arcname="score")
            shutil.rmtree(f / "score")

    def audited(self):
        """run_checks under an audit: (lines by check, ok, score.tar members extracted, whole reads of score.tar)."""
        members, whole = [], []
        extract, read_bytes = tarfile.TarFile.extractfile, pathlib.Path.read_bytes

        def watch_extract(archive, member):
            members.append(getattr(member, "name", member))
            return extract(archive, member)

        def watch_read(path):
            if path.name == "score.tar":
                whole.append(str(path))
            return read_bytes(path)

        with mock.patch.object(tarfile.TarFile, "extractfile", watch_extract), \
                mock.patch.object(pathlib.Path, "read_bytes", watch_read):
            got, ok = self.run_checks()
        return got, ok, members, whole


class ScoreTarInPlace(Copies):
    def test_copies_as_the_driver_leaves_them_pass_with_score_tar_read_in_place(self):
        for slug in self.r.slugs:
            (self.folder(slug) / "score" / "score.log").write_text(f"{self.ids[slug][0]} PRJNA000002\n",
                                                                     encoding="utf-8")
        self.as_driver()
        got, ok, members, whole = self.audited()
        self.assertTrue(ok, "\n".join(got.values()))
        self.assertEqual(whole, [])
        self.assertTrue(members)
        self.assertFalse(any(m.endswith("score.log") for m in members), members)
        self.assertIn("0 genome ids and 0 BioProject accessions", got["genome ids and BioProjects"])

    def test_a_table_in_score_tar_that_differs_from_sha256sums_fails(self):
        slug = self.r.slugs[1]
        p = self.folder(slug) / "score" / "datasets.tsv"
        p.write_text(p.read_text(encoding="utf-8") + "\n", encoding="utf-8")
        self.as_driver()
        self.fails("SHA256SUMS")

    def test_a_genome_id_in_a_table_of_score_tar_fails(self):
        slug = self.r.slugs[3]
        p = self.folder(slug) / "score" / "checks.tsv"
        p.write_text(p.read_text(encoding="utf-8") + f"{self.ids[slug][2]}\tx\t1\t\n", encoding="utf-8")
        self.seal(slug)
        self.as_driver()
        got = self.fails("genome ids and BioProjects")
        self.assertIn("1 genome ids and 0 BioProject accessions", got["genome ids and BioProjects"])


class ScoreLogNeverHashed(Copies):
    def listed(self, name):
        f = self.folder(self.r.slugs[0])
        with open(f / "SHA256SUMS", "a", encoding="utf-8") as out:
            out.write(f"{checks.sha256(f / 'score' / 'score.log')}  {name}\n")
        opened = []
        real = checks.sha256

        def watch(path):
            opened.append(str(path))
            return real(path)

        with mock.patch.object(checks, "sha256", watch):
            got, ok = self.run_checks()
        self.assertFalse(ok)
        self.assertTrue(got["SHA256SUMS"].startswith("FAIL "), got["SHA256SUMS"])
        self.assertFalse(any(p.endswith("score.log") for p in opened), opened)
        return got

    def test_score_log_named_in_sha256sums_is_counted_and_never_hashed(self):
        got = self.listed("score/score.log")
        self.assertIn("score/score.log listed 1 times", got["SHA256SUMS"])

    def test_score_log_named_another_way_in_sha256sums_is_counted_and_never_hashed(self):
        got = self.listed("score/./score.log")
        self.assertIn("score/score.log listed 1 times", got["SHA256SUMS"])


class Times(Copies):
    def test_a_failed_with_no_time_fails(self):
        slug = self.r.slugs[5]
        (self.folder(slug) / "FAILED").write_text("", encoding="utf-8")
        got = self.fails("DONE")
        self.assertIn(f"{slug} FAILED with no time in date -u +%FT%TZ form", got["DONE"])

    def test_a_failed_time_in_another_form_fails(self):
        slug = self.r.slugs[5]
        (self.folder(slug) / "FAILED").write_text(f"FAILED {slug} at fit 2026-11-23 13:00:00\n", encoding="utf-8")
        self.fails("DONE")

    def test_a_done_time_in_another_form_fails_beside_a_later_failed(self):
        slug = self.r.slugs[5]
        self.seal(slug, time="2026-11-23T9:00:00Z")
        (self.folder(slug) / "FAILED").write_text(f"FAILED {slug} at fit 2026-11-23T10:00:00Z\n", encoding="utf-8")
        got = self.fails("DONE")
        self.assertNotIn("set aside", got["DONE"])


class GenomeIds(Copies):
    def test_a_genome_id_joined_to_a_word_fails(self):
        for i, text in enumerate(("genome_{}", "id{}", "fit-{}:")):
            with self.subTest(text=text):
                slug = self.r.slugs[i]
                p = self.folder(slug) / "logs" / "run.log"
                p.write_text(text.format(self.ids[slug][1]) + "\n", encoding="utf-8")
                got = self.fails("genome ids and BioProjects")
                self.assertIn("1 genome ids", got["genome ids and BioProjects"])
                p.unlink()

    def test_an_id_inside_a_longer_number_is_no_id(self):
        slug = self.r.slugs[2]
        (self.folder(slug) / "logs" / "run.log").write_text(f"1{self.ids[slug][1]}1\n", encoding="utf-8")
        got, ok = self.run_checks()
        self.assertTrue(ok, got["genome ids and BioProjects"])


class ListsTxt(Copies):
    def test_copies_without_lists_txt_take_it_rebuilt_and_held_to_sha256sums(self):
        for slug in self.r.slugs:
            (self.folder(slug) / "lists.txt").unlink()
        got, ok = self.run_checks()
        self.assertTrue(ok, "\n".join(got.values()))
        self.assertIn("lists.txt not in 8 copies, rebuilt", got["SHA256SUMS"])
        self.assertTrue(got["expect.py"].startswith("pass "), got["expect.py"])

    def test_a_rebuilt_lists_txt_that_is_not_the_runs_fails(self):
        slug = self.r.slugs[2]
        p = self.folder(slug) / "lists.txt"
        p.write_text(p.read_text(encoding="utf-8").replace("twin\t3", "twin\t2"), encoding="utf-8")
        self.seal(slug)
        p.unlink()
        self.fails("SHA256SUMS")

    def test_a_rebuilt_lists_txt_with_no_line_in_sha256sums_fails(self):
        slug = self.r.slugs[4]
        f = self.folder(slug)
        lines = (f / "SHA256SUMS").read_text(encoding="utf-8").splitlines(keepends=True)
        (f / "SHA256SUMS").write_text("".join(l for l in lines if not l.rstrip().endswith("  lists.txt")),
                                      encoding="utf-8")
        (f / "lists.txt").unlink()
        self.fails("SHA256SUMS")


if __name__ == "__main__":
    unittest.main()
