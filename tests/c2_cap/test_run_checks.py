"""SYNTHETIC tests of the run's checks before any score is read, scripts/c2-cap/run_checks.py.

Each run copy here is tests/c2_cap/page_common.py's made-up reading, completed with the files run_species.sh
writes beside it (lists.txt, roots.json, the study fold files' table, DONE, SHA256SUMS) and a made-up release
whose genome ids are made up. No figure, id or count comes from a run.
"""

import json
import shutil
import tempfile
import unittest
from pathlib import Path

from tests.c2_cap import common, page_common as pc

checks = common.load("run_checks")
STUDY_COLUMNS = ["dataset", "arm", "added", "train_added", "unused_bioproject", "unused_unknown", "unused_duplicate",
                 "written"]


class Run(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.r = pc.Reading(self.tmp).write()
        self.release = self.tmp / "release"
        self.ids = {}
        for i, slug in enumerate(self.r.slugs):
            self.ids[slug] = [f"{9000 + i}.{n}" for n in range(1, 4)]
            p = self.release / "species" / f"{slug}.json"
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(json.dumps({"genomes": [{"id": g} for g in self.ids[slug]]}), encoding="utf-8")
            self.complete(slug)
        self.selection = {s: (5000, 2000, 2096) for s in self.r.slugs}

    def tearDown(self):
        self._tmp.cleanup()

    def folder(self, slug):
        return self.r.runs / slug

    def complete(self, slug, time="2026-11-23T14:00:00Z"):
        f = self.folder(slug)
        study = [d for d in self.r.datasets(slug) if d in self.r.study]
        (f / "lists.txt").write_text(f"species\t{slug}\ndatasets\t{','.join(self.r.datasets(slug))}\ndraws\t3\n"
                                     "twin\t3\narms\tcap-all\nsplits\tlineage,study\nscore\tyes\n", encoding="utf-8")
        (f / "results" / "roots.json").write_text(json.dumps({"roots": ["cap-all", "cap-2000"],
                                                              "fits": {"cap-all": "cap-all"}, "same_as": {}}),
                                                  encoding="utf-8")
        if study:
            pc.tsv(f / "results" / "study-files-cap-all.tsv",
                   [{"dataset": d, "arm": "cap-all", "added": 400, "train_added": 300, "unused_bioproject": 50,
                     "unused_unknown": 40, "unused_duplicate": 10, "written": 1} for d in study], STUDY_COLUMNS)
        (f / "score" / "score.log").write_text("score.py's output, never read\n", encoding="utf-8")
        self.seal(slug, time)

    def seal(self, slug, time="2026-11-23T14:00:00Z"):
        f = self.folder(slug)
        (f / "SHA256SUMS").write_text(pc.sums(f), encoding="utf-8")
        b = checks.binary_of(f)
        (f / "DONE").write_text(f"DONE {slug} {time} binary {b['grm']} code {b['code']} fits 4 SHA256SUMS "
                                f"{checks.sha256(f / 'SHA256SUMS')}\n", encoding="utf-8")

    def run_checks(self):
        lines, ok = checks.run(self.r.runs, self.release, pc.read(self.r.lists), self.selection)
        return {line.split(" ", 1)[1].split(":", 1)[0]: line for line in lines}, ok

    def fails(self, name):
        got, ok = self.run_checks()
        self.assertFalse(ok)
        self.assertTrue(got[name].startswith("FAIL "), got[name])
        return got

    def test_a_clean_run_passes_every_check(self):
        got, ok = self.run_checks()
        self.assertTrue(ok, "\n".join(got.values()))
        self.assertEqual(list(got), ["species", "SHA256SUMS", "DONE", "one binary", "expect.py", "cohorts",
                                     "genome ids and BioProjects", "refused lineage draws"])
        self.assertIn("8 copies for the 8 species", got["species"])
        self.assertIn("score/score.log listed 0 times", got["SHA256SUMS"])
        self.assertIn(f"{self.r.refused}: 3", got["refused lineage draws"])
        self.assertIn("1 datasets with the lineage split unavailable, 3 refused default slots",
                      got["refused lineage draws"])

    def test_a_missing_species_copy_fails(self):
        shutil.rmtree(self.folder(self.r.slugs[4]))
        got = self.fails("species")
        self.assertIn(self.r.slugs[4], got["species"])

    def test_a_changed_file_fails_its_sha256sums(self):
        p = self.folder(self.r.slugs[2]) / "results" / "added-labelled.tsv"
        p.write_text(p.read_text() + "\n", encoding="utf-8")
        self.fails("SHA256SUMS")

    def test_score_log_named_in_sha256sums_fails(self):
        f = self.folder(self.r.slugs[0])
        log = f / "score" / "score.log"
        with open(f / "SHA256SUMS", "a", encoding="utf-8") as out:
            out.write(f"{checks.sha256(log)}  score/score.log\n")
        got = self.fails("SHA256SUMS")
        self.assertIn("score/score.log listed 1 times", got["SHA256SUMS"])

    def test_a_done_with_another_code_sha256_fails(self):
        f = self.folder(self.r.slugs[6])
        done = (f / "DONE").read_text().replace(pc.CODE, "f" * 64)
        (f / "DONE").write_text(done, encoding="utf-8")
        self.fails("DONE")

    def test_a_failed_beside_a_newer_done_is_named_and_set_aside(self):
        f = self.folder(self.r.slugs[5])
        (f / "FAILED").write_text(f"FAILED {self.r.slugs[5]} at fit 2026-11-23T13:00:00Z\n", encoding="utf-8")
        got, ok = self.run_checks()
        self.assertTrue(ok, got["DONE"])
        self.assertIn(f"set aside: {self.r.slugs[5]} FAILED at 2026-11-23T13:00:00Z, DONE at 2026-11-23T14:00:00Z",
                      got["DONE"])

    def test_a_failed_newer_than_its_done_fails(self):
        f = self.folder(self.r.slugs[5])
        (f / "FAILED").write_text(f"FAILED {self.r.slugs[5]} at fit 2026-11-23T15:00:00Z\n", encoding="utf-8")
        got = self.fails("DONE")
        self.assertIn(f"{self.r.slugs[5]} FAILED at 2026-11-23T15:00:00Z with no newer DONE", got["DONE"])

    def test_species_on_two_binaries_fail_and_each_binary_is_named(self):
        slug = self.r.slugs[7]
        self.r.binary[slug] = (pc.COMMIT, pc.GRM, "e" * 64)
        self.r.write_species(slug)
        self.seal(slug)
        got = self.fails("one binary")
        self.assertIn(f"code {'e' * 64}: 1 species ({slug})", got["one binary"])
        self.assertIn(f"code {pc.CODE}: 7 species", got["one binary"])

    def test_an_expect_tsv_that_is_not_expect_pys_fails(self):
        slug = self.r.slugs[1]
        self.r.edit(self.folder(slug) / "expect.tsv",
                    lambda rows: [dict(r, role="draw") if r["role"] == "none-added" else r for r in rows])
        self.seal(slug)
        self.fails("expect.py")

    def test_a_cohort_past_4096_fails(self):
        slug = self.r.slugs[3]
        p = self.folder(slug) / "results" / "cohorts.json"
        doc = json.loads(p.read_text())
        doc["arms"]["cap-all"].update(genomes=4097, lines=4097, added=2097)
        p.write_text(json.dumps(doc), encoding="utf-8")
        self.seal(slug)
        got = self.fails("cohorts")
        self.assertIn(f"{slug}: cap-all holds 4097 genomes, past 4,096", got["cohorts"])

    def test_a_cohort_short_of_the_selection_rule_fails_and_fetch_failures_are_named(self):
        slug = self.r.slugs[0]
        self.selection[slug] = (5000, 2000, 2100)
        p = self.folder(slug) / "results" / "fetch.json"
        p.write_text(json.dumps({"failed_added": 4, "left_wider_cohorts": {"cap-4000": 0, "cap-all": 4}}),
                     encoding="utf-8")
        self.seal(slug)
        got = self.fails("cohorts")
        self.assertIn(f"{slug}: cap-all holds 4096 genomes and adds 2096, not 4100 and 2100", got["cohorts"])
        self.assertIn(f"{slug}: 4 added genomes failed the fetch", got["cohorts"])

    def test_a_genome_id_in_a_published_file_fails(self):
        slug = self.r.slugs[2]
        (self.folder(slug) / "logs" / "run.log").write_text(f"fitted on {self.ids[slug][1]}\n", encoding="utf-8")
        got = self.fails("genome ids and BioProjects")
        self.assertIn("1 genome ids and 0 BioProject accessions", got["genome ids and BioProjects"])

    def test_a_bioproject_accession_in_a_published_file_fails(self):
        slug = self.r.slugs[4]
        (self.folder(slug) / "logs" / "run.log").write_text("held out PRJNA000001\n", encoding="utf-8")
        got = self.fails("genome ids and BioProjects")
        self.assertIn("0 genome ids and 1 BioProject accessions", got["genome ids and BioProjects"])

    def test_score_log_is_not_scanned_or_read(self):
        slug = self.r.slugs[2]
        (self.folder(slug) / "score" / "score.log").write_text(f"{self.ids[slug][0]} PRJEB000001\n", encoding="utf-8")
        got, ok = self.run_checks()
        self.assertTrue(ok, got["genome ids and BioProjects"])

    def test_an_id_list_given_is_scanned_too(self):
        slug = self.r.slugs[3]
        (self.folder(slug) / "logs" / "run.log").write_text("added 8888.12\n", encoding="utf-8")
        extra = self.tmp / "cap-all.txt"
        extra.write_text("8888.12\n", encoding="utf-8")
        lines, ok = checks.run(self.r.runs, self.release, pc.read(self.r.lists), self.selection, lists=[extra])
        self.assertFalse(ok)
        self.assertTrue(any(line.startswith("FAIL genome ids and BioProjects: 1 genome ids") for line in lines))


if __name__ == "__main__":
    unittest.main()
