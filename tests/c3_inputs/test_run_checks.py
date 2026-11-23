"""scripts/c3-inputs/run_checks.py: the checks of the C3 run's copies, on SYNTHETIC copies.

Every copy, result file, genome id, binary and score here is made up in the
test, from tests/c3_inputs/reading.py's three made-up species. Each copy is
laid out as a species' copy comes back from its run: results/, check/,
logs/, expect.tsv, lists.txt, binary.txt, SHA256SUMS, DONE, and its score/
packed in score.tar.
"""

import builtins
import hashlib
import io
import json
import shutil
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests.c3_inputs import common
from tests.c3_inputs import reading as R

checks = common.load("run_checks")
expect = common.load("expect")
ARMS = "program,k31-program,k31-counts,k31-aa,unitigs"
IDS = {"synthia_coli": ["90001.11", "90001.12"], "fakeria_maxima": ["90002.21"], "mockella_minor": ["90003.31"]}


def sums(copy):
    names = sorted(p.relative_to(copy).as_posix() for top in ("results", "check", "score") if (copy / top).is_dir()
                   for p in (copy / top).rglob("*") if p.is_file() and p.suffix in (".json", ".tsv", ".txt"))
    names += ["binary.txt", "expect.tsv", "lists.txt", "logs/timing.tsv"]
    return "".join(f"{hashlib.sha256((copy / n).read_bytes()).hexdigest()}  {n}\n" for n in sorted(names))


def done(copy, slug, code=R.CODE):
    s = hashlib.sha256((copy / "SHA256SUMS").read_bytes()).hexdigest()
    n = sum(1 for _ in (copy / "results" / "picks").rglob("*.pick.json"))
    (copy / "DONE").write_text(f"DONE {slug} 2026-11-23T10:00:00Z binary {R.GRM} code {code} results {n} SHA256SUMS {s}\n")


def pack(copy):
    with tarfile.open(copy / "score.tar", "w") as tar:
        tar.add(copy / "score", arcname="score")
    shutil.rmtree(copy / "score")


def seal(copy, slug):
    """SHA256SUMS and DONE of a copy with score/ unpacked, then score/ packed into score.tar."""
    (copy / "SHA256SUMS").write_text(sums(copy))
    done(copy, slug)
    pack(copy)


def build(top):
    """SYNTHETIC copies of the three species under top/back, the lists at top/lists.tsv, the release at top/release."""
    top = Path(top)
    reading, lists = R.build(top / "r")
    back = top / "back"
    for slug in R.SPECIES:
        copy = back / slug
        ds = [d for d in R.FIRST if d.startswith(slug + "__")]
        a = [d for d in ds if d in R.IN_A]
        st = [d for d in ds if d in R.IN_S]
        rows = expect.rows(slug, ds, a, st, [42, 43, 44], [45, 46, 47], ARMS.split(","))
        picks = copy / "results" / "picks"
        for r in rows:
            dataset, arm, split, name, inputs = r[0], r[2], r[3], r[7], r[8]
            path = picks / name
            if split == "study" and (dataset, arm) in R.REFUSED:
                path = picks / (name[: -len(".pick.json")] + ".refused.json")
                rec = {"format": "grm-c3-refused-1", "dataset": dataset, "arm": arm, "scheme": "study-whole", "seed": 42,
                       "reason": "no candidate", "fold": 2, "training_genomes": 30, "support_floor": 3}
            else:
                rec = {"format": "grm-c9-pick-1", "dataset": dataset}
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(rec))
            if inputs != "-":
                (picks / inputs).write_text(json.dumps({"format": "grm-c3-inputs-1", "dataset": dataset}))
        (copy / "results" / "builds").mkdir(parents=True)
        for kind in ("counts", "aa", "unitigs"):
            (copy / "results" / "builds" / f"{kind}.build.json").write_text(json.dumps(R.build_record(slug, kind)))
        (copy / "check").mkdir()
        (copy / "check" / "refit.tsv").write_text("dataset\tpassed\nmade-up\t1\n")
        (copy / "logs").mkdir()
        (copy / "logs" / "timing.tsv").write_text("step\twall_s\tpeak_kb\nbuild\t1.0\t10\n")
        (copy / "logs" / "run.log").write_text("made up: [genome id] scrubbed\n")
        (copy / "expect.tsv").write_text(checks.expect_text(rows))
        (copy / "lists.txt").write_text(f"species\t{slug}\ndatasets\t{','.join(ds)}\na\t{','.join(a)}\nseeds\t42,43,44\n"
                                        f"twin\t45,46,47\narms\t{ARMS}\ngraph\tgrm\n")
        shutil.copy(reading / "run" / slug / "binary.txt", copy / "binary.txt")
        shutil.copytree(reading / "score" / slug, copy / "score")
        (copy / "score" / "score.log").write_text("checks: made up\n")
        seal(copy, slug)
    for slug, ids in IDS.items():
        path = top / "release" / "species" / f"{slug}.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"genomes": [{"id": i} for i in ids]}))
    return back, lists, top / "release"


def lines(rep):
    return {line.split(":")[0]: line for line in rep.lines if line.startswith(("PASS", "FAIL"))}


class RunChecks(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.back, self.lists, self.release = build(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def run_checks(self, **kw):
        return checks.run(self.back, self.lists, self.release, **kw)

    def unpack(self, slug):
        with tarfile.open(self.back / slug / "score.tar") as tar:
            tar.extractall(self.back / slug, filter="data")

    def repack(self, slug):
        """Pack score/ again and seal the copy's SHA256SUMS and DONE."""
        copy = self.back / slug
        (copy / "score.tar").unlink()
        seal(copy, slug)

    def test_clean_copies_pass_every_check_with_score_in_its_tar(self):
        rep = self.run_checks()
        self.assertEqual(rep.failed, 0, "\n".join(rep.lines))
        self.assertEqual(rep.lines[-1], "CHECKS 6 of 6")
        got = lines(rep)
        self.assertIn("1 of 1 refusal records", got["PASS refusals"])
        self.assertIn("5 constant-model rows in draws.tsv", got["PASS refusals"])
        self.assertIn("0 genome ids in", got["PASS genome ids"])
        self.assertIn("against 4 release ids of 3 species lists", got["PASS genome ids"])

    def test_clean_copies_pass_with_score_unpacked(self):
        for slug in R.SPECIES:
            self.unpack(slug)
        rep = self.run_checks()
        self.assertEqual(rep.failed, 0, "\n".join(rep.lines))

    def test_the_printout_holds_counts_and_no_dataset(self):
        text = "\n".join(self.run_checks().lines)
        for d, *_ in R.L:
            self.assertNotIn(d, text)

    def test_a_missing_species_fails(self):
        shutil.rmtree(self.back / "mockella_minor")
        self.assertIn("FAIL species", lines(self.run_checks()))

    def test_a_file_changed_since_its_sum_fails(self):
        path = self.back / "fakeria_maxima" / "results" / "builds" / "counts.build.json"
        path.write_text(path.read_text().replace("5000", "5001"))
        got = lines(self.run_checks())
        self.assertIn("FAIL SHA256SUMS", got)
        self.assertIn("2 of 3 species", got["FAIL SHA256SUMS"])

    def test_a_score_file_changed_inside_score_tar_fails(self):
        self.unpack("synthia_coli")
        (self.back / "synthia_coli" / "score" / "checks.tsv").write_text("file\tcheck\tpassed\tdetail\n")
        copy = self.back / "synthia_coli"
        (copy / "score.tar").unlink()
        pack(copy)
        self.assertIn("FAIL SHA256SUMS", lines(self.run_checks()))

    def test_a_done_naming_another_code_fails(self):
        path = self.back / "synthia_coli" / "DONE"
        path.write_text(path.read_text().replace(R.CODE, "f" * 64))
        got = lines(self.run_checks())
        self.assertIn("FAIL DONE and binary.txt", got)

    def test_a_failed_file_not_set_aside_fails_and_one_set_aside_is_named(self):
        (self.back / "mockella_minor" / "FAILED").write_text("FAILED mockella_minor at build 2026-11-23T08:00:00Z\n")
        self.assertIn("FAIL DONE and binary.txt", lines(self.run_checks()))
        rep = self.run_checks(set_aside=("mockella_minor",))
        self.assertEqual(rep.failed, 0, "\n".join(rep.lines))
        self.assertIn("  mockella_minor: FAILED at build at 2026-11-23T08:00:00Z, before DONE at 2026-11-23T10:00:00Z: "
                      "named and set aside", rep.lines)

    def test_a_failed_file_after_its_done_is_not_set_aside(self):
        (self.back / "mockella_minor" / "FAILED").write_text("FAILED mockella_minor at build 2026-11-23T11:00:00Z\n")
        self.assertIn("FAIL DONE and binary.txt", lines(self.run_checks(set_aside=("mockella_minor",))))

    def test_an_expect_list_that_is_not_expect_pys_fails(self):
        path = self.back / "fakeria_maxima" / "expect.tsv"
        path.write_text("".join(path.read_text().splitlines(keepends=True)[:-1]))
        self.repack_after_unpack("fakeria_maxima")
        self.assertIn("FAIL expect.py", lines(self.run_checks()))

    def test_a_result_file_more_than_its_rows_fails(self):
        extra = self.back / "synthia_coli" / "results" / "picks" / "k31" / "synthia_coli__drug_a.lineage-s48.pick.json"
        extra.write_text("{}")
        self.repack_after_unpack("synthia_coli")
        got = lines(self.run_checks())
        self.assertIn("FAIL expect.py", got)
        self.assertIn(", 1 more", got["FAIL expect.py"])

    def test_a_refused_default_fails(self):
        picks = self.back / "synthia_coli" / "results" / "picks"
        pick = picks / "k31" / "synthia_coli__drug_a.study-whole-s42.pick.json"
        pick.rename(pick.with_name("synthia_coli__drug_a.study-whole-s42.refused.json"))
        self.repack_after_unpack("synthia_coli")
        self.assertIn("FAIL expect.py", lines(self.run_checks()))

    def test_a_refusal_for_another_reason_fails(self):
        path = self.back / "synthia_coli" / "results" / "picks" / "program" / "synthia_coli__drug_a.study-whole-s42.refused.json"
        rec = json.loads(path.read_text())
        rec["reason"] = "out of memory"
        path.write_text(json.dumps(rec))
        self.repack_after_unpack("synthia_coli")
        self.assertIn("FAIL refusals", lines(self.run_checks()))

    def test_constant_rows_with_no_refusal_record_fail(self):
        path = self.back / "synthia_coli" / "results" / "picks" / "program" / "synthia_coli__drug_a.study-whole-s42.refused.json"
        path.rename(path.with_name("synthia_coli__drug_a.study-whole-s42.pick.json"))
        self.repack_after_unpack("synthia_coli")
        self.assertIn("FAIL refusals", lines(self.run_checks()))

    def test_a_genome_id_in_score_tar_is_counted(self):
        self.unpack("fakeria_maxima")
        path = self.back / "fakeria_maxima" / "score" / "checks.tsv"
        path.write_text(path.read_text() + "made up: 90002.21\n")
        self.repack("fakeria_maxima")
        got = lines(self.run_checks())
        self.assertIn("FAIL genome ids", got)
        self.assertIn("1 genome ids in", got["FAIL genome ids"])

    def test_a_genome_id_in_logs_or_a_path_named_with_ids_also_is_counted(self):
        (self.back / "synthia_coli" / "logs" / "run.log").write_text("made up: 90001.12 and 90003.31\n")
        got = lines(self.run_checks())
        self.assertIn("2 genome ids in", got["FAIL genome ids"])
        extra = Path(self.tmp.name) / "page.md"
        extra.write_text("made up: 90001.11 but not 90001.111\n")
        (self.back / "synthia_coli" / "logs" / "run.log").write_text("made up: [genome id]\n")
        got = lines(self.run_checks(ids_also=(str(extra),)))
        self.assertIn("1 genome ids in", got["FAIL genome ids"])

    def test_score_log_is_never_opened_in_score_or_taken_out_of_score_tar(self):
        for packed in (True, False):
            self.tearDown()
            self.setUp()
            self.unpack("fakeria_maxima")
            (self.back / "fakeria_maxima" / "score" / "score.log").write_text("made up: 90002.21\n")
            self.repack("fakeria_maxima")
            if not packed:
                for slug in R.SPECIES:
                    self.unpack(slug)
            opened, taken = [], []
            real_open, real_io_open, real_extract = builtins.open, io.open, tarfile.TarFile.extractfile

            def guard(real):
                def wrapped(file, *args, **kwargs):
                    if str(file).endswith("score/score.log"):
                        opened.append(str(file))
                    return real(file, *args, **kwargs)
                return wrapped

            def extract(archive, member):
                if getattr(member, "name", member).endswith("score/score.log"):
                    taken.append(member)
                return real_extract(archive, member)

            with mock.patch("builtins.open", guard(real_open)), mock.patch("io.open", guard(real_io_open)), \
                    mock.patch.object(tarfile.TarFile, "extractfile", extract):
                rep = self.run_checks()
            self.assertEqual(rep.failed, 0, "\n".join(rep.lines))
            self.assertEqual((opened, taken), ([], []), packed)
        copy = self.back / "fakeria_maxima"
        (copy / "SHA256SUMS").write_text((copy / "SHA256SUMS").read_text() + f"{'0' * 64}  score/score.log\n")
        done(copy, "fakeria_maxima")
        self.assertIn("FAIL SHA256SUMS", lines(self.run_checks()))

    def test_a_genome_id_after_a_letter_or_underscore_is_counted(self):
        for text in ("made up: genome_90001.11\n", "made up: id90001.12\n"):
            (self.back / "synthia_coli" / "logs" / "run.log").write_text(text)
            got = lines(self.run_checks())
            self.assertIn("1 genome ids in", got.get("FAIL genome ids", ""), text)
        (self.back / "synthia_coli" / "logs" / "run.log").write_text("made up: 190001.11 and 0.9000111\n")
        self.assertIn("PASS genome ids", lines(self.run_checks()))

    def test_a_failed_with_no_time_of_run_species_form_is_not_set_aside(self):
        for text in ("FAILED mockella_minor at build 08:00\n", "FAILED mockella_minor at build 1\n"):
            (self.back / "mockella_minor" / "FAILED").write_text(text)
            rep = self.run_checks(set_aside=("mockella_minor",))
            self.assertIn("FAIL DONE and binary.txt", lines(rep), text)
            self.assertFalse([line for line in rep.lines if "named and set aside" in line], text)

    def test_a_failed_at_the_second_of_its_done_is_not_set_aside(self):
        (self.back / "mockella_minor" / "FAILED").write_text("FAILED mockella_minor at build 2026-11-23T10:00:00Z\n")
        rep = self.run_checks(set_aside=("mockella_minor",))
        self.assertIn("FAIL DONE and binary.txt", lines(rep))
        self.assertFalse([line for line in rep.lines if "named and set aside" in line])

    def test_a_sha256sums_listing_a_file_beyond_the_runs_rule_fails(self):
        slug = "mockella_minor"
        copy = self.back / slug
        self.unpack(slug)
        (copy / "score.tar").unlink()
        (copy / "results" / "extra.bin").write_bytes(b"SYNTHETIC")
        (copy / "SHA256SUMS").write_text(sums(copy) + f"{hashlib.sha256(b'SYNTHETIC').hexdigest()}  results/extra.bin\n")
        done(copy, slug)
        pack(copy)
        self.assertIn("FAIL SHA256SUMS", lines(self.run_checks()))

    def test_a_species_run_without_one_arm_fails(self):
        slug = "synthia_coli"
        copy = self.back / slug
        self.unpack(slug)
        ds = [d for d in R.FIRST if d.startswith(slug + "__")]
        a = [d for d in ds if d in R.IN_A]
        st = [d for d in ds if d in R.IN_S]
        arms = [x for x in ARMS.split(",") if x != "k31-aa"]
        rows = expect.rows(slug, ds, a, st, [42, 43, 44], [45, 46, 47], arms)
        kept = {name for r in rows for name in (r[7], r[8])}
        for r in expect.rows(slug, ds, a, st, [42, 43, 44], [45, 46, 47], ARMS.split(",")):
            for name in (r[7], r[8]):
                path = copy / "results" / "picks" / name
                if name != "-" and name not in kept and path.is_file():
                    path.unlink()
        (copy / "expect.tsv").write_text(checks.expect_text(rows))
        text = (copy / "lists.txt").read_text()
        (copy / "lists.txt").write_text(text.replace(f"arms\t{ARMS}\n", f"arms\t{','.join(arms)}\n"))
        self.repack(slug)
        self.assertIn("FAIL expect.py", lines(self.run_checks()))

    def test_a_refusal_whose_constant_rows_miss_a_learner_fails(self):
        slug = "synthia_coli"
        self.unpack(slug)
        path = self.back / slug / "score" / "draws.tsv"
        rows = R.page.load(path)
        row = next(r for r in rows if r["constant"] == "1")
        row["constant"], row["why"] = "0", ""
        R.score.write_tsv(path, rows, R.score.DRAW_COLUMNS)
        self.repack(slug)
        self.assertIn("FAIL refusals", lines(self.run_checks()))

    def test_a_refusal_whose_mark_moved_to_a_lineage_row_fails(self):
        # The record's study run keeps four marks and its lineage draw takes the fifth: the rows still number five
        # a record, so only the count per refused study run can refuse it.
        slug = "synthia_coli"
        self.unpack(slug)
        path = self.back / slug / "score" / "draws.tsv"
        rows = R.page.load(path)
        row = next(r for r in rows if r["constant"] == "1")
        row["constant"], row["why"] = "0", ""
        other = next(r for r in rows if (r["dataset"], r["arm"], r["split"]) == (row["dataset"], row["arm"], "lineage"))
        other["constant"] = "1"
        R.score.write_tsv(path, rows, R.score.DRAW_COLUMNS)
        self.repack(slug)
        got = lines(self.run_checks())
        self.assertIn("FAIL refusals", got)
        self.assertIn("5 constant-model rows in draws.tsv", got["FAIL refusals"])

    def repack_after_unpack(self, slug):
        self.unpack(slug)
        self.repack(slug)


if __name__ == "__main__":
    unittest.main()
