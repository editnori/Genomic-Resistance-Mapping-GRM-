"""The run script: its receipts, its refusals and every stage, driven by a
stand-in for grm (fake_grm.py) inside a temporary committed kit."""

import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
import marshal
import os
import shutil
import stat
import struct
import sys
import tempfile
import unittest
from pathlib import Path

from tests.outside_study import fixtures as F

FAKE_SWITCHES = ("FAKE_MISSING", "FAKE_TRANSPORT", "FAKE_SILENT", "FAKE_VERDICTS", "FAKE_UNCOUNTED",
                 "FAKE_DISTANCES", "FAKE_UNJUDGED", "FAKE_REFERENCE", "FAKE_KMERS", "FAKE_GENES",
                 "FAKE_DIGESTS", "FAKE_CALLS_VERSION", "FAKE_PREDICT_WRONG", "FAKE_PREDICT_K")
# The world's run: one genome the database answers with no sequence, one
# failing the species check, one that does not count, one borderline.
WORLD_RUN = {"FAKE_MISSING": "1001.121", "FAKE_UNCOUNTED": "1001.123",
             "FAKE_VERDICTS": json.dumps({"1001.122": "fail", F.BORDERLINE: "borderline"})}


def tree_digest(root):
    h = hashlib.sha256()
    for p in sorted(Path(root).rglob("*")):
        h.update(str(p.relative_to(root)).encode())
        if p.is_file():
            h.update(p.read_bytes())
    return h.hexdigest()


def recorded_genome(g, verdict="pass"):
    """A genome of the master root's recorded species check, with the fields
    `grm species-check` writes and the stand-in's default values."""
    return {"genome_id": g, "species": F.NAME, "verdict": verdict, "gc": 0.33, "gc_z": 0.0, "distance": 0.01,
            "ani": 0.99, "nearest_species": None, "nearest_distance": None, "inside_species": None,
            "reason": "far from the reference" if verdict == "fail" else None}


def genome_contents(w):
    """What each genome holds, for the stand-in: k-mers, gene calls and
    k-mer set digests, the cohort's from the reports' published bits and
    calls, the new genomes' from the world."""
    kmers = {g: [F.KMER[f] for f in (1, 2, 3) if g in F.COHORT_F[f]] for g in F.COHORT}
    kmers.update({g: [F.KMER[f] for f in (1, 2) if found.get(f)] for g, found in w["found"].items()})
    genes = dict(F.COHORT_HITS)
    for g in w["genomes"]:
        genes[g["id"]] = [F.BLAZ] if g["known"] else []
    genes.update(w["known_gene"]["hits"]["new"])
    digests = {g: F.digest_of(g) for g in list(F.COHORT) + list(w["availability"])}
    return kmers, genes, digests


class RunCase(unittest.TestCase):
    """A committed kit holding the world's tables and pins, a master root
    holding the cohort, and the stand-in grm."""

    world_args = {}

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = self.root = Path(self.tmp.name)
        self.w = self.make_world()
        self.reports = root / "reports"
        self.kit = F.kit_repo(root / "kit", self.w, self.reports)
        self.R = self.kit.R
        m = self.master = root / "master"
        for sub in ("fasta", "catalog", "curated", "species-check/sketches", f"species/{F.SPECIES}"):
            (m / sub).mkdir(parents=True)
        for g in F.COHORT:
            (m / "fasta" / f"{g}.fna.gz").write_bytes(b"cohort " + g.encode())
            (m / "species-check" / "sketches" / f"{g}.json").write_text("{}")
        for name, text in (("AMR_CDS.fa", ">x\nACGT\n"), ("ReferenceGeneCatalog.txt", "x\n"),
                           ("version.txt", "2026-08-07.1\n")):
            (m / "catalog" / name).write_text(text)
        (m / "curated/species.json").write_text(json.dumps([{"slug": F.SPECIES, "species": F.NAME}]))
        (m / f"species/{F.SPECIES}/genomes.txt").write_text("\n".join(F.COHORT) + "\n")
        (m / "species-check/species_check.json").write_text(json.dumps({
            "species": [{"species": F.NAME, "slug": F.SPECIES, "reference": "999.1"}],
            "genomes": [recorded_genome(g) for g in F.COHORT]}))
        grm = self.grm = root / "grm"
        grm.write_text(f"#!{sys.executable} -I\n" + (F.ROOT / "tests/outside_study/fake_grm.py").read_text())
        grm.chmod(grm.stat().st_mode | stat.S_IXUSR)
        self.work = root / "work"
        self.env = {k: os.environ.get(k) for k in FAKE_SWITCHES}
        for k in FAKE_SWITCHES:
            os.environ.pop(k, None)
        kmers, genes, digests = genome_contents(self.w)
        for name, doc in (("FAKE_KMERS", kmers), ("FAKE_GENES", genes), ("FAKE_DIGESTS", digests)):
            p = root / f"{name.lower()}.json"
            p.write_text(json.dumps(doc))
            os.environ[name] = str(p)

    def make_world(self):
        return F.world(**self.world_args)

    def tearDown(self):
        for k, v in self.env.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        self.tmp.cleanup()

    def argv(self, stage, *extra, work=None, commit=True, reports=None):
        argv = [stage, "--grm", str(self.grm), "--master", str(self.master), "--work", str(work or self.work),
                "--species", F.SPECIES, "--population", "primary", "--reports", str(reports or self.reports),
                "--jobs", "2"]
        if commit:
            argv += ["--expect-commit", self.kit.commit]
        return argv + list(extra)

    def run_stage(self, stage, *extra, R=None, **kw):
        err = io.StringIO()
        with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
            code = (R or self.R).main(self.argv(stage, *extra, **kw))
        return code, err.getvalue()

    def run_ok(self, *stages, **kw):
        for stage in stages:
            code, err = self.run_stage(stage, **kw)
            self.assertEqual(code, 0, f"{stage}: {err}")

    def assertRefused(self, result, pattern):
        code, err = result
        self.assertEqual(code, 2, err)
        self.assertRegex(err, pattern)

    def receipt(self, stage):
        return json.loads((self.work / "receipts" / f"{stage}.json").read_text())

    def attempts(self):
        return [json.loads(x) for x in (self.work / "receipts/attempts.jsonl").read_text().splitlines()]

    def run_args(self):
        return argparse.Namespace(grm=str(self.grm), master=str(self.master), work=str(self.work),
                                  species=F.SPECIES, population="primary", reports=str(self.reports), jobs=2,
                                  rehearsal=False, expect_commit=self.kit.commit, tables=None, ids=None,
                                  datasets=None)


class Setup(RunCase):
    def test_setup_derives_the_lists_and_never_writes_the_master(self):
        before = tree_digest(self.master)
        self.run_ok("setup")
        self.assertEqual(tree_digest(self.master), before)
        r = self.receipt("setup")
        want = sorted([g["id"] for g in self.w["genomes"]] + list(F.EXTRA), key=self.R.S.genome_key)
        self.assertEqual(len(want), 126)
        self.assertEqual(r["frozen_ids"], want)
        self.assertEqual(r["datasets"], [F.A, F.B])
        self.assertEqual((self.work / "ids.txt").read_text().split(), want)
        self.assertEqual((self.work / "datasets.txt").read_text().split(), [F.A, F.B])
        self.assertEqual(r["counts"]["frozen ids requested"], 126)
        self.assertEqual(r["counts"]["cohort genomes for the replay checks"], 8)
        self.assertEqual(r["counts"]["catalog files copied"], 3)
        self.assertEqual(r["grm"]["sha256"], self.R.S.sha256_file(self.grm))
        self.assertEqual(r["kit"]["commit"], self.kit.commit)
        self.assertEqual(r["kit"]["expect_commit"], self.kit.commit)
        self.assertEqual(r["kit"]["frozen_sha256"], self.R.S.sha256_file(self.kit.root / self.R.S.FROZEN_REL))
        self.assertEqual(sorted(r["kit"]["files"]), sorted(self.R.S.FROZEN_FILES))
        self.assertIn("work/ids.txt", r["outputs"])
        self.assertTrue(all(len(v) == 64 for v in r["inputs"].values()))
        root = self.work / "root"
        self.assertTrue((root / "fasta" / "999.1.fna.gz").is_symlink())
        self.assertEqual((root / f"species/outside__{F.SPECIES}/genomes.txt").read_text().split(), want)
        self.assertEqual((root / f"species/replay__{F.SPECIES}/genomes.txt").read_text().split(), F.COHORT)
        self.assertEqual([(a["stage"], a["event"]) for a in self.attempts()], [("setup", "start"), ("setup", "done")])

    def test_a_real_run_needs_the_expected_commit(self):
        self.assertRefused(self.run_stage("setup", commit=False), "a real run needs --expect-commit")
        self.assertRefused(self.run_stage("setup", "--expect-commit", "abc", commit=False), "is not 40 lower-case hex")
        self.assertRefused(self.run_stage("setup", "--expect-commit", "0" * 40, commit=False),
                           "not the expected commit " + "0" * 40)
        self.assertFalse((self.work / "receipts/setup.json").exists())

    def test_rehearsal_options_are_refused_in_a_real_run(self):
        for opt in ("--tables", "--ids", "--datasets"):
            self.assertRefused(self.run_stage("setup", opt, str(self.root)), f"{opt} is for a rehearsal only")

    def test_ids_and_datasets_go_together(self):
        self.assertRefused(self.run_stage("setup", "--rehearsal", "--ids", str(self.root)),
                           "--ids and --datasets go together")

    def test_no_option_moves_the_pins_or_the_kit_root(self):
        for flag in ("--pins", "--kit-root", "--frozen"):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                self.R.main(self.argv("setup", flag, "x"))

    def test_a_rehearsal_takes_cohort_genomes_only(self):
        ids, ds = self.root / "ids.txt", self.root / "datasets.txt"
        ids.write_text("1001.1\n")
        ds.write_text(f"{F.A}\n")
        self.assertRefused(self.run_stage("setup", "--rehearsal", "--ids", str(ids), "--datasets", str(ds),
                                          commit=False), "1001.1 is not a cohort genome")

    def test_a_bad_ids_file_is_refused(self):
        ids, ds = self.root / "ids.txt", self.root / "datasets.txt"
        ds.write_text(f"{F.A}\n")
        for text, pattern in (("999.1\n999.1\n", "listed twice"), ("999.1\nGCF_000001\n", "not a BV-BRC genome id")):
            ids.write_text(text)
            self.assertRefused(self.run_stage("setup", "--rehearsal", "--ids", str(ids), "--datasets", str(ds),
                                              commit=False), pattern)

    def test_a_dataset_of_another_species_is_refused(self):
        ids, ds = self.root / "ids.txt", self.root / "datasets.txt"
        ids.write_text("999.1\n")
        ds.write_text("other_species__drugA\n")
        self.assertRefused(self.run_stage("setup", "--rehearsal", "--ids", str(ids), "--datasets", str(ds),
                                          commit=False), "is not a testus_bacterium dataset")

    def test_an_id_already_in_the_master_fasta_is_refused(self):
        (self.master / "fasta" / "1001.2.fna.gz").write_bytes(b"x")
        self.assertRefused(self.run_stage("setup"), "1001.2 already has a FASTA file in the master root")

    def test_the_work_directory_inside_the_master_is_refused(self):
        self.assertRefused(self.run_stage("setup", work=self.master / "w"), "inside the master root")

    def test_no_catalog_is_refused_outside_a_rehearsal(self):
        shutil.rmtree(self.master / "catalog")
        self.assertRefused(self.run_stage("setup"), "no catalog")

    def test_a_report_not_matching_its_pin_is_refused(self):
        p = self.reports / f"{F.A}.json"
        p.write_bytes(p.read_bytes() + b"\n")
        self.assertRefused(self.run_stage("setup"), "is not the pinned")

    def test_a_non_empty_scratch_root_is_refused(self):
        (self.work / "root").mkdir(parents=True)
        (self.work / "root" / "left-over").write_text("x")
        self.assertRefused(self.run_stage("setup"), "exists and is not empty; start a new work directory")


class CohortInTheTables(RunCase):
    def make_world(self):
        w = F.world()
        w["candidates"].append(F.candidate("999.3", F.PRIMARY, "P1", "Atlantis", "", True))
        return w

    def test_a_cohort_genome_in_the_frozen_tables_is_refused(self):
        self.assertRefused(self.run_stage("setup"), "999.3 is a cohort genome")


class Receipts(RunCase):
    def test_a_stage_without_the_earlier_receipt_is_refused(self):
        self.assertRefused(self.run_stage("fetch"), "the setup stage has no receipt")
        self.assertRefused(self.run_stage("count"), "the setup stage has no receipt")

    def test_a_stage_never_runs_twice_in_one_work_directory(self):
        self.run_ok("setup", "fetch")
        before = self.receipt("fetch")
        log = (self.work / "receipts/attempts.jsonl").read_bytes()
        self.assertRefused(self.run_stage("fetch"), "fetch: this stage already has a receipt .* "
                                                    "start a new work directory")
        self.assertEqual(self.receipt("fetch"), before)
        self.assertRefused(self.run_stage("setup"), "setup: this stage already has a receipt")
        self.assertRefused(self.run_stage("all"), "setup: this stage already has a receipt")
        # A try that is refused for its receipt writes no line to the log.
        self.assertEqual((self.work / "receipts/attempts.jsonl").read_bytes(), log)
        events = [(a["stage"], a["event"]) for a in self.attempts()]
        self.assertEqual(events, [("setup", "start"), ("setup", "done"), ("fetch", "start"), ("fetch", "done")])

    def test_a_changed_setup_output_is_refused(self):
        self.run_ok("setup")
        ids = self.work / "ids.txt"
        ids.write_text("\n".join(ids.read_text().split()[:3]) + "\n")
        self.assertRefused(self.run_stage("fetch"), "work/ids.txt changed since the setup stage")

    def test_a_changed_output_is_refused(self):
        self.run_ok("setup", "fetch")
        av = self.work / "fetch" / "availability.tsv"
        av.write_text(av.read_text().replace("1001.4\tyes", "1001.4\tno"))
        self.assertRefused(self.run_stage("species-check"), "availability.tsv changed since the fetch stage")

    def test_another_grm_binary_is_refused(self):
        self.run_ok("setup")
        with open(self.grm, "a") as f:
            f.write("\n# changed\n")
        self.assertRefused(self.run_stage("fetch"), "grm binary differs")

    def test_other_arguments_are_refused(self):
        self.run_ok("setup")
        other = self.root / "reports2"
        shutil.copytree(self.reports, other)
        self.assertRefused(self.run_stage("fetch", reports=other), "arguments differ")

    def test_every_stage_passes_the_gate(self):
        self.run_ok("setup")
        (self.kit.root / "docs/outside-study-test.md").write_text("changed plan\n")
        for stage in ("fetch", "known-gene", "predict"):
            self.assertRefused(self.run_stage(stage), "outside-study-test.md has sha256")

    def test_a_new_commit_stops_the_run(self):
        self.run_ok("setup")
        old = self.kit.commit
        (self.kit.root / "docs/outside-study-test.md").write_text("changed plan\n")
        F.repin(self.kit.root)
        F.commit_all(self.kit)
        self.assertRefused(self.run_stage("fetch", "--expect-commit", old, commit=False),
                           "not the expected commit")
        self.assertRefused(self.run_stage("fetch", "--expect-commit", self.kit.commit, commit=False),
                           "arguments differ|kit commit changed")

    def test_the_loaded_score_script_must_be_the_kits_own(self):
        other = F.kit_repo(self.root / "other", self.w, self.root / "other-reports")
        self.R.S = other.S
        self.assertRefused(self.run_stage("setup"), "the loaded score.py .* is not")

    def test_a_planted_bytecode_cache_is_never_run(self):
        # A __pycache__ file for score.py whose header carries the source's
        # time and size, so the import system takes it for the compiled
        # source. Its code is score.py with another class minimum.
        kit = self.kit.root / "scripts/outside-study"
        src = kit / "score.py"
        code = compile(src.read_text().replace("\nMIN_CLASS = 25\n", "\nMIN_CLASS = 61\n"), str(src), "exec")
        st = src.stat()
        pyc = Path(importlib.util.cache_from_source(str(src)))
        pyc.parent.mkdir(exist_ok=True)
        pyc.write_bytes(importlib.util.MAGIC_NUMBER + struct.pack("<III", 0, int(st.st_mtime) & 0xFFFFFFFF,
                                                                  st.st_size & 0xFFFFFFFF) + marshal.dumps(code))
        self.assertEqual(F.load("score", kit).MIN_CLASS, 61, "the import system does run the planted file")
        R = F.load("run", kit)
        self.assertEqual(R.S.MIN_CLASS, 25)
        self.assertEqual(R.S.LOADED_SHA256, hashlib.sha256(src.read_bytes()).hexdigest())
        self.run_ok("setup", R=R)
        self.assertEqual(self.receipt("setup")["kit"]["score.py"], R.S.LOADED_SHA256)

    def test_a_copy_of_the_run_script_elsewhere_is_refused(self):
        d = self.root / "elsewhere"
        d.mkdir()
        shutil.copyfile(self.kit.root / "scripts/outside-study/run.py", d / "run.py")
        shutil.copyfile(self.kit.root / "scripts/outside-study/score.py", d / "score.py")
        with self.assertRaises(FileNotFoundError):
            F.load("run", d)  # it looks for score.py in its own kit, two folders up
        kit = d / "a" / "scripts" / "outside-study"
        kit.mkdir(parents=True)
        for n in ("run.py", "score.py"):
            shutil.copyfile(d / n, kit / n)
        self.assertRefused(self.run_stage("setup", R=F.load("run", kit)), "not in a git work tree")


class Fetch(RunCase):
    def test_a_database_answer_of_no_sequence_is_the_availability_result(self):
        os.environ["FAKE_MISSING"] = "1001.2"
        self.run_ok("setup", "fetch")
        av = {r["genome_id"]: r for r in self.R.S.read_tsv(self.work / "fetch/availability.tsv")}
        self.assertEqual(av["1001.2"]["fasta"], "no")
        self.assertEqual(av["1001.2"]["detail"], "1001.2: giving up after 5 attempts: empty FASTA response")
        self.assertEqual(sum(r["fasta"] == "yes" for r in av.values()), 125)
        self.assertEqual(self.receipt("fetch")["counts"]["fetched"], 125)

    def assertStops(self, pattern):
        self.run_ok("setup")
        self.assertRefused(self.run_stage("fetch"), pattern)
        self.assertFalse((self.work / "receipts/fetch.json").exists())
        self.assertFalse((self.work / "fetch/availability.tsv").exists())
        for k in ("FAKE_TRANSPORT", "FAKE_SILENT"):
            os.environ.pop(k, None)
        self.run_ok("fetch")
        events = [(a["stage"], a["event"]) for a in self.attempts() if a["stage"] == "fetch"]
        self.assertEqual(events, [("fetch", "start"), ("fetch", "refused"), ("fetch", "start"), ("fetch", "done")])

    def test_a_connection_failure_stops_the_stage_without_a_receipt(self):
        os.environ["FAKE_TRANSPORT"] = "1001.3"
        self.assertStops("1 genomes have no FASTA and no answer of no sequence .*1001.3.*connection refused")

    def test_a_genome_neither_fetched_nor_recorded_stops_the_stage(self):
        os.environ["FAKE_SILENT"] = "1001.3"
        self.assertStops("1001.3: no FASTA file after fetch, and no failure recorded")

    def test_the_no_sequence_answer(self):
        self.assertTrue(self.R.no_sequence_answer("1001.2: giving up after 5 attempts: empty FASTA response"))
        self.assertTrue(self.R.no_sequence_answer("empty FASTA response"))
        for m in ("HTTP status 503", "error sending request: connection refused",
                  "giving up after 5 attempts: empty FASTA response from the proxy", ""):
            self.assertFalse(self.R.no_sequence_answer(m), m)


class SpeciesCheckAndCount(RunCase):
    def test_fetch_species_check_and_count(self):
        os.environ.update(WORLD_RUN)
        self.run_ok("setup", "fetch", "species-check", "count")
        av = {r["genome_id"]: r["fasta"] for r in self.R.S.read_tsv(self.work / "fetch/availability.tsv")}
        self.assertEqual(av["1001.121"], "no")
        v = {r["genome_id"]: r["verdict"] for r in self.R.S.read_tsv(self.work / "species-check/verdicts.tsv")}
        self.assertEqual((v["1001.122"], v[F.BORDERLINE], v["1001.2"]), ("fail", "borderline", "pass"))
        c = {r["genome_id"]: r for r in self.R.S.read_tsv(self.work / "count/counted.tsv")}
        self.assertEqual(c["1001.123"]["counted"], "no")
        self.assertEqual(c["1001.123"]["kset_digest"], "")
        self.assertEqual(c["1001.5"]["kset_digest"], F.digest_of("999.2"))
        cd = {r["genome_id"]: r["kset_digest"] for r in self.R.S.read_tsv(self.work / "count/cohort-digests.tsv")}
        self.assertEqual(cd, {g: F.digest_of(g) for g in F.COHORT})
        self.assertEqual(self.receipt("fetch")["counts"]["fetched"], 125)
        self.assertEqual(self.receipt("species-check")["counts"]["fail (excluded)"], 1)
        self.assertEqual(self.receipt("count")["counts"][f"replay__{F.SPECIES}: counted"], 8)
        survivors, excluded = self.R.Run(self.run_args()).survivors()
        self.assertEqual(survivors, self.w["survivors"])
        self.assertEqual({g: e[0] for g, e in excluded.items()},
                         {"1001.121": "no sequence", "1001.122": "species check failed", "1001.123": "not counted"})

    def test_a_missing_kset_digest_is_refused(self):
        p = Path(os.environ["FAKE_DIGESTS"])
        p.write_text(json.dumps(dict(json.loads(p.read_text()), **{"1001.9": ""})))
        self.run_ok("setup", "fetch", "species-check")
        self.assertRefused(self.run_stage("count"), "records no k-mer set digest for 1001.9")

    def species_check_refused(self, pattern, rows):
        self.run_ok("setup", "fetch")
        self.assertRefused(self.run_stage("species-check"), pattern)
        changes = self.R.S.read_tsv(self.work / "species-check/cohort-changes.tsv")
        self.assertEqual([(r["genome_id"], r["recorded"], r["now"]) for r in changes], rows)
        self.assertFalse((self.work / "receipts/species-check.json").exists())
        self.assertFalse((self.work / "species-check/cohort-gc-moves.tsv").exists())
        self.assertFalse((self.work / "species-check/verdicts.tsv").exists())

    def test_a_recorded_gc_flag_that_moves_is_listed_and_the_run_goes_on(self):
        os.environ["FAKE_VERDICTS"] = json.dumps({"999.6": "borderline"})
        self.run_ok("setup", "fetch", "species-check")
        moves = self.R.S.read_tsv(self.work / "species-check/cohort-gc-moves.tsv")
        self.assertEqual([(r["genome_id"], r["recorded"], r["now"]) for r in moves],
                         [("999.6", "pass", "borderline")])
        self.assertFalse((self.work / "species-check/cohort-changes.tsv").exists())
        counts = self.receipt("species-check")["counts"]
        self.assertEqual(counts["recorded genomes with unchanged verdicts"], len(F.COHORT) - 1)
        self.assertEqual(counts["recorded genomes whose GC flag alone moved (listed)"], 1)
        self.assertIn("work/species-check/cohort-gc-moves.tsv", self.receipt("species-check")["outputs"])

    def test_no_recorded_flag_moving_writes_an_empty_list(self):
        self.run_ok("setup", "fetch", "species-check")
        self.assertEqual(self.R.S.read_tsv(self.work / "species-check/cohort-gc-moves.tsv"), [])
        counts = self.receipt("species-check")["counts"]
        self.assertEqual(counts["recorded genomes with unchanged verdicts"], len(F.COHORT))
        self.assertEqual(counts["recorded genomes whose GC flag alone moved (listed)"], 0)

    def test_a_recorded_borderline_that_moves_back_to_pass_is_listed(self):
        check = self.master / "species-check/species_check.json"
        doc = json.loads(check.read_text())
        doc["genomes"] = [recorded_genome(g, "borderline" if g == "999.6" else "pass") for g in F.COHORT]
        check.write_text(json.dumps(doc))
        self.run_ok("setup", "fetch", "species-check")
        moves = self.R.S.read_tsv(self.work / "species-check/cohort-gc-moves.tsv")
        self.assertEqual([(r["genome_id"], r["recorded"], r["now"]) for r in moves],
                         [("999.6", "borderline", "pass")])

    def test_a_recorded_verdict_that_moves_to_fail_stops_the_run(self):
        os.environ["FAKE_VERDICTS"] = json.dumps({"999.6": "fail"})
        self.species_check_refused("1 recorded species-check results changed", [("999.6", "pass", "fail")])

    def test_a_recorded_fail_that_no_longer_fails_stops_the_run(self):
        check = self.master / "species-check/species_check.json"
        doc = json.loads(check.read_text())
        doc["genomes"] = [recorded_genome(g, "fail" if g == "999.6" else "pass") for g in F.COHORT]
        check.write_text(json.dumps(doc))
        os.environ["FAKE_VERDICTS"] = json.dumps({"999.6": "borderline"})
        self.species_check_refused("1 recorded species-check results changed", [("999.6", "fail", "borderline")])

    def test_a_recorded_flag_that_moves_with_another_distance_stops_the_run(self):
        os.environ["FAKE_VERDICTS"] = json.dumps({"999.6": "borderline"})
        os.environ["FAKE_DISTANCES"] = json.dumps({"999.6": 0.06})
        self.species_check_refused("1 recorded species-check results changed", [("999.6", "pass", "borderline")])

    def test_a_recorded_genome_the_new_check_leaves_out_stops_the_run(self):
        os.environ["FAKE_UNJUDGED"] = "999.6"
        self.species_check_refused("1 recorded species-check results changed", [("999.6", "pass", "")])

    def test_a_changed_reference_stops_the_run(self):
        os.environ["FAKE_REFERENCE"] = "999.2"
        self.species_check_refused("1 recorded species-check results changed",
                                   [("999.2", "reference 999.1", "reference 999.2")])

    def test_an_earlier_refusal_leaves_no_stale_list(self):
        os.environ["FAKE_VERDICTS"] = json.dumps({"999.6": "fail"})
        self.species_check_refused("1 recorded species-check results changed", [("999.6", "pass", "fail")])
        os.environ["FAKE_VERDICTS"] = json.dumps({"999.6": "borderline"})
        self.run_ok("species-check")
        self.assertFalse((self.work / "species-check/cohort-changes.tsv").exists())
        self.assertTrue((self.work / "species-check/cohort-gc-moves.tsv").exists())

    def test_cohort_genomes_the_replay_needs_must_all_count(self):
        os.environ["FAKE_UNCOUNTED"] = "999.4"
        self.run_ok("setup", "fetch", "species-check")
        self.assertRefused(self.run_stage("count"), "cohort genomes for the replay checks were not counted")


class KnownGeneAndPredict(RunCase):
    def setUp(self):
        super().setUp()
        os.environ.update(WORLD_RUN)

    def test_known_gene_calls(self):
        self.run_ok("setup", "fetch", "species-check", "count", "known-gene")
        ds = {r["dataset"]: r for r in self.R.S.read_tsv(self.work / "known-gene/datasets.tsv")}
        self.assertEqual((ds[F.A]["replay"], ds[F.B]["replay"]), ("pass", "no rule published"))
        calls = self.R.S.read_tsv(self.work / "known-gene/calls.tsv")
        self.assertEqual(len(calls), len(self.w["survivors"]))
        got = {r["genome_id"]: r["catalog_resistant"] for r in calls}
        want = {g: r["catalog_resistant"] for (d, g), r in self.w["known_gene"]["calls"].items()}
        self.assertEqual(got, want)

    def test_gene_calls_from_another_catalog_are_refused(self):
        os.environ["FAKE_CALLS_VERSION"] = "2026-01-01.1|group-calls-v2"
        self.run_ok("setup", "fetch", "species-check", "count")
        self.assertRefused(self.run_stage("known-gene"), "gene calls from catalog '2026-01-01.1|group-calls-v2'")

    def test_a_report_on_another_catalog_is_refused(self):
        (self.master / "catalog/version.txt").write_text("2026-09-01.1\n")
        self.run_ok("setup", "fetch", "species-check", "count")
        self.assertRefused(self.run_stage("known-gene"), "the report's catalog is 2026-08-07.1")

    def test_a_known_gene_replay_mismatch_stops_the_run(self):
        p = Path(os.environ["FAKE_GENES"])
        p.write_text(json.dumps(dict(json.loads(p.read_text()), **{"999.4": [F.BLAZ]})))
        self.run_ok("setup", "fetch", "species-check", "count")
        self.assertRefused(self.run_stage("known-gene"), "the known-gene replay failed on the cohort")
        self.assertIn("999.4", (self.work / "known-gene/replay-problems.txt").read_text())
        self.assertFalse((self.work / "receipts/known-gene.json").exists())
        self.assertFalse((self.work / "known-gene/calls.tsv").exists())

    def test_predict_replays_first_then_scores_the_survivors(self):
        self.run_ok("all")
        r = self.receipt("predict")
        self.assertEqual(r["counts"]["genomes scored"], len(self.w["survivors"]))
        self.assertEqual(r["counts"]["cohort genome-report replays"], 16)
        self.assertEqual((self.work / "predict/survivors.txt").read_text().split(), self.w["survivors"])
        for d in (F.A, F.B):
            rows = self.R.S.read_rows(self.work / "predict" / f"{d}.jsonl")
            self.assertEqual(len(rows), len(self.w["survivors"]) * len(F.LEARNERS[d]))
            self.assertIn(f"work/predict/{d}.jsonl", r["outputs"])
            self.assertIn(f"work/predict-replay/{d}.jsonl", r["outputs"])

    def assertReplayStops(self, pattern):
        self.run_ok("setup", "fetch", "species-check", "count", "known-gene")
        self.assertRefused(self.run_stage("predict"), pattern)
        self.assertTrue((self.work / "predict-replay/replay-problems.txt").exists())
        self.assertFalse((self.work / "predict").exists())
        self.assertFalse((self.work / "receipts/predict.json").exists())

    def test_a_replay_mismatch_stops_before_any_new_genome_is_predicted(self):
        os.environ["FAKE_PREDICT_WRONG"] = "999.1"
        self.assertReplayStops("the final models do not replay on the cohort")
        self.assertIn("999.1", (self.work / "predict-replay/replay-problems.txt").read_text())

    def test_rows_with_another_k_stop_the_run(self):
        os.environ["FAKE_PREDICT_K"] = "6"
        self.assertReplayStops("the final models do not replay on the cohort")
        self.assertIn("another k", (self.work / "predict-replay/replay-problems.txt").read_text())

    def score(self, out):
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            code = self.kit.S.main(["--work", str(self.work), "--reports", str(self.reports), "--species", F.SPECIES,
                                    "--population", "primary", "--expect-commit", self.kit.commit,
                                    "--out", str(out)])
        return code, err.getvalue()

    def test_a_refused_rerun_leaves_a_completed_run_as_it_was(self):
        # After the last stage ended, one more try of a stage is refused.
        # It must not add a line to the log: score.py refuses a log that
        # goes on after the last stage's end.
        self.run_ok("all")
        log = self.work / "receipts/attempts.jsonl"
        before = log.read_bytes()
        for stage in ("predict", "all", "known-gene"):
            self.assertRefused(self.run_stage(stage), "this stage already has a receipt")
        self.assertEqual(log.read_bytes(), before)
        code, err = self.score(self.root / "out")
        self.assertEqual(code, 0, err)
        res = json.loads((self.root / "out/results.json").read_text())
        self.assertEqual(res["attempts"], {s: 1 for s in self.R.STAGES})

    def test_the_run_scores_end_to_end(self):
        self.run_ok("all")
        out = self.root / "out"
        code, err = self.score(out)
        self.assertEqual(code, 0, err)
        res = json.loads((out / "results.json").read_text())
        want = self.kit.S.score(F.inputs(self.w))
        for key in ("metrics", "paired", "groups", "leave_one_study_out", "dataset_counts", "duplicates",
                    "duplicate_counts"):
            self.assertEqual(res[key], json.loads(json.dumps(want[key])), key)
        self.assertEqual(res["attempts"], {s: 1 for s in self.R.STAGES})
        self.assertTrue(res["table_pins_checked"])
        self.assertEqual(res["kit"]["commit"], self.kit.commit)
        self.assertEqual(res["grm_sha256"], self.R.S.sha256_file(self.grm))


if __name__ == "__main__":
    unittest.main()
