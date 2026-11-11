"""The run script: its receipts, its refusals and the stages that do not
need the real grm, driven by a stand-in (fake_grm.py)."""

import contextlib
import hashlib
import io
import json
import os
import shutil
import stat
import sys
import tempfile
import unittest
from pathlib import Path

from tests.outside_study import fixtures as F

R = F.load("run")
NEW = ["1001.1", "1001.2", "1001.3", "1001.4"]


def tree_digest(root):
    h = hashlib.sha256()
    for p in sorted(Path(root).rglob("*")):
        h.update(str(p.relative_to(root)).encode())
        if p.is_file():
            h.update(p.read_bytes())
    return h.hexdigest()


FAKE_SWITCHES = ("FAKE_MISSING", "FAKE_VERDICTS", "FAKE_UNCOUNTED", "FAKE_DISTANCES", "FAKE_UNJUDGED",
                 "FAKE_REFERENCE")


def recorded_genome(g, verdict="pass"):
    """A genome of the master root's recorded species check, with the fields
    `grm species-check` writes and the stand-in's default values."""
    return {"genome_id": g, "species": F.NAME, "verdict": verdict, "gc": 0.33, "gc_z": 0.0, "distance": 0.01,
            "ani": 0.99, "nearest_species": None, "nearest_distance": None, "inside_species": None,
            "reason": "far from the reference" if verdict == "fail" else None}


class RunScript(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = self.root = Path(self.tmp.name)
        self.paths = F.write_world(F.world(), root / "world")
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
        self.ids = root / "ids.txt"
        self.ids.write_text("\n".join(NEW) + "\n")
        self.datasets = root / "datasets.txt"
        self.datasets.write_text(f"{F.A}\n{F.B}\n")
        self.work = root / "work"
        self.kit = F.frozen_kit(root / "kit")
        self.env = {k: os.environ.get(k) for k in FAKE_SWITCHES}
        for k in FAKE_SWITCHES:
            os.environ.pop(k, None)

    def tearDown(self):
        for k, v in self.env.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v
        self.tmp.cleanup()

    def run_stage(self, stage, *extra, kit=None, work=None, ids=None):
        argv = [stage, "--grm", str(self.grm), "--master", str(self.master), "--work", str(work or self.work),
                "--species", F.SPECIES, "--ids", str(ids or self.ids), "--reports", str(self.paths["reports"]),
                "--datasets", str(self.datasets), "--jobs", "2", *extra]
        err = io.StringIO()
        with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
            code = R.main(argv, overrides={"pins": self.paths["pins"], "kit_root": kit or self.kit})
        return code, err.getvalue()

    def assertRefused(self, result, pattern):
        code, err = result
        self.assertEqual(code, 2, err)
        self.assertRegex(err, pattern)

    def receipt(self, stage):
        return json.loads((self.work / "receipts" / f"{stage}.json").read_text())

    # ---------------------------------------------------------------- setup

    def test_setup_builds_the_scratch_root_and_never_writes_the_master(self):
        before = tree_digest(self.master)
        code, err = self.run_stage("setup")
        self.assertEqual(code, 0, err)
        self.assertEqual(tree_digest(self.master), before)
        r = self.receipt("setup")
        self.assertEqual(r["counts"]["frozen ids requested"], 4)
        self.assertEqual(r["counts"]["cohort genomes for the replay checks"], 8)
        self.assertEqual(r["counts"]["catalog files copied"], 3)
        self.assertEqual(r["grm"]["sha256"], R.S.sha256_file(self.grm))
        self.assertTrue(all(len(v) == 64 for v in r["inputs"].values()))
        root = self.work / "root"
        self.assertTrue((root / "fasta" / "999.1.fna.gz").is_symlink())
        self.assertEqual((root / f"species/outside__{F.SPECIES}/genomes.txt").read_text().split(), NEW)
        self.assertEqual((root / f"species/replay__{F.SPECIES}/genomes.txt").read_text().split(), F.COHORT)

    def test_a_cohort_genome_in_the_ids_is_refused(self):
        self.ids.write_text("1001.1\n999.3\n")
        self.assertRefused(self.run_stage("setup"), "999.3 is a cohort genome")

    def test_a_rehearsal_takes_cohort_genomes_only(self):
        self.assertRefused(self.run_stage("setup", "--rehearsal"), "1001.1 is not a cohort genome")

    def test_an_id_already_in_the_master_fasta_is_refused(self):
        (self.master / "fasta" / "1001.2.fna.gz").write_bytes(b"x")
        self.assertRefused(self.run_stage("setup"), "1001.2 already has a FASTA file in the master root")

    def test_the_work_directory_inside_the_master_is_refused(self):
        self.assertRefused(self.run_stage("setup", work=self.master / "w"), "inside the master root")

    def test_no_catalog_is_refused_outside_a_rehearsal(self):
        shutil.rmtree(self.master / "catalog")
        self.assertRefused(self.run_stage("setup"), "no catalog")

    def test_a_bad_ids_file_is_refused(self):
        self.ids.write_text("1001.1\n1001.1\n")
        self.assertRefused(self.run_stage("setup"), "listed twice")
        self.ids.write_text("1001.1\nGCF_000001\n")
        self.assertRefused(self.run_stage("setup"), "not a BV-BRC genome id")

    def test_a_dataset_of_another_species_is_refused(self):
        self.datasets.write_text("other_species__drugA\n")
        self.assertRefused(self.run_stage("setup"), "is not a testus_bacterium dataset")

    def test_a_report_not_matching_its_pin_is_refused(self):
        p = self.paths["reports"] / f"{F.A}.json"
        p.write_bytes(p.read_bytes() + b"\n")
        self.assertRefused(self.run_stage("setup"), "is not the pinned")

    # ---------------------------------------------------------------- receipts

    def test_a_stage_without_the_earlier_receipt_is_refused(self):
        self.assertRefused(self.run_stage("fetch"), "the setup stage has no receipt")
        self.assertRefused(self.run_stage("count"), "the setup stage has no receipt")

    def test_a_changed_input_is_refused(self):
        self.assertEqual(self.run_stage("setup")[0], 0)
        self.ids.write_text("\n".join(NEW[:3]) + "\n")
        self.assertRefused(self.run_stage("fetch"), "ids.txt changed since the setup stage")

    def test_a_changed_output_is_refused(self):
        self.assertEqual(self.run_stage("setup")[0], 0)
        self.assertEqual(self.run_stage("fetch")[0], 0)
        av = self.work / "fetch" / "availability.tsv"
        av.write_text(av.read_text().replace("1001.4\tyes", "1001.4\tno"))
        self.assertRefused(self.run_stage("species-check"), "availability.tsv changed since the fetch stage")

    def test_another_grm_binary_is_refused(self):
        self.assertEqual(self.run_stage("setup")[0], 0)
        with open(self.grm, "a") as f:
            f.write("\n# changed\n")
        self.assertRefused(self.run_stage("fetch"), "grm binary differs")

    def test_other_arguments_are_refused(self):
        self.assertEqual(self.run_stage("setup")[0], 0)
        other = self.root / "ids2.txt"
        other.write_text(self.ids.read_text())
        self.assertRefused(self.run_stage("fetch", ids=other), "arguments differ")

    def test_the_gated_stages_refuse_without_the_frozen_files(self):
        self.assertEqual(self.run_stage("setup")[0], 0)
        # Earlier receipts in place (copies of setup's, whose files are
        # unchanged), so only the gate stands in the way.
        for stage in ("fetch", "species-check", "count", "known-gene"):
            shutil.copyfile(self.work / "receipts/setup.json", self.work / f"receipts/{stage}.json")
        nokit = self.root / "nokit"
        nokit.mkdir()
        for stage in ("known-gene", "predict"):
            self.assertRefused(self.run_stage(stage, kit=nokit), "not frozen")
        (self.kit / "docs/outside-study-test.md").write_text("changed plan\n")
        for stage in ("known-gene", "predict"):
            self.assertRefused(self.run_stage(stage), "outside-study-test.md has sha256")

    # ---------------------------------------------------------------- stages

    def test_fetch_species_check_and_count(self):
        os.environ["FAKE_MISSING"] = "1001.2"
        os.environ["FAKE_VERDICTS"] = json.dumps({"1001.3": "fail", "1001.4": "borderline"})
        os.environ["FAKE_UNCOUNTED"] = ""
        for stage in ("setup", "fetch", "species-check", "count"):
            code, err = self.run_stage(stage)
            self.assertEqual(code, 0, f"{stage}: {err}")
        av = {r["genome_id"]: r for r in R.S.read_tsv(self.work / "fetch/availability.tsv")}
        self.assertEqual({g: r["fasta"] for g, r in av.items()},
                         {"1001.1": "yes", "1001.2": "no", "1001.3": "yes", "1001.4": "yes"})
        self.assertEqual(av["1001.2"]["detail"], "no contigs returned")
        v = {r["genome_id"]: r["verdict"] for r in R.S.read_tsv(self.work / "species-check/verdicts.tsv")}
        self.assertEqual(v, {"1001.1": "pass", "1001.3": "fail", "1001.4": "borderline"})
        c = {r["genome_id"]: r["counted"] for r in R.S.read_tsv(self.work / "count/counted.tsv")}
        self.assertEqual(c["1001.2"], "no")
        self.assertEqual(self.receipt("fetch")["counts"]["fetched"], 3)
        self.assertEqual(self.receipt("species-check")["counts"]["fail (excluded)"], 1)
        self.assertEqual(self.receipt("count")["counts"][f"replay__{F.SPECIES}: counted"], 8)
        survivors, excluded = R.Run(self.args()).survivors()
        self.assertEqual(survivors, ["1001.1", "1001.4"])
        self.assertEqual({g: e[0] for g, e in excluded.items()},
                         {"1001.2": "no sequence", "1001.3": "species check failed"})

    def args(self):
        import argparse
        return argparse.Namespace(grm=str(self.grm), master=str(self.master), work=str(self.work),
                                  species=F.SPECIES, ids=str(self.ids), reports=str(self.paths["reports"]),
                                  datasets=str(self.datasets), jobs=2, rehearsal=False,
                                  pins=self.paths["pins"], kit_root=self.kit)

    def species_check_refused(self, pattern, rows):
        self.assertEqual(self.run_stage("setup")[0], 0)
        self.assertEqual(self.run_stage("fetch")[0], 0)
        self.assertRefused(self.run_stage("species-check"), pattern)
        changes = R.S.read_tsv(self.work / "species-check/cohort-changes.tsv")
        self.assertEqual([(r["genome_id"], r["recorded"], r["now"]) for r in changes], rows)
        self.assertFalse((self.work / "receipts/species-check.json").exists())
        self.assertFalse((self.work / "species-check/cohort-gc-moves.tsv").exists())
        self.assertFalse((self.work / "species-check/verdicts.tsv").exists())

    def test_a_recorded_gc_flag_that_moves_is_listed_and_the_run_goes_on(self):
        os.environ["FAKE_VERDICTS"] = json.dumps({"999.6": "borderline"})
        for stage in ("setup", "fetch", "species-check"):
            code, err = self.run_stage(stage)
            self.assertEqual(code, 0, f"{stage}: {err}")
        moves = R.S.read_tsv(self.work / "species-check/cohort-gc-moves.tsv")
        self.assertEqual([(r["genome_id"], r["recorded"], r["now"]) for r in moves],
                         [("999.6", "pass", "borderline")])
        self.assertFalse((self.work / "species-check/cohort-changes.tsv").exists())
        counts = self.receipt("species-check")["counts"]
        self.assertEqual(counts["recorded genomes with unchanged verdicts"], len(F.COHORT) - 1)
        self.assertEqual(counts["recorded genomes whose GC flag alone moved (listed)"], 1)
        self.assertIn(str(self.work / "species-check/cohort-gc-moves.tsv"), self.receipt("species-check")["outputs"])

    def test_no_recorded_flag_moving_writes_an_empty_list(self):
        for stage in ("setup", "fetch", "species-check"):
            self.assertEqual(self.run_stage(stage)[0], 0)
        self.assertEqual(R.S.read_tsv(self.work / "species-check/cohort-gc-moves.tsv"), [])
        counts = self.receipt("species-check")["counts"]
        self.assertEqual(counts["recorded genomes with unchanged verdicts"], len(F.COHORT))
        self.assertEqual(counts["recorded genomes whose GC flag alone moved (listed)"], 0)

    def test_a_recorded_borderline_that_moves_back_to_pass_is_listed(self):
        check = self.master / "species-check/species_check.json"
        doc = json.loads(check.read_text())
        doc["genomes"] = [recorded_genome(g, "borderline" if g == "999.6" else "pass") for g in F.COHORT]
        check.write_text(json.dumps(doc))
        for stage in ("setup", "fetch", "species-check"):
            self.assertEqual(self.run_stage(stage)[0], 0)
        moves = R.S.read_tsv(self.work / "species-check/cohort-gc-moves.tsv")
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
        code, err = self.run_stage("species-check")
        self.assertEqual(code, 0, err)
        self.assertFalse((self.work / "species-check/cohort-changes.tsv").exists())
        self.assertTrue((self.work / "species-check/cohort-gc-moves.tsv").exists())

    def test_cohort_genomes_the_replay_needs_must_all_count(self):
        os.environ["FAKE_UNCOUNTED"] = "999.4"
        for stage in ("setup", "fetch", "species-check"):
            self.assertEqual(self.run_stage(stage)[0], 0)
        self.assertRefused(self.run_stage("count"), "cohort genomes for the replay checks were not counted")

    def test_no_option_moves_the_pins_or_the_kit_root(self):
        for flag in ("--pins", "--kit-root"):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                R.main(["setup", "--grm", "g", "--master", "m", "--work", "w", "--species", "s", "--ids", "i",
                        "--reports", "r", "--datasets", "d", flag, "x"])


if __name__ == "__main__":
    unittest.main()
