"""Tests for what H-1 needs of scripts/input-comparison/run.py (H1PLAN's gaps): kset sharing (G1), `run.py combine`
(G2), resume (G6) and the H-1 spec files (G10). With GRM_INPUT_E2E_BINARY, a second synthetic dataset runs the
comparison whole and in two dataset parts through the real `grm`, and combine must give the whole run's J and K.

usage: python3 -B -m unittest discover -s tests/input_comparison -t .
"""
import argparse
import contextlib
import copy
import hashlib
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests.input_comparison import test_run as base

run = base.run
HERE = base.HERE
SPEC = base.SPEC
DS = "salmonella_enterica__meropenem"
# DECISIONS D2 after H1PLAN: revision 12's measured kset bytes per genome, rounded up.
D2 = {"salmonella_enterica": 29147870, "mycobacterium_tuberculosis": 25881860, "staphylococcus_aureus": 16812689,
      "neisseria_gonorrhoeae": 12750902}

# base.SpecTests' stub grm with two changes: `grm arm` writes a complete arm output (its one plan tagged with
# STUB_PLAN_TAG; with STUB_ARM_FAILS set it fails before writing it) and records each call in WORK/arm-calls.txt;
# with STUB_RECOUNT set, `grm count` replaces every kset already in the root, as a recount would.
STUB_GRM = base.SpecTests.STUB_GRM.replace(
    """    (out / "seen.json").write_text(json.dumps({"plans_present": (out / "plans").exists()}))
    sys.exit(0)""",
    """    (out / "seen.json").write_text(json.dumps({"plans_present": (out / "plans").exists()}))
    if os.environ.get("STUB_ARM_FAILS"):
        sys.exit(5)
    (out / "plans").mkdir(exist_ok=True)
    (out / "plans" / "lineage-0.json").write_text(json.dumps({"tag": os.environ.get("STUB_PLAN_TAG", "")}))
    (out / "runs.json").write_text("[]")
    (out / "arm.json").write_text("{}")
    (out / "predictions.jsonl").write_text("")
    with open(root.parents[1] / "arm-calls.txt", "a") as f:
        f.write(argv[argv.index("--arm") + 1] + "\\n")
    sys.exit(0)""").replace(
    """    sys.exit(3 if os.environ.get("STUB_COUNT_FAILS") else 0)""",
    """    if argv[0] == "count" and os.environ.get("STUB_RECOUNT"):
        for kset in list((root / "species" / "salmonella_enterica" / "ksets").glob("*.kset")):
            (kset.parent / "new.part").write_bytes(kset.read_bytes())
            os.replace(kset.parent / "new.part", kset)
    sys.exit(3 if os.environ.get("STUB_COUNT_FAILS") else 0)""")


class Stub(unittest.TestCase):
    """A WORK curated by the stub grm (base.SpecTests' spec: default cap 20, cap-10, cap-30 and k21; 50 eligible
    genomes g000..g049), and a Runner on it."""

    def setUp(self):
        self.assertNotEqual(STUB_GRM, base.SpecTests.STUB_GRM)
        self.assertIn("arm-calls.txt", STUB_GRM)
        self.assertIn("STUB_RECOUNT", STUB_GRM)
        patcher = mock.patch.object(run.Runner, "check_synthetic", lambda runner: None)
        patcher.start()
        self.addCleanup(patcher.stop)
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmp = Path(tmp.name)
        self.grm = self.tmp / "grm"
        self.grm.write_text(f"#!{sys.executable}\n" + STUB_GRM)
        self.grm.chmod(0o755)

    def spec_path(self, **changes):
        return base.SpecTests.stub_spec(self, self.tmp, **changes)

    def estimate(self, *flags, spec=None):
        out = io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()), \
                self.assertRaises(SystemExit) as e:
            run.main([str(spec or self.spec_path()), "--raw", str(self.tmp / "raw"), "--work", str(self.tmp / "work"),
                      "--grm", str(self.grm), "--source-commit", "0" * 40, "--below-panel-floor", "--estimate-only",
                      "--free-disk-bytes", str(10 ** 15), *flags])
        if isinstance(e.exception.code, str):
            raise AssertionError(e.exception.code)
        return json.loads(out.getvalue())

    def runner(self, *flags):
        """A Runner on the curated stub WORK, as run.main sets it up."""
        spec, default, arms = run.load_spec(self.spec_path())
        args = argparse.Namespace(grm=str(self.grm), raw=str(self.tmp / "raw"), work=str(self.tmp / "work"),
                                  no_share_ksets="--no-share-ksets" in flags,
                                  trace_training="--trace-training" in flags)
        r = run.Runner(args, spec)
        r.default, r.all_arms = default, [default, *arms]
        r.arms = {a["name"]: a for a in arms}
        return r

    def fake_pool(self, r, ids):
        """Pool ksets for `ids` (a kset and its sidecar each), as if the pool were counted; returns their folder."""
        union = r.union(r.all_arms)
        ksets = r.kset_dir(r.pool_root(union))
        ksets.mkdir(parents=True)
        for g in ids:
            (ksets / f"{g}.kset").write_bytes(f"kset {g}".encode())
            (ksets / f"{g}.json").write_text(json.dumps({"genome_id": g}))
        (ksets.parent / "structure.json").write_text(json.dumps({"genomes": union, "k": 31}))
        r.pool_ksets = ksets
        return ksets


class KsetSharingTests(Stub):
    """G1: the k31 roots link the pool's ksets; the disk check and disk_bytes count each once."""

    def test_disk_bytes_counts_a_hard_linked_file_once(self):
        a = self.tmp / "d" / "a"
        a.parent.mkdir()
        a.write_bytes(b"x" * 100000)
        one = run.disk_bytes(self.tmp / "d")
        os.link(a, self.tmp / "d" / "b")
        self.assertEqual(run.disk_bytes(self.tmp / "d"), one)
        shutil.copy(a, self.tmp / "d" / "c")
        self.assertEqual(run.disk_bytes(self.tmp / "d"), 2 * one)

    def test_a_k31_root_takes_the_pool_s_ksets_and_a_k21_root_its_own(self):
        self.estimate()
        r = self.runner()
        pool = self.fake_pool(r, [f"g{i:03d}" for i in range(30)])
        default = r.root(r.default)
        mine = r.kset_dir(default)
        ids = r.genomes_txt(default).read_text().split()
        self.assertEqual(len(ids), 20)
        for g in ids:
            for ext in ("kset", "json"):
                self.assertTrue(os.path.samefile(mine / f"{g}.{ext}", pool / f"{g}.{ext}"), (g, ext))
        self.assertEqual(json.loads((default / run.SHARED).read_text()),
                         {"genomes": 20, "shared": 20, "own": 0, "linked_now": 20, "pool": pool.parents[2].name})
        self.assertTrue((default / "built.json").exists())
        k21 = r.root(r.arms["k21"])
        self.assertFalse(list(r.kset_dir(k21).glob("*.kset")))
        self.assertEqual(json.loads((k21 / run.SHARED).read_text())["shared"], 0)
        # A built root is not linked again, and a root's own kset is kept.
        cap10 = r.curated(r.arms["cap-10"])
        r.kset_dir(cap10).mkdir(parents=True)
        (r.kset_dir(cap10) / "g000.kset").write_bytes(b"own")
        r.root(r.arms["cap-10"])
        self.assertEqual((r.kset_dir(cap10) / "g000.kset").read_bytes(), b"own")
        self.assertFalse((r.kset_dir(cap10) / "g000.json").exists())
        self.assertEqual(json.loads((cap10 / run.SHARED).read_text())["shared"], 9)

    def test_without_sharing_every_root_counts_its_own(self):
        self.estimate("--no-share-ksets")
        r = self.runner("--no-share-ksets")
        self.fake_pool(r, [f"g{i:03d}" for i in range(30)])
        default = r.root(r.default)
        self.assertFalse(list(r.kset_dir(default).glob("*")))
        self.assertEqual(json.loads((default / run.SHARED).read_text())["own"], 20)

    def test_a_recount_of_a_linked_kset_is_refused(self):
        self.estimate()
        r = self.runner()
        self.fake_pool(r, [f"g{i:03d}" for i in range(30)])
        os.environ["STUB_RECOUNT"] = "1"
        try:
            with self.assertRaises(SystemExit) as e:
                r.root(r.default)
        finally:
            del os.environ["STUB_RECOUNT"]
        self.assertIn("grm count rewrote 20 of the 20 ksets linked from the species pool", str(e.exception))
        self.assertFalse((r.curated(r.default) / "built.json").exists())

    def test_the_disk_check_counts_a_shared_kset_once(self):
        per = SPEC["pool_budget"]["default_kset_bytes"]
        shared = self.estimate()
        shutil.rmtree(self.tmp / "work")
        own = self.estimate("--no-share-ksets")
        k31 = {n: r for n, r in shared["roots"].items() if r["k"] == 31}
        # Roots: the default's (20), cap-10's and cap-30's at k 31, k21's.
        self.assertEqual(sorted(r["genomes"] for r in k31.values()), [10, 20, 30])
        for n, r in k31.items():
            self.assertEqual((r["pending"], r["shared"]), (0, r["genomes"]), n)
            self.assertEqual((own["roots"][n]["pending"], own["roots"][n]["shared"]), (r["genomes"], 0), n)
        (k21,) = [r for r in shared["roots"].values() if r["k"] == 21]
        self.assertEqual((k21["pending"], k21["shared"]), (20, 0))
        self.assertEqual(own["kset_bytes"] - shared["kset_bytes"], round(60 * per))
        self.assertEqual(shared["pool_bytes"], own["pool_bytes"])
        # With the pool counted, only the pool's finished ksets are shared.
        r = self.runner()
        self.fake_pool(r, [f"g{i:03d}" for i in range(25)])
        roots = r.counted_roots(r.all_arms)
        (cap30,) = [v for v in roots.values() if v["genomes"] == 30]
        self.assertEqual((cap30["shared"], cap30["pending"]), (25, 5))

    def test_estimate_only_reports_the_union_s_missing_fastas(self):
        disk = self.estimate()
        self.assertEqual(disk["union_without_fasta"], 30)
        (self.tmp / "raw" / "fasta").mkdir(parents=True, exist_ok=True)
        for i in range(29):
            (self.tmp / "raw" / "fasta" / f"g{i:03d}.fna.gz").write_bytes(b"")
        shutil.rmtree(self.tmp / "work")
        disk = self.estimate()
        self.assertEqual((disk["union_without_fasta"], disk["union_without_fasta_examples"]), (1, ["g029"]))


class ResumeTests(Stub):
    """G6: a reused WORK skips an arm whose outputs are complete and were made from the same inputs, spec, grouping
    and reference plans; anything else reruns it."""

    def setUp(self):
        super().setUp()
        self.estimate()
        self.groups = self.tmp / "structure.json"
        self.groups.write_text(json.dumps({"genomes": [], "k": 31}))

    def run_arms(self, *flags, spec=None):
        """Runs the default and cap-10 on DS; returns the arms grm ran this time."""
        r = self.runner(*flags)
        if spec is not None:
            r.spec = spec
        r.groups = self.groups
        calls = self.tmp / "work" / "arm-calls.txt"
        before = calls.read_text().split() if calls.exists() else []
        r.arm(r.default, DS)
        r.arm(r.arms["cap-10"], DS)
        after = calls.read_text().split()
        return after[len(before):], r

    def out(self, arm):
        return self.tmp / "work" / "arms" / arm / DS

    def test_a_finished_arm_is_resumed_and_records_its_outputs(self):
        ran, _ = self.run_arms()
        self.assertEqual(ran, ["default", "cap-10"])
        record = json.loads((self.out("cap-10") / run.RESUME).read_text())
        self.assertEqual(set(record["outputs"]),
                         {"arm.json", "runs.json", "predictions.jsonl", "seen.json", "plans/lineage-0.json"})
        self.assertEqual(record["reference_plans"],
                         {"lineage-0.json": run.file_sha256(self.out("default") / "plans" / "lineage-0.json")})
        self.assertEqual(record["inputs_sha256"], run.file_sha256(self.tmp / "work" / run.INPUTS))
        ran, r = self.run_arms()
        self.assertEqual((ran, sorted(r.resumed)), ([], [f"cap-10/{DS}", f"default/{DS}"]))

    def test_a_changed_input_reruns_the_arm(self):
        self.run_arms()
        spec = copy.deepcopy(json.loads(self.spec_path().read_text()))
        spec["margin"] = {"log_loss": 0.03}
        self.assertEqual(self.run_arms(spec=spec)[0], ["default", "cap-10"], "spec")
        self.assertEqual(self.run_arms()[0], ["default", "cap-10"], "spec back")
        self.assertEqual(self.run_arms("--trace-training")[0], ["default", "cap-10"], "trace")
        self.assertEqual(self.run_arms()[0], ["default", "cap-10"])
        self.groups.write_text(json.dumps({"genomes": ["g000"], "k": 31}))
        self.assertEqual(self.run_arms()[0], ["default", "cap-10"], "grouping")
        record = self.tmp / "work" / run.INPUTS
        record.write_text(record.read_text().replace('"raw_files":', '"raw_files": 0, "was_raw_files":'))
        self.assertEqual(self.run_arms()[0], ["default", "cap-10"], "inputs.json")
        # The default's plans change: cap-10's reference plans no longer match.
        os.environ["STUB_PLAN_TAG"] = "other"
        try:
            (self.out("default") / run.RESUME).unlink()
            self.assertEqual(self.run_arms()[0], ["default", "cap-10"], "reference plans")
        finally:
            del os.environ["STUB_PLAN_TAG"]
        # Nothing changed since: both resume, cap-10 on the default's new plans.
        self.assertEqual(self.run_arms()[0], [])

    def test_an_incomplete_or_changed_output_reruns_the_arm(self):
        self.run_arms()
        out = self.out("cap-10")
        for name, damage in (("predictions.jsonl", lambda p: p.write_text("{}\n")),
                             ("arm.json", lambda p: p.unlink()),
                             ("extra.json", lambda p: p.write_text("{}")),
                             (run.RESUME, lambda p: p.unlink()),
                             (run.RESUME, lambda p: p.write_text("{")),
                             ("plans/lineage-0.json", lambda p: p.write_text("{}"))):
            damage(out / name)
            ran, _ = self.run_arms()
            self.assertEqual(ran, ["cap-10"], name)
            (out / "extra.json").unlink(missing_ok=True)
            self.run_arms()
        # A record without a required output never resumes.
        record = json.loads((out / run.RESUME).read_text())
        del record["outputs"]["runs.json"]
        (out / "runs.json").unlink()
        (out / run.RESUME).write_text(run.canonical(record))
        self.assertEqual(self.run_arms()[0], ["cap-10"])

    def test_only_a_finished_arm_has_a_record(self):
        """A rerun removes the old resume.json first, so an arm whose rerun failed has none until it finishes."""
        self.run_arms()
        out = self.out("cap-10")
        (out / "predictions.jsonl").write_text("{}\n")
        os.environ["STUB_ARM_FAILS"] = "1"
        try:
            with self.assertRaises(RuntimeError):
                self.run_arms()
        finally:
            del os.environ["STUB_ARM_FAILS"]
        self.assertFalse((out / run.RESUME).exists())
        self.assertTrue((self.out("default") / run.RESUME).exists())
        self.assertEqual(self.run_arms()[0], ["cap-10"])
        self.assertTrue((out / run.RESUME).exists())

    def test_a_linked_output_is_never_resumed(self):
        self.run_arms()
        out = self.out("cap-10")
        target = self.tmp / "elsewhere.jsonl"
        shutil.move(out / "predictions.jsonl", target)
        (out / "predictions.jsonl").symlink_to(target)
        self.assertEqual(self.run_arms()[0], ["cap-10"])


def fake_part(top, name, spec, datasets, raw=None, **changes):
    """A finished part's WORK as `combine` reads it, with no arm outputs: spec.json (the spec with `datasets`),
    inputs.json, raw-manifest.txt, K's manifest, species-groups.json and summary.json. `changes` replaces one of
    inputs / manifest / groups / summary, or the spec (`spec_changes`)."""
    work = Path(top) / name
    (work / "k").mkdir(parents=True)
    s = {**spec, "datasets": [{"id": d, "partition": "confirmatory"} for d in datasets], **changes.get("spec_changes", {})}
    lines = raw if raw is not None else [f"{'a' * 64}  3  bvbrc/genomes.jsonl", f"{'b' * 64}  5  fasta/g1.fna.gz"]
    text = "\n".join(lines)
    inputs = {"raw_manifest_sha256": run.sha256_bytes(text.encode()), "raw_files": len(lines), "grm_sha256": "c" * 64,
              "grm": "/x/grm", "vocabulary_sha256": "d" * 64, **changes.get("inputs", {})}
    manifest = {"protocol_sha256": run.sha256_bytes(run.canonical(s).encode()), "source_commit": "e" * 40,
                "binary_sha256": inputs["grm_sha256"], "synthetic": True, "revision12_gate": "not_required",
                "compute": {"actual_seconds": 1.0}, **changes.get("manifest", {})}
    groups = {"structure": f"{work}/roots/pool/structure.json", "sha256": "f" * 64, "genomes": 30, "lineages": 3,
              "arms": ["default"], **changes.get("groups", {})}
    summary = {"eligible_genomes": 50, "panel_floor": 2001, "seconds": 1.0, **changes.get("summary", {})}
    for path, value in ((run.SPEC_COPY, run.canonical(s)), (run.INPUTS, run.canonical(inputs)),
                        ("k/manifest.json", json.dumps(manifest)), ("species-groups.json", run.canonical(groups)),
                        ("summary.json", json.dumps(summary)), (run.RAW_MANIFEST, text)):
        (work / path).write_text(value)
    return work


class CombineRefusalTests(unittest.TestCase):
    """G2: `combine` refuses parts that are not one spec split by datasets on one pool, grm, vocabulary, source
    commit and master manifest. A refusal comes before any arm output is read; valid parts get as far as that."""

    A, B = DS, "salmonella_enterica__ciprofloxacin"

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmp = Path(tmp.name)
        self.spec = {**copy.deepcopy(SPEC), "datasets": [{"id": self.A, "partition": "confirmatory"},
                                                         {"id": self.B, "partition": "confirmatory"}]}
        (self.tmp / "spec.json").write_text(json.dumps(self.spec))
        self.master = self.tmp / "MASTER.sha256"
        self.master.write_text(f"{'a' * 64}  bvbrc/genomes.jsonl\n{'b' * 64}  ./fasta/g1.fna.gz\n{'0' * 64}  other\n")

    def combine(self, *works, out=None):
        argv = ["combine", str(self.tmp / "spec.json"), "--master-manifest", str(self.master),
                "--out", str(out or self.tmp / f"out-{len(list(self.tmp.glob('out-*')))}")]
        for w in works:
            argv += ["--work-from", str(w)]
        with self.assertRaises(SystemExit) as e, contextlib.redirect_stdout(io.StringIO()):
            run.main(argv)
        return str(e.exception.code)

    def part(self, name, datasets, **changes):
        return fake_part(self.tmp, name, self.spec, datasets, **changes)

    def test_valid_parts_pass_every_check(self):
        a, b = self.part("a", [self.A]), self.part("b", [self.B])
        self.assertEqual(self.combine(a, b), f"part {a} has no output of arm default on {self.A}")
        self.assertEqual(self.combine(self.part("both", [self.A, self.B])),
                         f"part {self.tmp / 'both'} has no output of arm default on {self.A}")

    def test_parts_must_be_disjoint_and_cover_the_spec(self):
        a, both = self.part("a", [self.A]), self.part("both", [self.A, self.B])
        self.assertEqual(self.combine(a, both), f"dataset {self.A} is in two parts: {a} and {both}")
        self.assertEqual(self.combine(a), f"no part holds {self.B} of {self.tmp / 'spec.json'}")
        self.assertEqual(self.combine(a, a), f"dataset {self.A} is in two parts: {a} and {a}")
        c = self.part("c", ["salmonella_enterica__ampicillin"])
        self.assertIn("dataset salmonella_enterica__ampicillin with partition confirmatory is not in",
                      self.combine(a, c))
        d = fake_part(self.tmp, "d", self.spec, [], spec_changes={"datasets": [{"id": self.B, "partition": "development"}]})
        self.assertIn(f"dataset {self.B} with partition development is not in", self.combine(a, d))

    def test_a_part_spec_may_differ_only_in_its_datasets(self):
        a = self.part("a", [self.A])
        for changes, key in (({"margin": {"log_loss": 0.03}}, "margin"), ({"seed": 1}, "seed"),
                             ({"arms": SPEC["arms"][:-1]}, "arms"), ({"species": "other"}, "species")):
            b = self.part(f"b-{key}", [self.B], spec_changes=changes)
            self.assertEqual(self.combine(a, b), f"part {b}: its spec differs from {self.tmp / 'spec.json'} in {key}, "
                                                 "not only in its datasets")
        b = self.part("b-k", [self.B], manifest={"protocol_sha256": "0" * 64})
        self.assertEqual(self.combine(a, b), f"part {b}: K's manifest was not made from its spec.json")

    def test_parts_must_share_the_pool_grm_vocabulary_and_commit(self):
        a = self.part("a", [self.A])
        for name, changes, label in (
                ("pool", {"groups": {"sha256": "0" * 64}}, "the species pool"),
                ("lineages", {"groups": {"lineages": 4}}, "the species pool"),
                ("grm", {"inputs": {"grm_sha256": "0" * 64}, "manifest": {"binary_sha256": "0" * 64}},
                 "the grm binary"),
                ("vocabulary", {"inputs": {"vocabulary_sha256": None}}, "the vocabulary"),
                ("commit", {"manifest": {"source_commit": "1" * 40}}, "the source commit"),
                ("gate", {"manifest": {"revision12_gate": "passed"}}, "synthetic and the revision-12 gate"),
                ("eligible", {"summary": {"eligible_genomes": 49}}, "eligible genomes")):
            b = self.part(f"b-{name}", [self.B], **changes)
            self.assertIn(f"the parts differ in {label}", self.combine(a, b), name)
        # The pool's path may differ (another WORK); its content may not.
        b = self.part("b-path", [self.B], groups={"structure": "/elsewhere/structure.json"})
        self.assertIn("has no output of arm default", self.combine(a, b))
        b = self.part("b-binary", [self.B], manifest={"binary_sha256": "0" * 64})
        c = self.part("c-binary", [self.A], manifest={"binary_sha256": "0" * 64})
        self.assertEqual(self.combine(c, b), f"part {c}: K's binary_sha256 is not the grm its inputs.json hashed")

    def test_every_part_s_raw_must_be_in_the_master_manifest(self):
        a = self.part("a", [self.A])
        b = self.part("b", [self.B], raw=[f"{'a' * 64}  3  bvbrc/genomes.jsonl", f"{'9' * 64}  5  fasta/g1.fna.gz",
                                         f"{'b' * 64}  5  fasta/g2.fna.gz"])
        self.assertEqual(self.combine(a, b), f"part {b}: 2 of its 3 RAW files are not in the master manifest "
                                             f"{self.master} with that sha256, e.g. fasta/g1.fna.gz, fasta/g2.fna.gz")
        (b / run.RAW_MANIFEST).write_text(f"{'a' * 64}  3  bvbrc/genomes.jsonl")
        self.assertEqual(self.combine(a, b), f"part {b}: raw-manifest.txt is not the RAW manifest its inputs.json "
                                             "hashed")
        self.master.write_text("not a manifest\n")
        self.assertIn("line 1: not a sha256sum line", self.combine(a, self.part("c", [self.B])))

    def test_out_must_be_new_or_empty_and_parts_finished(self):
        a, b = self.part("a", [self.A]), self.part("b", [self.B])
        (self.tmp / "full" / "x").mkdir(parents=True)
        self.assertEqual(self.combine(a, b, out=self.tmp / "full"), f"--out {self.tmp / 'full'} is not an empty folder")
        (b / "summary.json").unlink()
        self.assertIn(f"part {b}: no readable summary.json", self.combine(a, b))


class H1SpecTests(unittest.TestCase):
    """G10: scripts/input-comparison/h1/ holds the pilot (D1) and the four species' stage-1 specs with the Salmonella
    parts; each is spec.json with only its species, datasets, kset size (D2) and arms (no k arm, M-2) changed."""

    H1 = HERE / "h1"
    COUNTS = {"salmonella_enterica": 25, "mycobacterium_tuberculosis": 22, "staphylococcus_aureus": 18,
              "neisseria_gonorrhoeae": 7}

    def specs(self):
        return {p.name: json.loads(p.read_text()) for p in sorted(self.H1.glob("*.json"))}

    def test_the_folder_holds_the_pilot_the_four_species_and_the_salmonella_parts(self):
        self.assertEqual(sorted(self.specs()), sorted(
            ["pilot.json", "h1-salmonella_enterica-part1.json", "h1-salmonella_enterica-part2.json",
             *(f"h1-{s}.json" for s in self.COUNTS)]))

    def test_each_spec_changes_only_its_species_datasets_kset_size_and_arms(self):
        for name, s in self.specs().items():
            self.assertEqual(set(s), set(SPEC), name)
            for key in set(SPEC) - {"decisions", "species", "datasets", "pool_budget", "arms"}:
                self.assertEqual(s[key], SPEC[key], (name, key))
            self.assertEqual({k: v for k, v in s["pool_budget"].items() if k not in ("default_kset_bytes", "reason")},
                             {k: v for k, v in SPEC["pool_budget"].items()
                              if k not in ("default_kset_bytes", "reason")}, name)
            self.assertEqual(s["pool_budget"]["default_kset_bytes"], D2[s["species"]], name)
            self.assertIn(f"{D2[s['species']]:,} B per genome", s["pool_budget"]["reason"], name)
            self.assertIn("M-2", s["decisions"], name)
            self.assertNotIn("synthetic", s, name)
            for arm in s["arms"]:
                self.assertIn(arm, SPEC["arms"], name)
            _, default, arms = run.load_spec(self.H1 / name)
            self.assertEqual(default["k"], 31)
            self.assertFalse([a["name"] for a in arms if a["factor"] == "k"], name)
            self.assertEqual(arms[-1]["name"], run.SENTINEL["name"])
            for d in s["datasets"]:
                self.assertEqual(d["partition"], "confirmatory", (name, d))
                self.assertTrue(d["id"].startswith(s["species"] + "__"), (name, d))

    def test_stage_1_runs_spec_json_s_arms_without_k21_and_k25(self):
        want = [a for a in SPEC["arms"] if a["name"] not in ("k21", "k25")]
        for species, n in self.COUNTS.items():
            s = self.specs()[f"h1-{species}.json"]
            self.assertEqual((s["species"], len(s["datasets"]), s["arms"]), (species, n, want), species)

    def test_the_salmonella_parts_split_its_datasets_13_and_12(self):
        specs = self.specs()
        full = specs["h1-salmonella_enterica.json"]
        parts = [specs[f"h1-salmonella_enterica-part{i}.json"] for i in (1, 2)]
        self.assertEqual([len(p["datasets"]) for p in parts], [13, 12])
        ids = [d["id"] for p in parts for d in p["datasets"]]
        self.assertEqual(sorted(ids), sorted(d["id"] for d in full["datasets"]))
        self.assertEqual(len(set(ids)), len(ids))
        for p in parts:
            self.assertEqual({k: v for k, v in p.items() if k != "datasets"},
                             {k: v for k, v in full.items() if k != "datasets"})

    def test_the_pilot_is_d1_s(self):
        s = self.specs()["pilot.json"]
        self.assertEqual((s["species"], s["datasets"]),
                         ("salmonella_enterica", [{"id": "salmonella_enterica__ciprofloxacin",
                                                   "partition": "confirmatory"}]))
        _, _, arms = run.load_spec(self.H1 / "pilot.json")
        self.assertEqual([(a["name"], a["factor"]) for a in arms],
                         [("grr-v1", "cohort"), ("k-c-h-p", "blocks"), ("budget-16000", "budget"),
                          ("sentinel-seed-twin", "calibration")])
        full = self.specs()["h1-salmonella_enterica.json"]
        self.assertEqual(s["pool_budget"], full["pool_budget"])


class SyntheticRootTests(unittest.TestCase):
    def test_the_second_drug_adds_its_records_and_changes_nothing_else(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name, flags in (("one", []), ("two", ["--second-drug"])):
                subprocess.run([sys.executable, "-I", "-B", str(HERE / "synthetic_root.py"), str(Path(tmp) / name),
                                "--genomes", "40", "--clades", "4", *flags], check=True, capture_output=True)
            one, two = Path(tmp) / "one", Path(tmp) / "two"
            files = sorted(str(p.relative_to(one)) for p in one.rglob("*") if p.is_file())
            self.assertEqual(files, sorted(str(p.relative_to(two)) for p in two.rglob("*") if p.is_file()))
            for f in files:
                if f not in ("bvbrc/genome_amr.jsonl", "truth.json"):
                    self.assertEqual((one / f).read_bytes(), (two / f).read_bytes(), f)
            amr_one = (one / "bvbrc/genome_amr.jsonl").read_text().splitlines()
            amr_two = (two / "bvbrc/genome_amr.jsonl").read_text().splitlines()
            self.assertEqual(amr_two[:len(amr_one)], amr_one)
            truth = json.loads((two / "truth.json").read_text())
            added = [json.loads(line) for line in amr_two[len(amr_one):]]
            self.assertEqual([r["genome_id"] for r in added], [g["genome_id"] for g in truth["genomes"]])
            for r, g in zip(added, truth["genomes"]):
                self.assertEqual(r["antibiotic"], "ciprofloxacin")
                self.assertEqual(r["resistant_phenotype"] == "Resistant", g["resistant_second"])
            plain = json.loads((one / "truth.json").read_text())
            self.assertNotIn("second_drug", plain["config"])
            self.assertNotIn("resistant_second", plain["genomes"][0])
            del truth["config"]["second_drug"]
            for g in truth["genomes"]:
                del g["resistant_second"]
            self.assertEqual(truth, plain)


timeless = base.timeless


def master_manifest(raw, path):
    """A sha256sum-format master manifest of every file under RAW's input parts."""
    lines = []
    for part in run.RAW_PARTS:
        for p in sorted((Path(raw) / part).rglob("*")):
            if p.is_file():
                lines.append(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(raw)}")
    Path(path).write_text("\n".join(lines) + "\n")
    return path


@unittest.skipUnless(os.environ.get("GRM_INPUT_E2E_BINARY"), "set GRM_INPUT_E2E_BINARY to a grm release binary")
class CombineEndToEndTests(unittest.TestCase):
    """G2 on a synthetic root with two datasets (meropenem and, with --second-drug, ciprofloxacin): the comparison
    run whole, and run as two parts of one dataset each, through the real grm; combine of the parts must give the
    whole run's J records and K report, after timings."""

    DATASETS = (DS, "salmonella_enterica__ciprofloxacin")

    @classmethod
    def setUpClass(cls):
        out = os.environ.get("GRM_INPUT_E2E_OUT")
        cls.tmp = None if out else tempfile.TemporaryDirectory()
        cls.out = Path(out or cls.tmp.name) / "combine"
        shutil.rmtree(cls.out, ignore_errors=True)
        cls.out.mkdir(parents=True)
        cls.raw = cls.out / "raw"
        subprocess.run([sys.executable, "-B", str(HERE / "synthetic_root.py"), str(cls.raw), "--genomes", "1200",
                        "--second-drug"], check=True)
        spec = json.loads((HERE / "spec.json").read_text())
        spec.pop("revision12_gate")
        spec["synthetic"] = True
        spec["pool_budget"]["default_kset_bytes"] = 2 ** 20
        # Block C's salmonella_qrdr assay calls on both drugs' genomes; H has no ciprofloxacin panel, so k-h's
        # ciprofloxacin contrasts would not be produced.
        spec["arms"] = [{"name": "cap-600", "cap": 600}, {"name": "k-c", "blocks": ["K", "C"]}]
        spec["datasets"] = [{"id": d, "partition": "confirmatory"} for d in cls.DATASETS]
        cls.spec = cls.out / "spec.json"
        cls.spec.write_text(json.dumps(spec, indent=1))
        cls.parts = []
        for i, d in enumerate(cls.DATASETS, 1):
            path = cls.out / f"part{i}.json"
            path.write_text(json.dumps({**spec, "datasets": [{"id": d, "partition": "confirmatory"}]}, indent=1))
            cls.parts.append(path)
        cls.whole = cls.run_spec(cls.spec, "whole")
        cls.works = [cls.run_spec(p, p.stem) for p in cls.parts]
        cls.master = master_manifest(cls.raw, cls.out / "MASTER.sha256")
        cls.combined = cls.out / "combined"
        cls.combine(cls.combined, *cls.works)

    @classmethod
    def run_spec(cls, spec, name):
        work = cls.out / name
        subprocess.run([sys.executable, "-B", str(HERE / "run.py"), str(spec), "--raw", str(cls.raw), "--work",
                        str(work), "--grm", os.environ["GRM_INPUT_E2E_BINARY"], "--source-commit", "0" * 40,
                        "--jobs", os.environ.get("GRM_INPUT_E2E_JOBS", "3"), "--below-panel-floor",
                        "--sigdrug-vocabulary", str(cls.raw / "vocabulary")], check=True)
        return work

    @classmethod
    def combine(cls, out, *works):
        argv = [sys.executable, "-B", str(HERE / "run.py"), "combine", str(cls.spec), "--master-manifest",
                str(cls.master), "--out", str(out)]
        for w in works:
            argv += ["--work-from", str(w)]
        subprocess.run(argv, check=True)

    @classmethod
    def tearDownClass(cls):
        if cls.tmp:
            cls.tmp.cleanup()

    def test_both_datasets_are_decided_and_each_part_ran_one(self):
        report = json.loads((self.combined / "report" / "ablation-report.json").read_text())
        self.assertEqual(report["contrast_counts"]["planned"], 2 * 3 * 5 * 2)  # 2 datasets, 2 arms and the sentinel
        self.assertEqual(report["contrast_counts"]["not_produced"], 0, report["not_produced"])
        for work, d in zip(self.works, self.DATASETS):
            self.assertEqual(sorted(p.name for p in (work / "arms" / "default").iterdir()), [d])
        groups = [json.loads((w / "species-groups.json").read_text()) for w in (self.whole, *self.works)]
        self.assertEqual(len({g["sha256"] for g in groups}), 1)

    def test_two_parts_combined_equal_the_whole_run(self):
        whole, combined = timeless(self.whole), timeless(self.combined)
        for name, data in whole["j"].items():
            self.assertEqual(combined["j"][name], data, name)
        self.assertEqual(combined["manifest"], whole["manifest"])
        self.assertEqual(combined["report"], whole["report"])
        self.assertEqual(combined["md"], whole["md"])
        for path in sorted((self.whole / "arms").rglob("predictions.jsonl")):
            self.assertEqual((self.combined / path.relative_to(self.whole)).read_bytes(), path.read_bytes(), path)
        record = json.loads((self.combined / "combine.json").read_text())
        self.assertEqual([p["datasets"] for p in record["parts"]], [[d] for d in self.DATASETS])
        self.assertEqual(record["spec_sha256"], whole["manifest"]["protocol_sha256"])
        self.assertEqual(json.loads((self.combined / "summary.json").read_text())["arms"],
                         json.loads((self.whole / "summary.json").read_text())["arms"])

    def test_a_single_work_combined_equals_itself(self):
        """H1-PLAN §7 acceptance 11: combine of a species' single WORK reproduces its own J and K."""
        out = self.out / "combined-whole"
        shutil.rmtree(out, ignore_errors=True)
        self.combine(out, self.whole)
        self.assertEqual(timeless(out), timeless(self.whole))
        for name in ("contrasts", "planned", "missing"):
            self.assertEqual((out / "j" / f"{name}.json").read_bytes(), (self.whole / "j" / f"{name}.json").read_bytes())


if __name__ == "__main__":
    unittest.main()
