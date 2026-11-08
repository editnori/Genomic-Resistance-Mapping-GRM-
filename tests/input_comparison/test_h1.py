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
    """    if argv[0] == "count" and os.environ.get("STUB_INPLACE"):
        # As species.rs's reuse branch rewrites a sidecar whose catalog version differs: in place, same inode.
        for side in list((root / "species" / "salmonella_enterica" / "ksets").glob("*.json")):
            with open(side, "r+") as f:
                f.write(json.dumps({"rewritten_in_place_by": root.name}))
                f.truncate()
    if argv[0] == "count" and os.environ.get("STUB_RECOUNT"):
        for kset in list((root / "species" / "salmonella_enterica" / "ksets").glob("*.kset")):
            (kset.parent / "new.part").write_bytes(kset.read_bytes())
            os.replace(kset.parent / "new.part", kset)
    sys.exit(3 if os.environ.get("STUB_COUNT_FAILS") else 0)""").replace(
    """"reported": "R" if i % 3 else "S"})""",
    """"reported": ("S" if os.environ.get("STUB_FLIP") and g == "g001" else "R") if i % 3 else
    ("R" if os.environ.get("STUB_UNSELECTED") and g == "g048" else "S")})""").replace(
    """"antibiotic": "meropenem"}]))""",
    """"antibiotic": "meropenem", **({"n_resistant": 7} if os.environ.get("STUB_COUNTS") else {})}]))""").replace(
    """(root / "species" / "salmonella_enterica" / "genomes.txt").write_text("\\n".join(ids[:cap]) + "\\n")""",
    """drawn = ids[:cap][::-1] if os.environ.get("STUB_ORDER") else ids[:cap]
if os.environ.get("STUB_SWAP"):
    drawn = drawn[1:] + [ids[-1]]
(root / "species" / "salmonella_enterica" / "genomes.txt").write_text("\\n".join(drawn) + "\\n")""").replace(
    """(root / "species" / "salmonella_enterica" / "genomes.txt").write_text("\\n".join(drawn) + "\\n")""",
    """(root / "species" / "salmonella_enterica" / "genomes.txt").write_text("\\n".join(drawn) + "\\n")
datasets = json.loads((root / "curated" / "datasets.json").read_text())
if os.environ.get("STUB_EXTRA"):
    datasets.append({"id": "salmonella_enterica__extra", "species": "Salmonella enterica",
                     "species_slug": "salmonella_enterica", "antibiotic": "extra"})
if os.environ.get("STUB_OTHER"):
    # A second panel species: 12 eligible genomes o000..o011, the first 8 selected (o001..o008 with
    # STUB_OTHER_SWAP), one dataset; STUB_OTHER_FLIP changes selected genome o002's label.
    other = [f"o{i:03d}" for i in range(12)]
    (root / "species" / "other_sp").mkdir(parents=True, exist_ok=True)
    picked = other[1:9] if os.environ.get("STUB_OTHER_SWAP") else other[:8]
    (root / "species" / "other_sp" / "genomes.txt").write_text("\\n".join(picked) + "\\n")
    datasets.append({"id": "other_sp__drug", "species": "Other sp", "species_slug": "other_sp", "antibiotic": "drug"})
    with open(root / "curated" / "phenotypes.jsonl", "a") as f:
        for i, g in enumerate(other):
            call = "R" if i % 2 else "S"
            if os.environ.get("STUB_OTHER_FLIP") and g == "o002":
                call = "R"
            f.write(json.dumps({"genome_id": g, "species": "Other sp", "antibiotic": "drug", "reported": call}) + "\\n")
(root / "curated" / "datasets.json").write_text(json.dumps(datasets))""")
# What the substance gate tests switch on (main's ruling after H1DATA); each must reach the stub.
GATE_HOOKS = ("STUB_FLIP", "STUB_UNSELECTED", "STUB_COUNTS", "STUB_ORDER", "STUB_SWAP", "STUB_EXTRA", "STUB_OTHER",
              "STUB_OTHER_SWAP", "STUB_OTHER_FLIP")


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
            # The kset is the pool's file; the sidecar a copy of the pool's (REVIEW-H1CODE L1).
            self.assertTrue(os.path.samefile(mine / f"{g}.kset", pool / f"{g}.kset"), g)
            self.assertFalse(os.path.samefile(mine / f"{g}.json", pool / f"{g}.json"), g)
            self.assertEqual((mine / f"{g}.json").read_bytes(), (pool / f"{g}.json").read_bytes(), g)
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
        self.assertIn("grm count rewrote the kset or sidecar of 20 of the 20 genomes shared with the species pool",
                      str(e.exception))
        self.assertFalse((r.curated(r.default) / "built.json").exists())

    def test_a_sidecar_rewritten_in_place_is_refused_and_never_reaches_the_pool(self):
        """REVIEW-H1CODE L1 (and R4): species.rs rewrites a reused sidecar in place when its catalog version differs.
        The root holds a copy, so the write stays in the root, and the changed copy refuses the root."""
        self.estimate()
        r = self.runner()
        pool = self.fake_pool(r, [f"g{i:03d}" for i in range(30)])
        before = {p.name: p.read_bytes() for p in pool.iterdir()}
        os.environ["STUB_INPLACE"] = "1"
        try:
            with self.assertRaises(SystemExit) as e:
                r.root(r.default)
        finally:
            del os.environ["STUB_INPLACE"]
        self.assertIn("grm count rewrote the kset or sidecar of 20 of the 20 genomes shared", str(e.exception))
        self.assertEqual({p.name: p.read_bytes() for p in pool.iterdir()}, before)
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


class GateCase(Stub):
    """The revision-12 gate on the stub: revision 12 is the stub's plain curate (with BASE_HOOKS set), pinned for
    SLUGS; `gated` reruns --estimate-only on a fresh WORK with one hook more."""
    BASE_HOOKS = ()
    SLUGS = (SPEC["species"],)

    def setUp(self):
        super().setUp()
        for hook in GATE_HOOKS:
            self.assertEqual(STUB_GRM.count(f'"{hook}"'), 1, hook)
        self.addCleanup(lambda: [os.environ.pop(h, None) for h in GATE_HOOKS])
        for hook in self.BASE_HOOKS:
            os.environ[hook] = "1"
        self.estimate()
        self.pins = base.gate_of(base.default_root(self.tmp / "work"), list(self.SLUGS))

    def gated(self, hook=None):
        """--estimate-only on a fresh WORK with revision 12's pins, curate changed by `hook`; (exit code, stderr,
        stdout)."""
        shutil.rmtree(self.tmp / "work")
        if hook:
            os.environ[hook] = "1"
        out, err = io.StringIO(), io.StringIO()
        try:
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err), self.assertRaises(SystemExit) as e:
                run.main([str(self.spec_path(revision12_gate=self.pins)), "--raw", str(self.tmp / "raw"), "--work",
                          str(self.tmp / "work"), "--grm", str(self.grm), "--source-commit", "0" * 40,
                          "--below-panel-floor", "--estimate-only", "--free-disk-bytes", str(10 ** 15)])
        finally:
            os.environ.pop(hook or "", None)
        return e.exception.code, err.getvalue(), out.getvalue()


class MultiSpeciesGateTests(GateCase):
    """REVIEW-H1CODE M1: in a run for the stub's Salmonella, the gate checks the other pinned panel species too: a
    changed genome set or a changed selected label of other_sp refuses."""
    BASE_HOOKS = ("STUB_OTHER",)
    SLUGS = (SPEC["species"], "other_sp")

    def test_both_species_as_pinned_pass(self):
        self.assertEqual(sorted(self.pins["genome_set"]), sorted(self.SLUGS))
        self.assertEqual(sorted(self.pins["datasets"]["other_sp"]), ["other_sp__drug"])
        code, err, out = self.gated()
        self.assertEqual(code, 0, err)
        self.assertIn("the default's genome sets for other_sp, salmonella_enterica", err)

    def test_another_species_genome_set_changed_refuses(self):
        code, _, _ = self.gated("STUB_OTHER_SWAP")
        self.assertIn("the default arm's genome set for other_sp is", str(code))

    def test_another_species_selected_label_changed_refuses(self):
        code, _, _ = self.gated("STUB_OTHER_FLIP")
        self.assertIn("the default arm's dataset other_sp__drug is", str(code))
        self.assertIn("a selected genome's phenotype row differs", str(code))


class SubstanceGateTests(GateCase):
    """Main's ruling after H1DATA: the revision-12 gate passes a default with revision 12's genome set and, per
    dataset, its spec outside the eligible-pool counts and its selected genomes' phenotype rows, whatever order curate
    drew; a swapped genome, a changed selected label or an extra dataset refuses (the last by main's ruling after
    REVIEW-H1DATA)."""

    def test_an_extra_dataset_refuses(self):
        code, _, _ = self.gated("STUB_EXTRA")
        self.assertIn("the default arm has 1 dataset(s) of salmonella_enterica that revision 12 hasn't: "
                      "salmonella_enterica__extra", str(code))

    def test_a_duplicate_phenotype_row_refuses(self):
        """REVIEW-H1CODE L5 (R5): two rows for one genome and drug can't be told apart, so substance refuses."""
        root = base.default_root(self.tmp / "work")
        rows = root / "curated" / "phenotypes.jsonl"
        rows.write_text(rows.read_text() + rows.read_text().splitlines()[0] + "\n")
        with self.assertRaises(SystemExit) as e:
            run.substance(root, [SPEC["species"]])
        self.assertIn("two phenotype rows for genome g000 and meropenem", str(e.exception))

    def test_revision_12_itself_passes_as_recorded(self):
        code, err, out = self.gated()
        self.assertEqual(code, 0)
        self.assertIn("revision-12 gate: PASSED", err)
        self.assertTrue(json.loads(out)["revision12_gate"]["recorded"]["as_the_spec_records"])

    def test_the_same_set_in_another_order_passes_and_the_order_is_recorded(self):
        code, err, out = self.gated("STUB_ORDER")
        self.assertEqual(code, 0, err)
        record = json.loads(out)["revision12_gate"]
        self.assertEqual(record["genome_set"], self.pins["genome_set"])
        self.assertNotEqual(record["recorded"]["genome_order"], self.pins["recorded"]["genome_order"])
        self.assertFalse(record["recorded"]["as_the_spec_records"])
        self.assertIn("NOT as the spec records them", err)

    def test_substance_compares_a_species_as_a_sorted_set(self):
        """Another panel species' genomes.txt keeps its drawn order in the default root (only the run's own species is
        rewritten sorted), and revision 12's pins come from drawn-order files too, so the set is sorted first; a
        duplicate id still changes it."""
        def root(name, genomes):
            r = self.tmp / name
            (r / "curated").mkdir(parents=True)
            (r / "species" / "other_sp").mkdir(parents=True)
            (r / "species" / "other_sp" / "genomes.txt").write_text("".join(g + "\n" for g in genomes))
            (r / "curated" / "datasets.json").write_text(json.dumps([{"id": "other_sp__drug", "species": "Other sp",
                                                                      "species_slug": "other_sp", "antibiotic": "drug"}]))
            (r / "curated" / "phenotypes.jsonl").write_text("".join(json.dumps(
                {"genome_id": g, "species": "Other sp", "antibiotic": "drug", "reported": "R"}) + "\n" for g in set(genomes)))
            return run.substance(r, ["other_sp"])

        drawn, other, doubled = root("a", ["x2", "x1", "x3"]), root("b", ["x1", "x3", "x2"]), root("c", ["x1", "x2", "x2", "x3"])
        self.assertEqual(drawn, other)
        self.assertEqual(drawn["genome_set"]["other_sp"], run.sha256_bytes(b"x1\nx2\nx3\n"))
        self.assertNotEqual(doubled["genome_set"], drawn["genome_set"])

    def test_one_genome_swapped_refuses(self):
        code, _, _ = self.gated("STUB_SWAP")
        self.assertIn("the default arm's genome set for salmonella_enterica is", str(code))
        self.assertFalse((self.tmp / "work" / "revision12-gate.json").exists())

    def test_one_selected_label_changed_refuses(self):
        code, _, _ = self.gated("STUB_FLIP")
        self.assertIn("the default arm's dataset salmonella_enterica__meropenem is", str(code))
        self.assertIn("a selected genome's phenotype row differs", str(code))

    def test_eligible_counts_and_an_unselected_genome_s_row_decide_nothing(self):
        code, err, _ = self.gated("STUB_COUNTS")
        self.assertEqual(code, 0, err)
        root = base.default_root(self.tmp / "work")
        self.assertNotEqual(run.datasets_view(root / "curated" / "datasets.json"),
                            self.pins["recorded"]["datasets_json_sha256"])
        code, err, _ = self.gated("STUB_UNSELECTED")
        self.assertEqual(code, 0, err)
        root = base.default_root(self.tmp / "work")
        rows = [json.loads(line) for line in (root / "curated" / "phenotypes.jsonl").read_text().splitlines()]
        self.assertEqual([r["reported"] for r in rows if r["genome_id"] == "g048"], ["R"])
        self.assertNotIn("g048", root.joinpath("species", "salmonella_enterica", "genomes.txt").read_text().split())


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
        identical = "--run-identical-cohorts" in flags
        r.arm(r.default, DS, identical)
        r.arm(r.arms["cap-10"], DS, identical)
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
        # REVIEW-H1CODE L5 (R3): --run-identical-cohorts changes what an identical cohort writes, so it reruns too.
        self.assertEqual(self.run_arms("--run-identical-cohorts")[0], ["default", "cap-10"], "identical cohorts")
        self.assertTrue(json.loads((self.out("cap-10") / run.RESUME).read_text())["run_identical_cohorts"])
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

    def test_a_linked_arm_folder_or_arm_parent_refuses_the_run(self):
        """REVIEW-H1CODE-B L-a: a dataset folder, or its arm folder, that is a link to another arm's folder inside
        arms/ is refused before anything is resumed or written through it."""
        self.run_arms()
        default = self.out("default")
        before = {p.name: p.read_bytes() for p in default.iterdir() if p.is_file()}
        shutil.rmtree(self.out("cap-10"))
        self.out("cap-10").symlink_to(default, target_is_directory=True)
        with self.assertRaises(SystemExit) as e:
            self.run_arms()
        self.assertIn(f"{self.out('cap-10')} is a symlink", str(e.exception))
        self.assertEqual({p.name: p.read_bytes() for p in default.iterdir() if p.is_file()}, before)
        self.out("cap-10").unlink()
        arm = self.out("cap-10").parent
        shutil.rmtree(arm)
        (arm.parent / "elsewhere").mkdir()
        arm.symlink_to(arm.parent / "elsewhere", target_is_directory=True)
        with self.assertRaises(SystemExit) as e:
            self.run_arms()
        self.assertIn(f"{arm} is a symlink", str(e.exception))
        self.assertEqual(list((arm.parent / "elsewhere").iterdir()), [])

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


SENTINEL_NAMES = [s["name"] for s in run.SENTINELS]


def arm_folder(out, name, ds, record, plans=None):
    """A finished arm folder as a run leaves it: runs.json, arm.json, predictions.jsonl, `plans` (name -> text) under
    plans/, and resume.json with `record` and every file's sha256."""
    out.mkdir(parents=True)
    for f, text in (("runs.json", "[]"), ("arm.json", "{}"), ("predictions.jsonl", f'{{"arm": "{name}"}}\n')):
        (out / f).write_text(text)
    for p, text in (plans or {}).items():
        (out / "plans").mkdir(exist_ok=True)
        (out / "plans" / p).write_text(text)
    outputs = {str(p.relative_to(out)): run.file_sha256(p) for p in sorted(out.rglob("*")) if p.is_file()}
    (out / run.RESUME).write_text(run.canonical({"arm": name, "dataset": ds, "trace_training": False,
                                                 "run_identical_cohorts": False, **record, "outputs": outputs}))


class ReferenceTests(Stub):
    """G2(b): a k add-on (--reference-from) takes a finished run's default arm, sentinels and grouping, checked and
    copied, and runs only its k arms, with no species pool; a reference run that differs in what the add-on would
    rest on refuses before anything is counted."""

    def setUp(self):
        super().setUp()
        self.full = json.loads(self.spec_path().read_text())   # default cap 20; cap-10, cap-30 and k21
        self.addon_spec = self.tmp / "addon.json"
        self.addon_spec.write_text(json.dumps({**self.full, "arms": [{"name": "k21", "k": 21}]}))
        self.ref = self.reference(self.tmp / "main")

    def reference(self, work, arms=None, raw="", inputs=None, manifest=None, summary=None, trace=False, spec=None):
        """A finished run of the full spec without k21 (or with `arms`, and `spec`'s changes) on DS, as `combine`
        takes a part: its records and the default's, cap arms' and sentinels' outputs, the default with one plan."""
        spec = {**self.full, "arms": arms if arms is not None else [a for a in self.full["arms"] if a["name"] != "k21"],
                **(spec or {})}
        (work / "k").mkdir(parents=True)
        ins = {"raw_manifest_sha256": run.sha256_bytes(raw.encode()), "raw_files": len(raw.split("\n")) if raw else 0,
               "grm_sha256": run.file_sha256(self.grm), "grm": str(self.grm), "vocabulary_sha256": None,
               **(inputs or {})}
        man = {"protocol_sha256": run.sha256_bytes(run.canonical(spec).encode()), "source_commit": "0" * 40,
               "binary_sha256": ins["grm_sha256"], "synthetic": True, "revision12_gate": "not_required",
               "compute": {"actual_seconds": 1.0}, **(manifest or {})}
        self.groups_record = {"structure": f"{work}/roots/pool/structure.json", "sha256": "f" * 64, "genomes": 30,
                              "lineages": 3, "arms": ["default", *(a["name"] for a in spec["arms"]), *SENTINEL_NAMES]}
        for path, value in ((run.SPEC_COPY, run.canonical(spec)), (run.INPUTS, run.canonical(ins)),
                            ("k/manifest.json", json.dumps(man)), ("species-groups.json", run.canonical(self.groups_record)),
                            ("summary.json", json.dumps({"eligible_genomes": 50, "panel_floor": 2001, "seconds": 1.0,
                                                         **(summary or {})})),
                            (run.RAW_MANIFEST, raw)):
            (work / path).write_text(value)
        record = {"inputs_sha256": run.file_sha256(work / run.INPUTS), "spec_sha256": man["protocol_sha256"],
                  "trace_training": trace}
        arm_folder(work / "arms" / "default" / DS, "default", DS, record, {"lineage-0.json": '{"tag": "main"}'})
        for name in [a["name"] for a in spec["arms"]] + SENTINEL_NAMES:
            arm_folder(work / "arms" / name / DS, name, DS, record)
        return work

    def addon(self, *flags, spec=None, ref=None):
        """run.main of the add-on with the stub grm and J and K left out; returns the refusal, or None."""
        def report_run(*args):
            return {"arms": {}, "no_decision_datasets": [], "refused_plans": []}
        argv = [str(spec or self.addon_spec), "--raw", str(self.tmp / "raw"), "--work", str(self.tmp / "work"),
                "--grm", str(self.grm), "--source-commit", "0" * 40, "--below-panel-floor", "--free-disk-bytes",
                str(10 ** 15), "--reference-from", str(ref or self.ref), *flags]
        with mock.patch.object(run, "report_run", report_run), contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()):
            try:
                run.main(argv)
            except SystemExit as e:
                return str(e.code)
        return None

    def work(self, *parts):
        return self.tmp.joinpath("work", *parts)

    def calls(self):
        p = self.work("arm-calls.txt")
        return p.read_text().split() if p.exists() else []

    def test_the_add_on_takes_the_default_and_sentinels_and_runs_only_its_k_arm(self):
        self.assertIsNone(self.addon())
        self.assertEqual(self.calls(), ["k21"])
        for name in ["default", *SENTINEL_NAMES]:
            self.assertTrue(run.same_folder(self.work("arms", name, DS), self.ref / "arms" / name / DS), name)
        # No species pool and no k31 root: only the k21 root is counted.
        counted = self.work("counted.txt").read_text().splitlines()
        self.assertEqual([c.split()[0] for c in counted], ["count", "structure", "matrix"])
        self.assertIn("--k 21", counted[0])
        self.assertEqual([p.name for p in self.work("roots").glob("pool-*")], [])
        disk = json.loads(self.work("disk-estimate.json").read_text())
        self.assertEqual((disk["pool_bytes"], [r["k"] for r in disk["roots"].values()]), (0, [21]))
        ref = json.loads(self.work(run.REFERENCE).read_text())
        self.assertEqual({k: ref[k] for k in ("work", "spec_sha256", "inputs_sha256", "species_groups")},
                         {"work": str(self.ref.resolve()), "spec_sha256": run.sha256_bytes(
                             (self.ref / run.SPEC_COPY).read_bytes()),
                          "inputs_sha256": run.file_sha256(self.ref / run.INPUTS), "species_groups": self.groups_record})
        self.assertEqual(sorted(ref["taken"]), sorted(f"{n}/{DS}" for n in ["default", *SENTINEL_NAMES]))
        self.assertEqual(json.loads(self.work("species-groups.json").read_text()), self.groups_record)
        k21 = json.loads(self.work("arms", "k21", DS, run.RESUME).read_text())
        self.assertEqual((k21["groups_sha256"], k21["reference_plans"]),
                         ("f" * 64, {"lineage-0.json": run.file_sha256(self.ref / "arms/default" / DS / "plans/lineage-0.json")}))
        self.assertEqual(k21["spec_sha256"], run.sha256_bytes(run.canonical(json.loads(self.addon_spec.read_text())).encode()))
        self.assertEqual(json.loads(self.work("summary.json").read_text())["reference_from"], ref)

    def test_a_rerun_resumes_the_k_arm_on_a_fresh_copy(self):
        self.assertIsNone(self.addon())
        (self.work("arms", "default", DS, "predictions.jsonl")).write_text("damaged copy\n")
        self.assertIsNone(self.addon())
        self.assertEqual(self.calls(), ["k21"])
        self.assertTrue(run.same_folder(self.work("arms", "default", DS), self.ref / "arms/default" / DS))

    def test_estimate_only_plans_the_k_roots_alone(self):
        disk = self.estimate("--reference-from", str(self.tmp / "not-read"), spec=self.addon_spec)
        self.assertEqual((disk["pool_bytes"], disk["genomes"], [(r["k"], r["genomes"]) for r in disk["roots"].values()]),
                         (0, 20, [(21, 20)]))

    def test_a_reference_that_differs_refuses_before_anything_is_counted(self):
        full = self.full
        k21 = {"name": "k21", "k": 21}
        other_partition = [{"id": DS, "partition": "exploratory"}]
        cases = (
            ("a cohort arm", {"spec": {**full, "arms": [{"name": "cap-10", "cap": 10}]}},
             "--reference-from runs k arms only (the k add-on, G2(b)); this spec's arms change cap-10: cap"),
            ("no arm", {"spec": {**full, "arms": []}}, "this spec's arms change nothing: it has none"),
            ("other margin", {"spec": {**full, "arms": [k21], "margin": {"log_loss": 0.03}}},
             "its spec differs from {spec} in margin, not only in its arms and datasets"),
            ("other partition", {"ref": {"spec": {"datasets": other_partition}}},
             f"did not run {DS} with this spec's partition"),
            ("ran k21", {"ref": {"arms": full["arms"]}}, "already ran k21"),
            ("grm", {"ref": {"inputs": {"grm_sha256": "0" * 64}}}, "ran another grm binary"),
            ("vocabulary", {"ref": {"inputs": {"vocabulary_sha256": "d" * 64}}}, "ran another vocabulary"),
            ("binary", {"ref": {"manifest": {"binary_sha256": "0" * 64}}},
             "K's binary_sha256 is not the grm its inputs.json hashed"),
            ("commit", {"ref": {"manifest": {"source_commit": "1" * 40}}}, f"ran source commit {'1' * 40}, not {'0' * 40}"),
            ("synthetic", {"ref": {"manifest": {"synthetic": False}}}, "K's manifest records synthetic False"),
            ("gate", {"ref": {"manifest": {"revision12_gate": "passed"}}}, "and gate passed"),
            ("eligible", {"ref": {"summary": {"eligible_genomes": 49}}}, "counted 49 eligible genomes"),
            ("floor", {"ref": {"summary": {"panel_floor": 7}}}, "against a floor of 7"),
            ("protocol", {"ref": {"manifest": {"protocol_sha256": "0" * 64}}},
             "K's manifest was not made from its spec.json"),
            ("raw", {"raw": "bvbrc/genomes.jsonl"},
             "1 of this run's 1 RAW files are not in its raw-manifest.txt with that size and sha256, "
             "e.g. bvbrc/genomes.jsonl"),
            ("trace", {"ref": {"trace": True}}, "ran with --trace-training True and --run-identical-cohorts False"),
            ("add-on", {"mark": run.REFERENCE}, "is itself a k add-on"),
            ("add-on by its summary", {"ref": {"summary": {"reference_from": {"work": "/elsewhere"}}}},
             "is itself a k add-on"),
            ("damaged", {"damage": "default"}, f"arm default on {DS}: predictions.jsonl changed after the arm finished"),
            ("missing sentinel", {"drop": "sentinel-seed-twin"}, f"has no output of arm sentinel-seed-twin on {DS}"),
            ("unfinished", {"drop_file": "summary.json"}, "no readable summary.json"),
        )
        for label, case, message in cases:
            with self.subTest(label):
                shutil.rmtree(self.tmp / "work", ignore_errors=True)
                shutil.rmtree(self.tmp / "raw", ignore_errors=True)
                ref = self.reference(self.tmp / f"ref-{label.replace(' ', '-')}", **case.get("ref", {}))
                spec = None
                if "spec" in case:
                    spec = self.tmp / f"addon-{label.replace(' ', '-')}.json"
                    spec.write_text(json.dumps(case["spec"]))
                if "raw" in case:
                    (self.tmp / "raw" / "bvbrc").mkdir(parents=True)
                    (self.tmp / "raw" / case["raw"]).write_text("{}\n")
                if "mark" in case:
                    (ref / case["mark"]).write_text("{}")
                if "damage" in case:
                    (ref / "arms" / case["damage"] / DS / "predictions.jsonl").write_text("")
                if "drop" in case:
                    shutil.rmtree(ref / "arms" / case["drop"])
                if "drop_file" in case:
                    (ref / case["drop_file"]).unlink()
                code = self.addon(spec=spec, ref=ref)
                self.assertIsNotNone(code)
                self.assertIn(message.format(spec=spec), code)
                self.assertFalse(self.work("counted.txt").exists(), code)
                self.assertEqual(self.calls(), [], code)

    def test_a_k_arm_whose_cohort_is_not_the_default_s_refuses(self):
        """A k arm's cohort is the default's (load_spec lets a k arm change nothing else); one that held another would
        change the union the reference run's grouping was made from."""
        real = run.Runner.genome_set

        def genome_set(runner, arm):
            s = real(runner, arm)
            return s - {sorted(s)[0]} if arm["name"] == "k21" else s
        with mock.patch.object(run.Runner, "genome_set", genome_set):
            code = self.addon()
        self.assertIn("arm k21's cohort is not the default's, so it would change the union the grouping of "
                      f"--reference-from {self.ref} was made from", code)
        self.assertFalse(self.work("counted.txt").exists(), code)
        self.assertEqual(self.calls(), [])

    def test_a_copy_that_differs_from_the_reference_refuses(self):
        """Each copied folder is checked against its source after the copy, byte for byte."""
        real = shutil.copytree

        def copytree(src, dst, *args, **kw):   # shutil's own recursion calls it too
            out = real(src, dst, *args, **kw)
            if Path(src).parent.name == "default":
                (Path(dst) / "predictions.jsonl").write_text("changed in the copy\n")
            return out
        with mock.patch.object(run.shutil, "copytree", copytree):
            code = self.addon()
        self.assertIn(f"the copy of arm default on {DS} from {self.ref} differs from it, or one of them holds a link",
                      code)
        self.assertEqual(self.calls(), [])

    def test_a_reference_arm_folder_holding_a_link_refuses(self):
        """A link in the reference run's arm folder refuses before the copy, as `combine` refuses a part's."""
        (self.tmp / "target.txt").write_text("x")
        (self.ref / "arms" / "default" / DS / "extra").symlink_to(self.tmp / "target.txt")
        code = self.addon()
        self.assertIn(f"part {self.ref}: arm default on {DS}: its folder holds a symlink", code)
        self.assertEqual(self.calls(), [])

    def test_a_linked_copy_destination_refuses(self):
        """The copy never writes through a link in the add-on's WORK."""
        (self.tmp / "elsewhere").mkdir()
        self.work("arms").mkdir(parents=True)
        self.work("arms", "default").symlink_to(self.tmp / "elsewhere", target_is_directory=True)
        code = self.addon()
        self.assertIn("is a symlink", code)
        self.assertEqual(list((self.tmp / "elsewhere").iterdir()), [])


class CombineRefusalTests(unittest.TestCase):
    """G2: `combine` refuses parts that are not one spec split by datasets on one pool, grm, vocabulary, source
    commit and master manifest. A refusal comes before any arm output is read; valid parts get as far as that."""

    A, B = DS, "salmonella_enterica__ciprofloxacin"

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmp = Path(tmp.name)
        self.master = self.tmp / "MASTER.sha256"
        self.master.write_text(f"{'a' * 64}  bvbrc/genomes.jsonl\n{'b' * 64}  ./fasta/g1.fna.gz\n{'0' * 64}  other\n")
        self.spec = {**copy.deepcopy(SPEC), "datasets": [{"id": self.A, "partition": "confirmatory"},
                                                         {"id": self.B, "partition": "confirmatory"}]}
        self.pin()

    def pin(self):
        """Pins the master manifest file as it is now in the spec (REVIEW-H1CODE L6) and writes spec.json."""
        self.spec["data_root"] = {"path": str(self.tmp), "manifest": self.master.name,
                                  "manifest_sha256": run.file_sha256(self.master)}
        (self.tmp / "spec.json").write_text(json.dumps(self.spec))

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
                             ({"species": "other"}, "species")):
            b = self.part(f"b-{key}", [self.B], spec_changes=changes)
            self.assertEqual(self.combine(a, b), f"part {b}: its spec differs from {self.tmp / 'spec.json'} in {key}, "
                                                 "not only in its datasets")
        # G2(b): a part may leave out k arms only, and keeps SPEC's order and statement of the rest.
        for label, arms in (("budget-16000 left out", SPEC["arms"][:-1]), ("reordered", SPEC["arms"][::-1]),
                            ("restated", [*SPEC["arms"][:-1], {"name": "budget-16000", "screened_features": 8000}]),
                            ("extra", [*SPEC["arms"], {"name": "cap-500", "cap": 500}])):
            b = self.part(f"b-{label.replace(' ', '-')}", [self.B], spec_changes={"arms": arms})
            self.assertEqual(self.combine(a, b), f"part {b}: its spec differs from {self.tmp / 'spec.json'} in arms, "
                                                 "not only in its datasets (a part may leave out only k arms, which a "
                                                 "k add-on part runs)", label)
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
        # The pool's path may differ (another WORK), and so may its arm list (a main part without the k arms); its
        # grouping may not.
        b = self.part("b-path", [self.B], groups={"structure": "/elsewhere/structure.json"})
        self.assertIn("has no output of arm default", self.combine(a, b))
        b = self.part("b-arms", [self.B], groups={"arms": ["default", "cap-1000"]})
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
        self.pin()
        self.assertIn("line 1: not a sha256sum line", self.combine(self.part("a2", [self.A]), self.part("c", [self.B])))

    def test_the_master_manifest_must_be_the_one_the_spec_pins(self):
        """REVIEW-H1CODE L6: any other file given as --master-manifest refuses, even one holding every RAW line."""
        a, b = self.part("a", [self.A]), self.part("b", [self.B])
        pinned = run.file_sha256(self.master)
        self.master.write_text(self.master.read_text() + f"{'1' * 64}  extra\n")
        self.assertEqual(self.combine(a, b), f"master manifest {self.master} is {run.file_sha256(self.master)}, not "
                                             f"the {pinned} {self.tmp / 'spec.json'} pins")
        del self.spec["data_root"]
        (self.tmp / "spec.json").write_text(json.dumps(self.spec))
        a, b = self.part("a-nopin", [self.A]), self.part("b-nopin", [self.B])
        self.assertIn("pins no data_root.manifest_sha256", self.combine(a, b))

    def test_a_spec_with_an_unpinned_data_root_never_loads(self):
        """Main's ruling after REVIEW-H1DATA: data_root.manifest_sha256 must be a sha256."""
        for value in ("PLACEHOLDER", None, "7dd00779", "A" * 64):
            self.spec["data_root"]["manifest_sha256"] = value
            (self.tmp / "spec.json").write_text(json.dumps(self.spec))
            with self.assertRaises(SystemExit, msg=value) as e:
                run.load_spec(self.tmp / "spec.json")
            self.assertIn(f"data_root.manifest_sha256 is {value!r}, not a sha256: pin the data root's manifest",
                          str(e.exception))

    def arm_outputs(self, work, datasets, damaged=None, **record):
        """A finished arm folder with its resume.json, for every arm of the spec on `datasets`, as a run writes them;
        `record` overrides the record's fields; `damaged` = (arm, dataset) truncates its predictions afterwards."""
        _, _, arms = run.load_spec(self.tmp / "spec.json")
        manifest = json.loads((Path(work) / "k" / "manifest.json").read_text())
        for name in [run.DEFAULT_ARM, *(a["name"] for a in arms)]:
            for ds in datasets:
                out = Path(work) / "arms" / name / ds
                out.mkdir(parents=True)
                for f, text in (("runs.json", "[]"), ("arm.json", "{}"), ("predictions.jsonl", '{"p": 0.5}\n')):
                    (out / f).write_text(text)
                outputs = {f: run.file_sha256(out / f) for f in ("arm.json", "predictions.jsonl", "runs.json")}
                (out / run.RESUME).write_text(run.canonical({
                    "arm": name, "dataset": ds, "inputs_sha256": run.file_sha256(Path(work) / run.INPUTS),
                    "spec_sha256": manifest["protocol_sha256"], "trace_training": False,
                    "run_identical_cohorts": False, **record, "outputs": outputs}))
                if damaged == (name, ds):
                    (out / "predictions.jsonl").write_text("")

    def reaches_the_report(self, *works):
        def report_run(*args):
            raise SystemExit("REPORT REACHED")
        with mock.patch.object(run, "report_run", report_run):
            return self.combine(*works)

    def test_combine_takes_only_arm_outputs_its_parts_recorded(self):
        """REVIEW-H1CODE L2: each copied arm folder must match its part's resume.json, and the parts' arms must share
        --trace-training and --run-identical-cohorts. Clean parts reach the report."""
        a, b = self.part("a", [self.A]), self.part("b", [self.B])
        self.arm_outputs(a, [self.A])
        self.arm_outputs(b, [self.B])
        self.assertEqual(self.reaches_the_report(a, b), "REPORT REACHED")
        cap = a / "arms" / "cap-1000" / self.A
        refusals = (
            ("predictions.jsonl", lambda: (cap / "predictions.jsonl").write_text(""),
             "predictions.jsonl changed after the arm finished"),
            ("extra file", lambda: (cap / "extra.json").write_text("{}"), "its files differ from resume.json's list: extra.json"),
            ("missing record", lambda: (cap / run.RESUME).unlink(), "no readable resume.json"),
            ("link", lambda: (cap / "plans").symlink_to(self.tmp, target_is_directory=True), "its folder holds a symlink"))
        for label, damage, message in refusals:
            shutil.rmtree(a / "arms")
            self.arm_outputs(a, [self.A])
            damage()
            code = self.reaches_the_report(a, b)
            self.assertIn(f"part {a}: arm cap-1000 on {self.A}: {message}", code, label)
        for label, record, message in (
                ("other arm", {"arm": "k-h"}, "resume.json is the record of k-h/"),
                ("other inputs", {"inputs_sha256": "0" * 64}, "was not made from the part's inputs.json"),
                ("other spec", {"spec_sha256": "0" * 64}, "was not made from the part's spec")):
            shutil.rmtree(a / "arms")
            self.arm_outputs(a, [self.A], **record)
            self.assertIn(message, self.reaches_the_report(a, b), label)
        for flag in ("run_identical_cohorts", "trace_training"):
            shutil.rmtree(a / "arms")
            shutil.rmtree(b / "arms")
            self.arm_outputs(a, [self.A])
            self.arm_outputs(b, [self.B], **{flag: True})
            self.assertIn("the parts' arms ran with different --trace-training / --run-identical-cohorts",
                          self.reaches_the_report(a, b), flag)
        # REVIEW-H1CODE-B L-a: the dataset folder, or its arm folder, a link to a valid copy of the same output.
        for level in ("dataset", "arm"):
            shutil.rmtree(a / "arms")
            self.arm_outputs(a, [self.A])
            folder = cap if level == "dataset" else cap.parent
            shutil.move(folder, self.tmp / f"moved-{level}")
            folder.symlink_to(self.tmp / f"moved-{level}", target_is_directory=True)
            self.assertIn(f"part {a}: arm cap-1000 on {self.A}: {folder} is a symlink", self.reaches_the_report(a, b),
                          level)
        # REVIEW-H1CODE-B L-b: a record without a flag refuses, naming the file.
        for flag in ("trace_training", "run_identical_cohorts"):
            shutil.rmtree(a / "arms")
            self.arm_outputs(a, [self.A])
            record = json.loads((cap / run.RESUME).read_text())
            del record[flag]
            (cap / run.RESUME).write_text(run.canonical(record))
            self.assertIn(f"part {a}: arm cap-1000 on {self.A}: {cap / run.RESUME} records no {flag}",
                          self.reaches_the_report(a, b), flag)
        # The same flags in both parts pass.
        shutil.rmtree(a / "arms")
        shutil.rmtree(b / "arms")
        self.arm_outputs(a, [self.A], run_identical_cohorts=True)
        self.arm_outputs(b, [self.B], run_identical_cohorts=True)
        self.assertEqual(self.reaches_the_report(a, b), "REPORT REACHED")

    # G2(b): main parts that leave out k21 and k25, and k add-on parts that ran them against a main part.
    K_ARMS = ("k21", "k25")

    def main_part(self, name, datasets, **changes):
        """A main part without k21 and k25, whose species-groups.json lists its own arms, as its run writes it."""
        arms = [x for x in SPEC["arms"] if x["name"] not in self.K_ARMS]
        own = {"arms": ["default", *(x["name"] for x in arms), *SENTINEL_NAMES]}
        return self.part(name, datasets, spec_changes={"arms": arms, **changes.pop("spec_changes", {})},
                         groups={**own, **changes.pop("groups", {})}, **changes)

    def addon_part(self, name, datasets, main, k=K_ARMS, arms=None, reference=None, **changes):
        """A k add-on part of `k` (or `arms` as stated) on `datasets`, with the reference.json a run with
        --reference-from `main` writes."""
        stated = arms if arms is not None else [x for x in SPEC["arms"] if x["name"] in k]
        # The add-on's species-groups.json is its main part's (take_reference).
        groups = {**json.loads((main / "species-groups.json").read_text()), **changes.pop("groups", {})}
        work = self.part(name, datasets, spec_changes={"arms": stated}, groups=groups, **changes)
        man = json.loads((main / "k" / "manifest.json").read_text())
        record = {"work": str(main), "spec_sha256": man["protocol_sha256"],
                  "inputs_sha256": run.file_sha256(main / run.INPUTS),
                  "species_groups": json.loads((main / "species-groups.json").read_text()), "taken": {}}
        (work / run.REFERENCE).write_text(reference if isinstance(reference, str) else
                                          run.canonical({**record, **(reference or {})}))
        return work

    def outputs(self, work, names, datasets, plans=None, **record):
        """`names`' finished arm folders on `datasets`, made from WORK's inputs.json and spec (`record` overrides)."""
        man = json.loads((work / "k" / "manifest.json").read_text())
        rec = {"inputs_sha256": run.file_sha256(work / run.INPUTS), "spec_sha256": man["protocol_sha256"], **record}
        for name in names:
            for ds in datasets:
                arm_folder(work / "arms" / name / ds, name, ds, rec, plans if name == run.DEFAULT_ARM else None)

    def finish_main(self, work, datasets):
        names = [x["name"] for x in json.loads((work / run.SPEC_COPY).read_text())["arms"]]
        self.outputs(work, [run.DEFAULT_ARM, *names, *SENTINEL_NAMES], datasets, plans={"lineage-0.json": '{"m": 1}'})

    def finish_addon(self, work, main, datasets, plans=None):
        """The add-on's copies of `main`'s default and sentinels, and its k arms scored on `main`'s default plans
        (or `plans`)."""
        for ds in datasets:
            for name in (run.DEFAULT_ARM, *SENTINEL_NAMES):
                shutil.copytree(main / "arms" / name / ds, work / "arms" / name / ds)
            ref = plans if plans is not None else {p.name: run.file_sha256(p) for p in
                                                   sorted((main / "arms" / run.DEFAULT_ARM / ds / "plans").glob("*.json"))}
            names = [x["name"] for x in json.loads((work / run.SPEC_COPY).read_text())["arms"]]
            self.outputs(work, names, [ds], reference_plans=ref)

    def addon_set(self):
        """Main parts a (A) and b (B) without k21 and k25; add-on ka runs both on A against a, kb21 and kb25 one
        each on B against b; all finished."""
        a, b = self.main_part("a", [self.A]), self.main_part("b", [self.B])
        self.finish_main(a, [self.A])
        self.finish_main(b, [self.B])
        ka = self.addon_part("ka", [self.A], a)
        kb21, kb25 = self.addon_part("kb21", [self.B], b, k=("k21",)), self.addon_part("kb25", [self.B], b, k=("k25",))
        self.finish_addon(ka, a, [self.A])
        self.finish_addon(kb21, b, [self.B])
        self.finish_addon(kb25, b, [self.B])
        return a, b, ka, kb21, kb25

    def test_k_add_on_parts_complete_main_parts(self):
        a, b, ka, kb21, kb25 = self.addon_set()
        self.assertEqual(self.reaches_the_report(a, b, ka, kb21, kb25), "REPORT REACHED")
        self.assertEqual(self.reaches_the_report(kb25, ka, b, kb21, a), "REPORT REACHED")
        spec = self.tmp / "spec.json"
        self.assertEqual(self.combine(a, b, ka), f"no part holds arm k21 on {self.B}, k25 on {self.B} of {spec}")
        self.assertEqual(self.combine(a, b, ka, kb21), f"no part holds arm k25 on {self.B} of {spec}")
        self.assertEqual(self.combine(a, b), f"no part holds arm k21 on {self.A}, k25 on {self.A}, k21 on {self.B}, "
                                             f"k25 on {self.B} of {spec}")
        # An add-on alone holds no dataset: its datasets belong to the main part it ran against.
        self.assertEqual(self.combine(ka, kb21, kb25), f"no part holds {self.A}, {self.B} of {spec}")

    def test_an_arm_on_a_dataset_is_in_one_part(self):
        a, b, ka, kb21, kb25 = self.addon_set()
        ka21 = self.addon_part("ka21", [self.A], a, k=("k21",))
        self.assertEqual(self.combine(a, b, ka, kb21, kb25, ka21), f"arm k21 on {self.A} is in two parts: {ka} and {ka21}")
        whole_a = self.part("whole-a", [self.A], groups={"arms": ["default", *(x["name"] for x in SPEC["arms"]),
                                                                 *SENTINEL_NAMES]})   # it ran k21 and k25 itself
        self.assertEqual(self.combine(whole_a, b, ka, kb21, kb25), f"arm k21 on {self.A} is in two parts: {whole_a} and {ka}")

    def test_an_add_on_holds_only_k_arms_as_spec_states_them(self):
        a, b, ka, kb21, kb25 = self.addon_set()
        spec = self.tmp / "spec.json"
        for label, arms, names in (("cohort arm", [{"name": "cap-1000", "cap": 1000}], "cap-1000"),
                                   ("restated", [{"name": "k21", "k": 23}], "k21"), ("none", [], "none"),
                                   ("twice", [{"name": "k21", "k": 21}] * 2, "k21, k21")):
            x = self.addon_part(f"x-{label}", [self.A], a, arms=arms)
            self.assertEqual(self.combine(a, b, x, kb21, kb25), f"part {x}: a k add-on part (--reference-from) holds "
                                                                f"only k arms of {spec} as it states them, not {names}",
                             label)

    def test_an_add_on_names_the_main_part_that_holds_its_datasets(self):
        a, b, ka, kb21, kb25 = self.addon_set()
        both = self.addon_part("both", [self.A, self.B], a, k=("k21",))
        ka25 = self.addon_part("ka25", [self.A], a, k=("k25",))
        self.assertEqual(self.combine(a, b, both, ka25, kb25),
                         f"part {both}: a k add-on part takes its datasets from one run, not 2")
        for label, ref in (("other run", {"spec_sha256": json.loads((b / "k/manifest.json").read_text())["protocol_sha256"],
                                          "inputs_sha256": run.file_sha256(b / run.INPUTS)}),
                           ("other inputs", {"inputs_sha256": "0" * 64}), ("no record", {"spec_sha256": None})):
            x = self.addon_part(f"x-{label.replace(' ', '-')}", [self.A], a, reference=ref)
            self.assertEqual(self.combine(a, b, x, kb21, kb25), f"part {x}: its reference.json names another run than "
                                                                f"{a}, the part that holds its datasets", label)
        x = self.addon_part("x-list", [self.A], a, reference="[]")
        self.assertEqual(self.combine(a, b, x, kb21, kb25), f"part {x}: reference.json is not a reference record")
        x = self.addon_part("x-pool", [self.A], a, groups={"sha256": "0" * 64})
        self.assertIn("the parts differ in the species pool", self.combine(a, b, x, kb21, kb25))

    def test_an_add_on_ran_against_its_main_part_s_default_arm(self):
        a, b, ka, kb21, kb25 = self.addon_set()
        default = ka / "arms" / run.DEFAULT_ARM / self.A
        for label, damage in (
                ("predictions", lambda: (default / "predictions.jsonl").write_text("")),
                ("plan", lambda: (default / "plans" / "lineage-0.json").write_text('{"m": 2}')),
                ("extra file", lambda: (default / "extra.json").write_text("{}")),
                # the same bytes, through a link to the main part's file
                ("file linked", lambda: ((default / "predictions.jsonl").unlink(), (default / "predictions.jsonl")
                                         .symlink_to(a / "arms" / run.DEFAULT_ARM / self.A / "predictions.jsonl"))),
                ("sentinel gone", lambda: shutil.rmtree(ka / "arms" / "sentinel-tuning-twin" / self.A)),
                ("linked", lambda: (shutil.rmtree(default), default.symlink_to(a / "arms" / run.DEFAULT_ARM / self.A,
                                                                               target_is_directory=True)))):
            shutil.rmtree(ka / "arms")
            self.finish_addon(ka, a, [self.A])
            damage()
            name = "sentinel-tuning-twin" if label == "sentinel gone" else run.DEFAULT_ARM
            self.assertEqual(self.reaches_the_report(a, b, ka, kb21, kb25),
                             f"part {ka}: arm {name} on {self.A} is not the copy of {a}'s it ran against", label)
        shutil.rmtree(ka / "arms")
        self.finish_addon(ka, a, [self.A], plans={"lineage-0.json": "0" * 64})
        self.assertEqual(self.reaches_the_report(a, b, ka, kb21, kb25),
                         f"part {ka}: arm k21 on {self.A} was scored on other reference plans than the default arm "
                         f"of {a}")
        # Its own k arms are checked as any part's.
        shutil.rmtree(ka / "arms")
        self.finish_addon(ka, a, [self.A])
        (ka / "arms" / "k25" / self.A / "predictions.jsonl").write_text("")
        self.assertEqual(self.reaches_the_report(a, b, ka, kb21, kb25),
                         f"part {ka}: arm k25 on {self.A}: predictions.jsonl changed after the arm finished")

    def test_out_must_be_new_or_empty_and_parts_finished(self):
        a, b = self.part("a", [self.A]), self.part("b", [self.B])
        (self.tmp / "full" / "x").mkdir(parents=True)
        self.assertEqual(self.combine(a, b, out=self.tmp / "full"), f"--out {self.tmp / 'full'} is not an empty folder")
        (b / "summary.json").unlink()
        self.assertIn(f"part {b}: no readable summary.json", self.combine(a, b))


class H1SpecTests(unittest.TestCase):
    """G10: scripts/input-comparison/h1/ holds the pilot (D1), the four species' stage-1 specs with all of spec.json's
    arms (Layth lifted register section 13 G2 on 8 Nov 2026), and the parts each runs as: main parts with every arm but
    k21 and k25, and one k add-on part (G2(b)) per main part and k arm. Each is spec.json with only its species,
    datasets, kset size (D2) and arms changed."""

    K_ARMS = ("k21", "k25")
    MAIN_PARTS = {"salmonella_enterica": 2, "mycobacterium_tuberculosis": 1, "staphylococcus_aureus": 1,
                  "neisseria_gonorrhoeae": 1}

    H1 = HERE / "h1"
    COUNTS = {"salmonella_enterica": 25, "mycobacterium_tuberculosis": 22, "staphylococcus_aureus": 18,
              "neisseria_gonorrhoeae": 7}

    def specs(self):
        return {p.name: json.loads(p.read_text()) for p in sorted(self.H1.glob("*.json"))}

    def load(self, name):
        return run.load_spec(self.H1 / name)

    def test_each_spec_names_the_h1_data_root_and_pins_its_manifest(self):
        """Main's ruling after REVIEW-H1DATA: H-1 runs from its own read-only data root, and combine pins its
        manifest, MASTER-MANIFEST-h1.sha256, final at 7dd00779."""
        for name, s in self.specs().items():
            self.assertEqual({k: v for k, v in s["data_root"].items() if k != "reason"},
                             {"path": "/home/layth/review/grm-h1-data-root-20261107",
                              "manifest": "/home/layth/review/grm-h1data-20261107/stage/MASTER-MANIFEST-h1.sha256",
                              "manifest_sha256": "7dd00779bc4b538b325d947137c525afb4e723ca494bd14b607a8730b7902887"},
                             name)
            self.load(name)

    def test_the_folder_holds_the_pilot_the_four_species_and_their_parts(self):
        parts = [f"h1-{s}-part{n}{k}.json" for s, m in self.MAIN_PARTS.items() for n in range(1, m + 1)
                 for k in ("", *(f"-{k}" for k in self.K_ARMS))]
        self.assertEqual(len(parts), 15)
        self.assertEqual(sorted(self.specs()), sorted(["pilot.json", *(f"h1-{s}.json" for s in self.COUNTS), *parts]))

    def test_each_spec_changes_only_its_species_datasets_kset_size_and_arms(self):
        for name, s in self.specs().items():
            self.assertEqual(set(s), set(SPEC) | {"data_root"}, name)
            for key in set(SPEC) - {"decisions", "species", "datasets", "pool_budget", "arms", "revision12_gate"}:
                self.assertEqual(s[key], SPEC[key], (name, key))
            self.assertEqual({k: v for k, v in s["pool_budget"].items() if k not in ("default_kset_bytes", "reason")},
                             {k: v for k, v in SPEC["pool_budget"].items()
                              if k not in ("default_kset_bytes", "reason")}, name)
            self.assertEqual(s["pool_budget"]["default_kset_bytes"], D2[s["species"]], name)
            self.assertIn(f"{D2[s['species']]:,} B per genome", s["pool_budget"]["reason"], name)
            self.assertIn("M-2", s["decisions"], name)
            self.assertIn("Layth lifted register section 13 G2 on 8 Nov 2026", s["decisions"], name)
            self.assertNotIn("left out", s["decisions"], name)
            self.assertNotIn("synthetic", s, name)
            self.assertEqual(s["arms"], [a for a in SPEC["arms"] if a in s["arms"]], name)
            _, default, arms = self.load(name)
            self.assertEqual(default["k"], 31)
            self.assertEqual([a["name"] for a in arms if a.get("sentinel")], [x["name"] for x in run.SENTINELS])
            for d in s["datasets"]:
                self.assertEqual(d["partition"], "confirmatory", (name, d))
                self.assertTrue(d["id"].startswith(s["species"] + "__"), (name, d))

    def test_each_spec_gates_on_revision_12_and_records_the_h1_species_check(self):
        """Main's ruling after H1DATA: the pins are revision 12's (spec.json's, from its master root); the recorded
        order and view are revision 12's cohort with the H-1 species check (Salmonella c60015cf, view 2a64c19c)."""
        gate = SPEC["revision12_gate"]
        self.assertEqual(sorted(gate["genome_set"]), sorted(D2))
        self.assertEqual({s: len(v) for s, v in gate["datasets"].items()},
                         {"salmonella_enterica": 25, "mycobacterium_tuberculosis": 22, "staphylococcus_aureus": 19,
                          "neisseria_gonorrhoeae": 7})
        self.assertEqual(gate["recorded"]["genome_order"]["salmonella_enterica"][:8], "5864bd5d")
        self.assertEqual(gate["recorded"]["datasets_json_sha256"][:8], "5519c61d")
        self.assertEqual(gate["reference"]["root"], "/home/layth/review/grm-r12-run-20261105/data-root")
        self.assertEqual(sorted(gate["reference"]["files"]), sorted(
            ["curated/datasets.json", "curated/phenotypes.jsonl", *(f"species/{s}/genomes.txt" for s in D2)]))
        for name, s in self.specs().items():
            g = s["revision12_gate"]
            self.assertEqual({k: v for k, v in g.items() if k != "recorded"},
                             {k: v for k, v in gate.items() if k != "recorded"}, name)
            self.assertEqual(g["recorded"]["genome_order"],
                             {**gate["recorded"]["genome_order"], "salmonella_enterica":
                              "c60015cf9e8d69d4e3189be0830dd6ac40b96351974d76c8c4b83bf1e6fdd805"}, name)
            self.assertEqual(g["recorded"]["datasets_json_sha256"],
                             "2a64c19ca796e013897235047c83e0760276ddb1927bd167c556ff8d1ac48af2", name)
            self.assertIn("revision 12's cohort with the H-1 species check", g["recorded"]["note"], name)
            # REVIEW-H1CODE L3: the final staged species check, not round 1's.
            self.assertIn("c2e7c389779e6aa009800a9d327a3b036a14a01153429c62b51cdd7b66c942ef", g["recorded"]["note"], name)
            self.assertNotIn("4af81f94", g["recorded"]["note"], name)

    def test_stage_1_runs_all_of_spec_json_s_arms(self):
        for species, n in self.COUNTS.items():
            s = self.specs()[f"h1-{species}.json"]
            self.assertEqual((s["species"], len(s["datasets"]), s["arms"]), (species, n, SPEC["arms"]), species)
            self.assertEqual(len(s["arms"]), 11)

    def test_each_species_parts_are_main_parts_and_k_add_ons_that_combine_takes(self):
        """The main parts split the species' datasets and run every arm but k21 and k25; each main part has one k
        add-on part per k arm on its datasets. Fed to `check_parts` as finished parts, they pass every check up to
        the master manifest, which this machine doesn't hold."""
        specs = self.specs()
        for species, m in self.MAIN_PARTS.items():
            full_name = f"h1-{species}.json"
            full = specs[full_name]
            with self.subTest(species), tempfile.TemporaryDirectory() as tmp:
                parts, mains = [], {}
                for n in range(1, m + 1):
                    for k in ("", *self.K_ARMS):
                        name = f"h1-{species}-part{n}{'-' + k if k else ''}.json"
                        s = specs[name]
                        self.assertEqual({x: v for x, v in s.items() if x not in ("arms", "datasets")},
                                         {x: v for x, v in full.items() if x not in ("arms", "datasets")}, name)
                        if k:
                            self.assertEqual(s["arms"], [a for a in SPEC["arms"] if a["name"] == k], name)
                            self.assertEqual(s["datasets"], mains[n]["spec"]["datasets"], name)
                            _, _, arms = self.load(name)
                            self.assertEqual({a["factor"] for a in arms if not a.get("sentinel")}, {"k"}, name)
                        else:
                            self.assertEqual(s["arms"], [a for a in SPEC["arms"] if a["name"] not in self.K_ARMS],
                                             name)
                        work = Path(tmp) / name
                        work.mkdir()
                        (work / run.INPUTS).write_text("{}")
                        protocol = run.sha256_bytes(run.canonical(s).encode())
                        part = {"work": work, "spec": s, "raw": [],
                                "inputs": {"grm_sha256": "c" * 64, "vocabulary_sha256": "d" * 64},
                                "manifest": {"protocol_sha256": protocol, "source_commit": "e" * 40,
                                             "binary_sha256": "c" * 64, "synthetic": False,
                                             "revision12_gate": "passed"},
                                "groups": {"sha256": "f" * 64}, "summary": {"eligible_genomes": 5000,
                                                                          "panel_floor": 2001},
                                "reference": None if not k else {
                                    "spec_sha256": mains[n]["manifest"]["protocol_sha256"],
                                    "inputs_sha256": run.file_sha256(mains[n]["work"] / run.INPUTS)}}
                        if not k:
                            mains[n] = part
                        parts.append(part)
                ids = [d["id"] for n in mains for d in mains[n]["spec"]["datasets"]]
                self.assertEqual(sorted(ids), sorted(d["id"] for d in full["datasets"]))
                _, _, arms = self.load(full_name)
                master = Path(tmp) / "MASTER.sha256"
                master.write_text("")
                with self.assertRaises(SystemExit) as e:
                    run.check_parts(self.H1 / full_name, full, arms, parts, master)
                self.assertIn(f"master manifest {master} is {run.file_sha256(master)}, not the "
                              f"{full['data_root']['manifest_sha256']}", str(e.exception))

    def test_the_salmonella_parts_split_its_datasets_13_and_12(self):
        specs = self.specs()
        full = specs["h1-salmonella_enterica.json"]
        parts = [specs[f"h1-salmonella_enterica-part{i}.json"] for i in (1, 2)]
        self.assertEqual([len(p["datasets"]) for p in parts], [13, 12])
        ids = [d["id"] for p in parts for d in p["datasets"]]
        self.assertEqual(sorted(ids), sorted(d["id"] for d in full["datasets"]))
        self.assertEqual(len(set(ids)), len(ids))
        for p in parts:
            self.assertEqual({k: v for k, v in p.items() if k not in ("datasets", "arms")},
                             {k: v for k, v in full.items() if k not in ("datasets", "arms")})

    def test_the_pilot_is_d1_s(self):
        s = self.specs()["pilot.json"]
        self.assertEqual((s["species"], s["datasets"]),
                         ("salmonella_enterica", [{"id": "salmonella_enterica__ciprofloxacin",
                                                   "partition": "confirmatory"}]))
        _, _, arms = self.load("pilot.json")
        self.assertEqual([(a["name"], a["factor"]) for a in arms],
                         [("grr-v1", "selection"), ("k21", "k"), ("k-c-h-p", "blocks"), ("budget-16000", "budget"),
                          ("sentinel-seed-twin", "calibration"), ("sentinel-tuning-twin", "calibration")])
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
    run whole, and run as parts through the real grm: part1 runs every arm on meropenem; part2 runs every arm but
    k21 on ciprofloxacin, and part2-k21, a k add-on (G2(b)), runs k21 on ciprofloxacin against part2's default arm.
    Combine of the three must give the whole run's J records and K report, after timings. part2 goes first, so the
    pool's arm list in K's manifest and summary.json must come from SPEC, not from the first part's."""

    DATASETS = (DS, "salmonella_enterica__ciprofloxacin")

    @classmethod
    def setUpClass(cls):
        out = os.environ.get("GRM_INPUT_E2E_OUT")
        cls.tmp = None if out else tempfile.TemporaryDirectory()
        if cls.tmp:
            # A class cleanup also runs when setUpClass fails, which tearDownClass does not.
            cls.addClassCleanup(cls.tmp.cleanup)
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
        spec["arms"] = [{"name": "cap-600", "cap": 600}, {"name": "k-c", "blocks": ["K", "C"]}, {"name": "k21", "k": 21}]
        spec["datasets"] = [{"id": d, "partition": "confirmatory"} for d in cls.DATASETS]
        # The parts' RAW is this root, and combine takes only the master manifest the spec pins (REVIEW-H1CODE L6).
        cls.master = master_manifest(cls.raw, cls.out / "MASTER.sha256")
        spec["data_root"] = {"path": str(cls.raw), "manifest": str(cls.master),
                             "manifest_sha256": run.file_sha256(cls.master)}
        cls.spec = cls.out / "spec.json"
        cls.spec.write_text(json.dumps(spec, indent=1))
        cls.parts = []
        k21 = [a for a in spec["arms"] if a["name"] == "k21"]
        for i, (d, arms) in enumerate(((cls.DATASETS[0], spec["arms"]),
                                       (cls.DATASETS[1], [a for a in spec["arms"] if a["name"] != "k21"]),
                                       (cls.DATASETS[1], k21))):
            path = cls.out / ("part2-k21.json" if i == 2 else f"part{i + 1}.json")
            path.write_text(json.dumps({**spec, "arms": arms, "datasets": [{"id": d, "partition": "confirmatory"}]},
                                       indent=1))
            cls.parts.append(path)
        cls.whole = cls.run_spec(cls.spec, "whole")
        cls.works = [cls.run_spec(p, p.stem) for p in cls.parts[:2]]
        cls.addon = cls.run_spec(cls.parts[2], cls.parts[2].stem, "--reference-from", str(cls.works[1]))
        cls.combined = cls.out / "combined"
        cls.combine(cls.combined, cls.works[1], cls.addon, cls.works[0])

    @classmethod
    def run_spec(cls, spec, name, *flags):
        work = cls.out / name
        subprocess.run([sys.executable, "-B", str(HERE / "run.py"), str(spec), "--raw", str(cls.raw), "--work",
                        str(work), "--grm", os.environ["GRM_INPUT_E2E_BINARY"], "--source-commit", "0" * 40,
                        "--jobs", os.environ.get("GRM_INPUT_E2E_JOBS", "3"), "--below-panel-floor",
                        "--sigdrug-vocabulary", str(cls.raw / "vocabulary"), *flags], check=True)
        return work

    @classmethod
    def combine(cls, out, *works):
        argv = [sys.executable, "-B", str(HERE / "run.py"), "combine", str(cls.spec), "--master-manifest",
                str(cls.master), "--out", str(out)]
        for w in works:
            argv += ["--work-from", str(w)]
        subprocess.run(argv, check=True)

    def test_both_datasets_are_decided_and_each_part_ran_one(self):
        report = json.loads((self.combined / "report" / "ablation-report.json").read_text())
        # 2 datasets, 3 arms and the two sentinels.
        self.assertEqual(report["contrast_counts"]["planned"], 2 * 5 * 5 * 2)
        self.assertEqual(report["contrast_counts"]["not_produced"], 0, report["not_produced"])
        for work, d in zip((*self.works, self.addon), (*self.DATASETS, self.DATASETS[1])):
            self.assertEqual(sorted(p.name for p in (work / "arms" / "default").iterdir()), [d])
        groups = [json.loads((w / "species-groups.json").read_text()) for w in (self.whole, *self.works, self.addon)]
        self.assertEqual(len({g["sha256"] for g in groups}), 1)

    def test_the_k_add_on_ran_k21_alone_on_part2_s_default_arm(self):
        """G2(b): the add-on counted no species pool and ran neither the default arm nor the sentinels: it took
        part2's, byte for byte, and its k21 output is the whole run's."""
        part2, addon, d = self.works[1], self.addon, self.DATASETS[1]
        self.assertEqual(sorted(p.name for p in (addon / "arms").iterdir()),
                         sorted(["default", "k21", *SENTINEL_NAMES]))
        for name in ("default", *SENTINEL_NAMES):
            self.assertTrue(run.same_folder(addon / "arms" / name / d, part2 / "arms" / name / d), name)
        self.assertEqual([p.name for p in (addon / "roots").glob("pool-*")], [])
        self.assertTrue(any((part2 / "roots").glob("pool-*")))
        self.assertEqual(json.loads((addon / "disk-estimate.json").read_text())["pool_bytes"], 0)
        ref = json.loads((addon / run.REFERENCE).read_text())
        self.assertEqual(ref["work"], str(part2.resolve()))
        self.assertEqual(json.loads((addon / "summary.json").read_text())["reference_from"], ref)
        k21 = json.loads((addon / "arms" / "k21" / d / run.RESUME).read_text())
        self.assertEqual(k21["groups_sha256"], json.loads((part2 / "species-groups.json").read_text())["sha256"])
        whole = json.loads((self.whole / "arms" / "k21" / d / run.RESUME).read_text())
        self.assertEqual(k21["reference_plans"], whole["reference_plans"])
        self.assertEqual(k21["groups_sha256"], whole["groups_sha256"])
        self.assertEqual((addon / "arms" / "k21" / d / "predictions.jsonl").read_bytes(),
                         (self.whole / "arms" / "k21" / d / "predictions.jsonl").read_bytes())

    def test_the_parts_and_the_k_add_on_combined_equal_the_whole_run(self):
        whole, combined = timeless(self.whole), timeless(self.combined)
        for name, data in whole["j"].items():
            self.assertEqual(combined["j"][name], data, name)
        self.assertEqual(combined["manifest"], whole["manifest"])
        self.assertEqual(combined["report"], whole["report"])
        self.assertEqual(combined["md"], whole["md"])
        for path in sorted((self.whole / "arms").rglob("predictions.jsonl")):
            self.assertEqual((self.combined / path.relative_to(self.whole)).read_bytes(), path.read_bytes(), path)
        record = json.loads((self.combined / "combine.json").read_text())
        # combine's order: part2, its add-on (both ciprofloxacin), then part1 (meropenem)
        self.assertEqual([p["datasets"] for p in record["parts"]],
                         [[self.DATASETS[1]], [self.DATASETS[1]], [self.DATASETS[0]]])
        self.assertEqual([p["arms"] for p in record["parts"]],
                         [["default", "cap-600", "k-c", *SENTINEL_NAMES], ["k21"],
                          ["default", "cap-600", "k-c", "k21", *SENTINEL_NAMES]])
        self.assertEqual([p["reference"] and p["reference"]["work"] for p in record["parts"]],
                         [None, str(self.works[1].resolve()), None])
        self.assertEqual(record["spec_sha256"], whole["manifest"]["protocol_sha256"])
        summary, whole_summary = (json.loads((w / "summary.json").read_text()) for w in (self.combined, self.whole))
        self.assertEqual(summary["arms"], whole_summary["arms"])
        # The first part left k21 out, so its own species-groups.json lists fewer arms than the pool was made over.
        self.assertNotIn("k21", json.loads((self.works[1] / "species-groups.json").read_text())["arms"])
        self.assertEqual({k: v for k, v in summary["species_groups"].items() if k != "structure"},
                         {k: v for k, v in whole_summary["species_groups"].items() if k != "structure"})

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
