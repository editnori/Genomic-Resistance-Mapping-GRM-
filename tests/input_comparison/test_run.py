"""Tests for scripts/input-comparison: the spec rules, and with GRM_INPUT_E2E_BINARY the
synthetic end-to-end comparison (all 11 arms plus the default through the real `grm`).

usage: python3 -B -m unittest discover -s tests/input_comparison -t .
       GRM_INPUT_E2E_BINARY=target/release/grm GRM_INPUT_E2E_OUT=DIR python3 -B -m unittest ...
"""
import argparse
import contextlib
import copy
import importlib.util
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

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT / "scripts" / "input-comparison"


def load(name):
    spec = importlib.util.spec_from_file_location(f"input_comparison_{name}", HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


run = load("run")
null_study = load("null_study")
SPEC = json.loads((HERE / "spec.json").read_text())
# The synthetic e2e default's revision-12-gate values (synthetic_root.py defaults, the default curated with the
# species check and the vocabulary), from the real grm: tools/synthpins.py in the lane folder.
SYNTHETIC_GATE = {"genome_order": {
                      "salmonella_enterica": "2299eda73ef338aa2bcd312c58d948f651e3dea9503d3185ba8a3271ce1a87af"},
                  "datasets_json_sha256": "fa73a6223da0f58abc9a8c098b67747ee48889211e5b3504c85417eb13250260"}
SYNTHETIC_CHECK = run.Runner.check_synthetic


def timeless(work):
    """A finished comparison's J records, K manifest and report without what only timing changes: each run's
    seconds and compute.actual_seconds in the manifest (also inside the report), the report's sha256 of that
    manifest, and the Markdown's seconds and input-hash lines."""
    def strip(m):
        for u in m["runs"]:
            u.pop("seconds")
        m["compute"].pop("actual_seconds")
        return m

    work = Path(work)
    report = json.loads((work / "report" / "ablation-report.json").read_text())
    strip(report["manifest"])
    report["inputs_sha256"].pop("manifest")
    md = [line for line in (work / "report" / "ablation-report.md").read_text().splitlines()
          if not line.startswith(("- Actual seconds:", "- Input sha256"))]
    return {"j": {n: (work / "j" / f"{n}.json").read_bytes() for n in ("contrasts", "planned", "missing")},
            "manifest": strip(json.loads((work / "k" / "manifest.json").read_text())), "report": report, "md": md}


class SpecTests(unittest.TestCase):
    def setUp(self):
        # The stub RAWs below are not synthetic_root.py's, so their synthetic stub specs would stop at the
        # synthetic-root check. These tests are about the steps after it; the check runs for real in
        # test_a_synthetic_spec_refuses_a_raw_that_is_not_a_synthetic_root and the e2e.
        patcher = mock.patch.object(run.Runner, "check_synthetic", lambda runner: None)
        patcher.start()
        self.addCleanup(patcher.stop)

    def write(self, spec):
        f = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
        json.dump(spec, f)
        f.close()
        self.addCleanup(os.unlink, f.name)
        return f.name

    def test_the_decided_spec_has_eleven_one_factor_arms_and_the_sentinel(self):
        _, default, arms = run.load_spec(HERE / "spec.json")
        *arms, sentinel = arms
        self.assertEqual((sentinel["name"], sentinel["learner_seed"], sentinel["pairing"], sentinel["sentinel"]),
                         ("sentinel-seed-twin", 1001, "same_plan", True))
        self.assertEqual(len(arms), 11)
        self.assertEqual((default["cap"], default["selection"], default["k"], default["blocks"],
                          default["screened_features"], default["min_class"]),
                         (2000, "legacy", 31, ["K"], 4000, 25))
        by_factor = {}
        for a in arms:
            by_factor.setdefault(a["factor"], []).append(a["name"])
        self.assertEqual(by_factor, {"cohort": ["cap-1000", "cap-4000", "grr-v1"], "k": ["k21", "k25"],
                                     "blocks": ["k-c", "k-h", "k-p", "k-c-h-p"],
                                     "budget": ["budget-1000", "budget-16000"]})
        self.assertEqual({a["name"]: a["pairing"] for a in arms if a["factor"] in ("cohort", "k")},
                         {"cap-1000": "common_outer_test_set_balanced", "cap-4000": "common_outer_test_set_balanced",
                          "grr-v1": "common_outer_test_set_balanced", "k21": "k31_groups", "k25": "k31_groups"})
        self.assertTrue(all(a["pairing"] == "same_plan" for a in arms if a["factor"] in ("blocks", "budget")))
        blocks = {b for a in arms for b in a["blocks"]}
        self.assertEqual(blocks, {"K", "C", "H", "P"})

    def test_an_arm_must_change_exactly_one_input(self):
        for arm in ({"name": "same"}, {"name": "same-cap", "cap": 2000},
                    {"name": "two", "cap": 1000, "k": 21}, {"name": "odd", "min_class": 10}):
            spec = copy.deepcopy(SPEC)
            spec["arms"] = [arm]
            with self.assertRaises(SystemExit, msg=arm["name"]):
                run.load_spec(self.write(spec))
        spec = copy.deepcopy(SPEC)
        spec["arms"].append(dict(spec["arms"][0]))
        with self.assertRaises(SystemExit):
            run.load_spec(self.write(spec))

    def test_the_spec_sets_the_h1_panel_floor(self):
        self.assertEqual(SPEC["panel_floor"]["min_eligible_genomes"], 2001)

    def test_each_cohort_arm_carries_its_kind_and_that_kind_s_base_rate(self):
        _, default, arms = run.load_spec(HERE / "spec.json")
        kinds = {a["name"]: a["kind"] for a in arms if a["factor"] == "cohort"}
        self.assertEqual(kinds, {"cap-1000": "cap-down", "cap-4000": "cap-up", "grr-v1": "selection-swap"})
        spec = {"base_rates": {"cap-up": {"rate": 0.1, "ci95": [0.0, 0.2], "seeds": 5, "source": "twin"}}}
        cohort = {a["name"]: a for a in arms if a["factor"] == "cohort"}
        self.assertEqual(run.base_rate(spec, cohort["cap-4000"]),
                         {"arm": "cap-4000", "kind": "cap-up", "rate": 0.1, "ci95": [0.0, 0.2], "seeds": 5,
                          "source": "twin"})
        self.assertEqual(run.base_rate(spec, cohort["cap-1000"]),
                         {"arm": "cap-1000", "kind": "cap-down", "rate": None, "ci95": None, "seeds": None,
                          "source": "not measured"})
        bad = copy.deepcopy(SPEC)
        bad["base_rates"] = {"cohort": {"rate": 0.14, "ci95": None, "seeds": 5, "source": "one rate for all"}}
        with self.assertRaises(SystemExit):
            run.load_spec(self.write(bad))

    def runner(self, tmp, spec=None):
        grm = Path(tmp) / "grm-not-called"
        if not grm.exists():
            grm.write_text("#!/bin/sh\necho grm-not-called \"$@\" >&2\nexit 99\n")
            grm.chmod(0o755)
        args = argparse.Namespace(grm=str(grm), raw=str(Path(tmp) / "raw"), work=str(Path(tmp) / "work"))
        r = run.Runner(args, spec or SPEC)
        r.default = run.load_spec(HERE / "spec.json")[1]
        return r

    def test_the_pool_refuses_a_genome_of_the_union_without_a_fasta(self):
        """REVIEW-INPUT4 M1: a missing FASTA refuses with its id, never shrinks the pool."""
        with tempfile.TemporaryDirectory() as tmp:
            r = self.runner(tmp)
            (Path(tmp) / "raw" / "fasta").mkdir(parents=True)
            for g in ("g1", "g3"):
                (Path(tmp) / "raw" / "fasta" / f"{g}.fna.gz").write_bytes(b"")
            r.union = lambda arms: ["g1", "g2", "g3", "g4"]
            r.counted_roots = lambda arms: {}
            with self.assertRaises(SystemExit) as e:
                r.pool([r.default])
            self.assertIn("2 of 4 have none: g2, g4", str(e.exception))

    def test_the_pool_refuses_a_structure_that_is_not_the_union(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = self.runner(tmp)
            union = ["g1", "g2", "g3"]
            r.union = lambda arms: union
            r.counted_roots = lambda arms: {}
            root = Path(tmp) / "work" / "roots" / f"pool-{run.sha256_bytes(run.canonical(union).encode())[:12]}"
            (root / "species" / SPEC["species"]).mkdir(parents=True)
            structure = root / "species" / SPEC["species"] / "structure.json"
            structure.write_text(json.dumps({"genomes": ["g1", "g3"], "k": 31}))
            with self.assertRaises(SystemExit) as e:
                r.pool([r.default])
            self.assertIn("holds 2 genomes, not the arms' union of 3; missing g2", str(e.exception))
            structure.write_text(json.dumps({"genomes": union, "k": 21}))
            with self.assertRaises(SystemExit) as e:
                r.pool([r.default])
            self.assertIn("counted at k = 21", str(e.exception))
            structure.write_text(json.dumps({"genomes": union, "k": 31}))
            self.assertEqual(r.pool([r.default]), structure)

    def test_the_pool_s_kset_estimate_refuses_more_than_its_disk_share(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = self.runner(tmp)
            per = SPEC["pool_budget"]["default_kset_bytes"]
            work = SPEC["pool_budget"]["work_files_bytes"]
            ok = r.disk_check(["g"] * 100, free=int((100 * per + work) / 0.8) + 1, measure=False)
            self.assertTrue(ok["ok"], ok)
            self.assertIn("stated default", ok["source"])
            self.assertFalse(r.disk_check(["g"] * 100, free=int((100 * per + work) / 0.8) - 100, measure=False)["ok"])
            # Every data root's ksets count beside the pool's, with the matrix margin on each root and the arm
            # outputs as a share of every root's ksets.
            roots = {"a-k31": {"genomes": 60, "k": 31, "arms": ["default"]},
                     "b-k21": {"genomes": 40, "k": 21, "arms": ["k21"]}}
            both = r.disk_check(["g"] * 100, free=10 ** 15, roots=roots, measure=False)
            self.assertEqual(both["roots"]["a-k31"]["bytes"], round(60 * per * 1.10))
            self.assertEqual(both["arm_output_bytes"], round(100 * per * 0.05))
            self.assertEqual(both["work_files_bytes"], SPEC["pool_budget"]["work_files_bytes"])
            self.assertEqual(both["kset_bytes"], round(100 * per + 100 * per * 1.10 + 100 * per * 0.05
                                                       + both["work_files_bytes"]))
            self.assertFalse(r.disk_check(["g"] * 100, free=int((100 * per + work) / 0.8) + 1, roots=roots,
                                          measure=False)["ok"])
            # A reused WORK: a built root writes nothing but its arm outputs; a root with 3 of 10 ksets finished
            # writes the other 7 and its matrix; the pool's structure.json leaves nothing pending.
            reused = {"built": {"genomes": 10, "k": 31, "arms": ["default"], "built": True, "pending": 0},
                      "partial": {"genomes": 10, "k": 27, "arms": ["k27"], "built": False, "pending": 7}}
            left = r.disk_check(["g"] * 100, free=10 ** 15, roots=reused, measure=False, pool_pending=0)
            self.assertEqual((left["pool_bytes"], left["roots"]["built"]["bytes"]), (0, 0))
            self.assertEqual(left["roots"]["partial"]["bytes"], round(7 * per + 10 * per * 0.10))
            self.assertEqual(left["kset_bytes"], round(7 * per + 10 * per * 0.10 + 20 * per * 0.05
                                                       + left["work_files_bytes"]))
            # On the box the smaller of the given and the measured free disk counts.
            self.assertLessEqual(r.disk_check(["g"], free=10 ** 18)["free_bytes"],
                                 shutil.disk_usage(Path(tmp) / "work").free)
            # The species' own counted ksets in RAW set the per-genome size when there are any.
            ksets = Path(tmp) / "raw" / "species" / SPEC["species"] / "ksets"
            ksets.mkdir(parents=True)
            for i, size in enumerate((1000, 3000)):
                (ksets / f"g{i}.kset").write_bytes(b"x" * size)
            own = r.disk_check(["g"] * 10, free=100000, measure=False)
            self.assertEqual((own["kset_bytes"] - own["work_files_bytes"], own["per_genome_bytes"]), (20000, 2000))
            self.assertIn("2 counted ksets", own["source"])
            r.union = lambda arms: ["g1"]
            r.counted_roots = lambda arms: {}
            r.pool_curated = lambda union: 0
            with self.assertRaises(SystemExit) as e:
                r.check_disk([r.default], free=100)
            self.assertIn("over 80% of", str(e.exception))
            self.assertFalse(json.loads((Path(tmp) / "work" / "disk-estimate.json").read_text())["ok"])

    def test_the_box_takes_the_smaller_of_the_given_and_the_measured_free_disk(self):
        """REVIEW-INPUT6 L2 (R9): with --free-disk-bytes given, the box still reads WORK's volume."""
        with tempfile.TemporaryDirectory() as tmp:
            r = self.runner(tmp)
            r.union = lambda arms: ["g1"]
            r.counted_roots = lambda arms: {}
            r.pool_curated = lambda union: 0
            real = run.shutil.disk_usage
            run.shutil.disk_usage = lambda path: real(path)._replace(free=1000)
            try:
                with self.assertRaises(SystemExit) as e:
                    r.check_disk([r.default], free=10 ** 15)
            finally:
                run.shutil.disk_usage = real
            self.assertIn("of 0.0 GiB free", str(e.exception))
            self.assertEqual(json.loads((Path(tmp) / "work" / "disk-estimate.json").read_text())["free_bytes"], 1000)

    STUB_GRM = """import json, os, sys
from pathlib import Path
# A stand-in for grm curate only: 50 eligible genomes g000..g049, the first --cap of them selected.
root, argv = Path(sys.argv[2]), sys.argv[3:]
if argv[0] == "arm":
    # Records whether an earlier run's plans were still there when the arm started.
    out = Path(argv[argv.index("--out") + 1])
    out.mkdir(parents=True, exist_ok=True)
    (out / "seen.json").write_text(json.dumps({"plans_present": (out / "plans").exists()}))
    sys.exit(0)
if argv[0] in ("count", "structure", "matrix"):
    # Every count is recorded beside WORK; with STUB_COUNT_FAILS set, a count fails the run.
    with open(root.parents[1] / "counted.txt", "a") as f:
        f.write(" ".join(argv) + "\\n")
    sys.exit(3 if os.environ.get("STUB_COUNT_FAILS") else 0)
if argv[0] != "curate":
    sys.exit(f"stub grm: {argv[0]} is not stubbed")
cap = int(argv[argv.index("--cap") + 1])
(root / "curate-env.json").write_text(json.dumps({"GRM_SIGDRUG_VOCABULARY": os.environ.get("GRM_SIGDRUG_VOCABULARY")}))
ids = [f"g{i:03d}" for i in range(50)]
(root / "curated").mkdir(parents=True, exist_ok=True)
(root / "curated" / "datasets.json").write_text(json.dumps([{"species_slug": "salmonella_enterica",
                                                              "antibiotic": "meropenem"}]))
(root / "curated" / "phenotypes.jsonl").write_text("".join(json.dumps({
    "genome_id": g, "antibiotic": "meropenem", "reported": "R" if i % 3 else "S"}) + "\\n" for i, g in enumerate(ids)))
# grm curate's provenance side file, which run.py drops from every root.
(root / "curated" / "provenance.jsonl").write_text("".join(json.dumps({
    "genome_id": g, "antibiotic": "meropenem", "pmids": ["1"]}) + "\\n" for g in ids))
(root / "species" / "salmonella_enterica").mkdir(parents=True, exist_ok=True)
(root / "species" / "salmonella_enterica" / "genomes.txt").write_text("\\n".join(ids[:cap]) + "\\n")
"""

    def stub_spec(self, tmp, **changes):
        """The stub runs' spec: the e2e's arms cut to cap-10, cap-30 and k21 around a default of cap 20, without
        the revision-12 gate and with neither curate input (the stub RAW has neither) unless `changes` sets them."""
        spec = copy.deepcopy(SPEC)
        spec["default"]["cap"] = 20
        spec["arms"] = [{"name": "cap-10", "cap": 10}, {"name": "cap-30", "cap": 30}, {"name": "k21", "k": 21}]
        spec["curate_inputs"] = {"species_check": False, "sigdrug_vocabulary": False}
        spec.pop("revision12_gate")
        spec["synthetic"] = True
        spec.update(changes)
        path = Path(tmp) / "spec.json"
        path.write_text(json.dumps(spec))
        return path

    def stub_grm(self, tmp):
        grm = Path(tmp) / "grm"
        grm.write_text(f"#!{sys.executable}\n" + self.STUB_GRM)
        grm.chmod(0o755)
        return grm

    def main(self, tmp, *flags, spec=None):
        """run.main with the stub grm; returns the SystemExit's code (a refusal's message) and stdout."""
        tmp = Path(tmp)
        spec = spec or self.stub_spec(tmp)
        out = io.StringIO()
        with contextlib.redirect_stdout(out), self.assertRaises(SystemExit) as e:
            run.main([str(spec), "--raw", str(tmp / "raw"), "--work", str(tmp / "work"), "--grm",
                      str(tmp / "grm"), "--source-commit", "0" * 40, "--below-panel-floor", *flags])
        return e.exception.code, out.getvalue()

    def test_main_checks_disk_before_it_counts_anything(self):
        """REVIEW-INPUT7 L1: run.main past --estimate-only, with a free disk that refuses and a grm whose count
        fails: the refusal comes first, and nothing is counted (kills a check after the pool or a root)."""
        with tempfile.TemporaryDirectory() as tmp:
            self.stub_grm(tmp)
            os.environ["STUB_COUNT_FAILS"] = "1"
            try:
                code, _ = self.main(tmp, "--free-disk-bytes", "1000")
            finally:
                del os.environ["STUB_COUNT_FAILS"]
            self.assertIn("the comparison still has to write about", str(code))
            self.assertFalse((Path(tmp) / "work" / "counted.txt").exists())

    def test_the_revision12_gate_refuses_another_default_cohort(self):
        """REVIEW-INPUT7 and REVIEW-INPUT8 M1: before anything is counted, the default's drawn genome order and its
        datasets.json without selection records must be revision 12's; --estimate-only checks it too (L1)."""
        with tempfile.TemporaryDirectory() as tmp:
            self.stub_grm(tmp)
            self.main(tmp, "--free-disk-bytes", "1000")
            root = next(r for r in (Path(tmp) / "work" / "roots").iterdir()
                        if json.loads((r / "curated.json").read_text())["arm"] == "default")
            slug = SPEC["species"]
            order = json.loads((root / "curated.json").read_text())["drawn_order_sha256"]
            view = run.datasets_view(root / "curated" / "datasets.json")
            ours = {"genome_order": {slug: order}, "datasets_json_sha256": view}

            def late(h):
                # The right value with only its last hex digit changed (REVIEW-INPUT9 L2, g13).
                return h[:-1] + ("0" if h[-1] != "0" else "1")

            for wrong, message in (
                    ({**ours, "genome_order": {slug: late(order)}},
                     f"the default arm's genome order for {slug} is {order}, not revision 12's {late(order)}"),
                    ({**ours, "datasets_json_sha256": late(view)},
                     f"the default arm's datasets.json is {view}, not revision 12's {late(view)}")):
                for flags in (("--free-disk-bytes", "1000"), ("--free-disk-bytes", str(10 ** 15), "--estimate-only")):
                    code, _ = self.main(tmp, *flags, spec=self.stub_spec(tmp, revision12_gate=wrong))
                    self.assertIn(message, str(code))
            # The run's species must be named; every named species is compared with its own root file.
            code, _ = self.main(tmp, "--free-disk-bytes", "1000",
                                spec=self.stub_spec(tmp, revision12_gate={**ours, "genome_order": {"other_sp": order}}))
            self.assertIn(f"names no genome order for {slug}", str(code))
            code, _ = self.main(tmp, "--free-disk-bytes", "1000", spec=self.stub_spec(
                tmp, revision12_gate={**ours, "genome_order": {slug: order, "other_sp": order}}))
            self.assertIn("has no species/other_sp/genomes.txt", str(code))
            (root / "species" / "other_sp").mkdir()
            (root / "species" / "other_sp" / "genomes.txt").write_text("x1\nx2\n")
            other = run.sha256_bytes(b"x1\nx2\n")
            code, _ = self.main(tmp, "--free-disk-bytes", "1000", spec=self.stub_spec(
                tmp, revision12_gate={**ours, "genome_order": {slug: order, "other_sp": late(other)}}))
            self.assertIn(f"genome order for other_sp is {other}, not revision 12's {late(other)}", str(code))
            self.assertFalse((Path(tmp) / "work" / "counted.txt").exists())
            ours = {**ours, "genome_order": {slug: order, "other_sp": other}}
            code, _ = self.main(tmp, "--free-disk-bytes", "1000", spec=self.stub_spec(tmp, revision12_gate=ours))
            self.assertIn("the comparison still has to write about", str(code))
            code, out = self.main(tmp, "--free-disk-bytes", str(10 ** 15), "--estimate-only",
                                  spec=self.stub_spec(tmp, revision12_gate=ours))
            self.assertEqual((code, json.loads(out)["revision12_gate"]), (0, ours))

    def test_a_real_run_without_the_gate_refuses(self):
        """REVIEW-INPUT9 L1: no revision12_gate passes only when the spec says "synthetic": true."""
        with tempfile.TemporaryDirectory() as tmp:
            self.stub_grm(tmp)
            for flags in (("--free-disk-bytes", "1000"), ("--free-disk-bytes", str(10 ** 15), "--estimate-only")):
                for synthetic in (False, None, "true"):
                    spec = json.loads(self.stub_spec(tmp).read_text())
                    if synthetic is None:
                        spec.pop("synthetic")
                    else:
                        spec["synthetic"] = synthetic
                    (Path(tmp) / "spec.json").write_text(json.dumps(spec))
                    code, _ = self.main(tmp, *flags, spec=Path(tmp) / "spec.json")
                    self.assertIn("a real run is never ungated", str(code), (flags, synthetic))
            code, out = self.main(tmp, "--free-disk-bytes", str(10 ** 15), "--estimate-only")
            self.assertEqual((code, json.loads(out)["revision12_gate"]), (0, None))

    def synthetic_raw(self, tmp, ids=None, fmt=run.SYNTHETIC_FORMAT, truth=True):
        """A stub RAW whose bvbrc/genomes.jsonl holds the stub's 50 genomes, with a truth.json of format `fmt`
        listing `ids` (those 50 by default), or none."""
        raw = Path(tmp) / "raw"
        held = [f"g{i:03d}" for i in range(50)]
        (raw / "bvbrc").mkdir(parents=True, exist_ok=True)
        (raw / "bvbrc" / "genomes.jsonl").write_text("".join(json.dumps({"genome_id": g}) + "\n" for g in held))
        marker = raw / "truth.json"
        if truth:
            marker.write_text(json.dumps({"format": fmt, "genomes": [{"genome_id": g} for g in (ids or held)]}))
        elif marker.exists():
            marker.unlink()

    def test_a_synthetic_spec_refuses_a_raw_that_is_not_a_synthetic_root(self):
        """REVIEW-INTEG7 L2: "synthetic": true runs only on a RAW whose truth.json is synthetic_root.py's and lists
        exactly RAW's genomes; otherwise it refuses in both modes, before any root is curated."""
        with mock.patch.object(run.Runner, "check_synthetic", SYNTHETIC_CHECK), tempfile.TemporaryDirectory() as tmp:
            self.stub_grm(tmp)
            cases = (({"truth": False}, "no readable truth.json"),
                     ({"fmt": "grm-ablation-synthetic"}, f"truth.json's format is not {run.SYNTHETIC_FORMAT}"),
                     # A marker copied from another synthetic root: right format, other genomes.
                     ({"ids": [f"x{i:03d}" for i in range(50)]},
                      "truth.json lists 50 genomes and bvbrc/genomes.jsonl holds 50, 100 of them in only one"),
                     ({"ids": [f"g{i:03d}" for i in range(49)]}, "1 of them in only one"))
            for kwargs, message in cases:
                self.synthetic_raw(tmp, **kwargs)
                for flags in (("--free-disk-bytes", "1000"), ("--free-disk-bytes", str(10 ** 15), "--estimate-only")):
                    code, _ = self.main(tmp, *flags)
                    self.assertIn("is not a synthetic root", str(code), (kwargs, flags))
                    self.assertIn(message, str(code), (kwargs, flags))
                    self.assertFalse((Path(tmp) / "work" / "roots").exists(), (kwargs, flags))
                    # REVIEW-INTEG7B L1: the check comes before RAW is hashed (check_inputs writes inputs.json).
                    self.assertFalse((Path(tmp) / "work" / run.INPUTS).exists(), (kwargs, flags))
            # A spec without "synthetic" isn't checked here; it needs its gate instead.
            spec = self.stub_spec(tmp)
            plain = json.loads(spec.read_text())
            plain.pop("synthetic")
            spec.write_text(json.dumps(plain))
            self.assertIn("a real run is never ungated", str(self.main(tmp, "--free-disk-bytes", "1000", spec=spec)[0]))
            self.synthetic_raw(tmp)
            code, out = self.main(tmp, "--free-disk-bytes", str(10 ** 15), "--estimate-only")
            self.assertEqual((code, json.loads(out)["revision12_gate"]), (0, None))

    def test_the_gate_decision_is_printed_on_stderr(self):
        """REVIEW-INTEG7 L2: the run says at its start, on stderr, whether the gate passed or wasn't required."""
        with mock.patch.object(run.Runner, "check_synthetic", SYNTHETIC_CHECK), tempfile.TemporaryDirectory() as tmp:
            self.stub_grm(tmp)
            self.synthetic_raw(tmp)
            err = io.StringIO()
            with contextlib.redirect_stderr(err):
                code, _ = self.main(tmp, "--free-disk-bytes", str(10 ** 15), "--estimate-only")
            self.assertEqual(code, 0)
            self.assertIn("revision-12 gate: NOT REQUIRED: the spec says \"synthetic\": true and RAW is a synthetic "
                          "root (truth.json); this run is ungated", err.getvalue())
            root = next(r for r in (Path(tmp) / "work" / "roots").iterdir()
                        if json.loads((r / "curated.json").read_text())["arm"] == "default")
            slug = SPEC["species"]
            gate = {"genome_order": {slug: json.loads((root / "curated.json").read_text())["drawn_order_sha256"]},
                    "datasets_json_sha256": run.datasets_view(root / "curated" / "datasets.json")}
            err = io.StringIO()
            with contextlib.redirect_stderr(err):
                code, _ = self.main(tmp, "--free-disk-bytes", str(10 ** 15), "--estimate-only",
                                    spec=self.stub_spec(tmp, revision12_gate=gate))
            self.assertEqual(code, 0)
            self.assertIn(f"revision-12 gate: PASSED: the default's genome order for {slug} and its datasets.json "
                          f"view {gate['datasets_json_sha256']} equal revision 12's (the spec also says synthetic; "
                          "the gate ran)", err.getvalue())
            self.assertNotIn("NOT REQUIRED", err.getvalue())
        # A real spec's line names no synthetic flag.
        self.assertTrue(run.gate_line(gate, {}).endswith("equal revision 12's"))

    def test_the_manifest_records_whether_the_gate_ran(self):
        """REVIEW-INTEG7B L2: an ungated synthetic run's manifest records revision12_gate "not_required" and a gated
        one "passed". run.main runs past the manifest with the counting, arm and contrast steps stubbed, and stops
        at K's report."""
        def manifest_after(tmp, spec=None):
            pool = Path(tmp) / "structure.json"
            pool.write_text(json.dumps({"genomes": [f"g{i:03d}" for i in range(50)], "lineages": {"n_clusters": 1}}))
            with mock.patch.object(run.Runner, "check_synthetic", SYNTHETIC_CHECK), \
                    mock.patch.object(run.Runner, "check_disk", lambda *a: None), \
                    mock.patch.object(run.Runner, "pool", lambda *a: pool), \
                    mock.patch.object(run.Runner, "root", lambda *a: None), \
                    mock.patch.object(run.Runner, "arm", lambda *a: None), \
                    mock.patch.object(run, "contrasts", lambda *a: ([], [], {})), \
                    mock.patch.object(run, "manifest", lambda *a: {}), \
                    mock.patch.object(run.k_report, "build_report", side_effect=SystemExit("stopped after the manifest")), \
                    contextlib.redirect_stderr(io.StringIO()):
                code, _ = self.main(tmp, "--free-disk-bytes", str(10 ** 15), spec=spec)
            self.assertEqual(code, "stopped after the manifest")
            m = json.loads((Path(tmp) / "work" / "k" / "manifest.json").read_text())
            return m["synthetic"], m["revision12_gate"]

        with tempfile.TemporaryDirectory() as tmp:
            self.stub_grm(tmp)
            self.synthetic_raw(tmp)
            self.assertEqual(manifest_after(tmp), (True, "not_required"))
            root = next(r for r in (Path(tmp) / "work" / "roots").iterdir()
                        if json.loads((r / "curated.json").read_text())["arm"] == "default")
            gate = {"genome_order": {SPEC["species"]: json.loads((root / "curated.json").read_text())["drawn_order_sha256"]},
                    "datasets_json_sha256": run.datasets_view(root / "curated" / "datasets.json")}
            self.assertEqual(manifest_after(tmp, self.stub_spec(tmp, revision12_gate=gate)), (True, "passed"))

    def test_the_gate_s_datasets_view_drops_only_selection_records(self):
        """REVIEW-INPUT8 M1: this tree's curate adds `selection` to every dataset; the view serializes the rest as
        grm writes it, so an older datasets.json hashes to its own bytes."""
        with tempfile.TemporaryDirectory() as tmp:
            old = [{"id": "s__d", "n_resistant": 3, "breakpoint": None, "min_class": 25}]
            text = json.dumps(old, indent=2)
            (Path(tmp) / "old.json").write_text(text)
            self.assertEqual(run.datasets_view(Path(tmp) / "old.json"), run.sha256_bytes(text.encode()))
            new = [{**old[0], "selection": {"rule": "legacy", "cap": 2000, "seed": None, "always_included": None}}]
            (Path(tmp) / "new.json").write_text(json.dumps(new, indent=2))
            self.assertEqual(run.datasets_view(Path(tmp) / "new.json"), run.sha256_bytes(text.encode()))
            (Path(tmp) / "other.json").write_text(json.dumps([{**old[0], "n_resistant": 4}], indent=2))
            self.assertNotEqual(run.datasets_view(Path(tmp) / "other.json"), run.sha256_bytes(text.encode()))
            # Every entry, not only the first, loses its selection record (REVIEW-INPUT9 L2, g9).
            many = [{**old[0], "id": f"s__d{i}"} for i in range(3)]
            many_text = json.dumps(many, indent=2)
            (Path(tmp) / "many.json").write_text(json.dumps([{**e, "selection": new[0]["selection"]} for e in many],
                                                           indent=2))
            self.assertEqual(run.datasets_view(Path(tmp) / "many.json"), run.sha256_bytes(many_text.encode()))

    def test_a_spec_must_state_its_curate_inputs(self):
        """REVIEW-INPUT8 L3."""
        with tempfile.TemporaryDirectory() as tmp:
            self.stub_grm(tmp)
            spec = json.loads(self.stub_spec(tmp).read_text())
            spec.pop("curate_inputs")
            (Path(tmp) / "spec.json").write_text(json.dumps(spec))
            self.assertIn("states no curate_inputs", str(self.main(tmp, spec=Path(tmp) / "spec.json")[0]))

    def test_grm_children_never_see_the_vocabulary_variable(self):
        """REVIEW-INPUT8 L3: GRM_SIGDRUG_VOCABULARY would give curate a vocabulary inputs.json doesn't hash."""
        with tempfile.TemporaryDirectory() as tmp:
            self.stub_grm(tmp)
            os.environ["GRM_SIGDRUG_VOCABULARY"] = str(Path(tmp) / "a-vocabulary")
            try:
                self.main(tmp, "--free-disk-bytes", str(10 ** 15), "--estimate-only")
            finally:
                del os.environ["GRM_SIGDRUG_VOCABULARY"]
            envs = list((Path(tmp) / "work" / "roots").glob("*/curate-env.json"))
            self.assertTrue(envs)
            for e in envs:
                self.assertEqual(json.loads(e.read_text()), {"GRM_SIGDRUG_VOCABULARY": None}, e)

    def test_curate_inputs_are_required_when_the_spec_asks(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.stub_grm(tmp)
            spec = self.stub_spec(tmp, curate_inputs={"species_check": True, "sigdrug_vocabulary": False})
            self.assertIn("has no species-check/species_check.json", str(self.main(tmp, spec=spec)[0]))
            (Path(tmp) / "raw" / "species-check").mkdir(parents=True)
            (Path(tmp) / "raw" / "species-check" / "species_check.json").write_text("{}")
            spec = self.stub_spec(tmp, curate_inputs={"species_check": True, "sigdrug_vocabulary": True})
            self.assertIn("pass --sigdrug-vocabulary", str(self.main(tmp, spec=spec)[0]))

    def test_arm_and_dataset_names_are_plain_and_unique(self):
        """REVIEW-INPUT7 L2 and REVIEW-INPUT8 I2: names become folders under WORK; whole-string match, at most 64."""
        spec = copy.deepcopy(SPEC)
        spec["arms"] = [{"name": "c" * 64, "cap": 1000}]
        run.load_spec(self.write(spec))
        for bad in ("../../escape", "K+H", "a/b", ".hidden", "", "cap-1\n", "c" * 65):
            spec = copy.deepcopy(SPEC)
            spec["arms"] = [{"name": bad, "cap": 1000}]
            with self.assertRaises(SystemExit, msg=bad):
                run.load_spec(self.write(spec))
            spec = copy.deepcopy(SPEC)
            spec["datasets"] = [{"id": bad, "partition": "confirmatory"}]
            with self.assertRaises(SystemExit, msg=bad):
                run.load_spec(self.write(spec))
        spec = copy.deepcopy(SPEC)
        spec["datasets"] = spec["datasets"] * 2
        with self.assertRaises(SystemExit) as e:
            run.load_spec(self.write(spec))
        self.assertIn("listed twice", str(e.exception))

    def test_an_arm_folder_outside_work_or_a_linked_plans_folder_is_refused(self):
        """REVIEW-INPUT7 L2: arm() never clears a folder outside WORK/arms, and never skips a plans link."""
        ds = "salmonella_enterica__meropenem"
        with tempfile.TemporaryDirectory() as tmp:
            r = self.runner(tmp)
            outside = Path(tmp) / "escape" / ds / "plans"
            outside.mkdir(parents=True)
            (Path(tmp) / "work" / "arms").mkdir(parents=True)
            with self.assertRaises(SystemExit) as e:
                r.arm({**r.default, "name": "../../escape", "factor": "blocks"}, ds)
            self.assertIn("is outside", str(e.exception))
            self.assertTrue(outside.exists())
            target = Path(tmp) / "elsewhere"
            (target / "lineage-1.json").parent.mkdir(parents=True)
            (target / "lineage-1.json").write_text("{}")
            out = Path(tmp) / "work" / "arms" / "cap-10" / ds
            out.mkdir(parents=True)
            (out / "plans").symlink_to(target, target_is_directory=True)
            with self.assertRaises(SystemExit) as e:
                r.arm({**r.default, "name": "cap-10", "factor": "cohort"}, ds)
            self.assertIn("is a symlink", str(e.exception))
            self.assertTrue((target / "lineage-1.json").exists())

    def test_an_identical_cohort_s_old_plans_are_cleared_too(self):
        """REVIEW-INPUT7 L2 (C2): the stub selects grr-v1's cohort as the default's, so it is an identical cohort."""
        ds = "salmonella_enterica__meropenem"
        with tempfile.TemporaryDirectory() as tmp:
            self.stub_grm(tmp)
            r = self.runner(tmp)
            r.grm = str(Path(tmp) / "grm")
            r.default = {**r.default, "cap": 20}
            arm = {**r.default, "name": "grr-v1", "factor": "cohort", "selection": "group-round-robin-v1",
                   "selection_seed": 1, "pairing": "common_outer_test_set_balanced"}
            plans = Path(tmp) / "work" / "arms" / "grr-v1" / ds / "plans"
            plans.mkdir(parents=True)
            (plans / "lineage-1.json").write_text(json.dumps({"plan": {"unavailable_reason": "stale"}}))
            out = r.arm(arm, ds)
            self.assertIn("identical_cohort", (out / "arm.json").read_text())
            self.assertFalse(plans.exists())

    def test_a_symlink_anywhere_in_raw_is_refused(self):
        """REVIEW-INPUT8 L5: a linked top-level part, a linked file inside a part, and a link outside the hashed
        parts each refuse."""
        def raw_with(tmp):
            raw = Path(tmp) / "raw"
            for p in ("bvbrc/genomes.jsonl", "fasta/g000.fna.gz", "species/s/genomes.txt"):
                (raw / p).parent.mkdir(parents=True, exist_ok=True)
                (raw / p).write_text("one")
            return raw

        for case in ("part", "file", "outside", "outside-dir"):
            with tempfile.TemporaryDirectory() as tmp:
                raw = raw_with(tmp)
                elsewhere = Path(tmp) / "elsewhere"
                elsewhere.mkdir()
                (elsewhere / "x.json").write_text("x")
                if case == "part":
                    (raw / "catalog").symlink_to(elsewhere, target_is_directory=True)
                elif case == "file":
                    (raw / "fasta" / "g001.fna.gz").symlink_to(elsewhere / "x.json")
                elif case == "outside":
                    (raw / "species" / "s" / "linked.txt").symlink_to(elsewhere / "x.json")
                else:
                    # A linked folder outside the hashed parts (REVIEW-INPUT9 L2).
                    (raw / "species" / "linked").symlink_to(elsewhere, target_is_directory=True)
                with self.assertRaises(AssertionError, msg=case) as e:
                    self.estimate(tmp, "--free-disk-bytes", str(10 ** 15))
                self.assertIn("symlink", str(e.exception), case)

    def test_a_stratum_file_is_hashed_and_must_lie_inside_raw(self):
        """REVIEW-INPUT8 L5 (H5): a changed stratum refuses a reuse; a stratum outside RAW refuses."""
        with tempfile.TemporaryDirectory() as tmp:
            raw = Path(tmp) / "raw"
            (raw / "strata").mkdir(parents=True)
            (raw / "strata" / "s.json").write_text("[]")
            self.stub_grm(tmp)
            arms = [{"name": "cap-30", "cap": 30},
                    {"name": "grr-s", "selection": "group-round-robin-v1", "selection_seed": 1,
                     "stratum": "strata/s.json"}]
            spec = json.loads(self.stub_spec(tmp).read_text())
            spec["arms"] = arms
            path = Path(tmp) / "spec.json"
            path.write_text(json.dumps(spec))
            self.assertEqual(self.main(tmp, "--free-disk-bytes", str(10 ** 15), "--estimate-only", spec=path)[0], 0)
            (raw / "strata" / "s.json").write_text('["g000"]')
            self.assertIn("another raw_manifest_sha256",
                          str(self.main(tmp, "--free-disk-bytes", str(10 ** 15), "--estimate-only", spec=path)[0]))
            spec["arms"][1]["stratum"] = "../outside.json"
            (Path(tmp) / "outside.json").write_text("[]")
            path.write_text(json.dumps(spec))
            shutil.rmtree(Path(tmp) / "work")
            self.assertIn("is outside RAW",
                          str(self.main(tmp, "--free-disk-bytes", str(10 ** 15), "--estimate-only", spec=path)[0]))

    def test_the_input_record_covers_non_fasta_inputs_and_refuses_links(self):
        """REVIEW-INPUT7 L3: a changed annotation, species check or vocabulary file refuses a reuse; a symlink
        inside RAW refuses outright."""
        with tempfile.TemporaryDirectory() as tmp:
            raw = Path(tmp) / "raw"
            for p in ("bvbrc/features/g000.json", "species-check/species_check.json", "fasta/g000.fna.gz"):
                (raw / p).parent.mkdir(parents=True, exist_ok=True)
                (raw / p).write_text("one")
            vocab = Path(tmp) / "vocabulary"
            vocab.mkdir()
            (vocab / "drug_name_typo_folds.tsv").write_text("observed\tcanonical\n")
            flags = ("--free-disk-bytes", str(10 ** 15), "--sigdrug-vocabulary", str(vocab))
            self.assertEqual(self.estimate(tmp, *flags)[0], 0)
            record = json.loads((Path(tmp) / "work" / "inputs.json").read_text())
            self.assertEqual(record["raw_files"], 3)
            self.assertIsNotNone(record["vocabulary_sha256"])
            for path, label in ((raw / "bvbrc/features/g000.json", "raw_manifest_sha256"),
                                (raw / "species-check/species_check.json", "raw_manifest_sha256"),
                                (vocab / "drug_name_typo_folds.tsv", "vocabulary_sha256")):
                path.write_text("two")
                with self.assertRaises(AssertionError, msg=path) as e:
                    self.estimate(tmp, *flags)
                self.assertIn(f"another {label}", str(e.exception))
                path.write_text("one" if label != "vocabulary_sha256" else "observed\tcanonical\n")
            self.assertEqual(self.estimate(tmp, *flags)[0], 0)
            (raw / "bvbrc" / "linked").symlink_to(raw / "fasta", target_is_directory=True)
            with self.assertRaises(AssertionError) as e:
                self.estimate(tmp, *flags)
            self.assertIn("holds a symlink", str(e.exception))

    def test_a_bare_grm_name_is_hashed_as_path_finds_it(self):
        """REVIEW-INPUT7 L3: --grm NAME is looked up on PATH, as subprocess runs it, not read from the cwd."""
        with tempfile.TemporaryDirectory() as tmp:
            bin_dir = Path(tmp) / "bin"
            bin_dir.mkdir()
            grm = bin_dir / "grm-stub"
            grm.write_text(f"#!{sys.executable}\n" + self.STUB_GRM)
            grm.chmod(0o755)
            decoy = Path(tmp) / "grm-stub"
            decoy.write_text("#!/bin/sh\nexit 9\n")
            decoy.chmod(0o755)
            path, cwd = os.environ["PATH"], os.getcwd()
            os.environ["PATH"] = f"{bin_dir}{os.pathsep}{path}"
            os.chdir(tmp)
            try:
                r = run.Runner(argparse.Namespace(grm="grm-stub", raw=str(Path(tmp) / "raw"),
                                                  work=str(Path(tmp) / "work")), SPEC)
            finally:
                os.environ["PATH"] = path
                os.chdir(cwd)
            self.assertEqual(r.grm, str(grm.resolve()))
            r.all_arms = [r.default if hasattr(r, "default") else run.load_spec(HERE / "spec.json")[1]]
            self.assertEqual(r.check_inputs()["grm_sha256"], run.file_sha256(grm))

    def estimate(self, tmp, *flags, stub=True):
        """run.main --estimate-only with the stub grm, on a spec whose union (cap 30) is smaller than the 50
        eligible genomes; returns (exit code, the printed estimate). A refusal message fails an assertion."""
        tmp = Path(tmp)
        grm = tmp / "grm"
        if stub:
            grm.write_text(f"#!{sys.executable}\n" + self.STUB_GRM)
        grm.chmod(0o755)
        path = self.stub_spec(tmp)
        out = io.StringIO()
        with contextlib.redirect_stdout(out), self.assertRaises(SystemExit) as e:
            run.main([str(path), "--raw", str(tmp / "raw"), "--work", str(tmp / "work"), "--grm", str(grm),
                      "--source-commit", "0" * 40, "--below-panel-floor", "--estimate-only", *flags])
        if isinstance(e.exception.code, str):
            raise AssertionError(e.exception.code)
        return e.exception.code, json.loads(out.getvalue())

    def test_estimate_only_sizes_the_box_it_is_given_and_exits_1_when_it_does_not_fit(self):
        """REVIEW-INPUT5 L1: --estimate-only uses --free-disk-bytes, not this machine's disk (Mh), and exits 1
        when the estimate exceeds 80% of it (Mi)."""
        per = SPEC["pool_budget"]["default_kset_bytes"]
        with tempfile.TemporaryDirectory() as tmp:
            code, disk = self.estimate(tmp, "--free-disk-bytes", str(10 ** 15))
        self.assertEqual((code, disk["free_bytes"], disk["ok"]), (0, 10 ** 15, True))
        with tempfile.TemporaryDirectory() as tmp:
            code, disk = self.estimate(tmp, "--free-disk-bytes", "1000")
        self.assertEqual((code, disk["free_bytes"], disk["ok"]), (1, 1000, False))
        # 30 pool genomes, and the default's (20, shared with the sentinel), cap-10's, cap-30's and k21's roots.
        self.assertEqual(disk["pool_bytes"], round(30 * per))
        self.assertEqual({r["genomes"] for r in disk["roots"].values()}, {20, 10, 30})
        self.assertEqual(sorted(len(r["arms"]) for r in disk["roots"].values()), [1, 1, 1, 2])
        self.assertEqual(disk["kset_bytes"], round(30 * per) + sum(r["bytes"] for r in disk["roots"].values())
                         + disk["arm_output_bytes"] + disk["pool_curated_bytes"] + disk["work_files_bytes"])
        self.assertGreater(disk["pool_curated_bytes"], 0)

    def test_a_reused_work_refuses_another_raw_or_another_grm(self):
        """REVIEW-INPUT6 L3: WORK records RAW's input manifest and the grm binary, and refuses a reuse when
        either changed, or when it has roots but no record."""
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "raw" / "fasta").mkdir(parents=True)
            fasta = Path(tmp) / "raw" / "fasta" / "g000.fna.gz"
            fasta.write_bytes(b"one")
            self.assertEqual(self.estimate(tmp, "--free-disk-bytes", str(10 ** 15))[0], 0)
            record = json.loads((Path(tmp) / "work" / "inputs.json").read_text())
            self.assertEqual(record["raw_files"], 1)
            self.assertEqual(self.estimate(tmp, "--free-disk-bytes", str(10 ** 15))[0], 0)
            fasta.write_bytes(b"two")
            with self.assertRaises(AssertionError) as e:
                self.estimate(tmp, "--free-disk-bytes", str(10 ** 15))
            self.assertIn("another raw_manifest_sha256", str(e.exception))
            fasta.write_bytes(b"one")
            (Path(tmp) / "grm").write_text((Path(tmp) / "grm").read_text() + "\n# another build\n")
            with self.assertRaises(AssertionError) as e:
                self.estimate(tmp, "--free-disk-bytes", str(10 ** 15), stub=False)
            self.assertIn("another grm_sha256", str(e.exception))
            (Path(tmp) / "work" / "inputs.json").unlink()
            with self.assertRaises(AssertionError) as e:
                self.estimate(tmp, "--free-disk-bytes", str(10 ** 15))
            self.assertIn("has data roots but no inputs.json", str(e.exception))

    def test_an_arm_s_plans_from_an_earlier_run_are_cleared_before_it_runs(self):
        """REVIEW-INPUT6 L3: a stale plan in a current arm's out/plans never reaches this run."""
        ds = "salmonella_enterica__meropenem"
        with tempfile.TemporaryDirectory() as tmp:
            self.estimate(tmp, "--free-disk-bytes", str(10 ** 15))
            r = self.runner(tmp)
            r.grm, r.groups = str(Path(tmp) / "grm"), Path(tmp) / "structure.json"
            arm = next(a for a in run.load_spec(Path(tmp) / "spec.json")[2] if a["name"] == "cap-10")
            plans = Path(tmp) / "work" / "arms" / "cap-10" / ds / "plans"
            plans.mkdir(parents=True)
            (plans / "lineage-1.json").write_text(json.dumps({"plan": {"unavailable_reason": "stale"}}))
            out = r.arm(arm, ds)
            self.assertEqual(json.loads((out / "seen.json").read_text()), {"plans_present": False})
            self.assertFalse(plans.exists())

    def test_curate_s_provenance_file_leaves_every_root_and_is_never_copied_or_sized(self):
        """PROV round 2: grm curate's curated/provenance.jsonl (about 310 MB on the real root, which no arm reads)
        leaves each root once it curates; a root curated before that rule still has one, and the pool neither
        copies it nor counts it in the disk estimate."""
        with tempfile.TemporaryDirectory() as tmp:
            self.estimate(tmp, "--free-disk-bytes", str(10 ** 15))
            curated = sorted((Path(tmp) / "work" / "roots").glob("*/curated"))
            self.assertEqual(len(curated), 4)
            for c in curated:
                self.assertTrue((c / "phenotypes.jsonl").exists(), c)
                self.assertFalse((c / "provenance.jsonl").exists(), c)
            r = self.runner(tmp, spec=json.loads(self.stub_spec(tmp).read_text()))
            r.grm = str(Path(tmp) / "grm")
            default = r.curated(r.default) / "curated"
            size = run.disk_bytes(default)
            (default / "provenance.jsonl").write_bytes(b"x" * 1_000_000)
            union = ["g000", "g001"]
            self.assertEqual(r.pool_curated(union), size)
            (Path(tmp) / "raw" / "fasta").mkdir(parents=True, exist_ok=True)
            for g in union:
                (Path(tmp) / "raw" / "fasta" / f"{g}.fna.gz").write_bytes(b"")
            r.union = lambda arms: union

            def call(name, root, *argv):
                if argv[0] == "structure":
                    path = root / "species" / r.slug / "structure.json"
                    path.write_text(json.dumps({"genomes": union, "k": r.default["k"]}))
            r.call = call
            pool = r.pool([r.default]).parents[2] / "curated"
            self.assertEqual(sorted(p.name for p in pool.iterdir()),
                             sorted(p.name for p in default.iterdir() if p.name != "provenance.jsonl"))

    def test_the_genome_dates_side_files_stay_out_of_raw_s_input_manifest(self):
        """PROV round 2: `grm fetch genome-dates` writing into RAW's bvbrc/ doesn't change RAW's input manifest,
        so a WORK made before the fetch is still reused; any other new file under bvbrc/ still refuses."""
        with tempfile.TemporaryDirectory() as tmp:
            raw = Path(tmp) / "raw"
            (raw / "bvbrc").mkdir(parents=True)
            (raw / "bvbrc" / "genomes.jsonl").write_text("one")
            flags = ("--free-disk-bytes", str(10 ** 15))
            self.assertEqual(self.estimate(tmp, *flags)[0], 0)
            was = json.loads((Path(tmp) / "work" / "inputs.json").read_text())
            self.assertEqual(was["raw_files"], 1)
            self.assertEqual(run.RAW_SIDE_FILES, ("bvbrc/genome_dates.jsonl", "bvbrc/genome_dates.manifest.json"))
            for name in run.RAW_SIDE_FILES:
                (raw / name).write_text("{}")
            self.assertEqual(self.estimate(tmp, *flags)[0], 0)
            self.assertEqual(json.loads((Path(tmp) / "work" / "inputs.json").read_text()), was)
            (raw / "bvbrc" / "genome_dates.jsonl.bak").write_text("{}")
            with self.assertRaises(AssertionError) as e:
                self.estimate(tmp, *flags)
            self.assertIn("another raw_manifest_sha256", str(e.exception))

    def test_the_pool_is_the_arms_union_not_the_eligible_set(self):
        """REVIEW-INPUT5 L2: 50 genomes are eligible, the arms select 30 (Me)."""
        with tempfile.TemporaryDirectory() as tmp:
            code, disk = self.estimate(tmp, "--free-disk-bytes", str(10 ** 15))
            r = self.runner(tmp)
            r.grm = str(Path(tmp) / "grm")
            self.assertEqual(len(r.eligible(r.default)), 50)
        self.assertEqual(disk["genomes"], 30)

    def test_a_salmonella_sized_comparison_refuses_on_one_xlarge(self):
        """REVIEW-INPUT5 M1: union 4,000-6,000 at 27.5 MiB per k31 kset, the spec's roots, an xlarge's 259.3 GiB."""
        roots = {"default": 2000, "cap-1000": 1000, "cap-4000": 4000, "grr-v1": 2000, "k21": 2000, "k25": 2000}
        roots = {name: {"genomes": n, "k": 31, "arms": [name]} for name, n in roots.items()}
        with tempfile.TemporaryDirectory() as tmp:
            r = self.runner(tmp)
            for union in (4000, 6000):
                disk = r.disk_check(["g"] * union, free=278_429_159_424, roots=roots, measure=False)
                self.assertFalse(disk["ok"], disk)
                self.assertGreater(disk["kset_bytes"], 0.8 * 278_429_159_424 * 2)

    def test_refused_plans_lists_an_unavailable_plan_with_its_reason(self):
        """REVIEW-INPUT4 L2: a plan folder with one refused fold plan and one refused cohort plan."""
        ds = "salmonella_enterica__meropenem"
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            for arm, name, body in (
                    ("default", "lineage-0.json", {"plan": {"unavailable_reason": "the gbt holdout lacks a class"}}),
                    ("default", "random-0.json", {"plan": {"unavailable_reason": None}}),
                    ("cap-1000", "lineage-0.cohort.json", {"unavailable_reason": "outer fold 2 lacks a class"}),
                    ("cap-1000", "random-0.cohort.json", {"unavailable_reason": None})):
                d = work / "arms" / arm / ds / "plans"
                d.mkdir(parents=True, exist_ok=True)
                (d / name).write_text(json.dumps(body))
            (work / "arms" / "dropped-arm" / ds / "plans").mkdir(parents=True)
            (work / "arms" / "dropped-arm" / ds / "plans" / "lineage-0.json").write_text(
                json.dumps({"plan": {"unavailable_reason": "an arm no longer in the spec"}}))
            self.assertEqual(run.refused_plans(work, [ds], ["default", "cap-1000"]), [
                {"dataset": ds, "arm": "cap-1000", "plan": "lineage-0.cohort.json",
                 "reason": "outer fold 2 lacks a class"},
                {"dataset": ds, "arm": "default", "plan": "lineage-0.json",
                 "reason": "the gbt holdout lacks a class"}])

    def null_records(self, arm, lineage):
        """J records of one arm in one seed: lineage contrasts with the given decisions (else inconclusive), and
        random-split contrasts in family D."""
        out = []
        for learner in run.LEARNERS:
            for split in ("lineage", "random"):
                family = "C" if split == "lineage" else "D"
                out.append({"id": f"ds/{family}:{learner}:{split}:default->{arm}:log_loss", "dataset": "ds",
                            "treatment": arm, "learner": learner, "split": split, "family": family,
                            "decision": lineage.get(learner, "inconclusive") if split == "lineage" else "inconclusive"})
        return out

    def seed_folder(self, top, s):
        """An empty seed folder with the cohort.json null_study.one_seed writes."""
        d = Path(top) / f"seed-{s}"
        d.mkdir(parents=True)
        (d / "cohort.json").write_text(json.dumps({"seed": s, "genomes": 3000, "resistant": 600}))
        return d

    def null_seeds(self, top, arms, contrasts, unreliable=(), seeds=20, missing=None, start=0):
        """Seed folders start..start+seeds-1, each summary with its own `seconds`, as every real run's has."""
        for s in range(start, start + seeds):
            d = self.seed_folder(top, s)
            summary = {"arms": run.summarize(SPEC, arms, contrasts, missing or {}, list(unreliable)),
                       "unreliable": list(unreliable), "seconds": s + 0.5}
            (d / "summary.json").write_text(json.dumps(summary))

    def null_unproduced(self, arm):
        """One arm's records with every lineage contrast not produced: the random-split rows, and the lineage ids
        with a reason, as run.contrasts() gives them to summarize()."""
        rows = self.null_records(arm, {})
        missing = {x["id"]: "the gbt holdout lacks a class" for x in rows if x["split"] == "lineage"}
        return [x for x in rows if x["split"] == "random"], missing

    def calibration_arms(self):
        """calibration.json's arms as summarize() takes them."""
        return [{"name": a["name"], "factor": "calibration", "pairing": a["pairing"]}
                for a in run.load_spec(HERE / "calibration.json")[2] if not a.get("sentinel")]

    def test_the_calibration_judges_every_arm_its_spec_names(self):
        """REVIEW-INPUT17 M1: the gate judged only the arms it found, so an empty folder printed "seeds 0; pass True"
        and a calibration arm absent from every seed passed. Now cohort-path and seed-twin are judged by name; no
        seeds or a missing arm fails and is named, a mistyped folder is refused, and the null study names a missing
        arm instead of raising KeyError."""
        self.assertEqual(null_study.CALIBRATION_ARMS, ("cohort-path", "seed-twin"))
        self.assertEqual([a["name"] for a in self.calibration_arms()], list(null_study.CALIBRATION_ARMS))
        with tempfile.TemporaryDirectory() as tmp:
            r = null_study.report([tmp], calibration=True)
            self.assertEqual((r["seeds"], r["missing_arms"], r["pass"]), (0, ["cohort-path", "seed-twin"], False))
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                null_study.main(["report", tmp, "--calibration"])
            self.assertIn(f"calibration: seeds 0; pass False (no seed summaries in {tmp}; arm cohort-path is in no "
                          "seed; arm seed-twin is in no seed)", out.getvalue())
            # A mistyped folder is refused, in both studies.
            for flags in ([], ["--calibration"]):
                with self.assertRaises(SystemExit) as e:
                    null_study.main(["report", str(Path(tmp) / "seeds-typo"), *flags])
                self.assertEqual(str(e.exception), f"{Path(tmp) / 'seeds-typo'}: no such folder")
            # Empty, the null study names the arms it can't find.
            with self.assertRaises(SystemExit) as e:
                null_study.report([tmp])
            self.assertIn("null study: arm cap-1000, cap-4000, grr-v1, k-h, k-p, k-c-h-p is in no seed's summary.json "
                          "(0 summaries)", str(e.exception))
        # seed-twin absent from all 20 seeds; cohort-path clean in each.
        with tempfile.TemporaryDirectory() as tmp:
            self.null_seeds(tmp, self.calibration_arms()[:1], self.null_records("cohort-path", {}))
            r = null_study.report([tmp], calibration=True)
            self.assertEqual((r["arms"]["cohort-path"]["seeds"], r["missing_arms"], r["pass"]),
                             (20, ["seed-twin"], False))
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                null_study.main(["report", tmp, "--calibration"])
            self.assertIn("calibration: seeds 20; pass False (arm seed-twin is in no seed)", out.getvalue())
        # The null study with grr-v1 in no seed.
        arms = [{"name": a, "factor": "blocks", "pairing": "same_plan"} for a in null_study.PLANTED_ARMS]
        arms += [{"name": a, "factor": "cohort", "pairing": "common_outer_test_set_balanced"}
                 for a in null_study.COHORT_ARMS[:-1]]
        con = [r for a in null_study.PLANTED_ARMS for r in self.null_records(a, {l: "improves" for l in run.LEARNERS})]
        con += [r for a in null_study.COHORT_ARMS[:-1] for r in self.null_records(a, {})]
        with tempfile.TemporaryDirectory() as tmp:
            self.null_seeds(tmp, arms, con)
            with self.assertRaises(SystemExit) as e:
                null_study.report([tmp])
        self.assertIn("null study: arm grr-v1 is in no seed's summary.json (20 summaries)", str(e.exception))

    def test_every_judged_arm_needs_20_seeds_not_any(self):
        """REVIEW-INPUT17 L1 (mutants N3 and N4): with two calibration arms, one in 20 seeds and one in 19, the
        calibration fails; in the null study, a planted arm in 19 seeds fails it as a cohort arm does."""
        both = self.calibration_arms()
        random_only, missing = self.null_unproduced("cohort-path")
        twin = self.null_records("seed-twin", {})
        with tempfile.TemporaryDirectory() as tmp:
            self.null_seeds(tmp, both, self.null_records("cohort-path", {}) + twin, seeds=19)
            self.null_seeds(tmp, both, random_only + twin, seeds=1, missing=missing, start=19)
            r = null_study.report([tmp], calibration=True)
        self.assertEqual((r["arms"]["seed-twin"]["seeds"], r["arms"]["cohort-path"]["seeds"]), (20, 19))
        self.assertEqual((r["missing_arms"], r["thin_arms"], r["pass"]), ([], ["cohort-path"], False))
        arms = [{"name": a, "factor": "blocks", "pairing": "same_plan"} for a in null_study.PLANTED_ARMS]
        arms += [{"name": a, "factor": "cohort", "pairing": "common_outer_test_set_balanced"}
                 for a in null_study.COHORT_ARMS]
        cohort = [r for a in null_study.COHORT_ARMS for r in self.null_records(a, {})]
        planted = {a: self.null_records(a, {l: "improves" for l in run.LEARNERS}) for a in null_study.PLANTED_ARMS}
        short, gone = self.null_unproduced("k-p")
        with tempfile.TemporaryDirectory() as tmp:
            self.null_seeds(tmp, arms, cohort + [r for rows in planted.values() for r in rows], seeds=19)
            self.null_seeds(tmp, arms, cohort + planted["k-h"] + planted["k-c-h-p"] + short, seeds=1, missing=gone,
                            start=19)
            r = null_study.report([tmp])
        self.assertEqual((r["arms"]["k-p"]["seeds"], r["arms"]["k-h"]["seeds"], r["arms"]["grr-v1"]["seeds"]),
                         (19, 20, 20))
        self.assertEqual(r["planted_improves_share"], {a: 1.0 for a in null_study.PLANTED_ARMS})
        self.assertFalse(r["pass"])

    def test_the_calibration_judges_the_runner_s_sentinel_too(self):
        """REVIEW-INPUT17B L1 (mutant N3e): the calibration judges every arm its summaries hold, the runner's sentinel
        included, but no test held a sentinel. With both named arms in 20 seeds and the sentinel in 19, it fails and
        names the sentinel; with the sentinel in 20 too, it passes."""
        arms = run.load_spec(HERE / "calibration.json")[2]
        self.assertEqual([a["name"] for a in arms], [*null_study.CALIBRATION_ARMS, run.SENTINEL["name"]])
        named = [r for a in null_study.CALIBRATION_ARMS for r in self.null_records(a, {})]
        random_only, missing = self.null_unproduced(run.SENTINEL["name"])
        for produced, passes in ((19, False), (20, True)):
            with tempfile.TemporaryDirectory() as tmp:
                self.null_seeds(tmp, arms, named + self.null_records(run.SENTINEL["name"], {}), seeds=produced)
                self.null_seeds(tmp, arms, named + random_only, seeds=20 - produced, missing=missing, start=produced)
                r = null_study.report([tmp], calibration=True)
                out = io.StringIO()
                with contextlib.redirect_stdout(out):
                    null_study.main(["report", tmp, "--calibration"])
            self.assertEqual([r["arms"][a]["seeds"] for a in null_study.CALIBRATION_ARMS], [20, 20])
            self.assertEqual((r["arms"][run.SENTINEL["name"]]["seeds"], r["missing_arms"], r["pass"]),
                             (produced, [], passes))
            self.assertEqual(r["thin_arms"], [] if passes else [run.SENTINEL["name"]])
            self.assertIn("calibration: seeds 20; pass True" if passes else
                          "calibration: seeds 20; pass False (arm sentinel-seed-twin has a lineage contrast in 19 "
                          "seeds, under 20)", out.getvalue())
        # REVIEW-INPUT18 L3 (mutant V2): the sentinel's rates are judged as the named arms' are. In all 20 seeds, a
        # sentinel that improves, or worsens, on every lineage contrast fails the calibration.
        for decision in ("improves", "worse"):
            with tempfile.TemporaryDirectory() as tmp:
                self.null_seeds(tmp, arms, named + self.null_records(run.SENTINEL["name"],
                                                                     {l: decision for l in run.LEARNERS}))
                r = null_study.report([tmp], calibration=True)
            self.assertEqual((r["missing_arms"], r["thin_arms"], r["pass"]), ([], [], False), decision)
            self.assertEqual(r["arms"][run.SENTINEL["name"]][decision], 100)
        # REVIEW-INPUT18 L3 (mutant V6): every arm found is checked for 20 seeds, the first included: cohort-path in
        # 19 seeds fails beside a full seed-twin and sentinel.
        short, gone = self.null_unproduced("cohort-path")
        rest = [r for a in ("seed-twin", run.SENTINEL["name"]) for r in self.null_records(a, {})]
        with tempfile.TemporaryDirectory() as tmp:
            self.null_seeds(tmp, arms, self.null_records("cohort-path", {}) + rest, seeds=19)
            self.null_seeds(tmp, arms, short + rest, seeds=1, missing=gone, start=19)
            r = null_study.report([tmp], calibration=True)
        self.assertEqual(list(r["arms"])[0], "cohort-path")
        self.assertEqual((r["thin_arms"], r["pass"]), (["cohort-path"], False))

    def test_a_seed_is_the_one_its_cohort_json_records(self):
        """REVIEW-INPUT17B L2 (its cases B1-B3): seeds were keyed by folder name, so 10 seeds passed as 20 beside
        `seed-0101` copies, `.bak` copies or copies renamed to seeds 111-120. Now a seed is cohort.json's: its copies
        count once when their summaries are byte-identical and are refused, naming both, when they differ; a summary
        whose seed can't be read is refused."""
        both = self.calibration_arms()
        con = self.null_records("cohort-path", {}) + self.null_records("seed-twin", {})
        with tempfile.TemporaryDirectory() as tmp:
            a = Path(tmp) / "a"
            self.null_seeds(a, both, con, seeds=10, start=101)
            # B1 and B2: seed-0101 and seed-101.bak beside seed-101.
            for copy in ("seed-0101", "seed-101.bak"):
                shutil.copytree(a / "seed-101", a / copy)
            # B3: seeds 101-110 copied under the names of seeds 111-120.
            renamed = Path(tmp) / "renamed"
            for s in range(101, 111):
                shutil.copytree(a / f"seed-{s}", renamed / f"seed-{s + 10}")
            for dirs in ([a], [a, renamed], [renamed]):
                r = null_study.report(dirs, calibration=True)
                self.assertEqual((r["seeds"], r["arms"]["seed-twin"]["seeds"], r["thin_arms"], r["pass"]),
                                 (10, 10, ["cohort-path", "seed-twin"], False), dirs)
            # The control: ten more seeds of their own make 20, and pass.
            self.null_seeds(a, both, con, seeds=10, start=111)
            r = null_study.report([a, renamed], calibration=True)
            self.assertEqual((r["seeds"], r["arms"]["seed-twin"]["seeds"], r["pass"]), (20, 20, True))
            # A copy under another name with other content is refused, naming both files.
            (a / "seed-0101" / "summary.json").write_text(json.dumps({"arms": {}}))
            with self.assertRaises(SystemExit) as e:
                null_study.report([a], calibration=True)
            self.assertEqual(str(e.exception), f"seed 101 has two different summaries: {a / 'seed-0101' / 'summary.json'}"
                                               f" and {a / 'seed-101' / 'summary.json'}")
            # REVIEW-INPUT18 L3 (mutant W4): a copy of the same length is compared byte for byte, as two runs of one
            # seed that differ only in `seconds` would be.
            first = (a / "seed-101" / "summary.json").read_text()
            same_length = first.replace('"seconds": 101.5', '"seconds": 102.5')
            self.assertEqual((len(same_length), same_length != first), (len(first), True))
            (a / "seed-0101" / "summary.json").write_text(same_length)
            with self.assertRaises(SystemExit) as e:
                null_study.report([a], calibration=True)
            self.assertIn("seed 101 has two different summaries", str(e.exception))
            shutil.rmtree(a / "seed-0101")
            # REVIEW-INPUT18 F1: a copy whose cohort.json names another seed is refused, naming both files, whether the
            # seed is new to the study (the reviewer's F1 recorded copies of 101-110 as 111-120) or negative (-5, their
            # N1).
            for s, folder in ((130, "seed-130"), (-5, "seed-x")):
                forged = Path(tmp) / f"forged{s}"
                shutil.copytree(a / "seed-101", forged / folder)
                (forged / folder / "cohort.json").write_text(json.dumps({"seed": s, "genomes": 3000, "resistant": 600}))
                with self.assertRaises(SystemExit) as e:
                    null_study.report([a, forged], calibration=True)
                self.assertEqual(str(e.exception), f"seeds 101 and {s} have byte-identical summaries: "
                                                   f"{a / 'seed-101' / 'summary.json'} and {forged / folder / 'summary.json'}")
            # REVIEW-INPUT18 L3 (mutant W7): seeds 100 apart are two seeds; 101-110 and 201-210 make 20, and pass.
            hundred = Path(tmp) / "hundred"
            self.null_seeds(hundred, both, con, seeds=10, start=101)
            self.null_seeds(hundred, both, con, seeds=10, start=201)
            r = null_study.report([hundred], calibration=True)
            self.assertEqual((r["seeds"], r["arms"]["seed-twin"]["seeds"], r["pass"]), (20, 20, True))
            # A summary whose seed can't be read is refused, in both studies.
            cohort = a / "seed-105" / "cohort.json"
            for text, why in (("", "JSONDecodeError"), ("[105]", "TypeError"), ('{"genomes": 3000}', "KeyError"),
                              ('{"seed": "105"}', "'105', not an integer"), ('{"seed": true}', "True, not an integer"),
                              ('{"seed": 105.0}', "105.0, not an integer"), (None, "FileNotFoundError")):
                cohort.unlink(missing_ok=True)
                if text is not None:
                    cohort.write_text(text)
                for flags in ([], ["--calibration"]):
                    with self.assertRaises(SystemExit) as e:
                        null_study.main(["report", str(a), *flags])
                    self.assertTrue(str(e.exception).startswith(f"{a / 'seed-105' / 'summary.json'}: "), e.exception)
                    self.assertIn(why, str(e.exception))

    def test_a_seed_killed_between_its_writes_reruns(self):
        """REVIEW-INPUT18 L2: one_seed wrote summary.json before cohort.json, and a rerun skips a folder with a
        summary, so a seed killed between the two was never rerun and report refused it. Killed after each of its
        writes in turn, a seed with a summary always has its cohort.json, and a rerun leaves one seed report reads.
        synthetic_root.py and run.py are stood in for, writing a truth.json, a summary and contrasts."""
        summary = {"arms": run.summarize(SPEC, self.calibration_arms(), self.null_records("cohort-path", {})
                                         + self.null_records("seed-twin", {}), {}, []), "unreliable": [], "seconds": 1.5}

        def stand_in(cmd, **_):
            if Path(cmd[2]).name == "synthetic_root.py":
                Path(cmd[3]).mkdir(parents=True)
                (Path(cmd[3]) / "truth.json").write_bytes(b'{"genomes": [{"resistant": true}, {"resistant": false}]}')
            else:
                work = Path(cmd[cmd.index("--work") + 1])
                (work / "j").mkdir(parents=True)
                (work / "summary.json").write_bytes(json.dumps(summary).encode())
                (work / "j" / "contrasts.json").write_bytes(b"[]")
            return subprocess.CompletedProcess(cmd, 0)

        class Killed(BaseException):
            pass

        writes = []

        def killing(real, at, target):
            """`real`, which writes the path in its argument `target`, raising Killed after the at-th write."""
            def write(*args, **kwargs):
                out = real(*args, **kwargs)
                writes.append(str(args[target]))
                if len(writes) == at:
                    raise Killed
                return out
            return write

        def one_seed(out, at=0):
            writes.clear()
            with mock.patch.object(null_study.subprocess, "run", stand_in), \
                    mock.patch.object(null_study.shutil, "copy", killing(shutil.copy, at, 1)), \
                    mock.patch.object(null_study.os, "replace", killing(os.replace, at, 1)), \
                    mock.patch.object(Path, "write_text", killing(Path.write_text, at, 0)):
                return null_study.one_seed(7, out, "grm", 1)

        with tempfile.TemporaryDirectory() as tmp:
            one_seed(Path(tmp) / "clean")
            done = len(writes)
            self.assertGreaterEqual(done, 4, writes)
            self.assertEqual(writes[-1], str(Path(tmp) / "clean" / "seed-7" / "summary.json"), writes)
            for at in range(1, done + 1):
                out = Path(tmp) / f"killed-{at}"
                with self.assertRaises(Killed):
                    one_seed(out, at)
                seed = out / "seed-7"
                if (seed / "summary.json").exists():
                    self.assertEqual(json.loads((seed / "cohort.json").read_text())["seed"], 7, at)
                one_seed(out)
                r = null_study.report([out], calibration=True)
                self.assertEqual((r["seeds"], r["arms"]["seed-twin"]["seeds"]), (1, 1), at)
                self.assertEqual((seed / "summary.json").read_bytes(), json.dumps(summary).encode(), at)

    def test_a_seed_counts_once_however_many_folders_hold_it(self):
        """REVIEW-INPUT17 L2: the same 10-seed folder given twice, or two folders holding the same seeds, counted 20
        seeds and passed. Byte-identical summaries of one seed now count once; different ones are refused."""
        both = self.calibration_arms()
        con = self.null_records("cohort-path", {}) + self.null_records("seed-twin", {})
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp) / "a", Path(tmp) / "b"
            self.null_seeds(a, both, con, seeds=10)
            self.null_seeds(b, both, con, seeds=10)
            for dirs in ([a, a], [a, b]):
                r = null_study.report(dirs, calibration=True)
                self.assertEqual((r["seeds"], r["arms"]["seed-twin"]["seeds"], r["thin_arms"], r["pass"]),
                                 (10, 10, ["cohort-path", "seed-twin"], False), dirs)
            # Seeds 10-19 in a third folder make 20 distinct seeds.
            c = Path(tmp) / "c"
            self.null_seeds(c, both, con, seeds=10, start=10)
            r = null_study.report([a, b, c], calibration=True)
            self.assertEqual((r["seeds"], r["arms"]["seed-twin"]["seeds"], r["pass"]), (20, 20, True))
            # The same seed with other content is refused, naming both files.
            (b / "seed-3" / "summary.json").write_text(json.dumps({"arms": {}}))
            with self.assertRaises(SystemExit) as e:
                null_study.report([a, b], calibration=True)
            self.assertEqual(str(e.exception), f"seed 3 has two different summaries: {a / 'seed-3' / 'summary.json'} "
                                               f"and {b / 'seed-3' / 'summary.json'}")

    def test_null_study_passes_only_with_a_lineage_contrast_in_20_seeds_per_arm(self):
        """REVIEW-NULL16 L1 (its case R4): a calibration arm with lineage contrasts in 1 of 20 seeds, 95 not produced,
        passed on "improves 0/5" because the gate counted summary files. Now each arm needs 20 seeds that gave it a
        lineage contrast, and the not-produced count is printed."""
        twin = [{"name": "seed-twin", "factor": "calibration", "pairing": "same_plan"}]
        random_only, missing = self.null_unproduced("seed-twin")
        with tempfile.TemporaryDirectory() as tmp:
            self.null_seeds(tmp, twin, self.null_records("seed-twin", {}), seeds=1)
            self.null_seeds(tmp, twin, random_only, seeds=19, missing=missing, start=1)
            r = null_study.report([tmp], calibration=True)
            a = r["arms"]["seed-twin"]
            self.assertEqual((r["seeds"], a["seeds"], a["contrasts"], a["not_produced"], a["improves"]),
                             (20, 1, 5, 95, 0))
            self.assertFalse(r["pass"])
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                null_study.main(["report", tmp, "--calibration"])
            self.assertIn("not produced 95; seeds with a lineage contrast 1/20", out.getvalue())
            self.assertIn("lineage contrasts not produced: 95", out.getvalue())
            self.assertIn("calibration: seeds 20; pass False", out.getvalue())
        # 25 summaries of both calibration arms: 20 seeds with seed-twin's lineage contrasts pass, 19 don't.
        both = self.calibration_arms()
        path = [r for r in self.null_records("cohort-path", {})]
        for produced, passes in ((20, True), (19, False)):
            with tempfile.TemporaryDirectory() as tmp:
                self.null_seeds(tmp, both, path + self.null_records("seed-twin", {}), seeds=produced)
                self.null_seeds(tmp, both, path + random_only, seeds=25 - produced, missing=missing, start=produced)
                r = null_study.report([tmp], calibration=True)
            self.assertEqual((r["seeds"], r["arms"]["seed-twin"]["seeds"], r["pass"]), (25, produced, passes))
            self.assertEqual(r["thin_arms"], [] if passes else ["seed-twin"])
        # The null study: one cohort arm short of 20 such seeds fails it, all else passing.
        arms = [{"name": a, "factor": "blocks", "pairing": "same_plan"} for a in null_study.PLANTED_ARMS]
        arms += [{"name": a, "factor": "cohort", "pairing": "common_outer_test_set_balanced"}
                 for a in null_study.COHORT_ARMS]
        planted = [r for a in null_study.PLANTED_ARMS for r in self.null_records(a, {l: "improves" for l in run.LEARNERS})]
        full = planted + [r for a in null_study.COHORT_ARMS for r in self.null_records(a, {})]
        short, gone = self.null_unproduced("grr-v1")
        thin = planted + [r for a in null_study.COHORT_ARMS[:-1] for r in self.null_records(a, {})] + short
        with tempfile.TemporaryDirectory() as tmp:
            self.null_seeds(tmp, arms, full, seeds=19)
            self.null_seeds(tmp, arms, thin, seeds=1, missing=gone, start=19)
            r = null_study.report([tmp])
        self.assertEqual((r["arms"]["grr-v1"]["seeds"], r["arms"]["grr-v1"]["not_produced"]), (19, 5))
        self.assertEqual(r["arms"]["cap-1000"]["seeds"], 20)
        self.assertFalse(r["pass"])

    def test_null_study_rates_count_decision_bearing_lineage_contrasts_only(self):
        """PR #69's review (null_study.py:144): today's summaries mark random-split rows descriptive; they must not
        enter a rate, so a planted arm improving on every lineage contrast scores 1.0 and the study can pass, and a
        calibration arm's lineage false-positive rate isn't halved."""
        arms = [{"name": a, "factor": "blocks", "pairing": "same_plan"} for a in null_study.PLANTED_ARMS]
        arms += [{"name": a, "factor": "cohort", "pairing": "common_outer_test_set_balanced"}
                 for a in null_study.COHORT_ARMS]
        contrasts = [r for a in null_study.PLANTED_ARMS for r in self.null_records(a, {l: "improves" for l in run.LEARNERS})]
        contrasts += [r for a in null_study.COHORT_ARMS for r in self.null_records(a, {})]
        with tempfile.TemporaryDirectory() as tmp:
            self.null_seeds(tmp, arms, contrasts)
            r = null_study.report([tmp])
        self.assertEqual(r["planted_improves_share"], {a: 1.0 for a in null_study.PLANTED_ARMS})
        self.assertTrue(r["pass"])
        self.assertEqual(r["arms"]["k-h"]["contrasts"], 5 * 20)
        with tempfile.TemporaryDirectory() as tmp:
            self.null_seeds(tmp, [{"name": "seed-twin", "factor": "calibration", "pairing": "same_plan"}],
                            self.null_records("seed-twin", {"gbt": "improves"}))
            r = null_study.report([tmp], calibration=True)
        self.assertAlmostEqual(r["arms"]["seed-twin"]["per_seed_improves"]["mean"], 0.2)
        self.assertFalse(r["pass"])

    def test_null_study_calibration_sees_what_the_sentinel_masks(self):
        """PR #69's review (null_study.py:143, P1): a sentinel firing on every seed turns its learners' rows
        unreliable; calibration counts J's decisions, fails, and reports the firings."""
        sentinel = [{"name": "sentinel-seed-twin", "factor": "calibration", "pairing": "same_plan", "sentinel": True}]
        con = self.null_records("sentinel-seed-twin", {"gbt": "improves", "rf": "worse"})
        unreliable = [{"species": "Salmonella enterica", "learner": l, "sentinel": "sentinel-seed-twin", "reason": "fired",
                       "contrasts": [x["id"] for x in con if x["learner"] == l]} for l in ("gbt", "rf")]
        with tempfile.TemporaryDirectory() as tmp:
            self.null_seeds(tmp, sentinel, con, unreliable)
            masked = json.loads((Path(tmp) / "seed-0" / "summary.json").read_text())["arms"]["sentinel-seed-twin"]
            self.assertEqual(masked["decisions"]["ds:gbt:lineage"], "unreliable")
            r = null_study.report([tmp], calibration=True)
        a = r["arms"]["sentinel-seed-twin"]
        self.assertEqual((a["improves"], a["worse"], a["contrasts"]), (20, 20, 100))
        self.assertFalse(r["pass"])
        self.assertEqual(r["sentinel_firings"], {"seeds_with_a_firing": 20, "firings": 40, "learners": {"gbt": 20, "rf": 20}})

    def test_null_study_reads_summaries_written_before_the_mask(self):
        """REVIEW-NULL16 L2 (mutant B2): the stored seeds' summaries have no raw_decisions, and both splits in
        decisions as J decided them; the study takes the lineage rows from decisions."""
        decisions = {f"ds:{l}:{split}": "improves" if l == "gbt" or split == "random" else "inconclusive"
                     for l in run.LEARNERS for split in ("lineage", "random")}
        with tempfile.TemporaryDirectory() as tmp:
            for s in range(20):
                (self.seed_folder(tmp, s) / "summary.json").write_text(json.dumps({"arms": {"seed-twin": {
                    "factor": "calibration", "pairing": "same_plan", "decisions": decisions}}, "seconds": s}))
            r = null_study.report([tmp], calibration=True)
        a = r["arms"]["seed-twin"]
        self.assertEqual((a["improves"], a["contrasts"], a["seeds"]), (20, 100, 20))
        self.assertAlmostEqual(a["per_seed_improves"]["mean"], 0.2)
        self.assertEqual(r["sentinel_firings"], {"seeds_with_a_firing": 0, "firings": 0, "learners": {}})

    def test_null_study_reads_the_arm_names_from_before_the_rename(self):
        """REVIEW-NULL16 I3: the stored null seeds name the block arms K+C, K+H, K+P and K+C+H+P (before round 12);
        `report` reads them under today's names instead of raising KeyError: 'k-h'."""
        def decisions(value):
            return {f"ds:{l}:{split}": value for l in run.LEARNERS for split in ("lineage", "random")}

        old = {"K+C": "inconclusive", "K+H": "improves", "K+P": "improves", "K+C+H+P": "improves",
               **{a: "inconclusive" for a in null_study.COHORT_ARMS}}
        with tempfile.TemporaryDirectory() as tmp:
            for s in range(20):
                (self.seed_folder(tmp, s) / "summary.json").write_text(json.dumps({"arms": {
                    arm: {"factor": "blocks", "pairing": "same_plan", "decisions": decisions(v)}
                    for arm, v in old.items()}, "seconds": s}))
            r = null_study.report([tmp])
        self.assertEqual(sorted(r["arms"]), sorted(["k-c", *null_study.PLANTED_ARMS, *null_study.COHORT_ARMS]))
        self.assertEqual(r["planted_improves_share"], {a: 1.0 for a in null_study.PLANTED_ARMS})
        self.assertEqual((r["arms"]["k-c"]["contrasts"], r["arms"]["k-c"]["improves"]), (100, 0))
        self.assertTrue(r["pass"])

    def test_null_study_counts_a_contrast_not_produced_apart_and_refuses_an_arm_without_any(self):
        """REVIEW-NULL16 L2 (mutants C1-C3): a lineage contrast not produced is counted apart, never as a decision in
        a denominator; an arm with no lineage contrast in any seed is refused."""
        twin = [{"name": "seed-twin", "factor": "calibration", "pairing": "same_plan"}]
        rows = self.null_records("seed-twin", {"gbt": "improves"})
        gap = {x["id"]: "outer fold 2 lacks a class" for x in rows if x["split"] == "lineage" and x["learner"] == "rf"}
        with tempfile.TemporaryDirectory() as tmp:
            self.null_seeds(tmp, twin, [x for x in rows if x["id"] not in gap], missing=gap)
            shown = json.loads((Path(tmp) / "seed-0" / "summary.json").read_text())["arms"]["seed-twin"]["decisions"]
            self.assertTrue(shown["ds:rf:lineage"].startswith("not produced"))
            r = null_study.report([tmp], calibration=True)
        a = r["arms"]["seed-twin"]
        self.assertEqual((a["contrasts"], a["improves"], a["not_produced"], a["seeds"]), (80, 20, 20, 20))
        self.assertAlmostEqual(a["per_seed_improves"]["mean"], 0.25)
        random_only, missing = self.null_unproduced("seed-twin")
        with tempfile.TemporaryDirectory() as tmp:
            self.null_seeds(tmp, twin, random_only, missing=missing)
            with self.assertRaises(SystemExit) as e:
                null_study.report([tmp], calibration=True)
        self.assertEqual(str(e.exception), "arm seed-twin: no lineage contrast was produced in any seed")

    def test_a_spec_with_repeats_and_a_cohort_arm_is_refused(self):
        """PR #69's review (stats.py:300): with repeats > 1 a cohort arm's subsample can differ by repeat, which J
        can't pair; no decided spec uses repeats > 1, so such a spec is refused up front."""
        self.assertEqual((SPEC["repeats"], json.loads((HERE / "calibration.json").read_text())["repeats"]), (1, 1))
        spec = copy.deepcopy(SPEC)
        spec["repeats"] = 2
        with self.assertRaises(SystemExit) as e:
            run.load_spec(self.write(spec))
        # REVIEW-NULL16 I4: the refusal names the arms.
        self.assertIn("repeats 2 with a cohort-paired arm (cap-1000, cap-4000, grr-v1)", str(e.exception))
        spec["arms"] = [{"name": "k-h", "blocks": ["K", "H"]}, {"name": "cap-4000", "cap": 4000}]
        with self.assertRaises(SystemExit) as e:
            run.load_spec(self.write(spec))
        self.assertIn("repeats 2 with a cohort-paired arm (cap-4000):", str(e.exception))
        spec["arms"] = [{"name": "k-h", "blocks": ["K", "H"]}]
        run.load_spec(self.write(spec))

    def test_summary_keeps_unreliable_and_descriptive_rows_out_of_its_decisions(self):
        """REVIEW-INPUT M14: summary.json's override of flagged contrasts."""
        arm = {"name": "k-h", "factor": "blocks", "pairing": "same_plan"}

        def x(learner, split, decision, family="C"):
            return {"id": f"ds/{family}:{learner}:{split}:default->k-h:log_loss", "dataset": "ds",
                    "treatment": "k-h", "learner": learner, "split": split, "family": family, "decision": decision}

        rows = [x("gbt", "lineage", "improves"), x("lasso", "lineage", "worse"),
                x("gbt", "random", "inconclusive", "D")]
        out = run.summarize(SPEC, [arm], rows, {}, [{"contrasts": [rows[0]["id"]]}])["k-h"]
        self.assertEqual(out["decisions"], {"ds:gbt:lineage": "unreliable", "ds:lasso:lineage": "worse",
                                            "ds:gbt:random": "descriptive"})
        self.assertEqual(out["raw_decisions"]["ds:gbt:lineage"], "improves")
        self.assertEqual((out["any_improves"], out["any_worse"]), (False, True))
        quiet = run.summarize(SPEC, [arm], rows, {}, [])["k-h"]
        self.assertTrue(quiet["any_improves"])

    def test_the_calibration_arms_change_no_input(self):
        _, default, arms = run.load_spec(HERE / "calibration.json")
        self.assertEqual([(a["name"], a["factor"], a["pairing"]) for a in arms],
                         [("cohort-path", "calibration", "common_outer_test_set_balanced"),
                          ("seed-twin", "calibration", "same_plan"),
                          ("sentinel-seed-twin", "calibration", "same_plan")])
        for a in arms:
            self.assertEqual(run.root_key(a), run.root_key(default))
            self.assertEqual((a["blocks"], a["screened_features"], a["k"]), (["K"], 4000, 31))
        for bad in ({"name": "x", "calibration": "cohort_path"}, {"name": "x", "calibration": "other", "seed": 1},
                    {"name": "x", "calibration": "learner_seed", "learner_seed": 1, "cap": 1000}):
            spec = copy.deepcopy(SPEC)
            spec["arms"] = [bad]
            with self.assertRaises(SystemExit, msg=bad):
                run.load_spec(self.write(spec))


@unittest.skipUnless(os.environ.get("GRM_INPUT_E2E_BINARY"), "set GRM_INPUT_E2E_BINARY to a grm release binary")
class EndToEndTests(unittest.TestCase):
    """The synthetic comparison: a planted effect (H and P call a carbapenemase family K
    cannot see) that the k-h, k-p and k-c-h-p arms must find, and a null C block that
    must not improve anything."""

    @classmethod
    def setUpClass(cls):
        out = os.environ.get("GRM_INPUT_E2E_OUT")
        cls.tmp = None if out else tempfile.TemporaryDirectory()
        cls.out = Path(out or cls.tmp.name)
        raw, work = cls.out / "raw", cls.out / "work"
        subprocess.run([sys.executable, "-B", str(HERE / "synthetic_root.py"), str(raw)], check=True)
        commit = os.environ.get("GRM_INPUT_E2E_COMMIT", "0" * 40)
        # The spec's 27.5 MiB per kset is a real Salmonella genome's; at that size this comparison (a 3,000-genome
        # pool and 12,000 genome-ksets in its roots) needs about 419 GiB and refuses on an xlarge before any count.
        # Synthetic genomes are 23 kb, their ksets about 0.17 MB: the e2e states 1 MiB, an upper bound.
        spec = json.loads((HERE / "spec.json").read_text())
        spec["pool_budget"]["default_kset_bytes"] = 2 ** 20
        # The gate on the synthetic default's own values, pinned here: a curate whose output or record changes
        # fails this run before H-1 does (REVIEW-INPUT8 M1).
        spec["revision12_gate"] = dict(SYNTHETIC_GATE)
        spec["synthetic"] = True
        cls.out.mkdir(parents=True, exist_ok=True)
        (cls.out / "spec-e2e.json").write_text(json.dumps(spec, indent=1))
        # A bare --grm, found on PATH, with a decoy named grm in the run's cwd (REVIEW-INPUT8 L2).
        binary = Path(os.environ["GRM_INPUT_E2E_BINARY"]).resolve()
        (cls.out / "grm").write_text("not the binary\n")
        env = {**os.environ, "PATH": f"{binary.parent}{os.pathsep}{os.environ['PATH']}"}
        subprocess.run([sys.executable, "-B", str(HERE / "run.py"), str(cls.out / "spec-e2e.json"), "--raw", str(raw),
                        "--work", str(work), "--grm", binary.name,
                        "--source-commit", commit, "--jobs", os.environ.get("GRM_INPUT_E2E_JOBS", "3"),
                        "--trace-training", "--sigdrug-vocabulary", str(raw / "vocabulary")],
                       check=True, cwd=cls.out, env=env)
        cls.summary = json.loads((work / "summary.json").read_text())
        cls.report = json.loads((work / "report" / "ablation-report.json").read_text())
        cls.work = work

    @classmethod
    def tearDownClass(cls):
        if cls.tmp:
            cls.tmp.cleanup()

    def decisions(self, arm):
        return self.summary["arms"][arm]["decisions"]

    def test_every_planned_contrast_is_produced(self):
        counts = self.report["contrast_counts"]
        self.assertEqual(counts["planned"], 12 * 5 * 2)  # 11 arms and the sentinel
        self.assertEqual(counts["not_produced"], 0, self.report["not_produced"])
        self.assertTrue((self.work / "report" / "ablation-report.md").exists())

    def test_the_planted_effect_is_found(self):
        for arm in ("k-h", "k-p", "k-c-h-p"):
            d = self.decisions(arm)
            self.assertTrue(self.summary["arms"][arm]["any_improves"], (arm, d))
            self.assertEqual(d["salmonella_enterica__meropenem:lasso:lineage"], "improves", (arm, d))
            # The random split is descriptive: an estimate, no p, no decision.
            self.assertEqual(d["salmonella_enterica__meropenem:lasso:random"], "descriptive", (arm, d))

    def test_the_null_block_improves_nothing(self):
        self.assertFalse(self.summary["arms"]["k-c"]["any_improves"], self.decisions("k-c"))

    def test_each_contrast_names_its_arms_pairing_rule(self):
        want = {a["name"]: a["pairing"] for a in run.load_spec(HERE / "spec.json")[2]}
        for x in self.report["contrasts"]:
            self.assertEqual(x["pairing"], want[x["treatment"]], x["id"])

    def test_random_split_contrasts_are_descriptive_and_the_sentinel_is_checked(self):
        for x in self.report["contrasts"]:
            if x["split"] == "random":
                self.assertEqual((x["family"], x["p"], x["decision"]), ("D", None, "inconclusive"), x["id"])
                self.assertIsNotNone(x["estimate"], x["id"])
            else:
                self.assertEqual(x["family"], "C", x["id"])
        self.assertEqual(self.report["manifest"]["sentinel_arms"], ["sentinel-seed-twin"])
        flagged = {(u["learner"]) for u in self.report["unreliable"]}
        sentinel = self.summary["arms"]["sentinel-seed-twin"]["raw_decisions"]
        for learner in ("scm", "cart", "rf", "gbt", "lasso"):
            quiet = sentinel[f"salmonella_enterica__meropenem:{learner}:lineage"] == "inconclusive"
            self.assertEqual(learner in flagged, not quiet, learner)
        self.assertIn("### Random split: descriptive only",
                      (self.work / "report" / "ablation-report.md").read_text())

    def default_root(self):
        (root,) = [r for r in (self.work / "roots").iterdir()
                   if json.loads((r / "curated.json").read_text())["arm"] == "default"]
        return root

    def test_the_report_prints_each_cohort_arm_s_kind_rate_or_not_measured(self):
        kinds = {"cap-1000": "cap-down", "cap-4000": "cap-up", "grr-v1": "selection-swap"}
        expect = {}
        for arm, kind in kinds.items():
            rate = SPEC["base_rates"].get(kind)
            expect[arm] = (f"No-signal base rate for this arm type ({kind}): not measured." if rate is None else
                           f"No-signal base rate for this arm type ({kind}): {rate['rate'] * 100:.1f}%")
        seen = set()
        md = (self.work / "report" / "ablation-report.md").read_text()
        for line in md.splitlines():
            if line.startswith("- `") and "/C:" in line:
                arm = line.split("->", 1)[1].split(":", 1)[0]
                if arm in expect:
                    self.assertIn(expect[arm], line)
                    seen.add(arm)
                else:
                    self.assertNotIn("base rate", line, line[:120])
        self.assertEqual(seen, set(kinds))

    def test_summary_lists_every_refused_plan_with_its_reason(self):
        """REVIEW-INPUT3 I-a: a keyed plan with a fold or holdout side lacking a class is refused, not repaired."""
        listed = {(r["arm"], r["plan"]) for r in self.summary["refused_plans"]}
        ds = "salmonella_enterica__meropenem"
        for path in (self.work / "arms").glob(f"*/{ds}/plans/*.json"):
            p = json.loads(path.read_text())
            reason = p["plan"]["unavailable_reason"] if "plan" in p else p["unavailable_reason"]
            self.assertEqual((path.parent.parent.parent.name, path.name) in listed, reason is not None, path)

    def test_every_arm_splits_on_one_species_grouping(self):
        """REVIEW-INPUT3 N2 and REVIEW-INPUT4 M1: one structure.json over exactly the union of the arms' cohorts,
        built once; the default and the cohort arms split on it, the others take the default's plans."""
        groups = json.loads((self.work / "species-groups.json").read_text())
        self.assertEqual(self.summary["species_groups"], groups)
        pool = json.loads(Path(groups["structure"]).read_text())["genomes"]
        union = set()
        for root in (self.work / "roots").iterdir():
            if not root.name.startswith("pool-"):
                union |= set((root / "species" / "salmonella_enterica" / "genomes.txt").read_text().split())
        self.assertEqual(pool, sorted(union))
        self.assertEqual(groups["genomes"], len(pool))
        _, _, arms = run.load_spec(HERE / "spec.json")
        split_on = {"default"} | {a["name"] for a in arms if a["pairing"] == "common_outer_test_set_balanced"}
        for log in (self.work / "logs").glob("arm-*.log"):
            arm = next(a for a in sorted(split_on | {a["name"] for a in arms}, key=len, reverse=True)
                       if log.name.startswith(f"arm-{a}-"))
            self.assertEqual(f"--groups {groups['structure']}" in log.read_text(), arm in split_on, log.name)
        manifest = json.loads((self.work / "k" / "manifest.json").read_text())
        self.assertEqual(manifest["species_pool"], {"genomes": len(pool), "arms": groups["arms"],
                                                    "structure_sha256": groups["sha256"]})
        md = (self.work / "report" / "ablation-report.md").read_text()
        self.assertIn(f"- Species pool: the union of the cohorts of this run's {len(groups['arms'])} arms", md)
        # REVIEW-INPUT9 L1: the e2e's spec says synthetic, and K records it.
        self.assertIs(manifest["synthetic"], True)
        # REVIEW-INTEG7 L2: the spec carried the pinned gate too, and K records that it passed.
        self.assertEqual(manifest["revision12_gate"], "passed")
        self.assertIn("- Synthetic run: yes.", md)
        self.assertIn("adding an arm means rerunning the whole comparison", md)

    def test_the_disk_estimate_covers_every_root_the_run_counted(self):
        """REVIEW-INPUT5 M1(a): the on-box estimate names exactly the data roots the run counted, and, at the
        per-genome size this run's pool measured, it is at least the disk the run used and within 1.5 times it."""
        disk = json.loads((self.work / "disk-estimate.json").read_text())
        roots = {r.name: r for r in (self.work / "roots").iterdir()}
        counted = {n for n, r in roots.items() if not n.startswith("pool-")
                   and (r / "species" / "salmonella_enterica" / "ksets").is_dir()}
        self.assertEqual(set(disk["roots"]), counted)

        # In allocated blocks, as df and the check see disk (REVIEW-INPUT7 L4).
        used = run.disk_bytes

        pool = next(r for n, r in roots.items() if n.startswith("pool-"))
        ksets = list((pool / "species" / "salmonella_enterica" / "ksets").glob("*.kset"))
        per = sum(p.stat().st_size for p in ksets) / len(ksets)
        budget = SPEC["pool_budget"]
        root_ksets = sum(r["genomes"] for r in disk["roots"].values()) * per
        # The k31 roots share the pool's ksets (G1): only the k21 and k25 roots write their own.
        pending = sum(r["pending"] for r in disk["roots"].values()) * per
        self.assertEqual({r["k"] for r in disk["roots"].values() if r["pending"]}, {21, 25})
        estimate = (disk["genomes"] * per + pending + root_ksets * budget["root_overhead_share"]
                    + root_ksets * budget["arm_output_share"] + disk["pool_curated_bytes"] + budget["work_files_bytes"])
        # The whole WORK: roots, pool, arm outputs, logs, j/, k/ and the report (RAW is outside it).
        measured = used(self.work)
        report = {"per_genome_bytes": per, "pool_measured": used(pool / "species"), "measured_work": measured,
                  "arms_measured": used(self.work / "arms"), "estimate": round(estimate),
                  "pool_curated_measured": used(pool / "curated"), "pool_curated_estimate": disk["pool_curated_bytes"],
                  "work_files_measured": measured - used(self.work / "roots") - used(self.work / "arms"),
                  "arm_output_estimate": round(root_ksets * budget["arm_output_share"]),
                  "ratio_estimate_to_measured": estimate / measured,
                  "roots": {n: {"measured": used(roots[n] / "species"), "genomes": disk["roots"][n]["genomes"],
                                "shared": disk["roots"][n]["shared"],
                                "estimate": round((disk["roots"][n]["pending"] + disk["roots"][n]["genomes"]
                                                   * budget["root_overhead_share"]) * per)} for n in sorted(counted)}}
        (self.work / "disk-measured.json").write_text(json.dumps(report, indent=1, sort_keys=True))
        self.assertLessEqual(measured, estimate, report)
        self.assertLessEqual(estimate, 1.5 * measured, report)

    def test_a_reused_work_checks_disk_before_counting_a_new_root(self):
        """REVIEW-INPUT6 M1: the finished e2e WORK (a link copy) plus a k27 arm at 10^15 B per kset refuses, and
        no k27 kset is written."""
        copy_ = self.out / "work-k27"
        shutil.rmtree(copy_, ignore_errors=True)
        shutil.copytree(self.work, copy_, symlinks=True, copy_function=os.link)
        # The run writes WORK-level files (disk-estimate.json); copies, not links, keep the e2e WORK unchanged.
        for f in copy_.iterdir():
            if f.is_file() and not f.is_symlink():
                f.unlink()
                shutil.copy2(self.work / f.name, f)
        spec = json.loads((self.out / "spec-e2e.json").read_text())
        spec["arms"].append({"name": "k27", "k": 27})
        spec["pool_budget"]["default_kset_bytes"] = 10 ** 15
        (self.out / "spec-k27.json").write_text(json.dumps(spec))
        p = subprocess.run([sys.executable, "-B", str(HERE / "run.py"), str(self.out / "spec-k27.json"), "--raw",
                            str(self.out / "raw"), "--work", str(copy_), "--grm", os.environ["GRM_INPUT_E2E_BINARY"],
                            "--source-commit", "0" * 40, "--sigdrug-vocabulary", str(self.out / "raw" / "vocabulary")],
                           capture_output=True, text=True)
        try:
            self.assertEqual(p.returncode, 1, p.stderr[-2000:])
            self.assertIn("the comparison still has to write about", p.stderr)
            disk = json.loads((copy_ / "disk-estimate.json").read_text())
            (k27,) = [n for n in disk["roots"] if n.endswith("-k27")]
            self.assertEqual((disk["roots"][k27]["pending"], disk["pool_bytes"]), (2000, 0))
            self.assertTrue(all(r["pending"] == 0 for n, r in disk["roots"].items() if n != k27))
            self.assertEqual(list((copy_ / "roots" / k27).rglob("*.kset")), [])
        finally:
            shutil.rmtree(copy_, ignore_errors=True)

    def test_the_k31_roots_share_the_pool_s_ksets_and_count_none_anew(self):
        """H-1 gap G1 on REVIEW-INPUT5 M1(b): grm count reuses a kset only from its own root's kset folder, so the
        k31 roots get the pool's ksets as hard links before their count. After it, every kset and sidecar of every
        k31 root is still the pool's file (a recount writes a new one): 0 counted anew. The k21 and k25 roots
        count their own."""
        roots = {r.name: r for r in (self.work / "roots").iterdir()}
        pool = next(r for n, r in roots.items() if n.startswith("pool-"))
        pool_ksets = pool / "species" / "salmonella_enterica" / "ksets"
        k31 = 0
        for name, root in roots.items():
            if root is pool or not (root / "built.json").exists():
                continue
            ids = (root / "species" / "salmonella_enterica" / "genomes.txt").read_text().split()
            shared = json.loads((root / run.SHARED).read_text())
            mine = root / "species" / "salmonella_enterica" / "ksets"
            if name.endswith("-k31"):
                k31 += 1
                self.assertEqual(shared, {"genomes": len(ids), "shared": len(ids), "own": 0, "linked_now": len(ids),
                                          "pool": pool.name}, name)
                for g in ids:
                    for ext in ("kset", "json"):
                        self.assertFalse((mine / f"{g}.{ext}").is_symlink())
                        self.assertTrue(os.path.samefile(mine / f"{g}.{ext}", pool_ksets / f"{g}.{ext}"), (name, g))
            else:
                self.assertEqual((shared["shared"], shared["own"]), (0, len(ids)), name)
                g = ids[0]
                self.assertNotEqual(os.stat(mine / f"{g}.kset").st_ino, os.stat(pool_ksets / f"{g}.kset").st_ino)
        self.assertEqual(k31, 4)  # the default's, cap-1000's, cap-4000's and grr-v1's

    def test_sharing_changes_no_output(self):
        """G1's identity: the same comparison with --no-share-ksets (every root counts its own) gives byte-identical
        predictions, plans and J records, and the same K manifest and report after timings."""
        other = self.out / "work-noshare"
        shutil.rmtree(other, ignore_errors=True)
        subprocess.run([sys.executable, "-B", str(HERE / "run.py"), str(self.out / "spec-e2e.json"), "--raw",
                        str(self.out / "raw"), "--work", str(other), "--grm", os.environ["GRM_INPUT_E2E_BINARY"],
                        "--source-commit", os.environ.get("GRM_INPUT_E2E_COMMIT", "0" * 40), "--jobs",
                        os.environ.get("GRM_INPUT_E2E_JOBS", "3"), "--trace-training", "--sigdrug-vocabulary",
                        str(self.out / "raw" / "vocabulary"), "--no-share-ksets"], check=True, capture_output=True)
        self.assertEqual(timeless(other), timeless(self.work))
        files = sorted(str(p.relative_to(self.work)) for p in (self.work / "arms").rglob("*")
                       if p.is_file() and (p.name in ("predictions.jsonl", "training.jsonl") or p.parent.name == "plans"))
        self.assertGreater(len(files), 13 * 2)
        for name in files:
            self.assertEqual((other / name).read_bytes(), (self.work / name).read_bytes(), name)
        for root in (other / "roots").iterdir():
            if (root / run.SHARED).exists():
                self.assertEqual(json.loads((root / run.SHARED).read_text())["shared"], 0, root.name)
        self.assertLess(run.disk_bytes(self.work / "roots"), run.disk_bytes(other / "roots"))

    def test_a_reused_work_resumes_every_finished_arm_and_reruns_a_damaged_one(self):
        """G6 through the real grm: the finished e2e WORK (a copy, ksets linked) run again skips every arm but the
        one whose predictions were damaged, and gives the same J and K after timings."""
        copy_ = self.out / "work-resume"
        shutil.rmtree(copy_, ignore_errors=True)
        shutil.copytree(self.work, copy_, symlinks=True, copy_function=lambda s, d: (
            os.link(s, d) if Path(s).parent.name == "ksets" else shutil.copy2(s, d)))
        damaged = copy_ / "arms" / "budget-1000" / "salmonella_enterica__meropenem" / "predictions.jsonl"
        damaged.write_text(damaged.read_text()[:-200])
        p = subprocess.run([sys.executable, "-B", str(HERE / "run.py"), str(self.out / "spec-e2e.json"), "--raw",
                            str(self.out / "raw"), "--work", str(copy_), "--grm", os.environ["GRM_INPUT_E2E_BINARY"],
                            "--source-commit", os.environ.get("GRM_INPUT_E2E_COMMIT", "0" * 40), "--trace-training",
                            "--sigdrug-vocabulary", str(self.out / "raw" / "vocabulary")],
                           capture_output=True, text=True)
        try:
            self.assertEqual(p.returncode, 0, p.stderr[-2000:])
            summary = json.loads((copy_ / "summary.json").read_text())
            self.assertEqual(len(summary["resumed_arms"]), 12)
            self.assertNotIn("budget-1000/salmonella_enterica__meropenem", summary["resumed_arms"])
            self.assertEqual(damaged.read_bytes(), (self.work / damaged.relative_to(copy_)).read_bytes())
            self.assertEqual(timeless(copy_), timeless(self.work))
            twice = [log.name for log in (copy_ / "logs").glob("arm-*.log")
                     if log.read_text().count(" arm --species ") > 1]
            self.assertEqual(twice, ["arm-budget-1000-salmonella_enterica__meropenem.log"])
        finally:
            shutil.rmtree(copy_, ignore_errors=True)

    def test_the_gate_held_the_pinned_synthetic_values(self):
        """REVIEW-INPUT8 M1: the e2e ran with the gate on, so its default matched both pins through the real binary."""
        root = self.default_root()
        self.assertEqual(json.loads((root / "curated.json").read_text())["drawn_order_sha256"],
                         SYNTHETIC_GATE["genome_order"]["salmonella_enterica"])
        self.assertEqual(run.datasets_view(root / "curated" / "datasets.json"), SYNTHETIC_GATE["datasets_json_sha256"])
        entries = json.loads((root / "curated" / "datasets.json").read_text())
        self.assertTrue(all("selection" in e for e in entries))

    def test_the_manifest_hashes_the_binary_the_run_ran(self):
        """REVIEW-INPUT8 L2: the run's --grm was a bare name with a decoy in its cwd."""
        inputs = json.loads((self.work / "inputs.json").read_text())
        manifest = json.loads((self.work / "k" / "manifest.json").read_text())
        binary = Path(os.environ["GRM_INPUT_E2E_BINARY"]).resolve()
        self.assertEqual(inputs["grm"], str(binary))
        self.assertEqual(manifest["binary_sha256"], inputs["grm_sha256"])
        self.assertEqual(inputs["grm_sha256"], run.file_sha256(binary))

    def test_every_root_curates_with_the_species_check_and_the_vocabulary(self):
        """REVIEW-INPUT8 L4: every curated root, cohort arms included, links species-check/, curates with
        --sigdrug-vocabulary, and holds no genome the check fails."""
        (failed,) = json.loads((self.out / "raw" / "truth.json").read_text())["curate_inputs"]["species_check_fails"]
        roots = [r for r in (self.work / "roots").iterdir() if not r.name.startswith("pool-")]
        self.assertGreaterEqual(len(roots), 6)
        for root in roots:
            self.assertTrue((root / "species-check").is_symlink(), root.name)
            curate = [line for line in (self.work / "logs" / f"root-{root.name}.log").read_text().splitlines()
                      if " curate " in line]
            self.assertEqual(len(curate), 1, root.name)
            self.assertIn(f"--sigdrug-vocabulary {self.out / 'raw' / 'vocabulary'}", curate[0])
            self.assertNotIn(failed, (root / "species" / "salmonella_enterica" / "genomes.txt").read_text().split(),
                             root.name)

    def test_the_default_curates_as_revision_12_did(self):
        """REVIEW-INPUT7 M1: the default root links RAW's species check and curates with the vocabulary. It drops
        the genome the check fails, renames the misspelled records to meropenem, and equals a plain grm curate
        run as S0 ran it (species-check/ in the root, --sigdrug-vocabulary)."""
        truth = json.loads((self.out / "raw" / "truth.json").read_text())["curate_inputs"]
        default = self.default_root()
        self.assertTrue((default / "species-check").is_symlink())
        ids = set((default / "species" / "salmonella_enterica" / "genomes.txt").read_text().split())
        (failed,) = truth["species_check_fails"]
        self.assertNotIn(failed, ids)
        labels = {}
        for line in (default / "curated" / "phenotypes.jsonl").read_text().splitlines():
            p = json.loads(line)
            labels.setdefault(p["genome_id"], set()).add(p["antibiotic"])
        self.assertTrue(truth["misspelled"])
        for g in truth["misspelled"]:
            self.assertEqual(labels.get(g), {"meropenem"}, g)
        self.assertNotIn(failed, labels)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for part in ("bvbrc", "fasta", "catalog", "species-check"):
                (root / part).symlink_to(self.out / "raw" / part, target_is_directory=True)
            subprocess.run([os.environ["GRM_INPUT_E2E_BINARY"], "--data", str(root), "curate", "--cap", "2000",
                            "--min-class", "25", "--sigdrug-vocabulary", str(self.out / "raw" / "vocabulary")],
                           check=True, capture_output=True)
            self.assertEqual((root / "curated" / "datasets.json").read_bytes(),
                             (default / "curated" / "datasets.json").read_bytes())
            plain = set((root / "species" / "salmonella_enterica" / "genomes.txt").read_text().split())
            self.assertEqual(plain, ids)
        inputs = json.loads((self.work / "inputs.json").read_text())
        self.assertIsNotNone(inputs["vocabulary_sha256"])

    def test_an_arm_refuses_a_grouping_counted_at_another_k(self):
        """REVIEW-INPUT4 I5: the k21 arm's root with the k = 31 grouping."""
        groups = json.loads((self.work / "species-groups.json").read_text())["structure"]
        (root,) = [r for r in (self.work / "roots").iterdir()
                   if json.loads((r / "curated.json").read_text())["arm"] == "k21"]
        ds = "salmonella_enterica__meropenem"
        with tempfile.TemporaryDirectory() as tmp:
            p = subprocess.run([os.environ["GRM_INPUT_E2E_BINARY"], "--data", str(root), "arm", "--species",
                                "salmonella_enterica", "--dataset", ds, "--arm", "k21", "--blocks", "K", "--groups",
                                groups, "--out", tmp], capture_output=True, text=True)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("was counted at k = 31, this arm at k = 21", p.stderr)

    def test_every_cohort_arm_keeps_the_default_folds_of_its_shared_genomes(self):
        """REVIEW-INPUT3 N2: a genome both cohorts hold keeps its inner fold and both holdout memberships in every
        outer fold it trains in, and its outer fold whenever it is scored; whichever lowest ids the arm holds."""
        ds = "salmonella_enterica__meropenem"
        kept = 0
        for arm in ("cap-1000", "cap-4000", "grr-v1"):
            trace = (self.work / "arms" / arm / ds / "training.jsonl").read_text().splitlines()
            arm_ids = json.loads(trace[0])["sample_ids"]
            for split in ("random", "lineage"):
                ref = json.loads((self.work / "arms" / "default" / ds / "plans" / f"{split}-0.json").read_text())
                coh = json.loads((self.work / "arms" / arm / ds / "plans" / f"{split}-0.cohort.json").read_text())
                plan, at = ref["plan"], {g: i for i, g in enumerate(ref["sample_ids"])}
                if coh["unavailable_reason"] is not None or plan["unavailable_reason"] is not None:
                    continue
                for j, g in enumerate(arm_ids):
                    i = at.get(g)
                    if i is None:
                        continue
                    self.assertIn(coh["outer"][j], (plan["outer"][i], None), f"{arm} {split} {g}")
                    for f in range(plan["k"]):
                        if coh["inner"][f][j] is None or plan["inner"][f][i] is None:
                            continue
                        self.assertEqual(coh["inner"][f][j], plan["inner"][f][i], f"{arm} {split} {g} fold {f}")
                        for h in ("gbt", "lasso"):
                            self.assertEqual(j in coh["holdouts"][f][h], i in plan["holdouts"][f][h],
                                             f"{arm} {split} {g} {h} fold {f}")
                        kept += 1
        self.assertGreater(kept, 0)

    def test_genomes_are_in_canonical_order_and_no_cohort_arm_is_identical(self):
        for root in (self.work / "roots").iterdir():
            ids = (root / "species" / "salmonella_enterica" / "genomes.txt").read_text().split()
            self.assertEqual(ids, sorted(ids), root.name)
        self.assertEqual(self.summary["identical_cohorts"], [])
        self.assertGreater(self.summary["eligible_genomes"], self.summary["panel_floor"] - 1)

    def test_the_default_arm_folds_match_an_independent_id_keyed_reference(self):
        """The comparison keys every fold on genome identity (main's ruling on REVIEW-INPUT2 N1); this
        recomputes the default arm's outer folds, inner folds and holdouts from the ids alone, in Python."""
        import hashlib

        def keyed(key, seed, purpose):
            return int.from_bytes(hashlib.sha256(f"grm-keyed:{purpose}:{seed}:{key}".encode()).digest()[:8], "big")

        ds = "salmonella_enterica__meropenem"
        plans = self.work / "arms" / "default" / ds / "plans"
        rows = [json.loads(l) for l in (self.work / "arms" / "default" / ds / "predictions.jsonl").read_text()
                .splitlines() if l.strip()]
        for split in ("random", "lineage"):
            pf = json.loads((plans / f"{split}-0.json").read_text())
            ids, plan = pf["sample_ids"], pf["plan"]
            self.assertEqual(plan["format"], "grm-fold-plan-keyed-v1")
            unit = list(ids)
            if split == "lineage":
                low = {}
                for i, c in enumerate(pf["cluster"]):
                    low[c] = min(low.get(c, ids[i]), ids[i])
                unit = [low[c] for c in pf["cluster"]]
            else:
                for group in pf["duplicates"]:
                    low = min(ids[i] for i in group)
                    for i in group:
                        unit[i] = low
            outer = [keyed(u, plan["repeat_seed"], "outer") % plan["k"] for u in unit]
            self.assertEqual(outer, plan["outer"], split)
            folds = {r["sample_id"]: r["fold"] for r in rows if r["split"] == split}
            self.assertEqual(folds, dict(zip(ids, outer)), split)
            for f in range(plan["k"]):
                seed = (plan["repeat_seed"] + 101 * f) % 2 ** 64
                inner = [keyed(u, seed, "inner") % plan["inner_k"] if outer[i] != f else None
                         for i, u in enumerate(unit)]
                self.assertEqual(inner, plan["inner"][f], (split, f))
                for purpose, m in (("gbt", 7), ("lasso", 5)):
                    held = [i for i, u in enumerate(unit) if outer[i] != f and keyed(u, seed, purpose) % m == 0]
                    self.assertEqual(held, plan["holdouts"][f][purpose], (split, f, purpose))

    def test_run_arm_never_trains_a_cohort_arm_on_its_held_out_or_scored_genomes(self):
        """REVIEW-INPUT2 L4 (M17 at run_arm): what `run_arm` handed the learners, from training.jsonl,
        against the cohort plans it wrote."""
        ds = "salmonella_enterica__meropenem"
        checked = 0
        for arm in ("cap-1000", "cap-4000", "grr-v1"):
            out = self.work / "arms" / arm / ds
            head, *lines = (out / "training.jsonl").read_text().splitlines()
            cohort_ids = json.loads(head)["sample_ids"]
            self.assertEqual(cohort_ids, sorted(cohort_ids))
            trained = {}
            for line in lines:
                t = json.loads(line)
                trained[(t["split"], t["repeat"], t["fold"])] = set(t["train"])
            for path in sorted((out / "plans").glob("*.cohort.json")):
                plan = json.loads(path.read_text())
                split, repeat = path.name.split(".")[0].rsplit("-", 1)
                self.assertEqual(len(plan["outer"]), len(cohort_ids), path.name)
                for f in range(plan["k"]):
                    train = trained[(split, int(repeat), f)]
                    out_of_training = {cohort_ids[i] for i in plan["held_out"][f]}
                    out_of_training |= {cohort_ids[i] for i, o in enumerate(plan["outer"]) if o == f}
                    self.assertFalse(train & out_of_training, (arm, path.name, f))
                    checked += 1
        self.assertGreater(checked, 0)

    def test_a_shuffled_genome_list_gives_byte_identical_arm_output(self):
        import argparse
        import random
        import tempfile
        default = self.default_root()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "root"
            root.mkdir()
            for part in ("bvbrc", "fasta", "catalog"):
                (root / part).symlink_to((default / part).resolve(), target_is_directory=True)
            shutil.copytree(default / "curated", root / "curated")
            ids = (default / "species" / "salmonella_enterica" / "genomes.txt").read_text().split()
            random.Random(7).shuffle(ids)
            (root / "species" / "salmonella_enterica").mkdir(parents=True)
            (root / "species" / "salmonella_enterica" / "genomes.txt").write_text("\n".join(ids) + "\n")
            args = argparse.Namespace(grm=os.environ["GRM_INPUT_E2E_BINARY"], raw=str(default), work=tmp)
            runner = run.Runner(args, SPEC)
            runner.canonicalize(root)
            runner.build(root, 31)
            ds = "salmonella_enterica__meropenem"
            out = Path(tmp) / "arm"
            groups = json.loads((self.work / "species-groups.json").read_text())["structure"]
            subprocess.run([args.grm, "--data", str(root), "arm", "--species", "salmonella_enterica", "--dataset", ds,
                            "--arm", "default", "--blocks", "K", "--groups", groups, "--out", str(out)],
                           check=True, capture_output=True)
            mine = self.work / "arms" / "default" / ds
            self.assertEqual((out / "predictions.jsonl").read_bytes(), (mine / "predictions.jsonl").read_bytes())
            for plan in sorted((mine / "plans").iterdir()):
                self.assertEqual((out / "plans" / plan.name).read_bytes(), plan.read_bytes(), plan.name)

    def test_the_reloaded_report_rebuilds(self):
        from scripts.ablation import report as k_report
        k_report.load_report(self.work / "report" / "ablation-report.json", prediction_root=self.work)


if __name__ == "__main__":
    unittest.main()
