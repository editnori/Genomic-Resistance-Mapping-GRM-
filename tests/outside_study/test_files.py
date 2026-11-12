"""The scoring script end to end on files, inside a temporary committed
kit: the receipts it checks, the pins, the frozen gate and every way tried
around it, the outputs and byte-identical reruns."""

import contextlib
import gzip
import io
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests.outside_study import fixtures as F
from tests.outside_study.fixtures import S

OUTPUTS = S.OUTPUT_NAMES
PLAIN = tuple(n for n in OUTPUTS if n != "results.json")


class KitCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.w = self.make_world()
        self.reports = self.root / "reports"
        self.kit = F.kit_repo(self.root / "kit", self.w, self.reports)
        self.S = self.kit.S
        self.work = F.write_work(self.w, self.root / "work", self.kit)

    def make_world(self):
        return F.world()

    def tearDown(self):
        self.tmp.cleanup()

    def argv(self, out, *extra, commit=True):
        argv = ["--work", str(self.work), "--reports", str(self.reports), "--species", F.SPECIES,
                "--population", "primary", "--out", str(out)]
        return argv + (["--expect-commit", self.kit.commit] if commit else []) + list(extra)

    def main(self, out, *extra, commit=True):
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            code = self.S.main(self.argv(out, *extra, commit=commit))
        return code, err.getvalue()

    def refused(self, pattern, *extra, commit=True):
        out = self.root / "refused"
        code, err = self.main(out, *extra, commit=commit)
        self.assertEqual(code, 2, err)
        self.assertRegex(err, pattern)
        self.assertFalse(out.exists(), "a refusal writes nothing")

    def edit_work(self, rel, old, new):
        p = self.work / rel
        text = p.read_text()
        self.assertIn(old, text)
        p.write_text(text.replace(old, new, 1))


class Files(KitCase):
    def test_scores_and_writes_every_output(self):
        out = self.root / "out"
        code, err = self.main(out)
        self.assertEqual(code, 0, err)
        self.assertEqual(sorted(p.name for p in out.iterdir()), sorted(OUTPUTS))
        self.assertIn("123 of 126 genomes scored on 2 datasets", err)
        res = json.loads((out / "results.json").read_text())
        self.assertEqual(res["kit"]["commit"], self.kit.commit)
        self.assertEqual(res["kit"]["expect_commit"], self.kit.commit)
        self.assertEqual(res["kit"]["frozen_sha256"], S.sha256_file(self.kit.root / S.FROZEN_REL))
        self.assertEqual(sorted(res["kit"]["files"]), sorted(S.FROZEN_FILES))
        self.assertTrue(res["table_pins_checked"])
        self.assertEqual(res["attempts"], {s: 1 for s in S.STAGES})
        self.assertEqual(res["grm_sha256"], F.GRM_SHA)
        self.assertIn("work/receipts/predict.json", res["inputs"])
        self.assertIn("tables/pairs.tsv", res["inputs"])
        metrics = S.read_tsv(out / "metrics.tsv")
        scm = next(r for r in metrics if r["dataset"] == F.A and r["learner"] == "scm")
        self.assertEqual((scm["tp"], scm["fn"], scm["fp"], scm["tn"]), ("48", "12", "7", "53"))
        self.assertEqual((scm["sensitivity"], scm["sensitivity_lo"], scm["reading_sensitivity"]),
                         ("0.800000", "0.676700", "below"))
        self.assertEqual(S.read_tsv(out / "reading_counts.tsv")[0],
                         {"learner": "scm", "below": "1", "consistent": "1", "above": "0", "not_read": "2"})
        self.assertEqual(len(S.read_tsv(out / "paired.tsv")), 3)
        self.assertEqual(len(S.read_tsv(out / "duplicates.tsv")), 3)

    def test_the_attempts_are_counted(self):
        log = self.work / "receipts/attempts.jsonl"
        with open(log, "a") as f:
            f.write(json.dumps({"stage": "fetch", "event": "start"}) + "\n")
        out = self.root / "out"
        self.assertEqual(self.main(out)[0], 0)
        res = json.loads((out / "results.json").read_text())
        self.assertEqual(res["attempts"]["fetch"], 2)
        self.assertIn("Stage attempts: setup 1, fetch 2,", (out / "reading.txt").read_text())

    def test_reruns_are_byte_identical_in_separate_processes(self):
        # Each `python3 -I` process draws its own random hash seed.
        script = self.kit.root / "scripts/outside-study/score.py"
        outs = []
        for i in range(3):
            out = self.root / f"run{i}"
            p = subprocess.run([sys.executable, "-I", str(script), *self.argv(out)],
                               capture_output=True, text=True, env={"PATH": "/usr/bin:/bin"})
            self.assertEqual(p.returncode, 0, p.stderr)
            outs.append(out)
        for out in outs[1:]:
            for name in OUTPUTS:
                self.assertEqual((outs[0] / name).read_bytes(), (out / name).read_bytes(), name)

    def test_shuffled_input_rows_give_the_same_outputs(self):
        a = self.root / "a"
        self.assertEqual(self.main(a)[0], 0)
        rng = random.Random(7)
        w = self.w
        for d in w["rows"]:
            rng.shuffle(w["rows"][d])
        for key in ("availability", "verdicts", "counted"):
            items = list(w[key].items())
            rng.shuffle(items)
            w[key] = dict(items)
        calls = list(w["known_gene"]["calls"].items())
        rng.shuffle(calls)
        w["known_gene"]["calls"] = dict(calls)
        shutil.rmtree(self.work)
        F.write_work(w, self.work, self.kit)
        b = self.root / "b"
        self.assertEqual(self.main(b)[0], 0)
        for name in PLAIN:
            self.assertEqual((a / name).read_bytes(), (b / name).read_bytes(), name)
        ra, rb = (json.loads((x / "results.json").read_text()) for x in (a, b))
        self.assertNotEqual(ra["inputs"], rb["inputs"], "the work files' hashes differ")
        ra.pop("inputs"), rb.pop("inputs")
        self.assertEqual(ra, rb)

    def test_a_rehearsal_prefixes_every_output(self):
        w = F.rehearsal_world()
        tables = F.write_tables(w, self.root / "tables")
        shutil.rmtree(self.work)
        F.write_work(w, self.work, self.kit, rehearsal=True)
        out = self.root / "out"
        code, err = self.main(out, "--rehearsal", "--tables", str(tables), commit=False)
        self.assertEqual(code, 0, err)
        self.assertEqual(sorted(p.name for p in out.iterdir()), sorted("rehearsal-" + n for n in OUTPUTS))
        res = json.loads((out / "rehearsal-results.json").read_text())
        self.assertFalse(res["table_pins_checked"])
        self.assertIsNone(res["kit"]["expect_commit"])
        self.assertTrue((out / "rehearsal-reading.txt").read_text().startswith("REHEARSAL ON COHORT GENOMES"))

    def test_tables_need_a_rehearsal(self):
        self.refused("--tables is for a rehearsal only", "--tables", str(self.root))

    def test_no_option_moves_the_pins_or_the_kit_root(self):
        for flag in ("--pins", "--kit-root", "--table-pins", "--frozen"):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                self.S.main(self.argv(self.root / "x") + [flag, "elsewhere"])

    def test_reports_must_match_their_pins(self):
        p = self.reports / f"{F.A}.json"
        p.write_bytes(p.read_bytes().replace(b'"schema": 2', b'"schema":  2'))
        self.refused("is not the pinned")

    def test_known_gene_unavailable_only_in_a_rehearsal(self):
        (self.work / "known-gene" / "unavailable.txt").write_text("no catalog\n")
        r = json.loads((self.work / "receipts/known-gene.json").read_text())
        r["outputs"]["work/known-gene/unavailable.txt"] = S.sha256_file(self.work / "known-gene/unavailable.txt")
        (self.work / "receipts/known-gene.json").write_text(json.dumps(r))
        self.refused("only a rehearsal may go on without them")

    def test_rows_naming_another_report(self):
        rel = f"predict/{F.A}.jsonl"
        self.edit_work(rel, self.w["shas"][F.A], "0" * 64)
        F.write_receipts(self.work, self.kit, self.receipt_files())
        self.refused("names report sha256")

    def receipt_files(self):
        out = {}
        for stage in S.STAGES:
            r = json.loads((self.work / "receipts" / f"{stage}.json").read_text())
            out[stage] = [k[len("work/"):] for k in r["outputs"]]
        return out


class CohortInTheFrozenList(KitCase):
    def make_world(self):
        w = F.world()
        w["candidates"][0]["genome_id"] = "999.3"
        return w

    def test_a_cohort_genome_in_the_frozen_list_is_refused(self):
        self.refused("999.3 is a cohort genome")


class ReceiptChecks(KitCase):
    """score.py takes nothing from the work directory that the run's six
    receipts do not vouch for."""

    def receipt(self, stage):
        return json.loads((self.work / "receipts" / f"{stage}.json").read_text())

    def put(self, stage, r):
        (self.work / "receipts" / f"{stage}.json").write_text(json.dumps(r))

    def test_a_missing_receipt(self):
        (self.work / "receipts/count.json").unlink()
        self.refused("the count stage has no receipt")

    def test_a_receipt_from_another_commit(self):
        r = self.receipt("fetch")
        r["kit"]["commit"] = "0" * 40
        self.put("fetch", r)
        self.refused("the fetch receipt names commit 0{40}")

    def test_a_receipt_from_other_scripts(self):
        r = self.receipt("predict")
        r["kit"]["run.py"] = "0" * 64
        self.put("predict", r)
        self.refused("the predict stage ran a run.py with sha256 0{64}")

    def test_receipts_from_two_grm_binaries(self):
        r = self.receipt("count")
        r["grm"]["sha256"] = "8" * 64
        self.put("count", r)
        self.refused("more than one grm binary")

    def test_a_receipt_for_other_arguments(self):
        for key, value, pattern in (("species", "other_species", "species 'other_species'"),
                                    ("population", "same-study", "population 'same-study'"),
                                    ("rehearsal", True, "rehearsal True")):
            r = self.receipt("setup")
            saved = json.dumps(r)
            r["args"][key] = value
            self.put("setup", r)
            self.refused(pattern)
            self.put("setup", json.loads(saved))

    def test_a_work_file_changed_after_its_stage(self):
        self.edit_work("fetch/availability.tsv", "1001.4\tyes", "1001.4\tno")
        self.refused("work/fetch/availability.tsv has sha256 .* the fetch stage recorded")

    def test_a_work_file_no_receipt_names(self):
        r = self.receipt("count")
        del r["outputs"]["work/count/cohort-digests.tsv"]
        self.put("count", r)
        self.refused("work/count/cohort-digests.tsv is in no stage's receipt")

    def test_a_work_file_that_is_missing(self):
        (self.work / "count/counted.tsv").unlink()
        self.refused("work/count/counted.tsv is missing")

    def test_a_report_that_is_missing(self):
        (self.reports / f"{F.A}.json").unlink()
        self.refused(f"{F.A}: no report file at ")

    def test_an_extra_prediction_file(self):
        shutil.copyfile(self.work / f"predict/{F.A}.jsonl", self.work / f"predict/{F.C}.jsonl")
        self.refused("the predict folder holds")

    def test_survivors_that_are_not_the_records(self):
        self.edit_work("predict/survivors.txt", "1001.2\n", "")
        r = self.receipt("predict")
        r["outputs"]["work/predict/survivors.txt"] = S.sha256_file(self.work / "predict/survivors.txt")
        self.put("predict", r)
        self.refused("predict/survivors.txt lists 122 genomes, not the 123 survivors")

    def test_a_prediction_file_that_is_not_json(self):
        with open(self.work / f"predict/{F.A}.jsonl", "a") as f:
            f.write("{not json\n")
        r = self.receipt("predict")
        r["outputs"][f"work/predict/{F.A}.jsonl"] = S.sha256_file(self.work / f"predict/{F.A}.jsonl")
        self.put("predict", r)
        self.refused(f"work/predict/{F.A}.jsonl: not JSON")

    def test_no_attempts_log(self):
        (self.work / "receipts/attempts.jsonl").unlink()
        self.refused("left no receipts/attempts.jsonl")


class CommittedPins(unittest.TestCase):
    """The repository's own pin files against the repository's own files,
    so a change to a frozen file that leaves its pin behind fails here and
    not first where the study runs."""

    def test_the_pin_file_names_the_frozen_files_and_matches_them(self):
        pins = S.read_pins(F.ROOT / S.FROZEN_REL)
        self.assertEqual(sorted(pins), sorted(S.FROZEN_FILES))
        self.assertEqual(len(S.FROZEN_FILES), 6)
        for rel in S.FROZEN_FILES:
            self.assertEqual(pins[rel], S.sha256_file(F.ROOT / rel), rel)

    def test_the_frozen_tables_match_their_pins(self):
        pins = S.read_pins(F.KIT / "tables.sha256")
        self.assertEqual(sorted(pins), sorted(S.TABLE_FILES))
        for name in S.TABLE_FILES:
            data = (F.KIT / "frozen" / S.TABLE_STORED[name]).read_bytes()
            if name.endswith(".tsv") and S.TABLE_STORED[name].endswith(".gz"):
                data = gzip.decompress(data)
            self.assertEqual(pins[name], S.sha256_bytes(data), name)


class FrozenGate(KitCase):
    """Each way tried around the gate, in a temporary committed kit. A real
    run passes the published commit; every attack here leaves the kit able
    to score with code or pins the commit does not hold."""

    def gate(self, expect=None, **kw):
        return self.S.frozen_gate(self.kit.commit if expect is None else expect, **kw)

    def assertGateRefuses(self, pattern, expect=None):
        with self.assertRaisesRegex(self.S.Refusal, pattern):
            self.gate(expect)

    def path(self, rel):
        return self.kit.root / rel

    def edit_score(self):
        p = self.path(S.SCORE_REL)
        p.write_text(p.read_text().replace("\nMIN_CLASS = 25\n", "\nMIN_CLASS = 20\n"))

    def git(self, *args, env=None):
        return F.git(self.kit.root, *args, env=env)

    def test_passes_on_the_committed_kit(self):
        gate = self.gate()
        self.assertEqual(gate["commit"], self.kit.commit)
        self.assertEqual(gate["expect_commit"], self.kit.commit)
        self.assertEqual(gate["frozen_sha256"], S.sha256_file(self.path(S.FROZEN_REL)))
        self.assertEqual(sorted(gate["files"]), sorted(S.FROZEN_FILES))

    def test_a_real_run_needs_the_expected_commit(self):
        with self.assertRaisesRegex(self.S.Refusal, "a real run needs --expect-commit"):
            self.S.frozen_gate(None)
        for bad in ("abc", self.kit.commit.upper(), self.kit.commit + "0", self.kit.commit[:12]):
            with self.assertRaisesRegex(self.S.Refusal, "is not 40 lower-case hex"):
                self.S.frozen_gate(bad)
        self.assertGateRefuses("not the expected commit", expect="0" * 40)
        self.assertIsNone(self.S.frozen_gate(None, rehearsal=True)["expect_commit"])

    def test_an_edited_frozen_file(self):
        for rel in S.FROZEN_FILES:
            if rel == S.FROZEN_REL:
                continue
            p = self.path(rel)
            before = p.read_bytes()
            p.write_bytes(before + b"\n")
            self.assertGateRefuses(f"{rel} has sha256")
            p.write_bytes(before)
        self.gate()

    def test_edited_and_repinned_but_not_committed(self):
        self.edit_score()
        F.repin(self.kit.root)
        self.assertGateRefuses("FROZEN.sha256 differs from the one committed at HEAD")

    def test_assume_unchanged_and_skip_worktree(self):
        self.edit_score()
        F.repin(self.kit.root)
        for flag in ("--assume-unchanged", "--skip-worktree"):
            self.git("update-index", flag, S.SCORE_REL, S.FROZEN_REL)
            self.assertEqual(self.git("status", "--short", "--untracked-files=no").stdout, "")
            self.assertGateRefuses("FROZEN.sha256 differs")
            self.git("update-index", flag.replace("--", "--no-"), S.SCORE_REL, S.FROZEN_REL)

    def test_an_edit_staged_but_not_committed(self):
        self.edit_score()
        F.repin(self.kit.root)
        self.git("add", "-A")
        self.assertGateRefuses("FROZEN.sha256 differs from the one committed at HEAD")

    def test_an_edit_committed_is_another_commit(self):
        old = self.kit.commit
        self.edit_score()
        F.repin(self.kit.root)
        F.commit_all(self.kit, "retune")
        self.assertGateRefuses(f"HEAD is {self.kit.commit}, not the expected commit {old}", expect=old)

    def test_a_committed_file_its_pin_does_not_match(self):
        self.edit_score()
        F.commit_all(self.kit, "edit without repin")
        self.assertGateRefuses("score.py at HEAD has sha256")

    def test_edited_pin_files(self):
        for rel in ("scripts/outside-study/reports.sha256", "scripts/outside-study/tables.sha256"):
            p = self.path(rel)
            before = p.read_text()
            p.write_text("0" + before[1:])
            self.assertGateRefuses(f"{rel} has sha256")
            p.write_text(before)

    def test_the_rules_text(self):
        with open(self.path("scripts/outside-study/frozen/RULES.md"), "a") as f:
            f.write("one more rule\n")
        self.assertGateRefuses("RULES.md has sha256")

    def test_git_variables_naming_another_repository_or_index(self):
        self.edit_score()
        F.repin(self.kit.root)
        side = self.root / "side"
        subprocess.run(["git", "clone", "-q", "--no-checkout", str(self.kit.root), str(side)], check=True)
        env = dict(os.environ, GIT_DIR=str(side / ".git"), GIT_WORK_TREE=str(self.kit.root))
        F.git(side, "add", "-A", env=env)
        F.git(side, "commit", "-q", "-m", "side", env=env)
        index = self.root / "index"
        shutil.copyfile(self.kit.root / ".git/index", index)
        F.git(self.kit.root, "update-index", "--assume-unchanged", S.SCORE_REL, S.FROZEN_REL,
              env=dict(os.environ, GIT_INDEX_FILE=str(index)))
        for extra in ({"GIT_DIR": str(side / ".git"), "GIT_WORK_TREE": str(self.kit.root)},
                      {"GIT_INDEX_FILE": str(index)},
                      {"GIT_CONFIG_PARAMETERS": "'core.worktree'='/'"}):
            saved = {k: os.environ.get(k) for k in extra}
            os.environ.update(extra)
            try:
                self.assertGateRefuses("FROZEN.sha256 differs")
            finally:
                for k, v in saved.items():
                    os.environ.pop(k) if v is None else os.environ.__setitem__(k, v)

    def test_a_stand_in_git_on_the_path(self):
        self.edit_score()
        F.repin(self.kit.root)
        fake = self.root / "bin"
        fake.mkdir()
        (fake / "git").write_text(f"#!/bin/sh\ncase \"$*\" in *rev-parse*) echo {self.kit.commit};; esac\nexit 0\n")
        (fake / "git").chmod(0o755)
        saved = os.environ["PATH"]
        os.environ["PATH"] = f"{fake}:{saved}"
        try:
            self.assertGateRefuses("FROZEN.sha256 differs")
        finally:
            os.environ["PATH"] = saved
        self.assertTrue(all(not p.startswith(str(fake)) for p in S.GIT_PATHS))

    def test_replace_objects(self):
        self.edit_score()
        F.repin(self.kit.root)
        self.git("add", "-A")
        tree = self.git("write-tree").stdout.strip()
        other = self.git("commit-tree", tree, "-m", "replacement").stdout.strip()
        self.git("reset", "-q")
        self.git("replace", self.kit.commit, other)
        self.assertEqual(self.git("show", f"HEAD:{S.FROZEN_REL}").stdout, self.path(S.FROZEN_REL).read_text(),
                         "plain git now shows the edited pins at HEAD")
        self.assertGateRefuses("FROZEN.sha256 differs")

    def test_a_work_tree_elsewhere(self):
        self.git("config", "core.worktree", str(self.root))
        self.assertGateRefuses("git puts the top of the work tree at")

    def test_a_copy_of_the_kit_without_git(self):
        copy = self.root / "copy"
        shutil.copytree(self.kit.root, copy, ignore=shutil.ignore_patterns(".git"))
        S2 = F.load("score", copy / "scripts/outside-study")
        with self.assertRaisesRegex(S2.Refusal, "is not in a git work tree"):
            S2.frozen_gate(self.kit.commit)

    def test_a_copy_of_the_kit_inside_another_repository(self):
        copy = self.root / "outer" / "inner"
        shutil.copytree(self.kit.root, copy, ignore=shutil.ignore_patterns(".git"))
        subprocess.run(["git", "init", "-q", str(self.root / "outer")], check=True)
        S2 = F.load("score", copy / "scripts/outside-study")
        with self.assertRaisesRegex(S2.Refusal, "not at the kit root"):
            S2.frozen_gate(self.kit.commit)

    def test_a_modified_copy_of_score_py_elsewhere(self):
        d = self.root / "x" / "y"
        d.mkdir(parents=True)
        shutil.copyfile(self.path(S.SCORE_REL), d / "score.py")
        S2 = F.load("score", d)
        with self.assertRaisesRegex(S2.Refusal, "not in a git work tree|is not the kit's own"):
            S2.frozen_gate(self.kit.commit)
        d2 = self.kit.root / "scripts" / "other"
        d2.mkdir()
        shutil.copyfile(self.path(S.SCORE_REL), d2 / "score.py")
        S3 = F.load("score", d2)
        with self.assertRaisesRegex(S3.Refusal, "is not the kit's own"):
            S3.frozen_gate(self.kit.commit)

    def test_symbolic_links(self):
        p = self.path(S.SCORE_REL)
        real = self.root / "score-real.py"
        shutil.copyfile(p, real)
        p.unlink()
        p.symlink_to(real)
        self.assertGateRefuses("is a symbolic link|is not the kit's own")
        p.unlink()
        shutil.copyfile(real, p)
        plan = self.path("docs/outside-study-test.md")
        shutil.move(plan, self.root / "plan-real.md")
        plan.symlink_to(self.root / "plan-real.md")
        self.assertGateRefuses("outside-study-test.md is a symbolic link")
        plan.unlink()
        shutil.move(self.root / "plan-real.md", plan)
        self.gate()
        docs = self.path("docs")
        shutil.move(docs, self.root / "docs-real")
        docs.symlink_to(self.root / "docs-real")
        self.assertGateRefuses("docs is a symbolic link")

    def test_pins_for_files_the_gate_does_not_know_or_missing_pins(self):
        p = self.path(S.FROZEN_REL)
        lines = p.read_text().splitlines(keepends=True)
        for text, pattern in (("".join(lines) + f"{'0' * 64}  extra.txt\n", "pins extra.txt, which is not a frozen file"),
                              ("".join(lines[1:]), "does not pin")):
            p.write_text(text)
            F.commit_all(self.kit, "pins")
            self.assertGateRefuses(pattern)

    def test_a_frozen_file_missing_from_the_commit(self):
        self.git("rm", "-q", "--cached", "docs/outside-study-test.md")
        self.git("commit", "-q", "-m", "drop")
        self.assertGateRefuses("docs/outside-study-test.md is not committed at HEAD")

    def test_no_pin_file(self):
        self.path(S.FROZEN_REL).unlink()
        self.assertGateRefuses("FROZEN.sha256 is missing")

    def test_the_frozen_tables_are_read_from_the_kit_and_checked(self):
        p = self.path("scripts/outside-study/frozen/pairs.tsv.gz")
        data = gzip.decompress(p.read_bytes()) + "\t".join(["1001.1", F.C, "S"] + [""] * 9).encode() + b"\n"
        p.write_bytes(gzip.compress(data, mtime=0))
        with self.assertRaisesRegex(self.S.Refusal, "the frozen table pairs.tsv has sha256"):
            self.S.load_tables(False)
        with self.assertRaisesRegex(self.S.Refusal, "--tables is for a rehearsal only"):
            self.S.load_tables(False, self.root)


if __name__ == "__main__":
    unittest.main()
