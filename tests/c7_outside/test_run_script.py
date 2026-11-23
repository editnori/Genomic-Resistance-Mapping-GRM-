"""The run script's own rules (scripts/c7-outside/run_species.sh), on SYNTHETIC files.

No genome, call or result is read: the expect lists, species files, logs and the stand-in for the runner are
written by the test.
"""

import hashlib
import json
import re
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests.c7_outside import common as C

SCRIPT = (C.KIT / "run_species.sh").read_text()
expect = C.load("expect")


def function(name):
    m = re.search(rf"^{name}\(\) \{{[^\n]*\n.*?^\}}\n", SCRIPT, re.S | re.M)
    if m is None:
        raise AssertionError(f"run_species.sh defines no {name}()")
    return m.group(0)


def fits_listed(expect_tsv):
    """The fit lines the script's awk takes from expect.tsv."""
    m = re.search(r"^awk (-F'\\t' '[^']*') \"\$OUT/expect.tsv\"", SCRIPT, re.M)
    if m is None:
        raise AssertionError("run_species.sh lists no fits from expect.tsv")
    prog = m.group(1).split(" ", 1)[1].strip("'")
    out = subprocess.run(["awk", "-F\t", prog, str(expect_tsv)], capture_output=True, text=True, check=True).stdout
    return [line.split() for line in out.splitlines()]


class Fits(unittest.TestCase):
    def test_the_fits_are_each_arm_at_its_seeds_and_the_study_never_the_default_or_twin(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "expect.tsv"
            rows = expect.rows(C.SPECIES, [C.DATASET, C.OTHER], [C.DATASET], [42, 43, 44], [45, 46, 47], [42])
            path.write_text("\t".join(expect.COLUMNS) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))
            fits = fits_listed(path)
        self.assertEqual(len(fits), 3 * 2 + 3)
        self.assertEqual({f[0] for f in fits}, {"xgboost", "lightgbm", "elasticnet"})
        self.assertEqual(sorted(f[1:] for f in fits if f[2] == "study-whole"), [[C.DATASET, "study-whole", "42"]] * 3)
        self.assertTrue(all(f[3] == "42" for f in fits))

    def test_a_finished_fit_is_skipped_and_a_failed_one_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            calls, fail = tmp / "calls", tmp / "fail"
            # A stand-in for the runner's python: it counts its calls and writes its --out, or fails.
            fake = tmp / "python"
            fake.write_text(f'#!/bin/bash\necho run >> "{calls}"\n[ -e "{fail}" ] && exit 1\n'
                            'while [ $# -gt 0 ]; do\n  [ "$1" = --out ] && { echo "{}" > "$2"; exit 0; }\n  shift\ndone\n'
                            'exit 3\n')
            fake.chmod(fake.stat().st_mode | stat.S_IXUSR)
            (tmp / "release" / "species").mkdir(parents=True)
            (tmp / "release" / "species" / f"{C.SPECIES}.json").write_text(json.dumps({"genomes": [{"id": "1234.5"}]}))
            out = tmp / "out"
            (out / "logs" / "fits").mkdir(parents=True)
            env = {"PATH": "/usr/bin:/bin", "HOME": str(tmp), "OUT": str(out), "W": str(tmp / "work"), "PY": str(fake),
                   "HERE": str(C.KIT), "UNITIGS": str(tmp / "u"), "PLANS": str(tmp / "plans"), "R": str(tmp / "release"),
                   "SP": C.SPECIES}
            cmd = function("scrub") + function("fit") + 'fit xgboost "$1" lineage 42'

            def go():
                return subprocess.run(["bash", "-c", cmd, "x", C.DATASET], env=env, capture_output=True, text=True)

            self.assertEqual(go().returncode, 0)
            self.assertEqual(go().returncode, 0)
            self.assertEqual(calls.read_text().count("run"), 1)
            report = out / "results" / "outside" / "xgboost" / f"{C.DATASET}.lineage-s42.outside.json"
            self.assertTrue(report.exists())
            report.unlink()
            fail.write_text("")
            self.assertEqual(go().returncode, 1)
            self.assertEqual(calls.read_text().count("run"), 2)
            self.assertTrue((out / "logs" / "fits" / f"xgboost-{C.DATASET}-lineage-s42.time").exists())

    def test_the_elastic_net_runs_after_the_boosters_as_many_at_once_as_memory_holds(self):
        m = re.search(r"^STEP=fits\n.*?^done\n", SCRIPT, re.S | re.M)
        self.assertIsNotNone(m, "run_species.sh has no fits block")
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            rows = expect.rows(C.SPECIES, [C.DATASET, C.OTHER], [], [42], [45], [42])
            (tmp / "expect.tsv").write_text("\t".join(expect.COLUMNS) + "\n"
                                            + "".join("\t".join(map(str, r)) + "\n" for r in rows))
            unitigs = tmp / "u"
            unitigs.mkdir()
            env = {"PATH": "/usr/bin:/bin", "OUT": str(tmp), "W": str(tmp), "UNITIGS": str(unitigs), "JOBS": "16"}
            cmd = ('say() { echo "$*"; }\nfit() { echo "$1" >> "$OUT/calls"; }\nexport -f fit\n'
                   'ENET_JOBS=$1\n' + m.group(0))

            def go(size, given=""):
                with open(unitigs / "patterns.bin", "wb") as f:
                    f.truncate(size)  # a sparse file: only its size is read
                (tmp / "calls").unlink(missing_ok=True)
                p = subprocess.run(["bash", "-c", cmd, "x", given], env=env, capture_output=True, text=True, check=True)
                return p.stdout, (tmp / "calls").read_text().split()

            out, calls = go(1 << 43)
            self.assertIn("the boosters 16 at once, the elastic net 1 at once", out)
            self.assertEqual(sorted(calls[:4]), ["lightgbm", "lightgbm", "xgboost", "xgboost"])
            self.assertEqual(calls[4:], ["elasticnet", "elasticnet"])
            self.assertIn("the elastic net 16 at once", go(8)[0])
            self.assertIn("the elastic net 3 at once", go(1 << 43, "3")[0])

    def test_the_elastic_net_can_be_left_out(self):
        m = re.search(r"^STEP=fits\n.*?^done\n", SCRIPT, re.S | re.M)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            rows = expect.rows(C.SPECIES, [C.DATASET], [C.DATASET], [42, 43, 44], [45], [42, 43, 44],
                               ("xgboost", "lightgbm"))
            (tmp / "expect.tsv").write_text("\t".join(expect.COLUMNS) + "\n"
                                            + "".join("\t".join(map(str, r)) + "\n" for r in rows))
            env = {"PATH": "/usr/bin:/bin", "OUT": str(tmp), "W": str(tmp), "UNITIGS": str(tmp / "none"), "JOBS": "4",
                   "ENET": "leave-out"}
            cmd = 'say() { echo "$*"; }\nfit() { echo "$1" >> "$OUT/calls"; }\nexport -f fit\nENET_JOBS=\n' + m.group(0)
            p = subprocess.run(["bash", "-c", cmd], env=env, capture_output=True, text=True, check=True)
            calls = (tmp / "calls").read_text().split()
        self.assertIn("the elastic net left out", p.stdout)
        self.assertEqual(sorted(set(calls)), ["lightgbm", "xgboost"])
        self.assertEqual(len(calls), 2 * 4)
        self.assertIn('  leave-out) ARMS=xgboost,lightgbm ;;', SCRIPT)
        self.assertIn('"$SEEDS" "$TWIN" "$ARM_SEEDS" \\\n  "$ARMS" > "$OUT/expect.tsv.part"', SCRIPT)

    def test_a_start_with_the_elastic_net_left_out_restarts_with_it_and_nothing_else_changes(self):
        # The lead's two passes on one OUT: the boosters with ENET=leave-out, then a restart with ENET=run.
        self.assertIn('lists_check "$OUT/lists.txt" "$OUT/lists.txt.part"\nmv "$OUT/lists.txt.part" "$OUT/lists.txt"\n',
                      SCRIPT)
        self.assertIn("printf 'species\\t%s\\ndatasets\\t%s\\nseeds\\t%s\\ntwin\\t%s\\narm_seeds\\t%s\\narms\\t%s\\n'", SCRIPT)

        def lists(arms, seeds="42,43,44", datasets=C.DATASET):
            return (f"species\t{C.SPECIES}\ndatasets\t{datasets}\nseeds\t{seeds}\ntwin\t45,46,47\n"
                    f"arm_seeds\t{seeds}\narms\t{arms}\n")

        boosters, all_three = "xgboost,lightgbm", "xgboost,lightgbm,elasticnet"
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)

            def go(first, this):
                (tmp / "first").unlink(missing_ok=True)
                if first is not None:
                    (tmp / "first").write_text(first)
                (tmp / "this").write_text(this)
                p = subprocess.run(["bash", "-c", 'say() { echo "$*"; }\n' + function("lists_check")
                                    + 'lists_check "$1" "$2"', "x", str(tmp / "first"), str(tmp / "this")],
                                   env={"PATH": "/usr/bin:/bin"}, capture_output=True, text=True)
                return p.returncode, p.stdout

            self.assertEqual(go(None, lists(boosters)), (0, ""))
            self.assertEqual(go(lists(all_three), lists(all_three)), (0, ""))
            rc, out = go(lists(boosters), lists(all_three))
            self.assertEqual(rc, 0)
            self.assertIn("the elastic net, left out at the first start, is listed now", out)
            for first, this in ((lists(all_three), lists(boosters)), (lists(boosters), lists(all_three, seeds="42,43")),
                                (lists(boosters), lists(all_three, datasets=f"{C.DATASET},{C.OTHER}")),
                                (lists(boosters), lists("xgboost,lightgbm,elasticnet,other"))):
                rc, out = go(first, this)
                self.assertEqual(rc, 1, this)
                self.assertIn("differ from the first start's; use a fresh OUT", out)

    def test_the_timing_table_is_written_from_each_fits_timing_file(self):
        # The step ran after every fit; before this test it stopped every run that fitted anything with a SyntaxError.
        m = re.search(r"^STEP=timing\n.*?(?=^STEP=score\n)", SCRIPT, re.S | re.M)
        self.assertIsNotNone(m, "run_species.sh has no timing block")
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "logs" / "fits").mkdir(parents=True)
            for name, t in (("xgboost-a-lineage-s42", (1.5, 0.72, 129552)), ("elasticnet-a-lineage-s43", (610.0, 600.63, 318480))):
                (tmp / "logs" / "fits" / f"{name}.json").write_text(
                    json.dumps({"wall_s": t[0], "cpu_s": t[1], "peak_kb": t[2]}))
            (tmp / "logs" / "fits" / "xgboost-a-lineage-s42.log").write_text("not a timing file\n")
            p = subprocess.run(["bash", "-c", "set -Eeuo pipefail\n" + m.group(0)], env={"PATH": "/usr/bin:/bin", "OUT": str(tmp)},
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertEqual((tmp / "logs" / "timing.tsv").read_text(),
                             "fit\twall_s\tcpu_s\tpeak_kb\nelasticnet-a-lineage-s43\t610.0\t600.63\t318480\n"
                             "xgboost-a-lineage-s42\t1.5\t0.72\t129552\n")

    def test_a_fits_log_comes_back_with_each_genome_id_masked(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "release" / "species").mkdir(parents=True)
            (tmp / "release" / "species" / f"{C.SPECIES}.json").write_text(
                json.dumps({"genomes": [{"id": "1234.5"}, {"id": "99.01"}]}))
            log = tmp / "logs" / "a.log"
            log.parent.mkdir()
            log.write_text("fold 1234.5 and 99.01, not 1234.56 or v99.01\n")
            env = {"PATH": "/usr/bin:/bin", "W": str(tmp / "work"), "R": str(tmp / "release"), "SP": C.SPECIES}
            subprocess.run(["bash", "-c", function("scrub") + 'scrub "$1"', "x", str(log)], env=env, check=True)
            self.assertEqual(log.read_text(), "fold [genome id] and [genome id], not 1234.56 or v99.01\n")
            self.assertIn("1234.5 ", (tmp / "work" / "raw-logs" / "a.log").read_text())


class Stage1(unittest.TestCase):
    """c3_stage on SYNTHETIC C3 run folders: read and checked on the C3 machine, run only on a new box."""

    COMMIT = "8bfd4ca8a64678c69431a6bc22b869af6d5cb81b"

    def c3_run(self, top, done=True, commit=COMMIT, score_changed=False, other_changed=False):
        run = top / "c3-run"
        (run / "results").mkdir(parents=True)
        (run / "score").mkdir()
        (run / "results" / "a.json").write_text("{}")
        (run / "score" / "count.tsv").write_text("sealed")
        grm = top / "grm"
        grm.write_text("binary")
        sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()  # noqa: E731
        (run / "binary.txt").write_text(f"commit {commit}\ngrm {sha(grm)}\ncode {'0' * 64}\n")
        (run / "SHA256SUMS").write_text(f"{sha(run / 'results' / 'a.json')}  results/a.json\n"
                                        f"{sha(run / 'score' / 'count.tsv')}  score/count.tsv\n"
                                        f"{sha(run / 'binary.txt')}  binary.txt\n")
        if done:
            (run / "DONE").write_text(f"DONE {C.SPECIES} 2026-11-23T00:00:00Z binary x code y results 1 SHA256SUMS "
                                      f"{sha(run / 'SHA256SUMS')}\n")
        if score_changed:
            (run / "score" / "count.tsv").write_text("changed")
        if other_changed:
            (run / "results" / "a.json").write_text("[]")
        return run, grm

    def go(self, top, run, grm, link=True):
        out = top / "out"
        out.mkdir()
        c3 = out / "c3"
        if link:
            c3.symlink_to(run)
        else:
            run.rename(c3)
        calls = top / "calls"
        repo = top / "repo"
        (repo / "scripts" / "c3-inputs").mkdir(parents=True)
        (repo / "scripts" / "c3-inputs" / "run_species.sh").write_text(f'echo stage1 >> "{calls}"\nexit 1\n')
        env = {"PATH": "/usr/bin:/bin", "C3": str(c3), "GRM": str(grm), "REPO": str(repo), "C3_COMMIT": self.COMMIT,
               "SP": C.SPECIES}
        cmd = 'say() { echo "$*"; }\n' + function("c3_stage") + "c3_stage"
        p = subprocess.run(["bash", "-e", "-c", cmd], env=env, capture_output=True, text=True)
        return p.returncode, p.stdout, calls.exists()

    def test_on_the_c3_machine_a_done_run_at_8bfd4ca_is_read_and_stage_one_never_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            top = Path(tmp)
            rc, out, ran = self.go(top, *self.c3_run(top, score_changed=True))
            self.assertEqual((rc, ran), (0, False), out)
            self.assertIn("read and not run", out)
            self.assertIn("2 files outside score/ held to SHA256SUMS", out)

    def test_on_the_c3_machine_a_run_with_no_done_is_refused_and_not_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            top = Path(tmp)
            rc, out, ran = self.go(top, *self.c3_run(top, done=False))
            self.assertEqual((rc != 0, ran), (True, False), out)
            self.assertIn("has no DONE", out)

    def test_a_c3_run_of_another_commit_or_a_changed_file_is_refused(self):
        for kw, msg in (({"commit": "0" * 40}, "was built at"), ({"other_changed": True}, "outside score/ is not")):
            with tempfile.TemporaryDirectory() as tmp:
                top = Path(tmp)
                rc, out, ran = self.go(top, *self.c3_run(top, **kw))
                self.assertEqual((rc != 0, ran), (True, False), out)
                self.assertIn(msg, out)

    def test_on_a_new_box_stage_one_is_c3s_run_with_the_unitig_arm_and_grms_graph(self):
        self.assertRegex(function("c3_stage"), r'bash "\$REPO/scripts/c3-inputs/run_species.sh" "\$KIT" "\$SP" '
                                               r'"\$DATASETS" "\$C3" "\$JOBS" "\$SEEDS" "\$TWIN" unitigs "" grm')
        with tempfile.TemporaryDirectory() as tmp:
            top = Path(tmp)
            rc, out, ran = self.go(top, *self.c3_run(top, done=False), link=False)
            self.assertEqual((rc != 0, ran), (True, True), out)


class PlansCheck(unittest.TestCase):
    """plans.py on SYNTHETIC exports, default results and release reports."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.top = Path(self._tmp.name)
        self.rows = expect.rows(C.SPECIES, [C.DATASET], [], [42, 43], [45], [42, 43])
        self.expect = self.top / "expect.tsv"
        self.expect.write_text("\t".join(expect.COLUMNS) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in self.rows))
        for s in (42, 43):
            C.write_export(self.top / "work" / f"export-s{s}", seed=s)
            self.default(s)
        rel = self.top / "release" / "datasets"
        rel.mkdir(parents=True)
        (rel / f"{C.DATASET}.json").write_text(json.dumps({"input_fingerprint": "synthetic-" + C.DATASET}))

    def tearDown(self):
        self._tmp.cleanup()

    def default(self, s, change=None):
        """The default's SYNTHETIC result at seed s, on the export's lineage folds by class, counted here."""
        meta = json.loads((self.top / "work" / f"export-s{s}" / C.DATASET / "dataset.json").read_text())
        plan = next(x["plan"] for x in meta["splits"] if x["plan"]["scheme"] == "lineage")
        y = meta["resistant"]
        folds = []
        for f in range(plan["k"]):
            test = [o == f for o in plan["outer"]]
            folds.append({"fold": f, "train": test.count(False), "train_resistant": sum(r and not t for r, t in zip(y, test)),
                          "test": test.count(True), "test_resistant": sum(r and t for r, t in zip(y, test))})
        doc = {"folds": folds, "tested": len(y), "tested_resistant": sum(y), "unavailable_reason": None}
        if change:
            change(doc)
        path = self.top / expect.picked("k31", C.DATASET, "lineage", s)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(doc))
        self.report(s, meta["samples"], plan["outer"])

    def report(self, s, ids, test_fold, seed=None):
        """The default's SYNTHETIC grm learn report at seed s where grm keeps it: its samples and lineage test folds."""
        species = self.top / "c3" / "work" / "root" / "species" / C.SPECIES
        path = species / ("results" if s == 42 else f"draws/seed-{s}") / f"{C.DATASET}.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps({"eval_config": {"seed": s if seed is None else seed},
                                    "samples": [{"id": i} for i in ids],
                                    "splits": [{"scheme": "lineage", "test_fold": list(test_fold)},
                                               {"scheme": "random", "test_fold": [0] * len(ids)}]}))

    def go(self):
        p = subprocess.run([sys.executable, "-I", "-B", str(C.KIT / "plans.py"), str(self.expect), str(self.top),
                            str(self.top / "release"), str(self.top / "work" / "export-s{seed}")],
                           capture_output=True, text=True)
        return p.returncode, p.stdout + p.stderr

    def test_exports_on_the_defaults_folds_pass(self):
        rc, out = self.go()
        self.assertEqual(rc, 0, out)
        self.assertIn("plans: 7 of 7 checks pass, 2 dataset seeds", out)

    def test_a_plan_that_differs_from_the_defaults_by_class_is_refused(self):
        def change(doc):
            doc["folds"][0]["test_resistant"] += 1
            doc["folds"][0]["train_resistant"] -= 1
        self.default(43, change)
        rc, out = self.go()
        self.assertEqual(rc, 1, out)
        self.assertIn(f"failed: {C.DATASET} seed 43: the plan's folds by class are the default's", out)

    def test_two_genomes_of_one_class_swapped_between_folds_pass_by_class_and_are_refused_by_genome(self):
        meta = json.loads((self.top / "work" / "export-s43" / C.DATASET / "dataset.json").read_text())
        outer = next(x["plan"] for x in meta["splits"] if x["plan"]["scheme"] == "lineage")["outer"]
        y = meta["resistant"]
        a = next(i for i in range(len(y)) if outer[i] == 0)
        b = next(i for i in range(len(y)) if outer[i] == 1 and y[i] == y[a])
        swapped = list(outer)
        swapped[a], swapped[b] = outer[b], outer[a]
        self.report(43, meta["samples"], swapped)
        rc, out = self.go()
        self.assertEqual(rc, 1, out)
        self.assertIn("plans: 6 of 7 checks pass", out)
        self.assertIn(f"failed: {C.DATASET} seed 43: the plan's folds by genome are the default's", out)
        self.assertNotIn("by class are the default's", out)
        self.assertFalse(any(i in out for i in meta["samples"]), "plans.py printed a genome id")

    def test_a_default_report_of_another_seed_or_none_is_refused(self):
        meta = json.loads((self.top / "work" / "export-s43" / C.DATASET / "dataset.json").read_text())
        outer = next(x["plan"] for x in meta["splits"] if x["plan"]["scheme"] == "lineage")["outer"]
        self.report(43, meta["samples"], outer, seed=42)
        rc, out = self.go()
        self.assertEqual(rc, 1, out)
        self.assertIn("the report is of seed 42, not 43", out)
        (self.top / "c3" / "work" / "root" / "species" / C.SPECIES / "draws" / "seed-43" / f"{C.DATASET}.json").unlink()
        rc, out = self.go()
        self.assertEqual(rc, 1, out)
        self.assertIn(f"failed: {C.DATASET} seed 43: the default's report at seed 43 is read", out)

    def test_an_export_of_another_seed_is_refused(self):
        import shutil
        shutil.rmtree(self.top / "work" / "export-s43")
        shutil.copytree(self.top / "work" / "export-s42", self.top / "work" / "export-s43")
        rc, out = self.go()
        self.assertEqual(rc, 1, out)
        self.assertIn("seed 42, not 43", out)

    def test_a_seed_42_fingerprint_that_is_not_the_releases_is_refused(self):
        (self.top / "release" / "datasets" / f"{C.DATASET}.json").write_text(json.dumps({"input_fingerprint": "other"}))
        rc, out = self.go()
        self.assertEqual(rc, 1, out)
        self.assertIn("the input fingerprint is the release report's", out)


class Rules(unittest.TestCase):
    def test_the_arms_run_at_the_defaults_seeds_and_the_export_at_each(self):
        self.assertIn("ARM_SEEDS=${8:-$SEEDS}", SCRIPT)
        self.assertIn('export --species "$SP" --only "$DATASETS" --seed "$s"', function("export_step"))
        self.assertIn('"$EXPORT_GRM" --data', function("export_step"))

    def test_the_plans_are_checked_before_any_fit(self):
        self.assertLess(SCRIPT.index('"$HERE/plans.py" "$OUT/expect.tsv" "$OUT" "$R" "$W/export-s{seed}"'),
                        SCRIPT.index("STEP=fits"))

    def test_the_libraries_are_installed_by_hash_and_each_fit_has_one_thread(self):
        self.assertIn('pip" install --quiet --no-deps --require-hashes -r "$HERE/requirements.txt"', SCRIPT)
        self.assertIn("OMP_NUM_THREADS=1", function("fit"))
        self.assertIn("OPENBLAS_NUM_THREADS=1", function("fit"))

    def test_the_scorers_output_stays_in_score_and_the_end_scans_every_folder_copied_back(self):
        self.assertIn('"$OUT/score" > "$OUT/score/score.log" 2>&1', SCRIPT)
        self.assertIn('scan.py" "$R/species/$SP.json" results score logs', SCRIPT)
        self.assertNotIn("score.log\" >> \"$OUT/logs", SCRIPT)


if __name__ == "__main__":
    unittest.main()
