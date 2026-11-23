"""The run script's own rules (scripts/c7-outside/run_species.sh), on SYNTHETIC files.

No genome, call or result is read: the expect lists, species files, logs and the stand-in for the runner are
written by the test.
"""

import json
import re
import stat
import subprocess
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


class Rules(unittest.TestCase):
    def test_the_arms_run_at_seed_42_alone_until_grm_export_takes_a_seed(self):
        self.assertIn('[ "$ARM_SEEDS" = 42 ] || { say "grm export writes seed 42\'s plans only (question 1)', SCRIPT)

    def test_stage_one_is_c3s_run_with_the_unitig_arm_and_grms_graph(self):
        self.assertRegex(SCRIPT, r'bash "\$REPO/scripts/c3-inputs/run_species.sh" "\$KIT" "\$SP" "\$DATASETS" "\$C3" '
                                 r'"\$JOBS" "\$SEEDS" "\$TWIN" unitigs "" grm')
        self.assertIn('if [ -e "$C3/DONE" ]; then', SCRIPT)

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
