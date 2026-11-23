"""The run script's own rules (scripts/c3-inputs/run_species.sh), on SYNTHETIC files.

No genome, call or result is read: the time records here are written by the
test, or by GNU time over a shell that the test ends with SIGKILL.
"""

import os
import re
import shutil
import subprocess
import tempfile
import unittest

from tests.c3_inputs.common import KIT

SCRIPT = (KIT / "run_species.sh").read_text()


def function(name):
    m = re.search(rf"^{name}\(\) \{{[^\n]*\n.*?^\}}\n", SCRIPT, re.S | re.M)
    if m is None:
        raise AssertionError(f"run_species.sh defines no {name}()")
    return m.group(0)


def killed(record):
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "step.time")
        with open(path, "w") as f:
            f.write(record)
        return subprocess.run(["bash", "-c", function("killed_for_memory") + 'killed_for_memory "$1"', "x", path]).returncode == 0


class MemoryKill(unittest.TestCase):
    def test_a_step_ended_by_sigkill_is_named_as_killed_for_memory(self):
        # GNU time over bash -c "<function>; <name>": the shell exits 137, and time writes that status.
        self.assertTrue(killed("Command exited with non-zero status 137\nprodigal 1.20 s wall, 3412 KB peak\n"))
        # GNU time -v over grm itself, as a fit runs it.
        self.assertTrue(killed("\tCommand terminated by signal 9\n\tMaximum resident set size (kbytes): 9\n"))
        self.assertFalse(killed("Command exited with non-zero status 1\nprodigal 1.20 s wall, 3412 KB peak\n"))
        self.assertFalse(killed("prodigal 1.20 s wall, 3412 KB peak\n"))

    @unittest.skipUnless(shutil.which("time", path="/usr/bin"), "GNU time is not at /usr/bin/time")
    def test_gnu_times_record_of_a_killed_step_is_named_as_killed_for_memory(self):
        # As step() runs a step: /usr/bin/time -o over bash -euo pipefail -c; the step's program ends by SIGKILL.
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "victim.time")
            r = subprocess.run(["/usr/bin/time", "-f", "victim %e s wall, %M KB peak", "-o", path,
                                "bash", "-euo", "pipefail", "-c", 'victim() { sh -c "kill -9 \\$\\$"; echo after; }; victim'],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            self.assertEqual(r.returncode, 137)
            with open(path) as f:
                record = f.read()
        self.assertNotIn("signal 9", record)
        self.assertTrue(killed(record))

    def test_every_step_and_fit_failure_asks_whether_memory_ran_out(self):
        self.assertIn('! killed_for_memory "$OUT/logs/$1.time"', function("step"))
        self.assertIn('! killed_for_memory "$OUT/logs/$name.time"', function("fit"))
        self.assertNotIn("grep -q 'signal 9'", SCRIPT)


class ScoreOutput(unittest.TestCase):
    def test_the_scorers_output_stays_in_score_and_out_of_the_logs_read_during_a_run(self):
        # logs/ and run.log are read while the run goes on and on every failure; score/ is opened only
        # when every species has run.
        call = [line for line in SCRIPT.splitlines() if 'score.py"' in line and not line.lstrip().startswith("#")]
        self.assertEqual(len(call), 1)
        self.assertIn('> "$OUT/score/score.log" 2>&1', call[0])
        self.assertNotIn("logs/score.log", SCRIPT)
        self.assertIsNone(re.search(r"(cat|tail|head)[^\n|;]*score\.log", SCRIPT))


class Builds(unittest.TestCase):
    def test_a_builds_sha256_list_holds_every_file_of_the_build(self):
        # Every start holds each built matrix to this list, so a file it leaves out could change unseen.
        with tempfile.TemporaryDirectory() as d:
            inputs, out = os.path.join(d, "inputs"), os.path.join(d, "out")
            os.makedirs(os.path.join(inputs, "counts", "matrix"))
            os.makedirs(os.path.join(out, "logs"))
            for name, text in (("build.json", '{"kind": "counts", "genomes": 8}'), ("genomes.tsv", "made-up\n"),
                               ("matrix/matrix.json", "{}"), ("matrix/patterns.bin", "made-up patterns")):
                with open(os.path.join(inputs, "counts", name), "w") as f:
                    f.write(text)
            r = subprocess.run(["bash", "-euo", "pipefail", "-c", function("keep_build") + "keep_build counts"],
                               env={"PATH": os.environ.get("PATH", "/usr/bin:/bin"), "INPUTS": inputs, "OUT": out},
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            self.assertEqual(r.returncode, 0)
            with open(os.path.join(out, "logs", "inputs-counts.sha256")) as f:
                listed = sorted(line.split("  ", 1)[1].strip() for line in f)
        self.assertEqual(listed, ["counts/build.json", "counts/genomes.tsv", "counts/matrix/matrix.json",
                                  "counts/matrix/patterns.bin"])

    def test_the_end_counts_genome_ids_in_every_folder_that_is_copied_back(self):
        self.assertIn('\' "$R/species/$SP.json" results check score logs\n', SCRIPT)


if __name__ == "__main__":
    unittest.main()
