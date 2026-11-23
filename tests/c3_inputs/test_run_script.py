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


if __name__ == "__main__":
    unittest.main()
