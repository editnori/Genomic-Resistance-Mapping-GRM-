"""The run script's own rules (scripts/c4-screen/run_species.sh), read as text.

No genome, call or result is read: the tests read the scripts, and one evaluates the C4 script's sed of
C3's functions in bash, which defines them and runs none.
"""

import re
import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = (ROOT / "scripts" / "c4-screen" / "run_species.sh").read_text()
C3_SCRIPT = (ROOT / "scripts" / "c3-inputs" / "run_species.sh").read_text()
DEFINED = re.compile(r"^([a-z_][a-z0-9_]*)\(\) \{", re.M)


def taken():
    """The names of C3's functions that the C4 script reads with its sed and holds in SHARED."""
    shared = re.search(r'^SHARED="([^"]*)"', SCRIPT, re.M).group(1).split()
    sed = re.search(r"^eval \"\$\(sed .*?run_species\.sh\"\)\"$", SCRIPT, re.S | re.M).group(0)
    names = set()
    for group in re.findall(r"/\^(?:\\\()?([a-z0-9_\\|]+?)(?:\\\))?\(\) \{", sed):
        names.update(group.split("\\|"))
    return shared, names


def body(name):
    """A function of C3's script as the C4 script's sed reads it: to a line "}" when its first line ends at the
    brace (or a comment after it), else to the first line that ends "; }"."""
    lines = C3_SCRIPT[re.search(rf"^{name}\(\) \{{", C3_SCRIPT, re.M).start():].splitlines()
    if re.match(r"^\S+\(\) \{\s*(#.*)?$", lines[0]):
        end = lines.index("}")
    else:
        end = next(i for i, line in enumerate(lines) if line.rstrip().endswith("; }"))
    return "\n".join(lines[: end + 1])


class ScoreOutput(unittest.TestCase):
    def test_the_scorers_output_stays_in_score_and_out_of_the_logs_read_during_a_run(self):
        # logs/ and run.log are read while the run goes on, by the driver's poll and on every failure; score/ is
        # opened only when every species has run, as scripts/c3-inputs/run_species.sh keeps it.
        call = [line for line in SCRIPT.splitlines() if 'score.py"' in line and not line.lstrip().startswith("#")]
        self.assertEqual(len(call), 1)
        self.assertIn('> "$OUT/score/score.log" 2>&1', call[0])
        self.assertNotIn("logs/score.log", SCRIPT)
        self.assertIsNone(re.search(r"(cat|tail|head)[^\n|;]*score\.log", SCRIPT))


class SharedSteps(unittest.TestCase):
    def test_the_sed_and_shared_name_the_same_functions_of_c3s_script(self):
        shared, names = taken()
        self.assertEqual(sorted(shared), sorted(names))
        self.assertTrue(set(names) <= set(DEFINED.findall(C3_SCRIPT)), set(names) - set(DEFINED.findall(C3_SCRIPT)))

    def test_every_c3_function_a_taken_function_calls_is_taken_or_the_c4_scripts_own(self):
        # step() and fit() call killed_for_memory on a failure and step() hands keep_build to its shell: a function
        # left out is "command not found" there, and the advice to restart with a smaller JOBS is lost.
        _, names = taken()
        c3, own = set(DEFINED.findall(C3_SCRIPT)), set(DEFINED.findall(SCRIPT))
        missing = sorted(
            f"{n} calls {other}"
            for n in names
            for other in c3 - names - own
            if re.search(rf"(?<![\w-]){other}(?![\w-])", body(n))
        )
        self.assertEqual(missing, [])

    def test_bash_defines_every_shared_function_from_the_sed(self):
        lines = re.search(r'^SHARED="[^"]*"\n^eval .*?run_species\.sh"\)"$', SCRIPT, re.S | re.M).group(0)
        prog = f'C3="{ROOT / "scripts" / "c3-inputs"}"\n{lines}\n' \
            'for f in $SHARED; do declare -F "$f" > /dev/null || echo "ABSENT $f"; done\n'
        out = subprocess.run(["bash", "-c", prog], capture_output=True, text=True, stdin=subprocess.DEVNULL)
        self.assertEqual((out.returncode, out.stdout, out.stderr), (0, "", ""))


if __name__ == "__main__":
    unittest.main()
