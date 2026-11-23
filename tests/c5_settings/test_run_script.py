"""The run script's own rules (scripts/c5-settings/run_species.sh), read as text.

No file is run and no genome, call or result is read: the test reads the script alone.
"""

import re
import unittest
from pathlib import Path

SCRIPT = (Path(__file__).resolve().parents[2] / "scripts" / "c5-settings" / "run_species.sh").read_text()


class ScoreOutput(unittest.TestCase):
    def test_the_scorers_output_stays_in_score_and_out_of_the_logs_read_during_a_run(self):
        # logs/ and run.log are read while the run goes on, by the driver's poll and on every failure; score/ is
        # opened only when every species has run, as scripts/c3-inputs/run_species.sh keeps it.
        call = [line for line in SCRIPT.splitlines() if 'score.py"' in line and not line.lstrip().startswith("#")]
        self.assertEqual(len(call), 1)
        self.assertIn('> "$OUT/score/score.log" 2>&1', call[0])
        self.assertNotIn("logs/score.log", SCRIPT)
        self.assertIsNone(re.search(r"(cat|tail|head)[^\n|;]*score\.log", SCRIPT))


if __name__ == "__main__":
    unittest.main()
