"""The committed copy of stage 1's export and the pins of scripts/matched-accuracy/stage2.py: the values ruling
R59 and the brief give, and the proof that no argument reaches them.

The release run here is SYNTHETIC (made_up2.py); its row lists are made up, so the real pins stop it."""

import ast
import contextlib
import hashlib
import io
import subprocess
import sys
import unittest

from tests.matched_accuracy import stage2_common as c
from tests.matched_accuracy.stage2_common import stage2

# The sha256 of each file of stage 1's export, as the lead's brief names them.
STAGE1_EXPORT_SHA256 = {
    "checks.tsv": "cde56a8503c24852ead865b4a4b5bc3008265e3eb985e68ef7bd38276463fade",
    "counts.tsv": "d686f07a3d3d34985b15a7be86e61c13aeac6c026fb925eba25c53f1496ced3e",
    "inputs.tsv": "55b6e4d945841872d6ef82e2c98de85d43d2b31f3e025d1c8c48035f34f8929b",
    "shares.tsv": "6b285b4f37123f412fb70752e49d2b2c56d39ebc6b73663a3e268b1819f9e921",
    "typing.tsv": "d25d1fedfbda087f40f91d4ca026705fa8d1d00808de309e06280629f90982fb",
}
# The sha256 stage 1 printed for each row list in its one real run (ruling R59's Record line).
R59_ROW_LISTS = {"A": "b135a8ccd24f3c9b87113498c58bb76181fc8c143135260df01eafbd8ad7dd26",
                 "B": "e3558b1b5a7a4a9f8ff23c79aca502c25586fcd6f889e4dac8d788cdd0b45d66"}
OPTIONS = ["--reports", "--rows-a", "--rows-b", "--hc20", "--out"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Committed(unittest.TestCase):
    def test_the_five_files_of_stage_1s_export_are_the_bytes_the_brief_names(self):
        self.assertEqual(sorted(p.name for p in c.STAGE1_COPY.iterdir()), sorted(STAGE1_EXPORT_SHA256))
        for name, digest in STAGE1_EXPORT_SHA256.items():
            self.assertEqual(sha(c.STAGE1_COPY / name), digest, name)

    def test_stage_2_reads_the_count_and_typing_tables_from_the_committed_copy(self):
        self.assertEqual(stage2.STAGE1_EXPORT, c.STAGE1_COPY)
        self.assertEqual(stage2.COUNTS, c.STAGE1_COPY / "counts.tsv")
        self.assertEqual(stage2.COUNTS_SHA256, STAGE1_EXPORT_SHA256["counts.tsv"])
        self.assertEqual(stage2.TYPING_TABLE, c.STAGE1_COPY / "typing.tsv")
        self.assertEqual(stage2.TYPING_TABLE_SHA256, STAGE1_EXPORT_SHA256["typing.tsv"])

    def test_both_committed_tables_pass_their_checks(self):
        self.assertEqual(len(stage2.check_count_table()), 596)
        typing = stage2.check_typing_table()
        self.assertEqual(len(typing), 1168)
        self.assertEqual({k[3] for k in typing}, {"mlst", "mlst_hc20"})
        self.assertEqual({k[2] for k in typing}, {"A", "B"})


class Pins(unittest.TestCase):
    def test_the_pins_are_the_values_of_the_brief_and_ruling_r59(self):
        self.assertEqual(stage2.ROW_LIST_SHA256, R59_ROW_LISTS)
        self.assertEqual(stage2.stage1.REPORT_SUMS, c.ROOT / "reference" / "lineage-control" /
                         "rev14-274-20261112.reports.sha256")
        self.assertEqual(sha(stage2.stage1.REPORT_SUMS), stage2.stage1.REPORT_SUMS_SHA256)
        self.assertEqual(stage2.NUMPY_VERSION, "1" + ".26" + ".4")
        self.assertEqual((stage2.SEED, stage2.DRAWS, stage2.MIN_SPECIES, stage2.MIN_PAIRS, stage2.PPM),
                         (42, 10000, 6, 10, 1000000))
        self.assertEqual(stage2.PERCENTILES, (2.5, 97.5))
        self.assertEqual(stage2.FIVE, ("scm", "cart", "rf", "gbt", "lasso"))
        self.assertEqual(stage2.OPENED, "scores opened: yes")


class Route(unittest.TestCase):
    """The made-up pins of the tests are set as module constants in the test's own process. The command line takes
    the five options alone, none of them a pin, and the script reads no environment variable."""

    def test_the_command_line_takes_the_five_options_alone(self):
        tree = ast.parse(c.SCRIPT.read_text(encoding="utf-8"))
        parse = next(f for f in tree.body if isinstance(f, ast.FunctionDef) and f.name == "parse")
        options = [n.args[0].value for n in ast.walk(parse) if isinstance(n, ast.Call)
                   and isinstance(n.func, ast.Attribute) and n.func.attr == "add_argument"]
        self.assertEqual(options, OPTIONS)
        text = c.SCRIPT.read_text(encoding="utf-8")
        for word in ("environ", "getenv", "argv[", "sys.argv"):
            self.assertEqual(text.count(word), 1 if word in ("argv[", "sys.argv") else 0, word)
        self.assertIn("sys.exit(main(sys.argv[1:]))", text)

    def test_an_option_naming_a_pin_is_refused(self):
        for flag in ("--row-list-sha256", "--counts", "--typing-table", "--seed"):
            with self.subTest(flag=flag):
                with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
                    stage2.parse(["--reports", "r", "--rows-a", "a", "--rows-b", "b", "--hc20", "h", "0" * 64,
                                  "--out", "o", flag, "x"])

    def test_made_up_row_lists_through_the_command_line_stop_at_the_real_pins(self):
        run = c.Run()
        self.addCleanup(run.close)
        self.assertEqual(run.code, 0, run.lines)
        argv = list(run.argv)
        out = run.out.parent / "by-command-line"
        argv[argv.index("--out") + 1] = str(out)
        done = subprocess.run([sys.executable, "-I", str(c.SCRIPT)] + argv, capture_output=True, text=True,
                              timeout=300, check=False)
        self.assertEqual(done.returncode, 1, done.stdout + done.stderr)
        self.assertEqual(done.stdout, c.STOP + "the row list of reading A differs from the sha256 stage 1 printed\n")
        self.assertEqual(done.stderr, "")
        self.assertFalse(out.exists())


if __name__ == "__main__":
    unittest.main()
