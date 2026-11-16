"""The lasso layout check's wiring: the workflow's job runs check.sh with no argument, the
script names the out-of-line `column_sum` of src/learn/lasso.rs and the page's three layouts,
and docs/ci.md and docs/lasso-layout.md say so. It reads files of the tree and builds nothing."""

import pathlib
import re
import unittest

from tests.lasso_layout import common
from tests.lasso_layout.common import CI_DOC, COMMAND, JOB, LASSO, PAGE, ROOT, SCRIPT, SUITE, WORKFLOW, job, jobs, section, steps, text

# The sentence of docs/lasso-layout.md, part Check, on when the check runs.
SENTENCE = "CI runs it on every push and pull request, in the `lasso-layout` job (docs/ci.md); the timing below is run by hand."
# The sentence it replaced, which ran the check by hand after a change to lasso.rs only.
OLD_SENTENCE = "Run it after a change to `lasso.rs`, to the release profile or to the toolchain"
LAYOUTS = ["profile", "one-unit", "sixteen-units"]
CACHE_PATHS = ["~/.cargo/registry/index/", "~/.cargo/registry/cache/", "~/.cargo/git/db/"]


def suite_count():
    loader = unittest.TestLoader()
    return loader.discover(start_dir=str(ROOT / SUITE), top_level_dir=str(ROOT)).countTestCases()


def named_step(job_text, name):
    found = [s for s in steps(job_text) if s[0] == f"      - name: {name}"]
    return found[0] if len(found) == 1 else None


def script_symbol():
    m = re.search(r"^symbol=(\S+)$", text(SCRIPT), re.M)
    return m.group(1) if m else None


def ok_line():
    """The script's line for a layout that passes, as a pattern with the layout as a group."""
    m = re.search(r'results\+=\("(\$layout: ok, [^"]*)"\)', text(SCRIPT))
    if m is None:
        return None
    pattern = ""
    for part in re.split(r"(\$layout|\$symbol|\$\(\(16#\$size\)\)|\$count)", m.group(1)):
        if part == "$layout":
            pattern += r"(?P<layout>[a-z-]+)"
        elif part == "$symbol":
            pattern += re.escape(script_symbol())
        elif part in ("$((16#$size))", "$count"):
            pattern += r"[1-9][0-9]*"
        else:
            pattern += re.escape(part)
    return pattern


class Job(unittest.TestCase):
    def setUp(self):
        self.job = job(JOB)
        self.assertIsNotNone(self.job, f"ci.yml has no single job {JOB}")

    def test_each_job_once(self):
        keys = [name for name, _ in jobs()]
        self.assertEqual(len(keys), len(set(keys)), keys)
        self.assertEqual(keys.count(JOB), 1)
        names = [re.search(r"\n    name: (.+)\n", chunk).group(1) for _, chunk in jobs()]
        self.assertEqual(len(names), len(set(names)), names)

    def test_runner_and_timeout(self):
        self.assertIn("\n    runs-on: ubuntu-24.04\n", self.job)
        m = re.search(r"\n    timeout-minutes: (\d+)\n", self.job)
        self.assertIsNotNone(m, "the job has no timeout")
        minutes = int(m.group(1))
        self.assertTrue(5 <= minutes <= 45 and minutes % 5 == 0, minutes)

    def test_checkout_as_rust_test(self):
        first = steps(self.job)[0]
        self.assertEqual(first, steps(job("rust-test"))[0])
        self.assertRegex(first[0], r"^      - uses: actions/checkout@[0-9a-f]{40} # v")
        self.assertEqual(first[1:], ["        with:", "          persist-credentials: false"])

    def test_install_rust_as_rust_test(self):
        step = named_step(self.job, "Install Rust")
        self.assertIsNotNone(step, "the job has no Install Rust step")
        self.assertEqual(step, named_step(job("rust-test"), "Install Rust"))
        self.assertTrue(any("grep -qw avx2 /proc/cpuinfo" in line for line in step))

    def test_checks_nm_and_objdump_first(self):
        found = [i for i, s in enumerate(steps(self.job)) if any(line.strip() == "for t in nm objdump; do" for line in s)]
        self.assertEqual(len(found), 1, "no single step checks for nm and objdump")
        step = steps(self.job)[found[0]]
        self.assertTrue(any("::error::$t not found" in line and "exit 1" in line for line in step), step)
        run = [i for i, s in enumerate(steps(self.job)) if f"        run: {COMMAND}" in s]
        self.assertEqual(len(run), 1)
        self.assertLess(found[0], run[0])

    def test_caches_registry_and_git_only(self):
        caches = [s for s in steps(self.job) if s[0].startswith("      - uses: actions/cache@")]
        self.assertEqual(len(caches), 1)
        self.assertEqual(caches[0][0], steps(job("rust-test"))[2][0], "not the pinned cache action")
        paths = [line.strip() for line in caches[0] if line.startswith("            ")]
        self.assertEqual(paths, CACHE_PATHS)
        key = [line.strip() for line in caches[0] if line.strip().startswith("key: ")]
        self.assertEqual(len(key), 1)
        self.assertTrue(key[0].startswith("key: cargo-lasso-layout-"), key[0])
        others = "".join(chunk for name, chunk in jobs() if name != JOB)
        self.assertNotIn("cargo-lasso-layout-", others)

    def test_runs_the_script_with_no_argument(self):
        found = [s for s in steps(self.job) if any(line.startswith("        run:") and "scripts/lasso-layout/check.sh" in line for line in s)]
        self.assertEqual(len(found), 1, "no single step runs check.sh")
        step = found[0]
        self.assertEqual(step[-1], f"        run: {COMMAND}")
        self.assertEqual(step, steps(self.job)[-1], "check.sh is not the job's last step")
        work = [line.strip() for line in step if line.strip().startswith("GRM_LAYOUT_WORK:")]
        self.assertEqual(len(work), 1, "GRM_LAYOUT_WORK is not set for the run")
        self.assertRegex(work[0], r"^GRM_LAYOUT_WORK: \$\{\{ runner\.temp \}\}/[a-z-]+$")

    def test_nothing_lets_it_pass_on_failure(self):
        self.assertNotRegex(self.job, r"\n *(if|continue-on-error):")

    def test_every_push_and_pull_request(self):
        on = text(WORKFLOW).split("\non:\n", 1)[1].split("\n\n", 1)[0]
        self.assertEqual(on, "  push:\n  pull_request:")
        self.assertNotRegex(self.job, r"\n *(paths|paths-ignore|branches):")

    def test_the_layouts_are_the_script_s(self):
        head = text(WORKFLOW).split("\njobs:\n", 1)[0]
        for name in ("RUSTFLAGS", "CARGO_PROFILE_", "CARGO_TARGET_DIR", "CARGO_BUILD_TARGET"):
            self.assertNotIn(name, self.job)
            self.assertNotIn(name, head)


class Symbol(unittest.TestCase):
    def test_column_sum_out_of_line(self):
        lines = text(LASSO).split("\n")
        found = [i for i, line in enumerate(lines) if re.match(r"\s*(pub(\([a-z]+\))? )?fn column_sum\(", line)]
        self.assertEqual(len(found), 1, "src/learn/lasso.rs has no single fn column_sum")
        self.assertEqual(lines[found[0] - 1].strip(), "#[inline(never)]")

    def test_script_names_the_symbol(self):
        crate = re.search(r'^\[package\]\nname = "([a-z_]+)"', text(common.CARGO), re.M)
        self.assertIsNotNone(crate)
        self.assertNotIn("\n[lib]\n", text(common.CARGO), "a [lib] name would change the symbol's crate")
        module = LASSO.relative_to(ROOT / "src").with_suffix("").parts
        self.assertEqual(script_symbol(), "::".join((crate.group(1), *module, "column_sum")))

    def test_layouts_are_the_page_s(self):
        script = text(SCRIPT)
        self.assertIn(f"\nlayouts=({' '.join(LAYOUTS)})\n", script)
        self.assertIn('\ndeclare -A units=([profile]="" [one-unit]=1 [sixteen-units]=16)\n', script)
        check = section(text(PAGE), "## Check")
        self.assertIn("with the release profile, with one codegen unit and with 16 codegen units", check)
        self.assertEqual(script.count("cargo build --release --locked --bin grm)"), 2)

    def test_no_argument_times_nothing(self):
        script = text(SCRIPT)
        self.assertRegex(script, r"\ncase \$# in\n    0 \| 3 \| 4\) ;;\n")
        stop = script.index('\nif [ "$#" -eq 0 ] || [ "$failed" -ne 0 ]; then\n    exit "$failed"\nfi\n')
        self.assertLess(stop, script.index(" learn --species "))
        self.assertLess(stop, script.index("\nseconds() {\n"))

    def test_each_result_names_its_layout(self):
        script = text(SCRIPT)
        results = re.findall(r'results\+=\("([^"]*)"', script)
        self.assertEqual(len(results), 4, results)
        self.assertTrue(all(r.startswith("$layout: ") for r in results), results)
        self.assertIn('results+=("$layout: FAIL, $symbol moves a floating-point register to or from the stack:" "$spills")', script)
        printed = script.index('\nprintf \'%s\\n\' "${results[@]}"\n')
        self.assertGreater(printed, script.rindex(" results+=("))
        self.assertLess(printed, script.index('\nif [ "$#" -eq 0 ] || [ "$failed" -ne 0 ]; then\n'))


class Documents(unittest.TestCase):
    def test_ci_doc_row(self):
        jobs_part = section(text(CI_DOC), "## Jobs")
        table = [line for line in jobs_part.split("\n") if line.startswith("| ")]
        rows = [line for line in table if line.startswith(f"| `{JOB}` |")]
        self.assertEqual(len(rows), 1, "the Jobs table has no single lasso-layout row")
        cells = rows[0].split(" | ")
        self.assertEqual(len(cells), 3, cells)
        self.assertTrue(cells[1].startswith(f"`{COMMAND}`"), cells[1])
        pattern = ok_line()
        self.assertIsNotNone(pattern, "check.sh has no line for a layout that passes")
        shown = [m.group("layout") for m in re.finditer(f"`{pattern}`", cells[2])]
        self.assertEqual(shown, LAYOUTS, "the row does not give the result line of each layout")

    def test_ci_doc_keeps_timing_outside_ci(self):
        outside = section(text(CI_DOC), "## Not in CI")
        items = [line for line in outside.split("\n") if line.startswith("- ") and "scripts/lasso-layout/check.sh" in line]
        self.assertEqual(len(items), 1, "Not in CI has no single item for the lasso timing")
        self.assertIn("`scripts/lasso-layout/check.sh DATA_ROOT SPECIES DATASET [ROUNDS]`", items[0])
        self.assertIn(f"`{JOB}`", items[0])

    def test_page_says_ci_runs_it(self):
        page = text(PAGE)
        self.assertEqual(page.count(SENTENCE), 1)
        self.assertIn(SENTENCE, section(page, "## Check"))
        self.assertNotIn(OLD_SENTENCE, page)


class Wiring(unittest.TestCase):
    def test_in_the_python_job(self):
        python = job("python")
        self.assertIsNotNone(python)
        self.assertRegex(python, rf"\n          - suite: {re.escape(SUITE)}\n            top: \.\n")
        self.assertEqual(python.count(f"suite: {SUITE}\n"), 1)

    def test_python_row_gives_the_count(self):
        rows = [line for line in text(CI_DOC).split("\n") if line.startswith("| `python`")]
        self.assertEqual(len(rows), 1)
        m = re.search(rf"`{re.escape(SUITE)}` (\d+)(,| \|)", rows[0])
        self.assertIsNotNone(m, "the python row does not give the suite's count")
        self.assertEqual(int(m.group(1)), suite_count())

    def test_ci_doc_names_the_suite_once_in_its_text(self):
        self.assertEqual(text(CI_DOC).count(f"The `{SUITE}` suite"), 1)

    def test_paths_inside_the_repository(self):
        for name in ("WORKFLOW", "CI_DOC", "PAGE", "SCRIPT", "LASSO", "CARGO"):
            path = getattr(common, name)
            self.assertTrue(path.is_relative_to(ROOT), name)
            self.assertTrue(path.exists(), name)
        for path in sorted(pathlib.Path(ROOT / SUITE).glob("*.py")):
            source = path.read_text(encoding="utf-8")
            self.assertIsNone(re.search(r"[\"'](/(home|Users|mnt|tmp|root|opt|srv)/|[A-Za-z]:\\)", source), path.name)
            self.assertIsNone(re.search(r"[\"']\.\./", source), path.name)
            self.assertIsNone(re.search(r"^ *(import|from) +(subprocess|shutil|urllib)", source, re.M), path.name)


if __name__ == "__main__":
    unittest.main()
