#!/usr/bin/env python3
"""Every plan entry point runs as the plan writes it: `python3 -I SCRIPT`, as a subprocess, from
another working directory (M1 of the plan-validators review).

usage: FOLDGEN=<foldgen binary> python3 -I scripts/mtb-4k/tests/test_entry_points.py

On Python 3.11 and later `-I` implies `-P`: the script's own directory is not on sys.path, so a
script that imports a sibling module must put that directory back itself. The other test files
insert the path for themselves, which is how write_run_manifest.py, run_revision8_checks.py and
compare_cohorts.py passed their suites while failing as the plan ran them. This file adds no path:

- it reads the plan's bash blocks, finds every `python3 -I` command, and requires that set of
  in-repo scripts to be exactly the reviewed one (a new plan command fails here until it is
  listed), and that no in-repo script is run without -I;
- it runs each with --help from a directory holding planted sibling-named modules (and with
  PYTHONPATH pointing at them): each must import its real siblings, never the planted ones;
- it runs each with real, plan-shaped arguments on fixtures, far enough to reach its own checks;
- a pin that moves stops each of the three report_contract CLIs with exit 2 (I1), not a traceback;
- no file in scripts/mtb-4k shadows a standard-library module, and every script that imports a
  sibling puts its own directory on sys.path first.
"""
import json, os, re, shutil, subprocess, sys, tempfile, unittest
from pathlib import Path

LANE = Path(__file__).resolve().parents[1]
SCRIPTS = LANE.parent
PLAN = LANE / "4k-execution-plan.md"
ENTRY = {"prepare_execution_roots.py", "fingerprint.py", "pilot_interrupt.py", "pilot_invalidation.py", "validate_matrix.py",
         "screen_duplicates.py", "write_run_manifest.py", "run_revision8_checks.py", "diff_revision7_prefix.py", "compare_cohorts.py",
         "metric-audit/audit.py"}
EXTERNAL = {"/root/grm-research/grm-research-mtb-4k-cohort-r2/scripts/verify_evidence.py"}
SIBLINGS = ("report_contract", "fingerprint", "prepare_execution_roots", "run_revision8_checks", "pilot_interrupt")
INVOKE = re.compile(r'python3\s+-I\s+"?([^"\s]+\.py)"?')
UNISOLATED = re.compile(r'python3\s+(?!-I\b|-c\b)"?([^"\s]+\.py)"?')


def bash_blocks(text):
    return re.findall(r"```bash\n(.*?)```", text, re.S)


def plan_scripts():
    """{in-repo script relative to scripts/mtb-4k or scripts/}, {external paths}, [unisolated]."""
    found, external, unisolated = set(), set(), []
    for block in bash_blocks(PLAN.read_text()):
        for path in INVOKE.findall(block):
            if path.startswith("$P/"):
                found.add(path[3:])
            elif path.startswith("$WS/scripts/"):
                found.add(path[len("$WS/scripts/"):])
            else:
                external.add(path)
        unisolated += [p for p in UNISOLATED.findall(block) if p.startswith(("$P/", "$WS/"))]
    return found, external, unisolated


def path_of(rel):
    return LANE / rel if (LANE / rel).is_file() else SCRIPTS / rel


class EntryPoints(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.cwd = Path(cls.tmp.name) / "elsewhere"
        cls.cwd.mkdir()
        for name in SIBLINGS + ("json",):
            (cls.cwd / f"{name}.py").write_text("import sys\nprint('PLANTED " + name + "'); sys.exit(99)\n")
        cls.env = {**os.environ, "PYTHONPATH": str(cls.cwd), "PYTHONSAFEPATH": ""}
        print(f"python {sys.version.split()[0]} ({sys.executable}); -I implies -P: {sys.version_info >= (3, 11)}", file=sys.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_isolated(self, rel, *args, cwd=None):
        return subprocess.run([sys.executable, "-I", "-B", str(path_of(rel)), *map(str, args)], cwd=cwd or self.cwd, env=self.env,
                              capture_output=True, text=True, timeout=600)

    def assertClean(self, p):
        for text in (p.stdout, p.stderr):
            self.assertNotIn("PLANTED", text)
            self.assertNotIn("Traceback", text)
            self.assertNotIn("ModuleNotFoundError", text)

    def test_the_plan_runs_exactly_these_entry_points_and_all_with_minus_I(self):
        found, external, unisolated = plan_scripts()
        self.assertEqual(found, ENTRY)
        self.assertEqual(external, EXTERNAL)
        self.assertEqual(unisolated, [])
        for rel in found:
            self.assertTrue(path_of(rel).is_file(), rel)

    def test_each_entry_point_loads_its_own_siblings_under_minus_I_from_another_directory(self):
        for rel in sorted(ENTRY):
            with self.subTest(rel):
                p = self.run_isolated(rel, "--help")
                self.assertEqual(p.returncode, 0, p.stderr[-2000:])
                self.assertIn("usage", p.stdout.lower())
                self.assertClean(p)

    def test_the_three_report_contract_clis_reach_their_own_checks(self):
        # M1's three: with plan-shaped arguments on a fixture cohort, each loads, runs, and stops
        # at the frozen cohort rule (the fixture's genome list is not the frozen one): exit 2,
        # from report_contract, with nothing written.
        sys.path.insert(0, str(LANE / "tests"))
        import fixtures as fx
        base = Path(self.tmp.name)
        root = fx.build_cohort(base / "prefix", "prefix", fx.make_rules())
        manifest = root / "validation-manifest.json"
        frozen = "cohort genome-list hash differs from frozen rule"
        held = base / "held-manifest.json"
        manifest.rename(held)
        p = self.run_isolated("write_run_manifest.py", root, "--cohort", "prefix", "--source-commit", fx.SOURCE,
                              "--binary", root / "synthetic-binary", "--log", root / "pipeline.log")
        self.assertEqual(p.returncode, 2, p.stderr)
        self.assertIn(frozen, p.stderr)
        self.assertClean(p)
        self.assertFalse(manifest.exists() or (root / "validation-manifest.json.part").exists())
        held.rename(manifest)
        audit = base / "audit"
        audit.mkdir()
        p = self.run_isolated("run_revision8_checks.py", audit / "prefix-audit.json", "--root", root, "--run-manifest", manifest,
                              "--expected-source-commit", fx.SOURCE, "--binary", root / "synthetic-binary", "--backup", root / "backup",
                              "--log", root / "pipeline.log", "--launch", "2000-01-01T00:00:00Z", "--foldgen", fx.foldgen_path(),
                              "--expected-foldgen-sha256", fx.foldgen_sha256(), "--checker-out", audit / "prefix-check-results.json",
                              "--summary-out", audit / "prefix-summary.json", "--receipt", audit / "prefix-revision8-receipt.json")
        self.assertEqual(p.returncode, 2, p.stderr)
        self.assertIn(frozen, p.stderr)
        self.assertClean(p)
        self.assertEqual(sorted(os.listdir(audit)), [])
        cohorts = [f"--cohort={n}={root}" for n in ("prefix", "added", "full", "matched2k")]
        receipts = [f"--checker-receipt={n}={audit}/{n}-revision8-receipt.json" for n in ("prefix", "added", "full", "matched2k")]
        p = self.run_isolated("compare_cohorts.py", *cohorts, *receipts, "--duplicates", audit / "duplicate-screen.json",
                              "--expected-source-commit", fx.SOURCE, "--expected-foldgen-sha256", fx.foldgen_sha256(),
                              "--out", audit / "compare-cohorts.json")
        self.assertEqual(p.returncode, 2, p.stderr)
        self.assertIn(frozen, p.stderr)
        self.assertClean(p)
        self.assertFalse((audit / "compare-cohorts.json").exists())
        # fingerprint.py, which the others import, runs to a full pass on the same root.
        p = self.run_isolated("fingerprint.py", root / "data", "--check", "--expect-reports", "22")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertClean(p)

    def test_the_other_entry_points_run_to_their_own_decisions(self):
        base = Path(self.tmp.name) / "others"
        base.mkdir()
        # pilot_interrupt.py refuses a launcher that is not its sibling job, and signals nothing.
        p = self.run_isolated("pilot_interrupt.py", "--launcher", "1", "--guard-resources", base / "guard.json", "--data", base,
                              "--grm", "/nonexistent/grm", "--folded", "128", "--receipt", base / "pilot-interrupt.json")
        self.assertEqual(p.returncode, 1)
        self.assertIn("not a running sibling job", p.stderr)
        self.assertClean(p)
        # pilot_invalidation.py stops on missing copies.
        p = self.run_isolated("pilot_invalidation.py", "--copies", base / "pilot-copies", "--grm", "/nonexistent/grm", "--width", "128",
                              "--receipt", base / "pilot-invalidation.json")
        self.assertEqual(p.returncode, 1)
        self.assertIn("missing", p.stderr)
        self.assertClean(p)
        # prepare_execution_roots.py --check-lane decides on a directory that is not a lane.
        (base / "data").mkdir()
        p = self.run_isolated("prepare_execution_roots.py", "--check-lane", base / "data")
        self.assertIn(p.returncode, (0, 1), p.stderr)
        self.assertClean(p)

    def test_a_moved_pin_stops_each_report_contract_cli_with_exit_2(self):
        with tempfile.TemporaryDirectory() as tmp:
            scripts = Path(tmp) / "scripts"
            shutil.copytree(SCRIPTS / "mtb-4k", scripts / "mtb-4k", ignore=shutil.ignore_patterns("tests", "__pycache__"))
            for rel in ("check-results.py", "result_scope.py", "metric-audit/audit.py", "metric-audit/summarize.py"):
                (scripts / rel).parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(SCRIPTS / rel, scripts / rel)
            clis = ("write_run_manifest.py", "run_revision8_checks.py", "compare_cohorts.py")
            for cli in clis:
                p = subprocess.run([sys.executable, "-I", "-B", str(scripts / "mtb-4k" / cli), "--help"], cwd=self.cwd, capture_output=True, text=True)
                self.assertEqual(p.returncode, 0, (cli, p.stderr))
            with (scripts / "check-results.py").open("a") as f:
                f.write("\n# edited\n")
            for cli in clis:
                with self.subTest(cli):
                    p = subprocess.run([sys.executable, "-I", "-B", str(scripts / "mtb-4k" / cli), "--help"], cwd=self.cwd, capture_output=True, text=True)
                    self.assertEqual(p.returncode, 2, p.stderr)
                    self.assertIn("STOP: cannot load the validators", p.stderr)
                    self.assertIn("differs from the reviewed pin", p.stderr)
                    self.assertNotIn("Traceback", p.stderr)

    def test_no_lane_file_shadows_the_standard_library(self):
        names = {p.stem for p in LANE.glob("*.py")} | {p.stem for p in (LANE / "tests").glob("*.py")}
        self.assertEqual(names & set(sys.stdlib_module_names), set())

    def test_every_script_that_imports_a_sibling_puts_its_directory_first(self):
        siblings = {p.stem for p in LANE.glob("*.py")}
        for path in sorted(LANE.glob("*.py")):
            text = path.read_text()
            uses = [m.start() for name in siblings - {path.stem}
                    for m in re.finditer(rf'^(?:import {name}\b|from {name} import|.*import_module\("{name}"\)|.*_sibling\("{name}"\))', text, re.M)]
            if not uses:
                continue
            with self.subTest(path.name):
                boot = text.find("sys.path.insert(0,") if "sys.path.insert(0," in text else text.find("sys.path.insert(0, ")
                self.assertGreaterEqual(boot, 0, "no sys.path bootstrap")
                self.assertLess(boot, min(uses))
                self.assertIn("Path(__file__).resolve().parent", text[:boot])


if __name__ == "__main__":
    result = unittest.main(exit=False, verbosity=2).result
    print(json.dumps({"python": sys.version.split()[0], "tests": result.testsRun, "ok": result.wasSuccessful()}))
    sys.exit(0 if result.wasSuccessful() else 1)
