#!/usr/bin/env python3
"""Every plan entry point runs as the plan writes it: `python3 -I SCRIPT`, as a subprocess, from
another working directory (M1 of the plan-validators review).

usage: FOLDGEN=<foldgen binary> python3 -I scripts/mtb-4k/tests/test_entry_points.py

On Python 3.11 and later `-I` implies `-P`: the script's own directory is not on sys.path, so a
script that imports a sibling module must put that directory back itself. The other test files
insert the path for themselves, which is how write_run_manifest.py, run_revision8_checks.py and
compare_cohorts.py passed their suites while failing as the plan ran them. This file adds no path:

- it reads every fenced block of the plan, finds every Python command (any interpreter name or
  path, an interpreter held in a variable such as `G="python3 $GUARD"`, flag clusters, `-c` and
  `-m`), and requires that set of in-repo scripts to be exactly the reviewed one (a new plan
  command fails here until it is listed), and every command to be exactly `python3 -I`, with
  one exception: the wiring lane's guard, run as `python3 $GUARD` (N3 and N7 of the focused
  re-check);
- it runs the plan's fold-rebuild `python3 -I -c` check as written, beside a planted json.py;
- it runs each with --help from a directory holding planted sibling-named modules (and with
  PYTHONPATH pointing at them): each must import its real siblings, never the planted ones;
- it runs each with real, plan-shaped arguments on fixtures, far enough to reach its own checks;
- a pin that moves stops each of the three report_contract CLIs with exit 2 (I1), not a traceback;
- so does --profile rN while revision_n.py holds its placeholders, naming each, before any
  argument is parsed and with nothing written; revision_n.py itself runs from another directory;
- no file in scripts/mtb-4k shadows a standard-library module, and every script that imports a
  sibling puts its own directory on sys.path first.
"""
import json, os, re, shlex, shutil, subprocess, sys, tempfile, unittest
from pathlib import Path

LANE = Path(__file__).resolve().parents[1]
SCRIPTS = LANE.parent
PLAN = LANE / "4k-execution-plan.md"
# Plan revision 6: the reference diff is diff_reference_prefix.py, and revision_n.py and cleaned_rerun.py gate the run.
ENTRY = {"prepare_execution_roots.py", "fingerprint.py", "pilot_interrupt.py", "pilot_invalidation.py", "validate_matrix.py",
         "screen_duplicates.py", "write_run_manifest.py", "run_revision8_checks.py", "diff_reference_prefix.py", "compare_cohorts.py",
         "metric-audit/audit.py", "revision_n.py", "cleaned_rerun.py", "verify_evidence.py"}
# Plan revision 7: verify_evidence.py runs from the reviewed checkout, so the plan runs no script from elsewhere.
EXTERNAL = set()
SIBLINGS = ("report_contract", "fingerprint", "prepare_execution_roots", "run_revision8_checks", "pilot_interrupt", "revision_n",
            "cleaned_rerun")
# The wiring lane's guard runs as `python3 $GUARD`, read-only and as written; pilot_interrupt.py requires the launcher's
# command line to be exactly that (N1), so it takes no -I.
GUARD = "$GUARD"
# The three CLIs that load report_contract run under the suite's profile (tests/profile_boot.py reads the same
# variable); their argparse default, r8, cannot load on a revision-N tree.
PROFILE = os.environ.get("MTB4K_TEST_PROFILE", "r8")
CONTRACT_CLIS = {"write_run_manifest.py", "run_revision8_checks.py", "compare_cohorts.py"}
FENCE = re.compile(r"```[^\n]*\n(.*?)```", re.S)  # every fenced block, whatever its language
INTERPRETER = re.compile(r"""(?<![^\s;|&({"'`])((?:/[\w.-]+)*/)?(python[0-9.]*)(?=[\s"']|$)""")
ALIAS = re.compile(r"""(?<![\w$])(\w+)=(?:"([^"]*)"|'([^']*)'|([^\s;"']+))""")
COMMENT = re.compile(r"(?:^|\s)#.*$", re.M)


def invocation(prefix, name, rest):
    """One Python command: the interpreter word, its flag letters, and what it runs ("script", "c", "m" or None)."""
    lexer = shlex.shlex(rest, posix=True, punctuation_chars=True)
    lexer.whitespace_split = True
    try:
        tokens = list(lexer)
    except ValueError as exc:
        return {"interpreter": (prefix or "") + name, "flags": None, "kind": None, "target": f"unparsable: {exc}"}
    letters, i = [], 0
    while i < len(tokens) and re.fullmatch(r"-[A-Za-z]+", tokens[i]):
        cluster = tokens[i][1:]
        for j, letter in enumerate(cluster):
            if letter in "cm":
                target = cluster[j + 1:] or (tokens[i + 1] if i + 1 < len(tokens) else None)
                return {"interpreter": (prefix or "") + name, "flags": letters, "kind": letter, "target": target}
            letters.append(letter)
        i += 1
    target = tokens[i] if i < len(tokens) and not re.fullmatch(r"[;&|()<>]+", tokens[i]) else None
    if target is not None:
        target = re.sub(r"\$\{(\w+)\}", r"$\1", target)
    return {"interpreter": (prefix or "") + name, "flags": letters, "kind": "script" if target else None, "target": target}


def parse_plan(text):
    """Every Python command in every fenced block, with interpreter variables (`G="python3 $GUARD"`, `PY=python3`)
    expanded where they are used, continuation lines joined and comments dropped."""
    found, aliases = [], {}
    for block in FENCE.findall(text):
        for line in COMMENT.sub("", block.replace("\\\n", " ")).splitlines():
            for m in list(ALIAS.finditer(line)):
                value = next(v for v in m.groups()[1:] if v is not None)
                if INTERPRETER.match(value):
                    aliases[m.group(1)] = value
                    line = line.replace(m.group(0), " ")
            for name, value in aliases.items():
                line = re.sub(rf"\$\{{{name}\}}|\${name}\b", lambda _: value, line)
            for m in INTERPRETER.finditer(line):
                found.append({**invocation(m.group(1), m.group(2), line[m.end():]), "line": line.strip()})
    return found


def plan_scripts(text=None):
    """{in-repo script relative to scripts/mtb-4k or scripts/}, {external paths}, [problems], [-c invocations]."""
    found, external, problems, dash_c = set(), set(), [], []
    for cmd in parse_plan(PLAN.read_text() if text is None else text):
        kind, target = cmd["kind"], cmd["target"]
        if kind == "script" and target.startswith("$P/"):
            found.add(target[3:])
        elif kind == "script" and target.startswith("$WS/scripts/"):
            found.add(target[len("$WS/scripts/"):])
        elif kind == "script" and target != GUARD:
            external.add(target)
        elif kind == "c":
            dash_c.append(cmd)
        want = [] if (kind, target) == ("script", GUARD) else ["I"]
        if cmd["interpreter"] != "python3" or cmd["flags"] != want or kind not in ("script", "c"):
            problems.append(cmd)
    return found, external, problems, dash_c


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
        found, external, problems, dash_c = plan_scripts()
        self.assertEqual(found, ENTRY)
        self.assertEqual(external, EXTERNAL)
        self.assertEqual(problems, [])
        self.assertEqual(len(dash_c), 2)  # the disk/RAM gate and the fold-rebuild check, both -I (N3)
        for rel in found:
            self.assertTrue(path_of(rel).is_file(), rel)

    def test_the_parser_sees_every_shape_the_recheck_used(self):
        # N7: each line, added to the plan, must change what the parser reports. The control must not.
        plan = PLAN.read_text()
        self.assertEqual(plan_scripts(plan + '\n```bash\npython3 -I "$P/compare_cohorts.py" --help\n```\n')[:3], plan_scripts(plan)[:3])
        shapes = ['python3 -B "$P/compare_cohorts.py" --help', 'python3 -I -B "$P/unreviewed.py"', 'python3 -IB "$P/unreviewed.py"',
                  'PY=python3; $PY -I "$P/unreviewed.py"', 'python3.12 -I "$P/unreviewed.py"', '/usr/bin/python3 -I "$P/unreviewed.py"',
                  'python3 -I "${P}/unreviewed.py"', 'python3 -I "$P"/unreviewed.py', 'python3 -I -m unreviewed',
                  'python3 "$P/compare_cohorts.py" --help', 'python3 -B "$P/prepare_execution_roots.py" --verify-frozen "$RUN"',
                  'python3 -I -u "$P/compare_cohorts.py" --help', 'x=$(python3 -c "import json")',
                  'G2="python3 -B $GUARD"; $G2 --log l --resources r -- true', 'python3 "$P/pilot_invalidation.py" \\\n  --copies c']
        for shape in shapes:
            for fence in ("bash", "sh", ""):
                with self.subTest(shape=shape, fence=fence):
                    found, external, problems, _ = plan_scripts(plan + f"\n```{fence}\n{shape}\n```\n")
                    self.assertTrue(found != ENTRY or external != EXTERNAL or problems, shape)

    def test_the_fold_rebuild_check_runs_isolated_as_written(self):
        # N3: the plan's -c check, beside a planted json.py whose load() reports a clean rebuild.
        [line] = [c["line"] for c in plan_scripts()[3] if "rN-fold-rebuild.json" in c["line"]]
        with tempfile.TemporaryDirectory() as tmp:
            run, cwd = Path(tmp) / "run", Path(tmp) / "planted"
            (run / "audit").mkdir(parents=True)
            cwd.mkdir()
            (cwd / "json.py").write_text('def load(f): return {"records": [{}] * 22, "summary": {"errors": []}}\n')
            audit = run / "audit" / "rN-fold-rebuild.json"
            issue = {"issues": [{"where": "splits[0].test_fold", "message": "differs"}]}
            error = {"file": "x.json", "error": "KeyError('splits')"}
            # M1 (REVIEW-MTB12): an error record or a summary error stops it too, and so does a count other than N
            cases = {"clean": ([{"issues": []}] * 22, [], "22", 0, ["22", "0", "0"]),
                     "test_fold issue": ([issue] * 22, [], "22", 1, ["22", "22", "0"]),
                     "error records": ([error] * 22, [], "22", 1, ["22", "0", "22"]),
                     "one error record": ([error] + [{"issues": []}] * 21, [], "22", 1, ["22", "0", "1"]),
                     "summary error": ([{"issues": []}] * 22, [error], "22", 1, ["22", "0", "1"]),
                     "count is not N": ([{"issues": []}] * 21, [], "22", 1, ["21", "0", "0"])}
            for name, (records, errors, n, want, out) in cases.items():
                with self.subTest(case=name):
                    audit.write_text(json.dumps({"records": records, "summary": {"errors": errors}}))
                    p = subprocess.run(["bash", "-c", line], cwd=cwd, env={**os.environ, "RUN": str(run), "N": n}, capture_output=True, text=True)
                    self.assertEqual(p.returncode, want, p.stderr)
                    self.assertEqual(p.stdout.split(), out)

    def test_each_entry_point_loads_its_own_siblings_under_minus_I_from_another_directory(self):
        for rel in sorted(ENTRY):
            with self.subTest(rel):
                p = self.run_isolated(rel, "--help", *(["--profile", PROFILE] if rel in CONTRACT_CLIS else []))
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
                              "--binary", root / "synthetic-binary", "--log", root / "pipeline.log", "--profile", PROFILE)
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
                              "--summary-out", audit / "prefix-summary.json", "--receipt", audit / "prefix-revision8-receipt.json",
                              "--profile", PROFILE)
        self.assertEqual(p.returncode, 2, p.stderr)
        self.assertIn(frozen, p.stderr)
        self.assertClean(p)
        self.assertEqual(sorted(os.listdir(audit)), [])
        cohorts = [f"--cohort={n}={root}" for n in ("prefix", "added", "full", "matched2k")]
        receipts = [f"--checker-receipt={n}={audit}/{n}-revision8-receipt.json" for n in ("prefix", "added", "full", "matched2k")]
        p = self.run_isolated("compare_cohorts.py", *cohorts, *receipts, "--duplicates", audit / "duplicate-screen.json",
                              "--expected-source-commit", fx.SOURCE, "--expected-foldgen-sha256", fx.foldgen_sha256(),
                              "--out", audit / "compare-cohorts.json", "--profile", PROFILE)
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
        # pilot_interrupt.py --list-runs lists nothing on a root nothing runs against (N2).
        p = self.run_isolated("pilot_interrupt.py", "--list-runs", "--data", base, "--grm", "/nonexistent/grm")
        self.assertEqual((p.returncode, json.loads(p.stdout)), (0, []), p.stderr)
        self.assertClean(p)
        # prepare_execution_roots.py --check-lane decides on a directory that is not a lane.
        (base / "data").mkdir()
        p = self.run_isolated("prepare_execution_roots.py", "--check-lane", base / "data")
        self.assertIn(p.returncode, (0, 1), p.stderr)
        self.assertClean(p)
        # Missing inputs are a STOP line, not a traceback (N13).
        p = self.run_isolated("prepare_execution_roots.py", "--verify-frozen", base / "no-run-root")
        self.assertEqual(p.returncode, 1, p.stderr)
        self.assertIn("STOP: --verify-frozen cannot read", p.stderr)
        self.assertClean(p)
        p = self.run_isolated("screen_duplicates.py", "--prefix-root", base / "no-prefix", "--added-root", base / "no-added",
                              "--out", base / "duplicate-screen.json")
        self.assertEqual(p.returncode, 1, p.stderr)
        self.assertIn("screen_duplicates: STOP: FileNotFoundError", p.stderr)
        self.assertClean(p)
        self.assertFalse((base / "duplicate-screen.json").exists())

    def test_a_moved_pin_stops_each_report_contract_cli_with_exit_2(self):
        with tempfile.TemporaryDirectory() as tmp:
            scripts = Path(tmp) / "scripts"
            shutil.copytree(SCRIPTS / "mtb-4k", scripts / "mtb-4k", ignore=shutil.ignore_patterns("tests", "__pycache__"))
            for rel in ("check-results.py", "result_scope.py", "metric-audit/audit.py", "metric-audit/summarize.py"):
                (scripts / rel).parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(SCRIPTS / rel, scripts / rel)
            clis = ("write_run_manifest.py", "run_revision8_checks.py", "compare_cohorts.py")
            for cli in clis:
                p = subprocess.run([sys.executable, "-I", "-B", str(scripts / "mtb-4k" / cli), "--help", "--profile", PROFILE], cwd=self.cwd,
                                   capture_output=True, text=True)
                self.assertEqual(p.returncode, 0, (cli, p.stderr))
            with (scripts / "check-results.py").open("a") as f:
                f.write("\n# edited\n")
            for cli in clis:
                with self.subTest(cli):
                    p = subprocess.run([sys.executable, "-I", "-B", str(scripts / "mtb-4k" / cli), "--help", "--profile", PROFILE], cwd=self.cwd,
                                       capture_output=True, text=True)
                    self.assertEqual(p.returncode, 2, p.stderr)
                    self.assertIn("STOP: cannot load the validators", p.stderr)
                    self.assertIn("differs from the reviewed pin", p.stderr)
                    self.assertNotIn("Traceback", p.stderr)

    def test_profile_rn_refuses_in_each_report_contract_cli_while_revision_n_is_unfilled(self):
        base = Path(self.tmp.name) / "profile-rn"
        base.mkdir()
        lane = Path(self.tmp.name) / "unfilled" / "scripts" / "mtb-4k"
        shutil.copytree(LANE, lane, ignore=shutil.ignore_patterns("tests", "__pycache__"))
        for rel in ("check-results.py", "result_scope.py", "metric-audit/audit.py", "metric-audit/summarize.py"):
            (lane.parent / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(SCRIPTS / rel, lane.parent / rel)
        text = (lane / "revision_n.py").read_text()
        start = text.index("\nVALUES = {\n") + 1
        end = text.index("\n}\n", start) + 3
        names = re.findall(r'^    "([A-Z0-9_]+)": ', text[start:end], flags=re.M)
        self.assertEqual(len(names), 8)
        block = "VALUES = {\n" + "".join(f'    "{name}": "@@{name}@@",\n' for name in names) + "}\n"
        (lane / "revision_n.py").write_text(text[:start] + block + text[end:])
        text = (lane / "cleaned_rerun.py").read_text()
        start = text.index("\nVALUES = {\n") + 1
        end = text.index("\n}\n", start) + 3
        names = re.findall(r'^    "([A-Z0-9_]+)": ', text[start:end], flags=re.M)
        (lane / "cleaned_rerun.py").write_text(text[:start] + "VALUES = {\n" + "".join(f'    "{n}": "@@{n}@@",\n' for n in names) + "}\n" + text[end:])

        def run_copy(rel, *args):
            return subprocess.run([sys.executable, "-I", "-B", str(lane / rel), *map(str, args)], cwd=self.cwd, env=self.env,
                                  capture_output=True, text=True, timeout=600)
        cases = {"write_run_manifest.py": [base, "--cohort", "prefix", "--source-commit", "1" * 40, "--binary", base / "grm"],
                 "run_revision8_checks.py": [base / "audit.json", "--root", base, "--run-manifest", base / "validation-manifest.json",
                                             "--expected-source-commit", "1" * 40, "--binary", base / "grm", "--backup", base,
                                             "--log", base / "log", "--launch", "2000-01-01T00:00:00Z", "--foldgen", base / "fg",
                                             "--expected-foldgen-sha256", "0" * 64, "--checker-out", base / "checker.json",
                                             "--summary-out", base / "summary.json", "--receipt", base / "receipt.json"],
                 "compare_cohorts.py": ["--cohort", f"prefix={base}", "--duplicates", base / "screen.json", "--expected-source-commit",
                                        "1" * 40, "--expected-foldgen-sha256", "0" * 64, "--out", base / "compare.json"]}
        for cli, args in cases.items():
            for profile in (["--profile", "rN"], ["--profile=rN"]):
                with self.subTest(cli=cli, profile=profile):
                    for extra in ([], ["--help"]):
                        p = run_copy(cli, *args, *profile, *extra)
                        self.assertEqual(p.returncode, 2, p.stderr)
                        self.assertIn("STOP: cannot load the validators: profile rN: unfilled placeholder(s): @@REVISION_N@@", p.stderr)
                        self.assertIn("@@RN_CONVERGENCE_STATUS_RULE@@", p.stderr)
                        self.assertNotIn("usage", p.stdout.lower())
                        self.assertClean(p)
        for cli, args in cases.items():
            with self.subTest(cli=cli, profile="rN-cleaned"):
                p = run_copy(cli, *args, "--profile", "rN-cleaned")
                self.assertEqual(p.returncode, 2, p.stderr)
                self.assertIn("STOP: cannot load the validators: profile rN-cleaned: unfilled placeholder(s): @@REVISION_N@@", p.stderr)
                self.assertIn("@@CLEANED_PREFIX_SHA256@@", p.stderr)
                self.assertClean(p)
        p = run_copy("prepare_execution_roots.py", base / "run", "--cohorts", "cleaned")
        self.assertEqual(p.returncode, 1, p.stderr)
        self.assertIn("refusing: unfilled placeholder(s): @@SPECIES_CHECK_PATH@@", p.stderr)
        self.assertClean(p)
        p = run_copy("cleaned_rerun.py", "status", "--group", "prepare")
        self.assertEqual(p.returncode, 2)
        self.assertIn("cleaned_rerun: STOP: unfilled placeholder(s)", p.stderr)
        self.assertEqual(sorted(os.listdir(base)), [])
        p = self.run_isolated("revision_n.py", "--help")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertClean(p)
        p = run_copy("revision_n.py", "status")
        self.assertEqual(p.returncode, 2)
        self.assertIn("revision_n: STOP: unfilled placeholder(s)", p.stderr)
        self.assertClean(p)

    def test_a_sibling_already_loaded_from_elsewhere_is_a_stop(self):
        # N10: the origin check after the import, with a look-alike already in sys.modules (where sys.path cannot help).
        code = ("import runpy, sys, types; m = types.ModuleType(sys.argv[1]); m.__file__ = '/elsewhere/' + sys.argv[1] + '.py'; "
                "sys.modules[sys.argv[1]] = m; script = sys.argv[2]; sys.argv = [script, '--help']; runpy.run_path(script, run_name='__main__')")
        cases = {"write_run_manifest.py": "report_contract", "run_revision8_checks.py": "report_contract",
                 "compare_cohorts.py": "report_contract", "pilot_invalidation.py": "pilot_interrupt"}
        for script, module in cases.items():
            with self.subTest(script):
                p = subprocess.run([sys.executable, "-I", "-B", "-c", code, module, str(LANE / script)], cwd=self.cwd,
                                   capture_output=True, text=True)
                self.assertNotEqual(p.returncode, 0, p.stdout)
                self.assertIn(f"{module} loaded from /elsewhere/{module}.py", p.stderr)
                self.assertIn("STOP", p.stderr)
                self.assertNotIn("usage", p.stdout.lower())

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
