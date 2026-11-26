"""The run script of the run on the datasets of M (scripts/c3-inputs/run_m.sh): its order and its own rules.

The script's text is read, and its own shell functions run in bash on
SYNTHETIC folders: made-up results, reports and score files, no genome. No
step of a run is started.
"""

import csv
import re
import subprocess
import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import common, m_common

expect = common.load("expect")
checks = common.load("m_checks")
SCRIPT = (common.KIT / "run_m.sh").read_text()
SPECIES_SCRIPT = (common.KIT / "run_species.sh").read_text()
CODE = "\n".join(line for line in SCRIPT.splitlines() if not line.lstrip().startswith("#"))
D = "synthia_coli__drug_m"
# The records file of ruling R26: the 1,285,111 lab records the pinned table was rebuilt from on 21 Nov 2026.
RECORDS_SHA256 = "d231fa768a31fc4ce4c8a6bb7024d05d632621846070688adb892fbfe340fd78"


def function(name, text=SCRIPT):
    """A shell function's text, to the first line that is `}` alone. m_missing ends `; }`, so its text runs on
    through keep_reports, which bash reads as two whole functions."""
    m = re.search(rf"(?ms)^{re.escape(name)}\(\) \{{.*?^\}}$", text)
    if m is None:
        raise AssertionError(f"no function {name}()")
    return m.group(0) + "\n"


def bash(code, env):
    r = subprocess.run(["bash", "-euo", "pipefail", "-c", code], env=dict(env, PATH="/usr/bin:/bin"),
                       capture_output=True, text=True)
    return r.returncode, r.stdout, r.stderr


def at(text):
    i = CODE.find(text)
    if i < 0:
        raise AssertionError(f"run_m.sh does not run: {text}")
    return i


class Order(unittest.TestCase):
    def test_the_three_checks_run_before_any_fit_of_m(self):
        steps = ["\nchanged_files\n", "step same_build first_build", 'm_checks.py" build', "step build changed_build",
                 '"$KIT/m-checked.txt" "$CHECK_SPECIES" "$CODE" "$GRM_COMMIT"', "step matrix matrix",
                 "step m_list m_list", 'm_fits "$OUT/expect-first.tsv"', 'm_checks.py" results',
                 'same_counts "$OUT/expect-first.tsv" first', 'touch "$OUT/flags/m_checked.done"',
                 '[ "$SP" != "$CHECK_SPECIES" ] || [ -e "$OUT/flags/m_checked.done" ]', 'm_fits "$OUT/expect.tsv"',
                 'same_counts "$OUT/expect.tsv" m', 'score.py" "$OUT/expect.tsv"', 'second.py" "$OUT/expect.tsv"',
                 "xargs sha256sum > SHA256SUMS", 'pack.py" check', 'pack.py" tar', 'echo "$line" > DONE']
        places = [at(s) for s in steps]
        self.assertEqual(places, sorted(places))
        self.assertEqual(CODE.count('m_fits "$OUT/expect'), 2)

    def test_each_check_names_what_differs_and_stops(self):
        for check in ("Same build differs", "Same results differs", "Same counts differs"):
            line = next(l for l in CODE.splitlines() if check in l)
            self.assertTrue(line.rstrip().endswith("false; }"), line)
        self.assertIn('echo "FAILED $SP $1 $STEP', function("fail"))

    def test_the_check_species_is_m_checks_and_any_other_needs_its_record(self):
        self.assertIn('CHECK_SPECIES=$(python3 -I -B "$HERE/m_checks.py" species "$HERE/lists.tsv"', CODE)
        self.assertIn('[ -f "$KIT/m-checked.txt" ] || { say "$SP is not the check species', CODE)


class Checked(unittest.TestCase):
    def test_another_species_runs_only_on_the_check_species_record_of_its_code(self):
        snippet = re.search(r"(?s)if \[ \"\$SP\" != \"\$CHECK_SPECIES\" \]; then\n  python3 -I -c '(.*?)' \"\$KIT/m-checked.txt\"",
                            SCRIPT).group(1)
        code, commit = "c" * 64, "f" * 40
        with tempfile.TemporaryDirectory() as d:
            rec = Path(d) / "m-checked.txt"
            # As the check species writes it: printf 'species %s\ncode %s\ncommit %s\nchecked %s\n'.
            rec.write_text(f"species neisseria_gonorrhoeae\ncode {code}\ncommit {commit}\nchecked 2026-11-26T00:00:00Z\n")
            run = lambda *a: subprocess.run(["python3", "-I", "-c", snippet, str(rec), *a], capture_output=True).returncode
            self.assertEqual(run("neisseria_gonorrhoeae", code, commit), 0)
            self.assertEqual((run("neisseria_gonorrhoeae", "d" * 64, commit), run("escherichia_coli", code, commit),
                              run("neisseria_gonorrhoeae", code, "e" * 40)), (1, 1, 1))


class Fits(unittest.TestCase):
    def test_the_one_lineage_fit_writes_scores_and_the_default_and_twin_keep_reports(self):
        calls = [l for l in CODE.splitlines() if re.match(r'\s*fit "fit-', l)]
        self.assertEqual(len(calls), 1)
        body = function("m_lineage_fit")
        self.assertIn('--exp-scores "$SCORES/$arm"', body)
        self.assertIn('case $arm in k31 | twin) keep_reports "$arm" "$seed" "$todo" ;; esac', body)
        self.assertNotRegex(CODE, r"(?m)^\s*lineage_fit ")
        self.assertIn('study_fit "$arm" "$d"', function("m_study"))

    def test_the_fits_take_the_first_runs_seeds_twin_arms_and_graph(self):
        self.assertIn("SEEDS=42,43,44; TWIN=45,46,47; ARMS=program,k31-program,k31-counts,k31-aa,unitigs; GRAPH=grm",
                      CODE)

    def test_m_missing_asks_for_the_pick_inputs_score_file_and_kept_report(self):
        rows = expect.rows("synthia_coli", [D], [D], [], [42], [45], ["program"])
        with tempfile.TemporaryDirectory() as d:
            top = Path(d)
            exp = top / "expect.tsv"
            with open(exp, "w") as f:
                f.write("\t".join(common.COLUMNS) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))
            env = {"PICKS": str(top / "picks"), "SCORES": str(top / "scores"), "KEPT": str(top / "reports")}
            files = [f"picks/k31/{D}.lineage-s42.pick.json", f"scores/k31/{D}.lineage-s42.scores.json",
                     f"reports/k31/{D}.lineage-s42.json", f"picks/program/{D}.lineage-s42.pick.json",
                     f"picks/program/{D}.drawn.inputs.json", f"scores/program/{D}.lineage-s42.scores.json",
                     f"picks/twin/{D}.lineage-s45.pick.json", f"scores/twin/{D}.lineage-s45.scores.json",
                     f"reports/twin/{D}.lineage-s45.json"]
            for f in files:
                (top / f).parent.mkdir(parents=True, exist_ok=True)
                (top / f).write_text("{}")
            ask = lambda arm, seed: bash(function("m_missing") + f'm_missing "{exp}" {arm} lineage {seed}', env)[1].strip()
            self.assertEqual((ask("k31", 42), ask("program", 42), ask("twin", 45)), ("", "", ""))
            for gone, arm, seed in ((files[2], "k31", 42), (files[4], "program", 42), (files[5], "program", 42),
                                    (files[8], "twin", 45), (files[0], "k31", 42)):
                (top / gone).rename(top / "aside")
                self.assertEqual(ask(arm, seed), D, gone)
                (top / "aside").rename(top / gone)

    def test_keep_reports_copies_each_report_where_its_seed_writes_it(self):
        with tempfile.TemporaryDirectory() as d:
            top = Path(d)
            env = {"REPORTS": str(top / "root/species/synthia_coli/results"), "ROOT": str(top / "root"),
                   "SP": "synthia_coli", "KEPT": str(top / "kept")}
            for folder, text in (("results", "seed 42"), ("draws/seed-43", "seed 43")):
                for x in (D, "synthia_coli__drug_n"):
                    p = top / "root/species/synthia_coli" / folder / f"{x}.json"
                    p.parent.mkdir(parents=True, exist_ok=True)
                    p.write_text(f"{text} {x}")
            for seed in (42, 43):
                rc, _, err = bash(function("keep_reports") + f'keep_reports twin {seed} "{D},synthia_coli__drug_n"', env)
                self.assertEqual(rc, 0, err)
            got = sorted(p.name for p in (top / "kept" / "twin").iterdir())
            self.assertEqual(got, sorted(f"{x}.lineage-s{s}.json" for x in (D, "synthia_coli__drug_n") for s in (42, 43)))
            self.assertEqual((top / "kept/twin" / f"{D}.lineage-s43.json").read_text(), f"seed 43 {D}")
            rc, _, _ = bash(function("keep_reports") + 'keep_reports k31 44 "synthia_coli__drug_q"', env)
            self.assertNotEqual(rc, 0)


class Shared(unittest.TestCase):
    def test_run_species_functions_are_read_from_it_27_of_27(self):
        block = re.search(r"(?ms)^SHARED=.*?^done\n", SCRIPT).group(0)
        code = 'say() { echo "$*"; }\n' + block + 'echo "imported $(echo $SHARED | wc -w)"\n'
        rc, out, err = bash(code, {"HERE": str(common.KIT)})
        self.assertEqual((rc, out.strip()), (0, "imported 27"), err)
        names = re.search(r'(?m)^SHARED="([^"]*)"$', SCRIPT).group(1).split()
        self.assertEqual([n for n in names if not re.search(rf"(?m)^{n}\(\) \{{", SPECIES_SCRIPT)], [])

    def test_run_species_sh_still_refuses_a_dataset_of_m(self):
        snippet = re.search(r"(?s)# Only datasets that are not M are fitted.*?python3 -I -c '(.*?)' \"\$KIT/marks.tsv\" \"\$DATASETS\"",
                            SPECIES_SCRIPT).group(1)
        with tempfile.TemporaryDirectory() as d:
            marks = Path(d) / "marks.tsv"
            marks.write_text("dataset\tpairs_S+I\tpairs_I+R\nx__a\t0\t0\nx__m\t2\t0\nx__n\t0\t1\n")
            run = lambda ds: subprocess.run(["python3", "-I", "-c", snippet, str(marks), ds], capture_output=True).returncode
            self.assertEqual((run("x__a"), run("x__a,x__m"), run("x__n")), (0, 1, 1))
        own = re.search(r"(?s)# The kit's label count must mark each of them M.*?python3 -I -c '(.*?)' \"\$KIT/marks.tsv\"",
                        SCRIPT).group(1)
        with tempfile.TemporaryDirectory() as d:
            marks = Path(d) / "marks.tsv"
            marks.write_text("dataset\tpairs_S+I\tpairs_I+R\nx__a\t0\t0\nx__m\t2\t0\n")
            run = lambda ds: subprocess.run(["python3", "-I", "-c", own, str(marks), ds], capture_output=True).returncode
            self.assertEqual((run("x__m"), run("x__m,x__a")), (0, 1))

    def test_the_pins_are_run_species_sh_s(self):
        pin = re.compile(r'(?m)^\s*pinned "(\$[A-Z]+[^"]*)" ([0-9a-f]{64})$')
        mine, theirs = dict(pin.findall(SCRIPT)), dict(pin.findall(SPECIES_SCRIPT))
        self.assertEqual(len(mine), 11)
        # The records file is run_m.sh's own pin (ruling R26); run_species.sh reads no records file.
        self.assertEqual(mine.pop("$KIT/genome_amr.jsonl"), RECORDS_SHA256)
        for path, sha in mine.items():
            self.assertEqual(theirs.get(path), sha, path)
        # bcalm2 runs only when GRAPH is not grm; the test genome is pinned inside program_env, which is read.
        self.assertEqual(set(theirs) - set(mine),
                         {"$KIT/bcalm-binaries-v2.2.3-Linux.tar.gz", "$AMR/testdna/test_dna.fa"})

    def test_the_records_file_is_pinned_by_its_own_sha256(self):
        # The 1,285,111 lab records the pinned table was rebuilt from (ruling R26), not the master list's file.
        self.assertIn(f'pinned "$KIT/genome_amr.jsonl" {RECORDS_SHA256}', CODE)
        self.assertNotIn("bvbrc/genome_amr.jsonl", CODE)
        self.assertIn('m_list.py" "$ROOT" "$KIT/genome_amr.jsonl"', function("m_list"))
        line = next(l for l in CODE.splitlines() if l.startswith('pinned "$KIT/genome_amr.jsonl"'))
        with tempfile.TemporaryDirectory() as d:
            kit = Path(d)
            (kit / "genome_amr.jsonl").write_text('{"made": "up"}\n')
            sha = subprocess.run(["sha256sum", str(kit / "genome_amr.jsonl")], capture_output=True, text=True).stdout[:64]
            code = 'say() { echo "$*"; }\n' + function("pinned", SPECIES_SCRIPT)
            self.assertEqual(bash(code + line.replace(RECORDS_SHA256, sha), {"KIT": str(kit)})[0], 0)
            rc, out, _ = bash(code + line, {"KIT": str(kit)})
            self.assertEqual((rc, "sha256 differs from the pinned value" in out), (1, True))

    def test_scrub_reads_an_id_as_the_pack_check_does(self):
        with tempfile.TemporaryDirectory() as d:
            top = Path(d)
            (top / "release" / "species").mkdir(parents=True)
            (top / "release" / "species" / "synthia_coli.json").write_text(
                '{"genomes": [{"id": "90007.3"}, {"id": "90007.4"}]}')
            log = top / "out" / "fit.log"
            log.parent.mkdir()
            log.write_text("wrote 90007.3_scores and 0.90007.4, read 90007.33 and v90007.4\n")
            env = {"W": str(top / "work"), "R": str(top / "release"), "SP": "synthia_coli", "HERE": str(common.KIT)}
            rc, _, err = bash(function("scrub") + f'scrub "{log}"', env)
            self.assertEqual(rc, 0, err)
            self.assertEqual(log.read_text(), "wrote [genome id]_scores and 0.[genome id], read 90007.33 and v90007.4\n")
            self.assertIn("90007.3_scores", (top / "work" / "raw-logs" / "fit.log").read_text())


class Builds(unittest.TestCase):
    """Both builds run with cargo's home at the first run's path and write it to their record; stubs stand in for
    git, tar, cargo and objcopy."""

    STUBS = (
        'git() { :; }\ntar() { :; }\n'
        'cargo() { mkdir -p "$CARGO_TARGET_DIR/release"; echo "home $CARGO_HOME" > "$CARGO_TARGET_DIR/release/grm"; '
        'echo "$CARGO_HOME" >> "$OUT/homes"; }\n'
        'objcopy() { cp "$4" "$5"; }\n')

    def run_build(self, name, export=None):
        export = re.search(r"(?m)^export CARGO_HOME=.*$", SCRIPT).group(0) if export is None else export
        with tempfile.TemporaryDirectory() as d:
            top = Path(d)
            env = {"W": str(top / "w"), "OUT": str(top / "out"), "REPO": str(top), "FIRST_COMMIT": "a" * 40,
                   "GRM_COMMIT": "b" * 40, "GRM": str(top / "w" / "target" / "release" / "grm")}
            (top / "out").mkdir()
            rc, out, err = bash(export + "\n" + self.STUBS + function(name) + name, env)
            record = "binary-first.txt" if name == "first_build" else "binary.txt"
            text = (top / "out" / record).read_text() if (top / "out" / record).exists() else ""
            homes = (top / "out" / "homes").read_text() if (top / "out" / "homes").exists() else ""
            return rc, text, homes, err

    def test_both_builds_use_the_first_runs_cargo_home(self):
        self.assertIn("export CARGO_HOME=/home/user/.cargo\n", SCRIPT)
        for name in ("first_build", "changed_build"):
            rc, _, homes, err = self.run_build(name)
            self.assertEqual((rc, homes), (0, "/home/user/.cargo\n"), err)

    def test_each_builds_record_names_its_cargo_home(self):
        for name in ("first_build", "changed_build"):
            rc, text, _, err = self.run_build(name)
            self.assertEqual(rc, 0, err)
            self.assertIn("\ncargo_home /home/user/.cargo\n", text)
        self.assertIn('[ "$(sed -n \'s/^cargo_home //p\' "$OUT/binary.txt")" = "$CARGO_HOME" ]', CODE)

    def test_same_build_names_the_cargo_home_when_it_differs(self):
        line = next(l for l in CODE.splitlines() if "Same build differs" in l)
        self.assertIn("cargo's home $CARGO_HOME", line)
        call = CODE[at('m_checks.py" build'):at("Same build differs")]
        self.assertIn('"$(sed -n \'s/^cargo_home //p\' "$OUT/binary-first.txt")"', call)


class SameCountsRun(unittest.TestCase):
    """same_counts itself, run in bash on SYNTHETIC fits: a difference ends it with a non-zero exit."""

    def lay_out(self, top, change=None):
        rows = expect.rows(m_common.SPECIES, [m_common.M_SET], [m_common.M_SET], [], [42], [], ["k31-counts"])
        gs = m_common.genomes()
        exp, picks, scores, kept, study = m_common.lay_out(top, rows, lambda e: gs)
        path = Path(top) / "expect.tsv"
        with open(path, "w") as f:
            f.write("\t".join(common.COLUMNS) + "\n" + "".join("\t".join(e[c] for c in common.COLUMNS) + "\n" for e in exp))
        if change:
            change(picks)
        (Path(top) / "out" / "check").mkdir(parents=True)
        env = {"HERE": str(common.KIT), "PICKS": str(picks), "SCORES": str(scores), "KEPT": str(kept),
               "STUDY": str(study), "OUT": str(Path(top) / "out"), "W": str(Path(top) / "work"),
               "R": str(Path(top) / "release"), "SP": m_common.SPECIES}
        code = 'say() { echo "said: $*"; }\n' + function("scrub") + function("same_counts") + f'same_counts "{path}" m\n'
        return bash(code, env)

    def test_same_counts_ends_non_zero_on_a_synthetic_difference(self):
        def change(picks):
            path = picks / "k31-counts" / f"{m_common.M_SET}.lineage-s42.pick.json"
            doc = checks.load(path)
            doc["folds"][1]["outer"][0]["tp"] += 1
            doc["folds"][1]["outer"][0]["fn"] -= 1
            m_common.write_json(path, doc)
        with tempfile.TemporaryDirectory() as d:
            rc, out, err = self.lay_out(d, change)
            self.assertNotEqual(rc, 0, out + err)
            self.assertIn("said: Same counts differs", out)
            self.assertIn("stopped", (Path(d) / "out" / "check" / "same-counts-m.txt").read_text())

    def test_same_counts_passes_on_synthetic_fits_that_agree(self):
        with tempfile.TemporaryDirectory() as d:
            rc, out, err = self.lay_out(d)
            self.assertEqual(rc, 0, out + err)
            self.assertIn("said: Same counts: ", out)


class Commits(unittest.TestCase):
    first = re.search(r"(?m)^FIRST_COMMIT=([0-9a-f]{40})$", SCRIPT).group(1)
    grm = re.search(r"(?m)^GRM_COMMIT=([0-9a-f]{40})$", SCRIPT).group(1)

    def test_the_first_runs_commit_is_every_page_binaries_row(self):
        with open(common.ROOT / "reference" / "input-arms" / "page-binaries.tsv", newline="") as f:
            commits = {r["commit"] for r in csv.DictReader(f, delimiter="\t")}
        self.assertEqual(commits, {self.first})
        self.assertNotEqual(self.grm, self.first)

    def test_the_changed_commit_touches_the_changed_files_alone(self):
        have = lambda c: subprocess.run(["git", "-C", str(common.ROOT), "cat-file", "-e", c + "^{commit}"],
                                        capture_output=True).returncode == 0
        if not (have(self.first) and have(self.grm)):
            self.skipTest("the clone does not hold both commits")
        moved = subprocess.run(["git", "-C", str(common.ROOT), "diff", "--name-only", self.first, self.grm],
                               capture_output=True, text=True, check=True).stdout.split()
        changed = re.search(r'(?m)^CHANGED="([^"]*)"$', SCRIPT).group(1).split()
        self.assertEqual(sorted(moved), sorted(changed))
        self.assertIn('git -C "$REPO" diff --name-only "$FIRST_COMMIT" "$GRM_COMMIT"', function("changed_files"))

    def test_a_diff_touching_a_file_off_the_changed_list_stops(self):
        # changed_files itself, in bash, with a stub git: the order of the diff does not matter, a file more or less
        # than CHANGED ends it with a non-zero exit. It runs where the clone holds neither commit.
        changed = re.search(r'(?m)^CHANGED="([^"]*)"$', SCRIPT).group(1).split()
        env = {"REPO": "/nowhere", "FIRST_COMMIT": self.first, "GRM_COMMIT": self.grm, "CHANGED": " ".join(changed)}
        for diff, stops in ((changed[::-1], False), (changed + ["src/synthetic.rs"], True), (changed[1:], True)):
            stub = "git() { printf '%s\\n' " + " ".join(diff) + "; }\n"
            rc, out, err = bash('say() { echo "said: $*"; }\n' + stub + function("changed_files") + "changed_files\n", env)
            self.assertEqual(rc != 0, stops, out + err)
            self.assertEqual("said: the changed commit's diff from the first run's touches:" in out, stops, out)


class End(unittest.TestCase):
    def test_the_end_counts_genome_ids_then_packs_with_pack_py(self):
        self.assertIn('\' "$R/species/$SP.json" results check score logs\n', SCRIPT)
        self.assertIn('pack.py" tar "$OUT" "$R/species/$SP.json" "$OUT/pack.tar"', CODE)
        self.assertIn("pack $pack\"", CODE)
        self.assertNotIn("work/scores", "\n".join(l for l in CODE.splitlines() if "sha256sum" in l))


if __name__ == "__main__":
    unittest.main()
