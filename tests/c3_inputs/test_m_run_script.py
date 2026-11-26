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

from tests.c3_inputs import common

expect = common.load("expect")
SCRIPT = (common.KIT / "run_m.sh").read_text()
SPECIES_SCRIPT = (common.KIT / "run_species.sh").read_text()
CODE = "\n".join(line for line in SCRIPT.splitlines() if not line.lstrip().startswith("#"))
D = "synthia_coli__drug_m"


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
        steps = ["step same_build first_build", 'm_checks.py" build', "step build changed_build",
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
        self.assertEqual(len(mine), 10)
        for path, sha in mine.items():
            self.assertEqual(theirs.get(path), sha, path)
        # bcalm2 runs only when GRAPH is not grm; the test genome is pinned inside program_env, which is read.
        self.assertEqual(set(theirs) - set(mine),
                         {"$KIT/bcalm-binaries-v2.2.3-Linux.tar.gz", "$AMR/testdna/test_dna.fa"})

    def test_the_records_file_is_held_to_its_line_in_the_master_list(self):
        self.assertIn('master bvbrc/genome_amr.jsonl "$KIT/genome_amr.jsonl"', CODE)
        self.assertIn('m_list.py" "$ROOT" "$KIT/genome_amr.jsonl"', function("m_list"))
        with tempfile.TemporaryDirectory() as d:
            top = Path(d)
            (top / "genome_amr.jsonl").write_text('{"made": "up"}\n')
            sha = subprocess.run(["sha256sum", str(top / "genome_amr.jsonl")], capture_output=True, text=True).stdout[:64]
            (top / "master").write_text(f"{'0' * 64}  bvbrc/other.jsonl\n{sha}  bvbrc/genome_amr.jsonl\n")
            code = function("master", SPECIES_SCRIPT) + f'master bvbrc/genome_amr.jsonl "{top}/genome_amr.jsonl"'
            self.assertEqual(bash(code, {"M": str(top / "master")})[0], 0)
            (top / "genome_amr.jsonl").write_text('{"made": "up", "changed": 1}\n')
            self.assertNotEqual(bash(code, {"M": str(top / "master")})[0], 0)


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
        self.assertIn('git -C "$REPO" diff --name-only "$FIRST_COMMIT" "$GRM_COMMIT"', function("changed_build"))


class End(unittest.TestCase):
    def test_the_end_counts_genome_ids_then_packs_with_pack_py(self):
        self.assertIn('\' "$R/species/$SP.json" results check score logs\n', SCRIPT)
        self.assertIn('pack.py" tar "$OUT" "$R/species/$SP.json" "$OUT/pack.tar"', CODE)
        self.assertIn("pack $pack\"", CODE)
        self.assertNotIn("work/scores", "\n".join(l for l in CODE.splitlines() if "sha256sum" in l))


if __name__ == "__main__":
    unittest.main()
