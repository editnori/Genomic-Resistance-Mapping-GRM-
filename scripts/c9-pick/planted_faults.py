"""Plant each fault of docs/served-pick.md's table and run the test meant to catch it.

  python3 -I planted_faults.py <out dir> [fault ...]

Copies this checkout's last commit to <out dir>/tree (git archive) and, one
fault at a time (every fault, or those named), changes one piece of text there, runs the one test named
for it, and puts the text back. A fault is caught when its test fails and
the tree still builds; each test must also pass with no fault planted.
Rust tests always build in <out dir>/target, whatever CARGO_TARGET_DIR
says: the copy has the same crate and relative paths as the checkout, so in
a shared target folder a planted build left there looks fresh to cargo in
the checkout, whose files are older. Writes
<out dir>/faults.tsv (fault, what it plants, file, test, result) and
<out dir>/logs/, prints one line per fault, and exits 1 when a fault is not
caught. Every test it runs is SYNTHETIC; it reads no data.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]

FAULTS = [
    ("1", "an outer test genome in an inner fold", "src/learn/inner_pick.rs",
     "train.iter_ones().filter(|&i| inner[i] != j)",
     "(0..n).filter(|&i| inner[i] != j)",
     "rust", "learn::inner_pick::tests::no_outer_test_genome_reaches_an_inner_fit"),
    ("2", "an inner split that ignores lineage groups", "src/learn/inner_pick.rs",
     "        Some(lineage),\n        k,",
     "        None,\n        k,",
     "rust", "learn::inner_pick::tests::inner_folds_keep_each_lineage_whole_and_are_the_release_tuning_partition"),
    ("3a", "inner folds that ignore the run's seed", "src/learn/inner_pick.rs",
     "crate::ablate::plan::fold_seed(seed, f),",
     "crate::ablate::plan::fold_seed(42, f),",
     "rust", "learn::inner_pick::tests::the_seed_keys_every_inner_fit_and_partition"),
    ("3b", "inner fits that ignore the run's seed", "src/learn/inner_pick.rs",
     'format!("c9-inner/{seed}/{fold}/{j}")',
     'format!("c9-inner/42/{fold}/{j}")',
     "rust", "learn::inner_pick::tests::the_seed_keys_every_inner_fit_and_partition"),
    ("4a", "a tie sent to the earlier learner under the later rule", "src/report.rs",
     "PickTie::Later => *tied.last()?,",
     "PickTie::Later => *tied.first()?,",
     "rust", "learn::inner_pick::tests::a_tie_goes_where_the_rule_says_and_is_listed"),
    ("4b", "the scorer's re-pick sending a tie the wrong way", "scripts/c9-pick/score.py",
     '    if tie == "later":\n        return tied[-1], tied',
     '    if tie == "later":\n        return tied[0], tied',
     "py", "tests.c9_pick.test_score.Served.test_a_tie_is_picked_again_by_the_scorers_rule"),
    ("5", "the default pick read from the lineage split", "src/report.rs",
     "    let rows = learner_rows(splits);",
     "    let mut rows = learner_rows(splits);\n    for l in rows.iter_mut() {\n        l.random = l.lineage;\n    }",
     "rust", "learn::inner_pick::tests::the_default_pick_reads_the_random_split_as_the_release_does"),
    ("6", "the picked learner scored on its inner folds", "scripts/c9-pick/score.py",
     '        counts["picks_" + learner] += 1\n'
     '        c = add(c, f["outer"][LEARNERS.index(learner)])',
     '        counts["picks_" + learner] += 1\n'
     '        inner = {l["learner"]: l for l in f["inner"]["learners"]}\n'
     '        k4 = ("tp", "fp", "tn", "fn")\n'
     '        c = add(c, dict(n=sum(inner[learner][k] for k in k4), loss_sum=0.0, **{k: inner[learner][k] for k in k4})\n'
     '                if learner in inner else f["outer"][LEARNERS.index(learner)])',
     "py", "tests.c9_pick.test_score.Served.test_the_served_figure_is_each_folds_pick_on_its_outer_test_fold"),
    ("7", "a twin row read as the inner pick's", "scripts/c9-pick/page_tables.py",
     'inner = cells[(d, split, "inner")]',
     'inner = cells.get((d, split, "twin"), cells[(d, split, "inner")])',
     "py", "tests.c9_pick.test_page_tables.Tables.test_the_inner_figure_is_the_inner_picks_and_the_twin_stands_apart"),
    ("8", "a beside answer read as the read one", "scripts/c9-pick/page_tables.py",
     'READ = "read"',
     'READ = "tie-earlier"',
     "py", "tests.c9_pick.test_page_tables.Tables.test_the_read_answer_is_the_reading_and_the_beside_answers_stand_apart"),
    ("8b", "a read folder scored under a beside answer's options passed", "scripts/c9-pick/page_tables.py",
     'if answers[a][slug]["options"] != options_text(opts):',
     'if False:',
     "py", "tests.c9_pick.test_page_tables.Gates.test_an_answer_scored_under_other_options_is_refused"),
    ("9", "a species dropped", "scripts/c9-pick/page_tables.py",
     "        if got != species:",
     "        if not set(got) <= set(species):",
     "py", "tests.c9_pick.test_page_tables.Gates.test_a_species_missing_from_the_reading_is_refused"),
    ("10", "a dataset of M counted", "scripts/c9-pick/page_tables.py",
     'return [r for r in lrows if r["first_run"] == "1"]',
     'return [r for r in lrows if r["first_run"] in ("0", "1")]',
     "py", "tests.c9_pick.test_page_tables.Gates.test_a_dataset_of_m_is_refused"),
    ("11", "a study row counted on the lineage split", "scripts/c9-pick/page_tables.py",
     '    mine = [r for r in rows if r["split"] == split]',
     '    mine = [r for r in rows if r["split"] != "random"]',
     "py", "tests.c9_pick.test_page_tables.Tables.test_a_study_row_stays_on_the_study_split"),
    ("12", "a missing draw read as zero", "scripts/c9-pick/page_tables.py",
     '    return None if x in ("", None) else float(x)',
     '    return 0.0 if x in ("", None) else float(x)',
     "py", "tests.c9_pick.test_page_tables.Tables.test_a_draw_not_scored_leaves_its_dataset_unscored_not_zero"),
]


def run_test(tree, kind, name, log, env):
    if kind == "rust":
        cmd = ["cargo", "test", "--release", "--locked", "--lib", "--", "--exact", name]
    else:
        # As CI runs the suite: discover from the top, so `tests` imports under -P.
        suite = "/".join(name.split(".")[:2])
        cmd = ["python3", "-s", "-P", "-B", "-W", "error", "-m", "unittest", "discover", "-s", suite, "-t", ".", "-k", name]
        env = {"PATH": "/usr/bin:/bin", "HOME": env.get("HOME", "/tmp"), "PYTHONHASHSEED": "0"}
    with open(log, "w") as f:
        p = subprocess.run(cmd, cwd=tree, env=env, stdout=f, stderr=subprocess.STDOUT)
    text = Path(log).read_text(errors="replace")
    built = not any(line.startswith("error[E") or line.startswith("error: could not compile") for line in text.splitlines())
    ran = (" 1 passed" in text or " 0 passed" in text) if kind == "rust" else "Ran 1 test" in text
    return p.returncode == 0, built and ran


def main(argv):
    if not argv:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    chosen = argv[1:]
    unknown = sorted(set(chosen) - {f[0] for f in FAULTS})
    if unknown:
        raise SystemExit(f"no fault {', '.join(unknown)}")
    out = Path(argv[0]).resolve()
    (out / "logs").mkdir(parents=True, exist_ok=True)
    tree = out / "tree"
    shutil.rmtree(tree, ignore_errors=True)
    tree.mkdir()
    archive = subprocess.run(["git", "-C", str(REPO), "archive", "HEAD"], check=True, capture_output=True).stdout
    subprocess.run(["tar", "-x", "-C", str(tree)], input=archive, check=True)
    env = dict(os.environ)
    env["CARGO_TARGET_DIR"] = str(out / "target")
    rows, bad = [], 0
    for n, what, path, old, new, kind, test in FAULTS:
        if chosen and n not in chosen:
            continue
        target = tree / path
        text = target.read_text()
        if text.count(old) != 1:
            raise SystemExit(f"fault {n}: the text to change occurs {text.count(old)} times in {path}")
        passed, ok = run_test(tree, kind, test, out / "logs" / f"{n}-clean.log", env)
        if not (passed and ok):
            result = "its test fails with no fault"
        else:
            target.write_text(text.replace(old, new))
            try:
                passed, ok = run_test(tree, kind, test, out / "logs" / f"{n}.log", env)
            finally:
                target.write_text(text)
            result = "caught" if not passed and ok else ("does not build" if not ok else "not caught")
        bad += result != "caught"
        rows.append((n, what, path, test, result))
        print(f"{n}: {result} by {test}", flush=True)
    with open(out / "faults.tsv", "w") as f:
        f.write("fault\tplants\tfile\ttest\tresult\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"caught {sum(r[4] == 'caught' for r in rows)} of {len(rows)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
