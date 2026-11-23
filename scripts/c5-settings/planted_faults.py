"""Plant each fault of docs/settings-arms.md's table and run the test meant to catch it.

  python3 -I planted_faults.py <out dir>

Block C9's runner (scripts/c9-pick/planted_faults.py) with C5's table, run as
block C3's runner runs it (scripts/c3-inputs/planted_faults.py, whose test
runner it takes for the binary's argument tests): it copies this checkout's
last commit to <out dir>/tree, plants one fault at a time, runs the one test
named for it, and puts the text back. A fault is caught when its test fails
and the tree still builds; each test must also pass with no fault planted. A
fault that changes two places in one file is planted as the whole file with
both changes, each change's text checked to occur once. Writes
<out dir>/faults.tsv and <out dir>/logs/, prints one line per fault, and exits
1 when a fault is not caught. Every test it runs is SYNTHETIC; it reads no
data.
"""

import importlib.util
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DRIVE = "scripts/c5-settings/drive.py"
SCORE = "scripts/c5-settings/score.py"
SETTINGS = "src/learn/settings.rs"
TD = "tests.c5_settings.test_drive.Drive"
TS = "tests.c5_settings.test_score"


def both(path, *changes):
    """(old, new) for a fault that makes several changes in one file: the file as committed, and with each
    change made. Each change's text must occur once."""
    text = subprocess.run(["git", "-C", str(REPO), "show", f"HEAD:{path}"], check=True, capture_output=True,
                          text=True).stdout
    new = text
    for old, to in changes:
        if new.count(old) != 1:
            raise SystemExit(f"{path}: the text to change occurs {new.count(old)} times")
        new = new.replace(old, to)
    return text, new


TIE = both(DRIVE, (
    '    if c["picked"] != 0:\n'
    '        trials.append({"id": c["picked"], "learners": {learner: setting(learner, c["picked"])}})',
    '    for t in (c["tied"] if c["picked"] != 0 else []):\n'
    '        trials.append({"id": t, "learners": {learner: setting(learner, t)}})'), (
    '            search = outer_calls(outer_docs[l], c["picked"], l, HALF)',
    '            search = max((outer_calls(outer_docs[l], t, l, HALF) for t in c["tied"]),\n'
    '                         key=lambda x: balanced(x) or 0)'))

FAULTS = [
    ("1", "a setting picked on the outer test fold: the inner run tests it as its third inner fold", DRIVE,
     '            "test": [[ids[i] for i in train if doc["inner"][fold][i] == j] for j in range(INNER_FOLDS)],\n'
     '            "train": [],\n'
     '            "unused": [ids[i] for i in range(len(ids)) if i not in at]}',
     '            "test": [[ids[i] for i in train if doc["inner"][fold][i] == j] for j in range(INNER_FOLDS - 1)]\n'
     '            + [[ids[i] for i in range(len(ids)) if doc["outer"][i] == fold]],\n'
     '            "train": [ids[i] for i in train if doc["inner"][fold][i] == INNER_FOLDS - 1],\n'
     '            "unused": [ids[i] for i in range(len(ids)) if i not in at and doc["outer"][i] != fold]}',
     "py", f"{TD}.test_no_outer_test_genome_is_read_before_the_pick"),
    ("2a", "an inner fold that splits a lineage cluster: inner folds drawn by genome", SETTINGS,
     "training_folds(y, &train, Some(lineage), inner_folds, seed_f)",
     "training_folds(y, &train, None, inner_folds, seed_f)",
     "rust", "learn::settings::tests::inner_folds_keep_lineages_whole_inside_the_outer_training_genomes"),
    ("2b", "the driver taking inner folds that split a lineage cluster", DRIVE,
     "        if lineage_fold.setdefault(doc[\"lineage\"][i], inner[i]) != inner[i]:",
     "        if lineage_fold.setdefault(doc[\"lineage\"][i], inner[i]) != inner[i] and False:",
     "py", f"{TD}.test_a_folds_file_that_splits_a_lineage_is_refused"),
    ("3", "the screen fitted on the fold's test genomes too", "src/learn/eval.rs",
     "    let feats = screen(x, y, train, cfg.screened_features);",
     "    let feats = screen(x, y, &Bits::ones(x.cols), cfg.screened_features);",
     "rust", "learn::settings::tests::no_test_genomes_label_reaches_its_folds_fit"),
    ("4", "the cutoff chosen on outer scores", DRIVE,
     '"cutoff": outer_calls(outer_docs[l], 0, l, c["cutoff"])},',
     '"cutoff": outer_calls(outer_docs[l], 0, l, choose_cutoff(grid_calls(outer_docs[l][0], l))[0])},',
     "py", f"{TD}.test_the_cutoff_is_chosen_on_inner_scores_alone"),
    ("5", "a tie broken by the outer score: every tied setting refitted, the best on the outer fold served",
     DRIVE, TIE[0], TIE[1], "py", f"{TD}.test_a_tie_goes_to_the_lowest_trial_by_inner_scores_alone"),
    ("6a", "the default refitted with a searched setting", DRIVE,
     '    trials = [{"id": 0, "learners": {learner: "release"}, "cutoffs": {learner: c["cutoff"]}}]',
     '    trials = [{"id": 0, "learners": {learner: setting(learner, c["picked"])}, "cutoffs": {learner: c["cutoff"]}}]',
     "py", f"{TD}.test_the_default_arm_is_the_defaults_own_fit"),
    ("6b", "the scorer taking a default that is not the default's own fit", SCORE,
     "same_calls(r[\"outer\"][\"default\"], want))",
     "True)",
     "py", f"{TS}.Score.test_a_choices_file_whose_default_is_not_the_defaults_fit_fails"),
    ("7a", "the twin scored as an arm: the search and the cutoff listed at the twin's seeds",
     "scripts/c5-settings/expect.py",
     '        if arm == DEFAULT and role in ("draw", "study"):',
     '        if arm in (DEFAULT, TWIN):',
     "py", "tests.c5_settings.test_lists.Expect.test_no_arm_runs_at_a_twin_seed"),
    ("7b", "the scorer taking an arm at a twin seed", SCORE,
     '        if e["role"] == "twin" or int(e["seed"]) in twin_seeds:',
     '        if e["role"] == "twin":',
     "py", f"{TS}.Score.test_an_arm_at_a_twin_seed_fails"),
    ("8", "a missing draw read as a zero", SCORE,
     '"scored": 0, "constant": 0,\n             "why": "missing draw"} for l in LEARNERS]',
     '"scored": 1, "constant": 0,\n             "why": "", **{m: 0.0 for m in C3.MEASURES}} for l in LEARNERS]',
     "py", f"{TS}.Score.test_a_missing_draw_is_not_scored_and_never_read_as_zero"),
    ("9", "Holm over the two arms alone, the twin left out of the family", SCORE,
     "module.DEFAULT, module.TWIN, module.FAMILY = DEFAULT, TWIN, list(FAMILY)",
     "module.DEFAULT, module.TWIN, module.FAMILY = DEFAULT, TWIN, list(ARMS)",
     "py", f"{TS}.Count.test_holm_runs_over_the_fifteen_lineage_rows_of_the_primary_measure"),
    ("10", "a trial fitted at another seed than the default's fold", SETTINGS,
     "    let seed = crate::ablate::plan::fold_seed(cfg.seed, fold);",
     "    let seed = crate::ablate::plan::fold_seed(cfg.seed, fold + 1);",
     "rust", "learn::settings::tests::the_release_setting_on_one_lineage_fold_reproduces_the_drawn_fold"),
    ("11", "--exp-settings without --force, so a kept report could skip it", "src/main.rs",
     'requires_all = ["fold_file", "force"]',
     'requires_all = ["fold_file"]',
     "rust-bin", "tests::learn_takes_the_settings_flags_only_with_force_and_alone"),
]


def main(argv):
    spec = importlib.util.spec_from_file_location("c3_faults_for_c5", HERE.parent / "c3-inputs" / "planted_faults.py")
    c3 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(c3)
    # C3's main hands its module's FAULTS to C9's runner, with its own test runner for the binary's tests.
    c3.FAULTS = FAULTS
    return c3.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
