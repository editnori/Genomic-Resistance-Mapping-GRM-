"""Plant each fault of docs/cap-arms.md's table and run the test meant to catch it.

  python3 -I planted_faults.py <out dir>

Copies this checkout's last commit to <out dir>/tree (git archive) and, one
fault at a time, changes one piece of text there, runs the one test named
for it, and puts the text back, through scripts/c9-pick/planted_faults.py's
runner, as scripts/c3-inputs/planted_faults.py does. A fault is caught when
its test fails; each test must also pass with no fault planted. Writes
<out dir>/faults.tsv (fault, what it plants, file, test, result) and
<out dir>/logs/, prints one line per fault, and exits 1 when a fault is not
caught. Every test it runs is SYNTHETIC; it reads no data.
"""

import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FOLDS = "scripts/c2-cap/check_folds.py"
STUDY = "scripts/c2-cap/study.py"
COHORTS = "scripts/c2-cap/cohorts.py"
SCORE = "scripts/c2-cap/score.py"
LISTS = "scripts/c2-cap/lists_check.py"
EXPECT = "scripts/c2-cap/expect.py"
FETCH = "scripts/c2-cap/bioprojects.py"
T = "tests.c2_cap"

FAULTS = [
    ("1", "an added genome of a held-out lineage in training, passed by the fold check", FOLDS,
     'bad.add("an added genome of a held-out lineage in training", sum(lin[g] in lineages for g in t_added))',
     'bad.add("an added genome of a held-out lineage in training", 0)',
     "py", f"{T}.test_folds.Lineage.test_an_added_genome_of_a_held_out_lineage_in_training_breaches"),
    ("2a", "an added genome of a held-out BioProject put in a study fold file's training", STUDY,
     "        elif bp[g] in test_bp:",
     "        elif False:",
     "py", f"{T}.test_folds.Study.test_an_added_genome_of_a_held_out_bioproject_stays_out_of_training"),
    ("2b", "an added genome of a held-out BioProject in training, passed by the fold check", FOLDS,
     "                    sum(g in bp and bp[g] in test_bp for g in train_added))",
     "                    0)",
     "py", f"{T}.test_folds.Study.test_the_study_check_breaches_on_a_held_out_bioproject_in_training"),
    ("3a", "an outer test fold that differs between arms, passed by the fold check", FOLDS,
     'bad.add("an outer test fold that differs between arms", int(s != test))',
     'bad.add("an outer test fold that differs between arms", 0)',
     "py", f"{T}.test_folds.Lineage.test_an_outer_test_fold_that_differs_between_arms_breaches"),
    ("3b", "an outer test fold that differs between arms, passed by the scorer", SCORE,
     "ok = base is not None and all(same_scored(got[l], base[l]) for l in LEARNERS)",
     "ok = base is not None",
     "py", f"{T}.test_score.Checks.test_an_arm_scored_on_other_genomes_than_the_defaults_outer_test_folds_fails"),
    ("4", "an added genome in an outer test fold, passed by the fold check", FOLDS,
     'bad.add("an added genome in an outer test fold", len(s & added))',
     'bad.add("an added genome in an outer test fold", 0)',
     "py", f"{T}.test_folds.Lineage.test_an_added_genome_in_an_outer_test_fold_breaches"),
    ("5", "a genome drawn that the rebuilt revision 14 table does not hold", COHORTS,
     "stray = [g for g in ids if g not in ranked]",
     "stray = []",
     "py", f"{T}.test_cohorts.Gate.test_a_genome_the_rebuilt_table_does_not_hold_is_refused"),
    ("6", "a cohort past 4,096", COHORTS,
     "        if cap > GATE:",
     "        if cap > 2 * GATE:",
     "py", f"{T}.test_cohorts.Gate.test_a_cohort_past_4096_is_refused"),
    ("7", "a wider cohort that drops a genome of the 2,000", COHORTS,
     "dropped = base - set(ids)",
     "dropped = set()",
     "py", f"{T}.test_cohorts.Gate.test_a_wider_cohort_that_drops_a_default_genome_is_refused"),
    ("8", "the twin scored as an arm: at a draw's slot or from a wider arm's fit", SCORE,
     'ok = e["fit"] == DEFAULT and int(e["slot"]) not in draw_slots',
     "ok = True",
     "py", f"{T}.test_score.Checks.test_a_twin_at_a_draw_slot_or_from_a_wider_fit_stops_the_scorer"),
    ("9", "a missing draw read as a zero", SCORE,
     '    checks.append((where, "the plan and the prediction rows exist", have, ""))\n'
     "    if not have:\n        return None",
     "    if not have:\n        return {learner: {} for learner in LEARNERS}",
     "py", f"{T}.test_score.Checks.test_a_missing_draw_stops_the_scorer_and_is_never_a_zero"),
    ("10a", "a dataset of M counted in the first run by the runner's list check", LISTS,
     '"M": sum(marked.get(d, True) for d in listed)}',
     '"M": 0}',
     "py", f"{T}.test_lists.ListsCheck.test_a_dataset_of_m_is_not_counted_in_the_first_run"),
    ("10b", "a dataset of M counted in the first run by the scorer", SCORE,
     'if r["first_run"] == "1" and r["in_M"] == "0"}',
     'if r["first_run"] == "1" or r["in_M"] == "1"}',
     "py", f"{T}.test_score.Checks.test_a_dataset_of_m_or_outside_the_first_run_is_refused"),
    ("11a", "the BioProject fetch asks another field of the genome table", FETCH,
     "body = f\"in(genome_id,({','.join(ids)}))&select(genome_id,bioproject_accession)&limit({len(ids)})\"",
     "body = f\"in(genome_id,({','.join(ids)}))&select(genome_id,bioproject_accession,collection_year)&limit({len(ids)})\"",
     "py", f"{T}.test_bioprojects.Fetch.test_the_fetch_asks_the_genome_table_for_genome_id_and_bioproject_accession_only"),
    ("11b", "the BioProject fetch asks another table", FETCH,
     'return f"{API}/genome/", body',
     'return f"{API}/genome_amr/", body',
     "py", f"{T}.test_bioprojects.Fetch.test_the_fetch_asks_the_genome_table_for_genome_id_and_bioproject_accession_only"),
    ("12", "a wider arm with no added genome labelled fitted and counted on the lineage split", EXPECT,
     'role = "none-added" if (a, d) in empty else "same-cohort" if a in same else "draw"',
     'role = "same-cohort" if a in same else "draw"',
     "py", f"{T}.test_lists.Expect.test_an_arm_with_no_added_genome_labelled_is_none_added_on_both_splits"),
    ("13", "a withdrawn arm's conditions read as a fitted arm's", SCORE,
     'gone = r["arm"] in withdrawn',
     "gone = False",
     "py", f"{T}.test_score.Count.test_a_withdrawn_arm_keeps_its_five_pairs_at_p_one_and_the_tables_name_it_withdrawn"),
    ("14", "Enterococcus faecium's one wider cohort scored under both wider arms", SCORE,
     'if e["role"] in ("none-added", "same-cohort"):',
     'if e["role"] == "none-added":',
     "py", f"{T}.test_score.Checks.test_a_none_added_or_same_cohort_lineage_row_is_not_scored"),
]


def _c9():
    spec = importlib.util.spec_from_file_location("c9_faults_for_c2", HERE.parent / "c9-pick" / "planted_faults.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(argv):
    if len(argv) != 1:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    c9 = _c9()
    c9.FAULTS = FAULTS
    c9.REPO = HERE.parents[1]
    return c9.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
