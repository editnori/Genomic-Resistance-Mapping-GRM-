"""Plant each fault of the result page's table script and run the test meant to catch it.

  python3 -I page_faults.py <out dir> [<fault> ...]

As scripts/c2-cap/planted_faults.py does for the run's scripts, through
scripts/c9-pick/planted_faults.py's runner: copies this checkout's last
commit to <out dir>/tree (git archive) and, one fault at a time, changes one
piece of text in scripts/c2-cap/page_tables.py, runs the one test named for
it, and puts the text back. A fault is caught when its test fails; each test
must also pass with no fault planted. Writes <out dir>/faults.tsv and
<out dir>/logs/, prints one line per fault, and exits 1 when a fault is not
caught. Fault 0 turns off gate 1's hold on run_checks.py's printout, and
faults 1 to 13 each turn off one gate of page_tables.GATES (2 to 14); 14 to
18 break a rule of the page's tables. Every test it runs is
SYNTHETIC (tests/c2_cap/page_common.py); it reads no data.
"""

import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAGE = "scripts/c2-cap/page_tables.py"
G = "tests.c2_cap.test_page_tables.Gates"
T = "tests.c2_cap.test_page_tables.Tables"

FAULTS = [
    ("0", "gate 1 inputs: a run checks printout with a FAIL line, passed", PAGE,
     '    if len(lines) != len(RUN_CHECKS) or passed != list(RUN_CHECKS) or any(line.startswith("FAIL") for line in lines):',
     "    if False:",
     "py", f"{G}.test_a_run_checks_printout_with_a_fail_line_is_refused"),
    ("1", "gate 2 species: a species folder missing from the runs, passed", PAGE,
     "    if found != species:",
     "    if False:",
     "py", f"{G}.test_a_species_folder_missing_from_the_runs_is_refused"),
    ("2", "gate 3 first run: a dataset of M in a species' scored tables, passed", PAGE,
     '                if r["dataset"] not in by_species[slug]:',
     "                if False:",
     "py", f"{G}.test_a_dataset_of_m_in_a_species_scored_tables_is_refused"),
    ("3", "gate 4 cohorts: a cohort past 4,096 genomes, passed", PAGE,
     '            if int(a["genomes"]) > CAP_GATE:',
     '            if int(a["genomes"]) > 2 * CAP_GATE:',
     "py", f"{G}.test_a_cohort_json_past_4096_genomes_is_refused"),
    ("4", "gate 5 binary: a binary.txt with no code sha256, passed", PAGE,
     '        if not run["binary"].get("commit") or not run["binary"].get("code"):',
     '        if not run["binary"].get("commit"):',
     "py", f"{G}.test_a_binary_txt_with_no_code_sha256_is_refused"),
    ("5", "gate 6 none-added: a none-added dataset scored under cap-all, passed", PAGE,
     '                if (r["dataset"], r["split"], r["arm"]) in gone:',
     "                if False:",
     "py", f"{G}.test_a_none_added_dataset_scored_under_cap_all_is_refused"),
    ("6", "gate 7 withdrawn: a condition of cap-4000 read as of an arm that ran, passed", PAGE,
     '    if not cond_mid or any(r["withdrawn"] != "1" or r["met"] != "" or r["learners"] != "" for r in cond_mid):',
     "    if not cond_mid:",
     "py", f"{G}.test_a_withdrawn_arms_condition_read_as_an_arm_that_ran_is_refused"),
    ("7", "gate 8 twin: a twin row from cap-all's fit, passed", PAGE,
     '(r["fit"] != DEFAULT or r["role"] != "twin" or r["slot"] in draw_slots',
     '(r["role"] != "twin" or r["slot"] in draw_slots',
     "py", f"{G}.test_a_twin_row_from_a_wider_arms_fit_is_refused"),
    ("8", "gate 9 study: a study draw put on the lineage split, passed", PAGE,
     '            if (r["split"] == "study") != (r["role"] == "study") or (r["split"] == "study") != (r["slot"] == "-"):',
     "            if False:",
     "py", f"{G}.test_a_study_draw_put_on_the_lineage_split_is_refused"),
    ("9", "gate 10 unavailable: a lineage split unavailable on one arm and scored on another, passed", PAGE,
     "        for d in gone:\n            lin = ",
     "        for d in []:\n            lin = ",
     "py", f"{G}.test_an_unavailable_lineage_split_scored_on_one_arm_is_refused"),
    ("10", "gate 11 missing draw: a draw with four of its five learners, passed", PAGE,
     '            if sorted(have.get((e["dataset"], e["split"], e["arm"], e["slot"]), [])) != sorted(LEARNERS):',
     '            if not have.get((e["dataset"], e["split"], e["arm"], e["slot"])):',
     "py", f"{G}.test_a_draw_missing_from_the_scorers_tables_is_refused"),
    ("11", "gate 12 scored genomes: a cap-all draw with no check that its scored genomes are the default's, passed", PAGE,
     "            if (where, check) not in passed:",
     "            if False:",
     "py", f"{G}.test_a_cap_all_row_with_no_check_of_its_scored_genomes_is_refused"),
    ("12", "gate 13 few added: the count without the datasets that add fewer than 20 labelled genomes, passed", PAGE,
     "        if same_rows(got_count, beside, C3.COUNT_COLUMNS):",
     "        if False:",
     "py", f"{G}.test_a_count_taken_without_the_few_added_datasets_is_refused"),
    ("13", "gate 14 count: a count cell that is not the scorer's, passed", PAGE,
     "    if not same_rows(got_count, want_count, C3.COUNT_COLUMNS):",
     "    if False:",
     "py", f"{G}.test_a_count_cell_that_is_not_the_scorers_is_refused"),
    ("14", "a missing figure read as zero", PAGE,
     '    return None if x in ("", None) else float(x)',
     '    return 0.0 if x in ("", None) else float(x)',
     "py", f"{T}.test_an_unavailable_lineage_split_is_named_with_its_reason_and_holds_no_figure"),
    ("15", "condition 5, which reads the study split, read as deciding", PAGE,
     '    if r["condition"] == "5":\n        return "decides nothing: it reads the study split"',
     '    if False:\n        return "decides nothing: it reads the study split"',
     "py", f"{T}.test_the_study_split_decides_nothing_and_condition_5_is_not_read"),
    ("16", "the count without the few added datasets written as the count", PAGE,
     '    write_tsv(out / "page-count.tsv", page_count(count_rows), COUNT_PAGE_COLUMNS)',
     '    write_tsv(out / "page-count.tsv", page_count(beside), COUNT_PAGE_COLUMNS)',
     "py", f"{T}.test_the_few_added_datasets_stay_in_the_count_and_the_count_beside_decides_nothing"),
    ("17", "cap-4000 read for conditions 1 and 2 on the page", PAGE,
     '        if r["split"] == "lineage" and r["measure"] == PRIMARY and r["withdrawn"] != "1":',
     '        if r["split"] == "lineage" and r["measure"] == PRIMARY:',
     "py", f"{T}.test_cap_4000_is_named_withdrawn_and_holds_no_figure"),
    ("18", "species on two binaries named as one binary", PAGE,
     '        groups.setdefault((b["commit"], b["code"]), [])',
     '        groups.setdefault((b["commit"], runs[min(runs)]["binary"]["code"]), [])',
     "py", f"{T}.test_species_on_two_binaries_are_each_named_with_their_binary"),
]


def _c9():
    spec = importlib.util.spec_from_file_location("c9_faults_for_c2_page", HERE.parent / "c9-pick" / "planted_faults.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(argv):
    if not argv:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    c9 = _c9()
    c9.FAULTS = FAULTS
    c9.REPO = HERE.parents[1]
    return c9.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
