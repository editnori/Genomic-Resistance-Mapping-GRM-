"""Plant each fault of docs/outside-learners-result.md's gates and run checks, and run the test meant to catch it.

  python3 -I page_faults.py <out dir> [fault ...]

Run it with the python of the venv that holds requirements.txt's wheels.
Copies this checkout's scripts/ and tests/ to <out dir>/tree, a plain copy,
so it runs the same in a copy with no git history. One fault at a time, it
changes one piece of text there, runs the one test named for it, and puts
the text back. A fault is caught when its test fails and ran; each test must
also pass with no fault planted. The tests run as CI runs them, through
scripts/c7-outside/planted_faults.py's run_test. Writes <out dir>/faults.tsv
(fault, gate, what it plants, file, test, result) and <out dir>/logs/,
prints one line per fault, and exits 1 when a fault is not caught. Every
test it runs is SYNTHETIC; it reads no data.
"""

import importlib.util
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PAGE = "scripts/c7-outside/page_tables.py"
CHECKS = "scripts/c7-outside/run_checks.py"
G = "tests.c7_outside.test_page_tables.Gates"
T = "tests.c7_outside.test_page_tables.Tables"
C = "tests.c7_outside.test_run_checks.RunChecks"
L = "tests.c7_outside.test_run_checks.ScoreLog"

# fault, gate of the brief, what it plants, file, text, planted text, test
FAULTS = [
    ("1", "reading", "a folder outside the reading passed", PAGE,
     '    extra = sorted(p.name for p in top.iterdir() if p.is_dir() and p.name not in ("count", "score", "run"))',
     "    extra = []",
     f"{G}.test_a_folder_outside_the_reading_is_refused"),
    ("2a", "species", "a first copy set aside read as its species", PAGE,
     "        got = sorted(p.name for p in (top / name).iterdir() if p.is_dir())\n"
     "        for folder_name in got:",
     '        got = sorted({p.name.split(".")[0] for p in (top / name).iterdir() if p.is_dir()})\n'
     "        for folder_name in []:",
     f"{G}.test_a_first_copy_set_aside_is_never_read_as_a_species"),
    ("2b", "species", "a species outside the first run passed", PAGE,
     "        if got != species:",
     "        if not set(species) <= set(got):",
     f"{G}.test_a_species_outside_the_first_run_is_refused"),
    ("3", "M", "a dataset of M counted", PAGE,
     "            if ds & in_m:\n"
     '                raise Stop(f"gate 3 M: score/{slug}/{f} holds {len(ds & in_m)} datasets of M")\n'
     "            if ds - first:",
     "            if ds & in_m and not in_m:\n"
     '                raise Stop(f"gate 3 M: score/{slug}/{f} holds {len(ds & in_m)} datasets of M")\n'
     "            if ds - first - in_m:",
     f"{G}.test_a_dataset_of_m_is_never_counted"),
    ("4a", "pairs", "the elastic net read on gbt", PAGE,
     "                    if learner != LEARNER_OF[arm]:",
     "                    if learner not in PAIRED:",
     f"{G}.test_an_elastic_net_read_against_gbt_is_refused"),
    ("4b", "pairs", "an arm read against the default's every draw, not at its own seeds", PAGE,
     "    return DEFAULT if arm == TWIN else ARM_BASE",
     "    return DEFAULT",
     f"{T}.test_an_arm_is_read_against_the_default_at_its_own_seeds"),
    ("5", "twin", "the twin read as an arm in conditions.tsv", PAGE,
     "    if TWIN in cond_arms or sorted(set(cond_arms)) != sorted(ARMS):",
     "    if sorted(set(cond_arms) - {TWIN}) != sorted(ARMS):",
     f"{G}.test_conditions_that_read_the_twin_as_an_arm_are_refused"),
    ("6a", "study", "a study row counted on the lineage split by the gate", PAGE,
     '            if (r["role"] == "study") != (r["split"] == "study"):',
     '            if r["split"] == "lineage" and r["role"] == "study":\n'
     "                continue\n"
     '            if (r["role"] == "study") != (r["split"] == "study"):',
     f"{G}.test_a_study_run_moved_onto_the_lineage_split_is_refused"),
    ("6b", "study", "a study row counted on the lineage split in the species table", PAGE,
     '        groups.setdefault((r["dataset"].split("__")[0], r["split"], r["pair"]), []).append(r)',
     '        groups.setdefault((r["dataset"].split("__")[0], "lineage", r["pair"]), []).append(r)',
     f"{T}.test_a_study_row_stays_on_the_study_split"),
    ("7", "family", "a family short of a pair passed", PAGE,
     '        if keys != sorted(FAMILY) or any(r["holm_p"] == "" for r in fam):',
     "        if not set(keys) <= set(FAMILY):",
     f"{G}.test_a_family_short_of_a_pair_is_refused"),
    ("8a", "withdrawn", "a withdrawn arm's reports let through", PAGE,
     '            folder = top / "run" / slug / "outside" / arm\n'
     "            if folder.is_dir() and any(folder.iterdir()):",
     '            folder = top / "run" / slug / "outside" / arm\n'
     "            if False:",
     f"{G}.test_a_withdrawn_arms_report_is_refused"),
    ("8b", "withdrawn", "a withdrawn arm's pair read off p = 1", PAGE,
     '                                                   or num(r["sign_p"]) != 1.0):',
     "                                                   ):",
     f"{G}.test_a_withdrawn_pair_off_p_one_is_refused"),
    ("9", "passes", "a first pass alone read as a whole run", PAGE,
     '            if "elasticnet" not in data.withdrawn:',
     "            if False:",
     f"{G}.test_a_first_pass_alone_is_refused_unless_the_elastic_net_is_withdrawn"),
    ("10", "rows", "a first-run dataset left out of datasets.tsv", PAGE,
     "        if len(set(keys)) != len(keys) or set(keys) != want:",
     "        if len(set(keys)) != len(keys):",
     f"{G}.test_a_dataset_left_out_of_datasets_tsv_is_refused"),
    ("11a", "missing", "a missing draw's figures let through", PAGE,
     '            if r["scored"] != "1" and (not r["why"] or any(r[m] != "" for m in MEASURES)):',
     '            if r["scored"] != "1" and not r["why"]:',
     f"{G}.test_a_missing_draw_with_figures_is_refused"),
    ("11b", "missing", "a dataset not scored read as zero in the dataset table", PAGE,
     "                row[diff_col] = r12(va - vb) if both and va is not None and vb is not None else None",
     "                row[diff_col] = r12((va or 0) - (vb or 0))",
     f"{T}.test_a_dataset_not_scored_holds_no_difference"),
    ("12", "draws", "a dataset figure that is not its draws' mean", PAGE,
     "                if not same_cell(r[c], C9.fmt(g.get(c))):",
     '                if c == "dataset" and not same_cell(r[c], C9.fmt(g.get(c))):',
     f"{G}.test_a_dataset_figure_that_is_not_its_draws_mean_is_refused"),
    ("13a", "reports", "a report's fold counts not held to its folds", PAGE,
     "            if any(doc.get(k) != v for k, v in sums.items()):",
     "            if False:",
     f"{G}.test_a_report_whose_fold_count_is_not_its_folds_sum_is_refused"),
    ("13b", "reports", "a report with no draw let through", PAGE,
     "        if on_disk != sorted(marks):",
     "        if not set(marks) <= set(on_disk):",
     f"{G}.test_a_report_with_no_draw_is_refused"),
    ("14", "checks", "a failed score check let through", PAGE,
     "        if not checks or failed:",
     "        if not checks:",
     f"{G}.test_a_failed_score_check_is_refused"),
    ("15a", "binary", "library versions that differ written as one row", PAGE,
     "    found = versions(data)\n"
     "    if len(found) == 1:",
     "    found = versions(data)\n"
     "    if True:",
     f"{T}.test_versions_that_differ_are_named_with_each_species"),
    ("15b", "binary", "a DONE's result count not held to the reports", PAGE,
     '        if int(d["results"]) != len(data.reports[slug]):',
     "        if False:",
     f"{G}.test_a_done_with_another_result_count_is_refused"),
    ("16a", "beside", "the no-constant and no-unconverged counts written as deciding", PAGE,
     '            row = {"file": name, "decides": int(name == DECIDING_COUNT), **r,',
     '            row = {"file": name, "decides": 1, **r,',
     f"{T}.test_the_count_rows_decide_for_count_tsv_alone"),
    ("16b", "beside", "the arm against its own learner's twin pair alone written as deciding", PAGE,
     '                    "condition": r["condition"], "decides": int(r["decides"] == "1"),',
     '                    "condition": r["condition"], "decides": 1,',
     f"{T}.test_the_beside_condition_decides_nothing"),
    ("17", "count", "count files score.py did not write let through", PAGE,
     "        if data.count[name] != rows_want:",
     "        if False:",
     f"{G}.test_a_count_figure_score_py_did_not_write_is_refused"),
    ("R1", "run: species", "a first copy newer than its DONE set aside", CHECKS,
     '        if m["slug"] in species and done is not None and iso(m["time"]) < done["time"]:',
     '        if m["slug"] in species:',
     f"{C}.test_a_first_copy_newer_than_its_done_fails"),
    ("R2", "run: SHA256SUMS", "a file that differs from its SHA256SUMS line let through", CHECKS,
     "            if data is None or sha(data) != want:",
     "            if data is None:",
     f"{C}.test_a_changed_result_fails_its_sha256sums"),
    ("R3", "run: DONE", "a first pass alone beside second passes let through", CHECKS,
     '    passes_ok = set(by["elasticnet"]) == {want_pass}',
     "    passes_ok = True",
     f"{C}.test_a_first_pass_beside_second_passes_fails_done"),
    ("R4", "run: DONE", "a FAILED newer than its DONE set aside", CHECKS,
     '            if f is not None and f["species"] == s and m is not None and f["time"] < m["time"]:',
     '            if f is not None and f["species"] == s:',
     f"{C}.test_a_failed_newer_than_its_done_fails"),
    ("R5", "run: expect.py", "an expect.tsv that is not expect.py's let through", CHECKS,
     "        ok_species += same_lists and same_expect and on_disk == owed",
     "        ok_species += same_lists and on_disk == owed",
     f"{C}.test_an_expect_tsv_that_is_not_expect_pys_fails"),
    ("R6", "run: genome ids", "a genome id in a log let through", CHECKS,
     '    rep.check("genome ids", hits == 0 and lists_found == len(species) and len(ids) > 0,',
     '    rep.check("genome ids", lists_found == len(species) and len(ids) > 0,',
     f"{C}.test_a_genome_id_in_a_log_fails"),
    ("R7", "run: features", "a matrix column name in a log let through", CHECKS,
     '    rep.check("features", sum(features) == 0,',
     '    rep.check("features", True,',
     f"{C}.test_a_column_name_in_a_log_fails"),
    ("R9", "run: species", "a second-pass copy read when the elastic net is withdrawn after some second passes ran",
     CHECKS, '    late = "elasticnet" in withdrawn', "    late = False",
     f"{C}.test_late_withdrawal_reads_each_first_copy_and_never_opens_a_second_pass"),
    ("R10", "reading", "the count of species that ran a second pass left off the page", PAGE,
     '    data.second_passes = int(found[0]["n"]) if found else 0', "    data.second_passes = 0",
     f"{G}.test_the_second_pass_count_goes_to_the_page"),
    ("R8", "run: score.log", "score/score.log read with score.py's tables", CHECKS,
     "            if name in NOT_READ or name.startswith(NOT_READ_PREFIX):",
     "            if name.startswith(NOT_READ_PREFIX):",
     f"{L}.test_score_log_is_neither_read_nor_copied"),
]


def _planted():
    spec = importlib.util.spec_from_file_location("c7_faults_for_c7_page", HERE / "planted_faults.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(argv):
    if not argv:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    chosen = argv[1:]
    unknown = sorted(set(chosen) - {f[0] for f in FAULTS})
    if unknown:
        raise SystemExit(f"no fault {', '.join(unknown)}")
    run_test = _planted().run_test
    out = Path(argv[0]).resolve()
    (out / "logs").mkdir(parents=True, exist_ok=True)
    tree = out / "tree"
    shutil.rmtree(tree, ignore_errors=True)
    tree.mkdir()
    for top in ("scripts", "tests"):
        shutil.copytree(REPO / top, tree / top, ignore=shutil.ignore_patterns("__pycache__"))
    env = {"HOME": "/tmp"}
    rows, bad = [], 0
    for n, gate, what, path, old, new, test in FAULTS:
        if chosen and n not in chosen:
            continue
        target = tree / path
        text = target.read_text()
        if text.count(old) != 1:
            raise SystemExit(f"fault {n}: the text to change occurs {text.count(old)} times in {path}")
        passed, ran = run_test(tree, "py", test, out / "logs" / f"{n}-clean.log", env)
        if not (passed and ran):
            result = "its test fails with no fault"
        else:
            target.write_text(text.replace(old, new))
            try:
                passed, ran = run_test(tree, "py", test, out / "logs" / f"{n}.log", env)
            finally:
                target.write_text(text)
            result = "caught" if not passed and ran else ("did not run" if not ran else "not caught")
        bad += result != "caught"
        rows.append((n, gate, what, path, test, result))
        print(f"{n}: {result} by {test}", flush=True)
    with open(out / "faults.tsv", "w") as f:
        f.write("fault\tgate\tplants\tfile\ttest\tresult\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"caught {sum(r[5] == 'caught' for r in rows)} of {len(rows)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
