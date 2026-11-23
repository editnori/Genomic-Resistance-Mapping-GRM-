"""Plant each fault of docs/input-arms-result.md's tables and run checks, and run the test meant to catch it.

  python3 -I page_faults.py <out dir> [fault ...]

Copies this checkout's scripts/ and tests/ to <out dir>/tree, a plain copy,
so it runs the same in a copy with no git history. One fault at a time, it
changes one piece of text there, runs the one test named for it, and puts
the text back. A fault is caught when its test fails and ran; each test must
also pass with no fault planted. The tests run as CI runs them, through
scripts/c9-pick/planted_faults.py's run_test. Writes <out dir>/faults.tsv
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
PAGE = "scripts/c3-inputs/page_tables.py"
CHECKS = "scripts/c3-inputs/run_checks.py"
T = "tests.c3_inputs.test_page_tables"
C = "tests.c3_inputs.test_run_checks.RunChecks"

# fault, gate of the brief, what it plants, file, text, planted text, test
FAULTS = [
    ("1", "reading", "a folder outside the reading passed", PAGE,
     '    extra = sorted(p.name for p in top.iterdir() if p.is_dir() and p.name not in ("count", "score", "run"))',
     "    extra = []",
     f"{T}.Gates.test_a_folder_outside_the_reading_is_refused"),
    ("1b", "reading", "a reading written though the run's checks did not all pass", PAGE,
     "    if (not tail or tail.group(1) != tail.group(2) or int(tail.group(2)) != len(RUN_CHECKS)\n"
     "            or passed != list(RUN_CHECKS) or any(line.startswith(\"FAIL\") for line in lines)):",
     "    if not tail:",
     f"{T}.Gates.test_a_checks_printout_with_a_failed_check_is_refused"),
    ("2", "species", "a species read twice", PAGE,
     "        if got != species:",
     "        if not set(species) <= set(got):",
     f"{T}.Gates.test_a_species_twice_is_refused"),
    ("3", "twin", "the twin read as an arm in conditions.tsv", PAGE,
     "    if TWIN in cond_arms or sorted(set(cond_arms)) != sorted(ARMS):",
     "    if sorted(set(cond_arms) - {TWIN}) != sorted(ARMS):",
     f"{T}.Gates.test_conditions_that_read_the_twin_as_an_arm_are_refused"),
    ("4a", "study", "a study row counted on the lineage split by the gate", PAGE,
     '            if (r["role"] == "study") != (r["split"] == "study"):',
     '            if r["role"] != "study" and r["split"] == "study":',
     f"{T}.Gates.test_a_study_run_on_the_lineage_split_is_refused"),
    ("4b", "study", "a study row counted on the lineage split in the species table", PAGE,
     '        groups.setdefault((slug, r["split"], r["pair"], r["learner"]), []).append(r)',
     '        groups.setdefault((slug, "lineage", r["pair"], r["learner"]), []).append(r)',
     f"{T}.Tables.test_a_study_row_stays_on_the_study_split"),
    ("4c", "study", "conditions 1 and 2 read on the study rows", PAGE,
     '            if r["split"] == "lineage" and r["measure"] == PRIMARY:',
     '            if r["measure"] == PRIMARY:',
     f"{T}.Tables.test_the_count_rows_carry_conditions_one_and_two_on_the_lineage_primary_rows"),
    ("5", "M", "a dataset of M counted", PAGE,
     '            if ds & in_m:\n'
     '                raise Stop(f"gate 5 M: score/{slug}/{f} holds {len(ds & in_m)} datasets of M")\n'
     '            if ds - first:',
     '            if ds & in_m and not in_m:\n'
     '                raise Stop(f"gate 5 M: score/{slug}/{f} holds {len(ds & in_m)} datasets of M")\n'
     '            if ds - first - in_m:',
     f"{T}.Gates.test_a_dataset_of_m_is_refused"),
    ("6a", "own datasets", "an arm's rows passed on datasets not its own", PAGE,
     "        if len(set(keys)) != len(keys) or set(keys) != want:",
     "        if len(set(keys)) != len(keys) or not want <= set(keys):",
     f"{T}.Gates.test_a_program_arm_outside_a_is_refused"),
    ("6b", "own datasets", "a program arm counted on L's datasets in the arms table", PAGE,
     "        mine = [r for r in data.run if runs_on(arm, r)]",
     "        mine = list(data.run)",
     f"{T}.Tables.test_each_arm_stands_on_its_own_datasets"),
    ("7a", "draws", "a value passed for a dataset with a draw not scored", PAGE,
     '            if r["scored_draws"] != r["draws"] and (r["scored"] != "0" or any(r[m] != "" for m in MEASURES)):',
     '            if r["scored_draws"] != r["draws"] and r["scored"] not in "01":',
     f"{T}.Gates.test_a_value_for_a_dataset_with_a_draw_not_scored_is_refused"),
    ("7b", "draws", "a missing draw read as zero", PAGE,
     "                va, vb = num(a[m]), num(b[m])",
     "                va, vb = num(a[m] or 0), num(b[m] or 0)",
     f"{T}.Tables.test_a_draw_not_scored_leaves_its_dataset_unscored_not_zero"),
    ("8a", "constant", "a constant-model draw passed without its reason", PAGE,
     "            if marked != why:",
     "            if why and not marked:",
     f"{T}.Gates.test_a_constant_model_draw_without_its_reason_is_refused"),
    ("8b", "constant", "a constant-model draw passed off an arm's study run", PAGE,
     '            if marked and not (r["arm"] in ARMS and r["split"] == "study" and r["scored"] == "1"):',
     '            if marked and not (r["arm"] in ARMS and r["scored"] == "1"):',
     f"{T}.Gates.test_a_constant_model_draw_off_an_arms_study_run_is_refused"),
    ("8c", "constant", "a study run refused for no candidate left unmarked", PAGE,
     '                   "study_refused": int(split == "study" and int(a["constant_draws"] or 0) > 0)}',
     '                   "study_refused": 0}',
     f"{T}.Tables.test_a_constant_model_row_is_kept_and_marked"),
    ("9", "checks", "a failed check of the scorer passed", PAGE,
     '        failed = sum(1 for r in checks if r["passed"] != "1")',
     "        failed = 0",
     f"{T}.Gates.test_a_failed_check_is_refused"),
    ("10", "builds", "a build record holding more than its counts passed", PAGE,
     "            if not counts_only(rec) or any(k not in rec for k in BUILD_COUNTS[:8] + BUILD_NEEDS[kind]):",
     "            if any(k not in rec for k in BUILD_COUNTS[:8] + BUILD_NEEDS[kind]):",
     f"{T}.Gates.test_a_build_record_beyond_its_counts_is_refused"),
    ("11a", "binary", "a DONE line naming another code passed", PAGE,
     '        if m is None or m["species"] != slug or (m["binary"], m["code"]) != (b["grm"], b["code"]):',
     '        if m is None or m["species"] != slug:',
     f"{T}.Gates.test_a_done_line_naming_another_code_is_refused"),
    ("11b", "binary", "a page that leaves a species' binary out passed", PAGE,
     '    if named != {(slug, data.binary[slug]["code"]) for slug in data.species}:',
     '    if not named <= {(slug, data.binary[slug]["code"]) for slug in data.species}:',
     f"{T}.Gates.test_a_page_that_leaves_a_species_binary_out_is_refused"),
    ("11c", "binary", "two binaries written as one", PAGE,
     "    if len(codes) == 1:",
     "    if len(codes) >= 1:",
     f"{T}.Tables.test_two_binaries_are_each_named_with_their_species"),
    ("12a", "conditions", "the stricter reading written as deciding", PAGE,
     '        if (r["decides"] == 1) != (r["reading"] == DECIDING):',
     '        if r["decides"] == 0 and r["reading"] == DECIDING:',
     f"{T}.Gates.test_a_page_that_writes_the_stricter_reading_as_deciding_is_refused"),
    ("12b", "conditions", "conditions.tsv passed without reading it again off count.tsv", PAGE,
     '            if k not in mine or whole(mine[k]["learners"]) != n:',
     "            if k not in mine:",
     f"{T}.Gates.test_conditions_that_differ_from_the_count_are_refused"),
    ("13", "count", "a count passed that is not the scorer's", PAGE,
     "        if data.count[name] != rows_want:",
     "        if len(data.count[name]) != len(rows_want):",
     f"{T}.Gates.test_a_count_that_is_not_the_scorers_is_refused"),
    ("R1", "run: SHA256SUMS", "a file changed since its sum passed", CHECKS,
     "            if data is None or sha(data) != want:",
     "            if data is None:",
     f"{C}.test_a_file_changed_since_its_sum_fails"),
    ("R2", "run: DONE and binary.txt", "a FAILED file not set aside passed", CHECKS,
     '            if s in set_aside and f is not None and f["species"] == s and m is not None and f["time"] < m["time"]:',
     '            if f is not None and f["species"] == s and m is not None and f["time"] < m["time"]:',
     f"{C}.test_a_failed_file_not_set_aside_fails_and_one_set_aside_is_named"),
    ("R3", "run: expect.py", "a result file more than expect.py's rows passed", CHECKS,
     "        extra = len(on_disk - listed) + both",
     "        extra = both",
     f"{C}.test_a_result_file_more_than_its_rows_fails"),
    ("R4", "run: expect.py", "a refused default passed", CHECKS,
     '            may_refuse = arm not in (DEFAULT, TWIN_ARM) and split == "study"',
     '            may_refuse = split == "study"',
     f"{C}.test_a_refused_default_fails"),
    ("R5", "run: refusals", "a refusal for another reason passed", CHECKS,
     '                         rec.get("reason")) == (REFUSED_FORMAT, dataset, arm, STUDY_SCHEME, STUDY_SEED, NO_CANDIDATE)):',
     '                         NO_CANDIDATE) == (REFUSED_FORMAT, dataset, arm, STUDY_SCHEME, STUDY_SEED, NO_CANDIDATE)):',
     f"{C}.test_a_refusal_for_another_reason_fails"),
    ("R6", "run: genome ids", "score.tar left out of the genome-id count", CHECKS,
     "        for data in scores[s].values():",
     "        for data in []:",
     f"{C}.test_a_genome_id_in_score_tar_is_counted"),
]


def _c9():
    spec = importlib.util.spec_from_file_location("c9_faults_for_c3_page", HERE.parent / "c9-pick" / "planted_faults.py")
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
    run_test = _c9().run_test
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
