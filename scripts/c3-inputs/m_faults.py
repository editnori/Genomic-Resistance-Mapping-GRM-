"""Plant each fault of the run on the datasets of M and run the test meant to catch it (docs/input-arms.md).

  python3 -I m_faults.py <out dir> [fault ...]

The run on the datasets of M's own rules: the list and the floor
(m_list.py), the three checks (m_checks.py), the second reading
(second.py), the three counts and the page tables (m_reading.py), the pack
(pack.py), the cost (cost.py) and the run script (run_m.sh). Copies this
checkout's scripts/ and tests/ to <out dir>/tree, a plain copy, so it runs
the same in a copy with no git history. One fault at a time, it changes one
piece of text there, runs the one test named for it, and puts the text back.
A fault is caught when its test fails and ran; each test must also pass with
no fault planted. The tests run as CI runs them, through
scripts/c9-pick/planted_faults.py's run_test. Writes <out dir>/faults.tsv
(fault, rule, what it plants, file, test, result) and <out dir>/logs/,
prints one line per fault, and exits 1 when a fault is not caught. Every
test it runs is SYNTHETIC or reads committed files; it reads no data.
"""

import importlib.util
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
LIST = "scripts/c3-inputs/m_list.py"
CHECKS = "scripts/c3-inputs/m_checks.py"
SECOND = "scripts/c3-inputs/second.py"
READING = "scripts/c3-inputs/m_reading.py"
PACK = "scripts/c3-inputs/pack.py"
COST = "scripts/c3-inputs/cost.py"
RUN = "scripts/c3-inputs/run_m.sh"
T = "tests.c3_inputs.test_m_"

# fault, rule, what it plants, file, text, planted text, test
FAULTS = [
    ("L1", "list", "a genome with I and R left off the list", LIST,
     '        held = [p for p in cohort if "I" in calls.get((p["genome_id"], spec["antibiotic"]), set())]',
     '        held = [p for p in cohort if calls.get((p["genome_id"], spec["antibiotic"]), set()) == {"S", "I"}]',
     f"{T}list.List.test_the_list_holds_s_with_i_and_i_with_r_and_not_i_alone_or_r_with_s"),
    ("L2", "list", "an intermediate call winning over S", LIST,
     '    for c in ("R", "S", "I"):',
     '    for c in ("I", "R", "S"):',
     f"{T}list.Calls.test_the_label_rule_r_or_s_wins"),
    ("L3", "list", "a genome labelled from an MIC kept in the cohort", LIST,
     '                         and p.get("reported") in ("R", "S") and not p.get("conflict") and not p.get("label_from_mic")),',
     '                         and p.get("reported") in ("R", "S") and not p.get("conflict")),',
     f"{T}list.List.test_the_list_holds_s_with_i_and_i_with_r_and_not_i_alone_or_r_with_s"),
    ("L4", "list", "records that consolidate to another call passed", LIST,
     "            bad += 1\n",
     "            bad += 0\n",
     f"{T}list.List.test_records_that_consolidate_to_another_call_stop"),
    ("L5", "list", "records whose name counts differ passed", LIST,
     "    if seen != want_records:",
     "    if set(seen) != set(want_records):",
     f"{T}list.List.test_records_whose_name_counts_differ_from_drug_names_stop"),
    ("L6", "list", "a cohort unlike lists.tsv's classes passed", LIST,
     '        if (len(cohort), r, s) != (int(lrow["samples"]), int(lrow["resistant"]), int(lrow["susceptible"])):',
     '        if len(cohort) != int(lrow["samples"]):',
     f"{T}list.List.test_a_cohort_unlike_lists_tsv_stops"),
    ("L7", "list", "a dataset outside M listed", LIST,
     '        if d not in listed_rows or listed_rows[d]["species_slug"] != slug or listed_rows[d]["in_M"] != "1":',
     '        if d not in listed_rows or listed_rows[d]["species_slug"] != slug:',
     f"{T}list.List.test_a_dataset_outside_m_is_refused"),
    ("F1", "floor", "a class under 25 entering the second reading", LIST,
     '        row["enters"] = int(row["second_resistant"] >= FLOOR and row["second_susceptible"] >= FLOOR)',
     '        row["enters"] = int(row["second_resistant"] >= FLOOR or row["second_susceptible"] >= FLOOR)',
     f"{T}list.Floor.test_twenty_five_of_each_class_enter_and_twenty_four_do_not"),
    ("F2", "floor", "a floor of 24", LIST,
     "FLOOR = 25\n",
     "FLOOR = 24\n",
     f"{T}list.Floor.test_twenty_five_of_each_class_enter_and_twenty_four_do_not"),
    ("C1", "check species", "the most first-run datasets taken for the check species", CHECKS,
     "    both = sorted((n, s) for s, n in first.items() if s in m)",
     "    both = sorted(((n, s) for s, n in first.items() if s in m), reverse=True)",
     f"{T}checks.CheckSpecies.test_the_check_species_on_the_committed_list_is_neisseria_gonorrhoeae"),
    ("C2", "Same build", "a build of another code sha256 passed", CHECKS,
     '    if row["code"] != code:',
     "    if not code:",
     f"{T}checks.SameBuild.test_a_build_whose_code_sha256_differs_stops"),
    ("C3", "Same results", "a cell compared at 9 decimals, not 12", CHECKS,
     "                equal = page.r12(float(a)) == page.r12(float(b))",
     "                equal = abs(float(a) - float(b)) < 1e-9",
     f"{T}checks.SameResults.test_one_cell_off_in_the_twelfth_decimal_stops"),
    ("C4", "Same results", "a recorded row not given back passed", CHECKS,
     "    if missing or extra:",
     "    if extra:",
     f"{T}checks.SameResults.test_a_missing_draw_or_a_missing_row_stops"),
    ("C5", "Same results", "a record unlike its SHA256SUMS read", CHECKS,
     '    for name in ("page-datasets.tsv",):',
     "    for name in ():",
     f"{T}checks.SameResults.test_a_record_that_differs_from_its_sha256sums_stops"),
    ("C6", "Same counts", "a fold's calls compared on n alone", CHECKS,
     '            for k in ("n", "tp", "fp", "tn", "fn"):',
     '            for k in ("n",):',
     f"{T}checks.SameCounts.test_one_fold_count_off_by_one_genome_stops"),
    ("C7", "Same counts", "a log loss sum compared at 1e-3", CHECKS,
     '            if abs(c["loss_sum"] - w["loss_sum"]) > LOSS_TOLERANCE * max(1.0, abs(w["loss_sum"])):',
     '            if abs(c["loss_sum"] - w["loss_sum"]) > 1e-3 * max(1.0, abs(w["loss_sum"])):',
     f"{T}checks.SameCounts.test_a_log_loss_sum_off_stops"),
    ("C8", "Same counts", "a kept report held to its score file by length alone", CHECKS,
     "            if from_rep != genomes:",
     "            if len(from_rep) != len(genomes):",
     f"{T}checks.SameCounts.test_a_kept_report_unlike_its_score_file_stops"),
    ("C9", "Same counts", "a score of exactly 0.5 called S", CHECKS,
     "        called = s >= CALL_CUTOFF",
     "        called = s > CALL_CUTOFF",
     f"{T}checks.SameCounts.test_the_folds_counts_from_scores_are_grms_cutoff_rule"),
    ("S1", "second reading", "the listed genomes kept in the test folds", SECOND,
     '    per = CHECKS.fold_calls(genomes, len(doc["folds"]), left_out)',
     '    per = CHECKS.fold_calls(genomes, len(doc["folds"]))',
     f"{T}second.Again.test_only_the_listed_genomes_leave_and_training_stays_as_fitted"),
    ("S2", "second reading", "the training counts changed with no refit", SECOND,
     "            rec = dict(rec, test=test, test_resistant=test_r,",
     '            rec = dict(rec, train=rec["train"] - (rec["test"] - test), test=test, test_resistant=test_r,',
     f"{T}second.Again.test_only_the_listed_genomes_leave_and_training_stays_as_fitted"),
    ("S3", "second reading", "a dataset under the floor scored", SECOND,
     '    rows_in = [e for e in expect if e["dataset"] in enters]',
     "    rows_in = list(expect)",
     f"{T}second.Reading.test_a_dataset_under_the_floor_has_no_row"),
    ("S4", "second reading", "a refused study fit's constant model on the folds as fitted", SECOND,
     '        refused[e["file"]] = docs[base]',
     '        refused[e["file"]] = CHECKS.load(picks / base)',
     f"{T}second.Reading.test_a_refused_study_fit_is_the_constant_model_on_the_defaults_folds_less_the_list"),
    ("S5", "second reading", "counts that do not come back from the scores read", SECOND,
     '        if diffs:\n            raise Stop(f"{len(diffs)} counts do not come back from the scores; "',
     '        if len(diffs) > 99:\n            raise Stop(f"{len(diffs)} counts do not come back from the scores; "',
     f"{T}second.Stops.test_counts_that_do_not_come_back_from_the_scores_stop"),
    ("S6", "second reading", "a listed genome outside the cohort passed", SECOND,
     "            if not left_out <= ids:",
     "            if not left_out & ids:",
     f"{T}second.Stops.test_a_listed_genome_outside_the_cohort_stops"),
    ("S7", "second reading", "kept counts held to floor.tsv by their sum alone", SECOND,
     '            if got != (int(r["second_resistant"]), int(r["second_susceptible"])):',
     '            if sum(got) != int(r["second_resistant"]) + int(r["second_susceptible"]):',
     f"{T}second.Stops.test_kept_counts_unlike_the_floor_tables_stop"),
    ("S8", "second reading", "a floor table's entry read without its rule", SECOND,
     '        if int(r["floor"]) != 25 or enters != int(r["enters"]):',
     '        if int(r["floor"]) != 25:',
     f"{T}second.Stops.test_a_floor_table_unlike_its_own_rule_stops"),
    ("R1", "first run", "a published count cell past the sixth column not compared", READING,
     "        for c in SCORE.COUNT_COLUMNS:\n            if not same_cell(have[k][c], want[k][c]):",
     "        for c in SCORE.COUNT_COLUMNS[:6]:\n            if not same_cell(have[k][c], want[k][c]):",
     f"{T}reading.FirstRun.test_a_published_count_cell_unlike_the_records_count_stops"),
    ("R2", "first run", "a published condition's met not compared", READING,
     "            if fmt(chave[k][c]) != cwant[k][c]:",
     '            if c != "met" and fmt(chave[k][c]) != cwant[k][c]:',
     f"{T}reading.FirstRun.test_a_published_condition_unlike_the_records_conditions_stops"),
    ("R3", "first run", "the published count read unlike its SHA256SUMS", READING,
     "    for name in names:\n        if sums.get(name) != CHECKS.sha256(Path(reference) / name):",
     "    for name in names[:1]:\n        if sums.get(name) != CHECKS.sha256(Path(reference) / name):",
     f"{T}reading.FirstRun.test_a_record_file_unlike_its_sha256sums_stops_before_any_count"),
    ("R4", "Holm", "Holm taken over both readings at once", READING,
     "            out[(name, reading)] = (counted, module.conditions(counted))",
     '            if reading == "second":\n'
     '                fam = [r for c in (out[(name, "first")][0], counted) for r in c\n'
     '                       if r["split"] == "lineage" and r["measure"] == PRIMARY]\n'
     '                for r, p in zip(fam, SCORE.C9.holm([r["sign_p"] for r in fam])):\n'
     '                    r["holm_p"] = p\n'
     "            out[(name, reading)] = (counted, module.conditions(counted))",
     f"{T}reading.Counts.test_holm_is_taken_apart_for_each_count_and_each_reading"),
    ("R5", "counts", "the datasets of M counted in the first run as read", READING,
     '        for name, module, rows in (("first run", two, record), ("M", six, m_rows[reading]),',
     '        for name, module, rows in (("first run", two, record + m_rows[reading]), ("M", six, m_rows[reading]),',
     f"{T}reading.Counts.test_three_counts_in_each_of_two_readings"),
    ("R6", "counts", "M alone counted on two measures", READING,
     "    two, six = scorer(RECORD_MEASURES), scorer(SCORE.MEASURES)",
     "    two, six = scorer(RECORD_MEASURES), scorer(RECORD_MEASURES)",
     f"{T}reading.Counts.test_the_m_count_reads_six_measures_and_the_record_counts_two"),
    ("R7", "tables", "the left-out table naming the datasets that enter", READING,
     '    return [{c: r[c] for c in LEFT_OUT_COLUMNS} for r in floor if r["enters"] != "1"]',
     '    return [{c: r[c] for c in LEFT_OUT_COLUMNS} for r in floor if r["enters"] != "0"]',
     f"{T}reading.Counts.test_the_left_out_table_names_each_dataset_under_the_floor_with_its_counts"),
    ("R8", "tables", "the second reading's size counts taken from the first", READING,
     '        for r in taken[("all datasets", reading)][0]:',
     '        for r in taken[("all datasets", "first")][0]:',
     f"{T}reading.Counts.test_the_size_table_gives_entered_better_and_worse_per_reading_arm_and_learner"),
    ("R9", "tables", "a condition met in the first reading alone written as met in both", READING,
     '                    "met_both": int(bool(a) and bool(b))})',
     '                    "met_both": int(bool(a))})',
     f"{T}reading.Counts.test_the_both_table_meets_a_condition_only_when_both_readings_meet_it"),
    ("R10", "tables", "a second reading holding a dataset under the floor read", READING,
     '        if {r["dataset"] for r in second} != enters:',
     '        if not {r["dataset"] for r in second} <= named:',
     f"{T}reading.SpeciesFolders.test_a_second_reading_holding_a_dataset_under_the_floor_stops"),
    ("R11", "tables", "a dataset of M read twice", READING,
     "        if named & seen:",
     "        if named == seen:",
     f"{T}reading.SpeciesFolders.test_each_dataset_of_m_is_read_once_and_nothing_else"),
    ("P1", "pack", "a score file passed by its name", PACK,
     '        if rel.endswith(".scores.json"):',
     '        if rel.endswith(".scores.json") and parts[0] == "work":',
     f"{T}pack.Check.test_a_score_file_is_refused_by_its_name_and_by_its_format"),
    ("P2", "pack", "a file of genomes one by one passed by its kind", PACK,
     '        why = genome_doc(text) if rel.endswith(".json") else None',
     "        why = None",
     f"{T}pack.Check.test_a_kept_report_is_refused"),
    ("P3", "pack", "a genome id passed", PACK,
     "        hits = sum(1 for m in GENOME_ID.finditer(text) if m.group(0) in ids)",
     "        hits = 0",
     f"{T}pack.Check.test_a_genome_id_in_a_table_or_a_log_is_refused"),
    ("P4", "pack", "an id read inside a longer token", PACK,
     'GENOME_ID = re.compile(r"(?<![\\w.])\\d+\\.\\d+(?!\\d)")',
     'GENOME_ID = re.compile(r"\\d+\\.\\d+")',
     f"{T}pack.Check.test_figures_and_versions_of_the_form_of_an_id_that_name_no_genome_pass"),
    ("P5", "pack", "a link to a file passed", PACK,
     "        if path.is_symlink():",
     "        if path.is_symlink() and path.is_dir():",
     f"{T}pack.Check.test_a_link_is_refused_to_a_file_or_a_folder"),
    ("P6", "pack", "a list file of genomes passed off work/", PACK,
     '        if any(parts[i:i + 2] == ("m", "list") for i in range(len(parts) - 1)):',
     '        if parts[:2] == ("m", "list"):',
     f"{T}pack.Check.test_a_list_file_of_genomes_is_refused_even_empty"),
    ("P7", "pack", "work/ walked into the pack", PACK,
     "            dirs[:] = [x for x in dirs if x not in LEFT]",
     '            dirs[:] = [x for x in dirs if x != "flags"]',
     f"{T}pack.Check.test_a_pack_of_counts_and_tables_passes_and_leaves_work_and_flags"),
    ("P8", "pack", "a file off the list written into the tar", PACK,
     "        for rel in files:\n            data = (out / rel).read_bytes()",
     '        for rel in files + ["run.lock"]:\n            data = (out / rel).read_bytes()',
     f"{T}pack.Tar.test_the_tar_holds_the_listed_files_alone_with_no_owner_or_time"),
    ("P9", "pack", "a tar written into a folder of OUT", PACK,
     '    if out.resolve() in tar_path.parents and tar_path != out.resolve() / "pack.tar":',
     '    if tar_path.parent == out.resolve() and tar_path != out.resolve() / "pack.tar":',
     f"{T}pack.Tar.test_a_tar_path_inside_out_other_than_pack_tar_is_refused"),
    ("K1", "cost", "the program on every species of A counted twice", COST,
     "    data = [t - sum(c[i] for c in ARM_COST.values()) - PROGRAM_ON_SPECIES[i] for i, t in enumerate(RUN_TOTAL)]",
     "    data = [t - sum(c[i] for c in ARM_COST.values()) for i, t in enumerate(RUN_TOTAL)]",
     f"{T}cost.Model.test_the_data_steps_are_the_runs_total_less_its_arm_rows"),
    ("K2", "cost", "a matrix arm taking the smaller share", COST,
     "            share = max(share, build_share)",
     "            share = min(share, build_share)",
     f"{T}cost.Model.test_a_matrix_arm_takes_the_larger_of_its_samples_share_and_the_genomes_share"),
    ("K3", "cost", "an arm's figure unlike the entry's", COST,
     'ARM_COST = {"k31": (12.12, 36.37),',
     'ARM_COST = {"k31": (12.21, 36.37),',
     f"{T}cost.Model.test_the_figures_are_the_entry_c3_arms_cost_models"),
    ("K4", "cost", "a fit's measured time not scaled to the datasets of M", COST,
     "                fit += wall * m / first if first else 0.0",
     "                fit += wall",
     f"{T}cost.Measured.test_each_species_costs_its_first_run_box_with_its_fits_scaled_to_its_datasets_of_m"),
    ("K5", "cost", "a timing table read unlike its SHA256SUMS", COST,
     '    if sums.get("logs/timing.tsv") != hashlib.sha256(path.read_bytes()).hexdigest():',
     '    if "logs/timing.tsv" not in sums:',
     f"{T}cost.Measured.test_a_timing_table_unlike_its_sha256sums_stops"),
    ("X1", "run script", "a lineage fit that writes no scores", RUN,
     '    --exp-pick-inner 0 --exp-scores "$SCORES/$arm" $(arm_args "$arm")',
     '    --exp-pick-inner 0 $(arm_args "$arm")',
     f"{T}run_script.Fits.test_the_one_lineage_fit_writes_scores_and_the_default_and_twin_keep_reports"),
    ("X2", "run script", "a fit of the default with no kept report taken as done", RUN,
     '          && { [ "$2" != k31 ] && [ "$2" != twin ] || [ -e "$KEPT/$2/$d.lineage-s$4.json" ]; }; } || echo "$d"',
     '          ; } || echo "$d"',
     f"{T}run_script.Fits.test_m_missing_asks_for_the_pick_inputs_score_file_and_kept_report"),
    ("X3", "run script", "a drawn seed's report kept from the seed-42 folder", RUN,
     '  [ "$2" = 42 ] || from=$ROOT/species/$SP/draws/seed-$2',
     '  [ "$2" = 42 ] || from=$REPORTS',
     f"{T}run_script.Fits.test_keep_reports_copies_each_report_where_its_seed_writes_it"),
    ("X4", "run script", "Same build not run", RUN,
     "step same_build first_build\n",
     "",
     f"{T}run_script.Order.test_the_three_checks_run_before_any_fit_of_m"),
    ("X5", "run script", "the fits of M run before the check species passed", RUN,
     '[ "$SP" != "$CHECK_SPECIES" ] || [ -e "$OUT/flags/m_checked.done" ] || { say "the check species has not passed its checks"; false; }\n',
     "",
     f"{T}run_script.Order.test_the_three_checks_run_before_any_fit_of_m"),
    ("X6", "run script", "Same counts not run over the fits of M", RUN,
     'same_counts "$OUT/expect.tsv" m\n',
     "",
     f"{T}run_script.Order.test_the_three_checks_run_before_any_fit_of_m"),
    ("X7", "run script", "the pack not checked", RUN,
     'python3 -I -B "$HERE/pack.py" check "$OUT" "$R/species/$SP.json"\n',
     "",
     f"{T}run_script.Order.test_the_three_checks_run_before_any_fit_of_m"),
    ("X8", "run script", "a pin unlike run_species.sh's", RUN,
     'pinned "$KIT/prodigal.linux" 80a3d7acb08ba9f63dde6b6ef4e95ec356128ac355f36960ea722f381150f246',
     'pinned "$KIT/prodigal.linux" 80a3d7acb08ba9f63dde6b6ef4e95ec356128ac355f36960ea722f381150f247',
     f"{T}run_script.Shared.test_the_pins_are_run_species_sh_s"),
    ("X9", "run script", "the records file not held to the master list", RUN,
     'master bvbrc/genome_amr.jsonl "$KIT/genome_amr.jsonl" \\',
     'true bvbrc/genome_amr.jsonl "$KIT/genome_amr.jsonl" \\',
     f"{T}run_script.Shared.test_the_records_file_is_held_to_its_line_in_the_master_list"),
    ("X10", "run script", "another species' record read on its species alone", RUN,
     'want = {"species": sys.argv[2], "code": sys.argv[3], "commit": sys.argv[4]}',
     'want = {"species": sys.argv[2]}',
     f"{T}run_script.Checked.test_another_species_runs_only_on_the_check_species_record_of_its_code"),
]


def _c9():
    spec = importlib.util.spec_from_file_location("c9_faults_for_c3_m", HERE.parent / "c9-pick" / "planted_faults.py")
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
    for top in ("scripts", "tests", "reference", "docs"):
        shutil.copytree(REPO / top, tree / top, ignore=shutil.ignore_patterns("__pycache__"))
    env = {"HOME": "/tmp"}
    rows, bad = [], 0
    for n, rule, what, path, old, new, test in FAULTS:
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
        rows.append((n, rule, what, path, test, result))
        print(f"{n}: {result} by {test}", flush=True)
    with open(out / "faults.tsv", "w") as f:
        f.write("fault\trule\tplants\tfile\ttest\tresult\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"caught {sum(r[5] == 'caught' for r in rows)} of {len(rows)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
