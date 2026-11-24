"""page_faults.py: the faults of plan block C5's table script (page_tables.py) and the run's check script
(run_checks.py), each planted in turn, and the SYNTHETIC test that must catch it.

    python3 -I scripts/c5-settings/page_faults.py [--out faults.tsv]

Run from the repository root. A fault is a text of its file and the text that replaces it: the text must occur
once in the file. For each fault the named test runs once as the script stands, where it must pass, and once with
the fault planted, where it must fail. A fault is planted in memory only: the file is read, the text replaced, the
result loaded as a module of its own and put in the test module's place for the script, so the test runs the
planted script. Exit 1 unless every fault is caught.

The scripts it plants into need numpy, as score.py does. Run with -I.
"""

import argparse
import contextlib
import importlib
import importlib.util
import io
import os
import sys
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

PAGE = "scripts/c5-settings/page_tables.py"
CHECKS = "scripts/c5-settings/run_checks.py"
RESCORE = "scripts/c5-settings/rescore.py"
# The test module and the name it holds each script under.
HOLDER = {PAGE: ("tests.c5_settings.test_page_tables", "page"), CHECKS: ("tests.c5_settings.test_run_checks", "checks"),
          RESCORE: ("tests.c5_settings.test_run_checks", "rescore")}
G = "tests.c5_settings.test_page_tables.Gates."
T = "tests.c5_settings.test_page_tables.Tables."
R = "tests.c5_settings.test_run_checks."
U = "tests.c5_settings.test_page_tables.UnavailableDraw."

# (n, what is planted, file, text, replacement, gate or table, test that must catch it)
FAULTS = [
    ("P1", "a score.log in the reading passes", PAGE,
     "    if logs:\n", "    if False:\n", "1 files", G + "test_a_score_log_in_the_reading_is_refused"),
    ("P2", "the species are those found, not those of the first run", PAGE,
     '    species = sorted({r["species_slug"] for r in run})\n',
     '    species = sorted(p.name for p in (reading / "score").iterdir() if p.is_dir())\n', "2 species",
     G + "test_a_species_missing_from_the_reading_is_refused"),
    ("P3", "a dataset of M on the first run's list is counted", PAGE,
     '    return [r for r in lrows if r["first_run"] == "1" and r["in_M"] == "0"]\n',
     '    return [r for r in lrows if r["first_run"] == "1"]\n', "3 datasets",
     G + "test_a_dataset_of_m_on_the_list_is_never_counted"),
    ("P4", "a study row on the lineage split passes", PAGE,
     '            if (r["split"] == "study") != (r["role"] == "study"):',
     '            if r["split"] == "study" and r["role"] != "study":', "4 study",
     G + "test_a_study_row_on_the_lineage_split_is_refused"),
    ("P5", "the twin read as an arm passes", PAGE,
     "    if named != sorted(ARMS):", "    if not set(ARMS) <= set(named):", "5 twin",
     G + "test_the_twin_read_as_an_arm_is_refused"),
    ("P6", "the plain best in a Holm family passes", PAGE,
     '        if r["pair"] != f"{PLAIN}-{DEFAULT}" or r["holm_p"] != "":',
     '        if r["pair"] != f"{PLAIN}-{DEFAULT}":', "6 plain",
     G + "test_the_plain_best_in_a_holm_family_is_refused"),
    ("P7", "a draw not scored keeps a figure", PAGE,
     '            if r["scored"] != "1" and any(r[m] != "" for m in MEASURES):',
     '            if r["scored"] != "1" and all(r[m] != "" for m in MEASURES):', "7 missing",
     G + "test_a_missing_draw_read_as_zero_is_refused"),
    ("P8", "a constant-model row left unmarked passes", PAGE,
     '            if marked != (r["why"] == REFUSED_WHY) or (marked and',
     "            if (marked and", "8 constant",
     G + "test_an_unmarked_constant_model_row_is_refused"),
    ("P9", "a draw listed twice passes", PAGE,
     "        if len(set(keys)) != len(keys) or set(keys) != want:\n            stop(9,",
     "        if set(keys) != want:\n            stop(9,", "9 draws", G + "test_a_draw_twice_is_refused"),
    ("P10", "a failed check passes", PAGE,
     '        bad = [r for r in rows if r["passed"] != "1"]', "        bad = []", "10 checks",
     G + "test_a_failed_check_is_refused"),
    ("P11", "the cutoff arm's log loss apart from the default's passes", PAGE,
     "        if passed != made:\n            stop(11,", "        if passed > made:\n            stop(11,",
     "11 cutoff", G + "test_the_cutoff_arms_log_loss_apart_from_the_defaults_is_refused"),
    ("P11b", "scm's cutoff moved from 0.5 passes", PAGE,
     '                f_pass += (r.get("cutoff") == HALF and outer.get(CUTOFF) is not None',
     "                f_pass += (outer.get(CUTOFF) is not None", "11 cutoff",
     G + "test_scms_cutoff_moved_from_half_is_refused"),
    ("P12", "a choices file missing inside a scored draw passes", PAGE,
     "        if folds != list(range(len(folds))) or len(folds) < 2:", "        if len(folds) < 2:", "12 choices",
     G + "test_a_choices_file_missing_inside_a_scored_draw_is_refused"),
    ("P12c", "a scored draw short of its default's outer folds passes", PAGE,
     "        if len(folds) != want[key]:", "        if False:", "12 choices",
     G + "test_a_scored_draw_short_of_its_defaults_folds_is_refused"),
    ("P12b", "a pick against the default-anchored rule passes", PAGE,
     "        if (picked == 0) != (lead is None or lead <= score.GAMMA) or (picked != 0 and picked != plain):",
     "        if picked != 0 and picked != plain:", "12 choices",
     G + "test_a_pick_against_the_default_anchored_rule_is_refused"),
    ("P13", "a setting in the counts of what the search chose passes", PAGE,
     "            if not ok:\n                stop(13,", "            if False:\n                stop(13,", "13 settings",
     G + "test_a_setting_in_the_counts_of_what_the_search_chose_is_refused"),
    ("P14", "a binary.txt without its code passes", PAGE,
     '                and re.fullmatch(r"[0-9a-f]{64}", fields.get("code", ""))):',
     '                and True):', "14 binary", G + "test_a_binary_txt_without_its_code_is_refused"),
    ("P15", "conditions.tsv apart from the count passes", PAGE,
     '        if key in read and r["met"] != str(int(read[key])):',
     '        if key in read and r["met"] not in ("0", "1"):', "15 conditions",
     G + "test_conditions_tsv_apart_from_the_count_is_refused"),
    ("P16", "a count that is not score.py's passes", PAGE,
     "        if counts[name] != want[name]:", "        if len(counts[name]) != len(want[name]):", "16 count",
     G + "test_a_count_that_is_not_score_pys_is_refused"),
    ("T1", "the stricter reading decides", PAGE,
     '"decides": int(r["reading"] == "ruling 15"),', '"decides": 1,', "table",
     T + "test_ruling_15_decides_and_the_stricter_reading_stands_beside"),
    ("T2", "the cutoff arm's log loss rows decide", PAGE,
     '        row = {"decides": int(not equal), **r,', '        row = {"decides": 1, **r,', "table",
     T + "test_the_cutoff_arms_log_loss_rows_decide_nothing"),
    ("T3", "the plain best decides", PAGE,
     '            out.append({"beside": beside, "decides": 0, **r, "family": 0,',
     '            out.append({"beside": beside, "decides": int(beside == "plain-best"), **r, "family": 0,', "table",
     T + "test_the_plain_best_and_the_count_without_constant_models_decide_nothing"),
    ("T4", "a missing value read as zero in the datasets table", PAGE,
     '                            va = num(mine[m]) if mine["scored"] == "1" else None',
     '                            va = num(mine[m]) if mine["scored"] == "1" else 0.0', "table",
     T + "test_a_missing_draw_is_never_read_as_zero_in_the_tables"),
    ("T5", "a missing figure printed as zero", PAGE,
     '    if x in ("", None):\n        return "none"\n    value = Decimal',
     '    if x in ("", None):\n        return "0.000"\n    value = Decimal', "table",
     T + "test_a_missing_draw_is_never_read_as_zero_in_the_tables"),
    ("T6", "a species named by its slug", PAGE,
     '            out.append({"species": species_name(found["by_species"], slug), "species_slug": slug,',
     '            out.append({"species": slug, "species_slug": slug,', "table",
     T + "test_each_species_is_named_in_full_and_decides_nothing"),
    ("T7", "a trial's settings ride with the counts", PAGE,
     '        a["failed_trials"] += int(r.get("failed_trials") or 0)',
     '        a["failed_trials"] = r.get("picked_setting")', "table",
     T + "test_the_counts_of_what_the_search_chose_hold_no_setting"),
    ("T8", "two binaries named as one", PAGE,
     "    if len(codes) == 1:\n        return (f\"All", "    if len(codes) >= 1:\n        return (f\"All", "table",
     T + "test_one_binary_is_named_once_and_two_with_their_species"),
    ("R1", "a species more than the first run's passes", CHECKS,
     '    ok("the first run\'s species, each once", found == species,',
     '    ok("the first run\'s species, each once", set(species) <= set(found),', "check",
     R + "Check.test_a_species_more_than_the_first_runs_is_refused"),
    ("R2", "a FAILED newer than its DONE is set aside", CHECKS,
     '        if failed < done["time"]:', "        if True:", "check",
     R + "Check.test_a_failed_file_newer_than_done_is_refused"),
    ("R3", "a FAILED older than its DONE is not named", CHECKS,
     "            return done, failed, f\"{slug}/FAILED at {failed}, older than its DONE at {done['time']}\", None",
     "            return done, failed, None, None", "check",
     R + "Check.test_a_failed_file_older_than_done_is_named_and_set_aside"),
    ("R4", "a file changed after its SHA256SUMS passes", CHECKS,
     "                 if n != SCORE_LOG and (not (out / n).is_file() or sha256(out / n) != sha)]",
     "                 if n != SCORE_LOG and not (out / n).is_file()]", "check",
     R + "Check.test_a_file_changed_after_its_sha256sums_is_refused"),
    ("R5", "score.log in SHA256SUMS passes", CHECKS,
     '    ok("score/score.log stays out of SHA256SUMS", all(',
     '    ok("score/score.log stays out of SHA256SUMS", True or all(', "check",
     R + "Check.test_score_log_in_sha256sums_is_refused"),
    ("R6", "two code sha256 pass", CHECKS,
     '    ok("every species ran one commit and one code sha256", len(commits) == 1 and len(codes) == 1,',
     '    ok("every species ran one commit and one code sha256", len(commits) == 1,', "check",
     R + "Check.test_two_code_sha256_are_refused"),
    ("R7", "an expect.tsv that is not expect.py's passes", CHECKS,
     "        if not got or got[0] != list(expect.COLUMNS) or sorted(got[1:]) != sorted(want):",
     "        if not got or got[0] != list(expect.COLUMNS):", "check",
     R + "Check.test_an_expect_tsv_that_is_not_expect_pys_is_refused"),
    ("R8", "a missing choices file passes", CHECKS,
     "        if have - refused != needed:", "        if not (have - refused) <= needed:", "check",
     R + "Check.test_a_missing_choices_file_is_refused"),
    ("R9", "the genome id count skips the choices files", CHECKS,
     "            if name == SCORE_LOG or not (out / name).is_file():",
     '            if name.endswith(".json") or not (out / name).is_file():', "check",
     R + "Check.test_a_genome_id_in_a_choices_file_is_refused"),
    ("R10", "a refusal record with no constant-model rows passes", CHECKS,
     "            if n_constant != LEARNERS * len(refused):", "            if n_constant > LEARNERS * len(refused):",
     "check", R + "Check.test_a_refusal_record_with_no_constant_model_rows_is_refused"),
    ("R11", "gather copies score.log", CHECKS,
     '          ("binary.txt", "binary.txt"))',
     '          ("binary.txt", "binary.txt"), ("score/score.log", "score.log"))', "check",
     R + "Gather.test_gather_copies_four_files_a_species_and_the_choices_and_never_score_log"),
    # Part 2: an unavailable default draw, the columns of 0, the allowance, the times and the rescored tables.
    ("P12d", "an unavailable default that owes choices files passes", PAGE,
     '            if marks == {("0", UNAVAILABLE_WHY)} and want.get(key) == 0:',
     '            if marks == {("0", UNAVAILABLE_WHY)}:', "12 choices",
     U + "test_an_unavailable_default_that_owes_choices_files_is_refused"),
    ("T9", "a column of 0 in every row printed", PAGE,
     '    keep = [i for i in range(len(head)) if not rows or any(str(row[i]) != "0" for row in rows)]',
     "    keep = list(range(len(head)))", "table", U + "test_a_column_of_0_in_every_row_is_not_printed"),
    ("R12", "a species with no DONE passes as stopped at score without the allowance", CHECKS,
     "    allowed = failed_at if failed_at in present else None\n",
     '    allowed = failed_at if failed_at in present else next((s for s in present if not (copies / s / "DONE").exists()), None)\n',
     "check", R + "Allowance.test_a_species_stopped_at_score_fails_without_the_allowance"),
    ("R13", "SHA256SUMS.lead of another sha256 passes", CHECKS,
     "    lines.append((lead.is_file() and listed > 0 and wrong == 0 and named_log == 0 and own == lead_sha,",
     "    lines.append((lead.is_file() and listed > 0 and wrong == 0 and named_log == 0,", "check",
     R + "Allowance.test_the_allowance_refuses_a_lead_sums_of_another_sha256"),
    ("R14", "a drive ended otherwise passes the allowance", CHECKS,
     "    lines.append((runs > 0 and to_run == len(done_ok) == len(set(done_ok)) and other == 0,",
     "    lines.append((runs > 0 and to_run == len(done_ok) == len(set(done_ok)),", "check",
     R + "Allowance.test_the_allowance_refuses_a_drive_ended_otherwise"),
    ("R15", "an id on any line of a drive's .time file counted as a time", CHECKS,
     "        if label and value >= 0 and m.start() - start >= value + 2:",
     "        if value >= 0 and m.start() - start >= value + 2:", "check",
     R + "Times.test_an_id_on_another_line_of_a_time_file_still_fails"),
    ("R16", "a seconds_inner value equal to an id fails as an id", CHECKS,
     'TIME_KEYS = ("seconds_inner", "seconds_outer")', 'TIME_KEYS = ("seconds_outer",)', "check",
     R + "Times.test_a_seconds_inner_value_equal_to_an_id_is_counted_apart"),
    ("R17", "gather with tables copies the copy's tables", CHECKS,
     "                src_path = Path(tables) / slug / dst",
     "                src_path = copies / slug / src", "check",
     R + "Tables.test_gather_with_tables_takes_the_rescored_tables"),
    ("R18", "a rescored table unlike its SHA256SUMS line passes", RESCORE,
     '        agree = [t for t in TABLES if listed.get(f"score/{t}") == sha256(dest / t)]',
     "        agree = list(TABLES)", "check",
     R + "Tables.test_rescore_refuses_a_table_unlike_its_sha256sums_line"),
]


def planted_module(path, old, new):
    """The script at path with old replaced by new, loaded as a module of its own."""
    full = os.path.join(ROOT, path)
    with open(full, encoding="utf-8") as f:
        text = f.read()
    if text.count(old) != 1:
        raise ValueError(f"{path}: the fault's text occurs {text.count(old)} times, not once")
    name = "c5_settings_planted_" + os.path.basename(path)[:-3]
    spec = importlib.util.spec_from_file_location(name, full)
    module = importlib.util.module_from_spec(spec)
    exec(compile(text.replace(old, new), full, "exec"), module.__dict__)
    return module


@contextlib.contextmanager
def plant(path, old, new):
    holder, attr = HOLDER[path]
    test_module = importlib.import_module(holder)
    with mock.patch.object(test_module, attr, planted_module(path, old, new)):
        yield


def outcome(test_id, planting=None):
    """pass, fail or not run: the test's result, with the fault planted if planting is given."""
    suite = unittest.defaultTestLoader.loadTestsFromName(test_id)
    result = unittest.TestResult()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        with planting() if planting else contextlib.nullcontext():
            suite.run(result)
    if result.testsRun == 0:
        return "not run"
    return "pass" if result.wasSuccessful() and not result.skipped else "fail"


def run_all():
    rows = []
    for n, what, path, old, new, gate, test_id in FAULTS:
        clean = outcome(test_id)
        with_fault = outcome(test_id, lambda: plant(path, old, new))
        rows.append({"fault": n, "planted": what, "file": path, "gate": gate,
                     "test": test_id.rsplit(".", 2)[-2] + "." + test_id.rsplit(".", 1)[-1],
                     "as_written": clean, "with_fault": with_fault,
                     "caught": "yes" if clean == "pass" and with_fault == "fail" else "no"})
    return rows


HEADER = ["fault", "planted", "file", "gate", "test", "as_written", "with_fault", "caught"]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    rows = run_all()
    lines = ["\t".join(HEADER)] + ["\t".join(str(r[h]) for h in HEADER) for r in rows]
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    caught = sum(r["caught"] == "yes" for r in rows)
    print(f"planted faults: {caught} of {len(rows)} caught")
    return 0 if caught == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
