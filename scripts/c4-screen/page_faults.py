"""page_faults.py: the faults of plan block C4's table script (page_tables.py) and the run's check script
(run_checks.py), each planted in turn, and the SYNTHETIC test that must catch it.

    python3 -I scripts/c4-screen/page_faults.py [--out faults.tsv]

Run from the repository root. A fault is a text of its file and the text that replaces it: the text must occur
once in the file. For each fault the named test runs once as the script stands, where it must pass, and once with
the fault planted, where it must fail. A fault is planted in memory only: the file is read, the text replaced, the
result loaded as a module of its own and put in the test module's place for the script, so a class fixture is
built again under the fault. Exit 1 unless every fault is caught.

Python 3.12; the scripts it plants into need numpy, as score.py does. Run with -I.
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

PAGE = "scripts/c4-screen/page_tables.py"
CHECKS = "scripts/c4-screen/run_checks.py"
# The test module and the name it holds each script under.
HOLDER = {PAGE: ("tests.c4_screen.test_page_tables", "page"), CHECKS: ("tests.c4_screen.test_run_checks", "checks")}
T = "tests.c4_screen.test_page_tables."
R = "tests.c4_screen.test_run_checks."

# (n, what is planted, file, text, replacement, gate or table, test that must catch it)
FAULTS = [
    ("P1", "a score.log in the reading passes", PAGE,
     "    if logs:\n", "    if False:\n", "1 files", T + "Gates.test_a_score_log_in_the_reading_is_refused"),
    ("P24", "a printout of fewer checks than 11 passes", PAGE,
     '    if not lines or lines[-1] != f"CHECKS {RUN_CHECKS} of {RUN_CHECKS}" or',
     '    if not lines or not re.fullmatch(r"CHECKS (\\d+) of \\1", lines[-1]) or', "1 files",
     T + "Gates.test_a_printout_of_fewer_checks_than_eleven_is_refused"),
    ("P2", "the species are those found, not those of the first run", PAGE,
     '    species = sorted({r["species_slug"] for r in run})\n',
     '    species = sorted(p.name for p in (reading / "score").iterdir() if p.is_dir())\n', "2 species",
     T + "Gates.test_a_species_missing_from_the_reading_is_refused"),
    ("P3", "a dataset more than the first run's, as one of M, passes", PAGE,
     "        if len(set(keys)) != len(keys) or set(keys) != want:\n            more",
     "        if len(set(keys)) != len(keys) or not want <= set(keys):\n            more", "3 datasets",
     T + "Gates.test_a_dataset_of_m_in_the_reading_is_refused"),
    ("P4", "a study row on the lineage split passes", PAGE,
     '            if (r["split"] == "study") != (r["role"] == "study"):',
     '            if r["split"] == "study" and r["role"] != "study":', "4 study",
     T + "Gates.test_a_study_row_on_the_lineage_split_is_refused"),
    ("P5", "the twin read as an arm passes", PAGE,
     '    named = sorted({r["arm"] for r in counts["conditions.tsv"]} - set(ARMS))',
     '    named = sorted({r["arm"] for r in counts["conditions.tsv"]} - set(ARMS) - {TWIN})', "5 twin",
     T + "Gates.test_the_twin_read_as_an_arm_is_refused"),
    ("P6", "a draw not scored keeps a figure", PAGE,
     '            if r["scored"] != "1" and any(r[m] != "" for m in MEASURES):',
     '            if r["scored"] != "1" and r["why"] != "split unavailable" and any(r[m] != "" for m in MEASURES):',
     "6 missing", T + "Gates.test_a_draw_not_scored_with_a_figure_is_refused"),
    ("P7", "a constant-model draw not marked passes", PAGE,
     '            if marked != (r["why"] == REFUSED_WHY) or (marked and',
     '            if (marked and r["why"] != REFUSED_WHY) or (marked and', "7 constant",
     T + "Gates.test_a_constant_model_draw_not_marked_is_refused"),
    ("P8", "datasets.tsv held to its draws on the names alone", PAGE,
     "            for c in score.DATASET_COLUMNS:", "            for c in score.DATASET_COLUMNS[:5]:", "8 draws",
     T + "Gates.test_a_dataset_figure_that_is_not_the_mean_of_its_draws_is_refused"),
    ("P9", "a check of the run that failed passes", PAGE,
     '        bad = [r for r in rows if r["passed"] != "1"]',
     '        bad = [r for r in rows if r["passed"] not in ("0", "1")]', "9 checks",
     T + "Gates.test_a_check_of_the_run_that_failed_is_refused"),
    ("P10", "a screen arm's cart held to the default's on whether it is scored alone", PAGE,
     '            return a is not None and b is not None and all(a[c] == b[c] for c in a if c != "arm")',
     '            return a is not None and b is not None and a["scored"] == b["scored"]', "9 checks",
     "tests.c4_screen.test_page_tables.CartOff.test_a_screen_arm_whose_cart_differs_from_the_default_is_refused"),
    ("P11", "a screen arm counted on any learner", PAGE,
     '            if score.arm_of(r["pair"]) in SCREEN_ARMS and r["learner"] not in SCREENED:',
     '            if score.arm_of(r["pair"]) in SCREEN_ARMS and r["learner"] not in LEARNERS:', "10 screen",
     T + "Gates.test_a_count_row_of_a_screen_arm_on_scm_is_refused"),
    ("P12", "the floor arm counted over every dataset passes", PAGE,
     '            if r["set"] != MOVES or int(r["datasets"]) != want:',
     '            if r["set"] not in (ALL, MOVES):', "11 floor",
     T + "Gates.test_the_floor_arm_counted_over_every_dataset_is_refused"),
    ("P13", "a withdrawn arm's pairs dropped from the family pass", PAGE,
     '    if got != want or any(r["split"] != "lineage" or r["measure"] != PRIMARY for r in family):',
     '    if not set(got) <= set(want):', "12 withdrawn",
     "tests.c4_screen.test_page_tables.Withdrawn.test_a_withdrawn_arm_dropped_from_the_family_is_refused"),
    ("P14", "a binary.txt with no code sha256 passes", PAGE,
     '                and re.fullmatch(r"[0-9a-f]{64}", fields.get("code", ""))):',
     '                and re.fullmatch(r"[0-9a-f]{64}", fields.get("code", "c" * 64))):', "13 binary",
     T + "Gates.test_a_binary_with_no_code_sha256_is_refused"),
    ("P15", "two binaries written as one", PAGE,
     "    if len(codes) == 1:", "    if len(codes) >= 1:", "13 binary",
     "tests.c4_screen.test_page_tables.Binaries.test_two_binaries_are_named_each_with_its_species"),
    ("P16", "a condition the count does not bear passes", PAGE,
     '        if key in read and r["met"] != str(int(read[key])):',
     '        if key in read and r["met"] not in ("0", "1"):', "14 conditions",
     T + "Gates.test_a_condition_that_the_count_does_not_bear_is_refused"),
    ("P17", "the count held to score.py's by its length alone", PAGE,
     "        if counts[name] != want[name]:", "        if len(counts[name]) != len(want[name]):", "15 count",
     T + "Gates.test_a_count_figure_that_score_py_does_not_write_is_refused"),
    ("P18", "fit times of the twin pass", PAGE,
     '        if arms != sorted(arms_run) or any(r["learner"] not in LEARNERS for r in rows):',
     '        if not set(arms_run) <= set(arms) or any(r["learner"] not in LEARNERS for r in rows):', "16 fits",
     T + "Gates.test_fit_times_of_the_twin_are_refused"),
    ("P19", "a figure not scored read as zero in the dataset table", PAGE,
     '                            va = num(mine[m]) if mine["scored"] == "1" else None',
     "                            va = num(mine[m]) or 0.0", "table",
     T + "Tables.test_a_missing_draw_is_never_read_as_zero"),
    ("P20", "a refused study run not marked in the dataset table", PAGE,
     '                               "refused": int(int(mine["constant_draws"] or 0) > 0)}',
     '                               "refused": 0}', "table",
     T + "Tables.test_a_refused_study_run_is_kept_and_marked"),
    ("P21", "the counts beside it marked as deciding", PAGE,
     '        out += [{"beside": beside, "decides": 0, **r} for r in counts[name]]',
     '        out += [{"beside": beside, "decides": 1, **r} for r in counts[name]]', "table",
     T + "Tables.test_the_floor_arm_decides_over_800_genomes_or_more_and_floor_all_stands_beside"),
    ("P22", "the stricter reading marked as deciding", PAGE,
     '"decides": int(not r["reading"].startswith("stricter")),', '"decides": 1,', "table",
     T + "Tables.test_the_stricter_reading_decides_nothing_and_no_twin_arm_is_read"),
    ("P23", "a species named by its slug", PAGE,
     '                    out.append({"species": mine[0]["species"], "species_slug": slug,',
     '                    out.append({"species": slug, "species_slug": slug,', "table",
     T + "Tables.test_species_are_named_in_full_and_their_counts_decide_nothing"),
    ("R1", "a species folder more passes", CHECKS,
     '    ok("the first run\'s species, each once", found == species,',
     '    ok("the first run\'s species, each once", set(species) <= set(found),', "check",
     R + "Check.test_a_species_folder_more_is_refused"),
    ("R2", "a FAILED newer than its DONE set aside", CHECKS,
     '        if failed < done["time"]:', "        if failed is not None:", "check",
     R + "Check.test_a_failed_newer_than_done_is_refused"),
    ("R3", "SHA256SUMS held to the files' names alone", CHECKS,
     "                 if n != SCORE_LOG and (not (out / n).is_file() or sha256(out / n) != sha)]",
     "                 if n != SCORE_LOG and not (out / n).is_file()]", "check",
     R + "Check.test_a_file_that_differs_from_its_sha256_is_refused"),
    ("R4", "score.log in SHA256SUMS passes", CHECKS,
     '    ok("score/score.log stays out of SHA256SUMS", all(',
     '    ok("score/score.log stays out of SHA256SUMS", True or all(', "check",
     R + "Check.test_score_log_in_sha256sums_is_refused"),
    ("R5", "DONE's SHA256SUMS sha256 not read", CHECKS,
     '                or not (out / "SHA256SUMS").is_file() or done["sums"] != sha256(out / "SHA256SUMS"):',
     '                or not (out / "SHA256SUMS").is_file():', "check",
     R + "Check.test_a_done_naming_other_sha256sums_is_refused"),
    ("R6", "two code sha256 pass", CHECKS,
     '    ok("every species ran one commit and one code sha256", len(commits) == 1 and len(codes) == 1,',
     '    ok("every species ran one commit and one code sha256", len(commits) == 1,', "check",
     "tests.c4_screen.test_run_checks.TwoBinaries.test_two_code_sha256_are_refused"),
    ("R7", "expect.tsv held to its header alone", CHECKS,
     "        if not got or got[0] != list(expect.COLUMNS) or got[1:] != want:",
     "        if not got or got[0] != list(expect.COLUMNS):", "check",
     R + "Check.test_an_expect_tsv_that_is_not_expect_py_s_is_refused"),
    ("R8", "a result more than listed passes", CHECKS,
     "        if have != needed:", "        if not needed <= have:", "check",
     R + "Check.test_a_result_more_than_listed_is_refused"),
    ("R9", "a refusal record with no constant-model rows passes", CHECKS,
     "        if set(marked) != refused or any(n != LEARNERS for n in marked.values()):",
     "        if any(n != LEARNERS for n in marked.values()):", "check",
     R + "Check.test_a_refusal_record_with_no_constant_model_rows_is_refused"),
    ("R10", "the genome id count skips the tables", CHECKS,
     "            if name == SCORE_LOG or not (out / name).is_file():",
     '            if name.endswith(".tsv") or not (out / name).is_file():', "check",
     R + "Check.test_a_genome_id_in_a_table_is_refused"),
    ("R11", "gather copies score.log", CHECKS,
     '          ("logs/fits.tsv", "fits.tsv"), ("binary.txt", "binary.txt"))',
     '          ("logs/fits.tsv", "fits.tsv"), ("binary.txt", "binary.txt"), ("score/score.log", "score.log"))',
     "check", R + "Gather.test_gather_copies_five_files_a_species_and_never_score_log"),
]


def planted_module(path, old, new):
    """The script at path with old replaced by new, loaded as a module of its own."""
    full = os.path.join(ROOT, path)
    with open(full, encoding="utf-8") as f:
        text = f.read()
    if text.count(old) != 1:
        raise ValueError(f"{path}: the fault's text occurs {text.count(old)} times, not once")
    name = "c4_screen_planted_" + os.path.basename(path)[:-3]
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
