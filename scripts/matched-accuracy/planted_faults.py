"""Plant each fault of stage 1 of the matched accuracy protocol and run the test meant to fail with it.

  python3 -I planted_faults.py <out dir> [fault ...]

Copies this checkout's scripts/, tests/, reference/ and docs/ to <out dir>/tree, a plain copy, so it runs the
same in a copy with no git history. One fault at a time, it changes one piece of text of
scripts/matched-accuracy/stage1.py there, runs the one test named for it, and puts the text back. Each test must
pass with no fault planted. A fault counts only when its test ran alone and failed by an assertion: unittest's
last line is FAILED with failures and no errors. A test that fails by an error does not count. The tests run as
CI runs them, through scripts/c9-pick/planted_faults.py's run_test. Writes <out dir>/faults.tsv (fault, rule,
what it plants, file, test, result) and <out dir>/logs/, prints one line per fault and a last line with the
count, and exits 1 when a fault does not fail its test by an assertion. Every test it runs is SYNTHETIC or reads
committed files; it reads no data and no score.
"""

import importlib.util
import re
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
S = "scripts/matched-accuracy/stage1.py"
T = "tests.matched_accuracy."
LOADER = T + "test_loader.Loader."
STOPS = T + "test_stops."
ROWS = T + "test_rows."
EXPORT = T + "test_export.Export."
SCAN = T + "test_scan."
MADE_UP = T + "test_made_up.MadeUp."

# fault, rule, what it plants, file, text, planted text, test
FAULTS = [
    ("P01", "loader", "the loader keeps each learner's scores", S,
     '"scores_length": len(e["scores"]) if isinstance(e.get("scores"), list) else None}',
     '"scores_length": len(e["scores"]) if isinstance(e.get("scores"), list) else None, "scores": e.get("scores")}',
     LOADER + "test_no_score_or_value_of_roc_metrics_conformal_or_no_call_is_reachable"),
    ("P02", "loader", "the loader keeps every key of a genome row", S,
     'out["samples"] = [{k: s.get(k) for k in ROW_KEYS} if isinstance(s, dict) else None for s in doc["samples"]]',
     'out["samples"] = [dict(s) if isinstance(s, dict) else None for s in doc["samples"]]',
     LOADER + "test_the_loader_returns_the_named_keys_alone"),
    ("P03", "loader", "the loader keeps each learner's roc", S,
     'kept["learners"] = [{"learner": e.get("learner"),',
     'kept["learners"] = [{"learner": e.get("learner"), "roc": e.get("roc"),',
     LOADER + "test_no_score_or_value_of_roc_metrics_conformal_or_no_call_is_reachable"),
    ("P04", "loader", "the loader reads the random split", S,
     'lineage = [s for s in doc["splits"] if isinstance(s, dict) and s.get("scheme") == "lineage"]',
     'lineage = [s for s in doc["splits"] if isinstance(s, dict) and s.get("scheme") in ("lineage", "random")][-1:]',
     LOADER + "test_the_loader_keeps_the_lineage_split_and_drops_the_random_one"),
    ("P05", "seal", "a lab record inserted at the seal passes", S,
     "        if head >= SEAL:\n",
     "        if head > SEAL:\n",
     STOPS + "Records.test_a_record_inserted_at_the_seal_stops"),
    ("P06", "seal", "a lab record with a blank insertion time passes", S,
     "        if head is None:\n            raise Stop(\"a lab record's date_inserted",
     "        if head is None:\n            head = \"\"\n        if head is None:\n"
     "            raise Stop(\"a lab record's date_inserted",
     STOPS + "Records.test_a_record_with_a_blank_insertion_date_stops"),
    ("P07", "line count", "fewer lab records than the pinned count pass", S,
     "    if n != RECORDS_LINES:\n",
     "    if n > RECORDS_LINES:\n",
     STOPS + "Records.test_a_line_count_other_than_the_pinned_one_stops"),
    ("P08", "folds", "a test_fold of another length passes", S,
     "    if not isinstance(fold, list) or len(fold) != n:\n",
     "    if not isinstance(fold, list):\n",
     STOPS + "Reports.test_a_test_fold_shorter_than_the_genome_rows_stops"),
    ("P09", "folds", "a scores list of another length passes", S,
     '        if e["scores_length"] != n:\n',
     '        if e["scores_length"] is None:\n',
     STOPS + "Reports.test_a_scores_list_longer_than_the_genome_rows_stops"),
    ("P10", "folds", "a lineage split with an empty fold passes", S,
     'if split["folds"] != FOLDS or not all(is_int(f) for f in fold) or set(fold) != set(range(FOLDS)):',
     'if split["folds"] != FOLDS or not all(is_int(f) for f in fold):',
     STOPS + "Reports.test_a_lineage_split_with_other_than_five_folds_stops"),
    ("P11", "lineage cluster", "a lineage cluster in two folds passes", S,
     'if seen.setdefault(s["lineage"], f) != f:',
     'if seen.setdefault((s["lineage"], f), f) != f:',
     STOPS + "Reports.test_a_lineage_cluster_in_two_folds_stops"),
    ("P12", "copy group", "a copy group in two folds passes", S,
     "if len({fold[i] for i in g}) != 1:",
     "if len({fold[i] for i in g}) > 2:",
     STOPS + "Reports.test_a_copy_group_in_two_folds_stops"),
    ("P13", "learners", "a lineage split with four of the five learners passes", S,
     "    if five and five != set(FIVE):\n",
     "    if five and len(five) < 3:\n",
     STOPS + "Reports.test_a_report_with_four_of_the_five_learners_stops"),
    ("P14", "learners", "a report without the five learners enters", S,
     "        if enters(rep):\n",
     "        if enters(rep) or True:\n",
     MADE_UP + "test_a_report_without_learners_is_counted_and_left_out"),
    ("P15", "lab call", "any R or S of the label table passes", S,
     'if table.get((s["id"], rep["antibiotic"])) != ("R" if s["resistant"] else "S"):',
     'if table.get((s["id"], rep["antibiotic"])) not in ("R", "S"):',
     STOPS + "Calls.test_a_row_whose_call_differs_from_the_label_table_stops"),
    ("P16", "lab call", "a row with from_mic true is checked against the table", S,
     "                mic += 1\n                continue\n",
     "                mic += 1\n",
     ROWS + "Labels.test_a_row_with_from_mic_true_is_counted_apart_and_not_checked"),
    ("P17", "lab call", "the label table's label stands for its reported call", S,
     'table[key] = None if p.get("conflict") or p.get("label_from_mic") else p.get("reported")',
     'table[key] = p.get("label")',
     STOPS + "Calls.test_a_row_whose_table_pair_is_labelled_from_an_mic_stops"),
    ("P18", "reasons", "a copy group row without a BioProject leaves under reason 2", S,
     "            if i in copies:\n",
     '            if i in copies and project[s["id"]] is not None:\n',
     ROWS + "Reasons.test_a_row_leaves_under_the_first_reason_that_applies"),
    ("P19", "BioProject", "a value with two BioProjects gives one of them", S,
     "one = rows is not None and len(set(rows)) == 1 and len(rows[0]) == 1",
     "one = rows is not None and len(set(rows)) == 1 and len(rows[0]) >= 1",
     ROWS + "Genomes.test_a_genome_gets_one_bioproject_or_none"),
    ("P20", "BioProject", "a genome the second copy lacks is not counted apart", S,
     'c["bioproject_not_in_table"] += absent[s["id"]]',
     'c["bioproject_not_in_table"] += 0',
     ROWS + "Reasons.test_a_genome_the_second_copy_lacks_is_counted_apart_under_reason_2"),
    ("P21", "drug names", "a name missing from the drug names file reads as a drug", S,
     "kind, canonical = resolve.get(raw, (None, None))",
     'kind, canonical = resolve.get(raw, ("drug", raw))',
     STOPS + "Records.test_a_record_name_missing_from_the_drug_names_file_stops"),
    ("P22", "drug names", "a record's own name stands for its antibiotic", S,
     'if kind != "drug" or (genome, canonical) not in wanted:',
     'if kind != "drug" or (genome, raw) not in wanted:',
     ROWS + "Records.test_a_record_reaches_its_antibiotic_through_the_drug_names_file"),
    ("P23", "tuple", "a blank standard passes", S,
     "        if not methods or not standards:\n",
     "        if not methods:\n",
     ROWS + "Records.test_a_blank_method_or_standard_marks_the_pair"),
    ("P24", "tuple", "a standard is not uppercased", S,
     'standards = {s.upper() for s in values(rec.get("testing_standard"))}',
     'standards = values(rec.get("testing_standard"))',
     ROWS + "Records.test_a_standard_is_trimmed_and_uppercased_as_grm_curate_does"),
    ("P25", "tuple", "a field keeps its outer whitespace", S,
     "        return v.strip() or None\n",
     "        return v or None\n",
     ROWS + "Fields.test_value_text_trims_and_joins_as_grm_reads_a_field"),
    ("P26", "tuple", "the year leaves the tuple", S,
     'p["tuples"] |= {(m, s, y) for m in methods for s in standards for y in years}',
     'p["tuples"] |= {(m, s, None) for m in methods for s in standards}',
     ROWS + "Records.test_two_years_of_one_method_and_standard_are_two_tuples"),
    ("P27", "tuple", "two years of one method and standard pass reason 5", S,
     '            elif len(pair["tuples"]) != 1:\n',
     '            elif len({t[:2] for t in pair["tuples"]}) != 1:\n',
     ROWS + "Reasons.test_a_pair_with_two_tuples_leaves_under_reason_5"),
    ("P28", "intermediate", "only the literal word Intermediate marks a pair", S,
     'if parse_call(value_text(rec.get("resistant_phenotype"))) == "I":',
     'if value_text(rec.get("resistant_phenotype")) == "Intermediate":',
     ROWS + "Records.test_an_intermediate_call_by_curates_words_marks_the_pair"),
    ("P29", "groups", "the BioProject leaves the group key", S,
     'return (r["dataset"], r["fold"], r["lineage"], r["bioproject"], r["method"], r["standard"], r["year"]) + \\',
     'return (r["dataset"], r["fold"], r["lineage"], r["method"], r["standard"], r["year"]) + \\',
     ROWS + "Groups.test_a_group_is_keyed_by_its_seven_fields"),
    ("P30", "groups", "a group with one lab outcome stays", S,
     'both = {k: g for k, g in by.items() if len({r["resistant"] for r in g}) == 2}',
     'both = {k: g for k, g in by.items() if len({r["resistant"] for r in g}) >= 1}',
     ROWS + "Groups.test_a_group_with_one_lab_outcome_leaves"),
    ("P31", "pairs", "a group's pairs are the smaller of its two outcome counts", S,
     "    return r * (len(group) - r)\n",
     "    return min(r, len(group) - r)\n",
     ROWS + "Groups.test_the_pairs_of_a_group_are_its_resistant_rows_times_its_susceptible_rows"),
    ("P32", "reading B", "reading B takes out resistant rows with an intermediate record alone", S,
     '            if r["intermediate"]:\n',
     '            if r["intermediate"] and r["resistant"]:\n',
     ROWS + "Readings.test_reading_b_takes_out_every_pair_with_an_intermediate_record_before_grouping"),
    ("P33", "no group", "a reading B without a group passes the check", S,
     "        if not groups:\n",
     '        if not groups and name == "A":\n',
     STOPS + "Groups.test_no_group_in_reading_b_stops"),
    ("P34", "weights", "every group of the export weighs the same", S,
     "return {k: fractions.Fraction(1, by_dataset[k[0]] * datasets_in[species_of[k[0]]] * n_species) for k in groups}",
     "return {k: fractions.Fraction(1, len(groups)) for k in groups}",
     ROWS + "Weights.test_groups_datasets_and_species_weigh_equally_and_exactly"),
    ("P35", "shares", "the largest BioProject's weight is the largest weight", S,
     "    its = max(weight[b] for b in rows if rows[b] == top)\n",
     "    its = max(weight.values())\n",
     ROWS + "Weights.test_the_largest_bioprojects_share_is_two_pairs_of_integers"),
    ("P36", "typing", "the primary match beside a narrow one counts every row", S,
     "        primary = support(groups_of(surviving)[0])\n",
     "        primary = support(groups_of(rows)[0])\n",
     ROWS + "Typing.test_a_narrow_match_is_counted_beside_the_primary_match_on_its_surviving_rows"),
    ("P37", "typing", "a row lacks a typing value only when it lacks both", S,
     "            if any(r2[k] is None for k in extra):\n",
     "            if all(r2[k] is None for k in extra):\n",
     MADE_UP + "test_the_typing_counts_equal_the_hand_count"),
    ("P38", "row list", "the row list is not sorted", S,
     '    return "".join(sorted(lines)).encode("utf-8")\n',
     '    return "".join(lines).encode("utf-8")\n',
     EXPORT + "test_the_row_list_is_sorted_utf8_lf_and_tab_separated"),
    ("P39", "row list", "the row lists are readable by every user", S,
     "        os.chmod(path, 0o600)\n",
     "        os.chmod(path, 0o644)\n",
     EXPORT + "test_the_row_lists_are_mode_0600_in_a_directory_of_mode_0700"),
    ("P40", "scan", "the positive control passes a scan that finds nothing", S,
     "    if hits(planted, ids, accessions) == 0:\n",
     "    if hits(planted, ids, accessions) < 0:\n",
     SCAN + "Scan.test_the_positive_control_stops_a_scan_that_finds_nothing"),
    ("P41", "scan", "the scan looks up genome ids alone", S,
     "    n = sum(text.count(x) for x in ids if x in text) + sum(text.count(x) for x in accessions if x in text)\n",
     "    n = sum(text.count(x) for x in ids if x in text)\n",
     SCAN + "Hits.test_an_accession_inside_a_file_is_a_hit"),
    ("P42", "scan", "the scan counts no id-shaped token", S,
     "    return n + len(ID_SHAPED.findall(text)) + len(ACCESSION_SHAPED.findall(text))\n",
     "    return n\n",
     SCAN + "Hits.test_a_token_with_the_shape_of_a_real_genome_id_is_a_hit"),
    ("P43", "scan", "the scan looks up the ids of the reports alone", S,
     "    ids = {s[\"id\"] for rep in entering + left_out for s in rep[\"samples\"]} | record_ids | label_ids | genome_ids\n",
     "    ids = {s[\"id\"] for rep in entering + left_out for s in rep[\"samples\"]}\n",
     SCAN + "RunWithAHit.test_a_hit_prints_each_file_and_its_count_and_writes_nothing"),
    ("P44", "export", "the export is written on a hit", S,
     '            print(f"{name} {n}")\n        return 1\n',
     '            print(f"{name} {n}")\n        write(Path(args.out), lists, export)\n        return 1\n',
     SCAN + "RunWithAHit.test_a_hit_prints_each_file_and_its_count_and_writes_nothing"),
    ("P45", "dropped figure", "the dropped records are named as measured", S,
     'f"figure of ruling R49, not measured here"',
     'f"figure, measured here"',
     EXPORT + "test_the_dropped_figure_is_named_as_recorded_never_as_measured"),
    ("P46", "sha256", "the sums list of the reports is not checked", S,
     "    if sha256(REPORT_SUMS) != REPORT_SUMS_SHA256:\n",
     "    if not REPORT_SUMS.is_file():\n",
     STOPS + "Sums.test_a_sums_list_other_than_its_sha256_stops"),
    ("P47", "sha256", "a report is not checked against its line", S,
     "        if sha256(directory / name) != want[name]:\n",
     "        if not (directory / name).is_file():\n",
     STOPS + "Sums.test_a_report_other_than_its_line_of_the_sums_list_stops"),
    ("P48", "sha256", "the lab records are checked against their own sha256", S,
     '("lab records", args.records, RECORDS_SHA256)',
     '("lab records", args.records, sha256(args.records))',
     STOPS + "Sums.test_a_lab_records_file_other_than_its_sha256_stops"),
    ("P49", "sha256", "the second genome copy is checked against its own sha256", S,
     '("genome table, second copy", args.genomes, GENOMES_SHA256)',
     '("genome table, second copy", args.genomes, sha256(args.genomes))',
     STOPS + "Sums.test_a_second_genome_copy_other_than_its_sha256_stops"),
    ("P50", "sha256", "the drug names file is checked against its own sha256", S,
     '("drug names", args.drug_names[0], args.drug_names[1])',
     '("drug names", args.drug_names[0], sha256(args.drug_names[0]))',
     STOPS + "Sums.test_a_drug_names_file_other_than_its_sha256_stops"),
    ("P51", "sha256", "the label table is not checked", S,
     "        if got != want:\n",
     '        if got != want and what != "label table":\n',
     STOPS + "Sums.test_a_label_table_other_than_its_sha256_stops"),
    ("P52", "sha256", "the 5 Nov table is checked against its own sha256", S,
     'pins.append(("genome table, 5 Nov", args.hc20[0], args.hc20[1]))',
     'pins.append(("genome table, 5 Nov", args.hc20[0], sha256(args.hc20[0])))',
     STOPS + "Sums.test_a_5_nov_genome_table_other_than_its_sha256_stops"),
    ("P53", "sha256", "the lab records are read before the sha256 checks", S,
     "    sums = check_sums(args)\n",
     "    read_records(args.records, drug_names(args.drug_names[0]), set())\n    sums = check_sums(args)\n",
     STOPS + "Sums.test_every_sha256_is_checked_before_any_other_read"),
    ("P54", "seal", "a row of the second genome copy inserted at the seal passes", S,
     "        if head is None or head >= SEAL:\n",
     "        if head is None:\n",
     STOPS + "Tables.test_a_genome_table_row_inserted_at_the_seal_stops"),
]

def _c9():
    spec = importlib.util.spec_from_file_location("c9_faults_for_matched", HERE.parent / "c9-pick" / "planted_faults.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verdict(log):
    """unittest's last line in the log: (failures, errors), or None when it gives neither OK nor FAILED."""
    lines = re.findall(r"^(?:OK|FAILED)\b.*$", Path(log).read_text(errors="replace"), re.M)
    if not lines:
        return None
    line = lines[-1]
    failures = re.search(r"failures=(\d+)", line)
    errors = re.search(r"errors=(\d+)", line)
    return int(failures.group(1)) if failures else 0, int(errors.group(1)) if errors else 0


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
            counts = verdict(out / "logs" / f"{n}.log")
            if not ran:
                result = "did not run"
            elif passed:
                result = "passed its test"
            elif counts and counts[0] > 0 and counts[1] == 0:
                result = "failed by an assertion"
            else:
                result = "failed by an error"
        bad += result != "failed by an assertion"
        rows.append((n, rule, what, path, test, result))
        print(f"{n}: {result}: {test}", flush=True)
    with open(out / "faults.tsv", "w") as f:
        f.write("fault\trule\tplants\tfile\ttest\tresult\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"{sum(r[5] == 'failed by an assertion' for r in rows)} of {len(rows)} faults failed their test by an "
          f"assertion")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
