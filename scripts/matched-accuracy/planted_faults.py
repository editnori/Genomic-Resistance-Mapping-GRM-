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
     "        if max(heads) >= SEAL:\n",
     "        if max(heads) > SEAL:\n",
     STOPS + "Records.test_a_record_inserted_at_the_seal_stops"),
    ("P06", "seal", "a lab record with a blank insertion time passes", S,
     "        if heads is None:\n            raise Stop(\"a lab record's date_inserted",
     "        if heads is None:\n            heads = [\"\"]\n        if heads is None:\n"
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
     '                    c["bioproject_" + why[s["id"]]] += 1\n',
     '                    c["bioproject_" + why[s["id"]]] += why[s["id"]] == "no_accession"\n',
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
     "        if any(got.get(name, 0) <= found.get(name, 0) for name in files):\n",
     "        if any(got.get(name, 0) < found.get(name, 0) for name in files):\n",
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
     "    ids = every_genome | record_ids | label_ids | genome_ids\n",
     "    ids = set(every_genome)\n",
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
     "        if heads is None or max(heads) >= SEAL:\n",
     "        if heads is None:\n",
     STOPS + "Tables.test_a_genome_table_row_inserted_at_the_seal_stops"),
    ("P55", "seal", "a lab record's list of insertion times is read as its joined text", S,
     '        heads = stamp_heads(rec.get("date_inserted"))\n',
     '        heads = stamp_heads(value_text(rec.get("date_inserted")))\n',
     STOPS + "Records.test_a_list_of_insertion_times_with_one_at_the_seal_stops"),
    ("P56", "seal", "a blank or malformed text of a list of insertion times is passed over", S,
     "        return None if not heads or None in heads else [h for hs in heads for h in hs]\n",
     "        return None if all(h is None for h in heads) else [h for hs in heads if hs for h in hs]\n",
     STOPS + "Records.test_a_list_of_insertion_times_with_one_blank_or_malformed_stops"),
    ("P57", "seal", "a second copy row's list of insertion times is read as its joined text", S,
     '        heads = stamp_heads(row.get("date_inserted"))\n',
     '        heads = stamp_heads(value_text(row.get("date_inserted")))\n',
     STOPS + "Tables.test_a_genome_table_row_with_a_list_of_times_one_at_the_seal_stops"),
    ("P58", "seal", "an insertion time is not read as a date, so month 13 passes", S,
     '        datetime.datetime.strptime(head, "%Y-%m-%dT%H:%M:%S")\n',
     '        head.encode("ascii")\n',
     STOPS + "Records.test_an_impossible_insertion_date_stops"),
    ("P59", "seal", "an insertion time in other than ASCII digits is read as a time", S,
     'STAMP = re.compile(r"\\d{4}-\\d\\d-\\d\\dT\\d\\d:\\d\\d:\\d\\d", re.ASCII)',
     'STAMP = re.compile(r"\\d{4}-\\d\\d-\\d\\dT\\d\\d:\\d\\d:\\d\\d")',
     STOPS + "Records.test_an_insertion_time_in_other_than_ascii_digits_stops"),
    ("P60", "input shape", "a report that is not JSON ends in a Python error", S,
     '    except ValueError:\n        raise Stop(f"{Path(path).name} is not a JSON report") from None\n',
     '    except KeyError:\n        raise Stop(f"{Path(path).name} is not a JSON report") from None\n',
     STOPS + "Inputs.test_a_report_that_is_not_json_stops"),
    ("P61", "input shape", "a drug names file without a list of names ends in a Python error", S,
     "    if not isinstance(names, list) or not all(\n",
     "    if isinstance(names, str) and not all(\n",
     STOPS + "Inputs.test_a_drug_names_file_without_a_list_of_names_stops"),
    ("P62", "input shape", "a label row whose antibiotic is a list ends in a Python error", S,
     "        if antibiotic is not None and not isinstance(antibiotic, str):\n",
     "        if isinstance(antibiotic, int):\n",
     STOPS + "Inputs.test_a_label_row_whose_antibiotic_is_not_text_stops"),
    ("P63", "input shape", "a null genome row of a report left out ends in a Python error", S,
     '            if not all(s is not None and isinstance(s["id"], str) and s["id"].strip() for s in rep["samples"]):\n',
     '            if not all(s is None or isinstance(s["id"], str) for s in rep["samples"]):\n',
     STOPS + "Inputs.test_a_genome_row_of_a_report_left_out_without_a_text_id_stops"),
    ("P64", "input shape", "an output path that is a file is not checked before the write", S,
     "    if out.exists() and not out.is_dir():\n",
     "    if out.is_dir() and not out.exists():\n",
     STOPS + "Inputs.test_an_output_path_that_is_a_file_stops_and_leaves_the_file"),
    ("P65", "write", "an OSError of the write ends in a Python error", S,
     "    except OSError as e:\n        left = 0\n",
     "    except KeyError as e:\n        left = 0\n",
     STOPS + "Inputs.test_a_umask_that_bars_the_output_stops_and_leaves_nothing"),
    ("P66", "write", "a write that fails removes nothing it made", S,
     "                m.rmdir() if m.is_dir() else m.unlink()\n",
     "                m.is_dir()\n",
     STOPS + "Inputs.test_a_write_that_fails_after_a_file_removes_what_it_made"),
    ("P67", "BioProject", "a value of one part without an accession's shape reads as a BioProject", S,
     "        project[g] = part if part is not None and is_accession(part) else None\n",
     "        project[g] = part\n",
     ROWS + "Genomes.test_a_value_of_one_part_without_an_accessions_shape_reads_as_no_bioproject"),
    ("P68", "BioProject", "a BioProject accession is matched with case", S,
     'BIOPROJECT = re.compile(r"PRJ[A-Z]{2}[0-9]+", re.IGNORECASE | re.ASCII)',
     'BIOPROJECT = re.compile(r"PRJ[A-Z]{2}[0-9]+", re.ASCII)',
     ROWS + "Genomes.test_a_value_of_one_part_without_an_accessions_shape_reads_as_no_bioproject"),
    ("P69", "BioProject", "a part that holds an accession's shape inside it reads as an accession", S,
     "    return BIOPROJECT.fullmatch(part) is not None\n",
     "    return BIOPROJECT.search(part) is not None\n",
     ROWS + "Genomes.test_a_value_of_one_part_without_an_accessions_shape_reads_as_no_bioproject"),
    ("P70", "BioProject", "a genome whose one part is no accession is not counted apart", S,
     '        why[g] = "not_in_table" if g not in ids else "no_accession" if part is not None and project[g] is None \\\n',
     '        why[g] = "not_in_table" if g not in ids else None if part is not None and project[g] is None \\\n',
     ROWS + "Genomes.test_a_value_of_one_part_without_an_accessions_shape_reads_as_no_bioproject"),
    ("P71", "scan", "the scan looks up every BioProject part, a dash too", S,
     "        accessions |= {p for p in parts if is_accession(p)}\n",
     "        accessions |= parts\n",
     ROWS + "Genomes.test_a_value_of_one_part_without_an_accessions_shape_reads_as_no_bioproject"),
    ("P72", "BioProject", "a blank row beside a row naming a BioProject gives the genome none", S,
     "        if genome in wanted and parts:\n",
     "        if genome in wanted:\n",
     ROWS + "Genomes.test_a_blank_row_beside_a_row_naming_a_bioproject_is_skipped"),
    ("P73", "BioProject", "the genomes the second copy lacks are counted over the reports that enter alone", S,
     '             "not_in_copy_every": len(every_genome - genome_ids), ',
     '             "not_in_copy_every": len(wanted_genomes - genome_ids), ',
     MADE_UP + "test_the_inputs_and_checks_equal_the_hand_count"),
    ("P74", "typing", "the rows lacking a value among the eligible rows count every row entering the groups", S,
     '                    if (r["dataset"], r["genome"]) in eligible:\n',
     '                    if r["dataset"]:\n',
     ROWS + "Typing.test_a_row_lacking_a_value_is_counted_among_the_rows_entering_groups_and_among_the_eligible_rows"),
    ("P75", "typing", "the two counts of rows lacking a value swap columns", S,
     "lacking_value_of_rows_entering_groups=lacking[level][n],",
     "lacking_value_of_rows_entering_groups=lacking_eligible[level][n],",
     ROWS + "Typing.test_a_row_lacking_a_value_is_counted_among_the_rows_entering_groups_and_among_the_eligible_rows"),
    ("P76", "export", "the 5 Nov table's note says nothing of its insertion time", S,
     'f"rows; the table has no insertion time, so no row of it is checked against "',
     'f"rows; no row of it is checked against "',
     MADE_UP + "test_the_inputs_and_checks_equal_the_hand_count"),
    ("P77", "scan", "the positive control reads the first export file alone", S,
     "        if any(got.get(name, 0) <= found.get(name, 0) for name in files):\n",
     "        if got.get(min(files), 0) <= found.get(min(files), 0):\n",
     SCAN + "Scan.test_a_scan_cut_to_the_first_file_fails_the_control"),
    ("P78", "scan", "the positive control plants no token of each shape", S,
     "    return plants + [(what, token, ids) for what, token in CONTROL_TOKENS]\n",
     "    return plants\n",
     SCAN + "Scan.test_a_scan_cut_to_one_shape_fails_the_control"),
    ("P79", "scan", "the positive control plants a read genome id that may have a shape", S,
     "    shapeless = sorted(i for i in ids if not shaped(i))\n",
     "    shapeless = sorted(ids)\n",
     SCAN + "Scan.test_a_scan_that_looks_up_no_read_id_fails_the_control"),
    ("P80", "scan", "the positive control plants nothing when every read id has a shape", S,
     "        plants = [(\"made-up genome id without an id's shape\", CONTROL_ID, ids | {CONTROL_ID})]\n",
     "        plants = []\n",
     SCAN + "Scan.test_when_every_read_id_has_a_shape_a_made_up_id_stands_in"),
    ("P81", "scan", "the accession shapes match with case", S,
     'ACCESSION_SHAPED = re.compile(r"(?:PRJ[A-Z]{1,2}|SAM[A-Z]{1,2}|GC[AF]_|[SED]R[RSXP])\\d+", re.IGNORECASE)',
     'ACCESSION_SHAPED = re.compile(r"(?:PRJ[A-Z]{1,2}|SAM[A-Z]{1,2}|GC[AF]_|[SED]R[RSXP])\\d+")',
     SCAN + "Hits.test_an_accession_without_case_inside_a_word_or_of_a_run_or_sample_is_a_hit"),
    ("P82", "scan", "the accession shapes start at a word boundary", S,
     'ACCESSION_SHAPED = re.compile(r"(?:PRJ',
     'ACCESSION_SHAPED = re.compile(r"\\b(?:PRJ',
     SCAN + "Hits.test_an_accession_without_case_inside_a_word_or_of_a_run_or_sample_is_a_hit"),
    ("P83", "scan", "the accession shapes hold no run or sample accession", S,
     "|GC[AF]_|[SED]R[RSXP])",
     "|GC[AF]_)",
     SCAN + "Hits.test_an_accession_without_case_inside_a_word_or_of_a_run_or_sample_is_a_hit"),
    ("P84", "seal", "a row of the second copy with a blank insertion time passes", S,
     "        if heads is None or max(heads) >= SEAL:\n",
     "        if heads is not None and max(heads) >= SEAL:\n",
     STOPS + "Tables.test_a_genome_table_row_with_a_blank_insertion_time_stops"),
    ("P85", "seal", "an insertion time with a space for the T passes", S,
     '    if not STAMP.fullmatch(head):\n        return None\n    try:\n'
     '        datetime.datetime.strptime(head, "%Y-%m-%dT%H:%M:%S")\n',
     '    if not STAMP.fullmatch(head.replace(" ", "T", 1)):\n        return None\n    try:\n'
     '        datetime.datetime.strptime(head.replace(" ", "T", 1), "%Y-%m-%dT%H:%M:%S")\n',
     STOPS + "Records.test_an_insertion_time_with_a_space_for_the_t_stops"),
    ("P86", "input shape", "a genome id repeated in a report passes", S,
     'if len({s["id"] for s in rows}) != n:',
     'if len({s["id"] for s in rows}) > n:',
     STOPS + "Reports.test_a_genome_id_repeated_in_a_report_stops"),
    ("P87", "copy group", "a genome in two copy groups passes", S,
     "if len(members) != len(set(members)):",
     "if len(members) < len(set(members)):",
     STOPS + "Reports.test_a_genome_in_two_copy_groups_stops"),
    ("P88", "learners", "a learner name outside the six passes", S,
     "or not set(names) <= set(LEARNER_NAMES):",
     "or False:",
     STOPS + "Reports.test_a_learner_name_outside_the_six_stops"),
    ("P89", "scan", "the scan's ids leave out the 5 Nov table's", S,
     "    if hc20_ids:\n        ids |= hc20_ids\n",
     "",
     SCAN + "RunWithAHit.test_a_genome_id_read_from_the_5_nov_table_alone_is_a_hit"),
    ("P90", "drug names", "a record name of a kind other than drug reaches its pair", S,
     'if kind != "drug" or (genome, canonical) not in wanted:',
     "if (genome, canonical) not in wanted:",
     T + "test_made_up.OneChange.test_a_record_name_of_another_kind_reaches_no_pair_though_its_canonical_name_is_a_drug"),
    ("P91", "folds", "a test_fold value true passes as fold 1", S,
     "not all(is_int(f) for f in fold)",
     "not all(isinstance(f, int) for f in fold)",
     STOPS + "Reports.test_a_test_fold_value_of_true_stops"),
    ("P92", "copy group", "a copy group index past the last row is not stopped", S,
     "all(is_int(i) and 0 <= i < n for i in g)",
     "all(is_int(i) and 0 <= i for i in g)",
     STOPS + "Reports.test_a_copy_group_index_past_the_last_row_stops"),
    ("P93", "lab call", "a label table pair given twice passes, the last standing", S,
     "        if key in table:\n",
     "        if False:\n",
     STOPS + "Calls.test_a_label_table_pair_given_twice_stops"),
    ("P94", "tuple", "two methods in one record's list are joined into one method", S,
     'methods = values(rec.get("laboratory_typing_method"))',
     'methods = {value_text(rec.get("laboratory_typing_method"))} - {None}',
     ROWS + "Records.test_a_list_of_two_methods_in_one_record_gives_two_tuples"),
    ("P95", "shares", "among BioProjects tied on rows the share gives the smallest weight", S,
     "its = max(weight[b] for b in rows if rows[b] == top)",
     "its = min(weight[b] for b in rows if rows[b] == top)",
     ROWS + "Weights.test_among_bioprojects_tied_on_rows_the_share_gives_the_largest_weight"),
    ("P96", "typing", "an HC20 value that differs between rows gives the first", S,
     "if rows and len(set(rows)) == 1 and len(rows[0]) == 1 else None",
     "if rows and len(rows[0]) == 1 else None",
     ROWS + "HC20.test_a_genome_gets_one_hc20_group_or_none"),
    ("P97", "scan", "the scan counts no SAM-shaped token", S,
     "(?:PRJ[A-Z]{1,2}|SAM[A-Z]{1,2}|GC[AF]_",
     "(?:PRJ[A-Z]{1,2}|GC[AF]_",
     SCAN + "Hits.test_a_biosample_accession_is_a_hit"),
    ("P98", "copy group", "the first member of each copy group stays", S,
     'copies = {i for g in rep["kset_duplicates"] or [] for i in g}',
     'copies = {i for g in rep["kset_duplicates"] or [] for i in g[1:]}',
     MADE_UP + "test_every_member_of_the_listed_copy_group_leaves_under_reason_1_in_both_readings"),
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
