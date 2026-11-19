#!/usr/bin/env python3
"""Stand-in population tables for a rehearsal of the outside-study kit on cohort
genomes (docs/outside-study-test.md, "Rehearsal").

Writes candidates.tsv, pairs.tsv and datasets.tsv in the frozen tables' columns for a list
of cohort genome ids. Labels are the reports' own `samples[].resistant`. Each
genome is put in a stand-in study (REHEARSAL-A, -B or -C, by its place in id
order) so the by-study and leave-one-study-out tables have something to
group. Where a genome is not a sample of a listed dataset, a pair with call
I (first such dataset) or conflict (second) is added, so the run shows that
those carry no label. Measures nothing: these genomes trained the models.
"""

import argparse
import json
import sys
import types
from pathlib import Path

KIT = Path(__file__).resolve().parent
# score.py compiled from its source, never from a bytecode cache (as run.py does).
S = types.ModuleType("outside_study_score")
S.__file__ = str(KIT / "score.py")
exec(compile((KIT / "score.py").read_bytes(), S.__file__, "exec"), S.__dict__)

CANDIDATE_COLS = (
    "genome_id species species_slug salmonella_tag genome_name biosample sra_accession assembly_accession "
    "strain bioproject publication label_pmids collection_date collection_year isolation_country host_name "
    "isolation_source mlst date_inserted_genome_record genome_status genome_quality contigs checkm_completeness "
    "checkm_contamination genome_length species_median_length_held length_ratio older_lab_record qc_verdict "
    "qc_verdict_cohort_median cross_dup_fields cross_dup_matched_ids cross_dup_matches_cohort "
    "repeats_older_pool_sample within_batch_group within_batch_shared within_batch_kept n_labelled_datasets "
    "labelled_datasets same_study_cohort study_keys_in_cohort study_in_held_set_only population "
    "eligible_datasets on_eligible_dataset drop_stage drop_reason"
).split()
PAIR_COLS = ("genome_id dataset_id call records_R records_S records_I records_no_call genome_population "
             "genome_drop_reason dataset_eligible_in_population in_frozen_test in_same_study_set").split()
DATASET_COLS = ("dataset_id species_slug antibiotic cohort_R cohort_S eligible_primary eligible_same_study").split()


def build(ids, datasets, reports_dir, pins):
    reports = {d: S.load_report(reports_dir, d, pins)[0] for d in datasets}
    species = reports[datasets[0]]["species"]
    slug = reports[datasets[0]]["species_slug"]
    meta = {}
    for d in datasets:
        for s in reports[d]["samples"]:
            meta.setdefault(s["id"], s)
    stray = [g for g in ids if g not in meta]
    if stray:
        S.refuse(f"{stray[0]} is not a sample of the listed reports")
    order = sorted(ids, key=S.genome_key)
    cands, pairs = [], []
    for i, g in enumerate(order):
        labelled = []
        odd = 0
        for d in datasets:
            sample = next((s for s in reports[d]["samples"] if s["id"] == g), None)
            if sample is None:
                if odd < 2:
                    call = "I" if odd == 0 else "conflict"
                    pairs.append({"genome_id": g, "dataset_id": d, "call": call,
                                  "records_R": 1 if call == "conflict" else 0,
                                  "records_S": 1 if call == "conflict" else 0,
                                  "records_I": 1 if call == "I" else 0, "records_no_call": 0,
                                  "genome_population": "primary (other-study)", "genome_drop_reason": "",
                                  "dataset_eligible_in_population": "yes", "in_frozen_test": "no",
                                  "in_same_study_set": "no"})
                    odd += 1
                continue
            labelled.append(d)
            call = "R" if sample["resistant"] else "S"
            pairs.append({"genome_id": g, "dataset_id": d, "call": call,
                          "records_R": int(call == "R"), "records_S": int(call == "S"),
                          "records_I": 0, "records_no_call": 0,
                          "genome_population": "primary (other-study)", "genome_drop_reason": "",
                          "dataset_eligible_in_population": "yes", "in_frozen_test": "yes",
                          "in_same_study_set": "no"})
        m = meta[g]
        cands.append({"genome_id": g, "species": species, "species_slug": slug,
                      "bioproject": f"REHEARSAL-{'ABC'[i % 3]}",
                      "isolation_country": m.get("country") or "", "mlst": m.get("mlst") or "",
                      "n_labelled_datasets": len(labelled), "labelled_datasets": ";".join(labelled),
                      "population": "primary (other-study)", "eligible_datasets": ";".join(labelled),
                      "on_eligible_dataset": "yes", "same_study_cohort": "no"})
    drows = [{"dataset_id": d, "species_slug": slug, "antibiotic": reports[d]["antibiotic"],
              "cohort_R": reports[d]["counts"]["resistant"], "cohort_S": reports[d]["counts"]["susceptible"],
              "eligible_primary": "yes", "eligible_same_study": "no"} for d in datasets]
    return cands, pairs, drows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ids", required=True)
    ap.add_argument("--datasets", required=True)
    ap.add_argument("--reports", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)
    try:
        ids = [x.strip() for x in Path(args.ids).read_text().split() if x.strip()]
        datasets = [x.strip() for x in Path(args.datasets).read_text().split() if x.strip()]
        cands, pairs, drows = build(ids, datasets, args.reports, S.read_pins(KIT / "reports.sha256"))
    except S.Refusal as e:
        print(f"rehearsal.py: refused: {e}", file=sys.stderr)
        return 2
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    S.write_tsv(out / "candidates.tsv", CANDIDATE_COLS, cands)
    S.write_tsv(out / "pairs.tsv", PAIR_COLS, pairs)
    S.write_tsv(out / "datasets.tsv", DATASET_COLS, drows)
    print(json.dumps({"genomes": len(cands), "pairs": len(pairs), "datasets": len(drows)}), file=sys.stderr)
    return 0


if __name__ == "__main__":
    S.isolated_or_exit("rehearsal.py")
    sys.exit(main())
