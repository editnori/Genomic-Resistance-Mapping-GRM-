"""SYNTHETIC inputs for scripts/kover-comparison/record.py: a made-up comparison and review directory.

make(root) writes both under root and returns their paths. Three made-up datasets of 20 genomes, 5 folds
of 4 genomes each (2 resistant, 2 susceptible), with the confusion counts chosen so that each cell the
script derives can be worked by hand:

- alpha: in fold 0 Kover calls one susceptible genome resistant (specificity 0.5, balanced accuracy 0.75)
  where GRM calls all 4 right; in folds 1 to 4 both call all 4 right. Kover has one rule in every fold,
  and in fold 0 two kept rules tie, one with Kover's 0.75 and one with 1.0.
- beta: in fold 2 Kover misses one resistant genome (recall 0.5); elsewhere both are right. Kover has
  two rules in every fold, GRM one, so no fold has a tie row.
- gamma: both miss the same resistant genome in every fold (balanced accuracy 0.75 on both sides), and
  Kover's one rule ties with a second rule of the same test result.

None of it is data.
"""

import hashlib
import json
import math
import pathlib

from tests.kover_comparison.common import load_script

record = load_script()

PLAN_TIME = "2026-11-15T12:00:00Z"
FIT_TIME = "2026-11-15T12:10:00Z"
REVIEW_TIME = "2026-11-15T16:00:00Z"
TIES_TIME = "2026-11-15T16:40:00Z"

PERFECT = (2, 0, 2, 0)
# dataset: (species, per fold (grm tp fn tn fp, kover tp fn tn fp, grm rules, kover rules, agreeing), tie rows)
DESIGN = {
    "genus_alpha__drugone": ("Genus alpha", [((2, 0, 2, 0), (2, 0, 1, 1), 1, 1, 3)] + [(PERFECT, PERFECT, 1, 1, 4)] * 4,
                             {0: (3, 2, [(0.75, 1), (1.0, 1)])} | {k: (1, 1, [(1.0, 1)]) for k in range(1, 5)}),
    "genus_beta__drugtwo": ("Genus beta", [(PERFECT, PERFECT, 1, 2, 4)] * 2 + [((2, 0, 2, 0), (1, 1, 2, 0), 1, 2, 3)]
                            + [(PERFECT, PERFECT, 1, 2, 4)] * 2, {}),
    "genus_gamma__drugthree": ("Genus gamma", [((1, 1, 2, 0), (1, 1, 2, 0), 2, 1, 4)] * 5,
                               {k: (2, 2, [(0.75, 2)]) for k in range(5)}),
}


def tsv(head, rows):
    return "\t".join(head) + "\n" + "".join("\t".join(str(c) for c in r) + "\n" for r in rows)


def plan_text():
    """A plan that gives every setting in the words the script looks for, and the settings table."""
    lines = ["# Plan (SYNTHETIC)", ""]
    for key, wanted in record.SETTINGS_IN_PLAN.items():
        for text in (wanted,) if isinstance(wanted, str) else wanted:
            lines.append("- %s: %s" % (key, text))
    lines += ["", "| Setting | GRM | Kover | Same |", "|---|---|---|---|"]
    for setting, grm, kover, same in record.MAPPING:
        grm = record.MAPPING_REWORDED.get(setting, (grm, None))[0]
        lines.append("| %s | `%s` | %s | %s |" % (setting, grm, kover, same))
    return "\n".join(lines) + "\n"


def sha(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def measures(tp, fn, tn, fp):
    return record.measures(tp, fn, tn, fp)


def make(root, plan=None, datasets_after_hash=None):
    """Write the made-up comparison and review directories under root; return (kover, review).

    plan replaces the plan's text; datasets_after_hash, if given, is written over datasets.tsv after
    the plan hashed it.
    """
    root = pathlib.Path(root)
    k, r = root / "kover", root / "review"
    (k / "results" / "review").mkdir(parents=True)
    (k / "scripts").mkdir()
    (r / "results").mkdir(parents=True)

    rows = []
    for i, (ds, (species, folds, _)) in enumerate(DESIGN.items(), 1):
        rows.append([i, ds, species, 20, 10, 10, 7, 5, 20, "4,4,4,4,4"])
    datasets = tsv(record.COLUMNS["datasets.tsv"], rows)
    plan = plan if plan is not None else plan_text()
    files = {"PLAN.md": plan, "results/datasets.tsv": datasets,
             "scripts/select_datasets.py": "# SYNTHETIC\n", "scripts/build_kmer_tools.sh": "# SYNTHETIC\n",
             "scripts/smoke_synthetic.sh": "# SYNTHETIC\n"}
    for name, text in files.items():
        (k / name).write_text(text, encoding="utf-8")
    (k / "PLAN.sha256").write_text(
        "# sha256 of PLAN.md and of the inputs and scripts it names, recorded at %s, before any Kover fit\n" % PLAN_TIME
        + "".join("%s  %s\n" % (sha(files[name]), name) for name in record.PLAN_FILES), encoding="utf-8")
    if datasets_after_hash is not None:
        (k / "results" / "datasets.tsv").write_text(datasets_after_hash, encoding="utf-8")
    (k / "results" / "run_all.log").write_text("run_all start %s orders 1..3\n" % FIT_TIME, encoding="utf-8")

    head = ["dataset", "fold", "n_test", "resistant", "susceptible", "agreement"]
    for m in ("balanced_accuracy", "recall", "specificity", "log_loss"):
        head += ["grm_" + m, "kover_" + m]
    head += ["grm_rules", "kover_rules", "kover_kind", "kover_p"]
    per_fold = []
    for ds, (_, folds, ties) in DESIGN.items():
        for fold, (g, kv, g_rules, k_rules, agree) in enumerate(folds):
            mg, mk = measures(*g), measures(*kv)
            row = [ds, fold, 4, 2, 2, "%.6f" % (agree / 4)]
            for m in ("balanced_accuracy", "recall", "specificity", "log_loss"):
                row += ["%.6f" % mg[m], "%.6f" % mk[m]]
            row += ["%.1f" % g_rules, k_rules, "disjunction" if ds.startswith("genus_gamma") else "conjunction", "1.0"]
            per_fold.append(row)
        tie_folds = []
        for fold in range(5):
            if fold not in ties:
                tie_folds.append({"fold": fold, "skipped": "Kover model has 2 rules"})
                continue
            tied, kept, values = ties[fold]
            tie_folds.append({"fold": fold, "rules_tied_at_best_utility": tied, "rules_after_risk_tiebreak": kept,
                              "distinct_training_patterns_among_kept": kept,
                              "first_index_reproduces_kover_calls": True,
                              "first_index_test_ba": measures(*folds[fold][1])["balanced_accuracy"],
                              "test_ba_of_kept_rules": [[b, c] for b, c in values]})
        (k / "results" / "review" / ("fold_ties_%s.json" % ds)).write_text(
            json.dumps({"dataset": ds, "folds": tie_folds}), encoding="utf-8")
    (k / "results" / "per_fold.tsv").write_text(tsv(head, per_fold), encoding="utf-8")

    check_head = ["dataset", "kind", "p", "step", "grm_covered", "grm_lost", "kover_covered", "kover_lost",
                  "counts_equal", "pattern_equal"]
    check = [["genus_alpha__drugone", "conjunction", "1.0", 1, 9, 1, 9, 1, "True", "True"],
             ["genus_alpha__drugone", "conjunction", "1.0", 2, 3, 0, 3, 0, "True", "False"],
             ["genus_gamma__drugthree", "disjunction", "0.5", 1, 8, 2, 8, 2, "True", "False"]]
    tied_head = ["dataset", "step", "grm_listed_tied", "kover_equals_grm_chosen", "kover_equals_a_listed_tied"]
    tied = [["genus_alpha__drugone", 1, 1, "True", "False"], ["genus_alpha__drugone", 2, 4, "False", "True"],
            ["genus_gamma__drugthree", 1, 3, "False", "True"]]
    (k / "results" / "final_check.tsv").write_text(tsv(check_head, check), encoding="utf-8")
    (k / "results" / "final_tied.tsv").write_text(tsv(tied_head, tied), encoding="utf-8")

    (k / "results" / "review" / "fold_ties_summary_exploratory.txt").write_text(
        "%s\nSYNTHETIC tie summary\n" % TIES_TIME, encoding="utf-8")
    (k / "results" / "review" / "imbalance.txt").write_text("\n".join([
        "chosen12 datasets 3 | resistant share median 0.500 range 0.5-0.5 | resistant-majority 0 | lineages median 7.0",
        "eligible datasets 9 | resistant share median 0.400 range 0.1-0.9 | resistant-majority 2 | lineages median 12.5",
        "genome-dataset entries: eligible 900, chosen 60 (6.67%)",
        "n median: eligible 80.0, chosen 20.0",
        "species: eligible 5, with a chosen dataset 3, without 2",
        "GRM lineage scm balanced accuracy median: chosen 0.900, eligible 0.850",
    ]) + "\n", encoding="utf-8")

    review_plan = "# Review plan (SYNTHETIC)\n"
    (r / "PLAN.md").write_text(review_plan, encoding="utf-8")
    (r / "PLAN.sha256").write_text("# sha256 of PLAN.md, recorded at %s\n%s  PLAN.md\n" % (REVIEW_TIME, sha(review_plan)),
                                   encoding="utf-8")
    return str(k), str(r)


# Cells worked by hand from DESIGN (Kover minus GRM).
CLIP_LOSS = -math.log(1e-4)
RIGHT_LOSS = -math.log(1 - 1e-4)
ONE_ERROR_IN_4 = (CLIP_LOSS + 3 * RIGHT_LOSS) / 4 - RIGHT_LOSS
BY_HAND = {
    # per dataset: balanced accuracy, recall, specificity, error rate, rules, agreeing share
    "genus_alpha__drugone": (-0.25 / 5, 0.0, -0.5 / 5, 0.25 / 5, 0.0, 19 / 20),
    "genus_beta__drugtwo": (-0.25 / 5, -0.5 / 5, 0.0, 0.25 / 5, 1.0, 19 / 20),
    "genus_gamma__drugthree": (0.0, 0.0, 0.0, 0.0, -1.0, 1.0),
}
