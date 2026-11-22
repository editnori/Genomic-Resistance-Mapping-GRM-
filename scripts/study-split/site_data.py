"""site_data.py ROOT [--check FILE]: the study-split summaries the web app bundles.

ROOT is reference/study-split. The script reads aggregate TSV files only
(summary.tsv, datasets.tsv, contrasts.tsv, first.tsv, rules.tsv) and writes,
as JSON, the dataset/species counts and per-learner contrast summaries for the
site's study-split section. It holds no genome row and estimates nothing.

  python3 -I scripts/study-split/site_data.py reference/study-split > web/src/data/study-split.json
  python3 -I scripts/study-split/site_data.py reference/study-split --check web/src/data/study-split.json
"""

import csv
import json
import pathlib
import sys

LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
FAMILIES = ("SR", "SM", "MR")


def load_tsv(path):
    with pathlib.Path(path).open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def species_name(slug):
    genus, _, rest = slug.partition("_")
    return f"{genus.capitalize()} {rest.replace('_', ' ')}"


def contrast_obj(row):
    return {
        "median": float(row["median"]),
        "low": float(row["low"]),
        "high": float(row["high"]),
        "reading": row["reading"],
        "datasets": int(row["datasets"]),
        "worse": int(row["worse"]),
        "equal": int(row["equal"]),
        "better": int(row["better"]),
    }


def build(root):
    root = pathlib.Path(root)
    summary = {r["key"]: r["value"] for r in load_tsv(root / "summary.tsv")}
    datasets = load_tsv(root / "datasets.tsv")
    contrasts = load_tsv(root / "contrasts.tsv")
    first_rows = load_tsv(root / "first.tsv")
    rule_rows = load_tsv(root / "rules.tsv")

    sp_counts = {}
    for d in datasets:
        sp = d["species"]
        sp_counts[sp] = sp_counts.get(sp, 0) + 1
    species_list = [
        {"slug": slug, "name": species_name(slug), "datasets": n}
        for slug, n in sorted(sp_counts.items(), key=lambda kv: (-kv[1], kv[0]))
    ]

    by_key = {
        (r["experiment"], r["learner"], r["family"], r["instance"], r["measure"]): r
        for r in contrasts
    }
    first_index = {(r["learner"], r["contrast"], r["measure"]): r for r in first_rows}
    rules_index = {(r["learner"], r["family"], r["rule"]): r for r in rule_rows}

    first = {}
    second = {}
    vs_lineage = {}
    for learner in LEARNERS:
        measures = ("recall", "specificity", "balanced_accuracy", "log_loss", "error_rate") if learner == "scm" else ("recall", "specificity", "balanced_accuracy", "log_loss")
        f_entry = {m: contrast_obj(by_key[("first", learner, "study", "-", m)]) for m in measures}
        f_rule = first_index[(learner, "study", "rule")]
        f_entry["rule"] = f_rule["reading"]
        f_entry["keep"] = int(f_rule["keep"])
        first[learner] = f_entry

        loss_measure = "error_rate" if learner == "scm" else "log_loss"
        s_entry = {}
        for fam in FAMILIES:
            r_either = rules_index[(learner, fam, "either")]
            r_recall = rules_index[(learner, fam, "recall")]
            r_loss = rules_index[(learner, fam, "loss")]
            s_entry[fam] = {
                "recall": contrast_obj(by_key[("second", learner, fam, "mean", "recall")]),
                "specificity": contrast_obj(by_key[("second", learner, fam, "mean", "specificity")]),
                "loss": contrast_obj(by_key[("second", learner, fam, "mean", loss_measure)]),
                "loss_measure": loss_measure,
                "pairs": int(r_either["pairs"]),
                "cost": int(r_either["cost"]),
                "recall_worse": int(r_recall["cost"]),
                "loss_worse": int(r_loss["cost"]),
                "answer": r_either["answer"],
            }
        second[learner] = s_entry

        lin_entry = {}
        for m in ("recall", "specificity", "log_loss"):
            obj = contrast_obj(by_key[("first", learner, "study-lineage", "-", m)])
            obj["keep"] = int(first_index[(learner, "study-lineage", m)]["keep"])
            lin_entry[m] = obj
        vs_lineage[learner] = lin_entry

    out = {
        "source": "reference/study-split",
        "evaluation_revision": 14,
        "datasets": int(summary["datasets"]),
        "species": int(summary["species"]),
        "pairs": int(summary["pairs"]),
        "resistant_pairs": int(summary["resistant_pairs"]),
        "no_study_pairs": int(summary["no_study_pairs"]),
        "distinct_genomes": int(summary["distinct_genomes"]),
        "distinct_studies": int(summary["distinct_studies"]),
        "lineage_datasets": int(summary["lineage_datasets"]),
        "species_counts": species_list,
        "first": first,
        "second": second,
        "vs_lineage": vs_lineage,
    }
    text = json.dumps(out, separators=(",", ":"), ensure_ascii=True) + "\n"
    assert json.loads(text) == out
    return text


def main(argv):
    args = argv[1:]
    if len(args) not in (1, 3) or (len(args) == 3 and args[1] != "--check"):
        raise SystemExit(__doc__)
    text = build(args[0])
    if len(args) == 3:
        if pathlib.Path(args[2]).read_text(encoding="utf-8") == text:
            return 0
        print("the bundled file differs from what the record gives", file=sys.stderr)
        return 1
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
