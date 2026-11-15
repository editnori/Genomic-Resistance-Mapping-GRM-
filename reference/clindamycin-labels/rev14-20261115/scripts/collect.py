#!/usr/bin/env python3
"""Gather the scored datasets into the report's tables (PLAN.md sections 6 and 7).

Reads results/scores/<d>/{views,effects,caveat}.tsv and summary.json for every scored dataset, and
results/chosen.tsv, and writes:

  results/refit-table.tsv   per dataset, split and row: V1 (released cohort) and V3 (without the study), recall
                            beside log loss and balanced accuracy, the total effect with its lineage-cluster interval
  results/caveat.tsv        per dataset: the existence and dependence parts of the rule, rows passing of rows tested

Run: python3 -I scripts/collect.py
"""
import csv
import json
import os

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
ORDER = ["scm", "cart", "rf", "gbt", "lasso", "catalog"]


def tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    chosen = {r["dataset"]: r for r in tsv(os.path.join(OUT, "chosen.tsv"))}
    scored = sorted(d for d in os.listdir(os.path.join(OUT, "scores")) if d in chosen)
    head = ["dataset", "split", "row", "n_with", "n_without", "recall_with", "recall_with_cp", "recall_without",
            "recall_without_cp", "recall_total", "recall_total_ci", "log_loss_with", "log_loss_without",
            "log_loss_total", "log_loss_total_ci", "ba_with", "ba_with_ci", "ba_without", "ba_without_ci",
            "ba_total", "ba_total_ci", "ba_composition", "ba_composition_ci", "ba_training", "ba_training_ci",
            "scores_binary"]
    rows = []
    for d in scored:
        views = {(r["split"], r["row"], r["view"]): r for r in tsv(os.path.join(OUT, "scores", d, "views.tsv"))}
        eff = {(r["split"], r["row"], r["measure"], r["effect"]): r
               for r in tsv(os.path.join(OUT, "scores", d, "effects.tsv"))}
        splits = sorted({k[0] for k in views}, key=lambda s: (s != "random", s))
        for sp in splits:
            for lr in ORDER:
                v1, v3 = views[(sp, lr, "V1")], views[(sp, lr, "V3")]

                def e(m, kind):
                    x = eff[(sp, lr, m, kind)]
                    return x["estimate"], f"[{x['cluster_low']}, {x['cluster_high']}]"
                rt, rtc = e("recall", "total")
                lt, ltc = e("log_loss", "total")
                bt, btc = e("balanced_accuracy", "total")
                bc, bcc = e("balanced_accuracy", "composition")
                btr, btrc = e("balanced_accuracy", "training")
                rows.append([d, sp, lr, v1["n"], v3["n"], v1["recall"],
                             f"[{v1['recall_cp_low']}, {v1['recall_cp_high']}]", v3["recall"],
                             f"[{v3['recall_cp_low']}, {v3['recall_cp_high']}]", rt, rtc, v1["log_loss"],
                             v3["log_loss"], lt, ltc, v1["balanced_accuracy"],
                             f"[{v1['ba_cluster_low']}, {v1['ba_cluster_high']}]", v3["balanced_accuracy"],
                             f"[{v3['ba_cluster_low']}, {v3['ba_cluster_high']}]", bt, btc, bc, bcc, btr, btrc,
                             v1["scores_binary"]])
    with open(os.path.join(OUT, "refit-table.tsv"), "w") as fh:
        fh.write("\t".join(head) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    with open(os.path.join(OUT, "caveat.tsv"), "w") as fh:
        fh.write("dataset\tchosen_pmid\texistence\tdependence\trows_passing\trows_tested\tcaveat\tnote\n")
        for d, c in sorted(chosen.items()):
            existence = c["refit"].startswith("no")
            if d in scored:
                s = json.load(open(os.path.join(OUT, "scores", d, "summary.json")))
                dep, npass, ntest = s["dependence_caveat"], s["rows_passing"], s["rows_tested"]
                note = f"lineage clusters {s['lineage_clusters']}; placebo arms {s['placebo_arms']}"
            else:
                dep, npass, ntest = None, "NA", "NA"
                note = (f"not refit: removal leaves reported R {c['reported_R_left']}, S {c['reported_S_left']}, "
                        f"under min_class {c['min_class']}")
            cav = existence or bool(dep)
            fh.write("\t".join([d, c["chosen_pmid"], "yes" if existence else "no",
                                "NA" if dep is None else ("yes" if dep else "no"), str(npass), str(ntest),
                                "yes" if cav else "no", note]) + "\n")
    print(f"scored: {scored}")
    print(open(os.path.join(OUT, "caveat.tsv")).read())


if __name__ == "__main__":
    main()
