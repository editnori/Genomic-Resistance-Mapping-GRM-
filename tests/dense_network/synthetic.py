"""SYNTHETIC inputs for scripts/dense-network/record.py: a made-up experiment and review directory.

make(root) writes both under root and returns their paths. Four made-up datasets in two made-up species,
two in the largest-gap group and two drawn. On every dataset and both splits, each comparator's measure
is the network's minus a fixed offset (DELTA), so every difference over datasets is that offset, with an
interval of width 0, and each derived cell can be worked by hand:

- balanced accuracy: rf is the best of the five on every dataset (the network is 0.08 below it), and the
  site's best is cart on both splits, so it is never the best of the five;
- the lineage baseline is 0.20 below the network, so rule 2 compares 0.20 with 0.28;
- the catalog rule is scored on 3 of the 4 datasets.

lift=True turns every offset of the five around, so the network is 0.05 above the best of the five and
the rule adopts it. None of it is data.
"""

import hashlib
import json
import pathlib
import statistics

from tests.dense_network.common import load_script

record = load_script()

PLAN_TIME = "2026-11-15T12:00:00Z"
FITTING = "12:02:00Z"
FIRST = "12:05:00Z"
LOOK = "12:10:00Z"
LAST = "12:30:00Z"

DATASETS = [  # dataset, species, group, genomes, resistant, gap
    ("genus_alpha__drugone", "genus_alpha", "largest_gap", 40, 20, 0.30),
    ("genus_beta__drugtwo", "genus_beta", "largest_gap", 60, 15, 0.25),
    ("genus_alpha__drugthree", "genus_alpha", "drawn", 50, 25, 0.10),
    ("genus_beta__drugfour", "genus_beta", "drawn", 80, 40, 0.06),
]
NETWORK = {  # balanced accuracy, recall, log loss
    "genus_alpha__drugone": (0.70, 0.60, 0.60),
    "genus_beta__drugtwo": (0.80, 0.70, 0.50),
    "genus_alpha__drugthree": (0.60, 0.50, 0.70),
    "genus_beta__drugfour": (0.75, 0.65, 0.40),
}
# The network minus each comparator, per measure: balanced accuracy, recall, log loss.
DELTA = {
    "scm": (0.02, 0.0, -0.5), "cart": (-0.04, -0.05, 0.1), "rf": (-0.08, -0.10, 0.2),
    "gbt": (-0.02, 0.05, 0.0), "lasso": (0.04, 0.10, -0.1),
    "mean4": (-0.03, -0.01, 0.05), "baseline": (0.20, None, None), "catalog": (0.10, 0.02, -0.2),
}
LIFT = {"scm": (0.09, 0.0, -0.5), "cart": (0.07, -0.005, 0.1), "rf": (0.05, 0.0, 0.2),
        "gbt": (0.11, 0.05, 0.0), "lasso": (0.13, 0.10, -0.1)}
NO_CATALOG = "genus_beta__drugfour"
SITE_BEST = "cart"
KEYS = ["ba", "recall", "log_loss"]


def plan_text():
    """A plan that gives every setting in the words the script looks for."""
    lines = ["# Plan (SYNTHETIC)", ""]
    for key, wanted in record.SETTINGS_IN_PLAN.items():
        for text in (wanted,) if isinstance(wanted, str) else wanted:
            lines.append("- %s: %s" % (key, text))
    return "\n".join(lines) + "\n"


def delta(lift):
    return dict(DELTA, **LIFT) if lift else dict(DELTA)


def learners(ds, lift):
    """{learner: {ba, recall, log_loss[, learner]}} for one dataset and split."""
    d = delta(lift)
    net = dict(zip(KEYS, NETWORK[ds]))
    out = {"candidate": dict(net)}
    for who, offsets in d.items():
        if who == "catalog" and ds == NO_CATALOG:
            continue
        out[who] = {k: net[k] - o for k, o in zip(KEYS, offsets) if o is not None}
    five = record.FIVE
    best = max(five, key=lambda w: out[w]["ba"])
    out["best5"] = dict(out[best], learner=best)
    out["site_best"] = dict(out[SITE_BEST], learner=SITE_BEST)
    return out


def summary(lift=False):
    """The experiment's summary file: rows per dataset, and the figures it derives from them."""
    rows = []
    for ds, species, group, genomes, resistant, _ in DATASETS:
        per = learners(ds, lift)
        rows.append({"id": ds, "species": species, "n": genomes, "n_resistant": resistant, "gate_b": True,
                     "lineage": per, "random": json.loads(json.dumps(per))})
    paired, medians = {}, {}
    for split in record.SPLITS:
        paired[split], medians[split] = {}, {}
        for who in record.COMPARATORS:
            paired[split][who], medians[split][who] = {}, {}
            for k in KEYS:
                have = [r for r in rows if who in r[split] and k in r[split][who]]
                diffs = [r[split]["candidate"][k] - r[split][who][k] for r in have]
                if not diffs:
                    continue
                med = statistics.median(diffs)
                paired[split][who][k] = {"median": med, "lo": min(diffs), "hi": max(diffs), "n": len(diffs),
                                         "species": len({r["species"] for r in have}),
                                         "positive": sum(1 for x in diffs if x > 0),
                                         "negative": sum(1 for x in diffs if x < 0)}
                medians[split][who][k] = statistics.median(r[split][who][k] for r in have)
    line = [r["lineage"] for r in rows]
    cost = {"fits": 24, "processes": 4, "cpu_s": 3600.0, "peak_rss_kb": 1000}
    margin_c = statistics.median(x["candidate"]["ba"] - x["baseline"]["ba"] for x in line)
    margin_b = statistics.median(x["best5"]["ba"] - x["baseline"]["ba"] for x in line)
    gain = paired["lineage"]["best5"]["ba"]
    adopt = (gain["median"] >= 0.02 and gain["lo"] > 0 and margin_c >= margin_b
             and paired["lineage"]["best5"]["recall"]["median"] >= -0.01)
    return {"excluded": [], "rows": rows, "cost": cost, "paired": paired, "medians": medians,
            "primary": {"margin_candidate": margin_c, "margin_best5": margin_b, "cpu_hours": cost["cpu_s"] / 3600},
            "outcome": "adopt" if adopt else "no gain established"}


def sha(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def make(root, plan=None, lift=False, stated=None, review=None, change=None):
    """Write the made-up experiment and review directories under root; return (run, review).

    plan replaces the plan's text; stated replaces the sha256 the report gives for it; review replaces the
    review's report; change(summary, stage1, gates) may alter the three result files before they are written.
    """
    root = pathlib.Path(root)
    run, rev = root / "run", root / "review"
    (run / "results").mkdir(parents=True)
    rev.mkdir()
    plan = plan if plan is not None else plan_text()
    (run / "PLAN.md").write_text(plan, encoding="utf-8")
    (run / "REPORT.md").write_text(
        "# Report (SYNTHETIC)\n\n`PLAN.md` closed all six before any candidate result; sha256 `%s` (file time %s; the "
        "first candidate result was written at %s).\n" % (stated or sha(plan), PLAN_TIME, FIRST), encoding="utf-8")
    review = review if review is not None else (
        "# Review (SYNTHETIC)\n\nThe experiment ran its analysis on the outcomes of the first 2 datasets at %s. The "
        "analysis script was last changed at %s.\n\nOrder: plan; `ba2_candidate.py` last change %s; first result.\n"
        % (LOOK, LAST, FITTING))
    (rev / "REPORT.md").write_text(review, encoding="utf-8")
    s = summary(lift)
    stage1 = {"seed": record.SEED, "pool": 113, "top12": [d[0] for d in DATASETS if d[2] == "largest_gap"],
              "drawn12": [d[0] for d in DATASETS if d[2] == "drawn"], "gap": {d[0]: d[5] for d in DATASETS}}
    gates = [{"id": d[0], "gate_a": True, "gate_b": True} for d in DATASETS]
    if change is not None:
        change(s, stage1, gates)
    for name, obj in (("ba2-summary.json", s), ("stage1-datasets.json", stage1), ("gates.json", gates)):
        (run / "results" / name).write_text(json.dumps(obj), encoding="utf-8")
    return str(run), str(rev)


# Cells worked by hand from DELTA (the network minus the comparator, lineage split).
BY_HAND = {
    "best5_balanced_accuracy": -0.08, "best5_recall": -0.10, "best5_log_loss": 0.2,
    "site_best_balanced_accuracy": -0.04, "baseline_balanced_accuracy": 0.20, "catalog_datasets": 3,
    "margin_candidate": 0.20, "margin_best5": 0.28,
    # how far each of the five sits below the best of the five
    "below_best5": {"scm": 0.10, "cart": 0.04, "rf": 0.0, "gbt": 0.06, "lasso": 0.12},
    "network_median_balanced_accuracy": 0.725,  # median of 0.70, 0.80, 0.60, 0.75
}
