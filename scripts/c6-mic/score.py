"""Scores the results of `grm learn --exp-mic` (format grm-mic-plan-result-v1, docs/mic-plan.md).

  python3 -I score.py <expect.tsv> <results dir> <release dir> <out dir>

<expect.tsv> lists every result the run must hold, one per line, with the
columns dataset, species_slug, split (lineage or study), plan (the result's
scheme: lineage-s42, study-whole, ...), seed, arms (released,ordered or
released) and file (the result's path under <results dir>).

First every result is checked; the checks read no score:
  - the file exists (a dataset grm learn skips writes no file and exits 0);
  - its format, evaluation revision, dataset, plan and arms are the ones
    expected, and its input fingerprint equals that of the dataset's report
    in <release dir>/datasets/;
  - each fold's test counts equal the rows it tests;
  - each fold's training median is micbase's B2 recounted from the other
    folds' rows, and every row of the fold carries it, when the plan tests
    every genome with an MIC (else the check says it was not possible);
  - an arm that is available has a prediction for every tested row.
They go to <out dir>/checks.tsv. When a file is missing or a check fails,
nothing else is written and the exit code is 2 or 1.

Then <out dir>/scores.tsv: one row per dataset, split, seed and arm
(released, ordered, median), over the genomes with an MIC that the plan
tests, with each figure's count and denominator:
  - EA-exact: exact MICs whose predicted step is within one dilution, over
    the exact MICs tested. It reads MIC values only.
  - Category agreement: calls equal to the reported call, over the MIC
    genomes tested; the very major error rate (resistant called S) over the
    resistant ones; the major error rate (susceptible called R) over the
    susceptible ones. They read the R and S labels, through the call and
    through the breakpoint each fold fits on its training genomes.
A genome is called R when its predicted step is above its fold's r_above.
The released arm's step is its margin rounded half away from zero, as Rust
rounds. An arm the result marks unavailable, or a fold with MIC test genomes
and no breakpoint, leaves its figures empty, with the folds counted; no
figure is ever taken over part of the folds. The file holds aggregates
only: no genome id and no row per genome.
"""

import csv
import json
import math
import os
import sys

FORMAT = "grm-mic-plan-result-v1"
REVISION = 14
ARMS = ["released", "ordered", "median"]
EXPECT_COLUMNS = ["dataset", "species_slug", "split", "plan", "seed", "arms", "file"]
CHECK_COLUMNS = [
    "dataset", "split", "plan", "seed", "file", "exists", "format", "fingerprint",
    "plan_and_arms", "fold_counts", "median_check", "predictions", "ok",
]
SCORE_COLUMNS = [
    "dataset", "species_slug", "antibiotic", "split", "plan", "seed", "draw", "arm",
    "unavailable", "folds", "folds_scored", "folds_undefined", "folds_without_breakpoint",
    "ea_exact", "ea_exact_count", "ea_exact_n",
    "category_agreement", "category_count", "category_n",
    "very_major_rate", "very_major_count", "very_major_n",
    "major_rate", "major_count", "major_n",
    "ea_exact_reads", "category_reads",
]
EA_READS = "MIC values only"
CATEGORY_READS = "R and S labels and the training breakpoint"


def rround(x: float) -> int:
    """Rust f64::round: nearest integer, half away from zero (micbase)."""
    a = abs(x)
    f = math.floor(a)
    r = f + 1 if a - f >= 0.5 else f  # a - f is exact for finite doubles
    return int(r) if x >= 0 else -int(r)


def agrees(p: int, mic: dict) -> bool:
    """Essential agreement of a predicted step with an MIC (micbase)."""
    v, c = mic["log2"], mic["censor"]
    if c == "Exact":
        return abs(p - v) <= 1
    if c == "AtMost":
        return p <= v + 1
    if c == "Above":
        return p >= v
    raise ValueError(c)


def lower_median(values):
    """micbase's B2: the lower median, the value at position ceil(m/2)."""
    exact = sorted(values)
    if not exact:
        return None
    return exact[math.ceil(len(exact) / 2) - 1]


def step(arm: str, value):
    """An arm's predicted step: the released margin rounded, else as given."""
    if value is None:
        return None
    return rround(value) if arm == "released" else int(value)


def read_expect(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if rows and list(rows[0].keys()) != EXPECT_COLUMNS:
        raise ValueError(f"{path}: columns {list(rows[0].keys())}, want {EXPECT_COLUMNS}")
    return rows


def plan_matches(result, exp):
    plan = result.get("plan") or {}
    if exp["split"] == "lineage":
        ok = plan.get("kind") == "lineage" and exp["plan"] == f"lineage-s{plan.get('seed')}"
        ok = ok and str(plan.get("seed")) == exp["seed"]
    else:
        ok = plan.get("kind") == "fold-file" and plan.get("scheme") == exp["plan"]
    return ok and sorted(result.get("arms") or []) == sorted(exp["arms"].split(","))


def fold_counts_equal(result):
    rows, folds = result["tested"], result["folds"]
    for f, fold in enumerate(folds):
        mine = [r for r in rows if r["fold"] == f]
        if (
            fold["test_mic"] != len(mine)
            or fold["test_exact"] != sum(1 for r in mine if r["censor"] == "Exact")
            or fold["test_resistant"] != sum(1 for r in mine if r["resistant"])
        ):
            return False
    return all(0 <= r["fold"] < len(folds) for r in rows)


def median_check(result):
    """'equal', 'differs' or 'not possible': each fold's training median
    against micbase's B2 on the other folds' tested exact MICs."""
    rows, folds = result["tested"], result["folds"]
    if result.get("unavailable") or not folds:
        return "not possible"
    if sum(f["test_mic"] for f in folds) != result["n_mic"]:
        return "not possible"  # the plan trains genomes it never tests
    for f, fold in enumerate(folds):
        train = [r["log2"] for r in rows if r["fold"] != f and r["censor"] == "Exact"]
        if fold["train_exact"] != len(train) or fold["train_mic"] != sum(1 for r in rows if r["fold"] != f):
            return "differs"
        # A fold that tests no genome with an MIC is not fitted, so it holds no median.
        want = lower_median(train) if fold["test_mic"] > 0 else None
        if fold["median"] != want:
            return "differs"
        mine = [r for r in rows if r["fold"] == f]
        if result.get("median_unavailable") is None and any(r["median"] != fold["median"] for r in mine):
            return "differs"
    return "equal"


def arm_unavailable(result, arm):
    return result.get("unavailable") or result.get(f"{arm}_unavailable")


def predictions_whole(result):
    """Every available arm predicts every tested row, and an unavailable
    arm predicts none."""
    for arm in [a for a in ARMS if a == "median" or a in result["arms"]]:
        have = [r[arm] is not None for r in result["tested"]]
        if arm_unavailable(result, arm):
            if any(have):
                return False
        elif not all(have):
            return False
    return True


def fold_has_fit(fold, arm):
    return {
        "released": fold.get("released_trees") is not None,
        "ordered": fold.get("ordered") is not None,
        "median": fold.get("median") is not None,
    }[arm]


def score_arm(result, exp, arm):
    rows, folds = result["tested"], result["folds"]
    scored = [f for f, fold in enumerate(folds) if fold["test_mic"] > 0]
    why = arm_unavailable(result, arm)
    exact = [r for r in rows if r["censor"] == "Exact"]
    resistant = [r for r in rows if r["resistant"]]
    susceptible = [r for r in rows if not r["resistant"]]
    no_bp = [f for f in scored if folds[f].get("breakpoint") is None]
    out = {
        "dataset": exp["dataset"],
        "species_slug": exp["species_slug"],
        "antibiotic": result.get("antibiotic", ""),
        "split": exp["split"],
        "plan": exp["plan"],
        "seed": exp["seed"],
        "draw": "draw" if int(exp["seed"]) in (42, 43, 44) else "twin",
        "arm": arm,
        "unavailable": why or "",
        "folds": len(folds),
        "folds_scored": len(scored),
        "folds_undefined": sum(1 for f in scored if not fold_has_fit(folds[f], arm)),
        "folds_without_breakpoint": len(no_bp),
        "ea_exact_n": len(exact),
        "category_n": len(rows),
        "very_major_n": len(resistant),
        "major_n": len(susceptible),
        "ea_exact_reads": EA_READS,
        "category_reads": CATEGORY_READS,
    }
    for key in ["ea_exact", "ea_exact_count", "category_agreement", "category_count",
                "very_major_rate", "very_major_count", "major_rate", "major_count"]:
        out[key] = ""
    if why:
        return out
    pred = [step(arm, r[arm]) for r in rows]
    if any(p is None for p in pred):
        raise ValueError(f"{exp['file']}: arm {arm} lacks a prediction")  # predictions_whole stops this first
    ea = sum(1 for r, p in zip(rows, pred) if r["censor"] == "Exact" and agrees(p, r))
    out["ea_exact_count"] = ea
    out["ea_exact"] = ea / len(exact) if exact else ""
    if no_bp or not rows:
        return out
    call = [p > folds[r["fold"]]["breakpoint"]["r_above"] for r, p in zip(rows, pred)]
    same = sum(1 for r, c in zip(rows, call) if c == r["resistant"])
    vme = sum(1 for r, c in zip(rows, call) if r["resistant"] and not c)
    me = sum(1 for r, c in zip(rows, call) if not r["resistant"] and c)
    out["category_count"], out["category_agreement"] = same, same / len(rows)
    out["very_major_count"] = vme
    out["very_major_rate"] = vme / len(resistant) if resistant else ""
    out["major_count"] = me
    out["major_rate"] = me / len(susceptible) if susceptible else ""
    return out


def check_one(exp, results_dir, release_dir):
    """(the check row, the result or None)."""
    row = {c: "" for c in CHECK_COLUMNS}
    row.update({c: exp[c] for c in ["dataset", "split", "plan", "seed", "file"]})
    path = os.path.join(results_dir, exp["file"])
    row["exists"] = "yes" if os.path.isfile(path) else "no"
    if row["exists"] == "no":
        row["ok"] = "no"
        return row, None
    with open(path) as f:
        result = json.load(f)
    row["format"] = "ok" if (
        result.get("format") == FORMAT
        and result.get("evaluation_revision") == REVISION
        and result.get("id") == exp["dataset"]
        and result.get("species_slug") == exp["species_slug"]
    ) else "differs"
    report_path = os.path.join(release_dir, "datasets", exp["dataset"] + ".json")
    if os.path.isfile(report_path):
        with open(report_path) as f:
            want = json.load(f).get("input_fingerprint")
        row["fingerprint"] = "equal" if want and result.get("input_fingerprint") == want else "differs"
    else:
        row["fingerprint"] = "no report"
    row["plan_and_arms"] = "ok" if plan_matches(result, exp) else "differs"
    row["fold_counts"] = "equal" if fold_counts_equal(result) else "differs"
    row["median_check"] = median_check(result)
    row["predictions"] = "whole" if predictions_whole(result) else "partial"
    row["ok"] = "yes" if (
        row["format"] == "ok"
        and row["fingerprint"] == "equal"
        and row["plan_and_arms"] == "ok"
        and row["fold_counts"] == "equal"
        and row["median_check"] != "differs"
        and row["predictions"] == "whole"
    ) else "no"
    return row, result


def write_tsv(path, columns, rows):
    def cell(v):
        return f"{v:.6f}" if isinstance(v, float) else str(v)

    with open(path, "w") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            f.write("\t".join(cell(r[c]) for c in columns) + "\n")


def main(argv):
    if len(argv) != 4:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 3
    expect_path, results_dir, release_dir, out = argv
    os.makedirs(out, exist_ok=True)
    expect = read_expect(expect_path)
    checks, scored = [], []
    for exp in expect:
        row, result = check_one(exp, results_dir, release_dir)
        checks.append(row)
        if result is not None and row["ok"] == "yes":
            arms = [a for a in ARMS if a == "median" or a in exp["arms"].split(",")]
            scored.extend(score_arm(result, exp, arm) for arm in arms)
    write_tsv(os.path.join(out, "checks.tsv"), CHECK_COLUMNS, checks)
    missing = sum(1 for r in checks if r["exists"] == "no")
    failed = sum(1 for r in checks if r["ok"] == "no")
    print(f"checks: {len(checks)} results expected, {missing} missing, {failed} failed")
    if os.path.exists(os.path.join(out, "scores.tsv")):
        os.remove(os.path.join(out, "scores.tsv"))
    if missing:
        return 2
    if failed:
        return 1
    write_tsv(os.path.join(out, "scores.tsv"), SCORE_COLUMNS, scored)
    print(f"scores: {len(scored)} rows written")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
