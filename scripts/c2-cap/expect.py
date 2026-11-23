"""The rows one species' C2 run must score (docs/cap-arms.md), as the list score.py reads.

  python3 -I expect.py <species_slug> <datasets> <plans dir> <draws> <twin> <arms> <splits> <roots.json>
                       <added-labelled.tsv>

<datasets> is a comma list; <draws> and <twin> are counts of grm arm's
keyed lineage slots; <arms> and <splits> comma lists; <plans dir> holds
study.txt (the datasets with a reviewed study plan) and each wider arm's
study fold files, <arm>/<dataset>.study-whole.json; <roots.json> is
roots.py's plan (`fits` and `same_as`); <added-labelled.tsv> is
cohorts.py's count of each wider arm's added genomes with a call per
dataset.

Per dataset, on the lineage split: the default (cap-2000) at slots 0 to
draws-1 (role draw), the twin, the default's slots draws to
draws+twin-1 (role twin, fitted by cap-2000), and each wider arm at slots
0 to draws-1 (role draw). On the study split, per study dataset, one row
for the default and one per wider arm (role study, slot -); a wider arm
with no fold file adds no training genome there (role none-added). A wider
arm with no added genome labelled on the dataset is none-added on both
splits. A none-added row is not fitted and stays out of the count's
denominator. An arm that `same_as` names shares another arm's cohort and
fit (Enterococcus faecium's cap-4000 and cap-all): its rows say so in
`fit`, take role same-cohort, and are counted under the other arm only.

Columns: dataset, species_slug, arm, split, slot, role, fit (the arm whose
grm run holds the row).
"""

import json
import sys
from pathlib import Path

DEFAULT, TWIN = "cap-2000", "twin"
WIDER = ("cap-4000", "cap-all")
STUDY_SCHEME = "study-whole"
COLUMNS = ["dataset", "species_slug", "arm", "split", "slot", "role", "fit"]


def comma(text):
    return [] if text in ("", "-") else text.split(",")


def none_added(path):
    """{(arm, dataset)} with no added genome labelled, from cohorts.py's added-labelled.tsv."""
    with open(path, encoding="utf-8") as f:
        head = f.readline().rstrip("\n").split("\t")
        rows = [dict(zip(head, line.rstrip("\n").split("\t"))) for line in f if line.strip()]
    return {(r["arm"], r["dataset"]) for r in rows if r["none_added"] == "1"}


def rows(slug, datasets, study, draws, twin, arms, splits, roots, have_file, empty=frozenset()):
    if draws < 1 or twin < 0:
        raise ValueError("draws is 1 or more and twin 0 or more")
    if any(a not in WIDER for a in arms) or len(set(arms)) != len(arms):
        raise ValueError(f"the wider arms are {', '.join(WIDER)}")
    if not set(study) <= set(datasets):
        raise ValueError("a study dataset is not listed")
    same = roots.get("same_as", {})
    fit_of = {a: same.get(a, a) for a in arms}
    if any(fit_of[a] not in roots["fits"] for a in arms):
        raise ValueError("an arm has no fit in roots.json")
    out = []
    for d in datasets:
        if "lineage" in splits:
            out += [[d, slug, DEFAULT, "lineage", str(s), "draw", DEFAULT] for s in range(draws)]
            out += [[d, slug, TWIN, "lineage", str(s), "twin", DEFAULT] for s in range(draws, draws + twin)]
            for a in arms:
                role = "none-added" if (a, d) in empty else "same-cohort" if a in same else "draw"
                out += [[d, slug, a, "lineage", str(s), role, fit_of[a]] for s in range(draws)]
        if "study" in splits and d in study:
            out.append([d, slug, DEFAULT, "study", "-", "study", DEFAULT])
            for a in arms:
                role = ("none-added" if (a, d) in empty or not have_file(fit_of[a], d)
                        else "same-cohort" if a in same else "study")
                out.append([d, slug, a, "study", "-", role, fit_of[a]])
    return out


def main(argv):
    if len(argv) != 10:
        raise SystemExit(__doc__)
    slug, datasets, plans, draws, twin, arms, splits, roots_json, labelled = argv[1:]
    plans = Path(plans)
    study = [line.strip() for line in (plans / "study.txt").read_text(encoding="utf-8").splitlines() if line.strip()]
    with open(roots_json, encoding="utf-8") as f:
        roots = json.load(f)
    try:
        out = rows(slug, comma(datasets), study, int(draws), int(twin), comma(arms), comma(splits), roots,
                   lambda a, d: (plans / a / f"{d}.{STUDY_SCHEME}.json").is_file(), none_added(labelled))
    except ValueError as e:
        raise SystemExit(str(e)) from None
    print("\t".join(COLUMNS))
    for r in out:
        print("\t".join(r))


if __name__ == "__main__":
    main(sys.argv)
