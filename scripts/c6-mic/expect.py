"""The results one species' MIC run must hold (docs/mic-plan.md), as the list score.py reads.

  python3 -I expect.py <species_slug> <lineage datasets> <study datasets> > expect.tsv

Each dataset list is a comma list of dataset ids, or - for none. Every
dataset id must begin with <species_slug>__.

Lineage split: seeds 42 to 47, each its own plan, lineage-s<seed>. Study
split: six fold files, one per seed: the five reviewed study plans
(study-whole at 42, study-d1 to study-d4 at 43 to 46) and c6mic-s47 at 47.
On both, seeds 42, 43 and 44 fit the released and the ordered model, and
seeds 45, 46 and 47 (the twin) the released model only. Every result also
holds the training median.

Columns: dataset, species_slug, split, plan, seed, arms, file (the result's
path under the results directory: lineage/<dataset>.<plan>.mic.json or
study/<dataset>.<plan>.mic.json, the names grm learn writes).
"""

import sys

SEEDS = [42, 43, 44, 45, 46, 47]
STUDY_PLANS = {42: "study-whole", 43: "study-d1", 44: "study-d2", 45: "study-d3", 46: "study-d4", 47: "c6mic-s47"}
COLUMNS = ["dataset", "species_slug", "split", "plan", "seed", "arms", "file"]


def arms(seed):
    return "released,ordered" if seed in (42, 43, 44) else "released"


def datasets(text, species):
    if text == "-":
        return []
    out = [d for d in text.split(",") if d]
    for d in out:
        if not d.startswith(species + "__"):
            raise ValueError(f"{d} is not a dataset of {species}")
    if len(set(out)) != len(out):
        raise ValueError("a dataset is listed twice")
    return out


def rows(species, lineage, study):
    out = []
    for split, names in (("lineage", lineage), ("study", study)):
        for d in names:
            for seed in SEEDS:
                plan = f"lineage-s{seed}" if split == "lineage" else STUDY_PLANS[seed]
                out.append({
                    "dataset": d, "species_slug": species, "split": split, "plan": plan,
                    "seed": seed, "arms": arms(seed), "file": f"{split}/{d}.{plan}.mic.json",
                })
    return out


def main(argv):
    if len(argv) != 3:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    species = argv[0]
    table = rows(species, datasets(argv[1], species), datasets(argv[2], species))
    if not table:
        print("no dataset given", file=sys.stderr)
        return 2
    sys.stdout.write("\t".join(COLUMNS) + "\n")
    for r in table:
        sys.stdout.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
