"""Holds every C2 fit's folds to the design (docs/cap-arms.md, The held-out genomes).

  python3 -I check_folds.py lineage <fits dir> <dataset> <cohorts dir> <structure.json> <roots.json> <draws> <twin> <out tsv>
                                   <added-labelled.tsv>
  python3 -I check_folds.py study <plans dir> <cohorts dir> <genome meta .jsonl.gz[,...]> <structure.json> <arms> <out tsv>

Reads fold plans, cohort plans, `grm arm --trace-training` files and fold
files only: no prediction, no score and no label beyond what a plan holds.
Prints the count of each breach and of each checked fold, never a genome
id, and exits 1 on any breach.

lineage: the default, <fits dir>/cap-2000/<dataset>, has a plan
(plans/lineage-<slot>.json) and a training trace for every slot 0 to
draws+twin-1; each wider arm of roots.json's `fits` has a cohort plan
(plans/lineage-<slot>.cohort.json) for slots 0 to draws-1, made from the
default's plan of that slot, unless cohorts.py's added-labelled.tsv marks
the arm none-added on the dataset: then the arm has no fit folder there (a
fit of a none-added arm breaches). Per slot and outer fold f:

- the default tests every one of its genomes once and trains on the rest;
- a wider arm's scored genomes are the default's test genomes of fold f,
  the same set (an outer test fold that differs between arms breaches);
- no added genome (in the arm, not in the default) is scored;
- no added genome trains in f when its lineage in the species grouping is
  the lineage of a genome scored in f, or its k-mer digest a scored
  genome's;
- the default's genomes train in the arm exactly as in the default;
- the cohort plan's lineage_dropped counts the added genomes of a scored
  lineage, and its balance_dropped is 0.

Writes one row per dataset, arm, split, slot and fold: scored, trained,
trained_added, kept_out_lineage, kept_out_duplicate, balance_dropped.

study: the BioProjects come from one or more genome row files, a comma
list (the master genome rows and bioprojects.py's answers for the added
genomes); a genome two files give different BioProjects is unknown. Per
dataset of <plans dir>/study.txt and wider arm with a fold file
<plans dir>/<arm>/<dataset>.study-whole.json, against the reviewed plan
<plans dir>/reviewed/folds/<dataset>.study-whole.json:

- the arm's test folds are the default's, fold by fold;
- no added genome is in a test fold;
- every added genome in train has a known BioProject that no tested genome
  has, and a k-mer digest no tested genome has;
- the default's train and unused genomes stay where they were.

Writes one row per dataset and arm: folds, tested, train_default,
train_added, unused_added.
"""

import gzip
import json
import sys
from pathlib import Path

DEFAULT = "cap-2000"
COLUMNS = ["dataset", "arm", "split", "slot", "fold", "scored", "trained", "trained_added", "kept_out_lineage",
           "kept_out_duplicate", "balance_dropped"]
STUDY_COLUMNS = ["dataset", "arm", "folds", "tested", "train_default", "train_added", "unused_added"]


def load_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def read_ids(path):
    return [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]


def grouping(structure):
    """{genome: lineage}, {genome: kset digest} of the species grouping."""
    s = load_json(structure)
    return dict(zip(s["genomes"], s["lineages"]["labels"])), dict(zip(s["genomes"], s["kset_digests"]))


def trace(path):
    """(sample ids, {(split, slot, fold): set of training ids}). grm arm writes a trace line's `repeat` as the
    slot, while a plan's `repeat` is the keyed repeat, which differs from the slot when a refused draw was
    replaced (src/ablate/arm.rs, arm_plans), so traces are read by slot."""
    with open(path, encoding="utf-8") as f:
        ids = json.loads(f.readline())["sample_ids"]
        out = {}
        for line in f:
            if line.strip():
                r = json.loads(line)
                out[(r["split"], r["repeat"], r["fold"])] = set(r["train"])
    return ids, out


class Breaches:
    def __init__(self):
        self.counts = {}

    def add(self, what, n=1):
        if n:
            self.counts[what] = self.counts.get(what, 0) + n

    def report(self):
        for what, n in sorted(self.counts.items()):
            print(f"breach: {what}: {n}")
        return 1 if self.counts else 0


def lineage(fits, dataset, cohorts, structure, roots, draws, twin, bad, empty=frozenset()):
    lin, digest = grouping(structure)
    base_dir = Path(fits) / DEFAULT / dataset
    base_ids, base_train = trace(base_dir / "training.jsonl")
    if not set(base_ids) <= set(read_ids(Path(cohorts) / f"{DEFAULT}.txt")):
        bad.add("a default genome outside the default's cohort")
    base_set = set(base_ids)
    rows, tests = [], {}
    for slot in range(draws + twin):
        p = base_dir / "plans" / f"lineage-{slot}.json"
        if not p.is_file():
            bad.add("a default slot with no plan")
            continue
        doc = load_json(p)
        plan = doc["plan"]
        if doc["sample_ids"] != base_ids:
            bad.add("a default plan over other genomes than its trace")
            continue
        folds = {}
        for i, f in enumerate(plan["outer"]):
            folds.setdefault(f, set()).add(base_ids[i])
        if sorted(folds) != list(range(plan["k"])):
            bad.add("a default plan whose outer folds are not 0 to k-1")
        tests[slot] = (plan, folds)
        for f, test in sorted(folds.items()):
            train = base_train.get(("lineage", slot, f))
            if train is None:
                bad.add("a default fold with no training trace")
                continue
            bad.add("a default fold that trains on other genomes than the rest", int(train != base_set - test))
            rows.append([dataset, DEFAULT, "lineage", slot, f, len(test), len(train), 0, 0, 0, 0])
    for arm, _root in roots["fits"].items():
        arm_dir = Path(fits) / arm / dataset
        if (arm, dataset) in empty:
            bad.add(f"{arm}: a none-added arm fitted", int(arm_dir.exists()))
            continue
        ids, arm_train = trace(arm_dir / "training.jsonl")
        if not set(ids) <= set(read_ids(Path(cohorts) / f"{arm}.txt")) and arm in roots["genomes"]:
            bad.add(f"{arm}: a genome outside its cohort")
        added = set(ids) - base_set
        for slot in range(draws):
            p = arm_dir / "plans" / f"lineage-{slot}.cohort.json"
            if not p.is_file() or slot not in tests:
                bad.add(f"{arm}: a slot with no cohort plan")
                continue
            c = load_json(p)
            plan, folds = tests[slot]
            if c["reference_sha256"] != plan["sha256"] or c["repeat"] != plan["repeat"]:
                bad.add(f"{arm}: a cohort plan made from another plan than the default's slot")
            scored = {}
            for i, f in enumerate(c["outer"]):
                if f is not None:
                    scored.setdefault(f, set()).add(ids[i])
            for f in range(plan["k"]):
                s = scored.get(f, set())
                test = folds.get(f, set())
                bad.add("an outer test fold that differs between arms", int(s != test))
                bad.add("an added genome in an outer test fold", len(s & added))
                train = arm_train.get(("lineage", slot, f))
                if train is None:
                    bad.add(f"{arm}: a fold with no training trace")
                    continue
                bad.add("a scored genome in training", len(train & s))
                lineages = {lin[g] for g in s}
                digests = {digest[g] for g in s}
                t_added = train & added
                bad.add("an added genome of a held-out lineage in training", sum(lin[g] in lineages for g in t_added))
                bad.add("an added copy of a scored genome in training", sum(digest[g] in digests for g in t_added))
                base_f = base_train.get(("lineage", slot, f))
                if base_f is not None:
                    bad.add("a default genome that trains otherwise than in the default", int((train & base_set) != base_f))
                kept_lin = sum(lin[g] in lineages for g in added)
                kept_dup = sum(digest[g] in digests and lin[g] not in lineages for g in added)
                bad.add("a cohort plan whose lineage count differs", int(c["lineage_dropped"][f] != kept_lin))
                bad.add("a balance drop", c["balance_dropped"][f])
                rows.append([dataset, arm, "lineage", slot, f, len(s), len(train), len(t_added), kept_lin, kept_dup,
                             c["balance_dropped"][f]])
    return rows


def bioprojects(meta):
    """{genome: BioProject} from a comma list of genome row files; a genome two files disagree on is left out."""
    out, differ = {}, set()
    for path in [p for p in str(meta).split(",") if p]:
        with gzip.open(path, "rt", encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    r = json.loads(line)
                    g, p = r["genome_id"], r.get("bioproject_accession")
                    if not p:
                        continue
                    if out.get(g, p) != p:
                        differ.add(g)
                    out[g] = p
    for g in differ:
        del out[g]
    return out


def none_added(path):
    """{(arm, dataset)} with no added genome labelled (cohorts.py's added-labelled.tsv)."""
    with open(path, encoding="utf-8") as f:
        head = f.readline().rstrip("\n").split("\t")
        rows = [dict(zip(head, line.rstrip("\n").split("\t"))) for line in f if line.strip()]
    return {(r["arm"], r["dataset"]) for r in rows if r["none_added"] == "1"}


def study(plans, cohorts, meta, structure, arms, bad):
    bp = bioprojects(meta)
    _, digest = grouping(structure)
    base = set(read_ids(Path(cohorts) / f"{DEFAULT}.txt"))
    rows = []
    for d in read_ids(Path(plans) / "study.txt"):
        f0 = load_json(Path(plans) / "reviewed" / "folds" / f"{d}.study-whole.json")
        tested = set().union(*map(set, f0["test"])) if f0["test"] else set()
        test_bp = {bp[g] for g in tested if g in bp}
        test_digest = {digest[g] for g in tested if g in digest}
        for arm in arms:
            p = Path(plans) / arm / f"{d}.study-whole.json"
            if not p.is_file():
                continue
            fa = load_json(p)
            same = len(fa["test"]) == len(f0["test"]) and all(set(a) == set(b) for a, b in zip(fa["test"], f0["test"]))
            bad.add("an outer test fold that differs between arms", int(not same))
            arm_tested = set().union(*map(set, fa["test"])) if fa["test"] else set()
            bad.add("an added genome in an outer test fold", len(arm_tested - base))
            train_added = set(fa["train"]) - base
            bad.add("an added genome of no known BioProject in training", sum(g not in bp for g in train_added))
            bad.add("an added genome of a held-out BioProject in training",
                    sum(g in bp and bp[g] in test_bp for g in train_added))
            bad.add("an added copy of a tested genome in training",
                    sum(digest.get(g) in test_digest for g in train_added))
            bad.add("a default genome moved", int(not (set(f0["train"]) <= set(fa["train"])
                                                     and set(f0["unused"]) <= set(fa["unused"]))))
            unused_added = set(fa["unused"]) - base
            rows.append([d, arm, len(fa["test"]), len(arm_tested), len(set(fa["train"]) & base), len(train_added),
                         len(unused_added)])
    return rows


def write(path, columns, rows):
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")


def main(argv):
    bad = Breaches()
    if len(argv) == 11 and argv[1] == "lineage":
        fits, dataset, cohorts, structure, roots_json, draws, twin, out, labelled = argv[2:]
        rows = lineage(fits, dataset, cohorts, structure, load_json(roots_json), int(draws), int(twin), bad,
                       none_added(labelled))
        write(out, COLUMNS, rows)
        print(f"{dataset}: lineage folds checked {len(rows)}")
    elif len(argv) == 8 and argv[1] == "study":
        plans, cohorts, meta, structure, arms, out = argv[2:]
        rows = study(plans, cohorts, meta, structure, [a for a in arms.split() if a], bad)
        write(out, STUDY_COLUMNS, rows)
        print(f"study fold files checked {len(rows)}")
    else:
        raise SystemExit(__doc__)
    return bad.report()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
