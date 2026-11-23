"""Shared loaders and SYNTHETIC fixtures for the cap-arm scripts (scripts/c2-cap, docs/cap-arms.md).

Every genome, lineage, BioProject, plan, trace, prediction and pick file here
is made up in the test, in the shapes grm arm and grm learn write. No genome
or call comes from a real dataset, and the species and dataset names are made
up, except that the scorer's first-run list is the committed
scripts/c2-cap/lists.tsv, so the scorer's fixtures borrow one dataset name
from it.
"""

import contextlib
import importlib.util
import io
import json
import pathlib
import random

ROOT = pathlib.Path(__file__).resolve().parents[2]
KIT = ROOT / "scripts" / "c2-cap"
SPECIES = "synthia_coli"
K = 5
LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]


def load(name):
    spec = importlib.util.spec_from_file_location(f"c2_cap_{name}", KIT / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def quiet(fn, *args):
    """(return value, printed text) of fn(*args)."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = fn(*args)
    return out, buf.getvalue()


def write_json(path, doc):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc), encoding="utf-8")


def write_ids(path, ids):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(g + "\n" for g in ids), encoding="utf-8")


class Species:
    """A SYNTHETIC species: 60 default genomes in 12 lineages and 30 added genomes, 6 of them in lineages
    of their own; two added genomes are k-mer copies of default genomes."""

    def __init__(self, seed=7):
        rng = random.Random(seed)
        self.base = [f"9{i:03d}.1" for i in range(60)]
        self.added = [f"9{i:03d}.2" for i in range(30)]
        self.lineage = {g: i % 12 for i, g in enumerate(self.base)}
        for i, g in enumerate(self.added):
            self.lineage[g] = 12 + i % 6 if i < 6 else i % 12
        self.digest = {g: f"d{g}" for g in self.base + self.added}
        self.digest[self.added[10]] = self.digest[self.base[3]]
        self.digest[self.added[11]] = self.digest[self.base[4]]
        self.y = {g: rng.random() < 0.4 for g in self.base + self.added}
        self.y[self.base[0]], self.y[self.base[1]] = True, False

    def structure(self, path):
        genomes = sorted(self.base + self.added)
        write_json(path, {"genomes": genomes, "k": 31, "kset_digests": [self.digest[g] for g in genomes],
                          "lineages": {"labels": [self.lineage[g] for g in genomes]}})

    def cohorts(self, folder):
        write_ids(folder / "cap-2000.txt", sorted(self.base))
        write_ids(folder / "cap-4000.txt", sorted(self.base + self.added))
        write_ids(folder / "cap-all.txt", sorted(self.base + self.added))


def lineage_fits(sp, fits, dataset, slots, draws, arms=("cap-all",), breach=None, keyed=None, refused=()):
    """Writes the default's and each arm's plans and training traces as grm arm does for a correct fit, then
    plants `breach`: 'lineage' (an added genome of a scored lineage trains), 'differs' (a default test genome
    left out of an arm's scored set), 'added-scored' (an added genome scored beside the default's test genomes),
    'copy' (an added copy of a scored genome trains). `keyed` maps a slot to the keyed repeat that replaced its
    refused draw: the plans and cohort plans hold that repeat, the trace lines the slot, as grm arm writes them.
    `refused` lists slots whose default plan grm arm refused past its replacements: those plans hold an
    unavailable_reason, the default traces no lineage fold, and each arm has no cohort plan at a refused slot and
    traces only the slots before the first one, as grm arm does when the default's lineage split is unavailable."""
    keyed = keyed or {}
    cut = min(refused) if refused else None
    base = sorted(sp.base)
    lineages = sorted({sp.lineage[g] for g in base})
    for slot in range(slots):
        rng = random.Random(100 + slot)
        order = list(lineages)
        rng.shuffle(order)
        fold_of_lineage = {c: i % K for i, c in enumerate(order)}
        outer = [fold_of_lineage[sp.lineage[g]] for g in base]
        plan = {"repeat": keyed.get(slot, slot), "k": K, "outer": outer, "sha256": f"plan{slot}"}
        if slot in refused:
            plan["unavailable_reason"] = "SYNTHETIC: refused past the replacements"
        write_json(fits / "cap-2000" / dataset / "plans" / f"lineage-{slot}.json",
                   {"sample_ids": base, "y": [int(sp.y[g]) for g in base], "plan": plan})
    with_trace(fits / "cap-2000" / dataset, base,
               {} if refused else {(slot, f): [g for g, o in zip(base, plan_outer(fits, dataset, slot)) if o != f]
                                   for slot in range(slots) for f in range(K)})
    for arm in arms:
        ids = sorted(sp.base + sp.added)
        train = {}
        for slot in range(draws):
            if slot in refused:
                continue
            outer_base = dict(zip(base, plan_outer(fits, dataset, slot)))
            outer = [outer_base.get(g) for g in ids]
            lineage_dropped, balance = [], []
            for f in range(K):
                scored = {g for g in base if outer_base[g] == f}
                lin = {sp.lineage[g] for g in scored}
                dig = {sp.digest[g] for g in scored}
                kept = [g for g in sp.added if sp.lineage[g] in lin]
                lineage_dropped.append(len(kept))
                balance.append(0)
                train[(slot, f)] = [g for g in ids if outer_base.get(g, -1) != f and g not in kept
                                    and not (g in sp.added and sp.digest[g] in dig)]
            if breach == "lineage" and slot == 0:
                f = 0
                scored_lin = {sp.lineage[g] for g in base if outer_base[g] == f}
                leak = next(g for g in sp.added if sp.lineage[g] in scored_lin)
                train[(slot, f)].append(leak)
            if breach == "copy" and slot == 0:
                copy = sp.added[10]
                f = outer_base[sp.base[3]]
                if copy not in train[(slot, f)]:
                    train[(slot, f)].append(copy)
            if breach == "differs" and slot == 0:
                i = ids.index(next(g for g in base if outer_base[g] == 0))
                outer[i] = None
            if breach == "added-scored" and slot == 0:
                i = ids.index(sp.added[0])
                outer[i] = 0
            write_json(fits / arm / dataset / "plans" / f"lineage-{slot}.cohort.json",
                       {"reference_sha256": f"plan{slot}", "repeat": keyed.get(slot, slot), "k": K, "outer": outer,
                        "lineage_dropped": lineage_dropped, "balance_dropped": balance})
        with_trace(fits / arm / dataset, ids, {k: v for k, v in train.items() if cut is None or k[0] < cut})


def plan_outer(fits, dataset, slot):
    p = fits / "cap-2000" / dataset / "plans" / f"lineage-{slot}.json"
    return json.loads(p.read_text(encoding="utf-8"))["plan"]["outer"]


def with_trace(folder, ids, train):
    folder.mkdir(parents=True, exist_ok=True)
    lines = [json.dumps({"sample_ids": ids})]
    for (slot, f), t in sorted(train.items()):
        lines.append(json.dumps({"split": "lineage", "repeat": slot, "fold": f, "train": sorted(t)}))
    (folder / "training.jsonl").write_text("\n".join(lines) + "\n", encoding="utf-8")


def predictions(fits, fit, dataset, slots, arm_name=None, drop_slot=None):
    """predictions.jsonl for a fit from its plans: every scored genome, every learner, with the slot as `repeat`
    (grm arm writes the slot there, not the plan's keyed repeat)."""
    rows = []
    folder = fits / fit / dataset
    for slot in slots:
        if slot == drop_slot:
            continue
        if fit == "cap-2000":
            doc = json.loads((folder / "plans" / f"lineage-{slot}.json").read_text(encoding="utf-8"))
            ids, outer = doc["sample_ids"], doc["plan"]["outer"]
            y = doc["y"]
        else:
            doc = json.loads((folder / "plans" / f"lineage-{slot}.cohort.json").read_text(encoding="utf-8"))
            ids = json.loads((folder / "training.jsonl").read_text(encoding="utf-8").splitlines()[0])["sample_ids"]
            outer = doc["outer"]
            y = [int(g.endswith("1") or g.endswith("3")) for g in ids]
        for learner in LEARNERS:
            for i, g in enumerate(ids):
                if outer[i] is None:
                    continue
                score = ((i * 7 + slot) % 10) / 10
                rows.append(json.dumps({"dataset": dataset, "arm": arm_name or fit, "learner": learner,
                                        "split": "lineage", "repeat": slot, "sample_id": g, "y": int(y[i]),
                                        "fold": outer[i], "score": score, "called": score >= 0.5,
                                        "cluster": 0, "fold_plan_sha256": f"plan{slot}"}))
    (folder / "predictions.jsonl").write_text("\n".join(rows) + "\n", encoding="utf-8")


def fold_file(sp, dataset):
    """The default's SYNTHETIC study-whole fold file: the 60 default genomes, BioProjects P0 to P4 one per test
    fold, and two genomes of no BioProject in unused."""
    tests = [[g for i, g in enumerate(sorted(sp.base)) if i % K == f and i >= 2] for f in range(K)]
    return {"format": "grm-fold-file-v1", "dataset": dataset, "scheme": "study-whole", "test": tests, "train": [],
            "unused": sorted(sp.base)[:2]}


def bioprojects(sp, f0):
    """{genome: BioProject}: test fold f's genomes in P<f>; the added genomes in P0 (held out) for the first 8,
    in P9 (new) for the next 12, and the last 10 with none."""
    out = {}
    for f, t in enumerate(f0["test"]):
        for g in t:
            out[g] = f"PRJ{f}"
    for i, g in enumerate(sp.added):
        if i < 8:
            out[g] = "PRJ0"
        elif i < 20:
            out[g] = "PRJ9"
    return out


def meta(path, bp):
    import gzip
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as f:
        for g, p in sorted(bp.items()):
            f.write(json.dumps({"genome_id": g, "bioproject_accession": p}) + "\n")
