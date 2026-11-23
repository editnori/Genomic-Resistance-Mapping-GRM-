"""SYNTHETIC readings for scripts/c5-settings/page_tables.py and SYNTHETIC run copies for run_checks.py.

Three made-up species, "Synthia alpha" with two datasets and "Synthia beta"
and "Synthia gamma" with one each, and a made-up dataset of M beside them.
Every draw, choices file, count and genome id here is made up. Each figure
is a multiple of 1/64, and each arm's difference from the default is one
value over the datasets, so every count row's differences are one value and
the count needs no resampling.
"""

import hashlib
import json
import pathlib

from tests.c5_settings import common

score = common.load("score")
C3 = score.C3
LEARNERS = list(common.LEARNERS)
SPECIES = [("synthia_alpha", "Synthia alpha", ["synthetic drug 1", "synthetic drug 2"]),
           ("synthia_beta", "Synthia beta", ["synthetic drug 3"]),
           ("synthia_gamma", "Synthia gamma", ["synthetic drug 4"])]
M_DATASET = ("synthia_alpha", "Synthia alpha", "synthetic drug 9")
LIST_COLUMNS = ["dataset", "species_slug", "species", "antibiotic", "samples", "resistant", "susceptible", "lineages",
                "in_M", "first_run", "pairs_S+R", "labels_depend_on"]
SEEDS, TWIN, STUDY = (42, 43, 44), (45, 46, 47), 42
FOLDS = 5
BASE = {"balanced_accuracy": 48, "recall": 40, "specificity": 56, "vme": 24, "me": 8, "log_loss": 20}
# Each arm's difference from the default, in 1/64, on the measures where higher is better; lower-is-better
# measures move the other way, and the cutoff arm's log loss not at all.
DELTA = {"default": 0, "search": 4, "cutoff": 2, "twin": 0, "plain-best": 4}
SETTING = {"max_depth": 4321}
COMMIT = "c" * 40
GRM = "a" * 64
CODE = "b" * 64


def slug_of(name):
    return name.replace(" ", "_")


def dataset_name(slug, antibiotic):
    return f"{slug}__{slug_of(antibiotic)}"


def lists_rows(m_in_run=False):
    rows = []
    for slug, species, drugs in SPECIES:
        for i, ab in enumerate(drugs):
            rows.append({"dataset": dataset_name(slug, ab), "species_slug": slug, "species": species,
                         "antibiotic": ab, "samples": 120 + 10 * i, "resistant": 60, "susceptible": 60 + 10 * i,
                         "lineages": 20, "in_M": 0, "first_run": 1, "pairs_S+R": 0, "labels_depend_on": "none"})
    slug, species, ab = M_DATASET
    rows.append({"dataset": dataset_name(slug, ab), "species_slug": slug, "species": species, "antibiotic": ab,
                 "samples": 90, "resistant": 30, "susceptible": 60, "lineages": 15, "in_M": 1,
                 "first_run": int(m_in_run),
                 "pairs_S+R": 3, "labels_depend_on": "the intermediate call"})
    return rows


def expect_rows(slug, datasets):
    exp = common.load("expect")
    return [dict(zip(exp.COLUMNS, r)) for r in exp.rows(slug, datasets, datasets, list(SEEDS), list(TWIN))]


def values(arm, learner):
    """{measure: value} of one arm and learner; scm's cutoff arm is the default."""
    d = 0 if (learner == "scm" and arm == "cutoff") else DELTA[arm]
    out = {}
    for m, v in BASE.items():
        if m == "log_loss":
            step = 0 if arm == "cutoff" else -d // 2
        elif m in C3.LOWER_IS_BETTER:
            step = -d
        else:
            step = d
        out[m] = (v + step) / 64
    return out


def draw_row(e, learner, species, antibiotic, scored=1, constant=0, why=""):
    r = {"dataset": e["dataset"], "species": species, "antibiotic": antibiotic, "arm": e["arm"], "split": e["split"],
         "seed": int(e["seed"]), "role": e["role"], "learner": learner, "scored": scored, "constant": constant,
         "why": why}
    if scored:
        r.update(n=64, resistant=32, susceptible=32, tp=20, fp=4, tn=28, fn=12)
        r.update(values(e["arm"], learner))
    return r


def calls():
    return {"n": 13, "tp": 6, "fp": 1, "tn": 5, "fn": 1, "loss_sum": 5.0}


def choice_doc(dataset, scheme, seed, fold):
    """A choices file as the driver writes it, its five learners each in a different place: scm keeps trial 0,
    cart keeps trial 0 with a plain best of trial 3 that differs, rf picks trial 7, gbt's trial 0 failed inside and
    lasso's every trial failed."""
    c = calls()
    learners = [
        dict(learner="scm", status="picked", picked_trial=0, plain_best_trial=0, release_inner_mean=0.75,
             best_inner_mean=0.75, cutoff=0.5, plain_best_refitted=False, plain_best_note=None),
        dict(learner="cart", status="picked", picked_trial=0, plain_best_trial=3, release_inner_mean=0.75,
             best_inner_mean=0.78125, cutoff=0.45, plain_best_refitted=True, plain_best_note=None),
        dict(learner="rf", status="picked", picked_trial=7, plain_best_trial=7, release_inner_mean=0.75,
             best_inner_mean=0.8125, cutoff=0.55, plain_best_refitted=False, plain_best_note=None),
        dict(learner="gbt", status="release failed inside", picked_trial=0, plain_best_trial=5,
             release_inner_mean=None, best_inner_mean=0.6875, cutoff=0.5, plain_best_refitted=True,
             plain_best_note=None),
        dict(learner="lasso", status="release failed inside", picked_trial=0, plain_best_trial=None,
             release_inner_mean=None, best_inner_mean=None, cutoff=0.5, plain_best_refitted=False,
             plain_best_note="no trial scored inside; the release's is served"),
    ]
    for r in learners:
        r.update(trials=31, failed_trials=31 if r["plain_best_trial"] is None else 0, picked_setting=dict(SETTING),
                 tied_at_best=1, trials_above_gamma=0, trials_better=0, trials_same=0, trials_worse=0,
                 cutoff_inner_balanced_accuracy=0.75, half_inner_balanced_accuracy=0.75, search_note=None,
                 outer={"default": dict(c), "search": dict(c), "cutoff": dict(c), "plain-best": dict(c)},
                 seconds_inner=1.0, seconds_outer=1.0)
    return {"format": score.CHOICE_FORMAT, "evaluation_revision": 14, "dataset": dataset, "scheme": scheme,
            "seed": seed, "fold": fold, "tuning": "synthetic", "inner_seed": 1, "outer_seed": 2,
            "sampler_seed": score.SAMPLER_SEED, "gamma": score.GAMMA, "grid": list(score.GRID),
            "random_trials": score.TRIALS, "counts": {"train": 52, "train_resistant": 26, "test": 13,
                                                     "test_resistant": 7},
            "inner_folds": 3, "min_support_inner": 2, "min_support_outer": 2,
            "learners": learners, "timing": {}}


def write_tsv(path, rows, columns):
    C3.write_tsv(path, rows, columns)


def build(top, missing=False, constant=False, change_draws=None, change_choice=None, species=None, m_in_run=False):
    """top/lists.tsv, top/checks.txt and top/reading: a reading of every species that passes every gate.

    missing: Synthia alpha's first dataset has no choices files at lineage seed 43, and its three arms are a
    missing draw there. constant: Synthia beta's search study run is a constant model, refused for no candidate.
    m_in_run: lists.tsv marks the dataset of M as on the first run, and the reading holds it.
    change_draws(slug, draws) changes the draws before the tables are made from them; change_choice(key, doc)
    returns a choices file's doc, or None for no file.
    """
    top = pathlib.Path(top)
    lrows = lists_rows(m_in_run)
    write_tsv(top / "lists.tsv", lrows, LIST_COLUMNS)
    (top / "checks.txt").write_text("PASS the first run's species, each once: 3 folders for 3 species\n"
                                    "CHECKS 1 of 1\n", encoding="utf-8")
    reading = top / "reading"
    (reading / "count").mkdir(parents=True)
    (reading / "reading.txt").write_text("score.py --count: SYNTHETIC\n", encoding="utf-8")
    every = []
    for slug, name, drugs in SPECIES if species is None else [s for s in SPECIES if s[0] in species]:
        if m_in_run and slug == M_DATASET[0]:
            drugs = drugs + [M_DATASET[2]]
        folder = reading / "score" / slug
        (folder / "choices").mkdir(parents=True)
        datasets = [dataset_name(slug, ab) for ab in drugs]
        antibiotic = {dataset_name(slug, ab): ab for ab in drugs}
        draws, gone = [], set()
        for e in expect_rows(slug, datasets):
            ab = antibiotic[e["dataset"]]
            scheme = e["scheme"]
            lost = missing and e["dataset"] == datasets[0] and slug == "synthia_alpha" and e["split"] == "lineage" \
                and int(e["seed"]) == 43 and e["arm"] in score.ARMS
            refused = constant and slug == "synthia_beta" and e["split"] == "study" and e["arm"] == "search"
            for learner in LEARNERS:
                if lost:
                    draws.append(draw_row(e, learner, name, ab, scored=0, why="missing draw"))
                    gone.add((e["dataset"], scheme, int(e["seed"])))
                elif refused:
                    draws.append(draw_row(e, learner, name, ab, constant=1,
                                          why="constant model: refused for no candidate"))
                else:
                    draws.append(draw_row(e, learner, name, ab))
            if e["arm"] == "default" and e["role"] in ("draw", "study") \
                    and (e["dataset"], scheme, int(e["seed"])) not in gone and not lost:
                for k in range(FOLDS):
                    key = (e["dataset"], scheme, int(e["seed"]), k)
                    doc = choice_doc(*key)
                    if change_choice:
                        doc = change_choice(key, doc)
                    if doc is not None:
                        path = folder / "choices" / f"{e['dataset']}.{scheme}-s{e['seed']}.f{k}.choices.json"
                        path.write_text(json.dumps(doc), encoding="utf-8")
        if missing and slug == "synthia_alpha":
            for k in range(FOLDS):
                p = folder / "choices" / f"{datasets[0]}.lineage-s43.f{k}.choices.json"
                if p.exists():
                    p.unlink()
        if change_draws:
            change_draws(slug, draws)
        write_tsv(folder / "draws.tsv", draws, C3.DRAW_COLUMNS)
        rows = C3.dataset_rows(draws)
        write_tsv(folder / "datasets.tsv", rows, C3.DATASET_COLUMNS)
        write_tsv(folder / "checks.tsv", [{"file": "default/x.pick.json", "check": "fold 0 scm: SYNTHETIC",
                                           "passed": True, "detail": ""},
                                          {"file": "default/x.pick.json", "check": "fold 1 scm: SYNTHETIC",
                                           "passed": True, "detail": ""}], ["file", "check", "passed", "detail"])
        (folder / "binary.txt").write_text(f"commit {COMMIT}\ngrm {GRM}\ncode {CODE}\n", encoding="utf-8")
        every += C3.read_tsv(folder / "datasets.tsv")
    count_files(reading, every)
    return reading, top / "lists.tsv", top / "checks.txt"


def count_files(reading, rows):
    """The four count tables as score.py --count writes them from datasets.tsv rows."""
    counted = C3.count(rows)
    write_tsv(reading / "count" / "count.tsv", counted, C3.COUNT_COLUMNS)
    write_tsv(reading / "count" / "count-no-constant.tsv", C3.count(rows, keep_constant=False), C3.COUNT_COLUMNS)
    write_tsv(reading / "count" / "conditions.tsv", C3.conditions(counted), C3.CONDITION_COLUMNS)
    write_tsv(reading / "count" / "count-plain-best.tsv", score.plain_count(rows), C3.COUNT_COLUMNS)


def recount(reading):
    """The count tables again from the reading's datasets.tsv files, after a test changed them."""
    rows = [r for p in sorted((reading / "score").glob("*/datasets.tsv")) for r in C3.read_tsv(p)]
    count_files(reading, rows)


# ---------- run copies for run_checks.py ----------

def genome(slug, i):
    """A made-up id in the form of a BV-BRC genome id."""
    return f"{1000 + 100 * [s[0] for s in SPECIES].index(slug) + i}.{7 + i}"


def sha256(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def copies(top, change=None):
    """top/copies/<slug>/ as run_species.sh leaves its OUT, top/release/species/<slug>.json and top/lists.tsv.

    change(slug, out) changes a species' copy before its SHA256SUMS and DONE are written.
    """
    top = pathlib.Path(top)
    write_tsv(top / "lists.tsv", lists_rows(), LIST_COLUMNS)
    exp = common.load("expect")
    for slug, name, drugs in SPECIES:
        out = top / "copies" / slug
        picks = out / "results" / "picks"
        datasets = [dataset_name(slug, ab) for ab in drugs]
        want = expect_rows(slug, datasets)
        for e in want:
            if e["arm"] in ("default", "twin"):
                p = picks / e["file"]
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text(json.dumps({"format": "grm-c9-pick-1", "folds": [{"fold": k} for k in range(FOLDS)]}),
                             encoding="utf-8")
            else:
                for k in range(FOLDS):
                    p = picks / e["file"].replace("<k>", str(k))
                    p.parent.mkdir(parents=True, exist_ok=True)
                    if not p.exists():
                        p.write_text(json.dumps(choice_doc(e["dataset"], e["scheme"], int(e["seed"]), k)),
                                     encoding="utf-8")
        (out / "check").mkdir(parents=True)
        (out / "check" / "refit.tsv").write_text("check\tpassed\nSYNTHETIC\t1\n", encoding="utf-8")
        (out / "score").mkdir()
        write_tsv(out / "score" / "checks.tsv", [{"file": "x", "check": "SYNTHETIC", "passed": True, "detail": ""}],
                  ["file", "check", "passed", "detail"])
        write_tsv(out / "score" / "draws.tsv", [draw_row(dict(e, seed=e["seed"]), l, name, drugs[0])
                                                for e in want[:2] for l in LEARNERS], C3.DRAW_COLUMNS)
        write_tsv(out / "score" / "datasets.tsv", [], C3.DATASET_COLUMNS)
        (out / "score" / "score.log").write_text(f"SYNTHETIC log naming {genome(slug, 0)}\n", encoding="utf-8")
        (out / "logs").mkdir()
        (out / "logs" / "timing.tsv").write_text("step\tseconds\nfit\t1\n", encoding="utf-8")
        (out / "logs" / "drive-timing.tsv").write_text("drive\tseconds\n1\t1\n", encoding="utf-8")
        (out / "binary.txt").write_text(f"commit {COMMIT}\ngrm {GRM}\ncode {CODE}\n", encoding="utf-8")
        (out / "lists.txt").write_text(f"species\t{slug}\ndatasets\t{','.join(datasets)}\nseeds\t42,43,44\n"
                                       "twin\t45,46,47\n", encoding="utf-8")
        with open(out / "expect.tsv", "w", encoding="utf-8") as handle:
            handle.write("\t".join(exp.COLUMNS) + "\n")
            for r in exp.rows(slug, datasets, datasets, list(SEEDS), list(TWIN)):
                handle.write("\t".join(str(x) for x in r) + "\n")
        if change:
            change(slug, out)
        seal(out, slug)
        species_json = top / "release" / "species" / f"{slug}.json"
        species_json.parent.mkdir(parents=True, exist_ok=True)
        species_json.write_text(json.dumps({"genomes": [{"id": genome(slug, i)} for i in range(4)]}),
                                encoding="utf-8")
    return top / "copies", top / "lists.tsv", top / "release"


def seal(out, slug, time="2026-11-23T12:00:00Z"):
    """SHA256SUMS and DONE as run_species.sh writes them."""
    names = sorted({p.relative_to(out).as_posix() for top in ("results", "check", "score") if (out / top).is_dir()
                    for p in (out / top).rglob("*") if p.is_file() and p.suffix in (".json", ".tsv", ".txt")}
                   | {"binary.txt", "expect.tsv", "lists.txt", "logs/timing.tsv", "logs/drive-timing.tsv"})
    (out / "SHA256SUMS").write_text("".join(f"{sha256(out / n)}  {n}\n" for n in names if (out / n).is_file()),
                                    encoding="utf-8")
    picks = len(list((out / "results" / "picks").rglob("*.pick.json")))
    n_choices = len(list((out / "results" / "picks" / "choices").rglob("*.choices.json")))
    named = dict(line.split(" ", 1) for line in (out / "binary.txt").read_text(encoding="utf-8").splitlines())
    (out / "DONE").write_text(f"DONE {slug} {time} binary {named['grm']} code {named['code']} results {picks} "
                              f"choices {n_choices} "
                              f"SHA256SUMS {sha256(out / 'SHA256SUMS')}\n", encoding="utf-8")
