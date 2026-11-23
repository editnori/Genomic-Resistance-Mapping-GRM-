"""The tables of docs/served-pick-result.md: plan block C9's one reading, as files and as the page prints them.

  python3 -I page_tables.py write <reading dir> <l.tsv> <out dir>
  python3 -I page_tables.py print <out dir>

The reading directory holds reading.txt (the reading's commands, each with
the sha256 of what it wrote) and one folder per answer of the design's
questions 3 to 5 (docs/served-pick.md, The scorer):

  read/            --tie later --no-pick default --default per-draw: the answers read
  tie-earlier/     --tie earlier: scored beside, decides nothing
  tie-default/     --tie default: scored beside, decides nothing
  default-seed42/  --default seed42: scored beside, decides nothing

Each answer's folder holds count/ (score.py --count's count.tsv and
conditions.tsv) and score/<species slug>/ (score.py's datasets.tsv,
draws.tsv, picks.tsv and options.txt). <l.tsv> is scripts/c9-pick/l.tsv,
the list of L: the first run is its rows with first_run 1, which are the
rows not in M.

`write` passes every gate of GATES before it writes a file, and stops at
the first that fails, naming it. It then writes <out dir> (absent or
empty): the scorer's count files under count/<answer>/, reading.txt, the
page's tables as page-*.tsv, gates.txt, README.md and SHA256SUMS. Every row
is an aggregate over a dataset, a species, a split or a learner: no file
holds a genome, a sequence or a feature.

`print` prints the page's tables from <out dir> alone, figures rounded to
3 decimals, and tests/c9_pick checks that each line stands in the page.
"""

import csv
import hashlib
import importlib.util
import math
import sys
from decimal import ROUND_HALF_EVEN, Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _score():
    spec = importlib.util.spec_from_file_location("c9_pick_score_for_page", HERE / "score.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


score = _score()

LEARNERS = score.LEARNERS
MEASURES = score.MEASURES
PRIMARY = score.PRIMARY
MARGIN = score.MARGIN
READ = "read"
ANSWERS = {
    "read": {"--tie": "later", "--no-pick": "default", "--default": "per-draw"},
    "tie-earlier": {"--tie": "earlier", "--no-pick": "default", "--default": "per-draw"},
    "tie-default": {"--tie": "default", "--no-pick": "default", "--default": "per-draw"},
    "default-seed42": {"--tie": "later", "--no-pick": "default", "--default": "seed42"},
}
# The answers whose default and twin rows must equal the read answer's: only the inner pick's tie differs.
SAME_DEFAULT = ("tie-earlier", "tie-default")
SEEDS, TWIN, STUDY_SEED = (42, 43, 44), (45, 46, 47), 42
SPLITS = ("lineage", "study")
ARMS = {"lineage": ("default", "inner", "twin"), "study": ("default", "inner")}
PAIRS = [("lineage", "inner"), ("lineage", "twin"), ("study", "inner")]
SCORE_FILES = ("datasets.tsv", "draws.tsv", "picks.tsv", "options.txt")
PICK_COLUMNS = ["dataset", "split", "seed", "folds", "ties", "no_pick", "same_as_default"] + ["picks_" + l for l in LEARNERS]
GATES = [
    (1, "answers", "the reading holds reading.txt and the four answers, each with count/count.tsv, count/conditions.tsv "
                   "and score/"),
    (2, "species", "each answer's score/ holds exactly the species of l.tsv's first run, each with its four files"),
    (3, "options", "each species' options.txt holds its answer's options"),
    (4, "datasets", "each datasets.tsv holds exactly its species' first-run datasets, none in M: the default, the inner "
                    "pick and the twin on the lineage split, and the default and the inner pick on the study split for "
                    "those in S; species and antibiotic equal l.tsv's"),
    (5, "draws", "each draws.tsv holds each dataset's draws at seeds 42, 43 and 44, its twin at 45, 46 and 47 and its "
                 "study run at seed 42, and nothing else; datasets.tsv is their mean; picks.tsv is their inner rows"),
    (6, "beside", "the tie answers' default and twin rows equal the read answer's"),
    (7, "count", "each answer's count.tsv and conditions.tsv are what score.py --count writes from its datasets.tsv files"),
    (8, "conditions", "the inner pick's conditions 1 and 2 as read here off count.tsv, and condition 3 off the twin's "
                      "row, equal conditions.tsv's"),
]
MEASURE_WORDS = {"balanced_accuracy": "balanced accuracy", "recall": "recall", "specificity": "specificity",
                 "vme": "very major error rate", "me": "major error rate", "log_loss": "log loss"}
ANSWER_WORDS = {"read": "read: ties to the later learner, each draw's own default",
                "tie-earlier": "beside: ties to the earlier learner",
                "tie-default": "beside: ties to the default's learner",
                "default-seed42": "beside: the seed-42 default at every draw"}


class Stop(Exception):
    pass


# ---------- reading files ----------

def load(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def num(x):
    return None if x in ("", None) else float(x)


def whole(x):
    return None if x in ("", None) else int(x)


def r12(x):
    return None if x is None else round(x, 12)


def first_run(lrows):
    return [r for r in lrows if r["first_run"] == "1"]


def cell_of(rows):
    return {(r["dataset"], r["split"], r["arm"]): r for r in rows}


# ---------- gates ----------

def expected_rows(datasets):
    """(dataset, split, arm) of datasets.tsv for the first-run rows of one species."""
    out = set()
    for r in datasets:
        for arm in ARMS["lineage"]:
            out.add((r["dataset"], "lineage", arm))
        if r["in_S"] == "1":
            for arm in ARMS["study"]:
                out.add((r["dataset"], "study", arm))
    return out


def expected_draws(datasets):
    """(dataset, split, seed, role, arm) of draws.tsv for the first-run rows of one species."""
    out = set()
    for r in datasets:
        for s in SEEDS:
            out.add((r["dataset"], "lineage", s, "draw", "default"))
            out.add((r["dataset"], "lineage", s, "draw", "inner"))
        for s in TWIN:
            out.add((r["dataset"], "lineage", s, "twin", "twin"))
        if r["in_S"] == "1":
            out.add((r["dataset"], "study", STUDY_SEED, "study", "default"))
            out.add((r["dataset"], "study", STUDY_SEED, "study", "inner"))
    return out


def typed_draw(r):
    """A draws.tsv row as score.dataset_rows takes it."""
    out = dict(r)
    out["seed"] = int(r["seed"])
    out["scored"] = int(r["scored"])
    for m in MEASURES:
        out[m] = num(r.get(m))
    for k in ("n", "tp", "fp", "tn", "fn"):
        out[k] = whole(r.get(k)) or 0
    return out


def same_cell(a, b):
    if a == b:
        return True
    try:
        return abs(float(a) - float(b)) <= 2e-12
    except ValueError:
        return False


def read_answer(folder, species):
    """{slug: {file: rows}} of one answer's score/ and its count files."""
    out = {}
    for slug in species:
        sp = folder / "score" / slug
        out[slug] = {name: load(sp / name) for name in ("datasets.tsv", "draws.tsv", "picks.tsv")}
        out[slug]["options"] = (sp / "options.txt").read_text()
    return out


def options_text(opts):
    return "".join(f"{k} {v}\n" for k, v in opts.items())


def check_reading(reading, lrows):
    """Every gate of GATES, in order. Returns (answers, count rows, condition rows); raises Stop naming the gate."""
    reading = Path(reading)
    run = first_run(lrows)
    species = sorted({r["species_slug"] for r in run})
    by_species = {s: [r for r in run if r["species_slug"] == s] for s in species}
    names = {r["dataset"]: (r["species"], r["antibiotic"]) for r in lrows}

    # 1 answers
    if not (reading / "reading.txt").is_file():
        raise Stop("gate 1 answers: no reading.txt")
    for a in ANSWERS:
        for need in ("count/count.tsv", "count/conditions.tsv"):
            if not (reading / a / need).is_file():
                raise Stop(f"gate 1 answers: {a}/{need} is missing")
        if not (reading / a / "score").is_dir():
            raise Stop(f"gate 1 answers: {a}/score is missing")
    extra = sorted(p.name for p in reading.iterdir() if p.is_dir() and p.name not in ANSWERS)
    if extra:
        raise Stop(f"gate 1 answers: folders that are no answer: {', '.join(extra)}")

    # 2 species; from here on every gate and table works over the species found.
    found = {}
    for a in ANSWERS:
        got = sorted(p.name for p in (reading / a / "score").iterdir() if p.is_dir())
        if got != species:
            missing, more = sorted(set(species) - set(got)), sorted(set(got) - set(species))
            raise Stop(f"gate 2 species: {a} has {len(got)} species of {len(species)}; missing {missing}, more {more}")
        for slug in got:
            for name in SCORE_FILES:
                if not (reading / a / "score" / slug / name).is_file():
                    raise Stop(f"gate 2 species: {a}/score/{slug}/{name} is missing")
        found[a] = got
    species = found[READ]
    if any(found[a] != species for a in ANSWERS):
        raise Stop("gate 2 species: the answers hold different species")
    answers = {a: read_answer(reading / a, species) for a in ANSWERS}

    # 3 options
    for a, opts in ANSWERS.items():
        for slug in species:
            if answers[a][slug]["options"] != options_text(opts):
                raise Stop(f"gate 3 options: {a}/score/{slug}/options.txt is not {' '.join(f'{k} {v}' for k, v in opts.items())}")

    # 4 datasets
    for a in ANSWERS:
        for slug in species:
            rows = answers[a][slug]["datasets.tsv"]
            keys = [(r["dataset"], r["split"], r["arm"]) for r in rows]
            want = expected_rows(by_species[slug])
            if len(set(keys)) != len(keys) or set(keys) != want:
                more = sorted({k[0] for k in set(keys) - want})
                less = sorted({k[0] for k in want - set(keys)})
                raise Stop(f"gate 4 datasets: {a}/score/{slug}/datasets.tsv holds other rows than the first run's; "
                           f"more {more}, less {less}")
            for r in rows:
                if (r["species"], r["antibiotic"]) != names[r["dataset"]] or r["species_slug"] != slug:
                    raise Stop(f"gate 4 datasets: {a}/{slug}: {r['dataset']} is named otherwise than in l.tsv")

    # 5 draws
    for a in ANSWERS:
        for slug in species:
            files = answers[a][slug]
            draws = files["draws.tsv"]
            keys = [(r["dataset"], r["split"], int(r["seed"]), r["role"], r["arm"]) for r in draws]
            if len(set(keys)) != len(keys) or set(keys) != expected_draws(by_species[slug]):
                raise Stop(f"gate 5 draws: {a}/score/{slug}/draws.tsv holds other draws than the first run's")
            again = score.dataset_rows([typed_draw(r) for r in draws])
            again_at = {(r["dataset"], r["split"], r["arm"]): r for r in again}
            for r in files["datasets.tsv"]:
                g = again_at[(r["dataset"], r["split"], r["arm"])]
                for c in score.DATASET_COLUMNS:
                    if not same_cell(r[c], score.fmt(g.get(c))):
                        raise Stop(f"gate 5 draws: {a}/{slug}: {r['dataset']} {r['split']} {r['arm']}: {c} is not the "
                                   "mean of its draws")
            inner = [{k: r.get(k, "") for k in PICK_COLUMNS} for r in draws if r["arm"] == "inner"]
            if inner != files["picks.tsv"]:
                raise Stop(f"gate 5 draws: {a}/score/{slug}/picks.tsv is not the inner rows of draws.tsv")

    # 6 beside
    for a in SAME_DEFAULT:
        for slug in species:
            mine = cell_of(answers[a][slug]["datasets.tsv"])
            read = cell_of(answers[READ][slug]["datasets.tsv"])
            for k, r in read.items():
                if k[2] != "inner" and mine[k] != r:
                    raise Stop(f"gate 6 beside: {a}: {k[0]} {k[1]} {k[2]} differs from the read answer's")

    # 7 count
    counts, conds = {}, {}
    for a in ANSWERS:
        rows = [r for slug in species for r in answers[a][slug]["datasets.tsv"]]
        counted = score.count(rows)
        want_count = [{c: score.fmt(r.get(c)) for c in score.COUNT_COLUMNS} for r in counted]
        want_cond = [{c: score.fmt(r.get(c)) for c in ("condition", "reads", "met")} for r in score.conditions(counted)]
        got_count = load(reading / a / "count" / "count.tsv")
        got_cond = load(reading / a / "count" / "conditions.tsv")
        if got_count != want_count:
            raise Stop(f"gate 7 count: {a}/count/count.tsv is not score.py's count of its datasets.tsv files")
        if got_cond != want_cond:
            raise Stop(f"gate 7 count: {a}/count/conditions.tsv is not score.py's conditions of its count")
        counts[a], conds[a] = got_count, got_cond

    # 8 conditions
    for a in ANSWERS:
        at = {(r["split"], r["pair"], r["measure"]): r for r in counts[a]}
        met = {int(r["condition"]): r["met"] for r in conds[a]}
        inner = meets(at[("lineage", "inner-default", PRIMARY)])
        twin = meets(at[("lineage", "twin-default", PRIMARY)])
        if (str(int(inner[0])), str(int(inner[1])), str(int(not (twin[0] or twin[1])))) != (met[1], met[2], met[3]):
            raise Stop(f"gate 8 conditions: {a}: conditions 1 to 3 read here differ from conditions.tsv")
    return answers, counts, conds


def meets(row):
    """(condition 1, condition 2) of one lineage primary row of count.tsv, as score.conditions reads them."""
    lo, holm = num(row["interval_lo"]), num(row["holm_p"])
    c1 = lo is not None and lo > 0
    c2 = holm is not None and holm < score.HOLM_ALPHA and int(row["better"]) > int(row["worse"])
    return c1, c2


# ---------- the page's tables ----------

def page_conditions(conds):
    rows = []
    for i, r in enumerate(conds[READ]):
        row = {"condition": r["condition"], "reads": r["reads"]}
        for a in ANSWERS:
            row[a] = conds[a][i]["met"]
        rows.append(row)
    return rows


def page_count(counts):
    rows = []
    for a in ANSWERS:
        for r in counts[a]:
            row = {"answer": a, **r, "meets_1": "", "meets_2": ""}
            if r["split"] == "lineage" and r["measure"] == PRIMARY:
                c1, c2 = meets(r)
                row["meets_1"], row["meets_2"] = int(c1), int(c2)
            rows.append(row)
    return rows


def diffs(rows, split, arm, measure=PRIMARY):
    """[(dataset, difference)] of one split and pair, rounded to 12 decimals as score.count takes them; and the
    count of datasets not scored."""
    mine = [r for r in rows if r["split"] == split]
    by = {(r["dataset"], r["arm"]): r for r in mine}
    out, unscored = [], 0
    for d in sorted({r["dataset"] for r in mine if r["arm"] == arm}):
        a, b = by[(d, arm)], by.get((d, "default"))
        va, vb = num(a[measure]), (num(b[measure]) if b else None)
        if b is None or a["scored"] != "1" or b["scored"] != "1" or va is None or vb is None:
            unscored += 1
            continue
        out.append((d, r12(va - vb)))
    return out, unscored


def tally(values, measure=PRIMARY):
    up = -1 if measure in score.LOWER_IS_BETTER else 1
    better = sum(1 for x in values if up * x > MARGIN)
    worse = sum(1 for x in values if up * x < -MARGIN)
    return better, worse, len(values) - better - worse


def median(values):
    import numpy as np
    return r12(float(np.median(values))) if values else None


def page_species(answers, lrows):
    run = first_run(lrows)
    out = []
    for slug in sorted(answers[READ]):
        mine = [r for r in run if r["species_slug"] == slug]
        rows = answers[READ][slug]["datasets.tsv"]
        row = {"species": mine[0]["species"], "species_slug": slug, "datasets": len(mine),
               "study_datasets": sum(1 for r in mine if r["in_S"] == "1")}
        for split, arm in PAIRS:
            got, unscored = diffs(rows, split, arm)
            values = [x for _, x in got]
            b, w, s = tally(values)
            key = f"{split}_{arm}"
            row.update({f"{key}_entered": len(values), f"{key}_not_scored": unscored, f"{key}_better": b,
                        f"{key}_worse": w, f"{key}_same": s, f"{key}_median": median(values)})
        out.append(row)
    return out


SPECIES_COLUMNS = ["species", "species_slug", "datasets", "study_datasets"] + [
    f"{split}_{arm}_{c}" for split, arm in PAIRS for c in ("entered", "not_scored", "better", "worse", "same", "median")]


def page_datasets(answers):
    out = []
    for slug in sorted(answers[READ]):
        files = answers[READ][slug]
        cells = cell_of(files["datasets.tsv"])
        draws = files["draws.tsv"]
        for d in sorted({r["dataset"] for r in files["datasets.tsv"]}):
            for split in SPLITS:
                if (d, split, "default") not in cells:
                    continue
                dflt = cells[(d, split, "default")]
                inner = cells[(d, split, "inner")]
                mine = sorted((r for r in draws if r["dataset"] == d and r["split"] == split),
                              key=lambda r: int(r["seed"]))
                defaults = [r["learner"] for r in mine if r["arm"] == "default" and r["scored"] == "1"]
                picks = [r for r in mine if r["arm"] == "inner" and r["scored"] == "1"]
                folds = sum(int(r["folds"]) for r in picks)
                no_pick = sum(int(r["no_pick"]) for r in picks)
                same = sum(int(r["same_as_default"]) for r in picks)
                row = {"dataset": d, "species": dflt["species"], "antibiotic": dflt["antibiotic"], "split": split,
                       "draws": dflt["draws"], "default_draws_scored": dflt["scored_draws"],
                       "inner_draws_scored": inner["scored_draws"], "default_learners": ",".join(defaults),
                       "default": num(dflt[PRIMARY]), "inner": num(inner[PRIMARY]),
                       "default_recall": num(dflt["recall"]), "inner_recall": num(inner["recall"]),
                       "folds": folds, "folds_no_pick": no_pick,
                       "folds_tied": sum(int(r["ties"]) for r in picks), "folds_same_as_default": same,
                       "folds_differing": folds - same - no_pick}
                row["difference"] = (r12(row["inner"] - row["default"])
                                     if row["inner"] is not None and row["default"] is not None else None)
                if split == "lineage":
                    twin = cells[(d, split, "twin")]
                    row["twin"] = num(twin[PRIMARY])
                    row["twin_difference"] = (r12(row["twin"] - row["default"])
                                              if row["twin"] is not None and row["default"] is not None else None)
                out.append(row)
    return out


DATASET_COLUMNS = ["dataset", "species", "antibiotic", "split", "draws", "default_draws_scored", "inner_draws_scored",
                   "default_learners", "default", "inner", "difference", "twin", "twin_difference", "default_recall",
                   "inner_recall", "folds", "folds_no_pick", "folds_tied", "folds_same_as_default", "folds_differing"]


def page_picks(answers):
    """Per answer and split: folds per learner the inner pick serves, and the default's, over the scored draws."""
    out = []
    for a in ANSWERS:
        for split in SPLITS:
            inner = {"answer": a, "split": split, "arm": "inner", "draws": 0, "folds": 0, "no_pick": 0, "tied": 0,
                     "same_as_default": 0, **{l: 0 for l in LEARNERS}}
            dflt = {"answer": a, "split": split, "arm": "default", "draws": 0, "folds": 0, **{l: 0 for l in LEARNERS}}
            for slug in sorted(answers[a]):
                draws = answers[a][slug]["draws.tsv"]
                learner_of = {(r["dataset"], r["seed"]): r["learner"] for r in draws
                              if r["split"] == split and r["arm"] == "default" and r["scored"] == "1"}
                for r in draws:
                    if r["split"] != split or r["arm"] != "inner" or r["scored"] != "1":
                        continue
                    f = int(r["folds"])
                    inner["draws"] += 1
                    inner["folds"] += f
                    inner["no_pick"] += int(r["no_pick"])
                    inner["tied"] += int(r["ties"])
                    inner["same_as_default"] += int(r["same_as_default"])
                    for l in LEARNERS:
                        inner[l] += int(r["picks_" + l])
                    dflt["draws"] += 1
                    dflt["folds"] += f
                    dflt[learner_of[(r["dataset"], r["seed"])]] += f
            inner["differs"] = inner["folds"] - inner["same_as_default"] - inner["no_pick"]
            out += [dflt, inner]
    return out


PICK_TABLE_COLUMNS = ["answer", "split", "arm", "draws", "folds"] + LEARNERS + ["no_pick", "tied", "same_as_default",
                                                                               "differs"]


def page_pick_pairs(answers):
    """The read answer: per split and default's learner, the folds the inner pick serves to each learner."""
    out = []
    for split in SPLITS:
        rows = {l: {"split": split, "default_learner": l, "draws": 0, "folds": 0, "no_pick": 0, **{k: 0 for k in LEARNERS}}
                for l in LEARNERS}
        for slug in sorted(answers[READ]):
            draws = answers[READ][slug]["draws.tsv"]
            learner_of = {(r["dataset"], r["seed"]): r["learner"] for r in draws
                          if r["split"] == split and r["arm"] == "default" and r["scored"] == "1"}
            for r in draws:
                if r["split"] != split or r["arm"] != "inner" or r["scored"] != "1":
                    continue
                row = rows[learner_of[(r["dataset"], r["seed"])]]
                row["draws"] += 1
                row["folds"] += int(r["folds"])
                row["no_pick"] += int(r["no_pick"])
                for k in LEARNERS:
                    row[k] += int(r["picks_" + k])
        out += [rows[l] for l in LEARNERS]
    return out


PAIR_COLUMNS = ["split", "default_learner", "draws", "folds", "no_pick"] + LEARNERS


def page_pooled(answers):
    """The read answer: per split and arm, the calls summed over the scored datasets and their draws, with Rules'
    measures and Wilson intervals. A genome is counted once a draw."""
    out = []
    for split in SPLITS:
        for arm in ARMS[split]:
            c = {"n": 0, "tp": 0, "fp": 0, "tn": 0, "fn": 0, "loss_sum": 0.0}
            datasets = scored = 0
            for slug in sorted(answers[READ]):
                for r in answers[READ][slug]["datasets.tsv"]:
                    if r["split"] != split or r["arm"] != arm:
                        continue
                    datasets += 1
                    if r["scored"] != "1":
                        continue
                    scored += 1
                    for k in ("n", "tp", "fp", "tn", "fn"):
                        c[k] += int(r[k])
            m = score.measures(c)
            row = {"split": split, "arm": arm, "datasets": datasets, "scored": scored}
            row.update({k: m[k] for k in ("n", "tp", "fp", "tn", "fn", "resistant", "susceptible", "balanced_accuracy")})
            for k in ("recall", "specificity", "vme", "me"):
                row.update({k: m[k], f"{k}_lo": m[f"{k}_lo"], f"{k}_hi": m[f"{k}_hi"]})
            out.append(row)
    return out


POOLED_COLUMNS = ["split", "arm", "datasets", "scored", "n", "resistant", "susceptible", "tp", "fp", "tn", "fn",
                  "balanced_accuracy"] + [c for k in ("recall", "specificity", "vme", "me") for c in (k, f"{k}_lo", f"{k}_hi")]


# ---------- writing ----------

README = """# The served-learner reading of 23 November 2026: the record

This directory holds the tables behind [docs/served-pick-result.md](../../docs/served-pick-result.md): plan block C9, the learner picked inside each outer training fold against the release's pick, read once on the first run's {datasets} datasets in {species} species, {study} of them with a study split. Every row is an aggregate over a dataset, a species, a split or a learner. No file holds a genome, a sequence or a feature.

`scripts/c9-pick/score.py` wrote the files under `count/`; `scripts/c9-pick/page_tables.py write` wrote every other file here from the run's own `score/` folders and the reading's three scorings beside them, after the eight gates of `gates.txt` passed. `reading.txt` lists the reading's commands, each with the sha256 of what it wrote. The rules of the reading were fixed at 20:36Z on 22 Nov 2026, before any C9 score was read (`docs/assumption-plan-corrections.md`, entry "C9 design").

| File | Holds |
|---|---|
| `gates.txt` | the eight gates the program passed before it wrote a table |
| `reading.txt` | the reading's commands and the sha256 of each one's output |
| `count/read/` | the count that decides: `--tie later --no-pick default --default per-draw` |
| `count/tie-earlier/`, `count/tie-default/`, `count/default-seed42/` | the counts beside it, each deciding nothing |
| `page-conditions.tsv` | the plan's change conditions 1 to 5 under each answer, as the scorer read them; 6 is the maintainer's |
| `page-count.tsv` | every count row of the four answers, with conditions 1 and 2 read on the lineage primary rows |
| `page-species.tsv` | per species, named in full: datasets, and better, worse and same on the primary measure per split and pair |
| `page-datasets.tsv` | per dataset and split, species and antibiotic in full: the default's learner per draw and figure, the inner pick's figure, the difference, the twin's, the draws scored and the folds |
| `page-picks.tsv` | per answer and split: the folds each learner serves under the default and under the inner pick, folds with no pick, tied and differing from the default's learner |
| `page-pick-pairs.tsv` | the read answer, per split and default's learner: the folds the inner pick serves to each learner |
| `page-pooled.tsv` | the read answer, per split and arm: the calls summed over datasets and draws, with Rules' measures and Wilson intervals |

A dataset's figure is the mean of its draws: seeds 42, 43 and 44 on the lineage split, the twin's 45, 46 and 47, and one study run at the release's seeds. The pooled calls count each genome once a draw. The default's pick reads the random split over every genome, the outer test folds' genomes included, so the default's served figure can carry the pick's optimism and the inner pick's cannot. A same or a worse is never read as showing that the release's pick is right.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/c9-pick/page_tables.py print` prints the page's tables from these files, and `tests/c9_pick/test_page.py` checks the page against them.
"""


def write_tsv(path, rows, columns):
    with open(path, "w", newline="") as handle:
        handle.write("\t".join(columns) + "\n")
        for r in rows:
            handle.write("\t".join(score.fmt(r.get(c)) for c in columns) + "\n")


def write(reading, l_path, out):
    reading, out = Path(reading), Path(out)
    lrows = load(l_path)
    if out.exists() and any(out.iterdir()):
        raise Stop(f"{out} is not empty")
    answers, counts, conds = check_reading(reading, lrows)
    out.mkdir(parents=True, exist_ok=True)
    for a in ANSWERS:
        (out / "count" / a).mkdir(parents=True)
        for name in ("count.tsv", "conditions.tsv"):
            (out / "count" / a / name).write_bytes((reading / a / "count" / name).read_bytes())
    (out / "reading.txt").write_bytes((reading / "reading.txt").read_bytes())
    write_tsv(out / "page-conditions.tsv", page_conditions(conds), ["condition", "reads"] + list(ANSWERS))
    write_tsv(out / "page-count.tsv", page_count(counts), ["answer"] + score.COUNT_COLUMNS + ["meets_1", "meets_2"])
    write_tsv(out / "page-species.tsv", page_species(answers, lrows), SPECIES_COLUMNS)
    write_tsv(out / "page-datasets.tsv", page_datasets(answers), DATASET_COLUMNS)
    write_tsv(out / "page-picks.tsv", page_picks(answers), PICK_TABLE_COLUMNS)
    write_tsv(out / "page-pick-pairs.tsv", page_pick_pairs(answers), PAIR_COLUMNS)
    write_tsv(out / "page-pooled.tsv", page_pooled(answers), POOLED_COLUMNS)
    with open(out / "gates.txt", "w") as handle:
        for n, name, words in GATES:
            handle.write(f"{n} {name}: {words}\n")
    run = first_run(lrows)
    (out / "README.md").write_text(README.format(datasets=len(run), species=len({r["species_slug"] for r in run}),
                                                 study=sum(1 for r in run if r["in_S"] == "1")))
    files = sorted(p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file())
    with open(out / "SHA256SUMS", "w") as handle:
        for name in files:
            handle.write(f"{hashlib.sha256((out / name).read_bytes()).hexdigest()}  {name}\n")
    return files


# ---------- printing ----------

def d3(x):
    if x in ("", None):
        return "none"
    value = Decimal(str(x)).quantize(Decimal("0.001"), rounding=ROUND_HALF_EVEN)
    if value == 0:
        value = Decimal("0.000")
    return str(value)


def pval(x):
    if x in ("", None):
        return "none"
    p = float(x)
    return "< 0.001" if p < 0.001 else d3(p)


def count_words(x):
    return f"{int(x):,}"


def of(part, whole_):
    return f"{count_words(part)} of {count_words(whole_)}"


def interval(lo, hi):
    return "none" if lo in ("", None) else f"{d3(lo)} to {d3(hi)}"


def yes(x):
    return {"1": "met", "0": "not met", "": "the maintainer's"}[x]


def table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return lines


PAIR_WORDS = {"inner-default": "inner pick less default", "twin-default": "twin less default"}


def print_conditions(rows):
    return table(["Condition", "Reads", "Read", "Ties to the earlier", "Ties to the default's", "Seed-42 default"],
                 [[r["condition"], r["reads"]] + [yes(r[a]) for a in ANSWERS] for r in rows])


def print_twin(rows):
    out = []
    for r in rows:
        if r["answer"] == READ and r["split"] == "lineage" and r["measure"] == PRIMARY:
            out.append([PAIR_WORDS[r["pair"]], yes(str(r["meets_1"])), yes(str(r["meets_2"])),
                        of(r["better"], r["entered"]), of(r["worse"], r["entered"]), of(r["same"], r["entered"]),
                        d3(r["median_difference"]), interval(r["interval_lo"], r["interval_hi"]), pval(r["sign_p"]),
                        pval(r["holm_p"])])
    return table(["Pair", "Condition 1", "Condition 2", "Better", "Worse", "Same", "Median", "Interval", "Sign p",
                  "Holm p"], out)


def print_count(rows):
    out = []
    for r in rows:
        if r["answer"] != READ:
            continue
        out.append([r["split"], PAIR_WORDS[r["pair"]], MEASURE_WORDS[r["measure"]], r["orientation"],
                    of(r["entered"], r["datasets"]), r["better"], r["worse"], r["same"], d3(r["median_difference"]),
                    interval(r["interval_lo"], r["interval_hi"]), r["species"]])
    return table(["Split", "Pair", "Measure", "Better is", "Entered", "Better", "Worse", "Same", "Median", "Interval",
                  "Species"], out)


def print_beside(rows, conds):
    out = []
    for a in ANSWERS:
        met = ", ".join(r["condition"] for r in conds if r[a] == "1") or "none"
        for r in rows:
            if r["answer"] != a or r["measure"] not in (PRIMARY, "recall"):
                continue
            out.append([ANSWER_WORDS[a], r["split"], PAIR_WORDS[r["pair"]], MEASURE_WORDS[r["measure"]],
                        r["better"], r["worse"], r["same"], d3(r["median_difference"]),
                        interval(r["interval_lo"], r["interval_hi"]), met])
    return table(["Answer", "Split", "Pair", "Measure", "Better", "Worse", "Same", "Median", "Interval",
                  "Conditions met"], out)


def print_pooled(rows):
    out = []
    for r in rows:
        out.append([r["split"], r["arm"], of(r["scored"], r["datasets"]), d3(r["balanced_accuracy"]),
                    f"{d3(r['recall'])} ({of(r['tp'], r['resistant'])})", interval(r["recall_lo"], r["recall_hi"]),
                    f"{d3(r['specificity'])} ({of(r['tn'], r['susceptible'])})",
                    f"{d3(r['vme'])} ({of(r['fn'], r['resistant'])})", f"{d3(r['me'])} ({of(r['fp'], r['susceptible'])})"])
    return table(["Split", "Arm", "Datasets scored", "Balanced accuracy", "Recall", "Recall interval", "Specificity",
                  "Very major error rate", "Major error rate"], out)


def print_picks(rows):
    out = []
    for r in rows:
        tail = ([r["no_pick"], r["tied"], r["same_as_default"], of(r["differs"], r["folds"])] if r["arm"] == "inner"
                else ["", "", "", ""])
        out.append([ANSWER_WORDS[r["answer"]], r["split"], r["arm"], r["draws"], count_words(r["folds"])]
                   + [count_words(r[l]) for l in LEARNERS] + tail)
    return table(["Answer", "Split", "Serves", "Draws", "Folds"] + [f"`{l}`" for l in LEARNERS]
                 + ["No pick", "Tied", "Same as default", "Differs"], out)


def print_pairs(rows):
    out = []
    for r in rows:
        if r["draws"] == "0":
            continue
        out.append([r["split"], f"`{r['default_learner']}`", r["draws"], count_words(r["folds"]), r["no_pick"]]
                   + [count_words(r[l]) for l in LEARNERS])
    return table(["Split", "Default's learner", "Draws", "Folds", "No pick"] + [f"Inner `{l}`" for l in LEARNERS], out)


def print_species(rows):
    out = []
    for r in rows:
        cells = [f"*{r['species']}*", r["datasets"], r["study_datasets"]]
        for split, arm in PAIRS:
            k = f"{split}_{arm}"
            cells.append(f"{r[k + '_better']} / {r[k + '_worse']} / {r[k + '_same']}")
            cells.append(d3(r[k + "_median"]))
        out.append(cells)
    return table(["Species", "Datasets", "In S", "Lineage inner better / worse / same", "Lineage inner median",
                  "Twin better / worse / same", "Twin median", "Study inner better / worse / same", "Study inner median"],
                 out)


def print_datasets(rows):
    out = []
    for r in rows:
        learners = ", ".join(f"`{l}`" for l in r["default_learners"].split(",")) if r["default_learners"] else "none"
        out.append([f"*{r['species']}*", r["antibiotic"], r["split"], learners, d3(r["default"]), d3(r["inner"]),
                    d3(r["difference"]), d3(r["twin_difference"]) if r["split"] == "lineage" else "",
                    of(r["inner_draws_scored"], r["draws"]), r["folds_no_pick"], of(r["folds_differing"], r["folds"])])
    return table(["Species", "Antibiotic", "Split", "Default's learner", "Default", "Inner pick", "Difference",
                  "Twin difference", "Draws scored", "No pick", "Folds differing"], out)


def tables(folder):
    folder = Path(folder)
    count = load(folder / "page-count.tsv")
    conds = load(folder / "page-conditions.tsv")
    return {
        "conditions": print_conditions(conds),
        "twin": print_twin(count),
        "count": print_count(count),
        "beside": print_beside(count, conds),
        "pooled": print_pooled(load(folder / "page-pooled.tsv")),
        "picks": print_picks(load(folder / "page-picks.tsv")),
        "pairs": print_pairs(load(folder / "page-pick-pairs.tsv")),
        "species": print_species(load(folder / "page-species.tsv")),
        "datasets": print_datasets(load(folder / "page-datasets.tsv")),
    }


def main(argv):
    if len(argv) == 4 and argv[0] == "write":
        try:
            files = write(argv[1], argv[2], argv[3])
        except Stop as e:
            print(f"stopped: {e}", file=sys.stderr)
            return 1
        print(f"gates: {len(GATES)} of {len(GATES)} pass; wrote {len(files) + 1} files")
        return 0
    if len(argv) == 2 and argv[0] == "print":
        for name, lines in tables(argv[1]).items():
            print(f"[{name}]")
            print("\n".join(lines))
            print()
        return 0
    sys.stderr.write(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
