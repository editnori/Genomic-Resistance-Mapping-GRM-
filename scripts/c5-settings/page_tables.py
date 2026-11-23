"""The tables of docs/settings-arms-result.md: plan block C5's one reading, as files and as the page prints them.

  python3 -I page_tables.py write <reading dir> <lists.tsv> <checks.txt> <out dir>
  python3 -I page_tables.py print <out dir>

<lists.tsv> is scripts/c5-settings/lists.tsv, the list of S: the first run
is its rows with first_run 1, the rows of S in L that are not in M.
<checks.txt> is the printout of scripts/c5-settings/run_checks.py over the
run's copies, made before any score was read. The reading directory holds:

  reading.txt            the reading's commands, each with the sha256 of what it wrote
  count/                 score.py --count's count.tsv, count-no-constant.tsv,
                         conditions.tsv and count-plain-best.tsv
  score/<species slug>/  one species' run as run_checks.py gather copies it:
                         score.py's checks.tsv, draws.tsv and datasets.tsv, the
                         run's binary.txt, folds.tsv (each default draw's and
                         study run's outer fold count), and choices/, the
                         driver's choices files; never score.log

`write` passes every gate of GATES before it writes a file, and stops at
the first that fails, naming it. It then writes <out dir> (absent or
empty): score.py's four count tables under count/, copied once the count
gate has held them to score.py's own count of the datasets.tsv files;
reading.txt and checks.txt, copied; the page's tables as page-*.tsv;
gates.txt, README.md and SHA256SUMS. Where score.py's tables already are
the page's (the count, the two counts beside it and the conditions), this
script checks and copies them, and its page-*.tsv files add only marks:
which rows decide, which stand in the Holm family, conditions 1 and 2 read
on each family row, and the cutoff arm's log loss rows marked equal by
construction. Of the choices files it writes counts alone: per learner the
folds that kept trial 0, a trial id's folds, a grid point's folds and the
folds where a trial failed; never a trial's settings. Every row is an
aggregate over a dataset, a species, an arm, a fold count or a learner: no
file holds a genome, a sequence or a feature.

`print` prints the page's tables from <out dir> alone, figures rounded to
3 decimals.
"""

import csv
import hashlib
import importlib.util
import json
import re
import statistics
import sys
from decimal import ROUND_HALF_EVEN, Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


score = _load("c5_settings_score_for_page", HERE / "score.py")
EXPECT = _load("c5_settings_expect_for_page", HERE / "expect.py")

C3 = score.C3
LEARNERS = list(score.LEARNERS)
MEASURES = list(C3.MEASURES)
PRIMARY = C3.PRIMARY
MARGIN = C3.MARGIN
DEFAULT, TWIN, PLAIN = score.DEFAULT, score.TWIN, score.PLAIN
CUTOFF = "cutoff"
ARMS = ("search", CUTOFF)
FAMILY = tuple(score.FAMILY)
# The arms the choices files hold beside the default: the family's two and the plain best.
CHOICE_ARMS = tuple(score.ARMS)
SEEDS, TWIN_SEEDS, STUDY_SEED = (42, 43, 44), (45, 46, 47), 42
GRID = list(score.GRID)
TRIALS = score.TRIALS
HALF = score.HALF
COUNT_FILES = ("count.tsv", "count-no-constant.tsv", "conditions.tsv", "count-plain-best.tsv")
SPECIES_FILES = ("checks.tsv", "draws.tsv", "datasets.tsv", "binary.txt", "folds.tsv")
REFUSED_WHY = "constant model: refused for no candidate"
MISSING_WHY = "missing draw"
CHOICE_NAME = re.compile(r"(?P<dataset>[a-z0-9_]+__[a-z0-9_]+)\.(?P<scheme>lineage|study-whole)-s(?P<seed>\d+)"
                         r"\.f(?P<fold>\d+)\.choices\.json")
SYNTHETIC_LOSS = (0.0126, 0.0088, 0.0163)
EQUAL_BY_CONSTRUCTION = "equal by construction: the cutoff arm moves no score"
TWIN_NOTE = "the twin: read for condition 3 alone, never as an arm"
GATES = [
    (1, "files", "the reading holds reading.txt, count/ with score.py --count's four tables, and score/; no score.log "
                 "is in it; the run's checks.txt says every check passed"),
    (2, "species", "score/ holds each species of lists.tsv's first run once and no other, each with its five files and "
                   "choices/, and their datasets sum to the first run's"),
    (3, "datasets", "each datasets.tsv holds exactly its species' first-run datasets, none in M: the default, the "
                    "search, the cutoff arm, the plain best and the twin on the lineage split, and the default, the "
                    "search, the cutoff arm and the plain best on the study split, each for the five learners; species "
                    "and antibiotic equal lists.tsv's"),
    (4, "study", "every study row is a study run at seed 42, and no study row stands on the lineage split"),
    (5, "twin", "the twin stands on the lineage split at seeds 45, 46 and 47 alone, no choices file is at a twin seed, "
                "and conditions.tsv reads the search and the cutoff arm as the arms and the twin as none"),
    (6, "plain", "the plain best, no anchor, stands in count-plain-best.tsv alone, with no Holm p, and in no condition "
                 "and no other count"),
    (7, "missing", "a draw not scored has no figure, and a dataset with a draw not scored is not scored and has no "
                   "figure: a missing draw is never read as zero"),
    (8, "constant", "each constant-model draw is an arm's study run refused for no candidate, kept, scored and marked, "
                    "and count-no-constant.tsv counts none"),
    (9, "draws", "each draws.tsv holds each dataset's draws at seeds 42, 43 and 44, the twin's at 45, 46 and 47 and "
                 "the study run at seed 42, each for the five learners, and nothing else; datasets.tsv is their mean"),
    (10, "checks", "every check of each checks.tsv passed"),
    (11, "cutoff", "the cutoff arm's log loss is the default's on every dataset, split and learner, so its log loss "
                   "rows count the same everywhere; scm's cutoff arm is the default on every fold: its cutoff 0.5 and "
                   "its calls the default's in each choices file, and its figures the default's in datasets.tsv"),
    (12, "choices", "each draw and study run of the default whose arms are scored has one choices file per outer fold, "
                    "folds 0 on, as many as folds.tsv gives its default's result, of its dataset, split, seed and fold, with the sampler seed, gamma, grid, 30 random "
                    "trials and the five learners; each pick a trial id from 0 to 30 by the default-anchored rule, "
                    "each cutoff a grid point; and a draw with no choices file is a missing draw"),
    (13, "settings", "the counts of what the search chose hold trial ids, grid points and fold counts alone: no "
                     "trial's settings and no row per genome"),
    (14, "binary", "each binary.txt names a commit, a binary and a code sha256; the page names each binary with its "
                   "species unless every species ran one code sha256"),
    (15, "conditions", "conditions 1 and 2, 3 and 4 of ruling 15 for the search and the cutoff arm, as read here off "
                       "count.tsv, equal conditions.tsv's"),
    (16, "count", "count.tsv, count-no-constant.tsv, conditions.tsv and count-plain-best.tsv are what score.py --count "
                  "writes from the datasets.tsv files"),
]
GATE_NAMES = {n: name for n, name, _ in GATES}
MEASURE_WORDS = {"balanced_accuracy": "balanced accuracy", "recall": "recall", "specificity": "specificity",
                 "vme": "very major error rate", "me": "major error rate", "log_loss": "log loss"}
ARM_WORDS = {"search": "search", CUTOFF: "cutoff", TWIN: "twin", PLAIN: "plain best, no anchor"}
# The columns of the tables of what the search chose: trial ids, grid points and fold counts alone.
CHOICE_COLUMNS = ["split", "learner", "folds", "kept_trial_0", "picked_other", "plain_best_trial_0",
                  "plain_best_is_pick", "plain_best_differs", "plain_best_refitted", "trial_0_failed",
                  "every_trial_failed", "failed_trials", "search_served_default", "plain_served_default",
                  "cutoff_half", "cutoff_moved"]
TRIAL_COLUMNS = ["split", "learner", "trial", "picked", "plain_best"]
CUTOFF_COLUMNS = ["split", "learner", "cutoff", "folds"]


class Stop(Exception):
    pass


def stop(n, text):
    raise Stop(f"gate {n} {GATE_NAMES[n]}: {text}")


# ---------- reading files ----------

def load(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def num(x):
    return None if x in ("", None) else float(x)


def r12(x):
    return None if x is None else round(x, 12)


def fmt(x):
    return score.C9.fmt(x)


def first_run(lrows):
    """The first run: S's rows with first_run 1, none of them in M."""
    return [r for r in lrows if r["first_run"] == "1" and r["in_M"] == "0"]


def in_m(lrow):
    return lrow["in_M"] == "1"


def split_of(scheme):
    return "lineage" if scheme == "lineage" else "study"


def expected(rows, slug):
    """expect.py's rows for one species' first-run rows, as dicts: every first-run dataset is in S and has a study
    run."""
    datasets = [r["dataset"] for r in rows]
    out = EXPECT.rows(slug, datasets, datasets, list(SEEDS), list(TWIN_SEEDS))
    return [dict(zip(EXPECT.COLUMNS, map(str, r))) for r in out]


def typed_draw(r):
    """A draws.tsv row as score.py's dataset_rows takes it."""
    out = dict(r)
    out["seed"] = int(r["seed"])
    out["scored"] = int(r["scored"])
    out["constant"] = int(r["constant"] or 0)
    for m in MEASURES:
        out[m] = num(r.get(m))
    for k in ("n", "tp", "fp", "tn", "fn"):
        out[k] = int(r[k]) if r.get(k) not in ("", None) else 0
    return out


def same_cell(a, b):
    if a == b:
        return True
    try:
        return abs(float(a) - float(b)) <= 2e-12
    except ValueError:
        return False


def fmt_rows(rows, columns):
    return [{c: fmt(r.get(c)) for c in columns} for r in rows]


def arm_of(pair):
    return pair[: -len(f"-{DEFAULT}")]


def read_choices(folder):
    """{(dataset, scheme, seed, fold): doc} of one species' choices/, each file's name read as the driver writes it."""
    out = {}
    for p in sorted(folder.iterdir()):
        m = CHOICE_NAME.fullmatch(p.name)
        if not m or not p.is_file():
            raise Stop(f"gate 12 choices: {folder.name}/{p.name} is not a choices file's name")
        out[(m["dataset"], m["scheme"], int(m["seed"]), int(m["fold"]))] = load_json(p)
    return out


# ---------- gates ----------

def check_reading(reading, lrows, checks_path):
    """Every gate of GATES, in order. Returns what the tables read; raises Stop naming the first gate that fails."""
    reading = Path(reading)
    run = first_run(lrows)
    lists_of = {r["dataset"]: r for r in lrows}

    # 1 files
    if not (reading / "reading.txt").is_file():
        stop(1, "no reading.txt")
    for name in COUNT_FILES:
        if not (reading / "count" / name).is_file():
            stop(1, f"count/{name} is missing")
    if not (reading / "score").is_dir():
        stop(1, "score/ is missing")
    logs = sorted(p.relative_to(reading).as_posix() for p in reading.rglob("score.log"))
    if logs:
        stop(1, f"score.log stays out of the record, and the reading holds {', '.join(logs)}")
    lines = Path(checks_path).read_text(encoding="utf-8").splitlines()
    tail = re.fullmatch(r"CHECKS (\d+) of (\d+)", lines[-1]) if lines else None
    if not tail or tail.group(1) != tail.group(2) or any(line.startswith("FAIL") for line in lines):
        stop(1, "the run's checks.txt does not say that every check passed")

    # 2 species
    species = sorted({r["species_slug"] for r in run})
    found = sorted(p.name for p in (reading / "score").iterdir() if p.is_dir())
    strays = sorted(p.name for p in (reading / "score").iterdir() if not p.is_dir())
    if found != species or strays:
        stop(2, f"score/ holds {len(found)} species of the first run's {len(species)}; missing "
                f"{sorted(set(species) - set(found))}, more {sorted(set(found) - set(species))}, files {strays}")
    for slug in species:
        for name in SPECIES_FILES:
            if not (reading / "score" / slug / name).is_file():
                stop(2, f"score/{slug}/{name} is missing")
        if not (reading / "score" / slug / "choices").is_dir():
            stop(2, f"score/{slug}/choices/ is missing")
    files = {slug: {name: load(reading / "score" / slug / name) for name in SPECIES_FILES if name.endswith(".tsv")}
             for slug in species}
    by_species = {slug: [r for r in run if r["species_slug"] == slug] for slug in species}
    held = sum(len({r["dataset"] for r in files[slug]["datasets.tsv"]}) for slug in species)
    if held != len(run):
        stop(2, f"the species' datasets.tsv files hold {held} datasets, the first run {len(run)}")
    counts = {name: load(reading / "count" / name) for name in COUNT_FILES}
    expect = {slug: expected(by_species[slug], slug) for slug in species}

    # 3 datasets
    for slug in species:
        rows = files[slug]["datasets.tsv"]
        keys = [(r["dataset"], r["split"], r["arm"], r["learner"]) for r in rows]
        want = {(e["dataset"], e["split"], e["arm"], l) for e in expect[slug] for l in LEARNERS}
        if len(set(keys)) != len(keys) or set(keys) != want:
            more = sorted({k[0] for k in set(keys) - want})
            less = sorted({k[0] for k in want - set(keys)})
            stop(3, f"score/{slug}/datasets.tsv holds other rows than the first run's; more {more}, less {less}")
        for r in rows:
            mine = lists_of.get(r["dataset"])
            if mine is None or in_m(mine):
                stop(3, f"{r['dataset']} is a dataset of M or of no list")
            if (r["species_slug"], r["species"], r["antibiotic"]) != (slug, mine["species"], mine["antibiotic"]):
                stop(3, f"{r['dataset']} is named otherwise than in lists.tsv")

    # 4 study
    for slug in species:
        for r in files[slug]["draws.tsv"]:
            if (r["split"] == "study") != (r["role"] == "study"):
                stop(4, f"score/{slug}/draws.tsv: {r['dataset']} {r['arm']} seed {r['seed']}: a {r['role']} row on "
                        f"the {r['split']} split")
            if r["role"] == "study" and int(r["seed"]) != STUDY_SEED:
                stop(4, f"score/{slug}/draws.tsv: {r['dataset']} {r['arm']}: a study run at seed {r['seed']}")

    # 5 twin
    for slug in species:
        for r in files[slug]["draws.tsv"]:
            is_twin = r["arm"] == TWIN
            if is_twin != (r["role"] == "twin") or (is_twin and (r["split"] != "lineage"
                                                                 or int(r["seed"]) not in TWIN_SEEDS)):
                stop(5, f"score/{slug}/draws.tsv: {r['dataset']} {r['arm']} as {r['role']} on the {r['split']} split "
                        f"at seed {r['seed']}")
    choices = {slug: read_choices(reading / "score" / slug / "choices") for slug in species}
    for slug in species:
        at_twin = sorted({k[2] for k in choices[slug] if k[2] in TWIN_SEEDS})
        if at_twin:
            stop(5, f"score/{slug}/choices/ holds a choices file at the twin's seed {at_twin[0]}")
    named = sorted({r["arm"] for r in counts["conditions.tsv"]})
    if named != sorted(ARMS):
        stop(5, f"conditions.tsv reads {', '.join(named)} as the arms, not the search and the cutoff arm")
    for name in ("count.tsv", "count-no-constant.tsv"):
        for r in counts[name]:
            if arm_of(r["pair"]) == TWIN and r["split"] != "lineage":
                stop(5, f"{name} counts the twin on the {r['split']} split")

    # 6 plain
    for name in ("count.tsv", "count-no-constant.tsv"):
        pairs = sorted({r["pair"] for r in counts[name]} - {f"{a}-{DEFAULT}" for a in FAMILY})
        if pairs:
            stop(6, f"{name} counts {', '.join(pairs)}, which stands in no count of the family")
    for r in counts["count-plain-best.tsv"]:
        if r["pair"] != f"{PLAIN}-{DEFAULT}" or r["holm_p"] != "":
            stop(6, f"count-plain-best.tsv holds {r['pair']} with Holm p '{r['holm_p']}'")
    if not counts["count-plain-best.tsv"]:
        stop(6, "count-plain-best.tsv is empty")
    if any(r["arm"] == PLAIN for r in counts["conditions.tsv"]):
        stop(6, "conditions.tsv reads the plain best as an arm")

    # 7 missing
    for slug in species:
        for r in files[slug]["draws.tsv"]:
            if r["scored"] != "1" and any(r[m] != "" for m in MEASURES):
                stop(7, f"score/{slug}/draws.tsv: {r['dataset']} {r['split']} {r['arm']} seed {r['seed']} "
                        f"{r['learner']}: a draw not scored has a figure")
        for r in files[slug]["datasets.tsv"]:
            whole_draws = r["scored_draws"] == r["draws"]
            if r["scored"] != str(int(whole_draws)) or (r["scored"] != "1" and any(r[k] != "" for k in MEASURES)):
                stop(7, f"score/{slug}/datasets.tsv: {r['dataset']} {r['split']} {r['arm']} {r['learner']}: "
                        f"{r['scored_draws']} of {r['draws']} draws scored, read as scored {r['scored']}")

    # 8 constant
    for slug in species:
        for r in files[slug]["draws.tsv"]:
            marked = r["constant"] == "1"
            if marked != (r["why"] == REFUSED_WHY) or (marked and (r["split"] != "study" or r["arm"] in (DEFAULT, TWIN)
                                                                   or r["scored"] != "1")):
                stop(8, f"score/{slug}/draws.tsv: {r['dataset']} {r['split']} {r['arm']} {r['learner']}: a constant "
                        f"model marked {r['constant']} for '{r['why']}'")
    for r in counts["count-no-constant.tsv"]:
        if r["constant"] not in ("", "0"):
            stop(8, f"count-no-constant.tsv counts {r['constant']} constant-model datasets for {r['pair']}")

    # 9 draws
    for slug in species:
        draws = files[slug]["draws.tsv"]
        keys = [(r["dataset"], r["split"], r["arm"], int(r["seed"]), r["role"], r["learner"]) for r in draws]
        want = {(e["dataset"], e["split"], e["arm"], int(e["seed"]), e["role"], l) for e in expect[slug]
                for l in LEARNERS}
        if len(set(keys)) != len(keys) or set(keys) != want:
            stop(9, f"score/{slug}/draws.tsv holds other draws than the first run's: {len(set(keys) - want)} more, "
                    f"{len(want - set(keys))} less, {len(keys) - len(set(keys))} twice")
        again = C3.dataset_rows([typed_draw(r) for r in draws])
        again_at = {(r["dataset"], r["split"], r["arm"], r["learner"]): r for r in again}
        for r in files[slug]["datasets.tsv"]:
            g = again_at[(r["dataset"], r["split"], r["arm"], r["learner"])]
            for c in C3.DATASET_COLUMNS:
                if not same_cell(r[c], fmt(g.get(c))):
                    stop(9, f"score/{slug}: {r['dataset']} {r['split']} {r['arm']} {r['learner']}: {c} is not the "
                            "mean of its draws")

    # 10 checks
    run_checks = []
    for slug in species:
        rows = files[slug]["checks.tsv"]
        bad = [r for r in rows if r["passed"] != "1"]
        if not rows or bad:
            stop(10, f"score/{slug}/checks.tsv: {len(bad)} of {len(rows)} checks did not pass")
        run_checks += rows

    # 11 cutoff
    agree = cutoff_checks(files, choices)
    for name, passed, made in agree:
        if passed != made:
            stop(11, f"{name}: {passed} of {made}")
    for r in counts["count.tsv"]:
        if arm_of(r["pair"]) == CUTOFF and r["measure"] == "log_loss" and (r["better"], r["worse"]) != ("0", "0"):
            stop(11, f"count.tsv counts the cutoff arm's log loss on the {r['split']} split for {r['learner']} as "
                     f"{r['better']} better and {r['worse']} worse")

    # 12 choices
    for slug in species:
        check_choices(slug, expect[slug], files[slug]["draws.tsv"], choices[slug], files[slug]["folds.tsv"])

    # 13 settings
    tables = {"page-choices.tsv": (page_choices(choices), CHOICE_COLUMNS),
              "page-trials.tsv": (page_trials(choices), TRIAL_COLUMNS),
              "page-cutoffs.tsv": (page_cutoffs(choices), CUTOFF_COLUMNS)}
    for name, (rows, columns) in tables.items():
        check_counts_only(name, rows, columns)

    # 14 binary
    bins = []
    for slug in species:
        fields, restarts = {}, 0
        for line in (reading / "score" / slug / "binary.txt").read_text(encoding="utf-8").splitlines():
            key, _, value = line.partition(" ")
            if key == "restart":
                restarts += 1
            elif key in ("commit", "grm", "code"):
                if key in fields:
                    stop(14, f"score/{slug}/binary.txt names its {key} twice")
                fields[key] = value.strip()
        if not (re.fullmatch(r"[0-9a-f]{40}", fields.get("commit", ""))
                and re.fullmatch(r"[0-9a-f]{64}", fields.get("grm", ""))
                and re.fullmatch(r"[0-9a-f]{64}", fields.get("code", ""))):
            stop(14, f"score/{slug}/binary.txt does not name a commit, a binary and a code sha256")
        bins.append({"species": by_species[slug][0]["species"], "species_slug": slug, **fields,
                     "restarts": restarts})
    if len(binary_codes(bins)) > 1 and binaries_sentence(bins).startswith("All "):
        stop(14, "the species ran more than one binary, and the page would name one")

    # 15 conditions
    read = read_conditions(counts["count.tsv"])
    for r in counts["conditions.tsv"]:
        key = (r["arm"], r["reading"], r["condition"])
        if key in read and r["met"] != str(int(read[key])):
            stop(15, f"{r['arm']}, {r['reading']}, condition {r['condition']}: conditions.tsv says {r['met']}, read "
                     f"here {int(read[key])}")
    if len([r for r in counts["conditions.tsv"] if (r["arm"], r["reading"], r["condition"]) in read]) != len(read):
        stop(15, "conditions.tsv lacks a condition read here")

    # 16 count
    rows = [r for slug in species for r in files[slug]["datasets.tsv"]]
    try:
        counted = C3.count(rows)
        beside = C3.count(rows, keep_constant=False)
        conds = C3.conditions(counted)
        plain = score.plain_count(rows)
    except SystemExit as e:
        stop(16, f"score.py refuses the count: {e}")
    want = {"count.tsv": fmt_rows(counted, C3.COUNT_COLUMNS),
            "count-no-constant.tsv": fmt_rows(beside, C3.COUNT_COLUMNS),
            "conditions.tsv": fmt_rows(conds, C3.CONDITION_COLUMNS),
            "count-plain-best.tsv": fmt_rows(plain, C3.COUNT_COLUMNS)}
    for name in COUNT_FILES:
        if counts[name] != want[name]:
            stop(16, f"count/{name} is not what score.py --count writes from the datasets.tsv files")

    return {"run": run, "species": species, "by_species": by_species, "files": files, "counts": counts,
            "choices": choices, "run_checks": run_checks, "agree": agree, "binaries": bins, "tables": tables}


def cutoff_checks(files, choices):
    """[(check, passed, made)]: the cutoff arm's log loss equals the default's in datasets.tsv; scm's cutoff arm is
    the default's in datasets.tsv, its figures and calls, where both are scored; and in each choices file scm's
    cutoff is 0.5 and its cutoff calls are the default's."""
    l_pass = l_made = s_pass = s_made = f_pass = f_made = 0
    for slug in sorted(files):
        cell = {(r["dataset"], r["split"], r["arm"], r["learner"]): r for r in files[slug]["datasets.tsv"]}
        for (d, split, arm, learner), base in sorted(cell.items()):
            if arm != DEFAULT:
                continue
            mine = cell.get((d, split, CUTOFF, learner))
            if mine is None:
                continue
            if mine["scored"] != "1" or base["scored"] != "1":
                # A missing draw leaves the dataset unscored in one arm: there is nothing to compare.
                continue
            l_made += 1
            l_pass += same_cell(mine["log_loss"], base["log_loss"])
            if learner == "scm":
                s_made += 1
                s_pass += all(same_cell(mine[c], base[c]) for c in MEASURES + ["n", "tp", "fp", "tn", "fn"])
        for key, doc in sorted(choices[slug].items()):
            for r in doc.get("learners", []):
                if r.get("learner") != "scm":
                    continue
                f_made += 1
                outer = r.get("outer") or {}
                f_pass += (r.get("cutoff") == HALF and outer.get(CUTOFF) is not None and outer.get(DEFAULT) is not None
                           and score.same_calls(outer[CUTOFF], outer[DEFAULT]))
    return [("the cutoff arm's log loss equals the default's, dataset by dataset", l_pass, l_made),
            ("scm's cutoff arm equals the default in datasets.tsv", s_pass, s_made),
            ("scm's cutoff is 0.5 and its calls the default's in each choices file", f_pass, f_made)]


def check_choices(slug, expect, draws, choices, folds_rows):
    """Gate 12 for one species. folds_rows is folds.tsv: each default draw's and study run's outer fold count, as
    run_checks.py gather read it from the default's result."""
    want = {}
    for r in folds_rows:
        key = (r["dataset"], r["scheme"], int(r["seed"]))
        if key in want or not r["folds"].isdigit():
            stop(12, f"score/{slug}/folds.tsv: {key[0]} {key[1]} seed {key[2]}: no single outer fold count")
        want[key] = int(r["folds"])
    runs = sorted({(e["dataset"], e["scheme"], int(e["seed"])) for e in expect
                   if e["arm"] == DEFAULT and e["role"] in ("draw", "study")})
    scored = {}
    for r in draws:
        if r["arm"] in CHOICE_ARMS:
            key = (r["dataset"], "lineage" if r["split"] == "lineage" else "study-whole", int(r["seed"]))
            scored.setdefault(key, set()).add((r["scored"], r["why"]))
    held = {}
    for (d, scheme, seed, fold) in choices:
        held.setdefault((d, scheme, seed), []).append(fold)
    more = sorted(set(held) - set(runs))
    if more:
        stop(12, f"score/{slug}/choices/ holds {more[0][0]} {more[0][1]} at seed {more[0][2]}, no draw or study run "
                 "of the default")
    for key in runs:
        folds = sorted(held.get(key, []))
        marks = scored.get(key, set())
        if not folds:
            if marks != {("0", MISSING_WHY)}:
                stop(12, f"score/{slug}: {key[0]} {key[1]} seed {key[2]} has no choices file and its arms are not "
                         "a missing draw")
            continue
        if folds != list(range(len(folds))) or len(folds) < 2:
            stop(12, f"score/{slug}: {key[0]} {key[1]} seed {key[2]} has choices files for folds {folds}")
        if any(s != "1" for s, _ in marks):
            # A draw with a fold missing is not scored; its files may stand short of the default's folds.
            continue
        if key not in want:
            stop(12, f"score/{slug}/folds.tsv gives {key[0]} {key[1]} seed {key[2]} no outer fold count")
        if len(folds) != want[key]:
            stop(12, f"score/{slug}: {key[0]} {key[1]} seed {key[2]} has choices files for {len(folds)} folds, its "
                     f"default's result {want[key]}")
        for fold in folds:
            doc = choices[key + (fold,)]
            check_choice_doc(slug, key, fold, doc)


def check_choice_doc(slug, key, fold, doc):
    where = f"score/{slug}: {key[0]} {key[1]} seed {key[2]} fold {fold}"
    if (doc.get("format"), doc.get("dataset"), doc.get("scheme"), doc.get("seed"), doc.get("fold")) \
            != (score.CHOICE_FORMAT, key[0], key[1], key[2], fold):
        stop(12, f"{where}: the file is not the choices file of its name")
    if (doc.get("sampler_seed"), doc.get("gamma"), doc.get("grid"), doc.get("random_trials")) \
            != (score.SAMPLER_SEED, score.GAMMA, GRID, TRIALS):
        stop(12, f"{where}: the sampler seed, gamma, grid or random trials are not the design's")
    learners = doc.get("learners") or []
    if [r.get("learner") for r in learners] != LEARNERS:
        stop(12, f"{where}: the learners are not the five in order")
    for r in learners:
        picked, plain = r.get("picked_trial"), r.get("plain_best_trial")
        if not (isinstance(picked, int) and 0 <= picked <= TRIALS):
            stop(12, f"{where} {r['learner']}: the pick {picked} is no trial id")
        if plain is not None and not (isinstance(plain, int) and 0 <= plain <= TRIALS):
            stop(12, f"{where} {r['learner']}: the plain best {plain} is no trial id")
        rel, best = r.get("release_inner_mean"), r.get("best_inner_mean")
        lead = None if rel is None or best is None else r12(best - rel)
        if (picked == 0) != (lead is None or lead <= score.GAMMA) or (picked != 0 and picked != plain):
            stop(12, f"{where} {r['learner']}: the pick {picked} does not follow the default-anchored rule")
        if r.get("cutoff") not in GRID:
            stop(12, f"{where} {r['learner']}: the cutoff {r.get('cutoff')} is no grid point")


def check_counts_only(name, rows, columns):
    """Gate 13 for one table of what the search chose."""
    allowed = {"page-choices.tsv": CHOICE_COLUMNS, "page-trials.tsv": TRIAL_COLUMNS,
               "page-cutoffs.tsv": CUTOFF_COLUMNS}[name]
    if columns != allowed:
        stop(13, f"{name} has the columns {columns}, not the counts' own")
    for r in rows:
        if set(r) != set(allowed):
            stop(13, f"{name}: a row holds the fields {sorted(set(r) - set(allowed))} beside the counts")
        for c in allowed:
            v = r[c]
            if c == "split":
                ok = v in ("lineage", "study")
            elif c == "learner":
                ok = v in LEARNERS
            elif c == "cutoff":
                ok = v in GRID
            elif c == "trial":
                ok = isinstance(v, int) and 0 <= v <= TRIALS
            else:
                ok = isinstance(v, int) and not isinstance(v, bool) and v >= 0
            if not ok:
                stop(13, f"{name}: {c} holds {v!r}, not a count, a trial id or a grid point")


def meets(row):
    """(condition 1, condition 2) of one lineage primary row of count.tsv, as score.py reads them."""
    lo, holm = num(row["interval_lo"]), num(row["holm_p"])
    c1 = lo is not None and lo > 0
    c2 = holm is not None and holm < C3.HOLM_ALPHA and int(row["better"]) > int(row["worse"])
    return c1, c2


def read_conditions(count_rows):
    """{(arm, reading, condition): met}: ruling 15's conditions 1 and 2, 3 and 4, read here off count.tsv's rows."""
    at = {(r["split"], r["pair"], r["learner"], r["measure"]): r for r in count_rows}

    def row(arm, learner, measure=PRIMARY):
        key = ("lineage", f"{arm}-{DEFAULT}", learner, measure)
        if key not in at:
            stop(15, f"count.tsv has no lineage row for {arm} with {learner} on {measure}")
        return at[key]

    def low(r):
        hi = num(r["interval_hi"])
        return hi is not None and hi < -MARGIN

    twin_most = max(sum(meets(row(TWIN, l))[0] for l in LEARNERS), sum(meets(row(TWIN, l))[1] for l in LEARNERS))
    out = {}
    for arm in ARMS:
        carry = {l for l in LEARNERS if all(meets(row(arm, l)))}
        out.update({(arm, "ruling 15", "1 and 2"): len(carry) >= C3.TOUCHED_AT_LEAST,
                    (arm, "ruling 15", "3"): twin_most < C3.TOUCHED_AT_LEAST,
                    (arm, "ruling 15", "4"): not any(low(row(arm, l, "recall")) for l in LEARNERS)})
    return out


# ---------- the page's tables ----------

def page_conditions(conds):
    """conditions.tsv's rows, each marked: ruling 15 decides (1), the stricter reading stands beside (0)."""
    return [{"arm": r["arm"], "reading": r["reading"], "decides": int(r["reading"] == "ruling 15"),
             "condition": r["condition"], "reads": r["reads"], "met": r["met"], "learners": r["learners"]}
            for r in conds]


CONDITION_COLUMNS = ["arm", "reading", "decides", "condition", "reads", "met", "learners"]


def by_construction(r):
    return arm_of(r["pair"]) == CUTOFF and r["measure"] == "log_loss"


def page_count(counts):
    """count.tsv's rows, the count that decides, each marked: in the Holm family or not, conditions 1 and 2 read on
    each family row, the cutoff arm's log loss rows equal by construction, deciding nothing, and the twin's rows
    read for condition 3 alone."""
    out = []
    for r in counts["count.tsv"]:
        equal = by_construction(r)
        row = {"decides": int(not equal), **r, "family": int(r["holm_p"] != ""), "meets_1": "", "meets_2": "",
               "note": EQUAL_BY_CONSTRUCTION if equal else TWIN_NOTE if arm_of(r["pair"]) == TWIN else ""}
        if row["family"]:
            c1, c2 = meets(r)
            row["meets_1"], row["meets_2"] = int(c1), int(c2)
        out.append(row)
    return out


COUNT_COLUMNS = ["decides"] + C3.COUNT_COLUMNS + ["family", "meets_1", "meets_2", "note"]


def page_beside(counts):
    """count-plain-best.tsv's and count-no-constant.tsv's rows, each beside the count, in no Holm family and
    deciding nothing."""
    out = []
    for beside, name in (("plain-best", "count-plain-best.tsv"), ("no-constant", "count-no-constant.tsv")):
        for r in counts[name]:
            out.append({"beside": beside, "decides": 0, **r, "family": 0,
                        "note": EQUAL_BY_CONSTRUCTION if by_construction(r) else ""})
    return out


BESIDE_COLUMNS = ["beside", "decides"] + C3.COUNT_COLUMNS + ["family", "note"]


def page_scm(counts):
    """SCM's search and plain best against the default on the primary measure, per split, beside the synthetic loss
    the plan names for the default-anchored rule: 0.0126 (0.0088 to 0.0163), the default's value less the rule's on
    the lineage split. The synthetic figure decides nothing."""
    out = []
    for name, pair in (("count.tsv", f"search-{DEFAULT}"), ("count-plain-best.tsv", f"{PLAIN}-{DEFAULT}")):
        for r in counts[name]:
            if r["pair"] != pair or r["learner"] != "scm" or r["measure"] != PRIMARY:
                continue
            loss, lo, hi = SYNTHETIC_LOSS
            out.append({"split": r["split"], "pair": r["pair"], "decides": int(name == "count.tsv"),
                        "datasets": r["datasets"], "entered": r["entered"], "better": r["better"],
                        "worse": r["worse"], "same": r["same"], "median_difference": r["median_difference"],
                        "interval_lo": r["interval_lo"], "interval_hi": r["interval_hi"], "holm_p": r["holm_p"],
                        "synthetic_loss": loss if r["split"] == "lineage" else None,
                        "synthetic_loss_lo": lo if r["split"] == "lineage" else None,
                        "synthetic_loss_hi": hi if r["split"] == "lineage" else None})
    return out


SCM_COLUMNS = ["split", "pair", "decides", "datasets", "entered", "better", "worse", "same", "median_difference",
               "interval_lo", "interval_hi", "holm_p", "synthetic_loss", "synthetic_loss_lo", "synthetic_loss_hi"]


def choice_rows(choices):
    """[(split, learner row)] over every species' choices files."""
    out = []
    for slug in sorted(choices):
        for (d, scheme, seed, fold), doc in sorted(choices[slug].items()):
            for r in doc.get("learners", []):
                out.append((split_of(scheme), r))
    return out


def page_choices(choices):
    """Per split and learner, over the outer folds searched: those that kept trial 0, those that picked another,
    the plain best's place, and the folds where trial 0 or every trial failed inside. Counts only."""
    acc = {}
    for split, r in choice_rows(choices):
        a = acc.setdefault((split, r["learner"]), {c: 0 for c in CHOICE_COLUMNS[2:]})
        plain, picked = r.get("plain_best_trial"), r.get("picked_trial")
        a["folds"] += 1
        a["kept_trial_0"] += picked == 0
        a["picked_other"] += picked != 0
        a["plain_best_trial_0"] += plain == 0
        a["plain_best_is_pick"] += plain is not None and plain == picked
        a["plain_best_differs"] += plain is not None and plain != picked
        a["plain_best_refitted"] += bool(r.get("plain_best_refitted"))
        a["trial_0_failed"] += r.get("status") == "release failed inside"
        a["every_trial_failed"] += plain is None
        a["failed_trials"] += int(r.get("failed_trials") or 0)
        a["search_served_default"] += r.get("search_note") is not None
        a["plain_served_default"] += r.get("plain_best_note") is not None
        a["cutoff_half"] += r.get("cutoff") == HALF
        a["cutoff_moved"] += r.get("cutoff") != HALF
    out = []
    for split in ("lineage", "study"):
        for learner in LEARNERS:
            if (split, learner) in acc:
                out.append({"split": split, "learner": learner,
                            **{k: int(v) for k, v in acc[(split, learner)].items()}})
    return out


def page_trials(choices):
    """Per split, learner and trial id 0 to 30: the folds the search picked it on, and the folds it was the plain
    best on. A trial is named by its id alone."""
    picked, plain, seen = {}, {}, set()
    for split, r in choice_rows(choices):
        seen.add((split, r["learner"]))
        k = (split, r["learner"], r.get("picked_trial"))
        picked[k] = picked.get(k, 0) + 1
        if r.get("plain_best_trial") is not None:
            k = (split, r["learner"], r["plain_best_trial"])
            plain[k] = plain.get(k, 0) + 1
    out = []
    for split in ("lineage", "study"):
        for learner in LEARNERS:
            if (split, learner) not in seen:
                continue
            for t in range(TRIALS + 1):
                out.append({"split": split, "learner": learner, "trial": t,
                            "picked": picked.get((split, learner, t), 0), "plain_best": plain.get((split, learner, t), 0)})
    return out


def page_cutoffs(choices):
    """Per split, learner and grid point: the folds whose cutoff it was."""
    acc, seen = {}, set()
    for split, r in choice_rows(choices):
        seen.add((split, r["learner"]))
        k = (split, r["learner"], r.get("cutoff"))
        acc[k] = acc.get(k, 0) + 1
    out = []
    for split in ("lineage", "study"):
        for learner in LEARNERS:
            if (split, learner) not in seen:
                continue
            for c in GRID:
                out.append({"split": split, "learner": learner, "cutoff": c, "folds": acc.get((split, learner, c), 0)})
    return out


def cells(rows):
    return {(r["dataset"], r["split"], r["arm"], r["learner"]): r for r in rows}


def diffs(rows, split, arm, learner, datasets, measure=PRIMARY):
    """[(dataset, difference)] of one split, pair and learner over `datasets`, the arm's value less the default's at
    12 decimals as score.py takes them; the datasets not scored; the datasets with a constant-model draw."""
    cell = cells(rows)
    out, unscored, constant = [], 0, 0
    for d in sorted(datasets):
        a, b = cell.get((d, split, arm, learner)), cell.get((d, split, DEFAULT, learner))
        if a is None:
            continue
        va, vb = num(a[measure]), (num(b[measure]) if b else None)
        if b is None or a["scored"] != "1" or b["scored"] != "1" or va is None or vb is None:
            unscored += 1
            continue
        out.append((d, r12(va - vb)))
        constant += int(a["constant_draws"] or 0) > 0
    return out, unscored, constant


def tally(values, measure=PRIMARY):
    up = -1 if measure in C3.LOWER_IS_BETTER else 1
    better = sum(1 for x in values if up * x > MARGIN)
    worse = sum(1 for x in values if up * x < -MARGIN)
    return better, worse, len(values) - better - worse


def median(values):
    return r12(statistics.median(values)) if values else None


def species_name(by_species, slug):
    return by_species[slug][0]["species"]


def page_species(found):
    """Per species, split, pair and learner, the primary measure's count over the species' datasets; it decides
    nothing."""
    out = []
    for slug in found["species"]:
        rows = found["files"][slug]["datasets.tsv"]
        pool = [r["dataset"] for r in found["by_species"][slug]]
        for split, arms in (("lineage", FAMILY + (PLAIN,)), ("study", ARMS + (PLAIN,))):
            for arm in arms:
                for learner in LEARNERS:
                    got, unscored, constant = diffs(rows, split, arm, learner, pool)
                    values = [x for _, x in got]
                    b, w, s = tally(values)
                    out.append({"species": species_name(found["by_species"], slug), "species_slug": slug,
                                "split": split, "pair": f"{arm}-{DEFAULT}", "learner": learner, "decides": 0,
                                "datasets": len(pool), "entered": len(values), "not_scored": unscored,
                                "constant": constant, "better": b, "worse": w, "same": s, "median": median(values)})
    return out


SPECIES_COLUMNS = ["species", "species_slug", "split", "pair", "learner", "decides", "datasets", "entered",
                   "not_scored", "constant", "better", "worse", "same", "median"]


def page_datasets(found):
    """Per dataset, split, pair and learner: the default's value and the arm's on balanced accuracy, recall and log
    loss, the differences, the draws scored, the genome count, the constant-model draws and a study run refused."""
    out = []
    for slug in found["species"]:
        cell = cells(found["files"][slug]["datasets.tsv"])
        for lrow in found["by_species"][slug]:
            d = lrow["dataset"]
            for split, arms in (("lineage", FAMILY + (PLAIN,)), ("study", ARMS + (PLAIN,))):
                for arm in arms:
                    for learner in LEARNERS:
                        base, mine = cell[(d, split, DEFAULT, learner)], cell[(d, split, arm, learner)]
                        row = {"dataset": d, "species": lrow["species"], "antibiotic": lrow["antibiotic"],
                               "genomes": lrow["samples"], "split": split, "pair": f"{arm}-{DEFAULT}",
                               "learner": learner, "draws": mine["draws"], "default_draws_scored": base["scored_draws"],
                               "draws_scored": mine["scored_draws"], "constant_draws": mine["constant_draws"],
                               "study_refused": int(split == "study" and int(mine["constant_draws"] or 0) > 0),
                               "note": EQUAL_BY_CONSTRUCTION if arm == CUTOFF else ""}
                        for m, names in ((PRIMARY, ("default_value", "value", "difference")),
                                         ("recall", ("default_recall", "recall_value", "recall_difference")),
                                         ("log_loss", ("default_log_loss", "log_loss_value", "log_loss_difference"))):
                            va = num(mine[m]) if mine["scored"] == "1" else None
                            vb = num(base[m]) if base["scored"] == "1" else None
                            row[names[0]], row[names[1]] = vb, va
                            row[names[2]] = r12(va - vb) if va is not None and vb is not None else None
                        out.append(row)
    return out


DATASET_COLUMNS = ["dataset", "species", "antibiotic", "genomes", "split", "pair", "learner", "draws",
                   "default_draws_scored", "draws_scored", "constant_draws", "study_refused", "default_value", "value",
                   "difference", "default_recall", "recall_value", "recall_difference", "default_log_loss",
                   "log_loss_value", "log_loss_difference", "note"]


def page_checks(found):
    """score.py's checks over every species, as checks passed of checks made per check, and the page's own checks of
    the cutoff arm."""
    by = {}
    for r in found["run_checks"]:
        # A check's name may carry its fold: "fold 3 rf: ..." is counted under "fold <k> rf: ...".
        name = re.sub(r"^fold \d+", "fold <k>", r["check"])
        n = by.setdefault(name, [0, 0])
        n[0] += r["passed"] == "1"
        n[1] += 1
    out = [{"source": "score.py checks.tsv", "check": c, "passed": p, "made": m} for c, (p, m) in sorted(by.items())]
    out += [{"source": "page_tables.py", "check": c, "passed": p, "made": m} for c, p, m in found["agree"]]
    return out


CHECK_COLUMNS = ["source", "check", "passed", "made"]
BINARY_COLUMNS = ["species", "species_slug", "commit", "grm", "code", "restarts"]


def binary_codes(bins):
    """The code sha256 values the species ran on."""
    return sorted({b["code"] for b in bins})


# ---------- writing ----------

README = """# The settings reading: the record

This directory holds the tables behind [docs/settings-arms-result.md](../../docs/settings-arms-result.md): plan block C5, a nested search of each learner's settings among 31 trials under the default-anchored rule at 0.04 (`search`), and a call cutoff chosen on inner scores from a 19-point grid (`cutoff`), each against the release's rule, beside a twin and the plain best with no anchor (`plain-best`), read once on the first run's {datasets} datasets in {species} species, each with a lineage and a study split. Every row is an aggregate over a dataset, a species, an arm, a fold count or a learner. No file holds a genome, a sequence or a feature, and no file holds a trial's settings: a trial is named by its id alone.

`scripts/c5-settings/score.py --count` wrote the four tables under `count/`. `checks.txt` is the printout of `scripts/c5-settings/run_checks.py`, the checks of the run's copies made before any score was read. `scripts/c5-settings/page_tables.py write` checked the reading against the {gates} gates of `gates.txt`, and wrote every other file here once they passed. `reading.txt` lists the reading's commands, each with the sha256 of what it wrote. The rules of the reading were fixed at 05:52Z on 23 Nov 2026, before any C5 score was read (`docs/assumption-plan-corrections.md`, entry "C5 design").

| File | Holds |
|---|---|
| `checks.txt` | the checks of the run's copies before the reading, counts only |
| `gates.txt` | the gates the table script passed before it wrote a table |
| `reading.txt` | the reading's commands and the sha256 of each one's output |
| `count/count.tsv`, `count/conditions.tsv` | the count that decides, and the plan's conditions per arm as `score.py` read them |
| `count/count-plain-best.tsv`, `count/count-no-constant.tsv` | the counts beside it, each in no Holm family and deciding nothing |
| `page-conditions.tsv` | the conditions per arm, `decides` 1 for ruling 15 and 0 for the stricter reading beside it; condition 6 is the maintainer's |
| `page-count.tsv` | every row of the count that decides, with the 15 rows of the Holm family marked and conditions 1 and 2 read on each; the cutoff arm's log loss rows marked equal by construction, deciding nothing |
| `page-beside.tsv` | every row of the two counts beside it, each marked as deciding nothing |
| `page-scm.tsv` | `scm`'s search and plain best against the default on balanced accuracy, beside the synthetic loss of 0.0126 (0.0088 to 0.0163) the plan names for the default-anchored rule |
| `page-choices.tsv` | per split and learner, the outer folds searched: those that kept trial 0, the plain best's place, the folds where trial 0 or every trial failed inside, and the cutoffs moved from 0.5 |
| `page-trials.tsv` | per split, learner and trial id, the folds the search picked it on and the folds it was the plain best on; ids alone |
| `page-cutoffs.tsv` | per split, learner and grid point, the folds whose cutoff it was |
| `page-species.tsv` | per species, named in full, split, pair and learner: better, worse and the same on the primary measure; it decides nothing |
| `page-datasets.tsv` | per dataset, species and antibiotic in full, split, pair and learner: the default's value and the arm's on balanced accuracy, recall and log loss, the differences, the draws scored, the genome count, the constant-model draws and a study run refused |
| `page-checks.tsv` | `score.py`'s checks of every result as checks passed of checks made, and the table script's own checks of the cutoff arm |
| `page-binaries.tsv` | each species' commit, binary and code sha256, and its restarts |

{binaries}

A dataset's figure is the mean of its draws: seeds 42, 43 and 44 on the lineage split, the twin's 45, 46 and 47, and one study run at the release's seeds. The cutoff arm moves no score, so its log loss equals the default's on every dataset by construction, and its log loss rows decide nothing. `scm` scores 0 or 1, so its cutoff arm is the default on every fold. A split by lineage or study is a retrospective split inside one public database, never external validation.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/c5-settings/page_tables.py print` prints the page's tables from these files.
"""


def binaries_sentence(bins):
    codes = binary_codes(bins)
    commits = sorted({b["commit"] for b in bins})
    if len(codes) == 1:
        return (f"All {len(bins)} species ran one binary by its code sha256, `{codes[0]}`, built from "
                f"{len(commits)} commit{'s' if len(commits) != 1 else ''}: "
                + ", ".join(f"`{c[:7]}`" for c in commits) + ".")
    return (f"The {len(bins)} species ran {len(codes)} binaries by their code sha256; `page-binaries.tsv` names each "
            "binary with its species.")


def write_tsv(path, rows, columns):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for r in rows:
            handle.write("\t".join(fmt(r.get(c)) for c in columns) + "\n")


def write(reading, l_path, checks_path, out):
    reading, out = Path(reading), Path(out)
    lrows = load(l_path)
    if out.exists() and any(out.iterdir()):
        raise Stop(f"{out} is not empty")
    found = check_reading(reading, lrows, checks_path)
    (out / "count").mkdir(parents=True, exist_ok=True)
    for name in COUNT_FILES:
        (out / "count" / name).write_bytes((reading / "count" / name).read_bytes())
    (out / "reading.txt").write_bytes((reading / "reading.txt").read_bytes())
    (out / "checks.txt").write_bytes(Path(checks_path).read_bytes())
    write_tsv(out / "page-conditions.tsv", page_conditions(found["counts"]["conditions.tsv"]), CONDITION_COLUMNS)
    write_tsv(out / "page-count.tsv", page_count(found["counts"]), COUNT_COLUMNS)
    write_tsv(out / "page-beside.tsv", page_beside(found["counts"]), BESIDE_COLUMNS)
    write_tsv(out / "page-scm.tsv", page_scm(found["counts"]), SCM_COLUMNS)
    for name, (rows, columns) in found["tables"].items():
        write_tsv(out / name, rows, columns)
    write_tsv(out / "page-species.tsv", page_species(found), SPECIES_COLUMNS)
    write_tsv(out / "page-datasets.tsv", page_datasets(found), DATASET_COLUMNS)
    write_tsv(out / "page-checks.tsv", page_checks(found), CHECK_COLUMNS)
    write_tsv(out / "page-binaries.tsv", found["binaries"], BINARY_COLUMNS)
    with open(out / "gates.txt", "w", encoding="utf-8") as handle:
        for n, name, words in GATES:
            handle.write(f"{n} {name}: {words}\n")
    (out / "README.md").write_text(README.format(
        datasets=len(found["run"]), species=len(found["species"]), gates=len(GATES),
        binaries=binaries_sentence(found["binaries"])), encoding="utf-8")
    files = sorted(p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file())
    with open(out / "SHA256SUMS", "w", encoding="utf-8") as handle:
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


def d4(x):
    return "none" if x in ("", None) else str(Decimal(str(x)).quantize(Decimal("0.0001"), rounding=ROUND_HALF_EVEN))


def pval(x):
    if x in ("", None):
        return "none"
    return "< 0.001" if float(x) < 0.001 else d3(x)


def count_words(x):
    return f"{int(x):,}"


def of(part, whole_):
    return f"{count_words(part)} of {count_words(whole_)}"


def interval(lo, hi):
    return "none" if lo in ("", None) else f"{d3(lo)} to {d3(hi)}"


def yes(x):
    return {"1": "met", "0": "not met", "": "the maintainer's"}[str(x)]


def table(head, rows):
    lines = ["| " + " | ".join(head) + " |", "|" + "|".join("---" for _ in head) + "|"]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return lines


def pair_words(pair):
    return f"{ARM_WORDS[arm_of(pair)]} less default"


def print_conditions(rows):
    return table(["Arm", "Reading", "Decides", "Condition", "Reads", "Met", "Learners"],
                 [[ARM_WORDS[r["arm"]], r["reading"], "yes" if r["decides"] == "1" else "no, beside", r["condition"],
                   r["reads"], yes(r["met"]), r["learners"] or "none"] for r in rows])


def print_family(rows):
    out = []
    for r in rows:
        if r["family"] != "1":
            continue
        out.append([pair_words(r["pair"]), f"`{r['learner']}`", of(r["better"], r["entered"]), r["worse"], r["same"],
                    d3(r["median_difference"]), interval(r["interval_lo"], r["interval_hi"]), pval(r["sign_p"]),
                    pval(r["holm_p"]), yes(r["meets_1"]), yes(r["meets_2"])])
    return table(["Pair", "Learner", "Better", "Worse", "Same", "Median", "Interval", "Sign p", "Holm p",
                  "Condition 1", "Condition 2"], out)


def print_count(rows, measures=(PRIMARY, "recall", "log_loss")):
    out = []
    for r in rows:
        if r["measure"] not in measures:
            continue
        out.append([r["split"], pair_words(r["pair"]), f"`{r['learner']}`", MEASURE_WORDS[r["measure"]],
                    r["orientation"], of(r["entered"], r["datasets"]), r["better"], r["worse"], r["same"],
                    r["constant"], d3(r["median_difference"]), interval(r["interval_lo"], r["interval_hi"]),
                    "no, equal by construction" if r["decides"] == "0" else "yes"])
    return table(["Split", "Pair", "Learner", "Measure", "Better is", "Entered", "Better", "Worse", "Same",
                  "Constant model", "Median", "Interval", "Decides"], out)


def print_measures(rows):
    """Per split, pair and learner: the median difference of each measure, with better and worse, recall beside log
    loss."""
    at = {(r["split"], r["pair"], r["learner"], r["measure"]): r for r in rows}
    keys = []
    for r in rows:
        k = (r["split"], r["pair"], r["learner"])
        if k not in keys:
            keys.append(k)
    out = []
    for k in keys:
        cells_ = [k[0], pair_words(k[1]), f"`{k[2]}`"]
        for m in MEASURES:
            r = at[k + (m,)]
            cells_.append(f"{d3(r['median_difference'])} ({r['better']} / {r['worse']})")
        out.append(cells_)
    return table(["Split", "Pair", "Learner"] + [MEASURE_WORDS[m].capitalize() for m in MEASURES], out)


def print_beside(rows):
    out = []
    for r in rows:
        if r["measure"] not in (PRIMARY, "recall", "log_loss"):
            continue
        out.append([r["beside"], r["split"], pair_words(r["pair"]), f"`{r['learner']}`", MEASURE_WORDS[r["measure"]],
                    of(r["entered"], r["datasets"]), r["better"], r["worse"], r["same"], d3(r["median_difference"]),
                    interval(r["interval_lo"], r["interval_hi"]), pval(r["sign_p"])])
    return table(["Beside", "Split", "Pair", "Learner", "Measure", "Entered", "Better", "Worse", "Same", "Median",
                  "Interval", "Sign p"], out)


def print_scm(rows):
    out = []
    for r in rows:
        synthetic = ("none" if r["synthetic_loss"] == "" else
                     f"{d4(r['synthetic_loss'])} ({d4(r['synthetic_loss_lo'])} to {d4(r['synthetic_loss_hi'])})")
        out.append([r["split"], pair_words(r["pair"]), of(r["entered"], r["datasets"]), r["better"], r["worse"],
                    r["same"], d3(r["median_difference"]), interval(r["interval_lo"], r["interval_hi"]),
                    pval(r["holm_p"]), synthetic])
    return table(["Split", "Pair", "Entered", "Better", "Worse", "Same", "Median", "Interval", "Holm p",
                  "Synthetic loss of the rule, default less rule"], out)


def print_choices(rows):
    return table(["Split", "Learner", "Folds", "Kept trial 0", "Picked another", "Plain best trial 0",
                  "Plain best is the pick", "Plain best differs", "Trial 0 failed", "Every trial failed",
                  "Cutoff moved from 0.5"],
                 [[r["split"], f"`{r['learner']}`", count_words(r["folds"]), of(r["kept_trial_0"], r["folds"]),
                   count_words(r["picked_other"]), count_words(r["plain_best_trial_0"]),
                   count_words(r["plain_best_is_pick"]), of(r["plain_best_differs"], r["folds"]),
                   count_words(r["trial_0_failed"]), count_words(r["every_trial_failed"]),
                   of(r["cutoff_moved"], r["folds"])] for r in rows])


def print_cutoffs(rows):
    keys, at = [], {}
    for r in rows:
        k = (r["split"], r["learner"])
        if k not in keys:
            keys.append(k)
        at[k + (r["cutoff"],)] = r["folds"]
    grid = [g for g in sorted({r["cutoff"] for r in rows}, key=float)]
    return table(["Split", "Learner"] + grid, [[k[0], f"`{k[1]}`"] + [at.get(k + (g,), "0") for g in grid]
                                               for k in keys])


def print_checks(rows):
    return table(["Source", "Check", "Passed"], [[r["source"], r["check"], of(r["passed"], r["made"])] for r in rows])


def print_species(rows):
    """Per species, split and pair: better, worse and the same summed over the five learners, as dataset and learner
    pairs of the denominator."""
    groups = {}
    for r in rows:
        groups.setdefault((r["species"], r["split"], r["pair"]), []).append(r)
    out = []
    for (sp, split, pair), rs in groups.items():
        pairs = sum(int(r["datasets"]) for r in rs)
        b, w, s = (sum(int(r[k]) for r in rs) for k in ("better", "worse", "same"))
        out.append([f"*{sp}*", split, pair_words(pair), rs[0]["datasets"], f"{b} / {w} / {s} of {pairs}",
                    sum(int(r["not_scored"]) for r in rs)])
    return table(["Species", "Split", "Pair", "Datasets", "Better / worse / same, dataset and learner pairs",
                  "Not scored"], out)


def print_datasets(rows):
    groups = {}
    for r in rows:
        groups.setdefault((r["species"], r["antibiotic"], r["genomes"], r["split"]), {}).setdefault(
            r["pair"], []).append(r)
    out = []
    for (sp, ab, genomes, split), pairs in groups.items():
        cells_ = [f"*{sp}*", ab, count_words(genomes), split]
        for arm in FAMILY + (PLAIN,):
            rs = pairs.get(f"{arm}-{DEFAULT}")
            cells_.append(" / ".join(d3(r["difference"]) for r in rs) if rs else "")
        refused = any(r["study_refused"] == "1" for rs in pairs.values() for r in rs)
        cells_.append("yes" if refused else "no")
        out.append(cells_)
    return table(["Species", "Antibiotic", "Genomes", "Split", "Search: five learners", "Cutoff: five learners",
                  "Twin: five learners", "Plain best: five learners", "Study run refused"], out)


def print_binaries(rows):
    codes = binary_codes(rows)
    if len(codes) < 2:
        return [f"All {len(rows)} species ran one binary by its code sha256, `{codes[0][:16]}`."]
    return table(["Species", "Commit", "Code sha256", "Restarts"],
                 [[f"*{r['species']}*", f"`{r['commit'][:7]}`", f"`{r['code'][:16]}`", r["restarts"]] for r in rows])


def tables(folder):
    folder = Path(folder)
    count = load(folder / "page-count.tsv")
    return {
        "conditions": print_conditions(load(folder / "page-conditions.tsv")),
        "family": print_family(count),
        "count": print_count(count),
        "measures": print_measures(count),
        "beside": print_beside(load(folder / "page-beside.tsv")),
        "scm": print_scm(load(folder / "page-scm.tsv")),
        "choices": print_choices(load(folder / "page-choices.tsv")),
        "cutoffs": print_cutoffs(load(folder / "page-cutoffs.tsv")),
        "checks": print_checks(load(folder / "page-checks.tsv")),
        "binaries": print_binaries(load(folder / "page-binaries.tsv")),
        "species": print_species(load(folder / "page-species.tsv")),
        "datasets": print_datasets(load(folder / "page-datasets.tsv")),
    }


def main(argv):
    if len(argv) == 5 and argv[0] == "write":
        try:
            files = write(argv[1], argv[2], argv[3], argv[4])
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
