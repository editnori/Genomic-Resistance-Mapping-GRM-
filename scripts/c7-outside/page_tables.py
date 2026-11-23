"""The tables of docs/outside-learners-result.md: plan block C7's one reading, as files and as the page prints them.

  python3 -I page_tables.py write <reading dir> <lists.tsv> <out dir>
  python3 -I page_tables.py print <out dir>

The reading directory holds what run_checks.py gather copied and what the
count wrote, and nothing else:

  reading.txt          the reading's commands, each with the sha256 of what it wrote
  checks.txt           run_checks.py check's printout, every check passed
  withdrawn.tsv        arm, reason, time: each library arm the lead withdrew before
                       it ran, with the reason and the time the lead gave (a header
                       alone when none was)
  count/               score.py --count's count.tsv, count-no-constant.tsv,
                       count-no-unconverged.tsv and conditions.tsv
  score/<species>/     each species' score/ as score.py wrote it: checks.tsv,
                       draws.tsv and datasets.tsv
  run/<species>/       from each species' copy: DONE, libraries.txt, timing.tsv
                       (logs/timing.tsv), plans.log (logs/plans.log) and outside/,
                       the run's results/outside/ (one grm-c7-outside-1 report per
                       library arm, dataset, split and seed: counts and calls)

<lists.tsv> is scripts/c7-outside/lists.tsv: the first run is its rows with
first_run 1, L's rows not in M; the study split runs on those in S.

`write` passes every gate of GATES before it writes a file, and stops at the
first that fails, naming it. The four count files are score.py's own
tables: `write` counts them again from the datasets.tsv files (gate 17),
copies them under count/ as they are, and sums them in SHA256SUMS with every
other file. score.py's checks.tsv per species is summed by check name into
page-checks.tsv, with plans.py's checks from its count line. It then writes
the page's tables as page-*.tsv, reading.txt, checks.txt, gates.txt,
README.md and SHA256SUMS into <out dir> (absent or empty). Every row is an
aggregate over a dataset, a species, a split, an arm, a pair or a fold
count: no file holds a genome, a k-mer, a unitig, a sequence, a position, a
feature name, a weight or an importance.

`print` prints the page's tables from <out dir> alone, figures rounded to 3
decimals. Standard library only, and numpy where score.py needs it; it
imports none of XGBoost, LightGBM or scikit-learn.
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


score = _load("c7_outside_score_for_page", HERE / "score.py")
O, C9 = score.O, score.C9

DEFAULT, TWIN, ARM_BASE = score.DEFAULT, score.TWIN, score.ARM_BASE
ARMS = ("xgboost", "lightgbm", "elasticnet")
LEARNER_OF = dict(score.ARMS)
PAIRED = tuple(score.PAIRED)
FAMILY = [tuple(p) for p in score.FAMILY]
MEASURES, PRIMARY, MARGIN = score.MEASURES, score.PRIMARY, score.MARGIN
SEEDS, TWIN_SEEDS, STUDY_SEED = (42, 43, 44), (45, 46, 47), 42
ARM_SEEDS = SEEDS
SPLITS = ("lineage", "study")
STUDY_SCHEME = "study-whole"
SCORE_FILES = ("checks.tsv", "draws.tsv", "datasets.tsv")
COUNT_FILES = ("count.tsv", "count-no-constant.tsv", "count-no-unconverged.tsv", "conditions.tsv")
DECIDING_COUNT = "count.tsv"
RUN_FILES = ("DONE", "libraries.txt", "timing.tsv", "plans.log")
LIBRARIES = ("python", "numpy", "scipy", "sklearn", "xgboost", "lightgbm")
DONE_LINE = re.compile(r"DONE (?P<species>\S+) (?P<time>\S+) grm (?P<grm>[0-9a-f]{64}) code (?P<code>[0-9a-f]{64}) "
                       r"export-grm (?P<export>[0-9a-f]{64}) elasticnet (?P<elasticnet>run|leave-out) "
                       r"results (?P<results>\d+) SHA256SUMS (?P<sums>[0-9a-f]{64})")
# A species' first copy, moved aside by the driver when the species restarts (`<slug>.<%Y%m%dT%H%M%SZ>`).
ASIDE = re.compile(r"(?P<slug>[a-z0-9_]+)\.(?P<time>\d{8}T\d{6}Z)")
PLANS_LINE = re.compile(r"plans: (?P<passed>\d+) of (?P<made>\d+) checks pass, (?P<seeds>\d+) dataset seeds")
CHECKS_LINE = re.compile(r"CHECKS (?P<passed>\d+) of (?P<made>\d+)")
# run_checks.py makes six checks: species, SHA256SUMS, DONE, expect.py, genome ids and features.
RUN_CHECKS = ("species", "SHA256SUMS", "DONE", "expect.py", "genome ids", "features")
TIME = re.compile(r"\d{4}-\d\d-\d\dT\d\d:\d\d(:\d\d)?Z")
# run_checks.py's species line when the elastic net is withdrawn: how many species ran a second pass, set aside.
SECOND_PASSES = re.compile(r"; the elastic net withdrawn: (?P<n>\d+) species ran a second pass, \d+ second-pass copies "
                           r"set aside")
# plans.py's checks when every one passes: three at each lineage dataset seed of the arms, and one more at seed 42.
PLANS_EXPORT = "the export at seed {seed} is read"
PLANS_CHECKS = ("the plan's folds by class are the default's", "the plan's folds by genome are the default's")
PLANS_AT_42 = "the input fingerprint is the release report's"
GATES = [
    (1, "reading", "the reading holds reading.txt, checks.txt with every check passed, withdrawn.tsv, count/ with the "
                   "four files, score/ and run/, and nothing else; withdrawn.tsv names, with a reason and a time, each "
                   "library arm count.tsv marks withdrawn, and no other"),
    (2, "species", "score/ and run/ each hold every species of lists.tsv's first run once, with its files; a folder "
                   "named <species slug>.<UTC time> is a first copy set aside and never read as a species"),
    (3, "M", "no dataset of M, or outside the first run, stands in any table; species and antibiotic are lists.tsv's"),
    (4, "pairs", "xgboost and lightgbm stand on gbt only and elasticnet on lasso only, each read against the "
                 "default's own learner at the arm's seeds (k31-at-arm-seeds); count.tsv holds the family's pairs"),
    (5, "twin", "the twin's rows are lineage draws at seeds 45, 46 and 47 only, no arm's row is at a twin seed, and "
                "conditions.tsv reads the three arms and never the twin"),
    (6, "study", "a study row is a study run at seed 42 on the study split, and every lineage row is a draw at 42, 43 "
                 "or 44 or the twin's: a study row is never counted on the lineage split"),
    (7, "family", "each count holds the five lineage pairs of the family on the primary measure, withdrawn or not, "
                  "each with its Holm p"),
    (8, "withdrawn", "a withdrawn arm holds no row, no report and no fold, enters no dataset of its pair, and its pair "
                     "stays in the family at p = 1"),
    (9, "passes", "a species whose DONE says the elastic net was left out (its first pass alone) holds no elastic net "
                  "row or report, and unless the elastic net is withdrawn every species' DONE says it ran"),
    (10, "rows", "each species' datasets.tsv and draws.tsv hold exactly the rows of its first-run datasets for the "
                 "default, the twin and each arm not withdrawn, the study split on those in S"),
    (11, "missing", "a draw not scored holds no figure, and a dataset with a draw not scored is not scored: a "
                    "missing draw is never read as zero"),
    (12, "draws", "datasets.tsv is the mean of draws.tsv's draws, k31-at-arm-seeds the default's at the arms' seeds"),
    (13, "reports", "run/<species>/outside/ holds one report per library arm draw of draws.tsv and no other, each of "
                    "its listed keys only, its fold counts the sums over its folds, draws.tsv's marks its own"),
    (14, "checks", "each species' checks.tsv holds checks and every one passed, and plans.log's count line says each "
                   "plans check passed at every lineage dataset seed of the arms"),
    (15, "binary", "each DONE names its species and its result count; one code sha256 by binary.txt's code line (the "
                   "DONE's code) and one version of each library in every species, or the page names each with its "
                   "species"),
    (16, "beside", "count-no-constant.tsv and count-no-unconverged.tsv are written as deciding nothing, never as the "
                   "count, and the arm against its own learner's twin pair alone decides nothing"),
    (17, "count", "the four count files are what score.py --count writes from the datasets.tsv files, with the arms "
                  "withdrawn.tsv names"),
]
MEASURE_WORDS = {"balanced_accuracy": "balanced accuracy", "recall": "recall", "specificity": "specificity",
                 "vme": "very major error rate", "me": "major error rate", "log_loss": "log loss"}
ARM_WORDS = {"xgboost": "XGBoost 3.0.5", "lightgbm": "LightGBM 4.6.0", "elasticnet": "scikit-learn 1.7.2 elastic net",
             "twin": "the twin, the default at seeds 45, 46 and 47"}


class Stop(Exception):
    pass


# ---------- reading files ----------

def load(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def num(x):
    return None if x in ("", None) else float(x)


def whole(x):
    return None if x in ("", None) else int(x)


def r12(x):
    return None if x is None else round(x, 12)


def pair_arm(pair):
    """The arm of a count.tsv pair, `<arm>-k31`."""
    return pair[: -len(f"-{DEFAULT}")]


def pair_name(arm, learner):
    """The pair as the page names it: xgboost-gbt, lightgbm-gbt, elasticnet-lasso, twin-gbt, twin-lasso."""
    return f"{arm}-{learner}"


def against(arm):
    return DEFAULT if arm == TWIN else ARM_BASE


def scheme_of(split):
    return "lineage" if split == "lineage" else STUDY_SCHEME


def same_cell(a, b):
    if a == b:
        return True
    try:
        return abs(float(a) - float(b)) <= 2e-12
    except ValueError:
        return False


class Reading:
    """What `write` reads: the count files, each species' score/ tables and run files, and the lists."""

    def __init__(self, folder, lrows):
        self.folder = Path(folder)
        self.lrows = lrows
        self.run = [r for r in lrows if r["first_run"] == "1"]
        self.species = sorted({r["species_slug"] for r in self.run})
        self.by_species = {s: [r for r in self.run if r["species_slug"] == s] for s in self.species}
        self.names = {r["dataset"]: (r["species"], r["antibiotic"]) for r in lrows}
        self.genomes = {r["dataset"]: int(r["samples"]) for r in lrows}
        self.species_name = {r["species_slug"]: r["species"] for r in lrows}
        self.count, self.score, self.done, self.libraries, self.timing, self.plans = {}, {}, {}, {}, {}, {}
        self.reports = {}
        self.withdrawn = {}
        self.second_passes = 0
        self.checks_txt = ""


# ---------- gates ----------

def expected_rows(rows, withdrawn):
    """(dataset, split, arm, learner) of datasets.tsv for the first-run rows of one species."""
    out = set()
    arms = [a for a in ARMS if a not in withdrawn]
    for r in rows:
        study = r["in_S"] == "1"
        for learner in PAIRED:
            out.add((r["dataset"], "lineage", DEFAULT, learner))
            out.add((r["dataset"], "lineage", TWIN, learner))
            if study:
                out.add((r["dataset"], "study", DEFAULT, learner))
        for arm in arms:
            learner = LEARNER_OF[arm]
            for split in ("lineage", "study") if study else ("lineage",):
                out.add((r["dataset"], split, arm, learner))
        if arms:
            # score.dataset_rows reads the default again at every split and seed an arm ran, on both its learners.
            for split in ("lineage", "study") if study else ("lineage",):
                for learner in PAIRED:
                    out.add((r["dataset"], split, ARM_BASE, learner))
    return out


def expected_draws(rows, withdrawn):
    """(dataset, split, seed, role, arm, learner) of draws.tsv for the first-run rows of one species."""
    out = set()
    for d, split, arm, learner in expected_rows(rows, withdrawn):
        if arm == ARM_BASE:
            continue
        if split == "study":
            out.add((d, split, STUDY_SEED, "study", arm, learner))
        elif arm == TWIN:
            out.update((d, split, s, "twin", arm, learner) for s in TWIN_SEEDS)
        else:
            out.update((d, split, s, "draw", arm, learner) for s in (ARM_SEEDS if arm in ARMS else SEEDS))
    return out


def typed_draw(r):
    """A draws.tsv row as score.dataset_rows takes it."""
    out = dict(r)
    out["seed"] = int(r["seed"])
    out["scored"] = int(r["scored"])
    out["constant"] = int(r["constant"] or 0)
    out["unconverged"] = int(r["unconverged"] or 0)
    for m in MEASURES:
        out[m] = num(r.get(m))
    for k in ("n", "tp", "fp", "tn", "fn"):
        out[k] = whole(r.get(k)) or 0
    return out


def report_name(arm, dataset, split, seed):
    return f"{arm}/{dataset}.{scheme_of(split)}-s{seed}.outside.json"


def fold_sums(doc):
    """The report's fold counts as outside.py sums them over its folds."""
    folds = doc.get("folds", [])
    return {"constant_folds": sum(int(f.get("constant", 0)) for f in folds),
            "one_tree_folds": sum(1 for f in folds if f.get("rounds_kept") == 1),
            "null_fit_folds": sum(int(f.get("null_fit", 0)) for f in folds),
            "unconverged_refits": sum(1 for f in folds if f.get("refit_converged") == 0)}


def read_libraries(text):
    out = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 2:
            out[parts[0]] = parts[1]
    return out


def check_reading(folder, lrows):
    """Gates 1 to 15, in order, on what the reading holds. Returns the Reading; raises Stop naming the gate."""
    data = Reading(folder, lrows)
    top = data.folder

    # 1 reading
    for name in ("reading.txt", "checks.txt", "withdrawn.tsv"):
        if not (top / name).is_file():
            raise Stop(f"gate 1 reading: no {name}")
    for name in COUNT_FILES:
        if not (top / "count" / name).is_file():
            raise Stop(f"gate 1 reading: count/{name} is missing")
    for name in ("score", "run"):
        if not (top / name).is_dir():
            raise Stop(f"gate 1 reading: {name}/ is missing")
    extra = sorted(p.name for p in top.iterdir() if p.is_dir() and p.name not in ("count", "score", "run"))
    if extra:
        raise Stop(f"gate 1 reading: folders that are no part of the reading: {', '.join(extra)}")
    data.checks_txt = (top / "checks.txt").read_text(encoding="utf-8")
    lines = [x for x in data.checks_txt.splitlines() if x.strip()]
    last = CHECKS_LINE.fullmatch(lines[-1]) if lines else None
    if last is None or last["passed"] != last["made"] or any(x.startswith("FAIL") for x in lines):
        raise Stop("gate 1 reading: checks.txt does not say that every check of the run passed")
    passed = [x.split(":")[0][len("PASS "):] for x in lines if x.startswith("PASS ")]
    if int(last["made"]) != len(RUN_CHECKS) or passed != list(RUN_CHECKS):
        raise Stop(f"gate 1 reading: checks.txt does not hold run_checks.py's {len(RUN_CHECKS)} checks, each passed "
                   "once, in order")
    data.count = {name: load(top / "count" / name) for name in COUNT_FILES}
    rows = load(top / "withdrawn.tsv")
    if rows and list(rows[0]) != ["arm", "reason", "time"]:
        raise Stop("gate 1 reading: withdrawn.tsv's columns are not arm, reason and time")
    for r in rows:
        if r["arm"] not in ARMS or not (r.get("reason") or "").strip() or not TIME.fullmatch(r.get("time") or ""):
            raise Stop(f"gate 1 reading: withdrawn.tsv: {r['arm']} is no library arm or has no reason or UTC time")
        if r["arm"] in data.withdrawn:
            raise Stop(f"gate 1 reading: withdrawn.tsv names {r['arm']} twice")
        data.withdrawn[r["arm"]] = (r["reason"].strip(), r["time"])
    marked = {pair_arm(r["pair"]) for r in data.count["count.tsv"] if r["withdrawn"] == "1"}
    if marked != set(data.withdrawn):
        raise Stop(f"gate 1 reading: count.tsv marks withdrawn {sorted(marked)}, withdrawn.tsv names "
                   f"{sorted(data.withdrawn)}")
    found = [m for x in lines if x.startswith("PASS species:") for m in SECOND_PASSES.finditer(x)]
    if len(found) != int("elasticnet" in data.withdrawn):
        raise Stop("gate 1 reading: checks.txt's species line does not give, once, how many species ran a second "
                   "pass when the elastic net is withdrawn, and only then")
    data.second_passes = int(found[0]["n"]) if found else 0

    # 2 species
    species = data.species
    for name, need in (("score", SCORE_FILES), ("run", RUN_FILES)):
        got = sorted(p.name for p in (top / name).iterdir() if p.is_dir())
        for folder_name in got:
            m = ASIDE.fullmatch(folder_name)
            if m is not None:
                raise Stop(f"gate 2 species: {name}/{folder_name} is a first copy of {m['slug']} set aside at "
                           f"{m['time']}, never read as a species")
        if got != species:
            missing, more = sorted(set(species) - set(got)), sorted(set(got) - set(species))
            raise Stop(f"gate 2 species: {name}/ has {len(got)} species of {len(species)}; missing {missing}, more {more}")
        for slug in got:
            for f in need:
                if not (top / name / slug / f).is_file():
                    raise Stop(f"gate 2 species: {name}/{slug}/{f} is missing")
    for slug in species:
        data.score[slug] = {f: load(top / "score" / slug / f) for f in SCORE_FILES}
        for f in ("draws.tsv", "datasets.tsv"):
            other = sorted({r["dataset"].split("__")[0] for r in data.score[slug][f]} - {slug})
            if other or any(r.get("species_slug", slug) != slug for r in data.score[slug][f]):
                raise Stop(f"gate 2 species: score/{slug}/{f} holds rows of another species: {other}")

    # 3 M
    first = {r["dataset"] for r in data.run}
    in_m = {r["dataset"] for r in lrows if r["in_M"] == "1"}
    for slug in species:
        for f in ("draws.tsv", "datasets.tsv"):
            ds = {r["dataset"] for r in data.score[slug][f]}
            if ds & in_m:
                raise Stop(f"gate 3 M: score/{slug}/{f} holds {len(ds & in_m)} datasets of M")
            if ds - first:
                raise Stop(f"gate 3 M: score/{slug}/{f} holds {len(ds - first)} datasets outside the first run")
            for r in data.score[slug][f]:
                if (r["species"], r["antibiotic"]) != data.names[r["dataset"]]:
                    raise Stop(f"gate 3 M: score/{slug}/{f}: {r['dataset']} is named otherwise than in lists.tsv")

    # 4 pairs
    for slug in species:
        for f in ("draws.tsv", "datasets.tsv"):
            for r in data.score[slug][f]:
                arm, learner = r["arm"], r["learner"]
                if arm in ARMS:
                    if learner != LEARNER_OF[arm]:
                        raise Stop(f"gate 4 pairs: score/{slug}/{f}: {arm} on {r['dataset']} stands on {learner}, "
                                   f"not {LEARNER_OF[arm]}")
                elif arm in (DEFAULT, TWIN) or (arm == ARM_BASE and f == "datasets.tsv"):
                    if learner not in PAIRED:
                        raise Stop(f"gate 4 pairs: score/{slug}/{f}: {arm} stands on {learner}")
                else:
                    raise Stop(f"gate 4 pairs: score/{slug}/{f}: a row of {arm}, which no pair reads")
                if arm in ARMS and f == "draws.tsv":
                    seeds = ARM_SEEDS if r["split"] == "lineage" else (STUDY_SEED,)
                    if int(r["seed"]) not in seeds:
                        raise Stop(f"gate 4 pairs: score/{slug}/draws.tsv: {arm} on {r['dataset']} at seed "
                                   f"{r['seed']}, not an arm seed")
        cells = {(r["dataset"], r["split"], r["arm"], r["learner"]) for r in data.score[slug]["datasets.tsv"]}
        for d, split, arm, learner in sorted(cells):
            if arm in ARMS and (d, split, ARM_BASE, learner) not in cells:
                raise Stop(f"gate 4 pairs: score/{slug}/datasets.tsv: {arm} on {d} {split} has no {ARM_BASE} row of "
                           f"{learner} to be read against")
    for name in ("count.tsv", "count-no-constant.tsv", "count-no-unconverged.tsv"):
        got = {(r["split"], r["pair"], r["learner"]) for r in data.count[name]}
        want = {("lineage", f"{a}-{DEFAULT}", l) for a, l in FAMILY} | {
            ("study", f"{a}-{DEFAULT}", l) for a, l in FAMILY if a != TWIN}
        if got != want:
            raise Stop(f"gate 4 pairs: count/{name} holds other pairs than the family's: "
                       f"{sorted(got ^ want)}")

    # 5 twin
    for slug in species:
        for r in data.score[slug]["draws.tsv"]:
            seed, is_twin = int(r["seed"]), r["arm"] == TWIN
            if is_twin and not (r["split"] == "lineage" and r["role"] == "twin" and seed in TWIN_SEEDS):
                raise Stop(f"gate 5 twin: score/{slug}/draws.tsv: a twin row of {r['dataset']} is no lineage draw at "
                           "a twin seed")
            if not is_twin and (r["role"] == "twin" or seed in TWIN_SEEDS):
                raise Stop(f"gate 5 twin: score/{slug}/draws.tsv: {r['arm']} has a row of {r['dataset']} at a twin "
                           "seed or role")
        if any(r["arm"] == TWIN and r["split"] != "lineage" for r in data.score[slug]["datasets.tsv"]):
            raise Stop(f"gate 5 twin: score/{slug}/datasets.tsv has a twin row off the lineage split")
    cond_arms = [r["arm"] for r in data.count["conditions.tsv"]]
    if TWIN in cond_arms or sorted(set(cond_arms)) != sorted(ARMS):
        raise Stop("gate 5 twin: conditions.tsv does not read the three arms alone")

    # 6 study
    for slug in species:
        for r in data.score[slug]["draws.tsv"]:
            seed = int(r["seed"])
            if (r["role"] == "study") != (r["split"] == "study"):
                raise Stop(f"gate 6 study: score/{slug}/draws.tsv: a {r['role']} row of {r['dataset']} stands on "
                           f"the {r['split']} split")
            if r["split"] == "study" and seed != STUDY_SEED:
                raise Stop(f"gate 6 study: score/{slug}/draws.tsv: a study run of {r['dataset']} at seed {seed}")
            if r["split"] == "lineage" and r["arm"] != TWIN and not (r["role"] == "draw" and seed in SEEDS):
                raise Stop(f"gate 6 study: score/{slug}/draws.tsv: a lineage row of {r['dataset']} is no draw")
            if r["split"] not in SPLITS:
                raise Stop(f"gate 6 study: score/{slug}/draws.tsv: a row on the {r['split']} split")
        if any(r["split"] not in SPLITS for r in data.score[slug]["datasets.tsv"]):
            raise Stop(f"gate 6 study: score/{slug}/datasets.tsv has a row on a split that is not counted")

    # 7 family
    for name in ("count.tsv", "count-no-constant.tsv", "count-no-unconverged.tsv"):
        fam = [r for r in data.count[name] if r["split"] == "lineage" and r["measure"] == PRIMARY]
        keys = sorted((pair_arm(r["pair"]), r["learner"]) for r in fam)
        if keys != sorted(FAMILY) or any(r["holm_p"] == "" for r in fam):
            raise Stop(f"gate 7 family: count/{name} holds {len(fam)} lineage pairs on {PRIMARY}, not the family's "
                       f"{len(FAMILY)} each with its Holm p")

    # 8 withdrawn
    for arm in sorted(data.withdrawn):
        for slug in species:
            for f in ("draws.tsv", "datasets.tsv"):
                if any(r["arm"] == arm for r in data.score[slug][f]):
                    raise Stop(f"gate 8 withdrawn: score/{slug}/{f} holds rows of {arm}, which was withdrawn")
            folder = top / "run" / slug / "outside" / arm
            if folder.is_dir() and any(folder.iterdir()):
                raise Stop(f"gate 8 withdrawn: run/{slug}/outside/{arm} holds reports of {arm}, which was withdrawn")
        for name in ("count.tsv", "count-no-constant.tsv", "count-no-unconverged.tsv"):
            for r in data.count[name]:
                if pair_arm(r["pair"]) == arm and (r["entered"] != "0" or r["datasets"] != "0"
                                                   or num(r["sign_p"]) != 1.0):
                    raise Stop(f"gate 8 withdrawn: count/{name}: the withdrawn {arm}'s pair enters datasets or is not "
                               "at p = 1")

    # 9 passes
    for slug in species:
        m = DONE_LINE.fullmatch((top / "run" / slug / "DONE").read_text(encoding="utf-8").strip())
        if m is None or m["species"] != slug:
            raise Stop(f"gate 9 passes: run/{slug}/DONE is no DONE line of {slug}")
        data.done[slug] = m.groupdict()
        if m["elasticnet"] == "leave-out":
            held = any(r["arm"] == "elasticnet" for f in ("draws.tsv", "datasets.tsv") for r in data.score[slug][f])
            folder = top / "run" / slug / "outside" / "elasticnet"
            if held or (folder.is_dir() and any(folder.iterdir())):
                raise Stop(f"gate 9 passes: {slug}'s run left the elastic net out, and its reading holds elastic net "
                           "results")
            if "elasticnet" not in data.withdrawn:
                raise Stop(f"gate 9 passes: {slug}'s DONE is of the first pass alone, which left the elastic net out, "
                           "and the elastic net is not withdrawn")

    # 10 rows
    for slug in species:
        rows = data.score[slug]["datasets.tsv"]
        keys = [(r["dataset"], r["split"], r["arm"], r["learner"]) for r in rows]
        want = expected_rows(data.by_species[slug], data.withdrawn)
        if len(set(keys)) != len(keys) or set(keys) != want:
            more = sorted({(k[0], k[2]) for k in set(keys) - want})
            less = sorted({(k[0], k[2]) for k in want - set(keys)})
            raise Stop(f"gate 10 rows: score/{slug}/datasets.tsv holds other rows than its first-run datasets'; "
                       f"more {more}, less {less}")
        draws = data.score[slug]["draws.tsv"]
        dkeys = [(r["dataset"], r["split"], int(r["seed"]), r["role"], r["arm"], r["learner"]) for r in draws]
        if len(set(dkeys)) != len(dkeys) or set(dkeys) != expected_draws(data.by_species[slug], data.withdrawn):
            raise Stop(f"gate 10 rows: score/{slug}/draws.tsv holds other draws than its first-run datasets'")

    # 11 missing
    for slug in species:
        for r in data.score[slug]["draws.tsv"]:
            if r["scored"] != "1" and (not r["why"] or any(r[m] != "" for m in MEASURES)):
                raise Stop(f"gate 11 missing: score/{slug}/draws.tsv: {r['dataset']} {r['arm']} seed {r['seed']}: a "
                           "draw not scored holds a figure")
        for r in data.score[slug]["datasets.tsv"]:
            if r["scored_draws"] != r["draws"] and (r["scored"] != "0" or any(r[m] != "" for m in MEASURES)):
                raise Stop(f"gate 11 missing: score/{slug}: {r['dataset']} {r['split']} {r['arm']} {r['learner']}: a "
                           "draw not scored is read as a value")

    # 12 draws
    for slug in species:
        files = data.score[slug]
        again = score.dataset_rows([typed_draw(r) for r in files["draws.tsv"]])
        again_at = {(r["dataset"], r["split"], r["arm"], r["learner"]): r for r in again}
        for r in files["datasets.tsv"]:
            g = again_at.get((r["dataset"], r["split"], r["arm"], r["learner"]))
            if g is None:
                raise Stop(f"gate 12 draws: score/{slug}: {r['dataset']} {r['split']} {r['arm']} {r['learner']} has "
                           "no draws")
            for c in score.DATASET_COLUMNS:
                if not same_cell(r[c], C9.fmt(g.get(c))):
                    raise Stop(f"gate 12 draws: score/{slug}: {r['dataset']} {r['split']} {r['arm']} {r['learner']}: "
                               f"{c} is not the mean of its draws")

    # 13 reports
    for slug in species:
        folder = top / "run" / slug / "outside"
        on_disk = sorted(p.relative_to(folder).as_posix() for p in folder.rglob("*") if p.is_file()) \
            if folder.is_dir() else []
        marks = {}
        for r in data.score[slug]["draws.tsv"]:
            if r["arm"] in ARMS and r["why"] != "draw missing":
                marks[report_name(r["arm"], r["dataset"], r["split"], int(r["seed"]))] = r
        if on_disk != sorted(marks):
            raise Stop(f"gate 13 reports: run/{slug}/outside holds {len(on_disk)} reports, its draws.tsv's library arm "
                       f"draws are {len(marks)}; more {len(set(on_disk) - set(marks))}, less "
                       f"{len(set(marks) - set(on_disk))}")
        data.reports[slug] = {}
        for name, r in sorted(marks.items()):
            with open(folder / name, encoding="utf-8") as handle:
                doc = json.load(handle)
            try:
                O.check_report(doc)
            except (O.Refuse, KeyError, TypeError) as e:
                raise Stop(f"gate 13 reports: run/{slug}/outside/{name}: {e}") from None
            if (doc.get("format"), doc.get("arm"), doc.get("pairs_with"), doc.get("dataset"), doc.get("scheme"),
                    doc.get("seed")) != (O.FORMAT, r["arm"], LEARNER_OF[r["arm"]], r["dataset"],
                                         scheme_of(r["split"]), int(r["seed"])):
                raise Stop(f"gate 13 reports: run/{slug}/outside/{name} is not the report of its draw")
            sums = fold_sums(doc)
            if any(doc.get(k) != v for k, v in sums.items()):
                raise Stop(f"gate 13 reports: run/{slug}/outside/{name}: a fold count is not the sum over its folds")
            if r["scored"] == "1" and (int(r["constant"]) != int(sums["constant_folds"] > 0)
                                       or int(r["unconverged"]) != int(sums["unconverged_refits"] > 0)):
                raise Stop(f"gate 13 reports: score/{slug}/draws.tsv marks {name} otherwise than its report")
            data.reports[slug][name] = doc

    # 14 checks
    for slug in species:
        checks = data.score[slug]["checks.tsv"]
        failed = sum(1 for r in checks if r["passed"] != "1")
        if not checks or failed:
            raise Stop(f"gate 14 checks: score/{slug}/checks.tsv: {failed} of {len(checks)} checks failed")
        text = (top / "run" / slug / "plans.log").read_text(encoding="utf-8").strip().splitlines()
        m = PLANS_LINE.fullmatch(text[-1]) if text else None
        # plans.py checks each lineage dataset seed a library arm ran: three checks at each, and one more at 42.
        lineage = {(r["dataset"], int(r["seed"])) for r in data.score[slug]["draws.tsv"]
                   if r["arm"] in ARMS and r["split"] == "lineage"}
        at42 = sum(1 for _, s in lineage if s == STUDY_SEED)
        if m is None or m["passed"] != m["made"] or int(m["seeds"]) != len(lineage) \
                or int(m["made"]) != 3 * len(lineage) + at42:
            raise Stop(f"gate 14 checks: run/{slug}/plans.log does not say that every plans check passed at each of "
                       f"its {len(lineage)} lineage dataset seeds")
        data.plans[slug] = {"seeds": len(lineage), "at42": at42,
                            "by_seed": {s: sum(1 for _, x in lineage if x == s) for s in ARM_SEEDS}}

    # 15 binary
    for slug in species:
        d = data.done[slug]
        if int(d["results"]) != len(data.reports[slug]):
            raise Stop(f"gate 15 binary: run/{slug}/DONE names {d['results']} results, its reports are "
                       f"{len(data.reports[slug])}")
        libs = read_libraries((top / "run" / slug / "libraries.txt").read_text(encoding="utf-8"))
        if sorted(libs) != sorted(LIBRARIES):
            raise Stop(f"gate 15 binary: run/{slug}/libraries.txt does not name {', '.join(LIBRARIES)} each once")
        data.libraries[slug] = libs
        data.timing[slug] = load(top / "run" / slug / "timing.tsv")
    return data


def codes(data):
    """{code sha256: [species slug]}: the binaries the species ran, by binary.txt's code line as each DONE names it."""
    out = {}
    for slug in data.species:
        out.setdefault(data.done[slug]["code"], []).append(slug)
    return out


def versions(data):
    """{(python, numpy, scipy, sklearn, xgboost, lightgbm): [species slug]}: the library versions each run read back."""
    out = {}
    for slug in data.species:
        out.setdefault(tuple(data.libraries[slug][k] for k in LIBRARIES), []).append(slug)
    return out


def meets(row):
    """(condition 1, condition 2) of one lineage primary row of count.tsv, as score.conditions reads them."""
    lo, holm = num(row["interval_lo"]), num(row["holm_p"])
    c1 = lo is not None and lo > 0
    c2 = holm is not None and holm < score.HOLM_ALPHA and int(row["better"]) > int(row["worse"])
    return c1, c2


def check_tables(data, tables):
    """Gates 15 and 16 on the page's own tables. Raises Stop naming the gate."""
    # 15 binary: one code and one version each, or each species named with its own.
    named = {(r["species_slug"], r["code"]) for r in tables["page-binaries.tsv"]}
    if named != {(slug, data.done[slug]["code"]) for slug in data.species}:
        raise Stop("gate 15 binary: page-binaries.tsv does not name each species with its code sha256")
    libs = tables["page-libraries.tsv"]
    if len(versions(data)) == 1:
        if len(libs) != 1 or libs[0]["species_count"] != len(data.species):
            raise Stop("gate 15 binary: the libraries agree in every species and are not written as one row")
    elif {(r["species_slug"], tuple(r[k] for k in LIBRARIES)) for r in libs} != {
            (slug, tuple(data.libraries[slug][k] for k in LIBRARIES)) for slug in data.species}:
        raise Stop("gate 15 binary: the library versions differ and page-libraries.tsv does not name each species")
    # 16 beside
    for r in tables["page-count.tsv"]:
        if (r["decides"] == 1) != (r["file"] == DECIDING_COUNT):
            raise Stop(f"gate 16 beside: page-count.tsv writes {r['file']} as "
                       f"{'deciding' if r['decides'] == 1 else 'deciding nothing'}")
    for r in tables["page-conditions.tsv"]:
        beside = r["reads"].startswith("beside")
        if str(r["decides"]) != {"1": "1", "0": "0"}.get(r["conditions_decides"], "x") or beside == (r["decides"] == 1):
            raise Stop(f"gate 16 beside: page-conditions.tsv: {r['arm']}: condition {r['condition']} is written "
                       "otherwise than conditions.tsv decides it")
    if any(r["decides"] != 0 for r in tables["page-species.tsv"]):
        raise Stop("gate 16 beside: page-species.tsv writes a species count as deciding")


def check_count(data):
    """Gate 17: the count files are score.py --count's from the datasets.tsv files. Raises Stop."""
    rows = [r for slug in data.species for r in data.score[slug]["datasets.tsv"]]
    withdrawn = sorted(data.withdrawn)
    counted = score.count(rows, withdrawn)
    want = {
        "count.tsv": counted,
        "count-no-constant.tsv": score.count(rows, withdrawn, keep_constant=False),
        "count-no-unconverged.tsv": score.count(rows, withdrawn, keep_unconverged=False),
    }
    want = {k: [{c: C9.fmt(r.get(c)) for c in score.COUNT_COLUMNS} for r in v] for k, v in want.items()}
    want["conditions.tsv"] = [{c: C9.fmt(r.get(c)) for c in score.CONDITION_COLUMNS} for r in score.conditions(counted)]
    for name, rows_want in want.items():
        if data.count[name] != rows_want:
            raise Stop(f"gate 17 count: count/{name} is not what score.py --count writes from the datasets.tsv files")


# ---------- the page's tables ----------

def page_conditions(data):
    out = []
    for r in data.count["conditions.tsv"]:
        out.append({"arm": r["arm"], "pair": pair_name(r["arm"], r["learner"]), "learner": r["learner"],
                    "condition": r["condition"], "decides": int(r["decides"] == "1"),
                    "conditions_decides": r["decides"], "reads": r["reads"], "met": r["met"],
                    "withdrawn": int(r["arm"] in data.withdrawn)})
    return out


CONDITION_TABLE_COLUMNS = ["arm", "pair", "learner", "condition", "decides", "conditions_decides", "reads", "met",
                           "withdrawn"]


def page_count(data):
    out = []
    for name in ("count.tsv", "count-no-constant.tsv", "count-no-unconverged.tsv"):
        for r in data.count[name]:
            arm = pair_arm(r["pair"])
            row = {"file": name, "decides": int(name == DECIDING_COUNT), **r,
                   "pair": pair_name(arm, r["learner"]), "against": against(arm), "meets_1": "", "meets_2": ""}
            if name == DECIDING_COUNT and r["split"] == "lineage" and r["measure"] == PRIMARY:
                c1, c2 = meets(r)
                row["meets_1"], row["meets_2"] = int(c1), int(c2)
            out.append(row)
    return out


COUNT_TABLE_COLUMNS = ["file", "decides", "split", "pair", "against", "learner", "measure", "orientation", "withdrawn",
                       "datasets", "entered", "not_scored", "constant", "unconverged", "better", "worse", "same",
                       "median_difference", "species", "interval_lo", "interval_hi", "sign_p", "holm_p", "meets_1",
                       "meets_2"]


def page_withdrawn(data):
    return [{"arm": arm, "pair": pair_name(arm, LEARNER_OF[arm]), "reason": reason, "time": time,
             "second_passes": data.second_passes if arm == "elasticnet" else 0}
            for arm, (reason, time) in sorted(data.withdrawn.items())]


WITHDRAWN_COLUMNS = ["arm", "pair", "reason", "time", "second_passes"]


def family_order(arm, learner):
    return FAMILY.index((arm, learner))


def page_datasets(data):
    """Per dataset, split and pair: the arm's balanced accuracy and recall, the default learner's at the arm's seeds
    (the twin's against the default's every draw), the differences, the draws scored and the dataset's genomes."""
    out = []
    for slug in data.species:
        cells = {(r["dataset"], r["split"], r["arm"], r["learner"]): r for r in data.score[slug]["datasets.tsv"]}
        keys = sorted((k for k in cells if k[2] in ARMS or k[2] == TWIN),
                      key=lambda k: (k[0], SPLITS.index(k[1]), family_order(k[2], k[3])))
        for d, split, arm, learner in keys:
            a, b = cells[(d, split, arm, learner)], cells[(d, split, against(arm), learner)]
            row = {"dataset": d, "species": a["species"], "antibiotic": a["antibiotic"], "genomes": data.genomes[d],
                   "split": split, "pair": pair_name(arm, learner), "against": against(arm), "learner": learner,
                   "draws": a["draws"], "arm_draws_scored": a["scored_draws"],
                   "default_draws_scored": b["scored_draws"], "constant_draws": a["constant_draws"],
                   "unconverged_draws": a["unconverged_draws"]}
            both = a["scored"] == "1" and b["scored"] == "1"
            for m, arm_col, default_col, diff_col in ((PRIMARY, "arm_value", "default_value", "difference"),
                                                      ("recall", "arm_recall", "default_recall", "recall_difference")):
                va, vb = num(a[m]), num(b[m])
                row[arm_col], row[default_col] = va, vb
                row[diff_col] = r12(va - vb) if both and va is not None and vb is not None else None
            out.append(row)
    return out


DATASET_TABLE_COLUMNS = ["dataset", "species", "antibiotic", "genomes", "split", "pair", "against", "learner", "draws",
                         "arm_draws_scored", "default_draws_scored", "constant_draws", "unconverged_draws", "arm_value",
                         "default_value", "difference", "arm_recall", "default_recall", "recall_difference"]


def tally(values, measure=PRIMARY):
    up = -1 if measure in score.LOWER_IS_BETTER else 1
    better = sum(1 for x in values if up * x > MARGIN)
    worse = sum(1 for x in values if up * x < -MARGIN)
    return better, worse, len(values) - better - worse


def page_species(data, dataset_rows):
    """Per species, split and pair, on the primary measure: datasets better, worse and the same of those entered.
    A row for every species and pair, a withdrawn pair's with no dataset. Decides nothing."""
    groups = {}
    for r in dataset_rows:
        groups.setdefault((r["dataset"].split("__")[0], r["split"], r["pair"]), []).append(r)
    out = []
    for slug in data.species:
        for split in SPLITS:
            for arm, learner in FAMILY:
                if split == "study" and arm == TWIN:
                    continue
                rs = groups.get((slug, split, pair_name(arm, learner)), [])
                values = [r["difference"] for r in rs if r["difference"] is not None]
                b, w, s = tally(values)
                out.append({"species": data.species_name[slug], "species_slug": slug, "split": split,
                            "pair": pair_name(arm, learner), "learner": learner, "decides": 0,
                            "withdrawn": int(arm in data.withdrawn), "datasets": len(rs), "entered": len(values),
                            "not_scored": len(rs) - len(values), "better": b, "worse": w, "same": s})
    return out


SPECIES_COLUMNS = ["species", "species_slug", "split", "pair", "learner", "decides", "withdrawn", "datasets", "entered",
                   "not_scored", "better", "worse", "same"]


def page_folds(data):
    """Per arm and split, the fold counts of its reports as sums over the folds, with the folds as denominator."""
    out = []
    for arm in ARMS:
        for split in SPLITS:
            docs = [doc for slug in data.species for name, doc in data.reports[slug].items()
                    if name.startswith(f"{arm}/") and doc["scheme"] == scheme_of(split)]
            row = {"arm": arm, "learner": LEARNER_OF[arm], "split": split, "withdrawn": int(arm in data.withdrawn),
                   "reports": len(docs), "reports_unavailable": sum(1 for d in docs if d.get("unavailable_reason")),
                   "folds": sum(len(d["folds"]) for d in docs)}
            for k in ("one_tree_folds", "constant_folds", "null_fit_folds", "unconverged_refits"):
                row[k] = sum(int(d[k]) for d in docs)
            row["reports_unconverged"] = sum(1 for d in docs if int(d["unconverged_refits"]) > 0)
            out.append(row)
    return out


FOLD_COLUMNS = ["arm", "learner", "split", "withdrawn", "reports", "reports_unavailable", "folds", "one_tree_folds",
                "constant_folds", "null_fit_folds", "unconverged_refits", "reports_unconverged"]
# What the printed page shows per arm: a column that is 0 by construction stays in the file and off the page.
PRINTED_FOLDS = {"xgboost": ("one_tree_folds", "constant_folds"), "lightgbm": ("one_tree_folds", "constant_folds"),
                 "elasticnet": ("unconverged_refits", "null_fit_folds", "constant_folds")}


def page_libraries(data):
    found = versions(data)
    if len(found) == 1:
        v, slugs = next(iter(found.items()))
        return [{"species": "every species", "species_slug": "all", "species_count": len(slugs),
                 **dict(zip(LIBRARIES, v))}]
    return [{"species": data.species_name[slug], "species_slug": slug, "species_count": 1,
             **{k: data.libraries[slug][k] for k in LIBRARIES}} for slug in data.species]


LIBRARY_COLUMNS = ["species", "species_slug", "species_count"] + list(LIBRARIES)


def page_binaries(data):
    return [{"species": data.species_name[slug], "species_slug": slug, "grm": data.done[slug]["grm"],
             "code": data.done[slug]["code"], "export_grm": data.done[slug]["export"],
             "elasticnet": data.done[slug]["elasticnet"], "results": data.done[slug]["results"]}
            for slug in data.species]


BINARY_COLUMNS = ["species", "species_slug", "grm", "code", "export_grm", "elasticnet", "results"]


def page_checks(data):
    """score.py's checks per name, passed of made over every species, and plans.py's from its count lines."""
    tally_ = {}
    for slug in data.species:
        for r in data.score[slug]["checks.tsv"]:
            p, m = tally_.get(("score.py", r["check"]), (0, 0))
            tally_[("score.py", r["check"])] = (p + int(r["passed"] == "1"), m + 1)
    seeds = sum(data.plans[s]["seeds"] for s in data.species)
    for s in ARM_SEEDS:
        n = sum(data.plans[slug]["by_seed"][s] for slug in data.species)
        tally_[("plans.py", PLANS_EXPORT.format(seed=s))] = (n, n)
    for name in PLANS_CHECKS:
        tally_[("plans.py", name)] = (seeds, seeds)
    at42 = sum(data.plans[s]["at42"] for s in data.species)
    tally_[("plans.py", PLANS_AT_42)] = (at42, at42)
    return [{"source": src, "check": name, "passed": p, "made": m} for (src, name), (p, m) in sorted(tally_.items())]


CHECK_COLUMNS = ["source", "check", "passed", "made"]


def fit_arm(fit):
    arm = fit.split("-", 1)[0]
    return arm if arm in ARMS else None


def page_times(data):
    """Per arm: the fits timed, and the sums and medians of their wall and CPU seconds."""
    by = {arm: [] for arm in ARMS}
    for slug in data.species:
        for r in data.timing[slug]:
            arm = fit_arm(r["fit"])
            if arm is None:
                raise Stop(f"gate 13 reports: run/{slug}/timing.tsv times a fit of no library arm")
            by[arm].append((float(r["wall_s"]), float(r["cpu_s"])))
    out = []
    for arm in ARMS:
        walls, cpus = [w for w, _ in by[arm]], [c for _, c in by[arm]]
        reports = sum(1 for slug in data.species for name in data.reports[slug] if name.startswith(f"{arm}/"))
        out.append({"arm": arm, "withdrawn": int(arm in data.withdrawn), "reports": reports, "fits": len(walls),
                    "wall_s_sum": r12(sum(walls)) if walls else None,
                    "wall_s_median": r12(statistics.median(walls)) if walls else None,
                    "cpu_s_sum": r12(sum(cpus)) if cpus else None,
                    "cpu_s_median": r12(statistics.median(cpus)) if cpus else None})
    return out


TIME_COLUMNS = ["arm", "withdrawn", "reports", "fits", "wall_s_sum", "wall_s_median", "cpu_s_sum", "cpu_s_median"]


def tables_of(data):
    datasets = page_datasets(data)
    return {
        "page-conditions.tsv": page_conditions(data),
        "page-count.tsv": page_count(data),
        "page-withdrawn.tsv": page_withdrawn(data),
        "page-species.tsv": page_species(data, datasets),
        "page-datasets.tsv": datasets,
        "page-folds.tsv": page_folds(data),
        "page-libraries.tsv": page_libraries(data),
        "page-binaries.tsv": page_binaries(data),
        "page-checks.tsv": page_checks(data),
        "page-times.tsv": page_times(data),
    }


TABLE_COLUMNS = {
    "page-conditions.tsv": CONDITION_TABLE_COLUMNS,
    "page-count.tsv": COUNT_TABLE_COLUMNS,
    "page-withdrawn.tsv": WITHDRAWN_COLUMNS,
    "page-species.tsv": SPECIES_COLUMNS,
    "page-datasets.tsv": DATASET_TABLE_COLUMNS,
    "page-folds.tsv": FOLD_COLUMNS,
    "page-libraries.tsv": LIBRARY_COLUMNS,
    "page-binaries.tsv": BINARY_COLUMNS,
    "page-checks.tsv": CHECK_COLUMNS,
    "page-times.tsv": TIME_COLUMNS,
}


# ---------- writing ----------

README = """# The outside-learner reading of plan block C7: the record

This directory holds the tables behind [docs/outside-learners-result.md](../../docs/outside-learners-result.md): plan block C7, XGBoost 3.0.5 and LightGBM 4.6.0 on the 4,000 rows `grm export` screens, each read against GRM's boosted trees (`gbt`), and a scikit-learn 1.7.2 elastic net on unitigs read against GRM's lasso (`lasso`), beside a twin, read once on the first run's {datasets} datasets in {species} species, {study} of them with a study split. Every row is an aggregate over a dataset, a species, a split, an arm, a pair or a fold count. No file holds a genome, a k-mer, a unitig, a sequence, a position, a feature name, an elastic net weight or a booster's feature importance.

`scripts/c7-outside/score.py --count` wrote the four files under `count/`. They are score.py's own tables: `scripts/c7-outside/page_tables.py write` counted them again from the species' `datasets.tsv` files (gate 17), copied them here as they are, and sums them in `SHA256SUMS`. It wrote every other file here from the run's own `score/` tables and each species' run files, after the {gates} gates of `gates.txt` passed. `checks.txt` is the printout of `scripts/c7-outside/run_checks.py check`. `reading.txt` lists the reading's commands, each with the sha256 of what it wrote. The rules of the reading were fixed before any C7 score was read (`docs/assumption-plan-corrections.md`, entry "C7 design").

{withdrawn}

{binaries}

{libraries}

| File | Holds |
|---|---|
| `gates.txt` | the gates the program passed before it wrote a table |
| `reading.txt` | the reading's commands and the sha256 of each one's output |
| `checks.txt` | the checks of the run, made before any score was read |
| `count/count.tsv` | the count that decides, per split, pair, learner and measure |
| `count/count-no-constant.tsv` | the count with the datasets that hold a constant fold left out, which decides nothing |
| `count/count-no-unconverged.tsv` | the count with the datasets whose elastic net draws hold a refit that reached 2,000 epochs left out, which decides nothing |
| `count/conditions.tsv` | the plan's conditions per arm; the arm against its own learner's twin pair alone stands beside with `decides` 0 |
| `page-conditions.tsv` | `conditions.tsv` per arm, with its pair and whether the arm was withdrawn |
| `page-count.tsv` | every row of the three counts, with `decides` 1 for `count.tsv` alone, the pair as the page names it, the row it is read against, and conditions 1 and 2 on the deciding lineage rows of the primary measure |
| `page-withdrawn.tsv` | each arm withdrawn, with the reason and the time the lead gave, and how many species had run a second pass with it |
| `page-species.tsv` | per species, named in full, split and pair: datasets better, worse and the same of those entered, which decides nothing |
| `page-datasets.tsv` | per dataset, split and pair, species and antibiotic in full: the arm's balanced accuracy and recall, the default learner's at the arm's seeds, the differences, the draws scored and the dataset's genomes |
| `page-folds.tsv` | per arm and split: the folds of its reports, and of them the boosters' one-tree folds and the elastic net's unconverged refits, null fits and constant folds |
| `page-libraries.tsv` | the library versions each species' run read back: one row when they agree, a row per species when they do not |
| `page-binaries.tsv` | per species: the binary, code and export binary sha256 its `DONE` names, the elastic net's pass, and its result count |
| `page-checks.tsv` | `score.py`'s checks per name and `plans.py`'s, passed of made |
| `page-times.tsv` | per arm: its reports, the fits timed, and the sums and medians of their wall and CPU seconds |

A dataset's figure is the mean of its draws: seeds 42, 43 and 44 on the lineage split, and one study run at seed 42. Each arm is read against the default's own learner at the arm's seeds (`k31-at-arm-seeds`); the twin, the default at seeds 45, 46 and 47, against the default. A difference is always the arm's value less the default's. A dataset with a draw not scored is not scored. "No difference shown" is the rule's reading and no statement that two learners are equal.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/c7-outside/page_tables.py print` prints the page's tables from these files.
"""


def withdrawn_sentence(data):
    if not data.withdrawn:
        return "No arm was withdrawn: all three arms ran."
    parts = "; ".join(f"`{arm}` ({ARM_WORDS[arm]}), at {time}: {reason}" for arm, (reason, time) in
                      sorted(data.withdrawn.items()))
    if not data.second_passes:
        return (f"Withdrawn before it ran, with no score read: {parts}. No species ran a second pass with it. Its pair "
                "stays in the family of five at p = 1, and it holds no result.")
    return (f"Withdrawn with no score read, after {data.second_passes} species had run a second pass with it: "
            f"{parts}. Each of those species is read from its first copy, whose `DONE` says the elastic net was "
            "left out; its second-pass copy is named in `checks.txt`, set aside and never opened. Its pair stays in "
            "the family of five at p = 1, and it holds no result.")


def binaries_sentence(data):
    found = codes(data)
    if len(found) == 1:
        code = next(iter(found))
        return (f"One binary: every species' `DONE` gives the code sha256 `{code}` of the C3 run's `binary.txt`, "
                f"{len(data.species)} of {len(data.species)} species.")
    parts = "; ".join(f"`{c}`: {', '.join(data.species_name[s] for s in sp)}" for c, sp in sorted(found.items()))
    return (f"{len(found)} binaries by `binary.txt`'s code sha256, each named with its species in "
            f"`page-binaries.tsv`: {parts}.")


def libraries_sentence(data):
    vs = versions(data)
    if len(vs) == 1:
        v = dict(zip(LIBRARIES, next(iter(vs))))
        return (f"One version of each library in every species: XGBoost {v['xgboost']}, LightGBM {v['lightgbm']}, "
                f"scikit-learn {v['sklearn']}, numpy {v['numpy']} and scipy {v['scipy']}, on python {v['python']}.")
    return (f"{len(vs)} sets of library versions, each species named with its own in `page-libraries.tsv`.")


def write_tsv(path, rows, columns):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for r in rows:
            handle.write("\t".join(C9.fmt(r.get(c)) for c in columns) + "\n")


def write(reading, l_path, out):
    reading, out = Path(reading), Path(out)
    lrows = load(l_path)
    if out.exists() and any(out.iterdir()):
        raise Stop(f"{out} is not empty")
    data = check_reading(reading, lrows)
    tables = tables_of(data)
    check_tables(data, tables)
    check_count(data)
    out.mkdir(parents=True, exist_ok=True)
    (out / "count").mkdir()
    for name in COUNT_FILES:
        (out / "count" / name).write_bytes((reading / "count" / name).read_bytes())
    for name in ("reading.txt", "checks.txt"):
        (out / name).write_bytes((reading / name).read_bytes())
    for name, rows in tables.items():
        write_tsv(out / name, rows, TABLE_COLUMNS[name])
    with open(out / "gates.txt", "w", encoding="utf-8") as handle:
        for n, name, words in GATES:
            handle.write(f"{n} {name}: {words}\n")
    run = data.run
    (out / "README.md").write_text(README.format(
        datasets=len(run), species=len(data.species), study=sum(1 for r in run if r["in_S"] == "1"), gates=len(GATES),
        withdrawn=withdrawn_sentence(data), binaries=binaries_sentence(data), libraries=libraries_sentence(data)),
        encoding="utf-8")
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


def yes(x, withdrawn="0"):
    if withdrawn == "1":
        return "withdrawn"
    return {"1": "met", "0": "not met", "": "the maintainer's"}[str(x)]


def table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return lines


def print_conditions(rows):
    out = [[f"`{r['pair']}`", r["condition"], "decides" if r["decides"] == "1" else "beside, decides nothing",
            r["reads"], yes(r["met"], r["withdrawn"] if r["condition"] != "6" else "0")] for r in rows]
    return table(["Pair", "Condition", "Role", "Reads", "Read"], out)


def print_count(rows, file="count.tsv", measure=PRIMARY, split="lineage"):
    out = []
    for r in rows:
        if r["file"] != file or r["split"] != split or r["measure"] != measure:
            continue
        cells = [f"`{r['pair']}`", f"`{r['against']}`", of(r["better"], r["entered"]), of(r["worse"], r["entered"]),
                 of(r["same"], r["entered"]), d3(r["median_difference"]), interval(r["interval_lo"], r["interval_hi"]),
                 pval(r["sign_p"]), pval(r["holm_p"])]
        if r["withdrawn"] == "1":
            cells[2:] = ["withdrawn"] * 3 + ["none", "none", pval(r["sign_p"]), pval(r["holm_p"])]
        if file == "count.tsv" and split == "lineage" and measure == PRIMARY:
            cells += [yes(r["meets_1"], r["withdrawn"]), yes(r["meets_2"], r["withdrawn"])]
        out.append(cells)
    header = ["Pair", "Against", "Better", "Worse", "Same", "Median", "Interval", "Sign p", "Holm p"]
    if file == "count.tsv" and split == "lineage" and measure == PRIMARY:
        header += ["Condition 1", "Condition 2"]
    return table(header, out)


def print_measures(rows):
    """count.tsv's other rows: every measure on the lineage split, and the study split, recall beside log loss."""
    out = []
    for r in rows:
        if r["file"] != "count.tsv" or (r["split"], r["measure"]) == ("lineage", PRIMARY):
            continue
        out.append([r["split"], f"`{r['pair']}`", MEASURE_WORDS[r["measure"]], r["orientation"],
                    "withdrawn" if r["withdrawn"] == "1" else of(r["entered"], r["datasets"]), r["better"],
                    r["worse"], r["same"], d3(r["median_difference"]), interval(r["interval_lo"], r["interval_hi"])])
    return table(["Split", "Pair", "Measure", "Better is", "Entered", "Better", "Worse", "Same", "Median", "Interval"],
                 out)


def print_withdrawn(rows):
    if not rows:
        return ["No arm was withdrawn."]
    return table(["Arm", "Pair", "Withdrawn at", "Reason", "Species that ran a second pass"],
                 [[ARM_WORDS[r["arm"]], f"`{r['pair']}`", r["time"], r["reason"], r["second_passes"]] for r in rows])


def print_species(rows):
    by, names = {}, {}
    for r in rows:
        if r["split"] != "lineage":
            continue
        names[r["species_slug"]] = r["species"]
        by[(r["species_slug"], r["pair"])] = r
    pairs = [pair_name(a, l) for a, l in FAMILY]
    out = []
    for slug in sorted(names):
        cells = [f"*{names[slug]}*"]
        for p in pairs:
            r = by[(slug, p)]
            cells.append("withdrawn" if r["withdrawn"] == "1" else
                         f"{r['better']} / {r['worse']} / {r['same']} of {r['entered']}")
        out.append(cells)
    return table(["Species"] + [f"`{p}`: better / worse / same" for p in pairs], out)


def print_folds(rows):
    out = []
    for r in rows:
        if r["withdrawn"] == "1":
            out.append([ARM_WORDS[r["arm"]], r["split"], "withdrawn", ""])
            continue
        counts = "; ".join(f"{k.replace('_', ' ')} {of(r[k], r['folds'])}" for k in PRINTED_FOLDS[r["arm"]])
        out.append([ARM_WORDS[r["arm"]], r["split"], count_words(r["folds"]), counts])
    return table(["Arm", "Split", "Folds", "Of the folds"], out)


def print_libraries(rows):
    return table(["Species", "python", "numpy", "scipy", "scikit-learn", "XGBoost", "LightGBM"],
                 [[r["species"] if r["species_slug"] == "all" else f"*{r['species']}*", r["python"], r["numpy"],
                   r["scipy"], r["sklearn"], r["xgboost"], r["lightgbm"]] for r in rows])


def print_binaries(rows):
    if len({r["code"] for r in rows}) == 1:
        r = rows[0]
        return [f"One code sha256 in {len(rows)} of {len(rows)} species: `{r['code'][:12]}`."]
    return table(["Species", "Code sha256", "Elastic net"],
                 [[f"*{r['species']}*", f"`{r['code'][:12]}`", r["elasticnet"]] for r in rows])


def print_checks(rows):
    return table(["Script", "Check", "Passed"], [[f"`{r['source']}`", r["check"], of(r["passed"], r["made"])]
                                                  for r in rows])


def print_times(rows):
    out = []
    for r in rows:
        if r["withdrawn"] == "1":
            out.append([ARM_WORDS[r["arm"]], "withdrawn", "", "", "", ""])
            continue
        out.append([ARM_WORDS[r["arm"]], count_words(r["fits"]), d3(r["wall_s_sum"]), d3(r["wall_s_median"]),
                    d3(r["cpu_s_sum"]), d3(r["cpu_s_median"])])
    return table(["Arm", "Fits", "Wall s, sum", "Wall s, median", "CPU s, sum", "CPU s, median"], out)


def print_datasets(rows):
    out = []
    for r in rows:
        if r["split"] != "lineage" and r["split"] != "study":
            continue
        out.append([f"*{r['species']}*", r["antibiotic"], count_words(r["genomes"]), r["split"], f"`{r['pair']}`",
                    f"{r['arm_draws_scored']} of {r['draws']}", d3(r["arm_value"]), d3(r["default_value"]),
                    d3(r["difference"]), d3(r["arm_recall"]), d3(r["default_recall"]), d3(r["recall_difference"])])
    return table(["Species", "Antibiotic", "Genomes", "Split", "Pair", "Draws scored", "Balanced accuracy",
                  "Against", "Difference", "Recall", "Against", "Difference"], out)


# The tables `print` prints, in its order; the page's part 2 placeholders name them.
PRINTS = ("conditions", "count", "measures", "no-constant", "no-unconverged", "withdrawn", "species", "folds",
          "libraries", "binaries", "checks", "times", "datasets")


def tables(folder):
    folder = Path(folder)
    count = load(folder / "page-count.tsv")
    return {
        "conditions": print_conditions(load(folder / "page-conditions.tsv")),
        "count": print_count(count),
        "measures": print_measures(count),
        "no-constant": print_count(count, file="count-no-constant.tsv"),
        "no-unconverged": print_count(count, file="count-no-unconverged.tsv"),
        "withdrawn": print_withdrawn(load(folder / "page-withdrawn.tsv")),
        "species": print_species(load(folder / "page-species.tsv")),
        "folds": print_folds(load(folder / "page-folds.tsv")),
        "libraries": print_libraries(load(folder / "page-libraries.tsv")),
        "binaries": print_binaries(load(folder / "page-binaries.tsv")),
        "checks": print_checks(load(folder / "page-checks.tsv")),
        "times": print_times(load(folder / "page-times.tsv")),
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
