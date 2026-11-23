"""The tables of docs/input-arms-result.md: plan block C3's one reading, as files and as the page prints them.

  python3 -I page_tables.py write <reading dir> <lists.tsv> <out dir>
  python3 -I page_tables.py print <out dir>

The reading directory holds:

  reading.txt          the reading's commands, each with the sha256 of what it wrote
  count/               score.py --count's count.tsv, count-no-constant.tsv and
                       conditions.tsv, taken once over the species' datasets.tsv
  score/<species>/     each species' run's own score/ (its score.tar): checks.tsv,
                       draws.tsv and datasets.tsv, as score.py wrote them
  run/<species>/       from each species' copy: binary.txt, DONE, and builds/,
                       the copy's results/builds/<kind>.build.json

<lists.tsv> is scripts/c3-inputs/lists.tsv: the first run is its rows with
first_run 1, which are L's rows not in M; the two program arms run on those
with first_run_A 1, the rows of A not in M; the study split runs on the
first run's rows in S.

`write` passes every gate of GATES before it writes a file, and stops at the
first that fails, naming it. The count files are score.py's own tables:
`write` checks them by counting again from the datasets.tsv files (gate 13),
copies them under count/ as they are, and sums them in SHA256SUMS with every
other file. It then writes the page's tables as page-*.tsv, reading.txt,
gates.txt, README.md and SHA256SUMS into <out dir> (absent or empty). Every
row is an aggregate over a dataset, a species, a split, an arm or a learner:
no file holds a genome, a k-mer, a unitig, a sequence, a position or a
column name of a matrix.

`print` prints the page's tables from <out dir> alone, figures rounded to 3
decimals.
"""

import csv
import hashlib
import importlib.util
import json
import re
import sys
from decimal import ROUND_HALF_EVEN, Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _score():
    spec = importlib.util.spec_from_file_location("c3_inputs_score_for_page", HERE / "score.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


score = _score()
C9 = score.C9

LEARNERS = score.LEARNERS
MEASURES = score.MEASURES
PRIMARY = score.PRIMARY
MARGIN = score.MARGIN
DEFAULT, TWIN = score.DEFAULT, score.TWIN
FAMILY = list(score.FAMILY)
# The five arms; the twin is the default fitted again at seeds of its own and is never one of them.
ARMS = [a for a in FAMILY if a != TWIN]
PROGRAM_ARMS = ("program", "k31-program")
BUILD_OF = dict(score.BUILD_OF)
SEEDS, TWIN_SEEDS, STUDY_SEED = (42, 43, 44), (45, 46, 47), 42
SPLITS = ("lineage", "study")
DECIDING, BESIDE = "ruling 15", "stricter"
CONSTANT_WHY = "constant model: refused for no candidate"
SCORE_FILES = ("checks.tsv", "draws.tsv", "datasets.tsv")
COUNT_FILES = ("count.tsv", "count-no-constant.tsv", "conditions.tsv")
BUILD_FORMAT = "grm-c3-build-1"
REVISION = score.REVISION
BUILD_COUNTS = ("genomes", "space_k", "matrix_support", "values_summed", "matrix_patterns", "matrix_values_kept",
                "matrix_values_rare", "matrix_values_core", "proteins", "graph_kmers", "unitigs", "pairs_held",
                "pairs_in_part", "kmers_without_unitig")
EQUALITY = ("kmers_grm", "kmers_outside", "kmers_both", "kmers_grm_only", "kmers_outside_only", "unitigs_grm",
            "unitigs_outside", "unitigs_both", "unitigs_grm_only", "unitigs_outside_only")
BUILD_NEEDS = {"counts": (), "aa": ("proteins",), "unitigs": ("graph", "graph_kmers", "unitigs", "pairs_held",
                                                               "pairs_in_part")}
GRAPHS = ("grm", "outside")
DONE_LINE = re.compile(r"DONE (?P<species>\S+) (?P<time>\S+) binary (?P<binary>[0-9a-f]{64}) "
                       r"code (?P<code>[0-9a-f]{64}) results (?P<results>\d+) SHA256SUMS (?P<sums>[0-9a-f]{64})")
GATES = [
    (1, "reading", "the reading holds reading.txt, count/ with count.tsv, count-no-constant.tsv and conditions.tsv, "
                   "score/ and run/, and nothing else; an arm count.tsv marks withdrawn is one of the five"),
    (2, "species", "score/ and run/ each hold exactly the species of the first run of lists.tsv, each once, with its "
                   "files, and each species' rows name it alone"),
    (3, "twin", "the twin's rows are lineage draws at seeds 45, 46 and 47 only, no arm's row is at a twin seed, and "
                "conditions.tsv reads the five arms and never the twin"),
    (4, "study", "a study row is a study run at seed 42 on the study split, and every lineage row is a draw at 42, 43 "
                 "or 44 or the twin's"),
    (5, "M", "no dataset of M, or outside the first run, stands in any table; species and antibiotic are lists.tsv's"),
    (6, "own datasets", "each arm's rows are on its own set's datasets only: the default, the twin and the three "
                        "matrix arms on L's first run, the two program arms on A's, the study split on those in S, "
                        "and a withdrawn arm on none"),
    (7, "draws", "datasets.tsv is the mean of draws.tsv's draws, and a dataset with a draw not scored is not scored, "
                 "never read as zero"),
    (8, "constant", "a constant-model draw is an arm's study run refused for no candidate, marked in draws.tsv and "
                    "counted in datasets.tsv, and its dataset is kept"),
    (9, "checks", "each species' checks.tsv holds checks and every one passed"),
    (10, "builds", "each species' builds are grm-c3-build-1 records of its species, one per kind of the matrix arms "
                   "run, holding counts only"),
    (11, "binary", "each species' binary.txt and DONE name the same binary and code; one code sha256 in every "
                   "species, or the page names each binary with its species"),
    (12, "conditions", "ruling 15's rows decide and the stricter reading's decide nothing, and conditions 1 and 2 as "
                       "read here off count.tsv equal conditions.tsv's"),
    (13, "count", "count.tsv, count-no-constant.tsv and conditions.tsv are what score.py --count writes from the "
                  "datasets.tsv files, with the arms count.tsv marks withdrawn"),
]
MEASURE_WORDS = {"balanced_accuracy": "balanced accuracy", "recall": "recall", "specificity": "specificity",
                 "vme": "very major error rate", "me": "major error rate", "log_loss": "log loss"}
ARM_WORDS = {"k31": "K31, the default", "unitigs": "unitigs", "k31-counts": "K31 plus k-mer counts",
             "k31-aa": "K31 plus amino-acid k-mers", "program": "the AMRFinderPlus program's calls alone",
             "k31-program": "K31 plus the AMRFinderPlus program's calls",
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


def first_run(lrows):
    return [r for r in lrows if r["first_run"] == "1"]


def in_a(r):
    return r["first_run_A"] == "1"


def arm_set(arm):
    return "A" if arm in PROGRAM_ARMS else "L"


def runs_on(arm, r):
    """Whether arm `arm` runs on the first-run row `r` of lists.tsv."""
    return arm not in PROGRAM_ARMS or in_a(r)


def pair_arm(pair):
    return pair[: -len(f"-{DEFAULT}")]


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
        self.run = first_run(lrows)
        self.species = sorted({r["species_slug"] for r in self.run})
        self.by_species = {s: [r for r in self.run if r["species_slug"] == s] for s in self.species}
        self.names = {r["dataset"]: (r["species"], r["antibiotic"]) for r in lrows}
        self.species_name = {r["species_slug"]: r["species"] for r in lrows}
        self.count = {}
        self.score = {}
        self.binary = {}
        self.done = {}
        self.builds = {}
        self.withdrawn = set()


# ---------- gates ----------

def expected_rows(rows, withdrawn):
    """(dataset, split, arm, learner) of datasets.tsv for the first-run rows of one species."""
    out = set()
    for r in rows:
        arms = [DEFAULT, TWIN] + [a for a in ARMS if a not in withdrawn and runs_on(a, r)]
        for arm in arms:
            for learner in LEARNERS:
                out.add((r["dataset"], "lineage", arm, learner))
                if r["in_S"] == "1" and arm != TWIN:
                    out.add((r["dataset"], "study", arm, learner))
    return out


def expected_draws(rows, withdrawn):
    """(dataset, split, seed, role, arm, learner) of draws.tsv for the first-run rows of one species."""
    out = set()
    for d, split, arm, learner in expected_rows(rows, withdrawn):
        if split == "study":
            out.add((d, split, STUDY_SEED, "study", arm, learner))
        elif arm == TWIN:
            out.update((d, split, s, "twin", arm, learner) for s in TWIN_SEEDS)
        else:
            out.update((d, split, s, "draw", arm, learner) for s in SEEDS)
    return out


def typed_draw(r):
    """A draws.tsv row as score.dataset_rows takes it."""
    out = dict(r)
    out["seed"] = int(r["seed"])
    out["scored"] = int(r["scored"])
    out["constant"] = int(r["constant"] or 0)
    for m in MEASURES:
        out[m] = num(r.get(m))
    for k in ("n", "tp", "fp", "tn", "fn"):
        out[k] = whole(r.get(k)) or 0
    return out


def read_binary(text):
    """{commit, grm, code, restarts, prodigal, bcalm} of one binary.txt."""
    out = {"restarts": 0, "prodigal": "", "bcalm": ""}
    for line in text.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] in ("commit", "grm", "code") and len(parts) == 2:
            out[parts[0]] = parts[1]
        elif parts[0] == "restart":
            out["restarts"] += 1
        elif parts[0] in ("prodigal", "bcalm") and len(parts) == 3:
            out[parts[0]] = f"{parts[1]} {parts[2]}"
    return out


def counts_only(rec):
    """Whether a build record holds the known fields only, each a count, a known word or the equality counts."""
    known = {"format", "evaluation_revision", "kind", "species", "graph", "equality"} | set(BUILD_COUNTS)
    if not set(rec) <= known:
        return False
    for k in BUILD_COUNTS:
        if k in rec and (isinstance(rec[k], bool) or not isinstance(rec[k], int) or rec[k] < 0):
            return False
    if "graph" in rec and rec["graph"] not in GRAPHS:
        return False
    if "equality" in rec:
        e = rec["equality"]
        if not isinstance(e, dict) or set(e) != set(EQUALITY):
            return False
        if any(isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in e.values()):
            return False
    return True


def check_reading(folder, lrows):
    """Gates 1 to 11, in order, on what the reading holds. Returns the Reading; raises Stop naming the gate."""
    data = Reading(folder, lrows)
    top = data.folder

    # 1 reading
    if not (top / "reading.txt").is_file():
        raise Stop("gate 1 reading: no reading.txt")
    for name in COUNT_FILES:
        if not (top / "count" / name).is_file():
            raise Stop(f"gate 1 reading: count/{name} is missing")
    for name in ("score", "run"):
        if not (top / name).is_dir():
            raise Stop(f"gate 1 reading: {name}/ is missing")
    extra = sorted(p.name for p in top.iterdir() if p.is_dir() and p.name not in ("count", "score", "run"))
    if extra:
        raise Stop(f"gate 1 reading: folders that are no part of the reading: {', '.join(extra)}")
    data.count = {name: load(top / "count" / name) for name in COUNT_FILES}
    data.withdrawn = {pair_arm(r["pair"]) for r in data.count["count.tsv"] if r["withdrawn"] == "1"}
    if not data.withdrawn <= set(ARMS):
        raise Stop(f"gate 1 reading: count.tsv marks withdrawn what is none of the five arms: "
                   f"{', '.join(sorted(data.withdrawn - set(ARMS)))}")

    # 2 species
    species = data.species
    for name, need in (("score", SCORE_FILES), ("run", ("binary.txt", "DONE"))):
        got = sorted(p.name for p in (top / name).iterdir() if p.is_dir())
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
    seen = {}
    for slug in species:
        for r in data.score[slug]["datasets.tsv"]:
            seen.setdefault(r["dataset"], set()).add(slug)
    if any(len(v) > 1 for v in seen.values()):
        raise Stop("gate 2 species: a dataset stands in two species' tables")

    # 3 twin
    for slug in species:
        for r in data.score[slug]["draws.tsv"]:
            seed, is_twin = int(r["seed"]), r["arm"] == TWIN
            if is_twin and not (r["split"] == "lineage" and r["role"] == "twin" and seed in TWIN_SEEDS):
                raise Stop(f"gate 3 twin: score/{slug}/draws.tsv: a twin row of {r['dataset']} is no lineage draw at "
                           "a twin seed")
            if not is_twin and (r["role"] == "twin" or seed in TWIN_SEEDS):
                raise Stop(f"gate 3 twin: score/{slug}/draws.tsv: {r['arm']} has a row of {r['dataset']} at a twin "
                           "seed or role")
        if any(r["arm"] == TWIN and r["split"] != "lineage" for r in data.score[slug]["datasets.tsv"]):
            raise Stop(f"gate 3 twin: score/{slug}/datasets.tsv has a twin row off the lineage split")
    cond_arms = [r["arm"] for r in data.count["conditions.tsv"]]
    if TWIN in cond_arms or sorted(set(cond_arms)) != sorted(ARMS):
        raise Stop("gate 3 twin: conditions.tsv does not read the five arms alone")

    # 4 study
    for slug in species:
        for r in data.score[slug]["draws.tsv"]:
            seed = int(r["seed"])
            if (r["role"] == "study") != (r["split"] == "study"):
                raise Stop(f"gate 4 study: score/{slug}/draws.tsv: a {r['role']} row of {r['dataset']} stands on "
                           f"the {r['split']} split")
            if r["split"] == "study" and seed != STUDY_SEED:
                raise Stop(f"gate 4 study: score/{slug}/draws.tsv: a study run of {r['dataset']} at seed {seed}")
            if r["split"] == "lineage" and r["arm"] != TWIN and not (r["role"] == "draw" and seed in SEEDS):
                raise Stop(f"gate 4 study: score/{slug}/draws.tsv: a lineage row of {r['dataset']} is no draw")
            if r["split"] not in SPLITS:
                raise Stop(f"gate 4 study: score/{slug}/draws.tsv: a row on the {r['split']} split")
        if any(r["split"] not in SPLITS for r in data.score[slug]["datasets.tsv"]):
            raise Stop(f"gate 4 study: score/{slug}/datasets.tsv has a row on a split that is not counted")

    # 5 M
    first = {r["dataset"] for r in data.run}
    in_m = {r["dataset"] for r in lrows if r["in_M"] == "1"}
    for slug in species:
        for f in ("draws.tsv", "datasets.tsv"):
            ds = {r["dataset"] for r in data.score[slug][f]}
            if ds & in_m:
                raise Stop(f"gate 5 M: score/{slug}/{f} holds {len(ds & in_m)} datasets of M")
            if ds - first:
                raise Stop(f"gate 5 M: score/{slug}/{f} holds {len(ds - first)} datasets outside the first run")
            for r in data.score[slug][f]:
                if (r["species"], r["antibiotic"]) != data.names[r["dataset"]]:
                    raise Stop(f"gate 5 M: score/{slug}/{f}: {r['dataset']} is named otherwise than in lists.tsv")

    # 6 own datasets
    for slug in species:
        rows = data.score[slug]["datasets.tsv"]
        keys = [(r["dataset"], r["split"], r["arm"], r["learner"]) for r in rows]
        want = expected_rows(data.by_species[slug], data.withdrawn)
        if len(set(keys)) != len(keys) or set(keys) != want:
            more = sorted({(k[0], k[2]) for k in set(keys) - want})
            less = sorted({(k[0], k[2]) for k in want - set(keys)})
            raise Stop(f"gate 6 own datasets: score/{slug}/datasets.tsv holds other rows than its arms' own datasets; "
                       f"more {len(more)} dataset and arm pairs, less {len(less)}")
        draws = data.score[slug]["draws.tsv"]
        dkeys = [(r["dataset"], r["split"], int(r["seed"]), r["role"], r["arm"], r["learner"]) for r in draws]
        if len(set(dkeys)) != len(dkeys) or set(dkeys) != expected_draws(data.by_species[slug], data.withdrawn):
            raise Stop(f"gate 6 own datasets: score/{slug}/draws.tsv holds other draws than its arms' own datasets'")

    # 7 draws
    for slug in species:
        files = data.score[slug]
        for r in files["datasets.tsv"]:
            if r["scored_draws"] != r["draws"] and (r["scored"] != "0" or any(r[m] != "" for m in MEASURES)):
                raise Stop(f"gate 7 draws: score/{slug}: {r['dataset']} {r['split']} {r['arm']} {r['learner']}: a "
                           "draw not scored is read as a value")
        again = score.dataset_rows([typed_draw(r) for r in files["draws.tsv"]])
        again_at = {(r["dataset"], r["split"], r["arm"], r["learner"]): r for r in again}
        for r in files["datasets.tsv"]:
            g = again_at[(r["dataset"], r["split"], r["arm"], r["learner"])]
            for c in score.DATASET_COLUMNS:
                if not same_cell(r[c], C9.fmt(g.get(c))):
                    raise Stop(f"gate 7 draws: score/{slug}: {r['dataset']} {r['split']} {r['arm']} {r['learner']}: "
                               f"{c} is not the mean of its draws")

    # 8 constant
    for slug in species:
        for r in data.score[slug]["draws.tsv"]:
            marked, why = r["constant"] == "1", r["why"] == CONSTANT_WHY
            if marked != why:
                raise Stop(f"gate 8 constant: score/{slug}/draws.tsv: {r['dataset']} {r['arm']}: a constant-model "
                           "draw and its reason do not go together")
            if marked and not (r["arm"] in ARMS and r["split"] == "study" and r["scored"] == "1"):
                raise Stop(f"gate 8 constant: score/{slug}/draws.tsv: {r['dataset']} {r['arm']} {r['split']}: a "
                           "constant-model draw that is no arm's scored study run")
        for r in data.score[slug]["datasets.tsv"]:
            if int(r["constant_draws"] or 0) > 0 and r["scored"] != "1":
                raise Stop(f"gate 8 constant: score/{slug}/datasets.tsv: {r['dataset']} {r['arm']}: a dataset with a "
                           "constant-model draw is not kept")

    # 9 checks
    for slug in species:
        checks = data.score[slug]["checks.tsv"]
        failed = sum(1 for r in checks if r["passed"] != "1")
        if not checks or failed:
            raise Stop(f"gate 9 checks: score/{slug}/checks.tsv: {failed} of {len(checks)} checks failed")

    # 10 builds
    kinds = sorted(BUILD_OF[a] for a in BUILD_OF if a not in data.withdrawn)
    for slug in species:
        folder = top / "run" / slug / "builds"
        got = sorted(p.name for p in folder.iterdir()) if folder.is_dir() else []
        if got != sorted(f"{k}.build.json" for k in kinds):
            raise Stop(f"gate 10 builds: run/{slug}/builds holds {got}, not one record per kind of {kinds}")
        data.builds[slug] = {}
        for kind in kinds:
            with open(folder / f"{kind}.build.json", encoding="utf-8") as handle:
                rec = json.load(handle)
            if (rec.get("format"), rec.get("evaluation_revision"), rec.get("kind"), rec.get("species")) != (
                    BUILD_FORMAT, REVISION, kind, slug):
                raise Stop(f"gate 10 builds: run/{slug}/builds/{kind}.build.json is not the {kind} build of {slug}")
            if not counts_only(rec) or any(k not in rec for k in BUILD_COUNTS[:8] + BUILD_NEEDS[kind]):
                raise Stop(f"gate 10 builds: run/{slug}/builds/{kind}.build.json holds more or less than its counts")
            data.builds[slug][kind] = rec

    # 11 binary
    for slug in species:
        b = read_binary((top / "run" / slug / "binary.txt").read_text(encoding="utf-8"))
        m = DONE_LINE.fullmatch((top / "run" / slug / "DONE").read_text(encoding="utf-8").strip())
        if not all(b.get(k) for k in ("commit", "grm", "code")):
            raise Stop(f"gate 11 binary: run/{slug}/binary.txt lacks its commit, binary or code line")
        if m is None or m["species"] != slug or (m["binary"], m["code"]) != (b["grm"], b["code"]):
            raise Stop(f"gate 11 binary: run/{slug}/DONE does not name {slug} and binary.txt's binary and code")
        data.binary[slug], data.done[slug] = b, m.groupdict()
    return data


def binaries(data):
    """{code sha256: [species slug]}: the binaries the species ran, by binary.txt's code line."""
    out = {}
    for slug in data.species:
        out.setdefault(data.binary[slug]["code"], []).append(slug)
    return out


def meets(row):
    """(condition 1, condition 2) of one lineage primary row of count.tsv, as score.conditions reads them."""
    lo, holm = num(row["interval_lo"]), num(row["holm_p"])
    c1 = lo is not None and lo > 0
    c2 = holm is not None and holm < score.HOLM_ALPHA and int(row["better"]) > int(row["worse"])
    return c1, c2


def check_tables(data, tables):
    """Gates 11 and 12 on the page's own tables. Raises Stop naming the gate."""
    # 11 binary: one code, or each binary named with its species.
    named = {(r["species_slug"], r["code"]) for r in tables["page-binaries.tsv"]}
    if named != {(slug, data.binary[slug]["code"]) for slug in data.species}:
        raise Stop("gate 11 binary: page-binaries.tsv does not name each species with its binary")
    # 12 conditions
    conds = data.count["conditions.tsv"]
    if {r["reading"] for r in conds} != {DECIDING, BESIDE}:
        raise Stop("gate 12 conditions: conditions.tsv holds other readings than ruling 15 and the stricter one")
    for r in tables["page-conditions.tsv"]:
        if (r["decides"] == 1) != (r["reading"] == DECIDING):
            raise Stop(f"gate 12 conditions: page-conditions.tsv: {r['arm']}: the {r['reading']} reading's condition "
                       f"{r['condition']} is written as {'deciding' if r['decides'] == 1 else 'deciding nothing'}")
    at = {(r["split"], r["pair"], r["learner"], r["measure"]): r for r in data.count["count.tsv"]}

    def met(arm):
        rows = [at[("lineage", f"{arm}-{DEFAULT}", l, PRIMARY)] for l in LEARNERS]
        return [l for l, r in zip(LEARNERS, rows) if meets(r)[0]], [l for l, r in zip(LEARNERS, rows) if meets(r)[1]]

    twin_one, twin_two = met(TWIN)
    for arm in ARMS:
        one, two = met(arm)
        mine = {(r["reading"], r["condition"]): r for r in conds if r["arm"] == arm}
        want = {(DECIDING, "1 and 2"): len(set(one) & set(two)), (DECIDING, "3"): max(len(twin_one), len(twin_two)),
                (BESIDE, "1"): len(one), (BESIDE, "2"): len(two)}
        for k, n in want.items():
            if k not in mine or whole(mine[k]["learners"]) != n:
                raise Stop(f"gate 12 conditions: {arm}: the {k[0]} reading's condition {k[1]} read here off count.tsv "
                           "differs from conditions.tsv")


def check_count(data):
    """Gate 13: the count files are score.py --count's from the datasets.tsv files. Raises Stop."""
    rows = [r for slug in data.species for r in data.score[slug]["datasets.tsv"]]
    withdrawn = sorted(data.withdrawn)
    counted = score.count(rows, withdrawn)
    want = {
        "count.tsv": [{c: C9.fmt(r.get(c)) for c in score.COUNT_COLUMNS} for r in counted],
        "count-no-constant.tsv": [{c: C9.fmt(r.get(c)) for c in score.COUNT_COLUMNS}
                                  for r in score.count(rows, withdrawn, keep_constant=False)],
        "conditions.tsv": [{c: C9.fmt(r.get(c)) for c in score.CONDITION_COLUMNS} for r in score.conditions(counted)],
    }
    for name, rows_want in want.items():
        if data.count[name] != rows_want:
            raise Stop(f"gate 13 count: count/{name} is not what score.py --count writes from the datasets.tsv files")


# ---------- the page's tables ----------

def page_conditions(data):
    return [{"arm": r["arm"], "reading": r["reading"], "decides": int(r["reading"] == DECIDING),
             "condition": r["condition"], "reads": r["reads"], "met": r["met"], "learners": r["learners"]}
            for r in data.count["conditions.tsv"]]


CONDITION_TABLE_COLUMNS = ["arm", "reading", "decides", "condition", "reads", "met", "learners"]


def page_count(data):
    out = []
    for name, decides in (("count.tsv", 1), ("count-no-constant.tsv", 0)):
        for r in data.count[name]:
            row = {"file": name, "decides": decides, **r, "meets_1": "", "meets_2": ""}
            if r["split"] == "lineage" and r["measure"] == PRIMARY:
                c1, c2 = meets(r)
                row["meets_1"], row["meets_2"] = int(c1), int(c2)
            out.append(row)
    return out


COUNT_TABLE_COLUMNS = ["file", "decides"] + score.COUNT_COLUMNS + ["meets_1", "meets_2"]


def arm_order(arm):
    return ([DEFAULT] + FAMILY).index(arm)


def page_datasets(data):
    """Per dataset, split, arm row and learner: the arm's figure and the default's, the difference, the draws."""
    out = []
    for slug in data.species:
        cells = {(r["dataset"], r["split"], r["arm"], r["learner"]): r for r in data.score[slug]["datasets.tsv"]}
        keys = sorted((k for k in cells if k[2] != DEFAULT),
                      key=lambda k: (k[0], SPLITS.index(k[1]), arm_order(k[2]), LEARNERS.index(k[3])))
        for d, split, arm, learner in keys:
            a, b = cells[(d, split, arm, learner)], cells[(d, split, DEFAULT, learner)]
            row = {"dataset": d, "species": a["species"], "antibiotic": a["antibiotic"], "split": split,
                   "pair": f"{arm}-{DEFAULT}", "learner": learner, "draws": a["draws"],
                   "arm_draws_scored": a["scored_draws"], "default_draws_scored": b["scored_draws"],
                   "constant_draws": a["constant_draws"],
                   "study_refused": int(split == "study" and int(a["constant_draws"] or 0) > 0)}
            both = a["scored"] == "1" and b["scored"] == "1"
            for m, arm_col, default_col, diff_col in ((PRIMARY, "arm_value", "default_value", "difference"),
                                                      ("recall", "arm_recall", "default_recall", "recall_difference")):
                va, vb = num(a[m]), num(b[m])
                row[arm_col], row[default_col] = va, vb
                row[diff_col] = r12(va - vb) if both and va is not None and vb is not None else None
            out.append(row)
    return out


DATASET_TABLE_COLUMNS = ["dataset", "species", "antibiotic", "split", "pair", "learner", "draws", "arm_draws_scored",
                         "default_draws_scored", "constant_draws", "study_refused", "arm_value", "default_value",
                         "difference", "arm_recall", "default_recall", "recall_difference"]


def tally(values, measure=PRIMARY):
    up = -1 if measure in score.LOWER_IS_BETTER else 1
    better = sum(1 for x in values if up * x > MARGIN)
    worse = sum(1 for x in values if up * x < -MARGIN)
    return better, worse, len(values) - better - worse


def median(values):
    import numpy as np
    return r12(float(np.median(values))) if values else None


def page_species_count(data, dataset_rows):
    """Per species, split, arm row and learner, on the primary measure: better, worse and same, and the median."""
    out = []
    groups = {}
    for r in dataset_rows:
        slug = r["dataset"].split("__")[0]
        groups.setdefault((slug, r["split"], r["pair"], r["learner"]), []).append(r)
    order = lambda k: (k[0], SPLITS.index(k[1]), arm_order(pair_arm(k[2])), LEARNERS.index(k[3]))
    for key in sorted(groups, key=order):
        slug, split, pair, learner = key
        rs = groups[key]
        values = [r["difference"] for r in rs if r["difference"] is not None]
        b, w, s = tally(values)
        out.append({"species": data.species_name[slug], "species_slug": slug, "split": split, "pair": pair,
                    "learner": learner, "datasets": len(rs), "entered": len(values), "not_scored": len(rs) - len(values),
                    "constant": sum(1 for r in rs if r["difference"] is not None and int(r["constant_draws"] or 0) > 0),
                    "better": b, "worse": w, "same": s, "median": median(values)})
    return out


SPECIES_COUNT_COLUMNS = ["species", "species_slug", "split", "pair", "learner", "datasets", "entered", "not_scored",
                         "constant", "better", "worse", "same", "median"]


def page_species(data):
    out = []
    for slug in data.species:
        mine = data.by_species[slug]
        rows = data.score[slug]["datasets.tsv"]
        out.append({"species": data.species_name[slug], "species_slug": slug, "datasets": len(mine),
                    "a_datasets": sum(1 for r in mine if in_a(r)),
                    "study_datasets": sum(1 for r in mine if r["in_S"] == "1"),
                    "study_a_datasets": sum(1 for r in mine if r["in_S"] == "1" and in_a(r)),
                    "m_datasets": sum(1 for r in data.lrows if r["species_slug"] == slug and r["in_M"] == "1"),
                    "constant_datasets": len({(r["dataset"], r["arm"]) for r in rows if int(r["constant_draws"] or 0) > 0}),
                    "default_lineage_not_scored": len({r["dataset"] for r in rows if r["arm"] == DEFAULT
                                                       and r["split"] == "lineage" and r["scored"] != "1"}),
                    "code": data.binary[slug]["code"]})
    return out


SPECIES_COLUMNS = ["species", "species_slug", "datasets", "a_datasets", "study_datasets", "study_a_datasets",
                   "m_datasets", "constant_datasets", "default_lineage_not_scored", "code"]


def page_arms(data):
    """Per arm row: its set, its first-run datasets and species, its study datasets, and its refused study runs."""
    out = []
    for arm in [DEFAULT] + FAMILY:
        mine = [r for r in data.run if runs_on(arm, r)]
        withdrawn = arm in data.withdrawn
        study = [r for r in mine if r["in_S"] == "1"] if arm != TWIN else []
        rows = [r for slug in data.species for r in data.score[slug]["datasets.tsv"] if r["arm"] == arm]
        out.append({"arm": arm, "words": ARM_WORDS[arm], "set": arm_set(arm), "withdrawn": int(withdrawn),
                    "datasets": 0 if withdrawn else len(mine),
                    "species": 0 if withdrawn else len({r["species_slug"] for r in mine}),
                    "study_datasets": 0 if withdrawn else len(study),
                    "study_species": 0 if withdrawn else len({r["species_slug"] for r in study}),
                    "build": BUILD_OF.get(arm, ""),
                    "study_refused": len({r["dataset"] for r in rows if r["split"] == "study"
                                          and int(r["constant_draws"] or 0) > 0})})
    return out


ARM_COLUMNS = ["arm", "words", "set", "withdrawn", "datasets", "species", "study_datasets", "study_species", "build",
               "study_refused"]


def page_builds(data):
    out = []
    for slug in data.species:
        for kind in sorted(data.builds[slug]):
            rec = data.builds[slug][kind]
            arm = [a for a, k in BUILD_OF.items() if k == kind][0]
            row = {"species": data.species_name[slug], "species_slug": slug, "kind": kind, "arm": arm,
                   "graph": rec.get("graph")}
            row.update({k: rec.get(k) for k in BUILD_COUNTS})
            row.update({f"equality_{k}": (rec.get("equality") or {}).get(k) for k in EQUALITY})
            out.append(row)
    return out


BUILD_COLUMNS = ["species", "species_slug", "kind", "arm", "graph"] + list(BUILD_COUNTS) + [f"equality_{k}" for k in
                                                                                           EQUALITY]


def page_binaries(data):
    return [{"species": data.species_name[slug], "species_slug": slug, **{k: data.binary[slug][k] for k in
            ("commit", "grm", "code", "restarts", "prodigal", "bcalm")}, "results": data.done[slug]["results"]}
            for slug in data.species]


BINARY_COLUMNS = ["species", "species_slug", "commit", "grm", "code", "restarts", "prodigal", "bcalm", "results"]


def tables_of(data):
    datasets = page_datasets(data)
    return {
        "page-conditions.tsv": page_conditions(data),
        "page-count.tsv": page_count(data),
        "page-arms.tsv": page_arms(data),
        "page-species.tsv": page_species(data),
        "page-species-count.tsv": page_species_count(data, datasets),
        "page-datasets.tsv": datasets,
        "page-builds.tsv": page_builds(data),
        "page-binaries.tsv": page_binaries(data),
    }


TABLE_COLUMNS = {
    "page-conditions.tsv": CONDITION_TABLE_COLUMNS,
    "page-count.tsv": COUNT_TABLE_COLUMNS,
    "page-arms.tsv": ARM_COLUMNS,
    "page-species.tsv": SPECIES_COLUMNS,
    "page-species-count.tsv": SPECIES_COUNT_COLUMNS,
    "page-datasets.tsv": DATASET_TABLE_COLUMNS,
    "page-builds.tsv": BUILD_COLUMNS,
    "page-binaries.tsv": BINARY_COLUMNS,
}


# ---------- writing ----------

README = """# The input-arm reading of plan block C3: the record

This directory holds the tables behind [docs/input-arms-result.md](../../docs/input-arms-result.md): plan block C3, five other inputs beside the presence of 31-mers (K31), each against K31 refitted in the same run, read once on the first run's {datasets} datasets in {species} species, {study} of them with a study split. The two program arms run on A's first run, {a_datasets} datasets in {a_species} species; the three matrix arms, the default and the twin on L's. Every row is an aggregate over a dataset, a species, a split, an arm or a learner. No file holds a genome, a k-mer, a unitig, a sequence, a position or a column name of a matrix.

`scripts/c3-inputs/score.py --count` wrote the three files under `count/`. They are score.py's own tables: `scripts/c3-inputs/page_tables.py write` counted them again from the species' `datasets.tsv` files (gate 13), copied them here as they are, and sums them in `SHA256SUMS`. It wrote every other file here from the run's own `score/` tables and each species' run files, after the {gates} gates of `gates.txt` passed. `reading.txt` lists the reading's commands, each with the sha256 of what it wrote. The rules of the reading were fixed before any C3 score was read (`docs/assumption-plan-corrections.md`, entries "C3 design" and "C3 arms").

{binaries}

| File | Holds |
|---|---|
| `gates.txt` | the gates the program passed before it wrote a table |
| `reading.txt` | the reading's commands and the sha256 of each one's output |
| `count/count.tsv` | the count that decides, per split, arm row, learner and measure |
| `count/count-no-constant.tsv` | the count with the constant-model datasets left out, which decides nothing |
| `count/conditions.tsv` | the plan's conditions per arm, read by ruling 15, which decides, and by the stricter reading, which decides nothing |
| `page-conditions.tsv` | `conditions.tsv` with a column `decides`: 1 for ruling 15's rows, 0 for the stricter reading's |
| `page-count.tsv` | every row of both counts, with `decides`, and conditions 1 and 2 read on the lineage primary rows |
| `page-arms.tsv` | per arm row: its set, its first-run datasets and species, its study datasets, whether it was withdrawn, its build and its study runs refused for no candidate |
| `page-species.tsv` | per species, named in full: datasets in L, A and S, datasets of M waiting, constant-model datasets and the binary's code sha256 |
| `page-species-count.tsv` | per species, split, arm row and learner: better, worse and same on the primary measure, and the median difference |
| `page-datasets.tsv` | per dataset, split, arm row and learner, species and antibiotic in full: the arm's balanced accuracy and recall, the default's, the differences, the draws scored, constant-model draws and a study run refused |
| `page-builds.tsv` | per species and matrix build, the counts its `build.json` gives, and nothing else |
| `page-binaries.tsv` | per species: the commit, binary and code sha256 of its run, restarts, and the outside programs' versions and sha256 |

A dataset's figure is the mean of its draws: seeds 42, 43 and 44 on the lineage split, the twin's 45, 46 and 47, and one study run at the release's seeds. A difference is always the arm's value less the default's. A dataset with a draw not scored is not scored. An arm's study run refused for no candidate is scored as the constant model on the default's folds (ruling 16), enters the count, and is marked in `constant`. "No difference shown" is the rule's reading and no statement that two inputs are equal.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/c3-inputs/page_tables.py print` prints the page's tables from these files.
"""


def binaries_sentence(data):
    codes = binaries(data)
    if len(codes) == 1:
        code = next(iter(codes))
        return f"One binary: every species' `binary.txt` gives the code sha256 `{code}`, {len(data.species)} of {len(data.species)} species."
    parts = "; ".join(f"`{c}`: {', '.join(data.species_name[s] for s in sp)}" for c, sp in sorted(codes.items()))
    return f"{len(codes)} binaries by `binary.txt`'s code sha256, each named with its species in `page-binaries.tsv`: {parts}."


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
    (out / "reading.txt").write_bytes((reading / "reading.txt").read_bytes())
    for name, rows in tables.items():
        write_tsv(out / name, rows, TABLE_COLUMNS[name])
    with open(out / "gates.txt", "w", encoding="utf-8") as handle:
        for n, name, words in GATES:
            handle.write(f"{n} {name}: {words}\n")
    run = data.run
    a_run = [r for r in run if in_a(r)]
    (out / "README.md").write_text(README.format(
        datasets=len(run), species=len(data.species), study=sum(1 for r in run if r["in_S"] == "1"),
        a_datasets=len(a_run), a_species=len({r["species_slug"] for r in a_run}), gates=len(GATES),
        binaries=binaries_sentence(data)), encoding="utf-8")
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


def yes(x):
    return {"1": "met", "0": "not met", "": "the maintainer's"}[str(x)]


def table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return lines


def pair_words(pair):
    return ARM_WORDS[pair_arm(pair)]


def print_conditions(rows, reading):
    out = [[ARM_WORDS[r["arm"]], r["condition"], r["reads"], yes(r["met"]), r["learners"] or ""]
           for r in rows if r["reading"] == reading]
    return table(["Arm", "Condition", "Reads", "Read", "Learners"], out)


def print_count(rows, file="count.tsv"):
    out = []
    for r in rows:
        if r["file"] != file or r["split"] != "lineage" or r["measure"] != PRIMARY:
            continue
        cells = [pair_words(r["pair"]), f"`{r['learner']}`", of(r["better"], r["entered"]), of(r["worse"], r["entered"]),
                 of(r["same"], r["entered"]), d3(r["median_difference"]), interval(r["interval_lo"], r["interval_hi"]),
                 pval(r["sign_p"]), pval(r["holm_p"]), r["constant"]]
        if file == "count.tsv":
            cells += [yes(r["meets_1"]), yes(r["meets_2"])]
        out.append(cells)
    header = ["Arm", "Learner", "Better", "Worse", "Same", "Median", "Interval", "Sign p", "Holm p", "Constant"]
    return table(header + (["Condition 1", "Condition 2"] if file == "count.tsv" else []), out)


def print_beside(rows):
    out = []
    for r in rows:
        if r["file"] != "count.tsv" or (r["split"], r["measure"]) == ("lineage", PRIMARY):
            continue
        if r["split"] == "study" and r["measure"] != PRIMARY:
            continue
        out.append([r["split"], pair_words(r["pair"]), f"`{r['learner']}`", MEASURE_WORDS[r["measure"]],
                    r["orientation"], of(r["entered"], r["datasets"]), r["better"], r["worse"], r["same"],
                    d3(r["median_difference"]), interval(r["interval_lo"], r["interval_hi"])])
    return table(["Split", "Arm", "Learner", "Measure", "Better is", "Entered", "Better", "Worse", "Same", "Median",
                  "Interval"], out)


def print_arms(rows):
    out = [[r["words"], r["set"], "withdrawn" if r["withdrawn"] == "1" else "ran", count_words(r["datasets"]),
            r["species"], r["study_datasets"], r["study_refused"]] for r in rows]
    return table(["Arm", "Set", "State", "Datasets", "Species", "Study datasets", "Study runs refused"], out)


def print_species(rows, counts):
    by = {}
    for r in counts:
        if r["split"] == "lineage":
            k = (r["species_slug"], r["pair"])
            b, w, s = by.get(k, (0, 0, 0))
            by[k] = (b + int(r["better"]), w + int(r["worse"]), s + int(r["same"]))
    out = []
    for r in rows:
        cells = [f"*{r['species']}*", r["datasets"], r["a_datasets"], r["study_datasets"], r["m_datasets"]]
        for arm in FAMILY:
            b, w, s = by.get((r["species_slug"], f"{arm}-{DEFAULT}"), (0, 0, 0))
            cells.append(f"{b} / {w} / {s}")
        out.append(cells)
    return table(["Species", "Datasets", "In A", "In S", "M waiting"] + [f"{ARM_WORDS[a]}: better / worse / same"
                                                                        for a in FAMILY], out)


def print_builds(rows):
    out = [[f"*{r['species']}*", r["kind"], count_words(r["genomes"]), count_words(r["values_summed"]),
            count_words(r["matrix_patterns"]), count_words(r["proteins"]) if r["proteins"] else "",
            count_words(r["unitigs"]) if r["unitigs"] else "",
            count_words(r["pairs_in_part"]) if r["pairs_in_part"] else ""] for r in rows]
    return table(["Species", "Build", "Genomes", "Values summed", "Patterns", "Proteins", "Unitigs",
                  "Pairs held in part"], out)


def print_binaries(rows):
    out = [[f"*{r['species']}*", f"`{r['commit'][:7]}`", f"`{r['code'][:12]}`", r["restarts"]] for r in rows]
    return table(["Species", "Commit", "Code sha256", "Restarts"], out)


def tables(folder):
    folder = Path(folder)
    count = load(folder / "page-count.tsv")
    conds = load(folder / "page-conditions.tsv")
    return {
        "conditions": print_conditions(conds, DECIDING),
        "stricter": print_conditions(conds, BESIDE),
        "count": print_count(count),
        "no-constant": print_count(count, file="count-no-constant.tsv"),
        "beside": print_beside(count),
        "arms": print_arms(load(folder / "page-arms.tsv")),
        "species": print_species(load(folder / "page-species.tsv"), load(folder / "page-species-count.tsv")),
        "builds": print_builds(load(folder / "page-builds.tsv")),
        "binaries": print_binaries(load(folder / "page-binaries.tsv")),
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
