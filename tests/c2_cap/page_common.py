"""SYNTHETIC readings for the result page's scripts (scripts/c2-cap/page_tables.py and run_checks.py).

Every figure, call, count and checksum here is made up in the test. The 8 species slugs and names are those of
scripts/c2-cap/species.tsv, so that the reading has the run's shape; their datasets are made up (drug_a, drug_b,
drug_m), and no figure comes from a run.
"""

import hashlib
import json

from tests.c2_cap import common

score = common.load("score")
expect = common.load("expect")
C3, C9 = score.C3, score.C9
LEARNERS = score.LEARNERS
SPECIES = [line.split("\t")[:2] for line in (common.KIT / "species.tsv").read_text(encoding="utf-8").splitlines()[1:]]
LIST_COLUMNS = ["dataset", "species_slug", "species", "antibiotic", "samples", "resistant", "susceptible", "lineages",
                "in_M", "in_S", "first_run", "labels_depend_on"]
ADDED_COLUMNS = ["dataset", "species", "antibiotic", "arm", "labelled", "added", "added_labelled", "added_resistant",
                 "added_susceptible", "none_added"]
EXPECT_COLUMNS = ["dataset", "species_slug", "arm", "split", "slot", "role", "fit"]
COMMIT, GRM, CODE = "c" * 40, "a" * 64, "b" * 64


def tsv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            f.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in columns) + "\n")


def read(path):
    return C9.read_tsv(path)


def calls(arm, learner):
    """Made-up summed calls: the default and the twin at 0.8 balanced accuracy, cap-all at 0.84."""
    up = 2 if arm == "cap-all" else 0
    c = {"n": 100, "tp": 40 + up, "fn": 10 - up, "tn": 40 + up, "fp": 10 - up, "loss_sum": 40.0 - up}
    return C9.measures(c)


class Reading:
    """A SYNTHETIC reading: lists.tsv, runs/<slug>/ in run_species.sh's OUT layout and count/.

    Each species has two first-run datasets, drug_a and drug_b, and one in M, drug_m. The first two species'
    drug_a is in S. The second species' drug_b adds no labelled genome under cap-all (none-added); the third
    species' drug_a adds 5 (fewer than 20); the fourth species' drug_b has its lineage split unavailable.
    """

    def __init__(self, tmp):
        self.tmp = tmp
        self.lists = tmp / "lists.tsv"
        self.runs = tmp / "runs"
        self.count = tmp / "count"
        self.slugs = [s for s, _ in SPECIES]
        self.name = dict(SPECIES)
        self.none_added = {(self.slugs[1], f"{self.slugs[1]}__drug_b")}
        self.few = f"{self.slugs[2]}__drug_a"
        self.refused = f"{self.slugs[3]}__drug_b"
        self.study = {f"{s}__drug_a" for s in self.slugs[:2]}
        self.binary = {s: (COMMIT, GRM, CODE) for s in self.slugs}

    def list_rows(self):
        out = []
        for slug in self.slugs:
            for drug, m in (("drug_a", 0), ("drug_b", 0), ("drug_m", 1)):
                d = f"{slug}__{drug}"
                out.append({"dataset": d, "species_slug": slug, "species": self.name[slug],
                            "antibiotic": drug.replace("_", " "), "samples": 300, "resistant": 100,
                            "susceptible": 200, "lineages": 20, "in_M": m, "in_S": int(d in self.study),
                            "first_run": 1 - m, "labels_depend_on": "none of the four"})
        return out

    def datasets(self, slug):
        return [f"{slug}__drug_a", f"{slug}__drug_b"]

    def added_rows(self, slug):
        out = []
        for drug in ("drug_a", "drug_b", "drug_m"):
            d = f"{slug}__{drug}"
            for arm in ("cap-4000", "cap-all"):
                n = 0 if (slug, d) in self.none_added else 5 if d == self.few else 400
                out.append({"dataset": d, "species": self.name[slug], "antibiotic": drug.replace("_", " "), "arm": arm,
                            "labelled": 300 + n, "added": 2096 if arm == "cap-all" else 2000, "added_labelled": n,
                            "added_resistant": n // 2, "added_susceptible": n - n // 2, "none_added": int(n == 0)})
        return out

    def expect_rows(self, slug):
        study = [d for d in self.datasets(slug) if d in self.study]
        empty = {("cap-all", d) for s, d in self.none_added if s == slug}
        roots = {"fits": {"cap-all": "cap-all"}, "same_as": {}}
        rows = expect.rows(slug, self.datasets(slug), study, 3, 3, ["cap-all"], ["lineage", "study"], roots,
                           lambda a, d: True, empty)
        return [dict(zip(EXPECT_COLUMNS, r)) for r in rows]

    def draw_rows(self, slug):
        out = []
        for e in self.expect_rows(slug):
            if e["role"] in ("none-added", "same-cohort"):
                continue
            sp, ab = self.name[slug], e["dataset"].split("__")[1].replace("_", " ")
            base = {"dataset": e["dataset"], "species": sp, "antibiotic": ab, "arm": e["arm"], "fit": e["fit"],
                    "split": e["split"], "slot": e["slot"], "role": e["role"], "constant": 0}
            for learner in LEARNERS:
                if e["split"] == "lineage" and e["dataset"] == self.refused:
                    out.append(dict(base, learner=learner, scored=0, why="lineage split unavailable"))
                else:
                    out.append(dict(base, learner=learner, scored=1, why="", **calls(e["arm"], learner)))
        return out

    def check_rows(self, slug):
        out = [{"row": e["dataset"], "check": "the plan and the prediction rows exist", "passed": 1, "detail": ""}
               for e in self.expect_rows(slug) if e["role"] == "draw"]
        for e in self.expect_rows(slug):
            if e["arm"] != "cap-all" or e["role"] in ("none-added", "same-cohort"):
                continue
            if e["split"] == "lineage" and e["dataset"] != self.refused:
                out.append({"row": f"{e['fit']}/{e['dataset']} lineage slot {e['slot']}",
                            "check": "the scored genomes of each outer fold are the default's", "passed": 1, "detail": ""})
            elif e["split"] == "study":
                out.append({"row": f"{e['fit']}/{e['dataset']} study slot {e['slot']}",
                            "check": "the test folds are the default's, with at least its training genomes",
                            "passed": 1, "detail": ""})
        return out

    def cohorts(self, slug):
        arms = {"cap-2000": {"genomes": 2000, "added": 0, "holds_default": 2000},
                "cap-4000": {"genomes": 4000, "added": 2000, "holds_default": 2000},
                "cap-all": {"genomes": 4096, "added": 2096, "holds_default": 2000}}
        for a in arms.values():
            a.update(lines=a["genomes"], sha256="d" * 64)
        return {"species": self.name[slug], "slug": slug, "ranked": 5000, "arms": arms}

    def write_species(self, slug, draws=None, checks=None):
        f = self.runs / slug
        draws = self.draw_rows(slug) if draws is None else draws
        tsv(f / "results" / "added-labelled.tsv", self.added_rows(slug), ADDED_COLUMNS)
        (f / "results" / "cohorts.json").write_text(json.dumps(self.cohorts(slug)), encoding="utf-8")
        tsv(f / "expect.tsv", self.expect_rows(slug), EXPECT_COLUMNS)
        commit, grm, code = self.binary[slug]
        (f / "binary.txt").write_text(f"commit {commit}\ngrm {grm}\ncode {code}\n", encoding="utf-8")
        (f / "score").mkdir(parents=True, exist_ok=True)
        C9.write_tsv(f / "score" / "draws.tsv", draws, score.DRAW_COLUMNS)
        C9.write_tsv(f / "score" / "datasets.tsv", C3.dataset_rows(draws), C3.DATASET_COLUMNS)
        tsv(f / "score" / "checks.tsv", self.check_rows(slug) if checks is None else checks,
            ["row", "check", "passed", "detail"])
        steps = [{"step": f"fit-cap-2000-{d}", "wall_s": "10.00", "cpu_s": "40.00", "peak_kb": "1000"}
                 for d in self.datasets(slug)]
        steps += [{"step": f"fit-cap-all-{d}", "wall_s": "20.00", "cpu_s": "80.00", "peak_kb": "2000"}
                  for d in self.datasets(slug) if (slug, d) not in self.none_added]
        tsv(f / "logs" / "timing.tsv", steps, ["step", "wall_s", "cpu_s", "peak_kb"])
        if slug == self.slugs[3]:
            (f / "logs" / f"folds-{self.refused}.log").write_text(
                f"{self.refused}: lineage split unavailable: 3 default slot plans refused past grm arm's replacements; "
                "no lineage run is produced\n", encoding="utf-8")

    def write(self, count_rows=None):
        tsv(self.lists, self.list_rows(), LIST_COLUMNS)
        for slug in self.slugs:
            self.write_species(slug)
        self.recount(count_rows)
        return self

    def all_rows(self):
        return [r for slug in self.slugs for r in read(self.runs / slug / "score" / "datasets.tsv")]

    def recount(self, rows=None):
        """count/ as score.py --count --withdrawn cap-4000 writes it, from every datasets.tsv or from `rows`."""
        rows = self.all_rows() if rows is None else rows
        counted = C3.count(rows, ["cap-4000"])
        self.count.mkdir(parents=True, exist_ok=True)
        C9.write_tsv(self.count / "count.tsv", counted, C3.COUNT_COLUMNS)
        C9.write_tsv(self.count / "conditions.tsv", score.conditions(counted, ["cap-4000"]), score.CONDITION_COLUMNS)

    def edit(self, path, fn):
        """Rewrites a TSV of the reading through fn(rows) -> rows, keeping its columns."""
        with open(path, encoding="utf-8") as f:
            columns = f.readline().rstrip("\n").split("\t")
        tsv(path, fn(read(path)), columns)


def sums(folder):
    """SHA256SUMS text over a species folder's results/, score/ (less score.log), binary.txt, expect.tsv,
    lists.txt and logs/timing.tsv, as run_species.sh writes it."""
    names = [p.relative_to(folder).as_posix() for sub in ("results", "score") for p in sorted((folder / sub).rglob("*"))
             if p.is_file() and p.suffix in (".json", ".tsv", ".txt")]
    names += [n for n in ("binary.txt", "expect.tsv", "lists.txt", "logs/timing.tsv") if (folder / n).is_file()]
    return "".join(f"{hashlib.sha256((folder / n).read_bytes()).hexdigest()}  {n}\n" for n in sorted(names))
