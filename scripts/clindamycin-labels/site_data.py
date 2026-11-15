#!/usr/bin/env python3
"""site_data.py REFERENCE [--check FILE]: the clindamycin kinds the web app bundles.

REFERENCE is reference/clindamycin-labels/rev14-20261115, the files behind
docs/clindamycin-labels.md. The script reads those files alone and writes, as
JSON, what the dataset report needs to show each of the six clindamycin
datasets' kind with the counts behind it:

- the release the check read and its evaluation revision (PLAN.md and
  results/count-checks.txt, which must agree), the plan's threshold for a
  change (PLAN.md), the class minimum (results/chosen.tsv), and the sha256 of
  each file read;
- for each dataset: its kind, the PMID the plan picks, the dataset's samples
  and the study's, the disagreements between GRM's known-gene rule and the
  reported phenotype in the dataset and in the study, the samples with no
  PMID, the reported resistant and susceptible samples left without the
  study, whether the dataset was scored without it, and where it was, the
  rows and the distinct tests that pass of those tested and whether the
  passing changes are rises or falls; and the released report's sha256 where
  results/logs/released-reports.sha256 holds it, null where it does not.

It holds no row per sample and no genome id.

The kind is not in any file: results/caveat.tsv gives the plan's verdict,
and the page adds the kind. The script derives it by this rule:

1. A dataset not scored without its study (chosen.tsv's `refit` is not
   "yes") is "one_study" when no reported resistant and no reported
   susceptible sample is left without the study, and "class_balance" when
   some are left but a class falls under the class minimum.
2. A dataset scored without its study is "measured_dependence" when at least
   one distinct test passes the plan's rule, and "no_dependence" when none
   does.

A distinct test is a split, a row and a measure of
results/scores/<dataset>/caveat.tsv, with the known-gene rule's rows counted
once, because they are the same on both splits; the script refuses a file
where they differ. The passing changes are "rises" when every one is above 0
and "falls" when every one is below; the script refuses a dataset whose
passing changes go both ways, since the site has no words for it. It also
refuses files that do not agree with each other: a dataset scored without its
study must leave both classes at the class minimum or more, and the caveat
table's study, existence part, dependence part, verdict and rows passing must
be what the other files give.

It tests nothing beyond that: the counts are the files'.

  python3 -I scripts/clindamycin-labels/site_data.py reference/clindamycin-labels/rev14-20261115 > web/src/data/clindamycin-labels.json
  python3 -I scripts/clindamycin-labels/site_data.py reference/clindamycin-labels/rev14-20261115 --check web/src/data/clindamycin-labels.json

With --check it prints nothing and exits 1 unless FILE is exactly the text
the script prints.
"""
import csv
import hashlib
import io
import json
import pathlib
import re
import sys

PAGE = "docs/clindamycin-labels.md"
KINDS = ("measured_dependence", "one_study", "class_balance", "no_dependence")
# The known-gene rule's row of the score tables, the same on both splits.
RULE_ROW = "catalog"
SPLITS = ("random", "lineage")

PLAN_RELEASE = re.compile(r"Release read: `([^`]+)` \(evaluation revision (\d+)\)\.")
PLAN_THRESHOLD = re.compile(r"is at least (0\.\d+) in absolute value and its 95% interval excludes 0")
INDEX_LINE = re.compile(r"^index: evaluation_revision=(\d+) datasets=(\d+) clindamycin datasets=(\d+)$", re.M)


def kind_of(scored, left_resistant, left_susceptible, distinct_passing):
    """The page's kind of a dataset, by the rule of the module's docstring."""
    if not scored:
        return "one_study" if left_resistant == 0 and left_susceptible == 0 else "class_balance"
    return "measured_dependence" if distinct_passing > 0 else "no_dependence"


class Reader:
    """Reads files of one reference directory and keeps the sha256 of each."""

    def __init__(self, root):
        self.root = pathlib.Path(root)
        self.sources = {}

    def text(self, rel):
        data = (self.root / rel).read_bytes()
        self.sources[rel] = hashlib.sha256(data).hexdigest()
        return data.decode("utf-8")

    def table(self, rel):
        return list(csv.DictReader(io.StringIO(self.text(rel), newline=""), delimiter="\t"))


def refuse(where, message):
    raise SystemExit(f"{where}: {message}")


def count(cell, where):
    if not re.fullmatch(r"\d+", cell):
        refuse(where, f"{cell!r} is not a count")
    return int(cell)


def tests_of(rows, dataset, where):
    """The distinct tests of one score table: (passes, change) for each, the
    known-gene rule's rows once. Also the rows passing and the rows tested."""
    distinct = {}
    for r in rows:
        if r["dataset"] != dataset:
            refuse(where, f"a row of {r['dataset']}")
        if r["split"] not in SPLITS or r["passes"] not in ("yes", "no"):
            refuse(where, f"row {r['split']} {r['row']} {r['measure']} has split {r['split']!r} and passes {r['passes']!r}")
        key = (r["row"], r["measure"]) if r["row"] == RULE_ROW else (r["split"], r["row"], r["measure"])
        value = (r["passes"] == "yes", float(r["total"]), r["cluster_low"], r["cluster_high"])
        if key in distinct and distinct[key] != value:
            refuse(where, f"the known-gene rule's {r['measure']} differs between the splits")
        if key in distinct and r["row"] != RULE_ROW:
            refuse(where, f"row {r['split']} {r['row']} {r['measure']} appears twice")
        distinct[key] = value
    passing_rows = sum(r["passes"] == "yes" for r in rows)
    passing = [change for passes, change, _, _ in distinct.values() if passes]
    return passing_rows, len(rows), len(passing), len(distinct), passing


def direction_of(changes, where):
    if not changes:
        return None
    if all(c > 0 for c in changes):
        return "rises"
    if all(c < 0 for c in changes):
        return "falls"
    return refuse(where, "the passing changes go both ways")


def build(reference):
    """The bundle's text for one reference directory."""
    read = Reader(reference)
    plan = read.text("PLAN.md")
    m = PLAN_RELEASE.search(plan)
    if m is None:
        refuse("PLAN.md", "names no release read")
    release, revision = m.group(1), int(m.group(2))
    t = PLAN_THRESHOLD.search(plan)
    if t is None:
        refuse("PLAN.md", "states no threshold for the dependence part")
    threshold = float(t.group(1))
    i = INDEX_LINE.search(read.text("results/count-checks.txt"))
    if i is None or int(i.group(1)) != revision:
        refuse("results/count-checks.txt", f"does not read an index of evaluation revision {revision}")

    chosen = read.table("results/chosen.tsv")
    caveats = {r["dataset"]: r for r in read.table("results/caveat.tsv")}
    by_study = read.table("results/by-study.tsv")
    hashes = {}
    for line in read.text("results/logs/released-reports.sha256").splitlines():
        digest, name = line.split("  ")
        if not (re.fullmatch(r"[0-9a-f]{64}", digest) and name.endswith(".json") and "/" not in name):
            refuse("results/logs/released-reports.sha256", f"{line!r} is not a report's sha256 line")
        hashes[name[: -len(".json")]] = digest

    if len(chosen) != int(i.group(3)) or sorted(caveats) != sorted(r["dataset"] for r in chosen):
        refuse("results/chosen.tsv", "does not hold the clindamycin datasets of the index and of results/caveat.tsv")
    if sorted(hashes) != sorted(set(hashes) & set(caveats)):
        refuse("results/logs/released-reports.sha256", "names a report outside the datasets")
    min_classes = {r["min_class"] for r in chosen}
    if len(min_classes) != 1:
        refuse("results/chosen.tsv", f"class minimums {sorted(min_classes)}")
    min_class = count(min_classes.pop(), "results/chosen.tsv")

    datasets = {}
    for r in sorted(chosen, key=lambda r: r["dataset"]):
        d = r["dataset"]
        where = f"results/chosen.tsv {d}"
        totals = [s for s in by_study if s["dataset"] == d and s["study"] == "TOTAL"]
        no_pmid = [s for s in by_study if s["dataset"] == d and s["study"] == "no-pmid"]
        if len(totals) != 1 or len(no_pmid) > 1:
            refuse(f"results/by-study.tsv {d}", "needs one TOTAL row and at most one no-pmid row")
        samples = count(r["samples_in_dataset"], where)
        if count(totals[0]["samples"], f"results/by-study.tsv {d}") != samples:
            refuse(f"results/by-study.tsv {d}", "its TOTAL is not the dataset's samples")
        left_resistant = count(r["reported_R_left"], where)
        left_susceptible = count(r["reported_S_left"], where)
        scored = r["refit"] == "yes"
        if scored != (min(left_resistant, left_susceptible) >= min_class):
            refuse(where, f"refit {r['refit']!r} with {left_resistant} resistant and {left_susceptible} susceptible left")

        c = caveats[d]
        if c["chosen_pmid"] != r["chosen_pmid"] or c["existence"] != ("no" if scored else "yes"):
            refuse(f"results/caveat.tsv {d}", "its study or existence part is not chosen.tsv's")
        if scored:
            rel = f"results/scores/{d}/caveat.tsv"
            rows_passing, rows_tested, distinct_passing, distinct_tested, changes = tests_of(read.table(rel), d, rel)
            direction = direction_of(changes, rel)
            if (c["dependence"], c["rows_passing"], c["rows_tested"]) != ("yes" if distinct_passing else "no", str(rows_passing), str(rows_tested)):
                refuse(f"results/caveat.tsv {d}", f"its dependence part is not what {rel} gives")
        else:
            rows_passing = rows_tested = distinct_passing = distinct_tested = direction = None
            if (c["dependence"], c["rows_passing"], c["rows_tested"]) != ("NA", "NA", "NA"):
                refuse(f"results/caveat.tsv {d}", "a dataset not scored has a dependence part")
        if c["caveat"] != ("yes" if (not scored or distinct_passing) else "no"):
            refuse(f"results/caveat.tsv {d}", f"its verdict {c['caveat']!r} is not the rule's")

        datasets[d] = {
            "kind": kind_of(scored, left_resistant, left_susceptible, distinct_passing),
            "pmid": r["chosen_pmid"],
            "samples": samples,
            "study_samples": count(r["samples_in_study"], where),
            "disagreements": count(r["disagreements_in_dataset"], where),
            "study_disagreements": count(r["disagreements_in_study"], where),
            "no_pmid_samples": count(no_pmid[0]["samples"], f"results/by-study.tsv {d}") if no_pmid else 0,
            "left_resistant": left_resistant,
            "left_susceptible": left_susceptible,
            "scored_without": scored,
            "rows_passing": rows_passing,
            "rows_tested": rows_tested,
            "distinct_passing": distinct_passing,
            "distinct_tested": distinct_tested,
            "direction": direction,
            "report_sha256": hashes.get(d),
        }

    data = {
        "page": PAGE,
        "reference": f"reference/clindamycin-labels/{pathlib.Path(reference).resolve().name}",
        "release": release,
        "evaluation_revision": revision,
        "min_class": min_class,
        "threshold": threshold,
        "sources": dict(sorted(read.sources.items())),
        "datasets": datasets,
    }
    return json.dumps(data, separators=(",", ":")) + "\n"


def main(argv):
    if len(argv) == 2:
        sys.stdout.write(build(argv[1]))
        return 0
    if len(argv) == 4 and argv[2] == "--check":
        return 0 if pathlib.Path(argv[3]).read_text(encoding="utf-8") == build(argv[1]) else 1
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
