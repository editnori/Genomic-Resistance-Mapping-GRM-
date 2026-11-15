#!/usr/bin/env python3
"""site_data.py REFERENCE [--check FILE]: the clindamycin kinds the web app bundles.

REFERENCE is reference/clindamycin-labels/rev14-20261115, the files behind
docs/clindamycin-labels.md. The script reads those files and one file beside
that directory, the release's checksum line of each report
(reference/lineage-control/<release>.reports.sha256), and writes, as JSON,
what the dataset report needs to show each of the six clindamycin datasets'
kind with the counts behind it:

- the release the check read and its evaluation revision (PLAN.md and
  results/count-checks.txt, which must agree), the plan's threshold for a
  change (PLAN.md), the class minimum (results/chosen.tsv), and the sha256 of
  each file read;
- for each dataset: its kind, the PMID the plan picks, the dataset's samples
  and the study's, the disagreements between GRM's known-gene rule and the
  reported phenotype in the dataset and in the study, the samples with no
  PMID and their disagreements, the samples and disagreements under other
  PMIDs (results/by-study.tsv), the reported resistant and susceptible
  samples left without the study, whether the dataset was scored without
  it, and where it was, the rows and the distinct tests that pass of those
  tested and whether the passing changes are rises or falls; and the
  released report's sha256, from the release's checksum lines. The script
  refuses a line of results/logs/released-reports.sha256 that differs from
  the release's.

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
where they differ, and a row that appears twice within one split. A row
passes when its change is at least the plan's threshold in absolute value
and its interval (cluster_low to cluster_high) excludes 0; the script
recomputes this for every row and refuses a `passes` cell that differs. The
passing changes are "rises" when every one is above 0 and "falls" when every
one is below; the script refuses a dataset whose passing changes go both
ways, since the site has no words for it. It refuses a dataset where only
the known-gene rule's rows pass, since the page's words for measured
dependence are about the learners' figures. It also refuses files that do
not agree with each other: a dataset scored without its study must leave
both classes at the class minimum or more, the caveat table's study,
existence part, dependence part, verdict and rows passing must be what the
other files give, and the study, the samples with no PMID and the other
PMIDs of results/by-study.tsv must add up to its TOTAL row.

The dataset report says that the study's share of the disagreements stands
out only against the samples with no PMID when all three of these hold, in
the counts of the file (no_pmid_contrast below):

1. most samples carry no PMID: 2 * no_pmid_samples > samples;
2. the study disagrees more often than the samples with no PMID:
   study_disagreements * no_pmid_samples > no_pmid_disagreements * study_samples;
3. among the samples that carry a PMID, the study holds at most one
   disagreement more than its share of those samples gives (the page's
   "close or reversed"): (study_disagreements - 1) * (study_samples +
   other_pmid_samples) <= (study_disagreements + other_pmid_disagreements) *
   study_samples.

For release rev14-274-20261112 the rule holds for
campylobacter_coli__clindamycin and campylobacter_jejuni__clindamycin and
for no other dataset; tests/clindamycin_labels checks that.

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


def no_pmid_contrast(d):
    """Whether the study's share of the disagreements stands out only against
    the samples with no PMID, by the rule of the module's docstring; `d` is a
    dataset entry of the bundle."""
    pmid_samples = d["study_samples"] + d["other_pmid_samples"]
    pmid_disagreements = d["study_disagreements"] + d["other_pmid_disagreements"]
    return (
        2 * d["no_pmid_samples"] > d["samples"]
        and d["study_disagreements"] * d["no_pmid_samples"] > d["no_pmid_disagreements"] * d["study_samples"]
        and (d["study_disagreements"] - 1) * pmid_samples <= pmid_disagreements * d["study_samples"]
    )


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


def passes_rule(change, low, high, threshold):
    """The plan's rule for one row: the change is at least the threshold in
    absolute value and its interval excludes 0."""
    return abs(change) >= threshold and (low > 0 or high < 0)


def tests_of(rows, dataset, where, threshold):
    """The distinct tests of one score table: (passes, change) for each, the
    known-gene rule's rows once. Also the rows passing and the rows tested,
    and whether a row other than the known-gene rule's passes."""
    distinct = {}
    seen = set()
    for r in rows:
        if r["dataset"] != dataset:
            refuse(where, f"a row of {r['dataset']}")
        if r["split"] not in SPLITS or r["passes"] not in ("yes", "no"):
            refuse(where, f"row {r['split']} {r['row']} {r['measure']} has split {r['split']!r} and passes {r['passes']!r}")
        if (r["split"], r["row"], r["measure"]) in seen:
            refuse(where, f"row {r['split']} {r['row']} {r['measure']} appears twice")
        seen.add((r["split"], r["row"], r["measure"]))
        if (r["passes"] == "yes") != passes_rule(float(r["total"]), float(r["cluster_low"]), float(r["cluster_high"]), threshold):
            refuse(where, f"row {r['split']} {r['row']} {r['measure']} has passes {r['passes']!r}, and the plan's rule gives the other")
        key = (r["row"], r["measure"]) if r["row"] == RULE_ROW else (r["split"], r["row"], r["measure"])
        value = (r["passes"] == "yes", float(r["total"]), r["cluster_low"], r["cluster_high"])
        if key in distinct and distinct[key] != value:
            refuse(where, f"the known-gene rule's {r['measure']} differs between the splits")
        distinct[key] = value
    passing_rows = sum(r["passes"] == "yes" for r in rows)
    passing = [change for passes, change, _, _ in distinct.values() if passes]
    learner_passes = any(r["passes"] == "yes" and r["row"] != RULE_ROW for r in rows)
    return passing_rows, len(rows), len(passing), len(distinct), passing, learner_passes


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
    # The release's own checksum line of each report, in the directory beside this one.
    release_rel = f"../../lineage-control/{release}.reports.sha256"
    released = {}
    for line in read.text(release_rel).splitlines():
        m = re.fullmatch(r"([0-9a-f]{64})  \./datasets/([a-z0-9_]+)\.json", line)
        if m is None:
            refuse(release_rel, f"{line!r} is not a report's sha256 line")
        released[m.group(2)] = m.group(1)
    for d, digest in hashes.items():
        if released.get(d) != digest:
            refuse("results/logs/released-reports.sha256", f"the sha256 of {d} is not its line in {release_rel}")

    if len(chosen) != int(i.group(3)) or sorted(caveats) != sorted(r["dataset"] for r in chosen):
        refuse("results/chosen.tsv", "does not hold the clindamycin datasets of the index and of results/caveat.tsv")
    if sorted(hashes) != sorted(set(hashes) & set(caveats)):
        refuse("results/logs/released-reports.sha256", "names a report outside the datasets")
    min_classes = {r["min_class"] for r in chosen}
    if len(min_classes) != 1:
        refuse("results/chosen.tsv", f"class minimums {sorted(min_classes)}")
    min_class = count(min_classes.pop(), "results/chosen.tsv")
    missing = sorted(r["dataset"] for r in chosen if r["dataset"] not in released)
    if missing:
        refuse(release_rel, f"has no line for {', '.join(missing)}")

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
        study_rows = [s for s in by_study if s["dataset"] == d and s["study"] == r["chosen_pmid"]]
        others = [s for s in by_study if s["dataset"] == d and s["study"] not in ("TOTAL", "no-pmid", r["chosen_pmid"])]
        if len(study_rows) != 1 or any(not re.fullmatch(r"\d+(\+\d+)*", s["study"]) for s in others):
            refuse(f"results/by-study.tsv {d}", "needs one row of the study the plan picks, and one or more PMIDs in every other row")
        parts = {}
        for name, group in (("study", study_rows), ("no_pmid", no_pmid), ("other", others)):
            parts[name] = tuple(sum(count(s[k], f"results/by-study.tsv {d}") for s in group) for k in ("samples", "disagreements"))
        if tuple(map(sum, zip(*parts.values()))) != (samples, count(totals[0]["disagreements"], f"results/by-study.tsv {d}")):
            refuse(f"results/by-study.tsv {d}", "the study, the samples with no PMID and the other PMIDs do not add up to its TOTAL")
        if parts["study"] != (count(r["samples_in_study"], where), count(r["disagreements_in_study"], where)) or parts["study"][1] + parts["no_pmid"][1] + parts["other"][1] != count(r["disagreements_in_dataset"], where):
            refuse(f"results/by-study.tsv {d}", "its study or its disagreements are not chosen.tsv's")
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
            rows_passing, rows_tested, distinct_passing, distinct_tested, changes, learner_passes = tests_of(read.table(rel), d, rel, threshold)
            direction = direction_of(changes, rel)
            if distinct_passing and not learner_passes:
                refuse(rel, "only the known-gene rule's rows pass, and the page's words for measured dependence are about the learners' figures")
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
            "no_pmid_disagreements": parts["no_pmid"][1],
            "other_pmid_samples": parts["other"][0],
            "other_pmid_disagreements": parts["other"][1],
            "left_resistant": left_resistant,
            "left_susceptible": left_susceptible,
            "scored_without": scored,
            "rows_passing": rows_passing,
            "rows_tested": rows_tested,
            "distinct_passing": distinct_passing,
            "distinct_tested": distinct_tested,
            "direction": direction,
            "report_sha256": released[d],
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
