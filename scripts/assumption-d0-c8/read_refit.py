"""Reads D0 and C8 of the refit over its species into the tables of a public page.

    python3 -I read_refit.py --out DIR --code-sha256 SHA \
        --species N --datasets N --d0-datasets N --d0-species N --d0-pairs N \
        --marks MARKS.tsv --c8-one-fold zero|undefined \
        --c8-over-five pairs|dataset-median --c8-bands none|LOW,HIGH \
        --order-rest first-only|default-after SPECIES_DIR...

Each SPECIES_DIR is one species' copy of scripts/refit/run_species.sh's OUT:
DONE, SHA256SUMS, binary.txt, expect.txt, d0.txt, check/ and results/
(docs/assumption-refit.md, "Copy"). The reading's rules are fixed by the
lead before any species' count is read; the rules the plan leaves open are
the arguments above, with no default.

It stops, prints `stop: <reason>` and writes nothing unless:

  1. every directory holds DONE in the run script's form and no FAILED,
     with the code sha256 given and the binary.txt it names;
  2. every file SHA256SUMS lists has its sha256, SHA256SUMS has the one
     DONE names, and every file read here is listed;
  3. every row of the check tables named by --equal-tables reads equal
     (by default the two tables of scripts/refit/refit_check.py);
  4. the species number --species, none twice, and the datasets number
     --datasets, none twice, each with its five learners in C8; with
     --marks they are exactly the datasets that are not M, species by
     species;
  5. D0's datasets are the plan's MISS set among each species' listed
     datasets (scripts/assumption-d0-c8/miss_sets.py), number --d0-datasets
     in --d0-species species, and hold --d0-pairs miss pairs (dataset and
     learner), each read from fold models; `d0.py combine` adds them up,
     with its own gates on the dataset and pair counts.

Then it writes the page's tables: the species, D0 per learner (datasets in
each cause, each cause's share of the gap to the AMRFinderPlus program,
ties, pairs with no missed fold), D0 per dataset, the ordering line as the
corrections page's entry "D0 causes" reads it, and C8 per learner and over
the five (`c8.py combine`). Species and antibiotics are named in full. D0's
sets are chosen by their result: its tables order comparisons and test
nothing. Standard library only.
"""

import argparse
import contextlib
import csv
import hashlib
import io
import os
import re
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c8  # noqa: E402
import d0  # noqa: E402
import foldread as fr  # noqa: E402
import miss_sets  # noqa: E402

SHA_RE = re.compile(r"[0-9a-f]{64}\Z")
DONE_RE = re.compile(
    r"DONE (?P<species>[a-z][a-z0-9_]*) (?P<time>\S+) binary (?P<binary>[0-9a-f]{64}) "
    r"code (?P<code>[0-9a-f]{64}) (?:d0 (?P<d0>\d+) datasets|d0 skipped: no MISS dataset among the listed) "
    r"SHA256SUMS (?P<sums>[0-9a-f]{64})\Z")
SUMS_RE = re.compile(r"(?P<sha>[0-9a-f]{64})  (?P<path>[A-Za-z0-9_.\-/]+)\Z")
CHECK_DATASETS = ["dataset", "fitted_fields", "fields_that_differ", "samples_with_other_descriptive_fields",
                  "fold_models", "samples_with_rule_groups", "unexpected_new_keys", "text_equal_less_allowed"]
# The check's later column (the correction "Refit cards"): feature cards whose `kmer` or `locus` differs from the
# release's. A table without it comes from the earlier check, which let no card differ, so its count is 0.
CHECK_CARDS = "cards_with_other_locus"
CHECK_LEARNERS = ["dataset", "split", "learner", "verdict", "fields_that_differ", "largest_abs_difference"]
SPECS = ["dataset", "breakpoint", "root_r_above", "release_r_above", "root_n", "release_n"]
FEATURE_LISTS = ["fetched", "equal_to_master"]
TABLES = ("datasets.tsv", "learners.tsv", "release-specs.tsv", "feature-lists.tsv")
NEED_EQUAL = ("datasets.tsv", "learners.tsv")
D0_FIELDS = {"folds_offered", "folds_inside_floor"}
FIRST = {"absent": "C3 first", "removed": "C4's floor arm first", "lineage": "C4's adjusted arms first",
         "unused": "C5 first"}
AFTER = {"absent": "C3, then C4 and C5", "removed": "C4 with its floor arm first, then C3 and C5",
         "lineage": "C4 with its adjusted arms first, then C3 and C5", "unused": "C5, then C4 and C3"}
DEFAULT_ORDER = "C4, C3, C5 (the plan's default)"
D0_NOTE = "D0's sets are chosen by their result: it orders comparisons and tests nothing"
POOLED_SHARE = "-"
POOLED_NOTE = "the gap share pooled over the five is not on the page (ruling 10); it is in d0/d0-summary.tsv"
OUTPUTS = ["gates.txt", "summary.txt", "page-species.tsv", "page-d0-learners.tsv", "page-d0-order.tsv",
           "page-d0-datasets.tsv", "page-c8-summary.tsv", "page-c8-datasets.tsv", "d0", "c8"]


class Stop(Exception):
    """An input the reading does not read on from."""


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_lines(path):
    with open(path, encoding="utf-8") as f:
        text = f.read()
    if text and not text.endswith("\n"):
        raise Stop(f"{path}: the last line has no newline")
    return text.splitlines()


def tsv(path, header, *more):
    """The table's rows by column name. Its header must be `header` or one of `more`."""
    with open(path, encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE))
    if not rows or rows[0] not in (header,) + more:
        raise Stop(f"{path}: its columns are not the run's")
    header = rows[0]
    for i, r in enumerate(rows[1:], start=2):
        if len(r) != len(header):
            raise Stop(f"{path}: line {i} has {len(r)} cells, want {len(header)}")
    return [dict(zip(header, r)) for r in rows[1:]]


def dataset_list(path, species):
    """expect.txt or d0.txt: dataset ids of the species, sorted, once each."""
    out = read_lines(path)
    for d in out:
        if not c8.DATASET_RE.match(d) or not d.startswith(species + "__"):
            raise Stop(f"{path}: a line is not a dataset of {species}")
    if out != sorted(set(out)):
        raise Stop(f"{path}: not sorted, or a dataset twice")
    return out


# ---------------------------------------------------------------- one species

def read_done(path, code_sha):
    """Gate 1: DONE, no FAILED, the pinned code sha256."""
    if os.path.exists(os.path.join(path, "FAILED")):
        first = (read_lines(os.path.join(path, "FAILED")) or [""])[0]
        raise Stop(f"{path}: FAILED is present ({first})")
    done_path = os.path.join(path, "DONE")
    if not os.path.isfile(done_path):
        raise Stop(f"{path}: no DONE")
    lines = read_lines(done_path)
    m = DONE_RE.match(lines[0]) if len(lines) == 1 else None
    if not m:
        raise Stop(f"{done_path}: not one line in the run script's form")
    if m["code"] != code_sha:
        raise Stop(f"{m['species']}: DONE has code sha256 {m['code']}, the pinned one is {code_sha}")
    return m


def verify_sums(path, done):
    """Gate 2: SHA256SUMS is DONE's, and every file it lists has its sha256."""
    sums_path = os.path.join(path, "SHA256SUMS")
    if not os.path.isfile(sums_path) or sha256(sums_path) != done["sums"]:
        raise Stop(f"{path}: SHA256SUMS is missing or is not the one DONE names")
    listed = {}
    for line in read_lines(sums_path):
        m = SUMS_RE.match(line)
        if not m or m["path"].startswith("/") or ".." in m["path"].split("/") or m["path"] in listed:
            raise Stop(f"{sums_path}: a line is not '<sha256>  <relative path>', or a path is twice")
        listed[m["path"]] = m["sha"]
    for rel, want in sorted(listed.items()):
        p = os.path.join(path, rel)
        if not os.path.isfile(p) or sha256(p) != want:
            raise Stop(f"{path}: {rel} is missing or differs from its sha256 in SHA256SUMS")
    return listed


def read_check(path, listed, expect, equal_tables):
    """Gate 3: the check tables named in equal_tables read equal on every
    row. All present tables are read and counted."""
    out = {}
    rows = tsv(os.path.join(path, "check", "datasets.tsv"), CHECK_DATASETS, CHECK_DATASETS + [CHECK_CARDS])
    got = [r["dataset"] for r in rows]
    if sorted(got) != expect or len(set(got)) != len(got):
        raise Stop(f"{path}: check/datasets.tsv does not hold each listed dataset once")
    for r in rows:
        if not r["samples_with_other_descriptive_fields"].isdigit() or not r.get(CHECK_CARDS, "0").isdigit():
            raise Stop(f"{path}: check/datasets.tsv: {r['dataset']}: a count is not a count")
        if "datasets.tsv" in equal_tables and (
                r["fitted_fields"] != "equal" or r["fields_that_differ"] != "-" or r["fold_models"] != "complete"
                or r["unexpected_new_keys"] != "-" or r["text_equal_less_allowed"] != "yes"):
            raise Stop(f"{path}: check/datasets.tsv: {r['dataset']} does not read equal")
    out["check_datasets"] = len(rows)
    out["descriptive"] = sum(int(r["samples_with_other_descriptive_fields"]) for r in rows)
    cards = [int(r.get(CHECK_CARDS, "0")) for r in rows]
    out["cards"], out["cards_datasets"] = sum(cards), sum(n > 0 for n in cards)
    rows = tsv(os.path.join(path, "check", "learners.tsv"), CHECK_LEARNERS)
    keys = [(r["dataset"], r["split"], r["learner"]) for r in rows]
    if len(set(keys)) != len(keys) or any(r["dataset"] not in expect for r in rows):
        raise Stop(f"{path}: check/learners.tsv holds a row twice or a dataset not listed")
    for d in expect:
        splits = {s for x, s, _ in keys if x == d}
        if not splits:
            raise Stop(f"{path}: check/learners.tsv has no row of {d}")
        for s in splits:
            if not set(fr.LEARNERS) <= {l for x, y, l in keys if x == d and y == s}:
                raise Stop(f"{path}: check/learners.tsv: {d} {s} lacks one of the five learners")
    for r in rows:
        if "learners.tsv" in equal_tables and (
                r["verdict"] != "equal" or r["fields_that_differ"] != "0" or r["largest_abs_difference"] != "0"):
            raise Stop(f"{path}: check/learners.tsv: {r['dataset']} {r['split']} {r['learner']} does not read equal")
    out["check_learners"] = len(rows)
    sp = os.path.join(path, "check", "release-specs.tsv")
    out["specs_differ"] = None
    if "check/release-specs.tsv" in listed:
        rows = tsv(sp, SPECS)
        if sorted(r["dataset"] for r in rows) != expect:
            raise Stop(f"{path}: check/release-specs.tsv does not hold each listed dataset once")
        if any(r["breakpoint"] not in ("equal", "differs") for r in rows):
            raise Stop(f"{path}: check/release-specs.tsv: a breakpoint cell is neither equal nor differs")
        out["specs_differ"] = sum(r["breakpoint"] == "differs" for r in rows)
        if "release-specs.tsv" in equal_tables and out["specs_differ"]:
            raise Stop(f"{path}: check/release-specs.tsv: {out['specs_differ']} breakpoints differ")
    elif "release-specs.tsv" in equal_tables:
        raise Stop(f"{path}: check/release-specs.tsv is not listed in SHA256SUMS")
    out["features_fetched"] = out["features_equal"] = None
    if "check/feature-lists.tsv" in listed:
        rows = tsv(os.path.join(path, "check", "feature-lists.tsv"), FEATURE_LISTS)
        if len(rows) != 1 or not all(v.isdigit() for v in rows[0].values()):
            raise Stop(f"{path}: check/feature-lists.tsv is not one row of two counts")
        out["features_fetched"], out["features_equal"] = int(rows[0]["fetched"]), int(rows[0]["equal_to_master"])
        if "feature-lists.tsv" in equal_tables and out["features_fetched"] != out["features_equal"]:
            raise Stop(f"{path}: check/feature-lists.tsv: {out['features_equal']} of {out['features_fetched']} "
                       "feature lists equal the master list")
    elif "feature-lists.tsv" in equal_tables:
        raise Stop(f"{path}: check/feature-lists.tsv is not listed in SHA256SUMS")
    return out


def read_d0_table(path, species, d0_list):
    """One species' results/d0/d0-datasets.tsv: the D0 datasets exactly, each
    with its five learners, read from fold models, no cause text but the
    plan's and the entry's."""
    p = os.path.join(path, "results", "d0", "d0-datasets.tsv")
    rows = tsv(p, d0.DATA_HEADER)
    seen = set()
    causes = re.compile(r"(?:%s)(?:\+(?:%s))*\Z" % ("|".join(d0.CAUSES), "|".join(d0.CAUSES)))
    for i, r in enumerate(rows, start=2):
        where = f"{p}: line {i}"
        if r["dataset"] not in d0_list or r["species"] != species or r["learner"] not in fr.LEARNERS:
            raise Stop(f"{where}: not a D0 dataset of {species} with one of the five learners")
        if (r["dataset"], r["learner"]) in seen:
            raise Stop(f"{where}: {r['dataset']} {r['learner']} twice")
        seen.add((r["dataset"], r["learner"]))
        if r["source"] != "fold models" or r["reading"] != d0.READING:
            raise Stop(f"{where}: not read from fold models, or another reading")
        if any(r[k] != "NA" for k in D0_FIELDS):
            raise Stop(f"{where}: holds a cause read without fold models")
        if r["in_miss"] not in ("yes", "no") or (r["cause"] != "NA" and not causes.match(r["cause"])):
            raise Stop(f"{where}: a mark or a cause that is not D0's")
        if r["cause"] != "NA":
            parts = r["cause"].split("+")
            if parts != [c for c in d0.CAUSES if c in parts] or r["in_miss"] != "yes":
                raise Stop(f"{where}: causes out of D0's order, or a cause outside the miss set")
        for k in ("program_ba", "learner_lineage_ba", "gap"):
            try:
                float(r[k])
            except ValueError:
                raise Stop(f"{where}: {k} is not a number") from None
        for k in ["missed_folds"] + [f"folds_{c}" for c in d0.CAUSES]:
            if not r[k].isdigit():
                raise Stop(f"{where}: {k} is not a count")
    if seen != {(d, l) for d in d0_list for l in fr.LEARNERS}:
        raise Stop(f"{p}: does not hold the five learners of each D0 dataset")
    return rows


def read_species(path, code_sha, equal_tables):
    done = read_done(path, code_sha)
    listed = verify_sums(path, done)
    species = done["species"]
    need = ["binary.txt", "expect.txt", "d0.txt", "check/datasets.tsv", "check/learners.tsv",
            "results/c8/c8-datasets.tsv"]
    for rel in need:
        if rel not in listed:
            raise Stop(f"{path}: {rel} is not listed in SHA256SUMS")
    binary = dict(line.split(" ", 1) for line in read_lines(os.path.join(path, "binary.txt")) if " " in line)
    if binary.get("code") != done["code"] or binary.get("grm") != done["binary"]:
        raise Stop(f"{path}: binary.txt does not name DONE's binary and code")
    expect = dataset_list(os.path.join(path, "expect.txt"), species)
    if not expect:
        raise Stop(f"{path}: expect.txt lists no dataset")
    d0_list = dataset_list(os.path.join(path, "d0.txt"), species)
    if not set(d0_list) <= set(expect):
        raise Stop(f"{path}: d0.txt lists a dataset expect.txt does not")
    if (done["d0"] is None) != (not d0_list) or (done["d0"] is not None and int(done["d0"]) != len(d0_list)):
        raise Stop(f"{path}: DONE's D0 count is not d0.txt's")
    check = read_check(path, listed, expect, equal_tables)
    c8_rows = c8.read_table(os.path.join(path, "results", "c8", "c8-datasets.tsv"))
    if {r["species"] for r in c8_rows} != {species}:
        raise Stop(f"{path}: the C8 table holds another species")
    if sorted({r["dataset"] for r in c8_rows}) != expect:
        raise Stop(f"{path}: the C8 table does not hold the listed datasets exactly")
    if any(r["m_mark"] is not False for r in c8_rows):
        raise Stop(f"{path}: the C8 table holds a dataset marked M or unmarked")
    d0_rows = []
    if d0_list:
        if "results/d0/d0-datasets.tsv" not in listed:
            raise Stop(f"{path}: results/d0/d0-datasets.tsv is not listed in SHA256SUMS")
        d0_rows = read_d0_table(path, species, d0_list)
    elif any(rel.startswith("results/d0/") for rel in listed):
        raise Stop(f"{path}: D0 was skipped and results/d0 holds tables")
    return {"path": path, "species": species, "expect": expect, "d0": d0_list, "check": check,
            "c8": c8_rows, "d0_rows": d0_rows, "time": done["time"]}


# ---------------------------------------------------------------- over species

def gates(spp, a, names):
    """Gates 4 and 5 over species. Returns lines for gates.txt."""
    lines = []
    slugs = [s["species"] for s in spp]
    for sp in sorted(set(slugs)):
        if slugs.count(sp) > 1:
            raise Stop(f"species {sp} twice: {', '.join(s['path'] for s in spp if s['species'] == sp)}")
    if len(spp) != a.species:
        raise Stop(f"{len(spp)} species, want {a.species}")
    datasets = [d for s in spp for d in s["expect"]]
    if len(datasets) != len(set(datasets)):
        raise Stop("a dataset twice over species")
    if len(datasets) != a.datasets:
        raise Stop(f"{len(datasets)} datasets, want {a.datasets}")
    for d in datasets:
        if d not in names:
            raise Stop(f"{d}: no species and antibiotic name in the names table")
    lines.append(f"species: {len(spp)} of {a.species}, each with DONE and code sha256 {a.code_sha256}")
    lines.append(f"datasets: {len(datasets)} of {a.datasets}, none twice; check tables read equal: "
                 f"{', '.join(sorted(a.equal_tables))}")
    if a.marks:
        marks = fr.read_marks(a.marks)
        notm = sorted(d for d, m in marks.items() if not m)
        missing = sorted(set(notm) - set(datasets))
        extra = sorted(set(datasets) - set(notm))
        if missing or extra:
            raise Stop(f"the datasets are not those that are not M: {len(missing)} missing"
                       f"{' (the first ' + missing[0] + ')' if missing else ''}, {len(extra)} not expected"
                       f"{' (the first ' + extra[0] + ')' if extra else ''}")
        lines.append(f"marks: the {len(datasets)} datasets are exactly those of {a.marks} that are not M")
    _, miss, _ = miss_sets.sets(a.record)
    union = {d for v in miss.values() for d in v}
    for s in spp:
        want = sorted(union & set(s["expect"]))
        if s["d0"] != want:
            raise Stop(f"{s['species']}: d0.txt lists {len(s['d0'])} datasets, the MISS set among its listed "
                       f"holds {len(want)}")
        for r in s["d0_rows"]:
            if (r["in_miss"] == "yes") != (r["dataset"] in miss[r["learner"]]):
                raise Stop(f"{s['species']}: {r['dataset']} {r['learner']}: in the miss set by D0's table and "
                           "not by the record, or the reverse")
    a.plan_datasets, a.plan_pairs = len(union), sum(len(v) for v in miss.values())
    d0_species = [s for s in spp if s["d0"]]
    n_d0 = sum(len(s["d0"]) for s in d0_species)
    if n_d0 != a.d0_datasets or len(d0_species) != a.d0_species:
        raise Stop(f"D0: {n_d0} datasets in {len(d0_species)} species, want {a.d0_datasets} in {a.d0_species}")
    lines.append(f"D0: {n_d0} datasets in {len(d0_species)} species, each species' list the MISS set among its "
                 f"listed datasets")
    return lines


def run_d0(spp, a, tmp):
    """`d0.py combine` over the species' D0 tables, with its own gates."""
    dirs = [os.path.join(s["path"], "results", "d0") for s in spp if s["d0"]]
    out = os.path.join(tmp, "d0")
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            d0.main(["combine", "--out", out, "--expect-datasets", str(a.d0_datasets),
                     "--expect-pairs", str(a.d0_pairs)] + dirs)
    except SystemExit as e:
        raise Stop(f"d0.py combine: {e}") from None
    summary = tsv(os.path.join(out, "d0-summary.tsv"), d0.SUMMARY_HEADER)
    line = read_lines(os.path.join(out, "d0-ordering.txt"))[0]
    every = {r["learner"]: r for r in summary if r["subset"] == "all"}
    for r in every.values():
        if any(r[f"datasets_{c}"] != "0" for c in d0.OPEN_ONLY):
            raise Stop("d0.py combine counted a cause read without fold models")
    if a.d0_learner_pairs:
        for item in a.d0_learner_pairs.split(","):
            learner, _, n = item.partition("=")
            if learner not in every or not n.isdigit():
                raise Stop(f"--d0-learner-pairs: {item} is not learner=count")
            if int(every[learner]["miss_datasets"]) != int(n):
                raise Stop(f"D0: {learner} has {every[learner]['miss_datasets']} miss pairs, want {n}")
    return out, every, line


def setter(line):
    """The cause `d0.ordering` says sets the order, or None."""
    m = re.match(r"(absent|removed|lineage|unused) holds ", line)
    return m.group(1) if m else None


def stability(counts, best, plan_pairs, pairs):
    """Ruling 6: the order is set on the miss pairs read, `pairs` of the
    plan's `plan_pairs`; m = plan_pairs - pairs are the MISS datasets' pairs
    that are M and not read. With x the leading cause's pairs and y the next
    cause's among the four, a named first block is stable when x / plan_pairs
    is 0.40 or more and x is above y + m; the default order is stable when no
    cause of the four reaches 0.40 of plan_pairs with m pairs added.
    `counts` holds the four causes' pairs. Returns (stable, words)."""
    m = plan_pairs - pairs
    if m < 0:
        raise Stop(f"D0: {pairs} miss pairs read, more than the plan's {plan_pairs}")
    ranked = sorted(((counts[c], c) for c in d0.FOUR), reverse=True)
    if best is not None:
        x = counts[best]
        y = max(counts[c] for c in d0.FOUR if c != best)
        stable = 5 * x >= 2 * plan_pairs and x > y + m
        words = (f"{best} holds {x} of the plan's {plan_pairs} miss pairs ({x / plan_pairs:.3f}) and the next cause "
                 f"of the four {y}; ")
        if stable:
            return True, words + f"{x} is 0.40 or more of {plan_pairs} and above {y} + {m}: the first block is stable"
        return False, words + f"the {m} pairs of the MISS datasets that are M could change the first block"
    most = ranked[0][0]
    stable = all(5 * (n + m) < 2 * plan_pairs for n, _ in ranked)
    words = f"the largest cause of the four holds {most}, {most + m} with the {m} pairs added, of {plan_pairs}; "
    if stable:
        return True, words + f"no cause of the four reaches 0.40 of {plan_pairs} with {m} pairs added: " \
                             "the default order is stable"
    return False, words + f"the {m} pairs of the MISS datasets that are M could change the order"


def page_tables(spp, a, names, every, line, c8_long, c8_summary):
    by_d0 = {}
    for s in spp:
        for r in s["d0_rows"]:
            by_d0[(r["dataset"], r["learner"])] = r
    species_rows = []
    for s in sorted(spp, key=lambda s: s["species"]):
        abx = [names[d][1] for d in s["expect"]]
        d0abx = [names[d][1] for d in s["d0"]]
        c = s["check"]
        species_rows.append([c8.species_name(s["species"]), len(s["expect"]), "; ".join(abx), len(s["d0"]),
                             "; ".join(d0abx) if d0abx else "none", c["check_datasets"], c["check_learners"],
                             c["specs_differ"], c["features_fetched"], c["features_equal"], c["descriptive"],
                             c["cards_datasets"], c["cards"]])
    species_header = ["species", "datasets", "antibiotics", "d0_datasets", "d0_antibiotics", "check_dataset_rows",
                      "check_learner_rows", "breakpoints_set_to_the_release", "feature_lists_fetched",
                      "feature_lists_equal_to_master", "samples_with_other_descriptive_fields",
                      "datasets_with_other_locus_cards", "cards_with_other_locus"]
    learners_header = ["learner", "miss_pairs"] + [f"datasets_{c}" for c in d0.CAUSES] + \
        [f"share_of_gap_to_the_amrfinderplus_program_{c}" for c in d0.CAUSES] + \
        ["ties", "pairs_with_no_missed_fold", "reading"]
    learners = []
    for l in list(fr.LEARNERS) + ["all"]:
        r = every[l]
        if l == "all":
            # Ruling 10: the plan's share is "of the learner's gap", so the share pooled over the five is not on the
            # page. It stays in the reader's file, d0/d0-summary.tsv.
            shares, note = [POOLED_SHARE] * len(d0.CAUSES), D0_NOTE + "; " + POOLED_NOTE
        else:
            shares = [None if r[f"gap_share_{c}"] == "NA" else float(r[f"gap_share_{c}"]) for c in d0.CAUSES]
            note = D0_NOTE
        learners.append([l if l != "all" else "all five (dataset and learner pairs)", int(r["miss_datasets"])] +
                        [int(r[f"datasets_{c}"]) for c in d0.CAUSES] + shares +
                        [int(r["ties"]), int(r["no_missed_fold"]), note])
    pairs = int(every["all"]["miss_datasets"])
    best = setter(line)
    most = max(int(every["all"][f"datasets_{c}"]) for c in d0.FOUR)
    order_rows = []
    for c in d0.CAUSES:
        n = int(every["all"][f"datasets_{c}"])
        order_rows.append([c, n, pairs, n / pairs if pairs else None,
                           "outside the four" if c == "chosen" else ("yes" if n == most and n else "no"),
                           "never (entry D0 causes)" if c == "chosen" else ("yes" if c == best else "no")])
    if best is None:
        meaning = DEFAULT_ORDER
    else:
        meaning = FIRST[best] if a.order_rest == "first-only" else AFTER[best]
    order_header = ["cause", "miss_pairs_in_cause", "miss_pairs", "share_of_miss_pairs", "largest_of_the_four",
                    "sets_the_order"]
    d0_rows = []
    for (d, l), r in sorted(by_d0.items()):
        if r["in_miss"] != "yes":
            continue
        sp, ab = names[d]
        d0_rows.append([sp, ab, d, l, float(r["program_ba"]), float(r["learner_lineage_ba"]), float(r["gap"]),
                        int(r["missed_folds"]), r["fallback"], r["cause"] if r["cause"] != "NA" else "no missed fold"])
    d0_header = ["species", "antibiotic", "dataset", "learner", "amrfinderplus_program_lineage_ba",
                 "learner_lineage_ba", "gap_amrfinderplus_program_minus_learner", "missed_folds", "fallback", "cause"]
    return {
        "page-species.tsv": (species_header, species_rows),
        "page-d0-learners.tsv": (learners_header, learners),
        "page-d0-order.tsv": (order_header, order_rows),
        "page-d0-datasets.tsv": (d0_header, d0_rows),
        "page-c8-summary.tsv": (c8.COMBINED_SUMMARY_HEADER, c8_summary),
        "page-c8-datasets.tsv": (c8.COMBINED_HEADER, c8_long),
    }, meaning


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--code-sha256", required=True, help="the build's code sha256 every DONE must hold")
    ap.add_argument("--species", type=int, required=True)
    ap.add_argument("--datasets", type=int, required=True)
    ap.add_argument("--d0-datasets", type=int, required=True)
    ap.add_argument("--d0-species", type=int, required=True)
    ap.add_argument("--d0-pairs", type=int, required=True)
    ap.add_argument("--d0-learner-pairs", help="learner=count,...: the miss pairs of each learner")
    ap.add_argument("--marks", help="the label count's table (datasets.tsv): the datasets must be those not M")
    ap.add_argument("--record", default=miss_sets.REF, help="the comparison record giving the MISS set")
    ap.add_argument("--names", default=c8.NAMES, help="dataset_id, species_slug and antibiotic per dataset")
    ap.add_argument("--equal-tables", default=",".join(NEED_EQUAL),
                    help=f"check tables every row of which must read equal, of {', '.join(TABLES)}")
    ap.add_argument("--c8-one-fold", required=True, choices=("zero", "undefined"))
    ap.add_argument("--c8-over-five", required=True, choices=("pairs", "dataset-median"))
    ap.add_argument("--c8-bands", required=True, help="none, or LOW,HIGH")
    ap.add_argument("--order-rest", required=True, choices=("first-only", "default-after"),
                    help="the order line names the first block only, or the first block then the default order")
    ap.add_argument("species_dirs", nargs="+")
    a = ap.parse_args(argv)
    a.equal_tables = set(t for t in a.equal_tables.split(",") if t)
    try:
        if not SHA_RE.match(a.code_sha256):
            raise Stop("--code-sha256 is not a sha256")
        if not a.equal_tables <= set(TABLES) or not set(NEED_EQUAL) <= a.equal_tables:
            raise Stop(f"--equal-tables must hold {' and '.join(NEED_EQUAL)}, and only tables of {', '.join(TABLES)}")
        for name in OUTPUTS:
            if os.path.exists(os.path.join(a.out, name)):
                raise Stop(f"{os.path.join(a.out, name)} exists")
        names = c8.read_names(a.names)
        spp = [read_species(os.path.normpath(p), a.code_sha256, a.equal_tables) for p in a.species_dirs]
        lines = gates(spp, a, names)
        with tempfile.TemporaryDirectory() as tmp:
            d0_out, every, line = run_d0(spp, a, tmp)
            lines.append(f"D0 pairs: {every['all']['miss_datasets']} of {a.d0_pairs}, by d0.py combine")
            c8_long, c8_summary = c8.combine_tables([(s["path"], s["c8"]) for s in spp], a.c8_one_fold,
                                                    a.c8_over_five, a.c8_bands, names)
            tables, meaning = page_tables(spp, a, names, every, line, c8_long, c8_summary)
            text = summary_text(spp, a, lines, every, line, meaning)
            os.makedirs(a.out, exist_ok=True)
            shutil.copytree(d0_out, os.path.join(a.out, "d0"))
    except (Stop, c8.CombineError, fr.ReportError) as e:
        print(f"stop: {e}", file=sys.stderr)
        return 2
    os.makedirs(os.path.join(a.out, "c8"), exist_ok=True)
    fr.write_tsv(os.path.join(a.out, "c8", "c8-combined-datasets.tsv"), c8.COMBINED_HEADER, c8_long)
    fr.write_tsv(os.path.join(a.out, "c8", "c8-combined-summary.tsv"), c8.COMBINED_SUMMARY_HEADER, c8_summary)
    for name, (header, rows) in tables.items():
        fr.write_tsv(os.path.join(a.out, name), header, rows)
    with open(os.path.join(a.out, "gates.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    with open(os.path.join(a.out, "summary.txt"), "w", encoding="utf-8") as f:
        f.write(text)
    print(text, end="")
    return 0


def summary_text(spp, a, lines, every, line, meaning):
    out = ["D0 and C8 of the refit over its species", ""] + ["gate passed: " + x for x in lines] + [
        "",
        f"rules: C8 one fold that chose: {a.c8_one_fold}; C8 over the five: {a.c8_over_five}; "
        f"C8 bands: {a.c8_bands}; order line: {a.order_rest}",
        "",
        f"D0 ordering ({D0_NOTE}): {line}",
        f"order of the input group: {meaning}",
        f"coverage: D0 reads {a.d0_datasets} of the plan's {a.plan_datasets} MISS datasets and "
        f"{every['all']['miss_datasets']} of its {a.plan_pairs} miss pairs; the MISS datasets that are M wait for "
        "the maintainer's label rule, and the line is read again then",
        "ruling 6: " + stability({c: int(every["all"][f"datasets_{c}"]) for c in d0.FOUR}, setter(line),
                                 a.plan_pairs, int(every["all"]["miss_datasets"]))[1],
        "",
        "C8 on the release decides nothing alone (the plan's row C8, Release).",
    ]
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    sys.exit(main())
