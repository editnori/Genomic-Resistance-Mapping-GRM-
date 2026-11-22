"""Reads the small-cohort run of plan block C1 part B into the tables of a public page (docs/cohort-draw.md).

  python3 -I read.py --out DIR --code-sha256 SHA --scorer SCORE_PY --scorer-sha256 SHA
      --species N --datasets N --results N --per-dataset N --names DATASETS_TSV
      --measures-over available|shared --line-50 beside|release-row --twin-rule beside|void
      --dataset-differences no|yes
      --pair FIRST SECOND [--pair FIRST SECOND ...]

Each --pair is one species: FIRST is the copy of run_species.sh's OUT (DONE,
binary.txt, expect.tsv, lists.txt, SHA256SUMS, results/), whose score/ this
script never opens; SECOND is the copy of the second scoring of the same
fitted results (score/ with checks.tsv, scores.tsv, lines.tsv and count.tsv,
logs/score.log, SHA256SUMS). Only SECOND's scores are read.

The reading stops, prints `stop: <reason>`, exits 2 and writes nothing unless:
  1. the scorer file has the sha256 given;
  2. every FIRST holds DONE in the run script's form and no FAILED, DONE's
     code sha256 is the one given, binary.txt names DONE's binary and code,
     SHA256SUMS has the sha256 DONE names, and the files read here
     (binary.txt, expect.tsv, lists.txt) and every result file expect.tsv
     names are listed in it, the read ones with their sha256;
  3. lists.txt, expect.tsv and DONE name one species; expect.tsv holds, for
     each dataset of lists.txt and no other, exactly the results of its lines,
     draws and twin, which are --per-dataset results;
  4. every SECOND lists exactly score/checks.tsv, score/count.tsv,
     score/lines.tsv, score/scores.tsv and logs/score.log in SHA256SUMS, each
     with its sha256, and score/ holds no other file;
  5. logs/score.log's checks: line says every result of expect.tsv passes
     every check, and checks.tsv holds that many results, each with the same
     checks, every one ok;
  6. scores.tsv and lines.tsv hold exactly the rows expect.tsv's results
     give, by key, for the species, each dataset with one antibiotic equal
     to the names file's;
  7. --species species, none twice; --datasets datasets, none twice, each
     named in the names file under its species; --results results;
  8. each species' count.tsv is byte for byte what the scorer's --count
     writes from that species' lines.tsv.
Then the scorer's own second form (`score.py --count`) counts over every
species' lines.tsv. This script writes no second count: its tables copy the
scorer's cells, and the entry of each dataset it lists is checked against the
scorer's entered, not-scored and no-line-at-25 counts (and, with
--dataset-differences yes, its median difference) before anything is written.

Arguments with no default are the lead's answers to READING-draft.md's
questions: --measures-over (Question 1), --line-50 (2), --twin-rule (3),
--dataset-differences (7).

Writes into --out (refused if any output exists): gates.txt, summary.txt,
page-deciding.tsv, page-twin.tsv, page-other.tsv, page-species.tsv,
page-datasets.tsv, page-entry.tsv, page-measures.tsv and count/count.tsv
(the scorer's). Outputs are aggregates: a row per species, dataset, line and
learner, never per genome; no genome id, k-mer, sequence, position or
feature. Standard library only; the scorer it calls needs numpy.
"""

import argparse
import csv
import hashlib
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
from fractions import Fraction

LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
TRUTHS = ["all", "new_lineage"]
BASE_LINE = "25"
TWIN_SEED = "45"
SEED = "42"
MIN_DRAWS = 3
DECIDING = ("lineage", "all", "levelled")
EXPECT_COLUMNS = ["dataset", "species_slug", "r", "draw", "seed", "study", "cohort", "cohort_sha256", "file"]
SECOND_FILES = ["logs/score.log", "score/checks.tsv", "score/count.tsv", "score/lines.tsv", "score/scores.tsv"]
DONE_FORM = re.compile(r"DONE (\S+) (\S+) binary ([0-9a-f]{64}) code ([0-9a-f]{64}) SHA256SUMS ([0-9a-f]{64})")
CHECKS_LINE = re.compile(r"checks: (\d+) of (\d+) results pass every check; (\d+) failed checks of (\d+)")
RATES = ["recall", "specificity", "very_major", "major"]
LABELS = ("not M: the same labels under label rules (a), (b) and (d); under rule (c), the most resistant call, "
          "a dataset holding a pair with an R and an S record is checked again")
CANDIDATE = "candidate under the plan's release row; the maintainer then decides whether to fetch and fit"
NOT_CANDIDATE = "not a candidate"
# Ruling 3 of 20:03Z, the twin beside: its reading, and the words a line that
# is not a candidate takes when a change of fold seed alone fails the row.
TWIN_HOLDS = "the row holds against a change of fold seed alone"
TWIN_FAILS = "a change of fold seed alone fails the row"
NOT_SHOWN = "not shown to be within the margin"
OUTPUTS = ["gates.txt", "summary.txt", "page-deciding.tsv", "page-twin.tsv", "page-other.tsv", "page-species.tsv",
           "page-datasets.tsv", "page-entry.tsv", "page-measures.tsv", "count"]


class Stop(Exception):
    """An input the reading does not read on from."""


# ---------------------------------------------------------------- files

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path, need=None):
    if not os.path.isfile(path):
        raise Stop(f"{path}: missing")
    with open(path, newline="") as f:
        rows = list(csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE))
    if not rows:
        raise Stop(f"{path}: empty")
    head = rows[0]
    if need is not None and [c for c in need if c not in head]:
        raise Stop(f"{path}: lacks the columns {[c for c in need if c not in head]}")
    out = []
    for i, cells in enumerate(rows[1:], start=2):
        if len(cells) != len(head):
            raise Stop(f"{path}: line {i} has {len(cells)} cells, want {len(head)}")
        out.append(dict(zip(head, cells)))
    return head, out


def sums(path):
    """{path: sha256} of a sha256sum list; a path twice stops."""
    if not os.path.isfile(path):
        raise Stop(f"{path}: missing")
    out = {}
    with open(path) as f:
        for i, line in enumerate(f, start=1):
            m = re.fullmatch(r"([0-9a-f]{64}) [ *](\S+)\n?", line)
            if not m:
                raise Stop(f"{path}: line {i} is not a sha256sum line")
            if m.group(2) in out:
                raise Stop(f"{path}: {m.group(2)} twice")
            out[m.group(2)] = m.group(1)
    return out


def species_name(slug):
    text = slug.replace("_", " ")
    return text[:1].upper() + text[1:]


def line_name(r, seed):
    return f"{r}-s{TWIN_SEED}" if seed == TWIN_SEED else str(r)


def line_order(name):
    return (int(name.split("-")[0]), name)


# ---------------------------------------------------------------- gates

def read_first(first, code_sha, per_dataset):
    """Gates 2 and 3 for one FIRST: (species, lists, expect rows, study datasets, DONE time)."""
    if os.path.exists(os.path.join(first, "FAILED")):
        raise Stop(f"{first}: FAILED is present")
    done_path = os.path.join(first, "DONE")
    if not os.path.isfile(done_path):
        raise Stop(f"{first}: no DONE")
    with open(done_path) as f:
        text = f.read().strip()
    m = DONE_FORM.fullmatch(text)
    if not m:
        raise Stop(f"{first}: DONE is not in the run script's form")
    sp, when, binary, code, sums_sha = m.groups()
    if code != code_sha:
        raise Stop(f"{sp}: DONE has code sha256 {code}, the given one is {code_sha}")
    listed = sums(os.path.join(first, "SHA256SUMS"))
    if sha256(os.path.join(first, "SHA256SUMS")) != sums_sha:
        raise Stop(f"{first}: SHA256SUMS is not the one DONE names")
    for name in ("binary.txt", "expect.tsv", "lists.txt"):
        p = os.path.join(first, name)
        if name not in listed or not os.path.isfile(p) or sha256(p) != listed[name]:
            raise Stop(f"{first}: {name} is missing, not listed or differs from its sha256 in SHA256SUMS")
    with open(os.path.join(first, "binary.txt")) as f:
        info = dict(line.rstrip("\n").split(" ", 1) for line in f if " " in line)
    if info.get("grm") != binary or info.get("code") != code:
        raise Stop(f"{first}: binary.txt does not name DONE's binary and code")
    lists = {}
    with open(os.path.join(first, "lists.txt")) as f:
        for line in f:
            k, _, v = line.rstrip("\n").partition("\t")
            lists[k] = v
    if lists.get("species") != sp:
        raise Stop(f"{first}: lists.txt names species {lists.get('species')}, DONE {sp}")
    datasets = lists.get("datasets", "").split(",")
    if not datasets or datasets == [""] or len(set(datasets)) != len(datasets):
        raise Stop(f"{first}: lists.txt's datasets are empty or hold one twice")
    for d in datasets:
        if not d.startswith(sp + "__"):
            raise Stop(f"{first}: {d} is not a dataset of {sp}")
    try:
        lines = [int(x) for x in lists["lines"].split(",")]
        draws = int(lists["draws"])
        twin = lists["twin"]
    except (KeyError, ValueError):
        raise Stop(f"{first}: lists.txt lacks lines, draws or twin") from None
    if int(BASE_LINE) not in lines or twin not in ("0", "1") or len(set(lines)) != len(lines):
        raise Stop(f"{first}: lists.txt's lines lack {BASE_LINE}, hold one twice, or its twin is not 0 or 1")
    grid = {(str(r), str(j), SEED) for r in lines for j in range(draws)}
    if twin == "1":
        grid |= {(BASE_LINE, str(j), TWIN_SEED) for j in range(draws)}
    if len(grid) != per_dataset:
        raise Stop(f"{sp}: lists.txt gives {len(grid)} results a dataset, want {per_dataset}")
    head, expect = read_tsv(os.path.join(first, "expect.tsv"))
    if head != EXPECT_COLUMNS:
        raise Stop(f"{first}/expect.tsv: columns {head}, want {EXPECT_COLUMNS}")
    seen, by_ds, study = set(), {}, {}
    for e in expect:
        key = (e["dataset"], e["r"], e["draw"], e["seed"])
        if key in seen:
            raise Stop(f"{first}/expect.tsv: {key} twice")
        seen.add(key)
        if e["species_slug"] != sp or not e["dataset"].startswith(sp + "__"):
            raise Stop(f"{first}/expect.tsv: {e['dataset']} is not of species {sp}")
        if e["file"] != f"{e['dataset']}.r{e['r']}.d{e['draw']}.s{e['seed']}.cohort.json":
            raise Stop(f"{first}/expect.tsv: {e['dataset']} r {e['r']} draw {e['draw']}: the result's name differs")
        if f"results/cohorts/{e['file']}" not in listed:
            raise Stop(f"{first}: SHA256SUMS does not list the result {e['file']}")
        by_ds.setdefault(e["dataset"], set()).add(key[1:])
        study.setdefault(e["dataset"], set()).add(e["study"])
    if sorted(by_ds) != sorted(datasets):
        raise Stop(f"{first}: expect.tsv's datasets are not lists.txt's: "
                   f"{len(set(datasets) - set(by_ds))} missing, {len(set(by_ds) - set(datasets))} not listed")
    for d, keys in sorted(by_ds.items()):
        if keys != grid:
            raise Stop(f"{sp}: {d} holds {len(keys)} results of the {len(grid)} its lines, draws and twin give")
        if len(study[d]) != 1 or study[d] - {"0", "1"}:
            raise Stop(f"{first}/expect.tsv: {d} has study marks {sorted(study[d])}")
    lines_all = sorted({line_name(r, s) for r, _, s in grid}, key=line_order)
    return {"species": sp, "done": when, "datasets": datasets, "expect": expect, "lines": lines_all,
            "draws": draws, "study": sorted(d for d in datasets if study[d] == {"1"})}


def read_second(second, sp, expect, study, lines, draws, names):
    """Gates 4 to 6 for one SECOND: its checks line, checks, scores and lines rows."""
    listed = sums(os.path.join(second, "SHA256SUMS"))
    if sorted(listed) != SECOND_FILES:
        raise Stop(f"{second}: SHA256SUMS lists {sorted(listed)}, want {SECOND_FILES}")
    for name in SECOND_FILES:
        p = os.path.join(second, name)
        if not os.path.isfile(p) or sha256(p) != listed[name]:
            raise Stop(f"{second}: {name} is missing or differs from its sha256 in SHA256SUMS")
    extra = sorted(set(os.listdir(os.path.join(second, "score"))) - {n[6:] for n in SECOND_FILES if n.startswith("score/")})
    if extra:
        raise Stop(f"{second}: score/ holds {extra}, which SHA256SUMS does not list")
    with open(os.path.join(second, "logs/score.log"), errors="replace") as f:
        found = [line.strip() for line in f if line.startswith("checks:")]
    if len(found) != 1:
        raise Stop(f"{second}/logs/score.log: {len(found)} checks: lines, want 1")
    m = CHECKS_LINE.fullmatch(found[0])
    if not m:
        raise Stop(f"{second}/logs/score.log: the checks: line is not in the scorer's form")
    passed, total, failed, checks_n = (int(x) for x in m.groups())
    if not (passed == total == len(expect)) or failed != 0:
        raise Stop(f"{second}/logs/score.log: {passed} of {total} results pass every check with {failed} failed "
                   f"checks; want {len(expect)} of {len(expect)} and 0")
    _, checks = read_tsv(os.path.join(second, "score/checks.tsv"), ["file", "check", "ok", "detail"])
    if len(checks) != checks_n:
        raise Stop(f"{second}/score/checks.tsv: {len(checks)} rows, the checks: line says {checks_n}")
    by_file = {}
    for c in checks:
        names_ = by_file.setdefault(c["file"], [])
        if c["check"] in names_:
            raise Stop(f"{second}/score/checks.tsv: {c['file']} holds the check {c['check']!r} twice")
        names_.append(c["check"])
        if c["ok"] != "1":
            raise Stop(f"{second}/score/checks.tsv: {c['file']} fails the check {c['check']!r}")
    files = {e["file"] for e in expect}
    if set(by_file) != files:
        raise Stop(f"{second}/score/checks.tsv: {len(files - set(by_file))} results of expect.tsv missing, "
                   f"{len(set(by_file) - files)} not in it")
    if len({tuple(sorted(v)) for v in by_file.values()}) != 1:
        raise Stop(f"{second}/score/checks.tsv: the results do not all hold the same checks")

    key_cols = ["dataset", "species_slug", "antibiotic", "r", "draw", "seed", "split", "truth_set", "learner", "available",
                "error"]
    need = key_cols + [f"{w}_{n}_{s}" for w in ("estimate", "truth") for n in RATES for s in ("k", "n")] + \
        [f"{n}_covers" for n in RATES] + [f"{w}_log_loss{s}" for w in ("estimate", "truth") for s in ("", "_n")]
    _, scores = read_tsv(os.path.join(second, "score/scores.tsv"), need)
    want = set()
    for e in expect:
        splits = ["lineage", "random"] + (["study"] if e["dataset"] in study else [])
        for s in splits:
            for t in TRUTHS:
                for l in LEARNERS:
                    want.add((e["dataset"], e["r"], e["draw"], e["seed"], s, t, l))
    got, antibiotic = set(), {}
    for r in scores:
        k = (r["dataset"], r["r"], r["draw"], r["seed"], r["split"], r["truth_set"], r["learner"])
        if k in got:
            raise Stop(f"{second}/score/scores.tsv: {k} twice")
        got.add(k)
        if r["species_slug"] != sp:
            raise Stop(f"{second}/score/scores.tsv: {r['dataset']} is marked species {r['species_slug']}")
        antibiotic.setdefault(r["dataset"], set()).add(r["antibiotic"])
    if got != want:
        raise Stop(f"{second}/score/scores.tsv: {len(want - got)} rows of expect.tsv's results missing, "
                   f"{len(got - want)} rows not in them")
    for d, a in sorted(antibiotic.items()):
        if d not in names:
            raise Stop(f"the names file does not name {d}")
        if a != {names[d]["antibiotic"]}:
            raise Stop(f"{second}/score/scores.tsv: {d} has antibiotic {sorted(a)}, the names file "
                       f"{names[d]['antibiotic']!r}")

    need = ["dataset", "species_slug", "split", "truth_set", "learner", "line", "draws", "available",
            "available_draws", "errors"]
    _, lrows = read_tsv(os.path.join(second, "score/lines.tsv"), need)
    datasets = sorted({e["dataset"] for e in expect})
    want = {(d, s, t, l, ln) for d in datasets for s in ["lineage", "random"] + (["study"] if d in study else [])
            for t in TRUTHS for l in LEARNERS for ln in lines}
    got = set()
    for r in lrows:
        k = (r["dataset"], r["split"], r["truth_set"], r["learner"], r["line"])
        if k in got:
            raise Stop(f"{second}/score/lines.tsv: {k} twice")
        got.add(k)
        if r["species_slug"] != sp or r["draws"] != str(draws):
            raise Stop(f"{second}/score/lines.tsv: {k} is marked species {r['species_slug']} with {r['draws']} draws")
    if got != want:
        raise Stop(f"{second}/score/lines.tsv: {len(want - got)} rows missing, {len(got - want)} rows not expected")
    # The entry and the count read lines.tsv, Rules' measures scores.tsv: each lines row's available draws and
    # errors must be, as text, what the scorer writes from scores.tsv, so both speak of the same results.
    avail = {}
    for r in scores:
        if r["available"] == "1":
            avail.setdefault((r["dataset"], r["split"], r["truth_set"], r["learner"], line_name(r["r"], r["seed"])),
                             []).append((int(r["draw"]), r["error"]))
    for r in lrows:
        k = (r["dataset"], r["split"], r["truth_set"], r["learner"], r["line"])
        av = sorted(avail.get(k, []))
        if r["available_draws"] != ",".join(str(d) for d, _ in av) or r["errors"] != ",".join(e for _, e in av):
            raise Stop(f"{second}/score/lines.tsv: {k}: its available draws or errors are not scores.tsv's")
    return scores, lrows


def run_count(python, scorer, out, line_files):
    """The scorer's own --count over line_files into out; its count.tsv's bytes."""
    r = subprocess.run([python, "-I", "-B", scorer, "--count", out] + list(line_files),
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise Stop(f"the scorer's --count ended with {r.returncode}: {(r.stderr or r.stdout).strip()[-300:]}")
    with open(os.path.join(out, "count.tsv"), "rb") as f:
        return f.read()


# ---------------------------------------------------------------- entry, measures

def draws_of(row):
    """{draw: error text} of a lines.tsv row, read as score.py's per_draw reads it."""
    ds = row["available_draws"].split(",") if row["available_draws"] else []
    es = row["errors"].split(",") if row["errors"] else []
    if len(ds) != len(es):
        raise Stop(f"lines.tsv: {row['dataset']} {row['line']}: {len(es)} errors for {len(ds)} available draws")
    return {int(d): float(e) for d, e in zip(ds, es)}


def entry(lrows, lines):
    """Per (split, truth, learner, line) not 25 and dataset: the line's and 25's draws, the shared, the reason."""
    cell = {(r["split"], r["truth_set"], r["learner"], r["line"], r["dataset"]): r for r in lrows}
    out = {}
    for (s, t, l, ln, d), row in sorted(cell.items()):
        if ln == BASE_LINE:
            continue
        own = draws_of(row)
        base_row = cell.get((s, t, l, BASE_LINE, d))
        base = draws_of(base_row) if base_row is not None else {}
        shared = sorted(set(own) & set(base))
        if base_row is None:
            why = "no line at 25"
        elif len(shared) < MIN_DRAWS:
            why = f"{len(shared)} shared draws, under {MIN_DRAWS}"
        else:
            why = ""
        diff = round(sum(own[j] - base[j] for j in shared) / len(shared), 12) if not why else None
        out[(s, t, l, ln, d)] = {"line_draws": sorted(own), "base_draws": sorted(base), "shared": shared,
                                 "why": why, "difference": diff, "species": row["species_slug"]}
    return out


def check_entry(ent, count_rows, differences):
    """The entry listed here against the scorer's counts, for every learner row of count.tsv."""
    for c in count_rows:
        if c["row"] != "learner":
            continue
        key = (c["split"], c["truth_set"], c["learner"], c["line"])
        mine = [v for k, v in ent.items() if k[:4] == key]
        entered = [v for v in mine if not v["why"]]
        got = (len(mine), len(entered), len(mine) - len(entered), sum(1 for v in mine if v["why"] == "no line at 25"))
        want = tuple(int(c[k]) for k in ("datasets", "entered", "not_scored", "no_base_line"))
        if got != want:
            raise Stop(f"the entry listed here, {got}, differs from the scorer's {want} at {key}")
        if differences:
            med = repr(round(float(statistics.median([v["difference"] for v in entered])), 12)) if entered else ""
            if med != c["median_difference"]:
                raise Stop(f"the median of the differences listed here differs from the scorer's at {key}")


def measures(scores, ent, over, lines):
    """Rules' measures beside the error, per split, truth set, learner, line and side, with denominators."""
    by = {}
    for r in scores:
        # Every split, truth set, learner and line gets its rows, with 0 results
        # where none is available (ruling 6: the file holds every row).
        rows = by.setdefault((r["split"], r["truth_set"], r["learner"], line_name(r["r"], r["seed"])), [])
        if r["available"] == "1":
            rows.append(r)
    groups = []
    for (s, t, l, ln), rows in sorted(by.items(), key=lambda kv: (kv[0][:3], line_order(kv[0][3]))):
        if over == "available":
            groups.append(((s, t, l, ln, "line"), rows))
        elif ln != BASE_LINE:
            keep = {(d, j) for (s2, t2, l2, ln2, d), v in ent.items() if (s2, t2, l2, ln2) == (s, t, l, ln)
                    and not v["why"] for j in v["shared"]}
            groups.append(((s, t, l, ln, "line"), [r for r in rows if (r["dataset"], int(r["draw"])) in keep]))
            base = by.get((s, t, l, BASE_LINE), [])
            groups.append(((s, t, l, ln, BASE_LINE), [r for r in base if (r["dataset"], int(r["draw"])) in keep]))
    out = []
    for (s, t, l, ln, side), rows in groups:
        datasets = sorted({r["dataset"] for r in rows})
        for who in ("estimate", "truth"):
            for name in RATES + ["log_loss"]:
                per_ds, k_sum, n_sum, loss_sum = {}, 0, 0, 0.0
                for r in rows:
                    if name == "log_loss":
                        v, n = r[f"{who}_log_loss"], int(r[f"{who}_log_loss_n"] or 0)
                        if v == "" or n == 0:
                            continue
                        loss_sum += float(v) * n
                        n_sum += n
                        per_ds.setdefault(r["dataset"], []).append(float(v))
                    else:
                        k, n = int(r[f"{who}_{name}_k"]), int(r[f"{who}_{name}_n"])
                        if n == 0:
                            continue
                        k_sum, n_sum = k_sum + k, n_sum + n
                        per_ds.setdefault(r["dataset"], []).append(k / n)
                pooled = (loss_sum / n_sum if name == "log_loss" else k_sum / n_sum) if n_sum else None
                means = [sum(v) / len(v) for v in per_ds.values()]
                cov = [r[f"{name}_covers"] for r in rows if name != "log_loss" and r[f"{name}_covers"] != ""]
                out.append({
                    "split": s, "truth_set": t, "learner": l, "line": ln, "side": side,
                    "over": "own available results" if over == "available" else "shared draws of the datasets that enter",
                    "who": who, "measure": name, "results": len(rows),
                    "pooled_count": "" if name == "log_loss" else k_sum, "pooled_denominator": n_sum,
                    "pooled_value": f6(pooled), "datasets": len(datasets), "datasets_defined": len(per_ds),
                    "median_over_datasets": f6(statistics.median(means)) if means else "",
                    "truth_inside_estimate_wilson": "" if who == "estimate" or name == "log_loss" else sum(1 for c in cov if c == "1"),
                    "of_results": "" if who == "estimate" or name == "log_loss" else len(cov),
                    "decides": "no: beside the error", "labels": LABELS})
    return out


# ---------------------------------------------------------------- tables

def f6(x):
    if x is None:
        return ""
    s = "%.6f" % float(x)
    return "0.000000" if s == "-0.000000" else s


def frac(k, n):
    return f6(Fraction(int(k), int(n))) if n not in ("", "0", 0) else ""


def why_not(c, line_50):
    why = []
    if c["split"] == "random":
        why.append("the random split")
    if c["split"] == "study":
        why.append("the study split")
    if c["truth_set"] == "new_lineage":
        why.append("truth on new lineages")
    if c["coverage"] == "fold_range":
        why.append("the fold range's coverage")
    if c["line"].endswith(f"-s{TWIN_SEED}"):
        why.append("the twin")
    if c["line"] == "50" and line_50 == "beside":
        why.append("line 50 is not a lower line")
    return why


DECIDING_COLUMNS = ["line", "learner", "row", "datasets", "entered", "not_scored", "no_base_line", "better", "worse",
                    "same", "median_difference", "interval_lo", "interval_hi", "species", "shared_draws", "covered",
                    "coverage_n", "coverage_at_line", "base_covered", "base_coverage_n", "coverage_at_25",
                    "coverage_minus_25", "line_covered", "line_coverage_n", "base_covered_own",
                    "base_coverage_n_own", "median_width", "base_median_width", "meets", "learners_meeting",
                    "meeting", "candidate", "reading", "labels"]


# The scorer's candidate cell is only "3 or more learners meet". Off the
# deciding rows that is no candidate (the twin never is, ruling of 19:01Z), so
# the twin's and the other rows' tables name the cell for what it counts.
MEET3 = "three_or_more_meet"
TWIN_COLUMNS = [MEET3 if k == "candidate" else k for k in DECIDING_COLUMNS]


def learner_reading(c):
    """Ruling 10: "meets", or which part of the release row a learner row fails, from the scorer's own cells.

    The parts are the interval's upper end (not wholly below 0.02, read at
    the 12 decimals the file writes), the coverage (not within 0.05 of 25's,
    in exact fractions) and no dataset entered. The scorer's meets cell must
    agree with them; nothing is counted again.
    """
    if c["entered"] in ("", "0"):
        parts = ["no dataset entered"]
    else:
        parts = []
        if c["interval_hi"] == "" or Fraction(c["interval_hi"]) >= Fraction(1, 50):
            parts.append("the interval's upper end")
        n, bn = c["coverage_n"], c["base_coverage_n"]
        if n in ("", "0") or bn in ("", "0") or abs(Fraction(int(c["covered"]), int(n))
                                                   - Fraction(int(c["base_covered"]), int(bn))) > Fraction(1, 20):
            parts.append("the coverage")
    if (c["meets"] == "1") != (not parts):
        raise Stop(f"count.tsv: {c['line']} {c['learner']}: the meets cell {c['meets']!r} differs from its cells")
    return "meets" if not parts else "does not meet: " + " and ".join(parts)


def deciding_row(c, reading):
    out = {k: c.get(k, "") for k in DECIDING_COLUMNS}
    out["learner"] = c["learner"]
    if c["row"] == "learner":
        out["coverage_at_line"] = frac(c["covered"], c["coverage_n"])
        out["coverage_at_25"] = frac(c["base_covered"], c["base_coverage_n"])
        if c["coverage_n"] not in ("", "0") and c["base_coverage_n"] not in ("", "0"):
            out["coverage_minus_25"] = f6(Fraction(int(c["covered"]), int(c["coverage_n"]))
                                          - Fraction(int(c["base_covered"]), int(c["base_coverage_n"])))
    out["reading"] = reading
    out["labels"] = LABELS
    return out


def write_tsv(path, columns, rows):
    with open(path, "w") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in columns) + "\n")


# ---------------------------------------------------------------- main

def parse(argv):
    ap = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--code-sha256", required=True)
    ap.add_argument("--scorer", required=True)
    ap.add_argument("--scorer-sha256", required=True)
    for k in ("species", "datasets", "results", "per-dataset"):
        ap.add_argument(f"--{k}", type=int, required=True)
    ap.add_argument("--names", required=True)
    ap.add_argument("--measures-over", choices=["available", "shared"], required=True)
    ap.add_argument("--line-50", choices=["beside", "release-row"], required=True)
    ap.add_argument("--twin-rule", choices=["beside", "void"], required=True)
    ap.add_argument("--dataset-differences", choices=["no", "yes"], required=True)
    ap.add_argument("--pair", nargs=2, action="append", required=True, metavar=("FIRST", "SECOND"))
    ap.add_argument("--python", default=sys.executable, help="the python that runs the scorer (needs numpy)")
    return ap.parse_args(argv)


def read(a, work):
    """Every gate, the scorer's count and the tables, in memory. Raises Stop."""
    gates = []
    if not os.path.isfile(a.scorer) or sha256(a.scorer) != a.scorer_sha256:
        raise Stop(f"{a.scorer}: missing or its sha256 is not the given {a.scorer_sha256}")
    gates.append(f"1 scorer: {os.path.basename(a.scorer)} has sha256 {a.scorer_sha256}")
    _, name_rows = read_tsv(a.names, ["dataset_id", "species_slug", "antibiotic"])
    names = {}
    for r in name_rows:
        if r["dataset_id"] in names:
            raise Stop(f"{a.names}: {r['dataset_id']} twice")
        names[r["dataset_id"]] = r
    species, seen = [], {}
    for first, second in a.pair:
        f = read_first(first, a.code_sha256, a.per_dataset)
        sp = f["species"]
        if sp in seen:
            raise Stop(f"species {sp} twice: {seen[sp]}, {first}")
        seen[sp] = first
        for d in f["datasets"]:
            if d not in names or names[d]["species_slug"] != sp:
                raise Stop(f"the names file does not name {d} under {sp}")
        f["scores"], f["line_rows"] = read_second(second, sp, f["expect"], set(f["study"]), f["lines"], f["draws"], names)
        f["second"] = second
        species.append(f)
    gates.append(f"2 first copies: DONE with code sha256 {a.code_sha256}, binary.txt and SHA256SUMS agree, "
                 f"in {len(species)} of {len(species)}")
    gates.append(f"3 expect.tsv: every dataset holds its {a.per_dataset} results of lists.txt's lines, draws and twin")
    gates.append(f"4 second scorings: the 5 files of SHA256SUMS hold their sha256 in {len(species)} of {len(species)}")
    gates.append("5 checks: every result passes every check in every species")
    gates.append("6 scores.tsv and lines.tsv hold exactly expect.tsv's rows; antibiotics equal the names file")
    datasets = [d for f in species for d in f["datasets"]]
    results = sum(len(f["expect"]) for f in species)
    if len(species) != a.species:
        raise Stop(f"{len(species)} species, want {a.species}")
    if len(set(datasets)) != len(datasets):
        raise Stop("a dataset is listed twice")
    if len(datasets) != a.datasets:
        raise Stop(f"{len(datasets)} datasets, want {a.datasets}")
    if results != a.results:
        raise Stop(f"{results} results, want {a.results}")
    lines = species[0]["lines"]
    if any(f["lines"] != lines for f in species):
        raise Stop("the species do not share one set of lines")
    gates.append(f"7 counts: {len(species)} species, {len(datasets)} datasets, {results} results, none twice")
    for f in species:
        got = run_count(a.python, a.scorer, os.path.join(work, "species", f["species"]),
                        [os.path.join(f["second"], "score/lines.tsv")])
        with open(os.path.join(f["second"], "score/count.tsv"), "rb") as fh:
            if fh.read() != got:
                raise Stop(f"{f['second']}/score/count.tsv is not what the scorer's --count writes from its lines.tsv")
    gates.append(f"8 count.tsv: the scorer's --count of each species' lines.tsv, in {len(species)} of {len(species)}")
    run_count(a.python, a.scorer, os.path.join(work, "count"), [os.path.join(f["second"], "score/lines.tsv") for f in species])
    _, count_rows = read_tsv(os.path.join(work, "count", "count.tsv"),
                             ["row", "split", "truth_set", "learner", "line", "coverage", "decides", "candidate",
                              "learners_meeting", "meeting", "meets", "datasets", "entered", "not_scored",
                              "no_base_line", "median_difference"])
    lrows = [r for f in species for r in f["line_rows"]]
    ent = entry(lrows, lines)
    check_entry(ent, count_rows, a.dataset_differences == "yes")
    gates.append("the entry listed here equals the scorer's datasets, entered, not scored and no line at 25 "
                 "on every learner row" + (", and its median difference" if a.dataset_differences == "yes" else ""))
    return species, gates, count_rows, ent, lines


def tables(a, species, count_rows, ent, lines):
    lower = [ln for ln in lines if ln != BASE_LINE and not ln.endswith(f"-s{TWIN_SEED}")
             and (ln != "50" or a.line_50 == "release-row")]
    twin_lines = [ln for ln in lines if ln.endswith(f"-s{TWIN_SEED}")]
    view = lambda c: (c["split"], c["truth_set"], c["coverage"]) == DECIDING
    twin_row = [c for c in count_rows if view(c) and c["row"] == "line" and c["line"] in twin_lines]
    twin_meets = any(int(c["learners_meeting"]) >= 3 for c in twin_row)
    void = a.twin_rule == "void" and twin_meets
    # Ruling 3 (beside): the twin's reading fixes the words a line is read in.
    twin_fails = bool(twin_row) and not twin_meets
    beside_twin = "".join(f"; the twin: {c['learners_meeting']} of {len(LEARNERS)} learners meet" for c in twin_row)
    deciding, twin, other = [], [], []
    for c in count_rows:
        if view(c) and c["line"] in lower:
            if c["row"] == "line":
                reading = ("decides nothing: the twin meets the release row on 3 or more learners" if void
                           else CANDIDATE + beside_twin if c["candidate"] == "1"
                           else NOT_CANDIDATE + (": " + NOT_SHOWN if twin_fails else ""))
            else:
                reading = learner_reading(c)
            deciding.append(deciding_row(c, reading))
        elif view(c) and c["line"] in twin_lines:
            reading = ("the count with no change; never a candidate"
                       + ("; it meets the release row on 3 or more learners, so the block decides nothing"
                          if void else "; " + (TWIN_HOLDS if int(c["learners_meeting"]) >= 3 else TWIN_FAILS)
                          if c["row"] == "line" else ""))
            t = deciding_row(c, reading if c["row"] == "line" else learner_reading(c))
            t[MEET3] = t.pop("candidate")
            twin.append(t)
        else:
            r = {(MEET3 if k == "candidate" else k): v for k, v in c.items()}
            r["decides"] = "no: " + ", ".join(why_not(c, a.line_50))
            r["labels"] = LABELS
            other.append(r)
    deciding.sort(key=lambda r: (line_order(r["line"]), r["row"] == "line",
                                 LEARNERS.index(r["learner"]) if r["learner"] in LEARNERS else 9))
    names = {}
    _, rows = read_tsv(a.names, ["dataset_id", "species_slug", "antibiotic"])
    for r in rows:
        names[r["dataset_id"]] = r
    sp_rows, ds_rows = [], []
    for f in sorted(species, key=lambda f: f["species"]):
        ds = sorted(f["datasets"], key=lambda d: names[d]["antibiotic"])
        sp_rows.append({"species": species_name(f["species"]), "species_slug": f["species"], "datasets": len(ds),
                        "results": len(f["expect"]), "study_datasets": len(f["study"]),
                        "antibiotics": ", ".join(names[d]["antibiotic"] for d in ds), "labels": LABELS})
        for d in ds:
            ds_rows.append({"species": species_name(f["species"]), "antibiotic": names[d]["antibiotic"], "dataset": d,
                            "results": sum(1 for e in f["expect"] if e["dataset"] == d),
                            "study_plan": "yes" if d in f["study"] else "no", "labels": LABELS})
    entry_rows = []
    for (s, t, l, ln, d), v in sorted(ent.items(), key=lambda kv: (kv[0][4], line_order(kv[0][3]), LEARNERS.index(kv[0][2]))):
        if (s, t) != DECIDING[:2]:
            continue
        entry_rows.append({"species": species_name(v["species"]), "antibiotic": names[d]["antibiotic"], "dataset": d,
                           "line": ln, "learner": l, "line_draws": ",".join(map(str, v["line_draws"])),
                           "draws_at_25": ",".join(map(str, v["base_draws"])), "shared_draws": len(v["shared"]),
                           "entered": "yes" if not v["why"] else "no: " + v["why"],
                           "difference": "" if a.dataset_differences == "no" or v["difference"] is None
                           else repr(v["difference"]), "labels": LABELS})
    return deciding, twin, other, sp_rows, ds_rows, entry_rows, void, twin_row


def main(argv=None):
    a = parse(argv)
    try:
        for name in OUTPUTS:
            if os.path.exists(os.path.join(a.out, name)):
                raise Stop(f"{os.path.join(a.out, name)} exists")
        with tempfile.TemporaryDirectory() as work:
            species, gates, count_rows, ent, lines = read(a, work)
            deciding, twin, other, sp_rows, ds_rows, entry_rows, void, twin_row = tables(a, species, count_rows, ent, lines)
            meas = measures([r for f in species for r in f["scores"]], ent, a.measures_over, lines)
            os.makedirs(os.path.join(a.out, "count"), exist_ok=True)
            shutil.copyfile(os.path.join(work, "count", "count.tsv"), os.path.join(a.out, "count", "count.tsv"))
    except Stop as e:
        print(f"stop: {e}", file=sys.stderr)
        return 2
    write_tsv(os.path.join(a.out, "page-deciding.tsv"), DECIDING_COLUMNS, deciding)
    write_tsv(os.path.join(a.out, "page-twin.tsv"), TWIN_COLUMNS, twin)
    other_cols = list(other[0].keys()) if other else ["decides", "labels"]
    write_tsv(os.path.join(a.out, "page-other.tsv"), other_cols, other)
    write_tsv(os.path.join(a.out, "page-species.tsv"),
              ["species", "species_slug", "datasets", "results", "study_datasets", "antibiotics", "labels"], sp_rows)
    write_tsv(os.path.join(a.out, "page-datasets.tsv"),
              ["species", "antibiotic", "dataset", "results", "study_plan", "labels"], ds_rows)
    write_tsv(os.path.join(a.out, "page-entry.tsv"),
              ["species", "antibiotic", "dataset", "line", "learner", "line_draws", "draws_at_25", "shared_draws",
               "entered", "difference", "labels"], entry_rows)
    meas_cols = list(meas[0].keys()) if meas else ["labels"]
    write_tsv(os.path.join(a.out, "page-measures.tsv"), meas_cols, meas)
    with open(os.path.join(a.out, "gates.txt"), "w") as f:
        f.write("\n".join(gates) + "\n")
    text = summary(a, species, gates, deciding, twin_row, void)
    with open(os.path.join(a.out, "summary.txt"), "w") as f:
        f.write(text)
    print(text, end="")
    return 0


def summary(a, species, gates, deciding, twin_row, void):
    lines = ["Small cohorts against the 25-per-class line, plan block C1 part B: the reading's tables",
             f"species: {len(species)}: " + ", ".join(species_name(f["species"]) for f in sorted(species, key=lambda f: f["species"])),
             "gates passed:"] + [f"  {g}" for g in gates]
    lines += ["rules chosen: " + ", ".join(f"--{k.replace('_', '-')} {getattr(a, k)}" for k in
                                           ("measures_over", "line_50", "twin_rule", "dataset_differences")),
              "the row that decides: the lineage split, the truth on all scored genomes, the levelled interval"]
    for r in deciding:
        if r["row"] == "line":
            lines.append(f"line {r['line']}: {r['learners_meeting']} of 5 learners meet ({r['meeting'] or 'none'}): "
                         f"{r['reading']}")
    for r in twin_row:
        lines.append(f"twin {r['line']}: {r['learners_meeting']} of 5 learners meet; the count with no change, "
                     f"never a candidate" + ("; under --twin-rule void the block decides nothing" if void else
                                             "; " + (TWIN_HOLDS if int(r["learners_meeting"]) >= 3 else TWIN_FAILS)))
    lines += [f"labels: {LABELS}",
              "A candidate is a line the maintainer may then decide to fetch and fit. Nothing changes in the release."]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    sys.exit(main())
