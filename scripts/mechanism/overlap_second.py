"""overlap_second.py: the learner overlap's second pass (docs/learner-overlap-result.md, Plan, "Second pass").

    python3 -I overlap_second.py plan
    python3 -I overlap_second.py read REPORTS_DIR --first FIRST_RECORD --out READING_DIR [--expect 274]
    python3 -I overlap_second.py write READING_DIR --first FIRST_RECORD --out RECORD_DIR
    python3 -I overlap_second.py print RECORD_DIR

The first pass (learner_overlap.py) counted the five learners one by one. The release ships their mean, so this
pass counts the resistant genomes the five-learner mean misses at the 0.5 cutoff and the susceptible genomes it
calls resistant, each by how many of the five miss or call it. It also works the first pass's yardstick and ratio
again with the genomes of every outer fold where one of the five had a null fit left out.

`read` runs every check of the first pass again through learner_overlap.read_reports, compares each report's
first-pass columns, counted again, with FIRST_RECORD/datasets.tsv, and writes reading.json: per report, counts
only. `write` writes vote.tsv, vote-pooled.tsv, gates.txt and SHA256SUMS after the shape scan finds no genome id.
`print` prints the page's tables from a record.

Nothing is fitted, trained or scored again: the mean of five stored scores is arithmetic. Every key named `id` is
dropped as a report is parsed (learner_overlap.load_report), so no genome id is held, printed or written. A failed
stop check writes nothing and exits 1.

Python 3.12, standard library only; run with -I.
"""

import argparse
import json
import math
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import learner_overlap as lo  # noqa: E402

FIVE = lo.FIVE
VOTE_SUM = lo.CUTOFF * len(FIVE)    # 2.5: the mean of the five at the release's 0.5 cutoff
NEAR = 1e-9
READING = "reading.json"
RECORD_FILES = ("vote.tsv", "vote-pooled.tsv", "gates.txt")
FIRST_COLUMNS = (["genomes", "resistant", "susceptible"] + [f"miss_{k}" for k in lo.KS]
                 + [f"wrong_{k}" for k in lo.KS] + [f"miss_{n}" for n in FIVE] + ["miss_all5", "miss_all5_in_null_fold"])
COUNT_COLUMNS = (["mean_miss", "mean_miss_by_all5", "mean_miss_by_1to4", "mean_miss_by_none",
                  "mean_wrong", "mean_wrong_by_all5", "mean_wrong_by_1to4", "mean_wrong_by_none",
                  "nonull_resistant"] + [f"nonull_miss_{n}" for n in FIVE] + ["nonull_miss_all5"])
VOTE_COLUMNS = ["species", "antibiotic", "genomes", "resistant", "susceptible"] + COUNT_COLUMNS + [
    "nonull_expected_all5", "nonull_ratio_all5"]
POOLED_COLUMNS = (["scope", "species", "datasets", "genomes", "resistant", "susceptible"] + COUNT_COLUMNS
                  + ["nonull_expected_all5", "mean_share_median", "nonull_datasets_with_ratio", "nonull_ratio_median",
                     "nonull_ratio_q1", "nonull_ratio_q3", "nonull_datasets_ratio_over_1", "paired_datasets",
                     "paired_ratio_median", "paired_nonull_ratio_median"])

Stop = lo.Stop


# ---- one report ---------------------------------------------------------------------------------------------------

def vote(scores):
    """The five-learner mean at the cutoff: the sum by math.fsum is 2.5 or more."""
    return math.fsum(scores) >= VOTE_SUM


def plain_vote(scores):
    total = 0.0
    for x in scores:
        total += x
    return total / len(scores) >= lo.CUTOFF


def null_union(learners):
    return set().union(*(lo.null_folds(lr) for lr in learners.values()))


def bucket(k):
    """'all5', '1to4' or 'none': how many of the five miss (or call) the genome."""
    return "all5" if k == len(FIVE) else ("none" if k == 0 else "1to4")


def second_counts(labels, folds, learners):
    """Counts of one report; `learners` maps each of the five to its learner entry."""
    nulls = null_union(learners)
    out = {c: 0 for c in COUNT_COLUMNS}
    edge = {"plain_differs": 0, "near": 0}
    for i, t in enumerate(labels):
        scores = [learners[n]["scores"][i] for n in FIVE]
        called = sum(1 for x in scores if lo.is_call(x))
        mean_call = vote(scores)
        edge["plain_differs"] += plain_vote(scores) != mean_call
        edge["near"] += abs(math.fsum(scores) - VOTE_SUM) < NEAR
        if t:
            if not mean_call:
                out["mean_miss"] += 1
                out[f"mean_miss_by_{bucket(len(FIVE) - called)}"] += 1
            if folds[i] not in nulls:
                out["nonull_resistant"] += 1
                for n, x in zip(FIVE, scores):
                    out[f"nonull_miss_{n}"] += 0 if lo.is_call(x) else 1
                out["nonull_miss_all5"] += 1 if called == 0 else 0
        elif mean_call:
            out["mean_wrong"] += 1
            out[f"mean_wrong_by_{bucket(called)}"] += 1
    return out, edge


def first_agrees(row, first):
    """The first-pass columns counted again equal the committed datasets.tsv row."""
    return all(str(row.get(c, "")) == first.get(c) for c in FIRST_COLUMNS)


def arithmetic_holds(counts, row):
    return (counts["mean_miss_by_none"] == 0 and counts["mean_wrong_by_none"] == 0
            and counts["mean_miss_by_all5"] == row["miss_5"] and counts["mean_wrong_by_all5"] == row["wrong_5"])


def null_cut_holds(counts, row):
    return counts["nonull_miss_all5"] + row["miss_all5_in_null_fold"] == row["miss_all5"]


def scan(paths):
    return lo.id_hits(paths)


# ---- the read -----------------------------------------------------------------------------------------------------

def first_rows(first_record):
    return {(r["species"], r["antibiotic"]): r for r in lo.read_tsv(os.path.join(first_record, "datasets.tsv"))}


def read_second(directory, first_record, expect=lo.EXPECT_FILES):
    """Run the first pass's checks, then the second pass's; return the reading or raise Stop."""
    reading = lo.read_reports(directory, expect)          # raises Stop on any first-pass check
    if reading["found"] != expect:
        raise Stop(f"check 1: {reading['found']} files found, {expect} expected: STOP")
    first = first_rows(first_record)
    used = {(d["species"], d["antibiotic"]): d for d in reading["datasets"]}
    fails = {2: [], 3: [], 4: []}
    if set(used) != set(first):
        fails[2].append(f"{len(set(used) ^ set(first))} reports stand in one of the read and the record only")
    edges = {"plain_differs": 0, "near": 0}
    rows = []
    for name in lo.report_files(directory):
        report = lo.load_report(os.path.join(directory, name))
        key = (report.get("species"), report.get("antibiotic"))
        if key not in used:
            continue
        d = used[key]
        row1 = lo.dataset_row(d)
        who = f"{key[0]} with {key[1]}"
        if key in first and not first_agrees(row1, first[key]):
            fails[2].append(who)
        row = {"species": key[0], "antibiotic": key[1], "genomes": d["genomes"], "resistant": d["resistant"],
               "susceptible": d["susceptible"], "five": lo.is_five(d)}
        if row["five"]:
            split = lo.lineage_splits(report)[0]
            learners = {lr["learner"]: lr for lr in lo.learner_entries(split)}
            labels = [s.get("resistant") for s in report.get("samples") or []]
            counts, edge = second_counts(labels, split["test_fold"], learners)
            for k in edges:
                edges[k] += edge[k]
            if not arithmetic_holds(counts, row1):
                fails[3].append(who)
            if not null_cut_holds(counts, row1):
                fails[4].append(who)
            row.update(counts)
        del report
        rows.append(row)
    five = sum(1 for r in rows if r["five"])
    lines = [
        f"second check 1 first pass: {reading['found']} files found, every first-pass check passed, "
        f"{reading['used']} reports used: pass",
        f"second check 2 first-pass columns: {len(rows)} reports read again, {len(first)} rows in datasets.tsv, "
        f"{len(fails[2])} that differ: " + ("pass" if not fails[2] else "STOP"),
        f"second check 3 arithmetic of the mean: {five} reports with the five, {len(fails[3])} where the mean misses "
        f"a genome none of the five miss, calls one none call, or parts from all five: "
        + ("pass" if not fails[3] else "STOP"),
        f"second check 4 null-fit cut: {len(fails[4])} reports where the all-five misses outside and inside null-fit "
        f"folds do not add up: " + ("pass" if not fails[4] else "STOP"),
        f"second check 5 the vote's edge: {edges['plain_differs']} genomes where a plain sum over 5 gives another "
        f"call than math.fsum; {edges['near']} genomes whose sum lies within 1e-9 of 2.5",
    ]
    for check in sorted(fails):
        lines += [f"  second check {check} failed: {d}" for d in fails[check]]
    if any("STOP" in line for line in lines):
        raise Stop("\n".join(reading["checks"] + lines))
    for r in rows:
        r.pop("five")
    return {"cutoff": lo.CUTOFF, "scheme": lo.SCHEME, "five": list(FIVE), "vote_sum": VOTE_SUM,
            "first_checks": reading["checks"], "checks": lines, "edges": edges,
            "datasets": sorted(rows, key=lambda r: (r["species"], r["antibiotic"]))}


# ---- the record ---------------------------------------------------------------------------------------------------

def vote_row(r):
    row = {c: r.get(c, "") for c in VOTE_COLUMNS[:-2]}
    if "mean_miss" in r:
        exp = lo.expected_all5(r["nonull_resistant"], [r[f"nonull_miss_{n}"] for n in FIVE])
        row["nonull_expected_all5"] = lo.d12(exp)
        row["nonull_ratio_all5"] = lo.d12(r["nonull_miss_all5"] / exp) if exp > 0 else ""
    else:
        row["nonull_expected_all5"] = row["nonull_ratio_all5"] = ""
    return row


def pooled_keys(first_record):
    """The first pass's pooled reports, by its own rule."""
    out = set()
    for r in lo.read_tsv(os.path.join(first_record, "datasets.tsv")):
        r = dict(r, learner_count=int(r["learner_count"]), resistant=int(r["resistant"]))
        if lo.pooled_in(r):
            out.add((r["species"], r["antibiotic"]))
    return out


def nonull_ratios(rows):
    return [float(r["nonull_ratio_all5"]) for r in rows if r["nonull_ratio_all5"] != ""]


def pool(rows, firsts, scope, species=""):
    out = {"scope": scope, "species": species, "datasets": len(rows)}
    for c in ["genomes", "resistant", "susceptible"] + COUNT_COLUMNS:
        out[c] = sum(int(r[c]) for r in rows)
    out["nonull_expected_all5"] = lo.d12(sum(float(r["nonull_expected_all5"]) for r in rows))
    out["mean_share_median"] = lo.d12(lo.median_of([int(r["mean_miss"]) / int(r["resistant"]) for r in rows]))
    ratios = nonull_ratios(rows)
    out["nonull_datasets_with_ratio"] = len(ratios)
    if ratios:
        q1, q3 = lo.quartiles(sorted(ratios))
        out.update({"nonull_ratio_median": lo.d12(lo.median_of(ratios)), "nonull_ratio_q1": lo.d12(q1),
                    "nonull_ratio_q3": lo.d12(q3), "nonull_datasets_ratio_over_1": sum(1 for x in ratios if x > 1)})
    else:
        out.update({"nonull_ratio_median": "", "nonull_ratio_q1": "", "nonull_ratio_q3": "",
                    "nonull_datasets_ratio_over_1": 0})
    paired = [(float(firsts[(r["species"], r["antibiotic"])]["ratio_all5"]), float(r["nonull_ratio_all5"]))
              for r in rows if r["nonull_ratio_all5"] != "" and firsts[(r["species"], r["antibiotic"])]["ratio_all5"] != ""]
    out["paired_datasets"] = len(paired)
    out["paired_ratio_median"] = lo.d12(lo.median_of([a for a, _ in paired])) if paired else ""
    out["paired_nonull_ratio_median"] = lo.d12(lo.median_of([b for _, b in paired])) if paired else ""
    return out


def pooled_rows(rows, first_record):
    keys = pooled_keys(first_record)
    firsts = first_rows(first_record)
    kept = [r for r in rows if (r["species"], r["antibiotic"]) in keys and r["mean_miss"] != ""]
    if not kept:
        raise Stop("no report enters the pooled rows")
    out = [pool(kept, firsts, "all")]
    for sp in sorted({r["species"] for r in kept}):
        out.append(pool([r for r in kept if r["species"] == sp], firsts, "species", sp))
    return out


def write_record(reading_dir, first_record, out):
    with open(os.path.join(reading_dir, READING), encoding="utf-8") as f:
        reading = json.load(f)
    if (reading.get("cutoff"), reading.get("scheme"), reading.get("five"), reading.get("vote_sum")) != (
            lo.CUTOFF, lo.SCHEME, list(FIVE), VOTE_SUM):
        raise Stop("the reading was made under other rules")
    lo.fresh_out(out)
    rows = [vote_row(r) for r in reading["datasets"]]
    pooled = pooled_rows(rows, first_record)
    rule = lo.reading_rule(float(pooled[0]["nonull_ratio_median"])) if pooled[0]["nonull_ratio_median"] else "none shown"
    parent = os.path.dirname(os.path.abspath(out))
    tmp = tempfile.mkdtemp(prefix=".overlap2-", dir=parent)
    try:
        lo.write_tsv(os.path.join(tmp, "vote.tsv"), VOTE_COLUMNS, rows)
        lo.write_tsv(os.path.join(tmp, "vote-pooled.tsv"), POOLED_COLUMNS, pooled)
        scanned = lo.files_under(reading_dir) + [os.path.join(tmp, n) for n in ("vote.tsv", "vote-pooled.tsv")]
        hits = scan(scanned)
        line6 = (f"second check 6 ids: {len(scanned) + 2} files scanned, {hits} id-shaped tokens or barred keys: "
                 + ("pass" if hits == 0 else "STOP"))
        with open(os.path.join(tmp, "gates.txt"), "w", encoding="utf-8") as f:
            f.write("\n".join(reading["first_checks"] + reading["checks"] + [line6]) + "\n")
            f.write(f"reading rule outside null-fit folds: {rule}\n")
        with open(os.path.join(tmp, "SHA256SUMS"), "w", encoding="utf-8") as f:
            for name in RECORD_FILES:
                f.write(f"{lo.sha256(os.path.join(tmp, name))}  {name}\n")
        hits += scan([os.path.join(tmp, "gates.txt"), os.path.join(tmp, "SHA256SUMS")])
        if hits:
            raise Stop(f"second check 6 ids: {hits} id-shaped tokens or barred keys: STOP")
        lo.publish(tmp, out)
    finally:
        if os.path.isdir(tmp):
            shutil.rmtree(tmp)
    return rows, pooled, line6, rule


def write_reading(reading, out):
    lo.fresh_out(out)
    parent = os.path.dirname(os.path.abspath(out))
    tmp = tempfile.mkdtemp(prefix=".overlap2-", dir=parent)
    try:
        with open(os.path.join(tmp, READING), "w", encoding="utf-8") as f:
            json.dump(reading, f, indent=1, sort_keys=True)
            f.write("\n")
        hits = scan([os.path.join(tmp, READING)])
        if hits:
            raise Stop(f"second check 6 ids: {hits} id-shaped tokens or barred keys in the reading: STOP")
        lo.publish(tmp, out)
    finally:
        if os.path.isdir(tmp):
            shutil.rmtree(tmp)


# ---- the page's tables --------------------------------------------------------------------------------------------

def tables(record, first_record):
    a = lo.read_tsv(os.path.join(record, "vote-pooled.tsv"))
    f = lo.read_tsv(os.path.join(first_record, "pooled.tsv"))[0]
    p, c, d = a[0], lo.count, lo.dec
    out = {}
    out["mean vote"] = ["| Five-learner mean at 0.5 | Genomes |", "|---|---:|",
                        f"| resistant genomes | {c(p['resistant'])} |",
                        f"| missed by the mean | {c(p['mean_miss'])} |",
                        f"| of those, missed by all five | {c(p['mean_miss_by_all5'])} |",
                        f"| of those, missed by 1 to 4 | {c(p['mean_miss_by_1to4'])} |",
                        f"| of those, missed by none | {c(p['mean_miss_by_none'])} |",
                        f"| median share of a report's resistant genomes the mean misses | "
                        f"{d(p['mean_share_median'], 3)} |",
                        f"| susceptible genomes | {c(p['susceptible'])} |",
                        f"| called resistant by the mean | {c(p['mean_wrong'])} |",
                        f"| of those, called resistant by all five | {c(p['mean_wrong_by_all5'])} |",
                        f"| of those, called resistant by 1 to 4 | {c(p['mean_wrong_by_1to4'])} |",
                        f"| of those, called resistant by none | {c(p['mean_wrong_by_none'])} |"]
    out["null-fit cut"] = ["| Measure | All genomes | Outside null-fit folds |", "|---|---:|---:|",
                           f"| resistant genomes | {c(f['resistant'])} | {c(p['nonull_resistant'])} |",
                           f"| missed by all five | {c(f['miss_all5'])} | {c(p['nonull_miss_all5'])} |",
                           f"| yardstick, summed | {d(f['expected_all5'], 2)} | {d(p['nonull_expected_all5'], 2)} |",
                           f"| reports with a ratio | {c(f['datasets_with_ratio'])} "
                           f"| {c(p['nonull_datasets_with_ratio'])} |",
                           f"| median ratio | {d(f['ratio_median'], 1)} | {d(p['nonull_ratio_median'], 1)} |",
                           f"| quartiles of the ratio | {d(f['ratio_q1'], 1)} to {d(f['ratio_q3'], 1)} "
                           f"| {d(p['nonull_ratio_q1'], 1)} to {d(p['nonull_ratio_q3'], 1)} |",
                           f"| reports with a ratio over 1 | {c(f['datasets_ratio_over_1'])} "
                           f"| {c(p['nonull_datasets_ratio_over_1'])} |",
                           f"| reports with both ratios | {c(p['paired_datasets'])} | {c(p['paired_datasets'])} |",
                           f"| median ratio over those | {d(p['paired_ratio_median'], 1)} "
                           f"| {d(p['paired_nonull_ratio_median'], 1)} |"]
    out["mean vote by species"] = (
        ["| Species | Reports | Resistant | Missed by the mean | Of those, by all five | Of those, by 1 to 4 "
         "| Median ratio outside null-fit folds |", "|---|---:|---:|---:|---:|---:|---:|"]
        + [f"| *{s['species']}* | {c(s['datasets'])} | {c(s['resistant'])} | {c(s['mean_miss'])} "
           f"| {c(s['mean_miss_by_all5'])} | {c(s['mean_miss_by_1to4'])} | {d(s['nonull_ratio_median'], 1)} |"
           for s in a[1:]])
    return out


# ---- the command line ---------------------------------------------------------------------------------------------

def plan():
    return {"cutoff": lo.CUTOFF, "vote": "math.fsum(five scores) >= 2.5", "vote_sum": VOTE_SUM, "near": NEAR,
            "scheme": lo.SCHEME, "five": list(FIVE), "first_columns": FIRST_COLUMNS, "vote_columns": VOTE_COLUMNS,
            "pooled_columns": POOLED_COLUMNS}


def main(argv=None):
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("plan")
    r = sub.add_parser("read")
    r.add_argument("reports")
    r.add_argument("--first", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--expect", type=int, default=lo.EXPECT_FILES)
    w = sub.add_parser("write")
    w.add_argument("reading")
    w.add_argument("--first", required=True)
    w.add_argument("--out", required=True)
    p = sub.add_parser("print")
    p.add_argument("record")
    p.add_argument("--first", required=True)
    args = ap.parse_args(argv)
    try:
        if args.cmd == "plan":
            print(json.dumps(plan(), indent=1))
        elif args.cmd == "read":
            lo.fresh_out(args.out)
            reading = read_second(args.reports, args.first, args.expect)
            write_reading(reading, args.out)
            print("\n".join(reading["first_checks"] + reading["checks"]))
            print(f"read: {len(reading['datasets'])} reports read again; wrote {READING}")
        elif args.cmd == "write":
            rows, pooled, line6, rule = write_record(args.reading, args.first, args.out)
            print(line6)
            print(f"wrote vote.tsv ({len(rows)} rows), vote-pooled.tsv ({len(pooled)} rows), gates.txt, SHA256SUMS")
            print(f"pooled: {pooled[0]['datasets']} reports, median ratio outside null-fit folds "
                  f"{pooled[0]['nonull_ratio_median'] or 'none'}; reading rule: {rule}")
        else:
            for name, lines in tables(args.record, args.first).items():
                print(f"## {name}\n" + "\n".join(lines) + "\n")
    except Stop as e:
        print(str(e))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
