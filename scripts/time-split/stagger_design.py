#!/usr/bin/env python3
"""The staggered split by time: cut points, windows, test sets, arms, pairs and fold files, per dataset and clock.

It applies the frozen plan of 25 Nov 2026 (sha256 PLAN_SHA256 below) with the readings of ruling R33 of 26 Nov
2026. It reads metadata alone: the data root's curated datasets and labels, its species lists and
bvbrc/genomes.jsonl, each genome's BV-BRC `date_inserted` and `bioproject_accession` (the batch-*.json files of a
fetch directory), and, once a species is built, its matrix's genome order and kset digests. It reads no score and
no result, and fits nothing.

  python3 -I stagger_design.py count   <data root> <fetch dir> <built root> <out dir>
  python3 -I stagger_design.py folds   <data root> <fetch dir> <built root> <out dir>
  python3 -I stagger_design.py check60 <data root> <fetch dir> <built root> <datasets.tsv> <out dir>

count: for each clock, share set, test minimum, contrast and horizon, the datasets and species with a pair, in
<out dir>/counts.tsv, and the share set and minimum the plan's rule takes for each clock, step by step, in
<out dir>/choice.tsv.

folds: the same count, then, for the share set and minimum taken on each clock, one fold file per dataset, origin,
arm and window that a pair uses, in the form `grm learn --fold-file` reads (docs/fold-file.md), under
<out dir>/fold-files/, each named with its sha256 in <out dir>/design.json. Both fits of a pair test the same
genomes. recent(1) is all(1), so one file serves both arms at origin 1, and the design says so.

check60: the design at the single share 60% on the public-day clock, with the rules of 14 Nov, against the cut day
of each dataset of the 14 Nov record's datasets.tsv. It stops, writing nothing and naming no cut day, unless the
table has 24 rows, each a different dataset and each a dataset of the metadata, matched letter for letter (ruling
R38; the 24 is CHECK_DATASETS, and no option changes it). A row with fewer cells than the header stops it too, with
one line (R41). It writes <out dir>/check60.tsv and exits 1 unless every cut day is the record's.

Every path is checked before it is opened: a path that holds `latetest`, `sealed` or `held-out` stops the script,
which names the path's role and that word and nothing more. An out dir that holds a file stops it too.

The rules, from the plan's Design items 1 to 8 and R33 points 1 to 6:

- Cut point at share s of a dated cohort of n genomes: the date of the genome at place ceil(s x n / 100) in date
  order (R33 point 1). A variant's origins are its shares, whatever its cut points do (R33 point 2).
- Window j runs from after c(j) to c(j+1), the last to the end. Its test set, for j of 2 or more: every cohort
  genome in it with a recorded study that no genome of the species list carries by c(j). A genome with no study is
  tested nowhere. An undated genome, and on the collection-year clock a genome with no four-digit year, is in no
  window, no training set and no test set, and makes no study known.
- Copies, per window (R33 point 3): a genome whose kset digest equals a test genome's of window j is in no training
  set scored on window j.
- all(k), for window j: every cohort genome dated by c(k), less window j's copies. recent(k): per class, the latest
  of all(k), as many as window j's all(1) holds; a tie in date goes to the smaller keyed hash
  sha256("stagger-20261125|" + clock + "|" + dataset + "|" + genome) (R33 point 4, the plan's item 6).
- A pair: origins k < j on window j whose test set holds the test minimum of each class and whose two training sets
  each hold 25 of each class (R33 point 5). Two training sets of the same genomes make no pair: two shares that
  give one cut point leave nothing to compare (a reading the plan and R33 leave open; see the hand-in).
- The plan's rule takes, per clock, the variant with the most origins among those with 12 datasets in 3 species
  at horizon 1 on stale all (R33 point 6), then the highest test minimum, then the most datasets, then the first in
  the plan's order; if none reaches 12 in 3, the one with the most datasets, and the clock is descriptive only.
"""

import hashlib
import json
import re
import sys
from pathlib import Path

PLAN_SHA256 = "ceac8b5bf286fb5e94f05dad69d6ade622ad6234ee0a6fd0ac2469dfc182e39c"
KEY = "stagger-20261125"
CLOCKS = ("public-day", "collection-year")
SHARE_SETS = ((40, 60, 80), (50, 70), (40, 60), (50, 75))
TEST_MINIMUMS = (20, 15, 10)
MIN_TRAIN_CLASS = 25
CONTRASTS = ("stale-all", "stale-same-size")
ARM_OF = {"stale-all": "all", "stale-same-size": "recent"}
CHOICE_CONTRAST, CHOICE_HORIZON = "stale-all", 1
MIN_DATASETS, MIN_SPECIES = 12, 3
CHECK_SHARE = 60
CHECK_DATASETS = 24
BARRED = ("latetest", "sealed", "held-out")
FOLD_FORMAT = "grm-fold-file-v1"
SCHEME = re.compile(r"[a-z0-9-]{1,48}")
DAY = re.compile(r"\d{4}-\d\d-\d\d")
YEAR = re.compile(r"\d{4}")

COUNT_COLUMNS = ["clock", "shares", "test_minimum", "origins", "contrast", "horizon", "datasets", "species", "pairs",
                 "provisional_datasets"]
CHOICE_COLUMNS = ["clock", "step", "rule", "kept", "variants"]
CHECK_COLUMNS = ["dataset", "cohort", "reference_cohort", "undated", "cut_day", "reference_cut_day",
                 "share_public_by_cut", "reference_share_public_by_cut", "same_cut_day"]


class Stop(Exception):
    pass


def need(cond, msg):
    if not cond:
        raise Stop(msg)


# ----------------------------------------------------------------- paths

def guard(path, role):
    """The path, once neither it nor where it leads holds a barred word. A stop names the role and the word alone."""
    p = Path(path)
    texts = [str(p).lower()]
    try:
        texts.append(str(p.resolve()).lower())
    except OSError:
        pass
    for word in BARRED:
        if any(word in t for t in texts):
            raise Stop(f"the {role} path holds '{word}'; nothing under it is opened")
    return p


def read_text(path, role):
    return guard(path, role).read_text(encoding="utf-8")


def read_json(path, role):
    try:
        return json.loads(read_text(path, role))
    except json.JSONDecodeError:
        raise Stop(f"the {role} is not JSON")


def read_jsonl(path, role):
    with open(guard(path, role), encoding="utf-8") as f:
        for line in f:
            if line.strip():
                yield json.loads(line)


def empty_out(path):
    out = guard(path, "out dir")
    need(not out.exists() or (out.is_dir() and not any(out.iterdir())), "the out dir holds a file; it is not written")
    return out


# ----------------------------------------------------------------- metadata

def value_text(v):
    """A genome field as grm reads it (src/bvbrc.rs, value_text): trimmed text, a number as its text, None if blank."""
    if v is None:
        return None
    if isinstance(v, str):
        return v.strip() or None
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, list):
        parts = [p for p in (value_text(x) for x in v) if p is not None]
        return ";".join(parts) if parts else None
    return repr(v) if isinstance(v, float) else str(v)


def year_of(v):
    """The collection year: the field's text when it is four digits, else None."""
    text = value_text(v)
    return int(text) if text is not None and YEAR.fullmatch(text) else None


def reported_label(p):
    """The evaluator's label (docs/fold-file.md): a reported R or S, not a conflict, not inferred from a MIC."""
    if p.get("reported") in ("R", "S") and not p.get("conflict") and not p.get("label_from_mic"):
        return p["reported"]
    return None


def load_fetch(fetch):
    """(day, study) per genome from the fetch's batch files. A missing day is None; a day that is not a date stops."""
    day, study, bad = {}, {}, 0
    folder = guard(fetch, "fetch dir")
    for path in sorted(folder.glob("batch-*.json")):
        for r in read_json(path, "fetch file"):
            g = str(r["genome_id"])
            text = str(r.get("date_inserted") or "").strip()[:10]
            if text and not DAY.fullmatch(text):
                bad += 1
            day[g] = text or None
            s = r.get("bioproject_accession")
            need(s is None or isinstance(s, str), "a fetched study accession is not text")
            study[g] = (s or "").strip().upper() or None
    need(bad == 0, f"{bad} fetched days are not dates")
    return day, study


def load_inputs(root, fetch, built):
    """Everything count, folds and check60 read, each species once."""
    root, built = guard(root, "data root"), guard(built, "built root")
    specs = read_json(root / "curated" / "datasets.json", "dataset list")
    labels = {}
    for p in read_jsonl(root / "curated" / "phenotypes.jsonl", "label table"):
        label = reported_label(p)
        if label:
            labels.setdefault((p["species"], p["antibiotic"]), {})[str(p["genome_id"])] = label
    day, study = load_fetch(fetch)
    year = {}
    for g in read_jsonl(root / "bvbrc" / "genomes.jsonl", "genome table"):
        year[str(g["genome_id"])] = year_of(g.get("collection_year"))
    species = {}
    for slug in sorted({s["species_slug"] for s in specs}):
        text = read_text(root / "species" / slug / "genomes.txt", "species list")
        listed = [x.strip() for x in text.splitlines() if x.strip()]
        matrix = guard(built / "species" / slug / "matrix" / "matrix.json", "matrix")
        if matrix.exists():
            m = read_json(matrix, "matrix")
            need("kset_digests" in m and len(m["kset_digests"]) == len(m["genomes"]),
                 f"{slug}: the matrix has no kset digest per genome")
            genomes = [str(g) for g in m["genomes"]]
            species[slug] = {"genomes": genomes, "digests": dict(zip(genomes, m["kset_digests"])), "built": True,
                             "listed": listed}
        else:
            species[slug] = {"genomes": listed, "digests": None, "built": False, "listed": listed}
    return {"specs": specs, "labels": labels, "dates": {"public-day": day, "collection-year": year},
            "study": study, "species": species}


def cohort_of(inputs, spec):
    """The dataset's cohort as `grm learn` evaluates it: the species' genomes with the evaluator's label."""
    labels = inputs["labels"].get((spec["species"], spec["antibiotic"]), {})
    sp = inputs["species"][spec["species_slug"]]
    cohort = [g for g in sp["genomes"] if g in labels]
    return cohort, labels, sp


# ----------------------------------------------------------------- the rules

def tie_key(clock, dataset, genome):
    """The plan's keyed hash for a tie in recent(k), item 6, to the letter."""
    return hashlib.sha256(("stagger-20261125|" + clock + "|" + dataset + "|" + genome).encode()).hexdigest()


def cut_points(dates, shares):
    """R33 point 1: per share s, the date at place ceil(s x n / 100) of the n dates in order, in whole numbers."""
    order = sorted(dates)
    n = len(order)
    return [order[-(-n * s // 100) - 1] for s in shares]


def window_of(date, cuts):
    """The window of a date, from 1: after c(j) up to c(j+1), the last to the end. None for undated or by c(1)."""
    if date is None or date <= cuts[0]:
        return None
    for j in range(1, len(cuts)):
        if date <= cuts[j]:
            return j
    return len(cuts)


def counts(ids, labels):
    r = sum(labels[g] == "R" for g in ids)
    return r, len(ids) - r


def recent(pool, size, labels, date, key):
    """R33 point 4: per class, the latest size[c] genomes of the pool; a tie in date goes to the smaller key."""
    chosen = set()
    for c in "RS":
        members = sorted((g for g in pool if labels[g] == c), key=key)
        members.sort(key=lambda g: date[g], reverse=True)
        chosen.update(members[:size[c]])
    return [g for g in pool if g in chosen]


def plan_dataset(dataset, clock, shares, cohort, labels, date, study, listed, digests):
    """One dataset, clock and share set: its cut points and, per window from 2, the test set, the copies and the
    arms all(k) and recent(k) for k = 1 to j, each in cohort order. No test minimum enters here."""
    dated = [g for g in cohort if date.get(g) is not None]
    out = {"dated": len(dated), "undated": len(cohort) - len(dated), "cuts": None, "windows": {}}
    if not dated:
        return out
    cuts = cut_points([date[g] for g in dated], shares)
    out["cuts"] = cuts
    key = {}
    for j in range(2, len(shares) + 1):
        # A study is known by c(j) when a dated genome of the species list carries it by then.
        known = {study[g] for g in listed if study.get(g) and date.get(g) is not None and date[g] <= cuts[j - 1]}
        test = [g for g in cohort if window_of(date.get(g), cuts) == j and study.get(g) and study[g] not in known]
        tested = set(test)
        copies = set()
        if digests is not None:
            seen = {digests[g] for g in test}
            copies = {g for g in cohort if g not in tested and digests[g] in seen}
        arms = {}
        for k in range(1, j + 1):
            arms[("all", k)] = [g for g in dated if date[g] <= cuts[k - 1] and g not in copies]
        size = dict(zip("RS", counts(arms[("all", 1)], labels)))
        for k in range(1, j + 1):
            pool = arms[("all", k)]
            for g in pool:
                if g not in key:
                    key[g] = tie_key(clock, dataset, g)
            arms[("recent", k)] = recent(pool, size, labels, date, key.__getitem__)
        out["windows"][j] = {"test": test, "copies": copies, "arms": arms}
    return out


def pairs_of(plan, labels, test_minimum, contrast):
    """R33 point 5: the pairs (k, j) of one contrast. Two training sets of the same genomes make no pair."""
    arm = ARM_OF[contrast]
    out = []
    for j, w in sorted(plan["windows"].items()):
        if min(counts(w["test"], labels)) < test_minimum:
            continue
        fresh = w["arms"][(arm, j)]
        if min(counts(fresh, labels)) < MIN_TRAIN_CLASS:
            continue
        for k in range(1, j):
            older = w["arms"][(arm, k)]
            if min(counts(older, labels)) < MIN_TRAIN_CLASS:
                continue
            if set(older) == set(fresh):
                continue
            out.append((k, j))
    return out


def variant_name(shares, minimum):
    return "-".join(map(str, shares)) + "/" + str(minimum)


def variants():
    """The variants in the plan's order: share sets as listed, each with the test minimums 20, 15 and 10."""
    return [(shares, m) for shares in SHARE_SETS for m in TEST_MINIMUMS]


def choose(entries):
    """The plan's rule over one clock's variants, given in the plan's order as dicts with shares, minimum, datasets
    and species at horizon 1 on stale all. (taken, descriptive only, steps); each step names what it keeps."""
    steps = []
    kept = [e for e in entries if e["datasets"] >= MIN_DATASETS and e["species"] >= MIN_SPECIES]
    steps.append(("1", f"at least {MIN_DATASETS} datasets in {MIN_SPECIES} species at horizon 1, stale all", kept))
    descriptive = not kept
    if kept:
        most = max(len(e["shares"]) for e in kept)
        kept = [e for e in kept if len(e["shares"]) == most]
        steps.append(("2", "the most origins", kept))
        highest = max(e["minimum"] for e in kept)
        kept = [e for e in kept if e["minimum"] == highest]
        steps.append(("3", "the highest test minimum", kept))
        most = max(e["datasets"] for e in kept)
        kept = [e for e in kept if e["datasets"] == most]
        steps.append(("4", "the most datasets", kept))
    else:
        most = max(e["datasets"] for e in entries)
        kept = [e for e in entries if e["datasets"] == most]
        steps.append(("2", f"none reaches {MIN_DATASETS} in {MIN_SPECIES}: the most datasets; descriptive only", kept))
    taken = kept[0]
    steps.append((str(len(steps) + 1), "a tie goes to the first in the plan's order", [taken]))
    return taken, descriptive, steps


# ----------------------------------------------------------------- count

def plan_all(inputs, clock, shares):
    """plan_dataset for every dataset with a cohort, in the dataset list's order."""
    out = []
    for spec in inputs["specs"]:
        cohort, labels, sp = cohort_of(inputs, spec)
        if not cohort:
            continue
        plan = plan_dataset(spec["id"], clock, shares, cohort, labels, inputs["dates"][clock], inputs["study"],
                            sp["listed"], sp["digests"])
        out.append((spec, cohort, labels, sp, plan))
    return out


def count(inputs):
    """(count rows, choice per clock, plans per clock and share set)."""
    rows, choices, plans = [], {}, {}
    for clock in CLOCKS:
        for shares in SHARE_SETS:
            plans[(clock, shares)] = plan_all(inputs, clock, shares)
        for shares, minimum in variants():
            for contrast in CONTRASTS:
                by_h = {h: [] for h in range(1, len(shares))}
                for spec, _, labels, sp, plan in plans[(clock, shares)]:
                    found = pairs_of(plan, labels, minimum, contrast)
                    for h in by_h:
                        n = sum(1 for k, j in found if j - k == h)
                        if n:
                            by_h[h].append((spec, sp, n))
                for h, hits in by_h.items():
                    row = {"clock": clock, "shares": "-".join(map(str, shares)), "test_minimum": minimum,
                           "origins": len(shares), "contrast": contrast, "horizon": h, "datasets": len(hits),
                           "species": len({s["species_slug"] for s, _, _ in hits}),
                           "pairs": sum(n for _, _, n in hits),
                           "provisional_datasets": sum(1 for _, sp, _ in hits if not sp["built"])}
                    rows.append(row)
        choices[clock] = choose(choice_entries(rows, clock))
    return rows, choices, plans


def choice_entries(rows, clock):
    """R33 point 6: what the rule reads, each variant's datasets and species on stale all at horizon 1, in the
    plan's order. The same-size contrast's rows stand beside and choose nothing."""
    by = {(r["shares"], r["test_minimum"]): r for r in rows
          if r["clock"] == clock and r["contrast"] == CHOICE_CONTRAST and r["horizon"] == CHOICE_HORIZON}
    out = []
    for shares, minimum in variants():
        r = by[("-".join(map(str, shares)), minimum)]
        out.append({"shares": shares, "minimum": minimum, "datasets": r["datasets"], "species": r["species"]})
    return out


def tsv(columns, rows):
    return "\t".join(columns) + "\n" + "".join("\t".join(str(r[c]) for c in columns) + "\n" for r in rows)


def choice_rows(choices):
    rows = []
    for clock in CLOCKS:
        _, _, steps = choices[clock]
        for step, rule, kept in steps:
            rows.append({"clock": clock, "step": step, "rule": rule, "kept": len(kept),
                         "variants": "; ".join(variant_name(e["shares"], e["minimum"]) for e in kept)})
    return rows


def choice_lines(choices):
    lines = []
    for clock in CLOCKS:
        taken, descriptive, steps = choices[clock]
        for step, rule, kept in steps:
            lines.append(f"{clock} step {step}, {rule}: {len(kept)} kept"
                         + (": " + ", ".join(variant_name(e["shares"], e["minimum"]) for e in kept) if kept else ""))
        lines.append(f"{clock}: the rule takes shares {'-'.join(map(str, taken['shares']))} with test minimum "
                     f"{taken['minimum']}, {taken['datasets']} datasets in {taken['species']} species at horizon 1, "
                     f"stale all" + ("; descriptive only" if descriptive else ""))
    return lines


# ----------------------------------------------------------------- folds

def scheme_of(clock, j, k, arm):
    s = f"stagger-{clock}-w{j}-o{k}-{arm}"
    need(SCHEME.fullmatch(s) is not None, f"{s}: not a scheme name")
    return s


def fold_body(dataset, scheme, cohort, test, train):
    tested, trained = set(test), set(train)
    need(not tested & trained, f"{dataset} {scheme}: a genome both tests and trains")
    return {"format": FOLD_FORMAT, "dataset": dataset, "scheme": scheme, "test": [list(test)],
            "train": list(train), "unused": [g for g in cohort if g not in tested and g not in trained]}


def fold_files(inputs, choices, plans):
    """(design, {file name: bytes}, lines) for the share set and minimum taken on each clock."""
    design = {"plan_sha256": PLAN_SHA256, "key": KEY, "min_train_class": MIN_TRAIN_CLASS, "clocks": {},
              "datasets": []}
    files, lines = {}, []
    for clock in CLOCKS:
        taken, descriptive, steps = choices[clock]
        shares, minimum = taken["shares"], taken["minimum"]
        design["clocks"][clock] = {"shares": list(shares), "test_minimum": minimum, "descriptive_only": descriptive,
                                   "datasets_at_horizon_1": taken["datasets"], "species_at_horizon_1": taken["species"]}
        shared = 0
        for spec, cohort, labels, sp, plan in plans[(clock, shares)]:
            row = {"id": spec["id"], "species_slug": spec["species_slug"], "antibiotic": spec["antibiotic"],
                   "clock": clock, "built": sp["built"], "cohort": len(cohort),
                   "cohort_resistant": counts(cohort, labels)[0], "dated": plan["dated"], "undated": plan["undated"],
                   "cuts": plan["cuts"], "windows": [], "pairs": {c: [] for c in CONTRASTS}, "why_not": ""}
            found = {c: pairs_of(plan, labels, minimum, c) for c in CONTRASTS}
            for c in CONTRASTS:
                row["pairs"][c] = [{"older": k, "fresh": j, "window": j, "horizon": j - k} for k, j in found[c]]
            if not any(found.values()):
                row["why_not"] = "no pair at any horizon"
            elif not sp["built"]:
                row["why_not"] = "species not built: its cohort is provisional and its copies are not known"
            design["datasets"].append(row)
            if row["why_not"]:
                continue
            for j, w in sorted(plan["windows"].items()):
                need_fits = set()
                for c in CONTRASTS:
                    for k, jj in found[c]:
                        if jj == j:
                            need_fits.update({(ARM_OF[c], k), (ARM_OF[c], j)})
                if not need_fits:
                    continue
                test_r, test_s = counts(w["test"], labels)
                wrow = {"window": j, "test": len(w["test"]), "test_resistant": test_r, "copies": len(w["copies"]),
                        "fits": []}
                need(set(w["arms"][("recent", 1)]) == set(w["arms"][("all", 1)]),
                     f"{spec['id']} {clock} window {j}: recent(1) is not all(1)")
                written = {}
                for arm, k in sorted(need_fits, key=lambda x: (x[1], x[0])):
                    train = w["arms"][(arm, k)]
                    fit = {"origin": k, "arm": arm, "train": len(train), "train_resistant": counts(train, labels)[0]}
                    base = (arm, k)
                    if k == 1 and arm == "recent":
                        fit["same_file_as"] = "all"
                        fit["note"] = "recent(1) is all(1): one file"
                        base = ("all", 1)
                        shared += 1
                    if base not in written:
                        scheme = scheme_of(clock, j, k, base[0])
                        body = fold_body(spec["id"], scheme, cohort, w["test"], w["arms"][base])
                        data = (json.dumps(body, separators=(",", ":")) + "\n").encode()
                        name = f"{spec['id']}.{scheme}.json"
                        need(name not in files, f"{name}: written twice")
                        files[name] = data
                        written[base] = {"scheme": scheme, "file": name, "sha256": hashlib.sha256(data).hexdigest()}
                    fit.update(written[base])
                    wrow["fits"].append(fit)
                row["windows"].append(wrow)
            lines.append(f"{clock} {spec['id']}: cuts {', '.join(map(str, plan['cuts']))}; "
                         + "; ".join(f"window {w['window']} tests {w['test']} ({w['test_resistant']} R), "
                                     f"{len(w['fits'])} fits" for w in row["windows"])
                         + "; pairs " + ", ".join(f"{c} {len(found[c])}" for c in CONTRASTS))
        lines.append(f"{clock}: {sum(1 for r in design['datasets'] if r['clock'] == clock and not r['why_not'])} "
                     f"datasets with fold files; {shared} files serve both all(1) and recent(1)")
    lines.append(f"fold files: {len(files)}")
    return design, files, lines


# ----------------------------------------------------------------- check60

def check60(inputs, table_path, datasets=CHECK_DATASETS):
    """Rows of the 60% check against the 14 Nov record's datasets.tsv, and whether every cut day is the record's.

    It stops before it cuts any day unless the table has `datasets` rows, each a different dataset and each a
    dataset of the metadata. main() gives the plan's 24 and nothing else; a test passes its own count."""
    text = read_text(table_path, "reference table")
    lines = text.rstrip("\n").split("\n")
    head = lines[0].split("\t")
    cells = [line.split("\t") for line in lines[1:]]
    reference = [dict(zip(head, row)) for row in cells]
    need(reference and {"dataset", "cohort", "cut_day", "share_public_by_cut"} <= set(head),
         "the reference table has no dataset rows of the 14 Nov record's columns")
    ragged = [str(i) for i, row in enumerate(cells, 1) if len(row) < len(head)]
    need(not ragged, f"{len(ragged)} of the reference table's {len(cells)} rows hold fewer cells than its header's "
                     f"{len(head)}: row {', '.join(ragged)}")
    names = [ref["dataset"] for ref in reference]
    need(len(reference) == datasets,
         f"the 60% check needs {datasets} dataset rows, each a different dataset of the metadata; the reference "
         f"table has {len(reference)}")
    twice = sorted({n for n in names if names.count(n) > 1})
    need(not twice, f"the 60% check needs {datasets} different datasets; the reference table's {len(reference)} rows "
                    f"name {len(set(names))}, and these more than once: {', '.join(twice)}")
    specs = {s["id"]: s for s in inputs["specs"]}
    lacking = sorted(set(names) - set(specs))
    need(not lacking, f"the 60% check needs {datasets} datasets of the metadata; the metadata lacks "
                      f"{len(lacking)} of the reference table's datasets: {', '.join(lacking)}")
    day = inputs["dates"]["public-day"]
    rows = []
    for ref in reference:
        row = {"dataset": ref["dataset"], "reference_cohort": ref["cohort"], "reference_cut_day": ref["cut_day"],
               "reference_share_public_by_cut": ref["share_public_by_cut"], "cohort": 0, "undated": 0,
               "cut_day": "", "share_public_by_cut": "", "same_cut_day": 0}
        spec = specs.get(ref["dataset"])
        if spec is not None:
            cohort, _, _ = cohort_of(inputs, spec)
            dated = [day[g] for g in cohort if day.get(g) is not None]
            row["cohort"], row["undated"] = len(cohort), len(cohort) - len(dated)
            if dated:
                cut = cut_points(dated, (CHECK_SHARE,))[0]
                row["cut_day"] = cut
                row["share_public_by_cut"] = str(round(sum(d <= cut for d in dated) / len(cohort), 4))
                row["same_cut_day"] = int(cut == ref["cut_day"] and row["undated"] == 0)
        rows.append(row)
    return rows, all(r["same_cut_day"] for r in rows)


# ----------------------------------------------------------------- main

def write(out, files):
    out.mkdir(parents=True, exist_ok=True)
    for name, data in files.items():
        path = out / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data if isinstance(data, bytes) else data.encode())


def main(argv):
    usage = __doc__.strip().split("\n\n")[1]
    if not argv or argv[0] not in ("count", "folds", "check60") or len(argv) != (6 if argv[0] == "check60" else 5):
        print(usage, file=sys.stderr)
        return 2
    try:
        mode = argv[0]
        out = empty_out(argv[-1])
        inputs = load_inputs(argv[1], argv[2], argv[3])
        if mode == "check60":
            rows, ok = check60(inputs, argv[4])
            write(out, {"check60.tsv": tsv(CHECK_COLUMNS, rows)})
            for r in rows:
                print(f"{r['dataset']}: cut day {r['cut_day'] or 'none'}, the record's {r['reference_cut_day']}, "
                      f"cohort {r['cohort']} against {r['reference_cohort']}, undated {r['undated']}")
            print(f"check60: {sum(r['same_cut_day'] for r in rows)} of {len(rows)} cut days are the record's")
            return 0 if ok else 1
        rows, choices, plans = count(inputs)
        files = {"counts.tsv": tsv(COUNT_COLUMNS, rows), "choice.tsv": tsv(CHOICE_COLUMNS, choice_rows(choices))}
        lines = choice_lines(choices)
        if mode == "folds":
            design, folds, more = fold_files(inputs, choices, plans)
            files["design.json"] = json.dumps(design, indent=1) + "\n"
            files.update({"fold-files/" + n: d for n, d in folds.items()})
            lines += more
        write(out, files)
        print("\n".join(lines))
        print(f"count: {len(rows)} rows; it read no score and no result")
        return 0
    except Stop as e:
        print(f"stagger_design: stopped: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
