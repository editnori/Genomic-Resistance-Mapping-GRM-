#!/usr/bin/env python3
"""Reviewer's own recount of the clindamycin table by study (brief point 2).

Written without reuse of lane clinda's code. Reads the six released reports, the curated phenotypes and every
raw BV-BRC AMR record (no text prefilter: each line is parsed), and writes under results/counts/:

  by-key.tsv        per dataset and study key, under two study fields:
                      label   = PMIDs of the raw records whose call equals the sample's label (the plan's field)
                      all     = PMIDs of every raw record of the pair (what src/curate.rs keeps as pair PMIDs)
  by-pmid.tsv       per dataset and PMID, membership counting (a sample with two PMIDs counts under both),
                    exclusive counting (only samples whose key is that PMID alone) and fractional counting
                    (a sample with k PMIDs counts 1/k)
  chosen.tsv        the plan's choice rule under each field and counting, with shares and exact intervals
  multi.tsv         samples with more than one PMID, per dataset
  checks.txt        input checks

Run: python3 -I scripts/recount.py
"""
import collections
import json
import math
import os

RELEASE = "/home/user/grm-ops/owner2/work/rev14-full/rev14-274-20261112"
ROOT = "/home/user/grm-h1raw"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "counts")
DRUG = "clindamycin"


def call_of(text):
    if text is None:
        return None
    t = text.strip().lower()
    if t in ("resistant", "nonsusceptible", "non-susceptible"):
        return "R"
    if t == "susceptible":
        return "S"
    if t in ("intermediate", "susceptible-dose dependent", "reduced susceptibility"):
        return "I"
    return None


def split_ids(v):
    if not v:
        return []
    return [x.strip() for x in str(v).split(";") if x.strip()]


# Exact binomial tail sums in log space.
def log_pmf(k, n, p):
    if p <= 0:
        return 0.0 if k == 0 else -math.inf
    if p >= 1:
        return 0.0 if k == n else -math.inf
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
            + k * math.log(p) + (n - k) * math.log1p(-p))


def cdf(k, n, p):  # P(X <= k)
    return min(1.0, sum(math.exp(log_pmf(i, n, p)) for i in range(0, k + 1)))


def sf(k, n, p):  # P(X >= k)
    return min(1.0, sum(math.exp(log_pmf(i, n, p)) for i in range(k, n + 1)))


def bisect(f, lo, hi, it=200):
    for _ in range(it):
        mid = (lo + hi) / 2
        if f(mid):
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def exact_ci(k, n, a=0.05):
    lo = 0.0 if k == 0 else bisect(lambda p: sf(k, n, p) >= a / 2, 0.0, 1.0)
    hi = 1.0 if k == n else bisect(lambda p: cdf(k, n, p) <= a / 2, 0.0, 1.0)
    return lo, hi


def fisher(a, b, c, d):
    r1, r2, c1, n = a + b, c + d, a + c, a + b + c + d

    def lp(x):
        return (math.lgamma(r1 + 1) - math.lgamma(x + 1) - math.lgamma(r1 - x + 1)
                + math.lgamma(r2 + 1) - math.lgamma(c1 - x + 1) - math.lgamma(r2 - c1 + x + 1)
                - math.lgamma(n + 1) + math.lgamma(c1 + 1) + math.lgamma(n - c1 + 1))
    lo, hi = max(0, c1 - r2), min(r1, c1)
    p0 = lp(a)
    return min(1.0, sum(math.exp(lp(x)) for x in range(lo, hi + 1) if lp(x) <= p0 + 1e-7))


def main():
    os.makedirs(OUT, exist_ok=True)
    checks = []
    idx = json.load(open(os.path.join(RELEASE, "index.json")))
    dsets = sorted(d["id"] for d in idx["datasets"] if d["antibiotic"] == DRUG)
    checks.append(f"index datasets={len(idx['datasets'])} clindamycin={dsets}")
    reps = {d: json.load(open(os.path.join(RELEASE, "datasets", d + ".json"))) for d in dsets}
    want = {s["id"] for r in reps.values() for s in r["samples"]}

    cur = collections.defaultdict(list)
    for line in open(os.path.join(ROOT, "curated", "phenotypes.jsonl")):
        r = json.loads(line)
        if r["antibiotic"] == DRUG and r["genome_id"] in want:
            cur[r["genome_id"]].append(r)

    names = json.load(open(os.path.join(ROOT, "curated", "drug_names.json")))["names"]
    raw_ok = {n["raw"] for n in names if n["canonical"] == DRUG and n["kind"] == "drug"}
    checks.append(f"raw names resolved to the drug: {sorted(raw_ok)}")
    raw = collections.defaultdict(list)
    n_lines = 0
    name_variants = collections.Counter()
    for line in open(os.path.join(ROOT, "bvbrc", "genome_amr.jsonl")):
        n_lines += 1
        r = json.loads(line)
        ab = r.get("antibiotic") or ""
        if "lindamycin" in ab.lower():
            name_variants[ab] += 1
        if ab in raw_ok and r["genome_id"] in want:
            raw[r["genome_id"]].append(r)
    checks.append(f"raw lines parsed={n_lines}; antibiotic names holding 'lindamycin' (any case): "
                  f"{dict(name_variants)}")

    key_rows, pmid_rows, chosen_rows, multi_rows = [], [], [], []
    for d in dsets:
        S = reps[d]["samples"]
        bad = 0
        mem = []
        for s in S:
            lab = "R" if s["resistant"] else "S"
            rule = "R" if s["catalog_resistant"] else "S"
            c = cur.get(s["id"], [])
            if not (len(c) == 1 and c[0]["reported"] == lab and not c[0]["conflict"] and not c[0]["label_from_mic"]):
                bad += 1
            recs = raw.get(s["id"], [])
            lab_p = sorted({p for r in recs if call_of(r.get("resistant_phenotype")) == lab
                            for p in split_ids(r.get("pmid"))}, key=lambda x: (len(x), x))
            all_p = sorted({p for r in recs for p in split_ids(r.get("pmid"))}, key=lambda x: (len(x), x))
            mem.append(dict(id=s["id"], lab=lab, rule=rule, label=lab_p, all=all_p, lineage=s["lineage"],
                            nrec=len(recs)))
        checks.append(f"{d}: samples={len(S)} curated mismatches={bad} samples with no raw record="
                      f"{sum(1 for m in mem if m['nrec'] == 0)} label-field differs from all-field="
                      f"{sum(1 for m in mem if m['label'] != m['all'])}")

        def tally(ms, w=None):
            c = collections.Counter()
            for m in ms:
                x = 1.0 if w is None else w(m)
                c["n"] += x
                c["rR"] += x * (m["lab"] == "R")
                c["rS"] += x * (m["lab"] == "S")
                c["gR"] += x * (m["rule"] == "R")
                c["gS"] += x * (m["rule"] == "S")
                c["fr"] += x * (m["rule"] == "R" and m["lab"] == "S")
                c["fs"] += x * (m["rule"] == "S" and m["lab"] == "R")
            c["dis"] = c["fr"] + c["fs"]
            return c

        tot = tally(mem)
        for field in ("label", "all"):
            keys = collections.defaultdict(list)
            for m in mem:
                keys["+".join(m[field]) or "no-pmid"].append(m)
            for k, ms in sorted(keys.items(), key=lambda kv: (-tally(kv[1])["dis"], -len(kv[1]), kv[0])):
                t = tally(ms)
                key_rows.append([field, d, k] + [int(t[x]) for x in ("n", "rR", "rS", "gR", "gS", "fr", "fs", "dis")])
            key_rows.append([field, d, "TOTAL"] + [int(tot[x]) for x in ("n", "rR", "rS", "gR", "gS", "fr", "fs", "dis")])
            multi = [m for m in mem if len(m[field]) > 1]
            multi_rows.append([field, d, len(multi), sum(1 for m in multi if m["lab"] != m["rule"]),
                               ";".join(sorted({"+".join(m[field]) for m in multi}))])
            pm = sorted({p for m in mem for p in m[field]}, key=lambda x: (len(x), x))
            counts = {}
            for p in pm:
                inm = [m for m in mem if p in m[field]]
                exc = [m for m in mem if m[field] == [p]]
                frac = tally([m for m in mem if p in m[field]], w=lambda m: 1.0 / len(m[field]))
                counts[p] = {"membership": tally(inm), "exclusive": tally(exc), "fractional": frac}
                for how in ("membership", "exclusive", "fractional"):
                    t = counts[p][how]
                    pmid_rows.append([field, how, d, p] + [round(t[x], 3) for x in
                                                           ("n", "rR", "rS", "gR", "gS", "fr", "fs", "dis")])
            for how in ("membership", "exclusive", "fractional"):
                cand = [p for p in pm if counts[p][how]["dis"] > 0]
                if not cand:
                    chosen_rows.append([field, how, d, "none"] + ["NA"] * 15)
                    continue
                cand.sort(key=lambda p: (-counts[p][how]["dis"], -counts[p][how]["fr"], -counts[p][how]["n"],
                                         int(p) if p.isdigit() else 10 ** 12))
                T = cand[0]
                runner = f"{cand[1]}:{counts[cand[1]][how]['dis']:g}" if len(cand) > 1 else "none"
                t = counts[T]["membership"]  # the study's own samples, whole
                k, n = int(t["dis"]), int(tot["dis"])
                lo, hi = exact_ci(k, n)
                rest_dis, rest_n = int(tot["dis"] - t["dis"]), int(tot["n"] - t["n"])
                p_f = fisher(k, int(t["n"]) - k, rest_dis, rest_n - rest_dis)
                left_R, left_S = int(tot["rR"] - t["rR"]), int(tot["rS"] - t["rS"])
                # Shares among samples that carry a PMID.
                with_p = [m for m in mem if m[field]]
                wp = tally(with_p)
                chosen_rows.append([field, how, d, T, runner, k, n, f"{k / n:.4f}", f"{lo:.4f}", f"{hi:.4f}",
                                    int(t["n"]), int(tot["n"]), f"{t['n'] / tot['n']:.4f}",
                                    f"{int(t['fr'])}/{int(tot['fr'])}", f"{int(t['fs'])}/{int(tot['fs'])}",
                                    f"{p_f:.3g}", f"{left_R}/{left_S}",
                                    f"{k}/{int(wp['dis'])}" + (f"={k / wp['dis']:.4f}" if wp["dis"] else ""),
                                    f"{int(t['n'])}/{int(wp['n'])}={t['n'] / wp['n']:.4f}"])

    def write(name, head, rows):
        with open(os.path.join(OUT, name), "w") as fh:
            fh.write("\t".join(head) + "\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")
    cols = ["samples", "reported_R", "reported_S", "rule_R", "rule_S", "false_resistant", "false_susceptible",
            "disagreements"]
    write("by-key.tsv", ["field", "dataset", "study_key"] + cols, key_rows)
    write("by-pmid.tsv", ["field", "counting", "dataset", "pmid"] + cols, pmid_rows)
    write("chosen.tsv", ["field", "counting", "dataset", "chosen", "runner_up", "dis_in_study", "dis_total",
                         "share_dis", "ci_low", "ci_high", "samples_in_study", "samples_total", "share_samples",
                         "fr_in/total", "fs_in/total", "fisher_p", "left_R/left_S",
                         "share_dis_among_samples_with_pmid", "share_samples_among_samples_with_pmid"], chosen_rows)
    write("multi.tsv", ["field", "dataset", "samples_with_2plus_pmids", "of_which_disagree", "keys"], multi_rows)
    with open(os.path.join(OUT, "checks.txt"), "w") as fh:
        fh.write("\n".join(checks) + "\n")
    print("\n".join(checks))


if __name__ == "__main__":
    main()
