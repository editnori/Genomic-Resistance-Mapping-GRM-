"""revlasso: compare my three-dataset replication with the lane and the release
(PLAN.md measures 1 to 4).

usage: python3 -I rep_compare.py RUNS_DIR LANE_EVENTS_DIR LANE_DUMP_DIR RELEASE_DIR OUT_JSON ID...
"""
import gzip
import json
import math
import sys
from pathlib import Path


def events(path):
    p = Path(path)
    text = p.read_text() if p.exists() else gzip.open(str(p) + ".gz", "rt").read()
    return [json.loads(x) for x in text.splitlines() if x.strip()]


def strip(e):
    return {k: v for k, v in e.items() if k != "seconds"}


def first_diff(a, b, path=""):
    if type(a) is not type(b):
        return f"{path}: type {type(a).__name__} vs {type(b).__name__}"
    if isinstance(a, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                return f"{path}.{k}: missing on one side"
            r = first_diff(a[k], b[k], f"{path}.{k}")
            if r:
                return r
        return None
    if isinstance(a, list):
        if len(a) != len(b):
            return f"{path}: length {len(a)} vs {len(b)}"
        for i, (x, y) in enumerate(zip(a, b)):
            r = first_diff(x, y, f"{path}[{i}]")
            if r:
                return r
        return None
    return None if a == b else f"{path}: {a!r} vs {b!r}"


def main(argv):
    runs, lane_ev, lane_dump, rel, outp = (Path(a) for a in argv[:5])
    ids = argv[5:]
    out = []
    for did in ids:
        r = json.loads(Path(rel, "datasets", did + ".json").read_text())
        y = [bool(s["resistant"]) for s in r["samples"]]
        mine = events(runs / "plain" / f"{did}.events.jsonl")
        mdump = events(runs / "dump" / f"{did}.events.jsonl")
        lane = events(lane_ev / f"{did}.events.jsonl")
        inp = [e for e in mine if e["event"] == "inputs"][0]
        row = {"id": did,
               "inputs": {k: inp[k] for k in ("input_fingerprint", "samples", "resistant", "features")},
               "release": {"input_fingerprint": r["input_fingerprint"], "samples": r["counts"]["samples"],
                           "resistant": r["counts"]["resistant"], "features": r["counts"]["features"]},
               "events_mine": len(mine), "events_lane": len(lane)}
        row["inputs_equal_release"] = all(row["inputs"][k] == row["release"][k] for k in row["inputs"])
        lane_inp = [e for e in lane if e["event"] == "inputs"][0]
        row["inputs_equal_lane"] = strip(inp) == strip(lane_inp)
        sm, sl, sd = [strip(e) for e in mine], [strip(e) for e in lane], [strip(e) for e in mdump]
        row["events_equal_lane"] = sm == sl
        row["events_first_diff_vs_lane"] = None if sm == sl else first_diff(sm, sl)
        row["dump_binary_events_equal_plain"] = sm == sd
        row["dump_binary_first_diff"] = None if sm == sd else first_diff(sm, sd)
        splits = []
        part_diffs = []
        for sp in [s for s in r["splits"] if not s.get("unavailable_reason")]:
            L = [l for l in sp["learners"] if l["learner"] == "lasso"][0]
            n = len(y)
            oof = [None] * n
            same_scores = 0
            compared = 0
            for fi in range(sp["folds"]):
                name = f"{did}__{sp['scheme']}__{fi}.json"
                a = json.loads((runs / "dumps" / name).read_text())
                b = json.loads(Path(lane_dump, name).read_text())
                for k in ("train", "fit", "val", "y", "features", "n", "seed", "audit", "complexity", "path", "refit"):
                    if a[k] != b[k]:
                        part_diffs.append(f"{sp['scheme']} {fi} {k}")
                if (a["scores"] is None) != (b["scores"] is None):
                    part_diffs.append(f"{sp['scheme']} {fi} scores presence")
                tr = set(a["train"])
                test = [i for i in range(n) if i not in tr]
                if test != [i for i in range(n) if sp["test_fold"][i] == fi]:
                    part_diffs.append(f"{sp['scheme']} {fi} test fold vs release")
                if set(a["y"]) != {i for i in range(n) if y[i]}:
                    part_diffs.append(f"{sp['scheme']} {fi} labels vs release")
                if a["scores"] is None:
                    prev = sum(1 for i in a["train"] if y[i]) / len(a["train"])
                    for i in test:
                        oof[i] = prev
                else:
                    for i in test:
                        oof[i] = a["scores"][i]
                        compared += 1
                        same_scores += int(b["scores"] is not None and a["scores"][i] == b["scores"][i])
            rs = L["scores"]
            splits.append({"scheme": sp["scheme"], "n": n,
                           "binary_scores_compared_with_lane": compared, "binary_scores_equal_lane": same_scores,
                           "pooled_equal_release": sum(1 for p, q in zip(oof, rs) if p == q),
                           "calls_changed_vs_release": sum(1 for p, q in zip(oof, rs) if (p >= 0.5) != (q >= 0.5))})
        row["splits"] = splits
        row["partition_and_fit_diffs_vs_lane"] = part_diffs
        fm = [e for e in mine if e.get("event") == "fit" and e["scheme"] == "final"]
        row["final_model_equal_release"] = bool(fm) and fm[0]["model"] == r["models"]["lasso"] and fm[0]["audit"] == r["models"]["lasso_fit"]
        row["seconds_mine_outer"] = sum(e["seconds"] for e in mine if e.get("event") == "fit" and e["scheme"] != "final")
        out.append(row)
        print(json.dumps({k: v for k, v in row.items() if k not in ("splits",)}, indent=None))
        for s in splits:
            print("  ", s)
    Path(outp).write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
