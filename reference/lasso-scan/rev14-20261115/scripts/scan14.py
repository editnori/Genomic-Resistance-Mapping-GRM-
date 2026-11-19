"""The large-species lasso scan at revision 14 (lane lassoscan, PLAN.md).

Runs the `lasso-scan` example binary built from main (852ebc9) on one dataset
per call, in the order of PLAN.md, and applies the scan's own per-fit test
(`unconverged` in scripts/lasso-scan/scan.py, LassoFitAudit::has_unconverged).
The input check is against the revision 14 released report's input
fingerprint, not revision 7's. Unlike scan.py it does not stop at the first
unconverged fit: it records where the registered rule would have stopped and
goes on. Each dataset has a 30-minute operational cap.

usage: python3 -I scan14.py --binary BIN --data ROOT --release DIR --list LIST --out OUT
           [--threads 16] [--only ID[,ID...]] [--dump DIR [--dump-full ID|SCHEME|FOLD]]
"""
import argparse
import datetime
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ORDER = ["mycobacterium_tuberculosis", "salmonella_enterica", "staphylococcus_aureus",
         "klebsiella_pneumoniae", "escherichia_coli"]
CAP_SECONDS = 1800


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def unconverged(audit):
    return audit["path_unconverged"] > 0 or any(
        audit[k] is not None and not audit[k]["converged"] for k in ("selected_convergence", "refit_convergence"))


def ordered(datasets):
    out = []
    for s in ORDER:
        out += [d for d in datasets if d["species"] == s]
    assert len(out) == len(datasets)
    return out


def fit_summary(ev):
    tr = ev["trace"]
    pts = tr["path"] + tr["refit"]
    a = ev["audit"]
    return {
        "scheme": ev["scheme"], "outer_fold": ev["outer_fold"], "seconds": ev["seconds"],
        "error": ev.get("error"), "audit": a,
        "converged": ev.get("error") is None and a is not None and not unconverged(a),
        "path_points": len(tr["path"]), "refit_points": len(tr["refit"]),
        "irls_total": sum(p["irls_steps"] for p in pts),
        "irls_max": max((p["irls_steps"] for p in pts), default=0),
        "capped_points": sum(p["sweep_cap_hits"] > 0 for p in pts),
        "weights": None if ev.get("model") is None else len(ev["model"]["coef"]),
    }


def run_dataset(args, entry, out):
    rel = json.loads(Path(args.release, "datasets", entry["id"] + ".json").read_text())
    cmd = [str(args.binary), str(args.data), entry["species"], entry["id"]]
    env = dict(os.environ, RAYON_NUM_THREADS=str(args.threads))
    if args.dump:
        env["GRM_LASSO_DUMP"] = str(args.dump)
        if args.dump_full:
            env["GRM_LASSO_DUMP_FULL"] = args.dump_full
    load0 = Path("/proc/loadavg").read_text().split()[:3]
    started_utc, t0 = now(), time.monotonic()
    with open(out / f"{entry['id']}.events.jsonl", "wb") as so, open(out / f"{entry['id']}.stderr.log", "wb") as se:
        p = subprocess.Popen(cmd, stdout=so, stderr=se, stdin=subprocess.DEVNULL, env=env)
        try:
            rc = p.wait(timeout=CAP_SECONDS)
        except subprocess.TimeoutExpired:
            p.kill()
            p.wait()
            rc = "killed at the 30-minute cap"
    wall = time.monotonic() - t0
    load1 = Path("/proc/loadavg").read_text().split()[:3]
    events = [json.loads(l) for l in (out / f"{entry['id']}.events.jsonl").read_text().splitlines() if l.strip()]
    inputs = [e for e in events if e.get("event") == "inputs"]
    fits = [fit_summary(e) for e in events if e.get("event") == "fit"]
    unavailable = [e["scheme"] for e in events if e.get("event") == "unavailable"]
    outer = [f for f in fits if f["scheme"] != "final"]
    want_outer = sum(s["folds"] for s in rel["splits"] if not s.get("unavailable_reason"))
    r7_fit = entry["lasso_seconds"] / entry["outer_fits"]
    fit_budget = max(60.0, 20.0 * r7_fit)
    dataset_budget = 5.0 * (entry["outer_fits"] + 1) * r7_fit
    rel_lasso_s = sum(l["seconds"] for s in rel["splits"] for l in s["learners"] if l["learner"] == "lasso")
    inp = inputs[0] if inputs else {}
    row = {
        "id": entry["id"], "species": entry["species"], "rank_rev7": entry["rank"],
        "started_utc": started_utc, "wall_seconds": wall, "exit": rc, "threads": args.threads,
        "load_before": load0, "load_after": load1,
        "input_fingerprint": inp.get("input_fingerprint"), "release_fingerprint": rel["input_fingerprint"],
        "fingerprint_match": inp.get("input_fingerprint") == rel["input_fingerprint"],
        "samples": inp.get("samples"), "features": inp.get("features"),
        "release_samples": rel["counts"]["samples"], "release_features": rel["counts"]["features"],
        "load_seconds": inp.get("seconds"),
        "unavailable_splits": unavailable,
        "fits": fits, "outer_fits": len(outer), "release_outer_fits": want_outer,
        "has_final": bool(fits) and fits[-1]["scheme"] == "final",
        "lasso_seconds_outer": sum(f["seconds"] for f in outer),
        "final_fit_seconds": fits[-1]["seconds"] if fits and fits[-1]["scheme"] == "final" else None,
        "slowest_fit_seconds": max((f["seconds"] for f in fits), default=None),
        "points": sum(f["path_points"] + f["refit_points"] for f in fits),
        "capped_points": sum(f["capped_points"] for f in fits),
        "irls_max": max((f["irls_max"] for f in fits), default=None),
        "rev7_lasso_seconds": entry["lasso_seconds"], "release_lasso_seconds": rel_lasso_s,
        "fit_budget_h1": fit_budget, "dataset_budget_h1": dataset_budget,
    }
    row["ratio_rev7_h1"] = row["lasso_seconds_outer"] / entry["lasso_seconds"]
    row["ratio_release"] = row["lasso_seconds_outer"] / rel_lasso_s if rel_lasso_s else None
    row["over_fit_budget_h1"] = [f"{f['scheme']} {f['outer_fold']}" for f in fits if f["seconds"] > fit_budget]
    row["over_dataset_budget_h1"] = sum(f["seconds"] for f in fits) > dataset_budget
    row["converged"] = (rc == 0 and row["fingerprint_match"] and len(outer) == want_outer and row["has_final"]
                        and all(f["converged"] for f in fits))
    reasons = []
    if rc != 0:
        reasons.append(f"exit {rc}")
    if not row["fingerprint_match"]:
        reasons.append("input fingerprint differs from the release")
    if len(outer) != want_outer or not row["has_final"]:
        reasons.append(f"{len(outer)} outer fits for {want_outer}, final {row['has_final']}")
    reasons += [f"{f['scheme']} {f['outer_fold']} unconverged or failed" for f in fits if not f["converged"]]
    row["registered_rule_would_stop"] = reasons + (["fit budget at h=1"] if row["over_fit_budget_h1"] else []) + (
        ["dataset budget at h=1"] if row["over_dataset_budget_h1"] else [])
    row["not_converged_reasons"] = reasons
    return row


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", type=Path, required=True)
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--release", type=Path, required=True)
    ap.add_argument("--list", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--only", default="")
    ap.add_argument("--dump", type=Path)
    ap.add_argument("--dump-full", default="")
    args = ap.parse_args(argv)
    args.binary = args.binary.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.dump:
        args.dump.mkdir(parents=True, exist_ok=True)
    datasets = ordered(json.loads(args.list.read_text())["datasets"])
    if args.only:
        keep = set(args.only.split(","))
        datasets = [d for d in datasets if d["id"] in keep]
    meta = {"started_utc": now(), "binary": str(args.binary),
            "binary_sha256": subprocess.run(["sha256sum", str(args.binary)], capture_output=True, text=True).stdout.split()[0],
            "nproc": os.cpu_count(), "threads": args.threads,
            "cpu": next(l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo") if l.startswith("model name")),
            "order": [d["id"] for d in datasets], "dump": str(args.dump) if args.dump else None}
    (args.out / "meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    rows = []
    for entry in datasets:
        row = run_dataset(args, entry, args.out)
        rows.append(row)
        print(f"{now()} {entry['id']}: converged={row['converged']} exit={row['exit']} "
              f"lasso {row['lasso_seconds_outer']:.1f} s, wall {row['wall_seconds']:.1f} s", flush=True)
        (args.out / "summary.json").write_text(json.dumps({"meta": meta, "datasets": rows}, indent=1) + "\n")
    meta["ended_utc"] = now()
    (args.out / "summary.json").write_text(json.dumps({"meta": meta, "datasets": rows}, indent=1) + "\n")
    (args.out / "scan.done").write_text(now() + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
