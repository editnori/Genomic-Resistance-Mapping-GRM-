#!/usr/bin/env python3
"""Interrupt the pilot's disk-sparse run at a checked point (stdlib, Linux).

usage: python3 -I pilot_interrupt.py --data ROOT --grm GRM --folded N
           --receipt OUT.json [--slug mycobacterium_tuberculosis]
           [--timeout 7200] [--poll 0.2]

The plan starts `GRM --data ROOT matrix ... --disk-sparse` under the
resource guard and /usr/bin/time, then runs this. It:

1. finds the one process whose command line starts with exactly
   `GRM --data ROOT matrix`, by reading /proc. The guard's and time's
   command lines start with python3 and /usr/bin/time, so neither is ever
   signalled (a `pkill -f` pattern matches all three);
2. polls the run's `run.json` until `folded_columns` reaches N, so the
   blocks before column N are folded: verified partials, their owned ksets
   deleted;
3. stops the process (SIGSTOP), records the state it stopped in, and
   requires it to be an interrupted run: `folded_columns` exactly N, no
   saved matrix, the folded blocks' partials and receipts present, and no
   finished matrix.json;
4. kills it (SIGKILL) and waits until it is gone.

It exits 0 only when all of that holds. When the run ends by itself first,
when the condition is not reached within --timeout, when the process is not
found, or when two processes match, it exits non-zero and says why, so the
plan never resumes from an unknown state. The receipt records every state
it saw.
"""

import argparse
import json
import os
import signal
import sys
import time
from pathlib import Path


def cmdline(pid):
    try:
        raw = Path(f"/proc/{pid}/cmdline").read_bytes()
    except OSError:
        return None
    return [a.decode(errors="replace") for a in raw.split(b"\0")[:-1]]


def find_run(grm, data, pids=None, read=cmdline):
    """PIDs whose command line starts with exactly `grm --data data matrix`."""
    want = [grm, "--data", data, "matrix"]
    found = []
    for entry in os.listdir("/proc") if pids is None else pids:
        if str(entry).isdigit() and int(entry) != os.getpid():
            argv = read(entry)
            if argv and argv[:4] == want:
                found.append(int(entry))
    return sorted(found)


def alive(pid):
    try:
        status = Path(f"/proc/{pid}/status").read_text()
    except OSError:
        return False
    state = next((line.split()[1] for line in status.splitlines() if line.startswith("State:")), "")
    return state != "Z"  # a zombie has exited; its parent has not reaped it yet


def run_state(species: Path):
    """(run dir, parsed run.json) of the single sparse run, or (None, None)."""
    partials = species / "matrix-partials"
    try:
        dirs = [d for d in partials.iterdir() if d.is_dir() and not d.is_symlink()]
    except OSError:
        return None, None
    if len(dirs) != 1:
        return None, None
    try:
        return dirs[0], json.loads((dirs[0] / "run.json").read_text())
    except (OSError, ValueError):
        return dirs[0], None  # absent or mid-rename: poll again


def snapshot(species: Path, run_dir, run):
    files = sorted(p.name for p in run_dir.iterdir()) if run_dir else []
    ksets = species / "ksets"
    return {
        "folded_columns": None if run is None else run.get("folded_columns"),
        "active": None if run is None else run.get("active"),
        "saved": None if run is None else run.get("saved") is not None,
        "block_width": None if run is None else run.get("config", {}).get("block_width"),
        "genomes": None if run is None else len(run.get("config", {}).get("genomes", [])),
        "run_dir": None if run_dir is None else str(run_dir),
        "run_files": files,
        "ksets": len([p for p in ksets.iterdir() if p.name.endswith(".kset")]) if ksets.is_dir() else 0,
        "matrix_json": (species / "matrix" / "matrix.json").exists(),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", required=True)
    ap.add_argument("--grm", required=True)
    ap.add_argument("--folded", required=True, type=int)
    ap.add_argument("--receipt", required=True, type=Path)
    ap.add_argument("--slug", default="mycobacterium_tuberculosis")
    ap.add_argument("--timeout", type=float, default=7200)
    ap.add_argument("--poll", type=float, default=0.2)
    args = ap.parse_args()
    species = Path(args.data) / "species" / args.slug
    receipt = {"data": args.data, "grm": args.grm, "folded_target": args.folded, "started": time.time(),
               "interrupted": False, "error": None, "seen": []}

    def finish(error=None):
        receipt["error"] = error
        receipt["finished"] = time.time()
        with args.receipt.open("x") as f:
            json.dump(receipt, f, indent=1)
            f.write("\n")
        if error:
            print(f"pilot_interrupt: STOP: {error}", file=sys.stderr)
            return 1
        print(json.dumps({k: receipt[k] for k in ("pid", "interrupted", "stopped_state")}, indent=1))
        return 0

    deadline = time.time() + args.timeout
    pids = []
    while time.time() < deadline and not pids:
        pids = find_run(args.grm, args.data)
        if not pids:
            time.sleep(args.poll)
    if len(pids) != 1:
        return finish(f"expected one `{args.grm} --data {args.data} matrix` process, found {pids}")
    pid = receipt["pid"] = pids[0]
    receipt["cmdline"] = cmdline(pid)
    last = None
    while True:
        if not alive(pid):
            receipt["final_state"] = snapshot(species, *run_state(species))
            return finish("the run ended before the kill point; nothing was interrupted, so the resume is not exercised")
        if time.time() > deadline:
            return finish(f"folded_columns did not reach {args.folded} within {args.timeout} s")
        run_dir, run = run_state(species)
        folded = None if run is None else run.get("folded_columns")
        if folded != last:
            receipt["seen"].append({"t": time.time(), "folded_columns": folded})
            last = folded
        if folded is not None and folded >= args.folded:
            break
        time.sleep(args.poll)
    os.kill(pid, signal.SIGSTOP)
    receipt["stopped_at"] = time.time()
    run_dir, run = run_state(species)
    state = receipt["stopped_state"] = snapshot(species, run_dir, run)
    width = state["block_width"] or 0
    want_files = {f"block-{c:05d}.partial" for c in range(0, args.folded, width)} if width else set()
    problems = []
    if state["folded_columns"] != args.folded:
        problems.append(f"folded_columns is {state['folded_columns']}, not {args.folded}")
    if state["saved"]:
        problems.append("the merged matrix is already saved")
    if state["genomes"] is not None and args.folded >= state["genomes"]:
        problems.append("every column is folded; no block is left to resume")
    missing = sorted(f for f in want_files | {f + ".receipt.json" for f in want_files} if f not in state["run_files"])
    if missing:
        problems.append(f"folded partials or receipts missing: {missing}")
    os.kill(pid, signal.SIGKILL)
    gone_by = time.time() + 60
    while alive(pid) and time.time() < gone_by:
        time.sleep(0.05)
    receipt["killed_at"] = time.time()
    receipt["gone"] = not alive(pid)
    if not receipt["gone"]:
        problems.append("the process did not exit after SIGKILL")
    receipt["after_kill_state"] = snapshot(species, *run_state(species))
    if problems:
        return finish("; ".join(problems))
    receipt["interrupted"] = True
    return finish()


if __name__ == "__main__":
    sys.exit(main())
