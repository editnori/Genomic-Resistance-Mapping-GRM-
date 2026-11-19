#!/usr/bin/env python3
"""Interrupt the pilot's disk-sparse run at a checked point (stdlib, Linux).

usage: python3 -I pilot_interrupt.py --launcher PID --guard-resources GUARD.json
           --data ROOT --grm GRM --folded N --receipt OUT.json
           [--slug mycobacterium_tuberculosis] [--timeout 7200] [--poll 0.2]
       python3 -I pilot_interrupt.py --list-runs --data ROOT --grm GRM

The plan starts the resource guard in the background (guard -> /usr/bin/time
-> `GRM --data ROOT matrix ... --disk-sparse`), then runs this with the
guard's PID (`$!`) as --launcher and the guard's --resources path. It:

1. reserves --receipt (an existing one is never overwritten, and nothing is
   signalled before it is reserved), and requires the launcher to be a
   running sibling of this process (a child of the same shell) whose
   command line is the resource guard of this run: `python3
   .../resource_guard.py`, `--resources` equal to --guard-resources, and
   `GRM --data ROOT matrix` after its `--` (N1 of the re-check). Any other
   launcher is refused before anything is signalled. It also requires the
   guard's resources record not to exist yet;
2. binds to the one process below the launcher whose command line starts
   with exactly `GRM --data ROOT matrix`, read from /proc (L1). Only the
   launcher's descendants can be bound, so a run left from an earlier
   attempt is never bound or signalled; and since two runs would share the
   root, any matching process outside the launcher's tree is a stop. The
   guard's and time's command lines start with python3 and /usr/bin/time,
   so neither is signalled on the way to a good interrupt;
3. polls the run's `run.json` until `folded_columns` reaches N, so the
   blocks before column N are folded: verified partials, their owned ksets
   deleted;
4. stops the bound process (SIGSTOP), scans again and requires the same
   single match, by PID and start time; records the state it stopped in and
   requires an interrupted run: `folded_columns` exactly N, no saved matrix,
   no matrix.json (I5), and the folded blocks' partials and receipts;
5. kills it (SIGKILL), waits until it is gone and the launcher has exited,
   and requires the guard's record to show the run killed by SIGKILL (rc 137
   from /usr/bin/time, or -9 without it) with no guard reason, for this
   command, and the root to be still in the stopped state.

It exits 0 only when all of that holds. On any stop (the run ends first, the
point is not reached within --timeout, no match, a second or foreign match,
a failed check, an error or a signal) it stops the launched run before it
exits 1: it freezes the launcher and everything below it, kills what is
below it, deepest first, and lets the launcher record that and exit (the
launcher is killed too if nothing ran below it or it does not exit). So the
plan never leaves a run going, and never resumes from an unknown state. The
receipt records every state it saw and what it signalled.

After any stop it lists every `GRM --data ROOT matrix` process still running,
in the receipt and on stderr, and never signals one (N2): the guard starts
its command in a new session, so a guard that dies before the bind leaves
GRM running outside the launch. --list-runs prints that list alone and exits
1 when it is not empty; the plan runs it after every step-a `wait`.
"""

import argparse
import json
import os
import signal
import sys
import time
from pathlib import Path

KILLED = (137, -9)  # /usr/bin/time exits 128 + 9 for a SIGKILLed command; Popen reports -9
GUARD = "resource_guard.py"


class Stopped(Exception):
    """A signal to this helper: stop the launched run, then exit."""


def cmdline(pid):
    try:
        raw = Path(f"/proc/{pid}/cmdline").read_bytes()
    except OSError:
        return None
    return [a.decode(errors="replace") for a in raw.split(b"\0")[:-1]]


def proc_stat(pid):
    """(state, ppid, start time in clock ticks) of a process, or None."""
    try:
        raw = Path(f"/proc/{pid}/stat").read_text()
    except OSError:
        return None
    fields = raw[raw.rindex(")") + 2:].split()  # the command name may hold spaces or parentheses
    return fields[0], int(fields[1]), int(fields[19])


def sibling(pid):
    """The launcher's identity when it is a live sibling of this helper (a child of the same shell,
    as the plan's `$!` is), else None. A wrong --launcher (1, the shell, a process in a subshell) is
    refused and never signalled; launcher_identity() also refuses a sibling that is not the guard."""
    if pid <= 1 or pid in (os.getpid(), os.getppid()):
        return None
    st = proc_stat(pid)
    if st is None or st[0] in ("Z", "X") or st[1] != os.getppid():
        return None
    return (pid, st[2])


def guard_problem(argv, resources, grm, data):
    """Why a command line is not the plan's resource guard for this run, else None: `python3
    .../resource_guard.py`, `--resources` exactly once before its `--` and equal to --guard-resources,
    and `GRM --data ROOT matrix` after it."""
    if not argv or len(argv) < 2 or not Path(argv[0]).name.startswith("python") or Path(argv[1]).name != GUARD:
        return f"its command line {argv} is not `python3 .../{GUARD} ...`"
    if "--" not in argv:
        return "its command line has no `--` before the guarded command"
    cut = argv.index("--")
    values = [argv[i + 1] for i in range(2, cut - 1) if argv[i] == "--resources"]
    if values != [resources]:
        return f"its --resources is {values}, not [{resources!r}]"
    if not contains(argv[cut + 1:], [grm, "--data", data, "matrix"]):
        return f"it does not guard `{grm} --data {data} matrix`"
    return None


def launcher_identity(args, wait=5.0):
    """(pid, start, argv) of the plan's own `$!` and None, else (None, why); nothing is signalled
    either way. The launcher must be a live sibling (sibling()) whose command line is this run's
    guard (guard_problem()). bash forks before it execs the guard, so the command line gets `wait`
    seconds to become the guard's."""
    ident = sibling(args.launcher)
    if ident is None:
        return None, (f"the launcher {args.launcher} is not a running sibling job of this helper (a child of the same "
                      "shell); pass the plan's own `$!`. Nothing was signalled")
    deadline = time.time() + wait
    while True:
        argv = cmdline(ident[0])  # read before alive(): a live identity means argv was this process's
        why = guard_problem(argv, args.guard_resources, args.grm, args.data)
        if not alive(ident):
            return None, f"the launcher {args.launcher} exited before it was identified. Nothing was signalled"
        if why is None:
            return ident + (argv,), None
        if time.time() > deadline:
            return None, (f"the launcher {args.launcher} is not the resource guard of this run: {why}; pass the plan's "
                          "own `$!`. Nothing was signalled")
        time.sleep(0.05)


def alive(ident):
    st = proc_stat(ident[0])
    return st is not None and st[2] == ident[1] and st[0] not in ("Z", "X")  # a zombie has exited


def descendants(launcher):
    """{pid: (depth, start)} of every live process below the launcher."""
    children = {}
    for entry in os.listdir("/proc"):
        if entry.isdigit():
            st = proc_stat(entry)
            if st is not None and st[0] not in ("Z", "X"):
                children.setdefault(st[1], []).append((int(entry), st[2]))
    out, todo = {}, [(launcher[0], 0)]
    while todo:
        parent, depth = todo.pop()
        for pid, start in children.get(parent, []):
            if pid not in out:
                out[pid] = (depth + 1, start)
                todo.append((pid, depth + 1))
    return out


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


def split_matches(matches, tree):
    """The matches below the launcher, as identities, and those outside its tree (L1)."""
    own = [(pid, tree[pid][1]) for pid in matches if pid in tree]
    return own, [pid for pid in matches if pid not in tree]


def signal_ident(ident, sig):
    """Signal a process only while it is still the one identified."""
    if not alive(ident):
        return False
    try:
        os.kill(ident[0], sig)
    except ProcessLookupError:
        return False
    return True


def wait_gone(idents, seconds):
    deadline = time.time() + seconds
    while any(alive(i) for i in idents) and time.time() < deadline:
        time.sleep(0.05)
    return not any(alive(i) for i in idents)


def stop_launched(launcher, grace=60, settle=5):
    """Stop the launched run: freeze the launcher and everything below it (repeating until no new
    process appears, so a fork in flight cannot escape), kill what is below it deepest first, then
    let the launcher record that and exit; kill the launcher too if nothing ran below it (it may
    not have started the run yet) or if it does not exit. A launcher with nothing below it first
    gets `settle` seconds to exit by itself: the guard of a run that has just ended is still
    writing its record. Returns what was signalled."""
    record = {"frozen": [], "killed": [], "launcher_killed": False}
    if launcher is None or not alive(launcher):
        record["launcher_exited"] = launcher is not None
        return record
    st = proc_stat(launcher[0])
    if st is None or st[2] != launcher[1]:
        record["launcher_exited"] = True
        return record
    if st[1] != os.getppid() or cmdline(launcher[0]) != launcher[2]:  # checked again before any signal
        record["refused"] = "the launcher is no longer this run's guard and a sibling of this helper"
        return record
    deadline = time.time() + settle
    while alive(launcher) and not descendants(launcher) and time.time() < deadline:
        time.sleep(0.05)
    if not alive(launcher):
        record["launcher_exited"] = True
        return record
    signal_ident(launcher, signal.SIGSTOP)
    record["frozen"].append(launcher[0])
    frozen = {}
    for _ in range(100):
        new = {pid: v for pid, v in descendants(launcher).items() if pid not in frozen}
        if not new:
            break
        for pid, (depth, start) in new.items():
            signal_ident((pid, start), signal.SIGSTOP)
            frozen[pid] = (depth, start)
            record["frozen"].append(pid)
        time.sleep(0.05)
    below = sorted(frozen.items(), key=lambda kv: -kv[1][0])
    for pid, (_, start) in below:
        if signal_ident((pid, start), signal.SIGKILL):
            record["killed"].append(pid)
    record["below_gone"] = wait_gone([(pid, start) for pid, (_, start) in below], grace)
    if below:
        signal_ident(launcher, signal.SIGCONT)  # the guard reaps its child, writes its record, exits
        if not wait_gone([launcher], grace):
            record["launcher_killed"] = signal_ident(launcher, signal.SIGKILL)
    else:
        record["launcher_killed"] = signal_ident(launcher, signal.SIGKILL)
    record["launcher_exited"] = wait_gone([launcher], grace)
    return record


def still_running(grm, data):
    """Every `GRM --data ROOT matrix` process still running, listed and never signalled (N2)."""
    return [{"pid": pid, "cmdline": cmdline(pid)} for pid in find_run(grm, data)]


def left_message(left, grm, data):
    return (f"; `{grm} --data {data} matrix` still running, not signalled: {[r['pid'] for r in left]}; "
            "stop them before anything else" if left else "")


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


def interrupted_problems(state, folded):
    """Why a snapshot is not the interrupted run the resume must start from; [] when it is."""
    width = state["block_width"] or 0
    want = {f"block-{c:05d}.partial" for c in range(0, folded, width)} if width else set()
    problems = []
    if state["folded_columns"] != folded:
        problems.append(f"folded_columns is {state['folded_columns']}, not {folded}")
    if state["saved"] or state["saved"] is None:
        problems.append("the merged matrix is already saved" if state["saved"] else "no run.json was readable")
    if state["matrix_json"]:
        problems.append("a matrix.json exists")
    if state["genomes"] is not None and folded >= state["genomes"]:
        problems.append("every column is folded; no block is left to resume")
    missing = sorted(f for f in want | {f + ".receipt.json" for f in want} if f not in state["run_files"])
    if not want or missing:
        problems.append(f"folded partials or receipts missing: {missing or 'no block width'}")
    return problems


def contains(seq, part):
    return any(seq[i:i + len(part)] == part for i in range(len(seq) - len(part) + 1))


def interrupt(args, species, launcher, receipt):
    """Steps 1-5 for an identified launcher; returns None on a checked interrupt, else the reason to stop."""
    receipt["launcher"] = {"pid": launcher[0], "start": launcher[1], "cmdline": launcher[2]}
    if os.path.lexists(args.guard_resources):
        return f"{args.guard_resources} exists before this run ended: a record from another launch"
    deadline = time.time() + args.timeout
    while True:
        if not alive(launcher):
            return "the launcher exited before the run was found"
        own, foreign = split_matches(find_run(args.grm, args.data), descendants(launcher))
        if foreign:
            return (f"`{args.grm} --data {args.data} matrix` process(es) {foreign} not launched by this run are running "
                    "against the same root; they were not signalled; stop them and start again")
        if len(own) > 1:
            return f"expected one `{args.grm} --data {args.data} matrix` process below the launcher, found {[p for p, _ in own]}"
        if own:
            break
        if time.time() > deadline:
            return f"no `{args.grm} --data {args.data} matrix` process appeared below the launcher within {args.timeout} s"
        time.sleep(args.poll)
    run = own[0]
    receipt["pid"], receipt["pid_start"], receipt["cmdline"] = run[0], run[1], cmdline(run[0])
    last = None
    while True:
        if not alive(run):
            receipt["final_state"] = snapshot(species, *run_state(species))
            return "the run ended before the kill point; nothing was interrupted, so the resume is not exercised"
        if time.time() > deadline:
            return f"folded_columns did not reach {args.folded} within {args.timeout} s"
        run_dir, state = run_state(species)
        folded = None if state is None else state.get("folded_columns")
        if folded != last:
            receipt["seen"].append({"t": time.time(), "folded_columns": folded})
            last = folded
        if folded is not None and folded >= args.folded:
            break
        time.sleep(args.poll)
    if not signal_ident(run, signal.SIGSTOP):
        return "the run ended at the kill point, before it could be stopped"
    receipt["stopped_at"] = time.time()
    problems = []
    own, foreign = split_matches(find_run(args.grm, args.data), descendants(launcher))
    receipt["matches_at_stop"] = {"below_launcher": [p for p, _ in own], "outside": foreign}
    if own != [run] or foreign:
        problems.append(f"at the stop the matching processes are {[p for p, _ in own]} below the launcher and {foreign} "
                        f"outside it, not only {run[0]}")
    state = receipt["stopped_state"] = snapshot(species, *run_state(species))
    problems += interrupted_problems(state, args.folded)
    signal_ident(run, signal.SIGKILL)
    receipt["killed_at"] = time.time()
    receipt["gone"] = wait_gone([run], 60)
    if not receipt["gone"]:
        problems.append("the process did not exit after SIGKILL")
    receipt["launcher_exited"] = wait_gone([launcher], 120)
    if not receipt["launcher_exited"]:
        problems.append("the launcher did not exit after the run was killed")
    else:
        try:
            guard = json.loads(Path(args.guard_resources).read_text())
            receipt["guard_record"] = {k: guard.get(k) for k in ("rc", "guard_reason", "command")}
            if guard.get("rc") not in KILLED or guard.get("guard_reason") is not None:
                problems.append(f"the guard recorded rc {guard.get('rc')} and reason {guard.get('guard_reason')!r}, "
                                f"not a SIGKILLed run (rc in {KILLED}, no reason)")
            if not contains(guard.get("command") or [], [args.grm, "--data", args.data, "matrix"]):
                problems.append("the guard's record is for another command")
        except (OSError, ValueError, AttributeError) as exc:
            problems.append(f"the guard's record {args.guard_resources} is unreadable: {exc}")
    after = receipt["after_kill_state"] = snapshot(species, *run_state(species))
    if after["run_dir"] != state["run_dir"] or interrupted_problems(after, args.folded):
        problems.append("after the kill the root is not in the stopped state: "
                        + "; ".join(interrupted_problems(after, args.folded) or ["another run directory"]))
    if problems:
        return "; ".join(problems)
    receipt["interrupted"] = True
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--list-runs", action="store_true",
                    help="only list every `GRM --data ROOT matrix` process still running; never signal; exit 1 if any")
    ap.add_argument("--launcher", type=int)
    ap.add_argument("--guard-resources")
    ap.add_argument("--data", required=True)
    ap.add_argument("--grm", required=True)
    ap.add_argument("--folded", type=int)
    ap.add_argument("--receipt", type=Path)
    ap.add_argument("--slug", default="mycobacterium_tuberculosis")
    ap.add_argument("--timeout", type=float, default=7200)
    ap.add_argument("--poll", type=float, default=0.2)
    args = ap.parse_args(argv)
    if args.list_runs:
        left = still_running(args.grm, args.data)
        print(json.dumps(left, indent=1))
        if left:
            print(f"pilot_interrupt: STOP{left_message(left, args.grm, args.data)}", file=sys.stderr)
        return 1 if left else 0
    missing = [name for name in ("launcher", "guard_resources", "folded", "receipt") if getattr(args, name) is None]
    if missing:
        ap.error("the following arguments are required: " + ", ".join("--" + n.replace("_", "-") for n in missing))
    species = Path(args.data) / "species" / args.slug

    def stop(signum, frame):
        raise Stopped(f"signal {signum}")

    def ignore_signals():  # finish stopping the run before anything else (N5)
        for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
            signal.signal(sig, signal.SIG_IGN)
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(sig, stop)
    try:
        out = args.receipt.open("x")
    except OSError as exc:
        ignore_signals()
        launcher, refusal = launcher_identity(args)
        stopped = stop_launched(launcher) if refusal is None else {"refused": refusal}
        left = still_running(args.grm, args.data)
        print(f"pilot_interrupt: STOP: cannot reserve the receipt {args.receipt}: {exc}; launched run stopped: {stopped}"
              + left_message(left, args.grm, args.data), file=sys.stderr)
        return 1
    receipt = {"data": args.data, "grm": args.grm, "folded_target": args.folded, "started": time.time(),
               "interrupted": False, "error": None, "seen": []}
    with out:
        launcher = None
        try:
            launcher, error = launcher_identity(args)
            error = error or interrupt(args, species, launcher, receipt)
        except BaseException as exc:  # any error or signal still stops the launched run
            error = f"{type(exc).__name__}: {exc}"
        if error:
            receipt["interrupted"] = False
            ignore_signals()
            receipt["launched_run_stopped"] = stop_launched(launcher)
            receipt["still_running"] = still_running(args.grm, args.data)
            error += left_message(receipt["still_running"], args.grm, args.data)
        receipt["error"] = error
        receipt["finished"] = time.time()
        json.dump(receipt, out, indent=1)
        out.write("\n")
    if error:
        print(f"pilot_interrupt: STOP: {error}", file=sys.stderr)
        return 1
    print(json.dumps({k: receipt[k] for k in ("pid", "interrupted", "stopped_state", "guard_record")}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
