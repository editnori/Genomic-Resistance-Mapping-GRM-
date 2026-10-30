"""Run the large-species lasso scan under its registered stop rule.

The rule is registered in docs/lasso-stopping.md ("Registered scan stop
rule"); the constants below are that rule, and test_scan.py pins them. The
scan list, scan-list.json beside this file, names the 18 datasets and their
revision-7 lasso times.

For each dataset, most expensive first, this runs the `lasso-scan` example
binary (main.rs beside this file), which fits only the lasso of the shipped
evaluator on the dataset's revision-7 inputs. It times every fit and stops at
once, killing the binary, on any of:

- a lasso error, an unconverged path point, or an invalid evaluation audit;
- a fit, dataset, or whole-scan wall-clock budget overrun;
- free disk that would fall below the floor before the next dataset;
- inputs that are not the revision-7 report's (its input fingerprint).

It never retries. It writes only OUT: plan.json, events.jsonl, one stderr log
per dataset, summary.json, and stop.json when stopped. The stderr log holds the
binary's progress line for each penalty point as it is solved, so a stopped
fit shows how far it got.

usage: python3 -I -B scan.py --input-root ROOT --out OUT [--binary BIN]
           [--list LIST] [--host-factor H] [--commit SHA] [--plan | --inputs-only]

--plan checks the list against ROOT's revision-7 reports, prints the budgets
and the disk estimate, and exits without running anything. --inputs-only also
runs BIN to load each dataset and compare its input fingerprint with the
revision-7 report's, without fitting. --host-factor h scales every budget: the
scan host's time for a fixed piece of work divided by the revision-7 host's,
so 1 on the Vast host class that ran revision 7.

Exit status: 0 finished, go for the 275-report run; 1 finished, no-go; 2
stopped by the rule (OUT/stop.json says why); 3 refused to start.
"""
import argparse
import ctypes
import datetime
import hashlib
import json
import math
import os
import selectors
import shutil
import signal
import statistics
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Rule:
    """The registered stop rule. `r7_fit` is a dataset's revision-7 lasso
    seconds per outer fit, and h the host factor."""

    fit_floor_seconds: float = 60.0
    # Per fit (tuning path and refit): max(floor, fit_factor * r7_fit * h).
    fit_factor: float = 20.0
    # Per dataset (outer fits and the final fit):
    # dataset_factor * (outer fits + 1) * r7_fit * h. The whole scan's
    # wall-clock budget is the sum over the list.
    dataset_factor: float = 5.0
    # Review before the 275-report run when more of the solved penalty
    # points (tuning and refit) than this hit the sweep cap in some IRLS step.
    cap_review_share: float = 0.05
    # Go only if the worst and median ratio of scan to revision-7 lasso
    # seconds and the projected 275-report run stay within these.
    max_ratio: float = 5.0
    max_median_ratio: float = 3.0
    max_projected_hours: float = 16.0
    # The projection's allowance for the untimed final fit: 11 fits per 10.
    final_fit_share: float = 1.1
    # Free disk that the scan must leave on OUT's filesystem.
    disk_floor_bytes: int = 30 * 2**30
    # Revision 7 ran with 64 rayon threads; the scan matches it.
    rayon_threads: int = 64


REGISTERED = Rule()

# Disk the scan writes, bounded from above. A fit line holds at most 30
# tuning and 30 refit statuses (about 4 KiB) and a refit with at most one
# weight per screened feature: 4,000 weights of at most 41 bytes each, about
# 160 KiB. The stderr log holds one progress line of under 140 bytes per
# penalty point: 11 fits of at most 60 points, under 100 KiB. The binary and
# the driver write no temporary files (the driver runs with -B).
FIT_LINE_BYTES = 256 * 2**10
START_LINE_BYTES = 2**10
INPUTS_LINE_BYTES = 4 * 2**10
STDERR_BYTES = 2**20
FIXED_BYTES = 4 * 2**20


class Stop(Exception):
    """The stop rule fired: `kind` names the condition."""

    def __init__(self, kind, detail, **facts):
        super().__init__(f"{kind}: {detail}")
        self.kind, self.detail, self.facts = kind, detail, facts


class Refuse(Exception):
    """A precondition failed before anything ran."""


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def species_of(dataset_id):
    return dataset_id.split("__", 1)[0]


def revision7_timing(root, entry):
    """The revision-7 report's lasso seconds and outer fits for `entry`,
    checked against the list."""
    path = Path(root, "species", entry["species"], "results", entry["id"] + ".json")
    if not path.is_file():
        raise Refuse(f"{entry['id']}: no revision-7 report at {path}")
    report = json.loads(path.read_text())
    seconds, fits = 0.0, 0
    for split in report["splits"]:
        if split.get("unavailable_reason"):
            continue
        rows = [r for r in split["learners"] if r["learner"] == "lasso"]
        if len(rows) != 1:
            raise Refuse(f"{entry['id']}: the {split['scheme']} split of {path} has {len(rows)} lasso rows")
        seconds += rows[0]["seconds"]
        fits += split["folds"]
    facts = {
        "revision": report["evaluation"]["revision"],
        "input_fingerprint": report.get("input_fingerprint"),
        "lasso_seconds": seconds,
        "outer_fits": fits,
    }
    want = {"revision": 7, "input_fingerprint": entry["input_fingerprint"], "outer_fits": entry["outer_fits"]}
    if any(facts[k] != v for k, v in want.items()) or not math.isclose(seconds, entry["lasso_seconds"], rel_tol=1e-9, abs_tol=1e-9):
        raise Refuse(f"{entry['id']}: the report at {path} is not the listed revision-7 report: {facts}")
    return facts


def budgets(entry, rule, host):
    r7_fit = entry["lasso_seconds"] / entry["outer_fits"]
    return {
        "r7_fit_seconds": r7_fit,
        "fit_seconds": max(rule.fit_floor_seconds, rule.fit_factor * r7_fit * host),
        "dataset_seconds": rule.dataset_factor * (entry["outer_fits"] + 1) * r7_fit * host,
    }


def disk_need(entry):
    """Bytes the scan writes for one dataset, bounded from above."""
    fits = entry["outer_fits"] + 1
    return INPUTS_LINE_BYTES + fits * (FIT_LINE_BYTES + START_LINE_BYTES) + STDERR_BYTES


def inputs_read(root, datasets):
    """Bytes of the inputs the binary reads; it writes none of them."""
    paths = {Path(root, "curated/datasets.json"), Path(root, "curated/phenotypes.jsonl"), Path(root, "bvbrc/genomes.jsonl")}
    for slug in {e["species"] for e in datasets}:
        species = Path(root, "species", slug)
        paths |= {species / "counted.json", species / "structure.json"}
        paths |= set((species / "matrix").glob("*"))
    paths |= set(Path(root, "catalog").rglob("*"))
    return sum(p.stat().st_size for p in paths if p.is_file())


def free_bytes(path):
    """Free bytes on the filesystem that holds, or will hold, `path`."""
    path = Path(path).resolve()
    while not path.exists():
        path = path.parent
    return shutil.disk_usage(path).free


def die_with_parent():
    """Run in the child: if the driver dies, even by SIGKILL, so does the
    binary, so no fit outlives its budget."""
    PR_SET_PDEATHSIG = 1
    ctypes.CDLL(None, use_errno=True).prctl(PR_SET_PDEATHSIG, signal.SIGKILL)


def plan(args, rule, scan_list):
    host = args.host_factor
    if not host > 0:
        raise Refuse("--host-factor must be positive")
    datasets = scan_list["datasets"]
    if any(e["lasso_seconds"] <= 0 or e["outer_fits"] <= 0 for e in datasets):
        raise Refuse("every listed dataset needs positive revision-7 lasso seconds and outer fits")
    rows = []
    for entry in datasets:
        if species_of(entry["id"]) != entry["species"]:
            raise Refuse(f"{entry['id']} is not a {entry['species']} dataset")
        revision7_timing(args.input_root, entry)
        rows.append({"id": entry["id"], "species": entry["species"], "rank": entry.get("rank"), **budgets(entry, rule, host), "disk_bytes": disk_need(entry)})
    need = FIXED_BYTES + sum(r["disk_bytes"] for r in rows)
    return {
        "input_root": str(Path(args.input_root).resolve()),
        "out": str(Path(args.out).resolve()),
        "host_factor": host,
        "rule": asdict(rule),
        "list": {"path": str(Path(args.list).resolve()), "sha256": sha256(args.list)},
        "revision7": scan_list["revision7"],
        "datasets": rows,
        "total_seconds": sum(r["dataset_seconds"] for r in rows),
        "disk": {
            "write_bytes": need,
            "temp_bytes": 0,
            "inputs_read_bytes": inputs_read(args.input_root, datasets),
            "floor_bytes": rule.disk_floor_bytes,
        },
    }


def check_out(out, root):
    out, root = Path(out).resolve(), Path(root).resolve()
    if out == root or root in out.parents:
        raise Refuse(f"OUT {out} is inside the input root {root}")
    if out.exists() and any(out.iterdir()):
        raise Refuse(f"OUT {out} is not empty; the scan writes only to a new directory")


def disk_report(p, free):
    gib = 2**30
    d = p["disk"]
    return (f"disk: writes at most {d['write_bytes'] / 2**20:.1f} MiB, no temporary files; reads "
            f"{d['inputs_read_bytes'] / gib:.2f} GiB of inputs in place; free {free / gib:.2f} GiB, "
            f"floor {d['floor_bytes'] / gib:.0f} GiB, headroom after the scan {(free - d['write_bytes'] - d['floor_bytes']) / gib:.2f} GiB")


def check_disk(path, need, rule, when):
    free = free_bytes(path)
    if free - need < rule.disk_floor_bytes:
        raise Stop("disk floor", f"{when}: {free} bytes free minus {need} still to write is under the {rule.disk_floor_bytes}-byte floor",
                   free_bytes=free, need_bytes=need)
    return free


class Scan:
    def __init__(self, args, rule, scan_list, p, out):
        self.args, self.rule, self.list, self.plan, self.out = args, rule, scan_list, p, out
        self.events = open(out / "events.jsonl", "a", buffering=1)
        self.started = time.monotonic()
        self.deadline_total = self.started + p["total_seconds"]
        self.rows = []
        self.proc = None

    def log(self, record):
        record = dict(record, t=round(time.monotonic() - self.started, 3))
        self.events.write(json.dumps(record) + "\n")

    def kill(self):
        proc, self.proc = self.proc, None
        if proc is None:
            return
        for sig, wait in ((signal.SIGTERM, 5), (signal.SIGKILL, 30)):
            if proc.poll() is not None:
                break
            try:
                os.killpg(proc.pid, sig)
                proc.wait(timeout=wait)
            except (ProcessLookupError, subprocess.TimeoutExpired):
                pass
        proc.stdout.close()

    def run_binary(self, entry, inputs_only, handle, deadline):
        """Runs the binary on one dataset and passes each event to `handle`.
        `deadline()` returns (when, kind, detail) for the nearest budget."""
        cmd = [self.args.binary] + (["--inputs-only"] if inputs_only else []) + [str(self.args.input_root), entry["species"], entry["id"]]
        env = dict(os.environ, RAYON_NUM_THREADS=str(self.rule.rayon_threads))
        with open(self.out / f"{entry['id']}.stderr.log", "wb") as err:
            self.proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=err, stdin=subprocess.DEVNULL,
                                         start_new_session=True, bufsize=0, env=env, preexec_fn=die_with_parent)
            fd = self.proc.stdout.fileno()
            sel = selectors.DefaultSelector()
            sel.register(fd, selectors.EVENT_READ)
            buf = b""
            try:
                while True:
                    when, kind, detail = deadline()
                    now = time.monotonic()
                    if now >= when:
                        raise Stop(kind, detail, dataset=entry["id"])
                    if not sel.select(timeout=when - now):
                        continue
                    chunk = os.read(fd, 1 << 16)
                    if not chunk:
                        break
                    buf += chunk
                    while b"\n" in buf:
                        line, buf = buf.split(b"\n", 1)
                        try:
                            event = json.loads(line)
                        except json.JSONDecodeError:
                            raise Stop("unexpected output", repr(line[:200]), dataset=entry["id"])
                        if event.get("dataset") != entry["id"]:
                            raise Stop("unexpected output", f"event for {event.get('dataset')}", dataset=entry["id"])
                        self.log(event)
                        handle(event)
                if buf:
                    raise Stop("unexpected output", f"unterminated line {buf[:200]!r}", dataset=entry["id"])
                rc = self.proc.wait()
            finally:
                sel.close()
        self.proc.stdout.close()
        self.proc = None
        return rc

    def check_inputs_event(self, entry, event):
        if event.get("event") != "inputs":
            raise Stop("unexpected output", f"expected the inputs line first, got {event.get('event')}", dataset=entry["id"])
        if event["input_fingerprint"] != entry["input_fingerprint"]:
            raise Stop("inputs differ from revision 7", f"input fingerprint {event['input_fingerprint']}, revision 7 {entry['input_fingerprint']}", dataset=entry["id"])

    def dataset(self, entry, b):
        rule, h = self.rule, self.args.host_factor
        state = {"inputs": None, "running": None, "spent": 0.0, "fits": [], "unavailable": []}

        def deadline():
            options = [(self.deadline_total, "scan budget overrun", f"the scan passed its {self.plan['total_seconds']:.0f} s budget")]
            if state["running"] is not None:
                start, site = state["running"]
                options.append((start + b["fit_seconds"], "fit budget overrun", f"{site} passed its {b['fit_seconds']:.0f} s budget"))
                options.append((start + b["dataset_seconds"] - state["spent"], "dataset budget overrun",
                                f"{entry['id']} passed its {b['dataset_seconds']:.0f} s budget during {site}"))
            return min(options)

        def handle(event):
            kind = event.get("event")
            if state["inputs"] is None:
                self.check_inputs_event(entry, event)
                state["inputs"] = event
                return
            if kind == "unavailable":
                state["unavailable"].append(event["scheme"])
            elif kind == "start" and state["running"] is None:
                site = f"{entry['id']} {event['scheme']} fold {event['outer_fold']}"
                state["running"] = (time.monotonic(), site)
            elif kind == "fit" and state["running"] is not None:
                start, site = state["running"]
                elapsed = time.monotonic() - start
                state["running"] = None
                fit = summarize_fit(event)
                state["fits"].append(fit)
                if event.get("error"):
                    raise Stop("lasso error", event["error"], dataset=entry["id"], fit=site)
                audit = event["audit"]
                if audit is None or unconverged(audit):
                    raise Stop("unconverged point", f"{site}: {audit}", dataset=entry["id"], fit=site)
                if max(event["seconds"], elapsed) > b["fit_seconds"]:
                    raise Stop("fit budget overrun", f"{site} took {elapsed:.1f} s of {b['fit_seconds']:.0f}", dataset=entry["id"])
                state["spent"] += event["seconds"]
                if state["spent"] > b["dataset_seconds"]:
                    raise Stop("dataset budget overrun", f"{entry['id']} spent {state['spent']:.1f} s of {b['dataset_seconds']:.0f}", dataset=entry["id"])
            else:
                raise Stop("unexpected output", f"{kind} event out of order", dataset=entry["id"])

        rc = self.run_binary(entry, False, handle, deadline)
        if rc != 0:
            raise Stop("scan binary failed", f"exit status {rc}; see {entry['id']}.stderr.log", dataset=entry["id"])
        fits = state["fits"]
        outer = [f for f in fits if f["scheme"] != "final"]
        if state["inputs"] is None or len(outer) != entry["outer_fits"] or [f["scheme"] for f in fits[-1:]] != ["final"]:
            raise Stop("unexpected output", f"{len(fits)} fits for {entry['outer_fits']} outer fits and the final fit", dataset=entry["id"])
        # The report's audit is valid exactly when every fit has its record
        # and none has an unconverged point; each was checked as it arrived.
        valid = all(not unconverged(f["audit"]) for f in fits)
        if not valid:
            raise Stop("invalid evaluation audit", entry["id"], dataset=entry["id"])
        scan_seconds = sum(f["seconds"] for f in outer)
        return {
            "id": entry["id"],
            "rank": entry.get("rank"),
            "samples": state["inputs"]["samples"],
            "features": state["inputs"]["features"],
            "unavailable_splits": state["unavailable"],
            "r7_lasso_seconds": entry["lasso_seconds"],
            "lasso_seconds": scan_seconds,
            "final_fit_seconds": fits[-1]["seconds"],
            "ratio": scan_seconds / (entry["lasso_seconds"] * h),
            "slowest_fit_seconds": max(f["seconds"] for f in fits),
            "fit_budget_seconds": b["fit_seconds"],
            "dataset_budget_seconds": b["dataset_seconds"],
            "points": sum(f["points"] for f in fits),
            "capped_points": sum(f["capped_points"] for f in fits),
            "valid": valid,
            "fits": fits,
        }

    def inputs_only(self, entry):
        seen = []

        def handle(event):
            if seen:
                raise Stop("unexpected output", "more than one line in --inputs-only", dataset=entry["id"])
            self.check_inputs_event(entry, event)
            seen.append(event)

        limit = time.monotonic() + 900
        rc = self.run_binary(entry, True, handle, lambda: (min(limit, self.deadline_total), "inputs check overrun", "15 minutes per dataset"))
        if rc != 0 or not seen:
            raise Stop("scan binary failed", f"exit status {rc}; see {entry['id']}.stderr.log", dataset=entry["id"])
        return {k: seen[0][k] for k in ("dataset", "samples", "resistant", "features", "lineages", "min_support", "input_fingerprint", "seconds")}


def unconverged(audit):
    """`LassoFitAudit::has_unconverged`: the evaluation audit's test."""
    return audit["path_unconverged"] > 0 or any(
        audit[k] is not None and not audit[k]["converged"] for k in ("selected_convergence", "refit_convergence"))


def summarize_fit(event):
    trace = event["trace"]
    points = trace["path"] + trace["refit"]
    return {
        "scheme": event["scheme"],
        "outer_fold": event["outer_fold"],
        "seconds": event["seconds"],
        "audit": event["audit"],
        "points": len(points),
        "capped_points": sum(p["sweep_cap_hits"] > 0 for p in points),
        "max_irls_steps": max((p["irls_steps"] for p in points), default=0),
    }


def decide(rows, rule, revision7):
    ratios = [r["ratio"] for r in rows]
    worst, median = max(ratios), statistics.median(ratios)
    projected = (revision7["report_seconds"] + (median - 1) * rule.final_fit_share * revision7["lasso_seconds"]) / 3600
    points = sum(r["points"] for r in rows)
    capped = sum(r["capped_points"] for r in rows)
    share = capped / points if points else 0.0
    reasons = []
    if worst > rule.max_ratio:
        reasons.append(f"worst ratio {worst:.2f} is above {rule.max_ratio}")
    if median > rule.max_median_ratio:
        reasons.append(f"median ratio {median:.2f} is above {rule.max_median_ratio}")
    if projected > rule.max_projected_hours:
        reasons.append(f"the 275-report run projects to {projected:.1f} h, above {rule.max_projected_hours} h")
    if share > rule.cap_review_share:
        reasons.append(f"{capped} of {points} solved points ({share:.1%}) hit the sweep cap, above {rule.cap_review_share:.0%}: review first")
    return {
        "decision": "no-go" if reasons else "go",
        "reasons": reasons,
        "worst_ratio": worst,
        "median_ratio": median,
        "projected_hours": projected,
        "points": points,
        "capped_points": capped,
        "capped_share": share,
    }


def now_utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def write_json(path, value):
    path.write_text(json.dumps(value, indent=1) + "\n")


def prepare(args, rule):
    """The checked plan, or Refuse. Nothing has run and nothing is written."""
    try:
        scan_list = json.loads(Path(args.list).read_text())
        check_out(args.out, args.input_root)
        p = plan(args, rule, scan_list)
    except Refuse:
        raise
    except Exception as e:  # a malformed list or report must not read as a no-go
        raise Refuse(f"{type(e).__name__}: {e}") from e
    p["commit"] = args.commit
    return scan_list, p


def run(args, rule=REGISTERED):
    """The scan, or one of its checks. Returns the exit status."""
    try:
        scan_list, p = prepare(args, rule)
        free = free_bytes(args.out)
        print(disk_report(p, free), flush=True)
        if args.plan:
            for r in p["datasets"]:
                print(f"  {r['id']:52} fit {r['fit_seconds']:7.0f} s  dataset {r['dataset_seconds']:7.0f} s", flush=True)
            print(f"whole scan {p['total_seconds']:.0f} s ({p['total_seconds'] / 3600:.2f} h)", flush=True)
            return 0
        if free - p["disk"]["write_bytes"] < rule.disk_floor_bytes:
            raise Refuse(f"{free} bytes free minus the {p['disk']['write_bytes']}-byte estimate is under the {rule.disk_floor_bytes}-byte floor")
        if not (Path(args.binary).is_file() and os.access(args.binary, os.X_OK)):
            raise Refuse(f"{args.binary} is not an executable")
        p["binary"] = {"path": str(Path(args.binary).resolve()), "sha256": sha256(args.binary)}
    except Refuse as e:
        print(f"refused: {e}", file=sys.stderr, flush=True)
        return 3

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    p.update(mode="inputs-only" if args.inputs_only else "scan", launched_utc=now_utc(), free_bytes_at_launch=free)
    write_json(out / "plan.json", p)
    scan = Scan(args, rule, scan_list, p, out)
    by_id = {e["id"]: e for e in scan_list["datasets"]}
    remaining = p["disk"]["write_bytes"]
    summary = {"mode": p["mode"], "launched_utc": p["launched_utc"], "host_factor": args.host_factor}
    try:
        try:
            for r in p["datasets"]:
                entry = by_id[r["id"]]
                check_disk(out, remaining, rule, f"before {entry['id']}")
                if args.inputs_only:
                    scan.rows.append(scan.inputs_only(entry))
                    print(f"{entry['id']}: inputs match revision 7", flush=True)
                else:
                    row = scan.dataset(entry, r)
                    scan.rows.append(row)
                    print(f"{entry['id']}: {row['lasso_seconds']:.1f} s of lasso, ratio {row['ratio']:.2f}, "
                          f"slowest fit {row['slowest_fit_seconds']:.1f} s, {row['capped_points']}/{row['points']} points capped", flush=True)
                remaining -= r["disk_bytes"]
        except Stop:
            raise
        except BaseException as e:  # a driver fault or a signal still stops the scan on the record
            raise Stop("driver stopped", f"{type(e).__name__}: {e}") from e
    except Stop as stop:
        scan.kill()
        record = {"kind": stop.kind, "detail": stop.detail, **stop.facts, "stopped_utc": now_utc(),
                  "elapsed_seconds": round(time.monotonic() - scan.started, 3), "completed": [row.get("id", row.get("dataset")) for row in scan.rows]}
        write_json(out / "stop.json", record)
        summary.update(status="stopped", stop=record, datasets=scan.rows)
        write_json(out / "summary.json", summary)
        print(f"STOPPED: {stop}", file=sys.stderr, flush=True)
        return 2
    finally:
        scan.kill()
        scan.events.close()

    summary.update(status="complete", elapsed_seconds=round(time.monotonic() - scan.started, 3), datasets=scan.rows)
    if args.inputs_only:
        write_json(out / "summary.json", summary)
        print("all inputs match revision 7", flush=True)
        return 0
    summary.update(decide(scan.rows, rule, scan_list["revision7"]))
    write_json(out / "summary.json", summary)
    print(f"{summary['decision'].upper()}: worst ratio {summary['worst_ratio']:.2f}, median {summary['median_ratio']:.2f}, "
          f"projected {summary['projected_hours']:.1f} h, {summary['capped_share']:.1%} of points capped", flush=True)
    for reason in summary["reasons"]:
        print(f"  {reason}", flush=True)
    return 0 if summary["decision"] == "go" else 1


def parse(argv):
    parser = argparse.ArgumentParser(description="Run the large-species lasso scan under its registered stop rule.")
    parser.add_argument("--input-root", required=True, type=Path, help="data root holding the revision-7 inputs and reports (read only)")
    parser.add_argument("--out", required=True, type=Path, help="new or empty output directory outside the input root")
    parser.add_argument("--binary", type=Path, help="the lasso-scan example binary")
    parser.add_argument("--list", type=Path, default=HERE / "scan-list.json")
    parser.add_argument("--host-factor", type=float, default=1.0)
    parser.add_argument("--commit", help="the source commit the binary was built from, recorded in plan.json")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--plan", action="store_true", help="check the list and print budgets and disk; run nothing")
    mode.add_argument("--inputs-only", action="store_true", help="load each dataset and check its fingerprint; fit nothing")
    args = parser.parse_args(argv)
    if not args.plan and args.binary is None:
        parser.error("--binary is required unless --plan")
    # Path("./lasso-scan") prints as "lasso-scan", which subprocess would look
    # up on PATH instead of in the working directory.
    if args.binary is not None:
        args.binary = args.binary.resolve()
    return args


def on_sigterm(signum, frame):
    raise SystemExit(f"signal {signum}")


if __name__ == "__main__":
    signal.signal(signal.SIGTERM, on_sigterm)
    sys.exit(run(parse(sys.argv[1:])))
