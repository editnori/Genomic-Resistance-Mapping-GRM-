"""Write the hyperparameter selection report for one or more trial histories.

Each run is one history (one outer training fold or repeat) with its protocol.
select.py, beside this file, checks each history and selects a config per
family by the one-standard-error rule; this writes a JSON report and a short
Markdown summary. Both list every trial, excluded ones included, with its
status and reason, and give each family's selection, candidate set and
budget, the budget flags, and, across runs, the stability of the selected
hyperparameters. Both say the SE is descriptive (select.SE_LABEL).

usage: python3 -I -B scripts/hyperopt/report.py
           --run LABEL PROTOCOL.json HISTORY.jsonl [--run ...]
           --out-json REPORT.json --out-md REPORT.md

PROTOCOL.json is Task B's protocol (protocol.py). Its protocol_hash must
match its content, so a protocol edited after the search is refused. The
flat test-only protocol shape is refused here. See select.Protocol.from_dict.

Exit status: 0 written; 2 refused, with every problem on stderr and nothing
written.
"""
import argparse
import importlib.util
import json
import os
import sys
from pathlib import Path

SCHEMA = "grm-hyperopt-report/1"


def _load_select():
    if __package__:
        from . import select
        return select
    name = "grm_hyperopt_select"
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(
            name, Path(__file__).resolve().with_name("select.py"))
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


select = _load_select()
HistoryError = select.HistoryError


def build_report(runs, trusted_flat=False):
    """runs: a list of {"label", "protocol", "trials"}, one per history.
    Returns the report dict; raises HistoryError if any history or protocol
    does not match, or the runs cannot be pooled for stability (see
    select.stability). trusted_flat accepts flat protocols, for tests only."""
    if not runs:
        raise HistoryError(["no_runs: the report needs at least one run"])
    out, selections, problems = [], [], []
    for run in runs:
        try:
            protocol = select.as_protocol(run["protocol"], trusted_flat)
            s = select.select_history(run["trials"], protocol)
        except HistoryError as e:
            problems += [f"{p} [run {run['label']}]" for p in e.problems]
            continue
        selections.append(s)
        trials = [dict(row, family=fam)
                  for fam, f in s["families"].items() for row in f["trials"]]
        out.append({
            "label": run["label"],
            "protocol": protocol.to_dict(),
            "history_sha256": s["history_sha256"],
            "folds": s["folds"],
            "n_trials": len(trials),
            "families": {fam: {k: v for k, v in f.items()
                               if k not in ("trials", "se_label")}
                         for fam, f in s["families"].items()},
            "budgets": s["budgets"],
            "trials": trials,
        })
    if problems:
        raise HistoryError(problems)
    stab = select.stability(selections, [r["label"] for r in runs])
    return {
        "schema": SCHEMA,
        "objective": selections[0]["objective"],
        "direction": selections[0]["direction"],
        "se_label": select.SE_LABEL,
        "rule": select.RULE,
        "runs": out,
        "stability": stab,
    }


def render_markdown(report):
    """The plain Markdown summary of a report."""
    lines = [
        "# Hyperparameter selection report",
        "",
        f"Objective: {report['objective']} ({report['direction']}).",
        f"Rule: {report['rule']}.",
        "",
        f"SE note: {report['se_label']}",
    ]
    for run in report["runs"]:
        lines += _run_lines(run)
    lines += _stability_lines(report["stability"])
    return "\n".join(lines) + "\n"


def _run_lines(run):
    p = run["protocol"]
    statuses = {}
    for t in run["trials"]:
        statuses[t["status"]] = statuses.get(t["status"], 0) + 1
    lines = [
        "", f"## Run {run['label']}", "",
        f"protocol_hash {p['protocol_hash']}, fold_plan_hash "
        f"{p['fold_plan_hash']}, history sha256 {run['history_sha256']}, "
        f"inner folds {_folds(run['folds'])}; "
        f"{run['n_trials']} trials ("
        + ", ".join(f"{n} {s}" for s, n in sorted(statuses.items())) + ").",
        "", "Constraints: " + ("; ".join(
            f"{c['metric']} {c['op']} {c['bound']}" for c in p["constraints"])
            or "none") + ".",
        "", "### Selection", "",
        "| family | selected trial | complexity | default | mean | best trial "
        "| best mean | best SE (descriptive) | candidate set |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for fam, f in run["families"].items():
        if f["selected"] is None:
            lines.append(f"| {fam} | none ({f['status']}) |  |  |  |  |  |  |  |")
            continue
        s, b = f["selected"], f["best"]
        lines.append(
            f"| {fam} | {s['trial_id']} | {_num(s['complexity'])} | "
            f"{_yes(s['is_default'])} | {_num(s['mean'])} | {b['trial_id']} | "
            f"{_num(b['mean'])} | {_num(b['se'])} | "
            + ", ".join(str(c["trial_id"]) for c in f["candidate_set"]) + " |")
    for fam, f in run["families"].items():
        if f["selected"] is not None:
            lines.append(f"\n{fam} selected config: "
                         f"`{json.dumps(f['selected']['config'], sort_keys=True)}`")

    lines += ["", "### Budgets", "",
              "| family | planned | used | ok | failed | unconverged "
              "| feasible | infeasible | fold set mismatch | fits | seconds "
              "| cost unreported |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for fam, b in run["budgets"]["families"].items():
        lines.append(f"| {fam} | {b['planned']} | {b['used']} | {b['ok']} | "
                     f"{b['failed']} | {b['unconverged']} | {b['feasible']} | "
                     f"{b['infeasible']} | {b['fold_set_mismatch']} | "
                     f"{b['fit_count']} | "
                     f"{_num(b['seconds'])} | {b['cost_unreported']} |")
    flags = run["budgets"]["flags"]
    lines += ["", "Budget flags:" + (" none." if not flags else "")]
    lines += [f"- {f['code']}: {f['message']}" for f in flags]

    excluded = [t for t in run["trials"] if t["excluded"]]
    lines += ["", f"### Excluded trials ({len(excluded)})", ""]
    if excluded:
        lines += ["| family | trial | status | reason | excluded because |",
                  "|---|---|---|---|---|"]
        lines += [f"| {t['family']} | {t['trial_id']} | {t['status']} | "
                  f"{_cell(t['reason'], t['status'] != 'ok')} | "
                  f"{_cell('; '.join(t['excluded_because']))} |"
                  for t in excluded]
    else:
        lines.append("None.")

    lines += ["", f"### All trials ({run['n_trials']})", "",
              "| family | trial | status | reason | excluded because "
              "| complexity | default | mean | SE (descriptive) | candidate "
              "| selected |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
    for t in run["trials"]:
        lines.append(
            f"| {t['family']} | {t['trial_id']} | {t['status']} | "
            f"{_cell(t['reason'], t['status'] != 'ok')} | "
            f"{_cell('; '.join(t['excluded_because']), False)} | "
            f"{_num(t['complexity'])} | "
            f"{_yes(t['is_default'])} | {_num(t['mean'])} | {_num(t['se'])} | "
            f"{_yes(t['in_candidate_set'])} | {_yes(t['selected'])} |")

    unused = [t for t in run["trials"]
              if t["fold_scores"] and not t["scores_used"]]
    lines += ["", f"### Fold scores not used ({len(unused)})", "",
              "Reported by failed, unconverged and fold-set-mismatch trials; "
              "kept as a record, used for nothing.", ""]
    if unused:
        lines += ["| family | trial | status | fold: score |", "|---|---|---|---|"]
        lines += [f"| {t['family']} | {t['trial_id']} | {t['status']} | "
                  + ", ".join(f"{s['fold']}: {_num(s['score']) or '-'}"
                              for s in t["fold_scores"]) + " |"
                  for t in unused]
    else:
        lines.append("None.")
    return lines


def _stability_lines(stab):
    lines = ["", f"## Stability across runs ({len(stab['runs'])})"]
    for fam, f in stab["families"].items():
        lines += ["", f"### {fam}", "",
                  f"Selected in {f['selected_in']} of {f['n_runs']} runs; the "
                  f"default was selected {f['default_selected']} times."]
        if f["runs_without_selection"]:
            lines.append("No selection in: "
                         + ", ".join(f["runs_without_selection"]) + ".")
        if not f["selected_in"]:
            continue
        lines += ["", "| selected config | count | runs |", "|---|---|---|"]
        lines += [f"| `{json.dumps(c['config'], sort_keys=True)}` | "
                  f"{c['count']} | {', '.join(c['runs'])} |"
                  for c in f["configs"]]
        lines += ["", "Min, median and max cover the numeric values; "
                  "skipped counts the others (such as an internal-search "
                  "default's value lists).", "",
                  "| hyperparameter | present | distinct | modal fraction "
                  "| min | median | max | skipped |",
                  "|---|---|---|---|---|---|---|---|"]
        rows = dict(f["params"], **{"(complexity)": {
            "present": f["selected_in"], "distinct": None,
            "modal_fraction": None, "numeric": f["complexity"]}})
        for name, prm in rows.items():
            n = prm["numeric"] or {}
            lines.append(
                f"| {name} | {prm['present']} | {_num(prm['distinct'])} | "
                f"{_num(prm['modal_fraction'])} | {_num(n.get('min'))} | "
                f"{_num(n.get('median'))} | {_num(n.get('max'))} | "
                f"{_num(n.get('skipped'))} |")
    return lines


def _folds(folds):
    return ", ".join(str(f) for f in folds) if folds else "none"


def _num(x):
    if x is None:
        return ""
    return str(x) if isinstance(x, int) else f"{x:.6g}"


def _yes(b):
    return "yes" if b else ""


def _cell(text, required=True):
    """Text as one table cell; a missing required text says so."""
    if text is None:
        return "(none given)" if required else ""
    return " ".join(str(text).split()).replace("|", "\\|")


def _write_all(outputs):
    """Write each (path, bytes) to a temporary file beside its path, and
    rename them into place only once every one is written. Any temporary
    file left by a failure is removed."""
    tmps = []
    try:
        for path, data in outputs:
            path = Path(path)
            tmps.append((path.with_name(f".{path.name}.{os.getpid()}.tmp"), path))
            tmps[-1][0].write_bytes(data)
        for tmp, path in tmps:
            os.replace(tmp, path)
    finally:
        for tmp, _ in tmps:
            tmp.unlink(missing_ok=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", nargs=3, action="append", required=True,
                    metavar=("LABEL", "PROTOCOL", "HISTORY"))
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-md", required=True)
    args = ap.parse_args(argv)
    try:
        runs = []
        for label, protocol, history in args.run:
            runs.append({"label": label,
                         "protocol": select.load_protocol(protocol),
                         "trials": select.load_history(history)})
        report = build_report(runs)
        # Encode both before writing either: a run label or family name that
        # is not valid UTF-8 (a lone surrogate) fails here, not mid-write.
        json_bytes = (json.dumps(report, indent=1, allow_nan=False)
                      + "\n").encode("utf-8")
        md_bytes = render_markdown(report).encode("utf-8")
        _write_all([(args.out_json, json_bytes), (args.out_md, md_bytes)])
    except HistoryError as e:
        print("refused:", file=sys.stderr)
        for p in e.problems:
            print(f"  {p}", file=sys.stderr)
        return 2
    except UnicodeEncodeError as e:
        print(f"refused: bad_name: a run label or family name is not valid "
              f"UTF-8 ({e})", file=sys.stderr)
        return 2
    except (OSError, ValueError, ArithmeticError, RecursionError) as e:
        print(f"refused: {type(e).__name__}: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
