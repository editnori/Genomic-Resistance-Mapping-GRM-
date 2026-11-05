"""Write the hyperparameter selection report for one or more trial histories.

Each run is one history (one outer training fold or repeat) with its protocol.
select.py, beside this file, checks each history and selects a config per
family (select.RULE: the fixed default, DECISIONS F4's fallback);
this writes a JSON report and a short
Markdown summary. Both list every trial, excluded ones included, with its
status and reason, and give each family's selection, candidate set and
budget, the budget flags, and, across runs, the stability of the selected
hyperparameters. Both say the SE is descriptive (select.SE_LABEL).

Every mean and SE in the report is an inner-tuning score, never test
performance (select.SCORE_LABEL), and both outputs say so. This package
does not refit a selected config on the outer training fold or score it
on the untouched outer test fold; both outputs name that as a remaining
dependency (select.OUTER_EVALUATION).

The JSON also says so in its key names. The `mean`, `sd` and `se` that
select.py gives each trial, and the `mean` and `se` of each family's
`best`, `selected` and `candidate_set` entries, are
`inner_tuning_mean`, `inner_tuning_sd` and `inner_tuning_se` here
(INNER_TUNING), so a value read on its own keeps its label. Schema 2 made
this rename and added each run's `warnings` (below); a schema 1 report used
the bare keys. Schema 3 selects each family's fixed default, not the
one-SE rule's pick: a family's `selected` and only candidate is its
default, `best` and `default_gain` are kept for the record, a family has
`default_feasible` and no `threshold`, a history without a default trial
gives `no_default`, and the selected entry's inner mean is labelled a
selection-biased tuning score (select.SELECTED_SCORE_LABEL).
Stability's `mean` is the mean of selected hyperparameter values, not a
score, and keeps its name.

usage: python3 -I -B scripts/hyperopt/report.py
           --run LABEL PROTOCOL.json HISTORY.jsonl FOLD_PLAN.json [--run ...]
           --out-json REPORT.json --out-md REPORT.md

PROTOCOL.json is Task B's protocol (protocol.py). Its protocol_hash must
match its content, so a protocol edited after the search is refused. The
flat test-only protocol shape is refused here. See select.Protocol.from_dict.
FOLD_PLAN.json is the fold plan the protocol pins. The history must be the
one the protocol's search writes on it: every planned trial with its
planned config, in run order, each scored on the plan's inner folds (see
select.validate_history).

Each run gives two digests of its history. `history_sha256` is the sha256
of its records as canonical JSON lines (select.history_sha256), the bytes
Task B's search writes; `history_file_sha256` is the sha256 of the file's
bytes, read once and parsed from those same bytes, so the records and the
digest describe one file. They differ exactly when the file is not in that
canonical form.

A run whose protocol searches a parameter production has no knob for, such
as a lasso protocol from before review F1, is selected as usual and
labelled in both outputs (`warnings`; select.NO_PRODUCTION_KNOB).

Exit status: 0 written; 2 refused, with every problem on stderr and nothing
written.
"""
import argparse
import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path

SCHEMA = "grm-hyperopt-report/3"
# The report's names for select.py's inner-tuning scores (review F3).
INNER_TUNING = {"mean": "inner_tuning_mean", "sd": "inner_tuning_sd",
                "se": "inner_tuning_se"}


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
    """runs: a list of {"label", "protocol", "fold_plan"}, one per history,
    each with its history as exactly one of
    - "history_bytes", the bytes of its file. This parses them
      (select.parse_history) and hashes the same bytes, so
      history_file_sha256 is the digest of the bytes the records came from;
    - "trials", records built in memory; history_file_sha256 is then null.
    The digest is never taken from the caller: a run that supplies
    "history_file_sha256" is refused (untrusted_digest).
    "fold_plan" is required with Task B's protocol. Returns the report dict;
    raises HistoryError if any history, protocol or fold plan does not
    match, or the runs cannot be pooled for stability (see
    select.stability). trusted_flat accepts flat protocols, for tests only."""
    if not runs:
        raise HistoryError(["no_runs: the report needs at least one run"])
    out, selections, problems = [], [], []
    for run in runs:
        try:
            records, file_sha256 = _history(run)
            protocol = select.as_protocol(run["protocol"], trusted_flat)
            s = select.select_history(records, protocol, run.get("fold_plan"))
        except HistoryError as e:
            problems += [f"{p} [run {run['label']}]" for p in e.problems]
            continue
        selections.append(s)
        trials = [_inner_tuning(dict(row, family=fam))
                  for fam, f in s["families"].items() for row in f["trials"]]
        out.append({
            "label": run["label"],
            "protocol": protocol.to_dict(),
            "history_sha256": s["history_sha256"],
            "history_file_sha256": file_sha256,
            "folds": s["folds"],
            "warnings": s["warnings"],
            "n_trials": len(trials),
            "families": {fam: _family(f) for fam, f in s["families"].items()},
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
        "score_label": select.SCORE_LABEL,
        "outer_evaluation": select.OUTER_EVALUATION,
        "rule": select.RULE,
        "runs": out,
        "stability": stab,
    }


def _history(run):
    """A run's records, and the sha256 of the bytes they were parsed from
    (None for records given as such)."""
    if "history_file_sha256" in run:
        raise HistoryError(["untrusted_digest: history_file_sha256 is the "
                            "sha256 of history_bytes, computed here, never "
                            "taken from the caller"])
    if ("trials" in run) == ("history_bytes" in run):
        raise HistoryError(["bad_run: give the history as exactly one of "
                            "trials and history_bytes"])
    if "trials" in run:
        return run["trials"], None
    data = run["history_bytes"]
    return select.parse_history(data), hashlib.sha256(data).hexdigest()


def _inner_tuning(row):
    """row with its mean, SD and SE under their INNER_TUNING names."""
    return {INNER_TUNING.get(k, k): v for k, v in row.items()}


def _family(f):
    """A family's selection without its trial rows, which the run lists,
    and with every score named for the inner tuning it comes from."""
    out = {k: v for k, v in f.items() if k not in ("trials", "se_label")}
    for key in ("best", "selected"):
        if out[key] is not None:
            out[key] = _inner_tuning(out[key])
    out["candidate_set"] = [_inner_tuning(c) for c in out["candidate_set"]]
    return out


def render_markdown(report):
    """The plain Markdown summary of a report."""
    lines = [
        "# Hyperparameter selection report",
        "",
        f"Objective: {report['objective']} ({report['direction']}).",
        f"Rule: {report['rule']}.",
        "",
        f"Scores: {report['score_label']}",
        "",
        f"Outer evaluation: {report['outer_evaluation']}",
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
        f"{p['fold_plan_hash']}, canonical history sha256 "
        f"{run['history_sha256']}, history file sha256 "
        f"{run['history_file_sha256'] or 'not given'}, "
        f"inner folds {_folds(run['folds'])}; "
        f"{run['n_trials']} trials ("
        + ", ".join(f"{n} {s}" for s, n in sorted(statuses.items())) + ").",
    ]
    if run["warnings"]:
        lines += ["", "Warnings:"]
        lines += [f"- {w['code']}: {w['message']}" for w in run["warnings"]]
    lines += [
        "", "Constraints: " + ("; ".join(
            f"{c['metric']} {c['op']} {c['bound']}" for c in p["constraints"])
            or "none") + ".",
        "", "### Selection", "",
        "| family | selected trial | complexity | default "
        "| inner-tuning mean (selection-biased) | best trial "
        "| best inner-tuning mean "
        "| best inner-tuning SE (descriptive) | candidate set |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for fam, f in run["families"].items():
        if f["selected"] is None:
            lines.append(f"| {fam} | none ({f['status']}) |  |  |  |  |  |  |  |")
            continue
        s, b = f["selected"], f["best"] or {}
        lines.append(
            f"| {fam} | {s['trial_id']} | {_num(s['complexity'])} | "
            f"{_yes(s['is_default'])} | {_num(s['inner_tuning_mean'])} | "
            f"{b.get('trial_id', '')} | {_num(b.get('inner_tuning_mean'))} | "
            f"{_num(b.get('inner_tuning_se'))} | "
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
              "| complexity | default | inner-tuning mean "
              "| inner-tuning SE (descriptive) | candidate | selected |",
              "|---|---|---|---|---|---|---|---|---|---|---|"]
    for t in run["trials"]:
        lines.append(
            f"| {t['family']} | {t['trial_id']} | {t['status']} | "
            f"{_cell(t['reason'], t['status'] != 'ok')} | "
            f"{_cell('; '.join(t['excluded_because']), False)} | "
            f"{_num(t['complexity'])} | "
            f"{_yes(t['is_default'])} | {_num(t['inner_tuning_mean'])} | "
            f"{_num(t['inner_tuning_se'])} | "
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
    ap.add_argument("--run", nargs=4, action="append", required=True,
                    metavar=("LABEL", "PROTOCOL", "HISTORY", "FOLD_PLAN"))
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-md", required=True)
    args = ap.parse_args(argv)
    try:
        runs = []
        for label, protocol, history, fold_plan in args.run:
            runs.append({"label": label,
                         "protocol": select.load_protocol(protocol),
                         # One read; build_report parses and hashes it.
                         "history_bytes": Path(history).read_bytes(),
                         "fold_plan": select.load_fold_plan(fold_plan)})
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
