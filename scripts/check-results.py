"""Audit completed GRM reports on Vast without changing the running pipeline.

Revisions 3 to 6 are schema 1: every learner scores every isolate. Revision
7 is schema 2: a gbt or lasso score is null for each isolate of a fold
whose fit was null, and the row's no_call counts those isolates by class;
metrics cover the called isolates and are null when none is called; and
each split that ran saves test_fold, the outer fold that held out each
isolate.

The best learner follows the schema. In schema 2 it is the highest
random-split balanced accuracy among scm, cart, rf, gbt, and lasso rows
that call every isolate on every split that ran, or none. Schema-1 runs
logged the highest random-split balanced accuracy of any row, the catalog
included. A tie goes to the later row in both. With --log, each report's
last done line in that pipeline log must name that learner, and every
revision-7 report needs one. --skip-log supplies only class-minimum skip
evidence and does not require done lines; it defaults to --log.

usage: python3 -I check-results.py [DATA_ROOT] [--log PIPELINE_LOG] [--skip-log PIPELINE_LOG] [--reshaped]
--reshaped accepts only schema 2, for revision-5 and revision-6 reports
reshaped to schema 2 for checking.
"""
import argparse
import datetime
import json
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from result_scope import select_results

LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
NO_WINNER = "no learner calls every isolate"
DONE = re.compile(r"\s+(\S+__\S+) done in \S+ \((.*)\)$")
WINNER = re.compile(r"best random-split learner among those that call every isolate: (\S+) (\S+)")
ANY_ROW = re.compile(r"best random-split balanced accuracy: (\S+) (\S+)")


def done_lines(log):
    """The last done line for each dataset: (form, learner, balanced accuracy text)."""
    done = {}
    for line in log.read_text(errors="replace").splitlines():
        m = DONE.match(line)
        if not m:
            continue
        words = m.group(2)
        if words == NO_WINNER:
            done[m.group(1)] = ("full coverage", None, None)
        elif w := WINNER.fullmatch(words):
            done[m.group(1)] = ("full coverage", w.group(1), w.group(2))
        elif w := ANY_ROW.fullmatch(words):
            done[m.group(1)] = ("any row", w.group(1), w.group(2))
        else:
            done[m.group(1)] = ("unrecognized", words, None)
    return done


def balanced_accuracy(truth, scores):
    """Balanced accuracy at the 0.5 cutoff; a class with no isolate scores 0, as in grm."""
    resistant = [s >= 0.5 for t, s in zip(truth, scores) if t]
    susceptible = [s < 0.5 for t, s in zip(truth, scores) if not t]
    sensitivity = sum(resistant) / len(resistant) if resistant else 0.0
    specificity = sum(susceptible) / len(susceptible) if susceptible else 0.0
    return (sensitivity + specificity) / 2


def full_coverage_best(report):
    """(learner, random-split balanced accuracy) by the schema-2 rule, or None."""
    rows = {(split["scheme"], learner["learner"]): learner for split in report["splits"] for learner in split["learners"]}

    def calls_every_isolate(learner):
        no_call = learner.get("no_call")
        return learner["metrics"] is not None and (no_call is None or no_call["resistant"] + no_call["susceptible"] == 0)

    lineage_ran = any(split["scheme"] == "lineage" and split.get("unavailable_reason") is None for split in report["splits"])
    best = None
    for name in LEARNERS:
        random, lineage = rows.get(("random", name)), rows.get(("lineage", name))
        if random is None or not calls_every_isolate(random) or (lineage_ran and (lineage is None or not calls_every_isolate(lineage))):
            continue
        if best is None or random["metrics"]["balanced_accuracy"] >= best[1]:
            best = (name, random["metrics"]["balanced_accuracy"])
    return best


def any_row_best(report):
    """(learner, balanced accuracy) of the best random-split row, the catalog included: the schema-1 log rule."""
    best = None
    for learner in next(s for s in report["splits"] if s["scheme"] == "random")["learners"]:
        if best is None or learner["metrics"]["balanced_accuracy"] >= best[1]:
            best = (learner["learner"], learner["metrics"]["balanced_accuracy"])
    return best


def main():
    parser = argparse.ArgumentParser(description="Audit completed GRM dataset reports.")
    parser.add_argument("root", nargs="?", default=Path(__file__).resolve().parents[1] / "data", type=Path)
    parser.add_argument("--out", type=Path, help="receipt output path (default ROOT/logs/result-validation-receipt.json)")
    parser.add_argument("--log", type=Path, help="pipeline log whose done lines must name each report's best learner")
    parser.add_argument("--skip-log", type=Path, help="class-minimum skip evidence only; no done-line validation (defaults to --log)")
    parser.add_argument("--reshaped", action="store_true", help="accept only schema 2, for earlier revisions reshaped to it")
    args = parser.parse_args()
    root = args.root
    logged = done_lines(args.log) if args.log else None
    skip_log = args.skip_log or args.log
    files, scope = select_results(root, log=skip_log)
    assert not scope["missing_curated_results"], ("missing results for curated dataset", scope["missing_curated_results"])
    checked = []
    for path in files:
        report = json.loads(path.read_text())
        assert report["id"] == path.stem and report["species_slug"] == path.parent.parent.name, ("report identity differs from its curated path", path)
        audit = report["evaluation"]
        counts = report["counts"]
        samples = report["samples"]
        schema = report.get("schema")
        n = len(samples)
        truth = [sample["resistant"] for sample in samples]
        nr = sum(truth)
        ns = n - nr
        assert audit["revision"] in (3, 4, 5, 6, 7) and audit["valid"] is True, ("revision", path)
        assert schema == (2 if args.reshaped or audit["revision"] >= 7 else 1), ("schema does not match the revision", path, schema)
        assert "not external or clinical validation" in audit["scope"], path
        assert nr > 0 and ns > 0, path
        assert len({sample["id"] for sample in samples}) == n, path
        assert not any(sample["from_mic"] for sample in samples), path
        assert (counts["samples"], counts["resistant"], counts["susceptible"], counts["labels_from_mic"]) == (n, nr, ns, 0), path
        entries = 0
        no_call_rows = 0
        for split in report["splits"]:
            if split.get("unavailable_reason"):
                assert split["scheme"] == "lineage" and not split["learners"] and "test_fold" not in split, ("unavailable split", path)
                continue
            folds = split.get("test_fold")
            if schema == 1:
                assert folds is None, ("test_fold in a schema-1 report", path, split["scheme"])
            else:
                assert folds is not None and len(folds) == n and all(0 <= f < split["folds"] for f in folds), ("test_fold", path, split["scheme"])
                for k in range(split["folds"]):
                    held = [truth[i] for i in range(n) if folds[i] == k]
                    assert any(held) and not all(held), ("test fold lacks a class", path, split["scheme"], k)
                if split["scheme"] == "lineage":
                    lineage_folds = {}
                    for sample, f in zip(samples, folds):
                        lineage_folds.setdefault(sample["lineage"], set()).add(f)
                    assert all(len(v) == 1 for v in lineage_folds.values()), ("a lineage spans test folds", path)
            names = [learner["learner"] for learner in split["learners"]]
            assert names == LEARNERS + (["catalog"] if report.get("catalog") else []), ("learner rows", path, split["scheme"], names)
            for learner in split["learners"]:
                scores = learner["scores"]
                metrics = learner["metrics"]
                no_call = learner.get("no_call")
                identity = (str(path), split["scheme"], learner["learner"])
                assert len(scores) == n and all(v is None or (isinstance(v, (int, float)) and math.isfinite(v)) for v in scores), ("scores",) + identity
                uncalled = {"resistant": sum(1 for t, v in zip(truth, scores) if v is None and t), "susceptible": sum(1 for t, v in zip(truth, scores) if v is None and not t)}
                if learner["learner"] == "catalog" or schema == 1:
                    assert no_call is None and None not in scores, ("null score or no_call on a row that calls every isolate",) + identity
                else:
                    assert no_call == uncalled, ("no_call differs from the null scores",) + identity
                    assert learner["learner"] in ("gbt", "lasso") or no_call == {"resistant": 0, "susceptible": 0}, ("only gbt and lasso fits can be null",) + identity
                called = [i for i, v in enumerate(scores) if v is not None]
                no_call_rows += len(called) < n
                if learner["learner"] != "catalog":
                    assert len(learner["fold_balanced_accuracy"]) == split["folds"], ("fold balanced accuracy length",) + identity
                    assert learner["conformal"] == [], identity
                if not called:
                    assert metrics is None and learner["roc"] == [] and all(b is None for b in learner["fold_balanced_accuracy"]), ("a row without calls needs null metrics, no ROC, and null fold balanced accuracy",) + identity
                    entries += 1
                    continue
                assert metrics is not None, ("null metrics on a row with calls",) + identity
                predicted = {i: scores[i] >= 0.5 for i in called}
                confusion = {
                    "tp": sum(truth[i] and predicted[i] for i in called),
                    "tn": sum(not truth[i] and not predicted[i] for i in called),
                    "fp": sum(not truth[i] and predicted[i] for i in called),
                    "fn": sum(truth[i] and not predicted[i] for i in called),
                }
                assert metrics["n"] == len(called) and all(metrics[k] == v for k, v in confusion.items()), ("confusion of the called isolates",) + identity
                assert math.isclose(metrics["balanced_accuracy"], balanced_accuracy([truth[i] for i in called], [scores[i] for i in called]), abs_tol=1e-12), ("balanced accuracy",) + identity
                if learner["learner"] != "catalog" and folds is not None:
                    fold_ba = learner["fold_balanced_accuracy"]
                    for k in range(split["folds"]):
                        held = [i for i in range(n) if folds[i] == k]
                        fold = [scores[i] for i in held]
                        if all(v is None for v in fold):
                            assert fold_ba[k] is None, ("a fold without calls needs null balanced accuracy",) + identity + (k,)
                            continue
                        assert None not in fold and fold_ba[k] is not None, ("a fold is called in part",) + identity + (k,)
                        assert math.isclose(fold_ba[k], balanced_accuracy([truth[i] for i in held], fold), abs_tol=1e-12), ("fold balanced accuracy",) + identity + (k,)
                entries += 1
        mic = report.get("mic")
        if mic:
            assert len(mic["predicted"]) == n, path
            assert all(v is None or math.isfinite(v) for v in mic["predicted"]), path
            assert sum(v is not None for v in mic["predicted"]) == mic["n"], path
            assert 0 <= mic["n_exact"] <= mic["n"] <= n, path
            if mic["categorical"]:
                assert mic["categorical"]["n"] == mic["n"], path
        best = full_coverage_best(report) if schema == 2 else any_row_best(report)
        if logged is not None:
            line = logged.get(report["id"])
            assert line is not None or audit["revision"] < 7, ("revision-7 report without a done line", path)
            assert line is None or line[0] == ("full coverage" if schema == 2 else "any row"), ("done line form does not match the report schema", path, line)
            if line is not None:
                want = (None, None) if best is None else (best[0], "%.3f" % best[1])
                assert line[1:] == want, ("done line names another best learner", path, line, want)
        checked.append({"id": report["id"], "revision": audit["revision"], "schema": schema, "n": n, "resistant": nr, "metric_entries": entries, "rows_with_no_call": no_call_rows,
                        "best_learner": None if best is None else best[0], "log_checked": logged is not None and report["id"] in logged,
                        "unavailable_splits": {s["scheme"]: s["unavailable_reason"] for s in report["splits"] if s.get("unavailable_reason")}})

    receipt = {"checked_at": datetime.datetime.now(datetime.timezone.utc).isoformat(), "log": None if args.log is None else str(args.log), "skip_log": None if skip_log is None else str(skip_log), "reports": checked, **scope}
    out = args.out or root / "logs/result-validation-receipt.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"reports": len(checked), "revisions": sorted({r["revision"] for r in checked}), "schemas": sorted({r["schema"] for r in checked}),
                      "rows_with_no_call": sum(r["rows_with_no_call"] for r in checked), "without_best_learner": sum(r["best_learner"] is None for r in checked),
                      "log_checked": sum(r["log_checked"] for r in checked), **scope}))


if __name__ == "__main__":
    main()
