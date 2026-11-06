"""Build the ablation report (JSON plus a Markdown summary) from a run manifest,
the contrast, pooled and feature records, and the protocol's planned contrasts.

usage: python3 -I -B scripts/ablation/report.py MANIFEST CONTRASTS POOLED FEATURES PLANNED PLANNED_POOLED
           OUT_DIR [--not-produced FILE] [--predictions-root DIR]
Each input is a JSON file (UTF-8). CONTRASTS, POOLED and FEATURES hold lists
of records, PLANNED the list of every planned per-dataset contrast id,
PLANNED_POOLED the list of every planned pooled key (families A, B and D),
and the optional --not-produced file maps a planned id or pooled key with no
record to its reason. With --predictions-root, every ok run's predictions
file (a path relative to DIR) is read and its rows validated against the
manifest. Renders OUT_DIR/ablation-report.json and OUT_DIR/ablation-report.md,
writes both as UTF-8 to temporary names (mode 0o666 & ~umask), moves any earlier
report in OUT_DIR to backup names, renames both into place and deletes the
backups. On any problem it prints "<path>: <message>" lines and exits 1, with
OUT_DIR as it was: no new file, an earlier report restored byte for byte, no
temporary or backup file left, and any folder this run created for OUT_DIR
removed again. If restoring a backup fails too, that backup is kept, never
deleted, and the message ends "earlier report kept at <backup path>".
"""
import argparse
import contextlib
import copy
import hashlib
import json
import os
import re
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

if __package__:
    from . import schema
else:  # run as a script; -I leaves the script's directory off sys.path
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import schema  # noqa: E402

JSON_NAME = "ablation-report.json"
MARKDOWN_NAME = "ablation-report.md"
RULES_HEADING = "## Reading rules"
FEATURE_SENTENCE = ("A block attribution is descriptive only: selection frequency and sign stability "
                    "say what the fitted models used on these data, not why resistance occurs.")

# Contract section 15: the within-margin note (descriptive, cluster contrasts only) and what the pooled p tests.
WITHIN_MARGIN_MIN_CLUSTERS = 60
WITHIN_MARGIN_SENTENCE = ("The 90% interval lies within ±margin. This is descriptive only: GRM's calibration did not "
                          "hold the false-equivalence rate at 5%, so it is not an equivalence finding.")
POOLED_P_SENTENCE = ("Pooled p tests 'no effect in any species'. It is not a test that the average species effect "
                     "is nonzero.")
FAMILY_TITLES = {
    "A": "Family A: confirmatory nested-chain increments, pooled",
    "B": "Family B: leave-one-block-out, pooled",
    "C": "Family C: per-dataset contrasts of the A and B comparisons",
    "D": "Family D: exploratory",
}
FAMILY_NOTES = {
    "A": "Each line pools, across species on r2_ll, an arm of the confirmatory nested chain against the arm "
         "before it, on held-out predictions, with Holm-adjusted p-values within family A. An improvement is "
         "a predictive increment on these datasets only.",
    "B": "Each line pools, across species on r2_ll, the removal of one block (or C+K against C), with "
         "Holm-adjusted p-values within family B. Where a line removes a block, a small leave-one-block-out "
         "drop means the block is redundant given the other blocks, not unimportant.",
    "C": "Per-dataset contrasts of the A and B comparisons on log_loss, with Benjamini-Hochberg-adjusted "
         "p-values across all of them. Where a line removes a block, a small leave-one-block-out drop "
         "means the block is redundant given the other blocks, not unimportant.",
    "D": "Exploratory contrasts: random and source splits, other learners, secondary endpoints and "
         "sensitivity analyses. Family D carries no p-values and decides no direction; development datasets "
         "may appear. Where a line compares splits, the random-vs-lineage gap is a diagnostic, not causal "
         "biology.",
}
RULES = """\
1. "Not significant" means these data did not show a difference at the stated level. It is never read as "no value".
2. A small leave-one-block-out drop means the block is redundant given the other blocks, not unimportant.
3. A block attribution is descriptive only: it says what the fitted models used, not why resistance occurs.
4. The random-vs-lineage gap is a diagnostic, not causal biology. Family D is exploratory and decides no direction.
5. Every interval is printed with its label and means only what its label says. The 95% interval and \
adjusted p-value decide improves or worse; anything else is inconclusive. There is no equivalence decision: GRM's \
calibration did not hold the false-equivalence rate at 5%. The 90% TOST interval is printed and is descriptive \
only.
6. A claim level is printed only beside an improves decision, and contrasts never claim more than \
predictive_increment. Features claim at most association_stable from the pipeline; mechanism_candidate and \
causal_supported come only from named manual annotations with evidence references.
7. This is software output for review. It makes no release or clinical claim.
8. A log-loss improves can come from better calibration or less overfitting, not only from better \
discrimination. The AUROC delta printed beside each log-loss decision is the discrimination change; it is \
descriptive only and decides nothing.
9. Random-split contrasts are descriptive only (family D, no p, no direction). With five outer folds a \
fold-aware test cannot reach p < 0.05, and the cluster-level test is anti-conservative on this split: in GRM's \
calibration over 20 seeds, a learner-seed twin of the default was non-inconclusive for gbt in 9 of 20 on the \
random split, against 1 of 20 for gbt and 0 of 40 for rf on the lineage split (2.0% improves and 3.0% worse \
over all its contrasts). Decisions come from the lineage split only.
10. A sentinel arm, the learner-seed twin of the default, is exchangeable with it. When its lineage-split \
contrast for a learner is anything but inconclusive, or was not produced, every decision of that learner on \
that species is marked unreliable: kept and listed, but left out of every summary of decisions. The sentinel \
probes only the learners that draw randomness, rf and gbt; scm, cart and the lasso draw nothing, so their twin \
contrast is zero by construction. The input comparison sorts every arm's genomes by id, which removes \
genome-order effects only between arms with equal genome sets. It keys every fold on genome identity, through \
the default's cohort and one grouping of the species' eligible genomes that every arm splits on, so a genome both \
arms hold is scored only in its default outer fold and has the default's inner fold and holdout memberships in \
every outer fold it trains in. Neither makes a cohort arm's decision free of the change in its training samples. \
Beside each cohort arm's contrasts is the no-signal base rate measured for its kind (cap down, cap up or \
selection swap), or "not measured".
11. An arm whose genome set equals its baseline's is not run: it is reported as an identical cohort, with its \
reason, and has no contrasts.
"""
REPORT_KEYS = ("schema", "schema_version", "inputs_sha256", "predictions_sha256", "manifest", "planned_contrasts",
               "planned_pooled", "not_produced_reasons", "attempts", "contrast_counts", "pooled_counts",
               "not_produced", "pooled_not_produced", "undefined_contrasts", "unavailable_runs",
               "unavailable_blocks", "identical_cohorts", "exclusions", "unreliable", "contrasts", "pooled",
               "features")


COMPARISON = ("learner", "split", "baseline", "treatment", "endpoint")
# K's records are a few levels deep; deeper input is refused before copy.deepcopy can exhaust the stack.
MAX_JSON_DEPTH = 100


def _comparison(parts):
    """What links a pooled record to its per-dataset inputs (contract section 13): learner,
    split, arms and endpoint, not the family. An A or B pool on r2_ll draws on the family D
    r2_ll contrasts of the same comparison."""
    return tuple(parts[k] for k in COMPARISON)


def _confirmatory_tests(parts):
    """The family C log_loss comparison reported beside an A or B pool (section 13)."""
    return _comparison({**parts, "endpoint": schema.PRIMARY_CONTRAST_ENDPOINT})


def _check_records(contrasts, pooled, features, planned, planned_pooled, not_produced):
    c = schema._Check()
    for name, kind, rows in (("contrasts", "contrast", contrasts), ("pooled", "pooled", pooled),
                             ("features", "feature", features)):
        if not isinstance(rows, list):
            c.fail(name, "must be a list")
            continue
        for i, row in enumerate(rows):
            schema.SPECS[kind](c, f"{name}[{i}]", row)
    schema.strings(pattern=schema.CONTRAST_ID)(c, "planned_contrasts", planned)
    schema.strings(pattern=schema.KEY)(c, "planned_pooled", planned_pooled)
    schema.mapping(schema.text())(c, "not_produced", not_produced)
    c.done()


def _check_adjustment(c, contrasts, pooled, planned, planned_pooled):
    """Section 17: recompute p_adj in families A to C from every planned hypothesis of the family,
    Holm in A and B and BH in C, with a missing or undefined hypothesis at p = 1."""
    by_name = {x["id"]: (f"contrasts[{i}]", x) for i, x in enumerate(contrasts)}
    by_name.update({p["contrast"]: (f"pooled[{j}]", p) for j, p in enumerate(pooled)})
    family_of = {pid: schema.parse_key(schema.contrast_key(pid))["family"] for pid in planned}
    family_of.update({key: schema.parse_key(key)["family"] for key in planned_pooled})
    for family in schema.CONFIRMATORY_FAMILIES:
        method = schema.FAMILY_ADJUST[family]
        names = [name for name, f in family_of.items() if f == family]
        found = [by_name.get(name, (None, None)) for name in names]
        pvalues = [1.0 if rec is None or rec["p"] is None else rec["p"] for _, rec in found]
        for (path, rec), expected in zip(found, schema.ADJUSTERS[method](pvalues)):
            if rec is not None and rec["p_adj"] is not None and abs(rec["p_adj"] - expected) > schema.ADJUST_TOLERANCE:
                c.fail(f"{path}.p_adj", f"is {rec['p_adj']!r}, but {method} over the {len(names)} planned hypotheses "
                                        f"of family {family}, missing or undefined at p = 1, gives {expected!r}")


def _check_pairing(c, path, x, base, treat):
    """A contrast's pairing rule is its treatment run's (same_plan when the run names none),
    and a run paired by another rule pairs only with the reference arm it names. Two arms
    each paired to a reference, or a paired arm against another baseline, are refused."""
    if base is None or treat is None:
        return
    rule = (treat.get("pairing") or {}).get("rule", schema.SAME_PLAN)
    if x["pairing"] != rule:
        c.fail(f"{path}.pairing", f"is {x['pairing']}, but run {treat['id']} pairs by {rule}")
    if base.get("pairing") is not None:
        c.fail(path, f"baseline arm {x['baseline']} is paired to {base['pairing']['reference_arm']}; "
                     "a contrast compares an arm with its own reference arm")
    if treat.get("pairing") is not None and treat["pairing"]["reference_arm"] != x["baseline"]:
        c.fail(path, f"arm {x['treatment']} is paired to {treat['pairing']['reference_arm']}, not {x['baseline']}")


def _check_links(manifest, contrasts, pooled, features, planned, planned_pooled, not_produced):
    """Cross-record checks; returns the planned contrasts and the planned pooled keys with no
    record, each with its reason."""
    c = schema._Check()
    datasets = {d["id"]: d for d in manifest["datasets"]}
    blocks = {b["id"] for b in manifest["blocks"]}
    runs = {(u["dataset"], u["arm"], u["learner"], u["split"]): u for u in manifest["runs"]}
    planned_set = set(planned)
    planned_by_comparison = defaultdict(set)
    for i, pid in enumerate(planned):
        dataset, parts = schema.contrast_dataset(pid), schema.parse_key(schema.contrast_key(pid))
        planned_by_comparison[_comparison(parts)].add(dataset)
        if dataset not in datasets:
            c.fail(f"planned_contrasts[{i}]", f"unknown dataset {dataset}")
        if parts["family"] not in schema.CONTRAST_FAMILIES:
            c.fail(f"planned_contrasts[{i}]", "families A and B hold pooled records only; plan their per-dataset "
                                              "contrasts as family C")
        elif parts["family"] == "C" and dataset in datasets and datasets[dataset]["partition"] != "confirmatory":
            c.fail(f"planned_contrasts[{i}]", f"development dataset {dataset} is refused in family C")

    schema._unique(c, "contrasts", contrasts, lambda x: x["id"], "contrast id")
    members = defaultdict(dict)
    # A development dataset in family C is refused on its planned id, which every contrast must have.
    for i, x in enumerate(contrasts):
        path, dataset = f"contrasts[{i}]", x["dataset"]
        members[_comparison(x)][dataset] = x
        if x["id"] not in planned_set:
            c.fail(f"{path}.id", "is not in planned_contrasts")
        if dataset not in datasets:
            c.fail(f"{path}.dataset", f"unknown dataset {dataset}")
            continue
        for arm in (x["baseline"], x["treatment"]):
            run = runs.get((dataset, arm, x["learner"], x["split"]))
            if run is None or run["status"] != "ok":
                c.fail(path, f"arm {arm} has no ok run for {dataset}, {x['learner']}, {x['split']}")
        _check_pairing(c, path, x, runs.get((dataset, x["baseline"], x["learner"], x["split"])),
                       runs.get((dataset, x["treatment"], x["learner"], x["split"])))

    pooled_planned_set = set(planned_pooled)
    for i, key in enumerate(planned_pooled):
        parts = schema.parse_key(key)
        if parts["family"] not in schema.POOLED_FAMILIES:
            c.fail(f"planned_pooled[{i}]", "family C holds per-dataset contrasts only; pooled keys are A, B or D")
        elif parts["family"] != "D" and parts["endpoint"] != schema.PRIMARY_POOLED_ENDPOINT:
            c.fail(f"planned_pooled[{i}]", f"family {parts['family']} pools the primary endpoint "
                                           f"{schema.PRIMARY_POOLED_ENDPOINT} only")

    produced = {x["id"] for x in contrasts}
    produced_pooled = {p["contrast"] for p in pooled}
    for pid in not_produced:
        if pid not in planned_set and pid not in pooled_planned_set:
            c.fail(f"not_produced.{pid}", "is not a planned contrast or pooled key")
        elif pid in produced or pid in produced_pooled:
            c.fail(f"not_produced.{pid}", "has a record")
    pooled_reasons = {(_comparison(schema.parse_key(p["contrast"])), e["dataset"]): e["reason"]
                      for p in pooled for e in p["excluded"]}
    missing = []
    for i, pid in enumerate(planned):
        if pid in produced:
            continue
        comparison = _comparison(schema.parse_key(schema.contrast_key(pid)))
        reason = not_produced.get(pid) or pooled_reasons.get((comparison, schema.contrast_dataset(pid)))
        if reason is None:
            c.fail(f"planned_contrasts[{i}]", f"{pid} has no record and no reason; give one in not_produced")
        else:
            missing.append({"id": pid, "reason": reason})

    pooled_missing = []
    for i, key in enumerate(planned_pooled):
        if key in produced_pooled:
            continue
        if key in not_produced:
            pooled_missing.append({"contrast": key, "reason": not_produced[key]})
        else:
            c.fail(f"planned_pooled[{i}]", f"{key} has no pooled record and no reason; give one in not_produced")

    schema._unique(c, "pooled", pooled, lambda p: p["contrast"], "contrast")
    for j, p in enumerate(pooled):
        path, confirmatory = f"pooled[{j}]", p["family"] in schema.CONFIRMATORY_FAMILIES
        if p["contrast"] not in pooled_planned_set:
            c.fail(f"{path}.contrast", "is not in planned_pooled")
        source = _comparison(schema.parse_key(p["contrast"]))
        scope, rows = planned_by_comparison.get(source, set()), members.get(source, {})
        if confirmatory:  # development datasets never enter an A or B denominator (section 2)
            scope = {d for d in scope if d in datasets and datasets[d]["partition"] == "confirmatory"}
        excluded = {e["dataset"] for e in p["excluded"]}
        if not scope:
            c.fail(f"{path}.contrast", f"pools a comparison with no planned per-dataset contrasts ({':'.join(source)})")
            continue
        if p["n_datasets"] != len(scope):
            c.fail(f"{path}.n_datasets", f"must count every planned dataset of the comparison, {len(scope)}")
        for d in sorted(excluded - scope):
            if confirmatory and d in datasets and datasets[d]["partition"] != "confirmatory":
                c.fail(f"{path}.excluded", f"development dataset {d} is refused in family {p['family']}")
            else:
                c.fail(f"{path}.excluded", f"dataset {d} is not planned for this comparison")
        included = scope - excluded
        for d in sorted(included):
            if d not in rows:
                c.fail(f"{path}.excluded", f"dataset {d} has no per-dataset record; exclude it with a reason")
            elif rows[d]["estimate"] is None:
                c.fail(f"{path}.excluded", f"undefined contrast {rows[d]['id']} must be excluded with its reason")
        species = {datasets[d]["species"] for d in included if d in datasets}
        if p["n_units"] != len(species):
            c.fail(f"{path}.n_units", f"must equal the {len(species)} contributing species")

    schema._unique(c, "features", features, lambda f: (f["dataset"], f["column"]), "(dataset, column)")
    for i, f in enumerate(features):
        if f["dataset"] not in datasets:
            c.fail(f"features[{i}].dataset", f"unknown dataset {f['dataset']}")
        for b in f["blocks"]:
            if b not in blocks:
                c.fail(f"features[{i}].blocks", f"unknown block {b}")
    _check_adjustment(c, contrasts, pooled, planned, planned_pooled)
    c.done()
    return missing, pooled_missing


def check_predictions(manifest, root):
    """Read every ok run's predictions file (JSON Lines, relative to root) and
    validate the rows as one set against the manifest. Each file holds rows of
    the ok runs that name it, and each such run has rows for every repeat.
    Returns {path: sha256 of the file bytes}."""
    c = schema._Check()
    root = Path(root).resolve()
    by_path = defaultdict(list)
    for u in manifest["runs"]:
        if u["status"] == "ok":
            by_path[u["predictions"]].append(u)
    rows, digests = [], {}
    for name, runs in sorted(by_path.items()):
        path = (root / name).resolve()
        if root not in path.parents:
            c.fail(name, "must stay inside the prediction root")
            continue
        try:
            data = path.read_bytes()
            lines = [_loads(line) for line in data.decode("utf-8").splitlines() if line.strip()]
        except (OSError, ValueError, RecursionError) as e:
            c.fail(name, str(e))
            continue
        digests[name] = hashlib.sha256(data).hexdigest()
        keys = {(u["dataset"], u["arm"], u["learner"], u["split"]): u for u in runs}
        seen = defaultdict(set)
        for n, row in enumerate(lines):
            before = len(c.errors)
            schema.PREDICTION(c, f"{name}[{n}]", row)
            if len(c.errors) != before:
                continue
            key = (row["dataset"], row["arm"], row["learner"], row["split"])
            if key not in keys:
                c.fail(f"{name}[{n}]", "belongs to no ok run that names this file")
            else:
                seen[key].add(row["repeat"])
                rows.append(row)
        for key, u in keys.items():
            if seen[key] != set(range(u["repeats"])):
                c.fail(name, f"run {u['id']} needs rows for repeats 0 to {u['repeats'] - 1}")
    c.done()
    schema.validate_predictions(rows, manifest)
    return digests


def build_report(manifest, contrasts, pooled, features, planned_contrasts, planned_pooled, not_produced=None,
                 out_dir=None, prediction_root=None):
    """Validate every input; report every planned contrast and planned pooled
    key as a record or as not produced with a reason; recompute p_adj in
    families A to C; refuse development datasets in families A to C; with
    prediction_root, also check every ok run's prediction rows. Return the
    report, and with out_dir also write its JSON and Markdown.

    planned_contrasts and planned_pooled are trusted input from the protocol:
    K checks every record, denominator and adjustment against them, but cannot
    tell a list with a hypothesis removed from the protocol's own list. Compare
    inputs_sha256.planned_contrasts and inputs_sha256.planned_pooled with the
    protocol's."""
    arguments = {"manifest": manifest, "contrasts": contrasts, "pooled": pooled, "features": features,
                 "planned_contrasts": planned_contrasts, "planned_pooled": planned_pooled,
                 "not_produced": not_produced}
    deep = [name for name, value in arguments.items() if _too_deep(value)]
    if deep:
        raise schema.SchemaError([f"build_report: {name} nested too deeply (more than {MAX_JSON_DEPTH} levels)"
                                  for name in deep])
    try:  # _too_deep counts dicts, lists and tuples; another deep container (a deque, a frozenset) gets here
        copied = copy.deepcopy(tuple(arguments.values()))
    except RecursionError:
        too_deep = f"build_report: input nested too deeply (more than {MAX_JSON_DEPTH} levels)"
        raise schema.SchemaError([too_deep]) from None
    manifest, contrasts, pooled, features, planned, planned_pooled, not_produced = copied
    manifest = schema.validate_manifest(manifest)
    not_produced = not_produced if not_produced is not None else {}
    _check_records(contrasts, pooled, features, planned, planned_pooled, not_produced)
    missing, pooled_missing = _check_links(manifest, contrasts, pooled, features, planned, planned_pooled,
                                           not_produced)
    predictions = None if prediction_root is None else check_predictions(manifest, prediction_root)
    runs = manifest["runs"]
    inputs = {"manifest": manifest, "contrasts": contrasts, "pooled": pooled, "features": features,
              "planned_contrasts": planned, "planned_pooled": planned_pooled, "not_produced": not_produced}
    undefined = [{"id": x["id"], "reason": x["undefined_reason"]}
                 for x in contrasts if x["undefined_reason"] is not None]
    report = {
        "schema": schema.SCHEMA,
        "schema_version": schema.SCHEMA_VERSION,
        "inputs_sha256": {name: schema.sha256(value, name) for name, value in inputs.items()},
        "predictions_sha256": predictions,
        "manifest": manifest,
        "planned_contrasts": planned,
        "planned_pooled": planned_pooled,
        "not_produced_reasons": not_produced,
        "attempts": {
            "runs": len(runs),
            "ok": sum(u["status"] == "ok" for u in runs),
            "unavailable": sum(u["status"] == "unavailable" for u in runs),
            "identical_cohort": sum(u["status"] == "identical_cohort" for u in runs),
        },
        "contrast_counts": {
            "planned": len(planned),
            "produced": len(contrasts),
            "defined": len(contrasts) - len(undefined),
            "undefined": len(undefined),
            "not_produced": len(missing),
        },
        "pooled_counts": {
            "planned": len(planned_pooled),
            "produced": len(pooled),
            "not_produced": len(pooled_missing),
        },
        "not_produced": missing,
        "pooled_not_produced": pooled_missing,
        "undefined_contrasts": undefined,
        "unavailable_runs": [
            {k: u[k] for k in ("id", "dataset", "arm", "learner", "split", "reason")}
            for u in runs if u["status"] == "unavailable"
        ],
        "identical_cohorts": [
            {k: u[k] for k in ("id", "dataset", "arm", "learner", "split", "reason")}
            for u in runs if u["status"] == "identical_cohort"
        ],
        "unavailable_blocks": [
            {"block": b["id"], "reason": b["unavailable_reason"]}
            for b in manifest["blocks"] if b["unavailable_reason"] is not None
        ],
        "exclusions": [
            {"pooled": p["contrast"], "dataset": e["dataset"], "reason": e["reason"]}
            for p in pooled for e in p["excluded"]
        ],
        "unreliable": _unreliable(manifest, contrasts, planned, missing),
        "contrasts": contrasts,
        "pooled": pooled,
        "features": features,
    }
    if out_dir is not None:
        write_report(report, out_dir)
    return report


def validate_report(report, prediction_root=None):
    """Check a loaded report: its inputs validate, and every section, compared
    as canonical JSON so true is not 1 and 16.0 is not 16, matches what
    build_report makes from them. A consistent edit of the inputs and every
    derived section together passes; compare inputs_sha256 with the canonical
    hashes of the original input files to catch that.

    predictions_sha256 is re-derived only with prediction_root. Without it, a
    recorded value is checked for shape (one hash per ok run's file), not
    against the files."""
    if not isinstance(report, dict) or set(report) != set(REPORT_KEYS):
        raise schema.SchemaError(["report: must have exactly the fields " + ", ".join(REPORT_KEYS)])
    if _too_deep(report):
        raise schema.SchemaError([f"report: nested too deeply (more than {MAX_JSON_DEPTH} levels)"])
    if type(report["schema_version"]) is not int or (report["schema"], report["schema_version"]) != (
            schema.SCHEMA, schema.SCHEMA_VERSION):
        raise schema.SchemaError([f"report: must be {schema.SCHEMA} schema {schema.SCHEMA_VERSION}"])
    rebuilt = build_report(report["manifest"], report["contrasts"], report["pooled"], report["features"],
                           report["planned_contrasts"], report["planned_pooled"], report["not_produced_reasons"],
                           prediction_root=prediction_root)
    recorded = report["predictions_sha256"]
    if prediction_root is None and recorded is not None:
        c = schema._Check()
        schema.mapping(schema.text(schema.HEX64))(c, "report.predictions_sha256", recorded)
        files = {u["predictions"] for u in rebuilt["manifest"]["runs"] if u["status"] == "ok"}
        if not c.errors and set(recorded) != files:
            c.fail("report.predictions_sha256", "must name every ok run's predictions file")
        c.done()
        rebuilt["predictions_sha256"] = recorded
    for k in REPORT_KEYS:
        if schema.canonical(rebuilt[k], f"report.{k}") != schema.canonical(report[k], f"report.{k}"):
            raise schema.SchemaError([f"report.{k}: does not match its records"])
    return report


def dumps(report):
    return json.dumps(report, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False) + "\n"


def _umask():
    mask = os.umask(0o022)
    os.umask(mask)
    return mask


def _reserve(out, final, suffix):
    """A fresh empty file named after final in out, closed; mkstemp picks a name nothing else holds. out is
    resolved first, as the kernel resolves it for the renames: mkstemp alone resolves ".." lexically, which
    after a symlinked folder names another folder (RKC-L1)."""
    fd, name = tempfile.mkstemp(dir=os.path.realpath(out), prefix=f".{final.name}.", suffix=suffix)
    os.close(fd)
    return Path(name)


def write_report(report, out_dir):
    """Section 17: render both files, then write them as UTF-8 to temporary names in out_dir,
    with the mode a plain open() would give (0o666 & ~umask). Move any earlier final files to
    backup names in out_dir, rename both new files into place, then delete the backups. On any
    failure remove the new files, restore the backups, so an earlier report survives byte for
    byte, remove the folders this call made for out_dir while empty, and raise
    SchemaError(["<path>: <message>"]). A final path that exists but isn't a regular file is
    refused before anything is written."""
    out = Path(out_dir)
    finals = [out / JSON_NAME, out / MARKDOWN_NAME]
    current, temps, reserved, backups, placed, created = out, [], [], [], [], []
    try:
        texts = [dumps(report), render_markdown(report)]
        payloads = []
        for current, text in zip(finals, texts):
            payloads.append((current, text.encode("utf-8")))
        for current in finals:
            if (current.exists() or current.is_symlink()) and not current.is_file():
                raise IsADirectoryError(f"exists and is not a regular file; {current.name} not replaced")
        current = out
        _make_dirs(out, created)
        mode = 0o666 & ~_umask()
        for current, data in payloads:
            temps.append(_reserve(out, current, ".tmp"))
            with open(temps[-1], "wb") as handle:
                os.fchmod(handle.fileno(), mode)
                handle.write(data)
        for current in finals:
            if current.exists() or current.is_symlink():
                reserved.append(_reserve(out, current, ".bak"))
                os.replace(current, reserved[-1])
                backups.append((current, reserved[-1]))
        for (current, _), temp in zip(payloads, temps):
            os.replace(temp, current)
            placed.append(current)
    except (OSError, UnicodeError) as e:
        for path in temps + placed:
            with contextlib.suppress(OSError):
                path.unlink()
        kept = []
        for final, backup in backups:
            try:
                os.replace(backup, final)
            except OSError:
                kept.append(backup)  # the only copy of that earlier file: never delete it
        for path in reserved:
            if path not in kept:
                with contextlib.suppress(OSError):
                    path.unlink()
        for folder in reversed(created):
            with contextlib.suppress(OSError):  # only while empty; through "..", the next may be a sibling
                folder.rmdir()
        message = f"{current}: {e}"
        if kept:
            message += "; earlier report kept at " + ", ".join(str(path) for path in kept)
        raise schema.SchemaError([message]) from None
    for _, backup in backups:
        with contextlib.suppress(OSError):
            backup.unlink()


def _make_dirs(path, created):
    """path.mkdir(parents=True, exist_ok=True) without recursion, appending each folder its os.mkdir calls
    made to created, shallowest first. Only those are ever removed again: a path through "..", or a folder
    that already existed, is never recorded."""
    pending = [(path, False)]
    while pending:
        folder, parent_made = pending.pop()
        try:
            os.mkdir(folder)
        except FileNotFoundError:
            if parent_made or folder.parent == folder:
                raise
            pending += [(folder, True), (folder.parent, False)]
        except OSError:
            if not folder.is_dir():
                raise
        else:
            created.append(folder)


def _unique_keys(pairs):
    keys = [k for k, _ in pairs]
    repeated = sorted({k for k in keys if keys.count(k) > 1})
    if repeated:
        raise ValueError(f"duplicate keys {repeated}")
    return dict(pairs)


def _too_deep(value):
    """True when value nests lists and objects more than MAX_JSON_DEPTH levels (no recursion). A dict's keys
    count as nested in it, as its values do. Other containers aren't walked; build_report's deep copy
    refuses those when deep."""
    stack = [(value, 1)]
    while stack:
        v, depth = stack.pop()
        if isinstance(v, (dict, list, tuple)):
            if depth > MAX_JSON_DEPTH:
                return True
            stack.extend((x, depth + 1) for x in ((*v.keys(), *v.values()) if isinstance(v, dict) else v))
    return False


def _loads(text):
    """json.loads that refuses a repeated key in any object (RFC 8259 leaves them undefined) and
    input nested more than MAX_JSON_DEPTH levels."""
    too_deep = f"input nested too deeply (more than {MAX_JSON_DEPTH} levels)"
    try:
        value = json.loads(text, object_pairs_hook=_unique_keys)
    except RecursionError:  # the decoder's own limit, near 10,000 levels
        raise ValueError(too_deep) from None
    if _too_deep(value):
        raise ValueError(too_deep)
    return value


def load_report(path, prediction_root=None):
    return validate_report(_load(path), prediction_root)


# Markdown summary. Free text is escaped; ids go in code spans.

_BREAKS = re.compile(r"[\r\n\v\f\x1c-\x1e\x85  ]+")
_SPECIAL = re.compile(r"([\\`*_\[\]<>|#~&])")


def _md(value):
    return _SPECIAL.sub(r"\\\1", _BREAKS.sub(" ", str(value)))


def _code(value):
    """A code span fenced with one more backtick than the longest run inside it."""
    value = _BREAKS.sub(" ", str(value))
    if "`" not in value:
        return f"`{value}`"
    fence = "`" * (max(len(run) for run in re.findall(r"`+", value)) + 1)
    return f"{fence} {value} {fence}"


def _plural(n, word):
    return f"{n} {word}" + ("" if n == 1 else "s")


def _num(v):
    text = f"{v:+.4f}"
    return "0.0000" if text in ("+0.0000", "-0.0000") else text


def _p(v):
    return f"{v:.3g}"


def _units(n, unit):
    return f"{n} {_md(unit)}" + ("" if n == 1 or unit.endswith("s") else "s")


def _intervals(rec):
    ci, tost = rec["ci"], rec["ci_tost"]
    return (f"95% CI [{_num(ci[0])}, {_num(ci[1])}], 90% TOST interval [{_num(tost[0])}, {_num(tost[1])}] "
            f"({_md(rec['interval_label'])})")


def _tests(rec, family):
    if rec["p"] is None:
        return f"no p-value (family {family} is exploratory, adjust none); margin ±{rec['margin']:g}"
    return f"p = {_p(rec['p'])}, p_adj ({rec['adjust']}) = {_p(rec['p_adj'])}; margin ±{rec['margin']:g}"


def _removes_block(arm_blocks, treatment, baseline):
    """True when the treatment arm's blocks are a strict subset of the baseline's (from the manifest's runs)."""
    if treatment not in arm_blocks or baseline not in arm_blocks:
        return False
    return arm_blocks[treatment] < arm_blocks[baseline]


def within_margin(rec):
    """Contract section 15: a defined inconclusive cluster contrast on at least 60 clusters whose
    90% interval lies strictly inside the margin. Pooled records (unit species) never qualify."""
    tost, margin = rec["ci_tost"], rec["margin"]
    return (rec["undefined_reason"] is None and rec["unit"] == "cluster"
            and rec["n_units"] >= WITHIN_MARGIN_MIN_CLUSTERS and rec["decision"] == "inconclusive"
            and -margin < tost[0] and tost[1] < margin)


def _decision(rec, removes, treatment, baseline):
    d, margin = rec["decision"], rec["margin"]
    if d in ("improves", "worse"):
        text = f"{treatment} improves on {baseline}" if d == "improves" else f"{treatment} is worse than {baseline}"
        if abs(rec["estimate"]) < margin:
            text += f", significant, but smaller than the margin of {margin:g}"
    elif rec["ci"][0] > 0 or rec["ci"][1] < 0:
        if rec["p_adj"] is None:
            text = "inconclusive: the 95% interval excludes zero, but exploratory family D decides no direction"
        else:
            text = (f"inconclusive: the 95% interval excludes zero, but the adjusted p-value "
                    f"({_p(rec['p_adj'])}) does not meet 0.05")
    else:
        text = "not significant: the 95% interval includes zero"
    if removes and d != "worse":
        text += "; a small drop means the block is redundant given the others, not unimportant"
    text += "."
    if d == "improves":
        text += f" Claim level {schema.CONTRAST_LEVEL}."
    return text


def _design(x):
    return (f"{_md(x['treatment'])} vs {_md(x['baseline'])}, {_md(x['learner'])}, {_md(x['split'])} split, "
            f"{_md(x['endpoint'])}")


RANDOM_SPLIT_SENTENCE = ("Descriptive only, with no p and no direction: with five outer folds a fold-aware test "
                         "cannot reach p < 0.05, and the cluster-level test is anti-conservative on this split "
                         "(a learner-seed twin's gbt came out non-inconclusive in 9 of 20 seeds).")
UNRELIABLE_SENTENCE = ("The learner-seed sentinel, which is exchangeable with its baseline, was not inconclusive "
                       "on the lineage split, so every decision of these learners on these species is unreliable: "
                       "kept below and flagged, and left out of every summary of decisions.")


def _unreliable(manifest, contrasts, planned, missing):
    """Each (species, learner) a sentinel arm flags, with the sentinel contrast, the reason and the
    contrasts it makes unreliable (rule 10)."""
    sentinels = set(manifest.get("sentinel_arms", []))
    if not sentinels:
        return []
    species = {d["id"]: d["species"] for d in manifest["datasets"]}
    by_id = {x["id"]: x for x in contrasts}
    absent = {m["id"]: m["reason"] for m in missing}
    flagged = {}
    for pid in planned:
        parts = schema.parse_key(schema.contrast_key(pid))
        if (parts["treatment"] not in sentinels or parts["split"] != schema.SENTINEL_SPLIT
                or parts["endpoint"] != schema.PRIMARY_CONTRAST_ENDPOINT):
            continue
        x = by_id.get(pid)
        if x is None:
            reason = f"the sentinel contrast was not produced: {absent.get(pid, 'no record')}"
        elif x["decision"] != "inconclusive":
            reason = f"the sentinel contrast is {x['decision']}"
        else:
            continue
        key = (species[schema.contrast_dataset(pid)], parts["learner"])
        flagged.setdefault(key, (pid, reason))
    out = []
    for (sp, learner), (pid, reason) in sorted(flagged.items()):
        ids = sorted(x["id"] for x in contrasts
                     if x["learner"] == learner and species[x["dataset"]] == sp)
        out.append({"species": sp, "learner": learner, "sentinel": pid, "reason": reason, "contrasts": ids})
    return out


def _contrast_line(x, datasets, arm_blocks, base_rates=None):
    base_rates = base_rates or {}
    dataset = datasets[x["dataset"]]
    head = f"- {_code(x['id'])} ({_md(dataset['species'])}, {dataset['partition']}): {_design(x)}: "
    method = f"{_md(x['ci_method'])}, n = {_units(x['n_units'], x['unit'])}"
    if x["undefined_reason"] is not None:
        line = (head + f"undefined: {_md(x['undefined_reason'])}. No estimate, interval or test ({method}). "
                f"Interval label: {_md(x['interval_label'])}.")
    else:
        removes = _removes_block(arm_blocks, x["treatment"], x["baseline"])
        line = (head + f"estimate {_num(x['estimate'])}, {_intervals(x)} ({method}); {_tests(x, x['family'])}. "
                f"Decision: {_decision(x, removes, _md(x['treatment']), _md(x['baseline']))}")
        if within_margin(x):
            line += f" {WITHIN_MARGIN_SENTENCE}"
    if x["few_clusters"]:
        line += (f" Few clusters (n < {schema.FEW_CLUSTERS}): the intervals are approximate; 95% coverage fell "
                 "below nominal in calibration.")
    if x["undefined_resamples"]:
        line += f" Undefined resamples: {x['undefined_resamples']}."
    rate = base_rates.get(x["treatment"])
    if rate is not None:
        kind = f" ({_md(rate['kind'])})" if rate.get("kind") else ""
        if rate["rate"] is None:
            line += f" No-signal base rate for this arm type{kind}: not measured."
        else:
            interval = (f" (95% interval [{rate['ci95'][0] * 100:.1f}, {rate['ci95'][1] * 100:.1f}]%)"
                        if rate["ci95"] else "")
            line += (f" No-signal base rate for this arm type{kind}: {rate['rate'] * 100:.1f}% of lineage-split "
                     f"contrasts non-inconclusive{interval} over {rate['seeds']} seeds ({_md(rate['source'])}).")
    if x["auroc_delta"] is not None:
        lo, hi = x["auroc_ci"]
        line += (f" AUROC delta {_num(x['auroc_delta'])} (95% interval [{_num(lo)}, {_num(hi)}]), descriptive "
                 "only.")
    if x["pairing"] != schema.SAME_PLAN:
        line += f" Pairing: {_code(x['pairing'])}."
    return line


def _few_species(p):
    """Contract section 13: with S < 6 species the exact sign-flip p cannot reach 0.05."""
    s = p["n_units"]
    if s >= schema.MIN_TESTABLE_SPECIES:
        return ""
    smallest = f"; the smallest p is 2/2^S = {_p(2 / 2 ** s)}" if s else ""
    return (f" A pooled p below 0.05 is not reachable with fewer than {schema.MIN_TESTABLE_SPECIES} species "
            f"(S = {s}{smallest}).")


def _pooled_line(p, arm_blocks, tests):
    """tests: the family C log_loss contrasts of the same comparison, beside an A or B pool."""
    included = p["n_datasets"] - p["n_excluded"]
    design = schema.parse_key(p["contrast"])
    removes = _removes_block(arm_blocks, design["treatment"], design["baseline"])
    head = f"- {_code(p['contrast'])} (family {p['family']}: {_design(design)})"
    if included == 0:
        body = (f"no pooled estimate; all {_plural(p['n_datasets'], 'dataset')} in the denominator were excluded "
                f"(adjust {p['adjust']}; margin ±{p['margin']:g}). Interval label: {_md(p['interval_label'])}.")
    else:
        body = (f"estimate {_num(p['estimate'])}, {_intervals(p)}, n = {_units(p['n_units'], 'species')}; "
                f"{_tests(p, p['family'])}. "
                f"Decision: {_decision(p, removes, _md(design['treatment']), _md(design['baseline']))} "
                f"Denominator: {_plural(p['n_datasets'], 'dataset')}, {included} contributing, "
                f"{p['n_excluded']} excluded.")
    if p["p"] is not None:
        body += f" {POOLED_P_SENTENCE}"
    body += _few_species(p)
    if p["family"] in ("A", "B"):
        listed = "; ".join(f"{_code(x['id'])} {x['decision']}" for x in tests) or "none produced"
        body += f" Per-dataset confirmatory tests (family C, log_loss): {listed}."
    for e in p["excluded"]:
        body += f"\n  - excluded {_code(e['dataset'])}: {_md(e['reason'])}"
    return f"{head}: {body}"


def _feature_claim(f):
    text = f["claim_level"]
    notes = f.get("annotations")
    if notes:
        steps = "; ".join(f"{_md(n['by'])} on {n['date']}, {n['previous_level']} to {n['level']}" for n in notes)
        text += f" (manual annotations: {steps}; evidence: {', '.join(_md(r) for r in f['evidence_refs'])})"
    return text


def _synthetic_line(m):
    """The header's synthetic line, saying what the revision-12 gate did (main's ruling on REVIEW-INTEG7 L2)."""
    gate = m.get("revision12_gate")
    if not m["synthetic"]:
        return "- Synthetic run: no. The revision-12 gate passed."
    if gate == "passed":
        return "- Synthetic run: yes. The spec also carries a revision-12 gate, and it passed."
    if gate == "not_required":
        return ("- Synthetic run: yes, on a synthetic data root. The spec carries no revision-12 gate, so none ran; "
                "a real run needs it.")
    return ("- Synthetic run: yes. The spec says it is synthetic, so the revision-12 gate a real run needs was not "
            "required; this manifest doesn't record whether one ran.")


def render_markdown(report):
    m = report["manifest"]
    datasets = {d["id"]: d for d in m["datasets"]}
    arm_blocks = {u["arm"]: set(u["blocks"]) for u in m["runs"]}
    confirmatory_tests = defaultdict(list)
    for x in report["contrasts"]:
        if x["family"] == "C":
            confirmatory_tests[_comparison(x)].append(x)
    a, n, pc = report["attempts"], report["contrast_counts"], report["pooled_counts"]
    out = [
        f"# Ablation report ({report['schema']} schema {report['schema_version']})",
        "",
        "Software output for review. It makes no release or clinical claim.",
        "",
        f"- Protocol sha256: `{m['protocol_sha256']}`",
        f"- Source commit: `{m['source_commit']}`",
        f"- Binary sha256: `{m['binary_sha256']}`",
        *([_synthetic_line(m)] if "synthetic" in m else []),
        *([f"- Species pool: the union of the cohorts of this run's {len(m['species_pool']['arms'])} arms "
           f"({m['species_pool']['genomes']} genomes; structure sha256 `{m['species_pool']['structure_sha256']}`). "
           "Lineages and the fold keys of units outside the default depend on this union, so adding an arm "
           "means rerunning the whole comparison."] if m.get("species_pool") else []),
        "- Input sha256 (canonical JSON): " + ", ".join(f"{k} `{v}`" for k, v in report["inputs_sha256"].items()),
        "- Prediction files: " + ("not checked" if report["predictions_sha256"] is None else
                                  f"{_plural(len(report['predictions_sha256']), 'file')} checked against the "
                                  "manifest, sha256 " + ", ".join(f"{_code(k)} `{v}`" for k, v in
                                                                  sorted(report["predictions_sha256"].items()))),
        "",
        "## Attempts",
        "",
        f"{a['runs']} runs attempted: {a['ok']} ok, {a['unavailable']} unavailable.",
        "",
    ]
    if report["identical_cohorts"]:
        out += ["Identical cohorts, not run (rule 11):", ""]
        out += [f"- {_code(u['id'])} ({_md(u['dataset'])}, arm {_md(u['arm'])}, {_md(u['learner'])}, "
                f"{_md(u['split'])} split): {_md(u['reason'])}" for u in report["identical_cohorts"]]
        out.append("")
    if report["unavailable_runs"]:
        out += ["Unavailable runs:", ""]
        out += [f"- {_code(u['id'])} ({_code(u['dataset'])}, arm {_md(u['arm'])}, {_md(u['learner'])}, "
                f"{_md(u['split'])} split): {_md(u['reason'])}" for u in report["unavailable_runs"]]
        out.append("")
    if report["unavailable_blocks"]:
        out += ["Unavailable blocks:", ""]
        out += [f"- {_code(b['block'])}: {_md(b['reason'])}" for b in report["unavailable_blocks"]]
        out.append("")
    out += [f"{n['planned']} contrasts planned: {n['produced']} produced ({n['defined']} with an estimate, "
            f"{n['undefined']} undefined), {n['not_produced']} not produced.", "",
            f"{_plural(pc['planned'], 'pooled record')} planned: {pc['produced']} produced, {pc['not_produced']} "
            "not produced (each listed in its family's section below).", ""]
    if report["not_produced"]:
        out += ["Not produced:", ""]
        out += [f"- {_code(x['id'])}: not produced: {_md(x['reason'])}" for x in report["not_produced"]]
        out.append("")
    if report["undefined_contrasts"]:
        out += ["Undefined:", ""]
        out += [f"- {_code(x['id'])}: {_md(x['reason'])}" for x in report["undefined_contrasts"]]
        out.append("")
    if report["exclusions"]:
        out += ["Pooling exclusions (each still counts in its denominator):", ""]
        out += [f"- {_code(e['pooled'])}: {_code(e['dataset'])}: {_md(e['reason'])}" for e in report["exclusions"]]
        out.append("")

    out += ["## Datasets", "", "| Dataset | Species | Partition | Samples | Groups |", "|---|---|---|---:|---|"]
    for d in m["datasets"]:
        groups = ", ".join(f"{_md(k)} {v}" for k, v in sorted(d["group_counts"].items())) or "none"
        out.append(f"| {_code(d['id'])} | {_md(d['species'])} | {d['partition']} | {d['n_samples']} | {groups} |")
    out.append("")

    out += ["## Pooled across species", ""]
    for family in schema.POOLED_FAMILIES:
        rows = [p for p in report["pooled"] if p["family"] == family]
        absent = [x for x in report["pooled_not_produced"] if schema.parse_key(x["contrast"])["family"] == family]
        if rows or absent:
            out += [f"### {FAMILY_TITLES[family]}", "", FAMILY_NOTES[family], ""]
            out += [_pooled_line(p, arm_blocks, confirmatory_tests.get(
                _confirmatory_tests(schema.parse_key(p["contrast"])), [])) for p in rows]
            out += [f"- {_code(x['contrast'])} (family {family}: {_design(schema.parse_key(x['contrast']))}): "
                    f"not produced: {_md(x['reason'])}." for x in absent]
            out.append("")
    if not report["planned_pooled"]:
        out += ["No pooled records planned.", ""]

    if report["unreliable"]:
        out += ["## Unreliable learners", "", UNRELIABLE_SENTENCE, ""]
        out += [f"- {_md(u['species'])}, {_code(u['learner'])}: sentinel {_code(u['sentinel'])}: "
                f"{_md(u['reason'])}; {_plural(len(u['contrasts']), 'contrast')} flagged."
                for u in report["unreliable"]]
        out.append("")
    unreliable = {i for u in report["unreliable"] for i in u["contrasts"]}
    rates = {b["arm"]: b for b in report["manifest"].get("base_rates", [])}

    def line(x):
        text = _contrast_line(x, datasets, arm_blocks, rates)
        return text + " Unreliable (rule 10)." if x["id"] in unreliable else text

    descriptive = [x for x in report["contrasts"]
                   if x["family"] == "D" and x["split"] in schema.DESCRIPTIVE_SPLITS]
    out += ["## Per-dataset contrasts", ""]
    for family in schema.CONTRAST_FAMILIES:
        rows = [x for x in report["contrasts"] if x["family"] == family and x not in descriptive]
        if rows:
            out += [f"### {FAMILY_TITLES[family]}", "", FAMILY_NOTES[family], ""]
            out += [line(x) for x in rows]
            out.append("")
    if descriptive:
        out += ["### Random split: descriptive only", "", RANDOM_SPLIT_SENTENCE, ""]
        out += [line(x) for x in descriptive]
        out.append("")
    if not report["contrasts"]:
        out += ["No contrasts.", ""]

    out += ["## Feature associations (descriptive)", "", FEATURE_SENTENCE, ""]
    for f in report["features"]:
        d = datasets[f["dataset"]]
        sign = "stable" if f["sign_stable"] else "not stable"
        out.append(f"- {_code(f['dataset'])} ({d['partition']}) {_code(f['column'])} from "
                   f"{', '.join(f['blocks'])}: selected in {f['selection_freq']:.0%} of fits, sign {sign}; "
                   f"claim level {_feature_claim(f)}.")
    if not report["features"]:
        out.append("No feature records.")
    out.append("")

    out += ["## Data-quality audits", "", "Missingness by block:", ""]
    for r in m["missingness"]:
        if r["indicator_auroc"] is None:
            auroc = "undefined"
        else:
            auroc = f"{r['indicator_auroc']:.3f} (oriented {r['indicator_auroc_oriented']:.3f})"
        flag = " **flagged**" if r["flag"] else ""
        rate_p = "not computed" if r["rate_diff_p"] is None else _p(r["rate_diff_p"])
        out.append(f"- {_code(r['dataset'])} / {_code(r['block'])}: missing in {r['rate_R']:.1%} of R and "
                   f"{r['rate_S']:.1%} of S (rate-difference p {rate_p}); missingness-indicator AUROC {auroc}{flag}")
    if not m["missingness"]:
        out.append("- none recorded")
    out += ["", "Group-proxy audit:", ""]
    for r in m["group_proxy_audit"]:
        value = "undefined" if r["value"] is None else f"{r['value']:g}"
        flag = " **flagged**" if r["flag"] else ""
        out.append(f"- {_code(r['dataset'])} / {_code(r['block'])}: {_md(r['metric'])} {value}{flag}")
    if not m["group_proxy_audit"]:
        out.append("- none recorded")
    out.append("")

    comp = m["compute"]
    out += ["## Compute", "",
            f"- Planned fit units: {comp['planned_fit_units']}; fits recorded in runs: "
            f"{sum(u['fit_count'] for u in m['runs'])}",
            f"- Actual seconds: {comp['actual_seconds']:g}; run seconds: {sum(u['seconds'] for u in m['runs']):g}",
            f"- Tiers completed: {', '.join(_md(t) for t in comp['tiers_completed']) or 'none'}",
            "", RULES_HEADING, "", RULES]
    return "\n".join(out)


def _load(path):
    try:
        return _loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError, RecursionError) as e:
        raise schema.SchemaError([f"{path}: {e}"]) from None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for name in ("manifest", "contrasts", "pooled", "features", "planned", "planned_pooled"):
        parser.add_argument(name, type=Path)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--not-produced", type=Path)
    parser.add_argument("--predictions-root", type=Path)
    args = parser.parse_args(argv)
    try:
        parts = [_load(getattr(args, name))
                 for name in ("manifest", "contrasts", "pooled", "features", "planned", "planned_pooled")]
        not_produced = _load(args.not_produced) if args.not_produced else None
        report = build_report(*parts, not_produced, prediction_root=args.predictions_root)
        write_report(report, args.out_dir)
    except schema.SchemaError as e:
        print("\n".join(e.errors), file=sys.stderr)
        return 1
    except (OSError, UnicodeError) as e:
        print(f"{args.out_dir}: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
