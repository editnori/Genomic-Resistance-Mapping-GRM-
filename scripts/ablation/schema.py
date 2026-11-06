"""Validators for the ablation-report schema, version 1 (stdlib only).

Follows the GRM ablation data contract v1.10 (main orchestration ruling,
2026-10-30); where this file and the contract disagree, the contract wins.
Every validator raises SchemaError listing each problem as "<path>: <message>".
Records are closed: a missing or unknown key is an error, so a renamed field
fails loudly instead of being dropped. docs/ablation-report.md is the field
reference.
"""
import datetime
import hashlib
import json
import math
import re
import reprlib
from collections import defaultdict

SCHEMA = "ablation-report"
SCHEMA_VERSION = 1

CLAIM_LEVELS = ("predictive_increment", "association_stable",
                "mechanism_candidate", "causal_supported")
# Contrasts always carry this level (contract section 1).
CONTRAST_LEVEL = "predictive_increment"
# The highest level pipeline output may set on a feature.
FEATURE_PIPELINE_MAX = "association_stable"
# Levels only the manual-annotation path (annotate_claim) may set, on features.
MANUAL_LEVELS = ("mechanism_candidate", "causal_supported")

# Family -> its p-value adjustment; D is exploratory and carries no p-values.
FAMILY_ADJUST = {"A": "holm", "B": "holm", "C": "bh", "D": "none"}
FAMILIES = tuple(FAMILY_ADJUST)
CONFIRMATORY_FAMILIES = ("A", "B", "C")
# Contract section 10: A and B hold pooled records only, per-dataset contrasts are C or D.
CONTRAST_FAMILIES = ("C", "D")
POOLED_FAMILIES = ("A", "B", "D")
# Contract section 9: the confirmatory families use the primary endpoint only.
PRIMARY_CONTRAST_ENDPOINT = "log_loss"
PRIMARY_POOLED_ENDPOINT = "r2_ll"
# Contract section 12: ci_method is exactly one of these; "none" only on an undefined record.
CI_METHODS = ("cluster_jackknife_t", "none")
# Contract section 13: few_clusters is n_units < FEW_CLUSTERS; the exact sign-flip test cannot give
# a pooled p below 2 / 2^S, so with fewer than MIN_TESTABLE_SPECIES species p < 0.05 is unreachable.
FEW_CLUSTERS = 20
MIN_TESTABLE_SPECIES = 6
# Contract section 17: a defined pooled record needs MIN_POOLED_SPECIES contributing species, and its p
# can't fall below the exact sign-flip floor 2 / 2^S by more than FLOOR_TOLERANCE. A recomputed
# p_adj may differ from the record's by at most ADJUST_TOLERANCE.
MIN_POOLED_SPECIES = 2
FLOOR_TOLERANCE = 1e-12
ADJUST_TOLERANCE = 1e-12
# Contract section 13: J's undefined reasons under the cluster jackknife (an endpoint undefined on
# the full data keeps its own reason).
ONE_CLUSTER_REASON = "one cluster: a jackknife interval needs at least 2"
UNDEFINED_SHARE_REASON = "more than 10% of resamples undefined"
MAX_UNDEFINED_SHARE = 0.10
SOME_UNDEFINED_REASON = re.compile(r"(\d+) of (\d+) leave-one-cluster-out estimates undefined; the jackknife needs all")
PARTITIONS = ("development", "confirmatory")
ADJUSTMENTS = ("holm", "bh", "none")
# Contract section 15 withdraws equivalent: no floor held the false-equivalence rate at 5%.
DECISIONS = ("improves", "worse", "inconclusive")
# identical_cohort: an input-comparison arm whose genome set equals its baseline's, so it is not run.
RUN_STATUSES = ("ok", "unavailable", "identical_cohort")
# J's pairing rules (DECISIONS V-1). A run without `pairing` is paired by
# same_plan; common_outer_test_set runs record their own cohort plans.
PAIRINGS = ("same_plan", "k31_groups", "common_outer_test_set_balanced", "common_outer_test_set")
SAME_PLAN = "same_plan"
# Rules whose treatment scores a subset of the reference's samples and records
# its own plans; the unbalanced one is kept for older records.
SUBSET_PAIRINGS = ("common_outer_test_set_balanced", "common_outer_test_set")
ALPHA = 0.05
# H's BlockKind -> collapse priority (C=0 < H=1 < P=2 < N=3 < X=4 < K=5).
H_KIND_PRIORITY = {"CuratedExact": 0, "Hierarchical": 1, "Presence": 2, "Count": 3,
                   "Context": 4, "Anonymous": 5}
# The missingness flag is (oriented AUROC > 0.60) or (rate_diff_p < 0.01) (contract v1.1 section 8).
MISSINGNESS_AUROC_FLAG = 0.60
MISSINGNESS_P_FLAG = 0.01

ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9._+-]*")
BLOCK_ID = re.compile(r"[A-Za-z][A-Za-z0-9_-]*")
HEX64 = re.compile(r"[0-9a-f]{64}")
COMMIT = re.compile(r"[0-9a-f]{40}|[0-9a-f]{64}")
# A contrast key is "<family>:<learner>:<split>:<baseline>-><treatment>:<endpoint>" (contract
# section 12). Every part is non-empty and holds no whitespace, "/", ":" or backtick; the arms
# also hold no "->".
_PART = r"[^\s/:`]+"
_ARM = r"(?:(?!->)[^\s/:`])+"
_KEY = rf"([ABCD]):({_PART}):({_PART}):({_ARM})->({_ARM}):({_PART})"
KEY = re.compile(_KEY)
KEY_PARTS = ("family", "learner", "split", "baseline", "treatment", "endpoint")
# A contrast id is "<dataset id>/<contrast key>"; pooled records name the key.
CONTRAST_ID = re.compile(rf"({ID.pattern})/([ABCD]:{_PART}:{_PART}:{_ARM}->{_ARM}:{_PART})")
COLUMN = re.compile(r"[^`]+")
DATE = re.compile(r"\d{4}-\d{2}-\d{2}")
# Integers beyond this are not exact in JSON readers that use doubles.
MAX_INT = 2 ** 53
# Float slack for the oriented AUROC.
TOLERANCE = 1e-12


class SchemaError(ValueError):
    def __init__(self, errors):
        self.errors = list(errors)
        super().__init__("; ".join(self.errors))


class _Check:
    def __init__(self):
        self.errors = []

    def fail(self, path, message):
        self.errors.append(f"{path}: {message}")

    def done(self):
        if self.errors:
            raise SchemaError(self.errors)


# Canonical JSON and hashes (contract section 7).

def _string_keys(value, path):
    if isinstance(value, dict):
        for k, v in value.items():
            if not isinstance(k, str):
                # reprlib keeps a deep key's text short, and its repr from exhausting the stack
                raise SchemaError([f"{path}: key {reprlib.repr(k)} is not a string"])
            _string_keys(v, f"{path}.{k}")
    elif isinstance(value, (list, tuple)):
        for i, v in enumerate(value):
            _string_keys(v, f"{path}[{i}]")


def canonical(value, path="value"):
    _string_keys(value, path)
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)
    except (TypeError, ValueError) as e:
        raise SchemaError([f"{path}: not canonical JSON: {e}"]) from None


def sha256(value, path="value"):
    return hashlib.sha256(canonical(value, path).encode("utf-8")).hexdigest()


# Field specs: each is a function (check, path, value) that records problems.

def string(c, path, v):
    if not isinstance(v, str):
        c.fail(path, "must be a string")


def text(pattern=None):
    def spec(c, path, v):
        if not isinstance(v, str) or not v.strip():
            c.fail(path, "must be a non-empty string")
        elif pattern is not None and not pattern.fullmatch(v):
            c.fail(path, f"must match {pattern.pattern}")
    return spec


def nullable(inner):
    def spec(c, path, v):
        if v is not None:
            inner(c, path, v)
    return spec


def integer(lo=None, hi=None):
    def spec(c, path, v):
        if type(v) is not int:
            c.fail(path, "must be an integer")
        elif abs(v) > MAX_INT:
            c.fail(path, "must be within ±2^53")
        elif (lo is not None and v < lo) or (hi is not None and v > hi):
            c.fail(path, f"must be in [{lo}, {hi}]")
    return spec


def _finite(v):
    if type(v) is int:
        return abs(v) <= MAX_INT
    return type(v) is float and math.isfinite(v)


def number(lo=None, hi=None, positive=False):
    def spec(c, path, v):
        if not _finite(v):
            c.fail(path, "must be a finite number within ±2^53")
        elif (lo is not None and v < lo) or (hi is not None and v > hi):
            c.fail(path, f"must be in [{lo}, {hi}]")
        elif positive and not v > 0:
            c.fail(path, "must be > 0")
    return spec


def boolean(c, path, v):
    if type(v) is not bool:
        c.fail(path, "must be true or false")


def binary(c, path, v):
    if type(v) is not int or v not in (0, 1):
        c.fail(path, "must be the integer 0 or 1")


def enum(options):
    def spec(c, path, v):
        if not isinstance(v, (str, int)) or isinstance(v, bool) or v not in options:
            c.fail(path, f"must be one of {', '.join(map(str, options))}")
    return spec


def _parse_date(v):
    try:
        return isinstance(v, str) and DATE.fullmatch(v) and datetime.date.fromisoformat(v)
    except ValueError:
        return None


def date(c, path, v):
    if not _parse_date(v):
        c.fail(path, "must be an ISO date YYYY-MM-DD")


def past_date(c, path, v):
    day = _parse_date(v)
    if not day:
        c.fail(path, "must be an ISO date YYYY-MM-DD")
    elif day > datetime.date.today():
        c.fail(path, "must not be in the future")


def strings(nonempty=False, pattern=None):
    def spec(c, path, v):
        if not isinstance(v, list):
            return c.fail(path, "must be a list of strings")
        if nonempty and not v:
            c.fail(path, "must not be empty")
        item = text(pattern)
        for i, s in enumerate(v):
            item(c, f"{path}[{i}]", s)
        if all(isinstance(s, str) for s in v) and len(set(v)) != len(v):
            c.fail(path, "must not repeat entries")
    return spec


def interval(c, path, v):
    if not isinstance(v, list) or len(v) != 2 or not all(_finite(x) for x in v):
        c.fail(path, "must be [lo, hi] with finite numbers")
    elif v[0] > v[1]:
        c.fail(path, "lo must not exceed hi")


def mapping(value_spec):
    def spec(c, path, v):
        if not isinstance(v, dict):
            return c.fail(path, "must be an object")
        for k, x in v.items():
            if not isinstance(k, str) or not k.strip():
                c.fail(path, f"keys must be non-empty strings, not {reprlib.repr(k)}")
                continue  # str() of a deep key would exhaust the stack
            value_spec(c, f"{path}.{k}", x)
    return spec


def record(fields, optional=None, after=None):
    """A closed object: every field required, optional ones may be absent,
    and the cross-field `after` check runs only when the fields themselves pass."""
    optional = optional or {}

    def spec(c, path, v):
        if not isinstance(v, dict):
            return c.fail(path, "must be an object")
        before = len(c.errors)
        for k in sorted(set(fields) - set(v)):
            c.fail(f"{path}.{k}", "is missing")
        for k in sorted(set(v) - set(fields) - set(optional), key=str):
            c.fail(f"{path}.{k}", "is not a field of this record")
        for k, field in {**fields, **optional}.items():
            if k in v:
                field(c, f"{path}.{k}", v[k])
        if after is not None and len(c.errors) == before:
            after(c, path, v)
    return spec


def records(item, nonempty=False):
    def spec(c, path, v):
        if not isinstance(v, list):
            return c.fail(path, "must be a list")
        if nonempty and not v:
            c.fail(path, "must not be empty")
        for i, x in enumerate(v):
            item(c, f"{path}[{i}]", x)
    return spec


def _unique(c, path, rows, key, what):
    seen = set()
    for i, row in enumerate(rows):
        k = key(row)
        if k in seen:
            c.fail(f"{path}[{i}]", f"repeats {what} {k}")
        seen.add(k)


# Prediction rows (contract section 4; J's field names).

def _cluster(c, path, v):
    if isinstance(v, bool) or not ((isinstance(v, str) and v.strip()) or (type(v) is int and 0 <= v <= MAX_INT)):
        c.fail(path, "must be a non-empty string or an integer >= 0")


def _prediction_after(c, path, r):
    if r["called"] and r["score"] < 0.5:
        c.fail(f"{path}.called", "a called row needs score >= 0.5")


PREDICTION = record({
    "dataset": text(ID), "arm": text(), "learner": text(), "split": text(),
    "repeat": integer(0), "sample_id": text(), "y": binary, "fold": integer(0),
    "score": number(0, 1), "called": boolean, "cluster": _cluster, "fold_plan_sha256": text(HEX64),
}, after=_prediction_after)


def _prediction_key(r):
    return (r["dataset"], r["arm"], r["learner"], r["split"], r["repeat"], r["sample_id"])


def validate_prediction(rec):
    c = _Check()
    PREDICTION(c, "prediction", rec)
    c.done()
    return rec


def validate_predictions(rows, manifest=None):
    """Validate prediction rows as one set: one row per (dataset, arm, learner,
    split, repeat, sample_id); a sample keeps one label and one cluster; one
    fold per sample and one fold plan per (dataset, split, repeat), across arms
    and learners; with a manifest, every row belongs to an ok run and carries
    that run's plan for its repeat."""
    c = _Check()
    rows = list(rows)
    for i, r in enumerate(rows):
        PREDICTION(c, f"predictions[{i}]", r)
    c.done()
    _unique(c, "predictions", rows, _prediction_key, "key")
    seen = {}
    for i, r in enumerate(rows):
        checks = (
            ("y", (r["dataset"], r["sample_id"]), "label"),
            ("cluster", (r["dataset"], r["sample_id"]), "cluster"),
            ("fold", (r["dataset"], r["split"], r["repeat"], r["sample_id"]), "fold"),
            ("fold_plan_sha256", (r["dataset"], r["split"], r["repeat"]), "fold plan"),
        )
        for field, key, what in checks:
            first = seen.setdefault((field, key), r[field])
            if first != r[field]:
                c.fail(f"predictions[{i}].{field}", f"differs from another row's {what} for {key}")
    if manifest is not None:
        runs = {(u["dataset"], u["arm"], u["learner"], u["split"]): u
                for u in validate_manifest(manifest)["runs"] if u["status"] == "ok"}
        for i, r in enumerate(rows):
            run = runs.get((r["dataset"], r["arm"], r["learner"], r["split"]))
            if run is None:
                c.fail(f"predictions[{i}]", "belongs to no ok run in the manifest")
            elif r["repeat"] >= run["repeats"]:
                c.fail(f"predictions[{i}].repeat", f"run {run['id']} has {run['repeats']} repeats")
            elif r["fold_plan_sha256"] != run["fold_plan_sha256"][r["repeat"]]:
                c.fail(f"predictions[{i}].fold_plan_sha256", f"differs from run {run['id']}'s plan for this repeat")
        _check_paired_samples(c, rows, runs)
    c.done()
    return rows


def _check_paired_samples(c, rows, runs):
    """A run paired to a reference arm scores the reference's samples (same_plan,
    k31_groups) or a subset of them (SUBSET_PAIRINGS), per repeat."""
    samples = defaultdict(set)
    for r in rows:
        samples[(r["dataset"], r["arm"], r["learner"], r["split"], r["repeat"])].add(r["sample_id"])
    for (dataset, arm, learner, split), run in sorted(runs.items()):
        pairing = run.get("pairing")
        if pairing is None:
            continue
        for repeat in range(run["repeats"]):
            own = samples.get((dataset, arm, learner, split, repeat), set())
            ref = samples.get((dataset, pairing["reference_arm"], learner, split, repeat), set())
            fits = own <= ref if pairing["rule"] in SUBSET_PAIRINGS else own == ref
            if not fits:
                c.fail(f"manifest.runs.{run['id']}",
                       f"repeat {repeat} scores {len(own - ref)} samples outside, and misses {len(ref - own)} of, "
                       f"reference arm {pairing['reference_arm']}'s samples; {pairing['rule']} needs "
                       + ("a subset of them" if pairing["rule"] in SUBSET_PAIRINGS else "the same samples"))


# Run manifest.

def _block_after(c, path, b):
    if b["unavailable_reason"] is None:
        if not b["version"].strip():
            c.fail(f"{path}.version", "must be non-empty for an available block")
        if not HEX64.fullmatch(b["source_sha256"]):
            c.fail(f"{path}.source_sha256", "must be 64 lowercase hex characters for an available block")
    elif b["source_sha256"] and not HEX64.fullmatch(b["source_sha256"]):
        c.fail(f"{path}.source_sha256", "must be empty or 64 lowercase hex characters")
    if b["priority"] != H_KIND_PRIORITY[b["kind"]]:
        c.fail(f"{path}.priority", f"kind {b['kind']} has priority {H_KIND_PRIORITY[b['kind']]}")


def _run_after(c, path, u):
    plans = u["fold_plan_sha256"]
    if plans is not None and len(plans) != u["repeats"]:
        c.fail(f"{path}.fold_plan_sha256", f"needs one plan per repeat ({u['repeats']})")
    if u["status"] == "ok":
        if u["reason"] is not None:
            c.fail(f"{path}.reason", "must be null for an ok run")
        for k in ("predictions", "fold_plan_sha256"):
            if u[k] is None:
                c.fail(f"{path}.{k}", "is required for an ok run")
        if u["fit_count"] < 1:
            c.fail(f"{path}.fit_count", "an ok run has at least one fit")
    elif u["reason"] is None:
        c.fail(f"{path}.reason", "an unavailable run must give its reason")
    pairing = u.get("pairing")
    if pairing is not None:
        own = pairing["arm_plan_sha256"]
        if (own is not None) != (pairing["rule"] in SUBSET_PAIRINGS) and u["status"] == "ok":
            c.fail(f"{path}.pairing.arm_plan_sha256", f"is set exactly when the rule is one of {SUBSET_PAIRINGS}")
        elif own is not None and len(own) != u["repeats"]:
            c.fail(f"{path}.pairing.arm_plan_sha256", f"needs one plan per repeat ({u['repeats']})")
        if pairing["reference_arm"] == u["arm"]:
            c.fail(f"{path}.pairing.reference_arm", "names the run's own arm")


def _missingness_after(c, path, r):
    raw, oriented, p = r["indicator_auroc"], r["indicator_auroc_oriented"], r["rate_diff_p"]
    if (raw is None) != (oriented is None):
        return c.fail(path, "indicator_auroc and indicator_auroc_oriented are both null or both numbers")
    if raw is not None and abs(oriented - max(raw, 1 - raw)) > TOLERANCE:
        return c.fail(f"{path}.indicator_auroc_oriented", "must equal max(indicator_auroc, 1 - indicator_auroc)")
    expected = ((oriented is not None and oriented > MISSINGNESS_AUROC_FLAG)
                or (p is not None and p < MISSINGNESS_P_FLAG))
    if r["flag"] != expected:
        c.fail(f"{path}.flag", f"must be {str(expected).lower()}: the flag is (indicator_auroc_oriented > "
                               f"{MISSINGNESS_AUROC_FLAG}) or (rate_diff_p < {MISSINGNESS_P_FLAG})")


BLOCK = record({
    "id": text(BLOCK_ID), "kind": enum(tuple(H_KIND_PRIORITY)), "version": string, "source_sha256": string,
    "priority": integer(0, 5), "screen_cap": nullable(integer(1)), "catalog_date": nullable(date),
    "circularity_risk": boolean, "unavailable_reason": nullable(text()),
}, after=_block_after)
DATASET = record({
    "id": text(ID), "species": text(), "partition": enum(PARTITIONS),
    "n_samples": integer(0), "group_counts": mapping(integer(0)),
})
MISSINGNESS = record({
    "dataset": text(ID), "block": text(BLOCK_ID), "rate_R": number(0, 1), "rate_S": number(0, 1),
    "rate_diff_p": nullable(number(0, 1)), "indicator_auroc": nullable(number(0, 1)),
    "indicator_auroc_oriented": nullable(number(0.5, 1)), "flag": boolean,
}, after=_missingness_after)
GROUP_PROXY = record({
    "dataset": text(ID), "block": text(BLOCK_ID), "metric": text(), "value": nullable(number()), "flag": boolean,
})
RUN = record({
    "id": text(ID), "dataset": text(ID), "arm": text(), "blocks": strings(pattern=BLOCK_ID),
    "learner": text(), "split": text(), "repeats": integer(1),
    "fold_plan_sha256": nullable(records(text(HEX64), nonempty=True)), "status": enum(RUN_STATUSES),
    "reason": nullable(text()), "metrics": mapping(nullable(number())),
    "fit_count": integer(0), "seconds": number(0), "predictions": nullable(text()),
}, optional={
    # An arm paired to a reference arm by another rule than identical plans
    # (DECISIONS V-1); absent means same_plan with any arm of the same plan.
    "pairing": record({
        "rule": enum(PAIRINGS), "reference_arm": text(),
        "arm_plan_sha256": nullable(records(text(HEX64), nonempty=True)),
    }),
}, after=_run_after)
COMPUTE = record({
    "planned_fit_units": integer(0), "actual_seconds": number(0), "tiers_completed": strings(),
})
# Splits whose contrasts are descriptive only (main's ruling after the calibration study, 5 Nov 2026):
# planned in family D, printed in their own block.
DESCRIPTIVE_SPLITS = ("random",)
# The split a sentinel arm is judged on.
SENTINEL_SPLIT = "lineage"


def _manifest_after(c, path, m):
    arms = {u["arm"] for u in m["runs"]}
    for i, arm in enumerate(m.get("sentinel_arms", [])):
        if arm not in arms:
            c.fail(f"{path}.sentinel_arms[{i}]", f"arm {arm} has no run")
    _unique(c, f"{path}.base_rates", m.get("base_rates", []), lambda b: b["arm"], "arm")
    for i, b in enumerate(m.get("base_rates", [])):
        if b["arm"] not in arms:
            c.fail(f"{path}.base_rates[{i}].arm", f"arm {b['arm']} has no run")
        if (b["rate"] is None) != (b["seeds"] is None) or (b["rate"] is None and b["ci95"] is not None):
            c.fail(f"{path}.base_rates[{i}]", "a rate needs its seeds; a rate not measured has no seeds or interval")


MANIFEST = record({
    "schema_version": enum((SCHEMA_VERSION,)),
    "protocol_sha256": text(HEX64), "source_commit": text(COMMIT), "binary_sha256": text(HEX64),
    "blocks": records(BLOCK), "datasets": records(DATASET, nonempty=True),
    "missingness": records(MISSINGNESS), "group_proxy_audit": records(GROUP_PROXY),
    "runs": records(RUN, nonempty=True), "compute": COMPUTE,
}, optional={
    # Arms exchangeable with their baseline (the input comparison's learner-seed twin): any decision on
    # one makes every decision of that learner on that species unreliable.
    "sentinel_arms": strings(),
    # The input comparison's spec said "synthetic": true, so it ran without the revision-12 gate a real run needs
    # (main's ruling on REVIEW-INPUT9 L1).
    "synthetic": boolean,
    # A measured no-signal decision rate printed beside each contrast of the arm (main's rulings on
    # REVIEW-INPUT2 and REVIEW-INPUT3 I-c: the rate of the arm's kind, or null when that kind has none).
    # The input comparison's species pool: the union of the run's arms' cohorts, which its lineages and the
    # fold keys of units outside the default depend on (main's ruling on REVIEW-INPUT4 M1).
    "species_pool": record({"genomes": integer(1), "arms": strings(nonempty=True), "structure_sha256": text(HEX64)}),
    "base_rates": records(record({
        "arm": text(), "rate": nullable(number(0, 1)), "ci95": nullable(interval), "seeds": nullable(integer(1)),
        "source": text(),
    }, optional={"kind": text()})),
}, after=_manifest_after)


def validate_manifest(m):
    c = _Check()
    MANIFEST(c, "manifest", m)
    c.done()
    _unique(c, "manifest.blocks", m["blocks"], lambda b: b["id"], "block id")
    _unique(c, "manifest.datasets", m["datasets"], lambda d: d["id"], "dataset id")
    _unique(c, "manifest.runs", m["runs"], lambda u: u["id"], "run id")
    _unique(c, "manifest.runs", m["runs"],
            lambda u: (u["dataset"], u["arm"], u["learner"], u["split"]), "(dataset, arm, learner, split)")
    _unique(c, "manifest.missingness", m["missingness"], lambda r: (r["dataset"], r["block"]), "(dataset, block)")
    _unique(c, "manifest.group_proxy_audit", m["group_proxy_audit"],
            lambda r: (r["dataset"], r["block"], r["metric"]), "(dataset, block, metric)")
    blocks = {b["id"]: b for b in m["blocks"]}
    datasets = {d["id"] for d in m["datasets"]}
    for part in ("missingness", "group_proxy_audit"):
        for i, r in enumerate(m[part]):
            if r["dataset"] not in datasets:
                c.fail(f"manifest.{part}[{i}].dataset", f"unknown dataset {r['dataset']}")
            if r["block"] not in blocks:
                c.fail(f"manifest.{part}[{i}].block", f"unknown block {r['block']}")
    arm_blocks, fold_plans = {}, {}
    for i, u in enumerate(m["runs"]):
        path = f"manifest.runs[{i}]"
        if u["dataset"] not in datasets:
            c.fail(f"{path}.dataset", f"unknown dataset {u['dataset']}")
        for b in u["blocks"]:
            if b not in blocks:
                c.fail(f"{path}.blocks", f"unknown block {b}")
            elif u["status"] == "ok" and blocks[b]["unavailable_reason"] is not None:
                c.fail(f"{path}.status", f"ok run uses unavailable block {b}")
        if arm_blocks.setdefault(u["arm"], u["blocks"]) != u["blocks"]:
            c.fail(f"{path}.blocks", f"arm {u['arm']} has different blocks in another run")
        for repeat, plan in enumerate(u["fold_plan_sha256"] or []):
            if fold_plans.setdefault((u["dataset"], u["split"], repeat), plan) != plan:
                c.fail(f"{path}.fold_plan_sha256[{repeat}]",
                       "differs from another run's plan for the same dataset, split and repeat")
    keyed = {(u["dataset"], u["arm"], u["learner"], u["split"]): u for u in m["runs"]}
    for i, u in enumerate(m["runs"]):
        pairing = u.get("pairing")
        if pairing is None:
            continue
        ref = keyed.get((u["dataset"], pairing["reference_arm"], u["learner"], u["split"]))
        if ref is None:
            c.fail(f"manifest.runs[{i}].pairing.reference_arm",
                   f"arm {pairing['reference_arm']} has no run for {u['dataset']}, {u['learner']}, {u['split']}")
        elif ref.get("pairing") is not None:
            c.fail(f"manifest.runs[{i}].pairing.reference_arm",
                   f"arm {pairing['reference_arm']} is itself paired to another arm")
    c.done()
    return m


# The decision rule (contract section 1), applied to contrasts and pooled records.

def expected_decision(ci, p_adj):
    """Section 1's rule without the equivalent branch (section 15); with no p (family D) it is
    always inconclusive."""
    if ci is not None and p_adj is not None and p_adj < ALPHA:
        if ci[0] > 0:
            return "improves"
        if ci[1] < 0:
            return "worse"
    return "inconclusive"


def p_floor(n_species):
    """Section 13: the smallest exact sign-flip p over S species, 2 / 2^S."""
    return math.ldexp(1.0, 1 - n_species)


def _order(pvalues):
    return sorted(range(len(pvalues)), key=lambda i: pvalues[i])  # stable, like J's mergesort


def holm(pvalues):
    """Holm step-down adjusted p-values, in input order (J's stats.holm in plain Python)."""
    m, out, running = len(pvalues), [0.0] * len(pvalues), 0.0
    for rank, i in enumerate(_order(pvalues)):
        running = max(running, (m - rank) * pvalues[i])
        out[i] = min(running, 1.0)
    return out


def bh(pvalues):
    """Benjamini-Hochberg adjusted p-values, in input order (J's stats.bh in plain Python). The
    running minimum starts from the largest p times 1, so no value exceeds 1 and J's clip is a no-op."""
    m, out, running = len(pvalues), [0.0] * len(pvalues), math.inf
    order = _order(pvalues)
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        running = min(running, pvalues[i] * (m / (rank + 1)))
        out[i] = running
    return out


ADJUSTERS = {"holm": holm, "bh": bh}


def _estimates(c, path, v, defined):
    """Null shapes: an undefined or all-excluded record has no estimate,
    intervals or tests; a defined one has all three, nested, and p-values
    unless it is family D. Then the decision must follow the rule."""
    values = ("estimate", "ci", "ci_tost", "p", "p_adj")
    before = len(c.errors)
    if v["adjust"] != FAMILY_ADJUST[v["family"]]:
        c.fail(f"{path}.adjust", f"family {v['family']} uses {FAMILY_ADJUST[v['family']]}")
    if not defined:
        for k in values:
            if v[k] is not None:
                c.fail(f"{path}.{k}", "must be null when there is no estimate")
    else:
        for k in ("estimate", "ci", "ci_tost"):
            if v[k] is None:
                c.fail(f"{path}.{k}", "is required when the estimate is defined")
        for k in ("p", "p_adj"):
            if v["family"] == "D" and v[k] is not None:
                c.fail(f"{path}.{k}", "must be null in family D")
            if v["family"] != "D" and v[k] is None:
                c.fail(f"{path}.{k}", f"is required in family {v['family']}")
        if v["p"] is not None and v["p_adj"] is not None and v["p_adj"] < v["p"]:
            c.fail(f"{path}.p_adj", "an adjusted p-value cannot be below p")
        if v["ci"] is not None and v["ci_tost"] is not None and not (
                v["ci"][0] <= v["ci_tost"][0] and v["ci_tost"][1] <= v["ci"][1]):
            c.fail(f"{path}.ci_tost", "must nest inside ci: ci[0] <= ci_tost[0] <= ci_tost[1] <= ci[1]")
    if len(c.errors) == before:
        expected = expected_decision(v["ci"], v["p_adj"])
        if v["decision"] != expected:
            c.fail(f"{path}.decision", f"is {v['decision']} but the contract rule gives {expected}")


def parse_key(key):
    """The parts of a contrast key, or None when it does not follow the grammar."""
    match = KEY.fullmatch(key) if isinstance(key, str) else None
    return dict(zip(KEY_PARTS, match.groups())) if match else None


def _jackknife_reason(c, path, v):
    """Contract section 13: any undefined leave-one-cluster-out replicate makes the contrast
    undefined, and J's three jackknife reasons must agree with n_units and undefined_resamples."""
    reason, u, k = v["undefined_reason"], v["undefined_resamples"], v["n_units"]
    if reason is None:
        if u:
            c.fail(f"{path}.undefined_resamples", "must be 0 on a defined record; any undefined replicate "
                                                  "makes the contrast undefined")
    elif reason == ONE_CLUSTER_REASON:
        if k != 1:
            c.fail(f"{path}.undefined_reason", "says one cluster but n_units is not 1")
    elif reason == UNDEFINED_SHARE_REASON:
        if not u > MAX_UNDEFINED_SHARE * k:
            c.fail(f"{path}.undefined_reason", "says more than 10% but undefined_resamples is at most 10% of n_units")
    elif match := SOME_UNDEFINED_REASON.fullmatch(reason):
        if (int(match.group(1)), int(match.group(2))) != (u, k):
            c.fail(f"{path}.undefined_reason", "must name undefined_resamples of n_units")
        elif not 0 < u <= MAX_UNDEFINED_SHARE * k:
            c.fail(f"{path}.undefined_reason", "is for 1 to 10% undefined replicates; above that use "
                                               f"'{UNDEFINED_SHARE_REASON}'")


def _contrast_after(c, path, v):
    before = len(c.errors)
    if v["dataset"] != contrast_dataset(v["id"]):
        c.fail(f"{path}.dataset", f"must equal the id prefix {contrast_dataset(v['id'])}")
    parts = parse_key(contrast_key(v["id"]))
    for k in KEY_PARTS:
        if parts[k] != v[k]:
            c.fail(f"{path}.{k}", f"disagrees with the id's key, which gives {parts[k]}")
    if v["family"] not in CONTRAST_FAMILIES:
        c.fail(f"{path}.family", "families A and B hold pooled records only; per-dataset contrasts are C or D")
    elif v["family"] == "C" and v["endpoint"] != PRIMARY_CONTRAST_ENDPOINT:
        c.fail(f"{path}.endpoint", f"family C uses the primary endpoint {PRIMARY_CONTRAST_ENDPOINT} only; "
                                   "secondary endpoints are family D")
    if v["claim_level"] != CONTRAST_LEVEL:
        c.fail(f"{path}.claim_level", f"contrasts always carry {CONTRAST_LEVEL}; manual levels apply to features")
    if v["undefined_reason"] is None and v["ci_method"] == "none":
        c.fail(f"{path}.ci_method", "a defined record names its interval method")
    if v["few_clusters"] != (v["n_units"] < FEW_CLUSTERS):
        c.fail(f"{path}.few_clusters", f"must be n_units < {FEW_CLUSTERS}")
    _jackknife_reason(c, path, v)
    if (v["auroc_delta"] is None) != (v["auroc_ci"] is None):
        c.fail(f"{path}.auroc_ci", "auroc_delta and auroc_ci are both null or both set")
    elif v["auroc_delta"] is not None:
        if v["endpoint"] != PRIMARY_CONTRAST_ENDPOINT:
            c.fail(f"{path}.auroc_delta", f"is carried only beside a {PRIMARY_CONTRAST_ENDPOINT} decision")
        elif not v["auroc_ci"][0] <= v["auroc_delta"] <= v["auroc_ci"][1]:
            c.fail(f"{path}.auroc_ci", "must contain auroc_delta")
    if len(c.errors) == before:
        _estimates(c, path, v, v["undefined_reason"] is None)


CONTRAST = record({
    "id": text(CONTRAST_ID), "dataset": text(ID), "family": enum(FAMILIES),
    "baseline": text(), "treatment": text(), "learner": text(), "split": text(), "endpoint": text(),
    "estimate": nullable(number()), "ci": nullable(interval), "ci_tost": nullable(interval),
    "ci_method": enum(CI_METHODS), "unit": text(), "n_units": integer(1), "few_clusters": boolean,
    "p": nullable(number(0, 1)), "p_adj": nullable(number(0, 1)), "adjust": enum(ADJUSTMENTS),
    "margin": number(positive=True), "decision": enum(DECISIONS), "claim_level": enum(CLAIM_LEVELS),
    "undefined_resamples": integer(0), "undefined_reason": nullable(text()), "interval_label": text(),
    "pairing": enum(PAIRINGS),
    # J's descriptive AUROC delta beside a log-loss decision; it decides nothing.
    "auroc_delta": nullable(number()), "auroc_ci": nullable(interval),
}, after=_contrast_after)


def _pooled_after(c, path, v):
    before = len(c.errors)
    if v["n_excluded"] != len(v["excluded"]):
        c.fail(f"{path}.n_excluded", "must equal the number of excluded entries")
    included = v["n_datasets"] - v["n_excluded"]
    if included < 0:
        c.fail(f"{path}.n_datasets", "is the full denominator and must count every excluded dataset")
    if v["n_units"] > max(included, 0) or (v["n_units"] == 0) != (included <= 0):
        c.fail(f"{path}.n_units", "species must come from the included datasets")
    _unique(c, f"{path}.excluded", v["excluded"], lambda e: e["dataset"], "dataset")
    parts = parse_key(v["contrast"])
    if parts["family"] != v["family"]:
        c.fail(f"{path}.family", f"must equal the family its key implies, {parts['family']}")
    if v["family"] not in POOLED_FAMILIES:
        c.fail(f"{path}.family", "family C holds per-dataset contrasts only; pooled records are A, B or D")
    elif v["family"] != "D" and parts["endpoint"] != PRIMARY_POOLED_ENDPOINT:
        c.fail(f"{path}.contrast", f"family {v['family']} pools the primary endpoint {PRIMARY_POOLED_ENDPOINT} only")
    if len(c.errors) == before and included > 0:
        if v["n_units"] < MIN_POOLED_SPECIES:
            c.fail(f"{path}.n_units", f"a defined pooled record needs at least {MIN_POOLED_SPECIES} contributing "
                                      "species (section 17); with fewer, exclude every dataset with its reason")
        elif v["p"] is not None and v["p"] < p_floor(v["n_units"]) - FLOOR_TOLERANCE:
            c.fail(f"{path}.p", f"is below the exact sign-flip floor 2/2^S = {p_floor(v['n_units']):.6g} for "
                                f"S = {v['n_units']} species (section 17)")
    if len(c.errors) == before:
        _estimates(c, path, v, included > 0)


EXCLUSION = record({"dataset": text(ID), "reason": text()})
POOLED = record({
    "contrast": text(KEY), "family": enum(FAMILIES), "unit": enum(("species",)), "n_units": integer(0),
    "n_datasets": integer(0), "n_excluded": integer(0), "excluded": records(EXCLUSION),
    "estimate": nullable(number()), "ci": nullable(interval), "ci_tost": nullable(interval),
    "p": nullable(number(0, 1)), "p_adj": nullable(number(0, 1)), "adjust": enum(ADJUSTMENTS),
    "margin": number(positive=True), "decision": enum(DECISIONS), "interval_label": text(),
}, after=_pooled_after)


# Features and claim levels.

ANNOTATION = record({
    "by": text(), "date": past_date, "level": enum(MANUAL_LEVELS), "previous_level": enum(CLAIM_LEVELS),
})


def _feature_after(c, path, v):
    level, refs, notes = v["claim_level"], v["evidence_refs"], v.get("annotations")
    rank = CLAIM_LEVELS.index
    if level == "causal_supported" and not refs:
        c.fail(f"{path}.claim_level", "causal_supported requires non-empty evidence_refs")
    if rank(level) >= rank("association_stable") and not v["sign_stable"]:
        c.fail(f"{path}.claim_level", f"{level} requires sign_stable")
    if notes is None:
        if rank(level) > rank(FEATURE_PIPELINE_MAX):
            c.fail(f"{path}.claim_level", f"pipeline output may set at most {FEATURE_PIPELINE_MAX} on a feature; "
                                          f"{level} needs annotate_claim")
        return
    if level not in MANUAL_LEVELS:
        c.fail(f"{path}.annotations", "only record a manual mechanism_candidate or causal_supported")
    if not refs:
        c.fail(f"{path}.evidence_refs", "a manual claim level requires non-empty evidence_refs")
    if rank(notes[0]["previous_level"]) > rank(FEATURE_PIPELINE_MAX):
        c.fail(f"{path}.annotations[0].previous_level", f"the first annotation follows a pipeline level, "
                                                         f"at most {FEATURE_PIPELINE_MAX}")
    for i in range(1, len(notes)):
        if notes[i]["previous_level"] != notes[i - 1]["level"]:
            c.fail(f"{path}.annotations[{i}].previous_level", "must equal the level the annotation before set")
        if notes[i]["date"] < notes[i - 1]["date"]:
            c.fail(f"{path}.annotations[{i}].date", "must not precede the annotation before")
    if notes[-1]["level"] != level:
        c.fail(f"{path}.claim_level", "must equal the last annotation's level")


FEATURE = record({
    "dataset": text(ID), "column": text(COLUMN), "blocks": strings(nonempty=True, pattern=BLOCK_ID),
    "selection_freq": number(0, 1), "sign_stable": boolean,
    "claim_level": enum(CLAIM_LEVELS), "evidence_refs": strings(),
}, optional={"annotations": records(ANNOTATION, nonempty=True)}, after=_feature_after)

SPECS = {"contrast": CONTRAST, "pooled": POOLED, "feature": FEATURE}


def _validate(kind, rec):
    c = _Check()
    SPECS[kind](c, kind, rec)
    c.done()
    return rec


def validate_contrast(rec):
    return _validate("contrast", rec)


def validate_pooled(rec):
    return _validate("pooled", rec)


def validate_feature(rec):
    return _validate("feature", rec)


def contrast_dataset(contrast_id):
    return CONTRAST_ID.fullmatch(contrast_id).group(1)


def contrast_key(contrast_id):
    return CONTRAST_ID.fullmatch(contrast_id).group(2)


def annotate_claim(rec, kind, level, evidence_refs, by, on):
    """The manual-annotation path, the only way to set mechanism_candidate or
    causal_supported. Features only: the contract fixes contrasts at
    predictive_increment. Returns a new record whose annotations chain records
    who set each level, when, and the level it replaced."""
    if kind != "feature":
        raise SchemaError([f"annotate_claim: only features take manual levels; contrasts always carry "
                           f"{CONTRAST_LEVEL}"])
    if not isinstance(evidence_refs, list) or not all(isinstance(r, str) for r in evidence_refs):
        raise SchemaError(["annotate_claim: evidence_refs must be a list of strings"])
    if level not in MANUAL_LEVELS:
        raise SchemaError([f"annotate_claim: level must be one of {', '.join(MANUAL_LEVELS)}"])
    _validate(kind, rec)
    out = json.loads(canonical(rec))
    out["evidence_refs"] = list(dict.fromkeys([*rec["evidence_refs"], *evidence_refs]))
    out["claim_level"] = level
    out["annotations"] = [*out.get("annotations", []),
                          {"by": by, "date": on, "level": level, "previous_level": rec["claim_level"]}]
    return _validate(kind, out)
