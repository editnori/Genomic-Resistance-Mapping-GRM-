"""Hyperparameter protocol: frozen search spaces, budgets, seeds and objective.

A protocol is one JSON object, pinned by `protocol_hash`. The search runner
(`search.py`) refuses a protocol whose stored hash does not match its
content, so editing any field, including one search-space bound, needs a new
protocol and a new history.

Canonical JSON
--------------
`canonical_json(obj)` is `json.dumps(obj, sort_keys=True,
separators=(",", ":"), ensure_ascii=True, allow_nan=False)`, after refusing
any object key that is not a string. Floats are written with CPython's
shortest round-trip `repr`, and NaN or infinity is refused. `protocol_hash`
is the sha256 hex digest of the canonical JSON of the protocol without its
`protocol_hash` key. History lines use the same canonical form, so equal
records are equal bytes.

Fold plans
----------
The fold plan is the object the runner hands to the objective. It must be
plain JSON: objects with string keys, lists, strings, finite numbers,
booleans and null. Tuples are written as lists; numpy values, sets and
other types are refused, so convert them first. `fold_plan_hash(plan)` is
the sha256 hex digest of `canonical_json(plan)` encoded as ASCII, the same
bytes as Task D's `canaries.fold_plan_hash` for any plain-JSON plan.

The plan must be an object with an `inner_folds` key, either
- an integer k >= 2, meaning inner folds 0..k-1; or
- a list of per-sample inner fold labels, each a non-negative integer or
  null for a sample outside the inner split (an outer test sample). The
  inner folds are the distinct labels, and there must be at least 2.
Other keys are free and are covered by the hash.

Protocol fields (schema_version 1)
----------------------------------
- `schema_version`: the integer 1.
- `fold_plan_hash`: `fold_plan_hash` of the plan this protocol tunes on.
  The runner hashes the plan it is given and refuses a mismatch.
- `objective`: `{"name", "direction", "constraints"}`.
  `name` is an inner metric key, such as `balanced_accuracy`.
  `direction` is `"minimize"` or `"maximize"`. `constraints` is a list of
  `{"metric", "op", "bound"}` with `op` in `"<="`, `">="`. A constraint on
  `complexity`, `seconds` or `fit_count` applies to that trial field; any
  other metric applies to the mean of that metric over the trial's inner
  folds. The runner records constraints and checks their metrics are
  reported; selection enforces them.
- `seeds`: `{"sampler": int, "fit": int}`. `sampler` drives config draws.
  `fit` is passed to every objective call, so every trial in a protocol
  sees the same learner randomness and configs are compared on paired fits.
- `trials_per_family`: one budget shared by every family, default included.
- `families`: `{family: {"default": config, "space": {param: spec}}}`.
  Family and parameter names are non-empty strings, and no parameter may be
  called `internal_search`. A spec is one of
  `{"type": "int", "low": int, "high": int, "log": bool}`,
  `{"type": "float", "low": float, "high": float, "log": bool}` or
  `{"type": "choice", "values": [scalar, ...]}`.
  Bounds are inclusive. A sampled config holds exactly the space's
  parameters, each fixed to one value: an int parameter takes a JSON
  integer, a float parameter a JSON float, and a choice parameter one of
  its listed values.
- `protocol_hash`: as above.

Defaults
--------
Trial 0 of each family is the production default, exactly. A default is
one of two kinds:
- a fixed default: an ordinary config inside the space (rf, gbt, lasso).
  The objective fits it like any sampled config.
- an internal-search default: `{"internal_search": true, param: [value,
  ...], ...}`, holding every space parameter as a non-empty list of
  distinct in-space values (scm, cart). It means "run the production code
  path, whose own inner-CV search picks among these values". The future
  objective adapter must run that production path for such a config (for
  scm `scm::select` then `scm::fit`; for cart `tree::select_depth` then
  `tree::fit_tree`), not fit one of the listed values.
Sampled trials never carry `internal_search`; they fix every value, and the
adapter fits exactly that config with no inner search of its own.

History JSONL (one line per trial)
----------------------------------
Each line is the canonical JSON of one object, ended by a bare `\\n` (LF;
a `\\r\\n` ending is refused), with
these 12 fields:

- `protocol_hash` (str): the protocol that produced the trial.
- `trial_id` (int): index within the family, from 0. `(family, trial_id)`
  is unique within a history.
- `family` (str): a key of `families`.
- `config` (object): the parameters passed to the objective; for trial 0
  the family default, which may be an internal-search default.
- `complexity` (number or null): objective-reported model size.
- `is_default` (bool): true exactly when `trial_id` is 0.
- `seed` (int): `seeds.fit`, the seed the objective received.
- `fold_plan_hash` (str): equal to the protocol's.
- `status` (str): `"ok"`, `"failed"` or `"unconverged"`.
- `inner_scores` (list): `[{"fold": int, "metrics": {name: number}}]`,
  sorted by fold, folds unique. An `ok` trial reports exactly the fold
  plan's inner folds (so at least 2), and every fold reports the objective
  metric and every constrained inner metric. A `failed` or `unconverged`
  trial reports a subset of the plan's inner folds, possibly none.
- `fit_count` (int or null): learner fits the objective ran.
- `seconds` (number or null): objective-reported compute time. The runner
  does not time trials, so a deterministic objective writes a byte-identical
  history.

plus a 13th field present exactly when `status` is not `"ok"`:

- `reason` (str, never empty): why the trial is not ok. A reason starting
  `objective raised` or `invalid objective result` was written by the
  runner; any other came from the objective.

`complexity`, `fit_count` and `seconds` are null only on a `failed` trial.
A line must be byte-equal to the canonical JSON of its record, so repeated
keys, whitespace and other hand edits are refused. Lines are ordered by
`trial_id`, then family name, so an interrupted run has spent the same
number of trials on every family, give or take one.
"""
import copy
import hashlib
import json
import math

SCHEMA_VERSION = 1
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
STATUSES = ("ok", "failed", "unconverged")
TRIAL_FIELDS = ("complexity", "fit_count", "seconds")
HISTORY_FIELDS = frozenset({
    "protocol_hash", "trial_id", "family", "config", "complexity", "is_default", "seed",
    "fold_plan_hash", "status", "inner_scores", "fit_count", "seconds",
})
OPTIONAL_HISTORY_FIELDS = frozenset({"reason"})
INTERNAL_SEARCH = "internal_search"

# Trial 0 is the production default. rf, gbt and lasso fit one fixed config:
# `ForestParams::default()`, `GbtParams::default()` and the lasso's 30-point
# path with its 1e-3 smallest-penalty ratio. scm and cart search internally
# (`scm::select` at `scm.rs:218-262`, `tree::select_depth` at
# `tree.rs:228-251`), so their defaults list production's grids: scm both
# kinds, `P_GRID` and 1..`scm_max_rules` (6) rules; cart depths
# 1..`cart_max_depth` (8) with `min_leaf` 1. `generic` exists for tests.
FAMILIES = {
    "scm": {
        "default": {INTERNAL_SEARCH: True, "kind": ["conjunction", "disjunction"],
                    "max_rules": [1, 2, 3, 4, 5, 6],
                    "p": [0.1, 0.316, 0.562, 1.0, 1.778, 3.162, 10.0]},
        "space": {
            "kind": {"type": "choice", "values": ["conjunction", "disjunction"]},
            "max_rules": {"type": "int", "low": 1, "high": 10, "log": False},
            "p": {"type": "float", "low": 0.1, "high": 10.0, "log": True},
        },
    },
    "cart": {
        "default": {INTERNAL_SEARCH: True, "max_depth": [1, 2, 3, 4, 5, 6, 7, 8],
                    "min_leaf": [1]},
        "space": {
            "max_depth": {"type": "int", "low": 1, "high": 30, "log": False},
            "min_leaf": {"type": "int", "low": 1, "high": 20, "log": True},
        },
    },
    "rf": {
        "default": {"max_depth": 30, "min_leaf": 1, "n_trees": 200, "subsample": 0.632},
        "space": {
            "max_depth": {"type": "int", "low": 2, "high": 30, "log": False},
            "min_leaf": {"type": "int", "low": 1, "high": 10, "log": True},
            "n_trees": {"type": "int", "low": 50, "high": 500, "log": True},
            "subsample": {"type": "float", "low": 0.3, "high": 1.0, "log": False},
        },
    },
    "gbt": {
        "default": {"colsample": 0.5, "lambda": 1.0, "learning_rate": 0.1, "max_depth": 4,
                    "min_child_weight": 1.0, "rounds": 400},
        "space": {
            "colsample": {"type": "float", "low": 0.1, "high": 1.0, "log": False},
            "lambda": {"type": "float", "low": 0.1, "high": 10.0, "log": True},
            "learning_rate": {"type": "float", "low": 0.01, "high": 0.3, "log": True},
            "max_depth": {"type": "int", "low": 1, "high": 8, "log": False},
            "min_child_weight": {"type": "float", "low": 0.1, "high": 10.0, "log": True},
            "rounds": {"type": "int", "low": 50, "high": 800, "log": True},
        },
    },
    "lasso": {
        "default": {"lambda_min_ratio": 0.001, "n_lambdas": 30},
        "space": {
            "lambda_min_ratio": {"type": "float", "low": 0.0001, "high": 0.1, "log": True},
            "n_lambdas": {"type": "int", "low": 10, "high": 100, "log": False},
        },
    },
    "generic": {
        "default": {"mode": "a", "n": 4, "x": 0.0},
        "space": {
            "mode": {"type": "choice", "values": ["a", "b", "c"]},
            "n": {"type": "int", "low": 1, "high": 8, "log": False},
            "x": {"type": "float", "low": -1.0, "high": 1.0, "log": False},
        },
    },
}


class ProtocolError(ValueError):
    """A protocol, fold plan or history that must not be used."""


def canonical_json(obj):
    _check_keys(obj)
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
                      allow_nan=False)


def _check_keys(obj):
    if isinstance(obj, dict):
        for key, value in obj.items():
            if not isinstance(key, str):
                raise TypeError(f"canonical JSON needs string keys, got {key!r}")
            _check_keys(value)
    elif isinstance(obj, (list, tuple)):
        for value in obj:
            _check_keys(value)


def compute_hash(protocol):
    """sha256 of the canonical protocol without its `protocol_hash` key."""
    body = {key: value for key, value in protocol.items() if key != "protocol_hash"}
    return hashlib.sha256(canonical_json(body).encode()).hexdigest()


def fold_plan_hash(plan):
    """sha256 of the plan's canonical JSON; refuses a plan that is not plain JSON."""
    try:
        text = canonical_json(plan)
    except (TypeError, ValueError) as err:
        raise ProtocolError(f"fold plan must be plain JSON: {err}") from err
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def inner_fold_ids(plan):
    """Sorted inner fold ids the plan defines; refuses fewer than 2."""
    _require(isinstance(plan, dict) and "inner_folds" in plan,
             "fold plan must be an object with an inner_folds key")
    spec = plan["inner_folds"]
    if _is_int(spec):
        folds = list(range(spec))
    else:
        _require(isinstance(spec, (list, tuple))
                 and all(label is None or (_is_int(label) and label >= 0) for label in spec),
                 "fold plan inner_folds must be an integer or a list of per-sample labels "
                 "(non-negative integers or null)")
        folds = sorted({label for label in spec if label is not None})
    _require(len(folds) >= 2, f"fold plan defines {len(folds)} inner folds; at least 2 needed")
    return folds


def build_protocol(fold_plan_hash, *, families=LEARNERS, trials_per_family=32,
                   objective_name="balanced_accuracy", direction="maximize", constraints=(),
                   sampler_seed=20261030, fit_seed=42):
    """Return a validated, hashed protocol for the named families."""
    protocol = {
        "schema_version": SCHEMA_VERSION,
        "fold_plan_hash": fold_plan_hash,
        "objective": {"name": objective_name, "direction": direction,
                      "constraints": [dict(c) for c in constraints]},
        "seeds": {"sampler": sampler_seed, "fit": fit_seed},
        "trials_per_family": trials_per_family,
        "families": {name: copy.deepcopy(FAMILIES[name]) for name in families},
    }
    validate_body(protocol)
    protocol["protocol_hash"] = compute_hash(protocol)
    return protocol


def verify_protocol(protocol):
    """Refuse a malformed protocol or one whose content no longer matches its hash."""
    validate_body(protocol)
    stored = protocol.get("protocol_hash")
    actual = compute_hash(protocol)
    if stored != actual:
        raise ProtocolError(f"protocol hash mismatch: stored {stored}, content {actual}")


def save_protocol(protocol, path):
    verify_protocol(protocol)
    with open(path, "w", encoding="ascii") as stream:
        stream.write(json.dumps(protocol, sort_keys=True, indent=2, allow_nan=False) + "\n")


def load_protocol(path):
    try:
        with open(path, encoding="ascii") as stream:
            protocol = json.load(stream, parse_constant=_refuse_constant,
                                 object_pairs_hook=_unique_keys)
    except UnicodeDecodeError as err:
        raise ProtocolError(f"{path}: not ASCII JSON: {err}") from err
    verify_protocol(protocol)
    return protocol


def _refuse_constant(name):
    raise ProtocolError(f"non-finite number {name} in JSON")


def _unique_keys(pairs):
    keys = [key for key, _ in pairs]
    repeated = sorted({key for key in keys if keys.count(key) > 1})
    if repeated:
        raise ProtocolError(f"repeated JSON keys {repeated}")
    return dict(pairs)


def _require(condition, message):
    if not condition:
        raise ProtocolError(message)


def _is_int(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _is_number(value):
    """A finite JSON number; an int too large for a float is not one."""
    if not (_is_int(value) or isinstance(value, float)):
        return False
    try:
        return math.isfinite(value)
    except OverflowError:
        return False


def _is_name(value):
    return isinstance(value, str) and value != ""


def validate_body(protocol):
    _require(isinstance(protocol, dict), "protocol must be an object")
    keys = {"schema_version", "fold_plan_hash", "objective", "seeds", "trials_per_family",
            "families"}
    extra = set(protocol) - keys - {"protocol_hash"}
    _require(keys <= set(protocol) and not extra,
             f"protocol fields must be {sorted(keys)} plus protocol_hash")
    _require(_is_int(protocol["schema_version"]) and protocol["schema_version"] == SCHEMA_VERSION,
             f"unsupported schema_version {protocol['schema_version']!r}")
    _require(_is_name(protocol["fold_plan_hash"]), "fold_plan_hash must be a non-empty string")
    _require(_is_int(protocol["trials_per_family"]) and protocol["trials_per_family"] >= 1,
             "trials_per_family must be an integer >= 1")
    seeds = protocol["seeds"]
    _require(isinstance(seeds, dict) and set(seeds) == {"sampler", "fit"}
             and all(_is_int(seeds[k]) and seeds[k] >= 0 for k in seeds),
             "seeds must be {sampler, fit} non-negative integers")
    _validate_objective(protocol["objective"])
    families = protocol["families"]
    _require(isinstance(families, dict) and families, "families must be a non-empty object")
    for name, family in families.items():
        _require(_is_name(name), f"family name {name!r} must be a non-empty string")
        _require(isinstance(family, dict) and set(family) == {"default", "space"},
                 f"family {name} must have exactly default and space")
        space = family["space"]
        _require(isinstance(space, dict) and space, f"family {name} has an empty space")
        for param, spec in space.items():
            _require(_is_name(param) and param != INTERNAL_SEARCH,
                     f"family {name} parameter {param!r} must be a non-empty string other "
                     f"than {INTERNAL_SEARCH}")
            _validate_spec(f"{name}.{param}", spec)
        problem = default_problem(family["default"], space)
        _require(problem is None, f"family {name} default: {problem}")


def _validate_objective(objective):
    _require(isinstance(objective, dict)
             and set(objective) == {"name", "direction", "constraints"},
             "objective must have exactly name, direction and constraints")
    _require(_is_name(objective["name"]), "objective name must be a non-empty string")
    _require(objective["direction"] in ("minimize", "maximize"),
             "objective direction must be minimize or maximize")
    _require(isinstance(objective["constraints"], list), "constraints must be a list")
    for c in objective["constraints"]:
        _require(isinstance(c, dict) and set(c) == {"metric", "op", "bound"}
                 and _is_name(c["metric"])
                 and c["op"] in ("<=", ">=") and _is_number(c["bound"]),
                 f"constraint {c!r} must be {{metric, op: <= or >=, bound: number}}")


def _validate_spec(where, spec):
    _require(isinstance(spec, dict), f"{where}: spec must be an object")
    kind = spec.get("type")
    if kind == "choice":
        values = spec.get("values")
        _require(set(spec) == {"type", "values"} and isinstance(values, list) and values,
                 f"{where}: choice needs a non-empty values list")
        _require(all(isinstance(v, (str, bool)) or _is_number(v) for v in values),
                 f"{where}: choice values must be strings, booleans or finite numbers")
        _require(len({canonical_json(v) for v in values}) == len(values),
                 f"{where}: choice values repeat")
        return
    _require(kind in ("int", "float") and set(spec) == {"type", "low", "high", "log"},
             f"{where}: spec must be int, float or choice with its exact fields")
    check = _is_int if kind == "int" else (lambda v: isinstance(v, float) and math.isfinite(v))
    _require(check(spec["low"]) and check(spec["high"]) and spec["low"] <= spec["high"],
             f"{where}: {kind} bounds must be {kind}s with low <= high")
    _require(isinstance(spec["log"], bool), f"{where}: log must be a boolean")
    _require(not spec["log"] or spec["low"] > 0, f"{where}: a log scale needs low > 0")


def _value_ok(value, spec):
    if spec["type"] == "choice":
        return canonical_json(value) in {canonical_json(v) for v in spec["values"]}
    typed = _is_int(value) if spec["type"] == "int" else isinstance(value, float)
    return typed and _is_number(value) and spec["low"] <= value <= spec["high"]


def config_problem(config, space):
    """Return why a sampled `config` is outside `space`, or None when it lies inside."""
    if not isinstance(config, dict):
        return "config must be an object"
    if set(config) != set(space):
        return f"parameters {sorted(config)} differ from space {sorted(space)}"
    for param, spec in space.items():
        if not _value_ok(config[param], spec):
            return f"{param}={config[param]!r} is outside {canonical_json(spec)}"
    return None


def is_internal_search(config):
    """True for an internal-search default: run the production code path."""
    return isinstance(config, dict) and config.get(INTERNAL_SEARCH) is True


def default_problem(default, space):
    """Return why `default` is neither a fixed nor an internal-search default."""
    if not is_internal_search(default):
        return config_problem(default, space)
    if set(default) != set(space) | {INTERNAL_SEARCH}:
        return (f"internal-search parameters {sorted(default)} differ from "
                f"{sorted(set(space) | {INTERNAL_SEARCH})}")
    for param, spec in space.items():
        grid = default[param]
        if not (isinstance(grid, list) and grid):
            return f"{param} must list the values the internal search covers"
        if len({canonical_json(v) for v in grid}) != len(grid):
            return f"{param} grid repeats a value"
        outside = [v for v in grid if not _value_ok(v, spec)]
        if outside:
            return f"{param} grid values {outside!r} are outside {canonical_json(spec)}"
    return None


def validate_record(record, protocol=None, folds=None):
    """Refuse a history record that breaks the documented JSONL schema.

    With `protocol`, the record must also belong to it. With `folds`, the
    fold plan's inner fold ids, its inner scores must cover exactly those
    folds when ok, and a subset of them otherwise.
    """
    _require(isinstance(record, dict), "history record must be an object")
    fields = set(record)
    missing = HISTORY_FIELDS - fields
    extra = fields - HISTORY_FIELDS - OPTIONAL_HISTORY_FIELDS
    _require(not missing and not extra,
             f"history record fields: missing {sorted(missing)}, unknown {sorted(extra)}")
    status = record["status"]
    _require(status in STATUSES, f"unknown status {status!r}")
    _require(_is_int(record["trial_id"]) and record["trial_id"] >= 0,
             "trial_id must be a non-negative integer")
    _require(record["is_default"] is (record["trial_id"] == 0),
             "is_default must be true exactly for trial 0")
    _require(_is_int(record["seed"]), "seed must be an integer")
    for key in ("protocol_hash", "family", "fold_plan_hash"):
        _require(_is_name(record[key]), f"{key} must be a non-empty string")
    _require(isinstance(record["config"], dict), "config must be an object")
    _require(("reason" in record) is (status != "ok"),
             "reason must be present exactly when status is not ok")
    if status != "ok":
        _require(_is_name(record["reason"]), "reason must be a non-empty string")
    for key in TRIAL_FIELDS:
        value = record[key]
        valid = _is_int(value) if key == "fit_count" else _is_number(value)
        _require((value is None and status == "failed") or (valid and value >= 0),
                 f"{key} must be a non-negative number (null only when failed)")
    problem = inner_scores_problem(record["inner_scores"])
    _require(problem is None, f"inner_scores: {problem}")
    if status == "ok":
        _require(record["inner_scores"], "an ok trial needs inner scores")
    if folds is not None:
        problem = folds_problem(record, folds)
        _require(problem is None, f"inner_scores: {problem}")
    if protocol is not None:
        _check_against_protocol(record, protocol)


def folds_problem(record, folds):
    """Return why the record's folds do not fit the plan's inner `folds`, or None."""
    seen = [entry["fold"] for entry in record["inner_scores"]]
    if record["status"] == "ok":
        if len(seen) < 2:
            return f"an ok trial needs at least 2 inner folds, got {seen}"
        if seen != list(folds):
            return (f"an ok trial must report exactly the plan's inner folds "
                    f"{list(folds)}, got {seen}")
    elif not set(seen) <= set(folds):
        return f"folds {sorted(set(seen) - set(folds))} are not in the plan's inner folds"
    return None


def _check_against_protocol(record, protocol):
    family = protocol["families"].get(record["family"])
    _require(family is not None, f"family {record['family']!r} is not in the protocol")
    _require(record["protocol_hash"] == protocol["protocol_hash"],
             f"record protocol_hash {record['protocol_hash']} is not "
             f"{protocol['protocol_hash']}")
    _require(record["fold_plan_hash"] == protocol["fold_plan_hash"],
             "record fold_plan_hash differs from the protocol")
    _require(record["seed"] == protocol["seeds"]["fit"], "record seed is not seeds.fit")
    _require(record["trial_id"] < protocol["trials_per_family"],
             "trial_id exceeds trials_per_family")
    if record["is_default"]:
        _require(canonical_json(record["config"]) == canonical_json(family["default"]),
                 "trial 0 must be the default config")
    else:
        problem = config_problem(record["config"], family["space"])
        _require(problem is None, f"config: {problem}")
    if record["status"] == "ok":
        problem = missing_metrics(record["inner_scores"], required_metrics(protocol))
        _require(problem is None, problem)


def inner_scores_problem(scores):
    """Return why `scores` is not a valid inner_scores list, or None."""
    if not isinstance(scores, list):
        return "must be a list"
    folds = []
    for entry in scores:
        if not (isinstance(entry, dict) and set(entry) == {"fold", "metrics"}):
            return "each entry must be exactly {fold, metrics}"
        if not (_is_int(entry["fold"]) and entry["fold"] >= 0):
            return "fold must be a non-negative integer"
        metrics = entry["metrics"]
        if not (isinstance(metrics, dict) and all(_is_name(k) for k in metrics)):
            return f"fold {entry['fold']} metrics must be an object keyed by name"
        bad = sorted(k for k, v in metrics.items() if not _is_number(v))
        if bad:
            return f"fold {entry['fold']} metrics {bad} are not finite numbers"
        folds.append(entry["fold"])
    if folds != sorted(set(folds)):
        return "folds must be unique and sorted"
    return None


def required_metrics(protocol):
    """Inner metrics every fold of an ok trial must report."""
    objective = protocol["objective"]
    names = {objective["name"]}
    names.update(c["metric"] for c in objective["constraints"] if c["metric"] not in TRIAL_FIELDS)
    return names


def missing_metrics(scores, names):
    for entry in scores:
        absent = sorted(names - set(entry["metrics"]))
        if absent:
            return f"fold {entry['fold']} is missing metrics {absent}"
    return None


def history_line(record):
    return canonical_json(record) + "\n"


def read_history(path, protocol=None, folds=None):
    """Parse and validate a history JSONL; refuse mixed protocol hashes.

    Every line must be canonical and end in a bare `\\n`; the file is read
    without newline translation, so a `\\r\\n` ending is refused. A last
    line without a newline was cut short by a crash; it is refused, not
    resumed, and the operator must delete the partial line before
    rerunning. With `protocol`, every record must also belong to it; with
    `folds`, its inner scores must fit those inner fold ids (see
    `validate_record`).
    """
    try:
        with open(path, encoding="ascii", newline="") as stream:
            text = stream.read()
    except UnicodeDecodeError as err:
        raise ProtocolError(f"{path}: history is not ASCII: {err}") from err
    _require(not text or text.endswith("\n"),
             f"{path}: last line is incomplete, probably cut short by a crash; delete "
             "everything after the last newline, then rerun to resume")
    records = []
    for number, line in enumerate(text.split("\n")[:-1], 1):
        try:
            _require("\r" not in line,
                     "carriage return in line; history lines must end in a bare \\n")
            record = json.loads(line, parse_constant=_refuse_constant,
                                object_pairs_hook=_unique_keys)
            validate_record(record, folds=folds)
            _require(canonical_json(record) == line, "line is not canonical JSON")
        except (ValueError, TypeError) as err:
            raise ProtocolError(f"{path}:{number}: {err}") from err
        records.append(record)
    hashes = sorted({r["protocol_hash"] for r in records})
    _require(len(hashes) <= 1, f"{path}: history mixes protocol hashes {hashes}")
    keys = [(r["family"], r["trial_id"]) for r in records]
    _require(len(set(keys)) == len(keys), f"{path}: a (family, trial_id) repeats")
    if protocol is not None:
        for number, record in enumerate(records, 1):
            try:
                _check_against_protocol(record, protocol)
            except ProtocolError as err:
                raise ProtocolError(f"{path}:{number}: {err}") from err
    return records
