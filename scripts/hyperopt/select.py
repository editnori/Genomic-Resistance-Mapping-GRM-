"""Hyperparameter selection, budget accounting and stability over trial
histories.

A trial history is the JSONL file the search (Task B) writes for one outer
training fold: one line per trial, each scored on the same inner folds. This
module checks the protocol against its own protocol_hash, checks the history
against the protocol and Task B's record rules, refuses one that does not
match, and then, for each learner family:

- with Task B's protocol, refuses a history its search did not write. The
  fold plan the protocol pins is required and must hash to its
  fold_plan_hash; every ok trial must score exactly that plan's inner folds
  and every failed or unconverged trial a subset of them. The records must
  be the protocol's planned trials (planned_trials, search.py's
  deterministic plan), every one, each with its planned config, in run
  order, so an off-plan config, a fold subset or a truncated history is
  refused;
- excludes the failed and unconverged trials, recording each one, its status
  and its reason, and keeping whatever fold scores they reported;
- with a flat test-only protocol and no fold plan, excludes an ok trial
  scored on a different fold set from the history's inner folds
  (fold_set_mismatch), and flags the family. The history's inner folds are
  then the fold set the most ok trials report; a tie goes to the larger
  set, then to the lexicographically smallest fold list. A failed or
  unconverged trial reporting a fold outside that set is refused;
- scores every ok trial by its mean objective over the inner folds, and the
  standard error of those fold scores;
- enforces the protocol's constraints, each {metric, op: "<=" or ">=",
  bound}, inclusive. A constraint on complexity, seconds or fit_count applies
  to that trial field; any other applies to the mean of that metric over the
  trial's inner folds. An ok trial that fails one is excluded as infeasible,
  with the constraint and its value recorded;
- applies the one-standard-error rule to the feasible ok trials. The best
  trial has the best mean (ties go to the preference order below). The
  candidates are the feasible trials whose mean is within one SE of the best
  mean, where the SE is the best trial's.
  Among the candidates it prefers lower complexity, then the default config,
  then the lowest trial id. The candidates are also returned as the
  indistinguishable-candidate set.

The SE is descriptive. Inner folds share training data, so the fold scores are
not independent samples and the SE is not the basis of any test, interval or
p-value. SE_LABEL says so, and every output that carries an SE carries it.

Every mean, SD and SE is an inner-tuning score, not test performance
(SCORE_LABEL). Nothing here refits a selected config on the outer training
fold or scores it on the untouched outer test fold (OUTER_EVALUATION); that
remains a dependency outside this package.

It also counts each family's budget (planned against used trials, ok, failed
and unconverged trials, fits and seconds), flags unequal budgets, and, given
the selections from several histories (outer folds or repeats), reports how
often each config and each hyperparameter value was selected and how spread
the selected values are. Stability pools only runs that share everything in
their protocols but the fold plan, and never the same history twice.

This file is named select.py, which shadows the standard library's `select`
module for any program that puts this directory first on sys.path. Load it by
path (report.py and the tests do) or run Python with -I, which keeps a
script's own directory off sys.path.
"""
import copy
import hashlib
import io
import json
import math
from collections import Counter
from dataclasses import asdict, dataclass

STATUSES = ("ok", "failed", "unconverged")
DIRECTIONS = ("minimize", "maximize")
OPS = ("<=", ">=")
# Constraints on these name a trial field; any other names an inner metric.
TRIAL_FIELDS = ("complexity", "fit_count", "seconds")
# Task B's protocol schema version, and the key that marks a family default
# running the production code path's own internal search.
SCHEMA_VERSION = 1
INTERNAL_SEARCH = "internal_search"

SE_LABEL = (
    "SE is descriptive: the sample standard deviation of the best trial's "
    "inner-fold scores divided by the square root of the number of folds. "
    "Inner folds share training data, so the fold scores are not independent "
    "samples and this SE is not an independent-sample inference."
)

SCORE_LABEL = (
    "Every mean, SD and SE here is an inner-tuning score: the objective over "
    "the inner folds of one outer training fold, used only to choose a "
    "config. None is test performance. Choosing the best inner mean biases "
    "it upward, so a selected trial's inner-tuning mean overstates how its "
    "config scores on samples it was not tuned on."
)

OUTER_EVALUATION = (
    "Not in this package: refitting the selected config on the whole outer "
    "training fold and scoring it on the untouched outer test fold. That "
    "outer refit and evaluation is a remaining dependency; until a driver "
    "does it, no score here estimates performance."
)

RULE = (
    "one-standard-error rule: candidates are the ok trials that meet every "
    "constraint and whose mean objective is within one SE (the best "
    "trial's) of the best such mean; among "
    "them prefer lower complexity, then the default config, then the lowest "
    "trial id"
)

# Searched parameters production has no knob for. A protocol may still list
# one: before review F1, Task B's lasso searched lambda_min_ratio. Such a
# protocol verifies and selects as before; select_history labels it
# (`warnings`), so its selection isn't read as a config production can fit.
NO_PRODUCTION_KNOB = {
    ("lasso", "lambda_min_ratio"): (
        "production fixes the lasso path's smallest-penalty ratio at 1e-3 "
        "(lasso.rs path_problem) and has no knob for it, so a config that "
        "varies it is not one production can fit (review F1)"),
}

# Relative slack on the one-SE boundary, so a trial whose mean lies exactly on
# it in real arithmetic is not dropped by floating-point rounding.
BOUNDARY_RTOL = 1e-12

# The trial record schema: field name and a check of its value. The fields in
# TRIAL_FIELDS may also be null on a failed trial. `reason` is present exactly
# when status is not ok, and no other field is allowed.
FIELDS = {
    "protocol_hash": ("a string", lambda v: isinstance(v, str) and v != ""),
    "trial_id": ("a non-negative integer", lambda v: _is_int(v) and v >= 0),
    "family": ("a string", lambda v: isinstance(v, str) and v != ""),
    "config": ("an object with finite numbers",
               lambda v: isinstance(v, dict) and _canonical_ok(v)),
    "complexity": ("a finite non-negative number",
                   lambda v: _is_finite(v) and v >= 0),
    "is_default": ("a boolean", lambda v: isinstance(v, bool)),
    "seed": ("an integer", lambda v: _is_int(v)),
    "fold_plan_hash": ("a string", lambda v: isinstance(v, str) and v != ""),
    "status": ("one of " + ", ".join(STATUSES), lambda v: v in STATUSES),
    "inner_scores": ("a list", lambda v: isinstance(v, list)),
    "fit_count": ("a non-negative integer", lambda v: _is_int(v) and v >= 0),
    "seconds": ("a finite non-negative number",
                lambda v: _is_finite(v) and v >= 0),
}


class HistoryError(ValueError):
    """A history or protocol that does not match. `problems` lists every
    mismatch found, each as "code: detail"; `codes` is the set of codes."""

    def __init__(self, problems):
        self.problems = list(problems)
        super().__init__("; ".join(self.problems))

    @property
    def codes(self):
        return {p.split(":", 1)[0] for p in self.problems}


@dataclass(frozen=True)
class Protocol:
    """What a history must match: the objective metric and its direction, the
    constraints, the expected hashes, and the planned number of trials for
    each family. From Task B's protocol, also the fit seed every record must
    carry, each family's default config and search space, and the sampler
    seed that draws the planned configs; all three are None on a flat
    protocol. `search_hash` is the sha256 of the protocol without
    protocol_hash and fold_plan_hash: runs that differ only in their fold
    plan share it."""
    objective: str
    direction: str
    protocol_hash: str
    fold_plan_hash: str
    budgets: dict
    constraints: tuple = ()
    fit_seed: object = None
    families: object = None
    search_hash: str = ""
    sampler_seed: object = None

    @classmethod
    def from_dict(cls, d, trusted_flat=False):
        """Check a protocol mapping in either of two shapes.

        Task B's protocol (protocol.py): objective {name, direction,
        constraints}, seeds {fit, ...}, trials_per_family and families
        {family: {default, space}}, plus protocol_hash and fold_plan_hash.
        Every family gets trials_per_family. protocol_hash must equal
        content_hash(d), so a protocol edited after it was hashed is refused
        (protocol_hash_mismatch).

        The flat shape: objective (the metric name), direction, budgets
        {family: planned trials}, optional constraints, protocol_hash and
        fold_plan_hash. It has no content hash to check, so it is refused
        unless trusted_flat is true. Only tests set trusted_flat; report.py
        never does."""
        problems = []
        if not isinstance(d, dict):
            raise HistoryError(["bad_protocol: the protocol is not an object"])
        fit_seed = families = sampler_seed = None
        if isinstance(d.get("objective"), dict):
            obj = d["objective"]
            objective, direction = obj.get("name"), obj.get("direction")
            constraints = obj.get("constraints")
            families, n = d.get("families"), d.get("trials_per_family")
            budgets = None
            if not (isinstance(families, dict) and families):
                problems.append("bad_protocol: families must be a non-empty "
                                "object")
                families = None
            elif not (_is_int(n) and n >= 1):
                problems.append("bad_protocol: trials_per_family must be a "
                                "positive integer")
            else:
                budgets = {fam: n for fam in families}
            for fam, spec in (families or {}).items():
                problems += [f"bad_protocol: family {fam}: {p}"
                             for p in _family_problems(spec)]
            if not (_is_int(d.get("schema_version"))
                    and d["schema_version"] == SCHEMA_VERSION):
                problems.append(f"bad_protocol: schema_version must be the "
                                f"integer {SCHEMA_VERSION}")
            seeds = d.get("seeds")
            fit_seed = seeds.get("fit") if isinstance(seeds, dict) else None
            if not _is_int(fit_seed):
                problems.append("bad_protocol: seeds.fit must be an integer")
            sampler_seed = (seeds.get("sampler") if isinstance(seeds, dict)
                            else None)
            if not _is_int(sampler_seed):
                problems.append("bad_protocol: seeds.sampler must be an "
                                "integer")
            stored = d.get("protocol_hash")
            try:
                actual = content_hash(d)
            except (TypeError, ValueError) as e:
                problems.append(f"bad_protocol: the protocol is not canonical "
                                f"JSON ({e})")
            else:
                if isinstance(stored, str) and stored and stored != actual:
                    problems.append(
                        f"protocol_hash_mismatch: the protocol stores "
                        f"protocol_hash {stored}, but its content hashes to "
                        f"{actual}, so it was edited after it was hashed")
        elif not trusted_flat:
            raise HistoryError([
                "untrusted_protocol: this is not Task B's protocol (objective "
                "must be an object), and a flat protocol has no content hash "
                "to check, so it is accepted only in tests"])
        else:
            objective, direction = d.get("objective"), d.get("direction")
            constraints = d.get("constraints", [])
            budgets = d.get("budgets")
            if not (isinstance(budgets, dict) and budgets):
                problems.append("bad_protocol: budgets must be a non-empty "
                                "object of family to planned trials")
                budgets = None
        if not (isinstance(objective, str) and objective):
            problems.append("bad_protocol: the objective must be a non-empty "
                            "metric name")
        for key in ("protocol_hash", "fold_plan_hash"):
            if not (isinstance(d.get(key), str) and d.get(key)):
                problems.append(f"bad_protocol: {key} must be a non-empty string")
        if direction not in DIRECTIONS:
            problems.append("bad_protocol: direction must be one of "
                            + ", ".join(DIRECTIONS))
        for fam, n in (budgets or {}).items():
            if not (isinstance(fam, str) and fam and _is_int(n) and n >= 1):
                problems.append(f"bad_protocol: budget for {fam!r} must be "
                                "a positive integer")
        if not isinstance(constraints, list):
            problems.append("bad_protocol: constraints must be a list")
            constraints = []
        for c in constraints:
            if not (isinstance(c, dict) and set(c) == {"metric", "op", "bound"}
                    and isinstance(c["metric"], str) and c["metric"]
                    and c["op"] in OPS and _is_finite(c["bound"])):
                problems.append(f"bad_protocol: constraint {c!r} must be "
                                "{metric, op: <= or >=, bound: finite number}")
        if not problems:
            try:
                search = _sha256(_canonical(
                    {k: v for k, v in d.items()
                     if k not in ("protocol_hash", "fold_plan_hash")}))
            except (TypeError, ValueError) as e:
                problems.append(f"bad_protocol: the protocol is not canonical "
                                f"JSON ({e})")
        if problems:
            raise HistoryError(problems)
        return cls(objective, direction, d["protocol_hash"],
                   d["fold_plan_hash"], dict(sorted(budgets.items())),
                   tuple(dict(c) for c in constraints), fit_seed, families,
                   search, sampler_seed)

    def to_dict(self):
        d = asdict(self)
        d["constraints"] = list(d["constraints"])
        return d

    def metric_names(self):
        """The inner metrics every ok trial must report on every fold."""
        return {self.objective} | {c["metric"] for c in self.constraints
                                   if c["metric"] not in TRIAL_FIELDS}


def as_protocol(p, trusted_flat=False):
    if isinstance(p, Protocol):
        return p
    return Protocol.from_dict(p, trusted_flat=trusted_flat)


def content_hash(protocol):
    """Task B's protocol_hash: the sha256 hex digest of the canonical JSON
    (sort_keys, separators "," and ":", ensure_ascii, no NaN or infinity) of
    the protocol without its protocol_hash key, encoded as UTF-8."""
    return _sha256(_canonical({k: v for k, v in protocol.items()
                               if k != "protocol_hash"}))


def history_sha256(trials):
    """sha256 of the history as canonical JSON lines. Task B writes each
    record as one canonical line, so for a file its search wrote this is the
    sha256 of the file; for any other file the two differ, and report.py
    gives both."""
    return _sha256("".join(_canonical(t) + "\n" for t in trials))


def planned_trials(protocol):
    """Task B's deterministic plan (search.plan_trials) for a Task B
    protocol: every trial in run order, trial id then family name, as
    (family, trial_id, config). Trial 0 is the family default; a later
    trial draws each parameter as search.py documents, from seeds.sampler.
    A flat protocol has no plan, so None."""
    p = as_protocol(protocol)
    if p.families is None:
        return None
    n, = set(p.budgets.values())    # Task B gives every family one budget
    return [(fam, tid, copy.deepcopy(p.families[fam]["default"]) if tid == 0
             else {param: _sample_value(spec, _uniform(p.sampler_seed, fam,
                                                       tid, param))
                   for param, spec in sorted(p.families[fam]["space"].items())})
            for tid in range(n) for fam in sorted(p.families)]


def load_protocol(path):
    """Read a protocol JSON file. Refuses bad JSON, NaN or Infinity and
    repeated keys; from_dict checks the content."""
    return _load_json(path, "bad_protocol")


def load_fold_plan(path):
    """Read a fold plan JSON file, refusing what load_protocol refuses;
    validate_history checks it against the protocol."""
    return _load_json(path, "bad_fold_plan")


def _load_json(path, code):
    with open(path, encoding="utf-8") as f:
        text = f.read()
    try:
        return json.loads(text, object_pairs_hook=_no_dup_keys,
                          parse_constant=_no_constant)
    except (ValueError, RecursionError) as e:
        raise HistoryError([f"{code}: {path}: {e}"]) from None


def load_history(path):
    """Read a JSONL trial history file; see parse_history."""
    with open(path, "rb") as f:
        return parse_history(f.read())


def parse_history(data):
    """Parse the bytes of a JSONL trial history, decoded as UTF-8 with
    universal newlines, as open() reads a text file. A caller that also
    hashes the file hashes these same bytes. Refuses a last line without
    its newline (cut short by a crash), blank lines, lines that are not
    JSON, NaN or Infinity, and objects with a repeated key."""
    trials, problems = [], []
    with io.TextIOWrapper(io.BytesIO(data), encoding="utf-8") as f:
        text = f.read()
    if text and not text.endswith("\n"):
        problems.append("bad_json: the last line has no newline, so it may "
                        "be incomplete")
    for n, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            problems.append(f"bad_json: line {n} is blank")
            continue
        try:
            trials.append(json.loads(line, object_pairs_hook=_no_dup_keys,
                                     parse_constant=_no_constant))
        except (ValueError, RecursionError) as e:
            problems.append(f"bad_json: line {n}: {e}")
    if problems:
        raise HistoryError(problems)
    return trials


def validate_history(trials, protocol, fold_plan=None):
    """Check every record against the schema, Task B's record rules and the
    protocol.

    With Task B's protocol, also check that the protocol's search wrote the
    history: fold_plan, required, must hash to the protocol's
    fold_plan_hash (bad_fold_plan, missing_fold_plan); an ok trial must
    score exactly the plan's inner folds and any other trial a subset of
    them (off_plan_folds); each record must be a planned trial with its
    planned config (off_plan); every planned trial must be present
    (incomplete_history); and the records must be in run order
    (off_plan_order). A flat protocol has no plan; a fold_plan given with
    one is checked the same way.

    Returns the history's sorted inner fold ids: the fold plan's, or
    without one the fold set the most ok trials report, ties to the larger
    set and then the lexicographically smallest list (none if no trial is
    ok). Raises HistoryError listing every problem otherwise."""
    p = as_protocol(protocol)
    problems = []

    def bad(code, detail):
        problems.append(f"{code}: {detail}")

    if not isinstance(trials, list) or not trials:
        raise HistoryError(["empty_history: no trial records"])

    plan_folds = None
    if fold_plan is not None or p.families is not None:
        plan_folds = _plan_folds(fold_plan, p, bad)
    planned = planned_trials(p)
    plan_configs = (None if planned is None else
                    {(fam, tid): config for fam, tid, config in planned})
    order = []             # (family, trial_id) of each well-formed record
    first_seen = {}
    protocol_hashes, plan_hashes = {}, {}
    ok_folds = {}          # fold tuple -> records of ok trials
    other_folds = []       # (where, folds) of failed and unconverged trials
    ok_trials = []         # (where, trial, folds) of ok trials
    unbudgeted = {}
    for n, t in enumerate(trials, 1):
        where = f"record {n}"
        if not isinstance(t, dict):
            bad("not_an_object", f"{where} is not a JSON object")
            continue
        missing = [k for k in FIELDS if k not in t]
        if missing:
            bad("missing_field", f"{where} lacks {', '.join(missing)}")
            continue
        unknown = sorted(set(t) - set(FIELDS) - {"reason"})
        if unknown:
            bad("unknown_field", f"{where} has fields outside the record "
                f"schema: {', '.join(map(str, unknown))}")
            continue
        failed = t["status"] == "failed"
        wrong = [f"{k} must be {what}" for k, (what, ok) in FIELDS.items()
                 if not (ok(t[k]) or (failed and k in TRIAL_FIELDS
                                      and t[k] is None))]
        if t["status"] == "ok" and "reason" in t:
            wrong.append("reason must be absent on an ok trial")
        if (t["status"] in ("failed", "unconverged")
                and not (isinstance(t.get("reason"), str) and t["reason"])):
            wrong.append("reason must be a non-empty string on a "
                         f"{t['status']} trial")
        if wrong:
            bad("bad_field", f"{where}: {'; '.join(wrong)}")
            continue

        key = (t["family"], t["trial_id"])
        order.append(key)
        where = f"record {n} ({t['family']} trial {t['trial_id']})"
        if key in first_seen:
            bad("duplicate_trial_id",
                f"{where} repeats the trial id of record {first_seen[key]}")
        else:
            first_seen[key] = n
        protocol_hashes.setdefault(t["protocol_hash"], []).append(n)
        plan_hashes.setdefault(t["fold_plan_hash"], []).append(n)
        if t["family"] not in p.budgets:
            unbudgeted.setdefault(t["family"], []).append(n)
        if t["is_default"] != (t["trial_id"] == 0):
            bad("bad_default", f"{where}: is_default is {t['is_default']}, "
                "but trial 0, and only trial 0, is a family's default")
        if p.fit_seed is not None and t["seed"] != p.fit_seed:
            bad("bad_seed", f"{where}: seed is {t['seed']}, not the "
                f"protocol's fit seed {p.fit_seed}")
        if p.families is not None and t["family"] in p.families:
            spec = p.families[t["family"]]
            if t["is_default"]:
                if _canonical(t["config"]) != _canonical(spec["default"]):
                    bad("bad_config", f"{where}: trial 0's config is not the "
                        "family default")
            else:
                problem = _config_problem(t["config"], spec["space"])
                if problem:
                    bad("bad_config", f"{where}: config {problem}")
                elif key not in plan_configs:
                    bad("off_plan", f"{where} is not a planned trial: the "
                        f"plan has trials 0..{p.budgets[t['family']] - 1}")
                elif _canonical(t["config"]) != _canonical(plan_configs[key]):
                    bad("off_plan", f"{where}: config "
                        f"{_canonical(t['config'])} is not the planned "
                        f"config {_canonical(plan_configs[key])}")

        folds = _check_scores(t, p.metric_names(), where, bad)
        if folds is None:
            continue
        if plan_folds is not None and (
                folds != plan_folds if t["status"] == "ok"
                else not set(folds) <= set(plan_folds)):
            bad("off_plan_folds", f"{where} is {t['status']} on folds "
                f"{list(folds)}; the fold plan's inner folds are "
                f"{list(plan_folds)}")
            continue
        if t["status"] == "ok":
            ok_folds.setdefault(folds, []).append(n)
            ok_trials.append((where, t, folds))
        else:
            other_folds.append((where, folds))

    for value, recs in protocol_hashes.items():
        if value != p.protocol_hash:
            bad("protocol_hash_mismatch", f"expected {p.protocol_hash}, found "
                f"{value} in {_records(recs)}")
    if len(protocol_hashes) > 1:
        bad("mixed_protocol_hash", "the history holds "
            + ", ".join(f"{v} ({_records(r)})"
                        for v, r in sorted(protocol_hashes.items())))
    for value, recs in plan_hashes.items():
        if value != p.fold_plan_hash:
            bad("fold_plan_hash_mismatch", f"expected {p.fold_plan_hash}, "
                f"found {value} in {_records(recs)}")
    if len(plan_hashes) > 1:
        bad("mixed_fold_plan_hash", "the history holds "
            + ", ".join(f"{v} ({_records(r)})"
                        for v, r in sorted(plan_hashes.items())))
    for fam, recs in sorted(unbudgeted.items()):
        bad("unbudgeted_family",
            f"family {fam} has no budget in the protocol ({_records(recs)})")
    # Counting needs every record well-formed.
    if planned is not None and len(order) == len(trials):
        want = [(fam, tid) for fam, tid, _ in planned]
        absent = [key for key in want if key not in first_seen]
        if absent:
            bad("incomplete_history", f"the plan has {len(want)} trials, "
                f"and {len(absent)} are missing: "
                + ", ".join(f"{fam} {tid}" for fam, tid in absent[:5])
                + (f" and {len(absent) - 5} more" if len(absent) > 5 else ""))
        elif order != want and len(order) == len(want):
            bad("off_plan_order", "the records are not in the plan's run "
                "order (trial id, then family name)")

    # The history's inner folds are the fold plan's; without one, the fold
    # set the most ok trials report, a tie going to the larger set, then to
    # the lexicographically smallest fold list. An ok trial on any other set
    # is excluded later (fold_set_mismatch); a failed or unconverged trial
    # must stay inside it.
    folds = ()
    if plan_folds is not None:
        folds = plan_folds
    elif ok_folds:
        folds = min(ok_folds, key=lambda f: (-len(ok_folds[f]), -len(f), f))
        if len(folds) < 2:
            bad("too_few_folds", f"ok trials have {len(folds)} inner fold; "
                "the SE needs at least 2")
        for where, f in other_folds:
            if not set(f) <= set(folds):
                bad("fold_mismatch", f"{where} reports folds {list(f)}, "
                    f"outside the history's inner folds {list(folds)}")
    if not problems:
        # Finite inputs can still overflow once summed or squared.
        for where, t, f in ok_trials:
            if f == folds:
                _check_range(t, p, where, bad)
        for fam in sorted({t["family"] for t in trials}):
            if not _finite_sum(t["seconds"] for t in trials
                               if t["family"] == fam and t["seconds"] is not None):
                bad("out_of_range", f"family {fam}'s seconds do not sum to a "
                    "finite number")

    if problems:
        raise HistoryError(problems)
    return list(folds)


def fold_stats(scores):
    """Mean, sample SD (n-1) and SE (SD / sqrt(n)) of fold scores. The SE is
    descriptive; see SE_LABEL."""
    k = len(scores)
    mean = math.fsum(scores) / k
    sd = math.sqrt(math.fsum((x - mean) ** 2 for x in scores) / (k - 1))
    return mean, sd, sd / math.sqrt(k)


def constraint_check(trial, protocol):
    """An ok trial's value for each constrained metric, and the constraints it
    fails, each as text. A trial field is read as is; any other metric is
    averaged over the trial's inner folds. Bounds are inclusive."""
    p = as_protocol(protocol)
    values, failed = {}, []
    for c in p.constraints:
        m = c["metric"]
        if m not in values:
            values[m] = trial[m] if m in TRIAL_FIELDS else math.fsum(
                s["metrics"][m] for s in trial["inner_scores"]) / len(
                trial["inner_scores"])
        v = values[m]
        if not (v <= c["bound"] if c["op"] == "<=" else v >= c["bound"]):
            failed.append(f"constraint {m} {c['op']} {c['bound']} not met "
                          f"({m} is {v})")
    return values, failed


def preference_key(trial):
    """Order among indistinguishable candidates: lower complexity, then the
    default config, then the lowest trial id."""
    return (trial["complexity"], not trial["is_default"], trial["trial_id"])


def select_history(trials, protocol, fold_plan=None):
    """Validate one history and select a config for each budgeted family.
    Task B's protocol needs the fold plan it pins (see validate_history).

    Returns the objective, direction, hashes, the protocol's search_hash,
    planned budgets and constraints, history_sha256, inner folds, SE_LABEL,
    SCORE_LABEL, RULE, `warnings` (each family that searches a parameter in
    NO_PRODUCTION_KNOB, as {code, family, message}; the selection is
    unchanged), `families` (one result per budgeted family; see
    _select_family) and `budgets` (see budget_accounting)."""
    p = as_protocol(protocol)
    folds = validate_history(trials, p, fold_plan)
    by_family = {fam: [] for fam in p.budgets}
    for t in trials:
        by_family[t["family"]].append(t)
    warnings = [
        _flag("no_production_knob", fam, f"{fam} searches {param}: {why}")
        for (fam, param), why in sorted(NO_PRODUCTION_KNOB.items())
        if param in (p.families or {}).get(fam, {}).get("space", {})]
    return {
        "objective": p.objective,
        "direction": p.direction,
        "protocol_hash": p.protocol_hash,
        "fold_plan_hash": p.fold_plan_hash,
        "search_hash": p.search_hash,
        "planned": p.budgets,
        "constraints": list(p.constraints),
        "history_sha256": history_sha256(trials),
        "folds": folds,
        "se_label": SE_LABEL,
        "score_label": SCORE_LABEL,
        "rule": RULE,
        "warnings": warnings,
        "families": {fam: _select_family(fam, ts, folds, p)
                     for fam, ts in by_family.items()},
        "budgets": _budgets(trials, p, folds),
    }


def _select_family(family, trials, folds, p):
    """The one-SE selection for one family.

    `trials` has one row per trial, excluded ones included, ordered by trial
    id: status, reason (the trial's own, null on ok trials), config,
    complexity, `fold_scores` (the objective on every fold the trial
    reported, null where it reported none), `scores_used`, scores in fold
    order, mean, SD and SE, constraint values and failed constraints,
    whether it is excluded and why, and whether it is in the candidate set
    and selected. Scores, mean, SD, SE and constraint values are set only on
    an ok trial scored on the history's folds (`scores_used`); a failed,
    unconverged or fold_set_mismatch trial keeps its fold_scores as
    reported, used for nothing. An infeasible ok trial keeps its mean and SE
    but is excluded.
    `excluded` lists the excluded trials and why. `best` is the best-mean
    feasible trial, `threshold` the one-SE boundary, `candidate_set` the
    indistinguishable candidates in preference order, and `selected` the
    first of them. `status` is "selected", "no_feasible_trials" (ok trials,
    none meets the constraints), "no_ok_trials" or "no_trials"."""
    rows = []
    for t in sorted(trials, key=lambda t: t["trial_id"]):
        row = {
            "trial_id": t["trial_id"],
            "status": t["status"],
            "reason": t.get("reason"),
            "excluded": t["status"] != "ok",
            "excluded_because": [f"status is {t['status']}"]
            if t["status"] != "ok" else [],
            "is_default": t["is_default"],
            "complexity": t["complexity"],
            "seed": t["seed"],
            "config": t["config"],
            "fit_count": t["fit_count"],
            "seconds": t["seconds"],
            "fold_scores": [{"fold": s["fold"],
                             "score": s["metrics"].get(p.objective)}
                            for s in t["inner_scores"]],
            "scores_used": False,
            "scores": None, "mean": None, "sd": None, "se": None,
            "constraint_values": None,
            "in_candidate_set": False,
            "selected": False,
        }
        if t["status"] == "ok":
            by_fold = {s["fold"]: s["metrics"][p.objective]
                       for s in t["inner_scores"]}
            if list(by_fold) != folds:
                row["excluded"] = True
                row["excluded_because"] = [
                    f"fold_set_mismatch: scored on folds {list(by_fold)}, not "
                    f"the history's inner folds {folds}"]
            else:
                row["scores_used"] = True
                row["scores"] = [by_fold[f] for f in folds]
                row["mean"], row["sd"], row["se"] = fold_stats(row["scores"])
                row["constraint_values"], failed = constraint_check(t, p)
                row["excluded"] = bool(failed)
                row["excluded_because"] = failed
        rows.append(row)

    excluded = [{"trial_id": r["trial_id"], "status": r["status"],
                 "reason": r["reason"],
                 "excluded_because": r["excluded_because"]}
                for r in rows if r["excluded"]]
    n_ok = sum(r["status"] == "ok" for r in rows)
    result = {
        "family": family,
        "n_trials": len(rows),
        "n_ok": n_ok,
        "n_feasible": len(rows) - len(excluded),
        "excluded": excluded,
        "best": None,
        "threshold": None,
        "candidate_set": [],
        "selected": None,
        "status": ("no_trials" if not rows else
                   "no_feasible_trials" if n_ok else "no_ok_trials"),
        "se_label": SE_LABEL,
        "trials": rows,
    }
    ok = [r for r in rows if not r["excluded"]]
    if not ok:
        return result

    sign = 1 if p.direction == "minimize" else -1
    best = min(ok, key=lambda r: (sign * r["mean"], preference_key(r)))
    threshold = best["mean"] + sign * best["se"]
    slack = BOUNDARY_RTOL * max(1.0, abs(threshold))
    candidates = sorted(
        (r for r in ok if sign * (r["mean"] - threshold) <= slack),
        key=preference_key)
    for r in candidates:
        r["in_candidate_set"] = True
    chosen = candidates[0]
    chosen["selected"] = True
    result.update(
        status="selected",
        best={"trial_id": best["trial_id"], "mean": best["mean"],
              "se": best["se"]},
        threshold=threshold,
        candidate_set=[_brief(r) for r in candidates],
        selected=dict(_brief(chosen), config=chosen["config"],
                      is_best=chosen is best),
    )
    return result


def budget_accounting(trials, protocol, fold_plan=None):
    """Validate a history and count its budget per family."""
    p = as_protocol(protocol)
    return _budgets(trials, p, validate_history(trials, p, fold_plan))


def _budgets(trials, p, folds):
    """Per family: planned and used trials; ok, failed and unconverged counts;
    feasible (ok, on the history's folds and meeting every constraint),
    infeasible (on the folds but failing a constraint) and fold_set_mismatch
    ok counts; total fit_count and seconds over the trials that report them, and
    `cost_unreported`, the failed trials with a null fit_count or seconds.
    `flags` lists, each with a code: families under or over budget or with no
    trials, trial ids that are not 0..used-1, families with fold_set_mismatch
    trials, unequal planned, used, ok or (with constraints) feasible
    counts across families, and unequal total fit_count or seconds, since
    equal trial counts can hide unequal compute."""
    families, flags = {}, []
    for fam, planned in p.budgets.items():
        ts = [t for t in trials if t["family"] == fam]
        status = Counter(t["status"] for t in ts)
        on_folds = [t for t in ts if t["status"] == "ok"
                    and [s["fold"] for s in t["inner_scores"]] == folds]
        infeasible = sum(1 for t in on_folds if constraint_check(t, p)[1])
        families[fam] = {
            "planned": planned,
            "used": len(ts),
            "ok": status["ok"],
            "failed": status["failed"],
            "unconverged": status["unconverged"],
            "feasible": len(on_folds) - infeasible,
            "infeasible": infeasible,
            "fold_set_mismatch": status["ok"] - len(on_folds),
            "fit_count": sum(t["fit_count"] for t in ts
                             if t["fit_count"] is not None),
            "seconds": math.fsum(t["seconds"] for t in ts
                                 if t["seconds"] is not None),
            "cost_unreported": sum(1 for t in ts if t["fit_count"] is None
                                   or t["seconds"] is None),
        }
        if not ts:
            flags.append(_flag("no_trials", fam,
                               f"{fam} planned {planned} trials but has none"))
        elif len(ts) < planned:
            flags.append(_flag("under_budget", fam,
                               f"{fam} used {len(ts)} of {planned} planned trials"))
        elif len(ts) > planned:
            flags.append(_flag("over_budget", fam,
                               f"{fam} used {len(ts)} trials, {planned} planned"))
        ids = sorted(t["trial_id"] for t in ts)
        if ids and ids != list(range(len(ids))):
            flags.append(_flag("trial_id_gap", fam,
                               f"{fam} trial ids are not 0..{len(ids) - 1}"))
        off = families[fam]["fold_set_mismatch"]
        if off:
            flags.append(_flag("fold_set_mismatch", fam,
                               f"{fam} has {off} ok trial(s) scored on folds "
                               f"other than the history's {folds}, excluded"))
    for what in ("planned", "used", "ok") + (
            ("feasible",) if p.constraints else ()):
        values = {fam: b[what] for fam, b in families.items()}
        if len(set(values.values())) > 1:
            flags.append(_flag(f"unequal_{what}", None, f"{what} trials differ "
                               "across families: " + ", ".join(
                                   f"{f} {v}" for f, v in values.items())))
    # An internal-search default runs production's whole grid in one trial,
    # so equal trial counts can still spend unequal compute (review F6).
    for what, name in (("fit_count", "fits"), ("seconds", "seconds")):
        values = {fam: b[what] for fam, b in families.items()}
        if len(set(values.values())) > 1:
            flags.append(_flag(f"unequal_{what}", None, f"{name} differ "
                               "across families: " + ", ".join(
                                   f"{f} {v}" for f, v in values.items())))
    return {"families": families, "flags": flags}


def stability(selections, labels=None):
    """How often each config and hyperparameter value was selected across
    several histories' selections (select_history results), and how spread
    the selected values are.

    Per family: `selected_in` of `n_runs`, the runs with no selection, each
    run's pick, `configs` (distinct selected configs with count and fraction),
    `default_selected` (count), `params` (per hyperparameter: runs where it is
    present, value counts, distinct values, modal fraction, and for numeric
    values min, max, range, mean, population SD and median) and the same
    numeric spread of the selected complexity. Spreads describe the selected
    values only.

    Refuses runs whose protocols differ in anything but the fold plan
    (mixed_objective, mixed_budgets, mixed_constraints, or mixed_protocol
    for other fields such as seeds or search spaces), two runs of one
    experiment (duplicate_experiment, by protocol_hash), and the same history
    given twice (duplicate_history, by history_sha256)."""
    if not selections:
        raise HistoryError(["no_selections: stability needs at least one run"])
    labels = list(labels) if labels is not None else [
        f"run {i}" for i in range(len(selections))]
    problems = []
    if len(labels) != len(selections):
        problems.append("bad_labels: one label per selection is needed")
    if len(set(labels)) != len(labels):
        problems.append("bad_labels: run labels must be unique")
    goals = {(s["objective"], s["direction"]) for s in selections}
    if len(goals) > 1:
        problems.append("mixed_objective: the runs select on "
                        + ", ".join(f"{o} ({d})" for o, d in sorted(goals)))
    for key, code in (("planned", "mixed_budgets"),
                      ("constraints", "mixed_constraints")):
        seen = {_canonical(s[key]) for s in selections}
        if len(seen) > 1:
            problems.append(f"{code}: the runs' protocols have different "
                            f"{key}: " + "; ".join(sorted(seen)))
    if not problems and len({s["search_hash"] for s in selections}) > 1:
        problems.append("mixed_protocol: the runs' protocols differ in more "
                        "than their fold plans (seeds, search spaces or "
                        "other fields)")
    # Pooled runs differ only in their fold plan, so two runs sharing a
    # protocol_hash share the fold plan, seeds and search: one experiment,
    # however its history was copied, cut, split or rerun.
    for key, code, what in (
            ("protocol_hash", "duplicate_experiment",
             "the same experiment (protocol_hash {})"),
            ("history_sha256", "duplicate_history",
             "the same history (sha256 {})")):
        runs_by = {}
        for label, s in zip(labels, selections):
            runs_by.setdefault(s[key], []).append(label)
        for value, runs in sorted(runs_by.items()):
            if len(runs) > 1:
                problems.append(f"{code}: runs {', '.join(runs)} hold "
                                + what.format(value))
    if problems:
        raise HistoryError(problems)

    families = sorted({f for s in selections for f in s["families"]})
    out = {}
    for fam in families:
        picks, missing = [], []
        for label, s in zip(labels, selections):
            chosen = s["families"].get(fam, {}).get("selected")
            if chosen is None:
                missing.append(label)
            else:
                picks.append((label, chosen))
        k = len(picks)
        configs = {}
        for label, c in picks:
            entry = configs.setdefault(_canonical(c["config"]), {
                "config": c["config"], "count": 0, "runs": []})
            entry["count"] += 1
            entry["runs"].append(label)
        params = {}
        for name in sorted({n for _, c in picks for n in c["config"]}):
            values = [c["config"][name] for _, c in picks if name in c["config"]]
            counts = Counter(_canonical(v) for v in values)
            params[name] = {
                "present": len(values),
                "distinct": len(counts),
                "modal_fraction": max(counts.values()) / len(values),
                "value_counts": [{"value": json.loads(v), "count": n,
                                  "fraction": n / len(values)}
                                 for v, n in _by_count(counts)],
                "numeric": _spread(values),
            }
        out[fam] = {
            "n_runs": len(selections),
            "selected_in": k,
            "runs_without_selection": missing,
            "picks": [{"run": label, "trial_id": c["trial_id"],
                       "complexity": c["complexity"],
                       "is_default": c["is_default"], "config": c["config"]}
                      for label, c in picks],
            "configs": [dict(configs[key], fraction=configs[key]["count"] / k)
                        for key, _ in _by_count(
                            {key: e["count"] for key, e in configs.items()})],
            "default_selected": sum(c["is_default"] for _, c in picks),
            "params": params,
            "complexity": _spread([c["complexity"] for _, c in picks]),
        }
    return {"runs": labels, "families": out}


def _check_scores(t, metrics, where, bad):
    """Check a trial's inner scores: sorted by fold, no fold twice, and, on an
    ok trial, a finite score for each of `metrics` on every fold. Returns the
    fold tuple, or None after reporting a problem."""
    folds = []
    for i, s in enumerate(t["inner_scores"]):
        if not (isinstance(s, dict) and set(s) == {"fold", "metrics"}
                and _is_int(s["fold"]) and s["fold"] >= 0
                and isinstance(s["metrics"], dict)
                and all(isinstance(k, str) and k and _is_finite(v)
                        for k, v in s["metrics"].items())):
            bad("bad_scores", f"{where}: inner_scores[{i}] must be exactly "
                "{fold: non-negative int, metrics: {name: finite number}}")
            return None
        folds.append(s["fold"])
        if t["status"] == "ok":
            absent = sorted(m for m in metrics
                            if not _is_finite(s["metrics"].get(m)))
            if absent:
                bad("missing_metric", f"{where}: fold {s['fold']} has no "
                    f"finite score for {', '.join(absent)}")
                return None
    if len(set(folds)) != len(folds):
        bad("duplicate_fold", f"{where} scores a fold twice: {folds}")
        return None
    if folds != sorted(folds):
        bad("unsorted_folds", f"{where} lists folds {folds}, not sorted")
        return None
    if t["status"] == "ok" and not folds:
        bad("missing_metric", f"{where} is ok but has no inner scores")
        return None
    return tuple(folds)


def _check_range(t, p, where, bad):
    """Report an ok trial whose finite scores still overflow: a mean, SD,
    SE, one-SE boundary or constraint mean that is not a finite float."""
    scores = [s["metrics"][p.objective] for s in t["inner_scores"]]
    try:
        mean, sd, se = fold_stats(scores)
        values = [mean, sd, se, mean + se, mean - se]
        values += list(constraint_check(t, p)[0].values())
        ok = all(math.isfinite(v) for v in values)
    except (ArithmeticError, ValueError):
        ok = False
    if not ok:
        bad("out_of_range", f"{where}: its scores are too large for a finite "
            "mean, SD, SE or constraint value")


def _finite_sum(values):
    try:
        return math.isfinite(math.fsum(values))
    except (ArithmeticError, ValueError):
        return False


def _plan_folds(fold_plan, p, bad):
    """The fold plan's sorted inner fold ids as a tuple, read as Task B's
    inner_fold_ids reads them, once the plan hashes to the protocol's
    fold_plan_hash (Task B's fold_plan_hash). None after reporting a
    problem."""
    if fold_plan is None:
        bad("missing_fold_plan", f"the protocol pins fold plan "
            f"{p.fold_plan_hash}; without it no trial's folds can be checked")
        return None
    try:
        digest = _sha256(_canonical(fold_plan))
    except (TypeError, ValueError) as e:
        bad("bad_fold_plan", f"the fold plan is not plain JSON ({e})")
        return None
    if digest != p.fold_plan_hash:
        bad("bad_fold_plan", f"the fold plan hashes to {digest}, not the "
            f"protocol's fold_plan_hash {p.fold_plan_hash}")
        return None
    spec = fold_plan.get("inner_folds") if isinstance(fold_plan, dict) else None
    if _is_int(spec):
        folds = tuple(range(spec))
    elif isinstance(spec, (list, tuple)) and all(
            v is None or (_is_int(v) and v >= 0) for v in spec):
        folds = tuple(sorted({v for v in spec if v is not None}))
    else:
        bad("bad_fold_plan", "the fold plan's inner_folds must be an integer "
            "or a list of per-sample labels (non-negative integers or null)")
        return None
    if len(folds) < 2:
        bad("bad_fold_plan", f"the fold plan defines {len(folds)} inner "
            "folds; at least 2 are needed")
        return None
    return folds


def _uniform(seed, family, trial_id, param):
    """search.py's draw: a uniform in [0, 1) from the sha256 of the
    canonical [seed, family, trial_id, param]."""
    digest = hashlib.sha256(_canonical([seed, family, trial_id, param])
                            .encode()).digest()
    return (int.from_bytes(digest[:8], "big") >> 11) / 2**53


def _sample_value(spec, u):
    """search.sample_value: the value spec draws at u. The tests compare
    planned_trials with search.plan_trials, so the two cannot drift."""
    if spec["type"] == "choice":
        values = spec["values"]
        return values[min(math.floor(u * len(values)), len(values) - 1)]
    low, high = spec["low"], spec["high"]
    if spec["type"] == "int":
        if spec["log"]:
            value = math.exp(math.log(low) + u * (
                math.log(high + 1) - math.log(low)))
            value = math.floor(float(f"{value:.10g}"))
        else:
            value = low + math.floor(u * (high - low + 1))
    else:
        if spec["log"]:
            value = math.exp(math.log(low) + u * (math.log(high)
                                                  - math.log(low)))
        else:
            value = low + u * (high - low)
        value = float(f"{value:.10g}")
    return min(max(value, low), high)


def _family_problems(spec):
    """Why a protocol family is not {default, space: {param: spec}} with each
    spec an int or float range or a choice list and the default one of Task
    B's two kinds, or []. A fixed default is a config inside the space. An
    internal-search default is {"internal_search": true, param: [value,
    ...]} for every space parameter, each list non-empty, distinct and
    inside the space; it stands for the production code path's own search.
    No space parameter may be called internal_search."""
    if not (isinstance(spec, dict) and isinstance(spec.get("default"), dict)
            and isinstance(spec.get("space"), dict) and spec["space"]):
        return ["must hold a default config object and a non-empty space"]
    problems = []
    if INTERNAL_SEARCH in spec["space"]:
        problems.append(f"{INTERNAL_SEARCH} is reserved and cannot be a "
                        "space parameter")
    for param, s in spec["space"].items():
        kind = s.get("type") if isinstance(s, dict) else None
        if kind == "choice":
            ok = (isinstance(s.get("values"), list) and s["values"]
                  and all(isinstance(v, (str, bool)) or _is_finite(v)
                          for v in s["values"]))
        elif kind in ("int", "float"):
            typed = _is_int if kind == "int" else (
                lambda v: isinstance(v, float) and _is_finite(v))
            ok = (typed(s.get("low")) and typed(s.get("high"))
                  and s["low"] <= s["high"] and isinstance(s.get("log"), bool)
                  and (not s["log"] or s["low"] > 0))
        else:
            ok = False
        if not ok:
            problems.append(f"space {param} must be {{type: int or float, "
                            "low <= high, log: boolean, low > 0 on a log "
                            "scale} or {type: choice, values: [...]}")
    if problems:
        return problems
    problem = _default_problem(spec["default"], spec["space"])
    return [f"default: {problem}"] if problem else []


def _default_problem(default, space):
    if default.get(INTERNAL_SEARCH) is not True:
        return _config_problem(default, space)
    if set(default) != set(space) | {INTERNAL_SEARCH}:
        return (f"internal-search parameters {sorted(default)} differ from "
                f"{sorted(set(space) | {INTERNAL_SEARCH})}")
    for param, s in space.items():
        grid = default[param]
        if not (isinstance(grid, list) and grid):
            return f"{param} must be a non-empty list of values"
        if len({_canonical(v) for v in grid}) != len(grid):
            return f"{param} lists a value twice"
        outside = [v for v in grid if not _value_ok(v, s)]
        if outside:
            return f"{param} values {outside!r} are outside {_canonical(s)}"
    return None


def _config_problem(config, space):
    """Why a config is outside a (checked) search space, or None. It must
    hold exactly the space's parameters; an int parameter takes an integer,
    a float parameter a float, both within the inclusive bounds, and a
    choice parameter one of its values. A sampled config never carries
    internal_search, which is not a space parameter."""
    if set(config) != set(space):
        return (f"parameters {sorted(config)} differ from the space's "
                f"{sorted(space)}")
    for param, s in space.items():
        if not _value_ok(config[param], s):
            return f"{param}={config[param]!r} is outside {_canonical(s)}"
    return None


def _value_ok(v, s):
    if s["type"] == "choice":
        return _canonical(v) in {_canonical(x) for x in s["values"]}
    typed = _is_int(v) if s["type"] == "int" else isinstance(v, float)
    return typed and _is_finite(v) and s["low"] <= v <= s["high"]


def _brief(r):
    return {"trial_id": r["trial_id"], "mean": r["mean"], "se": r["se"],
            "complexity": r["complexity"], "is_default": r["is_default"]}


def _flag(code, family, message):
    return {"code": code, "family": family, "message": message}


def _spread(values):
    """Min, max, range, mean, population SD and median of the numbers among
    values. Anything else (a string, a boolean, the value list of an
    internal-search default) is skipped and counted. With no number, n is 0
    and the statistics are null, but skipped is still counted."""
    nums = [v for v in values if _is_num(v)]
    if not nums:
        return dict({"n": 0, "skipped": len(values)}, **dict.fromkeys(
            ("min", "max", "range", "mean", "sd_population", "median")))
    try:
        mean = math.fsum(nums) / len(nums)
        s = sorted(nums)
        mid = len(s) // 2
        out = {
            "n": len(nums), "skipped": len(values) - len(nums),
            "min": s[0], "max": s[-1], "range": s[-1] - s[0],
            "mean": mean,
            "sd_population": math.sqrt(math.fsum((x - mean) ** 2 for x in nums)
                                       / len(nums)),
            "median": s[mid] if len(s) % 2 else (s[mid - 1] + s[mid]) / 2,
        }
        ok = all(math.isfinite(v) for v in out.values())
    except (ArithmeticError, ValueError):
        ok = False
    if not ok:
        raise HistoryError(["out_of_range: the selected values are too "
                            "large for a finite spread"])
    return out


def _by_count(counts):
    """(key, count) pairs, most frequent first, ties by key."""
    return sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))


def _canonical(value):
    """Task B's canonical JSON."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True, allow_nan=False)


def _sha256(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _canonical_ok(value):
    try:
        _canonical(value)
    except ValueError:
        return False
    return True


def _records(recs):
    shown = ", ".join(str(r) for r in recs[:5])
    more = f" and {len(recs) - 5} more" if len(recs) > 5 else ""
    return f"record{'s' if len(recs) > 1 else ''} {shown}{more}"


def _no_constant(name):
    raise ValueError(f"{name} is not allowed")


def _no_dup_keys(pairs):
    keys = [k for k, _ in pairs]
    if len(set(keys)) != len(keys):
        raise ValueError(f"repeated key in {keys}")
    return dict(pairs)


def _is_int(v):
    return isinstance(v, int) and not isinstance(v, bool)


def _is_num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _is_finite(v):
    """A JSON number that is a finite float; an integer too large to convert
    is not."""
    if not _is_num(v):
        return False
    try:
        return math.isfinite(v)
    except OverflowError:
        return False
