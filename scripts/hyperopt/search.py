"""Deterministic random search over a protocol, with an injected objective.

Sampling
--------
Trial 0 of every family is its default config, the production default (see
`protocol.py`). Each later draw is `u = (first 8 bytes of
sha256(canonical_json([seeds.sampler, family, trial_id, param]))
>> 11) / 2**53`, a uniform in [0, 1) that depends on no other parameter, so
adding a parameter or a family leaves the other draws unchanged.

- int: `low + floor(u * (high - low + 1))`; on a log scale
  `floor(exp(log(low) + u * (log(high + 1) - log(low))))`.
- float: `low + u * (high - low)`, or `exp` of the same on a log scale,
  then rounded to 10 significant digits so a last-bit difference in a
  platform's `exp` does not change the history bytes.
- choice: `values[floor(u * len(values))]`.

Every value is clamped into its bounds. Draws are not deduplicated: a small
discrete space can repeat a config, and the repeat still spends budget.

Objective contract
------------------
`objective(family, config, fold_plan, seed)` is called once per trial.

- `config` is a private copy. For trial 0 it is the family default; when
  that is an internal-search default (`config["internal_search"] is True`,
  scm and cart) the objective must run the production code path, whose own
  inner search picks among the listed values. Every other config fixes its
  values and is fitted as given.
- `fold_plan` is a private copy of the plan whose `fold_plan_hash` the
  protocol pins. It must be plain JSON with an `inner_folds` key; see
  "Fold plans" in `protocol.py`.
- `seed` is the protocol's `seeds.fit`, the value each record's `seed`
  states.

It returns an object with exactly `status`, `inner_scores`, `fit_count`,
`seconds` and `complexity`, plus `reason` when status is not `ok`, each as
documented for history lines in `protocol.py`. Inner scores may come in any
fold order; the runner sorts them. An `ok` result must score exactly the
plan's inner folds, and a `failed` or `unconverged` one a subset of them. A
result that breaks the contract, or an exception, is recorded as a
`failed` trial with nulls and the reason, and the search continues.

Resume
------
Each line is flushed and fsynced as it is written. KeyboardInterrupt and
other non-Exception errors stop the run; rerunning with the same history
file runs only the trials not yet written. A crash mid-write can leave a
partial last line without its newline. The runner refuses such a file
rather than guess; delete everything after the last newline, then rerun.
"""
import copy
import hashlib
import math
import os

from .protocol import (ProtocolError, canonical_json, history_line, inner_fold_ids,
                       read_history, validate_record, verify_protocol)
from .protocol import fold_plan_hash as compute_fold_plan_hash

RESULT_FIELDS = frozenset({"status", "inner_scores", "fit_count", "seconds", "complexity"})


def _uniform(seed, family, trial_id, param):
    digest = hashlib.sha256(canonical_json([seed, family, trial_id, param]).encode()).digest()
    return (int.from_bytes(digest[:8], "big") >> 11) / 2**53


def sample_value(spec, u):
    if spec["type"] == "choice":
        values = spec["values"]
        return values[min(math.floor(u * len(values)), len(values) - 1)]
    low, high = spec["low"], spec["high"]
    if spec["type"] == "int":
        if spec["log"]:
            value = math.floor(math.exp(math.log(low) + u * (math.log(high + 1) - math.log(low))))
        else:
            value = low + math.floor(u * (high - low + 1))
    else:
        if spec["log"]:
            value = math.exp(math.log(low) + u * (math.log(high) - math.log(low)))
        else:
            value = low + u * (high - low)
        value = float(f"{value:.10g}")
    return min(max(value, low), high)


def sample_config(protocol, family, trial_id):
    space = protocol["families"][family]["space"]
    seed = protocol["seeds"]["sampler"]
    return {param: sample_value(space[param], _uniform(seed, family, trial_id, param))
            for param in sorted(space)}


def plan_trials(protocol):
    """Every trial in run order: trial_id, then family name."""
    families = protocol["families"]
    return [
        {"family": family, "trial_id": trial_id, "is_default": trial_id == 0,
         "config": (copy.deepcopy(families[family]["default"]) if trial_id == 0
                    else sample_config(protocol, family, trial_id))}
        for trial_id in range(protocol["trials_per_family"])
        for family in sorted(families)
    ]


def check_fold_plan(protocol, fold_plan):
    """Return the plan's inner fold ids; refuse a plan the protocol does not pin."""
    digest = compute_fold_plan_hash(fold_plan)
    if digest != protocol["fold_plan_hash"]:
        raise ProtocolError(f"fold plan hashes to {digest}, not the protocol's "
                            f"{protocol['fold_plan_hash']}")
    return inner_fold_ids(fold_plan)


def _record(protocol, trial, **outcome):
    return {
        "protocol_hash": protocol["protocol_hash"],
        "fold_plan_hash": protocol["fold_plan_hash"],
        "seed": protocol["seeds"]["fit"],
        "family": trial["family"],
        "trial_id": trial["trial_id"],
        "is_default": trial["is_default"],
        "config": copy.deepcopy(trial["config"]),
        **outcome,
    }


def _failed(protocol, trial, reason):
    return _record(protocol, trial, status="failed", reason=reason, inner_scores=[],
                   fit_count=None, seconds=None, complexity=None)


def run_trial(protocol, objective, trial, fold_plan):
    """Run one trial and return its history record, failed or not.

    Refuses, without calling the objective, a fold plan whose hash is not
    the protocol's.
    """
    folds = check_fold_plan(protocol, fold_plan)
    fold_plan = copy.deepcopy(fold_plan)
    try:
        result = objective(trial["family"], copy.deepcopy(trial["config"]), fold_plan,
                           protocol["seeds"]["fit"])
    except Exception as err:
        return _failed(protocol, trial, f"objective raised {type(err).__name__}: {err}")
    try:
        if not isinstance(result, dict):
            raise ProtocolError(f"expected an object, got {type(result).__name__}")
        missing = sorted(RESULT_FIELDS - set(result))
        unknown = sorted(set(result) - RESULT_FIELDS - {"reason"})
        if missing or unknown:
            raise ProtocolError(f"missing fields {missing}, unknown fields {unknown}")
        outcome = copy.deepcopy(result)
        try:
            outcome["inner_scores"] = sorted(outcome["inner_scores"], key=lambda e: e["fold"])
        except (TypeError, KeyError):
            pass  # validate_record names the malformed entry
        record = _record(protocol, trial, **outcome)
        validate_record(record, protocol, folds)
    except (ValueError, TypeError) as err:
        return _failed(protocol, trial, f"invalid objective result: {err}")
    return record


def run_search(protocol, objective, *, fold_plan, history_path, fold_plan_hash=None):
    """Run every planned trial not yet in `history_path`, appending one line each.

    Refuses a protocol that fails its hash; a fold plan that is not plain
    JSON, defines fewer than 2 inner folds, or does not hash to the
    protocol's `fold_plan_hash`; a supplied `fold_plan_hash` other than the
    protocol's; and a history that mixes protocol hashes, belongs to another
    protocol or fold plan, has a partial last line, or is not a prefix of
    this protocol's plan. Returns every record in the history after the run.
    """
    verify_protocol(protocol)
    if fold_plan_hash is not None and fold_plan_hash != protocol["fold_plan_hash"]:
        raise ProtocolError(f"fold plan hash {fold_plan_hash!r} is not the protocol's "
                            f"{protocol['fold_plan_hash']!r}")
    folds = check_fold_plan(protocol, fold_plan)
    plan = plan_trials(protocol)
    done = read_history(history_path, protocol, folds) if os.path.exists(history_path) else []
    if len(done) > len(plan):
        raise ProtocolError(f"{history_path}: {len(done)} records exceed the plan's {len(plan)}")
    for number, (record, trial) in enumerate(zip(done, plan), 1):
        if any(canonical_json(record[key]) != canonical_json(trial[key])
               for key in ("family", "trial_id", "config")):
            raise ProtocolError(f"{history_path}:{number}: record is not planned trial "
                                f"{trial['family']}:{trial['trial_id']}")
    records = list(done)
    with open(history_path, "a", encoding="ascii", newline="") as stream:
        for trial in plan[len(done):]:
            record = run_trial(protocol, objective, trial, fold_plan)
            stream.write(history_line(record))
            stream.flush()
            os.fsync(stream.fileno())
            records.append(record)
    return records
