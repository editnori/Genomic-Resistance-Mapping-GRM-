"""The set covering machine against its two references.

Reference B is scm_reference.py, an independent implementation with GRM's
rule set. Reference A is pyscm-ml: its utility search is run on the samples
that remain before each step GRM took, and its own fit is run freely.

pyscm-ml's rules are stumps `x > t` and `x <= t`, t at each value a feature
takes among the remaining samples. On 0/1 features the stumps at t = 0 are
presence and absence. The stumps at a feature's largest value rule out every
sample (`x > max`) or none (`x <= max`); GRM has such a rule only through a
feature that is constant on the remaining samples.
"""

import math

import numpy as np
from pyscm import SetCoveringMachineClassifier
from pyscm._scm_utility import find_max

import panels as P
import scm_reference

TOLERANCE = 1e-7  # pyscm-ml's own tolerance on utilities (cpp_extensions/double_utils.h)


def dense(panel):
    x = np.zeros((panel["n"], len(panel["rows"])), dtype=np.double)
    for feature, row in enumerate(panel["rows"]):
        x[P.members(row), feature] = 1.0
    return x


def library_step(p, x_fit, y_working, argsort, remaining, feature, presence, grm_utility):
    """pyscm-ml's best rules on the remaining samples, read against GRM's rule.

    `feature` and `presence` are GRM's rule in the search's own terms: the
    rule that holds on the samples the search keeps.
    """
    utility, features, thresholds, kinds, covered, lost = find_max(p, x_fit, y_working, argsort, remaining)
    column = x_fit[remaining, feature]
    n_out = int((y_working[remaining] == 0).sum())
    n_keep = len(remaining) - n_out
    # A stump at its feature's largest remaining value rules out all or none.
    degenerate = [
        (int(k) == 0 and (int(c), int(l)) == (n_out, n_keep)) or (int(k) == 1 and (int(c), int(l)) == (0, 0))
        for k, c, l in zip(kinds, covered, lost)
    ]
    wanted_kind = 0 if presence else 1
    same_rule = any(int(f) == feature and t == 0.0 and int(k) == wanted_kind for f, t, k in zip(features, thresholds, kinds))
    # On a feature that is 1 on every remaining sample the library has no stump at 0.
    # Its stump `x > 1` on that feature then rules out the same remaining samples as absence.
    same_effect = (
        not same_rule
        and not presence
        and bool((column == 1.0).all())
        and any(int(f) == feature and t == 1.0 and int(k) == 0 for f, t, k in zip(features, thresholds, kinds))
    )
    if same_rule or same_effect:
        status = "certified"
    elif all(degenerate) and utility > grm_utility + TOLERANCE:
        status = "dominated"
    else:
        status = "other"
    return {
        "status": status,
        "same_effect_only": same_effect,
        "library_utility": float(utility),
        "library_best_rules": len(features),
        "library_best_degenerate": int(sum(degenerate)),
        "remaining_to_rule_out": n_out,
        "remaining_to_keep": n_keep,
    }


def library_fit(kind, p, max_rules, x, x_fit, y_fit, grm_rules, grm_calls):
    """pyscm-ml's own fit, run freely: reported, not judged."""
    try:
        model = SetCoveringMachineClassifier(model_type=kind, p=p, max_rules=max_rules, random_state=0).fit(x_fit, y_fit)
    except ValueError as error:
        return {"refused": str(error)}
    rules = []
    for stump in model.model_.rules:
        if stump.threshold == 0.0:
            rules.append([int(stump.feature_idx), stump.kind == "greater"])
        else:
            rules.append([int(stump.feature_idx), f"{stump.kind} {stump.threshold:g}"])
    calls = sum(1 << i for i in np.flatnonzero(model.predict(x) == 1).tolist())
    return {"rules": rules, "same_rules": rules == grm_rules, "same_calls": calls == grm_calls}


def check_fit(panel, x, fit, max_rules):
    n, rows, y, train = panel["n"], panel["rows"], panel["y"], panel["fit"]
    kind, p = fit["kind"], fit["p"]
    grm_rules = [[int(f), bool(s)] for f, s in fit["rules"]]
    grm_calls = [P.decode_bits(n, text) for text in fit["calls"]]
    out = {"kind": kind, "p": p, "rules": len(grm_rules)}

    # Reference B: the whole path.
    rules, steps = scm_reference.fit(rows, y, train, n, kind, p, max_rules)
    problems = []
    if [list(rule) for rule in rules] != grm_rules:
        problems.append("rules")
    else:
        for index, (step, mine) in enumerate(zip(fit["steps"], steps)):
            covered, lost, utility = mine
            if (step["covered"], step["lost"]) != (covered, lost):
                problems.append(f"counts of step {index + 1}")
            # GRM's float utility is the float of the exact one, and is the float expression itself.
            if step["utility"] != covered - p * lost or not math.isclose(step["utility"], float(utility), rel_tol=1e-12, abs_tol=1e-12):
                problems.append(f"utility of step {index + 1}")
        for k, calls in enumerate(grm_calls):
            if calls != scm_reference.predict(rows, n, kind, [tuple(rule) for rule in grm_rules[:k]]):
                problems.append(f"calls with {k} rules")
    out["independent"] = {"same": not problems, "differs": problems}

    # Reference A, step by step. The search keeps `keep` and rules out `rule_out`.
    members = P.members(train)
    x_fit = np.ascontiguousarray(x[members])
    y_fit = np.array([(y >> i) & 1 for i in members], dtype=np.int64)
    y_working = y_fit if kind == "conjunction" else 1 - y_fit
    one_class = len(set(y_fit.tolist())) < 2
    argsort = np.argsort(x_fit, axis=0).T.copy()
    remaining = np.arange(len(members))
    out["steps"] = []
    # A fit the library refuses is judged on Reference B alone.
    for (feature, stored_presence), step in zip([] if one_class else grm_rules, fit["steps"]):
        presence = stored_presence if kind == "conjunction" else not stored_presence
        result = library_step(p, x_fit, y_working, argsort, remaining, feature, presence, step["utility"])
        result["grm_utility"] = step["utility"]
        out["steps"].append(result)
        holds = x_fit[remaining, feature] == (1.0 if presence else 0.0)
        remaining = remaining[holds]
    out["library_refuses_one_class"] = one_class
    out["library_fit"] = library_fit(kind, p, max_rules, x, x_fit, y_fit, grm_rules, grm_calls[-1])

    statuses = [step["status"] for step in out["steps"]]
    if problems:
        out["outcome"], out["cause"] = "unexplained", "the independent implementation gives another path: " + ", ".join(problems)
    elif one_class or all(s == "certified" for s in statuses):
        out["outcome"], out["cause"] = "matched", ""
    elif all(s in ("certified", "dominated") for s in statuses):
        out["outcome"] = "explained"
        out["cause"] = "the library's stump at a feature's largest value, which rules out every remaining sample or none, has the higher utility"
    else:
        out["outcome"], out["cause"] = "unexplained", "a step is neither among the library's best rules nor below a stump at a largest value"
    return out


def check(panel, grm):
    x = dense(panel)
    fits = [check_fit(panel, x, fit, grm["max_rules"]) for fit in grm["fits"]]
    outcomes = [fit["outcome"] for fit in fits]
    steps = [step for fit in fits for step in fit["steps"]]
    differences = []
    for fit in fits:
        if fit["outcome"] == "matched":
            continue
        dominated = [i + 1 for i, step in enumerate(fit["steps"]) if step["status"] == "dominated"]
        other = [i + 1 for i, step in enumerate(fit["steps"]) if step["status"] == "other"]
        evidence = ""
        for index, step in enumerate(fit["steps"], 1):
            if step["status"] != "certified":
                evidence = (
                    f"step {index}: GRM utility {step['grm_utility']:.6g}, library utility {step['library_utility']:.6g} on "
                    f"{step['library_best_rules']} best rules, {step['library_best_degenerate']} of them at a largest value; "
                    f"{step['remaining_to_rule_out']} samples to rule out and {step['remaining_to_keep']} to keep remained"
                )
                break
        what = []
        if dominated:
            what.append(f"steps {dominated} are not the library's best rule")
        if other:
            what.append(f"steps {other} are not among the library's best rules")
        if not fit["independent"]["same"]:
            what.append("the independent implementation differs in " + ", ".join(fit["independent"]["differs"]))
        differences.append(
            {
                "unit": f"{fit['kind']} p={fit['p']:g}",
                "what": "; ".join(what),
                "cause": fit["cause"] if fit["outcome"] == "explained" else "",
                "evidence": evidence,
            }
        )
    refused = sum(fit["library_refuses_one_class"] for fit in fits)
    free = [fit["library_fit"] for fit in fits if "refused" not in fit["library_fit"]]
    return {
        "outcome": "matched" if all(o == "matched" for o in outcomes) else "unexplained" if "unexplained" in outcomes else "explained",
        "units": "fits",
        "compared": len(fits),
        "matched": outcomes.count("matched"),
        "explained": outcomes.count("explained"),
        "unexplained": outcomes.count("unexplained"),
        "detail": {
            "independent_same": sum(fit["independent"]["same"] for fit in fits),
            "rules": sum(fit["rules"] for fit in fits),
            "steps_read_against_library": len(steps),
            "steps_certified": sum(step["status"] == "certified" for step in steps),
            "steps_certified_same_effect_only": sum(step["same_effect_only"] for step in steps),
            "steps_dominated": sum(step["status"] == "dominated" for step in steps),
            "steps_other": sum(step["status"] == "other" for step in steps),
            "fits_library_refuses": refused,
            "library_free_fits": len(free),
            "library_free_same_rules": sum(f["same_rules"] for f in free),
            "library_free_same_calls": sum(f["same_calls"] for f in free),
        },
        "differences": differences,
        "fits": fits,
    }
