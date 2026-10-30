#!/usr/bin/env python3
"""Evaluation revision N's validator bindings: every value report_contract.py's profile rN needs
that is unknown until revision N's checker commit exists. Each starts as an @@NAME@@ placeholder
and fails closed until a reviewed commit fills it. Nothing here is a guess, and nothing here
depends on a cohort: profile rN checks the frozen cohorts on revision N (the revision-9 fallback
of a stopped 4k cohort), and the cleaned-cohort rerun adds its cohort pins elsewhere.

usage: python3 -I revision_n.py status
       python3 -I revision_n.py get NAME
       python3 -I revision_n.py verify-checkout WORKTREE

status lists every placeholder and whether it is filled; it exits 0 only when all are filled and
well formed. get prints one filled value (JSON for an object) and exits 2 while it is a
placeholder or malformed. verify-checkout requires WORKTREE to be a clean git checkout whose HEAD
descends from RN_CHECKER_COMMIT, with the four checker files report_contract.py pins holding
their pinned SHA-256 both in WORKTREE and at RN_CHECKER_COMMIT; it exits 1 on any problem.

To fill: replace placeholder strings in VALUES with the reviewed values, in a commit that changes
VALUES and nothing else. tests/test_revision_n.py holds for a wholly unfilled or a wholly filled
file (since e354f99), so no test changes with a fill.
Two values are decisions this file checks, not free text. RN_FIELD_CHANGES must be {} or, per
key set, a list that keeps every revision-8 key (an addition); report_contract.py refuses a
removal or rename, which needs a reviewed contract change. RN_CONVERGENCE_STATUS_RULE must be
"as-690815a" for the contract to run: revision N's check-results.py judges a convergence status
as 690815a's does (IRLS0, R1). Any other value means validate_convergence needs a reviewed change.
"""
import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

PLACEHOLDER = re.compile(r"@@[A-Z0-9_]+@@\Z")
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
COMMIT = re.compile(r"[0-9a-f]{40}\Z")
CONVERGENCE_AS_690815A = "as-690815a"
KEY_SETS = ("TOP_KEYS", "METRIC_KEYS", "CONFIG_KEYS", "FOREST_KEYS", "GBT_KEYS", "FIT_KEYS", "CONVERGENCE_KEYS")

# name: (kind, what it is). The plan's placeholder table lists the same names.
TABLE = {
    "REVISION_N": ("revision", "evaluation revision N, the revision the reports and audit.py --revision carry"),
    "RN_CHECKER_COMMIT": ("commit", "revision N's checker commit (the revision bump); the four pins are its blobs"),
    "RN_PIN_CHECK_RESULTS": ("sha256", "SHA-256 of scripts/check-results.py at RN_CHECKER_COMMIT"),
    "RN_PIN_METRIC_AUDIT": ("sha256", "SHA-256 of scripts/metric-audit/audit.py at RN_CHECKER_COMMIT"),
    "RN_PIN_SUMMARIZE": ("sha256", "SHA-256 of scripts/metric-audit/summarize.py at RN_CHECKER_COMMIT"),
    "RN_PIN_RESULT_SCOPE": ("sha256", "SHA-256 of scripts/result_scope.py at RN_CHECKER_COMMIT"),
    "RN_FIELD_CHANGES": ("field_changes", "report or status fields revision N adds, as {key set: full key list}; {} for none"),
    "RN_CONVERGENCE_STATUS_RULE": ("convergence", "how revision N's checker judges convergence statuses; only as-690815a runs"),
}

VALUES = {
    "REVISION_N": "@@REVISION_N@@",
    "RN_CHECKER_COMMIT": "@@RN_CHECKER_COMMIT@@",
    "RN_PIN_CHECK_RESULTS": "@@RN_PIN_CHECK_RESULTS@@",
    "RN_PIN_METRIC_AUDIT": "@@RN_PIN_METRIC_AUDIT@@",
    "RN_PIN_SUMMARIZE": "@@RN_PIN_SUMMARIZE@@",
    "RN_PIN_RESULT_SCOPE": "@@RN_PIN_RESULT_SCOPE@@",
    "RN_FIELD_CHANGES": "@@RN_FIELD_CHANGES@@",
    "RN_CONVERGENCE_STATUS_RULE": "@@RN_CONVERGENCE_STATUS_RULE@@",
}

# The checker files report_contract.py pins, and the placeholder that holds each pin.
PIN_NAMES = {"check-results.py": "RN_PIN_CHECK_RESULTS", "metric-audit/audit.py": "RN_PIN_METRIC_AUDIT",
             "metric-audit/summarize.py": "RN_PIN_SUMMARIZE", "result_scope.py": "RN_PIN_RESULT_SCOPE"}
CONTRACT_NAMES = ["REVISION_N", *PIN_NAMES.values(), "RN_FIELD_CHANGES", "RN_CONVERGENCE_STATUS_RULE"]


class Unfilled(ValueError):
    def __init__(self, names):
        self.names = list(names)
        super().__init__("unfilled placeholder(s): " + ", ".join(f"@@{n}@@" for n in self.names))


class Malformed(ValueError):
    pass


def is_placeholder(value):
    return isinstance(value, str) and bool(PLACEHOLDER.fullmatch(value))


def checked(name, value, table=None):
    """The value, or Malformed when it is not of its kind."""
    kind, what = (table or TABLE)[name]
    integer = isinstance(value, int) and not isinstance(value, bool)
    ok = {
        "revision": lambda v: integer and v >= 9,
        "commit": lambda v: isinstance(v, str) and bool(COMMIT.fullmatch(v)),
        "sha256": lambda v: isinstance(v, str) and bool(SHA256.fullmatch(v)),
        "field_changes": lambda v: isinstance(v, dict) and set(v) <= set(KEY_SETS)
        and all(isinstance(keys, list) and keys and all(isinstance(k, str) and k for k in keys) for keys in v.values()),
        "convergence": lambda v: isinstance(v, str) and bool(v),
    }[kind](value)
    if not ok:
        raise Malformed(f"{name} = {value!r} is not a valid {kind} ({what})")
    return value


def unfilled(values=None, table=None):
    values = VALUES if values is None else values
    return [n for n in (table or TABLE) if n not in values or is_placeholder(values[n])]


def require(names, values=None, table=None, check=None):
    """{name: value} for `names`; Unfilled naming every placeholder among them, else Malformed."""
    values = VALUES if values is None else values
    missing = [n for n in names if n not in values or is_placeholder(values[n])]
    if missing:
        raise Unfilled(missing)
    return {n: (check or checked)(n, values[n]) for n in names}


def sha256_file(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_checkout(worktree, values=None):
    """Problems with WORKTREE as revision N's source; [] when there are none."""
    v = require(["RN_CHECKER_COMMIT", *PIN_NAMES.values()], values)
    worktree = Path(worktree)

    def git(*args):
        return subprocess.run(["git", "-C", str(worktree), *args], capture_output=True)

    problems = []
    status = git("status", "--porcelain", "--untracked-files=all")
    if status.returncode or status.stdout:
        problems.append(f"{worktree} is not a clean git checkout")
    if git("merge-base", "--is-ancestor", v["RN_CHECKER_COMMIT"], "HEAD").returncode:
        problems.append(f"HEAD does not descend from RN_CHECKER_COMMIT {v['RN_CHECKER_COMMIT']}")
    for rel, name in PIN_NAMES.items():
        path = worktree / "scripts" / rel
        if not path.is_file() or sha256_file(path) != v[name]:
            problems.append(f"{path} does not hold the pinned SHA-256 {v[name]} ({name})")
        blob = git("show", f"{v['RN_CHECKER_COMMIT']}:scripts/{rel}")
        if blob.returncode or hashlib.sha256(blob.stdout).hexdigest() != v[name]:
            problems.append(f"scripts/{rel} at RN_CHECKER_COMMIT does not hold the pinned SHA-256 ({name})")
    return problems


def cli(module_name, table, values, argv=None, verify=None):
    """status / get NAME [/ verify-checkout WORKTREE] over one placeholder table."""
    ap = argparse.ArgumentParser(description=sys.modules[module_name].__doc__.split("\n\n")[0], allow_abbrev=False)
    sub = ap.add_subparsers(dest="command", required=True)
    sub.add_parser("status")
    sub.add_parser("get").add_argument("name", choices=sorted(table))
    if verify is not None:
        sub.add_parser("verify-checkout").add_argument("worktree", type=Path)
    args = ap.parse_args(argv)
    module = sys.modules[module_name]
    try:
        if args.command == "status":
            missing = unfilled(values(), table)
            print(json.dumps({n: n not in missing for n in table}, indent=1))
            if missing:
                raise Unfilled(missing)
            require(list(table), values(), table, module.checked)
            return 0
        if args.command == "get":
            value = require([args.name], values(), table, module.checked)[args.name]
            print(json.dumps(value, sort_keys=True) if isinstance(value, (dict, list)) else value)
            return 0
        problems = verify(args.worktree)
        print(json.dumps({"worktree": str(args.worktree), "problems": problems}, indent=1))
        return 1 if problems else 0
    except (Unfilled, Malformed) as exc:
        print(f"{Path(module.__file__).stem}: STOP: {exc}", file=sys.stderr)
        return 2


def main(argv=None):
    return cli(__name__, TABLE, lambda: VALUES, argv, verify=verify_checkout)


if __name__ == "__main__":
    sys.exit(main())
