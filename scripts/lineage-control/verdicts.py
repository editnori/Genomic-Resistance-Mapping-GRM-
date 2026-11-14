#!/usr/bin/env python3
"""A test of each SCM rule term's association with the lab calls inside lineage clusters.

  python3 -I scripts/lineage-control/verdicts.py RELEASE_DIR > OUT.tsv
  python3 -I scripts/lineage-control/verdicts.py RELEASE_DIR --check TABLE

RELEASE_DIR is an unpacked release: `datasets/*.json`, one report each. The
output has one row per rule term. Python 3.12, standard library only. It reads
reports and fits nothing. With --check it prints nothing and exits 1 unless
TABLE is exactly the text the script prints for RELEASE_DIR; that is how a
committed table is held to its release where the release is unpacked.

The test, fixed on 2026-11-13 before any result was read (the interpretability
review of that day), is run here with one addition, the fourth verdict named
after this list:

  Context     A term is tested on the isolates the earlier terms of its own
              model left undecided, the set the greedy search chose it on. In
              a conjunction those are the isolates every earlier term holds
              for; in a disjunction, the isolates no earlier term holds for.
  Strata      The report's lineage clusters (`samples[].lineage`). A stratum
              is informative when it holds both states of the term and both
              lab calls among at least two isolates.
  Statistic   Cochran-Mantel-Haenszel over the strata, with the continuity
              correction, one degree of freedom.
  Assessable  At least MIN_STRATA informative strata.
  Verdict     Benjamini-Hochberg over the assessable terms of the release at
              ALPHA. `seen`: adjusted p under ALPHA and more resistant
              isolates than expected where the term holds. `not_shown`:
              assessable and not under ALPHA. `cannot_test`: not assessable.

The rules fixed on 2026-11-13 give those three verdicts. This script adds a
fourth, `reversed`: adjusted p under ALPHA and fewer resistant isolates than
expected where the term holds. Without it such a term would have to be called
`seen` or `not_shown`, and it is neither. No term of release
rev14-274-20261112 has it.

A rule whose feature has no card in its report stops the run, with the
dataset and the feature named: the script will not test such a rule as one
that holds for no isolate.

A verdict is an association. `seen` says the association is not explained by
the cut into lineage clusters; it says nothing about cause. It is one test
that pools the clusters, so it does not say the association is there in each
cluster. The model and the order of its terms were chosen on all isolates, so
the test is in-sample with respect to that choice. The clusters are MinHash
clusters at one threshold per dataset; they are not a phylogeny.
"""
import base64
import hashlib
import json
import math
import sys
from pathlib import Path

MIN_STRATA = 5
ALPHA = 0.05
VERDICTS = ("seen", "not_shown", "cannot_test", "reversed")
COLUMNS = (
    "dataset", "revision", "position", "feature", "presence", "input_fingerprint", "report_sha256", "isolates",
    "context_isolates", "informative_strata", "excess", "p", "q", "verdict",
)


def bitset(b64, n):
    """A feature card's presence bits (base64 of little-endian words) as an int over n isolates."""
    return int.from_bytes(base64.b64decode(b64), "little") & ((1 << n) - 1)


def table(t, y, where):
    """Counts among `where`: term and resistant, term and susceptible, no term and resistant, neither."""
    return (
        (t & y & where).bit_count(),
        (t & ~y & where).bit_count(),
        (~t & y & where).bit_count(),
        (~t & ~y & where).bit_count(),
    )


def stratified(t, y, context, strata):
    """The Cochran-Mantel-Haenszel test of term against call over `strata`, within `context`.

    Returns the informative strata, the excess of resistant isolates where the
    term holds over what the strata's margins give, and the two-sided p with
    the continuity correction (None when no stratum is informative).
    """
    excess = variance = 0.0
    informative = 0
    for stratum in strata:
        a, b, c, d = table(t, y, context & stratum)
        n = a + b + c + d
        if n < 2:
            continue
        resistant, susceptible = a + c, b + d
        v = (a + b) * (c + d) * resistant * susceptible / (n * n * (n - 1))
        if v > 0:
            informative += 1
            excess += a - (a + b) * resistant / n
            variance += v
    p = None
    if variance > 0:
        chi_square = max(abs(excess) - 0.5, 0.0) ** 2 / variance
        p = math.erfc(math.sqrt(chi_square / 2))
    return informative, excess, p


def benjamini_hochberg(pvalues):
    """Adjusted p for each of `pvalues`, in their order."""
    order = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
    adjusted = [None] * len(pvalues)
    running = 1.0
    for rank in range(len(order), 0, -1):
        i = order[rank - 1]
        running = min(running, pvalues[i] * len(pvalues) / rank)
        adjusted[i] = running
    return adjusted


def report_terms(report, report_sha256):
    """One row per rule term of a report, without q and verdict."""
    samples = report["samples"]
    n = len(samples)
    everyone = (1 << n) - 1
    y = sum(1 << i for i, s in enumerate(samples) if s["resistant"])
    by_cluster = {}
    for i, s in enumerate(samples):
        by_cluster[s["lineage"]] = by_cluster.get(s["lineage"], 0) | 1 << i
    strata = list(by_cluster.values())
    cards = {f["feature"]: f for f in report["features"]}
    scm = report["models"]["scm"]
    if scm["kind"] not in ("conjunction", "disjunction"):
        raise SystemExit(f"{report['id']}: a set covering machine of kind {scm['kind']!r}")
    context = everyone
    rows = []
    for position, rule in enumerate(scm["rules"], start=1):
        card = cards.get(rule["feature"])
        if card is None:
            raise SystemExit(f"{report['id']}: rule {position} names feature {rule['feature']}, which has no feature card in the report")
        bits = bitset(card["bits"], n)
        t = bits if rule["presence"] else ~bits & everyone
        informative, excess, p = stratified(t, y, context, strata)
        rows.append({
            "dataset": report["id"],
            "revision": report["evaluation"]["revision"],
            "position": position,
            "feature": rule["feature"],
            "presence": int(rule["presence"]),
            "input_fingerprint": report.get("input_fingerprint", ""),
            "report_sha256": report_sha256,
            "isolates": n,
            "context_isolates": context.bit_count(),
            "informative_strata": informative,
            "excess": excess,
            "p": p,
        })
        context &= t if scm["kind"] == "conjunction" else ~t & everyone
    return rows


def judge(rows):
    """Adds q and verdict to every row: Benjamini-Hochberg over the assessable ones."""
    assessable = [r for r in rows if r["informative_strata"] >= MIN_STRATA and r["p"] is not None]
    for r, q in zip(assessable, benjamini_hochberg([r["p"] for r in assessable])):
        r["q"] = q
    for r in rows:
        q = r.setdefault("q", None)
        if q is None:
            r["verdict"] = "cannot_test"
        elif q < ALPHA and r["excess"] > 0:
            r["verdict"] = "seen"
        elif q < ALPHA:
            r["verdict"] = "reversed"
        else:
            r["verdict"] = "not_shown"
    return rows


def cell(value):
    if value is None:
        return ""
    return f"{value:.6g}" if isinstance(value, float) else str(value)


def build(release):
    """The table's text for one unpacked release: a header line, then one line per rule term."""
    paths = sorted((Path(release) / "datasets").glob("*.json"))
    if not paths:
        raise SystemExit(f"no report in {release}/datasets")
    rows = []
    for path in paths:
        raw = path.read_bytes()
        rows += report_terms(json.loads(raw), hashlib.sha256(raw).hexdigest())
    judge(rows)
    lines = ["\t".join(COLUMNS)] + ["\t".join(cell(r[c]) for c in COLUMNS) for r in rows]
    return "".join(line + "\n" for line in lines)


def main(argv):
    if len(argv) == 2:
        sys.stdout.write(build(argv[1]))
        return 0
    if len(argv) == 4 and argv[2] == "--check":
        return 0 if Path(argv[3]).read_bytes() == build(argv[1]).encode("utf-8") else 1
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
