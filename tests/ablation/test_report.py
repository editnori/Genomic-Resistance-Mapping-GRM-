"""build_report: planned-contrast completeness, family membership,
denominators, prediction files, summary wording, escaping and JSON round trips.

usage: python3 -I -B -m unittest discover -s tests/ablation -t .
"""
import contextlib
import copy
import errno
import io
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest
from collections import OrderedDict, UserList, deque
from pathlib import Path
from unittest import mock

from scripts.ablation import report as rp
from scripts.ablation import schema
from scripts.ablation.schema import SchemaError
from tests.ablation import test_schema as fx

ROOT = Path(__file__).resolve().parents[2]

# Phrases the summary must never use outside its reading rules.
FORBIDDEN = [
    r"\bno value\b", r"\bworthless\b", r"\buseless\b", r"\bvalueless\b",
    r"\bno (?:predictive )?(?:effect|contribution|benefit|signal)\b",
    r"(?<!not )\bunimportant\b", r"\bimportan(?:t|ce)\b",
    r"(?<!not )\bcaus(?:al|es|ed|ing|ation)\b",
    r"\b(?:drives|determines|explain\w*|mechanism|mechanisms)\b",
    r"does not settle",
    r"equivalen",
]
# Summary lines for pooled records and per-dataset contrasts, as render_markdown prints them.
POOLED_LINE = re.compile(r"- `[^`/]+` \(family [ABD]: ")
CONTRAST_LINE = re.compile(r"- `[^`/]+/[^`]+` \(")
# Held literally, so a change to the sentence in report.py fails a test.
FEATURE_SENTENCE = ("A block attribution is descriptive only: selection frequency and sign stability "
                    "say what the fitted models used on these data, not why resistance occurs.")
NOT_PRODUCED_REASON = "vme not computed: pairing refused (fixture)"
# Contract section 15's exact sentences, held literally.
WITHIN_MARGIN_SENTENCE = ("The 90% interval lies within ±margin. This is descriptive only: GRM's calibration did not "
                          "hold the false-equivalence rate at 5%, so it is not an equivalence finding.")
POOLED_P_SENTENCE = ("Pooled p tests 'no effect in any species'. It is not a test that the average species effect "
                     "is nonzero.")
LOBO = {"baseline": "K+C", "treatment": "K"}


def cid(dataset, **parts):
    return f"{dataset}/{fx.key(**parts)}"


FWD_A = cid("ecoli__amp")
FWD_C = cid("ecoli__cip")
LOBO_A = cid("ecoli__amp", **LOBO)
LOBO_C = cid("ecoli__cip", **LOBO)
D_MTB = cid("mtb__inh", family="D")
VME_C = cid("ecoli__cip", family="D", endpoint="vme")
VME_A = cid("ecoli__amp", family="D", endpoint="vme")
RAND_A = cid("ecoli__amp", family="D", split="random")
RAND_K = cid("kpneu__amp", family="D", split="random")
# Section 13: A and B pool the family D r2_ll contrasts of their comparison.
R2_A, R2_C, R2_K, R2_MTB = (cid(d, family="D", endpoint="r2_ll") for d in
                            ("ecoli__amp", "ecoli__cip", "kpneu__amp", "mtb__inh"))
R2L_A, R2L_C, R2L_K = (cid(d, family="D", endpoint="r2_ll", **LOBO) for d in ("ecoli__amp", "ecoli__cip", "kpneu__amp"))
# Section 17: a defined pool needs 2 contributing species, so two more confirmatory species feed A and B.
R2_S, R2_T = (cid(d, family="D", endpoint="r2_ll") for d in ("spneu__ery", "saur__oxa"))
R2L_S, R2L_T = (cid(d, family="D", endpoint="r2_ll", **LOBO) for d in ("spneu__ery", "saur__oxa"))
RF_A, RF_S = (cid(d, family="D", endpoint="r2_ll", learner="rf") for d in ("ecoli__amp", "spneu__ery"))
PLANNED = [FWD_A, FWD_C, cid("kpneu__amp"), LOBO_A, LOBO_C, cid("kpneu__amp", **LOBO), D_MTB, VME_C, VME_A,
           RAND_A, RAND_K, R2_A, R2_C, R2_K, R2_MTB, R2L_A, R2L_C, R2L_K, R2_S, R2_T, R2L_S, R2L_T, RF_A, RF_S]
POOL_A = fx.key("A", endpoint="r2_ll")
POOL_A_RF = fx.key("A", endpoint="r2_ll", learner="rf")
POOL_B = fx.key("B", endpoint="r2_ll", **LOBO)
POOL_B_NP = fx.key("B", "K+C", "C", "r2_ll")
POOL_D = fx.key("D", endpoint="r2_ll")
POOL_VME = fx.key("D", endpoint="vme")
POOL_RAND = fx.key("D", split="random")
PLANNED_POOLED = [POOL_A, POOL_A_RF, POOL_B, POOL_B_NP, POOL_D, POOL_VME, POOL_RAND]
RF_REASON = "rf runs not attempted (fixture)"
POOLED_NOT_PRODUCED_REASON = "C-only arm not run (fixture)"
NOT_PRODUCED = {VME_A: NOT_PRODUCED_REASON, cid("kpneu__amp"): "K+C run unavailable: one-class test part",
                cid("kpneu__amp", **LOBO): "K+C run unavailable", POOL_B_NP: POOLED_NOT_PRODUCED_REASON}


def good_manifest():
    """The schema fixture's manifest plus two confirmatory species with K and K+C lineage runs."""
    m = fx.good_manifest()
    m["datasets"] += [fx.good_dataset(id="spneu__ery", species="Streptococcus pneumoniae", n_samples=200,
                                      group_counts={"lineage": 8}),
                      fx.good_dataset(id="saur__oxa", species="Staphylococcus aureus", n_samples=220,
                                      group_counts={"lineage": 9})]
    n = 16
    for d in ("spneu__ery", "saur__oxa"):
        for arm, arm_blocks in (("K", ["K"]), ("K+C", ["K", "C"])):
            n += 1
            m["runs"].append(fx.good_run(id=f"r{n}", dataset=d, arm=arm, blocks=arm_blocks,
                                         predictions=f"predictions/r{n}.jsonl"))
    return m


def _family(name):
    parts = schema.parse_key(name.partition("/")[2] if "/" in name else name) if isinstance(name, str) else None
    return parts and parts["family"]


def _holm(ps):
    """Definition form: the k-th smallest gets min(1, max over j <= k of (m - j + 1) p_(j))."""
    m, ranked = len(ps), sorted(range(len(ps)), key=ps.__getitem__)
    out = [0.0] * m
    for k, i in enumerate(ranked):
        out[i] = min(1.0, max((m - j) * ps[ranked[j]] for j in range(k + 1)))
    return out


def _bh(ps):
    """Definition form: the k-th smallest gets min(1, min over j >= k of m p_(j) / j)."""
    m, ranked = len(ps), sorted(range(len(ps)), key=ps.__getitem__)
    out = [0.0] * m
    for k, i in enumerate(ranked):
        out[i] = min(1.0, min(ps[ranked[j]] * (m / (j + 1)) for j in range(k, m)))
    return out


def adjusted(contrasts, pooled, planned, planned_pooled):
    """Copies of contrasts and pooled with each family A to C p_adj set as section 17 recomputes it,
    for tests whose subject isn't the adjustment. Inputs a test made malformed are left alone."""
    if not all(isinstance(v, list) for v in (contrasts, pooled, planned, planned_pooled)):
        return contrasts, pooled
    contrasts, pooled = copy.deepcopy(contrasts), copy.deepcopy(pooled)
    recs = {x.get("id"): x for x in contrasts if isinstance(x, dict)}
    recs.update({x.get("contrast"): x for x in pooled if isinstance(x, dict)})
    for family, method in (("A", _holm), ("B", _holm), ("C", _bh)):
        names = [n for n in planned + planned_pooled if _family(n) == family]
        ps = [recs[n]["p"] if n in recs and isinstance(recs[n].get("p"), float) else 1.0 for n in names]
        for n, value in zip(names, method(ps)):
            if n in recs and recs[n].get("p_adj") is not None:
                recs[n]["p_adj"] = value
    return contrasts, pooled


def r2(ident, estimate, ci, tost, decision="inconclusive", **over):
    return fx.good_contrast(id=ident, estimate=estimate, ci=ci, ci_tost=tost, p=None, p_adj=None,
                            decision=decision, interval_label="95% jackknife-t CI, r2_ll", **over)


def good_inputs():
    """Seventeen contrasts (one undefined) of twenty-four planned; seven not produced, explained
    by a pooled exclusion or by not_produced. Pooled, seven planned: A and B over the family D
    r2_ll contrasts of their comparison on three species, with the family C log_loss contrasts
    reported beside them; an A pool with every dataset excluded; a B key not produced; and three
    family D pools, one with every dataset excluded. p_adj is BH over the six planned family C
    contrasts and Holm over the two planned A and two planned B keys, missing at p = 1 (section 17)."""
    contrasts = [
        fx.good_contrast(p_adj=0.024),
        fx.good_contrast(id=FWD_C, estimate=0.017, ci=[0.005, 0.030], ci_tost=[0.008, 0.027], p=0.03, p_adj=0.09,
                         decision="inconclusive", interval_label="95% lineage-cluster CI, cip"),
        fx.good_contrast(id=LOBO_A, estimate=-0.004, ci=[-0.012, 0.004], ci_tost=[-0.010, 0.002], p=0.30,
                         p_adj=0.6, decision="inconclusive", interval_label="95% lineage-cluster CI, leave C out"),
        fx.good_contrast(id=LOBO_C, estimate=-0.01, ci=[-0.030, 0.010], ci_tost=[-0.026, 0.006], p=0.5, p_adj=0.75,
                         decision="inconclusive", interval_label="95% lineage-cluster CI, leave C out, cip"),
        fx.good_contrast(id=D_MTB, estimate=0.06, ci=[0.03, 0.09], ci_tost=[0.035, 0.085], p=None, p_adj=None,
                         decision="inconclusive", few_clusters=True, n_units=4, unit="sublineage",
                         interval_label="95% sublineage CI, few clusters"),
        fx.good_contrast(id=VME_C, **fx.UNDEFINED, undefined_resamples=64, interval_label="95% CI, vme"),
        fx.good_contrast(id=RAND_A, estimate=0.01, ci=[-0.01, 0.03], ci_tost=[-0.005, 0.025], p=None, p_adj=None,
                         decision="inconclusive", interval_label="95% random-split CI"),
        fx.good_contrast(id=RAND_K, estimate=0.02, ci=[0.0, 0.04], ci_tost=[0.005, 0.035], p=None, p_adj=None,
                         decision="inconclusive", interval_label="95% random-split CI, kpneu"),
        r2(R2_A, 0.03, [0.01, 0.05], [0.015, 0.045]),
        r2(R2_C, 0.01, [-0.01, 0.03], [-0.005, 0.025]),
        r2(R2_MTB, 0.05, [0.02, 0.08], [0.025, 0.075], n_units=4, unit="sublineage", few_clusters=True),
        r2(R2L_A, -0.005, [-0.02, 0.01], [-0.015, 0.005]),
        r2(R2L_C, -0.01, [-0.03, 0.01], [-0.025, 0.005]),
        r2(R2_S, 0.02, [0.0, 0.04], [0.005, 0.035]),
        r2(R2_T, 0.025, [0.005, 0.045], [0.01, 0.04]),
        r2(R2L_S, -0.002, [-0.03, 0.016], [-0.025, 0.011]),
        r2(R2L_T, -0.004, [-0.026, 0.018], [-0.022, 0.014]),
    ]
    pooled = [
        fx.good_pooled(),
        fx.good_pooled(contrast=POOL_B, family="B", estimate=-0.007, ci=[-0.03, 0.01], ci_tost=[-0.025, 0.005],
                       p=0.375, p_adj=0.75, interval_label="95% species-level CI, leave C out",
                       excluded=[{"dataset": "kpneu__amp", "reason": "K+C run unavailable"}]),
        fx.good_pooled(contrast=POOL_D, family="D", n_units=4, n_datasets=6, n_excluded=1,
                       excluded=[{"dataset": "kpneu__amp", "reason": "K+C run unavailable: one-class test part"}],
                       estimate=0.04, ci=[0.01, 0.07], ci_tost=[0.015, 0.065], p=None, p_adj=None, adjust="none",
                       interval_label="95% species-level CI, all partitions"),
        fx.all_excluded_pooled(contrast=POOL_VME, family="D", adjust="none", n_datasets=2, n_excluded=2,
                               interval_label="95% species-level CI (not computed)",
                               excluded=[{"dataset": "ecoli__amp", "reason": NOT_PRODUCED_REASON},
                                         {"dataset": "ecoli__cip", "reason": fx.UNDEFINED["undefined_reason"]}]),
        fx.good_pooled(contrast=POOL_RAND, family="D", n_units=2, n_datasets=2, n_excluded=0, excluded=[],
                       estimate=0.015, ci=[0.0, 0.035], ci_tost=[0.005, 0.015], p=None, p_adj=None, adjust="none",
                       interval_label="95% species-level CI, random split"),
        fx.all_excluded_pooled(contrast=POOL_A_RF, n_datasets=2, n_excluded=2,
                               interval_label="95% species-level CI (not computed), rf",
                               excluded=[{"dataset": "ecoli__amp", "reason": RF_REASON},
                                         {"dataset": "spneu__ery", "reason": RF_REASON}]),
    ]
    features = [
        fx.good_feature(),
        fx.good_feature(dataset="ecoli__cip", column="kmer:ACGTTGCA", blocks=["K"], selection_freq=0.4,
                        sign_stable=False, claim_level="predictive_increment"),
        schema.annotate_claim(fx.good_feature(dataset="mtb__inh", column="katG_S315T"), "feature",
                              "causal_supported", ["PMID:00000000 (fixture)"], "Fixture Annotator", "2026-10-30"),
    ]
    return (good_manifest(), contrasts, pooled, features, list(PLANNED), list(PLANNED_POOLED),
            dict(NOT_PRODUCED))


NAMES = ("manifest", "contrasts", "pooled", "features", "planned_contrasts", "planned_pooled", "not_produced")


def build(adjust=True, **parts):
    """build_report on the fixture with parts replaced; adjust sets p_adj as section 17 recomputes it."""
    values = dict(zip(NAMES, good_inputs()))
    values.update(parts)
    if adjust:
        values["contrasts"], values["pooled"] = adjusted(values["contrasts"], values["pooled"],
                                                         values["planned_contrasts"], values["planned_pooled"])
    return rp.build_report(*(values[n] for n in NAMES), prediction_root=values.get("prediction_root"))


def body(markdown):
    """The summary without its reading rules, which quote the excluded phrasing."""
    head, rules = markdown.split("\n" + rp.RULES_HEADING + "\n")
    return head


def line_for(markdown, ident):
    (line,) = [x for x in markdown.splitlines() if x.startswith(f"- `{ident}` (")]
    return line


def write_predictions(root, manifest, tag=""):
    """One JSON Lines file per ok run: three samples, every repeat, the run's plan per repeat.
    tag is appended to each sample id; the files are raw UTF-8."""
    for u in manifest["runs"]:
        if u["status"] != "ok":
            continue
        rows = []
        for repeat in range(u["repeats"]):
            for i, y in enumerate((0, 1, 1)):
                rows.append(fx.good_prediction(dataset=u["dataset"], arm=u["arm"], learner=u["learner"],
                                               split=u["split"], repeat=repeat, sample_id=f"S{i}{tag}", y=y,
                                               fold=(i + repeat) % 3, score=0.8 if y else 0.2, called=bool(y),
                                               cluster=f"L{i}", fold_plan_sha256=u["fold_plan_sha256"][repeat]))
        path = Path(root) / u["predictions"]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("".join(json.dumps(r, ensure_ascii=False) + "\n" for r in rows), encoding="utf-8")


class BuildCases(unittest.TestCase):
    def setUp(self):
        self.m, self.c, self.p, self.f, self.planned, self.planned_pooled, self.np = good_inputs()
        self.report = rp.build_report(self.m, self.c, self.p, self.f, self.planned, self.planned_pooled, self.np)
        self.md = rp.render_markdown(self.report)

    def refuse(self, message, **parts):
        with self.assertRaisesRegex(SchemaError, message):
            build(**parts)

    def test_attempts(self):
        self.assertEqual(self.report["attempts"], {"runs": 21, "ok": 19, "unavailable": 2})
        self.assertIn("21 runs attempted: 19 ok, 2 unavailable.", self.md)

    def test_unavailable_runs_and_blocks_are_reported(self):
        reasons = {u["id"]: u["reason"] for u in self.report["unavailable_runs"]}
        self.assertEqual(reasons, {"r10": "one-class test part in fold 2", "r-plasmid": "block P unavailable"})
        for run_id, reason in reasons.items():
            self.assertRegex(self.md, rf"- `{run_id}` \(.*\): {re.escape(reason)}\n")
        self.assertIn("- `P`: replicon database not licensed for this run", self.md)

    def test_every_planned_contrast_is_reported(self):
        """M1: each planned id is a record or not produced with its reason."""
        self.assertEqual(self.report["contrast_counts"],
                         {"planned": 24, "produced": 17, "defined": 16, "undefined": 1, "not_produced": 7})
        self.assertEqual(self.report["not_produced"], [
            {"id": cid("kpneu__amp"), "reason": "K+C run unavailable: one-class test part"},
            {"id": cid("kpneu__amp", **LOBO), "reason": "K+C run unavailable"},
            {"id": VME_A, "reason": NOT_PRODUCED_REASON},
            {"id": R2_K, "reason": "K+C run unavailable: one-class test part"},  # from the pooled exclusions
            {"id": R2L_K, "reason": "K+C run unavailable"},
            {"id": RF_A, "reason": RF_REASON},
            {"id": RF_S, "reason": RF_REASON},
        ])
        self.assertIn("24 contrasts planned: 17 produced (16 with an estimate, 1 undefined), 7 not produced.",
                      self.md)
        for x in self.report["not_produced"]:
            self.assertIn(f"- `{x['id']}`: not produced: {x['reason']}", self.md)

    def test_dropped_contrast_without_reason_is_refused(self):
        contrasts = [x for x in self.c if x["id"] != LOBO_C]
        pooled = [p for p in self.p if p["contrast"] != POOL_B]
        self.refuse(f"{re.escape(LOBO_C)} has no record and no reason; give one in not_produced",
                    contrasts=contrasts, pooled=pooled)
        report = build(contrasts=contrasts, pooled=pooled,
                       not_produced={**self.np, LOBO_C: "dropped by fixture", R2L_K: "no run",
                                     POOL_B: "pool dropped by fixture"})
        self.assertIn({"id": LOBO_C, "reason": "dropped by fixture"}, report["not_produced"])
        self.assertIn({"contrast": POOL_B, "reason": "pool dropped by fixture"}, report["pooled_not_produced"])

    def test_planned_contrast_links(self):
        cases = [
            ({"contrasts": self.c + [fx.good_contrast(id=cid("ecoli__cip", split="random"))]},
             "is not in planned_contrasts"),
            ({"not_produced": {**self.np, FWD_A: "x"}}, "has a record"),
            ({"not_produced": {**self.np, cid("ecoli__amp", endpoint="zz"): "x"}},
             "is not a planned contrast or pooled key"),
            ({"not_produced": {**self.np, cid("kpneu__amp"): ""}}, "not_produced"),
            ({"planned_contrasts": self.planned + [cid("nope")]}, "unknown dataset nope"),
            ({"planned_contrasts": self.planned + [FWD_A]}, "must not repeat"),
            ({"planned_contrasts": self.planned + ["ecoli__amp/A-C"]}, r"planned_contrasts\[24\]: must match"),
            ({"planned_contrasts": "x"}, "must be a list"),
            ({"planned_contrasts": self.planned + [cid("ecoli__amp", family="A", endpoint="r2_ll")]},
             "families A and B hold pooled records only; plan their per-dataset contrasts as family C"),
        ]
        for parts, message in cases:
            with self.subTest(message=message):
                self.refuse(message, **parts)

    def test_development_dataset_planned_in_family_c(self):
        """I3: refused even when nothing was produced for it."""
        self.refuse(r"planned_contrasts\[24\]: development dataset mtb__inh is refused in family C",
                    planned_contrasts=self.planned + [cid("mtb__inh")],
                    not_produced={**self.np, cid("mtb__inh"): "not run"})

    def test_undefined_contrasts_are_printed_with_their_reason(self):
        reason = fx.UNDEFINED["undefined_reason"]
        self.assertEqual(self.report["undefined_contrasts"], [{"id": VME_C, "reason": reason}])
        line = line_for(self.md, VME_C)
        self.assertIn(f"undefined: {reason}. No estimate, interval or test", line)
        self.assertIn("Interval label: 95% CI, vme.", line)
        self.assertNotIn("Claim level", line)

    def test_exclusions_are_reported(self):
        self.assertEqual(len(self.report["exclusions"]), 7)
        for p in self.p:
            for e in p["excluded"]:
                self.assertIn({"pooled": p["contrast"], **e}, self.report["exclusions"])
                self.assertIn(f"- `{p['contrast']}`: `{e['dataset']}`: {e['reason']}", self.md)
                self.assertIn(f"  - excluded `{e['dataset']}`: {e['reason']}", self.md)

    def test_denominators(self):
        self.assertIn("Denominator: 5 datasets, 4 contributing, 1 excluded.", line_for(self.md, POOL_A))
        self.assertIn("Denominator: 6 datasets, 5 contributing, 1 excluded.", line_for(self.md, POOL_D))
        self.assertIn("Denominator: 2 datasets, 2 contributing, 0 excluded.", line_for(self.md, POOL_RAND))
        pooled = copy.deepcopy(self.p)
        pooled[4].update(n_datasets=1, n_excluded=1, n_units=0, excluded=[{"dataset": "ecoli__amp", "reason": "x"}],
                         estimate=None, ci=None, ci_tost=None)
        planned = [x for x in self.planned if x != RAND_K]
        contrasts = [x for x in self.c if x["id"] != RAND_K]
        md = rp.render_markdown(build(pooled=pooled, planned_contrasts=planned, contrasts=contrasts))
        self.assertIn("no pooled estimate; all 1 dataset in the denominator were excluded", line_for(md, POOL_RAND))
        self.assertIn(f"- `{POOL_VME}` (family D: K+C vs K, lasso, lineage split, vme): no pooled estimate; all 2 "
                      "datasets in the denominator were excluded (adjust none; margin ±0.02). Interval label: "
                      "95% species-level CI (not computed).", self.md)

    def test_pooled_records_link_by_comparison(self):
        """Section 13: a pool's inputs share its learner, split, arms and endpoint, whatever their family;
        A and B draw on the family D r2_ll contrasts, confirmatory datasets only."""
        pooled = copy.deepcopy(self.p)
        pooled[0].update(excluded=[], n_excluded=0)
        self.refuse("dataset kpneu__amp has no per-dataset record; exclude it with a reason", pooled=pooled)
        pooled = copy.deepcopy(self.p)
        pooled[0].update(n_datasets=4)
        self.refuse("must count every planned dataset of the comparison, 5", pooled=pooled)
        pooled = copy.deepcopy(self.p)
        pooled[0]["excluded"].append({"dataset": "mtb__inh", "reason": "x"})
        pooled[0].update(n_excluded=2, n_datasets=6)
        self.refuse("development dataset mtb__inh is refused in family A", pooled=pooled)
        # The same comparison pooled as family D counts the development dataset.
        self.assertEqual(self.p[2]["n_datasets"], 6)
        stray = fx.all_excluded_pooled(contrast=fx.key("A", endpoint="r2_ll", split="random"), excluded=[
            {"dataset": "ecoli__amp", "reason": "x"}])
        self.refuse(r"pools a comparison with no planned per-dataset contrasts \(lasso:random:K:K\+C:r2_ll\)",
                    pooled=self.p + [stray])
        # The family C log_loss contrasts are not a pool's inputs: without the D r2_ll ones, A has none.
        planned = [x for x in self.planned if "r2_ll" not in x]
        contrasts = [x for x in self.c if x["endpoint"] != "r2_ll"]
        self.refuse(r"pooled\[0\].contrast: pools a comparison with no planned per-dataset contrasts",
                    planned_contrasts=planned, contrasts=contrasts)

    def test_pooled_inputs_from_j_shape(self):
        """An undefined D r2_ll input (J's all-susceptible dataset) must be excluded from an A pool."""
        contrasts = copy.deepcopy(self.c)
        r2c = next(x for x in contrasts if x["id"] == R2_C)
        r2c.update(fx.UNDEFINED, undefined_reason="r2_ll undefined: no resistant (y=1) samples")
        self.refuse(f"undefined contrast {re.escape(R2_C)} must be excluded", contrasts=contrasts)

    def test_undefined_member_must_be_excluded_from_a_d_pool(self):
        contrasts = copy.deepcopy(self.c)
        next(x for x in contrasts if x["id"] == R2_MTB).update(
            fx.UNDEFINED, undefined_reason="r2_ll undefined: no resistant (y=1) samples")
        self.refuse(rf"pooled\[2\].excluded: undefined contrast {re.escape(R2_MTB)} must be excluded",
                    contrasts=contrasts)

    def test_pooled_units_equal_contributing_species(self):
        """L5: n_units below the species count is refused, as well as above it."""
        for n_units in (1, 0):
            pooled = copy.deepcopy(self.p)
            pooled[4]["n_units"] = n_units
            with self.subTest(n_units=n_units), self.assertRaisesRegex(
                    SchemaError, "n_units: (a defined pooled record needs at least 2|species must come)"):
                build(pooled=pooled)
        pooled = copy.deepcopy(self.p)
        pooled[0]["n_units"] = 4
        self.refuse("must equal the 3 contributing species", pooled=pooled)
        pooled = copy.deepcopy(self.p)
        pooled[2]["n_units"] = 3
        self.refuse("must equal the 4 contributing species", pooled=pooled)

    def test_development_datasets_refused_in_family_c(self):
        extra = fx.shaped("inside_margin", "C", id=cid("mtb__inh"))
        self.refuse("development dataset mtb__inh is refused in family C",
                    contrasts=self.c + [extra], planned_contrasts=self.planned + [extra["id"]])

    def test_development_datasets_allowed_in_family_d(self):
        self.assertIn(D_MTB, [x["id"] for x in self.report["contrasts"]])
        self.assertIn("(Mycobacterium tuberculosis, development)", self.md)

    def test_contrast_links(self):
        cases = [
            (fx.shaped("inside_margin", "D", id=cid("kpneu__amp", family="D")),
             "arm K\\+C has no ok run for kpneu__amp"),
            (fx.good_contrast(id=cid("nope")), "unknown dataset nope"),
            (fx.good_contrast(id=cid("ecoli__amp", learner="rf")), "arm K has no ok run"),
        ]
        for extra, message in cases:
            with self.subTest(message=message):
                self.refuse(message, contrasts=self.c + [extra], planned_contrasts=self.planned + [extra["id"]])
        self.refuse("repeats contrast id", contrasts=self.c + [fx.good_contrast()])

    def test_feature_links(self):
        cases = [
            (fx.good_feature(dataset="nope"), "unknown dataset nope"),
            (fx.good_feature(blocks=["nope"]), "unknown block nope"),
            (fx.good_feature(), r"repeats \(dataset, column\)"),
        ]
        for extra, message in cases:
            with self.subTest(message=message):
                self.refuse(message, features=self.f + [extra])

    def test_every_record_is_validated(self):
        cases = [
            ({"contrasts": self.c + [dict(self.c[0], claim_level="association_stable")]},
             r"contrasts\[17\].claim_level"),
            ({"pooled": [dict(self.p[0], unit="dataset")]}, r"pooled\[0\].unit"),
            ({"features": [dict(self.f[0], selection_freq=1.5)]}, r"features\[0\].selection_freq"),
            ({"manifest": dict(self.m, binary_sha256="x")}, "manifest.binary_sha256"),
            ({"contrasts": {}}, "contrasts: must be a list"),
        ]
        for parts, message in cases:
            with self.subTest(message=message):
                self.refuse(message, **parts)

    def test_wrong_decision_is_refused(self):
        contrasts = copy.deepcopy(self.c)
        contrasts[1]["decision"] = "improves"
        self.refuse(r"contrasts\[1\].decision: is improves but the contract rule gives inconclusive",
                    contrasts=contrasts)

    def test_non_finite_and_huge_values_raise_schema_errors(self):
        m = copy.deepcopy(self.m)
        m["compute"]["actual_seconds"] = float("nan")
        self.refuse("manifest.compute.actual_seconds: must be a finite number", manifest=m)
        contrasts = copy.deepcopy(self.c)
        contrasts[0]["ci"] = [0, 10 ** 400]
        self.refuse(r"contrasts\[0\].ci", contrasts=contrasts)
        self.refuse("not_produced: keys must be non-empty strings, not 1$", not_produced={1: "x"})

    def test_inputs_are_not_mutated_or_aliased(self):
        parts = good_inputs()
        before = copy.deepcopy(parts)
        report = rp.build_report(*parts)
        self.assertEqual(parts, before)
        report["contrasts"][0]["estimate"] = 9
        self.assertEqual(parts[1][0]["estimate"], before[1][0]["estimate"])

    def test_fixture_adjustment_is_section_17s(self):
        """The fixture's literal p_adj values are what section 17 recomputes (independent Holm and BH)."""
        self.assertEqual(adjusted(self.c, self.p, self.planned, self.planned_pooled), (self.c, self.p))
        self.assertEqual([x["p_adj"] for x in self.c if x["family"] == "C"], [0.024, 0.09, 0.6, 0.75])
        self.assertEqual([p["p_adj"] for p in self.p if p["family"] in "AB"], [0.5, 0.75, None])

    def test_pooled_keys_are_planned(self):
        """Section 17 (E2E-M1): every planned pooled key is a record or not produced with a reason."""
        self.assertEqual(self.report["planned_pooled"], PLANNED_POOLED)
        self.assertEqual(self.report["pooled_counts"], {"planned": 7, "produced": 6, "not_produced": 1})
        self.assertEqual(self.report["pooled_not_produced"],
                         [{"contrast": POOL_B_NP, "reason": POOLED_NOT_PRODUCED_REASON}])
        self.assertEqual(self.report["inputs_sha256"]["planned_pooled"], schema.sha256(PLANNED_POOLED))
        self.assertIn("7 pooled records planned: 6 produced, 1 not produced (each listed in its family's section "
                      "below).", self.md)
        section_b = self.md.split("### Family B: leave-one-block-out, pooled")[1].split("###")[0]
        self.assertIn(f"- `{POOL_B_NP}` (family B: C vs K+C, lasso, lineage split, r2\\_ll): not produced: "
                      f"{POOLED_NOT_PRODUCED_REASON}.", section_b)
        self.assertEqual(self.md.count(POOLED_NOT_PRODUCED_REASON), 1)
        drop = lambda *keys: [p for p in self.p if p["contrast"] not in keys]  # noqa: E731
        cases = [
            ({"pooled": drop(POOL_B)}, rf"planned_pooled\[2\]: {re.escape(POOL_B)} has no pooled record and no "
                                       r"reason; give one in not_produced"),
            ({"pooled": drop(POOL_A, POOL_A_RF, POOL_B)}, rf"planned_pooled\[0\]: {re.escape(POOL_A)} has no"),
            ({"pooled": []}, rf"planned_pooled\[6\]: {re.escape(POOL_RAND)} has no"),
            ({"planned_pooled": [k for k in PLANNED_POOLED if k != POOL_RAND]},
             r"pooled\[4\].contrast: is not in planned_pooled"),
            ({"planned_pooled": []}, r"pooled\[0\].contrast: is not in planned_pooled"),
            ({"planned_pooled": [], "not_produced": {k: v for k, v in self.np.items() if k != POOL_B_NP}},
             r"pooled\[5\].contrast: is not in planned_pooled"),
            ({"not_produced": {**self.np, fx.key("A", "C", "K+C", "r2_ll"): "x"}},
             "is not a planned contrast or pooled key"),
            ({"not_produced": {**self.np, POOL_A: "x"}}, rf"not_produced.{re.escape(POOL_A)}: has a record"),
            ({"not_produced": {**self.np, POOL_B_NP: ""}}, r"not_produced\..*: must be a non-empty string"),
            ({"planned_pooled": "x"}, "planned_pooled: must be a list"),
            ({"planned_pooled": PLANNED_POOLED + [POOL_A]}, "planned_pooled: must not repeat"),
            ({"planned_pooled": PLANNED_POOLED + ["ecoli__amp/" + POOL_A]}, r"planned_pooled\[7\]: must match"),
            ({"planned_pooled": PLANNED_POOLED + [fx.key("C")], "not_produced": {**self.np, fx.key("C"): "x"}},
             r"planned_pooled\[7\]: family C holds per-dataset contrasts only; pooled keys are A, B or D"),
            ({"planned_pooled": PLANNED_POOLED + [fx.key("A")], "not_produced": {**self.np, fx.key("A"): "x"}},
             r"planned_pooled\[7\]: family A pools the primary endpoint r2_ll only"),
        ]
        for parts, message in cases:
            with self.subTest(message=message):
                self.refuse(message, **parts)
        # A planned pool not produced, with its reason, is reported in its family's section.
        report = build(pooled=drop(POOL_B), not_produced={**self.np, POOL_B: "pool() refused (fixture)",
                                                          R2L_K: "K+C run unavailable"})
        self.assertEqual(report["pooled_counts"], {"planned": 7, "produced": 5, "not_produced": 2})
        md = rp.render_markdown(report)
        self.assertIn(f"- `{POOL_B}` (family B: K vs K+C, lasso, lineage split, r2\\_ll): not produced: pool() "
                      "refused (fixture).", md.split("### Family B")[1].split("###")[0])
        md = rp.render_markdown(build(pooled=[], planned_pooled=[], not_produced={
            k: v for k, v in self.np.items() if k != POOL_B_NP} | {R2_K: "x", R2L_K: "x", RF_A: "x", RF_S: "x"}))
        self.assertIn("0 pooled records planned: 0 produced, 0 not produced", md)
        self.assertIn("No pooled records planned.", md)

    def test_adjustment_is_recomputed(self):
        """Section 17: p_adj is Holm in A and B and BH in C over every planned hypothesis, missing at p = 1."""
        contrasts, pooled = copy.deepcopy(self.c), copy.deepcopy(self.p)
        contrasts[1]["p_adj"] = 0.07
        self.refuse(r"contrasts\[1\].p_adj: is 0.07, but bh over the 6 planned hypotheses of family C, missing or "
                    r"undefined at p = 1, gives 0.09", adjust=False, contrasts=contrasts)
        pooled[0]["p_adj"] = 0.25
        self.refuse(r"pooled\[0\].p_adj: is 0.25, but holm over the 2 planned hypotheses of family A",
                    adjust=False, pooled=pooled)
        pooled = copy.deepcopy(self.p)
        pooled[1]["p_adj"] = 0.375  # family B's p; Holm with the not-produced B key at p = 1 gives 0.75
        self.refuse(r"pooled\[1\].p_adj: is 0.375, but holm over the 2 planned hypotheses of family B, missing or "
                    r"undefined at p = 1, gives 0.75", adjust=False, pooled=pooled)
        contrasts = copy.deepcopy(self.c)
        contrasts[1]["p_adj"] = 0.09 + 1e-13
        build(adjust=False, contrasts=contrasts)
        contrasts[1]["p_adj"] = 0.09 + 1e-11
        self.refuse(r"contrasts\[1\].p_adj", adjust=False, contrasts=contrasts)
        # The not-produced C contrasts count at p = 1: with them left out of the plan, BH has four hypotheses.
        planned = [x for x in self.planned if x not in (cid("kpneu__amp"), cid("kpneu__amp", **LOBO))]
        self.refuse(r"contrasts\[0\].p_adj: is 0.024, but bh over the 4 planned hypotheses of family C",
                    adjust=False, planned_contrasts=planned)
        # The all-excluded A pool counts at p = 1; left out of the plan (with its inputs), Holm has one
        # hypothesis. K can't tell that plan from the protocol's (the documented trusted-plan limit).
        drop = dict(planned_contrasts=[x for x in self.planned if x not in (RF_A, RF_S)],
                    planned_pooled=[k for k in PLANNED_POOLED if k != POOL_A_RF],
                    pooled=[p for p in self.p if p["contrast"] != POOL_A_RF])
        self.refuse(r"pooled\[0\].p_adj: is 0.5, but holm over the 1 planned hypotheses of family A, missing or "
                    r"undefined at p = 1, gives 0.25", adjust=False, **drop)
        self.assertEqual(build(**drop)["pooled"][0]["p_adj"], 0.25)

    def test_inputs_are_hashed(self):
        self.assertEqual(self.report["inputs_sha256"]["contrasts"], schema.sha256(self.c))
        self.assertIn(self.report["inputs_sha256"]["manifest"], self.md)
        self.assertIsNone(self.report["predictions_sha256"])
        self.assertIn("- Prediction files: not checked", self.md)


class PredictionFileCases(unittest.TestCase):
    """I4: with a prediction root, every ok run's rows are validated as part of the build."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.m = good_inputs()[0]
        write_predictions(self.root, self.m)

    def tearDown(self):
        self.tmp.cleanup()

    def edit(self, run_id, change):
        u = next(u for u in self.m["runs"] if u["id"] == run_id)
        path = self.root / u["predictions"]
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        path.write_text("".join(json.dumps(r) + "\n" for r in change(rows)))
        return u["predictions"]

    def test_checked(self):
        report = build(prediction_root=self.root)
        self.assertEqual(len(report["predictions_sha256"]), 19)
        md = rp.render_markdown(report)
        self.assertIn("- Prediction files: 19 files checked against the manifest, sha256 `predictions/r1.jsonl`", md)
        self.assertEqual(rp.validate_report(json.loads(rp.dumps(report)), self.root), report)
        self.assertEqual(rp.validate_report(json.loads(rp.dumps(report))), report)

    def test_refusals(self):
        cases = [
            ("r1", lambda rows: [r for r in rows if r["repeat"] == 0], "run r1 needs rows for repeats 0 to 1"),
            ("r1", lambda rows: rows + [dict(rows[0], arm="K+C")], "belongs to no ok run that names this file"),
            ("r2", lambda rows: [dict(r, fold_plan_sha256="e" * 64) if r["repeat"] == 0 else r for r in rows],
             "fold plan"),
            ("r1", lambda rows: [dict(rows[0], score=1.5)] + rows[1:], r"predictions/r1.jsonl\[0\].score"),
            ("r1", lambda rows: [dict(rows[0], y=1 - rows[0]["y"])] + rows[1:], "label"),
        ]
        for run_id, change, message in cases:
            write_predictions(self.root, self.m)
            self.edit(run_id, change)
            with self.subTest(message=message), self.assertRaisesRegex(SchemaError, message):
                build(prediction_root=self.root)

    def test_rows_are_strict_utf8_json(self):
        """TK-L4: a repeated key in a prediction row is refused; rows are decoded as UTF-8 (see test_ascii_locale);
        E2E6-L1: a row nested too deeply is refused with the file's path."""
        path = self.root / "predictions/r1.jsonl"
        lines = path.read_text(encoding="utf-8").splitlines()
        path.write_text("\n".join([lines[0].replace('"y": ', '"y": 1, "y": ', 1)] + lines[1:]) + "\n",
                        encoding="utf-8")
        with self.assertRaisesRegex(SchemaError, r"^predictions/r1.jsonl: duplicate keys \['y'\]$"):
            build(prediction_root=self.root)
        for depth in (600, 10000):
            path.write_text("[" * depth + "]" * depth + "\n" + "\n".join(lines[1:]) + "\n", encoding="utf-8")
            with self.subTest(depth=depth), self.assertRaisesRegex(
                    SchemaError, r"^predictions/r1.jsonl: input nested too deeply \(more than 100 levels\)$"):
                build(prediction_root=self.root)
        write_predictions(self.root, self.m, tag="é")
        self.assertIn("S0é".encode("utf-8"), path.read_bytes())
        self.assertEqual(len(build(prediction_root=self.root)["predictions_sha256"]), 19)

    def test_missing_and_escaping_paths(self):
        (self.root / "predictions/r3.jsonl").unlink()
        with self.assertRaisesRegex(SchemaError, r"predictions/r3.jsonl: \[Errno 2\]"):
            build(prediction_root=self.root)
        write_predictions(self.root, self.m)
        m = copy.deepcopy(self.m)
        m["runs"][0]["predictions"] = "../outside.jsonl"
        with self.assertRaisesRegex(SchemaError, r"\.\./outside.jsonl: must stay inside the prediction root"):
            build(manifest=m, prediction_root=self.root)

    def test_reload(self):
        report = json.loads(rp.dumps(build(prediction_root=self.root)))
        self.edit("r1", lambda rows: rows[::-1])
        with self.assertRaisesRegex(SchemaError, "report.predictions_sha256: does not match"):
            rp.validate_report(report, self.root)
        bad = copy.deepcopy(report)
        bad["predictions_sha256"].popitem()
        with self.assertRaisesRegex(SchemaError, "must name every ok run's predictions file"):
            rp.validate_report(bad)
        bad["predictions_sha256"] = {"x": "zz"}
        with self.assertRaisesRegex(SchemaError, "predictions_sha256.x: must match"):
            rp.validate_report(bad)


class WordingCases(unittest.TestCase):
    def setUp(self):
        self.report = build()
        self.md = rp.render_markdown(self.report)

    def assert_clean(self, markdown):
        """Scan the summary outside its reading rules. Section 15's sentences are exempt only where they
        belong, once each: the pooled-p sentence (it quotes the null, "no effect in any species") on pooled
        lines, the within-margin sentence (it names equivalence) on per-dataset contrast lines."""
        lines = []
        for line in body(markdown).splitlines():
            if POOLED_LINE.match(line):
                line = line.replace(POOLED_P_SENTENCE, "<section 15>", 1)
            elif CONTRAST_LINE.match(line):
                line = line.replace(WITHIN_MARGIN_SENTENCE, "<section 15>", 1)
            lines.append(line)
        text = "\n".join(lines)
        for pattern in FORBIDDEN:
            self.assertIsNone(re.search(pattern, text, re.I), f"{pattern} in summary")

    def test_forbidden_phrases_absent(self):
        self.assert_clean(self.md)

    def test_every_family_and_decision(self):
        variants = (("C", cid("ecoli__amp")), ("C", LOBO_A), ("D", D_MTB),
                    ("D", cid("ecoli__amp", family="D", split="random", **LOBO)))
        for shape in fx.SHAPES:
            contrasts = [fx.shaped(shape, family, id=ident, interval_label=f"label {shape} {ident}")
                         for family, ident in variants]
            contrasts += [fx.good_contrast(id=cid("ecoli__cip"), **fx.UNDEFINED, interval_label="label undefined C"),
                          fx.good_contrast(id=cid("ecoli__cip", family="D", endpoint="vme"), **fx.UNDEFINED,
                                           interval_label="label undefined D")]
            planned = [x["id"] for x in contrasts]
            md = rp.render_markdown(build(contrasts=contrasts, pooled=[], planned_contrasts=planned,
                                          planned_pooled=[], not_produced={}))
            self.assert_clean(md)
            for x in contrasts:
                line = line_for(md, x["id"])
                with self.subTest(shape=shape, id=x["id"]):
                    self.assertIn(rp._md(x["interval_label"]), line)
                    self.assertEqual("Claim level" in line, x["decision"] == "improves", line)
                    self.assertNotIn(POOLED_P_SENTENCE, line)
                    if x["undefined_reason"] is not None:
                        self.assertIn(f"undefined: {x['undefined_reason']}", line)
                        continue
                    self.assertIn(f"90% TOST interval [{rp._num(x['ci_tost'][0])}, {rp._num(x['ci_tost'][1])}]",
                                  line)
                    off_zero = x["ci"][0] > 0 or x["ci"][1] < 0
                    if x["decision"] == "inconclusive" and off_zero:
                        self.assertIn("inconclusive: the 95% interval excludes zero, but", line)
                    if x["decision"] == "inconclusive" and not off_zero:
                        self.assertIn("not significant: the 95% interval includes zero", line)
                    removes = x["baseline"] == "K+C"
                    self.assertEqual("redundant given the others, not unimportant" in line,
                                     removes and x["decision"] != "worse", line)
                    small = x["decision"] in ("improves", "worse") and abs(x["estimate"]) < x["margin"]
                    self.assertEqual("significant, but smaller than the margin of 0.02" in line, small, line)
                    inside = (x["unit"] == "cluster" and x["n_units"] >= 60 and x["decision"] == "inconclusive"
                              and -0.02 < x["ci_tost"][0] and x["ci_tost"][1] < 0.02)
                    self.assertEqual(WITHIN_MARGIN_SENTENCE in line, inside, line)
                    if shape == "inside_margin":
                        self.assertTrue(inside, line)
                    if x["family"] == "D":
                        self.assertIn("no p-value (family D is exploratory, adjust none); margin ±0.02", line)

    def test_inconclusive_wording_follows_its_cause(self):
        """M2."""
        self.assertNotIn("does not settle", self.md)
        self.assertIn("Decision: inconclusive: the 95% interval excludes zero, but the adjusted p-value (0.09) "
                      "does not meet 0.05.", line_for(self.md, FWD_C))
        self.assertIn("Decision: not significant: the 95% interval includes zero; a small drop",
                      line_for(self.md, LOBO_C))
        self.assertIn("inconclusive: the 95% interval excludes zero, but exploratory family D decides no "
                      "direction.", line_for(self.md, D_MTB))

    def test_no_equivalence_wording(self):
        """Section 15 withdraws equivalent and section 12's floor wording."""
        self.assertIsNone(re.search(r"\bequivalent\b", self.md))
        for gone in ("fewer than 60", "not assessed", "not assessable", "TOST interval is not inside",
                     "TOST interval lies inside"):
            self.assertNotIn(gone, self.md)
        self.assertEqual(rp.WITHIN_MARGIN_SENTENCE, WITHIN_MARGIN_SENTENCE)
        self.assertEqual(rp.POOLED_P_SENTENCE, POOLED_P_SENTENCE)
        # CN-K-L2: outside the within-margin sentence and reading rule 5, no "equivalen" anywhere.
        rest = [line.replace(WITHIN_MARGIN_SENTENCE, "", 1) if CONTRAST_LINE.match(line) else line
                for line in body(self.md).splitlines()]
        self.assertEqual([line for line in rest if "equivalen" in line.lower()], [])
        rules = self.md.split(rp.RULES_HEADING)[1]
        self.assertEqual([x for x in rules.splitlines() if "equivalen" in x], [x for x in rules.splitlines()
                                                                                if x.startswith("5. ")])
        self.assertNotIn("equivalen", rp.dumps(self.report))

    def test_undefined_list_has_reasons_only(self):
        """CN-K-L2: the Undefined: summary list holds each id and its reason, nothing else."""
        block = self.md.split("\nUndefined:\n\n")[1].split("\n\n")[0]
        self.assertEqual(block.splitlines(), [f"- `{x['id']}`: {rp._md(x['reason'])}"
                                              for x in self.report["undefined_contrasts"]])

    def test_within_margin_note(self):
        """Section 15: cluster unit, n_units >= 60, inconclusive, -margin < ci_tost[0] and ci_tost[1] < margin."""
        for ident in (LOBO_A, R2L_A):
            self.assertIn("redundant given the others, not unimportant. " + WITHIN_MARGIN_SENTENCE,
                          line_for(self.md, ident))
        for ident in (FWD_A, FWD_C, LOBO_C, D_MTB, VME_C, R2_A, R2L_C):
            self.assertNotIn("90% interval lies within", line_for(self.md, ident))
        base = next(x for x in self.report["contrasts"] if x["id"] == LOBO_A)
        self.assertEqual((base["unit"], base["n_units"], base["margin"]), ("cluster", 64, 0.02))
        base = {**base, "ci": [-0.03, 0.03]}  # room for the TOST edges below
        cases = [
            ({}, True),
            ({"n_units": 60}, True),
            ({"n_units": 59}, False),
            ({"ci_tost": [-0.02, 0.002]}, False),
            ({"ci_tost": [-0.0199, 0.002]}, True),
            ({"ci_tost": [-0.010, 0.02]}, False),
            ({"ci_tost": [-0.010, 0.0199]}, True),
            ({"unit": "lineage"}, False),
            ({"unit": "species", "n_units": 400}, False),
            ({"margin": 0.01}, False),
            ({"margin": 0.0101}, True),
            ({**fx.SHAPES["improves_small"], "decision": "improves"}, False),
            ({**fx.SHAPES["improves_small"], "p": 0.04, "p_adj": 0.08, "decision": "inconclusive"}, True),
            ({**fx.UNDEFINED, "ci_method": "none"}, False),
        ]
        for change, expected in cases:
            rec = {**base, **change}
            with self.subTest(change=change):
                schema.validate_contrast(rec)
                self.assertEqual(rp.within_margin(rec), expected)
                contrasts = [rec if x["id"] == LOBO_A else x for x in self.report["contrasts"]]
                line = line_for(rp.render_markdown(build(contrasts=contrasts)), LOBO_A)
                self.assertEqual(WITHIN_MARGIN_SENTENCE in line, expected, line)
                self.assertEqual(line.count("90% interval lies within"), int(expected), line)

    def test_within_margin_never_on_pooled(self):
        """Section 15: never on a species-pooled record, even at 60 or more species inside the margin."""
        for p in self.report["pooled"]:
            self.assertNotIn("90% interval lies within", line_for(self.md, p["contrast"]))
        pooled = good_inputs()[2]
        pooled[1].update(ci_tost=[-0.015, 0.005])
        md = rp.render_markdown(build(pooled=pooled))
        self.assertNotIn("90% interval lies within", line_for(md, POOL_B))
        for n in (60, 400):
            line = rp._pooled_line({**self.report["pooled"][1], "ci_tost": [-0.015, 0.005], "n_units": n}, {}, [])
            with self.subTest(n=n):
                self.assertNotIn("90% interval lies within", line)

    def test_pooled_p_wording(self):
        """Section 15 (CN-K-L1): the pooled-p sentence is on every pooled line with a p, once, before the S < 6
        note; never on a pooled line without a p (family D, or A and B with every dataset excluded), never on
        a not-produced pooled line and never on a per-dataset line."""
        for pool in (POOL_A, POOL_B):
            line = line_for(self.md, pool)
            self.assertIn("margin ±0.02. Decision: ", line)
            self.assertEqual(line.count(POOLED_P_SENTENCE), 1, line)
            self.assertIn("excluded. " + POOLED_P_SENTENCE + " A pooled p below 0.05 is not reachable", line)
            # Nothing else is said between the S < 6 note and the confirmatory tests.
            self.assertIn("(S = 3; the smallest p is 2/2^S = 0.25). Per-dataset confirmatory tests", line)
        for pool in (POOL_D, POOL_VME, POOL_RAND, POOL_A_RF, POOL_B_NP):
            self.assertNotIn("Pooled p tests", line_for(self.md, pool))
        self.assertIn("no pooled estimate; all 2 datasets in the denominator were excluded (adjust holm",
                      line_for(self.md, POOL_A_RF))
        self.assertEqual(self.md.count(POOLED_P_SENTENCE), 2)
        for x in self.report["contrasts"]:
            self.assertNotIn("Pooled p tests", line_for(self.md, x["id"]))
        for s in (6, 12):
            line = rp._pooled_line({**self.report["pooled"][0], "n_units": s}, {}, [])
            with self.subTest(s=s):
                self.assertIn(POOLED_P_SENTENCE, line)
                self.assertNotIn("not reachable", line)

    def test_small_significant_effects(self):
        """Section 9: |estimate| below the margin."""
        for shape, verb in (("improves_small", "K+C improves on K"), ("worse", "K+C is worse than K")):
            values = dict(fx.SHAPES[shape])
            if shape == "worse":
                values.update(estimate=-0.01, ci=[-0.016, -0.004], ci_tost=[-0.015, -0.005])
            contrasts = copy.deepcopy(self.report["contrasts"])
            contrasts[0].update(values, decision=shape.split("_")[0])
            line = line_for(rp.render_markdown(build(contrasts=contrasts)), FWD_A)
            with self.subTest(shape=shape):
                self.assertIn(f"Decision: {verb}, significant, but smaller than the margin of 0.02.", line)
        self.assertNotIn("smaller than the margin", line_for(self.md, FWD_A))

    def test_claim_level_only_on_improves(self):
        """M3."""
        for x in self.report["contrasts"]:
            self.assertEqual("Claim level" in line_for(self.md, x["id"]), x["decision"] == "improves")
        self.assertIn("Decision: K+C improves on K. Claim level predictive_increment.", line_for(self.md, FWD_A))
        for p in self.report["pooled"]:
            self.assertNotIn("Claim level", line_for(self.md, p["contrast"]))

    def test_leave_one_block_out_is_redundancy(self):
        line = line_for(self.md, LOBO_A)
        self.assertIn("Decision: not significant: the 95% interval includes zero; a small drop means the block is "
                      "redundant given the others, not unimportant.", line)
        self.assertIn("redundant given the others, not unimportant", line_for(self.md, POOL_B))
        self.assertNotIn("redundant", line_for(self.md, FWD_A))
        self.assertNotIn("redundant", line_for(self.md, POOL_A))

    def test_no_redundancy_note_on_worse(self):
        """L5: a worse leave-one-block-out line does not call the block redundant. A pooled worse needs
        6 or more species (section 17), more than the fixture has, so the pooled line is rendered alone."""
        contrasts = copy.deepcopy(self.report["contrasts"])
        contrasts[3].update(fx.SHAPES["worse"], decision="worse")
        arm_blocks = {"K": {"K"}, "K+C": {"K", "C"}}
        worse = {**self.report["pooled"][1], "estimate": -0.03, "ci": [-0.05, -0.01], "ci_tost": [-0.045, -0.015],
                 "p": 0.001, "p_adj": 0.002, "n_units": 11, "decision": "worse"}
        for line in (line_for(rp.render_markdown(build(contrasts=contrasts)), LOBO_C),
                     rp._pooled_line(worse, arm_blocks, [])):
            self.assertIn("K is worse than K+C", line)
            self.assertNotIn("redundant", line)

    def test_attribution_is_descriptive(self):
        """L6."""
        self.assertIn("\n" + FEATURE_SENTENCE + "\n", body(self.md))
        self.assertEqual(rp.FEATURE_SENTENCE, FEATURE_SENTENCE)

    def test_family_sections(self):
        text = body(self.md)
        for title in ("### Family A: confirmatory nested-chain increments, pooled",
                      "### Family B: leave-one-block-out, pooled",
                      "### Family C: per-dataset contrasts of the A and B comparisons", "### Family D: exploratory"):
            self.assertIn(title, text)
        self.assertIn("random and source splits, other learners, secondary endpoints and sensitivity analyses",
                      text)
        self.assertIn("the random-vs-lineage gap is a diagnostic, not causal biology", text)
        self.assertLess(text.index("## Pooled across species"), text.index("## Per-dataset contrasts"))

    def test_every_interval_and_label_is_printed(self):
        """L8."""
        for x in self.report["contrasts"]:
            line = line_for(self.md, x["id"])
            self.assertIn(rp._md(x["interval_label"]), line)
            if x["ci"] is not None:
                self.assertIn(f"95% CI [{rp._num(x['ci'][0])}, {rp._num(x['ci'][1])}], 90% TOST interval "
                              f"[{rp._num(x['ci_tost'][0])}, {rp._num(x['ci_tost'][1])}] "
                              f"({rp._md(x['interval_label'])})", line)
        for p in self.report["pooled"]:
            self.assertIn(rp._md(p["interval_label"]), line_for(self.md, p["contrast"]))

    def test_pooled_lines_show_family_adjustment_and_margin(self):
        """I3 and I5: margin and adjustment on every pooled line."""
        line = line_for(self.md, POOL_A)
        self.assertIn("(family A: K+C vs K, lasso, lineage split, r2\\_ll)", line)
        self.assertIn("p = 0.25, p_adj (holm) = 0.5; margin ±0.02", line)
        self.assertIn("no p-value (family D is exploratory, adjust none); margin ±0.02", line_for(self.md, POOL_D))
        for p in self.report["pooled"]:
            self.assertIn("margin ±0.02", line_for(self.md, p["contrast"]))

    def test_confirmatory_tests_beside_a_and_b(self):
        """Section 13: the family C log_loss contrasts of the comparison are listed beside A and B."""
        self.assertIn(f"Per-dataset confirmatory tests (family C, log_loss): `{FWD_A}` improves; "
                      f"`{FWD_C}` inconclusive.", line_for(self.md, POOL_A))
        self.assertIn(f"Per-dataset confirmatory tests (family C, log_loss): `{LOBO_A}` inconclusive; "
                      f"`{LOBO_C}` inconclusive.", line_for(self.md, POOL_B))
        for pool in (POOL_D, POOL_VME, POOL_RAND):
            self.assertNotIn("confirmatory tests", line_for(self.md, pool))
        contrasts = [x for x in self.report["contrasts"] if x["family"] != "C"]
        planned = [x for x in self.report["planned_contrasts"] if "/C:" not in x]
        md = rp.render_markdown(build(contrasts=contrasts, planned_contrasts=planned, not_produced={
            VME_A: NOT_PRODUCED_REASON, POOL_B_NP: POOLED_NOT_PRODUCED_REASON}))
        self.assertIn("Per-dataset confirmatory tests (family C, log_loss): none produced.", line_for(md, POOL_A))

    def test_few_species(self):
        """Section 13: S < 6 contributing species cannot reach p < 0.05; said on every such pooled line."""
        self.assertIn("A pooled p below 0.05 is not reachable with fewer than 6 species (S = 3; the smallest p is "
                      "2/2^S = 0.25).", line_for(self.md, POOL_A))
        self.assertIn("(S = 4; the smallest p is 2/2^S = 0.125).", line_for(self.md, POOL_D))
        self.assertIn("not reachable with fewer than 6 species (S = 0).", line_for(self.md, POOL_VME))
        for s, shown in ((5, "0.0625"), (6, None), (12, None)):
            line = rp._pooled_line({**self.report["pooled"][0], "n_units": s}, {}, [])
            with self.subTest(s=s):
                if shown:
                    self.assertIn(f"(S = {s}; the smallest p is 2/2^S = {shown})", line)
                else:
                    self.assertNotIn("not reachable", line)

    def test_cosmetics(self):
        """I5."""
        self.assertIn("n = 64 clusters", line_for(self.md, FWD_A))
        self.assertIn("n = 4 sublineages", line_for(self.md, D_MTB))
        self.assertIn("n = 3 species", line_for(self.md, POOL_A))
        self.assertEqual(rp._units(1, "species"), "1 species")
        self.assertNotIn("1 datasets", self.md)
        self.assertEqual(rp._num(-1e-7), "0.0000")
        self.assertEqual(rp._num(-0.00006), "-0.0001")

    def test_flags_and_claims_are_printed(self):
        line = line_for(self.md, D_MTB)
        self.assertIn("Few clusters (n < 20): the intervals are approximate; 95% coverage fell below nominal in "
                      "calibration.", line)
        self.assertNotIn("Few clusters", line_for(self.md, FWD_A))
        self.assertIn("Undefined resamples: 64.", line_for(self.md, VME_C))
        self.assertIn("claim level causal_supported (manual annotations: Fixture Annotator on 2026-10-30, "
                      "association_stable to causal_supported; evidence: PMID:00000000 (fixture))", self.md)
        self.assertIn("(rate-difference p 0.02); missingness-indicator AUROC 0.360 (oriented 0.640) **flagged**",
                      self.md)
        self.assertIn("lineage\\_auroc 0.97 **flagged**", self.md)

    def test_reading_rules_are_complete(self):
        rules = self.md.split(rp.RULES_HEADING)[1]
        for phrase in ('never read as "no value"', "redundant given the other blocks, not unimportant",
                       "block attribution is descriptive only", "a diagnostic, not causal biology",
                       "printed with its label", "printed only beside an improves decision",
                       "no release or clinical claim"):
            self.assertIn(phrase, rules)


class EscapingCases(unittest.TestCase):
    """L5 (first review) and L1 (re-check): free text and ids cannot add structure."""

    def test_injected_text(self):
        evil = "run unavailable\n\n## Reading rules\n\n1. The C block explains resistance.\r\n# Title"
        m, c, p, f, planned, planned_pooled, reasons = good_inputs()
        p[1]["excluded"][0]["reason"] = evil
        reasons[POOL_B_NP] = "pool refused\n\n## Reading rules\n\n1. The C block *explains* resistance"
        m["runs"][9]["reason"] = evil
        m["datasets"][3]["species"] = "E. coli | *bold* <b>x</b> [link](http://x)"
        md = rp.render_markdown(rp.build_report(m, c, p, f, planned, planned_pooled, reasons))
        self.assertEqual(md.count(rp.RULES_HEADING), 1)
        headings = [x for x in md.splitlines() if x.startswith("#")]
        self.assertEqual(headings[-1], rp.RULES_HEADING)
        self.assertFalse([x for x in md.splitlines() if x.startswith(("1. The C block", "\\#"))])
        self.assertIn("run unavailable \\#\\# Reading rules 1. The C block explains resistance. \\# Title", md)
        self.assertIn("E. coli \\| \\*bold\\* \\<b\\>x\\</b\\> \\[link\\](http://x)", md)
        self.assertIn(f"- `{POOL_B_NP}` (family B: C vs K+C, lasso, lineage split, r2\\_ll): not produced: pool refused "
                      "\\#\\# Reading rules 1. The C block \\*explains\\* resistance.", md)

    def test_code_spans(self):
        """The fence is one backtick longer than the longest run inside, so the span cannot close early."""
        self.assertEqual(rp._code("plain"), "`plain`")
        self.assertEqual(rp._code("a`b"), "`` a`b ``")
        self.assertEqual(rp._code("a``b"), "``` a``b ```")
        self.assertEqual(rp._code("x``[click](x.invalid)**improves**``"),
                         "``` x``[click](x.invalid)**improves**`` ```")
        self.assertEqual(rp._code("a\nb"), "`a b`")
        self.assertEqual(rp._md("a b_c"), "a b\\_c")
        self.assertEqual(rp._md("#x & y"), "\\#x \\& y")


class RoundTripCases(unittest.TestCase):
    def setUp(self):
        self.report = build()

    def test_json_round_trip(self):
        text = rp.dumps(self.report)
        loaded = json.loads(text)
        self.assertEqual(rp.validate_report(loaded), self.report)
        self.assertEqual(rp.dumps(loaded), text)

    def test_files_round_trip(self):
        with tempfile.TemporaryDirectory() as tmp:
            rp.build_report(*good_inputs(), out_dir=Path(tmp) / "out")
            self.assertEqual(rp.load_report(Path(tmp) / "out" / rp.JSON_NAME), self.report)
            self.assertEqual((Path(tmp) / "out" / rp.MARKDOWN_NAME).read_text(encoding="utf-8"),
                             rp.render_markdown(self.report))
            self.assertEqual(sorted(x.name for x in (Path(tmp) / "out").iterdir()), [rp.JSON_NAME, rp.MARKDOWN_NAME])

    def test_load_report_refuses_a_tampered_file(self):
        """TJK-4: load_report validates what it reads; E2E-L2: a repeated key is refused."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / rp.JSON_NAME
            text = rp.dumps(self.report)
            report = json.loads(text)
            report["contrasts"][2]["decision"] = "equivalent"
            path.write_text(json.dumps(report), encoding="utf-8")
            with self.assertRaisesRegex(SchemaError, r"contrasts\[2\].decision: must be one of"):
                rp.load_report(path)
            path.write_text(text.replace('"decision": "inconclusive",', '"decision": "improves", '
                                         '"decision": "inconclusive",', 1), encoding="utf-8")
            with self.assertRaisesRegex(SchemaError, r"report.json: duplicate keys \['decision'\]"):
                rp.load_report(path)
            path.write_text(text, encoding="utf-8")
            self.assertEqual(rp.load_report(path), self.report)

    def test_tampered_reports_are_refused(self):
        cases = [
            (lambda r: r["attempts"].update(ok=16), "report.attempts: does not match"),
            (lambda r: r["attempts"].update(runs=21.0), "report.attempts: does not match"),
            (lambda r: r.update(schema_version=True), "must be ablation-report schema 1"),
            (lambda r: r["contrast_counts"].update(planned=10), "report.contrast_counts: does not match"),
            (lambda r: r["not_produced"].pop(), "report.not_produced: does not match"),
            (lambda r: r["undefined_contrasts"].clear(), "report.undefined_contrasts: does not match"),
            (lambda r: r["exclusions"].pop(), "report.exclusions: does not match"),
            (lambda r: r["inputs_sha256"].update(features="0" * 64), "report.inputs_sha256: does not match"),
            (lambda r: r.update(predictions_sha256={"x": "0" * 64}), "must name every ok run's predictions file"),
            (lambda r: r.update(extra=1), "must have exactly the fields"),
            (lambda r: r["contrasts"][1].update(decision="improves"), "the contract rule gives inconclusive"),
            (lambda r: r["contrasts"][1].update(ci_tost=[-0.01, 0.03]), "must nest"),
            (lambda r: r["contrasts"][2].update(decision="equivalent"), "decision: must be one of"),
            (lambda r: r["pooled"][1].update(decision="equivalent"), "decision: must be one of"),
            (lambda r: r["contrasts"][0].update(claim_level="causal_supported"), "claim_level"),
            (lambda r: r["pooled"][3].update(family="A", adjust="holm"), "must equal the family its key implies"),
            (lambda r: r["pooled_counts"].update(not_produced=0), "report.pooled_counts: does not match"),
            (lambda r: r["pooled_not_produced"].clear(), "report.pooled_not_produced: does not match"),
            (lambda r: r["planned_pooled"].remove(POOL_B_NP), r"not_produced.*is not a planned contrast or pooled key"),
            (lambda r: r.pop("planned_pooled"), "must have exactly the fields"),
        ]
        for mutate, message in cases:
            report = json.loads(rp.dumps(self.report))
            mutate(report)
            with self.subTest(message=message), self.assertRaisesRegex(SchemaError, message):
                rp.validate_report(report)

    def test_rehashed_p_adj_forgery_is_refused(self):
        """E2E-I1: with the input hash recomputed, a lowered p_adj and a matching decision still fail section 17."""
        report = json.loads(rp.dumps(self.report))
        report["contrasts"][1].update(p_adj=0.03, decision="improves")
        report["inputs_sha256"]["contrasts"] = schema.sha256(report["contrasts"])
        with self.assertRaisesRegex(SchemaError, r"contrasts\[1\].p_adj: is 0.03, but bh"):
            rp.validate_report(report)
        report = json.loads(rp.dumps(self.report))
        report["pooled"][1]["p_adj"] = 0.375
        report["inputs_sha256"]["pooled"] = schema.sha256(report["pooled"])
        with self.assertRaisesRegex(SchemaError, r"pooled\[1\].p_adj: is 0.375, but holm over the 2 planned "
                                                 r"hypotheses of family B"):
            rp.validate_report(report)

    def test_consistent_tamper_changes_the_input_hash(self):
        report = json.loads(rp.dumps(self.report))
        report["manifest"]["runs"] = [u for u in report["manifest"]["runs"] if u["id"] != "r-plasmid"]
        report["unavailable_runs"] = [u for u in report["unavailable_runs"] if u["id"] != "r-plasmid"]
        report["attempts"] = {"runs": 20, "ok": 19, "unavailable": 1}
        report["inputs_sha256"]["manifest"] = schema.sha256(report["manifest"])
        rp.validate_report(report)
        self.assertNotEqual(report["inputs_sha256"]["manifest"], schema.sha256(good_manifest()))


class WriteCases(unittest.TestCase):
    """Section 17 and the round-6 re-check: an earlier report survives a failed write byte for byte (E2E6-M1),
    files get the mode a plain open() gives (CN6-L3), and temporary files are made in OUT_DIR (TK-L3)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = Path(self.tmp.name) / "out"
        self.earlier = rp.build_report(*good_inputs(), out_dir=self.out)
        self.before = {f: (self.out / f).read_bytes() for f in (rp.JSON_NAME, rp.MARKDOWN_NAME)}
        values = good_inputs()
        values[-1][POOL_B_NP] = "a different reason (second run)"
        self.later = rp.build_report(*values)

    def tearDown(self):
        self.tmp.cleanup()

    def listing(self):
        return sorted(x.name for x in self.out.iterdir())

    def fail_when(self, test):
        real = os.replace

        def replace(src, dst):
            if test(Path(src), Path(dst)):
                raise OSError(errno.EIO, "injected rename failure")
            return real(src, dst)
        return mock.patch.object(rp.os, "replace", side_effect=replace)

    def test_failed_renames_restore_the_earlier_report(self):
        md, js = self.out / rp.MARKDOWN_NAME, self.out / rp.JSON_NAME
        cases = [
            ("second rename (Markdown) fails", lambda src, dst: src.suffix == ".tmp" and dst == md, md),
            ("first rename (JSON) fails", lambda src, dst: src.suffix == ".tmp" and dst == js, js),
            ("Markdown backup fails", lambda src, dst: src == md and dst.suffix == ".bak", md),
        ]
        for name, test, path in cases:
            with self.subTest(case=name):
                with self.fail_when(test), self.assertRaisesRegex(
                        SchemaError, rf"^{re.escape(str(path))}: \[Errno 5\] injected rename failure$"):
                    rp.write_report(self.later, self.out)
                self.assertEqual(self.listing(), [rp.JSON_NAME, rp.MARKDOWN_NAME])
                self.assertEqual({f: (self.out / f).read_bytes() for f in self.before}, self.before)
                self.assertEqual(rp.load_report(js), self.earlier)

    def test_failed_restore_keeps_the_backup(self):
        """W7-M1: when restoring a backup fails too, that backup is never deleted, holds the earlier file byte
        for byte, and is named in the message. Three triggers: two faults in a row; renames into the final
        names failing while deletes work; a directory created at the JSON path mid-write."""
        md, js = self.out / rp.MARKDOWN_NAME, self.out / rp.JSON_NAME
        real = os.replace

        def directory_appears(src, dst):
            if src == js and dst.suffix == ".bak":  # right after the JSON backup move, another process makes a dir
                real(src, dst)
                js.mkdir()
                return "done"
            return False

        cases = [
            ("Markdown rename, then the JSON restore", lambda src, dst: (src.suffix == ".tmp" and dst == md)
             or (src.suffix == ".bak" and dst == js), [rp.JSON_NAME]),
            ("every rename into a final name", lambda src, dst: dst in (js, md) and src.suffix in (".tmp", ".bak"),
             [rp.JSON_NAME, rp.MARKDOWN_NAME]),
            ("a directory appears at the JSON path", directory_appears, [rp.JSON_NAME]),
        ]
        for name, test, lost in cases:
            out = Path(self.tmp.name) / name.replace(" ", "-").replace(",", "")
            rp.build_report(*good_inputs(), out_dir=out)
            md, js = out / rp.MARKDOWN_NAME, out / rp.JSON_NAME

            def replace(src, dst, test=test):
                verdict = test(Path(src), Path(dst))
                if verdict == "done":
                    return None
                if verdict:
                    raise OSError(errno.EIO, "injected rename failure")
                return real(src, dst)
            with self.subTest(case=name):
                with mock.patch.object(rp.os, "replace", side_effect=replace), self.assertRaises(SchemaError) as e:
                    rp.write_report(self.later, out)
                (message,) = e.exception.errors
                kept = sorted(out.glob(".*.bak"))
                self.assertEqual([x.name.split(".")[1] + "." + x.name.split(".")[2] for x in kept], lost)
                for backup, final in zip(kept, lost):
                    self.assertEqual(backup.read_bytes(), self.before[final])
                self.assertTrue(message.endswith("; earlier report kept at " + ", ".join(str(x) for x in kept)),
                                message)
                for final in set(self.before) - set(lost):
                    self.assertEqual((out / final).read_bytes(), self.before[final])
                self.assertEqual(list(out.glob(".*.tmp")), [])

    def test_fresh_folder_failure_leaves_nothing(self):
        """W7-L4: an OUT_DIR this run created, with any parents it created, is removed again on failure."""
        base = Path(self.tmp.name)
        for out, gone in ((base / "fresh", base / "fresh"), (base / "absent" / "a" / "out", base / "absent")):
            with self.subTest(out=str(out.relative_to(base))):
                with self.fail_when(lambda src, dst, out=out: dst == out / rp.MARKDOWN_NAME), self.assertRaises(
                        SchemaError):
                    rp.write_report(self.later, out)
                self.assertFalse(gone.exists())
                self.assertTrue(base.exists())
        # A folder that existed before the run stays, empty, after a failure.
        (base / "empty").mkdir()
        with self.fail_when(lambda src, dst: dst == base / "empty" / rp.MARKDOWN_NAME), self.assertRaises(SchemaError):
            rp.write_report(self.later, base / "empty")
        self.assertEqual(list((base / "empty").iterdir()), [])
        # A mkdir error other than "no such file" on a folder that exists counts as existing, never as made.
        real_mkdir = os.mkdir

        def mkdir(path, *args, **kwargs):
            if Path(path) == base / "empty":
                raise OSError(errno.EACCES, "injected mkdir failure")
            return real_mkdir(path, *args, **kwargs)
        with mock.patch.object(rp.os, "mkdir", side_effect=mkdir), self.fail_when(
                lambda src, dst: dst == base / "empty" / rp.MARKDOWN_NAME), self.assertRaises(SchemaError):
            rp.write_report(self.later, base / "empty")
        self.assertTrue((base / "empty").is_dir())

    def test_a_failed_mkdir_removes_the_parents_it_made(self):
        """RK8-L1: when making OUT_DIR itself, or a parent, fails after the folders above it were made, those are
        removed again, and the message names OUT_DIR; a mkdir that still finds no parent after making it is tried
        once more, not again and again; a folder already gone at cleanup is skipped."""
        base = Path(self.tmp.name) / "m"
        base.mkdir()
        out = base / "absent" / "a" / "out"
        real_mkdir, real_rmdir = os.mkdir, os.rmdir
        for fails in (out, out.parent):
            def mkdir(path, *args, fails=fails, **kwargs):
                if Path(path) == fails:
                    raise OSError(errno.ENOSPC, "injected mkdir failure")
                return real_mkdir(path, *args, **kwargs)
            with self.subTest(fails=str(fails.relative_to(base))):
                with mock.patch.object(rp.os, "mkdir", side_effect=mkdir), self.assertRaisesRegex(
                        SchemaError, rf"^{re.escape(str(out))}: \[Errno 28\] injected mkdir failure$"):
                    rp.write_report(self.later, out)
                self.assertEqual(list(base.iterdir()), [])
        calls = []

        def vanishing(path, *args, **kwargs):  # its parent is made, yet OUT_DIR's mkdir still finds none
            if Path(path) == out:
                calls.append(path)
                if len(calls) > 20:
                    raise OSError(errno.ELOOP, "retried without end")
                raise FileNotFoundError(errno.ENOENT, "injected: parent vanished")
            return real_mkdir(path, *args, **kwargs)
        with mock.patch.object(rp.os, "mkdir", side_effect=vanishing), self.assertRaisesRegex(
                SchemaError, rf"^{re.escape(str(out))}: \[Errno 2\] injected: parent vanished$"):
            rp.write_report(self.later, out)
        self.assertEqual((len(calls), list(base.iterdir())), (2, []))

        def rmdir(path, *args, **kwargs):  # another process removes OUT_DIR first
            real_rmdir(path, *args, **kwargs)
            if Path(path) == out:
                raise FileNotFoundError(errno.ENOENT, "injected: already gone", str(path))
        with mock.patch.object(rp.os, "rmdir", side_effect=rmdir), self.fail_when(
                lambda src, dst: dst == out / rp.MARKDOWN_NAME), self.assertRaises(SchemaError):
            rp.write_report(self.later, out)
        self.assertEqual(list(base.iterdir()), [])

    def test_cleanup_goes_on_past_a_folder_that_stays(self):
        """RK9-I1: with OUT_DIR new/b/../c and another process's file in the new c, the empty new/b this run made
        is still removed; c, and new above it, stay."""
        base = Path(self.tmp.name) / "s"
        base.mkdir()
        out = base / "new" / "b" / ".." / "c"

        def theirs(src, dst):
            if dst == out / rp.MARKDOWN_NAME:
                (base / "new" / "c" / "theirs.txt").write_text("not ours")
                return True
            return False
        with self.fail_when(theirs), self.assertRaises(SchemaError):
            rp.write_report(self.later, out)
        self.assertEqual(sorted(p.relative_to(base).as_posix() for p in base.rglob("*")),
                         ["new", "new/c", "new/c/theirs.txt"])

    def test_out_dir_that_is_a_file(self):
        """RK9-I2: an OUT_DIR that is an existing file fails at its mkdir, naming OUT_DIR, and is left as it was."""
        base = Path(self.tmp.name) / "f"
        base.mkdir()
        for out in (base / "afile", base / "absent" / "afile"):
            out.parent.mkdir(exist_ok=True)
            out.write_text("a file")
            with self.subTest(out=out.relative_to(base).as_posix()), self.assertRaisesRegex(
                    SchemaError, rf"^{re.escape(str(out))}: \[Errno 17\] File exists"):
                rp.write_report(self.later, out)
            self.assertEqual(out.read_text(), "a file")
        self.assertEqual(sorted(p.relative_to(base).as_posix() for p in base.rglob("*")),
                         ["absent", "absent/afile", "afile"])

    def test_out_dir_through_dot_dot_or_very_deep(self):
        """RK8-I1: an OUT_DIR through ".." leaves nothing it made and keeps a folder that existed; an OUT_DIR
        1,100 new levels deep is made without RecursionError and removed again on failure."""
        base = Path(self.tmp.name) / "d"
        base.mkdir()
        out = base / "new" / ".." / "x"
        for existed in (False, True):
            if existed:
                (base / "x").mkdir()
            with self.subTest(x_existed=existed):
                with self.fail_when(lambda src, dst: dst == out / rp.MARKDOWN_NAME), self.assertRaises(SchemaError):
                    rp.write_report(self.later, out)
                self.assertEqual(sorted(p.name for p in base.iterdir()), ["x"] if existed else [])
        (base / "x").rmdir()
        deep = base.joinpath(*["n"] * 1100)
        try:
            with self.fail_when(lambda src, dst: dst == deep / rp.MARKDOWN_NAME), self.assertRaisesRegex(
                    SchemaError, rf"^{re.escape(str(deep / rp.MARKDOWN_NAME))}: \[Errno 5\]"):
                rp.write_report(self.later, deep)
            self.assertEqual(list(base.iterdir()), [])
            rp.write_report(self.later, deep)
            self.assertEqual(rp.load_report(deep / rp.JSON_NAME), self.later)
        finally:
            remove_tree(base)

    def test_success_replaces_the_earlier_report_and_deletes_the_backups(self):
        rp.write_report(self.later, self.out)
        self.assertEqual(self.listing(), [rp.JSON_NAME, rp.MARKDOWN_NAME])
        self.assertEqual(rp.load_report(self.out / rp.JSON_NAME), self.later)
        self.assertNotEqual((self.out / rp.JSON_NAME).read_bytes(), self.before[rp.JSON_NAME])

    def test_a_non_regular_final_path_is_refused_before_anything_changes(self):
        """A directory, a dangling symlink (TK7-L3) or a FIFO at a final path is refused, left as it was, and the
        earlier JSON is untouched."""
        md = self.out / rp.MARKDOWN_NAME
        cases = [("directory", md.mkdir, Path.is_dir),
                 ("dangling symlink", lambda: md.symlink_to(self.out / "nowhere"), Path.is_symlink),
                 ("FIFO", lambda: os.mkfifo(md), Path.is_fifo)]
        for name, make, kind in cases:
            with self.subTest(case=name):
                if md.is_dir() and not md.is_symlink():
                    md.rmdir()
                elif md.exists() or md.is_symlink():
                    md.unlink()
                make()
                with self.assertRaisesRegex(SchemaError, rf"^{re.escape(str(md))}: exists and is not a regular file"):
                    rp.write_report(self.later, self.out)
                self.assertTrue(kind(md))
                self.assertEqual((self.out / rp.JSON_NAME).read_bytes(), self.before[rp.JSON_NAME])
                self.assertEqual(self.listing(), [rp.JSON_NAME, rp.MARKDOWN_NAME])

    def test_file_mode_follows_the_umask(self):
        """CN6-L3: 0o666 & ~umask, the mode open() gives, on a fresh folder and over an earlier pair."""
        for mask in (0o022, 0o027, 0o077):
            old = os.umask(mask)
            try:
                for out in (Path(self.tmp.name) / f"fresh-{mask:o}", self.out):
                    rp.write_report(self.later, out)
                    with self.subTest(mask=oct(mask), out=out.name):
                        for f in (rp.JSON_NAME, rp.MARKDOWN_NAME):
                            self.assertEqual((out / f).stat().st_mode & 0o777, 0o666 & ~mask)
            finally:
                os.umask(old)

    def test_temporary_and_backup_files_are_made_in_out_dir(self):
        """TK-L3: every name mkstemp makes is in OUT_DIR (two temporary files, plus two backups over a report)."""
        made, real = [], tempfile.mkstemp

        def spy(*args, **kwargs):
            fd, name = real(*args, **kwargs)
            made.append(Path(name))
            return fd, name
        fresh = Path(self.tmp.name) / "fresh"
        for out, count in ((fresh, 2), (self.out, 4)):
            made.clear()
            with mock.patch.object(rp.tempfile, "mkstemp", side_effect=spy):
                rp.write_report(self.later, out)
            with self.subTest(out=out.name):
                self.assertEqual(len(made), count)
                self.assertEqual({x.parent for x in made}, {out})
                self.assertEqual(sorted(x.name for x in out.iterdir()), [rp.JSON_NAME, rp.MARKDOWN_NAME])


class CommandLineCases(unittest.TestCase):
    def write_inputs(self, tmp):
        paths = []
        for name, value in zip(NAMES, good_inputs()):
            path = Path(tmp) / f"{name}.json"
            path.write_text(json.dumps(value))
            paths.append(str(path))
        return paths[:6], paths[6]

    def run_main(self, args):
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            code = rp.main(args)
        return code, err.getvalue()

    def test_main(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths, reasons = self.write_inputs(tmp)
            write_predictions(Path(tmp) / "preds", good_inputs()[0])
            out = str(Path(tmp) / "out")
            self.assertEqual(self.run_main(paths + [out, "--not-produced", reasons, "--predictions-root",
                                                    str(Path(tmp) / "preds")]), (0, ""))
            report = rp.load_report(Path(out) / rp.JSON_NAME, Path(tmp) / "preds")
            self.assertEqual(len(report["predictions_sha256"]), 19)

    def test_errors_use_the_documented_format(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths, reasons = self.write_inputs(tmp)
            bad = Path(tmp) / "bad"
            contrasts = json.loads(Path(paths[1]).read_text())
            contrasts[0]["estimate"] = float("nan")
            Path(paths[1]).write_text(json.dumps(contrasts))
            code, err = self.run_main(paths + [str(bad), "--not-produced", reasons])
            self.assertEqual(code, 1)
            self.assertIn("contrasts[0].estimate: must be a finite number", err)
            Path(paths[1]).write_text("{not json")
            code, err = self.run_main(paths + [str(bad), "--not-produced", reasons])
            self.assertTrue(err.startswith(f"{paths[1]}: Expecting property name"), err)
            Path(paths[1]).write_text(json.dumps(good_inputs()[1]))
            code, err = self.run_main(paths[:5] + [str(Path(tmp) / "absent.json"), str(bad)])
            self.assertEqual(code, 1)
            self.assertIn("absent.json: [Errno 2]", err)
            code, err = self.run_main(paths + [str(bad), "--not-produced", reasons, "--predictions-root",
                                               str(Path(tmp) / "nowhere")])
            self.assertEqual(code, 1)
            self.assertIn("predictions/r1.jsonl: [Errno 2]", err)
            for line in err.strip().splitlines():
                self.assertRegex(line, r"^\S+: ")
            self.assertFalse(bad.exists())
            Path(paths[1]).write_text('[{"id": "a", "id": "b"}]')
            code, err = self.run_main(paths + [str(bad), "--not-produced", reasons])
            self.assertEqual((code, err), (1, f"{paths[1]}: duplicate keys ['id']\n"))

    def test_input_nested_too_deeply(self):
        """E2E6-L1: deep input is refused with "<path>: input nested too deeply", exit 1, nothing written."""
        with tempfile.TemporaryDirectory() as tmp:
            paths, reasons = self.write_inputs(tmp)
            out = Path(tmp) / "out"
            for depth, expected in ((100, f"{paths[3]}: must be a list"), (101, None), (600, None), (5000, None),
                                    (10000, None), (50000, None)):
                Path(paths[3]).write_text("[" * depth + "]" * depth)
                code, err = self.run_main(paths + [str(out), "--not-produced", reasons])
                with self.subTest(depth=depth):
                    self.assertEqual(code, 1)
                    if expected is None:
                        self.assertEqual(err, f"{paths[3]}: input nested too deeply (more than 100 levels)\n")
                    else:
                        self.assertNotIn("nested too deeply", err)
                        self.assertNotIn("Traceback", err)
            self.assertFalse(out.exists())
            report = json.loads(rp.dumps(build()))
            report["features"] = "@@"
            path = Path(tmp) / "deep.json"
            for depth in (600, 10000):
                path.write_text(json.dumps(report).replace('"@@"', "[" * depth + "]" * depth))
                with self.subTest(depth=depth), self.assertRaisesRegex(
                        SchemaError, rf"^{re.escape(str(path))}: input nested too deeply \(more than 100 levels\)$"):
                    rp.load_report(path)

    def test_deep_python_arguments_raise_schema_errors(self):
        """W7-M2 and RK8-M1: every build_report argument, and validate_report, refuse deep Python objects with
        SchemaError, never RecursionError, at depths where copy.deepcopy would exhaust the stack. Dicts (with
        subclasses), lists, tuples and deep tuple keys get the argument's name; another container (a deque, a
        frozenset, a UserList) is refused by the deep copy, and a deep frozenset key in a derived section by the
        string-key check."""
        makers = {"dict": lambda v: {"k": v}, "list": lambda v: [v], "tuple": lambda v: (v,),
                  "OrderedDict": lambda v: OrderedDict(k=v), "deque": lambda v: deque([v]),
                  "frozenset": lambda v: frozenset([v]), "UserList": lambda v: UserList([v])}

        def nest(kind, depth):
            if kind.endswith(" key"):
                return {nest(kind[:-4], depth): "x"}
            value = {"tuple": (), "frozenset": frozenset()}.get(kind, {})  # hashable all the way down, as a key
            for _ in range(depth):
                value = makers[kind](value)
            return value
        for i, name in enumerate(NAMES):
            for kind in ("deque", "frozenset", "UserList"):
                args = list(good_inputs())
                args[i] = nest(kind, 600)
                with self.subTest(argument=name, kind=kind), self.assertRaisesRegex(
                        SchemaError, r"^build_report: input nested too deeply \(more than 100 levels\)$"):
                    rp.build_report(*args)
        args = list(good_inputs())
        args[1], args[6] = nest("list", 600), nest("dict", 600)
        with self.assertRaises(SchemaError) as caught:
            rp.build_report(*args)
        self.assertEqual(caught.exception.errors, [f"build_report: {name} nested too deeply (more than 100 levels)"
                                                   for name in ("contrasts", "not_produced")])
        for kind, key in (("frozenset", "not_produced_reasons"), ("deque", "features"),
                          ("tuple key", "not_produced_reasons"), ("tuple key", "attempts"),
                          ("frozenset key", "attempts")):
            report = json.loads(rp.dumps(build()))
            report[key] = nest(kind, 5000)
            with self.subTest(report_key=key, kind=kind), self.assertRaises(SchemaError):
                rp.validate_report(report)
        for depth in (600, 5000, 50000):  # RK9-M1: checked by validate_report itself, outside the deep copy
            report = json.loads(rp.dumps(build()))
            report["predictions_sha256"] = {nest("frozenset", depth): "0" * 64}
            with self.subTest(predictions_sha256_key_depth=depth), self.assertRaisesRegex(
                    SchemaError, r"^report\.predictions_sha256: keys must be non-empty strings, not frozenset\("):
                rp.validate_report(report)
        for i, name in enumerate(NAMES):
            for kind in ("dict", "list", "tuple", "OrderedDict", "tuple key"):
                for depth in (600, 5000):
                    args = list(good_inputs())
                    args[i] = nest(kind, depth)
                    with self.subTest(argument=name, kind=kind, depth=depth), self.assertRaisesRegex(
                            SchemaError, rf"^build_report: {name} nested too deeply \(more than 100 levels\)$"):
                        rp.build_report(*args)
        for key in ("not_produced_reasons", "features", "attempts", "inputs_sha256"):
            report = json.loads(rp.dumps(build()))
            report[key] = nest("dict", 600)
            with self.subTest(report_key=key), self.assertRaisesRegex(
                    SchemaError, r"^report: nested too deeply \(more than 100 levels\)$"):
                rp.validate_report(report)

    def test_write_failures_leave_neither_file(self):
        """Section 17 (E2E-L3): a failed write prints "<path>: <message>", exits 1 and leaves neither file."""
        with tempfile.TemporaryDirectory() as tmp:
            paths, reasons = self.write_inputs(tmp)
            out = Path(tmp) / "out-mddir"
            (out / rp.MARKDOWN_NAME).mkdir(parents=True)
            code, err = self.run_main(paths + [str(out), "--not-produced", reasons])
            self.assertEqual(code, 1)
            self.assertEqual(err, f"{out / rp.MARKDOWN_NAME}: exists and is not a regular file; {rp.MARKDOWN_NAME} "
                                  "not replaced\n")
            self.assertEqual(sorted(x.name for x in out.iterdir()), [rp.MARKDOWN_NAME])
            # An encoding error (a lone surrogate in a reason) also leaves neither file.
            pooled = json.loads(Path(paths[2]).read_text())
            pooled[1]["excluded"][0]["reason"] = "bad \udc80 reason"
            Path(paths[2]).write_text(json.dumps(pooled))
            out = Path(tmp) / "out-encoding"
            code, err = self.run_main(paths + [str(out), "--not-produced", reasons])
            self.assertEqual(code, 1)
            self.assertTrue(err.startswith(f"{out / rp.MARKDOWN_NAME}: 'utf-8' codec can't encode"), err)
            self.assertEqual(err.count("\n"), 1)
            self.assertFalse(out.exists() and any(out.iterdir()))

    def test_script_entry_point(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths, reasons = self.write_inputs(tmp)
            done = subprocess.run([sys.executable, "-I", "-B", str(ROOT / "scripts/ablation/report.py"),
                                   *paths, str(Path(tmp) / "out"), "--not-produced", reasons],
                                  capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stderr)
            self.assertTrue((Path(tmp) / "out" / rp.MARKDOWN_NAME).is_file())

    def test_ascii_locale(self):
        """Section 17 (E2E-L3): with C-locale coercion and UTF-8 mode off, both files are still UTF-8 and
        byte-identical to a default run's."""
        with tempfile.TemporaryDirectory() as tmp:
            paths, reasons = self.write_inputs(tmp)
            features = good_inputs()[3]
            features[0]["column"] = "blaTEM-1é"  # a non-ASCII input, written as raw UTF-8
            Path(paths[3]).write_text(json.dumps(features, ensure_ascii=False), encoding="utf-8")
            write_predictions(Path(tmp) / "preds", good_inputs()[0], tag="é")  # non-ASCII sample ids (TK-L4)
            script = str(ROOT / "scripts/ablation/report.py")
            env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"), "LC_ALL": "C", "PYTHONCOERCECLOCALE": "0"}
            runs = {}
            for name, extra, run_env in (("ascii", ["-X", "utf8=0"], env), ("default", [], None)):
                out = Path(tmp) / name
                done = subprocess.run([sys.executable, *extra, "-I", "-B", script, *paths, str(out),
                                       "--not-produced", reasons, "--predictions-root", str(Path(tmp) / "preds")],
                                      capture_output=True, text=True, env=run_env)
                self.assertEqual((done.returncode, done.stderr), (0, ""), name)
                runs[name] = [(out / f).read_bytes() for f in (rp.JSON_NAME, rp.MARKDOWN_NAME)]
            self.assertEqual(runs["ascii"], runs["default"])
            self.assertIn("±".encode("utf-8"), runs["ascii"][1])
            self.assertIn("`blaTEM-1é`".encode("utf-8"), runs["ascii"][1])
            self.assertIn(b"19 files checked against the manifest", runs["ascii"][1])


def remove_tree(top):
    """shutil.rmtree, without its one recursion per level."""
    pending = [top]
    while pending:
        folder = pending[-1]
        below = [p for p in folder.iterdir() if p.is_dir() and not p.is_symlink()]
        if below:
            pending += below
            continue
        for p in folder.iterdir():
            p.unlink()
        folder.rmdir()
        pending.pop()


if __name__ == "__main__":
    unittest.main()
