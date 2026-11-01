"""Synthetic ablation cohorts with planted truths, for canaries (numpy only).

generate(config, seed, workers=1) builds genomes in lineages with a binary
phenotype y (1 = resistant, R; 0 = susceptible, S), the feature blocks the
ablation compares, and a truth manifest naming every planted effect. No real
data is read. to_json() gives the canonical bytes (contract section 7): the
same config and seed give the same bytes whatever the worker count, dict
order or PYTHONHASHSEED. The contract is
/home/layth/review/grm-ablation-contract-20261030/CONTRACT.md v1.7: block
specs and the missingness audit follow sections 5, 6 and 8, rate_diff_p and
the column format section 11, the single-class refusal's scope section 13,
and column names section 14.

usage: python3 -I -B scripts/ablation/synthetic.py --seed 7 [--config overrides.json]
           [--workers 4] [--out cohort.json | --sha256]

Cohort (a dict; every value is JSON-native):

    format, format_version   "grm-ablation-synthetic", 1
    seed                     the integer seed
    config                   the resolved config (DEFAULT_CONFIG plus overrides)
    samples                  [{sample_id, lineage, duplicate, y}] in sample order;
                             duplicate is a group id "D001" or null
    blocks                   BlockSpec records in BLOCK_ORDER, exactly
                             {id, kind, version, source_sha256, priority,
                             screen_cap, catalog_date, circularity_risk,
                             unavailable_reason}; below
    columns                  [{name, block, raw}], grouped by block in block order
    missingness              audit records, below
    truth                    {effects, disabled, null_columns}, below

Column names follow H's check_block rule (section 14): <block>:<name>, with
/ between levels inside the name (C:causal, H:fam001/allele), and no
backtick or whitespace. So integration can run H's check_block on synthetic
blocks; only check_spec's registry check is skipped (section 11).

Column format. raw is an interchange format (contract section 11), not
a serde image of H's types. H's RawColumn holds Bits and u16 codes and has no
serde derives, so the JSON does not deserialize into it. Integration (Q)
writes a test-side adapter that builds each column through H's public
constructors (RawColumn::binary, ::count, ::categorical) and refuses
out-of-range codes. One key names the kind:

    {"Binary": {"present": [bool], "assessable": [bool]}}
    {"Count": {"values": [int in [0, 2**32 - 1] | null]}}
    {"Categorical": {"codes": [int in [0, 65535] | null], "levels": [str]}}

Binary columns are JSON arrays of booleans. Counts are integers in
[0, 2**32 - 1] (H's u32), and validate_cohort() refuses a larger one.
Categorical codes are integers in [0, 65535] (H's u16), and generate() and
validate_cohort() refuse a larger code with ValueError. Every vector has one entry per sample, in sample order.
A not-assessable binary cell is present=false, assessable=false; a
not-assessable count or categorical cell is null. Codes index levels.

A cohort whose samples all have the same y is refused with ValueError
before any record is built (sections 11 and 13), and validate_cohort()
refuses one too.

Blocks (id, H's BlockKind, priority):

    C      CuratedExact  0  exact calls: C:causal (portable signal) and C:call### nulls
    H      Hierarchical  1  per family, a Binary H:fam### and its Categorical
                            H:fam###/allele (levels "absent", "a1".."aK"); they
                            share one assessable mask, and the binary is present
                            iff allele > 0
    P      Presence      2  gene presence/absence P:gene###
    N      Count         3  copy numbers N:copy### (Poisson)
    X      Context       4  always unavailable, no columns
    K      Anonymous     5  anonymous patterns K:pat####: lineage markers (group
                            proxy), then background patterns; always assessable
    noise  Anonymous     5  optional; noise:bin### and noise:cnt###, independent
                            of y, always assessable

Every BlockSpec has catalog_date null and circularity_risk false; screen_cap
is 4000 for K and null otherwise. An available block has version
"synthetic-v1" and source_sha256 = sha256(to_json(its column objects, in
column order)), the same objects that appear in columns. X has version "",
source_sha256 "" and unavailable_reason X_UNAVAILABLE.

C, H, P and N lose base_missing_rate of assessable units (a column, or an H
family) equally in both classes, unless the missingness shortcut targets the
block.

Missingness audit (contract section 6): one record per block among C, H, P,
N with base or planted missingness, in block order:

    block                     block id
    rate_R, rate_S            cell-level and descriptive: the not-assessable
                              fraction of the block's assessability cells
                              (sample x column, an H family counted once)
                              among y=1 and y=0 samples
    rate_diff_p               section 11: two-sided Fisher exact test on the
                              unit-level 2x2 table, class by whether the unit
                              has at least one missing cell in the block. A
                              duplicate group is one unit per group and class;
                              every other sample is its own unit. Two-sided
                              sums the tables no more likely than the observed
                              one (relative tolerance 1e-7); 6 significant
                              digits. Null when a class has no units, or when
                              no unit or every unit has a missing cell
    indicator_auroc           raw AUROC of the per-sample not-assessable unit
                              count predicting y = 1 (R), over all samples
    indicator_auroc_oriented  max(indicator_auroc, 1 - indicator_auroc)
    flag                      oriented > 0.60, or rate_diff_p < 0.01 (section 8)

A block with only base missingness (C, H and N at the defaults) is a null,
but it can still be flagged. Measured false flags over seeds 0-299 per row,
base-rate records only (3 per cohort). "records" is the share of base-rate
records flagged, split by which arm of the section 8 rule fired; "cohorts" is
the share of cohorts with at least one false-flagged block. The harness is
.tmp/taskl_fix1_flags.py of fix round 1.

    condition                               records  by p    by AUROC  cohorts
    defaults (n=400, duplicates 0.1)        0.89%    0.89%   0.00%     2.67%
    duplicates fraction 0                   0.67%    0.67%   0.00%     2.00%
    duplicates fraction 0.5                 0.44%    0.44%   0.00%     1.33%
    duplicates fraction 1.0                 0.44%    0.44%   0.00%     1.33%
    small class (~7% R), base rate 0.02     0.78%    0.22%   0.78%     2.33%
    small class (~7% R), base 0.002,
        duplicates off                      0.56%    0.56%   0.00%     1.67%
    n_samples 40, n_lineages 4              6.67%    0.00%   6.67%     19.00%
    n_samples 100, n_lineages 4, dup off    0.44%    0.44%   0.11%     1.33%

The p arm (Fisher exact on units) stays below 1% per record in every row,
with duplicates counted once per group. The AUROC arm (> 0.60) fires on noise
when n or a class is small: at n_samples 40 one cohort in five carries a
false flag. Canary consumers should pick n and seeds with this in mind, and
expect a false flag in about 1-3% of n=400 cohorts.

Scenarios (each has enabled in config; truth.effects carries one record per
enabled scenario with its columns, parameters, realized statistics and the
canaries below):

    portable_signal       C:causal ~ Bernoulli(frequency), independent of
                          lineage; logit P(y=1) gets +logit or -logit from it,
                          the same in every lineage.
        J: a C-adding contrast on the lineage split has a positive estimate
           with a CI above 0.
        K: C:causal may reach association_stable; the contrast claim_level
           stays predictive_increment.
        integration: C-only AUROC under lineage holdout stays within 0.05
           of its random-split AUROC when both are scored the same way, as
           the mean AUROC within each lineage. The generator alone spreads
           this gap at the defaults: over seeds 0-199, |gap| has mean
           0.013, 95th percentile 0.029 and maximum 0.0399 (seed 71), and
           over seeds 100-129 a maximum of 0.0375 (seed 110), so 0.05
           leaves headroom (.tmp/taskl_fix2_canary.py of fix round 2).
           Pooled leave-one-lineage-out
           scores run about 0.07 lower here from a fold-intercept artifact
           (holding out a high-prevalence lineage lowers its intercept), so
           read this canary on matched estimands.
    group_proxy           even-indexed lineages (L1, L3, ...) get
                          +lineage_logit, odd ones -lineage_logit, and K holds
                          markers_per_lineage patterns present in exactly one
                          lineage (each cell flipped with marker_flip). K
                          predicts y on a random split and not under lineage
                          holdout.
        J: a K contrast is large on the random split and near 0 on the
           lineage split; the splits are reported apart, never pooled.
        K: group_proxy_audit flags block K. realized.lineage_auroc (the
           planted parity offset as a score) and lineage_prevalence_auroc
           (each sample scored by its lineage's in-sample prevalence, the
           stronger) are reference scores for a pure lineage proxy, not
           ceilings: a fitted model can land above either.
           marker_lineage_auroc_mean and _min, and marker_agreement, measure
           whether the markers still tag their lineage (about 1 -
           marker_flip; 0.5 means the proxy is pure noise).
        integration: lineage folds hold each lineage whole; a leak shows as K
           keeping its random-split AUROC under lineage holdout.
    noise                 the noise block, drawn independently of y and lineage.
        J: adding it decides equivalent or inconclusive, not improves, at
           about the nominal false-positive rate.
        K: no noise column reaches association_stable.
        integration: noise-only runs score near AUROC 0.5 on both splits.
    missingness_shortcut  in one block (C, H, P or N) each unit is
                          not-assessable with rate_R for y=1 and rate_S for
                          y=0, so the not-assessable pattern alone predicts y.
        J: a missingness-indicator arm beats baseline though no call carries
           signal.
        K: the missingness record for the block is flagged, with
           indicator_auroc_oriented > 0.6; base-rate blocks are flagged
           only at the false-flag rates in the table above.
           realized is that block's audit record without its block field.
        integration: encoders (Tasks H, I) keep not-assessable apart from
           absent; collapsing them turns the shortcut into a fake call signal.
    duplicates            floor(fraction * n_samples / 2) pairs, each two
                          samples of one lineage with identical rows in every
                          column and identical y (k-mer-identical genomes).
        J: lineage-cluster resampling keeps each pair in one unit.
        K: no record of its own; prediction rows let the report check that
           a group's members share a fold.
        integration: the random split keeps each group in one outer test fold
           (docs/duplicate-isolation.md). truth groups use the
           kset_duplicates shape: ascending sample indices per group, ordered
           by first member.

Determinism: each column's values, parameter and assessable mask come from
their own PCG64 stream, seeded from (seed, sha256 of a stream name), so the
cohort is a pure function of config and seed. Stream names keep the dotted
form (C.causal, H.fam001.allele); display_name() maps one to its column
name, so the section 14 renaming left every value unchanged. Turning one scenario off leaves
unrelated columns byte-identical. Only Generator.random() uniforms are drawn,
and they are compared against probabilities computed in Python floats, so
the bytes do not depend on numpy's distribution samplers. Workers build
blocks in parallel; results are assembled in BLOCK_ORDER.
"""
import argparse
import concurrent.futures
import copy
import hashlib
import json
import math
import sys

import numpy as np

FORMAT = "grm-ablation-synthetic"
FORMAT_VERSION = 1

BLOCK_ORDER = ("C", "H", "P", "N", "X", "K", "noise")
# Kinds are H's Rust BlockKind names; priority and screen_cap follow contract
# sections 5 and 8.
BLOCK_KINDS = {"C": "CuratedExact", "H": "Hierarchical", "P": "Presence", "N": "Count",
               "X": "Context", "K": "Anonymous", "noise": "Anonymous"}
KIND_PRIORITY = {"CuratedExact": 0, "Hierarchical": 1, "Presence": 2, "Count": 3,
                 "Context": 4, "Anonymous": 5}
SCREEN_CAP = {"K": 4000}
BLOCK_VERSION = "synthetic-v1"
# Categorical codes are H's u16 and counts its u32 (contract section 11).
MAX_CODE = 65535
MAX_COUNT = 2 ** 32 - 1
BLOCK_SPEC_KEYS = ("id", "kind", "version", "source_sha256", "priority", "screen_cap",
                   "catalog_date", "circularity_risk", "unavailable_reason")
MISSINGNESS_KEYS = ("block", "rate_R", "rate_S", "rate_diff_p", "indicator_auroc",
                    "indicator_auroc_oriented", "flag")
# The missingness flag (contract section 8).
MISSINGNESS_AUROC_FLAG = 0.60
MISSINGNESS_P_FLAG = 0.01
SCENARIOS = ("portable_signal", "group_proxy", "noise", "missingness_shortcut", "duplicates")
MISSINGNESS_BLOCKS = ("C", "H", "P", "N")
X_UNAVAILABLE = "the synthetic generator plants no context features"
# P(allele == "absent") for every H family.
H_ABSENT = 0.4

DEFAULT_CONFIG = {
    "n_samples": 400,
    "n_lineages": 8,
    "base_missing_rate": 0.02,
    "blocks": {
        "C": {"n_null": 6},
        "H": {"n_families": 4, "n_alleles": 3},
        "P": {"n_columns": 20},
        "N": {"n_columns": 8},
        "K": {"n_background": 24},
    },
    "portable_signal": {"enabled": True, "frequency": 0.5, "logit": 2.0},
    "group_proxy": {"enabled": True, "lineage_logit": 2.0, "markers_per_lineage": 4,
                    "marker_flip": 0.02},
    "noise": {"enabled": True, "n_binary": 20, "n_count": 5},
    "missingness_shortcut": {"enabled": True, "block": "P", "rate_R": 0.05, "rate_S": 0.35},
    "duplicates": {"enabled": True, "fraction": 0.1},
}

CANARIES = {
    "portable_signal": [
        "J: a C-adding contrast on the lineage split has a positive estimate with a CI above 0",
        "K: C:causal may reach association_stable; the contrast claim_level stays predictive_increment",
        "integration: C-only within-lineage AUROC under lineage holdout stays within 0.05 of the random split's",
    ],
    "group_proxy": [
        "J: a K contrast is large on the random split and near 0 on the lineage split; report the splits apart",
        "K: group_proxy_audit flags block K; realized.lineage_auroc and lineage_prevalence_auroc are "
        "reference scores, not ceilings; marker_lineage_auroc_mean shows whether the markers still tag lineage",
        "integration: lineage folds hold each lineage whole; a leak shows as K keeping its random-split AUROC",
    ],
    "noise": [
        "J: adding the noise block decides equivalent or inconclusive, not improves, at about the nominal rate",
        "K: no noise column reaches association_stable",
        "integration: noise-only runs score near AUROC 0.5 on both splits",
    ],
    "missingness_shortcut": [
        "J: a missingness-indicator arm beats baseline though no call carries signal",
        "K: the missingness record for this block is flagged with indicator_auroc_oriented > 0.6; "
        "base-rate blocks only at the documented false-flag rates",
        "integration: encoders (Tasks H, I) keep not-assessable apart from absent",
    ],
    "duplicates": [
        "J: lineage-cluster resampling keeps each pair in one unit",
        "K: no record of its own; prediction rows let the report check that a group's members share a fold",
        "integration: the random split keeps each group in one outer test fold (docs/duplicate-isolation.md)",
    ],
}


# Config: overrides merge onto DEFAULT_CONFIG; unknown keys and bad values fail.

def _merge(base, over, path):
    if not isinstance(over, dict) or any(not isinstance(k, str) for k in over):
        raise ValueError(f"{path}: must be an object with string keys")
    for key in sorted(over):
        if key not in base:
            raise ValueError(f"{path}.{key}: is not a field")
        if isinstance(base[key], dict):
            _merge(base[key], over[key], f"{path}.{key}")
        else:
            base[key] = over[key]


def _int(v, path, lo, hi):
    if type(v) is not int or not lo <= v <= hi:
        raise ValueError(f"{path}: must be an integer in [{lo}, {hi}]")
    return v


def _num(v, path, lo, hi, open_ends=False):
    ok = type(v) in (int, float) and math.isfinite(v)
    ok = ok and ((lo < v < hi) if open_ends else (lo <= v <= hi))
    if not ok:
        raise ValueError(f"{path}: must be a number in {'(' if open_ends else '['}{lo}, {hi}{')' if open_ends else ']'}")
    return float(v)


def _bool(v, path):
    if type(v) is not bool:
        raise ValueError(f"{path}: must be true or false")
    return v


def resolve_config(config=None):
    """DEFAULT_CONFIG with config's overrides, checked; numbers documented as
    rates or logits come back as floats so 2 and 2.0 give the same bytes."""
    c = copy.deepcopy(DEFAULT_CONFIG)
    _merge(c, {} if config is None else config, "config")
    c["n_samples"] = _int(c["n_samples"], "config.n_samples", 8, 1_000_000)
    c["n_lineages"] = _int(c["n_lineages"], "config.n_lineages", 2, 10_000)
    if c["n_samples"] < 2 * c["n_lineages"]:
        raise ValueError("config.n_samples: must be at least 2 * n_lineages")
    c["base_missing_rate"] = _num(c["base_missing_rate"], "config.base_missing_rate", 0, 0.5)
    b = c["blocks"]
    b["C"]["n_null"] = _int(b["C"]["n_null"], "config.blocks.C.n_null", 0, 10_000)
    b["H"]["n_families"] = _int(b["H"]["n_families"], "config.blocks.H.n_families", 0, 10_000)
    b["H"]["n_alleles"] = _int(b["H"]["n_alleles"], "config.blocks.H.n_alleles", 1, 100)
    b["P"]["n_columns"] = _int(b["P"]["n_columns"], "config.blocks.P.n_columns", 0, 10_000)
    b["N"]["n_columns"] = _int(b["N"]["n_columns"], "config.blocks.N.n_columns", 0, 10_000)
    b["K"]["n_background"] = _int(b["K"]["n_background"], "config.blocks.K.n_background", 0, 100_000)
    for s in SCENARIOS:
        c[s]["enabled"] = _bool(c[s]["enabled"], f"config.{s}.enabled")
    p = c["portable_signal"]
    p["frequency"] = _num(p["frequency"], "config.portable_signal.frequency", 0, 1, open_ends=True)
    p["logit"] = _num(p["logit"], "config.portable_signal.logit", 0, 20)
    g = c["group_proxy"]
    g["lineage_logit"] = _num(g["lineage_logit"], "config.group_proxy.lineage_logit", 0, 20)
    g["markers_per_lineage"] = _int(g["markers_per_lineage"], "config.group_proxy.markers_per_lineage", 1, 1000)
    g["marker_flip"] = _num(g["marker_flip"], "config.group_proxy.marker_flip", 0, 0.5)
    z = c["noise"]
    z["n_binary"] = _int(z["n_binary"], "config.noise.n_binary", 0, 10_000)
    z["n_count"] = _int(z["n_count"], "config.noise.n_count", 0, 10_000)
    if z["enabled"] and z["n_binary"] + z["n_count"] == 0:
        raise ValueError("config.noise: an enabled noise block needs at least one column")
    m = c["missingness_shortcut"]
    if m["block"] not in MISSINGNESS_BLOCKS:
        raise ValueError(f"config.missingness_shortcut.block: must be one of {', '.join(MISSINGNESS_BLOCKS)}")
    m["rate_R"] = _num(m["rate_R"], "config.missingness_shortcut.rate_R", 0, 1)
    m["rate_S"] = _num(m["rate_S"], "config.missingness_shortcut.rate_S", 0, 1)
    if m["enabled"] and _n_units(c, m["block"]) == 0:
        raise ValueError(f"config.missingness_shortcut.block: block {m['block']} has no columns")
    d = c["duplicates"]
    d["fraction"] = _num(d["fraction"], "config.duplicates.fraction", 0, 1)
    return c


def _n_units(cfg, block):
    b = cfg["blocks"]
    return {"C": b["C"]["n_null"] + cfg["portable_signal"]["enabled"], "H": b["H"]["n_families"],
            "P": b["P"]["n_columns"], "N": b["N"]["n_columns"]}[block]


# Random streams: one per name, so no draw depends on another column's.

def _stream(seed, name):
    key = int.from_bytes(hashlib.sha256(name.encode("utf-8")).digest()[:8], "big")
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence([seed, key])))


def display_name(stream_name):
    """The section 14 column name for a dotted stream name: the first "."
    becomes ":" and any later one "/" (H.fam001.allele -> H:fam001/allele)."""
    return stream_name.replace(".", ":", 1).replace(".", "/")


def _name_ok(name, block):
    """Section 14: <block>:<name>, non-empty, with no backtick or whitespace."""
    return (isinstance(name, str) and isinstance(block, str) and name.startswith(block + ":")
            and len(name) > len(block) + 1 and not any(ch == "`" or ch.isspace() for ch in name))


def _sigmoid(x):
    return 1.0 / (1.0 + math.exp(-x))


def _poisson_cdf(lam):
    """Cumulative Poisson(lam) probabilities up to 1 - 1e-12, in Python floats."""
    cdf, p, k = [], math.exp(-lam), 0
    total = p
    while total < 1.0 - 1e-12 and k < 1000:
        cdf.append(total)
        k += 1
        p *= lam / k
        total += p
    return np.array(cdf)


class _Column:
    """kind is Binary, Count or Categorical; data holds present flags, counts
    or codes; missing marks not-assessable cells; role and tag feed truth."""
    __slots__ = ("name", "block", "kind", "data", "missing", "levels", "role", "tag", "unit")

    def __init__(self, name, block, kind, data, missing, levels=None, role="null", tag=None, unit=True):
        self.name, self.block, self.kind = name, block, kind
        self.data, self.missing, self.levels = data, missing, levels
        # unit is False for a column that shares another column's assessable
        # mask (H alleles), so the missingness audit counts each unit once.
        self.role, self.tag, self.unit = role, tag, unit

    def raw(self):
        miss = self.missing.tolist()
        if self.kind == "Binary":
            present = (self.data & ~self.missing).tolist()
            return {"Binary": {"present": present, "assessable": [not m for m in miss]}}
        values = [None if m else int(v) for v, m in zip(self.data.tolist(), miss)]
        if self.kind == "Count":
            return {"Count": {"values": values}}
        if any(v is not None and v > MAX_CODE for v in values):
            raise ValueError(f"{self.name}: categorical codes must be in [0, {MAX_CODE}] (H's u16)")
        return {"Categorical": {"codes": values, "levels": list(self.levels)}}


class _Context:
    def __init__(self, cfg, seed, lineage, c, y):
        self.cfg, self.seed, self.lineage, self.c, self.y = cfg, seed, lineage, c, y
        self.n = len(y)

    def uniform(self, name):
        return _stream(self.seed, "values:" + name).random(self.n)

    def param(self, name):
        return float(_stream(self.seed, "param:" + name).random())

    def missing(self, block, unit):
        u = _stream(self.seed, "missing:" + unit).random(self.n)
        m = self.cfg["missingness_shortcut"]
        if m["enabled"] and m["block"] == block:
            return u < np.where(self.y, m["rate_R"], m["rate_S"])
        return u < self.cfg["base_missing_rate"]

    def never_missing(self):
        return np.zeros(self.n, dtype=bool)


def _phenotype(cfg, seed):
    """Balanced lineage labels, the causal call c and y."""
    n, n_lin = cfg["n_samples"], cfg["n_lineages"]
    rank = np.empty(n, dtype=np.int64)
    rank[np.argsort(_stream(seed, "lineage").random(n), kind="stable")] = np.arange(n)
    lineage = rank % n_lin
    p, g = cfg["portable_signal"], cfg["group_proxy"]
    c = _stream(seed, "portable").random(n) < p["frequency"]
    offset = (np.where(lineage % 2 == 0, g["lineage_logit"], -g["lineage_logit"])
              if g["enabled"] else np.zeros(n))
    effect = np.where(c, p["logit"], -p["logit"]) if p["enabled"] else np.zeros(n)
    prob = np.array([_sigmoid(float(o) + float(e)) for o, e in zip(offset.tolist(), effect.tolist())])
    y = _stream(seed, "y").random(n) < prob
    return lineage, c, y


def _duplicate_pairs(cfg, seed, lineage):
    d = cfg["duplicates"]
    if not d["enabled"]:
        return []
    n = len(lineage)
    want = math.floor(d["fraction"] * n / 2)
    waiting, pairs = {}, []
    for i in np.argsort(_stream(seed, "duplicates").random(n), kind="stable").tolist():
        if len(pairs) == want:
            break
        lin = int(lineage[i])
        if lin in waiting:
            j = waiting.pop(lin)
            pairs.append((min(i, j), max(i, j)))
        else:
            waiting[lin] = i
    if len(pairs) < want:
        raise ValueError(f"config.duplicates.fraction: only {len(pairs)} same-lineage pairs fit, {want} asked")
    return sorted(pairs)


# Block builders: each returns its columns in a fixed order.

def _block_C(ctx):
    cols = []
    if ctx.cfg["portable_signal"]["enabled"]:
        cols.append(_Column("C.causal", "C", "Binary", ctx.c.copy(), ctx.missing("C", "C.causal"),
                            role="portable"))
    for j in range(ctx.cfg["blocks"]["C"]["n_null"]):
        name = f"C.call{j + 1:03d}"
        freq = 0.05 + 0.25 * ctx.param(name)
        cols.append(_Column(name, "C", "Binary", ctx.uniform(name) < freq, ctx.missing("C", name)))
    return cols


def _block_H(ctx):
    h, cols = ctx.cfg["blocks"]["H"], []
    levels = ["absent"] + [f"a{k + 1}" for k in range(h["n_alleles"])]
    for f in range(h["n_families"]):
        fam = f"H.fam{f + 1:03d}"
        weights = [w + 0.2 for w in _stream(ctx.seed, "param:" + fam).random(h["n_alleles"]).tolist()]
        total, acc, cdf = sum(weights), H_ABSENT, [H_ABSENT]
        for w in weights[:-1]:
            acc += (1.0 - H_ABSENT) * w / total
            cdf.append(acc)
        codes = np.searchsorted(np.array(cdf), ctx.uniform(fam), side="right")
        miss = ctx.missing("H", fam)
        cols.append(_Column(fam, "H", "Binary", codes > 0, miss))
        cols.append(_Column(fam + ".allele", "H", "Categorical", codes, miss.copy(), levels, unit=False))
    return cols


def _block_P(ctx):
    cols = []
    for j in range(ctx.cfg["blocks"]["P"]["n_columns"]):
        name = f"P.gene{j + 1:03d}"
        freq = 0.05 + 0.9 * ctx.param(name)
        cols.append(_Column(name, "P", "Binary", ctx.uniform(name) < freq, ctx.missing("P", name)))
    return cols


def _counts(ctx, name):
    return np.searchsorted(_poisson_cdf(0.5 + 3.5 * ctx.param(name)), ctx.uniform(name), side="right")


def _block_N(ctx):
    cols = []
    for j in range(ctx.cfg["blocks"]["N"]["n_columns"]):
        name = f"N.copy{j + 1:03d}"
        cols.append(_Column(name, "N", "Count", _counts(ctx, name), ctx.missing("N", name)))
    return cols


def _block_K(ctx):
    """Display names are anonymous and sequential; stream names are stable
    so turning the proxy off leaves the background patterns unchanged."""
    g, cols = ctx.cfg["group_proxy"], []
    if g["enabled"]:
        for lin in range(ctx.cfg["n_lineages"]):
            for m in range(g["markers_per_lineage"]):
                flip = ctx.uniform(f"K.marker.{lin + 1}.{m + 1}") < g["marker_flip"]
                cols.append(_Column(None, "K", "Binary", (ctx.lineage == lin) ^ flip, ctx.never_missing(),
                                    role="proxy", tag=lin))
    for j in range(ctx.cfg["blocks"]["K"]["n_background"]):
        key = f"K.background.{j + 1}"
        freq = 0.1 + 0.8 * ctx.param(key)
        cols.append(_Column(None, "K", "Binary", ctx.uniform(key) < freq, ctx.never_missing()))
    for k, col in enumerate(cols):
        col.name = f"K.pat{k + 1:04d}"
    return cols


def _block_noise(ctx):
    z, cols = ctx.cfg["noise"], []
    for j in range(z["n_binary"]):
        name = f"noise.bin{j + 1:03d}"
        freq = 0.1 + 0.8 * ctx.param(name)
        cols.append(_Column(name, "noise", "Binary", ctx.uniform(name) < freq, ctx.never_missing(),
                            role="noise"))
    for j in range(z["n_count"]):
        name = f"noise.cnt{j + 1:03d}"
        cols.append(_Column(name, "noise", "Count", _counts(ctx, name), ctx.never_missing(), role="noise"))
    return cols


BUILDERS = {"C": _block_C, "H": _block_H, "P": _block_P, "N": _block_N, "K": _block_K,
            "noise": _block_noise}


def _auroc(score, y):
    """Mann-Whitney AUROC with average ranks for ties, rounded to 6 places;
    None when a class is empty."""
    score, y = np.asarray(score, dtype=float), np.asarray(y, dtype=bool)
    n1 = int(y.sum())
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return None
    _, inv, cnt = np.unique(score, return_inverse=True, return_counts=True)
    ranks = (np.cumsum(cnt) - (cnt - 1) / 2.0)[inv]
    return round(float((ranks[y].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)), 6)


def _rate(x):
    return round(float(np.mean(x)), 6) if len(x) else None


def fisher_exact_p(a, b, c, d):
    """Two-sided Fisher exact p for the 2x2 table [[a, b], [c, d]]: the sum of
    the probabilities of all tables with the same margins that are no more
    likely than the observed one, with relative tolerance 1e-7. Probabilities
    are built relative to the modal table by the hypergeometric ratio, using
    only IEEE +, *, / (no libm), so the value is the same on every platform."""
    n1, n0, k = a + b, c + d, a + c
    lo, hi = max(0, k - n0), min(k, n1)

    def ratio(x):  # P(x + 1) / P(x)
        return ((k - x) * (n1 - x)) / ((x + 1) * (n0 - k + x + 1))

    mode = min(max((k + 1) * (n1 + 1) // (n1 + n0 + 2), lo), hi)
    rel = {mode: 1.0}
    for x in range(mode, hi):
        rel[x + 1] = rel[x] * ratio(x)
    for x in range(mode, lo, -1):
        rel[x - 1] = rel[x] / ratio(x - 1)
    cut = rel[a] * (1.0 + 1e-7)
    total = small = 0.0
    for x in range(lo, hi + 1):
        total += rel[x]
        if rel[x] <= cut:
            small += rel[x]
    return min(1.0, small / total)


def _rate_diff_p(miss_any, y, groups):
    """Contract section 11 rate_diff_p: Fisher exact on units (a duplicate
    group per class, or a lone sample) by whether any cell is missing; 6
    significant digits; None when a class has no units or no unit, or every
    unit, has a missing cell."""
    missing = {}  # (class, unit id) -> any missing cell
    for i, (cls, m) in enumerate(zip(y, miss_any)):
        key = (bool(cls), i if groups[i] is None else groups[i])
        missing[key] = missing.get(key, False) or bool(m)
    table = [[0, 0], [0, 0]]  # rows R, S; columns missing, complete
    for (cls, _), m in missing.items():
        table[0 if cls else 1][0 if m else 1] += 1
    (a, b), (c, d) = table
    if a + b == 0 or c + d == 0 or a + c == 0 or b + d == 0:
        return None
    return float(f"{fisher_exact_p(a, b, c, d):.6g}")


def missingness_flag(oriented, p):
    """The section 8 flag: oriented AUROC strictly above 0.60, or p below 0.01."""
    return bool((oriented is not None and oriented > MISSINGNESS_AUROC_FLAG)
                or (p is not None and p < MISSINGNESS_P_FLAG))


def _missingness_audit(cfg, columns, y, groups):
    """One contract section 6 record per block with base or planted
    missingness, in BLOCK_ORDER. Cells are assessability masks: a column, or
    an H family (its binary and allele columns share one mask)."""
    shortcut = cfg["missingness_shortcut"]
    out = []
    for b in MISSINGNESS_BLOCKS:
        units = [c.missing for c in columns if c.block == b and c.unit]
        planted = shortcut["enabled"] and shortcut["block"] == b
        if not units or not (planted or cfg["base_missing_rate"] > 0):
            continue
        miss = np.stack(units, axis=1)
        raw = _auroc(miss.sum(axis=1), y)
        oriented = None if raw is None else round(max(raw, 1.0 - raw), 6)
        p = _rate_diff_p(miss.any(axis=1).tolist(), y.tolist(), groups)
        out.append({
            "block": b, "rate_R": _rate(miss[y]), "rate_S": _rate(miss[~y]), "rate_diff_p": p,
            "indicator_auroc": raw, "indicator_auroc_oriented": oriented,
            "flag": missingness_flag(oriented, p),
        })
    return out


def _block_spec(block, column_json):
    """Contract section 5 BlockSpec; source_sha256 hashes the canonical JSON
    of the block's column objects, in column order."""
    available = block != "X"
    return {
        "id": block, "kind": BLOCK_KINDS[block],
        "version": BLOCK_VERSION if available else "",
        "source_sha256": hashlib.sha256(to_json(column_json).encode("utf-8")).hexdigest() if available else "",
        "priority": KIND_PRIORITY[BLOCK_KINDS[block]], "screen_cap": SCREEN_CAP.get(block),
        "catalog_date": None, "circularity_risk": False,
        "unavailable_reason": None if available else X_UNAVAILABLE,
    }


def _truth(cfg, columns, lineage_ids, lineage, y, pairs, missingness):
    effects = []
    by_role = {}
    for col in columns:
        by_role.setdefault(col.role, []).append(col)
    if cfg["portable_signal"]["enabled"]:
        col = by_role["portable"][0]
        ok = ~col.missing
        effects.append({
            "scenario": "portable_signal", "block": "C", "columns": [col.name],
            "params": dict(cfg["portable_signal"]),
            "realized": {
                "auroc": _auroc(col.data[ok], y[ok]),
                "auroc_by_lineage": {lid: _auroc(col.data[ok & (lineage == k)], y[ok & (lineage == k)])
                                     for k, lid in enumerate(lineage_ids)},
            },
            "canaries": CANARIES["portable_signal"],
        })
    if cfg["group_proxy"]["enabled"]:
        g = cfg["group_proxy"]
        markers = {lid: [] for lid in lineage_ids}
        for col in by_role["proxy"]:
            markers[lineage_ids[col.tag]].append(col.name)
        prevalence = np.array([float(np.mean(y[lineage == k])) for k in range(len(lineage_ids))])
        # How well each marker still tags its own lineage: about 1 - marker_flip.
        marker_auroc = [_auroc(col.data, lineage == col.tag) for col in by_role["proxy"]]
        effects.append({
            "scenario": "group_proxy", "block": "K", "columns": [c.name for c in by_role["proxy"]],
            "params": dict(g),
            "lineage_logit_offset": {lid: (g["lineage_logit"] if k % 2 == 0 else -g["lineage_logit"])
                                     for k, lid in enumerate(lineage_ids)},
            "markers_by_lineage": markers,
            "realized": {
                "prevalence_by_lineage": {lid: _rate(y[lineage == k]) for k, lid in enumerate(lineage_ids)},
                # Reference scores for a pure lineage proxy, not ceilings: the
                # planted parity offset, and each sample's lineage prevalence
                # (in-sample, so the stronger of the two).
                "lineage_auroc": _auroc(lineage % 2 == 0, y),
                "lineage_prevalence_auroc": _auroc(prevalence[lineage], y),
                "marker_lineage_auroc_mean": round(sum(marker_auroc) / len(marker_auroc), 6),
                "marker_lineage_auroc_min": min(marker_auroc),
                "marker_agreement": round(float(np.mean([np.mean(col.data == (lineage == col.tag))
                                                         for col in by_role["proxy"]])), 6),
            },
            "canaries": CANARIES["group_proxy"],
        })
    if cfg["noise"]["enabled"]:
        cols = by_role["noise"]
        effects.append({
            "scenario": "noise", "block": "noise", "columns": [c.name for c in cols],
            "params": dict(cfg["noise"]),
            "realized": {"max_abs_auroc_minus_half": round(max(abs(_auroc(c.data, y) - 0.5) for c in cols), 6)},
            "canaries": CANARIES["noise"],
        })
    if cfg["missingness_shortcut"]["enabled"]:
        m = cfg["missingness_shortcut"]
        audit = next(r for r in missingness if r["block"] == m["block"])
        effects.append({
            "scenario": "missingness_shortcut", "block": m["block"],
            "columns": [c.name for c in columns if c.block == m["block"]],
            "params": dict(m),
            "realized": {k: v for k, v in audit.items() if k != "block"},
            "canaries": CANARIES["missingness_shortcut"],
        })
    if cfg["duplicates"]["enabled"]:
        effects.append({
            "scenario": "duplicates", "block": None, "columns": [],
            "params": dict(cfg["duplicates"]),
            "groups": [list(p) for p in pairs],
            "realized": {"fraction": round(2 * len(pairs) / len(y), 6)},
            "canaries": CANARIES["duplicates"],
        })
    shortcut = cfg["missingness_shortcut"]["block"] if cfg["missingness_shortcut"]["enabled"] else None
    return {
        "effects": effects,
        "disabled": [s for s in SCENARIOS if not cfg[s]["enabled"]],
        # Columns with no planted relation to y or lineage (duplicate copying aside).
        "null_columns": [c.name for c in columns if c.role == "null" and c.block != shortcut],
    }


def generate(config=None, seed=0, workers=1):
    """The cohort dict for config (overrides of DEFAULT_CONFIG) and seed."""
    cfg = resolve_config(config)
    if type(seed) is not int or not 0 <= seed < 2 ** 63:
        raise ValueError("seed: must be an integer in [0, 2**63)")
    if type(workers) is not int or workers < 1:
        raise ValueError("workers: must be an integer >= 1")
    n, n_lin = cfg["n_samples"], cfg["n_lineages"]
    lineage, c, y = _phenotype(cfg, seed)
    pairs = _duplicate_pairs(cfg, seed, lineage)
    src = np.array([p[0] for p in pairs], dtype=np.int64)
    dst = np.array([p[1] for p in pairs], dtype=np.int64)
    # A copy takes its source's y; its columns are copied after they are built.
    y[dst] = y[src]
    if y.all() or not y.any():
        raise ValueError(f"cohort: every sample has y = {int(y[0])}; a single-class cohort is refused "
                         "(contract v1.7 sections 11 and 13)")
    ctx = _Context(cfg, seed, lineage, c, y)
    order = [b for b in BLOCK_ORDER if b in BUILDERS and (b != "noise" or cfg["noise"]["enabled"])]
    if workers == 1:
        built = [BUILDERS[b](ctx) for b in order]
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            built = list(pool.map(lambda b: BUILDERS[b](ctx), order))
    by_block = dict(zip(order, built))
    columns = [col for b in order for col in by_block[b]]
    for col in columns:
        col.name = display_name(col.name)

    for col in columns:
        col.data[dst] = col.data[src]
        col.missing[dst] = col.missing[src]

    width = len(str(n))
    lineage_ids = [f"L{k + 1:0{len(str(n_lin))}d}" for k in range(n_lin)]
    dup = {}
    for g, (i, j) in enumerate(pairs):
        dup[i] = dup[j] = f"D{g + 1:0{max(3, len(str(len(pairs))))}d}"
    samples = [{"sample_id": f"S{i + 1:0{width}d}", "lineage": lineage_ids[int(lineage[i])],
                "duplicate": dup.get(i), "y": int(y[i])} for i in range(n)]
    column_json = {b: [{"name": col.name, "block": col.block, "raw": col.raw()} for col in by_block.get(b, [])]
                   for b in BLOCK_ORDER if b in by_block or b == "X"}
    missingness = _missingness_audit(cfg, columns, y, [dup.get(i) for i in range(n)])
    return {
        "format": FORMAT, "format_version": FORMAT_VERSION, "seed": seed, "config": cfg,
        "samples": samples,
        "blocks": [_block_spec(b, column_json[b]) for b in BLOCK_ORDER if b in column_json],
        "columns": [col for b in BLOCK_ORDER if b in column_json for col in column_json[b]],
        "missingness": missingness,
        "truth": _truth(cfg, columns, lineage_ids, lineage, y, pairs, missingness),
    }


def _refuse_non_string_keys(obj, path="cohort"):
    if isinstance(obj, tuple):
        raise ValueError(f"{path}: tuples are not canonical; use a list")
    if isinstance(obj, dict):
        for k, v in obj.items():
            if not isinstance(k, str):
                raise ValueError(f"{path}: key {k!r} is not a string")
            _refuse_non_string_keys(v, f"{path}.{k}")
    elif isinstance(obj, list) and not set(map(type, obj)).isdisjoint((dict, list, tuple)):
        for i, v in enumerate(obj):
            if isinstance(v, (dict, list, tuple)):
                _refuse_non_string_keys(v, f"{path}[{i}]")


def to_json(obj):
    """Canonical JSON (contract section 7): sorted keys, no whitespace,
    ASCII, no NaN; non-string keys and tuples are refused."""
    _refuse_non_string_keys(obj)
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def sha256(obj):
    """sha256 over the UTF-8 bytes of to_json(obj)."""
    return hashlib.sha256(to_json(obj).encode("utf-8")).hexdigest()


def validate_cohort(cohort):
    """Structural checks a consumer can run before use; raises ValueError
    listing every problem."""
    errors = []
    keys = {"format", "format_version", "seed", "config", "samples", "blocks", "columns", "missingness", "truth"}
    if not isinstance(cohort, dict) or set(cohort) != keys:
        raise ValueError(f"cohort: must be an object with keys {', '.join(sorted(keys))}")
    if cohort["format"] != FORMAT or cohort["format_version"] != FORMAT_VERSION:
        errors.append(f"cohort: format must be {FORMAT} version {FORMAT_VERSION}")
    samples = cohort["samples"]
    n = len(samples)
    ids = [s.get("sample_id") for s in samples]
    if len(set(ids)) != n:
        errors.append("samples: sample_id values must be unique")
    groups = {}
    for i, s in enumerate(samples):
        if set(s) != {"sample_id", "lineage", "duplicate", "y"}:
            errors.append(f"samples[{i}]: keys must be sample_id, lineage, duplicate, y")
            continue
        if type(s["y"]) is not int or s["y"] not in (0, 1):
            errors.append(f"samples[{i}].y: must be 0 or 1")
        if s["duplicate"] is not None:
            groups.setdefault(s["duplicate"], []).append(s["lineage"])
    ys = {s.get("y") for s in samples if isinstance(s, dict)}
    if len(ys) == 1 and ys <= {0, 1}:
        errors.append(f"samples: every sample has y = {ys.pop()}; a single-class cohort is refused")
    for gid, lins in sorted(groups.items()):
        if len(lins) < 2 or len(set(lins)) != 1:
            errors.append(f"samples: duplicate group {gid} must have 2+ members in one lineage")
    names = [c.get("name") for c in cohort["columns"]]
    if len(set(names)) != len(names):
        errors.append("columns: names must be unique")
    block_ids = [b.get("id") for b in cohort["blocks"]]
    if len(set(block_ids)) != len(block_ids):
        errors.append("blocks: ids must be unique")
    order = [c.get("block") for c in cohort["columns"]]
    if order != sorted(order, key=lambda b: block_ids.index(b) if b in block_ids else -1):
        errors.append("columns: must be grouped by block, in block order")
    for b in cohort["blocks"]:
        path = f"blocks.{b.get('id')}"
        if set(b) != set(BLOCK_SPEC_KEYS):
            errors.append(f"{path}: keys must be {', '.join(BLOCK_SPEC_KEYS)}")
            continue
        if KIND_PRIORITY.get(b["kind"]) != b["priority"]:
            errors.append(f"{path}.priority: must match kind {b['kind']}")
        own = [c for c in cohort["columns"] if c.get("block") == b["id"]]
        if b["unavailable_reason"] is None:
            if b["version"] != BLOCK_VERSION or b["source_sha256"] != sha256(own):
                errors.append(f"{path}: version must be {BLOCK_VERSION} and source_sha256 must hash its columns")
        elif own or not b["unavailable_reason"]:
            errors.append(f"{path}: an unavailable block has a reason and no columns")
    for i, r in enumerate(cohort["missingness"]):
        path = f"missingness[{i}]"
        if set(r) != set(MISSINGNESS_KEYS) or r["block"] not in block_ids:
            errors.append(f"{path}: keys must be {', '.join(MISSINGNESS_KEYS)} and block a known block")
            continue
        raw, oriented, p = r["indicator_auroc"], r["indicator_auroc_oriented"], r["rate_diff_p"]
        if (raw is None) != (oriented is None) or (raw is not None and abs(oriented - max(raw, 1 - raw)) > 1e-9):
            errors.append(f"{path}.indicator_auroc_oriented: must be max(indicator_auroc, 1 - indicator_auroc)")
        expected = ((oriented is not None and oriented > MISSINGNESS_AUROC_FLAG)
                    or (p is not None and p < MISSINGNESS_P_FLAG))
        if r["flag"] is not expected:
            errors.append(f"{path}.flag: must be {str(expected).lower()}")
    available = {b.get("id") for b in cohort["blocks"] if b.get("unavailable_reason") is None}
    for c in cohort["columns"]:
        path = f"columns.{c.get('name')}"
        if c.get("block") not in available:
            errors.append(f"{path}.block: must be an available block")
        if not _name_ok(c.get("name"), c.get("block")):
            errors.append(f"{path}.name: must be <block>:<name> with no backtick or whitespace (section 14)")
        raw = c.get("raw")
        if not isinstance(raw, dict) or len(raw) != 1:
            errors.append(f"{path}.raw: must hold exactly one of Binary, Count, Categorical")
            continue
        (kind, body), = raw.items()
        if kind == "Binary":
            pr, asx = body.get("present"), body.get("assessable")
            if not (isinstance(pr, list) and isinstance(asx, list) and len(pr) == len(asx) == n):
                errors.append(f"{path}: present and assessable need {n} entries")
            elif any(type(v) is not bool for v in pr + asx) or any(p and not a for p, a in zip(pr, asx)):
                errors.append(f"{path}: entries must be booleans and present implies assessable")
        elif kind == "Count":
            vals = body.get("values")
            if not isinstance(vals, list) or len(vals) != n:
                errors.append(f"{path}.values: needs {n} entries")
            elif any(v is not None and (type(v) is not int or not 0 <= v <= MAX_COUNT) for v in vals):
                errors.append(f"{path}.values: entries must be integers in [0, {MAX_COUNT}] (H's u32) or null")
        elif kind == "Categorical":
            codes, levels = body.get("codes"), body.get("levels")
            if not isinstance(levels, list) or not levels or len(set(levels)) != len(levels) \
                    or any(not isinstance(lv, str) or not lv for lv in levels):
                errors.append(f"{path}.levels: must be distinct non-empty strings")
            elif not isinstance(codes, list) or len(codes) != n:
                errors.append(f"{path}.codes: needs {n} entries")
            elif any(v is not None and (type(v) is not int or v > MAX_CODE) for v in codes):
                errors.append(f"{path}.codes: entries must be integers in [0, {MAX_CODE}] (H's u16) or null")
            elif any(v is not None and not 0 <= v < len(levels) for v in codes):
                errors.append(f"{path}.codes: entries must index levels or be null")
        else:
            errors.append(f"{path}.raw: unknown kind {kind}")
    if errors:
        raise ValueError("; ".join(errors))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Write a synthetic ablation cohort as canonical JSON.")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--config", help="JSON file of overrides to DEFAULT_CONFIG")
    ap.add_argument("--workers", type=int, default=1)
    out = ap.add_mutually_exclusive_group()
    out.add_argument("--out", help="write the cohort here (default: stdout)")
    out.add_argument("--sha256", action="store_true", help="print only the sha256 of the canonical JSON")
    args = ap.parse_args(argv)
    config = None
    if args.config:
        with open(args.config, encoding="utf-8") as fh:
            config = json.load(fh)
    cohort = generate(config, args.seed, workers=args.workers)
    data = to_json(cohort).encode("ascii")
    if args.sha256:
        print(hashlib.sha256(data).hexdigest())
    elif args.out:
        with open(args.out, "wb") as fh:
            fh.write(data)
    else:
        sys.stdout.buffer.write(data)
    return 0


if __name__ == "__main__":
    sys.exit(main())
