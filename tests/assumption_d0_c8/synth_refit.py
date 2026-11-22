"""SYNTHETIC copies of the refit's species directories, for the tests of
scripts/assumption-d0-c8/read_refit.py and `c8.py combine`.

Three made-up species, six datasets that are not M and one that is M,
written as scripts/refit/run_species.sh writes its copy: DONE, SHA256SUMS,
binary.txt, expect.txt, d0.txt, check/ and results/. Every value is made up
here; none comes from a run.

D0 by hand (gap is the program's lineage balanced accuracy minus the
learner's):

    synthia_alpha__drugamycin   scm absent 0.20; rf removed+unused 0.15; lasso unused 0.30
    synthia_alpha__fakecillin   rf lineage 0.15
    synthia_beta__drugamycin    scm absent 0.15; cart no missed fold 0.15

So 3 D0 datasets in 2 species and 6 miss pairs: absent 2, removed 1,
lineage 1, chosen 0, unused 2; absent and unused tie for most, so no cause
sets the order. synthia_alpha__mdrug is in the record's MISS set and is M,
so the refit leaves it out.
"""

import hashlib
import os

import c8
import d0
import foldread as fr

CODE = "c0de" * 16
BINARY = "b1" * 32
SPECIES = {
    "synthia_alpha": ["synthia_alpha__drugamycin", "synthia_alpha__fakecillin",
                      "synthia_alpha__mockazole_trimethoprim"],
    "synthia_beta": ["synthia_beta__drugamycin", "synthia_beta__nalidixic_acid"],
    "synthia_gamma": ["synthia_gamma__fakecillin"],
}
M_DATASET = "synthia_alpha__mdrug"
ANTIBIOTIC = {
    "synthia_alpha__drugamycin": "drugamycin", "synthia_alpha__fakecillin": "fakecillin",
    "synthia_alpha__mockazole_trimethoprim": "mockazole/trimethoprim", "synthia_alpha__mdrug": "mdrug",
    "synthia_beta__drugamycin": "drugamycin", "synthia_beta__nalidixic_acid": "nalidixic acid",
    "synthia_gamma__fakecillin": "fakecillin",
}
PROGRAM = {"synthia_alpha__drugamycin": 0.9, "synthia_alpha__fakecillin": 0.8, "synthia_beta__drugamycin": 0.95,
           "synthia_alpha__mdrug": 0.9}
# (dataset, learner): (cause or None for no missed fold, gap). Every other
# pair of a D0 dataset is outside the miss set.
MISS = {
    ("synthia_alpha__drugamycin", "scm"): ("absent", 0.20),
    ("synthia_alpha__drugamycin", "rf"): ("removed+unused", 0.15),
    ("synthia_alpha__drugamycin", "lasso"): ("unused", 0.30),
    ("synthia_alpha__fakecillin", "rf"): ("lineage", 0.15),
    ("synthia_beta__drugamycin", "scm"): ("absent", 0.15),
    ("synthia_beta__drugamycin", "cart"): (None, 0.15),
    ("synthia_alpha__mdrug", "scm"): ("absent", 0.30),
}
D0 = {sp: sorted(d for (d, _) in MISS if d.startswith(sp + "__") and d != M_DATASET) for sp in SPECIES}
D0 = {sp: sorted(set(v)) for sp, v in D0.items()}
# The scm lineage overlap per dataset; synthia_gamma__fakecillin's one fold
# that chose gives 4 pairs at 0.
SCM_LINEAGE = {"synthia_alpha__drugamycin": 0.6, "synthia_alpha__fakecillin": 0.4,
               "synthia_alpha__mockazole_trimethoprim": 0.8, "synthia_beta__drugamycin": 0.2,
               "synthia_beta__nalidixic_acid": 0.3, "synthia_gamma__fakecillin": 0.0}


def write_record(directory):
    """reference/amrfinderplus-comparison's two tables, as miss_sets.py reads them."""
    os.makedirs(directory, exist_ok=True)
    every = [d for v in SPECIES.values() for d in v] + [M_DATASET]
    with open(os.path.join(directory, "datasets.tsv"), "w", encoding="utf-8") as f:
        f.write("dataset\tspecies\torganism_option\tlineage_split\n")
        for d in sorted(every):
            f.write(f"{d}\t{d.split('__')[0]}\tyes\tyes\n")
    with open(os.path.join(directory, "callers.tsv"), "w", encoding="utf-8") as f:
        f.write("dataset\tcaller\tbalanced_accuracy\n")
        for d in sorted(every):
            p = PROGRAM.get(d, 0.7)
            f.write(f"{d}\tprogram:tool_default\t{p:.6f}\n")
            for l in fr.LEARNERS:
                gap = MISS[(d, l)][1] if (d, l) in MISS else 0.0
                f.write(f"{d}\tlearner:{l}:lineage\t{p - gap:.6f}\n")
    return directory


def write_marks(path, extra=()):
    with open(path, "w", encoding="utf-8") as f:
        f.write("dataset\tpairs_S\tpairs_S+I\tpairs_I+R\n")
        for d in sorted([d for v in SPECIES.values() for d in v] + list(extra)):
            f.write(f"{d}\t5\t0\t0\n")
        f.write(f"{M_DATASET}\t5\t1\t0\n")
    return path


def write_names(path):
    with open(path, "w", encoding="utf-8") as f:
        f.write("dataset_id\tspecies_slug\tantibiotic\n")
        for d, ab in sorted(ANTIBIOTIC.items()):
            f.write(f"{d}\t{d.split('__')[0]}\t{ab}\n")
    return path


def c8_rows(slug):
    out = []
    for d in SPECIES[slug]:
        for l in fr.LEARNERS:
            row = [d, slug, False, l]
            for s in fr.SCHEMES:
                pattern = l in fr.PATTERN_LEARNERS
                if s == "random" and d == "synthia_beta__drugamycin":
                    row += [0, 0, None, None, None, None, None]
                elif s == "lineage" and l == "scm" and d == "synthia_gamma__fakecillin":
                    row += [5, 4, 0.0, 0.0, 1, 0, 0]
                elif s == "lineage" and l == "cart" and d == "synthia_beta__nalidixic_acid":
                    row += [5, 0, None, None, 0, 0, 0]
                else:
                    v = SCM_LINEAGE[d] if (s == "lineage" and l == "scm") else 0.5
                    row += [5, 10, v, v if pattern else None, 4, 2, 1]
            out.append(row)
    return out


def d0_rows(slug):
    rows = []
    for d in D0[slug]:
        for l in fr.LEARNERS:
            cause, gap = MISS.get((d, l), ("", 0.0))
            in_miss = (d, l) in MISS
            p = PROGRAM[d]
            tally = dict.fromkeys(d0.CAUSES, 0)
            parts = [c for c in (cause or "").split("+") if c]
            for c in parts:
                tally[c] += 1
            rows.append({"dataset": d, "species": slug, "learner": l, "source": "fold models", "program_ba": p,
                         "learner_lineage_ba": p - gap, "gap": gap, "in_miss": in_miss, "screen": 4000,
                         "missed_folds": len(parts) if in_miss else 0, "fallback": False,
                         **{f"folds_{c}": tally[c] for c in d0.CAUSES},
                         "cause": cause if in_miss and cause else None})
    return rows


def write_species(root, slug, code=CODE, name=None):
    """One species' copy, sealed; returns its directory."""
    out = os.path.join(root, name or slug)
    os.makedirs(os.path.join(out, "check"), exist_ok=True)
    os.makedirs(os.path.join(out, "results", "c8"), exist_ok=True)
    with open(os.path.join(out, "binary.txt"), "w", encoding="utf-8") as f:
        f.write(f"commit {'0' * 40}\ngrm {BINARY}\ncode {code}\n")
    with open(os.path.join(out, "expect.txt"), "w", encoding="utf-8") as f:
        f.write("".join(d + "\n" for d in SPECIES[slug]))
    with open(os.path.join(out, "d0.txt"), "w", encoding="utf-8") as f:
        f.write("".join(d + "\n" for d in D0[slug]))
    with open(os.path.join(out, "check", "datasets.tsv"), "w", encoding="utf-8") as f:
        f.write("dataset\tfitted_fields\tfields_that_differ\tsamples_with_other_descriptive_fields\tfold_models"
                "\tsamples_with_rule_groups\tunexpected_new_keys\ttext_equal_less_allowed\n")
        for d in SPECIES[slug]:
            f.write(f"{d}\tequal\t-\t40\tcomplete\t60\t-\tyes\n")
    with open(os.path.join(out, "check", "learners.tsv"), "w", encoding="utf-8") as f:
        f.write("dataset\tsplit\tlearner\tverdict\tfields_that_differ\tlargest_abs_difference\n")
        for d in SPECIES[slug]:
            for s in fr.SCHEMES:
                for l in list(fr.LEARNERS) + ["catalog"]:
                    f.write(f"{d}\t{s}\t{l}\tequal\t0\t0\n")
    with open(os.path.join(out, "check", "release-specs.tsv"), "w", encoding="utf-8") as f:
        f.write("dataset\tbreakpoint\troot_r_above\trelease_r_above\troot_n\trelease_n\n")
        for i, d in enumerate(SPECIES[slug]):
            f.write(f"{d}\t{'differs' if i == 0 else 'equal'}\t3\t3\t{101 if i == 0 else 100}\t100\n")
    with open(os.path.join(out, "check", "feature-lists.tsv"), "w", encoding="utf-8") as f:
        f.write("fetched\tequal_to_master\n11\t0\n")
    fr.write_tsv(os.path.join(out, "results", "c8", "c8-datasets.tsv"), c8.DATASET_HEADER, c8_rows(slug))
    if D0[slug]:
        dest = os.path.join(out, "results", "d0")
        d0.write_dataset(dest, [], d0_rows(slug))
    seal(out, slug, code)
    return out


def seal(out, slug=None, code=None):
    """SHA256SUMS over the files the run script lists, then DONE."""
    if slug is None or code is None:
        with open(os.path.join(out, "DONE"), encoding="utf-8") as f:
            parts = f.read().split()
        slug, code = parts[1], parts[6]
    names = []
    for base in ("check", "results"):
        for dirpath, _, files in os.walk(os.path.join(out, base)):
            for n in files:
                if n.endswith((".tsv", ".txt", ".json")):
                    names.append(os.path.relpath(os.path.join(dirpath, n), out))
    names += ["binary.txt", "expect.txt", "d0.txt"]
    lines = []
    for rel in sorted(names):
        with open(os.path.join(out, rel), "rb") as f:
            lines.append(f"{hashlib.sha256(f.read()).hexdigest()}  {rel}\n")
    with open(os.path.join(out, "SHA256SUMS"), "w", encoding="utf-8") as f:
        f.write("".join(lines))
    with open(os.path.join(out, "SHA256SUMS"), "rb") as f:
        sums = hashlib.sha256(f.read()).hexdigest()
    with open(os.path.join(out, "d0.txt"), encoding="utf-8") as f:
        n = len(f.read().split())
    d0line = f"d0 {n} datasets" if n else "d0 skipped: no MISS dataset among the listed"
    with open(os.path.join(out, "DONE"), "w", encoding="utf-8") as f:
        f.write(f"DONE {slug} 2026-11-22T00:00:00Z binary {BINARY} code {code} {d0line} SHA256SUMS {sums}\n")


def build(root):
    """Every input of one clean reading under `root`."""
    back = os.path.join(root, "back")
    dirs = [write_species(back, sp) for sp in SPECIES]
    return {
        "dirs": dirs,
        "record": write_record(os.path.join(root, "record")),
        "marks": write_marks(os.path.join(root, "marks.tsv")),
        "names": write_names(os.path.join(root, "names.tsv")),
    }


def args(paths, out, **over):
    """The reader's command line for `paths`, with options replaced by `over`
    (a value of None drops the option)."""
    opts = {"--out": out, "--code-sha256": CODE, "--species": "3", "--datasets": "6", "--d0-datasets": "3",
            "--d0-species": "2", "--d0-pairs": "6", "--marks": paths["marks"], "--record": paths["record"],
            "--names": paths["names"], "--c8-one-fold": "zero", "--c8-over-five": "pairs", "--c8-bands": "none",
            "--order-rest": "first-only"}
    dirs = over.pop("dirs", paths["dirs"])
    for k, v in over.items():
        opts["--" + k.replace("_", "-")] = v
    argv = []
    for k, v in opts.items():
        if v is not None:
            argv += [k, v]
    return argv + list(dirs)
