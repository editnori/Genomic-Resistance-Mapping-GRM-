"""tables.py ROOT [--check FILE]: the tables of docs/outside-study-results.md.

ROOT holds the outside-study kit's output files, one folder per run:
ROOT/<population>/<species>/ with `metrics.tsv`, `datasets.tsv`, `groups.tsv`
and `reading.txt` (docs/outside-study-test.md, Outputs).

The script estimates nothing. Each figure is a cell of those files, or a count
of their rows, or a cell shown as a percentage to one decimal. The readings
("below", "consistent", "above") are the kit's; the script counts them. The two
populations are never added together.

With --check FILE it prints nothing and exits 1 unless FILE holds, between the
lines `<!-- tables:begin -->` and `<!-- tables:end -->`, exactly the text this
script prints.
"""
import csv
import pathlib
import re
import sys

LEARNERS = ["scm", "cart", "lasso", "known_gene"]
NAMES = {"scm": "Set covering machine", "cart": "Decision tree", "lasso": "Lasso", "known_gene": "Known-gene rule"}
POPULATIONS = [("primary", "Primary population"), ("same-study", "Same-study population")]
READINGS = ["below", "consistent", "above"]
BEGIN, END = "<!-- tables:begin -->", "<!-- tables:end -->"
# Two drugs given together; their slug joins the two names with "_".
COMBINATIONS = {"amoxicillin_clavulanate", "ampicillin_sulbactam", "ceftazidime_avibactam", "piperacillin_tazobactam", "sulfamethoxazole_trimethoprim"}

# The rows that the document's text cites: population, species, drug, learner.
# The choice is by hand. Every figure in a printed row is the file's.
CITED = [
    ("primary", "staphylococcus_aureus", "oxacillin", "scm"),
    ("primary", "staphylococcus_aureus", "tetracycline", "known_gene"),
    ("primary", "streptococcus_pneumoniae", "chloramphenicol", "scm"),
    ("primary", "streptococcus_pneumoniae", "penicillin", "scm"),
    ("primary", "escherichia_coli", "amikacin", "scm"),
    ("primary", "escherichia_coli", "ceftazidime", "scm"),
    ("primary", "escherichia_coli", "sulfamethoxazole_trimethoprim", "scm"),
    ("primary", "mycobacterium_tuberculosis", "ethambutol", "scm"),
    ("primary", "mycobacterium_tuberculosis", "pyrazinamide", "scm"),
    ("primary", "mycobacterium_tuberculosis", "levofloxacin", "scm"),
    ("primary", "salmonella_enterica", "sulfamethoxazole", "scm"),
    ("primary", "salmonella_enterica", "sulfamethoxazole", "lasso"),
    ("primary", "salmonella_enterica", "sulfamethoxazole", "known_gene"),
    ("primary", "salmonella_enterica", "tetracycline", "scm"),
    ("primary", "salmonella_enterica", "tetracycline", "known_gene"),
    ("primary", "klebsiella_pneumoniae", "meropenem", "scm"),
    ("primary", "klebsiella_pneumoniae", "meropenem", "known_gene"),
    ("primary", "klebsiella_pneumoniae", "ertapenem", "cart"),
    ("primary", "klebsiella_pneumoniae", "tigecycline", "scm"),
    ("primary", "shigella_flexneri", "ciprofloxacin", "scm"),
    ("same-study", "klebsiella_pneumoniae", "meropenem", "scm"),
    ("same-study", "escherichia_coli", "ceftazidime", "scm"),
    ("same-study", "streptococcus_pneumoniae", "penicillin", "scm"),
    ("same-study", "streptococcus_pneumoniae", "penicillin", "lasso"),
    ("same-study", "enterococcus_faecium", "gentamicin", "scm"),
    ("same-study", "enterobacter_cloacae", "sulfamethoxazole_trimethoprim", "scm"),
    ("same-study", "salmonella_enterica", "amoxicillin_clavulanate", "scm"),
    ("same-study", "salmonella_enterica", "amoxicillin_clavulanate", "known_gene"),
    ("same-study", "salmonella_enterica", "trimethoprim", "scm"),
    ("same-study", "salmonella_enterica", "trimethoprim", "cart"),
]

# Cited rows that the text also reads study by study, from `groups.tsv`.
BY_STUDY = [
    ("primary", "streptococcus_pneumoniae", "chloramphenicol", "scm"),
    ("primary", "streptococcus_pneumoniae", "penicillin", "scm"),
]


def rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def title(slug):
    genus, _, rest = slug.partition("_")
    return f"*{genus[0].upper()}. {rest.replace('_', ' ')}*"


def drug(dataset):
    slug = dataset.split("__", 1)[1]
    return slug.replace("_", "/" if slug in COMBINATIONS else " ")


def num(n):
    return f"{int(n):,}"


def pct(x):
    return f"{100 * float(x):.1f}%"


def frozen_and_scored(reading):
    """The counts of the reading's line 'Frozen genomes: N. Scored: M (K of them borderline ...'."""
    m = re.search(r"^Frozen genomes: (\d+)\. Scored: (\d+) \((\d+) of them borderline", reading, re.M)
    if not m:
        raise SystemExit("no 'Frozen genomes' line in a reading.txt")
    return tuple(int(x) for x in m.groups())


def excluded(reading):
    return sum(int(x) for x in re.findall(r"^Genomes excluded, [^:]+: (\d+)\.", reading, re.M))


def same_kmer_set(reading):
    """The two counts of the reading's line on scored genomes with another genome's k-mer set."""
    m = re.search(r"^Scored genomes with the same k-mer set as a cohort genome: (\d+); as another scored genome: (\d+) ", reading, re.M)
    if not m:
        raise SystemExit("no 'same k-mer set' line in a reading.txt")
    return int(m.group(1)), int(m.group(2))


def attempts(reading):
    """'1' when every stage ran once, else the reading's own list of attempts per stage."""
    m = re.search(r"^Stage attempts: (.+)\.$", reading, re.M)
    if not m:
        raise SystemExit("no 'Stage attempts' line in a reading.txt")
    counts = [part.rsplit(" ", 1) for part in m.group(1).split(", ")]
    return "1" if all(n == "1" for _, n in counts) else m.group(1)


def tally(metrics, learner, column):
    got = [r[column] for r in metrics if r["learner"] == learner and r["reportable"] == "yes"]
    return [sum(1 for g in got if g == k) for k in READINGS] + [sum(1 for g in got if g not in READINGS)]


def cell(t):
    if sum(t) == 0:
        return "none scored"
    return " / ".join(str(x) for x in t[:3]) + (f" / {t[3]}" if t[3] else "")


def one_sided(metrics):
    """Known-gene rows with estimates where the rule calls no scored genome resistant, or every one."""
    rule = [r for r in metrics if r["learner"] == "known_gene" and r["reportable"] == "yes"]
    none = sum(1 for r in rule if int(r["tp"]) + int(r["fp"]) == 0)
    every = sum(1 for r in rule if int(r["tn"]) + int(r["fn"]) == 0)
    return len(rule), none, every


def holds(r):
    """True or False when both of a row's readings can be read: neither is "below". None otherwise."""
    if r["reportable"] != "yes" or r["reading_sensitivity"] not in READINGS or r["reading_specificity"] not in READINGS:
        return None
    return r["reading_sensitivity"] != "below" and r["reading_specificity"] != "below"


def load(root):
    """{population: {species: (reading, datasets, metrics, groups)}}, in a fixed order."""
    out = {}
    for pop, _ in POPULATIONS:
        base = root / pop
        if not base.is_dir():
            continue
        runs = {}
        for d in sorted(p for p in base.iterdir() if (p / "metrics.tsv").is_file()):
            reading = (d / "reading.txt").read_text(encoding="utf-8")
            runs[d.name] = (reading, rows(d / "datasets.tsv"), rows(d / "metrics.tsv"), rows(d / "groups.tsv"))
        if runs:
            out[pop] = runs
    return out


def population_tables(heading, runs):
    out = [f"**{heading}**", ""]
    out.append("| Species | Frozen genomes | Scored | Borderline in the species check | Excluded | Datasets | With estimates | Studies on a dataset |")
    out.append("|---|---:|---:|---:|---:|---:|---:|---|")
    totals = [0, 0, 0, 0, 0, 0]
    for s, (reading, datasets, _, _) in runs.items():
        frozen, scored, borderline = frozen_and_scored(reading)
        ex = excluded(reading)
        est = sum(1 for r in datasets if r["reportable"] == "yes")
        studies = sorted(int(r["studies"]) for r in datasets)
        span = str(studies[0]) if studies[0] == studies[-1] else f"{studies[0]} to {studies[-1]}"
        out.append(f"| {title(s)} | {num(frozen)} | {num(scored)} | {num(borderline)} | {num(ex)} | {len(datasets)} | {est} | {span} |")
        for i, v in enumerate((frozen, scored, borderline, ex, len(datasets), est)):
            totals[i] += v
    out.append(f"| All {len(runs)} | {num(totals[0])} | {num(totals[1])} | {num(totals[2])} | {num(totals[3])} | {totals[4]} | {totals[5]} | |")
    out.append("")

    out.append("Scored genomes whose k-mer set equals another genome's (listed, never excluded), and the attempts each stage of the run took.")
    out.append("")
    out.append("| Species | Same k-mer set as a cohort genome | Same k-mer set as another scored genome | Attempts per stage |")
    out.append("|---|---:|---:|---|")
    for s, (reading, _, _, _) in runs.items():
        cohort, scored = same_kmer_set(reading)
        out.append(f"| {title(s)} | {cohort} | {scored} | {attempts(reading)} |")
    out.append("")

    out.append("Where the scored genomes of a dataset with estimates come from. A study is a BioProject.")
    out.append("")
    out.append("| Species | Datasets with estimates | One study holds at least half of the resistant genomes | One study holds at least half of the susceptible genomes | Most resistant and most susceptible genomes come from different studies | No study can be left out |")
    out.append("|---|---:|---:|---:|---:|---:|")
    acc = [0, 0, 0, 0, 0]
    for s, (_, datasets, _, _) in runs.items():
        est = [r for r in datasets if r["reportable"] == "yes"]
        got = [
            len(est),
            sum(1 for r in est if float(r["top_study_resistant_share"]) >= 0.5),
            sum(1 for r in est if float(r["top_study_susceptible_share"]) >= 0.5),
            sum(1 for r in est if r["top_study_resistant"] != r["top_study_susceptible"]),
            sum(1 for r in est if int(r["loso_folds_reportable"]) == 0),
        ]
        acc = [a + b for a, b in zip(acc, got)]
        out.append(f"| {title(s)} | " + " | ".join(str(x) for x in got) + " |")
    out.append(f"| All {len(runs)} | " + " | ".join(str(x) for x in acc) + " |")
    out.append("")

    out.append("The study that most often holds the most resistant genomes of a dataset with estimates, and its share of them there.")
    out.append("")
    out.append("| Species | Study | Datasets where it holds the most resistant genomes | Its share of them, lowest to highest |")
    out.append("|---|---|---:|---|")
    for s, (_, datasets, _, _) in runs.items():
        est = [r for r in datasets if r["reportable"] == "yes"]
        count = {}
        for r in est:
            count[r["top_study_resistant"]] = count.get(r["top_study_resistant"], 0) + 1
        study = sorted(count, key=lambda k: (-count[k], k))[0]
        shares = sorted(float(r["top_study_resistant_share"]) for r in est if r["top_study_resistant"] == study)
        span = pct(shares[0]) if shares[0] == shares[-1] else f"{pct(shares[0])} to {pct(shares[-1])}"
        out.append(f"| {title(s)} | {study} | {count[study]} of {len(est)} | {span} |")
    out.append("")

    for column, word in (("reading_sensitivity", "Sensitivity"), ("reading_specificity", "Specificity")):
        out.append(f"{word} against the report's lineage-split figure, in datasets with estimates. Each cell counts datasets as below / consistent / above; a fourth number counts rows the plan leaves unread.")
        out.append("")
        out.append("| Species | " + " | ".join(NAMES[l] for l in LEARNERS) + " |")
        out.append("|---|" + "---|" * len(LEARNERS))
        sums = {l: [0, 0, 0, 0] for l in LEARNERS}
        for s, (_, _, metrics, _) in runs.items():
            cells = []
            for l in LEARNERS:
                t = tally(metrics, l, column)
                sums[l] = [a + b for a, b in zip(sums[l], t)]
                cells.append(cell(t))
            out.append(f"| {title(s)} | " + " | ".join(cells) + " |")
        out.append(f"| All {len(runs)} | " + " | ".join(cell(sums[l]) for l in LEARNERS) + " |")
        out.append("")

    out.append("The two readings of a row together, over all the species above. A row is one dataset with estimates and one learner.")
    out.append("")
    out.append("| Learner | Rows read | Neither reading below | Sensitivity below only | Specificity below only | Both below |")
    out.append("|---|---:|---:|---:|---:|---:|")
    for l in LEARNERS:
        both = [
            (r["reading_sensitivity"] == "below", r["reading_specificity"] == "below")
            for _, _, metrics, _ in runs.values()
            for r in metrics
            if r["learner"] == l and r["reportable"] == "yes" and r["reading_sensitivity"] in READINGS and r["reading_specificity"] in READINGS
        ]
        got = [len(both)] + [sum(1 for b in both if b == want) for want in ((False, False), (True, False), (False, True), (True, True))]
        out.append(f"| {NAMES[l]} | " + " | ".join(str(x) for x in got) + " |")
    out.append("")

    out.append("Datasets where a learner holds, of the datasets read for it. A learner holds when neither its sensitivity nor its specificity reads below.")
    out.append("")
    out.append("| Species | " + " | ".join(NAMES[l] for l in LEARNERS) + " |")
    out.append("|---|" + "---|" * len(LEARNERS))
    sums = {l: [0, 0] for l in LEARNERS}
    for s, (_, _, metrics, _) in runs.items():
        cells = []
        for l in LEARNERS:
            got = [holds(r) for r in metrics if r["learner"] == l and holds(r) is not None]
            sums[l] = [sums[l][0] + sum(got), sums[l][1] + len(got)]
            cells.append(f"{sum(got)} of {len(got)}" if got else "none read")
        out.append(f"| {title(s)} | " + " | ".join(cells) + " |")
    out.append(f"| All {len(runs)} | " + " | ".join(f"{sums[l][0]} of {sums[l][1]}" for l in LEARNERS) + " |")
    out.append("")

    out.append("Datasets with estimates where the known-gene rule calls every scored genome the same way.")
    out.append("")
    out.append("| Species | Rule rows with estimates | Calls no genome resistant | Calls every genome resistant |")
    out.append("|---|---:|---:|---:|")
    acc = [0, 0, 0]
    for s, (_, _, metrics, _) in runs.items():
        got = one_sided(metrics)
        acc = [a + b for a, b in zip(acc, got)]
        out.append(f"| {title(s)} | {got[0]} | {got[1]} | {got[2]} |")
    out.append(f"| All {len(runs)} | {acc[0]} | {acc[1]} | {acc[2]} |")
    out.append("")
    return out


def both_populations_table(data):
    """A dataset in both populations: one final model and one report figure, read on two sets of new genomes."""
    first, second = (data[pop] for pop, _ in POPULATIONS)
    out = ["**Datasets scored in both populations**", ""]
    out.append(
        "A dataset in both populations has one final model and one report figure, read on two sets of new genomes. "
        "A learner holds when neither its sensitivity nor its specificity reads below. No genome counts are added. "
        "The plan does not make this comparison; it was added after both results were read."
    )
    out.append("")
    out.append("| Learner | Rows read in both | Holds in both | Holds in the same-study population only | Holds in the primary population only | Holds in neither |")
    out.append("|---|---:|---:|---:|---:|---:|")
    for l in LEARNERS:
        pairs = []
        for s in first:
            if s not in second:
                continue
            other = {r["dataset"]: holds(r) for r in second[s][2] if r["learner"] == l}
            for r in first[s][2]:
                if r["learner"] == l and holds(r) is not None and other.get(r["dataset"]) is not None:
                    pairs.append((holds(r), other[r["dataset"]]))
        got = [len(pairs)] + [sum(1 for pair in pairs if pair == want) for want in ((True, True), (False, True), (True, False), (False, False))]
        out.append(f"| {NAMES[l]} | " + " | ".join(str(x) for x in got) + " |")
    out.append("")

    out.append(
        "The same rows by the estimate alone, with no interval. A small batch has a wide interval, so it reads below less often; "
        "this count does not depend on the batch's size. Each cell counts rows as higher in the same-study population / equal / higher in the primary population."
    )
    out.append("")
    out.append("| Learner | Rows read in both | Sensitivity | Specificity |")
    out.append("|---|---:|---|---|")
    for l in LEARNERS:
        cells, n = [], 0
        for k in ("sensitivity", "specificity"):
            got = [0, 0, 0]
            for s in first:
                if s not in second:
                    continue
                other = {r["dataset"]: r for r in second[s][2] if r["learner"] == l and holds(r) is not None}
                for r in first[s][2]:
                    if r["learner"] == l and holds(r) is not None and r["dataset"] in other:
                        a, b = float(other[r["dataset"]][k]), float(r[k])
                        got[0 if a > b else 1 if a == b else 2] += 1
            n = sum(got)
            cells.append(" / ".join(str(x) for x in got))
        out.append(f"| {NAMES[l]} | {n} | " + " | ".join(cells) + " |")
    out.append("")
    return out


def estimate(r, k):
    """'74.2% (55.4% to 88.1%), 23 of 31': the new estimate, its interval and its counts."""
    hit, miss = (int(r["tp"]), int(r["fn"])) if k == "sensitivity" else (int(r["tn"]), int(r["fp"]))
    return f"{pct(r[k])} ({pct(r[k + '_lo'])} to {pct(r[k + '_hi'])}), {num(hit)} of {num(hit + miss)}"


def reference(r, k):
    v = r["lineage_" + k]
    return "none" if v == "" else pct(v)


def cited_row(data, pop, species, d, learner):
    hit = [r for r in data[pop][species][2] if r["dataset"] == f"{species}__{d}" and r["learner"] == learner]
    if len(hit) != 1 or hit[0]["reportable"] != "yes":
        raise SystemExit(f"cited row not found with estimates: {pop} {species} {d} {learner}")
    return hit[0]


def by_study_table(data):
    out = ["**Cited rows, study by study**", ""]
    out.append("The same counts as the row above them in `metrics.tsv`, split by the genome's BioProject.")
    out.append("")
    out.append("| Population | Species | Drug | Learner | Study | Resistant genomes | Of them called resistant | Susceptible genomes | Of them called resistant |")
    out.append("|---|---|---|---|---|---:|---:|---:|---:|")
    for pop, species, d, learner in BY_STUDY:
        whole = cited_row(data, pop, species, d, learner)
        part = [r for r in data[pop][species][3] if r["dataset"] == whole["dataset"] and r["learner"] == learner and r["grouping"] == "bioproject"]
        for k in ("tp", "fn", "fp", "tn"):
            if sum(int(r[k]) for r in part) != int(whole[k]):
                raise SystemExit(f"the studies of {pop} {species} {d} {learner} do not add up to the row's {k}")
        for r in part:
            out.append(
                f"| {pop} | {title(species)} | {drug(whole['dataset'])} | {NAMES[learner]} | {r['group']} "
                f"| {num(r['scored_resistant'])} | {num(r['tp'])} | {num(r['scored_susceptible'])} | {num(r['fp'])} |"
            )
    out.append("")
    return out


def cited_table(data):
    out = ["**Rows cited in the text**", ""]
    out.append("Each estimate is followed by its 95% interval and its counts. The report's figure is the lineage-split figure for a learner, and the whole-cohort figure for the known-gene rule.")
    out.append("")
    out.append("| Population | Species | Drug | Learner | Sensitivity | Report | Reading | Specificity | Report | Reading |")
    out.append("|---|---|---|---|---|---|---|---|---|---|")
    for pop, species, d, learner in CITED:
        r = cited_row(data, pop, species, d, learner)
        out.append(
            f"| {pop} | {title(species)} | {drug(r['dataset'])} | {NAMES[learner]} | {estimate(r, 'sensitivity')} | {reference(r, 'sensitivity')} | {r['reading_sensitivity']} "
            f"| {estimate(r, 'specificity')} | {reference(r, 'specificity')} | {r['reading_specificity']} |"
        )
    out.append("")
    return out


def build(root):
    data = load(pathlib.Path(root))
    out = []
    for pop, heading in POPULATIONS:
        if pop in data:
            out += population_tables(heading, data[pop])
    if all(pop in data for pop, _ in POPULATIONS):
        out += both_populations_table(data)
    if any(key[0] in data and key[1] in data[key[0]] for key in CITED):
        out += cited_table(data)
        out += by_study_table(data)
    return "\n".join(out).rstrip("\n") + "\n"


def check(root, document):
    """True when the document's table block is exactly what build() prints."""
    doc = pathlib.Path(document).read_text(encoding="utf-8")
    m = re.search(re.escape(BEGIN) + r"\n(.*?)" + re.escape(END), doc, re.S)
    return bool(m) and m.group(1).strip("\n") == build(root).strip("\n")


def main(argv):
    if len(argv) not in (2, 4) or (len(argv) == 4 and argv[2] != "--check"):
        raise SystemExit(__doc__)
    if len(argv) == 4:
        if check(argv[1], argv[3]):
            return 0
        print("the tables in the document differ from the output files", file=sys.stderr)
        return 1
    sys.stdout.write(build(argv[1]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
