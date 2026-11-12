"""The committed outside-study results (reference/outside-study) and the two pages that read them.

The suite scores nothing. It checks that each page's tables are the ones the
output files give, that the sentences which carry a count agree with the files,
that the files are the listed ones, and that each run's record names the frozen
kit, its own population and the frozen genomes.
"""

import csv
import gzip
import hashlib
import importlib.util
import json
import pathlib
import re
import shutil
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
RESULTS = ROOT / "reference" / "outside-study"
DOCUMENT = ROOT / "docs" / "outside-study-results.md"
SPECIES_PAGE = ROOT / "docs" / "outside-study-results-by-species.md"
FROZEN = ROOT / "scripts" / "outside-study" / "frozen"
KIT = "e3e1c345a55412ac8f7d07b691ebb56573d2d3d8"
NAME = "a retrospective test on newly deposited studies from the same database"
FILES = sorted(
    [
        "datasets.tsv",
        "duplicates.tsv",
        "features.tsv",
        "genomes.tsv",
        "groups.tsv",
        "leave_one_study_out.tsv",
        "metrics.tsv",
        "paired.tsv",
        "reading_counts.tsv",
        "reading.txt",
        "results.json.gz",
    ]
)
# docs/outside-study-test.md, item 1: the value of `population` in the frozen
# candidates and the flag in the frozen dataset table, per population.
POPULATIONS = {"primary": ("primary (other-study)", "eligible_primary"), "same-study": ("same-study", "eligible_same_study")}
STAGES = ["count", "fetch", "known-gene", "predict", "setup", "species-check"]


def load_tables():
    spec = importlib.util.spec_from_file_location("outside_study_results_tables", ROOT / "scripts" / "outside-study-results" / "tables.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T = load_tables()


def tsv(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def record(run):
    with gzip.open(run / "results.json.gz", "rt", encoding="utf-8") as f:
        return json.load(f)


def runs():
    return sorted(p for pop in POPULATIONS for p in (RESULTS / pop).iterdir() if p.is_dir())


def blocks(reading):
    """{dataset: the indented lines under its header line} of a reading.txt.

    A dataset's block is its header line, `<dataset>: ...`, and the indented
    lines under it. No blank line separates two blocks."""
    out, name = {}, None
    for line in reading.split("\n"):
        if line.startswith(" "):
            if name is not None:
                out[name].append(line)
            continue
        head = line.split(": ", 1)[0]
        name = head if "__" in head and " " not in head and ": " in line else None
        if name is not None:
            if name in out:
                raise AssertionError(f"two blocks for {name}")
            out[name] = []
    return out


def prose(page=None):
    """A page outside its generated table block."""
    text = (page or DOCUMENT).read_text(encoding="utf-8")
    return text[: text.index(T.BEGIN)] + text[text.index(T.END) + len(T.END) :]


class Tables(unittest.TestCase):
    def test_the_document_holds_the_tables_the_files_give(self):
        self.assertTrue(T.check(RESULTS, DOCUMENT))
        self.assertTrue(T.check(RESULTS, SPECIES_PAGE, T.build_species))
        # Each page holds its own tables and not the other's.
        self.assertFalse(T.check(RESULTS, DOCUMENT, T.build_species))
        self.assertFalse(T.check(RESULTS, SPECIES_PAGE))

    def test_a_changed_figure_in_a_cited_row_fails_the_check(self):
        text = SPECIES_PAGE.read_text(encoding="utf-8")
        row = T.cited_row(T.load(RESULTS), *T.CITED[0])
        figure = T.estimate(row, "sensitivity")
        self.assertEqual(text.count(figure), 1)
        with tempfile.TemporaryDirectory() as tmp:
            changed = pathlib.Path(tmp) / "changed.md"
            changed.write_text(text.replace(figure, "0.0% (0.0% to 0.0%), 0 of 1"), encoding="utf-8")
            self.assertFalse(T.check(RESULTS, changed, T.build_species))

    def test_a_changed_cell_in_the_document_fails_the_check(self):
        text = DOCUMENT.read_text(encoding="utf-8")
        begin, end = text.index(T.BEGIN), text.index(T.END)
        block = text[begin:end]
        self.assertIn(" below ", block)
        with tempfile.TemporaryDirectory() as tmp:
            changed = pathlib.Path(tmp) / "changed.md"
            changed.write_text(text[:begin] + block.replace(" below ", " above ", 1) + text[end:], encoding="utf-8")
            self.assertFalse(T.check(RESULTS, changed))
            unmarked = pathlib.Path(tmp) / "unmarked.md"
            unmarked.write_text(text.replace(T.BEGIN, ""), encoding="utf-8")
            self.assertFalse(T.check(RESULTS, unmarked))

    def test_a_changed_reading_in_an_output_file_fails_the_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            copy = pathlib.Path(tmp) / "outside-study"
            shutil.copytree(RESULTS, copy)
            metrics = sorted(copy.glob("primary/*/metrics.tsv"))[0]
            text = metrics.read_text(encoding="utf-8")
            self.assertIn("\tbelow\t", text)
            metrics.write_text(text.replace("\tbelow\t", "\tabove\t", 1), encoding="utf-8")
            self.assertFalse(T.check(copy, DOCUMENT))

    def test_the_two_populations_are_printed_apart(self):
        text = T.build(RESULTS)
        self.assertEqual(text.count("**Primary population**"), 1)
        self.assertEqual(text.count("**Same-study population**"), 1)
        self.assertLess(text.index("**Primary population**"), text.index("**Same-study population**"))
        primary = text[: text.index("**Same-study population**")]
        frozen = sum(T.frozen_and_scored((run / "reading.txt").read_text(encoding="utf-8"))[0] for run in runs() if run.parent.name == "primary")
        self.assertIn(f"| All {len(list((RESULTS / 'primary').iterdir()))} | {frozen:,} |", primary)

    def test_a_dataset_in_both_populations_is_counted_once(self):
        text = T.build(RESULTS)
        after = text.split("**Datasets scored in both populations**", 1)[1]
        block = after.split("The same rows by the estimate alone", 1)
        self.assertEqual(len(block), 2)
        block = [block[0]] + block[1].split("The same-study genomes of the datasets with estimates in both populations", 1)
        self.assertEqual(len(block), 3)
        tables = [[[c.strip() for c in l.strip("|").split("|")] for l in part.split("\n") if l.startswith("| ")][1:] for part in block]
        for got in tables[:2]:
            self.assertEqual([g[0] for g in got], [T.NAMES[l] for l in T.LEARNERS])
        # The four ways a row can hold add up to the rows read in both.
        for g in tables[0]:
            self.assertEqual(int(g[1]), sum(int(x) for x in g[2:]), g)
        # Each measure's three counts add up to the same rows.
        for g, first in zip(tables[1], tables[0]):
            self.assertEqual(g[1], first[1])
            for cell in g[2:]:
                self.assertEqual(int(g[1]), sum(int(x) for x in cell.split(" / ")), g)
        # Each dataset in both populations is under one study, and the studies add up to the datasets.
        studies, total = tables[2][:-1], tables[2][-1]
        self.assertEqual(total[0], f"All {len(studies)}")
        self.assertEqual(int(total[1]), sum(int(g[1]) for g in studies))
        data = T.load(RESULTS)
        shared = 0
        for species, run in data["same-study"].items():
            if species in data["primary"]:
                here = {d["dataset"] for d in data["primary"][species][1] if d["reportable"] == "yes"}
                shared += sum(1 for d in run[1] if d["reportable"] == "yes" and d["dataset"] in here)
        self.assertEqual(int(total[1]), shared)
        # A learner is read on no more rows than there are datasets in both.
        for g in tables[0]:
            self.assertLessEqual(int(g[1]), shared, g)

    def test_the_different_studies_column_counts_the_datasets_the_kit_names(self):
        # docs/outside-study-test.md, Outputs, item 8: reading.txt says that the two
        # measures are measured mostly on different studies where more than half of
        # the resistant genomes are from one study and more than half of the
        # susceptible genomes from another. The column counts those datasets.
        data = T.load(RESULTS)
        named = 0
        for pop, runs_ in data.items():
            for species, (reading, datasets, _, _) in runs_.items():
                want = {d["dataset"] for d in datasets if d["reportable"] == "yes" and T.different_majority_studies(d)}
                got = {name for name, block in blocks(reading).items() if any("measured mostly on different studies" in l for l in block)}
                self.assertEqual(got, want, (pop, species))
                named += len(got)
        self.assertGreater(named, 0)

    def test_a_tied_or_shared_top_study_is_not_counted_as_different_studies(self):
        row = {"top_study_resistant": "A", "top_study_resistant_share": "0.600000", "top_study_susceptible": "B", "top_study_susceptible_share": "0.700000"}
        self.assertTrue(T.different_majority_studies(row))
        # A plurality is not a majority: the top studies differ and neither holds more than half.
        self.assertFalse(T.different_majority_studies(dict(row, top_study_resistant_share="0.403000", top_study_susceptible_share="0.403000")))
        # Exactly half can be a tie between two studies, so it does not count.
        self.assertFalse(T.different_majority_studies(dict(row, top_study_susceptible_share="0.500000")))
        self.assertFalse(T.different_majority_studies(dict(row, top_study_resistant_share="0.500000")))
        self.assertFalse(T.different_majority_studies(dict(row, top_study_susceptible="A")))

    def test_a_row_at_an_edge_is_counted_as_one(self):
        row = {"lineage_sensitivity": "0.000000", "lineage_specificity": "0.930000", "tp": "0", "fn": "30", "fp": "0", "tn": "40"}
        self.assertTrue(T.at_the_edge(row))
        self.assertTrue(T.at_the_edge(dict(row, lineage_sensitivity="0.800000", lineage_specificity="1.000000")))
        self.assertFalse(T.at_the_edge(dict(row, lineage_sensitivity="0.800000")))
        self.assertTrue(T.calls_one_way(row))
        self.assertTrue(T.calls_one_way(dict(row, tp="30", fn="0", fp="40", tn="0")))
        self.assertFalse(T.calls_one_way(dict(row, tp="1", fn="29")))

    def test_a_cited_row_that_is_missing_is_refused(self):
        data = T.load(RESULTS)
        with self.assertRaises(SystemExit):
            T.cited_row(data, "primary", "staphylococcus_aureus", "no_such_drug", "scm")

    def test_every_cited_figure_is_the_kit_s_own_text(self):
        # reading.txt gives each estimate in plain words. The percentages this
        # script formats from metrics.tsv must be the ones written there.
        data = T.load(RESULTS)
        for pop, species, drug, learner in T.CITED:
            row = T.cited_row(data, pop, species, drug, learner)
            block = blocks(data[pop][species][0])[row["dataset"]]
            line = [l for l in block if l.startswith(f"  {learner}: ")]
            self.assertEqual(len(line), 1, (pop, species, drug, learner))
            for k in ("sensitivity", "specificity"):
                want = f"{k.capitalize()} {T.pct(row[k])} [{T.pct(row[k + '_lo'])}, {T.pct(row[k + '_hi'])}]"
                self.assertIn(want, line[0], (pop, species, drug, learner))
            self.assertIn(f"TP {row['tp']}, FN {row['fn']}, FP {row['fp']}, TN {row['tn']}.", line[0])


class Prose(unittest.TestCase):
    """The sentences outside the table block that carry a count.

    Each is rebuilt here from the output files and must be in the text as
    written, so a count in a sentence cannot drift from the files. A sentence
    that says "each", "no" or "more often" is checked for every learner it
    covers."""

    @classmethod
    def setUpClass(cls):
        cls.data = T.load(RESULTS)
        cls.text = prose()

    def has(self, sentence):
        self.assertIn(sentence, self.text)

    def read(self, pop, learner):
        return [r for _, _, metrics, _ in self.data[pop].values() for r in metrics if r["learner"] == learner and T.holds(r) is not None]

    def held(self, pop, learner):
        got = [T.holds(r) for r in self.read(pop, learner)]
        return sum(got), len(got)

    def tally(self, pop, learner, column):
        return [sum(T.tally(metrics, learner, column)[i] for _, _, metrics, _ in self.data[pop].values()) for i in range(3)]

    def estimates(self, pop):
        return [d for _, datasets, _, _ in self.data[pop].values() for d in datasets if d["reportable"] == "yes"]

    def in_both(self, learner):
        """(primary row, same-study row) for each dataset the learner is read on in both populations."""
        other = {r["dataset"]: r for r in self.read("same-study", learner)}
        return [(r, other[r["dataset"]]) for r in self.read("primary", learner) if r["dataset"] in other]

    def test_the_headline_is_true_of_each_learner_it_names(self):
        self.has(
            "On genomes from studies the cohort does not hold, each learner's figures hold on fewer than half of the datasets read for it. "
            "On new genomes from studies the cohort already holds, each learner's figures hold on more than half."
        )
        self.has("Those two sentences are about the three learners: the set covering machine, the decision tree and the lasso.")
        for learner in ("scm", "cart", "lasso"):
            h, n = self.held("primary", learner)
            self.assertLess(2 * h, n, learner)
            h, n = self.held("same-study", learner)
            self.assertGreater(2 * h, n, learner)

    def test_the_counts_of_datasets_where_a_learner_holds(self):
        h, n = self.held("primary", "scm")
        self.has(f"- The set covering machine holds on {h} of the {n} datasets and reads below on the other {n - h}.")
        ch, cn = self.held("primary", "cart")
        self.has(f"- The decision tree holds on {ch} of the {cn} datasets.")
        lh, ln = self.held("primary", "lasso")
        self.has(f"- The lasso can be scored and read on {ln} of the {n} datasets, and holds on {lh} of them.")
        kh, kn = self.held("primary", "known_gene")
        self.has(f"It can be read on {kn} datasets. It holds on {kh} of them and reads below on at least one measure on the other {kn - kh}.")
        edge = [r for r in self.read("primary", "known_gene") if T.at_the_edge(r)]
        self.assertTrue(all(T.holds(r) for r in edge))
        self.has(f"On {len(edge)} of the {kh} its report figure is exactly 0% or 100% on a measure")
        flat = [r for r in self.read("primary", "known_gene") if T.calls_one_way(r)]
        self.has(
            f"In the primary population {len(edge)} of the {kn} known-gene rows read have such a figure, and the rule holds on all {len(edge)}. "
            f"On {len(flat)} of the {kn} it gives every scored genome the same call, and it holds on {sum(T.holds(r) for r in flat)} of those."
        )
        with_estimates, none, every = (sum(x) for x in zip(*(T.one_sided(metrics) for _, _, metrics, _ in self.data["primary"].values())))
        self.has(f"In the primary population it reads below on at least one measure on {kn - kh} of its {kn} datasets read.")
        self.has(f"The known-gene rule gives every scored genome the same call on {none + every} of its {with_estimates} primary datasets with estimates.")

        h, n = self.held("same-study", "scm")
        self.has(f"- The set covering machine holds on {h} of the {n} datasets.")
        ch, cn = self.held("same-study", "cart")
        self.has(f"- The decision tree holds on {ch} of the {cn} datasets.")
        lh, ln = self.held("same-study", "lasso")
        self.has(f"- The lasso can be scored and read on {ln} of the {n} datasets, and holds on {lh} of them.")
        kh, kn = self.held("same-study", "known_gene")
        self.has(f"- The known-gene rule holds on {kh} of the {kn} datasets.")

    def test_the_directions_the_figures_move(self):
        sens, spec = self.tally("primary", "scm", "reading_sensitivity"), self.tally("primary", "scm", "reading_specificity")
        read = self.read("primary", "scm")
        n = len(read)
        self.has(
            f"sensitivity reads above the report's figure on {sens[2]} of the {n} datasets and below it on {sens[0]}. "
            f"Its specificity reads below on {spec[0]} datasets and above on {spec[2]}."
        )
        below = [(r["reading_sensitivity"] == "below", r["reading_specificity"] == "below") for r in read]
        self.has(
            f"Only its specificity reads below on {below.count((False, True))} datasets, "
            f"only its sensitivity on {below.count((True, False))}, and both on {below.count((True, True))}."
        )
        flat = sum(1 for r in read if r["reading_sensitivity"] == r["reading_specificity"] == "consistent")
        self.has(f"both readings are consistent on {flat} of the {n} datasets. On the other {n - flat} at least one figure moved, up or down.")
        self.has(
            "In the primary population, for each learner and for the known-gene rule, sensitivity read above the report's figure more often than below it, "
            "and specificity read below more often than above."
        )
        self.has("In the same-study population specificity read above more often than below.")
        for learner in T.LEARNERS:
            sens, spec = self.tally("primary", learner, "reading_sensitivity"), self.tally("primary", learner, "reading_specificity")
            self.assertGreater(sens[2], sens[0], learner)
            self.assertGreater(spec[0], spec[2], learner)
            spec = self.tally("same-study", learner, "reading_specificity")
            self.assertGreater(spec[2], spec[0], learner)

    def test_the_comparison_of_the_two_populations(self):
        here = {d["dataset"] for d in self.estimates("primary")}
        shared = [d for d in self.estimates("same-study") if d["dataset"] in here]
        n = len(shared)
        self.has(f"{n} datasets are in both populations.")
        pairs = {l: [(T.holds(a), T.holds(b)) for a, b in self.in_both(l)] for l in T.LEARNERS}
        self.assertEqual(len(pairs["scm"]), n)
        first, second = sum(a for a, _ in pairs["scm"]), sum(b for _, b in pairs["scm"])
        self.has(f"The set covering machine holds on {second} of the {n} in the same-study population and on {first} of the {n} in the primary population.")
        self.has(f"the set covering machine holds on {second} with same-study genomes and on {first} with genomes from other studies.")
        alone = {l: pairs[l].count((True, False)) for l in T.LEARNERS}
        self.assertEqual((alone["scm"], alone["lasso"]), (0, 0))
        self.has("For the set covering machine, no dataset holds in the primary population alone.")
        self.has(f"The decision tree has {alone['cart']} such dataset and the known-gene rule has {alone['known_gene']}.")
        self.has("On these 34 datasets a learner's lineage-split figure held more often on new genomes from a study the cohort already holds than on genomes from other studies.".replace("34", str(n)))
        for learner in ("scm", "cart", "lasso"):
            self.assertGreater(sum(b for _, b in pairs[learner]), sum(a for a, _ in pairs[learner]), learner)

        def higher(learner, k):
            got = [0, 0, 0]
            for a, b in self.in_both(learner):
                got[0 if float(b[k]) > float(a[k]) else 1 if float(b[k]) == float(a[k]) else 2] += 1
            return got

        self.has("A count with no interval gives the same direction for the set covering machine and the decision tree.")
        for learner in ("scm", "cart"):
            for k in ("sensitivity", "specificity"):
                got = higher(learner, k)
                self.assertGreater(got[0], got[2], (learner, k))
        self.has(
            f"The set covering machine's sensitivity is higher in the same-study population on {higher('scm', 'sensitivity')[0]} of the {n} datasets, "
            f"and its specificity on {higher('scm', 'specificity')[0]}."
        )
        lasso = higher("lasso", "specificity")
        self.assertEqual((lasso[0], lasso[1]), (lasso[2], 0))
        self.has(f"The lasso is read on {len(pairs['lasso'])} of the {n}, and its specificity is higher in each population on {lasso[0]} of them.")
        rule = higher("known_gene", "sensitivity")
        self.assertEqual(len(pairs["known_gene"]), n)
        self.has(f"it is higher in the same-study population on {rule[0]} of the {n} and higher in the primary population on {rule[2]}.")

        count = {}
        for d in shared:
            self.assertGreaterEqual(float(d["top_study_resistant_share"]), 0.999, d["dataset"])
            count[d["top_study_resistant"]] = count.get(d["top_study_resistant"], 0) + 1
        order = sorted(count, key=lambda k: (-count[k], k))
        self.assertEqual(len(order), 3)
        self.has(
            f"one study has at least 99.9% of the resistant genomes on each of these {n} datasets, and there are three such studies: "
            f"{order[0]} on {count[order[0]]} datasets, {order[1]} on {count[order[1]]} and {order[2]} on {count[order[2]]}."
        )

    def test_the_limits_and_the_genome_counts(self):
        totals = {}
        for pop in POPULATIONS:
            frozen = [T.frozen_and_scored(reading) for reading, _, _, _ in self.data[pop].values()]
            totals[pop] = [sum(f[i] for f in frozen) for i in range(3)]
            datasets = sum(len(d) for _, d, _, _ in self.data[pop].values())
            name = "primary" if pop == "primary" else "same-study"
            self.has(f"The {name} population is {totals[pop][0]:,} genomes on {datasets} datasets in {len(self.data[pop])} species.")
        self.assertEqual(totals["primary"][0], totals["primary"][1])
        self.has("Every genome was scored and none was excluded.")
        self.has(f"{totals['same-study'][1]:,} of the {totals['same-study'][0]:,} genomes were scored.")

        def borderline(pop, species):
            return T.frozen_and_scored(self.data[pop][species][0])[2]

        self.has(
            f"The species check marks {totals['primary'][2]:,} of the {totals['primary'][0]:,} primary genomes as borderline, "
            f"{borderline('primary', 'escherichia_coli')} of them *E. coli*. "
            f"It marks {totals['same-study'][2]} of the {totals['same-study'][1]:,} scored same-study genomes, "
            f"{borderline('same-study', 'enterobacter_cloacae')} of them *E. cloacae*."
        )
        repeats = {s: T.same_kmer_set(reading)[1] for s, (reading, _, _, _) in self.data["primary"].items()}
        self.assertEqual(sorted(s for s, v in repeats.items() if v), ["salmonella_enterica", "streptococcus_pneumoniae"])
        self.has(
            f"{sum(repeats.values())} scored primary genomes have the k-mer set of another scored genome: "
            f"{repeats['salmonella_enterica']} *S. enterica* and {repeats['streptococcus_pneumoniae']} *S. pneumoniae*."
        )

        primary = self.estimates("primary")
        self.has(f"On {sum(1 for d in primary if float(d['top_study_resistant_share']) >= 0.5)} of the {len(primary)} datasets, a single study holds at least half of the resistant genomes.")
        self.has(
            f"On {sum(1 for d in primary if T.different_majority_studies(d))} of the {len(primary)} primary datasets, "
            "more than half of the resistant genomes come from one study and more than half of the susceptible genomes from another"
        )
        same = self.estimates("same-study")
        for d in same:
            self.assertGreaterEqual(float(d["top_study_resistant_share"]), 0.999, d["dataset"])
        self.has(
            f"On each of its {len(same)} datasets one study holds at least 99.9% of the resistant genomes, "
            f"and {sum(1 for d in same if int(d['studies']) == 1)} of the {len(same)} have one study only."
        )
        studies = {d["top_study_resistant"] for d in same}
        self.assertEqual(studies, {"PRJNA288601", "PRJNA771404", "PRJNA248792", "PRJEB3084"})
        self.has("Four studies fill that place: PRJNA288601, PRJNA771404, PRJNA248792 and PRJEB3084.")

        # The count of older genome records is the plan's own, and the plan is frozen.
        plan = (ROOT / "docs" / "outside-study-test.md").read_text(encoding="utf-8")
        self.assertIn(f"For 260 of the {totals['primary'][0]:,} genomes, the genome record was inserted before July 2025", plan)
        self.has(f"260 of the {totals['primary'][0]:,} primary genomes have a genome record from before July 2025.")


class SpeciesPage(unittest.TestCase):
    """The species page's own sentences: each species' heading and each figure set against a report's."""

    WRITTEN = {"salmonella_enterica": "Salmonella enterica"}

    @classmethod
    def setUpClass(cls):
        cls.data = T.load(RESULTS)
        text = prose(SPECIES_PAGE)
        first, second = text.index("## Primary population"), text.index("## Same-study population")
        cls.sections = {"primary": text[first:second], "same-study": text[second : text.index("## Rows")]}

    def paragraph(self, pop, species):
        name = self.WRITTEN.get(species, T.title(species).strip("*"))
        hit = [p for p in self.sections[pop].split("\n\n") if p.startswith(f"***{name}*, ")]
        self.assertEqual(len(hit), 1, (pop, species))
        return hit[0]

    def test_each_species_heading_gives_the_run_s_genomes_and_datasets(self):
        for pop, runs_ in self.data.items():
            for species, (reading, datasets, _, _) in runs_.items():
                frozen, scored, _ = T.frozen_and_scored(reading)
                genomes = f"{frozen:,} genomes" if scored == frozen else f"{scored:,} genomes scored of {frozen:,}"
                plural = "dataset" if len(datasets) == 1 else "datasets"
                name = self.WRITTEN.get(species, T.title(species).strip("*"))
                self.assertTrue(self.paragraph(pop, species).startswith(f"***{name}*, {genomes}, {len(datasets)} {plural}.**"), (pop, species))

    def test_no_species_of_a_population_is_left_out_or_added(self):
        for pop, runs_ in self.data.items():
            headed = [p for p in self.sections[pop].split("\n\n") if p.startswith("***")]
            self.assertEqual(len(headed), len(runs_), pop)

    def test_each_figure_set_against_a_report_s_is_a_cited_row_s(self):
        # "30.6% sensitivity against 87.0%" and "9.4% against 46.7%": the new
        # estimate and the report's figure of one measure of one cited row of
        # that population. A sentence that looks back at the primary population
        # ("where the primary population gave 22.9%") names a cited estimate there.
        pairs = {pop: set() for pop in self.data}
        estimates = {pop: set() for pop in self.data}
        for pop, species, drug, learner in T.CITED:
            row = T.cited_row(self.data, pop, species, drug, learner)
            for k in ("sensitivity", "specificity"):
                pairs[pop].add((T.pct(row[k]), T.reference(row, k)))
                estimates[pop].add(T.pct(row[k]))
        seen = 0
        for pop, text in self.sections.items():
            for new, old in re.findall(r"(\d+\.\d%)(?: sensitivity| specificity)? against (\d+\.\d%)", text):
                self.assertIn((new, old), pairs[pop], (pop, new, old))
                seen += 1
            for gave in re.findall(r"where the primary population gave (\d+\.\d%)", text):
                self.assertIn(gave, estimates["primary"], gave)
                seen += 1
        self.assertGreater(seen, 20)

    def test_the_datasets_where_the_set_covering_machine_holds(self):
        for species in ("enterococcus_faecium", "escherichia_coli", "klebsiella_pneumoniae", "salmonella_enterica", "streptococcus_pneumoniae"):
            got = [T.holds(r) for r in self.data["same-study"][species][2] if r["learner"] == "scm" and T.holds(r) is not None]
            self.assertIn(f"The set covering machine holds on {sum(got)} of the {len(got)} datasets", self.paragraph("same-study", species), species)


class Files(unittest.TestCase):
    def test_every_run_folder_holds_the_eleven_output_files(self):
        self.assertTrue(runs())
        for run in runs():
            self.assertEqual(sorted(p.name for p in run.iterdir()), FILES, run)

    def test_no_folder_but_the_two_populations(self):
        self.assertEqual(sorted(p.name for p in RESULTS.iterdir()), ["SHA256SUMS", "primary", "same-study"])

    def test_sha256sums_lists_every_file_with_its_hash(self):
        listed = {}
        for line in (RESULTS / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
            digest, name = line.split("  ", 1)
            self.assertNotIn(name, listed)
            listed[name] = digest
        on_disk = sorted(str(p.relative_to(RESULTS)) for p in RESULTS.rglob("*") if p.is_file() and p.name != "SHA256SUMS")
        self.assertEqual(sorted(listed), on_disk)
        for name, digest in listed.items():
            self.assertEqual(hashlib.sha256((RESULTS / name).read_bytes()).hexdigest(), digest, name)


class Records(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.candidates = tsv(FROZEN / "candidates.tsv.gz")
        cls.datasets = tsv(FROZEN / "datasets.tsv")

    def test_each_record_names_the_kit_and_its_own_run(self):
        binaries = set()
        for run in runs():
            r = record(run)
            self.assertEqual(r["kit"]["commit"], KIT, run)
            self.assertEqual(r["kit"]["expect_commit"], KIT, run)
            self.assertEqual(r["name"], NAME, run)
            self.assertEqual((r["population"], r["species"]), (run.parent.name, run.name))
            self.assertIs(r["rehearsal"], False, run)
            self.assertIs(r["table_pins_checked"], True, run)
            self.assertEqual(r["evaluation_revision"], 12, run)
            self.assertEqual(sorted(r["attempts"]), STAGES, run)
            binaries.add(r["grm_sha256"])
        self.assertEqual(len(binaries), 1, "every run used one grm binary")

    def test_each_run_covers_the_frozen_genomes_and_datasets(self):
        for run in runs():
            pop, species = run.parent.name, run.name
            value, flag = POPULATIONS[pop]
            frozen = sorted(c["genome_id"] for c in self.candidates if c["species_slug"] == species and c["population"] == value and c["on_eligible_dataset"] == "yes")
            eligible = sorted(d["dataset_id"] for d in self.datasets if d["species_slug"] == species and d[flag] == "yes")
            self.assertTrue(frozen and eligible, run)
            self.assertEqual(sorted(g["genome_id"] for g in tsv(run / "genomes.tsv")), frozen, run)
            self.assertEqual(sorted(d["dataset"] for d in tsv(run / "datasets.tsv")), eligible, run)
            self.assertEqual(record(run)["genomes_frozen"], len(frozen), run)

    def test_no_species_of_a_frozen_population_is_missing(self):
        for pop, (value, flag) in POPULATIONS.items():
            want = sorted({d["species_slug"] for d in self.datasets if d[flag] == "yes"})
            self.assertEqual(sorted(p.name for p in (RESULTS / pop).iterdir()), want, pop)

    def test_the_counts_in_metrics_are_the_record_s(self):
        for run in runs():
            by_key = {(m["dataset"], m["learner"]): m for m in record(run)["metrics"]}
            rows = tsv(run / "metrics.tsv")
            self.assertEqual(len(rows), len(by_key), run)
            for row in rows:
                m = by_key[(row["dataset"], row["learner"])]
                for k in ("tp", "fn", "fp", "tn"):
                    self.assertEqual(int(row[k]), m[k], (run, row["dataset"], row["learner"], k))
                self.assertEqual(row["reading_sensitivity"], m["reading_sensitivity"])
                self.assertEqual(row["reading_specificity"], m["reading_specificity"])


if __name__ == "__main__":
    unittest.main()
