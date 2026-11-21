"""docs/time-split.md against the committed record: the limit sentence, the section Reading, the exploratory
paragraphs under the study and time table, every table, the sentences that carry the record's counts, every other
paragraph and list item word for word, and every other number of the page."""

import re
import unittest

from tests.time_split import lines
from tests.time_split.common import ARMS, PAGE, RECORD
from tests.time_split.held import HELD
from tests.time_split.lines import n

LIMIT = ("This is a split by first-public day inside one database, BV-BRC. "
         "It is not a prospective test and not external validation.")
SECTIONS = ["Reading", "Plan and order", "Method", "Total contrast", "Study and time contrasts", "Per arm",
            "Test sets", "Review", "Limits", "Files"]
# Where each table the suite builds stands.
TABLE_SECTION = {"order": "Plan and order", "total": "Total contrast", "other": "Study and time contrasts",
                 "mix": "Study and time contrasts", "arms": "Per arm", "species": "Test sets",
                 "review": "Review"}
# The readings the review's fix changed, as the page names them, and what each reads now.
CHANGED = {("cart", "study recall"): "not shown", ("gbt", "study recall"): "not shown",
           ("gbt", "study specificity"): "not shown", ("cart", "time log_loss"): "not shown",
           ("lasso", "time log_loss"): "worse", ("lasso", "time balanced_accuracy"): "worse"}


def outside_code(text):
    """The page's lines outside fenced code blocks."""
    out, fenced = [], False
    for line in text.split("\n"):
        if line.startswith("```"):
            fenced = not fenced
            continue
        if not fenced:
            out.append(line)
    return out


def blocks(text):
    """The page's paragraphs and list items outside code, tables and headers, one string each, in page order."""
    out = []
    for para in "\n".join(outside_code(text)).split("\n\n"):
        para = para.strip("\n")
        if not para or para.startswith(("#", "|")):
            continue
        if para.startswith("- "):
            out += ["- " + x for x in para[2:].split("\n- ")]
        else:
            out.append(para)
    return out


def section(text, title):
    return text.split(f"\n## {title}\n", 1)[1].split("\n## ", 1)[0]


def page_tables(text):
    tables, cur = [], []
    for line in outside_code(text) + [""]:
        if line.startswith("|"):
            cur.append(line)
        elif cur:
            tables.append(cur)
            cur = []
    return tables


def counted(rec):
    """The sentences outside the section Reading that give the record's counts, rebuilt from it."""
    s, t, st = rec.summary, rec.totals, rec.settings
    pct = lambda v: "%g" % (100 * float(v))
    mix = {r["arm"]: r for r in rec.mix}
    species = ", ".join(f"{r['datasets']} of {lines.species_name(r['species'])}" for r in rec.species[:-1])
    species += f" and {rec.species[-1]['datasets']} of {lines.species_name(rec.species[-1]['species'])}"
    scm = rec.arms[("scm", "past")]
    return [
        f"Both are scored on the same later genomes of studies new after that day, in {s['datasets']} datasets of "
        f"{s['species']} species at evaluation revision {st['evaluation_revision']}.",
        f"The earliest day by which at least {pct(st['cut_share'])}% of the cohort was public.",
        f"A dataset enters when its test set holds at least {st['min_test_class']} genomes of each class, each arm "
        f"trains on at least {st['min_train_class']} of each class, and at least {st['min_reserve']} reserve genomes "
        f"and {st['min_later_other_studies']} later genomes of other studies remain in the pools.",
        f"Datasets. {s['datasets']} datasets in {s['species']} species enter: {species}.",
        f"with the default settings of `grm learn` (seed {st['seed']}, {n(st['screened_features'])} screened features, "
        f"{st['inner_folds']} inner folds).",
        f"A score of {st['cutoff']} or more calls resistant, and recall is the share of resistant genomes called "
        f"resistant. Log loss, unlike recall, reads the score itself: it clips each score to [{st['clip']}, 1 - "
        f"{st['clip']}], so one wrong call at the clip adds {float(st['loss_of_a_wrong_call_at_the_clip']):.2f} to a "
        f"sum over the test genomes. scm scores only 0 or 1, so its log loss is about "
        f"{float(st['loss_of_a_wrong_call_at_the_clip']):.2f} times its share of wrong calls,",
        f"Each contrast is read as its median over the {s['datasets']} datasets, with a {pct(st['interval_level'])}% "
        f"interval: {n(st['resamples'])} draws of the {s['species']} species with replacement (seed "
        f"{st['resample_seed']}), each the median over the datasets of the drawn species, from the {st['low_draw']}st "
        f"to the {n(st['high_draw'])}th sorted draw.",
        f"{s['species']} species give only {s['distinct_species_draws']} distinct draws ({s['ordered_species_draws']} "
        f"in order)",
        f"Over all {s['ordered_species_draws']} ordered draws the same two quantiles give the same reading in "
        f"{s['exact_readings_agree']} of {s['readings']} readings.",
        f"The experiment shows a cost when at least {st['min_learners_with_cost']} of {s['learners']} learners show "
        f"one and no learner reads better on recall or log loss, on at least {st['min_datasets']} datasets in "
        f"{st['min_species']} species.",
        f"Past minus random. Each cell gives the median over {s['datasets']} datasets, the interval and the reading.",
        f"raises the share of past genomes in training from {100 * float(mix['random']['past_share']):.2f}% in the "
        f"random arm to {100 * float(mix['study-out']['past_share']):.2f}% in the study-out arm.",
        f"Test genomes of the {s['datasets']} datasets added up per learner and arm",
        f"scm's share is {100 * float(scm['clip_share']):.1f}% at one decimal, as a right call at the clip adds about "
        f"{st['clip']} to the sum.",
        f"The test sets hold {n(s['tested_pairs'])} genome and dataset pairs, {n(s['resistant_pairs'])} of them "
        f"resistant, from {n(t['distinct_test_genomes'])} distinct genomes. A genome tested in several datasets of its "
        f"species is one pair per dataset and one distinct genome. A test set holds {s['test_min']} to "
        f"{s['test_max']} genomes (median {s['test_median']}), and each arm trains on {s['train_min']} to "
        f"{n(s['train_max'])} genomes (median {s['train_median']}), the same number in all three arms of a dataset.",
        f"The tested genomes come from {t['distinct_test_studies']} distinct studies. The per-species counts below add "
        f"up to {t['sum_of_per_species_study_counts']}, as {t['studies_tested_in_two_species']} studies have tested "
        f"genomes in two species. A test set holds genomes of {s['studies_with_a_tested_genome_min']} to "
        f"{s['studies_with_a_tested_genome_max']} studies ({s['drawn_test_studies_min']} to "
        f"{s['drawn_test_studies_max']} drawn), and {s['datasets_with_one_tested_study']} datasets hold one. In "
        f"{s['datasets_with_a_drawn_study_untested']} datasets a drawn study has no tested genome",
        f"the {t['result_files']} fold files ran again unchanged.",
        f"{s['review_fits_changed']} of the {s['review_fits']} fits changed a score ({n(s['review_scores_changed'])} "
        f"of {n(s['review_scores'])} scores), in {t['result_files_changed']} of the {t['result_files']} result files "
        f"and {t['datasets_with_a_changed_fit']} of the {s['datasets']} datasets. {s['review_readings_changed']} of "
        f"the {s['review_readings']} readings changed, none of them on the total contrast, and the rule's answer "
        f"stayed {s['answer']}.",
        f"the test sets hold genomes of {t['distinct_test_studies']} distinct studies where the first report said "
        f"{t['sum_of_per_species_study_counts']}, which counted each species' studies apart, and "
        f"{s['studies_with_a_tested_genome_min']} to {s['studies_with_a_tested_genome_max']} studies per dataset "
        f"where it said {s['drawn_test_studies_min']} to {s['drawn_test_studies_max']}.",
        f"One draw per arm, one cut per dataset and {s['species']} species:",
        f"holds the record: {s['record_files']} files that `scripts/time-split/record.py` writes",
        f"`checks.tsv` gives {dict((r['check'], r['count']) for r in lines.read_tsv(RECORD / 'checks.tsv'))['genome_ids_of_the_inputs_in_the_record']} "
        f"as the count of the input set's genome ids found in the record.",
    ]


class Page(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = PAGE.read_text(encoding="utf-8")
        cls.rec = lines.Record()

    def test_first_paragraph(self):
        paragraphs = [p for p in self.text.split("\n\n") if p.strip()]
        self.assertEqual(paragraphs[0], "# First-public-day split")
        first = paragraphs[1]
        self.assertIn(LIMIT, first)
        self.assertEqual(self.text.count(LIMIT), 1)
        self.assertTrue(first.startswith("This page compares two ways of choosing training genomes. The past arm trains "
                                         "only on genomes public in BV-BRC by each dataset's cut day. The random arm "
                                         "trains on a random draw of the same size and class balance. Both are scored on "
                                         "the same later genomes of studies new after that day"))

    def test_sections_in_order(self):
        body = outside_code(self.text)
        self.assertEqual([x[3:] for x in body if x.startswith("## ")], SECTIONS)
        self.assertEqual([x for x in body if x.startswith("# ")], ["# First-public-day split"])
        self.assertFalse([x for x in body if x.startswith("###")])

    def test_reading_is_rebuilt_from_the_record(self):
        body = section(self.text, "Reading").strip()
        self.assertEqual(body.split("\n\n"), lines.reading(self.rec))

    def test_reading_has_at_most_four_sentences(self):
        body = section(self.text, "Reading").strip()
        sentences = [x for p in body.split("\n\n") for x in re.split(r"(?<=\.)\s+(?=[A-Za-z])", p)]
        self.assertLessEqual(len(sentences), 4)
        for x in sentences:
            self.assertTrue(x.endswith("."), x)

    def test_reading_gives_each_effect_per_learner_with_its_denominator(self):
        body = section(self.text, "Reading")
        s = self.rec.summary
        for m in re.finditer(r"for (\d+) of (\d+)", body):
            self.assertEqual(m.group(2), s["learners"], m.group(0))
        self.assertIn(f"In {s['datasets']} datasets of {s['species']} species", body)
        for l in ("scm", "cart", "rf", "gbt", "lasso"):
            self.assertIn(f" for {l}", body)
        self.assertIn(f"of {n(s['tested_pairs'])} genome and dataset pairs", body)
        self.assertIn(f"of {n(s['resistant_pairs'])} resistant pairs", body)
        self.assertIn(f"the interval resamples only the {s['species']} species", body)

    def test_study_and_time_paragraph_is_rebuilt_from_the_record(self):
        body = section(self.text, "Study and time contrasts")
        built = lines.study_time(self.rec)
        for para in built:
            self.assertEqual(blocks(body).count(para), 1, para)
        at = blocks(body).index(built[0])
        self.assertEqual(blocks(body)[at:at + len(built)], built)

    def test_every_paragraph_and_list_item_is_built_or_held(self):
        known = set(lines.reading(self.rec)) | set(lines.study_time(self.rec)) | set(HELD)
        self.assertEqual([b for b in blocks(self.text) if b not in known], [])

    def test_held_lines_stand_once_and_in_order(self):
        self.assertEqual(len(HELD), len(set(HELD)))
        found = blocks(self.text)
        for line in HELD:
            self.assertEqual(found.count(line), 1, line)
        self.assertEqual([b for b in found if b in set(HELD)], HELD)
        built = set(lines.reading(self.rec)) | set(lines.study_time(self.rec))
        self.assertFalse(built & set(HELD), "a held line is also rebuilt")

    def test_every_table_is_built_from_the_record(self):
        built = lines.tables(self.rec)
        found = page_tables(self.text)
        self.assertEqual(len(found), len(built))
        for name, table in built.items():
            self.assertEqual(sum(1 for t in found if t == table), 1, name)
            self.assertIn("\n" + "\n".join(table) + "\n", section(self.text, TABLE_SECTION[name]), name)

    def test_sentences_that_give_the_records_counts(self):
        prose = " ".join(x for x in outside_code(self.text) if not x.startswith("|"))
        prose = prose.replace("- **", "").replace("**", "")
        for want in counted(self.rec):
            self.assertEqual(prose.count(want), 1, want)

    def test_every_number_is_in_the_record(self):
        found, _ = lines.values(RECORD)
        missing = [(x, line[:50]) for line in outside_code(self.text) for x in lines.numbers(line) if x not in found]
        self.assertEqual(missing, [])

    def test_hashes_are_the_records(self):
        st = self.rec.settings
        order = section(self.text, "Plan and order")
        self.assertIn(f"The plan (sha256 `{st['plan_sha256'][:8]}`) was hashed with its analysis script (sha256 "
                      f"`{st['analysis_sha256'][:8]}`) before any fold file, fit or result of the experiment existed, by "
                      "its hash list's own note.", order)
        self.assertIn(f"The rerun's plan (sha256 `{st['rerun_plan_sha256'][:8]}`)", order)
        _, hashes = lines.values(RECORD)
        for h in re.findall(r"`([0-9a-f]{6,64})`", self.text):
            self.assertTrue(len(h) >= 8 and any(x.startswith(h) for x in hashes), h)

    def test_the_review_names_the_readings_it_changed(self):
        changed = {(r["learner"], x) for r in self.rec.review for x in r["changed_readings"].split(";") if x}
        self.assertEqual(set(CHANGED), changed)
        for (learner, name), word in CHANGED.items():
            contrast, measure = name.split(" ")
            self.assertEqual(self.rec.readings[(learner, contrast, measure)]["reading"], word, (learner, name))
        self.assertEqual(self.rec.summary["review_total_contrast_readings_changed"], "0")
        self.assertIn("cart's and gbt's study recall, gbt's study specificity and cart's time log loss now read not "
                      "shown; lasso's time log loss and balanced accuracy now read worse.", self.text)

    def test_what_the_plan_fixed_and_what_is_exploratory(self):
        plan = {(r["contrast"], r["measure"]): r["plan"] for r in lines.read_tsv(RECORD / "readings.tsv")}
        self.assertEqual({k for k, v in plan.items() if v == "rule"}, {("total", "recall"), ("total", "log_loss")})
        self.assertEqual({k[1] for k, v in plan.items() if v == "supporting"}, {"specificity", "balanced_accuracy"})
        for phrase in ("It fixed the design, the rule that admits a dataset, recall and log loss per learner on three "
                       "contrasts, the species interval and the rule for a cost.",
                       "The plan named specificity and balanced accuracy as supporting figures and gave them no rule.",
                       "The plan gave the study and time contrasts no rule.",
                       "**Exploratory.** Every figure the plan did not fix is exploratory and enters no rule:",
                       "Exploratory, in the review's narrower form:",
                       "The last column is exploratory:",
                       "scm scores only 0 or 1, so its log loss is about 9.21 times its share of wrong calls, the "
                       "calls its recall and specificity count.",
                       "total is past minus random, the plan's question."):
            self.assertIn(phrase, self.text)

    def test_order_follows_the_record(self):
        order = [r["event"] for r in self.rec.order]
        self.assertEqual(order, list(lines.ORDER_LABELS))
        times = [r["time"] for r in self.rec.order]
        self.assertEqual(times, sorted(times))

    def test_arms_named_as_the_record_names_them(self):
        self.assertEqual([r["arm"] for r in self.rec.mix], list(ARMS))
        self.assertEqual(self.rec.mix[0]["from_later_other_studies"], "0")
        self.assertEqual(self.rec.mix[0]["from_later_test_studies"], "0")
        self.assertEqual(self.rec.mix[2]["from_later_test_studies"], "0")


if __name__ == "__main__":
    unittest.main()
