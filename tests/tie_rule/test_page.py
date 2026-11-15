"""docs/tie-rule.md and reference/tie-rule/README.md against the committed record.

Every line of the page that carries a number, outside code spans and a few
fixed phrases, is written again here from the record's files and must be a
line of the page, and every line the suite writes must be on the page once.
So a figure typed wrong, a figure left behind by a new record and a number
with no cell behind it all fail. The page's references into the code are
read at the tree's own files, and the words the brief rules out are looked
for. No count is typed here: each comes from the record.
"""

import importlib.util
import pathlib
import re
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "reference" / "tie-rule"
PAGE = ROOT / "docs" / "tie-rule.md"


def load_record():
    spec = importlib.util.spec_from_file_location("tie_rule_record", ROOT / "scripts" / "tie-rule" / "record.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = load_record()

# Where the present rule lives at the tree's commit: file, line, the text that line starts with.
CODE = (
    ("src/learn/tree.rs", 72, "pub fn fit_tree("),
    ("src/learn/tree.rs", 126, "let better = |a: (f64, u32), b: (f64, u32)| {"),
    ("src/learn/tree.rs", 308, "pub fn fit_forest("),
    ("src/learn/gbt.rs", 228, "fn grow("),
    ("src/learn/gbt.rs", 259, ".reduce_with(|a, b| {"),
)
COMPARISON = "if a.0 > b.0 || (a.0 == b.0 && a.1 < b.1) {"
# Phrases that carry a number and no figure of the record.
FIXED = ("evaluation revision 14", "15 November 2026", "64-bit", "95%", "Python 3.12")
NUMBER = re.compile(r"(?<![\w.])[+-]?\d[\d,]*(?:\.\d+)?(?:e[+-]?\d+)?%?(?![\w-])")


def summary():
    return {r["key"]: r["value"] for r in R.read_table(REFERENCE / "summary.tsv")}


def settings():
    return {r["key"]: str(r["value"]) for r in R.read_table(REFERENCE / "settings.tsv")}


def paired():
    return {(r["learner"], r["split"], r["measure"]): r for r in R.read_table(REFERENCE / "paired.tsv")}


def order():
    return {(r["split"], r["first"], r["second"]): r["order_changes"] for r in R.read_table(REFERENCE / "order.tsv")}


def n(x):
    return f"{x:,}"


def pct(a, b):
    return f"{100 * a / b:.1f}%"


def signed(x):
    if x != 0 and abs(x) < 0.00005:
        return f"{x:+.5f}"
    return f"{x:+.4f}"


def names(items):
    items = [f"`{i}`" for i in items]
    return ", ".join(items[:-1]) + " and " + items[-1] if len(items) > 1 else items[0]


def numbered(line):
    """The number tokens of a line outside code spans and the fixed phrases."""
    text = re.sub(r"`[^`]*`", "", line)
    text = re.sub(r"\]\([^)]*\)", "]", text)
    for phrase in FIXED:
        text = text.replace(phrase, "")
    return NUMBER.findall(text)


def page_lines():
    """Every number-bearing line of the page, written from the record."""
    s, p, o = summary(), paired(), order()
    run, release, lineage = s["datasets.run"], s["datasets.release"], s["datasets.lineage"]
    lines = []

    lines.append(f"- cart: `fit_tree` in `src/learn/tree.rs` (line {CODE[0][1]}), whose closure `better` (line {CODE[1][1]}) compares two candidates.")
    lines.append(f"- rf: `fit_forest` (`src/learn/tree.rs`, line {CODE[2][1]}) calls `fit_tree` with the features it draws at each node (`max_features`), so the same closure decides among them.")
    lines.append(f"- gbt: `grow` in `src/learn/gbt.rs` (line {CODE[3][1]}), whose `reduce_with` (line {CODE[4][1]}) compares two candidates the same way.")

    rep = [s[f"reproduction.{k}"] for k in ("main_equals_release", "present_equals_main", "traced_present_equals_present", "traced_exact_equals_exact")]
    assert len(set(rep)) == 1, rep
    lines.append(
        "Each dataset ran in five arms on the release's matrices, folds and seeds: A, `main`'s binary; B, the branch's binary with "
        "the flag off, which is the present rule; C, as B with a trace of every fit; D, the exact rule; E, as D with a trace. With "
        f"the `seconds` fields masked, A equals the released report, B equals A, C equals B and E equals D on {rep[0]} of {run} "
        "datasets (`reproduction.tsv`). The plan was fixed before any run (`PLAN.md`)."
    )
    lines.append(
        f"{run} of the {release} datasets and {s['species.run']} of the {s['species.release']} species were run. The plan takes the "
        "species in ascending order of the release's summed fitting time, and its stop rule left the "
        f"{s['species.not_run']} species with the largest sums unrun: {names(s['species.not_run.names'].split(','))}, "
        f"{s['datasets.not_run']} datasets. They differ from the datasets run: at the median, {n(s['genomes.median.not_run'])} "
        f"genomes against {n(s['genomes.median.run'])}, and {n(s['features.median.not_run'])} features against "
        f"{n(s['features.median.run'])}, {s['features.median.ratio']:.1f} times as many (`datasets.tsv`). Nothing on this page "
        "speaks for them."
    )
    genomes = {L: (s[f"{L}.lineage.genomes"], s[f"{L}.random.genomes"]) for L in R.LEARNERS}
    assert len(set(genomes.values())) == 1, genomes
    lines.append(
        f"The lineage split exists for {lineage} of the {run} datasets, so its denominators are {lineage}, and the random split's "
        f"are {run}. A held-out call is a score at or above {settings()['call_threshold']}. Each learner scores "
        f"{n(genomes['cart'][0])} held-out genomes on the lineage split and {n(genomes['cart'][1])} on the random split."
    )

    either = {L: s[f"{L}.either.metrics_changed"] for L in R.LEARNERS}
    assert s["trees.either.log_loss_changed"] == run
    lines.append(
        f"Each tree learner's released figures depend on the tie rule on part of the {run} datasets run (cart {either['cart']}, "
        f"rf {either['rf']}, gbt {either['gbt']}), and some tree learner's log loss changes on all {s['trees.either.log_loss_changed']}."
    )
    for L in R.LEARNERS:
        lines.append(
            f"| {L} | {s[f'{L}.lineage.metrics_changed']} of {lineage} / {s[f'{L}.random.metrics_changed']} of {run} / "
            f"{either[L]} of {run} | {s[f'{L}.either.log_loss_changed']} of {run} | {n(s[f'{L}.lineage.calls_changed'])} / "
            f"{n(s[f'{L}.random.calls_changed'])} | {s[f'{L}.all.recall.largest_abs_change']:.3f} / "
            f"{s[f'{L}.all.balanced_accuracy.largest_abs_change']:.3f} / {s[f'{L}.all.log_loss.largest_abs_change']:.3f} |"
        )
    assert either["scm"] == either["lasso"] == 0
    lines.append(
        f"Over the three tree learners and both splits, {n(s['trees.calls_changed'])} held-out calls change. scm and lasso, which "
        "grow no trees, change on no dataset."
    )
    lines.append(
        f"Datasets whose change, in absolute value, is over {R.SHARES[0]} / over {R.SHARES[1]}, out of {lineage} on the lineage "
        f"split and {run} on the random split:"
    )
    for L in R.TREES:
        for S in R.SPLITS:
            cells = [f"{s[f'{L}.{S}.{m}.over_{R.SHARES[0]}']} / {s[f'{L}.{S}.{m}.over_{R.SHARES[1]}']}" for m in ("recall", "log_loss", "balanced_accuracy")]
            lines.append(f"| {L}, {S} | " + " | ".join(cells) + " |")

    st = settings()
    lines.append(
        "For each tree learner and split, the change of each dataset, exact minus present, is averaged over the datasets, and its "
        f"95% interval comes from a bootstrap over species ({s['species.run']} species, {n(int(st['paired_draws']))} draws, seed "
        f"{st['paired_seed']}; `paired.tsv`)."
    )
    for L in R.TREES:
        for S in R.SPLITS:
            cells = []
            for m in ("recall", "log_loss", "balanced_accuracy", "specificity"):
                r = p[(L, S, m)]
                cells.append(f"{signed(r['mean'])} [{signed(r['low'])}, {signed(r['high'])}]")
            lines.append(f"| {L}, {S} | " + " | ".join(cells) + " |")
    excludes = sorted((L, S, m) for (L, S, m), r in p.items() if r["low"] > 0 or r["high"] < 0)
    for (L, S, m), r in p.items():
        if m in ("balanced_accuracy", "log_loss"):
            assert r["low"] <= 0 <= r["high"], (L, S, m)
    assert excludes == [("cart", "lineage", "recall"), ("cart", "lineage", "specificity"), ("cart", "random", "recall")], excludes
    rec, spe, ll = p[("cart", "lineage", "recall")], p[("cart", "lineage", "specificity")], p[("cart", "lineage", "log_loss")]
    assert rec["low"] > 0 and spe["high"] < 0
    lines.append(
        "A paired test over datasets, with a species bootstrap, finds no difference in balanced accuracy or log loss for any tree "
        f"learner on either split. For cart on the lineage split, the exact rule moves recall up by {rec['mean']:.3f} "
        f"[{rec['low']:.3f}, {rec['high']:.3f}] and specificity down by {-spe['mean']:.3f} [{-spe['high']:.4f}, {-spe['low']:.3f}], "
        f"beside log loss {ll['mean']:.3f} [{ll['low']:+.3f}, {ll['high']:+.3f}]."
    )
    lines.append(
        f"{len(excludes)} of the {len(p)} intervals exclude 0: cart's recall and specificity on the lineage split, and cart's recall "
        f"on the random split. About {0.05 * len(p):.1f} would at 95% by chance. The test shows no difference in balanced accuracy "
        "or log loss; it does not show that the two rules are equally accurate."
    )

    lines.append(
        f"On {s['order.random.best_changes']} of the {run} datasets run, the learner with the highest random-split balanced "
        "accuracy changes when only the tie rule changes; the order of any two learners by balanced accuracy changes on at most "
        f"{s['order.lineage.most_pair_changes']} of {lineage} datasets (lineage split) and {s['order.random.most_pair_changes']} "
        f"of {run} (random split)."
    )
    lines.append(f"| Pair of learners | Lineage split, of {lineage} | Random split, of {run} |")
    for a in range(len(R.LEARNERS)):
        for b in range(a + 1, len(R.LEARNERS)):
            x, y = R.LEARNERS[a], R.LEARNERS[b]
            lines.append(f"| {x}, {y} | {o[('lineage', x, y)]} | {o[('random', x, y)]} |")
    lines.append(
        "The tie effect of a dataset and split is the largest change in balanced accuracy over cart, rf and gbt. Under the present "
        f"rule, the gap from the best learner to the second has a median of {s['order.lineage.median_gap']:.4f} on the lineage "
        f"split and {s['order.random.median_gap']:.4f} on the random split, and the tie effect a median of "
        f"{s['order.lineage.median_tie_effect']:.4f} and {s['order.random.median_tie_effect']:.4f}. The tie effect is above 0 and "
        f"at least the gap on {s['order.lineage.tie_effect_at_least_gap']} of {lineage} datasets (lineage split) and "
        f"{s['order.random.tie_effect_at_least_gap']} of {run} (random split) (`splits.tsv`). There, the difference between the two "
        "best learners is no larger than what the tie rule alone moves."
    )

    t = {k: s[f"trace.{k[0]}.{k[1]}"] for k in ((L, c) for L in R.KINDS for c in (
        "searched_present", "exact_tied_present", "differ_present", "differ_tie_present", "differ_order_present",
        "fits", "fits_differing_present"))}
    for L in ("cart", "rf"):
        assert t[(L, "differ_order_present")] == 0 and t[(L, "differ_tie_present")] == t[(L, "differ_present")], L
    cart_tied = pct(t[("cart", "exact_tied_present")], t[("cart", "searched_present")])
    rf_tied = pct(t[("rf", "exact_tied_present")], t[("rf", "searched_present")])
    lines.append(
        f"At {n(t[('cart', 'exact_tied_present')])} of the {n(t[('cart', 'searched_present')])} cart nodes searched ({cart_tied}) "
        f"and {n(t[('rf', 'exact_tied_present')])} of the {n(t[('rf', 'searched_present')])} rf nodes ({rf_tied}), more than one "
        "candidate shares the best exact impurity. There the exact rule takes the lowest feature number, and the present rule "
        "takes it too unless rounding departs, so under both rules the split depends on the order of the matrix rows. Rounding "
        f"departs at {n(t[('cart', 'differ_present')])} cart nodes "
        f"({pct(t[('cart', 'differ_present')], t[('cart', 'searched_present')])}) and {n(t[('rf', 'differ_present')])} rf nodes "
        f"({pct(t[('rf', 'differ_present')], t[('rf', 'searched_present')])}). Each of them is an exact tie, and at none did "
        "rounding put a strictly worse split first (`traces.tsv`, the present trees)."
    )
    lines.append(
        f"For gbt, the exact rule moves a split only to a lower feature with the same two groups. It moves "
        f"{n(t[('gbt', 'differ_present')])} of the {n(t[('gbt', 'searched_present')])} nodes searched in the logistic models, in "
        f"{n(t[('gbt', 'fits_differing_present')])} of {n(t[('gbt', 'fits')])} fits, and {n(t[('gbt-mic', 'differ_present')])} of "
        f"{n(t[('gbt-mic', 'searched_present')])} in the MIC models, in {t[('gbt-mic', 'fits_differing_present')]} of "
        f"{t[('gbt-mic', 'fits')]} fits (`traces.tsv`)."
    )

    lines.append(
        "- **Keep the present rule and merge the flag.** Needs a Methods sentence that ties are broken on 64-bit gains, so the lower "
        "feature number wins only when two gains are equal to the last bit; tests for the code that passes the rule and the counts "
        "to the fits, where the review planted faults that every test passed and that change real reports; exact reports kept out "
        "of results directories, because the site and the index read neither `tie_rule` nor the input fingerprint; and the branch "
        "binary's lasso time solved first. With the flag off that binary writes `main`'s reports, but took a median of "
        f"{s['seconds.lasso.median_ratio_present_main']:.2f} times `main`'s lasso seconds per dataset, and "
        f"{s['seconds.cart.median_ratio_present_main']:.2f} times its cart seconds, on the run's shared machine (`seconds.tsv`). "
        "The review measured the lasso slowdown again on an idle machine."
    )
    lines.append(
        f"- **Adopt the exact rule at a later revision.** Needs `EVALUATION_REVISION` raised and all {release} datasets fitted "
        f"again, the {s['datasets.not_run']} not run included; a decision on gbt candidates with equal gains and different groups; "
        "the site and the index to show `tie_rule` wherever exact and present reports can be read side by side; and the time cost "
        "measured on an idle machine. On the run's shared machine the exact rule took "
        f"{s['seconds.cart.median_ratio_exact_present']:.2f} (cart), {s['seconds.rf.median_ratio_exact_present']:.2f} (rf) and "
        f"{s['seconds.gbt.median_ratio_exact_present']:.2f} (gbt) times the present rule's seconds, median per dataset."
    )
    lines.append(
        f"- **A rule that does not depend on feature numbers.** Both rules take the lowest number at {cart_tied} of cart and "
        f"{rf_tied} of rf nodes searched. Needs its own design and experiment."
    )
    lines.append(
        f"- {run} of the {release} datasets were run. The {s['datasets.not_run']} datasets of the {s['species.not_run']} largest "
        "species were not, and they have more features."
    )
    return lines


def readme_lines():
    s, st = summary(), settings()
    return [
        f"This folder holds the record behind [docs/tie-rule.md](../../docs/tie-rule.md). The experiment fitted {s['datasets.run']} "
        "datasets of release `rev14-274-20261112` (evaluation revision 14) twice, once with the present rule that the tree "
        "learners use to break ties between equal splits and once with an exact rule, and an independent review recomputed its "
        "figures on the same day. A change under the exact rule is a sensitivity of a figure to the tie rule. It does not show "
        "that either rule is more accurate. Every figure is cross-validated on public BV-BRC genomes.",
        f"A held-out call is a score at or above {st['call_threshold']}. Log loss is the mean over held-out genomes of the "
        f"negative log likelihood of the score clamped to [{st['log_loss_clamp']}, 1 - {st['log_loss_clamp']}]. The best learner "
        "is the index's rule (`src/report.rs`, `best_learner`): the highest balanced accuracy, a tie to the later learner of scm, "
        "cart, rf, gbt and lasso. The trace counts come from the branch's own counters, one line per fit: nodes searched, nodes "
        "with more than one candidate at the best exact impurity (`exact_tied`), and nodes where the exact pick differs from the "
        "present one, as an exact tie (`differ_tie`), a strictly better exact split (`differ_order`) or, for gbt, a lower "
        "candidate with the same two groups (`differ_partition`).",
        f"- {s['datasets.run']} of the {s['datasets.release']} datasets were run. The {s['datasets.not_run']} of the "
        f"{s['species.not_run']} largest species were not, and nothing here speaks for them.",
    ]


def text_lines(path):
    """The lines of a page outside fenced code blocks."""
    out, fenced = [], False
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("```"):
            fenced = not fenced
            continue
        if not fenced:
            out.append(line)
    return out


class Page(unittest.TestCase):
    def test_every_line_with_a_number_is_written_from_the_record(self):
        want = set(page_lines())
        for line in text_lines(PAGE):
            if numbered(line):
                self.assertIn(line, want, f"a line of the page with {numbered(line)} that the record does not give")

    def test_every_line_the_record_gives_is_on_the_page_once(self):
        lines = text_lines(PAGE)
        for line in page_lines():
            self.assertEqual(lines.count(line), 1, line)

    def test_the_code_the_page_cites_is_at_those_lines(self):
        for path, number, start in CODE:
            source = (ROOT / path).read_text(encoding="utf-8").splitlines()
            self.assertTrue(source[number - 1].strip().startswith(start), f"{path}:{number} is not {start!r}")
        # Both comparisons are the present rule: the larger gain, and a bit-equal gain to the lower feature.
        tree = (ROOT / "src/learn/tree.rs").read_text(encoding="utf-8").splitlines()
        gbt = (ROOT / "src/learn/gbt.rs").read_text(encoding="utf-8").splitlines()
        self.assertEqual(tree[CODE[1][1]].strip(), COMPARISON)
        self.assertEqual(gbt[CODE[4][1]].strip(), COMPARISON)
        # fit_forest's trees take a drawn feature set, so the same closure decides among them.
        forest = "\n".join(tree[CODE[2][1] - 1 : CODE[2][1] + 30])
        self.assertIn("fit_tree(", forest)
        self.assertIn("max_features: Some(max_features)", forest)
        # The index's best learner keeps the last of equal maxima, so a tie goes to the later learner.
        report = (ROOT / "src/report.rs").read_text(encoding="utf-8")
        self.assertIn(".max_by(|a, b| a.1.total_cmp(&b.1))", report)

    def test_the_page_names_the_commits_of_the_record(self):
        page = PAGE.read_text(encoding="utf-8")
        st = settings()
        for key in ("main_commit_runs", "branch_commit_runs", "branch_commit_head"):
            self.assertIn(f"`{st[key][:7]}`", page, key)
        self.assertIn(f"branch `{st['branch']}`, which is not on `main`", page)
        self.assertIn("The present rule stays the production rule of evaluation revision 14.", page)
        self.assertIn("recommends none", page)

    def test_the_words(self):
        for path in (PAGE, REFERENCE / "README.md"):
            text = path.read_text(encoding="utf-8")
            self.assertNotIn("—", text, path.name)
            for word in ("carrier", "provenance", "receipt"):
                self.assertNotIn(word, text.lower(), f"{path.name}: {word}")
        page = PAGE.read_text(encoding="utf-8")
        self.assertIn("A change under the exact rule is a sensitivity of a figure to the tie rule.", page)
        self.assertIn("It does not show that either rule is more accurate.", page)
        self.assertNotIn("rule is better", page)
        self.assertNotIn("more accurate rule", page)

    def test_the_readme_lines_with_a_number_are_written_from_the_record(self):
        want = readme_lines()
        lines = text_lines(REFERENCE / "README.md")
        for line in want:
            self.assertEqual(lines.count(line), 1, line)
        for line in lines:
            if numbered(line) and not line.startswith("| "):
                self.assertIn(line, want, f"a README line with {numbered(line)} that the record does not give")


class Wiring(unittest.TestCase):
    def test_the_suite_is_in_ci_and_its_document(self):
        workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
        self.assertTrue("          - suite: tests/tie_rule\n            top: .\n" in workflow, "tests/tie_rule is not in the python job")
        loader = unittest.defaultTestLoader
        count = 0
        for module in sorted(pathlib.Path(__file__).parent.glob("test_*.py")):
            spec = importlib.util.spec_from_file_location(f"tie_rule_count_{module.stem}", module)
            m = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(m)
            count += loader.loadTestsFromModule(m).countTestCases()
        self.assertTrue(f"`tests/tie_rule` {count} |" in (ROOT / "docs" / "ci.md").read_text(encoding="utf-8"),
                        f"docs/ci.md does not give tests/tie_rule its {count} tests")


if __name__ == "__main__":
    unittest.main()
