# What a rule term's words say

A set covering machine rule term is the presence or absence of one feature: a set of k-mers found in exactly the same genomes. The site prints each term in words. This page says where those words are written, what each part of them claims, and what none of them claims.

A term names a marker. Its k-mers go with the lab calls in the dataset's isolates. A gene name or a catalog family in the words says where the k-mers lie or what they match. It does not say why they predict, and annotation does not establish mechanism.

## Where the words are written

| Words | Written by | Changes when |
| --- | --- | --- |
| The Rule column of the datasets table, and its search text | `grm report`, into `scm_rules` of `index.json` (`report::rule_text` in `src/report.rs`) | a new release is built |
| The note under a rule that names several catalog families | the web app, from the index's words | the app is deployed |
| Each rule's heading and notes in the dataset inspector, and the rule slide's heading and notes | the web app, from the dataset's report (`web/src/data/derive/dataset.ts`) | the app is deployed |
| The within-lineage sentences under the inspector's rules | the web app, from `web/src/data/lineage-control.json` | the app is deployed, for the release the table was computed on |

A site that serves an older release keeps that release's index words until the next release, whatever app it runs.

The index's new words wait for a release under a new name. `grm report` stamps every index it writes with the time it ran (`generated_at`), so an index rebuilt in place from the same reports would differ from the published one in more than its words. The site's "Report built" date would move with no report changed, and `web/scripts/browser-check.cjs` would find no `reference/release-summary/*/guidance.json` for the tree, because it matches that file by the index's generation time. Do not rebuild the index of a served release to pick up a wording.

## The parts of a term's words

The index writes a term as up to three clauses, joined by `; `. The app's `indexRuleWords` builds the same words, and `tests/data/rule-text-cases.json` holds both to the same synthetic cases.

1. **The name.** `presence of` or `absence of`, then the first of these that applies:
   - the catalog family the k-mers match, or `k-mers matching N catalog families (…)` when they match several;
   - `k-mers N bp from the nearest annotated protein-coding gene, NAME` when the k-mers lie in no annotated protein-coding gene of the reference contig;
   - the gene, or else the product, the k-mers lie in;
   - `feature N` when nothing is named.

   "Gene" here means a protein-coding gene. `grm` reads only the coding genes of the reference genome's annotation (`fetch_features` in `src/bvbrc.rs` asks for CDS features), so a term worded by its distance can still lie in an rRNA, tRNA or other RNA gene, and the words claim nothing about those. For the same reason the inspector gives such a term's place as `outside protein-coding genes, N bp from NAME`, not as intergenic.
2. **The drug's rule.** For a term that names catalog families: which of their gene groups GRM's known-gene rule for this drug uses and which it does not, for example `GRM's known-gene rule for this drug uses blaNDM but not ble`, or `GRM's known-gene rule for this drug does not use armA`. GRM's rule is the subject and the gene group is its object. The rule reads the catalog's class and subclass for a group, so the clause describes GRM's own rule and says nothing of what the gene does. The clause reads the report's own gene-group rows, which exist only for groups GRM's catalog reader found in the dataset's genomes. A family with no row is left out, and a report without a catalog block gets no clause. The words then say nothing; they do not guess. One family can have a group with a subfamily and one without, and the second then carries the family's bare name; it is written `the group named aac(6')-Ib alone`, so that one name does not stand on both sides of the clause. The inspector and the rule slide give the same thing as a sentence: `Of the gene groups GRM found in these genomes, its known-gene rule for this drug does not use armA.`
3. **Ties.** `one of N features that separate these isolates equally well` when other features have the same covered and lost counts at the term's greedy step. For a later term the isolates are the ones the earlier rules left undecided, and the words say so.

A list of several families is one presence pattern. A feature matches a family when one of its sampled k-mers lies in any allele of it, so the k-mers may lie in any one of the families or span them. The list does not say which gene an isolate carries, and the app prints that under such a rule.

## Wording versions

The wording has changed without a change to the models, so the index numbers it in `rule_text_version` (`report::RULE_TEXT_VERSION`).

| Version | Words | Releases |
| --- | --- | --- |
| 1 | A term names one catalog family, the first in byte order | through evaluation revision 12; the index has no version field |
| 2 | A term names every catalog family its feature matches, after their count | `rev14-274-20261112`; the index has no version field |
| 3 | A term also gives its place outside a protein-coding gene, the drug's rule and its ties, as above | the next release |

A count of rule terms by the kind of gene they name is read from these words, so it holds for one wording only. Text that gives such a count names the evaluation revision and the wording version it read. The Methods page prints both under the set covering machine. Where the index has no version field the page prints `rule wording 2 or earlier: the index gives no number`, because two releases of one evaluation revision can differ in wording and the revision alone does not tell them apart.

The version is a number for the words only. It is not part of the evaluation: `EVALUATION_REVISION`, the report schema and the index schema do not change with it, and no dataset report changes.

## The within-lineage sentences

Under the inspector's rules the app keeps the caveat that a rule's k-mers can travel with the cause without being the cause, and adds one of three sentences for each rule when the table holds a verdict for this exact report:

- `the association is still seen when isolates are compared only inside their lineage cluster`
- `the test did not show the association inside lineages. That is not a finding of no association.`
- `too few lineage clusters hold both lab calls and both states of the rule, so the test could not be run`

After them the app says how the test was run and what it does not show. The counts in that text are read from the table's own counts (`lineageTested` in `web/src/data/derive/lineage.ts`), so they move with the table and are typed nowhere:

- One test pools the clusters. It does not show the association in each one, so the first sentence does not say the association holds in every lineage.
- How many of the release's rules could be tested, how many of those passed, and the level.
- The clusters are MinHash clusters, not a phylogeny.
- The rules were chosen on the same isolates, so the test is not held-out evidence.
- An association seen inside lineage clusters is still an association, not a cause.

`reference/lineage-control/README.md` gives the test, its counts and its limits. The second sentence says the test did not reach its level; it does not say the association is absent. The third sentence is not a negative result.
