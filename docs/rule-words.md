# What a rule term's words say

A set covering machine rule term is the presence or absence of one feature: a set of k-mers found in exactly the same genomes. The site prints each term in words. This page says where those words are written, what each part of them claims, and what none of them claims.

A term names a marker. Its k-mers go with the lab calls in the dataset's isolates. A gene name or a catalog family in the words says where the k-mers lie or what they match. It does not say why they predict, and annotation does not establish mechanism.

## Where the words are written

| Words | Written by | Changes when |
| --- | --- | --- |
| The Rule column of the datasets table, and its search text | `grm report`, into `scm_rules` of `index.json` (`report::rule_text` in `src/report.rs`) | a new release is built |
| The note under a rule that names several catalog families | the web app, from the index's words | the app is deployed |
| Each rule's heading and notes in the dataset inspector, and the rule slide | the web app, from the dataset's report (`web/src/data/derive/dataset.ts`) | the app is deployed |
| The within-lineage sentences under the inspector's rules | the web app, from `web/src/data/lineage-control.json` | the app is deployed, for the release the table was computed on |

A site that serves an older release keeps that release's index words until the next release, whatever app it runs.

## The parts of a term's words

The index writes a term as up to three clauses, joined by `; `. The app's `indexRuleWords` builds the same words, and `tests/data/rule-text-cases.json` holds both to the same synthetic cases.

1. **The name.** `presence of` or `absence of`, then the first of these that applies:
   - the catalog family the k-mers match, or `k-mers matching N catalog families (…)` when they match several;
   - `k-mers N bp from the nearest annotated gene, NAME` when the k-mers lie in no annotated gene of the reference contig;
   - the gene, or else the product, the k-mers lie in;
   - `feature N` when nothing is named.
2. **The drug's rule.** For a term that names catalog families: which of their gene groups count for this drug's known-gene rule and which do not, for example `blaNDM counts for this drug's known-gene rule, but ble does not`. The clause reads the report's own gene-group rows, which exist only for groups found in the dataset's genomes. A family with no row is left out, and a report without a catalog block gets no clause. The words then say nothing; they do not guess.
3. **Ties.** `one of N features that separate these isolates equally well` when other features have the same covered and lost counts at the term's greedy step. For a later term the isolates are the ones the earlier rules left undecided, and the words say so.

A list of several families is one presence pattern. A feature matches a family when one of its sampled k-mers lies in any allele of it, so the k-mers may lie in any one of the families or span them. The list does not say which gene an isolate carries, and the app prints that under such a rule.

## Wording versions

The wording has changed without a change to the models, so the index numbers it in `rule_text_version` (`report::RULE_TEXT_VERSION`).

| Version | Words | Releases |
| --- | --- | --- |
| 1 | A term names one catalog family, the first in byte order | through evaluation revision 12; the index has no version field |
| 2 | A term names every catalog family its feature matches, after their count | `rev14-274-20261112`; the index has no version field |
| 3 | A term also gives its place outside a gene, the drug's rule and its ties, as above | the next release |

A count of rule terms by the kind of gene they name is read from these words, so it holds for one wording only. Text that gives such a count names the evaluation revision and the wording version it read. The Methods page prints both under the set covering machine.

The version is a number for the words only. It is not part of the evaluation: `EVALUATION_REVISION`, the report schema and the index schema do not change with it, and no dataset report changes.

## The within-lineage sentences

Under the inspector's rules the app keeps the caveat that a rule's k-mers can travel with the cause without being the cause, and adds one of three sentences for each rule when the table holds a verdict for this exact report:

- `the association also holds inside lineages`
- `the test did not show the association inside lineages`
- `these genomes do not have enough mixed lineages to separate this from lineage`

`reference/lineage-control/README.md` gives the test, its counts and its limits. An association that holds inside lineages is still an association, and the rules were chosen on the same isolates, so the test is not held-out evidence. The second sentence says the test did not reach its level; it does not say the association is absent. The third sentence is not a negative result.
