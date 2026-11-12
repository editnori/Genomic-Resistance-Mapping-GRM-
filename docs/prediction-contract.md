# Prediction contract (grm-predict-1)

What a published dataset report's final model means on a genome outside its cohort, and what `grm predict` does with it. It scores genomes with a model exactly as `grm learn` fitted it: nothing is refitted, retuned or recalibrated. Every prediction row names this contract as `"contract": "grm-predict-1"`.

## The published model

A dataset report (schema 2, `REPORT_SCHEMA` in `src/analysis.rs:36`) holds the final models `grm learn` fits on every sample of the cohort (`final_models`, `src/learn/eval.rs:1244`). Their feature ids are rows of the dataset matrix: `final_plan` (`src/learn/eval.rs:1304`) requires that every row survives `training_candidates` on all samples, so the final models' candidates are the dataset matrix itself. A dataset row is every species-matrix pattern that is identical over the dataset's genomes (`Dataset::project`, `src/matrix.rs:1654`), and so every k-mer with one presence vector over the cohort.

The report lists some features under `features` (`model_features`, `src/analysis.rs:255`). Each entry has:

- `feature`, the dataset row;
- `n_kmers`, how many k-mers the row covers;
- `kmer`, one of them: the first k-mer placed on the reference contig, read on the contig's strand (`src/annotate.rs:365`), or, when no reference places the feature, the first of its k-mers in canonical form (`src/annotate.rs:173`). Either way it is a k-mer of the row;
- `bits`, the row's presence over the samples in report order: base64 of the little-endian u64 words, bit i for sample i (`base64`, `src/analysis.rs:232`; written at `src/analysis.rs:748`).

The served copy under `report/datasets/` is the file `grm learn` wrote, byte for byte (`src/report.rs:612`), so either file identifies the same model.

### Which learners are complete

| Learner | Published as | Complete? |
|---|---|---|
| SCM | `models.scm`: `kind` and `rules` (feature, presence) | Always. `model_features` lists every rule's feature. |
| CART | `models.cart.nodes` (value, split feature, children) and `models.cart_depth` | Always. `model_features` lists every split's feature, and the tree was grown to `cart_depth`. |
| Lasso | `models.lasso`: `intercept`, nonzero `coef`, `lambda` | Only when every nonzero coefficient's feature is listed. `model_features` lists the 30 largest coefficients by magnitude (`src/analysis.rs:270`), plus any feature listed for another learner. |
| Forest | `models.rf_importance`, the top 50 importances (`src/learn/eval.rs:1293`) | Never: no tree is published. |
| Boosted trees | `models.gbt_importance`, the top 50 importances | Never: no tree is published. |

Checked on 2026-11-11 against the 274 served revision-12 reports: every SCM rule and CART split names a listed feature with a valid 31-mer in all 274; the lasso is complete in 142 and incomplete in 132; none publishes a forest or boosted model beyond importances. `grm predict` on each report, with any genome, prints the learners it leaves out and why, which reproduces the lasso count. Decoding the stored `bits` of their own features, the published SCM rules reproduce every greedy step's `covered` and `lost`, and the published CART tree every node's `n` and `n_resistant`, in all 274 (a read-only script outside the tree; changing its bit order breaks the replay in 273 of them).

## A feature's value outside the cohort

Inside the cohort every k-mer of a row has the same presence vector, by construction. On a new genome they need not agree: the genome can carry some k-mers of a row and not others, for example part of a gene, or a variant at one position of a SNP pattern. The candidate definitions are the published k-mer alone, any k-mer of the row, every k-mer of the row, or a share of them. On the cohort they all give the row's bits, so the cohort cannot say which is right.

**Rule.** A feature is present on a genome when the genome carries the feature's published `kmer`, on either strand. A k-mer is found when it occurs in the genome's distinct canonical k-mers at the report's k, counted as `grm count` counts an assembly (`kset::count_file`, `src/kset.rs:28`, with `KmerSpace::for_each_canonical`, `src/dna.rs:103`). Bases other than A, C, G and T break the window (`src/dna.rs:109`), so a k-mer across an N or an IUPAC code is not found. Case is ignored.

**Why this one.** The published k-mer is the only sequence the report commits to. Any other rule needs the species matrix files (`kmers.bin` and the rest), which are not published, so it would make the served report alone unscorable, and it would make the model's input depend on files outside its identity. The published k-mer is also one fixed, exact sequence, so the value does not depend on a threshold.

**Disagreement is counted, not hidden.** When the species matrix the report was learned from is given (`--matrix`), every feature in a row also records `kmers_found`, how many of the row's `n_kmers` k-mers the genome carries, and the row records `partial_features`, the features with some but not all of their k-mers found. Those are exactly the features where the published k-mer, any k-mer and every k-mer can give different values. The model still reads the published k-mer's value.

`grm predict` accepts a matrix only if it rebuilds the report's dataset: the same k, every sample a column, the report's feature and pattern counts, the same `bits` for every model feature, `n_kmers` k-mers per feature and the published k-mer among them.

## A scorable genome

- An assembly, as FASTA (plain or gzip), or a kset `grm count` wrote from one at the report's k. A FASTQ file (reads) is refused: read errors would add k-mers no assembly holds, and `grm count` counts assemblies with every k-mer kept.
- At least one k-mer of length k.
- The genome is of the report's species and would pass the QC `grm curate` applies to cohort genomes (`QcConfig::default`, `src/curate.rs:192`): at most 1,000 contigs, CheckM completeness at least 90 and contamination at most 10, and a length within 25% of the species' median; and the species identity check (docs/species-check.md). `grm predict` cannot check CheckM scores, the species median or the species, so it does not decide this. For a FASTA input, each row records the assembly's length, contigs, GC, N50 and ambiguous bases, and whoever scores a genome applies the QC and the species check before reading its rows.

## A prediction row

One JSON line per genome and learner:

| Field | Meaning |
|---|---|
| `contract` | `grm-predict-1` |
| `dataset` | the report's `id` |
| `report_sha256` | SHA-256 of the report file's bytes: the model's identity |
| `input_fingerprint` | the report's own input fingerprint |
| `evaluation_revision` | the report's evaluation revision |
| `k` | the report's k |
| `genome` | the kset's name, or the FASTA file name without `.gz` and its last extension |
| `source` | the file read |
| `assembly` | length, contigs, gc, n50, ambiguous, for a FASTA input; null for a kset |
| `learner` | `scm`, `cart` or `lasso` |
| `call` | `resistant` or `susceptible`; null for a lasso with no coefficient |
| `score` | CART's leaf value, or the lasso's probability; null for the SCM and for a lasso with no coefficient |
| `features` | per feature the learner reads: `feature`, `kmer`, `found`, `n_kmers`, `kmers_found` (null without a matrix) |
| `partial_features` | features found in part; null without a matrix |

The SCM's call is its conjunction or disjunction of presence and absence terms (`ScmModel::predict`, `src/learn/scm.rs:80`). CART's score is the value of the leaf the genome reaches with nodes at `cart_depth` taken as leaves (`Tree::predict`, `src/learn/tree.rs:167`). The lasso's score is the logistic function of the intercept plus the coefficients of the features found (`LassoModel::predict_proba`, `src/learn/lasso.rs:104`); a lasso with no coefficient is the null fit (`lasso_scores`, `src/learn/eval.rs:447`), and `grm predict` makes no call for it: `call` and `score` are null. A score at or above 0.5 is a resistant call, the cutoff evaluation scores every learner with (`src/learn/eval.rs:756`). `grm predict` uses these functions themselves on a one-genome matrix of feature values, so a learner with a fitted model cannot be applied differently from how it was evaluated.

The null fit is the one case where the two differ. From evaluation revision 14, a cross-validation fold whose model is null scores each of its genomes at the share of resistant genomes in the fold's training set (`pool_folds`, `src/learn/eval.rs:728`), so that every learner is scored on the same genomes (docs/equal-scoring.md). That is a rule for scoring folds. It reads nothing from the genome, so `grm predict` does not apply it to a final model, at any revision.

## Refused

- A report whose `schema` is not 2, whose own evaluation audit is not valid, or whose k is outside 5..=32.
- A learner whose published form is not complete: a rule, split (above `cart_depth`) or nonzero coefficient on a feature the report does not list, or whose listed `kmer` is not k letters of A, C, G and T. Without `--learners`, such a learner is left out with its reason on standard error; named in `--learners`, it is refused. A report with no complete learner is refused.
- A kset counted at another k, a FASTQ file, a genome with no k-mer.
- A matrix that does not rebuild the report's dataset (above).

Rows are written only after every genome is scored, so a refusal leaves no partial output.

## Not claimed

- A call outside the cohort is what the frozen model gives for the published k-mers found and not found. It is a marker-level statement: it says which k-mers were found, not that any of them causes resistance.
- No calibration outside the cohort. Scores and calls carry the cross-validated metrics of their report only as a description of that cohort; on genomes from other studies, error rates are unmeasured until such genomes with laboratory phenotypes are scored.
- No clinical use.
- The forest and boosted models are not published, so they are not scored. A forest or boosted model refitted from the cohort is a new model, not the published one, and a row from it would need its own identity.
- The published k-mer was chosen by annotation for readability, not for robustness across studies; another k-mer of the same row could be found where it is not. `kmers_found` shows how often that happens; it does not change the call.

## Replay check

On the cohort the rule must give back the report's own model. `predict_replays_the_final_models_on_their_own_cohort` (`src/main.rs` tests) runs `grm count`, `structure`, `matrix` and `learn` on a synthetic species of 60 genomes, then `grm predict` on every cohort genome's FASTA with the matrix. Every feature's values over the cohort equal its published `bits`, every genome gets the final SCM's, CART's and lasso's fitted calls and scores, no feature is found in part, and the published SCM uses a presence and an absence term on a block that half its carriers hold reverse-complemented. The unit tests in `src/predict.rs` cover a feature found in part, a k-mer broken by a non-ACGT byte, both strands, and each refusal.
