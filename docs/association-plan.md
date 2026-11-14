# Association plan

This plan audits association evidence for known catalog elements and the stability of released rules. It does not search for, rank or prioritize a new resistance determinant or mutation, and it proposes no engineering intervention. A pattern that passes and that no known element explains is counted. Its sequence is not listed, annotated, ranked or followed up.

A k-mer pattern or catalog family in this plan is a marker with an association to the reported lab call. Annotation does not establish mechanism. Every threshold here is an adapted research screen. Only public BV-BRC and NCBI data are used.

The plan was fixed on 14 Nov 2026, before any association run on a revision 14 cohort. Nothing has been run under it. A run follows it as written, records its sha256 before reading any association statistic, and lists every deviation. Anything decided after a statistic is read is labelled exploratory. A change to the plan is a dated section added at its end, before the run it applies to reads a statistic, and the text above it stays as it is.

The first run under it reruns the October 2026 question on *E. coli* and tobramycin. That run does not reproduce the October numbers. The October row gives a cohort of 350 genomes and names no revision, no genome list of it survives, and no copy of its result file has been found (`docs/work-program-20261027.md`, entry of 2026-11-14 13:13 UTC).

## What the plan fixes

| Item | Fixed value for *E. coli* and tobramycin |
|---|---|
| Cohort | The 1,070 genomes of the revision 14 report `escherichia_coli__tobramycin`, 343 resistant and 727 susceptible |
| Unit tested | One distinct presence pattern of 31-mers over those 1,070 genomes, in at least 5 and at most 1,065 of them |
| Adjustment | A linear mixed model with a kinship matrix built from the patterns that match no catalog gene |
| Adequacy | Inflation factor against permuted labels, a permuted-label count of false passes, a known-element control, two sensitivity scans |
| Number of tests | 2,288,621 |
| Threshold | p below 0.05 / 2,288,621 = 2.184721716701892e-08, family-wise over one scan |
| Conditioning | The one catalog group that agrees best with the top pattern, then every aminoglycoside group, then every catalog group |
| Negative control | 200 label permutations inside lineage clusters, and 200 free shuffles |
| Positive control | The catalog group `aac(6')-Ib-cr`, and synthetic markers of fixed size |
| Replication | A second set of genomes from other studies with at least 25 genomes per class. None exists today |
| Nothing passes | A fixed table. No lower threshold, no second tier, no list of near misses |

## What was known when the plan was fixed

The plan is not blind to this cohort, and a reader should know what its author had seen.

- **Seen.** The audit of 13 Nov 2026 crossed three indicators that the revision 14 report already publishes with the lab call: the catalog rule's call, the families named in the first learned rule term, and their union. It found the rule in 158 of 343 resistant and 11 of 727 susceptible genomes. The choice of the positive control was made knowing that.
- **Seen.** Class counts by lineage cluster, the cohort's size, the count of patterns, and how many genomes hold each catalog group. Section 14 lists the commands.
- **Not seen.** No pattern of the 2,288,621 has been tested against the lab call. No kinship matrix, inflation factor, permutation or conditional scan exists for this cohort.

So the scan, the permutations and the conditional scans are fixed before read. The known-element tests of section 6 are fixed in form, but their direction is known from the audit, and the plan says so wherever it matters.

## 1. Cohort and source revision

| Item | Value |
|---|---|
| Release | `data-rev14-274-20261112`, evaluation revision 14, report schema 2, 274 reports. Its release notes give: built 2026-11-12 from commit `469a4f92266211cd0fd792018029f3d196f0d1a7` (tree `036215156abfae04c1a6d495e7888a5d4c355551`). Its `index.json` has sha256 `1e549f017271b127633ef0f75c7eb9bd88c01c6b361a875f82fe0b4f73ff9fbb` |
| Report | `datasets/escherichia_coli__tobramycin.json`, sha256 `002d36af44b494be71494ea73fedf7ba7e3e0fe91dd816128b6dd550b2ee63e9` |
| Input fingerprint | `9bb9b068c7d6a47c` |
| Genomes | 1,070, all with a FASTA file in the public data root |
| Genome list | The report's `samples[].id`, sorted, one per line. sha256 `604331f43da982c8990101838fdee9f0013e3cb76dd34eb14ef58585de4ff3c3` |
| Labelled list | The same lines as `id<TAB>R` or `id<TAB>S`. sha256 `3330ab75940a8af96e606258e435cb142fa27c5befa847acfedb862d56921ebb` |
| Label | The report's `samples[].resistant`: the reported R or S call. 343 resistant, 727 susceptible. None comes from a MIC (`labels_from_mic` 0) |
| Catalog | AMRFinderPlus reference gene catalog, version 2026-08-07.1. `catalog/AMR_CDS.fa` sha256 `bcc57417c5a170fb4a664de6898455990fcc296dbfe214b6874d761dad63962c` |

Rules for the cohort:

- No genome is added or dropped, for any reason found after the run starts. The sample order of every vector is the report's `samples` order.
- The two pairs of k-mer-identical genomes (`kset_duplicates`, sample indices 35 and 36, 193 and 407) stay in. In the kinship matrix a genome of such a pair is as close to its twin as to itself.
- A run on another dataset or revision needs its own table like this one, added to this plan as a dated section before that run reads a statistic.

**Why.** The October cohort cannot be rebuilt, and revision 14's is the cohort the site serves. A cohort named by a report hash and a list hash cannot be changed after the results are seen. Source: Nosek et al. 2018 (PNAS 115:2600), abstract: "An effective solution is to define the research questions and analysis plan before observing the research outcomes".

The cohort is a mix of collections. The 13 Nov audit found that country alone accounts for 0.43 of the label's variance and lineage cluster for 0.18. The report clusters the genomes into 65 lineage clusters (Mash distance threshold 0.01), 21 of them single genomes. The five largest hold 215, 158, 144, 118 and 49 genomes with 113, 40, 8, 46 and 21 resistant. 26 clusters hold both classes: 968 genomes, 340 resistant and 628 susceptible. The other 39 clusters hold 102 genomes of one class each.

## 2. Unit tested

The unit is a **feature** as `Dataset::project` builds it (`src/matrix.rs`): one distinct presence pattern over the cohort's genomes.

- The species matrix covers the 2,000 *E. coli* genomes of the data root's `species/escherichia_coli/genomes.txt` (sha256 `480ccb8db3a4f93a32257c35c46a41cb75c5895c1be2906fc161e8120f4a0fb7`). All 2,000 have a FASTA file there, and the cohort's 1,070 are among them. It collapses 31-mers that are in the same genomes into one pattern.
- Projected on the cohort, patterns that agree on all 1,070 genomes become one feature.
- A feature is kept when it is in at least 5 and at most 1,065 genomes. The floor is `max(n / 200, 3)` (`src/analysis.rs`), which is 5 for 1,070 genomes.
- A pattern and its complement are two features, as in the report. Counting both is conservative.

The report counts 2,288,621 features. They cover 3,020,033 patterns of the species matrix and 33,587,686 k-mers. The feature count is the number of distinct tests on this cohort.

**Building it.** `grm count` and `grm matrix` at k = 31 build the species matrix from the FASTA files. `grm matrix` must keep k-mers down to a support of 5 or lower, and its default of 2 does. `grm` has no command that writes a dataset's features, so the run calls the library from a small program of its own: `FeatureMatrix::load`, then `Dataset::project` with the report's samples as columns, in the report's order, and a `min_support` of 5. The program writes each feature as one row of bits and no k-mer sequence.

**Gate G1.** Before any statistic, the run shows that the matrix it built is the report's:

- its counts equal the report's `counts.kmers`, `counts.patterns` and `counts.features`;
- each of the 74 features the report describes is one of the run's features. `features[].bits` is base64 of 64-bit words written low byte first, and sample i is bit i mod 64 of word i div 64.

If either fails, the run stops. A count that is close is a different matrix.

**Why.** Patterns that mark the same genomes cannot be told apart by any test on those genomes, so they are one hypothesis. Source: Earle et al. 2016 (Nat Microbiol 1:16041), Methods: "Because the phenotypic contribution of multiple variants with identical phylogenetic patterns cannot be disentangled statistically, we found that pooling such variants improved the power by demanding a less conservative Bonferroni correction than correcting for the total number of variants".

## 3. The tests

All tests regress the label y (1 resistant, 0 susceptible) on one feature x (1 present, 0 absent) and on covariates C. C always holds the intercept. D is the number of columns of C, n is 1,070 and N = n − D.

**M0, unadjusted.** Ordinary least squares.

**M1, the mixed model (primary).** y = Cα + xβ + g + e, with g ~ N(0, σ²h²K) and e ~ N(0, σ²(1 − h²)I). K is the kinship matrix below. This is the FaST-LMM form that pyseer 1.4.2 runs for `--lmm`, and the plan takes its statistic as that source computes it (`pyseer/lmm.py`, `pyseer/fastlmm/lmm_cov.py` at commit `227c0b4247a80156075328b292cf5baea9d1b2f1`):

1. Scale K so that its diagonal sums to n.
2. Regress C out of y, of x and of both sides of K. Call the results ỹ, x̃ and K̃.
3. Take an orthonormal basis U of the space orthogonal to C in which K̃ is diagonal, with eigenvalues S. U has N columns. The source gets it by adding the identity to K before step 2, taking every eigenvector, dropping the D with the smallest eigenvalues and subtracting 1 from the rest. That keeps U orthogonal to C when K̃ has more zero eigenvalues than D, and this cohort's two identical pairs give it two more.
4. For a value h², write a·b for Σᵢ (Uᵀa)ᵢ (Uᵀb)ᵢ / (h²Sᵢ + 1 − h²).
5. Fit h² once, without x, by minimising ½ [Σᵢ log(h²Sᵢ + 1 − h²) + N (log(2π ỹ·ỹ / N) + 1)] over 0 ≤ h² ≤ 0.99999. The source evaluates ten values on a grid and then refines the minimum.
6. For each feature: β̂ = x̃·ỹ / x̃·x̃, r = ỹ·ỹ − β̂ (x̃·ỹ), and F = β̂² (N − 1) x̃·x̃ / r.
7. The p-value is the upper tail of F with 1 and N − 1 degrees of freedom. With the intercept alone, N − 1 is 1,068.

M0 is the same computation with h² = 0, where a·b is the plain inner product.

**Kinship matrix K.** Z holds one column per feature that matches no catalog gene (section 6), standardised: (x − p) / √(p(1 − p)) with p the feature's share of genomes. K = ZZᵀ / m_K, with m_K the number of such features. K is built from presence patterns alone and its sha256 is recorded before the labels are joined to it.

**Why a mixed model.** The October record's claims are about "the mixed model", and the model is the common choice for bacterial k-mer association.

- Lees et al. 2020 (mBio 11:e01344-20): "in exchange for reduced power, the linear mixed model consistently showed the best control of false positives in all settings, always <5%".
- Lees et al. 2018 (Bioinformatics 34:4310): "We have used the FaST-LMM implementation, which allows association in linear time (Lippert et al., 2011), using a kinship matrix estimated from a subset of variants or from a phylogeny."

**Why h² is fitted once.** It keeps one scan to one pass, and it is what the tool does. It is conservative for a strong marker. Source: Zhou and Stephens 2012 (Nat Genet 44:821): "using estimates of variance components or fitted residuals from the null model for testing may be expected to yield conservative p values, leading to a potential loss of power".

**Why K leaves out catalog-matching features.** A marker inside its own kinship matrix is partly explained away before it is tested.

- Listgarten et al. 2012 (Nat Methods 9:525): "when a given SNP is being tested, that SNP should be excluded from the computation of genetic similarity to avoid using it as a covariate. Including the SNP would make the log likelihood of the null model higher than it should be and lead to deflation of the test statistic and loss of power."
- Yang et al. 2014 (Nat Genet 46:100): "Pitfall: loss in power when candidate marker is included in GRM".
- Earle et al. 2016 built theirs from another variant class than the one tested: "For the analyses of SNPs, genes and kmers, we computed the relatedness matrix K from biallelic SNPs only."

The elements under audit are the catalog's, so those are the features left out. Patterns that travel with a catalog gene on the same mobile element stay in K. Section 7's positive controls measure what that costs.

**Why a linear model on a 0/1 label needs checking.** The cohort's collections differ widely in their resistant share (12% of the 313 United Kingdom genomes and 95% of the 110 from Thailand, in the 13 Nov audit). Source: Chen et al. 2016 (Am J Hum Genet 98:653), abstract: "LMMs are generally inappropriate for analyzing binary traits when population stratification leads to violation of the LMM's constant-residual variance assumption". Lees et al. 2020 add: "in highly clonal settings, the LMM can suffer from a very high false-positive rate unless the kinship matrix is carefully chosen". So M1's p-values are not taken on trust. Section 4 checks them on this cohort.

## 4. How the adjustment's adequacy is checked

**Inflation factor.** λ is the median chi-square statistic with 1 degree of freedom that corresponds to the p-values, divided by the chi-square median (0.4549). The run reports λ for every scan, over all features and over the **common features**: those in at least 54 and at most 1,016 genomes (5% of 1,070, rounded up, on each side). The checks below use the common features, because a pattern in a handful of genomes gives a statistic with few possible values, and its median is not the chi-square median. That restriction is this plan's own reasoning and has no source.

Source for λ: Price et al. 2010 (Nat Rev Genet 11:459): λ is "the median χ2 (1 dof) association statistic across SNPs divided by its theoretical median under the null distribution", and "Values of λGC <1.05 are generally considered benign, although inflation in λGC is proportional to sample size."

**λ alone settles nothing, in either direction.**

- A λ above 1 can come from real linked signal. Yang et al. 2011 (Eur J Hum Genet 19:807), abstract: "in the absence of population structure and other technical artefacts, but in the presence of polygenic inheritance, substantial genomic inflation is expected".
- A λ at or below 1 can come from a kinship matrix that absorbs the signal. Listgarten et al. 2012, Table 1, give λ 0.97 with 6 false negatives for the method that kept the tested markers in the matrix, against 1.09 with 2 for the one that removed them. The 13 Nov audit saw the same on simulated data: λ 0.84 with the marker found in 37% of replicates.

So adequacy is a four-part rule, and section 8 uses it for claim 4:

| Part | Passes when |
|---|---|
| A1 | Under M1, the median λ of the 200 lineage-cluster permutations (section 7) is between 0.95 and 1.05 |
| A2 | Under M1, at most 15 of those 200 permutations have any feature below the threshold |
| A3 | Under M1 on the real labels, λ lies between the 2.5th and 97.5th percentile of the 200 permutation values of A1 |
| A4 | The known-element positive control passes under M1 (section 7) |

**Two sensitivity scans.** They are reported beside M1 and decide no claim.

- **S1, lineage clusters as fixed effects.** M0 with one indicator for each of the 44 lineage clusters of two or more genomes. The 21 single-genome clusters share the reference level.
- **S2, country as fixed effects.** M1 with one indicator for each `country` value held by 25 or more genomes: United Kingdom (313), Ireland (180), USA (163), Thailand (110), blank (106), Kenya (45), India (45) and Mali (31). The other 77 genomes are the reference level. The report has no study field, so country is the nearest field to a collection.

For each, the run reports λ, the number of features that pass, and how many of M1's passing features also pass.

## 5. Number of tests and threshold

- **Family.** The 2,288,621 features of this cohort under one model. One scan is one family.
- **Rule.** Bonferroni at a family-wise 0.05. A feature **passes** when p < 0.05 / 2,288,621 = 2.184721716701892e-08. The threshold is an adapted research screen.
- **Every scan uses the same number.** M0, M1, the conditional scans and the sensitivity scans all divide by 2,288,621. Conditioning does not shrink the family.
- **M1 on the real labels is the primary scan.** The other scans are comparisons and are labelled as such. They are not added to the family.
- **The count is checked.** Gate G1 holds the run to 2,288,621. A run cannot count its tests after seeing them.

**Why Bonferroni over unique patterns.** It is the field's rule, and it is conservative here because patterns are correlated.

- Earle et al. 2016: a variant "was considered significant if its P value was smaller than α/n_p, where we took α = 0.05 to be the genome-wide false-positive rate and n_p to be the number of SNPs and genes, or kmers, with unique phylogenetic patterns".
- Lees et al. 2020: "Bonferroni corrections were calculated by using the number of unique unitig patterns as the number of tests."
- Johnson et al. 2010 (BMC Genomics 11:724): "Bonferroni adjustments can be made but are conservative due to the preponderance of linkage disequilibrium (LD) between genetic markers", and "Permutation tests are the gold standard for identifying appropriate significance thresholds".

**The permutation threshold is reported beside it.** The 10th smallest of the 200 smallest p-values from the lineage-cluster permutations estimates the 5% family-wise threshold for this cohort. The run reports it and the number of features that pass it. No claim uses it.

**Design floor.** This is arithmetic on the class counts, not a result. With 343 resistant and 727 susceptible genomes, an exact two-sided test with no adjustment passes the threshold only from these tables up:

| Susceptible genomes with the marker | Fewest resistant genomes with it | Share of the 343 |
|---:|---:|---:|
| 0 | 16 | 4.7% |
| 1 | 18 | 5.2% |
| 5 | 24 | 7.0% |
| 10 | 30 | 8.7% |
| 20 | 40 | 11.7% |

A marker in fewer resistant genomes than that cannot pass, whatever it marks. The mixed model needs more.

## 6. Known catalog elements and conditioning

**At pattern level.** A feature **matches the catalog** when at least one of its 31-mers, on either strand, occurs in a sequence of the pinned `catalog/AMR_CDS.fa`. This is the lookup GRM's own reader of the catalog uses (`Catalog::genes_containing`, `src/catalog.rs`). The run applies it to every k-mer of a feature (`Dataset::patterns_of`, then `FeatureMatrix::kmers_of`). The report's own annotation reads at most 400 k-mers of a feature (`src/annotate.rs`), so the run's flag can be wider.

**Gate G2.** Each of the 21 described features to which the report gives a catalog family must match the catalog in the run.

**At genome level.** The run takes the known-element indicators from the report and computes none of its own:

| Indicator | Source in the report | Count |
|---|---|---:|
| Each catalog group in at least 5 and at most 1,065 genomes | `samples[].groups` | 74 |
| Rare aminoglycoside groups, pooled: any of the 21 aminoglycoside-class groups in fewer than 5 genomes | the same | 1 (35 genomes) |
| Rare groups of other classes, pooled: any of the 38 | the same | 1 (53 genomes) |
| The catalog rule's call | `samples[].catalog_resistant` | 1 (169 genomes) |

The report lists 133 catalog groups in this cohort: 44 of class AMINOGLYCOSIDE (23 in 5 or more genomes) and 89 of other classes (51 in 5 or more). Two of the 77 indicators mark the same genomes (`mph(B) (Mph(B) family)` and `rdmCL`). They are one test and one covariate, and every table gives both names. That leaves 76 distinct indicators.

**Known-element table.** Each of the 76 is tested under M0 and M1 as if it were a feature. The table gives its genome count, the sign of β̂, both p-values, and whether it passes 0.05 / 76 = 6.58e-04 and whether it passes the scan threshold. These are tests of elements the catalog already lists, and the audit saw the rule's row before this plan. They are an audit of known elements and no discovery.

**Conditional scans.** Each is M1 with indicators added to C, so K is rebuilt in the residual space as section 3 describes.

| Scan | Covariates added | Answers |
|---|---|---|
| C1 | The one catalog group g* chosen in section 8 for claim 2 | The October statement as written |
| C2 | The 23 aminoglycoside groups and the pooled rare aminoglycoside indicator | What the drug class's known elements leave |
| C3 | All 74 groups and both pooled indicators, an identical pair entered once | What every known catalog element leaves |

The rule's call is in no conditional scan, because it is the union of groups that are already there. An indicator that is a linear combination of earlier ones is dropped, and the run reports how many were kept.

**Why.** Conditioning on a known marker is the standard way to ask whether a second signal is its own. Source: Yang et al. 2012 (Nat Genet 44:369): "Conditional analysis has been used as a tool to identify secondary association signals at a locus, involving association analysis conditioning on the primary associated SNP at the locus to test whether there are any other SNPs significantly associated."

**What a clean conditional scan shows.** It shows that the passing patterns and the covariate mark nearly the same genomes. The 13 Nov audit measured this on simulated data: a covariate that agreed with the marker on 95% of samples left a median of 177 of 200 linked patterns passing, and one that matched it left 0. So each conditional scan reports, beside its count, the share of the top feature's variance that its covariates explain. For C1 that is the squared correlation with g*, and the run adds the number of genomes where the two differ.

**Classes of a passing feature.** Each feature that passes M1 falls in one class:

| Class | Rule | Reported |
|---|---|---|
| A | It matches the catalog | Count, and the count per catalog family. Family names are the catalog's own, and a feature that matches several families is counted under each |
| B | It does not match the catalog, and it does not pass in C3 | Count only |
| C | It does not match the catalog, and it passes in C3 | Count only |

Class C is "a pattern that passes and that no known element explains". Section 11 binds what the run may do with it.

## 7. Controls

**Negative control: permuted labels.** Two sets of 200, both scanned under M0 and M1 with the real K, the real features and the same threshold. h² is refitted for each permuted label.

- **Lineage-cluster permutations.** Labels are shuffled among the genomes of one lineage cluster, so each cluster keeps its class counts. This keeps the between-lineage structure of the label and breaks any link to a marker inside a lineage. The cohort has about 10^212.6 such labelings.
- **Free shuffles.** Labels are shuffled over all 1,070 genomes. This removes all structure. It checks the code path and feeds claim 1.

Drawing them: numpy 2.5.4, `rng = numpy.random.default_rng([20261114, s, b])` with s = 1 for lineage-cluster permutations and s = 2 for free shuffles, and b from 0 to 199. Inside a cluster, taken in ascending cluster number with its genomes in sample order, the new labels are `y[idx][rng.permutation(len(idx))]`. A free shuffle is `y[rng.permutation(n)]`. The 400 label vectors are written and hashed before any statistic is read.

For each permutation the run keeps the smallest p-value, the number of features below the threshold, and λ.

**Reading it.** With a true family-wise rate of 0.05, the chance of more than 15 false-pass permutations in 200 is 0.044. The control flags a true rate of 0.11 with probability 0.93 and a rate of 0.07 with probability 0.33. So it can show a gross failure and cannot certify 0.05.

**Why inside lineage clusters.** A free shuffle is not a valid null for a structured cohort. Source: Abney 2015 (Genet Epidemiol 39:249): "naive application of permutations may result in an invalid test", and "Only when all individuals are equally related, as in an F2 cross, will a naive permutation approach obtain the correct type 1 error rate." The restricted shuffle assumes that genomes are exchangeable inside a cluster. Structure below the cluster, such as one hospital's clone or one country's collection, is broken by it, so the control says nothing about confounding at that level. S1 and S2 are the plan's only reading of that. Source for what a negative control can show: Lipsitch et al. 2010 (Epidemiology 21:383): the failure of a negative control "does not identify what form of bias is operating."

**Positive control 1: a known aminoglycoside-modifying family.** The catalog group `aac(6')-Ib-cr` (subclass AMIKACIN/KANAMYCIN/QUINOLONE/TOBRAMYCIN, in 135 genomes). It was chosen by a rule that reads no label: among the aminoglycoside groups that the tobramycin rule counts, the one in the most genomes. It passes when its indicator under M1 has p below the scan threshold and β̂ above 0.

If it fails, the mixed model has lost an element the catalog lists for this drug, and the run reports every scan count as unreadable for absence.

**Positive control 2: synthetic markers.** These are synthetic and are labelled so wherever they appear. A synthetic marker is a 0/1 vector drawn from the real labels with a fixed size. It is in no kinship matrix.

| Setting | Resistant genomes with it | Susceptible genomes with it |
|---|---:|---:|
| 10% and 0% | 34 | 0 |
| 15% and 0% | 51 | 0 |
| 20% and 1% | 69 | 7 |
| 30% and 5% | 103 | 36 |

Each setting is drawn 20 times in each of two placements: **spread**, drawn from all genomes of each class, and **one lineage**, drawn only from lineage cluster 2 (215 genomes, 113 resistant and 102 susceptible). A draw uses `rng = default_rng([20261114, 3, setting, placement, draw])`, with the setting 0 to 3 in the table's order, the placement 0 for spread and 1 for one lineage, and the draw 0 to 19. It takes the resistant genomes first, `rng.choice(r_idx, k_r, replace=False)`, then the susceptible ones, `rng.choice(s_idx, k_s, replace=False)`, with both index lists in sample order. The run reports, for M0 and M1, the share of the 20 draws that pass. This is the run's statement of what the screen can find, and it is reported whether or not anything real passes.

**Why known elements as the control.** The field checks an association method by whether it recovers what is known. Source: Earle et al. 2016: "Genuine resistance-conferring variants were detected in all but one study, demonstrating that the high accuracy attained in predicting antimicrobial resistance phenotypes from genotypes known from the literature is mirrored by good power to map the genotypes that confer antimicrobial resistance phenotypes using GWASs."

## 8. The four October claims

Each claim gets one of three words by the rule below: **confirmed**, **not confirmed**, or **not testable on this cohort**. A piece of the plan that did not run makes every claim that needs it not testable. A confirmed claim is confirmed on revision 14's cohort. No word here is a statement about the lost October files.

**Claim 1. "The earlier calibration failure was a printing-precision artefact."**

On the 200 free shuffles under M1, over the common features, compute two summaries: L, the median of the 200 values of λ, and T, the pooled share of p-values at or below 0.001. The calibration check passes when 0.95 ≤ L ≤ 1.05 and 0.0008 ≤ T ≤ 0.00125. Compute it on the 64-bit values, then again after writing every p-value as text and reading it back, in four formats: `%.17g`, `%.2E` (the format pyseer prints), `%.6f` and `%.3f`.

- Confirmed: the check passes on the 64-bit values and fails in at least one text format.
- Not confirmed: it fails on the 64-bit values, or it passes in every format.

Much of this is arithmetic and the run says so: `%.3f` turns the share at or below 0.001 into the share below 0.0015. What the run adds is whether full-precision calibration holds on this cohort. Which check failed in October is not recorded.

**Claim 2. "The strongest signal is explained by a known catalog gene family that the current tobramycin rule does not count."**

- T is the feature with the smallest p under M1. A tie goes to the feature in more genomes, then to the smaller pattern read as bytes.
- g* is the catalog group, among the 74, whose indicator has the largest squared correlation with T. A tie goes to the first by name.
- T is **explained** by g* when T does not pass in scan C1.
- g* is **not counted** when the report's `counts_for_rule` is false for it. The audit found that this flag comes from a group's first hit, so the run also prints the group's subclass.

Confirmed: T passes under M1, T is explained by g*, and g* is not counted. Not confirmed: T passes, and it is not explained or g* is counted. Not testable: no feature passes under M1.

The run also reports whether T matches the catalog and which families, the squared correlation of T with g*, and the number of genomes where they differ.

**Claim 3. "Conditioning on it removes all Bonferroni-significant patterns."**

Confirmed: no feature passes in C1. Not confirmed: at least one does. Not testable: claim 2 is not testable. The counts for C2 and C3 are reported beside it, by class.

**Claim 4. "Genomic control is adequate with the mixed model."**

Confirmed: A1, A2, A3 and A4 of section 4 all pass. Not confirmed: any of them fails, and the run names which. λ under M0 and M1 is printed beside the word, for the real labels and for both permutation sets. If only A3 fails, and on the high side, the run also prints λ of scan C3, since linked signal from known elements raises λ too.

## 9. Released rules

The report describes 74 features that its final models use, 21 of them with a catalog family. The scan tests them with every other feature, so this section adds no test.

The run reports how many of the 74 pass under M1, split by whether the feature matches the catalog, and how many still pass in S1 and S2. A described feature that matches the catalog is named by its catalog families. One that does not is counted.

These features were chosen by models fitted on these same 1,070 genomes, so a pass is after selection. It shows that a rule term keeps its association when lineage structure is modelled. It is not a new finding and says nothing about another cohort. Source: Kraft et al. 2009 (Stat Sci 24:561): "studying the marker in additional samples (beyond those used to initially identify the marker) will help produce more unbiased estimates of the genetic effect".

## 10. Replication

A result of this plan is **replicated** only when all of these hold:

1. A second set of genomes shares no genome with the cohort, and none of its studies (BioProject, or PubMed id where there is no BioProject) has a genome in the cohort.
2. Its labels follow the cohort's rule: a reported, unconflicted R or S call, none derived from a MIC.
3. It holds at least 25 genomes of each class, the project's own minimum for a dataset.
4. Its genome list and its sha256 are added to this plan before any of its statistics is read.
5. The same model M1 is run on it, with a kinship matrix built from that set alone.
6. Only known catalog elements are carried forward: the genome-level indicators that passed the scan threshold in the cohort. Each is tested once. It replicates when β̂ has the same sign and p is below 0.05 divided by the number carried forward.

Features of class B and C are never carried forward. A split of the 1,070 genomes is not replication, and neither is a set labelled from MICs.

**Today no such set exists.**

- The curated table of the public data root holds 488 *E. coli* genomes outside the cohort with a reported, unconflicted tobramycin call: 22 resistant and 466 susceptible. None of their FASTA files is in the data root.
- The frozen population of the outside-study test (`scripts/outside-study/frozen/datasets.tsv`) has 17 resistant and 60 susceptible genomes in scope for this dataset and marks it not eligible.
- Another 2,654 genomes outside the cohort have a label derived from a MIC (551 resistant, 2,103 susceptible). That is a different label rule.

So a run on the cohort alone is reported as one cohort, read once, not replicated.

**Why.** Source: Kraft et al. 2009 (Stat Sci 24:561): "replication across different designs and populations provides qualitative evidence that the association is not due to uncontrolled bias affecting a single study", and "Use the same analytic methods".

## 11. What the run may write about a pattern

This is a rule of the run, not a preference.

- **Class A.** The count, and the count per catalog family, with the family's catalog name, class and subclass. A catalog family named in a result is a marker. The run does not say what the gene does.
- **Classes B and C.** Counts only. The run's report, tables, figures, logs and output files give no k-mer, no sequence, no contig, no position, no gene or product name, no genome id and no feature number for any such feature.
- **No ranking.** The run publishes no ordered list of features and no list of features near the threshold.
- **No follow-up.** The run maps, aligns, annotates and searches nothing for a class B or C feature, and proposes no experiment, edit or construct for any feature.
- **Working files.** The per-feature p-values stay in the run's scratch, keyed by feature number and holding no sequence, so the counts can be recomputed. They are not copied to the run's output directory, a report, a repository or a release.
- **Wording.** A passing feature is "a pattern associated with the reported lab call in this cohort". It is never a determinant, a cause or a mechanism.

## 12. What is reported when nothing passes

If no feature passes under M1, the run reports this and nothing softer:

| Line | Content |
|---|---|
| Passing features | 0 of 2,288,621 at 2.184721716701892e-08 |
| Smallest p under M1 | The value. No feature is identified |
| Known-element table | As in section 6, complete |
| Positive control 1 | Passed or failed, with its p |
| Positive control 2 | The share of synthetic draws found, by setting and placement |
| Negative control | The count of false-pass permutations of 200, and the permutation threshold |
| λ | Under M0 and M1, real labels and permutations |
| Claims | Claims 2 and 3 are not testable. Claims 1 and 4 by their rules |

The threshold is not lowered. No "suggestive" tier, false discovery rate or top list is added afterwards. The sentence that goes with the table is: "No pattern passed an adapted research screen on this cohort. The screen cannot see a marker below the sizes in the synthetic-marker table, so this is weak evidence of absence."

If positive control 1 failed, the sentence is instead: "The screen did not find a catalog element known for this drug, so its empty result cannot be read."

**Why.** A threshold or a tier chosen after an empty result is a choice made with the result in hand. Source: Nosek et al. 2018, abstract: "Mistaking generation of postdictions with testing of predictions reduces the credibility of research findings." Johnson et al. 2010, quoted in section 5, is the source for reading an empty Bonferroni scan as weak evidence.

## 13. Precision, order of work and deviations

**Precision.**

- p-values are computed and compared in 64-bit floating point. The run also keeps log10 p from the log of the upper tail, and uses it wherever p underflows.
- Stored copies are 64-bit binary or text with 17 significant digits (`%.17g`). No decision is made on a printed value.

**Order of work.**

1. Copy this plan, record its sha256, add the run's specifics, and record the sha256 of that file.
2. Label-free steps: build the matrix, pass G1 and G2, flag catalog matches, build K and hash it, write the indicator table and hash it, draw and hash the 400 label vectors and the synthetic markers.
3. Check the implementation (below).
4. Only then the scans, in any order.

**Checking the implementation.** The run's code is checked against a second computation before step 4, on the 76 indicators and on 2,000 features: `default_rng([20261114, 4]).choice(2288621, 2000, replace=False)`, read as feature numbers in the order `Dataset::project` returns them.

- Preferred: pyseer 1.4.2 with `--lmm`, the same K as `--similarity`, the check rows as a `--pres` table named by number, the labels as `--phenotypes`, and no covariate or lineage option. Its defaults drop a row in fewer than 1% or more than 99% of genomes (`--min-af 0.01`, `--max-af 0.99`), so the run sets `--min-af 0.004 --max-af 0.996`. It also drops a row whose p-value is 1 (`--filter-pvalue 1`, `--lrt-pvalue 1`). A row the tool leaves out is listed with the tool's note and is not compared. For every other row the run's p, printed as `%.2E`, must equal the tool's `lrt-pvalue` or differ by one unit in the last digit.
- Fallback, if that release cannot be installed: generalised least squares in statsmodels 0.15.0 with covariance h²K̃ + (1 − h²)I at the run's h². p-values must agree to a relative 1e-8. The fallback is a deviation and is listed.

If the check fails, the run stops. It does not pick the implementation that gives the nicer result.

**Deviations.** Each is listed with its reason and with whether it was decided before or after a statistic was read. A scan that did not run is listed as not run, and section 8 then makes its claims not testable.

**What every run reports.** Whatever passes, the run's report holds these and each has a number or the words "not run":

| Group | Content |
|---|---|
| Gates | G1 and G2 with the counts compared, the sha256 of K, of the indicator table and of the 400 label vectors, and the implementation check |
| Scans | For M0, M1, C1, C2, C3, S1 and S2: h², λ over all features and over common features, the number of common features, and the number of features that pass |
| Classes | Under M1: the counts of classes A, B and C, and class A's count per catalog family |
| Known elements | The table of section 6 |
| Controls | Section 7: both permutation sets under M0 and M1, positive control 1, and the synthetic-marker shares |
| Released rules | Section 9 |
| Claims | One of the three words for each of the four claims, with the numbers its rule reads |
| Replication | "Not replicated: one cohort", until a set of section 10 exists |
| Deviations | The list |

## 14. Tools and the commands behind this plan

| Tool | Pin |
|---|---|
| `grm` | The repository at a commit whose `src/`, `Cargo.toml`, `Cargo.lock` and `.cargo/` equal the release commit's. `main` at `698db80` does. Rust 1.98.1 |
| Python | 3.12 |
| numpy | 2.5.4, wheel `numpy-2.5.4-cp312-cp312-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl`, sha256 `fbde6962867ee75b48b0ee29b2b9372ec5d617799dbaf38e82dc0596f2f7738a` |
| scipy | 1.18.1, wheel `scipy-1.18.1-cp312-cp312-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl`, sha256 `f55fa87b6c612ecd6b058f167c53231b1d14e412efe361d3d6e38b3631c73218` |
| statsmodels (fallback check) | 0.15.0, wheel `statsmodels-0.15.0-cp312-cp312-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl`, sha256 `c6f2bd48915cc9b3b95c3c7f6ac363b6f8aaa38172962ba5f7f6885f0dfe57f2` |
| pyseer (preferred check) | 1.4.2, tag commit `227c0b4247a80156075328b292cf5baea9d1b2f1`. The Python package index holds no file for it. The bioconda channel holds `noarch/pyseer-1.4.2-pyhdfd78af_0.conda`, sha256 `167125aad0ebf13efc14bc388716a977767843ea23daa72e567055f5e994bbf2` |

The hashes are the package index's and the channel's own, read on 14 Nov 2026. Nothing was installed to write this plan. The run records the versions and hashes it installed. If pyseer's own requirements force another numpy or scipy in its environment, pyseer gets an environment of its own, and the run's code keeps the pins above.

The numbers in this plan come from the report named in section 1 and from the public data root, by two scripts that read label-free fields and class counts and test nothing against the label:

```bash
R=rev14-274-20261112/datasets/escherichia_coli__tobramycin.json
python3 -I cohort_facts.py $R results    # cohort, list hashes, lineage clusters, counts, catalog groups
python3 -I plan_numbers.py $R results    # indicators, threshold, control cut, design floor, synthetic sizes
```

Both scripts, their outputs and the fetched source texts are kept with the plan's lane record, outside this repository.

## 15. Sources

Each passage quoted above was read in the full text or, where marked, the abstract, as PubMed Central serves it. Reference numbers and citations inside a passage are left out, except where the quote shows one.

| Choice | Source | Read |
|---|---|---|
| Plan fixed before outcomes | Nosek et al. 2018, PNAS 115:2600, PMC5856500 | abstract |
| Unit and test count: unique patterns | Earle et al. 2016, Nat Microbiol 1:16041, PMC5049680; Lees et al. 2020, mBio 11:e01344-20, PMC7343994 | full text |
| Mixed model for bacterial k-mer association | Lees et al. 2018, Bioinformatics 34:4310, PMC6289128; Lees et al. 2020 | full text |
| The statistic | pyseer 1.4.2 source, `pyseer/lmm.py` and `pyseer/fastlmm/lmm_cov.py`, commit `227c0b4` | source |
| Variance ratio fitted once | Zhou and Stephens 2012, Nat Genet 44:821, PMC3386377 | full text |
| Tested elements left out of the kinship matrix | Listgarten et al. 2012, Nat Methods 9:525, PMC3597090; Yang et al. 2014, Nat Genet 46:100, PMC3989144; Earle et al. 2016 | full text |
| A linear model on a 0/1 label under stratification | Chen et al. 2016, Am J Hum Genet 98:653, PMC4833218 | abstract |
| Inflation factor and its limits | Price et al. 2010, Nat Rev Genet 11:459, PMC2975875 (full text); Yang et al. 2011, Eur J Hum Genet 19:807, PMC3137506 (abstract) | as marked |
| Bonferroni is conservative; permutation as the reference | Johnson et al. 2010, BMC Genomics 11:724, PMC3023815 | full text |
| Conditional analysis | Yang et al. 2012, Nat Genet 44:369, PMC3593158 | full text |
| Permutation needs exchangeability | Abney 2015, Genet Epidemiol 39:249, PMC4634896 | full text |
| What a negative control shows | Lipsitch et al. 2010, Epidemiology 21:383, PMC3053408 | full text |
| Known elements as a positive control | Earle et al. 2016 | full text |
| Replication, and a pass after selection | Kraft et al. 2009, Stat Sci 24:561, PMC2865141 | full text |
| Nothing chosen after the result | Nosek et al. 2018 | abstract |
| The cohort's structure, the simulated limits of conditioning and of λ | The audit of 13 Nov 2026 of the October record, a project record outside this repository | full text |

## 16. Limits

- The plan was written by someone who had read the audit's crossing of the catalog rule with the lab call. The known-element tests are therefore not blind, and the plan treats them as an audit.
- Lineage clusters at a Mash distance of 0.01 are coarse: five clusters hold 684 of the 1,070 genomes. A kinship matrix adjusts more finely, and the permutations do not test that finer adjustment.
- Country stands in for collection. 106 genomes have none.
- The labels are reported R and S calls from many laboratories. 391 of the 1,070 genomes have a MIC.
- h² is fitted once for each scan, and the statistic treats a 0/1 label as continuous. Both are the tool's choices, checked here by permutation and not by theory.
- 200 permutations cannot certify a 0.05 error rate. They can show a gross failure.
- The restriction of λ to common features is the plan's own choice. So are the lower bound of 0.95 for λ, the band of 0.0008 to 0.00125 for the share of p-values at or below 0.001, the cut of 15 in 200 permutations, the four synthetic-marker sizes and the floor of 25 genomes for a country level. Only the upper bound of 1.05 has a source (Price et al. 2010).
- A pattern in fewer than 16 resistant genomes cannot pass. An empty scan says little about markers of that size.
- One species, one drug, one cohort. No replication set exists today.
- Bonferroni over 2,288,621 correlated tests is conservative, so a pattern that does not pass has not been shown to be unassociated.
