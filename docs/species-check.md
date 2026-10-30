# Species identity check (revision 9)

Status: implemented on branch `r9/species-check`; awaiting its independent review. Revision 9 includes it only if that review has passed by the time revision 9's rerun starts; otherwise it goes into revision 10. Revision 8 is unchanged.

## Why

A matrix review found two genomes labelled *Mycobacterium tuberculosis* in BV-BRC that are not *M. tuberculosis*: `1773.41701` (GC 64.2%, carries `blaMAB`, likely *M. abscessus*) and `1773.40949` (GC 68.4%, unidentified). Both have good CheckM scores, and both have Jaccard about 0.001 to the rest of the cohort. They are in the frozen MTB 2k cohort, so revisions 7 and 8 trained on them. CheckM scores a genome against the marker set of the lineage it was filed under, so a clean assembly of the wrong species can pass it.

The GRM admin ruling of 2026-10-29 added a species check for all 23 species. It compares each genome's GC and Mash distance with its labelled species' reference, and genomes that fail are dropped before curation. Later rulings that day set the nearer-species rule, its E. coli/Shigella exemption, and the handling of the MTB GC outliers below.

## How it runs

```
grm species-check          # judges each species' genomes.txt, writes the verdicts
grm curate                 # excludes the failed genomes, then selects as before
grm fetch fasta            # downloads replacement genomes the new selection picked
grm species-check          # judges the replacements and keeps every earlier verdict
grm curate                 # repeat the last three steps until nothing new fails
grm species-check --reset  # only to choose references and verdicts afresh
```

- `grm species-check` reads `curated/species.json` and each species' `genomes.txt`, plus every genome the existing `species-check/species_check.json` recorded for that species. Curation drops failed genomes from `genomes.txt`, so a check that read only the selection would no longer see them, and the next `curate` would let them back in. Earlier genomes are judged again from their cached sketches. It sketches every FASTA it finds into `species-check/sketches/<genome_id>.json` (GC, length and sketch, cached) and writes `species-check/species_check.json`. `--sketch-only` fills the cache without judging, so a machine holding some of the FASTA files can sketch them for a check that runs elsewhere.
- A species keeps the reference its first check chose, the medoid of its selection before any drop. Recomputed on later cohorts, the medoid moves once failures leave and replacements arrive: the independent review saw *A. baumannii* move from 470.2428 to 470.2437 and *E. cloacae* from 550.3494 to 550.1141. With the cohort only growing and the reference fixed, every genome keeps its distances, so a failed genome stays failed and `check -> curate` reaches a fixed point after the replacements are judged. Every recorded species is judged on every rerun, including one that has left `curated/species.json`, from its cached sketches. Its reference is the nearer species for the others: without the *K. quasipneumoniae* reference, the 67 *K. pneumoniae* that are *K. quasipneumoniae* would pass. If any recorded reference has no sketch or FASTA file on the machine, for example a species with no sketches at all, the check stops rather than judge without it. A recorded genome with neither a cached sketch nor a FASTA file keeps its earlier verdict.
- A genome recorded as failed stays failed. If a rerun judges it pass or borderline, its verdict stays `fail` with the recorded reason, `kept_fail` records what the rerun judged, and the command prints it. `grm species-check --reset` chooses references afresh and judges every genome without its recorded verdict. It still judges every genome the existing file recorded, because curation has dropped the failed ones from `genomes.txt`. Removing `species_check.json` instead would forget them, and the next `grm curate` would select them again. If a recorded failure has neither a sketch nor a FASTA file, `--reset` cannot judge it, so it stops and names the genome and species. A recorded pass or borderline genome in that state is listed as unchecked.
- `grm curate` reads `species-check/species_check.json` when it exists. Each failed genome is excluded before any other QC, and recorded in `curated/qc.json` with the reason `species check: sequence is another species: <evidence>` and in `curated/funnel.json` under that key. Its phenotypes leave every dataset, so class counts and the dataset gate recompute. `curated/species.json` gains `selected_unchecked` per species: the selected genomes the check has no verdict for.
- Selection is `select_genomes`, unchanged. The fresh selection is deterministic, but it is not the old selection minus the dropped genomes: dropped labels shrink class sizes, which shifts scores. For species above the cap, replacements fill the cap, and they need FASTA before the check can judge them. The admin accepted this on 2026-10-29.

### `species_check.json`

```
{ k: 21, sketch_size: 1000,
  config: { fail_distance: 0.13, species_distance: 0.05, gc_z: 6, gc_points: 1,
            genomic_groups: [["Escherichia coli", "Shigella *"]] },
  species: [{ species, slug, reference, gc_median, gc_sd, checked, failed, borderline,
              unchecked: [ids], distance_p50, distance_p99, distance_max_pass }],
  pairs:   [{ labelled, nearer, genomes, share, same_group }],
  genomes: [{ genome_id, species, verdict: "pass" | "borderline" | "fail", gc, gc_z,
              distance, ani, nearest_species, nearest_distance, inside_species, reason,
              kept_fail? }] }
```

Every checked genome is listed, sorted by species and then genome id. Every species an earlier check recorded appears with counts from this check, never copied from an earlier file. A genome without a sketch keeps its recorded entry unchanged. `kept_fail` is present only on a recorded failure that this check judged otherwise. On such an entry, `reason` and `inside_species` are the recorded check's, and the numeric fields (`gc`, `gc_z`, `distance`, `ani`, `nearest_*`) are this check's.

## Method

**Signal.** Each genome's sketch is the 1,000 smallest mixed canonical 21-mers of its assembly, the same values `KsetReader::sketch` reads from a kset counted at k = 21. k = 21 is Mash's default: at 95% ANI, two genomes still share about a fifth of their 21-mers, so the sketch resolves the species boundary. The Mash distance (`cluster::mash_distance`) approximates 1 - ANI near that boundary. GC is the G+C share of A, C, G and T in the assembly.

**Reference.** Each species' reference is its medoid: the genome whose median Mash distance to the other genomes of its cohort is smallest. Ties go to the lower genome id, and cohorts are sorted by id, so the choice is deterministic. The reference comes from the cohort itself, so no download is needed. A median is unaffected by a few foreign genomes, so contaminants become the reference only if they are most of the cohort. The medoid is chosen once, on the selection before any drop, and later checks keep it (see "How it runs"). The cost is that a species complex is measured from its largest member. For example, the *E. cloacae* reference 550.3494 sits in the cohort's main group, and the other complex members sit 0.085 to 0.111 from it.

**Verdicts.** A genome fails when
1. its distance to its own reference exceeds `fail_distance` = 0.13 (ANI estimate about 87%), or
2. another checked species' reference is nearer than its own and within that species' `species_distance` = 0.05 (about 95% ANI). The sequence then names that species. This rule never applies between two members of one genomic group (below).

A genome is borderline, recorded and never dropped, when it fails neither way but
- is outside its own 0.05 boundary,
- would fail rule 2 against a member of its own genomic group, or
- has GC more than 6 robust standard deviations (1.4826 x MAD, floored at 0.05 points) and at least 1 percentage point off its species median.

GC alone never fails a genome, since contaminating contigs and plasmids also move it. Genomes are never moved to another species in this revision.

**Thresholds.** 0.05 is the conventional 95% ANI species boundary. 0.13 is set from the distribution observed across the 23 species: species complexes reach 0.111 (the *E. cloacae* complex), and the next genome, a *Klebsiella* filed as *E. coli*, is at 0.146. Nothing falls between 0.112 and 0.145. A first cut at 0.10 split the *E. cloacae* complex's 0.08 to 0.12 band and was moved. Complexes that sit past 0.05 but under 0.13, with no other checked species nearer, stay borderline: *E. cloacae* (28 genomes, 0.085 to 0.111), *K. oxytoca* (51, 0.056 to 0.072), *E. faecium* (58, 0.050 to 0.055, likely clade B, *E. lactis*), *P. aeruginosa* (11, 0.050 to 0.055, likely *P. paraeruginosa*) and *A. baumannii* (8, 0.063 to 0.104).

**The E. coli/Shigella genomic group.** *E. coli* and every *Shigella* are one genomic species split by taxonomy. 532 of 2,000 *E. coli* (26.6%) are nearer the *S. sonnei* reference than their own, and 30 *S. sonnei* are nearer the *E. coli* reference. Rule 2 never applies between any two members of the group (`Escherichia coli` and `Shigella *`), in any direction. Rule 1 still applies to each against its own reference. Genomes rule 2 would have failed within the group are borderline, which is why *E. coli* shows 536 borderline genomes (532 plus 4 past its own boundary). The 20 *S. flexneri* nearer the *S. sonnei* reference sit 0.0151 to 0.0159 from it and 0.0193 to 0.0205 from their own reference 623.1307. For comparison, the own distances of all 1,984 *S. sonnei* to 624.3489 are min 0.0000, median 0.0013, p99 0.0201 and max 0.0305. 1,937 of its 1,984 genomes (97.6%) are within 0.015, so the 20 sit outside the bulk of the clone but inside its tail. The 20 also look like 10 duplicate pairs. Those go to the outer-fold duplicate-isolation lane, not this check.

**Coverage for revision 9.** 25,761 genomes checked, none unchecked: the 23,761 genomes of the frozen 23-species candidate (`/root/grm-extra-species-analysis-20261027/data/species/*/genomes.txt`, whose BV-BRC inputs match `frozen-inputs.json` of the MTB 4k pilot), with MTB checked on the 4k list (`grm-mtb-4k-pilot-opus-20261029/run/cohorts/full`, whose first 2,000 lines are the frozen 2k), which adds 2,000 genomes. The laptop sketched 18,909 genomes and Vast sketched the 6,852 held only there (16 threads, nice 10, read-only, 118 MB scratch in `/tmp/grm-r9-species-check`). `species_check.json` sha256 `25df721ee5d902a3a43a44bfb3f15bf67cb75a2ec6cc497cd69625f962fa2d71`.

## Results

| Species | Reference (medoid) | Checked | Failed | Borderline | GC median | Own distance p50 / p99 / max of passing |
|---|---|--:|--:|--:|--:|---|
| *Acinetobacter baumannii* | 470.2428 | 1705 | 3 | 8 | 39.0% | 0.0113 / 0.0226 / 0.0324 |
| *Campylobacter coli* | 195.2067 | 391 | 0 | 0 | 31.3% | 0.0072 / 0.0132 / 0.0210 |
| *Campylobacter jejuni* | 197.16128 | 508 | 0 | 0 | 30.4% | 0.0101 / 0.0212 / 0.0264 |
| *Clostridioides difficile* | 1496.2126 | 204 | 0 | 0 | 28.6% | 0.0072 / 0.0341 / 0.0342 |
| *Corynebacterium diphtheriae* | 1717.869 | 243 | 0 | 0 | 53.6% | 0.0144 / 0.0210 / 0.0216 |
| *Enterobacter cloacae* | 550.3494 | 271 | 1 | 28 | 55.1% | 0.0235 / 0.0391 / 0.0435 |
| *Enterococcus faecium* | 1352.10050 | 2000 | 0 | 58 | 37.7% | 0.0073 / 0.0487 / 0.0498 |
| *Escherichia coli* | 562.57275 | 2000 | 4 | 536 | 50.6% | 0.0277 / 0.0341 / 0.0371 |
| *Klebsiella aerogenes* | 548.511 | 93 | 0 | 0 | 55.0% | 0.0115 / 0.0333 / 0.0367 |
| *Klebsiella oxytoca* | 571.1931 | 112 | 0 | 51 | 55.2% | 0.0085 / 0.0110 / 0.0112 |
| *Klebsiella pneumoniae* | 573.31067 | 2000 | 68 | 2 | 57.2% | 0.0097 / 0.0447 / 0.0493 |
| *Klebsiella quasipneumoniae* | 1463164.366 | 109 | 0 | 0 | 57.5% | 0.0086 / 0.0344 / 0.0345 |
| *Mycobacterium tuberculosis* | 1267357.3 | 4000 | 2 | 6 | 65.6% | 0.0004 / 0.0012 / 0.0033 |
| *Neisseria gonorrhoeae* | 485.17768 | 2000 | 0 | 0 | 52.5% | 0.0027 / 0.0045 / 0.0052 |
| *Pseudomonas aeruginosa* | 1402501.3 | 1421 | 0 | 11 | 66.1% | 0.0082 / 0.0158 / 0.0212 |
| *Salmonella enterica* | 28901.24485 | 2000 | 1 | 0 | 52.1% | 0.0110 / 0.0184 / 0.0467 |
| *Shigella flexneri* | 623.1307 | 362 | 0 | 20 | 50.5% | 0.0012 / 0.0068 / 0.0072 |
| *Shigella sonnei* | 624.3489 | 1984 | 0 | 31 | 50.8% | 0.0013 / 0.0040 / 0.0202 |
| *Staphylococcus aureus* | 1280.15939 | 2000 | 0 | 0 | 32.7% | 0.0087 / 0.0199 / 0.0232 |
| *Staphylococcus haemolyticus* | 1283.472 | 126 | 0 | 0 | 32.6% | 0.0055 / 0.0091 / 0.0099 |
| *Streptococcus agalactiae* | 1311.2996 | 141 | 0 | 0 | 35.3% | 0.0007 / 0.0131 / 0.0134 |
| *Streptococcus pneumoniae* | 1313.46229 | 2000 | 0 | 0 | 39.6% | 0.0097 / 0.0118 / 0.0196 |
| *Vibrio cholerae* | 666.4954 | 91 | 0 | 0 | 47.5% | 0.0009 / 0.0367 / 0.0367 |
| **Total** | | 25761 | 79 | 751 | | |

### Nearer-species pairs

| Labelled | Nearer reference | Genomes | Share of labelled | Outcome |
|---|---|--:|--:|---|
| *Escherichia coli* | *Shigella sonnei* | 532 | 26.60% | borderline (E. coli/Shigella group) |
| *Klebsiella pneumoniae* | *Klebsiella quasipneumoniae* | 67 | 3.35% | fail |
| *Shigella flexneri* | *Shigella sonnei* | 20 | 5.52% | borderline (E. coli/Shigella group) |
| *Shigella sonnei* | *Escherichia coli* | 30 | 1.51% | borderline (E. coli/Shigella group) |

### Dropped genomes

| Genome | Labelled | Own distance (ANI est.) | GC (species median) | Nearest other reference | Rule |
|---|---|--:|--:|---|---|
| 1310623.3 | *Acinetobacter baumannii* | 0.1823 (81.8%) | 42.8% (39.0%) | *Vibrio cholerae* 0.2960 | own distance > 0.13 |
| 1310638.3 | *Acinetobacter baumannii* | 0.1868 (81.3%) | 43.0% (39.0%) | *Klebsiella aerogenes* 0.2960 | own distance > 0.13 |
| 470.7458 | *Acinetobacter baumannii* | 0.1823 (81.8%) | 42.5% (39.0%) | *Vibrio cholerae* 0.2630 | own distance > 0.13 |
| 1400153.4 | *Enterobacter cloacae* | 0.1542 (84.6%) | 57.3% (55.1%) | *Klebsiella pneumoniae* 0.0103 | own distance > 0.13 |
| 562.12961 | *Escherichia coli* | 0.1823 (81.8%) | 57.5% (50.6%) | *Klebsiella pneumoniae* 0.0083 | own distance > 0.13 |
| 562.12962 | *Escherichia coli* | 0.1477 (85.2%) | 57.4% (50.6%) | *Klebsiella pneumoniae* 0.0016 | own distance > 0.13 |
| 562.42818 | *Escherichia coli* | 0.1457 (85.4%) | 51.9% (50.6%) | *Shigella sonnei* 0.1420 | own distance > 0.13 |
| 562.42857 | *Escherichia coli* | 0.1591 (84.1%) | 54.9% (50.6%) | *Klebsiella aerogenes* 0.0022 | own distance > 0.13 |
| 573.13133 | *Klebsiella pneumoniae* | 0.0427 (95.7%) | 57.6% (57.2%) | *Klebsiella quasipneumoniae* 0.0124 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.13146 | *Klebsiella pneumoniae* | 0.0432 (95.7%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0127 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.13344 | *Klebsiella pneumoniae* | 0.0504 (95.0%) | 57.6% (57.2%) | *Klebsiella quasipneumoniae* 0.0102 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.13556 | *Klebsiella pneumoniae* | 0.0535 (94.6%) | 57.7% (57.2%) | *Klebsiella quasipneumoniae* 0.0344 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.13614 | *Klebsiella pneumoniae* | 0.0539 (94.6%) | 57.7% (57.2%) | *Klebsiella quasipneumoniae* 0.0335 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.13924 | *Klebsiella pneumoniae* | 0.0535 (94.6%) | 57.8% (57.2%) | *Klebsiella quasipneumoniae* 0.0334 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.14338 | *Klebsiella pneumoniae* | 0.0527 (94.7%) | 56.7% (57.2%) | *Klebsiella quasipneumoniae* 0.0120 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.15563 | *Klebsiella pneumoniae* | 0.0531 (94.7%) | 58.3% (57.2%) | *Klebsiella quasipneumoniae* 0.0324 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.15572 | *Klebsiella pneumoniae* | 0.0510 (94.9%) | 57.8% (57.2%) | *Klebsiella quasipneumoniae* 0.0104 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.17927 | *Klebsiella pneumoniae* | 0.0519 (94.8%) | 57.3% (57.2%) | *Klebsiella quasipneumoniae* 0.0104 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.24249 | *Klebsiella pneumoniae* | 0.0506 (94.9%) | 57.6% (57.2%) | *Klebsiella quasipneumoniae* 0.0113 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.24250 | *Klebsiella pneumoniae* | 0.0515 (94.8%) | 57.9% (57.2%) | *Klebsiella quasipneumoniae* 0.0098 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.24342 | *Klebsiella pneumoniae* | 0.0513 (94.9%) | 57.9% (57.2%) | *Klebsiella quasipneumoniae* 0.0096 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.24384 | *Klebsiella pneumoniae* | 0.0510 (94.9%) | 57.4% (57.2%) | *Klebsiella quasipneumoniae* 0.0103 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26078 | *Klebsiella pneumoniae* | 0.0496 (95.0%) | 57.8% (57.2%) | *Klebsiella quasipneumoniae* 0.0087 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26087 | *Klebsiella pneumoniae* | 0.0508 (94.9%) | 57.3% (57.2%) | *Klebsiella quasipneumoniae* 0.0105 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26093 | *Klebsiella pneumoniae* | 0.0473 (95.3%) | 57.8% (57.2%) | *Klebsiella quasipneumoniae* 0.0097 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26126 | *Klebsiella pneumoniae* | 0.0517 (94.8%) | 58.1% (57.2%) | *Klebsiella quasipneumoniae* 0.0096 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26131 | *Klebsiella pneumoniae* | 0.0521 (94.8%) | 58.0% (57.2%) | *Klebsiella quasipneumoniae* 0.0102 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26133 | *Klebsiella pneumoniae* | 0.0493 (95.1%) | 57.8% (57.2%) | *Klebsiella quasipneumoniae* 0.0094 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26135 | *Klebsiella pneumoniae* | 0.0500 (95.0%) | 57.7% (57.2%) | *Klebsiella quasipneumoniae* 0.0090 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26145 | *Klebsiella pneumoniae* | 0.0510 (94.9%) | 57.8% (57.2%) | *Klebsiella quasipneumoniae* 0.0093 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26171 | *Klebsiella pneumoniae* | 0.0489 (95.1%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0096 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26177 | *Klebsiella pneumoniae* | 0.0529 (94.7%) | 57.9% (57.2%) | *Klebsiella quasipneumoniae* 0.0344 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26179 | *Klebsiella pneumoniae* | 0.0489 (95.1%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0096 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26182 | *Klebsiella pneumoniae* | 0.0489 (95.1%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0096 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26183 | *Klebsiella pneumoniae* | 0.0489 (95.1%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0096 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26184 | *Klebsiella pneumoniae* | 0.0489 (95.1%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0096 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26185 | *Klebsiella pneumoniae* | 0.0489 (95.1%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0096 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26186 | *Klebsiella pneumoniae* | 0.0489 (95.1%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0096 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26187 | *Klebsiella pneumoniae* | 0.0489 (95.1%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0096 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26188 | *Klebsiella pneumoniae* | 0.0489 (95.1%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0096 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26189 | *Klebsiella pneumoniae* | 0.0537 (94.6%) | 58.2% (57.2%) | *Klebsiella quasipneumoniae* 0.0330 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26190 | *Klebsiella pneumoniae* | 0.0489 (95.1%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0096 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26191 | *Klebsiella pneumoniae* | 0.0489 (95.1%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0096 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26192 | *Klebsiella pneumoniae* | 0.0489 (95.1%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0096 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26193 | *Klebsiella pneumoniae* | 0.0489 (95.1%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0096 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26209 | *Klebsiella pneumoniae* | 0.0544 (94.6%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0323 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26210 | *Klebsiella pneumoniae* | 0.0498 (95.0%) | 57.6% (57.2%) | *Klebsiella quasipneumoniae* 0.0095 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26219 | *Klebsiella pneumoniae* | 0.0504 (95.0%) | 57.8% (57.2%) | *Klebsiella quasipneumoniae* 0.0103 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26223 | *Klebsiella pneumoniae* | 0.0519 (94.8%) | 57.8% (57.2%) | *Klebsiella quasipneumoniae* 0.0104 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26225 | *Klebsiella pneumoniae* | 0.0539 (94.6%) | 57.7% (57.2%) | *Klebsiella quasipneumoniae* 0.0327 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26241 | *Klebsiella pneumoniae* | 0.0486 (95.1%) | 57.9% (57.2%) | *Klebsiella quasipneumoniae* 0.0311 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26242 | *Klebsiella pneumoniae* | 0.0498 (95.0%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0089 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26243 | *Klebsiella pneumoniae* | 0.0498 (95.0%) | 57.6% (57.2%) | *Klebsiella quasipneumoniae* 0.0089 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26261 | *Klebsiella pneumoniae* | 0.0498 (95.0%) | 57.7% (57.2%) | *Klebsiella quasipneumoniae* 0.0096 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26272 | *Klebsiella pneumoniae* | 0.0479 (95.2%) | 58.1% (57.2%) | *Klebsiella quasipneumoniae* 0.0087 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26284 | *Klebsiella pneumoniae* | 0.0508 (94.9%) | 57.8% (57.2%) | *Klebsiella quasipneumoniae* 0.0089 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.26318 | *Klebsiella pneumoniae* | 0.0513 (94.9%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0101 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.31140 | *Klebsiella pneumoniae* | 0.0489 (95.1%) | 57.7% (57.2%) | *Klebsiella quasipneumoniae* 0.0103 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.31144 | *Klebsiella pneumoniae* | 0.0508 (94.9%) | 57.6% (57.2%) | *Klebsiella quasipneumoniae* 0.0106 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.31147 | *Klebsiella pneumoniae* | 0.0517 (94.8%) | 57.2% (57.2%) | *Klebsiella quasipneumoniae* 0.0024 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.46098 | *Klebsiella pneumoniae* | 0.0510 (94.9%) | 57.8% (57.2%) | *Klebsiella quasipneumoniae* 0.0080 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.46109 | *Klebsiella pneumoniae* | 0.0517 (94.8%) | 57.6% (57.2%) | *Klebsiella quasipneumoniae* 0.0102 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.46111 | *Klebsiella pneumoniae* | 0.1542 (84.6%) | 55.4% (57.2%) | *Enterobacter cloacae* 0.0134 | own distance > 0.13 |
| 573.5698 | *Klebsiella pneumoniae* | 0.0527 (94.7%) | 57.2% (57.2%) | *Klebsiella quasipneumoniae* 0.0110 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.5709 | *Klebsiella pneumoniae* | 0.0513 (94.9%) | 57.4% (57.2%) | *Klebsiella quasipneumoniae* 0.0334 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.5716 | *Klebsiella pneumoniae* | 0.0511 (94.9%) | 57.9% (57.2%) | *Klebsiella quasipneumoniae* 0.0099 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.5719 | *Klebsiella pneumoniae* | 0.0561 (94.4%) | 57.3% (57.2%) | *Klebsiella quasipneumoniae* 0.0366 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.5802 | *Klebsiella pneumoniae* | 0.0513 (94.9%) | 57.5% (57.2%) | *Klebsiella quasipneumoniae* 0.0119 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.5811 | *Klebsiella pneumoniae* | 0.0527 (94.7%) | 57.2% (57.2%) | *Klebsiella quasipneumoniae* 0.0110 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.5823 | *Klebsiella pneumoniae* | 0.0511 (94.9%) | 57.9% (57.2%) | *Klebsiella quasipneumoniae* 0.0099 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.5824 | *Klebsiella pneumoniae* | 0.0556 (94.4%) | 57.3% (57.2%) | *Klebsiella quasipneumoniae* 0.0367 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.66104 | *Klebsiella pneumoniae* | 0.0527 (94.7%) | 57.4% (57.2%) | *Klebsiella quasipneumoniae* 0.0341 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.66150 | *Klebsiella pneumoniae* | 0.0559 (94.4%) | 57.1% (57.2%) | *Klebsiella quasipneumoniae* 0.0349 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.67200 | *Klebsiella pneumoniae* | 0.0527 (94.7%) | 57.8% (57.2%) | *Klebsiella quasipneumoniae* 0.0339 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.67202 | *Klebsiella pneumoniae* | 0.0504 (95.0%) | 57.6% (57.2%) | *Klebsiella quasipneumoniae* 0.0336 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 573.67203 | *Klebsiella pneumoniae* | 0.0527 (94.7%) | 57.8% (57.2%) | *Klebsiella quasipneumoniae* 0.0339 | nearer *Klebsiella quasipneumoniae*, within 0.05 |
| 1773.40949 | *Mycobacterium tuberculosis* | 1.0000 (0.0%) | 68.4% (65.6%) | *Acinetobacter baumannii* 1.0000 | own distance > 0.13 |
| 1773.41701 | *Mycobacterium tuberculosis* | 0.2960 (70.4%) | 64.2% (65.6%) | *Acinetobacter baumannii* 1.0000 | own distance > 0.13 |
| 90370.5242 | *Salmonella enterica* | 0.1647 (83.5%) | 55.3% (52.1%) | *Enterobacter cloacae* 0.0962 | own distance > 0.13 |

### Borderline genomes outside the E. coli/Shigella group

| Genome | Labelled | Own distance | GC (z) | Nearest other reference | Why |
|---|---|--:|--:|---|---|
| 1310642.3 | *Acinetobacter baumannii* | 0.0950 | 38.6% (-4.0) | *Staphylococcus aureus* 0.2960 | past the 0.05 species boundary |
| 1310669.3 | *Acinetobacter baumannii* | 0.0703 | 38.6% (-3.5) | *Staphylococcus aureus* 0.2960 | past the 0.05 species boundary |
| 1310759.3 | *Acinetobacter baumannii* | 0.0932 | 39.1% (+1.3) | *Enterococcus faecium* 0.2960 | past the 0.05 species boundary |
| 1310816.3 | *Acinetobacter baumannii* | 0.0632 | 39.0% (+0.2) | *Staphylococcus aureus* 0.2960 | past the 0.05 species boundary |
| 470.11118 | *Acinetobacter baumannii* | 0.1037 | 38.6% (-3.4) | *Staphylococcus aureus* 0.2630 | past the 0.05 species boundary |
| 470.7417 | *Acinetobacter baumannii* | 0.0926 | 38.8% (-2.1) | *Staphylococcus aureus* 0.2960 | past the 0.05 species boundary |
| 470.7472 | *Acinetobacter baumannii* | 0.0684 | 38.7% (-2.7) | *Klebsiella quasipneumoniae* 0.2960 | past the 0.05 species boundary |
| 470.7930 | *Acinetobacter baumannii* | 0.0752 | 38.3% (-6.8) | *Staphylococcus aureus* 0.2960 | past the 0.05 species boundary |
| 1329845.3 | *Enterobacter cloacae* | 0.0995 | 55.5% (+1.2) | *Klebsiella pneumoniae* 0.1386 | past the 0.05 species boundary |
| 2077136.3 | *Enterobacter cloacae* | 0.1111 | 54.4% (-2.1) | *Klebsiella aerogenes* 0.1519 | past the 0.05 species boundary |
| 2077137.3 | *Enterobacter cloacae* | 0.1077 | 54.3% (-2.5) | *Klebsiella aerogenes* 0.1354 | past the 0.05 species boundary |
| 550.1668 | *Enterobacter cloacae* | 0.0988 | 55.4% (+1.1) | *Klebsiella pneumoniae* 0.1438 | past the 0.05 species boundary |
| 550.2116 | *Enterobacter cloacae* | 0.1030 | 55.1% (+0.2) | *Klebsiella aerogenes* 0.1457 | past the 0.05 species boundary |
| 550.2119 | *Enterobacter cloacae* | 0.1111 | 54.6% (-1.7) | *Klebsiella aerogenes* 0.1542 | past the 0.05 species boundary |
| 550.2139 | *Enterobacter cloacae* | 0.0975 | 55.6% (+1.5) | *Klebsiella pneumoniae* 0.1403 | past the 0.05 species boundary |
| 550.2148 | *Enterobacter cloacae* | 0.1060 | 55.0% (-0.4) | *Klebsiella aerogenes* 0.1420 | past the 0.05 species boundary |
| 550.2736 | *Enterobacter cloacae* | 0.0853 | 54.3% (-2.4) | *Klebsiella pneumoniae* 0.1369 | past the 0.05 species boundary |
| 550.2737 | *Enterobacter cloacae* | 0.0853 | 54.3% (-2.4) | *Klebsiella pneumoniae* 0.1369 | past the 0.05 species boundary |
| 550.3503 | *Enterobacter cloacae* | 0.1111 | 54.6% (-1.7) | *Klebsiella aerogenes* 0.1566 | past the 0.05 species boundary |
| 550.3506 | *Enterobacter cloacae* | 0.1111 | 54.6% (-1.5) | *Klebsiella aerogenes* 0.1566 | past the 0.05 species boundary |
| 550.3525 | *Enterobacter cloacae* | 0.1045 | 55.0% (-0.4) | *Klebsiella aerogenes* 0.1386 | past the 0.05 species boundary |
| 550.3530 | *Enterobacter cloacae* | 0.1030 | 55.0% (-0.3) | *Klebsiella aerogenes* 0.1457 | past the 0.05 species boundary |
| 550.3540 | *Enterobacter cloacae* | 0.0988 | 55.0% (-0.4) | *Klebsiella aerogenes* 0.1457 | past the 0.05 species boundary |
| 550.3551 | *Enterobacter cloacae* | 0.1053 | 55.1% (-0.1) | *Klebsiella pneumoniae* 0.1457 | past the 0.05 species boundary |
| 550.3552 | *Enterobacter cloacae* | 0.1016 | 54.7% (-1.2) | *Klebsiella pneumoniae* 0.1354 | past the 0.05 species boundary |
| 550.3573 | *Enterobacter cloacae* | 0.1023 | 54.7% (-1.2) | *Klebsiella pneumoniae* 0.1369 | past the 0.05 species boundary |
| 550.3581 | *Enterobacter cloacae* | 0.1068 | 54.7% (-1.4) | *Klebsiella pneumoniae* 0.1438 | past the 0.05 species boundary |
| 550.4370 | *Enterobacter cloacae* | 0.1045 | 55.0% (-0.2) | *Klebsiella pneumoniae* 0.1386 | past the 0.05 species boundary |
| 550.4573 | *Enterobacter cloacae* | 0.0858 | 54.4% (-2.1) | *Klebsiella pneumoniae* 0.1354 | past the 0.05 species boundary |
| 550.4581 | *Enterobacter cloacae* | 0.1053 | 54.5% (-2.0) | *Klebsiella aerogenes* 0.1457 | past the 0.05 species boundary |
| 550.4583 | *Enterobacter cloacae* | 0.1008 | 56.0% (+3.0) | *Klebsiella aerogenes* 0.1498 | past the 0.05 species boundary |
| 550.4710 | *Enterobacter cloacae* | 0.1045 | 55.0% (-0.2) | *Klebsiella pneumoniae* 0.1386 | past the 0.05 species boundary |
| 550.542 | *Enterobacter cloacae* | 0.0909 | 55.1% (+0.1) | *Klebsiella pneumoniae* 0.1438 | past the 0.05 species boundary |
| 550.544 | *Enterobacter cloacae* | 0.0893 | 54.4% (-2.3) | *Klebsiella pneumoniae* 0.1386 | past the 0.05 species boundary |
| 550.547 | *Enterobacter cloacae* | 0.1023 | 54.9% (-0.7) | *Klebsiella aerogenes* 0.1420 | past the 0.05 species boundary |
| 550.622 | *Enterobacter cloacae* | 0.0938 | 56.2% (+3.7) | *Klebsiella pneumoniae* 0.1477 | past the 0.05 species boundary |
| 1352.10031 | *Enterococcus faecium* | 0.0539 | 38.0% (+0.9) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.10078 | *Enterococcus faecium* | 0.0525 | 37.9% (+0.9) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.10265 | *Enterococcus faecium* | 0.0537 | 38.0% (+1.3) | *Neisseria gonorrhoeae* 0.2960 | past the 0.05 species boundary |
| 1352.10320 | *Enterococcus faecium* | 0.0517 | 38.0% (+1.2) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.10341 | *Enterococcus faecium* | 0.0508 | 38.4% (+2.5) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.10435 | *Enterococcus faecium* | 0.0502 | 38.3% (+2.2) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.10464 | *Enterococcus faecium* | 0.0544 | 37.9% (+0.7) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.10465 | *Enterococcus faecium* | 0.0502 | 38.0% (+1.0) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.10474 | *Enterococcus faecium* | 0.0504 | 38.1% (+1.4) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.10476 | *Enterococcus faecium* | 0.0523 | 37.9% (+0.7) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.10615 | *Enterococcus faecium* | 0.0541 | 38.2% (+1.8) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.10630 | *Enterococcus faecium* | 0.0533 | 38.0% (+1.1) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.10637 | *Enterococcus faecium* | 0.0519 | 38.1% (+1.3) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.11081 | *Enterococcus faecium* | 0.0548 | 38.4% (+2.5) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.11083 | *Enterococcus faecium* | 0.0535 | 38.4% (+2.5) | *Staphylococcus haemolyticus* 0.2960 | past the 0.05 species boundary |
| 1352.11093 | *Enterococcus faecium* | 0.0513 | 38.1% (+1.6) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.18073 | *Enterococcus faecium* | 0.0535 | 38.2% (+1.9) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.20821 | *Enterococcus faecium* | 0.0511 | 38.0% (+1.1) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.20940 | *Enterococcus faecium* | 0.0502 | 38.1% (+1.6) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.20943 | *Enterococcus faecium* | 0.0502 | 38.2% (+1.8) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.20944 | *Enterococcus faecium* | 0.0506 | 38.2% (+1.8) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.20948 | *Enterococcus faecium* | 0.0521 | 38.3% (+2.1) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.20949 | *Enterococcus faecium* | 0.0548 | 37.9% (+0.9) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.20957 | *Enterococcus faecium* | 0.0550 | 38.0% (+1.0) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.20958 | *Enterococcus faecium* | 0.0508 | 38.2% (+2.1) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.20963 | *Enterococcus faecium* | 0.0535 | 38.2% (+1.8) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.20965 | *Enterococcus faecium* | 0.0515 | 38.2% (+2.0) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.20967 | *Enterococcus faecium* | 0.0510 | 38.4% (+2.5) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.20975 | *Enterococcus faecium* | 0.0511 | 38.0% (+1.2) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.20984 | *Enterococcus faecium* | 0.0525 | 38.2% (+1.9) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.20988 | *Enterococcus faecium* | 0.0510 | 38.2% (+1.9) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.20989 | *Enterococcus faecium* | 0.0544 | 38.1% (+1.4) | *Streptococcus agalactiae* 0.2195 | past the 0.05 species boundary |
| 1352.20992 | *Enterococcus faecium* | 0.0506 | 38.1% (+1.5) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.20994 | *Enterococcus faecium* | 0.0529 | 38.1% (+1.3) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.21017 | *Enterococcus faecium* | 0.0552 | 38.1% (+1.4) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.21018 | *Enterococcus faecium* | 0.0550 | 38.2% (+1.9) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.21032 | *Enterococcus faecium* | 0.0533 | 38.2% (+1.9) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.21066 | *Enterococcus faecium* | 0.0541 | 38.0% (+1.1) | *Streptococcus agalactiae* 0.2195 | past the 0.05 species boundary |
| 1352.21069 | *Enterococcus faecium* | 0.0539 | 38.0% (+1.1) | *Streptococcus agalactiae* 0.2195 | past the 0.05 species boundary |
| 1352.21072 | *Enterococcus faecium* | 0.0552 | 38.0% (+1.0) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.21078 | *Enterococcus faecium* | 0.0519 | 37.9% (+0.8) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.21080 | *Enterococcus faecium* | 0.0517 | 38.0% (+1.2) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.21081 | *Enterococcus faecium* | 0.0502 | 38.2% (+1.9) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.21083 | *Enterococcus faecium* | 0.0523 | 38.1% (+1.4) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.21084 | *Enterococcus faecium* | 0.0500 | 38.2% (+2.0) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.21087 | *Enterococcus faecium* | 0.0504 | 38.2% (+2.0) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.21090 | *Enterococcus faecium* | 0.0533 | 38.1% (+1.6) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.21091 | *Enterococcus faecium* | 0.0535 | 38.0% (+1.2) | *Streptococcus agalactiae* 0.2195 | past the 0.05 species boundary |
| 1352.21113 | *Enterococcus faecium* | 0.0500 | 38.2% (+2.0) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.21123 | *Enterococcus faecium* | 0.0504 | 38.2% (+2.0) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.4398 | *Enterococcus faecium* | 0.0521 | 38.3% (+2.2) | *Acinetobacter baumannii* 1.0000 | past the 0.05 species boundary |
| 1352.9236 | *Enterococcus faecium* | 0.0521 | 38.2% (+1.8) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.9253 | *Enterococcus faecium* | 0.0537 | 38.3% (+2.1) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.9480 | *Enterococcus faecium* | 0.0517 | 38.1% (+1.7) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.9521 | *Enterococcus faecium* | 0.0513 | 38.1% (+1.3) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.9574 | *Enterococcus faecium* | 0.0515 | 38.0% (+1.2) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.9785 | *Enterococcus faecium* | 0.0515 | 38.2% (+1.7) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 1352.9899 | *Enterococcus faecium* | 0.0508 | 38.3% (+2.5) | *Shigella flexneri* 0.2960 | past the 0.05 species boundary |
| 562.144151 | *Escherichia coli* | 0.0803 | 50.2% (-3.0) | *Shigella sonnei* 0.0812 | past the 0.05 species boundary |
| 562.144184 | *Escherichia coli* | 0.0713 | 50.6% (-0.5) | *Shigella sonnei* 0.0741 | past the 0.05 species boundary |
| 562.144221 | *Escherichia coli* | 0.0791 | 50.5% (-1.1) | *Shigella sonnei* 0.0783 | past the 0.05 species boundary |
| 562.145036 | *Escherichia coli* | 0.0812 | 50.2% (-2.9) | *Shigella sonnei* 0.0803 | past the 0.05 species boundary |
| 571.1924 | *Klebsiella oxytoca* | 0.0710 | 55.3% (+0.4) | *Enterobacter cloacae* 0.1138 | past the 0.05 species boundary |
| 571.1925 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1926 | *Klebsiella oxytoca* | 0.0621 | 55.6% (+1.2) | *Klebsiella pneumoniae* 0.1245 | past the 0.05 species boundary |
| 571.1929 | *Klebsiella oxytoca* | 0.0666 | 55.1% (-0.0) | *Klebsiella aerogenes* 0.1310 | past the 0.05 species boundary |
| 571.1930 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1932 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1934 | *Klebsiella oxytoca* | 0.0716 | 55.8% (+1.7) | *Klebsiella pneumoniae* 0.1296 | past the 0.05 species boundary |
| 571.1936 | *Klebsiella oxytoca* | 0.0616 | 55.4% (+0.6) | *Enterobacter cloacae* 0.1188 | past the 0.05 species boundary |
| 571.1937 | *Klebsiella oxytoca* | 0.0687 | 55.1% (-0.2) | *Klebsiella aerogenes* 0.1324 | past the 0.05 species boundary |
| 571.1942 | *Klebsiella oxytoca* | 0.0596 | 55.9% (+2.0) | *Klebsiella aerogenes* 0.1324 | past the 0.05 species boundary |
| 571.1943 | *Klebsiella oxytoca* | 0.0584 | 55.8% (+1.8) | *Klebsiella aerogenes* 0.1310 | past the 0.05 species boundary |
| 571.1944 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1946 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1947 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1948 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.5) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1949 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1950 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1951 | *Klebsiella oxytoca* | 0.0723 | 55.4% (+0.6) | *Klebsiella pneumoniae* 0.1233 | past the 0.05 species boundary |
| 571.1953 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1954 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1955 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1956 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1957 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1959 | *Klebsiella oxytoca* | 0.0678 | 55.9% (+2.1) | *Klebsiella pneumoniae* 0.1270 | past the 0.05 species boundary |
| 571.1960 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1961 | *Klebsiella oxytoca* | 0.0723 | 55.4% (+0.5) | *Klebsiella pneumoniae* 0.1245 | past the 0.05 species boundary |
| 571.1962 | *Klebsiella oxytoca* | 0.0723 | 55.4% (+0.5) | *Klebsiella pneumoniae* 0.1233 | past the 0.05 species boundary |
| 571.1963 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1964 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1966 | *Klebsiella oxytoca* | 0.0723 | 55.4% (+0.5) | *Klebsiella pneumoniae* 0.1245 | past the 0.05 species boundary |
| 571.1968 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1969 | *Klebsiella oxytoca* | 0.0601 | 55.7% (+1.4) | *Klebsiella aerogenes* 0.1310 | past the 0.05 species boundary |
| 571.1971 | *Klebsiella oxytoca* | 0.0621 | 55.6% (+1.2) | *Klebsiella pneumoniae* 0.1245 | past the 0.05 species boundary |
| 571.1972 | *Klebsiella oxytoca* | 0.0603 | 55.7% (+1.4) | *Klebsiella aerogenes* 0.1324 | past the 0.05 species boundary |
| 571.1975 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1977 | *Klebsiella oxytoca* | 0.0700 | 55.7% (+1.5) | *Klebsiella pneumoniae* 0.1233 | past the 0.05 species boundary |
| 571.1978 | *Klebsiella oxytoca* | 0.0716 | 55.6% (+1.2) | *Klebsiella pneumoniae* 0.1233 | past the 0.05 species boundary |
| 571.1985 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1990 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1991 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1992 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1994 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.1995 | *Klebsiella oxytoca* | 0.0559 | 56.0% (+2.2) | *Klebsiella aerogenes* 0.1270 | past the 0.05 species boundary |
| 571.2000 | *Klebsiella oxytoca* | 0.0621 | 55.7% (+1.4) | *Klebsiella aerogenes* 0.1324 | past the 0.05 species boundary |
| 571.2001 | *Klebsiella oxytoca* | 0.0621 | 55.6% (+1.2) | *Klebsiella pneumoniae* 0.1245 | past the 0.05 species boundary |
| 571.2003 | *Klebsiella oxytoca* | 0.0723 | 55.4% (+0.5) | *Klebsiella pneumoniae* 0.1245 | past the 0.05 species boundary |
| 571.2005 | *Klebsiella oxytoca* | 0.0723 | 55.7% (+1.4) | *Klebsiella pneumoniae* 0.1257 | past the 0.05 species boundary |
| 571.2007 | *Klebsiella oxytoca* | 0.0574 | 55.9% (+1.9) | *Klebsiella aerogenes* 0.1270 | past the 0.05 species boundary |
| 571.217 | *Klebsiella oxytoca* | 0.0598 | 55.5% (+1.0) | *Klebsiella aerogenes* 0.1296 | past the 0.05 species boundary |
| 571.318 | *Klebsiella oxytoca* | 0.0598 | 55.8% (+1.8) | *Klebsiella aerogenes* 0.1310 | past the 0.05 species boundary |
| 571.445 | *Klebsiella oxytoca* | 0.0611 | 55.9% (+2.1) | *Klebsiella aerogenes* 0.1296 | past the 0.05 species boundary |
| 573.15469 | *Klebsiella pneumoniae* | 0.0119 | 55.7% (-7.8) | *Klebsiella quasipneumoniae* 0.0581 | GC -7.8 SD off |
| 573.31070 | *Klebsiella pneumoniae* | 0.0095 | 58.3% (+6.1) | *Klebsiella quasipneumoniae* 0.0561 | GC +6.1 SD off |
| 1773.18177 | *Mycobacterium tuberculosis* | 0.0023 | 63.8% (-36.2) | *Acinetobacter baumannii* 1.0000 | GC -36.2 SD off |
| 1773.18213 | *Mycobacterium tuberculosis* | 0.0029 | 63.0% (-51.1) | *Streptococcus pneumoniae* 0.2195 | GC -51.1 SD off |
| 1773.18404 | *Mycobacterium tuberculosis* | 0.0014 | 64.6% (-20.6) | *Staphylococcus haemolyticus* 0.2630 | GC -20.6 SD off |
| 1773.3676 | *Mycobacterium tuberculosis* | 0.0022 | 64.1% (-30.5) | *Streptococcus pneumoniae* 0.1420 | GC -30.5 SD off |
| 1773.40722 | *Mycobacterium tuberculosis* | 0.0024 | 62.7% (-58.2) | *Acinetobacter baumannii* 1.0000 | GC -58.2 SD off |
| 1773.5309 | *Mycobacterium tuberculosis* | 0.0033 | 63.6% (-39.4) | *Enterococcus faecium* 0.2960 | GC -39.4 SD off |
| 287.1032 | *Pseudomonas aeruginosa* | 0.0550 | 66.2% (+0.2) | *Enterobacter cloacae* 0.2960 | past the 0.05 species boundary |
| 287.12442 | *Pseudomonas aeruginosa* | 0.0537 | 66.5% (+1.5) | *Enterobacter cloacae* 0.2630 | past the 0.05 species boundary |
| 287.12724 | *Pseudomonas aeruginosa* | 0.0539 | 66.4% (+1.1) | *Acinetobacter baumannii* 0.2960 | past the 0.05 species boundary |
| 287.3105 | *Pseudomonas aeruginosa* | 0.0511 | 66.4% (+1.4) | *Enterobacter cloacae* 0.2960 | past the 0.05 species boundary |
| 287.5689 | *Pseudomonas aeruginosa* | 0.0546 | 65.6% (-2.3) | *Acinetobacter baumannii* 0.2960 | past the 0.05 species boundary |
| 287.5751 | *Pseudomonas aeruginosa* | 0.0546 | 65.6% (-2.3) | *Acinetobacter baumannii* 0.2960 | past the 0.05 species boundary |
| 287.7779 | *Pseudomonas aeruginosa* | 0.0517 | 65.6% (-2.1) | *Klebsiella pneumoniae* 0.2630 | past the 0.05 species boundary |
| 287.7794 | *Pseudomonas aeruginosa* | 0.0517 | 65.6% (-2.1) | *Klebsiella pneumoniae* 0.2630 | past the 0.05 species boundary |
| 287.7805 | *Pseudomonas aeruginosa* | 0.0531 | 66.6% (+2.0) | *Enterobacter cloacae* 0.2960 | past the 0.05 species boundary |
| 287.876 | *Pseudomonas aeruginosa* | 0.0541 | 66.3% (+0.9) | *Enterobacter cloacae* 0.2630 | past the 0.05 species boundary |
| 287.992 | *Pseudomonas aeruginosa* | 0.0504 | 66.0% (-0.4) | *Enterobacter cloacae* 0.2960 | past the 0.05 species boundary |
| 624.4427 | *Shigella sonnei* | 0.0054 | 52.1% (+11.3) | *Shigella flexneri* 0.0214 | GC +11.3 SD off |

### What the drops are

- *M. tuberculosis*: exactly the two known genomes, both in the frozen 2k prefix, at positions 353 (`1773.41701`) and 527 (`1773.40949`) of the 4k list. None of the 2,000 genomes added in 4k fails.
- *K. pneumoniae*: 67 *K. quasipneumoniae* filed as *K. pneumoniae*, own distance 0.043 to 0.056 and 0.002 to 0.037 from the *K. quasipneumoniae* reference 1463164.366. The rule takes the whole cluster, including the 25 inside 0.05 of their own reference. One more, 573.46111, is an *Enterobacter* (0.013 from the *E. cloacae* reference).
- *E. coli*: three *Klebsiella* (0.002 to 0.008 from the *K. pneumoniae* or *K. aerogenes* reference) and one unidentified genome at 0.146.
- *E. cloacae*: one *K. pneumoniae* (0.010 from its reference).
- *A. baumannii*: three genomes at 0.18 to 0.19 with GC 42.5 to 43.0% against 39.0%: another species, none of the 23.
- *S. enterica*: one genome filed under taxon 90370 (serovar Typhi) at 0.165, GC 55.3% against 52.1%: unidentified, nearest the *E. cloacae* reference at 0.096.

### MTB GC and length outliers

Six *M. tuberculosis* genomes are *M. tuberculosis* by sequence (own distance 0.0014 to 0.0033) but carry GC 62.7 to 64.6% against a 65.6% median, and assemblies of 4.52 to 4.83 Mb against about 4.41 Mb. They likely carry foreign contigs. By the admin ruling of 2026-10-29 they stay in and are recorded here. Their 4k-list positions are 115, 436, 869, 1086, 1265 and 3630; only the last is in the added 2k.

## Revision-10 notes

- Add a contig-level contamination check: genomes of the right species carrying foreign contigs, such as the six MTB outliers above, pass a whole-genome identity check.
- The nearer-species rule measures against one medoid per species. Species complexes measured from their largest member (*E. cloacae*, *K. oxytoca*) stay borderline. Splitting them would need per-complex references.
