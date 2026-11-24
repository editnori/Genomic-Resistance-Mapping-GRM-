# Lineage curve result

This page reads block M4 of the mechanism split once, by the rules `docs/mechanism-split.md` fixed before any run (section "M3 and M4"). Among the resistant genomes of each antibiotic of *Mycobacterium tuberculosis*, it counts the distinct graded variants of the WHO catalogue that appear, how that count grows as lineages are added, how many resistant genomes carry no graded variant, and how many variants are seen in one lineage only. It fits nothing: it reads TBProfiler's output beside the release's reports.

"The WHO catalogue" here is the catalogue as TBProfiler 6.6.3 ships it, `tbdb` commit `3ed1cb99`, whose graded entries carry the source label "WHO catalogue v2". A graded variant of an antibiotic is a variant the tool reports for the antibiotic's drug, passing the tool's own filter, with an entry graded "Assoc w R" or "Assoc w R - Interim". A report's antibiotic meets the tool's drug by the TBProfiler comparison's name map: rifampin is rifampicin, and ciprofloxacin, ofloxacin, rifabutin and prothionamide read another drug's calls by a shared name stem (proxy).

A curve describes the deposits of one public database. A share is no prevalence. A split by lineage is a retrospective split inside one database, never external validation. A table names no cause and no new marker. Nothing here is a change comparison: no default, no released file and no sentence on the site changes.

## Reading

Over the resistant genomes of the 20 reports the tool can read, the graded variants number from 0 to 136 per antibiotic, and how fast the count fills as main lineages are added differs by antibiotic. For 7 of the 18 antibiotics with a graded variant, the median over 1,000 random orders reaches the full count by the third main lineage added: amikacin (4 variants among 548 resistant genomes), ciprofloxacin (8 among 73), ethambutol (12 among 1,338), kanamycin (6 among 617), linezolid (3 among 202), ofloxacin (12 among 355) and prothionamide (14 among 92). For 5 it reaches the full count only with the last: clofazimine (5 among 193), ethionamide (136 among 915), pyrazinamide (94 among 381), rifabutin (27 among 665) and rifampin (40 among 1,365).

Resistant genomes with no graded variant run from 20 of 665 for rifabutin (0.030) to 444 of 448 for delamanid (0.991), and to all of them for cycloserine (119 of 119) and para-aminosalicylic acid (429 of 429), whose curves are 0 (section Limits). The share is under 0.1 for rifabutin, rifampin (72 of 1,365, 0.053), isoniazid (77 of 1,388, 0.055), streptomycin (44 of 539, 0.082) and ofloxacin (34 of 355, 0.096). It is over 0.8 for linezolid (178 of 202, 0.881), clofazimine (188 of 193, 0.974), delamanid, cycloserine and para-aminosalicylic acid.

Graded variants seen in one main lineage only: 121 of ethionamide's 136, 78 of pyrazinamide's 94, 20 of rifampin's 40, 16 of rifabutin's 27, 13 of prothionamide's 14 and all 5 of clofazimine's. None of amikacin's 4, ethambutol's 12, kanamycin's 6 or linezolid's 3 is. The resistant genomes that carry one of them: 245 of 915 for ethionamide (0.268), 196 of 381 for pyrazinamide (0.514) and 26 of 1,365 for rifampin (0.019).

Clofazimine and isoniazid are in set M and are marked. A curve describes the deposits of one public database, and a split by lineage is a retrospective split inside one database, never external validation.

Per report: the graded variants among its resistant genomes, the main-lineage curve's first step (median and band over the random orders) and last, the resistant genomes with no graded variant, and the graded variants seen in one main lineage only. Main lineages counts the groups as the tool writes them; a call that joins two lineages is a group of its own.

| Antibiotic | Kind | Set | Resistant genomes | Graded variants | Main lineages | First lineage, median (band) | All lineages | No graded variant (share) | In one main lineage only |
|---|---|---|---|---|---|---|---|---|---|
| amikacin | primary | not M | 548 | 4 | 6 | 3 (0 to 4) | 4 | 116 (0.212) | 0 |
| capreomycin | primary | not M | 422 | 10 | 5 | 2 (0 to 7) | 10 | 75 (0.178) | 8 |
| ciprofloxacin | proxy | not M | 73 | 8 | 4 | 4 (1 to 7) | 8 | 8 (0.110) | 3 |
| clofazimine | primary | M | 193 | 5 | 4 | 1 (1 to 2) | 5 | 188 (0.974) | 5 |
| cycloserine | primary | not M | 119 | 0 | 3 | 0 (0 to 0) | 0 | 119 (1.000) | 0 |
| delamanid | primary | not M | 448 | 3 | 6 | 0 (0 to 2) | 3 | 444 (0.991) | 2 |
| ethambutol | primary | not M | 1,338 | 12 | 8 | 5 (0 to 12) | 12 | 347 (0.259) | 0 |
| ethionamide | primary | not M | 915 | 136 | 7 | 8 (1 to 101) | 136 | 259 (0.283) | 121 |
| isoniazid | primary | M | 1,388 | 17 | 9 | 2 (0 to 11) | 17 | 77 (0.055) | 9 |
| kanamycin | primary | not M | 617 | 6 | 6 | 1 (0 to 6) | 6 | 128 (0.207) | 0 |
| levofloxacin | primary | not M | 753 | 16 | 8 | 4 (0 to 12) | 16 | 311 (0.413) | 6 |
| linezolid | primary | not M | 202 | 3 | 5 | 3 (0 to 3) | 3 | 178 (0.881) | 0 |
| moxifloxacin | primary | not M | 733 | 15 | 8 | 1 (0 to 13) | 15 | 136 (0.186) | 4 |
| ofloxacin | proxy | not M | 355 | 12 | 4 | 9 (0 to 11) | 12 | 34 (0.096) | 4 |
| para-aminosalicylic acid | primary | not M | 429 | 0 | 6 | 0 (0 to 0) | 0 | 429 (1.000) | 0 |
| prothionamide | proxy | not M | 92 | 14 | 2 | 14 (1 to 14) | 14 | 26 (0.283) | 13 |
| pyrazinamide | primary | not M | 381 | 94 | 5 | 15 (1 to 53) | 94 | 60 (0.157) | 78 |
| rifabutin | proxy | not M | 665 | 27 | 6 | 9 (1 to 17) | 27 | 20 (0.030) | 16 |
| rifampin | primary | not M | 1,365 | 40 | 7 | 10 (1 to 26) | 40 | 72 (0.053) | 20 |
| streptomycin | primary | not M | 539 | 15 | 5 | 6 (1 to 12) | 15 | 44 (0.082) | 9 |

## Data

**Genomes.** TBProfiler 6.6.3 ran again on the 2,000 genomes of the species' revision 14 reports, from their assemblies (FASTA), one genome per run with one thread, as the TBProfiler comparison ran it on 14 Nov 2026. That run's per-genome output was not kept, so the tool ran again. All 2,000 have a result and none failed (section Method).

**Reports.** The species has 22 reports. The tool reads 20 of them: 16 by the antibiotic's own drug (primary) and 4 by a proxy, ciprofloxacin and ofloxacin by levofloxacin and moxifloxacin, prothionamide by ethionamide and rifabutin by rifampicin. Amoxicillin and nicotinamide are not read: the tool has no drug of either name. A resistant genome is one the report's label calls resistant. Each report's genomes are its own; a genome can be resistant to several antibiotics and is counted once in each. Every resistant genome of the 20 reports has a result.

**Lineages.** Three groupings stand side by side: the tool's main lineage, its sublineage, and the report's own lineage cluster. A lineage is named as TBProfiler prints it. A genome the tool gives no lineage would enter no group of that grouping and be counted in the Datasets section; every resistant genome of the 20 reports has a main lineage, a sublineage and a cluster. A mixed call is its own group, as the tool writes it. Over the 20 reports the resistant genomes fall in 10 main-lineage groups: lineage1, lineage2, lineage3, lineage4, lineage6, La3, and four calls that join two lineages, La1;lineage4, lineage2;lineage3, lineage2;lineage4 and lineage3;lineage4. Per report they fall in 2 to 9 main-lineage groups, 5 to 46 sublineages and 34 to 584 clusters.

**Set M.** The reports of clofazimine and isoniazid are in set M: their labels wait for the maintainer's label rule. They are marked in every table and read with the release's labels, as the method page fixed.

## Curves

The curve counts the distinct graded variants among the resistant genomes as lineage groups are added one at a time. Over 1,000 random orders of the groups, seeded by the antibiotic and the grouping, each step gives the median and the 2.5 and 97.5 percentiles (the band). Beside it, the order by group size, largest first. No curve is drawn over under 5 resistant genomes; every read report has 5 or more.

Main lineages, every step:

| Antibiotic | Groups added | Median (2.5 to 97.5 percentile) | Largest first |
|---|---|---|---|
| amikacin | 1 of 6 | 3 (0 to 4) | 4 |
| amikacin | 2 of 6 | 4 (1 to 4) | 4 |
| amikacin | 3 of 6 | 4 (1 to 4) | 4 |
| amikacin | 4 of 6 | 4 (3 to 4) | 4 |
| amikacin | 5 of 6 | 4 (4 to 4) | 4 |
| amikacin | 6 of 6 | 4 (4 to 4) | 4 |
| capreomycin | 1 of 5 | 2 (0 to 7) | 7 |
| capreomycin | 2 of 5 | 5 (1 to 10) | 10 |
| capreomycin | 3 of 5 | 7 (2 to 10) | 10 |
| capreomycin | 4 of 5 | 10 (5 to 10) | 10 |
| capreomycin | 5 of 5 | 10 (10 to 10) | 10 |
| ciprofloxacin | 1 of 4 | 4 (1 to 7) | 7 |
| ciprofloxacin | 2 of 4 | 7 (4 to 8) | 8 |
| ciprofloxacin | 3 of 4 | 8 (6 to 8) | 8 |
| ciprofloxacin | 4 of 4 | 8 (8 to 8) | 8 |
| clofazimine | 1 of 4 | 1 (1 to 2) | 2 |
| clofazimine | 2 of 4 | 3 (2 to 3) | 3 |
| clofazimine | 3 of 4 | 4 (3 to 4) | 4 |
| clofazimine | 4 of 4 | 5 (5 to 5) | 5 |
| cycloserine | 1 of 3 | 0 (0 to 0) | 0 |
| cycloserine | 2 of 3 | 0 (0 to 0) | 0 |
| cycloserine | 3 of 3 | 0 (0 to 0) | 0 |
| delamanid | 1 of 6 | 0 (0 to 2) | 1 |
| delamanid | 2 of 6 | 1 (0 to 3) | 3 |
| delamanid | 3 of 6 | 2 (0 to 3) | 3 |
| delamanid | 4 of 6 | 2 (1 to 3) | 3 |
| delamanid | 5 of 6 | 3 (2 to 3) | 3 |
| delamanid | 6 of 6 | 3 (3 to 3) | 3 |
| ethambutol | 1 of 8 | 5 (0 to 12) | 12 |
| ethambutol | 2 of 8 | 11 (1 to 12) | 12 |
| ethambutol | 3 of 8 | 12 (2 to 12) | 12 |
| ethambutol | 4 of 8 | 12 (5 to 12) | 12 |
| ethambutol | 5 of 8 | 12 (11 to 12) | 12 |
| ethambutol | 6 of 8 | 12 (11 to 12) | 12 |
| ethambutol | 7 of 8 | 12 (12 to 12) | 12 |
| ethambutol | 8 of 8 | 12 (12 to 12) | 12 |
| ethionamide | 1 of 7 | 8 (1 to 101) | 101 |
| ethionamide | 2 of 7 | 35 (2 to 127) | 127 |
| ethionamide | 3 of 7 | 46 (3 to 129) | 129 |
| ethionamide | 4 of 7 | 104 (10 to 135) | 135 |
| ethionamide | 5 of 7 | 110 (19 to 136) | 135 |
| ethionamide | 6 of 7 | 134 (49 to 136) | 136 |
| ethionamide | 7 of 7 | 136 (136 to 136) | 136 |
| isoniazid | 1 of 9 | 2 (0 to 11) | 10 |
| isoniazid | 2 of 9 | 8 (1 to 14) | 12 |
| isoniazid | 3 of 9 | 10 (2 to 15) | 16 |
| isoniazid | 4 of 9 | 12 (4 to 16) | 16 |
| isoniazid | 5 of 9 | 13 (8 to 17) | 17 |
| isoniazid | 6 of 9 | 14 (9 to 17) | 17 |
| isoniazid | 7 of 9 | 15 (10 to 17) | 17 |
| isoniazid | 8 of 9 | 17 (13 to 17) | 17 |
| isoniazid | 9 of 9 | 17 (17 to 17) | 17 |
| kanamycin | 1 of 6 | 1 (0 to 6) | 6 |
| kanamycin | 2 of 6 | 6 (1 to 6) | 6 |
| kanamycin | 3 of 6 | 6 (1 to 6) | 6 |
| kanamycin | 4 of 6 | 6 (4 to 6) | 6 |
| kanamycin | 5 of 6 | 6 (6 to 6) | 6 |
| kanamycin | 6 of 6 | 6 (6 to 6) | 6 |
| levofloxacin | 1 of 8 | 4 (0 to 12) | 12 |
| levofloxacin | 2 of 8 | 6 (1 to 15) | 13 |
| levofloxacin | 3 of 8 | 12 (1 to 15) | 16 |
| levofloxacin | 4 of 8 | 13 (4 to 16) | 16 |
| levofloxacin | 5 of 8 | 13 (6 to 16) | 16 |
| levofloxacin | 6 of 8 | 15 (7 to 16) | 16 |
| levofloxacin | 7 of 8 | 16 (13 to 16) | 16 |
| levofloxacin | 8 of 8 | 16 (16 to 16) | 16 |
| linezolid | 1 of 5 | 3 (0 to 3) | 3 |
| linezolid | 2 of 5 | 3 (0 to 3) | 3 |
| linezolid | 3 of 5 | 3 (3 to 3) | 3 |
| linezolid | 4 of 5 | 3 (3 to 3) | 3 |
| linezolid | 5 of 5 | 3 (3 to 3) | 3 |
| moxifloxacin | 1 of 8 | 1 (0 to 13) | 13 |
| moxifloxacin | 2 of 8 | 7 (1 to 15) | 15 |
| moxifloxacin | 3 of 8 | 12 (1 to 15) | 15 |
| moxifloxacin | 4 of 8 | 13 (4 to 15) | 15 |
| moxifloxacin | 5 of 8 | 13 (6 to 15) | 15 |
| moxifloxacin | 6 of 8 | 15 (7 to 15) | 15 |
| moxifloxacin | 7 of 8 | 15 (13 to 15) | 15 |
| moxifloxacin | 8 of 8 | 15 (15 to 15) | 15 |
| ofloxacin | 1 of 4 | 9 (0 to 11) | 9 |
| ofloxacin | 2 of 4 | 9 (1 to 12) | 12 |
| ofloxacin | 3 of 4 | 12 (9 to 12) | 12 |
| ofloxacin | 4 of 4 | 12 (12 to 12) | 12 |
| para-aminosalicylic acid | 1 of 6 | 0 (0 to 0) | 0 |
| para-aminosalicylic acid | 2 of 6 | 0 (0 to 0) | 0 |
| para-aminosalicylic acid | 3 of 6 | 0 (0 to 0) | 0 |
| para-aminosalicylic acid | 4 of 6 | 0 (0 to 0) | 0 |
| para-aminosalicylic acid | 5 of 6 | 0 (0 to 0) | 0 |
| para-aminosalicylic acid | 6 of 6 | 0 (0 to 0) | 0 |
| prothionamide | 1 of 2 | 14 (1 to 14) | 14 |
| prothionamide | 2 of 2 | 14 (14 to 14) | 14 |
| pyrazinamide | 1 of 5 | 15 (1 to 53) | 53 |
| pyrazinamide | 2 of 5 | 47.5 (4 to 80) | 80 |
| pyrazinamide | 3 of 5 | 57 (19 to 90) | 90 |
| pyrazinamide | 4 of 5 | 84 (57 to 93) | 93 |
| pyrazinamide | 5 of 5 | 94 (94 to 94) | 94 |
| rifabutin | 1 of 6 | 9 (1 to 17) | 15 |
| rifabutin | 2 of 6 | 15 (3 to 23) | 23 |
| rifabutin | 3 of 6 | 19 (9 to 25) | 25 |
| rifabutin | 4 of 6 | 21 (15 to 27) | 27 |
| rifabutin | 5 of 6 | 25 (21 to 27) | 27 |
| rifabutin | 6 of 6 | 27 (27 to 27) | 27 |
| rifampin | 1 of 7 | 10 (1 to 26) | 26 |
| rifampin | 2 of 7 | 21 (1 to 35) | 35 |
| rifampin | 3 of 7 | 26 (3 to 38) | 39 |
| rifampin | 4 of 7 | 30 (10 to 40) | 40 |
| rifampin | 5 of 7 | 34 (22 to 40) | 40 |
| rifampin | 6 of 7 | 39 (29 to 40) | 40 |
| rifampin | 7 of 7 | 40 (40 to 40) | 40 |
| streptomycin | 1 of 5 | 6 (1 to 12) | 6 |
| streptomycin | 2 of 5 | 8 (2 to 15) | 13 |
| streptomycin | 3 of 5 | 12 (6 to 15) | 15 |
| streptomycin | 4 of 5 | 15 (8 to 15) | 15 |
| streptomycin | 5 of 5 | 15 (15 to 15) | 15 |

Sublineages and clusters, the first step and the last. Under the clusters, the median at the first group added is 0 or 1 for every antibiotic.

| Antibiotic | Sublineages | First sublineage, median (band) | Largest sublineage | All sublineages | Clusters | First cluster, median (band) | Largest cluster | All clusters |
|---|---|---|---|---|---|---|---|---|
| amikacin | 31 | 1 (0 to 4) | 4 | 4 | 249 | 1 (0 to 2) | 2 | 4 |
| capreomycin | 19 | 1 (0 to 7) | 7 | 10 | 152 | 1 (0 to 2) | 2 | 10 |
| ciprofloxacin | 17 | 1 (0 to 4) | 4 | 8 | 46 | 1 (0 to 2) | 1 | 8 |
| clofazimine | 19 | 0 (0 to 2) | 2 | 5 | 129 | 0 (0 to 1) | 1 | 5 |
| cycloserine | 10 | 0 (0 to 0) | 0 | 0 | 53 | 0 (0 to 0) | 0 | 0 |
| delamanid | 28 | 0 (0 to 2) | 1 | 3 | 276 | 0 (0 to 0) | 0 | 3 |
| ethambutol | 44 | 1 (0 to 10) | 11 | 12 | 584 | 1 (0 to 3) | 8 | 12 |
| ethionamide | 39 | 1 (0 to 83) | 83 | 136 | 380 | 1 (0 to 4) | 31 | 136 |
| isoniazid | 46 | 2 (0 to 9) | 9 | 17 | 564 | 1 (0 to 3) | 4 | 17 |
| kanamycin | 33 | 1 (0 to 6) | 6 | 6 | 258 | 1 (0 to 3) | 5 | 6 |
| levofloxacin | 33 | 1 (0 to 12) | 12 | 16 | 381 | 0 (0 to 3) | 8 | 16 |
| linezolid | 21 | 0 (0 to 3) | 2 | 3 | 148 | 0 (0 to 1) | 3 | 3 |
| moxifloxacin | 36 | 1 (0 to 9) | 12 | 15 | 347 | 1 (0 to 4) | 8 | 15 |
| ofloxacin | 17 | 1 (0 to 9) | 8 | 12 | 118 | 1 (0 to 5) | 8 | 12 |
| para-aminosalicylic acid | 27 | 0 (0 to 0) | 0 | 0 | 209 | 0 (0 to 0) | 0 | 0 |
| prothionamide | 5 | 0 (0 to 14) | 14 | 14 | 34 | 1 (0 to 5) | 5 | 14 |
| pyrazinamide | 24 | 2 (0 to 51) | 51 | 94 | 184 | 1 (0 to 5) | 4 | 94 |
| rifabutin | 26 | 2 (0 to 12) | 12 | 27 | 277 | 1 (0 to 2) | 2 | 27 |
| rifampin | 41 | 2 (0 to 16) | 23 | 40 | 541 | 1 (0 to 3) | 4 | 40 |
| streptomycin | 27 | 1 (0 to 6) | 5 | 15 | 196 | 1 (0 to 2) | 1 | 15 |

## No graded variant

A resistant genome with no graded variant carries no variant of the antibiotic's drug that the catalogue grades. It is a count, and nobody looked into these genomes further. A share is written over 5 or more resistant genomes. Over all the resistant genomes of each report, then per main lineage; a group of under 5 resistant genomes gets its count and no share, and `none` marks a group with no resistant genome of the antibiotic.

| Antibiotic | No graded variant | Share | La1;lineage4 | La3 | lineage1 | lineage2 | lineage2;lineage3 | lineage2;lineage4 | lineage3 | lineage3;lineage4 | lineage4 | lineage6 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| amikacin | 116 of 548 | 0.212 | none | none | 11 of 15 (0.733) | 69 of 389 (0.177) | 1 of 1 | 0 of 1 | 15 of 48 (0.312) | none | 20 of 94 (0.213) | none |
| capreomycin | 75 of 422 | 0.178 | none | none | 1 of 1 | 56 of 272 (0.206) | none | 0 of 1 | 1 of 9 (0.111) | none | 17 of 139 (0.122) | none |
| ciprofloxacin | 8 of 73 | 0.110 | none | none | 3 of 5 (0.600) | 0 of 7 (0.000) | none | none | 3 of 30 (0.100) | none | 2 of 31 (0.065) | none |
| clofazimine | 188 of 193 | 0.974 | none | none | 31 of 32 (0.969) | 80 of 82 (0.976) | none | none | 59 of 60 (0.983) | none | 18 of 19 (0.947) | none |
| cycloserine | 119 of 119 | 1.000 | none | none | none | 85 of 85 (1.000) | none | none | 1 of 1 | none | 33 of 33 (1.000) | none |
| delamanid | 444 of 448 | 0.991 | none | none | 73 of 73 (1.000) | 187 of 188 (0.995) | none | 2 of 2 | 122 of 124 (0.984) | 1 of 1 | 59 of 60 (0.983) | none |
| ethambutol | 347 of 1,338 | 0.259 | none | none | 59 of 82 (0.720) | 101 of 760 (0.133) | 0 of 1 | 3 of 4 | 118 of 209 (0.565) | 1 of 1 | 65 of 280 (0.232) | 0 of 1 |
| ethionamide | 259 of 915 | 0.283 | none | none | 79 of 106 (0.745) | 71 of 498 (0.143) | 0 of 1 | 0 of 1 | 54 of 83 (0.651) | 0 of 1 | 55 of 225 (0.244) | none |
| isoniazid | 77 of 1,388 | 0.055 | none | 1 of 1 | 16 of 63 (0.254) | 14 of 845 (0.017) | 0 of 1 | 1 of 2 | 26 of 169 (0.154) | 0 of 1 | 19 of 303 (0.063) | 0 of 3 |
| kanamycin | 128 of 617 | 0.207 | none | none | 14 of 18 (0.778) | 62 of 363 (0.171) | 1 of 1 | 0 of 1 | 25 of 59 (0.424) | none | 26 of 175 (0.149) | none |
| levofloxacin | 311 of 753 | 0.413 | 0 of 1 | none | 91 of 97 (0.938) | 39 of 372 (0.105) | 0 of 1 | 1 of 3 | 125 of 173 (0.723) | 1 of 1 | 54 of 105 (0.514) | none |
| linezolid | 178 of 202 | 0.881 | 1 of 1 | none | 37 of 37 (1.000) | 38 of 53 (0.717) | none | none | 69 of 74 (0.932) | none | 33 of 37 (0.892) | none |
| moxifloxacin | 136 of 733 | 0.186 | 0 of 1 | none | 52 of 58 (0.897) | 27 of 485 (0.056) | 0 of 1 | 1 of 3 | 19 of 67 (0.284) | 1 of 1 | 36 of 117 (0.308) | none |
| ofloxacin | 34 of 355 | 0.096 | none | none | none | 21 of 224 (0.094) | none | 1 of 1 | 2 of 3 | none | 10 of 127 (0.079) | none |
| para-aminosalicylic acid | 429 of 429 | 1.000 | none | none | 19 of 19 (1.000) | 285 of 285 (1.000) | none | 1 of 1 | 53 of 53 (1.000) | 1 of 1 | 70 of 70 (1.000) | none |
| prothionamide | 26 of 92 | 0.283 | none | none | none | 21 of 79 (0.266) | none | none | none | none | 5 of 13 (0.385) | none |
| pyrazinamide | 60 of 381 | 0.157 | none | none | 2 of 6 (0.333) | 36 of 226 (0.159) | none | none | 8 of 30 (0.267) | none | 14 of 118 (0.119) | 0 of 1 |
| rifabutin | 20 of 665 | 0.030 | none | none | 7 of 26 (0.269) | 0 of 480 (0.000) | 0 of 1 | 0 of 1 | 9 of 85 (0.106) | none | 4 of 72 (0.056) | none |
| rifampin | 72 of 1,365 | 0.053 | none | none | 18 of 50 (0.360) | 3 of 858 (0.003) | 0 of 1 | 0 of 2 | 35 of 160 (0.219) | none | 16 of 293 (0.055) | 0 of 1 |
| streptomycin | 44 of 539 | 0.082 | none | none | 0 of 3 | 12 of 336 (0.036) | none | 0 of 1 | 4 of 22 (0.182) | none | 28 of 177 (0.158) | none |

## Spread over lineages

Each graded variant is counted by the number of lineages it is seen in: 1, 2, 3, or 4 or more. Beside each count, the resistant genomes that carry at least one variant of that count, and their share of the report's resistant genomes.

| Antibiotic | Grouping | In 1 (resistant covered, share) | In 2 | In 3 | In 4 or more |
|---|---|---|---|---|---|
| amikacin | main lineage | 0 (0, 0.000) | 1 (2, 0.004) | 2 (23, 0.042) | 1 (407, 0.743) |
| amikacin | sublineage | 0 (0, 0.000) | 1 (2, 0.004) | 0 (0, 0.000) | 3 (430, 0.785) |
| capreomycin | main lineage | 8 (23, 0.055) | 0 (0, 0.000) | 1 (7, 0.017) | 1 (317, 0.751) |
| capreomycin | sublineage | 8 (23, 0.055) | 0 (0, 0.000) | 0 (0, 0.000) | 2 (324, 0.768) |
| ciprofloxacin | main lineage | 3 (5, 0.068) | 3 (11, 0.151) | 1 (23, 0.315) | 1 (26, 0.356) |
| ciprofloxacin | sublineage | 3 (5, 0.068) | 1 (2, 0.027) | 1 (3, 0.041) | 3 (55, 0.753) |
| clofazimine | main lineage | 5 (5, 0.026) | 0 (0, 0.000) | 0 (0, 0.000) | 0 (0, 0.000) |
| clofazimine | sublineage | 5 (5, 0.026) | 0 (0, 0.000) | 0 (0, 0.000) | 0 (0, 0.000) |
| cycloserine | main lineage | 0 (0, 0.000) | 0 (0, 0.000) | 0 (0, 0.000) | 0 (0, 0.000) |
| cycloserine | sublineage | 0 (0, 0.000) | 0 (0, 0.000) | 0 (0, 0.000) | 0 (0, 0.000) |
| delamanid | main lineage | 2 (2, 0.004) | 1 (2, 0.004) | 0 (0, 0.000) | 0 (0, 0.000) |
| delamanid | sublineage | 2 (2, 0.004) | 1 (2, 0.004) | 0 (0, 0.000) | 0 (0, 0.000) |
| ethambutol | main lineage | 0 (0, 0.000) | 1 (8, 0.006) | 6 (69, 0.052) | 5 (919, 0.687) |
| ethambutol | sublineage | 0 (0, 0.000) | 0 (0, 0.000) | 2 (23, 0.017) | 10 (971, 0.726) |
| ethionamide | main lineage | 121 (245, 0.268) | 9 (27, 0.030) | 2 (69, 0.075) | 4 (370, 0.404) |
| ethionamide | sublineage | 117 (228, 0.249) | 13 (44, 0.048) | 0 (0, 0.000) | 6 (438, 0.479) |
| isoniazid | main lineage | 9 (10, 0.007) | 1 (9, 0.006) | 2 (92, 0.066) | 5 (1,292, 0.931) |
| isoniazid | sublineage | 9 (10, 0.007) | 0 (0, 0.000) | 0 (0, 0.000) | 8 (1,304, 0.939) |
| kanamycin | main lineage | 0 (0, 0.000) | 2 (50, 0.081) | 3 (64, 0.104) | 1 (381, 0.618) |
| kanamycin | sublineage | 0 (0, 0.000) | 1 (44, 0.071) | 1 (6, 0.010) | 4 (445, 0.721) |
| levofloxacin | main lineage | 6 (12, 0.016) | 5 (39, 0.052) | 2 (72, 0.096) | 3 (341, 0.453) |
| levofloxacin | sublineage | 5 (7, 0.009) | 2 (11, 0.015) | 3 (27, 0.036) | 6 (411, 0.546) |
| linezolid | main lineage | 0 (0, 0.000) | 0 (0, 0.000) | 3 (24, 0.119) | 0 (0, 0.000) |
| linezolid | sublineage | 0 (0, 0.000) | 0 (0, 0.000) | 1 (5, 0.025) | 2 (23, 0.114) |
| moxifloxacin | main lineage | 4 (10, 0.014) | 5 (27, 0.037) | 3 (138, 0.188) | 3 (447, 0.610) |
| moxifloxacin | sublineage | 3 (4, 0.005) | 3 (16, 0.022) | 2 (8, 0.011) | 7 (579, 0.790) |
| ofloxacin | main lineage | 4 (6, 0.017) | 7 (269, 0.758) | 1 (47, 0.132) | 0 (0, 0.000) |
| ofloxacin | sublineage | 4 (6, 0.017) | 0 (0, 0.000) | 3 (28, 0.079) | 5 (289, 0.814) |
| para-aminosalicylic acid | main lineage | 0 (0, 0.000) | 0 (0, 0.000) | 0 (0, 0.000) | 0 (0, 0.000) |
| para-aminosalicylic acid | sublineage | 0 (0, 0.000) | 0 (0, 0.000) | 0 (0, 0.000) | 0 (0, 0.000) |
| prothionamide | main lineage | 13 (53, 0.576) | 1 (19, 0.207) | 0 (0, 0.000) | 0 (0, 0.000) |
| prothionamide | sublineage | 13 (53, 0.576) | 1 (19, 0.207) | 0 (0, 0.000) | 0 (0, 0.000) |
| pyrazinamide | main lineage | 78 (196, 0.514) | 15 (111, 0.291) | 1 (14, 0.037) | 0 (0, 0.000) |
| pyrazinamide | sublineage | 76 (190, 0.499) | 14 (76, 0.199) | 2 (22, 0.058) | 2 (33, 0.087) |
| rifabutin | main lineage | 16 (21, 0.032) | 2 (6, 0.009) | 4 (33, 0.050) | 5 (595, 0.895) |
| rifabutin | sublineage | 15 (19, 0.029) | 3 (8, 0.012) | 1 (3, 0.005) | 8 (624, 0.938) |
| rifampin | main lineage | 20 (26, 0.019) | 9 (32, 0.023) | 5 (151, 0.111) | 6 (1,107, 0.811) |
| rifampin | sublineage | 19 (24, 0.018) | 8 (23, 0.017) | 2 (11, 0.008) | 11 (1,257, 0.921) |
| streptomycin | main lineage | 9 (71, 0.132) | 1 (8, 0.015) | 4 (396, 0.735) | 1 (23, 0.043) |
| streptomycin | sublineage | 8 (68, 0.126) | 1 (8, 0.015) | 2 (32, 0.059) | 4 (392, 0.727) |

## Genes

Per gene as the WHO catalogue names it, the count of its graded variants seen among the resistant genomes, and the resistant genomes with one. A gene seen in 1 to 4 resistant genomes of an antibiotic is not named: such genes share one row, `other genes`. A variant's own name enters no table. Of the two largest counts, *ethA* holds 129 of ethionamide's 136 graded variants and *pncA* all 94 of pyrazinamide's.

| Antibiotic | Gene | Graded variants seen | Resistant genomes with one |
|---|---|---|---|
| amikacin | *rrs* | 3 | 415 |
| amikacin | *eis* | 1 | 17 |
| capreomycin | *rrs* | 3 | 325 |
| capreomycin | *tlyA* | 7 | 22 |
| ciprofloxacin | *gyrA* | 7 | 64 |
| ciprofloxacin | other genes | 1 | 1 |
| clofazimine | *mmpR5* | 5 | 5 |
| delamanid | other genes | 3 | 4 |
| ethambutol | *embB* | 12 | 991 |
| ethionamide | *inhA* | 7 | 439 |
| ethionamide | *ethA* | 129 | 271 |
| isoniazid | *katG* | 10 | 1,192 |
| isoniazid | *inhA* | 7 | 528 |
| kanamycin | *rrs* | 2 | 388 |
| kanamycin | *eis* | 4 | 107 |
| levofloxacin | *gyrA* | 10 | 432 |
| levofloxacin | *gyrB* | 6 | 19 |
| linezolid | *rplC* | 1 | 18 |
| linezolid | *rrl* | 2 | 6 |
| moxifloxacin | *gyrA* | 10 | 585 |
| moxifloxacin | *gyrB* | 5 | 21 |
| ofloxacin | *gyrA* | 8 | 314 |
| ofloxacin | *gyrB* | 4 | 8 |
| prothionamide | *ethA* | 10 | 41 |
| prothionamide | *inhA* | 4 | 25 |
| pyrazinamide | *pncA* | 94 | 321 |
| rifabutin | *rpoB* | 27 | 645 |
| rifampin | *rpoB* | 40 | 1,293 |
| streptomycin | *rpsL* | 2 | 291 |
| streptomycin | *rrs* | 3 | 127 |
| streptomycin | *gid* | 10 | 84 |

## Seen after

These points were seen in the tables after the reading. The rules did not ask for them. They change nothing above and decide nothing.

**The share with no graded variant differs between main lineages.** For rifampin it is 3 of 858 in lineage2 (0.003), 18 of 50 in lineage1 (0.360) and 35 of 160 in lineage3 (0.219). For levofloxacin it is 39 of 372 (0.105), 91 of 97 (0.938) and 125 of 173 (0.723). Nobody looked into it further.

**Calls that join two lineages.** Four main-lineage groups join two lineages as the tool writes them, and each holds 1 to 4 resistant genomes of an antibiotic. An assembly is one consensus sequence (section Limits). These calls were not looked into.

## Limits

**What a table can say.** A table says which graded variants of the catalogue the resistant genomes carry. It does not say what causes resistance in a genome, and it names no new marker. A gene is named only as the WHO catalogue names it, in 5 or more resistant genomes of the antibiotic.

**Deposits.** The genomes are BV-BRC deposits chosen for having lab results, capped at 2,000 a species. A share is no prevalence, and a curve describes the deposits of one public database, not the variants of any population.

**One database.** A split by lineage is a retrospective split inside one database, never external validation.

**Assemblies.** The tool is built for sequencing reads, and the project holds assemblies. An assembly is one consensus sequence and cannot show a mixed population, and a target region the assembly does not cover in full gives no variant there. On these genomes the TBProfiler comparison (`docs/tbprofiler-comparison.md`) counted the genome-drug rows with no call where an assembly did not cover a target region in full.

**Drugs with no graded entry.** The shipped database grades no entry for cycloserine or para-aminosalicylic acid (`docs/tbprofiler-comparison.md`). For those two the curve is 0 at every step and every resistant genome counts as one with no graded variant, whatever the genomes hold.

**Proxies.** Ciprofloxacin, ofloxacin, rifabutin and prothionamide count another drug's graded variants among their own resistant genomes.

**The catalogue's version.** The grades are those of `tbdb` `3ed1cb99`. Another version of the catalogue can grade other variants. The catalogue was built from public isolates with recorded phenotypes, and some of these genomes may be among them; the overlap was not checked.

**No graded variant.** A resistant genome with no graded variant is counted, and nothing more is said of it.

**Labels.** The reports in set M wait for the maintainer's label rule. They are marked.

**Markers.** No marker is found, named or ranked. The page holds aggregates over antibiotics, lineage groups and genes only: no genome id, no BioProject accession, no variant name, no sequence and no position.

## Method

**Tool run.** On one machine, from the TBProfiler comparison's own scripts: the environment from its explicit package list, each package held to its hash, with micromamba held to its published sha256; the annotation file the package ships copied to the place the tool reads, as that comparison did; the reference genome as a control; then `run_tbprofiler.sh` over the cohort list `cohort.py` builds from the reports. The cohort list equals that comparison's byte for byte.

**Checks of the run.** Before any table was read, `scripts/mechanism/curve_tables.py check` read the run's folder and printed counts only. All six checks passed. `checks.txt` in `reference/lineage-curve/` holds its printout.

| Check | Count |
|---|---|
| planned | 2,000 |
| results | 2,000 |
| failed lines | 0 |
| failed, by reason | none |
| wall time | 2,039 s, from 2026-11-24T21:41:01Z to 2026-11-24T22:15:00Z |

The TBProfiler comparison's own adapter, run over this run's result files, wrote a calls file and a database grades file with the sha256 that `reference/tbprofiler-comparison/inputs.sha256` gives for `calls.tsv` and `database-grades-by-drug.tsv`: the run again gives the calls of 14 Nov 2026 byte for byte.

**Reading.** The reading ran once, with the command of `docs/mechanism-split.md`, then the scan of its output, then the table script, then the scan of the record. `tb_curve.py` did not stop, and its reader was not changed. `reading.txt` lists the first three commands with the sha256 of each output.

```
python3 -I scripts/mechanism/tb_curve.py --species mycobacterium_tuberculosis --reports REPORTS --lists scripts/c3-inputs/lists.tsv --tool-dir RUN_DIR --bed TBDB.bed --out m4
python3 -I scripts/mechanism/scan.py --ids SPECIES.json m4
python3 -I scripts/mechanism/scan.py --ids ALL_SPECIES.json m4
python3 -I scripts/mechanism/curve_tables.py write READING scripts/c3-inputs/lists.tsv TBDB.bed reference/lineage-curve
python3 -I scripts/mechanism/scan.py --ids SPECIES.json reference/lineage-curve
python3 -I scripts/mechanism/scan.py --ids ALL_SPECIES.json reference/lineage-curve
```

REPORTS is the folder of the release's reports, RUN_DIR the tool run's folder, TBDB.bed the database's gene file, READING the folder that holds `m4/`, `reading.txt` and `checks.txt`, SPECIES.json the release's genome ids of the species and ALL_SPECIES.json every genome id of the release. Each scan found no hit.

**Tables.** `scripts/mechanism/curve_tables.py write` checked the reading against seven gates and wrote `reference/lineage-curve/`. All seven passed:

| Gate | Holds |
|---|---|
| 1. layout | the reading holds reading.txt, checks.txt and tb_curve.py's five files under m4/, and nothing more under m4/ |
| 2. versions | checks.txt counts every result file naming TBProfiler 6.6.3 and tbdb 3ed1cb99 and none naming another, and m4-summary.json names the same two |
| 3. accounting | result files plus failed lines equal the genomes planned, with no planned genome left with neither, none with both and no result outside the plan; m4-summary.json gives the same counts |
| 4. datasets | the summary's endpoints and reports not scored are exactly the species' reports in lists.tsv, named as lists.tsv names them, and every row of every table carries in_M as lists.tsv gives it |
| 5. floors | a share stands only over 5 or more resistant genomes, and no curve over under 5 |
| 6. names | no file names a variant, a position or a sequence: the gene column holds only gene names of tbdb.bed or `other genes` |
| 7. genes | a gene is named only where 5 or more resistant genomes of the antibiotic carry it |

**Tests.** `tests/mechanism/test_curve_tables.py` holds 20 SYNTHETIC tests of the table script. `tests/mechanism/test_curve_faults.py` runs `scripts/mechanism/curve_faults.py`, which plants each fault below in the table script in turn and runs the test meant to catch it. `tests/mechanism/test_lineage_page.py` holds 12 tests of this page against `reference/lineage-curve/`.

| Fault planted | Test that caught it |
|---|---|
| C1. A share written over under 5 resistant genomes | `Gates.test_a_share_under_5_resistant_genomes_is_refused` |
| C2. A variant's own name in an output | `Gates.test_a_variants_own_name_in_an_output_is_refused` |
| C3. A gene seen in 1 to 4 resistant genomes named | `Gates.test_a_gene_in_1_to_4_resistant_genomes_named_is_refused` |
| C4. A dataset in M left unmarked | `Gates.test_a_dataset_of_m_left_unmarked_is_refused` |
| C5. A result file of another tool or database version read | `Gates.test_a_result_of_another_tool_version_is_refused` |
| C6. A planned genome with neither a result nor a line in failed.tsv passed | `Gates.test_a_planned_genome_with_neither_a_result_nor_a_failed_line_is_refused` |
| C7. A curve written over under 5 resistant genomes | `Gates.test_a_curve_over_under_5_is_refused` |
| C8. The run check reads a result file of another tool version as the pinned one | `Check.test_a_result_of_another_tool_version_fails_the_check` |
| C9. The run check passes a planned genome with neither a result nor a failed line | `Check.test_a_planned_genome_with_neither_fails_the_check` |
| C10. The run check reads a result file of another database commit as the pinned one | `Check.test_a_result_of_another_database_commit_fails_the_check` |
| C11. A summary of another database commit read | `Gates.test_a_summary_of_another_database_commit_is_refused` |
| C12. A convergence share written where under 5 resistant genomes have a group | `Gates.test_a_convergence_share_under_5_resistant_genomes_with_a_group_is_refused` |

All 12 were caught.

**Files.** The reading's tables are in `reference/lineage-curve/`. `scripts/mechanism/curve_tables.py print` prints this page's tables from those files, rounded to 3 decimals.

## Datasets

One row per report of the species, species and antibiotic as `scripts/c3-inputs/lists.tsv` names them.

| Species | Antibiotic | Set | Kind | Tool drugs | Resistant with a result | Resistant with no result | No main lineage | No cluster |
|---|---|---|---|---|---|---|---|---|
| *Mycobacterium tuberculosis* | amikacin | not M | primary | amikacin | 548 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | amoxicillin | not M | none | none | none | none | none | none |
| *Mycobacterium tuberculosis* | capreomycin | not M | primary | capreomycin | 422 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | ciprofloxacin | not M | proxy | levofloxacin, moxifloxacin | 73 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | clofazimine | M | primary | clofazimine | 193 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | cycloserine | not M | primary | cycloserine | 119 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | delamanid | not M | primary | delamanid | 448 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | ethambutol | not M | primary | ethambutol | 1,338 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | ethionamide | not M | primary | ethionamide | 915 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | isoniazid | M | primary | isoniazid | 1,388 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | kanamycin | not M | primary | kanamycin | 617 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | levofloxacin | not M | primary | levofloxacin | 753 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | linezolid | not M | primary | linezolid | 202 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | moxifloxacin | not M | primary | moxifloxacin | 733 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | nicotinamide | not M | none | none | none | none | none | none |
| *Mycobacterium tuberculosis* | ofloxacin | not M | proxy | levofloxacin, moxifloxacin | 355 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | para-aminosalicylic acid | not M | primary | para-aminosalicylic acid | 429 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | prothionamide | not M | proxy | ethionamide | 92 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | pyrazinamide | not M | primary | pyrazinamide | 381 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | rifabutin | not M | proxy | rifampicin | 665 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | rifampin | not M | primary | rifampicin | 1,365 | 0 | 0 | 0 |
| *Mycobacterium tuberculosis* | streptomycin | not M | primary | streptomycin | 539 | 0 | 0 | 0 |
