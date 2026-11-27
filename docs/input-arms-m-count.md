# Count tables of M

The script `scripts/c3-inputs/m_page.py page` writes this page whole from the record, `reference/input-arms-m/`, and no one edits it by hand.

Datasets of M read: 35 of 35 in the first reading, 12 of 35 in the second (3 left out under the floor, 20 for a short list). The left-out page, `docs/input-arms-m-left-out.md`, names each dataset left out.

**Both readings.** Per arm, conditions 1 to 5 of ruling 15 on the count over all datasets.

| Arm | Condition | First reading | Second reading | Both |
|---|---|---|---|---|
| unitigs | 1 and 2 | not met | not met | not met |
| unitigs | 3 | met | met | met |
| unitigs | 4 | met | met | met |
| unitigs | 5 | not met | not met | not met |
| K31 plus k-mer counts | 1 and 2 | not met | not met | not met |
| K31 plus k-mer counts | 3 | met | met | met |
| K31 plus k-mer counts | 4 | met | met | met |
| K31 plus k-mer counts | 5 | not met | not met | not met |
| K31 plus amino-acid k-mers | 1 and 2 | not met | not met | not met |
| K31 plus amino-acid k-mers | 3 | met | met | met |
| K31 plus amino-acid k-mers | 4 | met | met | met |
| K31 plus amino-acid k-mers | 5 | met | met | met |
| the AMRFinderPlus program's calls alone | 1 and 2 | met | met | met |
| the AMRFinderPlus program's calls alone | 3 | met | met | met |
| the AMRFinderPlus program's calls alone | 4 | met | met | met |
| the AMRFinderPlus program's calls alone | 5 | met | met | met |
| K31 plus the AMRFinderPlus program's calls | 1 and 2 | not met | not met | not met |
| K31 plus the AMRFinderPlus program's calls | 3 | met | met | met |
| K31 plus the AMRFinderPlus program's calls | 4 | met | met | met |
| K31 plus the AMRFinderPlus program's calls | 5 | not met | not met | not met |

**Datasets read.** The datasets of M by species: in the second reading, left out under the floor, or left out for a short list.

| Species | Datasets of M | In the second reading | Left out, under the floor | Left out, short list |
|---|---|---|---|---|
| *Acinetobacter baumannii* | 11 | 7 of 11 | 0 of 11 | 4 of 11 |
| *Enterococcus faecium* | 2 | 1 of 2 | 0 of 2 | 1 of 2 |
| *Escherichia coli* | 7 | 0 of 7 | 0 of 7 | 7 of 7 |
| *Klebsiella pneumoniae* | 6 | 0 of 6 | 0 of 6 | 6 of 6 |
| *Mycobacterium tuberculosis* | 2 | 1 of 2 | 0 of 2 | 1 of 2 |
| *Neisseria gonorrhoeae* | 1 | 1 of 1 | 0 of 1 | 0 of 1 |
| *Pseudomonas aeruginosa* | 1 | 0 of 1 | 0 of 1 | 1 of 1 |
| *Staphylococcus aureus* | 1 | 1 of 1 | 0 of 1 | 0 of 1 |
| *Streptococcus pneumoniae* | 4 | 1 of 4 | 3 of 4 | 0 of 4 |
| All 9 species | 35 | 12 of 35 | 3 of 35 | 20 of 35 |

**All datasets.** The first run's datasets and M's, on the lineage split and balanced accuracy. Condition 1 is 1 where the interval lies wholly above 0. Condition 2 is 1 where better beats worse and the Holm p, the sign test's p after Holm's correction, is under 0.05. Each mark comes from the record's full figures, not the 3 decimals shown.

| Reading | Arm | Learner | Better | Worse | Same | Median | Interval | Sign p | Holm p | Constant | Condition 1 | Condition 2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| first | unitigs | `scm` | 26 of 267 | 32 of 267 | 209 of 267 | 0.000 | -0.001 to 0.000 | 0.821 | 1.000 | 0 | 0 | 0 |
| first | unitigs | `cart` | 26 of 267 | 23 of 267 | 218 of 267 | 0.000 | 0.000 to 0.000 | 0.388 | 1.000 | 0 | 0 | 0 |
| first | unitigs | `rf` | 11 of 267 | 16 of 267 | 240 of 267 | 0.000 | -0.001 to 0.000 | 0.876 | 1.000 | 0 | 0 | 0 |
| first | unitigs | `gbt` | 37 of 267 | 16 of 267 | 214 of 267 | 0.000 | 0.000 to 0.001 | 0.003 | 0.058 | 0 | 0 | 0 |
| first | unitigs | `lasso` | 12 of 267 | 13 of 267 | 242 of 267 | 0.000 | -0.001 to 0.000 | 0.655 | 1.000 | 0 | 0 | 0 |
| first | K31 plus k-mer counts | `scm` | 23 of 267 | 21 of 267 | 223 of 267 | 0.000 | 0.000 to 0.000 | 0.440 | 1.000 | 0 | 0 | 0 |
| first | K31 plus k-mer counts | `cart` | 30 of 267 | 26 of 267 | 211 of 267 | 0.000 | 0.000 to 0.000 | 0.344 | 1.000 | 0 | 0 | 0 |
| first | K31 plus k-mer counts | `rf` | 18 of 267 | 14 of 267 | 235 of 267 | 0.000 | 0.000 to 0.001 | 0.298 | 1.000 | 0 | 0 | 0 |
| first | K31 plus k-mer counts | `gbt` | 40 of 267 | 21 of 267 | 206 of 267 | 0.001 | 0.000 to 0.002 | 0.010 | 0.204 | 0 | 1 | 0 |
| first | K31 plus k-mer counts | `lasso` | 21 of 267 | 13 of 267 | 233 of 267 | 0.000 | 0.000 to 0.000 | 0.115 | 1.000 | 0 | 0 | 0 |
| first | K31 plus amino-acid k-mers | `scm` | 58 of 267 | 26 of 267 | 183 of 267 | 0.002 | 0.000 to 0.004 | < 0.001 | 0.007 | 0 | 0 | 1 |
| first | K31 plus amino-acid k-mers | `cart` | 61 of 267 | 28 of 267 | 178 of 267 | 0.000 | 0.000 to 0.003 | < 0.001 | 0.007 | 0 | 0 | 1 |
| first | K31 plus amino-acid k-mers | `rf` | 45 of 267 | 10 of 267 | 212 of 267 | 0.002 | 0.001 to 0.005 | < 0.001 | < 0.001 | 0 | 1 | 1 |
| first | K31 plus amino-acid k-mers | `gbt` | 57 of 267 | 20 of 267 | 190 of 267 | 0.003 | 0.000 to 0.005 | < 0.001 | < 0.001 | 0 | 1 | 1 |
| first | K31 plus amino-acid k-mers | `lasso` | 50 of 267 | 20 of 267 | 197 of 267 | 0.001 | 0.000 to 0.003 | < 0.001 | 0.005 | 0 | 0 | 1 |
| first | the AMRFinderPlus program's calls alone | `scm` | 79 of 236 | 62 of 236 | 95 of 236 | 0.004 | 0.000 to 0.012 | 0.089 | 1.000 | 0 | 0 | 0 |
| first | the AMRFinderPlus program's calls alone | `cart` | 101 of 236 | 51 of 236 | 84 of 236 | 0.010 | 0.005 to 0.020 | < 0.001 | < 0.001 | 0 | 1 | 1 |
| first | the AMRFinderPlus program's calls alone | `rf` | 145 of 236 | 32 of 236 | 59 of 236 | 0.040 | 0.030 to 0.054 | < 0.001 | < 0.001 | 0 | 1 | 1 |
| first | the AMRFinderPlus program's calls alone | `gbt` | 127 of 236 | 28 of 236 | 81 of 236 | 0.022 | 0.015 to 0.033 | < 0.001 | < 0.001 | 0 | 1 | 1 |
| first | the AMRFinderPlus program's calls alone | `lasso` | 120 of 236 | 36 of 236 | 80 of 236 | 0.020 | 0.012 to 0.027 | < 0.001 | < 0.001 | 0 | 1 | 1 |
| first | K31 plus the AMRFinderPlus program's calls | `scm` | 7 of 236 | 9 of 236 | 220 of 236 | 0.000 | 0.000 to 0.000 | 0.773 | 1.000 | 0 | 0 | 0 |
| first | K31 plus the AMRFinderPlus program's calls | `cart` | 6 of 236 | 6 of 236 | 224 of 236 | 0.000 | 0.000 to 0.000 | 0.613 | 1.000 | 0 | 0 | 0 |
| first | K31 plus the AMRFinderPlus program's calls | `rf` | 5 of 236 | 1 of 236 | 230 of 236 | 0.000 | 0.000 to 0.000 | 0.109 | 1.000 | 0 | 0 | 0 |
| first | K31 plus the AMRFinderPlus program's calls | `gbt` | 13 of 236 | 11 of 236 | 212 of 236 | 0.000 | 0.000 to 0.000 | 0.419 | 1.000 | 0 | 0 | 0 |
| first | K31 plus the AMRFinderPlus program's calls | `lasso` | 1 of 236 | 2 of 236 | 233 of 236 | 0.000 | 0.000 to 0.000 | 0.875 | 1.000 | 0 | 0 | 0 |
| first | the twin, the default at seeds 45, 46 and 47 | `scm` | 33 of 267 | 31 of 267 | 203 of 267 | 0.000 | 0.000 to 0.001 | 0.450 | 1.000 | 0 | 0 | 0 |
| first | the twin, the default at seeds 45, 46 and 47 | `cart` | 33 of 267 | 34 of 267 | 200 of 267 | 0.000 | 0.000 to 0.000 | 0.596 | 1.000 | 0 | 0 | 0 |
| first | the twin, the default at seeds 45, 46 and 47 | `rf` | 20 of 267 | 24 of 267 | 223 of 267 | -0.001 | -0.003 to 0.000 | 0.774 | 1.000 | 0 | 0 | 0 |
| first | the twin, the default at seeds 45, 46 and 47 | `gbt` | 48 of 267 | 37 of 267 | 182 of 267 | 0.000 | 0.000 to 0.001 | 0.139 | 1.000 | 0 | 0 | 0 |
| first | the twin, the default at seeds 45, 46 and 47 | `lasso` | 32 of 267 | 21 of 267 | 214 of 267 | 0.000 | -0.001 to 0.000 | 0.084 | 1.000 | 0 | 0 | 0 |
| second | unitigs | `scm` | 21 of 244 | 30 of 244 | 193 of 244 | 0.000 | -0.001 to 0.000 | 0.920 | 1.000 | 0 | 0 | 0 |
| second | unitigs | `cart` | 24 of 244 | 19 of 244 | 201 of 244 | 0.000 | 0.000 to 0.000 | 0.271 | 1.000 | 0 | 0 | 0 |
| second | unitigs | `rf` | 11 of 244 | 13 of 244 | 220 of 244 | 0.000 | -0.001 to 0.000 | 0.729 | 1.000 | 0 | 0 | 0 |
| second | unitigs | `gbt` | 34 of 244 | 15 of 244 | 195 of 244 | 0.000 | 0.000 to 0.001 | 0.005 | 0.099 | 0 | 0 | 0 |
| second | unitigs | `lasso` | 12 of 244 | 12 of 244 | 220 of 244 | 0.000 | -0.001 to 0.000 | 0.581 | 1.000 | 0 | 0 | 0 |
| second | K31 plus k-mer counts | `scm` | 21 of 244 | 20 of 244 | 203 of 244 | 0.000 | 0.000 to 0.000 | 0.500 | 1.000 | 0 | 0 | 0 |
| second | K31 plus k-mer counts | `cart` | 28 of 244 | 26 of 244 | 190 of 244 | 0.000 | 0.000 to 0.000 | 0.446 | 1.000 | 0 | 0 | 0 |
| second | K31 plus k-mer counts | `rf` | 18 of 244 | 12 of 244 | 214 of 244 | 0.000 | 0.000 to 0.000 | 0.181 | 1.000 | 0 | 0 | 0 |
| second | K31 plus k-mer counts | `gbt` | 38 of 244 | 18 of 244 | 188 of 244 | 0.001 | 0.000 to 0.002 | 0.005 | 0.105 | 0 | 1 | 0 |
| second | K31 plus k-mer counts | `lasso` | 20 of 244 | 11 of 244 | 213 of 244 | 0.000 | 0.000 to 0.000 | 0.075 | 1.000 | 0 | 0 | 0 |
| second | K31 plus amino-acid k-mers | `scm` | 51 of 244 | 24 of 244 | 169 of 244 | 0.001 | 0.000 to 0.004 | 0.001 | 0.027 | 0 | 0 | 1 |
| second | K31 plus amino-acid k-mers | `cart` | 56 of 244 | 26 of 244 | 162 of 244 | 0.000 | 0.000 to 0.003 | < 0.001 | 0.014 | 0 | 0 | 1 |
| second | K31 plus amino-acid k-mers | `rf` | 41 of 244 | 9 of 244 | 194 of 244 | 0.002 | 0.001 to 0.005 | < 0.001 | < 0.001 | 0 | 1 | 1 |
| second | K31 plus amino-acid k-mers | `gbt` | 53 of 244 | 17 of 244 | 174 of 244 | 0.002 | 0.001 to 0.005 | < 0.001 | < 0.001 | 0 | 1 | 1 |
| second | K31 plus amino-acid k-mers | `lasso` | 48 of 244 | 17 of 244 | 179 of 244 | 0.001 | 0.000 to 0.004 | < 0.001 | 0.002 | 0 | 0 | 1 |
| second | the AMRFinderPlus program's calls alone | `scm` | 75 of 214 | 55 of 214 | 84 of 214 | 0.005 | 0.000 to 0.013 | 0.048 | 0.905 | 0 | 0 | 0 |
| second | the AMRFinderPlus program's calls alone | `cart` | 91 of 214 | 43 of 214 | 80 of 214 | 0.010 | 0.005 to 0.020 | < 0.001 | < 0.001 | 0 | 1 | 1 |
| second | the AMRFinderPlus program's calls alone | `rf` | 129 of 214 | 29 of 214 | 56 of 214 | 0.036 | 0.026 to 0.053 | < 0.001 | < 0.001 | 0 | 1 | 1 |
| second | the AMRFinderPlus program's calls alone | `gbt` | 113 of 214 | 27 of 214 | 74 of 214 | 0.021 | 0.013 to 0.031 | < 0.001 | < 0.001 | 0 | 1 | 1 |
| second | the AMRFinderPlus program's calls alone | `lasso` | 107 of 214 | 33 of 214 | 74 of 214 | 0.020 | 0.010 to 0.027 | < 0.001 | < 0.001 | 0 | 1 | 1 |
| second | K31 plus the AMRFinderPlus program's calls | `scm` | 6 of 214 | 9 of 214 | 199 of 214 | 0.000 | 0.000 to 0.000 | 0.849 | 1.000 | 0 | 0 | 0 |
| second | K31 plus the AMRFinderPlus program's calls | `cart` | 5 of 214 | 5 of 214 | 204 of 214 | 0.000 | 0.000 to 0.000 | 0.623 | 1.000 | 0 | 0 | 0 |
| second | K31 plus the AMRFinderPlus program's calls | `rf` | 4 of 214 | 2 of 214 | 208 of 214 | 0.000 | 0.000 to 0.000 | 0.344 | 1.000 | 0 | 0 | 0 |
| second | K31 plus the AMRFinderPlus program's calls | `gbt` | 11 of 214 | 12 of 214 | 191 of 214 | 0.000 | 0.000 to 0.000 | 0.661 | 1.000 | 0 | 0 | 0 |
| second | K31 plus the AMRFinderPlus program's calls | `lasso` | 1 of 214 | 1 of 214 | 212 of 214 | 0.000 | 0.000 to 0.000 | 0.750 | 1.000 | 0 | 0 | 0 |
| second | the twin, the default at seeds 45, 46 and 47 | `scm` | 30 of 244 | 30 of 244 | 184 of 244 | 0.000 | 0.000 to 0.001 | 0.551 | 1.000 | 0 | 0 | 0 |
| second | the twin, the default at seeds 45, 46 and 47 | `cart` | 28 of 244 | 32 of 244 | 184 of 244 | 0.000 | -0.001 to 0.000 | 0.741 | 1.000 | 0 | 0 | 0 |
| second | the twin, the default at seeds 45, 46 and 47 | `rf` | 18 of 244 | 24 of 244 | 202 of 244 | -0.001 | -0.003 to 0.000 | 0.860 | 1.000 | 0 | 0 | 0 |
| second | the twin, the default at seeds 45, 46 and 47 | `gbt` | 42 of 244 | 35 of 244 | 167 of 244 | 0.000 | 0.000 to 0.001 | 0.247 | 1.000 | 0 | 0 | 0 |
| second | the twin, the default at seeds 45, 46 and 47 | `lasso` | 28 of 244 | 20 of 244 | 196 of 244 | 0.000 | -0.001 to 0.000 | 0.156 | 1.000 | 0 | 0 | 0 |

**All datasets, recall.** The first run's datasets and M's, on the lineage split and recall.

| Reading | Arm | Learner | Better | Worse | Same | Median | Interval |
|---|---|---|---|---|---|---|---|
| first | unitigs | `scm` | 38 of 267 | 45 of 267 | 184 of 267 | 0.000 | 0.000 to 0.000 |
| first | unitigs | `cart` | 42 of 267 | 45 of 267 | 180 of 267 | 0.000 | 0.000 to 0.000 |
| first | unitigs | `rf` | 18 of 267 | 18 of 267 | 231 of 267 | 0.000 | 0.000 to 0.000 |
| first | unitigs | `gbt` | 49 of 267 | 36 of 267 | 182 of 267 | 0.000 | 0.000 to 0.000 |
| first | unitigs | `lasso` | 25 of 267 | 22 of 267 | 220 of 267 | 0.000 | 0.000 to 0.000 |
| first | K31 plus k-mer counts | `scm` | 38 of 267 | 38 of 267 | 191 of 267 | 0.000 | 0.000 to 0.000 |
| first | K31 plus k-mer counts | `cart` | 40 of 267 | 41 of 267 | 186 of 267 | 0.000 | 0.000 to 0.000 |
| first | K31 plus k-mer counts | `rf` | 38 of 267 | 17 of 267 | 212 of 267 | 0.000 | 0.000 to 0.000 |
| first | K31 plus k-mer counts | `gbt` | 52 of 267 | 36 of 267 | 179 of 267 | 0.000 | 0.000 to 0.001 |
| first | K31 plus k-mer counts | `lasso` | 35 of 267 | 25 of 267 | 207 of 267 | 0.000 | 0.000 to 0.000 |
| first | K31 plus amino-acid k-mers | `scm` | 64 of 267 | 51 of 267 | 152 of 267 | 0.000 | 0.000 to 0.005 |
| first | K31 plus amino-acid k-mers | `cart` | 77 of 267 | 42 of 267 | 148 of 267 | 0.000 | 0.000 to 0.001 |
| first | K31 plus amino-acid k-mers | `rf` | 53 of 267 | 23 of 267 | 191 of 267 | 0.001 | 0.000 to 0.005 |
| first | K31 plus amino-acid k-mers | `gbt` | 61 of 267 | 45 of 267 | 161 of 267 | 0.000 | 0.000 to 0.003 |
| first | K31 plus amino-acid k-mers | `lasso` | 58 of 267 | 40 of 267 | 169 of 267 | 0.000 | 0.000 to 0.005 |
| first | the AMRFinderPlus program's calls alone | `scm` | 76 of 236 | 75 of 236 | 85 of 236 | 0.000 | -0.001 to 0.005 |
| first | the AMRFinderPlus program's calls alone | `cart` | 98 of 236 | 71 of 236 | 67 of 236 | 0.007 | 0.000 to 0.019 |
| first | the AMRFinderPlus program's calls alone | `rf` | 123 of 236 | 42 of 236 | 71 of 236 | 0.024 | 0.013 to 0.045 |
| first | the AMRFinderPlus program's calls alone | `gbt` | 108 of 236 | 55 of 236 | 73 of 236 | 0.012 | 0.003 to 0.027 |
| first | the AMRFinderPlus program's calls alone | `lasso` | 101 of 236 | 55 of 236 | 80 of 236 | 0.014 | 0.005 to 0.019 |
| first | K31 plus the AMRFinderPlus program's calls | `scm` | 19 of 236 | 11 of 236 | 206 of 236 | 0.000 | 0.000 to 0.000 |
| first | K31 plus the AMRFinderPlus program's calls | `cart` | 11 of 236 | 9 of 236 | 216 of 236 | 0.000 | 0.000 to 0.000 |
| first | K31 plus the AMRFinderPlus program's calls | `rf` | 12 of 236 | 9 of 236 | 215 of 236 | 0.000 | 0.000 to 0.000 |
| first | K31 plus the AMRFinderPlus program's calls | `gbt` | 22 of 236 | 20 of 236 | 194 of 236 | 0.000 | 0.000 to 0.000 |
| first | K31 plus the AMRFinderPlus program's calls | `lasso` | 3 of 236 | 4 of 236 | 229 of 236 | 0.000 | 0.000 to 0.000 |
| first | the twin, the default at seeds 45, 46 and 47 | `scm` | 50 of 267 | 49 of 267 | 168 of 267 | 0.000 | 0.000 to 0.000 |
| first | the twin, the default at seeds 45, 46 and 47 | `cart` | 50 of 267 | 51 of 267 | 166 of 267 | 0.000 | 0.000 to 0.000 |
| first | the twin, the default at seeds 45, 46 and 47 | `rf` | 26 of 267 | 42 of 267 | 199 of 267 | 0.000 | -0.001 to 0.000 |
| first | the twin, the default at seeds 45, 46 and 47 | `gbt` | 57 of 267 | 53 of 267 | 157 of 267 | 0.000 | 0.000 to 0.000 |
| first | the twin, the default at seeds 45, 46 and 47 | `lasso` | 46 of 267 | 37 of 267 | 184 of 267 | 0.000 | 0.000 to 0.001 |
| second | unitigs | `scm` | 34 of 244 | 41 of 244 | 169 of 244 | 0.000 | 0.000 to 0.000 |
| second | unitigs | `cart` | 38 of 244 | 37 of 244 | 169 of 244 | 0.000 | 0.000 to 0.000 |
| second | unitigs | `rf` | 17 of 244 | 15 of 244 | 212 of 244 | 0.000 | 0.000 to 0.000 |
| second | unitigs | `gbt` | 46 of 244 | 33 of 244 | 165 of 244 | 0.000 | 0.000 to 0.000 |
| second | unitigs | `lasso` | 25 of 244 | 18 of 244 | 201 of 244 | 0.000 | 0.000 to 0.000 |
| second | K31 plus k-mer counts | `scm` | 36 of 244 | 35 of 244 | 173 of 244 | 0.000 | 0.000 to 0.000 |
| second | K31 plus k-mer counts | `cart` | 36 of 244 | 38 of 244 | 170 of 244 | 0.000 | 0.000 to 0.000 |
| second | K31 plus k-mer counts | `rf` | 35 of 244 | 14 of 244 | 195 of 244 | 0.000 | 0.000 to 0.000 |
| second | K31 plus k-mer counts | `gbt` | 47 of 244 | 32 of 244 | 165 of 244 | 0.000 | 0.000 to 0.001 |
| second | K31 plus k-mer counts | `lasso` | 32 of 244 | 22 of 244 | 190 of 244 | 0.000 | 0.000 to 0.000 |
| second | K31 plus amino-acid k-mers | `scm` | 59 of 244 | 48 of 244 | 137 of 244 | 0.000 | 0.000 to 0.005 |
| second | K31 plus amino-acid k-mers | `cart` | 71 of 244 | 38 of 244 | 135 of 244 | 0.000 | 0.000 to 0.001 |
| second | K31 plus amino-acid k-mers | `rf` | 50 of 244 | 21 of 244 | 173 of 244 | 0.001 | 0.000 to 0.004 |
| second | K31 plus amino-acid k-mers | `gbt` | 56 of 244 | 39 of 244 | 149 of 244 | 0.000 | 0.000 to 0.003 |
| second | K31 plus amino-acid k-mers | `lasso` | 55 of 244 | 34 of 244 | 155 of 244 | 0.000 | 0.000 to 0.005 |
| second | the AMRFinderPlus program's calls alone | `scm` | 72 of 214 | 67 of 214 | 75 of 214 | 0.000 | -0.001 to 0.007 |
| second | the AMRFinderPlus program's calls alone | `cart` | 86 of 214 | 65 of 214 | 63 of 214 | 0.005 | 0.000 to 0.019 |
| second | the AMRFinderPlus program's calls alone | `rf` | 110 of 214 | 40 of 214 | 64 of 214 | 0.023 | 0.012 to 0.043 |
| second | the AMRFinderPlus program's calls alone | `gbt` | 95 of 214 | 50 of 214 | 69 of 214 | 0.012 | 0.003 to 0.026 |
| second | the AMRFinderPlus program's calls alone | `lasso` | 88 of 214 | 52 of 214 | 74 of 214 | 0.012 | 0.004 to 0.018 |
| second | K31 plus the AMRFinderPlus program's calls | `scm` | 20 of 214 | 10 of 214 | 184 of 214 | 0.000 | 0.000 to 0.000 |
| second | K31 plus the AMRFinderPlus program's calls | `cart` | 9 of 214 | 7 of 214 | 198 of 214 | 0.000 | 0.000 to 0.000 |
| second | K31 plus the AMRFinderPlus program's calls | `rf` | 10 of 214 | 8 of 214 | 196 of 214 | 0.000 | 0.000 to 0.000 |
| second | K31 plus the AMRFinderPlus program's calls | `gbt` | 21 of 214 | 19 of 214 | 174 of 214 | 0.000 | 0.000 to 0.000 |
| second | K31 plus the AMRFinderPlus program's calls | `lasso` | 3 of 214 | 3 of 214 | 208 of 214 | 0.000 | 0.000 to 0.000 |
| second | the twin, the default at seeds 45, 46 and 47 | `scm` | 46 of 244 | 47 of 244 | 151 of 244 | 0.000 | 0.000 to 0.000 |
| second | the twin, the default at seeds 45, 46 and 47 | `cart` | 40 of 244 | 47 of 244 | 157 of 244 | 0.000 | 0.000 to 0.000 |
| second | the twin, the default at seeds 45, 46 and 47 | `rf` | 25 of 244 | 40 of 244 | 179 of 244 | 0.000 | -0.001 to 0.000 |
| second | the twin, the default at seeds 45, 46 and 47 | `gbt` | 54 of 244 | 51 of 244 | 139 of 244 | 0.000 | 0.000 to 0.000 |
| second | the twin, the default at seeds 45, 46 and 47 | `lasso` | 42 of 244 | 35 of 244 | 167 of 244 | 0.000 | 0.000 to 0.001 |

**All datasets, study split.** The first run's datasets and M's, on the study split and balanced accuracy.

| Reading | Arm | Learner | Better | Worse | Same | Median | Interval |
|---|---|---|---|---|---|---|---|
| first | unitigs | `scm` | 10 of 69 | 9 of 69 | 50 of 69 | 0.000 | 0.000 to 0.002 |
| first | unitigs | `cart` | 9 of 69 | 6 of 69 | 54 of 69 | 0.000 | 0.000 to 0.001 |
| first | unitigs | `rf` | 2 of 69 | 3 of 69 | 64 of 69 | -0.001 | -0.002 to 0.000 |
| first | unitigs | `gbt` | 6 of 69 | 6 of 69 | 57 of 69 | 0.000 | -0.001 to 0.001 |
| first | unitigs | `lasso` | 3 of 69 | 1 of 69 | 65 of 69 | 0.000 | -0.001 to 0.003 |
| first | K31 plus k-mer counts | `scm` | 10 of 69 | 10 of 69 | 49 of 69 | 0.000 | 0.000 to 0.005 |
| first | K31 plus k-mer counts | `cart` | 13 of 69 | 15 of 69 | 41 of 69 | 0.000 | 0.000 to 0.002 |
| first | K31 plus k-mer counts | `rf` | 5 of 69 | 10 of 69 | 54 of 69 | 0.000 | -0.003 to 0.004 |
| first | K31 plus k-mer counts | `gbt` | 8 of 69 | 11 of 69 | 50 of 69 | 0.001 | 0.000 to 0.004 |
| first | K31 plus k-mer counts | `lasso` | 8 of 69 | 10 of 69 | 51 of 69 | 0.000 | -0.001 to 0.002 |
| first | K31 plus amino-acid k-mers | `scm` | 24 of 69 | 13 of 69 | 32 of 69 | 0.007 | 0.000 to 0.020 |
| first | K31 plus amino-acid k-mers | `cart` | 22 of 69 | 14 of 69 | 33 of 69 | 0.004 | 0.002 to 0.011 |
| first | K31 plus amino-acid k-mers | `rf` | 9 of 69 | 4 of 69 | 56 of 69 | 0.002 | 0.000 to 0.008 |
| first | K31 plus amino-acid k-mers | `gbt` | 19 of 69 | 9 of 69 | 41 of 69 | 0.005 | 0.002 to 0.013 |
| first | K31 plus amino-acid k-mers | `lasso` | 15 of 69 | 9 of 69 | 45 of 69 | 0.002 | 0.000 to 0.015 |
| first | the AMRFinderPlus program's calls alone | `scm` | 28 of 68 | 13 of 68 | 27 of 68 | 0.007 | -0.001 to 0.024 |
| first | the AMRFinderPlus program's calls alone | `cart` | 27 of 68 | 11 of 68 | 30 of 68 | 0.010 | 0.002 to 0.021 |
| first | the AMRFinderPlus program's calls alone | `rf` | 45 of 68 | 9 of 68 | 14 of 68 | 0.041 | 0.027 to 0.070 |
| first | the AMRFinderPlus program's calls alone | `gbt` | 33 of 68 | 10 of 68 | 25 of 68 | 0.017 | 0.009 to 0.039 |
| first | the AMRFinderPlus program's calls alone | `lasso` | 40 of 68 | 5 of 68 | 23 of 68 | 0.032 | 0.009 to 0.050 |
| first | K31 plus the AMRFinderPlus program's calls | `scm` | 3 of 68 | 3 of 68 | 62 of 68 | 0.000 | 0.000 to 0.000 |
| first | K31 plus the AMRFinderPlus program's calls | `cart` | 1 of 68 | 1 of 68 | 66 of 68 | 0.000 | 0.000 to 0.000 |
| first | K31 plus the AMRFinderPlus program's calls | `rf` | 3 of 68 | 4 of 68 | 61 of 68 | 0.000 | -0.002 to 0.000 |
| first | K31 plus the AMRFinderPlus program's calls | `gbt` | 0 of 68 | 3 of 68 | 65 of 68 | 0.000 | 0.000 to 0.001 |
| first | K31 plus the AMRFinderPlus program's calls | `lasso` | 1 of 68 | 0 of 68 | 67 of 68 | 0.000 | 0.000 to 0.000 |
| second | unitigs | `scm` | 9 of 58 | 6 of 58 | 43 of 58 | 0.000 | 0.000 to 0.002 |
| second | unitigs | `cart` | 7 of 58 | 6 of 58 | 45 of 58 | 0.000 | 0.000 to 0.002 |
| second | unitigs | `rf` | 2 of 58 | 2 of 58 | 54 of 58 | -0.001 | -0.001 to 0.000 |
| second | unitigs | `gbt` | 5 of 58 | 6 of 58 | 47 of 58 | 0.000 | -0.001 to 0.001 |
| second | unitigs | `lasso` | 2 of 58 | 1 of 58 | 55 of 58 | 0.000 | 0.000 to 0.003 |
| second | K31 plus k-mer counts | `scm` | 10 of 58 | 7 of 58 | 41 of 58 | 0.000 | -0.001 to 0.005 |
| second | K31 plus k-mer counts | `cart` | 11 of 58 | 12 of 58 | 35 of 58 | 0.000 | 0.000 to 0.002 |
| second | K31 plus k-mer counts | `rf` | 5 of 58 | 7 of 58 | 46 of 58 | 0.000 | -0.002 to 0.004 |
| second | K31 plus k-mer counts | `gbt` | 8 of 58 | 10 of 58 | 40 of 58 | 0.001 | 0.000 to 0.010 |
| second | K31 plus k-mer counts | `lasso` | 8 of 58 | 10 of 58 | 40 of 58 | 0.000 | 0.000 to 0.002 |
| second | K31 plus amino-acid k-mers | `scm` | 19 of 58 | 11 of 58 | 28 of 58 | 0.004 | 0.000 to 0.020 |
| second | K31 plus amino-acid k-mers | `cart` | 18 of 58 | 12 of 58 | 28 of 58 | 0.005 | 0.001 to 0.010 |
| second | K31 plus amino-acid k-mers | `rf` | 6 of 58 | 4 of 58 | 48 of 58 | 0.001 | 0.000 to 0.008 |
| second | K31 plus amino-acid k-mers | `gbt` | 15 of 58 | 7 of 58 | 36 of 58 | 0.006 | 0.002 to 0.013 |
| second | K31 plus amino-acid k-mers | `lasso` | 12 of 58 | 9 of 58 | 37 of 58 | 0.002 | 0.000 to 0.015 |
| second | the AMRFinderPlus program's calls alone | `scm` | 25 of 57 | 11 of 57 | 21 of 57 | 0.015 | 0.002 to 0.031 |
| second | the AMRFinderPlus program's calls alone | `cart` | 22 of 57 | 11 of 57 | 24 of 57 | 0.007 | 0.002 to 0.027 |
| second | the AMRFinderPlus program's calls alone | `rf` | 37 of 57 | 7 of 57 | 13 of 57 | 0.039 | 0.027 to 0.060 |
| second | the AMRFinderPlus program's calls alone | `gbt` | 26 of 57 | 9 of 57 | 22 of 57 | 0.015 | 0.006 to 0.039 |
| second | the AMRFinderPlus program's calls alone | `lasso` | 34 of 57 | 5 of 57 | 18 of 57 | 0.030 | 0.009 to 0.054 |
| second | K31 plus the AMRFinderPlus program's calls | `scm` | 3 of 57 | 2 of 57 | 52 of 57 | 0.000 | 0.000 to 0.000 |
| second | K31 plus the AMRFinderPlus program's calls | `cart` | 1 of 57 | 1 of 57 | 55 of 57 | 0.000 | 0.000 to 0.000 |
| second | K31 plus the AMRFinderPlus program's calls | `rf` | 3 of 57 | 2 of 57 | 52 of 57 | 0.000 | -0.001 to 0.000 |
| second | K31 plus the AMRFinderPlus program's calls | `gbt` | 0 of 57 | 3 of 57 | 54 of 57 | 0.000 | -0.003 to 0.001 |
| second | K31 plus the AMRFinderPlus program's calls | `lasso` | 1 of 57 | 0 of 57 | 56 of 57 | 0.000 | 0.000 to 0.000 |

**M alone.** The datasets of M, on the lineage split and balanced accuracy. This count decides nothing.

| Reading | Arm | Learner | Better | Worse | Same | Median | Interval | Sign p | Holm p | Constant |
|---|---|---|---|---|---|---|---|---|---|---|
| first | unitigs | `scm` | 7 of 35 | 3 of 35 | 25 of 35 | 0.000 | -0.001 to 0.004 | 0.172 | 1.000 | 0 |
| first | unitigs | `cart` | 2 of 35 | 4 of 35 | 29 of 35 | -0.001 | -0.003 to 0.000 | 0.891 | 1.000 | 0 |
| first | unitigs | `rf` | 1 of 35 | 4 of 35 | 30 of 35 | -0.001 | -0.005 to 0.001 | 0.969 | 1.000 | 0 |
| first | unitigs | `gbt` | 3 of 35 | 1 of 35 | 31 of 35 | 0.000 | -0.003 to 0.004 | 0.312 | 1.000 | 0 |
| first | unitigs | `lasso` | 0 of 35 | 2 of 35 | 33 of 35 | 0.000 | -0.001 to 0.001 | 1.000 | 1.000 | 0 |
| first | K31 plus k-mer counts | `scm` | 4 of 35 | 3 of 35 | 28 of 35 | 0.000 | -0.002 to 0.002 | 0.500 | 1.000 | 0 |
| first | K31 plus k-mer counts | `cart` | 3 of 35 | 3 of 35 | 29 of 35 | 0.000 | -0.002 to 0.000 | 0.656 | 1.000 | 0 |
| first | K31 plus k-mer counts | `rf` | 0 of 35 | 2 of 35 | 33 of 35 | 0.001 | -0.001 to 0.002 | 1.000 | 1.000 | 0 |
| first | K31 plus k-mer counts | `gbt` | 5 of 35 | 3 of 35 | 27 of 35 | 0.001 | 0.001 to 0.004 | 0.363 | 1.000 | 0 |
| first | K31 plus k-mer counts | `lasso` | 2 of 35 | 2 of 35 | 31 of 35 | 0.000 | -0.001 to 0.001 | 0.688 | 1.000 | 0 |
| first | K31 plus amino-acid k-mers | `scm` | 10 of 35 | 3 of 35 | 22 of 35 | 0.004 | -0.001 to 0.014 | 0.046 | 1.000 | 0 |
| first | K31 plus amino-acid k-mers | `cart` | 7 of 35 | 3 of 35 | 25 of 35 | 0.003 | -0.003 to 0.010 | 0.172 | 1.000 | 0 |
| first | K31 plus amino-acid k-mers | `rf` | 6 of 35 | 2 of 35 | 27 of 35 | 0.004 | 0.000 to 0.006 | 0.145 | 1.000 | 0 |
| first | K31 plus amino-acid k-mers | `gbt` | 6 of 35 | 4 of 35 | 25 of 35 | 0.002 | -0.005 to 0.006 | 0.377 | 1.000 | 0 |
| first | K31 plus amino-acid k-mers | `lasso` | 4 of 35 | 4 of 35 | 27 of 35 | 0.000 | -0.002 to 0.003 | 0.637 | 1.000 | 0 |
| first | the AMRFinderPlus program's calls alone | `scm` | 9 of 33 | 9 of 33 | 15 of 33 | 0.003 | -0.019 to 0.014 | 0.593 | 1.000 | 0 |
| first | the AMRFinderPlus program's calls alone | `cart` | 16 of 33 | 9 of 33 | 8 of 33 | 0.019 | -0.005 to 0.038 | 0.115 | 1.000 | 0 |
| first | the AMRFinderPlus program's calls alone | `rf` | 24 of 33 | 4 of 33 | 5 of 33 | 0.064 | 0.044 to 0.094 | < 0.001 | 0.003 | 0 |
| first | the AMRFinderPlus program's calls alone | `gbt` | 22 of 33 | 2 of 33 | 9 of 33 | 0.047 | 0.021 to 0.056 | < 0.001 | < 0.001 | 0 |
| first | the AMRFinderPlus program's calls alone | `lasso` | 18 of 33 | 4 of 33 | 11 of 33 | 0.026 | 0.003 to 0.052 | 0.002 | 0.061 | 0 |
| first | K31 plus the AMRFinderPlus program's calls | `scm` | 1 of 33 | 0 of 33 | 32 of 33 | 0.000 | 0.000 to 0.000 | 0.500 | 1.000 | 0 |
| first | K31 plus the AMRFinderPlus program's calls | `cart` | 1 of 33 | 1 of 33 | 31 of 33 | 0.000 | -0.001 to 0.000 | 0.750 | 1.000 | 0 |
| first | K31 plus the AMRFinderPlus program's calls | `rf` | 1 of 33 | 0 of 33 | 32 of 33 | 0.000 | 0.000 to 0.002 | 0.500 | 1.000 | 0 |
| first | K31 plus the AMRFinderPlus program's calls | `gbt` | 3 of 33 | 1 of 33 | 29 of 33 | 0.000 | -0.001 to 0.003 | 0.312 | 1.000 | 0 |
| first | K31 plus the AMRFinderPlus program's calls | `lasso` | 0 of 33 | 1 of 33 | 32 of 33 | 0.000 | 0.000 to 0.000 | 1.000 | 1.000 | 0 |
| first | the twin, the default at seeds 45, 46 and 47 | `scm` | 4 of 35 | 4 of 35 | 27 of 35 | 0.000 | 0.000 to 0.006 | 0.637 | 1.000 | 0 |
| first | the twin, the default at seeds 45, 46 and 47 | `cart` | 6 of 35 | 4 of 35 | 25 of 35 | 0.000 | -0.001 to 0.014 | 0.377 | 1.000 | 0 |
| first | the twin, the default at seeds 45, 46 and 47 | `rf` | 2 of 35 | 2 of 35 | 31 of 35 | -0.002 | -0.003 to 0.000 | 0.688 | 1.000 | 0 |
| first | the twin, the default at seeds 45, 46 and 47 | `gbt` | 8 of 35 | 4 of 35 | 23 of 35 | 0.000 | -0.001 to 0.007 | 0.194 | 1.000 | 0 |
| first | the twin, the default at seeds 45, 46 and 47 | `lasso` | 4 of 35 | 2 of 35 | 29 of 35 | 0.000 | -0.002 to 0.001 | 0.344 | 1.000 | 0 |
| second | unitigs | `scm` | 2 of 12 | 1 of 12 | 9 of 12 | 0.000 | -0.030 to 0.015 | 0.500 | 1.000 | 0 |
| second | unitigs | `cart` | 0 of 12 | 0 of 12 | 12 of 12 | -0.001 | -0.008 to 0.005 | 1.000 | 1.000 | 0 |
| second | unitigs | `rf` | 1 of 12 | 1 of 12 | 10 of 12 | 0.000 | -0.017 to 0.002 | 0.750 | 1.000 | 0 |
| second | unitigs | `gbt` | 0 of 12 | 0 of 12 | 12 of 12 | 0.000 | -0.003 to 0.014 | 1.000 | 1.000 | 0 |
| second | unitigs | `lasso` | 0 of 12 | 1 of 12 | 11 of 12 | 0.000 | -0.005 to 0.003 | 1.000 | 1.000 | 0 |
| second | K31 plus k-mer counts | `scm` | 2 of 12 | 2 of 12 | 8 of 12 | 0.003 | -0.004 to 0.042 | 0.688 | 1.000 | 0 |
| second | K31 plus k-mer counts | `cart` | 1 of 12 | 3 of 12 | 8 of 12 | -0.003 | -0.015 to 0.001 | 0.938 | 1.000 | 0 |
| second | K31 plus k-mer counts | `rf` | 0 of 12 | 0 of 12 | 12 of 12 | 0.000 | -0.009 to 0.008 | 1.000 | 1.000 | 0 |
| second | K31 plus k-mer counts | `gbt` | 3 of 12 | 0 of 12 | 9 of 12 | 0.008 | 0.002 to 0.035 | 0.125 | 1.000 | 0 |
| second | K31 plus k-mer counts | `lasso` | 1 of 12 | 0 of 12 | 11 of 12 | 0.000 | -0.001 to 0.024 | 0.500 | 1.000 | 0 |
| second | K31 plus amino-acid k-mers | `scm` | 3 of 12 | 1 of 12 | 8 of 12 | 0.012 | -0.006 to 0.016 | 0.312 | 1.000 | 0 |
| second | K31 plus amino-acid k-mers | `cart` | 2 of 12 | 1 of 12 | 9 of 12 | 0.005 | -0.005 to 0.014 | 0.500 | 1.000 | 0 |
| second | K31 plus amino-acid k-mers | `rf` | 2 of 12 | 1 of 12 | 9 of 12 | 0.006 | -0.020 to 0.011 | 0.500 | 1.000 | 0 |
| second | K31 plus amino-acid k-mers | `gbt` | 2 of 12 | 1 of 12 | 9 of 12 | 0.001 | -0.005 to 0.006 | 0.500 | 1.000 | 0 |
| second | K31 plus amino-acid k-mers | `lasso` | 2 of 12 | 1 of 12 | 9 of 12 | 0.003 | -0.004 to 0.006 | 0.500 | 1.000 | 0 |
| second | the AMRFinderPlus program's calls alone | `scm` | 5 of 11 | 2 of 11 | 4 of 11 | 0.014 | -0.054 to 0.177 | 0.227 | 1.000 | 0 |
| second | the AMRFinderPlus program's calls alone | `cart` | 6 of 11 | 1 of 11 | 4 of 11 | 0.033 | -0.009 to 0.148 | 0.062 | 1.000 | 0 |
| second | the AMRFinderPlus program's calls alone | `rf` | 8 of 11 | 1 of 11 | 2 of 11 | 0.079 | -0.045 to 0.095 | 0.020 | 0.586 | 0 |
| second | the AMRFinderPlus program's calls alone | `gbt` | 8 of 11 | 1 of 11 | 2 of 11 | 0.064 | -0.039 to 0.175 | 0.020 | 0.586 | 0 |
| second | the AMRFinderPlus program's calls alone | `lasso` | 5 of 11 | 1 of 11 | 5 of 11 | 0.006 | -0.007 to 0.042 | 0.109 | 1.000 | 0 |
| second | K31 plus the AMRFinderPlus program's calls | `scm` | 0 of 11 | 0 of 11 | 11 of 11 | 0.000 | -0.009 to 0.014 | 1.000 | 1.000 | 0 |
| second | K31 plus the AMRFinderPlus program's calls | `cart` | 0 of 11 | 0 of 11 | 11 of 11 | 0.000 | -0.011 to 0.004 | 1.000 | 1.000 | 0 |
| second | K31 plus the AMRFinderPlus program's calls | `rf` | 0 of 11 | 1 of 11 | 10 of 11 | 0.000 | -0.026 to 0.003 | 1.000 | 1.000 | 0 |
| second | K31 plus the AMRFinderPlus program's calls | `gbt` | 1 of 11 | 2 of 11 | 8 of 11 | 0.000 | -0.011 to 0.062 | 0.875 | 1.000 | 0 |
| second | K31 plus the AMRFinderPlus program's calls | `lasso` | 0 of 11 | 0 of 11 | 11 of 11 | 0.000 | 0.000 to 0.002 | 1.000 | 1.000 | 0 |
| second | the twin, the default at seeds 45, 46 and 47 | `scm` | 1 of 12 | 3 of 12 | 8 of 12 | 0.001 | -0.032 to 0.015 | 0.938 | 1.000 | 0 |
| second | the twin, the default at seeds 45, 46 and 47 | `cart` | 1 of 12 | 2 of 12 | 9 of 12 | -0.005 | -0.017 to 0.018 | 0.875 | 1.000 | 0 |
| second | the twin, the default at seeds 45, 46 and 47 | `rf` | 0 of 12 | 2 of 12 | 10 of 12 | -0.004 | -0.012 to 0.000 | 1.000 | 1.000 | 0 |
| second | the twin, the default at seeds 45, 46 and 47 | `gbt` | 2 of 12 | 2 of 12 | 8 of 12 | -0.004 | -0.017 to 0.037 | 0.688 | 1.000 | 0 |
| second | the twin, the default at seeds 45, 46 and 47 | `lasso` | 0 of 12 | 1 of 12 | 11 of 12 | -0.001 | -0.002 to 0.002 | 1.000 | 1.000 | 0 |

**M alone, recall.** The datasets of M, on the lineage split and recall.

| Reading | Arm | Learner | Better | Worse | Same | Median | Interval |
|---|---|---|---|---|---|---|---|
| first | unitigs | `scm` | 5 of 35 | 6 of 35 | 24 of 35 | -0.001 | -0.005 to 0.000 |
| first | unitigs | `cart` | 4 of 35 | 9 of 35 | 22 of 35 | -0.001 | -0.014 to 0.000 |
| first | unitigs | `rf` | 1 of 35 | 4 of 35 | 30 of 35 | 0.000 | -0.002 to 0.002 |
| first | unitigs | `gbt` | 4 of 35 | 6 of 35 | 25 of 35 | 0.000 | -0.005 to 0.002 |
| first | unitigs | `lasso` | 0 of 35 | 5 of 35 | 30 of 35 | 0.000 | -0.001 to 0.003 |
| first | K31 plus k-mer counts | `scm` | 4 of 35 | 4 of 35 | 27 of 35 | 0.000 | 0.000 to 0.007 |
| first | K31 plus k-mer counts | `cart` | 5 of 35 | 8 of 35 | 22 of 35 | 0.000 | -0.006 to 0.000 |
| first | K31 plus k-mer counts | `rf` | 4 of 35 | 4 of 35 | 27 of 35 | 0.000 | 0.000 to 0.003 |
| first | K31 plus k-mer counts | `gbt` | 9 of 35 | 6 of 35 | 20 of 35 | 0.000 | -0.001 to 0.007 |
| first | K31 plus k-mer counts | `lasso` | 4 of 35 | 3 of 35 | 28 of 35 | 0.000 | -0.001 to 0.006 |
| first | K31 plus amino-acid k-mers | `scm` | 9 of 35 | 6 of 35 | 20 of 35 | 0.002 | -0.001 to 0.014 |
| first | K31 plus amino-acid k-mers | `cart` | 11 of 35 | 6 of 35 | 18 of 35 | 0.000 | -0.005 to 0.011 |
| first | K31 plus amino-acid k-mers | `rf` | 5 of 35 | 4 of 35 | 26 of 35 | 0.002 | 0.000 to 0.007 |
| first | K31 plus amino-acid k-mers | `gbt` | 6 of 35 | 9 of 35 | 20 of 35 | 0.000 | -0.006 to 0.003 |
| first | K31 plus amino-acid k-mers | `lasso` | 5 of 35 | 7 of 35 | 23 of 35 | 0.000 | -0.008 to 0.004 |
| first | the AMRFinderPlus program's calls alone | `scm` | 10 of 33 | 13 of 33 | 10 of 33 | -0.001 | -0.035 to 0.001 |
| first | the AMRFinderPlus program's calls alone | `cart` | 19 of 33 | 8 of 33 | 6 of 33 | 0.027 | -0.019 to 0.067 |
| first | the AMRFinderPlus program's calls alone | `rf` | 19 of 33 | 4 of 33 | 10 of 33 | 0.071 | 0.013 to 0.187 |
| first | the AMRFinderPlus program's calls alone | `gbt` | 20 of 33 | 6 of 33 | 7 of 33 | 0.051 | 0.005 to 0.120 |
| first | the AMRFinderPlus program's calls alone | `lasso` | 19 of 33 | 5 of 33 | 9 of 33 | 0.050 | -0.003 to 0.097 |
| first | K31 plus the AMRFinderPlus program's calls | `scm` | 0 of 33 | 1 of 33 | 32 of 33 | 0.000 | 0.000 to 0.000 |
| first | K31 plus the AMRFinderPlus program's calls | `cart` | 2 of 33 | 2 of 33 | 29 of 33 | 0.000 | 0.000 to 0.000 |
| first | K31 plus the AMRFinderPlus program's calls | `rf` | 2 of 33 | 3 of 33 | 28 of 33 | 0.001 | -0.002 to 0.001 |
| first | K31 plus the AMRFinderPlus program's calls | `gbt` | 3 of 33 | 2 of 33 | 28 of 33 | 0.000 | -0.003 to 0.000 |
| first | K31 plus the AMRFinderPlus program's calls | `lasso` | 0 of 33 | 1 of 33 | 32 of 33 | 0.000 | 0.000 to 0.000 |
| first | the twin, the default at seeds 45, 46 and 47 | `scm` | 6 of 35 | 6 of 35 | 23 of 35 | 0.000 | -0.004 to 0.013 |
| first | the twin, the default at seeds 45, 46 and 47 | `cart` | 10 of 35 | 7 of 35 | 18 of 35 | 0.000 | -0.005 to 0.033 |
| first | the twin, the default at seeds 45, 46 and 47 | `rf` | 1 of 35 | 3 of 35 | 31 of 35 | 0.000 | 0.000 to 0.002 |
| first | the twin, the default at seeds 45, 46 and 47 | `gbt` | 4 of 35 | 6 of 35 | 25 of 35 | 0.001 | -0.001 to 0.004 |
| first | the twin, the default at seeds 45, 46 and 47 | `lasso` | 5 of 35 | 5 of 35 | 25 of 35 | 0.000 | -0.003 to 0.001 |
| second | unitigs | `scm` | 1 of 12 | 2 of 12 | 9 of 12 | 0.000 | -0.068 to 0.002 |
| second | unitigs | `cart` | 0 of 12 | 1 of 12 | 11 of 12 | 0.001 | -0.007 to 0.009 |
| second | unitigs | `rf` | 0 of 12 | 1 of 12 | 11 of 12 | 0.001 | -0.035 to 0.007 |
| second | unitigs | `gbt` | 1 of 12 | 3 of 12 | 8 of 12 | 0.000 | -0.033 to 0.018 |
| second | unitigs | `lasso` | 0 of 12 | 1 of 12 | 11 of 12 | 0.000 | -0.001 to 0.011 |
| second | K31 plus k-mer counts | `scm` | 2 of 12 | 1 of 12 | 9 of 12 | 0.007 | 0.000 to 0.039 |
| second | K31 plus k-mer counts | `cart` | 1 of 12 | 5 of 12 | 6 of 12 | -0.009 | -0.029 to 0.003 |
| second | K31 plus k-mer counts | `rf` | 1 of 12 | 1 of 12 | 10 of 12 | 0.000 | -0.013 to 0.011 |
| second | K31 plus k-mer counts | `gbt` | 4 of 12 | 2 of 12 | 6 of 12 | 0.001 | -0.012 to 0.032 |
| second | K31 plus k-mer counts | `lasso` | 1 of 12 | 0 of 12 | 11 of 12 | -0.001 | -0.006 to 0.039 |
| second | K31 plus amino-acid k-mers | `scm` | 4 of 12 | 3 of 12 | 5 of 12 | 0.001 | -0.023 to 0.012 |
| second | K31 plus amino-acid k-mers | `cart` | 5 of 12 | 2 of 12 | 5 of 12 | 0.000 | -0.055 to 0.041 |
| second | K31 plus amino-acid k-mers | `rf` | 2 of 12 | 2 of 12 | 8 of 12 | 0.005 | -0.039 to 0.011 |
| second | K31 plus amino-acid k-mers | `gbt` | 1 of 12 | 3 of 12 | 8 of 12 | 0.000 | -0.016 to 0.001 |
| second | K31 plus amino-acid k-mers | `lasso` | 2 of 12 | 1 of 12 | 9 of 12 | -0.003 | -0.011 to 0.042 |
| second | the AMRFinderPlus program's calls alone | `scm` | 6 of 11 | 5 of 11 | 0 of 11 | 0.021 | -0.157 to 0.606 |
| second | the AMRFinderPlus program's calls alone | `cart` | 7 of 11 | 2 of 11 | 2 of 11 | 0.027 | -0.066 to 0.145 |
| second | the AMRFinderPlus program's calls alone | `rf` | 6 of 11 | 2 of 11 | 3 of 11 | 0.022 | -0.090 to 0.291 |
| second | the AMRFinderPlus program's calls alone | `gbt` | 7 of 11 | 1 of 11 | 3 of 11 | 0.031 | -0.112 to 0.479 |
| second | the AMRFinderPlus program's calls alone | `lasso` | 6 of 11 | 2 of 11 | 3 of 11 | 0.022 | -0.020 to 0.206 |
| second | K31 plus the AMRFinderPlus program's calls | `scm` | 1 of 11 | 0 of 11 | 10 of 11 | 0.000 | 0.000 to 0.020 |
| second | K31 plus the AMRFinderPlus program's calls | `cart` | 0 of 11 | 0 of 11 | 11 of 11 | 0.000 | 0.000 to 0.004 |
| second | K31 plus the AMRFinderPlus program's calls | `rf` | 0 of 11 | 2 of 11 | 9 of 11 | 0.000 | -0.053 to 0.000 |
| second | K31 plus the AMRFinderPlus program's calls | `gbt` | 2 of 11 | 1 of 11 | 8 of 11 | -0.001 | -0.019 to 0.133 |
| second | K31 plus the AMRFinderPlus program's calls | `lasso` | 0 of 11 | 0 of 11 | 11 of 11 | 0.000 | 0.000 to 0.003 |
| second | the twin, the default at seeds 45, 46 and 47 | `scm` | 2 of 12 | 4 of 12 | 6 of 12 | -0.001 | -0.067 to 0.021 |
| second | the twin, the default at seeds 45, 46 and 47 | `cart` | 0 of 12 | 3 of 12 | 9 of 12 | -0.003 | -0.023 to -0.001 |
| second | the twin, the default at seeds 45, 46 and 47 | `rf` | 0 of 12 | 1 of 12 | 11 of 12 | -0.001 | -0.019 to 0.001 |
| second | the twin, the default at seeds 45, 46 and 47 | `gbt` | 1 of 12 | 4 of 12 | 7 of 12 | -0.001 | -0.044 to 0.081 |
| second | the twin, the default at seeds 45, 46 and 47 | `lasso` | 1 of 12 | 3 of 12 | 8 of 12 | 0.000 | -0.004 to 0.005 |

**M alone, study split.** The datasets of M, on the study split and balanced accuracy.

| Reading | Arm | Learner | Better | Worse | Same | Median | Interval |
|---|---|---|---|---|---|---|---|
| first | unitigs | `scm` | 2 of 14 | 3 of 14 | 9 of 14 | 0.003 | 0.003 to 0.014 |
| first | unitigs | `cart` | 2 of 14 | 0 of 14 | 12 of 14 | 0.000 | -0.002 to 0.010 |
| first | unitigs | `rf` | 0 of 14 | 1 of 14 | 13 of 14 | -0.002 | -0.009 to 0.004 |
| first | unitigs | `gbt` | 2 of 14 | 0 of 14 | 12 of 14 | 0.003 | -0.002 to 0.006 |
| first | unitigs | `lasso` | 1 of 14 | 0 of 14 | 13 of 14 | 0.001 | -0.002 to 0.008 |
| first | K31 plus k-mer counts | `scm` | 1 of 14 | 3 of 14 | 10 of 14 | 0.000 | -0.015 to 0.019 |
| first | K31 plus k-mer counts | `cart` | 2 of 14 | 4 of 14 | 8 of 14 | 0.000 | -0.033 to 0.013 |
| first | K31 plus k-mer counts | `rf` | 0 of 14 | 4 of 14 | 10 of 14 | -0.008 | -0.016 to 0.008 |
| first | K31 plus k-mer counts | `gbt` | 0 of 14 | 1 of 14 | 13 of 14 | 0.000 | -0.001 to 0.004 |
| first | K31 plus k-mer counts | `lasso` | 0 of 14 | 1 of 14 | 13 of 14 | 0.002 | -0.002 to 0.008 |
| first | K31 plus amino-acid k-mers | `scm` | 7 of 14 | 2 of 14 | 5 of 14 | 0.019 | 0.008 to 0.040 |
| first | K31 plus amino-acid k-mers | `cart` | 5 of 14 | 3 of 14 | 6 of 14 | 0.003 | -0.017 to 0.050 |
| first | K31 plus amino-acid k-mers | `rf` | 3 of 14 | 1 of 14 | 10 of 14 | 0.007 | -0.002 to 0.008 |
| first | K31 plus amino-acid k-mers | `gbt` | 6 of 14 | 2 of 14 | 6 of 14 | 0.008 | 0.003 to 0.035 |
| first | K31 plus amino-acid k-mers | `lasso` | 4 of 14 | 1 of 14 | 9 of 14 | -0.002 | -0.010 to 0.035 |
| first | the AMRFinderPlus program's calls alone | `scm` | 5 of 14 | 3 of 14 | 6 of 14 | 0.002 | -0.009 to 0.004 |
| first | the AMRFinderPlus program's calls alone | `cart` | 6 of 14 | 0 of 14 | 8 of 14 | 0.013 | 0.001 to 0.032 |
| first | the AMRFinderPlus program's calls alone | `rf` | 9 of 14 | 2 of 14 | 3 of 14 | 0.065 | 0.012 to 0.070 |
| first | the AMRFinderPlus program's calls alone | `gbt` | 8 of 14 | 2 of 14 | 4 of 14 | 0.037 | -0.018 to 0.053 |
| first | the AMRFinderPlus program's calls alone | `lasso` | 8 of 14 | 0 of 14 | 6 of 14 | 0.044 | 0.024 to 0.057 |
| first | K31 plus the AMRFinderPlus program's calls | `scm` | 0 of 14 | 1 of 14 | 13 of 14 | 0.000 | 0.000 to 0.000 |
| first | K31 plus the AMRFinderPlus program's calls | `cart` | 0 of 14 | 0 of 14 | 14 of 14 | 0.000 | -0.001 to 0.001 |
| first | K31 plus the AMRFinderPlus program's calls | `rf` | 0 of 14 | 2 of 14 | 12 of 14 | -0.004 | -0.006 to 0.004 |
| first | K31 plus the AMRFinderPlus program's calls | `gbt` | 0 of 14 | 0 of 14 | 14 of 14 | 0.000 | -0.007 to 0.001 |
| first | K31 plus the AMRFinderPlus program's calls | `lasso` | 0 of 14 | 0 of 14 | 14 of 14 | 0.000 | 0.000 to 0.000 |
| second | unitigs | `scm` | 1 of 3 | 0 of 3 | 2 of 3 | 0.013 | 0.001 to 0.082 |
| second | unitigs | `cart` | 0 of 3 | 0 of 3 | 3 of 3 | 0.000 | 0.000 to 0.001 |
| second | unitigs | `rf` | 0 of 3 | 0 of 3 | 3 of 3 | 0.002 | 0.002 to 0.004 |
| second | unitigs | `gbt` | 1 of 3 | 0 of 3 | 2 of 3 | 0.003 | 0.003 to 0.015 |
| second | unitigs | `lasso` | 0 of 3 | 0 of 3 | 3 of 3 | 0.005 | 0.003 to 0.007 |
| second | K31 plus k-mer counts | `scm` | 1 of 3 | 0 of 3 | 2 of 3 | 0.016 | 0.008 to 0.067 |
| second | K31 plus k-mer counts | `cart` | 0 of 3 | 1 of 3 | 2 of 3 | 0.000 | -0.025 to 0.003 |
| second | K31 plus k-mer counts | `rf` | 0 of 3 | 1 of 3 | 2 of 3 | 0.009 | -0.007 to 0.014 |
| second | K31 plus k-mer counts | `gbt` | 0 of 3 | 0 of 3 | 3 of 3 | 0.011 | 0.008 to 0.016 |
| second | K31 plus k-mer counts | `lasso` | 0 of 3 | 1 of 3 | 2 of 3 | 0.001 | -0.011 to 0.009 |
| second | K31 plus amino-acid k-mers | `scm` | 2 of 3 | 0 of 3 | 1 of 3 | 0.024 | 0.010 to 0.048 |
| second | K31 plus amino-acid k-mers | `cart` | 1 of 3 | 1 of 3 | 1 of 3 | 0.001 | -0.054 to 0.013 |
| second | K31 plus amino-acid k-mers | `rf` | 0 of 3 | 1 of 3 | 2 of 3 | -0.004 | -0.019 to -0.004 |
| second | K31 plus amino-acid k-mers | `gbt` | 2 of 3 | 0 of 3 | 1 of 3 | 0.021 | 0.021 to 0.032 |
| second | K31 plus amino-acid k-mers | `lasso` | 1 of 3 | 1 of 3 | 1 of 3 | 0.012 | -0.020 to 0.049 |
| second | the AMRFinderPlus program's calls alone | `scm` | 2 of 3 | 1 of 3 | 0 of 3 | 0.025 | 0.025 to 0.050 |
| second | the AMRFinderPlus program's calls alone | `cart` | 1 of 3 | 0 of 3 | 2 of 3 | 0.008 | 0.000 to 0.071 |
| second | the AMRFinderPlus program's calls alone | `rf` | 1 of 3 | 0 of 3 | 2 of 3 | 0.019 | -0.015 to 0.079 |
| second | the AMRFinderPlus program's calls alone | `gbt` | 1 of 3 | 1 of 3 | 1 of 3 | -0.014 | -0.042 to 0.068 |
| second | the AMRFinderPlus program's calls alone | `lasso` | 2 of 3 | 0 of 3 | 1 of 3 | 0.058 | 0.032 to 0.058 |
| second | K31 plus the AMRFinderPlus program's calls | `scm` | 0 of 3 | 0 of 3 | 3 of 3 | 0.000 | 0.000 to 0.000 |
| second | K31 plus the AMRFinderPlus program's calls | `cart` | 0 of 3 | 0 of 3 | 3 of 3 | -0.001 | -0.006 to 0.000 |
| second | K31 plus the AMRFinderPlus program's calls | `rf` | 0 of 3 | 0 of 3 | 3 of 3 | 0.000 | 0.000 to 0.002 |
| second | K31 plus the AMRFinderPlus program's calls | `gbt` | 0 of 3 | 0 of 3 | 3 of 3 | -0.006 | -0.012 to -0.003 |
| second | K31 plus the AMRFinderPlus program's calls | `lasso` | 0 of 3 | 0 of 3 | 3 of 3 | 0.000 | 0.000 to 0.000 |

**Files.** The 8 tables above, 370 rows for the 9 species of M, come from the record, `reference/input-arms-m/`, which holds every other table: each count on every measure and split (`m-count.tsv`), each count's conditions (`m-conditions.tsv`), the size table (`m-size.tsv`) and each dataset's figures (`m-datasets.tsv`). `scripts/c3-inputs/m_page.py print` prints them.
