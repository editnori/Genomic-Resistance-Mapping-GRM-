GRM revision 14 | Matched prediction accuracy
Protocol 1.0 | 26 November 2026

Question

Within the same species, drug, held-out fold, recorded lineage, BioProject and assay, how often does the frozen five-learner mean score a lab-resistant isolate above a lab-susceptible one? Study release rev14-274-20261112 only.

Stage 1 | Keep scores closed

Verify the release's training inputs, labels, fold assignments and score-row order from archived records. Every underlying BV-BRC source must predate the seal on 5 November 2026 at 19:10:00, using the archive's insertion timestamp and time convention. Stop if this cannot be verified from the existing archive or if eligible samples, including known aliases, cannot be verified as held out from their scorer's training set.

Use the release's binary lab calls. Require one identifiable biological sample, fold, recorded lineage, BioProject, method and standard per row. Exclude ambiguous or missing required fields and all score rows for samples duplicated within a dataset. Group by species-drug dataset, fold, lineage, BioProject, method, standard and recorded standard version. Unknown version is a separate level. Keep groups containing both lab outcomes.

Freeze the eligible rows before reading scores. Report exclusions, contributing samples, rows, groups, datasets and species, and the largest project's share of eligible rows and of total analysis weight. Stop if no group qualifies.

Stage 2 | Primary measure

Average the five archived held-out scores in release order, retaining recorded null fits. Stop on a missing or nonfinite component. Compare every resistant-susceptible pair within its group. Count one for a higher resistant score, half for an exact tie and zero otherwise. Divide by the group's pair count. Average groups equally within dataset, datasets equally within species and species equally overall.

With six or more species, resample species means with replacement 10,000 times at the observed species count using NumPy PCG64 seed 42. Use the linearly interpolated 2.5th and 97.5th percentiles. With fewer species, report the estimate and counts alone. A lower bound above 0.5 supports within-group discrimination on this cohort. An interval spanning 0.5 leaves the claim unresolved.

Secondary measures | Descriptive

Report all four combinations of GRM's frozen catalogue call for the resistant and susceptible member, in that order. Count unevaluable calls and absent drug rules separately. Give eligible groups, pairs, datasets and species per cell. Leave empty cells unestimated.

Repeat the primary comparison matching exact MLST, then both MLST and the recorded cgMLST group at HC20. Exclude missing typing values in those reads. At each step recompute the broader comparison on the same surviving isolates and datasets. Report lost coverage. Separately repeat the primary comparison excluding folds with any null fit. Keep the primary weighting throughout.

Execution

The 13 November limit applies. Return aggregate counts and score summaries only. Each stage requires the GRM lead's express authorization. Record any amendment and whether scores had been opened.
