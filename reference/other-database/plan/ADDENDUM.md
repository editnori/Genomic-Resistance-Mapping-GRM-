# Lane outside2: addendum to the plan for the eight other species

Written 15 Nov 2026 before any susceptibility result of the eight species below was read. Its sha256 is recorded in `ADDENDUM.sha256` with the time, before the first label of these species is fetched. It adds to `outside/PLAN.md` (sha256 `0dc4583857534d539b8c357b29f0c826f1683d6680dcc082a1f69e933d4f542d`) and changes no rule of it: the funnel steps, the draw of 4,000 per species (P2, the same seed text `20261115:`), the labels, the measures, the readings, the pass rule per learner and species, and the words all hold as written. It fixes only what that plan left open for further species. Anything decided after a label of these species is read is labelled exploratory.

## A1. Species and order

The eight release species that meet the plan's size rule S8 in the first lane's survey (`outside/results/survey-by-species.tsv`, run of 12:28:58Z), besides the pilot species *Escherichia coli*, run in falling `candidate_genome_drug_pairs`:

1. *Salmonella enterica* (65,610)
2. *Klebsiella pneumoniae* (15,787)
3. *Staphylococcus aureus* (11,680)
4. *Campylobacter jejuni* (9,771)
5. *Campylobacter coli* (8,564)
6. *Neisseria gonorrhoeae* (4,238)
7. *Pseudomonas aeruginosa* (3,980)
8. *Acinetobacter baumannii* (3,866)

The pilot species is not run again. Its row in every table is copied from the first lane (`outside/results/pilot/`) and marked so.

## A2. One new funnel step, counted

- O1. Outside-study candidate: a candidate whose BioSample, assembly key or any SRA run matches a candidate of the outside-study test of the same species (`scripts/outside-study/frozen/candidates.tsv.gz` of the repository at `852ebc9`, its genome_id, species_slug, biosample, sra_accession and assembly_accession fields only) is removed, so the two tests share no genome. The outside-study test is a retrospective test on newly deposited studies from the same database, and this test is kept apart from it.
- O1 comes after D2 and before D4, so a within-batch group keeps a genome that is not an outside-study candidate when it has one. Its count is written in each species' `funnel-label-blind.tsv`. The survey shows 441 such isolates for *Neisseria gonorrhoeae* and none for the other seven.

## A3. Each species as the plan runs it

- Training genomes, the replay check, QC, the species check, the sequence duplicates, the reports and the measures are those of the plan, taken for each species in turn: the training genomes are the union of `samples[].id` over that species' release reports; Q1's median is that species' median over the data root's BV-BRC table; Q2 runs on a scratch root of its own per species, holding the release's recorded check (`species-check/species_check.json` of the data root, unchanged), the 23 recorded references, that species' training genomes and its candidates.
- The reports of a species are pinned by sha256 (`results/<species>/reports.sha256`) before its first outside genome is scored.
- Labels last, per species: the label fetch of a species starts only after every prediction row of that species is written.

## A4. How the eight species are read together

- Counts only. For each learner (SCM, CART, lasso): the number of species where it carries over, where it does not carry over, and where it is not assessable, each by the plan's pass rule applied to that species alone.
- No pooled estimate across species, no median across species and no test across species.
- A species that is not run, or that stops (rule A5), is listed with the reason and counted in none of the three groups. A species where `grm predict` scores no report for a learner has that learner not assessable there, with the reason.
- The first lane's pilot result is shown beside these counts, marked as the first lane's, and is not added into them.
- RF and GBT are not scored: they are published as importances only. Every table says so.
- The result carries only the words of the plan's section 7: "a test on data from another public database". It is not external validation and not the outside-study test.

## A5. Stop rules

- No new species starts at or after 20:30Z (clock). A species whose downloads have begun is finished.
- Before each species, `df -h /home` is read. If less than 60 GB is free, no further species starts.
- If a species' replay check fails, that species stops and the rest go on.
- If `grm predict` cannot score a species' reports, the reason is written and the rest go on.
- On HTTP 429, 401 or 503 from NCBI, that step stops, the exact answer and its time are written, and no request is repeated in a loop.

## A6. Operations that change no rule

- At most 3 downloads at a time and one request at a time to each NCBI interface across all of this lane's processes (a file lock per interface), the generic client name of the first lane, no BV-BRC request.
- `grm predict` with at most 5 processes; the whole lane at most 8 threads.
- The lane's scratch stays under 10 GB. A species' FASTA files are deleted once its prediction rows are written and every file's sha256 is saved in `results/<species>/fasta-sha256.tsv.gz`; its sketches stay. Measures can be made again from the saved predictions and labels without a new download or a new prediction.
