GRM revision 14 | Matched prediction accuracy
Amendment 1 to Protocol 1.0 | 26 November 2026
Release rev14-274-20261112
Protocol SHA-256 e4c124466f7b85c649ea4684004831beb80002174e51317ad0ee7476a84401b0
Ruling R49 SHA-256 77a2d5edf0e25a7a9b89b4065e87dcf6daa6d20940eb601a21a2f45fc0b3b232

Split

Use only the five held-out folds of the lineage split. Each lineage cluster sits wholly in one fold. All 274 reports contain both the lineage and random splits.

Reports and learners

Include the 268 reports with scm, cart, rf, gbt and lasso under the lineage split. Count and exclude the six whose lineage split contains none of the five. Identify learners by name. The catalogue entry is outside the five-learner mean. The four catalogue-call cells are aggregate counts.

Sources

Take the frozen lab call, catalogue call, lineage cluster, MLST type and held-out lineage fold from the report. Stage 2 alone uses learner scores, if authorised. Take the unique BioProject from GRM's study-split mapping of the cut BV-BRC tables.

In the cut bvbrc/genome_amr.jsonl, read the lab method (`laboratory_typing_method`), standard (`testing_standard`) and version (`testing_standard_year`) together on each applicable record for the genome-antibiotic pair. Require one distinct tuple with a documented method and standard. Repeated identical tuples qualify. A blank year is the unknown-version level. Exclude missing or ambiguous required assignments.

Use the HC20 field (`cgmlst_hc20`) from the cut bvbrc/genomes.jsonl for the MLST plus HC20 comparison. Count missing coverage there. Zero eligible pairs means zero support and an unestimated concordance.

Intermediate records under R2

Intermediate records affect 35 of the 268 datasets in nine species. Reading A uses eligible release-labelled rows. Reading B removes every genome-antibiotic pair with any intermediate lab record in the cut table. Reapply group eligibility after removal. Freeze both row lists in Stage 1, keep them local and report each SHA-256.

Report both readings for every figure that includes an affected dataset. An above-chance primary claim requires both lower interval bounds to exceed 0.5. An interval spanning 0.5 is unresolved. Interpret results retrospectively within this database.

Execution under R49

A second worker reviews the script, tests, planted faults, protocol and amendment before public-repository merge. A GRM worker runs the merged Stage 1 on a GRM machine after the 35-dataset run starts. The machine holds neither sealed records nor the 21 November download. Verify the cut record of 1,285,111 kept before 2026-11-05T19:10:00 and 252,148 dropped. Discard every score list on report load and prove this in a test.

Export counts by dataset, species and overall, with the two hashes. Keep genome, pair, group and BioProject rows and all accessions local. Run the genome-ID scan with its positive control before export. Pre-freeze and non-GRM runs remain exploratory under R1.

Stage 2 awaits a separate ruling on the Stage 1 counts. If authorised, use one execution of the merged script, verify both frozen row-list hashes and report both readings whatever they show. No cell names or ranks a gene, feature or group.
