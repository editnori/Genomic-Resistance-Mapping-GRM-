# Input configuration: today's values and what reaches them

Layth specified a full-pipeline comparison over the input data: genome and cohort count, which genomes are included, k-mer length, feature families and the retained-feature budget. It isn't built, integrated or run, so no difference has been measured. This page gives each of the five inputs as the code has it today: its value, whether a flag sets it, and what reaches it from `src/main.rs`. [input-comparison-runner.md](input-comparison-runner.md) lists the components a comparison would use and what is missing. Line numbers are at commit `b90be1ca`.

| Input | Today's value | Flag | Evaluated as a factor |
|---|---|---|---|
| Genome and cohort count | at most 2,000 genomes per species; at least 25 per class per dataset | `grm curate --cap`, `--min-class` | no |
| Which genomes are included | `select_genomes`, a label-weighted ranking | none | no |
| k-mer length | 31 | `grm count --k` only; `grm pipeline` recounts at 31 | no |
| Feature families | the K block (k-mer patterns) only | none | no |
| Retained-feature budget | 4,000 screened features per fit, with the support floors below | `grm matrix --min-support` only; `grm pipeline` rebuilds at 2 | no |

Learner settings are a separate thing. `scripts/hyperopt` searches learner knobs only (`scripts/hyperopt/protocol.py:141-196`), and none of these five inputs.

## Genome and cohort count

**Today's value.** `grm curate` selects at most 2,000 genomes per species (`--cap`, `src/main.rs:35-37`). A species × drug dataset exists only when each class has at least 25 genomes with an unconflicted reported call (`--min-class`, `main.rs:38-41`; `QcConfig`, `src/curate.rs:178-199`).

**Flags and what reaches them.** Both flags are on `grm curate` only. `main.rs:175-179` passes them to `run_curate` (`main.rs:508`), which:
- refuses a class minimum below the outer fold count, `EvalConfig::default().folds` (`main.rs:519-523`);
- puts the class minimum in `QcConfig` (`main.rs:526-529`), which `curate::curate` applies to each dataset (`curate.rs:544`) and records in it (`curate.rs:559`);
- passes the cap to `curate::select_genomes`, once per species (`main.rs:611`).

`grm learn` applies each dataset's recorded class minimum again after download (`src/analysis.rs:781`, `:919-920`).

**What is recorded.** The class minimum is written to `curated/qc.json` (`main.rs:592-595`), to each dataset in `curated/datasets.json` (`curate.rs:256-257`) and to each dataset report (`analysis.rs:746`). The cap isn't written anywhere: `genomes.txt` holds the selected ids (`main.rs:623-625`), and `species.json` holds counts (`main.rs:489-498`, `:636-639`).

**In the pipeline.** `grm pipeline` doesn't run `curate` (`main.rs:203-245`). It uses whatever `genomes.txt` the last `grm curate` wrote. A data root (`--data`, `main.rs:19-20`) holds one curated cohort, since `curated/` and each species' folder are single paths (`src/layout.rs:31`, `:45`).

**As a factor.** Nothing compares two caps or two class minimums in one evaluation. The work program's A09, "Move beyond the 2,000-genome cap", is partial (`docs/work-program-20261027.md:136`).

## Which genomes are included

**Today's rule.** `select_genomes` (`curate.rs:623-671`) scores each genome by summing 1 / (size of that call's class) over its reported R or S calls on the species' dataset drugs. It keeps the top `cap`, with ties broken by genome id. A genome whose only labels were inferred from MIC scores nothing and isn't picked. So the choice reads labels, and one drug's cohort depends on the other drugs' labels.

**Flag.** None. `main.rs:611` is the only call to a selection rule.

**The opt-in alternative.** F's `plan_group_sampling` (`curate.rs:798`) draws genomes round-robin across groups without reading a label, except through a declared always-included stratum. Its design is named `group-round-robin-v1` (`curate.rs:675`). It has no caller outside `src/curate.rs`. E's read-only audit (`scripts/sampling-audit/`) ports `select_genomes` to Python and reports what it selects; no `grm` command calls it.

**Recorded direction.** The work program's entry of 30 Oct 10:1x UTC (`docs/work-program-20261027.md` at commit `744e636`, line 596; this branch's copy of the file doesn't have the entry) says: keep a legacy selection that is byte-compatible with today's, and add a separately versioned group-stratified design. "No unreviewed switch." No switch has been approved.

**What a call from `run_curate` would have to supply.** `plan_group_sampling` takes:
- `eligible`, one `EligibleGenome` per genome: its id and its group id, for a representative already past QC and deduplication (`curate.rs:690-701`);
- `always`, the always-included stratum: its ids, a reason, and whether labels chose it (`curate.rs:703-713`);
- `cap` and `seed`.

It returns a `GroupSamplingPlan` (`curate.rs:724-742`), which records the design, the generator versions, a hash of the eligible list, the seed, the cap, the stratum, the selected ids in draw order and each genome's inclusion probability.

`run_curate` has the cap today, but no seed, no stratum and no group ids. Its lineage groups don't exist at curate time: `grm structure` builds them later from the selected genomes' ksets (`src/cluster.rs:1-4`), and `grm learn` reads them from `structure.json` (`analysis.rs:794-798`). The only group source in the tree that works before the FASTA download is E's audit script, which reads one BV-BRC metadata field such as `cgmlst_hc100` (`scripts/sampling-audit/groups_from_metadata.py`). In the revision-9 metadata, whole species have no cgMLST (`scripts/sampling-audit/README.md:102`).

The rules recorded for such a call:
- **A versioned option.** The legacy `select_genomes` stays the default and stays byte-compatible. Nothing names the legacy rule's version today; E's port reproduces its output (`scripts/sampling-audit/legacy_select.py`).
- **One species per call** (`curate.rs:748-753`; [ablation-runner.md](ablation-runner.md#group-sampling-one-species-per-call)). A mixed call is wrong even when the species' group ids are disjoint, because one call has one cap and one stratum for all of them.
- **The cap and the class minimum recorded with the cohort.** A plan records its cap (`curate.rs:735`). Today the class minimum is recorded and the cap isn't (see above).

## k-mer length

**Today's value.** 31, `PIPELINE_K` (`main.rs:259`).

**Flag.** `grm count --k` takes 5 to 32 (`main.rs:67-72`, `src/dna.rs:17-18`), with `PIPELINE_K` as its default (`main.rs:69`). `grm pipeline` has no `--k` (`main.rs:107-121`).

**Where it's fixed.** `grm pipeline` uses `PIPELINE_K` to check the catalog (`main.rs:209`) and to recount unreadable ksets (`:231`). `count_to_matrix` uses it to load the catalog (`:267`), to check which stages are stale (`:273`) and to count (`:284`).

**What refuses another k.** `species::stale_stages` marks counting stale when `counted.json` was counted at another k (`src/species.rs:416`), and `count_to_matrix` passes it `PIPELINE_K`. So `grm pipeline` recounts stages built at another k, at 31, and rebuilds structure and matrix (`README.md:118`; test `stages_built_by_hand_at_another_k_or_min_support_are_rebuilt`, `main.rs:1081`). Ksets on disk are reused only at the same k (`species.rs:102`).

`grm learn` on its own doesn't compare k with 31. `load_inputs` refuses a structure or matrix that wasn't built from `counted.json` as it is (`analysis.rs:279-313`), and then uses the matrix's own k (`analysis.rs:319`, `:534`). Each report records k (`analysis.rs:743`), and the input fingerprint changes with it (`analysis.rs:2180`, in `fingerprint_changes_with_every_input_the_evaluation_reads`).

**What a configurable k would have to reach.** These read k today, or read something built at k:
- the catalog gene calls, loaded in a k-mer space of that k (`main.rs:318-337`);
- the stale-stage check (`species.rs:321`, `:416`) and kset reuse (`species.rs:102`);
- the stage folders, which hold one k per species in a data root (`layout.rs:45-55`);
- structure's MinHash sketches and Mash distances (`cluster.rs:3`, `:43`), and so the lineage labels that the lineage split and G's fold plans read;
- the kset digests behind the k-mer-identical duplicate groups (`analysis.rs:799`, `:814-822`), which the random split keeps whole;
- the report's region probe sets, built at the matrix's k (`analysis.rs:529-539`).

So a change of k can change the lineage labels and the duplicate groups, and with them the folds. The hyperparameter protocol has no k (`protocol.py:141-196`).

**Real cohorts.** The work program records that synthetic k = 21 and 31 checks passed, and that the real-cohort comparison "was explicitly refused and must not be rerouted" (`docs/work-program-20261027.md:148`, `:339`). That refusal stands. k stays fixed at 31 for the product unless Layth decides otherwise.

## Feature families

**Today.** The learners read one family, the K block's k-mer patterns. `eval_inputs` projects the species `FeatureMatrix` (`analysis.rs:793`), and `grm learn` evaluates that projection (`analysis.rs:542`). Catalog gene calls (`analysis.rs:319`) feed only the catalog baseline (`analysis.rs:556`) and the gene panels; no learner reads them.

**Flag.** None.

**H's registry.** `REGISTRY` (`src/ablate/blocks.rs:209`) defines six blocks. 4 of 6 adapter types present; N and X unavailable.

| Block | What it is | Adapter | State |
|---|---|---|---|
| C | curated exact calls (region probes) | `ProbeSource` (`blocks.rs:999`) | present |
| H | hierarchical panel calls | `PanelSource` (`blocks.rs:1166`) | present |
| P | catalog gene presence | `PresenceSource` (`blocks.rs:1287`) | present |
| N | gene copy counts | `CountSource` (`blocks.rs:1414`) | unavailable: nothing in GRM counts gene copies (`blocks.rs:204`) |
| X | genome context | `ContextSource` (`blocks.rs:1446`) | unavailable: no reviewed genome-context source exists (`blocks.rs:205`) |
| K | k-mer patterns | `KmerSource` (`blocks.rs:1482`) | present |

`admit` refuses N and X with those reasons (`blocks.rs:848-850`). No production path admits any block: `ablate::*` has no caller outside `src/ablate/`.

**Conditions for a future N or X source,** from the module doc and `admit`:
- a reviewed source. The registry marks N and X `unavailable` because neither has one (`blocks.rs:195-196`, `:237-254`);
- a sealed adapter: only adapters in `blocks.rs` implement `sealed::Adapter`, and `admit` accepts only those (`blocks.rs:27-36`, `:797-800`, `:843`);
- label-free: `build` sees sample ids, never phenotypes (`blocks.rs:17-18`, `:166-180`);
- one typed genotype input, never a group id (`blocks.rs:31-33`);
- proxy refusal: a block whose id, version or column names match a source proxy, such as a lineage, study, country or collection date, is refused (`blocks.rs:20-25`, `:734`, `:779`);
- no filled rows: a sample the source can't build is left out, never given a default row (`blocks.rs:166-176`);
- a non-empty `version` and a 64-hex `source_sha256` (ablation contract §5; `check_spec`, `blocks.rs:342`).

Building such a source is biological feature engineering, which is outside the current work.

## Retained-feature budget

**Today.** Fixed, with no setting `grm learn` reads: it uses `EvalConfig::default()` (`main.rs:385`).

| Value | What it limits | Set at | Applied at | Flag |
|---|---|---|---|---|
| `screened_features` 4,000 | chi-square-screened features offered to the forest, boosting and the lasso, per fit, on that fit's own fitting samples | `src/learn/eval.rs:77` | forest `eval.rs:371`, boosting `:421`, lasso `:505`, MIC regression `:1304` | none |
| `K_SCREEN_CAP` 4,000 | K block columns kept by `screen_block` | `blocks.rs:202`; K's `screen_cap` at `:259` | `screen_block` (`src/ablate/encode.rs:562`), which has no production caller | none |
| matrix `min_support` 2 | k-mers kept in the species matrix: present in, and absent from, at least 2 genomes | `PIPELINE_MIN_SUPPORT`, `main.rs:262` | `count_to_matrix` (`main.rs:274`, `:307`) | `grm matrix --min-support` (`main.rs:88-91`); `grm pipeline` rebuilds a matrix built at another value (`species.rs:411`) |
| dataset `min_support` max(n / 200, 3) | each fit's candidate features, through `training_candidates` | `analysis.rs:785` | `eval.rs:831` | none |
| lasso `max_nonzero` 400 | the lasso path stops once more than 400 weights are active | literal at `eval.rs:507` | `lasso::fit_path_with_status` (`src/learn/lasso.rs:208-209`, `:228`) | none |

A test pins `K_SCREEN_CAP` to `EvalConfig::default().screened_features` (`blocks.rs:1734-1737`, in `the_registry_lists_six_blocks_in_collapse_order`). `grm learn` refuses a species matrix whose `min_support` is above the dataset's (`analysis.rs:788-792`). `scripts/hyperopt` searches none of these values (`protocol.py:141-196`).
