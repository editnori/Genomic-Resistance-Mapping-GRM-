use anyhow::{Context, Result};
use clap::builder::RangedU64ValueParser;
use clap::{Parser, Subcommand};
use grm::bvbrc::{self, AmrRecord, Client, GenomeRecord};
use grm::catalog::Catalog;
use grm::curate::{self, QcConfig};
use grm::dna::{MAX_K, MIN_K};
use grm::drugs;
use grm::layout::{Layout, ensure_parent};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::Write;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "grm", about = "Genomic Resistance Mapping")]
struct Cli {
    /// Data root for downloads, intermediate files, and results.
    #[arg(long, global = true, default_value = "data")]
    data: PathBuf,
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Download records and assemblies from BV-BRC.
    Fetch {
        #[command(subcommand)]
        what: Fetch,
    },
    /// Consolidate phenotypes, run genome QC, find datasets, and choose the
    /// genomes to download for each species.
    Curate {
        /// Most genomes to select per species.
        #[arg(long, default_value_t = 2000)]
        cap: usize,
        /// Minimum genomes per class for a species x drug dataset, counting
        /// unconflicted reported calls only.
        #[arg(long, default_value_t = 25)]
        min_class: usize,
        /// SigDrug vocabulary directory
        /// (`data/drug/compendium/facts/vocabulary` in a SigDrug checkout).
        /// Without it, names resolve through GRM's vocabulary alone.
        #[arg(long, env = "GRM_SIGDRUG_VOCABULARY")]
        sigdrug_vocabulary: Option<PathBuf>,
        /// How genomes are chosen up to the cap: `legacy`, the label-weighted
        /// ranking of `curate::select_genomes`, or `group-round-robin-v1`,
        /// the label-free group draw of `curate::plan_group_sampling` with
        /// each genome's BV-BRC MLST sequence type as its group (a genome
        /// without one is its own singleton group).
        #[arg(long, value_enum, default_value_t = SelectionRule::Legacy)]
        selection: SelectionRule,
        /// How a genome and drug's reported calls become one call when its
        /// records mix categories: `r-or-s-wins` (the default, every release
        /// so far: S with I is S, I with R is R, R with S no call),
        /// `intermediate-aside` (a pair with an Intermediate and an R or S
        /// record gets no call), `most-resistant` (R over I over S, so R with
        /// S is R) or `intermediate-class` (a pair with an Intermediate
        /// record is I unless it holds both R and S). A rule other than the
        /// default is written to qc.json.
        #[arg(long, value_enum, default_value_t = LabelRuleArg::ROrSWins)]
        label_rule: LabelRuleArg,
        /// Seed of the group draw; required with `group-round-robin-v1`.
        #[arg(long)]
        selection_seed: Option<u64>,
        /// The always-included stratum for `group-round-robin-v1`, a JSON
        /// object {genome_ids, reason, label_dependent}; none when omitted.
        #[arg(long)]
        stratum: Option<PathBuf>,
    },
    /// Compare each curated species' selected genomes with the species'
    /// reference genome (Mash distance and GC) and record the genomes whose
    /// sequence is another species, which `grm curate` then excludes.
    SpeciesCheck {
        /// Only sketch the FASTA files on hand into the sketch cache, for a
        /// check that runs elsewhere.
        #[arg(long)]
        sketch_only: bool,
        /// Choose references afresh and judge every genome without its
        /// recorded verdict. The genomes the existing species check
        /// recorded are still judged, so the ones curation dropped stay out.
        #[arg(long)]
        reset: bool,
    },
    /// Count k-mers for each downloaded genome of the species, with
    /// assembly stats and catalog gene calls.
    Count {
        #[arg(long, value_delimiter = ',', required = true)]
        species: Vec<String>,
        #[arg(
            long,
            default_value_t = PIPELINE_K,
            value_parser = RangedU64ValueParser::<usize>::new().range(MIN_K as u64..=MAX_K as u64)
        )]
        k: usize,
        /// Count without catalog gene calls when the AMRFinderPlus catalog
        /// does not load, instead of stopping. Reports then have no catalog
        /// baseline, and the gene panels detect nothing.
        #[arg(long)]
        allow_missing_catalog: bool,
    },
    /// Lineage clusters, dendrogram, and 2D map from MinHash sketches.
    Structure {
        #[arg(long, value_delimiter = ',', required = true)]
        species: Vec<String>,
    },
    /// Merge ksets into the species pattern matrix.
    Matrix {
        #[arg(long, value_delimiter = ',', required = true)]
        species: Vec<String>,
        /// Keep k-mers present in at least this many genomes (and absent
        /// from at least this many).
        #[arg(long, default_value_t = PIPELINE_MIN_SUPPORT)]
        min_support: u32,
        /// Merge this many genomes at a time through column-block partials
        /// instead of opening every kset at once. A multiple of 64. The
        /// matrix files are the same; the ksets are kept. Unlike the default
        /// merge, it refuses a kset whose genome or digest differs from
        /// counted.json (docs/matrix-scaling.md).
        #[arg(long, value_parser = parse_column_blocks)]
        column_blocks: Option<usize>,
        /// Count and merge one column block at a time, keeping 30 GiB free.
        #[arg(long, conflicts_with = "column_blocks")]
        disk_sparse: bool,
        /// Genomes per sparse block; a multiple of 64 (default 256).
        #[arg(long, requires = "disk_sparse", value_parser = parse_block_width)]
        block_width: Option<usize>,
        /// With --disk-sparse, count without catalog gene calls when the
        /// AMRFinderPlus catalog does not load, instead of stopping (see
        /// `grm count`). The sparse run counts, then deletes, its ksets.
        #[arg(long, requires = "disk_sparse")]
        allow_missing_catalog: bool,
    },
    /// Cross-validate every learner on each species x drug dataset, fit the
    /// final models, annotate their features, and write one report each.
    Learn {
        #[arg(long, value_delimiter = ',', required = true)]
        species: Vec<String>,
        /// Only these antibiotics or dataset ids.
        #[arg(long, value_delimiter = ',')]
        only: Vec<String>,
        /// Rerun datasets that already have a report.
        #[arg(long)]
        force: bool,
        /// Retained-feature budget: chi-square-screened features per fit
        /// (`EvalConfig::screened_features`).
        #[arg(long, default_value_t = DEFAULT_SCREENED_FEATURES, value_parser = RangedU64ValueParser::<usize>::new().range(1..))]
        screened_features: usize,
        /// Score the one dataset this fold file names on the outer folds it
        /// supplies, in place of the random and lineage splits
        /// (docs/fold-file.md). Writes <fold-out>/<dataset>.<scheme>.json,
        /// which records the file's SHA-256, and no report. Takes one
        /// species; --force replaces an existing result. It runs at the
        /// release's seeds and tie rule, so --seed, --learner-seed and the
        /// --exp-tie-* options are refused with it. --exp-null-refit is
        /// refused with it too: it records the folds grm learn draws.
        #[arg(
            long,
            requires = "fold_out",
            conflicts_with_all = [
                "only",
                "seed",
                "learner_seed",
                "exp_tie_rule",
                "exp_tie_trace",
                "exp_tie_out",
                "exp_null_refit",
            ]
        )]
        fold_file: Option<PathBuf>,
        /// Directory for the --fold-file result.
        #[arg(long, requires = "fold_file")]
        fold_out: Option<PathBuf>,
        /// Seed of the evaluation's random draws (`EvalConfig::seed`): both
        /// splits' outer folds, each fit's inner folds and tuning holdouts,
        /// and the forest's and boosting's own draws. The species matrix,
        /// the cohort, the lineages and the screening rule do not depend on
        /// it. Every release ran at the default. Another value scores
        /// another fold draw and writes its reports under
        /// <data>/species/<slug>/draws/seed-<seed>/, which `grm report`
        /// does not read (docs/fold-draws.md).
        #[arg(long, default_value_t = RELEASE_SEED)]
        seed: u64,
        /// Seed the forest's and boosting's own draws with this instead of
        /// --seed, every partition unchanged (`EvalConfig::learner_seed`).
        /// Reports go under draws/seed-<seed>-learner-<learner-seed>/.
        #[arg(long)]
        learner_seed: Option<u64>,
        /// Experiment, off by default: how the decision tree, the forest and
        /// the boosted trees break ties between equally good features
        /// (`grm::learn::tie`). `float` is production's rule. `exact` makes
        /// the stated tie rule true for cart and rf; for gbt, only between
        /// splits with the same two groups. A report fitted with `exact` says
        /// so in `eval_config.tie_rule` and its input fingerprint hashes it,
        /// so `grm learn` refits it and the metric audit flags it. The site
        /// reads neither field: `grm report` stops at such a report, and at
        /// the release's seeds an exact run needs --exp-tie-out.
        #[arg(long, value_enum, default_value_t = TieRuleArg::Float)]
        exp_tie_rule: TieRuleArg,
        /// Experiment: write one JSON line per tree-learner fit to this file,
        /// with its nodes searched, the nodes where the two tie rules pick
        /// different features, and a digest of each tree. The fits are the
        /// same as without it.
        #[arg(long)]
        exp_tie_trace: Option<PathBuf>,
        /// Experiment: write this run's reports to this folder in place of
        /// results/ or the draw's folder. An --exp-tie-rule exact run at the
        /// release's seeds needs it: it refuses results/, the folder
        /// production writes and `grm report` indexes.
        #[arg(long)]
        exp_tie_out: Option<PathBuf>,
        /// Experiment, off by default (docs/null-fold-refit.md): also refit
        /// every boosting and lasso outer fold whose fit is null under the
        /// fallback settings, and append each one's scores to this JSON-lines
        /// file. While the file takes every line, the reports equal those of
        /// a run without it. A line it does not take ends the learn with an
        /// error naming the file and the fold, before that dataset's report
        /// is written. Refused beside --fold-file, and beside the --exp-tie-*
        /// options: its refits take production's tie rule and would add
        /// fits to a tie trace.
        #[arg(long, conflicts_with_all = ["exp_tie_rule", "exp_tie_trace", "exp_tie_out"])]
        exp_null_refit: Option<PathBuf>,
        /// Experiment, off by default (docs/mic-plan.md): fit only these MIC
        /// models, a comma list of `released` (the report's censored MIC
        /// regression) and `ordered` (one fit per dilution threshold), on
        /// the lineage split at --seed, or with --fold-file on the file's
        /// folds. Writes one result per dataset, to --exp-mic-out or with
        /// --fold-file to --fold-out, and fits no classifier and no final
        /// model; no report is read or written. --force replaces a result.
        /// Refused beside --learner-seed, the --exp-tie-* options and
        /// --exp-null-refit.
        #[arg(
            long,
            value_enum,
            value_delimiter = ',',
            conflicts_with_all = [
                "learner_seed",
                "exp_tie_rule",
                "exp_tie_trace",
                "exp_tie_out",
                "exp_null_refit",
            ]
        )]
        exp_mic: Vec<MicArmArg>,
        /// Directory for the --exp-mic results on the lineage split:
        /// <dir>/<dataset>.lineage-s<seed>.mic.json. With --fold-file the
        /// result goes to --fold-out, so this is refused beside it.
        #[arg(long, requires = "exp_mic", conflicts_with = "fold_file")]
        exp_mic_out: Option<PathBuf>,
        /// Experiment, off by default (docs/cohort-draw.md, plan block C1
        /// part B): for each drawn cohort file in this folder
        /// (<dataset>.r<r>.d<draw>.cohort.json) of the species' datasets, or
        /// of --only's, run the release's rule on the cohort at --seed and
        /// score its final models on the dataset's other genomes. Writes one
        /// result of counts per cohort to --exp-cohort-out and no report; a
        /// result that exists is kept unless --force. Refused beside
        /// --fold-file, --learner-seed, the --exp-tie-* options,
        /// --exp-null-refit and --exp-mic.
        #[arg(
            long,
            requires = "exp_cohort_out",
            conflicts_with_all = [
                "fold_file",
                "learner_seed",
                "exp_tie_rule",
                "exp_tie_trace",
                "exp_tie_out",
                "exp_null_refit",
                "exp_mic",
            ]
        )]
        exp_cohort: Option<PathBuf>,
        /// Directory for the --exp-cohort results:
        /// <dir>/<dataset>.r<r>.d<draw>.s<seed>.cohort.json.
        #[arg(long, requires = "exp_cohort")]
        exp_cohort_out: Option<PathBuf>,
    },
    /// Run fetch, count, structure, matrix, and learn for each species, then
    /// build the report. Stages skip work that is already on disk.
    Pipeline {
        /// Species slugs; every curated species when omitted.
        #[arg(long, value_delimiter = ',')]
        species: Vec<String>,
        #[arg(long, default_value_t = 8)]
        workers: usize,
        /// k-mer length to count, cluster and learn at. Stages built at
        /// another k are recounted and rebuilt at this one.
        #[arg(
            long,
            default_value_t = PIPELINE_K,
            value_parser = RangedU64ValueParser::<usize>::new().range(MIN_K as u64..=MAX_K as u64)
        )]
        k: usize,
        /// Retained-feature budget (see `grm learn`).
        #[arg(long, default_value_t = DEFAULT_SCREENED_FEATURES, value_parser = RangedU64ValueParser::<usize>::new().range(1..))]
        screened_features: usize,
        /// Delete each species' ksets after its datasets are learned, to cap
        /// disk use at one species' ksets (re-creatable with `grm count`).
        #[arg(long, conflicts_with = "disk_sparse")]
        drop_ksets: bool,
        /// Count and merge one column block at a time, keeping 30 GiB free.
        #[arg(long)]
        disk_sparse: bool,
        /// Genomes per sparse block; a multiple of 64 (default 256).
        #[arg(long, requires = "disk_sparse", value_parser = parse_block_width)]
        block_width: Option<usize>,
        /// Count without catalog gene calls when the AMRFinderPlus catalog
        /// does not load, instead of stopping (see `grm count`).
        #[arg(long)]
        allow_missing_catalog: bool,
    },
    /// Score one arm of an input comparison on one dataset: its feature
    /// blocks, encoded and assembled per outer fold, and every learner on the
    /// reference arm's fold plans (`grm::ablate::arm`). Writes contract §4
    /// prediction rows, run records and plans under --out. The species must
    /// be counted and built (count, structure, matrix) in the data root.
    Arm {
        #[arg(long)]
        species: String,
        #[arg(long)]
        dataset: String,
        /// The arm's name in prediction rows and contrasts.
        #[arg(long)]
        arm: String,
        /// Feature blocks: K, C, H, P (N and X have no source).
        #[arg(long, value_delimiter = ',', default_value = "K")]
        blocks: Vec<String>,
        /// Retained-feature budget (see `grm learn`).
        #[arg(long, default_value_t = DEFAULT_SCREENED_FEATURES, value_parser = RangedU64ValueParser::<usize>::new().range(1..))]
        screened_features: usize,
        #[arg(long, default_value_t = 1, value_parser = RangedU64ValueParser::<u32>::new().range(1..=100))]
        repeats: u32,
        /// How this arm pairs with the reference arm.
        #[arg(long, value_enum, default_value_t = PairingArg::SamePlan)]
        pairing: PairingArg,
        /// The reference arm's --out folder; omit for the reference itself.
        #[arg(long)]
        reference: Option<PathBuf>,
        /// A cohort arm's balanced-subsample seed (default: the reference
        /// plan's repeat seed).
        #[arg(long)]
        balance_seed: Option<u64>,
        /// Seed the learners' own draws with this instead of the plan's
        /// repeat seed, every partition unchanged (the calibration twin).
        #[arg(long)]
        learner_seed: Option<u64>,
        /// Score the reference's outer folds with inner folds and tuning
        /// holdouts drawn from this seed (the retuned twin; same_plan only).
        #[arg(long)]
        tuning_seed: Option<u64>,
        /// Write training.jsonl: each outer fold's training genome ids.
        #[arg(long)]
        trace_training: bool,
        /// The species grouping every arm splits on: the structure.json of
        /// one `grm structure` over the species' eligible pool (a cohort arm
        /// needs it).
        #[arg(long)]
        groups: Option<PathBuf>,
        #[arg(long)]
        out: PathBuf,
    },
    /// Apply one dataset report's final models, frozen, to genomes, as
    /// docs/prediction-contract.md defines: one JSON line per genome and
    /// learner, with each model feature's published k-mer found or not
    /// found. Fits nothing and reads no network.
    Predict {
        /// A dataset report that `grm learn` wrote or `grm report` served.
        #[arg(long)]
        report: PathBuf,
        /// Assemblies (FASTA, gzip or plain), counted at the report's k.
        #[arg(long, value_delimiter = ',')]
        fasta: Vec<PathBuf>,
        /// Ksets that `grm count` wrote at the report's k.
        #[arg(long, value_delimiter = ',')]
        kset: Vec<PathBuf>,
        /// The species matrix the report was learned from
        /// (<data>/species/<slug>/matrix). With it, each row also counts how
        /// many of every feature's k-mers the genome carries.
        #[arg(long)]
        matrix: Option<PathBuf>,
        /// Learners to apply; every complete one when omitted. A learner
        /// whose published form is not a complete model is refused.
        #[arg(long, value_enum, value_delimiter = ',')]
        learners: Vec<LearnerArg>,
        /// Output file; standard output when omitted.
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Build the dashboard data under <data>/report.
    Report {
        /// Published Kover results to compare against. The default is
        /// relative to the working folder; a missing file is an error.
        #[arg(long, default_value = KOVER_RESULTS)]
        kover: PathBuf,
    },
    /// Write each dataset's projected matrix, its random and lineage fold
    /// plans, and the forest's screened rows of each outer fold under
    /// --out/<dataset id>/, for a learner outside GRM
    /// (docs/dataset-export.md). Fits nothing. The species must be counted
    /// and built (count, structure, matrix).
    Export {
        #[arg(long)]
        species: String,
        /// Only these antibiotics or dataset ids. A name that matches no
        /// dataset of the species stops the export.
        #[arg(long, value_delimiter = ',')]
        only: Vec<String>,
        /// Retained-feature budget (see `grm learn`).
        #[arg(long, default_value_t = DEFAULT_SCREENED_FEATURES, value_parser = RangedU64ValueParser::<usize>::new().range(1..))]
        screened_features: usize,
        #[arg(long)]
        out: PathBuf,
        /// Replace matrix.bin and dataset.json in dataset directories that
        /// already exist under --out, instead of stopping.
        #[arg(long)]
        force: bool,
    },
    /// Write offline 95% Clopper-Pearson intervals for each learner row's
    /// VME and ME under <data>/report/intervals/, one file per dataset report
    /// that index.json lists, each pinned to its report's SHA-256 and counts.
    /// No report file changes.
    Intervals,
    /// Write what two predictors that read no genome reach on each dataset
    /// report's own folds under <data>/report/baselines/, one file per
    /// dataset report that index.json lists, each pinned to its report's
    /// SHA-256 and evaluation revision: constants at the training folds'
    /// median and mode MIC beside the MIC model, and lineage identity beside
    /// each split; and <data>/report/baselines.json for the pages that show
    /// many datasets, pinned to index.json's SHA-256. A dataset whose rebuilt
    /// MIC folds contradict its report is refused alone, and the command
    /// then exits nonzero. Run it after `grm report`. Fits nothing; no report
    /// file changes.
    Baselines,
    /// Write each indexed dataset report's provenance under
    /// <data>/report/provenance/: its samples' BV-BRC dates (`grm fetch
    /// genome-dates`) and record details (curated/provenance.jsonl), pinned to
    /// the report's SHA-256 and evaluation revision. No report file changes.
    Provenance {
        /// Write nothing; fail unless every provenance file's pins match the
        /// report and the inputs as they are now.
        #[arg(long)]
        check: bool,
    },
    /// Write the dashboard's research file: the revision-5 MIC baseline and
    /// selective-calling results under <data>/research, each with a guard on
    /// the revision-5 cohort it describes.
    Research {
        /// The revision-5 dataset reports that
        /// <data>/research/RESEARCH-SNAPSHOT-SHA256SUMS lists.
        #[arg(long)]
        revision5_reports: PathBuf,
        /// Output file, served beside the report as research.json.
        #[arg(long)]
        out: PathBuf,
    },
    /// Serve the built dashboard (web/dist) with its report data.
    Serve {
        #[arg(long, default_value = "web/dist")]
        app: PathBuf,
        #[arg(long, default_value_t = 8640)]
        port: u16,
    },
}

#[derive(Subcommand)]
enum Fetch {
    /// Every laboratory-method AMR record.
    Amr,
    /// Metadata for every genome that has an AMR record.
    Genomes,
    /// The NCBI AMRFinderPlus reference gene catalog.
    Catalog,
    /// BV-BRC's date_inserted and completion_date for every genome in
    /// species/*/genomes.txt, into bvbrc/genome_dates.jsonl with a manifest
    /// that pins it. genomes.jsonl is not changed.
    GenomeDates,
    /// Assemblies for the genomes `curate` selected.
    Fasta {
        /// Species slugs (from curated/species.json); all when omitted.
        #[arg(long, value_delimiter = ',')]
        species: Vec<String>,
        #[arg(long, default_value_t = 8)]
        workers: usize,
    },
}

fn main() -> Result<()> {
    run(Cli::parse())
}

fn run(cli: Cli) -> Result<()> {
    let layout = Layout::new(&cli.data);
    match cli.command {
        Command::Fetch { what } => fetch(&layout, what),
        Command::Curate {
            cap,
            min_class,
            sigdrug_vocabulary,
            selection,
            label_rule,
            selection_seed,
            stratum,
        } => {
            let rule = selection_rule(selection, selection_seed, stratum)?;
            run_curate_with(
                &layout,
                cap,
                min_class,
                sigdrug_vocabulary,
                &rule,
                label_rule.into(),
            )
        }
        Command::SpeciesCheck { sketch_only, reset } => species_check(&layout, sketch_only, reset),
        Command::Count {
            species,
            k,
            allow_missing_catalog,
        } => {
            require_catalog(&layout, k, allow_missing_catalog)?;
            species.iter().try_for_each(|s| count(&layout, s, k))
        }
        Command::Structure { species } => species.iter().try_for_each(|s| structure(&layout, s)),
        Command::Matrix {
            species,
            min_support,
            column_blocks,
            disk_sparse,
            block_width,
            allow_missing_catalog,
        } => {
            if disk_sparse {
                let width = block_width.unwrap_or(256);
                species.iter().try_for_each(|s| {
                    matrix_sparse(&layout, s, min_support, width, allow_missing_catalog)
                })
            } else {
                species
                    .iter()
                    .try_for_each(|s| matrix_with(&layout, s, min_support, column_blocks))
            }
        }
        Command::Learn {
            species,
            only,
            force,
            screened_features,
            seed,
            fold_file,
            fold_out,
            exp_mic,
            exp_mic_out,
            ..
        } if !exp_mic.is_empty() => learn_mic(
            &layout,
            &species,
            &only,
            force,
            screened_features,
            seed,
            fold_file.as_deref().zip(fold_out.as_deref()),
            exp_mic_out.as_deref(),
            &exp_mic,
        ),
        Command::Learn {
            species,
            only,
            force,
            screened_features,
            seed,
            exp_cohort: Some(cohorts),
            exp_cohort_out: Some(out),
            ..
        } => learn_cohorts(
            &layout,
            &species,
            &only,
            force,
            screened_features,
            seed,
            &cohorts,
            &out,
        ),
        Command::Learn {
            species,
            force,
            screened_features,
            fold_file: Some(fold_file),
            fold_out: Some(fold_out),
            ..
        } => learn_fold_file(
            &layout,
            &species,
            &fold_file,
            &fold_out,
            force,
            screened_features,
        ),
        Command::Learn {
            species,
            only,
            force,
            screened_features,
            seed,
            learner_seed,
            exp_tie_rule,
            exp_tie_trace,
            exp_tie_out,
            exp_null_refit,
            ..
        } => {
            let cfg = learn_config(screened_features, seed, learner_seed, exp_tie_rule.into());
            // Every species is checked before the first fit or trace line.
            for s in &species {
                refuse_exact_in_results(&layout, s, &cfg, exp_tie_out.as_deref())?;
            }
            if let Some(path) = &exp_tie_trace {
                grm::learn::tie::start_trace(path)
                    .with_context(|| format!("create {}", path.display()))?;
            }
            if let Some(path) = &exp_null_refit {
                grm::learn::null_refit::enable(path)?;
            }
            species
                .iter()
                .try_for_each(|s| learn_in(&layout, s, &only, force, &cfg, exp_tie_out.as_deref()))
        }
        Command::Pipeline {
            species,
            workers,
            k,
            screened_features,
            drop_ksets,
            disk_sparse,
            block_width,
            allow_missing_catalog,
        } => {
            require_catalog(&layout, k, allow_missing_catalog)?;
            require_kover(std::path::Path::new(KOVER_RESULTS))?;
            let slugs = if species.is_empty() {
                curated_species(&layout)?
            } else {
                species
            };
            let sparse = disk_sparse.then_some(block_width.unwrap_or(256));
            // The pipeline is the release's command: it always draws at the
            // release's seed and breaks ties by production's rule.
            let cfg = grm::learn::eval::EvalConfig {
                screened_features,
                ..Default::default()
            };
            for slug in &slugs {
                // Annotation reads reference assemblies, so the FASTA files
                // are fetched even when the matrix is already built.
                fetch(
                    &layout,
                    Fetch::Fasta {
                        species: vec![slug.clone()],
                        workers,
                    },
                )?;
                count_to_matrix_with_mode(&layout, slug, k, sparse, allow_missing_catalog)?;
                let learned = if disk_sparse {
                    learn_disk_sparse(&layout, slug, &cfg)
                } else {
                    learn(&layout, slug, &[], false, &cfg)
                };
                if let Err(e) = learned {
                    if !e.is::<grm::probes::KsetsUnreadable>() {
                        return Err(e);
                    }
                    if disk_sparse {
                        return Err(e.context(
                            "disk-sparse ksets were deleted after the matrix save; this species has kset-backed probe assays, so its probe wiring still needs a reviewed cache rather than an unbounded stock recount",
                        ));
                    }
                    eprintln!("{slug}: recounting ksets: {e:#}");
                    count(&layout, slug, k)?;
                    count_to_matrix(&layout, slug, k)?;
                    learn(&layout, slug, &[], false, &cfg)?;
                }
                if drop_ksets {
                    let dir = layout.kset_dir(slug);
                    if dir.exists() {
                        std::fs::remove_dir_all(&dir)?;
                        eprintln!("{slug}: removed ksets ({})", dir.display());
                    }
                }
                // Publish completed species while the remaining species run.
                report(&layout, &PathBuf::from(KOVER_RESULTS))?;
            }
            Ok(())
        }
        Command::Arm {
            species,
            dataset,
            arm,
            blocks,
            screened_features,
            repeats,
            pairing,
            reference,
            balance_seed,
            learner_seed,
            tuning_seed,
            trace_training,
            groups,
            out,
        } => {
            let cfg = grm::ablate::arm::ArmConfig {
                arm,
                blocks,
                eval: grm::learn::eval::EvalConfig {
                    screened_features,
                    ..Default::default()
                },
                repeats,
                pairing: pairing.into(),
                reference,
                balance_seed,
                learner_seed,
                tuning_seed,
                trace_training,
                groups,
            };
            let started = std::time::Instant::now();
            grm::ablate::arm::run_arm(&layout, &species, &dataset, &cfg, &out)?;
            eprintln!(
                "{dataset}: arm {} scored in {:.0}s; written to {}",
                cfg.arm,
                started.elapsed().as_secs_f64(),
                out.display()
            );
            Ok(())
        }
        Command::Predict {
            report,
            fasta,
            kset,
            matrix,
            learners,
            out,
        } => predict(&report, &fasta, &kset, matrix.as_deref(), &learners, out),
        Command::Report { kover } => report(&layout, &kover),
        Command::Export {
            species,
            only,
            screened_features,
            out,
            force,
        } => {
            let cfg = grm::learn::eval::EvalConfig {
                screened_features,
                ..Default::default()
            };
            let done = grm::export::export_species(&layout, &species, &only, &out, &cfg, force)?;
            let replaced = done.iter().filter(|d| d.replaced).count();
            eprintln!(
                "{species}: {} datasets exported to {}{}",
                done.len(),
                out.display(),
                if replaced > 0 {
                    format!(", {replaced} of them replacing an earlier export")
                } else {
                    String::new()
                }
            );
            Ok(())
        }
        Command::Intervals => {
            let n = grm::intervals::write(&layout.root.join("report"))?;
            eprintln!("intervals written for {n} dataset reports");
            Ok(())
        }
        Command::Baselines => {
            let w = grm::baselines::write(&layout.root.join("report"))?;
            for r in &w.refused {
                eprintln!("{r}");
            }
            eprintln!(
                "baselines written for {} dataset reports: {} with an MIC model, {} of them with folds their held-out predictions confirm; {} exact-MIC test isolates on {} datasets have a k-mer-identical copy in their training fold",
                w.files, w.mic, w.confirmed, w.copies, w.copy_datasets
            );
            w.outcome()
        }
        Command::Provenance { check: false } => {
            let n = grm::provenance::write(&layout)?;
            eprintln!("provenance written for {n} dataset reports");
            Ok(())
        }
        Command::Provenance { check: true } => {
            let n = grm::provenance::check(&layout)?;
            eprintln!("provenance pins match for {n} dataset reports");
            Ok(())
        }
        Command::Research {
            revision5_reports,
            out,
        } => research(&layout, &revision5_reports, &out),
        Command::Serve { app, port } => grm::serve::serve(&app, &layout.root.join("report"), port),
    }
}

const KOVER_RESULTS: &str = "reference/kover-summary.json";

/// The `grm pipeline --k` and `grm count` default.
const PIPELINE_K: usize = 31;

/// min_support for `grm pipeline`, and the `grm matrix` default.
const PIPELINE_MIN_SUPPORT: u32 = 2;

/// The `--screened-features` default: `EvalConfig::default()`'s budget.
const DEFAULT_SCREENED_FEATURES: usize = 4000;

/// The `grm learn --seed` default: `EvalConfig::default()`'s seed, which
/// every release's folds were drawn from.
const RELEASE_SEED: u64 = 42;

/// The evaluation `grm learn` runs: the release's settings, with the
/// screening budget, the seeds and the tie rule its flags give.
fn learn_config(
    screened_features: usize,
    seed: u64,
    learner_seed: Option<u64>,
    tie_rule: grm::learn::tie::TieRule,
) -> grm::learn::eval::EvalConfig {
    grm::learn::eval::EvalConfig {
        screened_features,
        seed,
        learner_seed,
        tie_rule,
        ..Default::default()
    }
}

/// Bytes a disk-sparse run keeps free: `DISK_SPARSE_RESERVE` in every build
/// but this binary's own tests, which keep none, so that they pass with
/// little free space under `TMPDIR`.
const SPARSE_RESERVE: u64 = if cfg!(test) {
    0
} else {
    grm::species::DISK_SPARSE_RESERVE
};

fn parse_block_width(raw: &str) -> std::result::Result<usize, String> {
    let width: usize = raw
        .parse()
        .map_err(|_| format!("block width must be an integer, not {raw:?}"))?;
    if width >= 64 && width.is_multiple_of(64) && width <= grm::matrix::REVIEWED_COLUMNS {
        Ok(width)
    } else {
        Err(format!(
            "block width must be a multiple of 64 in 64..={}, not {width}",
            grm::matrix::REVIEWED_COLUMNS
        ))
    }
}

/// Count, structure, and matrix at `k`, each rerun only when `stale_stages`
/// finds its output out of date.
fn count_to_matrix(layout: &Layout, slug: &str, k: usize) -> Result<()> {
    count_to_matrix_with_mode(layout, slug, k, None, false)
}

/// `count_to_matrix`, or with `sparse` (`grm pipeline --disk-sparse`) the
/// disk-sparse matrix at that block width and then the structure.
fn count_to_matrix_with_mode(
    layout: &Layout,
    slug: &str,
    k: usize,
    sparse: Option<usize>,
    allow_missing_catalog: bool,
) -> Result<()> {
    if let Some(width) = sparse {
        let catalog = sparse_catalog(layout, k, allow_missing_catalog)?;
        matrix_sparse_with_catalog(
            layout,
            slug,
            PIPELINE_MIN_SUPPORT,
            k,
            catalog.as_ref(),
            width,
        )?;
        let stale =
            grm::species::stale_stages(layout, slug, catalog.as_ref(), k, PIPELINE_MIN_SUPPORT);
        if stale.structure {
            structure_disk_sparse(layout, slug, SPARSE_RESERVE)?;
        } else {
            eprintln!("{slug}: structure matches the counted genomes, k, and ksets; reusing it");
        }
        return Ok(());
    }
    let catalog = load_catalog(layout, k);
    let stale_stages =
        || grm::species::stale_stages(layout, slug, catalog.as_ref(), k, PIPELINE_MIN_SUPPORT);
    // Count reruns whenever any stage is stale (see `stale_stages`).
    if !stale_stages().count {
        eprintln!(
            "{slug}: counts match genomes.txt and the catalog, and structure and matrix match the counts; skipping to learn"
        );
        return Ok(());
    }
    count_with(layout, slug, k, catalog.as_ref())?;
    // `stale_stages` takes a genome whose count failed as left out on
    // purpose, so stages built now would leave it out of every later run. A
    // failed count, such as a kset write on a full disk, stops here instead:
    // structure and matrix then do not match these counts, and the next run
    // counts again.
    let counted = grm::species::load_counted(layout, slug)?;
    if let Some((id, why)) = counted.failed.first() {
        anyhow::bail!(
            "{slug}: {} of {} selected genomes failed to count, first {id}: {why}; fix the cause and rerun, which counts them again before structure, matrix, and learning",
            counted.failed.len(),
            counted.failed.len() + counted.genomes.len()
        );
    }
    // Structure and matrix each record the counts they were built from, so
    // comparing them with the new counts alone decides each rebuild.
    let stale = stale_stages();
    if stale.structure {
        structure(layout, slug)?;
    } else {
        eprintln!("{slug}: structure matches the counted genomes, k, and ksets; reusing it");
    }
    if stale.matrix {
        matrix(layout, slug, PIPELINE_MIN_SUPPORT)
    } else {
        eprintln!("{slug}: matrix matches the counted genomes, k, and ksets; reusing it");
        Ok(())
    }
}

fn count(layout: &Layout, slug: &str, k: usize) -> Result<()> {
    count_with(layout, slug, k, load_catalog(layout, k).as_ref())
}

fn load_catalog(layout: &Layout, k: usize) -> Option<Catalog> {
    let space = grm::dna::KmerSpace::new(k);
    match Catalog::load(&layout.catalog_dir(), &space) {
        Ok(c) => Some(c),
        Err(e) => {
            eprintln!("no catalog gene calls: {e:#}");
            None
        }
    }
}

/// Stops `grm count` and `grm pipeline` when the catalog does not load,
/// unless `allow_missing_catalog`: counting without it records no gene calls.
fn require_catalog(layout: &Layout, k: usize, allow_missing_catalog: bool) -> Result<()> {
    if !allow_missing_catalog {
        Catalog::load(&layout.catalog_dir(), &grm::dna::KmerSpace::new(k))
            .context(CATALOG_REQUIRED)?;
    }
    Ok(())
}

const CATALOG_REQUIRED: &str = "the AMRFinderPlus catalog does not load, and counting without it records no gene calls, which leaves reports without the catalog baseline and gene panel detections; pass --allow-missing-catalog to count anyway";

/// Stops `grm pipeline` before it counts anything when the Kover comparison
/// at `path` doesn't load: each species ends in `grm report`, which would
/// refuse only after the learning.
fn require_kover(path: &std::path::Path) -> Result<()> {
    grm::report::kover_benchmark(path, grm::report::PIPELINE_KOVER_REMEDY).map(|_| ())
}

fn count_with(layout: &Layout, slug: &str, k: usize, catalog: Option<&Catalog>) -> Result<()> {
    let started = std::time::Instant::now();
    let s = grm::species::count(layout, slug, k, catalog)?;
    eprintln!(
        "{slug}: {} counted, {} missing FASTA, {} failed in {:.0}s",
        s.genomes.len(),
        s.missing_fasta.len(),
        s.failed.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

fn structure(layout: &Layout, slug: &str) -> Result<()> {
    let started = std::time::Instant::now();
    let s = grm::species::structure(layout, slug)?;
    log_structure(slug, &s, started);
    Ok(())
}

/// `structure` in `grm pipeline --disk-sparse` (see
/// `grm::species::structure_disk_sparse`).
fn structure_disk_sparse(layout: &Layout, slug: &str, reserve: u64) -> Result<()> {
    let started = std::time::Instant::now();
    let s = grm::species::structure_disk_sparse(layout, slug, reserve)?;
    log_structure(slug, &s, started);
    Ok(())
}

fn log_structure(slug: &str, s: &grm::species::Structure, started: std::time::Instant) {
    eprintln!(
        "{slug}: {} genomes, {} lineage clusters at Mash distance {} (largest {:.0}%) in {:.0}s",
        s.genomes.len(),
        s.lineages.n_clusters,
        s.lineages.threshold,
        s.lineages.largest_share * 100.0,
        started.elapsed().as_secs_f64()
    );
}

/// `grm matrix --column-blocks`: a width `build_matrix` accepts.
fn parse_column_blocks(s: &str) -> std::result::Result<usize, String> {
    let width = s.parse::<usize>().map_err(|e| e.to_string())?;
    grm::species::check_column_blocks(width).map_err(|e| e.to_string())?;
    Ok(width)
}

fn matrix(layout: &Layout, slug: &str, min_support: u32) -> Result<()> {
    matrix_with(layout, slug, min_support, None)
}

fn matrix_with(
    layout: &Layout,
    slug: &str,
    min_support: u32,
    column_blocks: Option<usize>,
) -> Result<()> {
    let started = std::time::Instant::now();
    let s = grm::species::build_matrix(layout, slug, min_support, column_blocks)?;
    log_matrix(slug, &s, started);
    Ok(())
}

fn matrix_sparse(
    layout: &Layout,
    slug: &str,
    min_support: u32,
    width: usize,
    allow_missing_catalog: bool,
) -> Result<()> {
    // A fresh sparse matrix builds counted.json incrementally. Existing
    // manifests keep their k; only a new one uses the count/pipeline default.
    let k = grm::species::load_counted(layout, slug).map_or(PIPELINE_K, |c| c.k);
    // A sparse matrix counts, so it passes the catalog gate `grm count` does
    // before anything is counted, written or deleted.
    require_catalog(layout, k, allow_missing_catalog)?;
    let catalog = sparse_catalog(layout, k, allow_missing_catalog)?;
    matrix_sparse_with_catalog(layout, slug, min_support, k, catalog.as_ref(), width)
}

/// The catalog a disk-sparse count calls genes with. The sparse stages
/// delete the ksets they count, so a catalog that does not load is an error
/// here too, as in `require_catalog`, unless `allow_missing_catalog`: never
/// a silent count without gene calls.
fn sparse_catalog(
    layout: &Layout,
    k: usize,
    allow_missing_catalog: bool,
) -> Result<Option<Catalog>> {
    match Catalog::load(&layout.catalog_dir(), &grm::dna::KmerSpace::new(k)) {
        Ok(c) => Ok(Some(c)),
        Err(e) if allow_missing_catalog => {
            eprintln!("no catalog gene calls: {e:#}");
            Ok(None)
        }
        Err(e) => Err(e.context(CATALOG_REQUIRED)),
    }
}

fn matrix_sparse_with_catalog(
    layout: &Layout,
    slug: &str,
    min_support: u32,
    k: usize,
    catalog: Option<&Catalog>,
    width: usize,
) -> Result<()> {
    let started = std::time::Instant::now();
    let s = grm::species::build_matrix_disk_sparse(
        layout,
        slug,
        k,
        catalog,
        min_support,
        width,
        SPARSE_RESERVE,
    )?;
    log_matrix(slug, &s, started);
    Ok(())
}

fn log_matrix(slug: &str, s: &grm::matrix::MatrixSummary, started: std::time::Instant) {
    eprintln!(
        "{slug}: {} genomes, {} distinct k-mers, {} kept in {} patterns ({} rare, {} core) in {:.0}s",
        s.genomes.len(),
        s.kmers_distinct,
        s.kmers_kept,
        s.n_patterns,
        s.kmers_rare,
        s.kmers_core,
        started.elapsed().as_secs_f64()
    );
}

fn learn(
    layout: &Layout,
    slug: &str,
    only: &[String],
    force: bool,
    cfg: &grm::learn::eval::EvalConfig,
) -> Result<()> {
    learn_in(layout, slug, only, force, cfg, None)
}

/// `learn` as `grm pipeline --disk-sparse` runs it: the reports' directory
/// and each report are written only through real directories below the data
/// root (`grm::analysis::analyze_species_disk_sparse`).
fn learn_disk_sparse(
    layout: &Layout,
    slug: &str,
    cfg: &grm::learn::eval::EvalConfig,
) -> Result<()> {
    grm::species::ensure_selection_counted(layout, slug)?;
    let written = grm::analysis::analyze_species_disk_sparse(layout, slug, &[], false, cfg)?;
    eprintln!("{}", learn_end_line(layout, slug, written.len(), cfg));
    Ok(())
}

/// `learn`, with the reports in `out` when `grm learn --exp-tie-out` names
/// it. That run ends with a line of its own, which does not begin with the
/// release line `scripts/metric-audit/audit.py` reads.
fn learn_in(
    layout: &Layout,
    slug: &str,
    only: &[String],
    force: bool,
    cfg: &grm::learn::eval::EvalConfig,
    out: Option<&std::path::Path>,
) -> Result<()> {
    grm::species::ensure_selection_counted(layout, slug)?;
    let line = match out {
        None => {
            let written = grm::analysis::analyze_species(layout, slug, only, force, cfg)?;
            learn_end_line(layout, slug, written.len(), cfg)
        }
        Some(dir) => {
            let written = grm::analysis::analyze_species_in(layout, slug, only, force, cfg, dir)?;
            format!(
                "{slug}: {} experiment reports written to {}",
                written.len(),
                dir.display()
            )
        }
    };
    eprintln!("{line}");
    Ok(())
}

/// An `--exp-tie-rule exact` run stops before any fit when its reports would
/// go to `results/`, the folder production reports are written to and
/// `grm report` indexes: by default at the release's seeds, or because
/// `--exp-tie-out` names that folder.
fn refuse_exact_in_results(
    layout: &Layout,
    slug: &str,
    cfg: &grm::learn::eval::EvalConfig,
    out: Option<&std::path::Path>,
) -> Result<()> {
    if cfg.tie_rule.is_float() {
        return Ok(());
    }
    let results = layout.results_dir(slug);
    let dir = out.map_or_else(|| grm::analysis::reports_dir(layout, slug, cfg), Into::into);
    let same = resolved(&dir) == resolved(&results);
    anyhow::ensure!(
        !same,
        "--exp-tie-rule exact would write {slug}'s reports to {}, the folder production reports are written to and `grm report` indexes; name another folder with --exp-tie-out DIR",
        results.display()
    );
    Ok(())
}

/// `path` made absolute, each part that exists resolved through the file
/// system, and a `..` after a missing part taken off in the text. Two
/// spellings of one folder then compare equal before the folder exists,
/// when `canonicalize` has nothing to resolve.
fn resolved(path: &std::path::Path) -> PathBuf {
    let absolute = std::path::absolute(path).unwrap_or_else(|_| path.to_path_buf());
    let mut out = PathBuf::new();
    let mut exists = true;
    for part in absolute.components() {
        if !exists && part == std::path::Component::ParentDir {
            out.pop();
            exists = out.exists();
            continue;
        }
        out.push(part);
        if exists {
            match std::fs::canonicalize(&out) {
                Ok(real) => out = real,
                Err(_) => exists = false,
            }
        }
    }
    out
}

/// The line `grm learn` ends a species with. A run that writes `results/`
/// ends with the release's own line, which `scripts/metric-audit/audit.py`
/// reads as the species finished. A draw names its folder, in words that
/// line's pattern cannot match, so a log that mixes a draw with a release
/// run never marks the species finished on the draw's line.
fn learn_end_line(
    layout: &Layout,
    slug: &str,
    written: usize,
    cfg: &grm::learn::eval::EvalConfig,
) -> String {
    let dir = grm::analysis::reports_dir(layout, slug, cfg);
    if dir == layout.results_dir(slug) {
        format!("{slug}: {written} dataset reports written")
    } else {
        format!(
            "{slug}: {written} draw reports written to {}",
            dir.display()
        )
    }
}

/// `grm learn --fold-file`: one dataset of one species on the outer folds a
/// file supplies (`analysis::evaluate_fold_file`), at the release's seeds
/// (the parser refuses `--seed` and `--learner-seed` beside it). It ends
/// with its own line, not `learn_end_line`.
fn learn_fold_file(
    layout: &Layout,
    species: &[String],
    fold_file: &std::path::Path,
    fold_out: &std::path::Path,
    force: bool,
    screened_features: usize,
) -> Result<()> {
    let [slug] = species else {
        anyhow::bail!(
            "--fold-file names one dataset, so it takes one species, not {}",
            species.len()
        );
    };
    grm::species::ensure_selection_counted(layout, slug)?;
    let cfg = grm::learn::eval::EvalConfig {
        screened_features,
        ..Default::default()
    };
    let written =
        grm::analysis::evaluate_fold_file(layout, slug, fold_file, fold_out, force, &cfg)?;
    eprintln!("{slug}: fold file result written to {}", written.display());
    Ok(())
}

/// `grm learn --exp-mic`: the MIC models on the lineage split at `seed`
/// (`analysis::mic_plan_species`), or on a fold file's folds at the
/// release's seeds (`analysis::mic_plan_fold_file`; the parser refuses
/// `--seed` beside `--fold-file`). Ends with lines of its own, which do not
/// begin with the release line `learn_end_line` writes.
#[allow(clippy::too_many_arguments)]
fn learn_mic(
    layout: &Layout,
    species: &[String],
    only: &[String],
    force: bool,
    screened_features: usize,
    seed: u64,
    fold_file: Option<(&std::path::Path, &std::path::Path)>,
    exp_mic_out: Option<&std::path::Path>,
    arms: &[MicArmArg],
) -> Result<()> {
    let arms: Vec<grm::learn::micplan::Arm> = arms.iter().map(|&a| a.into()).collect();
    let cfg = learn_config(
        screened_features,
        seed,
        None,
        grm::learn::tie::TieRule::Float,
    );
    if let Some((file, out)) = fold_file {
        let [slug] = species else {
            anyhow::bail!(
                "--fold-file names one dataset, so it takes one species, not {}",
                species.len()
            );
        };
        grm::species::ensure_selection_counted(layout, slug)?;
        let written =
            grm::analysis::mic_plan_fold_file(layout, slug, file, out, force, &cfg, &arms)?;
        eprintln!("{slug}: MIC plan result written to {}", written.display());
        return Ok(());
    }
    let out = exp_mic_out
        .context("--exp-mic needs --exp-mic-out DIR, or --fold-file FILE with --fold-out DIR")?;
    for slug in species {
        grm::species::ensure_selection_counted(layout, slug)?;
    }
    for slug in species {
        let written = grm::analysis::mic_plan_species(layout, slug, only, force, &cfg, &arms, out)?;
        eprintln!(
            "{slug}: {} MIC plan results written to {}",
            written.len(),
            out.display()
        );
    }
    Ok(())
}

/// `grm learn --exp-cohort`: the drawn cohorts of each species' datasets at
/// `seed` (`analysis::cohort_draw_species`). Ends with lines of its own,
/// which do not begin with the release line `learn_end_line` writes.
#[allow(clippy::too_many_arguments)]
fn learn_cohorts(
    layout: &Layout,
    species: &[String],
    only: &[String],
    force: bool,
    screened_features: usize,
    seed: u64,
    cohorts: &std::path::Path,
    out: &std::path::Path,
) -> Result<()> {
    let cfg = learn_config(
        screened_features,
        seed,
        None,
        grm::learn::tie::TieRule::Float,
    );
    for slug in species {
        grm::species::ensure_selection_counted(layout, slug)?;
    }
    for slug in species {
        let (written, kept) =
            grm::analysis::cohort_draw_species(layout, slug, only, force, &cfg, cohorts, out)?;
        eprintln!(
            "{slug}: {} cohort results written to {}, {kept} kept",
            written.len(),
            out.display()
        );
    }
    Ok(())
}

fn report(layout: &Layout, kover: &std::path::Path) -> Result<()> {
    let out = grm::report::build(layout, kover)?;
    eprintln!("dashboard data written to {}", out.display());
    Ok(())
}

fn research(
    layout: &Layout,
    revision5_reports: &std::path::Path,
    out: &std::path::Path,
) -> Result<()> {
    let r = grm::research::build(
        &layout.root.join("research"),
        revision5_reports,
        grm::research::REVISION5_SNAPSHOT,
    )?;
    grm::research::write(&r, out)?;
    eprintln!(
        "{} MIC datasets, {} withheld, and {} selective-calling datasets from snapshot {} written to {}",
        r.mic.datasets.len(),
        r.mic.withheld.len(),
        r.selective.datasets.len(),
        r.snapshot.sha256,
        out.display()
    );
    Ok(())
}

fn fetch(layout: &Layout, what: Fetch) -> Result<()> {
    let client = Client::new();
    match what {
        Fetch::Amr => {
            let out = layout.amr();
            ensure_parent(&out)?;
            let n = client.fetch_amr(&out)?;
            eprintln!("wrote {n} AMR records to {}", out.display());
        }
        Fetch::Genomes => {
            let records: Vec<AmrRecord> = bvbrc::read_jsonl(&layout.amr())?;
            let ids: Vec<String> = records
                .into_iter()
                .map(|r| r.genome_id)
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            let out = layout.genomes();
            ensure_parent(&out)?;
            let n = client.fetch_genomes(&ids, &out)?;
            eprintln!(
                "wrote metadata for {n} of {} genomes to {}",
                ids.len(),
                out.display()
            );
        }
        Fetch::Catalog => {
            let dir = layout.catalog_dir();
            std::fs::create_dir_all(&dir)?;
            for file in grm::catalog::CATALOG_FILES {
                client.download(
                    &format!("{}/{file}", grm::catalog::CATALOG_URL),
                    &dir.join(file),
                )?;
                eprintln!("  {file}");
            }
            eprintln!(
                "catalog version {}",
                std::fs::read_to_string(dir.join("version.txt"))?.trim()
            );
        }
        Fetch::GenomeDates => {
            let lists = grm::provenance::genome_lists(layout)?;
            let ids: Vec<String> = lists
                .iter()
                .flat_map(|l| l.ids.iter().cloned())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            let rows = client.fetch_genome_dates(&ids)?;
            let manifest = grm::provenance::write_genome_dates(layout, &lists, &ids, rows)?;
            eprintln!(
                "wrote dates for {} of {} genomes to {} ({} missing, {} without date_inserted); manifest {}",
                manifest.returned,
                manifest.requested,
                layout.genome_dates().display(),
                manifest.missing.len(),
                manifest.without_date_inserted,
                layout.genome_dates_manifest().display()
            );
        }
        Fetch::Fasta { species, workers } => {
            let slugs = if species.is_empty() {
                curated_species(layout)?
            } else {
                species
            };
            for slug in slugs {
                let ids = read_lines(&layout.species_dir(&slug).join("genomes.txt"))?;
                eprintln!("{slug}: {} genomes", ids.len());
                let failures = client.fetch_fasta_many(&ids, &layout.fasta_dir(), workers)?;
                if !failures.is_empty() {
                    let log = layout.species_dir(&slug).join("fasta_failures.tsv");
                    let mut w = std::fs::File::create(&log)?;
                    for (id, err) in &failures {
                        writeln!(w, "{id}\t{err}")?;
                    }
                    eprintln!(
                        "{slug}: {} downloads failed; see {}",
                        failures.len(),
                        log.display()
                    );
                }
            }
        }
    }
    Ok(())
}

#[derive(Serialize)]
struct SpeciesSelection {
    species: String,
    slug: String,
    datasets: usize,
    labeled_genomes: usize,
    selected_genomes: usize,
    /// Selected genomes the species check has no verdict for, typically
    /// ones not downloaded when it ran; None without a species check.
    selected_unchecked: Option<usize>,
}

#[derive(Serialize)]
struct DrugNameRow<'a> {
    records: usize,
    class: String,
    #[serde(flatten)]
    resolution: &'a drugs::Resolution,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
enum PairingArg {
    #[value(name = "same_plan")]
    SamePlan,
    #[value(name = "k31_groups")]
    K31Groups,
    #[value(name = "common_outer_test_set_balanced")]
    CommonOuterTestSetBalanced,
}

impl From<PairingArg> for grm::ablate::arm::Pairing {
    fn from(p: PairingArg) -> Self {
        match p {
            PairingArg::SamePlan => Self::SamePlan,
            PairingArg::K31Groups => Self::K31Groups,
            PairingArg::CommonOuterTestSetBalanced => Self::CommonOuterTestSetBalanced,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
enum TieRuleArg {
    Float,
    Exact,
}

impl From<TieRuleArg> for grm::learn::tie::TieRule {
    fn from(t: TieRuleArg) -> Self {
        match t {
            TieRuleArg::Float => Self::Float,
            TieRuleArg::Exact => Self::Exact,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
enum MicArmArg {
    Released,
    Ordered,
}

impl From<MicArmArg> for grm::learn::micplan::Arm {
    fn from(a: MicArmArg) -> Self {
        match a {
            MicArmArg::Released => Self::Released,
            MicArmArg::Ordered => Self::Ordered,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
enum LearnerArg {
    Scm,
    Cart,
    Lasso,
}

impl From<LearnerArg> for grm::predict::Learner {
    fn from(l: LearnerArg) -> Self {
        match l {
            LearnerArg::Scm => Self::Scm,
            LearnerArg::Cart => Self::Cart,
            LearnerArg::Lasso => Self::Lasso,
        }
    }
}

/// `grm predict`: rows for every genome, FASTA files first, in the order
/// given. Rows are written once every genome is scored, so a refused genome
/// leaves no partial output.
fn predict(
    report: &std::path::Path,
    fasta: &[PathBuf],
    kset: &[PathBuf],
    matrix: Option<&std::path::Path>,
    learners: &[LearnerArg],
    out: Option<PathBuf>,
) -> Result<()> {
    anyhow::ensure!(
        !fasta.is_empty() || !kset.is_empty(),
        "give genomes with --fasta or --kset"
    );
    let mut model = grm::predict::Model::load(report)?;
    let requested: Vec<grm::predict::Learner> = learners.iter().map(|&l| l.into()).collect();
    let learners = model.learners(&requested)?;
    for (learner, why) in &model.incomplete {
        eprintln!(
            "{}: {} left out, its published form is not a complete model: {why}",
            model.dataset,
            learner.name()
        );
    }
    if let Some(dir) = matrix {
        model.attach_matrix(dir)?;
    }
    let mut lines = Vec::new();
    for (path, is_fasta) in fasta
        .iter()
        .map(|p| (p, true))
        .chain(kset.iter().map(|p| (p, false)))
    {
        let genome = if is_fasta {
            model.genome_from_fasta(path)?
        } else {
            model.genome_from_kset(path)?
        };
        for row in model.predict(&genome, &learners) {
            serde_json::to_writer(&mut lines, &row)?;
            lines.push(b'\n');
        }
    }
    match out {
        Some(path) => {
            ensure_parent(&path)?;
            std::fs::write(&path, &lines).with_context(|| format!("write {}", path.display()))?;
        }
        None => std::io::stdout().write_all(&lines)?,
    }
    eprintln!(
        "{}: {} genomes scored with {}",
        model.dataset,
        fasta.len() + kset.len(),
        learners
            .iter()
            .map(|l| l.name())
            .collect::<Vec<_>>()
            .join(", ")
    );
    Ok(())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
enum SelectionRule {
    /// `curate::select_genomes`.
    Legacy,
    /// `curate::plan_group_sampling` (`curate::GROUP_SAMPLING_DESIGN`).
    #[value(name = "group-round-robin-v1")]
    GroupRoundRobinV1,
}

/// `curate::LabelRule` on the command line, with the same names as in qc.json.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
enum LabelRuleArg {
    #[value(name = "r-or-s-wins")]
    ROrSWins,
    #[value(name = "intermediate-aside")]
    IntermediateAside,
    #[value(name = "most-resistant")]
    MostResistant,
    #[value(name = "intermediate-class")]
    IntermediateClass,
}

impl From<LabelRuleArg> for curate::LabelRule {
    fn from(arg: LabelRuleArg) -> Self {
        match arg {
            LabelRuleArg::ROrSWins => curate::LabelRule::ROrSWins,
            LabelRuleArg::IntermediateAside => curate::LabelRule::IntermediateAside,
            LabelRuleArg::MostResistant => curate::LabelRule::MostResistant,
            LabelRuleArg::IntermediateClass => curate::LabelRule::IntermediateClass,
        }
    }
}

/// A selection rule with its inputs, checked before curation reads anything.
enum Selector {
    Legacy,
    /// The group draw; each genome's group is its BV-BRC MLST sequence type
    /// (`group_of_mlst`).
    Groups {
        seed: u64,
        stratum: curate::AlwaysIncluded,
    },
}

fn selection_rule(
    rule: SelectionRule,
    seed: Option<u64>,
    stratum: Option<PathBuf>,
) -> Result<Selector> {
    match rule {
        SelectionRule::Legacy => {
            anyhow::ensure!(
                seed.is_none() && stratum.is_none(),
                "--selection-seed and --stratum apply only to --selection group-round-robin-v1"
            );
            Ok(Selector::Legacy)
        }
        SelectionRule::GroupRoundRobinV1 => {
            let seed = seed.context("--selection group-round-robin-v1 needs --selection-seed")?;
            let stratum = match stratum {
                Some(p) => serde_json::from_slice(
                    &std::fs::read(&p).with_context(|| format!("read {}", p.display()))?,
                )
                .with_context(|| format!("parse the stratum {}", p.display()))?,
                None => curate::AlwaysIncluded::default(),
            };
            Ok(Selector::Groups { seed, stratum })
        }
    }
}

/// A genome's group for the group draw: `mlst:<sequence type>` from its
/// BV-BRC genome record, or `singleton:<genome id>` when the record has no
/// sequence type, so every genome has a group (the rulings after lane INPUT).
fn group_of_mlst(genome_id: &str, mlst: Option<&str>) -> (String, bool) {
    match mlst.map(str::trim).filter(|m| !m.is_empty()) {
        Some(st) => (format!("mlst:{st}"), false),
        None => (format!("singleton:{genome_id}"), true),
    }
}

/// The genomes `selector` chooses for `species` at `cap`, in the order
/// genomes.txt keeps, with the record of how. A group draw also returns its
/// plan, which `run_curate_with` saves beside genomes.txt. `mlst` maps a
/// genome to its BV-BRC MLST field.
fn select_species(
    curated: &curate::Curated,
    species: &str,
    cap: usize,
    selector: &Selector,
    mlst: &HashMap<&str, &str>,
) -> Result<(
    Vec<String>,
    curate::Selection,
    Option<curate::GroupSamplingPlan>,
)> {
    match selector {
        Selector::Legacy => Ok((
            curate::select_genomes(&curated.phenotypes, &curated.datasets, species, cap),
            curate::Selection {
                rule: curate::LEGACY_SELECTION.to_string(),
                cap,
                seed: None,
                always_included: None,
                input_sha256: None,
                group_source: None,
                singletons: None,
            },
            None,
        )),
        Selector::Groups { seed, stratum } => {
            let ids = curate::labeled_genomes(&curated.phenotypes, &curated.datasets, species);
            let mut singletons = 0;
            let eligible: Vec<curate::EligibleGenome> = ids
                .iter()
                .map(|id| {
                    let (group_id, single) = group_of_mlst(id, mlst.get(id.as_str()).copied());
                    singletons += single as usize;
                    curate::EligibleGenome {
                        genome_id: id.clone(),
                        group_id,
                        always_included: None,
                    }
                })
                .collect();
            // One plan per species (`plan_group_sampling`), with the part of
            // the stratum that is this species'.
            let own: BTreeSet<&str> = ids.iter().map(String::as_str).collect();
            let always = curate::AlwaysIncluded {
                genome_ids: stratum
                    .genome_ids
                    .iter()
                    .filter(|id| own.contains(id.as_str()))
                    .cloned()
                    .collect(),
                ..stratum.clone()
            };
            let plan = curate::plan_group_sampling(&eligible, &always, cap, *seed)
                .with_context(|| format!("{species}: group-stratified selection"))?;
            let selection = curate::Selection {
                rule: plan.design.clone(),
                cap,
                seed: Some(*seed),
                always_included: Some(plan.always_included.clone()),
                input_sha256: Some(plan.input_sha256.clone()),
                group_source: Some(curate::MLST_GROUP_SOURCE.to_string()),
                singletons: Some(singletons),
            };
            Ok((plan.selected.clone(), selection, Some(plan)))
        }
    }
}

#[cfg(test)]
fn run_curate(
    layout: &Layout,
    cap: usize,
    min_class: usize,
    sigdrug: Option<PathBuf>,
) -> Result<()> {
    run_curate_with(
        layout,
        cap,
        min_class,
        sigdrug,
        &Selector::Legacy,
        curate::LabelRule::default(),
    )
}

fn run_curate_with(
    layout: &Layout,
    cap: usize,
    min_class: usize,
    sigdrug: Option<PathBuf>,
    selector: &Selector,
    label_rule: curate::LabelRule,
) -> Result<()> {
    // `grm learn` gates each dataset at this minimum; below the outer fold
    // count, a random-split fold lacks a class and `grm learn` stops with an
    // error. At or above it, moving duplicate groups whole
    // (`learn::isolate_groups`) can still empty a fold of one class in a
    // small cohort with many duplicates (docs/duplicate-isolation.md).
    let folds = grm::learn::eval::EvalConfig::default().folds;
    anyhow::ensure!(
        min_class >= folds,
        "--min-class {min_class} is below {folds}: each of the {folds} random-split folds needs a resistant and a susceptible genome"
    );
    let records: Vec<AmrRecord> = bvbrc::read_jsonl(&layout.amr())?;
    let genomes: Vec<GenomeRecord> = bvbrc::read_jsonl(&layout.genomes())?;
    let mlst: HashMap<&str, &str> = genomes
        .iter()
        .filter_map(|g| Some((g.genome_id.as_str(), g.mlst.as_deref()?)))
        .collect();
    let cfg = QcConfig {
        min_class,
        label_rule,
        ..Default::default()
    };
    let vocab = match &sigdrug {
        Some(dir) => drugs::Vocabulary::load_sigdrug(dir)?,
        None => {
            eprintln!(
                "no SigDrug vocabulary given; resolving antibiotic names with GRM's vocabulary only"
            );
            drugs::Vocabulary::default()
        }
    };
    let mut resolver = drugs::Resolver::new(vocab);
    let check_path = layout.species_check();
    let check = if check_path.exists() {
        Some(grm::species_check::read(&check_path)?)
    } else {
        eprintln!(
            "no species check at {}; curating without one (run `grm species-check` to drop genomes of another species)",
            check_path.display()
        );
        None
    };
    let wrong_species = check.as_ref().map(|c| c.failures()).unwrap_or_default();
    let mut curated = curate::curate(&records, &genomes, &cfg, &wrong_species, &mut resolver)?;

    let dir = layout.curated();
    std::fs::create_dir_all(&dir)?;
    let mut raw_counts: BTreeMap<&str, usize> = BTreeMap::new();
    for r in &records {
        *raw_counts.entry(r.antibiotic.as_str()).or_default() += 1;
    }
    let mut names: Vec<DrugNameRow> = raw_counts
        .iter()
        .map(|(raw, &n)| {
            let resolution = &resolver.resolutions[*raw];
            DrugNameRow {
                records: n,
                class: drugs::class_of(&resolution.canonical),
                resolution,
            }
        })
        .collect();
    names.sort_by_key(|a| std::cmp::Reverse(a.records));
    std::fs::write(
        dir.join("drug_names.json"),
        serde_json::to_vec_pretty(&serde_json::json!({
            "vocabulary": resolver.vocabulary_source(),
            "names": names,
        }))?,
    )?;
    let mut w = std::io::BufWriter::new(std::fs::File::create(dir.join("phenotypes.jsonl"))?);
    for p in &curated.phenotypes {
        serde_json::to_writer(&mut w, p)?;
        w.write_all(b"\n")?;
    }
    w.flush()?;
    // A side file for `grm provenance`; every other curated file is as before.
    let mut w = std::io::BufWriter::new(std::fs::File::create(dir.join("provenance.jsonl"))?);
    for row in &curated.provenance {
        serde_json::to_writer(&mut w, row)?;
        w.write_all(b"\n")?;
    }
    w.flush()?;
    std::fs::write(
        dir.join("funnel.json"),
        serde_json::to_vec_pretty(&curated.funnel)?,
    )?;
    /// qc.json's first element: the QC settings and the selection rule.
    #[derive(Serialize)]
    struct QcRecord<'a> {
        #[serde(flatten)]
        qc: &'a QcConfig,
        selection: curate::Selection,
    }
    let (seed, always_included, group_source) = match selector {
        Selector::Legacy => (None, None, None),
        Selector::Groups { seed, stratum } => (
            Some(*seed),
            Some(stratum.clone()),
            Some(curate::MLST_GROUP_SOURCE.to_string()),
        ),
    };
    let rule = match selector {
        Selector::Legacy => curate::LEGACY_SELECTION,
        Selector::Groups { .. } => curate::GROUP_SAMPLING_DESIGN,
    };
    let record = QcRecord {
        qc: &cfg,
        selection: curate::Selection {
            rule: rule.to_string(),
            cap,
            seed,
            always_included,
            input_sha256: None,
            group_source,
            singletons: None,
        },
    };
    std::fs::write(
        dir.join("qc.json"),
        serde_json::to_vec_pretty(&(&record, &curated.qc_failures))?,
    )?;

    let mut by_species: BTreeMap<String, usize> = BTreeMap::new();
    for d in &curated.datasets {
        *by_species.entry(d.species.clone()).or_default() += 1;
    }
    let mut selections = Vec::new();
    let mut stratum_used: BTreeSet<String> = BTreeSet::new();
    for (species, n_datasets) in by_species {
        let species = species.as_str();
        let slug = curate::slug(species);
        let labeled = curated
            .phenotypes
            .iter()
            .filter(|p| p.species == species && p.label.is_some())
            .map(|p| &p.genome_id)
            .collect::<BTreeSet<_>>()
            .len();
        let (ids, selection, plan) = select_species(&curated, species, cap, selector, &mlst)?;
        if let Some(a) = &selection.always_included {
            stratum_used.extend(a.genome_ids.iter().cloned());
        }
        let plan_path = layout.species_dir(&slug).join("selection.json");
        match &plan {
            Some(plan) => {
                ensure_parent(&plan_path)?;
                std::fs::write(&plan_path, serde_json::to_vec_pretty(plan)?)?;
            }
            None if plan_path.exists() => std::fs::remove_file(&plan_path)?,
            None => {}
        }
        for d in curated.datasets.iter_mut().filter(|d| d.species == species) {
            d.selection = Some(selection.clone());
        }
        let selected_unchecked = check.as_ref().map(|c| {
            let checked: BTreeSet<&str> = c
                .genomes
                .iter()
                .filter(|g| g.species == species)
                .map(|g| g.genome_id.as_str())
                .collect();
            ids.iter()
                .filter(|id| !checked.contains(id.as_str()))
                .count()
        });
        let path = layout.species_dir(&slug).join("genomes.txt");
        ensure_parent(&path)?;
        std::fs::write(&path, ids.join("\n") + "\n")?;
        selections.push(SpeciesSelection {
            species: species.to_string(),
            slug,
            datasets: n_datasets,
            labeled_genomes: labeled,
            selected_genomes: ids.len(),
            selected_unchecked,
        });
    }
    if let Selector::Groups { stratum, .. } = selector {
        let unused: Vec<&String> = stratum
            .genome_ids
            .iter()
            .filter(|id| !stratum_used.contains(*id))
            .collect();
        anyhow::ensure!(
            unused.is_empty(),
            "{} always-included genomes are eligible in no curated species, first {}",
            unused.len(),
            unused.first().map_or("", |s| s.as_str())
        );
    }
    std::fs::write(
        dir.join("datasets.json"),
        serde_json::to_vec_pretty(&curated.datasets)?,
    )?;
    selections.sort_by_key(|a| std::cmp::Reverse(a.labeled_genomes));
    std::fs::write(
        dir.join("species.json"),
        serde_json::to_vec_pretty(&selections)?,
    )?;

    eprintln!("{}", serde_json::to_string_pretty(&curated.funnel)?);
    eprintln!(
        "{} datasets across {} species",
        curated.datasets.len(),
        selections.len()
    );
    for s in &selections {
        eprintln!(
            "  {:<40} {:>3} datasets  {:>6} labeled  {:>5} selected",
            s.species, s.datasets, s.labeled_genomes, s.selected_genomes
        );
    }
    Ok(())
}

fn species_check(layout: &Layout, sketch_only: bool, reset: bool) -> Result<()> {
    #[derive(serde::Deserialize)]
    struct Row {
        species: String,
        slug: String,
    }
    let path = layout.curated().join("species.json");
    let rows: Vec<Row> = serde_json::from_slice(
        &std::fs::read(&path)
            .with_context(|| format!("run `grm curate` first ({})", path.display()))?,
    )?;
    let selected = rows
        .into_iter()
        .map(|r| {
            let ids = read_lines(&layout.species_dir(&r.slug).join("genomes.txt"))?;
            Ok((r.species, r.slug, ids))
        })
        .collect::<Result<Vec<_>>>()?;
    // Genomes an earlier check judged stay in, even after --reset: curation
    // has since dropped the failures from genomes.txt. Without --reset its
    // references and verdicts stay too.
    let earlier_path = layout.species_check();
    let earlier = if earlier_path.exists() {
        Some(grm::species_check::read(&earlier_path)?)
    } else {
        None
    };
    let lists = grm::species_check::cohort_lists(selected, earlier.as_ref());
    let started = std::time::Instant::now();
    let cohorts = grm::species_check::load_cohorts(layout, &lists)?;
    let sketched: usize = cohorts.iter().map(|c| c.sketches.len()).sum();
    let missing: usize = cohorts.iter().map(|c| c.unchecked.len()).sum();
    eprintln!(
        "{sketched} genomes sketched, {missing} without FASTA or sketch, in {:.0}s",
        started.elapsed().as_secs_f64()
    );
    if sketch_only {
        return Ok(());
    }
    if reset && let Some(e) = &earlier {
        grm::species_check::ensure_recorded_failures_judgeable(&cohorts, e)?;
    }
    let recorded = if reset { None } else { earlier.as_ref() };
    let check = grm::species_check::check(cohorts, &Default::default(), recorded)?;
    let out = grm::species_check::write(layout, &check)?;
    for s in &check.species {
        eprintln!(
            "  {:<40} {:>5} checked  {:>3} failed  {:>3} borderline  {:>4} unchecked  reference {}",
            s.species,
            s.checked,
            s.failed,
            s.borderline,
            s.unchecked.len(),
            s.reference
        );
    }
    for p in &check.pairs {
        let outcome = if p.same_group {
            "borderline, one genomic group"
        } else {
            "fail"
        };
        eprintln!(
            "  {} nearer {}: {} genomes ({:.2}%), {outcome}",
            p.labelled,
            p.nearer,
            p.genomes,
            p.share * 100.0
        );
    }
    for g in &check.genomes {
        if let Some(now) = &g.kept_fail {
            eprintln!(
                "  kept failed {} ({}): recorded as failed, this check judged {now}; --reset to rejudge",
                g.genome_id, g.species
            );
        }
        if g.verdict == grm::species_check::Verdict::Fail {
            let why = g.reason.as_deref().unwrap_or("");
            eprintln!("  fail {} ({}): {why}", g.genome_id, g.species);
        }
    }
    eprintln!("species check written to {}", out.display());
    Ok(())
}

fn curated_species(layout: &Layout) -> Result<Vec<String>> {
    #[derive(serde::Deserialize)]
    struct Row {
        slug: String,
    }
    let path = layout.curated().join("species.json");
    let rows: Vec<Row> = serde_json::from_slice(
        &std::fs::read(&path)
            .with_context(|| format!("run `grm curate` first ({})", path.display()))?,
    )?;
    Ok(rows.into_iter().map(|r| r.slug).collect())
}

fn read_lines(path: &std::path::Path) -> Result<Vec<String>> {
    Ok(std::fs::read_to_string(path)
        .with_context(|| format!("read {}", path.display()))?
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .map(String::from)
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use grm::matrix::FeatureMatrix;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;
    use std::path::Path;

    #[test]
    fn label_rule_names_on_the_command_line_are_the_names_qc_json_records() {
        use clap::ValueEnum;
        for arg in LabelRuleArg::value_variants() {
            let name = arg.to_possible_value().unwrap().get_name().to_string();
            let rule: curate::LabelRule = (*arg).into();
            assert_eq!(serde_json::to_value(rule).unwrap(), name);
        }
        let cli = Cli::try_parse_from(["grm", "curate"]).unwrap();
        let Command::Curate { label_rule, .. } = cli.command else {
            panic!("not curate")
        };
        assert!(curate::LabelRule::from(label_rule).is_default());
    }

    #[test]
    fn a_curate_run_under_another_label_rule_records_it_in_qc_json() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        species_check_fixture(dir.path());
        let folds = grm::learn::eval::EvalConfig::default().folds;
        let qc = |layout: &Layout| curated_json(layout, "qc.json");
        run_curate(&layout, 30, folds, None).unwrap();
        assert!(qc(&layout)[0].get("label_rule").is_none());
        run_curate_with(
            &layout,
            30,
            folds,
            None,
            &Selector::Legacy,
            curate::LabelRule::IntermediateAside,
        )
        .unwrap();
        assert_eq!(qc(&layout)[0]["label_rule"], "intermediate-aside");
    }

    #[test]
    fn count_takes_only_a_k_the_kmer_space_supports() {
        let parse = |k: &str| Cli::try_parse_from(["grm", "count", "--species", "s", "--k", k]);
        for k in ["4", "33", "40"] {
            let err = parse(k).err().unwrap();
            assert_eq!(err.kind(), clap::error::ErrorKind::ValueValidation, "k={k}");
        }
        for k in [MIN_K, PIPELINE_K, MAX_K] {
            let Command::Count { k: parsed, .. } = parse(&k.to_string()).unwrap().command else {
                panic!("parsed another command");
            };
            assert_eq!(parsed, k);
        }
    }

    #[test]
    fn matrix_takes_only_a_column_block_width_of_whole_words() {
        let parse = |w: &str| {
            Cli::try_parse_from(["grm", "matrix", "--species", "s", "--column-blocks", w])
        };
        for w in ["0", "32", "96", "100", "x"] {
            let err = parse(w).err().unwrap();
            assert_eq!(err.kind(), clap::error::ErrorKind::ValueValidation, "{w}");
        }
        // clap reads a leading minus as another flag.
        assert!(parse("-64").is_err());
        for w in [64, 128, 4096] {
            let Command::Matrix { column_blocks, .. } = parse(&w.to_string()).unwrap().command
            else {
                panic!("parsed another command");
            };
            assert_eq!(column_blocks, Some(w));
        }
        let plain = Cli::try_parse_from(["grm", "matrix", "--species", "s"]).unwrap();
        let Command::Matrix { column_blocks, .. } = plain.command else {
            panic!("parsed another command");
        };
        assert_eq!(column_blocks, None);
    }

    /// Without `--disk-sparse`, `grm matrix` and `grm pipeline` parse to the
    /// stock path. The sparse flags need `--disk-sparse`, and `--disk-sparse`
    /// refuses `--column-blocks` and `--drop-ksets`.
    #[test]
    fn matrix_and_pipeline_keep_the_stock_path_unless_disk_sparse_is_given() {
        use clap::error::ErrorKind::{ArgumentConflict, MissingRequiredArgument};
        let parse = |args: &[&str]| Cli::try_parse_from(args.iter().copied());
        let Command::Matrix {
            column_blocks,
            disk_sparse,
            block_width,
            allow_missing_catalog,
            ..
        } = parse(&["grm", "matrix", "--species", "s"]).unwrap().command
        else {
            panic!("parsed another command");
        };
        assert_eq!(
            (
                column_blocks,
                disk_sparse,
                block_width,
                allow_missing_catalog
            ),
            (None, false, None, false)
        );
        let Command::Pipeline {
            drop_ksets,
            disk_sparse,
            block_width,
            ..
        } = parse(&["grm", "pipeline"]).unwrap().command
        else {
            panic!("parsed another command");
        };
        assert_eq!((drop_ksets, disk_sparse, block_width), (false, false, None));
        let kind = |args: &[&str]| parse(args).err().map(|e| e.kind());
        for args in [
            &["grm", "matrix", "--species", "s", "--block-width", "128"][..],
            &["grm", "matrix", "--species", "s", "--allow-missing-catalog"][..],
            &["grm", "pipeline", "--block-width", "128"][..],
        ] {
            assert_eq!(kind(args), Some(MissingRequiredArgument), "{args:?}");
        }
        for args in [
            &[
                "grm",
                "matrix",
                "--species",
                "s",
                "--disk-sparse",
                "--column-blocks",
                "64",
            ][..],
            &["grm", "pipeline", "--disk-sparse", "--drop-ksets"][..],
        ] {
            assert_eq!(kind(args), Some(ArgumentConflict), "{args:?}");
        }
    }

    /// `grm matrix` from its parsed arguments through `run`: the switch and
    /// `--min-support` reach the block merge. Only `main`'s `Cli::parse` of
    /// the process arguments and the printing of an error are left out.
    #[test]
    fn matrix_column_blocks_reach_the_block_merge_through_run() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        let data = dir.path().to_str().unwrap();
        let run_grm = |args: &[&str]| {
            let cli = Cli::try_parse_from(["grm", "--data", data].iter().chain(args)).unwrap();
            run(cli).map_err(|e| format!("{e:#}"))
        };
        run_grm(&["count", "--species", SLUG]).unwrap();
        let mx = layout.matrix_dir(SLUG);
        let files = || {
            [
                "matrix.json",
                "patterns.bin",
                "support.bin",
                "kmer_offsets.bin",
                "kmers.bin",
            ]
            .map(|f| std::fs::read(mx.join(f)).unwrap())
        };
        // The block merge removes a partials folder it finds; the default
        // merge leaves it alone. So the folder shows which merge ran.
        let partials = layout.partials_dir(SLUG);
        let mut written = Vec::new();
        for min_support in ["1", "2"] {
            std::fs::create_dir_all(&partials).unwrap();
            std::fs::write(partials.join("left.txt"), "a stopped run").unwrap();
            let plain = ["matrix", "--species", SLUG, "--min-support", min_support];
            run_grm(&plain).unwrap();
            assert!(partials.join("left.txt").exists(), "{min_support}");
            let want = files();
            run_grm(&[&plain[..], &["--column-blocks", "64"]].concat()).unwrap();
            assert!(
                !partials.exists(),
                "--min-support {min_support} --column-blocks 64 ran the default merge"
            );
            assert!(files() == want, "--min-support {min_support}");
            written.push(want);
        }
        assert!(written[0] != written[1]);
    }

    const SLUG: &str = "fixture";
    const GENE: &str = "ATGAAACGTTACGGCATTGCAGGCTTACCGATAGCTTGACCTAGGCATCGATCGGATCCATGCAAGCTTGG";

    /// Four genomes on a shared core: g0 and g1 also carry one block and the
    /// catalog's only gene, g2 and g3 another block. The catalog is v1, and no
    /// dataset is curated, which is enough for `grm::analysis::load_inputs`.
    fn fixture(root: &Path) -> Layout {
        let layout = Layout::new(root);
        let mut rng = ChaCha8Rng::seed_from_u64(11);
        let mut random = |n: usize| -> String {
            (0..n)
                .map(|_| b"ACGT"[rng.random_range(0..4)] as char)
                .collect()
        };
        let (core, a, b) = (random(3000), random(400), random(400));
        std::fs::create_dir_all(layout.fasta_dir()).unwrap();
        for (id, seq) in [
            ("g0", format!("{core}NN{a}NN{GENE}")),
            ("g1", format!("{core}NN{a}NN{GENE}")),
            ("g2", format!("{core}NN{b}")),
            ("g3", format!("{core}NN{b}")),
        ] {
            let path = bvbrc::fasta_path(&layout.fasta_dir(), id);
            std::fs::write(path, format!(">c\n{seq}\n")).unwrap();
        }
        write_genome_list(&layout, &["g0", "g1", "g2", "g3"]);
        let catalog = layout.catalog_dir();
        std::fs::create_dir_all(&catalog).unwrap();
        std::fs::write(
            catalog.join("AMR_CDS.fa"),
            format!(">WP_1.1|NG_1.1|1|1|blaX-1|blaX|beta-lactamase\n{GENE}\n"),
        )
        .unwrap();
        std::fs::write(
            catalog.join("ReferenceGeneCatalog.txt"),
            "allele\tgene_family\tproduct_name\ttype\tclass\tsubclass\trefseq_protein_accession\tgenbank_protein_accession\n\
             blaX-1\tblaX\tbeta-lactamase\tAMR\tBETA-LACTAM\tCARBAPENEM\tWP_1.1\t\n",
        )
        .unwrap();
        set_catalog_version(&layout, "v1");
        std::fs::create_dir_all(layout.curated()).unwrap();
        std::fs::write(layout.curated().join("datasets.json"), "[]").unwrap();
        std::fs::write(layout.curated().join("phenotypes.jsonl"), "").unwrap();
        ensure_parent(&layout.genomes()).unwrap();
        std::fs::write(layout.genomes(), "").unwrap();
        layout
    }

    fn write_genome_list(layout: &Layout, ids: &[&str]) {
        let path = layout.species_dir(SLUG).join("genomes.txt");
        ensure_parent(&path).unwrap();
        std::fs::write(path, ids.join("\n") + "\n").unwrap();
    }

    /// A FASTA file for `id`, as `grm fetch fasta` saves it, with g3's
    /// sequence.
    fn download(layout: &Layout, id: &str) {
        let fasta = |id| bvbrc::fasta_path(&layout.fasta_dir(), id);
        std::fs::copy(fasta("g3"), fasta(id)).unwrap();
    }

    fn set_catalog_version(layout: &Layout, version: &str) {
        std::fs::write(layout.catalog_dir().join("version.txt"), version).unwrap();
    }

    /// Zeroes a field of structure.json and of matrix.json that no staleness
    /// check reads, so a rebuild shows as the real value coming back.
    fn mark_stages(layout: &Layout) {
        for (path, field) in stage_fields(layout) {
            let mut v: serde_json::Value =
                serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
            v[field] = 0.into();
            std::fs::write(&path, serde_json::to_vec(&v).unwrap()).unwrap();
        }
    }

    /// Whether structure.json and matrix.json still carry `mark_stages`' zero.
    fn still_marked(layout: &Layout) -> (bool, bool) {
        let [structure, matrix] = stage_fields(layout).map(|(path, field)| {
            let v: serde_json::Value =
                serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
            v[field] == 0
        });
        (structure, matrix)
    }

    fn stage_fields(layout: &Layout) -> [(PathBuf, &'static str); 2] {
        [
            (
                layout.species_dir(SLUG).join("structure.json"),
                "sketch_size",
            ),
            (
                layout.matrix_dir(SLUG).join("matrix.json"),
                "kmers_distinct",
            ),
        ]
    }

    /// counted.json as `grm count` wrote it before kset digests, with the
    /// ksets and their sidecars deleted as `--drop-ksets` does.
    fn drop_ksets_and_digests(layout: &Layout) {
        std::fs::remove_dir_all(layout.kset_dir(SLUG)).unwrap();
        counts_without_digests(layout);
    }

    fn counts_without_digests(layout: &Layout) {
        let path = layout.species_dir(SLUG).join("counted.json");
        let mut counted: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
        for g in counted["genomes"].as_array_mut().unwrap() {
            g.as_object_mut().unwrap().remove("kset_digest");
        }
        std::fs::write(&path, serde_json::to_vec(&counted).unwrap()).unwrap();
    }

    fn remove_fields(path: &Path, fields: &[&str]) {
        let mut v: serde_json::Value =
            serde_json::from_slice(&std::fs::read(path).unwrap()).unwrap();
        for field in fields {
            v.as_object_mut().unwrap().remove(*field);
        }
        std::fs::write(path, serde_json::to_vec(&v).unwrap()).unwrap();
    }

    /// counted.json, each kset's sidecar, structure.json, and matrix.json as
    /// GRM revision 6 wrote them, before kset digests and the structure's k.
    fn written_before_kset_digests(layout: &Layout) {
        counts_without_digests(layout);
        for id in ["g0", "g1", "g2", "g3"] {
            let sidecar = layout.kset_dir(SLUG).join(format!("{id}.json"));
            remove_fields(&sidecar, &["kset_digest"]);
        }
        let structure = layout.species_dir(SLUG).join("structure.json");
        remove_fields(&structure, &["k", "kset_digests"]);
        remove_fields(
            &layout.matrix_dir(SLUG).join("matrix.json"),
            &["kset_digests"],
        );
    }

    fn counted_digests(layout: &Layout) -> Vec<Option<String>> {
        let counted = grm::species::load_counted(layout, SLUG).unwrap();
        counted.genomes.into_iter().map(|g| g.kset_digest).collect()
    }

    /// Why `grm learn` refuses the species as it is on disk, if it does.
    fn learn_refusal(layout: &Layout) -> Option<String> {
        grm::analysis::load_inputs(layout, SLUG)
            .err()
            .map(|e| format!("{e:#}"))
    }

    #[test]
    fn a_new_catalog_redoes_gene_calls_and_keeps_structure_and_matrix() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        mark_stages(&layout);
        set_catalog_version(&layout, "v2");
        // `grm learn` alone refuses the old calls; the pipeline recounts first.
        let refusal = learn_refusal(&layout).unwrap();
        assert!(
            refusal.contains("has gene calls from catalog v1|group-calls-v2"),
            "{refusal}"
        );
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        assert_eq!(learn_refusal(&layout), None);
        let counted = grm::species::load_counted(&layout, SLUG).unwrap();
        for g in &counted.genomes {
            assert_eq!(g.catalog_version, "v2|group-calls-v2", "{}", g.genome_id);
        }
        let genes: Vec<usize> = counted.genomes.iter().map(|g| g.genes.len()).collect();
        assert_eq!(genes, [1, 1, 0, 0]);
        assert_eq!(still_marked(&layout), (true, true));
    }

    #[test]
    fn direct_learn_refuses_added_removed_and_reordered_genome_selections() {
        for selection in [
            &["g0", "g1", "g2", "g3", "g4"][..],
            &["g0", "g1", "g2"][..],
            &["g3", "g2", "g1", "g0"][..],
        ] {
            let dir = tempfile::tempdir().unwrap();
            let layout = fixture(dir.path());
            count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
            download(&layout, "g4");
            write_genome_list(&layout, selection);
            let err = format!(
                "{:#}",
                learn(&layout, SLUG, &[], false, &Default::default()).unwrap_err()
            );
            assert!(
                err.contains("genomes.txt adds, removes, or reorders"),
                "{err}"
            );
            assert!(err.contains("rerun `grm count --species fixture`"), "{err}");
        }
    }

    #[test]
    fn another_genome_selection_recounts_and_rebuilds_structure_and_matrix() {
        for selection in [
            // g3 replaced by g4, downloaded after the last count.
            &["g0", "g1", "g2", "g4"][..],
            &["g0", "g1", "g2", "g3", "g4"],
            &["g0", "g1", "g2"],
            &["g3", "g2", "g1", "g0"],
        ] {
            let dir = tempfile::tempdir().unwrap();
            let layout = fixture(dir.path());
            count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
            mark_stages(&layout);
            download(&layout, "g4");
            write_genome_list(&layout, selection);
            count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
            assert_eq!(still_marked(&layout), (false, false), "{selection:?}");
            let m = FeatureMatrix::load(&layout.matrix_dir(SLUG)).unwrap();
            assert_eq!(m.summary.genomes, selection);
            assert_eq!(learn_refusal(&layout), None, "{selection:?}");
        }
    }

    #[test]
    fn a_selected_genome_is_counted_once_its_fasta_is_downloaded() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        mark_stages(&layout);
        // `grm count` leaves out a genome without a FASTA file, so selecting
        // g4 before its download changes nothing, the counts included.
        let counted_json = layout.species_dir(SLUG).join("counted.json");
        let counts = std::fs::read(&counted_json).unwrap();
        write_genome_list(&layout, &["g0", "g1", "g2", "g3", "g4"]);
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        assert_eq!(std::fs::read(&counted_json).unwrap(), counts);
        assert_eq!(still_marked(&layout), (true, true));
        download(&layout, "g4");
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        assert_eq!(still_marked(&layout), (false, false));
        let m = FeatureMatrix::load(&layout.matrix_dir(SLUG)).unwrap();
        assert_eq!(m.summary.genomes, ["g0", "g1", "g2", "g3", "g4"]);
    }

    #[test]
    fn a_genome_that_fails_to_count_stops_the_pipeline_until_a_full_count() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        // A directory where `kset::write` creates g1's temporary file fails
        // that write, as a full disk would.
        let part = layout.kset_dir(SLUG).join("g1.kset.part");
        std::fs::create_dir_all(&part).unwrap();
        let err = format!(
            "{:#}",
            count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap_err()
        );
        assert!(
            err.contains("fixture: 1 of 4 selected genomes failed to count, first g1:"),
            "{err}"
        );
        assert!(!layout.species_dir(SLUG).join("structure.json").exists());
        assert!(!layout.matrix_dir(SLUG).join("matrix.json").exists());
        std::fs::remove_dir(&part).unwrap();
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        let m = FeatureMatrix::load(&layout.matrix_dir(SLUG)).unwrap();
        assert_eq!(m.summary.genomes, ["g0", "g1", "g2", "g3"]);
        assert_eq!(learn_refusal(&layout), None);
    }

    #[test]
    fn a_count_at_another_k_is_redone_and_rebuilds_only_stages_built_at_that_k() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        mark_stages(&layout);
        // `grm count --k 21` alone: the recount at k = 31 writes again the
        // ksets that structure and matrix were built from.
        count(&layout, SLUG, 21).unwrap();
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        assert_eq!(
            grm::species::load_counted(&layout, SLUG).unwrap().k,
            PIPELINE_K
        );
        assert_eq!(still_marked(&layout), (true, true));
        // Structure and matrix built from counts at k = 21, then `grm count`.
        count(&layout, SLUG, 21).unwrap();
        structure(&layout, SLUG).unwrap();
        matrix(&layout, SLUG, 2).unwrap();
        count(&layout, SLUG, PIPELINE_K).unwrap();
        mark_stages(&layout);
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        assert_eq!(still_marked(&layout), (false, false));
        let m = FeatureMatrix::load(&layout.matrix_dir(SLUG)).unwrap();
        assert_eq!(m.summary.k, PIPELINE_K);
        assert_eq!(learn_refusal(&layout), None);
    }

    #[test]
    fn stages_built_by_hand_at_another_k_or_min_support_are_rebuilt() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        // `grm count --k 21`, `grm structure`, and `grm matrix`: stages that
        // agree with each other, at another k than the pipeline's.
        count(&layout, SLUG, 21).unwrap();
        structure(&layout, SLUG).unwrap();
        matrix(&layout, SLUG, 2).unwrap();
        mark_stages(&layout);
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        let counted = grm::species::load_counted(&layout, SLUG).unwrap();
        assert_eq!(counted.k, PIPELINE_K);
        assert_eq!(still_marked(&layout), (false, false));
        // `grm matrix --min-support 1` keeps k-mers the pipeline's matrix drops.
        matrix(&layout, SLUG, 1).unwrap();
        mark_stages(&layout);
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        assert_eq!(still_marked(&layout), (true, false));
        let m = FeatureMatrix::load(&layout.matrix_dir(SLUG)).unwrap();
        assert_eq!((m.summary.k, m.summary.min_support), (PIPELINE_K, 2));
        assert_eq!(learn_refusal(&layout), None);
    }

    #[test]
    fn a_failed_matrix_resave_is_rebuilt_on_the_next_run() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        let mx = layout.matrix_dir(SLUG);
        let built = FeatureMatrix::load(&mx).unwrap();
        // The re-save rewrites three binaries, then cannot create kmers.bin.
        std::fs::remove_file(mx.join("kmers.bin")).unwrap();
        std::fs::create_dir(mx.join("kmers.bin")).unwrap();
        assert!(grm::species::build_matrix(&layout, SLUG, 2, None).is_err());
        std::fs::remove_dir(mx.join("kmers.bin")).unwrap();
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        let rebuilt = FeatureMatrix::load(&mx).unwrap();
        assert_eq!(rebuilt.kmers, built.kmers);
        assert_eq!(rebuilt.summary.genomes, built.summary.genomes);
    }

    #[test]
    fn counts_from_before_kset_digests_are_redone_and_keep_structure_and_matrix() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        let counted_with_digests = counted_digests(&layout);
        assert!(counted_with_digests.iter().all(Option::is_some));
        mark_stages(&layout);
        drop_ksets_and_digests(&layout);
        let refusal = learn_refusal(&layout).unwrap();
        assert!(refusal.contains("no kset digest for g0"), "{refusal}");
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        assert_eq!(counted_digests(&layout), counted_with_digests);
        assert_eq!(still_marked(&layout), (true, true));
        assert_eq!(learn_refusal(&layout), None);
    }

    #[test]
    fn data_written_before_kset_digests_is_recounted_and_rebuilt() {
        // The ksets are still on disk for a species that had not finished
        // learning, and `--drop-ksets` deleted them for the others.
        for drop_ksets in [false, true] {
            let dir = tempfile::tempdir().unwrap();
            let layout = fixture(dir.path());
            count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
            let counted_with_digests = counted_digests(&layout);
            mark_stages(&layout);
            written_before_kset_digests(&layout);
            if drop_ksets {
                std::fs::remove_dir_all(layout.kset_dir(SLUG)).unwrap();
            }
            assert_eq!(
                learn_refusal(&layout).unwrap(),
                "counted.json records no kset digest for g0 (it predates them); rerun `grm count --species fixture`"
            );
            count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
            let state = format!("drop_ksets = {drop_ksets}");
            assert_eq!(counted_digests(&layout), counted_with_digests, "{state}");
            assert_eq!(still_marked(&layout), (false, false), "{state}");
            assert_eq!(learn_refusal(&layout), None, "{state}");
        }
    }

    #[test]
    fn a_structure_from_before_kset_digests_is_rebuilt_beside_a_current_matrix() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        mark_stages(&layout);
        let structure_json = layout.species_dir(SLUG).join("structure.json");
        remove_fields(&structure_json, &["k", "kset_digests"]);
        assert_eq!(
            learn_refusal(&layout).unwrap(),
            "structure.json records no kset digests (it predates them); rerun `grm structure`"
        );
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        assert_eq!(still_marked(&layout), (false, true));
        assert_eq!(learn_refusal(&layout), None);
    }

    #[test]
    fn counting_and_building_the_matrix_by_hand_leave_an_old_structure_refused() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        written_before_kset_digests(&layout);
        // `grm count` and `grm matrix`, without `grm structure`.
        count(&layout, SLUG, PIPELINE_K).unwrap();
        matrix(&layout, SLUG, 2).unwrap();
        assert_eq!(
            learn_refusal(&layout).unwrap(),
            "structure.json records no kset digests (it predates them); rerun `grm structure`"
        );
        structure(&layout, SLUG).unwrap();
        assert_eq!(learn_refusal(&layout), None);
    }

    #[test]
    fn learn_refuses_a_structure_and_matrix_built_from_other_ksets() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        let mut counted = grm::species::load_counted(&layout, SLUG).unwrap();
        counted.genomes[0].kset_digest = Some("changed kset".into());
        std::fs::write(
            layout.species_dir(SLUG).join("counted.json"),
            serde_json::to_vec(&counted).unwrap(),
        )
        .unwrap();
        let other_ksets = "was built from other ksets than counted.json records";
        assert_eq!(
            learn_refusal(&layout).unwrap(),
            format!(
                "structure.json {other_ksets}; matrix.json {other_ksets}; rerun `grm structure` and `grm matrix`"
            )
        );
    }

    #[test]
    fn changed_ksets_rebuild_structure_and_matrix() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        mark_stages(&layout);
        let fasta = bvbrc::fasta_path(&layout.fasta_dir(), "g0");
        let mut sequence = std::fs::read_to_string(&fasta).unwrap();
        sequence.push('A');
        std::fs::write(fasta, sequence).unwrap();
        drop_ksets_and_digests(&layout);
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        assert_eq!(still_marked(&layout), (false, false));
        let counted = grm::species::load_counted(&layout, SLUG).unwrap();
        let matrix = FeatureMatrix::load(&layout.matrix_dir(SLUG)).unwrap();
        let expected: Vec<String> = counted
            .genomes
            .iter()
            .map(|g| g.kset_digest.clone().unwrap())
            .collect();
        assert_eq!(matrix.summary.kset_digests, expected);
        assert_eq!(learn_refusal(&layout), None);
    }

    #[test]
    fn a_catalog_that_no_longer_loads_stops_learning_and_recounting() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        std::fs::remove_dir_all(layout.catalog_dir()).unwrap();
        let refusal = learn_refusal(&layout).unwrap();
        assert!(refusal.contains("which no longer loads"), "{refusal}");
        // A recount, due here for counts from before kset digests, would
        // record no gene calls in place of the ones counted.json holds.
        drop_ksets_and_digests(&layout);
        let err = format!(
            "{:#}",
            count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap_err()
        );
        assert!(err.contains("would drop the gene calls"), "{err}");
        let counted = grm::species::load_counted(&layout, SLUG).unwrap();
        let genes: Vec<usize> = counted.genomes.iter().map(|g| g.genes.len()).collect();
        assert_eq!(genes, [1, 1, 0, 0]);
    }

    #[test]
    fn required_catalog_rejects_incomplete_inputs_unless_allowed() {
        for (file, contents) in [
            ("ReferenceGeneCatalog.txt", Some("type\n")),
            ("version.txt", None),
            ("version.txt", Some(" \t\n")),
            ("AMR_CDS.fa", Some("")),
        ] {
            let dir = tempfile::tempdir().unwrap();
            let layout = fixture(dir.path());
            require_catalog(&layout, PIPELINE_K, false).unwrap();
            let path = layout.catalog_dir().join(file);
            if let Some(contents) = contents {
                std::fs::write(&path, contents).unwrap();
            } else {
                std::fs::remove_file(&path).unwrap();
            }
            let err = require_catalog(&layout, PIPELINE_K, false).unwrap_err();
            assert!(
                err.to_string().contains("pass --allow-missing-catalog"),
                "{err:#}"
            );
            require_catalog(&layout, PIPELINE_K, true).unwrap();
        }
    }

    #[test]
    fn pipeline_refuses_a_missing_kover_file_with_its_own_remedy() {
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("reference/kover-summary.json");
        let err = format!("{:#}", require_kover(&missing).unwrap_err());
        assert!(err.contains(missing.to_str().unwrap()), "{err}");
        assert!(
            err.contains("run `grm pipeline` from the source folder"),
            "{err}"
        );
        // pipeline has no --kover, so its message must not offer one.
        assert!(!err.contains("--kover"), "{err}");
        std::fs::create_dir_all(missing.parent().unwrap()).unwrap();
        std::fs::write(&missing, b"[]").unwrap();
        require_kover(&missing).unwrap();
        // The source folder's file loads, as `cargo test` runs from there.
        require_kover(std::path::Path::new(KOVER_RESULTS)).unwrap();
    }

    #[test]
    fn counting_stops_when_the_catalog_does_not_load_unless_allowed() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        std::fs::remove_dir_all(layout.catalog_dir()).unwrap();
        let data = dir.path().to_str().unwrap();
        let run_grm = |args: &[&str]| {
            let cli = Cli::try_parse_from(["grm", "--data", data].iter().chain(args)).unwrap();
            run(cli).map_err(|e| format!("{e:#}"))
        };
        for args in [
            &["count", "--species", SLUG][..],
            &["pipeline", "--species", SLUG],
        ] {
            let err = run_grm(args).unwrap_err();
            assert!(err.contains("pass --allow-missing-catalog"), "{err}");
            assert!(err.contains("run `grm fetch catalog` first"), "{err}");
        }
        assert!(!layout.species_dir(SLUG).join("counted.json").exists());
        run_grm(&["count", "--species", SLUG, "--allow-missing-catalog"]).unwrap();
        let counted = grm::species::load_counted(&layout, SLUG).unwrap();
        assert_eq!(counted.genomes.len(), 4);
        assert!(
            counted
                .genomes
                .iter()
                .all(|g| g.genes.is_empty() && g.catalog_version.is_empty())
        );
    }

    const FOREIGN: [&str; 3] = ["562.103", "562.108", "562.115"];

    /// One species of 60 labelled genomes with FASTA files, three of them
    /// unrelated sequence filed under it. Resistant calls are the rarer
    /// class, so a cap of 50 selects all 20 resistant genomes, the foreign
    /// three among them, and 30 of the 40 susceptible ones.
    fn species_check_fixture(root: &Path) {
        let mut rng = ChaCha8Rng::seed_from_u64(9);
        let random = |rng: &mut ChaCha8Rng| -> Vec<u8> {
            (0..100_000)
                .map(|_| b"ACGT"[rng.random_range(0..4)])
                .collect()
        };
        let species_root = random(&mut rng);
        let (mut amr, mut meta) = (String::new(), String::new());
        std::fs::create_dir_all(root.join("fasta")).unwrap();
        std::fs::create_dir_all(root.join("bvbrc")).unwrap();
        for i in 0..60 {
            let (id, call) = if i < 20 {
                (format!("562.{}", 100 + i), "Resistant")
            } else {
                (format!("562.{}", 180 + i), "Susceptible")
            };
            let seq: Vec<u8> = if FOREIGN.contains(&id.as_str()) {
                random(&mut rng)
            } else {
                species_root
                    .iter()
                    .map(|&b| {
                        if rng.random_bool(0.005) {
                            b"ACGT"[rng.random_range(0..4)]
                        } else {
                            b
                        }
                    })
                    .collect()
            };
            let mut fasta = format!(">{id}\n").into_bytes();
            fasta.extend(&seq);
            fasta.push(b'\n');
            std::fs::write(grm::bvbrc::fasta_path(&root.join("fasta"), &id), fasta).unwrap();
            let record = serde_json::json!({
                "genome_id": id, "antibiotic": "meropenem", "resistant_phenotype": call,
                "laboratory_typing_method": "Broth dilution",
            });
            amr += &format!("{record}\n");
            let genome = serde_json::json!({
                "genome_id": id, "species": "Escherichia coli", "genome_length": "100000",
            });
            meta += &format!("{genome}\n");
        }
        std::fs::write(root.join("bvbrc/genome_amr.jsonl"), amr).unwrap();
        std::fs::write(root.join("bvbrc/genomes.jsonl"), meta).unwrap();
    }

    /// `curate` drops the failed genomes from genomes.txt and selects
    /// replacements. The check that follows must still judge the dropped
    /// genomes against the same reference, or the next `curate` lets them
    /// back in.
    #[test]
    fn species_check_and_curate_reach_a_fixed_point() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        species_check_fixture(dir.path());
        let folds = grm::learn::eval::EvalConfig::default().folds;
        let selection =
            || read_lines(&layout.species_dir("escherichia_coli").join("genomes.txt")).unwrap();
        let foreign: Vec<String> = FOREIGN.map(String::from).to_vec();
        let state = || {
            let check: serde_json::Value =
                serde_json::from_slice(&std::fs::read(layout.species_check()).unwrap()).unwrap();
            let failed: Vec<String> = check["genomes"]
                .as_array()
                .unwrap()
                .iter()
                .filter(|g| g["verdict"] == "fail")
                .map(|g| g["genome_id"].as_str().unwrap().to_string())
                .collect();
            let references: Vec<String> = check["species"]
                .as_array()
                .unwrap()
                .iter()
                .map(|s| s["reference"].as_str().unwrap().to_string())
                .collect();
            let (_, qc): (serde_json::Value, BTreeMap<String, String>) =
                serde_json::from_slice(&std::fs::read(layout.curated().join("qc.json")).unwrap())
                    .unwrap();
            let dropped: BTreeMap<String, String> = qc
                .into_iter()
                .filter(|(_, why)| why.starts_with(curate::SPECIES_CHECK_FAILED))
                .collect();
            (failed, references, dropped, selection())
        };

        run_curate(&layout, 50, folds, None).unwrap();
        let first = selection();
        assert!(foreign.iter().all(|id| first.contains(id)), "{first:?}");
        let mut passes = Vec::new();
        for _ in 0..3 {
            species_check(&layout, false, false).unwrap();
            run_curate(&layout, 50, folds, None).unwrap();
            passes.push(state());
        }
        for (failed, _, dropped, selected) in &passes {
            assert_eq!(failed, &foreign);
            assert_eq!(dropped.keys().cloned().collect::<Vec<_>>(), foreign);
            assert_eq!(selected.len(), 50);
            assert!(
                selected.iter().all(|id| !foreign.contains(id)),
                "{selected:?}"
            );
        }
        // The first pass judged the 50 first selected; later passes also
        // judge the three replacements, so the evidence text may move once.
        assert_eq!(passes[1].1, passes[0].1, "references moved");
        assert_eq!(passes[1].3, passes[0].3, "selection moved");
        assert_eq!(passes[2], passes[1]);
        let check = grm::species_check::read(&layout.species_check()).unwrap();
        assert_eq!(check.species[0].checked, 53);
    }

    #[test]
    fn pipeline_and_learn_take_a_k_and_a_retained_feature_budget() {
        assert_eq!(
            DEFAULT_SCREENED_FEATURES,
            grm::learn::eval::EvalConfig::default().screened_features
        );
        let parse = |args: &[&str]| Cli::try_parse_from([&["grm", "pipeline"][..], args].concat());
        let Command::Pipeline {
            k,
            screened_features,
            ..
        } = parse(&[]).unwrap().command
        else {
            panic!("parsed another command");
        };
        assert_eq!(
            (k, screened_features),
            (PIPELINE_K, DEFAULT_SCREENED_FEATURES)
        );
        let Command::Pipeline {
            k,
            screened_features,
            ..
        } = parse(&["--k", "21", "--screened-features", "16000"])
            .unwrap()
            .command
        else {
            panic!("parsed another command");
        };
        assert_eq!((k, screened_features), (21, 16000));
        for bad in [&["--k", "33"][..], &["--screened-features", "0"][..]] {
            assert!(parse(bad).is_err(), "{bad:?}");
        }
        let learn = Cli::try_parse_from([
            "grm",
            "learn",
            "--species",
            "s",
            "--screened-features",
            "1000",
        ]);
        let Command::Learn {
            screened_features,
            exp_tie_rule,
            exp_tie_trace,
            fold_file,
            fold_out,
            ..
        } = learn.unwrap().command
        else {
            panic!("parsed another command");
        };
        assert_eq!(screened_features, 1000);
        // The tie-rule experiment is off unless asked for.
        assert_eq!((exp_tie_rule, exp_tie_trace), (TieRuleArg::Float, None));
        let learn =
            Cli::try_parse_from(["grm", "learn", "--species", "s", "--exp-tie-rule", "exact"]);
        let Command::Learn { exp_tie_rule, .. } = learn.unwrap().command else {
            panic!("parsed another command");
        };
        assert_eq!(exp_tie_rule, TieRuleArg::Exact);
        // Without --fold-file the command is the one it was.
        assert_eq!((fold_file, fold_out), (None, None));
    }

    #[test]
    fn learn_takes_the_mic_experiment_only_with_an_output_and_no_other_experiment() {
        let parse = |rest: &[&str]| {
            Cli::try_parse_from(["grm", "learn", "--species", "s"].iter().chain(rest))
        };
        // Off unless asked for.
        let Command::Learn {
            exp_mic,
            exp_mic_out,
            ..
        } = parse(&[]).unwrap().command
        else {
            panic!("parsed another command");
        };
        assert_eq!((exp_mic, exp_mic_out), (vec![], None));
        let Command::Learn {
            exp_mic,
            exp_mic_out,
            seed,
            ..
        } = parse(&[
            "--exp-mic",
            "released,ordered",
            "--exp-mic-out",
            "o",
            "--seed",
            "43",
        ])
        .unwrap()
        .command
        else {
            panic!("parsed another command");
        };
        assert_eq!(exp_mic, vec![MicArmArg::Released, MicArmArg::Ordered]);
        assert_eq!((exp_mic_out, seed), (Some("o".into()), 43));
        assert!(
            parse(&[
                "--exp-mic",
                "ordered",
                "--fold-file",
                "f.json",
                "--fold-out",
                "o"
            ])
            .is_ok()
        );
        for bad in [
            &["--exp-mic", "median", "--exp-mic-out", "o"][..],
            &["--exp-mic-out", "o"],
            &[
                "--exp-mic",
                "released",
                "--exp-mic-out",
                "o",
                "--learner-seed",
                "1",
            ],
            &[
                "--exp-mic",
                "released",
                "--exp-mic-out",
                "o",
                "--exp-tie-rule",
                "exact",
            ],
            &[
                "--exp-mic",
                "released",
                "--exp-mic-out",
                "o",
                "--exp-tie-trace",
                "t",
            ],
            &[
                "--exp-mic",
                "released",
                "--exp-mic-out",
                "o",
                "--exp-tie-out",
                "t",
            ],
            &[
                "--exp-mic",
                "released",
                "--exp-mic-out",
                "o",
                "--exp-null-refit",
                "t",
            ],
            &[
                "--exp-mic",
                "released",
                "--fold-file",
                "f",
                "--fold-out",
                "o",
                "--exp-mic-out",
                "p",
            ],
            &[
                "--exp-mic",
                "released",
                "--fold-file",
                "f",
                "--fold-out",
                "o",
                "--seed",
                "43",
            ],
        ] {
            assert!(parse(bad).is_err(), "{bad:?}");
        }
        // Without an output it stops before reading the data root.
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path().join("data"));
        let species = ["s".to_string()];
        let error = learn_mic(
            &layout,
            &species,
            &[],
            false,
            10,
            42,
            None,
            None,
            &[MicArmArg::Released],
        )
        .unwrap_err();
        assert!(format!("{error:#}").contains("--exp-mic-out"));
        assert!(!dir.path().join("data").exists());
    }

    #[test]
    fn learn_takes_the_cohort_experiment_only_with_an_output_and_no_other_experiment() {
        let parse = |rest: &[&str]| {
            Cli::try_parse_from(["grm", "learn", "--species", "s"].iter().chain(rest))
        };
        let Command::Learn {
            exp_cohort,
            exp_cohort_out,
            ..
        } = parse(&[]).unwrap().command
        else {
            panic!("parsed another command");
        };
        assert_eq!((exp_cohort, exp_cohort_out), (None, None));
        let Command::Learn {
            exp_cohort,
            exp_cohort_out,
            seed,
            only,
            ..
        } = parse(&[
            "--exp-cohort",
            "c",
            "--exp-cohort-out",
            "o",
            "--seed",
            "45",
            "--only",
            "d",
        ])
        .unwrap()
        .command
        else {
            panic!("parsed another command");
        };
        assert_eq!(
            (exp_cohort, exp_cohort_out, seed, only),
            (
                Some("c".into()),
                Some("o".into()),
                45,
                vec!["d".to_string()]
            )
        );
        let with = |extra: &[&'static str]| {
            let mut v = vec!["--exp-cohort", "c", "--exp-cohort-out", "o"];
            v.extend_from_slice(extra);
            v
        };
        for bad in [
            vec!["--exp-cohort", "c"],
            vec!["--exp-cohort-out", "o"],
            with(&["--fold-file", "f", "--fold-out", "p"]),
            with(&["--learner-seed", "1"]),
            with(&["--exp-tie-rule", "exact"]),
            with(&["--exp-tie-trace", "t"]),
            with(&["--exp-tie-out", "t"]),
            with(&["--exp-null-refit", "t"]),
            with(&["--exp-mic", "released"]),
            with(&["--exp-mic", "released", "--exp-mic-out", "p"]),
        ] {
            assert!(parse(&bad).is_err(), "{bad:?}");
        }
    }

    #[test]
    fn learn_takes_a_fold_file_only_with_its_output_directory_and_one_species() {
        let parse = |rest: &[&str]| {
            Cli::try_parse_from(["grm", "learn", "--species", "s"].iter().chain(rest))
        };
        let Command::Learn {
            fold_file,
            fold_out,
            force,
            ..
        } = parse(&["--fold-file", "f.json", "--fold-out", "o"])
            .unwrap()
            .command
        else {
            panic!("parsed another command");
        };
        assert_eq!(
            (fold_file, fold_out, force),
            (Some("f.json".into()), Some("o".into()), false)
        );
        for bad in [
            &["--fold-file", "f.json"][..],
            &["--fold-out", "o"][..],
            &[
                "--fold-file",
                "f.json",
                "--fold-out",
                "o",
                "--only",
                "meropenem",
            ][..],
        ] {
            assert!(parse(bad).is_err(), "{bad:?}");
        }
        // A fold file runs at production's tie rule: the experiment's
        // options are refused beside it, `float` named explicitly too.
        for extra in [
            &["--exp-tie-rule", "exact"][..],
            &["--exp-tie-rule", "float"][..],
            &["--exp-tie-trace", "t.jsonl"][..],
            &["--exp-tie-out", "x"][..],
        ] {
            let args = [&["--fold-file", "f.json", "--fold-out", "o"][..], extra].concat();
            assert!(parse(&args).is_err(), "{extra:?}");
        }
        // The file names one dataset, so a second species is refused before
        // anything is read.
        let dir = tempfile::tempdir().unwrap();
        let error = learn_fold_file(
            &Layout::new(dir.path()),
            &["a".into(), "b".into()],
            &dir.path().join("f.json"),
            &dir.path().join("o"),
            false,
            DEFAULT_SCREENED_FEATURES,
        )
        .unwrap_err();
        assert_eq!(
            format!("{error:#}"),
            "--fold-file names one dataset, so it takes one species, not 2"
        );
        assert!(std::fs::read_dir(dir.path()).unwrap().next().is_none());
    }

    /// A fold file run stays at the release's seeds: the parser refuses
    /// `--seed` and `--learner-seed` beside `--fold-file`, even at the
    /// release's own value, and without them it parses the release's seeds.
    #[test]
    fn learn_refuses_a_fold_file_beside_either_seed_flag() {
        let parse = |rest: &[&str]| {
            let all = [
                &["grm", "learn", "--species", "s"][..],
                &["--fold-file", "f.json", "--fold-out", "o"][..],
                rest,
            ]
            .concat();
            Cli::try_parse_from(all)
        };
        for bad in [
            &["--seed", "43"][..],
            &["--seed", "42"][..],
            &["--learner-seed", "7"][..],
            &["--seed", "42", "--learner-seed", "7"][..],
        ] {
            let error = parse(bad).err().unwrap_or_else(|| panic!("{bad:?} parsed"));
            assert_eq!(
                error.kind(),
                clap::error::ErrorKind::ArgumentConflict,
                "{bad:?}"
            );
        }
        let Command::Learn {
            fold_file,
            seed,
            learner_seed,
            ..
        } = parse(&[]).unwrap().command
        else {
            panic!("parsed another command");
        };
        assert_eq!(
            (fold_file, seed, learner_seed),
            (Some("f.json".into()), RELEASE_SEED, None)
        );
    }

    /// The null-fold experiment records the folds `grm learn` draws, so the
    /// parser refuses `--exp-null-refit` beside `--fold-file`, before FILE is
    /// opened, in either order. Without `--fold-file` the flag parses as it
    /// did.
    #[test]
    fn learn_refuses_the_null_fold_refit_beside_a_fold_file() {
        let parse = |args: &[&str]| {
            let all = [&["grm", "learn", "--species", "s"][..], args].concat();
            Cli::try_parse_from(all)
        };
        let refit = ["--exp-null-refit", "r.jsonl"];
        let fold = ["--fold-file", "f.json", "--fold-out", "o"];
        for bad in [
            [&refit[..], &fold[..]].concat(),
            [&fold[..], &refit[..]].concat(),
        ] {
            let error = parse(&bad)
                .err()
                .unwrap_or_else(|| panic!("{bad:?} parsed"));
            assert_eq!(
                error.kind(),
                clap::error::ErrorKind::ArgumentConflict,
                "{bad:?}"
            );
            let message = error.to_string();
            assert!(
                message.contains("--exp-null-refit") && message.contains("--fold-file"),
                "{message}"
            );
        }
        let Command::Learn {
            fold_file,
            exp_null_refit,
            ..
        } = parse(&refit).unwrap().command
        else {
            panic!("parsed another command");
        };
        assert_eq!((fold_file, exp_null_refit), (None, Some("r.jsonl".into())));
    }

    /// The null-fold experiment fits its fallbacks under production's tie
    /// rule, and those fits would add lines to a tie trace, so the parser
    /// refuses `--exp-null-refit` beside each option of the tie rule
    /// experiment, in either order, `float` named explicitly too.
    #[test]
    fn learn_refuses_the_null_fold_refit_beside_the_tie_rule_experiment() {
        let parse = |args: &[&str]| {
            let all = [&["grm", "learn", "--species", "s"][..], args].concat();
            Cli::try_parse_from(all)
        };
        let refit = ["--exp-null-refit", "r.jsonl"];
        for tie in [
            &["--exp-tie-rule", "exact"][..],
            &["--exp-tie-rule", "float"][..],
            &["--exp-tie-trace", "t.jsonl"][..],
            &["--exp-tie-out", "x"][..],
        ] {
            for bad in [[&refit[..], tie].concat(), [tie, &refit[..]].concat()] {
                let error = parse(&bad)
                    .err()
                    .unwrap_or_else(|| panic!("{bad:?} parsed"));
                assert_eq!(
                    error.kind(),
                    clap::error::ErrorKind::ArgumentConflict,
                    "{bad:?}"
                );
                assert!(error.to_string().contains("--exp-null-refit"), "{bad:?}");
            }
        }
    }

    /// `grm learn` draws at the release's seed unless told otherwise, and
    /// `grm pipeline`, the release's command, takes no seed at all.
    #[test]
    fn learn_takes_a_seed_and_a_learner_seed_and_the_pipeline_takes_neither() {
        assert_eq!(RELEASE_SEED, grm::learn::eval::EvalConfig::default().seed);
        let parse = |args: &[&str]| {
            let all = [&["grm", "learn", "--species", "s"][..], args].concat();
            let Command::Learn {
                seed, learner_seed, ..
            } = Cli::try_parse_from(all)?.command
            else {
                panic!("parsed another command");
            };
            Ok::<_, clap::Error>((seed, learner_seed))
        };
        assert_eq!(parse(&[]).unwrap(), (42, None));
        assert_eq!(
            parse(&["--seed", "10622993013839509822"]).unwrap(),
            (10622993013839509822, None)
        );
        assert_eq!(
            parse(&["--learner-seed", "7"]).unwrap(),
            (RELEASE_SEED, Some(7))
        );
        assert!(parse(&["--seed", "-1"]).is_err());
        assert!(parse(&["--seed"]).is_err());
        for flag in ["--seed", "--learner-seed"] {
            assert!(
                Cli::try_parse_from(["grm", "pipeline", flag, "43"]).is_err(),
                "{flag}"
            );
        }
    }

    /// At the release's seeds `grm learn` ends with `main`'s line, byte for
    /// byte, whatever the budget; a draw's line names its folder and does not
    /// begin with `<slug>: <n> dataset reports written`, the release line
    /// `scripts/metric-audit/audit.py` reads (its suite tests the pattern
    /// on the draw lines below).
    #[test]
    fn learn_ends_with_the_release_line_and_a_draw_with_a_line_of_its_own() {
        let layout = Layout::new("root");
        let line = |screened_features, seed, learner_seed| {
            let cfg = learn_config(screened_features, seed, learner_seed, Default::default());
            learn_end_line(&layout, "fixture_species", 3, &cfg)
        };
        let release = "fixture_species: 3 dataset reports written";
        assert_eq!(line(DEFAULT_SCREENED_FEATURES, RELEASE_SEED, None), release);
        assert_eq!(line(16000, RELEASE_SEED, None), release);
        assert_eq!(
            line(DEFAULT_SCREENED_FEATURES, 43, None),
            "fixture_species: 3 draw reports written to root/species/fixture_species/draws/seed-43"
        );
        assert_eq!(
            line(DEFAULT_SCREENED_FEATURES, RELEASE_SEED, Some(7)),
            "fixture_species: 3 draw reports written to root/species/fixture_species/draws/seed-42-learner-7"
        );
    }

    /// The evaluation `grm learn` runs, built from its parsed flags as the
    /// command builds it: the release's at the defaults, and otherwise the
    /// release's with the budget and both seeds the flags give.
    #[test]
    fn learn_evaluates_with_the_budget_and_seeds_its_flags_give() {
        use grm::learn::eval::EvalConfig;
        let config = |args: &[&str]| {
            let all = [&["grm", "learn", "--species", "s"][..], args].concat();
            let Command::Learn {
                screened_features,
                seed,
                learner_seed,
                exp_tie_rule,
                ..
            } = Cli::try_parse_from(all).unwrap().command
            else {
                panic!("parsed another command");
            };
            serde_json::to_value(learn_config(
                screened_features,
                seed,
                learner_seed,
                exp_tie_rule.into(),
            ))
            .unwrap()
        };
        let release = || EvalConfig::default();
        let expected = |cfg: EvalConfig| serde_json::to_value(cfg).unwrap();
        assert_eq!(config(&[]), expected(release()));
        assert_eq!(
            config(&["--seed", "43"]),
            expected(EvalConfig {
                seed: 43,
                ..release()
            })
        );
        assert_eq!(
            config(&["--learner-seed", "7"]),
            expected(EvalConfig {
                learner_seed: Some(7),
                ..release()
            })
        );
        assert_eq!(
            config(&[
                "--seed",
                "18446744073709551615",
                "--learner-seed",
                "0",
                "--screened-features",
                "1000",
            ]),
            expected(EvalConfig {
                seed: u64::MAX,
                learner_seed: Some(0),
                screened_features: 1000,
                ..release()
            })
        );
        assert_eq!(
            config(&["--exp-tie-rule", "exact", "--seed", "43"]),
            expected(EvalConfig {
                seed: 43,
                tie_rule: grm::learn::tie::TieRule::Exact,
                ..release()
            })
        );
        assert_eq!(config(&["--exp-tie-rule", "float"]), expected(release()));
    }

    /// The pipeline's stages follow `--k`: counts, structure and matrix built
    /// at 31 are rebuilt at 21, and each records 21.
    #[test]
    fn the_pipeline_counts_and_builds_at_the_k_it_is_given() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        count_to_matrix(&layout, SLUG, 21).unwrap();
        assert_eq!(learn_refusal(&layout), None);
        let (input, _) = grm::analysis::load_inputs(&layout, SLUG).unwrap();
        assert_eq!((input.structure.k, input.matrix.summary.k), (21, 21));
        let counted = grm::species::load_counted(&layout, SLUG).unwrap();
        assert_eq!(counted.k, 21);
    }

    /// `grm pipeline --k 21` from its parsed arguments through `run`, on the
    /// stock path and with `--disk-sparse`: counts, structure and matrix
    /// each record 21. A dispatch that passed another k on (the merge's
    /// `src/main.rs` hunk) would count at 31. Synthetic fixture.
    #[test]
    fn the_pipeline_command_counts_at_the_k_it_is_given() {
        for sparse in [false, true] {
            let dir = tempfile::tempdir().unwrap();
            let layout = fixture(dir.path());
            let data = dir.path().to_str().unwrap();
            let mut args = vec![
                "grm",
                "--data",
                data,
                "pipeline",
                "--species",
                SLUG,
                "--k",
                "21",
            ];
            if sparse {
                args.extend(["--disk-sparse", "--block-width", "64"]);
            }
            let ran = run(Cli::try_parse_from(args).unwrap()).map_err(|e| format!("{e:#}"));
            let counted = grm::species::load_counted(&layout, SLUG).unwrap();
            assert_eq!(counted.k, 21, "sparse {sparse}: {ran:?}");
            let (input, _) = grm::analysis::load_inputs(&layout, SLUG).unwrap();
            assert_eq!(
                (input.structure.k, input.matrix.summary.k),
                (21, 21),
                "sparse {sparse}: {ran:?}"
            );
            // The fixture has no curated/funnel.json, so the run ends at
            // the report step, after count, structure, matrix and learn.
            assert!(!layout.curated().join("funnel.json").exists());
            assert_eq!(
                ran,
                Err("No such file or directory (os error 2)".to_string()),
                "sparse {sparse}"
            );
        }
    }

    /// `grm pipeline` with `species/<slug>/results` linked outside the data
    /// root: the stock run makes and writes it through the link as `main`
    /// does and goes on to the report step, while `--disk-sparse` stops at
    /// learn with nothing written behind the link. Synthetic fixture.
    #[test]
    fn only_the_disk_sparse_pipeline_refuses_a_linked_results_directory() {
        for sparse in [false, true] {
            let dir = tempfile::tempdir().unwrap();
            let layout = fixture(dir.path());
            let outside = tempfile::tempdir().unwrap();
            let results = layout.results_dir(SLUG);
            ensure_parent(&results).unwrap();
            std::os::unix::fs::symlink(outside.path(), &results).unwrap();
            let data = dir.path().to_str().unwrap();
            let mut args = vec![
                "grm",
                "--data",
                data,
                "pipeline",
                "--species",
                SLUG,
                "--k",
                "21",
            ];
            if sparse {
                args.extend(["--disk-sparse", "--block-width", "64"]);
            }
            let ran = run(Cli::try_parse_from(args).unwrap()).map_err(|e| format!("{e:#}"));
            if sparse {
                let err = ran.unwrap_err();
                assert!(err.contains("symlinked directory"), "{err}");
            } else {
                assert_eq!(
                    ran,
                    Err("No such file or directory (os error 2)".to_string())
                );
            }
            assert!(results.is_symlink());
            assert_eq!(std::fs::read_dir(outside.path()).unwrap().count(), 0);
        }
    }

    fn curated_json(layout: &Layout, file: &str) -> serde_json::Value {
        serde_json::from_slice(&std::fs::read(layout.curated().join(file)).unwrap()).unwrap()
    }

    #[test]
    fn curate_records_the_legacy_rule_and_its_cap_without_changing_its_choice() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        species_check_fixture(dir.path());
        let folds = grm::learn::eval::EvalConfig::default().folds;
        run_curate(&layout, 40, folds, None).unwrap();
        let selected =
            read_lines(&layout.species_dir("escherichia_coli").join("genomes.txt")).unwrap();
        let phenotypes: Vec<curate::Phenotype> =
            bvbrc::read_jsonl(&layout.curated().join("phenotypes.jsonl")).unwrap();
        let datasets: Vec<curate::DatasetSpec> =
            serde_json::from_value(curated_json(&layout, "datasets.json")).unwrap();
        assert_eq!(
            selected,
            curate::select_genomes(&phenotypes, &datasets, "Escherichia coli", 40)
        );
        let want = serde_json::json!({
            "rule": "legacy", "cap": 40, "seed": null, "always_included": null,
        });
        assert_eq!(curated_json(&layout, "qc.json")[0]["selection"], want);
        assert_eq!(curated_json(&layout, "qc.json")[0]["min_class"], folds);
        for d in curated_json(&layout, "datasets.json").as_array().unwrap() {
            assert_eq!(d["selection"], want);
        }
        assert!(
            !layout
                .species_dir("escherichia_coli")
                .join("selection.json")
                .exists()
        );
    }

    #[test]
    fn curate_draws_round_robin_across_mlst_groups_and_records_the_plan() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        species_check_fixture(dir.path());
        let folds = grm::learn::eval::EvalConfig::default().folds;
        // Give the genome records a BV-BRC MLST field: six sequence types, and
        // none on every fifth genome, which then is its own singleton group.
        let records: Vec<serde_json::Value> = bvbrc::read_jsonl(&layout.genomes()).unwrap();
        let mut text = String::new();
        let mut without = BTreeSet::new();
        for (i, mut g) in records.into_iter().enumerate() {
            if i % 5 == 0 {
                without.insert(g["genome_id"].as_str().unwrap().to_string());
            } else {
                g["mlst"] = format!("MLST.Escherichia_coli_1.{}", i % 6).into();
            }
            text += &format!("{g}\n");
        }
        std::fs::write(layout.genomes(), text).unwrap();
        let stratum = dir.path().join("stratum.json");
        let pinned = without.iter().next().unwrap().clone();
        std::fs::write(
            &stratum,
            serde_json::json!({"genome_ids": [pinned], "reason": "fixture", "label_dependent": false})
                .to_string(),
        )
        .unwrap();
        let rule =
            selection_rule(SelectionRule::GroupRoundRobinV1, Some(7), Some(stratum)).unwrap();
        run_curate_with(
            &layout,
            30,
            folds,
            None,
            &rule,
            curate::LabelRule::default(),
        )
        .unwrap();
        let species_dir = layout.species_dir("escherichia_coli");
        let plan: curate::GroupSamplingPlan =
            serde_json::from_slice(&std::fs::read(species_dir.join("selection.json")).unwrap())
                .unwrap();
        let selected = read_lines(&species_dir.join("genomes.txt")).unwrap();
        assert_eq!(selected, plan.selected);
        assert_eq!(selected.len(), 30);
        assert_eq!(selected[0], pinned);
        // Every eligible genome has a group: its sequence type, or itself.
        assert_eq!(plan.inclusion.len(), 60);
        for g in &plan.inclusion {
            let want = if without.contains(&g.genome_id) {
                format!("singleton:{}", g.genome_id)
            } else {
                assert!(g.group_id.starts_with("mlst:MLST.Escherichia_coli_1."));
                g.group_id.clone()
            };
            assert_eq!(g.group_id, want);
        }
        let qc = &curated_json(&layout, "qc.json")[0]["selection"];
        assert_eq!(qc["rule"], curate::GROUP_SAMPLING_DESIGN);
        assert_eq!(qc["group_source"], curate::MLST_GROUP_SOURCE);
        assert_eq!(
            (qc["cap"].clone(), qc["seed"].clone()),
            (30.into(), 7.into())
        );
        assert!(qc.get("input_sha256").is_none() && qc.get("singletons").is_none());
        for d in curated_json(&layout, "datasets.json").as_array().unwrap() {
            let s = &d["selection"];
            assert_eq!(s["input_sha256"], plan.input_sha256.as_str());
            assert_eq!(s["group_source"], curate::MLST_GROUP_SOURCE);
            assert_eq!(s["singletons"], without.len());
            assert_eq!(s["always_included"]["genome_ids"][0], pinned.as_str());
        }
        // A legacy run afterwards drops the stale plan and names no group source.
        run_curate(&layout, 30, folds, None).unwrap();
        assert!(!species_dir.join("selection.json").exists());
        let qc = &curated_json(&layout, "qc.json")[0]["selection"];
        assert!(qc.get("group_source").is_none());
    }

    #[test]
    fn selection_flags_are_refused_where_they_do_not_apply() {
        let legacy = |seed, stratum| selection_rule(SelectionRule::Legacy, seed, stratum);
        assert!(legacy(Some(1), None).is_err());
        assert!(legacy(None, Some(PathBuf::from("s.json"))).is_err());
        let err = selection_rule(SelectionRule::GroupRoundRobinV1, None, None)
            .err()
            .unwrap()
            .to_string();
        assert!(err.contains("--selection-seed"), "{err}");
        assert!(selection_rule(SelectionRule::GroupRoundRobinV1, Some(1), None).is_ok());
        let parsed = Cli::try_parse_from(["grm", "curate", "--selection", "group-round-robin-v1"]);
        assert!(parsed.is_ok());
        assert!(Cli::try_parse_from(["grm", "curate", "--selection", "random"]).is_err());
        assert!(Cli::try_parse_from(["grm", "curate", "--groups", "g.tsv"]).is_err());
    }

    /// `--reset` drops the recorded references and verdicts, not the
    /// recorded genomes. Curation has removed the failed ones from
    /// genomes.txt, so a reset check that forgot them would let the next
    /// curate select them again.
    #[test]
    fn a_reset_check_still_judges_the_genomes_curation_dropped() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        species_check_fixture(dir.path());
        let folds = grm::learn::eval::EvalConfig::default().folds;
        run_curate(&layout, 50, folds, None).unwrap();
        species_check(&layout, false, false).unwrap();
        run_curate(&layout, 50, folds, None).unwrap();
        species_check(&layout, false, true).unwrap();
        run_curate(&layout, 50, folds, None).unwrap();
        let check = grm::species_check::read(&layout.species_check()).unwrap();
        let failures = check.failures();
        let failed: Vec<&str> = failures.keys().map(String::as_str).collect();
        assert_eq!(failed, FOREIGN);
        assert_eq!(check.species[0].checked, 53);
        let selected =
            read_lines(&layout.species_dir("escherichia_coli").join("genomes.txt")).unwrap();
        assert!(
            selected.iter().all(|id| !FOREIGN.contains(&id.as_str())),
            "{selected:?}"
        );
    }

    /// A reset check cannot judge a recorded failure with no sketch and no
    /// FASTA file, so it stops instead of dropping it from the record. A
    /// check without --reset keeps its recorded verdict.
    #[test]
    fn a_reset_check_stops_on_a_recorded_failure_it_cannot_judge() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        species_check_fixture(dir.path());
        let folds = grm::learn::eval::EvalConfig::default().folds;
        run_curate(&layout, 50, folds, None).unwrap();
        species_check(&layout, false, false).unwrap();
        run_curate(&layout, 50, folds, None).unwrap();
        let gone = FOREIGN[1];
        std::fs::remove_file(grm::species_check::sketch_path(&layout, gone)).unwrap();
        std::fs::remove_file(grm::bvbrc::fasta_path(&layout.fasta_dir(), gone)).unwrap();

        let err = species_check(&layout, false, true).unwrap_err().to_string();
        assert!(
            err.contains(&format!("cannot judge {gone} (Escherichia coli)")),
            "{err}"
        );
        species_check(&layout, false, false).unwrap();
        let check = grm::species_check::read(&layout.species_check()).unwrap();
        let failures = check.failures();
        let failed: Vec<&str> = failures.keys().map(String::as_str).collect();
        assert_eq!(failed, FOREIGN);
    }

    #[test]
    fn disk_sparse_matrix_command_counts_a_fresh_cohort_without_learning() {
        let dir = tempfile::tempdir().unwrap();
        let stock = fixture(&dir.path().join("stock"));
        let sparse = fixture(&dir.path().join("sparse"));
        count_to_matrix(&stock, SLUG, PIPELINE_K).unwrap();
        let data = sparse.root.to_str().unwrap();
        let cli = Cli::try_parse_from([
            "grm",
            "--data",
            data,
            "matrix",
            "--species",
            SLUG,
            "--min-support",
            "2",
            "--disk-sparse",
            "--block-width",
            "64",
        ])
        .unwrap();
        run(cli).unwrap();
        for name in [
            "matrix.json",
            "patterns.bin",
            "support.bin",
            "kmer_offsets.bin",
            "kmers.bin",
        ] {
            assert_eq!(
                std::fs::read(sparse.matrix_dir(SLUG).join(name)).unwrap(),
                std::fs::read(stock.matrix_dir(SLUG).join(name)).unwrap(),
                "{name}"
            );
        }
        let counted = grm::species::load_counted(&sparse, SLUG).unwrap();
        assert!(counted.pending.is_empty());
        assert!(counted.failed.is_empty());
        assert!(!sparse.kset_dir(SLUG).join("g0.kset").exists());
        assert!(!sparse.species_dir(SLUG).join("structure.json").exists());
        assert!(!sparse.results_dir(SLUG).exists());
        assert!(
            Cli::try_parse_from([
                "grm",
                "matrix",
                "--species",
                SLUG,
                "--disk-sparse",
                "--block-width",
                "65",
            ])
            .is_err()
        );
    }

    #[test]
    fn disk_sparse_matrix_passes_the_catalog_gate_before_counting_or_deleting() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        let data = layout.root.to_str().unwrap().to_string();
        let parse = |extra: &[&str]| {
            let mut args = vec![
                "grm",
                "--data",
                &data,
                "matrix",
                "--species",
                SLUG,
                "--disk-sparse",
                "--block-width",
                "64",
            ];
            args.extend_from_slice(extra);
            Cli::try_parse_from(args).unwrap()
        };
        let species = layout.species_dir(SLUG);
        let entries = || {
            let mut names: Vec<String> = std::fs::read_dir(&species)
                .unwrap()
                .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
                .collect();
            names.sort();
            names
        };
        // The catalog stops loading, as when its link target moves. Stock
        // `grm count` refuses; so does the sparse matrix, before it writes.
        std::fs::remove_file(layout.catalog_dir().join("AMR_CDS.fa")).unwrap();
        let count = Cli::try_parse_from(["grm", "--data", &data, "count", "--species", SLUG]);
        assert!(run(count.unwrap()).is_err());
        let err = run(parse(&[])).unwrap_err();
        assert!(
            format!("{err:#}").contains("--allow-missing-catalog"),
            "{err:#}"
        );
        assert_eq!(entries(), ["genomes.txt"]);
        // As with `grm count --allow-missing-catalog`, the flag counts without
        // gene calls.
        run(parse(&["--allow-missing-catalog"])).unwrap();
        let counted = grm::species::load_counted(&layout, SLUG).unwrap();
        assert!(counted.genomes.iter().all(|g| g.catalog_version.is_empty()));
        assert!(
            Cli::try_parse_from([
                "grm",
                "matrix",
                "--species",
                SLUG,
                "--allow-missing-catalog"
            ])
            .is_err()
        );

        // Ksets counted with the catalog are never deleted for a count
        // without it, flag or no flag.
        let other = tempfile::tempdir().unwrap();
        let layout = fixture(other.path());
        count_to_matrix(&layout, SLUG, PIPELINE_K).unwrap();
        let ksets = layout.kset_dir(SLUG);
        let kept = std::fs::read_dir(&ksets).unwrap().count();
        let before = std::fs::read(layout.species_dir(SLUG).join("counted.json")).unwrap();
        std::fs::remove_file(layout.catalog_dir().join("AMR_CDS.fa")).unwrap();
        let data = layout.root.to_str().unwrap();
        for extra in [&[][..], &["--allow-missing-catalog"]] {
            let mut args = vec![
                "grm",
                "--data",
                data,
                "matrix",
                "--species",
                SLUG,
                "--disk-sparse",
                "--block-width",
                "64",
            ];
            args.extend_from_slice(extra);
            assert!(run(Cli::try_parse_from(args).unwrap()).is_err());
            assert_eq!(std::fs::read_dir(&ksets).unwrap().count(), kept);
            assert_eq!(
                std::fs::read(layout.species_dir(SLUG).join("counted.json")).unwrap(),
                before
            );
            assert!(!layout.species_dir(SLUG).join("matrix-partials").exists());
        }
    }

    #[test]
    fn curate_checks_the_class_minimum_before_reading_input() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let folds = grm::learn::eval::EvalConfig::default().folds;
        // The data root is empty, so only a guard that runs first can name
        // the class minimum.
        let below = run_curate(&layout, 2000, folds - 1, None).unwrap_err();
        let floor = format!("--min-class {} is below {folds}:", folds - 1);
        assert!(below.to_string().starts_with(&floor), "{below:#}");
        let at = run_curate(&layout, 2000, folds, None).unwrap_err();
        assert!(
            format!("{at:#}").contains("bvbrc/genome_amr.jsonl"),
            "{at:#}"
        );
    }

    /// Sixty genomes on a shared core for one fixture dataset. The 30
    /// resistant genomes carry block A and lack block B; of the susceptible
    /// ones, 15 carry both and 15 neither. Every odd genome carries A
    /// reverse-complemented, and eight 100-base blocks come and go at random.
    /// Annotation reads a cached empty feature list, so nothing is fetched.
    fn replay_fixture(root: &Path) -> (Layout, Vec<String>) {
        let layout = fixture(root);
        let mut rng = ChaCha8Rng::seed_from_u64(29);
        let mut random = |n: usize| -> String {
            (0..n)
                .map(|_| b"ACGT"[rng.random_range(0..4)] as char)
                .collect()
        };
        let (core, a, b) = (random(3000), random(200), random(200));
        let noise: Vec<String> = (0..8).map(|_| random(100)).collect();
        let a_reverse: String = a
            .bytes()
            .rev()
            .map(|c| match c {
                b'A' => 'T',
                b'C' => 'G',
                b'G' => 'C',
                _ => 'A',
            })
            .collect();
        let mut draw = ChaCha8Rng::seed_from_u64(31);
        let ids: Vec<String> = (0..60).map(|i| format!("g{i}")).collect();
        let mut phenotypes = Vec::new();
        for (i, id) in ids.iter().enumerate() {
            let mut seq = core.clone();
            let (carries_a, carries_b) = (i < 45, (30..45).contains(&i));
            if carries_a {
                seq += "N";
                seq += if i % 2 == 1 { &a_reverse } else { &a };
            }
            if carries_b {
                seq += "N";
                seq += &b;
            }
            for block in &noise {
                if draw.random_bool(0.5) {
                    seq += "N";
                    seq += block;
                }
            }
            let path = bvbrc::fasta_path(&layout.fasta_dir(), id);
            std::fs::write(path, format!(">c\n{seq}\n")).unwrap();
            let features = layout.features_dir().join(format!("{id}.json"));
            ensure_parent(&features).unwrap();
            std::fs::write(features, "[]").unwrap();
            let call = if i < 30 {
                curate::Call::R
            } else {
                curate::Call::S
            };
            phenotypes.push(curate::Phenotype {
                genome_id: id.clone(),
                species: "Fixture species".into(),
                antibiotic: "fixturemycin".into(),
                reported: Some(call),
                conflict: false,
                mic: None,
                mic_spread: 0,
                label: Some(call),
                label_from_mic: false,
                methods: vec![],
                standards: vec![],
                n_records: 1,
            });
        }
        let names: Vec<&str> = ids.iter().map(String::as_str).collect();
        write_genome_list(&layout, &names);
        let spec = curate::DatasetSpec {
            id: format!("{SLUG}__fixturemycin"),
            species: "Fixture species".into(),
            species_slug: SLUG.into(),
            antibiotic: "fixturemycin".into(),
            drug_class: "fixture class".into(),
            n_resistant: 30,
            n_susceptible: 30,
            n_intermediate: 0,
            n_with_mic: 0,
            n_labels_from_mic: 0,
            breakpoint: None,
            min_class: 25,
            selection: None,
        };
        std::fs::write(
            layout.curated().join("datasets.json"),
            serde_json::to_vec(&[spec]).unwrap(),
        )
        .unwrap();
        let lines: Vec<String> = phenotypes
            .iter()
            .map(|p| serde_json::to_string(p).unwrap())
            .collect();
        std::fs::write(
            layout.curated().join("phenotypes.jsonl"),
            lines.join("\n") + "\n",
        )
        .unwrap();
        (layout, ids)
    }

    /// `grm learn --fold-file` passes `--screened-features` on: the result
    /// records the value the command was given. The fold file repeats the
    /// random split the fixture's own report drew, so it is a fair partition.
    #[test]
    fn learn_with_a_fold_file_screens_as_many_features_as_it_is_told() {
        let dir = tempfile::tempdir().unwrap();
        let (layout, _) = replay_fixture(dir.path());
        let data = dir.path().to_str().unwrap();
        let grm = |args: &[&str]| {
            let all = [&["grm", "--data", data][..], args].concat();
            run(Cli::try_parse_from(all).unwrap())
        };
        grm(&["count", "--species", SLUG]).unwrap();
        grm(&["structure", "--species", SLUG]).unwrap();
        grm(&["matrix", "--species", SLUG]).unwrap();
        grm(&["learn", "--species", SLUG]).unwrap();
        let report: serde_json::Value = serde_json::from_slice(
            &std::fs::read(
                layout
                    .results_dir(SLUG)
                    .join(format!("{SLUG}__fixturemycin.json")),
            )
            .unwrap(),
        )
        .unwrap();
        let random = &report["splits"][0];
        assert_eq!(random["scheme"], "random");
        let mut test = vec![Vec::new(); random["folds"].as_u64().unwrap() as usize];
        for (sample, fold) in report["samples"]
            .as_array()
            .unwrap()
            .iter()
            .zip(random["test_fold"].as_array().unwrap())
        {
            test[fold.as_u64().unwrap() as usize].push(sample["id"].clone());
        }
        let file = dir.path().join("folds.json");
        let text = serde_json::json!({
            "format": "grm-fold-file-v1",
            "dataset": format!("{SLUG}__fixturemycin"),
            "scheme": "own-random",
            "test": test,
            "train": [],
            "unused": [],
        });
        std::fs::write(&file, serde_json::to_vec(&text).unwrap()).unwrap();
        let out = dir.path().join("fold-out");
        grm(&[
            "learn",
            "--species",
            SLUG,
            "--fold-file",
            file.to_str().unwrap(),
            "--fold-out",
            out.to_str().unwrap(),
            "--screened-features",
            "7",
        ])
        .unwrap();
        let result: serde_json::Value = serde_json::from_slice(
            &std::fs::read(out.join(format!("{SLUG}__fixturemycin.own-random.json"))).unwrap(),
        )
        .unwrap();
        assert_eq!(result["eval_config"]["screened_features"], 7);
        assert_eq!(result["counts"]["tested"], 60);
    }

    /// `grm learn --exp-tie-rule exact` through the command's own path. At
    /// the release's seeds it would write to results/, so it stops before
    /// any fit and names --exp-tie-out; naming results/ itself stops it too.
    /// With another folder the report lands there and records the rule the
    /// fits ran under, and results/ is not written. A run without the flag
    /// still writes results/ with no tie rule recorded.
    #[test]
    fn an_exact_learn_refuses_results_and_writes_the_rule_to_the_folder_it_names() {
        let dir = tempfile::tempdir().unwrap();
        let (layout, _) = replay_fixture(dir.path());
        let data = dir.path().to_str().unwrap();
        let grm = |args: &[&str]| {
            let all = [&["grm", "--data", data][..], args].concat();
            run(Cli::try_parse_from(all).unwrap()).map_err(|e| format!("{e:#}"))
        };
        grm(&["count", "--species", SLUG]).unwrap();
        grm(&["structure", "--species", SLUG]).unwrap();
        grm(&["matrix", "--species", SLUG]).unwrap();
        let file = format!("{SLUG}__fixturemycin.json");
        let results = layout.results_dir(SLUG);
        let exact = ["learn", "--species", SLUG, "--exp-tie-rule", "exact"];
        let err = grm(&exact).unwrap_err();
        assert_eq!(
            err,
            format!(
                "--exp-tie-rule exact would write {SLUG}'s reports to {}, the folder production reports are written to and `grm report` indexes; name another folder with --exp-tie-out DIR",
                results.display()
            )
        );
        assert!(!results.exists());
        // results/ is not made yet, and another spelling of it is refused all
        // the same: through a folder that is missing, and through a link.
        let species = layout.species_dir(SLUG);
        let mut spellings = vec![species.join("missing").join("..").join("results")];
        #[cfg(unix)]
        {
            let link = dir.path().join("link");
            std::os::unix::fs::symlink(&species, &link).unwrap();
            spellings.push(link.join("results"));
        }
        for spelling in spellings {
            let named = [&exact[..], &["--exp-tie-out", spelling.to_str().unwrap()]].concat();
            assert!(grm(&named).unwrap_err().contains("--exp-tie-out DIR"));
            assert!(!results.exists() && !species.join("missing").exists());
        }
        std::fs::create_dir_all(&results).unwrap();
        let named = [&exact[..], &["--exp-tie-out", results.to_str().unwrap()]].concat();
        assert!(grm(&named).unwrap_err().contains("--exp-tie-out DIR"));
        let out = dir.path().join("exact");
        let named = [&exact[..], &["--exp-tie-out", out.to_str().unwrap()]].concat();
        grm(&named).unwrap();
        let read = |p: std::path::PathBuf| -> serde_json::Value {
            serde_json::from_slice(&std::fs::read(p).unwrap()).unwrap()
        };
        assert_eq!(read(out.join(&file))["eval_config"]["tie_rule"], "exact");
        assert!(!results.join(&file).exists());
        grm(&["learn", "--species", SLUG]).unwrap();
        let production = read(results.join(&file));
        assert!(production["eval_config"].get("tie_rule").is_none());
        assert_ne!(
            production["input_fingerprint"],
            read(out.join(&file))["input_fingerprint"]
        );
    }

    /// `--exp-tie-out` may name a folder outside the data root.
    #[test]
    fn an_exact_learn_writes_to_a_folder_outside_the_data_root() {
        let dir = tempfile::tempdir().unwrap();
        let elsewhere = tempfile::tempdir().unwrap();
        replay_fixture(dir.path());
        let data = dir.path().to_str().unwrap();
        let grm = |args: &[&str]| {
            let all = [&["grm", "--data", data][..], args].concat();
            run(Cli::try_parse_from(all).unwrap()).map_err(|e| format!("{e:#}"))
        };
        grm(&["count", "--species", SLUG]).unwrap();
        grm(&["structure", "--species", SLUG]).unwrap();
        grm(&["matrix", "--species", SLUG]).unwrap();
        let out = elsewhere.path().join("exact");
        let exact = ["learn", "--species", SLUG, "--exp-tie-rule", "exact"];
        grm(&[&exact[..], &["--exp-tie-out", out.to_str().unwrap()]].concat()).unwrap();
        assert!(out.join(format!("{SLUG}__fixturemycin.json")).exists());
    }

    /// From the review of pull request 128: with a relative --data and
    /// results/ not yet made, an absolute --exp-tie-out naming results/
    /// passed the guard, since the two spellings differ and neither folder
    /// could be resolved. The guard now compares the folders themselves.
    #[test]
    fn an_exact_learn_refuses_results_under_another_spelling_before_it_exists() {
        let layout = Layout::new("no-such-data-root");
        let cfg = grm::learn::eval::EvalConfig {
            tie_rule: grm::learn::tie::TieRule::Exact,
            seed: 7,
            ..Default::default()
        };
        let results = std::env::current_dir()
            .unwrap()
            .join(layout.results_dir(SLUG));
        assert!(!results.exists());
        let refused = |out: &Path| refuse_exact_in_results(&layout, SLUG, &cfg, Some(out)).is_err();
        assert!(refused(&results));
        assert!(refused(&results.join("..").join("results")));
        assert!(refused(&layout.results_dir(SLUG)));
        assert!(!refused(&results.with_file_name("exact")));
        // At a seed of its own the run writes to draws/, which the guard
        // lets through with no folder named.
        refuse_exact_in_results(&layout, SLUG, &cfg, None).unwrap();
    }

    /// The acceptance check of docs/prediction-contract.md: `grm predict`,
    /// from the report `grm learn` wrote and each cohort genome's FASTA,
    /// gives every feature the genome's bit in the report and every genome
    /// the final model's fitted call, for the SCM, CART and lasso.
    #[test]
    fn predict_replays_the_final_models_on_their_own_cohort() {
        let dir = tempfile::tempdir().unwrap();
        let (layout, ids) = replay_fixture(dir.path());
        let data = dir.path().to_str().unwrap();
        let grm = |args: &[&str]| {
            let all = [&["grm", "--data", data][..], args].concat();
            run(Cli::try_parse_from(all).unwrap())
        };
        grm(&["count", "--species", SLUG]).unwrap();
        grm(&["structure", "--species", SLUG]).unwrap();
        grm(&["matrix", "--species", SLUG]).unwrap();
        grm(&["learn", "--species", SLUG]).unwrap();
        let report_path = layout
            .results_dir(SLUG)
            .join(format!("{SLUG}__fixturemycin.json"));
        let report: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&report_path).unwrap()).unwrap();
        let samples: Vec<String> = report["samples"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s["id"].as_str().unwrap().to_string())
            .collect();
        assert_eq!(samples.len(), ids.len());

        let fasta: Vec<String> = samples
            .iter()
            .map(|id| {
                let path = bvbrc::fasta_path(&layout.fasta_dir(), id);
                path.to_str().unwrap().to_string()
            })
            .collect();
        let out = dir.path().join("predictions.jsonl");
        let matrix_dir = layout.matrix_dir(SLUG);
        grm(&[
            "predict",
            "--report",
            report_path.to_str().unwrap(),
            "--fasta",
            &fasta.join(","),
            "--matrix",
            matrix_dir.to_str().unwrap(),
            "--learners",
            "scm,cart,lasso",
            "--out",
            out.to_str().unwrap(),
        ])
        .unwrap();
        let rows: Vec<serde_json::Value> = std::fs::read_to_string(&out)
            .unwrap()
            .lines()
            .map(|l| serde_json::from_str(l).unwrap())
            .collect();
        assert_eq!(rows.len(), 3 * samples.len());

        // The final models as `learn` fitted them, on the dataset matrix
        // they were fitted on.
        let matrix = FeatureMatrix::load(&matrix_dir).unwrap();
        let columns: Vec<usize> = samples
            .iter()
            .map(|id| matrix.summary.genomes.iter().position(|g| g == id).unwrap())
            .collect();
        let resistant: Vec<bool> = report["samples"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s["resistant"].as_bool().unwrap())
            .collect();
        let min_support = report["counts"]["min_support"].as_u64().unwrap() as u32;
        let ds = grm::matrix::Dataset::project(&matrix, columns, &resistant, min_support);
        let models = &report["models"];
        let scm: grm::learn::scm::ScmModel = serde_json::from_value(models["scm"].clone()).unwrap();
        let cart: grm::learn::tree::Tree = serde_json::from_value(models["cart"].clone()).unwrap();
        let depth = models["cart_depth"].as_u64().unwrap() as u32;
        let lasso: grm::learn::lasso::LassoModel =
            serde_json::from_value(models["lasso"].clone()).unwrap();
        let fitted_scm = scm.predict(&ds.x);
        let fitted_cart = cart.predict(&ds.x, depth);
        let fitted_lasso = lasso.predict_proba(&ds.x);
        // The fixture's models use presence and absence terms, and the
        // lasso is listed in full.
        let presence: BTreeSet<bool> = scm.rules.iter().map(|r| r.presence).collect();
        assert_eq!(presence, BTreeSet::from([false, true]), "{:?}", scm.rules);
        assert!(cart.nodes.iter().any(|n| n.split.is_some()));
        assert!(!lasso.coef.is_empty());

        let call = |r: bool| if r { "resistant" } else { "susceptible" };
        let mut found: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
        for (i, id) in samples.iter().enumerate() {
            let of = |learner: &str| {
                rows.iter()
                    .find(|r| r["genome"] == id.as_str() && r["learner"] == learner)
                    .unwrap()
            };
            let (s, c, l) = (of("scm"), of("cart"), of("lasso"));
            assert_eq!(s["call"], call(fitted_scm.contains(i)), "{id}");
            assert_eq!(s["call"], call(resistant[i]), "{id}");
            assert_eq!(c["score"].as_f64().unwrap(), fitted_cart[i], "{id}");
            assert_eq!(c["call"], call(fitted_cart[i] >= 0.5), "{id}");
            assert_eq!(l["score"].as_f64().unwrap(), fitted_lasso[i], "{id}");
            assert_eq!(l["call"], call(fitted_lasso[i] >= 0.5), "{id}");
            for row in [s, c, l] {
                assert_eq!(row["partial_features"], 0, "{id}");
                for f in row["features"].as_array().unwrap() {
                    let carried = f["kmers_found"].as_u64().unwrap();
                    let all = f["n_kmers"].as_u64().unwrap();
                    let present = f["found"].as_bool().unwrap();
                    assert_eq!(carried, if present { all } else { 0 }, "{id}: {f}");
                    let list = found.entry(f["feature"].as_u64().unwrap()).or_default();
                    if present && list.last() != Some(&i) {
                        list.push(i);
                    }
                }
            }
        }
        // Each model feature's values over the cohort are its published bits.
        let published: BTreeMap<u64, &str> = report["features"]
            .as_array()
            .unwrap()
            .iter()
            .map(|f| (f["feature"].as_u64().unwrap(), f["bits"].as_str().unwrap()))
            .collect();
        for (feature, carriers) in &found {
            let bits = grm::bits::Bits::from_indices(samples.len(), carriers.iter().copied());
            assert_eq!(
                grm::analysis::base64(bits.words()),
                published[feature],
                "feature {feature}"
            );
        }
        // A published k-mer reads on one strand of its reference contig, and
        // some carriers hold only the other: the presence feature on A.
        let a = scm.rules.iter().find(|r| r.presence).unwrap().feature as u64;
        let kmer = rows[0]["features"]
            .as_array()
            .unwrap()
            .iter()
            .find(|f| f["feature"] == a)
            .unwrap()["kmer"]
            .as_str()
            .unwrap()
            .to_string();
        let mut strands = BTreeSet::new();
        for &i in &found[&a] {
            let path = bvbrc::fasta_path(&layout.fasta_dir(), &samples[i]);
            let text = String::from_utf8(std::fs::read(path).unwrap()).unwrap();
            strands.insert(text.contains(&kmer));
        }
        assert_eq!(strands, BTreeSet::from([false, true]));
    }

    const EXPORT_SLUG: &str = "escherichia_coli";
    const MEROPENEM: &str = "escherichia_coli__meropenem";
    const CIPROFLOXACIN: &str = "escherichia_coli__ciprofloxacin";

    /// A synthetic species for `grm export`, curated and built as
    /// `grm curate` and the pipeline's count, structure and matrix leave it.
    /// Its 60 genomes ("562.400" to "562.459") sit in 20 lineages of three,
    /// each lineage with its own copy of one random core, 2% of its bases
    /// drawn again. Each genome also carries short random blocks, each on its
    /// own contig, so each block is one matrix row for the genomes that hold it:
    /// - 24 blocks in three genomes each, the datasets' support floor, so an
    ///   outer fold that holds out one of the three leaves the block in two
    ///   training genomes;
    /// - 12 blocks in one set of eight genomes and one other genome each,
    ///   which differ from each other only on that genome;
    /// - 40 blocks in random sets.
    ///
    /// Genome 2 is a copy of genome 1, and genome 5 of genome 4: two
    /// k-mer-identical duplicate groups. Meropenem calls one genome of each
    /// lineage resistant (20 of 60), ciprofloxacin every fourth (15 of 60).
    fn export_fixture(root: &Path) -> Layout {
        let layout = Layout::new(root);
        let mut rng = ChaCha8Rng::seed_from_u64(21);
        let random = |rng: &mut ChaCha8Rng, n: usize| -> Vec<u8> {
            (0..n).map(|_| b"ACGT"[rng.random_range(0..4)]).collect()
        };
        let n = 60;
        let core = random(&mut rng, 6000);
        let lineage_cores: Vec<Vec<u8>> = (0..n / 3)
            .map(|_| {
                core.iter()
                    .map(|&b| {
                        if rng.random_bool(0.02) {
                            b"ACGT"[rng.random_range(0..4)]
                        } else {
                            b
                        }
                    })
                    .collect()
            })
            .collect();
        let copy_of = |g: usize| match g {
            2 => 1,
            5 => 4,
            g => g,
        };
        let units: Vec<usize> = (0..n).filter(|&g| copy_of(g) == g).collect();
        // Genomes with no copy, so a block they alone carry keeps its count.
        let single: Vec<usize> = units
            .iter()
            .copied()
            .filter(|&g| g != 1 && g != 4)
            .collect();
        let mut holders: Vec<Vec<usize>> = Vec::new();
        for _ in 0..24 {
            let mut set = BTreeSet::new();
            while set.len() < 3 {
                set.insert(single[rng.random_range(0..single.len())]);
            }
            holders.push(set.into_iter().collect());
        }
        let shared: Vec<usize> = (1..=8).map(|i| 6 * i).collect();
        for &extra in single.iter().filter(|g| !shared.contains(g)).take(12) {
            holders.push(shared.iter().copied().chain([extra]).collect());
        }
        for _ in 0..40 {
            holders.push(
                units
                    .iter()
                    .copied()
                    .filter(|_| rng.random_bool(0.4))
                    .collect(),
            );
        }
        let blocks: Vec<Vec<u8>> = holders.iter().map(|_| random(&mut rng, 80)).collect();
        std::fs::create_dir_all(layout.fasta_dir()).unwrap();
        std::fs::create_dir_all(root.join("bvbrc")).unwrap();
        let (mut amr, mut meta) = (String::new(), String::new());
        for g in 0..n {
            let id = format!("562.{}", 400 + g);
            let unit = copy_of(g);
            let mut fasta = b">core\n".to_vec();
            fasta.extend(&lineage_cores[unit / 3]);
            for (j, (set, block)) in holders.iter().zip(&blocks).enumerate() {
                if set.contains(&unit) {
                    fasta.extend(format!("\n>b{j}\n").into_bytes());
                    fasta.extend(block);
                }
            }
            fasta.push(b'\n');
            std::fs::write(bvbrc::fasta_path(&layout.fasta_dir(), &id), fasta).unwrap();
            for (antibiotic, resistant) in
                [("meropenem", g % 3 == 0), ("ciprofloxacin", g % 4 == 0)]
            {
                let call = if resistant {
                    "Resistant"
                } else {
                    "Susceptible"
                };
                let record = serde_json::json!({
                    "genome_id": id, "antibiotic": antibiotic, "resistant_phenotype": call,
                    "laboratory_typing_method": "Broth dilution",
                });
                amr += &format!("{record}\n");
            }
            let genome = serde_json::json!({
                "genome_id": id, "species": "Escherichia coli", "genome_length": "6000",
            });
            meta += &format!("{genome}\n");
        }
        std::fs::write(layout.amr(), amr).unwrap();
        std::fs::write(layout.genomes(), meta).unwrap();
        run_curate(&layout, 100, 5, None).unwrap();
        count_to_matrix(&layout, EXPORT_SLUG, PIPELINE_K).unwrap();
        layout
    }

    /// `grm --data <root> export --species escherichia_coli --out <out>`
    /// followed by `args`, with an error as `{:#}` prints it.
    fn run_export(root: &Path, out: &Path, args: &[&str]) -> Result<(), String> {
        let (data, out) = (root.to_str().unwrap(), out.to_str().unwrap());
        let base = [
            "grm",
            "--data",
            data,
            "export",
            "--species",
            EXPORT_SLUG,
            "--out",
            out,
        ];
        let cli = Cli::try_parse_from(base.iter().chain(args)).unwrap();
        run(cli).map_err(|e| format!("{e:#}"))
    }

    /// The bytes of a dataset directory's two files.
    fn export_files(dir: &Path) -> [Vec<u8>; 2] {
        ["matrix.bin", "dataset.json"].map(|f| std::fs::read(dir.join(f)).unwrap())
    }

    /// Candidate rows of a `matrix.bin` on the training samples, counted
    /// apart from `learn::training_candidates`: the distinct training
    /// patterns whose training support lies in [floor, |train| - floor].
    fn candidates_in(bytes: &[u8], words: usize, train: &[bool], floor: usize) -> usize {
        let n_train = train.iter().filter(|&&t| t).count();
        let rows = bytes.len() / (words * 8);
        let mut seen = std::collections::HashSet::new();
        for r in 0..rows {
            let pattern: Vec<bool> = (0..train.len())
                .filter(|&c| train[c])
                .map(|c| {
                    let w = (r * words + c / 64) * 8;
                    let word = u64::from_le_bytes(bytes[w..w + 8].try_into().unwrap());
                    (word >> (c % 64)) & 1 == 1
                })
                .collect();
            let support = pattern.iter().filter(|&&b| b).count();
            if support >= floor && support + floor <= n_train {
                seen.insert(pattern);
            }
        }
        seen.len()
    }

    /// `grm export` on a synthetic species, read back from the files it
    /// wrote, against what `grm learn` reads and evaluates for the same
    /// dataset at a non-default `--screened-features`: the evaluation
    /// revision, the input fingerprint, the cohort, the matrix, the outer
    /// folds `evaluate` scores with the duplicate groups, and per outer fold
    /// the training samples and the candidate rows at the dataset's own
    /// support floor.
    #[test]
    fn export_files_hold_what_learn_evaluates_on_a_synthetic_species() {
        use grm::learn::eval::EvalConfig;
        let dir = tempfile::tempdir().unwrap();
        let layout = export_fixture(dir.path());
        let out = dir.path().join("export");
        let cfg = EvalConfig {
            screened_features: 7,
            ..Default::default()
        };
        let only = ["meropenem".to_string()];
        let done =
            grm::export::export_species(&layout, EXPORT_SLUG, &only, &out, &cfg, false).unwrap();
        assert_eq!(done.len(), 1);
        assert_eq!((done[0].id.as_str(), done[0].replaced), (MEROPENEM, false));
        let [bytes, json] = export_files(&out.join(MEROPENEM));
        let d: serde_json::Value = serde_json::from_slice(&json).unwrap();

        let (input, catalog) = grm::analysis::load_inputs(&layout, EXPORT_SLUG).unwrap();
        let spec = input.datasets.iter().find(|s| s.id == MEROPENEM).unwrap();
        let e = grm::analysis::arm_inputs(&input, spec).unwrap().unwrap();
        assert_eq!(d["format"], grm::export::EXPORT_FORMAT);
        assert_eq!(d["evaluation_revision"], grm::analysis::EVALUATION_REVISION);
        // The fingerprint `grm learn` computes at this budget, which differs
        // from the default budget's.
        let digest = grm::analysis::species_digest(&input);
        let fingerprint = |cfg: &EvalConfig| {
            grm::analysis::dataset_fingerprint(&input, spec, catalog.as_ref(), cfg, digest)
        };
        assert_eq!(d["input_fingerprint"], fingerprint(&cfg));
        assert_ne!(fingerprint(&cfg), fingerprint(&EvalConfig::default()));
        assert_eq!(d["screened_features"], 7);

        let n = e.sample_ids.len();
        assert_eq!(n, 60);
        assert_eq!(d["samples"], serde_json::json!(e.sample_ids));
        let resistant: Vec<bool> = (0..n).map(|i| e.y.contains(i)).collect();
        assert_eq!(d["resistant"], serde_json::json!(resistant));
        assert_eq!(d["lineage"], serde_json::json!(e.lineage));
        let groups: BTreeSet<BTreeSet<&str>> = e
            .duplicates
            .iter()
            .map(|g| g.iter().map(|&i| e.sample_ids[i].as_str()).collect())
            .collect();
        let want = [["562.401", "562.402"], ["562.404", "562.405"]];
        assert_eq!(groups, want.map(BTreeSet::from).into());
        assert_eq!(d["duplicates"], serde_json::json!(e.duplicates));
        assert_eq!(e.min_support, 3);
        assert_eq!(d["min_support"], 3);
        assert_eq!(bytes, grm::export::matrix_bytes(&e.x));
        assert_eq!(
            (d["rows"].as_u64(), d["cols"].as_u64(), d["words"].as_u64()),
            (
                Some(e.x.rows as u64),
                Some(n as u64),
                Some(e.x.words as u64)
            )
        );

        let scored = grm::learn::eval::evaluate(
            MEROPENEM,
            &e.x,
            &e.y,
            &e.lineage,
            &e.duplicates,
            e.min_support,
            None,
            &cfg,
        )
        .unwrap();
        let splits = d["splits"].as_array().unwrap();
        assert_eq!(splits.len(), scored.splits.len());
        let (mut below_floor, mut on_test_samples, mut folds) = (false, false, 0);
        for (split, s) in splits.iter().zip(&scored.splits) {
            let plan = &split["plan"];
            assert_eq!(plan["scheme"], s.scheme.as_str());
            assert!(plan["unavailable_reason"].is_null(), "{}", s.scheme);
            assert!(s.unavailable_reason.is_none(), "{}", s.scheme);
            let outer: Vec<usize> = serde_json::from_value(plan["outer"].clone()).unwrap();
            assert_eq!(outer, s.test_fold, "{} outer folds", s.scheme);
            for g in &e.duplicates {
                assert!(
                    g.iter().all(|&i| outer[i] == outer[g[0]]),
                    "{} {g:?}",
                    s.scheme
                );
            }
            let screens = split["screens"].as_array().unwrap();
            assert_eq!(screens.len(), cfg.folds);
            for (f, screen) in screens.iter().enumerate() {
                assert_eq!(screen["fold"], f);
                let train: Vec<bool> = outer.iter().map(|&o| o != f).collect();
                let mask = grm::bits::Bits::from_indices(n, (0..n).filter(|&i| train[i]));
                assert_eq!(screen["n_train"], train.iter().filter(|&&t| t).count());
                let (cx, rows) = grm::learn::training_candidates(&e.x, &mask, e.min_support);
                let candidates = candidates_in(&bytes, e.x.words, &train, 3);
                assert_eq!(screen["n_candidates"], cx.rows, "{} fold {f}", s.scheme);
                assert_eq!(screen["n_candidates"], candidates, "{} fold {f}", s.scheme);
                let want: Vec<u32> = grm::learn::screen(&cx, &e.y, &mask, 7)
                    .into_iter()
                    .map(|r| rows[r as usize])
                    .collect();
                assert_eq!(screen["screened_rows"], serde_json::json!(want));
                // The fixture is one where a floor of 2, or candidates taken
                // on every sample, give other counts.
                below_floor |= candidates_in(&bytes, e.x.words, &train, 2) != candidates;
                on_test_samples |= candidates_in(&bytes, e.x.words, &[true; 60], 3) != candidates;
                folds += 1;
            }
        }
        assert_eq!(folds, 2 * cfg.folds);
        assert!(below_floor && on_test_samples);
    }

    /// An export stops at a dataset directory that exists, names it and
    /// writes nothing for any dataset; `--force` replaces the two files.
    #[test]
    fn export_stops_at_a_dataset_directory_that_exists_unless_forced() {
        let dir = tempfile::tempdir().unwrap();
        let layout = export_fixture(dir.path());
        let out = dir.path().join("export");
        let (mero, cipro) = (out.join(MEROPENEM), out.join(CIPROFLOXACIN));
        run_export(dir.path(), &out, &["--only", "meropenem"]).unwrap();
        let first = export_files(&mero);
        std::fs::write(mero.join("dataset.json"), "earlier").unwrap();
        std::fs::write(mero.join("notes.txt"), "kept").unwrap();

        let err = run_export(dir.path(), &out, &[]).unwrap_err();
        assert!(err.contains(&mero.display().to_string()), "{err}");
        assert!(err.contains("--force"), "{err}");
        assert!(!cipro.exists(), "a refused export wrote another dataset");
        assert_eq!(
            std::fs::read(mero.join("dataset.json")).unwrap(),
            b"earlier"
        );

        run_export(dir.path(), &out, &["--force"]).unwrap();
        assert_eq!(export_files(&mero), first);
        assert_eq!(std::fs::read(mero.join("notes.txt")).unwrap(), b"kept");
        assert!(cipro.join("dataset.json").exists());
        let cfg = grm::learn::eval::EvalConfig::default();
        let done =
            grm::export::export_species(&layout, EXPORT_SLUG, &[], &out, &cfg, true).unwrap();
        assert_eq!(done.len(), 2);
        assert!(done.iter().all(|d| d.replaced), "{done:?}");
        assert_eq!(export_files(&mero), first);

        // Anything at a dataset's path stops the export before it writes
        // any dataset, a file or a broken link too, whichever dataset's
        // path it is.
        for (n, (id, rest, link)) in [
            (CIPROFLOXACIN, MEROPENEM, false),
            (MEROPENEM, CIPROFLOXACIN, false),
            (CIPROFLOXACIN, MEROPENEM, true),
            (MEROPENEM, CIPROFLOXACIN, true),
        ]
        .into_iter()
        .enumerate()
        {
            let other = dir.path().join(format!("other-{n}"));
            std::fs::create_dir_all(&other).unwrap();
            if link {
                std::os::unix::fs::symlink(dir.path().join("nowhere"), other.join(id)).unwrap();
            } else {
                std::fs::write(other.join(id), "a file").unwrap();
            }
            let err = run_export(dir.path(), &other, &[]).unwrap_err();
            let want = format!(
                "1 dataset directory already exists: {}.",
                other.join(id).display()
            );
            assert!(err.contains(&want), "{id} link {link}: {err}");
            assert!(
                !other.join(rest).exists(),
                "{id} link {link}: {rest} written"
            );
        }
    }

    /// Each file is written under a temporary name and renamed, after the
    /// earlier `dataset.json` is removed and before the new one is renamed
    /// into place: a write that fails names its file and leaves a directory
    /// without `dataset.json`, which the next export without `--force`
    /// refuses.
    #[test]
    fn an_export_that_fails_to_write_names_the_file_and_leaves_no_dataset_json() {
        let dir = tempfile::tempdir().unwrap();
        export_fixture(dir.path());
        let out = dir.path().join("export");
        let mero = out.join(MEROPENEM);
        let args = ["--only", "meropenem", "--force"];
        run_export(dir.path(), &out, &args).unwrap();
        let first = export_files(&mero);
        for blocked in ["matrix.bin.tmp", "dataset.json.tmp"] {
            // A directory at the temporary name makes its write fail.
            std::fs::create_dir(mero.join(blocked)).unwrap();
            let err = run_export(dir.path(), &out, &args).unwrap_err();
            assert!(
                err.contains(&mero.join(blocked).display().to_string()),
                "{err}"
            );
            assert!(
                !mero.join("dataset.json").exists(),
                "{blocked}: dataset.json left"
            );
            assert!(mero.join("matrix.bin").exists());
            std::fs::remove_dir(mero.join(blocked)).unwrap();
        }
        let err = run_export(dir.path(), &out, &["--only", "meropenem"]).unwrap_err();
        assert!(err.contains(&mero.display().to_string()), "{err}");
        run_export(dir.path(), &out, &args).unwrap();
        assert_eq!(export_files(&mero), first);
        let names: BTreeSet<String> = std::fs::read_dir(&mero)
            .unwrap()
            .map(|e| e.unwrap().file_name().into_string().unwrap())
            .collect();
        assert_eq!(
            names,
            BTreeSet::from(["dataset.json".into(), "matrix.bin".into()])
        );
    }

    /// `--only` names that match no dataset of the species stop the export
    /// before it writes anything, and the message lists them.
    #[test]
    fn export_stops_on_an_only_name_that_matches_no_dataset() {
        let dir = tempfile::tempdir().unwrap();
        export_fixture(dir.path());
        let out = dir.path().join("export");
        for (only, unmatched) in [
            ("meropenem,ciprofloxacn", "ciprofloxacn"),
            ("nitrofurantoin", "nitrofurantoin"),
            (
                "klebsiella_oxytoca__meropenem",
                "klebsiella_oxytoca__meropenem",
            ),
            (
                "tetracycline,meropenem,ampicilin",
                "tetracycline, ampicilin",
            ),
        ] {
            let err = run_export(dir.path(), &out, &["--only", only]).unwrap_err();
            let want = format!("--only names no dataset of {EXPORT_SLUG}: {unmatched}.");
            assert!(err.contains(&want), "{only}: {err}");
            assert!(
                err.contains(MEROPENEM) && err.contains(CIPROFLOXACIN),
                "{err}"
            );
            assert!(!out.exists(), "{only}: something was written");
        }
        // An antibiotic and a dataset id both match.
        run_export(
            dir.path(),
            &out,
            &["--only", &format!("meropenem,{CIPROFLOXACIN}")],
        )
        .unwrap();
        assert!(out.join(MEROPENEM).join("dataset.json").exists());
        assert!(out.join(CIPROFLOXACIN).join("dataset.json").exists());
    }

    /// A genomes.txt changed after counting stops the export with the
    /// message `grm learn` gives, before anything is written.
    #[test]
    fn export_stops_as_learn_does_when_genomes_txt_changed_after_counting() {
        let dir = tempfile::tempdir().unwrap();
        let layout = export_fixture(dir.path());
        let path = layout.species_dir(EXPORT_SLUG).join("genomes.txt");
        let mut ids = read_lines(&path).unwrap();
        ids.swap(0, 1);
        std::fs::write(&path, ids.join("\n") + "\n").unwrap();
        let data = dir.path().to_str().unwrap();
        let learn = ["grm", "--data", data, "learn", "--species", EXPORT_SLUG];
        let learn_err = format!(
            "{:#}",
            run(Cli::try_parse_from(learn).unwrap()).unwrap_err()
        );
        let out = dir.path().join("export");
        let err = run_export(dir.path(), &out, &[]).unwrap_err();
        assert_eq!(err, learn_err);
        assert!(
            err.contains("rerun `grm count --species escherichia_coli`"),
            "{err}"
        );
        assert!(!out.exists());
    }

    /// A species without a matrix: the message names the three commands
    /// that build it, and a slug with no species folder is called out too.
    #[test]
    fn export_of_a_species_without_a_matrix_names_the_commands_to_run() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        let data = dir.path().to_str().unwrap();
        let out = dir.path().join("export");
        let export = |slug: &str| {
            let args = ["grm", "--data", data, "export", "--species", slug, "--out"];
            let cli = Cli::try_parse_from(args.iter().copied().chain([out.to_str().unwrap()]));
            format!("{:#}", run(cli.unwrap()).unwrap_err())
        };
        let err = export(SLUG);
        for step in ["count", "structure", "matrix"] {
            assert!(
                err.contains(&format!("`grm {step} --species {SLUG}`")),
                "{err}"
            );
        }
        assert!(
            err.contains(&layout.matrix_dir(SLUG).display().to_string()),
            "{err}"
        );
        assert!(!err.contains("check the species slug"), "{err}");
        let err = export("fixtur");
        assert!(err.contains("`grm count --species fixtur`"), "{err}");
        assert!(err.contains("check the species slug"), "{err}");
        assert!(!out.exists());
    }

    /// Every file under `dir`, by its path below `dir`.
    fn files_below(dir: &Path) -> BTreeMap<PathBuf, Vec<u8>> {
        let mut out = BTreeMap::new();
        let mut todo = vec![dir.to_path_buf()];
        while let Some(d) = todo.pop() {
            for entry in std::fs::read_dir(&d).unwrap() {
                let path = entry.unwrap().path();
                if path.is_dir() {
                    todo.push(path);
                } else {
                    let bytes = std::fs::read(&path).unwrap();
                    out.insert(path.strip_prefix(dir).unwrap().to_path_buf(), bytes);
                }
            }
        }
        out
    }

    /// A report file with its run clocks set to zero: `generated_at` in
    /// `index.json` and each learner's `seconds` in a dataset report. Both
    /// differ between two runs of one build on one root.
    fn without_clocks(bytes: &[u8]) -> serde_json::Value {
        fn zero(v: &mut serde_json::Value) {
            match v {
                serde_json::Value::Object(map) => {
                    for (key, value) in map.iter_mut() {
                        if key == "generated_at" || key == "seconds" {
                            *value = serde_json::Value::Null;
                        } else {
                            zero(value);
                        }
                    }
                }
                serde_json::Value::Array(items) => items.iter_mut().for_each(zero),
                _ => {}
            }
        }
        let mut v: serde_json::Value = serde_json::from_slice(bytes).unwrap();
        zero(&mut v);
        v
    }

    /// Ledger action 11: a data root fetched before `GENOME_FIELDS` asked
    /// for `collection_date` reads as before, and adding the field to every
    /// genome's metadata changes no report byte. `grm learn` keeps the
    /// reports it fitted as current; a forced refit and `grm report` give
    /// the same files, apart from the two run clocks.
    #[test]
    fn collection_dates_leave_every_report_file_as_it_was() {
        let dir = tempfile::tempdir().unwrap();
        let (layout, ids) = replay_fixture(dir.path());
        let data = dir.path().to_str().unwrap();
        let grm = |args: &[&str]| {
            let all = [&["grm", "--data", data][..], args].concat();
            run(Cli::try_parse_from(all).unwrap())
        };
        // `grm report` finds a species' files under the slug of its name, so
        // the fixture's species is named for its folder.
        for file in ["datasets.json", "phenotypes.jsonl"] {
            let path = layout.curated().join(file);
            let text = std::fs::read_to_string(&path).unwrap();
            std::fs::write(&path, text.replace(r#""Fixture species""#, r#""Fixture""#)).unwrap();
        }
        std::fs::write(
            layout.curated().join("funnel.json"),
            serde_json::to_vec(&curate::Funnel::default()).unwrap(),
        )
        .unwrap();
        std::fs::write(layout.amr(), "").unwrap();
        let kover = dir.path().join("kover.json");
        std::fs::write(&kover, "[]").unwrap();
        let kover = kover.to_str().unwrap();
        // Synthetic metadata, in the form `grm fetch genomes` wrote before
        // the field (`bvbrc::tests::genome_lines_fetched_before_collection_date_read_and_write_as_before`).
        let genomes: Vec<GenomeRecord> = ids
            .iter()
            .enumerate()
            .map(|(i, id)| GenomeRecord {
                genome_id: id.clone(),
                collection_year: (i % 4 != 3).then(|| format!("{}", 2010 + i % 7)),
                isolation_country: Some(["Kenya", "Peru", "Viet Nam"][i % 3].into()),
                host_name: Some("Human, Homo sapiens".into()),
                mlst: Some(format!("MLST.Fixture.{}", i % 5)),
                ..Default::default()
            })
            .collect();
        let write_genomes = |records: &[GenomeRecord]| {
            let lines: Vec<String> = records
                .iter()
                .map(|g| serde_json::to_string(g).unwrap())
                .collect();
            std::fs::write(layout.genomes(), lines.join("\n") + "\n").unwrap();
        };
        write_genomes(&genomes);
        assert!(
            !std::fs::read_to_string(layout.genomes())
                .unwrap()
                .contains("collection_date")
        );
        grm(&["count", "--species", SLUG]).unwrap();
        grm(&["structure", "--species", SLUG]).unwrap();
        grm(&["matrix", "--species", SLUG]).unwrap();
        grm(&["learn", "--species", SLUG]).unwrap();
        grm(&["report", "--kover", kover]).unwrap();
        let report_dir = layout.root.join("report");
        let before = files_below(&report_dir);
        let results_before = files_below(&layout.results_dir(SLUG));
        assert!(before.contains_key(Path::new("index.json")));
        assert!(before.contains_key(Path::new("corpus.json")));
        assert!(before.contains_key(Path::new("species/fixture.json")));
        assert!(before.contains_key(Path::new("datasets/fixture__fixturemycin.json")));
        let species: serde_json::Value =
            serde_json::from_slice(&before[Path::new("species/fixture.json")]).unwrap();
        assert_eq!(species["genomes"][0]["year"], "2010");

        // The same genomes with a collection date of every kind.
        let texts = ["2015-03-12", "2014-05", "2013", "2015-2016", "2017-02-30"];
        let dated: Vec<GenomeRecord> = genomes
            .iter()
            .enumerate()
            .map(|(i, g)| GenomeRecord {
                collection_date: match texts.get(i % 6) {
                    Some(text) => bvbrc::CollectionDate::parse(text),
                    None => bvbrc::CollectionDate::Missing,
                },
                ..g.clone()
            })
            .collect();
        write_genomes(&dated);
        let read: Vec<GenomeRecord> = bvbrc::read_jsonl(&layout.genomes()).unwrap();
        let unparsed = read
            .iter()
            .filter(|g| g.collection_date.is_unparsed())
            .count();
        assert_eq!(unparsed, 20);

        // Learn keeps its reports as current, and the report is rebuilt
        // byte for byte apart from the index's build time.
        grm(&["learn", "--species", SLUG]).unwrap();
        assert_eq!(files_below(&layout.results_dir(SLUG)), results_before);
        grm(&["report", "--kover", kover]).unwrap();
        let after = files_below(&report_dir);
        assert_eq!(
            after.keys().collect::<Vec<_>>(),
            before.keys().collect::<Vec<_>>()
        );
        for (path, bytes) in &after {
            if path == Path::new("index.json") {
                assert_eq!(without_clocks(bytes), without_clocks(&before[path]));
            } else {
                assert!(bytes == &before[path], "{} changed", path.display());
            }
        }

        // A forced refit reads the dated metadata and fits the same models.
        grm(&["learn", "--species", SLUG, "--force"]).unwrap();
        grm(&["report", "--kover", kover]).unwrap();
        let refit = files_below(&report_dir);
        assert_eq!(
            refit.keys().collect::<Vec<_>>(),
            before.keys().collect::<Vec<_>>()
        );
        for (path, bytes) in &refit {
            assert_eq!(
                without_clocks(bytes),
                without_clocks(&before[path]),
                "{}",
                path.display()
            );
        }
    }

    /// Ledger action 11: `grm curate` reads `genomes.jsonl` for species, QC
    /// and MLST groups. Run on a synthetic raw root before and after every
    /// genome gets a collection date, with each selection rule, it writes the
    /// same curated files and species files. The one allowance is the key
    /// order of `qc.json`'s failure map, a `HashMap` written as it iterates,
    /// which differs between two runs of `main` too; it is compared as JSON.
    /// No curated file holds the word `collection` or a date text of the
    /// fixture.
    #[test]
    fn collection_dates_leave_every_curated_file_as_it_was() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let data = dir.path().to_str().unwrap();
        std::fs::create_dir_all(dir.path().join("bvbrc")).unwrap();
        // Synthetic: 44 genomes of one species with a meropenem call each,
        // 22 resistant, three that fail QC, and one AST record whose genome
        // has no metadata line. The label dates are 2019 to 2021, apart
        // from every collection date below.
        let ids: Vec<String> = (0..44).map(|i| format!("9990006.{}", 100 + i)).collect();
        let mut amr = String::new();
        for (i, id) in ids
            .iter()
            .chain(["9990006.999".to_string()].iter())
            .enumerate()
        {
            let record = serde_json::json!({
                "genome_id": id, "antibiotic": "meropenem",
                "resistant_phenotype": if i % 2 == 0 { "Resistant" } else { "Susceptible" },
                "laboratory_typing_method": "Broth dilution", "testing_standard": "CLSI",
                "testing_standard_year": 2019 + i % 3, "pmid": 30000000 + i % 4,
                "date_inserted": format!("{}-0{}-1{}T08:00:00Z", 2019 + i % 3, 1 + i % 9, i % 10),
            });
            amr += &format!("{record}\n");
        }
        std::fs::write(layout.amr(), amr).unwrap();
        let genomes: Vec<GenomeRecord> = ids
            .iter()
            .enumerate()
            .map(|(i, id)| GenomeRecord {
                genome_id: id.clone(),
                species: Some("Escherichia coli".into()),
                genome_quality: Some(if i == 5 { "Poor" } else { "Good" }.into()),
                genome_length: (i != 6).then(|| format!("{}", 100_000 + 37 * i)),
                contigs: Some(if i == 7 { "5000" } else { "80" }.into()),
                collection_year: (i % 4 != 3).then(|| format!("{}", 2010 + i % 7)),
                mlst: (i % 5 != 4).then(|| format!("MLST.Escherichia_coli_1.{}", i % 6)),
                ..Default::default()
            })
            .collect();
        let write_genomes = |records: &[GenomeRecord]| {
            let lines: Vec<String> = records
                .iter()
                .map(|g| serde_json::to_string(g).unwrap())
                .collect();
            std::fs::write(layout.genomes(), lines.join("\n") + "\n").unwrap();
        };
        // Every kind of collection date: day, month, year, a range, a day
        // that does not exist, a cut range, and no value.
        let texts = [
            "2015-03-12",
            "2014-05",
            "1987",
            "2015-2016",
            "2017-02-30",
            "2004-20",
        ];
        let dated: Vec<GenomeRecord> = genomes
            .iter()
            .enumerate()
            .map(|(i, g)| GenomeRecord {
                collection_date: match texts.get(i % 7) {
                    Some(text) => bvbrc::CollectionDate::parse(text),
                    None => bvbrc::CollectionDate::Missing,
                },
                ..g.clone()
            })
            .collect();
        let curated_and_species = || {
            let mut files = files_below(&layout.curated());
            for (path, bytes) in files_below(&dir.path().join("species")) {
                files.insert(Path::new("species").join(path), bytes);
            }
            files
        };
        let rules: [&[&str]; 2] = [
            &[],
            &[
                "--selection",
                "group-round-robin-v1",
                "--selection-seed",
                "7",
            ],
        ];
        for rule in rules {
            let curate = || {
                let args = [
                    &[
                        "grm",
                        "--data",
                        data,
                        "curate",
                        "--cap",
                        "30",
                        "--min-class",
                        "10",
                    ][..],
                    rule,
                ]
                .concat();
                run(Cli::try_parse_from(args).unwrap()).unwrap();
            };
            write_genomes(&genomes);
            assert!(
                !std::fs::read_to_string(layout.genomes())
                    .unwrap()
                    .contains("collection_date")
            );
            curate();
            let before = curated_and_species();
            write_genomes(&dated);
            let lines = std::fs::read_to_string(layout.genomes()).unwrap();
            assert!(lines.lines().all(|l| l.contains(r#""collection_date":"#)));
            curate();
            let after = curated_and_species();

            // The fixture reaches every curated file: one dataset, a
            // selection under the cap, and the three QC failures. The group
            // rule also writes its plan, `selection.json`.
            let plan = usize::from(!rule.is_empty());
            assert_eq!(after.len(), 7 + 1 + plan, "{:?}", after.keys());
            let datasets: serde_json::Value =
                serde_json::from_slice(&after[Path::new("datasets.json")]).unwrap();
            assert_eq!(datasets.as_array().unwrap().len(), 1);
            let selected = &after[Path::new("species/escherichia_coli/genomes.txt")];
            assert_eq!(selected.iter().filter(|&&b| b == b'\n').count(), 30);
            let (_, qc): (serde_json::Value, BTreeMap<String, String>) =
                serde_json::from_slice(&after[Path::new("qc.json")]).unwrap();
            assert_eq!(qc.len(), 3, "{qc:?}");

            assert_eq!(
                after.keys().collect::<Vec<_>>(),
                before.keys().collect::<Vec<_>>()
            );
            for (path, bytes) in &after {
                if path == Path::new("qc.json") {
                    let value = |b: &[u8]| serde_json::from_slice::<serde_json::Value>(b).unwrap();
                    assert_eq!(value(bytes), value(&before[path]), "qc.json");
                } else {
                    assert!(bytes == &before[path], "{} changed", path.display());
                }
                let text = String::from_utf8_lossy(bytes);
                for word in texts.iter().chain(&["collection"]) {
                    assert!(!text.contains(word), "{} holds {word}", path.display());
                }
            }
        }
    }
}
