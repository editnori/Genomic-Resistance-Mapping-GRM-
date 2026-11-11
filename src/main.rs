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
        #[arg(long)]
        drop_ksets: bool,
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
    /// Write offline 95% Clopper-Pearson intervals for each learner row's
    /// VME and ME under <data>/report/intervals/, one file per dataset report
    /// that index.json lists, each pinned to its report's SHA-256 and counts.
    /// No report file changes.
    Intervals,
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
            selection_seed,
            stratum,
        } => {
            let rule = selection_rule(selection, selection_seed, stratum)?;
            run_curate_with(&layout, cap, min_class, sigdrug_vocabulary, &rule)
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
        } => species
            .iter()
            .try_for_each(|s| matrix(&layout, s, min_support)),
        Command::Learn {
            species,
            only,
            force,
            screened_features,
        } => species
            .iter()
            .try_for_each(|s| learn(&layout, s, &only, force, screened_features)),
        Command::Pipeline {
            species,
            workers,
            k,
            screened_features,
            drop_ksets,
            allow_missing_catalog,
        } => {
            require_catalog(&layout, k, allow_missing_catalog)?;
            require_kover(std::path::Path::new(KOVER_RESULTS))?;
            let slugs = if species.is_empty() {
                curated_species(&layout)?
            } else {
                species
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
                count_to_matrix(&layout, slug, k)?;
                if let Err(e) = learn(&layout, slug, &[], false, screened_features) {
                    if !e.is::<grm::probes::KsetsUnreadable>() {
                        return Err(e);
                    }
                    eprintln!("{slug}: recounting ksets: {e:#}");
                    count(&layout, slug, k)?;
                    count_to_matrix(&layout, slug, k)?;
                    learn(&layout, slug, &[], false, screened_features)?;
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
        Command::Intervals => {
            let n = grm::intervals::write(&layout.root.join("report"))?;
            eprintln!("intervals written for {n} dataset reports");
            Ok(())
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

/// Count, structure, and matrix at `k`, each rerun only when `stale_stages`
/// finds its output out of date.
fn count_to_matrix(layout: &Layout, slug: &str, k: usize) -> Result<()> {
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
        Catalog::load(&layout.catalog_dir(), &grm::dna::KmerSpace::new(k)).context(
            "the AMRFinderPlus catalog does not load, and counting without it records no gene calls, which leaves reports without the catalog baseline and gene panel detections; pass --allow-missing-catalog to count anyway",
        )?;
    }
    Ok(())
}

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
    eprintln!(
        "{slug}: {} genomes, {} lineage clusters at Mash distance {} (largest {:.0}%) in {:.0}s",
        s.genomes.len(),
        s.lineages.n_clusters,
        s.lineages.threshold,
        s.lineages.largest_share * 100.0,
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

fn matrix(layout: &Layout, slug: &str, min_support: u32) -> Result<()> {
    let started = std::time::Instant::now();
    let s = grm::species::build_matrix(layout, slug, min_support)?;
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
    Ok(())
}

fn learn(
    layout: &Layout,
    slug: &str,
    only: &[String],
    force: bool,
    screened_features: usize,
) -> Result<()> {
    grm::species::ensure_selection_counted(layout, slug)?;
    let cfg = grm::learn::eval::EvalConfig {
        screened_features,
        ..Default::default()
    };
    let written = grm::analysis::analyze_species(layout, slug, only, force, &cfg)?;
    eprintln!("{slug}: {} dataset reports written", written.len());
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
    run_curate_with(layout, cap, min_class, sigdrug, &Selector::Legacy)
}

fn run_curate_with(
    layout: &Layout,
    cap: usize,
    min_class: usize,
    sigdrug: Option<PathBuf>,
    selector: &Selector,
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
                learn(&layout, SLUG, &[], false, DEFAULT_SCREENED_FEATURES).unwrap_err()
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
        assert!(grm::species::build_matrix(&layout, SLUG, 2).is_err());
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
            screened_features, ..
        } = learn.unwrap().command
        else {
            panic!("parsed another command");
        };
        assert_eq!(screened_features, 1000);
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
        run_curate_with(&layout, 30, folds, None, &rule).unwrap();
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
}
