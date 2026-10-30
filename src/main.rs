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
use std::collections::{BTreeMap, BTreeSet};
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
    },
    /// Compare each curated species' selected genomes with the species'
    /// reference genome (Mash distance and GC) and record the genomes whose
    /// sequence is another species, which `grm curate` then excludes.
    SpeciesCheck {
        /// Only sketch the FASTA files on hand into the sketch cache, for a
        /// check that runs elsewhere.
        #[arg(long)]
        sketch_only: bool,
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
    },
    /// Run fetch, count, structure, matrix, and learn for each species, then
    /// build the report. Stages skip work that is already on disk.
    Pipeline {
        /// Species slugs; every curated species when omitted.
        #[arg(long, value_delimiter = ',')]
        species: Vec<String>,
        #[arg(long, default_value_t = 8)]
        workers: usize,
        /// Delete each species' ksets after its datasets are learned, to cap
        /// disk use at one species' ksets (re-creatable with `grm count`).
        #[arg(long)]
        drop_ksets: bool,
        /// Count without catalog gene calls when the AMRFinderPlus catalog
        /// does not load, instead of stopping (see `grm count`).
        #[arg(long)]
        allow_missing_catalog: bool,
    },
    /// Build the dashboard data under <data>/report.
    Report {
        /// Published Kover results to compare against.
        #[arg(long, default_value = KOVER_RESULTS)]
        kover: PathBuf,
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
        } => run_curate(&layout, cap, min_class, sigdrug_vocabulary),
        Command::SpeciesCheck { sketch_only } => species_check(&layout, sketch_only),
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
        } => species
            .iter()
            .try_for_each(|s| learn(&layout, s, &only, force)),
        Command::Pipeline {
            species,
            workers,
            drop_ksets,
            allow_missing_catalog,
        } => {
            require_catalog(&layout, PIPELINE_K, allow_missing_catalog)?;
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
                count_to_matrix(&layout, slug)?;
                if let Err(e) = learn(&layout, slug, &[], false) {
                    if !e.is::<grm::probes::KsetsUnreadable>() {
                        return Err(e);
                    }
                    eprintln!("{slug}: recounting ksets: {e:#}");
                    count(&layout, slug, PIPELINE_K)?;
                    count_to_matrix(&layout, slug)?;
                    learn(&layout, slug, &[], false)?;
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
        Command::Report { kover } => report(&layout, &kover),
        Command::Research {
            revision5_reports,
            out,
        } => research(&layout, &revision5_reports, &out),
        Command::Serve { app, port } => grm::serve::serve(&app, &layout.root.join("report"), port),
    }
}

const KOVER_RESULTS: &str = "reference/kover-summary.json";

/// k for `grm pipeline`, and the `grm count` default.
const PIPELINE_K: usize = 31;

/// min_support for `grm pipeline`, and the `grm matrix` default.
const PIPELINE_MIN_SUPPORT: u32 = 2;

/// Count, structure, and matrix, each rerun only when `stale_stages` finds
/// its output out of date.
fn count_to_matrix(layout: &Layout, slug: &str) -> Result<()> {
    let catalog = load_catalog(layout, PIPELINE_K);
    let stale_stages = || {
        grm::species::stale_stages(
            layout,
            slug,
            catalog.as_ref(),
            PIPELINE_K,
            PIPELINE_MIN_SUPPORT,
        )
    };
    // Count reruns whenever any stage is stale (see `stale_stages`).
    if !stale_stages().count {
        eprintln!(
            "{slug}: counts match genomes.txt and the catalog, and structure and matrix match the counts; skipping to learn"
        );
        return Ok(());
    }
    count_with(layout, slug, PIPELINE_K, catalog.as_ref())?;
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

fn learn(layout: &Layout, slug: &str, only: &[String], force: bool) -> Result<()> {
    grm::species::ensure_selection_counted(layout, slug)?;
    let cfg = grm::learn::eval::EvalConfig::default();
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

fn run_curate(
    layout: &Layout,
    cap: usize,
    min_class: usize,
    sigdrug: Option<PathBuf>,
) -> Result<()> {
    // `grm learn` gates each dataset at this minimum; below the outer fold
    // count, a random-split fold lacks a class and `grm learn` stops with an
    // error.
    let folds = grm::learn::eval::EvalConfig::default().folds;
    anyhow::ensure!(
        min_class >= folds,
        "--min-class {min_class} is below {folds}: each of the {folds} random-split folds needs a resistant and a susceptible genome"
    );
    let records: Vec<AmrRecord> = bvbrc::read_jsonl(&layout.amr())?;
    let genomes: Vec<GenomeRecord> = bvbrc::read_jsonl(&layout.genomes())?;
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
    let curated = curate::curate(&records, &genomes, &cfg, &wrong_species, &mut resolver)?;

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
    std::fs::write(
        dir.join("datasets.json"),
        serde_json::to_vec_pretty(&curated.datasets)?,
    )?;
    std::fs::write(
        dir.join("funnel.json"),
        serde_json::to_vec_pretty(&curated.funnel)?,
    )?;
    std::fs::write(
        dir.join("qc.json"),
        serde_json::to_vec_pretty(&(&cfg, &curated.qc_failures))?,
    )?;

    let mut by_species: BTreeMap<&str, usize> = BTreeMap::new();
    for d in &curated.datasets {
        *by_species.entry(&d.species).or_default() += 1;
    }
    let mut selections = Vec::new();
    for (species, n_datasets) in by_species {
        let slug = curate::slug(species);
        let labeled = curated
            .phenotypes
            .iter()
            .filter(|p| p.species == species && p.label.is_some())
            .map(|p| &p.genome_id)
            .collect::<BTreeSet<_>>()
            .len();
        let ids = curate::select_genomes(&curated.phenotypes, &curated.datasets, species, cap);
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

fn species_check(layout: &Layout, sketch_only: bool) -> Result<()> {
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
    // Genomes an earlier check judged stay in, and so do its references:
    // curation has since dropped the failures from genomes.txt.
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
    let check = grm::species_check::check(cohorts, &Default::default(), earlier.as_ref())?;
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
        count_to_matrix(&layout, SLUG).unwrap();
        mark_stages(&layout);
        set_catalog_version(&layout, "v2");
        // `grm learn` alone refuses the old calls; the pipeline recounts first.
        let refusal = learn_refusal(&layout).unwrap();
        assert!(
            refusal.contains("has gene calls from catalog v1|group-calls-v2"),
            "{refusal}"
        );
        count_to_matrix(&layout, SLUG).unwrap();
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
            count_to_matrix(&layout, SLUG).unwrap();
            download(&layout, "g4");
            write_genome_list(&layout, selection);
            let err = format!("{:#}", learn(&layout, SLUG, &[], false).unwrap_err());
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
            count_to_matrix(&layout, SLUG).unwrap();
            mark_stages(&layout);
            download(&layout, "g4");
            write_genome_list(&layout, selection);
            count_to_matrix(&layout, SLUG).unwrap();
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
        count_to_matrix(&layout, SLUG).unwrap();
        mark_stages(&layout);
        // `grm count` leaves out a genome without a FASTA file, so selecting
        // g4 before its download changes nothing, the counts included.
        let counted_json = layout.species_dir(SLUG).join("counted.json");
        let counts = std::fs::read(&counted_json).unwrap();
        write_genome_list(&layout, &["g0", "g1", "g2", "g3", "g4"]);
        count_to_matrix(&layout, SLUG).unwrap();
        assert_eq!(std::fs::read(&counted_json).unwrap(), counts);
        assert_eq!(still_marked(&layout), (true, true));
        download(&layout, "g4");
        count_to_matrix(&layout, SLUG).unwrap();
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
        let err = format!("{:#}", count_to_matrix(&layout, SLUG).unwrap_err());
        assert!(
            err.contains("fixture: 1 of 4 selected genomes failed to count, first g1:"),
            "{err}"
        );
        assert!(!layout.species_dir(SLUG).join("structure.json").exists());
        assert!(!layout.matrix_dir(SLUG).join("matrix.json").exists());
        std::fs::remove_dir(&part).unwrap();
        count_to_matrix(&layout, SLUG).unwrap();
        let m = FeatureMatrix::load(&layout.matrix_dir(SLUG)).unwrap();
        assert_eq!(m.summary.genomes, ["g0", "g1", "g2", "g3"]);
        assert_eq!(learn_refusal(&layout), None);
    }

    #[test]
    fn a_count_at_another_k_is_redone_and_rebuilds_only_stages_built_at_that_k() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        count_to_matrix(&layout, SLUG).unwrap();
        mark_stages(&layout);
        // `grm count --k 21` alone: the recount at k = 31 writes again the
        // ksets that structure and matrix were built from.
        count(&layout, SLUG, 21).unwrap();
        count_to_matrix(&layout, SLUG).unwrap();
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
        count_to_matrix(&layout, SLUG).unwrap();
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
        count_to_matrix(&layout, SLUG).unwrap();
        let counted = grm::species::load_counted(&layout, SLUG).unwrap();
        assert_eq!(counted.k, PIPELINE_K);
        assert_eq!(still_marked(&layout), (false, false));
        // `grm matrix --min-support 1` keeps k-mers the pipeline's matrix drops.
        matrix(&layout, SLUG, 1).unwrap();
        mark_stages(&layout);
        count_to_matrix(&layout, SLUG).unwrap();
        assert_eq!(still_marked(&layout), (true, false));
        let m = FeatureMatrix::load(&layout.matrix_dir(SLUG)).unwrap();
        assert_eq!((m.summary.k, m.summary.min_support), (PIPELINE_K, 2));
        assert_eq!(learn_refusal(&layout), None);
    }

    #[test]
    fn a_failed_matrix_resave_is_rebuilt_on_the_next_run() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        count_to_matrix(&layout, SLUG).unwrap();
        let mx = layout.matrix_dir(SLUG);
        let built = FeatureMatrix::load(&mx).unwrap();
        // The re-save rewrites three binaries, then cannot create kmers.bin.
        std::fs::remove_file(mx.join("kmers.bin")).unwrap();
        std::fs::create_dir(mx.join("kmers.bin")).unwrap();
        assert!(grm::species::build_matrix(&layout, SLUG, 2).is_err());
        std::fs::remove_dir(mx.join("kmers.bin")).unwrap();
        count_to_matrix(&layout, SLUG).unwrap();
        let rebuilt = FeatureMatrix::load(&mx).unwrap();
        assert_eq!(rebuilt.kmers, built.kmers);
        assert_eq!(rebuilt.summary.genomes, built.summary.genomes);
    }

    #[test]
    fn counts_from_before_kset_digests_are_redone_and_keep_structure_and_matrix() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        count_to_matrix(&layout, SLUG).unwrap();
        let counted_with_digests = counted_digests(&layout);
        assert!(counted_with_digests.iter().all(Option::is_some));
        mark_stages(&layout);
        drop_ksets_and_digests(&layout);
        let refusal = learn_refusal(&layout).unwrap();
        assert!(refusal.contains("no kset digest for g0"), "{refusal}");
        count_to_matrix(&layout, SLUG).unwrap();
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
            count_to_matrix(&layout, SLUG).unwrap();
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
            count_to_matrix(&layout, SLUG).unwrap();
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
        count_to_matrix(&layout, SLUG).unwrap();
        mark_stages(&layout);
        let structure_json = layout.species_dir(SLUG).join("structure.json");
        remove_fields(&structure_json, &["k", "kset_digests"]);
        assert_eq!(
            learn_refusal(&layout).unwrap(),
            "structure.json records no kset digests (it predates them); rerun `grm structure`"
        );
        count_to_matrix(&layout, SLUG).unwrap();
        assert_eq!(still_marked(&layout), (false, true));
        assert_eq!(learn_refusal(&layout), None);
    }

    #[test]
    fn counting_and_building_the_matrix_by_hand_leave_an_old_structure_refused() {
        let dir = tempfile::tempdir().unwrap();
        let layout = fixture(dir.path());
        count_to_matrix(&layout, SLUG).unwrap();
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
        count_to_matrix(&layout, SLUG).unwrap();
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
        count_to_matrix(&layout, SLUG).unwrap();
        mark_stages(&layout);
        let fasta = bvbrc::fasta_path(&layout.fasta_dir(), "g0");
        let mut sequence = std::fs::read_to_string(&fasta).unwrap();
        sequence.push('A');
        std::fs::write(fasta, sequence).unwrap();
        drop_ksets_and_digests(&layout);
        count_to_matrix(&layout, SLUG).unwrap();
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
        count_to_matrix(&layout, SLUG).unwrap();
        std::fs::remove_dir_all(layout.catalog_dir()).unwrap();
        let refusal = learn_refusal(&layout).unwrap();
        assert!(refusal.contains("which no longer loads"), "{refusal}");
        // A recount, due here for counts from before kset digests, would
        // record no gene calls in place of the ones counted.json holds.
        drop_ksets_and_digests(&layout);
        let err = format!("{:#}", count_to_matrix(&layout, SLUG).unwrap_err());
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
            species_check(&layout, false).unwrap();
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
}
