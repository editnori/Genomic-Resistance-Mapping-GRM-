//! Runs every dataset of one species end to end and writes one report per
//! dataset: samples, cross-validated learner results under both split
//! schemes, the catalog baseline and panel simulations, MIC regression, the
//! final interpretable models, and annotation for every feature they use.

use crate::annotate::{self, FeatureNote};
use crate::bits::{BitMatrix, Bits};
use crate::bvbrc::{self, Client, GenomeRecord};
use crate::catalog::{self, Catalog, FamilyHit};
use crate::curate::{Breakpoint, Call, Censor, DatasetSpec, Mic, Phenotype};
use crate::dna::KmerSpace;
use crate::fnv::Fnv;
use crate::layout::{Layout, remove_if_present};
use crate::learn::eval::{
    self, EvalConfig, FinalModels, LassoFitAudit, LassoOnlyEvent, LassoProgress, LearnerResult,
    MicInput, MicResult, SplitResult,
};
use crate::learn::gbt::GbtParams;
use crate::learn::tree::ForestParams;
use crate::learn::{BinaryMetrics, NoCall, binary_metrics};
use crate::matrix::{Dataset, FeatureMatrix, MatrixSummary};
use crate::panels;
use crate::probes;
use crate::species::{self, Structure};
use anyhow::{Context, Result};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::time::Instant;

/// Schema 2 lets a learner row hold null scores (no call), null fold
/// balanced accuracy, and null metrics when it called no isolate. A schema-1
/// reader stops at its version check instead of reading them as numbers.
pub const REPORT_SCHEMA: u32 = 2;
/// Revision 8 added the lasso solver evidence; revision 9 changes only the
/// lasso's active-set solve (docs/lasso-stopping.md). Revision 10 adds the
/// species identity check (docs/species-check.md) and keeps k-mer-identical
/// duplicates whole in one outer fold (docs/duplicate-isolation.md).
pub const EVALUATION_REVISION: u32 = 10;

#[derive(Serialize)]
pub struct EvaluationAudit {
    pub revision: u32,
    pub valid: bool,
    /// Why the metrics are not valid; empty when they are.
    pub reasons: Vec<String>,
    /// What a valid audit covers.
    pub scope: String,
}

fn evaluation_audit(
    splits: &[SplitResult],
    final_lasso: Option<&LassoFitAudit>,
) -> EvaluationAudit {
    let mut reasons = Vec::new();
    for split in splits.iter().filter(|s| s.unavailable_reason.is_none()) {
        let lasso = split.learners.iter().find(|l| l.learner == "lasso");
        match lasso.and_then(|l| l.lasso_fits.as_deref()) {
            None => reasons.push(format!(
                "{} split has no per-fold lasso solver status",
                split.scheme
            )),
            Some(fits) if fits.len() != split.folds => reasons.push(format!(
                "{} split has {} lasso solver records for {} folds",
                split.scheme,
                fits.len(),
                split.folds
            )),
            Some(fits) => {
                for (fold, fit) in fits.iter().enumerate() {
                    if fit.has_unconverged() {
                        reasons.push(format!(
                            "{} split outer fold {fold} has an unconverged two-class lasso fit",
                            split.scheme
                        ));
                    }
                }
            }
        }
    }
    match final_lasso {
        None => reasons.push("final model has no lasso solver status".into()),
        Some(fit) if fit.has_unconverged() => {
            reasons.push("final model has an unconverged two-class lasso fit".into())
        }
        Some(_) => {}
    }
    EvaluationAudit {
        revision: EVALUATION_REVISION,
        valid: reasons.is_empty(),
        reasons,
        scope: "Internal cross-validation on public BV-BRC genomes scored against reported R/S labels. Feature selection and tuning use only each outer fold's training samples; inner tuning shares that outer-training feature set. MIC thresholds are fitted within training folds. This is not external or clinical validation.".into(),
    }
}

#[derive(Serialize)]
pub struct Counts {
    pub samples: usize,
    pub resistant: usize,
    pub susceptible: usize,
    /// Samples whose label was inferred from MIC: 0, since samples carry
    /// reported calls only.
    pub labels_from_mic: usize,
    /// Matrix genomes whose resistant or susceptible label for this drug was
    /// inferred from MIC, left out of training and testing. Out of
    /// `samples` + this, the matrix genomes with either label.
    pub labels_from_mic_excluded: usize,
    pub features: usize,
    pub patterns: usize,
    pub kmers: u64,
    /// k of the species matrix the features come from.
    pub k: usize,
    pub lineages: usize,
    pub min_support: u32,
    /// Fewest resistant and fewest susceptible samples the dataset needed.
    pub min_class: usize,
}

#[derive(Serialize)]
pub struct Sample {
    pub id: String,
    pub resistant: bool,
    pub from_mic: bool,
    pub mic: Option<Mic>,
    pub lineage: u32,
    pub country: Option<String>,
    pub year: Option<String>,
    pub host: Option<String>,
    pub mlst: Option<String>,
    /// Indices into `catalog_groups`.
    pub groups: Vec<u16>,
    pub catalog_resistant: Option<bool>,
}

#[derive(Serialize)]
pub struct GroupRow {
    pub group: String,
    pub class: String,
    pub subclass: String,
    pub counts_for_rule: bool,
    pub n_resistant: usize,
    pub n_susceptible: usize,
}

#[derive(Serialize)]
pub struct CatalogBaseline {
    pub version: String,
    pub rule: Vec<Vec<String>>,
    pub metrics: BinaryMetrics,
    pub groups: Vec<GroupRow>,
}

#[derive(Serialize)]
pub struct PanelResult {
    pub id: String,
    pub name: String,
    pub mirrors: String,
    pub targets: Vec<String>,
    pub metrics: BinaryMetrics,
    /// Target groups detected, with resistant / susceptible carrier counts.
    pub detected: Vec<(String, usize, usize)>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub no_call: Option<NoCall>,
}

#[derive(Debug, PartialEq, Eq, Serialize)]
pub struct UnsupportedAssay {
    pub id: String,
    pub name: String,
    pub mirrors: String,
    pub reason: String,
}

#[derive(Serialize)]
pub struct FeatureCard {
    #[serde(flatten)]
    pub note: FeatureNote,
    /// Presence over samples: base64 of little-endian u64 words.
    pub bits: String,
}

#[derive(Serialize)]
pub struct DatasetReport {
    pub schema: u32,
    /// Missing on legacy reports, which have not passed this validity audit.
    pub evaluation: EvaluationAudit,
    /// From `input_fingerprint`; missing on reports written before it.
    pub input_fingerprint: String,
    /// The cross-validation and learner settings behind `splits`, `mic`, and
    /// `models`. Missing on reports written before it was recorded.
    pub eval_config: EvalConfig,
    pub id: String,
    pub species: String,
    pub species_slug: String,
    pub antibiotic: String,
    pub drug_class: String,
    pub breakpoint: Option<Breakpoint>,
    pub counts: Counts,
    pub lineage_threshold: f32,
    pub catalog_groups: Vec<String>,
    pub samples: Vec<Sample>,
    /// The k-mer-identical duplicate groups (`kset_duplicates`) as sample
    /// indices. Each sits whole in one test fold of every entry in `splits`
    /// that ran; `mic` draws its own folds and may split one. Left out when
    /// there are none, so such a report is unchanged.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub kset_duplicates: Vec<Vec<usize>>,
    pub splits: Vec<SplitResult>,
    pub catalog: Option<CatalogBaseline>,
    pub panels: Vec<PanelResult>,
    pub unsupported_assays: Vec<UnsupportedAssay>,
    pub mic: Option<MicResult>,
    pub models: FinalModels,
    pub features: Vec<FeatureCard>,
    pub seconds: f64,
}

pub fn base64(words: &[u64]) -> String {
    const ALPHABET: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let bytes: Vec<u8> = words.iter().flat_map(|w| w.to_le_bytes()).collect();
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for chunk in bytes.chunks(3) {
        let b = [
            chunk[0],
            *chunk.get(1).unwrap_or(&0),
            *chunk.get(2).unwrap_or(&0),
        ];
        let n = ((b[0] as u32) << 16) | ((b[1] as u32) << 8) | b[2] as u32;
        for (i, shift) in [18, 12, 6, 0].into_iter().enumerate() {
            if i <= chunk.len() {
                out.push(ALPHABET[((n >> shift) & 63) as usize] as char);
            } else {
                out.push('=');
            }
        }
    }
    out
}

/// Feature ids the final models refer to, in a fixed order of importance.
fn model_features(m: &FinalModels) -> BTreeSet<u32> {
    let mut f = BTreeSet::new();
    for (rule, step) in m.scm.rules.iter().zip(&m.scm.steps) {
        f.insert(rule.feature);
        f.extend(step.equivalent.iter().take(5).map(|r| r.feature));
    }
    f.extend(
        m.cart
            .nodes
            .iter()
            .filter_map(|n| n.split.as_ref().map(|s| s.feature)),
    );
    if let Some(l) = &m.lasso {
        let mut c = l.coef.clone();
        c.sort_by(|a, b| b.1.abs().total_cmp(&a.1.abs()));
        f.extend(c.iter().take(30).map(|c| c.0));
    }
    f.extend(m.rf_importance.iter().take(30).map(|w| w.feature));
    f.extend(m.gbt_importance.iter().take(30).map(|w| w.feature));
    f
}

pub struct SpeciesInputs {
    pub matrix: FeatureMatrix,
    pub structure: Structure,
    pub hits: HashMap<String, Vec<FamilyHit>>,
    pub contigs: HashMap<String, u32>,
    pub phenotypes: Vec<Phenotype>,
    pub genomes: HashMap<String, GenomeRecord>,
    pub datasets: Vec<DatasetSpec>,
}

/// Reads what `grm learn` evaluates for a species, with the catalog of its
/// gene calls (see `species::catalog_of_calls`). Fails when the counts record
/// no kset digests, or when structure.json or matrix.json was not built from
/// counted.json as it is (`species::not_built_from`); `grm pipeline` recounts
/// and rebuilds those before it learns.
pub fn load_inputs(layout: &Layout, slug: &str) -> Result<(SpeciesInputs, Option<Catalog>)> {
    let t = Instant::now();
    let matrix = FeatureMatrix::load(&layout.matrix_dir(slug))?;
    let counted = species::load_counted(layout, slug)?;
    if let Some(g) = counted.undigested() {
        anyhow::bail!(
            "counted.json records no kset digest for {} (it predates them); rerun `grm count --species {slug}`",
            g.genome_id
        );
    }
    let structure: Structure = serde_json::from_slice(
        &std::fs::read(layout.species_dir(slug).join("structure.json"))
            .context("run `grm structure` first")?,
    )?;
    let (mut why, mut rerun) = (Vec::new(), Vec::new());
    for (stage, k, genomes, digests) in [
        (
            "structure",
            structure.k,
            &structure.genomes,
            &structure.kset_digests,
        ),
        (
            "matrix",
            matrix.summary.k,
            &matrix.summary.genomes,
            &matrix.summary.kset_digests,
        ),
    ] {
        if let Some(reason) = species::not_built_from(&counted, k, genomes, digests) {
            why.push(format!("{stage}.json {reason}"));
            rerun.push(format!("`grm {stage}`"));
        }
    }
    anyhow::ensure!(
        why.is_empty(),
        "{}; rerun {}",
        why.join("; "),
        rerun.join(" and ")
    );
    let catalog = species::catalog_of_calls(layout, &counted, &KmerSpace::new(matrix.summary.k))?;
    let datasets: Vec<DatasetSpec> =
        serde_json::from_slice(&std::fs::read(layout.curated().join("datasets.json"))?)?;
    let datasets: Vec<DatasetSpec> = datasets
        .into_iter()
        .filter(|d| d.species_slug == slug)
        .collect();
    let species = datasets
        .first()
        .map(|d| d.species.clone())
        .unwrap_or_default();
    let in_matrix: BTreeSet<&str> = matrix.summary.genomes.iter().map(String::as_str).collect();
    let phenotypes: Vec<Phenotype> =
        bvbrc::read_jsonl::<Phenotype>(&layout.curated().join("phenotypes.jsonl"))?
            .into_iter()
            .filter(|p| p.species == species && in_matrix.contains(p.genome_id.as_str()))
            .collect();
    let genomes: HashMap<String, GenomeRecord> =
        bvbrc::read_jsonl::<GenomeRecord>(&layout.genomes())?
            .into_iter()
            .filter(|g| in_matrix.contains(g.genome_id.as_str()))
            .map(|g| (g.genome_id.clone(), g))
            .collect();
    let hits = counted
        .genomes
        .iter()
        .map(|g| (g.genome_id.clone(), g.genes.clone()))
        .collect();
    let contigs = counted
        .genomes
        .iter()
        .map(|g| (g.genome_id.clone(), g.stats.contigs))
        .collect();
    eprintln!(
        "  {slug}: loaded {} patterns over {} genomes, {} phenotypes in {:.0}s",
        matrix.summary.n_patterns,
        matrix.summary.genomes.len(),
        phenotypes.len(),
        t.elapsed().as_secs_f64()
    );
    let input = SpeciesInputs {
        matrix,
        structure,
        hits,
        contigs,
        phenotypes,
        genomes,
        datasets,
    };
    Ok((input, catalog))
}

/// Probe set for `assay` from its verified wild-type sequences, with calls
/// for every genome in `genomes` order, or the reason the assay is not
/// scored; unreadable ksets are an error.
fn build_probe_set(
    layout: &Layout,
    genomes: &[String],
    k: usize,
    spec: &DatasetSpec,
    assay: &probes::Assay,
) -> Result<std::result::Result<probes::ProbeSet, String>> {
    let dir = layout.kset_dir(&spec.species_slug);
    let ksets: Vec<std::path::PathBuf> = genomes
        .iter()
        .map(|g| dir.join(format!("{g}.kset")))
        .collect();
    let space = KmerSpace::new(k);
    match probes::build_probes(assay, &ksets, &space) {
        Ok(set) => Ok(Ok(set)),
        Err(e) if e.is::<probes::KsetsUnreadable>() => {
            Err(e.context(format!("{}: probe assay {}", spec.id, assay.id)))
        }
        Err(e) => {
            eprintln!(
                "    {}: probe assay {} unavailable: {e:#}",
                spec.id, assay.id
            );
            Ok(Err(format!("{e:#}")))
        }
    }
}

fn probe_panel(
    assay: &probes::Assay,
    set: &probes::ProbeSet,
    columns: &[usize],
    truth: &[bool],
) -> PanelResult {
    let tally = |want: probes::ProbeCall| {
        columns
            .iter()
            .zip(truth)
            .filter(|&(&c, _)| set.calls[c] == want)
            .fold(
                (0, 0),
                |(r, s), (_, &t)| if t { (r + 1, s) } else { (r, s + 1) },
            )
    };
    let (called, calls): (Vec<bool>, Vec<f64>) = columns
        .iter()
        .zip(truth)
        .filter(|&(&c, _)| set.calls[c] != probes::ProbeCall::NoCall)
        .map(|(&c, &t)| {
            (
                t,
                (set.calls[c] == probes::ProbeCall::MutationDetected) as u8 as f64,
            )
        })
        .unzip();
    let (n_r, n_s) = tally(probes::ProbeCall::MutationDetected);
    let (x_r, x_s) = tally(probes::ProbeCall::NoCall);
    PanelResult {
        id: assay.id.into(),
        name: assay.name.into(),
        mirrors: assay.mirrors.into(),
        targets: set
            .regions
            .iter()
            .map(|r| {
                format!(
                    "{} {}..{} ({} bp, {} probe k-mers) on {}, {}",
                    r.gene,
                    r.start,
                    r.end,
                    r.end - r.start + 1,
                    r.probe_kmers,
                    r.record,
                    r.source
                )
            })
            .collect(),
        metrics: binary_metrics(&called, &calls, 0.5, false),
        detected: vec![("mutation detected".into(), n_r, n_s)],
        no_call: Some(NoCall {
            resistant: x_r,
            susceptible: x_s,
        }),
    }
}

type ProbeCache = HashMap<&'static str, std::result::Result<probes::ProbeSet, String>>;

fn probe_results(
    assays: &[&probes::Assay],
    probe_cache: &ProbeCache,
    columns: &[usize],
    truth: &[bool],
) -> (Vec<PanelResult>, Vec<UnsupportedAssay>) {
    let mut panels = Vec::new();
    let mut unsupported = Vec::new();
    for assay in assays {
        match &probe_cache[assay.id] {
            Ok(set) => panels.push(probe_panel(assay, set, columns, truth)),
            Err(reason) => unsupported.push(UnsupportedAssay {
                id: assay.id.into(),
                name: assay.name.into(),
                mirrors: assay.mirrors.into(),
                reason: reason.clone(),
            }),
        }
    }
    (panels, unsupported)
}

// Keep the public analysis API and its positional call contract unchanged.
#[allow(clippy::too_many_arguments)]
pub fn analyze_dataset(
    layout: &Layout,
    client: &Client,
    input: &SpeciesInputs,
    spec: &DatasetSpec,
    catalog: Option<&Catalog>,
    cfg: &EvalConfig,
    probe_cache: &mut ProbeCache,
    fingerprint: String,
) -> Result<Option<DatasetReport>> {
    let started = Instant::now();
    let Some(EvalInputs {
        cohort:
            Cohort {
                rows,
                columns,
                n_resistant: n_r,
                labels_from_mic_excluded,
                ..
            },
        resistant,
        min_support,
        ds,
        groups,
        duplicates,
    }) = eval_inputs(input, spec)?
    else {
        return Ok(None);
    };
    let n = rows.len();
    let mics: Vec<Option<Mic>> = rows.iter().map(|p| p.mic).collect();
    eprintln!(
        "    {}: {n} samples ({n_r} R), {} features",
        spec.id,
        ds.n_features()
    );

    // Region-probe assays (Xpert-style wild-type probes), built once per
    // species and assay, then read for this dataset's samples.
    let assays: Vec<&probes::Assay> = probes::for_species(&spec.species)
        .into_iter()
        .filter(|a| a.drugs.contains(&spec.antibiotic.as_str()))
        .collect();
    for assay in &assays {
        if !probe_cache.contains_key(assay.id) {
            let set = build_probe_set(
                layout,
                &input.matrix.summary.genomes,
                input.matrix.summary.k,
                spec,
                assay,
            )?;
            probe_cache.insert(assay.id, set);
        }
    }

    let mut out = eval::evaluate(
        &spec.id,
        &ds.x,
        &ds.y,
        &groups,
        &duplicates,
        min_support,
        Some(MicInput {
            mic: &mics,
            breakpoint: spec.breakpoint,
        }),
        cfg,
    )?;

    // Catalog baseline: the same calls in every split, since nothing trains.
    let truth = &resistant;
    let rule = catalog::drug_rule(&spec.antibiotic);
    let sample_hits: Vec<&[FamilyHit]> = rows
        .iter()
        .map(|p| {
            input
                .hits
                .get(&p.genome_id)
                .map(Vec::as_slice)
                .unwrap_or(&[])
        })
        .collect();
    let catalog_calls: Option<Vec<bool>> = rule.as_ref().map(|r| {
        sample_hits
            .iter()
            .map(|h| catalog::explain(r, h).0)
            .collect()
    });
    let catalog_baseline = match (&rule, &catalog_calls, catalog) {
        (Some(r), Some(calls), Some(cat)) => {
            let scores: Vec<f64> = calls.iter().map(|&c| c as u8 as f64).collect();
            let metrics = binary_metrics(truth, &scores, 0.5, false);
            for split in &mut out.splits {
                if split.unavailable_reason.is_some() {
                    continue;
                }
                split.learners.push(LearnerResult {
                    learner: "catalog".into(),
                    metrics: Some(metrics.clone()),
                    fold_balanced_accuracy: vec![],
                    roc: vec![],
                    scores: scores.clone(),
                    conformal: vec![],
                    fold_complexity: vec![],
                    seconds: 0.0,
                    no_call: None,
                    lasso_fits: None,
                });
            }
            let mut table: BTreeMap<String, GroupRow> = BTreeMap::new();
            for (h, &is_r) in sample_hits.iter().zip(truth) {
                for hit in h.iter() {
                    let row = table.entry(hit.group.clone()).or_insert_with(|| GroupRow {
                        group: hit.group.clone(),
                        class: hit.class.clone(),
                        subclass: hit.subclass.clone(),
                        counts_for_rule: catalog::explain(r, std::slice::from_ref(hit))
                            .1
                            .contains(&hit.group),
                        n_resistant: 0,
                        n_susceptible: 0,
                    });
                    if is_r {
                        row.n_resistant += 1
                    } else {
                        row.n_susceptible += 1
                    }
                }
            }
            let mut groups: Vec<GroupRow> = table.into_values().collect();
            groups.sort_by(|a, b| {
                (b.n_resistant + b.n_susceptible).cmp(&(a.n_resistant + a.n_susceptible))
            });
            Some(CatalogBaseline {
                version: cat.version.clone(),
                rule: r
                    .iter()
                    .map(|g| g.iter().map(|t| t.to_string()).collect())
                    .collect(),
                metrics,
                groups,
            })
        }
        _ => None,
    };

    let panels: Vec<PanelResult> = panels::for_drug(&spec.antibiotic)
        .into_iter()
        .map(|panel| {
            let found: Vec<Vec<String>> = sample_hits
                .iter()
                .map(|h| panels::detected(panel, h))
                .collect();
            let scores: Vec<f64> = found.iter().map(|f| (!f.is_empty()) as u8 as f64).collect();
            let mut detected: BTreeMap<String, (usize, usize)> = BTreeMap::new();
            for (f, &is_r) in found.iter().zip(truth) {
                for g in f {
                    let e = detected.entry(g.clone()).or_default();
                    if is_r { e.0 += 1 } else { e.1 += 1 }
                }
            }
            PanelResult {
                id: panel.id.into(),
                name: panel.name.into(),
                mirrors: panel.mirrors.into(),
                targets: panel.targets.iter().map(|t| t.to_string()).collect(),
                metrics: binary_metrics(truth, &scores, 0.5, false),
                detected: detected.into_iter().map(|(g, (r, s))| (g, r, s)).collect(),
                no_call: None,
            }
        })
        .collect();

    let mut panels = panels;
    let (probe_panels, unsupported_assays) = probe_results(&assays, probe_cache, &columns, truth);
    panels.extend(probe_panels);

    // Catalog group dictionary for the sample table.
    let mut dict: BTreeMap<String, u16> = BTreeMap::new();
    for h in &sample_hits {
        for hit in h.iter() {
            let next = dict.len() as u16;
            dict.entry(hit.group.clone()).or_insert(next);
        }
    }
    let samples: Vec<Sample> = rows
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let meta = input.genomes.get(&p.genome_id);
            Sample {
                id: p.genome_id.clone(),
                resistant: resistant[i],
                from_mic: p.label_from_mic,
                mic: p.mic,
                lineage: groups[i],
                country: meta.and_then(|m| m.isolation_country.clone()),
                year: meta.and_then(|m| m.collection_year.clone()),
                host: meta.and_then(|m| m.host_name.clone()),
                mlst: meta.and_then(|m| m.mlst.clone()),
                groups: sample_hits[i].iter().map(|h| dict[&h.group]).collect(),
                catalog_resistant: catalog_calls.as_ref().map(|c| c[i]),
            }
        })
        .collect();
    let mut catalog_groups = vec![String::new(); dict.len()];
    for (g, i) in dict {
        catalog_groups[i as usize] = g;
    }

    let wanted = model_features(&out.finals);
    let notes = annotate::annotate(
        layout,
        client,
        &input.matrix,
        &ds,
        &wanted,
        catalog,
        &input.contigs,
    )?;
    let features = notes
        .into_iter()
        .map(|note| {
            let bits = base64(ds.x.row(note.feature as usize));
            FeatureCard { note, bits }
        })
        .collect();

    let n_kmers: u64 = (0..ds.n_features())
        .flat_map(|f| {
            ds.patterns_of(f)
                .iter()
                .map(|&p| input.matrix.kmers_of(p as usize).len() as u64)
        })
        .sum();
    out.lineage_threshold = Some(input.structure.lineages.threshold as f64);
    Ok(Some(DatasetReport {
        schema: REPORT_SCHEMA,
        evaluation: evaluation_audit(&out.splits, out.finals.lasso_fit.as_ref()),
        input_fingerprint: fingerprint,
        eval_config: cfg.clone(),
        id: spec.id.clone(),
        species: spec.species.clone(),
        species_slug: spec.species_slug.clone(),
        antibiotic: spec.antibiotic.clone(),
        drug_class: spec.drug_class.clone(),
        breakpoint: spec.breakpoint,
        counts: Counts {
            samples: n,
            resistant: n_r,
            susceptible: n - n_r,
            labels_from_mic: rows.iter().filter(|p| p.label_from_mic).count(),
            labels_from_mic_excluded,
            features: ds.n_features(),
            patterns: ds.feature_patterns.len(),
            kmers: n_kmers,
            k: input.matrix.summary.k,
            lineages: groups.iter().collect::<BTreeSet<_>>().len(),
            min_support,
            min_class: spec.min_class,
        },
        lineage_threshold: input.structure.lineages.threshold,
        catalog_groups,
        samples,
        kset_duplicates: duplicates,
        splits: out.splits,
        catalog: catalog_baseline,
        panels,
        unsupported_assays,
        mic: out.mic,
        models: out.finals,
        features,
        seconds: started.elapsed().as_secs_f64(),
    }))
}

/// What `evaluate` reads for one dataset.
struct EvalInputs<'a> {
    cohort: Cohort<'a>,
    /// Each cohort row's reported call.
    resistant: Vec<bool>,
    min_support: u32,
    ds: Dataset,
    /// Each sample's lineage.
    groups: Vec<u32>,
    /// The cohort's k-mer-identical duplicate groups (`kset_duplicates`).
    duplicates: Vec<Vec<usize>>,
}

/// The dataset's cohort, its projection at the dataset's support floor, and
/// its lineages, as `analyze_dataset` evaluates them; None when a class falls
/// short of the dataset's class minimum.
fn eval_inputs<'a>(input: &'a SpeciesInputs, spec: &DatasetSpec) -> Result<Option<EvalInputs<'a>>> {
    let c = cohort(spec, &input.phenotypes, &input.matrix.summary.genomes);
    if !c.meets_min_class {
        return Ok(None);
    }
    let resistant: Vec<bool> = c.rows.iter().map(|p| p.reported == Some(Call::R)).collect();
    let min_support = (c.rows.len() as u32 / 200).max(3);
    // Per-fit candidates are a subset of the projection only when the species
    // matrix kept every k-mer that min_support could admit.
    anyhow::ensure!(
        input.matrix.summary.min_support <= min_support,
        "the species matrix was built with min_support {} above the dataset's {min_support}; rebuild it with `grm matrix --min-support {min_support}` or lower",
        input.matrix.summary.min_support
    );
    let ds = Dataset::project(&input.matrix, c.columns.clone(), &resistant, min_support);
    let groups = c
        .columns
        .iter()
        .map(|&col| input.structure.lineages.labels[col])
        .collect();
    let duplicates = kset_duplicates(&c.columns, &input.matrix.summary.kset_digests);
    Ok(Some(EvalInputs {
        cohort: c,
        resistant,
        min_support,
        ds,
        groups,
        duplicates,
    }))
}

/// The samples whose genomes have the same kset digest, that is the same
/// canonical k-mer set: one ascending list of sample indices per group of two
/// or more, ordered by first member. `columns` holds each sample's matrix
/// column and `digests` each column's digest.
fn kset_duplicates(columns: &[usize], digests: &[String]) -> Vec<Vec<usize>> {
    let mut by_digest: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, &col) in columns.iter().enumerate() {
        by_digest.entry(digests[col].as_str()).or_default().push(i);
    }
    let mut groups: Vec<Vec<usize>> = by_digest.into_values().filter(|g| g.len() > 1).collect();
    groups.sort_unstable();
    groups
}

/// Runs only the lasso fits of each named dataset's evaluation
/// (`eval::lasso_only`), on the inputs, cohort, projection, and lineages
/// `analyze_dataset` would evaluate, after reporting each dataset's inputs and
/// input fingerprint. With `fit` false it reports the inputs alone.
/// `progress` sees each penalty point as it is solved. It reads `layout` and
/// writes nothing. The large-species lasso scan (`scripts/lasso-scan`) runs
/// it.
pub fn lasso_scan(
    layout: &Layout,
    slug: &str,
    ids: &[String],
    cfg: &EvalConfig,
    fit: bool,
    observe: &mut dyn FnMut(&str, &LassoOnlyEvent) -> Result<()>,
    progress: &mut dyn FnMut(&str, &LassoProgress),
) -> Result<()> {
    let (input, catalog) = load_inputs(layout, slug)?;
    let digest = matrix_digest(&input.matrix);
    for id in ids {
        let started = Instant::now();
        let spec = input
            .datasets
            .iter()
            .find(|d| d.id == *id)
            .with_context(|| format!("{slug} has no curated dataset {id}"))?;
        let fingerprint = input_fingerprint(&input, spec, catalog.as_ref(), cfg, digest);
        let e = eval_inputs(&input, spec)?
            .with_context(|| format!("{id} is short of its class minimum and has no report"))?;
        observe(
            id,
            &LassoOnlyEvent::Inputs {
                samples: e.resistant.len(),
                resistant: e.cohort.n_resistant,
                features: e.ds.n_features(),
                lineages: e.groups.iter().collect::<BTreeSet<_>>().len(),
                min_support: e.min_support,
                input_fingerprint: fingerprint,
                seconds: started.elapsed().as_secs_f64(),
            },
        )?;
        if fit {
            eval::lasso_only(
                id,
                &e.ds.x,
                &e.ds.y,
                &e.groups,
                &e.duplicates,
                e.min_support,
                cfg,
                &mut |event| observe(id, event),
                &mut |point| progress(id, point),
            )?;
        }
    }
    Ok(())
}

/// One dataset's samples among the species' matrix genomes.
struct Cohort<'a> {
    /// The drug's phenotypes with an unconflicted reported R or S call
    /// (`Phenotype::reported_label`), in species-matrix column order.
    rows: Vec<&'a Phenotype>,
    /// Each row's matrix column.
    columns: Vec<usize>,
    n_resistant: usize,
    /// Matrix genomes whose label for the drug was inferred from MIC; no
    /// learner trains or is scored on them.
    labels_from_mic_excluded: usize,
    /// Both classes have at least `spec.min_class` rows. A dataset short of
    /// that gets no report.
    meets_min_class: bool,
}

/// Selects a dataset's rows from the phenotypes of the matrix genomes, whose
/// columns `genomes` lists. `analyze_dataset` evaluates these rows and
/// `input_fingerprint` hashes them, so the two always agree.
fn cohort<'a>(spec: &DatasetSpec, phenotypes: &'a [Phenotype], genomes: &[String]) -> Cohort<'a> {
    let column_of: HashMap<&str, usize> = genomes
        .iter()
        .enumerate()
        .map(|(i, g)| (g.as_str(), i))
        .collect();
    let drug: Vec<&Phenotype> = phenotypes
        .iter()
        .filter(|p| p.antibiotic == spec.antibiotic && column_of.contains_key(p.genome_id.as_str()))
        .collect();
    let mut rows: Vec<(&Phenotype, usize)> = drug
        .iter()
        .filter(|p| p.reported_label().is_some())
        .map(|&p| (p, column_of[p.genome_id.as_str()]))
        .collect();
    rows.sort_by_key(|&(_, c)| c);
    let (rows, columns): (Vec<&Phenotype>, Vec<usize>) = rows.into_iter().unzip();
    let n_resistant = rows.iter().filter(|p| p.reported == Some(Call::R)).count();
    Cohort {
        meets_min_class: n_resistant >= spec.min_class
            && rows.len() - n_resistant >= spec.min_class,
        labels_from_mic_excluded: drug.iter().filter(|p| p.label_from_mic).count(),
        rows,
        columns,
        n_resistant,
    }
}

/// Runs the species' datasets (or the named subset), largest minority class
/// first, skipping datasets whose saved report is still current (see
/// `stale_reason`) unless `force`.
pub fn analyze_species(
    layout: &Layout,
    slug: &str,
    only: &[String],
    force: bool,
    cfg: &EvalConfig,
) -> Result<Vec<String>> {
    let (input, catalog) = load_inputs(layout, slug)?;
    let digest = matrix_digest(&input.matrix);
    let client = Client::new();
    let dir = layout.results_dir(slug);
    std::fs::create_dir_all(&dir)?;
    let mut specs: Vec<&DatasetSpec> = input
        .datasets
        .iter()
        .filter(|d| only.is_empty() || only.iter().any(|o| *o == d.antibiotic || *o == d.id))
        .collect();
    specs.sort_by_key(|d| std::cmp::Reverse(d.n_resistant.min(d.n_susceptible)));
    let mut written = Vec::new();
    let mut probe_cache = HashMap::new();
    for spec in specs {
        let path = dir.join(format!("{}.json", spec.id));
        let fingerprint = input_fingerprint(&input, spec, catalog.as_ref(), cfg, digest);
        if path.exists() && !force {
            let existing: serde_json::Value = serde_json::from_slice(&std::fs::read(&path)?)?;
            match stale_reason(&existing, &fingerprint) {
                None => continue,
                Some(reason) => eprintln!("    {}: {reason}; recomputing report", spec.id),
            }
        }
        match analyze_dataset(
            layout,
            &client,
            &input,
            spec,
            catalog.as_ref(),
            cfg,
            &mut probe_cache,
            fingerprint,
        )? {
            Some(report) => {
                write_atomic(&path, &serde_json::to_vec(&report)?)?;
                eprintln!(
                    "{}",
                    done_line(&spec.id, report.seconds, &report.evaluation, &report.splits)?
                );
                written.push(spec.id.clone());
            }
            None => {
                eprintln!(
                    "    {}: too few labeled genomes after download, skipped",
                    spec.id
                );
                // `grm report` publishes any saved report for a curated
                // dataset, so one left from an earlier cohort has to go.
                if remove_if_present(&path).with_context(|| format!("remove {}", path.display()))? {
                    eprintln!("    {}: removed its earlier report", spec.id);
                }
            }
        }
    }
    Ok(written)
}

/// The log's done line for a report just written. A report whose evaluation
/// audit is invalid is an error instead: it stays on disk for review, the
/// species stops as it does on a lasso error, and `stale_reason` recomputes
/// it on resume rather than reusing it.
fn done_line(
    id: &str,
    seconds: f64,
    evaluation: &EvaluationAudit,
    splits: &[SplitResult],
) -> Result<String> {
    anyhow::ensure!(
        evaluation.valid,
        "{id}: report written, but its evaluation audit is not valid: {}",
        evaluation.reasons.join("; ")
    );
    Ok(format!(
        "    {id} done in {seconds:.0}s ({})",
        best_learner_words(splits)?
    ))
}

/// The log's words for the learner `grm report` will name best, picked by
/// the same rule from the saved splits.
fn best_learner_words(splits: &[SplitResult]) -> Result<String> {
    let rows = crate::report::learner_rows(&serde_json::to_value(splits)?);
    Ok(match crate::report::best_learner(&rows) {
        Some((best, ba)) => format!(
            "best random-split learner among those that call every isolate: {} {ba:.3}",
            best.learner
        ),
        None => "no learner calls every isolate".into(),
    })
}

/// Readers such as scheduled report rebuilds see the old file or the new
/// one, never a partial write.
fn write_atomic(path: &std::path::Path, bytes: &[u8]) -> Result<()> {
    let tmp = path.with_extension("json.tmp");
    std::fs::write(&tmp, bytes).with_context(|| format!("write {}", tmp.display()))?;
    std::fs::rename(&tmp, path).with_context(|| format!("replace {}", path.display()))
}

/// Bits helper for tests and callers that hold a raw row.
pub fn row_bits(x: &BitMatrix, f: usize) -> Bits {
    Bits::from_words(x.cols, x.row(f).to_vec())
}

/// Why a saved report has to be recomputed, or None when this evaluation
/// revision wrote it from the inputs `fingerprint` identifies, its
/// evaluation audit is valid, and every feature locus has its reference
/// annotation. An invalid report is never reused as a finished one.
fn stale_reason(report: &serde_json::Value, fingerprint: &str) -> Option<&'static str> {
    if report["evaluation"]["revision"].as_u64() != Some(EVALUATION_REVISION as u64) {
        return Some("evaluation audit revision changed");
    }
    // The fingerprint leaves annotation out, so a locus that records a failed
    // annotation download is what recomputes the report.
    let unannotated = report["features"].as_array().is_some_and(|features| {
        features
            .iter()
            .any(|f| f["locus"]["annotation_error"].is_string())
    });
    match report["input_fingerprint"].as_str() {
        None => Some("report has no input fingerprint"),
        Some(saved) if saved != fingerprint => Some("inputs changed since the report was written"),
        Some(_) if report["evaluation"]["valid"].as_bool() != Some(true) => {
            Some("saved report's evaluation audit is not valid")
        }
        Some(_) if unannotated => Some("a feature's reference annotation failed to download"),
        Some(_) => None,
    }
}

/// Identifies what a dataset's evaluation reads apart from code: the species
/// matrix (`matrix`, from `matrix_digest`, which covers every genome's kset
/// digest), the antibiotic, its class minimum and breakpoint, the lineage
/// threshold, the catalog calls version, every `cfg` setting, and, from
/// `cohort`, the count of MIC-inferred labels left out and each evaluation
/// row's genome, reported call, MIC, lineage number, and gene calls, in
/// column order. Genome metadata, contig counts, and the drug class change
/// only descriptive report fields, so they are left out. Code changes need a
/// new EVALUATION_REVISION.
fn input_fingerprint(
    input: &SpeciesInputs,
    spec: &DatasetSpec,
    catalog: Option<&Catalog>,
    cfg: &EvalConfig,
    matrix: u64,
) -> String {
    let mut h = Fnv::new();
    h.u64(matrix);
    h.str(&spec.antibiotic);
    h.u64(spec.min_class as u64);
    match spec.breakpoint {
        None => h.bytes(&[0]),
        Some(Breakpoint {
            r_above,
            agreement,
            n,
        }) => {
            h.bytes(&[1]);
            h.bytes(&r_above.to_le_bytes());
            h.u64(agreement.to_bits());
            h.u64(n as u64);
        }
    }
    h.bytes(&input.structure.lineages.threshold.to_bits().to_le_bytes());
    match catalog {
        None => h.bytes(&[0]),
        Some(c) => {
            h.bytes(&[1]);
            h.str(&c.calls_version());
        }
    }
    // Every setting is named, so a new one does not compile until it is
    // hashed here.
    let EvalConfig {
        folds,
        inner_folds,
        scm_max_rules,
        cart_max_depth,
        forest:
            ForestParams {
                n_trees,
                subsample,
                max_depth: forest_depth,
                min_leaf,
                seed: forest_seed,
            },
        gbt:
            GbtParams {
                rounds,
                learning_rate,
                max_depth: gbt_depth,
                lambda,
                min_child_weight,
                colsample,
                patience,
                seed: gbt_seed,
            },
        screened_features,
        lasso_lambdas,
        seed,
    } = cfg;
    for n in [
        *folds,
        *inner_folds,
        *scm_max_rules,
        *n_trees,
        *rounds,
        *patience,
        *screened_features,
        *lasso_lambdas,
    ] {
        h.u64(n as u64);
    }
    for n in [*cart_max_depth, *forest_depth, *min_leaf, *gbt_depth] {
        h.bytes(&n.to_le_bytes());
    }
    for x in [
        *subsample,
        *learning_rate,
        *lambda,
        *min_child_weight,
        *colsample,
    ] {
        h.u64(x.to_bits());
    }
    for n in [*seed, *forest_seed, *gbt_seed] {
        h.u64(n);
    }
    let Cohort {
        rows,
        columns,
        labels_from_mic_excluded,
        // Both follow from the rows and the class minimum hashed above.
        n_resistant: _,
        meets_min_class: _,
    } = cohort(spec, &input.phenotypes, &input.matrix.summary.genomes);
    h.u64(labels_from_mic_excluded as u64);
    for (p, c) in rows.into_iter().zip(columns) {
        h.str(&p.genome_id);
        h.bytes(&[(p.reported == Some(Call::R)) as u8]);
        match p.mic {
            None => h.bytes(&[0]),
            Some(Mic { log2, censor }) => {
                let censor: u8 = match censor {
                    Censor::Exact => 1,
                    Censor::AtMost => 2,
                    Censor::Above => 3,
                };
                h.bytes(&[censor]);
                h.bytes(&log2.to_le_bytes());
            }
        }
        h.bytes(&input.structure.lineages.labels[c].to_le_bytes());
        let hits = input.hits.get(&p.genome_id).map_or(&[][..], Vec::as_slice);
        h.u64(hits.len() as u64);
        for FamilyHit {
            family,
            group,
            allele,
            containment,
            class,
            subclass,
        } in hits
        {
            for text in [family, group, allele, class, subclass] {
                h.str(text);
            }
            h.bytes(&containment.to_bits().to_le_bytes());
        }
    }
    h.hex()
}

/// Digest of a loaded species matrix: its summary, with the kset digest
/// `grm count` recorded for each genome, then every pattern word, support
/// count, k-mer offset, and k-mer, read in place.
fn matrix_digest(matrix: &FeatureMatrix) -> u64 {
    let FeatureMatrix {
        summary,
        patterns,
        support,
        kmer_offsets,
        kmers,
    } = matrix;
    let MatrixSummary {
        k,
        genomes,
        min_support,
        n_patterns,
        kmers_distinct,
        kmers_kept,
        kmers_rare,
        kmers_core,
        kset_digests,
    } = summary;
    let mut h = Fnv::new();
    h.u64(*k as u64);
    h.u64(genomes.len() as u64);
    for g in genomes {
        h.str(g);
    }
    h.u64(kset_digests.len() as u64);
    for digest in kset_digests {
        h.str(digest);
    }
    h.bytes(&min_support.to_le_bytes());
    for n in [
        *n_patterns as u64,
        *kmers_distinct,
        *kmers_kept,
        *kmers_rare,
        *kmers_core,
    ] {
        h.u64(n);
    }
    for &w in &patterns.data {
        h.u64(w);
    }
    for s in support {
        h.bytes(&s.to_le_bytes());
    }
    for &o in kmer_offsets {
        h.u64(o);
    }
    for &x in kmers {
        h.u64(x);
    }
    h.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reports_are_reused_only_at_the_current_revision_and_fingerprint() {
        let key = "0123456789abcdef";
        let saved = |revision: u32, fingerprint: &str| {
            serde_json::json!({
                "evaluation": {"revision": revision, "valid": true},
                "input_fingerprint": fingerprint,
            })
        };
        assert_eq!(EVALUATION_REVISION, 10);
        assert_eq!(stale_reason(&saved(EVALUATION_REVISION, key), key), None);
        // Revision 8's solver could stop where revision 9's converges, and a
        // revision-9 report may have scored a genome with its own copy in
        // training. The input fingerprint does not see either change.
        for old in [7, 8, 9] {
            assert_eq!(
                stale_reason(&saved(old, key), key),
                Some("evaluation audit revision changed"),
                "a same-fingerprint revision-{old} result must not enter revision 10"
            );
        }
        // A current report whose audit is invalid, or says nothing, is not a
        // finished one: an unconverged point that was not selected must not
        // pass for done on resume.
        for valid in [serde_json::json!(false), serde_json::json!(null)] {
            let mut report = saved(EVALUATION_REVISION, key);
            report["evaluation"]["valid"] = valid;
            assert_eq!(
                stale_reason(&report, key),
                Some("saved report's evaluation audit is not valid")
            );
        }
        // Only the evaluation revision, its validity, and the inputs
        // recompute a report; a schema bump alone keeps a current schema-1
        // report.
        assert_eq!(
            stale_reason(
                &serde_json::json!({
                    "schema": 1,
                    "evaluation": {"revision": EVALUATION_REVISION, "valid": true},
                    "input_fingerprint": key,
                }),
                key
            ),
            None
        );
        assert_eq!(
            stale_reason(&saved(EVALUATION_REVISION, "fedcba9876543210"), key),
            Some("inputs changed since the report was written")
        );
        assert_eq!(
            stale_reason(
                &serde_json::json!({"evaluation": {"revision": EVALUATION_REVISION}}),
                key
            ),
            Some("report has no input fingerprint")
        );
        for old in [
            saved(EVALUATION_REVISION - 1, key),
            serde_json::json!({"evaluation": {"revision": 1}}),
            serde_json::json!({"schema": 1}),
        ] {
            assert_eq!(
                stale_reason(&old, key),
                Some("evaluation audit revision changed"),
                "{old}"
            );
        }
    }

    #[test]
    fn a_report_with_a_locus_missing_its_annotation_is_recomputed() {
        let key = "0123456789abcdef";
        let saved = |locus: serde_json::Value| {
            serde_json::json!({
                "evaluation": {"revision": EVALUATION_REVISION, "valid": true},
                "input_fingerprint": key,
                "features": [{"feature": 0, "locus": null}, {"feature": 1, "locus": locus}],
            })
        };
        let annotated = serde_json::json!({"genome_id": "g0", "gene": "blaX"});
        assert_eq!(stale_reason(&saved(annotated), key), None);
        let missing = serde_json::json!({
            "genome_id": "g0",
            "gene": null,
            "annotation_error": "genome_feature: giving up after 5 attempts",
        });
        assert_eq!(
            stale_reason(&saved(missing), key),
            Some("a feature's reference annotation failed to download")
        );
    }

    /// `phenotype` with the call replaced by a label inferred from MIC.
    fn inferred(genome: &str, antibiotic: &str, label: Call) -> Phenotype {
        Phenotype {
            reported: None,
            label_from_mic: true,
            ..phenotype(genome, antibiotic, label, None)
        }
    }

    fn at_min_class(min_class: usize) -> DatasetSpec {
        DatasetSpec {
            min_class,
            ..dataset("fixture")
        }
    }

    #[test]
    fn cohorts_gate_on_the_dataset_min_class_and_count_mic_labels_left_out() {
        let mut phenotypes: Vec<Phenotype> = (0..24)
            .map(|i| phenotype(&format!("r{i}"), "fixture", Call::R, None))
            .collect();
        phenotypes.extend((0..30).map(|i| phenotype(&format!("s{i}"), "fixture", Call::S, None)));
        phenotypes.extend((0..10).map(|i| inferred(&format!("m{i}"), "fixture", Call::R)));
        phenotypes.push(Phenotype {
            reported: None,
            label: None,
            conflict: true,
            ..phenotype("c0", "fixture", Call::R, None)
        });
        phenotypes.push(phenotype("i0", "fixture", Call::I, None));
        phenotypes.push(inferred("o0", "other", Call::S));
        // Matrix columns run opposite to the listing order.
        let genomes: Vec<String> = phenotypes
            .iter()
            .rev()
            .map(|p| p.genome_id.clone())
            .collect();

        assert!(!cohort(&at_min_class(25), &phenotypes, &genomes).meets_min_class);
        let c = cohort(&at_min_class(10), &phenotypes, &genomes);
        assert!(c.meets_min_class);
        assert_eq!((c.rows.len(), c.n_resistant), (54, 24));
        assert_eq!(c.labels_from_mic_excluded, 10);
        let labeled = phenotypes
            .iter()
            .filter(|p| p.antibiotic == "fixture" && matches!(p.label, Some(Call::R | Call::S)))
            .count();
        assert_eq!(c.rows.len() + c.labels_from_mic_excluded, labeled);
        assert_eq!(c.rows[0].genome_id, "s29");
        for (p, &column) in c.rows.iter().zip(&c.columns) {
            assert_eq!(genomes[column], p.genome_id);
        }
        assert!(c.columns.windows(2).all(|w| w[0] < w[1]));
    }

    #[test]
    fn cohorts_admit_a_class_of_exactly_min_class() {
        let mut phenotypes: Vec<Phenotype> = (0..24)
            .map(|i| phenotype(&format!("r{i}"), "fixture", Call::R, None))
            .collect();
        phenotypes.extend((0..24).map(|i| phenotype(&format!("s{i}"), "fixture", Call::S, None)));
        let genomes: Vec<String> = phenotypes.iter().map(|p| p.genome_id.clone()).collect();

        let c = cohort(&at_min_class(24), &phenotypes, &genomes);
        assert!(c.meets_min_class);
        assert_eq!((c.rows.len(), c.n_resistant), (48, 24));
        assert!(!cohort(&at_min_class(25), &phenotypes, &genomes).meets_min_class);
    }

    #[test]
    fn audit_limits_validity_to_internal_cross_validation() {
        let convergence = crate::learn::lasso::Convergence {
            irls_steps: 3,
            sweep_cap_hits: 0,
            converged: true,
        };
        let clean = LassoFitAudit {
            path_points: 1,
            path_unconverged: 0,
            selected_point: Some(0),
            selected_convergence: Some(convergence),
            refit_convergence: Some(convergence),
        };
        let report = serde_json::to_value(evaluation_audit(&[], Some(&clean))).unwrap();
        assert_eq!(report["revision"], EVALUATION_REVISION);
        assert_eq!(report["valid"], true);
        assert_eq!(report["reasons"], serde_json::json!([]));
        assert!(
            report["scope"]
                .as_str()
                .unwrap()
                .contains("not external or clinical validation")
        );

        for bad in [
            LassoFitAudit {
                path_unconverged: 1,
                ..clean.clone()
            },
            LassoFitAudit {
                selected_convergence: Some(crate::learn::lasso::Convergence {
                    converged: false,
                    ..convergence
                }),
                ..clean.clone()
            },
            LassoFitAudit {
                refit_convergence: Some(crate::learn::lasso::Convergence {
                    converged: false,
                    ..convergence
                }),
                ..clean.clone()
            },
        ] {
            let audit = evaluation_audit(&[], Some(&bad));
            assert!(!audit.valid);
            assert!(!audit.reasons.is_empty());
        }
        assert!(!evaluation_audit(&[], None).valid);

        let split = SplitResult {
            scheme: "random".into(),
            folds: 1,
            n_groups: 2,
            test_fold: vec![],
            learners: vec![LearnerResult {
                learner: "lasso".into(),
                metrics: None,
                fold_balanced_accuracy: vec![],
                roc: vec![],
                scores: vec![],
                conformal: vec![],
                fold_complexity: vec![0.0],
                seconds: 0.0,
                no_call: None,
                lasso_fits: Some(vec![LassoFitAudit {
                    path_unconverged: 1,
                    ..clean.clone()
                }]),
            }],
            unavailable_reason: None,
        };
        assert!(!evaluation_audit(&[split], Some(&clean)).valid);
    }

    fn clean_lasso_fit() -> LassoFitAudit {
        let convergence = crate::learn::lasso::Convergence {
            irls_steps: 3,
            sweep_cap_hits: 0,
            converged: true,
        };
        LassoFitAudit {
            path_points: 4,
            path_unconverged: 0,
            selected_point: Some(2),
            selected_convergence: Some(convergence),
            refit_convergence: Some(convergence),
        }
    }

    /// A two-fold random split whose lasso row holds `fits`.
    fn lasso_split(fits: Option<Vec<LassoFitAudit>>) -> SplitResult {
        SplitResult {
            scheme: "random".into(),
            folds: 2,
            n_groups: 4,
            test_fold: vec![],
            learners: vec![LearnerResult {
                learner: "lasso".into(),
                metrics: None,
                fold_balanced_accuracy: vec![],
                roc: vec![],
                scores: vec![],
                conformal: vec![],
                fold_complexity: vec![0.0, 0.0],
                seconds: 0.0,
                no_call: None,
                lasso_fits: fits,
            }],
            unavailable_reason: None,
        }
    }

    #[test]
    fn every_available_split_needs_one_clean_lasso_record_per_fold() {
        let clean = clean_lasso_fit();
        let audit = |split| evaluation_audit(&[split], Some(&clean));
        let valid = audit(lasso_split(Some(vec![clean.clone(); 2])));
        assert!(valid.valid, "{:?}", valid.reasons);

        let missing = audit(lasso_split(None));
        assert!(!missing.valid);
        assert_eq!(
            missing.reasons,
            ["random split has no per-fold lasso solver status"]
        );
        for records in [1, 3] {
            let wrong = audit(lasso_split(Some(vec![clean.clone(); records])));
            assert!(!wrong.valid);
            assert_eq!(
                wrong.reasons,
                [format!(
                    "random split has {records} lasso solver records for 2 folds"
                )]
            );
        }
        // A split with no lasso row at all has no status either.
        let mut no_row = lasso_split(None);
        no_row.learners.clear();
        assert_eq!(
            audit(no_row).reasons,
            ["random split has no per-fold lasso solver status"]
        );
        // An unavailable split fits nothing and needs no record.
        let mut unavailable = lasso_split(None);
        unavailable.learners.clear();
        unavailable.unavailable_reason = Some("too few lineages".into());
        assert!(audit(unavailable).valid);
    }

    #[test]
    fn an_invalid_report_is_kept_but_never_logged_done() {
        let id = "fixture__meropenem";
        let clean = clean_lasso_fit();
        let line = done_line(id, 12.4, &evaluation_audit(&[], Some(&clean)), &[]).unwrap();
        assert_eq!(
            line,
            "    fixture__meropenem done in 12s (no learner calls every isolate)"
        );

        // An unconverged point that was not selected leaves a written report
        // whose audit is invalid: an error naming the dataset and why, and no
        // line the result checker would read as done.
        let excluded = LassoFitAudit {
            path_unconverged: 1,
            ..clean
        };
        let audit = evaluation_audit(&[], Some(&excluded));
        let error = done_line(id, 12.4, &audit, &[]).unwrap_err().to_string();
        assert!(
            error.starts_with("fixture__meropenem: report written"),
            "{error}"
        );
        assert!(
            error.contains("final model has an unconverged two-class lasso fit"),
            "{error}"
        );
        assert!(!error.contains(" done in "), "{error}");
        let saved = serde_json::json!({
            "evaluation": serde_json::to_value(&audit).unwrap(),
            "input_fingerprint": "0123456789abcdef",
        });
        assert_eq!(
            stale_reason(&saved, "0123456789abcdef"),
            Some("saved report's evaluation audit is not valid")
        );
    }

    #[test]
    fn the_log_names_the_best_learner_that_calls_every_isolate() {
        let row = |learner: &str, ba: f64, uncalled: usize| LearnerResult {
            learner: learner.into(),
            metrics: Some(BinaryMetrics {
                balanced_accuracy: ba,
                ..Default::default()
            }),
            fold_balanced_accuracy: vec![],
            roc: vec![],
            scores: vec![],
            conformal: vec![],
            fold_complexity: vec![],
            seconds: 0.0,
            no_call: Some(crate::learn::NoCall {
                resistant: uncalled,
                susceptible: 0,
            }),
            lasso_fits: None,
        };
        let split = |scheme: &str, learners| SplitResult {
            scheme: scheme.into(),
            folds: 5,
            n_groups: 12,
            test_fold: vec![],
            learners,
            unavailable_reason: None,
        };
        // The lasso scores higher but leaves two isolates uncalled.
        let splits = [
            split("random", vec![row("rf", 0.8124, 0), row("lasso", 0.9, 2)]),
            split("lineage", vec![row("rf", 0.7, 0), row("lasso", 0.6, 0)]),
        ];
        assert_eq!(
            best_learner_words(&splits).unwrap(),
            "best random-split learner among those that call every isolate: rf 0.812"
        );
        let partial = [split("random", vec![row("rf", 0.8, 1)])];
        assert_eq!(
            best_learner_words(&partial).unwrap(),
            "no learner calls every isolate"
        );
    }

    #[test]
    fn report_writes_replace_the_whole_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("d.json");
        std::fs::write(&path, b"old report, longer than the new one").unwrap();
        write_atomic(&path, b"{}").unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"{}");
        assert!(!dir.path().join("d.json.tmp").exists());
    }

    #[test]
    fn kset_duplicates_group_samples_by_digest_in_sample_order() {
        let digests: Vec<String> = ["a", "b", "a", "c", "b", "a", "d"]
            .map(String::from)
            .to_vec();
        // Samples 0..7 sit in columns 6, 0, 2, 4, 5, 1, 3: digests d a a b a b c.
        assert_eq!(
            kset_duplicates(&[6, 0, 2, 4, 5, 1, 3], &digests),
            [vec![1, 2, 4], vec![3, 5]]
        );
        assert!(kset_duplicates(&[0, 1, 3, 6], &digests).is_empty());
    }

    #[test]
    fn base64_matches_known_encoding() {
        // 0x0000000000616263 little endian is "cba\0\0\0\0\0".
        assert_eq!(base64(&[0x0000_0000_0061_6263]), "Y2JhAAAAAAA=");
    }

    use crate::probes::ProbeCall::{MutationDetected, NoCall, WildType};

    const WILD: &str = "GATTACAGCTATGAAACGTCTGGCAGTTCCGATCAGCTAATGCATCGGAACTTCAGGTCA";
    const MUTANT: &str = "GATTACAGCTATGAAACTTCTGGCAGTTCCGATCAGCTAATGCATCGGAACTTCAGGTCA";
    const ABSENT: &str = "CCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACC";

    const PROBED: probes::Assay = probes::Assay {
        id: "fixture_gyra",
        name: "fixture probes",
        mirrors: "fixture",
        species: "Fixture species",
        drugs: &["fixture"],
        regions: &[probes::Region {
            gene: "gyrA",
            span: probes::Span::Codons(3, 4),
            wild_type: Some(probes::WildType {
                source: "fixture wild type",
                record: "c1",
                cds_start: 11,
                cds_end: 40,
                strand: "+",
                cds: "ATGAAACGTCTGGCAGTTCCGATCAGCTAA",
                residues: "RL",
            }),
        }],
        unsupported: None,
    };

    const UNSOURCED: probes::Assay = probes::Assay {
        id: "fixture_katg",
        regions: &[probes::Region {
            gene: "katG",
            span: probes::Span::Codons(3, 4),
            wild_type: None,
        }],
        ..PROBED
    };

    const UNSUPPORTED: probes::Assay = probes::Assay {
        id: "fixture_qrdr",
        name: "fixture unsupported probes",
        unsupported: Some("fixture synonymous changes read as mutations"),
        ..PROBED
    };

    fn write_kset(layout: &Layout, id: &str, seq: &str) {
        write_kset_at(layout, id, seq, 7);
    }

    fn write_kset_at(layout: &Layout, id: &str, seq: &str, k: usize) {
        let space = KmerSpace::new(k);
        let dir = layout.kset_dir("fixture");
        std::fs::create_dir_all(&dir).unwrap();
        let fna = layout.root.join(format!("{id}.fna"));
        std::fs::write(&fna, format!(">c1\n{seq}\n")).unwrap();
        let (values, _) = crate::kset::count_file(&fna, &space, 1).unwrap();
        crate::kset::write(&dir.join(format!("{id}.kset")), id, &space, &values).unwrap();
    }

    fn probe_set(
        layout: &Layout,
        genomes: &[&str],
        assay: &probes::Assay,
    ) -> Result<std::result::Result<probes::ProbeSet, String>> {
        let genomes: Vec<String> = genomes.iter().map(|g| g.to_string()).collect();
        build_probe_set(layout, &genomes, 7, &dataset("fixture"), assay)
    }

    #[test]
    fn probe_assays_fail_instead_of_dropping_their_panel_without_ksets() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let cohort = ["g1", "r", "m", "g2"];
        let err = probe_set(&layout, &cohort, &PROBED).unwrap_err();
        assert!(err.is::<probes::KsetsUnreadable>(), "{err:#}");
        for id in ["g1", "r", "m"] {
            write_kset(&layout, id, WILD);
        }
        let err = probe_set(&layout, &cohort, &PROBED).unwrap_err();
        assert!(err.is::<probes::KsetsUnreadable>(), "{err:#}");
        assert!(format!("{err:#}").contains("g2.kset"), "{err:#}");
    }

    #[test]
    fn probe_assays_fail_on_a_kset_counted_at_another_k() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        for id in ["g1", "r"] {
            write_kset(&layout, id, WILD);
        }
        // A recount at k = 8 replaced this genome's kset after the matrix.
        write_kset_at(&layout, "m", WILD, 8);
        let cohort = ["g1", "r", "m"];
        let err = probe_set(&layout, &cohort, &PROBED).unwrap_err();
        assert!(err.is::<probes::KsetsUnreadable>(), "{err:#}");
        assert!(
            format!("{err:#}").contains("m.kset holds 8-mers, not the 7-mers"),
            "{err:#}"
        );
        let reason = probe_set(&layout, &cohort, &UNSUPPORTED)
            .unwrap()
            .unwrap_err();
        assert_eq!(reason, "fixture synonymous changes read as mutations");
    }

    #[test]
    fn unsourced_and_unsupported_probe_assays_give_a_reason_without_ksets() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let cohort = ["g1", "r", "m", "g2"];
        let reason = probe_set(&layout, &cohort, &UNSOURCED)
            .unwrap()
            .unwrap_err();
        assert!(
            reason.contains("katG has no verified wild-type sequence"),
            "{reason}"
        );
        let reason = probe_set(&layout, &cohort, &UNSUPPORTED)
            .unwrap()
            .unwrap_err();
        assert_eq!(reason, "fixture synonymous changes read as mutations");
    }

    #[test]
    fn cached_ksets_give_probe_calls_in_genome_order() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        for id in ["g1", "r", "g2"] {
            write_kset(&layout, id, WILD);
        }
        write_kset(&layout, "m", MUTANT);
        write_kset(&layout, "x", ABSENT);
        let cohort = ["g1", "r", "m", "x", "g2"];
        let set = probe_set(&layout, &cohort, &PROBED).unwrap().unwrap();
        assert_eq!(
            set.calls,
            [WildType, WildType, MutationDetected, NoCall, WildType]
        );
        assert!(probe_set(&layout, &cohort, &UNSOURCED).unwrap().is_err());
        assert!(probe_set(&layout, &cohort, &UNSUPPORTED).unwrap().is_err());
    }

    #[test]
    fn mutant_majority_cohorts_keep_wild_type_calls() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        for id in ["g1", "g2"] {
            write_kset(&layout, id, WILD);
        }
        for id in ["r", "m1", "m2", "m3"] {
            write_kset(&layout, id, MUTANT);
        }
        let alone = probe_set(&layout, &["g1", "m1"], &PROBED).unwrap().unwrap();
        assert_eq!(alone.calls, [WildType, MutationDetected]);
        let cohort = ["r", "m3", "g2", "m2", "m1", "g1"];
        let set = probe_set(&layout, &cohort, &PROBED).unwrap().unwrap();
        assert_eq!(
            set.calls,
            [
                MutationDetected,
                MutationDetected,
                WildType,
                MutationDetected,
                MutationDetected,
                WildType
            ]
        );
    }

    fn fixture_set() -> probes::ProbeSet {
        probes::ProbeSet {
            assay: PROBED.id.into(),
            regions: vec![probes::ResolvedRegion {
                gene: "gyrA".into(),
                record: "c1".into(),
                source: "fixture wild type".into(),
                start: 17,
                end: 22,
                strand: "+".into(),
                sequence: "CGTCTG".into(),
                probe_kmers: 12,
            }],
            calls: vec![
                MutationDetected,
                WildType,
                WildType,
                NoCall,
                MutationDetected,
                NoCall,
                MutationDetected,
            ],
        }
    }

    #[test]
    fn probe_panels_score_called_samples_and_count_no_calls_apart() {
        let truth = [true, true, true, true, true, false, false];
        let panel = probe_panel(&PROBED, &fixture_set(), &[6, 5, 4, 3, 2, 1, 0], &truth);
        let m = &panel.metrics;
        assert_eq!((m.n, m.tp, m.fn_, m.tn, m.fp), (5, 2, 1, 1, 1));
        assert_eq!(panel.detected, [("mutation detected".to_string(), 2, 1)]);
        assert_eq!(
            panel.no_call,
            Some(crate::learn::NoCall {
                resistant: 2,
                susceptible: 0
            })
        );
        assert_eq!(
            panel.targets,
            ["gyrA 17..22 (6 bp, 12 probe k-mers) on c1, fixture wild type"]
        );
    }

    #[test]
    fn unsupported_assays_are_listed_with_their_reason_instead_of_a_panel() {
        let reason = "fixture synonymous changes read as mutations";
        let cache: ProbeCache = [
            (PROBED.id, Ok(fixture_set())),
            (UNSUPPORTED.id, Err(reason.to_string())),
        ]
        .into_iter()
        .collect();
        let truth = [true, true, true, true, true, false, false];
        let (panels, unsupported) = probe_results(
            &[&PROBED, &UNSUPPORTED],
            &cache,
            &[6, 5, 4, 3, 2, 1, 0],
            &truth,
        );
        assert_eq!(panels.len(), 1);
        assert_eq!(panels[0].id, PROBED.id);
        assert_eq!(
            unsupported,
            [UnsupportedAssay {
                id: "fixture_qrdr".into(),
                name: "fixture unsupported probes".into(),
                mirrors: "fixture".into(),
                reason: reason.into(),
            }]
        );
        let json = serde_json::to_value(&panels[0]).unwrap();
        assert_eq!(
            json["no_call"],
            serde_json::json!({"resistant": 2, "susceptible": 0})
        );
    }

    fn hit(group: &str) -> FamilyHit {
        FamilyHit {
            family: group.into(),
            group: group.into(),
            allele: format!("{group}-1"),
            containment: 1.0,
            class: "BETA-LACTAM".into(),
            subclass: "CARBAPENEM".into(),
        }
    }

    fn phenotype(genome: &str, antibiotic: &str, reported: Call, mic: Option<Mic>) -> Phenotype {
        Phenotype {
            genome_id: genome.into(),
            species: "Fixture species".into(),
            antibiotic: antibiotic.into(),
            reported: Some(reported),
            conflict: false,
            mic,
            mic_spread: 0,
            label: Some(reported),
            label_from_mic: false,
            methods: vec![],
            standards: vec![],
            n_records: 1,
        }
    }

    fn dataset(antibiotic: &str) -> DatasetSpec {
        DatasetSpec {
            id: format!("fixture__{antibiotic}"),
            species: "Fixture species".into(),
            species_slug: "fixture".into(),
            antibiotic: antibiotic.into(),
            drug_class: "fixture class".into(),
            n_resistant: 2,
            n_susceptible: 2,
            n_intermediate: 0,
            n_with_mic: 2,
            n_labels_from_mic: 0,
            breakpoint: Some(Breakpoint {
                r_above: 2,
                agreement: 0.97,
                n: 60,
            }),
            min_class: 25,
        }
    }

    /// Six genomes. Meropenem evaluates g0 R, g1 S, g3 R, and g5 S, listed
    /// out of column order. It skips g2, whose label is marked as inferred
    /// from its MIC, and g4, whose records conflict.
    fn inputs() -> SpeciesInputs {
        let genomes: Vec<String> = (0..6).map(|i| format!("g{i}")).collect();
        let mic = |log2, censor| Some(Mic { log2, censor });
        let mut inferred = phenotype("g2", "meropenem", Call::R, mic(6, Censor::Exact));
        inferred.label_from_mic = true;
        let mut conflicted = phenotype("g4", "meropenem", Call::S, mic(0, Censor::Exact));
        conflicted.conflict = true;
        let kset_digests: Vec<String> = genomes.iter().map(|g| format!("kset of {g}")).collect();
        SpeciesInputs {
            matrix: FeatureMatrix {
                summary: MatrixSummary {
                    k: 7,
                    genomes: genomes.clone(),
                    min_support: 1,
                    n_patterns: 2,
                    kmers_distinct: 9,
                    kmers_kept: 3,
                    kmers_rare: 4,
                    kmers_core: 2,
                    kset_digests: kset_digests.clone(),
                },
                patterns: BitMatrix {
                    rows: 2,
                    cols: 6,
                    words: 1,
                    data: vec![0b001001, 0b010110],
                },
                support: vec![2, 3],
                kmer_offsets: vec![0, 2, 3],
                kmers: vec![0x11, 0x12, 0x13],
            },
            structure: Structure {
                genomes: genomes.clone(),
                k: 7,
                kset_digests,
                lineages: crate::cluster::Lineages {
                    threshold: 0.005,
                    labels: vec![0, 0, 1, 1, 2, 2],
                    n_clusters: 3,
                    largest_share: 1.0 / 3.0,
                },
                merges: vec![],
                leaf_order: vec![],
                coords: vec![],
                nearest: vec![],
                sketch_size: species::SKETCH,
            },
            hits: [
                ("g0", vec![hit("blaKPC")]),
                ("g3", vec![hit("blaKPC"), hit("blaOXA")]),
            ]
            .into_iter()
            .map(|(g, h)| (g.to_string(), h))
            .collect(),
            contigs: genomes.iter().map(|g| (g.clone(), 40)).collect(),
            phenotypes: vec![
                phenotype("g3", "meropenem", Call::R, mic(4, Censor::Above)),
                phenotype("g0", "meropenem", Call::R, mic(3, Censor::Exact)),
                phenotype("g5", "meropenem", Call::S, mic(-1, Censor::AtMost)),
                phenotype("g1", "meropenem", Call::S, None),
                inferred,
                conflicted,
                phenotype("g1", "ciprofloxacin", Call::R, None),
                phenotype("g2", "ciprofloxacin", Call::S, mic(-2, Censor::Exact)),
            ],
            genomes: genomes
                .iter()
                .map(|g| {
                    let record = GenomeRecord {
                        genome_id: g.clone(),
                        isolation_country: Some("Fixture".into()),
                        collection_year: Some("2020".into()),
                        ..Default::default()
                    };
                    (g.clone(), record)
                })
                .collect(),
            datasets: vec![dataset("meropenem"), dataset("ciprofloxacin")],
        }
    }

    /// The reuse key `analyze_species` computes for `spec` with the default
    /// evaluation settings.
    fn key(input: &SpeciesInputs, spec: &DatasetSpec, catalog: Option<&Catalog>) -> String {
        input_fingerprint(
            input,
            spec,
            catalog,
            &EvalConfig::default(),
            matrix_digest(&input.matrix),
        )
    }

    /// The same genomes, lineage numbers, and presence bits in reversed
    /// matrix columns.
    fn reverse_columns(x: &mut SpeciesInputs) {
        let n = x.matrix.summary.genomes.len() as u32;
        x.matrix.summary.genomes.reverse();
        x.matrix.summary.kset_digests.reverse();
        x.structure.genomes.reverse();
        x.structure.kset_digests.reverse();
        x.structure.lineages.labels.reverse();
        x.structure.coords.reverse();
        for w in &mut x.matrix.patterns.data {
            *w = w.reverse_bits() >> (64 - n);
        }
    }

    fn row<'a>(input: &'a mut SpeciesInputs, genome: &str, antibiotic: &str) -> &'a mut Phenotype {
        input
            .phenotypes
            .iter_mut()
            .find(|p| p.genome_id == genome && p.antibiotic == antibiotic)
            .unwrap()
    }

    type Edit = fn(&mut SpeciesInputs, &mut DatasetSpec);

    fn meropenem_key_after(edit: Edit) -> String {
        let (mut input, mut spec) = (inputs(), dataset("meropenem"));
        edit(&mut input, &mut spec);
        key(&input, &spec, None)
    }

    #[test]
    fn fingerprint_changes_with_every_input_the_evaluation_reads() {
        let base = key(&inputs(), &dataset("meropenem"), None);
        assert_eq!(meropenem_key_after(|_, _| ()), base);
        let changes: &[(&str, Edit)] = &[
            ("reported call", |x, _| {
                row(x, "g3", "meropenem").reported = Some(Call::S)
            }),
            ("MIC value", |x, _| {
                row(x, "g0", "meropenem").mic.as_mut().unwrap().log2 = 4
            }),
            ("MIC censoring", |x, _| {
                row(x, "g0", "meropenem").mic.as_mut().unwrap().censor = Censor::Above
            }),
            ("MIC added", |x, _| {
                row(x, "g1", "meropenem").mic = Some(Mic {
                    log2: 0,
                    censor: Censor::Exact,
                })
            }),
            ("genome dropped", |x, _| {
                x.phenotypes.retain(|p| p.genome_id != "g5")
            }),
            ("evaluated genome swapped", |x, _| {
                // g4 takes g5's call and MIC; both share a lineage and have no
                // gene calls, so only the genome, and its matrix column, change.
                let g5 = row(x, "g5", "meropenem").clone();
                *row(x, "g4", "meropenem") = Phenotype {
                    genome_id: "g4".into(),
                    ..g5
                };
                row(x, "g5", "meropenem").conflict = true;
            }),
            ("genome order", |x, _| reverse_columns(x)),
            ("lineage label", |x, _| {
                // g5 is evaluation row 3 and matrix column 5.
                x.structure.lineages.labels[5] = 0
            }),
            ("lineage threshold", |x, _| {
                x.structure.lineages.threshold = 0.0025
            }),
            ("breakpoint", |_, s| {
                s.breakpoint.as_mut().unwrap().r_above = 3
            }),
            ("breakpoint removed", |_, s| s.breakpoint = None),
            ("antibiotic", |x, s| {
                // The same rows under another drug name.
                s.antibiotic = "imipenem".into();
                for p in &mut x.phenotypes {
                    if p.antibiotic == "meropenem" {
                        p.antibiotic = "imipenem".into();
                    }
                }
            }),
            ("matrix bit", |x, _| x.matrix.patterns.data[1] ^= 1),
            ("pattern support", |x, _| x.matrix.support[0] = 3),
            ("k-mer offset", |x, _| x.matrix.kmer_offsets[1] = 1),
            ("k-mer", |x, _| x.matrix.kmers[2] ^= 1),
            ("k", |x, _| x.matrix.summary.k = 31),
            ("matrix min_support", |x, _| {
                x.matrix.summary.min_support = 2
            }),
            ("matrix counts", |x, _| x.matrix.summary.kmers_rare = 5),
            ("gene call", |x, _| {
                x.hits.get_mut("g3").unwrap()[1].allele = "blaOXA-2".into()
            }),
            ("gene containment", |x, _| {
                x.hits.get_mut("g0").unwrap()[0].containment = 0.9
            }),
            ("gene call added", |x, _| {
                x.hits.insert("g1".into(), vec![hit("blaNDM")]);
            }),
            ("characters moved between gene-call fields", |x, _| {
                let h = &mut x.hits.get_mut("g0").unwrap()[0];
                h.family = "blaKP".into();
                h.group = "CblaKPC".into();
            }),
            ("kset content", |x, _| {
                // A recount of g3, column 3, with as many k-mers as before
                // but other ones.
                x.matrix.summary.kset_digests[3] = "another kset of g3".into()
            }),
            ("kset of a genome no row evaluates", |x, _| {
                // g2 is in the species matrix, which every dataset projects.
                x.matrix.summary.kset_digests[2] = "another kset of g2".into()
            }),
        ];
        for (what, edit) in changes {
            assert_ne!(meropenem_key_after(*edit), base, "{what}");
        }

        let dir = tempfile::tempdir().unwrap();
        let mut catalog = write_catalog(&Layout::new(dir.path()), "2026-01-01.1");
        let with_catalog = key(&inputs(), &dataset("meropenem"), Some(&catalog));
        assert_ne!(with_catalog, base, "catalog loaded");
        catalog.version = "2026-02-01.1".into();
        assert_ne!(
            key(&inputs(), &dataset("meropenem"), Some(&catalog)),
            with_catalog,
            "catalog version"
        );
    }

    #[test]
    fn fingerprint_tracks_the_dataset_class_minimum() {
        let base = key(&inputs(), &dataset("meropenem"), None);
        assert_ne!(
            meropenem_key_after(|_, spec| spec.min_class = 24),
            base,
            "min_class 25 to 24"
        );
    }

    #[test]
    fn fingerprint_tracks_mic_inferred_labels_added_on_a_matrix_genome() {
        let mut input = inputs();
        let id = "g6".to_string();
        input.matrix.summary.genomes.push(id.clone());
        input
            .matrix
            .summary
            .kset_digests
            .push(format!("kset of {id}"));
        input.matrix.patterns.cols += 1;
        input.structure.genomes.push(id.clone());
        input.structure.kset_digests.push(format!("kset of {id}"));
        input.structure.lineages.labels.push(0);
        let before = key(&input, &dataset("meropenem"), None);
        input.phenotypes.push(inferred(&id, "meropenem", Call::R));
        assert_ne!(key(&input, &dataset("meropenem"), None), before);
    }

    #[test]
    fn fingerprint_tracks_every_evaluation_config_field() {
        fn changed(edit: fn(&mut EvalConfig)) -> String {
            let input = inputs();
            let spec = dataset("meropenem");
            let mut cfg = EvalConfig::default();
            edit(&mut cfg);
            input_fingerprint(&input, &spec, None, &cfg, matrix_digest(&input.matrix))
        }

        let base = changed(|_| {});
        type ConfigChange<'a> = (&'a str, fn(&mut EvalConfig));
        let changes: &[ConfigChange<'_>] = &[
            ("folds", |c| c.folds += 1),
            ("inner_folds", |c| c.inner_folds += 1),
            ("scm_max_rules", |c| c.scm_max_rules += 1),
            ("cart_max_depth", |c| c.cart_max_depth += 1),
            ("forest.n_trees", |c| c.forest.n_trees += 1),
            ("forest.subsample", |c| c.forest.subsample += 0.1),
            ("forest.max_depth", |c| c.forest.max_depth += 1),
            ("forest.min_leaf", |c| c.forest.min_leaf += 1),
            ("forest.seed", |c| c.forest.seed += 1),
            ("gbt.rounds", |c| c.gbt.rounds += 1),
            ("gbt.learning_rate", |c| c.gbt.learning_rate += 0.1),
            ("gbt.max_depth", |c| c.gbt.max_depth += 1),
            ("gbt.lambda", |c| c.gbt.lambda += 1.0),
            ("gbt.min_child_weight", |c| c.gbt.min_child_weight += 1.0),
            ("gbt.colsample", |c| c.gbt.colsample += 0.1),
            ("gbt.patience", |c| c.gbt.patience += 1),
            ("gbt.seed", |c| c.gbt.seed += 1),
            ("screened_features", |c| c.screened_features += 1),
            ("lasso_lambdas", |c| c.lasso_lambdas += 1),
            ("seed", |c| c.seed += 1),
        ];
        for (field, edit) in changes {
            assert_ne!(changed(*edit), base, "{field}");
        }
    }

    #[test]
    fn fingerprint_ignores_rows_and_fields_the_evaluation_does_not_read() {
        let input = inputs();
        let meropenem = dataset("meropenem");
        let Cohort { rows, columns, .. } =
            cohort(&meropenem, &input.phenotypes, &input.matrix.summary.genomes);
        let ids: Vec<&str> = rows.iter().map(|p| p.genome_id.as_str()).collect();
        assert_eq!(
            (ids, columns),
            (vec!["g0", "g1", "g3", "g5"], vec![0, 1, 3, 5])
        );
        let base = key(&input, &meropenem, None);
        let ignored: &[(&str, Edit)] = &[
            ("phenotype file order", |x, _| x.phenotypes.reverse()),
            ("other antibiotic", |x, _| {
                // g1 is also a meropenem row.
                row(x, "g1", "ciprofloxacin").reported = Some(Call::S)
            }),
            ("MIC-inferred row", |x, _| {
                let p = row(x, "g2", "meropenem");
                p.reported = Some(Call::S);
                p.mic = None;
            }),
            ("conflict row", |x, _| x.phenotypes.retain(|p| !p.conflict)),
            ("genome metadata", |x, _| {
                x.genomes.get_mut("g0").unwrap().isolation_country = Some("Elsewhere".into())
            }),
            ("drug class", |_, s| s.drug_class = "carbapenem".into()),
            ("dataset counts", |_, s| {
                s.n_resistant = 30;
                s.n_susceptible = 40;
                s.n_with_mic = 0;
            }),
            ("contigs", |x, _| {
                x.contigs.insert("g0".into(), 900);
            }),
            ("genome map", |x, _| {
                x.structure.coords = vec![[0.5, -0.5]; 6]
            }),
            ("lineage of a genome no row evaluates", |x, _| {
                // Column 2 holds g2, while evaluation row 2 is g3.
                x.structure.lineages.labels[2] = 9
            }),
        ];
        for (what, edit) in ignored {
            assert_eq!(meropenem_key_after(*edit), base, "{what}");
        }
    }

    #[test]
    fn rows_follow_matrix_columns_rather_than_genome_ids() {
        let mut input = inputs();
        reverse_columns(&mut input);
        let Cohort { rows, columns, .. } = cohort(
            &dataset("meropenem"),
            &input.phenotypes,
            &input.matrix.summary.genomes,
        );
        let ids: Vec<&str> = rows.iter().map(|p| p.genome_id.as_str()).collect();
        assert_eq!(
            (ids, columns),
            (vec!["g5", "g3", "g1", "g0"], vec![0, 2, 4, 5])
        );
    }

    #[test]
    fn swapping_genomes_with_identical_columns_changes_the_matrix_digest() {
        let base = inputs();
        let mut swapped = inputs();
        // g1 and g2 carry the same patterns, so the swap moves no bit.
        swapped.matrix.summary.genomes.swap(1, 2);
        assert_eq!(swapped.matrix.patterns.data, base.matrix.patterns.data);
        assert_ne!(matrix_digest(&swapped.matrix), matrix_digest(&base.matrix));
    }

    #[test]
    fn renumbering_singleton_lineages_changes_the_key_and_not_the_matrix_digest() {
        let mut input = inputs();
        let meropenem = dataset("meropenem");
        // g3 and g5 are singleton lineages. Swapping their numbers keeps every
        // grouping, but `grouped_folds` orders lineages by number before its
        // seeded shuffle.
        input.structure.lineages.labels = vec![0, 0, 1, 7, 2, 9];
        let (digest, before) = (matrix_digest(&input.matrix), key(&input, &meropenem, None));
        input.structure.lineages.labels = vec![0, 0, 1, 9, 2, 7];
        assert_eq!(matrix_digest(&input.matrix), digest);
        assert_ne!(key(&input, &meropenem, None), before);
    }

    /// Writes `input` where `load_inputs` reads it, with gene calls recorded
    /// under `calls`: a `Catalog::calls_version`, or "" for none.
    fn write_inputs(layout: &Layout, input: &SpeciesInputs, calls: &str) {
        input.matrix.save(&layout.matrix_dir("fixture")).unwrap();
        let summary = &input.matrix.summary;
        let counted = species::CountSummary {
            k: summary.k,
            genomes: summary
                .genomes
                .iter()
                .zip(&summary.kset_digests)
                .map(|(g, digest)| species::Counted {
                    genome_id: g.clone(),
                    n_kmers: 0,
                    stats: crate::fasta::SeqStats {
                        contigs: input.contigs[g],
                        ..Default::default()
                    },
                    catalog_version: calls.into(),
                    kset_digest: Some(digest.clone()),
                    genes: input.hits.get(g).cloned().unwrap_or_default(),
                })
                .collect(),
            missing_fasta: vec![],
            failed: vec![],
        };
        let dir = layout.species_dir("fixture");
        std::fs::write(
            dir.join("counted.json"),
            serde_json::to_vec(&counted).unwrap(),
        )
        .unwrap();
        std::fs::write(
            dir.join("structure.json"),
            serde_json::to_vec(&input.structure).unwrap(),
        )
        .unwrap();
        std::fs::create_dir_all(layout.curated()).unwrap();
        std::fs::write(
            layout.curated().join("datasets.json"),
            serde_json::to_vec(&input.datasets).unwrap(),
        )
        .unwrap();
        let phenotypes: Vec<String> = input
            .phenotypes
            .iter()
            .map(|p| serde_json::to_string(p).unwrap())
            .collect();
        std::fs::write(
            layout.curated().join("phenotypes.jsonl"),
            phenotypes.join("\n"),
        )
        .unwrap();
        let genomes: Vec<String> = input
            .genomes
            .values()
            .map(|g| serde_json::to_string(g).unwrap())
            .collect();
        crate::layout::ensure_parent(&layout.genomes()).unwrap();
        std::fs::write(layout.genomes(), genomes.join("\n")).unwrap();
    }

    /// Applies `edit` to `file` in the fixture's species directory.
    fn edit_json(layout: &Layout, file: &str, edit: impl FnOnce(&mut serde_json::Value)) {
        let path = layout.species_dir("fixture").join(file);
        let mut value = serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
        edit(&mut value);
        std::fs::write(&path, serde_json::to_vec(&value).unwrap()).unwrap();
    }

    /// Saves a stand-in report under the key `analyze_species` computes for
    /// `antibiotic`, and returns its path and bytes. Neither fixture dataset
    /// has 25 genomes per class, so a recompute returns None and removes the
    /// file: a file that survives unchanged was reused.
    fn save_report(
        layout: &Layout,
        input: &SpeciesInputs,
        catalog: Option<&Catalog>,
        antibiotic: &str,
    ) -> (std::path::PathBuf, Vec<u8>) {
        let results = layout.results_dir("fixture");
        std::fs::create_dir_all(&results).unwrap();
        let path = results.join(format!("fixture__{antibiotic}.json"));
        let saved = serde_json::to_vec(&serde_json::json!({
            "evaluation": {"revision": EVALUATION_REVISION, "valid": true},
            "input_fingerprint": key(input, &dataset(antibiotic), catalog),
        }))
        .unwrap();
        std::fs::write(&path, &saved).unwrap();
        (path, saved)
    }

    fn learn(layout: &Layout) -> Result<Vec<String>> {
        analyze_species(layout, "fixture", &[], false, &EvalConfig::default())
    }

    #[test]
    fn learn_keeps_matching_reports_and_removes_those_of_datasets_that_no_longer_qualify() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let input = inputs();
        write_inputs(&layout, &input, "");
        let (current, saved) = save_report(&layout, &input, None, "meropenem");
        let stale = layout
            .results_dir("fixture")
            .join("fixture__ciprofloxacin.json");
        std::fs::write(
            &stale,
            serde_json::to_vec(&serde_json::json!({
                "evaluation": {"revision": 7},
                "input_fingerprint": key(&input, &dataset("ciprofloxacin"), None),
            }))
            .unwrap(),
        )
        .unwrap();
        // No kset is on disk, as after `--drop-ksets`; the key reads the
        // digests matrix.json records.
        assert!(!layout.kset_dir("fixture").exists());
        assert!(learn(&layout).unwrap().is_empty());
        assert_eq!(std::fs::read(&current).unwrap(), saved);
        assert!(!stale.exists());

        // One flipped pattern bit is an input change, with matrix.json as it was.
        let patterns = layout.matrix_dir("fixture").join("patterns.bin");
        let mut bytes = std::fs::read(&patterns).unwrap();
        bytes[0] ^= 1;
        std::fs::write(&patterns, bytes).unwrap();
        learn(&layout).unwrap();
        assert!(!current.exists());
    }

    #[test]
    fn learn_recomputes_a_current_report_whose_audit_is_invalid() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let input = inputs();
        write_inputs(&layout, &input, "");
        let (report, saved) = save_report(&layout, &input, None, "meropenem");
        let mut invalid: serde_json::Value = serde_json::from_slice(&saved).unwrap();
        invalid["evaluation"]["valid"] = false.into();
        invalid["evaluation"]["reasons"] =
            serde_json::json!(["random split outer fold 1 has an unconverged two-class lasso fit"]);
        std::fs::write(&report, serde_json::to_vec(&invalid).unwrap()).unwrap();
        // Same revision and fingerprint, so only the audit decides: the
        // report is recomputed, which for this fixture removes it.
        learn(&layout).unwrap();
        assert!(!report.exists());
    }

    /// 150 genomes in 30 lineages of five. Meropenem calls every third
    /// genome resistant; two features mark resistance and three are noise.
    fn scan_inputs() -> SpeciesInputs {
        let n = 150;
        let genomes: Vec<String> = (0..n).map(|i| format!("g{i}")).collect();
        let kset_digests: Vec<String> = genomes.iter().map(|g| format!("kset of {g}")).collect();
        let rules: [fn(usize) -> bool; 5] = [
            |i| (i % 3 == 0 && i % 10 != 0) || i % 20 == 1,
            |i| i % 3 == 0 && i % 7 != 0,
            |i| i % 5 == 0,
            |i| i % 7 == 1,
            |i| i % 4 < 2,
        ];
        let mut patterns = BitMatrix::with_cols(n);
        let mut support = Vec::new();
        for rule in rules {
            let row = Bits::from_indices(n, (0..n).filter(|&i| rule(i)));
            support.push(row.count());
            patterns.push_row(row.words());
        }
        SpeciesInputs {
            matrix: FeatureMatrix {
                summary: MatrixSummary {
                    k: 7,
                    genomes: genomes.clone(),
                    min_support: 1,
                    n_patterns: rules.len(),
                    kmers_distinct: 5,
                    kmers_kept: 5,
                    kmers_rare: 0,
                    kmers_core: 0,
                    kset_digests: kset_digests.clone(),
                },
                patterns,
                support,
                kmer_offsets: (0..=5).collect(),
                kmers: (0x21..0x26).collect(),
            },
            structure: Structure {
                genomes: genomes.clone(),
                k: 7,
                kset_digests,
                lineages: crate::cluster::Lineages {
                    threshold: 0.005,
                    labels: (0..n).map(|i| (i / 5) as u32).collect(),
                    n_clusters: 30,
                    largest_share: 1.0 / 30.0,
                },
                merges: vec![],
                leaf_order: vec![],
                coords: vec![],
                nearest: vec![],
                sketch_size: species::SKETCH,
            },
            hits: HashMap::new(),
            contigs: genomes.iter().map(|g| (g.clone(), 40)).collect(),
            phenotypes: genomes
                .iter()
                .enumerate()
                .map(|(i, g)| {
                    let call = if i % 3 == 0 { Call::R } else { Call::S };
                    phenotype(g, "meropenem", call, None)
                })
                .collect(),
            genomes: genomes
                .iter()
                .map(|g| {
                    let record = GenomeRecord {
                        genome_id: g.clone(),
                        ..Default::default()
                    };
                    (g.clone(), record)
                })
                .collect(),
            datasets: vec![dataset("meropenem")],
        }
    }

    /// Every file under `dir` with its size and modification time.
    fn tree(dir: &std::path::Path) -> Vec<(std::path::PathBuf, u64, std::time::SystemTime)> {
        let mut out = Vec::new();
        for entry in std::fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            let meta = std::fs::metadata(&path).unwrap();
            if meta.is_dir() {
                out.extend(tree(&path));
            } else {
                out.push((path, meta.len(), meta.modified().unwrap()));
            }
        }
        out.sort();
        out
    }

    #[test]
    fn the_lasso_scan_fits_the_report_inputs_as_the_evaluator_does() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let input = scan_inputs();
        write_inputs(&layout, &input, "");
        let before = tree(dir.path());
        let cfg = EvalConfig {
            forest: ForestParams {
                n_trees: 4,
                ..Default::default()
            },
            gbt: GbtParams {
                rounds: 5,
                ..Default::default()
            },
            lasso_lambdas: 10,
            ..Default::default()
        };
        let spec = dataset("meropenem");
        let ids = [spec.id.clone()];
        let mut events = Vec::new();
        let mut points = 0;
        lasso_scan(
            &layout,
            "fixture",
            &ids,
            &cfg,
            true,
            &mut |id, e| {
                assert_eq!(id, spec.id);
                events.push(e.clone());
                Ok(())
            },
            &mut |id, _| {
                assert_eq!(id, spec.id);
                points += 1;
            },
        )
        .unwrap();
        // The scan reads the data root and writes nothing to it.
        assert_eq!(tree(dir.path()), before);

        let LassoOnlyEvent::Inputs {
            samples,
            resistant,
            lineages,
            min_support,
            input_fingerprint: fingerprint,
            ..
        } = &events[0]
        else {
            panic!("{:?}", events[0]);
        };
        assert_eq!(
            (*samples, *resistant, *lineages, *min_support),
            (150, 50, 30, 3)
        );
        let digest = matrix_digest(&input.matrix);
        assert_eq!(
            *fingerprint,
            input_fingerprint(&input, &spec, None, &cfg, digest)
        );

        // Its fits are the evaluator's lasso fits on the report's inputs.
        let e = eval_inputs(&input, &spec).unwrap().unwrap();
        let out = eval::evaluate(
            &spec.id,
            &e.ds.x,
            &e.ds.y,
            &e.groups,
            &e.duplicates,
            e.min_support,
            None,
            &cfg,
        )
        .unwrap();
        let fits: Vec<&eval::LassoOnlyFit> = events
            .iter()
            .filter_map(|e| match e {
                LassoOnlyEvent::Fit(f) => Some(f),
                _ => None,
            })
            .collect();
        assert_eq!(fits.len(), 2 * cfg.folds + 1);
        let traced: usize = fits
            .iter()
            .map(|f| f.trace.path.len() + f.trace.refit.len())
            .sum();
        assert_eq!(points, traced);
        for split in &out.splits {
            assert!(split.unavailable_reason.is_none());
            let lasso = split
                .learners
                .iter()
                .find(|l| l.learner == "lasso")
                .unwrap();
            for (fold, audit) in lasso.lasso_fits.as_ref().unwrap().iter().enumerate() {
                let fit = fits
                    .iter()
                    .find(|f| f.scheme == split.scheme && f.outer_fold == Some(fold))
                    .unwrap();
                assert_eq!(fit.audit.as_ref(), Some(audit));
            }
        }
        let last = fits.last().unwrap();
        assert_eq!(last.audit, out.finals.lasso_fit);
        assert_eq!(
            last.model.as_ref().map(|m| &m.coef),
            out.finals.lasso.as_ref().map(|m| &m.coef)
        );

        // Without fits it reports the inputs alone, and it refuses a dataset
        // the curation does not list before fitting anything.
        let mut inputs_only = Vec::new();
        lasso_scan(
            &layout,
            "fixture",
            &ids,
            &cfg,
            false,
            &mut |_, e| {
                inputs_only.push(e.clone());
                Ok(())
            },
            &mut |_, _| unreachable!("nothing is fitted"),
        )
        .unwrap();
        assert!(matches!(&inputs_only[..], [LassoOnlyEvent::Inputs { .. }]));
        let missing = ["fixture__nothing".to_string()];
        let error = lasso_scan(
            &layout,
            "fixture",
            &missing,
            &cfg,
            true,
            &mut |_, _| unreachable!("nothing to report"),
            &mut |_, _| unreachable!("nothing is fitted"),
        )
        .unwrap_err();
        assert!(format!("{error:#}").contains("no curated dataset fixture__nothing"));
    }

    /// Writes `scan_inputs` as a data root, with a stand-in revision-7 report
    /// whose fingerprint matches, for the scan harness's end-to-end test
    /// (`scripts/lasso-scan/smoke.sh`).
    #[test]
    #[ignore = "writes the scan smoke root; set GRM_LASSO_SCAN_FIXTURE to a new directory"]
    fn write_lasso_scan_fixture_root() {
        let root = std::path::PathBuf::from(
            std::env::var("GRM_LASSO_SCAN_FIXTURE").expect("set GRM_LASSO_SCAN_FIXTURE"),
        );
        assert!(
            std::fs::read_dir(&root).map_or(true, |mut d| d.next().is_none()),
            "{} must be new or empty",
            root.display()
        );
        let layout = Layout::new(&root);
        let input = scan_inputs();
        write_inputs(&layout, &input, "");
        let spec = dataset("meropenem");
        let results = layout.results_dir("fixture");
        std::fs::create_dir_all(&results).unwrap();
        // The scan reads revision-7 lasso seconds from each split's lasso row:
        // far more than a debug build takes on this fixture, so a clean smoke
        // run is a go.
        let split = |scheme: &str| {
            serde_json::json!({
                "scheme": scheme,
                "folds": 5,
                "learners": [{"learner": "lasso", "seconds": 30.0}],
            })
        };
        let report = serde_json::json!({
            "evaluation": {"revision": 7, "valid": true},
            "id": spec.id,
            "species_slug": "fixture",
            "input_fingerprint": key(&input, &spec, None),
            "seconds": 3.0,
            "splits": [split("random"), split("lineage")],
        });
        std::fs::write(
            results.join(format!("{}.json", spec.id)),
            serde_json::to_vec(&report).unwrap(),
        )
        .unwrap();
    }

    #[test]
    fn a_rebuild_in_reversed_genome_order_recomputes_and_a_new_map_alone_does_not() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let mut input = inputs();
        input.structure.lineages.labels = vec![0, 0, 1, 7, 2, 9];
        input.structure.coords = (0..6).map(|i| [i as f32, 0.0]).collect();
        write_inputs(&layout, &input, "");
        let (report, saved) = save_report(&layout, &input, None, "meropenem");

        input.structure.coords = (0..6).map(|i| [0.0, i as f32]).collect();
        write_inputs(&layout, &input, "");
        learn(&layout).unwrap();
        assert_eq!(std::fs::read(&report).unwrap(), saved);

        // Structure and matrix rebuilt after a recount read the genome list
        // backwards: reversed columns, a new map, and the singleton lineages
        // of g3 and g5 numbered the other way round.
        reverse_columns(&mut input);
        input.structure.lineages.labels = vec![7, 2, 9, 1, 0, 0];
        input.structure.coords = (0..6).map(|i| [i as f32, i as f32]).collect();
        write_inputs(&layout, &input, "");
        learn(&layout).unwrap();
        assert!(!report.exists());
    }

    /// Why `grm learn` refuses the fixture after `edit`, which leaves the
    /// saved report as it was.
    fn refusal_after(edit: impl FnOnce(&Layout)) -> String {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let input = inputs();
        write_inputs(&layout, &input, "");
        let (report, saved) = save_report(&layout, &input, None, "meropenem");
        edit(&layout);
        let err = format!("{:#}", learn(&layout).unwrap_err());
        assert_eq!(std::fs::read(&report).unwrap(), saved);
        err
    }

    #[test]
    fn learn_refuses_counts_recorded_before_kset_digests() {
        let err = refusal_after(|layout| {
            edit_json(layout, "counted.json", |c| {
                for g in c["genomes"].as_array_mut().unwrap() {
                    g.as_object_mut().unwrap().remove("kset_digest");
                }
            })
        });
        assert_eq!(
            err,
            "counted.json records no kset digest for g0 (it predates them); rerun `grm count --species fixture`"
        );
    }

    #[test]
    fn learn_refuses_a_structure_or_matrix_not_built_from_the_counts() {
        // `grm count --k 21` after both were built at k = 7.
        let err = refusal_after(|l| edit_json(l, "counted.json", |c| c["k"] = 21.into()));
        assert_eq!(
            err,
            "structure.json was built at k = 7 and counted.json records k = 21; matrix.json was built at k = 7 and counted.json records k = 21; rerun `grm structure` and `grm matrix`"
        );
        let err = refusal_after(|l| {
            edit_json(l, "counted.json", |c| {
                c["genomes"].as_array_mut().unwrap().reverse()
            })
        });
        let other_order =
            "lists other genomes than counted.json, or the same genomes in another order";
        assert_eq!(
            err,
            format!(
                "structure.json {other_order}; matrix.json {other_order}; rerun `grm structure` and `grm matrix`"
            )
        );
        // Each file as `grm structure` or `grm matrix` wrote it before kset
        // digests, beside the other one current.
        let err = refusal_after(|l| {
            edit_json(l, "structure.json", |s| {
                let s = s.as_object_mut().unwrap();
                s.remove("k");
                s.remove("kset_digests");
            })
        });
        assert_eq!(
            err,
            "structure.json records no kset digests (it predates them); rerun `grm structure`"
        );
        let err = refusal_after(|l| {
            edit_json(l, "matrix/matrix.json", |m| {
                m.as_object_mut().unwrap().remove("kset_digests");
            })
        });
        assert_eq!(
            err,
            "matrix.json records no kset digests (it predates them); rerun `grm matrix`"
        );
        let err = refusal_after(|l| {
            edit_json(l, "matrix/matrix.json", |m| {
                m["kset_digests"][1] = "another kset of g1".into()
            })
        });
        assert_eq!(
            err,
            "matrix.json was built from other ksets than counted.json records; rerun `grm matrix`"
        );
    }

    /// Writes a catalog with one synthetic gene at `version` and loads it at the
    /// fixture's k.
    fn write_catalog(layout: &Layout, version: &str) -> Catalog {
        let dir = layout.catalog_dir();
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(
            dir.join("ReferenceGeneCatalog.txt"),
            "type\tclass\tsubclass\trefseq_protein_accession\tgenbank_protein_accession\tgene_family\tproduct_name\n\
             AMR\tfixture\tfixture\tREF\t\tfixture\tfixture\n",
        ).unwrap();
        std::fs::write(
            dir.join("AMR_CDS.fa"),
            format!(">REF\n{}\n", "ACGT".repeat(16)),
        )
        .unwrap();
        std::fs::write(dir.join("version.txt"), version).unwrap();
        Catalog::load(&dir, &KmerSpace::new(7)).unwrap()
    }

    #[test]
    fn learn_refuses_gene_calls_from_another_catalog_or_one_that_no_longer_loads() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let input = inputs();
        let v1 = write_catalog(&layout, "v1");
        write_inputs(&layout, &input, &v1.calls_version());
        let (report, saved) = save_report(&layout, &input, Some(&v1), "meropenem");
        learn(&layout).unwrap();
        assert_eq!(std::fs::read(&report).unwrap(), saved);

        write_catalog(&layout, "v2");
        let err = format!("{:#}", learn(&layout).unwrap_err());
        assert!(
            err.contains("g0 has gene calls from catalog v1|group-calls-v2, and"),
            "{err}"
        );
        assert!(
            err.contains("now holds catalog v2|group-calls-v2; rerun `grm count`"),
            "{err}"
        );
        assert_eq!(std::fs::read(&report).unwrap(), saved);

        std::fs::remove_dir_all(layout.catalog_dir()).unwrap();
        let err = format!("{:#}", learn(&layout).unwrap_err());
        assert!(
            err.contains("which no longer loads: run `grm fetch catalog` first"),
            "{err}"
        );
        assert_eq!(std::fs::read(&report).unwrap(), saved);
    }

    #[test]
    fn reports_use_the_reported_cohort_and_record_its_counts_and_k() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        // The default minimum admits 25 reported resistant and 25 reported
        // susceptible genomes. One MIC-only matrix row stays out.
        let genomes: Vec<String> = (0..51).map(|i| format!("g{i}")).collect();
        for g in &genomes {
            // Annotation reads a reference assembly and its cached features.
            let fasta = crate::bvbrc::fasta_path(&layout.fasta_dir(), g);
            crate::layout::ensure_parent(&fasta).unwrap();
            std::fs::write(fasta, format!(">c\n{WILD}\n")).unwrap();
            let features = layout.features_dir().join(format!("{g}.json"));
            crate::layout::ensure_parent(&features).unwrap();
            std::fs::write(features, "[]").unwrap();
        }
        let mut phenotypes: Vec<Phenotype> = genomes[..50]
            .iter()
            .enumerate()
            .map(|(i, g)| {
                let call = if i < 25 { Call::R } else { Call::S };
                phenotype(g, "meropenem", call, None)
            })
            .collect();
        phenotypes.push(inferred("g50", "meropenem", Call::R));
        let input = SpeciesInputs {
            matrix: FeatureMatrix {
                summary: MatrixSummary {
                    k: 21,
                    genomes: genomes.clone(),
                    min_support: 3,
                    n_patterns: 2,
                    kmers_distinct: 2,
                    kmers_kept: 2,
                    kmers_rare: 0,
                    kmers_core: 0,
                    kset_digests: genomes.iter().map(|g| format!("kset of {g}")).collect(),
                },
                patterns: BitMatrix {
                    rows: 2,
                    cols: 51,
                    words: 1,
                    // The resistant genomes, and every other reported genome.
                    data: vec![(1 << 25) - 1, 0x5555_5555_5555_5555 >> 14],
                },
                support: vec![25, 25],
                kmer_offsets: vec![0, 1, 2],
                kmers: vec![0x21, 0x22],
            },
            structure: Structure {
                genomes: genomes.clone(),
                k: 21,
                kset_digests: genomes.iter().map(|g| format!("kset of {g}")).collect(),
                lineages: crate::cluster::Lineages {
                    threshold: 0.005,
                    labels: (0..51).map(|i| i % 10).collect(),
                    n_clusters: 10,
                    largest_share: 11.0 / 51.0,
                },
                merges: vec![],
                leaf_order: vec![],
                coords: vec![],
                nearest: vec![],
                sketch_size: species::SKETCH,
            },
            hits: HashMap::new(),
            contigs: HashMap::new(),
            phenotypes,
            genomes: HashMap::new(),
            datasets: vec![],
        };
        let spec = DatasetSpec {
            breakpoint: None,
            ..dataset("meropenem")
        };
        let report = analyze_dataset(
            &layout,
            &Client::new(),
            &input,
            &spec,
            None,
            &EvalConfig::default(),
            &mut HashMap::new(),
            "fixture key".into(),
        )
        .unwrap()
        .unwrap();
        assert_eq!(report.counts.k, 21);
        assert_eq!(report.counts.min_class, 25);
        assert_eq!(report.counts.samples, 50);
        assert_eq!(report.counts.resistant, 25);
        assert_eq!(report.counts.susceptible, 25);
        assert_eq!(report.counts.labels_from_mic, 0);
        assert_eq!(report.counts.labels_from_mic_excluded, 1);
        assert!(report.samples.iter().all(|s| s.id != "g50" && !s.from_mic));
        assert_eq!(report.input_fingerprint, "fixture key");
    }
}
