//! Curation: raw BV-BRC records to one phenotype per genome and drug.
//!
//! Every exclusion is counted, so the report can show the path from raw
//! records to the rows a model trained on.
//!
//! Labels come from two places. A reported call (Resistant, Susceptible,
//! Intermediate) wins when present. Records that only carry an MIC get a
//! label from an empirical breakpoint: for each species and drug, the MIC
//! threshold that best reproduces the reported calls on records that have
//! both. The threshold is used only when it reproduces at least 95% of them.
//!
//! Only unconflicted reported R/S calls count toward a dataset's class
//! minimum and the download ranking. MIC-inferred labels stay in the curated
//! phenotypes, marked `label_from_mic`, and are counted per dataset.

use crate::bvbrc::{AmrRecord, GenomeRecord};
use crate::drugs::{self, NameKind, Resolver};
use anyhow::{Result, ensure};
use rand::SeedableRng;
use rand::seq::SliceRandom;
use rand_chacha::ChaCha8Rng;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, HashMap};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum Call {
    S,
    I,
    R,
}

pub fn parse_call(text: &str) -> Option<Call> {
    match text.trim().to_lowercase().as_str() {
        "resistant" | "nonsusceptible" | "non-susceptible" => Some(Call::R),
        "susceptible" => Some(Call::S),
        "intermediate" | "susceptible-dose dependent" | "reduced susceptibility" => Some(Call::I),
        _ => None,
    }
}

/// How an MIC is bounded. Values are log2 dilution steps (0 = 1 mg/L).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Censor {
    Exact,
    /// MIC <= value (the lowest concentration tested already inhibited).
    AtMost,
    /// MIC > value (the highest concentration tested did not inhibit).
    Above,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Mic {
    pub log2: i32,
    pub censor: Censor,
}

/// Parses an MIC in mg/L. Zone diameters (disk diffusion, mm) are rejected.
pub fn parse_mic(r: &AmrRecord) -> Option<Mic> {
    let unit = r.measurement_unit.as_deref().map(str::to_lowercase);
    let method = r
        .laboratory_typing_method
        .as_deref()
        .unwrap_or("")
        .to_lowercase();
    let disk = method.contains("disk") || method.contains("disc");
    // Some disk-diffusion submissions label zone diameters "mg/L". Under that
    // method, only a gradient strip (E-test) reports an MIC.
    let gradient_strip = r
        .laboratory_typing_platform
        .as_deref()
        .is_some_and(|p| p.eq_ignore_ascii_case("E-test"));
    let mic_unit = match unit.as_deref() {
        Some("mg/l") => !disk || gradient_strip,
        Some(_) => false,
        None => !method.is_empty() && !disk,
    };
    if !mic_unit {
        return None;
    }
    let raw = r
        .measurement_value
        .as_deref()
        .or(r.measurement.as_deref())?;
    let (sign_in_value, number) = split_sign(raw);
    let sign = r
        .measurement_sign
        .as_deref()
        .filter(|s| !s.is_empty())
        .unwrap_or(sign_in_value);
    // Combination drugs report "2/38": the first component carries the MIC.
    let first = number.split('/').next()?.trim();
    let value: f64 = first.parse().ok()?;
    if !(value > 0.0 && value.is_finite()) {
        return None;
    }
    let step = value.log2().round() as i32;
    Some(match sign {
        "<=" | "=<" => Mic {
            log2: step,
            censor: Censor::AtMost,
        },
        "<" => Mic {
            log2: step - 1,
            censor: Censor::AtMost,
        },
        ">" => Mic {
            log2: step,
            censor: Censor::Above,
        },
        ">=" | "=>" => Mic {
            log2: step - 1,
            censor: Censor::Above,
        },
        _ => Mic {
            log2: step,
            censor: Censor::Exact,
        },
    })
}

fn split_sign(raw: &str) -> (&str, &str) {
    let raw = raw.trim();
    for sign in ["<=", ">=", "=<", "=>", "<", ">", "="] {
        if let Some(rest) = raw.strip_prefix(sign) {
            return (sign, rest.trim());
        }
    }
    ("", raw)
}

/// MIC above `r_above` (log2) is resistant.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Breakpoint {
    pub r_above: i32,
    /// Share of records with both a call and an MIC that the threshold
    /// classifies the same way (ambiguous censored values count against it).
    pub agreement: f64,
    pub n: usize,
}

/// Predicted call for an MIC under a threshold; None when censoring leaves
/// the side of the threshold unknown.
pub fn call_from_mic(mic: Mic, r_above: i32) -> Option<Call> {
    match mic.censor {
        Censor::Exact => Some(if mic.log2 > r_above { Call::R } else { Call::S }),
        Censor::AtMost => (mic.log2 <= r_above).then_some(Call::S),
        Censor::Above => (mic.log2 >= r_above).then_some(Call::R),
    }
}

pub fn infer_breakpoint(pairs: &[(Mic, Call)]) -> Option<Breakpoint> {
    let n = pairs.len();
    let n_r = pairs.iter().filter(|(_, c)| *c == Call::R).count();
    if n < 50 || n_r < 10 || n - n_r < 10 {
        return None;
    }
    let lo = pairs.iter().map(|(m, _)| m.log2).min()?;
    let hi = pairs.iter().map(|(m, _)| m.log2).max()?;
    let mut best: Option<Breakpoint> = None;
    for t in lo - 1..=hi {
        let correct = pairs
            .iter()
            .filter(|(m, c)| call_from_mic(*m, t) == Some(*c))
            .count();
        let agreement = correct as f64 / n as f64;
        if best.is_none_or(|b| agreement > b.agreement) {
            best = Some(Breakpoint {
                r_above: t,
                agreement,
                n,
            });
        }
    }
    best.filter(|b| b.agreement >= 0.95)
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct QcConfig {
    pub max_contigs: u32,
    pub min_completeness: f64,
    pub max_contamination: f64,
    /// Allowed relative deviation from the species median genome length.
    pub length_tolerance: f64,
    /// Minimum genomes with a reported call in each class for a dataset to
    /// exist.
    pub min_class: usize,
}

impl Default for QcConfig {
    fn default() -> Self {
        QcConfig {
            max_contigs: 1000,
            min_completeness: 90.0,
            max_contamination: 10.0,
            length_tolerance: 0.25,
            min_class: 25,
        }
    }
}

/// One genome x antibiotic after consolidation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Phenotype {
    pub genome_id: String,
    pub species: String,
    pub antibiotic: String,
    /// Consolidated reported call; None when absent or contradictory.
    pub reported: Option<Call>,
    /// The records reported both Resistant and Susceptible.
    pub conflict: bool,
    pub mic: Option<Mic>,
    /// Spread in dilution steps among exact MICs for this pair.
    pub mic_spread: u32,
    /// Final label: the reported call, else the MIC under the breakpoint.
    pub label: Option<Call>,
    pub label_from_mic: bool,
    pub methods: Vec<String>,
    pub standards: Vec<String>,
    pub n_records: u32,
}

impl Phenotype {
    /// The unconflicted reported R or S call, the only label that counts
    /// toward a dataset, ranks genomes for download, or trains and scores a
    /// model. MIC-inferred labels come from a threshold fitted on the whole
    /// cohort, held-out outcomes included, so they are never targets.
    pub fn reported_label(&self) -> Option<Call> {
        match self.reported {
            Some(call @ (Call::R | Call::S)) if !self.conflict && !self.label_from_mic => {
                Some(call)
            }
            _ => None,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DatasetSpec {
    pub id: String,
    pub species: String,
    pub species_slug: String,
    pub antibiotic: String,
    pub drug_class: String,
    /// Genomes with an unconflicted reported call in each class. Curated
    /// files without `min_class` counted MIC-inferred labels here too.
    pub n_resistant: usize,
    pub n_susceptible: usize,
    pub n_intermediate: usize,
    pub n_with_mic: usize,
    /// MIC-inferred labels for this species and drug, kept as provenance and
    /// counted in no class.
    pub n_labels_from_mic: usize,
    pub breakpoint: Option<Breakpoint>,
    /// The class minimum this dataset was curated with, which the analysis
    /// gate applies again after download.
    #[serde(default = "unrecorded_min_class")]
    pub min_class: usize,
}

/// Curated files and dataset reports written before `min_class` was recorded
/// do not say which `--min-class` curation used, but `grm learn` gated every
/// dataset at 25 then, so a missing value reads as that gate.
pub fn unrecorded_min_class() -> usize {
    25
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Funnel {
    pub records: usize,
    pub records_class_label: usize,
    pub records_not_antibiotic: usize,
    pub records_unresolved_name: usize,
    pub records_no_measure: usize,
    pub genomes_with_records: usize,
    pub genomes_without_metadata: usize,
    pub genomes_failed_qc: BTreeMap<String, usize>,
    pub genomes_passing_qc: usize,
    pub pairs: usize,
    pub pairs_conflict: usize,
    pub pairs_intermediate_only: usize,
    pub pairs_reported: usize,
    pub pairs_labeled_from_mic: usize,
    pub pairs_unlabeled: usize,
}

pub struct Curated {
    pub phenotypes: Vec<Phenotype>,
    pub datasets: Vec<DatasetSpec>,
    pub funnel: Funnel,
    pub qc_failures: HashMap<String, String>,
}

pub fn slug(species: &str) -> String {
    species
        .to_lowercase()
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect::<String>()
        .split('_')
        .filter(|s| !s.is_empty())
        .collect::<Vec<_>>()
        .join("_")
}

fn num(v: &Option<String>) -> Option<f64> {
    v.as_deref().and_then(|s| s.parse().ok())
}

/// Two-word species from metadata, falling back to the genome name.
fn species_of(g: &GenomeRecord) -> Option<String> {
    let text = g.species.as_deref().or(g.genome_name.as_deref())?;
    let words: Vec<&str> = text.split_whitespace().take(2).collect();
    (words.len() == 2 && words[0].chars().next()?.is_uppercase()).then(|| words.join(" "))
}

/// The QC reason for genomes `grm species-check` failed.
pub const SPECIES_CHECK_FAILED: &str = "species check: sequence is another species";

fn qc_failure(g: &GenomeRecord, median_len: f64, cfg: &QcConfig) -> Option<&'static str> {
    if g.genome_quality.as_deref() == Some("Poor") {
        return Some("BV-BRC quality: Poor");
    }
    if num(&g.contigs).is_some_and(|c| c > cfg.max_contigs as f64) {
        return Some("too many contigs");
    }
    if num(&g.checkm_completeness).is_some_and(|c| c < cfg.min_completeness) {
        return Some("CheckM completeness low");
    }
    if num(&g.checkm_contamination).is_some_and(|c| c > cfg.max_contamination) {
        return Some("CheckM contamination high");
    }
    match num(&g.genome_length) {
        None => Some("no genome length"),
        Some(len) if median_len > 0.0 && (len / median_len - 1.0).abs() > cfg.length_tolerance => {
            Some("genome length off species median")
        }
        _ => None,
    }
}

/// `wrong_species` maps the genomes `grm species-check` failed to why. They
/// are excluded before any other QC, so no dataset counts them.
pub fn curate(
    records: &[AmrRecord],
    genomes: &[GenomeRecord],
    cfg: &QcConfig,
    wrong_species: &BTreeMap<String, String>,
    resolver: &mut Resolver,
) -> Result<Curated> {
    let mut funnel = Funnel {
        records: records.len(),
        ..Default::default()
    };
    let meta: HashMap<&str, &GenomeRecord> =
        genomes.iter().map(|g| (g.genome_id.as_str(), g)).collect();

    // Species medians over genomes that have metadata and a species.
    let mut lengths: HashMap<String, Vec<f64>> = HashMap::new();
    for g in genomes {
        if let (Some(sp), Some(len)) = (species_of(g), num(&g.genome_length)) {
            lengths.entry(sp).or_default().push(len);
        }
    }
    let medians: HashMap<String, f64> = lengths
        .into_iter()
        .map(|(sp, mut v)| {
            v.sort_by(f64::total_cmp);
            (sp, v[v.len() / 2])
        })
        .collect();

    let mut genome_ids: BTreeSet<&str> = BTreeSet::new();
    let mut passing: HashMap<&str, String> = HashMap::new();
    let mut qc_failures = HashMap::new();
    for r in records {
        genome_ids.insert(&r.genome_id);
    }
    funnel.genomes_with_records = genome_ids.len();
    for id in &genome_ids {
        let Some(g) = meta.get(id) else {
            funnel.genomes_without_metadata += 1;
            continue;
        };
        if let Some(why) = wrong_species.get(*id) {
            *funnel
                .genomes_failed_qc
                .entry(SPECIES_CHECK_FAILED.into())
                .or_default() += 1;
            qc_failures.insert(id.to_string(), format!("{SPECIES_CHECK_FAILED}: {why}"));
            continue;
        }
        let Some(sp) = species_of(g) else {
            *funnel
                .genomes_failed_qc
                .entry("no species".into())
                .or_default() += 1;
            qc_failures.insert(id.to_string(), "no species".to_string());
            continue;
        };
        match qc_failure(g, medians.get(&sp).copied().unwrap_or(0.0), cfg) {
            Some(reason) => {
                *funnel.genomes_failed_qc.entry(reason.into()).or_default() += 1;
                qc_failures.insert(id.to_string(), reason.to_string());
            }
            None => {
                passing.insert(id, sp);
            }
        }
    }
    funnel.genomes_passing_qc = passing.len();

    #[derive(Default)]
    struct Acc {
        calls: BTreeSet<Call>,
        mics: Vec<Mic>,
        methods: BTreeSet<String>,
        standards: BTreeSet<String>,
        n: u32,
    }
    let mut pairs: BTreeMap<(String, String), Acc> = BTreeMap::new();
    for r in records {
        let name = resolver.resolve(&r.antibiotic);
        match name.kind {
            NameKind::Drug => {}
            NameKind::Class => {
                funnel.records_class_label += 1;
                continue;
            }
            NameKind::Excluded => {
                funnel.records_not_antibiotic += 1;
                continue;
            }
            NameKind::Unresolved => {
                funnel.records_unresolved_name += 1;
                continue;
            }
        }
        let antibiotic = name.canonical.clone();
        if !passing.contains_key(r.genome_id.as_str()) {
            continue;
        }
        let call = r.resistant_phenotype.as_deref().and_then(parse_call);
        let mic = parse_mic(r);
        if call.is_none() && mic.is_none() {
            funnel.records_no_measure += 1;
            continue;
        }
        let acc = pairs.entry((r.genome_id.clone(), antibiotic)).or_default();
        acc.n += 1;
        acc.calls.extend(call);
        acc.mics.extend(mic);
        acc.methods.extend(r.laboratory_typing_method.clone());
        acc.standards
            .extend(r.testing_standard.clone().map(|s| s.to_uppercase()));
    }

    let mut phenotypes: Vec<Phenotype> = pairs
        .into_iter()
        .map(|((genome_id, antibiotic), acc)| {
            let conflict = acc.calls.contains(&Call::R) && acc.calls.contains(&Call::S);
            let reported = if conflict {
                None
            } else if acc.calls.contains(&Call::R) {
                Some(Call::R)
            } else if acc.calls.contains(&Call::S) {
                Some(Call::S)
            } else {
                acc.calls.iter().next().copied()
            };
            let (mic, mic_spread) = consolidate_mics(&acc.mics);
            Phenotype {
                species: passing[genome_id.as_str()].clone(),
                genome_id,
                antibiotic,
                reported,
                conflict,
                mic,
                mic_spread,
                label: reported,
                label_from_mic: false,
                methods: acc.methods.into_iter().collect(),
                standards: acc.standards.into_iter().collect(),
                n_records: acc.n,
            }
        })
        .collect();

    // Empirical breakpoints per species and drug.
    let mut both: HashMap<(String, String), Vec<(Mic, Call)>> = HashMap::new();
    for p in &phenotypes {
        if let (Some(mic), Some(call @ (Call::R | Call::S))) = (p.mic, p.reported) {
            both.entry((p.species.clone(), p.antibiotic.clone()))
                .or_default()
                .push((mic, call));
        }
    }
    let breakpoints: HashMap<(String, String), Breakpoint> = both
        .into_iter()
        .filter_map(|(key, v)| infer_breakpoint(&v).map(|b| (key, b)))
        .collect();

    for p in &mut phenotypes {
        funnel.pairs += 1;
        if p.conflict {
            funnel.pairs_conflict += 1;
            continue;
        }
        match p.reported {
            Some(Call::R | Call::S) => funnel.pairs_reported += 1,
            Some(Call::I) => funnel.pairs_intermediate_only += 1,
            None => {
                let bp = breakpoints.get(&(p.species.clone(), p.antibiotic.clone()));
                match (p.mic, bp) {
                    (Some(mic), Some(bp)) => match call_from_mic(mic, bp.r_above) {
                        Some(call) => {
                            p.label = Some(call);
                            p.label_from_mic = true;
                            funnel.pairs_labeled_from_mic += 1;
                        }
                        None => funnel.pairs_unlabeled += 1,
                    },
                    _ => funnel.pairs_unlabeled += 1,
                }
            }
        }
    }

    let mut grouped: BTreeMap<(String, String), Vec<&Phenotype>> = BTreeMap::new();
    for p in &phenotypes {
        grouped
            .entry((p.species.clone(), p.antibiotic.clone()))
            .or_default()
            .push(p);
    }
    let datasets = grouped
        .into_iter()
        .filter_map(|((species, antibiotic), rows)| {
            let reported = |c: Call| {
                rows.iter()
                    .filter(|p| p.reported_label() == Some(c))
                    .count()
            };
            let (n_r, n_s) = (reported(Call::R), reported(Call::S));
            if n_r < cfg.min_class || n_s < cfg.min_class {
                return None;
            }
            let species_slug = slug(&species);
            Some(DatasetSpec {
                id: format!("{species_slug}__{}", slug(&antibiotic)),
                drug_class: drugs::class_of(&antibiotic),
                n_resistant: n_r,
                n_susceptible: n_s,
                n_intermediate: rows.iter().filter(|p| p.label == Some(Call::I)).count(),
                n_with_mic: rows.iter().filter(|p| p.mic.is_some()).count(),
                n_labels_from_mic: rows.iter().filter(|p| p.label_from_mic).count(),
                breakpoint: breakpoints
                    .get(&(species.clone(), antibiotic.clone()))
                    .copied(),
                min_class: cfg.min_class,
                species,
                species_slug,
                antibiotic,
            })
        })
        .collect();

    Ok(Curated {
        phenotypes,
        datasets,
        funnel,
        qc_failures,
    })
}

/// Median of exact MICs when any exist; otherwise the tightest bound when
/// all censored values agree in direction.
fn consolidate_mics(mics: &[Mic]) -> (Option<Mic>, u32) {
    let mut exact: Vec<i32> = mics
        .iter()
        .filter(|m| m.censor == Censor::Exact)
        .map(|m| m.log2)
        .collect();
    if !exact.is_empty() {
        exact.sort_unstable();
        let spread = (exact[exact.len() - 1] - exact[0]) as u32;
        return (
            Some(Mic {
                log2: exact[exact.len() / 2],
                censor: Censor::Exact,
            }),
            spread,
        );
    }
    let at_most = mics
        .iter()
        .filter(|m| m.censor == Censor::AtMost)
        .map(|m| m.log2)
        .min();
    let above = mics
        .iter()
        .filter(|m| m.censor == Censor::Above)
        .map(|m| m.log2)
        .max();
    match (at_most, above) {
        (Some(v), None) => (
            Some(Mic {
                log2: v,
                censor: Censor::AtMost,
            }),
            0,
        ),
        (None, Some(v)) => (
            Some(Mic {
                log2: v,
                censor: Censor::Above,
            }),
            0,
        ),
        _ => (None, 0),
    }
}

/// Picks up to `cap` genomes per species for download. Each genome scores
/// the sum over its reported R/S calls of 1 / (size of that call's class),
/// so genomes carrying rare classes (often the resistant ones) and many
/// calls come first. A genome whose only labels were inferred from MIC
/// scores nothing and is not picked. Ties break on genome id, so the choice
/// is reproducible.
pub fn select_genomes(
    phenotypes: &[Phenotype],
    datasets: &[DatasetSpec],
    species: &str,
    cap: usize,
) -> Vec<String> {
    let targets: BTreeSet<&str> = datasets
        .iter()
        .filter(|d| d.species == species)
        .map(|d| d.antibiotic.as_str())
        .collect();
    let mut class_sizes: HashMap<(&str, Call), usize> = HashMap::new();
    let mut labels: BTreeMap<&str, Vec<(&str, Call)>> = BTreeMap::new();
    for p in phenotypes {
        if p.species != species || !targets.contains(p.antibiotic.as_str()) {
            continue;
        }
        if let Some(call) = p.reported_label() {
            *class_sizes.entry((&p.antibiotic, call)).or_default() += 1;
            labels
                .entry(&p.genome_id)
                .or_default()
                .push((&p.antibiotic, call));
        }
    }
    let mut scored: Vec<(f64, &str)> = labels
        .iter()
        .map(|(id, ls)| {
            (
                ls.iter()
                    .map(|(a, c)| 1.0 / class_sizes[&(*a, *c)] as f64)
                    .sum(),
                *id,
            )
        })
        .collect();
    scored.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(b.1)));
    scored
        .into_iter()
        .take(cap)
        .map(|(_, id)| id.to_string())
        .collect()
}

/// Names the design `plan_group_sampling` implements, recorded in each plan
/// so a saved inclusion probability can be traced to its derivation.
pub const GROUP_SAMPLING_DESIGN: &str = "group-round-robin-v1";

/// The generator, seeding and shuffle behind a plan's draw order, recorded in
/// each plan because the order is reproducible only on these versions.
/// `rand_core`'s `SeedableRng::seed_from_u64` expands the seed into the
/// ChaCha8 key, so its version matters as much as `rand_chacha`'s. A test
/// checks every version against Cargo.lock, and another pins literal draw
/// orders, so a dependency change that moves the draws fails rather than
/// drifts.
pub const GROUP_SAMPLING_RNG: &str = concat!(
    "rand 0.9.5 SliceRandom::shuffle, ",
    "rand_chacha 0.9.0 ChaCha8Rng, ",
    "rand_core 0.9.5 SeedableRng::seed_from_u64"
);

/// One genome a group-stratified plan may draw: a representative already
/// past QC and deduplication, with the lineage group it belongs to.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct EligibleGenome {
    pub genome_id: String,
    pub group_id: String,
    /// Whether the step that froze this list marked the genome as always
    /// included. When present it must agree with the declared stratum, so a
    /// frozen list and a separately written declaration cannot drift apart.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub always_included: Option<bool>,
}

/// The genomes a plan selects whatever the seed and cap, and the rule that
/// chose them. This is the only part of a plan that may come from phenotype
/// labels, and `label_dependent` records whether it did.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlwaysIncluded {
    pub genome_ids: Vec<String>,
    /// For example "every reported resistant carrier of a drug with fewer
    /// than 25 resistant genomes".
    pub reason: String,
    pub label_dependent: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Inclusion {
    pub genome_id: String,
    pub group_id: String,
    /// The chance under the design that this genome is selected at the
    /// plan's cap, derived in `plan_group_sampling`.
    pub probability: f64,
}

/// A label-free, group-stratified download plan: the opt-in alternative to
/// `select_genomes`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GroupSamplingPlan {
    pub design: String,
    /// `GROUP_SAMPLING_RNG` when the plan was made.
    pub rng: String,
    /// SHA-256 of the canonical eligible list: one `genome\tgroup\tflag\n`
    /// line per genome in genome id order, the flag `-`, `0` or `1`.
    pub input_sha256: String,
    pub seed: u64,
    pub cap: usize,
    /// The declared stratum, its ids sorted.
    pub always_included: AlwaysIncluded,
    /// The first `cap` genomes of the draw order, in that order.
    pub selected: Vec<String>,
    /// Every eligible genome, in genome id order.
    pub inclusion: Vec<Inclusion>,
}

/// Plans which `cap` genomes to download from `eligible` without reading a
/// phenotype label, so the selection cannot depend on any task's labels
/// except through the declared `always` stratum.
///
/// Plan one species per call, with only that species' genomes in `eligible`.
/// Group ids aren't tied to a species, and `EligibleGenome` carries none, so
/// a group id that two species share would make one group of both. Its
/// round-robin share and every `probability` in it would then be wrong, with
/// no error. The caller must check this; `docs/ablation-runner.md` states it
/// with the ablation runner's other obligations.
///
/// The draw order is every stratum genome in id order, then the remaining
/// genomes taken round-robin across their groups:
///
/// 1. Group the remaining genomes by `group_id`. Stratum genomes do not
///    count toward their group.
/// 2. Seed a ChaCha8 generator with `seed`. Shuffle each group's members,
///    starting from id order and taking groups in id order, then shuffle the
///    order of the groups. One stream serves every shuffle. Adding or
///    removing a group of two or more members moves every later group's
///    draws; any added or removed group can move the group order. A group of
///    one takes nothing from the stream. Resizing a group may or may not move
///    the others, depending on how `rand` spends random words on a shuffle.
/// 3. Round k takes the k-th member of every group with more than k
///    members, in the shuffled group order, until every group is exhausted.
///
/// The plan selects the first `cap` genomes of that order. The order does
/// not depend on `cap`, so under the same inputs and seed a smaller cap
/// selects a prefix of a larger cap's selection. Every group is drawn once
/// before any group is drawn twice: a huge group fills only the tail, and a
/// singleton is selected once the cap reaches every group's first draw.
/// When the cap leaves fewer slots than there are groups, the groups
/// represented are the first ones in the shuffled order.
///
/// `probability` is the chance a genome is selected at `cap` over the
/// design's randomization, a uniform group order and independent uniform
/// orders within groups, which the seeded shuffles stand in for. Stratum
/// genomes have 1. For the rest, let m = cap - |stratum| be the open slots.
/// A member of a group of n sits at within-group rank k with chance 1/n.
/// Before round k, C_k = sum_h min(n_h, k) genomes have been drawn and D_k
/// groups have more than k members. The member's group falls uniformly among
/// those D_k, so it is drawn with chance clamp(m - C_k, 0, D_k) / D_k, and
///
/// ```text
/// p = (1/n) sum_{k<n} clamp(m - C_k, 0, D_k) / D_k
/// ```
///
/// Members of one group share p, p summed over the non-stratum genomes is
/// min(m, their count), and a Horvitz-Thompson weight is 1/p.
///
/// Fails on a duplicate genome id, an empty genome or group id, one holding a
/// tab or line break, a stratum genome that is not eligible, a stratum
/// without a reason, a precomputed flag that disagrees with the stratum, or a
/// cap below the stratum's size.
pub fn plan_group_sampling(
    eligible: &[EligibleGenome],
    always: &AlwaysIncluded,
    cap: usize,
    seed: u64,
) -> Result<GroupSamplingPlan> {
    let mut by_id: BTreeMap<&str, &EligibleGenome> = BTreeMap::new();
    for g in eligible {
        for field in [&g.genome_id, &g.group_id] {
            ensure!(
                !field.is_empty() && !field.contains(['\t', '\n', '\r']),
                "eligible genome {:?} has an empty id or group, or one holding a tab or line break",
                g.genome_id
            );
        }
        ensure!(
            by_id.insert(&g.genome_id, g).is_none(),
            "genome {} is eligible twice",
            g.genome_id
        );
    }
    let stratum: BTreeSet<&str> = always.genome_ids.iter().map(String::as_str).collect();
    ensure!(
        stratum.len() == always.genome_ids.len(),
        "the always-included stratum lists a genome twice"
    );
    ensure!(
        stratum.is_empty() || !always.reason.trim().is_empty(),
        "the always-included stratum needs a reason"
    );
    for id in &stratum {
        ensure!(
            by_id.contains_key(id),
            "always-included genome {id} is not eligible"
        );
    }
    for g in by_id.values() {
        if let Some(flag) = g.always_included {
            ensure!(
                flag == stratum.contains(g.genome_id.as_str()),
                "genome {}'s always_included flag disagrees with the declared stratum",
                g.genome_id
            );
        }
    }
    ensure!(
        cap >= stratum.len(),
        "cap {cap} is below the {} always-included genomes",
        stratum.len()
    );

    let mut ctx = ring::digest::Context::new(&ring::digest::SHA256);
    for g in by_id.values() {
        let flag = match g.always_included {
            None => "-",
            Some(false) => "0",
            Some(true) => "1",
        };
        ctx.update(format!("{}\t{}\t{flag}\n", g.genome_id, g.group_id).as_bytes());
    }
    let input_sha256 = ctx
        .finish()
        .as_ref()
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect();

    let mut members: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for g in by_id.values() {
        if !stratum.contains(g.genome_id.as_str()) {
            members.entry(&g.group_id).or_default().push(&g.genome_id);
        }
    }
    let group_size: BTreeMap<&str, usize> = members.iter().map(|(g, m)| (*g, m.len())).collect();
    let mut groups: Vec<Vec<&str>> = members.into_values().collect();
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    for m in &mut groups {
        m.shuffle(&mut rng);
    }
    groups.shuffle(&mut rng);

    let mut order: Vec<&str> = stratum.iter().copied().collect();
    order.extend(round_robin(&groups));
    let by_size = round_robin_inclusion(
        &groups.iter().map(Vec::len).collect::<Vec<_>>(),
        cap - stratum.len(),
    );
    let inclusion = by_id
        .values()
        .map(|g| Inclusion {
            genome_id: g.genome_id.clone(),
            group_id: g.group_id.clone(),
            probability: if stratum.contains(g.genome_id.as_str()) {
                1.0
            } else {
                by_size[group_size[g.group_id.as_str()] - 1]
            },
        })
        .collect();
    Ok(GroupSamplingPlan {
        design: GROUP_SAMPLING_DESIGN.to_string(),
        rng: GROUP_SAMPLING_RNG.to_string(),
        input_sha256,
        seed,
        cap,
        always_included: AlwaysIncluded {
            genome_ids: stratum.iter().map(|id| id.to_string()).collect(),
            reason: always.reason.clone(),
            label_dependent: always.label_dependent,
        },
        selected: order.iter().take(cap).map(|id| id.to_string()).collect(),
        inclusion,
    })
}

/// Round k takes the k-th member of every group that has one, in group order.
fn round_robin<'a>(groups: &[Vec<&'a str>]) -> Vec<&'a str> {
    let mut order = Vec::with_capacity(groups.iter().map(Vec::len).sum());
    let mut active: Vec<&[&'a str]> = groups.iter().map(Vec::as_slice).collect();
    let mut k = 0;
    loop {
        active.retain(|m| m.len() > k);
        if active.is_empty() {
            return order;
        }
        order.extend(active.iter().map(|m| m[k]));
        k += 1;
    }
}

/// The design's inclusion probability for a non-stratum genome by the size of
/// its group: entry n - 1 is p for a group of n, given every group's size and
/// `open` slots after the stratum.
fn round_robin_inclusion(sizes: &[usize], open: usize) -> Vec<f64> {
    let deepest = sizes.iter().copied().max().unwrap_or(0);
    // D_k: the groups with more than k members.
    let mut deeper = vec![0usize; deepest];
    for &n in sizes {
        for d in &mut deeper[..n] {
            *d += 1;
        }
    }
    let mut by_size = Vec::with_capacity(deepest);
    let mut drawn = 0; // C_k
    let mut sum = 0.0;
    for (k, &d) in deeper.iter().enumerate() {
        sum += open.saturating_sub(drawn).min(d) as f64 / d as f64;
        drawn += d;
        by_size.push(sum / (k + 1) as f64);
    }
    by_size
}

#[cfg(test)]
mod tests {
    use super::*;

    fn amr(value: &str, sign: Option<&str>, unit: Option<&str>, method: &str) -> AmrRecord {
        AmrRecord {
            genome_id: "1.1".into(),
            antibiotic: "x".into(),
            measurement_value: Some(value.into()),
            measurement_sign: sign.map(Into::into),
            measurement_unit: unit.map(Into::into),
            laboratory_typing_method: Some(method.into()),
            ..Default::default()
        }
    }

    #[test]
    fn parses_mics_and_rejects_zone_diameters() {
        let m = parse_mic(&amr("0.25", Some("<="), Some("mg/L"), "Broth dilution")).unwrap();
        assert_eq!(
            m,
            Mic {
                log2: -2,
                censor: Censor::AtMost
            }
        );
        let m = parse_mic(&amr(">2/38", None, Some("mg/L"), "MIC")).unwrap();
        assert_eq!(
            m,
            Mic {
                log2: 1,
                censor: Censor::Above
            }
        );
        let m = parse_mic(&amr("0.06", None, None, "Agar dilution")).unwrap();
        assert_eq!(
            m,
            Mic {
                log2: -4,
                censor: Censor::Exact
            }
        );
        assert!(parse_mic(&amr("28", None, Some("mm"), "Disk diffusion")).is_none());
        assert!(parse_mic(&amr("28", None, None, "Disk diffusion")).is_none());
        // BV-BRC carries zone diameters from disk-diffusion studies as "mg/L".
        assert!(parse_mic(&amr("30", None, Some("mg/L"), "Disk diffusion")).is_none());
        let mut strip = amr("0.023", None, Some("mg/L"), "Disk diffusion");
        strip.laboratory_typing_platform = Some("E-test".into());
        assert_eq!(
            parse_mic(&strip),
            Some(Mic {
                log2: -5,
                censor: Censor::Exact
            })
        );
    }

    #[test]
    fn zone_diameters_do_not_become_mic_labels() {
        let record = |id: String, call: Option<&str>, value: &str, method: &str| AmrRecord {
            genome_id: id,
            antibiotic: "meropenem".into(),
            resistant_phenotype: call.map(Into::into),
            measurement_value: Some(value.into()),
            measurement_unit: Some("mg/L".into()),
            laboratory_typing_method: Some(method.into()),
            ..Default::default()
        };
        let mut records = Vec::new();
        for i in 0..30 {
            let (r, s) = (format!("562.{i}"), format!("562.{}", 100 + i));
            records.push(record(r, Some("Resistant"), "16", "Broth dilution"));
            records.push(record(s, Some("Susceptible"), "0.5", "Broth dilution"));
        }
        records.push(record("562.900".into(), None, "8", "Broth dilution"));
        // A 30 mm zone means susceptible; read as 30 mg/L it would label R.
        records.push(record("562.910".into(), None, "30", "Disk diffusion"));
        let genomes: Vec<GenomeRecord> = records
            .iter()
            .map(|r| GenomeRecord {
                genome_id: r.genome_id.clone(),
                species: Some("Escherichia coli".into()),
                genome_length: Some("5000000".into()),
                ..Default::default()
            })
            .collect();
        let mut resolver = Resolver::new(Default::default());
        let curated = curate(
            &records,
            &genomes,
            &QcConfig::default(),
            &BTreeMap::new(),
            &mut resolver,
        )
        .unwrap();

        let label = |id: &str| {
            curated
                .phenotypes
                .iter()
                .find(|p| p.genome_id == id)
                .map(|p| (p.label, p.label_from_mic))
        };
        assert_eq!(label("562.900"), Some((Some(Call::R), true)));
        assert_eq!(label("562.910"), None);
        assert_eq!(curated.funnel.records_no_measure, 1);
        let d = &curated.datasets[0];
        let counts = (d.n_resistant, d.n_susceptible, d.n_labels_from_mic);
        assert_eq!(counts, (30, 30, 1));
    }

    /// One species and drug: reported resistant calls at 16 mg/L, reported
    /// susceptible calls at 0.5 mg/L, and MIC-only records at 16 mg/L that
    /// the inferred breakpoint labels resistant.
    fn reported_and_inferred(
        resistant: usize,
        susceptible: usize,
        inferred: usize,
    ) -> (Vec<AmrRecord>, Vec<GenomeRecord>) {
        let record = |id: String, call: Option<&str>, value: &str| AmrRecord {
            genome_id: id,
            antibiotic: "meropenem".into(),
            resistant_phenotype: call.map(Into::into),
            measurement_value: Some(value.into()),
            measurement_unit: Some("mg/L".into()),
            laboratory_typing_method: Some("Broth dilution".into()),
            ..Default::default()
        };
        let mut records = Vec::new();
        records.extend((0..resistant).map(|i| record(format!("562.{i}"), Some("Resistant"), "16")));
        records.extend(
            (0..susceptible)
                .map(|i| record(format!("562.{}", 100 + i), Some("Susceptible"), "0.5")),
        );
        records.extend((0..inferred).map(|i| record(format!("562.{}", 200 + i), None, "16")));
        let genomes = records
            .iter()
            .map(|r| GenomeRecord {
                genome_id: r.genome_id.clone(),
                species: Some("Escherichia coli".into()),
                genome_length: Some("5000000".into()),
                ..Default::default()
            })
            .collect();
        (records, genomes)
    }

    fn curate_at(
        (records, genomes): &(Vec<AmrRecord>, Vec<GenomeRecord>),
        min_class: usize,
    ) -> Curated {
        let cfg = QcConfig {
            min_class,
            ..Default::default()
        };
        curate(
            records,
            genomes,
            &cfg,
            &BTreeMap::new(),
            &mut Resolver::new(Default::default()),
        )
        .unwrap()
    }

    #[test]
    fn mic_inferred_labels_do_not_admit_a_dataset() {
        let with_inferred = reported_and_inferred(24, 30, 10);

        // 34 resistant labels, 24 of them reported.
        let strict = curate_at(&with_inferred, 25);
        let inferred = strict.phenotypes.iter().filter(|p| p.label_from_mic);
        assert_eq!(inferred.filter(|p| p.mic.is_some()).count(), 10);
        assert!(strict.datasets.is_empty(), "{:?}", strict.datasets);

        let loose = curate_at(&with_inferred, 10);
        assert_eq!(loose.datasets.len(), 1);
        let d = &loose.datasets[0];
        let counts = (
            d.n_resistant,
            d.n_susceptible,
            d.n_labels_from_mic,
            d.n_with_mic,
        );
        assert_eq!(counts, (24, 30, 10, 64));
        assert_eq!(d.min_class, 10);

        // The reported calls alone make the same dataset.
        let reported = curate_at(&reported_and_inferred(24, 30, 0), 10);
        assert_eq!(reported.datasets.len(), 1);
        let d = &reported.datasets[0];
        let counts = (d.n_resistant, d.n_susceptible, d.n_labels_from_mic);
        assert_eq!(counts, (24, 30, 0));
    }

    #[test]
    fn species_check_failures_are_excluded_before_qc_and_from_class_counts() {
        let (records, genomes) = reported_and_inferred(30, 30, 0);
        let wrong = BTreeMap::from([("562.0".to_string(), "Mash distance 0.3".to_string())]);
        let curated = curate(
            &records,
            &genomes,
            &QcConfig::default(),
            &wrong,
            &mut Resolver::new(Default::default()),
        )
        .unwrap();
        let d = &curated.datasets[0];
        assert_eq!((d.n_resistant, d.n_susceptible), (29, 30));
        assert!(curated.phenotypes.iter().all(|p| p.genome_id != "562.0"));
        assert_eq!(
            curated.qc_failures["562.0"],
            format!("{SPECIES_CHECK_FAILED}: Mash distance 0.3")
        );
        assert_eq!(curated.funnel.genomes_failed_qc[SPECIES_CHECK_FAILED], 1);
        assert_eq!(curated.funnel.genomes_passing_qc, 59);
    }

    /// A failed genome is excluded for its species check even when its
    /// metadata names no species, which the next QC step would report.
    #[test]
    fn the_species_check_exclusion_comes_before_the_no_species_qc() {
        let (records, mut genomes) = reported_and_inferred(30, 30, 0);
        assert_eq!(genomes[0].genome_id, "562.0");
        genomes[0].species = None;
        let wrong = BTreeMap::from([("562.0".to_string(), "Mash distance 0.3".to_string())]);
        let curated = curate(
            &records,
            &genomes,
            &QcConfig::default(),
            &wrong,
            &mut Resolver::new(Default::default()),
        )
        .unwrap();
        assert!(curated.qc_failures["562.0"].starts_with(SPECIES_CHECK_FAILED));
        assert_eq!(curated.funnel.genomes_failed_qc.get("no species"), None);
        assert_eq!(curated.funnel.genomes_failed_qc[SPECIES_CHECK_FAILED], 1);
    }

    #[test]
    fn datasets_admit_a_class_of_exactly_min_class() {
        let exact = reported_and_inferred(24, 24, 0);
        let admitted = curate_at(&exact, 24);
        assert_eq!(admitted.datasets.len(), 1);
        let d = &admitted.datasets[0];
        assert_eq!((d.n_resistant, d.n_susceptible, d.min_class), (24, 24, 24));
        assert!(curate_at(&exact, 25).datasets.is_empty());
    }

    #[test]
    fn curated_files_without_min_class_read_as_the_old_learn_gate_of_25() {
        let old = serde_json::json!({
            "id": "a_b__x", "species": "A b", "species_slug": "a_b", "antibiotic": "x",
            "drug_class": "other", "n_resistant": 30, "n_susceptible": 40, "n_intermediate": 0,
            "n_with_mic": 0, "n_labels_from_mic": 0, "breakpoint": null,
        });
        let spec: DatasetSpec = serde_json::from_value(old.clone()).unwrap();
        assert_eq!(spec.min_class, 25);
        let mut custom = old;
        custom["min_class"] = 10.into();
        let spec: DatasetSpec = serde_json::from_value(custom).unwrap();
        assert_eq!(spec.min_class, 10);
    }

    #[test]
    fn evaluation_targets_require_unconflicted_reported_labels() {
        let mut p = Phenotype {
            genome_id: "fixture".into(),
            species: "fixture".into(),
            antibiotic: "fixture".into(),
            reported: None,
            conflict: false,
            mic: None,
            mic_spread: 0,
            label: Some(Call::R),
            label_from_mic: true,
            methods: vec![],
            standards: vec![],
            n_records: 1,
        };
        assert_eq!(p.reported_label(), None);
        p.reported = Some(Call::R);
        p.label_from_mic = false;
        assert_eq!(p.reported_label(), Some(Call::R));
        p.conflict = true;
        assert_eq!(p.reported_label(), None);
        p.conflict = false;
        p.reported = Some(Call::I);
        assert_eq!(p.reported_label(), None);
    }

    #[test]
    fn breakpoint_respects_censoring() {
        let mut pairs = Vec::new();
        for _ in 0..40 {
            pairs.push((
                Mic {
                    log2: -1,
                    censor: Censor::Exact,
                },
                Call::S,
            ));
            pairs.push((
                Mic {
                    log2: 0,
                    censor: Censor::AtMost,
                },
                Call::S,
            ));
            pairs.push((
                Mic {
                    log2: 3,
                    censor: Censor::Exact,
                },
                Call::R,
            ));
            pairs.push((
                Mic {
                    log2: 4,
                    censor: Censor::Above,
                },
                Call::R,
            ));
        }
        let bp = infer_breakpoint(&pairs).unwrap();
        assert!((0..=2).contains(&bp.r_above), "threshold {}", bp.r_above);
        assert_eq!(bp.agreement, 1.0);
        // "<= 4 mg/L" cannot be placed against a 2 mg/L threshold.
        assert_eq!(
            call_from_mic(
                Mic {
                    log2: 2,
                    censor: Censor::AtMost
                },
                1
            ),
            None
        );
        assert_eq!(
            call_from_mic(
                Mic {
                    log2: 2,
                    censor: Censor::Above
                },
                1
            ),
            Some(Call::R)
        );
    }

    #[test]
    fn selection_prefers_rare_reported_class_carriers() {
        let p = |id: &str, call, from_mic: bool| Phenotype {
            genome_id: id.into(),
            species: "A b".into(),
            antibiotic: "x".into(),
            reported: (!from_mic).then_some(call),
            conflict: false,
            mic: None,
            mic_spread: 0,
            label: Some(call),
            label_from_mic: from_mic,
            methods: vec![],
            standards: vec![],
            n_records: 1,
        };
        let mut ph: Vec<Phenotype> = (0..9)
            .map(|i| p(&format!("s{i}"), Call::S, false))
            .collect();
        // This resistant id sorts after every susceptible id, so equal weights
        // cannot accidentally produce the expected rare-class-first selection.
        ph.push(p("z-r", Call::R, false));
        // Counted, these MIC-inferred labels would grow the resistant class
        // to 11, so z-r would score 1/11, below each susceptible genome's 1/9.
        ph.extend((0..10).map(|i| p(&format!("m{i}"), Call::R, true)));
        let ds = vec![DatasetSpec {
            id: "a_b__x".into(),
            species: "A b".into(),
            species_slug: "a_b".into(),
            antibiotic: "x".into(),
            drug_class: "other".into(),
            n_resistant: 1,
            n_susceptible: 9,
            n_intermediate: 0,
            n_with_mic: 0,
            n_labels_from_mic: 10,
            breakpoint: None,
            min_class: 1,
        }];
        assert_eq!(select_genomes(&ph, &ds, "A b", 2), vec!["z-r", "s0"]);
        let all = select_genomes(&ph, &ds, "A b", 100);
        assert_eq!(all.len(), 10);
        assert!(all.iter().all(|id| !id.starts_with('m')), "{all:?}");
    }

    fn phenotype(
        id: &str,
        species: &str,
        antibiotic: &str,
        call: Call,
        conflict: bool,
        from_mic: bool,
    ) -> Phenotype {
        Phenotype {
            genome_id: id.into(),
            species: species.into(),
            antibiotic: antibiotic.into(),
            reported: (!from_mic).then_some(call),
            conflict,
            mic: None,
            mic_spread: 0,
            label: Some(call),
            label_from_mic: from_mic,
            methods: vec![],
            standards: vec![],
            n_records: 1,
        }
    }

    fn dataset(species: &str, antibiotic: &str) -> DatasetSpec {
        DatasetSpec {
            id: format!("{}__{antibiotic}", slug(species)),
            species: species.into(),
            species_slug: slug(species),
            antibiotic: antibiotic.into(),
            drug_class: "other".into(),
            n_resistant: 1,
            n_susceptible: 1,
            n_intermediate: 0,
            n_with_mic: 0,
            n_labels_from_mic: 0,
            breakpoint: None,
            min_class: 1,
        }
    }

    #[test]
    fn legacy_selection_writes_the_same_genomes_txt_bytes() {
        use super::Call::{I, R, S};
        let ph = vec![
            phenotype("g1", "A b", "x", R, false, false),
            phenotype("g1", "A b", "y", S, false, false),
            phenotype("g2", "A b", "x", S, false, false),
            phenotype("g2", "A b", "y", S, false, false),
            phenotype("g3", "A b", "x", S, false, false),
            phenotype("g4", "A b", "y", R, false, false),
            // A conflicted call, an intermediate call, a drug without a
            // dataset, an MIC-inferred label and another species: none counts.
            phenotype("g5", "A b", "x", R, true, false),
            phenotype("g6", "A b", "x", I, false, false),
            phenotype("g7", "A b", "z", R, false, false),
            phenotype("g8", "A b", "x", R, false, true),
            phenotype("g9", "C d", "x", S, false, false),
            phenotype("g10", "A b", "y", S, false, false),
            phenotype("g0", "A b", "x", S, false, false),
        ];
        let ds = vec![
            dataset("A b", "x"),
            dataset("A b", "y"),
            dataset("C d", "x"),
        ];
        // Classes: x/R {g1}, x/S {g0, g2, g3}, y/R {g4}, y/S {g1, g2, g10}.
        // Scores: g1 4/3, g4 1, g2 2/3, then g0, g10 and g3 tie at 1/3 and
        // break on the id as a string, so g10 sorts before g3.
        let bytes = |ids: Vec<String>| ids.join("\n") + "\n";
        assert_eq!(
            bytes(select_genomes(&ph, &ds, "A b", 100)),
            "g1\ng4\ng2\ng0\ng10\ng3\n"
        );
        assert_eq!(
            bytes(select_genomes(&ph, &ds, "A b", 4)),
            "g1\ng4\ng2\ng0\n"
        );
        assert_eq!(select_genomes(&ph, &ds, "C d", 100), vec!["g9"]);
        assert!(select_genomes(&ph, &ds, "A b", 0).is_empty());
    }

    #[test]
    fn legacy_selection_reads_only_the_requested_species() {
        use super::Call::{R, S};
        // Both species test drug x. Each A b genome scores 1, so the order
        // falls to the ids. Counting C d's two x/R calls in A b's classes
        // would drop g1 to 1/3 and move it last, and selecting them would
        // add h1 and h2.
        let ph = vec![
            phenotype("g1", "A b", "x", R, false, false),
            phenotype("g2", "A b", "y", R, false, false),
            phenotype("g3", "A b", "x", S, false, false),
            phenotype("h1", "C d", "x", R, false, false),
            phenotype("h2", "C d", "x", R, false, false),
        ];
        let ds = vec![
            dataset("A b", "x"),
            dataset("A b", "y"),
            dataset("C d", "x"),
        ];
        assert_eq!(select_genomes(&ph, &ds, "A b", 9), vec!["g1", "g2", "g3"]);
        assert_eq!(select_genomes(&ph, &ds, "A b", 1), vec!["g1"]);
        assert_eq!(select_genomes(&ph, &ds, "C d", 9), vec!["h1", "h2"]);
    }

    #[test]
    fn target_drugs_come_from_the_requested_species_datasets_only() {
        use super::Call::R;
        // F-1: A b has a dataset for x only; y is a C d dataset drug. g2's
        // y call is A b's own phenotype, but no A b dataset tests y, so it
        // scores nothing. Taking target drugs from every species' datasets
        // would select g2 as well.
        let ph = vec![
            phenotype("g1", "A b", "x", R, false, false),
            phenotype("g2", "A b", "y", R, false, false),
            phenotype("h1", "C d", "y", R, false, false),
        ];
        let ds = vec![dataset("A b", "x"), dataset("C d", "y")];
        assert_eq!(select_genomes(&ph, &ds, "A b", 9), vec!["g1"]);
        assert_eq!(select_genomes(&ph, &ds, "C d", 9), vec!["h1"]);
    }

    fn eligible(rows: &[(&str, &str)]) -> Vec<EligibleGenome> {
        rows.iter()
            .map(|&(genome, group)| EligibleGenome {
                genome_id: genome.into(),
                group_id: group.into(),
                always_included: None,
            })
            .collect()
    }

    fn stratum(ids: &[&str], label_dependent: bool) -> AlwaysIncluded {
        AlwaysIncluded {
            genome_ids: ids.iter().map(|id| id.to_string()).collect(),
            reason: "test stratum".into(),
            label_dependent,
        }
    }

    fn permutations<T: Clone>(items: &[T]) -> Vec<Vec<T>> {
        if items.is_empty() {
            return vec![vec![]];
        }
        (0..items.len())
            .flat_map(|i| {
                let mut rest = items.to_vec();
                let first = rest.remove(i);
                permutations(&rest).into_iter().map(move |mut p| {
                    p.insert(0, first.clone());
                    p
                })
            })
            .collect()
    }

    #[test]
    fn group_plans_are_deterministic_and_ignore_input_order() {
        let rows = eligible(&[
            ("a1", "A"),
            ("a2", "A"),
            ("a3", "A"),
            ("b1", "B"),
            ("c1", "C"),
            ("c2", "C"),
            ("d1", "D"),
            ("e1", "E"),
        ]);
        let plan = plan_group_sampling(&rows, &stratum(&[], false), 5, 7).unwrap();
        assert_eq!(
            plan,
            plan_group_sampling(&rows, &stratum(&[], false), 5, 7).unwrap()
        );
        let mut reversed = rows.clone();
        reversed.reverse();
        assert_eq!(
            plan,
            plan_group_sampling(&reversed, &stratum(&[], false), 5, 7).unwrap()
        );
        assert_eq!(
            (
                plan.design.as_str(),
                plan.rng.as_str(),
                plan.seed,
                plan.cap,
                plan.selected.len()
            ),
            (GROUP_SAMPLING_DESIGN, GROUP_SAMPLING_RNG, 7, 5, 5)
        );
        // sha256 of "a1\tA\t-\na2\tA\t-\n...e1\tE\t-\n", from sha256sum.
        assert_eq!(
            plan.input_sha256,
            "560226e4c19511b54baa61f56aee307f8732c63438f62a56d24e97d1c2599622"
        );
        let mut moved = rows.clone();
        moved[0].group_id = "B".into();
        assert_ne!(
            plan_group_sampling(&moved, &stratum(&[], false), 5, 7)
                .unwrap()
                .input_sha256,
            plan.input_sha256
        );
        let orders: BTreeSet<Vec<String>> = (0..8)
            .map(|seed| {
                plan_group_sampling(&rows, &stratum(&[], false), 8, seed)
                    .unwrap()
                    .selected
            })
            .collect();
        assert!(orders.len() > 1, "the seed never changed the order");
    }

    /// Groups a (a1-a3), b (b1-b6) and c (c1-c6), plus whatever `extra` adds.
    fn three_groups(extra: &[(&str, &str)]) -> Vec<EligibleGenome> {
        let mut rows: Vec<EligibleGenome> = [("a", 3), ("b", 6), ("c", 6)]
            .into_iter()
            .flat_map(|(group, n)| {
                (1..=n).map(move |i| EligibleGenome {
                    genome_id: format!("{group}{i}"),
                    group_id: group.into(),
                    always_included: None,
                })
            })
            .collect();
        rows.extend(eligible(extra));
        rows
    }

    #[test]
    fn group_plans_draw_a_pinned_order_on_the_recorded_generator() {
        // The independent review's order on rand 0.9.5, rand_chacha 0.9.0 and
        // rand_core 0.9.5.
        // A shuffle-order, stream, round-walk or dependency change that keeps
        // the design's distribution still moves these literals.
        let order = [
            "a1", "b2", "c1", "a3", "b3", "c3", "a2", "b5", "c2", "b6", "c5", "b4", "c4", "b1",
            "c6",
        ];
        let plan = plan_group_sampling(&three_groups(&[]), &stratum(&[], false), 15, 7).unwrap();
        assert_eq!(plan.selected, order);

        // a0 joins group a but is in the stratum, so group a's shuffle and
        // the stream are unchanged; s1's group holds only stratum genomes and
        // drops out. The stratum leads in id order and the cap of 10 cuts the
        // same order after eight draws.
        let rows = three_groups(&[("a0", "a"), ("s1", "s")]);
        let plan = plan_group_sampling(&rows, &stratum(&["s1", "a0"], false), 10, 7).unwrap();
        assert_eq!(plan.selected[..2], ["a0", "s1"]);
        assert_eq!(plan.selected[2..], order[..8]);
    }

    #[test]
    fn a_singleton_group_leaves_every_other_groups_draws_alone() {
        // The re-check's case: groups b and c of six, then a group a that
        // sorts first. `rand` shuffles one member without drawing, so a
        // singleton a can move only the group order, while a pair a takes a
        // word from the stream and moves b's shuffle.
        let rows = |a: usize| -> Vec<EligibleGenome> {
            [("a", a), ("b", 6), ("c", 6)]
                .into_iter()
                .flat_map(|(group, n)| {
                    (1..=n).map(move |i| EligibleGenome {
                        genome_id: format!("{group}{i}"),
                        group_id: group.into(),
                        always_included: None,
                    })
                })
                .collect()
        };
        let order = |a: usize, seed: u64| -> Vec<String> {
            plan_group_sampling(&rows(a), &stratum(&[], false), 12 + a, seed)
                .unwrap()
                .selected
                .into_iter()
                .filter(|id| !id.starts_with('a'))
                .collect()
        };
        // Each group's members in the order the plan draws them.
        let members = |order: &[String], group: char| -> Vec<String> {
            order
                .iter()
                .filter(|id| id.starts_with(group))
                .cloned()
                .collect()
        };
        let (mut group_order_moved, mut pair_moved_b) = (0, 0);
        for seed in 0..50 {
            let alone = order(0, seed);
            let singleton = order(1, seed);
            for group in ['b', 'c'] {
                assert_eq!(
                    members(&singleton, group),
                    members(&alone, group),
                    "a singleton moved group {group}'s draws at seed {seed}"
                );
            }
            group_order_moved += usize::from(singleton != alone);
            pair_moved_b += usize::from(members(&order(2, seed), 'b') != members(&alone, 'b'));
        }
        assert!(
            group_order_moved > 0,
            "a singleton never moved the group order"
        );
        assert!(pair_moved_b > 0, "a pair never moved group b's draws");
    }

    #[test]
    fn the_recorded_generator_matches_cargo_lock() {
        // A checkout with CRLF line endings must read the same entries.
        let lock = include_str!("../Cargo.lock").replace("\r\n", "\n");
        const CRATES_IO: &str = "registry+https://github.com/rust-lang/crates.io-index";
        for (name, version) in [
            ("rand", "0.9.5"),
            ("rand_chacha", "0.9.0"),
            ("rand_core", "0.9.5"),
        ] {
            assert!(
                lock.contains(&format!("name = \"{name}\"\nversion = \"{version}\"\n")),
                "Cargo.lock no longer pins {name} {version}"
            );
            assert!(GROUP_SAMPLING_RNG.contains(&format!("{name} {version} ")));
            // F-2: the entry found is the one GRM builds with. Cargo.toml asks
            // for 0.9 of rand and rand_chacha, and cargo resolves one
            // semver-compatible version of each, so the lock holds one entry
            // per crate, inside that requirement. That entry must be the
            // crates.io release: a path, git or patched source of the same
            // version could draw differently.
            assert!(
                lock.contains(&format!(
                    "name = \"{name}\"\nversion = \"{version}\"\nsource = \"{CRATES_IO}\"\n"
                )),
                "Cargo.lock's {name} {version} is not the crates.io release"
            );
            assert_eq!(
                lock.matches(&format!("name = \"{name}\"\n")).count(),
                1,
                "Cargo.lock holds more than one {name}"
            );
        }
        let manifest = include_str!("../Cargo.toml").replace("\r\n", "\n");
        for name in ["rand", "rand_chacha"] {
            assert!(
                manifest.contains(&format!("\n{name} = \"0.9\"\n")),
                "Cargo.toml no longer asks for {name} 0.9"
            );
        }
    }

    #[test]
    fn smaller_caps_select_a_prefix_of_larger_caps() {
        let rows = eligible(&[
            ("a1", "A"),
            ("a2", "A"),
            ("a3", "A"),
            ("a4", "A"),
            ("b1", "B"),
            ("b2", "B"),
            ("b3", "B"),
            ("c1", "C"),
            ("d1", "D"),
            ("e1", "E"),
        ]);
        let always = stratum(&["b2"], false);
        for seed in 0..4 {
            let full = plan_group_sampling(&rows, &always, rows.len(), seed)
                .unwrap()
                .selected;
            assert_eq!(full.len(), rows.len());
            for cap in 1..=rows.len() {
                let plan = plan_group_sampling(&rows, &always, cap, seed).unwrap();
                assert_eq!(plan.selected, full[..cap], "seed {seed}, cap {cap}");
            }
        }
    }

    #[test]
    fn permuting_labels_moves_legacy_selection_but_not_the_group_plan() {
        let ids = ["g0", "g1", "g2", "g3", "g4", "g5"];
        let resistant_at = |r: usize| -> Vec<Phenotype> {
            ids.iter()
                .enumerate()
                .map(|(i, id)| {
                    let call = if i == r { Call::R } else { Call::S };
                    phenotype(id, "A b", "x", call, false, false)
                })
                .collect()
        };
        // The planner sees only genome ids and groups, whatever labels sit
        // beside them; no label can reach it.
        let plan_for = |ph: &[Phenotype]| {
            let rows: Vec<EligibleGenome> = ph
                .iter()
                .map(|p| {
                    let group = if p.genome_id.as_str() < "g3" {
                        "A"
                    } else {
                        "B"
                    };
                    EligibleGenome {
                        genome_id: p.genome_id.clone(),
                        group_id: group.into(),
                        always_included: None,
                    }
                })
                .collect();
            plan_group_sampling(&rows, &stratum(&[], false), 2, 3).unwrap()
        };
        let ds = vec![dataset("A b", "x")];
        let (before, after) = (resistant_at(0), resistant_at(4));
        assert_eq!(select_genomes(&before, &ds, "A b", 1), vec!["g0"]);
        assert_eq!(select_genomes(&after, &ds, "A b", 1), vec!["g4"]);
        assert_eq!(plan_for(&before), plan_for(&after));
    }

    #[test]
    fn the_always_included_stratum_is_always_selected_and_recorded() {
        let mut rows = eligible(&[
            ("a1", "A"),
            ("a2", "A"),
            ("b1", "B"),
            ("r1", "A"),
            ("r2", "C"),
        ]);
        rows[0].always_included = Some(false);
        rows[3].always_included = Some(true);
        let always = AlwaysIncluded {
            genome_ids: vec!["r2".into(), "r1".into()],
            reason: "reported resistant to x".into(),
            label_dependent: true,
        };
        for seed in 0..4 {
            for cap in 2..=rows.len() {
                let plan = plan_group_sampling(&rows, &always, cap, seed).unwrap();
                assert_eq!(plan.selected[..2], ["r1", "r2"]);
                let distinct: BTreeSet<&String> = plan.selected.iter().collect();
                assert_eq!(distinct.len(), cap, "{:?}", plan.selected);
                assert_eq!(
                    plan.always_included,
                    AlwaysIncluded {
                        genome_ids: vec!["r1".into(), "r2".into()],
                        ..always.clone()
                    }
                );
                for inc in plan
                    .inclusion
                    .iter()
                    .filter(|i| i.genome_id.starts_with('r'))
                {
                    assert_eq!(inc.probability, 1.0);
                }
            }
        }

        // The flags reach the hash: sha256 of "a1\tA\t0\na2\tA\t-\nb1\tB\t-\n
        // r1\tA\t1\nr2\tC\t-\n", from sha256sum. Dropping a flag changes it.
        let hash = |rows: &[EligibleGenome]| {
            plan_group_sampling(rows, &always, 2, 0)
                .unwrap()
                .input_sha256
        };
        assert_eq!(
            hash(&rows),
            "d90ada261f475a5126661f96b1dcbd8c750efb26d21f1e1ebafe2279e417f371"
        );
        for i in [0, 3] {
            let mut unflagged = rows.clone();
            unflagged[i].always_included = None;
            assert_ne!(hash(&unflagged), hash(&rows));
        }

        let error = |rows: &[EligibleGenome], always: &AlwaysIncluded, cap| {
            plan_group_sampling(rows, always, cap, 0)
                .unwrap_err()
                .to_string()
        };
        assert!(error(&rows, &always, 1).contains("below"));
        let unknown = AlwaysIncluded {
            genome_ids: vec!["zz".into()],
            ..always.clone()
        };
        assert!(error(&rows, &unknown, 3).contains("not eligible"));
        let twice = AlwaysIncluded {
            genome_ids: vec!["r1".into(), "r1".into()],
            ..always.clone()
        };
        assert!(error(&rows, &twice, 3).contains("twice"));
        let unexplained = AlwaysIncluded {
            reason: " ".into(),
            ..always.clone()
        };
        assert!(error(&rows, &unexplained, 3).contains("reason"));
        let mut disagreeing = rows.clone();
        disagreeing[3].always_included = Some(false);
        assert!(error(&disagreeing, &always, 3).contains("disagrees"));
    }

    #[test]
    fn group_plans_reject_ids_the_canonical_list_cannot_record() {
        let error = |rows: &[(&str, &str)]| {
            plan_group_sampling(&eligible(rows), &stratum(&[], false), 1, 0)
                .unwrap_err()
                .to_string()
        };
        assert!(error(&[("a", "A"), ("a", "B")]).contains("twice"));
        for row in [("a\tb", "A"), ("a\nb", "A"), ("a\rb", "A"), ("a", "A\tB")] {
            assert!(error(&[row]).contains("tab or line break"), "{row:?}");
        }
        for row in [("", "A"), ("a", "")] {
            assert!(error(&[row]).contains("empty"), "{row:?}");
        }
    }

    #[test]
    fn inclusion_probabilities_match_an_exhaustive_enumeration() {
        // Group a holds three genomes besides stratum genome s, b one, c two.
        let rows = eligible(&[
            ("a0", "a"),
            ("a1", "a"),
            ("a2", "a"),
            ("b0", "b"),
            ("c0", "c"),
            ("c1", "c"),
            ("s", "a"),
        ]);
        let always = stratum(&["s"], false);
        let groups = [vec!["a0", "a1", "a2"], vec!["b0"], vec!["c0", "c1"]];
        // Every group order with every order within each group: 6 x (6 x 1 x
        // 2) = 72 equally likely draw orders.
        let mut outcomes = Vec::new();
        for order in permutations(&[0usize, 1, 2]) {
            for a in permutations(&groups[0]) {
                for b in permutations(&groups[1]) {
                    for c in permutations(&groups[2]) {
                        let within = [a.clone(), b.clone(), c.clone()];
                        let ordered: Vec<Vec<&str>> =
                            order.iter().map(|&g| within[g].clone()).collect();
                        outcomes.push(round_robin(&ordered));
                    }
                }
            }
        }
        assert_eq!(outcomes.len(), 72);
        for open in 0..=6 {
            let plan = plan_group_sampling(&rows, &always, 1 + open, 0).unwrap();
            for inc in &plan.inclusion {
                let expected = if inc.genome_id == "s" {
                    1.0
                } else {
                    let hits = outcomes
                        .iter()
                        .filter(|o| o[..open].contains(&inc.genome_id.as_str()))
                        .count();
                    hits as f64 / 72.0
                };
                assert!(
                    (inc.probability - expected).abs() < 1e-12,
                    "{} with {open} open slots: {} vs {expected}",
                    inc.genome_id,
                    inc.probability
                );
            }
            let total: f64 = plan.inclusion.iter().map(|i| i.probability).sum();
            assert!((total - (1 + open) as f64).abs() < 1e-12);
        }

        // With two open slots: 2/9 in the group of three, 2/3 for the
        // singleton, 1/3 in the pair.
        let plan = plan_group_sampling(&rows, &always, 3, 0).unwrap();
        let expected = |id: &str| match &id[..1] {
            "a" => 2.0 / 9.0,
            "b" => 2.0 / 3.0,
            "c" => 1.0 / 3.0,
            _ => 1.0,
        };
        for inc in &plan.inclusion {
            assert!((inc.probability - expected(&inc.genome_id)).abs() < 1e-12);
        }
        // The seeded shuffles draw at those rates.
        let seeds = 4000;
        let mut hits: BTreeMap<String, usize> = BTreeMap::new();
        for seed in 0..seeds {
            for id in plan_group_sampling(&rows, &always, 3, seed)
                .unwrap()
                .selected
            {
                *hits.entry(id).or_default() += 1;
            }
        }
        for inc in &plan.inclusion {
            let rate = hits.get(&inc.genome_id).copied().unwrap_or(0) as f64 / seeds as f64;
            assert!(
                (rate - inc.probability).abs() < 0.04,
                "{}: drawn at {rate}, designed {}",
                inc.genome_id,
                inc.probability
            );
        }
    }

    #[test]
    fn group_plans_handle_full_caps_one_group_huge_groups_and_no_input() {
        let none = stratum(&[], false);
        let empty = plan_group_sampling(&[], &none, 10, 1).unwrap();
        assert!(empty.selected.is_empty() && empty.inclusion.is_empty());
        assert_eq!(
            empty.input_sha256,
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );

        let rows = eligible(&[("a1", "A"), ("a2", "A"), ("b1", "B"), ("c1", "C")]);
        for cap in [rows.len(), rows.len() + 5] {
            let plan = plan_group_sampling(&rows, &none, cap, 2).unwrap();
            let mut selected = plan.selected.clone();
            selected.sort();
            assert_eq!(selected, ["a1", "a2", "b1", "c1"]);
            assert!(plan.inclusion.iter().all(|i| i.probability == 1.0));
        }

        // One group: a seeded permutation of it, each genome at cap / n.
        let one = eligible(&[("x1", "G"), ("x2", "G"), ("x3", "G"), ("x4", "G")]);
        let plan = plan_group_sampling(&one, &none, 3, 5).unwrap();
        assert_eq!(plan.selected.len(), 3);
        assert!(plan.inclusion.iter().all(|i| i.probability == 0.75));

        // Every group's first draw comes before any second draw, so a huge
        // group cannot crowd out singletons.
        let mut huge: Vec<EligibleGenome> = (0..100)
            .map(|i| EligibleGenome {
                genome_id: format!("h{i:03}"),
                group_id: "H".into(),
                always_included: None,
            })
            .collect();
        huge.extend(eligible(&[("s1", "S1"), ("s2", "S2")]));
        for seed in 0..8 {
            let plan = plan_group_sampling(&huge, &none, 3, seed).unwrap();
            let selected: BTreeSet<&str> = plan.selected.iter().map(String::as_str).collect();
            assert!(
                selected.contains("s1") && selected.contains("s2"),
                "{selected:?}"
            );
            for inc in &plan.inclusion {
                let expected = if inc.group_id == "H" { 0.01 } else { 1.0 };
                assert!((inc.probability - expected).abs() < 1e-12);
            }
        }

        // Fewer open slots than singleton groups: each group at slots / groups,
        // and no group drawn twice.
        let singles = eligible(&[("a", "A"), ("b", "B"), ("c", "C"), ("d", "D")]);
        let plan = plan_group_sampling(&singles, &none, 2, 4).unwrap();
        assert_eq!(plan.selected.len(), 2);
        assert!(plan.inclusion.iter().all(|i| i.probability == 0.5));
    }
}
