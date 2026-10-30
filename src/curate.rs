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
use anyhow::Result;
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
}
