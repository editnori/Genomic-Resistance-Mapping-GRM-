//! Dashboard data: aggregates curated records, species structure, and
//! per-dataset reports into the files the web app reads.
//!
//!   index.json              manifest: sources, funnel, species, dataset catalog,
//!                           cross-dataset panel summary, Kover benchmark
//!   corpus.json             record-level distributions, drug-name resolution,
//!                           inferred breakpoints
//!   species/\<slug>.json    genome map, dendrogram, resistance strings per
//!                           genome, cross-resistance pairs
//!   datasets/\<id>.json     the per-dataset report as `grm learn` wrote it

use crate::bvbrc::{self, AmrRecord, GenomeRecord};
use crate::curate::{Breakpoint, Call, DatasetSpec, Funnel, Phenotype};
use crate::drugs;
use crate::layout::Layout;
use crate::learn::NoCall;
use crate::species::{self, Structure};
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;

/// Schema 2 lets a catalog row have no best learner and a learner row have
/// null rates when it called no isolate. A schema-1 reader stops at its
/// version check instead of reading either as a number.
pub const SCHEMA: u32 = 2;

/// Replace a complete file so a dashboard request cannot read partial JSON.
pub(crate) fn publish(path: &Path, bytes: &[u8]) -> Result<()> {
    let pending = path.with_extension(format!("json.{}.tmp", std::process::id()));
    std::fs::write(&pending, bytes)?;
    std::fs::rename(&pending, path)?;
    Ok(())
}

/// Summary rates. Panel and catalog rates are always present; a learner
/// row's are None when it called no isolate, so an empty call set is written
/// as null rather than 0% error.
#[derive(Serialize, Deserialize, Clone, Copy, Default)]
pub struct Lite {
    pub balanced_accuracy: Option<f64>,
    pub sensitivity: Option<f64>,
    pub specificity: Option<f64>,
    pub vme: Option<f64>,
    pub me: Option<f64>,
    pub auroc: Option<f64>,
    /// A learner row's samples in null-fit folds. Before revision 14 the
    /// rates above cover the rest; see `null_fits_scored`. Absent on panel
    /// metrics and on rows written before no-calls.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub no_call: Option<NoCall>,
    /// True on a revision-14 learner row with a null-fit fold: its `no_call`
    /// samples are scored at their fold's training prevalence and counted in
    /// the rates above. Absent (false) on every other row.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub null_fits_scored: bool,
}

fn lite(m: &Value) -> Lite {
    let f = |k: &str| Some(m.get(k).and_then(Value::as_f64).unwrap_or(0.0));
    Lite {
        balanced_accuracy: f("balanced_accuracy"),
        sensitivity: f("sensitivity"),
        specificity: f("specificity"),
        vme: f("vme"),
        me: f("me"),
        auroc: m.get("auroc").and_then(Value::as_f64),
        no_call: None,
        null_fits_scored: false,
    }
}

/// A learner row's rates and no-calls. Null metrics mean a row written
/// before revision 14 called no isolate, so its rates stay None. A row that
/// carries `called_metrics` (revision 14) scored its null-fit samples.
fn learner_lite(l: &Value) -> Lite {
    let no_call = NoCall::deserialize(&l["no_call"]).ok();
    let null_fits_scored = l.get("called_metrics").is_some();
    match &l["metrics"] {
        Value::Null => Lite {
            no_call,
            null_fits_scored,
            ..Lite::default()
        },
        m => Lite {
            no_call,
            null_fits_scored,
            ..lite(m)
        },
    }
}

/// Whether a learner row's rates cover every isolate. A revision-14 row
/// always does: a null-fit fold is scored at its training prevalence. A row
/// of revisions 7 to 13 does only when it left no isolate uncalled. Rows
/// written before no-calls scored every isolate, so they pass.
fn scored_on_every_isolate(m: &Lite) -> bool {
    m.balanced_accuracy.is_some()
        && (m.null_fits_scored || m.no_call.is_none_or(|x| x.resistant + x.susceptible == 0))
}

/// Every learner's rates on both splits of one report's `splits`, in
/// pipeline order.
pub(crate) fn learner_rows(splits: &Value) -> Vec<LearnerRow> {
    let scheme = |name: &str| splits.as_array()?.iter().find(|s| s["scheme"] == name);
    let lineage = scheme("lineage");
    let lineage_ran = lineage.is_some_and(|s| s["unavailable_reason"].is_null());
    crate::learn::eval::LEARNERS
        .iter()
        .map(|&l| LearnerRow {
            learner: l.into(),
            random: learner_row(scheme("random"), l).map(learner_lite),
            lineage: learner_row(lineage, l).map(learner_lite),
            lineage_ran,
        })
        .collect()
}

/// The best learner and its random-split balanced accuracy: the highest
/// score among learners whose rates cover every isolate on every split that
/// ran. From revision 14 every learner's do, so all five compete; in an
/// older report a learner that left isolates uncalled is passed over. None
/// when no learner qualifies. A tie goes to the later learner.
pub(crate) fn best_learner(rows: &[LearnerRow]) -> Option<(&LearnerRow, f64)> {
    rows.iter()
        .filter(|l| !l.lineage_ran || l.lineage.as_ref().is_some_and(scored_on_every_isolate))
        .filter_map(|l| {
            let random = l.random.filter(scored_on_every_isolate)?;
            Some((l, random.balanced_accuracy?))
        })
        .max_by(|a, b| a.1.total_cmp(&b.1))
}

#[derive(Serialize, Deserialize, Clone, Copy)]
pub struct Confusion {
    pub n: usize,
    pub tp: usize,
    pub fp: usize,
    pub tn: usize,
    #[serde(rename = "fn")]
    pub fn_: usize,
}

fn learner_row<'a>(split: Option<&'a Value>, learner: &str) -> Option<&'a Value> {
    split?["learners"]
        .as_array()?
        .iter()
        .find(|l| l["learner"] == learner)
}

#[derive(Serialize)]
pub struct LearnerRow {
    pub learner: String,
    pub random: Option<Lite>,
    pub lineage: Option<Lite>,
    /// Whether the lineage split ran. A missing learner row then means the
    /// learner did not cover that split, rather than that the split was absent.
    #[serde(skip)]
    lineage_ran: bool,
}

#[derive(Clone, Serialize)]
pub struct KoverRow {
    pub sensitivity: f64,
    pub specificity: f64,
    pub balanced_accuracy: f64,
    pub n_rules: f64,
    pub n_examples: f64,
}

#[derive(Serialize)]
pub struct CatalogRow {
    pub id: String,
    pub species: String,
    pub species_slug: String,
    pub antibiotic: String,
    pub drug_class: String,
    pub evaluation: Value,
    pub n: usize,
    pub n_resistant: usize,
    pub labels_from_mic: usize,
    /// The report's `counts.labels_from_mic_excluded`; absent for reports
    /// that do not record it.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub labels_from_mic_excluded: Option<usize>,
    /// The report's `counts.min_class`; 25 for a report that does not record
    /// it, the gate `grm learn` applied before reports recorded one.
    pub min_class: usize,
    pub features: usize,
    pub lineages: usize,
    pub best_learner: Option<String>,
    pub best_random: Option<f64>,
    pub best_lineage: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lineage_unavailable_reason: Option<String>,
    pub learners: Vec<LearnerRow>,
    pub catalog: Option<Lite>,
    pub mic_essential_agreement: Option<f64>,
    pub scm_rules: Vec<String>,
    pub panels: Vec<(String, Lite)>,
    pub kover: Option<KoverRow>,
}

#[derive(Serialize)]
pub struct SpeciesRow {
    pub species: String,
    pub slug: String,
    pub datasets_curated: usize,
    pub datasets_reported: usize,
    pub genomes_labeled: usize,
    pub genomes_analyzed: usize,
    pub lineages: Option<usize>,
    pub lineage_threshold: Option<f32>,
}

#[derive(Serialize)]
pub struct PanelSummaryRow {
    pub panel: String,
    pub name: String,
    pub dataset: String,
    pub species: String,
    pub antibiotic: String,
    pub metrics: Lite,
    pub best_model: Lite,
    pub best_learner: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub counts: Option<Confusion>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub best_model_counts: Option<Confusion>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub no_call: Option<NoCall>,
}

#[derive(Serialize)]
pub struct UnsupportedPanelRow {
    pub panel: String,
    pub name: String,
    pub dataset: String,
    pub species: String,
    pub antibiotic: String,
    pub reason: String,
}

#[derive(Serialize)]
pub struct Manifest {
    pub schema: u32,
    pub generated_at: String,
    pub evaluation_revision: u32,
    /// The `counts.k` every published dataset report records, or None when
    /// two reports differ, one predates the field, or none is published.
    pub k: Option<usize>,
    pub amr_records: usize,
    pub amrfinder_version: String,
    pub name_vocabulary: Option<String>,
    pub funnel: Funnel,
    pub species: Vec<SpeciesRow>,
    pub datasets: Vec<CatalogRow>,
    pub panels: Vec<PanelSummaryRow>,
    pub unsupported_panels: Vec<UnsupportedPanelRow>,
}

fn rule_text(report: &Value) -> Vec<String> {
    let feats: HashMap<u64, &Value> = report["features"]
        .as_array()
        .map(|a| {
            a.iter()
                .filter_map(|f| f["feature"].as_u64().map(|id| (id, f)))
                .collect()
        })
        .unwrap_or_default();
    report["models"]["scm"]["rules"]
        .as_array()
        .map(|rules| {
            rules
                .iter()
                .map(|r| {
                    let id = r["feature"].as_u64().unwrap_or(0);
                    let presence = r["presence"].as_bool().unwrap_or(true);
                    let f = feats.get(&id);
                    let label = f
                        .and_then(|f| {
                            let families: Vec<&str> = f["catalog_families"]
                                .as_array()
                                .map(|a| a.iter().filter_map(Value::as_str).collect())
                                .unwrap_or_default();
                            crate::annotate::family_label(&families).or_else(|| {
                                let locus = &f["locus"];
                                locus["gene"]
                                    .as_str()
                                    .map(String::from)
                                    .or_else(|| locus["product"].as_str().map(String::from))
                            })
                        })
                        .unwrap_or_else(|| format!("feature {id}"));
                    format!(
                        "{} {label}",
                        if presence {
                            "presence of"
                        } else {
                            "absence of"
                        }
                    )
                })
                .collect()
        })
        .unwrap_or_default()
}

/// What a missing Kover file tells `grm report`'s user, who can pass `--kover`.
pub const REPORT_KOVER_REMEDY: &str = "run grm from the source folder, which holds \
     reference/kover-summary.json, or give `grm report --kover` the file's path";

/// What it tells `grm pipeline`'s, which reads only the default path.
pub const PIPELINE_KOVER_REMEDY: &str = "run `grm pipeline` from the source folder, which holds \
     reference/kover-summary.json; pipeline reads the Kover file only from there";

/// The historical Kover rows `index.json` shows beside each dataset, keyed by
/// species and canonical antibiotic. A file that is missing or doesn't parse
/// is an error: the default path is relative to the working folder, and
/// reading nothing would drop every dataset's `kover` field without a word.
/// `remedy` is what a missing file tells the user to do.
pub fn kover_benchmark(path: &Path, remedy: &str) -> Result<HashMap<(String, String), KoverRow>> {
    let shown = match std::path::absolute(path) {
        Ok(abs) if abs != path => format!("{} ({})", path.display(), abs.display()),
        _ => path.display().to_string(),
    };
    let bytes = std::fs::read(path)
        .with_context(|| format!("read the Kover comparison {shown}: {remedy}"))?;
    let rows: Vec<Value> = serde_json::from_slice(&bytes)
        .with_context(|| format!("parse the Kover comparison {shown}"))?;
    let vocab = drugs::Vocabulary::default();
    Ok(rows
        .iter()
        .filter_map(|r| {
            let species = r["species"].as_str()?.to_string();
            let drug = drugs::resolve(r["antibiotic"].as_str()?, &vocab);
            (drug.kind == drugs::NameKind::Drug).then(|| {
                let (se, sp) = (
                    r["sensitivity"].as_f64().unwrap_or(0.0),
                    r["specificity"].as_f64().unwrap_or(0.0),
                );
                (
                    (species, drug.canonical),
                    KoverRow {
                        sensitivity: se,
                        specificity: sp,
                        balanced_accuracy: (se + sp) / 2.0,
                        n_rules: r["n_rules"].as_f64().unwrap_or(0.0),
                        n_examples: r["ds_n_examples"].as_f64().unwrap_or(0.0),
                    },
                )
            })
        })
        .collect())
}

#[derive(Serialize)]
struct GenomeRow<'a> {
    id: &'a str,
    lineage: u32,
    x: f32,
    y: f32,
    nearest: f32,
    country: Option<&'a str>,
    year: Option<&'a str>,
    host: Option<&'a str>,
    mlst: Option<&'a str>,
    length: u64,
    contigs: u32,
    n50: u64,
    groups: Vec<u16>,
    /// One character per drug in `drugs`: R, S, I, or '.' (untested).
    profile: String,
}

#[derive(Serialize)]
struct CrossPair {
    a: String,
    b: String,
    n: usize,
    rr: usize,
    rs: usize,
    sr: usize,
    ss: usize,
    phi: f64,
}

fn species_file(
    layout: &Layout,
    slug: &str,
    species: &str,
    datasets: &[&DatasetSpec],
    phenotypes: &[&Phenotype],
    genomes: &HashMap<&str, &GenomeRecord>,
) -> Result<Option<(usize, usize, f32)>> {
    let path = layout.species_dir(slug).join("structure.json");
    let structure_bytes = match std::fs::read(&path) {
        Ok(bytes) => bytes,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(None),
        Err(e) => return Err(e).with_context(|| format!("read {}", path.display())),
    };
    let structure: Structure = serde_json::from_slice(&structure_bytes)?;
    let counted = species::load_counted(layout, slug)?;
    let drugs: Vec<String> = datasets.iter().map(|d| d.antibiotic.clone()).collect();
    let drug_index: HashMap<&str, usize> = drugs
        .iter()
        .enumerate()
        .map(|(i, d)| (d.as_str(), i))
        .collect();
    let mut profiles: HashMap<&str, Vec<u8>> = HashMap::new();
    for p in phenotypes {
        if let (Some(&j), Some(call)) = (drug_index.get(p.antibiotic.as_str()), p.label) {
            let row = profiles
                .entry(p.genome_id.as_str())
                .or_insert_with(|| vec![b'.'; drugs.len()]);
            row[j] = match call {
                Call::R => b'R',
                Call::S => b'S',
                Call::I => b'I',
            };
        }
    }
    let mut dict: BTreeMap<String, u16> = BTreeMap::new();
    for g in &counted.genomes {
        for h in &g.genes {
            let next = dict.len() as u16;
            dict.entry(h.group.clone()).or_insert(next);
        }
    }
    let counted_by: HashMap<&str, &species::Counted> = counted
        .genomes
        .iter()
        .map(|g| (g.genome_id.as_str(), g))
        .collect();
    let rows: Vec<GenomeRow> = structure
        .genomes
        .iter()
        .enumerate()
        .map(|(i, id)| {
            let meta = genomes.get(id.as_str());
            let c = counted_by.get(id.as_str());
            GenomeRow {
                id,
                lineage: structure.lineages.labels[i],
                x: structure.coords[i][0],
                y: structure.coords[i][1],
                nearest: structure.nearest[i],
                country: meta.and_then(|m| m.isolation_country.as_deref()),
                year: meta.and_then(|m| m.collection_year.as_deref()),
                host: meta.and_then(|m| m.host_name.as_deref()),
                mlst: meta.and_then(|m| m.mlst.as_deref()),
                length: c.map_or(0, |c| c.stats.length),
                contigs: c.map_or(0, |c| c.stats.contigs),
                n50: c.map_or(0, |c| c.stats.n50),
                groups: c.map_or(vec![], |c| c.genes.iter().map(|h| dict[&h.group]).collect()),
                profile: String::from_utf8(
                    profiles
                        .get(id.as_str())
                        .cloned()
                        .unwrap_or_else(|| vec![b'.'; drugs.len()]),
                )
                .unwrap(),
            }
        })
        .collect();
    let mut cross = Vec::new();
    for a in 0..drugs.len() {
        for b in a + 1..drugs.len() {
            let (mut rr, mut rs, mut sr, mut ss) = (0, 0, 0, 0);
            for row in &rows {
                let p = row.profile.as_bytes();
                match (p[a], p[b]) {
                    (b'R', b'R') => rr += 1,
                    (b'R', b'S') => rs += 1,
                    (b'S', b'R') => sr += 1,
                    (b'S', b'S') => ss += 1,
                    _ => {}
                }
            }
            let n = rr + rs + sr + ss;
            if n < 30 {
                continue;
            }
            let (rr_f, rs_f, sr_f, ss_f) = (rr as f64, rs as f64, sr as f64, ss as f64);
            let denom = ((rr_f + rs_f) * (sr_f + ss_f) * (rr_f + sr_f) * (rs_f + ss_f)).sqrt();
            let phi = if denom == 0.0 {
                0.0
            } else {
                (rr_f * ss_f - rs_f * sr_f) / denom
            };
            cross.push(CrossPair {
                a: drugs[a].clone(),
                b: drugs[b].clone(),
                n,
                rr,
                rs,
                sr,
                ss,
                phi,
            });
        }
    }
    let mut catalog_groups = vec![String::new(); dict.len()];
    for (g, i) in dict {
        catalog_groups[i as usize] = g;
    }
    let out = serde_json::json!({
        "schema": SCHEMA,
        "species": species,
        "slug": slug,
        "drugs": drugs,
        "drug_classes": datasets.iter().map(|d| d.drug_class.clone()).collect::<Vec<_>>(),
        "catalog_groups": catalog_groups,
        "lineage_threshold": structure.lineages.threshold,
        "n_lineages": structure.lineages.n_clusters,
        "leaf_order": structure.leaf_order,
        "merges": structure.merges,
        "genomes": rows,
        "cross": cross,
    });
    let path = layout
        .root
        .join("report/species")
        .join(format!("{slug}.json"));
    crate::layout::ensure_parent(&path)?;
    publish(&path, &serde_json::to_vec(&out)?)?;
    Ok(Some((
        structure.genomes.len(),
        structure.lineages.n_clusters,
        structure.lineages.threshold,
    )))
}

fn top_counts<'a>(items: impl Iterator<Item = &'a str>, k: usize) -> Vec<(String, usize)> {
    let mut c: HashMap<&str, usize> = HashMap::new();
    for i in items {
        *c.entry(i).or_default() += 1;
    }
    let mut v: Vec<(String, usize)> = c.into_iter().map(|(a, b)| (a.to_string(), b)).collect();
    v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    v.truncate(k);
    v
}

/// Builds every dashboard file under `<data>/report/` and returns its path.
pub fn build(layout: &Layout, legacy_kover: &Path) -> Result<std::path::PathBuf> {
    let kover = kover_benchmark(legacy_kover, REPORT_KOVER_REMEDY)?;
    let out_dir = layout.root.join("report");
    std::fs::create_dir_all(out_dir.join("datasets"))?;
    // `prune` removes every file this build did not publish, including one a
    // concurrent build has just published, so builds take turns: each holds
    // an exclusive lock on the report directory until it returns.
    let lock =
        std::fs::File::open(&out_dir).with_context(|| format!("open {}", out_dir.display()))?;
    lock.lock()
        .with_context(|| format!("lock {}", out_dir.display()))?;
    let datasets: Vec<DatasetSpec> = serde_json::from_slice(
        &std::fs::read(layout.curated().join("datasets.json")).context("run `grm curate` first")?,
    )?;
    let funnel: Funnel =
        serde_json::from_slice(&std::fs::read(layout.curated().join("funnel.json"))?)?;
    let phenotypes: Vec<Phenotype> = bvbrc::read_jsonl(&layout.curated().join("phenotypes.jsonl"))?;
    let genome_records: Vec<GenomeRecord> = bvbrc::read_jsonl(&layout.genomes())?;
    let genomes: HashMap<&str, &GenomeRecord> = genome_records
        .iter()
        .map(|g| (g.genome_id.as_str(), g))
        .collect();
    let records: Vec<AmrRecord> = bvbrc::read_jsonl(&layout.amr())?;
    let amrfinder_version = std::fs::read_to_string(layout.catalog_dir().join("version.txt"))
        .unwrap_or_default()
        .trim()
        .to_string();

    let mut species_names: BTreeMap<&str, Vec<&DatasetSpec>> = BTreeMap::new();
    for d in &datasets {
        species_names.entry(&d.species).or_default().push(d);
    }
    let mut by_species_pheno: HashMap<&str, Vec<&Phenotype>> = HashMap::new();
    for p in &phenotypes {
        by_species_pheno.entry(&p.species).or_default().push(p);
    }

    let mut catalog_rows = Vec::new();
    let mut species_rows = Vec::new();
    let mut panel_rows = Vec::new();
    let mut unsupported_rows = Vec::new();
    let mut report_ks: BTreeSet<Option<u64>> = BTreeSet::new();
    // File names this build publishes under datasets/ and species/.
    let (mut dataset_files, mut species_files) = (BTreeSet::new(), BTreeSet::new());
    for (species, specs) in &species_names {
        let slug = crate::curate::slug(species);
        let phen = by_species_pheno
            .get(species)
            .map(Vec::as_slice)
            .unwrap_or(&[]);
        let structure = species_file(layout, &slug, species, specs, phen, &genomes)?;
        if structure.is_some() {
            species_files.insert(format!("{slug}.json"));
        }
        let mut reported = 0;
        for spec in specs {
            let src = layout.results_dir(&slug).join(format!("{}.json", spec.id));
            let bytes = match std::fs::read(&src) {
                Ok(bytes) => bytes,
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => continue,
                Err(e) => return Err(e).with_context(|| format!("read {}", src.display())),
            };
            let report: Value = serde_json::from_slice(&bytes)
                .with_context(|| format!("parse {}", src.display()))?;
            let file = format!("{}.json", spec.id);
            publish(&out_dir.join("datasets").join(&file), &bytes)?;
            dataset_files.insert(file);
            reported += 1;
            let splits = report["splits"].as_array().cloned().unwrap_or_default();
            let scheme = |name: &str| splits.iter().find(|s| s["scheme"] == name);
            let learners = learner_rows(&report["splits"]);
            let best = best_learner(&learners);
            let counts = &report["counts"];
            report_ks.insert(counts["k"].as_u64());
            let panels: Vec<(String, Lite)> = report["panels"]
                .as_array()
                .map(|a| {
                    a.iter()
                        .map(|p| {
                            (
                                p["id"].as_str().unwrap_or("").to_string(),
                                lite(&p["metrics"]),
                            )
                        })
                        .collect()
                })
                .unwrap_or_default();
            if let Some((b, _)) = best {
                let best_model_counts = learner_row(scheme("random"), &b.learner)
                    .and_then(|l| Confusion::deserialize(&l["metrics"]).ok());
                for p in report["panels"].as_array().into_iter().flatten() {
                    panel_rows.push(PanelSummaryRow {
                        panel: p["id"].as_str().unwrap_or("").into(),
                        name: p["name"].as_str().unwrap_or("").into(),
                        dataset: spec.id.clone(),
                        species: spec.species.clone(),
                        antibiotic: spec.antibiotic.clone(),
                        metrics: lite(&p["metrics"]),
                        best_model: b.random.unwrap(),
                        best_learner: b.learner.clone(),
                        counts: Confusion::deserialize(&p["metrics"]).ok(),
                        best_model_counts,
                        no_call: NoCall::deserialize(&p["no_call"]).ok(),
                    });
                }
            }
            for u in report["unsupported_assays"]
                .as_array()
                .into_iter()
                .flatten()
            {
                unsupported_rows.push(UnsupportedPanelRow {
                    panel: u["id"].as_str().unwrap_or("").into(),
                    name: u["name"].as_str().unwrap_or("").into(),
                    dataset: spec.id.clone(),
                    species: spec.species.clone(),
                    antibiotic: spec.antibiotic.clone(),
                    reason: u["reason"].as_str().unwrap_or("").into(),
                });
            }
            catalog_rows.push(CatalogRow {
                id: spec.id.clone(),
                species: spec.species.clone(),
                species_slug: slug.clone(),
                antibiotic: spec.antibiotic.clone(),
                drug_class: spec.drug_class.clone(),
                n: counts["samples"].as_u64().unwrap_or(0) as usize,
                n_resistant: counts["resistant"].as_u64().unwrap_or(0) as usize,
                labels_from_mic: counts["labels_from_mic"].as_u64().unwrap_or(0) as usize,
                labels_from_mic_excluded: counts["labels_from_mic_excluded"]
                    .as_u64()
                    .map(|n| n as usize),
                min_class: counts["min_class"]
                    .as_u64()
                    .map_or_else(crate::curate::unrecorded_min_class, |n| n as usize),
                evaluation: report.get("evaluation").cloned().unwrap_or_else(|| {
                    serde_json::json!({
                        "revision": 0,
                        "valid": false,
                        "reasons": ["This report predates the evaluation audit; its metrics are provisional and must be recomputed."]
                    })
                }),
                features: counts["features"].as_u64().unwrap_or(0) as usize,
                lineages: counts["lineages"].as_u64().unwrap_or(0) as usize,
                best_learner: best.map(|(b, _)| b.learner.clone()),
                best_random: best.map(|(_, ba)| ba),
                best_lineage: best
                    .and_then(|(b, _)| b.lineage)
                    .and_then(|m| m.balanced_accuracy),
                lineage_unavailable_reason: scheme("lineage")
                    .and_then(|s| s["unavailable_reason"].as_str())
                    .map(String::from),
                catalog: report["catalog"].get("metrics").map(lite),
                mic_essential_agreement: report["mic"]
                    .get("essential_agreement")
                    .and_then(Value::as_f64),
                scm_rules: rule_text(&report),
                panels,
                kover: kover
                    .get(&(spec.species.clone(), spec.antibiotic.clone()))
                    .cloned(),
                learners,
            });
        }
        let labeled: BTreeSet<&str> = phen
            .iter()
            .filter(|p| p.label.is_some())
            .map(|p| p.genome_id.as_str())
            .collect();
        species_rows.push(SpeciesRow {
            species: species.to_string(),
            slug,
            datasets_curated: specs.len(),
            datasets_reported: reported,
            genomes_labeled: labeled.len(),
            genomes_analyzed: structure.map_or(0, |s| s.0),
            lineages: structure.map(|s| s.1),
            lineage_threshold: structure.map(|s| s.2),
        });
    }
    species_rows.sort_by_key(|a| std::cmp::Reverse(a.genomes_labeled));

    let drug_names: Value = serde_json::from_slice(
        &std::fs::read(layout.curated().join("drug_names.json")).unwrap_or_default(),
    )
    .unwrap_or(Value::Null);
    let manifest = Manifest {
        schema: SCHEMA,
        generated_at: now_iso(),
        evaluation_revision: crate::analysis::EVALUATION_REVISION,
        k: match report_ks.len() {
            1 => report_ks.pop_first().flatten().map(|k| k as usize),
            _ => None,
        },
        amr_records: records.len(),
        amrfinder_version,
        name_vocabulary: drug_names["vocabulary"]
            .as_str()
            .map(|_| "SigDrug reviewed drug vocabulary".to_string()),
        funnel,
        species: species_rows,
        datasets: catalog_rows,
        panels: panel_rows,
        unsupported_panels: unsupported_rows,
    };

    // Corpus-level distributions.
    let by_species: Vec<Value> = species_names
        .keys()
        .map(|sp| {
            let phen = by_species_pheno.get(sp).map(Vec::as_slice).unwrap_or(&[]);
            let r = phen.iter().filter(|p| p.label == Some(Call::R)).count();
            let s = phen.iter().filter(|p| p.label == Some(Call::S)).count();
            let g: BTreeSet<&str> = phen.iter().map(|p| p.genome_id.as_str()).collect();
            serde_json::json!({"species": sp, "pairs_resistant": r, "pairs_susceptible": s, "genomes": g.len()})
        })
        .collect();
    let labeled_genomes: BTreeSet<&str> = phenotypes
        .iter()
        .filter(|p| p.label.is_some())
        .map(|p| p.genome_id.as_str())
        .collect();
    let years = top_counts(
        labeled_genomes
            .iter()
            .filter_map(|g| genomes.get(g).and_then(|m| m.collection_year.as_deref())),
        200,
    );
    let countries = top_counts(
        labeled_genomes
            .iter()
            .filter_map(|g| genomes.get(g).and_then(|m| m.isolation_country.as_deref())),
        80,
    );
    let methods = top_counts(
        records
            .iter()
            .filter_map(|r| r.laboratory_typing_method.as_deref()),
        20,
    );
    let standards = top_counts(
        records.iter().filter_map(|r| r.testing_standard.as_deref()),
        20,
    );
    let breakpoints: Vec<Value> = datasets
        .iter()
        .filter_map(|d| {
            let b: Breakpoint = d.breakpoint?;
            Some(serde_json::json!({
                "species": d.species, "antibiotic": d.antibiotic, "r_above_log2": b.r_above,
                "r_above_mg_l": 2f64.powi(b.r_above), "agreement": b.agreement, "n": b.n,
                "labels_from_mic": d.n_labels_from_mic,
            }))
        })
        .collect();
    let corpus = serde_json::json!({
        "schema": SCHEMA,
        "by_species": by_species,
        "years": years,
        "countries": countries,
        "methods": methods,
        "standards": standards,
        "breakpoints": breakpoints,
        "drug_names": drug_names["names"],
    });
    publish(&out_dir.join("corpus.json"), &serde_json::to_vec(&corpus)?)?;
    // The manifest advertises datasets only after their files are published.
    publish(&out_dir.join("index.json"), &serde_json::to_vec(&manifest)?)?;
    // An earlier build's file for a dataset or species this build left out
    // copies a result or structure that is gone, such as the result `grm
    // learn` removes when a cohort falls below its class minimum.
    prune(&out_dir.join("datasets"), &dataset_files)?;
    prune(&out_dir.join("species"), &species_files)?;
    Ok(out_dir)
}

/// Removes each `*.json` file in `dir` that `keep` does not name. Pending
/// writes (`*.tmp`) and other files stay.
fn prune(dir: &Path, keep: &BTreeSet<String>) -> Result<()> {
    let entries = match std::fs::read_dir(dir) {
        Ok(entries) => entries,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(()),
        Err(e) => return Err(e).with_context(|| format!("read {}", dir.display())),
    };
    for entry in entries {
        let entry = entry?;
        let name = entry.file_name();
        if name
            .to_str()
            .is_some_and(|n| n.ends_with(".json") && !keep.contains(n))
        {
            std::fs::remove_file(entry.path())
                .with_context(|| format!("remove {}", entry.path().display()))?;
        }
    }
    Ok(())
}

pub(crate) fn now_iso() -> String {
    iso_from_secs(
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_secs()) as i64,
    )
}

fn iso_from_secs(secs: i64) -> String {
    let (days, rem) = (secs.div_euclid(86_400), secs.rem_euclid(86_400));
    // Civil-from-days (Howard Hinnant), UTC.
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + (month <= 2) as i64;
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}Z",
        rem / 3600,
        rem % 3600 / 60,
        rem % 60
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(unix)]
    #[test]
    fn publishing_replaces_the_file_without_truncating_existing_readers() {
        use std::io::Read;

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("index.json");
        publish(&path, br#"{"datasets":[]}"#).unwrap();
        let mut reader = std::fs::File::open(&path).unwrap();
        publish(&path, br#"{"datasets":["new"]}"#).unwrap();

        let mut previous = String::new();
        reader.read_to_string(&mut previous).unwrap();
        assert_eq!(previous, r#"{"datasets":[]}"#);
        assert_eq!(
            std::fs::read_to_string(&path).unwrap(),
            r#"{"datasets":["new"]}"#
        );
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 1);
    }

    #[test]
    fn iso_dates_match_known_days() {
        // 2026-10-23 is day 20749 since the epoch; 2024-02-29 is day 19782.
        assert_eq!(
            iso_from_secs(20_749 * 86_400 + 3_661),
            "2026-10-23T01:01:01Z"
        );
        assert_eq!(iso_from_secs(19_782 * 86_400), "2024-02-29T00:00:00Z");
        assert_eq!(iso_from_secs(0), "1970-01-01T00:00:00Z");
    }

    fn built_index(reports: &[(&str, Value)]) -> Value {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        write_reports(&layout, reports);
        let out = build(&layout, &dir.path().join("kover.json")).unwrap();
        serde_json::from_slice(&std::fs::read(out.join("index.json")).unwrap()).unwrap()
    }

    /// Curated E. coli datasets for `reports`' antibiotics, each with its
    /// report saved, empty record files, and a Kover comparison with no rows
    /// at `<root>/kover.json`.
    fn write_reports(layout: &Layout, reports: &[(&str, Value)]) {
        std::fs::create_dir_all(&layout.root).unwrap();
        std::fs::write(layout.root.join("kover.json"), b"[]").unwrap();
        let species = "Escherichia coli";
        let slug = crate::curate::slug(species);
        let specs: Vec<DatasetSpec> = reports
            .iter()
            .map(|(antibiotic, _)| DatasetSpec {
                id: format!("{slug}__{antibiotic}"),
                species: species.into(),
                species_slug: slug.clone(),
                antibiotic: antibiotic.to_string(),
                drug_class: "carbapenem".into(),
                n_resistant: 3,
                n_susceptible: 4,
                n_intermediate: 0,
                n_with_mic: 0,
                n_labels_from_mic: 0,
                breakpoint: None,
                min_class: 3,
                selection: None,
            })
            .collect();
        std::fs::create_dir_all(layout.curated()).unwrap();
        std::fs::create_dir_all(layout.results_dir(&slug)).unwrap();
        std::fs::create_dir_all(layout.root.join("bvbrc")).unwrap();
        std::fs::write(
            layout.curated().join("datasets.json"),
            serde_json::to_vec(&specs).unwrap(),
        )
        .unwrap();
        std::fs::write(
            layout.curated().join("funnel.json"),
            serde_json::to_vec(&Funnel::default()).unwrap(),
        )
        .unwrap();
        for path in [
            layout.curated().join("phenotypes.jsonl"),
            layout.genomes(),
            layout.amr(),
        ] {
            std::fs::write(path, b"").unwrap();
        }
        for (spec, (_, report)) in specs.iter().zip(reports) {
            std::fs::write(
                layout.results_dir(&slug).join(format!("{}.json", spec.id)),
                serde_json::to_vec(report).unwrap(),
            )
            .unwrap();
        }
    }

    #[test]
    fn index_rows_carry_the_kover_row_for_their_species_and_antibiotic() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let report = || serde_json::json!({"counts": {"samples": 7, "resistant": 3}});
        write_reports(&layout, &[("meropenem", report()), ("imipenem", report())]);
        let kover = dir.path().join("kover.json");
        std::fs::write(
            &kover,
            br#"[{"species":"Escherichia coli","antibiotic":"Meropenem","sensitivity":0.8,
                 "specificity":0.6,"n_rules":2.0,"ds_n_examples":40.0}]"#,
        )
        .unwrap();
        let out = build(&layout, &kover).unwrap();
        let index: Value =
            serde_json::from_slice(&std::fs::read(out.join("index.json")).unwrap()).unwrap();
        let row = |id: &str| {
            index["datasets"]
                .as_array()
                .unwrap()
                .iter()
                .find(|d| d["id"] == id)
                .unwrap()
                .clone()
        };
        let meropenem = row("escherichia_coli__meropenem");
        assert_eq!(meropenem["kover"]["sensitivity"], 0.8);
        assert_eq!(meropenem["kover"]["balanced_accuracy"], 0.7);
        assert_eq!(meropenem["kover"]["n_examples"], 40.0);
        assert!(row("escherichia_coli__imipenem")["kover"].is_null());
    }

    #[test]
    fn a_missing_or_malformed_kover_comparison_refuses_before_writing() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let report = serde_json::json!({"counts": {"samples": 7, "resistant": 3}});
        write_reports(&layout, &[("meropenem", report)]);
        let kover = dir.path().join("kover.json");
        std::fs::remove_file(&kover).unwrap();
        // Run from a folder without reference/, the default path is missing:
        // the build names the path it expected and leaves no report behind.
        let err = format!("{:#}", build(&layout, &kover).unwrap_err());
        assert!(err.contains(kover.to_str().unwrap()), "{err}");
        assert!(err.contains("reference/kover-summary.json"), "{err}");
        assert!(err.contains("--kover"), "{err}");
        assert!(!layout.root.join("report").exists());

        std::fs::write(&kover, b"{not json").unwrap();
        let err = format!("{:#}", build(&layout, &kover).unwrap_err());
        assert!(err.contains("parse the Kover comparison"), "{err}");
        assert!(err.contains(kover.to_str().unwrap()), "{err}");
        assert!(!layout.root.join("report").exists());

        // A relative path is shown as given and resolved against the
        // working folder.
        let relative = Path::new("no-such-folder-grm-test/kover-summary.json");
        let err = format!("{:#}", build(&layout, relative).unwrap_err());
        let abs = std::path::absolute(relative).unwrap();
        assert!(
            err.contains(&format!("{} ({})", relative.display(), abs.display())),
            "{err}"
        );
    }

    #[test]
    fn a_rebuild_removes_files_for_datasets_and_species_it_no_longer_publishes() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let report = || serde_json::json!({"counts": {"samples": 7, "resistant": 3}});
        write_reports(&layout, &[("meropenem", report()), ("imipenem", report())]);
        let kover = dir.path().join("kover.json");
        let out = build(&layout, &kover).unwrap();
        let copy =
            |antibiotic: &str| out.join(format!("datasets/escherichia_coli__{antibiotic}.json"));
        assert!(copy("meropenem").exists() && copy("imipenem").exists());
        // `grm learn` removed imipenem's result, as it does when a cohort
        // falls below its class minimum. An earlier build published a
        // species that is no longer curated, and a concurrent build's write
        // is pending.
        let results = layout.results_dir("escherichia_coli");
        std::fs::remove_file(results.join("escherichia_coli__imipenem.json")).unwrap();
        let old_species = out.join("species/vibrio_cholerae.json");
        crate::layout::ensure_parent(&old_species).unwrap();
        std::fs::write(&old_species, "{}").unwrap();
        let pending = out.join("datasets/escherichia_coli__imipenem.json.4242.tmp");
        std::fs::write(&pending, "{}").unwrap();
        build(&layout, &kover).unwrap();
        assert!(copy("meropenem").exists());
        assert!(!copy("imipenem").exists());
        assert!(!old_species.exists());
        assert!(pending.exists());
        let index: Value =
            serde_json::from_slice(&std::fs::read(out.join("index.json")).unwrap()).unwrap();
        assert_eq!(index["datasets"].as_array().unwrap().len(), 1);
    }

    #[test]
    fn unreadable_result_preserves_the_published_file_and_manifest() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let report = serde_json::json!({"counts": {"samples": 7, "resistant": 3}});
        write_reports(&layout, &[("meropenem", report)]);
        let kover = dir.path().join("kover.json");
        let out = build(&layout, &kover).unwrap();
        let file = "escherichia_coli__meropenem.json";
        let published = out.join("datasets").join(file);
        let previous = std::fs::read(&published).unwrap();
        let manifest = std::fs::read(out.join("index.json")).unwrap();
        let index: Value = serde_json::from_slice(&manifest).unwrap();
        assert_eq!(index["datasets"][0]["id"], "escherichia_coli__meropenem");

        let source = layout.results_dir("escherichia_coli").join(file);
        std::fs::remove_file(&source).unwrap();
        // A directory fails even when the test runs as root.
        std::fs::create_dir(&source).unwrap();
        let kind = std::fs::read(&source).unwrap_err().kind();
        assert_ne!(kind, std::io::ErrorKind::NotFound);
        let err = build(&layout, &kover).unwrap_err();
        assert_eq!(err.downcast_ref::<std::io::Error>().unwrap().kind(), kind);
        assert!(
            err.to_string().contains(source.to_str().unwrap()),
            "{err:#}"
        );
        assert!(source.is_dir());
        assert_eq!(std::fs::read(&published).unwrap(), previous);
        assert_eq!(std::fs::read(out.join("index.json")).unwrap(), manifest);
    }

    #[test]
    fn unreadable_structure_preserves_the_published_file_and_manifest() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let report = serde_json::json!({"counts": {"samples": 7, "resistant": 3}});
        write_reports(&layout, &[("meropenem", report)]);
        let species = layout.species_dir("escherichia_coli");
        let source = species.join("structure.json");
        std::fs::write(
            &source,
            serde_json::to_vec(&serde_json::json!({
                "genomes": [], "coords": [], "nearest": [], "merges": [],
                "leaf_order": [], "sketch_size": 1000,
                "lineages": {"threshold": 0.01, "labels": [], "n_clusters": 0, "largest_share": 0.0}
            }))
            .unwrap(),
        )
        .unwrap();
        std::fs::write(
            species.join("counted.json"),
            br#"{"k":31,"genomes":[],"missing_fasta":[],"failed":[]}"#,
        )
        .unwrap();
        let kover = dir.path().join("kover.json");
        let out = build(&layout, &kover).unwrap();
        let published = out.join("species/escherichia_coli.json");
        let previous = std::fs::read(&published).unwrap();
        let manifest = std::fs::read(out.join("index.json")).unwrap();

        std::fs::remove_file(&source).unwrap();
        std::fs::create_dir(&source).unwrap();
        let kind = std::fs::read(&source).unwrap_err().kind();
        assert_ne!(kind, std::io::ErrorKind::NotFound);
        let err = build(&layout, &kover).unwrap_err();
        assert_eq!(err.downcast_ref::<std::io::Error>().unwrap().kind(), kind);
        assert!(
            err.to_string().contains(source.to_str().unwrap()),
            "{err:#}"
        );
        assert!(source.is_dir());
        assert_eq!(std::fs::read(&published).unwrap(), previous);
        assert_eq!(std::fs::read(out.join("index.json")).unwrap(), manifest);
    }

    #[test]
    fn a_build_waits_until_the_running_build_finishes() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let report = serde_json::json!({"counts": {"samples": 7, "resistant": 3}});
        write_reports(&layout, &[("meropenem", report)]);
        let out = layout.root.join("report");
        std::fs::create_dir_all(&out).unwrap();
        // The lock a running build holds on the report directory.
        let running = std::fs::File::open(&out).unwrap();
        running.lock().unwrap();
        let (done, finished) = std::sync::mpsc::channel();
        let kover = dir.path().join("kover.json");
        let waiting = std::thread::spawn({
            let layout = layout.clone();
            move || done.send(build(&layout, &kover).map(|_| ())).unwrap()
        });
        let wait = std::time::Duration::from_millis(500);
        assert!(
            finished.recv_timeout(wait).is_err(),
            "built during another build"
        );
        assert!(!out.join("index.json").exists());
        drop(running);
        finished.recv().unwrap().unwrap();
        waiting.join().unwrap();
        assert!(out.join("index.json").exists());
    }

    #[test]
    fn manifest_rows_carry_report_panel_counts_and_lineage_reason() {
        use crate::analysis::PanelResult;
        use crate::learn::binary_metrics;
        use crate::learn::eval::{LearnerResult, SplitResult};

        let truth = [true, true, true, false, false, false, false];
        let learner = |name: &str, scores: &[f64]| LearnerResult {
            learner: name.into(),
            metrics: Some(binary_metrics(&truth, scores, 0.5, false)),
            called_metrics: None,
            fold_balanced_accuracy: vec![],
            roc: vec![],
            scores: scores.to_vec(),
            conformal: vec![],
            fold_complexity: vec![],
            seconds: 0.0,
            no_call: None,
            lasso_fits: None,
        };
        let reason = "fewer than twice as many lineages as requested folds";
        let random = SplitResult {
            scheme: "random".into(),
            folds: 5,
            n_groups: 4,
            test_fold: vec![],
            learners: vec![
                learner("scm", &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
                learner("rf", &[0.9, 0.8, 0.7, 0.6, 0.1, 0.2, 0.3]),
            ],
            unavailable_reason: None,
        };
        let lineage = SplitResult {
            scheme: "lineage".into(),
            folds: 5,
            n_groups: 4,
            test_fold: vec![],
            learners: vec![],
            unavailable_reason: Some(reason.into()),
        };
        let panel = PanelResult {
            id: "carbapenemase".into(),
            name: "Carbapenemase PCR".into(),
            mirrors: String::new(),
            targets: vec![],
            metrics: binary_metrics(&truth, &[1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0], 0.5, false),
            detected: vec![],
            no_call: None,
        };
        let report = serde_json::json!({
            "counts": {"samples": 7, "resistant": 3},
            "splits": [random, lineage],
            "panels": [panel],
        });

        let index = built_index(&[("meropenem", report)]);
        let (row, panel) = (&index["datasets"][0], &index["panels"][0]);
        assert_eq!(row["best_learner"], "rf");
        assert_eq!(row["best_lineage"], Value::Null);
        assert_eq!(row["lineage_unavailable_reason"], reason);
        assert_eq!(
            panel["counts"],
            serde_json::json!({"n": 7, "tp": 2, "fp": 1, "tn": 3, "fn": 1})
        );
        assert_eq!(
            panel["best_model_counts"],
            serde_json::json!({"n": 7, "tp": 3, "fp": 1, "tn": 3, "fn": 0})
        );
        let count = |k: &str| panel["counts"][k].as_u64().unwrap();
        assert_eq!(count("n"), row["n"].as_u64().unwrap());
        assert_eq!(
            count("tp") + count("fn"),
            row["n_resistant"].as_u64().unwrap()
        );
        assert_eq!(
            panel["metrics"]["sensitivity"].as_f64().unwrap(),
            count("tp") as f64 / (count("tp") + count("fn")) as f64
        );
    }

    #[test]
    fn manifest_rows_omit_counts_and_reasons_their_reports_lack() {
        let rates = serde_json::json!({
            "balanced_accuracy": 0.75, "sensitivity": 0.5, "specificity": 1.0,
            "vme": 0.5, "me": 0.0, "auroc": null
        });
        let report = serde_json::json!({
            "counts": {"samples": 7, "resistant": 3},
            "splits": [
                {"scheme": "random", "learners": [{"learner": "rf", "metrics": rates}]},
                {"scheme": "lineage", "learners": [{"learner": "rf", "metrics": rates}]},
            ],
            "panels": [{"id": "carbapenemase", "name": "Carbapenemase PCR", "metrics": rates}],
        });

        let index = built_index(&[("meropenem", report)]);
        let (row, panel) = (&index["datasets"][0], &index["panels"][0]);
        assert_eq!(row["best_lineage"], 0.75);
        assert_eq!(panel["metrics"]["sensitivity"], 0.5);
        for (object, key) in [
            (row, "lineage_unavailable_reason"),
            (row, "labels_from_mic_excluded"),
            (panel, "counts"),
            (panel, "best_model_counts"),
            (panel, "no_call"),
        ] {
            assert!(object.get(key).is_none(), "{key}");
        }
        assert_eq!(index["unsupported_panels"], serde_json::json!([]));
    }

    #[test]
    fn manifest_rows_carry_the_class_minimum_and_excluded_mic_labels() {
        let report = serde_json::json!({
            "counts": {
                "samples": 7, "resistant": 3, "labels_from_mic": 0,
                "labels_from_mic_excluded": 4, "min_class": 3
            },
        });
        let index = built_index(&[("meropenem", report)]);
        let row = &index["datasets"][0];
        assert_eq!(row["labels_from_mic"], 0);
        assert_eq!(row["labels_from_mic_excluded"], 4);
        assert_eq!(row["min_class"], 3);
    }

    #[test]
    fn manifest_rows_read_a_report_without_min_class_as_the_old_learn_gate_of_25() {
        let report = serde_json::json!({"counts": {"samples": 7, "resistant": 3}});
        let index = built_index(&[("meropenem", report)]);
        let row = &index["datasets"][0];
        assert_eq!(row["min_class"], 25);
        assert!(row.get("labels_from_mic_excluded").is_none());
    }

    #[test]
    fn manifest_rows_carry_probe_no_calls_and_unsupported_assays() {
        let rates = serde_json::json!({
            "n": 5, "tp": 2, "fp": 1, "tn": 1, "fn": 1,
            "balanced_accuracy": 0.58, "sensitivity": 0.67, "specificity": 0.5,
            "vme": 0.33, "me": 0.5, "auroc": null
        });
        let reason = "synonymous base changes read as mutations";
        let report = serde_json::json!({
            "counts": {"samples": 7, "resistant": 5},
            "splits": [{"scheme": "random", "learners": [{"learner": "rf", "metrics": rates}]}],
            "panels": [{
                "id": "gc_gyra", "name": "gyrA codons 91-95 probes", "metrics": rates,
                "no_call": {"resistant": 2, "susceptible": 0}
            }],
            "unsupported_assays": [{
                "id": "ecoli_qrdr", "name": "gyrA and parC QRDR probes",
                "mirrors": "QRDR", "reason": reason
            }],
        });

        let index = built_index(&[("ciprofloxacin", report)]);
        assert_eq!(
            index["panels"][0]["no_call"],
            serde_json::json!({"resistant": 2, "susceptible": 0})
        );
        let unsupported = &index["unsupported_panels"][0];
        assert_eq!(unsupported["panel"], "ecoli_qrdr");
        assert_eq!(unsupported["name"], "gyrA and parC QRDR probes");
        assert_eq!(unsupported["antibiotic"], "ciprofloxacin");
        assert_eq!(unsupported["reason"], reason);
        assert_eq!(index["unsupported_panels"].as_array().unwrap().len(), 1);
    }

    #[test]
    fn manifest_k_is_the_one_k_every_report_records() {
        let at = |k: Option<u64>| {
            let mut counts = serde_json::json!({"samples": 7, "resistant": 3});
            if let Some(k) = k {
                counts["k"] = k.into();
            }
            serde_json::json!({ "counts": counts })
        };
        let k_of = |a, b| built_index(&[("meropenem", at(a)), ("imipenem", at(b))])["k"].clone();
        assert_eq!(k_of(Some(21), Some(21)), 21);
        assert_eq!(k_of(Some(21), Some(31)), Value::Null);
        assert_eq!(k_of(Some(21), None), Value::Null);
        assert_eq!(k_of(None, None), Value::Null);
    }

    #[test]
    fn the_best_learner_calls_every_random_split_sample() {
        // Lasso has the highest balanced accuracy, but only over the samples
        // it called; the random forest called all of them.
        let rates = |ba: f64| {
            serde_json::json!({
                "n": 7, "tp": 3, "fp": 0, "tn": 4, "fn": 0, "balanced_accuracy": ba,
                "sensitivity": 1.0, "specificity": 1.0, "vme": 0.0, "me": 0.0, "auroc": null
            })
        };
        let row = |learner: &str, ba: f64, uncalled: usize| {
            serde_json::json!({
                "learner": learner, "metrics": rates(ba),
                "no_call": {"resistant": uncalled, "susceptible": 0}
            })
        };
        let report = serde_json::json!({
            "counts": {"samples": 7, "resistant": 3},
            "splits": [
                {"scheme": "random", "learners": [row("rf", 0.8, 0), row("lasso", 0.9, 2)]},
                {"scheme": "lineage", "learners": [row("rf", 0.7, 0), row("lasso", 0.6, 3)]},
            ],
            "panels": [{"id": "carbapenemase", "name": "Carbapenemase PCR", "metrics": rates(0.5)}],
        });

        let index = built_index(&[("meropenem", report)]);
        let row = &index["datasets"][0];
        assert_eq!(row["best_learner"], "rf");
        assert_eq!(row["best_random"], 0.8);
        assert_eq!(index["panels"][0]["best_learner"], "rf");
        let lasso = row["learners"]
            .as_array()
            .unwrap()
            .iter()
            .find(|l| l["learner"] == "lasso")
            .unwrap();
        assert_eq!(
            lasso["lineage"]["no_call"],
            serde_json::json!({"resistant": 3, "susceptible": 0})
        );
    }

    #[test]
    fn the_best_learner_needs_full_coverage_on_every_available_split() {
        // The forest and boosting score above the lasso on random splits and
        // call every isolate there. With lineages held out the forest leaves
        // 2 isolates uncalled and boosting calls none, as every lineage fold
        // of E. coli amikacin kept no tree in revision 5.
        let rates = |ba: f64| {
            serde_json::json!({
                "n": 7, "tp": 3, "fp": 0, "tn": 4, "fn": 0, "balanced_accuracy": ba,
                "sensitivity": 1.0, "specificity": 1.0, "vme": 0.0, "me": 0.0, "auroc": null
            })
        };
        let row = |learner: &str, ba: f64, uncalled: usize| {
            serde_json::json!({
                "learner": learner, "metrics": rates(ba),
                "no_call": {"resistant": uncalled, "susceptible": 0}
            })
        };
        let uncalled = serde_json::json!({
            "learner": "gbt", "metrics": null, "no_call": {"resistant": 3, "susceptible": 4}
        });
        let report = serde_json::json!({
            "counts": {"samples": 7, "resistant": 3},
            "splits": [
                {"scheme": "random", "learners": [row("rf", 0.9, 0), row("gbt", 0.95, 0), row("lasso", 0.8, 0)]},
                {"scheme": "lineage", "learners": [row("rf", 0.7, 2), uncalled, row("lasso", 0.6, 0)]},
            ],
            "panels": [{"id": "carbapenemase", "name": "Carbapenemase PCR", "metrics": rates(0.5)}],
        });
        let index = built_index(&[("meropenem", report)]);
        let row = &index["datasets"][0];
        assert_eq!(row["best_learner"], "lasso");
        assert_eq!(row["best_random"], 0.8);
        assert_eq!(row["best_lineage"], 0.6);
        assert_eq!(index["panels"][0]["best_learner"], "lasso");
    }

    #[test]
    fn revision_14_rows_with_null_fits_compete_for_best_learner() {
        // Boosting had null-fit folds on both splits. In a revision-14 row
        // those isolates are scored at their fold's training prevalence and
        // the row keeps the older figure in `called_metrics`, so its rates
        // cover every isolate and it can be best. The same counts on a row
        // without `called_metrics` (revisions 7 to 13) cover only the called
        // isolates, and that row is passed over.
        let rates = |ba: f64| {
            serde_json::json!({
                "n": 7, "tp": 3, "fp": 0, "tn": 4, "fn": 0, "balanced_accuracy": ba,
                "sensitivity": 1.0, "specificity": 1.0, "vme": 0.0, "me": 0.0, "auroc": null
            })
        };
        let row = |learner: &str, ba: f64, in_null: usize, scored: bool| {
            let mut r = serde_json::json!({
                "learner": learner, "metrics": rates(ba),
                "no_call": {"resistant": in_null, "susceptible": 0}
            });
            if scored && in_null > 0 {
                r["called_metrics"] = rates(ba + 0.05);
            }
            r
        };
        let report = |scored: bool| {
            serde_json::json!({
                "counts": {"samples": 7, "resistant": 3},
                "splits": [
                    {"scheme": "random", "learners": [row("rf", 0.8, 0, scored), row("gbt", 0.85, 2, scored)]},
                    {"scheme": "lineage", "learners": [row("rf", 0.7, 0, scored), row("gbt", 0.6, 3, scored)]},
                ],
                "panels": [],
            })
        };
        let index = built_index(&[("meropenem", report(true))]);
        let row = &index["datasets"][0];
        assert_eq!(row["best_learner"], "gbt");
        assert_eq!(row["best_random"], 0.85);
        assert_eq!(row["best_lineage"], 0.6);
        let gbt = row["learners"]
            .as_array()
            .unwrap()
            .iter()
            .find(|l| l["learner"] == "gbt")
            .unwrap();
        assert_eq!(gbt["random"]["null_fits_scored"], true);
        assert_eq!(gbt["random"]["no_call"]["resistant"], 2);
        let rf = row["learners"]
            .as_array()
            .unwrap()
            .iter()
            .find(|l| l["learner"] == "rf")
            .unwrap();
        assert!(rf["random"].get("null_fits_scored").is_none());
        let index = built_index(&[("meropenem", report(false))]);
        assert_eq!(index["datasets"][0]["best_learner"], "rf");
        assert!(
            index["datasets"][0]["learners"]
                .as_array()
                .unwrap()
                .iter()
                .all(|l| l["random"].get("null_fits_scored").is_none())
        );
    }

    #[test]
    fn a_missing_lineage_learner_row_is_not_full_coverage() {
        let rates = |ba: f64| {
            serde_json::json!({
                "n": 7, "tp": 3, "fp": 0, "tn": 4, "fn": 0, "balanced_accuracy": ba,
                "sensitivity": 1.0, "specificity": 1.0, "vme": 0.0, "me": 0.0, "auroc": null
            })
        };
        let row = |learner: &str, ba: f64| {
            serde_json::json!({
                "learner": learner, "metrics": rates(ba),
                "no_call": {"resistant": 0, "susceptible": 0}
            })
        };
        let report = serde_json::json!({
            "counts": {"samples": 7, "resistant": 3},
            "splits": [
                {"scheme": "random", "learners": [row("rf", 0.9), row("lasso", 0.8)]},
                {"scheme": "lineage", "learners": [row("lasso", 0.6)]},
            ],
            "panels": [],
        });
        let index = built_index(&[("meropenem", report)]);
        let row = &index["datasets"][0];
        assert_eq!(row["best_learner"], "lasso");
        assert_eq!(row["best_random"], 0.8);
        assert_eq!(row["best_lineage"], 0.6);
    }

    #[test]
    fn a_learner_that_called_no_isolate_has_null_rates_in_the_index() {
        // Revision 5 scored E. coli amikacin's lineage boosting folds 0.5 and
        // reported 0% VME; with no call there is no rate to report.
        let rates = |ba: f64| {
            serde_json::json!({
                "n": 7, "tp": 3, "fp": 0, "tn": 4, "fn": 0, "balanced_accuracy": ba,
                "sensitivity": 1.0, "specificity": 1.0, "vme": 0.0, "me": 0.0, "auroc": 0.9
            })
        };
        let all = serde_json::json!({"resistant": 0, "susceptible": 0});
        let report = serde_json::json!({
            "counts": {"samples": 7, "resistant": 3},
            "splits": [
                {"scheme": "random", "learners": [
                    {"learner": "rf", "metrics": rates(0.8), "no_call": all},
                    {"learner": "gbt", "metrics": rates(0.6), "no_call": all},
                ]},
                {"scheme": "lineage", "learners": [
                    {"learner": "rf", "metrics": rates(0.7), "no_call": all},
                    {"learner": "gbt", "metrics": null, "no_call": {"resistant": 3, "susceptible": 4}},
                ]},
            ],
            "panels": [],
        });
        let index = built_index(&[("amikacin", report)]);
        let gbt = index["datasets"][0]["learners"]
            .as_array()
            .unwrap()
            .iter()
            .find(|l| l["learner"] == "gbt")
            .unwrap();
        assert_eq!(
            gbt["lineage"],
            serde_json::json!({
                "balanced_accuracy": null, "sensitivity": null, "specificity": null,
                "vme": null, "me": null, "auroc": null,
                "no_call": {"resistant": 3, "susceptible": 4}
            })
        );
        assert_eq!(gbt["random"]["vme"], 0.0);
        assert_eq!(index["datasets"][0]["best_learner"], "rf");
    }

    #[test]
    fn the_index_is_schema_2_and_reads_schema_1_reports_beside_schema_2_ones() {
        let rates = |ba: f64| {
            serde_json::json!({
                "n": 7, "tp": 3, "fp": 0, "tn": 4, "fn": 0, "balanced_accuracy": ba,
                "sensitivity": 1.0, "specificity": 1.0, "vme": 0.0, "me": 0.0, "auroc": null
            })
        };
        // A revision-6 report: schema 1, every learner scored every isolate.
        let old = serde_json::json!({
            "schema": 1,
            "counts": {"samples": 7, "resistant": 3},
            "splits": [{"scheme": "random", "learners": [
                {"learner": "rf", "metrics": rates(0.8)},
                {"learner": "lasso", "metrics": rates(0.9)},
            ]}],
            "panels": [],
        });
        // Schema 2: the lasso leaves one isolate uncalled.
        let new = serde_json::json!({
            "schema": 2,
            "counts": {"samples": 7, "resistant": 3},
            "splits": [{"scheme": "random", "learners": [
                {"learner": "rf", "metrics": rates(0.8), "no_call": {"resistant": 0, "susceptible": 0}},
                {"learner": "lasso", "metrics": rates(0.9), "no_call": {"resistant": 1, "susceptible": 0}},
            ]}],
            "panels": [],
        });
        let index = built_index(&[("imipenem", old), ("meropenem", new)]);
        assert_eq!(index["schema"], 2);
        let (old, new) = (&index["datasets"][0], &index["datasets"][1]);
        assert_eq!(old["antibiotic"], "imipenem");
        assert_eq!(old["best_learner"], "lasso");
        assert_eq!(new["antibiotic"], "meropenem");
        assert_eq!(new["best_learner"], "rf");
        let lasso = old["learners"]
            .as_array()
            .unwrap()
            .iter()
            .find(|l| l["learner"] == "lasso")
            .unwrap();
        assert_eq!(lasso["random"]["balanced_accuracy"], 0.9);
        assert!(lasso["random"].get("no_call").is_none());
    }

    #[test]
    fn no_full_coverage_learner_means_no_headline_winner() {
        let rates = |ba: f64| {
            serde_json::json!({
                "n": 7, "tp": 3, "fp": 0, "tn": 4, "fn": 0, "balanced_accuracy": ba,
                "sensitivity": 1.0, "specificity": 1.0, "vme": 0.0, "me": 0.0, "auroc": null
            })
        };
        let report = serde_json::json!({
            "counts": {"samples": 7, "resistant": 3},
            "splits": [
                {"scheme": "random", "learners": [
                    {"learner": "rf", "metrics": rates(0.8), "no_call": {"resistant": 1, "susceptible": 0}}
                ]},
                {"scheme": "lineage", "learners": [
                    {"learner": "rf", "metrics": rates(0.7), "no_call": {"resistant": 1, "susceptible": 0}}
                ]}
            ],
            "panels": [{"id": "carbapenemase", "name": "Carbapenemase PCR", "metrics": rates(0.5)}],
        });
        let index = built_index(&[("meropenem", report)]);
        assert_eq!(
            index["datasets"][0]["best_learner"],
            serde_json::Value::Null
        );
        assert_eq!(index["datasets"][0]["best_random"], serde_json::Value::Null);
        assert_eq!(
            index["datasets"][0]["best_lineage"],
            serde_json::Value::Null
        );
        assert_eq!(index["panels"].as_array().unwrap().len(), 0);
    }

    #[test]
    fn manifest_states_the_current_revision_beside_each_report_audit() {
        use crate::analysis::EVALUATION_REVISION;

        let audited = |revision: u32| {
            serde_json::json!({
                "counts": {"samples": 7, "resistant": 3},
                "evaluation": {"revision": revision, "valid": true, "reasons": [], "scope": "internal"},
            })
        };
        let index = built_index(&[
            ("meropenem", audited(EVALUATION_REVISION)),
            ("imipenem", audited(EVALUATION_REVISION - 1)),
            (
                "ertapenem",
                serde_json::json!({"counts": {"samples": 7, "resistant": 3}}),
            ),
        ]);
        assert_eq!(index["evaluation_revision"], EVALUATION_REVISION);
        let audits: Vec<(u64, bool)> = index["datasets"]
            .as_array()
            .unwrap()
            .iter()
            .map(|d| {
                let a = &d["evaluation"];
                (
                    a["revision"].as_u64().unwrap(),
                    a["valid"].as_bool().unwrap(),
                )
            })
            .collect();
        let current = EVALUATION_REVISION as u64;
        assert_eq!(
            audits,
            vec![(current, true), (current - 1, true), (0, false)]
        );
    }

    #[test]
    fn a_rule_term_names_every_catalog_family_its_feature_matches() {
        let report = serde_json::json!({
            "features": [
                {"feature": 1, "catalog_families": ["blaKPC"], "locus": {"gene": "kpc"}},
                {"feature": 2, "catalog_families": ["blaNDM", "ble", "blaAFM"], "locus": null},
                {"feature": 3, "catalog_families": [], "locus": {"gene": null, "product": "porin"}},
                {"feature": 4, "catalog_families": [], "locus": null},
            ],
            "models": {"scm": {"rules": [
                {"feature": 1, "presence": true},
                {"feature": 2, "presence": true},
                {"feature": 2, "presence": false},
                {"feature": 3, "presence": false},
                {"feature": 4, "presence": true},
                {"feature": 5, "presence": true},
            ]}},
        });
        assert_eq!(
            rule_text(&report),
            [
                "presence of blaKPC",
                "presence of k-mers matching 3 catalog families (blaAFM, blaNDM, ble)",
                "absence of k-mers matching 3 catalog families (blaAFM, blaNDM, ble)",
                "absence of porin",
                "presence of feature 4",
                "presence of feature 5",
            ]
        );
    }
}
