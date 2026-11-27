//! Reusable cross-database deduplication engine (`grm dedup`).
//!
//! Screens an external or newly acquired genome batch against one or more
//! versioned reference indices (`grm-dedup-index-v1`) and performs within-batch
//! connected-component deduplication and longitudinal patient grouping.
//!
//! ## Identity and overlap layers
//!
//! Every layer is recorded independently per query so overlapping evidence is
//! never collapsed:
//!
//! - Accession and alias equivalence (through a pluggable `--synonym-map` for
//!   INSDC/ENA/NCBI links such as `SAMEA`<->`ERS`, `SAMN`<->`SRS`,
//!   `PRJEB`<->`ERP`, and run-to-sample links): `D1` assembly (`GCA_`/`GCF_`
//!   numeric core ignoring prefix and version), `D2` BioSample/SRA, and
//!   same-species normalized `STRAIN`.
//! - Normalized assembly content identity (`EXACT_ASSEMBLY`: SHA-256 over
//!   strand-canonical uppercase ACGTN contigs sorted lexicographically).
//! - Canonical k-mer-set identity (`EXACT_KSET`: `sha256_k31` when both sides
//!   carry it, else `fnv64_k31` matched against `counted.json` digests with
//!   k-mer count verification).
//! - Genomic Mash distance (`SketchIndex` inverted MinHash candidate filter
//!   verified with exact [`crate::cluster::mash_distance`]): near-duplicates
//!   (`NEAR_DUP`, `<= near_dup_threshold`) separated from close relatives
//!   (`CLOSE_RELATIVE`, in `(near_dup_threshold, close_relative_threshold]`),
//!   using species-specific thresholds (`--thresholds`).
//! - Study-overlap tags (`SHARED_BIOPROJECT`, `SHARED_PMID`) kept separate
//!   from isolate identity.
//! - Longitudinal patient grouping (`subject_id`) kept separate from
//!   within-batch duplicate collapse so serial patient isolates stay intact
//!   across splits.
//!
//! ## Epistemic limits
//!
//! Every screen report records `claim_scope: "zero_detected_overlap_under_stated_checks"`
//! and `sealed_set_status: "unqueried_unknown"`. Unopened sealed records are
//! never queried by default, and accession gaps in either reference or query
//! cohorts make metadata overlap counts lower bounds.

use crate::cluster;
use crate::curate;
use crate::dna::{MAX_K, MIN_K};
use crate::kset::{self, FastaFingerprint};
use crate::layout::ensure_parent;
use crate::research::sha256;
use crate::species::{self, CountSummary};
use crate::species_check::GenomeSketch;
use anyhow::{Context, Result, bail, ensure};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

/// Version tag for reference index manifests written by `grm dedup index`.
pub const INDEX_SCHEMA: &str = "grm-dedup-index-v1";
/// Version tag for screen summary manifests written by `grm dedup screen`.
pub const SCREEN_SCHEMA: &str = "grm-dedup-screen-v1";
/// Version tag for species-specific threshold configuration files.
pub const THRESHOLDS_SCHEMA: &str = "grm-dedup-thresholds-v1";
/// Explicit epistemic scope recorded in every screen report.
pub const CLAIM_SCOPE: &str = "zero_detected_overlap_under_stated_checks";
/// Explicit sealed-set status recorded in every index and screen report.
pub const SEALED_SET_STATUS: &str = "unqueried_unknown";

pub const DEFAULT_K_DIGEST: usize = 31;
pub const DEFAULT_K_SKETCH: usize = 21;
pub const DEFAULT_SKETCH_SIZE: usize = species::SKETCH;
pub const DEFAULT_INDEX_HASHES: usize = 64;
pub const DEFAULT_NEAR_DUP_THRESHOLD: f64 = 0.0005;
pub const DEFAULT_CLOSE_RELATIVE_THRESHOLD: f64 = 0.001;

// ---------------------------------------------------------------- ordering & normalization

#[derive(PartialEq, Eq, PartialOrd, Ord)]
enum GenomeIdKey<'a> {
    TaxonSuffix(u64, u64, &'a str),
    Lex(&'a str),
}

fn genome_id_key(id: &str) -> GenomeIdKey<'_> {
    if let Some((left, right)) = id.split_once('.')
        && !left.is_empty()
        && !right.is_empty()
        && left.bytes().all(|b| b.is_ascii_digit())
        && right.bytes().all(|b| b.is_ascii_digit())
        && let (Ok(taxon), Ok(suffix)) = (left.parse::<u64>(), right.parse::<u64>())
    {
        return GenomeIdKey::TaxonSuffix(taxon, suffix, id);
    }
    GenomeIdKey::Lex(id)
}

/// Deterministic total order over genome ids: R0 `(taxon, suffix)` integer-pair
/// order for BV-BRC ids, lexicographical order for all other ids, with
/// `taxon.suffix` ids ordered before non-numeric ids so transitivity holds on
/// mixed batches.
pub fn compare_genome_ids(a: &str, b: &str) -> Ordering {
    genome_id_key(a).cmp(&genome_id_key(b))
}

fn is_placeholder(upper: &str) -> bool {
    matches!(
        upper,
        "" | "-"
            | "."
            | "NA"
            | "N/A"
            | "NONE"
            | "NULL"
            | "MISSING"
            | "UNKNOWN"
            | "NOT APPLICABLE"
            | "NOT COLLECTED"
            | "NOT PROVIDED"
            | "UNAVAILABLE"
    )
}

fn split_multi_tokens(raw: &str) -> Vec<String> {
    raw.split(|c: char| matches!(c, ',' | ';' | '|') || c.is_whitespace())
        .map(|tok| tok.trim().to_ascii_uppercase())
        .filter(|tok| !is_placeholder(tok))
        .collect()
}

/// Normalizes a single assembly accession under R0: uppercase, trim, and for
/// `GCA_`/`GCF_` extract the numeric core ignoring prefix and `.version` so
/// GenBank and RefSeq copies share one key (`ASM:<digits>`).
pub fn normalize_assembly_token(raw: &str) -> Option<String> {
    let up = raw.trim().to_ascii_uppercase();
    if is_placeholder(&up) {
        return None;
    }
    if let Some(rest) = up.strip_prefix("GCA_").or_else(|| up.strip_prefix("GCF_")) {
        let core = rest.split('.').next().unwrap_or(rest);
        if !core.is_empty() && core.bytes().all(|b| b.is_ascii_digit()) {
            return Some(format!("ASM:{core}"));
        }
    }
    Some(up)
}

/// Splits and normalizes assembly accessions into a sorted deduplicated list.
pub fn normalize_assembly_keys(raw: &[String]) -> Vec<String> {
    let mut set = BTreeSet::new();
    for item in raw {
        for tok in split_multi_tokens(item) {
            if let Some(norm) = normalize_assembly_token(&tok) {
                set.insert(norm);
            }
        }
    }
    set.into_iter().collect()
}

/// Normalizes BioSample accessions (trim, uppercase, strip optional
/// `BIOSAMPLE:` prefix).
pub fn normalize_biosample_keys(raw: &[String]) -> Vec<String> {
    let mut set = BTreeSet::new();
    for item in raw {
        for tok in split_multi_tokens(item) {
            let clean = tok.strip_prefix("BIOSAMPLE:").unwrap_or(&tok);
            if !is_placeholder(clean) {
                set.insert(clean.to_string());
            }
        }
    }
    set.into_iter().collect()
}

/// Splits and normalizes SRA/ENA/DDBJ run/sample accessions under R0
/// (separators `,`, `;`, `|`, whitespace; each token uppercase).
pub fn normalize_sra_keys(raw: &[String]) -> Vec<String> {
    let mut set = BTreeSet::new();
    for item in raw {
        for tok in split_multi_tokens(item) {
            set.insert(tok);
        }
    }
    set.into_iter().collect()
}

/// Splits and normalizes BioProject / study accessions.
pub fn normalize_bioproject_keys(raw: &[String]) -> Vec<String> {
    let mut set = BTreeSet::new();
    for item in raw {
        for tok in split_multi_tokens(item) {
            let clean = tok.strip_prefix("BIOPROJECT:").unwrap_or(&tok);
            if !is_placeholder(clean) {
                set.insert(clean.to_string());
            }
        }
    }
    set.into_iter().collect()
}

/// Splits and normalizes PubMed identifiers (strips optional `PMID:` prefix).
pub fn normalize_pmid_keys(raw: &[String]) -> Vec<String> {
    let mut set = BTreeSet::new();
    for item in raw {
        for tok in split_multi_tokens(item) {
            let clean = tok.strip_prefix("PMID:").unwrap_or(&tok);
            if !is_placeholder(clean) {
                set.insert(clean.to_string());
            }
        }
    }
    set.into_iter().collect()
}

/// Normalizes a strain label under R0: trim, collapse runs of whitespace to one
/// ASCII space, and uppercase.
pub fn normalize_strain(raw: Option<&str>) -> Option<String> {
    let raw = raw?;
    let collapsed: Vec<&str> = raw.split_whitespace().collect();
    if collapsed.is_empty() {
        return None;
    }
    let up = collapsed.join(" ").to_uppercase();
    (!is_placeholder(&up)).then_some(up)
}

/// Normalizes a longitudinal patient / subject identifier.
pub fn normalize_subject_id(raw: Option<&str>) -> Option<String> {
    let raw = raw?;
    let collapsed: Vec<&str> = raw.split_whitespace().collect();
    if collapsed.is_empty() {
        return None;
    }
    let up = collapsed.join(" ").to_uppercase();
    (!is_placeholder(&up)).then_some(up)
}

/// Two-word species name under R2 (`src/curate.rs` `species_of`): the first two
/// whitespace-separated words of `species`, else `genome_name`, when there are
/// two words and the first starts with an uppercase letter.
pub fn normalize_species(species: Option<&str>, genome_name: Option<&str>) -> Option<String> {
    for candidate in [species, genome_name].into_iter().flatten() {
        let words: Vec<&str> = candidate.split_whitespace().take(2).collect();
        if words.len() == 2 && words[0].chars().next().is_some_and(char::is_uppercase) {
            return Some(words.join(" "));
        }
    }
    None
}

// ---------------------------------------------------------------- union-find & synonym map

struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
        }
    }

    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] != x {
            self.parent[x] = self.parent[self.parent[x]];
            x = self.parent[x];
        }
        x
    }

    fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra != rb {
            if ra < rb {
                self.parent[rb] = ra;
            } else {
                self.parent[ra] = rb;
            }
        }
    }
}

#[derive(Default)]
struct EquivDomain {
    to_id: FxHashMap<String, usize>,
    keys: Vec<String>,
    uf: Vec<usize>,
    canonical: FxHashMap<String, String>,
}

impl EquivDomain {
    fn intern(&mut self, key: String) -> usize {
        if let Some(&id) = self.to_id.get(&key) {
            return id;
        }
        let id = self.keys.len();
        self.keys.push(key.clone());
        self.uf.push(id);
        self.to_id.insert(key, id);
        id
    }

    fn find(&mut self, mut x: usize) -> usize {
        while self.uf[x] != x {
            self.uf[x] = self.uf[self.uf[x]];
            x = self.uf[x];
        }
        x
    }

    fn link(&mut self, a: String, b: String) {
        let ia = self.intern(a);
        let ib = self.intern(b);
        let ra = self.find(ia);
        let rb = self.find(ib);
        if ra != rb {
            if ra < rb {
                self.uf[rb] = ra;
            } else {
                self.uf[ra] = rb;
            }
        }
    }

    fn finalize(&mut self) {
        let n = self.keys.len();
        let mut groups: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
        for i in 0..n {
            let r = self.find(i);
            groups.entry(r).or_default().push(i);
        }
        self.canonical.clear();
        for members in groups.values() {
            // Prefer BioSample/BioProject canonical prefix if present, then
            // lexicographically smallest key for order-independent determinism.
            let rep = members
                .iter()
                .map(|&i| self.keys[i].as_str())
                .min_by(|a, b| {
                    let rank = |s: &str| {
                        if s.starts_with("SAM") || s.starts_with("PRJ") || s.starts_with("ASM:") {
                            0u8
                        } else {
                            1u8
                        }
                    };
                    rank(a).cmp(&rank(b)).then_with(|| a.cmp(b))
                })
                .unwrap()
                .to_string();
            for &i in members {
                self.canonical.insert(self.keys[i].clone(), rep.clone());
            }
        }
    }

    fn resolve<'a>(&'a self, key: &'a str) -> &'a str {
        self.canonical.get(key).map_or(key, String::as_str)
    }
}

/// Pluggable INSDC/ENA/NCBI accession alias equivalence classes.
///
/// Resolves `SAMEA`<->`ERS`, `SAMN`<->`SRS`, `SAMD`<->`DRS`, run-to-sample
/// links (`ERR`/`SRR`/`DRR` <-> BioSample), `PRJEB`<->`ERP`, `PRJNA`<->`SRP`,
/// `PRJDB`<->`DRP`, and assembly aliases via transitive Union-Find closure.
#[derive(Default)]
pub struct SynonymMap {
    sample: EquivDomain,
    study: EquivDomain,
    assembly: EquivDomain,
    pair_count: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AliasDomain {
    Sample,
    Study,
    Assembly,
}

impl SynonymMap {
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of raw equivalence links inserted into the map.
    pub fn pair_count(&self) -> usize {
        self.pair_count
    }

    /// Adds an equivalence link `(a, b)` in `domain`. Call [`Self::finalize`]
    /// before querying if building manually.
    pub fn add_pair(&mut self, domain: AliasDomain, a: &str, b: &str) {
        match domain {
            AliasDomain::Sample => {
                let na = a.trim().to_ascii_uppercase();
                let nb = b.trim().to_ascii_uppercase();
                let na = na.strip_prefix("BIOSAMPLE:").unwrap_or(&na);
                let nb = nb.strip_prefix("BIOSAMPLE:").unwrap_or(&nb);
                if !is_placeholder(na) && !is_placeholder(nb) && na != nb {
                    self.sample.link(na.to_string(), nb.to_string());
                    self.pair_count += 1;
                }
            }
            AliasDomain::Study => {
                let na = a.trim().to_ascii_uppercase();
                let nb = b.trim().to_ascii_uppercase();
                let na = na.strip_prefix("BIOPROJECT:").unwrap_or(&na);
                let nb = nb.strip_prefix("BIOPROJECT:").unwrap_or(&nb);
                if !is_placeholder(na) && !is_placeholder(nb) && na != nb {
                    self.study.link(na.to_string(), nb.to_string());
                    self.pair_count += 1;
                }
            }
            AliasDomain::Assembly => {
                if let (Some(na), Some(nb)) =
                    (normalize_assembly_token(a), normalize_assembly_token(b))
                    && na != nb
                {
                    self.assembly.link(na, nb);
                    self.pair_count += 1;
                }
            }
        }
    }

    /// Computes canonical representatives for all equivalence classes.
    pub fn finalize(&mut self) {
        self.sample.finalize();
        self.study.finalize();
        self.assembly.finalize();
    }

    pub fn resolve_sample<'a>(&'a self, key: &'a str) -> &'a str {
        self.sample.resolve(key)
    }

    pub fn resolve_study<'a>(&'a self, key: &'a str) -> &'a str {
        self.study.resolve(key)
    }

    pub fn resolve_assembly<'a>(&'a self, key: &'a str) -> &'a str {
        self.assembly.resolve(key)
    }

    /// Loads one or more TSV or JSONL synonym maps, returning the finalized
    /// map and a deterministic SHA-256 receipt over the input file bytes (or
    /// `None` when `paths` is empty).
    pub fn load(paths: &[PathBuf]) -> Result<(Self, Option<String>)> {
        if paths.is_empty() {
            return Ok((Self::new(), None));
        }
        let mut map = Self::new();
        let mut digest_ctx = ring::digest::Context::new(&ring::digest::SHA256);
        for path in paths {
            let bytes = std::fs::read(path)
                .with_context(|| format!("read synonym map {}", path.display()))?;
            digest_ctx.update(&(bytes.len() as u64).to_le_bytes());
            digest_ctx.update(&bytes);
            let text = std::str::from_utf8(&bytes)
                .with_context(|| format!("synonym map {} is not UTF-8", path.display()))?;
            let first_non_empty = text
                .lines()
                .map(str::trim)
                .find(|l| !l.is_empty() && !l.starts_with('#'));
            let is_jsonl =
                path.extension().and_then(|e| e.to_str()).is_some_and(|e| {
                    e.eq_ignore_ascii_case("jsonl") || e.eq_ignore_ascii_case("json")
                }) || first_non_empty.is_some_and(|l| l.starts_with('{'));
            if is_jsonl {
                map.parse_jsonl(text, path)?;
            } else {
                map.parse_tsv(text, path)?;
            }
        }
        map.finalize();
        let receipt = hex_bytes(digest_ctx.finish().as_ref());
        Ok((map, Some(receipt)))
    }

    fn parse_tsv(&mut self, text: &str, path: &Path) -> Result<()> {
        for (idx, raw_line) in text.lines().enumerate() {
            let line = raw_line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').map(str::trim).collect();
            if idx == 0 && looks_like_synonym_header(&cols) {
                continue;
            }
            match cols.as_slice() {
                [kind, a, b, ..] => {
                    let domain = parse_alias_domain(kind, a, b).with_context(|| {
                        format!(
                            "{}:{}: unknown synonym kind {kind:?}",
                            path.display(),
                            idx + 1
                        )
                    })?;
                    self.add_pair(domain, a, b);
                }
                [a, b] => {
                    let domain = infer_alias_domain(a, b);
                    self.add_pair(domain, a, b);
                }
                _ => bail!(
                    "{}:{}: expected 2 or 3 tab-separated columns",
                    path.display(),
                    idx + 1
                ),
            }
        }
        Ok(())
    }

    fn parse_jsonl(&mut self, text: &str, path: &Path) -> Result<()> {
        for (idx, raw_line) in text.lines().enumerate() {
            let line = raw_line.trim();
            if line.is_empty() {
                continue;
            }
            let v: Value = serde_json::from_str(line)
                .with_context(|| format!("{}:{}", path.display(), idx + 1))?;
            let kind = v.get("kind").and_then(Value::as_str).unwrap_or("");
            if let Some(arr) = v.get("accessions").and_then(Value::as_array) {
                let items: Vec<&str> = arr.iter().filter_map(Value::as_str).collect();
                if items.len() >= 2 {
                    let domain = if kind.is_empty() {
                        infer_alias_domain(items[0], items[1])
                    } else {
                        parse_alias_domain(kind, items[0], items[1]).with_context(|| {
                            format!(
                                "{}:{}: unknown synonym kind {kind:?}",
                                path.display(),
                                idx + 1
                            )
                        })?
                    };
                    for pair in items.windows(2) {
                        self.add_pair(domain, pair[0], pair[1]);
                    }
                }
            } else if let (Some(a), Some(b)) = (
                v.get("a")
                    .or_else(|| v.get("accession_a"))
                    .and_then(Value::as_str),
                v.get("b")
                    .or_else(|| v.get("accession_b"))
                    .and_then(Value::as_str),
            ) {
                let domain = if kind.is_empty() {
                    infer_alias_domain(a, b)
                } else {
                    parse_alias_domain(kind, a, b).with_context(|| {
                        format!(
                            "{}:{}: unknown synonym kind {kind:?}",
                            path.display(),
                            idx + 1
                        )
                    })?
                };
                self.add_pair(domain, a, b);
            } else {
                bail!(
                    "{}:{}: synonym JSONL line must carry `accessions` or `a`/`b` (`accession_a`/`accession_b`)",
                    path.display(),
                    idx + 1
                );
            }
        }
        Ok(())
    }
}

fn looks_like_synonym_header(cols: &[&str]) -> bool {
    let low: Vec<String> = cols.iter().map(|c| c.to_ascii_lowercase()).collect();
    matches!(
        low.first().map(String::as_str),
        Some("kind" | "domain" | "type" | "accession_a" | "alias" | "source" | "from")
    )
}

fn parse_alias_domain(kind: &str, a: &str, b: &str) -> Result<AliasDomain> {
    match kind.trim().to_ascii_lowercase().as_str() {
        "sample" | "biosample" | "sra" | "run" | "run_to_sample" | "sample_run" | "experiment" => {
            Ok(AliasDomain::Sample)
        }
        "study" | "bioproject" | "project" => Ok(AliasDomain::Study),
        "assembly" | "gca_gcf" => Ok(AliasDomain::Assembly),
        "" => Ok(infer_alias_domain(a, b)),
        other => bail!("unsupported alias domain {other:?}"),
    }
}

fn infer_alias_domain(a: &str, b: &str) -> AliasDomain {
    for tok in [a, b] {
        let up = tok.trim().to_ascii_uppercase();
        if up.starts_with("PRJ")
            || up.starts_with("ERP")
            || up.starts_with("SRP")
            || up.starts_with("DRP")
        {
            return AliasDomain::Study;
        }
        if up.starts_with("GCA_") || up.starts_with("GCF_") || up.starts_with("ASM:") {
            return AliasDomain::Assembly;
        }
    }
    AliasDomain::Sample
}

// ---------------------------------------------------------------- thresholds

/// Species-specific genomic Mash distance thresholds.
///
/// `near_dup_threshold` separates near-duplicate sequence copies from distinct
/// isolates; `close_relative_threshold` bounds the close-relative diagnostic
/// band (`near_dup_threshold < d <= close_relative_threshold`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct SpeciesThresholds {
    pub near_dup_threshold: f64,
    pub close_relative_threshold: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub basis: Option<String>,
}

impl SpeciesThresholds {
    pub fn validate(&self, label: &str) -> Result<()> {
        ensure!(
            self.near_dup_threshold.is_finite()
                && self.close_relative_threshold.is_finite()
                && self.near_dup_threshold >= 0.0
                && self.close_relative_threshold <= 1.0
                && self.near_dup_threshold <= self.close_relative_threshold,
            "invalid thresholds for {label}: require 0.0 <= near_dup ({}) <= close_relative ({}) <= 1.0",
            self.near_dup_threshold,
            self.close_relative_threshold
        );
        Ok(())
    }
}

/// On-disk JSON format for `--thresholds`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ThresholdFile {
    #[serde(default = "default_thresholds_schema")]
    pub schema: String,
    #[serde(default)]
    pub default: Option<SpeciesThresholds>,
    #[serde(default)]
    pub species: BTreeMap<String, SpeciesThresholds>,
}

fn default_thresholds_schema() -> String {
    THRESHOLDS_SCHEMA.to_string()
}

/// Resolved threshold policy combining species-specific overrides with a
/// validated fallback.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ThresholdPolicy {
    pub default: SpeciesThresholds,
    /// Map keyed by normalized species slug (`curate::slug`).
    pub species: BTreeMap<String, SpeciesThresholds>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub config_sha256: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ResolvedThresholds {
    pub near_dup_threshold: f64,
    pub close_relative_threshold: f64,
    /// `"species"` when matched in `--thresholds`, `"default"` otherwise.
    pub source: String,
}

impl ThresholdPolicy {
    pub fn from_defaults(near_dup: f64, close_relative: f64) -> Result<Self> {
        let default = SpeciesThresholds {
            near_dup_threshold: near_dup,
            close_relative_threshold: close_relative,
            basis: Some(
                "cli fallback; validate species-specific thresholds on permitted development records"
                    .to_string(),
            ),
        };
        default.validate("default")?;
        Ok(ThresholdPolicy {
            default,
            species: BTreeMap::new(),
            config_sha256: None,
        })
    }

    pub fn load(
        path: Option<&Path>,
        fallback_near_dup: f64,
        fallback_close_relative: f64,
    ) -> Result<Self> {
        let mut policy = Self::from_defaults(fallback_near_dup, fallback_close_relative)?;
        let Some(path) = path else {
            return Ok(policy);
        };
        let bytes =
            std::fs::read(path).with_context(|| format!("read thresholds {}", path.display()))?;
        let file: ThresholdFile = serde_json::from_slice(&bytes)
            .with_context(|| format!("parse thresholds {}", path.display()))?;
        ensure!(
            file.schema == THRESHOLDS_SCHEMA,
            "{}: unsupported thresholds schema {:?}, expected {THRESHOLDS_SCHEMA:?}",
            path.display(),
            file.schema
        );
        if let Some(d) = file.default {
            d.validate("default")?;
            policy.default = d;
        }
        for (name, t) in file.species {
            t.validate(&name)?;
            let slug = curate::slug(&name);
            ensure!(
                !slug.is_empty(),
                "{}: empty species key {name:?}",
                path.display()
            );
            policy.species.insert(slug, t);
        }
        policy.config_sha256 = Some(sha256(&bytes));
        Ok(policy)
    }

    pub fn for_species(&self, species_slug: Option<&str>) -> ResolvedThresholds {
        if let Some(slug) = species_slug
            && let Some(t) = self.species.get(slug)
        {
            return ResolvedThresholds {
                near_dup_threshold: t.near_dup_threshold,
                close_relative_threshold: t.close_relative_threshold,
                source: "species".to_string(),
            };
        }
        ResolvedThresholds {
            near_dup_threshold: self.default.near_dup_threshold,
            close_relative_threshold: self.default.close_relative_threshold,
            source: "default".to_string(),
        }
    }
}

// ---------------------------------------------------------------- fast inverted MinHash index

/// Inverted posting index over a deterministic prefix of each bottom-`s`
/// MinHash sketch for sub-millisecond candidate filtering before exact
/// [`cluster::mash_distance`] verification.
///
/// For any pair whose sketches both hold at least `index_hashes` values, every
/// shared hash within the `index_hashes` smallest elements of their union is
/// guaranteed to appear in both prefixes. When the maximum mismatch count
/// allowed by `max_distance` is below `index_hashes`, the hit-count prefilter
/// is exact (zero false negatives). Sketches shorter than `index_hashes` are
/// checked directly.
pub struct SketchIndex {
    k: usize,
    sketch_size: usize,
    index_hashes: usize,
    sketches: Vec<Vec<u64>>,
    postings: FxHashMap<u64, Vec<u32>>,
    short_refs: Vec<u32>,
}

/// Thread-local scratch space for zero-allocation [`SketchIndex::query_with_scratch`].
pub struct QueryScratch {
    counts: Vec<u32>,
    touched: Vec<u32>,
}

impl QueryScratch {
    pub fn new(n_refs: usize) -> Self {
        QueryScratch {
            counts: vec![0; n_refs],
            touched: Vec::new(),
        }
    }
}

impl SketchIndex {
    pub fn build(
        sketches: Vec<Vec<u64>>,
        k: usize,
        sketch_size: usize,
        index_hashes: usize,
    ) -> Self {
        let index_hashes = index_hashes.max(1).min(sketch_size.max(1));
        let mut postings: FxHashMap<u64, Vec<u32>> = FxHashMap::default();
        let mut short_refs = Vec::new();
        for (idx, sketch) in sketches.iter().enumerate() {
            if sketch.is_empty() {
                continue;
            }
            if sketch.len() < index_hashes {
                short_refs.push(idx as u32);
            }
            for &h in sketch.iter().take(index_hashes) {
                postings.entry(h).or_default().push(idx as u32);
            }
        }
        SketchIndex {
            k,
            sketch_size,
            index_hashes,
            sketches,
            postings,
            short_refs,
        }
    }

    pub fn len(&self) -> usize {
        self.sketches.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sketches.is_empty()
    }

    /// Returns all `(ref_index, mash_distance)` with `mash_distance <= max_distance`,
    /// sorted ascending by `(mash_distance, ref_index)`.
    pub fn query(&self, sketch: &[u64], max_distance: f64) -> Vec<(usize, f64)> {
        let mut scratch = QueryScratch::new(self.sketches.len());
        self.query_with_scratch(sketch, max_distance, &mut scratch)
    }

    /// [`Self::query`] reusing caller-owned scratch counters.
    pub fn query_with_scratch(
        &self,
        sketch: &[u64],
        max_distance: f64,
        scratch: &mut QueryScratch,
    ) -> Vec<(usize, f64)> {
        if sketch.is_empty()
            || self.sketches.is_empty()
            || !max_distance.is_finite()
            || max_distance < 0.0
        {
            return Vec::new();
        }
        if scratch.counts.len() < self.sketches.len() {
            scratch.counts.resize(self.sketches.len(), 0);
        }
        let mut hits = Vec::new();
        if max_distance >= 1.0 {
            for (i, ref_s) in self.sketches.iter().enumerate() {
                if ref_s.is_empty() {
                    continue;
                }
                let d = cluster::mash_distance(
                    cluster::jaccard(sketch, ref_s, self.sketch_size),
                    self.k,
                );
                if d <= max_distance {
                    hits.push((i, d));
                }
            }
            hits.sort_by(|a, b| a.1.total_cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
            return hits;
        }

        let b = self.index_hashes;
        let e = (-(self.k as f64) * max_distance).exp();
        let j_min = e / (2.0 - e);
        let m_max = ((1.0 - j_min).max(0.0) * (self.sketch_size as f64)).ceil() as usize;
        // Deterministic zero-false-negative threshold when both sketches have >= B
        // values and m_max + 1 < B; otherwise fall back to >= 1 shared prefix hash.
        let min_hits = if sketch.len() >= b && m_max + 1 < b {
            (b - (m_max + 1)) as u32
        } else {
            1u32
        };

        for &h in sketch.iter().take(b) {
            if let Some(list) = self.postings.get(&h) {
                for &ref_idx in list {
                    let c = &mut scratch.counts[ref_idx as usize];
                    if *c == 0 {
                        scratch.touched.push(ref_idx);
                    }
                    *c += 1;
                }
            }
        }

        let mut checked_short: FxHashSet<u32> = FxHashSet::default();
        for &ref_idx in &scratch.touched {
            let count = scratch.counts[ref_idx as usize];
            let idx = ref_idx as usize;
            let ref_s = &self.sketches[idx];
            if count >= min_hits || ref_s.len() < b {
                if ref_s.len() < b {
                    checked_short.insert(ref_idx);
                }
                let d = cluster::mash_distance(
                    cluster::jaccard(sketch, ref_s, self.sketch_size),
                    self.k,
                );
                if d <= max_distance {
                    hits.push((idx, d));
                }
            }
        }
        for &ref_idx in &scratch.touched {
            scratch.counts[ref_idx as usize] = 0;
        }
        scratch.touched.clear();

        for &ref_idx in &self.short_refs {
            if checked_short.contains(&ref_idx) {
                continue;
            }
            let idx = ref_idx as usize;
            let ref_s = &self.sketches[idx];
            let d =
                cluster::mash_distance(cluster::jaccard(sketch, ref_s, self.sketch_size), self.k);
            if d <= max_distance {
                hits.push((idx, d));
            }
        }

        hits.sort_by(|a, b| a.1.total_cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
        hits
    }
}

// ---------------------------------------------------------------- flexible metadata ingestion

/// One genome metadata/provenance input row (compatible with `genomes.jsonl`,
/// sidecar `genome_provenance.jsonl`, and external portal query JSONL).
#[derive(Clone, Debug, Default)]
pub struct MetadataRow {
    pub genome_id: String,
    pub species: Option<String>,
    pub species_slug: Option<String>,
    pub assembly_accessions: Vec<String>,
    pub biosample_accessions: Vec<String>,
    pub sra_accessions: Vec<String>,
    pub bioproject_accessions: Vec<String>,
    pub pmids: Vec<String>,
    pub strain: Option<String>,
    pub subject_id: Option<String>,
    pub fasta_path: Option<PathBuf>,
    pub fnv64_k31: Option<String>,
    pub sha256_k31: Option<String>,
    pub assembly_sha256: Option<String>,
    pub n_kmers: Option<u64>,
    pub contigs: Option<u32>,
    pub bp: Option<u64>,
    pub sketch_k21: Option<Vec<u64>>,
    pub sketch_k31: Option<Vec<u64>>,
}

fn value_strings(v: &Value) -> Vec<String> {
    match v {
        Value::Null => Vec::new(),
        Value::String(s) => {
            let t = s.trim();
            if t.is_empty() {
                Vec::new()
            } else {
                vec![t.to_string()]
            }
        }
        Value::Number(n) => vec![n.to_string()],
        Value::Bool(_) | Value::Object(_) => Vec::new(),
        Value::Array(items) => items.iter().flat_map(value_strings).collect(),
    }
}

fn first_string(obj: &serde_json::Map<String, Value>, keys: &[&str]) -> Option<String> {
    for &key in keys {
        if let Some(v) = obj.get(key)
            && let Some(s) = value_strings(v).into_iter().next()
        {
            return Some(s);
        }
    }
    None
}

fn collect_strings(obj: &serde_json::Map<String, Value>, keys: &[&str]) -> Vec<String> {
    let mut out = Vec::new();
    for &key in keys {
        if let Some(v) = obj.get(key) {
            out.extend(value_strings(v));
        }
    }
    out
}

fn optional_u64(obj: &serde_json::Map<String, Value>, key: &str) -> Option<u64> {
    obj.get(key).and_then(|v| match v {
        Value::Number(n) => n.as_u64(),
        Value::String(s) => s.trim().parse::<u64>().ok(),
        _ => None,
    })
}

fn optional_u64_vec(obj: &serde_json::Map<String, Value>, key: &str) -> Result<Option<Vec<u64>>> {
    let Some(v) = obj.get(key) else {
        return Ok(None);
    };
    if v.is_null() {
        return Ok(None);
    }
    let arr = v
        .as_array()
        .with_context(|| format!("{key} must be an array of integers"))?;
    let mut out = Vec::with_capacity(arr.len());
    for item in arr {
        let n = item
            .as_u64()
            .with_context(|| format!("{key} contains non-u64 element {item}"))?;
        out.push(n);
    }
    Ok(Some(out))
}

fn parse_metadata_line(line: &str) -> Result<MetadataRow> {
    let v: Value = serde_json::from_str(line)?;
    let obj = v.as_object().context("JSONL row must be a JSON object")?;
    let genome_id = first_string(obj, &["genome_id", "sample_id", "id"])
        .context("row has no genome_id, sample_id, or id")?;
    let raw_species = first_string(obj, &["species"]);
    let raw_genome_name = first_string(obj, &["genome_name"]);
    let species = normalize_species(raw_species.as_deref(), raw_genome_name.as_deref())
        .or_else(|| raw_species.filter(|s| !is_placeholder(&s.to_ascii_uppercase())));
    let species_slug = species
        .as_deref()
        .map(curate::slug)
        .filter(|s| !s.is_empty());

    let assembly_accessions = normalize_assembly_keys(&collect_strings(
        obj,
        &[
            "assembly_accession",
            "assembly_accessions",
            "assembly_aliases",
            "refseq_assembly_accession",
        ],
    ));
    let biosample_accessions = normalize_biosample_keys(&collect_strings(
        obj,
        &[
            "biosample_accession",
            "biosample_accessions",
            "biosample_aliases",
            "biosample",
        ],
    ));
    let sra_accessions = normalize_sra_keys(&collect_strings(
        obj,
        &[
            "sra_accession",
            "sra_accessions",
            "run_accession",
            "ena_run",
        ],
    ));
    let bioproject_accessions = normalize_bioproject_keys(&collect_strings(
        obj,
        &[
            "bioproject_accession",
            "bioproject_accessions",
            "bioproject_aliases",
            "study_accession",
        ],
    ));
    let pmids = normalize_pmid_keys(&collect_strings(
        obj,
        &["pmid", "pmids", "publication", "publications", "pubmed_id"],
    ));
    let strain = normalize_strain(first_string(obj, &["strain"]).as_deref());
    let subject_id =
        normalize_subject_id(first_string(obj, &["subject_id", "patient_id"]).as_deref());
    let fasta_path = first_string(obj, &["fasta_path"]).map(PathBuf::from);

    let fnv64_k31 =
        first_string(obj, &["fnv64_k31", "kset_digest"]).map(|s| s.to_ascii_lowercase());
    let sha256_k31 = first_string(obj, &["sha256_k31"]).map(|s| s.to_ascii_lowercase());
    let assembly_sha256 = first_string(obj, &["assembly_sha256"]).map(|s| s.to_ascii_lowercase());
    let n_kmers = optional_u64(obj, "n_kmers");
    let contigs = optional_u64(obj, "contigs").and_then(|c| u32::try_from(c).ok());
    let bp = optional_u64(obj, "bp").or_else(|| optional_u64(obj, "genome_length"));
    let sketch_k21 = optional_u64_vec(obj, "sketch_k21")?;
    let sketch_k31 = optional_u64_vec(obj, "sketch_k31")?;

    Ok(MetadataRow {
        genome_id,
        species,
        species_slug,
        assembly_accessions,
        biosample_accessions,
        sra_accessions,
        bioproject_accessions,
        pmids,
        strain,
        subject_id,
        fasta_path,
        fnv64_k31,
        sha256_k31,
        assembly_sha256,
        n_kmers,
        contigs,
        bp,
        sketch_k21,
        sketch_k31,
    })
}

fn merge_sorted_vec(dst: &mut Vec<String>, src: Vec<String>) {
    if src.is_empty() {
        return;
    }
    let mut set: BTreeSet<String> = dst.drain(..).collect();
    set.extend(src);
    *dst = set.into_iter().collect();
}

fn merge_metadata_row(dst: &mut MetadataRow, src: MetadataRow) {
    if dst.species.is_none() {
        dst.species = src.species;
        dst.species_slug = src.species_slug;
    }
    merge_sorted_vec(&mut dst.assembly_accessions, src.assembly_accessions);
    merge_sorted_vec(&mut dst.biosample_accessions, src.biosample_accessions);
    merge_sorted_vec(&mut dst.sra_accessions, src.sra_accessions);
    merge_sorted_vec(&mut dst.bioproject_accessions, src.bioproject_accessions);
    merge_sorted_vec(&mut dst.pmids, src.pmids);
    if dst.strain.is_none() {
        dst.strain = src.strain;
    }
    if dst.subject_id.is_none() {
        dst.subject_id = src.subject_id;
    }
    if dst.fasta_path.is_none() {
        dst.fasta_path = src.fasta_path;
    }
    if dst.fnv64_k31.is_none() {
        dst.fnv64_k31 = src.fnv64_k31;
    }
    if dst.sha256_k31.is_none() {
        dst.sha256_k31 = src.sha256_k31;
    }
    if dst.assembly_sha256.is_none() {
        dst.assembly_sha256 = src.assembly_sha256;
    }
    if dst.n_kmers.is_none() {
        dst.n_kmers = src.n_kmers;
    }
    if dst.contigs.is_none() {
        dst.contigs = src.contigs;
    }
    if dst.bp.is_none() {
        dst.bp = src.bp;
    }
    if dst.sketch_k21.is_none() {
        dst.sketch_k21 = src.sketch_k21;
    }
    if dst.sketch_k31.is_none() {
        dst.sketch_k31 = src.sketch_k31;
    }
}

fn read_metadata_jsonl(path: &Path) -> Result<Vec<MetadataRow>> {
    let reader = BufReader::new(
        File::open(path).with_context(|| format!("open metadata JSONL {}", path.display()))?,
    );
    let mut rows = Vec::new();
    for (idx, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        rows.push(
            parse_metadata_line(trimmed)
                .with_context(|| format!("{}:{}", path.display(), idx + 1))?,
        );
    }
    Ok(rows)
}

fn resolve_fasta_for(
    genome_id: &str,
    explicit: Option<&Path>,
    fasta_dir: Option<&Path>,
) -> Result<Option<PathBuf>> {
    if let Some(p) = explicit {
        if p.is_absolute() {
            ensure!(
                p.is_file(),
                "{genome_id}: explicit fasta_path {} does not exist",
                p.display()
            );
            return Ok(Some(p.to_path_buf()));
        }
        if let Some(dir) = fasta_dir {
            let joined = dir.join(p);
            if joined.is_file() {
                return Ok(Some(joined));
            }
        }
        ensure!(
            p.is_file(),
            "{genome_id}: explicit fasta_path {} does not exist",
            p.display()
        );
        return Ok(Some(p.to_path_buf()));
    }
    let Some(dir) = fasta_dir else {
        return Ok(None);
    };
    for ext in ["fna.gz", "fna", "fa.gz", "fa", "fasta.gz", "fasta"] {
        let candidate = dir.join(format!("{genome_id}.{ext}"));
        if candidate.is_file() {
            return Ok(Some(candidate));
        }
    }
    Ok(None)
}

// ---------------------------------------------------------------- reference index (`grm dedup index`)

/// One canonical reference record stored in `<index>/records.jsonl`.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct IndexRecord {
    pub genome_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub species: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub species_slug: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub assembly_accessions: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub biosample_accessions: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub sra_accessions: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub bioproject_accessions: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub pmids: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strain: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub subject_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fnv64_k31: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sha256_k31: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assembly_sha256: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub n_kmers: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub contigs: Option<u32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bp: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sketch_k21: Option<Vec<u64>>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sketch_k31: Option<Vec<u64>>,
}

impl IndexRecord {
    /// True when the record holds no assembly, BioSample, or SRA accession.
    pub fn missing_accessions(&self) -> bool {
        self.assembly_accessions.is_empty()
            && self.biosample_accessions.is_empty()
            && self.sra_accessions.is_empty()
    }

    /// True when the record carries at least one sequence-derived fingerprint
    /// (a non-empty MinHash sketch or canonical k-mer digest).
    pub fn has_sequence_fingerprint(&self) -> bool {
        self.sketch_k21.as_ref().is_some_and(|s| !s.is_empty())
            || self.sha256_k31.is_some()
            || self.fnv64_k31.is_some()
            || self.assembly_sha256.is_some()
    }
}

/// Summary coverage counts for a reference index.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct IndexCounts {
    pub records: usize,
    pub with_assembly_accession: usize,
    pub with_biosample_accession: usize,
    pub with_sra_accession: usize,
    pub with_bioproject: usize,
    pub with_pmid: usize,
    pub with_strain: usize,
    pub with_subject_id: usize,
    pub accession_gap: usize,
    pub with_fnv64_k31: usize,
    pub with_sha256_k31: usize,
    pub with_assembly_sha256: usize,
    pub with_sketch_k21: usize,
    pub with_sketch_k31: usize,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct SpeciesIndexCounts {
    pub records: usize,
    pub accession_gap: usize,
    pub with_sketch_k21: usize,
    pub with_kset_digest: usize,
}

/// Versioned reference index manifest (`<index>/manifest.json`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct IndexManifest {
    pub schema: String,
    pub source_label: String,
    pub sealed_set_status: String,
    pub k_digest: usize,
    pub k_sketch: usize,
    pub sketch_size: usize,
    pub records_file: String,
    pub records_bytes: u64,
    pub records_sha256: String,
    pub index_sha256: String,
    pub counts: IndexCounts,
    pub per_species: BTreeMap<String, SpeciesIndexCounts>,
}

/// Configuration for [`build_index`].
#[derive(Clone, Debug)]
pub struct IndexConfig {
    pub metadata: Vec<PathBuf>,
    pub provenance: Option<PathBuf>,
    pub genome_ids: Option<PathBuf>,
    pub fasta_dir: Option<PathBuf>,
    pub counted: Vec<PathBuf>,
    pub sketches_dir: Option<PathBuf>,
    pub source_label: String,
    pub k_digest: usize,
    pub k_sketch: usize,
    pub sketch_size: usize,
    pub out_dir: PathBuf,
}

impl Default for IndexConfig {
    fn default() -> Self {
        IndexConfig {
            metadata: Vec::new(),
            provenance: None,
            genome_ids: None,
            fasta_dir: None,
            counted: Vec::new(),
            sketches_dir: None,
            source_label: "reference".to_string(),
            k_digest: DEFAULT_K_DIGEST,
            k_sketch: DEFAULT_K_SKETCH,
            sketch_size: DEFAULT_SKETCH_SIZE,
            out_dir: PathBuf::new(),
        }
    }
}

/// A loaded and receipt-verified reference index.
#[derive(Clone, Debug)]
pub struct LoadedIndex {
    pub path: PathBuf,
    pub manifest: IndexManifest,
    pub records: Vec<IndexRecord>,
}

fn load_counted_paths(
    paths: &[PathBuf],
    expected_k: usize,
) -> Result<FxHashMap<String, CountedInfo>> {
    let mut files = Vec::new();
    for p in paths {
        if p.is_file() {
            files.push(p.clone());
        } else if p.is_dir() {
            let direct = p.join("counted.json");
            if direct.is_file() {
                files.push(direct);
            } else {
                let mut children = Vec::new();
                for entry in std::fs::read_dir(p)
                    .with_context(|| format!("read counted directory {}", p.display()))?
                {
                    let entry = entry?;
                    let c = entry.path().join("counted.json");
                    if c.is_file() {
                        children.push(c);
                    }
                }
                children.sort();
                files.extend(children);
            }
        } else {
            bail!("--counted path does not exist: {}", p.display());
        }
    }
    let mut map: FxHashMap<String, CountedInfo> = FxHashMap::default();
    for file in files {
        let bytes = std::fs::read(&file).with_context(|| format!("read {}", file.display()))?;
        let summary: CountSummary =
            serde_json::from_slice(&bytes).with_context(|| format!("parse {}", file.display()))?;
        ensure!(
            summary.k == expected_k,
            "{}: counted at k={}, index expects k_digest={expected_k}",
            file.display(),
            summary.k
        );
        for g in summary.genomes {
            let info = CountedInfo {
                fnv64_k31: g.kset_digest,
                n_kmers: g.n_kmers,
                contigs: g.stats.contigs,
                bp: g.stats.length,
                sketch_k31: g.sketch,
            };
            map.insert(g.genome_id, info);
        }
    }
    Ok(map)
}

struct CountedInfo {
    fnv64_k31: Option<String>,
    n_kmers: u64,
    contigs: u32,
    bp: u64,
    sketch_k31: Option<Vec<u64>>,
}

fn load_species_sketch(
    dir: &Path,
    genome_id: &str,
    k_sketch: usize,
    sketch_size: usize,
) -> Result<Option<GenomeSketch>> {
    let path = dir.join(format!("{genome_id}.json"));
    if !path.is_file() {
        return Ok(None);
    }
    let bytes = std::fs::read(&path).with_context(|| format!("read {}", path.display()))?;
    let s: GenomeSketch =
        serde_json::from_slice(&bytes).with_context(|| format!("parse {}", path.display()))?;
    ensure!(
        s.k == k_sketch && s.sketch_size == sketch_size && s.genome_id == genome_id,
        "{}: sketch metadata (id={}, k={}, s={}) does not match expected ({genome_id}, k={k_sketch}, s={sketch_size})",
        path.display(),
        s.genome_id,
        s.k,
        s.sketch_size
    );
    Ok(Some(s))
}

fn read_id_allowlist(path: &Path) -> Result<BTreeSet<String>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("read genome id list {}", path.display()))?;
    Ok(text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(String::from)
        .collect())
}

fn compute_index_receipt(
    schema: &str,
    source_label: &str,
    k_digest: usize,
    k_sketch: usize,
    sketch_size: usize,
    records_sha256: &str,
) -> String {
    let payload = format!(
        "{schema}\n{source_label}\n{k_digest}\n{k_sketch}\n{sketch_size}\n{records_sha256}\n"
    );
    sha256(payload.as_bytes())
}

/// Builds a compact versioned reference index (`grm-dedup-index-v1`) on disk
/// without modifying any input metadata file (`genomes.jsonl` stays untouched).
pub fn build_index(cfg: &IndexConfig) -> Result<IndexManifest> {
    ensure!(
        (MIN_K..=MAX_K).contains(&cfg.k_digest),
        "k_digest must be within {MIN_K}..={MAX_K}, got {}",
        cfg.k_digest
    );
    ensure!(
        (MIN_K..=MAX_K).contains(&cfg.k_sketch),
        "k_sketch must be within {MIN_K}..={MAX_K}, got {}",
        cfg.k_sketch
    );
    ensure!(
        !cfg.metadata.is_empty() || cfg.provenance.is_some(),
        "provide at least one --metadata JSONL or --provenance sidecar JSONL"
    );

    let allowlist = cfg
        .genome_ids
        .as_deref()
        .map(read_id_allowlist)
        .transpose()?;

    let mut rows_by_id: FxHashMap<String, MetadataRow> = FxHashMap::default();
    let mut input_files: Vec<&Path> = cfg.metadata.iter().map(PathBuf::as_path).collect();
    if let Some(prov) = cfg.provenance.as_deref() {
        input_files.push(prov);
    }
    for path in input_files {
        for row in read_metadata_jsonl(path)? {
            if let Some(ref allowed) = allowlist
                && !allowed.contains(&row.genome_id)
            {
                continue;
            }
            rows_by_id
                .entry(row.genome_id.clone())
                .and_modify(|existing| merge_metadata_row(existing, row.clone()))
                .or_insert(row);
        }
    }
    if let Some(allowed) = allowlist {
        for id in allowed {
            rows_by_id.entry(id.clone()).or_insert_with(|| MetadataRow {
                genome_id: id,
                ..Default::default()
            });
        }
    }

    let counted_map = load_counted_paths(&cfg.counted, cfg.k_digest)?;
    let mut rows: Vec<MetadataRow> = rows_by_id.into_values().collect();
    rows.sort_by(|a, b| compare_genome_ids(&a.genome_id, &b.genome_id));

    let records: Vec<IndexRecord> = rows
        .into_par_iter()
        .map_init(
            || (Vec::<u64>::new(), Vec::<u64>::new()),
            |(values, scratch), mut row| -> Result<IndexRecord> {
                let fasta = resolve_fasta_for(
                    &row.genome_id,
                    row.fasta_path.as_deref(),
                    cfg.fasta_dir.as_deref(),
                )?;
                let fp: Option<FastaFingerprint> = match fasta {
                    Some(path) => Some(
                        kset::fingerprint_fasta_into(
                            &path,
                            cfg.k_digest,
                            cfg.k_sketch,
                            cfg.sketch_size,
                            values,
                            scratch,
                        )
                        .with_context(|| format!("fingerprint {}", path.display()))?,
                    ),
                    None => None,
                };
                if let Some(fp) = fp {
                    if let Some(c) = counted_map.get(&row.genome_id)
                        && let Some(ref counted_digest) = c.fnv64_k31
                    {
                        ensure!(
                            counted_digest == &fp.fnv64_k31,
                            "{}: FASTA fnv64_k31 ({}) disagrees with counted.json kset_digest ({counted_digest})",
                            row.genome_id,
                            fp.fnv64_k31
                        );
                    }
                    row.fnv64_k31 = Some(fp.fnv64_k31);
                    row.sha256_k31 = Some(fp.sha256_k31);
                    row.assembly_sha256 = Some(fp.assembly_sha256);
                    row.n_kmers = Some(fp.n_kmers);
                    row.contigs = Some(fp.contigs);
                    row.bp = Some(fp.bp);
                    row.sketch_k21 = Some(fp.sketch_k21);
                    row.sketch_k31 = fp.sketch_k31;
                } else {
                    if let Some(c) = counted_map.get(&row.genome_id) {
                        if row.fnv64_k31.is_none() {
                            row.fnv64_k31.clone_from(&c.fnv64_k31);
                        }
                        if row.n_kmers.is_none() {
                            row.n_kmers = Some(c.n_kmers);
                        }
                        if row.contigs.is_none() {
                            row.contigs = Some(c.contigs);
                        }
                        if row.bp.is_none() {
                            row.bp = Some(c.bp);
                        }
                        if row.sketch_k31.is_none() {
                            row.sketch_k31.clone_from(&c.sketch_k31);
                        }
                    }
                    if row.sketch_k21.is_none()
                        && let Some(dir) = cfg.sketches_dir.as_deref()
                        && let Some(sk) =
                            load_species_sketch(dir, &row.genome_id, cfg.k_sketch, cfg.sketch_size)?
                    {
                        if row.bp.is_none() {
                            row.bp = Some(sk.length);
                        }
                        row.sketch_k21 = Some(sk.sketch);
                    }
                }
                Ok(IndexRecord {
                    genome_id: row.genome_id,
                    species: row.species,
                    species_slug: row.species_slug,
                    assembly_accessions: row.assembly_accessions,
                    biosample_accessions: row.biosample_accessions,
                    sra_accessions: row.sra_accessions,
                    bioproject_accessions: row.bioproject_accessions,
                    pmids: row.pmids,
                    strain: row.strain,
                    subject_id: row.subject_id,
                    fnv64_k31: row.fnv64_k31,
                    sha256_k31: row.sha256_k31,
                    assembly_sha256: row.assembly_sha256,
                    n_kmers: row.n_kmers,
                    contigs: row.contigs,
                    bp: row.bp,
                    sketch_k21: row.sketch_k21,
                    sketch_k31: row.sketch_k31,
                })
            },
        )
        .collect::<Result<Vec<_>>>()?;

    let mut counts = IndexCounts {
        records: records.len(),
        ..Default::default()
    };
    let mut per_species: BTreeMap<String, SpeciesIndexCounts> = BTreeMap::new();
    let mut records_bytes = Vec::new();
    for rec in &records {
        counts.with_assembly_accession += usize::from(!rec.assembly_accessions.is_empty());
        counts.with_biosample_accession += usize::from(!rec.biosample_accessions.is_empty());
        counts.with_sra_accession += usize::from(!rec.sra_accessions.is_empty());
        counts.with_bioproject += usize::from(!rec.bioproject_accessions.is_empty());
        counts.with_pmid += usize::from(!rec.pmids.is_empty());
        counts.with_strain += usize::from(rec.strain.is_some());
        counts.with_subject_id += usize::from(rec.subject_id.is_some());
        let gap = rec.missing_accessions();
        counts.accession_gap += usize::from(gap);
        counts.with_fnv64_k31 += usize::from(rec.fnv64_k31.is_some());
        counts.with_sha256_k31 += usize::from(rec.sha256_k31.is_some());
        counts.with_assembly_sha256 += usize::from(rec.assembly_sha256.is_some());
        let has_sk21 = rec.sketch_k21.as_ref().is_some_and(|s| !s.is_empty());
        counts.with_sketch_k21 += usize::from(has_sk21);
        counts.with_sketch_k31 +=
            usize::from(rec.sketch_k31.as_ref().is_some_and(|s| !s.is_empty()));

        let sp_key = rec
            .species_slug
            .clone()
            .unwrap_or_else(|| "unknown".to_string());
        let sp = per_species.entry(sp_key).or_default();
        sp.records += 1;
        sp.accession_gap += usize::from(gap);
        sp.with_sketch_k21 += usize::from(has_sk21);
        sp.with_kset_digest += usize::from(rec.sha256_k31.is_some() || rec.fnv64_k31.is_some());

        serde_json::to_writer(&mut records_bytes, rec)?;
        records_bytes.push(b'\n');
    }

    let records_sha256 = sha256(&records_bytes);
    let index_sha256 = compute_index_receipt(
        INDEX_SCHEMA,
        &cfg.source_label,
        cfg.k_digest,
        cfg.k_sketch,
        cfg.sketch_size,
        &records_sha256,
    );
    let manifest = IndexManifest {
        schema: INDEX_SCHEMA.to_string(),
        source_label: cfg.source_label.clone(),
        sealed_set_status: SEALED_SET_STATUS.to_string(),
        k_digest: cfg.k_digest,
        k_sketch: cfg.k_sketch,
        sketch_size: cfg.sketch_size,
        records_file: "records.jsonl".to_string(),
        records_bytes: records_bytes.len() as u64,
        records_sha256,
        index_sha256,
        counts,
        per_species,
    };

    std::fs::create_dir_all(&cfg.out_dir)
        .with_context(|| format!("create index dir {}", cfg.out_dir.display()))?;
    let records_path = cfg.out_dir.join("records.jsonl");
    let manifest_path = cfg.out_dir.join("manifest.json");
    write_atomic(&records_path, &records_bytes)?;
    write_atomic(&manifest_path, &serde_json::to_vec_pretty(&manifest)?)?;
    Ok(manifest)
}

/// Loads a versioned reference index from its directory or `manifest.json`,
/// verifying schema and SHA-256 receipts.
pub fn load_index(path: &Path) -> Result<LoadedIndex> {
    let (dir, manifest_path) = if path.is_dir() {
        (path.to_path_buf(), path.join("manifest.json"))
    } else {
        let dir = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."))
            .to_path_buf();
        (dir, path.to_path_buf())
    };
    let manifest_bytes = std::fs::read(&manifest_path)
        .with_context(|| format!("read index manifest {}", manifest_path.display()))?;
    let manifest: IndexManifest = serde_json::from_slice(&manifest_bytes)
        .with_context(|| format!("parse index manifest {}", manifest_path.display()))?;
    ensure!(
        manifest.schema == INDEX_SCHEMA,
        "{}: unsupported index schema {:?}, expected {INDEX_SCHEMA:?}",
        manifest_path.display(),
        manifest.schema
    );
    let records_path = dir.join(&manifest.records_file);
    let records_bytes = std::fs::read(&records_path)
        .with_context(|| format!("read index records {}", records_path.display()))?;
    ensure!(
        records_bytes.len() as u64 == manifest.records_bytes
            && sha256(&records_bytes) == manifest.records_sha256,
        "{}: records file SHA-256 or size does not match {}",
        records_path.display(),
        manifest_path.display()
    );
    let expected_receipt = compute_index_receipt(
        &manifest.schema,
        &manifest.source_label,
        manifest.k_digest,
        manifest.k_sketch,
        manifest.sketch_size,
        &manifest.records_sha256,
    );
    ensure!(
        manifest.index_sha256 == expected_receipt,
        "{}: index_sha256 receipt mismatch",
        manifest_path.display()
    );
    let text = std::str::from_utf8(&records_bytes)
        .with_context(|| format!("{} is not UTF-8", records_path.display()))?;
    let mut records = Vec::with_capacity(manifest.counts.records);
    for (idx, line) in text.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        records.push(
            serde_json::from_str::<IndexRecord>(trimmed)
                .with_context(|| format!("{}:{}", records_path.display(), idx + 1))?,
        );
    }
    ensure!(
        records.len() == manifest.counts.records,
        "{}: manifest lists {} records, file holds {}",
        records_path.display(),
        manifest.counts.records,
        records.len()
    );
    Ok(LoadedIndex {
        path: dir,
        manifest,
        records,
    })
}

// ---------------------------------------------------------------- screening (`grm dedup screen`)

/// Configuration for [`screen`].
#[derive(Clone, Debug)]
pub struct ScreenConfig {
    pub references: Vec<PathBuf>,
    pub query: PathBuf,
    pub query_fasta_dir: Option<PathBuf>,
    pub query_sketches_dir: Option<PathBuf>,
    pub synonym_maps: Vec<PathBuf>,
    pub thresholds: Option<PathBuf>,
    pub fallback_near_dup: f64,
    pub fallback_close_relative: f64,
    pub k_digest: usize,
    pub k_sketch: usize,
    pub sketch_size: usize,
    pub index_hashes: usize,
    pub out_dir: PathBuf,
}

impl Default for ScreenConfig {
    fn default() -> Self {
        ScreenConfig {
            references: Vec::new(),
            query: PathBuf::new(),
            query_fasta_dir: None,
            query_sketches_dir: None,
            synonym_maps: Vec::new(),
            thresholds: None,
            fallback_near_dup: DEFAULT_NEAR_DUP_THRESHOLD,
            fallback_close_relative: DEFAULT_CLOSE_RELATIVE_THRESHOLD,
            k_digest: DEFAULT_K_DIGEST,
            k_sketch: DEFAULT_K_SKETCH,
            sketch_size: DEFAULT_SKETCH_SIZE,
            index_hashes: DEFAULT_INDEX_HASHES,
            out_dir: PathBuf::new(),
        }
    }
}

/// Primary disposition assigned to each query genome in the sequential funnel.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Disposition {
    #[serde(rename = "certified_clean")]
    CertifiedClean,
    #[serde(rename = "excluded_C4_ref_metadata")]
    ExcludedRefMetadata,
    #[serde(rename = "excluded_EXACT_ASSEMBLY")]
    ExcludedExactAssembly,
    #[serde(rename = "excluded_EXACT_KSET")]
    ExcludedExactKset,
    #[serde(rename = "excluded_D3_mash_near_dup")]
    ExcludedMashNearDup,
    #[serde(rename = "excluded_D4_within_batch_dup")]
    ExcludedWithinBatchDup,
    #[serde(rename = "excluded_UNRESOLVED_OVERLAP")]
    ExcludedUnresolvedOverlap,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct RefKeyHit {
    pub ref_id: String,
    pub source_label: String,
    pub matched_key: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub canonical_key: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RefDistanceHit {
    pub ref_id: String,
    pub source_label: String,
    pub mash_distance: f64,
    pub ref_missing_accessions: bool,
}

/// Separate reference-matching evidence layers for one query genome.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct RefMatchLayers {
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub d1_assembly: Vec<RefKeyHit>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub d2_biosample_sra: Vec<RefKeyHit>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub strain_same_species: Vec<RefKeyHit>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub exact_assembly: Vec<RefKeyHit>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub exact_kset: Vec<RefKeyHit>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub near_dup: Vec<RefDistanceHit>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub close_relative: Vec<RefDistanceHit>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub shared_bioproject: Vec<RefKeyHit>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub shared_pmid: Vec<RefKeyHit>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub shared_subject: Vec<RefKeyHit>,
}

/// One row of `<out>/decisions.jsonl`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct QueryDecision {
    pub genome_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub species: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub species_slug: Option<String>,
    pub disposition: Disposition,
    pub certified_clean: bool,
    pub matched_layers: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub study_tags: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub diagnostic_tags: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub unresolved_reasons: Vec<String>,
    pub ref_matches: RefMatchLayers,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nearest_ref_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nearest_ref_distance: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub within_batch_component_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub within_batch_representative_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub subject_id: Option<String>,
    pub split_group_id: String,
    pub thresholds_used: ResolvedThresholds,
}

/// Exact funnel accounting.
///
/// The sequential exclusion stages satisfy:
/// `C4_ref_metadata + EXACT_ASSEMBLY + EXACT_KSET + D3_mash_near_dup +
/// D4_within_batch_dup + UNRESOLVED_OVERLAP + certified_clean == total_queries`.
///
/// `D1_ref_assembly` and `D2_ref_biosample_sra` report metadata match counts
/// across queries; `CLOSE_RELATIVE`, `SHARED_BIOPROJECT`, `SHARED_PMID`, and
/// `ACCESSION_GAP` are non-exclusive tag/diagnostic counts.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[allow(non_snake_case)]
pub struct FunnelCounts {
    pub total_queries: usize,
    pub C4_ref_metadata: usize,
    pub D1_ref_assembly: usize,
    pub D2_ref_biosample_sra: usize,
    pub EXACT_ASSEMBLY: usize,
    pub EXACT_KSET: usize,
    pub D3_mash_near_dup: usize,
    pub CLOSE_RELATIVE: usize,
    pub D4_within_batch_dup: usize,
    pub UNRESOLVED_OVERLAP: usize,
    pub SHARED_BIOPROJECT: usize,
    pub SHARED_PMID: usize,
    pub ACCESSION_GAP: usize,
    pub certified_clean: usize,
}

/// Non-exclusive marginal hit counts for every independent identity and tag
/// layer across all queries (regardless of sequential stage order).
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[allow(non_snake_case)]
pub struct LayerCounts {
    pub D1_ref_assembly: usize,
    pub D2_ref_biosample_sra: usize,
    pub STRAIN_ref_same_species: usize,
    pub STRAIN_ref_only: usize,
    pub EXACT_ASSEMBLY: usize,
    pub EXACT_KSET: usize,
    pub D3_mash_near_dup: usize,
    pub CLOSE_RELATIVE: usize,
    pub SHARED_BIOPROJECT: usize,
    pub SHARED_PMID: usize,
    pub SHARED_SUBJECT: usize,
}

/// Diagnostic and study tags counted strictly within `certified_clean`.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[allow(non_snake_case)]
pub struct CleanTagCounts {
    pub SHARED_BIOPROJECT: usize,
    pub SHARED_PMID: usize,
    pub CLOSE_RELATIVE: usize,
    pub ACCESSION_GAP: usize,
    pub with_subject_group: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ReferenceSummary {
    pub path: String,
    pub source_label: String,
    pub index_sha256: String,
    pub record_count: usize,
    pub accession_gap_count: usize,
    pub accession_gap_share: f64,
    pub with_sketch_k21: usize,
    pub with_kset_digest: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SpeciesScreenSummary {
    pub funnel: FunnelCounts,
    pub layer_counts: LayerCounts,
    pub certified_clean_tags: CleanTagCounts,
    pub thresholds_used: ResolvedThresholds,
}

/// One within-batch duplicate connected component (`<out>/components.json`).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct BatchComponent {
    pub component_id: String,
    pub representative_id: String,
    pub members: Vec<String>,
    pub dropped: Vec<String>,
    pub edge_reasons: Vec<String>,
}

/// One longitudinal patient group (`<out>/subject_groups.json`).
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct SubjectGroup {
    pub subject_id: String,
    pub split_group_id: String,
    pub all_members: Vec<String>,
    pub certified_clean_members: Vec<String>,
}

/// Top-level screen summary written to `<out>/summary.json`.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ScreenSummary {
    pub schema: String,
    pub claim_scope: String,
    pub sealed_set_status: String,
    pub k_digest: usize,
    pub k_sketch: usize,
    pub sketch_size: usize,
    pub index_hashes: usize,
    pub references: Vec<ReferenceSummary>,
    pub thresholds: ThresholdPolicy,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub synonym_map_sha256: Option<String>,
    pub synonym_pairs_loaded: usize,
    pub funnel: FunnelCounts,
    pub layer_counts: LayerCounts,
    pub certified_clean_tags: CleanTagCounts,
    pub within_batch_components: usize,
    pub subject_groups: usize,
    pub per_species: BTreeMap<String, SpeciesScreenSummary>,
}

/// Full in-memory result returned by [`screen`].
pub struct ScreenResult {
    pub summary: ScreenSummary,
    pub decisions: Vec<QueryDecision>,
    pub components: Vec<BatchComponent>,
    pub subject_groups: Vec<SubjectGroup>,
}

struct RefEntry {
    source_label: String,
    rec: IndexRecord,
}

struct RefTables {
    entries: Vec<RefEntry>,
    by_assembly: FxHashMap<String, Vec<(usize, String)>>,
    by_sample_sra: FxHashMap<String, Vec<(usize, String)>>,
    by_species_strain: FxHashMap<(String, String), Vec<usize>>,
    by_assembly_sha256: FxHashMap<String, Vec<usize>>,
    by_kset_sha256: FxHashMap<String, Vec<usize>>,
    by_kset_fnv64: FxHashMap<String, Vec<usize>>,
    by_bioproject: FxHashMap<String, Vec<(usize, String)>>,
    by_pmid: FxHashMap<String, Vec<usize>>,
    by_subject: FxHashMap<String, Vec<usize>>,
    sketch_index: SketchIndex,
    any_has_sequence: bool,
}

impl RefTables {
    fn build(
        loaded: Vec<LoadedIndex>,
        synonyms: &SynonymMap,
        k_sketch: usize,
        sketch_size: usize,
        index_hashes: usize,
    ) -> Self {
        let mut entries = Vec::new();
        for idx in loaded {
            let label = idx.manifest.source_label;
            for rec in idx.records {
                entries.push(RefEntry {
                    source_label: label.clone(),
                    rec,
                });
            }
        }
        let mut by_assembly: FxHashMap<String, Vec<(usize, String)>> = FxHashMap::default();
        let mut by_sample_sra: FxHashMap<String, Vec<(usize, String)>> = FxHashMap::default();
        let mut by_species_strain: FxHashMap<(String, String), Vec<usize>> = FxHashMap::default();
        let mut by_assembly_sha256: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        let mut by_kset_sha256: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        let mut by_kset_fnv64: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        let mut by_bioproject: FxHashMap<String, Vec<(usize, String)>> = FxHashMap::default();
        let mut by_pmid: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        let mut by_subject: FxHashMap<String, Vec<usize>> = FxHashMap::default();
        let mut sketches = Vec::with_capacity(entries.len());
        let mut any_has_sequence = false;

        for (i, entry) in entries.iter().enumerate() {
            let rec = &entry.rec;
            if rec.has_sequence_fingerprint() {
                any_has_sequence = true;
            }
            for raw in &rec.assembly_accessions {
                let canon = synonyms.resolve_assembly(raw).to_string();
                by_assembly.entry(canon).or_default().push((i, raw.clone()));
            }
            for raw in rec
                .biosample_accessions
                .iter()
                .chain(rec.sra_accessions.iter())
            {
                let canon = synonyms.resolve_sample(raw).to_string();
                by_sample_sra
                    .entry(canon)
                    .or_default()
                    .push((i, raw.clone()));
            }
            if let (Some(sp), Some(strain)) = (&rec.species_slug, &rec.strain) {
                by_species_strain
                    .entry((sp.clone(), strain.clone()))
                    .or_default()
                    .push(i);
            }
            if let Some(ref h) = rec.assembly_sha256 {
                by_assembly_sha256.entry(h.clone()).or_default().push(i);
            }
            if let Some(ref h) = rec.sha256_k31 {
                by_kset_sha256.entry(h.clone()).or_default().push(i);
            }
            if let Some(ref h) = rec.fnv64_k31 {
                by_kset_fnv64.entry(h.clone()).or_default().push(i);
            }
            for raw in &rec.bioproject_accessions {
                let canon = synonyms.resolve_study(raw).to_string();
                by_bioproject
                    .entry(canon)
                    .or_default()
                    .push((i, raw.clone()));
            }
            for pmid in &rec.pmids {
                by_pmid.entry(pmid.clone()).or_default().push(i);
            }
            if let Some(ref subj) = rec.subject_id {
                by_subject.entry(subj.clone()).or_default().push(i);
            }
            sketches.push(rec.sketch_k21.clone().unwrap_or_default());
        }

        let sketch_index = SketchIndex::build(sketches, k_sketch, sketch_size, index_hashes);
        RefTables {
            entries,
            by_assembly,
            by_sample_sra,
            by_species_strain,
            by_assembly_sha256,
            by_kset_sha256,
            by_kset_fnv64,
            by_bioproject,
            by_pmid,
            by_subject,
            sketch_index,
            any_has_sequence,
        }
    }

    fn match_query(
        &self,
        q: &IndexRecord,
        synonyms: &SynonymMap,
        thresholds: &ResolvedThresholds,
        scratch: &mut QueryScratch,
    ) -> RefMatchLayers {
        let mut out = RefMatchLayers::default();

        // D1: assembly accession matches
        let mut seen_d1: FxHashSet<(usize, String)> = FxHashSet::default();
        for raw in &q.assembly_accessions {
            let canon = synonyms.resolve_assembly(raw);
            if let Some(list) = self.by_assembly.get(canon) {
                for (ri, _) in list {
                    if seen_d1.insert((*ri, canon.to_string())) {
                        out.d1_assembly.push(RefKeyHit {
                            ref_id: self.entries[*ri].rec.genome_id.clone(),
                            source_label: self.entries[*ri].source_label.clone(),
                            matched_key: raw.clone(),
                            canonical_key: (canon != raw).then(|| canon.to_string()),
                        });
                    }
                }
            }
        }
        sort_key_hits(&mut out.d1_assembly);

        // D2: BioSample / SRA / INSDC alias matches
        let mut seen_d2: FxHashSet<(usize, String)> = FxHashSet::default();
        for raw in q.biosample_accessions.iter().chain(q.sra_accessions.iter()) {
            let canon = synonyms.resolve_sample(raw);
            if let Some(list) = self.by_sample_sra.get(canon) {
                for (ri, _) in list {
                    if seen_d2.insert((*ri, canon.to_string())) {
                        out.d2_biosample_sra.push(RefKeyHit {
                            ref_id: self.entries[*ri].rec.genome_id.clone(),
                            source_label: self.entries[*ri].source_label.clone(),
                            matched_key: raw.clone(),
                            canonical_key: (canon != raw).then(|| canon.to_string()),
                        });
                    }
                }
            }
        }
        sort_key_hits(&mut out.d2_biosample_sra);

        // STRAIN: same normalized species + normalized strain
        if let (Some(sp), Some(strain)) = (&q.species_slug, &q.strain)
            && let Some(list) = self.by_species_strain.get(&(sp.clone(), strain.clone()))
        {
            for &ri in list {
                out.strain_same_species.push(RefKeyHit {
                    ref_id: self.entries[ri].rec.genome_id.clone(),
                    source_label: self.entries[ri].source_label.clone(),
                    matched_key: strain.clone(),
                    canonical_key: None,
                });
            }
        }
        sort_key_hits(&mut out.strain_same_species);

        // EXACT_ASSEMBLY: normalized assembly content SHA-256
        if let Some(ref h) = q.assembly_sha256
            && let Some(list) = self.by_assembly_sha256.get(h)
        {
            for &ri in list {
                out.exact_assembly.push(RefKeyHit {
                    ref_id: self.entries[ri].rec.genome_id.clone(),
                    source_label: self.entries[ri].source_label.clone(),
                    matched_key: h.clone(),
                    canonical_key: Some("assembly_sha256".to_string()),
                });
            }
        }
        sort_key_hits(&mut out.exact_assembly);

        // EXACT_KSET: prefer sha256_k31 when both carry it; fall back to
        // fnv64_k31 (+ n_kmers check when available) when one side came from
        // counted.json without SHA-256.
        let mut seen_kset: FxHashSet<usize> = FxHashSet::default();
        if let Some(ref sha) = q.sha256_k31
            && let Some(list) = self.by_kset_sha256.get(sha)
        {
            for &ri in list {
                if seen_kset.insert(ri) {
                    out.exact_kset.push(RefKeyHit {
                        ref_id: self.entries[ri].rec.genome_id.clone(),
                        source_label: self.entries[ri].source_label.clone(),
                        matched_key: sha.clone(),
                        canonical_key: Some("sha256_k31".to_string()),
                    });
                }
            }
        }
        if let Some(ref fnv) = q.fnv64_k31
            && let Some(list) = self.by_kset_fnv64.get(fnv)
        {
            for &ri in list {
                if seen_kset.contains(&ri) {
                    continue;
                }
                let ref_rec = &self.entries[ri].rec;
                // If both carry sha256_k31 and they disagreed, treat as an FNV
                // collision and do NOT call EXACT_KSET.
                if let (Some(qa), Some(ra)) = (&q.sha256_k31, &ref_rec.sha256_k31)
                    && qa != ra
                {
                    continue;
                }
                if let (Some(qn), Some(rn)) = (q.n_kmers, ref_rec.n_kmers)
                    && qn != rn
                {
                    continue;
                }
                if let (Some(qs), Some(rs)) = (&q.sketch_k31, &ref_rec.sketch_k31)
                    && qs != rs
                {
                    continue;
                }
                if seen_kset.insert(ri) {
                    out.exact_kset.push(RefKeyHit {
                        ref_id: ref_rec.genome_id.clone(),
                        source_label: self.entries[ri].source_label.clone(),
                        matched_key: fnv.clone(),
                        canonical_key: Some("fnv64_k31".to_string()),
                    });
                }
            }
        }
        sort_key_hits(&mut out.exact_kset);

        // Genomic Mash distance: NEAR_DUP vs CLOSE_RELATIVE
        if let Some(ref sketch) = q.sketch_k21
            && !sketch.is_empty()
        {
            let hits = self.sketch_index.query_with_scratch(
                sketch,
                thresholds.close_relative_threshold,
                scratch,
            );
            for (ri, dist) in hits {
                let entry = &self.entries[ri];
                let hit = RefDistanceHit {
                    ref_id: entry.rec.genome_id.clone(),
                    source_label: entry.source_label.clone(),
                    mash_distance: dist,
                    ref_missing_accessions: entry.rec.missing_accessions(),
                };
                if dist <= thresholds.near_dup_threshold {
                    out.near_dup.push(hit);
                } else {
                    out.close_relative.push(hit);
                }
            }
        }

        // Study-overlap tags (never conflated with isolate identity)
        let mut seen_bp: FxHashSet<(usize, String)> = FxHashSet::default();
        for raw in &q.bioproject_accessions {
            let canon = synonyms.resolve_study(raw);
            if let Some(list) = self.by_bioproject.get(canon) {
                for (ri, _) in list {
                    if seen_bp.insert((*ri, canon.to_string())) {
                        out.shared_bioproject.push(RefKeyHit {
                            ref_id: self.entries[*ri].rec.genome_id.clone(),
                            source_label: self.entries[*ri].source_label.clone(),
                            matched_key: raw.clone(),
                            canonical_key: (canon != raw).then(|| canon.to_string()),
                        });
                    }
                }
            }
        }
        sort_key_hits(&mut out.shared_bioproject);

        let mut seen_pmid: FxHashSet<(usize, String)> = FxHashSet::default();
        for pmid in &q.pmids {
            if let Some(list) = self.by_pmid.get(pmid) {
                for &ri in list {
                    if seen_pmid.insert((ri, pmid.clone())) {
                        out.shared_pmid.push(RefKeyHit {
                            ref_id: self.entries[ri].rec.genome_id.clone(),
                            source_label: self.entries[ri].source_label.clone(),
                            matched_key: pmid.clone(),
                            canonical_key: None,
                        });
                    }
                }
            }
        }
        sort_key_hits(&mut out.shared_pmid);

        // Longitudinal patient overlap against reference
        if let Some(ref subj) = q.subject_id
            && let Some(list) = self.by_subject.get(subj)
        {
            for &ri in list {
                out.shared_subject.push(RefKeyHit {
                    ref_id: self.entries[ri].rec.genome_id.clone(),
                    source_label: self.entries[ri].source_label.clone(),
                    matched_key: subj.clone(),
                    canonical_key: None,
                });
            }
        }
        sort_key_hits(&mut out.shared_subject);

        out
    }
}

fn sort_key_hits(hits: &mut [RefKeyHit]) {
    hits.sort_by(|a, b| {
        compare_genome_ids(&a.ref_id, &b.ref_id)
            .then_with(|| a.source_label.cmp(&b.source_label))
            .then_with(|| a.matched_key.cmp(&b.matched_key))
    });
}

/// Screens a query batch against one or more reference indices, performs
/// within-batch connected-component deduplication and longitudinal patient
/// grouping, and writes `summary.json`, `decisions.jsonl`, `components.json`,
/// and `subject_groups.json` under `cfg.out_dir`.
pub fn screen(cfg: &ScreenConfig) -> Result<ScreenResult> {
    ensure!(
        (MIN_K..=MAX_K).contains(&cfg.k_digest),
        "k_digest must be within {MIN_K}..={MAX_K}, got {}",
        cfg.k_digest
    );
    ensure!(
        (MIN_K..=MAX_K).contains(&cfg.k_sketch),
        "k_sketch must be within {MIN_K}..={MAX_K}, got {}",
        cfg.k_sketch
    );

    let loaded_refs: Vec<LoadedIndex> = cfg
        .references
        .iter()
        .map(|p| load_index(p))
        .collect::<Result<_>>()?;
    let (k_digest, k_sketch, sketch_size) = if let Some(first) = loaded_refs.first() {
        for other in &loaded_refs[1..] {
            ensure!(
                other.manifest.k_digest == first.manifest.k_digest
                    && other.manifest.k_sketch == first.manifest.k_sketch
                    && other.manifest.sketch_size == first.manifest.sketch_size,
                "reference index {} (k_digest={}, k_sketch={}, sketch_size={}) disagrees with {} (k_digest={}, k_sketch={}, sketch_size={})",
                other.path.display(),
                other.manifest.k_digest,
                other.manifest.k_sketch,
                other.manifest.sketch_size,
                first.path.display(),
                first.manifest.k_digest,
                first.manifest.k_sketch,
                first.manifest.sketch_size,
            );
        }
        (
            first.manifest.k_digest,
            first.manifest.k_sketch,
            first.manifest.sketch_size,
        )
    } else {
        (cfg.k_digest, cfg.k_sketch, cfg.sketch_size)
    };

    let (synonyms, synonym_map_sha256) = SynonymMap::load(&cfg.synonym_maps)?;
    let thresholds = ThresholdPolicy::load(
        cfg.thresholds.as_deref(),
        cfg.fallback_near_dup,
        cfg.fallback_close_relative,
    )?;

    let ref_summaries: Vec<ReferenceSummary> = loaded_refs
        .iter()
        .map(|idx| {
            let n = idx.manifest.counts.records;
            let gap = idx.manifest.counts.accession_gap;
            ReferenceSummary {
                path: idx.path.display().to_string(),
                source_label: idx.manifest.source_label.clone(),
                index_sha256: idx.manifest.index_sha256.clone(),
                record_count: n,
                accession_gap_count: gap,
                accession_gap_share: if n == 0 { 0.0 } else { gap as f64 / n as f64 },
                with_sketch_k21: idx.manifest.counts.with_sketch_k21,
                with_kset_digest: idx
                    .records
                    .iter()
                    .filter(|r| r.sha256_k31.is_some() || r.fnv64_k31.is_some())
                    .count(),
            }
        })
        .collect();
    let ref_tables = RefTables::build(
        loaded_refs,
        &synonyms,
        k_sketch,
        sketch_size,
        cfg.index_hashes,
    );

    // Load and deduplicate query rows by genome_id, then sort deterministically.
    let mut query_by_id: FxHashMap<String, MetadataRow> = FxHashMap::default();
    for row in read_metadata_jsonl(&cfg.query)? {
        query_by_id
            .entry(row.genome_id.clone())
            .and_modify(|existing| merge_metadata_row(existing, row.clone()))
            .or_insert(row);
    }
    let mut query_rows: Vec<MetadataRow> = query_by_id.into_values().collect();
    query_rows.sort_by(|a, b| compare_genome_ids(&a.genome_id, &b.genome_id));

    // Fingerprint query assemblies in parallel.
    let queries: Vec<IndexRecord> = query_rows
        .into_par_iter()
        .map_init(
            || (Vec::<u64>::new(), Vec::<u64>::new()),
            |(values, scratch), mut row| -> Result<IndexRecord> {
                let fasta = resolve_fasta_for(
                    &row.genome_id,
                    row.fasta_path.as_deref(),
                    cfg.query_fasta_dir.as_deref(),
                )?;
                if let Some(path) = fasta {
                    let fp = kset::fingerprint_fasta_into(
                        &path,
                        k_digest,
                        k_sketch,
                        sketch_size,
                        values,
                        scratch,
                    )
                    .with_context(|| format!("fingerprint query {}", path.display()))?;
                    row.fnv64_k31 = Some(fp.fnv64_k31);
                    row.sha256_k31 = Some(fp.sha256_k31);
                    row.assembly_sha256 = Some(fp.assembly_sha256);
                    row.n_kmers = Some(fp.n_kmers);
                    row.contigs = Some(fp.contigs);
                    row.bp = Some(fp.bp);
                    row.sketch_k21 = Some(fp.sketch_k21);
                    row.sketch_k31 = fp.sketch_k31;
                } else if row.sketch_k21.is_none()
                    && let Some(dir) = cfg.query_sketches_dir.as_deref()
                    && let Some(sk) =
                        load_species_sketch(dir, &row.genome_id, k_sketch, sketch_size)?
                {
                    if row.bp.is_none() {
                        row.bp = Some(sk.length);
                    }
                    row.sketch_k21 = Some(sk.sketch);
                }
                Ok(IndexRecord {
                    genome_id: row.genome_id,
                    species: row.species,
                    species_slug: row.species_slug,
                    assembly_accessions: row.assembly_accessions,
                    biosample_accessions: row.biosample_accessions,
                    sra_accessions: row.sra_accessions,
                    bioproject_accessions: row.bioproject_accessions,
                    pmids: row.pmids,
                    strain: row.strain,
                    subject_id: row.subject_id,
                    fnv64_k31: row.fnv64_k31,
                    sha256_k31: row.sha256_k31,
                    assembly_sha256: row.assembly_sha256,
                    n_kmers: row.n_kmers,
                    contigs: row.contigs,
                    bp: row.bp,
                    sketch_k21: row.sketch_k21,
                    sketch_k31: row.sketch_k31,
                })
            },
        )
        .collect::<Result<Vec<_>>>()?;

    let resolved_thresholds: Vec<ResolvedThresholds> = queries
        .iter()
        .map(|q| thresholds.for_species(q.species_slug.as_deref()))
        .collect();

    // Match every query across all reference layers in parallel.
    let ref_matches_all: Vec<RefMatchLayers> = queries
        .par_iter()
        .enumerate()
        .map_init(
            || QueryScratch::new(ref_tables.entries.len()),
            |scratch, (i, q)| {
                ref_tables.match_query(q, &synonyms, &resolved_thresholds[i], scratch)
            },
        )
        .collect();

    // Stage 1..4: reference exclusion precedence.
    let mut dispositions: Vec<Option<Disposition>> = vec![None; queries.len()];
    let mut survivors: Vec<usize> = Vec::new();
    for (i, m) in ref_matches_all.iter().enumerate() {
        let has_ref_meta = !m.d1_assembly.is_empty()
            || !m.d2_biosample_sra.is_empty()
            || !m.strain_same_species.is_empty();
        if has_ref_meta {
            dispositions[i] = Some(Disposition::ExcludedRefMetadata);
        } else if !m.exact_assembly.is_empty() {
            dispositions[i] = Some(Disposition::ExcludedExactAssembly);
        } else if !m.exact_kset.is_empty() {
            dispositions[i] = Some(Disposition::ExcludedExactKset);
        } else if !m.near_dup.is_empty() {
            dispositions[i] = Some(Disposition::ExcludedMashNearDup);
        } else {
            survivors.push(i);
        }
    }

    // Stage 5: within-batch connected components over survivors (Union-Find
    // across accessions, exact assembly, exact k-mer digest, and near_dup_threshold).
    let (components, comp_assignment) = build_within_batch_components(
        &queries,
        &survivors,
        &synonyms,
        &resolved_thresholds,
        k_sketch,
        sketch_size,
        cfg.index_hashes,
    );
    for (qi, (_, rep_id)) in &comp_assignment {
        if &queries[*qi].genome_id != rep_id {
            dispositions[*qi] = Some(Disposition::ExcludedWithinBatchDup);
        }
    }

    // Stage 6: unresolved overlap gating + clean certification.
    let has_references = !ref_tables.entries.is_empty();
    let mut unresolved_reasons_all: Vec<Vec<String>> = vec![Vec::new(); queries.len()];
    for i in 0..queries.len() {
        if dispositions[i].is_some() {
            continue;
        }
        let q = &queries[i];
        let m = &ref_matches_all[i];
        let reasons = &mut unresolved_reasons_all[i];
        if !q.has_sequence_fingerprint() {
            reasons.push("missing_sequence_fingerprint".to_string());
        }
        if !m.shared_subject.is_empty() {
            reasons.push("shared_subject_with_reference".to_string());
        }
        if !m.close_relative.is_empty() {
            if q.missing_accessions() {
                reasons.push("missing_accessions_with_close_relative".to_string());
            }
            if m.close_relative.iter().any(|h| h.ref_missing_accessions) {
                reasons.push("close_relative_to_unaccessioned_reference".to_string());
            }
        }
        if has_references && !ref_tables.any_has_sequence && q.missing_accessions() {
            reasons.push("missing_accessions_against_metadata_only_reference".to_string());
        }
        if !reasons.is_empty() {
            dispositions[i] = Some(Disposition::ExcludedUnresolvedOverlap);
        } else {
            dispositions[i] = Some(Disposition::CertifiedClean);
        }
    }

    // Assemble per-query decisions and exact funnel/layer/species counts.
    let mut decisions = Vec::with_capacity(queries.len());
    let mut funnel = FunnelCounts {
        total_queries: queries.len(),
        ..Default::default()
    };
    let mut layer_counts = LayerCounts::default();
    let mut clean_tags = CleanTagCounts::default();
    let mut per_species: BTreeMap<String, SpeciesScreenSummary> = BTreeMap::new();
    let mut subject_map: BTreeMap<String, (Vec<String>, Vec<String>)> = BTreeMap::new();

    for i in 0..queries.len() {
        let q = &queries[i];
        let m = &ref_matches_all[i];
        let disp = dispositions[i].unwrap();
        let clean = disp == Disposition::CertifiedClean;
        let t_used = resolved_thresholds[i].clone();

        let mut matched_layers = Vec::new();
        if !m.d1_assembly.is_empty() {
            matched_layers.push("D1_ref_assembly".to_string());
        }
        if !m.d2_biosample_sra.is_empty() {
            matched_layers.push("D2_ref_biosample_sra".to_string());
        }
        if !m.strain_same_species.is_empty() {
            matched_layers.push("STRAIN".to_string());
        }
        if !m.exact_assembly.is_empty() {
            matched_layers.push("EXACT_ASSEMBLY".to_string());
        }
        if !m.exact_kset.is_empty() {
            matched_layers.push("EXACT_KSET".to_string());
        }
        if !m.near_dup.is_empty() {
            matched_layers.push("NEAR_DUP".to_string());
        }
        if !m.close_relative.is_empty() {
            matched_layers.push("CLOSE_RELATIVE".to_string());
        }
        if disp == Disposition::ExcludedWithinBatchDup {
            matched_layers.push("D4_within_batch_dup".to_string());
        }

        let mut study_tags = Vec::new();
        if !m.shared_bioproject.is_empty() {
            study_tags.push("SHARED_BIOPROJECT".to_string());
        }
        if !m.shared_pmid.is_empty() {
            study_tags.push("SHARED_PMID".to_string());
        }

        let mut diagnostic_tags = Vec::new();
        if !m.close_relative.is_empty() {
            diagnostic_tags.push("CLOSE_RELATIVE".to_string());
        }
        if q.missing_accessions() {
            diagnostic_tags.push("ACCESSION_GAP".to_string());
        }
        if !m.shared_subject.is_empty() {
            diagnostic_tags.push("SHARED_SUBJECT".to_string());
        }

        let nearest = m.near_dup.first().or_else(|| m.close_relative.first());
        let (comp_id, rep_id) = match comp_assignment.get(&i) {
            Some((cid, rid)) => (Some(cid.clone()), Some(rid.clone())),
            None => (None, None),
        };
        let split_group_id = match &q.subject_id {
            Some(subj) => format!("subject:{subj}"),
            None => format!("isolate:{}", q.genome_id),
        };
        if let Some(ref subj) = q.subject_id {
            let entry = subject_map.entry(subj.clone()).or_default();
            entry.0.push(q.genome_id.clone());
            if clean {
                entry.1.push(q.genome_id.clone());
            }
        }

        // Update global and per-species counters identically.
        let sp_key = q
            .species_slug
            .clone()
            .unwrap_or_else(|| "unknown".to_string());
        let sp_summary = per_species
            .entry(sp_key)
            .or_insert_with(|| SpeciesScreenSummary {
                funnel: FunnelCounts::default(),
                layer_counts: LayerCounts::default(),
                certified_clean_tags: CleanTagCounts::default(),
                thresholds_used: t_used.clone(),
            });
        sp_summary.funnel.total_queries += 1;

        for (f, lc, ct) in [
            (&mut funnel, &mut layer_counts, &mut clean_tags),
            (
                &mut sp_summary.funnel,
                &mut sp_summary.layer_counts,
                &mut sp_summary.certified_clean_tags,
            ),
        ] {
            tally_counts(f, lc, ct, disp, q, m);
        }

        decisions.push(QueryDecision {
            genome_id: q.genome_id.clone(),
            species: q.species.clone(),
            species_slug: q.species_slug.clone(),
            disposition: disp,
            certified_clean: clean,
            matched_layers,
            study_tags,
            diagnostic_tags,
            unresolved_reasons: unresolved_reasons_all[i].clone(),
            ref_matches: m.clone(),
            nearest_ref_id: nearest.map(|h| h.ref_id.clone()),
            nearest_ref_distance: nearest.map(|h| h.mash_distance),
            within_batch_component_id: comp_id,
            within_batch_representative_id: rep_id,
            subject_id: q.subject_id.clone(),
            split_group_id,
            thresholds_used: t_used,
        });
    }

    let subject_groups: Vec<SubjectGroup> = subject_map
        .into_iter()
        .map(
            |(subject_id, (all_members, certified_clean_members))| SubjectGroup {
                split_group_id: format!("subject:{subject_id}"),
                subject_id,
                all_members,
                certified_clean_members,
            },
        )
        .collect();

    let summary = ScreenSummary {
        schema: SCREEN_SCHEMA.to_string(),
        claim_scope: CLAIM_SCOPE.to_string(),
        sealed_set_status: SEALED_SET_STATUS.to_string(),
        k_digest,
        k_sketch,
        sketch_size,
        index_hashes: cfg.index_hashes,
        references: ref_summaries,
        thresholds,
        synonym_map_sha256,
        synonym_pairs_loaded: synonyms.pair_count(),
        funnel,
        layer_counts,
        certified_clean_tags: clean_tags,
        within_batch_components: components.len(),
        subject_groups: subject_groups.len(),
        per_species,
    };

    // Write outputs atomically under cfg.out_dir.
    std::fs::create_dir_all(&cfg.out_dir)
        .with_context(|| format!("create screen out dir {}", cfg.out_dir.display()))?;
    let mut decision_bytes = Vec::new();
    for d in &decisions {
        serde_json::to_writer(&mut decision_bytes, d)?;
        decision_bytes.push(b'\n');
    }
    write_atomic(&cfg.out_dir.join("decisions.jsonl"), &decision_bytes)?;
    write_atomic(
        &cfg.out_dir.join("components.json"),
        &serde_json::to_vec_pretty(&components)?,
    )?;
    write_atomic(
        &cfg.out_dir.join("subject_groups.json"),
        &serde_json::to_vec_pretty(&subject_groups)?,
    )?;
    write_atomic(
        &cfg.out_dir.join("summary.json"),
        &serde_json::to_vec_pretty(&summary)?,
    )?;

    Ok(ScreenResult {
        summary,
        decisions,
        components,
        subject_groups,
    })
}

fn tally_counts(
    funnel: &mut FunnelCounts,
    layers: &mut LayerCounts,
    clean_tags: &mut CleanTagCounts,
    disp: Disposition,
    q: &IndexRecord,
    m: &RefMatchLayers,
) {
    match disp {
        Disposition::ExcludedRefMetadata => funnel.C4_ref_metadata += 1,
        Disposition::ExcludedExactAssembly => funnel.EXACT_ASSEMBLY += 1,
        Disposition::ExcludedExactKset => funnel.EXACT_KSET += 1,
        Disposition::ExcludedMashNearDup => funnel.D3_mash_near_dup += 1,
        Disposition::ExcludedWithinBatchDup => funnel.D4_within_batch_dup += 1,
        Disposition::ExcludedUnresolvedOverlap => funnel.UNRESOLVED_OVERLAP += 1,
        Disposition::CertifiedClean => funnel.certified_clean += 1,
    }
    let has_d1 = !m.d1_assembly.is_empty();
    let has_d2 = !m.d2_biosample_sra.is_empty();
    let has_strain = !m.strain_same_species.is_empty();
    let has_close = !m.close_relative.is_empty();
    let has_bp = !m.shared_bioproject.is_empty();
    let has_pmid = !m.shared_pmid.is_empty();
    let has_gap = q.missing_accessions();

    funnel.D1_ref_assembly += usize::from(has_d1);
    funnel.D2_ref_biosample_sra += usize::from(has_d2);
    funnel.CLOSE_RELATIVE += usize::from(has_close);
    funnel.SHARED_BIOPROJECT += usize::from(has_bp);
    funnel.SHARED_PMID += usize::from(has_pmid);
    funnel.ACCESSION_GAP += usize::from(has_gap);

    layers.D1_ref_assembly += usize::from(has_d1);
    layers.D2_ref_biosample_sra += usize::from(has_d2);
    layers.STRAIN_ref_same_species += usize::from(has_strain);
    layers.STRAIN_ref_only += usize::from(has_strain && !has_d1 && !has_d2);
    layers.EXACT_ASSEMBLY += usize::from(!m.exact_assembly.is_empty());
    layers.EXACT_KSET += usize::from(!m.exact_kset.is_empty());
    layers.D3_mash_near_dup += usize::from(!m.near_dup.is_empty());
    layers.CLOSE_RELATIVE += usize::from(has_close);
    layers.SHARED_BIOPROJECT += usize::from(has_bp);
    layers.SHARED_PMID += usize::from(has_pmid);
    layers.SHARED_SUBJECT += usize::from(!m.shared_subject.is_empty());

    if disp == Disposition::CertifiedClean {
        clean_tags.SHARED_BIOPROJECT += usize::from(has_bp);
        clean_tags.SHARED_PMID += usize::from(has_pmid);
        clean_tags.CLOSE_RELATIVE += usize::from(has_close);
        clean_tags.ACCESSION_GAP += usize::from(has_gap);
        clean_tags.with_subject_group += usize::from(q.subject_id.is_some());
    }
}

/// Union-Find connected-component deduplication over surviving queries across:
/// - normalized/alias-resolved assembly accessions,
/// - normalized/alias-resolved BioSample and SRA accessions,
/// - exact normalized assembly SHA-256,
/// - exact canonical k-mer digest (`sha256_k31` or verified `fnv64_k31`),
/// - genomic Mash distance `<= min(near_dup_a, near_dup_b)` via [`SketchIndex`].
fn build_within_batch_components(
    queries: &[IndexRecord],
    survivors: &[usize],
    synonyms: &SynonymMap,
    thresholds: &[ResolvedThresholds],
    k_sketch: usize,
    sketch_size: usize,
    index_hashes: usize,
) -> (Vec<BatchComponent>, FxHashMap<usize, (String, String)>) {
    let n = survivors.len();
    if n < 2 {
        return (Vec::new(), FxHashMap::default());
    }
    let mut uf = UnionFind::new(n);
    let mut edge_reasons_by_pair: Vec<(usize, usize, &'static str)> = Vec::new();

    let mut link_bucket = |map: FxHashMap<String, Vec<usize>>, reason: &'static str| {
        for group in map.values() {
            if group.len() < 2 {
                continue;
            }
            let first = group[0];
            for &other in &group[1..] {
                uf.union(first, other);
                edge_reasons_by_pair.push((first, other, reason));
            }
        }
    };

    let mut by_asm: FxHashMap<String, Vec<usize>> = FxHashMap::default();
    let mut by_sample: FxHashMap<String, Vec<usize>> = FxHashMap::default();
    let mut by_asm_sha: FxHashMap<String, Vec<usize>> = FxHashMap::default();
    let mut by_kset_sha: FxHashMap<String, Vec<usize>> = FxHashMap::default();
    let mut by_kset_fnv: FxHashMap<String, Vec<usize>> = FxHashMap::default();

    for (pos, &qi) in survivors.iter().enumerate() {
        let q = &queries[qi];
        for raw in &q.assembly_accessions {
            by_asm
                .entry(synonyms.resolve_assembly(raw).to_string())
                .or_default()
                .push(pos);
        }
        for raw in q.biosample_accessions.iter().chain(q.sra_accessions.iter()) {
            by_sample
                .entry(synonyms.resolve_sample(raw).to_string())
                .or_default()
                .push(pos);
        }
        if let Some(ref h) = q.assembly_sha256 {
            by_asm_sha.entry(h.clone()).or_default().push(pos);
        }
        if let Some(ref h) = q.sha256_k31 {
            by_kset_sha.entry(h.clone()).or_default().push(pos);
        }
        if let Some(ref h) = q.fnv64_k31 {
            by_kset_fnv.entry(h.clone()).or_default().push(pos);
        }
    }

    link_bucket(by_asm, "assembly_accession");
    link_bucket(by_sample, "biosample_or_sra");
    link_bucket(by_asm_sha, "exact_assembly");
    link_bucket(by_kset_sha, "exact_kset_sha256");

    // FNV-only fallback when at least one side lacks sha256_k31, verifying
    // no sha256_k31 disagreement and matching n_kmers / sketch_k31.
    for group in by_kset_fnv.values() {
        if group.len() < 2 {
            continue;
        }
        for a_idx in 0..group.len() {
            let pa = group[a_idx];
            let qa = &queries[survivors[pa]];
            for &pb in &group[a_idx + 1..] {
                let qb = &queries[survivors[pb]];
                if let (Some(sa), Some(sb)) = (&qa.sha256_k31, &qb.sha256_k31)
                    && sa != sb
                {
                    continue;
                }
                if let (Some(na), Some(nb)) = (qa.n_kmers, qb.n_kmers)
                    && na != nb
                {
                    continue;
                }
                if let (Some(ska), Some(skb)) = (&qa.sketch_k31, &qb.sketch_k31)
                    && ska != skb
                {
                    continue;
                }
                uf.union(pa, pb);
                edge_reasons_by_pair.push((pa, pb, "exact_kset_fnv64"));
            }
        }
    }

    // Fast within-batch genomic near-duplicate links via SketchIndex.
    let sketches: Vec<Vec<u64>> = survivors
        .iter()
        .map(|&qi| queries[qi].sketch_k21.clone().unwrap_or_default())
        .collect();
    if sketches.iter().any(|s| !s.is_empty()) {
        let batch_index = SketchIndex::build(sketches, k_sketch, sketch_size, index_hashes);
        let near_edges: Vec<(usize, usize)> = (0..n)
            .into_par_iter()
            .map_init(
                || QueryScratch::new(n),
                |scratch, pos| {
                    let qi = survivors[pos];
                    let Some(ref sk) = queries[qi].sketch_k21 else {
                        return Vec::new();
                    };
                    if sk.is_empty() {
                        return Vec::new();
                    }
                    let t_a = thresholds[qi].near_dup_threshold;
                    batch_index
                        .query_with_scratch(sk, t_a, scratch)
                        .into_iter()
                        .filter_map(|(other_pos, dist)| {
                            if other_pos <= pos {
                                return None;
                            }
                            let t_b = thresholds[survivors[other_pos]].near_dup_threshold;
                            (dist <= t_a.min(t_b)).then_some((pos, other_pos))
                        })
                        .collect::<Vec<_>>()
                },
            )
            .flatten()
            .collect();
        for (a, b) in near_edges {
            uf.union(a, b);
            edge_reasons_by_pair.push((a, b, "mash_near_dup"));
        }
    }

    // Group survivors by Union-Find root.
    let mut by_root: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for (pos, &qi) in survivors.iter().enumerate() {
        let r = uf.find(pos);
        by_root.entry(r).or_default().push(qi);
    }
    let mut reasons_by_root: FxHashMap<usize, BTreeSet<&'static str>> = FxHashMap::default();
    for (a, _b, reason) in edge_reasons_by_pair {
        let r = uf.find(a);
        reasons_by_root.entry(r).or_default().insert(reason);
    }

    let mut raw_components: Vec<(Vec<usize>, Vec<String>)> = Vec::new();
    for (root, mut member_indices) in by_root {
        if member_indices.len() < 2 {
            continue;
        }
        member_indices
            .sort_by(|&a, &b| compare_genome_ids(&queries[a].genome_id, &queries[b].genome_id));
        let reasons: Vec<String> = reasons_by_root
            .remove(&root)
            .unwrap_or_default()
            .into_iter()
            .map(String::from)
            .collect();
        raw_components.push((member_indices, reasons));
    }
    raw_components
        .sort_by(|a, b| compare_genome_ids(&queries[a.0[0]].genome_id, &queries[b.0[0]].genome_id));

    let mut components = Vec::with_capacity(raw_components.len());
    let mut assignment: FxHashMap<usize, (String, String)> = FxHashMap::default();
    for (idx, (member_indices, edge_reasons)) in raw_components.into_iter().enumerate() {
        let component_id = format!("C{:06}", idx + 1);
        let rep_id = queries[member_indices[0]].genome_id.clone();
        let members: Vec<String> = member_indices
            .iter()
            .map(|&qi| queries[qi].genome_id.clone())
            .collect();
        let dropped: Vec<String> = members[1..].to_vec();
        for &qi in &member_indices {
            assignment.insert(qi, (component_id.clone(), rep_id.clone()));
        }
        components.push(BatchComponent {
            component_id,
            representative_id: rep_id,
            members,
            dropped,
            edge_reasons,
        });
    }
    (components, assignment)
}

fn hex_bytes(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn write_atomic(path: &Path, bytes: &[u8]) -> Result<()> {
    ensure_parent(path)?;
    let tmp = path.with_extension(format!("part.{}", std::process::id()));
    std::fs::write(&tmp, bytes).with_context(|| format!("write {}", tmp.display()))?;
    std::fs::rename(&tmp, path)
        .with_context(|| format!("rename {} to {}", tmp.display(), path.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::KmerSpace;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    fn random_dna(rng: &mut ChaCha8Rng, len: usize) -> Vec<u8> {
        (0..len).map(|_| b"ACGT"[rng.random_range(0..4)]).collect()
    }

    fn mutate_dna(rng: &mut ChaCha8Rng, seq: &[u8], rate: f64) -> Vec<u8> {
        seq.iter()
            .map(|&b| {
                if rng.random_bool(rate) {
                    let alt: Vec<u8> = b"ACGT".iter().copied().filter(|&x| x != b).collect();
                    alt[rng.random_range(0..3)]
                } else {
                    b
                }
            })
            .collect()
    }

    fn write_fasta(dir: &Path, name: &str, contigs: &[&[u8]]) -> PathBuf {
        let path = dir.join(format!("{name}.fna"));
        let mut buf = Vec::new();
        for (i, c) in contigs.iter().enumerate() {
            buf.extend_from_slice(format!(">c{i}\n").as_bytes());
            buf.extend_from_slice(c);
            buf.push(b'\n');
        }
        std::fs::write(&path, buf).unwrap();
        path
    }

    #[test]
    fn accession_normalization_and_genome_id_order_match_r0() {
        assert_eq!(
            normalize_assembly_token("gca_000123456.1").as_deref(),
            Some("ASM:000123456")
        );
        assert_eq!(
            normalize_assembly_token(" GCF_000123456.3 ").as_deref(),
            Some("ASM:000123456")
        );
        assert_eq!(
            normalize_sra_keys(&["err100, srr200; drr300  ERR100".into()]),
            vec!["DRR300", "ERR100", "SRR200"]
        );
        assert_eq!(
            normalize_strain(Some("  k-12\t MG1655  ")).as_deref(),
            Some("K-12 MG1655")
        );
        assert_eq!(normalize_strain(Some("N/A")), None);

        // Transitive R0 genome id order on mixed numeric and non-numeric ids.
        let mut ids = vec![
            "SAMN002",
            "562.100",
            "562.20",
            "10.2",
            "SAMN001",
            "2.10",
            "2.10_extra",
        ];
        ids.sort_by(|a, b| compare_genome_ids(a, b));
        assert_eq!(
            ids,
            vec![
                "2.10",
                "10.2",
                "562.20",
                "562.100",
                "2.10_extra",
                "SAMN001",
                "SAMN002",
            ]
        );
    }

    #[test]
    fn synonym_map_resolves_transitive_insdc_ena_ncbi_aliases() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("aliases.tsv");
        std::fs::write(
            &tsv,
            "kind\taccession_a\taccession_b\n\
             sample\tSAMEA1001\tERS2001\n\
             run_to_sample\tERR3001\tERS2001\n\
             sample\tSAMN4001\tSRS4001\n\
             study\tPRJEB5001\tERP005001\n",
        )
        .unwrap();
        let (map, receipt) = SynonymMap::load(&[tsv]).unwrap();
        assert!(receipt.is_some());
        assert_eq!(map.pair_count(), 4);
        assert_eq!(
            map.resolve_sample("ERR3001"),
            map.resolve_sample("SAMEA1001")
        );
        assert_eq!(
            map.resolve_sample("ERS2001"),
            map.resolve_sample("SAMEA1001")
        );
        assert_eq!(
            map.resolve_sample("SRS4001"),
            map.resolve_sample("SAMN4001")
        );
        assert_eq!(
            map.resolve_study("ERP005001"),
            map.resolve_study("PRJEB5001")
        );
    }

    #[test]
    fn sketch_index_matches_exact_brute_force_pairwise_mash_distances() {
        let mut rng = ChaCha8Rng::seed_from_u64(101);
        let space = KmerSpace::new(21);
        let sketch_of = |seq: &[u8]| -> Vec<u64> {
            let mut vals = Vec::new();
            space.for_each_canonical(seq, |x| vals.push(space.mix(x)));
            vals.sort_unstable();
            vals.dedup();
            vals.truncate(1000);
            vals
        };

        let root_a = random_dna(&mut rng, 100_000);
        let root_b = random_dna(&mut rng, 100_000);
        let tiny = random_dna(&mut rng, 60); // shorter than index_hashes=64
        let refs = vec![
            sketch_of(&root_a),
            sketch_of(&mutate_dna(&mut rng, &root_a, 0.0002)),
            sketch_of(&mutate_dna(&mut rng, &root_a, 0.0008)),
            sketch_of(&mutate_dna(&mut rng, &root_a, 0.004)),
            sketch_of(&root_b),
            sketch_of(&tiny),
        ];
        let queries = vec![
            sketch_of(&root_a),
            sketch_of(&mutate_dna(&mut rng, &root_a, 0.0003)),
            sketch_of(&mutate_dna(&mut rng, &root_b, 0.0005)),
            sketch_of(&tiny),
        ];

        let index = SketchIndex::build(refs.clone(), 21, 1000, 64);
        for d_max in [0.0005, 0.001, 0.01] {
            for q in &queries {
                let got = index.query(q, d_max);
                let mut want: Vec<(usize, f64)> = refs
                    .iter()
                    .enumerate()
                    .filter_map(|(i, r)| {
                        let d = cluster::mash_distance(cluster::jaccard(q, r, 1000), 21);
                        (d <= d_max).then_some((i, d))
                    })
                    .collect();
                want.sort_by(|a, b| a.1.total_cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
                assert_eq!(got, want, "d_max={d_max}");
            }
        }
    }

    #[test]
    fn end_to_end_index_and_screen_enforces_all_layers_and_funnel_gates() {
        let dir = tempfile::tempdir().unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(777);
        let ref_fasta_dir = dir.path().join("ref_fasta");
        let query_fasta_dir = dir.path().join("query_fasta");
        std::fs::create_dir_all(&ref_fasta_dir).unwrap();
        std::fs::create_dir_all(&query_fasta_dir).unwrap();

        // Reference genomes:
        // r1: E. coli with full accessions + BioProject PRJEB100 + PMID 111111 + subject SUBJ_REF
        // r2: E. coli with NO accessions (tests accession_gap + unresolved close relative)
        // r3: M. tuberculosis (tighter species-specific threshold)
        // r4: E. coli unrelated strain with strain "EC-REF-4"
        let seq_r1 = random_dna(&mut rng, 80_000);
        let seq_r2 = random_dna(&mut rng, 80_000);
        let seq_r3 = random_dna(&mut rng, 80_000);
        let seq_r4 = random_dna(&mut rng, 80_000);
        write_fasta(&ref_fasta_dir, "562.1", &[&seq_r1]);
        write_fasta(&ref_fasta_dir, "562.2", &[&seq_r2]);
        write_fasta(&ref_fasta_dir, "1773.1", &[&seq_r3]);
        write_fasta(&ref_fasta_dir, "562.4", &[&seq_r4]);

        // Base genomes.jsonl (left untouched!) + sidecar genome_provenance.jsonl.
        let genomes_jsonl = dir.path().join("genomes.jsonl");
        let genomes_text = concat!(
            "{\"genome_id\":\"562.1\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMEA100\",\"assembly_accession\":\"GCA_000100001.1\"}\n",
            "{\"genome_id\":\"562.2\",\"species\":\"Escherichia coli\"}\n",
            "{\"genome_id\":\"1773.1\",\"species\":\"Mycobacterium tuberculosis\",\"biosample_accession\":\"SAMN300\",\"assembly_accession\":\"GCA_000300001.1\"}\n",
            "{\"genome_id\":\"562.4\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN400\"}\n",
        );
        std::fs::write(&genomes_jsonl, genomes_text).unwrap();
        let genomes_sha_before = sha256(genomes_text.as_bytes());

        let prov_jsonl = dir.path().join("genome_provenance.jsonl");
        std::fs::write(
            &prov_jsonl,
            concat!(
                "{\"genome_id\":\"562.1\",\"bioproject_accession\":\"PRJEB100\",\"pmid\":\"111111\",\"subject_id\":\"PAT_REF\"}\n",
                "{\"genome_id\":\"562.4\",\"strain\":\"EC-REF-4\"}\n",
            ),
        )
        .unwrap();

        let ref_index_dir = dir.path().join("ref_index");
        let manifest = build_index(&IndexConfig {
            metadata: vec![genomes_jsonl.clone()],
            provenance: Some(prov_jsonl),
            fasta_dir: Some(ref_fasta_dir),
            source_label: "bvbrc-dev".to_string(),
            out_dir: ref_index_dir.clone(),
            ..Default::default()
        })
        .unwrap();
        assert_eq!(manifest.schema, INDEX_SCHEMA);
        assert_eq!(manifest.sealed_set_status, SEALED_SET_STATUS);
        assert_eq!(manifest.counts.records, 4);
        assert_eq!(manifest.counts.accession_gap, 1);
        assert_eq!(
            sha256(&std::fs::read(&genomes_jsonl).unwrap()),
            genomes_sha_before
        );

        // Synonym map: SAMEA100 <-> ERS100 <-> ERR100, PRJEB100 <-> ERP000100.
        let syn_path = dir.path().join("synonyms.tsv");
        std::fs::write(
            &syn_path,
            "kind\taccession_a\taccession_b\n\
             sample\tSAMEA100\tERS100\n\
             run_to_sample\tERR100\tERS100\n\
             study\tPRJEB100\tERP000100\n",
        )
        .unwrap();

        // Species-specific thresholds: M. tuberculosis tighter than E. coli.
        let thresh_path = dir.path().join("thresholds.json");
        std::fs::write(
            &thresh_path,
            r#"{
                "schema": "grm-dedup-thresholds-v1",
                "default": {"near_dup_threshold": 0.0005, "close_relative_threshold": 0.0015},
                "species": {
                    "escherichia_coli": {"near_dup_threshold": 0.0005, "close_relative_threshold": 0.0020},
                    "mycobacterium_tuberculosis": {"near_dup_threshold": 0.0001, "close_relative_threshold": 0.0005}
                }
            }"#,
        )
        .unwrap();

        // Build query batch exercising every single layer and funnel gate:
        let seq_clean_1 = random_dna(&mut rng, 80_000);
        let seq_clean_2 = random_dna(&mut rng, 80_000);
        let seq_clean_3 = random_dna(&mut rng, 80_000);
        let seq_wb_a = random_dna(&mut rng, 80_000);
        let seq_near_r1 = mutate_dna(&mut rng, &seq_r1, 0.0002); // <= 0.0005 (NEAR_DUP)
        let seq_close_r1 = mutate_dna(&mut rng, &seq_r1, 0.0010); // in (0.0005, 0.0020] (CLOSE_RELATIVE)
        let seq_close_r1_b = mutate_dna(&mut rng, &seq_r1, 0.0010); // distinct from seq_close_r1 (~0.0020 apart)
        let seq_close_r2 = mutate_dna(&mut rng, &seq_r2, 0.0010); // close to unaccessioned r2 -> UNRESOLVED_OVERLAP
        let seq_mtb_near = mutate_dna(&mut rng, &seq_r3, 0.00025); // in (0.0001, 0.0005] for Mtb -> CLOSE_RELATIVE, NOT NEAR_DUP!

        // q_d1: RefSeq GCF matching GenBank GCA of r1 (different version)
        write_fasta(&query_fasta_dir, "q_d1", &[&seq_clean_1]);
        // q_d2_alias: ENA run ERR100 -> ERS100 -> SAMEA100 of r1
        write_fasta(&query_fasta_dir, "q_d2_alias", &[&seq_clean_1]);
        // q_strain: same-species strain EC-REF-4 of r4
        write_fasta(&query_fasta_dir, "q_strain", &[&seq_clean_1]);
        // q_exact_asm: different accession, identical sequence to r1 (lowercase/re-headered)
        write_fasta(&query_fasta_dir, "q_exact_asm", &[&seq_r1]);
        // q_exact_kset: duplicated contig of r1 -> same canonical k-mer set, different assembly_sha256!
        write_fasta(&query_fasta_dir, "q_exact_kset", &[&seq_r1, &seq_r1]);
        // q_near_dup: Mash <= 0.0005 to r1
        write_fasta(&query_fasta_dir, "q_near_dup", &[&seq_near_r1]);
        // q_close_clean: Mash ~0.0010 to r1, HAS accessions -> certified_clean + CLOSE_RELATIVE
        write_fasta(&query_fasta_dir, "q_close_clean", &[&seq_close_r1]);
        // q_close_unres_qgap: Mash ~0.0010 to r1, NO accessions -> UNRESOLVED_OVERLAP
        write_fasta(&query_fasta_dir, "q_close_unres_qgap", &[&seq_close_r1_b]);
        // q_close_unres_rgap: Mash ~0.0010 to unaccessioned r2, HAS accessions -> UNRESOLVED_OVERLAP
        write_fasta(&query_fasta_dir, "q_close_unres_rgap", &[&seq_close_r2]);
        // q_missing_fasta: NO FASTA, unique accession -> UNRESOLVED_OVERLAP
        // q_study_tag: clean genome sharing BioProject ERP000100 (alias of PRJEB100) + PMID 111111
        write_fasta(&query_fasta_dir, "q_study_tag", &[&seq_clean_1]);
        // q_wb_1 and q_wb_2: within-batch duplicates via Mash/sequence (q_wb_1 kept, q_wb_2 dropped D4)
        write_fasta(&query_fasta_dir, "q_wb_1", &[&seq_wb_a]);
        write_fasta(&query_fasta_dir, "q_wb_2", &[&seq_wb_a]);
        // q_pat_1 and q_pat_2: distinct clean serial isolates from same patient SUBJ_EXT
        write_fasta(&query_fasta_dir, "q_pat_1", &[&seq_clean_2]);
        write_fasta(&query_fasta_dir, "q_pat_2", &[&seq_clean_3]);
        // q_mtb_close: M. tuberculosis genome at Mash ~0.00025 (would be NEAR_DUP under 0.0005, but CLOSE_RELATIVE under Mtb 0.0001!)
        write_fasta(&query_fasta_dir, "q_mtb_close", &[&seq_mtb_near]);
        // q_pat_ref: distinct genome from reference patient PAT_REF -> UNRESOLVED_OVERLAP (shared_subject_with_reference)
        let seq_pat_ref = random_dna(&mut rng, 80_000);
        write_fasta(&query_fasta_dir, "q_pat_ref", &[&seq_pat_ref]);

        let query_jsonl = dir.path().join("queries.jsonl");
        std::fs::write(
            &query_jsonl,
            concat!(
                "{\"genome_id\":\"q_d1\",\"species\":\"Escherichia coli\",\"assembly_accession\":\"GCF_000100001.2\"}\n",
                "{\"genome_id\":\"q_d2_alias\",\"species\":\"Escherichia coli\",\"sra_accession\":\"ERR100\"}\n",
                "{\"genome_id\":\"q_strain\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN900\",\"strain\":\"  ec-ref-4 \"}\n",
                "{\"genome_id\":\"q_exact_asm\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN901\"}\n",
                "{\"genome_id\":\"q_exact_kset\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN902\"}\n",
                "{\"genome_id\":\"q_near_dup\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN903\"}\n",
                "{\"genome_id\":\"q_close_clean\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN904\"}\n",
                "{\"genome_id\":\"q_close_unres_qgap\",\"species\":\"Escherichia coli\"}\n",
                "{\"genome_id\":\"q_close_unres_rgap\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN905\"}\n",
                "{\"genome_id\":\"q_missing_fasta\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN906\"}\n",
                "{\"genome_id\":\"q_study_tag\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN907\",\"bioproject_accession\":\"ERP000100\",\"pmid\":\"111111\"}\n",
                "{\"genome_id\":\"q_wb_1\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN908\"}\n",
                "{\"genome_id\":\"q_wb_2\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN909\"}\n",
                "{\"genome_id\":\"q_pat_1\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN910\",\"subject_id\":\"SUBJ_EXT\"}\n",
                "{\"genome_id\":\"q_pat_2\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN911\",\"subject_id\":\"SUBJ_EXT\"}\n",
                "{\"genome_id\":\"q_mtb_close\",\"species\":\"Mycobacterium tuberculosis\",\"biosample_accession\":\"SAMN912\"}\n",
                "{\"genome_id\":\"q_pat_ref\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN913\",\"subject_id\":\"PAT_REF\"}\n",
            ),
        )
        .unwrap();

        let screen_out = dir.path().join("screen_out");
        let res = screen(&ScreenConfig {
            references: vec![ref_index_dir],
            query: query_jsonl,
            query_fasta_dir: Some(query_fasta_dir),
            synonym_maps: vec![syn_path],
            thresholds: Some(thresh_path),
            out_dir: screen_out.clone(),
            ..Default::default()
        })
        .unwrap();

        let by_id: FxHashMap<&str, &QueryDecision> = res
            .decisions
            .iter()
            .map(|d| (d.genome_id.as_str(), d))
            .collect();

        assert_eq!(by_id["q_d1"].disposition, Disposition::ExcludedRefMetadata);
        assert!(
            by_id["q_d1"]
                .matched_layers
                .contains(&"D1_ref_assembly".to_string())
        );

        assert_eq!(
            by_id["q_d2_alias"].disposition,
            Disposition::ExcludedRefMetadata
        );
        assert!(
            by_id["q_d2_alias"]
                .matched_layers
                .contains(&"D2_ref_biosample_sra".to_string())
        );

        assert_eq!(
            by_id["q_strain"].disposition,
            Disposition::ExcludedRefMetadata
        );
        assert!(
            by_id["q_strain"]
                .matched_layers
                .contains(&"STRAIN".to_string())
        );

        assert_eq!(
            by_id["q_exact_asm"].disposition,
            Disposition::ExcludedExactAssembly
        );
        assert!(
            by_id["q_exact_asm"]
                .matched_layers
                .contains(&"EXACT_ASSEMBLY".to_string())
        );
        assert!(
            by_id["q_exact_asm"]
                .matched_layers
                .contains(&"EXACT_KSET".to_string())
        );

        assert_eq!(
            by_id["q_exact_kset"].disposition,
            Disposition::ExcludedExactKset
        );
        assert!(
            !by_id["q_exact_kset"]
                .matched_layers
                .contains(&"EXACT_ASSEMBLY".to_string())
        );
        assert!(
            by_id["q_exact_kset"]
                .matched_layers
                .contains(&"EXACT_KSET".to_string())
        );

        assert_eq!(
            by_id["q_near_dup"].disposition,
            Disposition::ExcludedMashNearDup
        );

        assert_eq!(
            by_id["q_close_clean"].disposition,
            Disposition::CertifiedClean
        );
        assert!(
            by_id["q_close_clean"]
                .diagnostic_tags
                .contains(&"CLOSE_RELATIVE".to_string())
        );

        assert_eq!(
            by_id["q_close_unres_qgap"].disposition,
            Disposition::ExcludedUnresolvedOverlap
        );
        assert!(
            by_id["q_close_unres_qgap"]
                .unresolved_reasons
                .contains(&"missing_accessions_with_close_relative".to_string())
        );

        assert_eq!(
            by_id["q_close_unres_rgap"].disposition,
            Disposition::ExcludedUnresolvedOverlap
        );
        assert!(
            by_id["q_close_unres_rgap"]
                .unresolved_reasons
                .contains(&"close_relative_to_unaccessioned_reference".to_string())
        );

        assert_eq!(
            by_id["q_missing_fasta"].disposition,
            Disposition::ExcludedUnresolvedOverlap
        );
        assert!(
            by_id["q_missing_fasta"]
                .unresolved_reasons
                .contains(&"missing_sequence_fingerprint".to_string())
        );

        assert_eq!(
            by_id["q_study_tag"].disposition,
            Disposition::CertifiedClean
        );
        assert!(
            by_id["q_study_tag"]
                .study_tags
                .contains(&"SHARED_BIOPROJECT".to_string())
        );
        assert!(
            by_id["q_study_tag"]
                .study_tags
                .contains(&"SHARED_PMID".to_string())
        );

        assert_eq!(by_id["q_wb_1"].disposition, Disposition::CertifiedClean);
        assert_eq!(
            by_id["q_wb_2"].disposition,
            Disposition::ExcludedWithinBatchDup
        );
        assert_eq!(
            by_id["q_wb_2"].within_batch_representative_id.as_deref(),
            Some("q_wb_1")
        );

        assert_eq!(by_id["q_pat_1"].disposition, Disposition::CertifiedClean);
        assert_eq!(by_id["q_pat_2"].disposition, Disposition::CertifiedClean);
        assert_eq!(by_id["q_pat_1"].split_group_id, "subject:SUBJ_EXT");
        assert_eq!(by_id["q_pat_2"].split_group_id, "subject:SUBJ_EXT");
        let subj_ext = res
            .subject_groups
            .iter()
            .find(|g| g.subject_id == "SUBJ_EXT")
            .unwrap();
        assert_eq!(subj_ext.certified_clean_members, vec!["q_pat_1", "q_pat_2"]);

        assert_eq!(
            by_id["q_pat_ref"].disposition,
            Disposition::ExcludedUnresolvedOverlap
        );
        assert!(
            by_id["q_pat_ref"]
                .unresolved_reasons
                .contains(&"shared_subject_with_reference".to_string())
        );

        // Species-specific threshold kept Mtb close relative separate from near_dup!
        assert_eq!(
            by_id["q_mtb_close"].disposition,
            Disposition::CertifiedClean
        );
        assert_eq!(by_id["q_mtb_close"].thresholds_used.source, "species");
        assert_eq!(
            by_id["q_mtb_close"].thresholds_used.near_dup_threshold,
            0.0001
        );
        assert!(
            by_id["q_mtb_close"]
                .diagnostic_tags
                .contains(&"CLOSE_RELATIVE".to_string())
        );

        // Funnel partition invariant.
        let f = &res.summary.funnel;
        assert_eq!(f.total_queries, 17);
        assert_eq!(
            f.C4_ref_metadata
                + f.EXACT_ASSEMBLY
                + f.EXACT_KSET
                + f.D3_mash_near_dup
                + f.D4_within_batch_dup
                + f.UNRESOLVED_OVERLAP
                + f.certified_clean,
            f.total_queries
        );
        assert_eq!(f.C4_ref_metadata, 3);
        assert_eq!(f.D1_ref_assembly, 1);
        assert_eq!(f.D2_ref_biosample_sra, 1);
        assert_eq!(f.EXACT_ASSEMBLY, 1);
        assert_eq!(f.EXACT_KSET, 1);
        assert_eq!(f.D3_mash_near_dup, 1);
        assert_eq!(f.D4_within_batch_dup, 1);
        assert_eq!(f.UNRESOLVED_OVERLAP, 4);
        assert_eq!(f.certified_clean, 6);
        assert_eq!(f.SHARED_BIOPROJECT, 1);
        assert_eq!(f.SHARED_PMID, 1);
        assert_eq!(res.summary.claim_scope, CLAIM_SCOPE);
        assert_eq!(res.summary.sealed_set_status, SEALED_SET_STATUS);
        assert_eq!(res.summary.layer_counts.STRAIN_ref_only, 1);
        assert_eq!(res.summary.layer_counts.SHARED_SUBJECT, 1);

        // Verify output files written cleanly.
        for name in [
            "summary.json",
            "decisions.jsonl",
            "components.json",
            "subject_groups.json",
        ] {
            assert!(screen_out.join(name).is_file(), "missing {name}");
        }
    }

    #[test]
    fn sha256_k31_prevents_false_exact_kset_match_on_fnv64_collision() {
        let dir = tempfile::tempdir().unwrap();
        let ref_meta = dir.path().join("ref.jsonl");
        std::fs::write(
            &ref_meta,
            "{\"genome_id\":\"r1\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN1\",\"fnv64_k31\":\"deadbeef01234567\",\"sha256_k31\":\"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\",\"n_kmers\":1000,\"sketch_k21\":[1,2,3]}\n",
        )
        .unwrap();
        let idx_dir = dir.path().join("idx");
        build_index(&IndexConfig {
            metadata: vec![ref_meta],
            out_dir: idx_dir.clone(),
            ..Default::default()
        })
        .unwrap();

        let query_meta = dir.path().join("query.jsonl");
        std::fs::write(
            &query_meta,
            "{\"genome_id\":\"q1\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN2\",\"fnv64_k31\":\"deadbeef01234567\",\"sha256_k31\":\"bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\",\"n_kmers\":1000,\"sketch_k21\":[10,20,30]}\n",
        )
        .unwrap();
        let res = screen(&ScreenConfig {
            references: vec![idx_dir],
            query: query_meta,
            out_dir: dir.path().join("out"),
            ..Default::default()
        })
        .unwrap();
        assert_eq!(res.decisions[0].disposition, Disposition::CertifiedClean);
        assert!(res.decisions[0].ref_matches.exact_kset.is_empty());
    }

    #[test]
    fn load_index_rejects_tampered_records_and_uses_counted_and_sketches() {
        let dir = tempfile::tempdir().unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(55);
        let seq = random_dna(&mut rng, 10_000);
        let fna = write_fasta(dir.path(), "562.99", &[&seq]);
        let fp = kset::fingerprint_fasta(&fna, 31, 21, 1000).unwrap();

        // Build synthetic counted.json and species-check/sketches/562.99.json without FASTA dir.
        let counted_path = dir.path().join("counted.json");
        let counted = CountSummary {
            k: 31,
            genomes: vec![species::Counted {
                genome_id: "562.99".into(),
                n_kmers: fp.n_kmers,
                stats: crate::fasta::SeqStats {
                    length: fp.bp,
                    contigs: fp.contigs,
                    gc: 0.5,
                    n50: fp.bp,
                    ambiguous: 0,
                },
                catalog_version: String::new(),
                kset_digest: Some(fp.fnv64_k31.clone()),
                sketch: fp.sketch_k31.clone(),
                genes: vec![],
            }],
            missing_fasta: vec![],
            failed: vec![],
            pending: vec![],
        };
        std::fs::write(&counted_path, serde_json::to_vec(&counted).unwrap()).unwrap();

        let sketches_dir = dir.path().join("sketches");
        std::fs::create_dir_all(&sketches_dir).unwrap();
        let sk = GenomeSketch {
            genome_id: "562.99".into(),
            k: 21,
            sketch_size: 1000,
            length: fp.bp,
            gc: 0.5,
            sketch: fp.sketch_k21.clone(),
        };
        std::fs::write(
            sketches_dir.join("562.99.json"),
            serde_json::to_vec(&sk).unwrap(),
        )
        .unwrap();

        let meta = dir.path().join("meta.jsonl");
        std::fs::write(
            &meta,
            "{\"genome_id\":\"562.99\",\"species\":\"Escherichia coli\",\"biosample_accession\":\"SAMN555\"}\n",
        )
        .unwrap();

        let idx_dir = dir.path().join("idx");
        build_index(&IndexConfig {
            metadata: vec![meta],
            counted: vec![counted_path],
            sketches_dir: Some(sketches_dir),
            out_dir: idx_dir.clone(),
            ..Default::default()
        })
        .unwrap();

        let loaded = load_index(&idx_dir).unwrap();
        assert_eq!(loaded.records.len(), 1);
        assert_eq!(
            loaded.records[0].fnv64_k31.as_deref(),
            Some(fp.fnv64_k31.as_str())
        );
        assert_eq!(
            loaded.records[0].sketch_k21.as_deref(),
            Some(fp.sketch_k21.as_slice())
        );

        // Tampering records.jsonl must fail SHA-256 verification.
        let mut bytes = std::fs::read(idx_dir.join("records.jsonl")).unwrap();
        bytes.push(b' ');
        std::fs::write(idx_dir.join("records.jsonl"), bytes).unwrap();
        assert!(load_index(&idx_dir).is_err());
    }
}
