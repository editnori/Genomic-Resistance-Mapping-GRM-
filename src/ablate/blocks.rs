//! Feature blocks for ablation: the registry of what each block is and where
//! it comes from, and the adapters that build each block from genotype calls
//! GRM already makes.
//!
//! - C, curated exact calls: the region probes in `probes`, for the species
//!   whose assays can call.
//! - H, hierarchical panel calls: each `panels` panel for the drug, and each
//!   of its targets, over the genome's `catalog::family_calls`.
//! - P, presence: catalog gene groups found by `Catalog::detect`.
//! - N, gene copy counts: unavailable, since nothing in GRM counts copies.
//! - X, genome context: unavailable until a reviewed source exists.
//! - K, anonymous: the k-mer patterns of the species `FeatureMatrix`.
//!
//! When two blocks carry the same column, the lower `priority` keeps it
//! (C=0 < H=1 < P=2 < N=3 < X=4 < K=5); the encoder does that collapse.
//!
//! Every adapter is label-free by signature. It is made from genotype calls
//! alone, and `BlockSource::build` sees only sample ids, never phenotypes.
//!
//! Source proxies. Lineage and cluster ids, BioProject, study, submitter or
//! centre, country, collection and insertion dates, sequencing platform and
//! assembly QC say where a genome came from, not what it carries, and a model
//! can learn resistance prevalence from them. `proxy_kind` names them, and
//! `admit` refuses any block whose id, version or column names match one.
//! They can only ever be reported as a named diagnostic.
//!
//! Splitting groups. Task G's `GroupIds`, and any other group id, decide which
//! genomes are held out together, so they must never also be a predictor.
//! This module gives them no path to a predictor block:
//! - `admit` is the only way to make a `PredictorBlock`, and its sealed bound
//!   accepts only this module's adapters, so a `BlockSource` written
//!   elsewhere, over group ids or anything else, cannot be admitted;
//! - no adapter reads a group id type: each reads one typed genotype source
//!   (`ProbeSet`, `FamilyHit`, `Catalog` or `FeatureMatrix`);
//! - group ids are lineage, cluster, BioProject or study ids, which the
//!   deny-list also refuses by name.
//!
//! That guarantee stops at the input types. `ProbeSet`, `FamilyHit` and
//! `FeatureMatrix` have public fields, so code outside this module can write
//! group membership into one (as a k-mer pattern, say) and have it admitted.
//! Nothing here can detect that deliberate misuse; whoever builds those
//! inputs must build them from genotype calls only.
//!
//! The deny-list is a name check, and a name check has limits. A one-hot
//! column named by its value names no field at all. `proxy_kind` refuses the
//! values with a fixed shape (sequence types, clonal complexes, lineage
//! numbers, bare years), but not place names such as `Germany` or other free
//! values; those are out of scope. The seal on `admit` is the real guard.

use crate::bits::{self, Bits};
use crate::catalog::{self, Catalog, FamilyHit};
use crate::dna::KmerSpace;
use crate::matrix::FeatureMatrix;
use crate::panels::{self, Panel};
use crate::probes::{self, Assay, ProbeCall, ProbeSet};
use anyhow::{Context, Result, bail, ensure};
use serde::{Deserialize, Serialize};
use std::collections::hash_map::Entry;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};

// ------------------------------------------------------- reviewed interface

pub type BlockId = String; // "C","H","P","N","X","K"

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum BlockKind {
    CuratedExact,
    Hierarchical,
    Presence,
    Count,
    Context,
    Anonymous,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BlockSpec {
    pub id: BlockId,
    pub kind: BlockKind,
    pub version: String,
    pub source_sha256: String,
    pub priority: u8, // collapse priority C=0 < H=1 < P=2 < N=3 < X=4 < K=5
    pub screen_cap: Option<usize>, // K: Some(4000)
    pub catalog_date: Option<String>,
    pub circularity_risk: bool,
    pub unavailable_reason: Option<String>, // Some(..) = block unavailable
}

#[derive(Clone)]
pub enum RawColumn {
    Binary {
        present: Bits,
        assessable: Bits,
    },
    Count {
        values: Vec<Option<u32>>,
    },
    Categorical {
        codes: Vec<Option<u16>>,
        levels: Vec<String>,
    },
}

/// Constructors from plain per-row values, for callers that hold columns
/// as data rather than `Bits` (the synthetic interchange's bool arrays and
/// integer codes). Each checks what `check_block` would.
impl RawColumn {
    /// A binary column: row i is present when `present[i]`, and assessable
    /// when `assessable[i]`. Present rows must be assessable.
    pub fn binary(present: &[bool], assessable: &[bool]) -> Result<RawColumn> {
        ensure!(
            present.len() == assessable.len(),
            "{} present calls for {} assessable calls",
            present.len(),
            assessable.len()
        );
        ensure!(
            present.iter().zip(assessable).all(|(&p, &a)| a || !p),
            "a row is present where it is not assessable"
        );
        let bits = |v: &[bool]| {
            Bits::from_indices(v.len(), v.iter().enumerate().filter(|e| *e.1).map(|e| e.0))
        };
        Ok(RawColumn::Binary {
            present: bits(present),
            assessable: bits(assessable),
        })
    }

    /// A binary column assessable in every row.
    pub fn binary_assessable(present: &[bool]) -> RawColumn {
        RawColumn::Binary {
            present: Bits::from_indices(
                present.len(),
                present.iter().enumerate().filter(|e| *e.1).map(|e| e.0),
            ),
            assessable: Bits::ones(present.len()),
        }
    }

    /// A count column; `None` is not assessable.
    pub fn count(values: Vec<Option<u32>>) -> RawColumn {
        RawColumn::Count { values }
    }

    /// A categorical column of codes into `levels`; `None` is not
    /// assessable. Every code must name a level.
    pub fn categorical(codes: Vec<Option<u16>>, levels: Vec<String>) -> Result<RawColumn> {
        ensure!(
            codes.iter().flatten().all(|&c| (c as usize) < levels.len()),
            "a code has no level among {}",
            levels.len()
        );
        Ok(RawColumn::Categorical { codes, levels })
    }
}

#[derive(Clone)]
pub struct RawBlock {
    pub spec: BlockSpec,
    pub sample_ids: Vec<String>,
    pub names: Vec<String>,
    pub cols: Vec<RawColumn>,
}

/// Label-free by signature: a source never sees phenotypes.
///
/// `build` returns one row per requested sample it can build, in the order
/// asked, and leaves out every sample it cannot: a genome missing from the
/// source's input, such as one with no catalog hits computed. The encoder
/// reads a sample with no row as "block unavailable for that sample", so a
/// source must never fill an empty or default row in its place. A value
/// the source cannot assess for a sample it does have, such as a probe
/// `NoCall`, stays in that sample's row as not assessable. A block that is
/// unavailable altogether sets `spec.unavailable_reason` instead, and its
/// `build` fails.
pub trait BlockSource {
    fn spec(&self) -> &BlockSpec;
    fn build(&self, sample_ids: &[String]) -> anyhow::Result<RawBlock>;
}

// ----------------------------------------------------------------- registry

/// What a block is, apart from the data it is built from. An adapter adds
/// the provenance of its input to make the block's `BlockSpec`.
#[derive(Clone, Copy, Debug)]
pub struct BlockDef {
    pub id: &'static str,
    pub kind: BlockKind,
    pub priority: u8,
    pub screen_cap: Option<usize>,
    pub circularity_risk: bool,
    /// Why `circularity_risk` is set as it is.
    pub circularity_note: &'static str,
    /// Set when GRM has no source for the block at all.
    pub unavailable: Option<&'static str>,
}

/// Patterns the K block may pass to a learner: the evaluation's default
/// `screened_features`. The screen ranks patterns on labels, so it runs
/// inside each training fold, downstream of this module.
pub const K_SCREEN_CAP: usize = 4000;

const N_UNAVAILABLE: &str = "no copy-number source: a kset holds each distinct k-mer once (presence only) and the catalog calls genes by containment, so nothing in GRM counts gene copies";
const X_UNAVAILABLE: &str = "no reviewed genome-context source exists; the context fields in genome metadata are source proxies and are refused";
const NO_SOURCE_CIRCULARITY: &str =
    "no source exists, so the risk is unknown and counted as present";

pub static REGISTRY: [BlockDef; 6] = [
    BlockDef {
        id: "C",
        kind: BlockKind::CuratedExact,
        priority: 0,
        screen_cap: None,
        circularity_risk: true,
        circularity_note: "probe regions follow published assay targets, chosen from resistance associations that may include these genomes",
        unavailable: None,
    },
    BlockDef {
        id: "H",
        kind: BlockKind::Hierarchical,
        priority: 1,
        screen_cap: None,
        circularity_risk: true,
        circularity_note: "panel targets and the AMRFinderPlus catalog were curated from published resistance associations that may include these genomes",
        unavailable: None,
    },
    BlockDef {
        id: "P",
        kind: BlockKind::Presence,
        priority: 2,
        screen_cap: None,
        circularity_risk: true,
        circularity_note: "the AMRFinderPlus catalog was curated from published resistance associations that may include these genomes",
        unavailable: None,
    },
    BlockDef {
        id: "N",
        kind: BlockKind::Count,
        priority: 3,
        screen_cap: None,
        circularity_risk: true,
        circularity_note: NO_SOURCE_CIRCULARITY,
        unavailable: Some(N_UNAVAILABLE),
    },
    BlockDef {
        id: "X",
        kind: BlockKind::Context,
        priority: 4,
        screen_cap: None,
        circularity_risk: true,
        circularity_note: NO_SOURCE_CIRCULARITY,
        unavailable: Some(X_UNAVAILABLE),
    },
    BlockDef {
        id: "K",
        kind: BlockKind::Anonymous,
        priority: 5,
        screen_cap: Some(K_SCREEN_CAP),
        circularity_risk: false,
        circularity_note: "patterns come from the genomes alone, with no phenotype or catalog input",
        unavailable: None,
    },
];

pub fn def(id: &str) -> Option<&'static BlockDef> {
    REGISTRY.iter().find(|d| d.id == id)
}

fn registered(id: &str) -> &'static BlockDef {
    def(id).expect("a registered block id")
}

/// Where a block's input came from.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Provenance {
    pub version: String,
    pub source_sha256: String,
    pub catalog_date: Option<String>,
}

impl Provenance {
    /// The AMRFinderPlus catalog in `dir`, as `grm fetch catalog` left it:
    /// its trimmed `version.txt`, a SHA-256 over the catalog files, and the
    /// release date the version starts with.
    pub fn catalog_dir(dir: &Path) -> Result<Provenance> {
        let mut digest = Digest::new();
        let mut version = String::new();
        for name in catalog::CATALOG_FILES {
            let path = dir.join(name);
            let bytes = std::fs::read(&path).with_context(|| format!("read {}", path.display()))?;
            digest.field(name.as_bytes()).field(&bytes);
            if name == "version.txt" {
                version = String::from_utf8_lossy(&bytes).trim().to_string();
            }
        }
        ensure!(
            !version.is_empty(),
            "{}/version.txt is blank",
            dir.display()
        );
        Ok(Provenance {
            catalog_date: date_prefix(&version),
            version,
            source_sha256: digest.finish(),
        })
    }

    /// For a block with no source: nothing to name or hash.
    fn none() -> Provenance {
        Provenance {
            version: String::new(),
            source_sha256: String::new(),
            catalog_date: None,
        }
    }
}

impl BlockDef {
    /// This block's spec over an input with `provenance`. A block with no
    /// source stays unavailable whatever the provenance.
    pub fn spec(&self, provenance: Provenance) -> BlockSpec {
        BlockSpec {
            id: self.id.to_string(),
            kind: self.kind,
            version: provenance.version,
            source_sha256: provenance.source_sha256,
            priority: self.priority,
            screen_cap: self.screen_cap,
            catalog_date: provenance.catalog_date,
            circularity_risk: self.circularity_risk,
            unavailable_reason: self.unavailable.map(str::to_string),
        }
    }
}

/// Refuses a spec that disagrees with the registry or names a source proxy,
/// and holds the ablation contract's provenance rules (CONTRACT.md §5): an
/// available block has a version and a 64-character lowercase hex
/// `source_sha256`; an unavailable one has a non-empty reason, and may
/// leave both empty.
pub fn check_spec(spec: &BlockSpec) -> Result<()> {
    let d = def(&spec.id).with_context(|| format!("block {} is not registered", spec.id))?;
    ensure!(
        spec.kind == d.kind
            && spec.priority == d.priority
            && spec.screen_cap == d.screen_cap
            && spec.circularity_risk == d.circularity_risk,
        "block {}: kind, priority, screen cap or circularity risk differs from the registry",
        spec.id
    );
    ensure!(
        d.unavailable.is_none() || spec.unavailable_reason.is_some(),
        "block {} has no source and must stay unavailable",
        spec.id
    );
    match &spec.unavailable_reason {
        Some(reason) => ensure!(
            !reason.trim().is_empty(),
            "block {} is unavailable without a reason",
            spec.id
        ),
        None => ensure!(
            !spec.version.trim().is_empty(),
            "block {} has no version",
            spec.id
        ),
    }
    ensure!(
        is_sha256(&spec.source_sha256)
            || (spec.unavailable_reason.is_some() && spec.source_sha256.is_empty()),
        "block {}: source_sha256 {:?} is not a lowercase SHA-256 hex digest",
        spec.id,
        spec.source_sha256
    );
    refuse_proxies(&spec.id, [spec.id.as_str(), spec.version.as_str()])
}

impl BlockSpec {
    /// This spec in the contract's canonical JSON (CONTRACT.md §7).
    pub fn canonical_json(&self) -> Result<String> {
        canonical_json(&serde_json::to_value(self)?)
    }
}

/// The bytes Python's `json.dumps(obj, sort_keys=True, separators=(",",
/// ":"), ensure_ascii=True, allow_nan=False)` writes for `value`
/// (CONTRACT.md §7). Integers only: Python's float repr is not serde's, and
/// no block spec carries a float.
pub fn canonical_json(value: &serde_json::Value) -> Result<String> {
    let mut out = String::new();
    write_canonical(value, &mut out)?;
    Ok(out)
}

/// SHA-256 of `value`'s canonical JSON, over its UTF-8 bytes.
pub fn canonical_sha256(value: &serde_json::Value) -> Result<String> {
    let json = canonical_json(value)?;
    let digest = ring::digest::digest(&ring::digest::SHA256, json.as_bytes());
    Ok(digest.as_ref().iter().map(|b| format!("{b:02x}")).collect())
}

fn write_canonical(value: &serde_json::Value, out: &mut String) -> Result<()> {
    use serde_json::Value;
    match value {
        Value::Null => out.push_str("null"),
        Value::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
        Value::Number(n) => {
            ensure!(
                n.is_i64() || n.is_u64(),
                "canonical JSON here carries integers only, not {n}"
            );
            out.push_str(&n.to_string());
        }
        Value::String(s) => write_canonical_str(s, out),
        Value::Array(items) => {
            out.push('[');
            for (i, item) in items.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                write_canonical(item, out)?;
            }
            out.push(']');
        }
        Value::Object(map) => {
            // String order is code-point order, as Python sorts keys.
            let mut keys: Vec<&String> = map.keys().collect();
            keys.sort();
            out.push('{');
            for (i, key) in keys.into_iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                write_canonical_str(key, out);
                out.push(':');
                write_canonical(&map[key.as_str()], out)?;
            }
            out.push('}');
        }
    }
    Ok(())
}

/// A JSON string as Python's `ensure_ascii` writes it: printable ASCII as
/// is, the short escapes, and every other character as lowercase `\uXXXX`
/// UTF-16 units.
fn write_canonical_str(s: &str, out: &mut String) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{8}' => out.push_str("\\b"),
            '\u{c}' => out.push_str("\\f"),
            ' '..='~' => out.push(c),
            _ => {
                let mut units = [0u16; 2];
                for unit in c.encode_utf16(&mut units).iter() {
                    out.push_str(&format!("\\u{unit:04x}"));
                }
            }
        }
    }
    out.push('"');
}

fn is_sha256(s: &str) -> bool {
    s.len() == 64 && s.bytes().all(|b| matches!(b, b'0'..=b'9' | b'a'..=b'f'))
}

/// The YYYY-MM-DD that `s` starts with, as in AMRFinderPlus versions
/// ("2026-08-07.1").
fn date_prefix(s: &str) -> Option<String> {
    let b = s.as_bytes();
    let digits =
        |r: std::ops::Range<usize>| b.get(r).is_some_and(|d| d.iter().all(u8::is_ascii_digit));
    (digits(0..4)
        && b.get(4) == Some(&b'-')
        && digits(5..7)
        && b.get(7) == Some(&b'-')
        && digits(8..10)
        && !b.get(10).is_some_and(u8::is_ascii_digit))
    .then(|| s[..10].to_string())
}

/// The first word of `s` that starts with a date.
fn find_date(s: &str) -> Option<String> {
    s.split_whitespace().find_map(date_prefix)
}

/// SHA-256 over length-prefixed fields, so field boundaries are part of it.
struct Digest(ring::digest::Context);

impl Digest {
    fn new() -> Self {
        Digest(ring::digest::Context::new(&ring::digest::SHA256))
    }

    fn field(&mut self, bytes: &[u8]) -> &mut Self {
        self.0.update(&(bytes.len() as u64).to_le_bytes());
        self.0.update(bytes);
        self
    }

    fn words(&mut self, words: &[u64]) -> &mut Self {
        self.0.update(&(words.len() as u64).to_le_bytes());
        for chunk in words.chunks(1 << 13) {
            let bytes: Vec<u8> = chunk.iter().flat_map(|w| w.to_le_bytes()).collect();
            self.0.update(&bytes);
        }
        self
    }

    fn finish(self) -> String {
        self.0
            .finish()
            .as_ref()
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect()
    }
}

// ------------------------------------------------------------ source proxies

/// A metadata field that says where a genome came from rather than what it
/// carries. Never a predictor; at most a named diagnostic.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize)]
pub enum ProxyKind {
    Lineage,
    Cluster,
    BioProject,
    Study,
    Submitter,
    Country,
    CollectionDate,
    InsertionDate,
    SequencingPlatform,
    AssemblyQc,
}

impl ProxyKind {
    pub const ALL: [ProxyKind; 10] = [
        ProxyKind::Lineage,
        ProxyKind::Cluster,
        ProxyKind::BioProject,
        ProxyKind::Study,
        ProxyKind::Submitter,
        ProxyKind::Country,
        ProxyKind::CollectionDate,
        ProxyKind::InsertionDate,
        ProxyKind::SequencingPlatform,
        ProxyKind::AssemblyQc,
    ];

    pub fn label(self) -> &'static str {
        match self {
            ProxyKind::Lineage => "lineage id",
            ProxyKind::Cluster => "cluster id",
            ProxyKind::BioProject => "BioProject",
            ProxyKind::Study => "study",
            ProxyKind::Submitter => "submitter or centre",
            ProxyKind::Country => "country",
            ProxyKind::CollectionDate => "collection date",
            ProxyKind::InsertionDate => "insertion date",
            ProxyKind::SequencingPlatform => "sequencing platform",
            ProxyKind::AssemblyQc => "assembly QC",
        }
    }
}

/// How a name matches a proxy. `words` match a whole word of the name (split
/// at punctuation, spaces and lower-to-upper case changes); `compounds`
/// match anywhere in the name with its separators removed; `accessions`
/// match a word made of the prefix and digits only (PRJNA12345).
struct ProxyRule {
    kind: ProxyKind,
    words: &'static [&'static str],
    compounds: &'static [&'static str],
    accessions: &'static [&'static str],
}

/// First match wins, so insertion dates come before the general date words.
const PROXY_RULES: &[ProxyRule] = &[
    ProxyRule {
        kind: ProxyKind::Lineage,
        words: &[
            "clade",
            "clades",
            "mlst",
            "cgmlst",
            "serovar",
            "serotype",
            "spoligotype",
        ],
        compounds: &[
            "lineage",
            "sequencetype",
            "clonalcomplex",
            "subclade",
            "phylogroup",
        ],
        accessions: &[],
    },
    ProxyRule {
        kind: ProxyKind::Cluster,
        words: &[],
        compounds: &["cluster"],
        // NCBI Pathogen Detection SNP clusters (PDS000012345.1).
        accessions: &["pds"],
    },
    ProxyRule {
        kind: ProxyKind::BioProject,
        words: &[],
        compounds: &["bioproject"],
        accessions: &["prjna", "prjeb", "prjdb"],
    },
    ProxyRule {
        kind: ProxyKind::Study,
        words: &["study", "studies"],
        compounds: &["studyid", "studyaccession", "studytitle", "studyname"],
        accessions: &["srp", "erp", "drp"],
    },
    ProxyRule {
        kind: ProxyKind::Submitter,
        words: &[
            "center",
            "centre",
            "centers",
            "centres",
            "institution",
            "institute",
            "laboratory",
        ],
        compounds: &["submitter", "submitting", "centername"],
        accessions: &[],
    },
    ProxyRule {
        kind: ProxyKind::Country,
        words: &["geo", "continent", "location", "latitude", "longitude"],
        compounds: &["country", "countries", "geoloc", "geographic"],
        accessions: &[],
    },
    ProxyRule {
        kind: ProxyKind::InsertionDate,
        words: &["inserted"],
        compounds: &["insertiondate", "dateinserted"],
        accessions: &[],
    },
    ProxyRule {
        kind: ProxyKind::CollectionDate,
        words: &["date", "dates", "year", "years", "collected"],
        compounds: &[
            "collectiondate",
            "collectionyear",
            "isolationdate",
            "isolationyear",
            "samplingdate",
            "sampledate",
            "datetime",
        ],
        accessions: &[],
    },
    ProxyRule {
        kind: ProxyKind::SequencingPlatform,
        words: &["platform", "platforms", "sequencer", "instrument"],
        compounds: &[
            "sequencingplatform",
            "sequencingtechnology",
            "sequencingmethod",
            "platformname",
            "illumina",
            "nanopore",
            "pacbio",
            "miseq",
            "hiseq",
            "nextseq",
            "novaseq",
            "iontorrent",
        ],
        accessions: &[],
    },
    ProxyRule {
        kind: ProxyKind::AssemblyQc,
        words: &[
            "assembly", "contig", "contigs", "coverage", "qc", "quality", "n50", "l50",
        ],
        compounds: &[
            "assemblyqc",
            "checkm",
            "completeness",
            "contamination",
            "gccontent",
            "genomelength",
            "genomesize",
            "numcontigs",
            "coveragedepth",
        ],
        accessions: &[],
    },
];

/// Lowercase words of `name`, split at anything but ASCII letters and
/// digits and where a lowercase letter meets an uppercase one.
fn words(name: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut word = String::new();
    let mut prev_lower = false;
    for c in name.chars() {
        if !c.is_ascii_alphanumeric() {
            if !word.is_empty() {
                out.push(std::mem::take(&mut word));
            }
            prev_lower = false;
            continue;
        }
        if c.is_ascii_uppercase() && prev_lower {
            out.push(std::mem::take(&mut word));
        }
        prev_lower = c.is_ascii_lowercase();
        word.push(c.to_ascii_lowercase());
    }
    if !word.is_empty() {
        out.push(word);
    }
    out
}

/// The source proxy `name` refers to, if any.
pub fn proxy_kind(name: &str) -> Option<ProxyKind> {
    let words = words(name);
    let joined = words.concat();
    let accession = |prefix: &str| {
        words.iter().any(|w| {
            w.strip_prefix(prefix)
                .is_some_and(|d| !d.is_empty() && d.bytes().all(|b| b.is_ascii_digit()))
        })
    };
    PROXY_RULES
        .iter()
        .find(|r| {
            r.words.iter().any(|t| words.iter().any(|w| w == t))
                || r.compounds.iter().any(|c| joined.contains(c))
                || r.accessions.iter().any(|p| accession(p))
        })
        .map(|r| r.kind)
        .or_else(|| value_kind(name))
}

/// The proxy a name's last `:` part is a value of, when that part is only a
/// sequence type, clonal complex or lineage number (`ST131`, `CC17`, `L4`)
/// or a bare 4-digit year.
fn value_kind(name: &str) -> Option<ProxyKind> {
    let value = name.rsplit(':').next().unwrap_or(name).trim();
    let digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    let lineage = ["st", "cc", "l"].iter().any(|p| {
        value
            .get(..p.len())
            .is_some_and(|head| head.eq_ignore_ascii_case(p))
            && digits(&value[p.len()..])
    });
    if lineage {
        Some(ProxyKind::Lineage)
    } else if value.len() == 4 && digits(value) {
        Some(ProxyKind::CollectionDate)
    } else {
        None
    }
}

/// Refuses `names` (a block's id, version or column names) when any one is
/// a source proxy or has a character outside printable ASCII. The ASCII rule
/// stops lookalikes (`c\u{43e}untry` with a Cyrillic o, fullwidth letters)
/// from splitting a denied word.
pub fn refuse_proxies<'a>(block: &str, names: impl IntoIterator<Item = &'a str>) -> Result<()> {
    for name in names {
        ensure!(
            name.bytes().all(|b| (b' '..=b'~').contains(&b)),
            "block {block}: {name:?} has a character outside printable ASCII"
        );
        if let Some(kind) = proxy_kind(name) {
            bail!(
                "block {block}: {name:?} is a {} (a source proxy); it can only be reported as a named diagnostic, never a predictor block",
                kind.label()
            );
        }
    }
    Ok(())
}

// ---------------------------------------------------------------- admission

mod sealed {
    /// Implemented only by this module's adapters; see `admit`.
    pub trait Adapter {}
}

/// A block admitted as predictor input: built by one of this module's
/// adapters and checked against the registry, the proxy deny-list and its
/// own shape.
#[derive(Clone)]
pub struct PredictorBlock {
    raw: RawBlock,
}

impl PredictorBlock {
    pub fn raw(&self) -> &RawBlock {
        &self.raw
    }

    pub fn spec(&self) -> &BlockSpec {
        &self.raw.spec
    }
}

/// Builds `source` over `sample_ids` and admits it as predictor input.
///
/// Only this module's adapters can be admitted. A source written anywhere
/// else, such as one over splitting-group ids, does not compile here:
///
/// ```compile_fail,E0277
/// use grm::ablate::blocks::{BlockSource, BlockSpec, RawBlock, admit};
///
/// struct FromGroupIds(BlockSpec);
///
/// impl BlockSource for FromGroupIds {
///     fn spec(&self) -> &BlockSpec {
///         &self.0
///     }
///     fn build(&self, _: &[String]) -> anyhow::Result<RawBlock> {
///         unimplemented!()
///     }
/// }
///
/// fn try_admit(source: &FromGroupIds) {
///     let _ = admit(source, &[]);
/// }
/// ```
pub fn admit<S: BlockSource + sealed::Adapter>(
    source: &S,
    sample_ids: &[String],
) -> Result<PredictorBlock> {
    let spec = source.spec();
    check_spec(spec)?;
    if let Some(reason) = &spec.unavailable_reason {
        bail!("block {} is unavailable: {reason}", spec.id);
    }
    let raw = source.build(sample_ids)?;
    ensure!(
        raw.spec == *spec,
        "block {}: built under another spec",
        spec.id
    );
    check_block(&raw, sample_ids)?;
    Ok(PredictorBlock { raw })
}

/// Checks a block built for `requested`: its rows are some of the requested
/// samples, in their order and each once; every column has one value per
/// row; names are unique, under the block's id with a non-blank part after
/// it, printable ASCII without a backtick (contract §10), and no source
/// proxy; and a binary column is present only where it is assessable.
/// Spaces stay allowed: §14's no-whitespace rule covers synthetic names only.
pub fn check_block(raw: &RawBlock, requested: &[String]) -> Result<()> {
    let id = &raw.spec.id;
    let mut asked = requested.iter();
    for sample in &raw.sample_ids {
        ensure!(
            asked.any(|r| r == sample),
            "block {id}: row {sample} was not requested, is out of order or appears twice"
        );
    }
    ensure!(
        raw.names.len() == raw.cols.len(),
        "block {id}: {} names for {} columns",
        raw.names.len(),
        raw.cols.len()
    );
    let prefix = format!("{id}:");
    let mut seen = HashSet::new();
    for name in &raw.names {
        ensure!(
            name.starts_with(&prefix),
            "block {id}: column {name:?} is not named under {prefix}"
        );
        ensure!(
            !name[prefix.len()..].trim().is_empty(),
            "block {id}: column {name:?} has an empty name after {prefix}"
        );
        ensure!(
            !name.contains('`'),
            "block {id}: column {name:?} has a backtick"
        );
        ensure!(
            seen.insert(name.as_str()),
            "block {id}: column {name:?} appears twice"
        );
    }
    refuse_proxies(id, raw.names.iter().map(String::as_str))?;
    let n = raw.sample_ids.len();
    for (name, col) in raw.names.iter().zip(&raw.cols) {
        let rows = match col {
            RawColumn::Binary {
                present,
                assessable,
            } => {
                ensure!(
                    present.len() == assessable.len(),
                    "block {id}: column {name:?} has present and assessable calls of different lengths"
                );
                ensure!(
                    present.and_not(assessable).count() == 0,
                    "block {id}: column {name:?} is present where it is not assessable"
                );
                present.len()
            }
            RawColumn::Count { values } => values.len(),
            RawColumn::Categorical { codes, levels } => {
                ensure!(
                    codes.iter().flatten().all(|&c| (c as usize) < levels.len()),
                    "block {id}: column {name:?} has a code with no level"
                );
                codes.len()
            }
        };
        ensure!(
            rows == n,
            "block {id}: column {name:?} has {rows} values for {n} rows"
        );
    }
    Ok(())
}

fn unavailable(spec: &BlockSpec) -> Result<()> {
    if let Some(reason) = &spec.unavailable_reason {
        bail!("block {} is unavailable: {reason}", spec.id);
    }
    Ok(())
}

/// Each genome's position in an adapter's input.
struct Genomes(HashMap<String, usize>);

impl Genomes {
    fn new(genomes: Vec<String>) -> Result<Self> {
        let mut at = HashMap::with_capacity(genomes.len());
        for (i, g) in genomes.into_iter().enumerate() {
            match at.entry(g) {
                Entry::Occupied(e) => bail!("genome {} is listed twice", e.key()),
                Entry::Vacant(e) => {
                    e.insert(i);
                }
            }
        }
        Ok(Genomes(at))
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    /// The requested samples this input holds, in their order, with each
    /// one's position in the input. The rest are left out (see
    /// `BlockSource`).
    fn rows(&self, sample_ids: &[String]) -> Result<(Vec<String>, Vec<usize>)> {
        let mut seen = HashSet::with_capacity(sample_ids.len());
        let (mut ids, mut at) = (Vec::new(), Vec::new());
        for id in sample_ids {
            ensure!(seen.insert(id.as_str()), "sample {id} is requested twice");
            if let Some(&i) = self.0.get(id) {
                ids.push(id.clone());
                at.push(i);
            }
        }
        Ok((ids, at))
    }
}

fn always_assessable(present: Bits) -> RawColumn {
    let n = present.len();
    RawColumn::Binary {
        present,
        assessable: Bits::ones(n),
    }
}

// ------------------------------------------------------------- C: probes

const PROBES_VERSION: &str = "grm-probes";

/// Block C: each probe assay's call per genome. Mutation detected is
/// present, wild type absent, and no call (the genome carries none of the
/// probed k-mers) not assessable. Only species with an assay that can call
/// have the block.
pub struct ProbeSource {
    spec: BlockSpec,
    genomes: Genomes,
    sets: Vec<ProbeSet>,
    skipped: Vec<(String, String)>,
}

impl ProbeSource {
    /// From probe sets already called over `genomes` (`probes::build_probes`
    /// keeps the order it was given).
    pub fn new(species: &str, genomes: Vec<String>, sets: Vec<ProbeSet>) -> Result<Self> {
        Self::with_skipped(species, genomes, sets, Vec::new())
    }

    /// Calls every `species` assay from the genomes' ksets. An assay that
    /// cannot call at all (unsupported, or without a verified wild type) is
    /// skipped with its reason; an unreadable kset is an error.
    pub fn from_ksets(
        species: &str,
        genomes: Vec<String>,
        ksets: &[PathBuf],
        space: &KmerSpace,
    ) -> Result<Self> {
        ensure!(
            genomes.len() == ksets.len(),
            "{} ksets for {} genomes",
            ksets.len(),
            genomes.len()
        );
        let (mut sets, mut skipped) = (Vec::new(), Vec::new());
        for assay in probes::for_species(species) {
            match probes::build_probes(assay, ksets, space) {
                Ok(set) => sets.push(set),
                Err(e) if e.is::<probes::KsetsUnreadable>() => return Err(e),
                Err(e) => skipped.push((assay.id.to_string(), format!("{e:#}"))),
            }
        }
        Self::with_skipped(species, genomes, sets, skipped)
    }

    /// Assays `from_ksets` could not call, with the reason.
    pub fn skipped(&self) -> &[(String, String)] {
        &self.skipped
    }

    fn with_skipped(
        species: &str,
        genomes: Vec<String>,
        mut sets: Vec<ProbeSet>,
        skipped: Vec<(String, String)>,
    ) -> Result<Self> {
        let assays = probes::for_species(species);
        let genomes = Genomes::new(genomes)?;
        let mut seen = HashSet::new();
        for set in &sets {
            ensure!(
                assays.iter().any(|a| a.id == set.assay),
                "assay {} is not a {species} assay",
                set.assay
            );
            ensure!(
                seen.insert(set.assay.clone()),
                "assay {} is given twice",
                set.assay
            );
            ensure!(
                set.calls.len() == genomes.len(),
                "assay {}: {} calls for {} genomes",
                set.assay,
                set.calls.len(),
                genomes.len()
            );
        }
        sets.sort_by_key(|s| assays.iter().position(|a| a.id == s.assay));
        // The block's source is the assays it calls, or, with none, the
        // species' assay definitions.
        let described: Vec<&Assay> = if sets.is_empty() {
            assays.clone()
        } else {
            assays
                .iter()
                .copied()
                .filter(|a| seen.contains(a.id))
                .collect()
        };
        let catalog_date = described
            .iter()
            .flat_map(|a| a.regions)
            .filter_map(|r| r.wild_type)
            .filter_map(|w| find_date(w.source))
            .max();
        let ids: Vec<&str> = described.iter().map(|a| a.id).collect();
        let mut digest = Digest::new();
        digest.field(&serde_json::to_vec(&described)?);
        let mut spec = registered("C").spec(Provenance {
            version: format!("{PROBES_VERSION}:{}", ids.join("+")),
            source_sha256: digest.finish(),
            catalog_date,
        });
        if sets.is_empty() {
            spec.unavailable_reason = Some(if assays.is_empty() {
                format!("no probe assay covers {species}")
            } else if skipped.is_empty() {
                format!("no {species} probe calls were given")
            } else {
                let why: Vec<String> = skipped
                    .iter()
                    .map(|(id, why)| format!("{id}: {why}"))
                    .collect();
                format!("no {species} probe assay can call: {}", why.join("; "))
            });
        }
        Ok(ProbeSource {
            spec,
            genomes,
            sets,
            skipped,
        })
    }
}

impl sealed::Adapter for ProbeSource {}

impl BlockSource for ProbeSource {
    fn spec(&self) -> &BlockSpec {
        &self.spec
    }

    fn build(&self, sample_ids: &[String]) -> Result<RawBlock> {
        unavailable(&self.spec)?;
        let (ids, at) = self.genomes.rows(sample_ids)?;
        let (mut names, mut cols) = (Vec::new(), Vec::new());
        for set in &self.sets {
            let mut present = Bits::zeros(at.len());
            let mut assessable = Bits::zeros(at.len());
            for (i, &g) in at.iter().enumerate() {
                match set.calls[g] {
                    ProbeCall::MutationDetected => {
                        present.insert(i);
                        assessable.insert(i);
                    }
                    ProbeCall::WildType => assessable.insert(i),
                    ProbeCall::NoCall => {}
                }
            }
            names.push(format!("C:{}", set.assay));
            cols.push(RawColumn::Binary {
                present,
                assessable,
            });
        }
        Ok(RawBlock {
            spec: self.spec.clone(),
            sample_ids: ids,
            names,
            cols,
        })
    }
}

// ------------------------------------------------------------- H: panels

const PANELS_VERSION: &str = "panels-v1";

/// Block H: for each panel that targets the drug, whether the panel calls
/// (any target detected) and, under it, whether each target is detected,
/// from the genome's catalog gene calls. Every call is assessable.
pub struct PanelSource {
    spec: BlockSpec,
    genomes: Genomes,
    panels: Vec<&'static Panel>,
    calls: Vec<Vec<FamilyHit>>,
}

impl PanelSource {
    /// `calls` are each genome's `catalog::family_calls`, as the count
    /// sidecars keep them, made under `calls_version`
    /// (`Catalog::calls_version`) of the catalog `provenance` describes.
    /// List only genomes whose calls were computed; the others get no row.
    pub fn new(
        drug: &str,
        provenance: &Provenance,
        calls_version: &str,
        genomes: Vec<String>,
        calls: Vec<Vec<FamilyHit>>,
    ) -> Result<Self> {
        ensure!(
            calls_version.starts_with(&format!("{}|", provenance.version)),
            "gene calls made under {calls_version} are not from catalog {}",
            provenance.version
        );
        ensure!(
            calls.len() == genomes.len(),
            "{} gene call lists for {} genomes",
            calls.len(),
            genomes.len()
        );
        let panels = panels::for_drug(drug);
        let mut digest = Digest::new();
        digest
            .field(provenance.source_sha256.as_bytes())
            .field(calls_version.as_bytes())
            .field(drug.as_bytes())
            .field(&serde_json::to_vec(&panels)?);
        let mut spec = registered("H").spec(Provenance {
            version: format!("{calls_version}|{PANELS_VERSION}:{drug}"),
            source_sha256: digest.finish(),
            catalog_date: provenance.catalog_date.clone(),
        });
        if panels.is_empty() {
            spec.unavailable_reason = Some(format!("no panel targets {drug}"));
        }
        Ok(PanelSource {
            spec,
            genomes: Genomes::new(genomes)?,
            panels,
            calls,
        })
    }

    /// Calls each genome's genes from its sorted k-mers (a kset's values)
    /// as `grm count` does, then builds as `new`.
    pub fn from_catalog(
        drug: &str,
        cat: &Catalog,
        provenance: &Provenance,
        genomes: Vec<String>,
        genome_kmers: impl IntoIterator<Item = Result<Vec<u64>>>,
    ) -> Result<Self> {
        same_catalog(cat, provenance)?;
        let calls = catalog_hits(cat, &genomes, genome_kmers)?
            .iter()
            .map(|hits| catalog::family_calls(cat, hits))
            .collect();
        Self::new(drug, provenance, &cat.calls_version(), genomes, calls)
    }
}

impl sealed::Adapter for PanelSource {}

impl BlockSource for PanelSource {
    fn spec(&self) -> &BlockSpec {
        &self.spec
    }

    fn build(&self, sample_ids: &[String]) -> Result<RawBlock> {
        unavailable(&self.spec)?;
        let (ids, at) = self.genomes.rows(sample_ids)?;
        let n = at.len();
        let (mut names, mut cols) = (Vec::new(), Vec::new());
        for panel in &self.panels {
            let mut any = Bits::zeros(n);
            let mut targets = vec![Bits::zeros(n); panel.targets.len()];
            for (i, &g) in at.iter().enumerate() {
                let hits = &self.calls[g];
                if !panels::detected(panel, hits).is_empty() {
                    any.insert(i);
                }
                for (t, target) in panel.targets.iter().enumerate() {
                    if hits.iter().any(|h| panels::target_hit(target, h)) {
                        targets[t].insert(i);
                    }
                }
            }
            names.push(format!("H:{}", panel.id));
            cols.push(always_assessable(any));
            for (target, found) in panel.targets.iter().zip(targets) {
                names.push(format!("H:{}/{target}", panel.id));
                cols.push(always_assessable(found));
            }
        }
        Ok(RawBlock {
            spec: self.spec.clone(),
            sample_ids: ids,
            names,
            cols,
        })
    }
}

// ----------------------------------------------------------- P: presence

const PRESENCE_VERSION: &str = "group-presence-v1";

/// Block P: one column per catalog gene group, present when any allele of
/// the group reaches the call containment in the genome. Unlike
/// `family_calls`, near-identical alleles of other groups are not
/// suppressed. Every call is assessable.
pub struct PresenceSource {
    spec: BlockSpec,
    genomes: Genomes,
    groups: Vec<String>,
    /// Per genome, the columns present.
    present: Vec<Vec<u32>>,
}

impl PresenceSource {
    /// `genome_kmers` yields each genome's sorted k-mers (a kset's values),
    /// in `genomes` order. List only genomes whose ksets exist; the others
    /// get no row.
    pub fn new(
        cat: &Catalog,
        provenance: &Provenance,
        genomes: Vec<String>,
        genome_kmers: impl IntoIterator<Item = Result<Vec<u64>>>,
    ) -> Result<Self> {
        same_catalog(cat, provenance)?;
        let groups: Vec<String> = cat
            .genes
            .iter()
            .map(|g| g.group.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let column: HashMap<&str, u32> = groups
            .iter()
            .enumerate()
            .map(|(i, g)| (g.as_str(), i as u32))
            .collect();
        let present = catalog_hits(cat, &genomes, genome_kmers)?
            .into_iter()
            .map(|hits| {
                let mut cols: Vec<u32> = hits
                    .iter()
                    .map(|&(g, _)| column[cat.genes[g as usize].group.as_str()])
                    .collect();
                cols.sort_unstable();
                cols.dedup();
                cols
            })
            .collect();
        let mut digest = Digest::new();
        digest
            .field(provenance.source_sha256.as_bytes())
            .field(PRESENCE_VERSION.as_bytes())
            .field(&catalog::CALL_CONTAINMENT.to_le_bytes());
        let spec = registered("P").spec(Provenance {
            version: format!("{}|{PRESENCE_VERSION}", cat.version),
            source_sha256: digest.finish(),
            catalog_date: provenance.catalog_date.clone(),
        });
        Ok(PresenceSource {
            spec,
            genomes: Genomes::new(genomes)?,
            groups,
            present,
        })
    }
}

impl sealed::Adapter for PresenceSource {}

impl BlockSource for PresenceSource {
    fn spec(&self) -> &BlockSpec {
        &self.spec
    }

    fn build(&self, sample_ids: &[String]) -> Result<RawBlock> {
        unavailable(&self.spec)?;
        let (ids, at) = self.genomes.rows(sample_ids)?;
        let mut present = vec![Bits::zeros(at.len()); self.groups.len()];
        for (i, &g) in at.iter().enumerate() {
            for &c in &self.present[g] {
                present[c as usize].insert(i);
            }
        }
        Ok(RawBlock {
            spec: self.spec.clone(),
            sample_ids: ids,
            names: self.groups.iter().map(|g| format!("P:{g}")).collect(),
            cols: present.into_iter().map(always_assessable).collect(),
        })
    }
}

fn same_catalog(cat: &Catalog, provenance: &Provenance) -> Result<()> {
    ensure!(
        cat.version == provenance.version,
        "the provenance describes catalog {} but the loaded catalog is {}",
        provenance.version,
        cat.version
    );
    Ok(())
}

/// Each genome's catalog hits at the call containment.
fn catalog_hits(
    cat: &Catalog,
    genomes: &[String],
    genome_kmers: impl IntoIterator<Item = Result<Vec<u64>>>,
) -> Result<Vec<Vec<(u32, f32)>>> {
    let mut out = Vec::with_capacity(genomes.len());
    for kmers in genome_kmers {
        let kmers = kmers?;
        let genome = genomes
            .get(out.len())
            .context("more k-mer sets than genomes")?;
        ensure!(
            kmers.windows(2).all(|w| w[0] < w[1]),
            "genome {genome}: k-mers must be sorted and distinct, as a kset's are"
        );
        out.push(cat.detect(&kmers, catalog::CALL_CONTAINMENT));
    }
    ensure!(
        out.len() == genomes.len(),
        "{} k-mer sets for {} genomes",
        out.len(),
        genomes.len()
    );
    Ok(out)
}

// -------------------------------------------------------- N, X: no source

/// Block N, gene copy counts: always unavailable (`N_UNAVAILABLE`).
pub struct CountSource {
    spec: BlockSpec,
}

impl CountSource {
    pub fn new() -> Self {
        CountSource {
            spec: registered("N").spec(Provenance::none()),
        }
    }
}

impl Default for CountSource {
    fn default() -> Self {
        Self::new()
    }
}

impl sealed::Adapter for CountSource {}

impl BlockSource for CountSource {
    fn spec(&self) -> &BlockSpec {
        &self.spec
    }

    fn build(&self, _: &[String]) -> Result<RawBlock> {
        unavailable(&self.spec)?;
        bail!("block N has no source")
    }
}

/// Block X, genome context: always unavailable (`X_UNAVAILABLE`).
pub struct ContextSource {
    spec: BlockSpec,
}

impl ContextSource {
    pub fn new() -> Self {
        ContextSource {
            spec: registered("X").spec(Provenance::none()),
        }
    }
}

impl Default for ContextSource {
    fn default() -> Self {
        Self::new()
    }
}

impl sealed::Adapter for ContextSource {}

impl BlockSource for ContextSource {
    fn spec(&self) -> &BlockSpec {
        &self.spec
    }

    fn build(&self, _: &[String]) -> Result<RawBlock> {
        unavailable(&self.spec)?;
        bail!("block X has no source")
    }
}

// ------------------------------------------------------------ K: k-mers

/// Block K: one column per pattern of the species matrix, present in the
/// genomes that carry it. Every pattern is built; the screen to
/// `K_SCREEN_CAP` ranks on labels and so runs inside each training fold.
pub struct KmerSource<'a> {
    spec: BlockSpec,
    matrix: &'a FeatureMatrix,
    genomes: Genomes,
}

impl<'a> KmerSource<'a> {
    /// `matrix` must come from k-mer counting. Its fields are public, so
    /// this checks only its shape, not that `patterns` hold k-mer presence
    /// rather than, say, group membership.
    pub fn new(matrix: &'a FeatureMatrix) -> Result<Self> {
        let (n, p) = (matrix.n_genomes(), &matrix.patterns);
        ensure!(
            p.cols == n && p.words == bits::words_for(n) && p.data.len() == p.rows * p.words,
            "matrix patterns are not {n} genomes wide"
        );
        ensure!(
            matrix.kmer_offsets.len() == p.rows + 1
                && matrix.kmer_offsets.last().map(|&e| e as usize) == Some(matrix.kmers.len()),
            "matrix k-mer offsets do not cover {} patterns and {} k-mers",
            p.rows,
            matrix.kmers.len()
        );
        let mut digest = Digest::new();
        digest
            .field(&serde_json::to_vec(&matrix.summary)?)
            .words(&p.data)
            .words(&matrix.kmer_offsets)
            .words(&matrix.kmers);
        let s = &matrix.summary;
        let spec = registered("K").spec(Provenance {
            version: format!(
                "k={} min_support={} patterns={} genomes={n}",
                s.k, s.min_support, p.rows
            ),
            source_sha256: digest.finish(),
            catalog_date: None,
        });
        Ok(KmerSource {
            spec,
            matrix,
            genomes: Genomes::new(s.genomes.clone())?,
        })
    }

    pub fn column_name(pattern: usize) -> String {
        format!("K:p{pattern}")
    }

    /// The matrix pattern a K column holds, for `FeatureMatrix::kmers_of`.
    pub fn pattern_of(name: &str) -> Option<usize> {
        name.strip_prefix("K:p")?.parse().ok()
    }
}

impl sealed::Adapter for KmerSource<'_> {}

impl BlockSource for KmerSource<'_> {
    fn spec(&self) -> &BlockSpec {
        &self.spec
    }

    fn build(&self, sample_ids: &[String]) -> Result<RawBlock> {
        unavailable(&self.spec)?;
        let (ids, at) = self.genomes.rows(sample_ids)?;
        let rows = self.matrix.patterns.rows;
        let cols = (0..rows)
            .map(|r| {
                let row = self.matrix.patterns.row(r);
                always_assessable(Bits::from_indices(
                    at.len(),
                    at.iter()
                        .enumerate()
                        .filter(|&(_, &c)| bits::get(row, c))
                        .map(|(i, _)| i),
                ))
            })
            .collect();
        Ok(RawBlock {
            spec: self.spec.clone(),
            sample_ids: ids,
            names: (0..rows).map(Self::column_name).collect(),
            cols,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bits::BitMatrix;
    use crate::matrix::MatrixSummary;

    fn ids(names: &[&str]) -> Vec<String> {
        names.iter().map(|s| s.to_string()).collect()
    }

    /// A binary column's present and assessable rows.
    fn binary(block: &RawBlock, name: &str) -> (Vec<usize>, Vec<usize>) {
        let at = block
            .names
            .iter()
            .position(|n| n == name)
            .unwrap_or_else(|| panic!("no column {name}"));
        match &block.cols[at] {
            RawColumn::Binary {
                present,
                assessable,
            } => (
                present.iter_ones().collect(),
                assessable.iter_ones().collect(),
            ),
            _ => panic!("{name} is not binary"),
        }
    }

    fn hit(family: &str, group: &str) -> FamilyHit {
        FamilyHit {
            family: family.into(),
            group: group.into(),
            allele: group.into(),
            containment: 1.0,
            class: "fixture".into(),
            subclass: "fixture".into(),
        }
    }

    fn probe_source() -> ProbeSource {
        let set = ProbeSet {
            assay: "gc_gyra".into(),
            regions: vec![],
            calls: vec![
                ProbeCall::WildType,
                ProbeCall::MutationDetected,
                ProbeCall::NoCall,
            ],
        };
        ProbeSource::new("Neisseria gonorrhoeae", ids(&["g0", "g1", "g2"]), vec![set]).unwrap()
    }

    /// A fixture catalog's provenance (not a real release).
    fn fixture_provenance() -> Provenance {
        Provenance {
            version: "2026-01-02.1".into(),
            source_sha256: "a".repeat(64),
            catalog_date: Some("2026-01-02".into()),
        }
    }

    const CALLS_VERSION: &str = "2026-01-02.1|group-calls-v2";

    fn panel_source(genomes: &[&str], calls: Vec<Vec<FamilyHit>>) -> PanelSource {
        PanelSource::new(
            "meropenem",
            &fixture_provenance(),
            CALLS_VERSION,
            ids(genomes),
            calls,
        )
        .unwrap()
    }

    // Two synthetic genes, in different catalog groups.
    const GENE_A: &str = "ATGAAACGTTACGGCATTGCAGGCTTACCGATAGCTTGACCTAGGCATCGATCGGATCCATGCAAGCTTGG";
    const GENE_B: &str =
        "GGTACCGAGCTCGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCG";
    const NEITHER: &str = "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGG";

    /// A catalog of GENE_A (blaKPC) and GENE_B (tet(X)) with a fixture version.
    fn catalog_fixture() -> tempfile::TempDir {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(
            dir.path().join("ReferenceGeneCatalog.txt"),
            "type\tclass\tsubclass\trefseq_protein_accession\tgenbank_protein_accession\tgene_family\tproduct_name\n\
             AMR\tBETA-LACTAM\tCARBAPENEM\tWP_A\t\tblaKPC\tcarbapenem-hydrolyzing class A beta-lactamase\n\
             AMR\tTETRACYCLINE\tTETRACYCLINE\tWP_B\t\ttet(X)\ttetracycline-inactivating monooxygenase\n",
        )
        .unwrap();
        std::fs::write(
            dir.path().join("AMR_CDS.fa"),
            format!(
                ">WP_A|NG_A|1|1|blaKPC-2|blaKPC|beta-lactamase\n{GENE_A}\n>WP_B|NG_B|1|1|tet(X)|tet(X)|monooxygenase\n{GENE_B}\n"
            ),
        )
        .unwrap();
        std::fs::write(dir.path().join("version.txt"), "2026-01-02.1\n").unwrap();
        dir
    }

    fn kmers_of(seqs: &[&str]) -> Vec<u64> {
        let space = KmerSpace::new(31);
        let mut out = Vec::new();
        for seq in seqs {
            space.for_each_canonical(seq.as_bytes(), |x| out.push(space.mix(x)));
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    fn ok_all(kmers: Vec<Vec<u64>>) -> impl Iterator<Item = Result<Vec<u64>>> {
        kmers.into_iter().map(Ok)
    }

    /// Five genomes, three patterns: p0 in g0 g1, p1 in g2 g4, p2 in g1 g2 g3.
    fn matrix_fixture() -> FeatureMatrix {
        FeatureMatrix {
            summary: MatrixSummary {
                k: 31,
                genomes: ids(&["g0", "g1", "g2", "g3", "g4"]),
                min_support: 1,
                n_patterns: 3,
                kmers_distinct: 6,
                kmers_kept: 6,
                kmers_rare: 0,
                kmers_core: 0,
                kset_digests: vec![],
            },
            patterns: BitMatrix {
                rows: 3,
                cols: 5,
                words: 1,
                data: vec![0b00011, 0b10100, 0b01110],
            },
            support: vec![2, 2, 3],
            kmer_offsets: vec![0, 1, 3, 6],
            kmers: vec![11, 21, 22, 31, 32, 33],
        }
    }

    #[test]
    fn the_registry_lists_six_blocks_in_collapse_order() {
        use BlockKind as B;
        let got: Vec<(&str, BlockKind, u8)> = REGISTRY
            .iter()
            .map(|d| (d.id, d.kind, d.priority))
            .collect();
        assert_eq!(
            got,
            [
                ("C", B::CuratedExact, 0),
                ("H", B::Hierarchical, 1),
                ("P", B::Presence, 2),
                ("N", B::Count, 3),
                ("X", B::Context, 4),
                ("K", B::Anonymous, 5),
            ]
        );
        for d in &REGISTRY {
            assert_eq!(d.screen_cap, (d.id == "K").then_some(K_SCREEN_CAP));
            assert!(!d.circularity_note.is_empty());
        }
        assert_eq!(
            K_SCREEN_CAP,
            crate::learn::eval::EvalConfig::default().screened_features
        );
        let mut spec = probe_source().spec().clone();
        spec.priority = 5;
        assert!(check_spec(&spec).is_err());
        spec = probe_source().spec().clone();
        spec.circularity_risk = false;
        assert!(check_spec(&spec).is_err());
        spec = probe_source().spec().clone();
        spec.id = "G".into();
        assert!(check_spec(&spec).is_err());
    }

    #[test]
    fn the_k_screen_cap_is_the_evaluations_screened_features() {
        // One retained-feature budget. The K block's cap reaches
        // `screen_block` through the admitted block's spec, and `evaluate`
        // screens to `screened_features`. If the two drift apart, a K-only
        // arm keeps a different number of patterns than the pipeline offers
        // boosting, the lasso and the forest.
        let budget = crate::learn::eval::EvalConfig::default().screened_features;
        assert_eq!(K_SCREEN_CAP, budget);
        let matrix = matrix_fixture();
        let block = admit(&KmerSource::new(&matrix).unwrap(), &matrix.summary.genomes).unwrap();
        assert_eq!(block.spec().screen_cap, Some(budget));
    }

    #[test]
    fn the_deny_list_refuses_each_source_proxy() {
        use ProxyKind::*;
        let cases: &[(ProxyKind, &[&str])] = &[
            (
                Lineage,
                &[
                    "lineage",
                    "MTB sublineage",
                    "Clade",
                    "mlst",
                    "sequence_type",
                    "clonalComplex",
                ],
            ),
            (Cluster, &["cluster", "snp_cluster_id", "PDS000123456.1"]),
            (
                BioProject,
                &[
                    "BioProject",
                    "bioproject_accession",
                    "PRJNA12345",
                    "PRJEB678",
                ],
            ),
            (Study, &["study", "Study Title", "SRP012345"]),
            (
                Submitter,
                &[
                    "submitter",
                    "Sequencing Centre",
                    "sequencing_center",
                    "Institution",
                ],
            ),
            (
                Country,
                &[
                    "country",
                    "Isolation Country",
                    "geo_loc_name",
                    "isolationCountry",
                ],
            ),
            (
                CollectionDate,
                &["Collection Date", "collection_year", "collectionDate"],
            ),
            (InsertionDate, &["Date Inserted", "insertion_date"]),
            (
                SequencingPlatform,
                &[
                    "Sequencing Platform",
                    "sequencing_platform",
                    "Illumina MiSeq",
                    "Oxford Nanopore",
                ],
            ),
            (
                AssemblyQc,
                &[
                    "contig_n50",
                    "Contigs",
                    "Genome Length",
                    "GC Content",
                    "CheckM completeness",
                    "coverage",
                    "assembly_qc",
                ],
            ),
        ];
        let covered: BTreeSet<ProxyKind> = cases.iter().map(|c| c.0).collect();
        assert_eq!(covered, ProxyKind::ALL.into_iter().collect());
        let samples = ids(&["g0"]);
        for &(kind, names) in cases {
            for &name in names {
                assert_eq!(proxy_kind(name), Some(kind), "{name}");
                // As a column name.
                let mut block = probe_source().build(&samples).unwrap();
                block.names[0] = format!("C:{name}");
                let err = check_block(&block, &samples).unwrap_err();
                assert!(format!("{err:#}").contains("named diagnostic"), "{err:#}");
                // As the block's source.
                let mut spec = probe_source().spec().clone();
                spec.version = name.to_string();
                let err = check_spec(&spec).unwrap_err();
                assert!(format!("{err:#}").contains("named diagnostic"), "{err:#}");
            }
        }
    }

    /// A block whose columns are named for the provenance side files' fields.
    /// No adapter reads those files, so this one exists only in tests: it
    /// asks whether `admit` would refuse such a block if one did.
    struct SideFileSource {
        spec: BlockSpec,
        names: Vec<String>,
    }

    impl sealed::Adapter for SideFileSource {}

    impl BlockSource for SideFileSource {
        fn spec(&self) -> &BlockSpec {
            &self.spec
        }

        fn build(&self, sample_ids: &[String]) -> Result<RawBlock> {
            let n = sample_ids.len();
            let cols = self
                .names
                .iter()
                .map(|_| RawColumn::binary(&vec![false; n], &vec![true; n]))
                .collect::<Result<_>>()?;
            Ok(RawBlock {
                spec: self.spec.clone(),
                sample_ids: sample_ids.to_vec(),
                names: self.names.clone(),
                cols,
            })
        }
    }

    #[test]
    fn admit_refuses_a_block_of_the_provenance_side_file_fields() {
        // Every field of bvbrc/genome_dates.jsonl and curated/provenance.jsonl
        // but the row keys (`report/provenance/<id>.json` carries the same).
        let mut fields: Vec<String> = [
            serde_json::to_value(crate::bvbrc::GenomeDate::default()).unwrap(),
            serde_json::to_value(crate::curate::PairProvenance::default()).unwrap(),
        ]
        .iter()
        .flat_map(|row| row.as_object().unwrap().keys().cloned())
        .filter(|k| k != "genome_id" && k != "antibiotic")
        .collect();
        fields.sort();
        // Dates and years: source proxies, each refused by name alone.
        let proxies = [
            "completion_date",
            "date_inserted",
            "date_inserted_first",
            "date_inserted_last",
            "testing_standard_years",
        ];
        // Record and measurement details the name check does not catch.
        // Nothing reads them into a block: no adapter reads the side files,
        // and only this module's adapters can be admitted.
        let unnamed = [
            "methods",
            "mic_spread",
            "n_records",
            "pmids",
            "sources",
            "standards",
        ];
        let mut classified: Vec<&str> = proxies.iter().chain(&unnamed).copied().collect();
        classified.sort();
        assert_eq!(fields, classified, "classify each new side-file field here");

        let samples = ids(&["g0", "g1"]);
        let source = |names: &[&str]| SideFileSource {
            spec: probe_source().spec().clone(),
            names: names.iter().map(|n| format!("C:{n}")).collect(),
        };
        let all: Vec<&str> = fields.iter().map(String::as_str).collect();
        let err = admit(&source(&all), &samples).err().unwrap();
        assert!(format!("{err:#}").contains("named diagnostic"), "{err:#}");
        for name in proxies {
            let err = admit(&source(&[name]), &samples).err().unwrap();
            assert!(format!("{err:#}").contains("named diagnostic"), "{err:#}");
            // As the block's version, too.
            let mut versioned = source(&[]);
            versioned.spec.version = name.to_string();
            let err = admit(&versioned, &samples).err().unwrap();
            assert!(format!("{err:#}").contains("named diagnostic"), "{err:#}");
        }
        for name in unnamed {
            assert!(admit(&source(&[name]), &samples).is_ok(), "{name}");
        }
    }

    #[test]
    fn the_deny_list_passes_genotype_names() {
        let mut names: Vec<String> = [
            "K:p17",
            "P:blaOXA (OXA-48 family)",
            "P:aac(6')-Ib-cr",
            "P:tet(M)",
            "IS26 insertion sequence",
        ]
        .map(String::from)
        .to_vec();
        names.extend(probes::ASSAYS.iter().map(|a| format!("C:{}", a.id)));
        for p in panels::PANELS {
            names.push(format!("H:{}", p.id));
            names.extend(p.targets.iter().map(|t| format!("H:{}/{t}", p.id)));
        }
        for name in &names {
            assert_eq!(proxy_kind(name), None, "{name}");
        }
    }

    /// The error `check_block` gives a probe block whose first column is
    /// renamed to `name`.
    fn column_error(name: &str) -> String {
        let samples = ids(&["g0"]);
        let mut block = probe_source().build(&samples).unwrap();
        block.names[0] = name.to_string();
        format!("{:#}", check_block(&block, &samples).unwrap_err())
    }

    #[test]
    fn word_only_terms_and_camel_case_words_are_refused() {
        use ProxyKind::*;
        // Each term is listed only as a word, with no compound behind it,
        // and the camelCase names reach it only through the case split.
        for (name, kind) in [
            ("instrument", SequencingPlatform),
            ("Sequencing Instrument", SequencingPlatform),
            ("sequencingInstrument", SequencingPlatform),
            ("studies", Study),
            ("latitude", Country),
            ("sampleLatitude", Country),
            ("laboratory", Submitter),
            ("hostLaboratory", Submitter),
        ] {
            assert_eq!(proxy_kind(name), Some(kind), "{name}");
            let err = column_error(&format!("C:{name}"));
            assert!(err.contains("named diagnostic"), "{err}");
        }
    }

    #[test]
    fn the_deny_list_refuses_lookalike_run_together_and_value_names() {
        use ProxyKind::*;
        let run_together = [
            ("centername", Submitter),
            ("studyname", Study),
            ("datetime", CollectionDate),
            ("isolationyear", CollectionDate),
            ("sampledate", CollectionDate),
            ("platformname", SequencingPlatform),
            ("numcontigs", AssemblyQc),
            ("coveragedepth", AssemblyQc),
            ("subclade", Lineage),
            ("phylo_group", Lineage),
        ];
        let values = [
            ("ST131", Lineage),
            ("st131", Lineage),
            ("CC17", Lineage),
            ("L4", Lineage),
            ("l2", Lineage),
            ("2019", CollectionDate),
        ];
        for (name, kind) in run_together.into_iter().chain(values) {
            assert_eq!(proxy_kind(name), Some(kind), "{name}");
            let err = column_error(&format!("C:{name}"));
            assert!(err.contains("named diagnostic"), "{err}");
        }
        // Value shapes match only the whole last part of a name.
        for name in [
            "K:p17",
            "C:rrs_1401",
            "P:stx2",
            "P:blaCTX-M-15",
            "C:ST",
            "C:12345",
            "H:mec/ST131x",
        ] {
            assert_eq!(proxy_kind(name), None, "{name}");
        }
        // Lookalikes and control characters fall to the printable-ASCII
        // rule, as column names and as the block's source.
        for name in [
            "c\u{43e}untry",
            "\u{ff43}\u{ff4f}\u{ff55}\u{ff4e}\u{ff54}\u{ff52}\u{ff59}",
            "bio\u{440}roject",
            "d\u{430}te",
            "a\tb",
            "a\nb",
        ] {
            let err = column_error(&format!("C:{name}"));
            assert!(err.contains("printable ASCII"), "{err}");
            let mut spec = probe_source().spec().clone();
            spec.version = name.to_string();
            let err = format!("{:#}", check_spec(&spec).unwrap_err());
            assert!(err.contains("printable ASCII"), "{err}");
        }
    }

    #[test]
    fn check_block_refuses_backticks_and_empty_names() {
        for (name, why) in [
            ("C:a`b", "backtick"),
            ("C:`", "backtick"),
            ("C:", "empty name"),
            ("C:  ", "empty name"),
        ] {
            let err = column_error(name);
            assert!(err.contains(why), "{name}: {err}");
        }
        // Spaces inside a name stay allowed.
        let samples = ids(&["g0"]);
        let mut block = probe_source().build(&samples).unwrap();
        block.names[0] = "C:a b".to_string();
        check_block(&block, &samples).unwrap();
    }

    #[test]
    fn every_available_block_carries_its_provenance() {
        let dir = catalog_fixture();
        let cat = Catalog::load(dir.path(), &KmerSpace::new(31)).unwrap();
        let prov = Provenance::catalog_dir(dir.path()).unwrap();
        assert_eq!(prov.version, "2026-01-02.1");
        assert_eq!(prov.catalog_date.as_deref(), Some("2026-01-02"));
        assert!(is_sha256(&prov.source_sha256));
        let genomes = ids(&["g0", "g1"]);
        let kmers = || ok_all(vec![kmers_of(&[GENE_A]), kmers_of(&[])]);
        let matrix = matrix_fixture();
        let specs = [
            probe_source().spec().clone(),
            PanelSource::from_catalog("meropenem", &cat, &prov, genomes.clone(), kmers())
                .unwrap()
                .spec()
                .clone(),
            PresenceSource::new(&cat, &prov, genomes.clone(), kmers())
                .unwrap()
                .spec()
                .clone(),
            KmerSource::new(&matrix).unwrap().spec().clone(),
        ];
        // C's date is the AMRFinderPlus release of the gc_gyra wild type.
        let dates = [
            Some("2026-08-07"),
            Some("2026-01-02"),
            Some("2026-01-02"),
            None,
        ];
        for (spec, date) in specs.iter().zip(dates) {
            check_spec(spec).unwrap();
            assert!(spec.unavailable_reason.is_none(), "{}", spec.id);
            assert!(!spec.version.is_empty(), "{spec:?}");
            assert!(is_sha256(&spec.source_sha256), "{spec:?}");
            assert_eq!(spec.catalog_date.as_deref(), date, "{}", spec.id);
            assert_eq!(
                spec.circularity_risk,
                def(&spec.id).unwrap().circularity_risk
            );
        }
        let distinct: HashSet<&str> = specs.iter().map(|s| s.source_sha256.as_str()).collect();
        assert_eq!(distinct.len(), specs.len());
        // The digest follows the catalog's bytes.
        std::fs::write(dir.path().join("AMR_CDS.fa"), ">edited\nACGT\n").unwrap();
        assert_ne!(
            Provenance::catalog_dir(dir.path()).unwrap().source_sha256,
            prov.source_sha256
        );
        // An available spec without provenance is refused.
        let mut bare = specs[3].clone();
        bare.source_sha256.clear();
        assert!(check_spec(&bare).is_err());
        let mut bare = specs[3].clone();
        bare.version = " ".into();
        assert!(check_spec(&bare).is_err());
    }

    #[test]
    fn counts_and_context_are_unavailable() {
        fn refused<S: BlockSource + sealed::Adapter>(source: &S, why: &str) {
            let samples = ids(&["g0"]);
            let spec = source.spec();
            let reason = spec.unavailable_reason.as_deref().unwrap();
            assert!(reason.contains(why), "{reason}");
            check_spec(spec).unwrap();
            assert!(source.build(&samples).is_err());
            assert!(admit(source, &samples).is_err());
            let mut forged = spec.clone();
            forged.unavailable_reason = None;
            forged.version = "forged".into();
            forged.source_sha256 = "0".repeat(64);
            let err = check_spec(&forged).unwrap_err();
            assert!(
                format!("{err:#}").contains("must stay unavailable"),
                "{err:#}"
            );
        }
        refused(&CountSource::new(), "copy-number");
        refused(&ContextSource::new(), "context");
    }

    #[test]
    fn specs_hold_the_contract_provenance_rules() {
        let available = KmerSource::new(&matrix_fixture()).unwrap().spec().clone();
        let none = CountSource::new().spec().clone();
        // An unavailable block may leave version and digest empty.
        assert_eq!(
            (none.version.as_str(), none.source_sha256.as_str()),
            ("", "")
        );
        check_spec(&none).unwrap();
        let cases: [(&str, BlockSpec, bool); 6] = [
            ("available", available.clone(), true),
            (
                "uppercase digest",
                BlockSpec {
                    source_sha256: available.source_sha256.to_uppercase(),
                    ..available.clone()
                },
                false,
            ),
            (
                "short digest",
                BlockSpec {
                    source_sha256: available.source_sha256[..63].into(),
                    ..available.clone()
                },
                false,
            ),
            (
                "empty version",
                BlockSpec {
                    version: String::new(),
                    ..available.clone()
                },
                false,
            ),
            (
                "empty reason",
                BlockSpec {
                    unavailable_reason: Some(String::new()),
                    ..none.clone()
                },
                false,
            ),
            (
                "unavailable with a bad digest",
                BlockSpec {
                    source_sha256: "not hex".into(),
                    ..none.clone()
                },
                false,
            ),
        ];
        for (case, spec, ok) in cases {
            assert_eq!(check_spec(&spec).is_ok(), ok, "{case}");
        }
    }

    #[test]
    fn specs_serialise_to_the_contract_canonical_json() {
        let spec = registered("K").spec(Provenance {
            version: "k=31 min_support=1 patterns=3 genomes=5".into(),
            source_sha256: "0123456789abcdef".repeat(4),
            catalog_date: None,
        });
        let json = spec.canonical_json().unwrap();
        assert_eq!(
            json,
            r#"{"catalog_date":null,"circularity_risk":false,"id":"K","kind":"Anonymous","priority":5,"screen_cap":4000,"source_sha256":"0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef","unavailable_reason":null,"version":"k=31 min_support=1 patterns=3 genomes=5"}"#
        );
        assert_eq!(serde_json::from_str::<BlockSpec>(&json).unwrap(), spec);
        // Python's ensure_ascii escapes, surrogate pairs included.
        let spec = BlockSpec {
            catalog_date: Some("2026-01-02".into()),
            unavailable_reason: Some("caf\u{e9} \"q\" \\ \n\t\u{1}\u{7f}\u{1F600}".into()),
            ..CountSource::new().spec().clone()
        };
        let json = spec.canonical_json().unwrap();
        assert_eq!(
            json,
            r#"{"catalog_date":"2026-01-02","circularity_risk":true,"id":"N","kind":"Count","priority":3,"screen_cap":null,"source_sha256":"","unavailable_reason":"caf\u00e9 \"q\" \\ \n\t\u0001\u007f\ud83d\ude00","version":""}"#
        );
        assert!(json.is_ascii());
        assert_eq!(serde_json::from_str::<BlockSpec>(&json).unwrap(), spec);
        // Nested keys are sorted too, and floats are refused.
        let nested = serde_json::json!({"b": [1, {"d": true, "c": null}], "a": -2});
        assert_eq!(
            canonical_json(&nested).unwrap(),
            r#"{"a":-2,"b":[1,{"c":null,"d":true}]}"#
        );
        assert!(canonical_json(&serde_json::json!({"x": 0.5})).is_err());
        assert!(is_sha256(&canonical_sha256(&nested).unwrap()));
    }

    #[test]
    fn a_probe_no_call_is_not_assessable() {
        let block = admit(&probe_source(), &ids(&["g0", "g1", "g2"])).unwrap();
        assert_eq!(block.raw().names, ["C:gc_gyra"]);
        assert_eq!(binary(block.raw(), "C:gc_gyra"), (vec![1], vec![0, 1]));
    }

    #[test]
    fn probes_are_unavailable_where_no_assay_can_call() {
        let space = KmerSpace::new(31);
        let missing = [PathBuf::from("/nonexistent/g0.kset")];
        let g = ids(&["g0"]);
        for (species, why) in [
            ("Staphylococcus aureus", "no probe assay covers"),
            (
                "Mycobacterium tuberculosis",
                "no verified wild-type sequence",
            ),
            ("Escherichia coli", "synonymous"),
        ] {
            let source = ProbeSource::from_ksets(species, g.clone(), &missing, &space).unwrap();
            let reason = source.spec().unavailable_reason.clone().unwrap();
            assert!(reason.contains(why), "{species}: {reason}");
            assert!(admit(&source, &g).is_err());
        }
        // A species whose probes resolve must read its ksets.
        let err = ProbeSource::from_ksets("Neisseria gonorrhoeae", g.clone(), &missing, &space)
            .err()
            .unwrap();
        assert!(err.is::<probes::KsetsUnreadable>(), "{err:#}");
        // Calls must come from the species' own assays, one per genome.
        let set = |n: usize| ProbeSet {
            assay: "gc_gyra".into(),
            regions: vec![],
            calls: vec![ProbeCall::WildType; n],
        };
        assert!(ProbeSource::new("Escherichia coli", g.clone(), vec![set(1)]).is_err());
        assert!(ProbeSource::new("Neisseria gonorrhoeae", g.clone(), vec![set(2)]).is_err());
        assert!(ProbeSource::new("Neisseria gonorrhoeae", g, vec![set(1), set(1)]).is_err());
    }

    #[test]
    fn panel_columns_follow_their_targets() {
        let genomes = ["g0", "g1", "g2"];
        let calls = vec![
            vec![hit("blaKPC", "blaKPC")],
            vec![hit("blaOXA", "blaOXA (OXA-1 family)")],
            vec![],
        ];
        let block = admit(&panel_source(&genomes, calls.clone()), &ids(&genomes)).unwrap();
        let raw = block.raw();
        assert_eq!(binary(raw, "H:carbapenemase"), (vec![0], vec![0, 1, 2]));
        assert_eq!(binary(raw, "H:carbapenemase/blaKPC").0, [0]);
        assert!(
            binary(raw, "H:carbapenemase/blaOXA (OXA-48 family)")
                .0
                .is_empty()
        );
        // Each panel column is the union of its targets.
        let meropenem = panels::for_drug("meropenem");
        for panel in &meropenem {
            let union: BTreeSet<usize> = panel
                .targets
                .iter()
                .flat_map(|t| binary(raw, &format!("H:{}/{t}", panel.id)).0)
                .collect();
            assert_eq!(
                binary(raw, &format!("H:{}", panel.id)).0,
                union.into_iter().collect::<Vec<_>>()
            );
        }
        let columns: usize = meropenem.iter().map(|p| p.targets.len() + 1).sum();
        assert_eq!(raw.names.len(), columns);
        let none = PanelSource::new(
            "isoniazid",
            &fixture_provenance(),
            CALLS_VERSION,
            ids(&genomes),
            calls.clone(),
        )
        .unwrap();
        assert!(
            none.spec()
                .unavailable_reason
                .as_deref()
                .unwrap()
                .contains("no panel")
        );
        // Calls from another catalog release are refused.
        assert!(
            PanelSource::new(
                "meropenem",
                &fixture_provenance(),
                "2025-12-01.1|group-calls-v2",
                ids(&genomes),
                calls
            )
            .is_err()
        );
    }

    #[test]
    fn presence_follows_catalog_detection() {
        let dir = catalog_fixture();
        let cat = Catalog::load(dir.path(), &KmerSpace::new(31)).unwrap();
        let prov = Provenance::catalog_dir(dir.path()).unwrap();
        let genomes = ids(&["g0", "g1", "g2"]);
        let kmers = vec![
            kmers_of(&[GENE_A]),
            kmers_of(&[GENE_A, GENE_B]),
            kmers_of(&[NEITHER]),
        ];
        let source =
            PresenceSource::new(&cat, &prov, genomes.clone(), ok_all(kmers.clone())).unwrap();
        let block = admit(&source, &genomes).unwrap();
        assert_eq!(block.raw().names, ["P:blaKPC", "P:tet(X)"]);
        assert_eq!(binary(block.raw(), "P:blaKPC"), (vec![0, 1], vec![0, 1, 2]));
        assert_eq!(binary(block.raw(), "P:tet(X)").0, [1]);
        // H reads the same catalog through family_calls.
        let panels = PanelSource::from_catalog(
            "meropenem",
            &cat,
            &prov,
            genomes.clone(),
            ok_all(kmers.clone()),
        )
        .unwrap();
        let block = admit(&panels, &genomes).unwrap();
        assert_eq!(binary(block.raw(), "H:carbapenemase/blaKPC").0, [0, 1]);
        // Detection reads sorted k-mers, one set per genome, from the
        // catalog the provenance describes.
        let mut unsorted = kmers.clone();
        unsorted[1].reverse();
        assert!(PresenceSource::new(&cat, &prov, genomes.clone(), ok_all(unsorted)).is_err());
        assert!(
            PresenceSource::new(&cat, &prov, genomes.clone(), ok_all(kmers[..2].to_vec())).is_err()
        );
        let other = Provenance {
            version: "other".into(),
            ..prov.clone()
        };
        assert!(PresenceSource::new(&cat, &other, genomes, ok_all(kmers)).is_err());
    }

    #[test]
    fn kmer_patterns_round_trip() {
        let matrix = matrix_fixture();
        let source = KmerSource::new(&matrix).unwrap();
        let block = admit(&source, &matrix.summary.genomes).unwrap();
        let raw = block.raw();
        assert_eq!(raw.names, ["K:p0", "K:p1", "K:p2"]);
        let mut rebuilt = BitMatrix::with_cols(5);
        for (name, col) in raw.names.iter().zip(&raw.cols) {
            let RawColumn::Binary {
                present,
                assessable,
            } = col
            else {
                panic!("{name} is not binary")
            };
            assert_eq!(assessable.count(), 5);
            let p = KmerSource::pattern_of(name).unwrap();
            assert_eq!(p, rebuilt.rows);
            assert_eq!(present.count(), matrix.support[p]);
            rebuilt.push_row(present.words());
        }
        assert_eq!(rebuilt.data, matrix.patterns.data);
        assert_eq!(
            matrix.kmers_of(KmerSource::pattern_of("K:p1").unwrap()),
            [21, 22]
        );
        // A subset in another order reads each genome's own column.
        let block = admit(&source, &ids(&["g4", "g1"])).unwrap();
        assert_eq!(binary(block.raw(), "K:p0").0, [1]);
        assert_eq!(binary(block.raw(), "K:p1").0, [0]);
        assert_eq!(binary(block.raw(), "K:p2").0, [1]);
        // A matrix whose rows are not its genomes wide is refused.
        let mut narrow = matrix_fixture();
        narrow.patterns.cols = 4;
        assert!(KmerSource::new(&narrow).is_err());
    }

    #[test]
    fn samples_outside_the_input_get_no_row() {
        let dir = catalog_fixture();
        let cat = Catalog::load(dir.path(), &KmerSpace::new(31)).unwrap();
        let prov = Provenance::catalog_dir(dir.path()).unwrap();
        let presence = PresenceSource::new(
            &cat,
            &prov,
            ids(&["g0", "g1", "g2"]),
            ok_all(vec![
                kmers_of(&[GENE_A]),
                kmers_of(&[GENE_B]),
                kmers_of(&[]),
            ]),
        )
        .unwrap();
        let panels = panel_source(&["g0", "g1", "g2"], vec![vec![], vec![], vec![]]);
        let matrix = matrix_fixture();
        let kmers = KmerSource::new(&matrix).unwrap();
        let probes = probe_source();
        // g7 is in no input: a genome without probe calls, catalog hits or a
        // matrix column.
        let asked = ids(&["g2", "g7", "g0"]);
        for raw in [
            admit(&probes, &asked).unwrap().raw().clone(),
            admit(&panels, &asked).unwrap().raw().clone(),
            admit(&presence, &asked).unwrap().raw().clone(),
            admit(&kmers, &asked).unwrap().raw().clone(),
        ] {
            assert_eq!(raw.sample_ids, ["g2", "g0"], "{}", raw.spec.id);
            check_block(&raw, &asked).unwrap();
        }
        // g2's probe no call keeps its row, as not assessable.
        let block = admit(&probes, &asked).unwrap();
        assert_eq!(binary(block.raw(), "C:gc_gyra"), (vec![], vec![1]));
        // A block of none of the requested samples has no rows.
        let block = admit(&kmers, &ids(&["g9"])).unwrap();
        assert!(block.raw().sample_ids.is_empty());
        assert_eq!(block.raw().names.len(), 3);
    }

    #[test]
    fn admission_checks_rows_and_shape() {
        let source = probe_source();
        assert!(admit(&source, &ids(&["g0", "g0"])).is_err());
        let samples = ids(&["g0", "g1"]);
        let good = source.build(&samples).unwrap();
        check_block(&good, &samples).unwrap();
        let mut short = good.clone();
        short.cols[0] = RawColumn::Binary {
            present: Bits::zeros(1),
            assessable: Bits::zeros(1),
        };
        let mut leak = good.clone();
        leak.cols[0] = RawColumn::Binary {
            present: Bits::ones(2),
            assessable: Bits::zeros(2),
        };
        let mut foreign = good.clone();
        foreign.names[0] = "K:p0".into();
        let mut twice = good.clone();
        twice.names.push(twice.names[0].clone());
        twice.cols.push(twice.cols[0].clone());
        let mut level = good.clone();
        level.cols[0] = RawColumn::Categorical {
            codes: vec![Some(0), Some(1)],
            levels: vec!["a".into()],
        };
        let mut counts = good.clone();
        counts.cols[0] = RawColumn::Count {
            values: vec![Some(1)],
        };
        let mut unasked = good.clone();
        unasked.sample_ids[1] = "g2".into();
        for bad in [short, leak, foreign, twice, level, counts, unasked] {
            assert!(check_block(&bad, &samples).is_err());
        }
        // Rows keep the requested order.
        assert!(check_block(&good, &ids(&["g1", "g0"])).is_err());
    }

    #[test]
    fn columns_build_from_plain_values() {
        let RawColumn::Binary {
            present,
            assessable,
        } = RawColumn::binary(&[true, false, false], &[true, true, false]).unwrap()
        else {
            panic!("not binary")
        };
        assert_eq!(present.iter_ones().collect::<Vec<_>>(), [0]);
        assert_eq!(assessable.iter_ones().collect::<Vec<_>>(), [0, 1]);
        assert_eq!(present.len(), 3);
        // Present must be assessable, and the two must be the same length.
        assert!(RawColumn::binary(&[true], &[false]).is_err());
        assert!(RawColumn::binary(&[true], &[true, true]).is_err());
        let RawColumn::Binary {
            present,
            assessable,
        } = RawColumn::binary_assessable(&[false, true])
        else {
            panic!("not binary")
        };
        assert_eq!(present.iter_ones().collect::<Vec<_>>(), [1]);
        assert_eq!(assessable.count(), 2);
        let levels = ids(&["a", "b"]);
        assert!(RawColumn::categorical(vec![Some(1), None, Some(0)], levels.clone()).is_ok());
        assert!(RawColumn::categorical(vec![Some(2)], levels).is_err());
        assert!(RawColumn::categorical(vec![Some(u16::MAX)], vec![]).is_err());
        // A block of constructed columns passes the admission shape check.
        let mut block = probe_source().build(&ids(&["g0", "g1"])).unwrap();
        block.cols = vec![RawColumn::binary(&[false, true], &[true, true]).unwrap()];
        check_block(&block, &ids(&["g0", "g1"])).unwrap();
        block.cols = vec![RawColumn::count(vec![Some(3), None])];
        check_block(&block, &ids(&["g0", "g1"])).unwrap();
    }

    #[test]
    fn building_takes_no_labels() {
        // The adapters' code reaches no module that holds phenotypes.
        let code = include_str!("blocks.rs")
            .split("#[cfg(test)]")
            .next()
            .unwrap();
        // Every path into the crate is spelled `crate::`, so the list below
        // sees it: no `super::` (which reaches `crate::` from here) and no
        // `self::` (which can chain to `super::`).
        for path in ["super::", "self::"] {
            let found: Vec<usize> = code
                .match_indices(path)
                .map(|(at, _)| at)
                .filter(|&at| {
                    !code[..at].ends_with(|c: char| c.is_ascii_alphanumeric() || c == '_')
                })
                .collect();
            assert!(found.is_empty(), "{path} at bytes {found:?}");
        }
        let used: BTreeSet<&str> = code
            .match_indices("crate::")
            .map(|(at, _)| {
                let rest = &code[at + "crate::".len()..];
                let end = rest
                    .find(|c: char| !(c.is_ascii_alphanumeric() || c == '_'))
                    .unwrap_or(rest.len());
                &rest[..end]
            })
            .collect();
        assert_eq!(
            used,
            ["bits", "catalog", "dna", "matrix", "panels", "probes"]
                .into_iter()
                .collect()
        );
        // A block is a function of sample ids alone: reversing them reverses
        // every column, and building again gives the same block.
        let dir = catalog_fixture();
        let cat = Catalog::load(dir.path(), &KmerSpace::new(31)).unwrap();
        let prov = Provenance::catalog_dir(dir.path()).unwrap();
        let forward = ids(&["g0", "g1", "g2"]);
        let backward: Vec<String> = forward.iter().rev().cloned().collect();
        let presence = PresenceSource::new(
            &cat,
            &prov,
            forward.clone(),
            ok_all(vec![
                kmers_of(&[GENE_A]),
                kmers_of(&[GENE_B]),
                kmers_of(&[]),
            ]),
        )
        .unwrap();
        let panels = panel_source(
            &["g0", "g1", "g2"],
            vec![
                vec![hit("blaKPC", "blaKPC")],
                vec![],
                vec![hit("mecA", "mecA")],
            ],
        );
        let matrix = matrix_fixture();
        let kmers = KmerSource::new(&matrix).unwrap();
        let probes = probe_source();
        let sources: [&dyn BlockSource; 4] = [&probes, &panels, &presence, &kmers];
        let calls = |raw: &RawBlock| -> Vec<(Vec<usize>, Vec<usize>)> {
            raw.names.iter().map(|n| binary(raw, n)).collect()
        };
        for source in sources {
            let a = source.build(&forward).unwrap();
            let b = source.build(&backward).unwrap();
            assert_eq!(a.names, b.names);
            let flip = |v: Vec<usize>| -> Vec<usize> {
                let mut v: Vec<usize> = v.into_iter().map(|i| 2 - i).collect();
                v.sort_unstable();
                v
            };
            let flipped: Vec<_> = calls(&a)
                .into_iter()
                .map(|(p, s)| (flip(p), flip(s)))
                .collect();
            assert_eq!(flipped, calls(&b), "{}", a.spec.id);
            assert_eq!(calls(&a), calls(&source.build(&forward).unwrap()));
        }
    }
}
