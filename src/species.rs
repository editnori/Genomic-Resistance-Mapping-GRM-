//! Per-species pipeline stages: count, structure, matrix.
//!
//!   count      FASTA -> kset per genome, plus assembly stats, sketch, and
//!              catalog gene calls (cached in a JSON sidecar per genome)
//!   structure  sketches -> Mash distances -> lineage clusters, dendrogram,
//!              PCoA map
//!   matrix     ksets -> pattern matrix (see `matrix`)

use crate::bvbrc::fasta_path;
use crate::catalog::{self, Catalog, FamilyHit};
use crate::cluster;
use crate::dna::KmerSpace;
use crate::fasta::SeqStats;
use crate::kset::{self, KsetReader};
use crate::layout::{Layout, remove_if_present};
use crate::matrix;
use anyhow::{Context, Result};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

pub const SKETCH: usize = 1000;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Counted {
    pub genome_id: String,
    pub n_kmers: u64,
    pub stats: SeqStats,
    pub catalog_version: String,
    /// `kset::digest` of the kset this count wrote; None in counts recorded
    /// before digests.
    #[serde(default)]
    pub kset_digest: Option<String>,
    pub genes: Vec<FamilyHit>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CountSummary {
    pub k: usize,
    pub genomes: Vec<Counted>,
    pub missing_fasta: Vec<String>,
    pub failed: Vec<(String, String)>,
}

impl CountSummary {
    /// The first genome counted before kset digests were recorded.
    pub fn undigested(&self) -> Option<&Counted> {
        self.genomes.iter().find(|g| g.kset_digest.is_none())
    }

    /// The first genome whose gene calls came from another catalog than
    /// `calls_version`: a `Catalog::calls_version`, or "" for none.
    pub fn called_with_other_catalog(&self, calls_version: &str) -> Option<&Counted> {
        self.genomes
            .iter()
            .find(|g| g.catalog_version != calls_version)
    }

    /// Each genome's kset digest in counted order, which `grm {stage}`
    /// records beside what it builds from the ksets.
    fn kset_digests(&self, stage: &str) -> Result<Vec<String>> {
        self.genomes
            .iter()
            .map(|g| g.kset_digest.clone())
            .collect::<Option<_>>()
            .with_context(|| format!("run `grm count` to record kset digests before `grm {stage}`"))
    }
}

fn read_ids(path: &Path) -> Result<Vec<String>> {
    Ok(std::fs::read_to_string(path)
        .with_context(|| format!("read {}", path.display()))?
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .map(String::from)
        .collect())
}

/// genomes.txt split into the genomes with a FASTA file, which `count`
/// counts in this order, and those without one.
fn selection(layout: &Layout, slug: &str) -> Result<(Vec<String>, Vec<String>)> {
    let fasta_dir = layout.fasta_dir();
    Ok(read_ids(&layout.species_dir(slug).join("genomes.txt"))?
        .into_iter()
        .partition(|id| fasta_path(&fasta_dir, id).exists()))
}

fn count_one(
    id: &str,
    fasta: &Path,
    kset_path: &Path,
    sidecar: &Path,
    space: &KmerSpace,
    catalog: Option<&Catalog>,
    (values, scratch): &mut (Vec<u64>, Vec<u64>),
) -> Result<Counted> {
    let version = catalog.map_or(String::new(), Catalog::calls_version);
    // A kset counted under another `--k` is not a cache hit.
    if KsetReader::open(kset_path).is_ok_and(|r| r.k == space.k)
        && let Ok(bytes) = std::fs::read(sidecar)
        && let Ok(c) = serde_json::from_slice::<Counted>(&bytes)
    {
        if c.catalog_version == version && c.kset_digest.is_some() {
            return Ok(c);
        }
        // Catalog changed, or the sidecar predates kset digests:
        // redo both from the stored kset, with its k-mer count. A
        // sidecar from before digests can describe a kset of another
        // k, since those counts wrote the kset without removing it.
        let values = KsetReader::open(kset_path)?.read_all()?;
        let genes = catalog.map_or(vec![], |cat| {
            catalog::family_calls(cat, &cat.detect(&values, catalog::CALL_CONTAINMENT))
        });
        let c = Counted {
            n_kmers: values.len() as u64,
            catalog_version: version,
            kset_digest: Some(kset::digest(&values)),
            genes,
            ..c
        };
        std::fs::write(sidecar, serde_json::to_vec(&c)?)?;
        return Ok(c);
    }
    let stats = kset::count_into(fasta, space, 1, values, scratch)?;
    // The sidecar describes the kset beside it. Removing it first means a run
    // stopped between the two writes leaves a kset with no sidecar, which the
    // next run recounts, instead of the new kset with the old k-mer count and
    // gene calls.
    remove_if_present(sidecar).with_context(|| format!("remove {}", sidecar.display()))?;
    kset::write(kset_path, id, space, values)?;
    let genes = catalog.map_or(vec![], |cat| {
        catalog::family_calls(cat, &cat.detect(values, catalog::CALL_CONTAINMENT))
    });
    let c = Counted {
        genome_id: id.to_string(),
        n_kmers: values.len() as u64,
        stats,
        catalog_version: version,
        kset_digest: Some(kset::digest(values)),
        genes,
    };
    std::fs::write(sidecar, serde_json::to_vec(&c)?)?;
    Ok(c)
}

pub fn count(
    layout: &Layout,
    slug: &str,
    k: usize,
    catalog: Option<&Catalog>,
) -> Result<CountSummary> {
    if catalog.is_none() {
        // Counting without a catalog records no gene calls, and learning
        // without them leaves every report without its catalog baseline.
        let previous = load_counted(layout, slug).ok();
        if let Some(g) = previous
            .as_ref()
            .and_then(|c| c.called_with_other_catalog(""))
        {
            anyhow::bail!(
                "{slug}: no catalog loaded, and counting without one would drop the gene calls counted.json holds from catalog {}; restore it with `grm fetch catalog`, or delete counted.json to count without one",
                g.catalog_version
            );
        }
    }
    let space = KmerSpace::new(k);
    let (present, missing) = selection(layout, slug)?;
    let kdir = layout.kset_dir(slug);
    std::fs::create_dir_all(&kdir)?;
    let fasta_dir = layout.fasta_dir();
    let done = AtomicUsize::new(0);
    let results: Vec<(String, Result<Counted>)> = present
        .par_iter()
        .map_init(
            || (Vec::new(), Vec::new()),
            |buffers, id| {
                let r = count_one(
                    id,
                    &fasta_path(&fasta_dir, id),
                    &kdir.join(format!("{id}.kset")),
                    &kdir.join(format!("{id}.json")),
                    &space,
                    catalog,
                    buffers,
                );
                let d = done.fetch_add(1, Ordering::Relaxed) + 1;
                if d.is_multiple_of(200) || d == present.len() {
                    eprintln!("  {slug}: counted {d} / {}", present.len());
                }
                (id.clone(), r)
            },
        )
        .collect();
    let mut genomes = Vec::new();
    let mut failed = Vec::new();
    for (id, r) in results {
        match r {
            Ok(c) => genomes.push(c),
            Err(e) => failed.push((id, format!("{e:#}"))),
        }
    }
    let summary = CountSummary {
        k,
        genomes,
        missing_fasta: missing,
        failed,
    };
    std::fs::write(
        layout.species_dir(slug).join("counted.json"),
        serde_json::to_vec(&summary)?,
    )?;
    Ok(summary)
}

pub fn load_counted(layout: &Layout, slug: &str) -> Result<CountSummary> {
    let path = layout.species_dir(slug).join("counted.json");
    Ok(serde_json::from_slice(
        &std::fs::read(&path)
            .with_context(|| format!("run `grm count` first ({})", path.display()))?,
    )?)
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Structure {
    pub genomes: Vec<String>,
    /// k the genomes were counted at; 0 in files written before it was
    /// recorded.
    #[serde(default)]
    pub k: usize,
    /// Each genome's kset digest from counted.json; empty in files written
    /// before digests.
    #[serde(default)]
    pub kset_digests: Vec<String>,
    pub lineages: cluster::Lineages,
    pub merges: Vec<cluster::Merge>,
    pub leaf_order: Vec<u32>,
    pub coords: Vec<[f32; 2]>,
    /// Distance to the closest other genome (Mash), per genome.
    pub nearest: Vec<f32>,
    pub sketch_size: usize,
}

pub fn structure(layout: &Layout, slug: &str) -> Result<Structure> {
    let counted = load_counted(layout, slug)?;
    let kset_digests = counted.kset_digests("structure")?;
    let kdir = layout.kset_dir(slug);
    let sketches: Vec<Vec<u64>> = counted
        .genomes
        .par_iter()
        .map(|g| KsetReader::open(&kdir.join(format!("{}.kset", g.genome_id)))?.sketch(SKETCH))
        .collect::<Result<_>>()?;
    let n = sketches.len();
    let dist = cluster::distance_matrix(&sketches, SKETCH, counted.k);
    let merges = cluster::average_linkage(&dist, n);
    let lineages = cluster::choose_lineages(&merges, n, 20, 0.25);
    let nearest = (0..n)
        .map(|i| {
            (0..n)
                .filter(|&j| j != i)
                .map(|j| dist[i * n + j])
                .fold(f32::INFINITY, f32::min)
        })
        .collect();
    let s = Structure {
        genomes: counted
            .genomes
            .iter()
            .map(|g| g.genome_id.clone())
            .collect(),
        k: counted.k,
        kset_digests,
        leaf_order: cluster::leaf_order(&merges, n),
        coords: cluster::pcoa(&dist, n),
        lineages,
        merges,
        nearest,
        sketch_size: SKETCH,
    };
    std::fs::write(
        layout.species_dir(slug).join("structure.json"),
        serde_json::to_vec(&s)?,
    )?;
    Ok(s)
}

/// Merges ksets in counted order, which becomes the matrix column order.
pub fn build_matrix(
    layout: &Layout,
    slug: &str,
    min_support: u32,
) -> Result<matrix::MatrixSummary> {
    let counted = load_counted(layout, slug)?;
    let digests = counted.kset_digests("matrix")?;
    let kdir = layout.kset_dir(slug);
    let paths: Vec<PathBuf> = counted
        .genomes
        .iter()
        .map(|g| kdir.join(format!("{}.kset", g.genome_id)))
        .collect();
    let mut m = matrix::build(&paths, min_support, rayon::current_num_threads() * 2)?;
    m.summary.kset_digests = digests;
    m.save(&layout.matrix_dir(slug))?;
    Ok(m.summary)
}

/// Stages to rerun before learning so that the gene calls, structure, and
/// matrix all come from the selected genomes, their k, and the catalog.
pub struct Stale {
    pub count: bool,
    pub structure: bool,
    pub matrix: bool,
}

/// Compares counted.json with genomes.txt and `catalog`, and structure.json
/// and matrix.json each with counted.json. A catalog that cannot be read
/// (None) is no reason to recount: counting without one records no calls.
pub fn stale_stages(layout: &Layout, slug: &str, catalog: Option<&Catalog>) -> Stale {
    fn read(path: PathBuf) -> Option<BuiltFrom> {
        serde_json::from_slice(&std::fs::read(path).ok()?).ok()
    }
    let with_fasta = selection(layout, slug).ok().map(|(present, _)| present);
    stale(
        load_counted(layout, slug).ok().as_ref(),
        with_fasta.as_deref(),
        read(layout.species_dir(slug).join("structure.json")).as_ref(),
        read(layout.matrix_dir(slug).join("matrix.json")).as_ref(),
        catalog.map(Catalog::calls_version).as_deref(),
    )
}

/// What structure.json and matrix.json record about the counts they were
/// built from.
#[derive(Deserialize)]
struct BuiltFrom {
    #[serde(default)]
    k: usize,
    genomes: Vec<String>,
    #[serde(default)]
    kset_digests: Vec<String>,
}

/// Structure and matrix are each stale when the file is missing or
/// `not_built_from` the counts. Counts are stale when they list other
/// genomes than `count` would count now: the genomes in genomes.txt with a
/// FASTA file (`with_fasta`, None when genomes.txt cannot be read), less
/// those whose count failed, in file order. They are also stale when a
/// genome's gene calls carry another catalog version than `calls_version`,
/// the loaded catalog's (None when none loaded). Count also reruns when
/// structure or matrix is stale, because their rebuild reads the ksets, which
/// `--drop-ksets` deletes and a recount writes again (reusing those still on
/// disk). That covers counts recorded before kset digests, which no
/// structure or matrix matches.
fn selection_matches(counted: &CountSummary, with_fasta: &[String]) -> bool {
    let failed: HashSet<&str> = counted.failed.iter().map(|(id, _)| id.as_str()).collect();
    with_fasta
        .iter()
        .filter(|id| !failed.contains(id.as_str()))
        .eq(counted.genomes.iter().map(|g| &g.genome_id))
}

/// Direct learning uses the same selected-genome freshness check as the
/// pipeline, but stops for a recount instead of running one implicitly.
pub fn ensure_selection_counted(layout: &Layout, slug: &str) -> Result<()> {
    let counted = load_counted(layout, slug)?;
    let (with_fasta, _) = selection(layout, slug)?;
    anyhow::ensure!(
        selection_matches(&counted, &with_fasta),
        "genomes.txt adds, removes, or reorders selected genomes relative to counted.json after excluding failed counts; rerun `grm count --species {slug}`, then `grm structure` and `grm matrix` before learning"
    );
    Ok(())
}

fn stale(
    counted: Option<&CountSummary>,
    with_fasta: Option<&[String]>,
    structure: Option<&BuiltFrom>,
    matrix: Option<&BuiltFrom>,
    calls_version: Option<&str>,
) -> Stale {
    let Some(counted) = counted else {
        return Stale {
            count: true,
            structure: true,
            matrix: true,
        };
    };
    let other_genomes = with_fasta.is_none_or(|ids| !selection_matches(counted, ids));
    let outdated = |file: Option<&BuiltFrom>| {
        file.is_none_or(|f| not_built_from(counted, f.k, &f.genomes, &f.kset_digests).is_some())
    };
    let (structure, matrix) = (outdated(structure), outdated(matrix));
    let other_catalog =
        calls_version.is_some_and(|v| counted.called_with_other_catalog(v).is_some());
    Stale {
        count: other_genomes || other_catalog || structure || matrix,
        structure,
        matrix,
    }
}

/// Why a structure or matrix recorded as built at `k` from ksets with
/// `digests`, for `genomes` in that order, does not come from `counted`, or
/// None when it does. Counts without digests match no digests.
pub(crate) fn not_built_from(
    counted: &CountSummary,
    k: usize,
    genomes: &[String],
    digests: &[String],
) -> Option<String> {
    if digests.is_empty() && !genomes.is_empty() {
        Some("records no kset digests (it predates them)".into())
    } else if k != counted.k {
        Some(format!(
            "was built at k = {k} and counted.json records k = {}",
            counted.k
        ))
    } else if !genomes
        .iter()
        .eq(counted.genomes.iter().map(|g| &g.genome_id))
    {
        Some("lists other genomes than counted.json, or the same genomes in another order".into())
    } else if !digests
        .iter()
        .map(Some)
        .eq(counted.genomes.iter().map(|g| g.kset_digest.as_ref()))
    {
        Some("was built from other ksets than counted.json records".into())
    } else {
        None
    }
}

/// The catalog `grm learn` uses beside the counted gene calls: the one in
/// the catalog directory when the calls came from it, or None when the
/// genomes were counted without a catalog and none loads. Calls made with
/// another catalog version, or without one while a catalog loads, are an
/// error, which `grm pipeline` avoids by recounting first (see `stale`).
/// So are calls from a catalog that no longer loads, rather than reports
/// without their catalog baseline.
pub fn catalog_of_calls(
    layout: &Layout,
    counted: &CountSummary,
    space: &KmerSpace,
) -> Result<Option<Catalog>> {
    let dir = layout.catalog_dir();
    let loaded = Catalog::load(&dir, space);
    let version = loaded
        .as_ref()
        .map_or(String::new(), Catalog::calls_version);
    let Some(g) = counted.called_with_other_catalog(&version) else {
        return Ok(loaded.ok());
    };
    match loaded {
        Ok(_) if g.catalog_version.is_empty() => anyhow::bail!(
            "{} was counted without a catalog, and {} now holds catalog {version}; rerun `grm count`",
            g.genome_id,
            dir.display()
        ),
        Ok(_) => anyhow::bail!(
            "{} has gene calls from catalog {}, and {} now holds catalog {version}; rerun `grm count`",
            g.genome_id,
            g.catalog_version,
            dir.display()
        ),
        Err(e) => Err(e.context(format!(
            "{} has gene calls from catalog {}, which no longer loads",
            g.genome_id, g.catalog_version
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reuses_cached_ksets_only_at_the_requested_k() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("g.fna");
        let seq = ">c\nACGTTGCAAGGCTTACCGATAGCTAGGCTTAACGGTACCATGCAGTTCAGGA\n";
        std::fs::write(&fasta, seq).unwrap();
        let (kset_path, json) = (dir.path().join("g.kset"), dir.path().join("g.json"));
        let mut buffers = (Vec::new(), Vec::new());
        let mut count = |k| {
            let space = KmerSpace::new(k);
            count_one("g", &fasta, &kset_path, &json, &space, None, &mut buffers).unwrap()
        };
        let k21 = count(21);
        let k31 = count(31);
        let (values, _) = kset::count_file(&fasta, &KmerSpace::new(31), 1).unwrap();
        assert_eq!(k31.n_kmers, values.len() as u64);
        assert_ne!(k31.n_kmers, k21.n_kmers);
        assert_eq!(KsetReader::open(&kset_path).unwrap().k, 31);
        // A cache at the requested k is still used: the FASTA is not reread.
        std::fs::remove_file(&fasta).unwrap();
        assert_eq!(count(31).n_kmers, k31.n_kmers);
    }

    #[test]
    fn a_recount_never_leaves_the_old_sidecar_beside_a_new_kset() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("g.fna");
        std::fs::write(&fasta, format!(">c\n{SEQ}\n")).unwrap();
        let (kset_path, json) = (dir.path().join("g.kset"), dir.path().join("g.json"));
        let mut buffers = (Vec::new(), Vec::new());
        let mut count = |k| {
            let space = KmerSpace::new(k);
            count_one("g", &fasta, &kset_path, &json, &space, None, &mut buffers)
        };
        let k21 = count(21).unwrap();
        // A directory where `kset::write` creates its temporary file stops a
        // recount at k = 31 at the kset write, for root as for anyone. The old
        // sidecar has to be gone by then: a run killed, or a disk filled, just
        // after the write would leave it beside the new kset.
        let part = dir.path().join("g.kset.part");
        std::fs::create_dir(&part).unwrap();
        assert!(count(31).is_err());
        assert!(!json.exists());
        std::fs::remove_dir(&part).unwrap();
        let k31 = count(31).unwrap();
        let (values, _) = kset::count_file(&fasta, &KmerSpace::new(31), 1).unwrap();
        assert_eq!(k31.n_kmers, values.len() as u64);
        assert_ne!(k31.n_kmers, k21.n_kmers);
    }

    #[test]
    fn a_sidecar_from_before_digests_beside_a_kset_of_another_k_takes_the_ksets_count() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("g.fna");
        std::fs::write(&fasta, format!(">c\n{SEQ}\n")).unwrap();
        let (kset_path, json) = (dir.path().join("g.kset"), dir.path().join("g.json"));
        let mut buffers = (Vec::new(), Vec::new());
        let mut count = |k| {
            let space = KmerSpace::new(k);
            count_one("g", &fasta, &kset_path, &json, &space, None, &mut buffers).unwrap()
        };
        let k21 = count(21);
        let mut old: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&json).unwrap()).unwrap();
        let k31 = count(31);
        assert_ne!(k21.n_kmers, k31.n_kmers);
        // The k = 21 sidecar as `grm count` wrote it before kset digests,
        // beside the k = 31 kset: a recount from before digests wrote the
        // kset without removing the sidecar, then stopped.
        old.as_object_mut().unwrap().remove("kset_digest");
        std::fs::write(&json, serde_json::to_vec(&old).unwrap()).unwrap();
        std::fs::remove_file(&fasta).unwrap();
        let again = count(31);
        assert_eq!(again.n_kmers, k31.n_kmers);
        assert_eq!(again.kset_digest, k31.kset_digest);
    }

    const SEQ: &str = "ACGTTGCAAGGCTTACCGATAGCTAGGCTTAACGGTACCATGCAGTTCAGGA";

    #[test]
    fn kset_digests_tell_ksets_of_one_size_apart_and_a_recount_repeats_them() {
        let dir = tempfile::tempdir().unwrap();
        let space = KmerSpace::new(21);
        let mut buffers = (Vec::new(), Vec::new());
        let mut count = |id: &str, seq: &str| {
            let fasta = dir.path().join(format!("{id}.fna"));
            std::fs::write(&fasta, format!(">c\n{seq}\n")).unwrap();
            let (kset_path, json) = (
                dir.path().join(format!("{id}.kset")),
                dir.path().join(format!("{id}.json")),
            );
            count_one(id, &fasta, &kset_path, &json, &space, None, &mut buffers).unwrap()
        };
        // One base changed: the 21 k-mers over it give way to 21 others.
        let mut other = SEQ.to_string();
        other.replace_range(25..26, "T");
        let (a, b) = (count("a", SEQ), count("b", &other));
        assert_eq!(a.n_kmers, b.n_kmers);
        assert_ne!(a.kset_digest, b.kset_digest);
        // As after `--drop-ksets`: a recount of the same assembly repeats it.
        std::fs::remove_file(dir.path().join("a.kset")).unwrap();
        assert_eq!(count("a", SEQ).kset_digest, a.kset_digest);
    }

    #[test]
    fn a_sidecar_without_a_kset_digest_takes_it_from_the_stored_kset() {
        let dir = tempfile::tempdir().unwrap();
        let fasta = dir.path().join("g.fna");
        std::fs::write(&fasta, format!(">c\n{SEQ}\n")).unwrap();
        let (kset_path, json) = (dir.path().join("g.kset"), dir.path().join("g.json"));
        let space = KmerSpace::new(21);
        let mut buffers = (Vec::new(), Vec::new());
        let mut count = || count_one("g", &fasta, &kset_path, &json, &space, None, &mut buffers);
        let counted = count().unwrap();
        // The sidecar as `grm count` wrote it before kset digests.
        let mut old: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&json).unwrap()).unwrap();
        old.as_object_mut().unwrap().remove("kset_digest");
        std::fs::write(&json, serde_json::to_vec(&old).unwrap()).unwrap();
        // Without the assembly, only the stored kset can supply it.
        std::fs::remove_file(&fasta).unwrap();
        let again = count().unwrap();
        assert!(counted.kset_digest.is_some());
        assert_eq!(again.kset_digest, counted.kset_digest);
    }

    fn counted(k: usize, genomes: &[&str], calls_version: &str) -> CountSummary {
        CountSummary {
            k,
            genomes: genomes
                .iter()
                .map(|id| Counted {
                    genome_id: id.to_string(),
                    n_kmers: 0,
                    stats: SeqStats::default(),
                    catalog_version: calls_version.into(),
                    kset_digest: Some(format!("kset of {id}")),
                    genes: vec![],
                })
                .collect(),
            missing_fasta: vec![],
            failed: vec![],
        }
    }

    /// What a structure or matrix built at `k` for `genomes` records, with
    /// the kset digests `counted` gives them.
    fn built(k: usize, genomes: &[String]) -> BuiltFrom {
        BuiltFrom {
            k,
            genomes: genomes.to_vec(),
            kset_digests: genomes.iter().map(|g| format!("kset of {g}")).collect(),
        }
    }

    fn ids(ids: &[&str]) -> Vec<String> {
        ids.iter().map(|g| g.to_string()).collect()
    }

    #[test]
    fn structure_and_matrix_each_go_stale_on_their_own_k_genomes_and_ksets() {
        let (abc, cba) = (ids(&["a", "b", "c"]), ids(&["c", "b", "a"]));
        let at_31 = counted(31, &["a", "b", "c"], "v2");
        let plan = |counted: Option<&CountSummary>,
                    structure: Option<BuiltFrom>,
                    matrix: Option<BuiltFrom>| {
            let s = stale(
                counted,
                Some(&abc),
                structure.as_ref(),
                matrix.as_ref(),
                Some("v2"),
            );
            (s.count, s.structure, s.matrix)
        };
        let current = || Some(built(31, &abc));
        let fresh = |calls| {
            let (structure, matrix) = (built(31, &abc), built(31, &abc));
            let s = stale(
                Some(&at_31),
                Some(&abc),
                Some(&structure),
                Some(&matrix),
                calls,
            );
            (s.count, s.structure, s.matrix)
        };
        assert_eq!(fresh(Some("v2")), (false, false, false));
        // A new catalog redoes the gene calls and keeps the k-mer stages.
        assert_eq!(fresh(Some("v3")), (true, false, false));
        // An unreadable catalog leaves the gene calls as they are.
        assert_eq!(fresh(None), (false, false, false));
        // Counts recorded before kset digests are redone, and only then can
        // structure and matrix be compared with them.
        let mut undigested = at_31.clone();
        undigested.genomes[1].kset_digest = None;
        assert_eq!(
            plan(Some(&undigested), current(), current()),
            (true, true, true)
        );
        // The genomes were counted at another k than both were built at.
        let at_21 = counted(21, &["a", "b", "c"], "v2");
        assert_eq!(plan(Some(&at_21), current(), current()), (true, true, true));
        // The same genomes in another order are not the same matrix columns.
        let reordered = plan(Some(&at_31), current(), Some(built(31, &cba)));
        assert_eq!(reordered, (true, false, true));
        let structure_only = plan(Some(&at_31), Some(built(31, &cba)), current());
        assert_eq!(structure_only, (true, true, false));
        let fewer = plan(Some(&at_31), current(), Some(built(31, &abc[..2])));
        assert_eq!(fewer, (true, false, true));
        let other_k = plan(Some(&at_31), Some(built(21, &abc)), current());
        assert_eq!(other_k, (true, true, false));
        let mut changed_kset = at_31.clone();
        changed_kset.genomes[1].kset_digest = Some("changed kset".into());
        assert_eq!(
            plan(Some(&changed_kset), current(), current()),
            (true, true, true)
        );
        // Files written before structure.json recorded k, and before either
        // recorded kset digests, match no counts.
        let old_structure = BuiltFrom {
            k: 0,
            kset_digests: vec![],
            ..built(31, &abc)
        };
        assert_eq!(
            plan(Some(&at_31), Some(old_structure), current()),
            (true, true, false)
        );
        let old_matrix = BuiltFrom {
            kset_digests: vec![],
            ..built(31, &abc)
        };
        assert_eq!(
            plan(Some(&at_31), current(), Some(old_matrix)),
            (true, false, true)
        );
        // A missing file is stale; without counted.json nothing can be checked.
        assert_eq!(plan(None, current(), current()), (true, true, true));
        assert_eq!(plan(Some(&at_31), None, current()), (true, true, false));
        assert_eq!(plan(Some(&at_31), current(), None), (true, false, true));
    }

    #[test]
    fn counts_go_stale_when_genomes_txt_selects_other_genomes_with_fasta() {
        let abc = ids(&["a", "b", "c"]);
        let at_31 = counted(31, &["a", "b", "c"], "v2");
        // Structure and matrix still match counted.json here; they are
        // compared again with the counts the rerun writes.
        let plan = |with_fasta: Option<&[&str]>, counted: &CountSummary| {
            let with_fasta = with_fasta.map(ids);
            let (structure, matrix) = (built(31, &abc), built(31, &abc));
            let s = stale(
                Some(counted),
                with_fasta.as_deref(),
                Some(&structure),
                Some(&matrix),
                Some("v2"),
            );
            (s.count, s.structure, s.matrix)
        };
        assert_eq!(plan(Some(&["a", "b", "c"]), &at_31), (false, false, false));
        for selection in [
            &["a", "b", "d"][..],
            &["a", "b", "c", "d"],
            &["a", "b"],
            &["b", "a", "c"],
        ] {
            assert_eq!(
                plan(Some(selection), &at_31),
                (true, false, false),
                "{selection:?}"
            );
        }
        assert_eq!(plan(None, &at_31), (true, false, false));
        // A genome whose count failed stays out, as `count` left it.
        let mut failed = at_31.clone();
        failed.failed.push(("x".into(), "unreadable FASTA".into()));
        assert_eq!(
            plan(Some(&["a", "x", "b", "c"]), &failed),
            (false, false, false)
        );
    }

    #[test]
    fn learning_takes_only_the_catalog_the_gene_calls_came_from() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let space = KmerSpace::new(31);
        let catalog_for = |calls: &str| {
            catalog_of_calls(&layout, &counted(31, &["a", "b"], calls), &space)
                .map(|c| c.map(|c| c.version))
                .map_err(|e| format!("{e:#}"))
        };
        assert_eq!(catalog_for(""), Ok(None));
        let err = catalog_for("v1|group-calls-v2").unwrap_err();
        assert!(
            err.contains("a has gene calls from catalog v1|group-calls-v2, which no longer loads: run `grm fetch catalog` first"),
            "{err}"
        );
        let catalog_dir = layout.catalog_dir();
        std::fs::create_dir_all(&catalog_dir).unwrap();
        std::fs::write(catalog_dir.join("ReferenceGeneCatalog.txt"), "type\n").unwrap();
        std::fs::write(catalog_dir.join("AMR_CDS.fa"), "").unwrap();
        std::fs::write(catalog_dir.join("version.txt"), "v1").unwrap();
        assert_eq!(catalog_for("v1|group-calls-v2"), Ok(Some("v1".into())));
        let err = catalog_for("v0|group-calls-v2").unwrap_err();
        assert!(
            err.contains("a has gene calls from catalog v0|group-calls-v2, and"),
            "{err}"
        );
        let err = catalog_for("").unwrap_err();
        assert!(
            err.contains("a was counted without a catalog, and"),
            "{err}"
        );
    }
}
