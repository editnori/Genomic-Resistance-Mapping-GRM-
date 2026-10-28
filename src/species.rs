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
use crate::confine::{self, sync_dir};
use crate::dna::KmerSpace;
use crate::fasta::SeqStats;
use crate::kset::{self, KsetReader};
use crate::layout::{Layout, remove_if_present};
use crate::matrix;
use anyhow::{Context, Result};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

pub const SKETCH: usize = 1000;
pub const DISK_SPARSE_RESERVE: u64 = 30 * 1024 * 1024 * 1024;

/// Free disk that writers running at once share. Each write claims its exact
/// length under one lock, against free space less what earlier claims still
/// hold, and keeps the claim until its file is in place. Without the ledger,
/// several counting threads could each pass the same check and write
/// together past the reserve.
struct Reserve {
    bytes: u64,
    claimed: Mutex<u64>,
    free: fn(&Path) -> Result<u64>,
}

/// Released when dropped, once the claimed file is in place (and so counted
/// by the filesystem's own free space).
struct Claim<'a> {
    reserve: &'a Reserve,
    bytes: u64,
}

impl Reserve {
    fn new(bytes: u64) -> Self {
        Reserve {
            bytes,
            claimed: Mutex::new(0),
            free: matrix::free_bytes,
        }
    }

    fn claim(&self, dir: &Path, need: u64) -> Result<Claim<'_>> {
        if self.bytes == 0 {
            return Ok(Claim {
                reserve: self,
                bytes: 0,
            });
        }
        let mut claimed = self
            .claimed
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let free = (self.free)(dir)?;
        anyhow::ensure!(
            free >= self.bytes.saturating_add(*claimed).saturating_add(need),
            "{} has {free} bytes free: writing {need} more beside the {} bytes other writers have claimed would cut into the {}-byte reserve",
            dir.display(),
            *claimed,
            self.bytes
        );
        *claimed += need;
        Ok(Claim {
            reserve: self,
            bytes: need,
        })
    }
}

impl Drop for Claim<'_> {
    fn drop(&mut self) {
        if self.bytes != 0 {
            let mut claimed = self
                .reserve
                .claimed
                .lock()
                .unwrap_or_else(std::sync::PoisonError::into_inner);
            *claimed -= self.bytes;
        }
    }
}

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
    /// Bottom-sketch of the same values. `None` only in old counts; `Some([])`
    /// is a valid sketch for an all-empty kset.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sketch: Option<Vec<u64>>,
    pub genes: Vec<FamilyHit>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CountSummary {
    pub k: usize,
    pub genomes: Vec<Counted>,
    pub missing_fasta: Vec<String>,
    pub failed: Vec<(String, String)>,
    /// Selected genomes not counted yet in an interrupted disk-sparse run.
    /// Empty in every complete manifest, including manifests from old GRM.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub pending: Vec<String>,
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
        anyhow::ensure!(
            self.pending.is_empty(),
            "counted.json has {} selected genomes pending; finish the disk-sparse run before `grm {stage}`",
            self.pending.len()
        );
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

/// Writes a sidecar through `<id>.json.part` and a rename, never through a
/// link (see `confine`).
fn write_sidecar(path: &Path, counted: &Counted, reserve: &Reserve) -> Result<()> {
    let bytes = serde_json::to_vec(counted)?;
    let _claim = reserve.claim(path.parent().unwrap_or(Path::new(".")), bytes.len() as u64)?;
    confine::write_atomic(path, &bytes, false)
}

#[allow(clippy::too_many_arguments)]
fn count_one(
    id: &str,
    fasta: &Path,
    kset_path: &Path,
    sidecar: &Path,
    space: &KmerSpace,
    catalog: Option<&Catalog>,
    buffers: &mut (Vec<u64>, Vec<u64>),
    reserve: u64,
) -> Result<Counted> {
    let reserve = Reserve::new(reserve);
    count_one_with(
        id, fasta, kset_path, sidecar, space, catalog, buffers, &reserve, false,
    )
}

/// `count_one` against a shared reserve. With `keep`, a kset already at
/// `kset_path` is an input that this count did not create: it is read, and
/// may get a new sidecar, but is never rewritten. One at another k, or one
/// that differs from a fresh count of its FASTA, is refused instead.
#[allow(clippy::too_many_arguments)]
fn count_one_with(
    id: &str,
    fasta: &Path,
    kset_path: &Path,
    sidecar: &Path,
    space: &KmerSpace,
    catalog: Option<&Catalog>,
    (values, scratch): &mut (Vec<u64>, Vec<u64>),
    reserve: &Reserve,
    keep: bool,
) -> Result<Counted> {
    let version = catalog.map_or(String::new(), Catalog::calls_version);
    // The sidecar is GRM's own file: never read or written through a link.
    let has_sidecar = confine::regular_or_absent(sidecar)?;
    let kept = keep && std::fs::symlink_metadata(kset_path).is_ok();
    if kept {
        let reader = KsetReader::open(kset_path).with_context(|| {
            format!(
                "read the kset {} that this disk-sparse run did not create; it never rewrites one, so move it aside or delete it to have it recounted",
                kset_path.display()
            )
        })?;
        anyhow::ensure!(
            reader.k == space.k,
            "{} holds {}-mers, not the {}-mers this run counts; disk-sparse mode never rewrites a kset it did not create, so move it aside or delete it to have it recounted",
            kset_path.display(),
            reader.k,
            space.k
        );
    }
    // A kset counted under another `--k` is not a cache hit.
    if KsetReader::open(kset_path).is_ok_and(|r| r.k == space.k)
        && has_sidecar
        && let Ok(bytes) = std::fs::read(sidecar)
        && let Ok(c) = serde_json::from_slice::<Counted>(&bytes)
    {
        if c.catalog_version == version && c.kset_digest.is_some() && c.sketch.is_some() {
            return Ok(c);
        }
        // Catalog changed, or the sidecar predates kset digests/sketches:
        // redo those from the stored kset, with its k-mer count. A sidecar
        // from before digests can describe a kset of another k.
        let values = KsetReader::open(kset_path)?.read_all()?;
        let genes = catalog.map_or(vec![], |cat| {
            catalog::family_calls(cat, &cat.detect(&values, catalog::CALL_CONTAINMENT))
        });
        let c = Counted {
            n_kmers: values.len() as u64,
            catalog_version: version,
            kset_digest: Some(kset::digest(&values)),
            sketch: Some(values.iter().copied().take(SKETCH).collect()),
            genes,
            ..c
        };
        write_sidecar(sidecar, &c, reserve)?;
        return Ok(c);
    }
    let stats = kset::count_into(fasta, space, 1, values, scratch)?;
    // The sidecar describes the kset beside it. Removing it first means a run
    // stopped between the two writes leaves a kset with no sidecar, which the
    // next run recounts, instead of the new kset with the old k-mer count and
    // gene calls.
    remove_if_present(sidecar).with_context(|| format!("remove {}", sidecar.display()))?;
    if kept {
        // Only the sidecar was missing or unreadable. The kept kset stays
        // byte for byte when it is exactly what this count gives.
        let stored = KsetReader::open(kset_path)?.read_all()?;
        anyhow::ensure!(
            stored == *values,
            "{} differs from a fresh count of {}; disk-sparse mode never rewrites a kset it did not create, so move it aside or delete it to have it recounted",
            kset_path.display(),
            fasta.display()
        );
    } else {
        let dir = kset_path.parent().unwrap_or(Path::new("."));
        kset::write_checked(kset_path, id, space, values, |_, need| {
            reserve.claim(dir, need)
        })?;
    }
    let genes = catalog.map_or(vec![], |cat| {
        catalog::family_calls(cat, &cat.detect(values, catalog::CALL_CONTAINMENT))
    });
    let c = Counted {
        genome_id: id.to_string(),
        n_kmers: values.len() as u64,
        stats,
        catalog_version: version,
        kset_digest: Some(kset::digest(values)),
        sketch: Some(values.iter().copied().take(SKETCH).collect()),
        genes,
    };
    write_sidecar(sidecar, &c, reserve)?;
    Ok(c)
}

pub fn count(
    layout: &Layout,
    slug: &str,
    k: usize,
    catalog: Option<&Catalog>,
) -> Result<CountSummary> {
    let previous = load_counted(layout, slug).ok();
    if let Some(c) = previous.as_ref().filter(|c| !c.pending.is_empty()) {
        anyhow::bail!(
            "{slug}: counted.json has {} selected genomes pending from an interrupted disk-sparse run; resume `grm pipeline --disk-sparse` rather than starting a stock count",
            c.pending.len()
        );
    }
    if catalog.is_none() {
        // Counting without a catalog records no gene calls, and learning
        // without them leaves every report without its catalog baseline.
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
                    0,
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
        pending: vec![],
    };
    confine::write_atomic(
        &layout.species_dir(slug).join("counted.json"),
        &serde_json::to_vec(&summary)?,
        false,
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
    structure_with_reserve(layout, slug, 0)
}

/// `structure` with a free-disk preflight for its JSON output.
pub fn structure_with_reserve(layout: &Layout, slug: &str, reserve: u64) -> Result<Structure> {
    let counted = load_counted(layout, slug)?;
    let kset_digests = counted.kset_digests("structure")?;
    let kdir = layout.kset_dir(slug);
    let sketches: Vec<Vec<u64>> = counted
        .genomes
        .par_iter()
        .map(|g| {
            let path = kdir.join(format!("{}.kset", g.genome_id));
            if std::fs::symlink_metadata(&path).is_ok() {
                let reader = KsetReader::open(&path)?;
                anyhow::ensure!(
                    reader.k == counted.k,
                    "{} holds {}-mers, not the counted {}-mers",
                    path.display(),
                    reader.k,
                    counted.k
                );
                return reader.sketch(SKETCH);
            }
            let sketch = g
                .sketch
                .as_ref()
                .with_context(|| format!("no kset or cached sketch for {}", g.genome_id))?;
            anyhow::ensure!(
                sketch.len() as u64 == g.n_kmers.min(SKETCH as u64),
                "counted.json has the wrong sketch length for {}",
                g.genome_id
            );
            Ok(sketch.clone())
        })
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
    let out = layout.species_dir(slug).join("structure.json");
    let bytes = serde_json::to_vec(&s)?;
    if reserve != 0 {
        matrix::ensure_free_space(
            out.parent().unwrap_or(Path::new(".")),
            bytes.len() as u64,
            reserve,
        )?;
    }
    confine::write_atomic(&out, &bytes, false)?;
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

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
struct SparseConfig {
    k: usize,
    genomes: Vec<String>,
    missing_fasta: Vec<String>,
    catalog_version: String,
    block_width: usize,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct ActiveBlock {
    first_col: usize,
    /// Kset paths that predated this block. Every other kset in the block is
    /// owned by this sparse run and may be deleted after partial verification.
    preexisting: Vec<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct SavedSparse {
    digest: String,
    summary: matrix::MatrixSummary,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct SparseRun {
    config: SparseConfig,
    folded_columns: usize,
    active: Option<ActiveBlock>,
    saved: Option<SavedSparse>,
}

fn ensure_single_component(value: &str, what: &str) -> Result<()> {
    let mut components = Path::new(value).components();
    anyhow::ensure!(
        matches!(components.next(), Some(std::path::Component::Normal(_)))
            && components.next().is_none(),
        "{what} is not one safe path component: {value:?}"
    );
    Ok(())
}

/// Atomic, durable JSON for the two manifests a sparse resume depends on,
/// written through a no-follow temporary (see `confine::write_atomic`).
fn write_json_durably<T: Serialize>(path: &Path, value: &T, reserve: u64) -> Result<()> {
    let bytes = serde_json::to_vec_pretty(value)?;
    let parent = path.parent().unwrap_or(Path::new("."));
    if reserve != 0 {
        matrix::ensure_free_space(parent, bytes.len() as u64, reserve)?;
    }
    confine::write_atomic(path, &bytes, true)
}

/// Removes a direct-child regular file, never a symlink, directory, or tree.
fn remove_owned_file(dir: &Path, name: &str) -> Result<bool> {
    let path = dir.join(name);
    anyhow::ensure!(
        path.parent() == Some(dir),
        "refusing non-child cleanup path {}",
        path.display()
    );
    match std::fs::symlink_metadata(&path) {
        Ok(metadata) => anyhow::ensure!(
            metadata.file_type().is_file(),
            "refusing to remove non-regular-file {}",
            path.display()
        ),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(false),
        Err(e) => return Err(e.into()),
    }
    std::fs::remove_file(&path).with_context(|| format!("remove {}", path.display()))?;
    Ok(true)
}

fn partial_path(run_dir: &Path, first_col: usize) -> PathBuf {
    run_dir.join(format!("block-{first_col:05}.partial"))
}

fn expected_run_files(config: &SparseConfig) -> HashSet<String> {
    let mut expected: HashSet<String> = ["run.json".into(), "run.json.part".into()].into();
    for first in (0..config.genomes.len()).step_by(config.block_width) {
        let partial = partial_path(Path::new(""), first);
        let name = partial.file_name().unwrap().to_string_lossy().into_owned();
        expected.insert(name.clone());
        expected.insert(format!("{name}.part"));
        expected.insert(format!("{name}.receipt.json"));
        expected.insert(format!("{name}.receipt.json.part"));
    }
    expected
}

fn validate_run_entries(run_dir: &Path, config: &SparseConfig) -> Result<()> {
    let expected = expected_run_files(config);
    for entry in std::fs::read_dir(run_dir)? {
        let entry = entry?;
        let name = entry.file_name().into_string().map_err(|_| {
            anyhow::anyhow!("non-UTF-8 filename in sparse run {}", run_dir.display())
        })?;
        anyhow::ensure!(
            expected.contains(&name),
            "unknown file {} in sparse run; refusing to touch or clean it",
            entry.path().display()
        );
        // A link or directory under a run file's name is never written
        // through or removed.
        anyhow::ensure!(
            entry.file_type()?.is_file(),
            "refusing non-regular entry {} in sparse run",
            entry.path().display()
        );
    }
    Ok(())
}

fn sparse_run_dir(layout: &Layout, slug: &str, config: &SparseConfig) -> Result<PathBuf> {
    let bytes = serde_json::to_vec(config)?;
    let digest = ring::digest::digest(&ring::digest::SHA256, &bytes);
    let hex = digest
        .as_ref()
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect::<String>();
    Ok(layout.species_dir(slug).join("matrix-partials").join(hex))
}

fn digest_bytes(ctx: &mut ring::digest::Context, name: &str, bytes: &[u8]) {
    ctx.update(name.as_bytes());
    ctx.update(&(bytes.len() as u64).to_le_bytes());
    ctx.update(bytes);
}

fn digest_u64s(ctx: &mut ring::digest::Context, name: &str, values: &[u64]) {
    ctx.update(name.as_bytes());
    ctx.update(&((values.len() * 8) as u64).to_le_bytes());
    let mut bytes = Vec::with_capacity(1 << 20);
    for chunk in values.chunks(1 << 15) {
        bytes.clear();
        for value in chunk {
            bytes.extend(value.to_le_bytes());
        }
        ctx.update(&bytes);
    }
}

fn digest_u32s(ctx: &mut ring::digest::Context, name: &str, values: &[u32]) {
    ctx.update(name.as_bytes());
    ctx.update(&((values.len() * 4) as u64).to_le_bytes());
    let mut bytes = Vec::with_capacity(1 << 20);
    for chunk in values.chunks(1 << 15) {
        bytes.clear();
        for value in chunk {
            bytes.extend(value.to_le_bytes());
        }
        ctx.update(&bytes);
    }
}

fn digest_hex(ctx: ring::digest::Context) -> String {
    ctx.finish()
        .as_ref()
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

fn sparse_matrix_digest(m: &matrix::FeatureMatrix) -> Result<String> {
    let mut ctx = ring::digest::Context::new(&ring::digest::SHA256);
    digest_bytes(
        &mut ctx,
        "matrix.json",
        &serde_json::to_vec_pretty(&m.summary)?,
    );
    digest_u64s(&mut ctx, "patterns.bin", &m.patterns.data);
    digest_u32s(&mut ctx, "support.bin", &m.support);
    digest_u64s(&mut ctx, "kmer_offsets.bin", &m.kmer_offsets);
    digest_u64s(&mut ctx, "kmers.bin", &m.kmers);
    Ok(digest_hex(ctx))
}

fn saved_matrix_digest(dir: &Path) -> Result<String> {
    let mut ctx = ring::digest::Context::new(&ring::digest::SHA256);
    for name in [
        "matrix.json",
        "patterns.bin",
        "support.bin",
        "kmer_offsets.bin",
        "kmers.bin",
    ] {
        let path = dir.join(name);
        let len = std::fs::metadata(&path)?.len();
        ctx.update(name.as_bytes());
        ctx.update(&len.to_le_bytes());
        let mut file = File::open(&path)?;
        let mut bytes = vec![0u8; 1 << 20];
        loop {
            let got = file.read(&mut bytes)?;
            if got == 0 {
                break;
            }
            ctx.update(&bytes[..got]);
        }
    }
    Ok(digest_hex(ctx))
}

fn verify_saved_sparse(dir: &Path, saved: &SavedSparse) -> Result<matrix::MatrixSummary> {
    anyhow::ensure!(
        saved_matrix_digest(dir)? == saved.digest,
        "saved matrix in {} does not match the digest recorded before partial cleanup",
        dir.display()
    );
    let loaded = matrix::FeatureMatrix::load(dir)?;
    anyhow::ensure!(
        serde_json::to_vec(&loaded.summary)? == serde_json::to_vec(&saved.summary)?,
        "saved matrix summary in {} does not match the sparse run",
        dir.display()
    );
    Ok(loaded.summary.clone())
}

fn cleanup_sparse_run(run_dir: &Path, run: &SparseRun) -> Result<()> {
    validate_run_entries(run_dir, &run.config)?;
    for first in (0..run.config.genomes.len()).step_by(run.config.block_width) {
        let partial = partial_path(run_dir, first);
        let name = partial.file_name().unwrap().to_string_lossy().into_owned();
        for suffix in ["", ".part", ".receipt.json", ".receipt.json.part"] {
            remove_owned_file(run_dir, &format!("{name}{suffix}"))?;
        }
    }
    remove_owned_file(run_dir, "run.json.part")?;
    remove_owned_file(run_dir, "run.json")?;
    sync_dir(run_dir)?;
    std::fs::remove_dir(run_dir)?;
    if let Some(parent) = run_dir.parent()
        && let Err(e) = std::fs::remove_dir(parent)
        && e.kind() != std::io::ErrorKind::DirectoryNotEmpty
        && e.kind() != std::io::ErrorKind::NotFound
    {
        return Err(e.into());
    }
    Ok(())
}

/// A sparse run directory's name: the SHA-256 of its config, in hex.
fn is_run_name(name: &str) -> bool {
    name.len() == 64
        && name
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}

/// Settles `matrix-partials` before a run starts. An empty run directory is
/// what a cleanup leaves when it stops between removing run.json and its
/// directory; it holds nothing, so it is removed (and `matrix-partials` with
/// it once empty). A non-empty directory of another run (another selection,
/// k, catalog or block width) holds the only copy of the blocks whose ksets
/// it deleted, so a new run is refused while it exists: it is neither
/// orphaned nor counted against the reserve as if it were free. Anything
/// that is not a run directory is refused too.
fn settle_other_runs(species_dir: &Path, partials_dir: &Path, run_dir: &Path) -> Result<()> {
    let entries = match std::fs::read_dir(partials_dir) {
        Ok(entries) => entries,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(()),
        Err(e) => return Err(e).with_context(|| format!("read {}", partials_dir.display())),
    };
    let mut removed = false;
    let mut left = 0;
    for entry in entries {
        let entry = entry?;
        let path = entry.path();
        let name = entry.file_name();
        let name = name.to_string_lossy();
        let metadata = std::fs::symlink_metadata(&path)?;
        anyhow::ensure!(
            metadata.is_dir() && is_run_name(&name),
            "unknown entry {} in {}; refusing to touch or clean it",
            path.display(),
            partials_dir.display()
        );
        if std::fs::read_dir(&path)?.next().is_none() {
            std::fs::remove_dir(&path).with_context(|| format!("remove {}", path.display()))?;
            removed = true;
            continue;
        }
        anyhow::ensure!(
            path.file_name() == run_dir.file_name(),
            "{} holds an unfinished disk-sparse run of another selection, k, catalog or block width, whose partials stand in for ksets it deleted; finish it by restoring those inputs, or review and remove that directory, before starting another run",
            path.display()
        );
        left += 1;
    }
    if removed {
        sync_dir(partials_dir)?;
    }
    if left == 0 && removed {
        std::fs::remove_dir(partials_dir)
            .with_context(|| format!("remove {}", partials_dir.display()))?;
        sync_dir(species_dir)?;
    }
    Ok(())
}

/// Refuses, before any work, a lane whose directories or files are links
/// or other non-regular entries where sparse writes go (see `confine`).
fn check_sparse_lane(
    lane: &confine::Lane,
    species_dir: &Path,
    kdir: &Path,
    matrix_dir: &Path,
    partials_dir: &Path,
    genomes: &[String],
) -> Result<()> {
    for name in [
        "counted.json",
        "counted.json.part",
        "structure.json",
        "structure.json.part",
    ] {
        lane.check_file(&species_dir.join(name))?;
    }
    if lane.check_dir_if_present(matrix_dir)? {
        for name in MATRIX_FILES {
            lane.check_file(&matrix_dir.join(name))?;
        }
        lane.check_file(&matrix_dir.join("matrix.json.part"))?;
    }
    lane.check_dir_if_present(partials_dir)?;
    if lane.check_dir_if_present(kdir)? {
        // The ksets themselves may be links to a shared store: they are only
        // read, and never deleted unless this run wrote them (regular files).
        for id in genomes {
            for name in [
                format!("{id}.json"),
                format!("{id}.json.part"),
                format!("{id}.kset.part"),
            ] {
                lane.check_file(&kdir.join(name))?;
            }
        }
    }
    Ok(())
}

/// Every existing kset of `genomes` outside `owned` opens at `k`.
fn check_kept_ksets(
    kdir: &Path,
    genomes: &[String],
    owned: &HashSet<&str>,
    k: usize,
) -> Result<()> {
    genomes
        .par_iter()
        .filter(|id| !owned.contains(id.as_str()))
        .try_for_each(|id| {
            let path = kdir.join(format!("{id}.kset"));
            if std::fs::symlink_metadata(&path).is_err() {
                return Ok(());
            }
            let reader = KsetReader::open(&path).with_context(|| {
                format!(
                    "read the kset {} that this disk-sparse run did not create; it never rewrites one, so move it aside or delete it to have it recounted",
                    path.display()
                )
            })?;
            anyhow::ensure!(
                reader.k == k,
                "{} holds {}-mers, not the {k}-mers this run counts; disk-sparse mode never rewrites a kset it did not create, so move it aside or delete it to have it recounted",
                path.display(),
                reader.k
            );
            Ok(())
        })
}

const MATRIX_FILES: [&str; 5] = [
    "matrix.json",
    "patterns.bin",
    "support.bin",
    "kmer_offsets.bin",
    "kmers.bin",
];

/// Syncs the saved matrix and the manifest that records it once more. A
/// crash after the saved run.json was renamed but before its directory was
/// synced leaves both on disk, yet not necessarily durable; the partials
/// they replace go only after this.
fn resync_saved(matrix_dir: &Path, run_path: &Path, run_dir: &Path) -> Result<()> {
    for name in MATRIX_FILES {
        let path = matrix_dir.join(name);
        File::open(&path)
            .and_then(|f| f.sync_all())
            .with_context(|| format!("sync {}", path.display()))?;
    }
    sync_dir(matrix_dir)?;
    File::open(run_path)
        .and_then(|f| f.sync_all())
        .with_context(|| format!("sync {}", run_path.display()))?;
    sync_dir(run_dir)
}

/// An upper bound on the k-mers a FASTA can hold: its decompressed length,
/// since every k-mer ends at a distinct base.
fn sequence_bytes(path: &Path) -> Result<u64> {
    let mut reader = crate::fasta::open(path)?;
    std::io::copy(&mut reader, &mut std::io::sink())
        .with_context(|| format!("read {}", path.display()))
}

/// The most new kset bytes a block's count can write: an existing file at
/// its length, a counted genome at the 8-bytes-per-k-mer ceiling of its
/// recorded k-mer count, and an uncounted one at that ceiling for as many
/// k-mers as its FASTA has bases (`kset::max_encoded_len`). An upper bound,
/// not an estimate: real MTB ksets of 26-30 MB come in under it.
fn block_kset_bound(
    kdir: &Path,
    fasta_dir: &Path,
    ids: &[String],
    counted: &[Counted],
) -> Result<u64> {
    ids.par_iter()
        .enumerate()
        .map(|(i, id)| {
            if let Ok(m) = std::fs::metadata(kdir.join(format!("{id}.kset"))) {
                return Ok(m.len());
            }
            let n_kmers = match counted.get(i) {
                Some(g) => g.n_kmers,
                None => sequence_bytes(&fasta_path(fasta_dir, id))?,
            };
            Ok(kset::max_encoded_len(id, n_kmers))
        })
        .try_reduce(|| 0, |a, b| Ok(a.saturating_add(b)))
}

fn count_sparse_block(
    layout: &Layout,
    slug: &str,
    ids: &[String],
    space: &KmerSpace,
    catalog: Option<&Catalog>,
    reserve: &Reserve,
    kept: &HashSet<&str>,
) -> Result<Vec<Counted>> {
    let fasta_dir = layout.fasta_dir();
    let kdir = layout.kset_dir(slug);
    let results: Vec<(String, Result<Counted>)> = ids
        .par_iter()
        .map_init(
            || (Vec::new(), Vec::new()),
            |buffers, id| {
                let result = count_one_with(
                    id,
                    &fasta_path(&fasta_dir, id),
                    &kdir.join(format!("{id}.kset")),
                    &kdir.join(format!("{id}.json")),
                    space,
                    catalog,
                    buffers,
                    reserve,
                    kept.contains(id.as_str()),
                );
                (id.clone(), result)
            },
        )
        .collect();
    results
        .into_iter()
        .map(|(id, result)| {
            result.with_context(|| format!("{slug}: count {id} for the current disk-sparse block"))
        })
        .collect()
}

/// Counts, folds, merges and saves a matrix while holding at most one new
/// block of ksets. Production always leaves `DISK_SPARSE_RESERVE` bytes free.
pub fn build_matrix_disk_sparse(
    layout: &Layout,
    slug: &str,
    k: usize,
    catalog: Option<&Catalog>,
    min_support: u32,
    block_width: usize,
) -> Result<matrix::MatrixSummary> {
    build_matrix_disk_sparse_with(
        layout,
        slug,
        k,
        catalog,
        min_support,
        block_width,
        DISK_SPARSE_RESERVE,
        |_| Ok(()),
    )
}

/// The block's expected columns and k-mer counts, from counted.json.
fn block_expectations(
    counted: &CountSummary,
    first: usize,
    end: usize,
) -> (Vec<matrix::Column>, Vec<u64>) {
    counted.genomes[first..end]
        .iter()
        .map(|g| {
            (
                matrix::Column {
                    name: g.genome_id.clone(),
                    kset_digest: g.kset_digest.clone().unwrap(),
                },
                g.n_kmers,
            )
        })
        .unzip()
}

#[allow(clippy::too_many_arguments)]
fn build_matrix_disk_sparse_with(
    layout: &Layout,
    slug: &str,
    k: usize,
    catalog: Option<&Catalog>,
    min_support: u32,
    block_width: usize,
    reserve: u64,
    mut hook: impl FnMut(&str) -> Result<()>,
) -> Result<matrix::MatrixSummary> {
    anyhow::ensure!(
        block_width >= 64
            && block_width.is_multiple_of(64)
            && block_width <= matrix::REVIEWED_COLUMNS,
        "disk-sparse block width must be a multiple of 64 in 64..={} not {block_width}",
        matrix::REVIEWED_COLUMNS
    );
    let space = KmerSpace::new(k);
    ensure_single_component(slug, "species slug")?;
    // Every write below stays in real directories under the data root: no
    // directory from the root down, and no file this run writes, may be a
    // link (see `confine`).
    let lane = confine::Lane::new(&layout.root)?;
    let species_dir = layout.species_dir(slug);
    lane.check_dir(&species_dir)?;
    let (present, missing) = selection(layout, slug)?;
    anyhow::ensure!(!present.is_empty(), "no selected genomes with FASTA files");
    for id in &present {
        ensure_single_component(id, "genome id")?;
    }
    let unique: HashSet<&str> = present.iter().map(String::as_str).collect();
    anyhow::ensure!(
        unique.len() == present.len(),
        "selected genome ids are not unique"
    );
    matrix::ensure_reviewed(present.len(), matrix::REVIEWED_COLUMNS)?;
    let config = SparseConfig {
        k,
        genomes: present,
        missing_fasta: missing,
        catalog_version: catalog.map_or(String::new(), Catalog::calls_version),
        block_width,
    };
    // An exclusive directory lock is released by the OS even on a crash; two
    // sparse callers cannot share and truncate one block's temporary files.
    let lock =
        File::open(&species_dir).with_context(|| format!("open {}", species_dir.display()))?;
    lock.lock()
        .with_context(|| format!("lock {}", species_dir.display()))?;
    let _lock = lock;
    let kdir = layout.kset_dir(slug);
    let matrix_dir = layout.matrix_dir(slug);
    let partials_dir = species_dir.join("matrix-partials");
    check_sparse_lane(
        &lane,
        &species_dir,
        &kdir,
        &matrix_dir,
        &partials_dir,
        &config.genomes,
    )?;
    // A crash before the first manifest rename can leave only our deterministic
    // temporary. No other writer holds the species lock.
    if remove_owned_file(&species_dir, "counted.json.part")? {
        sync_dir(&species_dir)?;
    }
    if reserve != 0 {
        matrix::ensure_free_space(&species_dir, 0, reserve)?;
    }
    let run_dir = sparse_run_dir(layout, slug, &config)?;
    settle_other_runs(&species_dir, &partials_dir, &run_dir)?;
    let run_path = run_dir.join("run.json");
    let mut run: Option<SparseRun> = if lane.check_file(&run_path)? {
        let run: SparseRun = serde_json::from_slice(&std::fs::read(&run_path)?)?;
        anyhow::ensure!(
            run.config == config,
            "sparse run config does not match its directory"
        );
        Some(run)
    } else {
        None
    };
    let previous = load_counted(layout, slug).ok();
    if catalog.is_none()
        && let Some(g) = previous
            .as_ref()
            .and_then(|c| c.called_with_other_catalog(""))
    {
        anyhow::bail!(
            "{slug}: counting without a catalog would drop gene calls from catalog {}; restore it with `grm fetch catalog`",
            g.catalog_version
        );
    }
    let compatible = previous.as_ref().is_some_and(|c| {
        c.k == k
            && c.failed.is_empty()
            && c.missing_fasta == config.missing_fasta
            && c.genomes.len() + c.pending.len() == config.genomes.len()
            && c.genomes
                .iter()
                .map(|g| &g.genome_id)
                .chain(c.pending.iter())
                .eq(config.genomes.iter())
            && c.genomes
                .iter()
                .all(|g| g.kset_digest.is_some() && g.catalog_version == config.catalog_version)
    });
    if run.is_some() {
        anyhow::ensure!(
            compatible,
            "counted.json no longer reconstructs the sparse run's exact selection, k, catalog and completed columns"
        );
    }
    let mut counted = if compatible {
        previous.unwrap()
    } else {
        CountSummary {
            k,
            genomes: vec![],
            missing_fasta: config.missing_fasta.clone(),
            failed: vec![],
            pending: config.genomes.clone(),
        }
    };

    // A current matrix needs no sparse work. A manifest, especially one in
    // the saved-but-not-cleaned phase, is still resumed below.
    if run.is_none()
        && counted.pending.is_empty()
        && counted.genomes.iter().all(|g| g.sketch.is_some())
        && let Ok(loaded) = matrix::FeatureMatrix::load(&matrix_dir)
        && loaded.summary.min_support == min_support
        && not_built_from(
            &counted,
            loaded.summary.k,
            &loaded.summary.genomes,
            &loaded.summary.kset_digests,
        )
        .is_none()
    {
        return Ok(loaded.summary.clone());
    }

    // Ksets this run did not create are inputs it never rewrites (see
    // `count_one_with`). One that does not read at this k is refused before
    // anything is written or recorded; the active block's own may be redone.
    let owned: HashSet<&str> = run
        .as_ref()
        .and_then(|r| r.active.as_ref())
        .map(|a| {
            let end = (a.first_col + block_width).min(config.genomes.len());
            config.genomes[a.first_col.min(end)..end]
                .iter()
                .map(String::as_str)
                .filter(|id| !a.preexisting.iter().any(|p| p == id))
                .collect()
        })
        .unwrap_or_default();
    check_kept_ksets(&kdir, &config.genomes, &owned, k)?;
    lane.create_dir_all(&run_dir, true)?;
    validate_run_entries(&run_dir, &config)?;
    if !compatible {
        write_json_durably(&species_dir.join("counted.json"), &counted, reserve)?;
    }
    if run.is_none() {
        // Only the manifest's own temporary is provably ours before the first
        // successful manifest rename. Any other file is left for diagnosis.
        let names: Vec<String> = std::fs::read_dir(&run_dir)?
            .map(|e| Ok(e?.file_name().to_string_lossy().into_owned()))
            .collect::<Result<_>>()?;
        anyhow::ensure!(
            names.iter().all(|name| name == "run.json.part"),
            "sparse run {} has files but no run.json; refusing to adopt or delete them",
            run_dir.display()
        );
        if !names.is_empty() {
            remove_owned_file(&run_dir, "run.json.part")?;
            sync_dir(&run_dir)?;
        }
        run = Some(SparseRun {
            config: config.clone(),
            folded_columns: 0,
            active: None,
            saved: None,
        });
        write_json_durably(&run_path, run.as_ref().unwrap(), reserve)?;
    }
    let mut run = run.take().unwrap();
    let n = config.genomes.len();
    let count_ceiling = run.folded_columns.saturating_add(block_width).min(n);
    // A run counts at most one block past its folded columns, unless the
    // whole selection was already counted: a stock `grm count`, or a
    // finished sparse run now rebuilt at another min_support. That count is
    // adopted block by block; ksets it no longer has are recounted, and
    // must reproduce its digests.
    let fully_counted = counted.pending.is_empty() && counted.genomes.len() == n;
    anyhow::ensure!(
        run.folded_columns <= n
            && (run.folded_columns == n || run.folded_columns.is_multiple_of(block_width))
            && run.folded_columns <= counted.genomes.len()
            && (counted.genomes.len() <= count_ceiling || fully_counted),
        "sparse run and counted.json disagree about completed columns"
    );
    if let Some(active) = &run.active {
        anyhow::ensure!(
            active.first_col == run.folded_columns && active.first_col < n,
            "sparse run has an invalid active block start"
        );
        let end = (active.first_col + block_width).min(n);
        let unique: HashSet<&str> = active.preexisting.iter().map(String::as_str).collect();
        anyhow::ensure!(
            unique.len() == active.preexisting.len()
                && active
                    .preexisting
                    .iter()
                    .all(|id| { config.genomes[active.first_col..end].contains(id) }),
            "sparse run has an invalid active block"
        );
    }

    if let Some(saved) = &run.saved {
        anyhow::ensure!(
            run.folded_columns == n && counted.pending.is_empty() && run.active.is_none(),
            "saved sparse run has incomplete counting or folding"
        );
        // min_support is not part of the run's config: blocks do not depend
        // on it. A saved matrix at another min_support is finished first.
        anyhow::ensure!(
            saved.summary.min_support == min_support,
            "the interrupted disk-sparse run already saved a matrix at min_support {}; rerun it with --min-support {} to finish its cleanup, then build at min_support {min_support}",
            saved.summary.min_support,
            saved.summary.min_support
        );
        let summary = verify_saved_sparse(&matrix_dir, saved)?;
        resync_saved(&matrix_dir, &run_path, &run_dir)?;
        cleanup_sparse_run(&run_dir, &run)?;
        return Ok(summary);
    }

    // Folded blocks have already lost their owned ksets. A bad one is an
    // integrity error, not something to merge or silently replace.
    for first in (0..run.folded_columns).step_by(block_width) {
        let end = (first + block_width).min(n);
        let (columns, n_kmers) = block_expectations(&counted, first, end);
        matrix::verify_partial(&partial_path(&run_dir, first), &columns, &n_kmers, first, k)?;
    }

    lane.create_dir_all(&kdir, false)?;
    let ledger = Reserve::new(reserve);
    while run.folded_columns < n {
        let first = run.folded_columns;
        let end = (first + block_width).min(n);
        let ids = &config.genomes[first..end];
        if run.active.is_none() {
            let preexisting = ids
                .iter()
                .filter(|id| std::fs::symlink_metadata(kdir.join(format!("{id}.kset"))).is_ok())
                .cloned()
                .collect();
            run.active = Some(ActiveBlock {
                first_col: first,
                preexisting,
            });
            write_json_durably(&run_path, &run, reserve)?;
        }
        let active = run.active.clone().unwrap();
        let preexisting: HashSet<&str> = active.preexisting.iter().map(String::as_str).collect();
        let partial = partial_path(&run_dir, first);
        let partial_name = partial.file_name().unwrap().to_string_lossy().into_owned();
        let verified = counted.genomes.len() >= end && {
            let (columns, n_kmers) = block_expectations(&counted, first, end);
            matrix::verify_partial(&partial, &columns, &n_kmers, first, k).is_ok()
        };
        if !verified {
            let needs_count = counted.genomes.len() < end
                || counted.genomes[first..counted.genomes.len().min(end)]
                    .iter()
                    .any(|g| g.sketch.is_none())
                || ids
                    .iter()
                    .any(|id| !std::fs::symlink_metadata(kdir.join(format!("{id}.kset"))).is_ok());
            if needs_count {
                if reserve != 0 {
                    let bound = block_kset_bound(
                        &kdir,
                        &layout.fasta_dir(),
                        ids,
                        counted.genomes.get(first..).unwrap_or(&[]),
                    )?;
                    matrix::ensure_free_space(&kdir, bound, reserve).with_context(|| {
                        format!(
                            "preflight the current block's ksets at their {bound}-byte upper bound (each write also claims its exact length against the reserve)"
                        )
                    })?;
                }
                hook(&format!("before_count:{first}"))?;
                lane.check_dir(&kdir)?;
                let new =
                    count_sparse_block(layout, slug, ids, &space, catalog, &ledger, &preexisting)?;
                // Recounting to recreate a missing cache must reproduce every
                // existing expectation; it may not silently change a matrix.
                for (i, c) in new.iter().enumerate() {
                    if let Some(old) = counted.genomes.get(first + i) {
                        anyhow::ensure!(
                            old.genome_id == c.genome_id
                                && old.kset_digest == c.kset_digest
                                && old.catalog_version == c.catalog_version,
                            "{} changed genome, kset digest or catalog while recreating a sparse block; start a new run from the changed inputs instead of rewriting expectations",
                            c.genome_id
                        );
                    }
                }
                counted.genomes.truncate(first);
                counted.genomes.extend(new);
                counted.pending = config.genomes[end..].to_vec();
                write_json_durably(&species_dir.join("counted.json"), &counted, reserve)?;
            }
            for suffix in ["", ".part", ".receipt.json", ".receipt.json.part"] {
                remove_owned_file(&run_dir, &format!("{partial_name}{suffix}"))?;
            }
            let paths: Vec<PathBuf> = ids
                .iter()
                .map(|id| kdir.join(format!("{id}.kset")))
                .collect();
            lane.check_dir(&run_dir)?;
            matrix::write_partial(
                &paths,
                &block_expectations(&counted, first, end).0,
                first,
                &partial,
                rayon::current_num_threads() * 2,
                reserve,
            )?;
        }
        // Re-open, hash and read back the pair even when write_partial just
        // returned: every column must hold its genome's counted k-mers. No
        // owned kset is deleted before this succeeds.
        let (columns, n_kmers) = block_expectations(&counted, first, end);
        matrix::verify_partial(&partial, &columns, &n_kmers, first, k)?;
        for suffix in [".part", ".receipt.json.part"] {
            remove_owned_file(&run_dir, &format!("{partial_name}{suffix}"))?;
        }
        sync_dir(&run_dir)?;
        hook(&format!("after_partial:{first}"))?;
        lane.check_dir(&kdir)?;
        for id in ids.iter().filter(|id| !preexisting.contains(id.as_str())) {
            remove_owned_file(&kdir, &format!("{id}.kset"))?;
        }
        sync_dir(&kdir)?;
        run.folded_columns = end;
        run.active = None;
        write_json_durably(&run_path, &run, reserve)?;
        hook(&format!("after_fold:{first}"))?;
        eprintln!("  {slug}: folded {end} / {n} genomes");
    }
    anyhow::ensure!(counted.pending.is_empty());

    let columns: Vec<matrix::Column> = counted
        .genomes
        .iter()
        .map(|g| matrix::Column {
            name: g.genome_id.clone(),
            kset_digest: g.kset_digest.clone().unwrap(),
        })
        .collect();
    let partials: Vec<PathBuf> = (0..n)
        .step_by(block_width)
        .map(|first| partial_path(&run_dir, first))
        .collect();
    let mut m = matrix::build_from_partials(
        &partials,
        &columns,
        k,
        min_support,
        rayon::current_num_threads() * 2,
    )?;
    m.summary.kset_digests = counted.kset_digests("matrix")?;
    let saved = SavedSparse {
        digest: sparse_matrix_digest(&m)?,
        summary: m.summary.clone(),
    };
    lane.create_dir_all(&matrix_dir, true)?;
    m.save_with_reserve(&matrix_dir, reserve)?;
    drop(m);
    let summary = verify_saved_sparse(&matrix_dir, &saved)?;
    run.saved = Some(saved);
    write_json_durably(&run_path, &run, reserve)?;
    hook("after_saved")?;
    cleanup_sparse_run(&run_dir, &run)?;
    Ok(summary)
}

/// Stages to rerun before learning so that the gene calls, structure, and
/// matrix all come from the selected genomes, their k, and the catalog.
pub struct Stale {
    pub count: bool,
    pub structure: bool,
    pub matrix: bool,
}

/// Compares counted.json with genomes.txt, `catalog`, and the pipeline's `k`,
/// structure.json and matrix.json each with counted.json, and matrix.json
/// with the pipeline's `min_support`. A catalog that cannot be read (None) is
/// no reason to recount: counting without one records no calls.
pub fn stale_stages(
    layout: &Layout,
    slug: &str,
    catalog: Option<&Catalog>,
    k: usize,
    min_support: u32,
) -> Stale {
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
        k,
        min_support,
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
    /// matrix.json only.
    #[serde(default)]
    min_support: Option<u32>,
}

/// Structure and matrix are each stale when the file is missing or
/// `not_built_from` the counts. Counts are stale when they list other
/// genomes than `count` would count now: the genomes in genomes.txt with a
/// FASTA file (`with_fasta`, None when genomes.txt cannot be read), less
/// those whose count failed, in file order. They are also stale when a
/// genome's gene calls carry another catalog version than `calls_version`,
/// the loaded catalog's (None when none loaded), or when counted at another k
/// than the pipeline's; the matrix is also stale when built with another
/// min_support. Stages built by hand with `grm count --k` or `grm matrix
/// --min-support` agree with each other, and these rules rebuild them. Count
/// also reruns when structure or matrix is stale, because their rebuild reads
/// the ksets, which `--drop-ksets` deletes and a recount writes again (reusing
/// those still on disk). That covers counts recorded before kset digests,
/// which no structure or matrix matches.
fn selection_matches(counted: &CountSummary, with_fasta: &[String]) -> bool {
    if !counted.pending.is_empty() {
        return false;
    }
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
    k: usize,
    min_support: u32,
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
    let other_min_support = matrix.is_some_and(|m| m.min_support != Some(min_support));
    let (structure, matrix) = (outdated(structure), outdated(matrix) || other_min_support);
    let other_catalog =
        calls_version.is_some_and(|v| counted.called_with_other_catalog(v).is_some());
    Stale {
        count: other_genomes || counted.k != k || other_catalog || structure || matrix,
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
            count_one(
                "g",
                &fasta,
                &kset_path,
                &json,
                &space,
                None,
                &mut buffers,
                0,
            )
            .unwrap()
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
            count_one(
                "g",
                &fasta,
                &kset_path,
                &json,
                &space,
                None,
                &mut buffers,
                0,
            )
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
            count_one(
                "g",
                &fasta,
                &kset_path,
                &json,
                &space,
                None,
                &mut buffers,
                0,
            )
            .unwrap()
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
            count_one(id, &fasta, &kset_path, &json, &space, None, &mut buffers, 0).unwrap()
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
        let mut count = || {
            count_one(
                "g",
                &fasta,
                &kset_path,
                &json,
                &space,
                None,
                &mut buffers,
                0,
            )
        };
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
                    sketch: Some(vec![]),
                    genes: vec![],
                })
                .collect(),
            missing_fasta: vec![],
            failed: vec![],
            pending: vec![],
        }
    }

    /// What a structure or matrix built at `k` for `genomes` records, with
    /// the kset digests `counted` gives them and the pipeline's min_support.
    fn built(k: usize, genomes: &[String]) -> BuiltFrom {
        BuiltFrom {
            k,
            genomes: genomes.to_vec(),
            kset_digests: genomes.iter().map(|g| format!("kset of {g}")).collect(),
            min_support: Some(2),
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
                31,
                2,
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
                31,
                2,
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
        // Stages built by hand at k = 21 agree with each other but not with
        // the pipeline's k, so the genomes are counted again; structure and
        // matrix are then compared with those counts.
        let by_hand = || Some(built(21, &abc));
        assert_eq!(
            plan(Some(&at_21), by_hand(), by_hand()),
            (true, false, false)
        );
        // A matrix built by hand with another min_support keeps other k-mers.
        let min_support_1 = BuiltFrom {
            min_support: Some(1),
            ..built(31, &abc)
        };
        assert_eq!(
            plan(Some(&at_31), current(), Some(min_support_1)),
            (true, false, true)
        );
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
                31,
                2,
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
        std::fs::write(
            catalog_dir.join("ReferenceGeneCatalog.txt"),
            "type\tclass\tsubclass\trefseq_protein_accession\tgenbank_protein_accession\tgene_family\tproduct_name\n\
             AMR\tfixture\tfixture\tREF\t\tfixture\tfixture\n",
        ).unwrap();
        std::fs::write(
            catalog_dir.join("AMR_CDS.fa"),
            format!(">REF\n{}\n", "ACGT".repeat(16)),
        )
        .unwrap();
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

#[cfg(test)]
mod sparse_tests;
