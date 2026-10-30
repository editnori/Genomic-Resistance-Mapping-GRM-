//! Species identity: genomes whose sequence belongs to another species than
//! their BV-BRC label.
//!
//! BV-BRC species names come from submitters. CheckM scores a genome against
//! the marker set of the lineage it was filed under, so a clean assembly of
//! the wrong species can still score well. This check reads the sequence
//! itself. Each genome's bottom-s MinHash sketch is compared with a reference
//! genome for its labelled species, and its GC content with the species' GC
//! distribution.
//!
//! The reference is the species medoid: the genome whose median Mash
//! distance to the rest of the species' genomes is smallest. It comes from
//! the genomes themselves, so no external download is needed. A median is
//! unaffected by a few foreign genomes, so contaminants cannot become the
//! reference unless they are most of the species.
//!
//! Mash distance approximates 1 - ANI near the species boundary of about 95%
//! ANI (distance 0.05). Species complexes (E. cloacae, K. oxytoca, E. faecium
//! with E. lactis, P. aeruginosa with P. paraeruginosa) spread past that
//! boundary to distance 0.12, so a distance past it is not by itself clear
//! evidence. A genome fails when
//! - its distance to its own reference exceeds `fail_distance`, or
//! - another checked species' reference is nearer than its own and within
//!   that species' boundary, so the sequence names that species.
//!
//! The second rule never applies between two members of a group in
//! `genomic_groups`. E. coli and every Shigella are one genomic species split
//! by taxonomy: hundreds of correctly labelled E. coli are nearer the
//! S. sonnei reference than their own, and S. flexneri genomes sit 0.015 from
//! the S. sonnei reference. A genome the rule would fail within its group is
//! borderline instead. Genomes are never moved to another species.
//!
//! A genome is also borderline when it fails neither way but is outside its
//! own boundary or has GC far off the species median. Borderline genomes are
//! recorded and never dropped. GC alone does not fail a genome, since
//! contaminating contigs and plasmids also move it.

use crate::bvbrc::fasta_path;
use crate::cluster;
use crate::dna::KmerSpace;
use crate::fasta::{self, StatsBuilder};
use crate::layout::{Layout, ensure_parent};
use anyhow::{Context, Result};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// Mash's default k: at 95% ANI, two genomes still share about a fifth of
/// their 21-mers, so the sketch resolves the species boundary.
pub const K: usize = 21;
pub const SKETCH: usize = crate::species::SKETCH;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CheckConfig {
    /// Mash distance to the reference past which a genome is another species
    /// outright (about 87% ANI).
    pub fail_distance: f64,
    /// The species boundary: about 95% ANI.
    pub species_distance: f64,
    /// GC deviation from the species median, in robust standard deviations
    /// (1.4826 x MAD), past which a genome is borderline...
    pub gc_z: f64,
    /// ...when it is also at least this many percentage points off.
    pub gc_points: f64,
    /// Groups of species the nearer-species rule never applies between, in
    /// any direction. A member "Genus *" stands for every species of the
    /// genus.
    pub genomic_groups: Vec<Vec<String>>,
}

impl CheckConfig {
    /// Whether two species share a genomic group.
    pub fn same_group(&self, a: &str, b: &str) -> bool {
        let member = |pattern: &str, species: &str| match pattern.strip_suffix(" *") {
            Some(genus) => species.split(' ').next() == Some(genus),
            None => pattern == species,
        };
        self.genomic_groups
            .iter()
            .any(|group| group.iter().any(|p| member(p, a)) && group.iter().any(|p| member(p, b)))
    }
}

impl Default for CheckConfig {
    fn default() -> Self {
        CheckConfig {
            // Across the 23 species of revision 9, species complexes reach
            // 0.111 (E. cloacae) and the next genome is at 0.146, a
            // Klebsiella filed as E. coli. 0.13 falls in that gap.
            fail_distance: 0.13,
            species_distance: 0.05,
            gc_z: 6.0,
            gc_points: 1.0,
            // Admin rulings of 2026-10-29.
            genomic_groups: vec![vec!["Escherichia coli".into(), "Shigella *".into()]],
        }
    }
}

/// A genome's GC and MinHash sketch, cached per genome so a later check only
/// sketches newly downloaded assemblies.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GenomeSketch {
    pub genome_id: String,
    pub k: usize,
    pub sketch_size: usize,
    pub length: u64,
    pub gc: f64,
    pub sketch: Vec<u64>,
}

/// GC and the `s` smallest mixed canonical k-mers of an assembly: the same
/// sketch `KsetReader::sketch` reads from a kset counted at `space.k`.
pub fn sketch_fasta(id: &str, path: &Path, space: &KmerSpace, s: usize) -> Result<GenomeSketch> {
    let mut stats = StatsBuilder::default();
    let mut kept: Vec<u64> = Vec::with_capacity(4 * s);
    let mut cut = u64::MAX;
    fasta::for_each_record(path, |_, seq| {
        stats.add(seq);
        space.for_each_canonical(seq, |x| {
            let h = space.mix(x);
            if h < cut {
                kept.push(h);
                if kept.len() >= 4 * s {
                    kept.sort_unstable();
                    kept.dedup();
                    kept.truncate(s);
                    if kept.len() == s {
                        cut = kept[s - 1];
                    }
                }
            }
        });
    })?;
    kept.sort_unstable();
    kept.dedup();
    kept.truncate(s);
    let stats = stats.finish();
    Ok(GenomeSketch {
        genome_id: id.to_string(),
        k: space.k,
        sketch_size: s,
        length: stats.length,
        gc: stats.gc,
        sketch: kept,
    })
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Verdict {
    Pass,
    Borderline,
    Fail,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GenomeCheck {
    pub genome_id: String,
    /// The labelled species.
    pub species: String,
    pub verdict: Verdict,
    /// GC fraction, and its deviation from the species median in robust
    /// standard deviations.
    pub gc: f64,
    pub gc_z: f64,
    /// Mash distance to the labelled species' reference, and 1 minus it as
    /// an ANI estimate.
    pub distance: f64,
    pub ani: f64,
    /// The other checked species whose reference is closest, and its
    /// distance.
    pub nearest_species: Option<String>,
    pub nearest_distance: Option<f64>,
    /// The species whose reference is nearer than the genome's own and
    /// within its boundary; one outside the genome's genomic group first.
    pub inside_species: Option<String>,
    pub reason: Option<String>,
}

/// Genomes of one labelled species inside another species' boundary and
/// nearer its reference: what the nearer-species rule fails, or records as
/// borderline within a genomic group.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PairCount {
    pub labelled: String,
    pub nearer: String,
    pub genomes: usize,
    /// Share of the labelled species' checked genomes.
    pub share: f64,
    /// The two share a genomic group, so none of these genomes fail.
    pub same_group: bool,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SpeciesSummary {
    pub species: String,
    pub slug: String,
    /// The medoid genome every distance is measured to.
    pub reference: String,
    pub gc_median: f64,
    pub gc_sd: f64,
    pub checked: usize,
    pub failed: usize,
    pub borderline: usize,
    /// Listed genomes with neither a FASTA file nor a cached sketch.
    pub unchecked: Vec<String>,
    /// Distance to the reference at the 50th, 99th and 100th percentile of
    /// the genomes that pass.
    pub distance_p50: f64,
    pub distance_p99: f64,
    pub distance_max_pass: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SpeciesCheck {
    pub k: usize,
    pub sketch_size: usize,
    pub config: CheckConfig,
    pub species: Vec<SpeciesSummary>,
    pub pairs: Vec<PairCount>,
    /// Every checked genome, by species then genome id.
    pub genomes: Vec<GenomeCheck>,
}

impl SpeciesCheck {
    /// Failed genome id -> why, the genomes `grm curate` excludes.
    pub fn failures(&self) -> BTreeMap<String, String> {
        self.genomes
            .iter()
            .filter(|g| g.verdict == Verdict::Fail)
            .map(|g| (g.genome_id.clone(), g.reason.clone().unwrap_or_default()))
            .collect()
    }
}

/// One labelled species and its genomes' sketches.
pub struct Cohort {
    pub species: String,
    pub slug: String,
    pub sketches: Vec<GenomeSketch>,
    pub unchecked: Vec<String>,
}

fn median(sorted: &[f64]) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    sorted[sorted.len() / 2]
}

fn quantile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    sorted[((sorted.len() - 1) as f64 * q).round() as usize]
}

fn distance(a: &GenomeSketch, b: &GenomeSketch) -> f64 {
    cluster::mash_distance(cluster::jaccard(&a.sketch, &b.sketch, a.sketch_size), a.k)
}

/// The genome with the smallest median distance to the others; ties go to
/// the earlier genome, and cohorts are sorted by id, so the choice is
/// reproducible.
fn medoid(sketches: &[GenomeSketch]) -> usize {
    let n = sketches.len();
    if n <= 2 {
        return 0;
    }
    let s: Vec<Vec<u64>> = sketches.iter().map(|g| g.sketch.clone()).collect();
    let dist = cluster::distance_matrix(&s, sketches[0].sketch_size, sketches[0].k);
    let medians: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut row: Vec<f64> = (0..n)
                .filter(|&j| j != i)
                .map(|j| dist[i * n + j] as f64)
                .collect();
            row.sort_by(f64::total_cmp);
            median(&row)
        })
        .collect();
    (0..n)
        .min_by(|&a, &b| medians[a].total_cmp(&medians[b]).then(a.cmp(&b)))
        .unwrap()
}

/// Checks every cohort against its own reference and the other cohorts'.
pub fn check(mut cohorts: Vec<Cohort>, cfg: &CheckConfig) -> SpeciesCheck {
    cohorts.retain(|c| !c.sketches.is_empty());
    cohorts.sort_by(|a, b| a.species.cmp(&b.species));
    for c in &mut cohorts {
        c.sketches.sort_by(|a, b| a.genome_id.cmp(&b.genome_id));
        c.unchecked.sort();
    }
    let refs: Vec<&GenomeSketch> = cohorts
        .iter()
        .map(|c| &c.sketches[medoid(&c.sketches)])
        .collect();
    let mut species = Vec::new();
    let mut pairs = Vec::new();
    let mut genomes = Vec::new();
    for (ci, c) in cohorts.iter().enumerate() {
        let mut gcs: Vec<f64> = c.sketches.iter().map(|g| g.gc).collect();
        gcs.sort_by(f64::total_cmp);
        let gc_median = median(&gcs);
        let mut dev: Vec<f64> = gcs.iter().map(|g| (g - gc_median).abs()).collect();
        dev.sort_by(f64::total_cmp);
        // A floor keeps a species of near-identical genomes from calling a
        // 0.1-point difference many standard deviations.
        let gc_sd = (1.4826 * median(&dev)).max(0.0005);
        let checked: Vec<GenomeCheck> = c
            .sketches
            .par_iter()
            .map(|g| {
                let d = distance(g, refs[ci]);
                let nearest = (0..cohorts.len())
                    .filter(|&o| o != ci)
                    .map(|o| (o, distance(g, refs[o])))
                    .min_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
                let gc_z = (g.gc - gc_median) / gc_sd;
                let gc_off = gc_z.abs() > cfg.gc_z && (g.gc - gc_median).abs() * 100.0 >= cfg.gc_points;
                let evidence = format!(
                    "Mash distance {d:.4} (ANI ~{:.1}%) to the {} reference {}, GC {:.1}% against a species median of {:.1}%",
                    (1.0 - d) * 100.0,
                    c.species,
                    refs[ci].genome_id,
                    g.gc * 100.0,
                    gc_median * 100.0
                );
                let inside = |grouped: bool| {
                    (0..cohorts.len())
                        .filter(|&o| o != ci && cfg.same_group(&c.species, &cohorts[o].species) == grouped)
                        .map(|o| (o, distance(g, refs[o])))
                        .min_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)))
                        .filter(|&(_, od)| od <= cfg.species_distance && od < d)
                };
                let inside_other = inside(false).or_else(|| inside(true));
                let (verdict, why) = if d > cfg.fail_distance {
                    (Verdict::Fail, format!("past the {} fail distance", cfg.fail_distance))
                } else if let Some((o, od)) = inside_other {
                    let other = &cohorts[o].species;
                    let why = format!(
                        "nearer the {other} reference {} ({od:.4}) and within its {} boundary",
                        refs[o].genome_id, cfg.species_distance
                    );
                    if cfg.same_group(&c.species, other) {
                        (Verdict::Borderline, format!("{why}; exempt, one genomic group with {}", c.species))
                    } else {
                        (Verdict::Fail, why)
                    }
                } else if d > cfg.species_distance {
                    (Verdict::Borderline, format!("past the {} species boundary", cfg.species_distance))
                } else if gc_off {
                    (Verdict::Borderline, format!("GC {gc_z:+.1} SD off"))
                } else {
                    (Verdict::Pass, String::new())
                };
                let reason = (verdict != Verdict::Pass).then(|| format!("{evidence}; {why}"));
                GenomeCheck {
                    genome_id: g.genome_id.clone(),
                    species: c.species.clone(),
                    verdict,
                    gc: g.gc,
                    gc_z,
                    distance: d,
                    ani: 1.0 - d,
                    nearest_species: nearest.map(|(o, _)| cohorts[o].species.clone()),
                    nearest_distance: nearest.map(|(_, od)| od),
                    inside_species: inside_other.map(|(o, _)| cohorts[o].species.clone()),
                    reason,
                }
            })
            .collect();
        let mut pass: Vec<f64> = checked
            .iter()
            .filter(|g| g.verdict == Verdict::Pass)
            .map(|g| g.distance)
            .collect();
        pass.sort_by(f64::total_cmp);
        let count = |v: Verdict| checked.iter().filter(|g| g.verdict == v).count();
        species.push(SpeciesSummary {
            species: c.species.clone(),
            slug: c.slug.clone(),
            reference: refs[ci].genome_id.clone(),
            gc_median,
            gc_sd,
            checked: checked.len(),
            failed: count(Verdict::Fail),
            borderline: count(Verdict::Borderline),
            unchecked: c.unchecked.clone(),
            distance_p50: quantile(&pass, 0.5),
            distance_p99: quantile(&pass, 0.99),
            distance_max_pass: pass.last().copied().unwrap_or(f64::NAN),
        });
        let mut inside: BTreeMap<&str, usize> = BTreeMap::new();
        for g in &checked {
            if let Some(other) = &g.inside_species
                && g.distance <= cfg.fail_distance
            {
                *inside.entry(other).or_default() += 1;
            }
        }
        for (other, n) in inside {
            pairs.push(PairCount {
                labelled: c.species.clone(),
                nearer: other.to_string(),
                genomes: n,
                share: n as f64 / checked.len() as f64,
                same_group: cfg.same_group(&c.species, other),
            });
        }
        genomes.extend(checked);
    }
    SpeciesCheck {
        k: refs.first().map_or(K, |r| r.k),
        sketch_size: refs.first().map_or(SKETCH, |r| r.sketch_size),
        config: cfg.clone(),
        species,
        pairs,
        genomes,
    }
}

pub fn sketch_path(layout: &Layout, genome_id: &str) -> PathBuf {
    layout
        .species_check_dir()
        .join("sketches")
        .join(format!("{genome_id}.json"))
}

/// The cached sketch when it was taken at `K` and `SKETCH`, else a new one
/// from the FASTA file, else None when neither exists.
pub fn load_or_sketch(layout: &Layout, genome_id: &str) -> Result<Option<GenomeSketch>> {
    let cache = sketch_path(layout, genome_id);
    if let Ok(bytes) = std::fs::read(&cache)
        && let Ok(s) = serde_json::from_slice::<GenomeSketch>(&bytes)
        && s.k == K
        && s.sketch_size == SKETCH
        && s.genome_id == genome_id
    {
        return Ok(Some(s));
    }
    let fasta = fasta_path(&layout.fasta_dir(), genome_id);
    if !fasta.exists() {
        return Ok(None);
    }
    let s = sketch_fasta(genome_id, &fasta, &KmerSpace::new(K), SKETCH)
        .with_context(|| format!("sketch {}", fasta.display()))?;
    ensure_parent(&cache)?;
    let tmp = cache.with_extension("json.tmp");
    std::fs::write(&tmp, serde_json::to_vec(&s)?)?;
    std::fs::rename(&tmp, &cache)?;
    Ok(Some(s))
}

/// Sketches (or reads cached sketches of) each species' listed genomes.
pub fn load_cohorts(
    layout: &Layout,
    species: &[(String, String, Vec<String>)],
) -> Result<Vec<Cohort>> {
    species
        .iter()
        .map(|(name, slug, ids)| {
            let found: Vec<(String, Option<GenomeSketch>)> = ids
                .par_iter()
                .map(|id| Ok((id.clone(), load_or_sketch(layout, id)?)))
                .collect::<Result<_>>()?;
            let mut cohort = Cohort {
                species: name.clone(),
                slug: slug.clone(),
                sketches: vec![],
                unchecked: vec![],
            };
            for (id, s) in found {
                match s {
                    Some(s) => cohort.sketches.push(s),
                    None => cohort.unchecked.push(id),
                }
            }
            Ok(cohort)
        })
        .collect()
}

pub fn write(layout: &Layout, check: &SpeciesCheck) -> Result<PathBuf> {
    let path = layout.species_check();
    ensure_parent(&path)?;
    std::fs::write(&path, serde_json::to_vec_pretty(check)?)?;
    Ok(path)
}

pub fn read(path: &Path) -> Result<SpeciesCheck> {
    serde_json::from_slice(
        &std::fs::read(path).with_context(|| format!("read {}", path.display()))?,
    )
    .with_context(|| format!("parse {}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    const LEN: usize = 200_000;

    fn random_genome(rng: &mut ChaCha8Rng, gc: f64) -> Vec<u8> {
        (0..LEN)
            .map(|_| match (rng.random_bool(gc), rng.random_bool(0.5)) {
                (true, true) => b'G',
                (true, false) => b'C',
                (false, true) => b'A',
                (false, false) => b'T',
            })
            .collect()
    }

    /// Substitutes each base with probability `rate`: about `rate` Mash
    /// distance, and 1 - `rate` ANI, from `seq`.
    fn mutate(rng: &mut ChaCha8Rng, seq: &[u8], rate: f64) -> Vec<u8> {
        seq.iter()
            .map(|&b| {
                if rng.random_bool(rate) {
                    let other: Vec<u8> = b"ACGT".iter().copied().filter(|&c| c != b).collect();
                    other[rng.random_range(0..3)]
                } else {
                    b
                }
            })
            .collect()
    }

    fn sketch(dir: &Path, id: &str, seq: &[u8]) -> GenomeSketch {
        let path = dir.join(format!("{id}.fna"));
        let mut text = format!(">{id}\n").into_bytes();
        text.extend_from_slice(seq);
        text.push(b'\n');
        std::fs::write(&path, text).unwrap();
        sketch_fasta(id, &path, &KmerSpace::new(K), SKETCH).unwrap()
    }

    /// `n` strains of a species around `root`, each 0.5% from it, so strains
    /// sit about 1% apart, as within a clonal species.
    fn strains(
        rng: &mut ChaCha8Rng,
        dir: &Path,
        prefix: &str,
        root: &[u8],
        n: usize,
    ) -> Vec<GenomeSketch> {
        (0..n)
            .map(|i| sketch(dir, &format!("{prefix}.{i}"), &mutate(rng, root, 0.005)))
            .collect()
    }

    fn cohort(species: &str, sketches: Vec<GenomeSketch>) -> Cohort {
        Cohort {
            species: species.into(),
            slug: species.to_lowercase().replace(' ', "_"),
            sketches,
            unchecked: vec![],
        }
    }

    fn verdicts(check: &SpeciesCheck, v: Verdict) -> Vec<(&str, &str)> {
        check
            .genomes
            .iter()
            .filter(|g| g.verdict == v)
            .map(|g| (g.genome_id.as_str(), g.species.as_str()))
            .collect()
    }

    #[test]
    fn the_sketch_is_the_bottom_s_of_every_mixed_canonical_kmer() {
        let dir = tempfile::tempdir().unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(1);
        let seq = random_genome(&mut rng, 0.6);
        let s = sketch(dir.path(), "g", &seq);
        let space = KmerSpace::new(K);
        let mut all = Vec::new();
        space.for_each_canonical(&seq, |x| all.push(space.mix(x)));
        all.sort_unstable();
        all.dedup();
        all.truncate(SKETCH);
        assert_eq!(s.sketch, all);
        assert_eq!(s.length, LEN as u64);
        assert!((s.gc - 0.6).abs() < 0.01, "{}", s.gc);
    }

    #[test]
    fn a_mixed_species_cohort_fails_exactly_its_foreign_genomes() {
        let dir = tempfile::tempdir().unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(2);
        let alpha_root = random_genome(&mut rng, 0.65);
        let gamma_root = random_genome(&mut rng, 0.40);
        let mut alpha = strains(&mut rng, dir.path(), "alpha", &alpha_root, 12);
        // Labelled Alpha: an unrelated genome at nearly Alpha's GC, as
        // M. abscessus is beside M. tuberculosis, one at a far GC, and a
        // strain of the other checked species.
        let unrelated = random_genome(&mut rng, 0.64);
        alpha.push(sketch(dir.path(), "foreign.1", &unrelated));
        let far_gc = random_genome(&mut rng, 0.30);
        alpha.push(sketch(dir.path(), "foreign.2", &far_gc));
        let stray = mutate(&mut rng, &gamma_root, 0.005);
        alpha.push(sketch(dir.path(), "foreign.3", &stray));
        let gamma = strains(&mut rng, dir.path(), "gamma", &gamma_root, 8);

        let check = check(
            vec![cohort("Alpha a", alpha), cohort("Gamma g", gamma)],
            &CheckConfig::default(),
        );
        assert_eq!(
            verdicts(&check, Verdict::Fail),
            vec![
                ("foreign.1", "Alpha a"),
                ("foreign.2", "Alpha a"),
                ("foreign.3", "Alpha a"),
            ],
        );
        assert!(verdicts(&check, Verdict::Borderline).is_empty());
        let stray = check.genomes.iter().find(|g| g.genome_id == "foreign.3");
        assert_eq!(stray.unwrap().nearest_species.as_deref(), Some("Gamma g"));
        assert!(check.species[0].reference.starts_with("alpha."));
        let failures = check.failures();
        assert_eq!(failures.len(), 3);
        assert!(
            failures["foreign.1"].contains("fail distance"),
            "{failures:?}"
        );
    }

    /// Shigella sits inside E. coli (about 97-98% ANI), and K. pneumoniae and
    /// K. quasipneumoniae are about 94% apart, just past the species
    /// boundary. Correctly labelled genomes of each pass.
    #[test]
    fn close_relatives_labelled_as_themselves_pass() {
        let dir = tempfile::tempdir().unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(3);
        let coli = random_genome(&mut rng, 0.51);
        let shigella = mutate(&mut rng, &coli, 0.025);
        let pneumoniae = random_genome(&mut rng, 0.57);
        let quasi = mutate(&mut rng, &pneumoniae, 0.06);
        let mut kp = strains(&mut rng, dir.path(), "kp", &pneumoniae, 12);
        // A K. quasipneumoniae filed as K. pneumoniae names its species.
        let filed_wrong = mutate(&mut rng, &quasi, 0.005);
        kp.push(sketch(dir.path(), "kq-as-kp", &filed_wrong));
        // A Shigella filed as E. coli shares its genomic group, so it is
        // recorded, not failed.
        let mut ec = strains(&mut rng, dir.path(), "ec", &coli, 12);
        let shigella_as_coli = mutate(&mut rng, &shigella, 0.005);
        ec.push(sketch(dir.path(), "ss-as-ec", &shigella_as_coli));
        let ss = strains(&mut rng, dir.path(), "ss", &shigella, 10);
        let kq = strains(&mut rng, dir.path(), "kq", &quasi, 10);
        let cohorts = vec![
            cohort("Escherichia coli", ec),
            cohort("Shigella sonnei", ss),
            cohort("Klebsiella pneumoniae", kp),
            cohort("Klebsiella quasipneumoniae", kq),
        ];
        let check = check(cohorts, &CheckConfig::default());
        assert_eq!(
            verdicts(&check, Verdict::Fail),
            vec![("kq-as-kp", "Klebsiella pneumoniae")]
        );
        assert_eq!(
            verdicts(&check, Verdict::Borderline),
            vec![("ss-as-ec", "Escherichia coli")]
        );
        let get = |id: &str| check.genomes.iter().find(|g| g.genome_id == id).unwrap();
        assert_eq!(
            get("kq-as-kp").nearest_species.as_deref(),
            Some("Klebsiella quasipneumoniae")
        );
        let shigella_as_coli = get("ss-as-ec");
        assert!(shigella_as_coli.distance < 0.05);
        assert_eq!(
            shigella_as_coli.nearest_species.as_deref(),
            Some("Shigella sonnei")
        );
    }

    /// Sketches of ten M. tuberculosis genomes from the frozen 2k cohort and
    /// the two genomes revisions 7 and 8 trained on as M. tuberculosis:
    /// 1773.41701 (GC 64.2%, carries blaMAB, likely M. abscessus) and
    /// 1773.40949 (GC 68.4%, unidentified).
    #[test]
    fn the_two_known_genomes_fail_the_tuberculosis_check() {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/species-check-mtb.json.gz");
        let mut text = String::new();
        std::io::Read::read_to_string(&mut fasta::open(&path).unwrap(), &mut text).unwrap();
        let sketches: Vec<GenomeSketch> = serde_json::from_str(&text).unwrap();
        assert_eq!(sketches.len(), 12);
        let check = check(
            vec![cohort("Mycobacterium tuberculosis", sketches)],
            &CheckConfig::default(),
        );
        let failed = verdicts(&check, Verdict::Fail);
        assert_eq!(
            failed,
            vec![
                ("1773.40949", "Mycobacterium tuberculosis"),
                ("1773.41701", "Mycobacterium tuberculosis"),
            ]
        );
        assert!(verdicts(&check, Verdict::Borderline).is_empty());
        for g in check.genomes.iter().filter(|g| g.verdict == Verdict::Fail) {
            assert!(g.distance > 0.25, "{g:?}");
        }
        let gc = |id: &str| check.genomes.iter().find(|g| g.genome_id == id).unwrap().gc;
        assert!((gc("1773.41701") - 0.642).abs() < 0.001);
        assert!((gc("1773.40949") - 0.684).abs() < 0.001);
    }

    #[test]
    fn e_coli_and_every_shigella_share_one_genomic_group() {
        let cfg = CheckConfig::default();
        assert!(cfg.same_group("Escherichia coli", "Shigella sonnei"));
        assert!(cfg.same_group("Shigella flexneri", "Escherichia coli"));
        assert!(cfg.same_group("Shigella flexneri", "Shigella sonnei"));
        assert!(cfg.same_group("Shigella boydii", "Shigella dysenteriae"));
        assert!(!cfg.same_group("Escherichia coli", "Klebsiella pneumoniae"));
        assert!(!cfg.same_group("Escherichia albertii", "Shigella sonnei"));
        assert!(!cfg.same_group("Klebsiella pneumoniae", "Klebsiella quasipneumoniae"));
    }

    /// Within E. coli and Shigella, and between two Shigella, a genome nearer
    /// another member's reference is borderline; without the group it would
    /// fail even inside its own boundary. Its own reference still fails it
    /// past `fail_distance`.
    #[test]
    fn the_genomic_group_exempts_e_coli_and_shigella_from_the_nearer_species_rule() {
        let dir = tempfile::tempdir().unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(5);
        let coli = random_genome(&mut rng, 0.51);
        let sonnei = mutate(&mut rng, &coli, 0.025);
        let flexneri = mutate(&mut rng, &coli, 0.025);
        let mut ec = strains(&mut rng, dir.path(), "ec", &coli, 10);
        let ss_as_ec = mutate(&mut rng, &sonnei, 0.005);
        ec.push(sketch(dir.path(), "ss-as-ec", &ss_as_ec));
        let mut ss = strains(&mut rng, dir.path(), "ss", &sonnei, 10);
        let ec_as_ss = mutate(&mut rng, &coli, 0.005);
        ss.push(sketch(dir.path(), "ec-as-ss", &ec_as_ss));
        let mut sf = strains(&mut rng, dir.path(), "sf", &flexneri, 10);
        let ss_as_sf = mutate(&mut rng, &sonnei, 0.005);
        sf.push(sketch(dir.path(), "ss-as-sf", &ss_as_sf));
        let unrelated = random_genome(&mut rng, 0.51);
        sf.push(sketch(dir.path(), "foreign-as-sf", &unrelated));
        let cohorts = || {
            vec![
                cohort("Escherichia coli", ec.clone()),
                cohort("Shigella sonnei", ss.clone()),
                cohort("Shigella flexneri", sf.clone()),
            ]
        };
        let grouped = check(cohorts(), &CheckConfig::default());
        let crossed = [
            ("ec-as-ss", "Shigella sonnei"),
            ("ss-as-ec", "Escherichia coli"),
            ("ss-as-sf", "Shigella flexneri"),
        ];
        assert_eq!(
            verdicts(&grouped, Verdict::Fail),
            vec![("foreign-as-sf", "Shigella flexneri")]
        );
        let mut borderline = verdicts(&grouped, Verdict::Borderline);
        borderline.sort();
        assert_eq!(borderline, crossed);
        assert!(
            grouped.pairs.iter().all(|p| p.same_group),
            "{:?}",
            grouped.pairs
        );

        let ungrouped = CheckConfig {
            genomic_groups: vec![],
            ..CheckConfig::default()
        };
        let split = check(cohorts(), &ungrouped);
        let mut failed = verdicts(&split, Verdict::Fail);
        failed.sort();
        let mut expected = crossed.to_vec();
        expected.insert(1, ("foreign-as-sf", "Shigella flexneri"));
        assert_eq!(failed, expected);
        let ss_as_ec = split
            .genomes
            .iter()
            .find(|g| g.genome_id == "ss-as-ec")
            .unwrap();
        assert!(ss_as_ec.distance < 0.05, "{ss_as_ec:?}");
        assert_eq!(ss_as_ec.inside_species.as_deref(), Some("Shigella sonnei"));
    }

    #[test]
    fn a_genome_past_the_boundary_with_no_closer_species_is_borderline() {
        let dir = tempfile::tempdir().unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(4);
        let root = random_genome(&mut rng, 0.5);
        let mut genomes = strains(&mut rng, dir.path(), "a", &root, 10);
        let distant = mutate(&mut rng, &root, 0.07);
        genomes.push(sketch(dir.path(), "distant", &distant));
        let check = check(vec![cohort("Alpha a", genomes)], &CheckConfig::default());
        assert!(verdicts(&check, Verdict::Fail).is_empty());
        assert_eq!(
            verdicts(&check, Verdict::Borderline),
            vec![("distant", "Alpha a")]
        );
        assert_eq!(check.species[0].borderline, 1);
    }
}
