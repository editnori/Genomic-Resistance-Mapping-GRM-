//! Species-level presence matrix.
//!
//! Merging every genome's kset gives each distinct k-mer a presence bitset
//! over the genomes. Many k-mers share the same bitset (every k-mer along an
//! accessory gene, every k-mer overlapping one SNP), so k-mers collapse into
//! *patterns*: one row per distinct bitset, plus the list of k-mers that
//! carry it. Learners see patterns; annotation goes back to the k-mers.
//!
//! Files in a matrix directory:
//!   matrix.json        summary and genome order (column i = genomes[i])
//!   patterns.bin       n_patterns rows of `words` little-endian u64
//!   support.bin        u32 per pattern: genomes carrying it
//!   kmer_offsets.bin   u64 per pattern + 1: slice bounds into kmers.bin
//!   kmers.bin          mixed k-mer values grouped by pattern

use crate::bits::{self, BitMatrix, Bits, words_for};
use crate::kset::{BUCKET_BITS, BUCKETS, KsetReader, radix_sort};
use crate::layout::remove_if_present;
use anyhow::{Context, Result, ensure};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MatrixSummary {
    pub k: usize,
    pub genomes: Vec<String>,
    pub min_support: u32,
    pub n_patterns: usize,
    pub kmers_distinct: u64,
    pub kmers_kept: u64,
    /// Present in fewer than `min_support` genomes (singletons included).
    pub kmers_rare: u64,
    /// Present in more than n - min_support genomes (near-core).
    pub kmers_core: u64,
    #[serde(default)]
    pub kset_digests: Vec<String>,
}

pub struct FeatureMatrix {
    pub summary: MatrixSummary,
    pub patterns: BitMatrix,
    pub support: Vec<u32>,
    pub kmer_offsets: Vec<u64>,
    pub kmers: Vec<u64>,
}

impl FeatureMatrix {
    pub fn n_genomes(&self) -> usize {
        self.summary.genomes.len()
    }

    pub fn kmers_of(&self, pattern: usize) -> &[u64] {
        &self.kmers[self.kmer_offsets[pattern] as usize..self.kmer_offsets[pattern + 1] as usize]
    }
}

struct BucketOut {
    rows: Vec<u64>,
    pairs: Vec<(u32, u64)>,
    distinct: u64,
    rare: u64,
    core: u64,
}

/// Merges ksets (one per genome, in column order) into a pattern matrix.
/// A k-mer is kept when min_support <= support <= n - min_support.
pub fn build(
    kset_paths: &[PathBuf],
    min_support: u32,
    threads_chunk: usize,
) -> Result<FeatureMatrix> {
    ensure!(!kset_paths.is_empty(), "no ksets to merge");
    let readers: Vec<KsetReader> = kset_paths
        .iter()
        .map(|p| KsetReader::open(p))
        .collect::<Result<_>>()?;
    let k = readers[0].k;
    ensure!(
        readers.iter().all(|r| r.k == k),
        "ksets mix different k values"
    );
    let n = readers.len();
    ensure!(
        2 * k as u32 - BUCKET_BITS + genome_bits(n) <= 64,
        "{n} genomes at k = {k} do not fit the packed merge; use k <= 31 and at most 4096 genomes"
    );
    let words = words_for(n);
    let max_support = (n as u32).saturating_sub(min_support);

    let mut index: FxHashMap<Box<[u64]>, u32> = FxHashMap::default();
    let mut support: Vec<u32> = Vec::new();
    let mut pairs: Vec<(u32, u64)> = Vec::new();
    let (mut distinct, mut rare, mut core) = (0u64, 0u64, 0u64);

    let chunk = threads_chunk.max(1);
    for start in (0..BUCKETS).step_by(chunk) {
        let end = (start + chunk).min(BUCKETS);
        let outs: Vec<BucketOut> = (start..end)
            .into_par_iter()
            .map(|b| merge_bucket(&readers, b, words, min_support, max_support))
            .collect::<Result<_>>()?;
        for out in outs {
            distinct += out.distinct;
            rare += out.rare;
            core += out.core;
            let mut local_to_global = Vec::with_capacity(out.rows.len() / words.max(1));
            for row in out.rows.chunks_exact(words) {
                let next = support.len() as u32;
                let id = *index.entry(row.into()).or_insert_with(|| {
                    support.push(bits::count(row));
                    next
                });
                local_to_global.push(id);
            }
            pairs.extend(
                out.pairs
                    .into_iter()
                    .map(|(local, kmer)| (local_to_global[local as usize], kmer)),
            );
        }
        if end.is_multiple_of(128) || end == BUCKETS {
            eprintln!(
                "  merged {end} / {BUCKETS} buckets: {distinct} distinct k-mers, {} kept, {} patterns",
                pairs.len(),
                support.len()
            );
        }
    }

    let n_patterns = support.len();
    let mut patterns = BitMatrix {
        rows: n_patterns,
        cols: n,
        words,
        data: vec![0; n_patterns * words],
    };
    for (row, &id) in index.iter() {
        patterns.data[id as usize * words..(id as usize + 1) * words].copy_from_slice(row);
    }
    drop(index);

    pairs.par_sort_unstable();
    let mut kmer_offsets = Vec::with_capacity(n_patterns + 1);
    let mut kmers = Vec::with_capacity(pairs.len());
    let mut cursor = 0;
    for id in 0..n_patterns as u32 {
        kmer_offsets.push(kmers.len() as u64);
        while cursor < pairs.len() && pairs[cursor].0 == id {
            kmers.push(pairs[cursor].1);
            cursor += 1;
        }
    }
    kmer_offsets.push(kmers.len() as u64);

    let summary = MatrixSummary {
        k,
        genomes: readers.iter().map(|r| r.name.clone()).collect(),
        min_support,
        n_patterns,
        kmers_distinct: distinct,
        kmers_kept: kmers.len() as u64,
        kmers_rare: rare,
        kmers_core: core,
        kset_digests: vec![],
    };
    Ok(FeatureMatrix {
        summary,
        patterns,
        support,
        kmer_offsets,
        kmers,
    })
}

/// Bits needed to hold a genome index below `n`.
fn genome_bits(n: usize) -> u32 {
    (usize::BITS - (n.max(2) - 1).leading_zeros()).max(1)
}

/// Merges one bucket across all genomes. Every value in a bucket shares its
/// top BUCKET_BITS bits, so the remaining bits and the genome index pack into
/// one u64; a single radix sort then groups each k-mer's genomes together.
fn merge_bucket(
    readers: &[KsetReader],
    bucket: usize,
    words: usize,
    min: u32,
    max: u32,
) -> Result<BucketOut> {
    let low_bits = 2 * readers[0].k as u32 - BUCKET_BITS;
    let gbits = genome_bits(readers.len());
    let low_mask = (1u64 << low_bits) - 1;
    let total: usize = readers.iter().map(|r| r.bucket_len(bucket)).sum();
    let mut packed: Vec<u64> = Vec::with_capacity(total);
    let (mut bytes, mut values) = (Vec::new(), Vec::new());
    for (g, r) in readers.iter().enumerate() {
        values.clear();
        r.read_bucket(bucket, &mut bytes, &mut values)?;
        packed.extend(values.iter().map(|&v| ((v & low_mask) << gbits) | g as u64));
    }
    radix_sort(&mut packed, &mut values);
    let base = (bucket as u64) << low_bits;
    let gmask = (1u64 << gbits) - 1;

    let mut local: FxHashMap<Box<[u64]>, u32> = FxHashMap::default();
    let mut out = BucketOut {
        rows: Vec::new(),
        pairs: Vec::new(),
        distinct: 0,
        rare: 0,
        core: 0,
    };
    let mut row = vec![0u64; words];

    let mut i = 0;
    while i < packed.len() {
        let key = packed[i] >> gbits;
        let value = base | key;
        row.iter_mut().for_each(|w| *w = 0);
        let mut support = 0u32;
        while i < packed.len() && packed[i] >> gbits == key {
            bits::set(&mut row, (packed[i] & gmask) as usize);
            support += 1;
            i += 1;
        }
        out.distinct += 1;
        if support < min {
            out.rare += 1;
            continue;
        }
        if support > max {
            out.core += 1;
            continue;
        }
        let next = local.len() as u32;
        let id = *local.entry(row.as_slice().into()).or_insert_with(|| {
            out.rows.extend_from_slice(&row);
            next
        });
        out.pairs.push((id, value));
    }
    Ok(out)
}

fn write_u64s(path: &Path, values: &[u64]) -> Result<()> {
    let mut w =
        BufWriter::new(File::create(path).with_context(|| format!("create {}", path.display()))?);
    for v in values {
        w.write_all(&v.to_le_bytes())?;
    }
    w.flush()?;
    Ok(())
}

fn read_u64s(path: &Path) -> Result<Vec<u64>> {
    let mut bytes = Vec::new();
    BufReader::new(File::open(path).with_context(|| format!("open {}", path.display()))?)
        .read_to_end(&mut bytes)?;
    Ok(bytes
        .as_chunks::<8>()
        .0
        .iter()
        .map(|c| u64::from_le_bytes(*c))
        .collect())
}

impl FeatureMatrix {
    pub fn save(&self, dir: &Path) -> Result<()> {
        std::fs::create_dir_all(dir)?;
        // A summary left from an earlier save would mark half-written
        // binaries as a complete matrix, so it goes before they are touched.
        let summary = dir.join("matrix.json");
        remove_if_present(&summary).with_context(|| format!("remove {}", summary.display()))?;
        write_u64s(&dir.join("patterns.bin"), &self.patterns.data)?;
        let support: Vec<u8> = self.support.iter().flat_map(|s| s.to_le_bytes()).collect();
        std::fs::write(dir.join("support.bin"), support)?;
        write_u64s(&dir.join("kmer_offsets.bin"), &self.kmer_offsets)?;
        write_u64s(&dir.join("kmers.bin"), &self.kmers)?;
        // Summary last: its presence marks a complete matrix.
        std::fs::write(&summary, serde_json::to_vec_pretty(&self.summary)?)?;
        Ok(())
    }

    pub fn load(dir: &Path) -> Result<Self> {
        let summary: MatrixSummary = serde_json::from_slice(
            &std::fs::read(dir.join("matrix.json"))
                .with_context(|| format!("no matrix in {}", dir.display()))?,
        )?;
        let n = summary.genomes.len();
        let words = words_for(n);
        let data = read_u64s(&dir.join("patterns.bin"))?;
        ensure!(
            data.len() == summary.n_patterns * words,
            "patterns.bin size does not match matrix.json"
        );
        let support: Vec<u32> = std::fs::read(dir.join("support.bin"))?
            .as_chunks::<4>()
            .0
            .iter()
            .map(|c| u32::from_le_bytes(*c))
            .collect();
        ensure!(
            support.len() == summary.n_patterns,
            "support.bin size does not match matrix.json"
        );
        let kmer_offsets = read_u64s(&dir.join("kmer_offsets.bin"))?;
        let kmers = read_u64s(&dir.join("kmers.bin"))?;
        ensure!(
            kmer_offsets.len() == summary.n_patterns + 1,
            "kmer_offsets.bin is incomplete"
        );
        ensure!(
            kmer_offsets[summary.n_patterns] == kmers.len() as u64,
            "kmers.bin size does not match kmer_offsets.bin"
        );
        Ok(FeatureMatrix {
            patterns: BitMatrix {
                rows: summary.n_patterns,
                cols: n,
                words,
                data,
            },
            summary,
            support,
            kmer_offsets,
            kmers,
        })
    }
}

/// A learning problem: features restricted to the labeled genomes, with
/// patterns that became identical (or constant) on that subset collapsed.
/// This uses every labeled genome, so it is only a candidate superset: each
/// fit applies `learn::training_candidates` to its own training samples.
pub struct Dataset {
    /// Features x samples.
    pub x: BitMatrix,
    /// Bit i set = sample i is resistant (the positive class).
    pub y: Bits,
    /// Column of each sample in the species matrix.
    pub columns: Vec<usize>,
    /// Feature f covers species patterns `feature_patterns[offsets[f]..offsets[f+1]]`.
    pub feature_offsets: Vec<u32>,
    pub feature_patterns: Vec<u32>,
}

impl Dataset {
    pub fn n_samples(&self) -> usize {
        self.x.cols
    }

    pub fn n_features(&self) -> usize {
        self.x.rows
    }

    pub fn patterns_of(&self, feature: usize) -> &[u32] {
        &self.feature_patterns
            [self.feature_offsets[feature] as usize..self.feature_offsets[feature + 1] as usize]
    }

    /// Projects `matrix` onto `columns` (with labels `resistant`), keeping
    /// features present in at least `min_support` and at most
    /// n - `min_support` of those samples.
    pub fn project(
        matrix: &FeatureMatrix,
        columns: Vec<usize>,
        resistant: &[bool],
        min_support: u32,
    ) -> Dataset {
        assert_eq!(columns.len(), resistant.len());
        let n = columns.len();
        let words = words_for(n);
        let max = (n as u32).saturating_sub(min_support);
        let rows = matrix.patterns.rows;

        // Gather in parallel chunks; collapse duplicates within each chunk
        // first so the sequential merge sees far fewer rows.
        let chunk = 1 << 15;
        let parts: Vec<(Vec<u64>, Vec<Vec<u32>>)> = (0..rows.div_ceil(chunk))
            .into_par_iter()
            .map(|c| {
                let mut local: FxHashMap<Box<[u64]>, usize> = FxHashMap::default();
                let mut out_rows = Vec::new();
                let mut members: Vec<Vec<u32>> = Vec::new();
                let mut row = vec![0u64; words];
                for p in c * chunk..((c + 1) * chunk).min(rows) {
                    let src = matrix.patterns.row(p);
                    row.iter_mut().for_each(|w| *w = 0);
                    for (i, &col) in columns.iter().enumerate() {
                        if bits::get(src, col) {
                            bits::set(&mut row, i);
                        }
                    }
                    let s = bits::count(&row);
                    if s < min_support || s > max {
                        continue;
                    }
                    let next = members.len();
                    let id = *local.entry(row.as_slice().into()).or_insert_with(|| {
                        out_rows.extend_from_slice(&row);
                        members.push(Vec::new());
                        next
                    });
                    members[id].push(p as u32);
                }
                (out_rows, members)
            })
            .collect();

        let mut index: FxHashMap<Box<[u64]>, usize> = FxHashMap::default();
        let mut x = BitMatrix::with_cols(n);
        let mut members: Vec<Vec<u32>> = Vec::new();
        for (part_rows, part_members) in parts {
            for (row, m) in part_rows.chunks_exact(words).zip(part_members) {
                let next = members.len();
                let id = *index.entry(row.into()).or_insert_with(|| {
                    x.push_row(row);
                    members.push(Vec::new());
                    next
                });
                members[id].extend(m);
            }
        }
        let mut feature_offsets = Vec::with_capacity(members.len() + 1);
        let mut feature_patterns = Vec::new();
        for m in &members {
            feature_offsets.push(feature_patterns.len() as u32);
            feature_patterns.extend_from_slice(m);
        }
        feature_offsets.push(feature_patterns.len() as u32);

        let y = Bits::from_indices(
            n,
            resistant
                .iter()
                .enumerate()
                .filter(|(_, r)| **r)
                .map(|(i, _)| i),
        );
        Dataset {
            x,
            y,
            columns,
            feature_offsets,
            feature_patterns,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::KmerSpace;
    use crate::kset;

    fn kset_of(dir: &Path, name: &str, seq: &str, space: &KmerSpace) -> PathBuf {
        let fna = dir.join(format!("{name}.fna"));
        std::fs::write(&fna, format!(">c\n{seq}\n")).unwrap();
        let (values, _) = kset::count_file(&fna, space, 1).unwrap();
        let path = dir.join(format!("{name}.kset"));
        kset::write(&path, name, space, &values).unwrap();
        path
    }

    #[test]
    fn merge_groups_kmers_by_presence_pattern() {
        let dir = tempfile::tempdir().unwrap();
        let space = KmerSpace::new(7);
        // A core block shared by all, plus an accessory block in g1 and g2.
        let core = "ACGTTGCAAGGCTTACCGATA";
        let accessory = "TTTGGGCCCAAATGTGTCACA";
        let paths = vec![
            kset_of(dir.path(), "g0", core, &space),
            kset_of(dir.path(), "g1", &format!("{core}NNNN{accessory}"), &space),
            kset_of(dir.path(), "g2", &format!("{accessory}NNNN{core}"), &space),
        ];
        let m = build(&paths, 1, 8).unwrap();
        assert_eq!(m.summary.genomes, vec!["g0", "g1", "g2"]);
        // Core k-mers (support 3 > n - 1) are dropped; accessory ones form a
        // single pattern {g1, g2}.
        assert_eq!(m.n_genomes(), 3);
        assert_eq!(m.summary.n_patterns, 1);
        assert_eq!(m.support, vec![2]);
        assert!(!m.patterns.get(0, 0) && m.patterns.get(0, 1) && m.patterns.get(0, 2));
        let mut distinct = std::collections::HashSet::new();
        space.for_each_canonical(accessory.as_bytes(), |x| {
            distinct.insert(x);
        });
        assert_eq!(m.kmers_of(0).len(), distinct.len());
        assert!(m.summary.kmers_core > 0);
        for &v in m.kmers_of(0) {
            let fwd = space.decode(space.unmix(v));
            let rc = space.decode(space.reverse_complement(space.unmix(v)));
            assert!(accessory.contains(&fwd) || accessory.contains(&rc));
        }

        m.save(&dir.path().join("mx")).unwrap();
        let loaded = FeatureMatrix::load(&dir.path().join("mx")).unwrap();
        assert_eq!(loaded.patterns.data, m.patterns.data);
        assert_eq!(loaded.kmers, m.kmers);

        // On columns {g0, g1} the accessory pattern is still variable.
        let ds = Dataset::project(&loaded, vec![0, 1], &[false, true], 1);
        assert_eq!(ds.n_features(), 1);
        assert_eq!(ds.patterns_of(0), &[0]);
        // On {g1, g2} it is constant and disappears.
        let ds = Dataset::project(&loaded, vec![1, 2], &[false, true], 1);
        assert_eq!(ds.n_features(), 0);
    }

    #[test]
    fn load_rejects_truncated_files() {
        let dir = tempfile::tempdir().unwrap();
        let space = KmerSpace::new(7);
        let core = "ACGTTGCAAGGCTTACCGATA";
        let accessory = "TTTGGGCCCAAATGTGTCACA";
        let paths = vec![
            kset_of(dir.path(), "g0", core, &space),
            kset_of(dir.path(), "g1", &format!("{core}NNNN{accessory}"), &space),
        ];
        let m = build(&paths, 1, 8).unwrap();
        assert_eq!(m.summary.n_patterns, 1);
        // An interrupted save or copy leaves a short file behind.
        for (file, cut) in [("support.bin", 4), ("kmers.bin", 8)] {
            let mx = dir.path().join(file);
            m.save(&mx).unwrap();
            FeatureMatrix::load(&mx).unwrap();
            let path = mx.join(file);
            let len = std::fs::metadata(&path).unwrap().len();
            let f = std::fs::OpenOptions::new().write(true).open(&path).unwrap();
            f.set_len(len - cut).unwrap();
            assert!(
                FeatureMatrix::load(&mx).is_err(),
                "{file} short by {cut} bytes"
            );
        }
    }

    #[test]
    fn a_failed_resave_leaves_no_summary() {
        let dir = tempfile::tempdir().unwrap();
        let space = KmerSpace::new(7);
        let core = "ACGTTGCAAGGCTTACCGATA";
        let accessory = "TTTGGGCCCAAATGTGTCACA";
        let paths = vec![
            kset_of(dir.path(), "g0", core, &space),
            kset_of(dir.path(), "g1", &format!("{core}NNNN{accessory}"), &space),
        ];
        let m = build(&paths, 1, 8).unwrap();
        let mx = dir.path().join("mx");
        m.save(&mx).unwrap();
        // The second save rewrites three binaries, then cannot create kmers.bin.
        std::fs::remove_file(mx.join("kmers.bin")).unwrap();
        std::fs::create_dir(mx.join("kmers.bin")).unwrap();
        assert!(m.save(&mx).is_err());
        assert!(!mx.join("matrix.json").exists());
        let err = FeatureMatrix::load(&mx).err().unwrap();
        assert!(format!("{err:#}").contains("no matrix in"), "{err:#}");
    }
}
