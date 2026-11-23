//! Experiment, off by default (plan block C3, docs/input-arms.md): the
//! label-free builds of the three matrix arms of `grm learn --exp-inputs`.
//!
//! `grm exp-input-build KIND --species S --out DIR` builds, once per species
//! over every genome of the species matrix and in its column order, one
//! pattern matrix in the release's own format (`DIR/matrix`), which
//! `grm learn --exp-inputs ARM --exp-input-matrix DIR` projects onto each
//! cohort. Each genome's values go into a kset (`kset::write`) and
//! `matrix::build` merges them at `MATRIX_SUPPORT`, so identical patterns
//! collapse as they do for the species matrix:
//! - `counts`: the canonical 31-mers a genome's assembly holds 2 or more
//!   times, by the release's own counter (`kset::count_file` with its
//!   abundance floor at 2);
//! - `aa`: the amino-acid 10-mers of the proteins Prodigal called, each
//!   protein first held to GRM's own table 11 translation of the contig at
//!   Prodigal's coordinates; stop and ambiguous residues end a window;
//! - `unitigs`: the unitigs of the compacted de Bruijn graph over every
//!   31-mer of the species' genomes (single-genome 31-mers included), built
//!   by GRM or read from an outside build's FASTA; a genome holds a unitig
//!   when its kset holds every 31-mer of it.
//!
//! Nothing here reads a label, a phenotype, a fold or a seed. The build's
//! `build.json` holds counts only; its ksets, matrix and `genomes.tsv` name
//! genomes and stay with the run.

use crate::dna::KmerSpace;
use crate::kset::{self, BUCKETS, KsetReader};
use crate::matrix::{self, FeatureMatrix, MatrixSummary};
use anyhow::{Context, Result, bail, ensure};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

/// The build file's format. A change to what it holds changes it.
pub const BUILD_FORMAT: &str = "grm-c3-build-1";
/// The support every added matrix is merged at: the species matrix's.
pub const MATRIX_SUPPORT: u32 = 2;
/// The nucleotide k of the counts and unitig arms.
pub const K: usize = 31;
/// A genome's 31-mer enters the counts arm at this many copies.
pub const MIN_COUNT: u32 = 2;
/// The amino-acid k (the lead's ruling on question 5).
pub const AA_K: usize = 10;
/// The 20 standard residues, each a 5-bit code by its place here.
pub const RESIDUES: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
/// An amino-acid 10-mer is 50 bits; a kset of this k mixes 50 bits.
pub const AA_SPACE_K: usize = 25;
/// A unitig number is 32 bits; a kset of this k mixes 32 bits.
pub const UNITIG_SPACE_K: usize = 16;

/// Which label-free build.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Counts,
    Aa,
    Unitigs,
}

impl Kind {
    pub fn name(self) -> &'static str {
        match self {
            Kind::Counts => "counts",
            Kind::Aa => "aa",
            Kind::Unitigs => "unitigs",
        }
    }

    /// The k each genome's kset is written at.
    pub fn space_k(self) -> usize {
        match self {
            Kind::Counts => K,
            Kind::Aa => AA_SPACE_K,
            Kind::Unitigs => UNITIG_SPACE_K,
        }
    }
}

/// Which graph the unitig arm reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Graph {
    /// GRM's own compaction of every 31-mer of the species' ksets.
    Grm,
    /// An outside build's unitigs FASTA (bcalm2 or Bifrost).
    Outside,
}

// ---------------------------------------------------------------------------
// Counts

/// A genome's canonical 31-mers seen `MIN_COUNT` or more times in its
/// assembly, mixed and sorted: the release's own count with its abundance
/// floor raised.
pub fn count_values(path: &Path, space: &KmerSpace) -> Result<Vec<u64>> {
    Ok(kset::count_file(path, space, MIN_COUNT)?.0)
}

// ---------------------------------------------------------------------------
// Amino acids

/// The 5-bit code of a standard residue; None for a stop, X or anything
/// else.
pub fn residue_code(b: u8) -> Option<u64> {
    RESIDUES
        .iter()
        .position(|&r| r == b.to_ascii_uppercase())
        .map(|i| i as u64)
}

/// Calls `f` with the 50-bit code of every amino-acid 10-mer of `protein`
/// whose 10 residues are all standard: a stop or an ambiguous residue ends
/// the window, so no 10-mer spans one.
pub fn for_each_aa_kmer(protein: &[u8], mut f: impl FnMut(u64)) {
    let mask = (1u64 << (5 * AA_K)) - 1;
    let (mut v, mut filled) = (0u64, 0usize);
    for &b in protein {
        let Some(c) = residue_code(b) else {
            filled = 0;
            continue;
        };
        v = ((v << 5) | c) & mask;
        filled += 1;
        if filled >= AA_K {
            f(v);
        }
    }
}

/// One protein of Prodigal's `-a` output, with its header's coordinates.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Called {
    pub contig: String,
    /// 1-based, inclusive, on the contig's forward strand.
    pub start: usize,
    pub end: usize,
    pub forward: bool,
    /// Prodigal's `partial=LR`: the left and right ends in contig order run
    /// off the sequence (true) or are a real start or stop (false).
    pub partial: (bool, bool),
    pub protein: Vec<u8>,
}

/// Parses a Prodigal header: `<contig>_<n> # <start> # <end> # <strand> #
/// ID=...;partial=LR;...`.
pub fn parse_header(header: &str) -> Result<Called> {
    let parts: Vec<&str> = header.split(" # ").collect();
    ensure!(
        parts.len() == 5,
        "a protein header with {} fields, not Prodigal's 5",
        parts.len()
    );
    let id = parts[0].split_whitespace().next().unwrap_or("");
    let (contig, n) = id
        .rsplit_once('_')
        .context("a protein id with no gene number")?;
    ensure!(
        !contig.is_empty() && n.parse::<u32>().is_ok(),
        "a protein id that is not <contig>_<n>"
    );
    let start: usize = parts[1].trim().parse().context("a protein start")?;
    let end: usize = parts[2].trim().parse().context("a protein end")?;
    let forward = match parts[3].trim() {
        "1" => true,
        "-1" => false,
        _ => bail!("a protein strand that is not 1 or -1"),
    };
    let partial = parts[4]
        .split(';')
        .find_map(|kv| kv.strip_prefix("partial="))
        .context("a protein header with no partial=")?
        .as_bytes();
    ensure!(
        partial.len() == 2 && partial.iter().all(|c| *c == b'0' || *c == b'1'),
        "a partial= that is not two of 0 and 1"
    );
    ensure!(
        start >= 1 && end >= start,
        "a protein's coordinates are out of order"
    );
    Ok(Called {
        contig: contig.to_string(),
        start,
        end,
        forward,
        partial: (partial[0] == b'1', partial[1] == b'1'),
        protein: Vec::new(),
    })
}

/// Reads Prodigal's protein FASTA: each record's full header and its
/// residues.
pub fn read_proteins(path: &Path) -> Result<Vec<Called>> {
    let file = std::fs::File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut out: Vec<Called> = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line?;
        let line = line.trim_end();
        if let Some(h) = line.strip_prefix('>') {
            out.push(
                parse_header(h)
                    .with_context(|| format!("protein {} of {}", out.len() + 1, path.display()))?,
            );
        } else if !line.is_empty() {
            out.last_mut()
                .context("residues before the first protein header")?
                .protein
                .extend_from_slice(line.as_bytes());
        }
    }
    Ok(out)
}

/// The standard genetic code (NCBI tables 1 and 11 read every codon alike;
/// they differ only in which codons may start a gene). Bases in TCAG order.
const CODE11: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

fn tcag(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'T' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// GRM's own translation at table 11 of a coding sequence read 5' to 3';
/// a codon with any other base than ACGT is X.
pub fn translate11(cds: &[u8]) -> Vec<u8> {
    cds.as_chunks::<3>()
        .0
        .iter()
        .map(|c| match (tcag(c[0]), tcag(c[1]), tcag(c[2])) {
            (Some(a), Some(b), Some(d)) => CODE11[a * 16 + b * 4 + d],
            _ => b'X',
        })
        .collect()
}

/// Whether Prodigal's protein is GRM's own table 11 translation of the
/// contig at its coordinates: the first residue M where the gene has a real
/// start, as Prodigal writes any start codon, and the stop kept.
pub fn protein_matches(contig: &[u8], c: &Called) -> bool {
    if c.end > contig.len() || !(c.end - c.start + 1).is_multiple_of(3) {
        return false;
    }
    let seg = &contig[c.start - 1..c.end];
    let cds: Vec<u8> = if c.forward {
        seg.to_vec()
    } else {
        seg.iter().rev().map(|&b| complement(b)).collect()
    };
    let mut want = translate11(&cds);
    // The 5' end is the left one on the forward strand, the right one on the
    // reverse strand.
    let start_real = if c.forward {
        !c.partial.0
    } else {
        !c.partial.1
    };
    if start_real && !want.is_empty() {
        want[0] = b'M';
    }
    want == c.protein
}

/// One genome's amino-acid 10-mers, mixed into the 50-bit space and sorted,
/// after every protein is held to the genome's contigs. A protein that GRM's
/// translation does not give stops the build (question 20 of the C3 arms
/// entry in docs/assumption-plan-corrections.md).
pub fn aa_values(genome: &Path, proteins: &Path, space: &KmerSpace) -> Result<(Vec<u64>, usize)> {
    let mut contigs: HashMap<String, Vec<u8>> = HashMap::new();
    crate::fasta::for_each_record(genome, |name, seq| {
        contigs.insert(name.to_string(), seq.to_vec());
    })?;
    let called = read_proteins(proteins)?;
    let bad = called
        .iter()
        .filter(|c| {
            contigs
                .get(&c.contig)
                .is_none_or(|s| !protein_matches(s, c))
        })
        .count();
    ensure!(
        bad == 0,
        "{bad} of {} proteins in {} are not GRM's table 11 translation of the genome at their coordinates",
        called.len(),
        proteins.display()
    );
    let mut values = Vec::new();
    for c in &called {
        for_each_aa_kmer(&c.protein, |v| values.push(space.mix(v)));
    }
    values.sort_unstable();
    values.dedup();
    Ok((values, called.len()))
}

// ---------------------------------------------------------------------------
// Unitigs

/// A sorted set of mixed 31-mers with a prefix index, and each one's unitig.
pub struct Table {
    pub values: Vec<u64>,
    /// `values[index[p]..index[p + 1]]` hold the values whose top `prefix`
    /// bits are p.
    index: Vec<u64>,
    shift: u32,
    /// Each value's unitig number, filled by the compaction or the reader.
    pub unitig: Vec<u32>,
    /// Each unitig's 31-mer count.
    pub sizes: Vec<u32>,
}

impl Table {
    /// Indexes sorted, distinct values of a `bits`-bit space.
    pub fn new(values: Vec<u64>, bits: u32) -> Table {
        let prefix = (64 - (values.len() as u64 / 8).max(1).leading_zeros()).clamp(1, bits);
        let shift = bits - prefix;
        let mut index = vec![0u64; (1usize << prefix) + 1];
        for &v in &values {
            index[(v >> shift) as usize + 1] += 1;
        }
        for p in 1..index.len() {
            index[p] += index[p - 1];
        }
        Table {
            values,
            index,
            shift,
            unitig: Vec::new(),
            sizes: Vec::new(),
        }
    }

    pub fn find(&self, v: u64) -> Option<usize> {
        let p = (v >> self.shift) as usize;
        let (lo, hi) = (self.index[p] as usize, self.index[p + 1] as usize);
        self.values[lo..hi].binary_search(&v).ok().map(|i| lo + i)
    }

    fn set_unitigs(&mut self, unitig: Vec<u32>, count: usize) {
        let mut sizes = vec![0u32; count];
        for &u in &unitig {
            sizes[u as usize] += 1;
        }
        self.unitig = unitig;
        self.sizes = sizes;
    }
}

/// Every distinct value of the ksets: the union over the species' genomes,
/// sorted. A value held by one genome is kept. Each bucket is read over
/// every genome, then only its distinct values are kept, so the memory is
/// the distinct values twice plus one bucket's reads a thread, not the sum
/// of every genome's values.
pub fn union_values(readers: &[KsetReader]) -> Result<Vec<u64>> {
    let parts: Vec<Vec<u64>> = (0..BUCKETS)
        .into_par_iter()
        .map(|b| {
            let mut bytes = Vec::new();
            let mut out = Vec::new();
            for r in readers {
                r.read_bucket(b, &mut bytes, &mut out)?;
            }
            out.sort_unstable();
            out.dedup();
            Ok(out.as_slice().to_vec())
        })
        .collect::<Result<_>>()?;
    let mut all = Vec::with_capacity(parts.iter().map(Vec::len).sum());
    for p in parts {
        all.extend_from_slice(&p);
    }
    Ok(all)
}

const NONE: u32 = u32::MAX;

/// The canonical form of an oriented k-mer, and whether it is the forward one.
fn canon(space: &KmerSpace, y: u64) -> (u64, bool) {
    let r = space.reverse_complement(y);
    if y <= r { (y, true) } else { (r, false) }
}

/// The oriented k-mers following `x` (its successors) that the table holds,
/// as (index, oriented value), at most 4.
fn successors(t: &Table, space: &KmerSpace, x: u64, out: &mut Vec<(usize, u64)>) {
    out.clear();
    let mask = (1u64 << space.bits()) - 1;
    for b in 0..4u64 {
        let y = ((x << 2) | b) & mask;
        if let Some(i) = t.find(space.mix(canon(space, y).0)) {
            out.push((i, y));
        }
    }
}

fn predecessors(t: &Table, space: &KmerSpace, x: u64, out: &mut Vec<(usize, u64)>) {
    out.clear();
    let shift = space.bits() - 2;
    for b in 0..4u64 {
        let y = (x >> 2) | (b << shift);
        if let Some(i) = t.find(space.mix(canon(space, y).0)) {
            out.push((i, y));
        }
    }
}

/// The link on `i`'s right (`right`) or left side: the one neighbour when
/// `i` has one there, that neighbour has `i` alone on its facing side, and
/// it is another 31-mer.
fn link(
    t: &Table,
    space: &KmerSpace,
    i: usize,
    right: bool,
    a: &mut Vec<(usize, u64)>,
    b: &mut Vec<(usize, u64)>,
) -> u32 {
    let x = space.unmix(t.values[i]);
    if right {
        successors(t, space, x, a);
    } else {
        predecessors(t, space, x, a);
    }
    if a.len() != 1 || a[0].0 == i {
        return NONE;
    }
    let (j, y) = a[0];
    if right {
        predecessors(t, space, y, b);
    } else {
        successors(t, space, y, b);
    }
    if b.len() == 1 { j as u32 } else { NONE }
}

fn find_root(parent: &mut [u32], mut i: u32) -> u32 {
    while parent[i as usize] != i {
        let p = parent[i as usize];
        parent[i as usize] = parent[p as usize];
        i = p;
    }
    i
}

/// GRM's own compaction: each 31-mer's unitig number in `t`, unitigs
/// numbered by their smallest index in the table. The unitigs are the
/// connected parts of the compacted edges, so a path or a cycle is one.
/// Returns the unitig count.
pub fn compact(t: &mut Table, space: &KmerSpace) -> Result<usize> {
    let n = t.values.len();
    ensure!(
        n < NONE as usize,
        "{n} 31-mers exceed the unitig table's u32 index"
    );
    let links: Vec<(u32, u32)> = (0..n)
        .into_par_iter()
        .map_init(
            || (Vec::with_capacity(4), Vec::with_capacity(4)),
            |(a, b), i| {
                (
                    link(t, space, i, true, a, b),
                    link(t, space, i, false, a, b),
                )
            },
        )
        .collect();
    let mut parent: Vec<u32> = (0..n as u32).collect();
    for (i, &(r, l)) in links.iter().enumerate() {
        for j in [r, l] {
            if j == NONE {
                continue;
            }
            let (a, b) = (find_root(&mut parent, i as u32), find_root(&mut parent, j));
            if a != b {
                let (lo, hi) = if a < b { (a, b) } else { (b, a) };
                parent[hi as usize] = lo;
            }
        }
    }
    drop(links);
    // A root is its part's smallest index, so numbering roots in index
    // order numbers the unitigs by their smallest 31-mer.
    let mut count = 0u32;
    let mut unitig = vec![0u32; n];
    for i in 0..n {
        let r = find_root(&mut parent, i as u32) as usize;
        if r == i {
            unitig[i] = count;
            count += 1;
        } else {
            unitig[i] = unitig[r];
        }
    }
    drop(parent);
    t.set_unitigs(unitig, count as usize);
    Ok(count as usize)
}

/// Reads an outside build's unitigs FASTA: each record is one unitig, and
/// each of its canonical 31-mers is that unitig's. A 31-mer in two unitigs
/// is refused.
pub fn read_unitigs(path: &Path, space: &KmerSpace) -> Result<(Table, usize)> {
    let mut pairs: Vec<(u64, u32)> = Vec::new();
    let mut count = 0u32;
    crate::fasta::for_each_record(path, |_, seq| {
        space.for_each_canonical(seq, |x| pairs.push((space.mix(x), count)));
        count += 1;
    })?;
    pairs.par_sort_unstable();
    pairs.dedup();
    let twice = pairs.windows(2).filter(|w| w[0].0 == w[1].0).count();
    ensure!(
        twice == 0,
        "{twice} 31-mers of {} are in two unitigs",
        path.display()
    );
    let (values, unitig): (Vec<u64>, Vec<u32>) = pairs.into_iter().unzip();
    let mut t = Table::new(values, space.bits());
    t.set_unitigs(unitig, count as usize);
    Ok((t, count as usize))
}

/// The equality check of question 1: each unitig is the set of its 31-mers.
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Equality {
    pub kmers_grm: u64,
    pub kmers_outside: u64,
    pub kmers_both: u64,
    pub kmers_grm_only: u64,
    pub kmers_outside_only: u64,
    pub unitigs_grm: u64,
    pub unitigs_outside: u64,
    /// Unitigs whose 31-mer sets are equal in both builds.
    pub unitigs_both: u64,
    pub unitigs_grm_only: u64,
    pub unitigs_outside_only: u64,
}

pub fn compare(grm: &Table, outside: &Table) -> Equality {
    const CONFLICT: u32 = u32::MAX - 1;
    let mut g_map = vec![NONE; grm.sizes.len()];
    let mut o_map = vec![NONE; outside.sizes.len()];
    let mark = |m: &mut [u32], at: u32, to: u32| {
        let cell = &mut m[at as usize];
        *cell = if *cell == NONE || *cell == to {
            to
        } else {
            CONFLICT
        };
    };
    let (mut i, mut j) = (0usize, 0usize);
    let (mut both, mut g_only, mut o_only) = (0u64, 0u64, 0u64);
    while i < grm.values.len() || j < outside.values.len() {
        let gv = grm.values.get(i).copied().unwrap_or(u64::MAX);
        let ov = outside.values.get(j).copied().unwrap_or(u64::MAX);
        if i < grm.values.len() && (j >= outside.values.len() || gv < ov) {
            g_map[grm.unitig[i] as usize] = CONFLICT;
            g_only += 1;
            i += 1;
        } else if j < outside.values.len() && (i >= grm.values.len() || ov < gv) {
            o_map[outside.unitig[j] as usize] = CONFLICT;
            o_only += 1;
            j += 1;
        } else {
            let (g, o) = (grm.unitig[i], outside.unitig[j]);
            mark(&mut g_map, g, o);
            mark(&mut o_map, o, g);
            both += 1;
            i += 1;
            j += 1;
        }
    }
    let equal = g_map
        .iter()
        .enumerate()
        .filter(|&(g, &o)| o < CONFLICT && o_map[o as usize] == g as u32)
        .count() as u64;
    Equality {
        kmers_grm: grm.values.len() as u64,
        kmers_outside: outside.values.len() as u64,
        kmers_both: both,
        kmers_grm_only: g_only,
        kmers_outside_only: o_only,
        unitigs_grm: grm.sizes.len() as u64,
        unitigs_outside: outside.sizes.len() as u64,
        unitigs_both: equal,
        unitigs_grm_only: grm.sizes.len() as u64 - equal,
        unitigs_outside_only: outside.sizes.len() as u64 - equal,
    }
}

/// One genome against the graph: the unitigs it holds whole, sorted, and
/// how many it holds in part and how many of its 31-mers no unitig has.
pub fn presence(t: &Table, kmers: &[u64]) -> (Vec<u32>, u64, u64) {
    let mut ids: Vec<u32> = Vec::with_capacity(kmers.len());
    let mut missing = 0u64;
    for &v in kmers {
        match t.find(v) {
            Some(i) => ids.push(t.unitig[i]),
            None => missing += 1,
        }
    }
    ids.sort_unstable();
    let mut held = Vec::new();
    let mut part = 0u64;
    let mut s = 0;
    while s < ids.len() {
        let mut e = s + 1;
        while e < ids.len() && ids[e] == ids[s] {
            e += 1;
        }
        if (e - s) as u32 == t.sizes[ids[s] as usize] {
            held.push(ids[s]);
        } else {
            part += 1;
        }
        s = e;
    }
    (held, part, missing)
}

// ---------------------------------------------------------------------------
// The build

/// What one build made: counts only.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Build {
    pub format: String,
    pub evaluation_revision: u32,
    pub kind: String,
    pub species: String,
    pub genomes: usize,
    /// The k of the added ksets: 31 for counts, 25 (a 50-bit space of
    /// amino-acid 10-mers) for aa, 16 (a 32-bit space of unitig numbers).
    pub space_k: usize,
    pub matrix_support: u32,
    /// Values over every genome's kset, summed.
    pub values_summed: u64,
    pub matrix_patterns: usize,
    pub matrix_values_kept: u64,
    pub matrix_values_rare: u64,
    pub matrix_values_core: u64,
    /// aa: proteins read, each equal to GRM's translation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub proteins: Option<u64>,
    /// unitigs: which graph, its 31-mers and unitigs, and the genome and
    /// unitig pairs held whole and in part.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub graph: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub graph_kmers: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unitigs: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pairs_held: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pairs_in_part: Option<u64>,
    /// Genome 31-mers no unitig of an outside graph holds.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kmers_without_unitig: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub equality: Option<Equality>,
}

/// Where a build's per-genome counts go: genome id, values, and for the
/// unitig arm the unitigs held whole and in part. It names genomes and stays
/// with the run.
pub const GENOMES_FILE: &str = "genomes.tsv";

/// The arguments of one build.
pub struct Request<'a> {
    pub kind: Kind,
    pub species: &'a str,
    /// The species matrix's genomes, in its column order.
    pub genomes: &'a [String],
    pub fasta_dir: &'a Path,
    pub kset_dir: &'a Path,
    /// aa: Prodigal's `<genome id>.faa`.
    pub proteins: Option<&'a Path>,
    pub graph: Graph,
    /// unitigs: the outside build's FASTA, read as the graph (`Graph::Outside`)
    /// or compared with GRM's (`Graph::Grm`).
    pub outside: Option<&'a Path>,
    pub out: &'a Path,
}

fn write_json<T: Serialize>(path: &Path, value: &T) -> Result<()> {
    let part = path.with_extension("json.part");
    std::fs::write(&part, serde_json::to_vec_pretty(value)?)?;
    std::fs::rename(&part, path)?;
    Ok(())
}

/// Builds the added matrix of `req.kind` into `req.out/matrix`, its counts
/// into `req.out/build.json`, and per-genome counts into `GENOMES_FILE`.
/// The added ksets are removed once the matrix is saved.
pub fn build(req: &Request) -> Result<Build> {
    let n = req.genomes.len();
    ensure!(n > 0, "no genomes");
    let kdir = req.out.join("ksets");
    std::fs::create_dir_all(&kdir)?;
    let space = KmerSpace::new(req.kind.space_k());
    let space31 = KmerSpace::new(K);
    let paths: Vec<PathBuf> = req
        .genomes
        .iter()
        .map(|g| kdir.join(format!("{g}.kset")))
        .collect();
    let mut b = Build {
        format: BUILD_FORMAT.to_string(),
        evaluation_revision: crate::analysis::EVALUATION_REVISION,
        kind: req.kind.name().to_string(),
        species: req.species.to_string(),
        genomes: n,
        space_k: space.k,
        matrix_support: MATRIX_SUPPORT,
        ..Default::default()
    };
    let mut rows: Vec<String> = Vec::with_capacity(n);
    match req.kind {
        Kind::Counts => {
            let got: Vec<usize> = req
                .genomes
                .par_iter()
                .zip(&paths)
                .map(|(g, p)| {
                    let v = count_values(&crate::bvbrc::fasta_path(req.fasta_dir, g), &space31)?;
                    kset::write(p, g, &space, &v)?;
                    Ok(v.len())
                })
                .collect::<Result<_>>()?;
            for (g, c) in req.genomes.iter().zip(&got) {
                rows.push(format!("{g}\t{c}"));
            }
            b.values_summed = got.iter().map(|&c| c as u64).sum();
        }
        Kind::Aa => {
            let dir = req
                .proteins
                .context("the aa build needs Prodigal's proteins")?;
            let got: Vec<(usize, usize)> = req
                .genomes
                .par_iter()
                .zip(&paths)
                .map(|(g, p)| {
                    let (v, proteins) = aa_values(
                        &crate::bvbrc::fasta_path(req.fasta_dir, g),
                        &dir.join(format!("{g}.faa")),
                        &space,
                    )?;
                    kset::write(p, g, &space, &v)?;
                    Ok((v.len(), proteins))
                })
                .collect::<Result<_>>()?;
            for (g, (c, p)) in req.genomes.iter().zip(&got) {
                rows.push(format!("{g}\t{c}\t{p}"));
            }
            b.values_summed = got.iter().map(|x| x.0 as u64).sum();
            b.proteins = Some(got.iter().map(|x| x.1 as u64).sum());
        }
        Kind::Unitigs => {
            let species_ksets: Vec<PathBuf> = req
                .genomes
                .iter()
                .map(|g| req.kset_dir.join(format!("{g}.kset")))
                .collect();
            let readers: Vec<KsetReader> = species_ksets
                .iter()
                .map(|p| KsetReader::open(p))
                .collect::<Result<_>>()?;
            ensure!(
                readers.iter().all(|r| r.k == K),
                "the species ksets are not all at k {K}"
            );
            let (table, equality) = match req.graph {
                Graph::Grm => {
                    let mut t = Table::new(union_values(&readers)?, space31.bits());
                    compact(&mut t, &space31)?;
                    let eq = match req.outside {
                        Some(path) => Some(compare(&t, &read_unitigs(path, &space31)?.0)),
                        None => None,
                    };
                    (t, eq)
                }
                Graph::Outside => {
                    let path = req
                        .outside
                        .context("an outside graph needs its unitigs FASTA")?;
                    (read_unitigs(path, &space31)?.0, None)
                }
            };
            let got: Vec<(usize, u64, u64)> = readers
                .par_iter()
                .zip(&paths)
                .map(|(r, p)| {
                    let (held, part, missing) = presence(&table, &r.read_all()?);
                    let mut v: Vec<u64> = held.iter().map(|&u| space.mix(u as u64)).collect();
                    v.sort_unstable();
                    kset::write(p, &r.name, &space, &v)?;
                    Ok((held.len(), part, missing))
                })
                .collect::<Result<_>>()?;
            for (g, (h, part, miss)) in req.genomes.iter().zip(&got) {
                rows.push(format!("{g}\t{h}\t{part}\t{miss}"));
            }
            b.values_summed = got.iter().map(|x| x.0 as u64).sum();
            b.graph = Some(
                match req.graph {
                    Graph::Grm => "grm",
                    Graph::Outside => "outside",
                }
                .to_string(),
            );
            b.graph_kmers = Some(table.values.len() as u64);
            b.unitigs = Some(table.sizes.len() as u64);
            b.pairs_held = b.values_summed.into();
            b.pairs_in_part = Some(got.iter().map(|x| x.1).sum());
            b.kmers_without_unitig = Some(got.iter().map(|x| x.2).sum());
            b.equality = equality;
        }
    }
    let header = match req.kind {
        Kind::Counts => "genome\tkmers_count_2",
        Kind::Aa => "genome\taa_10mers\tproteins",
        Kind::Unitigs => "genome\tunitigs_held\tunitigs_in_part\tkmers_without_unitig",
    };
    std::fs::write(
        req.out.join(GENOMES_FILE),
        format!("{header}\n{}\n", rows.join("\n")),
    )?;
    let m = matrix::build(&paths, MATRIX_SUPPORT, rayon::current_num_threads() * 2)?;
    ensure!(
        m.summary.genomes == req.genomes,
        "the added matrix's genomes are not the species matrix's, in order"
    );
    b.matrix_patterns = m.summary.n_patterns;
    b.matrix_values_kept = m.summary.kmers_kept;
    b.matrix_values_rare = m.summary.kmers_rare;
    b.matrix_values_core = m.summary.kmers_core;
    m.save(&req.out.join("matrix"))?;
    std::fs::remove_dir_all(&kdir)?;
    write_json(&req.out.join("build.json"), &b)?;
    Ok(b)
}

/// The species matrix's genomes in column order, from its `matrix.json`.
pub fn species_genomes(matrix_dir: &Path) -> Result<Vec<String>> {
    let s: MatrixSummary = serde_json::from_slice(
        &std::fs::read(matrix_dir.join("matrix.json"))
            .with_context(|| format!("no matrix in {}", matrix_dir.display()))?,
    )?;
    Ok(s.genomes)
}

/// Loads a build's matrix and its counts for `grm learn --exp-input-matrix`.
pub fn load(dir: &Path) -> Result<(FeatureMatrix, Build)> {
    let b: Build = serde_json::from_slice(
        &std::fs::read(dir.join("build.json"))
            .with_context(|| format!("no build.json in {}", dir.display()))?,
    )?;
    ensure!(
        b.format == BUILD_FORMAT,
        "{} is not a {BUILD_FORMAT} build",
        dir.display()
    );
    Ok((FeatureMatrix::load(&dir.join("matrix"))?, b))
}

/// The build's per-genome unitig counts (held whole, in part), by genome.
pub fn unitig_pairs(dir: &Path) -> Result<HashMap<String, (u64, u64)>> {
    let text = std::fs::read_to_string(dir.join(GENOMES_FILE))?;
    let mut out = HashMap::new();
    for line in text.lines().skip(1) {
        let cells: Vec<&str> = line.split('\t').collect();
        ensure!(cells.len() == 4, "a {GENOMES_FILE} row without 4 cells");
        out.insert(cells[0].to_string(), (cells[1].parse()?, cells[2].parse()?));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    //! SYNTHETIC tests: made-up sequences, proteins and graphs, no genome.
    use super::*;
    use std::collections::{BTreeMap, BTreeSet};

    fn rc(s: &str) -> String {
        s.bytes().rev().map(|b| complement(b) as char).collect()
    }

    /// A deterministic sequence of `n` bases.
    fn seq(n: usize, seed: u64) -> String {
        let mut x = seed.wrapping_mul(0x9e37_79b9_7f4a_7c15) | 1;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x >> 30) as usize & 3] as char
            })
            .collect()
    }

    fn tmp(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("grm-c3b-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn the_count_arm_keeps_the_31mers_seen_twice_as_a_recount_from_strings_does() {
        // SYNTHETIC: a contig with a 200-base repeat, a copy on the other
        // strand, an N run and a lower-case stretch. The recount below is
        // written apart: canonical strings, counted in a map.
        let a = seq(900, 1);
        let rep = seq(200, 2);
        let contig1 = format!("{a}{rep}{}NNNN{}", seq(300, 3), &rep[..120]);
        let contig2 = format!(
            "{}{}{}",
            seq(150, 4),
            rc(&rep[40..]),
            seq(80, 5).to_lowercase()
        );
        let dir = tmp("counts");
        let path = dir.join("g.fna");
        std::fs::write(&path, format!(">c1\n{contig1}\n>c2 note\n{contig2}\n")).unwrap();
        let space = KmerSpace::new(K);
        let got: BTreeSet<String> = count_values(&path, &space)
            .unwrap()
            .into_iter()
            .map(|v| space.decode(space.unmix(v)))
            .collect();
        let mut counts: BTreeMap<String, u32> = BTreeMap::new();
        for c in [&contig1, &contig2] {
            let up = c.to_uppercase();
            for w in up.as_bytes().windows(K) {
                if w.iter().all(|b| b"ACGT".contains(b)) {
                    let s = std::str::from_utf8(w).unwrap().to_string();
                    let r = rc(&s);
                    *counts.entry(if s <= r { s } else { r }).or_default() += 1;
                }
            }
        }
        let want: BTreeSet<String> = counts
            .iter()
            .filter(|(_, c)| **c >= 2)
            .map(|(s, _)| s.clone())
            .collect();
        assert!(want.len() > 150, "the fixture repeats k-mers");
        assert_eq!(got, want);
        let once = counts.len();
        assert!(once > want.len() + 1000, "a count column is not presence");
        std::fs::remove_dir_all(&dir).unwrap();
    }

    /// A SYNTHETIC gene: ATG, `codons` random sense codons, a stop.
    fn gene(codons: usize, seed: u64) -> String {
        let s = seq(codons * 3 + 60, seed);
        let mut out = String::from("ATG");
        for c in s.as_bytes().as_chunks::<3>().0 {
            if out.len() >= 3 + codons * 3 {
                break;
            }
            let c = std::str::from_utf8(&c[..]).unwrap();
            if !["TAA", "TAG", "TGA"].contains(&c) {
                out.push_str(c);
            }
        }
        out.push_str("TGA");
        out
    }

    #[test]
    fn a_protein_must_be_the_table_11_translation_at_its_coordinates() {
        // SYNTHETIC: one forward and one reverse gene on a contig, with
        // headers as Prodigal writes them.
        let g1 = gene(80, 7);
        let g2 = gene(60, 8);
        let pre = seq(50, 9);
        let mid = seq(40, 10);
        let contig = format!("{pre}{g1}{mid}{}{}", rc(&g2), seq(30, 11));
        let s1 = pre.len() + 1;
        let e1 = s1 + g1.len() - 1;
        let s2 = e1 + mid.len() + 1;
        let e2 = s2 + g2.len() - 1;
        let p1 = {
            let mut p = translate11(g1.as_bytes());
            p[0] = b'M';
            p
        };
        let p2 = translate11(g2.as_bytes());
        assert_eq!(p1.last(), Some(&b'*'));
        let h = |n: usize, s: usize, e: usize, st: &str, part: &str| {
            format!(
                "c1_{n} # {s} # {e} # {st} # ID=1_{n};partial={part};start_type=ATG;rbs_motif=None"
            )
        };
        let mut c1 = parse_header(&h(1, s1, e1, "1", "00")).unwrap();
        c1.protein = p1.clone();
        let mut c2 = parse_header(&h(2, s2, e2, "-1", "00")).unwrap();
        c2.protein = p2.clone();
        assert!(protein_matches(contig.as_bytes(), &c1));
        assert!(protein_matches(contig.as_bytes(), &c2));
        // Another frame.
        let mut f = c1.clone();
        f.start += 1;
        f.end += 1;
        assert!(!protein_matches(contig.as_bytes(), &f));
        // The other strand.
        let mut o = c1.clone();
        o.forward = false;
        assert!(!protein_matches(contig.as_bytes(), &o));
        // Table 4 reads TGA as W: a protein so translated is not table 11's.
        let mut t4 = c1.clone();
        t4.protein = p1
            .iter()
            .map(|&a| if a == b'*' { b'W' } else { a })
            .collect();
        assert!(!protein_matches(contig.as_bytes(), &t4));
        // The whole read: a .faa as Prodigal writes it, residues on lines of 60.
        let dir = tmp("aa");
        let wrap = |p: &[u8]| {
            p.chunks(60)
                .map(|c| String::from_utf8(c.to_vec()).unwrap())
                .collect::<Vec<_>>()
                .join("\n")
        };
        std::fs::write(dir.join("g.fna"), format!(">c1 synthetic\n{contig}\n")).unwrap();
        std::fs::write(
            dir.join("g.faa"),
            format!(
                ">{}\n{}\n>{}\n{}\n",
                h(1, s1, e1, "1", "00"),
                wrap(&p1),
                h(2, s2, e2, "-1", "00"),
                wrap(&p2)
            ),
        )
        .unwrap();
        let space = KmerSpace::new(AA_SPACE_K);
        let (v, n) = aa_values(&dir.join("g.fna"), &dir.join("g.faa"), &space).unwrap();
        assert_eq!(n, 2);
        let mut want = BTreeSet::new();
        for p in [&p1, &p2] {
            for w in p.windows(AA_K) {
                if !w.contains(&b'*') {
                    want.insert(w.to_vec());
                }
            }
        }
        assert_eq!(v.len(), want.len());
        // A frame-shifted protein file stops the build.
        std::fs::write(
            dir.join("bad.faa"),
            format!(">{}\n{}\n", h(1, s1 + 1, e1 + 1, "1", "00"), wrap(&p1)),
        )
        .unwrap();
        assert!(aa_values(&dir.join("g.fna"), &dir.join("bad.faa"), &space).is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn no_amino_acid_10mer_spans_a_stop_or_an_ambiguous_residue() {
        let mut got = Vec::new();
        for_each_aa_kmer(b"MKVLAAGHHHE*WWYYTTRRPPX", |v| got.push(v));
        // 11 residues before the stop give 2 windows; 10 after it and before
        // X give 1.
        let code = |s: &[u8]| {
            s.iter()
                .fold(0u64, |v, &b| (v << 5) | residue_code(b).unwrap())
        };
        assert_eq!(
            got,
            vec![
                code(b"MKVLAAGHHH"),
                code(b"KVLAAGHHHE"),
                code(b"WWYYTTRRPP")
            ]
        );
        let mut none = Vec::new();
        for_each_aa_kmer(b"MKVLA*AGHHHWBYZ", |v| none.push(v));
        assert!(none.is_empty());
    }

    /// SYNTHETIC ksets of `genomes`, at k 31, in a folder; returns readers.
    fn ksets(dir: &Path, genomes: &[(&str, String)]) -> Vec<KsetReader> {
        let space = KmerSpace::new(K);
        genomes
            .iter()
            .map(|(name, s)| {
                let path = dir.join(format!("{name}.fna"));
                std::fs::write(&path, format!(">c\n{s}\n")).unwrap();
                let (v, _) = kset::count_file(&path, &space, 1).unwrap();
                let kpath = dir.join(format!("{name}.kset"));
                kset::write(&kpath, name, &space, &v).unwrap();
                KsetReader::open(&kpath).unwrap()
            })
            .collect()
    }

    /// The unitigs of `t` as sets of decoded canonical 31-mers.
    fn parts(t: &Table, space: &KmerSpace) -> BTreeSet<BTreeSet<String>> {
        let mut by: BTreeMap<u32, BTreeSet<String>> = BTreeMap::new();
        for (i, &v) in t.values.iter().enumerate() {
            by.entry(t.unitig[i])
                .or_default()
                .insert(space.decode(space.unmix(v)));
        }
        by.into_values().collect()
    }

    #[test]
    fn the_graph_holds_single_genome_kmers_and_breaks_at_their_branch() {
        // SYNTHETIC: three genomes share a backbone; genome c alone has one
        // substitution in its middle, so the graph branches there, and its
        // 31 variant k-mers are single-genome k-mers that make a unitig.
        let back = seq(400, 21);
        let mut var = back.clone().into_bytes();
        var[200] = if var[200] == b'A' { b'C' } else { b'A' };
        let var = String::from_utf8(var).unwrap();
        let dir = tmp("graph");
        let readers = ksets(
            &dir,
            &[("a", back.clone()), ("b", back.clone()), ("c", var.clone())],
        );
        let space = KmerSpace::new(K);
        let mut t = Table::new(union_values(&readers).unwrap(), space.bits());
        let n_back = back.len() - K + 1;
        assert_eq!(
            t.values.len(),
            n_back + K,
            "every 31-mer, the 31 seen once included"
        );
        let count = compact(&mut t, &space).unwrap();
        // The backbone's left part, the branch, its right part, and the
        // variant path: 4 unitigs.
        assert_eq!(count, 4);
        let sizes: BTreeSet<u32> = t.sizes.iter().copied().collect();
        assert!(sizes.contains(&(K as u32)), "{sizes:?}");
        // Genome c holds the variant unitig whole; a and b hold none of it.
        let (held_a, part_a, miss_a) = presence(&t, &readers[0].read_all().unwrap());
        let (held_c, part_c, _) = presence(&t, &readers[2].read_all().unwrap());
        assert_eq!((held_a.len(), part_a, miss_a), (3, 0, 0));
        assert_eq!((held_c.len(), part_c), (3, 0));
        assert_ne!(held_a, held_c);
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn a_unitig_is_held_only_when_every_kmer_of_it_is() {
        // SYNTHETIC: genome b holds the first 300 bases of a's 500, so it
        // holds a's single unitig in part only.
        let a = seq(500, 31);
        let dir = tmp("held");
        let readers = ksets(&dir, &[("a", a.clone()), ("b", a[..300].to_string())]);
        let space = KmerSpace::new(K);
        let mut t = Table::new(union_values(&readers).unwrap(), space.bits());
        assert_eq!(compact(&mut t, &space).unwrap(), 1);
        let (held, part, _) = presence(&t, &readers[1].read_all().unwrap());
        assert_eq!((held.len(), part), (0, 1));
        let (held, part, _) = presence(&t, &readers[0].read_all().unwrap());
        assert_eq!((held.len(), part), (1, 0));
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn grm_and_an_outside_fasta_of_the_same_unitigs_are_equal_as_sets() {
        // SYNTHETIC: the outside unitigs written as sequences, one on the
        // reverse strand; the check finds every unitig in both. A FASTA that
        // merges two unitigs differs in both counts.
        let back = seq(300, 41);
        let mut var = back.clone().into_bytes();
        var[150] = if var[150] == b'G' { b'T' } else { b'G' };
        let var = String::from_utf8(var).unwrap();
        let dir = tmp("equal");
        let readers = ksets(&dir, &[("a", back.clone()), ("b", var.clone())]);
        let space = KmerSpace::new(K);
        let mut t = Table::new(union_values(&readers).unwrap(), space.bits());
        let count = compact(&mut t, &space).unwrap();
        assert_eq!(count, 4);
        let p = parts(&t, &space);
        // Each unitig's sequence, from the two genomes' paths.
        let unitigs = [
            back[..150].to_string(),
            back[150 - (K - 1)..151 + (K - 1)].to_string(),
            var[150 - (K - 1)..151 + (K - 1)].to_string(),
            back[151..].to_string(),
        ];
        let fasta: String = unitigs
            .iter()
            .enumerate()
            .map(|(i, u)| format!(">{i}\n{}\n", if i == 2 { rc(u) } else { u.clone() }))
            .collect();
        std::fs::write(dir.join("u.fa"), &fasta).unwrap();
        let (o, n) = read_unitigs(&dir.join("u.fa"), &space).unwrap();
        assert_eq!(n, 4);
        assert_eq!(parts(&o, &space), p);
        let eq = compare(&t, &o);
        assert_eq!(
            (
                eq.unitigs_both,
                eq.unitigs_grm_only,
                eq.unitigs_outside_only
            ),
            (4, 0, 0)
        );
        assert_eq!((eq.kmers_grm_only, eq.kmers_outside_only), (0, 0));
        let merged = format!(">0\n{back}\n>1\n{}\n", &var[150 - (K - 1)..151 + (K - 1)]);
        std::fs::write(dir.join("m.fa"), merged).unwrap();
        let (m, _) = read_unitigs(&dir.join("m.fa"), &space).unwrap();
        let eq = compare(&t, &m);
        assert_eq!(
            (
                eq.unitigs_both,
                eq.unitigs_grm_only,
                eq.unitigs_outside_only
            ),
            (1, 3, 1)
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
