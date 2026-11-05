//! Feature-block encoding for ablation arms. Each block's encoder is fitted on
//! the outer training samples and applied to every sample; each block is then
//! filtered and screened on its own, and columns of different blocks that are
//! identical on training samples collapse to one, which on each sample takes
//! the value of the highest-priority block that observes it (contract §18).
//!
//! `fit_encoder`, `transform` and `block_candidates` take no labels, so
//! nothing before `screen_block` can depend on phenotypes. Everything they fit
//! (cutpoints, medians, levels, which indicators exist, the support filter and
//! within-block duplicate collapse) comes from the outer training samples.
//! `screen_block` reads labels, so it runs once per fit, on that fit's rows.
//!
//! A sample can be missing from a block in three ways, each with its own
//! pattern:
//! - absent: the feature bit is 0;
//! - not assessable: the feature bits are 0 and the column's not-assessable
//!   indicator is 1;
//! - block unavailable (the block has no row for the sample; a `BlockSource`
//!   leaves such samples out of `build`): only the block's unavailable
//!   indicator is 1, and counts take the training median.
//!
//! An indicator column exists only when at least `min_support` training
//! samples are in its state. The one threshold covers both the per-column
//! not-assessable indicators and the block-unavailable indicator.
//!
//! Indicators are observations of availability (contract §20, worded as in
//! §21): the block-unavailable indicator observes every sample, and a
//! not-assessable indicator observes every sample where its block has a row.
//! On every sample an indicator observes, its value is an observation, not a
//! placeholder. Where its block has no row, a not-assessable indicator does
//! not observe: it reads 0 and gives way like a feature's placeholder. A
//! feature's value in a missing state is a placeholder: a 0, or the training
//! median's thermometer bits.
//!
//! A block with no encoding is unavailable, and `fit_encoder` returns an error
//! carrying the reason; the arm planner reports that arm unavailable with it.
//! That covers a block unavailable for every sample (`unavailable_reason`
//! set) and a block with no row for any training sample.

use crate::ablate::blocks::{BlockId, BlockSpec, RawBlock, RawColumn};
use crate::bits::{BitMatrix, Bits};
use crate::learn::{screen, training_candidates};
use anyhow::{Result, bail, ensure};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct EncodeConfig {
    /// Training samples a missing state needs before it gets an indicator
    /// column, and the support filter `block_candidates` applies.
    pub min_support: u32,
}

/// What an encoded column records about its raw column.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum ColumnRole {
    /// A binary feature, present and assessable.
    Present,
    /// A count at or above this cutpoint.
    AtLeast(u32),
    /// A categorical level seen in training.
    Level(String),
    /// The raw column could not be assessed. For a categorical column this
    /// includes a level training never saw.
    NotAssessable,
    /// The block has no row for the sample.
    BlockUnavailable,
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ColumnSource {
    pub block: BlockId,
    pub priority: u8,
    /// The raw column, or `None` for the block-unavailable indicator.
    pub column: Option<String>,
    pub role: ColumnRole,
    /// The samples whose value this source gives the encoded column: those
    /// it observes that no higher-priority source does. A within-block twin
    /// supplies none.
    pub supplied: u32,
}

/// Every block column an encoded column stands for, all identical on
/// training samples. `sources[0]` is the representative. After `assemble`,
/// each block's column follows in priority order, and each is followed by
/// the within-block twins `block_candidates` dropped, in block order.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ColumnProv {
    pub sources: Vec<ColumnSource>,
}

/// Encoded columns over the arm's samples: row r of `x` is described by
/// `prov[r]`, and `available[block]` marks the samples each source block has
/// a row for. `observed[block][column]` marks the samples where the block has
/// a row and that raw column is assessable, for the columns where that is
/// fewer than `available[block]`. Every block an `EncodedBlock` combines must
/// have been encoded over the same `samples`.
#[derive(Clone, Debug)]
pub struct EncodedBlock {
    pub x: BitMatrix,
    pub prov: Vec<ColumnProv>,
    pub available: BTreeMap<BlockId, Bits>,
    pub observed: BTreeMap<BlockId, BTreeMap<String, Bits>>,
}

impl EncodedBlock {
    fn empty(cols: usize) -> Self {
        EncodedBlock {
            x: BitMatrix::with_cols(cols),
            prov: Vec::new(),
            available: BTreeMap::new(),
            observed: BTreeMap::new(),
        }
    }

    /// The samples `source` observes: for a feature bit, where its block has
    /// a row and its raw column is assessable; for a not-assessable
    /// indicator, where its block has a row; for the block-unavailable
    /// indicator, every sample. An indicator's value is an observation of
    /// availability (contract §20), so each observes every sample it is
    /// defined on. Elsewhere a feature's value is a missing-state placeholder.
    pub fn observes(&self, source: &ColumnSource) -> Bits {
        let available = &self.available[&source.block];
        match (&source.role, &source.column) {
            (ColumnRole::BlockUnavailable, _) => Bits::ones(self.x.cols),
            (ColumnRole::NotAssessable, _) | (_, None) => available.clone(),
            (_, Some(column)) => self
                .observed
                .get(&source.block)
                .and_then(|columns| columns.get(column))
                .unwrap_or(available)
                .clone(),
        }
    }

    fn push(&mut self, row: &[u64], prov: ColumnProv) {
        self.x.push_row(row);
        self.prov.push(prov);
    }

    fn select(&self, rows: &[u32]) -> EncodedBlock {
        let mut out = EncodedBlock::empty(self.x.cols);
        out.available = self.available.clone();
        out.observed = self.observed.clone();
        for &r in rows {
            out.push(self.x.row(r as usize), self.prov[r as usize].clone());
        }
        out
    }
}

/// A block's encoding, fitted on training samples only.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BlockEncoder {
    pub block: BlockId,
    /// SHA-256 of the source `BlockSpec`'s canonical JSON; `transform`
    /// refuses a block whose spec hashes differently.
    pub spec_sha256: String,
    pub priority: u8,
    pub columns: Vec<ColumnFit>,
    /// Whether the block-unavailable indicator is emitted.
    pub unavailable_indicator: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ColumnFit {
    pub name: String,
    pub kind: FitKind,
    /// Whether the not-assessable indicator is emitted.
    pub indicator: bool,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum FitKind {
    Binary,
    /// Thermometer bits `>= t` for each cutpoint, strictly increasing. Tied
    /// quartiles merge, so a column can have fewer than three bits. The
    /// median imputes samples the block is unavailable for.
    Count {
        cutpoints: Vec<u32>,
        median: Option<u32>,
    },
    /// One bit per level seen in training, sorted by name.
    Categorical {
        levels: Vec<String>,
    },
}

impl BlockEncoder {
    fn prov(&self, column: Option<&str>, role: ColumnRole, supplied: u32) -> ColumnProv {
        ColumnProv {
            sources: vec![ColumnSource {
                block: self.block.clone(),
                priority: self.priority,
                column: column.map(str::to_string),
                role,
                supplied,
            }],
        }
    }
}

/// Where each requested sample sits in a raw block.
struct Align {
    rows: Vec<Option<usize>>,
    /// The block's samples are exactly the requested ones, in order.
    identity: bool,
    available: Bits,
}

impl Align {
    fn new(raw: &RawBlock, samples: &[String]) -> Result<Self> {
        let id = &raw.spec.id;
        if let Some(reason) = &raw.spec.unavailable_reason {
            bail!("block {id} is unavailable: {reason}");
        }
        ensure!(
            raw.names.len() == raw.cols.len(),
            "block {id}: {} names for {} columns",
            raw.names.len(),
            raw.cols.len()
        );
        let m = raw.sample_ids.len();
        for (name, col) in raw.names.iter().zip(&raw.cols) {
            let fits = match col {
                RawColumn::Binary {
                    present,
                    assessable,
                } => present.len() == m && assessable.len() == m,
                RawColumn::Count { values } => values.len() == m,
                RawColumn::Categorical { codes, levels } => {
                    codes.len() == m
                        && codes
                            .iter()
                            .flatten()
                            .all(|&c| usize::from(c) < levels.len())
                }
            };
            ensure!(
                fits,
                "block {id} column {name} does not match its {m} samples or its levels"
            );
        }
        let mut names: rustc_hash::FxHashSet<&str> = Default::default();
        for name in &raw.names {
            ensure!(
                names.insert(name.as_str()),
                "block {id} lists column {name} twice"
            );
        }
        let mut index: rustc_hash::FxHashMap<&str, usize> = Default::default();
        for (r, s) in raw.sample_ids.iter().enumerate() {
            ensure!(
                index.insert(s.as_str(), r).is_none(),
                "block {id} lists sample {s} twice"
            );
        }
        let mut requested: rustc_hash::FxHashSet<&str> = Default::default();
        for s in samples {
            ensure!(
                requested.insert(s.as_str()),
                "sample {s} is requested twice"
            );
        }
        let rows: Vec<Option<usize>> = samples
            .iter()
            .map(|s| index.get(s.as_str()).copied())
            .collect();
        let available = mask(samples.len(), |i| rows[i].is_some());
        Ok(Align {
            identity: raw.sample_ids == samples,
            rows,
            available,
        })
    }

    /// A raw bitset over the requested samples; 0 where the block has no row.
    fn gather(&self, b: &Bits) -> Bits {
        if self.identity {
            return b.clone();
        }
        mask(self.rows.len(), |i| {
            self.rows[i].is_some_and(|r| b.contains(r))
        })
    }

    /// Raw values of the masked samples that have a row.
    fn values<'a, T: Copy>(
        &'a self,
        of: &'a Bits,
        values: &'a [T],
    ) -> impl Iterator<Item = T> + 'a {
        of.iter_ones()
            .filter_map(move |i| self.rows[i].map(|r| values[r]))
    }
}

fn mask(n: usize, f: impl Fn(usize) -> bool) -> Bits {
    Bits::from_indices(n, (0..n).filter(|&i| f(i)))
}

/// Thermometer cutpoints from sorted training counts: the quartiles, strictly
/// increasing and above the minimum, so no bit is set on every assessable
/// training sample.
fn cutpoints(sorted: &[u32]) -> Vec<u32> {
    let Some(&min) = sorted.first() else {
        return Vec::new();
    };
    let mut cuts: Vec<u32> = (1..=3)
        .map(|k| sorted[k * sorted.len() / 4])
        .filter(|&t| t > min)
        .collect();
    cuts.dedup();
    cuts
}

/// The spec in the contract's canonical JSON (CONTRACT.md §7): H's
/// `BlockSpec::canonical_json`, which escapes every character outside
/// printable ASCII as Python's `ensure_ascii=True` does.
fn spec_json(spec: &BlockSpec) -> Result<String> {
    spec.canonical_json()
}

fn spec_sha256(spec: &BlockSpec) -> Result<String> {
    let digest = ring::digest::digest(&ring::digest::SHA256, spec_json(spec)?.as_bytes());
    Ok(digest.as_ref().iter().map(|b| format!("{b:02x}")).collect())
}

/// Fits `raw`'s encoding on the `train` samples of `samples`. Samples are
/// matched to the block's rows by ID; a sample with no row is
/// block-unavailable. Takes no labels, and reads no values outside `train`.
/// Errors if the whole block is unavailable, naming its reason, or if the
/// block has no row for any training sample.
pub fn fit_encoder(
    raw: &RawBlock,
    samples: &[String],
    train: &Bits,
    cfg: &EncodeConfig,
) -> Result<BlockEncoder> {
    ensure!(
        train.len() == samples.len(),
        "training mask covers {} samples, not {}",
        train.len(),
        samples.len()
    );
    let al = Align::new(raw, samples)?;
    let fit = train.and(&al.available);
    ensure!(
        fit.count() > 0,
        "block {} has no row for any training sample",
        raw.spec.id
    );
    let columns = raw
        .names
        .iter()
        .zip(&raw.cols)
        .map(|(name, col)| {
            let (kind, missing) = match col {
                RawColumn::Binary { assessable, .. } => {
                    (FitKind::Binary, fit.and_not(&al.gather(assessable)).count())
                }
                RawColumn::Count { values } => {
                    let mut seen: Vec<u32> = al.values(&fit, values).flatten().collect();
                    seen.sort_unstable();
                    let kind = FitKind::Count {
                        cutpoints: cutpoints(&seen),
                        median: seen.get(seen.len().saturating_sub(1) / 2).copied(),
                    };
                    (kind, fit.count() - seen.len() as u32)
                }
                RawColumn::Categorical { codes, levels } => {
                    let seen: std::collections::BTreeSet<&String> = al
                        .values(&fit, codes)
                        .flatten()
                        .map(|c| &levels[usize::from(c)])
                        .collect();
                    let missing = al.values(&fit, codes).filter(Option::is_none).count();
                    let levels = seen.into_iter().cloned().collect();
                    (FitKind::Categorical { levels }, missing as u32)
                }
            };
            ColumnFit {
                name: name.clone(),
                kind,
                indicator: missing >= cfg.min_support,
            }
        })
        .collect();
    Ok(BlockEncoder {
        block: raw.spec.id.clone(),
        spec_sha256: spec_sha256(&raw.spec)?,
        priority: raw.spec.priority,
        columns,
        unavailable_indicator: train.and_not(&al.available).count() >= cfg.min_support,
    })
}

/// Encodes `raw` over `samples` with a fitted encoder. `raw` must have the
/// fitted block's spec and columns, in order; its samples may differ from the
/// fit's. Per raw column the rows are its feature bits, then its
/// not-assessable indicator; the block-unavailable indicator comes last.
/// The result records the samples each raw column observes, for `assemble`.
pub fn transform(enc: &BlockEncoder, raw: &RawBlock, samples: &[String]) -> Result<EncodedBlock> {
    ensure!(
        raw.spec.id == enc.block,
        "encoder for block {} given block {}",
        enc.block,
        raw.spec.id
    );
    ensure!(
        spec_sha256(&raw.spec)? == enc.spec_sha256,
        "block {} spec differs from the one the encoder was fitted on",
        enc.block
    );
    ensure!(
        raw.names.len() == enc.columns.len()
            && raw
                .names
                .iter()
                .zip(&enc.columns)
                .all(|(n, c)| *n == c.name),
        "block {} columns differ from the fitted ones",
        enc.block
    );
    let al = Align::new(raw, samples)?;
    let n = samples.len();
    let mut out = EncodedBlock::empty(n);
    out.available
        .insert(enc.block.clone(), al.available.clone());
    for (fit, col) in enc.columns.iter().zip(&raw.cols) {
        let name = Some(fit.name.as_str());
        let (bits, not_assessable): (Vec<(Bits, ColumnRole)>, Bits) = match (&fit.kind, col) {
            (
                FitKind::Binary,
                RawColumn::Binary {
                    present,
                    assessable,
                },
            ) => {
                let assessable = al.gather(assessable);
                let bit = al.gather(present).and(&assessable);
                (
                    vec![(bit, ColumnRole::Present)],
                    al.available.and_not(&assessable),
                )
            }
            (FitKind::Count { cutpoints, median }, RawColumn::Count { values }) => {
                let value: Vec<Option<u32>> = al
                    .rows
                    .iter()
                    .map(|r| r.map_or(*median, |r| values[r]))
                    .collect();
                let bits = cutpoints
                    .iter()
                    .map(|&t| {
                        let bit = mask(n, |i| value[i].is_some_and(|v| v >= t));
                        (bit, ColumnRole::AtLeast(t))
                    })
                    .collect();
                (
                    bits,
                    mask(n, |i| al.rows[i].is_some() && value[i].is_none()),
                )
            }
            (FitKind::Categorical { levels: seen }, RawColumn::Categorical { codes, levels }) => {
                let level: Vec<Option<&str>> = al
                    .rows
                    .iter()
                    .map(|r| {
                        r.and_then(|r| codes[r])
                            .map(|c| levels[usize::from(c)].as_str())
                    })
                    .collect();
                let bits = seen
                    .iter()
                    .map(|l| {
                        let bit = mask(n, |i| level[i] == Some(l.as_str()));
                        (bit, ColumnRole::Level(l.clone()))
                    })
                    .collect();
                let unseen = |l: &str| seen.binary_search_by(|s| s.as_str().cmp(l)).is_err();
                (
                    bits,
                    mask(n, |i| al.rows[i].is_some() && level[i].is_none_or(unseen)),
                )
            }
            _ => bail!(
                "block {} column {} changed kind since the fit",
                enc.block,
                fit.name
            ),
        };
        let observed = al.available.and_not(&not_assessable);
        for (bit, role) in bits {
            out.push(bit.words(), enc.prov(name, role, observed.count()));
        }
        if fit.indicator {
            out.push(
                not_assessable.words(),
                enc.prov(name, ColumnRole::NotAssessable, al.available.count()),
            );
        }
        if observed != al.available {
            out.observed
                .entry(enc.block.clone())
                .or_default()
                .insert(fit.name.clone(), observed);
        }
    }
    if enc.unavailable_indicator {
        out.push(
            al.available.not().words(),
            enc.prov(None, ColumnRole::BlockUnavailable, n as u32),
        );
    }
    Ok(out)
}

/// The block's columns any fit in the outer fold may use, from
/// `learn::training_candidates` on the outer training mask: the support
/// filter and within-block duplicate collapse, first row wins. Reads no
/// labels. Each dropped twin is appended to its representative's provenance,
/// supplying no samples: the kept row is the representative's own.
///
/// Known limit (contract §19): within-block twins don't coalesce, so a
/// dropped twin's observation is lost where the kept column is not
/// assessable. The loss is the same in every arm that contains the block.
pub fn block_candidates(block: &EncodedBlock, train: &Bits, cfg: &EncodeConfig) -> EncodedBlock {
    let (_, rows) = training_candidates(&block.x, train, cfg.min_support);
    let mut out = block.select(&rows);
    let kept: rustc_hash::FxHashMap<Vec<u64>, usize> = rows
        .iter()
        .enumerate()
        .map(|(i, &r)| (train.and_slice(block.x.row(r as usize)).words().to_vec(), i))
        .collect();
    // A row dropped for support has no kept twin, since twins share support.
    for r in 0..block.x.rows {
        if let Some(&i) = kept.get(train.and_slice(block.x.row(r)).words())
            && rows[i] as usize != r
        {
            let twins = block.prov[r].sources.iter().map(|s| ColumnSource {
                supplied: 0,
                ..s.clone()
            });
            out.prov[i].sources.extend(twins);
        }
    }
    out
}

/// Keeps a candidate block's top `cap` columns by `learn::screen`, which
/// ranks them by chi-square against `y` on the `fit` samples. A block within
/// its cap is kept whole and `y` is not read. Survivors keep their block
/// order.
///
/// `fit` must be the rows of the exact fit that will use the screened block
/// (contract §16): each inner fit (the `fit_m` part of `tuning_holdout`, or
/// an inner CV training fold) gets its own call on its own rows, and only the
/// final outer fit uses the outer training mask. A block screened on the
/// outer training mask must never reach inner tuning, or validation labels
/// choose the capped block's columns.
pub fn screen_block(
    block: &EncodedBlock,
    y: &Bits,
    fit: &Bits,
    cap: Option<usize>,
) -> EncodedBlock {
    match cap {
        Some(cap) if block.x.rows > cap => {
            let mut top = screen(&block.x, y, fit, cap);
            top.sort_unstable();
            block.select(&top)
        }
        _ => block.clone(),
    }
}

/// Joins screened blocks into one arm matrix (contract §18). Columns of
/// different blocks that are identical on every training sample are twins,
/// and they collapse to one column whatever their blocks' availability. The
/// column from the lowest-priority-number block (C before H, P, N, X, K,
/// then block id) represents them, and the provenance lists every source in
/// that order.
///
/// On each sample the collapsed column takes the value of the
/// highest-priority source that observes it (`EncodedBlock::observes`).
/// Where none does, it keeps the representative's missing-state placeholder.
/// So an observed value gives way only to a higher-priority source's
/// observed value, never to a placeholder, whether the feature's block has
/// no row or its column is not assessable there. An indicator's 1 is an
/// observed value (contract §20) and wins by priority like any other, so a
/// higher-priority block's indicator can set a lower-priority twin's
/// test-sample value. Each source records how many samples it supplied; a
/// sample no source observes counts for none. Twins agree on training
/// samples, so training values, and everything fitted from them, are the
/// representative's.
///
/// Collapse reads training patterns, like `training_candidates`, and the
/// coalescing reads availability and assessability; neither carries labels.
/// Pass each block through `block_candidates` first: two training-identical
/// rows of one block would coalesce here, which contract §19 rules out. Call
/// this before any learner's own duplicate collapse, which keeps one twin's
/// row as it stands.
///
/// Never split a collapsed column back into per-block columns, for example to
/// report each block's importance: that puts the representative's placeholder
/// back on samples only a lower-priority twin observes (contract §18). The
/// provenance says which source supplied how many samples. The arm runner's
/// full order, `block_candidates`, then `screen_block` once per fit (§16),
/// then this, then the learners, is in `docs/ablation-runner.md`.
pub fn assemble(blocks: &[EncodedBlock], train: &Bits) -> EncodedBlock {
    let n = train.len();
    assert!(blocks.iter().all(|b| b.x.cols == n));
    let mut out = EncodedBlock::empty(n);
    for (id, available) in blocks.iter().flat_map(|b| &b.available) {
        let seen = out
            .available
            .entry(id.clone())
            .or_insert_with(|| available.clone());
        assert_eq!(seen, available, "block {id} comes with two availabilities");
    }
    for (id, columns) in blocks.iter().flat_map(|b| &b.observed) {
        for (column, observed) in columns {
            let seen = out
                .observed
                .entry(id.clone())
                .or_default()
                .entry(column.clone())
                .or_insert_with(|| observed.clone());
            assert_eq!(seen, observed, "block {id} column {column} comes twice");
        }
    }
    let mut order: Vec<(usize, usize)> = blocks
        .iter()
        .enumerate()
        .flat_map(|(b, block)| (0..block.x.rows).map(move |r| (b, r)))
        .collect();
    // Stable, so rows of one block keep their order.
    order.sort_by(|&(ab, ar), &(bb, br)| {
        let (a, b) = (
            &blocks[ab].prov[ar].sources[0],
            &blocks[bb].prov[br].sources[0],
        );
        (a.priority, &a.block).cmp(&(b.priority, &b.block))
    });
    // Twin groups by training pattern, each in priority order.
    let mut twins: Vec<Vec<(usize, usize)>> = Vec::new();
    let mut first: rustc_hash::FxHashMap<Vec<u64>, usize> = Default::default();
    for (b, r) in order {
        let key = train.and_slice(blocks[b].x.row(r)).words().to_vec();
        match first.entry(key) {
            std::collections::hash_map::Entry::Occupied(e) => twins[*e.get()].push((b, r)),
            std::collections::hash_map::Entry::Vacant(e) => {
                e.insert(twins.len());
                twins.push(vec![(b, r)]);
            }
        }
    }
    let row = |(b, r): (usize, usize)| Bits::from_words(n, blocks[b].x.row(r).to_vec());
    for group in twins {
        let mut value = row(group[0]);
        let mut filled = Bits::zeros(n);
        let mut prov = ColumnProv {
            sources: Vec::new(),
        };
        for (b, r) in group {
            let mut sources = blocks[b].prov[r].sources.clone();
            let observed = out.observes(&sources[0]);
            // What this source observes and no higher-priority one did.
            let supplies = observed.and_not(&filled);
            value = value.and_not(&supplies).or(&row((b, r)).and(&supplies));
            filled = filled.or(&observed);
            sources[0].supplied = supplies.count();
            prov.sources.extend(sources);
        }
        out.push(value.words(), prov);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ablate::blocks::BlockKind;
    use rand::seq::SliceRandom;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    const CFG: EncodeConfig = EncodeConfig { min_support: 2 };

    fn spec(id: &str, kind: BlockKind, priority: u8, screen_cap: Option<usize>) -> BlockSpec {
        BlockSpec {
            id: id.into(),
            kind,
            version: "test".into(),
            source_sha256: String::new(),
            priority,
            screen_cap,
            catalog_date: None,
            circularity_risk: false,
            unavailable_reason: None,
        }
    }

    fn ids(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("s{i}")).collect()
    }

    /// A block of always-assessable binary columns over `samples`.
    fn binary_block(spec: BlockSpec, samples: &[String], rows: &[Bits]) -> RawBlock {
        RawBlock {
            spec,
            sample_ids: samples.to_vec(),
            names: (0..rows.len()).map(|i| format!("f{i}")).collect(),
            cols: rows
                .iter()
                .map(|present| RawColumn::Binary {
                    present: present.clone(),
                    assessable: Bits::ones(present.len()),
                })
                .collect(),
        }
    }

    /// A block with one binary column, `f0`, over `samples`.
    fn binary(spec: BlockSpec, samples: &[String], present: Bits, assessable: Bits) -> RawBlock {
        RawBlock {
            spec,
            sample_ids: samples.to_vec(),
            names: vec!["f0".into()],
            cols: vec![RawColumn::Binary {
                present,
                assessable,
            }],
        }
    }

    /// Each source of encoded column `r`, with the samples it supplied.
    fn supplied(e: &EncodedBlock, r: usize) -> Vec<(&str, u32)> {
        e.prov[r]
            .sources
            .iter()
            .map(|s| (s.block.as_str(), s.supplied))
            .collect()
    }

    fn random_rows(rows: usize, n: usize, seed: u64) -> Vec<Bits> {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        (0..rows)
            .map(|_| Bits::from_indices(n, (0..n).filter(|_| rng.random_bool(0.5))))
            .collect()
    }

    fn encode(raw: &RawBlock, samples: &[String], train: &Bits) -> EncodedBlock {
        let enc = fit_encoder(raw, samples, train, &CFG).unwrap();
        transform(&enc, raw, samples).unwrap()
    }

    /// One binary, one count and one categorical column. Sample `i`'s values
    /// come from `seed(i)` alone; the block has rows for `rows`, in order.
    fn mixed_block(all: &[String], rows: &[usize], seed: impl Fn(usize) -> u64) -> RawBlock {
        let m = rows.len();
        let (mut present, mut assessable) = (Bits::zeros(m), Bits::zeros(m));
        let (mut counts, mut codes) = (Vec::new(), Vec::new());
        for (r, &i) in rows.iter().enumerate() {
            let mut rng = ChaCha8Rng::seed_from_u64(seed(i));
            if rng.random_bool(0.4) {
                present.insert(r);
            }
            if rng.random_bool(0.8) {
                assessable.insert(r);
            }
            counts.push(rng.random_bool(0.8).then(|| rng.random_range(0..12u32)));
            codes.push(rng.random_bool(0.8).then(|| rng.random_range(0..4u16)));
        }
        RawBlock {
            spec: spec("X", BlockKind::Context, 4, None),
            sample_ids: rows.iter().map(|&i| all[i].clone()).collect(),
            names: vec!["gene".into(), "copies".into(), "context".into()],
            cols: vec![
                RawColumn::Binary {
                    present,
                    assessable,
                },
                RawColumn::Count { values: counts },
                RawColumn::Categorical {
                    codes,
                    levels: ["chromosome", "plasmid", "phage", "novel"]
                        .map(String::from)
                        .to_vec(),
                },
            ],
        }
    }

    fn roles(e: &EncodedBlock) -> Vec<ColumnRole> {
        e.prov.iter().map(|p| p.sources[0].role.clone()).collect()
    }

    /// The encoded bits of sample `i`, one per row.
    fn pattern(e: &EncodedBlock, i: usize) -> Vec<bool> {
        (0..e.x.rows).map(|r| e.x.get(r, i)).collect()
    }

    #[test]
    fn the_fit_reads_training_samples_only() {
        let all = ids(60);
        let train = mask(60, |i| i % 4 != 0);
        // Samples 1, 2, 3 and 6 have no row in the block.
        let rows: Vec<usize> = (0..60).filter(|i| ![1, 2, 3, 6].contains(i)).collect();
        let fit = |raw: &RawBlock| {
            serde_json::to_vec(&fit_encoder(raw, &all, &train, &CFG).unwrap()).unwrap()
        };
        let base = fit(&mixed_block(&all, &rows, |i| i as u64));
        let enc: BlockEncoder = serde_json::from_slice(&base).unwrap();
        assert!(enc.unavailable_indicator && enc.columns.iter().all(|c| c.indicator));

        // New values for every test sample, test samples dropped from the
        // block, and the block's rows in another order.
        let test_seed = |i: usize| {
            if train.contains(i) {
                i as u64
            } else {
                1000 + i as u64
            }
        };
        let mut moved: Vec<usize> = rows
            .iter()
            .copied()
            .filter(|&i| train.contains(i) || i > 30)
            .collect();
        moved.shuffle(&mut ChaCha8Rng::seed_from_u64(1));
        assert_eq!(fit(&mixed_block(&all, &moved, test_seed)), base);

        // Training values do reach the fit.
        assert_ne!(fit(&mixed_block(&all, &rows, |i| 5000 + i as u64)), base);
    }

    #[test]
    fn indicators_cutpoints_and_the_median_are_decided_on_training_samples() {
        // s0..s5 train, s6..s9 test, min_support 2. Each missing state has
        // one training sample, so no indicator exists; counting test samples,
        // or the block-unavailable s5, would emit them. Training counts are
        // out of order.
        let all = ids(10);
        let rows = [0, 1, 2, 3, 4, 6, 7]; // no row for s5, s8 or s9
        let raw = RawBlock {
            spec: spec("X", BlockKind::Context, 4, None),
            sample_ids: rows.iter().map(|&i| all[i].clone()).collect(),
            names: vec!["gene".into(), "copies".into(), "context".into()],
            cols: vec![
                // Not assessable: s0, s6, s7.
                RawColumn::Binary {
                    present: Bits::from_indices(7, [3]),
                    assessable: Bits::from_indices(7, [1, 2, 3, 4]),
                },
                // Missing: s2, s6. Test s7 is 100.
                RawColumn::Count {
                    values: vec![Some(6), Some(2), None, Some(4), Some(1), None, Some(100)],
                },
                // Missing: s1, s6, s7.
                RawColumn::Categorical {
                    codes: vec![Some(0), None, Some(1), Some(0), Some(1), None, None],
                    levels: ["a", "b"].map(String::from).to_vec(),
                },
            ],
        };
        let train = mask(10, |i| i < 6);
        let enc = fit_encoder(&raw, &all, &train, &CFG).unwrap();
        assert!(!enc.unavailable_indicator);
        assert!(enc.columns.iter().all(|c| !c.indicator));
        // Training counts sort to 1, 2, 4, 6; with s7 the median would be 4.
        assert_eq!(
            enc.columns[1].kind,
            FitKind::Count {
                cutpoints: vec![2, 4, 6],
                median: Some(2)
            }
        );
    }

    #[test]
    fn a_count_without_training_values_keeps_no_row_apart_from_not_assessable() {
        // s0..s3 train. The block has no row for s2 or s3, and no count for
        // s0, s1 or s5, so there is no training median to impute.
        let all = ids(6);
        let raw = RawBlock {
            spec: spec("N", BlockKind::Count, 3, None),
            sample_ids: [0, 1, 4, 5].map(|i| all[i].clone()).to_vec(),
            names: vec!["copies".into()],
            cols: vec![RawColumn::Count {
                values: vec![None, None, Some(3), None],
            }],
        };
        let train = mask(6, |i| i < 4);
        let e = encode(&raw, &all, &train);
        assert_eq!(
            roles(&e),
            [ColumnRole::NotAssessable, ColumnRole::BlockUnavailable]
        );
        assert_eq!(pattern(&e, 0), [true, false]);
        assert_eq!(pattern(&e, 2), [false, true]);
        assert_eq!(pattern(&e, 4), [false, false]);
    }

    #[test]
    fn blocks_with_no_encoding_and_repeated_samples_are_refused() {
        let all = ids(10);
        let train = mask(10, |i| i < 8);
        let refused = |raw: &RawBlock, samples: &[String]| {
            fit_encoder(raw, samples, &mask(samples.len(), |i| i < 8), &CFG)
                .unwrap_err()
                .to_string()
        };
        let unavailable = RawBlock {
            spec: BlockSpec {
                unavailable_reason: Some("no count sidecar".into()),
                ..spec("N", BlockKind::Count, 3, None)
            },
            sample_ids: all.clone(),
            names: vec!["copies".into()],
            cols: vec![RawColumn::Count {
                values: vec![Some(1); 10],
            }],
        };
        assert_eq!(
            refused(&unavailable, &all),
            "block N is unavailable: no count sidecar"
        );
        // The reviewer's K block with no rows, and one with rows only for test
        // samples.
        let k = spec("K", BlockKind::Anonymous, 5, Some(4000));
        let empty = binary_block(k.clone(), &[], &[Bits::zeros(0)]);
        let test_only = binary_block(k, &all[8..], &[Bits::ones(2)]);
        for raw in [&empty, &test_only] {
            assert_eq!(
                refused(raw, &all),
                "block K has no row for any training sample"
            );
        }
        // One training row is enough.
        let one_row = binary_block(
            spec("P", BlockKind::Presence, 2, None),
            &[0, 8, 9].map(|i| all[i].clone()),
            &[Bits::from_indices(3, [0])],
        );
        assert!(fit_encoder(&one_row, &all, &train, &CFG).is_ok());
        let mut repeated = all.clone();
        repeated[1] = all[0].clone();
        let raw = binary_block(
            spec("P", BlockKind::Presence, 2, None),
            &all,
            &[Bits::ones(10)],
        );
        assert_eq!(refused(&raw, &repeated), "sample s0 is requested twice");
        assert!(fit_encoder(&raw, &all, &train, &CFG).is_ok());
        let mut names_twice = raw.clone();
        names_twice.names.push("f0".into());
        names_twice.cols.push(raw.cols[0].clone());
        assert_eq!(refused(&names_twice, &all), "block P lists column f0 twice");
    }

    #[test]
    fn the_encoder_pins_its_spec_by_canonical_hash() {
        let k = BlockSpec {
            version: "k=31 min_support=1 patterns=3 genomes=5".into(),
            source_sha256: "0123456789abcdef".repeat(4),
            ..spec("K", BlockKind::Anonymous, 5, Some(4000))
        };
        // The bytes H's `specs_serialise_to_the_contract_canonical_json`
        // expects for this spec.
        assert_eq!(
            spec_json(&k).unwrap(),
            r#"{"catalog_date":null,"circularity_risk":false,"id":"K","kind":"Anonymous","priority":5,"screen_cap":4000,"source_sha256":"0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef","unavailable_reason":null,"version":"k=31 min_support=1 patterns=3 genomes=5"}"#
        );
        let all = ids(6);
        let raw = binary_block(k.clone(), &all, &random_rows(3, 6, 8));
        let enc = fit_encoder(&raw, &all, &mask(6, |i| i < 4), &CFG).unwrap();
        assert_eq!(enc.spec_sha256, spec_sha256(&k).unwrap());
        assert!(transform(&enc, &raw, &all).is_ok());
        // Same id and columns, another matrix.
        let other = RawBlock {
            spec: BlockSpec {
                source_sha256: "f".repeat(64),
                ..k
            },
            ..raw.clone()
        };
        assert_eq!(
            transform(&enc, &other, &all).unwrap_err().to_string(),
            "block K spec differs from the one the encoder was fitted on"
        );
    }

    #[test]
    fn encoding_is_label_free_and_only_screening_reads_labels() {
        let all = ids(64);
        let train = mask(64, |i| i < 48);
        let y = mask(64, |i| i % 3 == 0);
        let mut perm: Vec<usize> = (0..64).collect();
        perm.shuffle(&mut ChaCha8Rng::seed_from_u64(2));
        let y_perm = mask(64, |i| y.contains(perm[i]));
        let k = binary_block(
            spec("K", BlockKind::Anonymous, 5, Some(10)),
            &all,
            &random_rows(200, 64, 3),
        );
        let x = mixed_block(&all, &(0..60).collect::<Vec<_>>(), |i| i as u64);
        // The arm's path up to screening: no step can be handed `y`.
        let arm = |_y: &Bits| {
            [&k, &x].map(|raw| block_candidates(&encode(raw, &all, &train), &train, &CFG))
        };
        for (a, b) in arm(&y).iter().zip(&arm(&y_perm)) {
            assert_eq!(a.x.data, b.x.data);
            assert_eq!(a.prov, b.prov);
        }
        // The permutation does change what screening keeps.
        let kept = |y: &Bits| screen_block(&arm(y)[0], y, &train, Some(10)).prov;
        assert_ne!(kept(&y), kept(&y_perm));
    }

    #[test]
    fn thermometer_cutpoints_are_training_quartiles_and_the_top_bin_is_open() {
        assert_eq!(cutpoints(&[1, 2, 3, 4, 5, 6, 7, 8]), [3, 5, 7]);
        // Tied quartiles merge, and a cutpoint at the minimum is dropped.
        assert_eq!(cutpoints(&[0, 0, 0, 0, 0, 0, 1, 5]), [1]);
        assert_eq!(cutpoints(&[4, 4, 4]), Vec::<u32>::new());
        assert_eq!(cutpoints(&[]), Vec::<u32>::new());
        for seed in 0..50 {
            let mut rng = ChaCha8Rng::seed_from_u64(seed);
            let mut v: Vec<u32> = (0..rng.random_range(1..40))
                .map(|_| rng.random_range(0..20))
                .collect();
            v.sort_unstable();
            let c = cutpoints(&v);
            assert!(c.windows(2).all(|w| w[0] < w[1]) && c.iter().all(|&t| t > v[0]));
            assert_eq!(cutpoints(&v), c);
        }

        // Training counts 1..=8 on s0..s7; test s8 is far above the range and
        // s9 below it.
        let all = ids(10);
        let values = vec![
            Some(1),
            Some(2),
            Some(3),
            Some(4),
            Some(5),
            Some(6),
            Some(7),
            Some(8),
            Some(500),
            Some(0),
        ];
        let raw = RawBlock {
            spec: spec("N", BlockKind::Count, 3, None),
            sample_ids: all.clone(),
            names: vec!["copies".into()],
            cols: vec![RawColumn::Count { values }],
        };
        let train = mask(10, |i| i < 8);
        let enc = fit_encoder(&raw, &all, &train, &CFG).unwrap();
        // Shuffling the block's rows gives the same fit.
        let mut order: Vec<usize> = (0..10).collect();
        order.shuffle(&mut ChaCha8Rng::seed_from_u64(4));
        let RawColumn::Count { values } = &raw.cols[0] else {
            unreachable!()
        };
        let shuffled = RawBlock {
            sample_ids: order.iter().map(|&i| all[i].clone()).collect(),
            cols: vec![RawColumn::Count {
                values: order.iter().map(|&i| values[i]).collect(),
            }],
            ..raw.clone()
        };
        assert_eq!(fit_encoder(&shuffled, &all, &train, &CFG).unwrap(), enc);

        let e = transform(&enc, &raw, &all).unwrap();
        assert_eq!(roles(&e), [3, 5, 7].map(ColumnRole::AtLeast));
        assert_eq!(pattern(&e, 7), [true; 3]); // 8, the training maximum
        assert_eq!(pattern(&e, 8), [true; 3]); // 500
        assert_eq!(pattern(&e, 3), [true, false, false]); // 4
        assert_eq!(pattern(&e, 9), [false; 3]); // 0
    }

    #[test]
    fn the_three_missing_states_differ_and_indicators_need_min_support() {
        // s0..s7 train, s8 and s9 test. s2 and s3 cannot be assessed; the
        // block has no row for s4 or s5.
        let all = ids(10);
        let rows = [9, 0, 1, 2, 3, 6, 7, 8];
        let at = |idx: &[usize]| {
            Bits::from_indices(
                rows.len(),
                (0..rows.len()).filter(|&r| idx.contains(&rows[r])),
            )
        };
        let raw = RawBlock {
            spec: spec("P", BlockKind::Presence, 2, None),
            sample_ids: rows.iter().map(|&i| all[i].clone()).collect(),
            names: vec!["gene".into()],
            cols: vec![RawColumn::Binary {
                present: at(&[0, 1, 2, 9]),
                assessable: at(&[0, 1, 6, 7, 8, 9]),
            }],
        };
        let train = mask(10, |i| i < 8);
        let e = encode(&raw, &all, &train);
        assert_eq!(
            roles(&e),
            [
                ColumnRole::Present,
                ColumnRole::NotAssessable,
                ColumnRole::BlockUnavailable
            ]
        );
        assert_eq!(pattern(&e, 0), [true, false, false]); // present
        assert_eq!(pattern(&e, 6), [false, false, false]); // absent
        assert_eq!(pattern(&e, 2), [false, true, false]); // not assessable, present bit ignored
        assert_eq!(pattern(&e, 4), [false, false, true]); // block unavailable
        assert_eq!(pattern(&e, 9), [true, false, false]); // a test sample

        // Two training samples per state: below a min_support of 3, neither
        // indicator exists and the states read as absent.
        let enc = fit_encoder(&raw, &all, &train, &EncodeConfig { min_support: 3 }).unwrap();
        let e = transform(&enc, &raw, &all).unwrap();
        assert_eq!(roles(&e), [ColumnRole::Present]);
        assert_eq!(pattern(&e, 2), pattern(&e, 6));
        assert_eq!(pattern(&e, 4), pattern(&e, 6));
    }

    #[test]
    fn unavailable_counts_take_the_training_median_and_unseen_levels_are_not_assessable() {
        // s0..s7 train, s8..s9 test; the block has no row for s6, s7 or s9.
        let all = ids(10);
        let levels = ["chromosome", "plasmid", "phage"]
            .map(String::from)
            .to_vec();
        let raw = RawBlock {
            spec: spec("X", BlockKind::Context, 4, None),
            sample_ids: ids(9)[..6].iter().chain([&all[8]]).cloned().collect(),
            names: vec!["copies".into(), "context".into()],
            cols: vec![
                RawColumn::Count {
                    values: vec![Some(1), Some(3), Some(3), Some(9), None, None, Some(30)],
                },
                // Training sees chromosome and plasmid; s8 has phage.
                RawColumn::Categorical {
                    codes: vec![Some(0), Some(1), Some(0), None, None, Some(1), Some(2)],
                    levels,
                },
            ],
        };
        let train = mask(10, |i| i < 8);
        let enc = fit_encoder(&raw, &all, &train, &CFG).unwrap();
        assert_eq!(
            enc.columns[0].kind,
            FitKind::Count {
                cutpoints: vec![3, 9],
                median: Some(3)
            }
        );
        let e = transform(&enc, &raw, &all).unwrap();
        let level = |l: &str| ColumnRole::Level(l.into());
        assert_eq!(
            roles(&e),
            [
                ColumnRole::AtLeast(3),
                ColumnRole::AtLeast(9),
                ColumnRole::NotAssessable,
                level("chromosome"),
                level("plasmid"),
                ColumnRole::NotAssessable,
                ColumnRole::BlockUnavailable,
            ]
        );
        // Unavailable: the median's thermometer bits (3), no column indicators.
        assert_eq!(
            pattern(&e, 6),
            [true, false, false, false, false, false, true]
        );
        // A missing count and a missing level.
        assert_eq!(
            pattern(&e, 4),
            [false, false, true, false, false, true, false]
        );
        // An out-of-range count and a level training never saw.
        assert_eq!(
            pattern(&e, 8),
            [true, true, false, false, false, true, false]
        );

        // A raw block whose level codes are numbered differently encodes by
        // name.
        let RawColumn::Categorical { codes, .. } = &raw.cols[1] else {
            unreachable!()
        };
        let mut renumbered = raw.clone();
        renumbered.cols[1] = RawColumn::Categorical {
            codes: codes.iter().map(|c| c.map(|c| 2 - c)).collect(),
            levels: ["phage", "plasmid", "chromosome"]
                .map(String::from)
                .to_vec(),
        };
        assert_eq!(transform(&enc, &renumbered, &all).unwrap().x.data, e.x.data);
    }

    #[test]
    fn screening_respects_caps_and_keeps_small_blocks_whole() {
        let all = ids(64);
        let train = mask(64, |i| i < 48);
        let y = mask(64, |i| i % 2 == 0);
        let mut rows = random_rows(50, 64, 5);
        // f50 is f3's twin on training samples only: test samples s60 and
        // s61 are flipped. The first wins and the twin joins its provenance.
        // f51 is present on test samples only.
        let mut twin = rows[3].clone();
        for i in [60, 61] {
            if twin.contains(i) {
                twin.remove(i);
            } else {
                twin.insert(i);
            }
        }
        assert_ne!(twin, rows[3]);
        assert_eq!(train.and(&twin), train.and(&rows[3]));
        rows.push(twin);
        rows.push(Bits::from_indices(64, [60, 61, 62]));
        let e = encode(
            &binary_block(spec("P", BlockKind::Presence, 2, None), &all, &rows),
            &all,
            &train,
        );
        let cands = block_candidates(&e, &train, &CFG);
        assert_eq!(cands.x.rows, 50);
        let columns = |p: &ColumnProv| -> Vec<String> {
            p.sources
                .iter()
                .map(|s| s.column.clone().unwrap())
                .collect()
        };
        assert_eq!(columns(&cands.prov[3]), ["f3", "f50"]);
        assert_eq!(supplied(&cands, 3), [("P", 64), ("P", 0)]);
        assert!(
            cands
                .prov
                .iter()
                .all(|p| p.sources.len() == 1 || columns(p)[0] == "f3")
        );
        assert!(
            !cands
                .prov
                .iter()
                .any(|p| columns(p).contains(&"f51".into()))
        );

        assert_eq!(screen_block(&cands, &y, &train, Some(5)).x.rows, 5);
        // A block within its cap is kept whole without reading `y`: an empty
        // `y` would fail `screen`'s length check against `fit`.
        for cap in [None, Some(50), Some(500)] {
            let kept = screen_block(&cands, &Bits::zeros(0), &train, cap);
            assert_eq!(kept.prov, cands.prov);
            assert_eq!(kept.x.data, cands.x.data);
        }

        // An all-blocks arm, K over its 4,000 cap.
        let k_rows = random_rows(5000, 64, 6);
        let mut blocks = vec![binary_block(
            spec("K", BlockKind::Anonymous, 5, Some(4000)),
            &all,
            &k_rows,
        )];
        for (i, (id, kind)) in [
            ("C", BlockKind::CuratedExact),
            ("H", BlockKind::Hierarchical),
            ("P", BlockKind::Presence),
            ("N", BlockKind::Count),
            ("X", BlockKind::Context),
        ]
        .into_iter()
        .enumerate()
        {
            // Every other block repeats some K columns.
            let mut own = random_rows(30, 64, 10 + i as u64);
            own.extend(k_rows[i * 7..i * 7 + 5].iter().cloned());
            blocks.push(binary_block(spec(id, kind, i as u8, Some(20)), &all, &own));
        }
        // Screened for the final outer fit, on the outer training mask.
        let screened: Vec<EncodedBlock> = blocks
            .iter()
            .map(|raw| {
                let cands = block_candidates(&encode(raw, &all, &train), &train, &CFG);
                screen_block(&cands, &y, &train, raw.spec.screen_cap)
            })
            .collect();
        assert_eq!(screened[0].x.rows, 4000);
        assert!(screened[1..].iter().all(|b| b.x.rows == 20));
        let arm = assemble(&screened, &train);
        let from_k = arm
            .prov
            .iter()
            .filter(|p| p.sources.iter().any(|s| s.block == "K"))
            .count();
        assert!(from_k <= 4000);
        assert_eq!(
            arm.prov
                .iter()
                .filter(|p| p.sources[0].block == "K")
                .count()
                + 5 * 20,
            arm.x.rows
        );
    }

    #[test]
    fn k_alone_matches_the_base_pipeline_and_each_fit_screens_on_its_own_rows() {
        let all = ids(64);
        let train = mask(64, |i| i % 4 != 1);
        // An inner fit's rows, as `tuning_holdout`'s `fit_m` would give.
        let fit_m = mask(64, |i| i % 4 != 1 && i % 5 != 0);
        let y = mask(64, |i| i % 3 != 0);
        let x_rows = random_rows(5000, 64, 7);
        let mut x = BitMatrix::with_cols(64);
        for r in &x_rows {
            x.push_row(r.words());
        }
        let raw = binary_block(
            spec("K", BlockKind::Anonymous, 5, Some(4000)),
            &all,
            &x_rows,
        );
        let e = encode(&raw, &all, &train);
        assert_eq!(e.x.data, x.data);
        let cands = block_candidates(&e, &train, &CFG);

        // The base learners: candidates on the outer training mask, then
        // `screen` on the rows of the fit being made.
        let (cx, rows) = training_candidates(&x, &train, CFG.min_support);
        assert!(cx.rows > 4000);
        assert_eq!(cands.x.data, cx.data);
        let mut kept = Vec::new();
        for fit in [&train, &fit_m] {
            let mut top: Vec<u32> = screen(&cx, &y, fit, 4000)
                .into_iter()
                .map(|f| rows[f as usize])
                .collect();
            top.sort_unstable();
            let mut want = BitMatrix::with_cols(64);
            for &r in &top {
                want.push_row(x.row(r as usize));
            }
            let ours = screen_block(&cands, &y, fit, raw.spec.screen_cap);
            assert_eq!(ours.x.data, want.data);
            let names: Vec<Option<String>> = top.iter().map(|r| Some(format!("f{r}"))).collect();
            assert_eq!(
                ours.prov
                    .iter()
                    .map(|p| p.sources[0].column.clone())
                    .collect::<Vec<_>>(),
                names
            );
            kept.push(ours.prov);
        }
        // The outer fit's screen and an inner fit's screen keep different
        // columns, so neither can stand in for the other.
        assert_ne!(kept[0], kept[1]);
    }

    #[test]
    fn cross_block_collapse_keeps_the_higher_priority_column_and_lists_every_source() {
        let all = ids(12);
        let train = mask(12, |i| i < 9);
        let shared = Bits::from_indices(12, [0, 1, 2, 9]);
        // Identical to `shared` on training samples, different on test.
        let shared_on_train = Bits::from_indices(12, [0, 1, 2, 10, 11]);
        let k_only = Bits::from_indices(12, [3, 4, 5]);
        let c = binary_block(
            spec("C", BlockKind::CuratedExact, 0, None),
            &all,
            std::slice::from_ref(&shared),
        );
        let p = binary_block(
            spec("P", BlockKind::Presence, 2, None),
            &all,
            &[Bits::from_indices(12, [6, 7]), shared_on_train],
        );
        // K's twin also differs from C's column on test samples, so the
        // representative's own row is the one kept.
        let k = binary_block(
            spec("K", BlockKind::Anonymous, 5, Some(4000)),
            &all,
            &[k_only.clone(), Bits::from_indices(12, [0, 1, 2, 11])],
        );
        // Input order is lowest priority first; priority still decides.
        let arm = assemble(&[&k, &p, &c].map(|raw| encode(raw, &all, &train)), &train);
        let sources = |r: usize| -> Vec<(String, String)> {
            arm.prov[r]
                .sources
                .iter()
                .map(|s| (s.block.clone(), s.column.clone().unwrap()))
                .collect()
        };
        let pair = |block: &str, col: &str| (block.to_string(), col.to_string());
        assert_eq!(arm.x.rows, 3);
        assert_eq!(
            sources(0),
            [pair("C", "f0"), pair("P", "f1"), pair("K", "f1")]
        );
        // C observes every sample, so its own row is the value.
        assert_eq!(supplied(&arm, 0), [("C", 12), ("P", 0), ("K", 0)]);
        assert_eq!(arm.x.row(0), Bits::from_indices(12, [0, 1, 2, 9]).words());
        assert_eq!(sources(1), [pair("P", "f0")]);
        assert_eq!(sources(2), [pair("K", "f0")]);
        assert_eq!(arm.x.row(2), k_only.words());
    }

    #[test]
    fn twins_with_different_availability_collapse_and_the_learner_reads_the_observed_value() {
        // Probe P1, the reviewer's case: s0..s7 train. H has no row for test
        // s8, and no training sample lacks a row, so H has no block
        // indicator. P observes the gene on s8.
        let all = ids(10);
        let train = mask(10, |i| i < 8);
        let h_ids: Vec<String> = all.iter().filter(|s| *s != "s8").cloned().collect();
        let h = binary(
            spec("H", BlockKind::Hierarchical, 1, None),
            &h_ids,
            Bits::from_indices(9, [0, 1, 2]),
            Bits::ones(9),
        );
        let p = binary(
            spec("P", BlockKind::Presence, 2, None),
            &all,
            Bits::from_indices(10, [0, 1, 2, 8]),
            Bits::ones(10),
        );
        let blocks = [&h, &p].map(|raw| encode(raw, &all, &train));
        assert_eq!(roles(&blocks[0]), [ColumnRole::Present]);
        let arm = assemble(&blocks, &train);
        assert_eq!(arm.x.rows, 1);
        assert_eq!(supplied(&arm, 0), [("H", 9), ("P", 1)]);
        assert_eq!(arm.available.keys().collect::<Vec<_>>(), ["H", "P"]);

        // What a base learner reads: its own candidates on the training mask.
        let (cx, rows) = training_candidates(&arm.x, &train, CFG.min_support);
        assert_eq!(rows, [0]);
        let p_only = assemble(&[encode(&p, &all, &train)], &train);
        let (px, _) = training_candidates(&p_only.x, &train, CFG.min_support);
        assert!(cx.get(0, 8), "s8 reads P's observed 1");
        assert_eq!(cx.row(0), px.row(0));
    }

    #[test]
    fn a_not_assessable_representative_takes_the_value_of_an_observing_twin() {
        // Probe P2: C cannot assess test s8, where P observes the gene. No
        // training sample is a no-call, so C has no not-assessable indicator.
        let all = ids(10);
        let train = mask(10, |i| i < 8);
        let mut c_assessable = Bits::ones(10);
        c_assessable.remove(8);
        let c = binary(
            spec("C", BlockKind::CuratedExact, 0, None),
            &all,
            Bits::from_indices(10, [0, 1, 2, 8]),
            c_assessable,
        );
        let p = binary(
            spec("P", BlockKind::Presence, 2, None),
            &all,
            Bits::from_indices(10, [0, 1, 2, 8]),
            Bits::ones(10),
        );
        let blocks = [&p, &c].map(|raw| encode(raw, &all, &train));
        assert_eq!(roles(&blocks[1]), [ColumnRole::Present]);
        assert!(!blocks[1].x.get(0, 8));
        let arm = assemble(&blocks, &train);
        assert_eq!(arm.x.rows, 1);
        assert_eq!(supplied(&arm, 0), [("C", 9), ("P", 1)]);
        assert_eq!(arm.x.row(0), Bits::from_indices(10, [0, 1, 2, 8]).words());
    }

    #[test]
    fn where_twins_both_observe_a_sample_the_higher_priority_value_wins() {
        // s0..s7 train; C, H and P agree on every training sample.
        // s8: C cannot assess it; H observes 0 and P observes 1.
        // s9: C has no row; H observes 1 and P observes 0.
        let all = ids(10);
        let train = mask(10, |i| i < 8);
        let mut c_assessable = Bits::ones(9);
        c_assessable.remove(8);
        let c = binary(
            spec("C", BlockKind::CuratedExact, 0, None),
            &all[..9],
            Bits::from_indices(9, [0, 1, 2, 8]),
            c_assessable,
        );
        let h = binary(
            spec("H", BlockKind::Hierarchical, 1, None),
            &all,
            Bits::from_indices(10, [0, 1, 2, 9]),
            Bits::ones(10),
        );
        let p = binary(
            spec("P", BlockKind::Presence, 2, None),
            &all,
            Bits::from_indices(10, [0, 1, 2, 8]),
            Bits::ones(10),
        );
        let blocks = [&p, &c, &h].map(|raw| encode(raw, &all, &train));
        assert!(blocks.iter().all(|b| b.x.rows == 1));
        let arm = assemble(&blocks, &train);
        assert_eq!(arm.x.rows, 1);
        assert!(!arm.x.get(0, 8), "H's 0 beats P's 1");
        assert!(arm.x.get(0, 9), "H's 1 beats P's 0");
        assert_eq!(arm.x.row(0), Bits::from_indices(10, [0, 1, 2, 9]).words());
        assert_eq!(supplied(&arm, 0), [("C", 8), ("H", 2), ("P", 0)]);
    }

    #[test]
    fn columns_that_differ_on_any_training_sample_never_collapse() {
        // s0..s7 train. P's f0 is C's twin and the only one to observe test
        // s8. P's f1 differs from it on training s7 alone, and f2 has its
        // training support but another pattern.
        let all = ids(10);
        let train = mask(10, |i| i < 8);
        let c = binary(
            spec("C", BlockKind::CuratedExact, 0, None),
            &all[..8],
            Bits::from_indices(8, [0, 1, 2]),
            Bits::ones(8),
        );
        let p_rows = [
            Bits::from_indices(10, [0, 1, 2, 8]),
            Bits::from_indices(10, [0, 1, 2, 7, 8]),
            Bits::from_indices(10, [0, 1, 3, 8]),
        ];
        let p = binary_block(spec("P", BlockKind::Presence, 2, None), &all, &p_rows);
        let arm = assemble(&[&c, &p].map(|raw| encode(raw, &all, &train)), &train);
        assert_eq!(arm.x.rows, 3);
        assert_eq!(supplied(&arm, 0), [("C", 8), ("P", 2)]);
        assert_eq!(arm.x.row(0), p_rows[0].words());
        for (r, want) in [(1, &p_rows[1]), (2, &p_rows[2])] {
            assert_eq!(supplied(&arm, r), [("P", 10)]);
            assert_eq!(arm.x.row(r), want.words());
        }
    }

    #[test]
    #[should_panic(expected = "block P comes with two availabilities")]
    fn one_block_with_two_availabilities_is_refused() {
        let all = ids(6);
        let train = mask(6, |i| i < 4);
        let col = |m: usize| [Bits::from_indices(m, [0, 1])];
        let p = spec("P", BlockKind::Presence, 2, None);
        let full = binary_block(p.clone(), &all, &col(6));
        let short = binary_block(p, &all[..5], &col(5));
        assemble(
            &[&full, &short].map(|raw| encode(raw, &all, &train)),
            &train,
        );
    }

    #[test]
    fn equal_priorities_break_ties_on_block_id() {
        // The synthetic noise block shares K's priority (contract §8).
        let all = ids(8);
        let train = mask(8, |i| i < 6);
        let col = [Bits::from_indices(8, [0, 1, 2])];
        let noise = binary_block(spec("noise", BlockKind::Anonymous, 5, None), &all, &col);
        let k = binary_block(spec("K", BlockKind::Anonymous, 5, Some(4000)), &all, &col);
        let arm = assemble(&[&noise, &k].map(|raw| encode(raw, &all, &train)), &train);
        let blocks: Vec<&str> = arm.prov[0]
            .sources
            .iter()
            .map(|s| s.block.as_str())
            .collect();
        assert_eq!(blocks, ["K", "noise"]);
    }

    #[test]
    fn a_sample_no_twin_observes_keeps_the_representatives_placeholder() {
        // Probe Q1: s0..s7 train. N's count has no row for s8 or s9, so its
        // AtLeast(5) bit takes the training median 5 there. K is N's twin on
        // training, observes s8 as 0 and cannot assess s9.
        let all = ids(10);
        let train = mask(10, |i| i < 8);
        let n = RawBlock {
            spec: spec("N", BlockKind::Count, 3, None),
            sample_ids: all[..8].to_vec(),
            names: vec!["copies".into()],
            cols: vec![RawColumn::Count {
                values: [5, 5, 5, 5, 5, 0, 0, 0].map(Some).to_vec(),
            }],
        };
        let mut k_assessable = Bits::ones(10);
        k_assessable.remove(9);
        let k = binary(
            spec("K", BlockKind::Anonymous, 5, Some(4000)),
            &all,
            Bits::from_indices(10, [0, 1, 2, 3, 4]),
            k_assessable,
        );
        let blocks = [&k, &n].map(|raw| encode(raw, &all, &train));
        assert_eq!(roles(&blocks[1]), [ColumnRole::AtLeast(5)]);
        assert!(blocks[1].x.get(0, 8) && blocks[1].x.get(0, 9));
        let arm = assemble(&blocks, &train);
        assert_eq!(arm.x.rows, 1);
        assert!(!arm.x.get(0, 8), "K observes s8 as 0");
        assert!(arm.x.get(0, 9), "nobody observes s9: N's median bit stands");
        assert_eq!(
            arm.x.row(0),
            Bits::from_indices(10, [0, 1, 2, 3, 4, 9]).words()
        );
        // s9 counts for neither source.
        assert_eq!(supplied(&arm, 0), [("N", 8), ("K", 1)]);
    }

    #[test]
    fn a_block_unavailable_indicator_observes_every_sample() {
        // Probe Q2a, contract §20: H has no row for s0, s1 or test s8, so its
        // block-unavailable indicator exists and is P's gene's twin on
        // training. P observes the gene absent on s8 and present on s9.
        let all = ids(10);
        let train = mask(10, |i| i < 8);
        let h = binary(
            spec("H", BlockKind::Hierarchical, 1, None),
            &[2, 3, 4, 5, 6, 7, 9].map(|i| all[i].clone()),
            Bits::from_indices(7, [0, 1]),
            Bits::ones(7),
        );
        let p = binary(
            spec("P", BlockKind::Presence, 2, None),
            &all,
            Bits::from_indices(10, [0, 1, 9]),
            Bits::ones(10),
        );
        let blocks = [&h, &p].map(|raw| encode(raw, &all, &train));
        assert_eq!(
            roles(&blocks[0]),
            [ColumnRole::Present, ColumnRole::BlockUnavailable]
        );
        let arm = assemble(&blocks, &train);
        assert_eq!(arm.x.rows, 2);
        assert_eq!(roles(&arm)[1], ColumnRole::BlockUnavailable);
        assert_eq!(supplied(&arm, 1), [("H", 10), ("P", 0)]);
        assert_eq!(arm.x.row(1), Bits::from_indices(10, [0, 1, 8]).words());
    }

    #[test]
    fn a_not_assessable_indicator_observes_where_its_block_has_a_row() {
        // Probe Q2b, contract §20: C cannot assess its column on s0, s1 or
        // test s8, so the not-assessable indicator exists and is P's gene's
        // twin on training. C has no row for test s9.
        let all = ids(10);
        let train = mask(10, |i| i < 8);
        let c = binary(
            spec("C", BlockKind::CuratedExact, 0, None),
            &all[..9],
            Bits::from_indices(9, [2, 3]),
            Bits::from_indices(9, [2, 3, 4, 5, 6, 7]),
        );
        let p = binary(
            spec("P", BlockKind::Presence, 2, None),
            &all,
            Bits::from_indices(10, [0, 1, 9]),
            Bits::ones(10),
        );
        let blocks = [&c, &p].map(|raw| encode(raw, &all, &train));
        assert_eq!(
            roles(&blocks[0]),
            [ColumnRole::Present, ColumnRole::NotAssessable]
        );
        let arm = assemble(&blocks, &train);
        assert_eq!(arm.x.rows, 2);
        assert_eq!(roles(&arm)[1], ColumnRole::NotAssessable);
        // C's 1 beats P's observed 0 on s8; P supplies s9, where C has no row.
        assert_eq!(supplied(&arm, 1), [("C", 9), ("P", 1)]);
        assert_eq!(arm.x.row(1), Bits::from_indices(10, [0, 1, 8, 9]).words());
    }

    #[test]
    fn a_within_block_twin_supplies_nothing_after_assemble() {
        // Probe Q3, contract §19: P's f0 and f1 agree on training, so
        // `block_candidates` keeps f0. C's twin has no row for test s8.
        let all = ids(10);
        let train = mask(10, |i| i < 8);
        let p = binary_block(
            spec("P", BlockKind::Presence, 2, None),
            &all,
            &[[0, 1, 2, 8], [0, 1, 2, 9]].map(|idx| Bits::from_indices(10, idx)),
        );
        let c = binary(
            spec("C", BlockKind::CuratedExact, 0, None),
            &all.iter()
                .filter(|s| *s != "s8")
                .cloned()
                .collect::<Vec<_>>(),
            Bits::from_indices(9, [0, 1, 2]),
            Bits::ones(9),
        );
        let p_cands = block_candidates(&encode(&p, &all, &train), &train, &CFG);
        assert_eq!(p_cands.x.rows, 1);
        let arm = assemble(&[encode(&c, &all, &train), p_cands], &train);
        assert_eq!(arm.x.rows, 1);
        let columns: Vec<(&str, Option<&str>)> = arm.prov[0]
            .sources
            .iter()
            .map(|s| (s.block.as_str(), s.column.as_deref()))
            .collect();
        assert_eq!(
            columns,
            [("C", Some("f0")), ("P", Some("f0")), ("P", Some("f1"))]
        );
        assert_eq!(supplied(&arm, 0), [("C", 9), ("P", 1), ("P", 0)]);
        assert_eq!(arm.x.row(0), Bits::from_indices(10, [0, 1, 2, 8]).words());
    }

    #[test]
    fn the_spec_hash_reads_hs_canonical_json_on_non_ascii_specs() {
        // H's canonical bytes escape é, DEL, ≥ and a character outside the
        // BMP as Python's json.dumps(..., ensure_ascii=True) does; serde_json
        // writes them raw, which is the form this module used before Q.
        let ascii = spec("K", BlockKind::Anonymous, 5, Some(4000));
        let wide = BlockSpec {
            version: "v\u{e9}\u{7f}\u{2265}\u{1f600}".into(),
            ..ascii.clone()
        };
        for s in [&ascii, &wide] {
            assert_eq!(spec_json(s).unwrap(), s.canonical_json().unwrap());
        }
        let json = spec_json(&wide).unwrap();
        assert!(
            json.contains(r#""version":"v\u00e9\u007f\u2265\ud83d\ude00""#),
            "{json}"
        );
        assert!(json.is_ascii());
        let raw = serde_json::to_string(&serde_json::to_value(&wide).unwrap()).unwrap();
        assert_ne!(json, raw);
        assert_eq!(
            serde_json::to_string(&serde_json::to_value(&ascii).unwrap()).unwrap(),
            spec_json(&ascii).unwrap()
        );
        let digest = ring::digest::digest(&ring::digest::SHA256, json.as_bytes());
        let hex: String = digest.as_ref().iter().map(|b| format!("{b:02x}")).collect();
        assert_eq!(spec_sha256(&wide).unwrap(), hex);
    }

    /// Main's close-check probe MP1 (C-1). Section 18's fallback where the
    /// representative's placeholder is 0 and a lower-priority twin's
    /// placeholder is 1. s0..s7 train. C (priority 0) has no row for s9 and
    /// cannot assess s8, so it reads 0 on both. N (priority 3) has no row for
    /// s8 or s9, so its AtLeast(5) bit takes the training median 5 there and
    /// reads 1. Nobody observes s8 or s9, so both take C's 0; a fallback that
    /// ORs the sources' rows (mutant MC03) gives 1.
    #[test]
    fn the_fallback_is_the_representatives_placeholder_where_it_is_0() {
        let all = ids(10);
        let train = mask(10, |i| i < 8);
        let n = RawBlock {
            spec: spec("N", BlockKind::Count, 3, None),
            sample_ids: all[..8].to_vec(),
            names: vec!["copies".into()],
            cols: vec![RawColumn::Count {
                values: [5, 5, 5, 5, 5, 0, 0, 0].map(Some).to_vec(),
            }],
        };
        let mut c_assessable = Bits::ones(9);
        c_assessable.remove(8);
        let c = binary(
            spec("C", BlockKind::CuratedExact, 0, None),
            &all[..9],
            Bits::from_indices(9, [0, 1, 2, 3, 4]),
            c_assessable,
        );
        let blocks = [&n, &c].map(|raw| encode(raw, &all, &train));
        assert_eq!(blocks[0].prov.len(), 1);
        assert_eq!(blocks[0].prov[0].sources[0].role, ColumnRole::AtLeast(5));
        assert!(blocks[0].x.get(0, 8) && blocks[0].x.get(0, 9));
        assert_eq!(blocks[1].prov.len(), 1, "C has no indicator");
        let arm = assemble(&blocks, &train);
        assert_eq!(arm.x.rows, 1);
        assert_eq!(supplied(&arm, 0), [("C", 8), ("N", 0)]);
        assert!(!arm.x.get(0, 8), "nobody observes s8: C's 0 stands");
        assert!(!arm.x.get(0, 9), "nobody observes s9: C's 0 stands");
        assert_eq!(
            arm.x.row(0),
            Bits::from_indices(10, [0, 1, 2, 3, 4]).words()
        );
    }
}
