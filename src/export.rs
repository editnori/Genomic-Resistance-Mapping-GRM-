//! Dataset export for a learner outside GRM: one dataset's projected matrix,
//! the fold plans `learn::eval::evaluate` scores it on, and each outer fold's
//! candidate and screened rows, so an outside model can be scored on the same
//! genomes, folds and screened columns as the five learners without
//! refitting them (docs/dataset-export.md).
//!
//! Nothing here fits a model. The export reads the same inputs `grm learn`
//! reads (`analysis::arm_inputs`), records the dataset's input fingerprint so
//! it can be matched to a published report, and writes:
//!
//! - `matrix.bin`: the projected matrix, row-major, each row `words` 64-bit
//!   little-endian words, bit `c % 64` of word `c / 64` set when feature row
//!   `r` is present in sample `c` (`bits::BitMatrix`);
//! - `dataset.json`: the samples, labels, lineages, duplicate groups, the
//!   plans of both schemes (`ablate::plan::fold_plan`, repeat 0 at the
//!   evaluation seed) and, per available plan and outer fold, the candidate
//!   rows count and the screened rows, as row indices of `matrix.bin` in
//!   screen order. The screen is the forest's: `learn::screen` over the
//!   outer fold's training candidates (`learn::training_candidates`) on its
//!   training samples.

use crate::ablate::plan::{FoldPlan, GroupIds, Scheme, fold_plan};
use crate::analysis;
use crate::bits::{BitMatrix, Bits};
use crate::layout::Layout;
use crate::learn::eval::EvalConfig;
use crate::learn::{screen, training_candidates};
use anyhow::{Context, Result, ensure};
use serde::Serialize;
use std::path::Path;

/// The `format` of every `dataset.json` this module writes.
pub const EXPORT_FORMAT: &str = "grm-dataset-export-v1";

/// One outer fold's columns as the forest sees them.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct FoldScreen {
    pub fold: usize,
    /// Training samples of this outer fold.
    pub n_train: usize,
    /// Rows left after `training_candidates` on this fold's training samples.
    pub n_candidates: usize,
    /// The screened rows, as indices into the projected matrix, best first.
    pub screened_rows: Vec<u32>,
}

/// One scheme's plan with its screens; `screens` is empty when the plan is
/// unavailable.
#[derive(Clone, Debug, Serialize)]
pub struct SplitExport {
    pub plan: FoldPlan,
    pub screens: Vec<FoldScreen>,
}

/// Everything `dataset.json` holds except the matrix bytes.
#[derive(Clone, Debug, Serialize)]
pub struct DatasetExport {
    pub format: String,
    pub id: String,
    pub species_slug: String,
    pub antibiotic: String,
    /// The fingerprint `grm learn` records as the report's
    /// `input_fingerprint` for the same inputs and settings.
    pub input_fingerprint: String,
    pub screened_features: usize,
    pub folds: usize,
    pub inner_folds: usize,
    pub seed: u64,
    pub min_support: u32,
    pub samples: Vec<String>,
    pub resistant: Vec<bool>,
    pub lineage: Vec<u32>,
    pub duplicates: Vec<Vec<usize>>,
    pub rows: usize,
    pub cols: usize,
    pub words: usize,
    /// Lowercase hex SHA-256 of `matrix.bin`.
    pub matrix_sha256: String,
    pub splits: Vec<SplitExport>,
}

/// The matrix as `matrix.bin` stores it.
pub fn matrix_bytes(x: &BitMatrix) -> Vec<u8> {
    let mut out = Vec::with_capacity(x.data.len() * 8);
    for w in &x.data {
        out.extend_from_slice(&w.to_le_bytes());
    }
    out
}

fn sha256_hex(bytes: &[u8]) -> String {
    let digest = ring::digest::digest(&ring::digest::SHA256, bytes);
    digest.as_ref().iter().map(|b| format!("{b:02x}")).collect()
}

/// The forest's screen of each outer fold of an available plan, as rows of
/// `x`. Only training samples are read: `training_candidates` and `screen`
/// both take the fold's training mask.
pub fn fold_screens(
    x: &BitMatrix,
    y: &Bits,
    plan: &FoldPlan,
    min_support: u32,
    screened_features: usize,
) -> Vec<FoldScreen> {
    if !plan.is_available() {
        return vec![];
    }
    (0..plan.k)
        .map(|f| {
            let train = Bits::from_indices(
                plan.outer.len(),
                (0..plan.outer.len()).filter(|&i| plan.outer[i] != f),
            );
            let (cx, rows) = training_candidates(x, &train, min_support);
            let screened_rows = screen(&cx, y, &train, screened_features)
                .into_iter()
                .map(|r| rows[r as usize])
                .collect();
            FoldScreen {
                fold: f,
                n_train: train.count() as usize,
                n_candidates: cx.rows,
                screened_rows,
            }
        })
        .collect()
}

/// Both schemes' plans, repeat 0 at the evaluation seed, with their screens.
pub fn dataset_splits(
    x: &BitMatrix,
    y: &Bits,
    lineage: &[u32],
    duplicates: &[Vec<usize>],
    min_support: u32,
    cfg: &EvalConfig,
) -> Result<Vec<SplitExport>> {
    ensure!(
        x.cols == y.len() && y.len() == lineage.len(),
        "export dimensions differ: {} columns, {} labels, {} lineages",
        x.cols,
        y.len(),
        lineage.len()
    );
    let groups = GroupIds::new(lineage.to_vec(), duplicates.to_vec(), None)?;
    [Scheme::Random, Scheme::Lineage]
        .into_iter()
        .map(|scheme| {
            let plan = fold_plan(y, &groups, scheme, cfg.folds, cfg.inner_folds, cfg.seed, 0)?;
            let screens = fold_screens(x, y, &plan, min_support, cfg.screened_features);
            Ok(SplitExport { plan, screens })
        })
        .collect()
}

/// Exports every dataset of `slug` (or the ones `only` names, by antibiotic
/// or id) under `out/<dataset id>/`. Returns the ids written.
pub fn export_species(
    layout: &Layout,
    slug: &str,
    only: &[String],
    out: &Path,
    cfg: &EvalConfig,
) -> Result<Vec<String>> {
    let (input, catalog) = analysis::load_inputs(layout, slug)?;
    let digest = analysis::species_digest(&input);
    let mut written = Vec::new();
    for spec in input
        .datasets
        .iter()
        .filter(|d| only.is_empty() || only.iter().any(|o| *o == d.antibiotic || *o == d.id))
    {
        let Some(e) = analysis::arm_inputs(&input, spec)? else {
            eprintln!("    {}: too few labeled genomes, not exported", spec.id);
            continue;
        };
        let fingerprint =
            analysis::dataset_fingerprint(&input, spec, catalog.as_ref(), cfg, digest);
        let splits = dataset_splits(&e.x, &e.y, &e.lineage, &e.duplicates, e.min_support, cfg)?;
        let bytes = matrix_bytes(&e.x);
        let export = DatasetExport {
            format: EXPORT_FORMAT.to_string(),
            id: spec.id.clone(),
            species_slug: spec.species_slug.clone(),
            antibiotic: spec.antibiotic.clone(),
            input_fingerprint: fingerprint,
            screened_features: cfg.screened_features,
            folds: cfg.folds,
            inner_folds: cfg.inner_folds,
            seed: cfg.seed,
            min_support: e.min_support,
            resistant: (0..e.y.len()).map(|i| e.y.contains(i)).collect(),
            samples: e.sample_ids,
            lineage: e.lineage,
            duplicates: e.duplicates,
            rows: e.x.rows,
            cols: e.x.cols,
            words: e.x.words,
            matrix_sha256: sha256_hex(&bytes),
            splits,
        };
        let dir = out.join(&spec.id);
        std::fs::create_dir_all(&dir).with_context(|| format!("create {}", dir.display()))?;
        std::fs::write(dir.join("matrix.bin"), &bytes)?;
        std::fs::write(dir.join("dataset.json"), serde_json::to_vec(&export)?)?;
        eprintln!(
            "    {}: exported {} rows x {} samples, fingerprint {}",
            spec.id, export.rows, export.cols, export.input_fingerprint
        );
        written.push(spec.id.clone());
    }
    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::learn::eval::evaluate;
    use rand::{Rng, SeedableRng};

    /// A synthetic dataset: `n` samples in `n / 4` lineages, row 0 present
    /// in no sample (dropped as a candidate), row 1 the label, the rest
    /// random. Labels alternate within lineages so both classes are in
    /// every fold.
    fn synthetic(n: usize, rows: usize, seed: u64) -> (BitMatrix, Bits, Vec<u32>) {
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(seed);
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
        let lineage: Vec<u32> = (0..n).map(|i| (i / 4) as u32).collect();
        let mut x = BitMatrix::with_cols(n);
        x.push_row(Bits::zeros(n).words());
        x.push_row(y.words());
        for _ in 2..rows {
            let row = Bits::from_indices(n, (0..n).filter(|_| rng.random_bool(0.4)));
            x.push_row(row.words());
        }
        (x, y, lineage)
    }

    #[test]
    fn matrix_bytes_keep_every_bit_in_row_and_column_order() {
        let (x, _, _) = synthetic(70, 9, 1);
        let bytes = matrix_bytes(&x);
        assert_eq!(bytes.len(), x.rows * x.words * 8);
        for r in 0..x.rows {
            for c in 0..x.cols {
                let w = r * x.words + c / 64;
                let word = u64::from_le_bytes(bytes[w * 8..w * 8 + 8].try_into().unwrap());
                assert_eq!((word >> (c % 64)) & 1 == 1, x.get(r, c), "row {r} col {c}");
            }
        }
    }

    #[test]
    fn plans_are_the_outer_folds_evaluate_scores() {
        let (raw, y, lineage) = synthetic(80, 30, 2);
        // `evaluate` takes a projection: filtered and collapsed on all samples.
        let (x, _) = training_candidates(&raw, &Bits::ones(80), 3);
        let cfg = EvalConfig {
            screened_features: 10,
            ..Default::default()
        };
        // Duplicate groups, each inside one lineage (`synthetic` puts four
        // samples in each): each must sit whole in one outer fold of the
        // random plan, as in `evaluate`'s.
        let duplicates = vec![vec![0, 1, 2], vec![5, 6], vec![40, 41, 42, 43]];
        let splits = dataset_splits(&x, &y, &lineage, &duplicates, 3, &cfg).unwrap();
        let out = evaluate("t", &x, &y, &lineage, &duplicates, 3, None, &cfg).unwrap();
        assert_eq!(splits.len(), out.splits.len());
        for (s, e) in splits.iter().zip(&out.splits) {
            assert_eq!(s.plan.scheme.as_str(), e.scheme);
            assert_eq!(s.plan.outer, e.test_fold, "{} outer folds", e.scheme);
            assert_eq!(s.plan.is_available(), e.unavailable_reason.is_none());
        }
        let random = &splits[0].plan;
        assert_eq!(random.scheme, Scheme::Random);
        for g in &duplicates {
            assert!(
                g.iter().all(|&i| random.outer[i] == random.outer[g[0]]),
                "duplicate group {g:?} split across outer folds"
            );
        }
    }

    #[test]
    fn screens_are_matrix_rows_ranked_on_training_samples_only() {
        let (x, y, lineage) = synthetic(80, 40, 3);
        let cfg = EvalConfig {
            screened_features: 6,
            ..Default::default()
        };
        let splits = dataset_splits(&x, &y, &lineage, &[], 3, &cfg).unwrap();
        for split in &splits {
            assert_eq!(split.screens.len(), cfg.folds);
            for s in &split.screens {
                let train =
                    Bits::from_indices(80, (0..80).filter(|&i| split.plan.outer[i] != s.fold));
                // Row 0 is no candidate, so a screen index that skipped the
                // mapping to matrix rows would name the wrong rows.
                assert!(!s.screened_rows.contains(&0));
                assert_eq!(s.screened_rows[0], 1, "the label row ranks first");
                let (cx, rows) = training_candidates(&x, &train, 3);
                let expect: Vec<u32> = screen(&cx, &y, &train, 6)
                    .into_iter()
                    .map(|r| rows[r as usize])
                    .collect();
                assert_eq!(s.screened_rows, expect);
                assert_eq!(s.n_train, train.count() as usize);
                assert_eq!(s.n_candidates, cx.rows);
                // Flipping every test label leaves the screen unchanged.
                let mut flipped = y.clone();
                for i in (0..80).filter(|&i| !train.contains(i)) {
                    if flipped.contains(i) {
                        flipped.remove(i);
                    } else {
                        flipped.insert(i);
                    }
                }
                let again = fold_screens(&x, &flipped, &split.plan, 3, 6);
                assert_eq!(again[s.fold].screened_rows, s.screened_rows);
            }
        }
    }

    #[test]
    fn an_unavailable_lineage_plan_exports_no_screens() {
        // Four lineages are fewer than twice the five folds.
        let (x, y, _) = synthetic(40, 12, 4);
        let lineage: Vec<u32> = (0..40).map(|i| (i % 4) as u32).collect();
        let cfg = EvalConfig::default();
        let splits = dataset_splits(&x, &y, &lineage, &[], 3, &cfg).unwrap();
        let lin = splits
            .iter()
            .find(|s| s.plan.scheme == Scheme::Lineage)
            .unwrap();
        assert!(!lin.plan.is_available());
        assert!(lin.screens.is_empty());
        let rnd = splits
            .iter()
            .find(|s| s.plan.scheme == Scheme::Random)
            .unwrap();
        assert_eq!(rnd.screens.len(), cfg.folds);
    }

    /// The chi-square of a row on the training samples, written out apart
    /// from `learn::screen`.
    fn chi_square(row: &Bits, y: &Bits, train: &Bits) -> f64 {
        let (mut a, mut b, mut c, mut d) = (0f64, 0f64, 0f64, 0f64);
        for i in (0..y.len()).filter(|&i| train.contains(i)) {
            match (row.contains(i), y.contains(i)) {
                (true, true) => a += 1.0,
                (true, false) => b += 1.0,
                (false, true) => c += 1.0,
                (false, false) => d += 1.0,
            }
        }
        let denom = (a + b) * (c + d) * (a + c) * (b + d);
        if denom == 0.0 {
            0.0
        } else {
            (a + b + c + d) * (a * d - b * c).powi(2) / denom
        }
    }

    /// A synthetic matrix of complementary row pairs: a row and its
    /// complement have the same chi-square on any training samples, so every
    /// score is tied at least twice, and an odd screen size cuts inside a
    /// tie. The complement comes first in some pairs and second in others.
    #[test]
    fn the_screen_keeps_the_lower_row_of_a_tie_at_the_cut() {
        let n = 60;
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(7);
        let y = Bits::from_indices(n, (0..n).filter(|i| i % 3 == 0));
        let lineage: Vec<u32> = (0..n).map(|i| (i / 2) as u32).collect();
        let mut x = BitMatrix::with_cols(n);
        let mut rows: Vec<Bits> = Vec::new();
        for p in 0..12 {
            // Rows that lean towards the label, by a different amount per pair.
            let lean = 0.05 + 0.03 * p as f64;
            let row = Bits::from_indices(
                n,
                (0..n).filter(|&i| {
                    rng.random_bool(if y.contains(i) {
                        0.5 + lean
                    } else {
                        0.5 - lean
                    })
                }),
            );
            let mut comp = Bits::ones(n);
            for i in 0..n {
                if row.contains(i) {
                    comp.remove(i);
                }
            }
            let (first, second) = if p % 2 == 0 { (comp, row) } else { (row, comp) };
            rows.push(first);
            rows.push(second);
        }
        for r in &rows {
            x.push_row(r.words());
        }
        let cfg = EvalConfig {
            screened_features: 5,
            ..Default::default()
        };
        let splits = dataset_splits(&x, &y, &lineage, &[], 3, &cfg).unwrap();
        let mut cuts_in_a_tie = 0;
        for split in splits.iter().filter(|s| s.plan.is_available()) {
            for s in &split.screens {
                let train =
                    Bits::from_indices(n, (0..n).filter(|&i| split.plan.outer[i] != s.fold));
                // Every row is a candidate here: distinct, and well inside the
                // support bounds on any fold's training samples.
                assert_eq!(s.n_candidates, x.rows);
                let mut order: Vec<(f64, u32)> = rows
                    .iter()
                    .enumerate()
                    .map(|(r, row)| (chi_square(row, &y, &train), r as u32))
                    .collect();
                order.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
                let want: Vec<u32> = order[..5].iter().map(|o| o.1).collect();
                assert_eq!(
                    s.screened_rows,
                    want,
                    "{} fold {}",
                    split.plan.scheme.as_str(),
                    s.fold
                );
                // The 5th and 6th rows tie: the cut falls inside the tie, and
                // the lower of the two rows is the one kept.
                let (kept, left) = (order[4], order[5]);
                assert_eq!(kept.0, left.0);
                assert!(kept.1 < left.1);
                cuts_in_a_tie += 1;
            }
        }
        assert!(cuts_in_a_tie >= 5, "{cuts_in_a_tie} folds checked");
    }
}
