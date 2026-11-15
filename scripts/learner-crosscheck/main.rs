//! The GRM side of the learner cross-check. `compare.py` beside this file
//! sets what this writes against the reference libraries; see README.md.
//!
//! usage: learner-crosscheck FIXTURES OUT_DIR
//!
//! It reads the saved test panels of FIXTURES
//! (tests/data/lasso-stopping-fixtures.txt), fits every learner but the L1
//! logistic regression on each with the public functions of `grm::learn`, and
//! writes OUT_DIR/<panel>.json and OUT_DIR/index.json. It changes nothing
//! under `src/`, so the learners are those of the commit it is built from.
//!
//! Every learner is fitted on the samples of the panel's `fit` mask and
//! scored on all samples. Settings are the defaults of `EvalConfig`,
//! `ForestParams` and `GbtParams`.
//!
//! No reference library can be given GRM's random stream, so this also
//! replays the draws GRM made, from the documented seeds, and writes them:
//! each forest tree's subsample and each node's drawn features, and each
//! boosting round's drawn features. It stops if a replayed draw does not
//! reproduce the fitted tree.
//!
//! The output holds no time and no path, and the same build writes the same
//! bytes at any thread count.

use anyhow::{Context, Result, bail, ensure};
use grm::bits::{BitMatrix, Bits};
use grm::learn::eval::EvalConfig;
use grm::learn::gbt::{self, GTree, Gbt, GbtParams, Objective};
use grm::learn::scm::{self, Kind};
use grm::learn::tree::{self, ForestParams, Tree, TreeParams};
use grm::learn::{sample_indices, screen, training_folds, tuning_holdout, validate_folds};
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use serde_json::{Value, json};
use std::path::Path;

const USAGE: &str = "usage: learner-crosscheck FIXTURES OUT_DIR";

/// `eval::GBT_HOLDOUT`, which is private to the crate.
const GBT_HOLDOUT: f64 = 0.15;

struct Panel {
    name: String,
    n: usize,
    declared_features: usize,
    y: Option<Bits>,
    fit: Option<Bits>,
    x: BitMatrix,
}

/// A bitset of the fixture file: hex nibbles, low bit first.
fn hex_bits(n: usize, s: &str) -> Result<Bits> {
    ensure!(
        s.len() == n.div_ceil(4),
        "a bitset of {n} samples has {} nibbles",
        s.len()
    );
    let mut bits = Bits::zeros(n);
    for (i, c) in s.chars().enumerate() {
        let v = c
            .to_digit(16)
            .with_context(|| format!("{c:?} is no hex digit"))?;
        for b in 0..4 {
            if (v >> b) & 1 == 1 {
                ensure!(4 * i + b < n, "a bit is set past sample {n}");
                bits.insert(4 * i + b);
            }
        }
    }
    Ok(bits)
}

fn bits_hex(bits: &Bits) -> String {
    (0..bits.len().div_ceil(4))
        .map(|i| {
            let v = (0..4)
                .filter(|b| 4 * i + b < bits.len() && bits.contains(4 * i + b))
                .fold(0u32, |v, b| v | (1 << b));
            char::from_digit(v, 16).expect("a nibble")
        })
        .collect()
}

fn read_panels(path: &Path) -> Result<Vec<Panel>> {
    let text =
        std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut panels: Vec<Panel> = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let t: Vec<&str> = line.split(' ').collect();
        let at = || format!("{} line {}", path.display(), i + 1);
        match t[0] {
            "fixture" => {
                ensure!(t.len() == 6, "{}: a fixture line has 6 fields", at());
                let n: usize = t[2].parse().with_context(at)?;
                panels.push(Panel {
                    name: t[1].to_string(),
                    n,
                    declared_features: t[3].parse().with_context(at)?,
                    y: None,
                    fit: None,
                    x: BitMatrix::with_cols(n),
                });
            }
            "y" | "fit" | "r" => {
                ensure!(t.len() == 2, "{}: a bitset line has 2 fields", at());
                let p = panels
                    .last_mut()
                    .with_context(|| format!("{}: no fixture yet", at()))?;
                let bits = hex_bits(p.n, t[1]).with_context(at)?;
                match t[0] {
                    "y" => p.y = Some(bits),
                    "fit" => p.fit = Some(bits),
                    _ => p.x.push_row(bits.words()),
                }
            }
            // The saved lasso paths, which this check does not use.
            "point" => {}
            _ if line.is_empty() || line.starts_with('#') => {}
            other => bail!("{}: unknown line type {other:?}", at()),
        }
    }
    for p in &panels {
        ensure!(
            p.y.is_some() && p.fit.is_some(),
            "panel {} lacks y or fit",
            p.name
        );
        ensure!(
            p.x.rows == p.declared_features,
            "panel {} declares {} features and has {}",
            p.name,
            p.declared_features,
            p.x.rows
        );
    }
    Ok(panels)
}

fn tree_json(t: &Tree) -> Value {
    let split = |f: &dyn Fn(&tree::Split) -> Value| -> Vec<Value> {
        t.nodes
            .iter()
            .map(|n| n.split.as_ref().map_or(Value::Null, f))
            .collect()
    };
    json!({
        "value": t.nodes.iter().map(|n| n.value).collect::<Vec<_>>(),
        "n": t.nodes.iter().map(|n| n.n).collect::<Vec<_>>(),
        "n_resistant": t.nodes.iter().map(|n| n.n_resistant).collect::<Vec<_>>(),
        "depth": t.nodes.iter().map(|n| n.depth).collect::<Vec<_>>(),
        "feature": split(&|s| s.feature.into()),
        "gain": split(&|s| s.gain.into()),
        "absent": split(&|s| s.absent.into()),
        "present": split(&|s| s.present.into()),
    })
}

fn gtree_json(t: &GTree) -> Value {
    let split = |f: &dyn Fn(&gbt::GNode) -> Value| -> Vec<Value> {
        t.nodes
            .iter()
            .map(|n| {
                if n.feature.is_some() {
                    f(n)
                } else {
                    Value::Null
                }
            })
            .collect()
    };
    json!({
        "feature": split(&|n| n.feature.into()),
        "absent": split(&|n| n.absent.into()),
        "present": split(&|n| n.present.into()),
        "gain": split(&|n| n.gain.into()),
        "value": t.nodes.iter().map(|n| n.value).collect::<Vec<_>>(),
    })
}

fn scm_fits(p: &Panel, y: &Bits, train: &Bits, cfg: &EvalConfig) -> Vec<Value> {
    let mut fits = Vec::new();
    for kind in [Kind::Conjunction, Kind::Disjunction] {
        for &penalty in &scm::P_GRID {
            let m = scm::fit(&p.x, y, train, kind, penalty, cfg.scm_max_rules, true);
            fits.push(json!({
                "kind": kind,
                "p": penalty,
                "rules": m.rules.iter().map(|r| json!([r.feature, r.presence])).collect::<Vec<_>>(),
                "steps": m.steps.iter().map(|s| json!({
                    "covered": s.covered,
                    "lost": s.lost,
                    "utility": s.utility,
                    "n_equivalent": s.n_equivalent,
                })).collect::<Vec<_>>(),
                // The calls of the model cut to its first 0, 1, ... rules.
                "calls": (0..=m.rules.len())
                    .map(|k| bits_hex(&m.truncated(k).predict(&p.x)))
                    .collect::<Vec<_>>(),
            }));
        }
    }
    fits
}

/// What production would select on this panel, recorded and not compared.
fn selected(p: &Panel, y: &Bits, train: &Bits, cfg: &EvalConfig) -> Value {
    let k = cfg.inner_folds;
    let folds = training_folds(y, train, None, k, cfg.seed);
    if let Err(e) = validate_folds(y, train, &folds, k) {
        return json!({ "error": format!("{e:#}") });
    }
    let c = scm::select(&p.x, y, train, cfg.scm_max_rules, &folds, k);
    let (depth, tree_error) = tree::select_depth(&p.x, y, train, &folds, k, cfg.cart_max_depth);
    json!({
        "inner_folds": k,
        "seed": cfg.seed,
        "scm": { "kind": c.kind, "p": c.p, "n_rules": c.n_rules, "cv_error": c.cv_error },
        "tree": { "depth": depth, "cv_error": tree_error },
    })
}

fn decision_tree(p: &Panel, y: &Bits, train: &Bits, cfg: &EvalConfig) -> Value {
    let t = tree::fit_tree(
        &p.x,
        y,
        train,
        &TreeParams {
            max_depth: cfg.cart_max_depth,
            min_leaf: 1,
            max_features: None,
            seed: 0,
        },
    );
    json!({
        "max_depth": cfg.cart_max_depth,
        "tree": tree_json(&t),
        "scores_by_depth": (1..=cfg.cart_max_depth).map(|d| t.predict(&p.x, d)).collect::<Vec<_>>(),
    })
}

/// The feature draws of one forest tree, replayed from its seed in the order
/// `fit_tree` visits nodes: a stack, the present child before the absent one,
/// and one draw for every node that is searched for a split.
fn node_draws(
    t: &Tree,
    n_features: usize,
    m: usize,
    params: &ForestParams,
    seed: u64,
) -> Result<Vec<Value>> {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut draws = Vec::new();
    let mut stack = vec![0usize];
    while let Some(id) = stack.pop() {
        let node = &t.nodes[id];
        let (n, pos) = (node.n, node.n_resistant);
        if node.depth >= params.max_depth || pos == 0 || pos == n || n < 2 * params.min_leaf {
            continue;
        }
        let draw = sample_indices(n_features, m, &mut rng);
        if let Some(s) = &node.split {
            ensure!(
                draw.contains(&s.feature),
                "node {id} splits on a feature outside its replayed draw"
            );
            stack.push(s.absent as usize);
            stack.push(s.present as usize);
        }
        draws.push(json!([id, draw]));
    }
    Ok(draws)
}

fn forest(p: &Panel, y: &Bits, train: &Bits, cfg: &EvalConfig) -> Result<Value> {
    let params = cfg.forest;
    // As `eval::fit_rf`: the forest sees the features in the screen's order
    // and its trees name positions in that order.
    let order = screen(&p.x, y, train, cfg.screened_features);
    let mut screened = BitMatrix::with_cols(p.x.cols);
    for &f in &order {
        screened.push_row(p.x.row(f as usize));
    }
    let fitted = tree::fit_forest(&screened, y, train, &params);

    let members: Vec<usize> = train.iter_ones().collect();
    let take = ((members.len() as f64) * params.subsample).round() as usize;
    let m = ((screened.rows as f64).sqrt().ceil() as usize).max(1);
    let mut trees = Vec::new();
    for (t, grown) in fitted.trees.iter().enumerate() {
        let seed = params.seed.wrapping_mul(1_000_003).wrapping_add(t as u64);
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        rng.set_stream(1);
        let chosen = sample_indices(members.len(), take, &mut rng);
        let mask = Bits::from_indices(train.len(), chosen.into_iter().map(|i| members[i as usize]));
        // The replayed subsample is the one the forest used only if one tree
        // fitted on it is the forest's tree.
        let again = tree::fit_tree(
            &screened,
            y,
            &mask,
            &TreeParams {
                max_depth: params.max_depth,
                min_leaf: params.min_leaf,
                max_features: Some(m),
                seed,
            },
        );
        ensure!(
            tree_json(&again) == tree_json(grown),
            "panel {} forest tree {t}: the replayed subsample does not give the forest's tree",
            p.name
        );
        let draws = node_draws(grown, screened.rows, m, &params, seed)
            .with_context(|| format!("panel {} forest tree {t}", p.name))?;
        trees.push(json!({
            "seed": seed,
            "subsample": bits_hex(&mask),
            "tree": tree_json(grown),
            "draws": draws,
        }));
    }
    Ok(json!({
        "params": params,
        "screen": order,
        "take": take,
        "max_features": m,
        "trees": trees,
        "scores": fitted.predict(&screened),
    }))
}

/// The features each boosting round drew, replayed from the seed as
/// `gbt::fit` draws them, for the first `rounds` rounds.
fn round_draws(features: &[u32], params: &GbtParams, rounds: usize) -> Vec<Vec<u32>> {
    let per_tree = ((features.len() as f64 * params.colsample).ceil() as usize).max(1);
    let mut rng = ChaCha8Rng::seed_from_u64(params.seed);
    (0..rounds)
        .map(|_| {
            sample_indices(features.len(), per_tree, &mut rng)
                .into_iter()
                .map(|i| features[i as usize])
                .collect()
        })
        .collect()
}

fn gbt_json(p: &Panel, m: &Gbt, features: &[u32], params: &GbtParams) -> Result<Value> {
    let draws = round_draws(features, params, m.trees.len());
    for (r, (t, draw)) in m.trees.iter().zip(&draws).enumerate() {
        ensure!(
            t.nodes
                .iter()
                .filter_map(|n| n.feature)
                .all(|f| draw.contains(&f)),
            "panel {} boosting round {r} splits on a feature outside its replayed draw",
            p.name
        );
    }
    Ok(json!({
        "params": params,
        "features": features,
        "kept": m.trees.len(),
        "draws": draws,
        "trees": m.trees.iter().map(gtree_json).collect::<Vec<_>>(),
        "margin": m.margin(&p.x),
        "probability": m.predict_proba(&p.x),
    }))
}

/// Boosting on the whole `fit` mask with the same mask as validation and
/// patience as long as the run, so every round is grown and early stopping
/// drops only the trailing rounds that no longer lower the training loss.
fn gbt_tree_by_tree(p: &Panel, y: &Bits, train: &Bits, cfg: &EvalConfig) -> Result<Value> {
    let features = screen(&p.x, y, train, cfg.screened_features);
    let params = GbtParams {
        patience: cfg.gbt.rounds,
        ..cfg.gbt
    };
    let m = gbt::fit(
        &p.x,
        &Objective::Logistic { y },
        train,
        train,
        &features,
        &params,
    );
    gbt_json(p, &m, &features, &params)
}

/// Boosting as `eval::fit_gbt` runs it: a tuning holdout, the screen on the
/// fitting part, early stopping on the holdout. Once with every feature
/// offered to every tree, which a library's own loop can follow, and once at
/// the production column share for its tree count.
fn gbt_end_to_end(p: &Panel, y: &Bits, train: &Bits, cfg: &EvalConfig) -> Result<Value> {
    let (fit_m, val) = match tuning_holdout(train, y, None, "group", GBT_HOLDOUT, cfg.seed) {
        Ok(parts) => parts,
        Err(e) => return Ok(json!({ "error": format!("{e:#}") })),
    };
    let features = screen(&p.x, y, &fit_m, cfg.screened_features);
    let objective = Objective::Logistic { y };
    let all = GbtParams {
        seed: cfg.seed,
        colsample: 1.0,
        ..cfg.gbt
    };
    let m_all = gbt::fit(&p.x, &objective, &fit_m, &val, &features, &all);
    let half = GbtParams {
        seed: cfg.seed,
        ..cfg.gbt
    };
    let m_half = gbt::fit(&p.x, &objective, &fit_m, &val, &features, &half);
    Ok(json!({
        "holdout": GBT_HOLDOUT,
        "seed": cfg.seed,
        "fit": bits_hex(&fit_m),
        "validation": bits_hex(&val),
        "colsample_1": gbt_json(p, &m_all, &features, &all)?,
        "colsample_production": {
            "params": half,
            "kept": m_half.trees.len(),
            "probability": m_half.predict_proba(&p.x),
        },
    }))
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [fixtures, out] = &args[..] else {
        bail!(USAGE);
    };
    let out = Path::new(out);
    std::fs::create_dir_all(out)?;
    let cfg = EvalConfig::default();
    let panels = read_panels(Path::new(fixtures))?;
    let mut index = Vec::new();
    for p in &panels {
        let (y, train) = (
            p.y.as_ref().expect("checked"),
            p.fit.as_ref().expect("checked"),
        );
        let counts = json!({
            "samples": p.n,
            "features": p.x.rows,
            "resistant": y.count(),
            "fit": train.count(),
            "fit_resistant": train.and(y).count(),
            "present": (0..p.x.rows).map(|f| grm::bits::count(p.x.row(f))).collect::<Vec<_>>(),
        });
        let doc = json!({
            "panel": p.name,
            "counts": counts,
            "scm": { "max_rules": cfg.scm_max_rules, "fits": scm_fits(p, y, train, &cfg) },
            "selected": selected(p, y, train, &cfg),
            "decision_tree": decision_tree(p, y, train, &cfg),
            "forest": forest(p, y, train, &cfg)?,
            "gbt_tree_by_tree": gbt_tree_by_tree(p, y, train, &cfg)?,
            "gbt_end_to_end": gbt_end_to_end(p, y, train, &cfg)?,
        });
        std::fs::write(out.join(format!("{}.json", p.name)), format!("{doc}\n"))?;
        eprintln!("{}: {} samples, {} features", p.name, p.n, p.x.rows);
        index.push(json!({ "panel": p.name, "counts": doc["counts"] }));
    }
    let index =
        json!({ "evaluation_revision": grm::analysis::EVALUATION_REVISION, "panels": index });
    std::fs::write(out.join("index.json"), format!("{index}\n"))?;
    Ok(())
}
