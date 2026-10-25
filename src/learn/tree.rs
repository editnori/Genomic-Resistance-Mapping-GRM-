//! CART decision trees on presence/absence features, and random forests.
//!
//! A node is a sample mask; a split on feature f sends the samples carrying
//! f to `present` and the rest to `absent`. With class-balanced weights the
//! impurity is weighted Gini, and each node's value is the weighted share of
//! resistant samples.
//!
//! A node splits only on positive gain (above 1e-12), never on zero gain. An
//! exact interaction such as resistance = f0 XOR f1, with each combination
//! on equally many samples, leaves both sides of either feature half
//! resistant, so no first split has gain and the tree stays one 0.5 leaf. A
//! forest node that draws `max_features` features becomes a leaf when none
//! of them can split; it does not draw again.
//!
//! Every node keeps its own value, so predicting with a smaller `max_depth`
//! than the tree was grown to gives exactly the tree a depth-limited fit
//! would have produced: greedy splits above depth d do not depend on the
//! limit. Cross-validation scores all depths from one fit.

use super::{fold_masks, sample_indices};
use crate::bits::{self, BitMatrix, Bits};
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Split {
    pub feature: u32,
    pub gain: f64,
    pub absent: u32,
    pub present: u32,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Node {
    pub value: f64,
    pub n: u32,
    pub n_resistant: u32,
    pub depth: u32,
    pub split: Option<Split>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Tree {
    pub nodes: Vec<Node>,
}

#[derive(Clone, Copy, Debug)]
pub struct TreeParams {
    pub max_depth: u32,
    pub min_leaf: u32,
    /// Features drawn at random per node; None evaluates every feature.
    pub max_features: Option<usize>,
    pub seed: u64,
}

struct Weights {
    pos: f64,
    neg: f64,
}

fn gini(p: f64, n: f64) -> f64 {
    let t = p + n;
    if t == 0.0 {
        return 0.0;
    }
    let q = p / t;
    2.0 * q * (1.0 - q)
}

pub fn fit_tree(x: &BitMatrix, y: &Bits, train: &Bits, params: &TreeParams) -> Tree {
    let n_pos = train.and(y).count() as f64;
    let n_neg = train.count() as f64 - n_pos;
    let total = n_pos + n_neg;
    let w = Weights {
        pos: total / (2.0 * n_pos.max(1.0)),
        neg: total / (2.0 * n_neg.max(1.0)),
    };
    let mut rng = ChaCha8Rng::seed_from_u64(params.seed);
    let mut tree = Tree { nodes: Vec::new() };
    let mut stack: Vec<(usize, Bits)> = Vec::new();

    let push = |tree: &mut Tree, mask: &Bits, depth: u32| -> usize {
        let pos = mask.and(y).count();
        let neg = mask.count() - pos;
        // The class-weighted share from exact count products: rounded class
        // weights can score a node with the training class ratio below 0.5.
        let (wp, wn) = (pos as f64 * n_neg.max(1.0), neg as f64 * n_pos.max(1.0));
        tree.nodes.push(Node {
            value: if wp + wn == 0.0 { 0.5 } else { wp / (wp + wn) },
            n: pos + neg,
            n_resistant: pos,
            depth,
            split: None,
        });
        tree.nodes.len() - 1
    };
    let root = push(&mut tree, train, 0);
    stack.push((root, train.clone()));

    while let Some((id, mask)) = stack.pop() {
        let node = &tree.nodes[id];
        let (n, pos) = (node.n, node.n_resistant);
        if node.depth >= params.max_depth || pos == 0 || pos == n || n < 2 * params.min_leaf {
            continue;
        }
        let mask_pos = mask.and(y);
        let mask_neg = mask.and_not(y);
        let (np, nn) = (pos as f64 * w.pos, (n - pos) as f64 * w.neg);
        let parent = (np + nn) * gini(np, nn);

        let evaluate = |f: u32| -> Option<(f64, u32)> {
            let row = x.row(f as usize);
            let pp = bits::and_count(row, mask_pos.words());
            let pn = bits::and_count(row, mask_neg.words());
            let present = pp + pn;
            if present < params.min_leaf || n - present < params.min_leaf {
                return None;
            }
            let (lp, ln) = (pp as f64 * w.pos, pn as f64 * w.neg);
            let (rp, rn) = (np - lp, nn - ln);
            let gain = parent - (lp + ln) * gini(lp, ln) - (rp + rn) * gini(rp, rn);
            (gain > 1e-12).then_some((gain, f))
        };
        let better = |a: (f64, u32), b: (f64, u32)| {
            if a.0 > b.0 || (a.0 == b.0 && a.1 < b.1) {
                a
            } else {
                b
            }
        };
        let best = match params.max_features {
            Some(m) => sample_indices(x.rows, m, &mut rng)
                .into_iter()
                .filter_map(evaluate)
                .reduce(better),
            None => (0..x.rows as u32)
                .into_par_iter()
                .filter_map(evaluate)
                .reduce_with(better),
        };
        let Some((gain, feature)) = best else {
            continue;
        };
        let row = x.row(feature as usize);
        let present_mask = mask.and_slice(row);
        let absent_mask = mask.and_not_slice(row);
        let depth = tree.nodes[id].depth + 1;
        let absent = push(&mut tree, &absent_mask, depth);
        let present = push(&mut tree, &present_mask, depth);
        tree.nodes[id].split = Some(Split {
            feature,
            gain: gain / total,
            absent: absent as u32,
            present: present as u32,
        });
        stack.push((absent, absent_mask));
        stack.push((present, present_mask));
    }
    tree
}

impl Tree {
    /// Resistant-class scores for every sample, treating nodes at
    /// `max_depth` as leaves.
    pub fn predict(&self, x: &BitMatrix, max_depth: u32) -> Vec<f64> {
        let mut out = vec![0.0; x.cols];
        let mut stack = vec![(0usize, Bits::ones(x.cols))];
        while let Some((id, mask)) = stack.pop() {
            let node = &self.nodes[id];
            match &node.split {
                Some(s) if node.depth < max_depth => {
                    let row = x.row(s.feature as usize);
                    stack.push((s.absent as usize, mask.and_not_slice(row)));
                    stack.push((s.present as usize, mask.and_slice(row)));
                }
                _ => {
                    for i in mask.iter_ones() {
                        out[i] = node.value;
                    }
                }
            }
        }
        out
    }

    pub fn depth(&self) -> u32 {
        self.nodes.iter().map(|n| n.depth).max().unwrap_or(0)
    }

    /// Copy of the tree with everything below `max_depth` removed.
    pub fn truncated(&self, max_depth: u32) -> Tree {
        let mut nodes = Vec::new();
        let mut stack = vec![(0usize, usize::MAX, false)];
        while let Some((old, parent, present)) = stack.pop() {
            let mut node = self.nodes[old].clone();
            let split = node.split.take().filter(|_| node.depth < max_depth);
            let new = nodes.len();
            nodes.push(node);
            if parent != usize::MAX {
                let s: &mut Split = nodes[parent].split.as_mut().unwrap();
                if present {
                    s.present = new as u32
                } else {
                    s.absent = new as u32
                }
            }
            if let Some(s) = split {
                stack.push((s.absent as usize, new, false));
                stack.push((s.present as usize, new, true));
                nodes[new].split = Some(s);
            }
        }
        Tree { nodes }
    }

    pub fn add_importance(&self, into: &mut std::collections::HashMap<u32, f64>, scale: f64) {
        for node in &self.nodes {
            if let Some(s) = &node.split {
                *into.entry(s.feature).or_default() += s.gain * scale;
            }
        }
    }
}

/// Picks the depth with the lowest mean balanced error across inner folds.
pub fn select_depth(
    x: &BitMatrix,
    y: &Bits,
    train: &Bits,
    inner_folds: &[usize],
    k: usize,
    max_depth: u32,
) -> (u32, f64) {
    let mut err = vec![0.0; max_depth as usize];
    let masks = fold_masks(inner_folds, k);
    for (fit_mask, val) in &masks {
        let (fit_mask, val) = (fit_mask.and(train), val.and(train));
        let tree = fit_tree(
            x,
            y,
            &fit_mask,
            &TreeParams {
                max_depth,
                min_leaf: 1,
                max_features: None,
                seed: 0,
            },
        );
        for d in 1..=max_depth {
            let scores = tree.predict(x, d);
            err[d as usize - 1] += balanced_error(&scores, y, &val) / masks.len() as f64;
        }
    }
    let (best, e) = err
        .iter()
        .enumerate()
        .fold((0, f64::INFINITY), |acc, (i, &e)| {
            if e < acc.1 - 1e-12 { (i, e) } else { acc }
        });
    (best as u32 + 1, e)
}

pub fn balanced_error(scores: &[f64], y: &Bits, mask: &Bits) -> f64 {
    let (mut fnr, mut fpr, mut np, mut nn) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for i in mask.iter_ones() {
        if y.contains(i) {
            np += 1.0;
            fnr += (scores[i] < 0.5) as u8 as f64;
        } else {
            nn += 1.0;
            fpr += (scores[i] >= 0.5) as u8 as f64;
        }
    }
    (fnr / np.max(1.0) + fpr / nn.max(1.0)) / 2.0
}

#[derive(Clone, Copy, Debug, Serialize)]
pub struct ForestParams {
    pub n_trees: usize,
    pub subsample: f64,
    pub max_depth: u32,
    pub min_leaf: u32,
    /// Evaluation replaces this per fit with a seed from `EvalConfig::seed`
    /// and the outer fold, so a recorded config leaves it out.
    #[serde(skip)]
    pub seed: u64,
}

impl Default for ForestParams {
    fn default() -> Self {
        ForestParams {
            n_trees: 200,
            subsample: 0.632,
            max_depth: 30,
            min_leaf: 1,
            seed: 17,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Forest {
    pub trees: Vec<Tree>,
}

pub fn fit_forest(x: &BitMatrix, y: &Bits, train: &Bits, params: &ForestParams) -> Forest {
    let members: Vec<usize> = train.iter_ones().collect();
    let take = ((members.len() as f64) * params.subsample).round() as usize;
    let max_features = ((x.rows as f64).sqrt().ceil() as usize).max(1);
    let trees = (0..params.n_trees)
        .into_par_iter()
        .map(|t| {
            let seed = params.seed.wrapping_mul(1_000_003).wrapping_add(t as u64);
            // fit_tree draws node features from stream 0 of this seed. The
            // subsample reads stream 1, so the two draws share no numbers.
            let mut rng = ChaCha8Rng::seed_from_u64(seed);
            rng.set_stream(1);
            let chosen = sample_indices(members.len(), take, &mut rng);
            let mask =
                Bits::from_indices(train.len(), chosen.into_iter().map(|i| members[i as usize]));
            fit_tree(
                x,
                y,
                &mask,
                &TreeParams {
                    max_depth: params.max_depth,
                    min_leaf: params.min_leaf,
                    max_features: Some(max_features),
                    seed,
                },
            )
        })
        .collect();
    Forest { trees }
}

impl Forest {
    pub fn predict(&self, x: &BitMatrix) -> Vec<f64> {
        let per_tree: Vec<Vec<f64>> = self
            .trees
            .par_iter()
            .map(|t| t.predict(x, u32::MAX))
            .collect();
        let mut out = vec![0.0; x.cols];
        for scores in &per_tree {
            for (o, s) in out.iter_mut().zip(scores) {
                *o += s;
            }
        }
        // Divide once: summing pre-divided votes can put an exact half
        // below 0.5.
        let n = per_tree.len().max(1) as f64;
        out.iter_mut().for_each(|o| *o /= n);
        out
    }

    /// Mean decrease in weighted Gini impurity per feature.
    pub fn importance(&self) -> std::collections::HashMap<u32, f64> {
        let mut imp = std::collections::HashMap::new();
        for t in &self.trees {
            t.add_importance(&mut imp, 1.0 / self.trees.len() as f64);
        }
        imp
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(rows: &[&[usize]], cols: usize) -> BitMatrix {
        let mut x = BitMatrix::with_cols(cols);
        for r in rows {
            x.push_row(Bits::from_indices(cols, r.iter().copied()).words());
        }
        x
    }

    #[test]
    fn learns_xor_with_depth_two_and_truncates() {
        // Resistant iff exactly one of f0, f1 is present. f2, on susceptible
        // sample 3 alone, is the only positive-gain first split; once it
        // unbalances the classes, f0 and f1 split with positive gain.
        let x = matrix(&[&[0, 1, 4, 5], &[0, 2, 4, 6], &[3]], 8);
        let y = Bits::from_indices(8, [1, 2, 5, 6]);
        let train = Bits::ones(8);
        let tree = fit_tree(
            &x,
            &y,
            &train,
            &TreeParams {
                max_depth: 3,
                min_leaf: 1,
                max_features: None,
                seed: 0,
            },
        );
        let scores = tree.predict(&x, 3);
        for (i, &score) in scores[..8].iter().enumerate() {
            assert_eq!(score >= 0.5, y.contains(i), "sample {i}");
        }
        // XOR has no single informative split, so a depth-1 tree cannot fit
        // it; truncation must reproduce the depth-limited prediction.
        let shallow = tree.truncated(1);
        assert!(shallow.depth() <= 1);
        assert_eq!(shallow.predict(&x, u32::MAX), tree.predict(&x, 1));
    }

    #[test]
    fn pure_xor_has_no_positive_gain_split_so_the_tree_is_one_leaf() {
        // Resistant iff exactly one of f0, f1 is present, each combination
        // on 20 of 80 samples. Either feature leaves both of its sides half
        // resistant, so its gain is zero and the root stays the only node.
        let f0: Vec<usize> = (40..80).collect();
        let f1: Vec<usize> = (0..80).filter(|i| i / 20 % 2 == 1).collect();
        let x = matrix(&[&f0, &f1], 80);
        let y = Bits::from_indices(80, (0..80).filter(|i| (*i >= 40) != (i / 20 % 2 == 1)));
        let tree = fit_tree(
            &x,
            &y,
            &Bits::ones(80),
            &TreeParams {
                max_depth: 8,
                min_leaf: 1,
                max_features: None,
                seed: 0,
            },
        );
        assert_eq!(tree.nodes.len(), 1);
        assert_eq!(tree.predict(&x, 8), vec![0.5; 80]);
    }

    #[test]
    fn forest_separates_a_clean_marker() {
        let x = matrix(&[&[0, 1, 2, 3, 4], &[0, 5, 7], &[1, 6]], 10);
        let y = Bits::from_indices(10, [0, 1, 2, 3, 4]);
        let train = Bits::ones(10);
        let forest = fit_forest(
            &x,
            &y,
            &train,
            &ForestParams {
                n_trees: 50,
                ..Default::default()
            },
        );
        let scores = forest.predict(&x);
        let auc = super::super::auroc(&(0..10).map(|i| y.contains(i)).collect::<Vec<_>>(), &scores)
            .unwrap();
        assert!(auc > 0.95, "auc {auc}");
        let imp = forest.importance();
        let top = imp.iter().max_by(|a, b| a.1.total_cmp(b.1)).unwrap().0;
        assert_eq!(*top, 0);
    }

    #[test]
    fn forest_nodes_stay_leaves_when_no_drawn_feature_can_split() {
        // f0 separates the classes; f1 to f3 are on no sample and cannot
        // split. Every tree holds all 8 samples and each node draws
        // ceil(sqrt(4)) = 2 of the 4 features, so some roots draw only empty
        // ones. Those roots stay impure 0.5 leaves: no node draws again.
        let x = matrix(&[&[0, 1, 2, 3], &[], &[], &[]], 8);
        let y = Bits::from_indices(8, [0, 1, 2, 3]);
        let forest = fit_forest(
            &x,
            &y,
            &Bits::ones(8),
            &ForestParams {
                n_trees: 40,
                subsample: 1.0,
                ..Default::default()
            },
        );
        let (stopped, split): (Vec<&Tree>, Vec<&Tree>) = forest
            .trees
            .iter()
            .partition(|t| t.nodes[0].split.is_none());
        assert!(!stopped.is_empty() && !split.is_empty());
        for t in stopped {
            let root = &t.nodes[0];
            assert_eq!(
                (t.nodes.len(), root.n, root.n_resistant, root.value),
                (1, 8, 4, 0.5)
            );
        }
        for t in split {
            assert_eq!(t.nodes[0].split.as_ref().unwrap().feature, 0);
        }
    }

    #[test]
    fn a_tree_draws_node_features_independently_of_its_subsample() {
        // Sixteen samples, sample 0 the only resistant one, and sixteen
        // copies of a feature present on it alone: any drawn feature splits
        // a root that holds sample 0, and the tie goes to the lowest drawn
        // index. Each tree subsamples 4 samples and each node draws 4
        // features. Drawn from one random stream, both draws would pick the
        // same four indices, so every tree holding sample 0 would split on
        // feature 0.
        let n = 16;
        let x = matrix(&vec![&[0usize][..]; n], n);
        let y = Bits::from_indices(n, [0]);
        let forest = fit_forest(
            &x,
            &y,
            &Bits::ones(n),
            &ForestParams {
                n_trees: 200,
                subsample: 0.25,
                ..Default::default()
            },
        );
        let roots: Vec<u32> = forest
            .trees
            .iter()
            .filter(|t| t.nodes[0].n_resistant == 1)
            .map(|t| t.nodes[0].split.as_ref().unwrap().feature)
            .collect();
        assert!(roots.len() > 20, "{} trees hold sample 0", roots.len());
        assert!(roots.iter().any(|&f| f != 0), "{roots:?}");
    }

    #[test]
    fn balanced_nodes_score_exactly_one_half() {
        // 11 R and 4 S with no candidate feature: the root is the only leaf,
        // and its class-weighted resistant share is exactly 1/2.
        let x = BitMatrix::with_cols(15);
        let y = Bits::from_indices(15, 0..11);
        let tree = fit_tree(
            &x,
            &y,
            &Bits::ones(15),
            &TreeParams {
                max_depth: 2,
                min_leaf: 1,
                max_features: None,
                seed: 0,
            },
        );
        assert_eq!(tree.predict(&x, 2), vec![0.5; 15]);
    }

    #[test]
    fn forest_scores_are_exact_vote_shares() {
        // Six of twelve single-leaf trees vote resistant.
        let leaf = |value| Tree {
            nodes: vec![Node {
                value,
                n: 1,
                n_resistant: 0,
                depth: 0,
                split: None,
            }],
        };
        let forest = Forest {
            trees: (0..12)
                .map(|t| leaf(if t < 6 { 1.0 } else { 0.0 }))
                .collect(),
        };
        assert_eq!(forest.predict(&BitMatrix::with_cols(3)), vec![0.5; 3]);
    }
}
