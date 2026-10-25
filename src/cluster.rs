//! Population structure from MinHash sketches.
//!
//! Each genome's sketch is the s smallest mixed k-mer values of its kset.
//! Mash distance between two sketches estimates divergence; average-linkage
//! clustering on those distances gives lineage groups for held-out-lineage
//! evaluation, a dendrogram for display, and PCoA gives a 2D map.

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Jaccard estimate from two bottom-s sketches (sorted ascending).
pub fn jaccard(a: &[u64], b: &[u64], s: usize) -> f64 {
    let (mut i, mut j, mut seen, mut shared) = (0, 0, 0, 0);
    while seen < s && i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                shared += 1;
                i += 1;
                j += 1;
            }
        }
        seen += 1;
    }
    let seen = (seen + (a.len() - i) + (b.len() - j)).min(s);
    if seen == 0 {
        0.0
    } else {
        shared as f64 / seen as f64
    }
}

/// Mash distance (Ondov et al. 2016): -ln(2J / (1 + J)) / k, capped at 1.
pub fn mash_distance(j: f64, k: usize) -> f64 {
    if j <= 0.0 {
        return 1.0;
    }
    (-(2.0 * j / (1.0 + j)).ln() / k as f64).min(1.0)
}

/// Condensed-to-square symmetric distance matrix, row-major n x n.
pub fn distance_matrix(sketches: &[Vec<u64>], s: usize, k: usize) -> Vec<f32> {
    let n = sketches.len();
    let rows: Vec<Vec<f32>> = (0..n)
        .into_par_iter()
        .map(|i| {
            (0..n)
                .map(|j| {
                    if i == j {
                        0.0
                    } else {
                        mash_distance(jaccard(&sketches[i], &sketches[j], s), k) as f32
                    }
                })
                .collect()
        })
        .collect();
    rows.into_iter().flatten().collect()
}

/// One agglomeration step: clusters `a` and `b` (dendrogram node ids: leaves
/// are 0..n, merge i creates node n + i) joined at `height`.
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct Merge {
    pub a: u32,
    pub b: u32,
    pub height: f32,
    pub size: u32,
}

/// Average-linkage (UPGMA) clustering by the nearest-neighbor chain
/// algorithm, O(n^2) time. Returns merges sorted by height.
pub fn average_linkage(dist: &[f32], n: usize) -> Vec<Merge> {
    if n < 2 {
        return vec![];
    }
    let mut d = dist.to_vec();
    let mut size = vec![1u32; n];
    // Height of the merge that formed each active cluster (-inf for a leaf).
    let mut height = vec![f32::NEG_INFINITY; n];
    let mut active = vec![true; n];
    let mut raw: Vec<(usize, usize, f32)> = Vec::with_capacity(n - 1);
    let mut chain: Vec<usize> = Vec::new();
    let mut remaining = n;
    while remaining > 1 {
        if chain.is_empty() {
            chain.push(active.iter().position(|&a| a).unwrap());
        }
        let c = *chain.last().unwrap();
        let prev = if chain.len() >= 2 {
            Some(chain[chain.len() - 2])
        } else {
            None
        };
        // Nearest active neighbor; prefer the previous chain element on ties
        // so the chain cannot cycle.
        let mut best = usize::MAX;
        let mut best_d = f32::INFINITY;
        for j in 0..n {
            if j == c || !active[j] {
                continue;
            }
            let dj = d[c * n + j];
            if dj < best_d || (dj == best_d && Some(j) == prev) {
                best = j;
                best_d = dj;
            }
        }
        if Some(best) == prev {
            chain.pop();
            chain.pop();
            let (keep, gone) = (c.min(best), c.max(best));
            // f32 rounding in the update below can put a merge a few ulps
            // under the merge that formed one of its clusters (tied distances
            // averaged together). The height sort would then order the parent
            // first and reuse a stale node, so clamp to the children.
            let h = [best_d, height[keep], height[gone]]
                .into_iter()
                .max_by(f32::total_cmp)
                .unwrap();
            raw.push((keep, gone, h));
            height[keep] = h;
            let (sk, sg) = (size[keep] as f32, size[gone] as f32);
            for j in 0..n {
                if active[j] && j != keep && j != gone {
                    let v = (sk * d[keep * n + j] + sg * d[gone * n + j]) / (sk + sg);
                    d[keep * n + j] = v;
                    d[j * n + keep] = v;
                }
            }
            size[keep] += size[gone];
            active[gone] = false;
            remaining -= 1;
        } else {
            chain.push(best);
        }
    }
    // Average linkage is reducible and the clamp above keeps that true in
    // f32, so heights are monotone along the tree. A stable sort by height
    // keeps each child before its parent and gives a valid agglomeration order.
    raw.sort_by(|a, b| a.2.total_cmp(&b.2));
    let mut node_of: Vec<u32> = (0..n as u32).collect();
    let mut sizes: Vec<u32> = vec![1; n];
    raw.into_iter()
        .enumerate()
        .map(|(i, (keep, gone, h))| {
            let (a, b) = (node_of[keep], node_of[gone]);
            let s = sizes[keep] + sizes[gone];
            node_of[keep] = (n + i) as u32;
            sizes[keep] = s;
            Merge {
                a,
                b,
                height: h,
                size: s,
            }
        })
        .collect()
}

/// Flat clusters from merges at or below `threshold`; labels are 0..k in
/// order of first appearance.
pub fn cut(merges: &[Merge], n: usize, threshold: f32) -> Vec<u32> {
    let mut parent: Vec<usize> = (0..n + merges.len()).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    for (i, m) in merges.iter().enumerate() {
        if m.height <= threshold {
            let node = n + i;
            let (ra, rb) = (
                find(&mut parent, m.a as usize),
                find(&mut parent, m.b as usize),
            );
            parent[ra] = node;
            parent[rb] = node;
        }
    }
    let mut label_of = std::collections::HashMap::new();
    (0..n)
        .map(|i| {
            let root = find(&mut parent, i);
            let next = label_of.len() as u32;
            *label_of.entry(root).or_insert(next)
        })
        .collect()
}

/// Leaf order of the dendrogram (left-to-right), for display.
pub fn leaf_order(merges: &[Merge], n: usize) -> Vec<u32> {
    if merges.is_empty() {
        return (0..n as u32).collect();
    }
    let mut out = Vec::with_capacity(n);
    let mut roots: Vec<u32> = Vec::new();
    // Nodes never used as a child are roots (a forest if merges stop early).
    let mut used = vec![false; n + merges.len()];
    for m in merges {
        used[m.a as usize] = true;
        used[m.b as usize] = true;
    }
    for id in (0..n + merges.len()).rev() {
        if !used[id] {
            roots.push(id as u32);
        }
    }
    for root in roots {
        let mut stack = vec![root];
        while let Some(id) = stack.pop() {
            if (id as usize) < n {
                out.push(id);
            } else {
                let m = merges[id as usize - n];
                stack.push(m.b);
                stack.push(m.a);
            }
        }
    }
    out
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Lineages {
    pub threshold: f32,
    pub labels: Vec<u32>,
    pub n_clusters: usize,
    pub largest_share: f64,
}

pub const THRESHOLDS: [f32; 7] = [0.01, 0.005, 0.0025, 0.001, 0.0005, 0.00025, 0.0001];

/// The coarsest threshold whose clusters are numerous and small enough to
/// spread across folds: at least `min_clusters`, none above `max_share`.
pub fn choose_lineages(
    merges: &[Merge],
    n: usize,
    min_clusters: usize,
    max_share: f64,
) -> Lineages {
    let mut last = None;
    for &t in &THRESHOLDS {
        let labels = cut(merges, n, t);
        let k = labels.iter().max().map_or(0, |m| *m as usize + 1);
        let mut counts = vec![0usize; k];
        for &l in &labels {
            counts[l as usize] += 1;
        }
        let largest_share = counts.iter().copied().max().unwrap_or(0) as f64 / n.max(1) as f64;
        let lin = Lineages {
            threshold: t,
            labels,
            n_clusters: k,
            largest_share,
        };
        if k >= min_clusters && largest_share <= max_share {
            return lin;
        }
        last = Some(lin);
    }
    last.unwrap()
}

/// Classical multidimensional scaling (PCoA) to two dimensions by power
/// iteration on the double-centered squared-distance matrix.
pub fn pcoa(dist: &[f32], n: usize) -> Vec<[f32; 2]> {
    if n < 3 {
        return vec![[0.0, 0.0]; n];
    }
    let mut b = vec![0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            let v = dist[i * n + j] as f64;
            b[i * n + j] = -0.5 * v * v;
        }
    }
    let row_mean: Vec<f64> = (0..n)
        .map(|i| b[i * n..(i + 1) * n].iter().sum::<f64>() / n as f64)
        .collect();
    let grand = row_mean.iter().sum::<f64>() / n as f64;
    for i in 0..n {
        for j in 0..n {
            b[i * n + j] += grand - row_mean[i] - row_mean[j];
        }
    }
    let tol = 1e-6 * b.iter().map(|x| x * x).sum::<f64>().sqrt();
    let mut coords = vec![[0f32; 2]; n];
    let mut found: Vec<(f64, Vec<f64>)> = Vec::new();
    let mut axis = 0;
    for trial in 0..n {
        if axis == 2 {
            break;
        }
        let v: Vec<f64> = (0..n)
            .map(|i| ((i * 7919 + trial * 104_729) % 1000) as f64 / 1000.0 - 0.5)
            .collect();
        let (lambda, rayleigh, v) = power_iteration(&b, n, &found, v, None);
        if rayleigh >= lambda - tol {
            for i in 0..n {
                coords[i][axis] = (v[i] * lambda.max(0.0).sqrt()) as f32;
            }
            found.push((lambda, v));
            axis += 1;
        } else {
            let (flipped, _, u) = power_iteration(&b, n, &found, v, Some(lambda));
            found.push((lambda - flipped, u));
        }
    }
    coords
}

fn power_iteration(
    b: &[f64],
    n: usize,
    found: &[(f64, Vec<f64>)],
    mut v: Vec<f64>,
    flip: Option<f64>,
) -> (f64, f64, Vec<f64>) {
    let mut lambda = 0.0;
    let mut rayleigh = 0.0;
    for _ in 0..300 {
        let mut w: Vec<f64> = (0..n)
            .into_par_iter()
            .map(|i| (0..n).map(|j| b[i * n + j] * v[j]).sum())
            .collect();
        for (l, u) in found {
            let dot: f64 = u.iter().zip(&v).map(|(a, b)| a * b).sum();
            for (wi, ui) in w.iter_mut().zip(u) {
                *wi -= l * dot * ui;
            }
        }
        if let Some(rho) = flip {
            for (wi, vi) in w.iter_mut().zip(&v) {
                *wi = rho * vi - *wi;
            }
        }
        let norm = w.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm == 0.0 {
            break;
        }
        rayleigh = w.iter().zip(&v).map(|(a, b)| a * b).sum();
        let next: Vec<f64> = w.iter().map(|x| x / norm).collect();
        let delta: f64 = next.iter().zip(&v).map(|(a, b)| (a - b).abs()).sum();
        v = next;
        lambda = norm;
        if delta < 1e-9 {
            break;
        }
    }
    (lambda, rayleigh, v)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jaccard_of_overlapping_sketches() {
        let a: Vec<u64> = (0..100).collect();
        let b: Vec<u64> = (50..150).collect();
        // Bottom-100 of the union is 0..100; b holds 50 of those.
        assert!((jaccard(&a, &b, 100) - 0.5).abs() < 1e-12);
        assert_eq!(mash_distance(1.0, 21), 0.0);
        assert_eq!(mash_distance(0.0, 21), 1.0);
    }

    #[test]
    fn short_sketches_count_the_whole_union() {
        assert_eq!(jaccard(&[1, 2], &[1, 5, 6, 7], 100), 0.2);
        assert_eq!(jaccard(&[1, 2], &[1, 5, 6, 7], 3), 1.0 / 3.0);
    }

    #[test]
    fn average_linkage_recovers_two_groups() {
        // Points on a line: {0, 1, 2} and {10, 11}.
        let pos = [0.0f32, 1.0, 2.0, 10.0, 11.0];
        let n = pos.len();
        let dist: Vec<f32> = (0..n * n)
            .map(|k| (pos[k / n] - pos[k % n]).abs())
            .collect();
        let merges = average_linkage(&dist, n);
        assert_eq!(merges.len(), n - 1);
        assert!(merges.windows(2).all(|w| w[0].height <= w[1].height));
        assert_eq!(merges.last().unwrap().size, 5);
        let labels = cut(&merges, n, 3.0);
        assert_eq!(labels, vec![0, 0, 0, 1, 1]);
        let order = leaf_order(&merges, n);
        let mut sorted = order.clone();
        sorted.sort();
        assert_eq!(sorted, vec![0, 1, 2, 3, 4]);
        // The two groups stay contiguous in leaf order.
        let pos_of = |leaf: u32| order.iter().position(|&l| l == leaf).unwrap();
        assert_eq!((pos_of(3) as i32 - pos_of(4) as i32).abs(), 1);
        // The final UPGMA height is the mean pairwise distance between groups.
        let between: f32 = [0.0f32, 1.0, 2.0]
            .iter()
            .flat_map(|a| [10.0f32, 11.0].map(move |b| b - a))
            .sum::<f32>()
            / 6.0;
        assert!((merges.last().unwrap().height - between).abs() < 1e-4);
    }

    #[test]
    fn tied_distances_still_form_one_tree() {
        // 14 genomes in 9 clonal groups, every pair across groups at the Mash
        // distance for J = 0.963, k = 21. Averaging those ties in f32 drifts
        // later merges one ulp below the merges that formed their clusters.
        let h = f32::from_bits(0x3a6d_88f2);
        let groups = [0, 1, 2, 3, 1, 4, 5, 6, 7, 8, 8, 3, 2, 8];
        let n = groups.len();
        let dist: Vec<f32> = (0..n * n)
            .map(|k| {
                if groups[k / n] == groups[k % n] {
                    0.0
                } else {
                    h
                }
            })
            .collect();
        let merges = average_linkage(&dist, n);
        let mut used = vec![false; n + merges.len()];
        for (i, m) in merges.iter().enumerate() {
            for child in [m.a as usize, m.b as usize] {
                assert!(
                    child < n + i && !used[child],
                    "merge {i} reuses node {child}"
                );
                used[child] = true;
            }
        }
        let mut order = leaf_order(&merges, n);
        order.sort();
        assert_eq!(order, (0..n as u32).collect::<Vec<_>>());
    }

    #[test]
    fn identical_genomes_keep_their_mash_height() {
        let d = mash_distance(1.0, 21) as f32;
        let merges = average_linkage(&[0.0, d, d, 0.0], 2);
        assert_eq!(merges[0].height.to_bits(), d.to_bits());
    }

    #[test]
    fn pcoa_ignores_negative_eigenvalues() {
        let n = 4;
        let mut dist = vec![0f32; n * n];
        for (i, j) in [(0, 1), (1, 0), (2, 3), (3, 2)] {
            dist[i * n + j] = 1.0;
        }
        let c = pcoa(&dist, n);
        for i in 0..n {
            for j in i + 1..n {
                let got = ((c[i][0] - c[j][0]).powi(2) + (c[i][1] - c[j][1]).powi(2)).sqrt();
                let want = if dist[i * n + j] == 1.0 {
                    1.0
                } else {
                    0.5f32.sqrt()
                };
                assert!((got - want).abs() < 1e-4, "{i}-{j}: {c:?}");
            }
        }
    }

    #[test]
    fn pcoa_preserves_a_line() {
        let pos = [0.0f32, 1.0, 2.0, 3.0];
        let n = pos.len();
        let dist: Vec<f32> = (0..n * n)
            .map(|k| (pos[k / n] - pos[k % n]).abs())
            .collect();
        let c = pcoa(&dist, n);
        let d01 = ((c[0][0] - c[1][0]).powi(2) + (c[0][1] - c[1][1]).powi(2)).sqrt();
        let d03 = ((c[0][0] - c[3][0]).powi(2) + (c[0][1] - c[3][1]).powi(2)).sqrt();
        assert!(
            (d01 - 1.0).abs() < 1e-3 && (d03 - 3.0).abs() < 1e-3,
            "{d01} {d03}"
        );
    }
}
