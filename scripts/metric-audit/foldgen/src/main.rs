//! Read-only fold reconstruction for the revision-6 metric audit.
//!
//! Copies `stratified_folds` and `grouped_folds` from
//! src/learn/mod.rs, with `Bits` replaced by `[bool]`, so the audit
//! can rebuild each outer fold from a saved report. Uses the same rand 0.9.5
//! and rand_chacha 0.9.0 the pipeline links.
//!
//! stdin, one request per line:
//!   S <k> <seed> <labels as 0/1 string>
//!   G <k> <seed> <labels as 0/1 string> <comma-separated group ids>
//! stdout: one line of comma-separated fold indices per request.

use rand::SeedableRng;
use rand::seq::SliceRandom;
use rand_chacha::ChaCha8Rng;
use std::io::{self, BufRead, Write};

fn stratified_folds(y: &[bool], k: usize, seed: u64) -> Vec<usize> {
    assert!(k >= 2);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mut fold = vec![0; y.len()];
    for class in [true, false] {
        let mut members: Vec<usize> = (0..y.len()).filter(|&i| y[i] == class).collect();
        members.shuffle(&mut rng);
        for (j, i) in members.into_iter().enumerate() {
            fold[i] = j % k;
        }
    }
    fold
}

fn grouped_folds(y: &[bool], groups: &[u32], k: usize, seed: u64) -> Vec<usize> {
    assert!(k >= 2);
    assert_eq!(y.len(), groups.len());
    let mut members: std::collections::BTreeMap<u32, Vec<usize>> = Default::default();
    for (i, &g) in groups.iter().enumerate() {
        members.entry(g).or_default().push(i);
    }
    let mut order: Vec<(u32, Vec<usize>)> = members.into_iter().collect();
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    order.shuffle(&mut rng);
    order.sort_by_key(|(_, m)| std::cmp::Reverse(m.len()));
    let total_r = y.iter().filter(|&&v| v).count() as f64;
    let total_s = y.len() as f64 - total_r;
    let mut load = vec![(0.0f64, 0.0f64); k];
    let mut fold = vec![0; y.len()];
    for (_, m) in order {
        let r = m.iter().filter(|&&i| y[i]).count() as f64;
        let s = m.len() as f64 - r;
        let best = (0..k)
            .min_by(|&a, &b| {
                let cost = |f: usize| {
                    let (lr, ls) = load[f];
                    ((lr + r).powi(2) - lr.powi(2)) / total_r.max(1.0).powi(2)
                        + ((ls + s).powi(2) - ls.powi(2)) / total_s.max(1.0).powi(2)
                };
                cost(a).total_cmp(&cost(b))
            })
            .unwrap();
        load[best].0 += r;
        load[best].1 += s;
        for i in m {
            fold[i] = best;
        }
    }
    fold
}

fn main() {
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());
    for line in stdin.lock().lines() {
        let line = line.expect("stdin");
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.is_empty() {
            continue;
        }
        let k: usize = parts[1].parse().expect("k");
        let seed: u64 = parts[2].parse().expect("seed");
        let y: Vec<bool> = parts[3].bytes().map(|b| b == b'1').collect();
        let folds = match parts[0] {
            "S" => stratified_folds(&y, k, seed),
            "G" => {
                let groups: Vec<u32> = parts[4].split(',').map(|g| g.parse().expect("group")).collect();
                grouped_folds(&y, &groups, k, seed)
            }
            other => panic!("unknown request {other}"),
        };
        let text: Vec<String> = folds.iter().map(|f| f.to_string()).collect();
        writeln!(out, "{}", text.join(",")).expect("stdout");
        out.flush().expect("flush");
    }
}
