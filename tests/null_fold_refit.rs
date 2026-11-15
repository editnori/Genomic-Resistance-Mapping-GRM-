//! The null-fold experiment flag (`grm learn --exp-null-refit`) changes no
//! evaluation output and writes one record per null fold. Its own test
//! binary, because the flag is set once per process.

use grm::bits::{BitMatrix, Bits};
use grm::learn::eval::{EvalConfig, EvalOutput, evaluate};
use grm::learn::gbt::GbtParams;
use grm::learn::null_refit;
use grm::learn::tree::ForestParams;

fn config() -> EvalConfig {
    EvalConfig {
        folds: 3,
        inner_folds: 2,
        scm_max_rules: 1,
        cart_max_depth: 1,
        forest: ForestParams {
            n_trees: 2,
            max_depth: 2,
            ..Default::default()
        },
        gbt: GbtParams {
            rounds: 30,
            max_depth: 1,
            patience: 3,
            ..Default::default()
        },
        lasso_lambdas: 8,
        ..Default::default()
    }
}

/// Every number of an evaluation but wall-clock seconds.
fn without_seconds(out: &EvalOutput) -> serde_json::Value {
    fn strip(v: &mut serde_json::Value) {
        match v {
            serde_json::Value::Object(m) => {
                m.remove("seconds");
                m.values_mut().for_each(strip);
            }
            serde_json::Value::Array(a) => a.iter_mut().for_each(strip),
            _ => {}
        }
    }
    let mut v = serde_json::to_value(out).unwrap();
    strip(&mut v);
    v
}

#[test]
fn the_flag_changes_no_output_and_records_each_null_fold_once() {
    // 120 samples in 30 groups of 4. Two datasets: no feature at all, where
    // every boosting and lasso fold is null, and one with a feature that
    // marks resistant samples in a third of the groups only.
    let n = 120;
    let y = Bits::from_indices(n, (0..n).filter(|i| i % 4 == 0 || [1, 5, 9].contains(i)));
    let groups: Vec<u32> = (0..n).map(|i| (i / 4) as u32).collect();
    let mut partial = BitMatrix::with_cols(n);
    partial.push_row(Bits::from_indices(n, (0..n).filter(|i| i % 4 == 0 && i % 3 == 0)).words());
    partial.push_row(Bits::from_indices(n, (0..n).filter(|i| i % 5 == 0)).words());
    let cfg = config();
    let run = |x: &BitMatrix| {
        evaluate("fixture__fixturemycin", x, &y, &groups, &[], 3, None, &cfg).unwrap()
    };
    let empty = BitMatrix::with_cols(n);
    assert!(!null_refit::enabled());
    let before: Vec<EvalOutput> = [&empty, &partial].map(run).into();

    let dir = tempfile::tempdir().unwrap();
    let file = dir.path().join("null-refit.jsonl");
    null_refit::enable(&file).unwrap();
    assert!(null_refit::enabled());
    assert!(null_refit::enable(&file).is_err());
    let after: Vec<EvalOutput> = [&empty, &partial].map(run).into();
    for (b, a) in before.iter().zip(&after) {
        assert_eq!(without_seconds(b), without_seconds(a));
    }

    // One record per null boosting or lasso outer fold, in evaluation order,
    // naming that fold's test samples.
    let lines: Vec<serde_json::Value> = std::fs::read_to_string(&file)
        .unwrap()
        .lines()
        .map(|l| serde_json::from_str(l).unwrap())
        .collect();
    let mut want = vec![];
    for out in &after {
        for split in &out.splits {
            for fold in 0..cfg.folds {
                for l in &split.learners {
                    if l.fold_complexity.get(fold) == Some(&0.0)
                        && (l.learner == "gbt" || l.learner == "lasso")
                    {
                        let test: Vec<usize> =
                            (0..n).filter(|&i| split.test_fold[i] == fold).collect();
                        want.push((split.scheme.clone(), fold, l.learner.clone(), test));
                    }
                }
            }
        }
    }
    let got: Vec<(String, usize, String, Vec<usize>)> = lines
        .iter()
        .map(|r| {
            (
                r["scheme"].as_str().unwrap().to_string(),
                r["fold"].as_u64().unwrap() as usize,
                r["learner"].as_str().unwrap().to_string(),
                serde_json::from_value(r["test"].clone()).unwrap(),
            )
        })
        .collect();
    assert_eq!(got, want);
    // The empty dataset alone gives 2 splits x 3 folds x 2 learners.
    assert!(want.len() >= 12, "{}", want.len());
    for r in &lines {
        let test = r["test"].as_array().unwrap().len();
        match r["learner"].as_str().unwrap() {
            "lasso" => {
                assert!(r["gbt"].is_null());
                for p in r["lasso"]["points"].as_array().unwrap() {
                    assert_eq!(p["test_scores"].as_array().unwrap().len(), test);
                }
            }
            _ => {
                assert!(r["lasso"].is_null());
                let g = &r["gbt"];
                assert_eq!(g["full_test_scores"].as_array().unwrap().len(), test);
                assert_eq!(
                    g["holdout_loss"].as_array().unwrap().len(),
                    cfg.gbt.rounds + 1
                );
            }
        }
    }
}
