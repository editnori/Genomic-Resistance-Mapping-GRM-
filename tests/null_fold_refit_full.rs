//! A FILE that takes no line stops the evaluation (`grm learn
//! --exp-null-refit`, docs/null-fold-refit.md): a skipped record would bias
//! every count taken from FILE. The error names FILE, the dataset, the
//! split, the learner and the fold. Its own test binary, because the flag is
//! set once per process. It writes to /dev/full, which Linux provides and
//! which fails every write with "No space left on device"; on a system
//! without it the test is ignored, so it neither passes nor fails there.

use grm::bits::{BitMatrix, Bits};
use grm::learn::eval::{EvalConfig, evaluate};
use grm::learn::gbt::GbtParams;
use grm::learn::null_refit;
use grm::learn::tree::ForestParams;

#[test]
#[cfg_attr(
    not(target_os = "linux"),
    ignore = "needs /dev/full, a file every write to which fails"
)]
fn a_file_that_takes_no_line_stops_the_evaluation_and_the_error_names_the_fold() {
    let cfg = EvalConfig {
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
    };
    // No feature at all, so every boosting and lasso fold is null; the first
    // to be refitted is boosting's on the random split's outer fold 0.
    let n = 120;
    let y = Bits::from_indices(n, (0..n).filter(|i| i % 4 == 0 || [1, 5, 9].contains(i)));
    let groups: Vec<u32> = (0..n).map(|i| (i / 4) as u32).collect();
    let empty = BitMatrix::with_cols(n);
    let run = || {
        evaluate(
            "fixture__fixturemycin",
            &empty,
            &y,
            &groups,
            &[],
            3,
            None,
            &cfg,
        )
    };
    assert!(run().is_ok());
    null_refit::enable(std::path::Path::new("/dev/full")).unwrap();
    let message = format!("{:#}", run().unwrap_err());
    assert!(
        message.starts_with(
            "appending the null-fold record of fixture__fixturemycin, random split, gbt, \
             outer fold 0 to /dev/full: "
        ),
        "{message}"
    );
    assert!(message.contains("os error 28"), "{message}");
}
