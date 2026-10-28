//! The fitting process of the large-species lasso scan. `scan.py` beside this
//! file runs it and applies the registered stop rule; see
//! docs/lasso-stopping.md.
//!
//! usage: lasso-scan [--inputs-only] DATA_ROOT SPECIES DATASET...
//!
//! It runs `grm::analysis::lasso_scan`, the shipped lasso fits of the
//! evaluator with the default evaluation settings, on each named dataset of
//! SPECIES in turn. It writes one JSON line per event to stdout: a dataset's
//! inputs and input fingerprint, then a start line and a fit line for each of
//! its lasso fits. `--inputs-only` stops after the inputs lines. It reads
//! DATA_ROOT and writes nothing else. A lasso error is written as a fit line
//! and then ends the run with exit status 1.

use anyhow::{Result, bail};
use grm::analysis::lasso_scan;
use grm::layout::Layout;
use grm::learn::eval::EvalConfig;
use std::io::Write;

const USAGE: &str = "usage: lasso-scan [--inputs-only] DATA_ROOT SPECIES DATASET...";

fn main() -> Result<()> {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let fit = args.first().map(String::as_str) != Some("--inputs-only");
    if !fit {
        args.remove(0);
    }
    let [root, species, datasets @ ..] = &args[..] else {
        bail!(USAGE);
    };
    if datasets.is_empty() || args.iter().any(|a| a.starts_with("--")) {
        bail!(USAGE);
    }
    let mut out = std::io::stdout().lock();
    lasso_scan(
        &Layout::new(root),
        species,
        datasets,
        &EvalConfig::default(),
        fit,
        &mut |dataset, event| {
            let mut line = serde_json::to_value(event)?;
            line["dataset"] = dataset.into();
            writeln!(out, "{line}")?;
            out.flush()?;
            Ok(())
        },
    )
}
