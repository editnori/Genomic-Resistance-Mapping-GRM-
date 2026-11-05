//! Offline confidence intervals for the error rates a report already holds.
//!
//! `grm intervals` reads each dataset report that `index.json` lists and
//! writes `intervals/<id>.json` beside it: for every learner row of every
//! split, a two-sided 95% Clopper-Pearson (exact binomial) interval for VME,
//! the false-susceptible share of the resistant isolates called (x = fn over
//! tp + fn), and for ME, the false-resistant share of the susceptible isolates
//! called (x = fp over tn + fp). It is arithmetic on the saved counts only:
//! nothing is fitted, and no report file changes. Each file records the
//! SHA-256 of the report it was computed from and the counts it used, so a
//! reader can tell when the report has changed since. A zero denominator
//! gives no interval (null), never [0, 1].
//!
//! These intervals describe sampling of these isolates under internal
//! cross-validation. They are not calibration, clinical bounds or an
//! external validation.

use anyhow::{Context, Result, bail};
use serde_json::{Value, json};
use std::collections::BTreeSet;
use std::path::Path;

/// Version of the intervals file layout.
pub const SCHEMA: u32 = 1;
/// Two-sided error rate: a 95% interval.
pub const ALPHA: f64 = 0.05;

/// ln Γ(x) for x > 0 (Lanczos, g = 7, nine terms; relative error near 1e-15).
fn ln_gamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // Reflection keeps the series in its accurate range.
        let pi = std::f64::consts::PI;
        return (pi / (pi * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut a = C[0];
    let t = x + G + 0.5;
    for (i, c) in C.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// Continued fraction for the incomplete beta function (modified Lentz).
fn beta_cf(a: f64, b: f64, x: f64) -> f64 {
    const TINY: f64 = 1e-300;
    let (qab, qap, qam) = (a + b, a + 1.0, a - 1.0);
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < TINY {
        d = TINY;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=10_000 {
        let m = m as f64;
        let m2 = 2.0 * m;
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        h *= d * c;
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-16 {
            break;
        }
    }
    h
}

/// The regularized incomplete beta function I_x(a, b), the CDF of Beta(a, b).
pub fn beta_cdf(x: f64, a: f64, b: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let front =
        (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        front * beta_cf(a, b, x) / a
    } else {
        1.0 - front * beta_cf(b, a, 1.0 - x) / b
    }
}

/// The p quantile of Beta(a, b), by bisection on the CDF: it is monotone,
/// so this converges to the last bit the CDF resolves.
pub fn beta_quantile(p: f64, a: f64, b: f64) -> f64 {
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            break;
        }
        if beta_cdf(mid, a, b) < p {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Two-sided Clopper-Pearson interval for x events in n trials at level
/// 1 - alpha, or None when n is 0. Lower is 0 at x = 0 and upper is 1 at x = n.
pub fn clopper_pearson(x: u64, n: u64, alpha: f64) -> Option<(f64, f64)> {
    if n == 0 || x > n {
        return None;
    }
    let (xf, nf) = (x as f64, n as f64);
    let lo = if x == 0 {
        0.0
    } else {
        beta_quantile(alpha / 2.0, xf, nf - xf + 1.0)
    };
    let hi = if x == n {
        1.0
    } else {
        beta_quantile(1.0 - alpha / 2.0, xf + 1.0, nf - xf)
    };
    Some((lo, hi))
}

fn rate(x: u64, n: u64) -> Value {
    match clopper_pearson(x, n, ALPHA) {
        None => Value::Null,
        Some((lo, hi)) => json!({ "x": x, "n": n, "lo": lo, "hi": hi }),
    }
}

fn count(metrics: &Value, key: &str, id: &str) -> Result<u64> {
    metrics[key]
        .as_u64()
        .with_context(|| format!("{id}: a learner row's metrics have no whole-number {key}"))
}

/// The intervals file for one dataset report: its id, the report's SHA-256,
/// its evaluation revision (null when it records none) and, per split and
/// learner row, the counts used with the VME and ME intervals. A row that
/// called no isolate (null metrics) is left out.
pub fn for_report(bytes: &[u8]) -> Result<Value> {
    let report: Value = serde_json::from_slice(bytes).context("parse a dataset report")?;
    let id = report["id"]
        .as_str()
        .context("a dataset report has no id")?;
    let mut splits = Vec::new();
    for split in report["splits"].as_array().into_iter().flatten() {
        let mut rows = Vec::new();
        for row in split["learners"].as_array().into_iter().flatten() {
            let m = &row["metrics"];
            if m.is_null() {
                continue;
            }
            let [n, tp, fp, tn, fn_] = ["n", "tp", "fp", "tn", "fn"].map(|k| count(m, k, id));
            let (n, tp, fp, tn, fn_) = (n?, tp?, fp?, tn?, fn_?);
            if tp + fp + tn + fn_ != n {
                bail!("{id}: a learner row's tp, fp, tn and fn do not sum to n");
            }
            rows.push(json!({
                "learner": row["learner"],
                "counts": { "n": n, "tp": tp, "fp": fp, "tn": tn, "fn": fn_ },
                "vme": rate(fn_, tp + fn_),
                "me": rate(fp, tn + fp),
            }));
        }
        splits.push(json!({ "scheme": split["scheme"], "learners": rows }));
    }
    Ok(json!({
        "schema": SCHEMA,
        "method": "clopper-pearson",
        "confidence": 1.0 - ALPHA,
        "dataset": id,
        "report_sha256": crate::research::sha256(bytes),
        "evaluation_revision": report["evaluation"]["revision"],
        "splits": splits,
    }))
}

/// Writes `intervals/<id>.json` under `report_dir` for every dataset that
/// its `index.json` lists, removes any other file there, and returns how
/// many it wrote. It holds the report directory's lock, as `grm report`
/// does, so it never reads a build half done.
pub fn write(report_dir: &Path) -> Result<usize> {
    let lock = std::fs::File::open(report_dir)
        .with_context(|| format!("open {}", report_dir.display()))?;
    lock.lock()
        .with_context(|| format!("lock {}", report_dir.display()))?;
    let index: Value = serde_json::from_slice(
        &std::fs::read(report_dir.join("index.json")).context("run `grm report` first")?,
    )?;
    let out = report_dir.join("intervals");
    std::fs::create_dir_all(&out)?;
    let mut kept = BTreeSet::new();
    for row in index["datasets"].as_array().into_iter().flatten() {
        let id = row["id"]
            .as_str()
            .context("index.json lists a dataset with no id")?;
        let path = report_dir.join("datasets").join(format!("{id}.json"));
        let bytes = std::fs::read(&path).with_context(|| format!("read {}", path.display()))?;
        let file = for_report(&bytes)?;
        if file["dataset"] != id {
            bail!("{} records the id {}", path.display(), file["dataset"]);
        }
        let name = format!("{id}.json");
        crate::report::publish(&out.join(&name), &serde_json::to_vec(&file)?)?;
        kept.insert(name);
    }
    for entry in std::fs::read_dir(&out)? {
        let entry = entry?;
        if !entry.file_name().to_str().is_some_and(|n| kept.contains(n)) {
            std::fs::remove_file(entry.path())
                .with_context(|| format!("remove {}", entry.path().display()))?;
        }
    }
    Ok(kept.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// (x, n, lower, upper) from the offline revision-7 run's SciPy 1.14.1
    /// `beta.ppf` intervals (data/research/grm-all23-saved-intervals-20261027-kimi,
    /// vme-me-clopper-pearson.csv), chosen across zero, full and large counts.
    const SCIPY: [(u64, u64, f64, f64); 30] = [
        (0, 6, 0.0, 0.45925812643990044),
        (0, 7, 0.0, 0.40961639722500337),
        (0, 21, 0.0, 0.1610976152190796),
        (0, 1886, 0.0, 0.00195401602039081),
        (3, 111, 0.005608683216675948, 0.0769555750675793),
        (4, 87, 0.012667375820641409, 0.11355344796666053),
        (5, 65, 0.025448467019043023, 0.17045629048892497),
        (5, 83, 0.019847321241740542, 0.1350082396842335),
        (5, 209, 0.0078124846568847685, 0.05494333888926265),
        (6, 110, 0.020277720969548958, 0.11494662512228843),
        (6, 280, 0.007903552986266925, 0.04605508562773042),
        (7, 1319, 0.0021362919591686753, 0.010903806367681847),
        (8, 101, 0.03481691039879477, 0.15011620100508763),
        (8, 317, 0.010956946291487912, 0.04911870180113821),
        (8, 1886, 0.0018330248422895048, 0.008340805876968928),
        (10, 10, 0.6915028921812392, 1.0),
        (11, 11, 0.7150858470818455, 1.0),
        (11, 120, 0.04665226965524853, 0.15809854119583328),
        (14, 57, 0.14126857773596005, 0.37761030480194585),
        (22, 175, 0.08048892787899473, 0.18412749654984673),
        (23, 451, 0.032599285708675686, 0.07554156137766942),
        (24, 1886, 0.008169922719802096, 0.01887546766118562),
        (26, 80, 0.22445958517745335, 0.4389370456044484),
        (41, 304, 0.09855207074221606, 0.1784867586358306),
        (53, 408, 0.0988407923817262, 0.16645003066433706),
        (63, 848, 0.05755743818575499, 0.09405440618587402),
        (75, 619, 0.09650626678664745, 0.14949904178259882),
        (78, 618, 0.10105919143253055, 0.15500146718542052),
        (83, 161, 0.43556769905190557, 0.5949061282319127),
        (208, 530, 0.3506387002748109, 0.43547319573686855),
    ];

    #[test]
    fn clopper_pearson_matches_scipy_beta_ppf() {
        for (x, n, lo, hi) in SCIPY {
            let (l, h) = clopper_pearson(x, n, ALPHA).unwrap();
            assert!((l - lo).abs() < 1e-12, "{x}/{n}: lower {l} vs {lo}");
            assert!((h - hi).abs() < 1e-12, "{x}/{n}: upper {h} vs {hi}");
        }
    }

    #[test]
    fn clopper_pearson_boundaries_have_closed_forms() {
        for n in [1u64, 2, 6, 25, 400, 1886, 10_000] {
            let (lo, hi) = clopper_pearson(0, n, ALPHA).unwrap();
            assert_eq!(lo, 0.0);
            assert!(
                (hi - (1.0 - (ALPHA / 2.0).powf(1.0 / n as f64))).abs() < 1e-12,
                "0/{n}"
            );
            let (lo, hi) = clopper_pearson(n, n, ALPHA).unwrap();
            assert_eq!(hi, 1.0);
            assert!(
                (lo - (ALPHA / 2.0).powf(1.0 / n as f64)).abs() < 1e-12,
                "{n}/{n}"
            );
        }
        // 1 of 2: Beta(1, 2) and Beta(2, 1) quantiles are 1 - sqrt(1 - p) and sqrt(p).
        let (lo, hi) = clopper_pearson(1, 2, ALPHA).unwrap();
        let root = (1.0 - ALPHA / 2.0).sqrt();
        assert!((lo - (1.0 - root)).abs() < 1e-12 && (hi - root).abs() < 1e-12);
        assert_eq!(clopper_pearson(0, 0, ALPHA), None);
        assert_eq!(clopper_pearson(3, 2, ALPHA), None);
    }

    #[test]
    fn clopper_pearson_tails_hold_alpha_over_two() {
        // The defining identities: P(X >= x | p = lower) = P(X <= x | p = upper) = alpha/2.
        for (x, n) in [(1u64, 9u64), (5, 83), (40, 300), (700, 1900)] {
            let (lo, hi) = clopper_pearson(x, n, ALPHA).unwrap();
            let (xf, nf) = (x as f64, n as f64);
            assert!((beta_cdf(lo, xf, nf - xf + 1.0) - ALPHA / 2.0).abs() < 1e-12);
            assert!((1.0 - beta_cdf(hi, xf + 1.0, nf - xf) - ALPHA / 2.0).abs() < 1e-12);
            assert!(lo < xf / nf && xf / nf < hi);
        }
    }

    fn report_json(rows: Value) -> Vec<u8> {
        serde_json::to_vec(&json!({
            "id": "sp__drug",
            "evaluation": { "revision": 11 },
            "splits": [{ "scheme": "random", "learners": rows }],
        }))
        .unwrap()
    }

    #[test]
    fn intervals_file_pins_the_report_and_its_counts() {
        let bytes = report_json(json!([
            { "learner": "scm", "metrics": { "n": 1278, "tp": 755, "fp": 41, "tn": 452, "fn": 30 } },
            { "learner": "gbt", "metrics": null },
            { "learner": "rf", "metrics": { "n": 20, "tp": 0, "fp": 2, "tn": 18, "fn": 0 } },
        ]));
        let file = for_report(&bytes).unwrap();
        assert_eq!(file["schema"], 1);
        assert_eq!(file["method"], "clopper-pearson");
        assert_eq!(file["confidence"], 0.95);
        assert_eq!(file["dataset"], "sp__drug");
        assert_eq!(file["evaluation_revision"], 11);
        assert_eq!(file["report_sha256"], crate::research::sha256(&bytes));
        let rows = file["splits"][0]["learners"].as_array().unwrap();
        // The row that called nothing has no interval row.
        assert_eq!(rows.len(), 2);
        // The offline run's first row: scm on acinetobacter_baumannii__amikacin.
        assert_eq!(rows[0]["vme"]["x"], 30);
        assert_eq!(rows[0]["vme"]["n"], 785);
        assert!((rows[0]["vme"]["lo"].as_f64().unwrap() - 0.02593041110157624).abs() < 1e-12);
        assert!((rows[0]["vme"]["hi"].as_f64().unwrap() - 0.054109136028033265).abs() < 1e-12);
        assert_eq!(rows[0]["me"]["n"], 493);
        assert!((rows[0]["me"]["lo"].as_f64().unwrap() - 0.0603395944989114).abs() < 1e-12);
        assert!((rows[0]["me"]["hi"].as_f64().unwrap() - 0.11113119109649895).abs() < 1e-12);
        assert_eq!(rows[0]["counts"]["tn"], 452);
        // No resistant isolate called: VME has no denominator, so no interval.
        assert!(rows[1]["vme"].is_null());
        assert_eq!(rows[1]["me"]["x"], 2);
    }

    #[test]
    fn intervals_file_refuses_counts_that_do_not_add_up() {
        let bytes = report_json(json!([
            { "learner": "scm", "metrics": { "n": 10, "tp": 5, "fp": 1, "tn": 1, "fn": 1 } },
        ]));
        assert!(for_report(&bytes).is_err());
        let bytes = report_json(json!([
            { "learner": "scm", "metrics": { "n": 10, "tp": 5.5, "fp": 1, "tn": 3, "fn": 1 } },
        ]));
        assert!(for_report(&bytes).is_err());
    }

    #[test]
    fn write_covers_the_indexed_datasets_and_leaves_reports_untouched() {
        let temp = tempfile::tempdir().unwrap();
        let dir = temp.path();
        std::fs::create_dir_all(dir.join("datasets")).unwrap();
        std::fs::create_dir_all(dir.join("intervals")).unwrap();
        let bytes = report_json(json!([
            { "learner": "scm", "metrics": { "n": 4, "tp": 2, "fp": 0, "tn": 2, "fn": 0 } },
        ]));
        std::fs::write(dir.join("datasets/sp__drug.json"), &bytes).unwrap();
        std::fs::write(
            dir.join("index.json"),
            br#"{"datasets":[{"id":"sp__drug"}]}"#,
        )
        .unwrap();
        std::fs::write(dir.join("intervals/gone__drug.json"), b"{}").unwrap();
        assert_eq!(write(dir).unwrap(), 1);
        assert_eq!(
            std::fs::read(dir.join("datasets/sp__drug.json")).unwrap(),
            bytes
        );
        assert!(!dir.join("intervals/gone__drug.json").exists());
        let file: Value =
            serde_json::from_slice(&std::fs::read(dir.join("intervals/sp__drug.json")).unwrap())
                .unwrap();
        assert_eq!(file, for_report(&bytes).unwrap());
    }
}
