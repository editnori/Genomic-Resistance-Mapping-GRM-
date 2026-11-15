//! Outer folds read from a file (`grm learn --fold-file`, docs/fold-file.md).
//!
//! `evaluate` draws its two outer partitions itself. A fold file supplies
//! one instead, for a question the random and lineage draws cannot ask, such
//! as a split by the day each genome became public. It names one dataset and
//! says, for every genome of that dataset's cohort, whether it is tested in
//! an outer fold, trains in every fold, or is left out:
//!
//! ```json
//! {
//!   "format": "grm-fold-file-v1",
//!   "dataset": "fixture__meropenem",
//!   "scheme": "first-public-past",
//!   "test": [["g7", "g8"]],
//!   "train": ["g0", "g1", "g2"],
//!   "unused": ["g3"]
//! }
//! ```
//!
//! `test` holds one list per outer fold. Fold f is scored on `test[f]` and
//! trains on `train` and on every other fold's test genomes. A genome under
//! `unused` is in no fold.
//!
//! The file is used as written or not at all. Every cohort genome must be
//! listed exactly once, so a file made for another cohort, or one that lost a
//! line, is refused instead of scored on fewer genomes. Its SHA-256 is of the
//! file's bytes, and the result records it.

use super::eval::{DRAWN_SCHEMES, Role};
use anyhow::{Context, Result, bail, ensure};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

/// The `format` of every fold file this module reads.
pub const FOLD_FILE_FORMAT: &str = "grm-fold-file-v1";

/// The longest `scheme` name a fold file may give its partition.
const SCHEME_MAX: usize = 48;

/// A fold file's contents. A key this type does not name is refused.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FoldFile {
    /// `FOLD_FILE_FORMAT`.
    pub format: String,
    /// The dataset id the partition is over.
    pub dataset: String,
    /// The partition's name in the result: lowercase letters, digits and
    /// hyphens, and neither of the names `evaluate` writes.
    pub scheme: String,
    /// Per outer fold, the genomes it scores.
    pub test: Vec<Vec<String>>,
    /// Genomes that train in every outer fold and are never scored.
    pub train: Vec<String>,
    /// Genomes in no fold.
    pub unused: Vec<String>,
}

/// A fold file with the SHA-256 of the bytes it was read from.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LoadedFoldFile {
    pub file: FoldFile,
    /// Lowercase hex SHA-256 of the file as it is on disk.
    pub sha256: String,
}

/// Reads and checks the fold file at `path`.
pub fn load(path: &Path) -> Result<LoadedFoldFile> {
    let bytes = std::fs::read(path).with_context(|| format!("read {}", path.display()))?;
    parse(&bytes).with_context(|| format!("fold file {}", path.display()))
}

/// Parses a fold file's bytes and checks its header (`FoldFile::check`).
pub fn parse(bytes: &[u8]) -> Result<LoadedFoldFile> {
    let file: FoldFile = serde_json::from_slice(bytes).context("parse")?;
    file.check()?;
    Ok(LoadedFoldFile {
        file,
        sha256: crate::research::sha256(bytes),
    })
}

impl FoldFile {
    /// Fails unless the format is `FOLD_FILE_FORMAT`, the scheme is a name a
    /// result can carry, and the file has an outer fold.
    pub fn check(&self) -> Result<()> {
        ensure!(
            self.format == FOLD_FILE_FORMAT,
            "format {:?} is not {FOLD_FILE_FORMAT}",
            self.format
        );
        ensure!(
            !self.scheme.is_empty()
                && self.scheme.len() <= SCHEME_MAX
                && self
                    .scheme
                    .bytes()
                    .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'-'),
            "scheme {:?} is not 1 to {SCHEME_MAX} lowercase letters, digits and hyphens",
            self.scheme
        );
        ensure!(
            !DRAWN_SCHEMES.contains(&self.scheme.as_str()),
            "scheme {:?} is a split `grm learn` draws itself; a fold file needs another name",
            self.scheme
        );
        ensure!(!self.test.is_empty(), "the file has no test fold");
        Ok(())
    }

    /// The outer folds the file holds.
    pub fn folds(&self) -> usize {
        self.test.len()
    }

    /// Each sample's role, in `sample_ids` order, for the cohort of
    /// `dataset`. Fails when the file is for another dataset, names a genome
    /// outside the cohort, lists a genome twice, leaves a cohort genome out,
    /// or has an empty test fold.
    pub fn roles(&self, dataset: &str, sample_ids: &[String]) -> Result<Vec<Role>> {
        self.check()?;
        ensure!(
            self.dataset == dataset,
            "the fold file is for dataset {:?}, not {dataset:?}",
            self.dataset
        );
        let at: HashMap<&str, usize> = sample_ids
            .iter()
            .enumerate()
            .map(|(i, id)| (id.as_str(), i))
            .collect();
        ensure!(
            at.len() == sample_ids.len(),
            "{dataset} lists a sample twice"
        );
        let mut roles: Vec<Option<Role>> = vec![None; sample_ids.len()];
        let lists = self
            .test
            .iter()
            .enumerate()
            .map(|(f, ids)| (Role::Test(f), ids))
            .chain([(Role::Train, &self.train), (Role::Unused, &self.unused)]);
        for (role, ids) in lists {
            ensure!(
                !ids.is_empty() || !matches!(role, Role::Test(_)),
                "{} is empty",
                role_words(role)
            );
            for id in ids {
                let Some(&i) = at.get(id.as_str()) else {
                    bail!(
                        "{} names genome {id}, which is not a sample of {dataset}",
                        role_words(role)
                    );
                };
                if let Some(first) = roles[i] {
                    bail!(
                        "genome {id} is listed twice, under {} and {}",
                        role_words(first),
                        role_words(role)
                    );
                }
                roles[i] = Some(role);
            }
        }
        let missing = roles.iter().filter(|r| r.is_none()).count();
        if let Some(i) = roles.iter().position(Option::is_none) {
            bail!(
                "the fold file leaves out {missing} of the {} samples of {dataset}, first {}; list every sample under test, train or unused",
                sample_ids.len(),
                sample_ids[i]
            );
        }
        Ok(roles.into_iter().flatten().collect())
    }
}

fn role_words(role: Role) -> String {
    match role {
        Role::Test(f) => format!("test fold {f}"),
        Role::Train => "train".into(),
        Role::Unused => "unused".into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Every fixture here is synthetic: invented genome ids and no real data.

    fn ids(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("g{i}")).collect()
    }

    fn file(test: &[&[&str]], train: &[&str], unused: &[&str]) -> FoldFile {
        let own = |ids: &[&str]| ids.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        FoldFile {
            format: FOLD_FILE_FORMAT.into(),
            dataset: "fixture__meropenem".into(),
            scheme: "first-public-past".into(),
            test: test.iter().map(|ids| own(ids)).collect(),
            train: own(train),
            unused: own(unused),
        }
    }

    fn refusal(f: &FoldFile, samples: &[String]) -> String {
        format!("{:#}", f.roles("fixture__meropenem", samples).unwrap_err())
    }

    #[test]
    fn roles_follow_the_sample_order_whatever_order_the_file_lists() {
        let f = file(&[&["g5", "g1"], &["g0"]], &["g4", "g2"], &["g3"]);
        assert_eq!(f.folds(), 2);
        assert_eq!(
            f.roles("fixture__meropenem", &ids(6)).unwrap(),
            [
                Role::Test(1),
                Role::Test(0),
                Role::Train,
                Role::Unused,
                Role::Train,
                Role::Test(0),
            ]
        );
    }

    #[test]
    fn the_hash_is_of_the_files_bytes_not_of_its_parsed_content() {
        let compact = br#"{"format":"grm-fold-file-v1","dataset":"d","scheme":"s","test":[["a"]],"train":["b"],"unused":[]}"#;
        let spaced = br#"{"format": "grm-fold-file-v1", "dataset": "d", "scheme": "s", "test": [["a"]], "train": ["b"], "unused": []}"#;
        let (a, b) = (parse(compact).unwrap(), parse(spaced).unwrap());
        assert_eq!(a.file, b.file);
        // `sha256sum` of each fixture line, taken outside this code.
        assert_eq!(
            a.sha256,
            "7cf08ed1ac3784c2543bcf6c9b323272a05bf82882b011bd22b0fa1afa43e47b"
        );
        assert_eq!(
            b.sha256,
            "a7dfc0517874110424cef06bd28fc50b3827164a83750856ba80d76a6e936400"
        );
    }

    #[test]
    fn a_file_is_refused_unless_it_lists_every_cohort_sample_exactly_once() {
        let samples = ids(6);
        let whole = file(&[&["g5", "g1"]], &["g0", "g2", "g4"], &["g3"]);
        assert!(whole.roles("fixture__meropenem", &samples).is_ok());

        let short = file(&[&["g5", "g1"]], &["g0", "g2"], &[]);
        assert_eq!(
            refusal(&short, &samples),
            "the fold file leaves out 2 of the 6 samples of fixture__meropenem, first g3; list every sample under test, train or unused"
        );
        let foreign = file(&[&["g5", "g9"]], &["g0", "g2", "g4", "g1"], &["g3"]);
        assert_eq!(
            refusal(&foreign, &samples),
            "test fold 0 names genome g9, which is not a sample of fixture__meropenem"
        );
        let twice = file(&[&["g5", "g1"]], &["g0", "g2", "g4", "g1"], &["g3"]);
        assert_eq!(
            refusal(&twice, &samples),
            "genome g1 is listed twice, under test fold 0 and train"
        );
        let repeated = file(&[&["g5", "g1"]], &["g0", "g2", "g4"], &["g3", "g3"]);
        assert_eq!(
            refusal(&repeated, &samples),
            "genome g3 is listed twice, under unused and unused"
        );
        let empty = file(&[&["g5", "g1"], &[]], &["g0", "g2", "g4"], &["g3"]);
        assert_eq!(refusal(&empty, &samples), "test fold 1 is empty");
    }

    #[test]
    fn a_file_for_another_dataset_or_with_a_drawn_split_name_is_refused() {
        let samples = ids(3);
        let good = file(&[&["g0"]], &["g1", "g2"], &[]);
        assert_eq!(
            format!(
                "{:#}",
                good.roles("fixture__ciprofloxacin", &samples).unwrap_err()
            ),
            "the fold file is for dataset \"fixture__meropenem\", not \"fixture__ciprofloxacin\""
        );
        for name in ["random", "lineage"] {
            let mut drawn = good.clone();
            drawn.scheme = name.into();
            assert!(
                format!("{:#}", drawn.check().unwrap_err())
                    .contains("is a split `grm learn` draws itself")
            );
        }
        for name in ["", "First-Public", "past only", &"x".repeat(49)] {
            let mut named = good.clone();
            named.scheme = name.to_string();
            assert!(
                format!("{:#}", named.check().unwrap_err()).contains("lowercase letters"),
                "{name:?}"
            );
        }
        let mut versioned = good.clone();
        versioned.format = "grm-fold-file-v2".into();
        assert_eq!(
            format!("{:#}", versioned.check().unwrap_err()),
            "format \"grm-fold-file-v2\" is not grm-fold-file-v1"
        );
        let mut foldless = good;
        foldless.test.clear();
        assert_eq!(
            format!("{:#}", foldless.check().unwrap_err()),
            "the file has no test fold"
        );
    }

    #[test]
    fn a_scheme_that_could_name_a_path_is_refused() {
        // The scheme becomes part of the result's file name
        // (`analysis::fold_result_path`), so nothing but a-z, 0-9 and the
        // hyphen may pass: no separator, no dot, no underscore, no NUL.
        let good = file(&[&["g0"]], &["g1", "g2"], &[]);
        for name in [
            "../escape",
            "a/b",
            "a.b",
            "a\\b",
            "past_only",
            "a\0b",
            "\u{e9}",
        ] {
            let mut named = good.clone();
            named.scheme = name.to_string();
            assert!(named.check().is_err(), "{name:?}");
        }
    }

    #[test]
    fn a_file_with_an_unknown_or_missing_key_is_refused() {
        let extra = br#"{"format":"grm-fold-file-v1","dataset":"d","scheme":"s","test":[["a"]],"train":["b"],"unused":[],"seed":7}"#;
        assert!(format!("{:#}", parse(extra).unwrap_err()).contains("unknown field `seed`"));
        let no_unused =
            br#"{"format":"grm-fold-file-v1","dataset":"d","scheme":"s","test":[["a"]],"train":["b"]}"#;
        assert!(format!("{:#}", parse(no_unused).unwrap_err()).contains("missing field `unused`"));
    }

    #[test]
    fn load_names_the_file_it_could_not_read_or_parse() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("folds.json");
        assert!(format!("{:#}", load(&path).unwrap_err()).starts_with("read "));
        std::fs::write(&path, b"{").unwrap();
        let error = format!("{:#}", load(&path).unwrap_err());
        assert!(
            error.starts_with(&format!("fold file {}: parse", path.display())),
            "{error}"
        );
        let good = br#"{"format":"grm-fold-file-v1","dataset":"d","scheme":"s","test":[["a"]],"train":["b"],"unused":[]}"#;
        std::fs::write(&path, good).unwrap();
        assert_eq!(load(&path).unwrap(), parse(good).unwrap());
    }
}
