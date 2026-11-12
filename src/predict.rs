//! `grm predict`: a published dataset report's final models, frozen, applied
//! to genomes. docs/prediction-contract.md is the contract this implements.
//!
//! A model feature is a dataset row: every k-mer with one presence pattern
//! over the cohort. The report lists one of them per feature (`kmer`), and
//! on a genome outside the cohort the feature's value is whether that k-mer
//! is found, on either strand. When the species matrix is given, each row
//! also counts how many of the feature's k-mers the genome carries, so a
//! feature found only in part is counted, not hidden.
//!
//! Only learners whose published form is a complete model are applied: the
//! SCM and CART always, the lasso when the report lists the feature of every
//! nonzero coefficient. Forests and boosted models are published as
//! importances only. Nothing is fitted, and calls use the cutoff the
//! evaluation scores with (`CALL_CUTOFF`).

use crate::analysis::{REPORT_SCHEMA, base64};
use crate::bits::BitMatrix;
use crate::dna::{KmerSpace, MAX_K, MIN_K};
use crate::fasta::{self, SeqStats};
use crate::kset::{self, KsetReader};
use crate::learn::lasso::LassoModel;
use crate::learn::scm::ScmModel;
use crate::learn::tree::Tree;
use crate::matrix::{Dataset, FeatureMatrix};
use crate::research::sha256;
use anyhow::{Context, Result, bail, ensure};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io::BufRead;
use std::path::{Path, PathBuf};

/// Names this contract in every prediction row.
pub const CONTRACT: &str = "grm-predict-1";

/// A score at or above this is a resistant call: the cutoff every learner's
/// metrics use (`binary_metrics(.., 0.5, ..)` in `learn::eval`).
pub const CALL_CUTOFF: f64 = 0.5;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Learner {
    Scm,
    Cart,
    Lasso,
}

impl Learner {
    pub const ALL: [Learner; 3] = [Learner::Scm, Learner::Cart, Learner::Lasso];

    pub fn name(self) -> &'static str {
        match self {
            Learner::Scm => "scm",
            Learner::Cart => "cart",
            Learner::Lasso => "lasso",
        }
    }
}

/// The parts of a dataset report (`analysis::DatasetReport`) a prediction
/// reads, its schema checked first. Other fields are ignored.
#[derive(Deserialize)]
struct Published {
    evaluation: PublishedEvaluation,
    input_fingerprint: String,
    id: String,
    counts: PublishedCounts,
    samples: Vec<PublishedSample>,
    models: PublishedModels,
    features: Vec<PublishedFeature>,
}

#[derive(Deserialize)]
struct PublishedEvaluation {
    revision: u32,
    valid: bool,
}

#[derive(Deserialize)]
struct PublishedCounts {
    features: usize,
    patterns: usize,
    k: usize,
    min_support: u32,
}

#[derive(Deserialize)]
struct PublishedSample {
    id: String,
    resistant: bool,
}

#[derive(Deserialize)]
struct PublishedModels {
    scm: ScmModel,
    cart: Tree,
    cart_depth: u32,
    lasso: Option<LassoModel>,
}

#[derive(Deserialize)]
struct PublishedFeature {
    feature: u32,
    n_kmers: usize,
    kmer: String,
    bits: String,
}

/// One feature a complete learner uses.
#[derive(Clone, Debug)]
pub struct Feature {
    /// Row of the report's dataset matrix, as the report numbers it.
    pub id: u32,
    pub kmer: String,
    /// `kmer` canonical and mixed, as ksets store k-mers.
    value: u64,
    pub n_kmers: usize,
    /// Every k-mer of the feature, mixed and sorted, once a matrix is given.
    pattern: Option<Vec<u64>>,
}

/// A report's final models with their features resolved, numbered by
/// position in `features`.
pub struct Model {
    pub dataset: String,
    pub report_sha256: String,
    pub input_fingerprint: String,
    pub evaluation_revision: u32,
    pub space: KmerSpace,
    pub features: Vec<Feature>,
    scm: ScmModel,
    cart: Tree,
    cart_depth: u32,
    lasso: Option<LassoModel>,
    /// Each learner's features as positions in `features`.
    uses: Vec<(Learner, Vec<usize>)>,
    /// Learners whose published form is not a complete model, and why.
    pub incomplete: Vec<(Learner, String)>,
    published: Published,
}

/// One genome's distinct canonical k-mers, mixed and sorted as in a kset.
pub struct Genome {
    pub name: String,
    pub source: PathBuf,
    /// Assembly statistics, for a genome counted from FASTA here.
    pub assembly: Option<SeqStats>,
    values: Vec<u64>,
}

#[derive(Debug, Serialize, PartialEq)]
pub struct FeatureValue {
    pub feature: u32,
    pub kmer: String,
    /// The value the model reads: the published k-mer is found.
    pub found: bool,
    pub n_kmers: usize,
    /// How many of the feature's k-mers the genome carries; None without a
    /// matrix.
    pub kmers_found: Option<usize>,
}

#[derive(Debug, Serialize)]
pub struct Prediction {
    pub contract: &'static str,
    pub dataset: String,
    pub report_sha256: String,
    pub input_fingerprint: String,
    pub evaluation_revision: u32,
    pub k: usize,
    pub genome: String,
    pub source: String,
    pub assembly: Option<SeqStats>,
    pub learner: Learner,
    /// "resistant" or "susceptible"; None for a null fit, which calls none.
    pub call: Option<&'static str>,
    /// CART's leaf value or the lasso's probability; the SCM has no score.
    pub score: Option<f64>,
    pub features: Vec<FeatureValue>,
    /// Features found in part: some of their k-mers but not all, where the
    /// published k-mer alone, any k-mer and every k-mer can give different
    /// values. None without a matrix.
    pub partial_features: Option<usize>,
}

impl Model {
    /// Reads a dataset report. Refuses a schema other than this build's, an
    /// evaluation the report itself marks invalid, and a k the k-mer space
    /// cannot hold. Learners that are not complete are recorded in
    /// `incomplete`, not refused here.
    pub fn load(path: &Path) -> Result<Model> {
        let bytes = std::fs::read(path).with_context(|| format!("read {}", path.display()))?;
        let schema: serde_json::Value = serde_json::from_slice(&bytes)
            .with_context(|| format!("{} is not JSON", path.display()))?;
        ensure!(
            schema["schema"].as_u64() == Some(REPORT_SCHEMA as u64),
            "{}: report schema {} is not {REPORT_SCHEMA}, the one grm predict reads",
            path.display(),
            schema["schema"]
        );
        let published: Published = serde_json::from_value(schema)
            .with_context(|| format!("{}: not a dataset report", path.display()))?;
        ensure!(
            published.evaluation.valid,
            "{}: the report's own evaluation audit is not valid",
            path.display()
        );
        let k = published.counts.k;
        ensure!(
            (MIN_K..=MAX_K).contains(&k),
            "{}: k {k} is outside {MIN_K}..={MAX_K}",
            path.display()
        );
        Model::resolve(published, sha256(&bytes))
    }

    fn resolve(published: Published, report_sha256: String) -> Result<Model> {
        let space = KmerSpace::new(published.counts.k);
        let mut listed: HashMap<u32, &PublishedFeature> = HashMap::new();
        for f in &published.features {
            ensure!(
                listed.insert(f.feature, f).is_none(),
                "the report lists feature {} twice",
                f.feature
            );
        }
        let m = &published.models;
        let wanted: [(Learner, Option<Vec<u32>>); 3] = [
            (
                Learner::Scm,
                Some(m.scm.rules.iter().map(|r| r.feature).collect()),
            ),
            (Learner::Cart, Some(cart_features(&m.cart, m.cart_depth))),
            (
                Learner::Lasso,
                m.lasso
                    .as_ref()
                    .map(|l| l.coef.iter().map(|c| c.0).collect()),
            ),
        ];
        let mut features: Vec<Feature> = Vec::new();
        let mut position: HashMap<u32, usize> = HashMap::new();
        let (mut uses, mut incomplete) = (Vec::new(), Vec::new());
        for (learner, ids) in wanted {
            let Some(ids) = ids else {
                incomplete.push((learner, "the report has no final model".to_string()));
                continue;
            };
            let unlisted = ids.iter().filter(|f| !listed.contains_key(f)).count();
            let unreadable = ids
                .iter()
                .filter_map(|f| listed.get(f))
                .filter(|f| space.canonical(&f.kmer).is_none())
                .count();
            if unlisted + unreadable > 0 {
                incomplete.push((
                    learner,
                    format!(
                        "{} of its {} features are not listed in the report and {unreadable} have no readable {}-mer",
                        unlisted,
                        ids.len(),
                        space.k
                    ),
                ));
                continue;
            }
            let mut rows = Vec::with_capacity(ids.len());
            for id in ids {
                let next = features.len();
                let at = *position.entry(id).or_insert(next);
                if at == next {
                    let f = listed[&id];
                    let canonical = space.canonical(&f.kmer).unwrap();
                    features.push(Feature {
                        id,
                        kmer: f.kmer.clone(),
                        value: space.mix(canonical),
                        n_kmers: f.n_kmers,
                        pattern: None,
                    });
                }
                rows.push(at);
            }
            uses.push((learner, rows));
        }
        let renumber = |f: &mut u32| *f = position.get(f).map_or(u32::MAX, |&p| p as u32);
        let mut scm = published.models.scm.clone();
        scm.rules.iter_mut().for_each(|r| renumber(&mut r.feature));
        let mut cart = published.models.cart.clone();
        for node in &mut cart.nodes {
            if let Some(s) = node.split.as_mut() {
                renumber(&mut s.feature);
            }
        }
        let mut lasso = published.models.lasso.clone();
        if let Some(l) = lasso.as_mut() {
            l.coef.iter_mut().for_each(|c| renumber(&mut c.0));
        }
        Ok(Model {
            dataset: published.id.clone(),
            report_sha256,
            input_fingerprint: published.input_fingerprint.clone(),
            evaluation_revision: published.evaluation.revision,
            space,
            features,
            cart_depth: published.models.cart_depth,
            scm,
            cart,
            lasso,
            uses,
            incomplete,
            published,
        })
    }

    /// The learners to apply: `requested`, or every complete one when it is
    /// empty. Refuses a requested learner that is not complete.
    pub fn learners(&self, requested: &[Learner]) -> Result<Vec<Learner>> {
        for l in requested {
            if let Some((_, why)) = self.incomplete.iter().find(|(i, _)| i == l) {
                bail!(
                    "{}: the published {} is not a complete model: {why}",
                    self.dataset,
                    l.name()
                );
            }
        }
        let chosen: Vec<Learner> = if requested.is_empty() {
            self.uses.iter().map(|(l, _)| *l).collect()
        } else {
            Learner::ALL
                .into_iter()
                .filter(|l| requested.contains(l))
                .collect()
        };
        ensure!(
            !chosen.is_empty(),
            "{}: no learner is complete",
            self.dataset
        );
        Ok(chosen)
    }

    /// Reads every k-mer of each model feature from the species matrix the
    /// report was learned from. The matrix must rebuild the report's dataset
    /// exactly: same k, every sample a column, the same feature and pattern
    /// counts, the same presence bits for every listed feature, and each
    /// published k-mer among its feature's k-mers.
    pub fn attach_matrix(&mut self, dir: &Path) -> Result<()> {
        let matrix = FeatureMatrix::load(dir)?;
        let p = &self.published;
        ensure!(
            matrix.summary.k == self.space.k,
            "the matrix in {} has k {}, the report k {}",
            dir.display(),
            matrix.summary.k,
            self.space.k
        );
        let column: HashMap<&str, usize> = matrix
            .summary
            .genomes
            .iter()
            .enumerate()
            .map(|(i, g)| (g.as_str(), i))
            .collect();
        let columns = p
            .samples
            .iter()
            .map(|s| {
                column.get(s.id.as_str()).copied().with_context(|| {
                    format!("sample {} of the report is not in {}", s.id, dir.display())
                })
            })
            .collect::<Result<Vec<usize>>>()?;
        let resistant: Vec<bool> = p.samples.iter().map(|s| s.resistant).collect();
        let ds = Dataset::project(&matrix, columns, &resistant, p.counts.min_support);
        ensure!(
            ds.n_features() == p.counts.features && ds.feature_patterns.len() == p.counts.patterns,
            "the matrix in {} gives {} features over {} patterns, the report {} over {}",
            dir.display(),
            ds.n_features(),
            ds.feature_patterns.len(),
            p.counts.features,
            p.counts.patterns
        );
        let bits: HashMap<u32, &str> = p
            .features
            .iter()
            .map(|f| (f.feature, f.bits.as_str()))
            .collect();
        for f in &mut self.features {
            let row = f.id as usize;
            ensure!(
                row < ds.n_features() && base64(ds.x.row(row)) == bits[&f.id],
                "feature {}: the matrix in {} gives other presence bits than the report",
                f.id,
                dir.display()
            );
            let mut pattern: Vec<u64> = ds
                .patterns_of(row)
                .iter()
                .flat_map(|&q| matrix.kmers_of(q as usize).iter().copied())
                .collect();
            pattern.sort_unstable();
            ensure!(
                pattern.len() == f.n_kmers && pattern.binary_search(&f.value).is_ok(),
                "feature {}: the matrix in {} holds {} k-mers without {}, the report {}",
                f.id,
                dir.display(),
                pattern.len(),
                f.kmer,
                f.n_kmers
            );
            f.pattern = Some(pattern);
        }
        Ok(())
    }

    /// Counts an assembly at the report's k, as `grm count` does. Refuses
    /// reads (FASTQ) and a file with no k-mer.
    pub fn genome_from_fasta(&self, path: &Path) -> Result<Genome> {
        ensure!(
            !is_fastq(path)?,
            "{} holds reads (FASTQ); only an assembly is scorable",
            path.display()
        );
        let (values, stats) = kset::count_file(path, &self.space, 1)?;
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .map(|n| n.strip_suffix(".gz").unwrap_or(n))
            .map(|n| n.rsplit_once('.').map_or(n, |(stem, _)| stem))
            .unwrap_or_default()
            .to_string();
        self.genome(name, path, Some(stats), values)
    }

    /// Reads a kset that `grm count` wrote. Refuses one counted at another k.
    pub fn genome_from_kset(&self, path: &Path) -> Result<Genome> {
        let reader = KsetReader::open(path)?;
        ensure!(
            reader.k == self.space.k,
            "{} was counted at k {}, the report at k {}",
            path.display(),
            reader.k,
            self.space.k
        );
        let values = reader.read_all()?;
        self.genome(reader.name.clone(), path, None, values)
    }

    fn genome(
        &self,
        name: String,
        path: &Path,
        assembly: Option<SeqStats>,
        values: Vec<u64>,
    ) -> Result<Genome> {
        ensure!(
            !values.is_empty(),
            "{} holds no {}-mer",
            path.display(),
            self.space.k
        );
        Ok(Genome {
            name,
            source: path.to_path_buf(),
            assembly,
            values,
        })
    }

    /// One row per learner in `learners` for `genome`.
    pub fn predict(&self, genome: &Genome, learners: &[Learner]) -> Vec<Prediction> {
        let found: Vec<bool> = self
            .features
            .iter()
            .map(|f| genome.values.binary_search(&f.value).is_ok())
            .collect();
        let carried: Vec<Option<usize>> = self
            .features
            .iter()
            .map(|f| {
                f.pattern.as_ref().map(|p| {
                    p.iter()
                        .filter(|v| genome.values.binary_search(v).is_ok())
                        .count()
                })
            })
            .collect();
        let mut x = BitMatrix::with_cols(1);
        for &f in &found {
            x.push_row(&[f as u64]);
        }
        learners
            .iter()
            .map(|&learner| {
                let rows = &self.uses.iter().find(|(l, _)| *l == learner).unwrap().1;
                let score = match learner {
                    Learner::Scm => None,
                    Learner::Cart => Some(self.cart.predict(&x, self.cart_depth)[0]),
                    Learner::Lasso => self
                        .lasso
                        .as_ref()
                        .filter(|l| !l.coef.is_empty())
                        .map(|l| l.predict_proba(&x)[0]),
                };
                let resistant = match learner {
                    Learner::Scm => Some(self.scm.predict(&x).contains(0)),
                    _ => score.map(|s| s >= CALL_CUTOFF),
                };
                let features: Vec<FeatureValue> = rows
                    .iter()
                    .map(|&r| FeatureValue {
                        feature: self.features[r].id,
                        kmer: self.features[r].kmer.clone(),
                        found: found[r],
                        n_kmers: self.features[r].n_kmers,
                        kmers_found: carried[r],
                    })
                    .collect();
                let partial_features = features
                    .iter()
                    .map(|f| f.kmers_found.map(|c| (0 < c && c < f.n_kmers) as usize))
                    .sum();
                Prediction {
                    contract: CONTRACT,
                    dataset: self.dataset.clone(),
                    report_sha256: self.report_sha256.clone(),
                    input_fingerprint: self.input_fingerprint.clone(),
                    evaluation_revision: self.evaluation_revision,
                    k: self.space.k,
                    genome: genome.name.clone(),
                    source: genome.source.display().to_string(),
                    assembly: genome.assembly.clone(),
                    learner,
                    call: resistant.map(|r| if r { "resistant" } else { "susceptible" }),
                    score,
                    features,
                    partial_features,
                }
            })
            .collect()
    }
}

/// Split features of the nodes a prediction at `depth` can reach, in node
/// order and without repeats (`Tree::predict` treats nodes at `depth` as
/// leaves).
fn cart_features(tree: &Tree, depth: u32) -> Vec<u32> {
    let mut out: Vec<u32> = Vec::new();
    for node in &tree.nodes {
        if let Some(s) = node.split.as_ref().filter(|_| node.depth < depth)
            && !out.contains(&s.feature)
        {
            out.push(s.feature);
        }
    }
    out
}

/// Whether the file's first record is a FASTQ record.
fn is_fastq(path: &Path) -> Result<bool> {
    let mut reader = fasta::open(path)?;
    let first = reader
        .fill_buf()
        .with_context(|| format!("read {}", path.display()))?;
    Ok(first.first() == Some(&b'@'))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;
    use serde_json::{Value, json};

    const K: usize = 15;

    fn random(rng: &mut ChaCha8Rng, n: usize) -> String {
        (0..n)
            .map(|_| b"ACGT"[rng.random_range(0..4)] as char)
            .collect()
    }

    fn reverse_complement(seq: &str) -> String {
        seq.bytes()
            .rev()
            .map(|b| match b {
                b'A' => 'T',
                b'C' => 'G',
                b'G' => 'C',
                b'T' => 'A',
                other => other as char,
            })
            .collect()
    }

    struct Fixture {
        dir: tempfile::TempDir,
        /// The block g0 and g1 carry, and the model's only feature.
        block: String,
        /// The block's published k-mer starts here.
        at: usize,
        core: String,
        report: Value,
    }

    impl Fixture {
        fn path(&self, name: &str) -> PathBuf {
            self.dir.path().join(name)
        }

        fn write_report(&self, report: &Value) -> PathBuf {
            let path = self.path("report.json");
            std::fs::write(&path, serde_json::to_vec(report).unwrap()).unwrap();
            path
        }

        fn model(&self) -> Model {
            Model::load(&self.write_report(&self.report)).unwrap()
        }

        fn fasta(&self, name: &str, seq: &str) -> PathBuf {
            let path = self.path(&format!("{name}.fna"));
            std::fs::write(&path, format!(">c\n{seq}\n")).unwrap();
            path
        }

        fn rows(&self, model: &Model, seq: &str) -> Vec<Prediction> {
            let genome = model.genome_from_fasta(&self.fasta("new", seq)).unwrap();
            model.predict(&genome, &model.learners(&[]).unwrap())
        }

        fn published_kmer(&self) -> &str {
            &self.block[self.at..self.at + K]
        }
    }

    /// Four genomes on a shared core: g0 and g1 (resistant) carry a 60-base
    /// block, g2 and g3 another. The species matrix is saved under
    /// `matrix`, and the report's three learners all call resistant on the
    /// first block alone, its k-mer published as read at base 20.
    fn fixture() -> Fixture {
        let dir = tempfile::tempdir().unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(3);
        let (core, block, other) = (
            random(&mut rng, 300),
            random(&mut rng, 60),
            random(&mut rng, 60),
        );
        let space = KmerSpace::new(K);
        let ksets: Vec<PathBuf> = [
            ("g0", format!("{core}N{block}")),
            ("g1", format!("{block}N{core}")),
            ("g2", format!("{core}N{other}")),
            ("g3", format!("{other}N{core}")),
        ]
        .iter()
        .map(|(id, seq)| {
            let fna = dir.path().join(format!("{id}.fna"));
            std::fs::write(&fna, format!(">c\n{seq}\n")).unwrap();
            let (values, _) = kset::count_file(&fna, &space, 1).unwrap();
            let path = dir.path().join(format!("{id}.kset"));
            kset::write(&path, id, &space, &values).unwrap();
            path
        })
        .collect();
        let m = matrix::build(&ksets, 1, 8).unwrap();
        m.save(&dir.path().join("matrix")).unwrap();
        let ds = Dataset::project(&m, vec![0, 1, 2, 3], &[true, true, false, false], 1);
        let at = 20;
        let value = space.mix(space.canonical(&block[at..at + K]).unwrap());
        let f = (0..ds.n_features())
            .find(|&f| {
                ds.patterns_of(f)
                    .iter()
                    .any(|&p| m.kmers_of(p as usize).contains(&value))
            })
            .unwrap();
        let n_kmers: usize = ds
            .patterns_of(f)
            .iter()
            .map(|&p| m.kmers_of(p as usize).len())
            .sum();
        assert_eq!(n_kmers, 60 - K + 1);
        let report = json!({
            "schema": REPORT_SCHEMA,
            "evaluation": {"revision": 12, "valid": true, "reasons": []},
            "input_fingerprint": "0123456789abcdef",
            "id": "fixture__drug",
            "counts": {
                "samples": 4,
                "features": ds.n_features(),
                "patterns": ds.feature_patterns.len(),
                "k": K,
                "min_support": 1,
            },
            "samples": [
                {"id": "g0", "resistant": true},
                {"id": "g1", "resistant": true},
                {"id": "g2", "resistant": false},
                {"id": "g3", "resistant": false},
            ],
            "models": {
                "scm": {
                    "kind": "conjunction",
                    "p": 1.0,
                    "rules": [{"feature": f, "presence": true}],
                    "steps": [],
                },
                "cart": {"nodes": [
                    {"value": 0.5, "n": 4, "n_resistant": 2, "depth": 0,
                     "split": {"feature": f, "gain": 0.5, "absent": 1, "present": 2}},
                    {"value": 0.0, "n": 2, "n_resistant": 0, "depth": 1, "split": null},
                    {"value": 1.0, "n": 2, "n_resistant": 2, "depth": 1, "split": null},
                ]},
                "cart_depth": 1,
                "lasso": {"intercept": -1.0, "coef": [[f, 2.5]], "lambda": 0.1},
            },
            "features": [{
                "feature": f,
                "n_kmers": n_kmers,
                "kmer": &block[at..at + K],
                "bits": base64(ds.x.row(f)),
            }],
        });
        Fixture {
            dir,
            block,
            at,
            core,
            report,
        }
    }

    fn calls(rows: &[Prediction]) -> Vec<Option<&'static str>> {
        rows.iter().map(|r| r.call).collect()
    }

    const RESISTANT: [Option<&str>; 3] = [Some("resistant"); 3];
    const SUSCEPTIBLE: [Option<&str>; 3] = [Some("susceptible"); 3];

    #[test]
    fn the_cohort_genomes_get_their_own_bits_from_fasta_and_kset_alike() {
        let fx = fixture();
        let model = fx.model();
        assert_eq!(model.learners(&[]).unwrap(), Learner::ALL);
        for (id, want) in [
            ("g0", RESISTANT),
            ("g1", RESISTANT),
            ("g2", SUSCEPTIBLE),
            ("g3", SUSCEPTIBLE),
        ] {
            let from_fasta = model.genome_from_fasta(&fx.path(&format!("{id}.fna")));
            let from_kset = model.genome_from_kset(&fx.path(&format!("{id}.kset")));
            for genome in [from_fasta.unwrap(), from_kset.unwrap()] {
                assert_eq!(genome.name, id);
                let rows = model.predict(&genome, &Learner::ALL);
                assert_eq!(calls(&rows), want, "{id}");
                assert_eq!(
                    rows[1].score,
                    Some(if want == RESISTANT { 1.0 } else { 0.0 })
                );
                assert_eq!(rows[0].score, None);
            }
        }
    }

    #[test]
    fn a_row_names_the_model_and_the_genome() {
        let fx = fixture();
        let path = fx.write_report(&fx.report);
        let model = Model::load(&path).unwrap();
        let rows = fx.rows(&model, &fx.block);
        let row = serde_json::to_value(&rows[2]).unwrap();
        assert_eq!(row["contract"], CONTRACT);
        assert_eq!(row["dataset"], "fixture__drug");
        assert_eq!(
            row["report_sha256"],
            sha256(&std::fs::read(&path).unwrap()).as_str()
        );
        assert_eq!(row["input_fingerprint"], "0123456789abcdef");
        assert_eq!(row["evaluation_revision"], 12);
        assert_eq!(row["k"], K);
        assert_eq!(row["genome"], "new");
        assert_eq!(row["learner"], "lasso");
        assert_eq!(row["call"], "resistant");
        assert_eq!(row["assembly"]["length"], 60);
        let p = 1.0 / (1.0 + (-1.5f64).exp());
        assert!((row["score"].as_f64().unwrap() - p).abs() < 1e-12, "{row}");
        assert_eq!(row["features"][0]["kmer"], fx.published_kmer());
        assert_eq!(row["features"][0]["found"], true);
        assert_eq!(row["features"][0]["kmers_found"], Value::Null);
        assert_eq!(row["partial_features"], Value::Null);
    }

    #[test]
    fn a_feature_found_in_part_is_counted_and_its_value_follows_the_published_kmer() {
        let fx = fixture();
        let mut model = fx.model();
        model.attach_matrix(&fx.path("matrix")).unwrap();
        let n = 60 - K + 1;
        // Whole, missing, the first 40 bases (with the published k-mer at
        // 20..35), and the last 30 bases (without it).
        for (seq, found, carried) in [
            (fx.block.clone(), true, n),
            (fx.core.clone(), false, 0),
            (fx.block[..40].to_string(), true, 40 - K + 1),
            (fx.block[30..].to_string(), false, 30 - K + 1),
        ] {
            let rows = fx.rows(&model, &seq);
            assert_eq!(calls(&rows), if found { RESISTANT } else { SUSCEPTIBLE });
            for row in &rows {
                let f = &row.features[0];
                assert_eq!(
                    (f.found, f.n_kmers, f.kmers_found),
                    (found, n, Some(carried))
                );
                let partial = (0 < carried && carried < n) as usize;
                assert_eq!(row.partial_features, Some(partial));
            }
        }
    }

    #[test]
    fn the_published_kmer_is_found_on_either_strand_and_not_across_a_non_acgt_byte() {
        let fx = fixture();
        let model = fx.model();
        let found = |seq: &str| fx.rows(&model, seq)[0].features[0].found;
        assert!(found(&fx.block));
        assert!(found(&reverse_complement(&fx.block)));
        assert!(found(&fx.block.to_lowercase()));
        // An N, or an IUPAC code, inside the published k-mer's window.
        for byte in ["N", "R"] {
            let cut = fx.at + K / 2;
            let broken = format!("{}{byte}{}", &fx.block[..cut], &fx.block[cut + 1..]);
            assert!(!found(&broken), "{byte}");
            assert!(!found(&reverse_complement(&broken)), "{byte}");
        }
    }

    #[test]
    fn a_report_of_another_schema_or_an_invalid_evaluation_is_refused() {
        let fx = fixture();
        let refusal = |edit: &dyn Fn(&mut Value)| {
            let mut report = fx.report.clone();
            edit(&mut report);
            let path = fx.write_report(&report);
            format!("{:#}", Model::load(&path).err().unwrap())
        };
        let err = refusal(&|r| r["schema"] = json!(1));
        assert!(err.contains("report schema 1 is not 2"), "{err}");
        let err = refusal(&|r| {
            r.as_object_mut().unwrap().remove("schema");
        });
        assert!(err.contains("report schema null"), "{err}");
        let err = refusal(&|r| r["evaluation"]["valid"] = json!(false));
        assert!(err.contains("evaluation audit is not valid"), "{err}");
        let err = refusal(&|r| r["counts"]["k"] = json!(33));
        assert!(err.contains("k 33 is outside"), "{err}");
    }

    #[test]
    fn a_learner_with_an_unresolved_feature_is_left_out_and_refused_when_named() {
        let fx = fixture();
        let f = fx.report["features"][0]["feature"].as_u64().unwrap();
        let load = |edit: &dyn Fn(&mut Value)| {
            let mut report = fx.report.clone();
            edit(&mut report);
            Model::load(&fx.write_report(&report)).unwrap()
        };
        // A lasso coefficient on a feature the report does not list.
        let model = load(&|r| r["models"]["lasso"]["coef"] = json!([[f, 2.5], [f + 1, 0.3]]));
        assert_eq!(model.learners(&[]).unwrap(), [Learner::Scm, Learner::Cart]);
        let err = model.learners(&[Learner::Lasso]).unwrap_err().to_string();
        assert!(
            err.contains("published lasso is not a complete model: 1 of its 2"),
            "{err}"
        );
        // No final lasso at all.
        let model = load(&|r| r["models"]["lasso"] = Value::Null);
        assert!(model.learners(&[Learner::Lasso]).is_err());
        // A published k-mer that is not k ACGT letters leaves out every
        // learner that reads it, and nothing is left to apply.
        for kmer in [&fx.published_kmer()[1..], "ACGTNACGTACGTAC"] {
            let model = load(&|r| r["features"][0]["kmer"] = json!(kmer));
            assert_eq!(model.incomplete.len(), 3, "{kmer}");
            let err = model.learners(&[]).unwrap_err().to_string();
            assert!(err.contains("no learner is complete"), "{err}");
        }
        // A CART split below the published depth is never read.
        let model = load(&|r| r["models"]["cart_depth"] = json!(0));
        assert_eq!(model.learners(&[Learner::Cart]).unwrap(), [Learner::Cart]);
        assert_eq!(
            calls(&fx.rows(&model, &fx.block))[1],
            Some("resistant"),
            "a 0.5 root is a resistant call"
        );
    }

    #[test]
    fn a_lasso_without_coefficients_makes_no_call() {
        let fx = fixture();
        let mut report = fx.report.clone();
        report["models"]["lasso"]["coef"] = json!([]);
        let model = Model::load(&fx.write_report(&report)).unwrap();
        let rows = fx.rows(&model, &fx.block);
        assert_eq!((rows[2].call, rows[2].score), (None, None));
        assert!(rows[2].features.is_empty());
    }

    #[test]
    fn a_genome_of_another_k_reads_or_no_kmer_is_refused() {
        let fx = fixture();
        let model = fx.model();
        let other = KmerSpace::new(K + 2);
        let (values, _) = kset::count_file(&fx.path("g0.fna"), &other, 1).unwrap();
        let path = fx.path("g0.k17.kset");
        kset::write(&path, "g0", &other, &values).unwrap();
        let err = model.genome_from_kset(&path).err().unwrap().to_string();
        assert!(err.contains("counted at k 17, the report at k 15"), "{err}");

        let reads = fx.path("reads.fq");
        std::fs::write(
            &reads,
            format!("@r1\n{}\n+\n{}\n", fx.block, "I".repeat(60)),
        )
        .unwrap();
        let err = model.genome_from_fasta(&reads).err().unwrap().to_string();
        assert!(err.contains("holds reads (FASTQ)"), "{err}");

        let short = fx.fasta("short", &fx.block[..K - 1]);
        let err = model.genome_from_fasta(&short).err().unwrap().to_string();
        assert!(err.contains("holds no 15-mer"), "{err}");
    }

    #[test]
    fn a_matrix_that_does_not_rebuild_the_report_is_refused() {
        let fx = fixture();
        let attach = |edit: &dyn Fn(&mut Value)| {
            let mut report = fx.report.clone();
            edit(&mut report);
            let mut model = Model::load(&fx.write_report(&report)).unwrap();
            model
                .attach_matrix(&fx.path("matrix"))
                .err()
                .map(|e| format!("{e:#}"))
        };
        assert_eq!(attach(&|_| {}), None);
        let err = attach(&|r| r["samples"][3]["id"] = json!("g9")).unwrap();
        assert!(err.contains("sample g9 of the report is not in"), "{err}");
        let err = attach(&|r| r["counts"]["features"] = json!(99)).unwrap();
        assert!(err.contains("the report 99 over"), "{err}");
        // g1 and g2 swap labels and ids: the same columns, other bits.
        let err = attach(&|r| {
            r["samples"][1]["id"] = json!("g2");
            r["samples"][2]["id"] = json!("g1");
        })
        .unwrap();
        assert!(err.contains("other presence bits than the report"), "{err}");
        let err = attach(&|r| r["features"][0]["kmer"] = json!(&fx.core[..K])).unwrap();
        assert!(err.contains("k-mers without"), "{err}");
        let err = attach(&|r| r["counts"]["k"] = json!(K + 2)).unwrap();
        assert!(err.contains("has k 15, the report k 17"), "{err}");
    }
}
