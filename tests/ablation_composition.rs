//! A toy composition of the ablation parts on a generic random bit matrix:
//! G's `fold_plan`; H's `KmerSource` over a toy `FeatureMatrix`, admitted;
//! I's `fit_encoder`, `transform`, `block_candidates`, `screen_block` and
//! `assemble` for each fit; boosting (`learn::gbt::fit`) on the assembled
//! columns; and one prediction row per outer test sample, in K's row shape
//! (contract §4, `docs/ablation-report.md`).
//!
//! This is a software composition test. It is not the arm runner, the input
//! optimizer or an experiment, and it measures nothing: one arm (K alone),
//! one learner, random labels and random bits. It checks that the public
//! parts compose in the order `docs/ablation-runner.md` gives, with every
//! mask taken from the plan; that no outer test sample reaches a fit; and
//! that the rows are the same at any thread count. What an arm runner owes
//! is in `docs/ablation-runner.md`.
use grm::ablate::blocks::{K_SCREEN_CAP, KmerSource, PredictorBlock, admit};
use grm::ablate::encode::{
    EncodeConfig, EncodedBlock, assemble, block_candidates, fit_encoder, screen_block, transform,
};
use grm::ablate::plan::{FoldPlan, GroupIds, Scheme, fold_plan, fold_seed};
use grm::bits::{BitMatrix, Bits};
use grm::learn::eval::EvalConfig;
use grm::learn::gbt::{self, GbtParams, Objective};
use grm::matrix::{FeatureMatrix, MatrixSummary};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use serde_json::{Value, json};
use std::collections::BTreeMap;

const SEED: u64 = 0x2026_1102_0005_0001;
/// Samples, in lineages of four.
const N: usize = 60;
/// More patterns than K's screen cap, so each fit's screen drops some.
const PATTERNS: usize = K_SCREEN_CAP + 200;
const MIN_SUPPORT: u32 = 2;
/// K's prediction row fields, sorted (`PREDICTION` in
/// `scripts/ablation/schema.py`).
const ROW_FIELDS: [&str; 12] = [
    "arm",
    "called",
    "cluster",
    "dataset",
    "fold",
    "fold_plan_sha256",
    "learner",
    "repeat",
    "sample_id",
    "score",
    "split",
    "y",
];
/// What `FoldRun::fitted` holds, in order.
const STAGES: [&str; 5] = [
    "encoder",
    "candidates",
    "screened block",
    "assembled arm",
    "model",
];

/// Generic synthetic samples: random labels and random pattern bits.
#[derive(Clone)]
struct Toy {
    samples: Vec<String>,
    y: Bits,
    groups: GroupIds,
    /// Each pattern's samples.
    rows: Vec<Bits>,
}

impl Toy {
    /// N samples in lineages of four, each lineage holding both classes.
    /// Every third lineage, from the second, holds a duplicate pair.
    fn new() -> Toy {
        let mut rng = ChaCha8Rng::seed_from_u64(SEED);
        let mut labels: Vec<bool> = (0..N).map(|_| rng.random_bool(0.4)).collect();
        for lineage in labels.chunks_mut(4) {
            if lineage.iter().all(|&l| l == lineage[0]) {
                let i = rng.random_range(0..4);
                lineage[i] = !lineage[i];
            }
        }
        let lineage = (0..N as u32).map(|i| i / 4).collect();
        let duplicate = (0..N / 4)
            .filter(|j| j % 3 == 1)
            .map(|j| vec![4 * j + 1, 4 * j + 2])
            .collect();
        let mut toy = Toy {
            samples: (0..N).map(|i| format!("s{i:02}")).collect(),
            y: Bits::from_indices(N, (0..N).filter(|&i| labels[i])),
            groups: GroupIds::new(lineage, duplicate, None).unwrap(),
            rows: vec![Bits::zeros(N); PATTERNS],
        };
        toy.redraw(&Bits::ones(N), &mut rng);
        toy
    }

    /// Redraws the pattern bits of the samples in `which`. Duplicates keep
    /// identical patterns, as k-mer-identical genomes do.
    fn redraw(&mut self, which: &Bits, rng: &mut ChaCha8Rng) {
        for row in &mut self.rows {
            for i in which.iter_ones() {
                if rng.random_bool(0.5) {
                    row.insert(i);
                } else {
                    row.remove(i);
                }
            }
            for group in self.groups.duplicate() {
                for &i in &group[1..] {
                    if row.contains(group[0]) {
                        row.insert(i);
                    } else {
                        row.remove(i);
                    }
                }
            }
        }
    }

    /// The toy's k-mer matrix, one k-mer per pattern.
    fn matrix(&self) -> FeatureMatrix {
        let mut patterns = BitMatrix::with_cols(N);
        for row in &self.rows {
            patterns.push_row(row.words());
        }
        let p = self.rows.len();
        FeatureMatrix {
            summary: MatrixSummary {
                k: 31,
                genomes: self.samples.clone(),
                min_support: MIN_SUPPORT,
                n_patterns: p,
                kmers_distinct: p as u64,
                kmers_kept: p as u64,
                kmers_rare: 0,
                kmers_core: 0,
                kset_digests: vec![],
            },
            patterns,
            support: self.rows.iter().map(Bits::count).collect(),
            kmer_offsets: (0..=p as u64).collect(),
            kmers: (0..p as u64).collect(),
        }
    }

    /// H's K block over every sample, admitted.
    fn block(&self) -> PredictorBlock {
        let matrix = self.matrix();
        admit(&KmerSource::new(&matrix).unwrap(), &self.samples).unwrap()
    }

    /// The groups `scheme` keeps whole: duplicate groups, or lineages.
    fn whole(&self, scheme: Scheme) -> Vec<Vec<usize>> {
        match scheme {
            Scheme::Random => self.groups.duplicate().to_vec(),
            Scheme::Lineage => {
                let mut members: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
                for (i, &l) in self.groups.lineage().iter().enumerate() {
                    members.entry(l).or_default().push(i);
                }
                members.into_values().collect()
            }
        }
    }
}

/// `y` with the labels of `which` flipped.
fn flip(y: &Bits, which: &Bits) -> Bits {
    y.and_not(which).or(&which.and_not(y))
}

/// A stage's provenance and its rows on the training samples alone.
fn training_view(block: &EncodedBlock, train: &Bits) -> Value {
    let rows: Vec<Vec<u64>> = (0..block.x.rows)
        .map(|r| train.and_slice(block.x.row(r)).words().to_vec())
        .collect();
    json!({"prov": block.prov, "rows": rows})
}

/// One outer fold: its masks, what its fit made, and its scores.
struct FoldRun {
    train: Bits,
    /// Boosting's tuning fit (`fit_m`) and its early-stopping holdout.
    fit: Bits,
    val: Bits,
    /// Candidate columns, then screened columns.
    widths: (usize, usize),
    /// Each of `STAGES`, as JSON. The encoded stages keep only their
    /// training rows, since every sample is encoded and scored.
    fitted: Vec<Value>,
    /// Every sample's score, or None when boosting kept no tree: the 0.5
    /// prior, a no-call.
    scores: Option<Vec<f64>>,
    /// The share of the outer training samples that is resistant: a no-call
    /// row's score.
    prevalence: f64,
}

/// Outer fold `f` of the K-only arm, in `docs/ablation-runner.md`'s order.
/// Boosting makes one fit, on the plan's gbt `fit_m`, stopping early on the
/// holdout, as `learn::eval` does; so the block is screened once, on `fit_m`.
fn run_fold(
    block: &PredictorBlock,
    samples: &[String],
    y: &Bits,
    plan: &FoldPlan,
    f: usize,
) -> FoldRun {
    let n = samples.len();
    let cfg = EncodeConfig {
        min_support: MIN_SUPPORT,
    };
    // 1. Every mask from the plan.
    let train = Bits::from_indices(n, (0..n).filter(|&i| plan.outer[i] != f));
    let val = Bits::from_indices(n, plan.holdouts[f].gbt.iter().copied());
    let fit = train.and_not(&val);
    // 2. The encoder, fitted on the outer training samples.
    let encoder = fit_encoder(block.raw(), samples, &train, &cfg).unwrap();
    let encoded = transform(&encoder, block.raw(), samples).unwrap();
    // 3. The candidates, on the outer training mask; no labels.
    let candidates = block_candidates(&encoded, &train, &cfg);
    // 4. The screen reads labels, so only the fit's own rows.
    let screened = screen_block(&candidates, y, &fit, block.spec().screen_cap);
    // 5. Assemble on the outer training mask, then the learner.
    let arm = assemble(std::slice::from_ref(&screened), &train);
    let features: Vec<u32> = (0..arm.x.rows as u32).collect();
    let params = GbtParams {
        seed: fold_seed(plan.repeat_seed, f),
        ..EvalConfig::default().gbt
    };
    let model = gbt::fit(
        &arm.x,
        &Objective::Logistic { y },
        &fit,
        &val,
        &features,
        &params,
    );
    let mut fitted_encoder = serde_json::to_value(&encoder).unwrap();
    // The spec hash covers the whole matrix, outer test columns included: it
    // names the input, and nothing is fitted from it.
    fitted_encoder
        .as_object_mut()
        .unwrap()
        .remove("spec_sha256");
    FoldRun {
        widths: (candidates.x.rows, screened.x.rows),
        fitted: vec![
            fitted_encoder,
            training_view(&candidates, &train),
            training_view(&screened, &train),
            training_view(&arm, &train),
            serde_json::to_value(&model).unwrap(),
        ],
        scores: (!model.trees.is_empty()).then(|| model.predict_proba(&arm.x)),
        prevalence: f64::from(train.and(y).count()) / f64::from(train.count()),
        train,
        fit,
        val,
    }
}

/// One split of the K-only arm: its plan, block, fold runs and rows.
struct Arm {
    plan: FoldPlan,
    block: PredictorBlock,
    runs: Vec<FoldRun>,
    rows: Vec<Value>,
}

fn arm(toy: &Toy, scheme: Scheme) -> Arm {
    let cfg = EvalConfig::default();
    let plan = fold_plan(
        &toy.y,
        &toy.groups,
        scheme,
        cfg.folds,
        cfg.inner_folds,
        cfg.seed,
        0,
    )
    .unwrap();
    assert!(
        plan.is_available(),
        "{scheme}: {:?}",
        plan.unavailable_reason
    );
    let block = toy.block();
    let runs: Vec<FoldRun> = (0..plan.k)
        .map(|f| run_fold(&block, &toy.samples, &toy.y, &plan, f))
        .collect();
    let rows = (0..N)
        .map(|i| {
            let run = &runs[plan.outer[i]];
            let (score, called) = match &run.scores {
                Some(scores) => (scores[i], scores[i] >= 0.5),
                None => (run.prevalence, false),
            };
            json!({
                "dataset": "toy",
                "arm": "K",
                "learner": "gbt",
                "split": scheme.as_str(),
                "repeat": plan.repeat,
                "sample_id": toy.samples[i],
                "y": u8::from(toy.y.contains(i)),
                "fold": plan.outer[i],
                "score": score,
                "called": called,
                "cluster": toy.groups.lineage()[i],
                "fold_plan_sha256": plan.sha256,
            })
        })
        .collect();
    Arm {
        plan,
        block,
        runs,
        rows,
    }
}

/// Step 1's masks. Fold f trains on the samples outside its outer test
/// fold, and boosting fits on those outside the plan's gbt holdout. No mask
/// holds an outer test sample, each holds both classes, and no group the
/// split keeps whole is cut. Each fit's screen keeps K's cap.
fn check_masks(toy: &Toy, arm: &Arm) {
    let plan = &arm.plan;
    let mut tested = Bits::zeros(N);
    for (f, run) in arm.runs.iter().enumerate() {
        let at = format!("{} fold {f}", plan.scheme);
        let test = Bits::from_indices(N, (0..N).filter(|&i| plan.outer[i] == f));
        assert_eq!(run.train, test.not(), "{at}: training samples");
        assert_eq!(
            run.val.iter_ones().collect::<Vec<_>>(),
            plan.holdouts[f].gbt,
            "{at}: holdout"
        );
        assert_eq!(run.fit.or(&run.val), run.train, "{at}: fit and holdout");
        assert_eq!(run.fit.and(&run.val).count(), 0, "{at}: fit and holdout");
        for (name, mask) in [
            ("training", &run.train),
            ("fit", &run.fit),
            ("holdout", &run.val),
        ] {
            assert_eq!(
                mask.and(&test).count(),
                0,
                "{at}: the {name} mask holds an outer test sample"
            );
            let resistant = mask.and(&toy.y).count();
            assert!(
                resistant > 0 && resistant < mask.count(),
                "{at}: the {name} mask lacks a class"
            );
        }
        for group in toy.whole(plan.scheme) {
            for (name, mask) in [("outer test fold", &test), ("holdout", &run.val)] {
                assert!(
                    group
                        .iter()
                        .all(|&i| mask.contains(i) == mask.contains(group[0])),
                    "{at}: the {name} cuts {group:?}"
                );
            }
        }
        assert!(
            run.widths.0 > K_SCREEN_CAP && run.widths.1 == K_SCREEN_CAP,
            "{at}: {:?} candidate and screened columns",
            run.widths
        );
        tested = tested.or(&test);
    }
    assert_eq!(tested, Bits::ones(N), "{}: an untested sample", plan.scheme);
}

/// One row per outer test sample in K's row shape, under the verified plan.
/// A called row scores 0.5 or more; a no-call row comes from a fold whose
/// fit kept no tree.
fn check_rows(toy: &Toy, arm: &Arm) {
    let plan = &arm.plan;
    plan.verify().unwrap();
    assert!(
        plan.sha256.len() == 64
            && plan
                .sha256
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    );
    assert_eq!(arm.rows.len(), N);
    for (i, row) in arm.rows.iter().enumerate() {
        let mut fields: Vec<&str> = row
            .as_object()
            .unwrap()
            .keys()
            .map(String::as_str)
            .collect();
        fields.sort_unstable();
        assert_eq!(fields, ROW_FIELDS);
        assert_eq!(row["sample_id"], toy.samples[i]);
        assert_eq!(row["fold"], plan.outer[i]);
        assert_eq!(row["fold_plan_sha256"], plan.sha256);
        let score = row["score"].as_f64().unwrap();
        assert!(
            (0.0..=1.0).contains(&score),
            "{}: score {score}",
            toy.samples[i]
        );
        if row["called"].as_bool().unwrap() {
            assert!(score >= 0.5, "{}: called at {score}", toy.samples[i]);
        }
    }
    // Boosting kept trees somewhere, so the checks compare real fits.
    assert!(
        arm.runs.iter().any(|r| r.scores.is_some()),
        "{}: every fit kept no tree",
        plan.scheme
    );
}

/// No outer test sample reaches a fit: with every outer test sample of fold
/// f relabelled and its pattern bits redrawn, fold f fits exactly what it
/// did. The plan stays as drawn: its outer folds are stratified on every
/// label by design, and G's tests show its inner folds and holdouts read
/// training labels only. Then §16: relabelling boosting's holdout leaves the
/// candidates and the screened block as they were.
fn check_no_test_sample_reaches_a_fit(toy: &Toy, arm: &Arm) {
    let mut rng = ChaCha8Rng::seed_from_u64(SEED + 1);
    for (f, run) in arm.runs.iter().enumerate() {
        let at = format!("{} fold {f}", arm.plan.scheme);
        let test = run.train.not();
        let mut changed = Toy {
            y: flip(&toy.y, &test),
            ..toy.clone()
        };
        changed.redraw(&test, &mut rng);
        assert!(changed.rows != toy.rows, "{at}: nothing redrawn");
        let again = run_fold(&changed.block(), &toy.samples, &changed.y, &arm.plan, f);
        for (stage, (a, b)) in STAGES.iter().zip(again.fitted.iter().zip(&run.fitted)) {
            assert!(a == b, "{at}: an outer test sample changed the {stage}");
        }
        let relabelled = flip(&toy.y, &run.val);
        let again = run_fold(&arm.block, &toy.samples, &relabelled, &arm.plan, f);
        let label_free = STAGES.iter().zip(again.fitted.iter().zip(&run.fitted));
        for (stage, (a, b)) in label_free.skip(1).take(2) {
            assert!(a == b, "{at}: the holdout's labels changed the {stage}");
        }
    }
}

/// A software composition test of G, H, I and one learner on a toy K-only
/// arm. It is not the arm runner, the input optimizer or an experiment; see
/// `docs/ablation-runner.md` for what an arm runner owes. Prints the rows as
/// JSON lines, for K's `validate_predictions`.
#[test]
fn a_toy_k_arm_composes_the_plan_blocks_encoder_and_boosting() {
    let toy = Toy::new();
    let in_pool = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| [Scheme::Random, Scheme::Lineage].map(|scheme| arm(&toy, scheme)))
    };
    let arms = in_pool(1);
    for arm in &arms {
        check_masks(&toy, arm);
        check_rows(&toy, arm);
        check_no_test_sample_reaches_a_fit(&toy, arm);
    }
    for threads in [3, 16] {
        for (a, b) in arms.iter().zip(in_pool(threads)) {
            assert!(
                a.rows == b.rows,
                "{}: rows differ at {threads} threads",
                a.plan.scheme
            );
            for (f, (x, y)) in a.runs.iter().zip(&b.runs).enumerate() {
                assert!(
                    x.fitted == y.fitted,
                    "{} fold {f}: fits differ at {threads} threads",
                    a.plan.scheme
                );
            }
        }
    }
    for row in arms.iter().flat_map(|a| &a.rows) {
        println!("{row}");
    }
}
