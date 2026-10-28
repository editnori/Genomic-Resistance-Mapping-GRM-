//! Synthetic assurance follow-up, 2026-10-27. Extends the fold-mask
//! conformance evidence past n = 0..=5 and k in {2, 3}: larger exhaustive
//! domains, seeded random vectors up to n = 4097 and k = 300, out-of-contract
//! k in {0, 1}, and the caller preconditions that training_folds builds and
//! validate_folds checks. Synthetic integers only. It calls the actual public
//! grm functions and writes evidence for an independent oracle and for Bend.
//!
//! `extended_foldmask_smoke` runs every tier and check at small counts for
//! routine `cargo test`. `extended_foldmask_vectors` is the full recorded run
//! (2,669,501 case IDs, about 0.9 GB of output) and is ignored by default; see
//! docs/synthetic-assurance.md for the release command and expected hashes.
use grm::bits::{Bits, and_count, words_for};
use grm::learn::{fold_masks, training_folds, validate_folds};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use serde_json::{Value, json};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

const SEED: u64 = 0x2026_1027_5eed_0001;
const RANDOM_CASES: usize = 100_000;
const BEND_RANDOM_CASES: usize = 1_500;
const CONTRACT_CASES: usize = 2_000;
const CONTRACT_BEND_CASES: usize = 40;
const VALIDATE_CASES: usize = 50_000;
const TRAINING_CASES: usize = 20_000;
const WORDS_INLINE_MAX_N: usize = 256;
const BOUNDARY_N: [usize; 27] = [
    63, 64, 65, 127, 128, 129, 191, 192, 193, 255, 256, 257, 319, 320, 321, 511, 512, 513, 1023,
    1024, 1025, 2047, 2048, 2049, 4095, 4096, 4097,
];
const BEND_BOUNDARY_N: [usize; 11] = [63, 64, 65, 127, 128, 129, 255, 256, 257, 258, 300];
/// Case IDs in the full run, as recorded by the 2026-10-27 receipts.
const FULL_CASE_IDS: usize = 2_669_501;

/// Case counts for one run. Only the counts differ between runs; every tier,
/// generator and assertion is shared.
struct Scale {
    /// (k, lowest n, highest n) domains of the exhaustive tier.
    exhaustive: &'static [(usize, u32, u32)],
    /// The ext45 Bend book covers every input with n up to this, k in {4, 5}.
    ext45_max_n: u32,
    random: usize,
    bend_random: usize,
    contract: usize,
    contract_bend: usize,
    validate: usize,
    training: usize,
}

/// The recorded evidence. Its output must stay byte-identical to the receipts.
const FULL: Scale = Scale {
    exhaustive: &[
        (2, 6, 6),
        (3, 6, 6),
        (4, 0, 5),
        (5, 0, 5),
        (6, 0, 4),
        (7, 0, 4),
        (8, 0, 4),
    ],
    ext45_max_n: 4,
    random: RANDOM_CASES,
    bend_random: BEND_RANDOM_CASES,
    contract: CONTRACT_CASES,
    contract_bend: CONTRACT_BEND_CASES,
    validate: VALIDATE_CASES,
    training: TRAINING_CASES,
};

/// Seconds in a debug build, for routine `cargo test`.
const SMOKE: Scale = Scale {
    exhaustive: &[
        (2, 0, 4),
        (3, 0, 3),
        (4, 0, 2),
        (5, 0, 2),
        (6, 0, 2),
        (7, 0, 2),
        (8, 0, 2),
    ],
    ext45_max_n: 2,
    random: 1_000,
    bend_random: 150,
    contract: 200,
    contract_bend: CONTRACT_BEND_CASES,
    validate: 2_000,
    training: 1_000,
};

impl Scale {
    /// Case IDs the run must assign, derived from the counts alone.
    fn case_ids(&self) -> usize {
        // alphabet^n assignments times 2^n eligibility masks for each n.
        let domains = |alphabet: usize, lo: u32, hi: u32| -> usize {
            (lo..=hi).map(|n| (alphabet * 2).pow(n)).sum()
        };
        let exhaustive: usize = self
            .exhaustive
            .iter()
            .map(|&(k, lo, hi)| domains(k + 2, lo, hi))
            .sum();
        let ext45: usize = [4, 5]
            .iter()
            .map(|&k| domains(k + 1, 0, self.ext45_max_n))
            .sum();
        exhaustive
            + ext45
            + self.random
            + self.bend_random
            + self.contract
            + self.validate
            + self.training
    }
}

struct FoldOut {
    raw_fit: Bits,
    raw_test: Bits,
    fit: Bits,
    val: Bits,
}

impl FoldOut {
    fn counts(&self) -> [u32; 4] {
        [
            self.raw_fit.count(),
            self.raw_test.count(),
            self.fit.count(),
            self.val.count(),
        ]
    }
}

fn check_bits(label: &str, what: &str, b: &Bits, expected: &[bool]) {
    assert_eq!(b.len(), expected.len(), "{label}: {what}: length");
    assert_eq!(
        b.words().len(),
        words_for(expected.len()),
        "{label}: {what}: word count"
    );
    let mut want = Vec::new();
    for (i, &on) in expected.iter().enumerate() {
        assert_eq!(b.contains(i), on, "{label}: {what}: membership at {i}");
        if on {
            want.push(i);
        }
    }
    assert_eq!(
        b.iter_ones().collect::<Vec<_>>(),
        want,
        "{label}: {what}: iter_ones including tail"
    );
    assert_eq!(b.count() as usize, want.len(), "{label}: {what}: count");
    let rem = b.len() % 64;
    if rem != 0 {
        assert_eq!(
            b.words().last().unwrap() >> rem,
            0,
            "{label}: {what}: cleared tail bits"
        );
    }
}

/// Runs the actual fold_masks and the callers' `.and(eligible)` restriction,
/// and checks every output against the specification computed here.
fn evaluate(label: &str, a: &[usize], eligible: &Bits, k: usize) -> (Vec<FoldOut>, bool) {
    assert_eq!(a.len(), eligible.len());
    let n = a.len();
    let membership: Vec<bool> = (0..n).map(|i| eligible.contains(i)).collect();
    check_bits(label, "eligible input", eligible, &membership);
    let actual = fold_masks(a, k);
    assert_eq!(actual.len(), k, "{label}: fold count");
    let admissible = (0..n).all(|i| !membership[i] || a[i] < k);
    let mut visits = vec![0usize; n];
    let mut outs = Vec::with_capacity(k);
    for (f, (raw_fit, raw_test)) in actual.into_iter().enumerate() {
        let want_test: Vec<bool> = a.iter().map(|&id| id == f).collect();
        let want_raw_fit: Vec<bool> = want_test.iter().map(|&b| !b).collect();
        let want_fit: Vec<bool> = (0..n).map(|i| membership[i] && want_raw_fit[i]).collect();
        let want_val: Vec<bool> = (0..n).map(|i| membership[i] && want_test[i]).collect();
        let fit = raw_fit.and(eligible);
        let val = raw_test.and(eligible);
        check_bits(label, "raw fit", &raw_fit, &want_raw_fit);
        check_bits(label, "raw test", &raw_test, &want_test);
        check_bits(label, "restricted fit", &fit, &want_fit);
        check_bits(label, "restricted validation", &val, &want_val);
        assert_eq!(fit.and(&val).count(), 0, "{label}: disjoint at fold {f}");
        assert_eq!(fit.or(&val), *eligible, "{label}: coverage at fold {f}");
        assert_eq!(
            eligible.and_not(&raw_test),
            fit,
            "{label}: and_not at fold {f}"
        );
        assert_eq!(and_count(raw_fit.words(), eligible.words()), fit.count());
        assert_eq!(and_count(raw_test.words(), eligible.words()), val.count());
        for i in val.iter_ones() {
            visits[i] += 1;
        }
        outs.push(FoldOut {
            raw_fit,
            raw_test,
            fit,
            val,
        });
    }
    for i in 0..n {
        assert_eq!(
            visits[i],
            usize::from(membership[i] && a[i] < k),
            "{label}: validation visits at {i}"
        );
        if k >= 2 && admissible && membership[i] {
            assert_eq!(visits[i], 1, "{label}: unique validation fold at {i}");
        }
    }
    (outs, admissible)
}

fn hex_words(b: &Bits) -> Vec<String> {
    b.words().iter().map(|w| format!("{w:016x}")).collect()
}

/// SHA-256 over the length (u64 LE) and the words (u64 LE) of all four masks,
/// raw fit, raw test, fit, validation. Tail bits are included, so a dirty tail
/// changes the digest.
fn digest_fold(o: &FoldOut) -> String {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&(o.raw_fit.len() as u64).to_le_bytes());
    for b in [&o.raw_fit, &o.raw_test, &o.fit, &o.val] {
        for w in b.words() {
            bytes.extend_from_slice(&w.to_le_bytes());
        }
    }
    ring::digest::digest(&ring::digest::SHA256, &bytes)
        .as_ref()
        .iter()
        .map(|x| format!("{x:02x}"))
        .collect()
}

fn folds_json(outs: &[FoldOut], inline: bool) -> Value {
    Value::Array(
        outs.iter()
            .map(|o| {
                if inline {
                    json!({"c": o.counts(), "w": [hex_words(&o.raw_fit), hex_words(&o.raw_test),
                        hex_words(&o.fit), hex_words(&o.val)]})
                } else {
                    json!({"c": o.counts(), "d": digest_fold(o)})
                }
            })
            .collect(),
    )
}

fn bool_list(b: &Bits) -> String {
    format!(
        "[{}]",
        (0..b.len())
            .map(|i| if b.contains(i) { "True{}" } else { "False{}" })
            .collect::<Vec<_>>()
            .join(",")
    )
}

/// Rust IDs below k become Fold(id); every ID at or above k is Outside.
fn assignment_list(a: &[usize], k: usize) -> String {
    format!(
        "[{}]",
        a.iter()
            .map(|&x| {
                if x < k {
                    format!("M.Fold{{{x}n}}")
                } else {
                    "M.Outside{}".to_owned()
                }
            })
            .collect::<Vec<_>>()
            .join(",")
    )
}

struct BendBook {
    dir: PathBuf,
    manifest: BufWriter<File>,
    shard: Option<BufWriter<File>>,
    per_shard: usize,
    in_shard: usize,
    vectors: usize,
    shards: usize,
}

impl BendBook {
    fn new(dir: PathBuf, per_shard: usize) -> Self {
        fs::create_dir_all(&dir).unwrap();
        let mut manifest = BufWriter::new(File::create(dir.join("VECTORS.bend")).unwrap());
        writeln!(
            manifest,
            "# Rust-produced concrete checks; separate from the universal proof."
        )
        .unwrap();
        Self {
            dir,
            manifest,
            shard: None,
            per_shard,
            in_shard: 0,
            vectors: 0,
            shards: 0,
        }
    }

    fn push(
        &mut self,
        id: usize,
        label: &str,
        a: &[usize],
        eligible: &Bits,
        k: usize,
        outs: &[FoldOut],
    ) {
        if self.shard.is_none() || self.in_shard == self.per_shard {
            if let Some(mut prior) = self.shard.take() {
                prior.flush().unwrap();
            }
            let name = format!("cases_{:04}.bend", self.shards);
            let mut shard = BufWriter::new(File::create(self.dir.join(&name)).unwrap());
            writeln!(
                shard,
                "import Base\nimport ../source/MODEL.bend as M\nimport ./CHECK.bend as C\n"
            )
            .unwrap();
            writeln!(self.manifest, "import ./{name} as V{}", self.shards).unwrap();
            self.shard = Some(shard);
            self.shards += 1;
            self.in_shard = 0;
        }
        // CHECK.all lists folds from k - 1 down to zero.
        let reports: Vec<String> = outs
            .iter()
            .rev()
            .map(|o| {
                let [a1, a2, a3, a4] = o.counts();
                format!(
                    "C.Masks{{{},{},{},{},{a1}n,{a2}n,{a3}n,{a4}n}}",
                    bool_list(&o.raw_fit),
                    bool_list(&o.raw_test),
                    bool_list(&o.fit),
                    bool_list(&o.val)
                )
            })
            .collect();
        writeln!(
            self.shard.as_mut().unwrap(),
            "# Rust input case {id}: {label} n={} k={k}\ndef case_{id}() -> {{C.all({k}n, {}, {}) == [{}] : List<&2, C.Report>}}:\n  {{==}}\n",
            a.len(),
            assignment_list(a, k),
            bool_list(eligible),
            reports.join(",")
        )
        .unwrap();
        self.vectors += 1;
        self.in_shard += 1;
    }

    fn finish(mut self) -> Value {
        if let Some(mut shard) = self.shard.take() {
            shard.flush().unwrap();
        }
        self.manifest.flush().unwrap();
        json!({"vectors": self.vectors, "shards": self.shards})
    }
}

#[derive(Default)]
struct IdStats {
    entries: u64,
    in_range: u64,
    outside: u64,
    equal_k: u64,
    low32_in_range_high: u64,
    at_least_2_63: u64,
    usize_max: u64,
}

impl IdStats {
    fn add(&mut self, a: &[usize], k: usize) {
        for &x in a {
            self.entries += 1;
            if x < k {
                self.in_range += 1;
                continue;
            }
            self.outside += 1;
            self.equal_k += u64::from(x == k);
            self.low32_in_range_high += u64::from(x >= 1usize << 32 && (x as u32 as usize) < k);
            self.at_least_2_63 += u64::from(x >= 1usize << 63);
            self.usize_max += u64::from(x == usize::MAX);
        }
    }

    fn json(&self) -> Value {
        json!({"entries": self.entries, "in_range": self.in_range, "outside": self.outside,
            "equal_k": self.equal_k, "low32_in_range_high": self.low32_in_range_high,
            "at_least_2_63": self.at_least_2_63, "usize_max": self.usize_max})
    }
}

#[derive(Default)]
struct Tier {
    cases: u64,
    fold_rows: u64,
    admissible: u64,
    min_n: Option<usize>,
    max_n: usize,
    min_k: Option<usize>,
    max_k: usize,
    ids: IdStats,
}

impl Tier {
    fn add(&mut self, a: &[usize], k: usize, admissible: bool) {
        self.cases += 1;
        self.fold_rows += k as u64;
        self.admissible += u64::from(admissible);
        let n = a.len();
        self.min_n = Some(self.min_n.map_or(n, |m| m.min(n)));
        self.max_n = self.max_n.max(n);
        self.min_k = Some(self.min_k.map_or(k, |m| m.min(k)));
        self.max_k = self.max_k.max(k);
        self.ids.add(a, k);
    }

    fn json(&self) -> Value {
        json!({"cases": self.cases, "fold_rows": self.fold_rows, "admissible": self.admissible,
            "min_n": self.min_n, "max_n": self.max_n, "min_k": self.min_k, "max_k": self.max_k,
            "ids": self.ids.json()})
    }
}

fn sentinel(rng: &mut ChaCha8Rng, k: usize) -> usize {
    let low = rng.random_range(0..k.max(1));
    match rng.random_range(0..10) {
        0 => k,
        1 => k + 1,
        2 => k + rng.random_range(0..1000),
        3 => rng.random_range(k.max(1)..=u32::MAX as usize),
        // Low 32 bits look like an in-range fold; catches truncating comparisons.
        4 => (rng.random_range(1..=0xffff_ffffusize) << 32) | low,
        5 => 1usize << 63,
        6 => (1usize << 63) | low,
        7 => usize::MAX - rng.random_range(0..1000),
        8 => usize::MAX,
        _ => rng.random_range(k..=usize::MAX),
    }
}

fn assignments(rng: &mut ChaCha8Rng, n: usize, k: usize) -> Vec<usize> {
    let pattern = rng.random_range(0..8);
    let p_out = [0.0, 0.0, 0.01, 0.1, 0.3, 1.0][rng.random_range(0..6)];
    (0..n)
        .map(|i| {
            if k == 0 || rng.random_bool(p_out) {
                return sentinel(rng, k);
            }
            match pattern {
                1 => i % k,
                2 => 0,
                3 => k - 1,
                4 => i * k / n,
                5 => {
                    if rng.random_bool(0.8) {
                        0
                    } else {
                        rng.random_range(0..k)
                    }
                }
                6 => rng.random_range(0..k.min(2)),
                _ => rng.random_range(0..k),
            }
        })
        .collect()
}

fn random_mask(rng: &mut ChaCha8Rng, n: usize) -> Bits {
    match rng.random_range(0..7) {
        0 => Bits::zeros(n),
        1 => Bits::ones(n),
        // Dirty tail words must be cleared by from_words.
        2 => Bits::from_words(n, (0..words_for(n)).map(|_| rng.random::<u64>()).collect()),
        _ => {
            let p = [0.05, 0.5, 0.95, 0.5][rng.random_range(0..4)];
            let idx: Vec<usize> = (0..n).filter(|_| rng.random_bool(p)).collect();
            Bits::from_indices(n, idx)
        }
    }
}

fn eligibility(rng: &mut ChaCha8Rng, a: &[usize], k: usize) -> Bits {
    let n = a.len();
    match rng.random_range(0..9) {
        0 => Bits::zeros(n),
        1 => Bits::ones(n),
        2 => Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0)),
        3 => Bits::from_indices(n, (0..n).filter(|&i| a[i] < k)),
        4 => {
            let idx: Vec<usize> = (0..n)
                .filter(|&i| a[i] < k || rng.random_bool(0.02))
                .collect();
            Bits::from_indices(n, idx)
        }
        5 => Bits::from_indices(n, n.checked_sub(1)),
        6 => Bits::from_indices(n, 0..n.min(1)),
        _ => random_mask(rng, n),
    }
}

fn draw_n(rng: &mut ChaCha8Rng) -> usize {
    match rng.random_range(0..20) {
        0..=6 => rng.random_range(0..=70),
        7..=11 => BOUNDARY_N[rng.random_range(0..BOUNDARY_N.len())],
        12..=16 => rng.random_range(71..=600),
        _ => rng.random_range(601..=4100),
    }
}

fn draw_k(rng: &mut ChaCha8Rng, n: usize) -> usize {
    let k = match rng.random_range(0..40) {
        0..=9 => 5,
        10..=19 => 3,
        20..=31 => rng.random_range(2..=16),
        32..=38 => rng.random_range(17..=64),
        _ => rng.random_range(65..=300),
    };
    if n > 1100 { k.min(16) } else { k }
}

fn draw_bend_n(rng: &mut ChaCha8Rng) -> usize {
    match rng.random_range(0..10) {
        0..=5 => rng.random_range(0..=40),
        6..=7 => BEND_BOUNDARY_N[rng.random_range(0..6)],
        8 => rng.random_range(41..=200),
        _ => BEND_BOUNDARY_N[rng.random_range(6..BEND_BOUNDARY_N.len())],
    }
}

fn draw_bend_k(rng: &mut ChaCha8Rng, n: usize) -> usize {
    match rng.random_range(0..10) {
        0..=2 => 5,
        3..=4 => 3,
        // Fold IDs and k above 256 exercise Bend's U32-backed Nat literals.
        9 if n <= 12 => rng.random_range(250..=262),
        _ => rng.random_range(2..=12),
    }
}

fn mask_json(b: &Bits) -> Value {
    json!({"n": b.len(), "w": hex_words(b)})
}

struct Sink {
    dir: PathBuf,
    next_case: usize,
}

impl Sink {
    fn file(&self, name: &str) -> BufWriter<File> {
        BufWriter::new(File::create(self.dir.join(name)).unwrap())
    }

    fn id(&mut self) -> usize {
        self.next_case += 1;
        self.next_case - 1
    }
}

fn line(out: &mut BufWriter<File>, value: &Value) {
    serde_json::to_writer(&mut *out, value).unwrap();
    writeln!(out).unwrap();
}

fn first_word(b: &Bits) -> String {
    assert!(b.words().len() <= 1);
    b.words()
        .first()
        .map_or_else(|| "-".to_owned(), |w| format!("{w:x}"))
}

/// Exhaustive domains past the prior n = 0..=5, k in {2, 3}. Alphabet:
/// 0..k, k itself, and usize::MAX; every eligibility mask.
fn exhaustive_tier(sink: &mut Sink, scale: &Scale) -> Value {
    let mut out = sink.file("exhaustive.txt");
    let mut tier = Tier::default();
    let mut domains = Vec::new();
    for &(k, lo, hi) in scale.exhaustive {
        for n in lo..=hi {
            let base = k + 2;
            let before = tier.cases;
            for code in 0..base.pow(n) {
                let mut c = code;
                let a: Vec<usize> = (0..n)
                    .map(|_| {
                        let digit = c % base;
                        c /= base;
                        if digit == k + 1 { usize::MAX } else { digit }
                    })
                    .collect();
                let a_csv = if a.is_empty() {
                    "-".to_owned()
                } else {
                    a.iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                };
                for mask in 0..(1usize << n) {
                    let m = Bits::from_indices(
                        n as usize,
                        (0..n as usize).filter(|&i| mask & (1 << i) != 0),
                    );
                    let (outs, admissible) = evaluate("exhaustive", &a, &m, k);
                    tier.add(&a, k, admissible);
                    let id = sink.id();
                    write!(
                        out,
                        "E {id} {k} {n} {code} {mask} {} {a_csv} {}",
                        u8::from(admissible),
                        first_word(&m)
                    )
                    .unwrap();
                    for o in &outs {
                        let [c1, c2, c3, c4] = o.counts();
                        write!(
                            out,
                            " {}:{}:{}:{}:{c1}:{c2}:{c3}:{c4}",
                            first_word(&o.raw_fit),
                            first_word(&o.raw_test),
                            first_word(&o.fit),
                            first_word(&o.val)
                        )
                        .unwrap();
                    }
                    writeln!(out).unwrap();
                }
            }
            domains.push(json!({"k": k, "n": n, "cases": tier.cases - before}));
        }
    }
    out.flush().unwrap();
    json!({"summary": tier.json(), "domains": domains})
}

fn random_tier(sink: &mut Sink, rng: &mut ChaCha8Rng, scale: &Scale) -> Value {
    let mut out = sink.file("random.jsonl");
    let mut tier = Tier::default();
    let mut digested = 0u64;
    for _ in 0..scale.random {
        let n = draw_n(rng);
        let k = draw_k(rng, n);
        let a = assignments(rng, n, k);
        let m = eligibility(rng, &a, k);
        let (outs, admissible) = evaluate("random", &a, &m, k);
        tier.add(&a, k, admissible);
        let inline = n <= WORDS_INLINE_MAX_N && k <= 16;
        digested += u64::from(!inline);
        let id = sink.id();
        line(
            &mut out,
            &json!({"case": id, "tier": "random", "k": k, "a": a, "elig": mask_json(&m),
                "adm": admissible, "folds": folds_json(&outs, inline)}),
        );
    }
    out.flush().unwrap();
    json!({"summary": tier.json(), "digest_only_cases": digested})
}

fn bend_tiers(sink: &mut Sink, rng: &mut ChaCha8Rng, bend_dir: &Path, scale: &Scale) -> Value {
    let mut out = sink.file("bend-cases.jsonl");
    let mut contract_out = sink.file("contract.jsonl");
    let mut ext45 = BendBook::new(bend_dir.join("ext45"), 1024);
    let mut rand_book = BendBook::new(bend_dir.join("rand"), 128);
    let mut ext_tier = Tier::default();
    let mut rand_tier = Tier::default();
    let mut contract_tier = Tier::default();
    let mut contract_emitted = 0u64;
    // Every abstract input for k in {4, 5}, n <= 4 (full run): in-range IDs plus
    // one Outside representative (usize::MAX), every eligibility mask.
    for k in [4usize, 5] {
        for n in 0..=scale.ext45_max_n {
            let base = k + 1;
            for code in 0..base.pow(n) {
                let mut c = code;
                let a: Vec<usize> = (0..n)
                    .map(|_| {
                        let digit = c % base;
                        c /= base;
                        if digit == k { usize::MAX } else { digit }
                    })
                    .collect();
                for mask in 0..(1usize << n) {
                    let m = Bits::from_indices(
                        n as usize,
                        (0..n as usize).filter(|&i| mask & (1 << i) != 0),
                    );
                    let (outs, admissible) = evaluate("ext45", &a, &m, k);
                    ext_tier.add(&a, k, admissible);
                    let id = sink.id();
                    ext45.push(id, "ext45", &a, &m, k, &outs);
                    line(
                        &mut out,
                        &json!({"case": id, "tier": "ext45", "book": "ext45", "k": k, "a": a,
                            "elig": mask_json(&m), "adm": admissible, "folds": folds_json(&outs, true)}),
                    );
                }
            }
        }
    }
    for _ in 0..scale.bend_random {
        let n = draw_bend_n(rng);
        let k = draw_bend_k(rng, n);
        let a = assignments(rng, n, k);
        let m = eligibility(rng, &a, k);
        let (outs, admissible) = evaluate("bend-random", &a, &m, k);
        rand_tier.add(&a, k, admissible);
        let id = sink.id();
        rand_book.push(id, "bend-random", &a, &m, k, &outs);
        line(
            &mut out,
            &json!({"case": id, "tier": "bend-random", "book": "rand", "k": k, "a": a,
                "elig": mask_json(&m), "adm": admissible, "folds": folds_json(&outs, true)}),
        );
    }
    // Out of contract: validate_folds rejects k < 2, so these record raw behavior only.
    for j in 0..scale.contract {
        let k = j % 2;
        let emit = j < scale.contract_bend;
        let n = if emit {
            rng.random_range(0..=20)
        } else {
            draw_bend_n(rng)
        };
        let a = assignments(rng, n, k);
        let m = eligibility(rng, &a, k);
        let (outs, admissible) = evaluate("contract", &a, &m, k);
        contract_tier.add(&a, k, admissible);
        let id = sink.id();
        if emit {
            rand_book.push(id, "out-of-contract-k", &a, &m, k, &outs);
            contract_emitted += 1;
        }
        line(
            &mut contract_out,
            &json!({"case": id, "tier": "contract", "book": if emit { "rand" } else { "" }, "k": k,
                "a": a, "elig": mask_json(&m), "adm": admissible,
                "folds": folds_json(&outs, n <= WORDS_INLINE_MAX_N)}),
        );
    }
    out.flush().unwrap();
    contract_out.flush().unwrap();
    json!({"ext45": {"summary": ext_tier.json(), "book": ext45.finish()},
        "bend_random": {"summary": rand_tier.json()},
        "contract": {"summary": contract_tier.json(), "emitted_to_rand_book": contract_emitted},
        "rand_book": rand_book.finish()})
}

fn validate_tier(sink: &mut Sink, rng: &mut ChaCha8Rng, scale: &Scale) -> Value {
    let mut out = sink.file("validate.jsonl");
    let mut ok = 0u64;
    let mut kinds = std::collections::BTreeMap::<String, u64>::new();
    for _ in 0..scale.validate {
        let structured = rng.random_bool(0.3);
        let k = if structured {
            rng.random_range(2..=6)
        } else {
            rng.random_range(0..=8)
        };
        let n = if structured {
            rng.random_range(4 * k..=4 * k + 120)
        } else {
            match rng.random_range(0..10) {
                0..=6 => rng.random_range(0..=40),
                7..=8 => rng.random_range(41..=200),
                _ => [63, 64, 65, 127, 128, 129, 255, 256, 257][rng.random_range(0..9)],
            }
        };
        let (y, folds, train) = if structured {
            // Class-wise round robin: each class has at least 2k members, so
            // every fold's training and validation sets start with both classes.
            let y = Bits::from_indices(n, (0..n).filter(|i| i % 2 == 0));
            let mut seen = [0usize; 2];
            let mut folds: Vec<usize> = (0..n)
                .map(|i| {
                    let c = usize::from(y.contains(i));
                    seen[c] += 1;
                    (seen[c] - 1) % k
                })
                .collect();
            let mut train = Bits::ones(n);
            if rng.random_bool(0.3) {
                let i = rng.random_range(0..n);
                folds[i] = sentinel(rng, k);
                if rng.random_bool(0.5) {
                    train.remove(i);
                }
            }
            if rng.random_bool(0.3) {
                let f = rng.random_range(0..k);
                let class = rng.random_bool(0.5);
                for (i, &fold) in folds.iter().enumerate() {
                    if fold == f && y.contains(i) == class && rng.random_bool(0.9) {
                        train.remove(i);
                    }
                }
            }
            (y, folds, train)
        } else {
            let nf = match rng.random_range(0..20) {
                0 => n + 1,
                1 => n.saturating_sub(1),
                _ => n,
            };
            let nt = match rng.random_range(0..20) {
                0 => n + 1,
                1 => n.saturating_sub(1),
                _ => n,
            };
            let y = random_mask(rng, n);
            let folds = assignments(rng, nf, k);
            let train = if nt == nf && rng.random_bool(0.4) {
                Bits::from_indices(nt, (0..nt).filter(|&i| folds[i] < k))
            } else {
                random_mask(rng, nt)
            };
            (y, folds, train)
        };
        let result = match validate_folds(&y, &train, &folds, k) {
            Ok(()) => {
                ok += 1;
                "ok".to_owned()
            }
            Err(e) => format!("error: {e:#}"),
        };
        let kind = result
            .split(' ')
            .take(if result.starts_with("error: fold") {
                4
            } else {
                3
            })
            .collect::<Vec<_>>()
            .join(" ");
        *kinds.entry(kind).or_default() += 1;
        let id = sink.id();
        line(
            &mut out,
            &json!({"case": id, "tier": "validate", "k": k, "y": mask_json(&y), "train": mask_json(&train),
                "folds": folds, "result": result}),
        );
    }
    out.flush().unwrap();
    json!({"cases": scale.validate, "ok": ok, "result_prefixes": kinds})
}

fn training_tier(sink: &mut Sink, rng: &mut ChaCha8Rng, scale: &Scale) -> Value {
    let mut out = sink.file("training.jsonl");
    let mut tier = Tier::default();
    let mut grouped = 0u64;
    let mut valid = 0u64;
    for _ in 0..scale.training {
        let n = draw_bend_n(rng);
        let k = rng.random_range(2..=8);
        let seed: u64 = rng.random();
        let y = random_mask(rng, n);
        let train = random_mask(rng, n);
        let groups: Option<Vec<u32>> = if rng.random_bool(0.5) {
            None
        } else {
            let g = [1, 2, 3, (n / 4).max(1), n.max(1)][rng.random_range(0..5)] as u32;
            Some((0..n).map(|_| rng.random_range(0..g)).collect())
        };
        grouped += u64::from(groups.is_some());
        let folds = training_folds(&y, &train, groups.as_deref(), k, seed);
        assert_eq!(folds.len(), n, "training: length");
        for (i, &fold) in folds.iter().enumerate() {
            if train.contains(i) {
                assert!(fold < k, "training: member {i} has no fold");
            } else {
                assert_eq!(fold, usize::MAX, "training: outside sentinel at {i}");
            }
        }
        if let Some(g) = &groups {
            let mut fold_of = std::collections::BTreeMap::new();
            for i in train.iter_ones() {
                let first = *fold_of.entry(g[i]).or_insert(folds[i]);
                assert_eq!(folds[i], first, "training: group {} split", g[i]);
            }
        }
        // Labels and groups outside `train` must not change the inner folds.
        let mut y2 = y.clone();
        let mut g2 = groups.clone();
        for i in 0..n {
            if !train.contains(i) && rng.random_bool(0.5) {
                if y2.contains(i) {
                    y2.remove(i);
                } else {
                    y2.insert(i);
                }
                if let Some(g) = g2.as_mut() {
                    g[i] = rng.random_range(0..1000);
                }
            }
        }
        assert_eq!(
            training_folds(&y2, &train, g2.as_deref(), k, seed),
            folds,
            "training: outside samples changed the inner folds"
        );
        let (_, admissible) = evaluate("training", &folds, &train, k);
        assert!(
            admissible,
            "training: built folds violate the fold-mask precondition"
        );
        tier.add(&folds, k, admissible);
        let result = match validate_folds(&y, &train, &folds, k) {
            Ok(()) => {
                valid += 1;
                "ok".to_owned()
            }
            Err(e) => format!("error: {e:#}"),
        };
        let id = sink.id();
        line(
            &mut out,
            &json!({"case": id, "tier": "training", "k": k, "seed": seed, "y": mask_json(&y),
                "train": mask_json(&train), "groups": groups, "folds": folds, "validate": result}),
        );
    }
    out.flush().unwrap();
    json!({"summary": tier.json(), "grouped": grouped, "validate_ok": valid})
}

/// The directory named by `var`, which must be empty so nothing is
/// overwritten, or else a temporary directory removed when the guard drops.
fn output_dir(var: &str) -> (PathBuf, Option<tempfile::TempDir>) {
    let Some(dir) = std::env::var_os(var) else {
        let tmp = tempfile::tempdir().unwrap();
        return (tmp.path().to_owned(), Some(tmp));
    };
    let dir = PathBuf::from(dir);
    fs::create_dir_all(&dir).unwrap();
    assert!(
        fs::read_dir(&dir).unwrap().next().is_none(),
        "{var} must be empty"
    );
    (dir, None)
}

/// Runs every tier at `scale`, writes the evidence into `dir`, and checks the
/// number of case IDs against the count derived from `scale`.
fn run(dir: &Path, scale: &Scale) {
    let bend_dir = dir.join("bend");
    let mut sink = Sink {
        dir: dir.to_owned(),
        next_case: 0,
    };
    let mut rng = ChaCha8Rng::seed_from_u64(SEED);
    let exhaustive = exhaustive_tier(&mut sink, scale);
    let random = random_tier(&mut sink, &mut rng, scale);
    let bend = bend_tiers(&mut sink, &mut rng, &bend_dir, scale);
    let validate = validate_tier(&mut sink, &mut rng, scale);
    let training = training_tier(&mut sink, &mut rng, scale);
    assert_eq!(sink.next_case, scale.case_ids(), "case ID count");
    let summary = json!({
        "seed": format!("{SEED:#x}"), "usize_bits": usize::BITS, "total_case_ids": sink.next_case,
        "debug_assertions": cfg!(debug_assertions),
        "exhaustive": exhaustive, "random": random, "bend": bend,
        "validate": validate, "training": training
    });
    fs::write(
        dir.join("summary.json"),
        serde_json::to_vec_pretty(&summary).unwrap(),
    )
    .unwrap();
    println!("{summary}");
}

/// Every tier and check at small counts, written to a temporary directory.
#[test]
fn extended_foldmask_smoke() {
    let tmp = tempfile::tempdir().unwrap();
    run(tmp.path(), &SMOKE);
}

/// The full recorded run: about 36 s in release and 317 s in debug, about
/// 0.9 GB written. GRM_SYNTH_OUT, if set, names an empty directory that keeps
/// the output for hashing; otherwise it goes to a temporary directory.
#[test]
#[ignore = "full run, ~0.9 GB; use --release, see docs/synthetic-assurance.md"]
fn extended_foldmask_vectors() {
    assert_eq!(FULL.case_ids(), FULL_CASE_IDS);
    let (dir, _tmp) = output_dir("GRM_SYNTH_OUT");
    run(&dir, &FULL);
}
