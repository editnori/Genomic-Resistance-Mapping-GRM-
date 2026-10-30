//! Synthetic assurance follow-up, 2026-10-27: multi-k differential vectors for
//! grm::dna::KmerSpace at every supported k (5..=32). Random synthetic byte
//! strings and integers only; an independent oracle re-derives every value.
//!
//! `multi_k_space_vectors` runs at full size in routine `cargo test` (under
//! 20 s in debug, about 30 MB of output). The k = 13, 14 bijection check is ignored by
//! default; see docs/synthetic-assurance.md.
use grm::dna::{KmerSpace, MAX_K, MIN_K};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use serde_json::{Value, json};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::PathBuf;

const SEED: u64 = 0x2026_1027_5eed_0002;
const SEQS_PER_K: usize = 250;
const VALUES_PER_K: usize = 4_000;
const TEXTS_PER_K: usize = 1_000;
/// Every value of the 2k-bit space is mixed for k <= 12 (2^24 values).
const EXHAUSTIVE_MAX_K: usize = 12;
/// Full mix tables are digested for k <= 10 so the oracle can recompute them.
const TABLE_DIGEST_MAX_K: usize = 10;

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn sha256_hex(bytes: &[u8]) -> String {
    hex(ring::digest::digest(&ring::digest::SHA256, bytes).as_ref())
}

fn random_seq(rng: &mut ChaCha8Rng, k: usize) -> Vec<u8> {
    let len = match rng.random_range(0..10) {
        0..=5 => rng.random_range(0..=80),
        6..=8 => rng.random_range(81..=400),
        _ => rng.random_range(401..=2000),
    };
    match rng.random_range(0..4) {
        0 => (0..len).map(|_| b"ACGT"[rng.random_range(0..4)]).collect(),
        1 | 2 => (0..len)
            .map(|_| match rng.random_range(0..100) {
                0..=69 => b"ACGT"[rng.random_range(0..4)],
                70..=84 => b"acgt"[rng.random_range(0..4)],
                85..=89 => b'N',
                90..=91 => b'n',
                92..=95 => b"RYKMSWBDHVU-*."[rng.random_range(0..14)],
                _ => rng.random::<u8>(),
            })
            .collect(),
        _ => {
            // Clean runs around k long separated by 1..=5 window-breaking bytes.
            let mut s = Vec::with_capacity(len + 8);
            while s.len() < len {
                for _ in 0..rng.random_range(0..=2 * k + 2) {
                    s.push(b"ACGTacgt"[rng.random_range(0..8)]);
                }
                for _ in 0..rng.random_range(1..=5) {
                    s.push(b"NnRY-.\n\r\0\xff"[rng.random_range(0..10)]);
                }
            }
            s.truncate(len);
            s
        }
    }
}

fn exhaustive_bijection(space: &KmerSpace, mask: u64, digest: bool) -> Value {
    let size = mask as usize + 1;
    let mut seen = vec![0u64; size.div_ceil(64)];
    let mut table = Vec::with_capacity(if digest { size * 8 } else { 0 });
    for x in 0..=mask {
        let y = space.mix(x);
        assert!(y <= mask, "k={}: mix escaped the space", space.k);
        let (w, b) = ((y >> 6) as usize, y & 63);
        assert_eq!((seen[w] >> b) & 1, 0, "k={}: mix collision at {y}", space.k);
        seen[w] |= 1u64 << b;
        assert_eq!(space.unmix(y), x, "k={}: unmix", space.k);
        if digest {
            table.extend_from_slice(&y.to_le_bytes());
        }
    }
    let covered: u64 = seen.iter().map(|w| u64::from(w.count_ones())).sum();
    assert_eq!(covered, mask + 1, "k={}: mix is not onto", space.k);
    json!({"k": space.k, "values": mask + 1, "permutation": true,
        "table_sha256": digest.then(|| sha256_hex(&table))})
}

fn space_mask(space: &KmerSpace) -> u64 {
    if space.bits() == 64 {
        u64::MAX
    } else {
        (1u64 << space.bits()) - 1
    }
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

#[test]
fn multi_k_space_vectors() {
    let (dir, _tmp) = output_dir("GRM_KMER_OUT");
    let mut rng = ChaCha8Rng::seed_from_u64(SEED);
    let mut seqs = BufWriter::new(File::create(dir.join("sequences.jsonl")).unwrap());
    let mut values = BufWriter::new(File::create(dir.join("values.jsonl")).unwrap());
    let mut texts = BufWriter::new(File::create(dir.join("texts.jsonl")).unwrap());
    let mut tables = Vec::new();
    let (mut n_seq, mut n_kmers, mut n_values, mut n_texts) = (0u64, 0u64, 0u64, 0u64);
    for k in MIN_K..=MAX_K {
        let space = KmerSpace::new(k);
        assert_eq!(space.k, k);
        assert_eq!(space.bits() as usize, 2 * k);
        let mask = space_mask(&space);
        for _ in 0..SEQS_PER_K {
            let seq = random_seq(&mut rng, k);
            let mut canon = Vec::new();
            space.for_each_canonical(&seq, |x| canon.push(x));
            let mixed: Vec<u64> = canon.iter().map(|&x| space.mix(x)).collect();
            for (&x, &y) in canon.iter().zip(&mixed) {
                assert!(x <= mask && y <= mask, "k={k}: value escaped the space");
                assert_eq!(space.unmix(y), x, "k={k}: unmix(mix(x))");
                assert_eq!(space.reverse_complement(space.reverse_complement(x)), x);
                assert!(
                    x <= space.reverse_complement(x),
                    "k={k}: canonical is the smaller strand"
                );
            }
            n_seq += 1;
            n_kmers += canon.len() as u64;
            serde_json::to_writer(
                &mut seqs,
                &json!({"k": k, "seq": hex(&seq), "canon": canon, "mixed": mixed}),
            )
            .unwrap();
            writeln!(seqs).unwrap();
        }
        let mut xs: Vec<u64> = vec![
            0,
            1,
            2,
            3,
            mask,
            mask - 1,
            mask >> 1,
            (mask >> 1) + 1,
            0x5555_5555_5555_5555 & mask,
            0xaaaa_aaaa_aaaa_aaaa & mask,
            0x0f0f_0f0f_0f0f_0f0f & mask,
            0x1b1b_1b1b_1b1b_1b1b & mask,
        ];
        xs.extend((0..VALUES_PER_K).map(|_| rng.random::<u64>() & mask));
        for x in xs {
            assert_eq!(space.mix(space.unmix(x)), x, "k={k}: mix(unmix(x))");
            serde_json::to_writer(
                &mut values,
                &json!({"k": k, "x": x, "mix": space.mix(x), "unmix": space.unmix(x),
                    "rc": space.reverse_complement(x), "decode": space.decode(x)}),
            )
            .unwrap();
            writeln!(values).unwrap();
            n_values += 1;
        }
        for _ in 0..TEXTS_PER_K {
            let len = match rng.random_range(0..10) {
                0 => k - 1,
                1 => k + 1,
                2 => 0,
                _ => k,
            };
            let text: String = (0..len)
                .map(|_| {
                    if rng.random_range(0..100) < 85 {
                        b"ACGTacgt"[rng.random_range(0..8)] as char
                    } else {
                        b"NnRYX-"[rng.random_range(0..6)] as char
                    }
                })
                .collect();
            serde_json::to_writer(
                &mut texts,
                &json!({"k": k, "text": text, "canonical": space.canonical(&text)}),
            )
            .unwrap();
            writeln!(texts).unwrap();
            n_texts += 1;
        }
        if k <= EXHAUSTIVE_MAX_K {
            tables.push(exhaustive_bijection(&space, mask, k <= TABLE_DIGEST_MAX_K));
        }
    }
    for w in [&mut seqs, &mut values, &mut texts] {
        w.flush().unwrap();
    }
    let summary = json!({"seed": format!("{SEED:#x}"), "k_range": [MIN_K, MAX_K],
        "sequences": n_seq, "canonical_kmers": n_kmers, "values": n_values, "texts": n_texts,
        "exhaustive_max_k": EXHAUSTIVE_MAX_K, "table_digest_max_k": TABLE_DIGEST_MAX_K,
        "tables": tables});
    fs::write(
        dir.join("summary.json"),
        serde_json::to_vec_pretty(&summary).unwrap(),
    )
    .unwrap();
    println!("{summary}");
}

/// Release-profile extension of the exhaustive permutation check to k = 13, 14
/// (2^26 and 2^28 values, a 32 MiB bitmap). It runs only when requested with
/// --ignored. GRM_KMER_LARGE_OUT, if set, names a new file for the result.
#[test]
#[ignore = "k = 13, 14 bijection; use --release, see docs/synthetic-assurance.md"]
fn multi_k_mixer_bijection_large() {
    let out = std::env::var_os("GRM_KMER_LARGE_OUT").map(|out| File::create_new(out).unwrap());
    let tables: Vec<Value> = (13..=14)
        .map(|k| {
            let space = KmerSpace::new(k);
            exhaustive_bijection(&space, space_mask(&space), false)
        })
        .collect();
    let result = json!({"tables": tables});
    if let Some(mut out) = out {
        out.write_all(&serde_json::to_vec_pretty(&result).unwrap())
            .unwrap();
    }
    println!("{result}");
}
