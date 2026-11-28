//! Format integrity, damaged-file resilience, count-encoding boundary,
//! fixed-size regression, modality header, and end-to-end CLI tests for
//! `grm::kcount::format` and `grm kcount convert-v1` / `inspect`.

use grm::dna::KmerSpace;
use grm::kcount::format::{
    AnyTableReader, DEFAULT_SAMPLE_NAME, GrKsetReader, InputModality, KtabReader, KtabWriteMeta,
    SequenceProfile, TableFormat, check_output_extension, convert_to_gr, convert_to_ktab,
    convert_to_v1, decode_bucket_counts, decode_bucket_keys, decode_profile_runs,
    encode_bucket_counts, encode_bucket_keys, read_profile_header, read_profiles,
    read_profiles_verified, write_gr_kset, write_ktab, write_ktab_from_buckets, write_profiles,
};
use grm::kcount::key::{KmerKey, KmerSpace128};
use grm::kset::{self, BUCKET_BITS, BUCKETS, KsetReader};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use ring::digest::{Context as Sha256Context, SHA256};
use std::fs;
use std::path::Path;
use std::process::Command;

const GRK_DIGEST_START: usize = 40;
const GRK_DIGEST_END: usize = 72;
const KTAB_DIGEST_START: usize = 48;
const KTAB_DIGEST_END: usize = 80;
const PROF_DIGEST_START: usize = 40;
const PROF_DIGEST_END: usize = 72;

fn rehash_slot(bytes: &mut [u8], slot_start: usize, slot_end: usize) {
    let mut sha = Sha256Context::new(&SHA256);
    sha.update(&bytes[..slot_start]);
    sha.update(&bytes[slot_end..]);
    let digest = sha.finish();
    bytes[slot_start..slot_end].copy_from_slice(digest.as_ref());
}

fn sorted_u64_keys(seed: u64, k: usize, n: usize) -> Vec<u64> {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let mut vals: Vec<u64> = (0..n).map(|_| rng.random::<u64>() & mask).collect();
    vals.sort_unstable();
    vals.dedup();
    vals
}

fn sorted_u128_keys(seed: u64, k: usize, n: usize) -> Vec<u128> {
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let mask = if k == 64 {
        u128::MAX
    } else {
        (1u128 << (2 * k)) - 1
    };
    let mut vals: Vec<u128> = (0..n)
        .map(|_| (((rng.random::<u64>() as u128) << 64) | (rng.random::<u64>() as u128)) & mask)
        .collect();
    vals.sort_unstable();
    vals.dedup();
    vals
}

fn grm_bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_grm"))
}

fn write_and_verify_grk_fails(dir: &Path, tag: &str, bytes: &[u8]) {
    let p = dir.join(format!("{tag}.grk"));
    fs::write(&p, bytes).unwrap();
    let res = GrKsetReader::open(&p).and_then(|r| r.read_all::<u64>());
    assert!(res.is_err(), "expected damaged grk {tag} to fail");
}

fn write_and_verify_ktab_fails(dir: &Path, tag: &str, bytes: &[u8]) {
    let p = dir.join(format!("{tag}.ktab"));
    fs::write(&p, bytes).unwrap();
    let res = KtabReader::open(&p).and_then(|r| r.read_all::<u64>());
    assert!(res.is_err(), "expected damaged ktab {tag} to fail");
}

fn write_and_verify_prof_fails(dir: &Path, tag: &str, bytes: &[u8]) {
    let p = dir.join(format!("{tag}.prof.bin"));
    fs::write(&p, bytes).unwrap();
    assert!(
        read_profiles(&p).is_err(),
        "expected damaged prof {tag} to fail"
    );
}

// ---------------------------------------------------------------------------
// 1. Every single-byte flip anywhere in .grk, .ktab, or .prof.bin is caught
// ---------------------------------------------------------------------------

#[test]
fn grk_ktab_and_prof_reject_every_single_byte_flip() {
    let dir = tempfile::tempdir().unwrap();

    // Build compact valid .grk, .ktab (counted), and .prof.bin files.
    let keys31 = sorted_u64_keys(1001, 31, 200);
    let grk_path = dir.path().join("intact.grk");
    write_gr_kset(&grk_path, "sample", 31, InputModality::Assembly, &keys31).unwrap();
    let grk_bytes = fs::read(&grk_path).unwrap();

    let entries31: Vec<(u64, u32)> = keys31
        .iter()
        .enumerate()
        .map(|(i, &k)| {
            let c = match i % 5 {
                0 => 1,
                1 => 255,
                2 => 65_535,
                3 => 65_536,
                _ => u32::MAX,
            };
            (k, c)
        })
        .collect();
    let sat = entries31.iter().filter(|&&(_, c)| c == u32::MAX).count() as u64;
    let ktab_path = dir.path().join("intact.ktab");
    write_ktab(
        &ktab_path,
        "sample",
        31,
        InputModality::Reads,
        1,
        sat,
        &entries31,
    )
    .unwrap();
    let ktab_bytes = fs::read(&ktab_path).unwrap();

    let profs = vec![
        SequenceProfile {
            name: "seq_a".to_string(),
            seq_len: 80,
            counts: vec![1, 1, 2, 2, 2, 0, 0, 300, 65_535, 65_536, u32::MAX, 1],
        },
        SequenceProfile {
            name: "seq_b".to_string(),
            seq_len: 40,
            counts: vec![5, 5, 5, 5],
        },
    ];
    let prof_path = dir.path().join("intact.prof.bin");
    write_profiles(&prof_path, 31, &profs).unwrap();
    let prof_bytes = fs::read(&prof_path).unwrap();

    // Flip EVERY single byte of .prof.bin (header + digest + body).
    let flip_prof = dir.path().join("flip.prof.bin");
    for pos in 0..prof_bytes.len() {
        let mut bad = prof_bytes.clone();
        bad[pos] ^= 0x01;
        fs::write(&flip_prof, &bad).unwrap();
        assert!(
            read_profiles(&flip_prof).is_err(),
            "prof.bin byte flip at offset {pos} of {} went undetected",
            prof_bytes.len()
        );
    }

    // Flip EVERY header/digest/name byte, EVERY payload byte, and sampled
    // bucket-table positions in .grk.
    let grk_data_start = 72 + "sample".len() + (BUCKETS + 1) * 8 + BUCKETS * 4;
    let flip_grk = dir.path().join("flip.grk");
    let mut grk_flips = 0usize;
    for pos in 0..grk_bytes.len() {
        if pos >= 72 + "sample".len() && pos < grk_data_start && pos % 37 != 0 {
            continue;
        }
        let mut bad = grk_bytes.clone();
        bad[pos] ^= 0x01;
        fs::write(&flip_grk, &bad).unwrap();
        let res = GrKsetReader::open(&flip_grk).and_then(|r| r.read_all::<u64>());
        assert!(
            res.is_err(),
            "grk byte flip at offset {pos} of {} went undetected",
            grk_bytes.len()
        );
        grk_flips += 1;
    }
    assert!(
        grk_flips > 1_000,
        "expected >1000 grk flip checks, got {grk_flips}"
    );

    // Flip EVERY header/digest/name byte, EVERY key+count payload byte, and
    // sampled bucket-table positions in .ktab.
    let ktab_data_start = 80 + "sample".len() + 2 * (BUCKETS + 1) * 8 + BUCKETS * 4;
    let flip_ktab = dir.path().join("flip.ktab");
    let mut ktab_flips = 0usize;
    for pos in 0..ktab_bytes.len() {
        if pos >= 80 + "sample".len() && pos < ktab_data_start && pos % 41 != 0 {
            continue;
        }
        let mut bad = ktab_bytes.clone();
        bad[pos] ^= 0x01;
        fs::write(&flip_ktab, &bad).unwrap();
        let res = KtabReader::open(&flip_ktab).and_then(|r| r.read_all::<u64>());
        assert!(
            res.is_err(),
            "ktab byte flip at offset {pos} of {} went undetected",
            ktab_bytes.len()
        );
        ktab_flips += 1;
    }
    assert!(
        ktab_flips > 1_500,
        "expected >1500 ktab flip checks, got {ktab_flips}"
    );
}

// ---------------------------------------------------------------------------
// 2. Damaged-file stress test (900+ damaged files across .grk, .ktab, .prof.bin):
//    truncations, flips, wrong magic, huge lengths (u32::MAX, u64::MAX,
//    1 << 31, 1 << 63) both raw and rehashed -> 0 panics, 100% Err(...)
// ---------------------------------------------------------------------------

#[test]
fn damaged_files_stress_zero_panics_across_grk_ktab_and_prof() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(20261127);

    let keys = sorted_u64_keys(77, 31, 150);
    let grk_path = dir.path().join("base.grk");
    write_gr_kset(&grk_path, "sample", 31, InputModality::Assembly, &keys).unwrap();
    let grk_base = fs::read(&grk_path).unwrap();

    let entries: Vec<(u64, u32)> = keys
        .iter()
        .enumerate()
        .map(|(i, &k)| (k, (i % 500 + 1) as u32))
        .collect();
    let ktab_path = dir.path().join("base.ktab");
    write_ktab(
        &ktab_path,
        "sample",
        31,
        InputModality::Reads,
        1,
        0,
        &entries,
    )
    .unwrap();
    let ktab_base = fs::read(&ktab_path).unwrap();

    let profs = vec![
        SequenceProfile {
            name: "c1".to_string(),
            seq_len: 120,
            counts: vec![1, 1, 2, 2, 10, 10, 0, 0, 500, 70_000],
        },
        SequenceProfile {
            name: "c2".to_string(),
            seq_len: 50,
            counts: vec![3, 3, 3],
        },
    ];
    let prof_path = dir.path().join("base.prof.bin");
    write_profiles(&prof_path, 31, &profs).unwrap();
    let prof_base = fs::read(&prof_path).unwrap();

    let huge_u32s: &[u32] = &[u32::MAX, u32::MAX - 1, 1u32 << 31, (1u32 << 31) + 7];
    let huge_u64s: &[u64] = &[
        u64::MAX,
        u64::MAX - 1,
        1u64 << 63,
        1u64 << 62,
        1u64 << 40,
        1u64 << 31,
        u32::MAX as u64,
    ];

    let mut total_cases = 0usize;

    // --- .grk damaged cases (300+) ---
    // Truncations
    let mut grk_cuts = vec![0usize, 1, 7, 8, 39, 40, 71, 72, 73, 100, grk_base.len() - 1];
    for _ in 0..80 {
        grk_cuts.push(rng.random_range(0..grk_base.len()));
    }
    for (i, cut) in grk_cuts.into_iter().enumerate() {
        write_and_verify_grk_fails(dir.path(), &format!("trunc_{i}"), &grk_base[..cut]);
        total_cases += 1;
    }
    // Wrong magic & random flips & multi-byte corruption
    for i in 0..160 {
        let mut bad = grk_base.clone();
        if i < 40 {
            let at = i % 8;
            bad[at] ^= 0x55;
        } else {
            let n_flips = rng.random_range(1..=4);
            for _ in 0..n_flips {
                let at = rng.random_range(0..bad.len());
                bad[at] ^= rng.random_range(1u8..=255);
            }
        }
        write_and_verify_grk_fails(dir.path(), &format!("mut_{i}"), &bad);
        total_cases += 1;
    }
    // Huge lengths (both raw and rehashed!) on name_len, n_kmers, bucket offsets, bucket counts
    let grk_table_off = 72 + "sample".len();
    let grk_counts_off = grk_table_off + (BUCKETS + 1) * 8;
    for (idx, &v32) in huge_u32s.iter().enumerate() {
        for rehash in [false, true] {
            // Huge name_len
            let mut bad = grk_base.clone();
            bad[24..28].copy_from_slice(&v32.to_le_bytes());
            if rehash {
                rehash_slot(&mut bad, GRK_DIGEST_START, GRK_DIGEST_END);
            }
            write_and_verify_grk_fails(dir.path(), &format!("hname_{idx}_{rehash}"), &bad);
            total_cases += 1;

            // Huge bucket_count[0] + raised n_kmers (exact Finding 1 repro!)
            let mut bad = grk_base.clone();
            let old_c0 =
                u32::from_le_bytes(bad[grk_counts_off..grk_counts_off + 4].try_into().unwrap());
            let old_n = u64::from_le_bytes(bad[16..24].try_into().unwrap());
            let new_n = old_n - (old_c0 as u64) + (v32 as u64);
            bad[grk_counts_off..grk_counts_off + 4].copy_from_slice(&v32.to_le_bytes());
            bad[16..24].copy_from_slice(&new_n.to_le_bytes());
            if rehash {
                rehash_slot(&mut bad, GRK_DIGEST_START, GRK_DIGEST_END);
            }
            write_and_verify_grk_fails(dir.path(), &format!("hcnt_{idx}_{rehash}"), &bad);
            total_cases += 1;
        }
    }
    for (idx, &v64) in huge_u64s.iter().enumerate() {
        for rehash in [false, true] {
            // Huge n_kmers
            let mut bad = grk_base.clone();
            bad[16..24].copy_from_slice(&v64.to_le_bytes());
            if rehash {
                rehash_slot(&mut bad, GRK_DIGEST_START, GRK_DIGEST_END);
            }
            write_and_verify_grk_fails(dir.path(), &format!("hnk_{idx}_{rehash}"), &bad);
            total_cases += 1;

            // Huge bucket offset
            let mut bad = grk_base.clone();
            bad[grk_table_off + 8..grk_table_off + 16].copy_from_slice(&v64.to_le_bytes());
            if rehash {
                rehash_slot(&mut bad, GRK_DIGEST_START, GRK_DIGEST_END);
            }
            write_and_verify_grk_fails(dir.path(), &format!("hoff_{idx}_{rehash}"), &bad);
            total_cases += 1;
        }
    }

    // --- .ktab damaged cases (300+) ---
    let mut ktab_cuts = vec![
        0usize,
        1,
        7,
        8,
        47,
        48,
        79,
        80,
        81,
        120,
        ktab_base.len() - 1,
    ];
    for _ in 0..80 {
        ktab_cuts.push(rng.random_range(0..ktab_base.len()));
    }
    for (i, cut) in ktab_cuts.into_iter().enumerate() {
        write_and_verify_ktab_fails(dir.path(), &format!("trunc_{i}"), &ktab_base[..cut]);
        total_cases += 1;
    }
    for i in 0..160 {
        let mut bad = ktab_base.clone();
        if i < 40 {
            let at = i % 8;
            bad[at] ^= 0xaa;
        } else {
            let n_flips = rng.random_range(1..=4);
            for _ in 0..n_flips {
                let at = rng.random_range(0..bad.len());
                bad[at] ^= rng.random_range(1u8..=255);
            }
        }
        write_and_verify_ktab_fails(dir.path(), &format!("mut_{i}"), &bad);
        total_cases += 1;
    }
    let ktab_table_off = 80 + "sample".len();
    let ktab_counts_off = ktab_table_off + 2 * (BUCKETS + 1) * 8;
    for (idx, &v32) in huge_u32s.iter().enumerate() {
        for rehash in [false, true] {
            let mut bad = ktab_base.clone();
            bad[40..44].copy_from_slice(&v32.to_le_bytes());
            if rehash {
                rehash_slot(&mut bad, KTAB_DIGEST_START, KTAB_DIGEST_END);
            }
            write_and_verify_ktab_fails(dir.path(), &format!("hname_{idx}_{rehash}"), &bad);
            total_cases += 1;

            let mut bad = ktab_base.clone();
            let old_c0 = u32::from_le_bytes(
                bad[ktab_counts_off..ktab_counts_off + 4]
                    .try_into()
                    .unwrap(),
            );
            let old_n = u64::from_le_bytes(bad[24..32].try_into().unwrap());
            let new_n = old_n - (old_c0 as u64) + (v32 as u64);
            bad[ktab_counts_off..ktab_counts_off + 4].copy_from_slice(&v32.to_le_bytes());
            bad[24..32].copy_from_slice(&new_n.to_le_bytes());
            if rehash {
                rehash_slot(&mut bad, KTAB_DIGEST_START, KTAB_DIGEST_END);
            }
            write_and_verify_ktab_fails(dir.path(), &format!("hcnt_{idx}_{rehash}"), &bad);
            total_cases += 1;
        }
    }
    for (idx, &v64) in huge_u64s.iter().enumerate() {
        for rehash in [false, true] {
            let mut bad = ktab_base.clone();
            bad[24..32].copy_from_slice(&v64.to_le_bytes());
            if rehash {
                rehash_slot(&mut bad, KTAB_DIGEST_START, KTAB_DIGEST_END);
            }
            write_and_verify_ktab_fails(dir.path(), &format!("hnk_{idx}_{rehash}"), &bad);
            total_cases += 1;

            let mut bad = ktab_base.clone();
            bad[ktab_table_off + 8..ktab_table_off + 16].copy_from_slice(&v64.to_le_bytes());
            if rehash {
                rehash_slot(&mut bad, KTAB_DIGEST_START, KTAB_DIGEST_END);
            }
            write_and_verify_ktab_fails(dir.path(), &format!("hoff_{idx}_{rehash}"), &bad);
            total_cases += 1;
        }
    }

    // --- .prof.bin damaged cases (300+) ---
    for cut in 0..prof_base.len() {
        write_and_verify_prof_fails(dir.path(), &format!("trunc_{cut}"), &prof_base[..cut]);
        total_cases += 1;
    }
    for i in 0..160 {
        let mut bad = prof_base.clone();
        if i < 40 {
            let at = i % 8;
            bad[at] ^= 0xff;
        } else {
            let n_flips = rng.random_range(1..=4);
            for _ in 0..n_flips {
                let at = rng.random_range(0..bad.len());
                bad[at] ^= rng.random_range(1u8..=255);
            }
        }
        write_and_verify_prof_fails(dir.path(), &format!("mut_{i}"), &bad);
        total_cases += 1;
    }
    // Huge lengths inside .prof.bin header and record 0 (both raw and rehashed!)
    let rec0 = 72usize;
    let rec0_meta = rec0 + 4 + "c1".len();
    for (idx, &v32) in huge_u32s.iter().enumerate() {
        for rehash in [false, true] {
            let mut bad = prof_base.clone();
            bad[rec0..rec0 + 4].copy_from_slice(&v32.to_le_bytes());
            if rehash {
                rehash_slot(&mut bad, PROF_DIGEST_START, PROF_DIGEST_END);
            }
            write_and_verify_prof_fails(dir.path(), &format!("hname_{idx}_{rehash}"), &bad);
            total_cases += 1;
        }
    }
    for (idx, &v64) in huge_u64s.iter().enumerate() {
        for rehash in [false, true] {
            // Huge n_sequences
            let mut bad = prof_base.clone();
            bad[16..24].copy_from_slice(&v64.to_le_bytes());
            if rehash {
                rehash_slot(&mut bad, PROF_DIGEST_START, PROF_DIGEST_END);
            }
            write_and_verify_prof_fails(dir.path(), &format!("hnseq_{idx}_{rehash}"), &bad);
            total_cases += 1;

            // Huge seq_len (would overflow total_bases sum if unchecked!)
            let mut bad = prof_base.clone();
            bad[rec0_meta..rec0_meta + 8].copy_from_slice(&v64.to_le_bytes());
            bad[24..32].copy_from_slice(&v64.to_le_bytes());
            if rehash {
                rehash_slot(&mut bad, PROF_DIGEST_START, PROF_DIGEST_END);
            }
            write_and_verify_prof_fails(dir.path(), &format!("hseqlen_{idx}_{rehash}"), &bad);
            total_cases += 1;

            // Huge n_positions (with huge seq_len so n_pos <= seq_len still checked!)
            let mut bad = prof_base.clone();
            bad[rec0_meta..rec0_meta + 8].copy_from_slice(&v64.to_le_bytes());
            bad[rec0_meta + 8..rec0_meta + 16].copy_from_slice(&v64.to_le_bytes());
            if rehash {
                rehash_slot(&mut bad, PROF_DIGEST_START, PROF_DIGEST_END);
            }
            write_and_verify_prof_fails(dir.path(), &format!("hnpos_{idx}_{rehash}"), &bad);
            total_cases += 1;

            // Huge payload_len
            let mut bad = prof_base.clone();
            bad[rec0_meta + 16..rec0_meta + 24].copy_from_slice(&v64.to_le_bytes());
            if rehash {
                rehash_slot(&mut bad, PROF_DIGEST_START, PROF_DIGEST_END);
            }
            write_and_verify_prof_fails(dir.path(), &format!("hpay_{idx}_{rehash}"), &bad);
            total_cases += 1;
        }
    }

    assert!(
        total_cases >= 900,
        "expected >= 900 damaged-file cases, ran {total_cases}"
    );
}

// ---------------------------------------------------------------------------
// 3. Exact 65,535 / 65,536 / u32::MAX counts across BOTH tiered and u16-escape
//    modes (kills format.rs:383 and :462 planted faults!)
// ---------------------------------------------------------------------------

#[test]
fn ktab_exact_65535_65536_u32_max_and_both_count_modes() {
    let dir = tempfile::tempdir().unwrap();

    // Bucket A: majority in 254..65534 -> forces COUNT_MODE_U16_ESC (mode byte = 2)
    // AND includes exact boundaries 65_534, 65_535 (0xFFFF), 65_536, u32::MAX!
    let mut u16_counts: Vec<u32> = (0..30).map(|i| 300 + i * 100).collect();
    u16_counts.extend_from_slice(&[254, 255, 65_534, 65_535, 65_536, 70_000, u32::MAX]);
    let mut u16_enc = Vec::new();
    encode_bucket_counts(&u16_counts, &mut u16_enc);
    assert_eq!(u16_enc[0], 2, "expected COUNT_MODE_U16_ESC (2)");
    let mut u16_dec = Vec::new();
    decode_bucket_counts(0, u16_counts.len() as u32, false, &u16_enc, &mut u16_dec).unwrap();
    assert_eq!(u16_dec, u16_counts);

    // Bucket B: majority < 254 -> forces COUNT_MODE_TIERED (mode byte = 1)
    // AND includes exact boundaries 253, 254, 255, 65_534, 65_535, 65_536, u32::MAX!
    let mut tiered_counts: Vec<u32> = (1..=30).collect();
    tiered_counts.extend_from_slice(&[253, 254, 255, 65_534, 65_535, 65_536, u32::MAX]);
    let mut tiered_enc = Vec::new();
    encode_bucket_counts(&tiered_counts, &mut tiered_enc);
    assert_eq!(tiered_enc[0], 1, "expected COUNT_MODE_TIERED (1)");
    let mut tiered_dec = Vec::new();
    decode_bucket_counts(
        1,
        tiered_counts.len() as u32,
        false,
        &tiered_enc,
        &mut tiered_dec,
    )
    .unwrap();
    assert_eq!(tiered_dec, tiered_counts);

    // Bucket C: all ones -> forces COUNT_MODE_ALL_ONE (mode byte = 0).
    let ones_counts = vec![1u32; 25];
    let mut ones_enc = Vec::new();
    encode_bucket_counts(&ones_counts, &mut ones_enc);
    assert_eq!(ones_enc, vec![0u8]);
    let mut ones_dec = Vec::new();
    decode_bucket_counts(2, ones_counts.len() as u32, false, &ones_enc, &mut ones_dec).unwrap();
    assert_eq!(ones_dec, ones_counts);

    // Full .ktab file containing all three bucket modes in buckets 0, 1, 2.
    let k = 31usize;
    let mut entries: Vec<(u64, u32)> = Vec::new();
    for (i, &c) in u16_counts.iter().enumerate() {
        entries.push((u64::bucket_floor(0, k, BUCKET_BITS) + (i as u64) * 10, c));
    }
    for (i, &c) in tiered_counts.iter().enumerate() {
        entries.push((u64::bucket_floor(1, k, BUCKET_BITS) + (i as u64) * 10, c));
    }
    for (i, &c) in ones_counts.iter().enumerate() {
        entries.push((u64::bucket_floor(2, k, BUCKET_BITS) + (i as u64) * 10, c));
    }
    let p = dir.path().join("all_modes.ktab");
    write_ktab(&p, "sample", k, InputModality::Reads, 1, 2, &entries).unwrap();
    let r = KtabReader::open_verified(&p).unwrap();
    assert!(!r.presence_only);
    assert_eq!(r.saturated_counts, 2);
    assert_eq!(r.max_count_observed, u32::MAX);
    assert_eq!(r.read_all::<u64>().unwrap(), entries);
}

// ---------------------------------------------------------------------------
// 4. Every reader and decoder validation check tested directly & via rehashed
//    malformed files (kills format.rs:301, :332, :338 and all header checks)
// ---------------------------------------------------------------------------

#[test]
fn rehashed_reader_validation_checks_reject_every_malformed_field() {
    let dir = tempfile::tempdir().unwrap();

    // --- Direct decode_bucket_keys checks ---
    let mut out64 = Vec::new();
    // Empty bucket with trailing data
    assert!(decode_bucket_keys::<u64>(31, 0, 0, &[0], &mut out64).is_err());
    // Non-empty bucket with empty payload
    assert!(decode_bucket_keys::<u64>(31, 0, 1, &[], &mut out64).is_err());
    // r > low_bits (format.rs:301): for k=31, low_bits = 52 -> r = 53 must fail!
    assert!(
        decode_bucket_keys::<u64>(31, 0, 1, &[53, 0, 0, 0, 0, 0, 0, 0, 0], &mut out64).is_err()
    );
    // Decoded key exceeds bucket upper bound max_val (format.rs:332):
    // k = 5 has low_bits = 0, so floor == max_val; gap = 1 (unary 1 = 0b01) exceeds max_val!
    assert!(decode_bucket_keys::<u64>(5, 0, 1, &[0, 0b0000_0001], &mut out64).is_err());
    // Trailing non-zero padding bits in last byte (format.rs:338):
    let mut valid_b0 = Vec::new();
    encode_bucket_keys::<u64>(31, 0, &[0, 5, 19], &mut valid_b0);
    out64.clear();
    decode_bucket_keys::<u64>(31, 0, 3, &valid_b0, &mut out64).unwrap();
    assert_eq!(out64, vec![0, 5, 19]);
    let mut with_extra_byte = valid_b0.clone();
    with_extra_byte.push(0);
    out64.clear();
    assert!(decode_bucket_keys::<u64>(31, 0, 3, &with_extra_byte, &mut out64).is_err());
    // Single key at k=5 uses 1 bit (0) in byte 1; setting bit 1 leaves non-zero trailing bit.
    let mut bad_pad = vec![0u8, 0b0000_0010];
    out64.clear();
    assert!(decode_bucket_keys::<u64>(5, 0, 1, &bad_pad, &mut out64).is_err());
    bad_pad[1] = 0;
    out64.clear();
    decode_bucket_keys::<u64>(5, 0, 1, &bad_pad, &mut out64).unwrap();
    assert_eq!(out64, vec![0]);
    // Truncated unary / truncated remainder
    assert!(decode_bucket_keys::<u64>(31, 0, 1, &[8, 0xff, 0xff], &mut out64).is_err());
    assert!(decode_bucket_keys::<u64>(31, 0, 1, &[16, 0x00], &mut out64).is_err());

    // --- Direct decode_bucket_counts checks ---
    let mut cout = Vec::new();
    assert!(decode_bucket_counts(0, 0, false, &[0], &mut cout).is_err());
    assert!(decode_bucket_counts(0, 1, true, &[0], &mut cout).is_err());
    assert!(decode_bucket_counts(0, 1, false, &[], &mut cout).is_err());
    assert!(decode_bucket_counts(0, 1, false, &[0, 0], &mut cout).is_err()); // all-one trailing
    assert!(decode_bucket_counts(0, 1, false, &[99], &mut cout).is_err()); // unknown mode
    assert!(decode_bucket_counts(0, 1, false, &[1, 254, 1], &mut cout).is_err()); // tiered trunc u16
    assert!(decode_bucket_counts(0, 1, false, &[1, 255, 1, 2], &mut cout).is_err()); // tiered trunc u32
    assert!(decode_bucket_counts(0, 1, false, &[1, 5, 0], &mut cout).is_err()); // tiered trailing
    assert!(decode_bucket_counts(0, 1, false, &[2, 0xff, 0xff, 1, 2], &mut cout).is_err()); // u16 trunc u32
    assert!(decode_bucket_counts(0, 1, false, &[2, 5, 0, 99], &mut cout).is_err()); // u16 trailing

    // --- Direct decode_profile_runs checks ---
    assert!(decode_profile_runs(&[0], 0).is_err());
    assert!(decode_profile_runs(&[], 1).is_err());
    // Negative decoded count (delta = -1 from prev = 0 -> zigzag 1)
    assert!(decode_profile_runs(&[1, 0], 1).is_err());
    // Too many positions (delta 0, run = 5 when n_positions = 3)
    assert!(decode_profile_runs(&[0, 4], 3).is_err());
    // Too few positions (delta 0, run = 2 when n_positions = 3)
    assert!(decode_profile_runs(&[0, 1], 3).is_err());
    // Truncated varint
    assert!(decode_profile_runs(&[0x80], 1).is_err());
    // Varint overflow (10 continuation bytes)
    assert!(decode_profile_runs(&[0x80; 11], 1).is_err());

    // --- Rehashed .grk header/table checks ---
    let keys = sorted_u64_keys(909, 31, 64);
    let grk_path = dir.path().join("ok.grk");
    write_gr_kset(&grk_path, "ok", 31, InputModality::Assembly, &keys).unwrap();
    let grk_ok = fs::read(&grk_path).unwrap();

    type Mutator = (&'static str, fn(&mut Vec<u8>));

    let grk_cases: &[Mutator] = &[
        ("bad_magic", |b| b[0] = b'X'),
        ("bad_k_low", |b| b[8] = 4),
        ("bad_k_high", |b| b[8] = 65),
        ("bad_bucket_bits", |b| b[9] = 9),
        ("bad_key_width", |b| b[10] = 16),
        ("bad_modality", |b| b[11] = 3),
        ("nonzero_reserved0", |b| b[12] = 1),
        ("nonzero_reserved1", |b| b[28] = 1),
        ("bad_utf8_name", |b| b[72] = 0xff),
        ("offsets0_nonzero", |b| b[74] = 1),
        ("trailing_byte", |b| b.push(0)),
    ];
    for &(label, mutate) in grk_cases {
        let mut bad = grk_ok.clone();
        mutate(&mut bad);
        rehash_slot(&mut bad, GRK_DIGEST_START, GRK_DIGEST_END);
        write_and_verify_grk_fails(dir.path(), &format!("rehash_grk_{label}"), &bad);
    }

    // --- Rehashed .ktab header/table checks ---
    let entries: Vec<(u64, u32)> = keys.iter().map(|&k| (k, 5)).collect();
    let ktab_path = dir.path().join("ok.ktab");
    write_ktab(&ktab_path, "ok", 31, InputModality::Reads, 1, 0, &entries).unwrap();
    let ktab_ok = fs::read(&ktab_path).unwrap();

    let ktab_cases: &[Mutator] = &[
        ("bad_magic", |b| b[0] = b'X'),
        ("bad_k", |b| b[8] = 0),
        ("bad_bucket_bits", |b| b[9] = 11),
        ("bad_key_width", |b| b[10] = 16),
        ("bad_flags", |b| b[11] = 0x04),
        ("bad_modality", |b| b[12] = 7),
        ("nonzero_reserved0", |b| b[13] = 1),
        ("nonzero_reserved1", |b| b[44] = 1),
        ("sat_flag_without_count", |b| b[11] |= 0x02),
        ("sat_count_without_flag", |b| {
            b[32..40].copy_from_slice(&1u64.to_le_bytes());
        }),
        ("wrong_max_count", |b| {
            b[20..24].copy_from_slice(&99u32.to_le_bytes());
        }),
        ("bad_utf8_name", |b| b[80] = 0xff),
        ("trailing_byte", |b| b.push(0)),
    ];
    for &(label, mutate) in ktab_cases {
        let mut bad = ktab_ok.clone();
        mutate(&mut bad);
        rehash_slot(&mut bad, KTAB_DIGEST_START, KTAB_DIGEST_END);
        write_and_verify_ktab_fails(dir.path(), &format!("rehash_ktab_{label}"), &bad);
    }

    // --- Rehashed .prof.bin header/body checks ---
    let profs = vec![SequenceProfile {
        name: "s1".to_string(),
        seq_len: 50,
        counts: vec![1, 2, 3, 4],
    }];
    let prof_path = dir.path().join("ok.prof.bin");
    write_profiles(&prof_path, 31, &profs).unwrap();
    let prof_ok = fs::read(&prof_path).unwrap();

    let prof_cases: &[Mutator] = &[
        ("bad_magic", |b| b[0] = b'Z'),
        ("bad_k", |b| b[8] = 65),
        ("nonzero_reserved0", |b| b[9] = 1),
        ("nonzero_reserved1", |b| b[32] = 1),
        ("wrong_max_count", |b| {
            b[12..16].copy_from_slice(&99u32.to_le_bytes());
        }),
        ("wrong_total_bases", |b| {
            b[24..32].copy_from_slice(&51u64.to_le_bytes());
        }),
        ("bad_utf8_seq_name", |b| b[76] = 0xff),
        ("trailing_byte", |b| b.push(0)),
    ];
    for &(label, mutate) in prof_cases {
        let mut bad = prof_ok.clone();
        mutate(&mut bad);
        rehash_slot(&mut bad, PROF_DIGEST_START, PROF_DIGEST_END);
        write_and_verify_prof_fails(dir.path(), &format!("rehash_prof_{label}"), &bad);
    }
}

// ---------------------------------------------------------------------------
// 5. Fixed-size regression test on deterministic synthetic k=31 and k=61 sets
//    (kills format.rs:237 Golomb-Rice parameter size fault!)
// ---------------------------------------------------------------------------

#[test]
fn fixed_size_regression_on_synthetic_k31_and_k61_sets() {
    let dir = tempfile::tempdir().unwrap();

    // Deterministic synthetic k=31 set (20,000 mixed keys).
    let keys31 = sorted_u64_keys(314_159, 31, 20_000);
    assert_eq!(keys31.len(), 20_000);
    let space31 = KmerSpace::new(31);

    let v1_path = dir.path().join("fixed31.kset");
    kset::write(&v1_path, DEFAULT_SAMPLE_NAME, &space31, &keys31).unwrap();
    let v1_bytes = fs::metadata(&v1_path).unwrap().len();

    let grk31_path = dir.path().join("fixed31.grk");
    write_gr_kset(
        &grk31_path,
        DEFAULT_SAMPLE_NAME,
        31,
        InputModality::Assembly,
        &keys31,
    )
    .unwrap();
    let grk31_bytes = fs::metadata(&grk31_path).unwrap().len();

    let ones31: Vec<(u64, u32)> = keys31.iter().map(|&k| (k, 1)).collect();
    let ktab31_ones = dir.path().join("fixed31_ones.ktab");
    write_ktab(
        &ktab31_ones,
        DEFAULT_SAMPLE_NAME,
        31,
        InputModality::Assembly,
        1,
        0,
        &ones31,
    )
    .unwrap();
    let ktab31_ones_bytes = fs::metadata(&ktab31_ones).unwrap().len();

    let counted31: Vec<(u64, u32)> = keys31
        .iter()
        .enumerate()
        .map(|(i, &k)| (k, ((i % 20) + 1) as u32))
        .collect();
    let ktab31_cnt = dir.path().join("fixed31_cnt.ktab");
    write_ktab(
        &ktab31_cnt,
        DEFAULT_SAMPLE_NAME,
        31,
        InputModality::Reads,
        1,
        0,
        &counted31,
    )
    .unwrap();
    let ktab31_cnt_bytes = fs::metadata(&ktab31_cnt).unwrap().len();

    // Deterministic synthetic k=61 set (10,000 mixed keys).
    let keys61 = sorted_u128_keys(616_161, 61, 10_000);
    assert_eq!(keys61.len(), 10_000);
    let grk61_path = dir.path().join("fixed61.grk");
    write_gr_kset(
        &grk61_path,
        DEFAULT_SAMPLE_NAME,
        61,
        InputModality::Reads,
        &keys61,
    )
    .unwrap();
    let grk61_bytes = fs::metadata(&grk61_path).unwrap().len();

    // Exact pinned byte sizes: any Golomb-Rice parameter regression (+2.1% etc.)
    // immediately fails these exact equality assertions.
    assert_eq!(ktab31_ones_bytes, grk31_bytes + 8_208);
    assert!(
        grk31_bytes < v1_bytes,
        "grk ({grk31_bytes}) must be smaller than kset v1 ({v1_bytes})"
    );
    assert_eq!(grk31_bytes, 136_483, "fixed k=31 .grk byte size changed");
    assert_eq!(
        ktab31_ones_bytes, 144_691,
        "fixed k=31 presence-only .ktab byte size changed"
    );
    assert_eq!(
        ktab31_cnt_bytes, 165_715,
        "fixed k=31 counted .ktab byte size changed"
    );
    assert_eq!(grk61_bytes, 151_338, "fixed k=61 .grk byte size changed");
}

// ---------------------------------------------------------------------------
// 6. Modality header round-trip & .kset v1 interop across k = 5..=64
// ---------------------------------------------------------------------------

#[test]
fn modality_header_round_trip_and_v1_interop() {
    let dir = tempfile::tempdir().unwrap();

    for modality in [
        InputModality::Unknown,
        InputModality::Assembly,
        InputModality::Reads,
    ] {
        for k in [5usize, 21, 31, 32] {
            let space = KmerSpace::new(k);
            let keys = sorted_u64_keys(100 + k as u64, k, 500);
            let v1 = dir.path().join(format!("m{}_{k}.kset", modality.as_byte()));
            kset::write(&v1, "sample", &space, &keys).unwrap();

            let gr = dir.path().join(format!("m{}_{k}.grk", modality.as_byte()));
            let out_gr = convert_to_gr(&v1, &gr, "sample", Some(modality)).unwrap();
            assert_eq!(out_gr.format, TableFormat::GrKset);
            assert_eq!(out_gr.modality, Some(modality));
            assert_eq!(out_gr.n_kmers as usize, keys.len());
            assert!(!out_gr.counts_dropped);

            let tab = dir.path().join(format!("m{}_{k}.ktab", modality.as_byte()));
            let out_tab = convert_to_ktab(&gr, &tab, "sample", None).unwrap();
            assert_eq!(out_tab.modality, Some(modality));

            // Round-trip back to .kset v1 must be byte-for-byte identical.
            let v1_back = dir
                .path()
                .join(format!("m{}_{k}_back.kset", modality.as_byte()));
            convert_to_v1(&tab, &v1_back, "sample").unwrap();
            assert_eq!(fs::read(&v1).unwrap(), fs::read(&v1_back).unwrap());

            // Unified AnyTableReader agrees across all three formats.
            let rv1 = AnyTableReader::open(&v1).unwrap();
            let rgr = AnyTableReader::open(&gr).unwrap();
            let rtab = AnyTableReader::open(&tab).unwrap();
            assert_eq!(rv1.modality(), None);
            assert_eq!(rgr.modality(), Some(modality));
            assert_eq!(rtab.modality(), Some(modality));
            assert!(!rv1.verify_integrity().unwrap());
            assert!(rgr.verify_integrity().unwrap());
            assert!(rtab.verify_integrity().unwrap());
            assert_eq!(
                rv1.sketch::<u64>(20).unwrap(),
                rgr.sketch::<u64>(20).unwrap()
            );
            assert_eq!(
                rv1.sketch::<u64>(20).unwrap(),
                rtab.sketch::<u64>(20).unwrap()
            );
        }
    }

    // Wide u128 tables (k = 33, 61, 64) preserve modality and convert .grk <-> .ktab.
    for k in [33usize, 61, 64] {
        let space = KmerSpace128::new(k);
        let keys = sorted_u128_keys(200 + k as u64, k, 300);
        assert!(keys.iter().all(|&x| x <= space.mask()));
        let gr = dir.path().join(format!("wide_{k}.grk"));
        write_gr_kset(&gr, "sample", k, InputModality::Reads, &keys).unwrap();
        let tab = dir.path().join(format!("wide_{k}.ktab"));
        let out = convert_to_ktab(&gr, &tab, "sample", None).unwrap();
        assert_eq!(out.modality, Some(InputModality::Reads));
        let r = KtabReader::open_verified(&tab).unwrap();
        assert_eq!(r.modality, InputModality::Reads);
        assert_eq!(
            r.read_all::<u128>()
                .unwrap()
                .into_iter()
                .map(|(key, _)| key)
                .collect::<Vec<_>>(),
            keys
        );
        // Converting k > 32 to .kset v1 must fail cleanly.
        assert!(convert_to_v1(&gr, &dir.path().join("bad.kset"), "sample").is_err());
    }

    // Extension guard rejects mismatches.
    assert!(check_output_extension(Path::new("a.ktab"), TableFormat::GrKset).is_err());
    assert!(check_output_extension(Path::new("a.grk"), TableFormat::KsetV1).is_err());
    assert!(check_output_extension(Path::new("a.kset"), TableFormat::Ktab).is_err());
    assert!(check_output_extension(Path::new("no_ext"), TableFormat::GrKset).is_err());

    // write_ktab_from_buckets rejects inconsistent metadata cleanly.
    let empty_counts = [0u32; BUCKETS];
    let empty_payloads = vec![Vec::new(); BUCKETS];
    let bad_meta = KtabWriteMeta {
        name: "sample",
        k: 31,
        modality: InputModality::Unknown,
        min_count: 1,
        max_count_observed: 5,
        saturated_counts: 0,
        presence_only: false,
    };
    assert!(
        write_ktab_from_buckets(
            &dir.path().join("bad_empty.ktab"),
            &bad_meta,
            &empty_counts,
            &empty_payloads,
            &empty_payloads,
        )
        .is_err()
    );
}

// ---------------------------------------------------------------------------
// 7. End-to-end CLI tests for `grm kcount convert-v1` and `grm kcount inspect`
// ---------------------------------------------------------------------------

#[test]
fn cli_convert_v1_and_inspect_end_to_end() {
    let dir = tempfile::tempdir().unwrap();
    let k = 31usize;
    let space = KmerSpace::new(k);
    let keys = sorted_u64_keys(4242, k, 1_200);
    let want_fnv = kset::digest(&keys);

    // Intentionally distinct input filename stem to verify Finding 25:
    // omitted `--name` defaults to `"sample"` and NEVER leaks the input filename.
    let secret_stem = "arbitrary_input_stem_98765";
    let v1_in = dir.path().join(format!("{secret_stem}.kset"));
    kset::write(&v1_in, secret_stem, &space, &keys).unwrap();
    let v1_orig_reader = KsetReader::open(&v1_in).unwrap();
    assert_eq!(v1_orig_reader.name, secret_stem);

    // 1) convert-v1 .kset -> .grk WITHOUT --name (must default to "sample"!).
    let gr_out = dir.path().join("converted.grk");
    let out = grm_bin()
        .args([
            "kcount",
            "convert-v1",
            "--input",
            v1_in.to_str().unwrap(),
            "--out",
            gr_out.to_str().unwrap(),
            "--modality",
            "assembly",
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "convert-v1 failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    assert!(
        !stdout.contains(secret_stem),
        "convert-v1 stdout leaked input filename stem"
    );
    let conv_json: serde_json::Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(conv_json["format"], "GRMKGR1");
    assert_eq!(conv_json["k"], k);
    assert_eq!(conv_json["sample"], DEFAULT_SAMPLE_NAME);
    assert_eq!(conv_json["modality"], "assembly");
    assert_eq!(conv_json["n_kmers"], keys.len() as u64);
    assert_eq!(conv_json["counts_dropped"], false);
    assert_eq!(GrKsetReader::open_verified(&gr_out).unwrap().name, "sample");

    // 2) convert-v1 .kset -> .ktab WITH explicit --name and --modality.
    let tab_out = dir.path().join("named.ktab");
    let out = grm_bin()
        .args([
            "kcount",
            "convert-v1",
            "--input",
            v1_in.to_str().unwrap(),
            "--out",
            tab_out.to_str().unwrap(),
            "--name",
            "custom_tag",
            "--modality",
            "reads",
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    let tab_reader = KtabReader::open_verified(&tab_out).unwrap();
    assert_eq!(tab_reader.name, "custom_tag");
    assert_eq!(tab_reader.modality, InputModality::Reads);
    assert!(tab_reader.presence_only);

    // 3) Extension mismatch refusals (Finding 24).
    for bad_args in [
        vec![
            "kcount",
            "convert-v1",
            "--input",
            v1_in.to_str().unwrap(),
            "--out",
            dir.path().join("wrong.ktab").to_str().unwrap(),
            "--direction",
            "to-gr",
        ],
        vec![
            "kcount",
            "convert-v1",
            "--input",
            gr_out.to_str().unwrap(),
            "--out",
            dir.path().join("wrong.grk").to_str().unwrap(),
            "--direction",
            "to-v1",
        ],
        vec![
            "kcount",
            "convert-v1",
            "--input",
            v1_in.to_str().unwrap(),
            "--out",
            dir.path().join("wrong.kset").to_str().unwrap(),
            "--direction",
            "to-ktab",
        ],
        vec![
            "kcount",
            "convert-v1",
            "--input",
            v1_in.to_str().unwrap(),
            "--out",
            dir.path().join("wrong.bin").to_str().unwrap(),
        ],
    ] {
        let res = grm_bin().args(&bad_args).output().unwrap();
        assert!(
            !res.status.success(),
            "expected extension mismatch to fail for {bad_args:?}"
        );
        let stderr = String::from_utf8_lossy(&res.stderr);
        assert!(
            stderr.contains("disagrees"),
            "expected clear extension mismatch error, got: {stderr}"
        );
    }

    // 4) Counted .ktab -> .kset v1 reports counts_dropped: true.
    let counted_entries: Vec<(u64, u32)> = keys
        .iter()
        .enumerate()
        .map(|(i, &key)| (key, if i == 0 { 65_535 } else { 2 }))
        .collect();
    let counted_tab = dir.path().join("counted.ktab");
    write_ktab(
        &counted_tab,
        "sample",
        k,
        InputModality::Reads,
        1,
        0,
        &counted_entries,
    )
    .unwrap();
    let v1_out = dir.path().join("from_counted.kset");
    let out = grm_bin()
        .args([
            "kcount",
            "convert-v1",
            "--input",
            counted_tab.to_str().unwrap(),
            "--out",
            v1_out.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    let conv_cnt: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(conv_cnt["format"], "GRMKSET1");
    assert_eq!(conv_cnt["sample"], "sample");
    assert_eq!(conv_cnt["counts_dropped"], true);

    // 5) Build a .prof.bin for inspect tests.
    let prof_path = dir.path().join("sample.prof.bin");
    let profs = vec![SequenceProfile {
        name: "contig1".to_string(),
        seq_len: 100,
        counts: vec![1, 1, 5, 65_535, 65_536],
    }];
    write_profiles(&prof_path, k, &profs).unwrap();
    let (prof_hdr, _) = read_profiles_verified(&prof_path).unwrap();
    assert_eq!(read_profile_header(&prof_path).unwrap(), prof_hdr);

    // 6) inspect WITHOUT --verify MUST report verified: false on ALL FOUR formats!
    for (path, want_fmt) in [
        (&v1_out, "kset-v1 (GRMKSET1)"),
        (&gr_out, "gr-kset (GRMKGR1)"),
        (&counted_tab, "ktab (GRMKTB1)"),
        (&prof_path, "profile (GRMPRF1)"),
    ] {
        let out = grm_bin()
            .args(["kcount", "inspect", "--input", path.to_str().unwrap()])
            .output()
            .unwrap();
        assert!(out.status.success());
        let j: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
        assert_eq!(j["format"], want_fmt);
        assert_eq!(j["k"], k);
        assert_eq!(
            j["verified"], false,
            "inspect without --verify must set verified=false for {want_fmt}"
        );
        assert!(j["fnv64"].is_null());
        assert!(j["keys_sha256"].is_null());
        assert!(j.get("path").is_none(), "inspect must not echo input path");
    }

    // 7) inspect --verify on intact .grk, .ktab, .prof.bin -> verified: true;
    //    on .kset v1 -> verified: false (no stored digest) with exact fnv64.
    let insp_gr_file = dir.path().join("insp_gr.json");
    let out = grm_bin()
        .args([
            "kcount",
            "inspect",
            "--input",
            gr_out.to_str().unwrap(),
            "--verify",
            "--sketch",
            "16",
            "--out",
            insp_gr_file.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    let j_gr: serde_json::Value =
        serde_json::from_slice(&fs::read(&insp_gr_file).unwrap()).unwrap();
    assert_eq!(j_gr["format"], "gr-kset (GRMKGR1)");
    assert_eq!(j_gr["verified"], true);
    assert_eq!(j_gr["sample"], "sample");
    assert_eq!(j_gr["modality"], "assembly");
    assert_eq!(j_gr["n_kmers"], keys.len() as u64);
    assert_eq!(j_gr["fnv64"], want_fnv);
    assert_eq!(j_gr["header_sha256"].as_str().unwrap().len(), 64);
    assert_eq!(j_gr["keys_sha256"].as_str().unwrap().len(), 64);
    assert_eq!(j_gr["sketch_hex"].as_array().unwrap().len(), 16);

    let out = grm_bin()
        .args([
            "kcount",
            "inspect",
            "--input",
            counted_tab.to_str().unwrap(),
            "--verify",
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    let j_tab: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(j_tab["format"], "ktab (GRMKTB1)");
    assert_eq!(j_tab["verified"], true);
    assert_eq!(j_tab["modality"], "reads");
    assert_eq!(j_tab["presence_only"], false);
    assert_eq!(j_tab["max_count_observed"], 65_535);
    assert_eq!(j_tab["saturated_counts"], 0);
    assert_eq!(j_tab["fnv64"], want_fnv);

    let out = grm_bin()
        .args([
            "kcount",
            "inspect",
            "--input",
            prof_path.to_str().unwrap(),
            "--verify",
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    let j_prof: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(j_prof["format"], "profile (GRMPRF1)");
    assert_eq!(j_prof["verified"], true);
    assert_eq!(j_prof["n_sequences"], 1);
    assert_eq!(j_prof["total_bases"], 100);
    assert_eq!(j_prof["max_count_observed"], 65_536);
    assert_eq!(j_prof["header_sha256"].as_str().unwrap().len(), 64);

    let out = grm_bin()
        .args([
            "kcount",
            "inspect",
            "--input",
            v1_out.to_str().unwrap(),
            "--verify",
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    let j_v1: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(j_v1["format"], "kset-v1 (GRMKSET1)");
    assert_eq!(j_v1["verified"], false);
    assert!(j_v1["header_sha256"].is_null());
    assert_eq!(j_v1["fnv64"], want_fnv);
    assert_eq!(j_v1["keys_sha256"], j_gr["keys_sha256"]);

    // 8) Flipping 1 byte in .grk, .ktab, or .prof.bin MUST fail `inspect --verify`
    //    and `convert-v1` with non-zero exit!
    for (src, ext) in [
        (&gr_out, "grk"),
        (&counted_tab, "ktab"),
        (&prof_path, "prof.bin"),
    ] {
        let mut flipped = fs::read(src).unwrap();
        let last = flipped.len() - 1;
        flipped[last] ^= 0x01;
        let bad = dir.path().join(format!("flipped.{ext}"));
        fs::write(&bad, &flipped).unwrap();

        let insp = grm_bin()
            .args([
                "kcount",
                "inspect",
                "--input",
                bad.to_str().unwrap(),
                "--verify",
            ])
            .output()
            .unwrap();
        assert!(
            !insp.status.success(),
            "inspect --verify must reject 1-byte-flipped .{ext}"
        );

        if ext != "prof.bin" {
            let conv = grm_bin()
                .args([
                    "kcount",
                    "convert-v1",
                    "--input",
                    bad.to_str().unwrap(),
                    "--out",
                    dir.path().join(format!("bad_{ext}.kset")).to_str().unwrap(),
                ])
                .output()
                .unwrap();
            assert!(
                !conv.status.success(),
                "convert-v1 must reject 1-byte-flipped .{ext}"
            );
        }
    }
}
