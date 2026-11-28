//! Checked 2-bit SPSS simplitig archive integrity, damaged-file stress,
//! fixed-size regression, byte-for-byte `.kset` v1 / `.grk` round-trip, and
//! end-to-end CLI tests for `grm::kcount::spss` and
//! `grm kcount archive` / `expand` / `inspect`.

use flate2::Compression;
use flate2::write::GzEncoder;
use grm::dna::KmerSpace;
use grm::kcount::format::{
    DEFAULT_SAMPLE_NAME, GrKsetReader, InputModality, KtabReader, write_gr_kset, write_ktab,
};
use grm::kcount::key::{KmerSpace64, KmerSpace128};
use grm::kcount::spss::{
    MAGIC_SPSS1, archive_from_fasta, archive_from_table, expand_spss_keys, expand_to_gr,
    expand_to_kset_v1, expand_to_ktab, read_spss_header, verify_spss, write_spss_from_mixed,
};
use grm::kset;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use ring::digest::{Context as Sha256Context, SHA256};
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use std::process::Command;

const SPSS_DIGEST_START: usize = 48;
const SPSS_DIGEST_END: usize = 80;

fn rehash_spss(bytes: &mut [u8]) {
    let mut sha = Sha256Context::new(&SHA256);
    sha.update(&bytes[..SPSS_DIGEST_START]);
    sha.update(&bytes[SPSS_DIGEST_END..]);
    let digest = sha.finish();
    bytes[SPSS_DIGEST_START..SPSS_DIGEST_END].copy_from_slice(digest.as_ref());
}

fn rand_dna(rng: &mut ChaCha8Rng, len: usize) -> Vec<u8> {
    (0..len).map(|_| b"ACGT"[rng.random_range(0..4)]).collect()
}

fn grm_bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_grm"))
}

fn write_and_verify_spss_fails(dir: &Path, tag: &str, bytes: &[u8], expect_header_fail: bool) {
    let p = dir.join(format!("{tag}.spss"));
    fs::write(&p, bytes).unwrap();
    if expect_header_fail {
        assert!(
            read_spss_header(&p).is_err(),
            "expected read_spss_header (inspect) to reject damaged spss {tag}"
        );
    }
    assert!(
        expand_spss_keys::<u64>(&p).is_err(),
        "expected expand_spss_keys (expand) to reject damaged spss {tag}"
    );
    assert!(
        verify_spss(&p).is_err(),
        "expected verify_spss (inspect --verify) to reject damaged spss {tag}"
    );
}

// ---------------------------------------------------------------------------
// 1. Every single-byte flip anywhere in .spss is caught by expand and verify
// ---------------------------------------------------------------------------

#[test]
fn spss_rejects_every_single_byte_flip() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(3141);
    let seq = rand_dna(&mut rng, 1_200);
    let fna = dir.path().join("mini.fna");
    fs::write(&fna, [b">c1\n".as_slice(), &seq, b"\n"].concat()).unwrap();

    let spss_path = dir.path().join("intact.spss");
    archive_from_fasta(&fna, &spss_path, 31, "sample", InputModality::Assembly).unwrap();
    let spss_bytes = fs::read(&spss_path).unwrap();
    assert!(spss_bytes.len() > 350);

    let flip_path = dir.path().join("flip.spss");
    for pos in 0..spss_bytes.len() {
        let mut bad = spss_bytes.clone();
        bad[pos] ^= 0x01;
        fs::write(&flip_path, &bad).unwrap();
        assert!(
            expand_spss_keys::<u64>(&flip_path).is_err(),
            "spss byte flip at offset {pos} of {} escaped expand",
            spss_bytes.len()
        );
        assert!(
            verify_spss(&flip_path).is_err(),
            "spss byte flip at offset {pos} of {} escaped verify",
            spss_bytes.len()
        );
    }
}

// ---------------------------------------------------------------------------
// 2. Damaged .spss stress test (>600 damaged files: truncations, flips,
//    wrong magic, huge u32::MAX / u64::MAX / 1<<31 / 1<<63 lengths raw and
//    rehashed) -> 0 panics, 0 OOM aborts, 100% clean Err in expand & inspect
// ---------------------------------------------------------------------------

#[test]
fn damaged_spss_stress_zero_panics_in_expand_and_inspect() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(20261127);

    let fna = dir.path().join("base.fna");
    let mut seq = rand_dna(&mut rng, 400);
    seq.extend_from_slice(b"NNNN");
    seq.extend_from_slice(&rand_dna(&mut rng, 350));
    fs::write(&fna, [b">c1\n".as_slice(), &seq, b"\n"].concat()).unwrap();

    let spss_path = dir.path().join("base.spss");
    archive_from_fasta(&fna, &spss_path, 31, "sample", InputModality::Assembly).unwrap();
    let base = fs::read(&spss_path).unwrap();

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

    // All byte-by-byte truncations up to base.len() (280+ cuts) must fail BOTH
    // read_spss_header (inspect) AND expand_spss_keys (expand).
    for cut in 0..base.len() {
        write_and_verify_spss_fails(dir.path(), &format!("trunc_{cut}"), &base[..cut], true);
        total_cases += 1;
    }

    // Wrong magic (40) + random single/multi-byte corruptions (200).
    for i in 0..240 {
        let mut bad = base.clone();
        let header_fail = i < 40;
        if i < 40 {
            bad[i % 8] ^= 0x5a;
        } else {
            let n_flips = rng.random_range(1..=4);
            for _ in 0..n_flips {
                let at = rng.random_range(0..bad.len());
                bad[at] ^= rng.random_range(1u8..=255);
            }
        }
        write_and_verify_spss_fails(dir.path(), &format!("mut_{i}"), &bad, header_fail);
        total_cases += 1;
    }

    // Huge lengths (both raw and rehashed) at every length/count header field:
    // name_len [12..16], n_kmers [16..24], n_simplitigs [24..32],
    // total_bases [32..40], and lengths_byte_len [80 + name_len .. +8].
    // EVERY one of these MUST fail in BOTH inspect (read_spss_header) AND expand!
    let len_field_off = 80 + "sample".len();

    for (idx, &v32) in huge_u32s.iter().enumerate() {
        for rehash in [false, true] {
            let mut bad = base.clone();
            bad[12..16].copy_from_slice(&v32.to_le_bytes());
            if rehash {
                rehash_spss(&mut bad);
            }
            write_and_verify_spss_fails(dir.path(), &format!("hname_{idx}_{rehash}"), &bad, true);
            total_cases += 1;
        }
    }

    for (idx, &v64) in huge_u64s.iter().enumerate() {
        for rehash in [false, true] {
            // Huge n_simplitigs alone (exact Finding 1 repro: 1 << 63!)
            let mut bad = base.clone();
            bad[24..32].copy_from_slice(&v64.to_le_bytes());
            if rehash {
                rehash_spss(&mut bad);
            }
            write_and_verify_spss_fails(dir.path(), &format!("hsimp_{idx}_{rehash}"), &bad, true);
            total_cases += 1;

            // Huge n_simplitigs + matching n_kmers + wrapped total_bases
            let mut bad = base.clone();
            let wrapped_bases = v64.wrapping_add(30u64.wrapping_mul(v64));
            bad[16..24].copy_from_slice(&v64.to_le_bytes());
            bad[24..32].copy_from_slice(&v64.to_le_bytes());
            bad[32..40].copy_from_slice(&wrapped_bases.to_le_bytes());
            if rehash {
                rehash_spss(&mut bad);
            }
            write_and_verify_spss_fails(
                dir.path(),
                &format!("hsimp_wrap_{idx}_{rehash}"),
                &bad,
                true,
            );
            total_cases += 1;

            // Huge n_kmers
            let mut bad = base.clone();
            bad[16..24].copy_from_slice(&v64.to_le_bytes());
            if rehash {
                rehash_spss(&mut bad);
            }
            write_and_verify_spss_fails(dir.path(), &format!("hnk_{idx}_{rehash}"), &bad, true);
            total_cases += 1;

            // Huge total_bases
            let mut bad = base.clone();
            bad[32..40].copy_from_slice(&v64.to_le_bytes());
            if rehash {
                rehash_spss(&mut bad);
            }
            write_and_verify_spss_fails(dir.path(), &format!("hbases_{idx}_{rehash}"), &bad, true);
            total_cases += 1;

            // Huge lengths_byte_len
            let mut bad = base.clone();
            bad[len_field_off..len_field_off + 8].copy_from_slice(&v64.to_le_bytes());
            if rehash {
                rehash_spss(&mut bad);
            }
            write_and_verify_spss_fails(dir.path(), &format!("hlenb_{idx}_{rehash}"), &bad, true);
            total_cases += 1;
        }
    }

    assert!(
        total_cases > 550,
        "expected > 550 damaged SPSS cases, ran {total_cases}"
    );
}

// ---------------------------------------------------------------------------
// 3. Rehashed malformed .spss files checking EVERY header/body/uniqueness gate
//    (including exact duplicate k-mer gate at spss.rs bucket uniqueness check!)
// ---------------------------------------------------------------------------

#[test]
fn rehashed_spss_validation_rejects_malformed_fields_and_duplicates() {
    let dir = tempfile::tempdir().unwrap();

    // Build a valid base .spss at k = 5 (2 simplitigs of length 5 -> 10 bases ->
    // 3 packed bytes with 4 unused trailing padding bits in the final byte).
    let space5 = KmerSpace64::new(5);
    let c1 = space5.canonical("AAAAA").unwrap();
    let c2 = space5.canonical("CCCCC").unwrap();
    let mut mixed2 = vec![space5.mix(c1), space5.mix(c2)];
    mixed2.sort_unstable();
    let ok_path = dir.path().join("ok5.spss");
    write_spss_from_mixed(&ok_path, "ok", 5, InputModality::Assembly, &mixed2).unwrap();
    let ok_bytes = fs::read(&ok_path).unwrap();

    type Mutator = (&'static str, bool, fn(&mut Vec<u8>));

    let cases: &[Mutator] = &[
        ("bad_magic", true, |b| b[0] = b'X'),
        ("bad_k_low", true, |b| b[8] = 4),
        ("bad_k_high", true, |b| b[8] = 65),
        ("bad_modality", true, |b| b[9] = 3),
        ("nonzero_reserved", true, |b| b[10] = 1),
        ("bad_utf8_name", true, |b| b[80] = 0xff),
        ("trailing_byte", true, |b| b.push(0)),
        ("nonzero_padding_bits", false, |b| {
            let last = b.len() - 1;
            b[last] |= 0xf0;
        }),
        ("wrong_fnv64", false, |b| {
            b[40] ^= 0x01;
        }),
    ];
    for &(label, header_fail, mutate) in cases {
        let mut bad = ok_bytes.clone();
        mutate(&mut bad);
        rehash_spss(&mut bad);
        write_and_verify_spss_fails(dir.path(), &format!("rehash_{label}"), &bad, header_fail);
    }

    // Craft a rehashed .spss with a malformed length stream:
    // n_simplitigs = 2, n_kmers = 3, total_bases = 3 + 4*2 = 11, lengths_byte_len = 2,
    // ULEB extras [2, 0] -> sum_bases = 7 + 5 = 12 != 11 (exceeds total_bases!).
    {
        let mut bad = Vec::new();
        bad.extend_from_slice(MAGIC_SPSS1);
        bad.push(5); // k
        bad.push(0); // modality
        bad.extend_from_slice(&[0, 0]);
        bad.extend_from_slice(&2u32.to_le_bytes()); // name_len
        bad.extend_from_slice(&3u64.to_le_bytes()); // n_kmers
        bad.extend_from_slice(&2u64.to_le_bytes()); // n_simplitigs
        bad.extend_from_slice(&11u64.to_le_bytes()); // total_bases
        bad.extend_from_slice(&0u64.to_le_bytes()); // fnv
        bad.extend_from_slice(&[0u8; 32]); // digest placeholder
        bad.extend_from_slice(b"ok");
        bad.extend_from_slice(&2u64.to_le_bytes()); // lengths_byte_len
        bad.extend_from_slice(&[2u8, 0u8]); // extras -> 7 + 5 = 12 != 11
        bad.extend_from_slice(&[0u8; 3]); // ceil(11/4) = 3 packed bytes
        rehash_spss(&mut bad);
        write_and_verify_spss_fails(dir.path(), "rehash_sum_bases_mismatch", &bad, false);
    }

    // Craft a rehashed .spss with leftover unconsumed length byte:
    // n_simplitigs = 1, n_kmers = 1, total_bases = 5, lengths_byte_len = 2,
    // len_bytes = [0, 0] (reads first 0, leaves second byte unconsumed!).
    {
        let mut bad = Vec::new();
        bad.extend_from_slice(MAGIC_SPSS1);
        bad.push(5);
        bad.push(0);
        bad.extend_from_slice(&[0, 0]);
        bad.extend_from_slice(&2u32.to_le_bytes());
        bad.extend_from_slice(&1u64.to_le_bytes());
        bad.extend_from_slice(&1u64.to_le_bytes());
        bad.extend_from_slice(&5u64.to_le_bytes());
        bad.extend_from_slice(&0u64.to_le_bytes());
        bad.extend_from_slice(&[0u8; 32]);
        bad.extend_from_slice(b"ok");
        bad.extend_from_slice(&2u64.to_le_bytes());
        bad.extend_from_slice(&[0u8, 0u8]);
        bad.extend_from_slice(&[0u8; 2]); // ceil(5/4) = 2
        rehash_spss(&mut bad);
        write_and_verify_spss_fails(dir.path(), "rehash_leftover_len_byte", &bad, false);
    }

    // Craft a rehashed .spss with a truncated ULEB128 varint (continuation bit set).
    {
        let mut bad = Vec::new();
        bad.extend_from_slice(MAGIC_SPSS1);
        bad.push(5);
        bad.push(0);
        bad.extend_from_slice(&[0, 0]);
        bad.extend_from_slice(&2u32.to_le_bytes());
        bad.extend_from_slice(&1u64.to_le_bytes());
        bad.extend_from_slice(&1u64.to_le_bytes());
        bad.extend_from_slice(&5u64.to_le_bytes());
        bad.extend_from_slice(&0u64.to_le_bytes());
        bad.extend_from_slice(&[0u8; 32]);
        bad.extend_from_slice(b"ok");
        bad.extend_from_slice(&1u64.to_le_bytes());
        bad.extend_from_slice(&[0x80u8]); // truncated ULEB
        bad.extend_from_slice(&[0u8; 2]);
        rehash_spss(&mut bad);
        write_and_verify_spss_fails(dir.path(), "rehash_trunc_uleb", &bad, false);
    }

    // Craft a rehashed .spss encoding TWO IDENTICAL k-mers ("AAAAA" twice) WITH
    // matching n_kmers = 2, matching n_simplitigs = 2, matching total_bases = 10,
    // valid zero-padded packed bases, matching FNV-1a64 over [mixed(AAAAA), mixed(AAAAA)],
    // AND valid SHA-256! ONLY the duplicate-k-mer check inside expand_spss_keys
    // catches this file (kills the spss.rs bucket uniqueness fault).
    {
        let m_a = space5.mix(c1);
        let fnv_dup = u64::from_str_radix(&kset::digest(&[m_a, m_a]), 16).unwrap();

        let mut dup = Vec::new();
        dup.extend_from_slice(MAGIC_SPSS1);
        dup.push(5); // k
        dup.push(InputModality::Assembly.as_byte());
        dup.extend_from_slice(&[0, 0]);
        dup.extend_from_slice(&2u32.to_le_bytes()); // name "ok"
        dup.extend_from_slice(&2u64.to_le_bytes()); // n_kmers = 2
        dup.extend_from_slice(&2u64.to_le_bytes()); // n_simplitigs = 2
        dup.extend_from_slice(&10u64.to_le_bytes()); // total_bases = 2 + 4*2 = 10
        dup.extend_from_slice(&fnv_dup.to_le_bytes()); // matching fnv64!
        dup.extend_from_slice(&[0u8; 32]); // digest slot
        dup.extend_from_slice(b"ok");
        dup.extend_from_slice(&2u64.to_le_bytes()); // lengths_byte_len = 2
        dup.extend_from_slice(&[0u8, 0u8]); // two length-5 simplitigs
        dup.extend_from_slice(&[0u8, 0u8, 0u8]); // 10 'A' bases (all 0 bits) + 0 padding
        rehash_spss(&mut dup);

        let p = dir.path().join("rehash_duplicate_kmer.spss");
        fs::write(&p, &dup).unwrap();
        // Header is structurally valid; expand MUST reject on duplicate k-mer!
        assert!(read_spss_header(&p).is_ok());
        let err = expand_spss_keys::<u64>(&p).unwrap_err();
        assert!(
            format!("{err:#}").contains("duplicate k-mer"),
            "expected duplicate k-mer error, got: {err:#}"
        );
        assert!(verify_spss(&p).is_err());
    }
}

// ---------------------------------------------------------------------------
// 4. Fixed-size regression on deterministic synthetic 500-contig input
//    (kills spss.rs right/left greedy simplitig extension size faults!)
// ---------------------------------------------------------------------------

#[test]
fn fixed_size_regression_on_synthetic_500_contig_genome() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(500_031);
    let k = 31usize;
    let space = KmerSpace::new(k);

    // Deterministic synthetic 500-contig assembly (~459 kb) with shared
    // repeats and N runs to exercise both left/right greedy simplitig extension
    // and branching splits.
    let shared_rep = rand_dna(&mut rng, 180);
    let fna = dir.path().join("synth_500_contigs.fna");
    {
        let mut w = std::io::BufWriter::new(File::create(&fna).unwrap());
        for c in 0..500usize {
            let mut contig = rand_dna(&mut rng, 900);
            if c % 25 == 0 {
                contig.extend_from_slice(&shared_rep);
            }
            if c % 10 == 0 {
                contig.extend_from_slice(b"NNNNN");
                contig.extend_from_slice(&rand_dna(&mut rng, 100));
            }
            writeln!(w, ">contig_{c}").unwrap();
            w.write_all(&contig).unwrap();
            writeln!(w).unwrap();
        }
        w.flush().unwrap();
    }

    let (values, _) = kset::count_file(&fna, &space, 1).unwrap();
    let v1_path = dir.path().join("synth_500.kset");
    kset::write(&v1_path, DEFAULT_SAMPLE_NAME, &space, &values).unwrap();
    let v1_bytes = fs::metadata(&v1_path).unwrap().len();

    let grk_path = dir.path().join("synth_500.grk");
    write_gr_kset(
        &grk_path,
        DEFAULT_SAMPLE_NAME,
        k,
        InputModality::Assembly,
        &values,
    )
    .unwrap();
    let grk_bytes = fs::metadata(&grk_path).unwrap().len();

    let spss_path = dir.path().join("synth_500.spss");
    let arch = archive_from_table(
        &v1_path,
        &spss_path,
        DEFAULT_SAMPLE_NAME,
        Some(InputModality::Assembly),
        false,
    )
    .unwrap();

    // Exact pinned counts and byte size: any regression in right/left greedy
    // extension increases n_simplitigs and archive_bytes and fails immediately.
    assert_eq!(arch.n_kmers, 439_226, "500-contig n_kmers changed");
    assert_eq!(arch.n_simplitigs, 550, "500-contig n_simplitigs changed");
    assert_eq!(arch.total_bases, 455_726, "500-contig total_bases changed");
    assert_eq!(
        arch.archive_bytes, 115_076,
        "500-contig .spss byte size changed"
    );
    assert_eq!(grk_bytes, 2_471_268, "500-contig .grk byte size changed");
    assert!(arch.bytes_per_kmer < 0.27);
    assert!(arch.archive_bytes * 20 < v1_bytes);

    // Exact byte-for-byte round-trip back to .kset v1.
    let rebuilt_v1 = dir.path().join("rebuilt_500.kset");
    let exp = expand_to_kset_v1(&spss_path, &rebuilt_v1).unwrap();
    assert_eq!(exp.fnv64, arch.fnv64);
    assert_eq!(fs::read(&rebuilt_v1).unwrap(), fs::read(&v1_path).unwrap());
}

// ---------------------------------------------------------------------------
// 5. Exact byte-for-byte .kset v1 -> .spss -> .kset v1 (k = 5, 21, 31, 32)
//    AND .grk -> .spss -> .grk (k = 33, 45, 61, 64) + direct FASTA/FASTQ.gz
// ---------------------------------------------------------------------------

#[test]
fn round_trip_byte_for_byte_kset_v1_u64_and_grk_u128() {
    let dir = tempfile::tempdir().unwrap();
    let mut rng = ChaCha8Rng::seed_from_u64(6432);

    // Multi-contig synthetic genome with repeats, N runs, and lowercase bases.
    let c1 = rand_dna(&mut rng, 6_000);
    let rep = c1[1_000..1_200].to_vec();
    let mut c2 = rand_dna(&mut rng, 3_000);
    c2.extend_from_slice(b"NNNN");
    c2.extend_from_slice(&rep);
    c2.extend_from_slice(b"acgtacgt");
    let fna = dir.path().join("multi.fna");
    fs::write(
        &fna,
        [b">c1\n".as_slice(), &c1, b"\n>c2\n", &c2, b"\n"].concat(),
    )
    .unwrap();
    let fna_gz = dir.path().join("multi.fna.gz");
    {
        let mut enc = GzEncoder::new(File::create(&fna_gz).unwrap(), Compression::fast());
        enc.write_all(&fs::read(&fna).unwrap()).unwrap();
        enc.finish().unwrap();
    }

    // 1) u64 widths: .kset v1 -> .spss -> .kset v1 byte-for-byte at k = 5, 21, 31, 32.
    for k in [5usize, 21, 31, 32] {
        let space = KmerSpace::new(k);
        let (values, _) = kset::count_file(&fna, &space, 1).unwrap();
        let v1_orig = dir.path().join(format!("orig_k{k}.kset"));
        kset::write(&v1_orig, "sample", &space, &values).unwrap();
        let want_v1 = fs::read(&v1_orig).unwrap();

        let spss = dir.path().join(format!("arch_k{k}.spss"));
        let arch = archive_from_table(&v1_orig, &spss, "sample", None, false).unwrap();
        assert_eq!(arch.k, k);
        assert_eq!(arch.n_kmers as usize, values.len());
        assert_eq!(arch.fnv64, kset::digest(&values));

        let v1_back = dir.path().join(format!("back_k{k}.kset"));
        let exp = expand_to_kset_v1(&spss, &v1_back).unwrap();
        assert_eq!(exp.fnv64, arch.fnv64);
        assert_eq!(
            fs::read(&v1_back).unwrap(),
            want_v1,
            ".kset v1 round-trip differed at k={k}"
        );

        // Direct FASTA.gz archive produces the exact same .spss bytes!
        let spss_gz = dir.path().join(format!("arch_gz_k{k}.spss"));
        archive_from_fasta(&fna_gz, &spss_gz, k, "sample", InputModality::Unknown).unwrap();
        assert_eq!(
            fs::read(&spss_gz).unwrap(),
            fs::read(&spss).unwrap(),
            "FASTA.gz archive differed from table archive at k={k}"
        );
    }

    // 2) u128 widths: .grk -> .spss -> .grk byte-for-byte at k = 33, 45, 61, 64.
    for (idx, k) in [33usize, 45, 61, 64].into_iter().enumerate() {
        let space = KmerSpace128::new(k);
        let modality = if idx % 2 == 0 {
            InputModality::Assembly
        } else {
            InputModality::Reads
        };
        let mut mixed: Vec<u128> = Vec::new();
        space.for_each_canonical(&c1, |x| mixed.push(space.mix(x)));
        space.for_each_canonical(&c2, |x| mixed.push(space.mix(x)));
        mixed.sort_unstable();
        mixed.dedup();

        let grk_orig = dir.path().join(format!("orig_k{k}.grk"));
        write_gr_kset(&grk_orig, "sample", k, modality, &mixed).unwrap();
        let want_grk = fs::read(&grk_orig).unwrap();

        let spss = dir.path().join(format!("arch_k{k}.spss"));
        let arch = archive_from_table(&grk_orig, &spss, "sample", None, false).unwrap();
        assert_eq!(arch.k, k);
        assert_eq!(arch.modality, modality);
        assert_eq!(arch.n_kmers as usize, mixed.len());

        let grk_back = dir.path().join(format!("back_k{k}.grk"));
        let exp = expand_to_gr(&spss, &grk_back).unwrap();
        assert_eq!(exp.fnv64, arch.fnv64);
        assert_eq!(exp.modality, Some(modality));
        assert_eq!(
            fs::read(&grk_back).unwrap(),
            want_grk,
            ".grk round-trip differed at k={k}"
        );

        // Expanding k > 32 to .kset v1 must fail cleanly.
        assert!(expand_to_kset_v1(&spss, &dir.path().join(format!("bad_{k}.kset"))).is_err());

        // Expanding to presence-only .ktab preserves all keys and modality with count = 1.
        let tab_back = dir.path().join(format!("back_k{k}.ktab"));
        expand_to_ktab(&spss, &tab_back).unwrap();
        let rtab = KtabReader::open_verified(&tab_back).unwrap();
        assert!(rtab.presence_only);
        assert_eq!(rtab.modality, modality);
        let tab_entries = rtab.read_all::<u128>().unwrap();
        assert_eq!(tab_entries.len(), mixed.len());
        assert!(tab_entries.iter().all(|&(_, c)| c == 1));
        assert_eq!(
            tab_entries
                .into_iter()
                .map(|(key, _)| key)
                .collect::<Vec<_>>(),
            mixed
        );

        // Direct FASTA.gz archive matches table archive byte-for-byte when modality matches.
        let spss_gz = dir.path().join(format!("arch_gz_k{k}.spss"));
        archive_from_fasta(&fna_gz, &spss_gz, k, "sample", modality).unwrap();
        assert_eq!(
            fs::read(&spss_gz).unwrap(),
            fs::read(&spss).unwrap(),
            "FASTA.gz u128 archive differed at k={k}"
        );
    }

    // 3) Verify k = 61 high-bit / low-bit pairs round-trip cleanly without collision.
    let s61 = KmerSpace128::new(61);
    let prefix32 = "ACGT".repeat(8);
    let suffix32 = "TGCA".repeat(8);
    let mut pair_keys = vec![
        s61.mix(
            s61.canonical(&format!("{prefix32}{}", "A".repeat(29)))
                .unwrap(),
        ),
        s61.mix(
            s61.canonical(&format!("{prefix32}{}C", "A".repeat(28)))
                .unwrap(),
        ),
        s61.mix(
            s61.canonical(&format!("{}{suffix32}", "C".repeat(29)))
                .unwrap(),
        ),
        s61.mix(
            s61.canonical(&format!("{}G{suffix32}", "C".repeat(28)))
                .unwrap(),
        ),
    ];
    pair_keys.sort_unstable();
    let gr61 = dir.path().join("pairs61.grk");
    write_gr_kset(&gr61, "sample", 61, InputModality::Reads, &pair_keys).unwrap();
    let sp61 = dir.path().join("pairs61.spss");
    archive_from_table(&gr61, &sp61, "sample", None, false).unwrap();
    let gr61_back = dir.path().join("pairs61_back.grk");
    expand_to_gr(&sp61, &gr61_back).unwrap();
    assert_eq!(
        GrKsetReader::open_verified(&gr61_back)
            .unwrap()
            .read_all::<u128>()
            .unwrap(),
        pair_keys
    );
}

// ---------------------------------------------------------------------------
// 6. End-to-end CLI tests:
//    - archive refuses counted .ktab without --drop-counts
//    - archive succeeds with --drop-counts and expands to counts = 1
//    - archive defaults --name to "sample" (zero filename stem leak)
//    - expand enforces output extensions
//    - inspect rejects crafted n_simplitigs = 1<<63 without --verify
//    - inspect without --verify -> verified: false; with --verify -> true
//    - inspect --verify rejects 1-byte-flipped .spss
// ---------------------------------------------------------------------------

#[test]
fn cli_archive_expand_inspect_drop_counts_and_no_filename_leak() {
    let dir = tempfile::tempdir().unwrap();
    let k = 31usize;
    let space = KmerSpace::new(k);
    let mut rng = ChaCha8Rng::seed_from_u64(777);
    let seq = rand_dna(&mut rng, 2_500);

    let secret_stem = "arbitrary_archive_stem_77777";
    let fna_in = dir.path().join(format!("{secret_stem}.fna"));
    fs::write(&fna_in, [b">c1\n".as_slice(), &seq, b"\n"].concat()).unwrap();

    let (values, _) = kset::count_file(&fna_in, &space, 1).unwrap();
    let want_fnv = kset::digest(&values);

    // Write a counted .ktab with counts > 1 AND a presence-only .ktab.
    let counted_entries: Vec<(u64, u32)> = values
        .iter()
        .enumerate()
        .map(|(i, &key)| (key, if i % 3 == 0 { 12 } else { 1 }))
        .collect();
    let counted_tab = dir.path().join(format!("{secret_stem}.ktab"));
    write_ktab(
        &counted_tab,
        secret_stem,
        k,
        InputModality::Reads,
        1,
        0,
        &counted_entries,
    )
    .unwrap();

    let presence_entries: Vec<(u64, u32)> = values.iter().map(|&key| (key, 1)).collect();
    let presence_tab = dir.path().join("presence_only.ktab");
    write_ktab(
        &presence_tab,
        secret_stem,
        k,
        InputModality::Assembly,
        1,
        0,
        &presence_entries,
    )
    .unwrap();

    // 1) Archive on counted .ktab WITHOUT --drop-counts MUST FAIL cleanly!
    let refused_spss = dir.path().join("refused.spss");
    let out = grm_bin()
        .args([
            "kcount",
            "archive",
            "--input",
            counted_tab.to_str().unwrap(),
            "--out",
            refused_spss.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "archive must refuse counted .ktab without --drop-counts"
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("--drop-counts"),
        "expected error mentioning --drop-counts, got: {stderr}"
    );
    assert!(!refused_spss.exists());

    // 2) Archive on counted .ktab WITH --drop-counts succeeds, defaults sample
    //    to "sample" (NO filename leak!), and reports counts_dropped: true.
    let dropped_spss = dir.path().join("dropped.spss");
    let out = grm_bin()
        .args([
            "kcount",
            "archive",
            "--input",
            counted_tab.to_str().unwrap(),
            "--out",
            dropped_spss.to_str().unwrap(),
            "--drop-counts",
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "archive --drop-counts failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    assert!(
        !stdout.contains(secret_stem),
        "archive stdout leaked input filename stem"
    );
    let arch_j: serde_json::Value = serde_json::from_str(&stdout).unwrap();
    assert_eq!(arch_j["sample"], DEFAULT_SAMPLE_NAME);
    assert_eq!(arch_j["k"], k);
    assert_eq!(arch_j["modality"], "reads");
    assert_eq!(arch_j["n_kmers"], values.len() as u64);
    assert_eq!(arch_j["counts_dropped"], true);
    assert_eq!(arch_j["fnv64"], want_fnv);
    assert_eq!(arch_j["header_sha256"].as_str().unwrap().len(), 64);

    // Expanding dropped.spss to .ktab produces presence-only counts == 1.
    let exp_tab = dir.path().join("expanded.ktab");
    let out = grm_bin()
        .args([
            "kcount",
            "expand",
            "--input",
            dropped_spss.to_str().unwrap(),
            "--out",
            exp_tab.to_str().unwrap(),
            "--format",
            "ktab",
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "expand --format ktab failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let exp_j: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(exp_j["format"], "GRMKTB1");
    assert_eq!(exp_j["sample"], DEFAULT_SAMPLE_NAME);
    assert_eq!(exp_j["fnv64"], want_fnv);
    let exp_reader = KtabReader::open_verified(&exp_tab).unwrap();
    assert!(exp_reader.presence_only);
    assert_eq!(exp_reader.max_count_observed, 1);
    assert_eq!(exp_reader.name, DEFAULT_SAMPLE_NAME);
    let exp_entries = exp_reader.read_all::<u64>().unwrap();
    assert!(exp_entries.iter().all(|&(_, c)| c == 1));
    assert_eq!(
        exp_entries
            .into_iter()
            .map(|(key, _)| key)
            .collect::<Vec<_>>(),
        values
    );

    // 3) Archive on presence-only .ktab succeeds WITHOUT --drop-counts.
    let pres_spss = dir.path().join("presence.spss");
    let out = grm_bin()
        .args([
            "kcount",
            "archive",
            "--input",
            presence_tab.to_str().unwrap(),
            "--out",
            pres_spss.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    let pres_j: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(pres_j["counts_dropped"], false);
    assert_eq!(pres_j["sample"], DEFAULT_SAMPLE_NAME);

    // 4) Archive directly from FASTA WITHOUT --name defaults to "sample" and
    //    never leaks the input filename stem into header or stdout.
    let fa_spss = dir.path().join("from_fa.spss");
    let out = grm_bin()
        .args([
            "kcount",
            "archive",
            "--input",
            fna_in.to_str().unwrap(),
            "--out",
            fa_spss.to_str().unwrap(),
            "--k",
            "31",
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    let fa_stdout = String::from_utf8(out.stdout).unwrap();
    assert!(!fa_stdout.contains(secret_stem));
    let fa_j: serde_json::Value = serde_json::from_str(&fa_stdout).unwrap();
    assert_eq!(fa_j["sample"], DEFAULT_SAMPLE_NAME);
    assert_eq!(
        read_spss_header(&fa_spss).unwrap().name,
        DEFAULT_SAMPLE_NAME
    );

    // 5) Expand enforces output extension agreement.
    let bad_ext = grm_bin()
        .args([
            "kcount",
            "expand",
            "--input",
            fa_spss.to_str().unwrap(),
            "--out",
            dir.path().join("wrong.ktab").to_str().unwrap(),
            "--format",
            "kset-v1",
        ])
        .output()
        .unwrap();
    assert!(!bad_ext.status.success());
    assert!(String::from_utf8_lossy(&bad_ext.stderr).contains("disagrees"));

    // 6) inspect WITHOUT --verify reports verified: false and no input path;
    //    inspect WITH --verify reports verified: true, fnv64, keys_sha256, sketch.
    let out = grm_bin()
        .args([
            "kcount",
            "inspect",
            "--input",
            dropped_spss.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    let insp_no: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(insp_no["format"], "spss-archive (GRMSPS1)");
    assert_eq!(insp_no["k"], k);
    assert_eq!(insp_no["sample"], DEFAULT_SAMPLE_NAME);
    assert_eq!(insp_no["verified"], false);
    assert!(insp_no["fnv64"].is_null());
    assert!(insp_no["keys_sha256"].is_null());
    assert_eq!(insp_no["header_sha256"], arch_j["header_sha256"]);
    assert!(insp_no["n_simplitigs"].as_u64().unwrap() >= 1);
    assert!(insp_no.get("path").is_none());

    let insp_file = dir.path().join("insp_verified.json");
    let out = grm_bin()
        .args([
            "kcount",
            "inspect",
            "--input",
            dropped_spss.to_str().unwrap(),
            "--verify",
            "--sketch",
            "8",
            "--out",
            insp_file.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(out.status.success());
    let insp_yes: serde_json::Value =
        serde_json::from_slice(&fs::read(&insp_file).unwrap()).unwrap();
    assert_eq!(insp_yes["format"], "spss-archive (GRMSPS1)");
    assert_eq!(insp_yes["verified"], true);
    assert_eq!(insp_yes["fnv64"], want_fnv);
    assert_eq!(insp_yes["keys_sha256"].as_str().unwrap().len(), 64);
    assert_eq!(insp_yes["sketch_hex"].as_array().unwrap().len(), 8);

    // 7) 1-byte-flipped .spss MUST fail `inspect --verify` AND `expand`!
    let mut flipped = fs::read(&dropped_spss).unwrap();
    let last = flipped.len() - 1;
    flipped[last] ^= 0x01;
    let flipped_spss = dir.path().join("flipped.spss");
    fs::write(&flipped_spss, &flipped).unwrap();

    let out = grm_bin()
        .args([
            "kcount",
            "inspect",
            "--input",
            flipped_spss.to_str().unwrap(),
            "--verify",
        ])
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "inspect --verify must reject 1-byte-flipped .spss"
    );

    let out = grm_bin()
        .args([
            "kcount",
            "expand",
            "--input",
            flipped_spss.to_str().unwrap(),
            "--out",
            dir.path().join("should_not_exist.kset").to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "expand must reject 1-byte-flipped .spss"
    );

    // 8) Crafted .spss with n_simplitigs = 1 << 63 MUST fail `inspect` even
    //    WITHOUT `--verify` (never returning exit code 0)!
    let mut crafted = fs::read(&dropped_spss).unwrap();
    crafted[24..32].copy_from_slice(&(1u64 << 63).to_le_bytes());
    let crafted_spss = dir.path().join("crafted_huge.spss");
    fs::write(&crafted_spss, &crafted).unwrap();

    let out = grm_bin()
        .args([
            "kcount",
            "inspect",
            "--input",
            crafted_spss.to_str().unwrap(),
        ])
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "inspect without --verify must reject n_simplitigs = 1<<63"
    );
}
