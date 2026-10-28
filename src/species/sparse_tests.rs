use super::*;
use crate::confine::suffixed;
use std::process::Command;

const SLUG: &str = "sparse-fixture";
const FILES: [&str; 5] = [
    "matrix.json",
    "patterns.bin",
    "support.bin",
    "kmer_offsets.bin",
    "kmers.bin",
];

fn sequence(g: usize) -> String {
    let mut seq = String::from("ACGTTGCAAGGCTTACCGATAGCTAGGCTTAACGGTACCATGCAGTTCAGGA");
    seq.push_str("NNNN");
    seq.extend((0..96usize).map(|i| {
        let x = g.wrapping_mul(29) + i.wrapping_mul(17) + (i / 11) * 23;
        b"ACGT"[(x + (x >> 5)) % 4] as char
    }));
    seq
}

fn fixture(root: &Path, n: usize, empty: bool) -> Layout {
    let layout = Layout::new(root);
    let ids: Vec<String> = (0..n).map(|g| format!("g{g:04}")).collect();
    std::fs::create_dir_all(layout.fasta_dir()).unwrap();
    std::fs::create_dir_all(layout.species_dir(SLUG)).unwrap();
    std::fs::write(
        layout.species_dir(SLUG).join("genomes.txt"),
        ids.join("\n") + "\n",
    )
    .unwrap();
    for (g, id) in ids.iter().enumerate() {
        let seq = if empty { "NNNN".into() } else { sequence(g) };
        std::fs::write(
            fasta_path(&layout.fasta_dir(), id),
            format!(">contig\n{seq}\n"),
        )
        .unwrap();
    }
    layout
}

fn sparse(layout: &Layout, k: usize) -> Result<matrix::MatrixSummary> {
    build_matrix_disk_sparse_with(layout, SLUG, k, None, 1, 64, 0, |_| Ok(()))
}

fn sparse_hook(layout: &Layout, k: usize, stop: &str) -> Result<matrix::MatrixSummary> {
    build_matrix_disk_sparse_with(layout, SLUG, k, None, 1, 64, 0, |event| {
        if event == stop {
            anyhow::bail!("injected stop at {event}")
        }
        Ok(())
    })
}

fn run_dir(layout: &Layout) -> PathBuf {
    let root = layout.species_dir(SLUG).join("matrix-partials");
    let dirs: Vec<PathBuf> = std::fs::read_dir(root)
        .unwrap()
        .map(|e| e.unwrap().path())
        .collect();
    assert_eq!(dirs.len(), 1, "expected one active sparse run");
    dirs[0].clone()
}

fn read_run(layout: &Layout) -> SparseRun {
    serde_json::from_slice(&std::fs::read(run_dir(layout).join("run.json")).unwrap()).unwrap()
}

fn matrix_files(layout: &Layout) -> Vec<Vec<u8>> {
    FILES
        .iter()
        .map(|name| std::fs::read(layout.matrix_dir(SLUG).join(name)).unwrap())
        .collect()
}

fn assert_same_matrix(a: &Layout, b: &Layout) {
    assert_eq!(matrix_files(a), matrix_files(b));
    let loaded = matrix::FeatureMatrix::load(&b.matrix_dir(SLUG)).unwrap();
    assert_eq!(loaded.summary.genomes.len(), loaded.patterns.cols);
}

fn ksets_remain(layout: &Layout, mut ids: std::ops::Range<usize>) -> bool {
    let dir = layout.kset_dir(SLUG);
    ids.any(|g| dir.join(format!("g{g:04}.kset")).exists())
}

#[test]
fn disk_sparse_matches_stock_and_structure_after_kset_deletion() {
    let dir = tempfile::tempdir().unwrap();
    let stock = fixture(&dir.path().join("stock"), 70, false);
    let sparse_layout = fixture(&dir.path().join("sparse"), 70, false);
    count(&stock, SLUG, 7, None).unwrap();
    build_matrix(&stock, SLUG, 1).unwrap();
    let stock_structure = structure(&stock, SLUG).unwrap();
    sparse(&sparse_layout, 7).unwrap();
    assert_same_matrix(&stock, &sparse_layout);
    assert!(!ksets_remain(&sparse_layout, 0..70));
    assert!(ksets_remain(&stock, 0..70));
    let counted = load_counted(&sparse_layout, SLUG).unwrap();
    assert!(counted.pending.is_empty());
    assert!(counted.genomes.iter().all(|g| g.sketch.is_some()));
    let sparse_structure = structure(&sparse_layout, SLUG).unwrap();
    assert_eq!(sparse_structure.genomes, stock_structure.genomes);
    assert_eq!(sparse_structure.kset_digests, stock_structure.kset_digests);
    assert_eq!(
        sparse_structure.lineages.labels,
        stock_structure.lineages.labels
    );
    assert!(
        !sparse_layout
            .species_dir(SLUG)
            .join("matrix-partials")
            .exists()
    );
}

#[test]
fn sparse_preserves_preexisting_ksets() {
    let dir = tempfile::tempdir().unwrap();
    let stock = fixture(&dir.path().join("stock"), 64, false);
    let existing = fixture(&dir.path().join("existing"), 64, false);
    count(&stock, SLUG, 7, None).unwrap();
    build_matrix(&stock, SLUG, 1).unwrap();
    count(&existing, SLUG, 7, None).unwrap();
    sparse(&existing, 7).unwrap();
    assert_same_matrix(&stock, &existing);
    assert!(ksets_remain(&existing, 0..64));
}

#[test]
fn sparse_resume_binds_genome_order_and_kset_digests() {
    let dir = tempfile::tempdir().unwrap();
    let layout = fixture(dir.path(), 64, false);
    assert!(sparse_hook(&layout, 7, "after_partial:0").is_err());
    let counted_path = layout.species_dir(SLUG).join("counted.json");
    let original = std::fs::read(&counted_path).unwrap();
    let mut counted = load_counted(&layout, SLUG).unwrap();
    counted.genomes.swap(0, 1);
    std::fs::write(&counted_path, serde_json::to_vec(&counted).unwrap()).unwrap();
    let err = sparse(&layout, 7).unwrap_err();
    assert!(
        format!("{err:#}").contains("no longer reconstructs"),
        "{err:#}"
    );
    assert!(ksets_remain(&layout, 0..64));

    std::fs::write(&counted_path, &original).unwrap();
    let mut counted = load_counted(&layout, SLUG).unwrap();
    counted.genomes[0].kset_digest = Some("wrong".into());
    std::fs::write(&counted_path, serde_json::to_vec(&counted).unwrap()).unwrap();
    let err = sparse(&layout, 7).unwrap_err();
    assert!(
        format!("{err:#}").contains("not the expected wrong"),
        "{err:#}"
    );
    assert!(ksets_remain(&layout, 0..64));
    std::fs::write(&counted_path, original).unwrap();
    sparse(&layout, 7).unwrap();
    assert!(!ksets_remain(&layout, 0..64));
    matrix::FeatureMatrix::load(&layout.matrix_dir(SLUG)).unwrap();
}

#[test]
fn sparse_resume_cleans_orphans_and_refuses_unknown_files() {
    let dir = tempfile::tempdir().unwrap();
    let layout = fixture(dir.path(), 65, false);
    assert!(sparse_hook(&layout, 7, "after_partial:0").is_err());
    let counted = load_counted(&layout, SLUG).unwrap();
    assert_eq!(counted.genomes.len(), 64);
    assert_eq!(counted.pending, ["g0064"]);
    assert!(count(&layout, SLUG, 7, None).is_err());
    assert!(ensure_selection_counted(&layout, SLUG).is_err());
    let run = run_dir(&layout);
    let partial = partial_path(&run, 0);
    std::fs::remove_file(matrix::receipt_path(&partial)).unwrap();
    std::fs::write(suffixed(&partial, ".part"), "orphan").unwrap();
    std::fs::write(suffixed(&matrix::receipt_path(&partial), ".part"), "orphan").unwrap();
    std::fs::write(run.join("run.json.part"), "orphan").unwrap();
    let unknown = run.join("do-not-delete");
    std::fs::write(&unknown, "keep me").unwrap();
    let err = sparse(&layout, 7).unwrap_err();
    assert!(format!("{err:#}").contains("unknown file"), "{err:#}");
    assert_eq!(std::fs::read(&unknown).unwrap(), b"keep me");
    assert!(ksets_remain(&layout, 0..64));
    std::fs::remove_file(unknown).unwrap();
    sparse(&layout, 7).unwrap();
    assert!(!ksets_remain(&layout, 0..65));
    matrix::FeatureMatrix::load(&layout.matrix_dir(SLUG)).unwrap();
}

#[test]
fn sparse_rejects_a_stale_folded_block_without_recreating_it() {
    let dir = tempfile::tempdir().unwrap();
    let layout = fixture(dir.path(), 65, false);
    assert!(sparse_hook(&layout, 7, "after_fold:0").is_err());
    assert_eq!(read_run(&layout).folded_columns, 64);
    assert!(!ksets_remain(&layout, 0..64));
    let receipt = matrix::receipt_path(&partial_path(&run_dir(&layout), 0));
    let original = std::fs::read(&receipt).unwrap();
    let mut forged: matrix::PartialReceipt = serde_json::from_slice(&original).unwrap();
    forged.columns[0].kset_digest = "stale".into();
    std::fs::write(&receipt, serde_json::to_vec(&forged).unwrap()).unwrap();
    let err = sparse(&layout, 7).unwrap_err();
    assert!(
        format!("{err:#}").contains("header and receipt disagree"),
        "{err:#}"
    );
    assert_eq!(read_run(&layout).folded_columns, 64);
    assert!(!layout.matrix_dir(SLUG).join("matrix.json").exists());
    std::fs::write(receipt, original).unwrap();
    sparse(&layout, 7).unwrap();
    assert!(!ksets_remain(&layout, 0..65));
}

#[test]
fn sparse_rewrites_a_mixed_k_all_empty_active_block() {
    let dir = tempfile::tempdir().unwrap();
    let stock = fixture(&dir.path().join("stock"), 65, true);
    let layout = fixture(&dir.path().join("sparse"), 65, true);
    count(&stock, SLUG, 31, None).unwrap();
    build_matrix(&stock, SLUG, 1).unwrap();
    assert!(sparse_hook(&layout, 31, "after_partial:0").is_err());
    let kdir = layout.kset_dir(SLUG);
    let paths: Vec<PathBuf> = (0..64)
        .map(|g| kdir.join(format!("g{g:04}.kset")))
        .collect();
    let original: Vec<Vec<u8>> = paths.iter().map(|p| std::fs::read(p).unwrap()).collect();
    let space32 = KmerSpace::new(32);
    for (g, path) in paths.iter().enumerate() {
        kset::write(path, &format!("g{g:04}"), &space32, &[]).unwrap();
    }
    let counted = load_counted(&layout, SLUG).unwrap();
    assert!(counted.genomes.iter().all(|g| g.n_kmers == 0));
    assert!(
        counted
            .genomes
            .iter()
            .all(|g| g.kset_digest == counted.genomes[0].kset_digest)
    );
    let columns: Vec<matrix::Column> = counted.genomes[..64]
        .iter()
        .map(|g| matrix::Column {
            name: g.genome_id.clone(),
            kset_digest: g.kset_digest.clone().unwrap(),
        })
        .collect();
    let partial = partial_path(&run_dir(&layout), 0);
    matrix::write_partial(&paths, &columns, 0, &partial, 8, 0).unwrap();
    for (path, bytes) in paths.iter().zip(original) {
        std::fs::write(path, bytes).unwrap();
    }
    sparse(&layout, 31).unwrap();
    assert_same_matrix(&stock, &layout);
    assert!(!ksets_remain(&layout, 0..65));
}

#[test]
fn sparse_recreates_old_missing_sketches_even_when_matrix_is_current() {
    let dir = tempfile::tempdir().unwrap();
    let layout = fixture(dir.path(), 64, false);
    count(&layout, SLUG, 7, None).unwrap();
    build_matrix(&layout, SLUG, 1).unwrap();
    let expected = matrix_files(&layout);
    let counted_path = layout.species_dir(SLUG).join("counted.json");
    let mut counted = load_counted(&layout, SLUG).unwrap();
    for g in &mut counted.genomes {
        g.sketch = None;
    }
    std::fs::write(&counted_path, serde_json::to_vec(&counted).unwrap()).unwrap();
    for g in 0..64 {
        std::fs::remove_file(layout.kset_dir(SLUG).join(format!("g{g:04}.kset"))).unwrap();
    }
    sparse(&layout, 7).unwrap();
    assert_eq!(matrix_files(&layout), expected);
    assert!(
        load_counted(&layout, SLUG)
            .unwrap()
            .genomes
            .iter()
            .all(|g| g.sketch.is_some())
    );
    assert!(!ksets_remain(&layout, 0..64));
    structure(&layout, SLUG).unwrap();
}

#[test]
fn structure_prefers_real_ksets_and_validates_cache_only_when_absent() {
    let dir = tempfile::tempdir().unwrap();
    let layout = fixture(dir.path(), 4, false);
    count(&layout, SLUG, 7, None).unwrap();
    let path = layout.species_dir(SLUG).join("counted.json");
    let mut counted = load_counted(&layout, SLUG).unwrap();
    counted.genomes[0].sketch = Some(vec![]);
    std::fs::write(&path, serde_json::to_vec(&counted).unwrap()).unwrap();
    // The wrong-length cache is irrelevant while the real kset exists.
    structure(&layout, SLUG).unwrap();
    std::fs::remove_file(layout.kset_dir(SLUG).join("g0000.kset")).unwrap();
    let err = structure(&layout, SLUG).unwrap_err();
    assert!(
        format!("{err:#}").contains("wrong sketch length"),
        "{err:#}"
    );
}

#[test]
fn sparse_cleans_an_initial_manifest_temp_but_not_unowned_files() {
    let dir = tempfile::tempdir().unwrap();
    let layout = fixture(dir.path(), 64, false);
    let config = SparseConfig {
        k: 7,
        genomes: (0..64).map(|g| format!("g{g:04}")).collect(),
        missing_fasta: vec![],
        catalog_version: String::new(),
        block_width: 64,
    };
    let run = sparse_run_dir(&layout, SLUG, &config).unwrap();
    std::fs::create_dir_all(&run).unwrap();
    std::fs::write(run.join("run.json.part"), "interrupted create").unwrap();
    sparse(&layout, 7).unwrap();
    assert!(!run.exists());
    matrix::FeatureMatrix::load(&layout.matrix_dir(SLUG)).unwrap();
}

#[test]
fn mtb_probe_panels_do_not_read_deleted_ksets() {
    let missing = vec![PathBuf::from("definitely-absent.kset")];
    let assays = crate::probes::for_species("Mycobacterium tuberculosis");
    assert!(!assays.is_empty());
    for assay in assays {
        let err = crate::probes::build_probes(assay, &missing, &KmerSpace::new(31)).unwrap_err();
        assert!(
            !err.is::<crate::probes::KsetsUnreadable>(),
            "{} unexpectedly required ksets: {err:#}",
            assay.id
        );
    }
}

#[test]
fn sparse_reserve_and_partial_failure_happen_before_kset_deletion() {
    let dir = tempfile::tempdir().unwrap();
    let refused = fixture(&dir.path().join("refused"), 64, false);
    assert_eq!(DISK_SPARSE_RESERVE, 30 * 1024 * 1024 * 1024);
    let err = build_matrix_disk_sparse_with(&refused, SLUG, 7, None, 1, 64, u64::MAX, |_| Ok(()))
        .unwrap_err();
    assert!(format!("{err:#}").contains("reserve"), "{err:#}");
    assert!(!ksets_remain(&refused, 0..64));

    let layout = fixture(&dir.path().join("failure"), 64, false);
    assert!(sparse_hook(&layout, 7, "before_count:0").is_err());
    let blocker = suffixed(&partial_path(&run_dir(&layout), 0), ".part");
    std::fs::create_dir(&blocker).unwrap();
    // A blocker already there is refused before anything is counted.
    let err = sparse(&layout, 7).unwrap_err();
    assert!(format!("{err:#}").contains("non-regular"), "{err:#}");
    assert!(!layout.kset_dir(SLUG).join("g0000.kset").exists());
    std::fs::remove_dir(&blocker).unwrap();
    // One that appears while the block is counted stops the partial write,
    // and every kset it would have covered stays.
    let err = build_matrix_disk_sparse_with(&layout, SLUG, 7, None, 1, 64, 0, |event| {
        if event == "before_count:0" {
            std::fs::create_dir(&blocker)?;
        }
        Ok(())
    })
    .unwrap_err();
    assert!(format!("{err:#}").contains("non-regular"), "{err:#}");
    assert!((0..64).all(|g| ksets_remain(&layout, g..g + 1)));
    std::fs::remove_dir(blocker).unwrap();
    sparse(&layout, 7).unwrap();
    assert!(!ksets_remain(&layout, 0..64));
}

#[test]
fn sparse_keeps_partials_when_matrix_save_fails() {
    let dir = tempfile::tempdir().unwrap();
    let layout = fixture(dir.path(), 64, false);
    let blocker = layout.matrix_dir(SLUG).join("kmers.bin");
    // A blocker already there is refused before anything is counted.
    std::fs::create_dir_all(&blocker).unwrap();
    let err = sparse(&layout, 7).unwrap_err();
    assert!(format!("{err:#}").contains("kmers.bin"), "{err:#}");
    assert!(!layout.kset_dir(SLUG).exists());
    std::fs::remove_dir(&blocker).unwrap();
    // One that appears after the last fold stops the save itself.
    let err = build_matrix_disk_sparse_with(&layout, SLUG, 7, None, 1, 64, 0, |event| {
        if event == "after_fold:0" {
            std::fs::create_dir_all(&blocker)?;
        }
        Ok(())
    })
    .unwrap_err();
    assert!(format!("{err:#}").contains("kmers.bin"), "{err:#}");
    assert!(!layout.matrix_dir(SLUG).join("matrix.json").exists());
    assert!(!ksets_remain(&layout, 0..64));
    let run = run_dir(&layout);
    assert!(partial_path(&run, 0).exists());
    assert!(matrix::receipt_path(&partial_path(&run, 0)).exists());
    assert!(read_run(&layout).saved.is_none());
    std::fs::remove_dir(layout.matrix_dir(SLUG).join("kmers.bin")).unwrap();
    sparse(&layout, 7).unwrap();
    matrix::FeatureMatrix::load(&layout.matrix_dir(SLUG)).unwrap();
    assert!(!run.exists());
}

#[test]
fn sparse_saved_digest_gates_cleanup_and_interrupted_cleanup_resumes() {
    let dir = tempfile::tempdir().unwrap();
    let layout = fixture(dir.path(), 64, false);
    assert!(sparse_hook(&layout, 7, "after_saved").is_err());
    assert!(read_run(&layout).saved.is_some());
    let expected = matrix_files(&layout);
    let run = run_dir(&layout);
    assert!(partial_path(&run, 0).exists());
    let support = layout.matrix_dir(SLUG).join("support.bin");
    let original = std::fs::read(&support).unwrap();
    std::fs::write(&support, [original.as_slice(), &[0]].concat()).unwrap();
    let err = sparse(&layout, 7).unwrap_err();
    assert!(
        format!("{err:#}").contains("does not match the digest"),
        "{err:#}"
    );
    assert!(partial_path(&run, 0).exists());
    std::fs::write(&support, original).unwrap();
    sparse(&layout, 7).unwrap();
    assert_eq!(matrix_files(&layout), expected);
    assert!(!run.exists());
}

#[test]
fn sparse_resume_after_kill_before_kset_delete() {
    let dir = tempfile::tempdir().unwrap();
    let stock = fixture(&dir.path().join("stock"), 65, false);
    let sparse_layout = fixture(&dir.path().join("sparse"), 65, false);
    count(&stock, SLUG, 7, None).unwrap();
    build_matrix(&stock, SLUG, 1).unwrap();
    let marker = dir.path().join("partial-verified");
    let exe = std::env::current_exe().unwrap();
    let mut child = Command::new(exe)
        .arg("--exact")
        .arg("species::sparse_tests::sparse_kill_after_partial_probe")
        .arg("--ignored")
        .arg("--nocapture")
        .arg("--test-threads=1")
        .env("GRM_SPARSE_KILL_ROOT", &sparse_layout.root)
        .env("GRM_SPARSE_KILL_MARKER", &marker)
        .env("RAYON_NUM_THREADS", "4")
        .spawn()
        .unwrap();
    let started = std::time::Instant::now();
    while !marker.exists() {
        if let Some(status) = child.try_wait().unwrap() {
            panic!("kill probe exited before its marker: {status}");
        }
        if started.elapsed() > std::time::Duration::from_secs(30) {
            child.kill().unwrap();
            child.wait().unwrap();
            panic!("kill probe did not reach its marker");
        }
        std::thread::sleep(std::time::Duration::from_millis(25));
    }
    child.kill().unwrap();
    let status = child.wait().unwrap();
    assert!(!status.success());
    let run = read_run(&sparse_layout);
    assert_eq!(run.folded_columns, 0);
    assert!(run.active.is_some());
    assert!(ksets_remain(&sparse_layout, 0..64));
    let partial = partial_path(&run_dir(&sparse_layout), 0);
    assert!(partial.exists());
    assert!(matrix::receipt_path(&partial).exists());
    sparse(&sparse_layout, 7).unwrap();
    assert_same_matrix(&stock, &sparse_layout);
    assert!(!ksets_remain(&sparse_layout, 0..65));
}

#[test]
#[ignore = "run by sparse_resume_after_kill_before_kset_delete"]
fn sparse_kill_after_partial_probe() {
    let layout = Layout::new(PathBuf::from(
        std::env::var_os("GRM_SPARSE_KILL_ROOT").unwrap(),
    ));
    let marker = PathBuf::from(std::env::var_os("GRM_SPARSE_KILL_MARKER").unwrap());
    build_matrix_disk_sparse_with(&layout, SLUG, 7, None, 1, 64, 0, |event| {
        if event == "after_partial:0" {
            std::fs::write(&marker, b"verified")?;
            std::thread::sleep(std::time::Duration::from_secs(120));
        }
        Ok(())
    })
    .unwrap();
}

#[test]
fn sparse_merge_under_an_open_file_limit_that_stops_stock_build() {
    let dir = tempfile::tempdir().unwrap();
    let stock = fixture(&dir.path().join("stock"), 200, false);
    let sparse_layout = fixture(&dir.path().join("sparse"), 200, false);
    count(&stock, SLUG, 7, None).unwrap();
    build_matrix(&stock, SLUG, 1).unwrap();
    let exe = std::env::current_exe().unwrap();
    let status = Command::new("/bin/sh")
        .arg("-c")
        .arg("ulimit -Sn 96 && exec \"$0\" --exact species::sparse_tests::sparse_open_file_limit_probe --ignored --nocapture --test-threads=1")
        .arg(exe)
        .env("GRM_SPARSE_FD_ROOT", &sparse_layout.root)
        .env("RAYON_NUM_THREADS", "4")
        .status()
        .unwrap();
    assert!(status.success());
    assert_same_matrix(&stock, &sparse_layout);
    assert!(!ksets_remain(&sparse_layout, 0..200));
}

#[test]
#[ignore = "run by sparse_merge_under_an_open_file_limit_that_stops_stock_build"]
fn sparse_open_file_limit_probe() {
    let root = std::env::var_os("GRM_SPARSE_FD_ROOT").unwrap();
    let layout = Layout::new(PathBuf::from(root));
    sparse(&layout, 7).unwrap();
}

fn sparse_at(layout: &Layout, min_support: u32, width: usize) -> Result<matrix::MatrixSummary> {
    build_matrix_disk_sparse_with(layout, SLUG, 7, None, min_support, width, 0, |_| Ok(()))
}

fn stock_matrix(root: &Path, n: usize, min_support: u32) -> Layout {
    let stock = fixture(root, n, false);
    count(&stock, SLUG, 7, None).unwrap();
    build_matrix(&stock, SLUG, min_support).unwrap();
    stock
}

/// Every entry below `dir`: file bytes, and links as their targets.
fn tree(dir: &Path) -> Vec<(PathBuf, Vec<u8>)> {
    let mut out = Vec::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(d) = stack.pop() {
        for entry in std::fs::read_dir(&d).unwrap() {
            let path = entry.unwrap().path();
            let m = std::fs::symlink_metadata(&path).unwrap();
            if m.file_type().is_symlink() {
                let target = std::fs::read_link(&path).unwrap();
                out.push((path, target.into_os_string().into_encoded_bytes()));
            } else if m.is_dir() {
                out.push((path.clone(), b"dir".to_vec()));
                stack.push(path);
            } else {
                out.push((path.clone(), std::fs::read(&path).unwrap()));
            }
        }
    }
    out.sort();
    out
}

fn kset_bytes(layout: &Layout, ids: std::ops::Range<usize>) -> Vec<Vec<u8>> {
    ids.map(|g| std::fs::read(layout.kset_dir(SLUG).join(format!("g{g:04}.kset"))).unwrap())
        .collect()
}

#[test]
fn sparse_refuses_links_where_it_writes_before_any_work() {
    use std::os::unix::fs::symlink;
    let dir = tempfile::tempdir().unwrap();
    type Plant = fn(&Layout, &Path, &Path);
    let cases: [(&str, Plant); 8] = [
        ("species-dir", |l, out, _| {
            let real = out.join("species-real");
            std::fs::rename(l.species_dir(SLUG), &real).unwrap();
            symlink(&real, l.species_dir(SLUG)).unwrap();
        }),
        ("ksets-dir", |l, out, _| {
            symlink(out, l.kset_dir(SLUG)).unwrap()
        }),
        ("matrix-dir", |l, out, _| {
            std::fs::write(out.join("matrix.json"), "frozen reference").unwrap();
            symlink(out, l.matrix_dir(SLUG)).unwrap();
        }),
        ("matrix-file", |l, _, victim| {
            std::fs::create_dir_all(l.matrix_dir(SLUG)).unwrap();
            symlink(victim, l.matrix_dir(SLUG).join("patterns.bin")).unwrap();
        }),
        ("matrix-partials", |l, out, _| {
            symlink(out, l.species_dir(SLUG).join("matrix-partials")).unwrap()
        }),
        ("counted-json", |l, _, victim| {
            symlink(victim, l.species_dir(SLUG).join("counted.json")).unwrap()
        }),
        ("sidecar", |l, _, victim| {
            std::fs::create_dir_all(l.kset_dir(SLUG)).unwrap();
            symlink(victim, l.kset_dir(SLUG).join("g0003.json")).unwrap();
        }),
        ("kset-temporary", |l, _, victim| {
            std::fs::create_dir_all(l.kset_dir(SLUG)).unwrap();
            symlink(victim, l.kset_dir(SLUG).join("g0003.kset.part")).unwrap();
        }),
    ];
    for (name, plant) in cases {
        let layout = fixture(&dir.path().join(name), 70, false);
        let out = dir.path().join(format!("{name}-outside"));
        std::fs::create_dir_all(&out).unwrap();
        let victim = out.join("victim");
        std::fs::write(&victim, "a shared file outside the lane").unwrap();
        plant(&layout, &out, &victim);
        let before = tree(&out);
        let err = sparse(&layout, 7).unwrap_err();
        assert!(format!("{err:#}").contains("symlink"), "{name}: {err:#}");
        assert_eq!(tree(&out), before, "{name}");
    }
}

#[test]
fn sparse_refuses_links_in_its_run_directory_on_resume() {
    use std::os::unix::fs::symlink;
    let dir = tempfile::tempdir().unwrap();
    let stock = stock_matrix(&dir.path().join("stock"), 70, 1);
    let out = dir.path().join("outside");
    std::fs::create_dir_all(&out).unwrap();
    let victim = out.join("victim");
    std::fs::write(&victim, "outside").unwrap();

    let layout = fixture(&dir.path().join("part"), 70, false);
    assert!(sparse_hook(&layout, 7, "after_fold:0").is_err());
    let run = run_dir(&layout);
    symlink(&victim, run.join("run.json.part")).unwrap();
    let err = sparse(&layout, 7).unwrap_err();
    assert!(format!("{err:#}").contains("non-regular entry"), "{err:#}");
    assert_eq!(std::fs::read(&victim).unwrap(), b"outside");
    std::fs::remove_file(run.join("run.json.part")).unwrap();
    sparse(&layout, 7).unwrap();
    assert_same_matrix(&stock, &layout);

    let layout = fixture(&dir.path().join("run-dir"), 70, false);
    assert!(sparse_hook(&layout, 7, "after_fold:0").is_err());
    let run = run_dir(&layout);
    let moved = out.join("run-elsewhere");
    std::fs::rename(&run, &moved).unwrap();
    symlink(&moved, &run).unwrap();
    let before = tree(&out);
    let err = sparse(&layout, 7).unwrap_err();
    assert!(format!("{err:#}").contains("unknown entry"), "{err:#}");
    assert_eq!(tree(&out), before);
    std::fs::remove_file(&run).unwrap();
    std::fs::rename(&moved, &run).unwrap();
    sparse(&layout, 7).unwrap();
    assert_same_matrix(&stock, &layout);
}

#[test]
fn a_linked_data_root_is_resolved_once_and_a_linked_kept_kset_is_only_read() {
    let dir = tempfile::tempdir().unwrap();
    let stock = stock_matrix(&dir.path().join("stock"), 70, 1);
    let real = fixture(&dir.path().join("real"), 70, false);
    let alias = dir.path().join("alias");
    std::os::unix::fs::symlink(&real.root, &alias).unwrap();
    // A kset in a shared store, linked in without its sidecar: it is counted
    // again in memory, found equal, and left as it is.
    let shared = dir.path().join("shared-g0003.kset");
    std::fs::copy(stock.kset_dir(SLUG).join("g0003.kset"), &shared).unwrap();
    let shared_bytes = std::fs::read(&shared).unwrap();
    std::fs::create_dir_all(real.kset_dir(SLUG)).unwrap();
    let linked = real.kset_dir(SLUG).join("g0003.kset");
    std::os::unix::fs::symlink(&shared, &linked).unwrap();
    sparse(&Layout::new(&alias), 7).unwrap();
    assert_same_matrix(&stock, &real);
    assert!(linked.is_symlink());
    assert_eq!(std::fs::read(&shared).unwrap(), shared_bytes);
    assert!(!ksets_remain(&real, 0..3) && !ksets_remain(&real, 4..70));
}

#[test]
fn sparse_never_rewrites_a_kset_it_did_not_create() {
    let dir = tempfile::tempdir().unwrap();
    let stock = stock_matrix(&dir.path().join("stock"), 70, 1);
    // Ksets left at another k, without the counted.json that went with them.
    let layout = fixture(&dir.path().join("sparse"), 70, false);
    count(&layout, SLUG, 9, None).unwrap();
    std::fs::remove_file(layout.species_dir(SLUG).join("counted.json")).unwrap();
    let before = kset_bytes(&layout, 0..70);
    let err = sparse(&layout, 7).unwrap_err();
    assert!(
        format!("{err:#}").contains("never rewrites a kset it did not create"),
        "{err:#}"
    );
    assert_eq!(kset_bytes(&layout, 0..70), before);
    assert!(!layout.species_dir(SLUG).join("counted.json").exists());
    assert!(!layout.species_dir(SLUG).join("matrix-partials").exists());

    // At the right k but with other contents and no sidecar, a kept kset is
    // refused once counted, and still not rewritten.
    for g in 0..70 {
        std::fs::remove_file(layout.kset_dir(SLUG).join(format!("g{g:04}.kset"))).unwrap();
    }
    let space = KmerSpace::new(7);
    let fasta = fasta_path(&layout.fasta_dir(), "g0003");
    let (mut values, _) = kset::count_file(&fasta, &space, 1).unwrap();
    values.pop();
    let odd = layout.kset_dir(SLUG).join("g0003.kset");
    kset::write(&odd, "g0003", &space, &values).unwrap();
    std::fs::remove_file(layout.kset_dir(SLUG).join("g0003.json")).unwrap();
    let odd_bytes = std::fs::read(&odd).unwrap();
    let err = sparse(&layout, 7).unwrap_err();
    assert!(
        format!("{err:#}").contains("differs from a fresh count"),
        "{err:#}"
    );
    assert_eq!(std::fs::read(&odd).unwrap(), odd_bytes);
    // The operator removes it; the run then counts it like the others.
    std::fs::remove_file(&odd).unwrap();
    sparse(&layout, 7).unwrap();
    assert_same_matrix(&stock, &layout);
}

#[test]
fn sparse_refuses_a_forged_partial_whose_receipt_matches_it() {
    let dir = tempfile::tempdir().unwrap();
    let stock = stock_matrix(&dir.path().join("stock"), 70, 1);
    let layout = fixture(&dir.path().join("sparse"), 70, false);
    assert!(sparse_hook(&layout, 7, "after_fold:0").is_err());
    let counted = load_counted(&layout, SLUG).unwrap();
    let (columns, _) = block_expectations(&counted, 0, 64);
    // The forger's block: g0000 lost one k-mer.
    let forge = dir.path().join("forge");
    std::fs::create_dir_all(&forge).unwrap();
    let space = KmerSpace::new(7);
    let paths: Vec<PathBuf> = (0..64)
        .map(|g| {
            let id = format!("g{g:04}");
            let fasta = fasta_path(&layout.fasta_dir(), &id);
            let (mut values, _) = kset::count_file(&fasta, &space, 1).unwrap();
            if g == 0 {
                values.remove(values.len() / 2);
            }
            let path = forge.join(format!("{id}.kset"));
            kset::write(&path, &id, &space, &values).unwrap();
            path
        })
        .collect();
    let mut forged_columns = columns.clone();
    let dropped = KsetReader::open(&paths[0]).unwrap().read_all().unwrap();
    forged_columns[0].kset_digest = kset::digest(&dropped);
    let forged = forge.join("block-00000.partial");
    matrix::write_partial(&paths, &forged_columns, 0, &forged, 8, 0).unwrap();
    // The header claims g0000's real digest (column 0's field follows its
    // 4-byte name length and name), and the receipt matches the edited file.
    let mut bytes = std::fs::read(&forged).unwrap();
    let at = 48 + 4 + "g0000".len();
    let real = u64::from_str_radix(&columns[0].kset_digest, 16).unwrap();
    bytes[at..at + 8].copy_from_slice(&real.to_le_bytes());
    let mut receipt: matrix::PartialReceipt =
        serde_json::from_slice(&std::fs::read(matrix::receipt_path(&forged)).unwrap()).unwrap();
    receipt.columns = columns.clone();
    receipt.bytes = bytes.len() as u64;
    receipt.sha256 = ring::digest::digest(&ring::digest::SHA256, &bytes)
        .as_ref()
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect();
    let partial = partial_path(&run_dir(&layout), 0);
    let originals = [
        std::fs::read(&partial).unwrap(),
        std::fs::read(matrix::receipt_path(&partial)).unwrap(),
    ];
    std::fs::write(&partial, &bytes).unwrap();
    std::fs::write(
        matrix::receipt_path(&partial),
        serde_json::to_vec(&receipt).unwrap(),
    )
    .unwrap();
    let err = sparse(&layout, 7).unwrap_err();
    assert!(
        format!("{err:#}").contains("does not hold that genome's kset"),
        "{err:#}"
    );
    assert!(!layout.matrix_dir(SLUG).join("matrix.json").exists());
    std::fs::write(&partial, &originals[0]).unwrap();
    std::fs::write(matrix::receipt_path(&partial), &originals[1]).unwrap();
    sparse(&layout, 7).unwrap();
    assert_same_matrix(&stock, &layout);
}

#[test]
fn sparse_refuses_to_orphan_another_runs_partials() {
    let dir = tempfile::tempdir().unwrap();
    let stock = stock_matrix(&dir.path().join("stock"), 130, 1);
    let layout = fixture(&dir.path().join("sparse"), 130, false);
    assert!(sparse_hook(&layout, 7, "after_fold:0").is_err());
    let counted_path = layout.species_dir(SLUG).join("counted.json");
    let counted_before = std::fs::read(&counted_path).unwrap();
    let old_run = run_dir(&layout);
    let err = sparse_at(&layout, 1, 128).unwrap_err();
    assert!(
        format!("{err:#}").contains("unfinished disk-sparse run"),
        "{err:#}"
    );
    assert_eq!(std::fs::read(&counted_path).unwrap(), counted_before);
    assert_eq!(run_dir(&layout), old_run);
    // A changed selection is refused the same way, so reverting it resumes.
    let genomes = layout.species_dir(SLUG).join("genomes.txt");
    let original = std::fs::read_to_string(&genomes).unwrap();
    let mut ids: Vec<&str> = original.lines().collect();
    ids.swap(128, 129);
    std::fs::write(&genomes, ids.join("\n") + "\n").unwrap();
    let err = sparse(&layout, 7).unwrap_err();
    assert!(
        format!("{err:#}").contains("unfinished disk-sparse run"),
        "{err:#}"
    );
    assert_eq!(std::fs::read(&counted_path).unwrap(), counted_before);
    std::fs::write(&genomes, original).unwrap();
    sparse(&layout, 7).unwrap();
    assert_same_matrix(&stock, &layout);
    assert!(!layout.species_dir(SLUG).join("matrix-partials").exists());
}

#[test]
fn an_interrupted_cleanups_empty_run_directory_is_removed() {
    let dir = tempfile::tempdir().unwrap();
    let layout = fixture(dir.path(), 64, false);
    sparse(&layout, 7).unwrap();
    let expected = matrix_files(&layout);
    let config = SparseConfig {
        k: 7,
        genomes: (0..64).map(|g| format!("g{g:04}")).collect(),
        missing_fasta: vec![],
        catalog_version: String::new(),
        block_width: 64,
    };
    // What a cleanup leaves when it stops after removing run.json.
    std::fs::create_dir_all(sparse_run_dir(&layout, SLUG, &config).unwrap()).unwrap();
    sparse(&layout, 7).unwrap();
    assert!(!layout.species_dir(SLUG).join("matrix-partials").exists());
    assert_eq!(matrix_files(&layout), expected);
}

#[test]
fn sparse_rebuilds_at_another_min_support_and_adopts_a_complete_stock_count() {
    let dir = tempfile::tempdir().unwrap();
    let stock = stock_matrix(&dir.path().join("stock"), 130, 2);
    // After a finished run deleted its ksets, another min_support recounts
    // them block by block, and each must reproduce its recorded digest.
    let layout = fixture(&dir.path().join("rebuilt"), 130, false);
    sparse_at(&layout, 1, 64).unwrap();
    sparse_at(&layout, 2, 64).unwrap();
    assert_same_matrix(&stock, &layout);
    assert!(!ksets_remain(&layout, 0..130));
    assert!(!layout.species_dir(SLUG).join("matrix-partials").exists());

    // A complete stock count wider than one block is adopted, ksets intact.
    let adopted = fixture(&dir.path().join("adopted"), 130, false);
    count(&adopted, SLUG, 7, None).unwrap();
    let before = kset_bytes(&adopted, 0..130);
    sparse_at(&adopted, 2, 64).unwrap();
    assert_same_matrix(&stock, &adopted);
    assert_eq!(kset_bytes(&adopted, 0..130), before);

    // A matrix already saved at one min_support is finished at it first.
    let saved = fixture(&dir.path().join("saved"), 70, false);
    assert!(sparse_hook(&saved, 7, "after_saved").is_err());
    let err = sparse_at(&saved, 2, 64).unwrap_err();
    assert!(format!("{err:#}").contains("--min-support 1"), "{err:#}");
    sparse_at(&saved, 1, 64).unwrap();
    assert_eq!(sparse_at(&saved, 2, 64).unwrap().min_support, 2);
    assert!(!saved.species_dir(SLUG).join("matrix-partials").exists());
}

#[test]
fn parallel_writes_share_one_reserve() {
    fn free(_: &Path) -> Result<u64> {
        Ok(1_000)
    }
    let reserve = Reserve {
        bytes: 400,
        claimed: Mutex::new(0),
        free,
    };
    let here = Path::new(".");
    let a = reserve.claim(here, 300).unwrap();
    // Free space alone would pass this; the first claim already holds 300.
    let err = reserve.claim(here, 301).err().unwrap();
    assert!(
        format!("{err:#}").contains("other writers have claimed"),
        "{err:#}"
    );
    let b = reserve.claim(here, 300).unwrap();
    drop(a);
    let c = reserve.claim(here, 300).unwrap();
    drop((b, c));
    assert_eq!(*reserve.claimed.lock().unwrap(), 0);
    assert!(Reserve::new(0).claim(here, u64::MAX).is_ok());
}

#[test]
fn the_block_preflight_bounds_the_ksets_it_writes() {
    let dir = tempfile::tempdir().unwrap();
    let layout = fixture(dir.path(), 64, false);
    let ids: Vec<String> = (0..64).map(|g| format!("g{g:04}")).collect();
    let kdir = layout.kset_dir(SLUG);
    let fresh = block_kset_bound(&kdir, &layout.fasta_dir(), &ids, &[]).unwrap();
    count(&layout, SLUG, 7, None).unwrap();
    let written: u64 = ids
        .iter()
        .map(|id| {
            std::fs::metadata(kdir.join(format!("{id}.kset")))
                .unwrap()
                .len()
        })
        .sum();
    let counted = load_counted(&layout, SLUG).unwrap();
    for id in &ids {
        std::fs::remove_file(kdir.join(format!("{id}.kset"))).unwrap();
    }
    let known = block_kset_bound(&kdir, &layout.fasta_dir(), &ids, &counted.genomes).unwrap();
    assert!(
        written <= known && known <= fresh,
        "written {written}, bound from counts {known}, from FASTA {fresh}"
    );
}
