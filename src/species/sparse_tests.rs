use super::*;
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
    let err = sparse(&layout, 7).unwrap_err();
    assert!(
        format!("{err:#}").contains("refusing to remove non-regular-file"),
        "{err:#}"
    );
    assert!(ksets_remain(&layout, 0..64));
    std::fs::remove_dir(blocker).unwrap();
    sparse(&layout, 7).unwrap();
    assert!(!ksets_remain(&layout, 0..64));
}

#[test]
fn sparse_keeps_partials_when_matrix_save_fails() {
    let dir = tempfile::tempdir().unwrap();
    let layout = fixture(dir.path(), 64, false);
    std::fs::create_dir_all(layout.matrix_dir(SLUG).join("kmers.bin")).unwrap();
    let err = sparse(&layout, 7).unwrap_err();
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
