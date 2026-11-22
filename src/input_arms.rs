//! Experiment, off by default (plan block C3, docs/input-arms.md): other
//! inputs beside 31-mer presence, as arms of `grm learn`.
//!
//! `grm learn --exp-inputs ARM --exp-program-tables DIR` turns it on for that
//! process; without the flag `active` is None, nothing here is called, and
//! `grm learn` writes what it wrote before, byte for byte. With it, the one
//! place both `grm learn` paths project a cohort (`analysis::project_cohort`)
//! hands the projected K31 dataset to `apply`, which returns the arm's
//! matrix, and `evaluate` then fits and scores it with the release's own
//! code: each fit's candidates at the dataset's support floor on its own
//! training genomes (`learn::training_candidates`) and the chi-square screen
//! on its own fitting genomes. An arm writes no report; its results are the
//! C9 writer's per-fold counts (`--exp-pick DIR --exp-pick-inner 0`) and one
//! `<dataset>.<tag>.inputs.json` of column counts in the same folder.
//!
//! The arms built here read the AMRFinderPlus program's per-genome tables
//! (the program itself, not GRM's own reader of its catalog):
//! - `program`: one binary column per element symbol (or subclass, or
//!   class) the program's tables name, over the cohort, and nothing else;
//! - `k31-program`: those columns first, then the dataset's K31 rows, in one
//!   matrix, so one candidate pass and one screen run over both blocks, and a
//!   program column identical on a fit's training genomes to a K31 pattern
//!   is kept as the program's.
//!
//! What this module builds reads no label: `program_block` sees each
//! genome's calls and the support floor only, and `arm_matrix` takes the
//! dataset's matrix, never its labels. The one cohort-wide filter is the
//! label-free range `Dataset::project` applies to K31, a column's cohort
//! support in [floor, n - floor]; every training subset implies it, so it
//! adds and drops no training candidate. A cohort genome with no table stops
//! the dataset; no row is ever filled. A result holds counts only: no
//! column name, genome id or row per genome.

use crate::bits::{BitMatrix, Bits};
use crate::matrix::Dataset;
use anyhow::{Context, Result, bail, ensure};
use serde::Serialize;
use std::collections::{BTreeSet, HashSet};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

/// The result format. A change to what a result holds changes it.
pub const FORMAT: &str = "grm-c3-inputs-1";

/// An arm of block C3 that this build can run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Arm {
    /// The program's columns alone.
    Program,
    /// The program's columns, then the dataset's K31 rows.
    K31Program,
}

impl Arm {
    pub fn name(self) -> &'static str {
        match self {
            Arm::Program => "program",
            Arm::K31Program => "k31-program",
        }
    }
}

/// What one program column is (the design's question 7).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Column {
    /// The table's `Element symbol`, as the program names it.
    Symbol,
    /// The table's `Subclass`.
    Subclass,
    /// The table's `Class`.
    Class,
}

impl Column {
    pub fn name(self) -> &'static str {
        match self {
            Column::Symbol => "symbol",
            Column::Subclass => "subclass",
            Column::Class => "class",
        }
    }

    fn header(self) -> &'static str {
        match self {
            Column::Symbol => "Element symbol",
            Column::Subclass => "Subclass",
            Column::Class => "Class",
        }
    }
}

/// Which table rows give a column (the design's question 8).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Rows {
    /// Rows of Type AMR, both scopes, every method.
    Amr,
    /// Every row.
    All,
}

impl Rows {
    pub fn name(self) -> &'static str {
        match self {
            Rows::Amr => "amr",
            Rows::All => "all",
        }
    }
}

/// The experiment's settings, as `enable` set them.
#[derive(Clone, Debug)]
pub struct Settings {
    pub arm: Arm,
    /// The species' program tables, `<genome id>.tsv`.
    pub tables: PathBuf,
    pub column: Column,
    pub rows: Rows,
    /// Where `<dataset>.<tag>.inputs.json` goes: the `--exp-pick` folder.
    pub out: PathBuf,
}

static SINK: OnceLock<Settings> = OnceLock::new();

/// Turns the experiment on for this process. A second call fails.
pub fn enable(settings: Settings) -> Result<()> {
    ensure!(
        settings.tables.is_dir(),
        "--exp-program-tables {} is not a folder",
        settings.tables.display()
    );
    std::fs::create_dir_all(&settings.out)
        .with_context(|| format!("create {}", settings.out.display()))?;
    SINK.set(settings)
        .map_err(|_| anyhow::anyhow!("the C3 input arm is already set"))
}

/// The arm this process runs, if any.
pub fn active() -> Option<&'static Settings> {
    SINK.get()
}

/// One genome's program calls: each name `column` gives on a row that
/// `rows` admits. A table with a header and no row is a genome with no call.
pub fn read_table(text: &str, column: Column, rows: Rows) -> Result<BTreeSet<String>> {
    let mut lines = text.lines();
    let header: Vec<&str> = lines
        .next()
        .context("an empty table")?
        .split('\t')
        .collect();
    let at = |name: &str| {
        header
            .iter()
            .position(|h| *h == name)
            .with_context(|| format!("the table has no column {name:?}"))
    };
    let col = at(column.header())?;
    let kind = at("Type")?;
    let mut out = BTreeSet::new();
    for (i, line) in lines.enumerate() {
        if line.is_empty() {
            continue;
        }
        let cells: Vec<&str> = line.split('\t').collect();
        ensure!(
            cells.len() == header.len(),
            "row {} has {} cells for {} columns",
            i + 2,
            cells.len(),
            header.len()
        );
        if rows == Rows::Amr && cells[kind] != "AMR" {
            continue;
        }
        let name = cells[col].trim();
        ensure!(!name.is_empty(), "row {} has no {}", i + 2, column.header());
        out.insert(name.to_string());
    }
    Ok(out)
}

/// The program's columns over a cohort, from each genome's calls and the
/// support floor alone: one row per name in byte order whose cohort support
/// lies in [floor, n - floor], identical patterns kept once, the first.
pub struct ProgramBlock {
    pub x: BitMatrix,
    /// Names seen in any genome's calls.
    pub seen: usize,
    /// Names whose cohort support lies outside [floor, n - floor].
    pub out_of_range: usize,
    /// Names dropped as an earlier name's identical pattern.
    pub collapsed: usize,
}

pub fn program_block(calls: &[BTreeSet<String>], min_support: u32) -> ProgramBlock {
    let n = calls.len();
    let names: BTreeSet<&String> = calls.iter().flatten().collect();
    let max = (n as u32).saturating_sub(min_support);
    let mut x = BitMatrix::with_cols(n);
    let mut seen_rows: HashSet<Vec<u64>> = HashSet::new();
    let (mut out_of_range, mut collapsed) = (0, 0);
    for name in &names {
        let row = Bits::from_indices(
            n,
            calls
                .iter()
                .enumerate()
                .filter(|(_, c)| c.contains(*name))
                .map(|(i, _)| i),
        );
        let s = row.count();
        if s < min_support || s > max {
            out_of_range += 1;
            continue;
        }
        if !seen_rows.insert(row.words().to_vec()) {
            collapsed += 1;
            continue;
        }
        x.push_row(row.words());
    }
    ProgramBlock {
        x,
        seen: names.len(),
        out_of_range,
        collapsed,
    }
}

/// The arm's matrix over the dataset's samples: the program's rows, then,
/// for `k31-program`, every K31 row of `ds` once. It reads `ds.x` only.
pub fn arm_matrix(arm: Arm, ds: &Dataset, program: &ProgramBlock) -> BitMatrix {
    assert_eq!(program.x.cols, ds.x.cols);
    let mut x = BitMatrix::with_cols(ds.x.cols);
    for r in 0..program.x.rows {
        x.push_row(program.x.row(r));
    }
    if arm == Arm::K31Program {
        for r in 0..ds.x.rows {
            x.push_row(ds.x.row(r));
        }
    }
    x
}

/// What one arm's matrix of one dataset held: counts only.
#[derive(Debug, Serialize)]
pub struct InputsResult {
    pub format: &'static str,
    pub evaluation_revision: u32,
    pub dataset: String,
    pub tag: String,
    pub arm: &'static str,
    pub column: &'static str,
    pub rows: &'static str,
    pub samples: usize,
    pub min_support: u32,
    pub tables_read: usize,
    pub program_names_seen: usize,
    pub program_names_out_of_range: usize,
    pub program_names_collapsed: usize,
    pub program_rows: usize,
    pub k31_rows: usize,
    pub arm_rows: usize,
}

/// Reads each cohort genome's table, `<dir>/<genome id>.tsv`. A genome with
/// no table stops the dataset: its row would have to be filled.
pub fn cohort_calls(
    dir: &Path,
    genome_ids: &[String],
    column: Column,
    rows: Rows,
) -> Result<Vec<BTreeSet<String>>> {
    let missing = genome_ids
        .iter()
        .filter(|id| !dir.join(format!("{id}.tsv")).is_file())
        .count();
    if missing > 0 {
        bail!(
            "{missing} of {} cohort genomes have no program table in {}; the arm fills no row",
            genome_ids.len(),
            dir.display()
        );
    }
    genome_ids
        .iter()
        .map(|id| {
            let path = dir.join(format!("{id}.tsv"));
            let text = std::fs::read_to_string(&path)
                .with_context(|| format!("read {}", path.display()))?;
            read_table(&text, column, rows).with_context(|| format!("table {}", path.display()))
        })
        .collect()
}

/// The arm's dataset in place of the K31 projection `ds`, and its counts.
/// The returned dataset keeps `ds`'s labels and columns; its rows carry no
/// k-mer patterns, so nothing may annotate them, and an arm writes no
/// report.
pub fn arm_dataset(
    settings: &Settings,
    ds: Dataset,
    dataset: &str,
    tag: &str,
    genome_ids: &[String],
    min_support: u32,
) -> Result<(Dataset, InputsResult)> {
    ensure!(
        genome_ids.len() == ds.x.cols,
        "{} genome ids for {} samples",
        genome_ids.len(),
        ds.x.cols
    );
    let calls = cohort_calls(&settings.tables, genome_ids, settings.column, settings.rows)?;
    let block = program_block(&calls, min_support);
    let x = arm_matrix(settings.arm, &ds, &block);
    let result = InputsResult {
        format: FORMAT,
        evaluation_revision: crate::analysis::EVALUATION_REVISION,
        dataset: dataset.to_string(),
        tag: tag.to_string(),
        arm: settings.arm.name(),
        column: settings.column.name(),
        rows: settings.rows.name(),
        samples: ds.x.cols,
        min_support,
        tables_read: calls.len(),
        program_names_seen: block.seen,
        program_names_out_of_range: block.out_of_range,
        program_names_collapsed: block.collapsed,
        program_rows: block.x.rows,
        k31_rows: ds.x.rows,
        arm_rows: x.rows,
    };
    let rows = x.rows;
    Ok((
        Dataset {
            x,
            y: ds.y,
            columns: ds.columns,
            feature_offsets: vec![0; rows + 1],
            feature_patterns: Vec::new(),
        },
        result,
    ))
}

/// `analysis::project_cohort`'s hook: the arm's dataset when an arm is
/// active, with its counts written to `<out>/<dataset>.<tag>.inputs.json`;
/// `ds` itself otherwise.
pub fn apply(
    ds: Dataset,
    dataset: &str,
    tag: &str,
    genome_ids: &[String],
    min_support: u32,
) -> Result<Dataset> {
    let Some(settings) = active() else {
        return Ok(ds);
    };
    let (ds, result) = arm_dataset(settings, ds, dataset, tag, genome_ids, min_support)?;
    let path = settings.out.join(format!("{dataset}.{tag}.inputs.json"));
    let bytes = serde_json::to_vec_pretty(&result)?;
    let part = path.with_extension("json.part");
    std::fs::write(&part, &bytes).with_context(|| format!("write {}", part.display()))?;
    std::fs::rename(&part, &path).with_context(|| format!("rename to {}", path.display()))?;
    eprintln!(
        "    {dataset}: input arm {}: {} program rows and {} K31 rows of {}, {} rows",
        settings.arm.name(),
        result.program_rows,
        result.k31_rows,
        if settings.arm == Arm::K31Program {
            "both"
        } else {
            "which none enter"
        },
        result.arm_rows
    );
    Ok(ds)
}

#[cfg(test)]
mod tests {
    //! SYNTHETIC tests: made-up tables and cohorts, no real genome.
    use super::*;
    use crate::learn::training_candidates;

    /// The 22 columns of a table the program 4.2.7 writes with `--plus`,
    /// as the AMRFinderPlus kit's fixture tables have them.
    const HEADER: &str = "Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\tHMM accession\tHMM description";

    fn row(symbol: &str, kind: &str, class: &str, subclass: &str) -> String {
        format!(
            "NA\tcontig1\t1\t90\t+\t{symbol}\tsynthetic element\tcore\t{kind}\t{kind}\t{class}\t{subclass}\tEXACTX\t30\t30\t100.00\t100.00\t30\tNA\tsynthetic\tNA\tNA"
        )
    }

    fn table(rows: &[String]) -> String {
        let mut t = String::from(HEADER);
        for r in rows {
            t.push('\n');
            t.push_str(r);
        }
        t.push('\n');
        t
    }

    /// A deterministic stream of bits.
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn below(&mut self, n: u64) -> bool {
            self.next() % 100 < n
        }
    }

    /// SYNTHETIC: n genomes, their calls over `names` program names, and a
    /// K31 projection of `k` rows; resistance is the first name's presence
    /// with 10% flipped.
    fn cohort(n: usize, names: usize, k: usize, seed: u64) -> (Vec<BTreeSet<String>>, Dataset) {
        let mut g = Lcg(seed);
        let calls: Vec<BTreeSet<String>> = (0..n)
            .map(|_| {
                (0..names)
                    .filter(|j| g.below(10 + 7 * (*j as u64 % 9)))
                    .map(|j| format!("synthetic_{j:02}"))
                    .collect()
            })
            .collect();
        let mut x = BitMatrix::with_cols(n);
        for _ in 0..k {
            let row = Bits::from_indices(n, (0..n).filter(|_| g.below(35)));
            x.push_row(row.words());
        }
        let y = Bits::from_indices(
            n,
            (0..n).filter(|&i| calls[i].contains("synthetic_00") ^ (i % 10 == 3)),
        );
        let ds = Dataset {
            x,
            y,
            columns: (0..n).collect(),
            feature_offsets: vec![0; k + 1],
            feature_patterns: Vec::new(),
        };
        (calls, ds)
    }

    fn rows_of(x: &BitMatrix) -> Vec<Vec<u64>> {
        (0..x.rows).map(|r| x.row(r).to_vec()).collect()
    }

    #[test]
    fn a_table_gives_the_named_column_of_the_admitted_rows() {
        let t = table(&[
            row("blaSYN-1", "AMR", "BETA-LACTAM", "CARBAPENEM"),
            row("gyrA_S00L", "AMR", "QUINOLONE", "QUINOLONE"),
            row("stressSYN", "STRESS", "METAL", "COPPER"),
            row("blaSYN-1", "AMR", "BETA-LACTAM", "CARBAPENEM"),
        ]);
        let amr = read_table(&t, Column::Symbol, Rows::Amr).unwrap();
        assert_eq!(
            amr.into_iter().collect::<Vec<_>>(),
            vec!["blaSYN-1", "gyrA_S00L"]
        );
        let all = read_table(&t, Column::Symbol, Rows::All).unwrap();
        assert_eq!(all.len(), 3);
        let sub = read_table(&t, Column::Subclass, Rows::Amr).unwrap();
        assert_eq!(
            sub.into_iter().collect::<Vec<_>>(),
            vec!["CARBAPENEM", "QUINOLONE"]
        );
        let class = read_table(&t, Column::Class, Rows::All).unwrap();
        assert_eq!(class.len(), 3);
        // A header with no row is a genome with no call.
        assert!(
            read_table(&table(&[]), Column::Symbol, Rows::Amr)
                .unwrap()
                .is_empty()
        );
        // A table without the named column, or a short row, is refused.
        assert!(read_table("Type\tClass\nAMR\tX\n", Column::Symbol, Rows::Amr).is_err());
        let short = format!("{HEADER}\nNA\tcontig1\n");
        assert!(read_table(&short, Column::Symbol, Rows::Amr).is_err());
    }

    #[test]
    fn the_program_block_keeps_the_cohort_range_and_each_pattern_once() {
        let s = |v: &[&str]| v.iter().map(|x| x.to_string()).collect::<BTreeSet<_>>();
        // SYNTHETIC: 8 genomes, floor 2.
        let calls = vec![
            s(&["a", "b", "c", "e"]),
            s(&["a", "b", "e"]),
            s(&["a", "e"]),
            s(&["a", "d", "e"]),
            s(&["a", "e"]),
            s(&["a", "e"]),
            s(&["a", "e"]),
            s(&["e"]),
        ];
        let block = program_block(&calls, 2);
        // a: support 7 > 8 - 2, out; b: 2, kept; c: 1, out; d: 1, out; e: 8, out.
        assert_eq!(block.seen, 5);
        assert_eq!(block.out_of_range, 4);
        assert_eq!(block.collapsed, 0);
        assert_eq!(block.x.rows, 1);
        assert_eq!(block.x.row(0), Bits::from_indices(8, [0, 1]).words());
        // Two names with one pattern: the first in byte order is kept.
        let twin = vec![s(&["x", "w"]), s(&["x", "w"]), s(&[]), s(&[]), s(&[])];
        let block = program_block(&twin, 2);
        assert_eq!((block.x.rows, block.collapsed), (1, 1));
    }

    #[test]
    fn the_joined_arm_holds_each_block_once_the_program_first() {
        let (calls, ds) = cohort(60, 12, 40, 7);
        let block = program_block(&calls, 3);
        assert!(block.x.rows >= 3, "the fixture gives program rows");
        let joined = arm_matrix(Arm::K31Program, &ds, &block);
        assert_eq!(joined.rows, block.x.rows + ds.x.rows);
        let rows = rows_of(&joined);
        assert_eq!(&rows[..block.x.rows], &rows_of(&block.x)[..]);
        assert_eq!(&rows[block.x.rows..], &rows_of(&ds.x)[..]);
        let alone = arm_matrix(Arm::Program, &ds, &block);
        assert_eq!(rows_of(&alone), rows_of(&block.x));
    }

    #[test]
    fn an_arm_matrix_does_not_change_with_the_labels() {
        let (calls, ds) = cohort(60, 12, 40, 11);
        let block = program_block(&calls, 3);
        let flipped = Dataset {
            x: ds.x.clone(),
            y: ds.y.not(),
            columns: ds.columns.clone(),
            feature_offsets: ds.feature_offsets.clone(),
            feature_patterns: Vec::new(),
        };
        for arm in [Arm::Program, Arm::K31Program] {
            assert_eq!(
                rows_of(&arm_matrix(arm, &ds, &block)),
                rows_of(&arm_matrix(arm, &flipped, &block))
            );
        }
    }

    /// The candidates a fit on `train` gets, as their training patterns.
    fn candidates(x: &BitMatrix, train: &Bits, floor: u32) -> Vec<Vec<u64>> {
        let (sub, _) = training_candidates(x, train, floor);
        (0..sub.rows)
            .map(|r| {
                sub.row(r)
                    .iter()
                    .zip(train.words())
                    .map(|(a, b)| a & b)
                    .collect()
            })
            .collect()
    }

    #[test]
    fn no_outer_test_genome_changes_a_fits_candidates() {
        // SYNTHETIC: 60 genomes in 5 outer folds of 12. For each fold, every
        // test genome's program calls are changed in turn (all removed, all
        // given); the fit's candidates on the joined matrix must not move.
        let floor = 3;
        for seed in [3, 5, 9] {
            let (calls, ds) = cohort(60, 16, 30, seed);
            for f in 0..5 {
                let test: Vec<usize> = (f * 12..f * 12 + 12).collect();
                let train = Bits::from_indices(60, (0..60).filter(|i| !test.contains(i)));
                let base = arm_matrix(Arm::K31Program, &ds, &program_block(&calls, floor));
                let want = candidates(&base, &train, floor);
                let every: BTreeSet<String> = calls.iter().flatten().cloned().collect();
                for fill in [BTreeSet::new(), every] {
                    let mut moved = calls.clone();
                    for &i in &test {
                        moved[i] = fill.clone();
                    }
                    let x = arm_matrix(Arm::K31Program, &ds, &program_block(&moved, floor));
                    assert_eq!(candidates(&x, &train, floor), want, "seed {seed}, fold {f}");
                }
            }
        }
    }

    #[test]
    fn a_cohort_genome_without_a_table_stops_the_dataset() {
        let dir = std::env::temp_dir().join(format!("grm-c3-tables-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(
            dir.join("g1.tsv"),
            table(&[row("blaSYN-1", "AMR", "B", "C")]),
        )
        .unwrap();
        std::fs::write(dir.join("g2.tsv"), table(&[])).unwrap();
        let ids = |v: &[&str]| v.iter().map(|x| x.to_string()).collect::<Vec<_>>();
        let got = cohort_calls(&dir, &ids(&["g1", "g2"]), Column::Symbol, Rows::Amr).unwrap();
        assert_eq!(got[0].len(), 1);
        assert!(got[1].is_empty());
        let err = cohort_calls(&dir, &ids(&["g1", "g3"]), Column::Symbol, Rows::Amr)
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("1 of 2 cohort genomes have no program table"),
            "{err}"
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn the_arm_dataset_keeps_labels_and_columns_and_counts_its_rows() {
        let (calls, ds) = cohort(40, 10, 20, 13);
        let dir = std::env::temp_dir().join(format!("grm-c3-arm-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let ids: Vec<String> = (0..40).map(|i| format!("syn{i}")).collect();
        for (id, c) in ids.iter().zip(&calls) {
            let rows: Vec<String> = c.iter().map(|s| row(s, "AMR", "B", "C")).collect();
            std::fs::write(dir.join(format!("{id}.tsv")), table(&rows)).unwrap();
        }
        for arm in [Arm::Program, Arm::K31Program] {
            let settings = Settings {
                arm,
                tables: dir.clone(),
                column: Column::Symbol,
                rows: Rows::Amr,
                out: dir.clone(),
            };
            let y = ds.y.clone();
            let k = ds.x.rows;
            let copy = Dataset {
                x: ds.x.clone(),
                y: ds.y.clone(),
                columns: ds.columns.clone(),
                feature_offsets: ds.feature_offsets.clone(),
                feature_patterns: Vec::new(),
            };
            let (got, result) =
                arm_dataset(&settings, copy, "synthetic__drug", "drawn", &ids, 3).unwrap();
            assert_eq!(got.y, y);
            assert_eq!(got.columns, ds.columns);
            assert_eq!(got.feature_offsets.len(), got.x.rows + 1);
            let p = program_block(&calls, 3).x.rows;
            assert_eq!(result.program_rows, p);
            assert_eq!(result.k31_rows, k);
            let want = if arm == Arm::K31Program { p + k } else { p };
            assert_eq!((result.arm_rows, got.x.rows), (want, want));
            assert_eq!(result.tables_read, 40);
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
