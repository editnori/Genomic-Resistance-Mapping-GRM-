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
//!   matrix, so one candidate pass and one screen run over both blocks. A
//!   K31 row whose cohort pattern equals a program column's is left out, so
//!   the matrix holds each cohort pattern once, as the final models need,
//!   and a program column identical on a fit's training genomes to a K31
//!   pattern is kept as the program's.
//!
//! The three matrix arms read a pattern matrix that `grm exp-input-build`
//! made once per species with no label (`input_build`), and project it
//! onto the cohort at the dataset's floor, with every label set to S:
//! - `k31-counts`: the 31-mers a genome holds 2 or more times, then K31;
//! - `k31-aa`: the amino-acid 10-mers of Prodigal's proteins, then K31;
//! - `unitigs`: the unitigs of the species' graph, in K31's place.
//!
//! A joined arm (`k31-program`, `k31-counts`, `k31-aa`) also writes, per
//! split and seed, how many of each block's candidates pass the 4,000
//! screen of each screened learner (`<dataset>.<scheme>-s<seed>.screen.json`).
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
use crate::input_build::{self, Build, Kind};
use crate::matrix::{Dataset, FeatureMatrix};
use anyhow::{Context, Result, bail, ensure};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::{Mutex, OnceLock};

/// The result format. A change to what a result holds changes it.
pub const FORMAT: &str = "grm-c3-inputs-1";

/// An arm of block C3 that this build can run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Arm {
    /// The program's columns alone.
    Program,
    /// The program's columns, then the dataset's K31 rows.
    K31Program,
    /// The 31-mers seen 2 or more times in a genome, then K31.
    K31Counts,
    /// The amino-acid 10-mers of called proteins, then K31.
    K31Aa,
    /// The unitigs of the species' graph, in K31's place.
    Unitigs,
}

impl Arm {
    pub fn name(self) -> &'static str {
        match self {
            Arm::Program => "program",
            Arm::K31Program => "k31-program",
            Arm::K31Counts => "k31-counts",
            Arm::K31Aa => "k31-aa",
            Arm::Unitigs => "unitigs",
        }
    }

    /// Whether the arm reads the program's tables.
    pub fn is_program(self) -> bool {
        matches!(self, Arm::Program | Arm::K31Program)
    }

    /// Whether the arm joins a block before K31's rows.
    pub fn joined(self) -> bool {
        matches!(self, Arm::K31Program | Arm::K31Counts | Arm::K31Aa)
    }

    /// The label-free build a matrix arm reads.
    pub fn build_kind(self) -> Option<Kind> {
        match self {
            Arm::K31Counts => Some(Kind::Counts),
            Arm::K31Aa => Some(Kind::Aa),
            Arm::Unitigs => Some(Kind::Unitigs),
            Arm::Program | Arm::K31Program => None,
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
    /// The species' program tables, `<genome id>.tsv`, for a program arm.
    pub tables: Option<PathBuf>,
    /// `grm exp-input-build`'s folder, for a matrix arm.
    pub matrix: Option<PathBuf>,
    pub column: Column,
    pub rows: Rows,
    /// Where `<dataset>.<tag>.inputs.json` goes: the `--exp-pick` folder.
    pub out: PathBuf,
}

static SINK: OnceLock<Settings> = OnceLock::new();

/// A matrix arm's build: its matrix, its counts, and for the unitig arm each
/// genome's unitigs held whole and in part.
struct Added {
    matrix: FeatureMatrix,
    build: Build,
    pairs: Option<HashMap<String, (u64, u64)>>,
}

static ADDED: OnceLock<Added> = OnceLock::new();

/// Per dataset, the rows of a joined arm's first block, for the screen count.
static JOINED: Mutex<BTreeMap<String, usize>> = Mutex::new(BTreeMap::new());

/// Turns the experiment on for this process. A second call fails.
pub fn enable(settings: Settings) -> Result<()> {
    if settings.arm.is_program() {
        let tables = settings
            .tables
            .as_deref()
            .context("a program arm needs --exp-program-tables")?;
        ensure!(
            tables.is_dir(),
            "--exp-program-tables {} is not a folder",
            tables.display()
        );
    } else {
        let dir = settings
            .matrix
            .as_deref()
            .context("a matrix arm needs --exp-input-matrix")?;
        let (matrix, build) = input_build::load(dir)?;
        let want = settings.arm.build_kind().expect("a matrix arm").name();
        ensure!(
            build.kind == want,
            "--exp-inputs {} reads a {want} build, and {} holds a {} build",
            settings.arm.name(),
            dir.display(),
            build.kind
        );
        let pairs = if settings.arm == Arm::Unitigs {
            Some(input_build::unitig_pairs(dir)?)
        } else {
            None
        };
        ADDED
            .set(Added {
                matrix,
                build,
                pairs,
            })
            .map_err(|_| anyhow::anyhow!("the C3 input matrix is already loaded"))?;
    }
    std::fs::create_dir_all(&settings.out)
        .with_context(|| format!("create {}", settings.out.display()))?;
    SINK.set(settings)
        .map_err(|_| anyhow::anyhow!("the C3 input arm is already set"))
}

/// A joined arm's first-block rows in `dataset`'s matrix, once `apply` made
/// it; None for any other arm or none at all.
pub fn joined_rows(dataset: &str) -> Option<usize> {
    JOINED.lock().ok()?.get(dataset).copied()
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

/// The arm's matrix over the dataset's samples, and how many K31 rows it
/// left out as a program row's pattern: the program's rows, then, for
/// `k31-program`, each K31 row of `ds` whose cohort pattern no program row
/// has. It reads `ds.x` only.
pub fn arm_matrix(arm: Arm, ds: &Dataset, program: &ProgramBlock) -> (BitMatrix, usize) {
    assert_eq!(program.x.cols, ds.x.cols);
    let mut x = BitMatrix::with_cols(ds.x.cols);
    let mut held: HashSet<&[u64]> = HashSet::new();
    for r in 0..program.x.rows {
        x.push_row(program.x.row(r));
        held.insert(program.x.row(r));
    }
    let mut as_program = 0;
    if arm == Arm::K31Program {
        for r in 0..ds.x.rows {
            if held.contains(ds.x.row(r)) {
                as_program += 1;
                continue;
            }
            x.push_row(ds.x.row(r));
        }
    }
    (x, as_program)
}

/// What one arm's matrix of one dataset held: counts only.
#[derive(Debug, Serialize)]
pub struct InputsResult {
    pub format: &'static str,
    pub evaluation_revision: u32,
    pub dataset: String,
    pub tag: String,
    pub arm: &'static str,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub column: Option<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rows: Option<&'static str>,
    pub samples: usize,
    pub min_support: u32,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tables_read: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub program_names_seen: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub program_names_out_of_range: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub program_names_collapsed: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub program_rows: Option<usize>,
    pub k31_rows: usize,
    /// K31 rows left out of a joined arm as a program row's cohort pattern.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub k31_rows_as_program: Option<usize>,
    pub arm_rows: usize,
    /// A matrix arm: its build's kind and species patterns, the added rows
    /// at the dataset's floor, and the K31 rows left out of a joined arm as
    /// an added row's cohort pattern.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub added_build: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub added_species_patterns: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub added_rows: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub k31_rows_as_added: Option<usize>,
    /// The unitig arm: the cohort's genome and unitig pairs held whole and
    /// held in part (some but not every 31-mer of the unitig).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unitig_pairs_held: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unitig_pairs_in_part: Option<u64>,
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
    let tables = settings
        .tables
        .as_deref()
        .context("a program arm needs --exp-program-tables")?;
    let calls = cohort_calls(tables, genome_ids, settings.column, settings.rows)?;
    let block = program_block(&calls, min_support);
    let (x, as_program) = arm_matrix(settings.arm, &ds, &block);
    let result = InputsResult {
        format: FORMAT,
        evaluation_revision: crate::analysis::EVALUATION_REVISION,
        dataset: dataset.to_string(),
        tag: tag.to_string(),
        arm: settings.arm.name(),
        column: Some(settings.column.name()),
        rows: Some(settings.rows.name()),
        samples: ds.x.cols,
        min_support,
        tables_read: Some(calls.len()),
        program_names_seen: Some(block.seen),
        program_names_out_of_range: Some(block.out_of_range),
        program_names_collapsed: Some(block.collapsed),
        program_rows: Some(block.x.rows),
        k31_rows: ds.x.rows,
        k31_rows_as_program: Some(as_program),
        arm_rows: x.rows,
        added_build: None,
        added_species_patterns: None,
        added_rows: None,
        k31_rows_as_added: None,
        unitig_pairs_held: None,
        unitig_pairs_in_part: None,
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

/// A matrix arm's block over the cohort: the build's matrix projected onto
/// the cohort's species-matrix columns at the dataset's floor, as
/// `Dataset::project` projects K31, with every label set to S so that no
/// label can reach it. Identical cohort patterns collapse to the first.
pub fn added_block(matrix: &FeatureMatrix, columns: &[usize], min_support: u32) -> BitMatrix {
    let susceptible = vec![false; columns.len()];
    Dataset::project(matrix, columns.to_vec(), &susceptible, min_support).x
}

/// A matrix arm's matrix over the dataset's samples, and how many K31 rows
/// it left out as an added row's pattern: the added rows, then, for a joined
/// arm, each K31 row of `ds` whose cohort pattern no added row has. It reads
/// `ds.x` only.
pub fn matrix_arm(joined: bool, ds: &Dataset, added: &BitMatrix) -> (BitMatrix, usize) {
    assert_eq!(added.cols, ds.x.cols);
    let mut x = BitMatrix::with_cols(ds.x.cols);
    let mut first: HashSet<&[u64]> = HashSet::new();
    for r in 0..added.rows {
        x.push_row(added.row(r));
        first.insert(added.row(r));
    }
    let mut as_added = 0;
    if joined {
        for r in 0..ds.x.rows {
            if first.contains(ds.x.row(r)) {
                as_added += 1;
                continue;
            }
            x.push_row(ds.x.row(r));
        }
    }
    (x, as_added)
}

/// `arm_dataset` for a matrix arm: the build's block over the cohort, joined
/// before K31's rows or in their place.
#[allow(clippy::too_many_arguments)]
fn matrix_dataset(
    settings: &Settings,
    added: &Added,
    ds: Dataset,
    dataset: &str,
    tag: &str,
    genome_ids: &[String],
    species_genomes: &[String],
    min_support: u32,
) -> Result<(Dataset, InputsResult)> {
    ensure!(
        genome_ids.len() == ds.x.cols && ds.columns.len() == ds.x.cols,
        "{} genome ids for {} samples",
        genome_ids.len(),
        ds.x.cols
    );
    ensure!(
        added.matrix.summary.genomes == species_genomes,
        "the {} build's genomes are not the species matrix's, in order",
        added.build.kind
    );
    let block = added_block(&added.matrix, &ds.columns, min_support);
    let (x, as_added) = matrix_arm(settings.arm.joined(), &ds, &block);
    let (held, in_part) = match &added.pairs {
        Some(pairs) => {
            let mut sum = (0u64, 0u64);
            for id in genome_ids {
                let (h, p) = pairs.get(id).with_context(|| {
                    format!("the unitig build has no row for a cohort genome of {dataset}")
                })?;
                sum = (sum.0 + h, sum.1 + p);
            }
            (Some(sum.0), Some(sum.1))
        }
        None => (None, None),
    };
    let result = InputsResult {
        format: FORMAT,
        evaluation_revision: crate::analysis::EVALUATION_REVISION,
        dataset: dataset.to_string(),
        tag: tag.to_string(),
        arm: settings.arm.name(),
        column: None,
        rows: None,
        samples: ds.x.cols,
        min_support,
        tables_read: None,
        program_names_seen: None,
        program_names_out_of_range: None,
        program_names_collapsed: None,
        program_rows: None,
        k31_rows: ds.x.rows,
        k31_rows_as_program: None,
        arm_rows: x.rows,
        added_build: Some(added.build.kind.clone()),
        added_species_patterns: Some(added.matrix.summary.n_patterns),
        added_rows: Some(block.rows),
        k31_rows_as_added: settings.arm.joined().then_some(as_added),
        unitig_pairs_held: held,
        unitig_pairs_in_part: in_part,
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
/// `ds` itself otherwise. `species_genomes` are the species matrix's, in
/// column order.
pub fn apply(
    ds: Dataset,
    dataset: &str,
    tag: &str,
    genome_ids: &[String],
    species_genomes: &[String],
    min_support: u32,
) -> Result<Dataset> {
    let Some(settings) = active() else {
        return Ok(ds);
    };
    let (ds, result) = if settings.arm.is_program() {
        arm_dataset(settings, ds, dataset, tag, genome_ids, min_support)?
    } else {
        let added = ADDED.get().context("the C3 input matrix is not loaded")?;
        matrix_dataset(
            settings,
            added,
            ds,
            dataset,
            tag,
            genome_ids,
            species_genomes,
            min_support,
        )?
    };
    if settings.arm.joined() {
        let first = result
            .program_rows
            .or(result.added_rows)
            .context("a joined arm with no first block")?;
        JOINED
            .lock()
            .map_err(|_| anyhow::anyhow!("the joined-rows lock is poisoned"))?
            .insert(dataset.to_string(), first);
    }
    let path = settings.out.join(format!("{dataset}.{tag}.inputs.json"));
    let bytes = serde_json::to_vec_pretty(&result)?;
    let part = path.with_extension("json.part");
    std::fs::write(&part, &bytes).with_context(|| format!("write {}", part.display()))?;
    std::fs::rename(&part, &path).with_context(|| format!("rename to {}", path.display()))?;
    if settings.arm.is_program() {
        eprintln!(
            "    {dataset}: input arm {}: {} program rows and {} K31 rows of {}, {} as a program row, {} rows",
            settings.arm.name(),
            result.program_rows.unwrap_or(0),
            result.k31_rows,
            if settings.arm == Arm::K31Program {
                "both"
            } else {
                "which none enter"
            },
            result.k31_rows_as_program.unwrap_or(0),
            result.arm_rows
        );
    } else {
        eprintln!(
            "    {dataset}: input arm {}: {} added rows and {} K31 rows of {}, {} as an added row, {} rows",
            settings.arm.name(),
            result.added_rows.unwrap_or(0),
            result.k31_rows,
            if settings.arm.joined() {
                "both"
            } else {
                "which none enter"
            },
            result.k31_rows_as_added.unwrap_or(0),
            result.arm_rows
        );
    }
    Ok(ds)
}

/// The screen file's format.
pub const SCREEN_FORMAT: &str = "grm-c3-screen-1";

/// One screened learner's outer fit on one fold: its fitting genomes, the
/// features its screen kept, and how many of them each block gave.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct ScreenLearner {
    pub learner: &'static str,
    pub fitting: usize,
    pub screened: usize,
    pub passed_first: usize,
    pub passed_k31: usize,
}

/// One outer fold of a joined arm: its training genomes, each block's
/// candidates (`training_candidates`), and the screened learners.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct ScreenFold {
    pub fold: usize,
    pub train: usize,
    pub candidates_first: usize,
    pub candidates_k31: usize,
    pub learners: Vec<ScreenLearner>,
}

#[derive(Serialize)]
struct ScreenResult<'a> {
    format: &'static str,
    evaluation_revision: u32,
    dataset: &'a str,
    scheme: &'a str,
    seed: u64,
    arm: &'static str,
    first_block: &'static str,
    screened_features: usize,
    folds: &'a [ScreenFold],
}

/// Writes a joined arm's screen count for one split and seed beside its
/// pick file: `<out>/<dataset>.<scheme>-s<seed>.screen.json`. Counts only.
pub fn write_screen(
    dataset: &str,
    scheme: &str,
    seed: u64,
    screened_features: usize,
    folds: &[ScreenFold],
) -> Result<()> {
    let settings = active().context("no C3 input arm is set")?;
    let doc = ScreenResult {
        format: SCREEN_FORMAT,
        evaluation_revision: crate::analysis::EVALUATION_REVISION,
        dataset,
        scheme,
        seed,
        arm: settings.arm.name(),
        first_block: match settings.arm {
            Arm::K31Program => "program",
            Arm::K31Counts => "counts",
            Arm::K31Aa => "aa",
            Arm::Program | Arm::Unitigs => bail!("{} joins no block", settings.arm.name()),
        },
        screened_features,
        folds,
    };
    let path = settings
        .out
        .join(format!("{dataset}.{scheme}-s{seed}.screen.json"));
    let part = path.with_extension("json.part");
    std::fs::write(&part, serde_json::to_vec_pretty(&doc)?)
        .with_context(|| format!("write {}", part.display()))?;
    std::fs::rename(&part, &path).with_context(|| format!("rename to {}", path.display()))?;
    Ok(())
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
        let (joined, none) = arm_matrix(Arm::K31Program, &ds, &block);
        assert_eq!(none, 0, "the fixture's K31 rows hold no program pattern");
        assert_eq!(joined.rows, block.x.rows + ds.x.rows);
        let rows = rows_of(&joined);
        assert_eq!(&rows[..block.x.rows], &rows_of(&block.x)[..]);
        assert_eq!(&rows[block.x.rows..], &rows_of(&ds.x)[..]);
        let (alone, zero) = arm_matrix(Arm::Program, &ds, &block);
        assert_eq!((rows_of(&alone), zero), (rows_of(&block.x), 0));
    }

    #[test]
    fn a_k31_row_equal_to_a_program_column_enters_once_as_the_programs() {
        // SYNTHETIC: two of the K31 rows are copies of program columns over
        // the cohort. The joined matrix holds each cohort pattern once, the
        // program's copy first, and so passes the final models' check that
        // every row is a candidate on all samples.
        let (calls, mut ds) = cohort(60, 12, 40, 7);
        let block = program_block(&calls, 3);
        assert!(block.x.rows >= 2, "the fixture gives program rows");
        let k = ds.x.rows;
        let mut x = BitMatrix::with_cols(ds.x.cols);
        for r in 0..k {
            x.push_row(ds.x.row(r));
            if r == 5 {
                x.push_row(block.x.row(0));
            }
        }
        x.push_row(block.x.row(1));
        ds.x = x;
        let (joined, as_program) = arm_matrix(Arm::K31Program, &ds, &block);
        assert_eq!(as_program, 2);
        assert_eq!(joined.rows, block.x.rows + k);
        let rows = rows_of(&joined);
        assert_eq!(&rows[..block.x.rows], &rows_of(&block.x)[..]);
        let distinct: HashSet<&Vec<u64>> = rows.iter().collect();
        assert_eq!(distinct.len(), rows.len());
        let all = Bits::ones(ds.x.cols);
        assert_eq!(training_candidates(&joined, &all, 3).1.len(), joined.rows);
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
                rows_of(&arm_matrix(arm, &ds, &block).0),
                rows_of(&arm_matrix(arm, &flipped, &block).0)
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
                let (base, _) = arm_matrix(Arm::K31Program, &ds, &program_block(&calls, floor));
                let want = candidates(&base, &train, floor);
                let every: BTreeSet<String> = calls.iter().flatten().cloned().collect();
                for fill in [BTreeSet::new(), every] {
                    let mut moved = calls.clone();
                    for &i in &test {
                        moved[i] = fill.clone();
                    }
                    let (x, _) = arm_matrix(Arm::K31Program, &ds, &program_block(&moved, floor));
                    assert_eq!(candidates(&x, &train, floor), want, "seed {seed}, fold {f}");
                }
            }
        }
    }

    #[test]
    fn the_folds_at_a_seed_are_the_same_genomes_with_and_without_an_arm() {
        // SYNTHETIC: one cohort's K31 projection and both arms' matrices,
        // with the same labels, lineages and duplicate group. `evaluate`
        // draws each split's outer folds from those and the seed alone, so
        // every arm's folds at seed s are K31's, genome by genome, which
        // the scorer can only check by class counts.
        use crate::learn::eval::{EvalConfig, evaluate};
        use crate::learn::gbt::GbtParams;
        use crate::learn::tree::ForestParams;
        let n = 120;
        let (calls, ds) = cohort(n, 12, 40, 21);
        let block = program_block(&calls, 3);
        assert!(block.x.rows >= 3, "the fixture gives program rows");
        let groups: Vec<u32> = (0..n).map(|i| (i / 6) as u32).collect();
        let duplicates = vec![vec![4, 5]];
        let mut available = 0;
        for seed in [42, 43, 44] {
            let cfg = EvalConfig {
                folds: 5,
                inner_folds: 2,
                scm_max_rules: 1,
                cart_max_depth: 1,
                forest: ForestParams {
                    n_trees: 2,
                    max_depth: 2,
                    ..Default::default()
                },
                gbt: GbtParams {
                    rounds: 2,
                    max_depth: 1,
                    ..Default::default()
                },
                screened_features: 4,
                lasso_lambdas: 2,
                seed,
                ..Default::default()
            };
            let run = |x: &BitMatrix| {
                evaluate(
                    "synthetic__drug",
                    x,
                    &ds.y,
                    &groups,
                    &duplicates,
                    3,
                    None,
                    &cfg,
                )
                .unwrap()
            };
            let k31 = run(&ds.x);
            for arm in [Arm::Program, Arm::K31Program] {
                let (x, _) = arm_matrix(arm, &ds, &block);
                assert_ne!(x.rows, ds.x.rows);
                let got = run(&x);
                assert_eq!(got.splits.len(), k31.splits.len());
                for (a, b) in k31.splits.iter().zip(&got.splits) {
                    assert_eq!(a.scheme, b.scheme);
                    assert_eq!(a.unavailable_reason, b.unavailable_reason, "{}", a.scheme);
                    assert_eq!(a.test_fold, b.test_fold, "{} at seed {seed}", a.scheme);
                    available += (a.unavailable_reason.is_none() && a.scheme == "lineage") as usize;
                }
            }
        }
        assert!(
            available > 0,
            "the fixture's lineage split runs at some seed"
        );
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

    /// SYNTHETIC: a species matrix of `rows` patterns over `n` genomes, as
    /// `grm exp-input-build` saves one, each genome in a pattern at 30%.
    fn species_matrix(n: usize, rows: usize, seed: u64) -> FeatureMatrix {
        species_matrix_at(n, &vec![30; rows], seed)
    }

    /// SYNTHETIC: as `species_matrix`, with pattern r holding each genome at
    /// `percent[r]`%, so a low one gives patterns held by a few genomes.
    fn species_matrix_at(n: usize, percent: &[u64], seed: u64) -> FeatureMatrix {
        let rows = percent.len();
        let mut g = Lcg(seed);
        let mut patterns = BitMatrix::with_cols(n);
        for &p in percent {
            patterns.push_row(Bits::from_indices(n, (0..n).filter(|_| g.below(p))).words());
        }
        FeatureMatrix {
            summary: crate::matrix::MatrixSummary {
                k: 31,
                genomes: (0..n).map(|i| format!("syn{i}")).collect(),
                min_support: 2,
                n_patterns: rows,
                kmers_distinct: rows as u64,
                kmers_kept: rows as u64,
                kmers_rare: 0,
                kmers_core: 0,
                kset_digests: Vec::new(),
            },
            patterns,
            support: vec![0; rows],
            kmer_offsets: (0..=rows as u64).collect(),
            kmers: (0..rows as u64).collect(),
        }
    }

    const MATRIX_ARMS: [Arm; 3] = [Arm::K31Counts, Arm::K31Aa, Arm::Unitigs];

    #[test]
    fn a_matrix_arm_holds_each_block_once_the_added_first() {
        // SYNTHETIC: a cohort's K31 rows and an added species matrix over 80
        // genomes, the cohort its first 60. Two K31 rows copy added rows
        // over the cohort; a joined arm keeps the added copy alone.
        let (_, mut ds) = cohort(60, 4, 40, 17);
        let m = species_matrix(80, 50, 19);
        let block = added_block(&m, &ds.columns, 3);
        assert!(block.rows >= 10, "the fixture gives added rows");
        let k = ds.x.rows;
        let mut x = BitMatrix::with_cols(ds.x.cols);
        for r in 0..k {
            x.push_row(ds.x.row(r));
            if r == 7 {
                x.push_row(block.row(2));
            }
        }
        x.push_row(block.row(5));
        ds.x = x;
        for arm in MATRIX_ARMS {
            let (got, as_added) = matrix_arm(arm.joined(), &ds, &block);
            let rows = rows_of(&got);
            assert_eq!(&rows[..block.rows], &rows_of(&block)[..], "{}", arm.name());
            if arm.joined() {
                assert_eq!(as_added, 2);
                assert_eq!(got.rows, block.rows + k, "{}", arm.name());
                let distinct: HashSet<&Vec<u64>> = rows.iter().collect();
                assert_eq!(distinct.len(), rows.len());
                let all = Bits::ones(ds.x.cols);
                assert_eq!(training_candidates(&got, &all, 3).1.len(), got.rows);
            } else {
                assert_eq!((got.rows, as_added), (block.rows, 0), "{}", arm.name());
            }
        }
    }

    #[test]
    fn a_matrix_arm_does_not_change_with_the_labels() {
        let (_, ds) = cohort(60, 4, 40, 23);
        let m = species_matrix(70, 60, 29);
        let flipped = Dataset {
            x: ds.x.clone(),
            y: ds.y.not(),
            columns: ds.columns.clone(),
            feature_offsets: ds.feature_offsets.clone(),
            feature_patterns: Vec::new(),
        };
        for arm in MATRIX_ARMS {
            let a = matrix_arm(arm.joined(), &ds, &added_block(&m, &ds.columns, 3)).0;
            let b = matrix_arm(
                arm.joined(),
                &flipped,
                &added_block(&m, &flipped.columns, 3),
            )
            .0;
            assert_eq!(rows_of(&a), rows_of(&b), "{}", arm.name());
        }
    }

    #[test]
    fn no_outer_test_genome_changes_a_matrix_arms_candidates() {
        // SYNTHETIC: 60 cohort genomes of 72 in 5 outer folds of 12. For each
        // fold, the outer test genomes' presence in the species' added
        // matrix is cleared, then filled; no fit's candidates may move. Half
        // the added patterns hold a genome at 8%, so some sit near the floor
        // on the training genomes, where a count over the whole cohort would
        // keep or drop them.
        let floor = 3;
        let percent: Vec<u64> = (0..80).map(|r| if r % 2 == 0 { 30 } else { 8 }).collect();
        for seed in [3, 5, 9] {
            let (_, ds) = cohort(60, 4, 30, seed);
            let m = species_matrix_at(72, &percent, seed + 100);
            for arm in MATRIX_ARMS {
                for f in 0..5 {
                    let test: Vec<usize> = (f * 12..f * 12 + 12).collect();
                    let train = Bits::from_indices(60, (0..60).filter(|i| !test.contains(i)));
                    let base =
                        matrix_arm(arm.joined(), &ds, &added_block(&m, &ds.columns, floor)).0;
                    let want = candidates(&base, &train, floor);
                    for fill in [false, true] {
                        let mut moved = species_matrix_at(72, &percent, seed + 100);
                        for r in 0..moved.patterns.rows {
                            let mut row = Bits::from_indices(
                                72,
                                (0..72).filter(|&c| crate::bits::get(m.patterns.row(r), c)),
                            );
                            for &i in &test {
                                let c = ds.columns[i];
                                row = if fill {
                                    row.or(&Bits::from_indices(72, [c]))
                                } else {
                                    row.and_not(&Bits::from_indices(72, [c]))
                                };
                            }
                            let w = moved.patterns.words;
                            moved.patterns.data[r * w..(r + 1) * w].copy_from_slice(row.words());
                        }
                        let x =
                            matrix_arm(arm.joined(), &ds, &added_block(&moved, &ds.columns, floor))
                                .0;
                        assert_eq!(
                            candidates(&x, &train, floor),
                            want,
                            "{}, seed {seed}, fold {f}",
                            arm.name()
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn the_folds_at_a_seed_are_the_same_genomes_with_a_matrix_arm() {
        // SYNTHETIC: as the program arms' test, for the three matrix arms.
        use crate::learn::eval::{EvalConfig, evaluate};
        use crate::learn::gbt::GbtParams;
        use crate::learn::tree::ForestParams;
        let n = 120;
        let (_, ds) = cohort(n, 4, 40, 31);
        let m = species_matrix(n, 60, 37);
        let block = added_block(&m, &ds.columns, 3);
        let groups: Vec<u32> = (0..n).map(|i| (i / 6) as u32).collect();
        let duplicates = vec![vec![4, 5]];
        let mut available = 0;
        for seed in [42, 43, 44] {
            let cfg = EvalConfig {
                folds: 5,
                inner_folds: 2,
                scm_max_rules: 1,
                cart_max_depth: 1,
                forest: ForestParams {
                    n_trees: 2,
                    max_depth: 2,
                    ..Default::default()
                },
                gbt: GbtParams {
                    rounds: 2,
                    max_depth: 1,
                    ..Default::default()
                },
                screened_features: 4,
                lasso_lambdas: 2,
                seed,
                ..Default::default()
            };
            let run = |x: &BitMatrix| {
                evaluate(
                    "synthetic__drug",
                    x,
                    &ds.y,
                    &groups,
                    &duplicates,
                    3,
                    None,
                    &cfg,
                )
                .unwrap()
            };
            let k31 = run(&ds.x);
            for arm in MATRIX_ARMS {
                let (x, _) = matrix_arm(arm.joined(), &ds, &block);
                assert_ne!(x.rows, ds.x.rows);
                let got = run(&x);
                assert_eq!(got.splits.len(), k31.splits.len());
                for (a, b) in k31.splits.iter().zip(&got.splits) {
                    assert_eq!(a.unavailable_reason, b.unavailable_reason, "{}", a.scheme);
                    assert_eq!(
                        a.test_fold,
                        b.test_fold,
                        "{} {} at seed {seed}",
                        arm.name(),
                        a.scheme
                    );
                    available += (a.unavailable_reason.is_none() && a.scheme == "lineage") as usize;
                }
            }
        }
        assert!(
            available > 0,
            "the fixture's lineage split runs at some seed"
        );
    }

    #[test]
    fn a_unitig_arm_sums_its_cohorts_pairs_and_refuses_a_build_of_other_genomes() {
        // SYNTHETIC: a unitig build over a species matrix of 50 genomes, the
        // cohort its first 30. Genome i holds i + 1 unitigs whole and 2(i + 1)
        // in part, so the cohort's sums are 465 and 930. A build whose genomes
        // are the species matrix's in another order is refused, as is one with
        // no row for a cohort genome.
        let (_, ds) = cohort(30, 4, 20, 37);
        let m = species_matrix(50, 40, 41);
        let species = m.summary.genomes.clone();
        let ids = species[..30].to_vec();
        let pairs: HashMap<String, (u64, u64)> = species
            .iter()
            .enumerate()
            .map(|(i, g)| (g.clone(), (i as u64 + 1, 2 * (i as u64 + 1))))
            .collect();
        let settings = Settings {
            arm: Arm::Unitigs,
            tables: None,
            matrix: None,
            column: Column::Symbol,
            rows: Rows::Amr,
            out: std::env::temp_dir(),
        };
        let build = Build {
            kind: "unitigs".to_string(),
            ..Build::default()
        };
        let copy = |ds: &Dataset| Dataset {
            x: ds.x.clone(),
            y: ds.y.clone(),
            columns: ds.columns.clone(),
            feature_offsets: ds.feature_offsets.clone(),
            feature_patterns: Vec::new(),
        };
        let added = Added {
            matrix: species_matrix(50, 40, 41),
            build: build.clone(),
            pairs: Some(pairs.clone()),
        };
        let (_, result) = matrix_dataset(
            &settings,
            &added,
            copy(&ds),
            "synthetic__drug",
            "drawn",
            &ids,
            &species,
            3,
        )
        .unwrap();
        assert_eq!(
            (result.unitig_pairs_held, result.unitig_pairs_in_part),
            (Some(465), Some(930))
        );
        let mut other = species_matrix(50, 40, 41);
        other.summary.genomes.reverse();
        let reordered = Added {
            matrix: other,
            build: build.clone(),
            pairs: Some(pairs.clone()),
        };
        assert!(
            matrix_dataset(
                &settings,
                &reordered,
                copy(&ds),
                "synthetic__drug",
                "drawn",
                &ids,
                &species,
                3
            )
            .is_err()
        );
        let mut short = pairs;
        short.remove(&ids[29]);
        let missing = Added {
            matrix: m,
            build,
            pairs: Some(short),
        };
        assert!(
            matrix_dataset(
                &settings,
                &missing,
                copy(&ds),
                "synthetic__drug",
                "drawn",
                &ids,
                &species,
                3
            )
            .is_err()
        );
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
                tables: Some(dir.clone()),
                matrix: None,
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
            assert_eq!(result.program_rows, Some(p));
            assert_eq!(result.k31_rows, k);
            assert_eq!(result.k31_rows_as_program, Some(0));
            let want = if arm == Arm::K31Program { p + k } else { p };
            assert_eq!((result.arm_rows, got.x.rows), (want, want));
            assert_eq!(result.tables_read, Some(40));
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
