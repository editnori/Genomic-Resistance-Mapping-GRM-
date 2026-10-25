//! research.json: revision-5 research results pinned for the dashboard.
//!
//! Two research lanes ran on the revision-5 snapshot whose hash list is
//! `RESEARCH-SNAPSHOT-SHA256SUMS`. The MIC-first lane scored GRM's MIC model
//! beside a no-genome baseline; the selective-calling lane asked on which
//! datasets the genome alone can make the R/S call for a confident subset.
//! `grm research` reads their saved outputs under `<data>/research/`, checks
//! them against the revision-5 dataset reports, and writes one file that the
//! dashboard serves beside the report.
//!
//! A dataset id is not a cohort: a later revision can keep the id and change
//! the genomes, their order, labels, MICs, or held-out scores. Each entry
//! therefore carries a guard of SHA-256 digests over exactly the report
//! values its lane read, taken from the revision-5 report file whose hash
//! the snapshot lists. The web app (web/src/data/research.ts) computes the
//! same digests from the report it loaded and uses an entry only when every
//! one matches. Each digest covers one line per value, each ending in "\n":
//!
//!   sample_ids      the sample id, in report order
//!   labels          R or S, then m when the label came from an MIC
//!   lineages        the lineage cluster as a decimal integer
//!   mic_values      - without an MIC, else "<censor> <log2>"
//!   mic_predicted   the saved predicted log2 MICs: - for null, else the
//!                   f64's bits as 16 lowercase hex digits; the guard holds
//!                   null when the report has no MIC model
//!   scores          for each split and learner in report order,
//!                   "#<scheme> <learner>" and then its held-out scores,
//!                   encoded as mic_predicted is
//!
//! The MIC lane read the samples, MICs, and saved MIC predictions; the
//! selective-calling lane read the samples and every split's held-out scores
//! (the lineage split for its verdicts, the random split for category C).
//! A MIC rerun's guard holds the MICs and random-split predictions the lane
//! scored, so a report matches it only where the rerun reproduced the
//! report's values.
//!
//! A MIC entry whose model scores disk-diffusion zone diameters as MICs
//! ships no numbers: `mic.withheld` lists its id and the reason, and the web
//! app shows that reason on every report revision.

use crate::curate::{Censor, Mic};
use anyhow::{Context, Result, bail, ensure};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::io::BufRead;
use std::path::Path;

pub const SCHEMA: u32 = 1;
/// SHA-256 of `RESEARCH-SNAPSHOT-SHA256SUMS`, the revision-5 snapshot both
/// lane receipts name.
pub const REVISION5_SNAPSHOT: &str =
    "75a724c1c362e6f5c9555886788c57cbf99f617d25f8a91dc8beb4c54b0b92fe";
const SNAPSHOT_SUMS: &str = "RESEARCH-SNAPSHOT-SHA256SUMS";
const MIC_LANE: &str = "grm-research-mic-first-r1";
const SELECTIVE_LANE: &str = "grm-research-selective-calling-r1";

#[derive(Serialize)]
pub struct Research {
    pub schema: u32,
    pub snapshot: Snapshot,
    /// Every lane file read, relative to the research directory.
    pub sources: Vec<Source>,
    pub mic: MicSection,
    pub selective: SelectiveSection,
}

#[derive(Serialize)]
pub struct Snapshot {
    pub revision: u32,
    pub sha256: String,
    /// Dataset reports the hash list covers.
    pub datasets: usize,
}

#[derive(Serialize)]
pub struct Source {
    pub path: String,
    pub sha256: String,
    /// Where a symlinked lane file points.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub target: Option<String>,
}

#[derive(Serialize, Clone, Debug, PartialEq)]
pub struct CohortGuard {
    pub samples: usize,
    pub sample_ids: String,
    pub labels: String,
    pub lineages: String,
}

#[derive(Serialize)]
pub struct MicGuard {
    #[serde(flatten)]
    pub cohort: CohortGuard,
    pub mic_values: String,
    pub mic_predicted: Option<String>,
}

#[derive(Serialize)]
pub struct SelectiveGuard {
    #[serde(flatten)]
    pub cohort: CohortGuard,
    pub scores: String,
}

#[derive(Serialize)]
pub struct MicSection {
    pub limits: MicLimits,
    pub datasets: Vec<MicDataset>,
    /// Datasets whose model scores zone diameters as MICs: no numbers, only
    /// the reason.
    pub withheld: Vec<MicWithheld>,
}

#[derive(Serialize, Debug, PartialEq)]
pub struct MicWithheld {
    pub id: String,
    pub reason: String,
}

#[derive(Serialize)]
pub struct MicLimits {
    pub scope: &'static str,
    pub baseline: &'static str,
    pub exact_minimum: &'static str,
    pub saved_only: &'static str,
    /// `{removed}` stands for the entry's removed zone diameters, which the
    /// web app writes in.
    pub zone_fix: &'static str,
}

// `MicInputs::dataset` refuses a rerun with a training fold that has no
// exact MIC, which keeps the baseline sentence true.
const MIC_LIMITS: MicLimits = MicLimits {
    scope: "Internal cross-validation on public BV-BRC genomes in the revision-5 snapshot. This is not external or clinical validation.",
    baseline: "The no-genome baseline predicts the median exact MIC of each test isolate's training fold.",
    exact_minimum: "GRM fits an MIC model only with at least 50 exact MICs, and this rerun scored fewer.",
    saved_only: "The research did not copy this species' feature matrix, which its rerun needs to refit the lineage split and rebuild the baseline's folds, so this row scores the report's saved random-split predictions.",
    zone_fix: "Revision 5 read {removed} as MICs; the rerun drops them and scores the MICs that remain.",
};

#[derive(Serialize)]
pub struct MicDataset {
    pub id: String,
    pub guard: MicGuard,
    pub zone_diameters: ZoneDiameters,
    pub model: MicModel,
    pub baseline: Option<MicBaseline>,
}

#[derive(Serialize, Debug, PartialEq)]
pub struct ZoneDiameters {
    /// Report MICs the rerun dropped because the zone-diameter fix removes them.
    pub removed: usize,
    /// Scored MICs that are zone diameters the fix would remove; `build`
    /// withholds an entry with any.
    pub scored: usize,
}

#[derive(Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ModelSource {
    Rerun,
    SavedReport,
}

#[derive(Serialize)]
pub struct MicModel {
    pub source: ModelSource,
    pub random: MicRun,
    pub lineage: Option<MicRun>,
}

#[derive(Serialize)]
pub struct MicBaseline {
    pub random: MicRun,
    pub lineage: MicRun,
}

/// One prediction set scored as in the lane's prediction_quality.csv.
#[derive(Serialize)]
pub struct MicRun {
    pub n_predicted: usize,
    pub n_exact: usize,
    pub ea_exact: Option<f64>,
    /// Lineage-split model runs: folds whose early-stopping holdout fell back
    /// to a random holdout.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub holdout_random_folds: Option<usize>,
}

#[derive(Serialize)]
pub struct SelectiveSection {
    pub categories: Vec<Category>,
    pub limits: SelectiveLimits,
    pub datasets: Vec<SelectiveDataset>,
}

#[derive(Serialize)]
pub struct Category {
    pub code: &'static str,
    /// What the category means for a dataset, in full sentences. A2's names
    /// its refit learner as `{learner}`, which the web app writes in.
    pub label: &'static str,
}

#[derive(Serialize)]
pub struct SelectiveLimits {
    pub scope: &'static str,
    /// The bar every category's pass or fail refers to (`fda_pass` in the
    /// lane's common.py).
    pub bar: &'static str,
    pub no_call: &'static str,
    pub optimism: &'static str,
    pub procedure: &'static str,
    pub per_100: &'static str,
}

// The counts in optimism are the lane summary's: 13 of the 228 datasets with
// a lineage split pass the primary rule and 15 each pass the gbt and rf
// refits, and A1 and A2 hold 13 and 17.
const SELECTIVE_LIMITS: SelectiveLimits = SelectiveLimits {
    scope: "Internal cross-validation on public genomes in the revision-5 snapshot, on the lineage split. This is not external or prospective validation, and there is no challenge set or essential agreement on MICs.",
    bar: "A rule passes when the isolates it calls meet the reference bar from FDA guidance for susceptibility test devices: category agreement above 89.9%, major error (ME, susceptible called resistant) under 3%, and a two-sided 95% Clopper-Pearson interval for very major error (VME, resistant called susceptible) with its lower bound under 1.5% and its upper bound under 7.5%.",
    no_call: "VME and ME count called isolates only, so a rule that calls fewer isolates can show lower rates, and a rule that calls none has no errors.",
    optimism: "A1 and A2 choose their rule inside the training folds, but trying three procedures makes A1 and A2 together mildly optimistic: each procedure alone passes 13 to 15 of the 228 datasets with a lineage split, and A1 and A2 together pass 30.",
    procedure: "The primary rule can pick a different learner in each fold, so an A1 result describes that procedure across folds.",
    per_100: "Per-100 values count every resistant or susceptible isolate in the dataset, so right calls, wrong calls, and no calls add up to 100 for each class. The resistant share comes from BV-BRC sampling, so it is not a hospital's resistance rate.",
};

/// Codes research.json uses for the lane's categories (species_drug.py):
/// the code, the lane's text after the code, and what the category means.
/// D1's 48 is the lane's power floor: with 47 called resistant isolates and
/// no misses the Clopper-Pearson upper bound is 7.55%, with 48 it is 7.40%.
const CATEGORIES: [(&str, &str, &str); 7] = [
    (
        "A1",
        "A1: genome-only calls on a confident subset (primary nested rule, lineage)",
        "The primary rule passes with lineages held out. It picks the learner and no-call band inside the training folds and was fixed before any held-out result was seen.",
    ),
    (
        "A2",
        "A2: genome-only calls on a confident subset (refit-nested gbt or rf, lineage)",
        "The primary rule fails, but {learner} passes with lineages held out when it is refit inside each training fold to choose its no-call band.",
    ),
    (
        "B",
        "B: passes only with a test-selected rule (unconfirmed)",
        "No rule chosen inside the training folds passes with lineages held out. A learner and no-call band picked on the held-out results pass, so the pass is optimistic until new isolates confirm it.",
    ),
    (
        "C",
        "C: only when close relatives are in training (random split)",
        "No rule passes with lineages held out, even one picked on the held-out results. A rule passes only on the random split, where test isolates have close relatives in training.",
    ),
    (
        "D1",
        "D1: never on this data - too few resistant isolates",
        "This dataset has fewer than 48 resistant isolates, too few for any rule to pass: calling every one of them with no misses still leaves the upper bound of the VME interval above 7.5%.",
    ),
    (
        "D2",
        "D2: never on this data - errors too frequent",
        "Selective calling ran on both splits, and no rule passes on either because errors among the called isolates are too frequent.",
    ),
    (
        "no_lineage",
        "no lineage split",
        "This report has no lineage split, so selective calling has no verdict for it.",
    ),
];

#[derive(Serialize)]
pub struct SelectiveDataset {
    pub id: String,
    pub guard: SelectiveGuard,
    pub category: &'static str,
    /// Present for A1 and A2 only.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub result: Option<SelectiveResult>,
}

/// Genome-only calls on the lineage split, from the summary's primary_
/// columns for A1 and refit_best_ columns for A2, which species_drug.py
/// copies from the lane's coverage_risk.csv or nested_refit.csv row. VME and
/// ME are among called isolates; their bounds are two-sided 95%
/// Clopper-Pearson intervals.
#[derive(Serialize)]
pub struct SelectiveResult {
    /// A2: the refit learner whose rule passes.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub learner: Option<String>,
    pub n: usize,
    pub n_r: usize,
    pub n_s: usize,
    pub tp: usize,
    #[serde(rename = "fn")]
    pub fn_: usize,
    pub fp: usize,
    pub tn: usize,
    pub coverage: f64,
    pub coverage_r: f64,
    pub coverage_s: f64,
    pub vme: f64,
    pub vme_lo: f64,
    pub vme_hi: f64,
    pub me: f64,
    pub me_lo: f64,
    pub me_hi: f64,
    pub ca: f64,
    pub per100r_caught: f64,
    pub per100r_missed: f64,
    pub per100r_lab: f64,
    pub per100s_cleared: f64,
    pub per100s_overcalled: f64,
    pub per100s_lab: f64,
}

// ---------------------------------------------------------------- digests

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn sha256(bytes: &[u8]) -> String {
    hex(ring::digest::digest(&ring::digest::SHA256, bytes).as_ref())
}

/// SHA-256 over each line followed by "\n".
fn digest_lines<S: AsRef<str>>(lines: impl IntoIterator<Item = S>) -> String {
    let mut ctx = ring::digest::Context::new(&ring::digest::SHA256);
    for line in lines {
        ctx.update(line.as_ref().as_bytes());
        ctx.update(b"\n");
    }
    hex(ctx.finish().as_ref())
}

fn value_line(v: Option<f64>) -> String {
    v.map_or_else(|| "-".to_string(), |x| format!("{:016x}", x.to_bits()))
}

fn mic_line(m: Option<Mic>) -> String {
    match m {
        None => "-".to_string(),
        Some(m) => {
            let censor = match m.censor {
                Censor::Exact => "Exact",
                Censor::AtMost => "AtMost",
                Censor::Above => "Above",
            };
            format!("{censor} {}", m.log2)
        }
    }
}

fn values_digest(values: &[Option<f64>]) -> String {
    digest_lines(values.iter().map(|v| value_line(*v)))
}

fn mics_digest(mics: impl Iterator<Item = Option<Mic>>) -> String {
    digest_lines(mics.map(mic_line))
}

/// The parts of a dataset report the guards cover.
#[derive(Deserialize)]
struct Report {
    id: String,
    antibiotic: String,
    counts: ReportCounts,
    samples: Vec<ReportSample>,
    mic: Option<SavedMic>,
    splits: Vec<ReportSplit>,
}

#[derive(Deserialize)]
struct ReportCounts {
    samples: usize,
    resistant: usize,
    susceptible: usize,
}

#[derive(Deserialize)]
struct ReportSample {
    id: String,
    resistant: bool,
    from_mic: bool,
    mic: Option<Mic>,
    lineage: u32,
}

#[derive(Deserialize)]
struct SavedMic {
    predicted: Vec<Option<f64>>,
}

#[derive(Deserialize)]
struct ReportSplit {
    scheme: String,
    learners: Vec<ReportLearner>,
}

#[derive(Deserialize)]
struct ReportLearner {
    learner: String,
    scores: Vec<Option<f64>>,
}

fn cohort_guard(samples: &[ReportSample]) -> Result<CohortGuard> {
    if let Some(s) = samples.iter().find(|s| s.id.contains('\n')) {
        bail!("sample id {:?} contains a newline", s.id);
    }
    Ok(CohortGuard {
        samples: samples.len(),
        sample_ids: digest_lines(samples.iter().map(|s| s.id.as_str())),
        labels: digest_lines(samples.iter().map(|s| match (s.resistant, s.from_mic) {
            (true, false) => "R",
            (true, true) => "Rm",
            (false, false) => "S",
            (false, true) => "Sm",
        })),
        lineages: digest_lines(samples.iter().map(|s| s.lineage.to_string())),
    })
}

fn scores_digest(report: &Report) -> String {
    digest_lines(report.splits.iter().flat_map(|split| {
        split.learners.iter().flat_map(|l| {
            std::iter::once(format!("#{} {}", split.scheme, l.learner))
                .chain(l.scores.iter().map(|v| value_line(*v)))
        })
    }))
}

// -------------------------------------------------------------- lane files

/// The research directory, recording the hash of every file read.
struct Inputs<'a> {
    root: &'a Path,
    sources: Vec<Source>,
}

impl Inputs<'_> {
    fn read(&mut self, rel: &str) -> Result<Vec<u8>> {
        let path = self.root.join(rel);
        let bytes = std::fs::read(&path).with_context(|| format!("read {}", path.display()))?;
        let target = if std::fs::symlink_metadata(&path)?.file_type().is_symlink() {
            Some(std::fs::read_link(&path)?.to_string_lossy().into_owned())
        } else {
            None
        };
        self.sources.push(Source {
            path: rel.to_string(),
            sha256: sha256(&bytes),
            target,
        });
        Ok(bytes)
    }
}

/// A CSV file as the lanes write it: a header row, quoted fields that may
/// hold commas, and no line breaks inside fields.
struct Table {
    name: String,
    columns: HashMap<String, usize>,
    rows: Vec<Vec<String>>,
}

fn csv_fields(line: &str) -> Result<Vec<String>> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut quoted = false;
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        match (quoted, c) {
            (true, '"') if chars.peek() == Some(&'"') => {
                field.push('"');
                chars.next();
            }
            (true, '"') => quoted = false,
            (false, '"') if field.is_empty() => quoted = true,
            (false, ',') => fields.push(std::mem::take(&mut field)),
            (_, c) => field.push(c),
        }
    }
    ensure!(!quoted, "unterminated quoted field");
    fields.push(field);
    Ok(fields)
}

impl Table {
    fn parse(name: &str, bytes: &[u8]) -> Result<Table> {
        let text = std::str::from_utf8(bytes).with_context(|| format!("{name} is not UTF-8"))?;
        let mut lines = text.lines();
        let header = csv_fields(lines.next().with_context(|| format!("{name} is empty"))?)?;
        let rows = lines
            .enumerate()
            .map(|(i, line)| {
                let row = csv_fields(line).with_context(|| format!("{name} line {}", i + 2))?;
                ensure!(
                    row.len() == header.len(),
                    "{name} line {} has {} fields for {} columns",
                    i + 2,
                    row.len(),
                    header.len()
                );
                Ok(row)
            })
            .collect::<Result<_>>()?;
        Ok(Table {
            name: name.to_string(),
            columns: header
                .into_iter()
                .enumerate()
                .map(|(i, c)| (c, i))
                .collect(),
            rows,
        })
    }

    fn cell<'a>(&self, row: &'a [String], column: &str) -> Result<&'a str> {
        let i = self
            .columns
            .get(column)
            .with_context(|| format!("{} has no column {column}", self.name))?;
        Ok(row[*i].as_str())
    }

    fn real(&self, row: &[String], column: &str) -> Result<Option<f64>> {
        match self.cell(row, column)? {
            "" | "nan" => Ok(None),
            s => {
                let v: f64 = s
                    .parse()
                    .with_context(|| format!("{} {column}: {s:?} is not a number", self.name))?;
                ensure!(v.is_finite(), "{} {column}: {s:?} is not finite", self.name);
                Ok(Some(v))
            }
        }
    }

    fn number(&self, row: &[String], column: &str) -> Result<f64> {
        self.real(row, column)?
            .with_context(|| format!("{} {column} is empty", self.name))
    }

    /// pandas writes counts as floats ("835.0") once a column has a gap.
    fn count(&self, row: &[String], column: &str) -> Result<usize> {
        let v = self.number(row, column)?;
        ensure!(
            v >= 0.0 && v.fract() == 0.0 && v < 1e15,
            "{} {column}: {v} is not a count",
            self.name
        );
        Ok(v as usize)
    }

    fn flag(&self, row: &[String], column: &str) -> Result<Option<bool>> {
        match self.cell(row, column)? {
            "True" => Ok(Some(true)),
            "False" => Ok(Some(false)),
            "" => Ok(None),
            s => bail!("{} {column}: {s:?} is not True or False", self.name),
        }
    }
}

fn parse_mic(log2: &str, censor: &str) -> Result<Option<Mic>> {
    let censor = match censor {
        "" => return Ok(None),
        "Exact" => Censor::Exact,
        "AtMost" => Censor::AtMost,
        "Above" => Censor::Above,
        c => bail!("unknown MIC censor {c:?}"),
    };
    let log2 = log2.parse().with_context(|| format!("MIC log2 {log2:?}"))?;
    Ok(Some(Mic { log2, censor }))
}

/// (snapshot MIC, post-fix MIC) for every genome-drug pair whose MIC the
/// zone-diameter fix changes, from the lane's pairs table.
type ZoneChanges = HashMap<(String, String), (Option<Mic>, Option<Mic>)>;

fn zone_changes(gz: &[u8]) -> Result<ZoneChanges> {
    let reader = std::io::BufReader::new(flate2::read::MultiGzDecoder::new(gz));
    let mut lines = reader.lines();
    let header: Vec<String> = lines
        .next()
        .context("pairs table is empty")??
        .split('\t')
        .map(String::from)
        .collect();
    let col = |name: &str| {
        header
            .iter()
            .position(|c| c == name)
            .with_context(|| format!("pairs table has no column {name}"))
    };
    let [genome, drug, snap_log2, snap_censor, log2, censor] = [
        "genome_id",
        "antibiotic",
        "snap_mic_log2",
        "snap_mic_censor",
        "mic_log2",
        "mic_censor",
    ]
    .map(col);
    let (genome, drug, snap_log2, snap_censor, log2, censor) =
        (genome?, drug?, snap_log2?, snap_censor?, log2?, censor?);
    let mut changed = HashMap::new();
    for line in lines {
        let line = line?;
        let c: Vec<&str> = line.split('\t').collect();
        ensure!(
            c.len() == header.len(),
            "pairs table row has {} fields",
            c.len()
        );
        let snap = parse_mic(c[snap_log2], c[snap_censor])?;
        let post = parse_mic(c[log2], c[censor])?;
        if snap != post {
            changed.insert((c[genome].to_string(), c[drug].to_string()), (snap, post));
        }
    }
    Ok(changed)
}

/// The parts of a micexp output file (`out/splits_final/<id>.json`) the
/// checks read.
#[derive(Deserialize)]
struct LaneSplits {
    id: String,
    samples: Vec<LaneSample>,
    runs: Vec<LaneRun>,
}

#[derive(Deserialize)]
struct LaneSample {
    id: String,
    resistant: bool,
    mic: Option<Mic>,
    lineage: u32,
    extra: bool,
}

#[derive(Deserialize)]
struct LaneRun {
    scheme: String,
    augmented: bool,
    predicted: Vec<Option<f64>>,
    holdout_kind: Vec<String>,
    /// Each sample's outer test fold; null for a sample the run did not test.
    folds: Vec<Option<usize>>,
}

// ---------------------------------------------------------------- MIC lane

struct MicInputs {
    quality: Table,
    /// prediction_quality.csv row indices per dataset id.
    rows: BTreeMap<String, Vec<usize>>,
    changed: ZoneChanges,
}

impl MicInputs {
    fn read(inputs: &mut Inputs) -> Result<MicInputs> {
        let quality = Table::parse(
            "prediction_quality.csv",
            &inputs.read(&format!("{MIC_LANE}/out/analysis/prediction_quality.csv"))?,
        )?;
        let mut rows: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (i, row) in quality.rows.iter().enumerate() {
            let source = quality.cell(row, "source")?;
            if matches!(
                source,
                "saved_report" | "splits_final" | "splits_final/null_median"
            ) {
                ensure!(
                    quality.flag(row, "augmented")? == Some(false),
                    "prediction_quality.csv: augmented {source} row"
                );
                rows.entry(quality.cell(row, "id")?.to_string())
                    .or_default()
                    .push(i);
            }
        }
        let changed = zone_changes(&inputs.read(&format!("{MIC_LANE}/out/pairs/pairs.tsv.gz"))?)?;
        Ok(MicInputs {
            quality,
            rows,
            changed,
        })
    }

    fn pick(&self, id: &str, source: &str, scheme: &str) -> Result<Vec<&[String]>> {
        let mut out = Vec::new();
        for &i in self.rows.get(id).map_or(&[][..], Vec::as_slice) {
            let row = &self.quality.rows[i];
            if self.quality.cell(row, "source")? == source
                && self.quality.cell(row, "scheme")? == scheme
            {
                out.push(row.as_slice());
            }
        }
        Ok(out)
    }

    /// One prediction set across the standards; essential agreement does not
    /// depend on the standard, so every row must carry the same values.
    fn run(&self, id: &str, source: &str, scheme: &str, holdout: Option<usize>) -> Result<MicRun> {
        let rows = self.pick(id, source, scheme)?;
        let q = &self.quality;
        let first = *rows
            .first()
            .with_context(|| format!("{id}: no {source} {scheme} rows"))?;
        let ea = |r: &[String]| -> Result<_> {
            Ok((
                q.count(r, "n_predicted")?,
                q.count(r, "n_exact")?,
                q.real(r, "ea_exact")?,
            ))
        };
        let (n_predicted, n_exact, ea_exact) = ea(first)?;
        for &r in &rows {
            ensure!(
                ea(r)? == (n_predicted, n_exact, ea_exact),
                "{id}: {source} {scheme} rows disagree on essential agreement"
            );
            if holdout.is_some() {
                let kinds = q.cell(r, "lineage_holdouts")?;
                ensure!(
                    Some(kinds.split(',').filter(|k| *k == "random").count()) == holdout,
                    "{id}: lineage_holdouts {kinds:?} differs from the micexp file"
                );
            }
        }
        Ok(MicRun {
            n_predicted,
            n_exact,
            ea_exact,
            holdout_random_folds: holdout,
        })
    }

    /// Scored MICs that are the pre-fix value of a pair the zone fix changes.
    fn zone_scored<'a>(
        &self,
        antibiotic: &str,
        scored: impl Iterator<Item = (&'a str, Option<Mic>)>,
    ) -> usize {
        scored
            .filter(|(id, m)| {
                m.is_some()
                    && self
                        .changed
                        .get(&(id.to_string(), antibiotic.to_string()))
                        .is_some_and(|(snap, _)| snap == m)
            })
            .count()
    }

    fn dataset(
        &self,
        inputs: &mut Inputs,
        report: &Report,
        cohort: &CohortGuard,
    ) -> Result<MicDataset> {
        let id = report.id.as_str();
        if self.pick(id, "splits_final", "random")?.is_empty() {
            return self.saved_only(report, cohort);
        }
        let rel = format!("{MIC_LANE}/out/splits_final/{id}.json");
        let lane: LaneSplits =
            serde_json::from_slice(&inputs.read(&rel)?).with_context(|| format!("parse {rel}"))?;
        let target = inputs.sources.last().and_then(|s| s.target.clone());
        let post_fix = match target.as_deref() {
            Some(t) if t.starts_with("../splits_postfix/") => true,
            Some(t) if t.starts_with("../splits/") => false,
            t => bail!("{rel} should link into splits/ or splits_postfix/, not {t:?}"),
        };
        ensure!(lane.id == id, "{rel} holds dataset {}", lane.id);
        ensure!(
            lane.samples.len() == report.samples.len(),
            "{id}: the lane scored {} samples, the report has {}",
            lane.samples.len(),
            report.samples.len()
        );
        let mut removed = 0;
        for (i, (a, b)) in lane.samples.iter().zip(&report.samples).enumerate() {
            ensure!(
                a.id == b.id && a.resistant == b.resistant && a.lineage == b.lineage && !a.extra,
                "{id}: lane sample {i} ({}) is not the report's sample {} with the same label and lineage",
                a.id,
                b.id
            );
            if a.mic == b.mic {
                continue;
            }
            let fix = self.changed.get(&(b.id.clone(), report.antibiotic.clone()));
            ensure!(
                post_fix && a.mic.is_none() && fix == Some(&(b.mic, None)),
                "{id}: sample {} has MIC {:?} in the lane and {:?} in the report, which the zone fix does not explain",
                b.id,
                a.mic,
                b.mic
            );
            removed += 1;
        }
        let run = |scheme: &str| {
            let run = lane
                .runs
                .iter()
                .find(|r| r.scheme == scheme && !r.augmented)
                .with_context(|| format!("{rel} has no {scheme} run"))?;
            ensure!(
                run.predicted.len() == lane.samples.len() && run.folds.len() == lane.samples.len(),
                "{rel}: short {scheme} run"
            );
            // MIC_LIMITS.baseline says the baseline predicts each training
            // fold's median exact MIC; micfirst.py falls back to censored
            // values for a training fold with none, so refuse that fold.
            let mut exact = BTreeMap::new();
            for (s, fold) in lane.samples.iter().zip(&run.folds) {
                if let Some(fold) = fold {
                    *exact.entry(*fold).or_insert(0) +=
                        usize::from(!s.extra && s.mic.is_some_and(|m| m.censor == Censor::Exact));
                }
            }
            let total: usize = exact.values().sum();
            ensure!(
                exact.values().all(|&n| n < total),
                "{rel}: a {scheme} training fold has no exact MIC, so its baseline is not a median exact MIC"
            );
            Ok(run)
        };
        let (random, lineage) = (run("random")?, run("lineage")?);
        let holdout = lineage
            .holdout_kind
            .iter()
            .filter(|k| *k == "random")
            .count();
        Ok(MicDataset {
            id: id.to_string(),
            guard: MicGuard {
                cohort: cohort.clone(),
                mic_values: mics_digest(lane.samples.iter().map(|s| s.mic)),
                // The predictions the lane scored. The web digests the
                // report's saved predictions, so this matches only where the
                // rerun reproduced them, and is null like the report's when
                // the report has no MIC model.
                mic_predicted: report
                    .mic
                    .as_ref()
                    .map(|_| values_digest(&random.predicted)),
            },
            zone_diameters: ZoneDiameters {
                removed,
                scored: self.zone_scored(
                    &report.antibiotic,
                    lane.samples.iter().map(|s| (s.id.as_str(), s.mic)),
                ),
            },
            model: MicModel {
                source: ModelSource::Rerun,
                random: self.run(id, "splits_final", "random", None)?,
                lineage: Some(self.run(id, "splits_final", "lineage", Some(holdout))?),
            },
            baseline: Some(MicBaseline {
                random: self.run(id, "splits_final/null_median", "random", None)?,
                lineage: self.run(id, "splits_final/null_median", "lineage", None)?,
            }),
        })
    }

    fn saved_only(&self, report: &Report, cohort: &CohortGuard) -> Result<MicDataset> {
        let id = report.id.as_str();
        for (source, scheme) in [
            ("splits_final", "lineage"),
            ("splits_final/null_median", "random"),
            ("splits_final/null_median", "lineage"),
            ("saved_report", "lineage"),
        ] {
            ensure!(
                self.pick(id, source, scheme)?.is_empty(),
                "{id}: {source} {scheme} rows without a splits_final random run"
            );
        }
        let saved = report
            .mic
            .as_ref()
            .with_context(|| format!("{id}: saved_report rows but the report has no MIC model"))?;
        Ok(MicDataset {
            id: id.to_string(),
            guard: MicGuard {
                cohort: cohort.clone(),
                mic_values: mics_digest(report.samples.iter().map(|s| s.mic)),
                mic_predicted: Some(values_digest(&saved.predicted)),
            },
            zone_diameters: ZoneDiameters {
                removed: 0,
                scored: self.zone_scored(
                    &report.antibiotic,
                    report.samples.iter().map(|s| (s.id.as_str(), s.mic)),
                ),
            },
            model: MicModel {
                source: ModelSource::SavedReport,
                random: self.run(id, "saved_report", "random", None)?,
                lineage: None,
            },
            baseline: None,
        })
    }
}

/// An entry whose model scores zone diameters as MICs keeps only its id and
/// the reason it ships no numbers. The reason carries no revision-5 count,
/// because the web app shows it on every report revision.
fn withhold(d: MicDataset) -> MicWithheld {
    let model = match d.model.source {
        ModelSource::SavedReport => "saved",
        ModelSource::Rerun => "rerun",
    };
    MicWithheld {
        reason: format!(
            "The revision-5 {model} MIC model scored disk-diffusion zone diameters as MICs, and there is no rerun or no-genome baseline without them."
        ),
        id: d.id,
    }
}

// ---------------------------------------------------------- selective lane

struct SelectiveInputs {
    summary: Table,
    /// Summary row index per dataset.
    rows: BTreeMap<String, usize>,
}

impl SelectiveInputs {
    fn read(inputs: &mut Inputs) -> Result<SelectiveInputs> {
        let summary = Table::parse(
            "species_drug_summary.csv",
            &inputs.read(&format!("{SELECTIVE_LANE}/out/species_drug_summary.csv"))?,
        )?;
        let rows = summary
            .rows
            .iter()
            .enumerate()
            .map(|(i, r)| Ok((summary.cell(r, "dataset")?.to_string(), i)))
            .collect::<Result<_>>()?;
        Ok(SelectiveInputs { summary, rows })
    }

    fn dataset(&self, report: &Report, cohort: &CohortGuard) -> Result<SelectiveDataset> {
        let id = report.id.as_str();
        let s = &self.summary;
        let row = &s.rows[self.rows[id]];
        let label = s.cell(row, "category")?;
        let code = CATEGORIES
            .iter()
            .find(|(_, l, _)| *l == label)
            .map(|(c, _, _)| *c)
            .with_context(|| format!("{id}: unknown category {label:?}"))?;
        let (n, n_r, n_s) = (
            s.count(row, "n")?,
            s.count(row, "n_r")?,
            s.count(row, "n_s")?,
        );
        let c = &report.counts;
        ensure!(
            (n, n_r, n_s) == (c.samples, c.resistant, c.susceptible),
            "{id}: the summary counts {n} isolates ({n_r} R, {n_s} S); the report {} ({} R, {} S)",
            c.samples,
            c.resistant,
            c.susceptible
        );
        let lineage = report
            .splits
            .iter()
            .any(|s| s.scheme == "lineage" && !s.learners.is_empty());
        ensure!(
            (code == "no_lineage") == !lineage,
            "{id}: category {code}, but the report has {} lineage split",
            if lineage { "a" } else { "no" }
        );
        let result = match code {
            "A1" => Some(self.result(id, row, "primary", None, (n, n_r, n_s))?),
            "A2" => {
                let learner = s.cell(row, "refit_best_learner")?;
                ensure!(
                    matches!(learner, "gbt" | "rf"),
                    "{id}: an A2 result needs refit_best_learner gbt or rf, not {learner:?}"
                );
                Some(self.result(id, row, "refit_best", Some(learner), (n, n_r, n_s))?)
            }
            _ => None,
        };
        Ok(SelectiveDataset {
            id: id.to_string(),
            guard: SelectiveGuard {
                cohort: cohort.clone(),
                scores: scores_digest(report),
            },
            category: code,
            result,
        })
    }

    /// The passing rule's values, which species_drug.py copied into the
    /// summary under `prefix`. Their counts must partition the isolates, as
    /// the dashboard's parser requires.
    fn result(
        &self,
        id: &str,
        row: &[String],
        prefix: &str,
        learner: Option<&str>,
        (n, n_r, n_s): (usize, usize, usize),
    ) -> Result<SelectiveResult> {
        let s = &self.summary;
        let count = |k: &str| s.count(row, &format!("{prefix}_{k}"));
        let x = |k: &str| s.number(row, &format!("{prefix}_{k}"));
        let r = SelectiveResult {
            learner: learner.map(String::from),
            n,
            n_r,
            n_s,
            tp: count("tp")?,
            fn_: count("fn")?,
            fp: count("fp")?,
            tn: count("tn")?,
            coverage: x("coverage")?,
            coverage_r: x("coverage_r")?,
            coverage_s: x("coverage_s")?,
            vme: x("vme")?,
            vme_lo: x("vme_lo")?,
            vme_hi: x("vme_hi")?,
            me: x("me")?,
            me_lo: x("me_lo")?,
            me_hi: x("me_hi")?,
            ca: x("ca")?,
            per100r_caught: x("per100r_caught")?,
            per100r_missed: x("per100r_missed")?,
            per100r_lab: x("per100r_lab")?,
            per100s_cleared: x("per100s_cleared")?,
            per100s_overcalled: x("per100s_overcalled")?,
            per100s_lab: x("per100s_lab")?,
        };
        ensure!(
            n == n_r + n_s && r.tp + r.fn_ <= n_r && r.fp + r.tn <= n_s,
            "{id}: the {prefix}_ counts do not partition the dataset"
        );
        Ok(r)
    }
}

// ------------------------------------------------------------------- build

/// Reads both lanes under `research` and the revision-5 reports in
/// `revision5_reports`, which must match the snapshot hash list, and
/// assembles research.json. `snapshot_sha256` pins the hash list.
pub fn build(research: &Path, revision5_reports: &Path, snapshot_sha256: &str) -> Result<Research> {
    let mut inputs = Inputs {
        root: research,
        sources: Vec::new(),
    };
    let sums = inputs.read(SNAPSHOT_SUMS)?;
    let sums_sha256 = sha256(&sums);
    ensure!(
        sums_sha256 == snapshot_sha256,
        "{SNAPSHOT_SUMS} has sha256 {sums_sha256}; the pinned revision-5 snapshot is {snapshot_sha256}"
    );
    let mut listed = BTreeMap::new();
    for line in std::str::from_utf8(&sums)?.lines() {
        let (sha, path) = line
            .split_once("  ")
            .with_context(|| format!("{SNAPSHOT_SUMS}: {line:?}"))?;
        if let Some(id) = path
            .strip_prefix("data/report/datasets/")
            .and_then(|p| p.strip_suffix(".json"))
        {
            listed.insert(id.to_string(), sha.to_string());
        }
    }
    for lane in [MIC_LANE, SELECTIVE_LANE] {
        let receipt = inputs.read(&format!("{lane}/receipt.md"))?;
        ensure!(
            String::from_utf8_lossy(&receipt).contains(snapshot_sha256),
            "{lane}/receipt.md does not name snapshot {snapshot_sha256}"
        );
    }
    let mic_inputs = MicInputs::read(&mut inputs)?;
    let selective_inputs = SelectiveInputs::read(&mut inputs)?;
    for id in mic_inputs.rows.keys().chain(selective_inputs.rows.keys()) {
        ensure!(
            listed.contains_key(id),
            "{id} has lane rows but no report in the snapshot"
        );
    }

    let (mut mic, mut withheld, mut selective) = (Vec::new(), Vec::new(), Vec::new());
    for (id, listed_sha256) in &listed {
        let in_mic = mic_inputs.rows.contains_key(id);
        let in_selective = selective_inputs.rows.contains_key(id);
        if !in_mic && !in_selective {
            continue;
        }
        let path = revision5_reports.join(format!("{id}.json"));
        let bytes = std::fs::read(&path).with_context(|| format!("read {}", path.display()))?;
        ensure!(
            sha256(&bytes) == *listed_sha256,
            "{} does not match the snapshot hash list",
            path.display()
        );
        let report: Report =
            serde_json::from_slice(&bytes).with_context(|| format!("parse {}", path.display()))?;
        ensure!(
            report.id == *id,
            "{} holds dataset {}",
            path.display(),
            report.id
        );
        ensure!(
            report.counts.samples == report.samples.len(),
            "{id}: counts.samples differs from the sample list"
        );
        let cohort = cohort_guard(&report.samples).with_context(|| id.clone())?;
        if in_mic {
            let d = mic_inputs.dataset(&mut inputs, &report, &cohort)?;
            if d.zone_diameters.scored > 0 {
                withheld.push(withhold(d));
            } else {
                mic.push(d);
            }
        }
        if in_selective {
            selective.push(selective_inputs.dataset(&report, &cohort)?);
        }
    }

    let categories = CATEGORIES
        .iter()
        .map(|&(code, _, label)| Category { code, label })
        .collect();
    Ok(Research {
        schema: SCHEMA,
        snapshot: Snapshot {
            revision: 5,
            sha256: sums_sha256,
            datasets: listed.len(),
        },
        sources: inputs.sources,
        mic: MicSection {
            limits: MIC_LIMITS,
            datasets: mic,
            withheld,
        },
        selective: SelectiveSection {
            categories,
            limits: SELECTIVE_LIMITS,
            datasets: selective,
        },
    })
}

pub fn write(research: &Research, out: &Path) -> Result<()> {
    crate::layout::ensure_parent(out)?;
    crate::report::publish(out, &serde_json::to_vec(research)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The same report and digests as web/src/data/research.test.ts; the
    /// digests were computed separately with Python's hashlib.
    const GOLDEN: &str = r#"{"id":"g","antibiotic":"x","counts":{"samples":3,"resistant":1,"susceptible":2},
        "samples":[
          {"id":"562.1","resistant":true,"from_mic":false,"mic":{"log2":3,"censor":"Above"},"lineage":0},
          {"id":"562.22","resistant":false,"from_mic":true,"mic":null,"lineage":2},
          {"id":"562.3","resistant":false,"from_mic":false,"mic":{"log2":-2,"censor":"Exact"},"lineage":0}],
        "mic":{"predicted":[2.75,null,-1.5]},
        "splits":[
          {"scheme":"random","learners":[{"learner":"gbt","scores":[0.5,0.5,0.5]}]},
          {"scheme":"lineage","learners":[
            {"learner":"gbt","scores":[0.9,0.1,0.25]},
            {"learner":"catalog","scores":[1,0,0]}]}]}"#;

    #[test]
    fn guard_digests_match_the_golden_vectors() {
        let r: Report = serde_json::from_str(GOLDEN).unwrap();
        let g = cohort_guard(&r.samples).unwrap();
        assert_eq!(g.samples, 3);
        assert_eq!(
            g.sample_ids,
            "04d55a6c086cda50b8411babb85468ef748eeffccff37e86131d6292312563b1"
        );
        assert_eq!(
            g.labels,
            "db5593b64018f1438d8c6813be1a0cfa5f306815b3b1f4f61614f1719599d794"
        );
        assert_eq!(
            g.lineages,
            "04ed82dae8b902d411accf765f4ded9af83a46b71ffb1fce5f02c611908e6bcd"
        );
        assert_eq!(
            mics_digest(r.samples.iter().map(|s| s.mic)),
            "3937937471419ea303f03db2232fe979f5b4fcd79ee03207ba1b6394273b8113"
        );
        assert_eq!(
            values_digest(&r.mic.as_ref().unwrap().predicted),
            "e9f5ea58df3eacbc7f595ba68475d7d67ae83d4d96bfc787bfad4bc25e66b577"
        );
        assert_eq!(
            scores_digest(&r),
            "40adb03707afe2d0fe731c7e8e48010a481daf5a4a8fc8efc00722ce5a245535"
        );
    }

    #[test]
    fn csv_fields_keep_quoted_commas() {
        assert_eq!(
            csv_fields(r#"a,"0,0,0",,"say ""hi""",b"#).unwrap(),
            ["a", "0,0,0", "", r#"say "hi""#, "b"]
        );
        assert!(csv_fields(r#"a,"open"#).is_err());
    }

    // A four-dataset research directory in the lanes' file formats:
    // E. coli ampicillin is rerun after the zone fix and passes the primary
    // rule (A1); K. pneumoniae meropenem is rerun unchanged and passes the
    // gbt refit (A2); K. pneumoniae ciprofloxacin is rerun, but its report
    // has no MIC model; E. faecium vancomycin has only its saved model,
    // scores one zone diameter, and passes only a test-selected band (B).
    use serde_json::{Value, json};
    use std::path::PathBuf;

    const AMP: &str = "escherichia_coli__ampicillin";
    const MERO: &str = "klebsiella_pneumoniae__meropenem";
    const CIP: &str = "klebsiella_pneumoniae__ciprofloxacin";
    const VAN: &str = "enterococcus_faecium__vancomycin";
    const AMP_REPORT_MICS: [Option<i32>; 4] = [Some(3), Some(5), None, Some(-1)];
    const AMP_LANE_MICS: [Option<i32>; 4] = [Some(3), None, None, Some(-1)];
    const MERO_MICS: [Option<i32>; 4] = [Some(3), Some(4), Some(0), Some(-1)];
    const MERO_SAVED: [Option<f64>; 4] = [Some(2.5), Some(3.25), Some(0.5), Some(-0.75)];
    const AMP_RERUN: [Option<f64>; 4] = [Some(2.0), None, None, Some(-1.25)];
    const CIP_MICS: [Option<i32>; 4] = [Some(2), Some(1), Some(-2), Some(-3)];
    const CIP_RERUN: [Option<f64>; 4] = [Some(1.5), Some(0.5), Some(-1.0), Some(-2.5)];

    /// The fields species_drug.py copies from a coverage or refit row into
    /// the summary under a primary_ or refit_best_ prefix.
    const RESULT_FIELDS: [&str; 20] = [
        "coverage",
        "coverage_r",
        "coverage_s",
        "vme",
        "vme_lo",
        "vme_hi",
        "me",
        "me_lo",
        "me_hi",
        "ca",
        "per100r_caught",
        "per100r_missed",
        "per100r_lab",
        "per100s_cleared",
        "per100s_overcalled",
        "per100s_lab",
        "tp",
        "fn",
        "fp",
        "tn",
    ];

    struct Fixture {
        _dir: tempfile::TempDir,
        research: PathBuf,
        reports: PathBuf,
        snapshot: String,
    }

    impl Fixture {
        fn mic_lane(&self, rel: &str) -> PathBuf {
            self.research.join(MIC_LANE).join(rel)
        }

        fn selective_lane(&self, rel: &str) -> PathBuf {
            self.research.join(SELECTIVE_LANE).join(rel)
        }

        fn build(&self) -> Result<Research> {
            build(&self.research, &self.reports, &self.snapshot)
        }

        fn error(&self) -> String {
            format!("{:#}", self.build().err().expect("build should fail"))
        }
    }

    fn mic_json(log2: Option<i32>) -> Value {
        log2.map_or(Value::Null, |l| json!({"log2": l, "censor": "Exact"}))
    }

    /// Four isolates, the first two resistant, in two lineages.
    fn samples(id: &str, mics: [Option<i32>; 4], lane: bool) -> Vec<Value> {
        (0..4)
            .map(|i| {
                let mut s = json!({"id": format!("{id}.{i}"), "resistant": i < 2,
                    "mic": mic_json(mics[i]), "lineage": i % 2});
                s[if lane { "extra" } else { "from_mic" }] = json!(false);
                s
            })
            .collect()
    }

    fn report_json(id: &str, mics: [Option<i32>; 4], predicted: Option<[Option<f64>; 4]>) -> Value {
        let learner = |name: &str, scores: [f64; 4]| json!({"learner": name, "scores": scores});
        json!({"id": id, "antibiotic": id.split("__").nth(1).unwrap(),
            "counts": {"samples": 4, "resistant": 2, "susceptible": 2},
            "samples": samples(id, mics, false),
            "mic": predicted.map(|p| json!({"predicted": p})),
            "splits": [
                {"scheme": "random", "learners": [learner("gbt", [0.9, 0.8, 0.2, 0.1])]},
                {"scheme": "lineage", "learners": [learner("gbt", [0.7, 0.6, 0.3, 0.2]),
                    learner("rf", [0.95, 0.55, 0.45, 0.05])]}]})
    }

    fn lane_json(id: &str, mics: [Option<i32>; 4], random: [Option<f64>; 4]) -> Value {
        // One fold per lineage. Each holds an isolate with an exact MIC in
        // every dataset here, so each training fold has one.
        let run = |scheme: &str, predicted: [Option<f64>; 4]| {
            json!({"scheme": scheme, "augmented": false, "predicted": predicted,
                "holdout_kind": vec![scheme; 5], "folds": [0, 1, 0, 1]})
        };
        json!({"id": id, "samples": samples(id, mics, true),
            "runs": [run("random", random), run("lineage", [Some(1.0), Some(2.0), None, Some(0.0)])]})
    }

    /// prediction_quality.csv rows for one prediction set: CLSI with a
    /// breakpoint and FDA without.
    fn quality_rows(id: &str, source: &str, scheme: &str, ea: f64) -> String {
        let holdouts = if source == "saved_report" {
            String::new()
        } else {
            format!("\"{}\"", [scheme; 5].join(","))
        };
        let lead = format!("{id},{source},{scheme},False");
        format!(
            "{lead},CLSI,True,3,2,{ea},0.75,3.0,1.0,0.0,2.0,0.6666666666666666,1.0,1.0,0.025,1.0,0.0,0.0,0.0,0.0,False,False,True,False,{holdouts}\n\
             {lead},FDA,False,3,2,{ea},0.75{}{holdouts}\n",
            ",".repeat(18)
        )
    }

    /// The 20 summary fields of a rule that calls tp + fn of 2 resistant and
    /// fp + tn of 2 susceptible isolates, computed as common.metrics does.
    fn result_cells(tp: usize, fn_: usize, fp: usize, tn: usize) -> String {
        let div = |a: usize, b: usize| a as f64 / b as f64;
        let (rc, sc) = (tp + fn_, fp + tn);
        [
            div(rc + sc, 4),
            div(rc, 2),
            div(sc, 2),
            div(fn_, rc),
            0.0,
            0.05,
            div(fp, sc),
            0.0,
            0.05,
            div(tp + tn, rc + sc),
            div(100 * tp, 2),
            div(100 * fn_, 2),
            div(100 * (2 - rc), 2),
            div(100 * tn, 2),
            div(100 * fp, 2),
            div(100 * (2 - sc), 2),
            tp as f64,
            fn_ as f64,
            fp as f64,
            tn as f64,
        ]
        .map(|v| v.to_string())
        .join(",")
    }

    fn write(path: PathBuf, bytes: impl AsRef<[u8]>) {
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        std::fs::write(path, bytes).unwrap();
    }

    /// The gzipped pairs table: the zone fix removes one ampicillin and one
    /// vancomycin MIC and leaves a meropenem MIC as it was; `extra` adds rows.
    fn pairs(extra: &str) -> Vec<u8> {
        use std::io::Write as _;
        let text = format!(
            "genome_id\tantibiotic\tsnap_mic_log2\tsnap_mic_censor\tmic_log2\tmic_censor\n\
             {AMP}.1\tampicillin\t5\tExact\t\t\n{VAN}.0\tvancomycin\t4\tExact\t\t\n{MERO}.0\tmeropenem\t3\tExact\t3\tExact\n{extra}"
        );
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        gz.write_all(text.as_bytes()).unwrap();
        gz.finish().unwrap()
    }

    #[cfg(unix)]
    fn fixture() -> Fixture {
        let dir = tempfile::tempdir().unwrap();
        let (research, reports) = (dir.path().join("research"), dir.path().join("reports"));
        let mut sums = String::from("00  data/curated/datasets.json\n");
        for (id, mics, predicted) in [
            (
                AMP,
                AMP_REPORT_MICS,
                Some([Some(1.5), Some(4.0), None, Some(-0.5)]),
            ),
            (MERO, MERO_MICS, Some(MERO_SAVED)),
            (CIP, CIP_MICS, None),
            (
                VAN,
                [Some(4), Some(2), Some(1), None],
                Some([Some(3.0), Some(2.0), Some(1.0), None]),
            ),
        ] {
            let bytes = serde_json::to_vec(&report_json(id, mics, predicted)).unwrap();
            sums += &format!("{}  data/report/datasets/{id}.json\n", sha256(&bytes));
            write(reports.join(format!("{id}.json")), bytes);
        }
        write(research.join(SNAPSHOT_SUMS), &sums);
        let snapshot = sha256(sums.as_bytes());
        let f = Fixture {
            _dir: dir,
            research,
            reports,
            snapshot,
        };
        for lane in [MIC_LANE, SELECTIVE_LANE] {
            write(
                f.research.join(lane).join("receipt.md"),
                format!("snapshot {}\n", f.snapshot),
            );
        }

        write(
            f.mic_lane(&format!("out/splits_postfix/{AMP}.json")),
            lane_json(AMP, AMP_LANE_MICS, AMP_RERUN).to_string(),
        );
        write(
            f.mic_lane(&format!("out/splits/{MERO}.json")),
            lane_json(MERO, MERO_MICS, MERO_SAVED).to_string(),
        );
        write(
            f.mic_lane(&format!("out/splits/{CIP}.json")),
            lane_json(CIP, CIP_MICS, CIP_RERUN).to_string(),
        );
        std::fs::create_dir_all(f.mic_lane("out/splits_final")).unwrap();
        for (id, dir) in [(AMP, "splits_postfix"), (MERO, "splits"), (CIP, "splits")] {
            std::os::unix::fs::symlink(
                format!("../{dir}/{id}.json"),
                f.mic_lane(&format!("out/splits_final/{id}.json")),
            )
            .unwrap();
        }
        let mut quality = String::from(
            "id,source,scheme,augmented,standard_applied,has_breakpoint,n_predicted,n_exact,ea_exact,ea_all,n_category,ref_R,ref_I,ref_S,ca,vme,vme_rate,vme_cp_lo,vme_cp_hi,me,me_rate,minor,minor_rate,fda_ea,fda_ca,fda_me,fda_vme,lineage_holdouts\n",
        );
        for id in [AMP, MERO, CIP] {
            for scheme in ["random", "lineage"] {
                quality += &quality_rows(id, "splits_final", scheme, 0.5);
                quality += &quality_rows(id, "splits_final/null_median", scheme, 0.4);
            }
        }
        // Reruns keep their saved_report rows, which the build does not read;
        // ampicillin's scored the zone diameter, so they differ from its rerun.
        quality += &quality_rows(AMP, "saved_report", "random", 0.25);
        quality += &quality_rows(MERO, "saved_report", "random", 0.5);
        quality += &quality_rows(VAN, "saved_report", "random", 0.5);
        write(f.mic_lane("out/analysis/prediction_quality.csv"), quality);
        write(f.mic_lane("out/pairs/pairs.tsv.gz"), pairs(""));

        let (pass_a1, pass_a2, refit_too, fail) = (
            result_cells(2, 0, 0, 2),
            result_cells(1, 0, 0, 2),
            result_cells(2, 0, 0, 1),
            result_cells(1, 1, 1, 1),
        );
        let prefixed = |p: &str| RESULT_FIELDS.map(|k| format!("{p}_{k}")).join(",");
        let empty = ",".repeat(RESULT_FIELDS.len() - 1);
        let label =
            |code: &str| format!("\"{}\"", CATEGORIES.iter().find(|c| c.0 == code).unwrap().1);
        // Ampicillin passes a refit too; its A1 result is the primary rule's.
        write(
            f.selective_lane("out/species_drug_summary.csv"),
            format!(
                "dataset,n,n_r,n_s,category,refit_best_learner,{},{}\n\
                 {AMP},4,2,2,{},gbt,{pass_a1},{refit_too}\n\
                 {MERO},4,2,2,{},gbt,{fail},{pass_a2}\n\
                 {VAN},4,2,2,{},,{fail},{empty}\n",
                prefixed("primary"),
                prefixed("refit_best"),
                label("A1"),
                label("A2"),
                label("B"),
            ),
        );
        f
    }

    #[cfg(unix)]
    #[test]
    fn builds_guarded_entries_from_lane_outputs() {
        let f = fixture();
        let r = f.build().unwrap();
        assert_eq!(r.snapshot.sha256, f.snapshot);
        assert_eq!(r.snapshot.datasets, 4);
        assert_eq!(r.sources.len(), 9);
        let mic = |id: &str| r.mic.datasets.iter().find(|d| d.id == id).unwrap();
        let mics = |m: [Option<i32>; 4]| {
            mics_digest(m.into_iter().map(|l| {
                l.map(|log2| Mic {
                    log2,
                    censor: Censor::Exact,
                })
            }))
        };

        // The zone-fix rerun is guarded by the MICs and predictions it scored.
        let amp = mic(AMP);
        assert_eq!(
            amp.zone_diameters,
            ZoneDiameters {
                removed: 1,
                scored: 0
            }
        );
        assert_eq!(amp.guard.mic_values, mics(AMP_LANE_MICS));
        assert_ne!(amp.guard.mic_values, mics(AMP_REPORT_MICS));
        assert_eq!(
            amp.guard.mic_predicted.as_deref(),
            Some(values_digest(&AMP_RERUN).as_str())
        );
        let random = &amp.model.random;
        assert_eq!(
            (random.n_predicted, random.n_exact, random.ea_exact),
            (3, 2, Some(0.5))
        );
        assert_eq!(
            amp.model.lineage.as_ref().unwrap().holdout_random_folds,
            Some(0)
        );
        assert_eq!(amp.baseline.as_ref().unwrap().lineage.ea_exact, Some(0.4));

        let mero = mic(MERO);
        assert_eq!(
            mero.guard.mic_predicted.as_deref(),
            Some(values_digest(&MERO_SAVED).as_str())
        );
        // A report without an MIC model gives the web guard no predictions.
        assert_eq!(mic(CIP).guard.mic_predicted, None);

        // The saved-only model scores a zone diameter, so only its reason
        // ships, with no revision-5 count in it.
        assert_eq!(r.mic.datasets.len(), 3);
        assert!(r.mic.datasets.iter().all(|d| d.id != VAN));
        assert_eq!(
            r.mic.withheld,
            [MicWithheld {
                id: VAN.to_string(),
                reason: "The revision-5 saved MIC model scored disk-diffusion zone diameters as MICs, and there is no rerun or no-genome baseline without them.".to_string(),
            }]
        );
        assert!(r.mic.limits.zone_fix.contains("{removed}"));

        let selective = |id: &str| r.selective.datasets.iter().find(|d| d.id == id).unwrap();
        let a1 = selective(AMP);
        assert_eq!(a1.guard.cohort, amp.guard.cohort);
        assert_eq!(a1.category, "A1");
        let a1_result = a1.result.as_ref().unwrap();
        assert_eq!(
            (
                a1_result.learner.as_deref(),
                a1_result.coverage,
                a1_result.tn
            ),
            (None, 1.0, 2)
        );
        let a2 = selective(MERO);
        assert_eq!(a2.category, "A2");
        let a2_result = a2.result.as_ref().unwrap();
        assert_eq!(a2_result.learner.as_deref(), Some("gbt"));
        assert_eq!((a2_result.coverage, a2_result.per100r_lab), (0.75, 50.0));
        let b = selective(VAN);
        assert_eq!((b.category, b.result.is_none()), ("B", true));
        let codes: Vec<&str> = r.selective.categories.iter().map(|c| c.code).collect();
        assert_eq!(codes, CATEGORIES.map(|c| c.0));
        // Only A2's label names a learner for the web app to write in.
        let named: Vec<&str> = r
            .selective
            .categories
            .iter()
            .filter(|c| c.label.contains("{learner}"))
            .map(|c| c.code)
            .collect();
        assert_eq!(named, ["A2"]);
        assert_eq!(
            r.selective.categories[6].label,
            "This report has no lineage split, so selective calling has no verdict for it."
        );
    }

    #[cfg(unix)]
    #[test]
    fn withholds_any_model_that_scores_zone_diameters() {
        let f = fixture();
        // Meropenem's rerun now scores an MIC that the zone fix removes.
        write(
            f.mic_lane("out/pairs/pairs.tsv.gz"),
            pairs(&format!("{MERO}.1\tmeropenem\t4\tExact\t\t\n")),
        );
        let r = f.build().unwrap();
        let ids: Vec<&str> = r.mic.datasets.iter().map(|d| d.id.as_str()).collect();
        assert_eq!(ids, [AMP, CIP]);
        assert_eq!(r.mic.withheld[0].id, VAN);
        assert_eq!(
            r.mic.withheld[1],
            MicWithheld {
                id: MERO.to_string(),
                reason: "The revision-5 rerun MIC model scored disk-diffusion zone diameters as MICs, and there is no rerun or no-genome baseline without them.".to_string(),
            }
        );
        assert_eq!(r.selective.datasets.len(), 3);
    }

    #[cfg(unix)]
    #[test]
    fn refuses_a_rerun_training_fold_without_an_exact_mic() {
        let f = fixture();
        // Meropenem's only lineage-split fold holds every exact MIC, so its
        // training fold has none.
        let mut lane = lane_json(MERO, MERO_MICS, MERO_SAVED);
        lane["runs"][1]["folds"] = json!([0, 0, 0, 0]);
        write(
            f.mic_lane(&format!("out/splits/{MERO}.json")),
            lane.to_string(),
        );
        assert!(f.error().contains(&format!(
            "splits_final/{MERO}.json: a lineage training fold has no exact MIC"
        )));
    }

    #[cfg(unix)]
    #[test]
    fn refuses_another_snapshot_and_unlisted_report_bytes() {
        let f = fixture();
        let other = build(&f.research, &f.reports, &"0".repeat(64))
            .err()
            .unwrap();
        assert!(format!("{other:#}").contains("pinned revision-5 snapshot"));
        let path = f.reports.join(format!("{MERO}.json"));
        let mut bytes = std::fs::read(&path).unwrap();
        bytes.push(b'\n');
        std::fs::write(&path, bytes).unwrap();
        assert!(f.error().contains("does not match the snapshot hash list"));
    }

    #[cfg(unix)]
    #[test]
    fn refuses_lane_rows_on_another_sample_order_or_genome() {
        let f = fixture();
        let path = f.mic_lane(&format!("out/splits/{MERO}.json"));
        let mut lane = lane_json(MERO, MERO_MICS, MERO_SAVED);
        lane["samples"].as_array_mut().unwrap().swap(0, 1);
        write(path.clone(), lane.to_string());
        assert!(f.error().contains("is not the report's sample"));
        let mut lane = lane_json(MERO, MERO_MICS, MERO_SAVED);
        lane["samples"][2]["id"] = json!("573.99999");
        write(path, lane.to_string());
        assert!(
            f.error()
                .contains("lane sample 2 (573.99999) is not the report's sample")
        );
    }

    #[cfg(unix)]
    #[test]
    fn refuses_mic_differences_the_zone_fix_does_not_explain() {
        let f = fixture();
        let mics = [Some(3), Some(4), Some(1), Some(-1)];
        write(
            f.mic_lane(&format!("out/splits/{MERO}.json")),
            lane_json(MERO, mics, MERO_SAVED).to_string(),
        );
        assert!(f.error().contains("which the zone fix does not explain"));
    }

    #[cfg(unix)]
    #[test]
    fn refuses_mic_rows_that_disagree_on_essential_agreement() {
        let f = fixture();
        let path = f.mic_lane("out/analysis/prediction_quality.csv");
        let quality = std::fs::read_to_string(&path).unwrap();
        // Meropenem's FDA row for the random split, with another value.
        let fda = format!("{MERO},splits_final,random,False,FDA,False,3,2,0.5,");
        assert_eq!(quality.matches(&fda).count(), 1);
        write(
            path,
            quality.replace(&fda, &fda.replace(",0.5,", ",0.625,")),
        );
        assert!(f.error().contains(&format!(
            "{MERO}: splits_final random rows disagree on essential agreement"
        )));
    }

    #[cfg(unix)]
    #[test]
    fn refuses_summary_rows_the_report_or_the_dashboard_would_reject() {
        let f = fixture();
        let path = f.selective_lane("out/species_drug_summary.csv");
        let summary = std::fs::read_to_string(&path).unwrap();
        let a2 = CATEGORIES[1].1;
        let mero = format!("{MERO},4,2,2,\"{a2}\",gbt,");
        assert_eq!(summary.matches(&mero).count(), 1);
        for (edited, message) in [
            (
                summary.replace(&mero, &format!("{MERO},4,1,3,\"{a2}\",gbt,")),
                "the summary counts 4 isolates (1 R, 3 S); the report 4 (2 R, 2 S)",
            ),
            (
                summary.replace(a2, "A3: a category the lane never wrote"),
                "unknown category",
            ),
            (
                summary.replace(a2, CATEGORIES[6].1),
                "category no_lineage, but the report has a lineage split",
            ),
            (
                summary.replace(&mero, &format!("{MERO},4,2,2,\"{a2}\",lasso,")),
                "an A2 result needs refit_best_learner gbt or rf, not \"lasso\"",
            ),
            // The refit calls 3 of meropenem's 2 resistant isolates.
            (
                summary.replace(",1,0,0,2\n", ",3,0,0,2\n"),
                "the refit_best_ counts do not partition the dataset",
            ),
        ] {
            write(path.clone(), edited);
            let error = f.error();
            assert!(error.contains(&format!("{MERO}: {message}")), "{error}");
        }
    }
}
