//! Provenance side files: where each sample's genome and AST records came
//! from, kept beside the pinned data without changing a byte of it.
//!
//! - `grm fetch genome-dates` writes `bvbrc/genome_dates.jsonl`, BV-BRC's
//!   `date_inserted` and `completion_date` for each genome in
//!   `species/*/genomes.txt` (`bvbrc::GenomeDate`), and
//!   `bvbrc/genome_dates.manifest.json`, which pins it: the genome lists it
//!   was asked for, how many came back, which did not, and the file's
//!   SHA-256. `genomes.jsonl` is not touched.
//! - `grm curate` writes `curated/provenance.jsonl`, one
//!   `curate::PairProvenance` per phenotype.
//! - `grm provenance` writes `provenance/<id>.json` beside each dataset
//!   report that `index.json` lists, as `grm intervals` does: the report's
//!   samples, in its order, each with its genome's dates and its pair's
//!   record details. Each file records the SHA-256 and evaluation revision
//!   of the report it describes and the SHA-256 of both inputs, and `grm
//!   provenance --check` fails when any of them no longer matches.
//!
//! Nothing here enters an input fingerprint, a feature block or a learner,
//! and no report file changes. The dates and years are source proxies that
//! `ablate::blocks` refuses by name; no block adapter reads these files, and
//! only its adapters can be admitted (docs/provenance-fields.md).

use crate::bvbrc::{GENOME_DATE_FIELDS, GenomeDate};
use crate::curate::PairProvenance;
use crate::layout::{Layout, ensure_parent};
use crate::research::sha256;
use anyhow::{Context, Result, bail, ensure};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;

/// Version of the genome-dates manifest and the provenance file layouts.
pub const SCHEMA: u32 = 1;

/// One species' `genomes.txt`: its genome ids and the file's SHA-256.
#[derive(Clone, Debug, PartialEq)]
pub struct GenomeList {
    pub slug: String,
    pub ids: Vec<String>,
    pub sha256: String,
}

/// The genome list a dates file was fetched for, as the manifest pins it.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct ListPin {
    pub genomes: usize,
    pub sha256: String,
}

/// `bvbrc/genome_dates.manifest.json`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DatesManifest {
    pub schema: u32,
    /// The BV-BRC collection and the fields selected from it.
    pub source: String,
    pub fields: String,
    pub fetched_at: String,
    /// Each species' genomes.txt, by slug.
    pub genome_lists: BTreeMap<String, ListPin>,
    /// The distinct genome ids over every list, and the SHA-256 of their
    /// sorted lines, each ended by "\n".
    pub requested: usize,
    pub requested_sha256: String,
    /// Rows written, one per genome BV-BRC returned.
    pub returned: usize,
    /// Requested genomes BV-BRC returned no record for.
    pub missing: Vec<String>,
    pub without_date_inserted: usize,
    pub without_completion_date: usize,
    /// The dates file, beside the manifest.
    pub file: String,
    pub bytes: u64,
    pub sha256: String,
}

/// Every `species/<slug>/genomes.txt` under the data root, by slug.
pub fn genome_lists(layout: &Layout) -> Result<Vec<GenomeList>> {
    let dir = layout.root.join("species");
    let mut slugs = Vec::new();
    for entry in std::fs::read_dir(&dir)
        .with_context(|| format!("run `grm curate` first ({})", dir.display()))?
    {
        let entry = entry?;
        if entry.path().join("genomes.txt").is_file() {
            let name = entry.file_name();
            let slug = name
                .to_str()
                .with_context(|| format!("{} is not UTF-8", entry.path().display()))?;
            slugs.push(slug.to_string());
        }
    }
    slugs.sort();
    ensure!(
        !slugs.is_empty(),
        "no species/*/genomes.txt under {}; run `grm curate` first",
        layout.root.display()
    );
    slugs
        .into_iter()
        .map(|slug| {
            let path = layout.species_dir(&slug).join("genomes.txt");
            let bytes = std::fs::read(&path).with_context(|| format!("read {}", path.display()))?;
            let ids = std::str::from_utf8(&bytes)
                .with_context(|| format!("{} is not UTF-8", path.display()))?
                .lines()
                .map(str::trim)
                .filter(|l| !l.is_empty())
                .map(String::from)
                .collect();
            Ok(GenomeList {
                sha256: sha256(&bytes),
                slug,
                ids,
            })
        })
        .collect()
}

/// Writes `rows`, the dates BV-BRC returned for `ids` (sorted and distinct,
/// over `lists`), as `bvbrc/genome_dates.jsonl` in genome id order, then the
/// manifest that pins it, and returns the manifest. A row for a genome not
/// asked for, or a second row for one, is an error, and nothing is written.
pub fn write_genome_dates(
    layout: &Layout,
    lists: &[GenomeList],
    ids: &[String],
    mut rows: Vec<GenomeDate>,
) -> Result<DatesManifest> {
    ensure!(
        ids.windows(2).all(|w| w[0] < w[1]),
        "the requested genome ids are not sorted and distinct"
    );
    let wanted: BTreeSet<&str> = ids.iter().map(String::as_str).collect();
    rows.sort_by(|a, b| a.genome_id.cmp(&b.genome_id));
    for (i, row) in rows.iter().enumerate() {
        ensure!(
            wanted.contains(row.genome_id.as_str()),
            "BV-BRC returned genome {}, which was not asked for",
            row.genome_id
        );
        ensure!(
            i == 0 || rows[i - 1].genome_id != row.genome_id,
            "BV-BRC returned genome {} twice",
            row.genome_id
        );
    }
    let returned: BTreeSet<&str> = rows.iter().map(|r| r.genome_id.as_str()).collect();
    let mut bytes = Vec::new();
    for row in &rows {
        serde_json::to_writer(&mut bytes, row)?;
        bytes.push(b'\n');
    }
    let requested_lines: String = ids.iter().map(|id| format!("{id}\n")).collect();
    let path = layout.genome_dates();
    let manifest = DatesManifest {
        schema: SCHEMA,
        source: "BV-BRC genome".into(),
        fields: GENOME_DATE_FIELDS.into(),
        fetched_at: crate::report::now_iso(),
        genome_lists: lists
            .iter()
            .map(|l| {
                let pin = ListPin {
                    genomes: l.ids.len(),
                    sha256: l.sha256.clone(),
                };
                (l.slug.clone(), pin)
            })
            .collect(),
        requested: ids.len(),
        requested_sha256: sha256(requested_lines.as_bytes()),
        returned: rows.len(),
        missing: ids
            .iter()
            .filter(|id| !returned.contains(id.as_str()))
            .cloned()
            .collect(),
        without_date_inserted: rows.iter().filter(|r| r.date_inserted.is_none()).count(),
        without_completion_date: rows.iter().filter(|r| r.completion_date.is_none()).count(),
        file: file_name(&path)?,
        bytes: bytes.len() as u64,
        sha256: sha256(&bytes),
    };
    ensure_parent(&path)?;
    crate::report::publish(&path, &bytes)?;
    crate::report::publish(
        &layout.genome_dates_manifest(),
        &serde_json::to_vec_pretty(&manifest)?,
    )?;
    Ok(manifest)
}

fn file_name(path: &Path) -> Result<String> {
    Ok(path
        .file_name()
        .and_then(|n| n.to_str())
        .with_context(|| format!("{} has no file name", path.display()))?
        .to_string())
}

/// Splits JSON lines into records.
fn parse_lines<T: serde::de::DeserializeOwned>(bytes: &[u8], what: &Path) -> Result<Vec<T>> {
    std::str::from_utf8(bytes)
        .with_context(|| format!("{} is not UTF-8", what.display()))?
        .lines()
        .enumerate()
        .filter(|(_, line)| !line.is_empty())
        .map(|(i, line)| {
            serde_json::from_str(line).with_context(|| format!("{}:{}", what.display(), i + 1))
        })
        .collect()
}

/// The two inputs of every provenance file, each with the SHA-256 it is
/// pinned by.
pub struct Inputs {
    pub dates: HashMap<String, GenomeDate>,
    pub dates_sha256: String,
    pub pairs: HashMap<(String, String), PairProvenance>,
    pub pairs_sha256: String,
}

impl Inputs {
    /// Reads the genome dates, refusing a file its manifest does not pin,
    /// and curated/provenance.jsonl.
    pub fn load(layout: &Layout) -> Result<Inputs> {
        let manifest_path = layout.genome_dates_manifest();
        let manifest: DatesManifest =
            serde_json::from_slice(&std::fs::read(&manifest_path).with_context(|| {
                format!(
                    "run `grm fetch genome-dates` first ({})",
                    manifest_path.display()
                )
            })?)
            .with_context(|| format!("parse {}", manifest_path.display()))?;
        ensure!(
            manifest.schema == SCHEMA,
            "{} has schema {}, not {SCHEMA}",
            manifest_path.display(),
            manifest.schema
        );
        let path = layout.genome_dates();
        let bytes = std::fs::read(&path).with_context(|| format!("read {}", path.display()))?;
        ensure!(
            manifest.file == file_name(&path)?
                && manifest.bytes == bytes.len() as u64
                && manifest.sha256 == sha256(&bytes),
            "{} does not match {}: rerun `grm fetch genome-dates`",
            path.display(),
            manifest_path.display()
        );
        let mut dates = HashMap::new();
        for row in parse_lines::<GenomeDate>(&bytes, &path)? {
            if let Some(old) = dates.insert(row.genome_id.clone(), row) {
                bail!("{} lists genome {} twice", path.display(), old.genome_id);
            }
        }
        let path = layout.curated().join("provenance.jsonl");
        let bytes = std::fs::read(&path)
            .with_context(|| format!("run `grm curate` first, which writes {}", path.display()))?;
        let mut pairs = HashMap::new();
        for row in parse_lines::<PairProvenance>(&bytes, &path)? {
            let key = (row.genome_id.clone(), row.antibiotic.clone());
            if pairs.insert(key, row).is_some() {
                bail!("{} lists a genome and antibiotic twice", path.display());
            }
        }
        Ok(Inputs {
            dates,
            dates_sha256: manifest.sha256,
            pairs,
            pairs_sha256: sha256(&bytes),
        })
    }
}

/// The provenance file for one dataset report: per sample, in the report's
/// order, its genome's dates (null when the dates file has no row for it)
/// and its pair's record details (null when curation has none), pinned to
/// the report's SHA-256 and evaluation revision and to both inputs.
pub fn for_report(bytes: &[u8], inputs: &Inputs) -> Result<Value> {
    let report: Value = serde_json::from_slice(bytes).context("parse a dataset report")?;
    let id = report["id"]
        .as_str()
        .context("a dataset report has no id")?;
    let antibiotic = report["antibiotic"]
        .as_str()
        .with_context(|| format!("{id}: the report has no antibiotic"))?;
    let (mut undated, mut unpaired) = (0, 0);
    let mut samples = Vec::new();
    for sample in report["samples"].as_array().into_iter().flatten() {
        let genome = sample["id"]
            .as_str()
            .with_context(|| format!("{id}: a sample has no id"))?;
        let dates = inputs.dates.get(genome).map(
            |d| json!({ "date_inserted": d.date_inserted, "completion_date": d.completion_date }),
        );
        let pair = match inputs
            .pairs
            .get(&(genome.to_string(), antibiotic.to_string()))
        {
            Some(p) => {
                let mut v = serde_json::to_value(p)?;
                let fields = v.as_object_mut().context("a pair row is not an object")?;
                fields.remove("genome_id");
                fields.remove("antibiotic");
                Some(v)
            }
            None => None,
        };
        undated += usize::from(dates.is_none());
        unpaired += usize::from(pair.is_none());
        samples.push(json!({ "genome_id": genome, "genome_dates": dates, "pair": pair }));
    }
    Ok(json!({
        "schema": SCHEMA,
        "dataset": id,
        "antibiotic": antibiotic,
        "report_sha256": sha256(bytes),
        "evaluation_revision": report["evaluation"]["revision"],
        "genome_dates_sha256": inputs.dates_sha256,
        "pairs_sha256": inputs.pairs_sha256,
        "counts": {
            "samples": samples.len(),
            "without_genome_dates": undated,
            "without_pair": unpaired,
        },
        "samples": samples,
    }))
}

/// The dataset ids `index.json` lists, under the report directory's lock,
/// which the caller holds while it reads the reports, as `grm report` does.
fn locked_index(report_dir: &Path) -> Result<(std::fs::File, Vec<String>)> {
    let lock = std::fs::File::open(report_dir)
        .with_context(|| format!("open {}", report_dir.display()))?;
    lock.lock()
        .with_context(|| format!("lock {}", report_dir.display()))?;
    let index: Value = serde_json::from_slice(
        &std::fs::read(report_dir.join("index.json")).context("run `grm report` first")?,
    )?;
    let ids = index["datasets"]
        .as_array()
        .into_iter()
        .flatten()
        .map(|row| {
            row["id"]
                .as_str()
                .map(String::from)
                .context("index.json lists a dataset with no id")
        })
        .collect::<Result<_>>()?;
    Ok((lock, ids))
}

fn read_report(report_dir: &Path, id: &str) -> Result<Vec<u8>> {
    let path = report_dir.join("datasets").join(format!("{id}.json"));
    std::fs::read(&path).with_context(|| format!("read {}", path.display()))
}

/// Writes `report/provenance/<id>.json` for every dataset `index.json`
/// lists, removes any other file there, and returns how many it wrote.
pub fn write(layout: &Layout) -> Result<usize> {
    let report_dir = layout.root.join("report");
    let (_lock, ids) = locked_index(&report_dir)?;
    let inputs = Inputs::load(layout)?;
    let out = report_dir.join("provenance");
    std::fs::create_dir_all(&out)?;
    let mut kept = BTreeSet::new();
    for id in &ids {
        let file = for_report(&read_report(&report_dir, id)?, &inputs)?;
        if file["dataset"] != id.as_str() {
            bail!("the report for {id} records the id {}", file["dataset"]);
        }
        let name = format!("{id}.json");
        crate::report::publish(&out.join(&name), &serde_json::to_vec(&file)?)?;
        kept.insert(name);
    }
    for entry in std::fs::read_dir(&out)? {
        let entry = entry?;
        if !entry.file_name().to_str().is_some_and(|n| kept.contains(n)) {
            std::fs::remove_file(entry.path())
                .with_context(|| format!("remove {}", entry.path().display()))?;
        }
    }
    Ok(kept.len())
}

/// Checks every provenance file against the reports and inputs as they are
/// now: one per indexed dataset and no other, each with its report naming
/// that dataset, each pinned to its report's
/// current SHA-256 and evaluation revision and to the inputs' SHA-256, with
/// the content those give. Returns how many it checked, or fails naming
/// each that does not hold.
pub fn check(layout: &Layout) -> Result<usize> {
    let report_dir = layout.root.join("report");
    let (_lock, ids) = locked_index(&report_dir)?;
    let inputs = Inputs::load(layout)?;
    let dir = report_dir.join("provenance");
    let mut problems = Vec::new();
    for id in &ids {
        let report = read_report(&report_dir, id)?;
        let path = dir.join(format!("{id}.json"));
        let Ok(bytes) = std::fs::read(&path) else {
            problems.push(format!("{id}: no {}", path.display()));
            continue;
        };
        let file: Value =
            serde_json::from_slice(&bytes).with_context(|| format!("parse {}", path.display()))?;
        let expected = for_report(&report, &inputs)?;
        let before = problems.len();
        // As `write` refuses: the report and its file belong to this id, not
        // another dataset's copied under its name.
        for (what, dataset) in [("report", &expected["dataset"]), ("file", &file["dataset"])] {
            if dataset != id.as_str() {
                problems.push(format!("{id}: the {what} records the dataset {dataset}"));
            }
        }
        for pin in [
            "report_sha256",
            "evaluation_revision",
            "genome_dates_sha256",
            "pairs_sha256",
        ] {
            if file[pin] != expected[pin] {
                problems.push(format!(
                    "{id}: {pin} is {}, now {}",
                    file[pin], expected[pin]
                ));
            }
        }
        // A file whose pins all hold can still have been edited in place.
        if problems.len() == before && file != expected {
            problems.push(format!("{id}: the file differs from what its inputs give"));
        }
    }
    let indexed: BTreeSet<String> = ids.iter().map(|id| format!("{id}.json")).collect();
    if let Ok(entries) = std::fs::read_dir(&dir) {
        for entry in entries {
            let name = entry?.file_name().to_string_lossy().into_owned();
            if !indexed.contains(&name) {
                problems.push(format!("{name}: not a dataset index.json lists"));
            }
        }
    }
    if !problems.is_empty() {
        bail!(
            "{} provenance problems:\n  {}",
            problems.len(),
            problems.join("\n  ")
        );
    }
    Ok(ids.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn date(id: &str, inserted: Option<&str>) -> GenomeDate {
        GenomeDate {
            genome_id: id.into(),
            date_inserted: inserted.map(Into::into),
            completion_date: Some("2014-07-29T00:00:00Z".into()),
        }
    }

    /// A data root with species lists a (g1, g2) and b (g3), one report over
    /// g1, g2 and g3 for meropenem, and a curated pair for g1 and g2.
    fn root() -> (tempfile::TempDir, Layout) {
        let temp = tempfile::tempdir().unwrap();
        let layout = Layout::new(temp.path());
        for (slug, ids) in [("a", "g2\ng1\n"), ("b", "g3\n")] {
            std::fs::create_dir_all(layout.species_dir(slug)).unwrap();
            std::fs::write(layout.species_dir(slug).join("genomes.txt"), ids).unwrap();
        }
        let report = layout.root.join("report");
        std::fs::create_dir_all(report.join("datasets")).unwrap();
        std::fs::write(
            report.join("index.json"),
            br#"{"datasets":[{"id":"sp__meropenem"}]}"#,
        )
        .unwrap();
        std::fs::write(report.join("datasets/sp__meropenem.json"), report_json(12)).unwrap();
        std::fs::create_dir_all(layout.curated()).unwrap();
        let pair = |g: &str| PairProvenance {
            genome_id: g.into(),
            antibiotic: "meropenem".into(),
            methods: vec!["Broth dilution".into()],
            n_records: 2,
            date_inserted_first: Some("2015-01-01T00:00:00Z".into()),
            pmids: vec!["26686880".into()],
            ..Default::default()
        };
        let rows: String = ["g1", "g2"]
            .iter()
            .map(|g| serde_json::to_string(&pair(g)).unwrap() + "\n")
            .collect();
        std::fs::write(layout.curated().join("provenance.jsonl"), rows).unwrap();
        (temp, layout)
    }

    fn report_json(revision: u32) -> Vec<u8> {
        serde_json::to_vec(&json!({
            "id": "sp__meropenem",
            "antibiotic": "meropenem",
            "evaluation": { "revision": revision },
            "samples": [{ "id": "g3" }, { "id": "g1" }, { "id": "g2" }],
        }))
        .unwrap()
    }

    fn fetch(layout: &Layout, rows: Vec<GenomeDate>) -> Result<DatesManifest> {
        let lists = genome_lists(layout)?;
        let ids: Vec<String> = lists
            .iter()
            .flat_map(|l| l.ids.iter().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        write_genome_dates(layout, &lists, &ids, rows)
    }

    #[test]
    fn genome_dates_are_sorted_and_pinned_with_what_was_missing() {
        let (_temp, layout) = root();
        let manifest = fetch(
            &layout,
            vec![
                date("g2", None),
                date("g1", Some("2014-12-08T22:10:02.394Z")),
            ],
        )
        .unwrap();
        let bytes = std::fs::read(layout.genome_dates()).unwrap();
        assert_eq!(
            String::from_utf8(bytes.clone()).unwrap(),
            "{\"genome_id\":\"g1\",\"date_inserted\":\"2014-12-08T22:10:02.394Z\",\"completion_date\":\"2014-07-29T00:00:00Z\"}\n\
             {\"genome_id\":\"g2\",\"date_inserted\":null,\"completion_date\":\"2014-07-29T00:00:00Z\"}\n"
        );
        assert_eq!(manifest.sha256, sha256(&bytes));
        assert_eq!(manifest.bytes, bytes.len() as u64);
        assert_eq!((manifest.requested, manifest.returned), (3, 2));
        assert_eq!(manifest.missing, ["g3"]);
        assert_eq!(manifest.without_date_inserted, 1);
        assert_eq!(manifest.without_completion_date, 0);
        assert_eq!(manifest.requested_sha256, sha256(b"g1\ng2\ng3\n"));
        assert_eq!(manifest.genome_lists["a"].genomes, 2);
        assert_eq!(manifest.genome_lists["b"].sha256, sha256(b"g3\n"));
        assert_eq!(manifest.fields, "genome_id,date_inserted,completion_date");
        let on_disk: DatesManifest =
            serde_json::from_slice(&std::fs::read(layout.genome_dates_manifest()).unwrap())
                .unwrap();
        assert_eq!(on_disk, manifest);
    }

    #[test]
    fn genome_dates_refuse_rows_not_asked_for_and_repeats() {
        let (_temp, layout) = root();
        let err = fetch(&layout, vec![date("g9", None)]).unwrap_err();
        assert!(format!("{err:#}").contains("not asked for"), "{err:#}");
        let err = fetch(&layout, vec![date("g1", None), date("g1", None)]).unwrap_err();
        assert!(format!("{err:#}").contains("twice"), "{err:#}");
        assert!(!layout.genome_dates().exists());
    }

    #[test]
    fn provenance_files_pin_the_report_and_leave_it_unchanged() {
        let (_temp, layout) = root();
        fetch(&layout, vec![date("g1", Some("2014-12-08T22:10:02.394Z"))]).unwrap();
        let report_path = layout.root.join("report/datasets/sp__meropenem.json");
        let before = std::fs::read(&report_path).unwrap();
        std::fs::create_dir_all(layout.root.join("report/provenance")).unwrap();
        std::fs::write(layout.root.join("report/provenance/gone__x.json"), b"{}").unwrap();
        assert_eq!(write(&layout).unwrap(), 1);
        assert_eq!(std::fs::read(&report_path).unwrap(), before);
        assert!(!layout.root.join("report/provenance/gone__x.json").exists());
        let file: Value = serde_json::from_slice(
            &std::fs::read(layout.root.join("report/provenance/sp__meropenem.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(file["report_sha256"], sha256(&before));
        assert_eq!(file["evaluation_revision"], 12);
        assert_eq!(
            file["genome_dates_sha256"],
            sha256(&std::fs::read(layout.genome_dates()).unwrap())
        );
        assert_eq!(
            file["pairs_sha256"],
            sha256(&std::fs::read(layout.curated().join("provenance.jsonl")).unwrap())
        );
        assert_eq!(
            file["counts"],
            json!({ "samples": 3, "without_genome_dates": 2, "without_pair": 1 })
        );
        // The report's order: g3 has neither, g1 both, g2 a pair only.
        let samples = file["samples"].as_array().unwrap();
        assert_eq!(
            samples[0],
            json!({ "genome_id": "g3", "genome_dates": null, "pair": null })
        );
        assert_eq!(
            samples[1]["genome_dates"]["date_inserted"],
            "2014-12-08T22:10:02.394Z"
        );
        assert_eq!(samples[1]["pair"]["pmids"], json!(["26686880"]));
        assert_eq!(samples[1]["pair"]["n_records"], 2);
        assert!(samples[1]["pair"].get("genome_id").is_none());
        assert!(samples[2]["genome_dates"].is_null());
        assert_eq!(check(&layout).unwrap(), 1);
    }

    #[test]
    fn check_fails_when_a_report_or_input_changes_after_the_pin() {
        let (_temp, layout) = root();
        fetch(&layout, vec![date("g1", None)]).unwrap();
        write(&layout).unwrap();
        let report_path = layout.root.join("report/datasets/sp__meropenem.json");
        let original = std::fs::read(&report_path).unwrap();

        // An edited report, same revision.
        let mut edited = original.clone();
        edited.extend(b" ");
        std::fs::write(&report_path, &edited).unwrap();
        let err = format!("{:#}", check(&layout).unwrap_err());
        assert!(err.contains("sp__meropenem: report_sha256"), "{err}");
        // A report of another evaluation revision.
        std::fs::write(&report_path, report_json(13)).unwrap();
        let err = format!("{:#}", check(&layout).unwrap_err());
        assert!(err.contains("evaluation_revision is 12, now 13"), "{err}");
        std::fs::write(&report_path, &original).unwrap();
        check(&layout).unwrap();

        // A provenance file edited in place keeps its pins but not its content.
        let path = layout.root.join("report/provenance/sp__meropenem.json");
        let mut file: Value = serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
        file["samples"][1]["pair"]["n_records"] = json!(9);
        std::fs::write(&path, serde_json::to_vec(&file).unwrap()).unwrap();
        let err = format!("{:#}", check(&layout).unwrap_err());
        assert!(err.contains("differs from what its inputs give"), "{err}");
        write(&layout).unwrap();
        check(&layout).unwrap();

        // New curated rows change the pin.
        let pairs = layout.curated().join("provenance.jsonl");
        let mut rows = std::fs::read(&pairs).unwrap();
        rows.extend(b"\n");
        std::fs::write(&pairs, rows).unwrap();
        let err = format!("{:#}", check(&layout).unwrap_err());
        assert!(err.contains("pairs_sha256"), "{err}");
        write(&layout).unwrap();

        // A dates file its manifest does not pin is refused.
        let mut dates = std::fs::read(layout.genome_dates()).unwrap();
        dates.extend(b"\n");
        std::fs::write(layout.genome_dates(), dates).unwrap();
        let err = format!("{:#}", check(&layout).unwrap_err());
        assert!(err.contains("does not match"), "{err}");
        let err = format!("{:#}", write(&layout).unwrap_err());
        assert!(err.contains("does not match"), "{err}");

        // A stray file, and a missing one.
        fetch(&layout, vec![date("g1", None)]).unwrap();
        write(&layout).unwrap();
        std::fs::write(layout.root.join("report/provenance/other__x.json"), b"{}").unwrap();
        std::fs::remove_file(&path).unwrap();
        let err = format!("{:#}", check(&layout).unwrap_err());
        assert!(err.contains("2 provenance problems"), "{err}");
        assert!(
            err.contains("sp__meropenem: no ") && err.contains("other__x.json: not"),
            "{err}"
        );
    }

    #[test]
    fn check_fails_when_one_dataset_s_report_and_file_sit_under_another_id() {
        // REVIEW-PROV L1: copy X's report and side file over Y's. Each still
        // pins its own bytes, so only the id comparison catches it.
        let (_temp, layout) = root();
        let report = layout.root.join("report");
        let mut other: Value = serde_json::from_slice(&report_json(12)).unwrap();
        other["id"] = json!("sp__other");
        std::fs::write(
            report.join("datasets/sp__other.json"),
            serde_json::to_vec(&other).unwrap(),
        )
        .unwrap();
        std::fs::write(
            report.join("index.json"),
            br#"{"datasets":[{"id":"sp__meropenem"},{"id":"sp__other"}]}"#,
        )
        .unwrap();
        fetch(&layout, vec![date("g1", None)]).unwrap();
        assert_eq!(write(&layout).unwrap(), 2);
        assert_eq!(check(&layout).unwrap(), 2);
        for dir in ["datasets", "provenance"] {
            std::fs::copy(
                report.join(format!("{dir}/sp__meropenem.json")),
                report.join(format!("{dir}/sp__other.json")),
            )
            .unwrap();
        }
        let err = format!("{:#}", check(&layout).unwrap_err());
        assert!(
            err.contains("sp__other: the report records the dataset \"sp__meropenem\""),
            "{err}"
        );
        assert!(
            err.contains("sp__other: the file records the dataset \"sp__meropenem\""),
            "{err}"
        );
        // `write` refuses the same tree.
        let err = format!("{:#}", write(&layout).unwrap_err());
        assert!(err.contains("records the id"), "{err}");
    }

    #[test]
    fn provenance_needs_both_inputs() {
        let (_temp, layout) = root();
        let err = format!("{:#}", write(&layout).unwrap_err());
        assert!(err.contains("grm fetch genome-dates"), "{err}");
        fetch(&layout, vec![]).unwrap();
        std::fs::remove_file(layout.curated().join("provenance.jsonl")).unwrap();
        let err = format!("{:#}", write(&layout).unwrap_err());
        assert!(err.contains("grm curate"), "{err}");
    }
}
