//! BV-BRC data access: laboratory AMR records, genome metadata, assemblies,
//! and gene annotations, all through the public REST API.
//!
//! The ftp.bvbrc.org HTTPS mirror times out from this network, so assemblies
//! come from the `genome_sequence` collection as FASTA instead.

use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::write::GzEncoder;
use serde::de::DeserializeOwned;
use serde::{Deserialize, Deserializer, Serialize};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Duration;

const API: &str = "https://www.bv-brc.org/api";
const PAGE: usize = 25_000;

/// One laboratory AST record, as stored in `genome_amr`.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct AmrRecord {
    #[serde(deserialize_with = "text")]
    pub genome_id: String,
    #[serde(default, deserialize_with = "opt_text")]
    pub genome_name: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub taxon_id: Option<String>,
    #[serde(deserialize_with = "text")]
    pub antibiotic: String,
    #[serde(default, deserialize_with = "opt_text")]
    pub resistant_phenotype: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub measurement: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub measurement_sign: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub measurement_value: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub measurement_unit: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub laboratory_typing_method: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub laboratory_typing_platform: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub vendor: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub testing_standard: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub testing_standard_year: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub source: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub pmid: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub date_inserted: Option<String>,
}

const AMR_FIELDS: &str = "genome_id,genome_name,taxon_id,antibiotic,resistant_phenotype,measurement,measurement_sign,\
measurement_value,measurement_unit,laboratory_typing_method,laboratory_typing_platform,vendor,testing_standard,\
testing_standard_year,source,pmid,date_inserted";

/// Genome-level metadata used for curation, QC, lineage, and context views.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct GenomeRecord {
    #[serde(deserialize_with = "text")]
    pub genome_id: String,
    #[serde(default, deserialize_with = "opt_text")]
    pub genome_name: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub species: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub genus: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub taxon_id: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub genome_status: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub genome_quality: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub genome_length: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub contigs: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub contig_n50: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub gc_content: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub checkm_completeness: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub checkm_contamination: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub isolation_country: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub geographic_group: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub collection_year: Option<String>,
    /// BV-BRC's `collection_date`, read strictly (`CollectionDate`). A
    /// `genomes.jsonl` line written before `GENOME_FIELDS` asked for it has
    /// no such key: it reads as `NotFetched` and is written back without the
    /// key, so the line and everything built from it stay as they were.
    #[serde(
        default,
        deserialize_with = "collection_date",
        skip_serializing_if = "CollectionDate::not_fetched"
    )]
    pub collection_date: CollectionDate,
    #[serde(default, deserialize_with = "opt_text")]
    pub host_name: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub isolation_source: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub mlst: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub cgmlst_hc5: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub cgmlst_hc20: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub cgmlst_hc100: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub sequencing_platform: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub assembly_method: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub biosample_accession: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub assembly_accession: Option<String>,
}

const GENOME_FIELDS: &str = "genome_id,genome_name,species,genus,taxon_id,genome_status,genome_quality,genome_length,\
contigs,contig_n50,gc_content,checkm_completeness,checkm_contamination,isolation_country,geographic_group,\
collection_year,collection_date,host_name,isolation_source,mlst,cgmlst_hc5,cgmlst_hc20,cgmlst_hc100,\
sequencing_platform,assembly_method,biosample_accession,assembly_accession";

/// When an isolate was collected, as BV-BRC's `collection_date` states it.
/// It is a genome field: it says when the sample was taken, not when a lab
/// result was recorded, reported or entered (`AmrRecord::date_inserted`,
/// `curate::PairProvenance`), nor when BV-BRC stored the genome
/// (`GenomeDate`). No lab record and no curated file carries it.
///
/// The text is read strictly. Only `YYYY-MM-DD` (a real calendar day),
/// `YYYY-MM` (month 01 to 12) and `YYYY` are dates, at day, month and year
/// precision. Any other text, such as a range (`2015-2016`), a truncated
/// range (`2004-20`), a time stamp, a slash or month-name form, or a day
/// that does not exist (`2015-02-30`), is `Unparsed`: kept as written and
/// never turned into a date. The year is not checked against any range;
/// whether a stated year is plausible is for the analysis that reads it.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum CollectionDate {
    /// The record comes from a file fetched before `GENOME_FIELDS` asked for
    /// `collection_date`, so nothing is known about it. Not the same as
    /// `Missing`.
    #[default]
    NotFetched,
    /// Asked for, and BV-BRC had no value.
    Missing,
    /// A day, month or year, at the precision the text states.
    Stated(StatedDate),
    /// Text that is none of the three forms, as BV-BRC wrote it (trimmed).
    Unparsed(String),
}

/// A collection date at the precision its text states. Nothing finer is
/// filled in: a month has no day and a year has no month.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StatedDate {
    Day { year: u16, month: u8, day: u8 },
    Month { year: u16, month: u8 },
    Year { year: u16 },
}

/// How finely a collection date is stated.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Precision {
    Day,
    Month,
    Year,
}

impl Precision {
    pub fn as_str(self) -> &'static str {
        match self {
            Precision::Day => "day",
            Precision::Month => "month",
            Precision::Year => "year",
        }
    }
}

impl StatedDate {
    pub fn precision(self) -> Precision {
        match self {
            StatedDate::Day { .. } => Precision::Day,
            StatedDate::Month { .. } => Precision::Month,
            StatedDate::Year { .. } => Precision::Year,
        }
    }

    pub fn year(self) -> u16 {
        match self {
            StatedDate::Day { year, .. }
            | StatedDate::Month { year, .. }
            | StatedDate::Year { year } => year,
        }
    }
}

/// The date in the one form it was read from, so a stated date is written
/// back exactly as BV-BRC wrote it.
impl std::fmt::Display for StatedDate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            StatedDate::Day { year, month, day } => write!(f, "{year:04}-{month:02}-{day:02}"),
            StatedDate::Month { year, month } => write!(f, "{year:04}-{month:02}"),
            StatedDate::Year { year } => write!(f, "{year:04}"),
        }
    }
}

impl CollectionDate {
    /// Reads one `collection_date` text, which fetch has already trimmed.
    pub fn parse(text: &str) -> Self {
        match stated_date(text) {
            Some(date) => CollectionDate::Stated(date),
            None => CollectionDate::Unparsed(text.to_string()),
        }
    }

    pub fn stated(&self) -> Option<StatedDate> {
        match self {
            CollectionDate::Stated(date) => Some(*date),
            _ => None,
        }
    }

    /// Day, month or year; None when there is no date or it did not parse.
    pub fn precision(&self) -> Option<Precision> {
        self.stated().map(StatedDate::precision)
    }

    /// The text BV-BRC wrote, whether or not it parsed.
    pub fn text(&self) -> Option<String> {
        match self {
            CollectionDate::Stated(date) => Some(date.to_string()),
            CollectionDate::Unparsed(text) => Some(text.clone()),
            CollectionDate::NotFetched | CollectionDate::Missing => None,
        }
    }

    pub fn is_unparsed(&self) -> bool {
        matches!(self, CollectionDate::Unparsed(_))
    }

    pub fn not_fetched(&self) -> bool {
        matches!(self, CollectionDate::NotFetched)
    }
}

/// `YYYY`, `YYYY-MM` or `YYYY-MM-DD` with ASCII digits only, a month from 01
/// to 12 and a day that exists in that month; anything else is None.
fn stated_date(text: &str) -> Option<StatedDate> {
    let b = text.as_bytes();
    let number = |range: std::ops::Range<usize>| -> Option<u16> {
        let digits = b.get(range)?;
        digits
            .iter()
            .all(u8::is_ascii_digit)
            .then(|| digits.iter().fold(0, |n, d| n * 10 + u16::from(d - b'0')))
    };
    let dash = |i: usize| b.get(i) == Some(&b'-');
    let year = number(0..4)?;
    if b.len() == 4 {
        return Some(StatedDate::Year { year });
    }
    if !dash(4) {
        return None;
    }
    let month = u8::try_from(number(5..7)?).ok()?;
    if !(1..=12).contains(&month) {
        return None;
    }
    if b.len() == 7 {
        return Some(StatedDate::Month { year, month });
    }
    if b.len() != 10 || !dash(7) {
        return None;
    }
    let day = u8::try_from(number(8..10)?).ok()?;
    (1..=days_in_month(year, month))
        .contains(&day)
        .then_some(StatedDate::Day { year, month, day })
}

fn days_in_month(year: u16, month: u8) -> u8 {
    match month {
        4 | 6 | 9 | 11 => 30,
        2 if year.is_multiple_of(4) && (!year.is_multiple_of(100) || year.is_multiple_of(400)) => {
            29
        }
        2 => 28,
        _ => 31,
    }
}

/// `CollectionDate` is stored as its text: a stated date in the form it was
/// read from, an unparsed one as written, and no value as null.
impl Serialize for CollectionDate {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        match self {
            CollectionDate::Stated(date) => s.collect_str(date),
            CollectionDate::Unparsed(text) => s.serialize_str(text),
            CollectionDate::NotFetched | CollectionDate::Missing => s.serialize_none(),
        }
    }
}

/// How many genomes of a fetch had each kind of collection date.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CollectionDateCounts {
    pub day: usize,
    pub month: usize,
    pub year: usize,
    pub unparsed: usize,
    pub missing: usize,
}

impl CollectionDateCounts {
    pub fn add(&mut self, date: &CollectionDate) {
        match date.precision() {
            Some(Precision::Day) => self.day += 1,
            Some(Precision::Month) => self.month += 1,
            Some(Precision::Year) => self.year += 1,
            None if date.is_unparsed() => self.unparsed += 1,
            None => self.missing += 1,
        }
    }
}

impl std::fmt::Display for CollectionDateCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "day {}, month {}, year {}, unparsed and kept as text {}, missing {}",
            self.day, self.month, self.year, self.unparsed, self.missing
        )
    }
}

/// Writes one page of BV-BRC's `genome` response as `genomes.jsonl` lines
/// and counts its collection dates. `GENOME_FIELDS` asks for
/// `collection_date`, and BV-BRC leaves out a field it has no value for, so
/// a record without the key has no date: it is written as null (`Missing`),
/// never left `NotFetched`.
fn write_genome_page(
    w: &mut impl Write,
    page: Vec<GenomeRecord>,
    counts: &mut CollectionDateCounts,
) -> Result<usize> {
    let n = page.len();
    for mut r in page {
        if r.collection_date.not_fetched() {
            r.collection_date = CollectionDate::Missing;
        }
        counts.add(&r.collection_date);
        serde_json::to_writer(&mut *w, &r)?;
        w.write_all(b"\n")?;
    }
    Ok(n)
}

/// When BV-BRC added a genome and when its assembly was completed, as BV-BRC
/// wrote them. `grm fetch genome-dates` keeps these in their own file
/// (`bvbrc/genome_dates.jsonl`), so `genomes.jsonl` and `GENOME_FIELDS` stay
/// as pinned. They are provenance for `grm provenance` only: both are source
/// proxies (`ablate::blocks::proxy_kind`), never predictors.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct GenomeDate {
    #[serde(deserialize_with = "text")]
    pub genome_id: String,
    #[serde(default, deserialize_with = "opt_text")]
    pub date_inserted: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub completion_date: Option<String>,
}

pub const GENOME_DATE_FIELDS: &str = "genome_id,date_inserted,completion_date";

/// A CDS from the PATRIC annotation of one genome.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Feature {
    #[serde(default, deserialize_with = "opt_text")]
    pub patric_id: Option<String>,
    #[serde(deserialize_with = "text")]
    pub accession: String,
    pub start: u64,
    pub end: u64,
    #[serde(deserialize_with = "text")]
    pub strand: String,
    #[serde(default, deserialize_with = "opt_text")]
    pub gene: Option<String>,
    #[serde(default, deserialize_with = "opt_text")]
    pub product: Option<String>,
}

pub struct Client {
    agent: ureq::Agent,
    /// Wait after the first failed attempt, doubled after each later one.
    backoff: Duration,
}

impl Default for Client {
    fn default() -> Self {
        Self::new()
    }
}

impl Client {
    pub fn new() -> Self {
        let agent = ureq::AgentBuilder::new()
            .timeout_connect(Duration::from_secs(30))
            .timeout_read(Duration::from_secs(600))
            .user_agent("grm/0.2 (genomic-resistance-mapping)")
            .build();
        Client {
            agent,
            backoff: Duration::from_secs(2),
        }
    }

    /// A client whose every request fails at once, without the network.
    #[cfg(test)]
    pub(crate) fn offline() -> Self {
        let agent = ureq::AgentBuilder::new()
            .resolver(|_: &str| -> std::io::Result<Vec<std::net::SocketAddr>> {
                Err(std::io::Error::other("offline test client"))
            })
            .build();
        Client {
            agent,
            backoff: Duration::ZERO,
        }
    }

    /// Runs `op` up to five times with exponential backoff.
    fn retry<T>(&self, what: &str, mut op: impl FnMut() -> Result<T>) -> Result<T> {
        let mut delay = self.backoff;
        let mut last = None;
        for attempt in 1..=5 {
            match op() {
                Ok(v) => return Ok(v),
                Err(e) => {
                    eprintln!("  {what}: attempt {attempt} failed: {e:#}");
                    last = Some(e);
                    std::thread::sleep(delay);
                    delay *= 2;
                }
            }
        }
        Err(last
            .unwrap()
            .context(format!("{what}: giving up after 5 attempts")))
    }

    fn get_json<T: DeserializeOwned>(&self, url: &str) -> Result<Vec<T>> {
        self.retry(url, || {
            let resp = self
                .agent
                .get(url)
                .set("Accept", "application/json")
                .call()?;
            Ok(serde_json::from_reader(BufReader::new(resp.into_reader()))?)
        })
    }

    fn post_rql<T: DeserializeOwned>(&self, collection: &str, query: &str) -> Result<Vec<T>> {
        let url = format!("{API}/{collection}/");
        self.retry(collection, || {
            let resp = self
                .agent
                .post(&url)
                .set("Content-Type", "application/rqlquery+x-www-form-urlencoded")
                .set("Accept", "application/json")
                .send_string(query)?;
            Ok(serde_json::from_reader(BufReader::new(resp.into_reader()))?)
        })
    }

    /// Downloads every laboratory-method AMR record to JSON lines. Records
    /// from BV-BRC's own computational predictions are excluded at the query.
    pub fn fetch_amr(&self, out: &Path) -> Result<usize> {
        let part = out.with_extension("part");
        let mut w = BufWriter::new(File::create(&part)?);
        let mut total = 0;
        loop {
            let url = format!(
                "{API}/genome_amr/?eq(evidence,%22Laboratory%20Method%22)&sort(+id)&limit({PAGE},{total})&select({AMR_FIELDS})"
            );
            let page: Vec<AmrRecord> = self.get_json(&url)?;
            for r in &page {
                serde_json::to_writer(&mut w, r)?;
                w.write_all(b"\n")?;
            }
            total += page.len();
            eprintln!("  amr records: {total}");
            if page.len() < PAGE {
                break;
            }
        }
        w.flush()?;
        drop(w);
        std::fs::rename(&part, out)?;
        Ok(total)
    }

    /// Fetches metadata for `ids` in POST batches and writes JSON lines.
    pub fn fetch_genomes(&self, ids: &[String], out: &Path) -> Result<usize> {
        let part = out.with_extension("part");
        let mut w = BufWriter::new(File::create(&part)?);
        let mut total = 0;
        let mut counts = CollectionDateCounts::default();
        for batch in ids.chunks(1000) {
            let query = format!(
                "in(genome_id,({}))&select({GENOME_FIELDS})&limit({PAGE})",
                batch.join(",")
            );
            let page: Vec<GenomeRecord> = self.post_rql("genome", &query)?;
            total += write_genome_page(&mut w, page, &mut counts)?;
            eprintln!("  genome metadata: {total} / {}", ids.len());
        }
        eprintln!("  collection_date: {counts}");
        w.flush()?;
        drop(w);
        std::fs::rename(&part, out)?;
        Ok(total)
    }

    /// The `GenomeDate` rows BV-BRC returns for `ids`, in POST batches, in
    /// the order it returns them. A genome it has no record for has no row.
    pub fn fetch_genome_dates(&self, ids: &[String]) -> Result<Vec<GenomeDate>> {
        let mut rows = Vec::with_capacity(ids.len());
        for batch in ids.chunks(1000) {
            let query = format!(
                "in(genome_id,({}))&select({GENOME_DATE_FIELDS})&limit({PAGE})",
                batch.join(",")
            );
            rows.extend(self.post_rql::<GenomeDate>("genome", &query)?);
            eprintln!("  genome dates: {} / {}", rows.len(), ids.len());
        }
        Ok(rows)
    }

    /// Downloads one assembly as gzipped FASTA. Returns the number of
    /// nucleotides written. Writes to `.part` and renames on success.
    pub fn fetch_fasta(&self, genome_id: &str, out: &Path) -> Result<u64> {
        let url = format!(
            "{API}/genome_sequence/?eq(genome_id,{genome_id})&limit(100000)&http_accept=application/dna%2Bfasta"
        );
        let part = out.with_extension("part");
        let bases = self.retry(genome_id, || {
            let resp = self.agent.get(&url).call()?;
            let mut reader = BufReader::with_capacity(1 << 20, resp.into_reader());
            let mut gz = GzEncoder::new(BufWriter::new(File::create(&part)?), Compression::fast());
            let mut line = Vec::new();
            let (mut bases, mut headers) = (0u64, 0u64);
            loop {
                line.clear();
                if reader.read_until(b'\n', &mut line)? == 0 {
                    break;
                }
                if line.first() == Some(&b'>') {
                    headers += 1;
                } else {
                    bases += line.iter().filter(|b| b.is_ascii_alphabetic()).count() as u64;
                }
                gz.write_all(&line)?;
            }
            gz.finish()?.flush()?;
            if headers == 0 || bases == 0 {
                bail!("empty FASTA response");
            }
            Ok(bases)
        })?;
        std::fs::rename(&part, out)?;
        Ok(bases)
    }

    /// Downloads assemblies for `ids` into `dir` with `workers` parallel
    /// connections, skipping files that already exist. Returns failures.
    pub fn fetch_fasta_many(
        &self,
        ids: &[String],
        dir: &Path,
        workers: usize,
    ) -> Result<Vec<(String, String)>> {
        std::fs::create_dir_all(dir)?;
        let todo: Vec<&String> = ids
            .iter()
            .filter(|id| !fasta_path(dir, id).exists())
            .collect();
        eprintln!(
            "  fasta: {} cached, {} to download",
            ids.len() - todo.len(),
            todo.len()
        );
        let next = AtomicUsize::new(0);
        let done = AtomicUsize::new(0);
        let failures = Mutex::new(Vec::new());
        std::thread::scope(|s| {
            for _ in 0..workers.max(1) {
                s.spawn(|| {
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        let Some(id) = todo.get(i) else { break };
                        if let Err(e) = self.fetch_fasta(id, &fasta_path(dir, id)) {
                            failures
                                .lock()
                                .unwrap()
                                .push(((*id).clone(), format!("{e:#}")));
                        }
                        let d = done.fetch_add(1, Ordering::Relaxed) + 1;
                        if d.is_multiple_of(100) || d == todo.len() {
                            eprintln!("  fasta: {d} / {}", todo.len());
                        }
                    }
                });
            }
        });
        Ok(failures.into_inner().unwrap())
    }

    /// PATRIC CDS annotations for one genome.
    pub fn fetch_features(&self, genome_id: &str) -> Result<Vec<Feature>> {
        let url = format!(
            "{API}/genome_feature/?eq(genome_id,{genome_id})&eq(annotation,PATRIC)&eq(feature_type,CDS)\
&select(patric_id,accession,start,end,strand,gene,product)&limit(25000)"
        );
        self.get_json(&url)
    }

    /// Downloads a URL to a file (used for the AMRFinderPlus catalog).
    pub fn download(&self, url: &str, out: &Path) -> Result<()> {
        let part = out.with_extension("part");
        self.retry(url, || {
            let resp = self.agent.get(url).call()?;
            let mut body = Vec::new();
            resp.into_reader().read_to_end(&mut body)?;
            std::fs::write(&part, &body)?;
            Ok(())
        })?;
        std::fs::rename(&part, out).with_context(|| format!("rename to {}", out.display()))?;
        Ok(())
    }
}

pub fn fasta_path(dir: &Path, genome_id: &str) -> PathBuf {
    dir.join(format!("{genome_id}.fna.gz"))
}

/// Reads JSON lines written by the fetchers.
pub fn read_jsonl<T: DeserializeOwned>(path: &Path) -> Result<Vec<T>> {
    let reader =
        BufReader::new(File::open(path).with_context(|| format!("open {}", path.display()))?);
    let mut out = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        out.push(
            serde_json::from_str(&line).with_context(|| format!("{}:{}", path.display(), i + 1))?,
        );
    }
    Ok(out)
}

// BV-BRC returns the same field as a string in one record and a number in
// another (taxon_id, testing_standard_year, measurement_value). Accept both.
fn value_text(v: serde_json::Value) -> Option<String> {
    match v {
        serde_json::Value::Null => None,
        serde_json::Value::String(s) => {
            let t = s.trim();
            (!t.is_empty()).then(|| t.to_string())
        }
        serde_json::Value::Array(items) => {
            let parts: Vec<String> = items.into_iter().filter_map(value_text).collect();
            (!parts.is_empty()).then(|| parts.join(";"))
        }
        other => Some(other.to_string()),
    }
}

fn text<'de, D: Deserializer<'de>>(d: D) -> Result<String, D::Error> {
    value_text(serde_json::Value::deserialize(d)?)
        .ok_or_else(|| serde::de::Error::custom("missing required text"))
}

fn opt_text<'de, D: Deserializer<'de>>(d: D) -> Result<Option<String>, D::Error> {
    Ok(value_text(serde_json::Value::deserialize(d)?))
}

/// Called only when the key is present: null or empty text is `Missing`.
fn collection_date<'de, D: Deserializer<'de>>(d: D) -> Result<CollectionDate, D::Error> {
    Ok(match value_text(serde_json::Value::deserialize(d)?) {
        Some(text) => CollectionDate::parse(&text),
        None => CollectionDate::Missing,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_mixed_type_fields() {
        let r: AmrRecord = serde_json::from_str(
            r#"{"genome_id":"1280.10002","taxon_id":1280,"antibiotic":"methicillin",
                "measurement_value":"","testing_standard_year":2021,"pmid":[26686880,1]}"#,
        )
        .unwrap();
        assert_eq!(r.taxon_id.as_deref(), Some("1280"));
        assert_eq!(r.measurement_value, None);
        assert_eq!(r.testing_standard_year.as_deref(), Some("2021"));
        assert_eq!(r.pmid.as_deref(), Some("26686880;1"));
    }

    #[test]
    fn genome_dates_keep_only_their_three_fields() {
        let r: GenomeDate = serde_json::from_str(
            r#"{"genome_id":562.1,"date_inserted":"2014-12-08T22:10:02.394Z",
                "completion_date":"","genome_name":"Escherichia coli"}"#,
        )
        .unwrap();
        assert_eq!(r.genome_id, "562.1");
        assert_eq!(r.date_inserted.as_deref(), Some("2014-12-08T22:10:02.394Z"));
        assert_eq!(r.completion_date, None);
        assert_eq!(
            serde_json::to_string(&r).unwrap(),
            r#"{"genome_id":"562.1","date_inserted":"2014-12-08T22:10:02.394Z","completion_date":null}"#
        );
        // The genome record and its field list are untouched.
        assert!(!GENOME_FIELDS.contains("date_inserted") && !GENOME_FIELDS.contains("completion"));
    }

    fn day(year: u16, month: u8, day: u8) -> Option<StatedDate> {
        Some(StatedDate::Day { year, month, day })
    }

    #[test]
    fn collection_dates_parse_at_their_stated_precision_and_nothing_else() {
        let month = |year, month| Some(StatedDate::Month { year, month });
        let year = |year| Some(StatedDate::Year { year });
        let cases: &[(&str, Option<StatedDate>)] = &[
            ("2015-03-12", day(2015, 3, 12)),
            ("2016-02-29", day(2016, 2, 29)),
            ("2000-02-29", day(2000, 2, 29)),
            ("2015-12-31", day(2015, 12, 31)),
            ("2015-03", month(2015, 3)),
            ("2004-12", month(2004, 12)),
            ("2015", year(2015)),
            ("1934", year(1934)),
            // A day that does not exist is not moved to one that does.
            ("1900-02-29", None),
            ("2015-02-29", None),
            ("2015-04-31", None),
            ("2015-03-00", None),
            ("2015-03-32", None),
            // A month out of range is not read as its year.
            ("2015-13", None),
            ("2015-00", None),
            ("2004-20", None),
            // Ranges, other layouts and words stay text.
            ("2015-2016", None),
            ("2015/2016", None),
            ("2015-03-12T00:00:00Z", None),
            ("2015-03-12 ", None),
            ("2015/03/12", None),
            ("12/03/2015", None),
            ("03-12-2015", None),
            ("Mar-2015", None),
            ("2015-3-1", None),
            ("2015-3", None),
            ("20150312", None),
            ("201", None),
            ("20151", None),
            ("+201", None),
            ("-2015", None),
            ("2015-+3", None),
            ("2015-03-+1", None),
            ("missing", None),
            ("not collected", None),
            ("\u{ff12}\u{ff10}\u{ff11}\u{ff15}", None),
            ("2015-\u{e9}", None),
            ("2015-0\u{e9}", None),
            ("2015-03-\u{e9}", None),
        ];
        for &(text, want) in cases {
            let date = CollectionDate::parse(text);
            assert_eq!(date.stated(), want, "{text:?}");
            assert_eq!(
                date.precision(),
                want.map(StatedDate::precision),
                "{text:?}"
            );
            assert_eq!(date.is_unparsed(), want.is_none(), "{text:?}");
            // Parsed or not, the text is what BV-BRC wrote.
            assert_eq!(date.text().as_deref(), Some(text), "{text:?}");
            assert_eq!(
                serde_json::to_string(&date).unwrap(),
                serde_json::to_string(text).unwrap(),
                "{text:?}"
            );
        }
        assert_eq!(
            CollectionDate::parse("2015-03").stated().unwrap().year(),
            2015
        );
        assert_eq!(Precision::Month.as_str(), "month");
    }

    /// Synthetic: a `genome` response in BV-BRC's shape, written for this
    /// test. It is not a BV-BRC record. As BV-BRC does, it leaves out a
    /// field with no value and sends some values as numbers or lists.
    const SYNTHETIC_GENOME_RESPONSE: &str = r#"[
        {"genome_id":"9990001.1","genome_name":"Synthetic fixture 1","collection_year":2015,"collection_date":"2015-03-12"},
        {"genome_id":"9990001.2","genome_name":"Synthetic fixture 2","collection_year":2014,"collection_date":"2014-05"},
        {"genome_id":"9990001.3","genome_name":"Synthetic fixture 3","collection_date":2013},
        {"genome_id":"9990001.4","genome_name":"Synthetic fixture 4","collection_year":2016,"collection_date":"2015-2016"},
        {"genome_id":"9990001.5","genome_name":"Synthetic fixture 5","collection_year":2012,"collection_date":["2012-07-04"]},
        {"genome_id":"9990001.6","genome_name":"Synthetic fixture 6","collection_year":2011},
        {"genome_id":"9990001.7","genome_name":"Synthetic fixture 7","collection_date":"  "},
        {"genome_id":"9990001.8","genome_name":"Synthetic fixture 8","collection_date":"2017-02-30"},
        {"genome_id":"9990001.9","genome_name":"Synthetic fixture 9","collection_date":null}
    ]"#;

    #[test]
    fn a_synthetic_genome_response_is_written_with_its_collection_dates() {
        assert!(GENOME_FIELDS.split(',').any(|f| f == "collection_date"));
        let page: Vec<GenomeRecord> = serde_json::from_str(SYNTHETIC_GENOME_RESPONSE).unwrap();
        let mut out = Vec::new();
        let mut counts = CollectionDateCounts::default();
        assert_eq!(write_genome_page(&mut out, page, &mut counts).unwrap(), 9);
        assert_eq!(
            counts,
            CollectionDateCounts {
                day: 2,
                month: 1,
                year: 1,
                unparsed: 2,
                missing: 3
            }
        );
        assert_eq!(
            counts.to_string(),
            "day 2, month 1, year 1, unparsed and kept as text 2, missing 3"
        );
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        let field = |line: &str, key: &str| {
            serde_json::from_str::<serde_json::Value>(line).unwrap()[key].clone()
        };
        let written: Vec<serde_json::Value> =
            lines.iter().map(|l| field(l, "collection_date")).collect();
        assert_eq!(
            written,
            [
                "2015-03-12".into(),
                "2014-05".into(),
                "2013".into(),
                "2015-2016".into(),
                "2012-07-04".into(),
                serde_json::Value::Null,
                serde_json::Value::Null,
                "2017-02-30".into(),
                serde_json::Value::Null,
            ]
        );
        // Every line has the key, so a genome with no date is told apart
        // from a file fetched before the field was asked for.
        assert!(lines.iter().all(|l| l.contains(r#""collection_date":"#)));
        // The year field is as BV-BRC sent it: the range's end year stays,
        // and no year is filled in from the date.
        let years: Vec<serde_json::Value> =
            lines.iter().map(|l| field(l, "collection_year")).collect();
        assert_eq!(years[3], "2016");
        assert!(years[2].is_null() && years[7].is_null());

        let path = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(path.path(), &text).unwrap();
        let read: Vec<GenomeRecord> = read_jsonl(path.path()).unwrap();
        let dates: Vec<&CollectionDate> = read.iter().map(|g| &g.collection_date).collect();
        assert_eq!(dates[0].stated(), day(2015, 3, 12));
        assert_eq!(dates[1].precision(), Some(Precision::Month));
        assert_eq!(dates[2].precision(), Some(Precision::Year));
        assert_eq!(dates[3], &CollectionDate::Unparsed("2015-2016".into()));
        assert_eq!(dates[4].stated(), day(2012, 7, 4));
        assert_eq!(dates[5], &CollectionDate::Missing);
        assert_eq!(dates[7], &CollectionDate::Unparsed("2017-02-30".into()));
        // Written back, each line is the line that was read.
        for (record, line) in read.iter().zip(&lines) {
            assert_eq!(&serde_json::to_string(record).unwrap(), line);
        }
    }

    /// Synthetic `genomes.jsonl` lines in the form `grm fetch genomes` wrote
    /// at `main` 852ebc9, before `GENOME_FIELDS` asked for `collection_date`:
    /// every field of the record then, in order, with no value as null.
    const SYNTHETIC_LINES_BEFORE_COLLECTION_DATE: [&str; 2] = [
        r#"{"genome_id":"9990002.1","genome_name":"Synthetic fixture 1","species":"Escherichia coli","genus":"Escherichia","taxon_id":"562","genome_status":"WGS","genome_quality":"Good","genome_length":"5000000","contigs":"80","contig_n50":"150000","gc_content":"50.6","checkm_completeness":"99.8","checkm_contamination":"0.4","isolation_country":"Kenya","geographic_group":"Africa","collection_year":"2015","host_name":"Human, Homo sapiens","isolation_source":"urine","mlst":"MLST.Escherichia_coli_1.131","cgmlst_hc5":null,"cgmlst_hc20":null,"cgmlst_hc100":null,"sequencing_platform":"Illumina","assembly_method":"SPAdes","biosample_accession":"SAMN00000001","assembly_accession":"GCA_000000001.1"}"#,
        r#"{"genome_id":"9990002.2","genome_name":null,"species":null,"genus":null,"taxon_id":null,"genome_status":null,"genome_quality":null,"genome_length":null,"contigs":null,"contig_n50":null,"gc_content":null,"checkm_completeness":null,"checkm_contamination":null,"isolation_country":null,"geographic_group":null,"collection_year":null,"host_name":null,"isolation_source":null,"mlst":null,"cgmlst_hc5":null,"cgmlst_hc20":null,"cgmlst_hc100":null,"sequencing_platform":null,"assembly_method":null,"biosample_accession":null,"assembly_accession":null}"#,
    ];

    #[test]
    fn genome_lines_fetched_before_collection_date_read_and_write_as_before() {
        for line in SYNTHETIC_LINES_BEFORE_COLLECTION_DATE {
            let record: GenomeRecord = serde_json::from_str(line).unwrap();
            assert_eq!(record.collection_date, CollectionDate::NotFetched);
            assert_eq!(record.collection_date.text(), None);
            assert_eq!(record.collection_date.precision(), None);
            assert!(!record.collection_date.is_unparsed());
            assert_eq!(serde_json::to_string(&record).unwrap(), line);
        }
        let first: GenomeRecord =
            serde_json::from_str(SYNTHETIC_LINES_BEFORE_COLLECTION_DATE[0]).unwrap();
        assert_eq!(first.collection_year.as_deref(), Some("2015"));
        assert_eq!(first.isolation_country.as_deref(), Some("Kenya"));
        // A record made in code without the field is written as before too.
        let made = GenomeRecord {
            genome_id: "g1".into(),
            ..Default::default()
        };
        assert!(
            !serde_json::to_string(&made)
                .unwrap()
                .contains("collection_date")
        );
    }

    #[test]
    fn collection_date_stays_off_every_label_record() {
        // Synthetic: a `genome_amr` record in BV-BRC's shape, with a stray
        // collection_date that no lab record should carry.
        let r: AmrRecord = serde_json::from_str(
            r#"{"genome_id":"9990003.1","antibiotic":"ampicillin","resistant_phenotype":"Resistant",
                "date_inserted":"2019-06-01T10:00:00Z","collection_date":"2015-03-12"}"#,
        )
        .unwrap();
        assert!(!serde_json::to_string(&r).unwrap().contains("collection"));
        assert!(!AMR_FIELDS.contains("collection"));
        assert_eq!(
            GENOME_DATE_FIELDS,
            "genome_id,date_inserted,completion_date"
        );
        let pair = serde_json::to_string(&crate::curate::PairProvenance::default()).unwrap();
        assert!(!pair.contains("collection"), "{pair}");
        let phenotype = serde_json::to_string(&crate::curate::Phenotype {
            genome_id: "9990003.1".into(),
            species: "Escherichia coli".into(),
            antibiotic: "ampicillin".into(),
            reported: None,
            conflict: false,
            mic: None,
            mic_spread: 0,
            label: None,
            label_from_mic: false,
            methods: vec![],
            standards: vec![],
            n_records: 1,
        })
        .unwrap();
        assert!(!phenotype.contains("collection"), "{phenotype}");
    }
}
