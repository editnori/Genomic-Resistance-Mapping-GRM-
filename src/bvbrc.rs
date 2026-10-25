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
collection_year,host_name,isolation_source,mlst,cgmlst_hc5,cgmlst_hc20,cgmlst_hc100,sequencing_platform,\
assembly_method,biosample_accession,assembly_accession";

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
        for batch in ids.chunks(1000) {
            let query = format!(
                "in(genome_id,({}))&select({GENOME_FIELDS})&limit({PAGE})",
                batch.join(",")
            );
            let page: Vec<GenomeRecord> = self.post_rql("genome", &query)?;
            for r in &page {
                serde_json::to_writer(&mut w, r)?;
                w.write_all(b"\n")?;
            }
            total += page.len();
            eprintln!("  genome metadata: {total} / {}", ids.len());
        }
        w.flush()?;
        drop(w);
        std::fs::rename(&part, out)?;
        Ok(total)
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
}
