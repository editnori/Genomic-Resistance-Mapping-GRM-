//! From selected features back to biology.
//!
//! A feature is a set of k-mers with one presence pattern. To make it
//! readable, each feature is:
//! - checked against the AMRFinderPlus catalog (is any k-mer inside a known
//!   resistance gene?);
//! - located on a reference assembly that carries it, then mapped onto the
//!   PATRIC CDS annotation of that assembly: gene, product, codon range, or
//!   the nearest gene when the k-mers sit between genes (promoters);
//! - given flanking sequence, which is what a PCR primer design would start
//!   from.
//!
//! Reference assemblies are chosen greedily to cover the most features with
//! the fewest genomes, so few annotation downloads are needed.

use crate::bits;
use crate::bvbrc::{self, Client, Feature};
use crate::catalog::Catalog;
use crate::dna::KmerSpace;
use crate::fasta;
use crate::layout::Layout;
use crate::matrix::{Dataset, FeatureMatrix};
use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashMap};

/// Most k-mers per feature used for locating (a large accessory gene can
/// carry tens of thousands; a sample places it just as well).
const LOCATE_KMERS: usize = 400;
const CONTEXT: usize = 60;
const MAX_REFERENCES: usize = 10;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Locus {
    pub genome_id: String,
    pub contig: String,
    /// 1-based inclusive span of the feature's k-mers on the contig.
    pub start: u64,
    pub end: u64,
    pub kmers_placed: usize,
    pub gene: Option<String>,
    pub product: Option<String>,
    pub patric_id: Option<String>,
    pub gene_strand: Option<String>,
    /// Codon range within the gene (1-based), when inside a CDS.
    pub codons: Option<(u64, u64)>,
    /// Nearest CDS and distance when the span is intergenic.
    pub nearest: Option<(String, i64)>,
    /// Contig sequence around the span (forward strand), span in uppercase.
    pub context: String,
    /// Why the genome has no PATRIC annotation here: its download failed, so
    /// gene, product, patric_id, gene_strand, codons, and nearest are None.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub annotation_error: Option<String>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FeatureNote {
    pub feature: u32,
    pub n_kmers: usize,
    /// One k-mer of the feature, as it reads on the reference contig.
    pub kmer: String,
    pub catalog_families: Vec<String>,
    pub catalog_products: Vec<String>,
    pub locus: Option<Locus>,
    pub n_resistant: u32,
    pub n_susceptible: u32,
}

fn feature_kmers(matrix: &FeatureMatrix, ds: &Dataset, f: u32) -> (usize, Vec<u64>) {
    let mut total = 0;
    let mut out = Vec::new();
    for &p in ds.patterns_of(f as usize) {
        let ks = matrix.kmers_of(p as usize);
        total += ks.len();
        out.extend_from_slice(ks);
    }
    if out.len() > LOCATE_KMERS {
        let step = out.len() as f64 / LOCATE_KMERS as f64;
        out = (0..LOCATE_KMERS)
            .map(|i| out[(i as f64 * step) as usize])
            .collect();
    }
    (total, out)
}

fn load_features(layout: &Layout, client: &Client, genome_id: &str) -> Result<Vec<Feature>> {
    let path = layout.features_dir().join(format!("{genome_id}.json"));
    if let Ok(bytes) = std::fs::read(&path)
        && let Ok(f) = serde_json::from_slice(&bytes)
    {
        return Ok(f);
    }
    let features = client.fetch_features(genome_id)?;
    std::fs::create_dir_all(layout.features_dir())?;
    std::fs::write(&path, serde_json::to_vec(&features)?)?;
    Ok(features)
}

struct Hit {
    contig: usize,
    pos: u64,
    kmer: String,
}

pub fn annotate(
    layout: &Layout,
    client: &Client,
    matrix: &FeatureMatrix,
    ds: &Dataset,
    features: &BTreeSet<u32>,
    catalog: Option<&Catalog>,
    contigs_of: &HashMap<String, u32>,
) -> Result<Vec<FeatureNote>> {
    let space = KmerSpace::new(matrix.summary.k);
    let mut kmers_of: HashMap<u32, Vec<u64>> = HashMap::new();
    let mut notes: HashMap<u32, FeatureNote> = HashMap::new();
    for &f in features {
        let (total, ks) = feature_kmers(matrix, ds, f);
        let mut families = BTreeSet::new();
        let mut products = BTreeSet::new();
        if let Some(cat) = catalog {
            for &k in &ks {
                for g in cat.genes_containing(k) {
                    let gene = &cat.genes[g as usize];
                    families.insert(gene.family.clone());
                    products.insert(gene.product.clone());
                }
            }
        }
        let row = ds.x.row(f as usize);
        notes.insert(
            f,
            FeatureNote {
                feature: f,
                n_kmers: total,
                kmer: ks
                    .first()
                    .map(|&k| space.decode(space.unmix(k)))
                    .unwrap_or_default(),
                catalog_families: families.into_iter().collect(),
                catalog_products: products.into_iter().collect(),
                locus: None,
                n_resistant: bits::and_count(row, ds.y.words()),
                n_susceptible: bits::count(row) - bits::and_count(row, ds.y.words()),
            },
        );
        kmers_of.insert(f, ks);
    }

    // Greedy choice of reference genomes: most uncovered features first,
    // fewer contigs breaking ties.
    let mut uncovered: BTreeSet<u32> = features.clone();
    let mut references: Vec<(usize, Vec<u32>)> = Vec::new();
    while !uncovered.is_empty() && references.len() < MAX_REFERENCES {
        let best = (0..ds.n_samples())
            .map(|s| {
                let covers: Vec<u32> = uncovered
                    .iter()
                    .copied()
                    .filter(|&f| ds.x.get(f as usize, s))
                    .collect();
                let genome = &matrix.summary.genomes[ds.columns[s]];
                (
                    covers.len(),
                    std::cmp::Reverse(contigs_of.get(genome).copied().unwrap_or(u32::MAX)),
                    s,
                    covers,
                )
            })
            .max_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)).then(b.2.cmp(&a.2)));
        let Some((n, _, s, covers)) = best else { break };
        if n == 0 {
            break;
        }
        for f in &covers {
            uncovered.remove(f);
        }
        references.push((s, covers));
    }

    for (sample, covers) in references {
        let genome_id = matrix.summary.genomes[ds.columns[sample]].clone();
        let mut wanted: HashMap<u64, Vec<u32>> = HashMap::new();
        for &f in &covers {
            for &k in &kmers_of[&f] {
                wanted.entry(k).or_default().push(f);
            }
        }
        let mut contig_names = Vec::new();
        let mut contig_seqs: Vec<Vec<u8>> = Vec::new();
        let mut hits: HashMap<u32, Vec<Hit>> = HashMap::new();
        fasta::for_each_record(
            &bvbrc::fasta_path(&layout.fasta_dir(), &genome_id),
            |name, seq| {
                let ci = contig_names.len();
                contig_names.push(name.trim_start_matches("accn|").to_string());
                let k = space.k as u64;
                // Track positions: for_each_canonical skips windows over N, so
                // walk the sequence to recover each k-mer's end coordinate.
                let mut end_positions = Vec::new();
                let mut filled = 0u64;
                for (i, &b) in seq.iter().enumerate() {
                    if crate::dna::CODE[b as usize] > 3 {
                        filled = 0;
                    } else {
                        filled += 1;
                        if filled >= k {
                            end_positions.push(i as u64);
                        }
                    }
                }
                let mut idx = 0;
                space.for_each_canonical(seq, |x| {
                    if let Some(fs) = wanted.get(&space.mix(x)) {
                        let pos = end_positions[idx] + 2 - k; // 1-based start
                        let text = String::from_utf8_lossy(
                            &seq[(pos - 1) as usize..(pos - 1 + k) as usize],
                        )
                        .to_uppercase();
                        for &f in fs {
                            hits.entry(f).or_default().push(Hit {
                                contig: ci,
                                pos,
                                kmer: text.clone(),
                            });
                        }
                    }
                    idx += 1;
                });
                contig_seqs.push(seq.to_vec());
            },
        )?;
        // A failed download leaves this genome's loci without gene or product.
        // Each records why, `grm learn` recomputes the report (see
        // `analysis::stale_reason`), and nothing is cached, so that run
        // downloads the annotation again.
        let (annotation, annotation_error) = match load_features(layout, client, &genome_id) {
            Ok(features) => (features, None),
            Err(e) => {
                eprintln!(
                    "    {genome_id}: no PATRIC annotation, so its loci have no gene or product: {e:#}"
                );
                (Vec::new(), Some(format!("{e:#}")))
            }
        };
        for (f, mut hs) in hits {
            // Densest contig, then the span of hits on it.
            let mut per_contig: HashMap<usize, usize> = HashMap::new();
            for h in &hs {
                *per_contig.entry(h.contig).or_default() += 1;
            }
            let (&contig, &count) = per_contig
                .iter()
                .max_by_key(|(c, n)| (**n, std::cmp::Reverse(**c)))
                .unwrap();
            hs.retain(|h| h.contig == contig);
            hs.sort_by_key(|h| h.pos);
            let start = hs.first().unwrap().pos;
            let end = hs.last().unwrap().pos + space.k as u64 - 1;
            let contig_name = &contig_names[contig];
            let seq = &contig_seqs[contig];
            let lo = (start as usize - 1).saturating_sub(CONTEXT);
            let hi = (end as usize + CONTEXT).min(seq.len());
            let context: String = (lo..hi)
                .map(|i| {
                    let c = seq[i] as char;
                    if i + 1 >= start as usize && i < end as usize {
                        c.to_ascii_uppercase()
                    } else {
                        c.to_ascii_lowercase()
                    }
                })
                .collect();
            let on_contig: Vec<&Feature> = annotation
                .iter()
                .filter(|a| &a.accession == contig_name)
                .collect();
            let inside = on_contig
                .iter()
                .filter(|a| a.start <= end && a.end >= start)
                .min_by_key(|a| a.end - a.start);
            let (gene, product, patric_id, gene_strand, codons, nearest) = match inside {
                Some(a) => {
                    let (c1, c2) = if a.strand == "+" {
                        (
                            (start.max(a.start) - a.start) / 3 + 1,
                            (end.min(a.end) - a.start) / 3 + 1,
                        )
                    } else {
                        (
                            (a.end - end.min(a.end)) / 3 + 1,
                            (a.end - start.max(a.start)) / 3 + 1,
                        )
                    };
                    (
                        a.gene.clone(),
                        a.product.clone(),
                        a.patric_id.clone(),
                        Some(a.strand.clone()),
                        Some((c1, c2)),
                        None,
                    )
                }
                None => {
                    let near = on_contig
                        .iter()
                        .map(|a| {
                            let d = if a.end < start {
                                start as i64 - a.end as i64
                            } else {
                                a.start as i64 - end as i64
                            };
                            (d.abs(), d, *a)
                        })
                        .min_by_key(|(abs, _, _)| *abs);
                    match near {
                        Some((_, d, a)) => (
                            a.gene.clone(),
                            a.product.clone(),
                            a.patric_id.clone(),
                            Some(a.strand.clone()),
                            None,
                            Some((a.gene.clone().or(a.product.clone()).unwrap_or_default(), d)),
                        ),
                        None => (None, None, None, None, None, None),
                    }
                }
            };
            if let Some(note) = notes.get_mut(&f) {
                note.kmer = hs[0].kmer.clone();
                note.locus = Some(Locus {
                    genome_id: genome_id.clone(),
                    contig: contig_name.clone(),
                    start,
                    end,
                    kmers_placed: count,
                    gene,
                    product,
                    patric_id,
                    gene_strand,
                    codons,
                    nearest,
                    context,
                    annotation_error: annotation_error.clone(),
                });
            }
        }
    }
    let mut out: Vec<FeatureNote> = notes.into_values().collect();
    out.sort_by_key(|n| n.feature);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kset;

    #[test]
    fn a_failed_annotation_download_is_recorded_on_the_locus_and_not_cached() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path());
        let space = KmerSpace::new(15);
        // g0 carries an accessory block, the only feature, that g1 lacks.
        let core = "ACGTTGCAAGGCTTACCGATAGCTAGGCTTAACGGTACCATGCAGTTCAGGA";
        let accessory = "TTTGGGCCCAAATGTGTCACAGGATCCATGCAAGCTTGGACTGACCTAGG";
        std::fs::create_dir_all(layout.fasta_dir()).unwrap();
        let ksets: Vec<_> = [("g0", format!("{core}NN{accessory}")), ("g1", core.into())]
            .into_iter()
            .map(|(id, seq)| {
                let fasta = bvbrc::fasta_path(&layout.fasta_dir(), id);
                std::fs::write(&fasta, format!(">c\n{seq}\n")).unwrap();
                let (values, _) = kset::count_file(&fasta, &space, 1).unwrap();
                let path = dir.path().join(format!("{id}.kset"));
                kset::write(&path, id, &space, &values).unwrap();
                path
            })
            .collect();
        let matrix = crate::matrix::build(&ksets, 1, 8).unwrap();
        let ds = Dataset::project(&matrix, vec![0, 1], &[true, false], 1);
        let locus = || {
            let notes = annotate(
                &layout,
                &Client::offline(),
                &matrix,
                &ds,
                &BTreeSet::from([0]),
                None,
                &HashMap::new(),
            )
            .unwrap();
            serde_json::to_value(&notes[0].locus).unwrap()
        };
        // No features are cached for g0, and BV-BRC cannot be reached.
        let missing = locus();
        assert_eq!(missing["genome_id"], "g0");
        assert_eq!(missing["gene"], serde_json::Value::Null);
        let error = missing["annotation_error"].as_str().unwrap_or_default();
        assert!(error.contains("giving up after 5 attempts"), "{missing}");
        assert!(!layout.features_dir().join("g0.json").exists());
        // Once the annotation is cached, the locus has its gene and no error.
        std::fs::create_dir_all(layout.features_dir()).unwrap();
        std::fs::write(
            layout.features_dir().join("g0.json"),
            r#"[{"accession": "c", "start": 1, "end": 200, "strand": "+", "gene": "blaX"}]"#,
        )
        .unwrap();
        let annotated = locus();
        assert_eq!(annotated["gene"], "blaX");
        assert!(annotated.get("annotation_error").is_none(), "{annotated}");
    }
}
