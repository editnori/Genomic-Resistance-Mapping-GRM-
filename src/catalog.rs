//! Known-mechanism baseline from the NCBI AMRFinderPlus reference catalog.
//!
//! Acquired resistance genes (AMR_CDS.fa, joined to ReferenceGeneCatalog.txt
//! for class and subclass) are detected in a genome by 31-mer containment:
//! the share of a reference allele's k-mers present in the genome's kset.
//! Containment c relates to nucleotide identity roughly as c = identity^k, so
//! the 0.5 call threshold corresponds to about 97.8% identity over the gene.
//!
//! This baseline sees acquired genes only. Point mutations (gyrA, rpoB, katG,
//! porin loss) are invisible to it, which is exactly where k-mer learning has
//! room to add something.

use crate::dna::KmerSpace;
use crate::fasta;
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

pub const CATALOG_URL: &str =
    "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest";
pub const CATALOG_FILES: [&str; 3] = ["AMR_CDS.fa", "ReferenceGeneCatalog.txt", "version.txt"];
pub const CALL_CONTAINMENT: f32 = 0.5;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Gene {
    pub symbol: String,
    pub family: String,
    /// The family, split by the subfamily its product names ("OXA-48
    /// family ..."), so blaOXA-48-like and blaOXA-1 are different groups.
    pub group: String,
    pub product: String,
    pub class: String,
    pub subclass: String,
    pub n_kmers: u32,
}

/// Group key for a catalog allele: `family` alone, or `family (X family)`
/// when the product names a subfamily X.
pub fn group_of(family: &str, product: &str) -> String {
    let words: Vec<&str> = product.split_whitespace().collect();
    match words.iter().position(|w| *w == "family") {
        Some(i) if i > 0 => format!("{family} ({} family)", words[i - 1]),
        _ => family.to_string(),
    }
}

pub struct Catalog {
    pub version: String,
    pub genes: Vec<Gene>,
    /// (mixed k-mer, gene index), sorted by k-mer.
    index: Vec<(u64, u32)>,
    /// Each gene's sorted k-mers, for overlap between candidate alleles.
    kmers: Vec<Vec<u64>>,
}

impl Catalog {
    pub fn load(dir: &Path, space: &KmerSpace) -> Result<Catalog> {
        let version_path = dir.join("version.txt");
        let version = std::fs::read_to_string(&version_path)
            .with_context(|| format!("run `grm fetch catalog` first ({})", version_path.display()))?
            .trim()
            .to_string();
        anyhow::ensure!(!version.is_empty(), "{} is blank", version_path.display());
        let table_path = dir.join("ReferenceGeneCatalog.txt");
        let table = std::fs::read_to_string(&table_path)
            .with_context(|| format!("run `grm fetch catalog` first ({})", dir.display()))?;
        let mut lines = table.lines();
        let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
        let col = |name: &str| {
            header.iter().position(|h| *h == name).with_context(|| {
                format!("{}: missing required column {name}", table_path.display())
            })
        };
        let (c_type, c_class, c_sub, c_ref, c_gb, c_family, c_product) = (
            col("type")?,
            col("class")?,
            col("subclass")?,
            col("refseq_protein_accession")?,
            col("genbank_protein_accession")?,
            col("gene_family")?,
            col("product_name")?,
        );
        let mut by_protein: HashMap<String, (String, String, String, String)> = HashMap::new();
        for line in lines {
            let f: Vec<&str> = line.split('\t').collect();
            let get = |c: usize| f.get(c).copied().unwrap_or("").to_string();
            if get(c_type) != "AMR" {
                continue;
            }
            let entry = (get(c_family), get(c_product), get(c_class), get(c_sub));
            for acc in [get(c_ref), get(c_gb)] {
                if !acc.is_empty() {
                    by_protein.insert(acc, entry.clone());
                }
            }
        }

        let mut genes = Vec::new();
        let mut index = Vec::new();
        let mut gene_kmers: Vec<Vec<u64>> = Vec::new();
        // Header: protein acc | nucleotide acc | . | . | allele | family | product
        fasta::for_each_record(&dir.join("AMR_CDS.fa"), |name, seq| {
            let fields: Vec<&str> = name.split('|').collect();
            let Some((family, product, class, subclass)) =
                fields.first().and_then(|acc| by_protein.get(*acc))
            else {
                return;
            };
            let mut kmers = Vec::new();
            space.for_each_canonical(seq, |x| kmers.push(space.mix(x)));
            kmers.sort_unstable();
            kmers.dedup();
            if kmers.is_empty() {
                return;
            }
            let gi = genes.len() as u32;
            genes.push(Gene {
                symbol: fields.get(4).unwrap_or(&"").to_string(),
                family: family.clone(),
                group: group_of(family, product),
                product: product.clone(),
                class: class.clone(),
                subclass: subclass.clone(),
                n_kmers: kmers.len() as u32,
            });
            index.extend(kmers.iter().map(|&k| (k, gi)));
            gene_kmers.push(kmers);
        })?;
        anyhow::ensure!(
            !genes.is_empty(),
            "{}: no usable genes after joining ReferenceGeneCatalog.txt and AMR_CDS.fa",
            dir.display()
        );
        index.sort_unstable();
        Ok(Catalog {
            version,
            genes,
            index,
            kmers: gene_kmers,
        })
    }

    /// Identifies both the catalog release and the shape of the gene calls
    /// derived from it, so cached calls are redone when either changes.
    pub fn calls_version(&self) -> String {
        format!("{}|group-calls-v2", self.version)
    }

    /// Share of gene `b`'s k-mers that gene `a` also has.
    fn shared_share(&self, a: u32, b: u32) -> f32 {
        let (x, y) = (&self.kmers[a as usize], &self.kmers[b as usize]);
        let (mut i, mut j, mut shared) = (0, 0, 0u32);
        while i < x.len() && j < y.len() {
            match x[i].cmp(&y[j]) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    shared += 1;
                    i += 1;
                    j += 1;
                }
            }
        }
        shared as f32 / y.len().max(1) as f32
    }

    /// Catalog genes whose sequence contains this mixed k-mer.
    pub fn genes_containing(&self, kmer: u64) -> impl Iterator<Item = u32> + '_ {
        let start = self.index.partition_point(|e| e.0 < kmer);
        self.index[start..]
            .iter()
            .take_while(move |e| e.0 == kmer)
            .map(|e| e.1)
    }

    /// Containment of every catalog gene with at least `min` of its k-mers
    /// in `genome` (sorted mixed k-mers), best first.
    pub fn detect(&self, genome: &[u64], min: f32) -> Vec<(u32, f32)> {
        let mut hits = vec![0u32; self.genes.len()];
        let (mut i, mut j) = (0, 0);
        while i < self.index.len() && j < genome.len() {
            let k = self.index[i].0;
            match k.cmp(&genome[j]) {
                std::cmp::Ordering::Less => i += 1,
                std::cmp::Ordering::Greater => j += 1,
                std::cmp::Ordering::Equal => {
                    while i < self.index.len() && self.index[i].0 == k {
                        hits[self.index[i].1 as usize] += 1;
                        i += 1;
                    }
                    j += 1;
                }
            }
        }
        let mut out: Vec<(u32, f32)> = hits
            .iter()
            .enumerate()
            .filter_map(|(g, &h)| {
                let c = h as f32 / self.genes[g].n_kmers as f32;
                (c >= min).then_some((g as u32, c))
            })
            .collect();
        out.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
        out
    }
}

/// A detected gene group in one genome: its best-matching allele.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FamilyHit {
    pub family: String,
    #[serde(default)]
    pub group: String,
    pub allele: String,
    pub containment: f32,
    pub class: String,
    pub subclass: String,
}

/// Calls genes from allele hits, then collapses them to one call per group.
///
/// Near-identical alleles from different catalog families (aac(6')-Ib and
/// aac(6')-Ib-cr differ by two SNPs) all clear the containment threshold
/// for the same locus. Hits are accepted best first (containment, then gene
/// length), and a hit is dropped when an accepted gene already holds most of
/// its k-mers: the locus is explained by the better-matching allele.
pub fn family_calls(catalog: &Catalog, hits: &[(u32, f32)]) -> Vec<FamilyHit> {
    let mut ranked: Vec<(u32, f32)> = hits
        .iter()
        .copied()
        .filter(|&(_, c)| c >= CALL_CONTAINMENT)
        .collect();
    ranked.sort_by(|a, b| {
        b.1.total_cmp(&a.1)
            .then(
                catalog.genes[b.0 as usize]
                    .n_kmers
                    .cmp(&catalog.genes[a.0 as usize].n_kmers),
            )
            .then(a.0.cmp(&b.0))
    });
    let mut accepted: Vec<(u32, f32)> = Vec::new();
    for (g, c) in ranked {
        if accepted
            .iter()
            .all(|&(a, _)| catalog.shared_share(a, g) <= 0.5)
        {
            accepted.push((g, c));
        }
    }
    let mut best: HashMap<&str, (u32, f32)> = HashMap::new();
    for (g, c) in accepted {
        let group = catalog.genes[g as usize].group.as_str();
        let e = best.entry(group).or_insert((g, c));
        if c > e.1 {
            *e = (g, c);
        }
    }
    let mut out: Vec<FamilyHit> = best
        .into_values()
        .map(|(g, c)| {
            let gene = &catalog.genes[g as usize];
            FamilyHit {
                family: gene.family.clone(),
                group: gene.group.clone(),
                allele: gene.symbol.clone(),
                containment: c,
                class: gene.class.clone(),
                subclass: gene.subclass.clone(),
            }
        })
        .collect();
    out.sort_by(|a, b| a.group.cmp(&b.group));
    out
}

/// Gene requirements under which the catalog calls a drug resistant: every
/// group must be met by some detected gene; a group lists alternatives as
/// `C:<class token>` or `S:<subclass token>`. None when the drug's
/// resistance is mainly mutational and no acquired-gene rule applies.
pub fn drug_rule(drug: &str) -> Option<Vec<Vec<&'static str>>> {
    let any = |t: &[&'static str]| Some(vec![t.to_vec()]);
    let (principal, inhibited) = crate::drugs::components(drug);
    let class = crate::drugs::class_of(drug);
    if inhibited {
        // An inhibitor restores the partner against narrow-spectrum enzymes,
        // so only inhibitor-resistant mechanisms count.
        return match class.as_str() {
            "penicillin + inhibitor" => any(&["S:CEPHALOSPORIN", "S:CARBAPENEM", "S:METHICILLIN"]),
            "cephalosporin + inhibitor" | "carbapenem + inhibitor" => any(&["S:CARBAPENEM"]),
            _ => None,
        };
    }
    if class == "folate pathway combination" {
        return Some(vec![vec!["C:TRIMETHOPRIM"], vec!["C:SULFONAMIDE"]]);
    }
    if class == "streptogramin" {
        return any(&["C:STREPTOGRAMIN"]);
    }
    let [drug] = principal.as_slice() else {
        return None;
    };
    match *drug {
        "meropenem" | "imipenem" | "ertapenem" | "doripenem" | "biapenem" => any(&["S:CARBAPENEM"]),
        "ceftriaxone" | "cefotaxime" | "ceftazidime" | "cefepime" | "cefpodoxime" | "ceftiofur"
        | "cefixime" | "cefoperazone" | "cefazolin" | "cephalothin" | "cefuroxime"
        | "cefotetan" | "cefmetazole" | "cefoxitin" => {
            any(&["S:CEPHALOSPORIN", "S:CARBAPENEM", "S:METHICILLIN"])
        }
        "aztreonam" => any(&["S:CEPHALOSPORIN"]),
        "ampicillin" | "amoxicillin" | "penicillin" | "piperacillin" | "ticarcillin"
        | "carbenicillin" => any(&["C:BETA-LACTAM"]),
        "oxacillin" | "methicillin" => any(&["S:METHICILLIN"]),
        "gentamicin" => any(&["S:GENTAMICIN", "S:AMINOGLYCOSIDE"]),
        "amikacin" => any(&["S:AMIKACIN", "S:AMINOGLYCOSIDE"]),
        "tobramycin" => any(&["S:TOBRAMYCIN", "S:AMINOGLYCOSIDE"]),
        "kanamycin" => any(&["S:KANAMYCIN", "S:AMINOGLYCOSIDE"]),
        "streptomycin" => any(&["S:STREPTOMYCIN"]),
        "spectinomycin" => any(&["S:SPECTINOMYCIN"]),
        "ciprofloxacin" | "levofloxacin" | "moxifloxacin" | "ofloxacin" | "norfloxacin"
        | "nalidixic acid" | "enrofloxacin" => any(&["C:QUINOLONE"]),
        "tetracycline" | "doxycycline" | "minocycline" => any(&["C:TETRACYCLINE"]),
        "tigecycline" => any(&["S:TIGECYCLINE"]),
        "azithromycin" | "erythromycin" | "clarithromycin" => any(&["C:MACROLIDE"]),
        "clindamycin" | "lincomycin" => any(&["C:LINCOSAMIDE"]),
        "trimethoprim" => any(&["C:TRIMETHOPRIM"]),
        "sulfisoxazole" | "sulfamethoxazole" | "sulfathiazole" => any(&["C:SULFONAMIDE"]),
        "chloramphenicol" => any(&["S:CHLORAMPHENICOL", "S:PHENICOL"]),
        "florfenicol" => any(&["S:FLORFENICOL", "S:PHENICOL"]),
        "colistin" | "polymyxin b" => any(&["C:COLISTIN"]),
        "vancomycin" => any(&["S:VANCOMYCIN"]),
        "teicoplanin" => any(&["S:TEICOPLANIN"]),
        "fosfomycin" => any(&["C:FOSFOMYCIN"]),
        "nitrofurantoin" => any(&["C:NITROFURAN"]),
        "fusidic acid" => any(&["C:FUSIDIC ACID"]),
        "mupirocin" => any(&["C:MUPIROCIN"]),
        "linezolid" => any(&["C:OXAZOLIDINONE"]),
        "daptomycin" => any(&["C:LIPOPEPTIDE"]),
        _ => None,
    }
}

fn token_match(token: &str, class: &str, subclass: &str) -> bool {
    let (kind, want) = token.split_at(2);
    let field = if kind == "C:" { class } else { subclass };
    field.split('/').any(|t| t.trim() == want)
}

/// Families among `hits` that satisfy one of the rule's groups, and whether
/// all groups are satisfied (the catalog's resistant call).
pub fn explain(rule: &[Vec<&str>], hits: &[FamilyHit]) -> (bool, Vec<String>) {
    let mut used = Vec::new();
    let mut all = true;
    for group in rule {
        let matching: Vec<&FamilyHit> = hits
            .iter()
            .filter(|h| group.iter().any(|t| token_match(t, &h.class, &h.subclass)))
            .collect();
        if matching.is_empty() {
            all = false;
        }
        used.extend(matching.into_iter().map(|h| h.group.clone()));
    }
    used.sort();
    used.dedup();
    (all, used)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn catalog_fixture() -> tempfile::TempDir {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(
            dir.path().join("ReferenceGeneCatalog.txt"),
            "type\tclass\tsubclass\trefseq_protein_accession\tgenbank_protein_accession\tgene_family\tproduct_name\n\
             AMR\tfixture class\tfixture subclass\tREF\tGB\tfixture_family\tfixture product\n",
        ).unwrap();
        std::fs::write(
            dir.path().join("AMR_CDS.fa"),
            format!(
                ">REF|nucleotide|1|1|fixture_gene|fixture_family|fixture product\n{}\n",
                "ACGT".repeat(16)
            ),
        )
        .unwrap();
        std::fs::write(dir.path().join("version.txt"), "fixture-v1\n").unwrap();
        dir
    }

    #[test]
    fn catalog_requires_each_column() {
        let dir = catalog_fixture();
        let path = dir.path().join("ReferenceGeneCatalog.txt");
        let table = std::fs::read_to_string(&path).unwrap();
        let rows: Vec<Vec<&str>> = table
            .lines()
            .map(|line| line.split('\t').collect())
            .collect();
        for (missing, name) in rows[0].iter().enumerate() {
            let truncated = rows
                .iter()
                .map(|row| {
                    row.iter()
                        .enumerate()
                        .filter_map(|(i, value)| (i != missing).then_some(*value))
                        .collect::<Vec<_>>()
                        .join("\t")
                })
                .collect::<Vec<_>>()
                .join("\n");
            std::fs::write(&path, truncated).unwrap();
            let err = Catalog::load(dir.path(), &KmerSpace::new(31))
                .err()
                .expect(name);
            let message = format!("{err:#}");
            assert!(message.contains("missing required column"), "{message}");
            assert!(message.contains(name), "{message}");
            assert!(message.contains(path.to_str().unwrap()), "{message}");
        }
    }

    #[test]
    fn catalog_requires_a_version_file() {
        let dir = catalog_fixture();
        let path = dir.path().join("version.txt");
        std::fs::remove_file(&path).unwrap();
        let err = Catalog::load(dir.path(), &KmerSpace::new(31))
            .err()
            .expect("missing version");
        assert_eq!(
            err.downcast_ref::<std::io::Error>().unwrap().kind(),
            std::io::ErrorKind::NotFound
        );
        assert!(err.to_string().contains(path.to_str().unwrap()), "{err:#}");
    }

    #[test]
    fn catalog_requires_a_nonblank_version() {
        let dir = catalog_fixture();
        let path = dir.path().join("version.txt");
        for blank in ["", " \t\r\n"] {
            std::fs::write(&path, blank).unwrap();
            let err = Catalog::load(dir.path(), &KmerSpace::new(31))
                .err()
                .expect("blank version");
            assert!(err.to_string().contains("blank"), "{err:#}");
            assert!(err.to_string().contains(path.to_str().unwrap()), "{err:#}");
        }
    }

    #[test]
    fn catalog_requires_usable_genes_after_joining() {
        let dir = catalog_fixture();
        let table_path = dir.path().join("ReferenceGeneCatalog.txt");
        let fasta_path = dir.path().join("AMR_CDS.fa");
        let table = std::fs::read_to_string(&table_path).unwrap();
        let fasta = std::fs::read_to_string(&fasta_path).unwrap();
        for (case, table, fasta) in [
            (
                "header only",
                table.lines().next().unwrap().to_string(),
                fasta.clone(),
            ),
            (
                "no AMR rows",
                table.replace("\nAMR\t", "\nOTHER\t"),
                fasta.clone(),
            ),
            ("empty FASTA", table.clone(), String::new()),
            (
                "unmatched accession",
                table.clone(),
                fasta.replace(">REF|", ">UNMATCHED|"),
            ),
            ("too short", table.clone(), ">REF\nACGT\n".into()),
            (
                "no canonical kmers",
                table,
                format!(">REF\n{}\n", "N".repeat(64)),
            ),
        ] {
            std::fs::write(&table_path, table).unwrap();
            std::fs::write(&fasta_path, fasta).unwrap();
            let err = Catalog::load(dir.path(), &KmerSpace::new(31))
                .err()
                .expect(case);
            assert!(
                err.to_string().contains("no usable genes"),
                "{case}: {err:#}"
            );
            assert!(
                err.to_string().contains(dir.path().to_str().unwrap()),
                "{case}: {err:#}"
            );
        }
    }

    #[test]
    fn catalog_loads_both_accessions_and_trims_the_version() {
        let dir = catalog_fixture();
        let fasta_path = dir.path().join("AMR_CDS.fa");
        let fasta = std::fs::read_to_string(&fasta_path).unwrap();
        std::fs::write(dir.path().join("version.txt"), " \tfixture-v1\r\n").unwrap();
        for accession in ["REF", "GB"] {
            std::fs::write(
                &fasta_path,
                fasta.replace(">REF|", &format!(">{accession}|")),
            )
            .unwrap();
            let catalog = Catalog::load(dir.path(), &KmerSpace::new(31)).unwrap();
            assert_eq!(catalog.version, "fixture-v1");
            assert_eq!(catalog.genes.len(), 1);
            assert_eq!(catalog.genes[0].symbol, "fixture_gene");
            assert_eq!(catalog.genes[0].family, "fixture_family");
            assert_eq!(catalog.genes[0].product, "fixture product");
            assert_eq!(catalog.genes[0].class, "fixture class");
            assert_eq!(catalog.genes[0].subclass, "fixture subclass");
            assert!(catalog.genes[0].n_kmers > 0);
        }
    }

    fn hit(family: &str, class: &str, subclass: &str) -> FamilyHit {
        FamilyHit {
            family: family.into(),
            group: family.into(),
            allele: family.into(),
            containment: 1.0,
            class: class.into(),
            subclass: subclass.into(),
        }
    }

    #[test]
    fn near_identical_alleles_resolve_to_the_better_match() {
        let dir = tempfile::tempdir().unwrap();
        let space = KmerSpace::new(31);
        let a = "ATGAAACGTTACGGCATTGCAGGCTTACCGATAGCTTGACCTAGGCATCGATCGGATCCATGCAAGCTTGGTACCGAGCTCGAATTCACTGGCCGTCGTTTTAC";
        // b differs from a by one substitution in the middle.
        let mut b = a.as_bytes().to_vec();
        b[52] = if b[52] == b'A' { b'C' } else { b'A' };
        let b = String::from_utf8(b).unwrap();
        std::fs::write(
            dir.path().join("AMR_CDS.fa"),
            format!(">WP_A|NG_A|1|1|aacA|aacA|acetyltransferase\n{a}\n>WP_B|NG_B|1|1|aacA-cr|aacA-cr|acetyltransferase variant\n{b}\n"),
        )
        .unwrap();
        std::fs::write(
            dir.path().join("ReferenceGeneCatalog.txt"),
            "gene_family\tproduct_name\ttype\tclass\tsubclass\trefseq_protein_accession\tgenbank_protein_accession\n\
             aacA\tacetyltransferase\tAMR\tAMINOGLYCOSIDE\tAMIKACIN\tWP_A\t\n\
             aacA-cr\tacetyltransferase variant\tAMR\tAMINOGLYCOSIDE/QUINOLONE\tAMIKACIN/QUINOLONE\tWP_B\t\n",
        )
        .unwrap();
        std::fs::write(dir.path().join("version.txt"), "test\n").unwrap();
        let cat = Catalog::load(dir.path(), &space).unwrap();
        let mut genome = Vec::new();
        space.for_each_canonical(a.as_bytes(), |x| genome.push(space.mix(x)));
        genome.sort_unstable();
        genome.dedup();
        let hits = cat.detect(&genome, 0.3);
        assert_eq!(
            hits.len(),
            2,
            "both alleles clear the threshold on their own"
        );
        let calls = family_calls(&cat, &hits);
        assert_eq!(
            calls.iter().map(|c| c.family.as_str()).collect::<Vec<_>>(),
            vec!["aacA"]
        );
        assert!(!explain(&drug_rule("ciprofloxacin").unwrap(), &calls).0);
    }

    #[test]
    fn groups_split_families_by_named_subfamily() {
        assert_eq!(
            group_of("blaOXA", "OXA-48 family class D beta-lactamase OXA-1012"),
            "blaOXA (OXA-48 family)"
        );
        assert_eq!(
            group_of(
                "blaKPC",
                "carbapenem-hydrolyzing class A beta-lactamase KPC-2"
            ),
            "blaKPC"
        );
    }

    #[test]
    fn rules_match_class_and_subclass_tokens() {
        let kpc = hit("blaKPC", "BETA-LACTAM", "CARBAPENEM");
        let tem = hit("blaTEM", "BETA-LACTAM", "BETA-LACTAM");
        let sul = hit("sul1", "SULFONAMIDE", "SULFONAMIDE");
        let dfr = hit("dfrA1", "TRIMETHOPRIM", "TRIMETHOPRIM");
        let aac = hit(
            "aac(6')-Ib",
            "AMINOGLYCOSIDE",
            "AMIKACIN/KANAMYCIN/TOBRAMYCIN",
        );

        assert_eq!(
            explain(
                &drug_rule("meropenem").unwrap(),
                &[tem.clone(), kpc.clone()]
            ),
            (true, vec!["blaKPC".into()])
        );
        assert!(!explain(&drug_rule("meropenem").unwrap(), std::slice::from_ref(&tem)).0);
        assert!(explain(&drug_rule("ampicillin").unwrap(), &[tem]).0);
        // Co-trimoxazole needs both a dfr and a sul gene.
        assert!(
            !explain(
                &drug_rule("sulfamethoxazole/trimethoprim").unwrap(),
                std::slice::from_ref(&sul)
            )
            .0
        );
        assert!(
            explain(
                &drug_rule("sulfamethoxazole/trimethoprim").unwrap(),
                &[sul, dfr]
            )
            .0
        );
        assert!(explain(&drug_rule("amikacin").unwrap(), std::slice::from_ref(&aac)).0);
        assert!(!explain(&drug_rule("gentamicin").unwrap(), &[aac]).0);
        assert!(drug_rule("isoniazid").is_none());
    }

    #[test]
    fn detects_a_gene_by_containment() {
        let dir = tempfile::tempdir().unwrap();
        let space = KmerSpace::new(31);
        let gene = "ATGAAACGTTACGGCATTGCAGGCTTACCGATAGCTTGACCTAGGCATCGATCGGATCCATGCAAGCTTGG";
        std::fs::write(
            dir.path().join("AMR_CDS.fa"),
            format!(">WP_1.1|NG_1.1|1|1|blaX-1|blaX|beta-lactamase\n{gene}\n"),
        )
        .unwrap();
        std::fs::write(
            dir.path().join("ReferenceGeneCatalog.txt"),
            "allele\tgene_family\tproduct_name\ttype\tclass\tsubclass\trefseq_protein_accession\tgenbank_protein_accession\n\
             blaX-1\tblaX\tbeta-lactamase\tAMR\tBETA-LACTAM\tCARBAPENEM\tWP_1.1\t\n",
        )
        .unwrap();
        std::fs::write(dir.path().join("version.txt"), "test\n").unwrap();
        let cat = Catalog::load(dir.path(), &space).unwrap();
        assert_eq!(cat.genes.len(), 1);
        let mut genome = Vec::new();
        space.for_each_canonical(format!("TTTTTTTTTT{gene}GGGGGGGGG").as_bytes(), |x| {
            genome.push(space.mix(x))
        });
        genome.sort_unstable();
        genome.dedup();
        let hits = cat.detect(&genome, 0.1);
        assert_eq!(hits.len(), 1);
        assert!((hits[0].1 - 1.0).abs() < 1e-6);
        let fams = family_calls(&cat, &hits);
        assert_eq!(fams[0].family, "blaX");
        assert!(explain(&drug_rule("meropenem").unwrap(), &fams).0);
    }
}
