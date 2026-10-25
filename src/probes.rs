//! Simulated mutation-detecting assays.
//!
//! Xpert MTB/RIF does not look for any one mutation. Its probes cover the
//! 81-bp rifampicin resistance-determining region (RRDR) of rpoB and report
//! resistance when the sample's sequence there differs from wild type. In
//! k-mer terms: take every 31-mer that overlaps the region on a verified
//! wild-type coding sequence, and call "mutation detected" when a genome
//! carries some of them but none that covers one of the probed bases; a
//! genome with none of them gets no call. A single SNP removes up to 31
//! consecutive k-mers, so the call does not depend on which SNP it is.
//!
//! The same mechanism covers the other region-targeted assays below.
//! Synonymous changes in a probed region also read as "mutation detected",
//! which is a known source of false resistance calls in the real assays.

use crate::dna::KmerSpace;
use crate::kset::{BUCKET_BITS, KsetReader};
use anyhow::{Context, Result, bail, ensure};
use rayon::prelude::*;
use serde::Serialize;
use std::collections::BTreeMap;

/// Where a probed region sits relative to its gene.
#[derive(Clone, Copy, Debug, Serialize)]
pub enum Span {
    /// Codons (1-based, inclusive) in the gene's own numbering.
    Codons(u32, u32),
    /// Bases upstream of the start codon (1-based, inclusive): a promoter.
    Upstream(u32, u32),
}

#[derive(Clone, Copy, Debug, Serialize)]
pub struct WildType {
    pub source: &'static str,
    pub record: &'static str,
    pub cds_start: u64,
    pub cds_end: u64,
    pub strand: &'static str,
    pub cds: &'static str,
    pub residues: &'static str,
}

#[derive(Clone, Copy, Debug, Serialize)]
pub struct Region {
    pub gene: &'static str,
    pub span: Span,
    pub wild_type: Option<WildType>,
}

#[derive(Clone, Debug, Serialize)]
pub struct Assay {
    pub id: &'static str,
    pub name: &'static str,
    pub mirrors: &'static str,
    pub species: &'static str,
    pub drugs: &'static [&'static str],
    pub regions: &'static [Region],
    pub unsupported: Option<&'static str>,
}

const ECOLI_SYNONYMOUS: &str = "In E. coli, synonymous base changes in gyrA codons 83-87 and parC codons 80-84 are common, and these nucleotide probes read them as mutations: in an 80-genome check, 17 of 35 genomes with wild-type amino acids in both regions read mutation detected.";

const ECOLI_GYRA: WildType = WildType {
    source: "AMRFinderPlus 2026-08-07.1 wild-type reference WP_001281243.1",
    record: "CP031609.1",
    cds_start: 2607950,
    cds_end: 2610586,
    strand: "+",
    cds: "ATGAGCGACCTTGCGAGAGAAATTACACCGGTCAACATTGAGGAAGAGCTGAAGAGCTCCTATCTGGATTATGCGATGTCGGTCATTGTTGGCCGTGCGCTGCCAGATGTCCGAGATGGCCTGAAGCCGGTACACCGTCGCGTACTTTACGCCATGAACGTACTAGGCAATGACTGGAACAAAGCCTATAAAAAATCTGCCCGTGTCGTTGGTGACGTAATCGGTAAATACCATCCCCATGGTGACTCGGCGGTTTATGACACGATCGTCCGTATGGCGCAGCCATTCTCGC",
    residues: "SAVYD",
};

const ECOLI_PARC: WildType = WildType {
    source: "AMRFinderPlus 2026-08-07.1 wild-type reference WP_001281881.1",
    record: "U00096.3",
    cds_start: 3163715,
    cds_end: 3165973,
    strand: "-",
    cds: "ATGAGCGATATGGCAGAGCGCCTTGCGCTACATGAATTTACGGAAAACGCCTACTTAAACTACTCCATGTACGTGATCATGGACCGTGCGTTGCCGTTTATTGGTGATGGTCTGAAACCTGTTCAGCGCCGCATTGTGTATGCGATGTCTGAACTGGGCCTGAATGCCAGCGCCAAATTTAAAAAATCGGCCCGTACCGTCGGTGACGTACTGGGTAAATACCATCCGCACGGCGATAGCGCCTGTTATGAAGCGATGGTCCTGATGGCGCAACCGTTCTCTT",
    residues: "SACYE",
};

const GONORRHOEAE_GYRA: WildType = WildType {
    source: "AMRFinderPlus 2026-08-07.1 wild-type reference WP_003688895.1",
    record: "CP012027.1",
    cds_start: 478958,
    cds_end: 481708,
    strand: "+",
    cds: "ATGACCGACGCAACCATCCGCCACGACCACAAATTCGCCCTCGAAACCCTGCCCGTCAGCCTTGAAGACGAAATGCGCAAAAGCTATCTCGACTACGCCATGAGCGTCATTGTCGGGCGCGCGCTGCCGGACGTTCGCGACGGCCTAAAGCCGGTGCACCGGCGCGTACTGTACGCGATGCACGAGCTGAAAAATAACTGGAATGCCGCCTACAAAAAATCGGCGCGCATCGTCGGCGACGTCATCGGTAAATACCACCCCCACGGCGATTCCGCAGTTTACGACACCATCGTCCGTATGGCGCAAAATTTCGCTA",
    residues: "SAVYD",
};

const SALMONELLA_GYRA: WildType = WildType {
    source: "AMRFinderPlus 2026-08-07.1 wild-type reference WP_001281271.1",
    record: "AAAEUU010000006.1",
    cds_start: 382271,
    cds_end: 384907,
    strand: "+",
    cds: "ATGAGCGACCTTGCGAGAGAAATTACACCGGTCAACATTGAGGAGGAGCTGAAGAGCTCCTATCTGGATTATGCGATGTCGGTCATTGTTGGCCGTGCGCTGCCGGATGTCCGAGATGGCCTGAAGCCGGTACACCGTCGCGTACTTTACGCCATGAACGTATTGGGCAATGACTGGAACAAAGCCTATAAAAAATCTGCCCGTGTCGTTGGTGACGTAATCGGTAAATACCATCCCCACGGCGATTCCGCAGTGTATGACACCATCGTTCGTATGGCGCAGCCATTCTCGC",
    residues: "SAVYD",
};

pub const ASSAYS: &[Assay] = &[
    Assay {
        id: "xpert_rif",
        name: "rpoB RRDR wild-type probes",
        mirrors: "Cepheid Xpert MTB/RIF (Ultra): loss of wild-type sequence across rpoB codons 426-452",
        species: "Mycobacterium tuberculosis",
        drugs: &["rifampin", "rifabutin"],
        regions: &[Region {
            gene: "rpoB",
            span: Span::Codons(426, 452),
            wild_type: None,
        }],
        unsupported: None,
    },
    Assay {
        id: "xpert_xdr_inh",
        name: "katG 315 and inhA promoter probes",
        mirrors: "Cepheid Xpert MTB/XDR isoniazid targets: katG codon 315 region and the fabG1-inhA promoter",
        species: "Mycobacterium tuberculosis",
        drugs: &["isoniazid"],
        regions: &[
            Region {
                gene: "katG",
                span: Span::Codons(311, 319),
                wild_type: None,
            },
            Region {
                gene: "fabG1",
                span: Span::Upstream(1, 32),
                wild_type: None,
            },
        ],
        unsupported: None,
    },
    Assay {
        id: "xpert_xdr_eto",
        name: "inhA promoter probes",
        mirrors: "Cepheid Xpert MTB/XDR ethionamide target: the fabG1-inhA promoter",
        species: "Mycobacterium tuberculosis",
        drugs: &["ethionamide"],
        regions: &[Region {
            gene: "fabG1",
            span: Span::Upstream(1, 32),
            wild_type: None,
        }],
        unsupported: None,
    },
    Assay {
        id: "xpert_xdr_fq",
        name: "gyrA and gyrB QRDR probes",
        mirrors: "Fluoroquinolone targets as in Xpert MTB/XDR: the gyrA QRDR (codons 88-94) and the gyrB QRDR (codons 461-501, H37Rv numbering)",
        species: "Mycobacterium tuberculosis",
        drugs: &["ofloxacin", "levofloxacin", "moxifloxacin", "ciprofloxacin"],
        regions: &[
            Region {
                gene: "gyrA",
                span: Span::Codons(88, 94),
                wild_type: None,
            },
            Region {
                gene: "gyrB",
                span: Span::Codons(461, 501),
                wild_type: None,
            },
        ],
        unsupported: None,
    },
    Assay {
        id: "xpert_xdr_amk",
        name: "rrs 1401 region and eis promoter probes",
        mirrors: "Cepheid Xpert MTB/XDR second-line injectable targets: eis promoter (rrs is an RNA gene and not probed here)",
        species: "Mycobacterium tuberculosis",
        drugs: &["amikacin", "kanamycin", "capreomycin"],
        regions: &[Region {
            gene: "eis",
            span: Span::Upstream(1, 40),
            wild_type: None,
        }],
        unsupported: None,
    },
    Assay {
        id: "gc_gyra",
        name: "gyrA codons 91-95 probes",
        mirrors: "Gonorrhoea ciprofloxacin genotyping (for example SpeeDx ResistancePlus GC targets gyrA S91)",
        species: "Neisseria gonorrhoeae",
        drugs: &["ciprofloxacin"],
        regions: &[Region {
            gene: "gyrA",
            span: Span::Codons(91, 95),
            wild_type: Some(GONORRHOEAE_GYRA),
        }],
        unsupported: None,
    },
    Assay {
        id: "ecoli_qrdr",
        name: "gyrA and parC QRDR probes",
        mirrors: "Quinolone resistance-determining regions: gyrA codons 83-87 and parC codons 80-84",
        species: "Escherichia coli",
        drugs: &["ciprofloxacin", "levofloxacin", "nalidixic acid"],
        regions: &[
            Region {
                gene: "gyrA",
                span: Span::Codons(83, 87),
                wild_type: Some(ECOLI_GYRA),
            },
            Region {
                gene: "parC",
                span: Span::Codons(80, 84),
                wild_type: Some(ECOLI_PARC),
            },
        ],
        unsupported: Some(ECOLI_SYNONYMOUS),
    },
    Assay {
        id: "salmonella_qrdr",
        name: "gyrA QRDR probes",
        mirrors: "Salmonella quinolone resistance-determining region: gyrA codons 83-87",
        species: "Salmonella enterica",
        drugs: &["ciprofloxacin", "nalidixic acid"],
        regions: &[Region {
            gene: "gyrA",
            span: Span::Codons(83, 87),
            wild_type: Some(SALMONELLA_GYRA),
        }],
        unsupported: None,
    },
];

pub fn for_species(species: &str) -> Vec<&'static Assay> {
    ASSAYS.iter().filter(|a| a.species == species).collect()
}

/// One region resolved on its verified wild-type sequence.
#[derive(Clone, Debug, Serialize)]
pub struct ResolvedRegion {
    pub gene: String,
    pub record: String,
    pub source: String,
    /// 1-based inclusive record coordinates of the probed bases.
    pub start: u64,
    pub end: u64,
    pub strand: String,
    /// The probed bases, read along the gene.
    pub sequence: String,
    pub probe_kmers: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ProbeCall {
    WildType,
    MutationDetected,
    NoCall,
}

#[derive(Clone, Debug, Serialize)]
pub struct ProbeSet {
    pub assay: String,
    pub regions: Vec<ResolvedRegion>,
    /// Per genome, in the order `build_probes` was given.
    #[serde(skip)]
    pub calls: Vec<ProbeCall>,
}

#[derive(Debug)]
pub struct KsetsUnreadable;

impl std::fmt::Display for KsetsUnreadable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("probe calls need every genome's kset, which `grm count` recreates")
    }
}

fn translate(codons: &[u8]) -> String {
    const AMINO: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    codons
        .chunks(3)
        .map(|codon| {
            codon
                .iter()
                .try_fold(0, |index, &base| {
                    Some(
                        index * 4
                            + match base {
                                b'T' => 0,
                                b'C' => 1,
                                b'A' => 2,
                                b'G' => 3,
                                _ => return None,
                            },
                    )
                })
                .map_or('X', |index| AMINO[index] as char)
        })
        .collect()
}

/// Resolves a region on its verified wild-type coding sequence and returns
/// the k-mers overlapping it, one per start position.
fn region_kmers(space: &KmerSpace, region: &Region) -> Result<(ResolvedRegion, Vec<u64>)> {
    let wild = region
        .wild_type
        .with_context(|| format!("gene {} has no verified wild-type sequence", region.gene))?;
    let Span::Codons(a, b) = region.span else {
        bail!("gene {} is probed outside its coding sequence", region.gene);
    };
    let cds = wild.cds.as_bytes();
    let (lo, hi, k) = (3 * (a as usize - 1), 3 * b as usize, space.k);
    if lo + 1 < k || hi + k - 1 > cds.len() {
        bail!(
            "the {} wild-type sequence does not extend {} bases past codons {a}-{b}",
            region.gene,
            k - 1
        );
    }
    let residues = translate(&cds[lo..hi]);
    if residues != wild.residues {
        bail!(
            "the {} wild-type sequence encodes {residues} at codons {a}-{b}, not {}",
            region.gene,
            wild.residues
        );
    }
    let mut kmers = Vec::new();
    space.for_each_canonical(&cds[lo + 1 - k..hi + k - 1], |x| kmers.push(space.mix(x)));
    if kmers.len() != hi - lo + k - 1 {
        bail!(
            "the {} wild-type sequence has a base other than ACGT",
            region.gene
        );
    }
    let (start, end) = if wild.strand == "+" {
        (wild.cds_start + lo as u64, wild.cds_start + hi as u64 - 1)
    } else {
        (wild.cds_end + 1 - hi as u64, wild.cds_end - lo as u64)
    };
    Ok((
        ResolvedRegion {
            gene: region.gene.to_string(),
            record: wild.record.to_string(),
            source: wild.source.to_string(),
            start,
            end,
            strand: wild.strand.to_string(),
            sequence: String::from_utf8_lossy(&cds[lo..hi]).into_owned(),
            probe_kmers: kmers.len(),
        },
        kmers,
    ))
}

fn region_call(carried: &[bool], k: usize) -> ProbeCall {
    if !carried.contains(&true) {
        ProbeCall::NoCall
    } else if carried.windows(k).all(|covering| covering.contains(&true)) {
        ProbeCall::WildType
    } else {
        ProbeCall::MutationDetected
    }
}

fn assay_call(per_region: &[ProbeCall]) -> ProbeCall {
    if per_region.contains(&ProbeCall::MutationDetected) {
        ProbeCall::MutationDetected
    } else if per_region.contains(&ProbeCall::NoCall) {
        ProbeCall::NoCall
    } else {
        ProbeCall::WildType
    }
}

/// Which of `kmers` (in `space`) each genome carries, read from the genomes'
/// ksets one bucket at a time. A kset counted at another k holds other
/// values, so it is an error rather than a genome without the probes.
pub fn presence(
    ksets: &[std::path::PathBuf],
    kmers: &[u64],
    space: &KmerSpace,
) -> Result<Vec<Vec<bool>>> {
    let mut by_bucket: BTreeMap<usize, Vec<(usize, u64)>> = BTreeMap::new();
    for (i, &k) in kmers.iter().enumerate() {
        by_bucket
            .entry((k >> (space.bits() - BUCKET_BITS)) as usize)
            .or_default()
            .push((i, k));
    }
    ksets
        .par_iter()
        .map(|path| {
            let reader = KsetReader::open(path)?;
            ensure!(
                reader.k == space.k,
                "{} holds {}-mers, not the {}-mers these probes use",
                path.display(),
                reader.k,
                space.k
            );
            let mut has = vec![false; kmers.len()];
            let (mut bytes, mut values) = (Vec::new(), Vec::new());
            for (&b, wanted) in &by_bucket {
                values.clear();
                reader.read_bucket(b, &mut bytes, &mut values)?;
                for &(i, k) in wanted {
                    has[i] = values.binary_search(&k).is_ok();
                }
            }
            Ok(has)
        })
        .collect()
}

/// Builds the probe set for `assay` from its verified wild-type sequences and
/// calls each genome of `ksets` from its own kset alone.
pub fn build_probes(
    assay: &Assay,
    ksets: &[std::path::PathBuf],
    space: &KmerSpace,
) -> Result<ProbeSet> {
    if let Some(reason) = assay.unsupported {
        bail!("{reason}");
    }
    let mut regions = Vec::new();
    let mut kmers = Vec::new();
    for r in assay.regions {
        let (resolved, region) = region_kmers(space, r)?;
        regions.push(resolved);
        kmers.push(region);
    }
    let has = presence(ksets, &kmers.concat(), space).context(KsetsUnreadable)?;
    let calls = has
        .iter()
        .map(|h| {
            let mut rest = h.as_slice();
            let per_region: Vec<ProbeCall> = kmers
                .iter()
                .map(|region| {
                    let (carried, tail) = rest.split_at(region.len());
                    rest = tail;
                    region_call(carried, space.k)
                })
                .collect();
            assay_call(&per_region)
        })
        .collect();
    Ok(ProbeSet {
        assay: assay.id.to_string(),
        regions,
        calls,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    const CDS: &str = "ATGAAACGTCTGGCAGTTCCGATCAGCTAA";

    fn wild(cds: &'static str, strand: &'static str) -> WildType {
        WildType {
            source: "fixture",
            record: "rec",
            cds_start: 101,
            cds_end: 400,
            strand,
            cds,
            residues: "RL",
        }
    }

    fn region(wild_type: Option<WildType>) -> Region {
        Region {
            gene: "gyrA",
            span: Span::Codons(3, 4),
            wild_type,
        }
    }

    fn call(space: &KmerSpace, probes: &[u64], genome: &str) -> ProbeCall {
        let mut carried = HashSet::new();
        space.for_each_canonical(genome.as_bytes(), |x| {
            carried.insert(space.mix(x));
        });
        let carried: Vec<bool> = probes.iter().map(|x| carried.contains(x)).collect();
        region_call(&carried, space.k)
    }

    #[test]
    fn probed_codons_map_to_record_coordinates_on_either_strand() {
        let space = KmerSpace::new(7);
        let (r, kmers) = region_kmers(&space, &region(Some(wild(CDS, "+")))).unwrap();
        assert_eq!((r.start, r.end, r.sequence.as_str()), (107, 112, "CGTCTG"));
        assert_eq!((r.probe_kmers, kmers.len()), (12, 12));
        let (r, _) = region_kmers(&space, &region(Some(wild(CDS, "-")))).unwrap();
        assert_eq!((r.start, r.end), (389, 394));
    }

    #[test]
    fn sources_must_encode_their_wild_type_residues() {
        let space = KmerSpace::new(7);
        let mutant = wild("ATGAAACTTCTGGCAGTTCCGATCAGCTAA", "+");
        let err = region_kmers(&space, &region(Some(mutant))).unwrap_err();
        assert!(format!("{err:#}").contains("encodes LL"), "{err:#}");
        assert!(region_kmers(&space, &region(None)).is_err());
    }

    #[test]
    fn sources_must_extend_a_k_mer_past_both_ends_of_the_probed_codons() {
        let space = KmerSpace::new(7);
        let upstream = Region {
            span: Span::Codons(2, 3),
            wild_type: Some(WildType {
                residues: "KR",
                ..wild(CDS, "+")
            }),
            ..region(None)
        };
        let downstream = region(Some(wild("ATGAAACGTCTGGCA", "+")));
        for short in [upstream, downstream] {
            let err = region_kmers(&space, &short).unwrap_err();
            assert!(
                format!("{err:#}").contains("does not extend 6 bases"),
                "{err:#}"
            );
        }
    }

    #[test]
    fn calls_follow_the_probed_bases_and_ignore_flanking_changes() {
        let space = KmerSpace::new(7);
        let (_, probes) = region_kmers(&space, &region(Some(wild(CDS, "+")))).unwrap();
        let genome = format!("GATTACAGCT{CDS}CTTCAGGTCA");
        assert_eq!(call(&space, &probes, &genome), ProbeCall::WildType);
        let cds_at = 10;
        let (lo, hi) = (cds_at + 6, cds_at + 12);
        for at in lo - 6..hi + 6 {
            for base in *b"ACGT" {
                let mut changed = genome.clone().into_bytes();
                if changed[at] == base {
                    continue;
                }
                changed[at] = base;
                let got = call(&space, &probes, std::str::from_utf8(&changed).unwrap());
                let want = if (lo..hi).contains(&at) {
                    ProbeCall::MutationDetected
                } else {
                    ProbeCall::WildType
                };
                assert_eq!(got, want, "base {at} changed to {}", base as char);
            }
        }
        assert_eq!(
            call(&space, &probes, "CCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGG"),
            ProbeCall::NoCall
        );
    }

    #[test]
    fn a_mutated_region_decides_the_assay_before_a_region_without_a_call() {
        let (wild, mutated, none) = (
            ProbeCall::WildType,
            ProbeCall::MutationDetected,
            ProbeCall::NoCall,
        );
        assert_eq!(assay_call(&[wild, wild]), wild);
        assert_eq!(assay_call(&[wild, none]), none);
        assert_eq!(assay_call(&[none, mutated]), mutated);
        assert_eq!(assay_call(&[mutated, wild]), mutated);
    }

    #[test]
    fn verified_sources_resolve_on_their_records() {
        let want = [
            ("gc_gyra", "gyrA", 479228, 479242, "TCCGCAGTTTACGAC"),
            ("ecoli_qrdr", "gyrA", 2608196, 2608210, "TCGGCGGTTTATGAC"),
            ("ecoli_qrdr", "parC", 3165722, 3165736, "AGCGCCTGTTATGAA"),
            ("salmonella_qrdr", "gyrA", 382517, 382531, "TCCGCAGTGTATGAC"),
        ];
        for k in [31, 32] {
            let space = KmerSpace::new(k);
            let mut got = Vec::new();
            for assay in ASSAYS {
                for r in assay.regions {
                    match region_kmers(&space, r) {
                        Ok((resolved, kmers)) => {
                            assert_eq!(kmers.len(), 15 + k - 1);
                            got.push((
                                assay.id,
                                r.gene,
                                resolved.start,
                                resolved.end,
                                resolved.sequence,
                            ));
                        }
                        Err(e) => assert_eq!(assay.species, "Mycobacterium tuberculosis", "{e:#}"),
                    }
                }
            }
            let got: Vec<_> = got
                .iter()
                .map(|(a, g, s, e, q)| (*a, *g, *s, *e, q.as_str()))
                .collect();
            assert_eq!(got, want);
        }
    }

    #[test]
    fn only_ecoli_qrdr_is_unsupported_and_it_reads_no_kset() {
        let space = KmerSpace::new(31);
        let missing = [std::path::PathBuf::from("/nonexistent/g.kset")];
        let sourced = ASSAYS
            .iter()
            .filter(|a| a.regions.iter().all(|r| r.wild_type.is_some()));
        for assay in sourced {
            let err = build_probes(assay, &missing, &space).unwrap_err();
            match assay.unsupported {
                Some(reason) => {
                    assert!(!err.is::<KsetsUnreadable>(), "{err:#}");
                    assert_eq!(format!("{err:#}"), reason);
                }
                None => assert!(err.is::<KsetsUnreadable>(), "{}: {err:#}", assay.id),
            }
        }
        let unsupported: Vec<&str> = ASSAYS
            .iter()
            .filter(|a| a.unsupported.is_some())
            .map(|a| a.id)
            .collect();
        assert_eq!(unsupported, ["ecoli_qrdr"]);
    }
}
