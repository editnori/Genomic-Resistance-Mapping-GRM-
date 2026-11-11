//! Simulated targeted molecular tests.
//!
//! A PCR panel reports a target only for the organisms its maker lists for
//! that target, and its maker reads each target against a stated drug class.
//! Here a target is a catalog gene group detected by k-mer containment, and a
//! panel call for a drug is positive when any target interpreted for that
//! drug, and reported for the dataset's species, is detected. Each panel
//! answers: if a lab ran an assay with these targets on these isolates, how
//! much of the phenotypic resistance would the panel report, and how often
//! would it report a susceptible isolate as carrying a marker?
//!
//! Target lists, drug readings and organism lists follow the manufacturers'
//! published documents, cited on each list below. The blood-culture panel
//! also reads its carbapenemase targets for the other drugs a primary source
//! says the enzyme hydrolyses, cited on those lists. Detection here is in silico,
//! so it measures what the targets can explain, not how a given cartridge
//! performs, and a detected marker is what the panel would report, not proof
//! that the gene causes the resistance.

use crate::catalog::FamilyHit;
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct Target {
    /// A catalog gene family or group; a family entry matches every group in
    /// it.
    pub gene: &'static str,
    /// The drugs the manufacturer reads this marker against.
    pub drugs: &'static [&'static str],
    /// The organisms the assay reports this marker for, each a genus or a
    /// species name as datasets spell it.
    pub organisms: &'static [&'static str],
}

#[derive(Clone, Debug, Serialize)]
pub struct Panel {
    pub id: &'static str,
    pub name: &'static str,
    /// The kind of assay this mirrors.
    pub mirrors: &'static str,
    pub targets: &'static [Target],
}

/// Why a panel with a target interpreted for a drug reports no result on a
/// dataset of that drug.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Unsupported {
    /// The assay reports none of the panel's targets for the species.
    Organism,
    /// The assay reports some of the panel's targets for the species, but
    /// none of those is interpreted for the drug.
    Drug,
}

impl Unsupported {
    pub fn reason(self) -> &'static str {
        match self {
            Unsupported::Organism => {
                "the panel does not cover this organism: the assay reports none of its targets for it"
            }
            Unsupported::Drug => {
                "no target of this panel is interpreted for this drug in this organism"
            }
        }
    }
}

const CARBAPENEMS: &[&str] = &["meropenem", "imipenem", "ertapenem", "doripenem"];
const EXTENDED_CEPHALOSPORINS: &[&str] = &[
    "ceftriaxone",
    "cefotaxime",
    "ceftazidime",
    "cefepime",
    "aztreonam",
    "cefpodoxime",
];
/// KPC, Bush-Jacoby group 2f: "Carbapenems, oxyimino-β-lactams,
/// cephamycins" (Sawa et al., J Intensive Care 2020;8:13, Table 1, after
/// Bush and Jacoby 2010); KPC-1 "hydrolyzed not only carbapenems but also
/// penicillins, cephalosporins, and monobactams" (Hammoudi Halat and Ayoub
/// Moubareck, Antibiotics 2020;9:186, section 3). So the carbapenems, the
/// extended-spectrum cephalosporins and aztreonam.
const KPC_DRUGS: &[&str] = &[
    "meropenem",
    "imipenem",
    "ertapenem",
    "doripenem",
    "ceftriaxone",
    "cefotaxime",
    "ceftazidime",
    "cefepime",
    "aztreonam",
    "cefpodoxime",
];
/// NDM, VIM and IMP, metallo-β-lactamases: they "hydrolyze all β-lactams
/// except aztreonam" (Nordmann et al., Emerg Infect Dis 2011;17:1791, p.
/// 1792, of VIM and IMP); NDM "confers resistance to all current β-lactam
/// antibiotics, with the exception of aztreonam" and IMP reaches "oxymino
/// cephalosporins" (BCID2 IFU pp. 9 and 8); MBLs hydrolyse "virtually all
/// β-lactam antibiotics including carbapenems with the exception of
/// monobactams" (Hammoudi Halat and Ayoub Moubareck 2020, section 4). So the
/// carbapenems and the extended-spectrum cephalosporins, not aztreonam.
const MBL_DRUGS: &[&str] = &[
    "meropenem",
    "imipenem",
    "ertapenem",
    "doripenem",
    "ceftriaxone",
    "cefotaxime",
    "ceftazidime",
    "cefepime",
    "cefpodoxime",
];
const POLYMYXINS: &[&str] = &["colistin", "polymyxin b"];
const METHICILLIN: &[&str] = &["oxacillin", "methicillin", "cefoxitin"];
const GLYCOPEPTIDES: &[&str] = &["vancomycin", "teicoplanin"];

/// Genera of the Enterobacterales that BCID2 reports as Enterobacterales,
/// from its on-panel list (BCID2 IFU RFIT-PRT-0841-05, Table 130, p. 102).
/// Shigella reads as E. coli there (footnote f).
const ENTEROBACTERALES: &[&str] = &[
    "Citrobacter",
    "Enterobacter",
    "Escherichia",
    "Klebsiella",
    "Morganella",
    "Proteus",
    "Providencia",
    "Salmonella",
    "Serratia",
    "Shigella",
];

/// BCID2's applicable bacteria for CTX-M, IMP, KPC, NDM and VIM: the
/// Enterobacterales, the A. calcoaceticus-baumannii complex and
/// P. aeruginosa (IFU Table 6, p. 29; complex members, Table 130, p. 102).
const BCID2_GRAM_NEGATIVES: &[&str] = &[
    "Citrobacter",
    "Enterobacter",
    "Escherichia",
    "Klebsiella",
    "Morganella",
    "Proteus",
    "Providencia",
    "Salmonella",
    "Serratia",
    "Shigella",
    "Acinetobacter baumannii",
    "Acinetobacter calcoaceticus",
    "Acinetobacter nosocomialis",
    "Acinetobacter pittii",
    "Acinetobacter seifertii",
    "Pseudomonas aeruginosa",
];

/// BCID2's applicable bacteria for mcr-1 (IFU Table 6, p. 28): E. cloacae
/// complex, E. coli, K. aerogenes, K. oxytoca, K. pneumoniae group and
/// Salmonella, with the species the IFU says each result covers (pp. 25,
/// 102). Shigella is reported as E. coli (Table 130, footnote f).
const MCR1_HOSTS: &[&str] = &[
    "Enterobacter asburiae",
    "Enterobacter cloacae",
    "Enterobacter hormaechei",
    "Enterobacter kobei",
    "Enterobacter ludwigii",
    "Enterobacter mori",
    "Escherichia coli",
    "Shigella",
    "Klebsiella aerogenes",
    "Klebsiella oxytoca",
    "Klebsiella grimontii",
    "Klebsiella michiganensis",
    "Klebsiella pneumoniae",
    "Klebsiella quasipneumoniae",
    "Klebsiella variicola",
    "Salmonella",
];

/// Xpert Carba-R colonies: "Enterobacteriaceae, Acinetobacter baumannii, or
/// Pseudomonas aeruginosa" (FDA K160901, section H).
const CARBA_R_HOSTS: &[&str] = &[
    "Citrobacter",
    "Enterobacter",
    "Escherichia",
    "Klebsiella",
    "Morganella",
    "Proteus",
    "Providencia",
    "Salmonella",
    "Serratia",
    "Shigella",
    "Acinetobacter baumannii",
    "Pseudomonas aeruginosa",
];

/// MRSA assays report S. aureus (FDA K162444 and K130894, section H); BCID2
/// also reports mecA/C with S. epidermidis and S. lugdunensis (IFU Table 6).
const STAPHYLOCOCCI: &[&str] = &[
    "Staphylococcus aureus",
    "Staphylococcus epidermidis",
    "Staphylococcus lugdunensis",
];

/// BCID2 reports vanA/B with E. faecalis and E. faecium (IFU Table 6); the
/// Xpert assays name vancomycin-resistant enterococci.
const ENTEROCOCCI: &[&str] = &["Enterococcus faecalis", "Enterococcus faecium"];

const fn target(
    gene: &'static str,
    drugs: &'static [&'static str],
    organisms: &'static [&'static str],
) -> Target {
    Target {
        gene,
        drugs,
        organisms,
    }
}

pub const PANELS: &[Panel] = &[
    // FDA K160901, section H: the five genes are "associated with
    // carbapenem-non-susceptibility".
    Panel {
        id: "carbapenemase",
        name: "Carbapenemase PCR",
        mirrors: "Cepheid Xpert Carba-R targets: KPC, NDM, VIM, IMP, OXA-48-like",
        targets: &[
            target("blaKPC", CARBAPENEMS, CARBA_R_HOSTS),
            target("blaNDM", CARBAPENEMS, CARBA_R_HOSTS),
            target("blaVIM", CARBAPENEMS, CARBA_R_HOSTS),
            target("blaIMP", CARBAPENEMS, CARBA_R_HOSTS),
            target("blaOXA (OXA-48 family)", CARBAPENEMS, CARBA_R_HOSTS),
        ],
    },
    // BCID2 IFU: CTX-M is read against ESBL activity (Table 78), the five
    // carbapenemases against carbapenem resistance (Table 79), mcr-1 against
    // colistin (Table 82); OXA-48-like is reported for Enterobacterales only
    // (Table 6). This panel asks what its targets explain, so KPC, NDM, VIM
    // and IMP are also read for the drugs the enzyme hydrolyses (KPC_DRUGS,
    // MBL_DRUGS). OXA-48-like stays on the carbapenems: "extremely weak
    // activity against expanded-spectrum cephalosporins" (IFU p. 9). The
    // `ctx_m` panel keeps the maker's CTX-M-only cephalosporin reading.
    Panel {
        id: "gram_negative_bcid",
        name: "Gram-negative blood-culture panel",
        mirrors: "BioFire BCID2 resistance targets: CTX-M, KPC, NDM, VIM, IMP, OXA-48-like, mcr-1",
        targets: &[
            target("blaCTX-M", EXTENDED_CEPHALOSPORINS, BCID2_GRAM_NEGATIVES),
            target("blaKPC", KPC_DRUGS, BCID2_GRAM_NEGATIVES),
            target("blaNDM", MBL_DRUGS, BCID2_GRAM_NEGATIVES),
            target("blaVIM", MBL_DRUGS, BCID2_GRAM_NEGATIVES),
            target("blaIMP", MBL_DRUGS, BCID2_GRAM_NEGATIVES),
            target("blaOXA (OXA-48 family)", CARBAPENEMS, ENTEROBACTERALES),
            target("mcr-1", POLYMYXINS, MCR1_HOSTS),
        ],
    },
    Panel {
        id: "ctx_m",
        name: "CTX-M PCR",
        mirrors: "ESBL screening by CTX-M detection (BCID2 and similar)",
        targets: &[target(
            "blaCTX-M",
            EXTENDED_CEPHALOSPORINS,
            BCID2_GRAM_NEGATIVES,
        )],
    },
    // BCID2 IFU Tables 80 and 81: methicillin (oxacillin/cefoxitin).
    Panel {
        id: "mec",
        name: "mecA/mecC PCR",
        mirrors: "MRSA molecular assays (Xpert MRSA, BCID2 mecA/C)",
        targets: &[
            target("mecA", METHICILLIN, STAPHYLOCOCCI),
            target("mecC", METHICILLIN, STAPHYLOCOCCI),
        ],
    },
    // Xpert vanA/vanB insert 301-0188 Rev. F, section 5: "genes for
    // vancomycin/teicoplanin resistance (vanA, vanB)".
    Panel {
        id: "van",
        name: "vanA/vanB PCR",
        mirrors: "VRE molecular assays (Xpert vanA/vanB, BCID2 vanA/B)",
        targets: &[
            target("vanA", GLYCOPEPTIDES, ENTEROCOCCI),
            target("vanB", GLYCOPEPTIDES, ENTEROCOCCI),
        ],
    },
    // No manufacturer document names this assay; the organisms are BCID2's
    // mcr-1 hosts, the one published organism list for an mcr target.
    Panel {
        id: "mcr",
        name: "mcr PCR",
        mirrors: "Plasmid colistin resistance (mcr-1 to mcr-10) assays",
        targets: &[
            target("mcr-1", POLYMYXINS, MCR1_HOSTS),
            target("mcr-2", POLYMYXINS, MCR1_HOSTS),
            target("mcr-3", POLYMYXINS, MCR1_HOSTS),
            target("mcr-4", POLYMYXINS, MCR1_HOSTS),
            target("mcr-5", POLYMYXINS, MCR1_HOSTS),
            target("mcr-6", POLYMYXINS, MCR1_HOSTS),
            target("mcr-7", POLYMYXINS, MCR1_HOSTS),
            target("mcr-8", POLYMYXINS, MCR1_HOSTS),
            target("mcr-9", POLYMYXINS, MCR1_HOSTS),
            target("mcr-10", POLYMYXINS, MCR1_HOSTS),
        ],
    },
];

/// Whether `species` is, or belongs to a genus, in `organisms`.
pub fn covers(organisms: &[&str], species: &str) -> bool {
    organisms.iter().any(|o| {
        species
            .strip_prefix(o)
            .is_some_and(|rest| rest.is_empty() || rest.starts_with(' '))
    })
}

impl Panel {
    /// The targets this panel would report for `drug` in `species`, or why
    /// it would report none.
    pub fn targets_for(&self, drug: &str, species: &str) -> Result<Vec<&Target>, Unsupported> {
        let reported: Vec<&Target> = self
            .targets
            .iter()
            .filter(|t| covers(t.organisms, species))
            .collect();
        if reported.is_empty() {
            return Err(Unsupported::Organism);
        }
        let used: Vec<&Target> = reported
            .into_iter()
            .filter(|t| t.drugs.contains(&drug))
            .collect();
        if used.is_empty() {
            return Err(Unsupported::Drug);
        }
        Ok(used)
    }
}

pub fn target_hit(target: &str, hit: &FamilyHit) -> bool {
    hit.group == target || hit.family == target
}

/// Targets detected among a genome's catalog hits, as the groups that
/// matched.
pub fn detected(targets: &[&Target], hits: &[FamilyHit]) -> Vec<String> {
    let mut out: Vec<String> = hits
        .iter()
        .filter(|h| targets.iter().any(|t| target_hit(t.gene, h)))
        .map(|h| h.group.clone())
        .collect();
    out.sort();
    out.dedup();
    out
}

/// A panel as it would run on one dataset, or why it would not.
pub type Applied = (&'static Panel, Result<Vec<&'static Target>, Unsupported>);

/// Every panel with a target interpreted for `drug`, each with the targets
/// it would report on a `species` dataset or the reason it reports none.
pub fn for_dataset(drug: &str, species: &str) -> Vec<Applied> {
    PANELS
        .iter()
        .filter(|p| p.targets.iter().any(|t| t.drugs.contains(&drug)))
        .map(|p| (p, p.targets_for(drug, species)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hit(family: &str, group: &str) -> FamilyHit {
        FamilyHit {
            family: family.into(),
            group: group.into(),
            allele: "x".into(),
            containment: 1.0,
            class: "BETA-LACTAM".into(),
            subclass: "CARBAPENEM".into(),
        }
    }

    fn panel(id: &str) -> &'static Panel {
        PANELS.iter().find(|p| p.id == id).unwrap()
    }

    fn genes(panel: &str, drug: &str, species: &str) -> Result<Vec<&'static str>, Unsupported> {
        panel_targets(panel, drug, species).map(|t| t.iter().map(|t| t.gene).collect())
    }

    fn panel_targets(
        id: &str,
        drug: &str,
        species: &str,
    ) -> Result<Vec<&'static Target>, Unsupported> {
        panel(id).targets_for(drug, species)
    }

    #[test]
    fn oxa48_target_does_not_match_other_oxa_groups() {
        let oxa = |group: &str| hit("blaOXA", group);
        let targets = panel_targets("carbapenemase", "meropenem", "Klebsiella pneumoniae").unwrap();
        assert_eq!(
            detected(&targets, &[oxa("blaOXA (OXA-48 family)")]),
            vec!["blaOXA (OXA-48 family)"]
        );
        assert!(detected(&targets, &[oxa("blaOXA (OXA-1 family)")]).is_empty());
        let meropenem = for_dataset("meropenem", "Klebsiella pneumoniae");
        assert_eq!(meropenem.len(), 2);
        assert!(meropenem.iter().all(|(_, t)| t.is_ok()));
    }

    /// Every (panel, drug) pair's targets in a species each target is
    /// reported for, and no drug outside the five classes.
    #[test]
    fn each_panel_and_drug_has_the_documented_targets() {
        let carbapenemases = [
            "blaKPC",
            "blaNDM",
            "blaVIM",
            "blaIMP",
            "blaOXA (OXA-48 family)",
        ];
        let mcr = [
            "mcr-1", "mcr-2", "mcr-3", "mcr-4", "mcr-5", "mcr-6", "mcr-7", "mcr-8", "mcr-9",
            "mcr-10",
        ];
        let mut expected: Vec<(&str, &str, Vec<&str>)> = Vec::new();
        for drug in CARBAPENEMS {
            expected.push(("carbapenemase", drug, carbapenemases.to_vec()));
            expected.push(("gram_negative_bcid", drug, carbapenemases.to_vec()));
        }
        // BCID2 reads KPC for aztreonam too, and the metallo-β-lactamases
        // for the cephalosporins only; the CTX-M panel reads CTX-M alone.
        for drug in EXTENDED_CEPHALOSPORINS {
            let bcid = if *drug == "aztreonam" {
                vec!["blaCTX-M", "blaKPC"]
            } else {
                vec!["blaCTX-M", "blaKPC", "blaNDM", "blaVIM", "blaIMP"]
            };
            expected.push(("gram_negative_bcid", drug, bcid));
            expected.push(("ctx_m", drug, vec!["blaCTX-M"]));
        }
        for drug in POLYMYXINS {
            expected.push(("gram_negative_bcid", drug, vec!["mcr-1"]));
            expected.push(("mcr", drug, mcr.to_vec()));
        }
        for drug in METHICILLIN {
            expected.push(("mec", drug, vec!["mecA", "mecC"]));
        }
        for drug in GLYCOPEPTIDES {
            expected.push(("van", drug, vec!["vanA", "vanB"]));
        }
        let mut pairs = 0;
        for p in PANELS {
            let drugs: std::collections::BTreeSet<&str> = p
                .targets
                .iter()
                .flat_map(|t| t.drugs.iter().copied())
                .collect();
            for drug in drugs {
                let wanted = expected
                    .iter()
                    .find(|(id, d, _)| *id == p.id && *d == drug)
                    .unwrap_or_else(|| panic!("{} is read for {drug}", p.id));
                let got: Vec<&str> = p
                    .targets
                    .iter()
                    .filter(|t| t.drugs.contains(&drug))
                    .map(|t| t.gene)
                    .collect();
                assert_eq!(got, wanted.2, "{} {drug}", p.id);
                pairs += 1;
            }
        }
        assert_eq!(pairs, expected.len());
        // In E. coli every interpreted target is reported.
        assert_eq!(
            genes("gram_negative_bcid", "meropenem", "Escherichia coli").unwrap(),
            carbapenemases
        );
    }

    #[test]
    fn ctx_m_is_not_read_as_a_carbapenem_call() {
        let ctx_m = [hit("blaCTX-M", "blaCTX-M")];
        let bcid = panel("gram_negative_bcid");
        for drug in CARBAPENEMS {
            let targets = bcid.targets_for(drug, "Klebsiella pneumoniae").unwrap();
            assert!(detected(&targets, &ctx_m).is_empty(), "{drug}");
        }
        let ceftriaxone = bcid
            .targets_for("ceftriaxone", "Klebsiella pneumoniae")
            .unwrap();
        assert_eq!(detected(&ceftriaxone, &ctx_m), ["blaCTX-M"]);
        // Nor is mcr-1 read for either.
        let mcr = [hit("mcr-1", "mcr-1")];
        assert!(detected(&ceftriaxone, &mcr).is_empty());
    }

    #[test]
    fn bcid2_reads_a_carbapenemase_for_the_drugs_the_enzyme_hydrolyses() {
        let species = "Klebsiella pneumoniae";
        let call = |id: &str, drug: &str, hits: &[FamilyHit]| {
            detected(&panel_targets(id, drug, species).unwrap(), hits)
        };
        // A KPC genome without CTX-M: BCID2 reports it on a cephalosporin,
        // the CTX-M panel does not.
        let kpc = [hit("blaKPC", "blaKPC")];
        for drug in EXTENDED_CEPHALOSPORINS {
            assert_eq!(call("gram_negative_bcid", drug, &kpc), ["blaKPC"], "{drug}");
            assert!(call("ctx_m", drug, &kpc).is_empty(), "{drug}");
        }
        // A metallo-β-lactamase is read for ceftazidime, not for aztreonam.
        for mbl in ["blaNDM", "blaVIM", "blaIMP"] {
            let hits = [hit(mbl, mbl)];
            assert_eq!(call("gram_negative_bcid", "ceftazidime", &hits), [mbl]);
            assert!(call("gram_negative_bcid", "aztreonam", &hits).is_empty());
        }
        // OXA-48-like and mcr-1 are no cephalosporin call, and Carba-R reads
        // no cephalosporin.
        let others = [
            hit("blaOXA", "blaOXA (OXA-48 family)"),
            hit("mcr-1", "mcr-1"),
        ];
        assert!(call("gram_negative_bcid", "cefotaxime", &others).is_empty());
        assert!(
            for_dataset("ceftriaxone", species)
                .iter()
                .all(|(p, _)| p.id != "carbapenemase")
        );
    }

    #[test]
    fn a_panel_reports_only_for_the_organisms_its_assay_names() {
        let organism = Err(Unsupported::Organism);
        assert_eq!(genes("mec", "cefoxitin", "Escherichia coli"), organism);
        assert_eq!(genes("mec", "cefoxitin", "Salmonella enterica"), organism);
        assert_eq!(
            genes("mec", "cefoxitin", "Staphylococcus aureus").unwrap(),
            ["mecA", "mecC"]
        );
        assert_eq!(
            genes("mec", "oxacillin", "Staphylococcus haemolyticus"),
            organism
        );
        assert_eq!(
            genes("van", "vancomycin", "Staphylococcus aureus"),
            organism
        );
        assert!(genes("van", "vancomycin", "Enterococcus faecium").is_ok());
        assert_eq!(
            genes("carbapenemase", "meropenem", "Streptococcus pneumoniae"),
            organism
        );
        assert_eq!(
            genes("gram_negative_bcid", "ceftriaxone", "Neisseria gonorrhoeae"),
            organism
        );
        assert_eq!(
            genes("ctx_m", "ceftriaxone", "Clostridioides difficile"),
            organism
        );
        // OXA-48-like is reported with the Enterobacterales only.
        assert_eq!(
            genes("gram_negative_bcid", "imipenem", "Acinetobacter baumannii").unwrap(),
            ["blaKPC", "blaNDM", "blaVIM", "blaIMP"]
        );
        assert_eq!(
            genes("carbapenemase", "imipenem", "Acinetobacter baumannii")
                .unwrap()
                .len(),
            5
        );
        // A genus entry covers its species, and only whole names match.
        assert!(covers(ENTEROBACTERALES, "Shigella sonnei"));
        assert!(covers(MCR1_HOSTS, "Shigella sonnei"));
        assert!(covers(MCR1_HOSTS, "Klebsiella quasipneumoniae"));
        assert!(!covers(MCR1_HOSTS, "Klebsiella"));
        assert!(!covers(
            &["Klebsiella pneumoniae"],
            "Klebsiella pneumoniaeX"
        ));
        assert!(!covers(ENTEROBACTERALES, "Pseudomonas aeruginosa"));
    }

    #[test]
    fn each_reason_is_recorded_apart() {
        // Covered organism, but its only polymyxin target is not reported
        // there.
        assert_eq!(
            genes("gram_negative_bcid", "colistin", "Pseudomonas aeruginosa"),
            Err(Unsupported::Drug)
        );
        assert_eq!(
            genes("gram_negative_bcid", "colistin", "Streptococcus pneumoniae"),
            Err(Unsupported::Organism)
        );
        assert_ne!(Unsupported::Organism.reason(), Unsupported::Drug.reason());
        let mec = for_dataset("cefoxitin", "Escherichia coli");
        assert_eq!(mec.len(), 1);
        assert_eq!(
            (mec[0].0.id, mec[0].1.as_ref().err()),
            ("mec", Some(&Unsupported::Organism))
        );
        // A panel with no target read for the drug is not listed at all.
        assert!(for_dataset("ciprofloxacin", "Escherichia coli").is_empty());
    }

    #[test]
    fn host_lists_extend_the_enterobacterales() {
        assert_eq!(
            &BCID2_GRAM_NEGATIVES[..ENTEROBACTERALES.len()],
            ENTEROBACTERALES
        );
        assert_eq!(&CARBA_R_HOSTS[..ENTEROBACTERALES.len()], ENTEROBACTERALES);
    }
}
