//! Simulated targeted molecular tests.
//!
//! A PCR panel calls resistance when any of its targets is detected. Here a
//! target is a catalog gene group detected by k-mer containment, so each
//! panel answers: if a lab ran an assay with these targets on these isolates,
//! how much of the phenotypic resistance would it catch, and how often would
//! it call a susceptible isolate resistant?
//!
//! Target lists follow the manufacturers' published target lists; detection
//! here is in silico, so it measures what the targets can explain, not how a
//! given cartridge performs.

use crate::catalog::FamilyHit;
use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct Panel {
    pub id: &'static str,
    pub name: &'static str,
    /// The kind of assay this mirrors.
    pub mirrors: &'static str,
    /// Gene families or groups; a family entry matches every group in it.
    pub targets: &'static [&'static str],
    pub drugs: &'static [&'static str],
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

pub const PANELS: &[Panel] = &[
    Panel {
        id: "carbapenemase",
        name: "Carbapenemase PCR",
        mirrors: "Cepheid Xpert Carba-R targets: KPC, NDM, VIM, IMP, OXA-48-like",
        targets: &[
            "blaKPC",
            "blaNDM",
            "blaVIM",
            "blaIMP",
            "blaOXA (OXA-48 family)",
        ],
        drugs: CARBAPENEMS,
    },
    Panel {
        id: "gram_negative_bcid",
        name: "Gram-negative blood-culture panel",
        mirrors: "BioFire BCID2 resistance targets: CTX-M, KPC, NDM, VIM, IMP, OXA-48-like, mcr-1",
        targets: &[
            "blaCTX-M",
            "blaKPC",
            "blaNDM",
            "blaVIM",
            "blaIMP",
            "blaOXA (OXA-48 family)",
            "mcr-1",
        ],
        drugs: &[
            "ceftriaxone",
            "cefotaxime",
            "ceftazidime",
            "cefepime",
            "aztreonam",
            "cefpodoxime",
            "meropenem",
            "imipenem",
            "ertapenem",
            "doripenem",
        ],
    },
    Panel {
        id: "ctx_m",
        name: "CTX-M PCR",
        mirrors: "ESBL screening by CTX-M detection (BCID2 and similar)",
        targets: &["blaCTX-M"],
        drugs: EXTENDED_CEPHALOSPORINS,
    },
    Panel {
        id: "mec",
        name: "mecA/mecC PCR",
        mirrors: "MRSA molecular assays (Xpert MRSA, BCID2 mecA/C)",
        targets: &["mecA", "mecC"],
        drugs: &["oxacillin", "methicillin", "cefoxitin"],
    },
    Panel {
        id: "van",
        name: "vanA/vanB PCR",
        mirrors: "VRE molecular assays (Xpert vanA/vanB, BCID2 vanA/B)",
        targets: &["vanA", "vanB"],
        drugs: &["vancomycin", "teicoplanin"],
    },
    Panel {
        id: "mcr",
        name: "mcr PCR",
        mirrors: "Plasmid colistin resistance (mcr-1 to mcr-10) assays",
        targets: &[
            "mcr-1", "mcr-2", "mcr-3", "mcr-4", "mcr-5", "mcr-6", "mcr-7", "mcr-8", "mcr-9",
            "mcr-10",
        ],
        drugs: &["colistin", "polymyxin b"],
    },
];

pub fn target_hit(target: &str, hit: &FamilyHit) -> bool {
    hit.group == target || hit.family == target
}

/// Panel targets detected among a genome's catalog hits.
pub fn detected(panel: &Panel, hits: &[FamilyHit]) -> Vec<String> {
    let mut out: Vec<String> = hits
        .iter()
        .filter(|h| panel.targets.iter().any(|t| target_hit(t, h)))
        .map(|h| h.group.clone())
        .collect();
    out.sort();
    out.dedup();
    out
}

pub fn for_drug(drug: &str) -> Vec<&'static Panel> {
    PANELS.iter().filter(|p| p.drugs.contains(&drug)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn oxa48_target_does_not_match_other_oxa_groups() {
        let hit = |group: &str| FamilyHit {
            family: "blaOXA".into(),
            group: group.into(),
            allele: "x".into(),
            containment: 1.0,
            class: "BETA-LACTAM".into(),
            subclass: "CARBAPENEM".into(),
        };
        let carba = &PANELS[0];
        assert_eq!(
            detected(carba, &[hit("blaOXA (OXA-48 family)")]),
            vec!["blaOXA (OXA-48 family)"]
        );
        assert!(detected(carba, &[hit("blaOXA (OXA-1 family)")]).is_empty());
        assert_eq!(for_drug("meropenem").len(), 2);
    }
}
