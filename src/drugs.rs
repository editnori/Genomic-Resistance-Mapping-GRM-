//! Antibiotic names: resolution of raw BV-BRC strings, and drug classes.
//!
//! BV-BRC antibiotic strings come from many submitters: Danish spellings
//! (tigecyklin), typos (trimotheprim), brands (Synercid, Timentin),
//! underscores, swapped combination order, class labels (carbapenem), stray
//! encoding bytes, and a few non-drugs. Each raw string resolves through:
//!
//! 1. SigDrug's reviewed vocabulary (exact ingredient names), its typo,
//!    brand-to-generic, and abbreviation folds, and its salt-to-base table,
//!    when a SigDrug vocabulary directory is supplied;
//! 2. GRM's antibiotic vocabulary (the class table below) and its own
//!    overrides for names SigDrug does not carry (veterinary and TB drugs,
//!    international spellings);
//! 3. edit distance of at most two to a GRM antibiotic name, unique best.
//!
//! Combinations resolve component by component; the canonical name lists the
//! principal drugs alphabetically and beta-lactamase inhibitors last. Every
//! resolution records the step that produced it.

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::path::Path;

const INHIBITORS: &[&str] = &[
    "clavulanate",
    "sulbactam",
    "tazobactam",
    "avibactam",
    "relebactam",
    "vaborbactam",
    "taniborbactam",
    "durlobactam",
    "enmetazobactam",
];

/// Raw-string fixes SigDrug does not carry. Values may name several
/// components ("a/b").
const OVERRIDES: &[(&str, &str)] = &[
    ("colomycin", "colistin"),
    ("synercid", "dalfopristin/quinupristin"),
    ("timentin", "ticarcillin/clavulanate"),
    ("pristimycin", "pristinamycin"),
    ("phosphomycin", "fosfomycin"),
    ("cefalotin", "cephalothin"),
    ("cefalexin", "cephalexin"),
    ("rifampicin", "rifampin"),
    ("clavulanic acid", "clavulanate"),
    ("clavulanat", "clavulanate"),
    ("sulphadimethoxine", "sulfadimethoxine"),
    ("penicillin g", "penicillin"),
    ("benzylpenicillin", "penicillin"),
];

const CLASS_LABELS: &[&str] = &[
    "beta-lactam",
    "extended spectrum beta lactamase",
    "cephalosporin",
    "carbapenem",
    "fluoroquinolones",
    "sulfonamides",
    "sulfa",
    "aminoglycosides",
    "aminogycosides",
    "macrolides",
];

const NOT_ANTIBIOTICS: &[&str] = &["instrument", "publicatiion", "azidothymidine"];

/// Close to more than one plausible drug: "geamycin" is two edits from
/// neomycin but may be a gentamicin typo; "sulfobactam" is two from
/// sulbactam but sits where co-trimoxazole's sulfamethoxazole belongs.
/// These stay unresolved rather than guessed.
const AMBIGUOUS: &[&str] = &["geamycin", "sulfobactam"];

/// Drug class of a single principal antibiotic. Its keys are GRM's
/// antibiotic vocabulary.
fn single_class(name: &str) -> Option<&'static str> {
    Some(match name {
        "ampicillin" | "amoxicillin" | "penicillin" | "piperacillin" | "ticarcillin"
        | "carbenicillin" | "oxacillin" | "methicillin" | "temocillin" | "nafcillin"
        | "dicloxacillin" | "cloxacillin" => "penicillin",
        "cefazolin" | "cephalothin" | "cephalexin" | "cefoxitin" | "cefotetan" | "cefmetazole"
        | "cefuroxime" | "cefaclor" => "cephalosporin (1st-2nd gen)",
        "ceftriaxone" | "cefotaxime" | "ceftazidime" | "cefixime" | "cefpodoxime" | "ceftiofur"
        | "cefoperazone" | "cefquinome" | "ceftizoxime" | "cefdinir" | "ceftibuten" => {
            "cephalosporin (3rd gen)"
        }
        "ceftolozane" => "cephalosporin (anti-pseudomonal)",
        "cefepime" | "cefpirome" | "cefozopran" => "cephalosporin (4th gen)",
        "ceftaroline" | "ceftobiprole" => "cephalosporin (anti-MRSA)",
        "cefiderocol" => "siderophore cephalosporin",
        "meropenem" | "imipenem" | "ertapenem" | "doripenem" | "biapenem" => "carbapenem",
        "aztreonam" => "monobactam",
        "gentamicin" | "amikacin" | "tobramycin" | "kanamycin" | "streptomycin" | "netilmicin"
        | "neomycin" | "spectinomycin" | "capreomycin" | "apramycin" | "plazomicin"
        | "arbekacin" => "aminoglycoside",
        "ciprofloxacin" | "levofloxacin" | "moxifloxacin" | "ofloxacin" | "norfloxacin"
        | "nalidixic acid" | "enrofloxacin" | "marbofloxacin" | "pefloxacin" | "gatifloxacin"
        | "delafloxacin" | "sparfloxacin" | "danofloxacin" | "gemifloxacin" | "trovafloxacin" => {
            "quinolone"
        }
        "azithromycin" | "erythromycin" | "clarithromycin" | "telithromycin" | "tylosin"
        | "tilmicosin" | "spiramycin" | "tulathromycin" => "macrolide",
        "clindamycin" | "lincomycin" | "pirlimycin" => "lincosamide",
        "quinupristin" | "dalfopristin" | "pristinamycin" | "virginiamycin" => "streptogramin",
        "tetracycline" | "doxycycline" | "minocycline" | "oxytetracycline"
        | "chlortetracycline" => "tetracycline",
        "eravacycline" => "fluorocycline",
        "tigecycline" => "glycylcycline",
        "colistin" | "polymyxin b" => "polymyxin",
        "trimethoprim" => "trimethoprim",
        "sulfisoxazole" | "sulfathiazole" | "sulfamethoxazole" | "sulfamethazine"
        | "sulfadimethoxine" => "sulfonamide",
        "chloramphenicol" | "florfenicol" => "phenicol",
        "vancomycin" | "teicoplanin" | "telavancin" | "dalbavancin" | "oritavancin" => {
            "glycopeptide"
        }
        "linezolid" | "tedizolid" => "oxazolidinone",
        "daptomycin" | "surotomycin" => "lipopeptide",
        "fosfomycin" => "fosfomycin",
        "nitrofurantoin" | "furazolidone" => "nitrofuran",
        "metronidazole" => "nitroimidazole",
        "fusidic acid" => "fusidane",
        "mupirocin" => "mupirocin",
        "novobiocin" => "aminocoumarin",
        "tiamulin" => "pleuromutilin",
        "avilamycin" => "orthosomycin",
        "bacitracin" => "bacitracin",
        "rifampin" | "rifabutin" | "rifaximin" => "rifamycin",
        "fidaxomicin" => "macrocycle",
        "chlorhexidine gluconate" => "antiseptic",
        "isoniazid" | "ethionamide" | "prothionamide" => "anti-TB (mycolic acid)",
        "ethambutol" => "anti-TB (arabinogalactan)",
        "pyrazinamide" | "nicotinamide" => "anti-TB (pyrazinamide)",
        "bedaquiline" => "anti-TB (ATP synthase)",
        "delamanid" | "pretomanid" => "anti-TB (nitroimidazole)",
        "clofazimine" => "anti-TB (riminophenazine)",
        "cycloserine" => "anti-TB (cycloserine)",
        "para-aminosalicylic acid" => "anti-TB (folate)",
        _ => return None,
    })
}

const GRM_VOCABULARY: &[&str] = &[
    "ampicillin",
    "amoxicillin",
    "penicillin",
    "piperacillin",
    "ticarcillin",
    "carbenicillin",
    "oxacillin",
    "methicillin",
    "temocillin",
    "nafcillin",
    "dicloxacillin",
    "cloxacillin",
    "cefazolin",
    "cephalothin",
    "cephalexin",
    "cefoxitin",
    "cefotetan",
    "cefmetazole",
    "cefuroxime",
    "cefaclor",
    "ceftriaxone",
    "cefotaxime",
    "ceftazidime",
    "cefixime",
    "cefpodoxime",
    "ceftiofur",
    "cefoperazone",
    "cefquinome",
    "ceftizoxime",
    "cefdinir",
    "ceftibuten",
    "ceftolozane",
    "arbekacin",
    "gemifloxacin",
    "trovafloxacin",
    "chlortetracycline",
    "eravacycline",
    "rifaximin",
    "fidaxomicin",
    "cefepime",
    "cefpirome",
    "cefozopran",
    "ceftaroline",
    "ceftobiprole",
    "cefiderocol",
    "meropenem",
    "imipenem",
    "ertapenem",
    "doripenem",
    "biapenem",
    "aztreonam",
    "gentamicin",
    "amikacin",
    "tobramycin",
    "kanamycin",
    "streptomycin",
    "netilmicin",
    "neomycin",
    "spectinomycin",
    "capreomycin",
    "apramycin",
    "plazomicin",
    "ciprofloxacin",
    "levofloxacin",
    "moxifloxacin",
    "ofloxacin",
    "norfloxacin",
    "nalidixic acid",
    "enrofloxacin",
    "marbofloxacin",
    "pefloxacin",
    "gatifloxacin",
    "delafloxacin",
    "sparfloxacin",
    "danofloxacin",
    "azithromycin",
    "erythromycin",
    "clarithromycin",
    "telithromycin",
    "tylosin",
    "tilmicosin",
    "spiramycin",
    "tulathromycin",
    "clindamycin",
    "lincomycin",
    "pirlimycin",
    "quinupristin",
    "dalfopristin",
    "pristinamycin",
    "virginiamycin",
    "tetracycline",
    "doxycycline",
    "minocycline",
    "oxytetracycline",
    "tigecycline",
    "colistin",
    "polymyxin b",
    "trimethoprim",
    "sulfisoxazole",
    "sulfathiazole",
    "sulfamethoxazole",
    "sulfamethazine",
    "sulfadimethoxine",
    "chloramphenicol",
    "florfenicol",
    "vancomycin",
    "teicoplanin",
    "telavancin",
    "dalbavancin",
    "oritavancin",
    "linezolid",
    "tedizolid",
    "daptomycin",
    "surotomycin",
    "fosfomycin",
    "nitrofurantoin",
    "furazolidone",
    "metronidazole",
    "fusidic acid",
    "mupirocin",
    "novobiocin",
    "tiamulin",
    "avilamycin",
    "bacitracin",
    "rifampin",
    "rifabutin",
    "isoniazid",
    "ethionamide",
    "prothionamide",
    "ethambutol",
    "pyrazinamide",
    "nicotinamide",
    "bedaquiline",
    "delamanid",
    "pretomanid",
    "clofazimine",
    "cycloserine",
    "para-aminosalicylic acid",
];

/// SigDrug's reviewed drug vocabulary, read from its
/// `data/drug/compendium/facts/vocabulary` directory.
#[derive(Default)]
pub struct Vocabulary {
    pub source: Option<String>,
    known: HashSet<String>,
    folds: HashMap<String, (String, String)>,
    salts: HashMap<String, String>,
}

fn read_tsv(path: &Path) -> Result<Vec<HashMap<String, String>>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("read {}", path.display()))?;
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    Ok(lines
        .map(|l| {
            header
                .iter()
                .zip(l.split('\t'))
                .map(|(h, v)| (h.to_string(), v.to_string()))
                .collect()
        })
        .collect())
}

impl Vocabulary {
    pub fn load_sigdrug(dir: &Path) -> Result<Vocabulary> {
        let mut v = Vocabulary {
            source: Some(dir.display().to_string()),
            ..Default::default()
        };
        for file in ["preferred_ingredients.tsv", "preferred_base_names.tsv"] {
            for row in read_tsv(&dir.join(file))? {
                if matches!(
                    row.get("classification").map(String::as_str),
                    Some("generic" | "ingredient" | "combination")
                ) {
                    for key in ["value", "preferred_name"] {
                        if let Some(name) = row.get(key) {
                            v.known.insert(name.to_lowercase());
                        }
                    }
                }
            }
        }
        for file in [
            "drug_name_typo_folds.tsv",
            "brand_generic_folds.tsv",
            "ingredient_abbrev_folds.tsv",
        ] {
            for row in read_tsv(&dir.join(file))? {
                if let (Some(o), Some(c)) = (row.get("observed"), row.get("canonical")) {
                    v.folds
                        .entry(o.to_lowercase())
                        .or_insert((c.to_lowercase(), format!("sigdrug:{file}")));
                }
            }
        }
        for row in read_tsv(&dir.join("salt_base_ingredients.tsv"))? {
            if let (Some(i), Some(b)) =
                (row.get("ingredient_name"), row.get("salt_base_ingredient"))
            {
                v.salts.insert(i.to_lowercase(), b.to_lowercase());
            }
        }
        Ok(v)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NameKind {
    Drug,
    /// Names a drug class or test, never a learning target.
    Class,
    /// Not an antibiotic (or not a drug at all).
    Excluded,
    Unresolved,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Component {
    pub text: String,
    pub canonical: String,
    pub source: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Resolution {
    pub raw: String,
    pub canonical: String,
    pub kind: NameKind,
    pub components: Vec<Component>,
}

fn clean(raw: &str) -> String {
    let ascii: String = raw.chars().filter(|c| c.is_ascii()).collect();
    let s = ascii.trim().to_lowercase().replace('_', " ");
    s.split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
        .replace(" / ", "/")
}

/// Optimal-string-alignment edit distance.
fn edit_distance(a: &str, b: &str) -> usize {
    let (a, b): (Vec<char>, Vec<char>) = (a.chars().collect(), b.chars().collect());
    let mut d = vec![vec![0usize; b.len() + 1]; a.len() + 1];
    for (i, row) in d.iter_mut().enumerate() {
        row[0] = i;
    }
    for j in 0..=b.len() {
        d[0][j] = j;
    }
    for i in 1..=a.len() {
        for j in 1..=b.len() {
            let cost = (a[i - 1] != b[j - 1]) as usize;
            d[i][j] = (d[i - 1][j] + 1)
                .min(d[i][j - 1] + 1)
                .min(d[i - 1][j - 1] + cost);
            if i > 1 && j > 1 && a[i - 1] == b[j - 2] && a[i - 2] == b[j - 1] {
                d[i][j] = d[i][j].min(d[i - 2][j - 2] + 1);
            }
        }
    }
    d[a.len()][b.len()]
}

fn grm_known(name: &str) -> bool {
    GRM_VOCABULARY.contains(&name) || INHIBITORS.contains(&name)
}

/// Resolves one component to canonical names (an override can expand to
/// several), or None.
fn resolve_component(text: &str, vocab: &Vocabulary) -> Option<(Vec<String>, String)> {
    let one = |s: &str, src: &str| Some((vec![s.to_string()], src.to_string()));
    let grm_spelling = |s: &str| -> String {
        // SigDrug may carry a salt or alternate spelling GRM files under
        // another name; fold those onto GRM's vocabulary.
        OVERRIDES
            .iter()
            .find(|(o, _)| *o == s)
            .map_or(s.to_string(), |(_, c)| c.to_string())
    };
    if AMBIGUOUS.contains(&text) {
        return None;
    }
    if vocab.known.contains(text) {
        return one(&grm_spelling(text), "sigdrug:vocabulary");
    }
    if let Some((c, src)) = vocab.folds.get(text) {
        // SigDrug writes combinations as "a-b"; split only when every part
        // is itself a known drug, so hyphenated single names stay whole.
        let canon = grm_spelling(c);
        let parts: Vec<String> = canon.split(['/', '-']).map(String::from).collect();
        let split = parts.len() > 1
            && parts
                .iter()
                .all(|p| vocab.known.contains(p) || grm_known(p));
        return Some((
            if split {
                parts.into_iter().map(|p| grm_spelling(&p)).collect()
            } else {
                vec![canon]
            },
            src.clone(),
        ));
    }
    if let Some(base) = vocab.salts.get(text) {
        return one(&grm_spelling(base), "sigdrug:salt_base_ingredients.tsv");
    }
    if grm_known(text) {
        return one(text, "grm:vocabulary");
    }
    if let Some((_, c)) = OVERRIDES.iter().find(|(o, _)| *o == text) {
        return Some((
            c.split('/').map(String::from).collect(),
            "grm:override".to_string(),
        ));
    }
    if text.len() >= 7 {
        let mut scored: Vec<(usize, &str)> = GRM_VOCABULARY
            .iter()
            .chain(INHIBITORS)
            .map(|&k| (edit_distance(text, k), k))
            .filter(|(d, _)| *d <= 2)
            .collect();
        scored.sort();
        if let Some(&(d, best)) = scored.first() {
            if scored.get(1).is_none_or(|(d2, _)| *d2 > d) {
                return one(best, &format!("grm:edit-distance-{d}"));
            }
        }
    }
    None
}

pub fn resolve(raw: &str, vocab: &Vocabulary) -> Resolution {
    let text = clean(raw);
    let make = |canonical: String, kind, components| Resolution {
        raw: raw.to_string(),
        canonical,
        kind,
        components,
    };
    if CLASS_LABELS.contains(&text.as_str()) {
        return make(text, NameKind::Class, vec![]);
    }
    if NOT_ANTIBIOTICS.contains(&text.as_str()) {
        return make(text, NameKind::Excluded, vec![]);
    }
    // Whole string first (names like "nalidixic acid" contain spaces), then
    // '/' components, then a two-way split at a space for "a b" combos.
    let mut parts: Vec<String> = text.split('/').map(|p| p.trim().to_string()).collect();
    if parts.len() == 1 && resolve_component(&text, vocab).is_none() {
        let words: Vec<&str> = text.split(' ').collect();
        for cut in 1..words.len() {
            let (a, b) = (words[..cut].join(" "), words[cut..].join(" "));
            if resolve_component(&a, vocab).is_some() && resolve_component(&b, vocab).is_some() {
                parts = vec![a, b];
                break;
            }
        }
    }
    let mut components = Vec::new();
    let mut names = Vec::new();
    let mut unresolved = false;
    for p in &parts {
        match resolve_component(p, vocab) {
            Some((canon, source)) => {
                for c in canon {
                    names.push(c.clone());
                    components.push(Component {
                        text: p.clone(),
                        canonical: c,
                        source: source.clone(),
                    });
                }
            }
            None => {
                unresolved = true;
                names.push(p.clone());
                components.push(Component {
                    text: p.clone(),
                    canonical: p.clone(),
                    source: "unresolved".into(),
                });
            }
        }
    }
    let (mut principal, mut inhibitors): (Vec<String>, Vec<String>) = names
        .into_iter()
        .partition(|n| !INHIBITORS.contains(&n.as_str()));
    principal.sort();
    principal.dedup();
    inhibitors.sort();
    inhibitors.dedup();
    let canonical = principal
        .into_iter()
        .chain(inhibitors)
        .collect::<Vec<_>>()
        .join("/");
    let kind = if unresolved {
        NameKind::Unresolved
    } else {
        NameKind::Drug
    };
    make(canonical, kind, components)
}

/// Resolution for every raw string, with a fast lookup by raw string.
pub struct Resolver {
    pub resolutions: HashMap<String, Resolution>,
    vocab: Vocabulary,
}

impl Resolver {
    pub fn new(vocab: Vocabulary) -> Self {
        Resolver {
            resolutions: HashMap::new(),
            vocab,
        }
    }

    pub fn resolve(&mut self, raw: &str) -> &Resolution {
        if !self.resolutions.contains_key(raw) {
            let r = resolve(raw, &self.vocab);
            self.resolutions.insert(raw.to_string(), r);
        }
        &self.resolutions[raw]
    }

    pub fn vocabulary_source(&self) -> Option<&str> {
        self.vocab.source.as_deref()
    }
}

pub fn components(canonical: &str) -> (Vec<&str>, bool) {
    let parts: Vec<&str> = canonical.split('/').collect();
    let has_inhibitor = parts.iter().any(|p| INHIBITORS.contains(p));
    (
        parts
            .into_iter()
            .filter(|p| !INHIBITORS.contains(p))
            .collect(),
        has_inhibitor,
    )
}

/// Drug class of a canonical name, combinations included.
pub fn class_of(canonical: &str) -> String {
    let (principal, inhibited) = components(canonical);
    let classes: Vec<&str> = principal.iter().filter_map(|p| single_class(p)).collect();
    let base = match classes.as_slice() {
        [] => "other".to_string(),
        [one] => one.to_string(),
        many if many
            .iter()
            .all(|c| matches!(*c, "trimethoprim" | "sulfonamide")) =>
        {
            "folate pathway combination".into()
        }
        many if many.iter().all(|c| *c == many[0]) => many[0].to_string(),
        many => many.join(" + "),
    };
    if inhibited {
        let family = if base.starts_with("cephalosporin") {
            "cephalosporin"
        } else {
            base.as_str()
        };
        format!("{family} + inhibitor")
    } else {
        base
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vocab() -> Vocabulary {
        let mut v = Vocabulary::default();
        v.known.insert("tigecycline".into());
        v.known.insert("clavulanic acid".into());
        v.folds.insert(
            "bactrim".into(),
            (
                "sulfamethoxazole-trimethoprim".into(),
                "sigdrug:brand_generic_folds.tsv".into(),
            ),
        );
        v.salts
            .insert("ceftiofur sodium".into(), "ceftiofur".into());
        v
    }

    #[test]
    fn resolves_submitter_variants() {
        let v = vocab();
        let r = resolve("Tigecyklin", &v);
        assert_eq!(
            (r.canonical.as_str(), r.kind),
            ("tigecycline", NameKind::Drug)
        );
        assert_eq!(r.components[0].source, "grm:edit-distance-2");
        assert_eq!(
            resolve("tigecycline", &v).components[0].source,
            "sigdrug:vocabulary"
        );
        assert_eq!(
            resolve("amoxicillin_clavulanat", &v).canonical,
            "amoxicillin/clavulanate"
        );
        assert_eq!(
            resolve("cefotaxime/clavulanic acid\u{c2}", &v).canonical,
            "cefotaxime/clavulanate"
        );
        assert_eq!(
            resolve("trimethoprim/sulfamethoxazole", &v).canonical,
            "sulfamethoxazole/trimethoprim"
        );
        assert_eq!(
            resolve("bactrim", &v).canonical,
            "sulfamethoxazole/trimethoprim"
        );
        assert_eq!(
            resolve("ceftiofur sodium", &v).components[0].source,
            "sigdrug:salt_base_ingredients.tsv"
        );
        assert_eq!(
            resolve("synercid", &v).canonical,
            "dalfopristin/quinupristin"
        );
        assert_eq!(resolve("geamycin", &v).kind, NameKind::Unresolved);
        assert_eq!(resolve("Carbapenem", &v).kind, NameKind::Class);
        assert_eq!(resolve("azidothymidine", &v).kind, NameKind::Excluded);
        assert_eq!(resolve("zzzqqq", &v).kind, NameKind::Unresolved);
    }

    #[test]
    fn classes_combinations_by_component() {
        assert_eq!(class_of("meropenem"), "carbapenem");
        assert_eq!(
            class_of("amoxicillin/clavulanate"),
            "penicillin + inhibitor"
        );
        assert_eq!(
            class_of("ceftazidime/avibactam"),
            "cephalosporin + inhibitor"
        );
        assert_eq!(
            class_of("sulfamethoxazole/trimethoprim"),
            "folate pathway combination"
        );
        assert_eq!(class_of("dalfopristin/quinupristin"), "streptogramin");
    }
}
