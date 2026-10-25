import type { Lite, Manifest, PanelSummaryRow } from "../../data/schema";

const lite = (balanced_accuracy: number, sensitivity: number, specificity: number, vme: number, me: number, auroc: number | null): Lite => ({
  balanced_accuracy,
  sensitivity,
  specificity,
  vme,
  me,
  auroc,
});

export const REAL_PROBE_ROWS: PanelSummaryRow[] = [
  {
    panel: "gc_gyra",
    name: "gyrA codons 91-95 probes",
    dataset: "neisseria_gonorrhoeae__ciprofloxacin",
    species: "Neisseria gonorrhoeae",
    antibiotic: "ciprofloxacin",
    metrics: lite(0.9782575906122264, 0.9924433249370277, 0.9640718562874252, 0.007556675062972292, 0.03592814371257485, null),
    best_model: lite(0.9842456145643222, 0.9924433249370277, 0.9760479041916168, 0.007556675062972292, 0.023952095808383235, 0.9765758156231618),
    best_learner: "cart",
    counts: { n: 961, tp: 788, fp: 6, tn: 161, fn: 6 },
    best_model_counts: { n: 961, tp: 788, fp: 4, tn: 163, fn: 6 },
    no_call: { resistant: 0, susceptible: 0 },
  },
  {
    panel: "salmonella_qrdr",
    name: "gyrA QRDR probes",
    dataset: "salmonella_enterica__ciprofloxacin",
    species: "Salmonella enterica",
    antibiotic: "ciprofloxacin",
    metrics: lite(0.875976468007973, 0.93048128342246, 0.8214716525934861, 0.06951871657754011, 0.17852834740651388, null),
    best_model: lite(0.9303669321719863, 0.8835978835978836, 0.9771359807460891, 0.1164021164021164, 0.02286401925391095, 0.9555374095085286),
    best_learner: "gbt",
    counts: { n: 1203, tp: 348, fp: 148, tn: 681, fn: 26 },
    best_model_counts: { n: 1209, tp: 334, fp: 19, tn: 812, fn: 44 },
    no_call: { resistant: 4, susceptible: 2 },
  },
  {
    panel: "salmonella_qrdr",
    name: "gyrA QRDR probes",
    dataset: "salmonella_enterica__nalidixic_acid",
    species: "Salmonella enterica",
    antibiotic: "nalidixic acid",
    metrics: lite(0.9136172783122413, 0.8621359223300971, 0.9650986342943855, 0.1378640776699029, 0.03490136570561457, null),
    best_model: lite(0.9136518791228331, 0.8651252408477842, 0.962178517397882, 0.1348747591522158, 0.037821482602118005, 0.9040383724082447),
    best_learner: "cart",
    counts: { n: 1174, tp: 444, fp: 23, tn: 636, fn: 71 },
    best_model_counts: { n: 1180, tp: 449, fp: 25, tn: 636, fn: 70 },
    no_call: { resistant: 4, susceptible: 2 },
  },
];

const REAL_PROBE_DATASETS = [
  { id: "neisseria_gonorrhoeae__ciprofloxacin", species: "Neisseria gonorrhoeae", species_slug: "neisseria_gonorrhoeae", antibiotic: "ciprofloxacin", n: 961, n_resistant: 794 },
  { id: "salmonella_enterica__ciprofloxacin", species: "Salmonella enterica", species_slug: "salmonella_enterica", antibiotic: "ciprofloxacin", n: 1209, n_resistant: 378 },
  { id: "salmonella_enterica__nalidixic_acid", species: "Salmonella enterica", species_slug: "salmonella_enterica", antibiotic: "nalidixic acid", n: 1180, n_resistant: 519 },
];

export const withRealProbeRows = (m: Manifest): Manifest => ({
  ...m,
  datasets: REAL_PROBE_DATASETS.map((d) => ({ ...m.datasets[0]!, ...d, drug_class: "quinolone" })),
  panels: REAL_PROBE_ROWS,
});
