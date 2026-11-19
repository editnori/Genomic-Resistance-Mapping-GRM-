// The evaluation audit: whether a report's metrics pass the internal
// evaluation audit, and what a pass covers. A report written before the
// audit carries none and reads exactly as report.rs labels it, provisional.

import type { EvaluationAudit, Manifest } from "../schema";

/** The evaluation revision whose reports leave null fits uncalled. Set it to
 *  the revision src/analysis.rs bumps to when it ships those reports. */
export const NO_CALL_REVISION = 7;

/** How a null fit counts under the manifest's evaluation revision. */
/** The first evaluation revision that scores a null fit's isolates at their
 *  fold's training prevalence (src/analysis.rs). */
export const EQUAL_SCORING_REVISION = 14;

/** The first evaluation revision whose gene-panel rows count only the
 *  targets a panel reads for the drug, in a species the assay reports them
 *  for (src/panels.rs). Earlier rows count every target of the panel. */
export const PANEL_TARGETS_REVISION = 13;

export interface DefWords {
  plain: string;
  detail: string;
}

const MARKER = "A detected gene is what the panel would report; it does not prove that the gene causes the isolate's resistance.";

/** How a simulated panel calls, and what the blood-culture and CTX-M rows
 *  count, under the manifest's evaluation revision. */
export const panelWords = (revision: number | undefined): { panels: DefWords; bloodCulture: DefWords } =>
  revision !== undefined && revision >= PANEL_TARGETS_REVISION
    ? {
        panels: {
          plain: "Commercial lab gene tests, imitated on the genome sequence: a panel calls resistant when it finds a gene that it reads for that drug in that species.",
          detail: `A gene panel calls resistant when the assembled genome carries one of its targets that the panel reads for the drug, in a species the assay reports that target for. Targets, drugs and species follow the manufacturers' published documents. A probe assay calls resistant when a probed region differs from the wild-type sequence. Detection is in silico, from assembled genomes. ${MARKER}`,
        },
        bloodCulture: {
          plain: "Two rows that read the same blood-culture gene test in two ways: every gene on it that breaks down the drug, or the CTX-M gene alone.",
          detail:
            "On a cephalosporin the two rows read the same panel two ways. The blood-culture row counts every gene of the panel that breaks down the drug: CTX-M, and also the carbapenemase genes KPC, NDM, VIM and IMP, as published studies of these enzymes describe. It shows what the panel's genes could explain. The CTX-M row counts CTX-M alone, as the maker compares the panel's cephalosporin result. For aztreonam the blood-culture row adds KPC only, since NDM, VIM and IMP do not break down aztreonam. OXA-48-like is read for carbapenems only.",
        },
      }
    : {
        panels: {
          plain: "Commercial lab gene tests, imitated on the genome sequence: a panel calls resistant when it finds any gene it looks for.",
          detail: `A gene panel calls resistant when the assembled genome carries any of its targets, which follow the manufacturers' published target lists. In these reports a panel's row counts every target of the panel on every drug the row is shown for, in every species. Reports from evaluation revision ${PANEL_TARGETS_REVISION} count only the targets the panel reads for that drug, in a species the assay reports them for. A probe assay calls resistant when a probed region differs from the wild-type sequence. Detection is in silico, from assembled genomes. ${MARKER}`,
        },
        bloodCulture: {
          plain: "Two rows that read the same blood-culture gene test in two ways: every gene on it, or the CTX-M gene alone.",
          detail: `On a cephalosporin or aztreonam the blood-culture row counts every target of the panel: CTX-M, the carbapenemase genes KPC, NDM, VIM, IMP and OXA-48-like, and mcr-1. The CTX-M row counts CTX-M alone. Reports from evaluation revision ${PANEL_TARGETS_REVISION} limit the blood-culture row to the genes the panel reads for the drug.`,
        },
      };

export const nullFitWords = (revision: number | undefined): string =>
  revision !== undefined && revision >= EQUAL_SCORING_REVISION
    ? `A boosted model that keeps no tree, or a lasso that keeps no weight, has learned nothing beyond the class balance. From evaluation revision ${EQUAL_SCORING_REVISION} it scores each isolate of that fold at the share of resistant isolates in the fold's training set, and the 0.5 cutoff calls that score resistant when the share is 0.5 or more. Balanced accuracy, VME, ME, and AUROC cover every isolate for every learner, and each dataset's learner table counts the isolates scored this way. Reports computed under revisions ${NO_CALL_REVISION} to ${EQUAL_SCORING_REVISION - 1} made no call for those isolates and scored the rest; earlier revisions scored them 0.5. A fitted model that scores an isolate exactly 0.5 calls it resistant.`
    : revision !== undefined && revision >= NO_CALL_REVISION
    ? `A boosted model that keeps no tree, or a lasso that keeps no weight, has learned nothing beyond the class balance. From evaluation revision ${NO_CALL_REVISION} it makes no call for the isolates of that fold: balanced accuracy, VME, ME, and AUROC cover the isolates called, and each dataset's learner table gives the share called. Reports computed under earlier revisions scored those isolates 0.5, which the 0.5 cutoff calls resistant. A fitted model that scores an isolate exactly 0.5 still calls it resistant.`
    : "A boosted model that keeps no tree, or a lasso that keeps no weight, has learned nothing beyond the class balance. It scores every isolate of that fold 0.5, which the 0.5 cutoff calls resistant.";

export const LEGACY_AUDIT: EvaluationAudit = {
  revision: 0,
  valid: false,
  reasons: ["This report predates the evaluation audit; its metrics are provisional and must be recomputed."],
};

const staleReason = (revision: number, current: number): string =>
  `Computed under evaluation audit revision ${revision}; the current revision is ${current}, so these metrics are provisional until the dataset is recomputed.`;

export const auditOf = (x: { evaluation?: EvaluationAudit }, current?: number): EvaluationAudit => {
  const a = x.evaluation ?? LEGACY_AUDIT;
  return current === undefined || a.revision === 0 || a.revision === current ? a : { ...a, valid: false, reasons: [staleReason(a.revision, current), ...a.reasons] };
};

export interface ManifestAudit {
  total: number;
  /** Datasets whose evaluation does not pass the audit. */
  invalid: number;
  /** Distinct reasons among those, most common first. */
  reasons: string[];
  /** Distinct scopes among the datasets that pass, most common first. */
  scopes: string[];
}

const byCount = (counts: Map<string, number>): string[] => [...counts.entries()].sort((a, b) => b[1] - a[1]).map(([k]) => k);

export const manifestAudit = (m: Manifest): ManifestAudit => {
  const reasons = new Map<string, number>();
  const scopes = new Map<string, number>();
  let invalid = 0;
  m.datasets.forEach((d) => {
    const a = auditOf(d, m.evaluation_revision);
    if (a.valid) {
      if (a.scope !== undefined && a.scope !== "") scopes.set(a.scope, (scopes.get(a.scope) ?? 0) + 1);
      return;
    }
    invalid += 1;
    a.reasons.forEach((r) => reasons.set(r, (reasons.get(r) ?? 0) + 1));
  });
  return { total: m.datasets.length, invalid, reasons: byCount(reasons), scopes: byCount(scopes) };
};
