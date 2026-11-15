// Synthetic baselines files for the tests: counts written by hand, whose
// figures are their own counts' as parseBaselines requires. Not BV-BRC data.

export const SHA = "a".repeat(64);

export const agreement = (agree: number, n: number) => ({ agree, n, essential_agreement: agree / Math.max(n, 1) });

/** A lineage entry whose rates are its counts' own. */
export const lineageEntry = (scheme: "random" | "lineage", n: number, tp: number, fp: number, tn: number, unseen: number, folds = 5) => {
  const fn = n - tp - fp - tn;
  const sensitivity = tp + fn === 0 ? 0 : tp / (tp + fn);
  const specificity = tn + fp === 0 ? 0 : tn / (tn + fp);
  return { scheme, folds, n, tp, fp, tn, fn, balanced_accuracy: (sensitivity + specificity) / 2, sensitivity, specificity, unseen, unseen_share: unseen / Math.max(n, 1) };
};

/** An MIC part over `nExact` exact MICs dealt into five folds. */
export const micPart = (nExact: number, model: number, median: number, mode: number, confirmed = true) => {
  const tests = [0, 1, 2, 3, 4].map((f) => Math.floor(nExact / 5) + (f < nExact % 5 ? 1 : 0));
  return {
    folds: 5,
    seed: 42,
    prediction_check: confirmed ? { confirmed: true, groups: 4, isolates: 9, spanning: 0, constant_folds: [] } : { confirmed: false, groups: 0, isolates: 0, spanning: 0, constant_folds: [0, 1, 2, 3, 4] },
    threshold_check: { confirmed: true, checked: 5 },
    exact_with_training_copy: 2,
    n: nExact + 12,
    n_exact: nExact,
    model: agreement(model, nExact),
    median: agreement(median, nExact),
    mode: agreement(mode, nExact),
    fold_constants: tests.map((t) => ({ train_exact: nExact - t, test_exact: t, median: 2, upper_median_differs: false, mode: 1, mode_tied: false })),
  };
};

export const file = (over: Record<string, unknown> = {}) => ({
  schema: 1,
  dataset: "synthetic__fixturemycin",
  report_sha256: SHA,
  evaluation_revision: 14,
  mic: micPart(96, 58, 81, 60),
  lineage: [lineageEntry("random", 120, 50, 10, 45, 8), lineageEntry("lineage", 120, 30, 30, 30, 120)],
  ...over,
});


/** One entry of the summary beside the index (baselines.json). */
export const summaryEntry = (dataset: string, model: number, constant: number, confirmed = true) => ({ confirmed, constant, dataset, model, report_sha256: SHA });

/** A summary pinned to the index bytes with SHA-256 `indexSha`. */
export const summary = (indexSha: string, datasets: readonly ReturnType<typeof summaryEntry>[], revision = 14) => ({
  datasets,
  evaluation_revision: revision,
  index_sha256: indexSha,
  schema: 1,
});
