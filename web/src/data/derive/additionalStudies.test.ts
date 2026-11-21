import { describe, expect, it } from "vitest";
import { ADDITIONAL_STUDIES, COMPARISON_LEARNERS } from "./additionalStudies";

describe("additionalStudies derived bundle", () => {
  it("loads all five revision 14 comparison records with expected learner keys", () => {
    expect(ADDITIONAL_STUDIES.evaluation_revision).toBe(14);
    expect(ADDITIONAL_STUDIES.time_split.datasets).toBe(24);
    expect(ADDITIONAL_STUDIES.time_split.total_log_loss_worse).toBe(5);
    expect(ADDITIONAL_STUDIES.kover.datasets).toBe(12);
    expect(ADDITIONAL_STUDIES.kover.balanced_accuracy.lower).toBe(8);
    expect(ADDITIONAL_STUDIES.dense_network.outcome).toBe("no gain established");
    expect(ADDITIONAL_STUDIES.amrfinderplus.datasets).toBe(235);
    expect(ADDITIONAL_STUDIES.amrfinderplus.program_vs_known_gene_rule.ahead_by_10pts).toBe(48);
    expect(ADDITIONAL_STUDIES.tbprofiler.graded_endpoints).toBe(14);
    for (const l of COMPARISON_LEARNERS) {
      expect(ADDITIONAL_STUDIES.time_split.by_learner[l]).toBeDefined();
      expect(ADDITIONAL_STUDIES.dense_network.by_learner[l]).toBeDefined();
      expect(ADDITIONAL_STUDIES.amrfinderplus.by_learner_lineage[l]).toBeDefined();
      expect(ADDITIONAL_STUDIES.tbprofiler.graded_vs_learners_lineage[l]).toBeDefined();
    }
  });
});
