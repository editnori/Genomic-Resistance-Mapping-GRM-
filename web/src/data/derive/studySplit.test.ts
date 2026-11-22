import { describe, expect, it } from "vitest";
import { STUDY_LEARNERS, STUDY_SPLIT, studySplitSummary } from "./studySplit";

describe("study-split bundled summary", () => {
  it("matches the 70-dataset, 9-species record counts and per-learner readings", () => {
    const s = studySplitSummary(STUDY_SPLIT);
    expect(s.datasets).toBe(70);
    expect(s.species).toBe(9);
    expect(s.pairs).toBe(57_914);
    expect(s.costShown).toBe(5);
    expect(s.logLossWorse).toBe(5);
    expect(s.specificityWorse).toBe(5);
    expect(s.recallWorseLearners).toEqual(["rf", "gbt"]);
    expect(s.recallNotShownLearners).toEqual(["scm", "cart", "lasso"]);
    expect(s.lineageLogLossNotShown).toBe(4);
  });

  it("carries second-experiment SR, SM and MR draw-pair counts for all five learners", () => {
    for (const l of STUDY_LEARNERS) {
      expect(STUDY_SPLIT.second[l].SR.pairs).toBe(25);
      expect(STUDY_SPLIT.second[l].SR.cost).toBe(25);
      expect(STUDY_SPLIT.second[l].SM.pairs).toBe(5);
      expect(STUDY_SPLIT.second[l].SM.cost).toBe(5);
    }
  });
});
