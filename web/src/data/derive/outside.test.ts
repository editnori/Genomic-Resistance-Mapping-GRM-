// The outside-study rows the site bundles, read as the documents read them.
// Each number here is one docs/outside-study-results.md states; the page
// computes it from the same rows. That the bundled file is exactly what
// reference/outside-study gives is the Python suite's check
// (tests/outside_study_results/test_site_data.py).

import { describe, expect, it } from "vitest";
import { OUTSIDE_STUDY, holds, outsidePoints, populationById, populationSummary } from "./outside";

const primary = populationSummary(populationById("primary"));
const same = populationSummary(populationById("same-study"));

describe("outside-study rows", () => {
  it("read revision 12's reports with 95% intervals", () => {
    expect(OUTSIDE_STUDY.evaluation_revision).toBe(12);
    expect(OUTSIDE_STUDY.confidence).toBe(0.95);
    expect(OUTSIDE_STUDY.populations.map((p) => p.id)).toEqual(["primary", "same-study"]);
  });

  it("give the primary population's size and each learner's holds over its own denominator", () => {
    expect([primary.frozen, primary.scored, primary.datasets, primary.species, primary.withEstimates]).toEqual([7258, 7258, 82, 7, 82]);
    expect(primary.holds).toEqual({ scm: { holds: 32, read: 81 }, cart: { holds: 23, read: 81 }, lasso: { holds: 7, read: 25 }, known_gene: { holds: 39, read: 74 } });
  });

  it("give the set covering machine's readings both ways on the primary population", () => {
    expect(primary.scm).toEqual({ sensitivityAbove: 36, sensitivityBelow: 27, specificityAbove: 25, specificityBelow: 35, onlySpecificity: 22, onlySensitivity: 14, both: 13 });
    expect(primary.oneStudyHalf).toBe(69);
  });

  it("give the same-study population apart", () => {
    expect([same.frozen, same.scored, same.datasets, same.species]).toEqual([3155, 3144, 49, 6]);
    expect(same.holds).toEqual({ scm: { holds: 36, read: 49 }, cart: { holds: 31, read: 49 }, lasso: { holds: 9, read: 15 }, known_gene: { holds: 30, read: 49 } });
  });

  it("draw one mark per set covering machine row read, with its reading", () => {
    for (const p of OUTSIDE_STUDY.populations) {
      const read = p.rows.filter((r) => r.learner === "scm" && holds(r) !== null).length;
      for (const m of ["sensitivity", "specificity"] as const) {
        const points = outsidePoints(p, m);
        expect(points).toHaveLength(read);
        // A mark reads below exactly when the report's figure is above its whole interval.
        for (const x of points) expect(x.reading === "below").toBe(x.report > x.figure.hi);
      }
    }
  });
});
