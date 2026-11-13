import { describe, expect, it } from "vitest";
import { INPUT_COMPARISON, armById, inputSummary, prediction, speciesIndex } from "./inputs";

const s = inputSummary();
const p4 = (x: number) => x.toFixed(4);

describe("the input comparison's bundled reading", () => {
  it("covers 72 datasets on four species, 49 with a lineage result", () => {
    expect(s.species).toBe(4);
    expect(s.datasets).toBe(72);
    expect(s.answered).toBe(49);
    expect(INPUT_COMPARISON.species.map((x) => [x.name, x.answered, x.datasets])).toEqual([
      ["Salmonella enterica", 18, 25],
      ["Staphylococcus aureus", 7, 18],
      ["Neisseria gonorrhoeae", 4, 7],
      ["Mycobacterium tuberculosis", 20, 22],
    ]);
  });

  it("gives the counts by dataset the document gives", () => {
    const cap1000 = armById("cap-1000").all;
    expect([cap1000.worse, cap1000.datasets, p4(cap1000.p)]).toEqual([25, 33, "0.0046"]);
    const budget1000 = armById("budget-1000").all;
    expect([budget1000.worse, budget1000.datasets, p4(budget1000.p)]).toEqual([37, 49, "0.0005"]);
    const cap4000 = armById("cap-4000");
    expect([cap4000.all.better, cap4000.all.datasets, p4(cap4000.all.p)]).toEqual([31, 49, "0.0854"]);
    expect(cap4000.species[speciesIndex("Salmonella enterica")]!.better).toBe(14);
  });

  it("reads the failed prediction on the two species it named", () => {
    expect(INPUT_COMPARISON.predictions.species).toEqual(["Staphylococcus aureus", "Mycobacterium tuberculosis"]);
    const t = prediction("cap-4000");
    expect([t.better, t.of]).toEqual([15, 27]);
    expect(t.p_one_sided).toBeGreaterThan(0.05);
  });

  it("counts the firm decisions: 23 of 1,945 tested, 5 better, 18 worse, 20 on the random forest", () => {
    expect([s.firm, s.tested, s.firmBetter, s.firmWorse, s.firmForest]).toEqual([23, 1945, 5, 18, 20]);
  });

  it("keeps the two runs of the default apart from the eleven changed inputs", () => {
    expect(s.changed).toHaveLength(11);
    expect(s.twins.map((a) => a.id)).toEqual(["sentinel-seed-twin", "sentinel-tuning-twin"]);
    for (const a of INPUT_COMPARISON.arms) {
      expect(a.species.reduce((n, x) => n + x.datasets, 0)).toBe(a.all.datasets);
      expect(a.all.better + a.all.worse + a.all.unchanged).toBe(a.all.datasets);
    }
  });
});
