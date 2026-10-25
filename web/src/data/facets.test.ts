import { describe, expect, it } from "vitest";
import { applyAxes, countBy, memberTest, type Axis } from "./facets";

interface Row {
  species: string;
  drugClass: string;
}

const rows: Row[] = [
  { species: "A", drugClass: "carbapenem" },
  { species: "A", drugClass: "quinolone" },
  { species: "B", drugClass: "carbapenem" },
  { species: "B", drugClass: "carbapenem" },
];

describe("predicting counts", () => {
  const axes = (species: string[], classes: string[]): Axis<Row>[] => [
    { id: "species", test: memberTest(new Set(species), (r) => r.species) },
    { id: "class", test: memberTest(new Set(classes), (r) => r.drugClass) },
  ];

  it("counts a facet against the other axes with its own selection lifted", () => {
    // Species B is chosen, so the class counts only see species B rows...
    expect(countBy(rows, axes(["B"], []), "class", (r) => r.drugClass)).toEqual(new Map([["carbapenem", 2]]));
    // ...while the species counts ignore the species choice itself.
    expect(countBy(rows, axes(["B"], ["carbapenem"]), "species", (r) => r.species)).toEqual(
      new Map([
        ["A", 1],
        ["B", 2],
      ]),
    );
  });

  it("treats an empty many-of selection as every value", () => {
    expect(applyAxes(rows, axes([], []))).toHaveLength(4);
    expect(applyAxes(rows, axes(["A"], ["carbapenem"]))).toHaveLength(1);
  });
});
