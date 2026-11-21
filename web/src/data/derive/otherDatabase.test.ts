import { describe, expect, it } from "vitest";
import { OTHER_DATABASE, otherDatabaseSummary } from "./otherDatabase";

describe("other-database bundled summary", () => {
  it("matches the 13,377-genome NCBI Pathogen Detection record", () => {
    const s = otherDatabaseSummary(OTHER_DATABASE);
    expect(s.totalGenomes).toBe(13_377);
    expect(s.totalSpecies).toBe(9);
    expect(s.pilotGenomes).toBe(3_897);
    expect(s.eightGenomes).toBe(9_480);
    expect(s.eightSpecies).toBe(8);
    expect(s.totals.scm.plan_meets).toBe(6);
    expect(s.totals.scm.holds).toBe(52);
    expect(s.totals.scm.estimable).toBe(84);
    expect(s.totals.cart.plan_meets).toBe(5);
    expect(s.totals.cart.holds).toBe(50);
    expect(s.totals.lasso.plan_meets).toBe(3);
    expect(s.totals.lasso.assessable).toBe(4);
    expect(s.totals.lasso.holds).toBe(17);
    expect(s.totals.lasso.estimable).toBe(24);
    expect(s.pilot.scm.holds).toBe(14);
    expect(s.pilot.scm.estimable).toBe(25);
    expect(s.pilot.cart.holds).toBe(12);
  });
});
