import { describe, expect, it } from "vitest";
import type { GenomeRow, SpeciesReport } from "../data/schema";
import { fixtureFile } from "../dev/fixtures/files";
import { parseSpecies } from "../data/validate";
import { MIN_JOINT, crossPairs, jointCalls, pairSummary, undefinedWhy, zeroMargin } from "./crossPairs";

const genome = (i: number, profile: string): GenomeRow => ({
  id: `G${i}`,
  lineage: 0,
  x: 0,
  y: 0,
  nearest: 0,
  country: null,
  year: null,
  host: null,
  mlst: null,
  length: 0,
  contigs: 0,
  n50: 0,
  groups: [],
  profile,
});

/** `count` genomes with each profile. */
const many = (spec: readonly [profile: string, count: number][]): GenomeRow[] => spec.flatMap(([p, c]) => Array.from({ length: c }, () => p)).map((p, i) => genome(i, p));

const species = (genomes: GenomeRow[], cross: SpeciesReport["cross"]): SpeciesReport => ({
  schema: 1,
  species: "Test species",
  slug: "test",
  drugs: ["amikacin", "colistin", "meropenem", "tetracycline"],
  drug_classes: ["aminoglycoside", "polymyxin", "carbapenem", "tetracycline"],
  catalog_groups: [],
  lineage_threshold: 0.001,
  n_lineages: 1,
  leaf_order: genomes.map((_, i) => i),
  merges: [],
  genomes,
  cross,
});

describe("crossPairs", () => {
  // Amikacin and colistin: 40 jointly called genomes; meropenem is called on
  // 20; tetracycline is susceptible everywhere; three genomes carry only an
  // intermediate call.
  const genomes = many([
    ["RRSS", 10],
    ["RSRS", 10],
    ["SR.S", 10],
    ["SS.S", 10],
    ["I...", 3],
  ]);
  const listed: SpeciesReport["cross"] = [
    { a: "amikacin", b: "colistin", n: 40, rr: 10, rs: 10, sr: 10, ss: 10, phi: 0 },
    // Written in the other order: its one-sided cells swap.
    { a: "tetracycline", b: "amikacin", n: 40, rr: 0, rs: 0, sr: 20, ss: 20, phi: 0 },
    { a: "colistin", b: "tetracycline", n: 40, rr: 0, rs: 20, sr: 0, ss: 20, phi: 0 },
  ];
  const rows = crossPairs(species(genomes, listed));

  it("lists every pair once, in drug order", () => {
    expect(rows.map((p) => `${p.a}|${p.b}`)).toEqual([
      "amikacin|colistin",
      "amikacin|meropenem",
      "amikacin|tetracycline",
      "colistin|meropenem",
      "colistin|tetracycline",
      "meropenem|tetracycline",
    ]);
  });

  it("keeps the file's counts for listed pairs, flipping a pair written the other way", () => {
    const at = rows.find((p) => p.key === "amikacin|colistin")!;
    expect(at).toMatchObject({ n: 40, rr: 10, rs: 10, sr: 10, ss: 10, phi: 0, status: "evaluated" });
    const flipped = rows.find((p) => p.key === "amikacin|tetracycline")!;
    expect(flipped).toMatchObject({ rr: 0, rs: 20, sr: 0, ss: 20 });
  });

  it("reads a zero margin as undefined, never as phi 0", () => {
    const flipped = rows.find((p) => p.key === "amikacin|tetracycline")!;
    expect(flipped.status).toBe("undefined");
    expect(flipped.phi).toBeNull();
    expect(undefinedWhy(flipped)).toBe("none is resistant to tetracycline");
    const both = rows.find((p) => p.key === "colistin|tetracycline")!;
    expect(both.status).toBe("undefined");
    expect(undefinedWhy(both)).toBe("none is resistant to tetracycline");
  });

  it("tallies unlisted pairs from the profiles, R and S only", () => {
    const mero = rows.find((p) => p.key === "amikacin|meropenem")!;
    // Only RRSS and RSRS are called on meropenem; I and untested never count.
    expect(mero).toMatchObject({ n: 20, rr: 10, rs: 10, sr: 0, ss: 0, status: "not-evaluated", phi: null });
    expect(mero.n).toBeLessThan(MIN_JOINT);
    expect(rows.find((p) => p.key === "colistin|meropenem")).toMatchObject({ n: 20, rr: 0, rs: 10, sr: 10, ss: 0 });
  });

  it("summarizes evaluated, undefined and not evaluated pairs and the highest defined phi", () => {
    const s = pairSummary(rows);
    expect(s).toMatchObject({ total: 6, evaluated: 3, undefined: 2, notEvaluated: 3 });
    expect(s.highest?.key).toBe("amikacin|colistin");
  });

  it("matches report.rs on the fixture: listed pairs keep their counts, the rest fall under the minimum", () => {
    for (const slug of ["klebsiella_pneumoniae", "staphylococcus_aureus"]) {
      const s = parseSpecies(fixtureFile(`/data/species/${slug}.json`));
      const all = crossPairs(s);
      expect(all).toHaveLength((s.drugs.length * (s.drugs.length - 1)) / 2);
      s.cross.forEach((c) => {
        const row = all.find((p) => p.key === `${c.a}|${c.b}`)!;
        expect([row.n, row.rr, row.rs, row.sr, row.ss]).toEqual([c.n, c.rr, c.rs, c.sr, c.ss]);
        expect(row.status).toBe(zeroMargin(c) ? "undefined" : "evaluated");
      });
      all
        .filter((p) => p.status === "not-evaluated")
        .forEach((p) => {
          expect(p.n).toBeLessThan(MIN_JOINT);
          expect(jointCalls(s, p.ia, p.ib)).toEqual({ n: p.n, rr: p.rr, rs: p.rs, sr: p.sr, ss: p.ss });
        });
    }
  });
});
