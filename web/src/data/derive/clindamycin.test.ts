// The clindamycin kinds the dataset report shows (docs/clindamycin-labels.md).
// That the bundled file is exactly what the reference files give, and that
// each count equals the page's, is the Python suite's check
// (tests/clindamycin_labels). Here: the file validates, a broken copy does
// not, a note is given only for a report the check read, and the words say
// what each kind can be read as. The reports are the synthetic fixture's,
// renamed to one of the six with the file's revision and sample count.

import { describe, expect, it } from "vitest";
import { fixtureFile } from "../../dev/fixtures/files";
import { parseDataset, parseManifest } from "../validate";
import type { CatalogRow, DatasetReport } from "../schema";
import {
  CLINDAMYCIN_KINDS,
  CLINDAMYCIN_LABELS,
  clindamycinEntry,
  clindamycinMarkTip,
  clindamycinMethods,
  clindamycinNote,
  clindamycinRowEntry,
  kindCounts,
  kindOf,
  noPmidContrast,
  parseClindamycinLabels,
  type ClindamycinDataset,
  type ClindamycinLabels,
} from "./clindamycin";
import bundled from "../clindamycin-labels.json";

const manifest = parseManifest(fixtureFile("/data/index.json"));
const base = parseDataset(fixtureFile(`/data/datasets/${manifest.datasets[0]!.id}.json`));
const file = CLINDAMYCIN_LABELS!;
const copy = (): Record<string, unknown> => JSON.parse(JSON.stringify(bundled)) as Record<string, unknown>;
const datasetsOf = (raw: Record<string, unknown>) => raw.datasets as Record<string, Record<string, unknown>>;

const SA = "staphylococcus_aureus__clindamycin";
const CC = "campylobacter_coli__clindamycin";
const CJ = "campylobacter_jejuni__clindamycin";
const SG = "streptococcus_agalactiae__clindamycin";
const EF = "enterococcus_faecium__clindamycin";
const SP = "streptococcus_pneumoniae__clindamycin";

/** The fixture's report as the check read dataset `id`: its id, the file's revision and sample count. */
const reportAs = (id: string, over: { revision?: number; samples?: number } = {}): DatasetReport => ({
  ...base,
  id,
  evaluation: { revision: over.revision ?? file.evaluation_revision, valid: true, reasons: [] },
  counts: { ...base.counts, samples: over.samples ?? file.datasets[id]!.samples },
});
const rowAs = (id: string, over: { revision?: number; n?: number } = {}): CatalogRow => ({
  ...manifest.datasets[0]!,
  id,
  evaluation: { revision: over.revision ?? file.evaluation_revision, valid: true, reasons: [] },
  n: over.n ?? file.datasets[id]!.samples,
});
const words = (id: string): string => {
  const note = clindamycinNote(reportAs(id))!;
  return [note.title, ...note.counts, note.reading].join(" ");
};

describe("the bundled clindamycin file", () => {
  it("validates, names one release and holds six datasets of the four kinds", () => {
    expect(CLINDAMYCIN_LABELS).not.toBeNull();
    expect([file.release, file.evaluation_revision, file.page]).toEqual(["rev14-274-20261112", 14, "docs/clindamycin-labels.md"]);
    expect(Object.keys(file.datasets).sort()).toEqual([CC, CJ, EF, SA, SG, SP]);
    expect(new Set(Object.values(file.datasets).map((d) => d.kind))).toEqual(new Set(CLINDAMYCIN_KINDS));
  });

  it("gives each dataset the kind its counts give", () => {
    for (const d of Object.values(file.datasets)) expect(d.kind).toBe(kindOf(d));
    expect(kindCounts(file).reduce((n, k) => n + k.datasets, 0)).toBe(Object.keys(file.datasets).length);
  });

  it("refuses two datasets' kinds swapped", () => {
    const raw = copy();
    const ds = datasetsOf(raw);
    [ds[SA]!.kind, ds[CJ]!.kind] = [ds[CJ]!.kind, ds[SA]!.kind];
    expect(() => parseClindamycinLabels(raw)).toThrow(/kind is/);
  });

  it("refuses a count changed so that the parts no longer add up", () => {
    const raw = copy();
    datasetsOf(raw)[SP]!.study_samples = (datasetsOf(raw)[SP]!.study_samples as number) + 1;
    expect(() => parseClindamycinLabels(raw)).toThrow(/add up/);
  });

  it("refuses a field per sample, in a dataset or at the top", () => {
    const raw = copy();
    datasetsOf(raw)[SA]!.sample_ids = ["made-up-1", "made-up-2"];
    expect(() => parseClindamycinLabels(raw)).toThrow(/no field per sample/);
    const top = { ...copy(), samples: [{ id: "made-up", resistant: true }] };
    expect(() => parseClindamycinLabels(top)).toThrow(/expected exactly/);
  });

  it("refuses a dataset scored without its study when a class falls under the class minimum, and one not scored when both keep it", () => {
    // jejuni keeps 8 reported resistant samples without its study: made out to be scored, with tests and a kind that fit them.
    const scored = copy();
    Object.assign(datasetsOf(scored)[CJ]!, { scored_without: true, rows_passing: 0, rows_tested: 36, distinct_passing: 0, distinct_tested: 33, direction: null, kind: "no_dependence" });
    expect(() => parseClindamycinLabels(scored)).toThrow(/scored_without should say whether both classes keep the class minimum/);
    // pneumoniae keeps 41 and 392: made out to be not scored, with no tests and the kind that fits.
    const unscored = copy();
    Object.assign(datasetsOf(unscored)[SP]!, { scored_without: false, rows_passing: null, rows_tested: null, distinct_passing: null, distinct_tested: null, direction: null, kind: "class_balance" });
    expect(() => parseClindamycinLabels(unscored)).toThrow(/scored_without should say whether both classes keep the class minimum/);
  });

  it("refuses a released report's sha256 that is null, and parts by PMID that do not add up", () => {
    const sha = copy();
    datasetsOf(sha)[CJ]!.report_sha256 = null;
    expect(() => parseClindamycinLabels(sha)).toThrow(/report_sha256 should be a sha256/);
    const parts = copy();
    datasetsOf(parts)[CC]!.other_pmid_disagreements = (datasetsOf(parts)[CC]!.other_pmid_disagreements as number) - 1;
    expect(() => parseClindamycinLabels(parts)).toThrow(/other PMIDs should add up/);
    expect(Object.values(file.datasets).every((d) => /^[0-9a-f]{64}$/.test(d.report_sha256))).toBe(true);
  });

  it("refuses a source without its sha256, a passing test with no direction, and a not-scored dataset with tests", () => {
    const sources = copy();
    (sources.sources as Record<string, string>)["results/chosen.tsv"] = "not a hash";
    expect(() => parseClindamycinLabels(sources)).toThrow(/sha256/);
    const direction = copy();
    datasetsOf(direction)[CC]!.direction = null;
    expect(() => parseClindamycinLabels(direction)).toThrow(/direction/);
    const tests = copy();
    datasetsOf(tests)[CJ]!.distinct_passing = 3;
    expect(() => parseClindamycinLabels(tests)).toThrow(/no tests/);
  });
});

describe("when a note shows", () => {
  it("shows for each of the six reports the check read", () => {
    for (const id of Object.keys(file.datasets)) expect(clindamycinEntry(reportAs(id))).toBe(file.datasets[id]);
  });

  it("shows nothing for a dataset outside the six", () => {
    expect(clindamycinNote(base)).toBeNull();
    expect(clindamycinNote({ ...reportAs(SA), id: "staphylococcus_aureus__erythromycin" })).toBeNull();
    expect(clindamycinMarkTip(manifest.datasets[0]!)).toBeNull();
  });

  it("shows nothing when the evaluation revision or the sample count differs", () => {
    expect(clindamycinNote(reportAs(SA, { revision: file.evaluation_revision + 1 }))).toBeNull();
    expect(clindamycinNote({ ...reportAs(SA), evaluation: undefined })).toBeNull();
    expect(clindamycinNote(reportAs(SA, { samples: file.datasets[SA]!.samples - 1 }))).toBeNull();
    expect(clindamycinRowEntry(rowAs(SA, { revision: file.evaluation_revision + 1 }))).toBeNull();
    expect(clindamycinRowEntry(rowAs(SA, { n: file.datasets[SA]!.samples + 1 }))).toBeNull();
    expect(clindamycinRowEntry(rowAs(SA))).toBe(file.datasets[SA]);
  });

  it("shows nothing when there is no valid file", () => {
    expect(clindamycinNote(reportAs(SA), null)).toBeNull();
    expect(clindamycinMarkTip(rowAs(SA), null)).toBeNull();
    expect(clindamycinMethods({ evaluation_revision: file.evaluation_revision }, null)).toBeNull();
  });
});

describe("the note's words", () => {
  it("names each kind as the page does", () => {
    expect(clindamycinNote(reportAs(SA))!.title).toBe("Clindamycin labels by study: measured dependence");
    expect(clindamycinNote(reportAs(EF))!.title).toBe("Clindamycin labels by study: one study");
    expect(clindamycinNote(reportAs(CJ))!.title).toBe("Clindamycin labels by study: class balance, dependence not measured");
    expect(clindamycinNote(reportAs(SP))!.title).toBe("Clindamycin labels by study: no dependence under the rule");
  });

  it("says dependence was not measured for the two class-balance datasets, and reads as no measured effect", () => {
    for (const id of [CJ, SG]) {
      const text = words(id);
      expect(text).toContain("dependence not measured");
      expect(text).toContain("was not measured");
      expect(text).toContain("cannot be scored without that study");
      expect(text).not.toMatch(/\b(rise|rises|fall|falls)\b|figures tested|measured dependence/);
    }
  });

  it("says the passing changes fall for the Campylobacter coli dataset and rise for the Staphylococcus aureus one", () => {
    expect(words(CC)).toMatch(/distinct figures tested fall by/);
    expect(words(CC)).not.toMatch(/\brises?\b/);
    expect(words(SA)).toMatch(/distinct figures tested rise by/);
    expect(words(SA)).not.toMatch(/\bfalls?\b/);
  });

  it("says the study's share stands out only against the samples with no PMID by the writer's rule, for the two Campylobacter datasets alone", () => {
    const said = "and the study's share of the disagreements stands out only against those samples.";
    for (const id of Object.keys(file.datasets)) expect(words(id).includes(said), id).toBe(noPmidContrast(file.datasets[id]!));
    expect(Object.keys(file.datasets).filter((id) => noPmidContrast(file.datasets[id]!)).sort()).toEqual([CC, CJ]);
    expect(words(CC)).toContain(`294 of the 389 samples carry no PMID, ${said}`);
    expect(words(CJ)).toContain(`475 of the 508 samples carry no PMID, ${said}`);
    // Each of the rule's three parts is needed.
    const coli = file.datasets[CC]!;
    expect(noPmidContrast({ ...coli, no_pmid_samples: Math.floor(coli.samples / 2) })).toBe(false);
    expect(noPmidContrast({ ...coli, no_pmid_disagreements: coli.no_pmid_samples })).toBe(false);
    expect(noPmidContrast({ ...coli, other_pmid_disagreements: 0 })).toBe(false);
  });

  it("names the figures that pass by their count, and says the others do not pass", () => {
    expect(clindamycinNote(reportAs(CC))!.reading).toBe(
      "Those 9 figures are associated with that one study's labels or genomes. The other 24 do not pass the rule, and neither result shows which reading is wrong.",
    );
    expect(clindamycinNote(reportAs(SA))!.reading).toBe(
      "Those 19 figures are associated with that one study's labels or genomes. The other 14 do not pass the rule, and neither result shows which reading is wrong.",
    );
    for (const id of [CC, SA]) expect(words(id)).not.toContain("The released figure");
  });

  it("does not read no passing figure as no dependence", () => {
    const reading = clindamycinNote(reportAs(SP))!.reading;
    expect(reading).toBe("This does not show that the figures do not depend on that study: some intervals are wide, and the fold plans of this dataset were not varied.");
    expect(words(SP)).toContain("none of the 33 distinct figures tested changes by 0.05 or more");
    expect(words(SP)).not.toMatch(/barely move|do not depend on that study\./);
  });

  it("says one passing figure, and one other, in the singular", () => {
    const one = (d: Partial<ClindamycinDataset>, id = CC): ClindamycinLabels => ({ ...file, datasets: { ...file.datasets, [id]: { ...file.datasets[id]!, ...d } } });
    const falls = clindamycinNote(reportAs(CC), one({ rows_passing: 1, distinct_passing: 1 }))!;
    expect(falls.counts.join(" ")).toContain("1 of the 33 distinct figures tested falls by 0.05 or more");
    expect(falls.reading).toBe("That figure is associated with that one study's labels or genomes. The other 32 do not pass the rule, and neither result shows which reading is wrong.");
    const rises = clindamycinNote(reportAs(SA), one({ rows_passing: 1, distinct_passing: 1 }, SA))!;
    expect(rises.counts.join(" ")).toContain("1 of the 33 distinct figures tested rises by 0.05 or more");
    expect(clindamycinNote(reportAs(SA), one({ rows_passing: 32, distinct_passing: 32 }, SA))!.reading).toContain("The other 1 does not pass the rule");
    expect(clindamycinNote(reportAs(SA), one({ rows_passing: 36, distinct_passing: 33 }, SA))!.reading).toBe(
      "Those 33 figures are associated with that one study's labels or genomes, and the result does not show which reading is wrong.",
    );
  });

  it("says the class minimum is the release's, with the file's number", () => {
    expect(words(CJ)).toContain(`Without it 8 reported resistant samples remain, under the class minimum of ${file.min_class} this release uses, so`);
    expect(words(SG)).toContain(`Without it 18 reported susceptible samples remain, under the class minimum of ${file.min_class} this release uses, so`);
    expect(words(CJ)).not.toContain("GRM needs");
  });

  it("holds the one-study reading to the page's sentences for that dataset, word for word", async () => {
    const fs = (await import(/* @vite-ignore */ `node:${"fs"}`)) as { readFileSync(path: string, encoding: string): string };
    const cwd = (globalThis as unknown as { process: { cwd(): string } }).process.cwd();
    const page = fs.readFileSync(`${cwd}/../${file.page}`, "utf8");
    const line = page.split("\n").find((l) => l.startsWith(`- **\`${EF}\`: one study.**`));
    expect(line).toBeDefined();
    // "... so the figures describe that study's genomes and labels alone. No figure without that study can exist."
    const [first, second] = clindamycinNote(reportAs(EF))!.reading.replace(/\.$/, "").split(", and ");
    expect(line).toContain(`${first!.charAt(0).toLowerCase()}${first!.slice(1)}.`);
    expect(line).toContain(`${second!.charAt(0).toUpperCase()}${second!.slice(1)}.`);
  });

  it("calls the threshold a research screen, and uses none of the barred words", () => {
    for (const id of Object.keys(file.datasets)) {
      const text = words(id);
      expect(text).not.toMatch(/carrier|provenance|AMRFinderPlus|validat|clinical|patient|error|wrong label|mislabel|caus/i);
      if (file.datasets[id]!.scored_without) expect(text).toContain(`(${file.threshold} is an adapted research screen chosen for this check)`);
    }
  });
});

describe("the Methods paragraph", () => {
  it("counts the datasets of each kind from the file, for the release the check read only", () => {
    const text = clindamycinMethods({ evaluation_revision: file.evaluation_revision })!;
    const counts = Object.fromEntries(kindCounts(file).map((k) => [k.kind, k.datasets]));
    expect(text).toContain(`measured dependence (${counts.measured_dependence} datasets); one study (${counts.one_study}); class balance, dependence not measured (${counts.class_balance}); and no dependence under the rule (${counts.no_dependence})`);
    expect(text).toContain(`release ${file.release}`);
    expect(text).toContain("neither is called the error");
    expect(text).toContain("The page sorts the 6 datasets into four kinds");
    expect(text).not.toContain("The check sorts");
    expect(clindamycinMethods({ evaluation_revision: file.evaluation_revision + 1 })).toBeNull();
    expect(clindamycinMethods({ evaluation_revision: undefined })).toBeNull();
  });

  it("would say one dataset in the singular", () => {
    const one: ClindamycinLabels = { ...file, datasets: { [SA]: file.datasets[SA]! } };
    expect(clindamycinMethods({ evaluation_revision: file.evaluation_revision }, one)).toContain("measured dependence (1 dataset); one study (0)");
  });
});
