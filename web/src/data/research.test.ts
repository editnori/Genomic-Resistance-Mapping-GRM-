// @vitest-environment happy-dom

import { act, createElement } from "react";
import { createRoot, type Root } from "react-dom/client";
import { afterEach, describe, expect, it, vi } from "vitest";
import { fixtureFile } from "../dev/fixtures/files";
import { micEntry, micRun, researchWith, SELECTIVE_CODES, SELECTIVE_RESULT, selectiveEntry, WITHHELD_REASON } from "../dev/fixtures/research";
import { PATHS, useDataset, useResearch, useResource, type Resource } from "./load";
import { researchFor, reportGuard, sha256Hex, type DatasetResearch } from "./research";
import type { DatasetReport, MicDataset, MicRun, MicWithheld, Research, SelectiveCategory } from "./schema";
import { parseDataset, parseResearch } from "./validate";

Object.assign(globalThis, { IS_REACT_ACT_ENVIRONMENT: true });

const ID = "klebsiella_pneumoniae__meropenem";
const report = (): DatasetReport => parseDataset(fixtureFile(`/data/datasets/${ID}.json`));
/** A fixture report without an MIC model, so its guard holds no predictions. */
const unmodeled = (): DatasetReport => parseDataset(fixtureFile("/data/datasets/klebsiella_pneumoniae__ceftriaxone.json"));
const copy = <T,>(v: T): T => JSON.parse(JSON.stringify(v)) as T;
/** research.json as fetched, typed loosely enough to break. */
const raw = <T,>(r: Research): T => copy(r) as unknown as T;
const ready = (data: Research): Resource<Research> => ({ status: "ready", data });

/** The failed guard parts, or the status when nothing failed. */
const failed = (s: DatasetResearch["mic"] | DatasetResearch["selective"]): string | string[] => (s.status === "mismatch" ? s.failed : s.status);

let root: Root | null = null;

afterEach(() => {
  act(() => root?.unmount());
  root = null;
  vi.unstubAllGlobals();
});

describe("research guard", () => {
  it("hashes with SHA-256", () => {
    expect(sha256Hex("")).toBe("e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    expect(sha256Hex("abc")).toBe("ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    // 56 bytes: the length no longer fits in the first block.
    expect(sha256Hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq")).toBe("248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
    expect(sha256Hex("a".repeat(1_000_000))).toBe("cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
    expect(sha256Hex("é")).toBe("4a99557e4033c3539de2eb65472017cad5f9557f7a0625a09f1c3f6e2ba69c4c");
  });

  it("computes the golden digests that src/research.rs and Python's hashlib give", () => {
    const golden = {
      samples: [
        { id: "562.1", resistant: true, from_mic: false, mic: { log2: 3, censor: "Above" }, lineage: 0 },
        { id: "562.22", resistant: false, from_mic: true, mic: null, lineage: 2 },
        { id: "562.3", resistant: false, from_mic: false, mic: { log2: -2, censor: "Exact" }, lineage: 0 },
      ],
      mic: { predicted: [2.75, null, -1.5] },
      splits: [
        { scheme: "random", learners: [{ learner: "gbt", scores: [0.5, 0.5, 0.5] }] },
        { scheme: "lineage", learners: [{ learner: "gbt", scores: [0.9, 0.1, 0.25] }, { learner: "catalog", scores: [1, 0, 0] }] },
      ],
    } as unknown as DatasetReport;
    expect(reportGuard(golden)).toEqual({
      samples: 3,
      sample_ids: "04d55a6c086cda50b8411babb85468ef748eeffccff37e86131d6292312563b1",
      labels: "db5593b64018f1438d8c6813be1a0cfa5f306815b3b1f4f61614f1719599d794",
      lineages: "04ed82dae8b902d411accf765f4ded9af83a46b71ffb1fce5f02c611908e6bcd",
      mic_values: "3937937471419ea303f03db2232fe979f5b4fcd79ee03207ba1b6394273b8113",
      mic_predicted: "e9f5ea58df3eacbc7f595ba68475d7d67ae83d4d96bfc787bfad4bc25e66b577",
      scores: "40adb03707afe2d0fe731c7e8e48010a481daf5a4a8fc8efc00722ce5a245535",
    });
  });

  it("applies an entry to the report whose cohort it guards", () => {
    const r = report();
    const state = researchFor(ready(parseResearch(raw(researchWith([micEntry(r)], [selectiveEntry(r, "A1")])))), r);
    expect(state.mic.status).toBe("current");
    expect(state.selective.status).toBe("current");
  });

  it("hashes a report object once however often its research is read", () => {
    const r = report();
    const research = ready(researchWith([micEntry(r)], [selectiveEntry(r, "A1")]));
    let hashes = 0;
    const counted = (target: DatasetReport): DatasetReport =>
      new Proxy(target, {
        get: (t, key) => {
          // reportGuard reads samples once per hash.
          if (key === "samples") hashes += 1;
          return Reflect.get(t, key);
        },
      });
    const shown = counted(r);
    // Both lanes, read by two inspector sections on each of three renders.
    for (let i = 0; i < 3; i += 1) {
      expect(researchFor(research, shown)).toMatchObject({ mic: { status: "current" }, selective: { status: "current" } });
      expect(researchFor(research, shown).selective.status).toBe("current");
    }
    expect(hashes).toBe(1);
    // Another object is hashed afresh, so an edited copy still fails its guard.
    const changed = copy(r);
    changed.samples[0]!.id = "573.99999";
    expect(researchFor(research, counted(changed)).mic).toEqual({ status: "mismatch", failed: ["sample_ids"], researchModel: true });
    expect(hashes).toBe(2);
  });

  it("rejects a changed sample order, id, label, MIC, saved MIC prediction, or held-out score", () => {
    const base = report();
    const research = ready(researchWith([micEntry(base)], [selectiveEntry(base, "A1")]));
    const change = (edit: (r: DatasetReport) => void) => {
      const r = copy(base);
      edit(r);
      const s = researchFor(research, r);
      return { mic: failed(s.mic), selective: failed(s.selective) };
    };
    const withMic = base.samples.findIndex((s) => s.mic !== null);
    const predicted = base.mic!.predicted.findIndex((p) => p !== null);
    const scheme = (r: DatasetReport, name: string) => r.splits.find((s) => s.scheme === name)!.learners[0]!.scores;

    expect(change((r) => r.samples.reverse())).toMatchObject({ mic: expect.arrayContaining(["sample_ids"]), selective: expect.arrayContaining(["sample_ids"]) });
    expect(change((r) => (r.samples[0]!.id = "573.99999"))).toEqual({ mic: ["sample_ids"], selective: ["sample_ids"] });
    expect(change((r) => (r.samples[0]!.resistant = !r.samples[0]!.resistant))).toEqual({ mic: ["labels"], selective: ["labels"] });
    expect(change((r) => (r.samples[0]!.from_mic = !r.samples[0]!.from_mic))).toEqual({ mic: ["labels"], selective: ["labels"] });
    // Each lane fails only on the values it read.
    expect(change((r) => (r.samples[withMic]!.mic!.log2 += 1))).toEqual({ mic: ["mic_values"], selective: "current" });
    expect(change((r) => (r.mic!.predicted[predicted]! += 1e-9))).toEqual({ mic: ["mic_predicted"], selective: "current" });
    expect(change((r) => (scheme(r, "lineage")[0]! += 1e-12))).toEqual({ mic: "current", selective: ["scores"] });
    expect(change((r) => (scheme(r, "random")[0]! += 1e-12))).toEqual({ mic: "current", selective: ["scores"] });
  });

  it("refuses research built from another snapshot or schema", () => {
    const other = raw<{ snapshot: { sha256: string } }>(researchWith([], []));
    other.snapshot.sha256 = "1".repeat(64);
    expect(() => parseResearch(other)).toThrow(/pins revision-5 snapshot 75a724c1/);
    expect(() => parseResearch({ ...copy(researchWith([], [])), schema: 2 })).toThrow(/research schema 1/);
  });

  it("gives a revision-6 report that keeps a revision-5 dataset id none of its numbers", () => {
    const rev5 = report();
    const research = ready(researchWith([micEntry(rev5)], [selectiveEntry(rev5, "A1")]));
    // Same id; two genomes replaced and every learner refitted.
    const rev6 = copy(rev5);
    rev6.samples[0]!.id = "573.900001";
    rev6.samples[1]!.id = "573.900002";
    rev6.splits.forEach((s) => s.learners.forEach((l) => (l.scores = l.scores.map((x) => x === null ? null : x / 2))));
    expect(researchFor(research, rev6)).toEqual({
      mic: { status: "mismatch", failed: ["sample_ids"], researchModel: true },
      selective: { status: "mismatch", failed: ["sample_ids", "scores"] },
    });
  });

  it("gives a report without an MIC model no MIC research, even when every guard part matches", () => {
    const mic = (entryFrom: DatasetReport, shown: DatasetReport) => researchFor(ready(researchWith([micEntry(entryFrom)], [])), shown).mic;
    const modeled = report();
    const without = (edit: (r: DatasetReport) => void = () => {}) => {
      const r = copy(modeled);
      r.mic = null;
      edit(r);
      return r;
    };
    // Neither the revision-5 report nor this one has a model, and the cohort
    // matches (K. pneumoniae ciprofloxacin on revisions 5 and 6).
    expect(mic(unmodeled(), unmodeled())).toEqual({ status: "no_model", failed: [] });
    // The revision-5 report had a model and this one has none (S. enterica
    // ampicillin on revision 6).
    expect(mic(modeled, without())).toEqual({ status: "no_model", failed: ["mic_predicted"] });
    // The cohort changed too (E. coli aztreonam on revision 6).
    expect(mic(modeled, without((r) => r.samples.pop()))).toEqual({
      status: "no_model",
      failed: ["samples", "sample_ids", "labels", "lineages", "mic_values", "mic_predicted"],
    });
    // A report with a model the revision-5 report lacked is a mismatch whose
    // research had no model.
    const gained = copy(unmodeled());
    gained.mic = copy(modeled.mic);
    gained.mic!.predicted = gained.samples.map(() => 1);
    expect(mic(unmodeled(), gained)).toEqual({ status: "mismatch", failed: ["mic_predicted"], researchModel: false });
  });

  it("withholds a zone-scored MIC model with its reason on every report revision", () => {
    const r = report();
    const research = ready(parseResearch(raw(researchWith([], [selectiveEntry(r, "A1")], [{ id: r.id, reason: WITHHELD_REASON }]))));
    expect(researchFor(research, r)).toMatchObject({ mic: { status: "withheld", reason: WITHHELD_REASON }, selective: { status: "current" } });
    // A later revision with the same id gets the reason too; selective
    // calling still checks its guard.
    const rev6 = copy(r);
    rev6.samples[0]!.id = "573.900001";
    expect(researchFor(research, rev6)).toEqual({ mic: { status: "withheld", reason: WITHHELD_REASON }, selective: { status: "mismatch", failed: ["sample_ids"] } });
    // The reason wins before any guard check, even over a matching entry.
    expect(researchFor(ready(researchWith([micEntry(r)], [], [{ id: r.id, reason: WITHHELD_REASON }])), r).mic).toEqual({ status: "withheld", reason: WITHHELD_REASON });
  });

  it("refuses a MIC entry that scores zone diameters and a withheld entry without a reason or with an entry", () => {
    const r = report();
    const scored = raw<{ mic: { datasets: { zone_diameters: { scored: number } }[] } }>(researchWith([micEntry(r, true)], []));
    scored.mic.datasets[0]!.zone_diameters.scored = 49;
    expect(() => parseResearch(scored)).toThrow(/scores 49 zone diameters as MICs, so it belongs in mic.withheld/);
    const withheld = (mic: MicDataset[], list: MicWithheld[]) => () => parseResearch(raw(researchWith(mic, [], list)));
    expect(withheld([], [{ id: r.id, reason: WITHHELD_REASON }])).not.toThrow();
    expect(withheld([micEntry(r)], [{ id: r.id, reason: WITHHELD_REASON }])).toThrow(/mic.withheld\[0\]: needs a reason and an id with no other MIC entry/);
    expect(withheld([], [{ id: r.id, reason: "" }])).toThrow(/mic.withheld\[0\]/);
    expect(withheld([], [{ id: r.id, reason: WITHHELD_REASON }, { id: r.id, reason: WITHHELD_REASON }])).toThrow(/mic.withheld\[1\]/);
    const missing = raw<{ mic: { withheld?: unknown } }>(researchWith([], []));
    delete missing.mic.withheld;
    expect(() => parseResearch(missing)).toThrow(/"withheld" should be an array/);
  });

  it("reads a loading research.json as loading and a failed one as absent", () => {
    const r = report();
    expect(researchFor({ status: "loading" }, r)).toEqual({ mic: { status: "loading" }, selective: { status: "loading" } });
    expect(researchFor({ status: "error", error: { kind: "missing", path: PATHS.research } }, r)).toEqual({ mic: { status: "absent" }, selective: { status: "absent" } });
  });

  it("keeps a saved MIC model with no lineage run or baseline usable", () => {
    const r = report();
    const state = researchFor(ready(parseResearch(raw(researchWith([micEntry(r, true)], [])))), r).mic;
    expect(state).toMatchObject({ status: "current", entry: { model: { source: "saved_report", lineage: null }, baseline: null } });
    const claimsBaseline = raw<{ mic: { datasets: Record<string, unknown>[] } }>(researchWith([micEntry(r, true)], []));
    claimsBaseline.mic.datasets[0]!.baseline = { random: micRun(0.7), lineage: micRun(0.5) };
    expect(() => parseResearch(claimsBaseline)).toThrow(/saved report model with no lineage run, no baseline/);
    const noExact = raw<{ mic: { datasets: { model: { random: Partial<MicRun> } }[] } }>(researchWith([micEntry(r, true)], []));
    delete noExact.mic.datasets[0]!.model.random.n_exact;
    expect(() => parseResearch(noExact)).toThrow(/model.random should be \{n_predicted, n_exact, ea_exact\}/);
  });

  it("takes essential agreement as a number or null and nothing else", () => {
    const r = report();
    const withEa = (ea: unknown) => {
      const doc = raw<{ mic: { datasets: { model: { random: Record<string, unknown> } }[] } }>(researchWith([micEntry(r)], []));
      if (ea === undefined) delete doc.mic.datasets[0]!.model.random.ea_exact;
      else doc.mic.datasets[0]!.model.random.ea_exact = ea;
      return () => parseResearch(doc);
    };
    expect(withEa(null)).not.toThrow();
    expect(withEa(0.8)).not.toThrow();
    expect(withEa("0.8")).toThrow(/model.random should be .* with ea_exact a number or null/);
    expect(withEa(undefined)).toThrow(/model.random should be .* with ea_exact a number or null/);
  });

  it("refuses a rerun whose baseline counts other isolates than its model", () => {
    const r = report();
    const edited = (edit: (d: MicDataset) => void) => {
      const entry = micEntry(r);
      edit(entry);
      return () => parseResearch(raw(researchWith([entry], [])));
    };
    expect(edited(() => {})).not.toThrow();
    expect(edited((d) => (d.baseline!.random.n_exact = 59))).toThrow(
      /mic.datasets\[0\]: on the random split the model scores 60 exact MICs of 80 predictions and the baseline 59 of 80, but both score the same isolates/,
    );
    expect(edited((d) => (d.model.lineage!.n_predicted = 79))).toThrow(/on the lineage split the model scores 60 exact MICs of 79 predictions and the baseline 60 of 80/);
  });

  it("refuses wording with no place for the removed zone diameters or the A2 learner", () => {
    const doc = () => raw<{ mic: { limits: { zone_fix: string } }; selective: { categories: { code: string; label: string }[] } }>(researchWith([], []));
    const noCount = doc();
    noCount.mic.limits.zone_fix = "Revision 5 read zone diameters as MICs.";
    expect(() => parseResearch(noCount)).toThrow(/zone_fix should hold \{removed\}/);
    const noLearner = doc();
    noLearner.selective.categories[1]!.label = "Category A2.";
    expect(() => parseResearch(noLearner)).toThrow(/only A2's holds \{learner\}/);
    const strayLearner = doc();
    strayLearner.selective.categories[2]!.label = "Category B with {learner}.";
    expect(() => parseResearch(strayLearner)).toThrow(/only A2's holds \{learner\}/);
  });

  it("lets only categories A1 and A2 carry numbers", () => {
    const r = report();
    SELECTIVE_CODES.forEach((c) => expect(() => parseResearch(raw(researchWith([], [selectiveEntry(r, c)])))).not.toThrow());
    const tampered = (category: SelectiveCategory, edit: (d: Record<string, unknown>) => void) => {
      const doc = raw<{ selective: { datasets: Record<string, unknown>[] } }>(researchWith([], [selectiveEntry(r, category)]));
      edit(doc.selective.datasets[0]!);
      return () => parseResearch(doc);
    };
    expect(tampered("B", (d) => (d.result = SELECTIVE_RESULT))).toThrow(/category B carries no numbers/);
    expect(tampered("C", (d) => (d.result = SELECTIVE_RESULT))).toThrow(/category C carries no numbers/);
    expect(tampered("A1", (d) => delete d.result)).toThrow(/category A1 needs its result/);
    expect(tampered("A2", (d) => delete (d.result as Record<string, unknown>).learner)).toThrow(/names its refit learner/);
    expect(tampered("A1", (d) => ((d.result as Record<string, unknown>).fn = 11))).toThrow(/do not partition/);
    const label = researchFor(ready(parseResearch(raw(researchWith([], [selectiveEntry(r, "B")])))), r).selective;
    expect(label).toMatchObject({ status: "current", entry: { category: "B" } });
    expect(label.status === "current" && "result" in label.entry).toBe(false);
    const reordered = raw<{ selective: { categories: unknown[] } }>(researchWith([], []));
    reordered.selective.categories.reverse();
    expect(() => parseResearch(reordered)).toThrow(/selective.categories does not list every category in order/);
  });

  it("loads the dataset when research.json is missing or invalid", async () => {
    const invalid = "/data/research-invalid.json";
    const json = (v: unknown) => new Response(JSON.stringify(v), { headers: { "Content-Type": "application/json" } });
    vi.stubGlobal("fetch", async (path: string) =>
      path === PATHS.research ? new Response("not found", { status: 404 }) : json(path === invalid ? { ...copy(researchWith([], [])), mic: [] } : fixtureFile(path)),
    );
    const seen: { research?: Resource<Research>; invalid?: Resource<Research>; dataset?: Resource<DatasetReport> } = {};
    const Probe = () => {
      seen.research = useResearch();
      seen.invalid = useResource(invalid, parseResearch);
      seen.dataset = useDataset(ID);
      return null;
    };
    const r = createRoot(document.createElement("div"));
    root = r;
    await act(async () => {
      r.render(createElement(Probe));
    });
    await act(async () => {
      await new Promise((done) => setTimeout(done, 0));
    });
    expect(seen.research).toEqual({ status: "error", error: { kind: "missing", path: PATHS.research } });
    expect(seen.invalid).toMatchObject({ status: "error", error: { kind: "shape" } });
    expect(seen.dataset?.status).toBe("ready");
    if (seen.dataset?.status === "ready") {
      const absent = { mic: { status: "absent" }, selective: { status: "absent" } };
      expect(researchFor(seen.research!, seen.dataset.data)).toEqual(absent);
      expect(researchFor(seen.invalid!, seen.dataset.data)).toEqual(absent);
    }
  });
});
