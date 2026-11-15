// The rule words of a published report, when GRM_REPORT_DIR names a report
// directory (the variable vite.config.ts reads): every rule term of every
// dataset, worded as the index words it, with the app's family note and the
// within-lineage verdict the inspector would attach.
//
// The second test is the cross-check of the two writers on real data. The
// index's words are written by `report::rule_text` (src/report.rs) and built
// again here by `indexRuleWords`. With GRM_RULE_WORDS_RUST naming the file the
// ignored Rust test `report::tests::the_rule_words_of_published_reports`
// wrote for the same directory (its GRM_RULE_WORDS_OUT), this suite compares
// the two itself, term by term. Without that file the test is reported as
// skipped: it never passes without making the comparison. docs/ci.md gives
// the two commands.
//
// With GRM_RULE_WORDS_OUT set, the app's words are also written there, one
// `dataset<TAB>position<TAB>words` line a term. Skipped without
// GRM_REPORT_DIR, so the fixture-only suite stays self-contained.

import { describe, expect, it } from "vitest";
import { parseDataset, parseManifest } from "../validate";
import { familiesNamed, featureIndex, indexRuleWords, RULE_TEXT_VERSION } from "./dataset";
import { LINEAGE_CONTROL, lineageVerdicts } from "./lineage";

const env = (globalThis as { process?: { env: Record<string, string | undefined> } }).process?.env ?? {};
const dir = env.GRM_REPORT_DIR;
const rustWords = env.GRM_RULE_WORDS_RUST;

// A computed specifier keeps Node's types out of the app's compile.
const files = async () => (await import(/* @vite-ignore */ ["node", "fs/promises"].join(":"))) as { readFile(p: string, encoding: string): Promise<string>; writeFile(p: string, text: string): Promise<void> };

/** Every rule term of the directory's reports in the app's words, as the
 *  lines the Rust test writes: datasets in the order of their file names
 *  (`<id>.json`), and terms in the model's order. */
const appLines = async (): Promise<string[]> => {
  const fs = await files();
  const read = async (path: string): Promise<unknown> => JSON.parse(await fs.readFile(`${dir}/${path}`, "utf8")) as unknown;
  const manifest = parseManifest(await read("index.json"));
  const byFile = manifest.datasets.map((row) => ({ id: row.id, file: `${row.id}.json` })).sort((a, b) => (a.file < b.file ? -1 : a.file > b.file ? 1 : 0));
  const lines: string[] = [];
  for (const { id } of byFile) {
    const r = parseDataset(await read(`datasets/${id}.json`));
    r.models.scm.rules.forEach((_, i) => lines.push(`${id}\t${i + 1}\t${indexRuleWords(r, i)}`));
  }
  return lines;
};

describe.skipIf(dir === undefined)("published rule words", () => {
  it("word every rule term as the index does, and attach a verdict only to the release the table was computed on", async () => {
    const fs = await files();
    const read = async (path: string): Promise<unknown> => JSON.parse(await fs.readFile(`${dir}/${path}`, "utf8")) as unknown;
    const manifest = parseManifest(await read("index.json"));
    const lines: string[] = [];
    let verdicts = 0;
    for (const row of manifest.datasets) {
      const r = parseDataset(await read(`datasets/${row.id}.json`));
      const cards = featureIndex(r);
      const words = r.models.scm.rules.map((_, i) => indexRuleWords(r, i));
      expect([row.id, row.scm_rules.length]).toEqual([row.id, words.length]);
      // An index of this wording holds exactly these words. The served index
      // of an older wording does not, so this line alone is no cross-check
      // of the two writers: the next test is.
      if (manifest.rule_text_version === RULE_TEXT_VERSION) expect([row.id, row.scm_rules]).toEqual([row.id, words]);
      r.models.scm.rules.forEach((rule, i) => {
        const families = new Set(cards.get(rule.feature)?.catalog_families ?? []).size;
        expect([row.id, i, familiesNamed(words[i]!)]).toEqual([row.id, i, families > 1 ? families : null]);
        lines.push(`${row.id}\t${i + 1}\t${words[i]}`);
      });
      const entry = Object.hasOwn(LINEAGE_CONTROL.datasets, row.id) ? LINEAGE_CONTROL.datasets[row.id] : undefined;
      const same = entry !== undefined && r.evaluation?.revision === LINEAGE_CONTROL.evaluation_revision && r.input_fingerprint === entry.fingerprint;
      const got = lineageVerdicts(r);
      if (same) {
        expect([row.id, got?.length]).toEqual([row.id, words.length]);
        verdicts += words.length;
      } else {
        expect([row.id, got]).toEqual([row.id, null]);
      }
    }
    // Every term of the table's own release has a verdict; another release has none.
    expect([0, LINEAGE_CONTROL.terms]).toContain(verdicts);
    if (env.GRM_RULE_WORDS_OUT !== undefined) await fs.writeFile(env.GRM_RULE_WORDS_OUT, `${lines.join("\n")}\n`);
  }, 600_000);

  it.skipIf(rustWords === undefined)("are the words the Rust index writer gives the same reports, term by term (GRM_RULE_WORDS_RUST)", async () => {
    const fs = await files();
    const rust = (await fs.readFile(rustWords!, "utf8")).split("\n");
    // The file ends with a newline, and holds nothing but the term lines.
    expect(rust.pop()).toBe("");
    const app = await appLines();
    expect(app.length).toBeGreaterThan(0);
    expect(rust.length).toBe(app.length);
    // Term by term first, so a failure names the dataset and the two wordings.
    app.forEach((line, i) => expect(rust[i]).toBe(line));
    expect(rust).toEqual(app);
  }, 600_000);
});
