// The rule words of a published report, when GRM_REPORT_DIR names a report
// directory (the variable vite.config.ts reads): every rule term of every
// dataset, worded as the index words it, with the app's family note and the
// within-lineage verdict the inspector would attach. With GRM_RULE_WORDS_OUT
// set, the words are also written there, one `dataset<TAB>position<TAB>words`
// line a term, the lines the ignored Rust test
// `report::tests::the_rule_words_of_published_reports` prints, so the two can
// be compared. Skipped without GRM_REPORT_DIR, so the fixture-only suite
// stays self-contained.

import { describe, expect, it } from "vitest";
import { parseDataset, parseManifest } from "../validate";
import { familiesNamed, featureIndex, indexRuleWords, RULE_TEXT_VERSION } from "./dataset";
import { LINEAGE_CONTROL, lineageVerdicts } from "./lineage";

const env = (globalThis as { process?: { env: Record<string, string | undefined> } }).process?.env ?? {};
const dir = env.GRM_REPORT_DIR;

// A computed specifier keeps Node's types out of the app's compile.
const files = async () => (await import(/* @vite-ignore */ ["node", "fs/promises"].join(":"))) as { readFile(p: string, encoding: string): Promise<string>; writeFile(p: string, text: string): Promise<void> };

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
      // An index of this wording holds exactly these words.
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
});
