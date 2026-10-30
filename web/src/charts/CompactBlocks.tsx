// A small population map. When a species' genomes make fewer than
// COMPACT_BLOCKS lineage blocks, the block columns keep a fixed width instead
// of stretching a handful of columns (or one) across the page, and this note
// and table say in words how many blocks there are, what each holds, and the
// map drug's calls in each. Nothing is added or split: the rows are the
// blocks lineageBlocks.ts makes, and their counts sum to the species totals.

import type { SpeciesReport } from "../data/schema";
import { capitalize, int, listWords, plural } from "../data/format";
import { BLOCK_MIN_GENOMES, blockSize, lineageSizes, rangeCalls, type CallCounts, type LineageBlock } from "./lineageBlocks";

/** Fewer blocks than this draw as a compact map. */
export const COMPACT_BLOCKS = 7;
/** The widest a compact map's column grows. */
export const COMPACT_BLOCK_PX = 72;

export type BlockLayout = "compact" | "spread";

export const blockLayout = (blocks: number): BlockLayout => (blocks > 0 && blocks < COMPACT_BLOCKS ? "compact" : "spread");

/** "126 genomes in 31 lineages make only 1 block: 1 pool of lineages under 20 genomes." */
export const compactBlocksText = (report: SpeciesReport, blocks: readonly LineageBlock[]): string => {
  const named = blocks.filter((b) => b.named).length;
  const pools = blocks.length - named;
  const parts = [
    ...(named > 0 ? [`${plural(named, "lineage")} of at least ${BLOCK_MIN_GENOMES} genomes`] : []),
    ...(pools > 0 ? [`${plural(pools, "pool")} of neighbouring lineages under ${BLOCK_MIN_GENOMES} genomes`] : []),
  ];
  const n = blocks.reduce((a, b) => a + blockSize(b), 0);
  return `${plural(n, "genome")} in ${plural(lineageSizes(report).size, "lineage")} make only ${plural(blocks.length, "block")}: ${listWords(parts)}. With this few, each column keeps a fixed width instead of stretching across the page.`;
};

const STATES = ["r", "i", "s", "untested"] as const;
const STATE_HEAD: Record<keyof CallCounts, [string, string]> = { r: ["R", "Resistant"], i: ["I", "Intermediate"], s: ["S", "Susceptible"], untested: ["–", "Not tested"] };

/** The note and the map drug's calls per block, for a compact map. */
export const CompactBlocks = ({ report, blocks, drug }: { report: SpeciesReport; blocks: readonly LineageBlock[]; drug: number }) => {
  const drugName = capitalize(report.drugs[drug] ?? "");
  const label = `${drugName} calls in each of the ${plural(blocks.length, "block")}`;
  return (
    <div data-compact-blocks={blocks.length} className="mt-3 space-y-3">
      <p className="prose-measure text-control text-ink">{compactBlocksText(report, blocks)}</p>
      <div className="scroll-thin max-w-full overflow-x-auto">
        <table aria-label={label} className="border-separate border-spacing-0 text-control tabular-nums">
          <thead>
            <tr className="text-tiny text-muted-text">
              <th scope="col" className="pb-1 pr-4 text-left font-medium">
                Block
              </th>
              <th scope="col" className="pb-1 pl-3 text-right font-medium">
                Genomes
              </th>
              {STATES.map((k) => (
                <th key={k} scope="col" className="pb-1 pl-3 text-right font-medium">
                  <span aria-hidden="true">{STATE_HEAD[k][0]}</span>
                  <span className="sr-only">{STATE_HEAD[k][1]}</span>
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {blocks.map((b) => {
              const c = rangeCalls(report, b.start, b.end, drug);
              return (
                <tr key={b.index} className="text-ink">
                  <th scope="row" className="border-t border-hair py-1 pr-4 text-left font-normal">
                    {`${b.index + 1}. ${b.lineages.length === 1 ? `Lineage ${b.lineages[0]}` : `${int(b.lineages.length)} lineages pooled`}`}
                  </th>
                  <td className="border-t border-hair py-1 pl-3 text-right">{int(blockSize(b))}</td>
                  {STATES.map((k) => (
                    <td key={k} className="border-t border-hair py-1 pl-3 text-right">
                      {int(c[k])}
                    </td>
                  ))}
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
      <p className="prose-measure text-note text-muted-text">{`${drugName} calls: R resistant, I intermediate, S susceptible, – not tested. Genomes shows each genome; Table lists every genome's calls.`}</p>
    </div>
  );
};
