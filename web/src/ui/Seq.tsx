import { contextRuns } from "../data/derive/dataset";
import { int } from "../data/format";

/** Contig context with the k-mer span in capitals and tinted, so the span
 *  reads by case and by tint; a long span shows its ends and a base count. */
export const Seq = ({ context, className }: { context: string; className?: string }) => (
  <p className={`seq rounded-md bg-surface-subtle px-2.5 py-2 text-control leading-relaxed text-muted-text ${className ?? ""}`}>
    {contextRuns(context).map((run, i) =>
      run.elided !== undefined ? (
        <span key={i} className="mx-1 whitespace-nowrap text-tiny font-medium tracking-normal text-ink">
          … {int(run.elided)} more bases …
        </span>
      ) : run.span ? (
        <mark key={i} className="rounded-sm bg-blue/15 px-px font-semibold text-ink">
          {run.text}
        </mark>
      ) : (
        <span key={i}>{run.text}</span>
      ),
    )}
  </p>
);
