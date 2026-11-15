// The clindamycin label check on the site (docs/clindamycin-labels.md): the
// note beside a dataset report's learner comparison, the short mark in the
// datasets table, and the Methods paragraph. Each shows only for the release
// the check read (data/derive/clindamycin.ts); otherwise it renders nothing,
// and the view is drawn as before.

import type { CatalogRow, DatasetReport, Manifest } from "../../data/schema";
import { CLINDAMYCIN_LABELS, clindamycinMarkTip, clindamycinMethods, clindamycinNote, type ClindamycinLabels } from "../../data/derive/clindamycin";

const PAGE_CODE = "font-mono text-[0.92em]";

/** The page, named by its path as the app names a repository document. */
const PageLink = ({ page }: { page: string }) => (
  <>
    The check and its limits are in <code className={PAGE_CODE}>{page}</code> in the GRM repository.
  </>
);

/** One note under the learner comparison of each of the six datasets: the
 *  kind, the counts behind it, what the figure can be read as, and the page.
 *  A caller that passes `digest`, the loaded report's sha256 (undefined
 *  until it is read, null when none exists), gets the note only when it is
 *  the file's report_sha256; a caller that passes none gets the pin by
 *  dataset, revision and sample count. */
export const ClindamycinNote = (props: { r: DatasetReport; file?: ClindamycinLabels | null; digest?: string | null | undefined }) => {
  const { r, file = CLINDAMYCIN_LABELS } = props;
  const note = clindamycinNote(r, file, Object.hasOwn(props, "digest") ? { digest: props.digest } : undefined);
  if (note === null) return null;
  return (
    <div role="note" data-clindamycin-note={note.kind} className="border-y border-hair py-3">
      <p className="prose-measure text-body text-ink">
        <strong className="font-semibold">{note.title}.</strong> {note.counts.join(" ")}
      </p>
      <p className="prose-measure type-note mt-1.5">
        {note.reading} <PageLink page={note.page} />
      </p>
    </div>
  );
};

/** The datasets table's mark beside the drug name, with the kind in its
 *  native tooltip and for a screen reader. */
export const ClindamycinMark = ({ row, file = CLINDAMYCIN_LABELS }: { row: CatalogRow; file?: ClindamycinLabels | null }) => {
  const tip = clindamycinMarkTip(row, file);
  if (tip === null) return null;
  return (
    <span data-clindamycin-mark="" title={tip} className="shrink-0 cursor-help text-tiny font-medium text-muted-text underline decoration-dotted underline-offset-2">
      Label check<span className="sr-only">. {tip}</span>
    </span>
  );
};

/** The Methods paragraph on the check, at the definitions' size. */
export const ClindamycinMethods = ({ manifest, file = CLINDAMYCIN_LABELS }: { manifest: Pick<Manifest, "evaluation_revision">; file?: ClindamycinLabels | null }) => {
  const text = clindamycinMethods(manifest, file);
  if (text === null || file === null) return null;
  return (
    <p data-clindamycin-methods="" className="prose-measure text-note text-ink">
      {text} <PageLink page={file.page} />
    </p>
  );
};
