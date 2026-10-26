import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { auditOf } from "../data/derive/evaluation";
import { parseManifest } from "../data/validate";
import { fixtureFile } from "../dev/fixtures/files";
import { LIMITS_HREF } from "../app/router";
import { ReportStatus } from "../app/ReportStatus";
import { AuditChip, DatasetNotice, datasetTip } from "./EvaluationNotice";

const scope = "Internal cross-validation against reported lab calls; not external or clinical validation.";
const audit = (revision: number, current?: number) => auditOf({ evaluation: { revision, valid: true, reasons: [], scope } }, current);
const notice = (revision: number, current?: number) => renderToStaticMarkup(createElement(DatasetNotice, { audit: audit(revision, current) }));
const chip = (revision: number, current?: number) => renderToStaticMarkup(createElement(AuditChip, { audit: audit(revision, current) }));
const tipText = (revision: number, current?: number) =>
  datasetTip(audit(revision, current))
    .map((l) => [l.term, l.text, l.source].filter((s) => s !== undefined).join(" "))
    .join(" ");
const text = (html: string) => html.replace(/<[^>]+>/g, "");

describe("dataset audit notice", () => {
  it("states a passing audit once, as plain header words, and marks an earlier revision provisional", () => {
    // A passing audit has no band; the header says it in plain words, and
    // its scope stays available to the view that prints it.
    expect(notice(6, 6)).toBe("");
    expect(chip(6, 6)).toContain("Internal cross-validation only");
    expect(chip(6, 6)).not.toContain('role="tooltip"');
    expect(tipText(6, 6)).toContain(scope);
    expect(tipText(6, 6)).toContain("passes the internal audit (revision 6)");

    // An earlier revision is provisional: a band with its reasons, bold words.
    const stale = notice(5, 6);
    expect(stale).toContain("Provisional");
    expect(stale).toContain("does not pass the current evaluation audit");
    expect(stale).toContain("current revision is 6");
    expect(stale).not.toContain("Internal cross-validation only");
    expect(stale).not.toContain(scope);
    expect(chip(5, 6)).toMatch(/<strong[^>]*>Provisional<\/strong>/);
    expect(chip(5, 6)).not.toContain("Internal cross-validation only");
  });
});

describe("report status line", () => {
  const fixture = parseManifest(fixtureFile("/data/index.json"));
  // Every evaluation passes at the manifest's revision.
  const manifest = { ...fixture, evaluation_revision: 7, datasets: fixture.datasets.map((d) => ({ ...d, evaluation: { revision: 7, valid: true, reasons: [], scope } })) };
  const line = (m = manifest) => renderToStaticMarkup(createElement(ReportStatus, { manifest: m }));

  it("prints the build date, the catalog version and the audit scope as plain words linked to the limits", () => {
    const html = line();
    expect(text(html)).toContain("Report built");
    expect(text(html)).toContain(`AMRFinderPlus ${manifest.amrfinder_version}`);
    expect(html).toMatch(new RegExp(`<a[^>]*href="${LIMITS_HREF}"[^>]*>Internal cross-validation only</a>`));
    expect(html).not.toContain('role="tooltip"');
    expect(text(html)).not.toContain("BV-BRC lab records");
  });

  it("says Provisional metrics in bold when an evaluation fails the audit", () => {
    const html = line({ ...manifest, evaluation_revision: 8 });
    expect(html).toMatch(new RegExp(`<a[^>]*href="${LIMITS_HREF}"[^>]*><strong[^>]*>Provisional metrics</strong></a>`));
  });

  it("ends each fact with its separator, so a wrapped line never starts with one", () => {
    const html = line();
    const separators = html.match(/<span aria-hidden="true">·<\/span>/g) ?? [];
    expect(separators).toHaveLength(2);
    expect(html.match(/<span aria-hidden="true">·<\/span><\/span>/g)).toHaveLength(2);
  });
});
