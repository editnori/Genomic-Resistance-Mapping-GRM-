import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { auditOf } from "../data/derive/evaluation";
import { AuditChip, DatasetNotice, datasetTip } from "./EvaluationNotice";

const scope = "Internal cross-validation against reported lab calls; not external or clinical validation.";
const audit = (revision: number, current?: number) => auditOf({ evaluation: { revision, valid: true, reasons: [], scope } }, current);
const notice = (revision: number, current?: number) => renderToStaticMarkup(createElement(DatasetNotice, { audit: audit(revision, current) }));
const chip = (revision: number, current?: number) => renderToStaticMarkup(createElement(AuditChip, { audit: audit(revision, current) }));
const tipText = (revision: number, current?: number) =>
  datasetTip(audit(revision, current))
    .map((l) => [l.term, l.text, l.source].filter((s) => s !== undefined).join(" "))
    .join(" ");

describe("dataset audit notice", () => {
  it("states a passing audit once, as the header chip with its scope on hover, and marks an earlier revision provisional", () => {
    // A passing audit has no band; the chip carries it, and its hover the scope.
    expect(notice(6, 6)).toBe("");
    expect(chip(6, 6)).toContain("Internal cross-validation only");
    expect(tipText(6, 6)).toContain(scope);
    expect(tipText(6, 6)).toContain("passes the internal audit (revision 6)");

    // An earlier revision is provisional: a band with its reasons, a coral chip.
    const stale = notice(5, 6);
    expect(stale).toContain("Provisional");
    expect(stale).toContain("does not pass the current evaluation audit");
    expect(stale).toContain("current revision is 6");
    expect(stale).not.toContain("Internal cross-validation only");
    expect(stale).not.toContain(scope);
    expect(chip(5, 6)).toContain("Provisional");
    expect(chip(5, 6)).not.toContain("Internal cross-validation only");
  });
});
