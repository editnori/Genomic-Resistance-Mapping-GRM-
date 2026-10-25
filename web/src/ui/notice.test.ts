import { createElement } from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { auditOf } from "../data/derive/evaluation";
import { DatasetNotice } from "./EvaluationNotice";

const scope = "Internal cross-validation against reported lab calls; not external or clinical validation.";
const notice = (revision: number, current?: number) =>
  renderToStaticMarkup(createElement(DatasetNotice, { audit: auditOf({ evaluation: { revision, valid: true, reasons: [], scope } }, current) }));

describe("dataset audit notice", () => {
  it("leaves passing audit details in Report and marks an earlier revision provisional", () => {
    const current = notice(6, 6);
    const stale = notice(5, 6);

    expect(current).toBe("");
    expect(stale).toContain("Provisional");
    expect(stale).toContain("does not pass the current evaluation audit");
    expect(stale).toContain("current revision is 6");
    expect(stale).not.toContain("Internal cross-validation only");
    expect(stale).not.toContain(scope);
  });
});
