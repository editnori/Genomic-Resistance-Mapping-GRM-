// Early-return body cancellation in fetchBytes. Every response here is
// synthetic so unit tests run offline without touching network or BV-BRC data.
import { afterEach, describe, expect, it, vi } from "vitest";
import { fetchBytes } from "./load";

function trackedResponse(bodyText: string, init: ResponseInit): {
  response: Response;
  isOpen: () => boolean;
  cancelCount: () => number;
} {
  let open = true;
  let cancels = 0;
  const encoded = new TextEncoder().encode(bodyText);
  const stream = new ReadableStream<Uint8Array>({
    pull(controller) {
      controller.enqueue(encoded);
      controller.close();
      open = false;
    },
    cancel() {
      cancels += 1;
      open = false;
    },
  });
  return {
    response: new Response(stream, init),
    isOpen: () => open,
    cancelCount: () => cancels,
  };
}

describe("fetchBytes", () => {
  afterEach(() => {
    vi.restoreAllMocks();
  });

  it("cancels the response body on 404 so no request stream stays open", async () => {
    const tracked = trackedResponse("Not Found", {
      status: 404,
      headers: { "Content-Type": "text/plain" },
    });
    vi.spyOn(globalThis, "fetch").mockResolvedValue(tracked.response);

    const result = await fetchBytes("/data/baselines.json");
    await Promise.resolve();

    expect(result).toEqual({ ok: false, error: { kind: "missing", path: "/data/baselines.json" } });
    expect(tracked.cancelCount()).toBe(1);
    expect(tracked.isOpen()).toBe(false);
  });

  it("cancels the response body on non-404 HTTP error statuses", async () => {
    const tracked = trackedResponse("Server Error", {
      status: 503,
      headers: { "Content-Type": "text/plain" },
    });
    vi.spyOn(globalThis, "fetch").mockResolvedValue(tracked.response);

    const result = await fetchBytes("/data/index.json");
    await Promise.resolve();

    expect(result).toEqual({ ok: false, error: { kind: "http", path: "/data/index.json", status: 503 } });
    expect(tracked.cancelCount()).toBe(1);
    expect(tracked.isOpen()).toBe(false);
  });

  it("cancels the response body on a 200 text/html fallback", async () => {
    const tracked = trackedResponse("<!doctype html><html></html>", {
      status: 200,
      headers: { "Content-Type": "text/html; charset=utf-8" },
    });
    vi.spyOn(globalThis, "fetch").mockResolvedValue(tracked.response);

    const result = await fetchBytes("/data/research.json");
    await Promise.resolve();

    expect(result).toEqual({ ok: false, error: { kind: "missing", path: "/data/research.json" } });
    expect(tracked.cancelCount()).toBe(1);
    expect(tracked.isOpen()).toBe(false);
  });

  it("reads JSON bytes without cancelling on 200 application/json", async () => {
    const tracked = trackedResponse('{"ok":true}', {
      status: 200,
      headers: { "Content-Type": "application/json" },
    });
    vi.spyOn(globalThis, "fetch").mockResolvedValue(tracked.response);

    const result = await fetchBytes("/data/index.json");

    expect(result.ok).toBe(true);
    if (result.ok) {
      expect(new TextDecoder().decode(result.bytes)).toBe('{"ok":true}');
    }
    expect(tracked.cancelCount()).toBe(0);
    expect(tracked.isOpen()).toBe(false);
  });
});
