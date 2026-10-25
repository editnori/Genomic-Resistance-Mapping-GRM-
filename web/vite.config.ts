import { createReadStream, statSync } from "node:fs";
import { resolve, sep } from "node:path";
import { fileURLToPath } from "node:url";
import type { Connect } from "vite";
import { defineConfig } from "vitest/config";
import react from "@vitejs/plugin-react";

// `grm report` writes the dashboard data under data/report/. In dev and
// preview, /data/* reads straight from it, the same mapping `grm serve` uses.
// GRM_REPORT_DIR points at a report built elsewhere (`grm report --data <root>`).
const REPORT_DIR = resolve(process.env.GRM_REPORT_DIR ?? resolve(fileURLToPath(new URL(".", import.meta.url)), "../data/report"));

const serveReport: Connect.NextHandleFunction = (req, res, next) => {
  const url = decodeURIComponent((req.url ?? "/").split("?")[0] ?? "/");
  const file = resolve(REPORT_DIR, `.${url}`);
  // Refuse any path that resolves outside the report directory.
  if (!file.startsWith(REPORT_DIR + sep)) {
    res.statusCode = 404;
    res.end("not found");
    return;
  }
  try {
    if (!statSync(file).isFile()) throw new Error("not a file");
  } catch {
    res.statusCode = 404;
    res.end("not found");
    return;
  }
  res.setHeader("Content-Type", file.endsWith(".json") ? "application/json" : "application/octet-stream");
  res.setHeader("Cache-Control", "no-store");
  createReadStream(file).on("error", next).pipe(res);
};

export default defineConfig({
  plugins: [
    react(),
    {
      name: "grm-report-data",
      configureServer(server) {
        server.middlewares.use("/data", serveReport);
      },
      configurePreviewServer(server) {
        server.middlewares.use("/data", serveReport);
      },
    },
  ],
  build: { outDir: "dist", emptyOutDir: true },
  test: { include: ["src/**/*.test.ts"], environment: "node" },
});
