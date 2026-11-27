// Tiny zero-data slug list shared by the router and flow preview so the main
// bundle never imports flow.json until /flow is opened.

export type FlowSlug = "genome-table" | "genome-path" | "bvbrc" | "piece-count" | "decision-tree";

export const FLOW_SLUGS: readonly { slug: FlowSlug; viewKey: "genome_table" | "genome_path" | "bvbrc" | "piece_count" | "decisions" }[] = [
  { slug: "genome-table", viewKey: "genome_table" },
  { slug: "genome-path", viewKey: "genome_path" },
  { slug: "bvbrc", viewKey: "bvbrc" },
  { slug: "piece-count", viewKey: "piece_count" },
  { slug: "decision-tree", viewKey: "decisions" },
] as const;
