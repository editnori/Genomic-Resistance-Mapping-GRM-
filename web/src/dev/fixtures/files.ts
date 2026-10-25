// The test fixture by /data path, round-tripped through JSON so the edge
// parsers see what a fetch would deliver.

import { buildFixture, type FixtureSet } from "./generate";

let built: FixtureSet | null = null;

export const fixtureFile = (path: string): unknown => {
  built ??= buildFixture();
  const name = decodeURIComponent(path.replace(/^\/data\//, "").replace(/\.json$/, ""));
  const value =
    name === "index"
      ? built.manifest
      : name === "corpus"
        ? built.corpus
        : name.startsWith("species/")
          ? built.species.get(name.slice("species/".length))
          : name.startsWith("datasets/")
            ? built.datasets.get(name.slice("datasets/".length))
            : undefined;
  return value === undefined ? null : (JSON.parse(JSON.stringify(value)) as unknown);
};
