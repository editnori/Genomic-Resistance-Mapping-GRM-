import { describe, expect, it } from "vitest";
import { datasetTitle, hasEverydayName, speciesName, speciesOption, speciesShort } from "./format";

describe("species and dataset names", () => {
  it("gives a species its everyday name only where one is in wide use", () => {
    expect(speciesName("Staphylococcus aureus")).toBe("Staph");
    expect(speciesName("Mycobacterium tuberculosis")).toBe("TB");
    expect(speciesName("Neisseria gonorrhoeae")).toBe("Gonorrhea");
    expect(speciesName("Salmonella enterica")).toBe("Salmonella");
    // E. coli's short scientific name is its everyday name.
    expect(speciesName("Escherichia coli")).toBe(speciesShort("Escherichia coli"));
    expect(hasEverydayName("Escherichia coli")).toBe(false);
  });

  it("falls back to the short scientific name, and to the name itself, for any other species", () => {
    expect(speciesName("Klebsiella pneumoniae")).toBe("K. pneumoniae");
    expect(speciesName("Synthetic species B2")).toBe("S. species B2");
    expect(speciesName("synthetic_species_b2")).toBe("synthetic_species_b2");
    expect(hasEverydayName("Klebsiella pneumoniae")).toBe(false);
    // Only own entries count: an inherited property name is not a species.
    expect(hasEverydayName("constructor")).toBe(false);
    expect(speciesName("constructor")).toBe("constructor");
  });

  it("lists a species by its scientific name, with the everyday name after it", () => {
    expect(speciesOption("Staphylococcus aureus")).toBe("Staphylococcus aureus (Staph)");
    expect(speciesOption("Klebsiella pneumoniae")).toBe("Klebsiella pneumoniae");
  });

  it("titles a dataset species first, with vs held to the drug", () => {
    expect(datasetTitle({ species: "Salmonella enterica", antibiotic: "ciprofloxacin" })).toBe("Salmonella vs ciprofloxacin");
    expect(datasetTitle({ species: "Klebsiella pneumoniae", antibiotic: "meropenem" })).toBe("K. pneumoniae vs meropenem");
    expect(datasetTitle({ species: "Synthetic species C", antibiotic: "sulfamethoxazole/trimethoprim" })).toBe("S. species C vs sulfamethoxazole/trimethoprim");
  });
});
