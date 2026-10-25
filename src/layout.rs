//! Where every artifact lives under the data root.

use std::path::{Path, PathBuf};

#[derive(Clone, Debug)]
pub struct Layout {
    pub root: PathBuf,
}

impl Layout {
    pub fn new(root: impl Into<PathBuf>) -> Self {
        Layout { root: root.into() }
    }

    pub fn amr(&self) -> PathBuf {
        self.root.join("bvbrc/genome_amr.jsonl")
    }

    pub fn genomes(&self) -> PathBuf {
        self.root.join("bvbrc/genomes.jsonl")
    }

    pub fn fasta_dir(&self) -> PathBuf {
        self.root.join("fasta")
    }

    pub fn catalog_dir(&self) -> PathBuf {
        self.root.join("catalog")
    }

    pub fn curated(&self) -> PathBuf {
        self.root.join("curated")
    }

    pub fn species_dir(&self, slug: &str) -> PathBuf {
        self.root.join("species").join(slug)
    }

    pub fn kset_dir(&self, slug: &str) -> PathBuf {
        self.species_dir(slug).join("ksets")
    }

    pub fn matrix_dir(&self, slug: &str) -> PathBuf {
        self.species_dir(slug).join("matrix")
    }

    pub fn results_dir(&self, slug: &str) -> PathBuf {
        self.species_dir(slug).join("results")
    }

    pub fn features_dir(&self) -> PathBuf {
        self.root.join("bvbrc/features")
    }
}

pub fn ensure_parent(path: &Path) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    Ok(())
}

/// Removes a file and says whether it was there; one that is already gone
/// is not an error.
pub fn remove_if_present(path: &Path) -> std::io::Result<bool> {
    match std::fs::remove_file(path) {
        Ok(()) => Ok(true),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(false),
        Err(e) => Err(e),
    }
}
