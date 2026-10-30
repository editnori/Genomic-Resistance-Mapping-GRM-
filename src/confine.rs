//! Writes that stay inside the data root.
//!
//! A data root is often assembled from links: FASTA files, the catalog and
//! the BV-BRC records may point into shared, frozen stores. Those are read
//! through their links. What GRM writes, it writes only below the root
//! itself, and only through real directories:
//!
//! - `Lane::check_dir` and `Lane::create_dir_all` walk a directory from the
//!   root one component at a time and refuse any component that is a
//!   symlink or not a directory; the directory's canonical path must then be
//!   its place below the canonical root.
//! - `create_part` makes a temporary with `O_CREAT | O_EXCL | O_NOFOLLOW`,
//!   so it never opens an existing file or follows a link; a stale regular
//!   temporary is unlinked first, anything else is refused.
//! - `replace_with` renames a finished temporary over its target only when
//!   the target is absent or a regular file. A rename replaces the name, not
//!   what a link points to, but a link where a GRM file belongs is refused
//!   rather than silently replaced.
//!
//! The data root itself may be a link; it is resolved once.

use anyhow::{Context, Result, bail, ensure};
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::os::unix::fs::OpenOptionsExt;
use std::path::{Component, Path, PathBuf};

/// The data root as given and as resolved.
#[derive(Clone, Debug)]
pub struct Lane {
    given: PathBuf,
    canonical: PathBuf,
}

impl Lane {
    /// The root must exist and be a directory.
    pub fn new(root: &Path) -> Result<Lane> {
        let canonical = std::fs::canonicalize(root)
            .with_context(|| format!("resolve the data root {}", root.display()))?;
        ensure!(
            canonical.is_dir(),
            "the data root {} is not a directory",
            root.display()
        );
        Ok(Lane {
            given: root.to_path_buf(),
            canonical,
        })
    }

    /// `dir` below the root, as plain components.
    fn relative<'a>(&self, dir: &'a Path) -> Result<&'a Path> {
        let rel = dir.strip_prefix(&self.given).with_context(|| {
            format!(
                "{} is not below the data root {}",
                dir.display(),
                self.given.display()
            )
        })?;
        ensure!(
            rel.components().all(|c| matches!(c, Component::Normal(_))),
            "{} is not a plain path below the data root",
            dir.display()
        );
        Ok(rel)
    }

    /// Walks `dir` from the root. A missing component is created when
    /// `create` is set (syncing its parent when `durable`), and otherwise
    /// ends the walk with Ok(false).
    fn walk(&self, dir: &Path, create: bool, durable: bool) -> Result<bool> {
        let rel = self.relative(dir)?;
        let mut at = self.canonical.clone();
        for component in rel.components() {
            at.push(component);
            match std::fs::symlink_metadata(&at) {
                Ok(m) if m.file_type().is_symlink() => bail!(
                    "refusing symlinked directory {} inside the data root {}: GRM writes only through real directories below it",
                    at.display(),
                    self.canonical.display()
                ),
                Ok(m) if m.is_dir() => {}
                Ok(_) => bail!("{} is not a directory", at.display()),
                Err(e) if e.kind() == std::io::ErrorKind::NotFound && create => {
                    // mkdir never follows a link at its last component: one
                    // planted since the lstat makes it fail with EEXIST.
                    std::fs::create_dir(&at).with_context(|| format!("create {}", at.display()))?;
                    if durable {
                        sync_dir(at.parent().unwrap_or(Path::new("/")))?;
                    }
                }
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(false),
                Err(e) => return Err(e).with_context(|| format!("inspect {}", at.display())),
            }
        }
        let resolved =
            std::fs::canonicalize(dir).with_context(|| format!("resolve {}", dir.display()))?;
        ensure!(
            resolved == at,
            "{} resolves to {}, not to its place {} below the data root",
            dir.display(),
            resolved.display(),
            at.display()
        );
        Ok(true)
    }

    /// `dir` exists below the root, through real directories only.
    pub fn check_dir(&self, dir: &Path) -> Result<()> {
        ensure!(
            self.walk(dir, false, false)?,
            "{} does not exist",
            dir.display()
        );
        Ok(())
    }

    /// As `check_dir`, but a missing `dir` (or ancestor) is Ok(false).
    pub fn check_dir_if_present(&self, dir: &Path) -> Result<bool> {
        self.walk(dir, false, false)
    }

    /// Creates `dir` and any missing ancestor below the root, refusing links
    /// on the way. `durable` syncs the parent of each directory it makes.
    pub fn create_dir_all(&self, dir: &Path, durable: bool) -> Result<()> {
        self.walk(dir, true, durable).map(drop)
    }

    /// `path`'s directory passes `check_dir_if_present`, and `path` itself is
    /// absent or a regular file. Ok(true) when it is a regular file.
    pub fn check_file(&self, path: &Path) -> Result<bool> {
        let parent = path
            .parent()
            .with_context(|| format!("{} has no directory", path.display()))?;
        if !self.check_dir_if_present(parent)? {
            return Ok(false);
        }
        regular_or_absent(path)
    }
}

/// Ok(true) for a regular file, Ok(false) when nothing is there, and an
/// error for a symlink, directory or anything else (never followed).
pub fn regular_or_absent(path: &Path) -> Result<bool> {
    match std::fs::symlink_metadata(path) {
        Ok(m) if m.file_type().is_file() => Ok(true),
        Ok(m) if m.file_type().is_symlink() => bail!(
            "refusing to write through symlink {}: a file GRM writes must be a regular file",
            path.display()
        ),
        Ok(_) => bail!("refusing to replace non-regular-file {}", path.display()),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(false),
        Err(e) => Err(e).with_context(|| format!("inspect {}", path.display())),
    }
}

/// Creates `path` for writing; fails if anything, a link included, is there.
pub fn create_new(path: &Path) -> Result<File> {
    OpenOptions::new()
        .write(true)
        .create_new(true)
        .custom_flags(libc::O_NOFOLLOW)
        .open(path)
        .with_context(|| format!("create {}", path.display()))
}

/// A fresh temporary at `part`: a stale regular file left by an interrupted
/// write is removed first; a link, directory or other entry is refused.
pub fn create_part(part: &Path) -> Result<File> {
    if regular_or_absent(part)? {
        std::fs::remove_file(part).with_context(|| format!("remove {}", part.display()))?;
    }
    create_new(part)
}

/// Renames the finished temporary `part` to `target`, which must be absent
/// or a regular file.
pub fn replace_with(part: &Path, target: &Path) -> Result<()> {
    regular_or_absent(target)?;
    std::fs::rename(part, target).with_context(|| format!("rename to {}", target.display()))
}

/// Makes renames, creations and removals in `dir` durable.
pub fn sync_dir(dir: &Path) -> Result<()> {
    let dir = if dir.as_os_str().is_empty() {
        Path::new(".")
    } else {
        dir
    };
    File::open(dir)
        .and_then(|d| d.sync_all())
        .with_context(|| format!("sync directory {}", dir.display()))
}

/// `path` with `suffix` appended to its file name.
pub fn suffixed(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

/// Writes `bytes` to `path` through `<path>.part` and a rename, so a reader
/// never sees a torn file and nothing outside the name is written. With
/// `durable`, the file and its directory are synced before this returns.
/// A failure removes the temporary.
pub fn write_atomic(path: &Path, bytes: &[u8], durable: bool) -> Result<()> {
    let part = suffixed(path, ".part");
    let mut file = create_part(&part)?;
    let result = (|| {
        file.write_all(bytes)
            .with_context(|| format!("write {}", part.display()))?;
        if durable {
            file.sync_all()
                .with_context(|| format!("sync {}", part.display()))?;
        }
        replace_with(&part, path)?;
        if durable {
            sync_dir(path.parent().unwrap_or(Path::new(".")))?;
        }
        Ok(())
    })();
    drop(file);
    if result.is_err() && regular_or_absent(&part).unwrap_or(false) {
        let _ = std::fs::remove_file(&part);
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::os::unix::fs::symlink;

    #[test]
    fn writes_never_follow_or_replace_links() {
        let dir = tempfile::tempdir().unwrap();
        let victim = dir.path().join("victim");
        std::fs::write(&victim, "outside").unwrap();
        let target = dir.path().join("target.json");
        symlink(&victim, &target).unwrap();
        let err = write_atomic(&target, b"lane", false).unwrap_err();
        assert!(format!("{err:#}").contains("symlink"), "{err:#}");
        assert_eq!(std::fs::read(&victim).unwrap(), b"outside");
        assert!(target.is_symlink());
        assert!(!suffixed(&target, ".part").exists());

        // A link at the temporary's name is refused, not followed or removed.
        let plain = dir.path().join("plain.json");
        symlink(&victim, suffixed(&plain, ".part")).unwrap();
        let err = write_atomic(&plain, b"lane", false).unwrap_err();
        assert!(format!("{err:#}").contains("symlink"), "{err:#}");
        assert_eq!(std::fs::read(&victim).unwrap(), b"outside");
        assert!(suffixed(&plain, ".part").is_symlink());
        assert!(create_new(&suffixed(&plain, ".part")).is_err());

        // A stale regular temporary is replaced, and the write lands.
        std::fs::remove_file(suffixed(&plain, ".part")).unwrap();
        std::fs::write(suffixed(&plain, ".part"), "stale").unwrap();
        write_atomic(&plain, b"lane", true).unwrap();
        assert_eq!(std::fs::read(&plain).unwrap(), b"lane");
        assert!(!suffixed(&plain, ".part").exists());
    }

    #[test]
    fn lane_directories_are_real_and_stay_below_the_root() {
        let dir = tempfile::tempdir().unwrap();
        let root = dir.path().join("data");
        let outside = dir.path().join("outside");
        std::fs::create_dir_all(&root).unwrap();
        std::fs::create_dir_all(&outside).unwrap();
        let lane = Lane::new(&root).unwrap();
        let deep = root.join("species/x/matrix");
        assert!(!lane.check_dir_if_present(&deep).unwrap());
        lane.create_dir_all(&deep, true).unwrap();
        lane.check_dir(&deep).unwrap();

        // A link anywhere below the root is refused, before and after it.
        symlink(&outside, root.join("species/x/ksets")).unwrap();
        for path in [root.join("species/x/ksets"), root.join("species/x/ksets/a")] {
            let err = lane.create_dir_all(&path, false).unwrap_err();
            assert!(
                format!("{err:#}").contains("symlinked directory"),
                "{err:#}"
            );
            assert!(lane.check_dir_if_present(&path).is_err());
        }
        assert_eq!(std::fs::read_dir(&outside).unwrap().count(), 0);
        let file = root.join("species/x/ksets/g.json");
        assert!(lane.check_file(&file).is_err());
        assert!(lane.check_dir(&outside).is_err());
        assert!(lane.check_dir(&root.join("species/../species")).is_err());

        // The root itself may be a link; it is resolved once.
        let alias = dir.path().join("alias");
        symlink(&root, &alias).unwrap();
        let lane = Lane::new(&alias).unwrap();
        lane.check_dir(&alias.join("species/x/matrix")).unwrap();
        assert!(
            !lane
                .check_file(&alias.join("species/x/matrix/m.bin"))
                .unwrap()
        );
    }
}
