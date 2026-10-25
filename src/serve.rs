//! Serves the built dashboard and its report data from one process.
//!
//!   /data/<path>   files under the report directory (`grm report` output)
//!   /<path>        files under the built web app; unknown paths get
//!                  index.html so client-side routes load

use anyhow::{Context, Result};
use std::path::{Component, Path, PathBuf};

fn content_type(path: &Path) -> &'static str {
    match path.extension().and_then(|e| e.to_str()) {
        Some("html") => "text/html; charset=utf-8",
        Some("js") | Some("mjs") => "text/javascript; charset=utf-8",
        Some("css") => "text/css; charset=utf-8",
        Some("json") => "application/json",
        Some("svg") => "image/svg+xml",
        Some("png") => "image/png",
        Some("woff2") => "font/woff2",
        Some("woff") => "font/woff",
        Some("ico") => "image/x-icon",
        _ => "application/octet-stream",
    }
}

/// Joins a URL path under `root`, refusing anything that would leave it.
fn resolve(root: &Path, url_path: &str) -> Option<PathBuf> {
    let rel = Path::new(url_path.trim_start_matches('/'));
    if rel.components().any(|c| !matches!(c, Component::Normal(_))) {
        return None;
    }
    Some(root.join(rel))
}

pub fn serve(web_dist: &Path, report: &Path, port: u16) -> Result<()> {
    let server = tiny_http::Server::http(("127.0.0.1", port))
        .map_err(|e| anyhow::anyhow!("bind port {port}: {e}"))?;
    eprintln!(
        "serving http://127.0.0.1:{port}  (app {}, data {})",
        web_dist.display(),
        report.display()
    );
    let index = web_dist.join("index.html");
    for request in server.incoming_requests() {
        let url = request.url().split('?').next().unwrap_or("/").to_string();
        let target = match url.strip_prefix("/data/") {
            Some(rest) => resolve(report, rest),
            None if url == "/" => Some(index.clone()),
            None => resolve(web_dist, &url)
                .filter(|p| p.is_file())
                .or_else(|| Some(index.clone())),
        };
        let response = match target.filter(|p| p.is_file()) {
            Some(path) => {
                let file = std::fs::File::open(&path)
                    .with_context(|| format!("open {}", path.display()))?;
                let header = tiny_http::Header::from_bytes(
                    &b"Content-Type"[..],
                    content_type(&path).as_bytes(),
                )
                .unwrap();
                request.respond(tiny_http::Response::from_file(file).with_header(header))
            }
            None => {
                request.respond(tiny_http::Response::from_string("not found").with_status_code(404))
            }
        };
        if let Err(e) = response {
            eprintln!("  {url}: {e}");
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn refuses_paths_that_leave_the_root() {
        let root = Path::new("/srv/app");
        assert_eq!(
            resolve(root, "/assets/a.js"),
            Some(PathBuf::from("/srv/app/assets/a.js"))
        );
        assert_eq!(resolve(root, "/../etc/passwd"), None);
        assert_eq!(resolve(root, "/assets/../../x"), None);
    }
}
