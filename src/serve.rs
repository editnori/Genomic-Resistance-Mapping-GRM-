//! Serves the built dashboard and its report data from one process.
//!
//!   /data/<path>   files under the report directory (`grm report` output)
//!   /<path>        files under the built web app; unknown paths get
//!                  index.html so client-side routes load

use anyhow::{Context, Result};
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::{Component, Path, PathBuf};
use std::sync::{Mutex, mpsc};
use tiny_http::{Header, Method, Request, Response, ResponseBox, StatusCode};

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

fn resolve_file(root: &Path, url_path: &str) -> Option<PathBuf> {
    let path = resolve(root, url_path)?.canonicalize().ok()?;
    if path.is_file() && path.starts_with(root) {
        Some(path)
    } else {
        None
    }
}

fn header(name: &'static [u8], value: &[u8]) -> Header {
    Header::from_bytes(name, value).unwrap()
}

fn text_response(status: u16, text: &str) -> ResponseBox {
    Response::from_string(text)
        .with_status_code(StatusCode(status))
        .with_header(header(b"Cache-Control", b"no-store"))
        .boxed()
}

fn internal_error(path: &Path, error: &std::io::Error) -> ResponseBox {
    eprintln!("  {}: {error}", path.display());
    text_response(500, "internal server error")
}

enum RangeSelection {
    Full,
    Partial(u64, u64),
    Invalid,
}

fn requested_range(headers: &[Header], file_length: u64) -> RangeSelection {
    let mut range_headers = headers.iter().filter(|header| header.field.equiv("Range"));
    let Some(range) = range_headers.next() else {
        return RangeSelection::Full;
    };
    if range_headers.next().is_some() {
        return RangeSelection::Invalid;
    }
    let Some(spec) = range.value.as_str().strip_prefix("bytes=") else {
        return RangeSelection::Invalid;
    };
    if spec.contains(',') {
        return RangeSelection::Invalid;
    }
    let Some((first, last)) = spec.split_once('-') else {
        return RangeSelection::Invalid;
    };
    if first.is_empty() {
        let Ok(suffix) = last.parse::<u64>() else {
            return RangeSelection::Invalid;
        };
        if suffix == 0 || file_length == 0 {
            return RangeSelection::Invalid;
        }
        let length = suffix.min(file_length);
        return RangeSelection::Partial(file_length - length, length);
    }
    let Ok(start) = first.parse::<u64>() else {
        return RangeSelection::Invalid;
    };
    if start >= file_length {
        return RangeSelection::Invalid;
    }
    let end = if last.is_empty() {
        file_length - 1
    } else {
        let Ok(end) = last.parse::<u64>() else {
            return RangeSelection::Invalid;
        };
        end.min(file_length - 1)
    };
    if end < start {
        return RangeSelection::Invalid;
    }
    RangeSelection::Partial(start, end - start + 1)
}

fn file_response(path: &Path, web_dist: &Path, url: &str, headers: &[Header]) -> ResponseBox {
    let mut file = match File::open(path) {
        Ok(file) => file,
        Err(error) => return internal_error(path, &error),
    };
    let file_length = match file.metadata() {
        Ok(metadata) => metadata.len(),
        Err(error) => return internal_error(path, &error),
    };
    let (status, range) = match requested_range(headers, file_length) {
        RangeSelection::Full => (200, None),
        RangeSelection::Partial(start, length) => (206, Some((start, length))),
        RangeSelection::Invalid => {
            let content_range = format!("bytes */{file_length}");
            return text_response(416, "range not satisfiable")
                .with_header(header(b"Accept-Ranges", b"bytes"))
                .with_header(header(b"Content-Range", content_range.as_bytes()));
        }
    };
    let mut response = match range {
        Some((start, length)) => {
            if let Err(error) = file.seek(SeekFrom::Start(start)) {
                return internal_error(path, &error);
            }
            Response::new(
                StatusCode(status),
                Vec::new(),
                file.take(length),
                Some(length as usize),
                None,
            )
            .boxed()
        }
        None => Response::from_file(file).boxed(),
    };
    let cache_control = if url.starts_with("/assets/")
        && path
            .strip_prefix(web_dist)
            .is_ok_and(|relative| relative.starts_with("assets"))
    {
        "public, max-age=31536000, immutable"
    } else {
        "no-cache"
    };
    response = response
        .with_header(header(b"Content-Type", content_type(path).as_bytes()))
        .with_header(header(b"Cache-Control", cache_control.as_bytes()))
        .with_header(header(b"Accept-Ranges", b"bytes"));
    if let Some((start, length)) = range {
        let content_range = format!("bytes {}-{}/{file_length}", start, start + length - 1);
        response = response.with_header(header(b"Content-Range", content_range.as_bytes()));
    }
    response
}

fn handle_request(request: Request, web_dist: &Path, report: &Path, index: &Path) {
    let url = request.url().split('?').next().unwrap_or("/").to_string();
    let response = if !matches!(request.method(), Method::Get | Method::Head) {
        text_response(405, "method not allowed").with_header(header(b"Allow", b"GET, HEAD"))
    } else {
        let target = match url.strip_prefix("/data/") {
            Some(rest) => resolve_file(report, rest),
            None if url == "/" => Some(index.to_path_buf()),
            None => resolve_file(web_dist, &url).or_else(|| Some(index.to_path_buf())),
        };
        match target {
            Some(path) => file_response(&path, web_dist, &url, request.headers()),
            None => text_response(404, "not found"),
        }
    };
    if let Err(error) = request.respond(response) {
        eprintln!("  {url}: {error}");
    }
}

pub fn serve(web_dist: &Path, report: &Path, port: u16) -> Result<()> {
    let web_dist = web_dist
        .canonicalize()
        .with_context(|| format!("app root {}", web_dist.display()))?;
    let report = report
        .canonicalize()
        .with_context(|| format!("report root {}", report.display()))?;
    let index = resolve_file(&web_dist, "index.html").context("app index.html")?;
    let server = tiny_http::Server::http(("127.0.0.1", port))
        .map_err(|e| anyhow::anyhow!("bind port {port}: {e}"))?;
    eprintln!(
        "serving http://127.0.0.1:{port}  (app {}, data {})",
        web_dist.display(),
        report.display()
    );
    let (sender, receiver) = mpsc::sync_channel(64);
    let receiver = Mutex::new(receiver);
    std::thread::scope(|scope| {
        for _ in 0..4 {
            let receiver = &receiver;
            let web_dist = &web_dist;
            let report = &report;
            let index = &index;
            scope.spawn(move || {
                loop {
                    let request = match receiver.lock() {
                        Ok(receiver) => receiver.recv(),
                        Err(_) => break,
                    };
                    match request {
                        Ok(request) => handle_request(request, web_dist, report, index),
                        Err(_) => break,
                    }
                }
            });
        }
        for request in server.incoming_requests() {
            if sender.send(request).is_err() {
                break;
            }
        }
        drop(sender);
    });
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;
    use std::fs;
    use std::io::{self, Read, Write};
    use std::net::{SocketAddr, TcpListener, TcpStream};
    use std::thread;
    use std::time::Duration;
    use tempfile::TempDir;

    struct HttpResponse {
        status: u16,
        headers: BTreeMap<String, String>,
        body: Vec<u8>,
    }

    fn fixture() -> TempDir {
        let temp = tempfile::tempdir().unwrap();
        fs::create_dir_all(temp.path().join("app/assets")).unwrap();
        fs::create_dir_all(temp.path().join("report")).unwrap();
        fs::write(temp.path().join("app/index.html"), "<html>app</html>").unwrap();
        fs::write(
            temp.path().join("app/assets/index-BfBf1GGn.js"),
            "console.log(1);",
        )
        .unwrap();
        fs::write(temp.path().join("report/index.json"), "{\"ok\":true}").unwrap();
        temp
    }

    fn start_server(web_dist: &Path, report: &Path) -> u16 {
        let listener = TcpListener::bind(("127.0.0.1", 0)).unwrap();
        let port = listener.local_addr().unwrap().port();
        drop(listener);
        let web_dist = web_dist.to_path_buf();
        let report = report.to_path_buf();
        thread::spawn(move || serve(&web_dist, &report, port).unwrap());
        let address: SocketAddr = ([127, 0, 0, 1], port).into();
        for _ in 0..100 {
            if let Ok(stream) = TcpStream::connect_timeout(&address, Duration::from_millis(20)) {
                drop(stream);
                return port;
            }
            thread::sleep(Duration::from_millis(20));
        }
        panic!("server did not start on port {port}");
    }

    fn request(port: u16, method: &str, path: &str, headers: &[(&str, &str)]) -> HttpResponse {
        request_with_timeout(port, method, path, headers, Duration::from_secs(3)).unwrap()
    }

    fn request_with_timeout(
        port: u16,
        method: &str,
        path: &str,
        headers: &[(&str, &str)],
        timeout: Duration,
    ) -> io::Result<HttpResponse> {
        let mut stream = TcpStream::connect(("127.0.0.1", port))?;
        stream.set_read_timeout(Some(timeout))?;
        stream.set_write_timeout(Some(timeout))?;
        let mut request =
            format!("{method} {path} HTTP/1.1\r\nHost: localhost\r\nConnection: close\r\n");
        for (name, value) in headers {
            request.push_str(&format!("{name}: {value}\r\n"));
        }
        request.push_str("\r\n");
        stream.write_all(request.as_bytes())?;
        let mut bytes = Vec::new();
        stream.read_to_end(&mut bytes)?;
        let header_end = bytes
            .windows(4)
            .position(|window| window == b"\r\n\r\n")
            .map(|position| position + 4)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "incomplete response"))?;
        let head = std::str::from_utf8(&bytes[..header_end])
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
        let mut lines = head.split("\r\n");
        let status = lines
            .next()
            .and_then(|line| line.split_whitespace().nth(1))
            .and_then(|status| status.parse().ok())
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid status"))?;
        let headers = lines
            .filter_map(|line| line.split_once(':'))
            .map(|(name, value)| (name.to_ascii_lowercase(), value.trim().to_string()))
            .collect();
        Ok(HttpResponse {
            status,
            headers,
            body: bytes[header_end..].to_vec(),
        })
    }

    fn get(port: u16, path: &str) -> HttpResponse {
        request(port, "GET", path, &[])
    }

    #[test]
    fn refuses_paths_that_leave_the_root() {
        let root = Path::new("/srv/app");
        assert_eq!(
            resolve(root, "/assets/a.js"),
            Some(PathBuf::from("/srv/app/assets/a.js"))
        );
        assert_eq!(resolve(root, "/../etc/passwd"), None);
        assert_eq!(resolve(root, "/assets/../../x"), None);
        for url in [
            "/%2e%2e/etc/passwd",
            "//etc/passwd",
            r"C:\etc\passwd",
            r"/..\etc\passwd",
            r"/assets\..\..\etc\passwd",
        ] {
            if let Some(path) = resolve(root, url) {
                assert!(path.starts_with(root), "{url} resolved to {path:?}");
            }
        }
    }

    #[test]
    fn content_types_cover_the_built_app_and_reports() {
        for (path, expected) in [
            ("index.html", "text/html; charset=utf-8"),
            ("index.js", "text/javascript; charset=utf-8"),
            ("index.mjs", "text/javascript; charset=utf-8"),
            ("index.css", "text/css; charset=utf-8"),
            ("index.json", "application/json"),
            ("favicon.svg", "image/svg+xml"),
            ("font.woff2", "font/woff2"),
            ("font.woff", "font/woff"),
            ("favicon.ico", "image/x-icon"),
            ("image.png", "image/png"),
        ] {
            assert_eq!(content_type(Path::new(path)), expected);
        }
    }

    #[cfg(unix)]
    #[test]
    fn symlinks_cannot_escape_either_root() {
        let temp = fixture();
        let outside = temp.path().join("outside.json");
        fs::write(&outside, "outside secret").unwrap();
        std::os::unix::fs::symlink(&outside, temp.path().join("app/leak.json")).unwrap();
        std::os::unix::fs::symlink(&outside, temp.path().join("report/leak.json")).unwrap();
        let port = start_server(&temp.path().join("app"), &temp.path().join("report"));
        let app = get(port, "/leak.json");
        let report = get(port, "/data/leak.json");
        assert_eq!(app.body, b"<html>app</html>");
        assert_eq!(report.status, 404);
        assert_ne!(report.body, b"outside secret");
    }

    #[test]
    fn head_matches_get_headers_without_a_body() {
        let temp = fixture();
        let port = start_server(&temp.path().join("app"), &temp.path().join("report"));
        let get = get(port, "/");
        let head = request(port, "HEAD", "/", &[]);
        assert_eq!(head.status, 200);
        assert!(head.body.is_empty());
        assert_eq!(
            head.headers.get("content-type"),
            get.headers.get("content-type")
        );
        assert_eq!(
            head.headers.get("content-length"),
            get.headers.get("content-length")
        );
    }

    #[test]
    fn serves_bounded_single_byte_ranges() {
        let temp = fixture();
        fs::write(temp.path().join("report/data.json"), "0123456789").unwrap();
        let port = start_server(&temp.path().join("app"), &temp.path().join("report"));
        let middle = request(port, "GET", "/data/data.json", &[("Range", "bytes=2-5")]);
        let open = request(port, "GET", "/data/data.json", &[("Range", "bytes=7-")]);
        let suffix = request(port, "GET", "/data/data.json", &[("Range", "bytes=-3")]);
        let invalid = request(port, "GET", "/data/data.json", &[("Range", "bytes=20-")]);
        let head = request(port, "HEAD", "/data/data.json", &[("Range", "bytes=2-5")]);
        assert_eq!(middle.status, 206);
        assert_eq!(middle.body, b"2345");
        assert_eq!(
            middle.headers.get("content-range").map(String::as_str),
            Some("bytes 2-5/10")
        );
        assert_eq!(
            middle.headers.get("content-length").map(String::as_str),
            Some("4")
        );
        assert_eq!(open.body, b"789");
        assert_eq!(suffix.body, b"789");
        assert_eq!(invalid.status, 416);
        assert_eq!(
            invalid.headers.get("content-range").map(String::as_str),
            Some("bytes */10")
        );
        assert_eq!(head.status, 206);
        assert!(head.body.is_empty());
        assert_eq!(
            head.headers.get("content-range"),
            middle.headers.get("content-range")
        );
    }

    #[test]
    fn applies_cache_policy_by_asset_lifetime() {
        let temp = fixture();
        let port = start_server(&temp.path().join("app"), &temp.path().join("report"));
        let asset = get(port, "/assets/index-BfBf1GGn.js");
        let index = get(port, "/");
        let report = get(port, "/data/index.json");
        let missing = get(port, "/data/missing.json");
        assert_eq!(
            asset.headers.get("cache-control").map(String::as_str),
            Some("public, max-age=31536000, immutable")
        );
        assert_eq!(
            index.headers.get("cache-control").map(String::as_str),
            Some("no-cache")
        );
        assert_eq!(
            report.headers.get("cache-control").map(String::as_str),
            Some("no-cache")
        );
        assert_eq!(missing.status, 404);
        assert_eq!(
            missing.headers.get("cache-control").map(String::as_str),
            Some("no-store")
        );
        assert_eq!(
            missing.headers.get("content-type").map(String::as_str),
            Some("text/plain; charset=UTF-8")
        );
    }

    #[test]
    fn rejects_methods_other_than_get_and_head() {
        let temp = fixture();
        let port = start_server(&temp.path().join("app"), &temp.path().join("report"));
        let response = request(port, "POST", "/data/index.json", &[]);
        assert_eq!(response.status, 405);
        assert_eq!(
            response.headers.get("allow").map(String::as_str),
            Some("GET, HEAD")
        );
        assert_eq!(get(port, "/data/index.json").status, 200);
    }

    #[test]
    fn a_slow_download_does_not_block_other_clients() {
        let temp = fixture();
        fs::File::create(temp.path().join("app/large.bin"))
            .unwrap()
            .set_len(32 * 1024 * 1024)
            .unwrap();
        let port = start_server(&temp.path().join("app"), &temp.path().join("report"));
        let mut slow = TcpStream::connect(("127.0.0.1", port)).unwrap();
        slow.set_read_timeout(Some(Duration::from_secs(3))).unwrap();
        slow.write_all(b"GET /large.bin HTTP/1.1\r\nHost: localhost\r\nConnection: close\r\n\r\n")
            .unwrap();
        let mut headers = Vec::new();
        while !headers.ends_with(b"\r\n\r\n") {
            let mut byte = [0];
            slow.read_exact(&mut byte).unwrap();
            headers.push(byte[0]);
            assert!(headers.len() < 16 * 1024);
        }
        assert!(
            String::from_utf8(headers)
                .unwrap()
                .starts_with("HTTP/1.1 200")
        );
        thread::sleep(Duration::from_millis(100));
        let second =
            request_with_timeout(port, "GET", "/data/index.json", &[], Duration::from_secs(2));
        drop(slow);
        let second = second.expect("a serial server blocked the second client");
        assert_eq!(second.status, 200);
        assert_eq!(second.body, b"{\"ok\":true}");
    }

    #[test]
    fn missing_asset_fallback_is_not_immutable() {
        let temp = fixture();
        let port = start_server(&temp.path().join("app"), &temp.path().join("report"));
        let response = get(port, "/assets/missing.js");
        assert_eq!(response.body, b"<html>app</html>");
        assert_eq!(
            response.headers.get("cache-control").map(String::as_str),
            Some("no-cache")
        );
    }

    #[test]
    fn file_open_errors_become_500_responses() {
        let temp = fixture();
        let missing = temp.path().join("report/missing.json");
        let response = file_response(
            &missing,
            &temp.path().join("app"),
            "/data/missing.json",
            &[],
        );
        assert_eq!(response.status_code(), StatusCode(500));
        assert_eq!(
            response
                .headers()
                .iter()
                .find(|header| header.field.equiv("Cache-Control"))
                .map(|header| header.value.as_str()),
            Some("no-store")
        );
    }
}
