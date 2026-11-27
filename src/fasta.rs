//! Streaming FASTA and FASTQ reading, gzip-aware, plus assembly statistics.

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Read, Seek};
use std::path::Path;

/// Opens a sequence file, transparently decompressing gzip by magic bytes.
pub fn open(path: &Path) -> Result<Box<dyn BufRead + Send>> {
    let mut file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut magic = [0u8; 2];
    let mut n = 0;
    while n < magic.len() {
        match file.read(&mut magic[n..]) {
            Ok(0) => break,
            Ok(got) => n += got,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e.into()),
        }
    }
    let is_gzip = n == 2 && magic == [0x1f, 0x8b];
    if file.rewind().is_ok() {
        if is_gzip {
            Ok(Box::new(BufReader::with_capacity(
                1 << 20,
                MultiGzDecoder::new(file),
            )))
        } else {
            Ok(Box::new(BufReader::with_capacity(1 << 20, file)))
        }
    } else {
        // Non-seekable inputs (e.g. FIFOs): prepend the sniffed bytes.
        let chained = Cursor::new(magic).take(n as u64).chain(file);
        if is_gzip {
            Ok(Box::new(BufReader::with_capacity(
                1 << 20,
                MultiGzDecoder::new(chained),
            )))
        } else {
            Ok(Box::new(BufReader::with_capacity(1 << 20, chained)))
        }
    }
}

/// Calls `f(name, sequence)` for each record. FASTA records may span lines;
/// FASTQ records are the standard four lines.
pub fn for_each_record(path: &Path, mut f: impl FnMut(&str, &[u8])) -> Result<()> {
    let mut reader = open(path)?;
    let mut line = Vec::with_capacity(256);
    let mut skip = Vec::new();
    let mut name = String::new();
    let mut seq: Vec<u8> = Vec::with_capacity(1 << 20);
    let mut in_record = false;

    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        trim_newline(&mut line);
        match line.first() {
            Some(b'>') => {
                if in_record {
                    f(&name, &seq);
                }
                set_header_name(&line[1..], &mut name);
                seq.clear();
                in_record = true;
            }
            Some(b'@') if !in_record || seq.is_empty() => {
                // FASTQ: header, sequence, '+', quality.
                set_header_name(&line[1..], &mut name);
                line.clear();
                reader.read_until(b'\n', &mut line)?;
                trim_newline(&mut line);
                f(&name, &line);
                skip.clear();
                reader.read_until(b'\n', &mut skip)?;
                skip.clear();
                reader.read_until(b'\n', &mut skip)?;
                in_record = false;
            }
            Some(_) if in_record => seq.extend_from_slice(&line),
            Some(_) => bail!("{}: sequence data before any header", path.display()),
            None => {}
        }
    }
    if in_record {
        f(&name, &seq);
    }
    Ok(())
}

fn trim_newline(line: &mut Vec<u8>) {
    while matches!(line.last(), Some(b'\n' | b'\r')) {
        line.pop();
    }
}

fn set_header_name(header: &[u8], out: &mut String) {
    out.clear();
    let text = String::from_utf8_lossy(header);
    if let Some(first) = text.split_whitespace().next() {
        out.push_str(first);
    }
}

/// Assembly-level facts used for genome quality control.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct SeqStats {
    pub length: u64,
    pub contigs: u32,
    pub gc: f64,
    pub n50: u64,
    /// Bases other than A, C, G, T (N runs, IUPAC codes).
    pub ambiguous: u64,
}

#[derive(Default)]
pub struct StatsBuilder {
    lengths: Vec<u64>,
    gc: u64,
    acgt: u64,
    ambiguous: u64,
}

impl StatsBuilder {
    pub fn add(&mut self, seq: &[u8]) {
        self.lengths.push(seq.len() as u64);
        for &b in seq {
            match b {
                b'G' | b'C' | b'g' | b'c' => {
                    self.gc += 1;
                    self.acgt += 1;
                }
                b'A' | b'T' | b'a' | b't' => self.acgt += 1,
                _ => self.ambiguous += 1,
            }
        }
    }

    /// Records pre-tallied base counts for a contig of length `len`, matching
    /// `add` without rescanning the sequence.
    #[inline]
    pub fn add_counts(&mut self, len: u64, gc: u64, ambiguous: u64) {
        debug_assert!(ambiguous <= len && gc <= len - ambiguous);
        self.lengths.push(len);
        self.gc += gc;
        self.acgt += len - ambiguous;
        self.ambiguous += ambiguous;
    }

    pub fn finish(mut self) -> SeqStats {
        let length: u64 = self.lengths.iter().sum();
        self.lengths.sort_unstable_by(|a, b| b.cmp(a));
        let mut running = 0;
        let mut n50 = 0;
        for &len in &self.lengths {
            running += len;
            if running * 2 >= length {
                n50 = len;
                break;
            }
        }
        SeqStats {
            length,
            contigs: self.lengths.len() as u32,
            gc: if self.acgt == 0 {
                0.0
            } else {
                self.gc as f64 / self.acgt as f64
            },
            n50,
            ambiguous: self.ambiguous,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn collect(path: &Path) -> Vec<(String, String)> {
        let mut out = Vec::new();
        for_each_record(path, |n, s| {
            out.push((n.to_string(), String::from_utf8_lossy(s).into_owned()))
        })
        .unwrap();
        out
    }

    #[test]
    fn reads_multiline_fasta_and_gzip() {
        let dir = tempfile::tempdir().unwrap();
        let text = ">c1 desc\nACGT\nAC\r\n>c2\nGGGG\n";
        let plain = dir.path().join("a.fna");
        std::fs::write(&plain, text).unwrap();
        let gz = dir.path().join("a.fna.gz");
        let mut enc =
            flate2::write::GzEncoder::new(File::create(&gz).unwrap(), flate2::Compression::fast());
        enc.write_all(text.as_bytes()).unwrap();
        enc.finish().unwrap();

        let want = vec![
            ("c1".to_string(), "ACGTAC".to_string()),
            ("c2".to_string(), "GGGG".to_string()),
        ];
        assert_eq!(collect(&plain), want);
        assert_eq!(collect(&gz), want);
    }

    #[test]
    fn reads_fastq() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("r.fq");
        std::fs::write(&path, "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\n@III\n").unwrap();
        assert_eq!(
            collect(&path),
            vec![
                ("r1".to_string(), "ACGT".to_string()),
                ("r2".to_string(), "TTTT".to_string())
            ]
        );
    }

    #[test]
    fn gzip_reads_every_member_and_rejects_truncation() {
        // bgzip and `cat a.gz b.gz` produce several gzip members; a cut-off
        // download must fail instead of reading as a shorter genome.
        let gz = |text: &str| {
            let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
            enc.write_all(text.as_bytes()).unwrap();
            enc.finish().unwrap()
        };
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("g.fna.gz");
        std::fs::write(&path, [gz(">a\nACGT\n"), gz("AC\n>b\nGGGG\n")].concat()).unwrap();
        assert_eq!(
            collect(&path),
            vec![
                ("a".to_string(), "ACGTAC".to_string()),
                ("b".to_string(), "GGGG".to_string())
            ]
        );

        let whole = gz(">c\nACGTACGTACGT\n");
        std::fs::write(&path, &whole[..whole.len() - 4]).unwrap();
        assert!(for_each_record(&path, |_, _| {}).is_err());
    }

    #[test]
    fn stats_report_n50_and_gc() {
        let mut b = StatsBuilder::default();
        b.add(b"GGGGGGGGGG"); // 10
        b.add(b"AAAAANNN"); // 8, 3 ambiguous
        b.add(b"AT"); // 2
        let s = b.finish();
        assert_eq!(s.length, 20);
        assert_eq!(s.contigs, 3);
        assert_eq!(s.n50, 10);
        assert_eq!(s.ambiguous, 3);
        assert!((s.gc - 10.0 / 17.0).abs() < 1e-12);
    }

    #[test]
    fn open_handles_empty_short_fastq_and_non_seekable_streams() {
        let dir = tempfile::tempdir().unwrap();

        // Empty and 1-byte files (< 2 magic bytes) open and rewind cleanly.
        let empty = dir.path().join("empty.fna");
        std::fs::write(&empty, b"").unwrap();
        assert!(collect(&empty).is_empty());

        let one = dir.path().join("one.fna");
        std::fs::write(&one, b">").unwrap();
        assert_eq!(collect(&one), vec![("".to_string(), "".to_string())]);

        // Multi-read FASTQ with quality lines starting with '@' and '+'.
        let fq = dir.path().join("multi.fq");
        let mut text = String::new();
        let mut want = Vec::new();
        for i in 0..32 {
            let seq = if i % 2 == 0 { "ACGTACGT" } else { "GGCCAT" };
            let qual = if i % 2 == 0 { "@+IIIIII" } else { "+@IIII" };
            text.push_str(&format!("@read_{i} extra info\n{seq}\n+\n{qual}\n"));
            want.push((format!("read_{i}"), seq.to_string()));
        }
        std::fs::write(&fq, &text).unwrap();
        assert_eq!(collect(&fq), want);

        // Non-seekable FIFO fallback for both plain and gzipped inputs.
        #[cfg(unix)]
        {
            use std::ffi::CString;
            use std::os::unix::ffi::OsStrExt;
            let fasta_text = ">c1 desc\nACGT\nAC\n>c2\nGGGG\n";
            let expected = vec![
                ("c1".to_string(), "ACGTAC".to_string()),
                ("c2".to_string(), "GGGG".to_string()),
            ];
            let mut gz_bytes = Vec::new();
            {
                let mut enc =
                    flate2::write::GzEncoder::new(&mut gz_bytes, flate2::Compression::fast());
                enc.write_all(fasta_text.as_bytes()).unwrap();
                enc.finish().unwrap();
            }
            for (name, payload) in [
                ("plain.fifo", fasta_text.as_bytes().to_vec()),
                ("gz.fifo", gz_bytes),
            ] {
                let fifo = dir.path().join(name);
                let c_fifo = CString::new(fifo.as_os_str().as_bytes()).unwrap();
                // SAFETY: c_fifo is a valid NUL-terminated path.
                assert_eq!(unsafe { libc::mkfifo(c_fifo.as_ptr(), 0o600) }, 0);
                let writer_path = fifo.clone();
                let handle = std::thread::spawn(move || {
                    let mut f = File::create(&writer_path).unwrap();
                    f.write_all(&payload).unwrap();
                });
                assert_eq!(collect(&fifo), expected, "{name}");
                handle.join().unwrap();
            }
        }
    }
}
