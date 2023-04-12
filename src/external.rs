use ahash::AHashMap;
use anyhow::bail;
use lazy_static::lazy_static;
use regex::Regex;
use seq_io::fasta::OwnedRecord;
use seq_io::BaseRecord;
use std::process::Stdio;
use std::{path::PathBuf, process::Command};
use tracing::debug;

use crate::config::ExternalContext;

pub fn hmmalign<'a, R>(hmm_path: &PathBuf, seqs: R) -> anyhow::Result<Vec<u8>>
where
    R: Iterator<Item = &'a OwnedRecord>,
{
    let mut child = Command::new("hmmalign")
        .arg("--informat")
        .arg("fasta")
        .arg("--outformat")
        .arg("afa")
        .arg(hmm_path)
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;
    if let Some(mut stdin) = child.stdin.take() {
        for s in seqs {
            s.write(&mut stdin)?;
        }
    } else {
        bail!("Failed to get stdin handle");
    }
    let output = child.wait_with_output()?;
    if !output.status.success() {
        bail!("hmmalign failed: {:?}", output);
    }
    Ok(output.stdout)
}

pub fn hmmbuild<'a, R>(seqs: R, name: &str, outpath: &PathBuf) -> anyhow::Result<()>
where
    R: Iterator<Item = &'a OwnedRecord>,
{
    let mut child = Command::new("hmmbuild")
        .arg("--cpu")
        .arg("0")
        .arg("--informat")
        .arg("afa")
        .arg("--ere")
        .arg("0.59")
        .arg("--symfrac")
        .arg("0.0")
        .arg("-n")
        .arg(name)
        .arg(outpath)
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;
    if let Some(mut stdin) = child.stdin.take() {
        for s in seqs {
            s.write(&mut stdin)?;
        }
    } else {
        bail!("Failed to get stdin handle");
    }
    let output = child.wait_with_output()?;
    if !output.status.success() {
        bail!("hmmbuild failed: {:?}", output);
    }
    Ok(())
}

pub fn hmmsearch<'a, R>(
    hmm_path: &PathBuf,
    seqs: R,
    seq_id: &AHashMap<String, u32>,
    config: &ExternalContext,
) -> anyhow::Result<Vec<(u32, f64)>>
where
    R: Iterator<Item = &'a OwnedRecord>,
{
    let mut child = Command::new("hmmsearch")
        .arg("--cpu")
        .arg(if config.io_bound { "1" } else { "0" })
        .arg("--noali")
        .arg("--max")
        .arg("-E")
        .arg("999999999")
        .arg(hmm_path)
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;
    let mut cnt = 0;
    if let Some(mut stdin) = child.stdin.take() {
        for s in seqs {
            s.write(&mut stdin)?;
            cnt += 1;
        }
    } else {
        bail!("Failed to get stdin handle");
    }
    debug!("{} sequences written to hmmsearch", cnt);
    let output = child.wait_with_output()?;
    if !output.status.success() {
        bail!("hmmsearch failed: {:?}", output);
    }
    lazy_static! {
        static ref RE: Regex = Regex::new(r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)").unwrap();
    }
    let raw_output = String::from_utf8(output.stdout)?;
    let mut res: Vec<(u32, f64)> = vec![];
    let mut start_reading = false;
    for l in raw_output.lines().map(|l| l.trim()) {
        if !start_reading && l.starts_with("E-value") {
            start_reading = true;
        } else if start_reading && l.is_empty() {
            break;
        } else if start_reading {
            if let Some(caps) = RE.captures(l) {
                let entire = caps.get(0).unwrap().as_str();
                if entire.contains("--") {
                    continue;
                }
                let seq_name = caps.get(9).unwrap().as_str();
                let seq_id = *seq_id.get(seq_name).unwrap();
                let bitscore = caps.get(2).unwrap().as_str().parse::<f64>()?;
                res.push((seq_id, bitscore));
            }
        }
    }
    Ok(res)
}
