use anyhow::bail;
use seq_io::fasta::OwnedRecord;
use seq_io::BaseRecord;
use std::fs::File;
use std::process::Stdio;
use std::{fs::rename, path::PathBuf, process::Command};

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
    // let stdout: &[u8] = &output.stdout;
    // let mut reader = seq_io::fasta::Reader::new(stdout);
    // let r: Result<Vec<_>, _> = reader.records().into_iter().into_iter().collect();
    // Ok(r?)
}
