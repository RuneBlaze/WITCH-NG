//! # Crucible
//!
//! `crucible` aims to be an efficient implementation of the WITCH algorithm
//! for aligning fragments to an existing alignment (called a "reference"
//! or "backbone" alignment).
mod adder;
mod combined;
mod compact_printer;
mod external;
mod matching;
mod melt;
mod score_calc;
mod structures;

use std::{path::PathBuf, time::Instant};
use anyhow::Ok;
use clap::{Parser, Subcommand};
use tracing::info;

#[derive(Parser, Debug, Hash, PartialEq)]
#[clap(author, version, about)]
struct Args {
    #[clap(subcommand)]
    cmd: SubCommand,
}

#[derive(Subcommand, Debug, PartialEq, Hash)]
enum SubCommand {
    Add {
        /// Path to query sequences (fragments) in FASTA format
        #[clap(short, long)]
        input: PathBuf,
        /// Either a path to a full-length MSA (in FASTA, also requiring a "tree") or a directory of eHMMs
        #[clap(short, long)]
        backbone: PathBuf,
        /// Output path for the intermediate eHMMs (only used if backbone is MSA, not eHMMs, in which case defaults to backbone MSA path with "ehmm" extension)
        #[clap(short, long)]
        ehmm_path: Option<PathBuf>,
        /// Path to backbone tree for the full-length MSA (when specified in "backbone")
        #[clap(short, long)]
        tree: Option<PathBuf>,
        /// Output path of the merged MSA
        #[clap(short, long)]
        output: PathBuf,
        /// Trim singleton columns in the output. Not implemented.
        #[clap(long)]
        trim: bool,
        /// Forgo outputting the backbone; must go with "--trim"
        #[clap(long)]
        only_queries: bool,
        /// Set level of parallelism; defaults to number of logical cores
        #[clap(long)]
        threads: Option<usize>,
    },
}

fn main() -> anyhow::Result<()> {
    let now = Instant::now();
    let args = Args::parse();
    tracing_subscriber::fmt::init();
    match args.cmd {
        SubCommand::Add {
            input,
            backbone,
            ehmm_path,
            tree,
            output,
            trim,
            only_queries,
            threads,
        } => {
            if let Some(t) = threads {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(t)
                    .build_global()?;
            }
            combined::combined_analysis(
                input,
                backbone,
                output,
                ehmm_path,
                tree,
                trim,
                only_queries,
            )?;
        }
    }
    info!("total elapsed time: {:?}", now.elapsed());
    Ok(())
}
