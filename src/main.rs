mod adder;
mod compact_printer;
mod external;
mod matching;
mod melt;
mod structures;

use std::{path::PathBuf, time::Instant};

use anyhow::Ok;
use clap::{Parser, Subcommand};
use melt::oneshot_melt;
use tracing::info;

use crate::adder::oneshot_add_queries;

#[derive(Parser, Debug, Hash, PartialEq)]
#[clap(author, version, about)]
struct Args {
    #[clap(subcommand)]
    cmd: SubCommand,
}

#[derive(Subcommand, Debug, PartialEq, Hash)]
enum SubCommand {
    /// Decompose input alignment by a tree into MSAs ready to become HMMs
    Melt {
        #[clap(short, long)]
        input: PathBuf,
        #[clap(short, long)]
        tree: PathBuf,
        #[clap(short, long)]
        outdir: PathBuf,
        #[clap(short = 's', long)]
        max_size: usize,
    },

    /// Receive payload from WITCH frontend and merges in the query sequences
    Dance {
        #[clap(short, long)]
        root: PathBuf,
    },
}

fn main() -> anyhow::Result<()> {
    let now = Instant::now();
    let args = Args::parse();
    tracing_subscriber::fmt::init();
    match args.cmd {
        SubCommand::Melt {
            input,
            tree,
            outdir,
            max_size,
        } => {
            oneshot_melt(&input, &tree, max_size, &outdir)?;
        }
        SubCommand::Dance { root } => {
            oneshot_add_queries(&root)?;
        }
    }
    info!("total elapsed time: {:?}", now.elapsed());
    Ok(())
}
