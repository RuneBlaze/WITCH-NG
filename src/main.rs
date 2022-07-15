mod melt;
mod structures;

use std::path::PathBuf;

use clap::{Parser, Subcommand};

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
    Dance {},
}

fn main() {
    let args = Args::parse();
    tracing_subscriber::fmt::init();
    println!("Hello, world!");
}
