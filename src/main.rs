//! # WITCH-NG
//!
//! WITCH-NG (code-name `crucible`) aims to be an efficient implementation of the WITCH algorithm
//! for aligning fragments to an existing alignment (called a "reference"
//! or "backbone" alignment).
mod adder;
mod combined;
mod compact_printer;
mod config;
mod external;
mod matching;
mod melt;
mod progress_reporter;
mod score_calc;
mod structures;

use anyhow::Ok;
use clap::{Parser, Subcommand};
use std::{path::PathBuf, time::Instant};
use tracing::{debug, info, warn};

use crate::config::ExternalContext;

#[derive(Parser, Debug, Hash, PartialEq)]
#[clap(author, version, about)]
struct Args {
    #[clap(subcommand)]
    cmd: SubCommand,
}

#[derive(Subcommand, Debug, PartialEq, Hash)]
enum SubCommand {
    /// Add query sequences to a reference alignment
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
        /// The HMM decomposition size lower bound; how many sequences must each HMM contain? Defaults to 10
        #[clap(long)]
        hmm_size_lb: Option<usize>,
        /// Specify to use a strategy for IO bound situations; make each worker use two threads, with one thread for IO
        #[clap(long)]
        io_bound: bool,
        /// Enable checkpointing in the hmmsearch stage; the checkpoint file will be a suffix of the output file
        #[clap(long)]
        checkpoint: bool,
        /// Log progress every ten seconds for the search phase
        #[clap(long)]
        progress: bool,
        /// Set level of parallelism; defaults to number of logical cores
        #[clap(long)]
        threads: Option<usize>,
    },
}

fn main() -> anyhow::Result<()> {
    let now = Instant::now();
    let args = Args::parse();
    tracing_subscriber::fmt::init();
    debug!("command: {:?}", &args.cmd);
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
            hmm_size_lb,
            io_bound,
            checkpoint,
            progress,
        } => {
            let checkpoint_path = output.with_extension("checkpoint");

            let nworkers = if let Some(t) = threads {
                t
            } else {
                if io_bound {
                    num_cpus::get() / 2
                } else {
                    num_cpus::get()
                }
            };

            let nthreads_per_worker = if io_bound { 2 } else { 1 };
            let external_context = ExternalContext {
                hmm_size_lb: hmm_size_lb.unwrap_or(10),
                show_progress: progress,
                io_bound,
                trim,
                only_queries,
                db: checkpoint.then(|| {
                    warn!("reloading a checkpoint is only guaranteed to work on the same input files and architecture, with same hmm-size-lb");
                    sled::Config::default()
                        .path(&checkpoint_path)
                        .flush_every_ms(Some(3000))
                        .use_compression(true)
                        .compression_factor(3)
                        .open().expect(
                        format!(
                            "failed to open checkpoint file at {:?}",
                            checkpoint_path.clone()
                        )
                        .as_str(),
                    )
                }),
                num_workers: nworkers,
                num_threads_per_worker: nthreads_per_worker,
            };

            rayon::ThreadPoolBuilder::new()
                .num_threads(nworkers)
                .build_global()?;
            info!("using {:?} workers", nworkers);
            if checkpoint {
                let num_entries = &external_context.db.as_ref().unwrap().len();
                if external_context.db.as_ref().unwrap().was_recovered() {
                    info!("recovered from checkpoint file at {:?}", &checkpoint_path);
                    debug!("checkpoint file contains {:?} entries", num_entries);
                }
            }
            combined::combined_analysis(
                input,
                backbone,
                output,
                ehmm_path,
                tree,
                &external_context,
            )?;
        }
    }
    info!("total elapsed time: {:?}", now.elapsed());
    Ok(())
}
