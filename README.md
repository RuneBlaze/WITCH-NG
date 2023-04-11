WITCH-NG
================

Fast, efficient sequence adding to existing alignments, a somewhat optimized implementation of [WITCH](https://github.com/c5shen/WITCH), and a more accurate version of HMMER's `hmmbuild` paired with `hmmalign` when the sequences are evolutionary diverse.

## Use cases

 - Instead of `hmmbuild+hmmalign`, the `witch-ng add` subcommand can add sequences to an existing alignment. When the sequences are quite diverse, the output alignment is likely more accurate.
 - For *de novo* sequence alignment with diverse input lengths (e.g., some sequences ~500NT in length and others ~1500NT), `witch-ng msa` can scalably align sequences. The output alignment should match or surpass the accuracy of `mafft` or [MAGUS](https://github.com/vlasmirnov/MAGUS), particularly with sequences that are evolutionarily distant from each other.

## Quick start

### Dependencies

 - [HMMER](http://hmmer.org/) (in particular `hmmbuild`, `hmmsearch`, `hmmalign`) must be installed. The easiest way might be through [package managers](http://hmmer.org/documentation.html). `v3.1b2` or newer versions of HMMER should all work.

If you can run `hmmsearch -h`, `hmmalign -h`, and `hmmbuild -h` from the command line, you should be good to go.

### Input files

 - WITCH-NG is a replacement for `hmmbuild` paired with `hmmalign`, but can also be used to align sequences from scratch. Same as `hmmbuild`+`hmmalign`, one must have a reference alignment ("backbone") and a set of query unaligned sequences ("queries") ready. The goal is to "add" the query sequences to the backbone alignment.
 - WITCH-NG, in addition to requiring the backbone MSA and the query sequences, requires a tree to be inferred on the backbone alignment. This "backbone tree" can be inferred either through FastTree or RAxML. The tree must be in the Newick format, and must be single-line

The rest of this README assumes the following files have been prepared and ready:

  - `queries.fa`: the path to the query sequences (in unaligned FASTA format)
  - `backbone.afa`: the path to the backbone alignment (aligned FASTA)
  - `backbone.tre`: the path to the backbone tree (single-line Newick tree inferred on `backbone.afa`)

We also provide an example version of the above files for the sake of trying WITCH-NG out. To download the example files:

```bash
# for example files
git clone https://gist.github.com/RuneBlaze/7d480dcedc032b4a476b06959410916b witch-ng-examples && cd witch-ng-examples

# one can also download them manually at https://gist.github.com/RuneBlaze/7d480dcedc032b4a476b06959410916b
```

### Downloading WITCH-NG

Prebuilt executables for this software, WITCH-NG, can be found in the [releases](https://github.com/RuneBlaze/WITCH-NG/releases/) section. If you prefer to compile the software from scratch, see [this section](#building-witch-ng-from-scratch). A subset of the provided executables can be found here:

| Binary Download                                                            | Operating System |
|----------------------------------------------------------------------------------|---------------------------------|
| [`witch-ng-x86_64-unknown-linux-musl.tar.gz`](https://github.com/RuneBlaze/WITCH-NG/releases/download/v0.0.2/witch-ng-x86_64-unknown-linux-musl.tar.gz)         | Linux (x86_64)                   |
| [`witch-ng-aarch64-apple-darwin.tar.gz`](https://github.com/RuneBlaze/WITCH-NG/releases/download/v0.0.2/witch-ng-aarch64-apple-darwin.tar.gz)              | macOS (Apple Silicon)            |
| [`witch-ng-x86_64-apple-darwin.tar.gz`](https://github.com/RuneBlaze/WITCH-NG/releases/download/v0.0.2/witch-ng-x86_64-apple-darwin.tar.gz)                | macOS (x86_64)                   |

For example, if you are on a Linux x86_64 machine, you can install WITCH-NG with the following command:

```bash
curl -L https://github.com/RuneBlaze/WITCH-NG/releases/download/v0.0.2/witch-ng-x86_64-unknown-linux-musl.tar.gz | tar -xz
./witch-ng --help
```


After checking that the executable `witch-ng` is installed under the current directory, the following command can be run

```bash
./witch-ng add -i queries.fa -b backbone.afa -t backbone.tre -o extended_alignment.afa
# using the example files, this command took less than 4 minute on my laptop
```

The on-screen logs from WITCH-NG will be something like this, with the output file `extended_alignment.afa`:

```bash
2023-01-25T02:59:28.499149Z  INFO witch_ng::melt: decomposed input tree num_subsets=275
2023-01-25T02:59:31.284642Z  INFO witch_ng::score_calc: read 13950 query sequences
2023-01-25T03:02:36.910727Z  INFO witch_ng::combined: all-against-all hmmsearch (with adjusted bitscore calculation) took 185.626686208s
2023-01-25T03:02:43.067799Z  INFO witch_ng::adder: output homologies formatted, output alignment will have 874 columns
2023-01-25T03:02:43.125048Z  INFO witch_ng: total elapsed time: 194.624788166s
```

The output `extended_alignment.afa` contains an aligned version of all sequences. We adopt the UPP convention of extended alignments, in which lower-case letters denote singleton insertion sites (i.e., anything lower-case is not homologous to anything). Future work will allow automatic masking of the singleton insertion sites.

## Options for `witch-ng add`

### `--threads <THREADS>`

Set level of parallelism; WITCH-NG by default uses all logical cores.

## Building WITCH-NG from Scratch

WITCH-NG is written in Rust. Install the [Rust toolchain](https://www.rust-lang.org/tools/install) and then compile the binary should build the standalone binary:

```shell
cargo build --release
```

## Misc

 - This project obvious reuses stuff from [WITCH](https://github.com/c5shen/WITCH).
 - This project uses [ogcat](https://github.com/RuneBlaze/ogcat), which in turn translates code from [TreeSwift](https://niema.net/TreeSwift/).

## Citation

If you use our software, please cite the following paper:

```bibtex
@article{witch-ng,
  title={{WITCH}-{NG}: efficient and accurate alignment of datasets with sequence length heterogeneity},
  author={Liu, Baqiao and Warnow, Tandy},
  journal={Bioinformatics Advances},
  volume={3},
  number={1},
  pages={vbad024},
  year={2023},
  publisher={Oxford University Press}
}
```

## License

WITCH-NG is available under the GPLv3 LICENSE.