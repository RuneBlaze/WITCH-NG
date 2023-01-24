WITCH-NG
================

[![shields.io](https://img.shields.io/badge/recommended_version-0.0.2-blue?style=for-the-badge)](https://github.com/RuneBlaze/WITCH-NG/releases/tag/v0.0.2) [![shields.io](https://img.shields.io/badge/research_paper-biorxiv_preprint-blue?style=for-the-badge)](https://www.biorxiv.org/content/10.1101/2022.08.08.503232v1)

Fast, efficient sequence adding to existing alignments, a somewhat optimized implementation of [WITCH](https://github.com/c5shen/WITCH).

## Quick start

 1. You can find prebuilt binaries for this software, WITCH-NG, in the [releases](https://github.com/RuneBlaze/WITCH-NG/releases/) section. If you prefer to build the software from scratch, instructions on how to do so can be found in [this section](#building-witch-ng-from-scratch). Please note that WITCH-NG requires HMMER, so WITCH-NG is not currently supported on Windows machines. However, it may be possible to use WITCH-NG on Windows Subsystem for Linux (WSL) by compiling it. If you're using Apple Silicon, you should prefer the `aarch64-apple-darwin` prebuilt binary. The following instructions assume that `witch-ng` is under the current working directory
 2. `HMMER` (in particular `hmmbuild`, `hmmsearch`, `hmmalign`) must also be installed and discoverable by PATH. The easiest way if you don't have them can be through [package managers](http://hmmer.org/documentation.html). This requirement will be lifted in the future, but right now `v3.1b2` or newer versions of HMMER should all work.
 3. A backbone ("reference") alignment and a backbone tree on this alignment must already exist to add the query sequences. This requirement will be lifted in the future for align-from-scratch analyses. For example, the files can be named like this:

    - `queries.fa`: the path to the query sequences (in unaligned FASTA format)
    - `backbone.afa`: the path to the backbone alignment (aligned FASTA)
    - `backbone.tre`: the path to the backbone tree (single-line Newick tree inferred on `backbone.afa`)

We also provide an example version of the above files for the sake of trying WITCH-NG out. To download the example files:

```bash
# for directly using the example files
git clone https://gist.github.com/RuneBlaze/272e086436190557b715dd980fd39903 witch-ng-examples && cd witch-ng-examples
```

After these input are gathered, the following command can be run:

```bash
# using the example files, this command took less than 1 minute on my laptop
witch-ng add -i queries.fa -b backbone.afa -t backbone.tre -o extended_alignment.afa
```

The output `extended_alignment.afa` then contains an aligned version of all sequences. Note that we adopt the UPP convention of extended alignments, in which lower-case letters denote singleton insertion sites (i.e., anything lower-case is not homologous to anything). Future work will allow automatic masking of the singleton insertion sites.

## Building WITCH-NG from Scratch

WITCH-NG is written in Rust. Install the [Rust toolchain](https://www.rust-lang.org/tools/install) and then compile the binary should build the standalone binary:

```shell
cargo build --release
```

## Misc

 - This project obvious reuses stuff from [WITCH](https://github.com/c5shen/WITCH).
 - This project uses [ogcat](https://github.com/RuneBlaze/ogcat), which in turn translates code from [TreeSwift](https://niema.net/TreeSwift/).

## License

WITCH-NG is available under the GPLv3 LICENSE.