WITCH-NG
================

[![shields.io](https://img.shields.io/badge/version-0.0.2-blue?style=for-the-badge)]()

Fast, efficient sequence adding to existing alignments. An implementation of [WITCH](https://github.com/c5shen/WITCH).

## Quick start

First some requirements:

 - `HMMER` (in particular `hmmbuild`, `hmmsearch`, `hmmalign`) must be in installed (and discoverable by PATH). This requirement will be lifted in the future, but right now `v3.1b2` or newer of HMMER should all work.
 - A backbone ("reference") alignment and a backbone tree on this alignment must be present to add the query sequences. Again, this requirement will be lifted in the future for full-stack analyses.

Prebuilt binaries of WITCH-NG can be found in [releases](https://github.com/RuneBlaze/WITCH-NG/releases/), or see the "Building WITCH-NG from Scratch" section for building instructions. After installing the standalone binary, and preparing the files such as these:

 - `queries.fa`: the path to the query sequences (in unaligned FASTA) format
 - `backbone.afa`: the path to the backbone alignment (aligned FASTA)
 - `backbone.tre`: the path to the backbone tree (single-line Newick tree on `backbone.afa`)

the following command can be run:

```shell
witch-ng add -i queries.fa -b backbone.afa -t backbone.tre -o extended_alignment.afa
```

after which `extended_alignment.afa` contains an aligned version of all sequences. Note that we adopt the UPP convention of extended alignments, in which lower-case letters denote singleton insertion sites (i.e., anything lower-case is not homologous to anything). Future work will allow automatic masking of the singleton insertion sites.

## Building WITCH-NG from Scratch

WITCH-NG is written in Rust. Therefore, install the [Rust toolchain](https://www.rust-lang.org/tools/install) and then compile the binary should work:

```shell
cargo build --release
```

## Misc

 - This project obvious reuses stuff from [WITCH](https://github.com/c5shen/WITCH).
 - This project uses `ogcat`, which in turn translates code from [TreeSwift](https://niema.net/TreeSwift/).