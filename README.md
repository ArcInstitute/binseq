# BINSEQ

[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE.md)
![actions status](https://github.com/arcinstitute/binseq/workflows/CI/badge.svg)
[![Crates.io](https://img.shields.io/crates/d/binseq?color=orange&label=crates.io)](https://crates.io/crates/binseq)
[![docs.rs](https://img.shields.io/docsrs/binseq?color=green&label=docs.rs)](https://docs.rs/binseq/latest/binseq/)

## Overview

BINSEQ is a binary file format family for efficient storage and processing of DNA sequences.
It uses two-bit encoding for nucleotides and is optimized for high-performance parallel processing.

The recommended variant is **CBQ** (`*.cbq`): a columnar, block-compressed format for variable-length records with optional quality scores and headers.
It is lossless by default (native `N` support), compresses well, and decodes fast.
For details on its structure see the [documentation](https://docs.rs/binseq/latest/binseq/cbq/).

Two earlier variants remain supported:

- **BQ** (`*.bq`): fixed-length records without quality scores. Minimal and fast for uniform reads.
- **VBQ** (`*.vbq`): variable-length records with optional quality scores and headers. Superseded by CBQ, which improves on its compression and decoding speed; new projects should use CBQ.

All variants support both single and paired sequences.

## Getting Started

This is a **library** for reading and writing BINSEQ files; for a **command-line interface** see [bqtools](https://github.com/arcinstitute/bqtools).

To get started please refer to our [documentation](https://docs.rs/binseq/latest/binseq/).
For example programs which make use of the library check out our [examples directory](https://github.com/arcinstitute/binseq/tree/main/examples).

For more information about the BINSEQ file family, please refer to our [preprint](https://www.biorxiv.org/content/10.1101/2025.04.08.647863v2).
