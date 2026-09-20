# MultiOmicsCancer

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**MultiOmicsCancer** is a modular, cancer-agnostic R package for integrative
multi-omics analysis in cancer research. It provides a unified framework for
differential methylation analysis, gene expression analysis, multi-omics
integration, survival analysis, and publication-quality visualization.

The package is designed to work with **any cancer type** — from TCGA and
cBioPortal to user-provided datasets — without requiring code modifications.

## Features

- **Differential methylation analysis** — Identify hyper/hypomethylated probes
  between tumor and normal samples using `limma`.
- **Differential gene expression analysis** — Support for both `limma` and
  `DESeq2` backends.
- **Methylation-expression integration** — Discover epigenetically regulated
  genes via correlation analysis.
- **Survival analysis** — Kaplan-Meier curves, log-rank tests, and univariate
  Cox proportional hazards models.
- **Publication-quality visualization** — Volcano plots, correlation
  histograms, and survival curves.

## Installation

You can install the development version of `MultiOmicsCancer` from GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install MultiOmicsCancer
devtools::install_github("maryam6770/MultiOmicsCancer")
