---
title: 'MultiOmicsCancer: An R Package for Integrative Multi-Omics Analysis in Pan-Cancer Research'
tags:
  - R
  - bioinformatics
  - multi-omics
  - cancer
  - epigenetics
  - DNA methylation
  - gene expression
  - survival analysis
authors:
  - name: Maryam Emami
    orcid: 0009-0008-7710-3274
    affiliation: 1
affiliations:
  - name: Royan Institute (ACECR) & Islamic Azad University, Marvdasht, Iran
    index: 1
date: 21 September 2026
bibliography: paper.bib
---

# Summary

`MultiOmicsCancer` is an open-source R package that provides a modular,
cancer-agnostic framework for integrative multi-omics analysis in cancer
research. The package unifies several analytical steps that are typically
performed with separate tools: differential methylation analysis, differential
gene expression analysis, methylation-expression integration, survival analysis,
and publication-quality visualization. All functions operate on a single
`MultiAssayExperiment` object, ensuring consistent sample handling and
interoperability across omics layers.

Unlike existing tools that are restricted to a specific cancer type or data
source, `MultiOmicsCancer` is designed to work with **any cancer type** and
supports data from cBioPortal, TCGA, and user-provided datasets. The package
implements an automatic profile discovery mechanism for cBioPortal, which
identifies the appropriate molecular profile IDs for each study without
requiring hard-coded mappings. This design makes the package immediately
applicable to new datasets without code modifications.

# Statement of Need

Cancer is a heterogeneous disease driven by the interplay of genetic and
epigenetic alterations. Understanding these alterations requires the integration
of multiple omics layers — DNA methylation, gene expression, and clinical
outcomes — into a unified analytical framework. While several R packages
address individual omics layers (e.g., `limma`, `DESeq2`, `minfi`), few
provide an end-to-end workflow that is both cancer-agnostic and accessible to
researchers without extensive bioinformatics training.

Existing multi-omics packages such as `ELMER` and `MOFA2` offer powerful
statistical frameworks but require substantial setup and are primarily designed
for users with advanced computational expertise. In contrast, `MultiOmicsCancer`
prioritizes usability: a single function call loads data, another performs
differential analysis, and another generates publication-ready figures. This
lowers the barrier to entry for experimental biologists and clinical
researchers who wish to perform integrative analyses on their own data.

The package is particularly suited for:

- **Cancer researchers** who wish to identify epigenetically regulated genes
  from paired methylation and expression data.
- **Clinical investigators** who need to perform survival analysis on
  candidate biomarkers without writing custom code.
- **Bioinformatics educators** who seek a teaching tool that demonstrates
  multi-omics integration in a transparent and reproducible manner.

`MultiOmicsCancer` fills a gap between complex statistical frameworks and
simple, single-purpose scripts, providing a middle ground that is both rigorous
and accessible.

# Acknowledgements

The author thanks the Royan Institute for its support during the development
of this package.

# References
