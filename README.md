# MultiOmicsCancer

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23467139.svg)](https://doi.org/10.5281/zenodo.23467139)
[![R-CMD-check](https://github.com/maryam6770/MultiOmicsCancer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/maryam6770/MultiOmicsCancer/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

An R package for integrative multi-omics analysis in cancer research, with a novel **Immune Methylation Score (IMS)** framework and built-in validation tools.

---

## Key Features

### Standard analyses
- Differential methylation (limma)
- Differential expression (DESeq2)
- Multi-omics integration
- Survival analysis (KM + Cox)
- Visualization

### Novel Q1-level tools
- `build_ims()` — Construct an Immune Methylation Score
- `bootstrap_validate_ims()` — Nested bootstrap validation
- `compare_ims_clinical()` — Compare IMS with clinical variables
- `estimate_purity()` — ESTIMATE tumor purity and immune score
- `validate_proteins()` — Protein-level validation via UniProt

---

## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("maryam6770/MultiOmicsCancer")
