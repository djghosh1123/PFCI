ReadMe
================
2026-03-14

# PFCI: Penalized Fast Causal Inference for High-Dimensional Structure Learning

PFCI implements **Penalized Fast Causal Inference (PFCI)**, a scalable
two-stage procedure for learning graphical structures in
high-dimensional settings with potential latent variables and selection
bias.

The method combines:

- **Graphical lasso screening** to obtain a sparse super-skeleton
- **Constrained Fast Causal Inference (FCI)** for orientation and
  refinement

This enables computationally efficient structure learning while
preserving theoretical guarantees under sparsity assumptions.

------------------------------------------------------------------------

## Installation

Install from CRAN:

``` r
install.packages("PFCI")
```

The development version is available on GitHub:

``` r
devtools::install_github("djghosh1123/PFCI")
```

Core functionality requires pcalg and graph from Bioconductor:

``` r
install.packages("BiocManager")
BiocManager::install(c("pcalg", "graph", "RBGL", "Rgraphviz"))
```

------------------------------------------------------------------------

## Basic usage

``` r
library(PFCI)

sim <- simulate_pfci_toy(p = 100, n = 100, edge_prob = 0.02, seed = 1)
fit <- pfci_fit(sim$X, alpha = 0.05)
met <- pfci_metrics(sim, fit)
met
plot_pag(fit)
```

------------------------------------------------------------------------

## Reference

Pal, S., Ghosh, D., & Yang, S. (2025). Penalized FCI for Causal
Structure Learning in a Sparse DAG for Biomarker Discovery in
Parkinson’s Disease. *Annals of Applied Statistics*.
<doi:10.48550/arXiv.2507.00173>
