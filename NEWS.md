# PFCI 0.1.1

* Removed `install.packages()` calls from `inst/scripts/` files
  to comply with CRAN policies.
* Fixed writing path in `inst/scripts/run_simulations_full.R`
  to use `tempdir()` instead of `getwd()`.

# PFCI 0.1.0

* Initial CRAN submission.
* Implements the two-stage Penalized Fast Causal Inference (PFCI) algorithm
  combining graphical lasso screening with FCI orientation.
* Core functions: `pfci_fit()`, `pfci_metrics()`, `plot_pag()`,
  `simulate_pfci_toy()`, `simulate_with_latent()`, `metrics_with_latent()`.
* Accompanies Pal, Ghosh, and Yang (2025), Annals of Applied Statistics,
  <doi:10.48550/arXiv.2507.00173>.
