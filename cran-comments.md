## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* Local: Windows, R 4.4.2
* GitHub Actions: windows-latest (R-release)
* GitHub Actions: macOS-latest (R-release)
* GitHub Actions: ubuntu-latest (R-devel, R-release, R-oldrel-1)
* win-builder: R-release and R-devel

## Resubmission notes (0.1.0 -> 0.1.1)

* Commented out `install.packages()` calls in `inst/scripts/` files
  to comply with CRAN policies.
* Replaced `getwd()` with `tempdir()` in `inst/scripts/run_simulations_full.R`.

## win-builder NOTE

* `'New submission'`: expected.
* `'Possibly misspelled words'`: FCI, PFCI, RFCI, and PAG are
  standard causal inference acronyms. Ghosh is an author surname.
  None are misspellings.

## Dependencies

This package uses `pcalg`, `graph`, `RBGL`, and `Rgraphviz` from
Bioconductor, listed under `Suggests` (not `Imports`). Functions that
require these packages check for their availability at runtime and
provide clear installation instructions if they are not found.

## Submission notes

* All examples that call `pfci_fit()` are wrapped in `\donttest{}`
  because they require Bioconductor packages that may not be
  present on all check systems.
* The vignette uses `eval=FALSE` throughout for the same reason.
