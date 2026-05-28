## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* Local: macOS [your version], R [your version]
* R-hub: Windows Server, R-devel
* R-hub: Ubuntu Linux, R-release

## Dependencies

This package uses `pcalg`, `graph`, `RBGL`, and `Rgraphviz` from
Bioconductor, listed under Suggests (not Imports). Functions that
require these packages check for their availability at runtime and
provide clear installation instructions if they are not found.

## Submission notes

* This is the first submission of PFCI to CRAN.
* All examples that call `pfci_fit()` are wrapped in `\donttest{}`
  because they require Bioconductor packages that may not be present
  on all check systems.
* The vignette uses `eval=FALSE` throughout for the same reason.
