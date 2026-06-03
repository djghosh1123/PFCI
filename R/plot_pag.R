#' Plot a PAG returned by PFCI
#'
#' Plots the Partial Ancestral Graph (PAG) estimated by \code{\link{pfci_fit}}
#' using the \pkg{pcalg} plot method. Requires \pkg{Rgraphviz} to be installed.
#'
#' @param fit A \code{pfci_fit} object returned by \code{\link{pfci_fit}},
#'   or a raw \code{pcalg} fci object.
#' @param ... Additional arguments passed to the \pkg{pcalg} plot method.
#' @return Invisibly returns \code{NULL}. Called for its side effect of
#'   producing a graph plot.
#' @seealso \code{\link{pfci_fit}}
#' @examples
#' \donttest{
#'   sim <- simulate_pfci_toy(p = 20, n = 100, edge_prob = 0.05, seed = 1)
#'   fit <- pfci_fit(sim$X, alpha = 0.05)
#'   plot_pag(fit)
#' }
#' @export

plot_pag <- function(fit, ...) {

  if (!requireNamespace("pcalg", quietly = TRUE)) {
    stop(
      "Package 'pcalg' is required for plot_pag() but is not installed.\n",
      "Please install 'pcalg' from Bioconductor.",
      call. = FALSE
    )
  }

  if (!requireNamespace("Rgraphviz", quietly = TRUE)) {
    stop(
      "Package 'Rgraphviz' is required for plot_pag() but is not installed.\n",
      "Please install Rgraphviz from Bioconductor.",
      call. = FALSE
    )
  }

  if (inherits(fit, "pfci_fit")) {
    pag_obj <- fit$pag
  } else {
    pag_obj <- fit
  }

  # Explicitly call pcalg plot method
  pcalg::plot(pag_obj, ...)
}
