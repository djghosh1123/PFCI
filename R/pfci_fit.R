#' Penalized FCI (PFCI): glasso screening + constrained FCI
#'
#' Runs a two-stage procedure:
#' (1) Graphical lasso screening to obtain a sparse undirected super-skeleton
#' (2) FCI on the restricted search space using fixedGaps and a gated CI test
#'
#' @param X Numeric matrix or data.frame of dimension n x p.
#' @param alpha Significance level for conditional independence tests in FCI.
#' @param rho Graphical lasso penalty. If NULL, uses a default depending on n.
#' @param approx Passed to glasso::glasso.
#' @param skel.method Skeleton method for pcalg::fci (default "stable").
#' @param doPdsep Logical; passed to pcalg::fci. Default FALSE.
#' @param labels Optional variable names (length p). If NULL uses colnames or X1..Xp.
#'
#' @return An object of class \code{pfci_fit}, a list containing:
#' \describe{
#'   \item{amat}{Adjacency matrix of the estimated PAG (integer codes: 0=none,
#'     1=circle, 2=arrowhead, 3=tail).}
#'   \item{pag}{The raw \code{fci} output object from \pkg{pcalg}.}
#'   \item{screen_adj}{Logical adjacency matrix from the glasso screening step.}
#'   \item{fixedGaps}{Logical matrix of fixed gaps passed to FCI.}
#'   \item{rho}{The glasso penalty used.}
#'   \item{alpha}{The significance level used.}
#'   \item{time}{A list with \code{glasso}, \code{fci}, and \code{total} runtimes in seconds.}
#' }
#' @references
#' Pal, S., Ghosh, D., and Yang, S. (2025). Penalized FCI for Causal Structure
#' Learning in a Sparse DAG for Biomarker Discovery in Parkinson's Disease.
#' \emph{Annals of Applied Statistics}. \doi{10.48550/arXiv.2507.00173}
#' @seealso \code{\link{pfci_metrics}}, \code{\link{plot_pag}},
#'   \code{\link{simulate_pfci_toy}}
#' @examples
#' \donttest{
#'   sim <- simulate_pfci_toy(p = 30, n = 100, edge_prob = 0.05, seed = 1)
#'   fit <- pfci_fit(sim$X, alpha = 0.05)
#'   print(fit)
#' }
#' @export
pfci_fit <- function(X, alpha = 0.05, rho = NULL, approx = TRUE,
                     skel.method = "stable", doPdsep = FALSE, labels = NULL) {

  if (!requireNamespace("pcalg", quietly = TRUE)) {
    stop(
      "Package 'pcalg' is required for pfci_fit() but is not installed.\n",
      " Please install pcalg from Bioconductor.",
      call. = FALSE
    )
  }

  if (is.data.frame(X)) X <- as.matrix(X)
  stopifnot(is.matrix(X))
  n <- nrow(X); p <- ncol(X)

  if (is.null(labels)) {
    labels <- colnames(X)
    if (is.null(labels)) labels <- paste0("X", seq_len(p))
  }
  colnames(X) <- labels

  if (is.null(rho)) {
    deg <- 2.5
    eps <- log(deg) / log(n) + 0.01
    rho <- n^(-(1 - eps) / 2)
  }

  t0 <- Sys.time()
  S <- stats::cov(X)
  gfit <- glasso::glasso(S, rho = rho, approx = approx)
  adj <- (gfit$wi != 0)
  diag(adj) <- FALSE
  screen_adj <- adj | t(adj)
  fixedGaps <- !screen_adj
  t1 <- Sys.time()

  suffStat <- list(C = stats::cor(X), n = n)

  gate_test <- function(x, y, Sset, suffStat) {
    if (!screen_adj[x, y]) return(1)
    pcalg::gaussCItest(x, y, Sset, suffStat)
  }

  f0 <- Sys.time()
  fit <- pcalg::fci(suffStat = suffStat,
                    indepTest = gate_test,
                    alpha = alpha,
                    labels = labels,
                    fixedGaps = fixedGaps,
                    skel.method = skel.method,
                    doPdsep = doPdsep)
  f1 <- Sys.time()

  out <- list(
    amat = fit@amat,
    pag = fit,
    screen_adj = screen_adj,
    fixedGaps = fixedGaps,
    rho = rho,
    alpha = alpha,
    time = list(
      glasso = as.numeric(difftime(t1, t0, units = "secs")),
      fci    = as.numeric(difftime(f1, f0, units = "secs")),
      total  = as.numeric(difftime(f1, t0, units = "secs"))
    )
  )
  class(out) <- "pfci_fit"
  out
}

#' @export
print.pfci_fit <- function(x, ...) {
  cat("PFCI fit\n")
  cat("  p =", nrow(x$amat), "\n")
  cat("  alpha =", x$alpha, "\n")
  cat("  rho   =", x$rho, "\n")
  cat("  time (sec): glasso =", sprintf("%.3f", x$time$glasso),
      ", fci =", sprintf("%.3f", x$time$fci),
      ", total =", sprintf("%.3f", x$time$total), "\n")
  invisible(x)
}
