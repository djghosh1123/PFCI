#' Metrics for latent simulation using oracle-FCI truth (skeleton only)
#'
#' Designed for the 3-line workflow:
#'   sim <- simulate_with_latent(...)
#'   fit <- pfci_fit(sim$X, ...)
#'   met <- metrics_with_latent(sim, fit)
#'
#' Returns only: SHD, F1_total, MCC, Time.
#'
#' @param sim Output from simulate_with_latent().
#' @param fit Output from pfci_fit() (must contain $amat and $time$total).
#'
#' @return A named list with SHD, F1_total, MCC, Time.
#' @seealso \code{\link{simulate_with_latent}}, \code{\link{pfci_fit}}
#' @examples
#' \donttest{
#'   sim <- simulate_with_latent(p_obs = 30, gamma = 0.05, n = 100, seed_graph = 1)
#'   fit <- pfci_fit(sim$X, alpha = 0.05)
#'   met <- metrics_with_latent(sim, fit)
#'   print(met)
#' }
#' @export

metrics_with_latent <- function(sim, fit) {
  stopifnot(is.list(sim), is.list(sim$truth), is.matrix(sim$truth$skel))
  stopifnot(is.list(fit), !is.null(fit$amat))

  true_skel <- sim$truth$skel
  est_skel <- .amat_to_skeleton(as(fit$amat, "matrix"))
  SHD <- .skel_shd(true_skel, est_skel)
  fm  <- .skel_f1_mcc(true_skel, est_skel)

  # est_skel  <- .amat_to_skeleton_latent(as(fit$amat, "matrix"))
  #
  # SHD <- .skel_shd_latent(true_skel, est_skel)
  # fm  <- .skel_f1_mcc_latent(true_skel, est_skel)

  list(
    SHD = as.numeric(SHD),
    F1_total = unname(fm["F1"]),
    MCC = unname(fm["MCC"]),
    Time = if (!is.null(fit$time$total)) as.numeric(fit$time$total) else NA_real_
  )
}
