#' Simulate data with latent variables and oracle-FCI truth skeleton
#'
#' This follows the exact latent SEM + oracle truth scheme:
#'  - Build a DAG over (observed + latent) nodes with:
#'      * observed->observed edges only for i<j (acyclic)
#'      * latent->observed edges (Poisson out-degree)
#'  - Simulate data from linear SEM with chosen error distribution
#'  - Construct "truth" by running FCI on the ORACLE correlation of observed nodes
#'    using a very large virtual sample size and alpha_truth (oracle-ish),
#'    with m.max controlling speed (e.g., m.max = 2)
#'
#' The returned truth is the *skeleton* implied by the oracle-FCI PAG (not marks).
#'
#' @param p_obs Number of observed variables.
#' @param gamma Latent ratio; p_lat = max(1, round(gamma * p_obs)).
#' @param n Sample size.
#' @param edge_prob_obs Edge probability among observed nodes (i<j only).
#' @param latent_out_deg Mean outgoing degree for each latent to observed (Poisson).
#' @param w_sd SD of nonzero edge weights.
#' @param errDist Error distribution for SEM noise: "normal", "t4", "mixt3".
#' @param noise_sd Noise SD multiplier.
#' @param mix Mixing proportion for "mixt3" heavy tail component.
#' @param seed_graph Seed controlling graph + weights.
#' @param seed_data Seed controlling data noise draws.
#' @param truth_alpha Alpha for oracle-truth FCI (typical: 0.9999).
#' @param truth_mmax Maximum conditioning set size in oracle FCI (speed knob; e.g., 2).
#' @param truth_verbose Logical; verbose output from oracle FCI.
#'
#' @return A list with elements: X, truth (skel + amat), meta, sem (A,W,indices).
#' @seealso \code{\link{pfci_fit}}, \code{\link{metrics_with_latent}}
#' @examples
#' \donttest{
#'   sim <- simulate_with_latent(p_obs = 30, gamma = 0.05, n = 100, seed_graph = 1)
#'   str(sim$truth)
#' }
#' @export
simulate_with_latent <- function(
    p_obs = 100, gamma = 0.05, n = 100,
    edge_prob_obs = 0.02, latent_out_deg = 3,
    w_sd = 0.8,
    errDist = c("normal", "t4", "mixt3"),
    noise_sd = 1,
    mix = 0.05,
    seed_graph = 1,
    seed_data  = 2,
    truth_alpha = 0.9999,
    truth_mmax  = 2,
    truth_verbose = FALSE
) {
  errDist <- match.arg(errDist)


  if (!requireNamespace("pcalg", quietly = TRUE)) {
    stop(
      "Package 'pcalg' is required for simulate_with_latent() but is not installed.\n",
      "Please install 'pcalg' from Bioconductor.",
      call. = FALSE
    )
  }


  sem <- .make_toy_sem(
    p_obs = p_obs, gamma = gamma,
    edge_prob_obs = edge_prob_obs,
    latent_out_deg = latent_out_deg,
    w_sd = w_sd,
    seed_graph = seed_graph
  )

  X <- .simulate_data_sem(
    sem = sem, n = n, errDist = errDist,
    noise_sd = noise_sd, mix = mix,
    seed_data = seed_data
  )

  # oracle truth PAG over observed nodes (amat with 0/1/2/3 marks)
  truth_amat <- .truth_pag_fci_oracle(
    sem = sem,
    alpha_truth = truth_alpha,
    noise_sd = noise_sd,
    m.max = truth_mmax,
    verbose = truth_verbose
  )

  truth_skel <- .amat_to_skeleton(truth_amat)

  out <- list(
    X = X,
    truth = list(
      amat = truth_amat,
      skel = truth_skel
    ),
    meta = list(
      p_obs = p_obs, gamma = gamma, n = n,
      edge_prob_obs = edge_prob_obs,
      latent_out_deg = latent_out_deg,
      w_sd = w_sd,
      errDist = errDist,
      noise_sd = noise_sd,
      mix = mix,
      seed_graph = seed_graph,
      seed_data = seed_data,
      truth_alpha = truth_alpha,
      truth_mmax = truth_mmax
    ),
    sem = sem
  )

  class(out) <- c("pfci_latent_sim", class(out))
  out
}

# -----------------------
# INTERNAL HELPERS (NOT EXPORTED)
# -----------------------

.make_toy_sem <- function(
    p_obs, gamma, edge_prob_obs, latent_out_deg, w_sd, seed_graph
) {
  set.seed(seed_graph)

  p_lat <- max(1L, as.integer(round(gamma * p_obs)))
  p_tot <- p_obs + p_lat

  obs <- seq_len(p_obs)
  lat <- (p_obs + 1L):p_tot

  A <- matrix(0L, p_tot, p_tot)

  # observed->observed edges only i<j => acyclic/topo by index
  if (p_obs >= 2L) {
    for (i in 1L:(p_obs - 1L)) for (j in (i + 1L):p_obs) {
      if (stats::runif(1) < edge_prob_obs) A[i, j] <- 1L
    }
  }

  # latent->observed edges
  for (L in lat) {
    d <- stats::rpois(1, lambda = latent_out_deg)
    if (d > 0L) {
      kids <- sample(obs, size = min(d, p_obs), replace = FALSE)
      A[L, kids] <- 1L
    }
  }

  W <- matrix(0, p_tot, p_tot)
  if (sum(A) > 0L) {
    W[A == 1L] <- stats::rnorm(sum(A), mean = 0, sd = w_sd)
  }

  list(
    A = A, W = W,
    obs = obs, lat = lat,
    p_obs = p_obs, p_lat = p_lat, p_tot = p_tot
  )
}

.rerr <- function(n, dist = c("normal", "t4", "mixt3"), mix = 0.05) {
  dist <- match.arg(dist)
  if (dist == "normal") {
    e <- stats::rnorm(n)
  } else if (dist == "t4") {
    # scaled to unit-ish variance: Var(t4)=2 so divide by sqrt(2)
    e <- stats::rt(n, df = 4) / sqrt(2)
  } else {
    is_out <- (stats::rbinom(n, 1, mix) == 1)
    e <- stats::rnorm(n)
    if (any(is_out)) e[is_out] <- stats::rt(sum(is_out), df = 3) / sqrt(3)
  }
  as.numeric(scale(e))
}

.simulate_data_sem <- function(
    sem, n, errDist, noise_sd, mix, seed_data
) {
  set.seed(seed_data)

  A <- sem$A
  W <- sem$W
  p_tot <- sem$p_tot

  X <- matrix(0, n, p_tot)
  topo <- seq_len(p_tot) # valid by construction: parents have smaller index

  for (j in topo) {
    pa <- which(A[, j] == 1L)
    mu <- if (length(pa) > 0L) X[, pa, drop = FALSE] %*% W[pa, j] else rep(0, n)
    eps <- noise_sd * .rerr(n, dist = errDist, mix = mix)
    X[, j] <- as.numeric(mu) + eps
  }

  X_obs <- X[, sem$obs, drop = FALSE]
  colnames(X_obs) <- paste0("X", seq_len(ncol(X_obs)))
  X_obs
}

.oracle_corr_obs <- function(sem, noise_sd = 1) {
  B <- sem$W
  p_tot <- sem$p_tot

  I <- diag(p_tot)
  M <- I - B
  Minv <- solve(M)

  D <- diag(rep(noise_sd^2, p_tot))
  Sigma <- t(Minv) %*% D %*% Minv

  Sobs <- Sigma[sem$obs, sem$obs, drop = FALSE]
  stats::cov2cor(Sobs)
}

.truth_pag_fci_oracle <- function(
    sem, alpha_truth = 0.9999, noise_sd = 1, m.max = 2, verbose = FALSE
) {
  Cobs <- .oracle_corr_obs(sem, noise_sd = noise_sd)
  labels <- paste0("X", seq_len(ncol(Cobs)))

  suff_oracle <- list(C = Cobs, n = 1e9) # virtual sample size (oracle)
  fit <- pcalg::fci(
    suffStat = suff_oracle,
    indepTest = pcalg::gaussCItest,
    alpha = alpha_truth,
    labels = labels,
    skel.method = "stable",
    doPdsep = FALSE,
    m.max = m.max,
    verbose = verbose
  )

  fit@amat
}
