## ============================================================================
## PFCI paper simulations (package-based)
## - Uses PFCI::simulate_pfci_toy() and PFCI::pfci_fit()
## - Compares to pcalg::rfci()
## - Optionally compares to GFCI (Tetrad) if rJava + rcausal are installed
## - Returns tidy results + summary table
## ============================================================================

suppressPackageStartupMessages({
  library(PFCI)
  library(pcalg)
  library(glasso)
  library(stats)
  library(progress)   # <- add this
})

## -------------------------------
## Optional GFCI support
## -------------------------------
DO_GFCI <- FALSE  # set TRUE only if you want the GFCI column

has_gfci <- function() {
  DO_GFCI &&
    requireNamespace("rJava", quietly = TRUE) &&
    requireNamespace("rcausal", quietly = TRUE)
}

## Convert Tetrad graph -> pcalg amat.pag encoding (0/1/2/3 ; 1=circle, 2=arrow, 3=tail)
tetrad_pag_to_amat_pag <- function(graph) {
  edf <- rcausal::graphToDataFrame(graph)  # node1,node2,endpoint1,endpoint2
  nodes <- unique(c(edf$node1, edf$node2))
  p <- length(nodes)
  idx <- seq_len(p); names(idx) <- nodes
  amat <- matrix(0L, p, p, dimnames = list(nodes, nodes))

  map <- function(s) {
    if (s == "ARROW") 2L else if (s == "TAIL") 3L else if (s == "CIRCLE") 1L else 0L
  }

  for (k in seq_len(nrow(edf))) {
    a <- edf$node1[k]; b <- edf$node2[k]
    ia <- idx[[a]]; ib <- idx[[b]]
    amat[ia, ib] <- map(edf$endpoint2[k])  # mark at column endpoint (b)
    amat[ib, ia] <- map(edf$endpoint1[k])  # mark at column endpoint (a)
  }
  diag(amat) <- 0L
  amat
}

gfci_tetrad <- function(X, alpha = 0.05,
                        score = "sem-bic", test = "fisher-z",
                        penaltyDiscount = 2.0, depth = 3,
                        maxPathLength = -1, faithfulnessAssumed = TRUE,
                        verbose = FALSE) {
  stopifnot(has_gfci())
  X <- as.data.frame(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))

  tdat <- rcausal::loadContinuousData(X)
  runner <- rcausal::TetradRunner()
  runner$setAlgorithm("gfci")
  runner$setData(tdat)
  runner$setScore(score)
  runner$setTest(test)
  runner$setParameters(list(
    alpha = alpha,
    penaltyDiscount = penaltyDiscount,
    depth = depth,
    maxPathLength = maxPathLength,
    faithfulnessAssumed = faithfulnessAssumed,
    verbose = verbose
  ))

  t0 <- Sys.time()
  runner$run()
  t1 <- Sys.time()

  list(
    graph = runner$getTetradGraph(),
    time_sec = as.numeric(difftime(t1, t0, units = "secs"))
  )
}

run_gfci <- function(X, alpha = 0.05) {
  gf <- gfci_tetrad(X, alpha = alpha, score = "sem-bic", test = "fisher-z",
                    penaltyDiscount = 2.0, depth = 3)
  amat <- tetrad_pag_to_amat_pag(gf$graph)

  want <- colnames(X)
  if (!all(want %in% colnames(amat))) stop("GFCI node labels do not match X colnames.")
  amat <- amat[want, want, drop = FALSE]

  list(amat = amat, time_sec = gf$time_sec)
}

## -------------------------------
## Metrics (PAG / skeleton / endpoints)
## -------------------------------

pair_type <- function(amat, i, j) {
  a <- amat[i, j]; b <- amat[j, i]
  if (a == 0L && b == 0L) return("none")
  if (a == 2L && b == 3L) return("dir")     # i -> j
  if (a == 1L && b == 3L) return("parDir")  # i o-> j
  if (a == 2L && b == 2L) return("bidir")   # i <-> j
  if (a == 1L && b == 1L) return("und")     # i o-o j
  "other"
}

edge_type_metrics <- function(true_amat, est_amat) {
  p <- nrow(true_amat)
  types <- c("dir", "parDir", "bidir", "und")
  out <- setNames(vector("list", length(types)), types)

  for (tp in types) {
    TP <- 0; FP <- 0; FN <- 0
    for (i in seq_len(p)) for (j in seq_len(p)) if (i < j) {
      ttp <- pair_type(true_amat, i, j)
      etp <- pair_type(est_amat, i, j)
      if (etp == tp && ttp == tp) TP <- TP + 1
      if (etp == tp && ttp != tp) FP <- FP + 1
      if (etp != tp && ttp == tp) FN <- FN + 1
    }
    prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    rec  <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) 2 * prec * rec / (prec + rec) else NA_real_
    out[[tp]] <- c(precision = prec, recall = rec, f1 = f1)
  }
  out
}

arrow_tail_metrics <- function(true_amat, est_amat) {
  p <- nrow(true_amat)

  agg <- function(mark_code) {
    TP <- 0; FP <- 0; FN <- 0
    for (i in seq_len(p)) for (j in seq_len(p)) if (i < j) {
      t_ij <- (true_amat[i, j] == mark_code)
      e_ij <- (est_amat [i, j] == mark_code)
      if (t_ij && e_ij) TP <- TP + 1 else if (!t_ij && e_ij) FP <- FP + 1 else if (t_ij && !e_ij) FN <- FN + 1

      t_ji <- (true_amat[j, i] == mark_code)
      e_ji <- (est_amat [j, i] == mark_code)
      if (t_ji && e_ji) TP <- TP + 1 else if (!t_ji && e_ji) FP <- FP + 1 else if (t_ji && !e_ji) FN <- FN + 1
    }
    prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    rec  <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) 2 * prec * rec / (prec + rec) else NA_real_
    c(precision = prec, recall = rec, f1 = f1)
  }

  list(arrow = agg(2L), tail = agg(3L))
}

skeleton_metrics <- function(true_amat, est_amat) {
  p <- nrow(true_amat)
  has_edge <- function(amat, i, j) (amat[i, j] != 0L || amat[j, i] != 0L)

  TP <- 0; FP <- 0; FN <- 0
  for (i in seq_len(p)) for (j in seq_len(p)) if (i < j) {
    t <- has_edge(true_amat, i, j)
    e <- has_edge(est_amat, i, j)
    if (e && t) TP <- TP + 1 else if (e && !t) FP <- FP + 1 else if (!e && t) FN <- FN + 1
  }

  prec <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  rec  <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  f1   <- if (!is.na(prec) && !is.na(rec) && (prec + rec) > 0) 2 * prec * rec / (prec + rec) else NA_real_

  c(precision = prec, recall = rec, f1 = f1)
}

## -------------------------------
## Estimators wrappers
## -------------------------------

run_rfci <- function(X, alpha = 0.05) {
  suff <- list(C = cor(X), n = nrow(X))
  t0 <- Sys.time()
  fit <- pcalg::rfci(suffStat = suff, indepTest = pcalg::gaussCItest,
                     alpha = alpha, labels = colnames(X),
                     skel.method = "stable")
  t1 <- Sys.time()
  list(amat = fit@amat, time_sec = as.numeric(difftime(t1, t0, units = "secs")))
}

run_pfci <- function(X, alpha = 0.05) {
  fit <- PFCI::pfci_fit(X, alpha = alpha)
  list(amat = fit$amat, time_sec = fit$time$total)
}

## -------------------------------
## One replicate, tidy output
## -------------------------------

one_replicate <- function(p_obs = 200, gamma = 0.05, n = 100,
                          edge_ratio = 1L, latent_out_deg = 3L,
                          seed = 1L, alpha = 0.05, do_gfci = FALSE) {

  sim <- PFCI::simulate_pfci_toy(
    p_obs = p_obs, gamma = gamma, n = n,
    edge_ratio = edge_ratio, latent_out_deg = latent_out_deg, seed = seed
  )

  X <- sim$X
  true_amat <- sim$truth_amat

  rf <- run_rfci(X, alpha = alpha)
  pf <- run_pfci(X, alpha = alpha)

  # Optional GFCI
  gf <- NULL
  if (do_gfci && has_gfci()) {
    gf <- run_gfci(X, alpha = alpha)
  }

  pack_method <- function(name, est_amat, time_sec) {
    sk <- skeleton_metrics(true_amat, est_amat)
    et <- edge_type_metrics(true_amat, est_amat)
    at <- arrow_tail_metrics(true_amat, est_amat)

    data.frame(
      method = name,
      p_obs = p_obs,
      gamma = gamma,
      n = n,
      edge_ratio = edge_ratio,
      latent_out_deg = latent_out_deg,
      seed = seed,
      F1_skel  = unname(sk["f1"]),
      F1_dir   = unname(et$dir["f1"]),
      F1_oDir  = unname(et$parDir["f1"]),
      F1_bi    = unname(et$bidir["f1"]),
      F1_und   = unname(et$und["f1"]),
      F1_arrow = unname(at$arrow["f1"]),
      F1_tail  = unname(at$tail["f1"]),
      time_sec = time_sec,
      stringsAsFactors = FALSE
    )
  }

  out <- rbind(
    pack_method("RFCI", rf$amat, rf$time_sec),
    pack_method("PFCI", pf$amat, pf$time_sec)
  )

  if (!is.null(gf)) {
    out <- rbind(out, pack_method("GFCI", gf$amat, gf$time_sec))
  }

  out
}

## -------------------------------
## Sweep runner
## -------------------------------

run_sweep <- function(p_grid = c(100, 200, 300, 400, 500),
                      gamma_set = c(0.05, 0.20),
                      n = 100,
                      R = 20,
                      edge_ratio_choices = c(1L, 2L),
                      latent_out_deg = 3L,
                      alpha = 0.05,
                      base_seed = 123,
                      do_gfci = FALSE) {

  total_runs <- length(p_grid) * length(gamma_set) * R
  pb <- progress::progress_bar$new(
    format = "  Running [:bar] :percent | p=:p gamma=:g rep=:r | eta: :eta",
    total = total_runs,
    clear = FALSE,
    width = 70
  )

  res <- vector("list", total_runs)
  k <- 1L

  for (g in gamma_set) {
    for (p in p_grid) {
      for (r in seq_len(R)) {

        seed_r <- base_seed + 100000L * which(gamma_set == g) +
          1000L * which(p_grid == p) + r

        er <- sample(edge_ratio_choices, 1)

        res[[k]] <- one_replicate(
          p_obs = p, gamma = g, n = n,
          edge_ratio = er,
          latent_out_deg = latent_out_deg,
          seed = seed_r,
          alpha = alpha,
          do_gfci = do_gfci
        )

        pb$tick(tokens = list(p = p, g = g, r = r))
        k <- k + 1L
      }
    }
  }

  do.call(rbind, res)
}


## -------------------------------
## Summaries for paper tables
## -------------------------------

mean_sd <- function(x) sprintf("%.3f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))

summarize_table <- function(df) {
  metric_cols <- c("F1_skel","F1_dir","F1_oDir","F1_bi","F1_und","F1_arrow","F1_tail","time_sec")

  agg <- aggregate(
    df[, metric_cols],
    by = list(method = df$method, p_obs = df$p_obs, gamma = df$gamma),
    FUN = mean_sd
  )

  # Sort for readability
  agg[order(agg$gamma, agg$p_obs, agg$method), ]
}

## -------------------------------
## Main execution
## -------------------------------

# Light sanity run
# tmp <- one_replicate(p_obs = 200, gamma = 0.05, n = 100, seed = 1, do_gfci = DO_GFCI)
# print(tmp)

# Full sweep for paper-like tables
results <- run_sweep(
  p_grid = c(100, 200, 300, 400, 500),
  gamma_set = c(0.05, 0.20),
  n = 100,
  R = 20,
  edge_ratio_choices = c(1L, 2L),
  latent_out_deg = 3L,
  alpha = 0.05,
  base_seed = 123,
  do_gfci = DO_GFCI
)

summary_tab <- summarize_table(results)

# Save outputs
out_dir <- file.path(getwd(), "simulation_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(results, file.path(out_dir, "pfci_sim_results_long.csv"), row.names = FALSE)
write.csv(summary_tab, file.path(out_dir, "pfci_sim_summary_table.csv"), row.names = FALSE)

message("Wrote:\n  ", file.path(out_dir, "pfci_sim_results_long.csv"),
        "\n  ", file.path(out_dir, "pfci_sim_summary_table.csv"))
