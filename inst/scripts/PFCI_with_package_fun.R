################################################################################
## PFCI-only simulation script
## - Reports: SHD (skeleton), F1 (skeleton), MCC (skeleton), runtime (glasso/fci/total)
################################################################################

## ================================
## Packages
## ================================
suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  # BiocManager::install("Rgraphviz") # only if you truly need it
  if (!requireNamespace("graph", quietly = TRUE))  install.packages("graph")
  if (!requireNamespace("glasso", quietly = TRUE)) install.packages("glasso")
  if (!requireNamespace("pcalg", quietly = TRUE))  install.packages("pcalg")
  if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")

  library(graph)
  library(glasso)
  library(pcalg)
  library(igraph)
})

## ================================
## DAG generator
## ================================
generate_general_dag <- function(n_nodes = 100,
                                 edge_prob = 0.2,
                                 sparsity = 100) {

  dag_matrix <- matrix(0, nrow = n_nodes, ncol = n_nodes)

  s_nodes <- sort(sample(n_nodes, sparsity, replace = FALSE))

  for (i in seq_along(s_nodes)) {
    for (j in (i + 1):length(s_nodes)) {
      if (runif(1) < edge_prob) {
        dag_matrix[s_nodes[i], s_nodes[j]] <- 1
      }
    }
  }

  true_dag <- as(dag_matrix, "graphNEL")
  return(list(true_dag = true_dag, adj_mat = dag_matrix))
}

## ================================
## PFCI function (from package)
## ================================
pfci_fit <- function(X, alpha = 0.05, rho = NULL, approx = TRUE,
                     skel.method = "stable", doPdsep = FALSE, labels = NULL) {

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

## ================================
## Helpers: skeleton extraction + metrics
## ================================
amat_to_skeleton <- function(amat) {
  # pcalg PAG amat uses {0,1,2,3}; any nonzero mark means adjacency
  A <- (amat != 0) * 1
  A <- ((A + t(A)) > 0) * 1
  diag(A) <- 0
  A
}

dag_adj_to_skeleton <- function(adj_mat) {
  A <- (adj_mat != 0) * 1
  A <- ((A + t(A)) > 0) * 1
  diag(A) <- 0
  A
}

skeleton_shd <- function(true_skel, est_skel) {
  # SHD for undirected skeletons (counts differing edges once)
  p <- nrow(true_skel)
  diff_count <- 0L
  for (i in 1:(p-1)) for (j in (i+1):p) {
    if (true_skel[i,j] != est_skel[i,j]) diff_count <- diff_count + 1L
  }
  diff_count
}

skeleton_f1_mcc <- function(true_skel, est_skel) {
  # Both are symmetric 0/1 with diag=0
  TP <- sum(true_skel == 1 & est_skel == 1) / 2
  FP <- sum(true_skel == 0 & est_skel == 1) / 2
  FN <- sum(true_skel == 1 & est_skel == 0) / 2
  TN <- sum(true_skel == 0 & est_skel == 0) / 2 - nrow(true_skel)/2  # remove diagonal contribution

  prec <- if ((TP+FP) > 0) TP/(TP+FP) else 0
  rec  <- if ((TP+FN) > 0) TP/(TP+FN) else 0
  f1   <- if ((prec+rec) > 0) 2*prec*rec/(prec+rec) else 0

  denom <- (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)
  mcc <- if (denom > 0) (TP*TN - FP*FN)/sqrt(denom) else 0

  c(F1 = f1, MCC = mcc, Precision = prec, Recall = rec,
    TP = TP, FP = FP, FN = FN, TN = TN)
}

## ================================
## One simulation run
## ================================
run_pfci_simulation_once <- function(p, s, n, edge_prob, errDist,
                                     alpha = 0.05,
                                     approx = TRUE,
                                     skel.method = "stable",
                                     doPdsep = FALSE,
                                     seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  true_dag <- generate_general_dag(n_nodes = p,
                                   sparsity = s,
                                   edge_prob = edge_prob)

  data <- rmvDAG(n, true_dag$true_dag, errDist = errDist)

  fit <- pfci_fit(data,
                  alpha = alpha,
                  approx = approx,
                  skel.method = skel.method,
                  doPdsep = doPdsep)

  true_skel <- dag_adj_to_skeleton(true_dag$adj_mat)
  est_skel  <- amat_to_skeleton(as(fit$amat, "matrix"))

  shd <- skeleton_shd(true_skel, est_skel)
  fm  <- skeleton_f1_mcc(true_skel, est_skel)

  data.frame(
    SHD = shd,
    F1_total = fm["F1"],
    MCC = fm["MCC"],
    time_total = fit$time$total,
    time_glasso = fit$time$glasso,
    time_fci = fit$time$fci,
    rho = fit$rho
  )
}

## ================================
## Run batches (3 scenarios)
## ================================
run_pfci_batch <- function(p, s, n, edge_prob, R = 20,
                           alpha = 0.05,
                           approx = TRUE,
                           skel.method = "stable",
                           doPdsep = FALSE,
                           seed0 = 12345) {

  scenarios <- data.frame(
    scenario = c("Gaussian errors", "t4 heavy tails", "mixt3"),
    errDist  = c("normal", "t4", "mixt3"),
    stringsAsFactors = FALSE
  )

  results <- list()

  counter <- 1

  for (k in seq_len(nrow(scenarios))) {

    sc <- scenarios$scenario[k]
    ed <- scenarios$errDist[k]

    for (r in seq_len(R)) {

      seed_use <- seed0 + 1000*k + r

      row <- run_pfci_simulation_once(
        p = p, s = s, n = n,
        edge_prob = edge_prob,
        errDist = ed,
        alpha = alpha,
        approx = approx,
        skel.method = skel.method,
        doPdsep = doPdsep,
        seed = seed_use
      )

      row$scenario <- sc
      row$rep <- r

      results[[counter]] <- row
      counter <- counter + 1
    }
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}

summarize_pfci <- function(df) {

  aggregate(
    cbind(SHD, F1_total, MCC, time_total) ~ scenario,
    data = df,
    FUN = function(x) {
      c(mean = mean(x), sd = sd(x))
    }
  )
}

## ================================
## Example:
## ================================

## (A) “mixt3” case from your first block: p=500, s=500, n=100, edge_prob=0.005
set.seed(1)
res_A <- run_pfci_batch(p = 500, s = 500, n = 100, edge_prob = 0.005, R = 20)

## (B) “t4” case from your second block: p=700, s=700, n=100, edge_prob=0.005
set.seed(2)
res_B <- run_pfci_batch(p = 700, s = 700, n = 100, edge_prob = 0.005, R = 20)

cat("\n=== Results A (p=500) ===\n")
print(summarize_pfci(res_A))

cat("\n=== Results B (p=700) ===\n")
print(summarize_pfci(res_B))

suppressPackageStartupMessages({
  library(pcalg)
  library(glasso)
})

###############################################################################
## 1) Generate a topologically ordered DAG
##    (edges only from lower index -> higher index)
###############################################################################
generate_general_dag_topo <- function(n_nodes = 100, edge_prob = 0.02, sparsity = 100, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  A <- matrix(0L, n_nodes, n_nodes)

  s_nodes <- sort(sample.int(n_nodes, size = min(sparsity, n_nodes), replace = FALSE))

  # only i<j edges among s_nodes => topo ordered by node index
  for (ii in seq_along(s_nodes)) {
    i <- s_nodes[ii]
    if (ii == length(s_nodes)) next
    for (jj in (ii + 1):length(s_nodes)) {
      j <- s_nodes[jj]
      if (runif(1) < edge_prob) A[i, j] <- 1L
    }
  }

  g <- as(A, "graphNEL")
  list(true_dag = g, adj_mat = A)
}

###############################################################################
## 2) Truth amat from CPDAG: directed edges become tail/arrow, undirected become o-o
## amat codes: 0 none, 1 circle, 2 arrow, 3 tail
###############################################################################
cpdag_graphNEL_to_amat <- function(cpdag) {
  labs <- nodes(cpdag)
  p <- length(labs)
  amat <- matrix(0L, p, p, dimnames = list(labs, labs))

  M <- as(cpdag, "matrix")
  M[M != 0] <- 1

  for (i in 1:(p - 1)) for (j in (i + 1):p) {
    ij <- (M[i, j] == 1)
    ji <- (M[j, i] == 1)
    if (!ij && !ji) next

    if (ij && ji) {
      amat[i, j] <- 1L
      amat[j, i] <- 1L
    } else if (ij && !ji) {
      amat[i, j] <- 2L
      amat[j, i] <- 3L
    } else {
      amat[j, i] <- 2L
      amat[i, j] <- 3L
    }
  }
  diag(amat) <- 0L
  amat
}

true_amat_from_dag <- function(dag_graphNEL) {
  cpdag <- pcalg::dag2cpdag(dag_graphNEL)
  cpdag_graphNEL_to_amat(cpdag)
}

###############################################################################
## 3) Metrics
###############################################################################
safe_prf <- function(TP, FP, FN) {
  prec <- if ((TP + FP) > 0) TP / (TP + FP) else 0
  rec  <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  f1   <- if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0
  c(precision = prec, recall = rec, f1 = f1)
}

has_edge_undirected <- function(amat, i, j) (amat[i, j] != 0L) || (amat[j, i] != 0L)

pair_class <- function(amat, i, j) {
  a <- amat[i, j]  # mark at j
  b <- amat[j, i]  # mark at i
  if (a == 0L && b == 0L) return("none")
  if ((a == 2L && b == 3L) || (a == 3L && b == 2L)) return("dir")   # ->
  if ((a == 2L && b == 1L) || (a == 1L && b == 2L)) return("oDir")  # o->
  if (a == 2L && b == 2L) return("bidir")                           # <->
  if (a == 1L && b == 1L) return("circ")                            # o-o
  "other"
}

skeleton_f1 <- function(true_amat, est_amat) {
  p <- nrow(true_amat)
  TP <- FP <- FN <- 0L
  for (i in 1:(p - 1)) for (j in (i + 1):p) {
    t <- has_edge_undirected(true_amat, i, j)
    e <- has_edge_undirected(est_amat , i, j)
    if (e && t)  TP <- TP + 1L
    if (e && !t) FP <- FP + 1L
    if (!e && t) FN <- FN + 1L
  }
  safe_prf(TP, FP, FN)[["f1"]]
}

class_f1 <- function(true_amat, est_amat, cls = c("dir", "oDir", "bidir", "circ")) {
  cls <- match.arg(cls)
  p <- nrow(true_amat)
  TP <- FP <- FN <- 0L
  for (i in 1:(p - 1)) for (j in (i + 1):p) {
    t <- (pair_class(true_amat, i, j) == cls)
    e <- (pair_class(est_amat , i, j) == cls)
    if (e && t)  TP <- TP + 1L
    if (e && !t) FP <- FP + 1L
    if (!e && t) FN <- FN + 1L
  }
  safe_prf(TP, FP, FN)[["f1"]]
}

endpoint_f1 <- function(true_amat, est_amat, mark = c("arrow", "tail")) {
  mark <- match.arg(mark)
  code <- if (mark == "arrow") 2L else 3L
  p <- nrow(true_amat)
  TP <- FP <- FN <- 0L
  for (i in 1:p) for (j in 1:p) if (i != j) {
    t <- (true_amat[i, j] == code)
    e <- (est_amat [i, j] == code)
    if (e && t)  TP <- TP + 1L
    if (e && !t) FP <- FP + 1L
    if (!e && t) FN <- FN + 1L
  }
  safe_prf(TP, FP, FN)[["f1"]]
}

###############################################################################
## 4) One replicate (PFCI only)
###############################################################################
run_one_pfci_marks <- function(p = 100, s = 100, n = 100, edge_prob = 0.02,
                               errDist = "normal", alpha = 0.05,
                               approx = TRUE, seed = 1) {
  set.seed(seed)

  dag_obj <- generate_general_dag_topo(n_nodes = p, edge_prob = edge_prob, sparsity = s)
  dag <- dag_obj$true_dag

  # DATA GENERATION: keep rmvDAG exactly (this requires topo order, now satisfied)
  X <- pcalg::rmvDAG(n, dag, errDist = errDist)

  # PFCI (your exact function)
  fit <- pfci_fit(X, alpha = alpha, approx = approx)

  # TRUTH
  amat_true <- true_amat_from_dag(dag)
  amat_est  <- fit$amat

  data.frame(
    scenario = errDist,
    F1_score = skeleton_f1(amat_true, amat_est),
    `->`     = class_f1(amat_true, amat_est, "dir"),
    `o->`    = class_f1(amat_true, amat_est, "oDir"),
    `<->`    = class_f1(amat_true, amat_est, "bidir"),
    `o-o`    = class_f1(amat_true, amat_est, "circ"),
    Arrow    = endpoint_f1(amat_true, amat_est, "arrow"),
    Tail     = endpoint_f1(amat_true, amat_est, "tail"),
    Time     = as.numeric(fit$time$total),
    check.names = FALSE
  )
}

###############################################################################
## 5) Table runner + summary
###############################################################################
fmt <- function(x) sprintf("%.3f (%.2f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))

run_pfci_table <- function(p = 100, s = 100, n = 100, edge_prob = 0.02,
                           R = 10, seed0 = 1, alpha = 0.05, approx = TRUE) {

  scenarios <- c(normal = "normal", t4 = "t4", mixt3 = "mixt3")

  res <- do.call(rbind, lapply(names(scenarios), function(sc_name) {
    err <- scenarios[[sc_name]]
    do.call(rbind, lapply(seq_len(R), function(r) {
      run_one_pfci_marks(
        p = p, s = s, n = n, edge_prob = edge_prob,
        errDist = err, alpha = alpha, approx = approx,
        seed = seed0 + 1000L * match(sc_name, names(scenarios)) + r
      )
    }))
  }))

  summ <- aggregate(
    res[, c("F1_score","->","o->","<->","o-o","Arrow","Tail","Time"), drop = FALSE],
    by = list(scenario = res$scenario),
    FUN = fmt
  )
  summ
}

###############################################################################
## 6) RUN
###############################################################################
set.seed(1)
tab <- run_pfci_table(p = 100, s = 100, n = 100, edge_prob = 0.02, R = 10, seed0 = 1)
print(tab)
