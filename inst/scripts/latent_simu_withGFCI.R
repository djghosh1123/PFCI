## ================================
## Packages
## ================================
suppressPackageStartupMessages({
  if (!requireNamespace("pcalg", quietly = TRUE)) install.packages("pcalg")
  if (!requireNamespace("glasso", quietly = TRUE)) install.packages("glasso")
  if (!requireNamespace("MASS", quietly = TRUE))   install.packages("MASS")
  if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
  library(pcalg); library(glasso); library(MASS); library(igraph)
})

## Optional (GFCI):
DO_GFCI <- TRUE
if (DO_GFCI) {
  if (!requireNamespace("rJava", quietly = TRUE)) install.packages("rJava")
  if (!requireNamespace("rcausal", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github("bd2kccd/r-causal")
  }
  library(rJava); library(rcausal)
}

## ================================
## Utilities
## ================================

# Convert Tetrad graph -> pcalg amat.pag (0/1/2/3 ; 1=circle, 2=arrow, 3=tail)
tetrad_pag_to_amat_pag <- function(graph) {
  edf <- rcausal::graphToDataFrame(graph)  # expects node1,node2,endpoint1,endpoint2
  nodes <- unique(c(edf$node1, edf$node2))
  p <- length(nodes); idx <- seq_len(p); names(idx) <- nodes
  amat <- matrix(0L, p, p, dimnames = list(nodes, nodes))
  map <- function(s) if (s=="ARROW") 2L else if (s=="TAIL") 3L else if (s=="CIRCLE") 1L else 0L
  for (k in seq_len(nrow(edf))) {
    a <- edf$node1[k]; b <- edf$node2[k]
    ia <- idx[[a]]; ib <- idx[[b]]
    amat[ia, ib] <- map(edf$endpoint2[k])  # mark at column-endpoint (b)
    amat[ib, ia] <- map(edf$endpoint1[k])  # mark at column-endpoint (a)
  }
  diag(amat) <- 0L
  amat
}

# GFCI (Tetrad) runner
gfci_tetrad <- function(X, alpha=0.05, score="sem-bic", test="fisher-z",
                        penaltyDiscount=2.0, depth=3, maxPathLength=-1,
                        faithfulnessAssumed=TRUE, verbose=FALSE) {
  X <- as.data.frame(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  tdat <- rcausal::loadContinuousData(X)
  runner <- rcausal::TetradRunner()
  runner$setAlgorithm("gfci")
  runner$setData(tdat)
  runner$setScore(score); runner$setTest(test)
  runner$setParameters(list(
    alpha=alpha, penaltyDiscount=penaltyDiscount, depth=depth,
    maxPathLength=maxPathLength, faithfulnessAssumed=faithfulnessAssumed,
    verbose=verbose
  ))
  t0 <- Sys.time(); runner$run(); t1 <- Sys.time()
  list(graph=runner$getTetradGraph(), time_sec=as.numeric(difftime(t1,t0,units="secs")))
}

# Decode a pair (i,j) from amat.pag into an "edge type label"
# returns one of: "dir" (i->j), "parDir" (o->), "bidir" (<->), "und" (o-o), or "none"
pair_type <- function(amat, i, j) {
  a <- amat[i,j]; b <- amat[j,i]
  if (a==0L && b==0L) return("none")
  if (a==2L && b==3L) return("dir")     # i -> j
  if (a==1L && b==3L) return("parDir")  # i o-> j
  if (a==2L && b==2L) return("bidir")   # i <-> j
  if (a==1L && b==1L) return("und")     # i o-o j
  # any other mixed marks are ambiguous; treat as "other" (won't be counted)
  return("other")
}

# Edge-type F1/prec/recall between true and est amat.pag (observed nodes only)
edge_type_metrics <- function(true_amat, est_amat, obs_idx) {
  stopifnot(all(dim(true_amat)==dim(est_amat)))
  m <- list(dir=0, parDir=0, bidir=0, und=0)
  types <- names(m)
  out <- lapply(types, function(tp) {
    TP=0; FP=0; FN=0
    for (i in obs_idx) for (j in obs_idx) if (i<j) {
      ttp <- pair_type(true_amat, i, j)
      etp <- pair_type(est_amat , i, j)
      if (etp == tp && ttp == tp) TP <- TP + 1
      if (etp == tp && ttp != tp) FP <- FP + 1
      if (etp != tp && ttp == tp) FN <- FN + 1
    }
    prec <- if ((TP+FP)>0) TP/(TP+FP) else NA
    rec  <- if ((TP+FN)>0) TP/(TP+FN) else NA
    f1   <- if (!is.na(prec) && !is.na(rec) && (prec+rec)>0) 2*prec*rec/(prec+rec) else NA
    c(precision=prec, recall=rec, f1=f1)
  })
  names(out) <- types
  out
}

# Arrowhead/Tail F1 (aggregate): treat arrowheads and tails presence regardless of the other end
arrow_tail_metrics <- function(true_amat, est_amat, obs_idx) {
  # For each unordered pair (i<j), count whether each end is an arrowhead or a tail
  # Arrowhead at j: true_amat[i,j]==2 ; Tail at j: true_amat[i,j]==3 ; Circle=1 ignored here.
  agg <- function(mark_code) {
    TP=0; FP=0; FN=0
    for (i in obs_idx) for (j in obs_idx) if (i<j) {
      t_ij <- (true_amat[i,j]==mark_code)
      e_ij <- (est_amat [i,j]==mark_code)
      if (t_ij && e_ij) TP <- TP+1 else if (!t_ij && e_ij) FP <- FP+1 else if (t_ij && !e_ij) FN <- FN+1
      t_ji <- (true_amat[j,i]==mark_code)
      e_ji <- (est_amat [j,i]==mark_code)
      if (t_ji && e_ji) TP <- TP+1 else if (!t_ji && e_ji) FP <- FP+1 else if (t_ji && !e_ji) FN <- FN+1
    }
    prec <- if ((TP+FP)>0) TP/(TP+FP) else NA
    rec  <- if ((TP+FN)>0) TP/(TP+FN) else NA
    f1   <- if (!is.na(prec) && !is.na(rec) && (prec+rec)>0) 2*prec*rec/(prec+rec) else NA
    c(precision=prec, recall=rec, f1=f1)
  }
  list(arrow = agg(2L), tail = agg(3L))
}

# Skeleton F1 (undirected presence/absence)
skeleton_f1 <- function(true_amat, est_amat, obs_idx) {
  TP=0; FP=0; FN=0
  has_edge <- function(amat,i,j) (amat[i,j]!=0L || amat[j,i]!=0L)
  for (i in obs_idx) for (j in obs_idx) if (i<j) {
    t <- has_edge(true_amat, i,j); e <- has_edge(est_amat, i,j)
    if (e && t) TP <- TP+1 else if (e && !t) FP <- FP+1 else if (!e && t) FN <- FN+1
  }
  prec <- if ((TP+FP)>0) TP/(TP+FP) else NA
  rec  <- if ((TP+FN)>0) TP/(TP+FN) else NA
  f1   <- if (!is.na(prec) && !is.na(rec) && (prec+rec)>0) 2*prec*rec/(prec+rec) else NA
  c(precision=prec, recall=rec, f1=f1)
}

## ================================
## Simulation design
## (fast true-PAG construction)
## ================================
# We build a DAG with:
# - observed nodes 1..p_obs
# - latent nodes (p_obs+1)..(p_obs+p_lat)
# - edges among observed chosen uniformly from upper triangle (acyclic)
# - latent nodes have ONLY outgoing edges to observed (random fan-out)
# Then "true PAG":
# - observed i->j if (i->j) is a directed edge in the DAG
# - observed i<->j if they share ≥1 latent parent
# - others: no edge

simulate_once <- function(p_obs=200, p_lat=40, n=100,
                          edge_ratio = 1L,    # total edges among observed = edge_ratio * p_tot
                          latent_out_deg = 3L,# avg #children per latent into observed
                          seed = 1L) {
  set.seed(seed)
  p_tot <- p_obs + p_lat
  obs_idx <- seq_len(p_obs)
  lat_idx <- (p_obs+1):p_tot
  
  # Topological order is 1..p_tot ; we only add edges i->j for i<j
  A <- matrix(0L, p_tot, p_tot)
  
  # (1) observed→observed edges
  # choose m edges from upper triangle among observed nodes
  up <- which(upper.tri(matrix(0,p_obs,p_obs)), arr.ind = TRUE)
  m  <- min(nrow(up), edge_ratio * p_tot)
  sel <- up[sample.int(nrow(up), m, replace = FALSE), , drop=FALSE]
  for (k in seq_len(nrow(sel))) A[sel[k,1], sel[k,2]] <- 1L
  
  # (2) latent→observed edges: each latent picks Poisson(latent_out_deg) children uniformly in observed
  for (L in lat_idx) {
    d <- rpois(1, lambda = latent_out_deg)
    if (d>0) {
      kids <- sample(obs_idx, size = min(d, p_obs), replace = FALSE)
      A[L, kids] <- 1L
    }
  }
  
  # Edge weights and Gaussian simulation in topo order
  W <- matrix(0, p_tot, p_tot)
  W[A==1L] <- rnorm(sum(A), mean=0, sd=0.8)
  topo <- 1:p_tot
  X <- matrix(0, n, p_tot)
  eps_sd <- rep(1, p_tot)
  for (j in topo) {
    pa <- which(A[,j]==1L)
    mu <- if (length(pa)) X[,pa,drop=FALSE] %*% W[pa,j] else 0
    X[,j] <- mu + rnorm(n, 0, eps_sd[j])
  }
  
  colnames(X) <- paste0("X", 1:p_tot)
  
  # Build "true PAG" amat.pag over ALL nodes first, then we’ll slice to observed later
  # Start as no edges
  amat_true <- matrix(0L, p_tot, p_tot)
  
  # observed directed edges kept as arrows
  for (i in obs_idx) for (j in obs_idx) if (A[i,j]==1L) {
    amat_true[i,j] <- 2L  # arrowhead at j
    amat_true[j,i] <- 3L  # tail at i
  }
  
  # bidirected among observed if they share at least one latent parent
  for (i in obs_idx) for (j in obs_idx) if (i<j) {
    # check if any latent L with A[L,i]==1 & A[L,j]==1
    shared <- which(A[lat_idx, i]==1L & A[lat_idx, j]==1L)
    if (length(shared)>0) {
      amat_true[i,j] <- 2L
      amat_true[j,i] <- 2L
    }
  }
  
  # Data returned: observed variables only
  list(
    X_obs = X[, obs_idx, drop=FALSE],
    amat_true_full = amat_true,
    obs_idx = obs_idx
  )
}

## ================================
## Estimators
## ================================

# RFCI (pcalg) on observed only
run_rfci <- function(X) {
  suff <- list(C = cor(X), n = nrow(X))
  t0 <- Sys.time()
  fit <- pcalg::rfci(suffStat = suff, indepTest = gaussCItest,
                     alpha = 0.05, labels = colnames(X),
                     skel.method = "stable")
  t1 <- Sys.time()
  list(amat = fit@amat, time_sec = as.numeric(difftime(t1,t0,units="secs")))
}

# PFCI = glasso skeleton + constrained FCI with gate
run_pfci <- function(X, rho = NULL) {
  S <- cov(X); n <- nrow(X)
  if (is.null(rho)) {
    # simple theory-guided default tuned by avg degree target ~ 2.5
    deg <- 2.5; eps <- log(deg)/log(n) + 0.01; rho <- n^(-(1 - eps)/2)
  }
  t0 <- Sys.time()
  gfit <- glasso::glasso(S, rho = rho, approx = TRUE)
  adj <- (gfit$wi != 0); diag(adj) <- FALSE
  adj_sym <- (adj | t(adj))
  fixedGaps <- !adj_sym  # forbid edges not in glasso adjacency (symmetric)
  suff <- list(C = cor(X), n = nrow(X))
  
  # gate indepTest: if glasso says no edge, return p=1 quickly; else do gaussCItest
  gate <- function(x, y, Sset, suffStat) {
    if (!adj_sym[x,y]) return(1) else return(gaussCItest(x,y,Sset,suffStat))
  }
  
  fit <- pcalg::fci(suffStat = suff, indepTest = gate, alpha = 0.05,
                    labels = colnames(X), fixedGaps = fixedGaps,
                    skel.method = "stable", doPdsep = FALSE)
  t1 <- Sys.time()
  list(amat = fit@amat, time_sec = as.numeric(difftime(t1,t0,units="secs")))
}

# GFCI (Tetrad) on observed only (if Java available)
run_gfci <- function(X) {
  gf <- gfci_tetrad(X, alpha=0.05, score="sem-bic", test="fisher-z",
                    penaltyDiscount=2.0, depth=3)
  amat <- tetrad_pag_to_amat_pag(gf$graph)
  # Ensure same node order/labels as X:
  want <- colnames(X)
  if (!all(want %in% colnames(amat))) {
    stop("GFCI graph nodes do not match data columns.")
  }
  amat <- amat[want, want]
  list(amat = amat, time_sec = gf$time_sec)
}

## ================================
## One-run wrapper (for a table row)
## ================================
run_one <- function(p_obs=200, gamma=0.05, n=100, edge_ratio = 1L,
                    latent_out_deg = 3L, seed=1L,
                    do_gfci = DO_GFCI) {
  
  p_lat <- round(gamma * p_obs)
  sim <- simulate_once(p_obs=p_obs, p_lat=p_lat, n=n,
                       edge_ratio=edge_ratio, latent_out_deg=latent_out_deg,
                       seed=seed)
  X   <- sim$X_obs
  tru <- sim$amat_true_full
  obs <- sim$obs_idx
  
  # truth restricted to observed block
  true_amat_obs <- tru[obs, obs]
  
  # RFCI
  rfci <- run_rfci(X)
  # PFCI
  pfci <- run_pfci(X)
  # GFCI (optional)
  gfci <- if (do_gfci) run_gfci(X) else NULL
  
  # Metrics
  sk_rfci <- skeleton_f1(true_amat_obs, rfci$amat, obs)
  sk_pfci <- skeleton_f1(true_amat_obs, pfci$amat, obs)
  sk_gfci <- if (!is.null(gfci)) skeleton_f1(true_amat_obs, gfci$amat, obs) else c(precision=NA, recall=NA, f1=NA)
  
  et_rfci <- edge_type_metrics(true_amat_obs, rfci$amat, obs)
  et_pfci <- edge_type_metrics(true_amat_obs, pfci$amat, obs)
  et_gfci <- if (!is.null(gfci)) edge_type_metrics(true_amat_obs, gfci$amat, obs) else NULL
  
  at_rfci <- arrow_tail_metrics(true_amat_obs, rfci$amat, obs)
  at_pfci <- arrow_tail_metrics(true_amat_obs, pfci$amat, obs)
  at_gfci <- if (!is.null(gfci)) arrow_tail_metrics(true_amat_obs, gfci$amat, obs) else NULL
  
  out_row <- list(
    p_obs = p_obs,
    gamma = gamma,
    edge_ratio = edge_ratio,
    
    RFCI = list(
      F1_skel = unname(sk_rfci["f1"]),
      F1_dir  = et_rfci$dir["f1"],
      F1_oDir = et_rfci$parDir["f1"],
      F1_bi   = et_rfci$bidir["f1"],
      F1_und  = et_rfci$und["f1"],
      F1_arrow= at_rfci$arrow["f1"],
      F1_tail = at_rfci$tail["f1"],
      time    = rfci$time_sec
    ),
    PFCI = list(
      F1_skel = unname(sk_pfci["f1"]),
      F1_dir  = et_pfci$dir["f1"],
      F1_oDir = et_pfci$parDir["f1"],
      F1_bi   = et_pfci$bidir["f1"],
      F1_und  = et_pfci$und["f1"],
      F1_arrow= at_pfci$arrow["f1"],
      F1_tail = at_pfci$tail["f1"],
      time    = pfci$time_sec
    ),
    GFCI = if (!is.null(gfci)) list(
      F1_skel = unname(sk_gfci["f1"]),
      F1_dir  = et_gfci$dir["f1"],
      F1_oDir = et_gfci$parDir["f1"],
      F1_bi   = et_gfci$bidir["f1"],
      F1_und  = et_gfci$und["f1"],
      F1_arrow= at_gfci$arrow["f1"],
      F1_tail = at_gfci$tail["f1"],
      time    = gfci$time_sec
    ) else NULL
  )
  out_row
}

## ================================
## Example: reproduce a row like your table
## ================================
# Change p_obs and gamma to sweep rows; replicate over seeds to get mean(sd)
set.seed(123)
res_ex <- run_one(p_obs = 200, gamma = 0.05, n = 100,
                  edge_ratio = 2L, latent_out_deg = 3L, seed = 42,
                  do_gfci = DO_GFCI)

str(res_ex, max.level = 1)

## ================================
## Sweep over settings
## ================================
# p_grid   <- c(100,200,300,400,500)
# gamma_set<- c(0.05, 0.20)
# R        <- 20
# tab      <- list()
# idx <- 1
# for (g in gamma_set) {
#   for (p in p_grid) {
#     rows <- replicate(R, run_one(p_obs=p, gamma=g, n=100,
#                                  edge_ratio=sample(c(1L,2L),1),
#                                  latent_out_deg=3L,
#                                  seed=sample.int(1e6,1),
#                                  do_gfci=DO_GFCI), simplify=FALSE)
#     # summarize to mean(sd)
#     pack <- function(method) {
#       m <- do.call(rbind, lapply(rows, function(z) unlist(z[[method]])))
#       apply(m, 2, function(col) sprintf("%.3f (%.2f)", mean(col, na.rm=TRUE), sd(col, na.rm=TRUE)))
#     }
#     RFCI <- pack("RFCI"); PFCI <- pack("PFCI")
#     GFCI <- if (DO_GFCI) pack("GFCI") else NULL
#     tab[[idx]] <- list(p_obs=p, gamma=g, RFCI=RFCI, PFCI=PFCI, GFCI=GFCI)
#     idx <- idx + 1
#   }
# }
# tab  # format into your LaTeX table as needed
