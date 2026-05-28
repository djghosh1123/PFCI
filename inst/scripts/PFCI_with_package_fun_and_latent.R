suppressPackageStartupMessages({
  library(pcalg)
  library(glasso)
})

## ----------------------------
## 1) PFCI function (from package)
## ----------------------------
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

## ----------------------------
## 2) Data generation
## ----------------------------
rerr <- function(n, dist = c("normal","t4","mixt3"), mix = 0.05) {
  dist <- match.arg(dist)
  if (dist == "normal") {
    e <- rnorm(n)
  } else if (dist == "t4") {
    e <- rt(n, df = 4) / sqrt(2)
  } else {
    is_out <- rbinom(n, 1, mix) == 1
    e <- rnorm(n)
    e[is_out] <- rt(sum(is_out), df = 3) / sqrt(3)
  }
  as.numeric(scale(e))
}

make_toy_sem <- function(p_obs = 100, gamma = 0.05,
                         edge_prob_obs = 0.02, latent_out_deg = 3,
                         w_sd = 0.8, seed = 1) {
  set.seed(seed)
  p_lat <- max(1, round(gamma * p_obs))
  p_tot <- p_obs + p_lat
  obs <- 1:p_obs
  lat <- (p_obs + 1):p_tot

  A <- matrix(0L, p_tot, p_tot)

  for (i in 1:(p_obs-1)) for (j in (i+1):p_obs) {
    if (runif(1) < edge_prob_obs) A[i, j] <- 1L
  }

  for (L in lat) {
    d <- rpois(1, lambda = latent_out_deg)
    if (d > 0) {
      kids <- sample(obs, size = min(d, p_obs), replace = FALSE)
      A[L, kids] <- 1L
    }
  }

  W <- matrix(0, p_tot, p_tot)
  W[A == 1L] <- rnorm(sum(A), mean = 0, sd = w_sd)

  list(A=A, W=W, obs=obs, lat=lat, p_obs=p_obs, p_lat=p_lat, p_tot=p_tot)
}

simulate_data <- function(sem, n = 100, errDist = c("normal","t4","mixt3"),
                          mix = 0.05, noise_sd = 1, seed = 1) {
  set.seed(seed)
  errDist <- match.arg(errDist)

  A <- sem$A; W <- sem$W
  p_tot <- sem$p_tot
  obs <- sem$obs

  X <- matrix(0, n, p_tot)
  topo <- 1:p_tot

  for (j in topo) {
    pa <- which(A[, j] == 1L)
    mu <- if (length(pa) > 0) X[, pa, drop=FALSE] %*% W[pa, j] else rep(0, n)
    eps <- noise_sd * rerr(n, dist = errDist, mix = mix)
    X[, j] <- as.numeric(mu) + eps
  }

  X_obs <- X[, obs, drop=FALSE]
  colnames(X_obs) <- paste0("X", obs)
  X_obs
}

## ----------------------------
## 3) Metrics
## ----------------------------
safe_prf <- function(TP, FP, FN) {
  prec <- if ((TP + FP) > 0) TP/(TP+FP) else 0
  rec  <- if ((TP + FN) > 0) TP/(TP+FN) else 0
  f1   <- if ((prec + rec) > 0) 2*prec*rec/(prec+rec) else 0
  c(precision=prec, recall=rec, f1=f1)
}

has_edge_undirected <- function(amat, i, j) {
  (amat[i,j] != 0L) || (amat[j,i] != 0L)
}

pair_type <- function(amat, i, j) {
  a <- amat[i,j]
  b <- amat[j,i]
  if (a==0L && b==0L) return("none")
  if (a==2L && b==3L) return("dir_ij")    # i -> j
  if (a==3L && b==2L) return("dir_ji")    # j -> i
  if (a==2L && b==1L) return("parDir_ij") # i o-> j
  if (a==1L && b==2L) return("parDir_ji") # j o-> i
  if (a==2L && b==2L) return("bidir")     # i <-> j
  if (a==1L && b==1L) return("und")       # i o-o j
  "other"
}

skeleton_f1 <- function(trueA, estA) {
  p <- nrow(trueA)
  TP <- FP <- FN <- 0L
  for (i in 1:(p-1)) for (j in (i+1):p) {
    t <- has_edge_undirected(trueA, i, j)
    e <- has_edge_undirected(estA , i, j)
    if (e && t) TP <- TP+1L
    if (e && !t) FP <- FP+1L
    if (!e && t) FN <- FN+1L
  }
  unname(safe_prf(TP, FP, FN)["f1"])
}

edge_class_f1 <- function(trueA, estA, cls=c("dir","parDir","bidir","und")) {
  cls <- match.arg(cls)
  p <- nrow(trueA)
  TP <- FP <- FN <- 0L

  is_cls <- function(amat, i, j) {
    tp <- pair_type(amat, i, j)
    if (cls=="dir")    return(tp %in% c("dir_ij","dir_ji"))
    if (cls=="parDir") return(tp %in% c("parDir_ij","parDir_ji"))
    if (cls=="bidir")  return(tp=="bidir")
    if (cls=="und")    return(tp=="und")
    FALSE
  }

  for (i in 1:(p-1)) for (j in (i+1):p) {
    t <- is_cls(trueA, i, j)
    e <- is_cls(estA , i, j)
    if (e && t)  TP <- TP+1L
    if (e && !t) FP <- FP+1L
    if (!e && t) FN <- FN+1L
  }
  unname(safe_prf(TP, FP, FN)["f1"])
}

endpoint_mark_f1 <- function(trueA, estA, mark=c("arrow","tail")) {
  mark <- match.arg(mark)
  code <- if (mark=="arrow") 2L else 3L
  p <- nrow(trueA)
  TP <- FP <- FN <- 0L
  for (i in 1:p) for (j in 1:p) if (i != j) {
    t <- (trueA[i,j] == code)
    e <- (estA [i,j] == code)
    if (e && t)  TP <- TP+1L
    if (e && !t) FP <- FP+1L
    if (!e && t) FN <- FN+1L
  }
  unname(safe_prf(TP, FP, FN)["f1"])
}

all_marks <- function(trueA, estA, time_sec) {
  c(
    F1_score = skeleton_f1(trueA, estA),
    `->`     = edge_class_f1(trueA, estA, "dir"),
    `o->`    = edge_class_f1(trueA, estA, "parDir"),
    `<->`    = edge_class_f1(trueA, estA, "bidir"),
    `o-o`    = edge_class_f1(trueA, estA, "und"),
    Arrow    = endpoint_mark_f1(trueA, estA, "arrow"),
    Tail     = endpoint_mark_f1(trueA, estA, "tail"),
    Time     = time_sec
  )
}

## ----------------------------
## 4) One replicate:
## ----------------------------
get_truth_amat <- function(X, alpha) {
  suff <- list(C = cor(X), n = nrow(X))

  ## Try RFCI
  rf <- try(pcalg::rfci(suffStat = suff, indepTest = pcalg::gaussCItest,
                        alpha = alpha, labels = colnames(X),
                        skel.method = "stable"),
            silent = TRUE)
  if (!inherits(rf, "try-error")) return(rf@amat)

  ## Fallback: FCI
  fc <- pcalg::fci(suffStat = suff, indepTest = pcalg::gaussCItest,
                   alpha = alpha, labels = colnames(X),
                   skel.method = "stable", doPdsep = FALSE)
  fc@amat
}

run_one_pfci <- function(p_obs=100, gamma=0.05, n=100,
                                 edge_prob_obs=0.02, latent_out_deg=3,
                                 errDist=c("normal","t4","mixt3"),
                                 alpha=0.05, approx=TRUE,
                                 w_sd=0.8, noise_sd=1,
                                 seed_graph=1, seed_data=2) {
  errDist <- match.arg(errDist)

  sem <- make_toy_sem(
    p_obs=p_obs, gamma=gamma,
    edge_prob_obs=edge_prob_obs,
    latent_out_deg=latent_out_deg,
    w_sd=w_sd, seed=seed_graph
  )
  X <- simulate_data(sem, n=n, errDist=errDist, noise_sd=noise_sd, seed=seed_data)

  truthA <- get_truth_amat(X, alpha=alpha)

  pf <- pfci_fit(X, alpha=alpha, approx=approx, labels=colnames(X))
  estA <- pf$amat

  all_marks(truthA, estA, time_sec = pf$time$total)
}

## ----------------------------
## 5) Table runner
## ----------------------------
run_pfci_table_vs_rfci <- function(p_obs=100, gamma=0.05, n=100,
                                   edge_prob_obs=0.02, latent_out_deg=3,
                                   R=30, seed0=1,
                                   alpha=0.05, approx=TRUE,
                                   w_sd=0.8, noise_sd=1) {

  scenarios <- c(normal="normal", t4="t4", mixt3="mixt3")
  metric_names <- c("F1_score","->","o->","<->","o-o","Arrow","Tail","Time")

  fmt <- function(x) sprintf("%.3f (%.2f)", mean(x), sd(x))

  raw <- list()

  for (sc_name in names(scenarios)) {
    err <- scenarios[[sc_name]]

    mat <- matrix(NA_real_, nrow=R, ncol=length(metric_names))
    colnames(mat) <- metric_names

    ## fixed graph (like a paper table); data seed varies across reps
    seed_graph <- seed0 + 999

    for (r in seq_len(R)) {
      seed_data <- seed0 + 1000L*r + match(sc_name, names(scenarios))*10L

      out <- run_one_pfci(
        p_obs=p_obs, gamma=gamma, n=n,
        edge_prob_obs=edge_prob_obs, latent_out_deg=latent_out_deg,
        errDist=err,
        alpha=alpha, approx=approx,
        w_sd=w_sd, noise_sd=noise_sd,
        seed_graph=seed_graph, seed_data=seed_data
      )

      mat[r, ] <- as.numeric(out[metric_names])
    }

    raw[[sc_name]] <- mat
  }

  tab <- do.call(rbind, lapply(names(raw), function(sc) {
    mat <- raw[[sc]]
    data.frame(
      scenario = sc,
      F1_score = fmt(mat[,"F1_score"]),
      `->`     = fmt(mat[,"Arrow"]),
      `o->`    = fmt(mat[,"o->"]),
      `<->`    = fmt(mat[,"<->"]),
      `o-o`    = fmt(mat[,"o-o"]),
      Time     = fmt(mat[,"Time"]),
      check.names = FALSE,
      row.names = NULL
    )
  }))

  attr(tab, "raw") <- raw
  tab
}

## ============================================================
## RUN
## ============================================================
set.seed(1)

tab <- run_pfci_table_vs_rfci(
  p_obs=200, gamma=0.2, n=100,
  edge_prob_obs=0.02, latent_out_deg=3,
  R=3, seed0=1,
  alpha=0.05, approx=TRUE,
  w_sd=1, noise_sd=1
)

print(tab)
