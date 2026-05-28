#' @importFrom methods as
#' @importFrom stats runif
NULL

# -------------------------------------------------------------------------
# Internal utility functions shared across PFCI
# NOT exported
# -------------------------------------------------------------------------

.amat_to_skeleton <- function(amat) {
  A <- (amat != 0L) * 1L
  A <- ((A + t(A)) > 0L) * 1L
  diag(A) <- 0L
  A
}

.skel_shd <- function(true_skel, est_skel) {
  p <- nrow(true_skel)
  d <- 0L
  for (i in 1:(p - 1L)) for (j in (i + 1L):p) {
    if (true_skel[i, j] != est_skel[i, j]) d <- d + 1L
  }
  d
}

.skel_f1_mcc <- function(true_skel, est_skel) {
  stopifnot(is.matrix(true_skel), is.matrix(est_skel))
  stopifnot(all(dim(true_skel) == dim(est_skel)))

  ut <- upper.tri(true_skel, diag = FALSE)
  tU <- true_skel[ut]
  eU <- est_skel[ut]

  TP <- sum(tU == 1 & eU == 1, na.rm = TRUE) * 1.0
  FP <- sum(tU == 0 & eU == 1, na.rm = TRUE) * 1.0
  FN <- sum(tU == 1 & eU == 0, na.rm = TRUE) * 1.0
  TN <- sum(tU == 0 & eU == 0, na.rm = TRUE) * 1.0

  prec <- if ((TP + FP) > 0) TP / (TP + FP) else 0
  rec  <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  f1   <- if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0

  denom <- (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)
  mcc   <- if (is.finite(denom) && denom > 0) {
    (TP * TN - FP * FN) / sqrt(denom)
  } else 0

  c(F1 = f1, MCC = mcc, Precision = prec, Recall = rec,
    TP = TP, FP = FP, FN = FN, TN = TN)
}

.safe_prf <- function(TP, FP, FN) {
  prec <- if ((TP + FP) > 0) TP / (TP + FP) else 0
  rec  <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  f1   <- if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0
  c(precision = prec, recall = rec, f1 = f1)
}
