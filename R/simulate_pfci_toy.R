#' Simulate toy data for PFCI using topo-ordered DAG + rmvDAG
#'
#' Workflow:
#'   sim <- simulate_pfci_toy(...)
#'   fit <- pfci_fun(sim$X, ...)
#'   met <- pfci_metrics(sim, fit)
#'
#' This simulator:
#'  - generates a *topologically ordered* DAG (edges only i -> j for i < j)
#'  - simulates data via pcalg::rmvDAG with requested errDist
#'  - returns truth skeleton (undirected) and an "amat-style" truth from dag2cpdag
#'
#' NOTE: The returned truth_amat is derived from the CPDAG of the generating DAG
#' (so it contains directed and o-o circle edges, but not latent-induced o-> / <->).
#'
#' Backward-compat: accepts old args p_obs/gamma (ignored) so old vignettes won't fail.
#'
#' @param p Number of observed variables (preferred).
#' @param sparsity Number of nodes eligible for edges (<= p). Default p.
#' @param n Sample size.
#' @param edge_prob Edge probability among eligible nodes.
#' @param errDist Error distribution for pcalg::rmvDAG ("normal","t4","mixt3").
#' @param seed Random seed.
#' @param p_obs (legacy) alias for p.
#' @param gamma (legacy) ignored (kept only for backward compatibility).
#'
#' @return A list: X, truth (true_dag, adj_mat, skel, amat), meta
#' @seealso \code{\link{pfci_fit}}, \code{\link{pfci_metrics}}
#' @examples
#' sim <- simulate_pfci_toy(p = 30, n = 100, edge_prob = 0.05, seed = 1)
#' str(sim$truth)
#' @export
simulate_pfci_toy <- function(p = NULL,
                              sparsity = NULL,
                              n = 100,
                              edge_prob = 0.02,
                              errDist = c("normal", "t4", "mixt3"),
                              seed = 1L,
                              p_obs = NULL,
                              gamma = 0.1) {

  errDist <- match.arg(errDist)

  if (!requireNamespace("pcalg", quietly = TRUE)) {
    stop(
      "Package 'pcalg' is required for simulate_pfci_toy() but is not installed.\n",
      "Please install 'pcalg' from Bioconductor.",
      call. = FALSE
    )
  }

  if (!requireNamespace("graph", quietly = TRUE)) {
    stop(
      "Package 'graph' is required for simulate_pfci_toy() but is not installed.\n",
      "Please install 'graph' from Bioconductor.",
      call. = FALSE
    )
  }

  # Backward compat: if p not provided, use p_obs
  if (is.null(p)) {
    if (!is.null(p_obs)) {
      p <- as.integer(p_obs)
    } else {
      p <- 100L
    }
  }
  p <- as.integer(p)

  if (is.null(sparsity)) sparsity <- p
  sparsity <- as.integer(sparsity)

  set.seed(seed)

  dag_obj <- .generate_general_dag_topo(
    n_nodes   = p,
    sparsity  = sparsity,
    edge_prob = edge_prob
  )

  # rmvDAG requires topo-sorted graphNEL; guaranteed by construction
  X <- pcalg::rmvDAG(n = n, dag_obj$true_dag, errDist = errDist)
  colnames(X) <- paste0("X", seq_len(ncol(X)))

  truth_skel <- .dag_adj_to_skeleton(dag_obj$adj_mat)
  truth_amat <- .true_amat_from_dag(dag_obj$true_dag)

  list(
    X = X,
    truth = list(
      true_dag = dag_obj$true_dag,
      adj_mat  = dag_obj$adj_mat,
      skel     = truth_skel,
      amat     = truth_amat
    ),
    meta = list(
      p = p,
      sparsity = sparsity,
      n = n,
      edge_prob = edge_prob,
      errDist = errDist,
      seed = seed,
      # legacy args kept for record only
      p_obs = p_obs,
      gamma = gamma
    )
  )
}

# -------------------------------------------------------------------------
# INTERNAL helpers (NOT exported)
# -------------------------------------------------------------------------

.generate_general_dag_topo <- function(n_nodes = 100L, sparsity = n_nodes, edge_prob = 0.02) {
  n_nodes <- as.integer(n_nodes)
  sparsity <- as.integer(min(sparsity, n_nodes))

  A <- matrix(0L, n_nodes, n_nodes)
  s_nodes <- sort(sample.int(n_nodes, size = sparsity, replace = FALSE))

  # only edges i -> j with i < j (topological order is node index)
  for (ii in seq_along(s_nodes)) {
    i <- s_nodes[ii]
    if (ii == length(s_nodes)) next
    for (jj in (ii + 1):length(s_nodes)) {
      j <- s_nodes[jj]
      if (runif(1) < edge_prob) A[i, j] <- 1L
    }
  }

  # Build graphNEL *without* relying on as(matrix,"graphNEL") coercion
  node_names <- as.character(seq_len(n_nodes))
  g <- graph::graphNEL(nodes = node_names, edgemode = "directed")

  # graph::addEdge is exported; addEdges is NOT (this was breaking your check)
  for (i in seq_len(n_nodes)) {
    kids <- which(A[i, ] == 1L)
    if (length(kids)) {
      for (j in kids) {
        g <- graph::addEdge(from = node_names[i], to = node_names[j], graph = g)
      }
    }
  }

  list(true_dag = g, adj_mat = A)
}

.dag_adj_to_skeleton <- function(adj_mat) {
  A <- (adj_mat != 0L) * 1L
  A <- ((A + t(A)) > 0L) * 1L
  diag(A) <- 0L
  A
}

# Convert CPDAG to pcalg-style amat codes:
# 0 none, 1 circle, 2 arrow, 3 tail
.cpdag_to_amat <- function(cpdag_graphNEL) {
  labs <- graph::nodes(cpdag_graphNEL)
  p <- length(labs)
  amat <- matrix(0L, p, p, dimnames = list(labs, labs))

  M <- as(cpdag_graphNEL, "matrix")
  M[M != 0] <- 1

  for (i in 1:(p - 1)) for (j in (i + 1):p) {
    ij <- (M[i, j] == 1)
    ji <- (M[j, i] == 1)
    if (!ij && !ji) next

    if (ij && ji) {
      # undirected edge => o-o (circle endpoints)
      amat[i, j] <- 1L
      amat[j, i] <- 1L
    } else if (ij && !ji) {
      # i -> j
      amat[i, j] <- 2L
      amat[j, i] <- 3L
    } else {
      # j -> i
      amat[j, i] <- 2L
      amat[i, j] <- 3L
    }
  }

  diag(amat) <- 0L
  amat
}

.true_amat_from_dag <- function(dag_graphNEL) {
  cpdag <- pcalg::dag2cpdag(dag_graphNEL)
  .cpdag_to_amat(cpdag)
}
