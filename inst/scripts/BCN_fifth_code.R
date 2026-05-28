install.packages("BiocManager")
BiocManager::install("Rgraphviz")
library(graph)
library(glasso)
library(pcalg)
library(igraph)

generate_general_dag <- function(n_nodes = 100, 
                                 edge_prob = 0.2,
                                 sparsity = 100) {
  
  # Create an empty adjacency matrix
  dag_matrix <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  
  # # Select nodes that Y (node 1) will point to
  # target_nodes <- sample(2:n_nodes, n_outgoing_from_Y, replace = FALSE)
  # dag_matrix[1, target_nodes] <- 1  # Y points to selected nodes
  
  # Sample a subset of nodes for sparse connections
  s_nodes <- sort(sample(n_nodes, sparsity, replace = FALSE))
  
  # Generate a sparse DAG structure among these selected nodes
  for (i in seq_along(s_nodes)) {
    for (j in (i + 1):length(s_nodes)) {
      if (runif(1) < edge_prob) {
        dag_matrix[s_nodes[i], s_nodes[j]] <- 1  # Directed edge from i to j
      }
    }
  }
  
  # Generate a sparse DAG structure among the rest of the nodes
  # for (i in 2:n_nodes) {
  #   for (j in (i+1):n_nodes) {
  #     if (runif(1) < edge_prob) {
  #       dag_matrix[i, j] <- 1  # Directed edge from i to j
  #     }
  #   }
  # }
  
  # Convert adjacency matrix to an DAG object
  true_dag <- as(dag_matrix, "graphNEL")
  
  return(list(true_dag = true_dag, adj_mat = dag_matrix))
}

# Example usage
set.seed(864898)
dag <- generate_general_dag(n_nodes = 700, edge_prob = 0.015, sparsity = 500)

# dag <- generate_general_dag(n_nodes = 700, n_outgoing_from_Y = 50,
#                             edge_prob = 0.005, sparsity = 500)
true_dag <- dag$true_dag
adj_mat <- dag$adj_mat
plot(true_dag)
sum(adj_mat)

run_simulation <- function() {
  # Define a known DAG with p groups of nodes
  p <- 500
  s <- 500
  true_dag <- generate_general_dag(n_nodes = p, sparsity = s, edge_prob = 0.005)  # Adjust group number as needed
  
  # Introduce Multicollinearity and generate the data
  n <- 100
  data <- rmvDAG(n, true_dag$true_dag, errDist = "mixt3")  # Generate data from the DAG
  
  # Compute the empirical covariance matrix
  S <- cov(data)
  
  # Apply graphical lasso
  # Tune lambda
  # lambdas <- seq(0, 0.5, by = 0.01)
  deg <- 1
  eps <- log(deg)/log(n) + 0.01
  best_lambda <- n^(-(1-eps)/2)
  start_pfci <- Sys.time()
  fit_glasso <- glasso(S, rho = best_lambda, approx = TRUE)
  end_pfci <- Sys.time()
  
  # Extract the adjacency matrix (skeleton)
  adj_glasso <- (fit_glasso$wi != 0)*1 # Non-zero precision matrix entries
  diag(adj_glasso) <- 0  # Remove diagonal
  
  adj_glasso_sym = (adj_glasso + t(adj_glasso))/2
  adj_complement = !adj_glasso_sym
  
  # Use the skeleton from graphical lasso as input
  suffStat <- list(C = cor(data), n = nrow(data))  # Sufficient statistics
  
  custom_gaussCItest <- function(x, y, S, suffStat) {
    if (adj_glasso[x, y] == 0) {
      return(1)  # Treat this edge as independent if glasso says no edge
    }
    return(gaussCItest(x, y, S, suffStat))  # Otherwise, perform the test
  }
  
  time1 = end_pfci - start_pfci
  # FCI
  start_pfci2 = Sys.time()
  fci_res <- fci(suffStat, indepTest = custom_gaussCItest, 
                 alpha = 0.05, labels = colnames(data), 
                 fixedGaps = adj_complement, skel.method = "stable",
                 doPdsep = FALSE)
  end_pfci2 = Sys.time()
  
  time2 = end_pfci2 - start_pfci2
  time_pfci = time1 + time2
  
  # RFCI
  start_rfci = Sys.time()
  rfci_res <- rfci(suffStat, indepTest = gaussCItest, alpha = 0.05, 
                   labels = colnames(data), skel.method = "stable")
  end_rfci = Sys.time()
  
  time_rfci = end_rfci - start_rfci
  
  # Compare the true DAG with the FCI result
  # oracle_dag <- generate_oracle_dag(true_dag)
  oracle_dag <- true_dag$true_dag
  true_skeleton <- as(oracle_dag, "matrix")
  
  # Convert any non-zero values in true_skeleton to 1 (making it binary)
  true_skeleton[true_skeleton != 0] <- 1
  
  # Ensure the skeleton is symmetric (since skeletons are undirected)
  true_skeleton <- (true_skeleton | t(true_skeleton)) * 1
  
  # Extract estimated skeletons
  estimated_skeleton_fci <- as(fci_res@amat, "matrix")
  estimated_skeleton_rfci <- as(rfci_res@amat, "matrix")
  
  # Convert any non-zero values in estimated_skeleton to 1
  estimated_skeleton_fci[estimated_skeleton_fci != 0] <- 1
  estimated_skeleton_rfci[estimated_skeleton_rfci != 0] <- 1
  
  # Ensure the skeleton is symmetric
  estimated_skeleton_fci <- (estimated_skeleton_fci | t(estimated_skeleton_fci)) * 1
  estimated_skeleton_rfci <- (estimated_skeleton_rfci | t(estimated_skeleton_rfci)) * 1
  
  # Convert your true_skeleton and estimated_skeleton to graphNEL format
  true_dag_graphNEL <- as(as(true_dag$adj_mat, "matrix"), "graphNEL")  # Assuming true_dag is already in a suitable format
  estimated_dag_fci <- as(as(fci_res@amat, "matrix"), "graphNEL")  # Use your estimated DAG
  estimated_dag_rfci <- as(as(rfci_res@amat, "matrix"), "graphNEL")  # Use your estimated DAG
  
  # Compute SHD between true DAG and estimated DAG
  shd_fci <- shd(true_dag_graphNEL, estimated_dag_fci)
  shd_rfci <- shd(true_dag_graphNEL, estimated_dag_rfci)
  
  # Calculate True Positives, False Positives, False Negatives, True Negatives for FCI
  TP_fci <- sum((true_skeleton == 1) & (estimated_skeleton_fci == 1))
  FP_fci <- sum((true_skeleton == 0) & (estimated_skeleton_fci == 1))
  FN_fci <- sum((true_skeleton == 1) & (estimated_skeleton_fci == 0))
  TN_fci <- sum((true_skeleton == 0) & (estimated_skeleton_fci == 0))
  
  # MCC calculation, converting values to numeric to avoid integer overflow
  MCC_fci <- (as.numeric(TP_fci) * as.numeric(TN_fci) - as.numeric(FP_fci) * as.numeric(FN_fci)) / 
    sqrt((as.numeric(TP_fci + FP_fci) * as.numeric(TP_fci + FN_fci) * as.numeric(TN_fci + FP_fci) * as.numeric(TN_fci + FN_fci)))
  
  # Calculate metrics for FCI
  precision_fci <- TP_fci / (TP_fci + FP_fci)
  recall_fci <- TP_fci / (TP_fci + FN_fci)
  
  # Calculate True Positives, False Positives, False Negatives, True Negatives for RFCI
  TP_rfci <- sum((true_skeleton == 1) & (estimated_skeleton_rfci == 1))
  FP_rfci <- sum((true_skeleton == 0) & (estimated_skeleton_rfci == 1))
  FN_rfci <- sum((true_skeleton == 1) & (estimated_skeleton_rfci == 0))
  TN_rfci <- sum((true_skeleton == 0) & (estimated_skeleton_rfci == 0))
  
  # Calculate MCC for RFCI
  MCC_rfci <- (as.numeric(TP_rfci) * as.numeric(TN_rfci) - as.numeric(FP_rfci) * as.numeric(FN_rfci)) / 
    sqrt((as.numeric(TP_rfci + FP_rfci) * as.numeric(TP_rfci + FN_rfci) * as.numeric(TN_rfci + FP_rfci) * as.numeric(TN_rfci + FN_rfci)))
  
  
  # Calculate metrics for RFCI
  precision_rfci <- TP_rfci / (TP_rfci + FP_rfci)
  recall_rfci <- TP_rfci / (TP_rfci + FN_rfci)
  
  F1_pfci = (2*precision_fci*recall_fci)/(precision_fci + recall_fci)
  F1_rfci = (2*precision_rfci*recall_rfci)/(precision_rfci + recall_rfci)
  
  # return(c(precision_fci = precision_fci, recall_fci = recall_fci,
  #          precision_rfci = precision_rfci, recall_rfci = recall_rfci,
  #          shd_fci = shd_fci, shd_rfci = shd_rfci,
  #          TP_fci = TP_fci, FP_fci = FP_fci, FN_fci = FN_fci,
  #          TP_rfci = TP_rfci, FP_rfci = FP_rfci, FN_rfci = FN_rfci))
  
  return(c(pfci_shd = shd_fci, rfci_shd = shd_rfci,
           pfci_f1 = F1_pfci, rfci_f1 = F1_rfci,
           pfci_time = time_pfci, rfci_time = time_rfci,
           pfci_mcc = MCC_fci, rfci_mcc = MCC_rfci))
}


num_replicates <- 20
results <- replicate(num_replicates, run_simulation())
rowMeans(results)


library(graph)
library(glasso)
library(pcalg)
library(igraph)

################################################################################
############ Non-normal                                            #############
################################################################################

run_simulation <- function() {
  # Define a known DAG with p groups of nodes
  p <- 700
  s <- 700
  true_dag <- generate_general_dag(n_nodes = p, sparsity = s, edge_prob = 0.005)  # Adjust group number as needed
  
  # Introduce Multicollinearity and generate the data
  n <- 100
  data <- rmvDAG(n, true_dag$true_dag, errDist = "t4")  # Generate data from the DAG
  
  # Compute the empirical covariance matrix
  S <- cov(data)
  
  # Apply graphical lasso
  # Tune lambda
  # lambdas <- seq(0, 0.5, by = 0.01)
  deg <- 1
  eps <- log(deg)/log(n) + 0.01
  best_lambda <- n^(-(1-eps)/2)
  start_pfci <- Sys.time()
  fit_glasso <- glasso(S, rho = best_lambda, approx = TRUE)
  end_pfci <- Sys.time()
  
  # Extract the adjacency matrix (skeleton)
  adj_glasso <- (fit_glasso$wi != 0)*1 # Non-zero precision matrix entries
  diag(adj_glasso) <- 0  # Remove diagonal
  
  adj_glasso_sym = (adj_glasso + t(adj_glasso))/2
  adj_complement = !adj_glasso_sym
  
  # Use the skeleton from graphical lasso as input
  suffStat <- list(C = cor(data), n = nrow(data))  # Sufficient statistics
  
  custom_gaussCItest <- function(x, y, S, suffStat) {
    if (adj_glasso[x, y] == 0) {
      return(1)  # Treat this edge as independent if glasso says no edge
    }
    return(gaussCItest(x, y, S, suffStat))  # Otherwise, perform the test
  }
  
  time1 = end_pfci - start_pfci
  # FCI
  start_pfci2 = Sys.time()
  fci_res <- fci(suffStat, indepTest = custom_gaussCItest, 
                 alpha = 0.05, labels = colnames(data), 
                 fixedGaps = adj_complement, skel.method = "stable",
                 doPdsep = FALSE)
  end_pfci2 = Sys.time()
  
  time2 = end_pfci2 - start_pfci2
  time_pfci = time1 + time2
  
  # RFCI
  start_rfci = Sys.time()
  rfci_res <- rfci(suffStat, indepTest = gaussCItest, alpha = 0.05, 
                   labels = colnames(data), skel.method = "stable")
  end_rfci = Sys.time()
  
  time_rfci = end_rfci - start_rfci
  
  # Compare the true DAG with the FCI result
  # oracle_dag <- generate_oracle_dag(true_dag)
  oracle_dag <- true_dag$true_dag
  true_skeleton <- as(oracle_dag, "matrix")
  
  # Convert any non-zero values in true_skeleton to 1 (making it binary)
  true_skeleton[true_skeleton != 0] <- 1
  
  # Ensure the skeleton is symmetric (since skeletons are undirected)
  true_skeleton <- (true_skeleton | t(true_skeleton)) * 1
  
  # Extract estimated skeletons
  estimated_skeleton_fci <- as(fci_res@amat, "matrix")
  estimated_skeleton_rfci <- as(rfci_res@amat, "matrix")
  
  # Convert any non-zero values in estimated_skeleton to 1
  estimated_skeleton_fci[estimated_skeleton_fci != 0] <- 1
  estimated_skeleton_rfci[estimated_skeleton_rfci != 0] <- 1
  
  # Ensure the skeleton is symmetric
  estimated_skeleton_fci <- (estimated_skeleton_fci | t(estimated_skeleton_fci)) * 1
  estimated_skeleton_rfci <- (estimated_skeleton_rfci | t(estimated_skeleton_rfci)) * 1
  
  # Convert your true_skeleton and estimated_skeleton to graphNEL format
  true_dag_graphNEL <- as(as(true_dag$adj_mat, "matrix"), "graphNEL")  # Assuming true_dag is already in a suitable format
  estimated_dag_fci <- as(as(fci_res@amat, "matrix"), "graphNEL")  # Use your estimated DAG
  estimated_dag_rfci <- as(as(rfci_res@amat, "matrix"), "graphNEL")  # Use your estimated DAG
  
  # Compute SHD between true DAG and estimated DAG
  shd_fci <- shd(true_dag_graphNEL, estimated_dag_fci)
  shd_rfci <- shd(true_dag_graphNEL, estimated_dag_rfci)
  
  # Calculate True Positives, False Positives, False Negatives, True Negatives for FCI
  TP_fci <- sum((true_skeleton == 1) & (estimated_skeleton_fci == 1))
  FP_fci <- sum((true_skeleton == 0) & (estimated_skeleton_fci == 1))
  FN_fci <- sum((true_skeleton == 1) & (estimated_skeleton_fci == 0))
  TN_fci <- sum((true_skeleton == 0) & (estimated_skeleton_fci == 0))
  
  # MCC calculation, converting values to numeric to avoid integer overflow
  MCC_fci <- (as.numeric(TP_fci) * as.numeric(TN_fci) - as.numeric(FP_fci) * as.numeric(FN_fci)) / 
    sqrt((as.numeric(TP_fci + FP_fci) * as.numeric(TP_fci + FN_fci) * as.numeric(TN_fci + FP_fci) * as.numeric(TN_fci + FN_fci)))
  
  # Calculate metrics for FCI
  precision_fci <- TP_fci / (TP_fci + FP_fci)
  recall_fci <- TP_fci / (TP_fci + FN_fci)
  
  # Calculate True Positives, False Positives, False Negatives, True Negatives for RFCI
  TP_rfci <- sum((true_skeleton == 1) & (estimated_skeleton_rfci == 1))
  FP_rfci <- sum((true_skeleton == 0) & (estimated_skeleton_rfci == 1))
  FN_rfci <- sum((true_skeleton == 1) & (estimated_skeleton_rfci == 0))
  TN_rfci <- sum((true_skeleton == 0) & (estimated_skeleton_rfci == 0))
  
  # Calculate MCC for RFCI
  MCC_rfci <- (as.numeric(TP_rfci) * as.numeric(TN_rfci) - as.numeric(FP_rfci) * as.numeric(FN_rfci)) / 
    sqrt((as.numeric(TP_rfci + FP_rfci) * as.numeric(TP_rfci + FN_rfci) * as.numeric(TN_rfci + FP_rfci) * as.numeric(TN_rfci + FN_rfci)))
  
  
  # Calculate metrics for RFCI
  precision_rfci <- TP_rfci / (TP_rfci + FP_rfci)
  recall_rfci <- TP_rfci / (TP_rfci + FN_rfci)
  
  F1_pfci = (2*precision_fci*recall_fci)/(precision_fci + recall_fci)
  F1_rfci = (2*precision_rfci*recall_rfci)/(precision_rfci + recall_rfci)
  
  # return(c(precision_fci = precision_fci, recall_fci = recall_fci,
  #          precision_rfci = precision_rfci, recall_rfci = recall_rfci,
  #          shd_fci = shd_fci, shd_rfci = shd_rfci,
  #          TP_fci = TP_fci, FP_fci = FP_fci, FN_fci = FN_fci,
  #          TP_rfci = TP_rfci, FP_rfci = FP_rfci, FN_rfci = FN_rfci))
  
  return(c(pfci_shd = shd_fci, rfci_shd = shd_rfci,
           pfci_f1 = F1_pfci, rfci_f1 = F1_rfci,
           pfci_time = time_pfci, rfci_time = time_rfci,
           pfci_mcc = MCC_fci, rfci_mcc = MCC_rfci))
}

num_replicates <- 20
results <- replicate(num_replicates, run_simulation())
rowMeans(results)
