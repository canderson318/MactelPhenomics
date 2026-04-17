# =============================================================================
# 04_network_flow_gtest.R
#
# This script tests whether the observed longitudinal cluster transitions are
# better explained by the inferred bifurcated MST trajectory than by several
# alternative graph structures. For each candidate graph (bifurcated MST,
# linear-by-pseudoprogression, linear-by-Chew-grade, fully unordered) an
# expected transition matrix is derived and compared to the observed matrix
# using the G-statistic (log-likelihood ratio test). A null distribution is
# built from semi-random matrices that preserve observed row sums.
#
# Pipeline position: runs after clustering (02_clustering_trajectory.R).
# Requires transition_counts from that script.
# =============================================================================

# Required packages ------------------------------------------------------------
library(igraph)   # graph construction and distance matrices
library(dplyr)    # data manipulation

# Expected input ---------------------------------------------------------------
# transition_counts : integer matrix (num_clusts x num_clusts) of observed
#                     cluster-to-cluster transition counts.
#                     Rows = destination cluster, cols = source cluster.
#                     Output of 02_clustering_trajectory.R.
#
# mst_edges         : data.frame with columns "from" and "to" (character/integer),
#                     listing edges of the MST from 02_clustering_trajectory.R.
#
# cluster_pp_means  : data.frame with columns "my_cluster" and "mean_pseudoprog",
#                     cluster mean pseudoprogression scores.
#
# cluster_chew_means: data.frame with columns "my_cluster" and "mean_chew",
#                     cluster mean Chew grade.

# Parameters -------------------------------------------------------------------
n_null_matrices <- 500    # number of semi-random null matrices for null distribution
null_seed       <- 4219   # reproducibility seed for null distribution

# =============================================================================
# STEP 1: Observed transition matrix (symmetric, diagonal = 0)
# =============================================================================

obs_mat <- t(transition_counts)   # rows = source, cols = destination
diag(obs_mat) <- 0L
num_clusts <- nrow(obs_mat)
clust_ids  <- as.character(seq_len(num_clusts))
dimnames(obs_mat) <- list(clust_ids, clust_ids)

# =============================================================================
# STEP 2: Helper — G-statistic
# Compares observed counts to expected counts derived from a graph structure.
# Expected counts are proportional to graph adjacency (connected = positive weight,
# disconnected = 0), scaled to match observed row sums.
# =============================================================================

g_stat <- function(observed, graph) {
  dist_mat   <- as.matrix(distances(graph, mode = "all"))
  adj_mat    <- (dist_mat == 1) * 1.0            # adjacency: 1 if directly connected
  diag(adj_mat) <- 0

  # for each source cluster (row), distribute row sum according to adjacency
  expected <- sweep(adj_mat, 1, rowSums(adj_mat), "/")   # row-normalise
  expected <- sweep(expected, 1, rowSums(observed), "*")  # scale to row counts

  # G = 2 * sum(O * log(O/E)), ignoring zeros
  mask <- observed > 0 & expected > 0
  2 * sum(observed[mask] * log(observed[mask] / expected[mask]))
}

# =============================================================================
# STEP 3: Bifurcated MST graph
# =============================================================================

bif_g <- graph_from_data_frame(mst_edges, directed = FALSE,
                                vertices = data.frame(name = clust_ids))
g_bif <- g_stat(obs_mat, bif_g)

# =============================================================================
# STEP 4: Linear graph ordered by pseudoprogression mean
# =============================================================================

pp_order <- cluster_pp_means$my_cluster[order(cluster_pp_means$mean_pseudoprog)]
pp_edges <- data.frame(
  from = as.character(pp_order[-length(pp_order)]),
  to   = as.character(pp_order[-1])
)
pp_g  <- graph_from_data_frame(pp_edges, directed = FALSE,
                                vertices = data.frame(name = clust_ids))
g_pp  <- g_stat(obs_mat, pp_g)

# =============================================================================
# STEP 5: Linear graph ordered by Chew grade mean
# =============================================================================

chew_order <- cluster_chew_means$my_cluster[order(cluster_chew_means$mean_chew)]
chew_edges <- data.frame(
  from = as.character(chew_order[-length(chew_order)]),
  to   = as.character(chew_order[-1])
)
chew_g <- graph_from_data_frame(chew_edges, directed = FALSE,
                                 vertices = data.frame(name = clust_ids))
g_chew <- g_stat(obs_mat, chew_g)

# =============================================================================
# STEP 6: Fully unordered (complete) graph
# =============================================================================

unord_edges <- as.data.frame(t(combn(clust_ids, 2)))
names(unord_edges) <- c("from", "to")
unord_g <- graph_from_data_frame(unord_edges, directed = FALSE,
                                  vertices = data.frame(name = clust_ids))
g_unord <- g_stat(obs_mat, unord_g)

# =============================================================================
# STEP 7: Null distribution — semi-random matrices preserving row sums
# =============================================================================

row_sums_obs <- rowSums(obs_mat)

set.seed(null_seed)
null_g_stats <- replicate(n_null_matrices, {
  # sample a random matrix row by row, preserving each row's sum
  rand_mat <- t(sapply(row_sums_obs, function(rs) {
    x <- rep(0, num_clusts)
    if (rs > 0) {
      idx <- sample(seq_len(num_clusts), min(rs, num_clusts), replace = rs > num_clusts)
      x[idx] <- x[idx] + 1
    }
    x
  }))
  dimnames(rand_mat) <- list(clust_ids, clust_ids)
  diag(rand_mat) <- 0
  g_stat(rand_mat, bif_g)   # compare random matrix against bifurcated structure
})

# =============================================================================
# STEP 8: Summary
# =============================================================================

results <- data.frame(
  graph       = c("bifurcated_MST", "linear_pseudoprog", "linear_chew", "unordered"),
  G_statistic = c(g_bif, g_pp, g_chew, g_unord)
)

null_quantiles <- quantile(null_g_stats, c(0.025, 0.5, 0.975))

cat("G-statistics for each candidate graph structure:\n")
print(results)
cat("\nNull distribution quantiles (bifurcated graph, random row-sum-matched matrices):\n")
print(null_quantiles)

# =============================================================================
# OUTPUT
# =============================================================================
# results        : data.frame of G-statistics per candidate graph
# null_g_stats   : numeric vector of G-statistics under the null
# null_quantiles : 2.5 / 50 / 97.5 percentiles of the null distribution
