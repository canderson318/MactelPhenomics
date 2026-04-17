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
n_null_matrices <- 1e4    # number of semi-random null matrices for null distribution
null_seed       <- 123    # reproducibility seed for null distribution

# =============================================================================
# STEP 1: Observed transition matrix (symmetric, diagonal = 0)
# =============================================================================

obs_mat <- t(transition_counts)   # rows = source, cols = destination
diag(obs_mat) <- 0L
num_clusts <- nrow(obs_mat)
clust_ids  <- as.character(seq_len(num_clusts))
dimnames(obs_mat) <- list(clust_ids, clust_ids)

# =============================================================================
# STEP 2: Helpers
# =============================================================================

# G-statistic: 2 * sum(O * log(O/E)), off-diagonal only, ignoring zero cells
compute_g_stat <- function(O, E) {
  valid <- (O > 0 & E > 0) & (row(O) != col(O))
  2 * sum(O[valid] * log(O[valid] / E[valid]))
}

# Convert a graph distance matrix to a row-normalised proximity-based expected
# matrix scaled to the observed row sums.
# Proximity = 1 / (1 + distance); diagonal forced to 0; then row-normalised.
make_expected <- function(graph, observed) {
  dist_mat <- as.matrix(distances(graph, mode = "all"))
  prox     <- 1 / (1 + dist_mat)
  diag(prox) <- 0
  prob     <- sweep(prox, 1, rowSums(prox), "/")   # row-normalise
  sweep(prob, 1, rowSums(observed), "*")            # scale to observed row sums
}

# =============================================================================
# STEP 3: Bifurcated MST graph
# =============================================================================

bif_g <- graph_from_data_frame(mst_edges, directed = FALSE,
                                vertices = data.frame(name = clust_ids))
g_bif <- compute_g_stat(obs_mat, make_expected(bif_g, obs_mat))

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
g_pp  <- compute_g_stat(obs_mat, make_expected(pp_g, obs_mat))

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
g_chew <- compute_g_stat(obs_mat, make_expected(chew_g, obs_mat))

# =============================================================================
# STEP 6: Fully unordered (complete) graph
# =============================================================================

unord_edges <- as.data.frame(t(combn(clust_ids, 2)))
names(unord_edges) <- c("from", "to")
unord_g <- graph_from_data_frame(unord_edges, directed = FALSE,
                                  vertices = data.frame(name = clust_ids))
g_unord <- compute_g_stat(obs_mat, make_expected(unord_g, obs_mat))

# =============================================================================
# STEP 7: Null distribution — semi-random matrices preserving row sums
# =============================================================================

set.seed(null_seed)
null_g_stats <- replicate(n_null_matrices, {
  # build random matrix preserving each row's sum via rmultinom with uniform probs
  rand_mat <- matrix(0L, nrow = num_clusts, ncol = num_clusts)
  for (i in seq_len(num_clusts)) {
    rs <- sum(obs_mat[i, ])
    if (rs == 0) next
    off_cols <- setdiff(seq_len(num_clusts), i)
    probs    <- runif(length(off_cols))
    probs    <- probs / sum(probs)
    rand_mat[i, off_cols] <- rmultinom(1, size = rs, prob = probs)
  }
  dimnames(rand_mat) <- list(clust_ids, clust_ids)
  # use random matrix as expected; compare against observed
  compute_g_stat(obs_mat, rand_mat)
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
