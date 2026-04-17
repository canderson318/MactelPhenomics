# =============================================================================
# 02_clustering_trajectory.R
#
# This script takes the imputed feature matrix from 01_imputation.R (stored in
# a SingleCellExperiment object) and performs dimensionality reduction,
# two-step clustering, cluster merging, and pseudotime inference via minimum
# spanning tree (MST). The result is a cluster label and continuous
# pseudoprogression score for each sample. Cluster-to-cluster transition
# frequencies (from longitudinal data) are also calculated.
#
# Pipeline position: runs after imputation (01_imputation.R), before XGBoost
# profiling (03_xgboost_shap.R) and network validation (04_network_flow_gtest.R).
# =============================================================================

# Required packages ------------------------------------------------------------
library(SingleCellExperiment)  # data container
library(scater)                # runPCA, runTSNE, runUMAP
library(bluster)               # clusterCells, TwoStepParam, KmeansParam, NNGraphParam
library(TSCAN)                 # createClusterMST, mapCellsToEdges, orderCells, averagePseudotime
library(scran)                 # aggregateAcrossCells

# Expected input ---------------------------------------------------------------
# sce : SingleCellExperiment object where:
#   - assay "pmm_imp_of_scaled" contains the imputed scaled feature matrix
#     (features x samples), output of 01_imputation.R
#   - colData contains at minimum a sample ID column

# Parameters -------------------------------------------------------------------
pca_ncomponents <- 50      # number of PCA components

# Two-step clustering params (seed + KmeansParam + NNGraphParam)
clust_seed    <- 1010
clust_centers <- 1440      # number of K-means micro-clusters
clust_k       <- 10        # k for nearest-neighbour graph step

# Cluster merges: similar clusters identified by visual/quantitative inspection
# Format: named vector where name → merges into value
cluster_merges <- c("4" = "7", "9" = "2", "12" = "3")

# TSNE params
tsne_params <- list(
  dimred              = "PCA",
  perplexity          = 50,
  final_momentum      = 0.9,
  momentum            = 0.5,
  eta                 = 100,
  exaggeration_factor = 4,
  seed                = 29209
)

# UMAP params
umap_params <- list(
  dimred      = "PCA",
  metric      = "euclidean",
  spread      = 1.2,
  min_dist    = 0.02,
  n_trees     = 50,
  n_neighbors = 20,
  search_k    = 2000,
  seed        = 29209
)

# =============================================================================
# STEP 1: PCA
# =============================================================================

set.seed(umap_params$seed)
sce <- runPCA(sce, exprs_values = "pmm_imp_of_scaled", ncomponents = pca_ncomponents)

# =============================================================================
# STEP 2: TSNE and UMAP on PCA
# =============================================================================

tsne_par <- c(tsne_params[names(tsne_params) != "seed"], list(x = sce))
set.seed(tsne_params$seed)
sce <- do.call(runTSNE, tsne_par)

umap_par <- c(umap_params[names(umap_params) != "seed"], list(x = sce))
set.seed(umap_params$seed)
sce <- do.call(runUMAP, umap_par)

# =============================================================================
# STEP 3: Two-step clustering (K-means micro-clusters → NNGraph)
# =============================================================================

set.seed(clust_seed)
km_nn <- clusterCells(
  sce,
  use.dimred = "PCA",
  BLUSPARAM  = TwoStepParam(
    first  = KmeansParam(centers = clust_centers),
    second = NNGraphParam(k = clust_k)
  ),
  full = TRUE
)

sce$km_nn_cluster <- km_nn$clusters

# =============================================================================
# STEP 4: Merge similar clusters
# =============================================================================

sce$my_cluster <- as.character(sce$km_nn_cluster)

for (from_clust in names(cluster_merges)) {
  sce$my_cluster[sce$my_cluster == from_clust] <- cluster_merges[[from_clust]]
}

# Re-factor sequentially
sce$my_cluster <- factor(
  as.numeric(factor(sce$my_cluster, levels = as.character(sort(unique(as.numeric(sce$my_cluster))))))
)

# =============================================================================
# STEP 5: MST on cluster centroids in PCA space
# =============================================================================

colLabels(sce) <- sce$my_cluster

# aggregate per cluster to get centroids
by_cluster <- aggregateAcrossCells(sce,
                                   use.assay.type = "pmm_imp_of_scaled",
                                   ids            = sce$my_cluster,
                                   statistics     = "mean")

centroids_pca <- reducedDim(by_cluster, "PCA")

# build MST
mst_pca <- TSCAN::createClusterMST(centroids_pca, clusters = NULL)

# report MST edges in PCA and TSNE space for plotting
pca_line_data  <- reportEdges(by_cluster, mst = mst_pca, clusters = NULL, use.dimred = "PCA")
tsne_line_data <- reportEdges(by_cluster, mst = mst_pca, clusters = NULL, use.dimred = "TSNE")

# =============================================================================
# STEP 6: Pseudoprogression — project each sample onto MST path
# =============================================================================

map_tscan       <- mapCellsToEdges(sce, mst = mst_pca, use.dimred = "PCA")
tscan_pseudo    <- orderCells(map_tscan, mst_pca)
sce$pseudoprog  <- averagePseudotime(tscan_pseudo)

# =============================================================================
# STEP 7: Cluster-to-cluster transition frequencies (longitudinal data)
# Requires colData to contain: patID_eye (eye-level ID), month (visit time),
# my_cluster (cluster assignment per visit).
# =============================================================================

num_clusts <- nlevels(sce$my_cluster)
clust_levels <- levels(sce$my_cluster)

pat_d <- as.data.frame(colData(sce)) |>
  dplyr::arrange(patID_eye, month) |>
  dplyr::select(patID_eye, month, my_cluster) |>
  dplyr::mutate(my_cluster = as.numeric(as.character(my_cluster))) |>
  dplyr::group_by(patID_eye) |>
  dplyr::filter(dplyr::n_distinct(my_cluster) > 1) |>
  dplyr::ungroup()

# matrix: rows = destination cluster, cols = source cluster
transition_counts <- matrix(
  0L,
  nrow = num_clusts, ncol = num_clusts,
  dimnames = list(
    paste0("to.",   clust_levels),
    paste0("from.", clust_levels)
  )
)

for (i in seq_len(num_clusts)) {
  from_i <- pat_d |>
    dplyr::group_by(patID_eye) |>
    dplyr::arrange(month) |>
    dplyr::filter(dplyr::lag(my_cluster) == i) |>
    dplyr::ungroup() |>
    dplyr::pull(my_cluster) |>
    table()

  for (nm in names(from_i)) {
    transition_counts[paste0("to.", nm), paste0("from.", i)] <- from_i[[nm]]
  }
}

# =============================================================================
# OUTPUT
# =============================================================================
# sce              : SCE with reducedDims PCA/TSNE/UMAP, colData columns
#                    km_nn_cluster, my_cluster, pseudoprog
# mst_pca          : igraph MST object
# pca_line_data    : data.frame of MST edge coordinates in PCA space
# tsne_line_data   : data.frame of MST edge coordinates in TSNE space
# transition_counts: integer matrix (num_clusts x num_clusts) of observed transitions
