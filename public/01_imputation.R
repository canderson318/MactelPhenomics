# =============================================================================
# 01_imputation.R
#
# This script takes a pre-cleaned numeric feature matrix (samples x features)
# and performs the imputation and derived-variable steps used in the manuscript.
# Key steps: feature and sample filtering by missingness and skewness; grouping
# variables by their missing-data pattern; parallel MICE PMM imputation with
# quality control; post-imputation constraint fixes; and calculation of derived
# features (parafoveal average, foveal slope). Output is an imputed matrix
# ready for dimensionality reduction.
#
# Pipeline position: runs after data cleaning, before clustering (02_clustering_trajectory.R).
# =============================================================================

# Required packages ------------------------------------------------------------
library(mice)       # multiple imputation
library(psych)      # skewness
library(parallel)   # mclapply

# Expected input ---------------------------------------------------------------
# feature_matrix : numeric matrix, rows = samples, cols = features
#                  Row names = sample IDs. NAs allowed (this script handles them).
#                  Assumed to contain only continuous/ordinal measurement variables;
#                  no ID columns, no demographic columns.
#
# mt_zone_cols   : character vector of macular thickness column names (e.g. "MT01":"MT09")
#                  used to calculate foveal slope and parafoveal average.
#
# isos_cols      : character vector of length 3: c("ISOS_BREAK", "ISOS_LOC_JF", "ISOS_LOC_CEN")
#                  used to enforce post-imputation logical constraints.

# Parameters -------------------------------------------------------------------
skew_threshold    <- 30      # remove features with |skewness| > this
miss_var_cutoff   <- 0.20    # remove features missing > this fraction
miss_samp_cutoff  <- NULL    # if NULL, uses boxplot outlier method on sample missingness
n_miss_clusters   <- 5       # number of missingness-pattern clusters for grouped imputation
imp_method        <- "pmm"   # mice imputation method
imp_m             <- 11      # number of parallel imputations to average
imp_seed          <- 23939   # seed for reproducibility
zone_distance_um  <- 500     # inter-zone distance in microns (for slope calculation)
imp_qc_ratio_na   <- 0.10    # fraction to mask per missingness-group during QC test

# =============================================================================
# STEP 1: Filter features by missingness and skewness
# =============================================================================

miss_ratio_vars <- colSums(is.na(feature_matrix)) / nrow(feature_matrix)

# remove features above missingness threshold
out_vars_miss <- names(miss_ratio_vars)[miss_ratio_vars >= miss_var_cutoff]

# remove features with extreme skewness
skw <- apply(feature_matrix, 2, function(x) abs(psych::skew(x, na.rm = TRUE)))
out_vars_skew <- names(skw)[skw > skew_threshold]

out_vars <- unique(c(out_vars_miss, out_vars_skew))
feature_matrix <- feature_matrix[, !colnames(feature_matrix) %in% out_vars]

# =============================================================================
# STEP 2: Filter samples by missingness (boxplot outlier method)
# =============================================================================

miss_ratio_samps <- rowSums(is.na(feature_matrix)) / ncol(feature_matrix)

if (is.null(miss_samp_cutoff)) {
  miss_samp_cutoff <- min(boxplot(miss_ratio_samps, plot = FALSE)$out)
}

out_samples <- rownames(feature_matrix)[miss_ratio_samps >= miss_samp_cutoff]
feature_matrix <- feature_matrix[!rownames(feature_matrix) %in% out_samples, ]

# =============================================================================
# STEP 3: Scale (mean-center + standardize)
# =============================================================================

scaled_matrix <- scale(feature_matrix)
scale_center <- attr(scaled_matrix, "scaled:center")
scale_sd     <- attr(scaled_matrix, "scaled:scale")

# =============================================================================
# STEP 4: Cluster variables by missingness pattern
# (groups co-missing variables so imputation respects data-collection structure)
# =============================================================================

miss_binary <- is.na(scaled_matrix) * 1L
var_miss_hc <- hclust(dist(t(miss_binary)), method = "ward.D2")
var_miss_clusters <- cutree(var_miss_hc, k = n_miss_clusters)
data_groups <- split(names(var_miss_clusters), var_miss_clusters)

# =============================================================================
# STEP 5: Imputation QC — mask 10% per missingness group, impute, compare to truth
# =============================================================================

truth <- scaled_matrix[rowSums(is.na(scaled_matrix)) == 0, ]
test  <- truth

set.seed(imp_seed)
seeds_qc <- round(runif(length(data_groups), 10, 1e6))
inds_pool <- seq(nrow(test))
masked_inds <- list()

for (i in seq_along(data_groups)) {
  grp_vars <- data_groups[[i]]
  set.seed(seeds_qc[i])
  row_inds <- sample(inds_pool, floor(nrow(test) * imp_qc_ratio_na))
  inds_pool <- inds_pool[-row_inds]
  test[row_inds, grp_vars] <- NA
  masked_inds[[i]] <- row_inds
}

inds_2d <- which(is.na(test), arr.ind = TRUE)

# run QC imputation
imp_test <- mclapply(seq(imp_m), function(n) {
  mice(test, m = 1, method = imp_method, printFlag = FALSE)
}, mc.cores = detectCores() - 1)

imp_test_avg <- Reduce("+", lapply(seq(imp_m), function(i) complete(imp_test[[i]], 1))) / imp_m

# summarise per-variable R between actual and predicted
qc_results <- data.frame(variable = character(), R = numeric(), p_value = numeric(),
                         stringsAsFactors = FALSE)
for (nm in colnames(truth)) {
  col_inds <- inds_2d[inds_2d[, 2] == which(colnames(truth) == nm), ]
  if (nrow(col_inds) == 0) next
  ct <- cor.test(truth[col_inds], imp_test_avg[col_inds], use = "complete.obs")
  qc_results <- rbind(qc_results, data.frame(variable = nm,
                                              R        = as.numeric(ct$estimate),
                                              p_value  = ct$p.value))
}

# =============================================================================
# STEP 6: Full imputation on the complete dataset
# =============================================================================

scaled_pred <- quickpred(scaled_matrix)

imp_list <- mclapply(seq(imp_m), function(n) {
  set.seed(imp_seed)
  mice(scaled_matrix, m = 1, method = imp_method,
       predictorMatrix = scaled_pred, printFlag = FALSE)
}, mc.cores = detectCores() - 1)

imp_avg <- Reduce("+", lapply(seq(imp_m), function(i) complete(imp_list[[i]], 1))) / imp_m

# =============================================================================
# STEP 7: Post-imputation constraint fix for ISOS_LOC variables
# (ISOS_LOC_JF and ISOS_LOC_CEN must be 0 when ISOS_BREAK == 0;
#  if ISOS_BREAK > 0 but both LOC columns are 0, default to JF = 1)
# =============================================================================

if (all(isos_cols %in% colnames(imp_avg))) {
  isos_raw <- imp_avg[, isos_cols]
  # rescale back to original units for rounding
  isos_rescaled <- sweep(sweep(isos_raw, 2, scale_sd[isos_cols], "*"),
                         2, scale_center[isos_cols], "+")
  isos_round <- apply(isos_rescaled, 2, round)

  # samples where ISOS_BREAK > 0 but both LOC columns == 0 → set JF = 1
  fix_samps <- rownames(isos_round)[
    isos_round[, isos_cols[1]] > 0 &
    isos_round[, isos_cols[2]] == 0 &
    isos_round[, isos_cols[3]] == 0
  ]
  val_1_scaled <- (1 - scale_center[isos_cols[2]]) / scale_sd[isos_cols[2]]
  imp_avg[fix_samps, isos_cols[2]] <- val_1_scaled
}

# =============================================================================
# STEP 8: Calculate derived features — parafoveal average and foveal slope
# =============================================================================

imp_df <- as.data.frame(imp_avg)

# parafoveal average (zones 2-5)
para_cols <- grep("MT0[2345]$", mt_zone_cols, value = TRUE)
if (length(para_cols) > 0) {
  imp_df$parafovea <- rowMeans(imp_df[, para_cols])
}

# foveal slope (degrees from zone 1 to each parafoveal zone)
calc_slope <- function(zone, zone1) {
  180 * atan((zone - zone1) / zone_distance_um) / pi
}
zone1_col <- grep("MT01$", mt_zone_cols, value = TRUE)
slope_cols <- grep("MT0[2345]$", mt_zone_cols, value = TRUE)

if (length(zone1_col) == 1 && length(slope_cols) > 0) {
  # rescale zone values back to microns before calculating slope
  zone1_um <- imp_df[[zone1_col]] * scale_sd[zone1_col] + scale_center[zone1_col]
  for (sc in slope_cols) {
    zone_um <- imp_df[[sc]] * scale_sd[sc] + scale_center[sc]
    imp_df[[paste0(sc, "_slope")]] <- calc_slope(zone_um, zone1_um)
  }
  # rescale new slope columns to z-scores
  for (nm in grep("_slope$", colnames(imp_df), value = TRUE)) {
    imp_df[[nm]] <- scale(imp_df[[nm]])[, 1]
  }
}

imp_avg <- as.matrix(imp_df)

# =============================================================================
# STEP 9: Remove poorly-imputed variables
# (non-significant QC correlation AND >5% missing in original data)
# =============================================================================

poor_vars <- qc_results$variable[
  !is.na(qc_results$p_value) &
  qc_results$p_value > 0.05 &
  (colSums(is.na(feature_matrix[, qc_results$variable, drop = FALSE])) /
     nrow(feature_matrix)) > 0.05
]
imp_avg <- imp_avg[, !colnames(imp_avg) %in% poor_vars]

# =============================================================================
# OUTPUT
# =============================================================================
# imp_avg : numeric matrix (samples x features), fully imputed and scaled,
#           with derived features appended and poorly-imputed features removed.
#           Ready for PCA in 02_clustering_trajectory.R.
