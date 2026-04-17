# =============================================================================
# 03_xgboost_shap.R
#
# This script trains a multiclass XGBoost model predicting cluster membership
# from phenotypic features, then computes SHAP (SHapley Additive exPlanations)
# scores to profile each cluster. SHAP scores reveal which features are most
# important for distinguishing each cluster from the rest. The script also
# identifies the most representative sample per cluster using predicted class
# probabilities, restricted to samples with low raw missingness.
#
# =============================================================================

# Required packages ------------------------------------------------------------
library(xgboost)   # gradient boosted trees
library(shapviz)   # SHAP visualisation

# Expected input ---------------------------------------------------------------
# feature_matrix : numeric matrix (samples x features), imputed and scaled —
#                  the "pmm_imp_of_scaled" assay from 02_clustering_trajectory.R.
#                  Row names = sample IDs.
#
# cluster_labels : named factor or character vector of cluster assignments,
#                  one entry per sample (same order as rows of feature_matrix).
#                  Corresponds to sce$my_cluster from 02_clustering_trajectory.R.
#
# missingness_ratio : named numeric vector (one value per sample) of original
#                     (pre-imputation) feature missingness fractions.
#                     Used to restrict representative selection to "clean" samples.

# Parameters -------------------------------------------------------------------
xgb_nrounds           <- 50     # boosting rounds
xgb_nthread           <- 1      # threads per model fit
max_display_shap      <- 20     # features shown in beeswarm importance plot
missingness_threshold <- 0.20   # samples above this fraction excluded from representative selection

# =============================================================================
# STEP 1: Prepare data
# =============================================================================

assay <- as.data.frame(feature_matrix)
assay$my_cluster <- cluster_labels

# XGBoost requires 0-indexed integer labels
n_classes  <- length(unique(cluster_labels))
label_vec  <- as.integer(factor(cluster_labels)) - 1L

dtrain <- xgb.DMatrix(
  data   = data.matrix(assay[, colnames(assay) != "my_cluster"]),
  label  = label_vec,
  nthread = xgb_nthread
)

# =============================================================================
# STEP 2: Train XGBoost multiclass model
# =============================================================================

params <- list(
  objective  = "multi:softprob",
  num_class  = n_classes,
  nthread    = xgb_nthread
)

fit <- xgb.train(params = params, data = dtrain, nrounds = xgb_nrounds)

# =============================================================================
# STEP 3: Compute SHAP scores
# =============================================================================

shap_obj <- shapviz(
  fit,
  X_pred = data.matrix(assay[, colnames(assay) != "my_cluster"]),
  X      = assay[, colnames(assay) != "my_cluster"]
)

# name each shapviz object by cluster label
names(shap_obj) <- levels(factor(cluster_labels))

# =============================================================================
# STEP 4: Visualise SHAP importance
# =============================================================================

sv_importance(shap_obj, max_display = max_display_shap, kind = "beeswarm")

# mean |SHAP| per feature per cluster as a matrix (features x clusters)
shap_scores <- sv_importance(shap_obj, kind = "no")

# =============================================================================
# STEP 5: Identify most representative sample per cluster
# (highest predicted class probability among low-missingness samples)
# =============================================================================

pred_vec  <- predict(fit, newdata = dtrain)
pred_mat  <- matrix(pred_vec, nrow = nrow(dtrain), ncol = n_classes, byrow = TRUE)
rownames(pred_mat) <- rownames(feature_matrix)

# restrict to samples with low raw missingness
clean_samps <- names(missingness_ratio)[missingness_ratio < missingness_threshold]
pred_mat_clean <- pred_mat[rownames(pred_mat) %in% clean_samps, ]

# sample with highest predicted probability for each cluster
cluster_reps <- setNames(
  rownames(pred_mat_clean)[apply(pred_mat_clean, 2, which.max)],
  levels(factor(cluster_labels))
)

# =============================================================================
# OUTPUT
# =============================================================================
# fit          : trained xgb.Booster object
# shap_obj     : mshapviz object (one shapviz per cluster)
# shap_scores  : matrix (features x clusters) of mean absolute SHAP values
# cluster_reps : named character vector — most representative sample per cluster
