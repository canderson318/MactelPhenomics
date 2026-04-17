# MacTel Disease Trajectory Analysis — Public Analysis Scripts

Companion code for the manuscript. Scripts reproduce the core pipeline without data loading or dataset-specific paths. Run sequentially; outputs of earlier scripts feed later ones.

---

## Scripts

### 01_imputation.R
Filter features (missingness, skewness) and samples (missingness outliers). Group variables by missing-data pattern (hclust ward.D2, k=5), run parallel MICE PMM (m=11). QC via mask-and-predict. Post-imputation constraint fixes. Append parafoveal average and foveal slope features.

**Input:** `feature_matrix` (samples × features, NAs allowed), `mt_zone_cols`, `isos_cols`  
**Output:** `imp_avg` (imputed scaled matrix), `qc_results`  
**Packages:** `mice`, `psych`, `parallel`

---

### 02_clustering_trajectory.R
PCA (50 components), TSNE, UMAP on imputed matrix. Two-step clustering via `bluster::clusterCells`: K-means to 1440 micro-clusters → NNGraph (k=10) on centroids. Manual cluster merges (12→3, 4→7, 9→2). MST on PCA centroids (`TSCAN::createClusterMST`); project samples onto MST path for pseudoprogression score.

**Input:** `sce` — SingleCellExperiment with assay `"pmm_imp_of_scaled"`  
**Output:** `sce` (with `km_nn_cluster`, `pseudoprog`), `mst_pca`, `pca_line_data`, `tsne_line_data`  
**Packages:** `SingleCellExperiment`, `scater`, `bluster`, `TSCAN`, `scran`

---

### 03_xgboost_shap.R
Train multiclass XGBoost (`multi:softprob`) predicting cluster from features. Compute SHAP scores per cluster via `shapviz`. Identify most representative sample per cluster by highest predicted class probability (restricted to low-missingness samples).

**Input:** `feature_matrix`, `cluster_labels`, `missingness_ratio`  
**Output:** `fit`, `shap_obj` (mshapviz), `shap_scores` (features × clusters), `cluster_reps`  
**Packages:** `xgboost`, `shapviz`

---

### 04_network_flow_gtest.R
G-test (log-likelihood ratio) of observed cluster transition matrix against four candidate graph structures: bifurcated MST, linear-by-pseudoprogression, linear-by-Chew-grade, fully unordered. Expected matrices derived from graph distances as proximity weights (1/(1+dist)), row-normalised, scaled to observed row sums. Null distribution from 10,000 row-sum-matched random matrices.

**Input:** `transition_counts`, `mst_edges`, `cluster_pp_means`, `cluster_chew_means`  
**Output:** `results` (G-statistics per graph), `null_g_stats`, `null_quantiles`  
**Packages:** `igraph`, `dplyr`

---

### 05_linear_modeling.R
Three modelling themes:

- **Trajectory (A vs B):** `glmer` (binomial) with age/sex covariates; looped over phenotypes and SNPs (BH-adjusted)
- **Progression rate:** per-eye slopes from `lm(score ~ year)`; `lmer` against covariates, phenotypes, SNPs
- **Stepwise time:** pairwise visit differences → `lm(prog_diff ~ year_diff * stage)`; invert interaction coefficients for cumulative years per stage

**Input:** `long_dat` (one row per visit: `patID`, `patID_eye`, `year`, `trajectory`, `pseudoprog`, `chew_grade`, `AGE`, `SEXF`, phenotype/SNP/PC columns)  
**Output:** covariate models, per-phenotype and per-SNP association tables (with q-values), `pp_years_per_stage`, `chew_years_per_stage`  
**Packages:** `lme4`, `lmerTest`, `dplyr`

---

### 06_bayesian_time_modeling.R
Bayesian non-linear mixed-effects model (brms) with logarithmic growth and individual-level temporal shifts:

```
pseudoprog ~ a + b * log(c + year + shift),  shift ~ 1 + (1 | patID_eye)
```

Fit for both pseudoprogression and Chew grade. Extract per-eye shifts; use population curve to estimate years between adjacent MST clusters.

**Input:** `long_dat`, `cluster_pp_means`, `mst_edges`  
**Output:** `fit_pp`, `fit_chew`, `shift_summ`, `newdat` (population curve), `mst_edges` (with `year_diff`)  
**Packages:** `brms`, `cmdstanr`, `bayesplot`, `dplyr`

---

## Installation

```r
# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "scater", "bluster", "TSCAN", "scran"))

# CRAN
install.packages(c("mice", "psych", "xgboost", "shapviz", "igraph",
                   "dplyr", "lme4", "lmerTest", "brms", "bayesplot"))

# cmdstanr (requires separate Stan installation)
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```
