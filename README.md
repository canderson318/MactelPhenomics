# MacTel Disease Trajectory Analysis — Public Analysis Scripts

Companion code for the manuscript. These six R scripts reproduce the core analytical pipeline: imputation, clustering, feature profiling, network validation, and progression modeling. No data loading, file paths, or dataset-specific code is included — the user supplies input objects as described in each script's header.

---

## Pipeline Overview

| Script | Purpose |
|--------|---------|
| `01_imputation.R` | Feature/sample filtering, grouped MICE PMM imputation, QC, derived variables |
| `02_clustering_trajectory.R` | PCA, two-step clustering (K-means → NNGraph), MST, pseudoprogression |
| `03_xgboost_shap.R` | XGBoost cluster profiling + SHAP feature importance |
| `04_network_flow_gtest.R` | G-test of observed cluster transitions vs. candidate graph structures |
| `05_linear_modeling.R` | Mixed-effects models: disease trajectory, progression rate, stepwise time |
| `06_bayesian_time_modeling.R` | Bayesian logarithmic growth model with individual-level temporal shifts |

Scripts are designed to run sequentially. Outputs of earlier scripts feed later ones.

---

## Script Outlines

### 01_imputation.R

Filters a raw feature matrix by missingness and skewness, groups variables by their missing-data pattern, and runs parallel MICE PMM imputation. Includes a QC step (mask-and-predict) and post-imputation constraint fixes for ordinal/binary variables. Appends derived features (parafoveal average, foveal slope).

**Input:**
- `feature_matrix` — numeric matrix (samples × features), NAs allowed
- `mt_zone_cols` — macular thickness column names (for derived features)
- `isos_cols` — ISOS break/location column names (for constraint fixes)

**Output:**
- `imp_avg` — fully imputed, scaled matrix ready for PCA
- `qc_results` — data.frame of per-variable imputation R and p-value

**Packages:** `mice`, `psych`, `parallel`

---

### 02_clustering_trajectory.R

Takes an imputed feature matrix (in a SingleCellExperiment container) and runs PCA, TSNE, and UMAP for dimensionality reduction. Clusters observations using a two-step approach: K-means to 1440 micro-clusters, followed by nearest-neighbour graph clustering on centroids. Merges three cluster pairs identified by inspection. Builds a minimum spanning tree (MST) on PCA centroids and projects each sample onto the MST path to derive a continuous pseudoprogression score.

**Input:**
- `sce` — SingleCellExperiment object with assay `"pmm_imp_of_scaled"` (features × samples)

**Output:**
- `sce` — updated SCE with `km_nn_cluster` and `pseudoprog` in `colData`; PCA/TSNE/UMAP in `reducedDims`
- `mst_pca` — igraph MST object
- `pca_line_data`, `tsne_line_data` — MST edge coordinates for plotting

**Packages:** `SingleCellExperiment`, `scater`, `bluster`, `TSCAN`, `scran`

---

### 03_xgboost_shap.R

Trains a multiclass XGBoost model predicting cluster membership from phenotypic features. Computes SHAP scores (via `shapviz`) to profile which features drive each cluster's identity. Identifies the most representative sample per cluster as the observation with the highest predicted class probability among low-missingness samples.

**Input:**
- `feature_matrix` — numeric matrix (samples × features), imputed and scaled
- `cluster_labels` — named character/factor vector of cluster assignments
- `missingness_ratio` — named numeric vector of per-sample raw missingness fractions

**Output:**
- `fit` — trained `xgb.Booster` object
- `shap_obj` — `mshapviz` object (one `shapviz` per cluster)
- `shap_scores` — matrix (features × clusters) of mean absolute SHAP values
- `cluster_reps` — named character vector of most representative sample per cluster

**Packages:** `xgboost`, `shapviz`

---

### 04_network_flow_gtest.R

Tests whether observed longitudinal cluster transitions are better explained by the bifurcated MST trajectory than by alternative graph structures (linear by pseudoprogression, linear by Chew grade, fully unordered). For each candidate graph, derives an expected transition matrix using proximity weights (1/(1+distance)), row-normalised and scaled to observed row sums. Computes the G-statistic (log-likelihood ratio) for each. Builds a null distribution from 10,000 semi-random matrices preserving observed row sums.

**Input:**
- `transition_counts` — integer matrix (clusters × clusters) of observed transition counts
- `mst_edges` — data.frame with columns `from` and `to` (MST edges)
- `cluster_pp_means` — data.frame with `my_cluster` and `mean_pseudoprog`
- `cluster_chew_means` — data.frame with `my_cluster` and `mean_chew`

**Output:**
- `results` — data.frame of G-statistics per candidate graph
- `null_g_stats` — numeric vector of G-statistics under the null
- `null_quantiles` — 2.5 / 50 / 97.5 percentiles of the null distribution

**Packages:** `igraph`, `dplyr`

---

### 05_linear_modeling.R

Three modelling themes applied to longitudinal data:

- **A. Trajectory** — `glmer` (binomial) of vascular vs. neurodegenerative route on age/sex covariates, then looped over phenotypes and SNPs (BH-adjusted).
- **B. Progression rate** — per-eye progression rates derived from `lm(score ~ year)` slopes, then `lmer` models relating rates to covariates, phenotypes, and SNPs.
- **C. Stepwise time** — pairwise forward-visit differences modelled as `lm(prog_diff ~ year_diff * stage)`; interaction coefficients inverted to yield cumulative years per disease stage.

**Input:**
- `long_dat` — data.frame (one row per visit) with `patID`, `patID_eye`, `year`, `trajectory`, `pseudoprog`, `chew_grade`, `AGE`, `SEXF`, phenotype columns, SNP columns, `PC1`–`PC5`

**Output:**
- `cov_traj_mod`, `pp_cov_mod`, `chew_cov_mod` — baseline covariate models
- `traj_pheno_results`, `pp_pheno_results` — per-phenotype association tables
- `snp_traj_results`, `snp_pheno_results` — per-SNP association tables
- `pp_time_mod`, `chew_time_mod` — stepwise time models
- `pp_years_per_stage`, `chew_years_per_stage` — cumulative years to each stage

**Packages:** `lme4`, `lmerTest`, `dplyr`

---

### 06_bayesian_time_modeling.R

Fits a Bayesian non-linear mixed-effects model (`brms`) with logarithmic growth form:

```
pseudoprog ~ a + b * log(c + year + shift),  shift ~ 1 + (1 | patID_eye)
```

Each eye receives an individual-level temporal shift, anchoring it on the population growth curve. The fitted population curve is then used to read off the estimated time separating adjacent clusters on the MST. The same model is fit for Chew grade.

**Input:**
- `long_dat` — data.frame with `patID_eye`, `year`, `pseudoprog`, `chew_grade`
- `cluster_pp_means` — data.frame with `my_cluster` and `mean_pseudoprog`
- `mst_edges` — data.frame with `from` and `to`

**Output:**
- `fit_pp`, `fit_chew` — `brmsfit` objects
- `shift_summ` — per-eye temporal shift estimates with 95% CI
- `newdat` — population-level prediction curve
- `mst_edges` — updated with `year_diff` column per edge
- `cluster_pp_means` — updated with `estimated_year` column

**Packages:** `brms`, `cmdstanr`, `bayesplot`, `dplyr`

---

## Dependencies

```r
# CRAN packages
install.packages(c(
  "mice", "psych", "parallel",
  "SingleCellExperiment", "scater", "bluster", "TSCAN", "scran",
  "xgboost", "shapviz",
  "igraph", "dplyr",
  "lme4", "lmerTest",
  "brms", "bayesplot"
))

# Bioconductor packages (SingleCellExperiment, scater, bluster, TSCAN, scran)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "scater", "bluster", "TSCAN", "scran"))

# cmdstanr requires a separate Stan installation
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```
