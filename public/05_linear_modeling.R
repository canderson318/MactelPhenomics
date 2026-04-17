# =============================================================================
# 05_linear_modeling.R
#
# This script applies mixed-effects linear and logistic regression models to
# relate phenotypic features, demographic covariates, and genetic variants
# to disease trajectory and progression rate. It also implements a stepwise
# time-change model that estimates the number of years that correspond to one
# unit of Chew grade or pseudoprogression change at each disease stage.
#
# Three modelling themes:
#   A. Trajectory (vascular vs. neurodegenerative route) ~ phenotypes / SNPs
#   B. Progression rate (Chew rate, pseudoprogression rate) ~ phenotypes / SNPs
#   C. Time modelling: lm(progression_diff ~ year_diff * stage)
#
# Pipeline position: runs after clustering (02_clustering_trajectory.R).
# Requires longitudinal data with trajectory labels and pseudoprogression scores.
# =============================================================================

# Required packages ------------------------------------------------------------
library(lme4)      # lmer, glmer
library(lmerTest)  # p-values for lmer
library(dplyr)     # data manipulation

# Expected input ---------------------------------------------------------------
# long_dat : data.frame with one row per visit, columns:
#   patID        — patient identifier
#   patID_eye    — eye-level identifier (patient + eye)
#   year         — visit time in years from first visit
#   trajectory   — binary factor: "A" (vascular) or "B" (neurodegenerative)
#   pseudoprog   — continuous pseudoprogression score
#   chew_grade   — ordinal Chew severity grade (0–6)
#   SEXF         — binary covariate (0/1)
#   AGE          — continuous covariate (age at baseline)
#   <phenotype columns> — continuous or binary phenotypic variables
#   <SNP columns>       — additive coded SNP dosage (0/1/2)
#   PC1:PC5      — genetic principal components (for SNP models)

# Parameters -------------------------------------------------------------------
pheno_cols <- NULL   # character vector of phenotype column names to loop over;
                     # set to NULL to skip phenotype loops
snp_cols   <- NULL   # character vector of SNP column names to loop over;
                     # set to NULL to skip genetic models
pc_cols    <- paste0("PC", 1:5)   # genetic PC column names for covariate adjustment
p_adjust_method <- "BH"           # multiple testing correction (Benjamini-Hochberg)

# =============================================================================
# THEME A: Trajectory models (binary outcome: route A vs B)
# =============================================================================

# Recode trajectory to 0/1 for binomial glmer
long_dat$traj_bin <- as.integer(long_dat$trajectory == "B")

# A1. Baseline covariates only
cov_traj_mod <- glmer(
  traj_bin ~ AGE + SEXF + (1 | patID),
  family = "binomial",
  data   = long_dat
)
summary(cov_traj_mod)

# A2. Each phenotype → trajectory (loop)
if (!is.null(pheno_cols)) {
  traj_pheno_results <- lapply(pheno_cols, function(pheno) {
    mod <- tryCatch(
      glmer(traj_bin ~ long_dat[[pheno]] + AGE + SEXF + (1 | patID),
            family = "binomial", data = long_dat),
      error = function(e) NULL
    )
    if (is.null(mod)) return(NULL)
    coef_row <- summary(mod)$coefficients["long_dat[[pheno]]", , drop = FALSE]
    data.frame(phenotype = pheno,
               estimate  = coef_row[, "Estimate"],
               p_value   = coef_row[, "Pr(>|z|)"])
  })
  traj_pheno_results <- do.call(rbind, Filter(Negate(is.null), traj_pheno_results))
  traj_pheno_results$q_value <- p.adjust(traj_pheno_results$p_value, method = p_adjust_method)
}

# A3. Each SNP → trajectory (loop)
if (!is.null(snp_cols)) {
  snp_traj_formula_rhs <- paste("AGE + SEXF +", paste(pc_cols, collapse = " + "), "+ (1 | patID)")
  snp_traj_results <- lapply(snp_cols, function(snp) {
    mod <- tryCatch(
      glmer(as.formula(paste("traj_bin ~", snp, "+", snp_traj_formula_rhs)),
            family = "binomial", data = long_dat),
      error = function(e) NULL
    )
    if (is.null(mod)) return(NULL)
    coef_row <- summary(mod)$coefficients[snp, , drop = FALSE]
    data.frame(snp      = snp,
               estimate = coef_row[, "Estimate"],
               p_value  = coef_row[, "Pr(>|z|)"])
  })
  snp_traj_results <- do.call(rbind, Filter(Negate(is.null), snp_traj_results))
  snp_traj_results$q_value <- p.adjust(snp_traj_results$p_value, method = p_adjust_method)
}

# =============================================================================
# THEME B: Progression rate models
# =============================================================================

# B1. Derive per-eye progression rates (slope of progression ~ year)
rate_dat <- long_dat |>
  dplyr::group_by(patID_eye) |>
  dplyr::filter(dplyr::n() > 1, !is.na(pseudoprog), !is.na(chew_grade)) |>
  dplyr::summarise(
    patID      = dplyr::first(patID),
    AGE        = dplyr::first(AGE),
    SEXF       = dplyr::first(SEXF),
    trajectory = dplyr::first(trajectory),
    pp_rate    = tryCatch(coef(lm(pseudoprog  ~ year))[["year"]], error = function(e) NA_real_),
    chew_rate  = tryCatch(coef(lm(chew_grade  ~ year))[["year"]], error = function(e) NA_real_),
    dplyr::across(all_of(c(if (!is.null(pheno_cols)) pheno_cols,
                            if (!is.null(snp_cols))   snp_cols,
                            pc_cols)), dplyr::first),
    .groups = "drop"
  )

# B2. Covariates → progression rates
pp_cov_mod   <- lmer(pp_rate   ~ AGE + SEXF + (1 | patID), data = rate_dat)
chew_cov_mod <- lmer(chew_rate ~ AGE + SEXF + (1 | patID), data = rate_dat)
summary(pp_cov_mod)
summary(chew_cov_mod)

# B3. Each phenotype → pseudoprogression rate (loop)
if (!is.null(pheno_cols)) {
  pp_pheno_results <- lapply(pheno_cols, function(pheno) {
    mod <- tryCatch(
      lmer(as.formula(paste("pp_rate ~", pheno, "+ AGE + SEXF + (1 | patID)")),
           data = rate_dat),
      error = function(e) NULL
    )
    if (is.null(mod)) return(NULL)
    coef_row <- summary(mod)$coefficients[pheno, , drop = FALSE]
    data.frame(phenotype = pheno,
               estimate  = coef_row[, "Estimate"],
               p_value   = coef_row[, "Pr(>|t|)"])
  })
  pp_pheno_results <- do.call(rbind, Filter(Negate(is.null), pp_pheno_results))
  pp_pheno_results$q_value <- p.adjust(pp_pheno_results$p_value, method = p_adjust_method)
}

# B4. Each SNP → phenotype (loop)
if (!is.null(snp_cols) && !is.null(pheno_cols)) {
  snp_pheno_formula_rhs <- paste(paste(pc_cols, collapse = " + "), "+ (1 | patID)")
  snp_pheno_results <- lapply(snp_cols, function(snp) {
    lapply(pheno_cols, function(pheno) {
      mod <- tryCatch(
        lmer(as.formula(paste(pheno, "~", snp, "+", snp_pheno_formula_rhs)),
             data = long_dat),
        error = function(e) NULL
      )
      if (is.null(mod)) return(NULL)
      coef_row <- summary(mod)$coefficients[snp, , drop = FALSE]
      data.frame(snp      = snp,
                 phenotype = pheno,
                 estimate  = coef_row[, "Estimate"],
                 p_value   = coef_row[, "Pr(>|t|)"])
    })
  })
  snp_pheno_results <- do.call(rbind, Filter(Negate(is.null),
                                              do.call(c, snp_pheno_results)))
  snp_pheno_results$q_value <- p.adjust(snp_pheno_results$p_value, method = p_adjust_method)
}

# =============================================================================
# THEME C: Stepwise time models — how much time does one unit of progression take?
# =============================================================================

# C1. Build pairwise forward-visit dataset per eye (expand.grid approach)
make_pairwise_diffs <- function(dat, id_col, time_col, prog_col) {
  out <- lapply(unique(dat[[id_col]]), function(eye) {
    sub <- dat[dat[[id_col]] == eye, c(time_col, prog_col)]
    sub <- sub[order(sub[[time_col]]), ]
    combos <- expand.grid(sub[[time_col]], sub[[time_col]])
    combos2 <- expand.grid(sub[[prog_col]], sub[[prog_col]])
    d <- cbind(setNames(combos, c("year", "next_year")),
               setNames(combos2, c("prog", "next_prog")))
    d[[id_col]] <- eye
    d[d$year < d$next_year, ]
  })
  do.call(rbind, out)
}

pp_diff_dat   <- make_pairwise_diffs(long_dat, "patID_eye", "year", "pseudoprog")
chew_diff_dat <- make_pairwise_diffs(long_dat, "patID_eye", "year", "chew_grade")

pp_diff_dat$year_diff <- pp_diff_dat$next_year - pp_diff_dat$year
pp_diff_dat$prog_diff <- pp_diff_dat$next_prog - pp_diff_dat$prog
pp_diff_dat <- pp_diff_dat[pp_diff_dat$prog <= pp_diff_dat$next_prog, ]

chew_diff_dat$year_diff <- chew_diff_dat$next_year - chew_diff_dat$year
chew_diff_dat$prog_diff <- chew_diff_dat$next_prog - chew_diff_dat$prog
chew_diff_dat <- chew_diff_dat[chew_diff_dat$prog <= chew_diff_dat$next_prog, ]

# quantile-bin pseudoprog for stage interaction
pp_quantiles        <- quantile(long_dat$pseudoprog, seq(0, 1, 0.1), na.rm = TRUE)
pp_diff_dat$pp_stage <- cut(pp_diff_dat$prog, pp_quantiles, labels = FALSE,
                             include.lowest = TRUE) - 1L

# C2. lm(prog_diff ~ year_diff * stage) for pseudoprogression and Chew
pp_time_mod   <- lm(prog_diff ~ year_diff * pp_stage, data = pp_diff_dat)
chew_time_mod <- lm(prog_diff ~ year_diff * prog,     data = chew_diff_dat)

summary(pp_time_mod)
summary(chew_time_mod)

# C3. Invert interaction coefficients → cumulative years per stage
invert_to_years <- function(mod, interaction_pattern = "year_diff") {
  coefs <- coef(summary(mod))[, "Estimate"]
  year_coefs <- coefs[grepl(interaction_pattern, names(coefs))]
  # baseline slope + each stage's additional slope
  full_effects <- c(year_coefs[1], year_coefs[1] + cumsum(year_coefs[-1]))
  recip        <- 1 / full_effects
  cumsum(recip)   # cumulative years to each stage
}

pp_years_per_stage   <- invert_to_years(pp_time_mod)
chew_years_per_stage <- invert_to_years(chew_time_mod)

cat("Estimated cumulative years to each pseudoprogression stage:\n")
print(round(pp_years_per_stage, 2))
cat("Estimated cumulative years to each Chew grade:\n")
print(round(chew_years_per_stage, 2))

# =============================================================================
# OUTPUT
# =============================================================================
# cov_traj_mod         : glmer, trajectory ~ age + sex
# traj_pheno_results   : data.frame, per-phenotype trajectory association
# snp_traj_results     : data.frame, per-SNP trajectory association
# pp_cov_mod           : lmer, pp_rate ~ covariates
# chew_cov_mod         : lmer, chew_rate ~ covariates
# pp_pheno_results     : data.frame, per-phenotype pp_rate association
# snp_pheno_results    : data.frame, per-SNP per-phenotype association
# pp_time_mod          : lm, pp_diff ~ year_diff * pp_stage
# chew_time_mod        : lm, chew_diff ~ year_diff * chew_stage
# pp_years_per_stage   : cumulative years to each pseudoprogression quantile
# chew_years_per_stage : cumulative years to each Chew grade
