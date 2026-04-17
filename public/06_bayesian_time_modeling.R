# =============================================================================
# 06_bayesian_time_modeling.R
#
# This script models the relationship between pseudoprogression and calendar
# time using a Bayesian non-linear mixed-effects model (brms). Each patient's
# eye is assigned an individual-level temporal shift (random intercept on the
# time axis), effectively asking: "given where this eye is in disease space,
# how long ago did their disease start?" The model is logarithmic growth:
#
#   pseudoprog ~ a + b * log(c + year + shift)
#   shift ~ 1 + (1 | patID_eye)
#
# Once fit, the population-level curve is used to infer how many years
# separate adjacent clusters on the MST. The same procedure is repeated for
# Chew grade as a complementary disease staging measure.
#
# =============================================================================

# Required packages ------------------------------------------------------------
library(brms)       # Bayesian regression models
library(cmdstanr)   # Stan backend (must be installed separately)
library(bayesplot)  # posterior diagnostic plots
library(dplyr)      # data manipulation

# Expected input ---------------------------------------------------------------
# long_dat : data.frame with one row per visit, columns:
#   patID_eye    — eye-level identifier
#   year         — visit time in years
#   pseudoprog   — continuous pseudoprogression score (will be z-scored internally)
#   chew_grade   — ordinal Chew severity grade (optional; for parallel Chew model)
#
# cluster_pp_means : data.frame with columns "my_cluster" and "mean_pseudoprog",
#                    cluster mean pseudoprogression scores (from 02_clustering_trajectory.R).
#
# mst_edges        : data.frame with columns "from" and "to", listing MST edges
#                    (from 02_clustering_trajectory.R), for reading off cluster sequence.

# Parameters -------------------------------------------------------------------
brm_chains  <- 4
brm_iter    <- 2000
brm_cores   <- parallel::detectCores() - 1
brm_backend <- "cmdstanr"

# =============================================================================
# STEP 1: Prepare pseudoprogression modeling data
# =============================================================================

raw_dat <- long_dat |>
  dplyr::select(patID_eye, year, pseudoprog) |>
  dplyr::group_by(patID_eye) |>
  dplyr::filter(dplyr::n() > 1) |>
  dplyr::mutate(pseudoprog = cummax(pseudoprog)) |>  # enforce monotonicity
  dplyr::ungroup() |>
  dplyr::arrange(year)

# z-score pseudoprog (store attributes for back-transformation)
pp_scaled        <- scale(raw_dat$pseudoprog)
pp_scale_center  <- attr(pp_scaled, "scaled:center")
pp_scale_sd      <- attr(pp_scaled, "scaled:scale")
raw_dat$pseudoprog <- as.numeric(pp_scaled)

modeling_dat <- raw_dat

# =============================================================================
# STEP 2: Define and fit Bayesian logarithmic growth model with temporal shifts
# =============================================================================

form <- brms::bf(
  pseudoprog ~ a + b * log(c + (year + shift)),
  a + b + c ~ 1,
  shift ~ 1 + (1 | patID_eye),   # individual temporal shifts
  nl = TRUE
)

priors <- c(
  prior(normal(0, 2),  nlpar = "a"),
  prior(normal(1, 1),  nlpar = "b"),
  prior(normal(1, 10), nlpar = "c", lb = 1),
  prior(normal(5, 10), nlpar = "shift")
)

fit_pp <- brm(
  form,
  prior   = priors,
  data    = modeling_dat,
  chains  = brm_chains,
  cores   = brm_cores,
  iter    = brm_iter,
  algorithm = "sampling",
  backend = brm_backend,
  silent  = 2
)

# =============================================================================
# STEP 3: Extract individual-level temporal shifts
# =============================================================================

post_summ    <- posterior_summary(fit_pp)
global_shift <- if ("b_shift_Intercept" %in% rownames(post_summ)) {
  post_summ["b_shift_Intercept", "Estimate"]
} else { 0 }

shift_summ <- post_summ[grep("r_patID_eye", rownames(post_summ)), , drop = FALSE] |>
  as.data.frame()
shift_summ[, c("Estimate", "Q2.5", "Q97.5")] <-
  shift_summ[, c("Estimate", "Q2.5", "Q97.5")] + global_shift
shift_summ$patID_eye <- sub(".*\\[([^,]+).*", "\\1", rownames(shift_summ))

modeling_dat <- dplyr::left_join(modeling_dat, shift_summ, by = "patID_eye")
modeling_dat$shift_year <- modeling_dat$year + modeling_dat$Estimate

# normalise so values are positive
abs_min <- abs(min(modeling_dat$shift_year, na.rm = TRUE))
modeling_dat$shift_year_norm <- modeling_dat$shift_year + abs_min

# =============================================================================
# STEP 4: Population-level prediction curve (line of best fit)
# =============================================================================

rescale_pp <- function(x) pp_scale_sd * x + pp_scale_center

newdat <- data.frame(
  shift_year_norm = seq(min(modeling_dat$shift_year_norm, na.rm = TRUE),
                        max(modeling_dat$shift_year_norm, na.rm = TRUE),
                        length.out = 1000)
)
newdat$year <- newdat$shift_year_norm - abs_min - global_shift
newdat$patID_eye <- NA   # marginal predictions (no random effect)

preds              <- fitted(fit_pp, newdata = newdat, re_formula = NA)
newdat$pseudoprog_fit <- rescale_pp(preds[, "Estimate"])
newdat$lower          <- rescale_pp(preds[, "Q2.5"])
newdat$upper          <- rescale_pp(preds[, "Q97.5"])

# linear approximation of curve slope (pp per year)
lobf <- lm(pseudoprog_fit ~ shift_year_norm, data = newdat)
pp_per_year <- coef(lobf)[["shift_year_norm"]]
yr_per_pp   <- 1 / pp_per_year

# =============================================================================
# STEP 5: Estimate cluster-to-cluster time from population curve
# =============================================================================

# function: look up year for a given pseudoprogression value along the curve
year_from_pp <- function(pp_val) {
  approx(x = newdat$pseudoprog_fit, y = newdat$shift_year_norm, xout = pp_val)$y
}

cluster_pp_means$estimated_year <- year_from_pp(cluster_pp_means$mean_pseudoprog)

# year difference along each MST edge
mst_edges$year_diff <- mapply(
  function(f, t) {
    year_from_pp(cluster_pp_means$mean_pseudoprog[cluster_pp_means$my_cluster == t]) -
    year_from_pp(cluster_pp_means$mean_pseudoprog[cluster_pp_means$my_cluster == f])
  },
  mst_edges$from,
  mst_edges$to
)

cat("Estimated years along each MST edge:\n")
print(mst_edges)

# =============================================================================
# STEP 6: Chew grade model (parallel structure)
# =============================================================================

chew_dat <- long_dat |>
  dplyr::select(patID_eye, year, chew_grade) |>
  dplyr::group_by(patID_eye) |>
  dplyr::filter(dplyr::n() > 1, !is.na(chew_grade)) |>
  dplyr::ungroup() |>
  dplyr::arrange(year)

chew_scaled      <- scale(chew_dat$chew_grade)
chew_scale_center <- attr(chew_scaled, "scaled:center")
chew_scale_sd     <- attr(chew_scaled, "scaled:scale")
chew_dat$chew_grade <- as.numeric(chew_scaled)

form_chew <- brms::bf(
  chew_grade ~ a + b * log(c + (year + shift)),
  a + b + c ~ 1,
  shift ~ 1 + (1 | patID_eye),
  nl = TRUE
)

fit_chew <- brm(
  form_chew,
  prior   = priors,
  data    = chew_dat,
  chains  = brm_chains,
  cores   = brm_cores,
  iter    = brm_iter,
  algorithm = "sampling",
  backend = brm_backend,
  silent  = 2
)

# =============================================================================
# OUTPUT
# =============================================================================
# fit_pp             : brmsfit object — pseudoprogression temporal shift model
# fit_chew           : brmsfit object — Chew grade temporal shift model
# shift_summ         : data.frame — per-eye temporal shift estimates with 95% CI
# newdat             : data.frame — population-level prediction curve
# pp_per_year        : numeric — estimated pseudoprogression units per year
# yr_per_pp          : numeric — estimated years per pseudoprogression unit
# mst_edges          : data.frame — MST edges with estimated year_diff column
# cluster_pp_means   : data.frame — cluster means with estimated_year column
