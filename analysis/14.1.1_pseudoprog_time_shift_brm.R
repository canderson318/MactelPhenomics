   
# SETUP ----
options(max.print = 100)
# Packages
pacman::p_load(tidyr,reshape2,purrr,tidyverse,ggplot2,gridExtra,ComplexHeatmap,viridis,zoo,mgcv,mice,ggrepel,htmltools,Metrics,psych,car,stats,lme4,lmerTest,ape,SingleCellExperiment,scater,dplyr,HDF5Array,parallel, brms,cmdstanr,txtplot, magrittr, bayesplot)


rm(list = ls())
gc()

# functions
# script with all functions i have made
source("/Users/canderson/Documents/LMRI/PRJ2024001/version015/analysis/00_MyFunctions.R")


# LOAD modeling_datA ----
wd_path = "~/Documents/LMRI/PRJ2024001/version015/"

how.long()
# load experiment
sce <-  loadHDF5SummarizedExperiment(paste0(wd_path, "processed_data/sce006.h5"))
how.long()


data.frame(colData(sce)) %>% 
  select(patID_eye, year, pseudoprog) %>% 
  arrange(patID_eye, year) %>% 
  mutate(pseudoprog = cummax(pseudoprog), .by = patID_eye) %>% 
  drop_na() %>% 
  ggplot(aes(year,  pseudoprog ))+
  geom_line(aes(group = patID_eye), alpha = .2)+
  scale_color_viridis_c()+
  ggtitle("Pseudoprog Over Time", "Each line an eye")+
  xlab("Pseudoprogression")+
  ylab("Year")+
  theme_classic()

sprintf("%sresults/clustering/14.1_progressionTimeModeling/pp_over_time.pdf", wd_path) %>% 
  {ggsave(., height = 5, width = 7); sprintf("cp %s %s", .,sprintf("%sresults/publication_figures/%s",wd_path, basename(.)) ) %>% system()}

# MAKE modeling_datA FOR MODELING----
raw_modeling_modeling_dat <- as.data.frame(colData(sce)) %>% 
  select(sample,patID_eye, month, pseudoprog,BL_eye_age) %>% 
  mutate(year = month/12) %>% 
  group_by(patID_eye) %>% 
  filter(n()>1) %>%
  mutate(pseudoprog = cummax(pseudoprog)) %>% 
  ungroup() %>% 
  arrange(year)

pseudoprog_scale <- scale(raw_modeling_modeling_dat$pseudoprog)
pseudoprog_scale_attr <- attributes(pseudoprog_scale)
raw_modeling_modeling_dat$pseudoprog <- as.numeric(pseudoprog_scale)

bl_pseudoprog <- raw_modeling_modeling_dat %>% group_by(patID_eye) %>% dplyr::slice(1) 

bl_qs <- quantile(bl_pseudoprog$pseudoprog  ,  seq(0,1,.2),   na.rm = TRUE)
# bl_qs <- seq(min(bl_pseudoprog$pseudoprog, na.rm = TRUE), max(bl_pseudoprog$pseudoprog, na.rm = TRUE), length.out = 6)
bl_pseudoprog$bl_q <- cut(bl_pseudoprog$pseudoprog, breaks = bl_qs,labels = FALSE, include.lowest = TRUE)

raw_modeling_modeling_dat <- left_join(raw_modeling_modeling_dat, bl_pseudoprog %>% select(patID_eye,bl_q), by = 'patID_eye')

raw_modeling_modeling_dat$bl_q <- raw_modeling_modeling_dat$bl_q-1
# raw_modeling_modeling_dat$bl_q <- factor(raw_modeling_modeling_dat$bl_q,levels = 0:n_distinct(raw_modeling_modeling_dat$bl_q))


raw_modeling_modeling_dat %>% 
  drop_na() %>% 
  ggplot(aes(year,  pseudoprog ))+
  geom_line(aes(group = patID_eye), alpha = .3)+
  # stat_smooth(aes(group = patID_eye),method = "lm", se = FALSE, color = rgb(0,0,0,.2))+
  geom_smooth(aes(group = bl_q, color = bl_q),method = "lm")+
  # facet_wrap(~bl_q)
  theme_bw()


# MODEL ----

modeling_dat <- raw_modeling_modeling_dat

# # test d
# set.seed(10303)
# x <- modeling_dat[modeling_dat$patID_eye%in% sample(unique(modeling_dat$patID_eye), 100), ]


# # estimate logistic growth ----
# form <- bf(
#   pseudoprog ~ a + b * exp(-c * (year + shift)),
#   a+b+c ~ 1,
#   # shift ~0+ (1|bl_q),
#   shift ~0+ (1|patID_eye),
#   nl = TRUE
# )
# 
# priors <- c(
#   prior(normal(0, 2), nlpar = "a"),
#   prior(normal(1, 2), nlpar = "b"),
#   prior(normal(1, 0.5), nlpar = "c", lb = 0),
#   # prior(student_t(0,3,5),class = "sd", nlpar = "shift")
#   prior(normal(0,5),class = "sd", nlpar = "shift")
# )
# 
# fit <- brm(
#   form,
#   data= modeling_dat,
#   # data= x,
#   prior = priors,
#   chains = 4,
#   cores = parallel::detectCores() - 1,
#   iter = 2000,
#   algorithm = "sampling",
#   backend = "cmdstanr",
#   silent = 2
# )
# saveRDS(fit, sprintf("%sresults/pseudoprog_time_shift_brm_wpriors.rds", wd_path))
# logistic_fit <- readRDS(sprintf("%sresults/pseudoprog_time_shift_brm_wpriors.rds", wd_path))

# estimate logarithmic growth with c ----
# form <- bf(
#   pseudoprog ~ a + b * log(c + (year + shift)),
#   a+b+c ~ 1,
#   shift ~0+ (1|patID_eye),
#   nl = TRUE
# )
# 
# priors <- c(
#   prior(normal(0, 2), nlpar = "a"),
#   prior(normal(1, 1), nlpar = "b"),
#   prior(normal(0, 10), nlpar = "c", lb = 1),
#   # prior(normal(0, 5), class = "sd", nlpar = "shift")
#   prior(normal(0, 5), class = "sd", nlpar = "shift", lb = 0)
# )
# 
# fit <- brm(
#   form,
#   prior = priors,
#   data= modeling_dat,
#   chains = 4,
#   cores = parallel::detectCores() - 1,
#   iter = 2000,
#   algorithm = "sampling",
#   backend = "cmdstanr",
#   silent = 2
# ); saveRDS(fit, sprintf("%sresults/logorithmic_pseudoprog_time_shift_brm_wpriors.rds", wd_path))
logorithmic_fit <- readRDS(sprintf("%sresults/logorithmic_pseudoprog_time_shift_brm_wpriors.rds", wd_path))

# estimate logarithmic growth with c and global shift intercept ----
# form <- bf(
#   # keep c because allows for flexibility; shift in when pseudoprog starts increasing
#   pseudoprog ~ a + b * log(c + (year + shift)),
#   a+b+c ~ 1,
#   shift ~ 1+ (1|patID_eye), # include global intercept, each patient intercept becomes difference from global
#   nl = TRUE
# )
# 
# priors <- c(
#   prior(normal(0, 1), nlpar = "a"),
#   prior(normal(1, 1), nlpar = "b"),
#   prior(normal(1, 10), nlpar = "c", lb = 1),
#   prior(normal(5, 10), nlpar = "shift") # global mean shift
# )
# 
# fit <- brm(
#   form,
#   prior = priors,
#   data= modeling_dat,
#   chains = 4,
#   cores = parallel::detectCores() - 1,
#   iter = 2000,
#   algorithm = "sampling",
#   backend = "cmdstanr",
#   silent = 1
# );
# saveRDS(fit, sprintf("%sresults/global_shift_intercept_logorithmic_pseudoprog_time_shift_brm_wpriors.rds", wd_path))
glbl_shft_intrcpt_logorithmic_fit <- readRDS(sprintf("%sresults/global_shift_intercept_logorithmic_pseudoprog_time_shift_brm_wpriors.rds", wd_path))


# estimate logarithmic growth with c and global shift intercept; new priors ----
# form <- bf(
#   # keep c because allows for flexibility; shift in when pseudoprog starts increasing
#   pseudoprog ~ a + b * log(c + (year + shift)),
#   a+b+c ~ 1,
#   shift ~ 1+ (1|patID_eye), # include global intercept, each patient intercept becomes difference from global
#   nl = TRUE
# )
# 
# priors <- c(
#   prior(normal(0, 3), nlpar = "a"),
#   prior(normal(3, 3), nlpar = "b"),
#   prior(normal(5, 20), nlpar = "c", lb = 1),
#   prior(normal(10, 20), nlpar = "shift")
# )
# 
# print(form)
# priors
# fit <- brm(
#   form,
#   prior = priors,
#   data= modeling_dat,
#   chains = 4,
#   cores = parallel::detectCores() - 1,
#   iter = 2000,
#   algorithm = "sampling",
#   backend = "cmdstanr",
#   silent = 2
# );
# saveRDS(fit, sprintf("%sresults/new_priors_global_shift_intercept_logorithmic_pseudoprog_time_shift_brm_wpriors.rds", wd_path))
new_priors_glbl_shft_intrcpt_logorithmic_fit <- readRDS(sprintf("%sresults/new_priors_global_shift_intercept_logorithmic_pseudoprog_time_shift_brm_wpriors.rds", wd_path))


# prior only estimate logarithmic growth with c and global shift intercept ----
# form <- bf(
#   # keep c because allows for flexibility; shift in when pseudoprog starts increasing
#   pseudoprog ~ a + b * log(c + (year + shift)),
#   a+b+c ~ 1,
#   shift ~ 1+ (1|patID_eye), # include global intercept, each patient intercept becomes difference from global
#   nl = TRUE
# )
# 
# priors <- c(
#   prior(normal(0, 2), nlpar = "a"),
#   prior(normal(1, 1), nlpar = "b"),
#   prior(normal(1, 1), nlpar = "c", lb = 1),
#   prior(normal(5, 10), nlpar = "shift"), # global mean shift
#   prior(normal(0, 10), nlpar = "shift", class = "sd", lb = 0), # global mean shift sd
# )
# 
# fit <- brm(
#   form,
#   prior = priors,
#   sample_prior = "only",
#   data= modeling_dat,
#   chains = 4,
#   cores = parallel::detectCores() - 1,
#   iter = 2000,
#   algorithm = "sampling",
#   backend = "cmdstanr",
#   silent = 1
# ); 
# saveRDS(fit, sprintf("%sresults/prior_only_global_shift_intercept_logorithmic_pseudoprog_time_shift_brm_wpriors.rds", wd_path))
prior_only_glbl_shft_intrcpt_logorithmic_fit <- readRDS(sprintf("%sresults/prior_only_global_shift_intercept_logorithmic_pseudoprog_time_shift_brm_wpriors.rds", wd_path))

# estimate logarithmic growth without c ----
# form <- bf(
#   pseudoprog ~ a + b * log(1+year + shift),
#   a+b+c ~ 1,
#   shift ~ 1+ (1|patID_eye), # include global intercept, each patient intercept becomes difference from global
#   nl = TRUE
# )
# 
# priors <- c(
#   prior(normal(0, 2), nlpar = "a"),
#   prior(normal(1, 1), nlpar = "b"),
#   prior(normal(5, 10), nlpar = "shift") # global mean shift
# )
# 
# fit <- brm(
#   form,
#   prior = priors,
#   data= modeling_dat,
#   chains = 4,
#   cores = parallel::detectCores() - 1,
#   iter = 2000,
#   algorithm = "sampling",
#   backend = "cmdstanr",
#   silent = 1
# );
# saveRDS(fit, sprintf("%sresults/no_c_global_intercept_logorithmic_pseudoprog_time_shift_brm_wpriors.rds", wd_path))
# no_c_logorithmic_fit <- readRDS(sprintf("%sresults/no_c_global_intercept_logorithmic_pseudoprog_time_shift_brm_wpriors.rds", wd_path))

# # estimate sqrt growth ----
# form <- bf(
#   # global intercept+global rate * log(year with patient level shift)
#   pseudoprog ~ a + b * sqrt(year + shift ),
#   a+b ~ 1,
#   shift ~0+ (1|patID_eye),
#   nl = TRUE
# )
# get_prior(form, modeling_dat)
# 
# priors <- c(
#   # Try to keep intercept at zero; when pp 0, year 0
#   prior(normal(0,.01), nlpar = "a"),
#   # bias for positive rate; we observe a positive relationship
#   prior(normal(.1, 1), nlpar = "b"),
#   # bias for positive shifts; patients should be 'later', i.e. their disease started before their first visit
#   prior(exponential(5), class = "sd", nlpar = "shift")
# )
# 
# fit <- brm(
#   form,
#   prior = priors,
#   data= modeling_dat,
#   chains = 4,
#   cores = parallel::detectCores() - 1,
#   iter = 2000,
#   algorithm = "sampling",
#   backend = "cmdstanr",
#   silent = 1
# ); saveRDS(fit, sprintf("%sresults/sqrt_pseudoprog_time_shift_brm_wpriors.rds", wd_path))
# 
# sqrt_fit <- readRDS(sprintf("%sresults/sqrt_pseudoprog_time_shift_brm_wpriors.rds", wd_path))

# # estimate linear pattern----
# form <- bf(
#   pseudoprog ~ b * (year + shift), # each individual curve is function of specific shift and time
#   b ~ 1, # estimate common slope
#   shift ~ 0 + (1|patID_eye), # make shift eye specific
#   nl = TRUE
# )
#
# priors <- c(
#   prior(normal(.2, .5), nlpar = "b"),
#   prior(normal(5,5), class = "sd", nlpar = "shift")
# )
#
# fit <- brm(
#   form,
#   data= modeling_dat,
#   prior = priors,
#   chains = 5,
#   cores = parallel::detectCores() - 1,
#   iter = 2000,
#   algorithm = "sampling",
#   # algorithm = "meanfield",
#   backend = "cmdstanr",
#   silent = 2
#   )
# saveRDS(linear_fit, sprintf("%sresults/linear_pseudoprog_time_shift_brm_wpriors.rds", wd_path))
linear_fit <- readRDS(sprintf("%sresults/linear_pseudoprog_time_shift_brm_wpriors.rds", wd_path))

# # estimate linear pattern with only positive shifts allowed ----
# form <- bf(
#   pseudoprog ~ b * (year + exp(shift)), # each individual curve is function of specific shift and time
#   b ~ 1, # estimate common slope
#   shift ~ 0 + (1|patID_eye), # make shift eye specific
#   nl = TRUE
# )
# 
# priors <- c(
#   prior(normal(.2, .5), nlpar = "b"),
#   prior(normal(0,1.5), class = "sd", nlpar = "shift")
# )
# 
# fit <- brm(
#   form,
#   data= modeling_dat,
#   prior = priors,
#   chains = 5,
#   cores = parallel::detectCores() - 1,
#   iter = 2000,
#   algorithm = "sampling",
#   backend = "cmdstanr",
#   silent = 2
#   )
# saveRDS(fit, sprintf("%sresults/pos_shift_linear_pseudoprog_time_shift_brm_wpriors.rds", wd_path))
# pos_shift_linear_fit <- readRDS(sprintf("%sresults/pos_shift_linear_pseudoprog_time_shift_brm_wpriors.rds", wd_path))

# # estimate linear pattern with quantile interaction ----
# form <- bf(
#   pseudoprog ~ b * (year + shift), # each individual curve is function of specific shift and time
#   b ~ 1+bl_q, # estimate common slope for each quintile
#   shift ~ 0 + (1|patID_eye), # make shift eye-specific
#   nl = TRUE
# )
#
# get_prior(form, modeling_dat)
#
# priors <- c(
#   prior(normal(.2, .5), nlpar = "b"),
#   prior(normal(-0.2, 0.5), nlpar = "b", coef = "bl_q"),
#   prior(normal(5,5), class = "sd", nlpar = "shift")
# )
#
# fit <- brm(
#   form,
#   data= modeling_dat,
#   prior = priors,
#   chains = 5,
#   cores = parallel::detectCores() - 1,
#   iter = 2000,
#   # algorithm = "sampling",
#   algorithm = "meanfield",
#   backend = "cmdstanr",
#   silent = 2
# )
# saveRDS(fit,sprintf("%sresults/blq_slope_linear_pseudoprog_time_shift_brm_wpriors.rds", wd_path))
# blq_linear_fit <- readRDS(#   sprintf("%sresults/blq_slope_linear_pseudoprog_time_shift_brm_wpriors.rds", wd_path))


# EVAL MODELS ----
## compare models 
compare_mods <- function(mod1, mod2) {
  mod1_name <- deparse(substitute(mod1))
  mod2_name <- deparse(substitute(mod2))
  
  loo_mod1 <- loo(mod1, moment_match = FALSE)
  loo_mod2 <- loo(mod2, moment_match = FALSE)
  
  comp <- loo_compare(loo_mod1, loo_mod2)
  rownames(comp) <- c(mod1_name, mod2_name)
  return(comp)
}

# '>' expected log probability density and '<' standard error
# compare_mods(glbl_shft_intrcpt_logorithmic_fit, new_priors_glbl_shft_intrcpt_logorithmic_fit)

#                              elpd_diff   se_diff
# glbl_shft_intrcpt_logorithmic_fit  0.0       0.0   
# logorithmic_fit                   -2.6       3.2   
#                                          elpd_diff   se_diff
# glbl_shft_intrcpt_logorithmic_fit              0.0       0.0  
# new_priors_glbl_shft_intrcpt_logorithmic_fit -16.0       3.9  

# assess how much prior biasses 

# Posterior draws
post1 <- as_draws_df(glbl_shft_intrcpt_logorithmic_fit)
# Prior-only draws
post2<- as_draws_df(new_priors_glbl_shft_intrcpt_logorithmic_fit)

list(mcmc_areas(cbind(post = post1, new=post2),regex_pars = "b_shift_"),
mcmc_areas(cbind(post = post1, new=post2),regex_pars = "b_b_Int"),
mcmc_areas(cbind(post = post1, new=post2),regex_pars = "b_c_Int"),
mcmc_areas(cbind(post = post1, new=post2),regex_pars = "sd_"),
mcmc_areas(cbind(post = post1, new=post2),regex_pars = "b_a_Int")) %>% 
  wrap_plots(ncol = 3)


# which fit to use for estimates?
fit <- glbl_shft_intrcpt_logorithmic_fit

# plot(fit)
# pp_check(fit)
# conditional_effects(fit)


# PLOT SHIFTS----
modeling_dat <- raw_modeling_modeling_dat

# Extract posterior summary
post_summ <- posterior_summary(fit)

# Detect and extract global shift intercept if present
global_shift <- if ("b_shift_Intercept" %in% rownames(post_summ)) {
  post_summ["b_shift_Intercept", "Estimate"]
} else {
  0
}

# Extract patient-level shift deviations
shift_summ <- post_summ[grep("r_patID_eye", rownames(post_summ)), , drop = FALSE] %>%
  as.data.frame()

# Add global intercept if applicable
shift_summ[, c("Estimate", "Q2.5", "Q97.5")] <- shift_summ[, c("Estimate", "Q2.5", "Q97.5")] + global_shift
shift_summ$patID_eye <- str_extract(rownames(shift_summ), "(?<=\\[)[^,]+")

# Join to modeling data
modeling_dat <- left_join(modeling_dat, shift_summ)
modeling_dat$shift_year <- modeling_dat$year + modeling_dat$Estimate

# Normalize time (to keep values positive)
abs_min_shift_year <- min(modeling_dat$shift_year, na.rm = TRUE)*-1
modeling_dat$shift_year_norm <- modeling_dat$shift_year + abs_min_shift_year

sce$shift_year_norm <- modeling_dat$shift_year_norm[match(sce$patID_eye, modeling_dat$patID_eye)]
sce$shift_year <- modeling_dat$shift_year[match(sce$patID_eye, modeling_dat$patID_eye)]


# # Save shift values
# select(modeling_dat, sample, patID_eye, shift_year, shift_year_norm) %>% 
#   write_csv(sprintf("%sprocessed_data/logorithmic_fit_eye_time_shifts.csv", wd_path))

# Rescale pseudoprog
rescale_pp <- function(x) {
  pseudoprog_scale_attr$`scaled:scale` * x + pseudoprog_scale_attr$`scaled:center`
}
modeling_dat$rescaled_pseudoprog <- rescale_pp(modeling_dat$pseudoprog)

# Prepare plotting data
plot_dat <- modeling_dat %>% drop_na()
plot_dat_ci <- plot_dat %>% dplyr::slice(1, .by = patID_eye)
plot_dat_points <- plot_dat %>%
  group_by(patID_eye) %>%
  summarize(
    rescaled_pseudoprog = mean(rescaled_pseudoprog, na.rm = TRUE),
    shift_year_norm = mean(shift_year_norm, na.rm = TRUE),
    .groups = "drop"
  )

# Prediction grid
newdat <- data.frame(
  shift_year_norm = seq(min(plot_dat$shift_year_norm), max(plot_dat$shift_year_norm), length.out = 1000)
)

# Convert to raw year (remove global shift and normalization)
newdat$year <- newdat$shift_year_norm - abs_min_shift_year - global_shift
newdat$patID_eye <- NA  # disable random effects

# Model predictions (marginal)
preds <- fitted(fit, newdata = newdat, re_formula = NA)
newdat$pseudoprog_fit <- preds[, "Estimate"]
newdat$lower <- preds[, "Q2.5"]
newdat$upper <- preds[, "Q97.5"]

# Rescale predictions
newdat <- newdat %>%
  mutate(across(matches("fit|low|up"), rescale_pp))

# Fit linear model to inferred line for text annotation
infered_line <- lm(pseudoprog_fit ~ shift_year_norm, newdat)

# Final plot
ggplot() +
  
  geom_segment(
    data = plot_dat_ci,
    aes(x = Q2.5 + abs_min_shift_year,
        xend = Q97.5 + abs_min_shift_year,
        y = rescaled_pseudoprog,
        yend = rescaled_pseudoprog,
        color = Estimate),
    alpha = 0.3, linewidth = 2
  ) +
  # scale_color_gradient2(low = "blue", high = "red") +
  scale_color_viridis_c(alpha = 1)+
  labs(color = "Shift") +
  
  geom_point(
    data = plot_dat_ci,
    aes(x = Estimate + abs_min_shift_year,
        y = rescaled_pseudoprog),
    shape = 15, alpha = 0.3
  ) +
  geom_line(
    data = plot_dat,
    aes(shift_year_norm, rescaled_pseudoprog, group = patID_eye),
    alpha = 0.3, color = "black", linewidth = 0.3
  ) +
  
  geom_ribbon(
    data = newdat,
    aes(x = shift_year_norm, ymin = lower, ymax = upper),
    alpha = 0.6, fill = hcl(240, 100, 80)
  ) +
  geom_line(
    data = newdat,
    aes(x = shift_year_norm, y = pseudoprog_fit),
    linewidth = 1.2, color = hcl(240,80, 50)
  ) +
  
  xlab("Disease age (years)") +
  ylab("Pseudoprogression") +
  ggtitle("Shifted Visit Time\nPoints show BL shifted time with CI",
          as.character(fit$formula)[1]) +
  theme_bw()
ggsave(sprintf("%sresults/clustering/14.1_progressionTimeModeling/patient_time_shift.pdf",wd_path), height = 8, width = 12)

# copy to pub plots
sprintf("%sresults/clustering/14.1_progressionTimeModeling/patient_time_shift.pdf",wd_path) %>% 
  sprintf("cp %s %s", ., sprintf("%sresults/publication_figures/%s",wd_path, basename(.))) %>% 
  system()

# COMPARE SHIFTS BETWEEN MODELS ---
files <- list.files(sprintf("%sprocessed_data/", wd_path), full.names = T) %>% 
  grep("eye_time_shifts", ., value = TRUE)

eye_time_shifts <- lapply(files, read.csv)
names(eye_time_shifts) <- str_extract(files, "(?<=//).*") %>% 
  str_remove(., "\\.csv")

eye_time_shift_d <- sapply(eye_time_shifts, `[[`, "shift_year_norm") %>% 
  as.data.frame() 

with(eye_time_shift_d, plot(glbl_shft_intrcpt_logorithmic_fit_eye_time_shifts,new_priors_glbl_shft_intrcpt_logorithmic_fit_eye_time_shifts)); add_lines()

posterior_summary(glbl_shft_intrcpt_logorithmic_fit)

# PREDICT TIME TO EACH CLUSTER ----
cluster_pp_means <- data.frame(colData(sce)) %>%
  select(pseudoprog, my_cluster) %>%
  summarize(mean_pseudoprog = mean(pseudoprog), .by = my_cluster)


# between PPs how much time passes?

pdf(sprintf("%sresults/clustering/14.1_progressionTimeModeling/time_to_pp_from_0.pdf",wd_path), height = 5, width = 7)
  with(newdat, {
    plot(pseudoprog_fit~shift_year_norm, type = "l",
         main = "Years to Pseudoprogression", 
         xlab = "Pseudoprogression", ylab = "Years")
    
    polygon(
      y = c(lower, rev(upper)),  # bounds along x-axis
      x = c(shift_year_norm, rev(shift_year_norm)),
      col = rgb(0, 0, 1, 0.2),  # semi-transparent blue
      border = NA
    )
  })
  # add_lines("xy")
  walk(seq(0,50,5), function(x) abline(v = x, col = rgb(0,0, 0,.2)))
  walk(seq(0,50,5), function(x) abline(h = x, col = rgb(0,0, 0,.2)))
dev.off()



#+ use newdata as proxy of reciprocal of line fit
# FUNCTION to estimate year from pp along logistic line
calc_year <- function(pp){
  sapply(pp, function(x){
    # newdat$shift_year_norm[which.min(abs(x-newdat$pseudoprog_fit))]
    approx(x = newdat$pseudoprog_fit, y = newdat$shift_year_norm, xout = pp)$y
  })
}

# FUNCTION to estimate year difference between two pps
calc_year_diff <- function(from_pp, to_pp){
  if(length(from_pp)==length(to_pp)){
    sapply(to_pp, calc_year)-  sapply(from_pp, calc_year)
  }else{
    warning("vectors not equal length")
  }
}

from_pp <- 1:9
to_pp <- 10:18
calc_year_diff(from_pp , to_pp)

newdat %>% 
  ggplot(aes(pseudoprog_fit,shift_year_norm))+ 
  geom_line()+
  ylim(c(0,6))+xlim(c(0,35))+scale_y_continuous(breaks = seq(0,6,1))+scale_x_continuous(breaks = seq(0,35,2))+
  labs(color = "From")+ylab("To")+xlab("Estimated Year Difference")+ggtitle("Estimated time from one Chew to another")+
  theme_bw()

pdf(sprintf("%sresults/clustering/14.1_progressionTimeModeling/time_to_pp_from_0.pdf",wd_path), height = 5, width = 7)
plot(newdat$shift_year_norm ~ newdat$pseudoprog_fit,
     type = "l",
     xlim = c(0,30), ylim = c(0,55),
     main = "Years To Each Pseudoprogression", xlab = "Pseudoprogression", ylab = "Years")
add_lines("xy",lty = 2)
walk(seq(-10,30,5), function(x) abline(v = x, col = rgb(0,0, 0,.1)))
walk(seq(-10,50,5), function(x) abline(h = x, col = rgb(0,0, 0,.1)))
dev.off()

# copy to publication figures
system(sprintf("cp '%sresults/clustering/14.1_progressionTimeModeling/time_to_pp_from_0.pdf' '%sresults/publication_figures/time_to_pp_from_0.pdf'",wd_path,wd_path))
system(sprintf("ls %sresults/publication_figures/",wd_path))


from_to <- matrix(c(3,6,6,4,4,9,9,2,2,7,2,13,4,1,1,5,1,10,1,11,1,8 ), byrow = TRUE, ncol = 2)

from_to_d <- as.data.frame(from_to) %>%
  rename_with(~c("from", "to")) %>%
  mutate_all(factor) %>%
  left_join(cluster_pp_means %>% dplyr::rename(from=my_cluster,from_pp = mean_pseudoprog)) %>%
  left_join(cluster_pp_means%>% dplyr::rename(to =my_cluster, to_pp = mean_pseudoprog)) %>%
  mutate(nonlinear_year_diff = calc_year_diff(from_pp,to_pp))

infered_line %>%
  summary()
# Call:
#   lm(formula = pseudoprog_fit ~ shift_year_norm, data = newdat)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.1351 -0.8345  0.3221  1.0204  1.2553 
# 
# Coefficients:
#             Estimate Std.    Error t   value   Pr(>|t|)    
# (Intercept)     -0.137932   0.072200   -1.91   0.0564 .  
# shift_year_norm  0.576775   0.002347  245.75   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.142 on 998 degrees of freedom
# Multiple R-squared:  0.9837,	Adjusted R-squared:  0.9837 
# F-statistic: 6.039e+04 on 1 and 998 DF,  p-value: < 2.2e-16
#^^ ~0.6 pp per year

pp_p_yr <- coef(infered_line)[2]
yr_p_pp <- 1/pp_p_yr

# linear years are those estimated from inverse of lm line
cluster_pp_means$linear_yrs_to_pp <- cluster_pp_means$mean_pseudoprog * yr_p_pp

# calculate year difference between clusters using linear estimate ((1/slope) * pp)
from_to_d$linear_year_diff <- NA
for(i in seq(nrow(from_to_d))){
  # origin cluster
  from = from_to_d$from[i]
  # destination cluster
  to =  from_to_d$to[i]
  # find difference between source and dest clusters times
  from_years <- cluster_pp_means$linear_yrs_to_pp[cluster_pp_means$my_cluster==from]
  to_years <- cluster_pp_means$linear_yrs_to_pp[cluster_pp_means$my_cluster==to]
  year_diff <- to_years-from_years
  from_to_d$linear_year_diff[i] <- year_diff
  # from_to_d$linear_year_diff[i] = year_to_clust[
  #   rownames(year_to_clust) == as.character(from_to_d$from[i]) ,
  #   colnames(year_to_clust) == as.character(from_to_d$to[i])]
  #
}

# compare linear years to clust (using lm 1/slope) to non linear (using newdat of logistic fit)
pdf(sprintf("%sresults/clustering/14.1_progressionTimeModeling/nonLinearclust_clust_year_diff.pdf",wd_path),height = 9, width = 9, bg = "white")
with(from_to_d, plot(linear_year_diff, nonlinear_year_diff,type = "n", 
     main = "Linear versus non linear cluster-cluster year estimates", 
     xlab = "Linear Estimate (1/ß(pp/year))",
     ylab = "Non-linear Estimate (predict(fit))")
       )
add_lines()
with(from_to_d, text(linear_year_diff, nonlinear_year_diff, labels = paste(from,to, sep = "--")))
dev.off()

#^^ non-linear estimates logner times at later pseudotimes and less time at earlier

# how  long from from clust to to clust? ----
# fill matrix with differences between these clusters
year_to_clust <- matrix(NA, ncol = n_distinct(sce$my_cluster), nrow = n_distinct(sce$my_cluster),
                        dimnames = list(as.character(seq(n_distinct(sce$my_cluster))),as.character(seq(n_distinct(sce$my_cluster)))))


#+ difference in times between clusters of fromcluster and tocluster
for(i in 1:nrow(from_to)){
  row <- from_to[i,]

  # origin cluster
  from = row[[1]]
  # destination cluster
  to = row[[2]]

  # year_to_clust[from,to] <- from_to_d$linear_year_diff[from_to_d$from==from& from_to_d$to==to]
  year_to_clust[from,to] <- from_to_d$nonlinear_year_diff[from_to_d$from==from& from_to_d$to==to]
}


# plot
## filter out all na row/cols
year_to_clust_simplified <- year_to_clust[rowSums(!is.na(year_to_clust))!=0 , colSums(!is.na(year_to_clust))!=0]
## reorder
year_to_clust_simplified <- year_to_clust_simplified[order(rowSums(year_to_clust_simplified,na.rm = TRUE), decreasing = FALSE),order(colSums(year_to_clust_simplified,na.rm = TRUE), decreasing = FALSE)]

# stepwise
year_to_clust_hm <- Heatmap(year_to_clust_simplified,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = circlize::colorRamp2(seq(min(year_to_clust_simplified, na.rm = TRUE),
                                       max(year_to_clust_simplified, na.rm = TRUE), length.out = 10),
                                   plasma(10, begin = .4)),
        heatmap_legend_param = list(title = "Time"),
        show_heatmap_legend = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f",  year_to_clust_simplified[i, j]), x, y, gp = gpar(fontsize = 12))
        },
        row_names_side = "left",
        border = TRUE,
        column_title = "BRM Estimated years from cluster to cluster")

pdf(sprintf("%sresults/clustering/14.1_progressionTimeModeling/year_to_clust_hm.pdf",wd_path),  height = 6, width = 10)
draw(year_to_clust_hm)
dev.off()


# Calculate years from cluster 3
year_from_3 <-  matrix(NA, ncol = n_distinct(sce$my_cluster), nrow = n_distinct(sce$my_cluster),
                       dimnames = list(as.character(seq(n_distinct(sce$my_cluster))),as.character(seq(n_distinct(sce$my_cluster)))))

year_from_3[3,6] <- year_to_clust[3,6]

year_from_3[6,4] <- year_from_3[3,6]+year_to_clust[6,4]
year_from_3[4,1] <- year_from_3[6,4]+year_to_clust[4,1]
year_from_3[1,5] <- year_from_3[4,1]+year_to_clust[1,5]
year_from_3[1,8] <- year_from_3[4,1]+year_to_clust[1,8]
year_from_3[1,10] <- year_from_3[4,1]+year_to_clust[1,10]
year_from_3[1,11] <- year_from_3[4,1]+year_to_clust[1,11]

year_from_3[4,9] <- year_from_3[6,4]+year_to_clust[4,9]
year_from_3[9,2] <- year_from_3[4,9] +year_to_clust[9,2]
year_from_3[2,7] <- year_from_3[9,2]+year_to_clust[2,7]
year_from_3[2,13] <- year_from_3[9,2]+year_to_clust[2,13]


# plot
## filter out all na row/cols
year_from_3 <- year_from_3[ as.numeric(rownames(year_to_clust)) , as.numeric(colnames(year_to_clust))]

## filter out all na row/cols
year_from_3_simplified <- year_from_3[rowSums(!is.na(year_from_3))!=0 , colSums(!is.na(year_from_3))!=0]
## reorder
year_from_3_simplified <- year_from_3_simplified[order(rowSums(year_from_3_simplified,na.rm = TRUE), decreasing = FALSE),order(colSums(year_from_3_simplified,na.rm = TRUE), decreasing = FALSE)]

rownames(year_from_3_simplified) <- rep("3", nrow(year_from_3_simplified))

# cummulative
year_from_3_hm <- Heatmap(year_from_3_simplified,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = circlize::colorRamp2(seq(min(year_from_3_simplified, na.rm = TRUE),
                                       max(year_from_3_simplified, na.rm = TRUE), length.out = 10),
                                   plasma(10, begin = .4)),
        heatmap_legend_param = list(title = "Time"),
        show_heatmap_legend = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f",  year_from_3_simplified[i, j]), x, y, gp = gpar(fontsize = 12))
        },
        row_names_side = "left",
        border = TRUE,
        column_title = "BRM Estimated years to each cluster from 3")
pdf(sprintf("%sresults/clustering/14.1_progressionTimeModeling/year_from_3_hm.pdf",wd_path),  height = 6, width = 10)
draw(year_from_3_hm)
dev.off()




# plot time from three on tsne ----
tsne <- as.data.frame(reducedDim(sce, "TSNE")) %>% 
  rownames_to_column("sample") %>% 
  left_join(
    data.frame(colData(sce)) %>% select(sample, my_cluster)
  ) 

from_3_1d <-setNames(as.vector(year_from_3),paste0(row(year_from_3), "--", col(year_from_3)))
from_3_1d <- sort(from_3_1d[!is.na(from_3_1d)])

from_3_d <- as.data.frame(from_3_1d) %>%
  rownames_to_column("transition") %>% 
  separate(transition, c("from", "to"), "--") 

tree <- metadata(sce)$trajectory_data$tsne_line_data %>% 
  mutate(x = edge) %>% 
  separate(x, c("from", "to"), "--") 

labels <- tree %>% 
  arrange(TSNE1,TSNE2) %>% 
  distinct(across(c(TSNE1,TSNE2)), .keep_all = T)

labels$my_cluster <- factor(c(5,10,3,11,8,1,6,2,4,9,13,7,12))
# edge      TSNE1      TSNE2 from to my_cluster
# 1   5--1 -30.351792 -16.752198    5  1          5
# 2  10--1 -29.712359  -1.922073   10  1         10
# 3   6--3 -17.975657  21.497607    6  3          3
# 4  11--1 -17.359004  -7.207370   11  1         11
# 5   8--1  -9.066364 -31.040676    8  1          8
# 6   4--1  -8.555592  -9.062096    4  1          1
# 7   6--3  -1.895921  15.512989    6  3          6
# 8   7--2  10.997100 -22.604083    7  2          2
# 9   4--1  11.998380   5.740602    4  1          4
# 10  9--2  15.691141  -7.191089    9  2          9
# 11 13--2  17.040655 -36.295400   13  2         13
# 12  7--2  24.064459 -28.294888    7  2          7
# 13 12--4  34.064702   6.634810   12  4         12

labels <- labels %>% 
  left_join(from_3_d %>% select(to, from_3_1d), by = join_by(my_cluster == to)) %>% 
  dplyr::rename(year_from_3=from_3_1d)

ggplot()+
  geom_point(data = tsne, aes(TSNE1, TSNE2, color = my_cluster), alpha = .4, show.legend = TRUE)+
  coord_equal()+
  geom_line(data = tree, aes(TSNE1, TSNE2, group = edge), alpha = .6, linewidth = 1)+
  # geom_label(data = labels, aes(TSNE1, TSNE2, label = my_cluster),  alpha = .5)+
  geom_label(data = labels, size = 5, show.legend = FALSE,
             aes(TSNE1, TSNE2, fill = my_cluster, label = round(year_from_3,1)), alpha = .9)+
  scale_color_manual(values =  cluster_cols)+
  scale_fill_manual(values =  cluster_cols,guide = "none")+
  ggtitle("Years to each cluster from cluster 3")+
  labs( color = "Cluster")+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  theme_bw()

ggsave(sprintf("%sresults/clustering/14.1_progressionTimeModeling/brm_TSNE_cluster_with_time_from_3.pdf",wd_path), height = 10, width = 11)
ggsave(sprintf("%sresults/publication_figures/brm_TSNE_cluster_with_time_from_3.pdf", wd_path), height = 10, width = 11)

tsne <- tsne %>% 
  mutate(my_cluster = as.character(my_cluster)) %>% 
  mutate(my_cluster = case_when(
    my_cluster == '3'~  '1',
    my_cluster == '6'~  '2',
    my_cluster == '4'~  '3',
    my_cluster == '9'~  '4',
    my_cluster == '1'~  '5',
    my_cluster == '2'~  '6',
    my_cluster == '12' ~ '7',
    my_cluster == '10' ~ '8',
    my_cluster == '11' ~ '9',
    my_cluster == '8'~  '10',
    my_cluster == '5'~  '11',
    my_cluster == '7'~  '12',
    my_cluster == '13' ~ '13',
    TRUE ~ NA
  ))

labels <- labels %>% 
  mutate(across(c(from,to,my_cluster), ~as.character(.))) %>% 
  mutate(across(c(from,to,my_cluster), 
                ~case_when(
                  . == '3'~  '1',
                  . == '6'~  '2',
                  . == '4'~  '3',
                  . == '9'~  '4',
                  . == '1'~  '5',
                  . == '2'~  '6',
                  . == '12' ~ '7',
                  . == '10' ~ '8',
                  . == '11' ~ '9',
                  . == '8'~  '10',
                  . == '5'~  '11',
                  . == '7'~  '12',
                  . == '13' ~ '13',
                  TRUE ~ NA
                )
                ))

ggplot()+
  geom_point(data = tsne, aes(TSNE1, TSNE2, color = my_cluster), alpha = .4, show.legend = TRUE)+
  coord_equal()+
  geom_line(data = tree, aes(TSNE1, TSNE2, group = edge), alpha = .6, linewidth = 1)+
  # geom_label(data = labels, aes(TSNE1, TSNE2, label = my_cluster),  alpha = .5)+
  geom_label(data = labels, size = 5, show.legend = FALSE,
             aes(TSNE1, TSNE2, fill = my_cluster, label = round(year_from_3,1)), alpha = .9)+
  scale_color_manual(values =  cluster_cols)+
  scale_fill_manual(values =  cluster_cols,guide = "none")+
  ggtitle("Years to each cluster from cluster 3")+
  labs( color = "Cluster")+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  theme_bw()

ggsave(sprintf("%sresults/clustering/14.1_progressionTimeModeling/new_clust_brm_TSNE_cluster_with_time_from_3.pdf",wd_path), height = 10, width = 11)
ggsave(sprintf("%sresults/publication_figures/new_clust_brm_TSNE_cluster_with_time_from_3.pdf", wd_path), height = 10, width = 11)






