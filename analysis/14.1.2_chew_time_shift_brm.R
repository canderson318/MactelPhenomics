
# SETUP ----
options(max.print = 100)
# Packages
pacman::p_load(tidyr,reshape2,purrr,tidyverse,ggplot2,gridExtra,ComplexHeatmap,viridis,zoo,mgcv,mice,ggrepel,htmltools,Metrics,psych,car,stats,lme4,lmerTest,ape,SingleCellExperiment,scater,dplyr,HDF5Array,parallel, brms,cmdstanr,txtplot,ggnewscale)


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


# PLOT CHEW OVER TIME----

data.frame(colData(sce)) %>% 
  select(patID_eye, year, chew = chew_grade_derived_mine_loose) %>% 
  arrange(patID_eye, year) %>% 
  mutate(chew = cummax(chew), .by = patID_eye) %>% 
  drop_na() %>% 
  ggplot(aes(year,  chew ))+
  geom_line(aes(group = patID_eye), alpha = .2)+
  scale_color_viridis_c()+
  ggtitle("Chew Grade Over Time", "Each line an eye")+
  xlab("Chew Grade")+
  ylab("Year")+
  theme_classic()

sprintf("%sresults/clustering/14.1_progressionTimeModeling/chew_over_time.pdf", wd_path) %>% 
  {ggsave(., height = 5, width = 7); sprintf("cp %s %s", .,sprintf("%sresults/publication_figures/%s",wd_path, basename(.)) ) %>% system()}


# MAKE modeling_datA FOR MODELING----
raw_modeling_modeling_dat <- as.data.frame(colData(sce)) %>% 
  select(sample,patID_eye, month, chew = chew_grade_derived_mine_loose,BL_eye_age,AGE) %>% 
  mutate(year = month/12) %>% 
  group_by(patID_eye) %>% 
  arrange(year) %>% 
  filter(n()>1) %>%
  mutate(chew = cummax(chew)) %>% 
  ungroup() 

chew_scale <- scale(raw_modeling_modeling_dat$chew)
chew_scale_attr <- attributes(chew_scale)
raw_modeling_modeling_dat$chew <- as.numeric(chew_scale)

# bl_chew <- raw_modeling_modeling_dat %>% group_by(patID_eye) %>% dplyr::slice(1) 
# 
# bl_qs <- quantile(bl_chew$chew  ,  seq(0,1,.3),   na.rm = TRUE)
# # bl_qs <- seq(min(bl_chew$chew, na.rm = TRUE), max(bl_chew$chew, na.rm = TRUE), length.out = 6)
# bl_chew$bl_q <- cut(bl_chew$chew, breaks = bl_qs,labels = FALSE, include.lowest = TRUE)
# 
# raw_modeling_modeling_dat <- left_join(raw_modeling_modeling_dat, bl_chew %>% select(patID_eye,bl_q), by = 'patID_eye')
# 
# raw_modeling_modeling_dat$bl_q <- raw_modeling_modeling_dat$bl_q-1
# # raw_modeling_modeling_dat$bl_q <- factor(raw_modeling_modeling_dat$bl_q,levels = 0:n_distinct(raw_modeling_modeling_dat$bl_q))

raw_modeling_modeling_dat <- raw_modeling_modeling_dat %>% 
  mutate(BL_chew = dplyr::first(chew), .by = patID_eye)


raw_modeling_modeling_dat %>%
  drop_na() %>%
  mutate(across(matches('chew'),
    ~((.*chew_scale_attr$`scaled:scale`) + chew_scale_attr$`scaled:center`) %>% round)
    ) %>% 
  ggplot(aes(month, chew)) +
  
  # First color scale: for individual lines
  geom_line(aes(group = patID_eye), alpha = .2) +
  scale_alpha(range = c(0.01,0.3))+
  # scale_color_viridis_c(option = "plasma", name = "BL Eye Age") +
  
  # Start new color scale
  new_scale_color() +
  
  # Second color scale: for group smooths
  geom_smooth(aes(group = BL_chew, color = factor(BL_chew)), method = "lm", se = TRUE,linewidth = 2) +
  geom_smooth(aes(group = 1), method = "loess", se = TRUE,linewidth = 2, color = "hotpink") +
  scale_color_viridis_d(option = "turbo")+
  ggtitle("Chew Grade Over Time") +
  xlab("Year") +
  ylab("Chew Grade") +
  theme_classic()


x <- raw_modeling_modeling_dat %>%
  drop_na() %>%
  mutate(across(matches('chew'),
  ~((.*chew_scale_attr$`scaled:scale`) + chew_scale_attr$`scaled:center`) %>% round)  ) %>% 
  select(patID_eye, year, chew, BL_chew)

yr_stat_fun <- function(x) {
  quantile(x, probs = seq(0,1,.25)) %>%
    {.[length(.)-1]} %>% # grab second to last quantile
    # {.[length(.)]} %>% # grab last quantile, the max
    # {.^(1.5)} %>%
    {}
}

stat_yrs <- aggregate(year ~ BL_chew,yr_stat_fun , data = x)
stat_yrs$year <- cumsum(stat_yrs$year)- min(cumsum(stat_yrs$year)) + .01

# stat_yrs <- x %>% 
#   group_by(BL_chew) %>% 
#   summarize(mean_max_year = mean(max(year))) %>% 
#   mutate(mean_max_year = cumsum(mean_max_year) %>% {.-min(.) + .01}) %>% 
#   dplyr::rename(year = mean_max_year)


x <- x %>% 
  left_join(stat_yrs %>% dplyr::rename(year_shift = year))

x$shifted_year <- x$year+x$year_shift

# fit <- lmer(chew~ log(shifted_year) + (1|patID_eye), data = x)
fit <- lm(chew~ sqrt(shifted_year), data = x)
(summ <- summary(fit))
# plot(fit)

std_err <- exp(summ$coefficients[2, 2])
line_d <- data.frame(shifted_year = seq(.001,max(x$shifted_year)+10,length.out = 500))
line_d$chew <- predict(fit, line_d, re.form = NA, type = "response")
line_d$lower <- line_d$chew-std_err
line_d$upper <- line_d$chew+std_err


x %>% 
  ggplot(aes(shifted_year, chew)) +
  geom_line(aes(group = patID_eye), alpha = .2) +
  geom_smooth(aes(group = BL_chew, color = factor(BL_chew)), 
              method = "lm", se = TRUE,linewidth = 2)+
  # geom_smooth(aes(group = 1),
  #             method = "lm", formula = y ~ I(log(x))
  #             )+
  geom_ribbon(data = line_d, aes(ymin = lower, ymax = upper),fill = hcl(245,100,60,.3))+
  geom_line(data = line_d, linewidth = 2, color = hcl(245,100,57,.8))+
  theme_classic()

yr_to_chew_d <- data.frame(chew = 0L:6L)
yr_to_chew_d$approx_yr <- approx(line_d$chew, line_d$shifted_year,yr_to_chew_d$chew)$y
yr_to_chew_d$approx_yr[is.na(yr_to_chew_d$approx_yr)] <- 0

ggplot(yr_to_chew_d, aes(approx_yr,chew))+ 
  geom_point( size = 4)+geom_line()+
  ylim(c(0,6))+xlim(c(0,35))+scale_y_continuous(breaks = seq(0,6,1))+scale_x_continuous(breaks = seq(0,100,2))+
  labs(color = "From")+ylab("Chew")+xlab("Estimated Year Difference")+ggtitle("Estimated time from one Chew to another")+
  theme_bw()



# MODEL ----

modeling_dat <- raw_modeling_modeling_dat

# estimate logarithmic growth with c ----
# form <- bf(
#   chew ~ a + b * log(c + (year + shift)),
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
# ); saveRDS(fit, sprintf("%sresults/logorithmic_chew_time_shift_brm_wpriors.rds", wd_path))
# logorithmic_fit <- readRDS(sprintf("%sresults/logorithmic_chew_time_shift_brm_wpriors.rds", wd_path))

# # estimatete logarithmic growth w/out c ----
# form <- bf(
#   chew ~ a + b * log( (year + shift)),
#   a+b ~ 1,
#   shift ~1+ (1|patID_eye),
#   nl = TRUE
# )
# 
# priors <- c(
#   prior(normal(0, 2), nlpar = "a"),
#   prior(normal(1, 1), nlpar = "b"),
#   # prior(normal(0, 5), class = "sd", nlpar = "shift")
#   prior(normal(0, 5), class = "sd", nlpar = "shift", lb = 0)
# )
# 
# priors
# form
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
# ); saveRDS(fit, sprintf("%sresults/no_c_logorithmic_chew_time_shift_brm_wpriors.rds", wd_path))

no_c_logorithmic_fit <- readRDS(sprintf("%sresults/no_c_logorithmic_chew_time_shift_brm_wpriors.rds", wd_path))

# # estimate logarithmic growth with c and global shift intercept ----
# form <- bf(
#   # keep c because allows for flexibility; shift in when chew starts increasing
#   chew ~ a + b * log(c + (year + shift)),
#   a+b+c ~ 1,
#   shift ~ 1+ (1|patID_eye), # include global intercept, each patient intercept becomes difference from global
#   nl = TRUE
# )
# 
# priors <- c(
#   prior(normal(0, 1), nlpar = "a"),
#   # prior(normal(-2.15223, 1), nlpar = "a"), # global curve at zero when chew min???
#   prior(normal(1, 1), nlpar = "b"),
#   prior(normal(1, 10), nlpar = "c", lb = 1),
#   prior(normal(5, 10), nlpar = "shift") # global mean shift
# )
# 
# priors
# form
# 
# fit <- brm(
#   form,
#   prior = priors,
#   data= modeling_dat ,
#   chains = 4,
#   cores = parallel::detectCores() - 1,
#   iter = 2000,
#   algorithm = "sampling",
#   backend = "cmdstanr",
#   silent = 2
# );
# saveRDS(fit, sprintf("%sresults/global_shift_intercept_logorithmic_chew_time_shift_brm_wpriors.rds", wd_path))
# stop()

glbl_shft_intrcpt_logorithmic_fit <- readRDS(sprintf("%sresults/global_shift_intercept_logorithmic_chew_time_shift_brm_wpriors.rds", wd_path))

# # estimate linear pattern----
# form <- bf(
#   chew ~ b * (year + shift), # each individual curve is function of specific shift and time
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
# saveRDS(linear_fit, sprintf("%sresults/linear_chew_time_shift_brm_wpriors.rds", wd_path))
# linear_fit <- readRDS(sprintf("%sresults/linear_chew_time_shift_brm_wpriors.rds", wd_path))


# # estimate linear pattern with quantile interaction ----
# form <- bf(
#   chew ~ b * (year + shift), # each individual curve is function of specific shift and time
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
# saveRDS(fit,sprintf("%sresults/blq_slope_linear_chew_time_shift_brm_wpriors.rds", wd_path))
# blq_linear_fit <- readRDS(#   sprintf("%sresults/blq_slope_linear_chew_time_shift_brm_wpriors.rds", wd_path))


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
# compare_mods(glbl_shft_intrcpt_logorithmic_fit, no_c_logorithmic_fit)
#                               elpd_diff   se_diff
# glbl_shft_intrcpt_logorithmic_fit   0.0       0.0  
# no_c_logorithmic_fit              -16.4       3.1  

fit <- glbl_shft_intrcpt_logorithmic_fit

# plot(fit)
pp_check(fit, ndraws = 10)
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

# Add global intercept if achewlicable
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


# Save shift values
# select(modeling_dat, sample, patID_eye, shift_year, shift_year_norm) %>%
#   write_csv(sprintf("%sprocessed_data/chew_logorithmic_fit_eye_time_shifts.csv", wd_path))

# Rescale chew
rescale_chew <- function(x) {
  chew_scale_attr$`scaled:scale` * x + chew_scale_attr$`scaled:center`
}
modeling_dat$rescaled_chew <- rescale_chew(modeling_dat$chew)

# Prepare plotting data
plot_dat <- modeling_dat %>% drop_na()
plot_dat_ci <- plot_dat %>% dplyr::slice(1, .by = patID_eye)
plot_dat_points <- plot_dat %>%
  group_by(patID_eye) %>%
  summarize(
    rescaled_chew = mean(rescaled_chew, na.rm = TRUE),
    shift_year_norm = mean(shift_year_norm, na.rm = TRUE),
    .groups = "drop"
  )

# Prediction grid
newdat <- data.frame(
  shift_year_norm = seq(min(plot_dat$shift_year_norm), max(plot_dat$shift_year_norm), length.out = 500)
)

# Convert to raw year (remove global shift and normalization)
newdat$year <- newdat$shift_year_norm - abs_min_shift_year - global_shift
newdat$patID_eye <- NA  # disable random effects

# Model predictions (marginal)
preds <- fitted(fit, newdata = newdat, re_formula = NA)
newdat$chew_fit <- preds[, "Estimate"]
newdat$lower <- preds[, "Q2.5"]
newdat$upper <- preds[, "Q97.5"]

# Rescale predictions
newdat <- newdat %>%
  mutate(across(matches("fit|low|up"), rescale_chew))

# Fit linear model to inferred line for text annotation
infered_line <- lm(chew_fit ~ shift_year_norm, newdat)

# Final plot
ggplot() +
  
  geom_segment(
    data = plot_dat_ci,
    aes(x = Q2.5 + abs_min_shift_year,
        xend = Q97.5 + abs_min_shift_year,
        y = rescaled_chew,
        yend = rescaled_chew,
        color = Estimate),
    alpha = 0.1, linewidth = 2
  ) +
  # scale_color_gradient2(low = "blue", high = "red") +
  scale_color_viridis_c(alpha = 1)+
  labs(color = "Shift") +
  
  geom_point(
    data = plot_dat_ci,
    aes(x = Estimate + abs_min_shift_year,
        y = rescaled_chew),
    shape = 15, alpha = 0.1
  ) +
  geom_line(
    data = plot_dat,
    aes(shift_year_norm, rescaled_chew, group = patID_eye),
    alpha = 0.1, color = "black", linewidth = 0.3
  ) +
  
  geom_ribbon(
    data = newdat,
    aes(x = shift_year_norm, ymin = lower, ymax = upper),
    alpha = 0.6, fill = hcl(240, 100, 80)
  ) +
  geom_line(
    data = newdat,
    aes(x = shift_year_norm, y = chew_fit),
    linewidth = 1.2, color = hcl(240,80, 50)
  ) +
  
  xlab("Disease age (years)") +
  ylab("Chew Grade") +
  ggtitle("Shifted Visit Time\nPoints show BL shifted time with CI",
          as.character(fit$formula)[1]) +
  theme_bw()

ggsave(sprintf("%sresults/clustering/14.1_progressionTimeModeling/chew_patient_time_shift.pdf",wd_path), height = 8, width = 12)

sprintf("%sresults/clustering/14.1_progressionTimeModeling/chew_patient_time_shift.pdf",wd_path) %>% 
  {ggsave(., height = 8, width = 12); sprintf("cp %s %s", .,sprintf("%sresults/publication_figures/%s",wd_path, basename(.)) ) %>% system()}

# # COMPARE SHIFTS BETWEEN MODELS ---
files <- list.files(sprintf("%sprocessed_data/", wd_path), full.names = T) %>%
  grep("eye_time_shifts", ., value = TRUE) %>% 
  grep("chew", ., value = TRUE)

eye_time_shifts <- lapply(files, read.csv)
names(eye_time_shifts) <- str_extract(files, "(?<=//).*") %>%
  str_remove(., "\\.csv")

eye_time_shift_d <- sapply(eye_time_shifts, `[[`, "shift_year_norm") %>%
  as.data.frame()

plot(eye_time_shift_d);add_lines()

# posterior_summary(glbl_shft_intrcpt_logorithmic_fit)

pp <- read.csv(sprintf("%sprocessed_data/%s", wd_path,"glbl_shft_intrcpt_logorithmic_fit_eye_time_shifts.csv"))
chew <- read.csv(sprintf("%sprocessed_data/%s", wd_path,"chew_logorithmic_fit_eye_time_shifts.csv"))

shift_time <- left_join(pp, chew, by = c("sample", 'patID_eye'), suffix = c(".pp", ".chew"))
pptime <- shift_time$shift_year_norm.pp
chewtime <- shift_time$shift_year_norm.chew


# Plot the scatter ----
# Fit linear model
mod <- lm(pptime ~ chewtime)
summ <- summary(mod)
coefs <- coef(summ)

lind <- data.frame(chewtime = seq(min(chewtime), max(chewtime), length.out = 200))
# Get predictions and confidence intervals
lind <-cbind(lind, predict(mod, newdata = lind, interval = "confidence"))
  
ggplot(shift_time, aes(shift_year_norm.chew, shift_year_norm.pp))+
  geom_hline(yintercept = 0, linetype = 2, col = "grey")+
  geom_vline(xintercept = 0, linetype = 2, col = "grey")+
  geom_abline(intercept = 0,slope = 1, linetype = 2, col = "grey")+
  geom_line(aes(group = patID_eye),show.legend = FALSE, alpha = .1)+
  geom_ribbon(data = lind, aes(x =chewtime,y = fit, ymin = lwr, ymax = upr), 
              fill = hcl(240,100,100), color = NA)+
  geom_abline(slope = coefs[2,1], intercept = coefs[1,1], 
              color = hcl(240,100,60), linewidth = 1)+
  geom_smooth()+
  ggtitle("Time Shifted With Pseudoprog and Chew Grade",
          sprintf("pseudoprogtime ~ %.1f + %.1f*chewtime; (p = %.1e)",coefs[1,1],coefs[2,1], coefs[2,4]))+
  xlab("Chew Shifted Time (years)")+
  ylab("Pseudoprogression Shifted Time (years)")+
  coord_equal()+
  theme_classic()


# PREDICT TIME TO EACH CLUSTER ----
cluster_chew_means <- data.frame(colData(sce)) %>%
  select(chew = chew_grade_derived_mine_loose, my_cluster) %>%
  summarize(mean_chew = mean(chew), .by = my_cluster)


# between chews how much time passes?
#+ use line data as proxy of the reciprocal of line fit

# ****** FUNCTION to estimate year from chew along logistic line ******
#+ using the line fitted to shifted time plot, estimate year shift for chew input
#+ return the  shifted time where chew closest to fitted chew value
calc_year <- function(chew){
  sapply(chew, function(x){
    # newdat$shift_year_norm[which.min(abs(x-newdat$chew_fit))]
    approx(x = newdat$chew_fit, y = newdat$shift_year_norm, xout = chew)$y
  })
}


# ****** FUNCTION to estimate year difference between two chews
#+ find difference in year estimate between to clusteres
calc_year_diff <- function(from_chew, to_chew){
  if(length(from_chew)==length(to_chew)){
    sapply(to_chew, calc_year)-  sapply(from_chew, calc_year)
  }else{
    warning("vectors not equal length")
  }
}

z <- expand.grid(from = 0:6, to = 0:6, KEEP.OUT.ATTRS = F)
z <- z[z[,1]<z[,2] , ]
b <- cbind(z, year_diff = calc_year_diff(z[,1],z[,2])) %>% 
  data.frame() 

# plot all combinations of from chew to chew year differences
ggplot(b, aes(year_diff,to))+ 
  geom_point(aes(color = factor(from)), size = 4)+geom_line(aes(color = factor(from)), size = 1)+
  ylim(c(0,6))+xlim(c(0,35))+scale_y_continuous(breaks = seq(0,6,1))+scale_x_continuous(breaks = seq(0,35,2))+
  labs(color = "From")+ylab("To")+xlab("Estimated Year Difference")+ggtitle("Estimated time from one Chew to another")+
  theme_bw()

b<- rbind(data.frame(from =  0, to = 0, year_diff = 0), b)
ggplot(b %>% filter(from == 0), aes(year_diff,to))+ 
  geom_point( size = 4)+geom_line()+
  ylim(c(0,6))+xlim(c(0,35))+scale_y_continuous(breaks = seq(0,6,1))+scale_x_continuous(breaks = seq(0,35,2))+
  labs(color = "From")+ylab("Chew")+xlab("Estimated Year Difference")+ggtitle("Estimated time from one Chew to another")+
  theme_bw()

ggsave(sprintf("%sresults/clustering/14.1_progressionTimeModeling/time_to_chew_from_0.pdf",wd_path), width = 6, height =5)


z <- data.frame(brm = b$year_diff[b$from==0],man= yr_to_chew_d$approx_yr[])

{
  plot(z,
       xlim = c(0,max(as.matrix(z))),
       ylim = c(0,max(as.matrix(z))),
       # asp = 1,
       type ="n",
       main = "Comparison between BRMS approx. and manually shifted year to chew" ,
       xlab = "BRMS estimated year to chew from 0",
       ylab = "Manual estimated year to chew from 0")
  add_lines()
  points(z, 
       type = "c")
  text(z, labels = 0:6)
  legend("topleft", 
         legend = sprintf("RMSE = %0.2f; R^2 = %0.4f",
                          rmse(z[,1], z[,2]),
                          cor(z[,1], z[,2])
                          )
         )
}

stop("STOPPED")

# yrs_from_0 <- b[b$from==0,]
# pdf(sprintf("%sresults/clustering/14.1_progressionTimeModeling/time_to_chew_from_0.pdf",wd_path), height = 5, width = 7)
# plot(yrs_from_0$to ,yrs_from_0$year_diff,
#      type = "b",#xlim = c(0,6), ylim = c(0,30),
#      main = "Years To Each Chew Grade", xlab = "Chew Grade", ylab = "Years")
# # add_lines("xy")
# walk(seq(0,6,1), function(x) abline(v = x, col = rgb(0,0, 0,.2)))
# walk(seq(0,30,5), function(x) abline(h = x, col = rgb(0,0, 0,.2)))
# dev.off()

system(sprintf("ls %sresults/",wd_path))
system(sprintf("ls %sresults/clustering/14.1_progressionTimeModeling/",wd_path))
system(sprintf("cp '%sresults/clustering/14.1_progressionTimeModeling/time_to_chew_from_0.pdf' '%sresults/publication_figures/time_to_chew_from_0.pdf'",wd_path,wd_path))
system(sprintf("ls %sresults/publication_figures/",wd_path))
  
  # from_to <- matrix(c(3,6,6,4,4,9,9,2,2,7,2,13,4,1,1,5,1,10,1,11,1,8 ), byrow = TRUE, ncol = 2)
# 
# from_to_d <- as.data.frame(from_to) %>%
#   rename_with(~c("from", "to")) %>%
#   mutate_all(factor) %>%
#   left_join(cluster_chew_means %>%select(my_cluster, mean_chew) %>%  rename(from=my_cluster,from_chew = mean_chew)) %>%
#   left_join(cluster_chew_means%>%select(my_cluster, mean_chew) %>%  rename(to =my_cluster, to_chew = mean_chew)) %>%
#   mutate(nonlinear_year_diff = calc_year_diff(from_chew,to_chew))
# 
# infered_line %>%
#   summary()
# # Call:
# #   lm(formula = chew_fit ~ shift_year_norm, data = newdat)
# # 
# # Residuals:
# #   Min      1Q  Median      3Q     Max 
# # -1.4209 -0.2929  0.1252  0.3815  0.4678 
# # 
# # Coefficients:
# #   Estimate Std. Error t value Pr(>|t|)    
# # (Intercept)     -0.09348    0.03931  -2.378   0.0178 *  
# #   shift_year_norm  0.20235    0.00201 100.668   <2e-16 ***
# #   ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# # 
# # Residual standard error: 0.4401 on 498 degrees of freedom
# # Multiple R-squared:  0.9532,	Adjusted R-squared:  0.9531 
# # F-statistic: 1.013e+04 on 1 and 498 DF,  p-value: < 2.2e-16
# #^^ ~0.2 chew per year
# 
# chew_p_yr <- coef(infered_line)[2]
# yr_p_chew <- 1/chew_p_yr
# 
# cluster_chew_means$linear_yrs_to_chew <- cluster_chew_means$mean_chew * yr_p_chew
# cluster_chew_means %>% ggplot(aes(mean_chew, linear_yrs_to_chew, color = my_cluster))+geom_point()+scale_color_manual(values = cluster_cols)+scale_y_continuous(breaks = seq(0,30,2))+ggtitle("years to cluster average chew")
# 
# # calculate year difference between clusters using linear estimate ((1/slope) * chew)
# from_to_d$linear_year_diff <- NA
# for(i in seq(nrow(from_to_d))){
#   # origin cluster
#   from = from_to_d$from[i]
#   # destination cluster
#   to =  from_to_d$to[i]
#   # find difference between source and dest clusters times
#   from_years <- cluster_chew_means$linear_yrs_to_chew[cluster_chew_means$my_cluster==from]
#   to_years <- cluster_chew_means$linear_yrs_to_chew[cluster_chew_means$my_cluster==to]
#   year_diff <- to_years-from_years
#   from_to_d$linear_year_diff[i] <- year_diff
# }
# 
# # compare linear years to clust (using lm 1/slope) to non linear (using newdat of logistic fit)
# pdf(sprintf("%sresults/clustering/14.1_progressionTimeModeling/chew_nonLinearclust_clust_year_diff.pdf",wd_path),height = 9, width = 9, bg = "white")
# with(from_to_d, plot(linear_year_diff, nonlinear_year_diff,type = "n", 
#                      main = "Linear versus non linear cluster-cluster year estimates", 
#                      xlab = "Linear Estimate (1/ß(chew/year))",
#                      ylab = "Non-linear Estimate (predict(fit))")
# )
# add_lines()
# with(from_to_d, text(linear_year_diff, nonlinear_year_diff, labels = paste(from,to, sep = "--")))
# dev.off()
# 
# 
# # how  long from from clust to to clust? ----
# # fill matrix with differences between these clusters
# year_to_clust <- matrix(NA, ncol = n_distinct(sce$my_cluster), nrow = n_distinct(sce$my_cluster),
#                         dimnames = list(as.character(seq(n_distinct(sce$my_cluster))),as.character(seq(n_distinct(sce$my_cluster)))))
# 
# 
# #+ difference in times between clusters of fromcluster and tocluster
# for(i in 1:nrow(from_to)){
#   row <- from_to[i,]
#   # origin cluster
#   from = row[[1]]
#   # destination cluster
#   to = row[[2]]
#   # year_to_clust[from,to] <- from_to_d$linear_year_diff[from_to_d$from==from& from_to_d$to==to]
#   year_to_clust[from,to] <- from_to_d$nonlinear_year_diff[from_to_d$from==from& from_to_d$to==to]
# }
# 
# # plot
# ## filter out all na row/cols
# year_to_clust_simplified <- year_to_clust[rowSums(!is.na(year_to_clust))!=0 , colSums(!is.na(year_to_clust))!=0]
# ## reorder
# year_to_clust_simplified <- year_to_clust_simplified[order(rowSums(year_to_clust_simplified,na.rm = TRUE), decreasing = FALSE),order(colSums(year_to_clust_simplified,na.rm = TRUE), decreasing = FALSE)]
# 
# # stepwise
# year_to_clust_hm <- Heatmap(year_to_clust_simplified,
#                             cluster_rows = FALSE,
#                             cluster_columns = FALSE,
#                             col = circlize::colorRamp2(seq(min(year_to_clust_simplified, na.rm = TRUE),
#                                                            max(year_to_clust_simplified, na.rm = TRUE), length.out = 10),
#                                                        plasma(10, begin = .4)),
#                             heatmap_legend_param = list(title = "Time"),
#                             show_heatmap_legend = FALSE,
#                             cell_fun = function(j, i, x, y, width, height, fill) {
#                               grid.text(sprintf("%.1f",  year_to_clust_simplified[i, j]), x, y, gp = gpar(fontsize = 12))
#                             },
#                             row_names_side = "left",
#                             border = TRUE,
#                             column_title = "BRM Estimated years from cluster to cluster")
# 
# pdf(sprintf("%sresults/clustering/14.1_progressionTimeModeling/chew_year_to_clust_hm.pdf",wd_path),  height = 6, width = 10)
# draw(year_to_clust_hm)
# dev.off()
# 
# 
# # Calculate years from cluster 3
# year_from_3 <-  matrix(NA, ncol = n_distinct(sce$my_cluster), nrow = n_distinct(sce$my_cluster),
#                        dimnames = list(as.character(seq(n_distinct(sce$my_cluster))),as.character(seq(n_distinct(sce$my_cluster)))))
# 
# year_from_3[3,6] <- year_to_clust[3,6]
# 
# year_from_3[6,4] <- year_from_3[3,6]+year_to_clust[6,4]
# year_from_3[4,1] <- year_from_3[6,4]+year_to_clust[4,1]
# year_from_3[1,5] <- year_from_3[4,1]+year_to_clust[1,5]
# year_from_3[1,8] <- year_from_3[4,1]+year_to_clust[1,8]
# year_from_3[1,10] <- year_from_3[4,1]+year_to_clust[1,10]
# year_from_3[1,11] <- year_from_3[4,1]+year_to_clust[1,11]
# 
# year_from_3[4,9] <- year_from_3[6,4]+year_to_clust[4,9]
# year_from_3[9,2] <- year_from_3[4,9] +year_to_clust[9,2]
# year_from_3[2,7] <- year_from_3[9,2]+year_to_clust[2,7]
# year_from_3[2,13] <- year_from_3[9,2]+year_to_clust[2,13]
# 
# 
# # plot
# ## filter out all na row/cols
# year_from_3 <- year_from_3[ as.numeric(rownames(year_to_clust)) , as.numeric(colnames(year_to_clust))]
# 
# ## filter out all na row/cols
# year_from_3_simplified <- year_from_3[rowSums(!is.na(year_from_3))!=0 , colSums(!is.na(year_from_3))!=0]
# ## reorder
# year_from_3_simplified <- year_from_3_simplified[order(rowSums(year_from_3_simplified,na.rm = TRUE), decreasing = FALSE),order(colSums(year_from_3_simplified,na.rm = TRUE), decreasing = FALSE)]
# 
# rownames(year_from_3_simplified) <- rep("3", nrow(year_from_3_simplified))
# 
# # cummulative
# year_from_3_hm <- Heatmap(year_from_3_simplified,
#                           cluster_rows = FALSE,
#                           cluster_columns = FALSE,
#                           col = circlize::colorRamp2(seq(min(year_from_3_simplified, na.rm = TRUE),
#                                                          max(year_from_3_simplified, na.rm = TRUE), length.out = 10),
#                                                      plasma(10, begin = .4)),
#                           heatmap_legend_param = list(title = "Time"),
#                           show_heatmap_legend = FALSE,
#                           cell_fun = function(j, i, x, y, width, height, fill) {
#                             grid.text(sprintf("%.1f",  year_from_3_simplified[i, j]), x, y, gp = gpar(fontsize = 12))
#                           },
#                           row_names_side = "left",
#                           border = TRUE,
#                           column_title = "BRM Estimated years to each cluster from 3")
# pdf(sprintf("%sresults/clustering/14.1_progressionTimeModeling/chew_year_from_3_hm.pdf",wd_path),  height = 6, width = 10)
# draw(year_from_3_hm)
# dev.off()
# 
