
options(max.print = 100)
cat(rep("\n",5))
cat("timeAX.R\n")
cat(rep("\n",5))

# Packages
pacman::p_load(tidyr,reshape2,purrr,tidyverse,ggplot2,gridExtra,ComplexHeatmap,viridis,
               zoo,mgcv,ggrepel,Metrics,car,stats,SingleCellExperiment,scater,dplyr,HDF5Array,
               txtplot,TimeAx, brms)

rm(list = ls())
gc()

#my functions
source("/Users/canderson/Documents/LMRI/PRJ2024001/version015/analysis/00_MyFunctions.R")


# LOAD modeling_datA ----
wd_path = "~/Documents/LMRI/PRJ2024001/version015/"

# results directory 
sprintf("%sresults/timeAX/", wd_path) %>% 
  {if(!dir.exists(.)) dir.create(.)}

# load experiment
sce <-  loadHDF5SummarizedExperiment(paste0(wd_path, "processed_data/sce006.h5"))

# load time shif data
time_shift <- read_csv(sprintf("%sprocessed_data/eye_time_shifts.csv",wd_path))

# PREP DAATA ----

# add time shift
all(sce$sample[match(time_shift$sample, sce$sample)] ==  time_shift$sample[na.omit(match(sce$sample,time_shift$sample))])
sce$shift_year_norm[na.omit(match(time_shift$sample, sce$sample))] <-   time_shift$shift_year_norm[na.omit(match(sce$sample,time_shift$sample))]
sce$shift_year[na.omit(match(time_shift$sample, sce$sample))] <-   time_shift$shift_year[na.omit(match(sce$sample,time_shift$sample))]



# mean and sd
scale <- rowSds(assays(sce)$measure_vars , na.rm = TRUE)
center <- rowMeans(assays(sce)$measure_vars , na.rm = TRUE)

# filter for >=3 visits
eyes_gr_3_vis <- (table(sce$patID_eye)>=3) %>% 
  {names(.)[.]} 
sub_sce <- sce[ , sce$patID_eye%in%eyes_gr_3_vis]

# filter out where no scale or mean
sub_sce <- sub_sce[!is.na(scale) & !is.na(center),]

# order so each patient has consecutive visits
sub_sce <- sub_sce[ , order(sub_sce$patID_eye, sub_sce$month)]

# txtplot(sub_sce$month)
# txtplot(sce$month)
# head(sub_sce$month,20)
# any(table(sub_sce$patID_eye)<3)

# grab assay
# assay <- assays(sub_sce)$pmm_imp_of_scaled %>%  
#   as.matrix() 

assay <-t(reducedDim(sub_sce, "PCA"))

colnames(assay) <- colnames(sub_sce)


# Make PCs positively correlated with diesease progression
corm <- cbind(t(assay), as.matrix(t(assays(sub_sce)$pmm_imp_of_scaled)) )%>% 
  cor(use = "pairwise.complete") 

# filter for pcs in columns
corm <- corm[-grep("PC\\d", rownames(corm)),grep("PC\\d", rownames(corm))]
Heatmap(t(corm))

# filter for features that correlate to progression
corm <- corm[grep("NEOVASC|IR|ISOS_BREAK|HR|RPE", rownames(corm)),]
Heatmap(t(corm))

# sign of correlation to disease
signs <- sign(colMeans(corm) )
  

# barplot(signs,las = 2)
# assay["PC1",] %>% density() %>% plot()
# sweep(assay, 1, signs,"*")["PC1",] %>% density() %>% lines(col = "blue")

# multiply rows of assay by respective sign
assay <-   sweep(assay, 1, signs,"*")


# Transpose to match TimeAx expectation: rows = features, cols = samples
# Recreate grouping
grouped <- split(as.data.frame(t(assay)), sub_sce$patID_eye)

# Check which features exist in ALL samples
#+ for each eye, check that all of their features have missingngess less than ncol
valid_features <- Reduce(intersect, lapply(grouped, function(x) {
  x <- t(x)
  rownames(x)[rowSums(is.na(x) | is.nan(x) | is.infinite(x)) < ncol(x)]
  })
)

setdiff_intersect(list(rownames(assay), valid_features), values = T)

stopifnot(length(sub_sce$patID_eye)==ncol(assay))

# set.seed(203)
# test_ids <-  colData(sub_sce)[ sub_sce$patID_eye %in% sample(sub_sce$patID_eye, 10, F) , c("patID_eye", "sample")] 
# test <- assay[, colnames(assay)%in% test_ids$sample]

# how.long(T)
# model = modelCreation(
#   # TEsT
#   # trainData = test, # data
#   # sampleNames = test_ids$patID_eye, # vector of ids for each sample, each id must have >3 obs
# 
#   trainData = assay, # data
#   sampleNames = sub_sce$patID_eye, # vector of ids for each sample, each id must have >3 obs
#   seed = valid_features, # what features to use in DTW
#   no_cores = 6,
#   ratio = FALSE
#   )
# how.long() # ~ 1hr
# saveRDS(model, sprintf("%sresults/timeAX/mod.rds", wd_path))
model <- readRDS(sprintf("%sresults/timeAX/mod.rds", wd_path))


# pseudotimeStats = predictByConsensus(model,assay,no_cores = 6)
# saveRDS(pseudotimeStats, sprintf("%sresults/timeAX/pt_stats.rds", wd_path))
pseudotimeStats <- readRDS(sprintf("%sresults/timeAX/pt_stats.rds", wd_path))


# robustnessStats = robustness(model,assay,sub_sce$patID_eye,no_cores = 6)
# saveRDS(robustnessStats, sprintf("%sresults/timeAX/robust_stats.rds", wd_path))
# robustnessStats <- readRDS( sprintf("%sresults/timeAX/robust_stats.rds", wd_path))


# Look at Model ----
trajs <- sapply(model$consensusList, `[[`,2)

trajs[[1]] %>% plot(type = "l", col = rgb(0, 0, 1, 0.1), xlim = c(0,2350), ylim = c(0,1.1))
walk(seq_along(trajs), ~ lines(trajs[[.x]], col = rgb(0, 0, 1, 0.1)))


pseudotime = pseudotimeStats$predictions
sub_sce$timeax <- pseudotime
uncertainty = pseudotimeStats$uncertainty

# robustnessPseudo = robustnessStats$robustnessPseudo
# robustnessScore = robustnessStats$score

# Fit GAM model
gamfit <- gam(uncertainty ~ s(pseudotime))

# Create new data for prediction
gamline <- data.frame(pseudotime = seq(min(pseudotime), max(pseudotime), length.out = 200))

# Predict with standard errors
preds <- predict(gamfit, newdata = gamline, se.fit = TRUE)
gamline$fit <- preds$fit
gamline$upper <- preds$fit + 2 * preds$se.fit
gamline$lower <- preds$fit - 2 * preds$se.fit

# Plot
plot(pseudotime, uncertainty, col = rgb(0, 0, 1, 0.2), pch = 16)
polygon(c(gamline$pseudotime, rev(gamline$pseudotime)),c(gamline$upper, rev(gamline$lower)),col = hcl(300,100,50,.4), border = NA)
lines(gamline$pseudotime, gamline$fit, col = hcl(300,100,45,1), lwd = 2)

# Prepare data
plot_data <- colData(sub_sce) %>%
  data.frame() %>%
  select(pseudoprog, timeax, year, shift_year_norm)


GGally::ggpairs(
  plot_data,
  lower = list(continuous = GGally::wrap("smooth", alpha = 0.1, color = rgb(0,0,0.4,.3))),
  upper = list(continuous = GGally::wrap("cor", size = 4)),  # default is cor anyway
  diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.5))
)+
  theme_bw()+
  theme(aspect.ratio = 1)

corm <- colData(sub_sce) %>%
  data.frame() %>%
  select(sample,pseudoprog, year, shift_year_norm, timeax, matches("RW"), matches("rs"), chew_grade_derived_mine_loose) %>%
  left_join(reducedDim(sce, "PCA") %>% as.data.frame() %>% rownames_to_column("sample")) %>%
  left_join(assays(sce)$pmm_imp_of_scaled %>% t() %>% as.data.frame() %>% rownames_to_column("sample")) %>%
  select(-sample) %>% 
  cor(use = "pairwise.complete",method = "pearson") 
  

corm <- corm[grep("timeax|pseudoprog|shift_year|norm",colnames(corm)),-grep("timeax|pseudoprog|shift_year|norm",colnames(corm))  ]

Heatmap(corm, 
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2"
        )


plotTSNE(sub_sce, color_by = 'timeax')+
  coord_equal()

dat <- colData(sub_sce) %>%
  data.frame() %>%
  left_join(assays(sce)$pmm_imp_of_scaled %>% t() %>% as.data.frame() %>% rownames_to_column("sample")) %>% 
  select(sample,patID_eye, year,chew_grade_derived_mine_loose,pseudoprog, shift_year_norm, timeax) %>% 
  group_by(patID_eye) %>% 
  arrange(year) %>% 
  # mutate(timeax = timeax - dplyr::first(timeax)) %>% 
  ungroup()

dat %>% 
  ggplot(aes(year, timeax))+
  geom_hline(yintercept = 0, linetype = 2, color = hcl(120,100, 80), linewidth = 1.5)+
  geom_smooth(aes(group = patID_eye), 
              show.legend = FALSE, linewidth = 1,
              se = FALSE,
              color = rgb(0,0,0,.1)
              )+
  geom_smooth()+
  theme_classic()


# fit <- brm(
#   bf(
#     timeax ~ b * log(a + year),
#     a + b ~ 1,
#     nl = TRUE
#   ),
#   family = gaussian(),
#   data = dat,
#   backend = "cmdstanr",
#   iter = 2000,
#   chains = 4,
#   cores = parallel::detectCores() - 1,
#   algorithm = "sampling",
#   silent = 2
# )
# 
# fit
# bayesplot::mcmc_areas(fit, regex_pars = "b_|sigma|sd")
# conditional_effects(fit)

qs <- data.frame(colData(sub_sce)) %>% 
  select(patID_eye, year, timeax) %>% 
  arrange(year) %>% 
  dplyr::slice(1,.by = patID_eye) %>% 
  mutate(quantile = cut(timeax, quantile(timeax,probs = seq(0,1,.2)), include.lowest = T, labels = FALSE) )

sub_sce$timeax_bl_q <- NA  
sub_sce$timeax_bl_q[match(sub_sce$patID_eye,qs$patID_eye )] <- qs$quantile[match(sub_sce$patID_eye ,qs$patID_eye)]

plotTSNE(sub_sce[,!is.na(sub_sce$timeax_bl_q)], color_by = "timeax_bl_q")+ coord_equal()


cbind(chew = sub_sce$chew_grade_derived_mine_loose, pp= sub_sce$pseudoprog, tmax = sub_sce$timeax, shiftyr = sub_sce$shift_year_norm) %>% 
  # cbind(., rowProds(.)) %>% 
  scale() %>% 
  prcomp() %>% 
  {print(.[["rotation"]]); print(.[["sdev"]] %>% {((.**2 / sum(.**2)) * 100) } ); .[["x"]]} %>% 
  as.data.frame() %>% 
  pairs(col = rgb(0,.3,.3, .1), pch = 16)

dat <- colData(sub_sce) %>%
  data.frame() %>%
  left_join(assays(sce)$pmm_imp_of_scaled %>% t() %>% as.data.frame() %>% rownames_to_column("sample"))


X <- select(dat, all_of(rownames(sce))) %>% 
  as.matrix()

timeax_mod <-  lm(RW064 ~ timeax, data = dat)
pp_mod <-      lm(RW064 ~ pseudoprog, data = dat)
chew_mod <-    lm(RW064 ~ chew_grade_derived_mine_loose, data = dat)
yrshift_mod <- lm(RW064 ~ shift_year_norm, data = dat)

compare_models <- function(mod1, mod2, mod1_name = "Model 1", mod2_name = "Model 2") {
  cat("=== ANOVA ===\n")
  print(anova(mod1, mod2))
  
  cat("\n=== BIC ===\n")
  bic <- BIC(mod1, mod2)
  rownames(bic) <- c(mod1_name, mod2_name)
  print(bic)
  
  cat("\n=== AIC ===\n")
  aic <- AIC(mod1, mod2)
  rownames(aic) <- c(mod1_name, mod2_name)
  print(aic)
  
  cat("\n=== Summary: ", mod1_name, " ===\n", sep = "")
  print(summary(mod1))
  
  cat("\n=== Summary: ", mod2_name, " ===\n", sep = "")
  print(summary(mod2))
}

compare_models(timeax_mod,pp_mod, 'timeax_mod','pp_mod')
# === ANOVA ===
#   Analysis of Variance Table
# 
# Model 1: RW064 ~ timeax
# Model 2: RW064 ~ pseudoprog
# Res.Df    RSS Df Sum of Sq F Pr(>F)
# 1   5173 6139.5                      
# 2   5173 4457.3  0    1682.2         
# 
# === BIC ===
#   df      BIC
# timeax_mod  3 15596.10
# pp_mod      3 13939.04
# 
# === AIC ===
#   df      AIC
# timeax_mod  3 15576.44
# pp_mod      3 13919.39
# 
# === Summary: timeax_mod ===
#   
#   Call:
#   lm(formula = RW064 ~ timeax, data = dat)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.31446  0.04107  0.62806  0.66936  0.71113 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept) -0.05457    0.03057  -1.785   0.0744 .
# timeax       0.11318    0.05575   2.030   0.0424 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.089 on 5173 degrees of freedom
# Multiple R-squared:  0.000796,	Adjusted R-squared:  0.0006029 
# F-statistic: 4.121 on 1 and 5173 DF,  p-value: 0.0424
# 
# 
# === Summary: pp_mod ===
#   
#   Call:
#   lm(formula = RW064 ~ pseudoprog, data = dat)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.6058 -0.2832  0.1849  0.6475  2.2618 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.939779   0.024863   37.80   <2e-16 ***
#   pseudoprog  -0.092014   0.002079  -44.25   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9282 on 5173 degrees of freedom
# Multiple R-squared:  0.2746,	Adjusted R-squared:  0.2744 
# F-statistic:  1958 on 1 and 5173 DF,  p-value: < 2.2e-16
# 

compare_models(chew_mod,pp_mod, 'chew_mod','pp_mod')
# === ANOVA ===
#   Analysis of Variance Table
# 
# Model 1: RW064 ~ chew_grade_derived_mine_loose
# Model 2: RW064 ~ pseudoprog
# Res.Df    RSS Df Sum of Sq F Pr(>F)
# 1   5173 5493.9                      
# 2   5173 4457.3  0    1036.6         
# 
# === BIC ===
#   df      BIC
# chew_mod  3 15021.09
# pp_mod    3 13939.04
# 
# === AIC ===
#   df      AIC
# chew_mod  3 15001.44
# pp_mod    3 13919.39
# 
# === Summary: chew_mod ===
#   
#   Call:
#   lm(formula = RW064 ~ chew_grade_derived_mine_loose, data = dat)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.0179 -0.2111  0.3196  0.7429  1.1662 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                    0.760247   0.033918   22.41   <2e-16 ***
#   chew_grade_derived_mine_loose -0.211645   0.008551  -24.75   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.031 on 5173 degrees of freedom
# Multiple R-squared:  0.1059,	Adjusted R-squared:  0.1057 
# F-statistic: 612.5 on 1 and 5173 DF,  p-value: < 2.2e-16
# 
# 
# === Summary: pp_mod ===
#   
#   Call:
#   lm(formula = RW064 ~ pseudoprog, data = dat)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.6058 -0.2832  0.1849  0.6475  2.2618 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.939779   0.024863   37.80   <2e-16 ***
#   pseudoprog  -0.092014   0.002079  -44.25   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9282 on 5173 degrees of freedom
# Multiple R-squared:  0.2746,	Adjusted R-squared:  0.2744 
# F-statistic:  1958 on 1 and 5173 DF,  p-value: < 2.2e-16
# 

compare_models(yrshift_mod,pp_mod, 'yrshift_mod','pp_mod')
# === ANOVA ===
#   Analysis of Variance Table
# 
# Model 1: RW064 ~ shift_year_norm
# Model 2: RW064 ~ pseudoprog
# Res.Df    RSS Df Sum of Sq F Pr(>F)
# 1   5173 4371.7                      
# 2   5173 4457.3  0   -85.548         
# 
# === BIC ===
#   df      BIC
# yrshift_mod  3 13838.76
# pp_mod       3 13939.04
# 
# === AIC ===
#   df      AIC
# yrshift_mod  3 13819.10
# pp_mod       3 13919.39
# 
# === Summary: yrshift_mod ===
#   
#   Call:
#   lm(formula = RW064 ~ shift_year_norm, data = dat)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.6757 -0.2261  0.2001  0.5560  2.9263 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      1.049000   0.026241   39.98   <2e-16 ***
#   shift_year_norm -0.063264   0.001381  -45.80   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9193 on 5173 degrees of freedom
# Multiple R-squared:  0.2885,	Adjusted R-squared:  0.2884 
# F-statistic:  2098 on 1 and 5173 DF,  p-value: < 2.2e-16
# 
# 
# === Summary: pp_mod ===
#   
#   Call:
#   lm(formula = RW064 ~ pseudoprog, data = dat)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.6058 -0.2832  0.1849  0.6475  2.2618 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.939779   0.024863   37.80   <2e-16 ***
#   pseudoprog  -0.092014   0.002079  -44.25   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9282 on 5173 degrees of freedom
# Multiple R-squared:  0.2746,	Adjusted R-squared:  0.2744 
# F-statistic:  1958 on 1 and 5173 DF,  p-value: < 2.2e-16


# CLEAR LINGERING PROCESSES ----
system("ps aux | grep 'R --no-echo' | grep -v grep ")
# # system("ps aux | grep 'R --no-echo' | grep -v grep ®| awk '{print $2}' | xargs kill")
