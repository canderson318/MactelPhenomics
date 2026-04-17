

options(max.print = 100)


# Packages
pacman::p_load(tidyr,reshape2,purrr,GGally,tidyverse,ggplot2,ComplexHeatmap,viridis,ggrepel,SingleCellExperiment,scater,skimr,dplyr,HDF5Array,parallel,brms)
rm(list = ls())
gc()

# script with all functions i have made
source("/Users/canderson/Documents/LMRI/PRJ2024001/version015/analysis/00_MyFunctions.R")

# Load Data   
## sce
wd_path = "~/Documents/LMRI/PRJ2024001/version015/"

how.long()
# load experiment
sce <-  loadHDF5SummarizedExperiment(paste0(wd_path, "processed_data/sce006.h5"))
how.long()


# rename clusters
## order of clusteres along pseudoprog
new_cluster_names <-  data.frame(colData(sce))%>% 
  select(my_cluster, pseudoprog) %>% 
  summarize(mean_pp = mean(pseudoprog),median_pp = median(pseudoprog), .by = my_cluster) 
new_cluster_names <- new_cluster_names[order(new_cluster_names$mean_pp),]
new_cluster_names$new_cluster <- seq(nrow(new_cluster_names))
## assign new cluster based on ppmean
sce$old_my_cluster <- sce$my_cluster
sce$my_cluster <- factor(new_cluster_names$new_cluster[match(as.character(sce$my_cluster), as.character(new_cluster_names$my_cluster))])


# ANALYSIS ----
assay <- assays(sce)$pmm_imp_of_scaled %>% 
  t() %>% 
  data.frame() %>% 
  cbind(cluster = sce$my_cluster) 

d <- assay%>% 
  pivot_longer(-cluster)

pars_base <- list(
  formula = value ~ 1,
  family = gaussian(),
  sample_prior = "no",
  prior = prior(normal(0,5), class = "Intercept"),
  backend = "cmdstanr",
  chains = 3,
  cores = 3,
  iter = 1000,
  algorithm = "sampling",
  silent = 2
)


combs <- select(d, cluster, name) %>% distinct()

results <- data.frame()
cat("\n")
for(i in seq(nrow(combs))){
  
  clust <- combs[[i,1]]
  var <- combs[[i,2]]
  
  df_group = d %>% 
    filter(cluster == clust, name == var)
  
  # # Fit BRM model
  # pars <- append(pars_base, list(data = df_group))
  # fit <- do.call(brm, pars)
  # # Extract posterior summary
  # post <- posterior_summary(fit)
  # mu_est = post[1,1]
  # mu_err = post[1,2]
  # mu_Q0.5 = post[1,3]
  # mu_Q99.5 = post[1,4]
  # sigma_est = post[2,1]
  # sigma_err = post[2,2]
  # sigma_Q0.5 = post[2,3]
  # sigma_Q99.5 = post[2,4]
  
  # Fit lm
  n <- nrow(df_group)
  if(n < 2) next  # skip underpowered groups

  fit <- lm(value ~ 1, data = df_group)
  mu_est <- coef(fit)[1]
  mu_err <- suppressWarnings(summary(fit)$coefficients[1,2])
  mu_Q0.5 <- confint(fit)[1]
  mu_Q99.5 <- confint(fit)[2]

  sigma_est <- sd(df_group$value)
  sigma_err <- sigma_est / sqrt(2 * (n - 1))  # std error of sample SD
  # Approximate CI for SD using chi-squared distribution
  chi2_lower <- qchisq(0.005, df = n - 1)
  chi2_upper <- qchisq(0.995, df = n - 1)
  sigma_Q0.5 <- sqrt((n - 1) * sigma_est^2 / chi2_lower)
  sigma_Q99.5 <- sqrt((n - 1) * sigma_est^2 / chi2_upper)

  
  # Return results as tibble
  res <- tibble(
    cluster = clust,
    name = var,
    mu_est = mu_est,
    mu_err = mu_err,
    mu_Q0.5 = mu_Q0.5,
    mu_Q99.5 = mu_Q99.5,
    sigma_est = sigma_est,
    sigma_err = sigma_err,
    sigma_Q0.5 = sigma_Q0.5, 
    sigma_Q99.5 = sigma_Q99.5
  )
  
  results <- bind_rows(results, res)
  
  cat(" | ")
}
cat("\n")

# Set up factor levels
results <- results %>%
  mutate(name = factor(name, levels = rev(sort(unique(name)))))

# ggplot(results, aes(y = factor(name, levels = rev(sort(unique(name)))), x = mu_est)) +
#   geom_point(size = .5) +
#   geom_errorbar(aes(xmin = mu_est - mu_err, xmax = mu_est + mu_err), width = 0.2) +
#   geom_ribbon(aes(xmin = mu_Q0.5, xmax = mu_Q99.5,group = 1),
#               fill = "blue", alpha = 0.2, inherit.aes = TRUE) +
#   facet_wrap(~ factor(cluster), scales = "free_x", nrow = 1) +
#   labs(title = "Cluster-level posterior means with SD bars and 95% credible interval ribbons") +
#   theme_minimal()


results %>% 
  group_by(name) %>%
  mutate(across(-c( cluster), ~ as.numeric(scale(.)))) %>%
  ungroup() %>% 
  ggplot(aes(y = tidytext::reorder_within(name,mu_est, cluster))) +
  # geom_vline(xintercept = c(-.5,.5), linetype  = 3, color = "darkgrey")+
  geom_vline(xintercept = c(-1,1), linetype  = 3, color = "darkgrey")+
  geom_vline(xintercept = 0, linetype  = 1, color = "darkgrey")+
  geom_hline(aes(yintercept = name), linetype =1, color = 'darkgrey', alpha =.3)+
  geom_ribbon(aes(xmin = mu_Q0.5,
                  xmax = mu_Q99.5,
                  fill = cluster,
                  group = cluster),
              show.legend = F,
              alpha = 0.5) +
  geom_segment(aes(x = mu_Q0.5,
                  xend = mu_Q99.5,
                  color = cluster,
                  group = cluster),
              show.legend = F,
              alpha = 1) +
  # geom_point(aes(x = mu_est, color = cluster, size = abs(mu_est)),show.legend = F)+
  scale_size_continuous(range = c(.1, 5))+
  tidytext::scale_y_reordered()+
  labs(
    x = "Mean",
    y = "Variable",
    title = "95% Confidence Interval Ribbons of Cluster Feature Means, Scaled Across Cluster",
    fill = "Cluster"
  ) +
  facet_wrap(~cluster, scales = "free", nrow = 3)+
  theme_classic()


# plot(abs(results$mu_est)+1,abs(results$mu_Q99.5-results$mu_Q0.5)+1 , log = "xy")
# lm((abs(results$mu_Q99.5-results$mu_Q0.5)+1 ) ~ (abs(results$mu_est)+1)) %>% summary()
