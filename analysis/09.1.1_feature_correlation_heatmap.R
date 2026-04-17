# Wed May  7 17:37:49 2025 ------------------------------

pacman::p_load(tidyr,reshape2,purrr,tidyverse,ggplot2,gridExtra,ComplexHeatmap,viridis,zoo,mgcv,mice, htmltools,Metrics,psych,car,stats,lme4,lmerTest,ape,SingleCellExperiment,scater,dplyr,HDF5Array,parallel,txtplot)

rm(list = ls())
gc()


# LOAD DATA ----
wd_path = "~/Documents/LMRI/PRJ2024001/version015/"


sce <- HDF5Array::loadHDF5SummarizedExperiment(sprintf("%sprocessed_data/sce001.h5", wd_path))

# FEATURE CORRELATIONS ----
assay <- assays(sce)$pmm_imp_of_scaled %>% 
  t() %>% 
  as.matrix()

cor <- cor(assay, use = "complete",method = "pearson")

hm <- Heatmap(cor, 
              clustering_method_columns =  "ward.D2",
              clustering_method_rows=  "ward.D2", 
              name = "R"
                )

pdf(sprintf("%sresults/featurewise_correlation_hm.pdf", wd_path), height = 15, width = 15)
draw(hm)
dev.off()

sprintf("cp %s %s",
        sprintf("%sresults/featurewise_correlation_hm.pdf", wd_path),
        sprintf("%sresults/publication_figures/featurewise_correlation_hm.pdf", wd_path)) %>% 
  system()
