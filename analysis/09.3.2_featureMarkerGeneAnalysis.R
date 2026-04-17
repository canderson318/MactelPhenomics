library(tidyr)
library(dplyr)
library(reshape2)
library(purrr)
library(tidyverse)
library(gridExtra)
library(ComplexHeatmap)
library(viridis)
library(zoo)
library(mgcv)
library(ggrepel)
library(car)
library(caret)
library(stats)
library(ape)
library(SingleCellExperiment)
library(scater)
library(scran)


rm(list = ls())
gc()

# FUNCTIONS----
source("/Users/canderson/Documents/LMRI/PRJ2024001/version015/analysis/00_MyFunctions.R")
how.long()


# Load Data ----

wd_path = "~/Documents/LMRI/PRJ2024001/version015/"

# load experiment
sce <-  HDF5Array::loadHDF5SummarizedExperiment(sprintf("%sprocessed_data/sce002.2.h5", wd_path))


a_path <- paste0(wd_path, "results/clustering/09_clustering/9.3.2_featureMarkerGeneAnalysis/")
if(!dir.exists(a_path)) {
  dir.create(a_path, recursive = TRUE)
  
}

# MARKER GENE ANALYSIS ----


marker.info <- scoreMarkers(sce, colLabels(sce), assay.type = "pmm_imp_of_scaled")

# order by mean AUC 
marker.info <- lapply(marker.info, function(a_DataFrame){
  a_DataFrame <- a_DataFrame[order(a_DataFrame$mean.AUC, decreasing = TRUE) , ]
  a_data.frame <- as.data.frame(a_DataFrame)
  return(a_data.frame)
})


marker.info %>% lapply(function(d){
  d %>% 
    arrange(rank.AUC) %>% 
    head(11)
})

all_markers <- lapply(names(marker.info), function(name){
  x <- marker.info[[name]]
  x$cluster <- name
  x$Feature <- rownames(x)
  return(x)
  }) %>%
  do.call(rbind,.)

# pc <- prcomp(all_markers %>% select(where(is.numeric)))
# 
# (pc$sdev^2/sum(pc$sdev^2)) %>% 
#    cumsum(.) %>%
#    txtplot(., height = 30)

all_markers %>%
  ggplot(aes(median.AUC, reorder_within(Feature, median.AUC, cluster) ,fill= cluster))+
  geom_bar(stat = 'identity', position = "dodge", show.legend  = FALSE)+
  scale_fill_manual(values = cluster_cols)+
  facet_wrap(~cluster, scales = "free", ncol = 6)+
  ggtitle("AUC for each cluster marker features")+
  theme_minimal()
