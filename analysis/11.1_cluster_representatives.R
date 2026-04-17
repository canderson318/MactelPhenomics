pacman::p_load(SingleCellExperiment, scran, scuttle, dplyr, ggplot2, scater, reshape2, stringr, tibble, xgboost)

rm(list=ls())
gc()



wd_path <- "~/Library/CloudStorage/Egnyte-lmri/Shared/LMRI/DS_Group/PRJ2024001/analysis_versions/version015/"

# functions
source(sprintf("%sanalysis/00_MyFunctions.R",wd_path))

# load data
sce <- HDF5Array::loadHDF5SummarizedExperiment(sprintf("%sprocessed_data/sce004.h5", wd_path))


# (((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))
# REDO CLUSTER ORDERING ----
# invert PC1
## rdims
pcs <- reducedDim(sce, "PCA")
pcs[,"PC1"] <- pcs[,"PC1"] * -1
reducedDim(sce, "PCA") <- pcs
## tree
pc_tree <- metadata(sce)$trajectory_data$pca_line_data
pc_tree$PC1 <- pc_tree$PC1 * -1
metadata(sce)$trajectory_data$pca_line_data <- pc_tree

# rename clusters
## order of clusteres along pseudoprog
new_cluster_names <-  data.frame(colData(sce))%>% 
  dplyr::select(my_cluster, pseudoprog) %>% 
  summarize(mean_pp = mean(pseudoprog),median_pp = median(pseudoprog), .by = my_cluster) 

with(new_cluster_names,{
  plot(mean_pp,median_pp, type = "n")
  text(mean_pp,median_pp, labels = my_cluster)
  add_lines()
})

new_cluster_names <- new_cluster_names[order(new_cluster_names$mean_pp),]
new_cluster_names$new_cluster <- seq(nrow(new_cluster_names))

## assign new cluster based on ppmean
sce$old_my_cluster <- sce$my_cluster
sce$my_cluster <- factor(new_cluster_names$new_cluster[match(as.character(sce$my_cluster), as.character(new_cluster_names$my_cluster))])

# ****** Function to replace each part of the edge ******
replace_edge_clusters <- function(edge_strs) {
  sapply(edge_strs, function(edge_str){
    parts <- strsplit(edge_str, "--")[[1]]
    new_parts <- sapply(parts, function(x) {
      new_val <- new_cluster_names$new_cluster[match(as.integer(x), new_cluster_names$my_cluster)]
      as.character(new_val)
    })
    paste(new_parts, collapse = "--")
  })
}

## replace in trees
metadata(sce)$trajectory_data$pca_line_data$old_edge<- metadata(sce)$trajectory_data$pca_line_data$edge
metadata(sce)$trajectory_data$pca_line_data$edge <- replace_edge_clusters(metadata(sce)$trajectory_data$pca_line_data$old_edge)

metadata(sce)$trajectory_data$tsne_line_data$old_edge<- metadata(sce)$trajectory_data$tsne_line_data$edge
metadata(sce)$trajectory_data$tsne_line_data$edge <- replace_edge_clusters(metadata(sce)$trajectory_data$tsne_line_data$old_edge)


# (((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))
# PREDICT CLUSTER GIVEN FEATURES ----
assay <- assays(sce)$pmm_imp_of_scaled %>% 
  t() %>% 
  as.data.frame()

assay <-cbind(assay, my_cluster = sce$my_cluster)

inds <- sample(1:nrow(assay), .7*nrow(assay))

stop("\nProbably unecessary to split 70/30 to get probs for each sample, consider revising")

x_train <- assay[inds,-which(colnames(assay) == "my_cluster")] %>% 
  data.matrix()
y_train <- assay$my_cluster[inds]
x_test <- assay[-inds,-which(colnames(assay) == "my_cluster")] %>% 
  data.matrix()
y_test <-  assay$my_cluster[-inds]


# TRAIN AND EXTRACT PROBS ----
params <- list(objective = "multi:softprob", num_class = n_distinct(assay$my_cluster), nthread = 1)
dtrain <- xgb.DMatrix( data.matrix(select(assay, -my_cluster)) , 
                       label = as.integer(assay$my_cluster) - 1,
                       nthread = 1)
fit <- xgb.train(params = params, data = dtrain, nrounds = 50)


# predict cluster for each patient; returns vector of each samplexcluster probability
pred <- predict(fit, newdata = dtrain)  # vector of length nrow * num_class

# Reshape into matrix: rows = observations, cols = classes
pred_mat <- matrix(pred, nrow = nrow(dtrain), ncol = fit$params$num_class, byrow = TRUE)
rownames(pred_mat) <- rownames(assay)

# find sample with highest probabilities for each cluster classification for where no missing in og data
txtplot::txtdensity(sce$msngns_ratio, pch = '•')
non_missy_samps <- colnames(sce)[sce$msngns_ratio<.2]

# filter pred mat for non missy samples
pred_mat_filt <- pred_mat[rownames(pred_mat)%in%non_missy_samps, ]


apply(pred_mat_filt,2, boxplot)

# find sample with highest prob for each cluster classification
cluster_reps <- setNames(rownames(pred_mat_filt)[apply(pred_mat_filt,2, which.max)], seq(fit$params$num_class))

# plot where subjects lie
ggplot()+
  geom_point(data = reducedDim(sce, "TSNE"),
             aes(TSNE1, TSNE2),
             color = "darkgrey", alpha = .4)+
  geom_point(data = reducedDim(sce[ , sce$sample%in%unlist(cluster_reps)], 'TSNE'),
             aes(TSNE1, TSNE2, color = sce$my_cluster[sce$sample%in%unlist(cluster_reps)]),
             size = 5)+
  scale_color_manual(values = cluster_cols)+
  ggtitle("Samples with highest SHAP in each cluster")+
  labs(color = "my_cluster")+
  coord_equal()+
  theme_bw()

# save subjects
cluster_reps_d <- as.matrix(cluster_reps) %>%
  as.data.frame() %>%
  rownames_to_column("cluster") %>%
  rename_with(~c('cluster', 'sample')) %>%
  mutate(sample = unlist(sample)) %>%
  left_join(colData(sce) %>% as.data.frame() ) %>%
  left_join(assays(sce)$measure_vars %>% t() %>% as.data.frame() %>% rownames_to_column("sample"), by = "sample")

