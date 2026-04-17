

# (((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))
# Setup ----
pacman::p_load(shapviz, patchwork, tidyr, dplyr, reshape2, purrr, GGally, tidyverse, ggplot2, gridExtra, ComplexHeatmap, viridis, zoo, mgcv, zoo, tools, ggrepel, car, caret, stats, lme4, lmerTest, ape, SingleCellExperiment, scater, velociraptor, igraph, Rtsne, xgboost)

rm(list = ls())
gc()


source("/Users/canderson/Documents/LMRI/PRJ2024001/version015/analysis/00_MyFunctions.R")
how.long()

# (((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))
# Load Data ----

wd_path = "~/Documents/LMRI/PRJ2024001/version015/"

# load experiment
sce <-  HDF5Array::loadHDF5SummarizedExperiment(sprintf("%sprocessed_data/sce002.2.h5", wd_path))


a_path <- paste0(wd_path, "results/clustering/09_clustering/9.3.1_featureClusterShap")
if(!dir.exists(a_path)) {
  dir.create(a_path, recursive = TRUE)
  
}


# (((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))
# PREDICT CLUSTER GIVEN FEATURES ----
assay <- assays(sce)$pmm_imp_of_scaled %>% 
  t() %>% 
  as.data.frame()

assay <-cbind(assay, my_cluster = sce$my_cluster)

inds <- sample(1:nrow(assay), .7*nrow(assay))

x_train <- assay[inds,-which(colnames(assay) == "my_cluster")] %>% 
  data.matrix()
y_train <- assay$my_cluster[inds]
x_test <- assay[-inds,-which(colnames(assay) == "my_cluster")] %>% 
  data.matrix()
y_test <-  assay$my_cluster[-inds]

# # Define parameters
# params <- list(
#   objective = "multi:softmax",  # For classification
#   num_class = n_distinct(sce$my_cluster), # Number of classes
#   eval_metric = "mlogloss"     # Multiclass log-loss
# )

# xgb <- xgboost(as.matrix(x_train),
#                label = y_train,
#                nrounds = 20
#                )
# # accuracy
# pred <- predict(xgb, x_test)
# confusion <- table(predicted = round(pred,0), actual = y_test)
# confusion <- confusion[as.character(1:11), as.character(1:11)]
# sum(diag(confusion))/ sum(confusion)
# 
# # importance
# gain <- setNames(xgb.importance(model = xgb)$Gain, xgb.importance(model = xgb)$Feature)
# par(mar = c(10,4,4,2))
# gain %>% sort  %>% barplot(las = 2, main = "XGB importance score")
# 
# 
# # look at proportion imputed
# orig_assay <- assays(sce)$measure_vars %>%
#   t() %>%
#   as.data.frame()
# # ratio imputed to total
# missingness <- colSums(is.na(orig_assay)) / nrow(orig_assay)
# 
# setdiff(names(missingness), names(gain))
# 
# # dataframe of missingness
# ratio_missing_d <- as.data.frame(missingness) %>%
#   rownames_to_column("Feature") %>%
#   rename_with(~c("Feature", "ratio_missing"))
# 
# # df of gain
# gain_d <- as.data.frame(gain) %>%
#   rownames_to_column("Feature") %>%
#   rename_with(~c("Feature", "gain_importance")) %>%
#   filter(Feature %in% names(missingness))
# 
# # join two together
# gain_and_missingness <- left_join(gain_d, ratio_missing_d)
# 
# # plot gain~missingness
# gain_and_missingness %>%
#   ggplot(aes(ratio_missing*100, gain_importance, color = Feature))+
#   geom_point(show.legend = FALSE)+
#   geom_label_repel(aes(label = Feature), show.legend = FALSE, max.overlaps = 20)+
#   ggtitle("Percent missing versus XGB importance (gain)")+
#   labs(x = "Percent Missing", y = "Gain")+
#   theme_bw()
# ggsave(paste0(wd_path, "results/clustering/09_clustering/9.3.1_featureClusterShap/feature_gain_by_missingness.pdf"),height = 10, width = 10)


# SHAP ----
params <- list(objective = "multi:softprob", num_class = n_distinct(assay$my_cluster), nthread = 1)
dtrain <- xgb.DMatrix( data.matrix(select(assay, -my_cluster)) , 
                       label = as.integer(assay$my_cluster) - 1,
                       nthread = 1)
fit <- xgb.train(params = params, data = dtrain, nrounds = 50)

# Create "mshapviz" object (logit scale)
(x <- shapviz(fit, X_pred = data.matrix(select(assay, -my_cluster)), X = assay))

# # Contains "shapviz" objects for all classes
# all.equal(x[[3]], shapviz(fit, X_pred = data.matrix(select(assay, -my_cluster)), X = assay, which_class =  n_distinct(assay$my_cluster) ))

# Better names
names(x) <- levels(assay$my_cluster)
x
#> 'mshapviz' object representing 3 'shapviz' objects:
#>    'setosa': 150 x 4 SHAP matrix
#>    'versicolor': 150 x 4 SHAP matrix
#>    'virginica': 150 x 4 SHAP matrix

bg_job <- parallel::mcparallel({
  beeswarm <- sv_importance(x, max_display = 10, kind = "beeswarm")
  ggsave(paste0(wd_path, "results/clustering/09_clustering/9.3.1_featureClusterShap/10var_shaps_beeswarm.png"),beeswarm, height = 10, width = 20)
}, mc.set.seed = FALSE)

bg_jobx <- parallel::mcparallel({
  beeswarm <- sv_importance(x, max_display = 20, kind = "beeswarm")
  ggsave(paste0(wd_path, "results/clustering/09_clustering/9.3.1_featureClusterShap/20var_shaps_beeswarm.pdf"),beeswarm, height = 20, width = 30)
}, mc.set.seed = FALSE)


# mccollect()

# # (((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))
# # FIND what values perfectly fit each cluster given shap values ----
# ##+ extract shaps, values, and baselines
# ##+ individual level shap = baseline + sum(shaps across values)
# ##+ 
# 
# values <- purrr::map(x, `[[`,which(names(x$`1`) =="X") )
# shaps <- purrr::map(x, `[[`, which(names(x$`1`) =="S"))
# baselines <- purrr::map(x, `[[`, 3)
# 
# # add rownames to shaps
# shaps <- Map(function(x,y){
#   rownames(x) <- rownames(y)
#   return(x)
# }, shaps, values)
# 
# 
# # calculate each sample's shap score: sum(shaps) + baseline = indiv. shap score
# samp_shaps <- Map(function(x,y) {
#   overall_shap <- rowSums(x) + y
#   overall_shap <- set_names(overall_shap,rownames(x))
#   return(overall_shap)
#   }, shaps, baselines)
# 
# # find sample with highest shap for each cluster classification
# cluster_reps <- purrr::map(samp_shaps, which.max)
# 
# # find sample with highest shap for each cluster classification for where no missing in og data
# txtplot::txtdensity(sce$msngns_ratio, pch = '•')
# non_missy_samps <- colnames(sce)[sce$msngns_ratio<.1]
# 
# # find sample with highest shap for each cluster classification for where no missing in og data
# non_missy_cluster_reps <- purrr::map(samp_shaps, ~{
#   x <- .
#   z <- x[names(x)%in%non_missy_samps]
#   names(which.max(z))
# }) 
# 
# 
# # plot where subjects lie
# ggplot()+
#   geom_point(data = reducedDim(sce, "TSNE"), 
#              aes(TSNE1, TSNE2), 
#              color = "darkgrey", alpha = .4)+
#   geom_point(data = reducedDim(sce[ , sce$sample%in%unlist(non_missy_cluster_reps)], 'TSNE'), 
#              aes(TSNE1, TSNE2, color = sce$my_cluster[sce$sample%in%unlist(non_missy_cluster_reps)]),
#              size = 5)+
#   scale_color_manual(values = cluster_cols)+
#   ggtitle("Samples with highest SHAP in each cluster")+
#   labs(color = "my_cluster")+
#   coord_equal()+
#   theme_bw()
# ggsave(sprintf('%sresults/clustering/09_clustering/9.3.1_featureClusterShap/shap_cluster_representatives.png', wd_path), width = 8, height = 8)
# 
# 
# # save subjects
# cluster_reps_d <- as.matrix(non_missy_cluster_reps) %>% 
#   as.data.frame() %>% 
#   rownames_to_column("cluster") %>% 
#   rename_with(~c('cluster', 'sample')) %>% 
#   mutate(sample = unlist(sample)) %>% 
#   left_join(colData(sce) %>% as.data.frame() ) %>% 
#   left_join(assays(sce)$measure_vars %>% t() %>% as.data.frame() %>% rownames_to_column("sample"), by = "sample")
# 
# write_csv(cluster_reps_d, sprintf("%sprocessed_data/shap_cluster_representatives.csv", wd_path))
# 
# # plot subjects heatmap
# row_order = metadata(sce)$other$cluster_agg_hm_ordering$row_order
# column_order = metadata(sce)$other$cluster_agg_hm_ordering$col_order
# 
# m <- cluster_reps_d %>% 
#   select(all_of(row_order)) %>% 
#   scale() %>% 
#   t() %>% 
#   na.omit() %>% 
#   t() 
# 
# 
# row_order <- row_order[row_order%in%colnames(m)]
# rownames(m) <- paste("samp:", cluster_reps_d$sample,":clust:", cluster_reps_d$my_cluster, sep = "")
# 
# column_order_inds <- match(column_order, str_extract(colnames(t(m)), "\\d+$"))
# row_order_inds <-  match(row_order, rownames(t(m)))
# 
# row_hclust <- hclust(dist(t(m)))
# col_hclust <- hclust(dist(m))
# 
# # Turn into dendrograms
# row_dend <- as.dendrogram(row_hclust)
# col_dend <- as.dendrogram(col_hclust)
# 
# dendextend::order.dendrogram(row_dend) <- row_order_inds
# dendextend::order.dendrogram(col_dend) <- column_order_inds
# 
# 
# # Plot with fixed order *and* visible dendrograms
# hm <- Heatmap(
#   t(m),
#   name = "Scaled feature value",
#   # cluster_rows = TRUE,
#   # cluster_columns = TRUE,
#   # clustering_method_columns = "ward.D2",
#   # clustering_method_rows = "ward.D2"
#   cluster_rows = row_dend,
#   cluster_columns = col_dend,
#   row_order = match(row_order, rownames(t(m))),
#   column_order = match(column_order, str_extract(colnames(t(m)), "\\d+$")),
#   row_dend_reorder = FALSE,
#   column_dend_reorder = FALSE
# )
# 
# pdf(sprintf('%sresults/clustering/09_clustering/9.3.1_featureClusterShap/shap_cluster_representatives_hm.pdf', wd_path), height = 15, width = 10)
# draw(hm,)
# dev.off()


# (((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))
# FIND what values perfectly fit each cluster given XGB probs ----
##+ extract shaps, values, and baselines
##+ individual level shap = baseline + sum(shaps across values)
##+


# predict cluster for each patient; returns vector of each samplexcluster probability
pred <- predict(fit, newdata = dtrain)  # vector of length nrow * num_class

# Reshape into matrix: rows = observations, cols = classes
pred_mat <- matrix(pred, nrow = nrow(dtrain), ncol = fit$params$num_class, byrow = TRUE)
rownames(pred_mat) <- rownames(assay)

# find sample with highest shap for each cluster classification for where no missing in og data
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

ggsave(sprintf('%sresults/clustering/09_clustering/9.3.1_featureClusterShap/shap_cluster_representatives.png', wd_path), width = 8, height = 8)


# save subjects
cluster_reps_d <- as.matrix(cluster_reps) %>%
  as.data.frame() %>%
  rownames_to_column("cluster") %>%
  rename_with(~c('cluster', 'sample')) %>%
  mutate(sample = unlist(sample)) %>%
  left_join(colData(sce) %>% as.data.frame() ) %>%
  left_join(assays(sce)$measure_vars %>% t() %>% as.data.frame() %>% rownames_to_column("sample"), by = "sample")

write_csv(cluster_reps_d, sprintf("%sprocessed_data/shap_cluster_representatives.csv", wd_path))

# plot subjects heatmap
row_order = metadata(sce)$other$cluster_agg_hm_ordering$row_order
column_order = metadata(sce)$other$cluster_agg_hm_ordering$col_order

m <- cluster_reps_d %>%
  select(all_of(as.character(row_order))) %>%
  scale() %>%
  t() %>%
  na.omit() %>%
  t()

row_order <- row_order[row_order%in%colnames(m)]
rownames(m) <- paste("samp:", cluster_reps_d$sample,":clust:", cluster_reps_d$my_cluster, sep = "")

column_order_inds <- match(column_order, str_extract(colnames(t(m)), "\\d+$"))
row_order_inds <-  match(row_order, rownames(t(m)))

row_hclust <- hclust(dist(t(m)))
col_hclust <- hclust(dist(m))

# Turn into dendrograms
row_dend <- as.dendrogram(row_hclust)
col_dend <- as.dendrogram(col_hclust)

dendextend::order.dendrogram(row_dend) <- row_order_inds
dendextend::order.dendrogram(col_dend) <- column_order_inds


# Plot with fixed order *and* visible dendrograms
hm_manual_order <- Heatmap(
  t(m),
  name = "Scaled feature value",
  cluster_rows = row_dend,
  cluster_columns = col_dend,
  row_order = match(row_order, rownames(t(m))),
  column_order = match(column_order, str_extract(colnames(t(m)), "\\d+$")),
  row_dend_reorder = FALSE,
  column_dend_reorder = FALSE
)

pdf(sprintf('%sresults/clustering/09_clustering/9.3.1_featureClusterShap/shap_cluster_representatives_hm_manual_order.pdf', wd_path), height = 15, width = 10)
draw(hm_manual_order)
dev.off()

# Plot 
hm_auto_order <- Heatmap(
  t(m),
  name = "Scaled feature value",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_method_columns = "ward.D2",
  clustering_method_rows = "ward.D2"
)

pdf(sprintf('%sresults/clustering/09_clustering/9.3.1_featureClusterShap/shap_cluster_representatives_hm_auto_order.pdf', wd_path), height = 15, width = 10)
draw(hm_auto_order)
dev.off()

# (((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))(((((((((((((((((((())))))))))))))))))))
# SHAP plots ----

# bg_job <- parallel::mcparallel({
#   beeswarm <- sv_importance(x, max_display = 20, kind = "beeswarm", ncol = 5)
#   ggsave(paste0(wd_path, "results/clustering/09_clustering/9.3.1_featureClusterShap/var_shaps_beeswarm.png"),beeswarm, height = 20, width = 30)
# }, mc.set.seed = FALSE)

# bar <- sv_importance(x, kind = "bar",  viridis_args = list(option = "turbo"))
# ggsave(paste0(wd_path, "results/clustering/09_clustering/9.3.1_featureClusterShap/var_shaps_bar.pdf"), bar,height = 30, width = 8)

shap_scores_orig <- sv_importance(x,kind = "no" ,  viridis_args = list(option = "turbo"))

pdf(paste0(wd_path, "results/clustering/09_clustering/9.3.1_featureClusterShap/feature_SHAPS_HM.pdf"), height = 20, width = 8)
Heatmap(shap_scores_orig, 
        column_title = "Feature SHAP scores", 
        # col =colorRampPalette(c("#02f2f5", "#fc2cff"))(100)
        col =colorRampPalette(c("white", "red"))(100)
        )
dev.off()

shap_scores <- log(shap_scores_orig)
shap_scores[is.infinite(shap_scores)] <- 0

pdf(paste0(wd_path, "results/clustering/09_clustering/9.3.1_featureClusterShap/feature_SHAPS_log_HM.pdf"), height = 20, width = 8)
Heatmap(shap_scores, 
        column_title = "Feature SHAP logged scores",
        # col =colorRampPalette(c("#02f2f5", "#fc2cff"))(100)
        col =colorRampPalette(c("white", "red"))(100)
        )
dev.off()

# which vars max for each cluster
top_shaps <- setNames(rownames(shap_scores_orig)[apply(shap_scores_orig, 2, which.max)], colnames(shap_scores_orig))
expressionPlot <- plotExpression(sce, x = "label", features = top_shaps, exprs_values = "pmm_imp_of_scaled", color_by = "my_cluster")+
  geom_hline(yintercept = 0, lty = 3, alpha = .5)
ggsave(paste0(wd_path, "results/clustering/09_clustering/9.3.1_featureClusterShap/top_shaps_expression_plot.pdf"), expressionPlot,
       height = 20, width = 20)


# check if process done
result <- parallel::mccollect(bg_job, wait = FALSE)
if(is.null(result)){
  cat("Awaiting Process")
  while(is.null(result)){
    cat(" . ")
    Sys.sleep(2)
    result <- parallel::mccollect(bg_job, wait = FALSE)
  }
  cat("\nProcess Completed\n")
  
}




