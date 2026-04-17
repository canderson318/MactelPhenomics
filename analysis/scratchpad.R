
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(tibble)
library(e1071)
library(rsvd)

 
# rm()
# gc()


####  OPTIONS ####
options(max.print = 2000)
options(max.print = 1000)
options(max.print = 450)
options(max.print = 300)
options(max.print = 150)
options(max.print = 100)

tools::pskill(jb$pid)
quartz()
dev.off();   options(device = "RStuzdioGD")


#### DATA MUNG ####

# 
# # WHICH FEATURES CAUSE ASYMETRY?----
# pacman::p_load(HDF5Array, SingleCellExperiment, dplyr, stringr,tidyr, tibble, ggplot2)
# 
# wd_path = "/Users/christiananderson/Library/CloudStorage/OneDrive-LowyMedicalResearchInstitute/Database_EDA/analysis_versions/version013/"
# 
# # load experiment
# sce <-  loadHDF5SummarizedExperiment(paste0(wd_path, "processed_data/sce006.h5"))

# d <-
#   t(assays(sce)$pmm_imp_of_scaled) %>% 
#   as.data.frame() %>% 
#   rownames_to_column("subject") %>% 
#   left_join(colData(sce) %>% as.data.frame() %>% 
#               select(subject,patID, patID_eye, month, pseudoprog)) %>% 
#   mutate(eye = stringr::str_extract(patID_eye,"\\d$")) %>% 
#   select(subject,patID, eye, patID_eye, month, pseudoprog,everything()) %>% 
#   mutate(month= factor(month)) %>% 
#   arrange(  patID, month, eye) %>% 
#   filter(month == 0 )
# 
# 
# # calculate difference of all vars between eyes
# d_diffs <- d %>% 
#   group_by(patID, month) %>% 
#   mutate(across(where(is.numeric), 
#                 ~abs(first(.)-last(.)), 
#                 .names = "diff_{.col}")) %>% 
#   select(!where(is.numeric), matches("diff") ) %>% 
#   ungroup() %>% 
#   as.data.frame()
# 
# d_diffs_distinct <- d_diffs %>% 
#   select(-c(eye, subject, patID_eye)) %>% 
#   distinct()
# 
# d_diffs_distinct %>% nrow()
# d_diffs %>% nrow()
# 
# X_diffs_distinct %>% 
#   select(matches("diff_")) %>% 
#   # mutate_all(~log(.)) %>%
#   reshape2::melt() %>% 
#   ggplot(aes(y = value, color = variable))+
#   geom_boxplot(show.legend = FALSE)
# 
# y <- as.matrix(d_diffs_distinct$diff_pseudoprog)
# 
# X <- d_diffs_distinct %>% 
#   select(matches("diff")) %>% 
#   select(-diff_pseudoprog) %>% 
#   as.matrix()
# 
# mods <- apply(X, 2, function(x) 
#   lm(y~x)
#   ) 
# summs <- lapply(mods, summary)
# 
# # names(summs) <- paste("diff_pseudoprog", names(summs), sep= " ~ ")
# 
# results_l <- lapply(summs, broom::tidy) 
#   
# results <- lapply(names(results_l), function(nm){
#   x <- results_l[[nm]]
#   x$term[x$term != "(Intercept)"] <- nm
#   return(x)
# }) %>% 
#   do.call(rbind,.) %>% 
#   filter(!grepl("Intercept", term)) %>% 
#   mutate(q = p.adjust(p.value, method = "BH")) 
# 
#   ggplot(results, aes(estimate, -log10(q), color = term))+
#   geom_hline(yintercept = -log10(.05), linetype = 3, color = 'darkgrey')+
#   geom_point(show.legend = FALSE)+
#   ggrepel::geom_label_repel(aes(label =term), show.legend = FALSE)+
#   ggtitle("What feature differences correalate to severity differences (|R-L|)?", "dPP~dx_i")+
#   theme_bw()

  

# # MIXTURE MODEL ON ASSAY----
# # Load necessary libraries
# library(mclust)
# 
# m <- assays(sce)$pmm_imp_of_scaled %>% 
#   t() %>%
#   data.matrix()
# 
# x <- sce$km_nn_cluster
# 
# # Fit a Gaussian Mixture Model
# gmm_model <- Mclust(m)
# 
# # Summary of the model
# summary(gmm_model)
# 
# means <- gmm_model$parameters$mean 
# colnames(means) <- as.character(seq(ncol(means)))
#   
# ComplexHeatmap::Heatmap(means)
#   
# shape <- gmm_model$parameters$variance$shape %>% 
#   setNames(., colnames(m))
#   
# barplot(shape, las = 2, cex.axis = 1,cex.names = .7)
# 
# sce$gmm_clust <- factor(gmm_model$classification)
# 
# scater::plotTSNE(sce, color_by = "gmm_clust")



# HEATMAP OF COLDATA----
# m <- as.data.frame(colData(sce)) %>% 
#   select(where(is.numeric))
# 
# sds <- colSds(as.matrix(m), na.rm = TRUE)
# 
# hm <- m[, sds>0] %>% 
#   cor(., use = "pairwise.complete.obs", method = "pearson") %>% 
#   ComplexHeatmap::Heatmap(.,cluster_columns = TRUE, cluster_rows = TRUE)
# 
# pdf('Desktop/hm.pdf', bg= "white", height = 15, width = 15)
# draw(hm)
# dev.off()


# assay <- assays(sce)$pmm_imp_of_scaled
# sce %>% colData() %>% names
#  
# boxplot(sce$AGE~sce$has_phgdh_rare_var)
# boxplot(sce$pseudoprog~sce$has_phgdh_rare_var)
# 
# 
# # plot ratio of stratus, cirrus, spectralis
# m <- metadata(sce)$d %>% 
#   select(patID_eye, month, matches("^spectralis.*01$|^stratus.*01$|^cirrus.*01$")) %>% 
#   arrange(month)
# 
# m_sub <- m %>% 
#   select(-patID_eye, -month) %>% 
#   mutate_all(~as.character(.))
# 
# for(j in 1:ncol(m_sub)){
#   
#   m_sub[!is.na(m_sub[,j]), j] <- colnames(m_sub)[j]
#   
# }
# 
# m_sub <- m_sub %>% 
#   cbind(m[,c("patID_eye","month")])
# 
# m_sub %>% 
#   mutate(machines = paste(spectralis_MT01,cirrus_MT01,stratus_MT01) %>% 
#            gsub("NA| ","",.)) %>% 
#   select(patID_eye, month, machines) %>% 
#   drop_na() %>% 
#   ggplot(aes(month, patID_eye))+
#   geom_point(aes(color = machines), show.legend = TRUE)
  
# assay <- assays(sce)$pmm_imp_of_scaled %>% 
#   as.matrix()
# library(trigger)
# trig <- trigger.build(marker=assay,)
# trigger.eigenR2(assay)


# # MODEL ROUTE SCORE WITH STRUCTURAL  FACTORS ----
# m_d <- colData(sce) %>%
#   as.data.frame() %>%
#   dplyr::select( patID_eye, patID,
#                  all_of(c(numeric_responses, categorical_responses, snps, covariates[covariates%in%names(colData(sce))] )),
#                  matches("PA\\d")
#                  ) %>%
#   unique() %>%
#   mutate(across(all_of(c(snps, covariates[covariates%in%names(colData(sce))])),
#                 ~ as.numeric(scale(.))))
# 
# fac_route_summs <- sapply(paste0("PA", 1:9), function(fac) {
#   form <- as.formula( paste(sep = " ",
#     "rf_route_score ~ ",fac," + (1|patID) +SEX+ BL_eye_age + ",paste(paste0("PC", 1:10), collapse = " + ")
#   ))
#   lmer(form, data = m_d)
#   }) %>%
#   map(summary)
# 
# names(fac_route_summs) <- paste("rf_route_score ~ ", names(fac_route_summs))
# 
# fac_route <- extract_results(fac_route_summs)
# fac_route %>%
#   filter(!Term %in% grep("SEX|PC|Intercept|age", fac_route$Term, value = TRUE)) %>%
#   volcano.plot

# # FA OF STRUCTURAL FEATURES ----
# wd_path = "/Users/christiananderson/Library/CloudStorage/OneDrive-LowyMedicalResearchInstitute/Database_EDA/analysis_versions/version013/"
# 
# how.long()
# # load experiment
# sce <-  loadHDF5SummarizedExperiment(paste0(wd_path, "processed_data/sce006.h5"))
# how.long()
# 
# assay <- assays(sce[grep("RW\\d.*4$|SNELL|ACUIT", rownames(sce), invert = TRUE) ,])$pmm_imp_of_scaled %>% 
#   t() %>% 
#   as.matrix()
# 
# set.seed(10329)
# fa.parallel(assay)
# library(psych)
# msa <- psych::KMO(assay)
# par(mar = c(8,4,2,2))
# msa$MSAi %>% sort() %>% barplot(las = 2, horiz = FALSE, cex.names = .5)
# # vars <- names(msa$MSAi)[msa$MSAi>.6]
# 
# # 
# nFacs <- 9
# FA <- psych::fa(assay, nfactors = nFacs, rotate = "varimax")
# 
# FA$loadings %>% Heatmap(.)
# 
# colData(sce) <- colData(sce) %>% 
#   cbind(FA$scores)
# # FA ordering decreasing latency =  MR1    MR2    MR4    MR3    MR6    MR5   
# 
# plots <- list(); for(f in paste("MR", 1:nFacs, sep = "")){
#   plots[[f]] <- plotReducedDim(sce, "TSNE", color_by  = f)
# }
  
# # ROTATE POINTS AROUND CENTER AND CLUSTER ELLIPSES ----
# source("OneDrive - Lowy Medical Research Institute/Database_EDA/analysis_versions/version013/analysis/MyFunctions.R")
# m <- matrix(c(-1,1, 
#               -.4,-3, 
#               -2,-4, 
#               1,3, 
#               3,-2,
#               1,4, 
#               0,3), ncol = 2, byrow = TRUE)
# 
# lim <- c(-10,10)
# plot(m, type = "n",xlim = lim , ylim = lim); lines(m)
# points <- m
# for(i in seq(1,360,5)) {
#   
#   new_points <- rotate_points(m,i)
#   # if(i%%5) new_points[,1] <- new_points[,1] + (10/i)
#   # if(i%%3) new_points[,2] <- new_points[,2] - (5/i)
#   lines(new_points)
#   points <- rbind(points, new_points)
# }
# plot(points, xlim = lim , ylim = lim)
# 
# tsne <- Rtsne::Rtsne(points)$Y
# plot(tsne)
# 
# umap <- umap::umap(tsne)$layout
# plot(umap)
# 
# clust <- kmeans(umap, 3)
# colors <-setNames(RColorBrewer::brewer.pal(n_distinct(clust$cluster), name = 'Set1'), as.character(unique(clust$cluster)))
# plot(umap, col = colors[ match(clust$cluster, names(colors))] )
# 
# library(dbscan)
# clust <- dbscan(umap, eps = 1.5, minPts = 1, borderPoints = TRUE)
# clust
# 
# colors <-setNames(RColorBrewer::brewer.pal(n_distinct(clust$cluster), name = 'Set1'), as.character(unique(clust$cluster)))
# plot(umap, col = colors[ match(clust$cluster, names(colors))] )


# library(mclust)
# clust <- Mclust(mds,G = 3)
# 
# colors <-setNames(RColorBrewer::brewer.pal(n_distinct(clust$classification), name = 'Set1'), as.character(unique(clust$classification)))
# plot(mds, col = colors[ match(clust$classification, names(colors))] )




# # PLOT PSEUDOPROG ~ TIME RELATIONSHIP ----
# # f <- function(x)   (x*13) - 84
# 
# f <- function(x){
#   # x= (x*13) - 84
#   x= (x+ 84)/13
#   x = log(x)
#   # x= x^(1/1.2)
#   # x = ((log(x) *  lag(x,default = 1))^(1))
#   return(x)
#   }
# 
# eyes <- sample(sce$patID_eye, 10)
# sce %>% 
#   colData() %>% 
#   as.data.frame() %>% 
#   filter(patID_eye%in% eyes) %>%
#   # filter(patID_eye=="005001_1") %>%
#   mutate(pseudoprog = cummax_na(pseudoprog),.by = patID_eye) %>% 
#   ggplot() +
#   # geom_smooth(show.legend = TRUE, se = FALSE)
#   geom_line( aes(month,pseudoprog, color = patID_eye), show.legend = TRUE)+
#   # geom_line( aes( month,pseudoprog), color = "black", show.legend = TRUE)+
#   # geom_line( aes(pseudoprog %>% f(), month), color = "blue", show.legend = TRUE)+
#   # geom_abline(slope = 1, intercept = 0, linetype = 2, col = "grey")+
#   # coord_cartesian(xlim = c(0,100), ylim = c(0,100))
#   theme_linedraw()
# 
# sce %>% 
#   colData() %>% 
#   as.data.frame() %>% 
#   # filter(patID_eye%in% eyes) %>%
#   filter(patID_eye=="005001_1") %>%
#   mutate(pseudoprog = cummax_na(pseudoprog),.by = patID_eye) %>% 
#   lm(month~pseudoprog, .)


# # WEIGHTED GRAPH FROM TRANSITION MATRIX ----
# # Extract non-zero entries as edge list
# edge_list <- which(transition_matrix > 0, arr.ind = TRUE)
# weights <- transition_matrix[edge_list]
# 
# # Adjust indices to fit igraph (since R uses 1-based indexing)
# edges <- as.matrix(edge_list)#  - 1  # Convert to 0-based indexing for igraph
# 
# # Create a directed graph
# g <- graph_from_edgelist(edges, directed = TRUE)
# 
# # Assign weights to the edges
# E(g)$weight <- weights
# 
# # Step 3: Layout the Graph as a Tree with Root at Node 1
# layout <- layout_as_tree(g, root = 1, circular = FALSE)
# 
# # Step 4: Adjust Node Positions Based on Weights (Optional)
# # Inverse weights so that larger weights bring nodes closer together
# scaled_weights <- sqrt(E(g)$weight)
# E(g)$weight <- scaled_weights  # Reassign for proximity
# 
# set.seed(10202340);plot(g, 
#                         edge.width = E(g)$weight %>% sqrt(), 
#                         edge.arrow.size = .8
#                         # layout = layout_as_tree(g, 
#                         #                         root = 1,
#                         #                         circular = TRUE,
#                         #                         mode = "all", 
#                         #                         flip.y = FALSE )
#                         )





 # library(tokenizers)
# text <- read_file("downloads/20211110.txt")
# tokens <- tokenize_words(text, strip_punct = TRUE)[[1]] %>% 
#         gsub("_|\\d", "", .)
# 
# freq <- table(tokens) %>%
#         sort() %>% 
#         as.data.frame() %>% 
#         mutate(token = as.character(tokens)) %>% 
#         mutate(length = nchar(token))
# 
# plot(freq$length, freq$Freq)
# freq %>% 
#         arrange(Freq, -length)
 
# sce %>% 
#         colData() %>% 
#         as.data.frame() %>% 
#         select(patID_eye, month, Normalised_PRS) %>% 
#         filter(n_distinct(month)>1, .by = patID_eye) %>% 
#         pull(patID_eye) %>% 
#         n_distinct()


# # TABLES OF ROUTE ASSIGNMENTS ----
# wd_path = "/Users/christiananderson/Library/CloudStorage/OneDrive-LowyMedicalResearchInstitute/Database_EDA/analysis_versions/version013/"
# sce <-  HDF5Array::loadHDF5SummarizedExperiment(paste0(wd_path, "processed_data/sce006.h5"))
# 
# tsne <- as.data.frame(colData(sce)) %>% 
#   select(subject, patID_eye, month, pseudotime, rf_route, km_nn_cluster) %>% 
#   left_join(
#     reducedDim(sce, "TSNE") %>% 
#       as.data.frame() %>% 
#       rownames_to_column("subject")
#   ) %>% 
#   mutate( num_visits = n(), .by = patID_eye)
# 
# 
# 
# tsne %>% 
#   mutate(rf_route = case_when(
#     rf_route == "C" ~ "A", 
#     rf_route == "D" ~ "B", 
#     TRUE ~ rf_route
#   )) %>% 
#   ggplot(.,aes(TSNE1, TSNE2))+
#   # geom_point(aes(color = rf_route))+scale_color_manual(values = c("A" = "#EF5FA7", "C" = "#B1467C",  "B" = "#F9B900", "D" = "#B38600"))+ ggtitle("TSNEs colored by rf_routes")
#   # geom_point(aes(color = num_visits))+scale_colour_viridis_c()+ ggtitle("TSNEs colored by eye's number of visits")
#   geom_point(aes(color = rf_route))+ scale_color_manual(values = c("A" = "#EF5FA7", "B" = "#F9B900", "X" = "#87B1E7"))+ ggtitle("TSNEs colored by rf_routes")
# 
# 
# 
# cluster_cols <- c("1" = "#3A6FA9",
#                   "2" = "#B2C6E5",
#                   "3" = "#EE8635",
#                   "7" = "#A12D27",
#                   "9" = "#77599B",
#                   "11" = "#6E4940",
#                   "4" = "#F2BF85",
#                   "6" = "#A8DD93",
#                   "8" = "#F09D99",
#                   "10" = "#C1B0D3",
#                   "5" = "#509E3E")
# # MST TREE
# mean_cluster_locs <- tsne %>% 
#   group_by(km_nn_cluster) %>% 
#   summarize(across(c(TSNE1,TSNE2), ~mean(.))) %>% 
#   mutate(km_nn_cluster=as.numeric(as.character((km_nn_cluster))) )
# 
# tree_structure <- data.frame(
#   parent = c(1, 2, 2,2, 3, 3, 3, 4, 4, 4),
#   child = c(2, 3, 4, 5,11, 9, 7, 10, 6, 8)
# )
# 
# # Join the parent and child mean coordinates for plotting
# tree_plot_data <- tree_structure %>%
#   left_join(mean_cluster_locs, by = c("parent" = "km_nn_cluster")) %>%
#   rename(from_x = TSNE1, from_y = TSNE2) %>%
#   left_join(mean_cluster_locs, by = c("child" = "km_nn_cluster")) %>%
#   rename(to_x = TSNE1, to_y = TSNE2)
# 
# 
# # transparent tsne of clusters
# tsne %>% 
#   ggplot(.,aes(TSNE1, TSNE2))+
#   scale_color_manual(values = cluster_cols)+ 
#   geom_point(aes(color = km_nn_cluster),size = 3, show.legend = FALSE)+
#   geom_segment(data = tree_plot_data, aes(x = from_x, y = from_y, xend = to_x, yend = to_y), size = 1, color = "white")+
#   geom_label(data = mean_cluster_locs, aes(label = km_nn_cluster, TSNE1, TSNE2),show.legend = FALSE,color = "black",  fill = "white", size = 10)+
#   theme(aspect.ratio = 1,
#     panel.background = element_rect(fill = "black", color = "black"),  
#     plot.background = element_rect(fill = "black", color = "black"),   
#     axis.text = element_blank(),                         
#     axis.title = element_blank(),                        
#     plot.title = element_blank(),                        
#     panel.grid = element_blank()                        
#   ) 
# ggsave("OneDrive - Lowy Medical Research Institute/Database_EDA/analysis_versions/version013/results/clustering/14_geneticProgressionModeling/black_bg_tree_label_cluster_tsne.png", height = 15, width = 15, dpi = 200)# transparent tsne
# 
# # transparent tsne of pseudotime
# tsne %>% 
#   ggplot(.,aes(TSNE1, TSNE2))+
#   geom_point(aes(color = pseudotime),size = 3, show.legend = FALSE)+
#   scale_color_viridis_c()+
#   theme(aspect.ratio = 1,
#     panel.background = element_rect(fill = "black", color = "black"),  
#     plot.background = element_rect(fill = "black", color = "black"),   
#     axis.text = element_blank(),                         
#     axis.title = element_blank(),                        
#     plot.title = element_blank(),                        
#     panel.grid = element_blank()                        
#   ) 
# ggsave("OneDrive - Lowy Medical Research Institute/Database_EDA/analysis_versions/version013/results/clustering/14_geneticProgressionModeling/black_bg_pt_tsne.png", height = 15, width = 15, dpi = 200)
# 
# 
# h <- colData(sce) %>% 
#   as.data.frame() %>% 
#   select(patID_eye, month, rf_route)
# 
# # par(mar = rep(5,4))
# # table of route categorized data points
# h %>%
#   mutate(rf_route =  ifelse(!is.na(rf_route), "Classified", "Unclassifiable")) %>% 
#   ggplot(aes(fill = rf_route, x = ''))+
#   geom_bar(stat = "count", position = "fill")+
#   ggtitle(paste("Proportion of data points \nclassified to route"))+
#   theme(axis.text = element_text(size = 20), 
#         axis.title = element_text(size = 25), 
#         title = element_text(size = 20), 
#         legend.text = element_text(size = 15))
#   theme_classic()
# 
# # table of route categorized eyes
# h %>% 
#   select(patID_eye, rf_route) %>% 
#   distinct() %>% 
#   mutate(rf_route =  ifelse(!is.na(rf_route), "Classified", "Unclassifiable")) %>% 
#   pull(rf_route) %>%
#   table(useNA = "ifany") %>% 
#   barplot(main = "Eyes with route" )
# # par(mfrow = c(1,1))





# A = c(6, 0, 18 ,4 ,10 )
# B = c(2, 0, 0 ,6 , 2 )
# C = c(16, 3 ,426 ,156, 34 )
# D = c(1, 3, 147, 434, 22 )
# X = c(6, 2, 44 ,35, 93 )
# 
# 
# tab <- matrix(c(A,B,C,D,X), byrow = TRUE, ncol = 5, 
#               dimnames= list(c(LETTERS[1:4], "X"),c(LETTERS[1:4], "X")))


# # RANDOM WALK ON GRAPH ----
# 
# edges <- which(transition_matrix>0, arr.ind = TRUE)
# g <- graph_from_edgelist(edges)
# 
# E(g)$weight <- transition_matrix[edges]
# set.seed(1312);plot(g, 
#                     edge.weight = E(g)$weight,
#                     layout = layout_as_tree(g, root = 1, circular = TRUE)
#                     )
# 
# P <- transition_matrix 
# 
# P <- P / rowSums(P)
# 
# P[is.nan(P)] <- 0
# 
# P[upper.tri(P)] <- P[upper.tri(P)]/ sum(P[upper.tri(P)])
# P[lower.tri(P)] <- P[lower.tri(P)]/ sum(P[lower.tri(P)])
# 
# my.map(P )
# 
# 
#    
# # Nsteps = 10; probability_matrix = P; start_node =1
# rwalk <-  function(probability_matrix, Nsteps = 5, start_node = NULL){
# 
#   # if no start node, make one
#   if(is.null(start_node)) start_node = sample(1:11,1)
# 
#   # initialize path
#   path= c(start_node) # start node
# 
#   # walk along nsteps times
#   for(i in 2:Nsteps){
#     current_node = path[[i-1]]
#     probs <- probability_matrix[current_node,]
# 
#     # end path if no next step
#     if(all(probs == 0)) break
# 
#     next_node = sample(1:11,1, prob =  probs)
#     path[[i]] <- next_node
#   }
# 
#   # make all same length
#   if(length(path) != Nsteps) path <- c(path,rep(NA, Nsteps - length(path)))
# 
#   return(path)
# }
# 
# # random walk on probability matrix
# iter <- 1.5e4
# Nsteps <- 10
# rwalks <- sapply(1:iter, function(i) rwalk(P, start_node = 1, Nsteps = Nsteps)) %>%
#   t()
# 
# # count of instances of each node for each step
# node_frequencies_at_each_step <- apply(rwalks, 2,table, simplify = FALSE)
# 
# # make matrix
# node_frequencies_at_each_step <- node_frequencies_at_each_step %>%
#   sapply(function(row){
#     row = as.matrix(row) %>% t()
#     if(ncol(row)!=11) row <- matrix(c(row, rep(NA, 11-ncol(row))), nrow = 1)
#     return(row)
#   })
# colnames(node_frequencies_at_each_step) <- paste("step",as.character(1:Nsteps), sep = "_")
# rownames(node_frequencies_at_each_step) <- paste("node",as.character(1:11), sep = "_")
# 
# # plot
# col_norm_node_frequencies_at_each_step <- node_frequencies_at_each_step%>%
#         sweep(., 2,colSums(., na.rm = TRUE), '/') 
# 
# col_norm_node_frequencies_at_each_step[is.na(col_norm_node_frequencies_at_each_step)] <- 0
# 
# col_norm_node_frequencies_at_each_step %>% 
#   # sweep(.,1,rowSums(., na.rm = TRUE), '/') %>% 
#   # sweep(., 2,colSums(., na.rm = TRUE), '/') %>%
#   sqrt() %>%
#   my.map


# COST FUNCTION TEST----
# cost_function <- function(matrix, target_row_sums, target_col_sums, distance_matrix) {
#   # Ensure the matrix is non-negative
#   matrix[matrix < 0] <- 0
#   
#   # Normalize the matrix so that row sums and column sums are 1
#   matrix <- sweep(matrix, 1, rowSums(matrix), "/") # Normalize row sums
#   matrix <- sweep(matrix, 2, colSums(matrix), "/") # Normalize column sums
#   
#   # Calculate row and column differences from the target
#   current_row_sums <- rowSums(matrix)
#   current_col_sums <- colSums(matrix)
#   
#   row_diff <- sum((current_row_sums - target_row_sums)^2)
#   col_diff <- sum((current_col_sums - target_col_sums)^2)
#   
#   # Fit the matrix to the distance matrix
#   fit_matrix <- 1 / matrix
#   fit_matrix[is.infinite(fit_matrix)] <- 0
#   
#   # Penalize matrices that don't fit the distance matrix
#   distance_penalty <- sum((matrix - distance_matrix)^2)
#   
#   # Return the total cost with penalties
#   return(row_diff + col_diff + distance_penalty)
# }
# 
# # cost_function(transition_matrix_prob_m, rep(1,11), rep(1,11), bif_distance_matrix)
# # cost_function(t(transition_matrix_prob_m), rep(1,11), rep(1,11), bif_distance_matrix)
# 
# # Example usage:
# result <- GenSA(par = as.numeric(a_matrix),
#                 fn = function(x) cost_function(matrix(x, nrow=nrow(a_matrix), ncol=ncol(a_matrix)),
#                                                rep(1, 11), rep(1, 11), bif_distance_matrix),
#                 lower = rep(0, length(as.numeric(a_matrix))),
#                 upper = rep(Inf, length(as.numeric(a_matrix)))) # ensure non-negativity
# 
# # Extract the resulting matrix
# annealed_matrix <- matrix(result$par, nrow=nrow(a_matrix), ncol=ncol(a_matrix))
# 
# # Normalize rows and columns to sum to 1
# annealed_matrix <- sweep(annealed_matrix, 1, rowSums(annealed_matrix), "/") 
# annealed_matrix <- sweep(annealed_matrix, 2, colSums(annealed_matrix), "/")



# #  LOOK AT AVERAGE PSEUDOTIME LINEAR PROGRESSION CLUSTER ORDERING ----
# cluster_pt <- sce %>% 
#      colData() %>% 
#      as.data.frame() %>% 
#      select(km_nn_cluster, pseudotime) %>% 
#      mutate(km_nn_cluster = factor(km_nn_cluster)) 
# 
# cluster_pt_means <- aggregate(pseudotime~ km_nn_cluster , mean, data = cluster_pt)
# cluster_pt_means <- cluster_pt_means[order(cluster_pt_means$pseudotime) , ]
# cluster_pt$km_nn_cluster <- factor(cluster_pt$km_nn_cluster, levels = cluster_pt_means$km_nn_cluster)
# boxplot(cluster_pt$pseudotime~x$km_nn_cluster, main = "Cluster Pseudotime ordered by mean")



 
# d <- colData(sce) %>% as.data.frame()
# d %>% 
#   select(patID, month, Normalised_PRS) %>% 
#   distinct() %>% 
#   filter(!is.na(Normalised_PRS)) %>% 
#   filter(n() > 1, .by = patID) %>% 
#   arrange(patID, month) %>% 
#   pull(patID) %>% 
#   n_distinct()


# # MAKE TREE DIAGRAM OF ROUTES ----
# library(igraph)
# # Define the edges according to the tree structure
# edges <- c(1, 2,
#            2,5,
#            2, 3, 
#            3, 11, 
#            3, 7, 
#            3, 9, 
#            2, 4, 
#            4, 10, 
#            4, 6, 
#            4, 8)
# 
# # Create the graph object
# g <- graph(edges = edges, directed = TRUE)
# 
# cols = c( "1" = "#3A6FA9", 
#           "2" = "#B2C6E5",
#           "3" = "#EE8635", 
#           "7" = "#A12D27",
#           "9" = "#77599B", 
#           "11" ="#6E4940", 
#           "4" = "#F2BF85",
#           "6" = "#A8DD93",
#           "8" = "#F09D99",
#           "10" ="#C1B0D3",
#           "5" = "#509E3E"
# )
# par(bg = "black")
# set.seed(12433)
# plot(g, 
#      layout = layout_as_tree(g, root = 1),  # Set 1 as the root
#      vertex.size = 30,                     # Size of nodes
#      vertex.label.cex = 1.5,               # Text size
#      vertex.label.color = "black",         # Label color
#      vertex.color = cols,
#      edge.arrow.size = 0.5,                # Arrow size
#      main = "Tree Diagram")   
# par(bg  = "white")








# # LOOK AT FDR CORRECTION ----
# ps <- at_quintile_PT_interaction_snp_results$`Pr(>|t|)` %>% 
#   sort() %>% 
#   cbind(1:length(at_quintile_PT_interaction_snp_results$`Pr(>|t|)`), .) %>% 
#   as.data.frame() %>% 
#   rename_with(~c("rank", "P"))
# 
# a <- 0.05
# ps$crit <- (ps$rank / nrow(ps)) * a
# ps$reject_null <- ps$P < ps$crit
# 
# # Calculate raw Q values
# ps$Q_raw <- (nrow(ps) * ps$P) / ps$rank
# 
# # Ensure Q values are monotonic (non-increasing)
# ps <- ps %>% 
#   arrange(desc(rank)) %>% 
#   mutate(Q = cummin(Q_raw)) %>% 
#   arrange(rank)
# 
# # Plot P vs Q
# f <- function(n) -log10(n)
# par(mar = c(4,4,2,2))
# plot(ps$P %>% f, ps$Q %>% f, main = "P-values vs Q-values", xlab = "P-values", ylab = "Q-values")
# abline(0, 1)

# # MAKE SLOPE QUINTILE PLOT AGAIN ----
# # uses impupted data so cleaner looking
# sce <- HDF5Array::loadHDF5SummarizedExperiment("OneDrive - Lowy Medical Research Institute/Database_EDA/analysis_versions/version013/processed_data/sce006.h5")
# 
# # select data and calculate parafovea_mean_excl_3
# slope_d <-
#   colData(sce) %>%
#   as.data.frame() %>%
#   select(patID, eye_num, patID_eye, subject, pseudotime, month) %>%
#   left_join(
#     assays(sce)$pmm_imp_of_scaled %>%
#       t() %>%
#       as.data.frame() %>%
#       rownames_to_column("subject") %>%
#       select(subject, matches("avg_MT_slope"))
#   ) 
# # make BL_parafovea_quintile
# slope_quantile_d <- slope_d %>%
#   arrange(month, patID_eye) %>%
#   group_by(patID_eye) %>%
#   filter(row_number() ==1) %>%
#   ungroup() %>%
#   mutate(
#     slope_quintile =
#       cut(avg_MT_slope_sans3, quantile(avg_MT_slope_sans3, seq(0,1,.2)), labels = 1:5)
#   ) %>%
#   select(patID_eye, slope_quintile) %>%
#   distinct()
# 
# slope_d <- slope_d %>%
#   left_join(slope_quantile_d)
# 
# # rescale parafovea to distr of mean/sd(avt_mt02,4,5)
# AVG_SD <- metadata(sce)$d %>%
#   select(avg_MT_slope_sans3) %>%
#   mutate(rowmean= rowMeans(across(everything()))) %>%
#   mutate(avg = mean(rowmean, na.rm = TRUE)) %>%
#   mutate(sd = sd(rowmean, na.rm = TRUE)) %>%
#   select(avg, sd) %>%
#   distinct
# 
# ## rescale
# slope_d$avg_MT_slope_sans3 <- (slope_d$avg_MT_slope_sans3 * AVG_SD$sd) + AVG_SD$avg
# 
# 
# 
# # plot parafovea over time faceted for quintile
# ggplot(slope_d %>% drop_na(), aes(month, avg_MT_slope_sans3))+
#   # geom_hline(yintercept = 0, linetype = 3)+
#   geom_smooth(method= "lm", aes(group = 1))+
#   facet_wrap(~factor(slope_quintile, levels = 5:1),ncol = 5, scales = "free_x")+
#   ggtitle("Average parafoveal slope (excluding zone 3) over time quantiled by BL thickness")+
#   theme_bw()
# ggsave("OneDrive - Lowy Medical Research Institute/Database_EDA/analysis_versions/version013/results/longitudinal_analysis/parafoveal_slope_mean_quintile_plot(imputedMTfromsce).png", height = 5, width = 20)

# # MAKE PARAFOVEA QUINTILE PLOT AGAIN ----
# # uses impupted data so cleaner looking
# sce <- HDF5Array::loadHDF5SummarizedExperiment("OneDrive - Lowy Medical Research Institute/Database_EDA/analysis_versions/version013/processed_data/sce006.h5")
# 
# # select data and calculate parafovea_mean_excl_3
# parafovea <-
#   colData(sce) %>%
#   as.data.frame() %>%
#   select(patID, eye_num, patID_eye, subject, pseudotime, month, FA1, FA2) %>%
#   left_join(
#     assays(sce)$pmm_imp_of_scaled %>%
#       t() %>%
#       as.data.frame() %>%
#       rownames_to_column("subject") %>%
#       select(subject, matches("avg_MT0[245]$"))
#   ) %>%
#   mutate(parafovea_excl_3 = rowMeans(across(matches("avg_MT"))))
# 
# # make BL_parafovea_quintile
# parafovea_quantile <- parafovea %>%
#   arrange(month, patID_eye) %>%
#   group_by(patID_eye) %>%
#   filter(row_number() ==1) %>%
#   ungroup() %>%
#   mutate(
#     parafovea_quintile =
#       cut(parafovea_excl_3, quantile(parafovea_excl_3, seq(0,1,.2)), labels = 1:5)
#   ) %>%
#   select(patID_eye, parafovea_quintile) %>%
#   distinct()
# 
# parafovea <- parafovea %>%
#   left_join(parafovea_quantile)
# 
# # rescale parafovea to distr of mean/sd(avt_mt02,4,5)
# AVG_SD <- metadata(sce)$d %>%
#   select(matches("avg_MT0[245]$")) %>%
#   mutate(rowmean= rowMeans(across(everything()))) %>%
#   mutate(avg = mean(rowmean, na.rm = TRUE)) %>%
#   mutate(sd = sd(rowmean, na.rm = TRUE)) %>%
#   select(avg, sd) %>%
#   distinct
# ## rescale
# parafovea$parafovea_excl_3 <- (parafovea$parafovea_excl_3 * AVG_SD$sd) + AVG_SD$avg
# 
# 
# # plot parafovea over time faceted for quintile
# ggplot(parafovea %>% drop_na(), aes(month, parafovea_excl_3))+
#   # geom_hline(yintercept = 0, linetype = 3)+
#   geom_smooth(method= "lm", aes(group = 1))+
#   facet_wrap(~factor(parafovea_quintile, levels = 5:1),ncol = 5, scales = "free_x")+
#   ggtitle("Average parafoveal macular thickness over time quantiled by BL thickness")+
#   theme_bw()
# ggsave("OneDrive - Lowy Medical Research Institute/Database_EDA/analysis_versions/version013/results/longitudinal_analysis/parafoveal_mean_quintile_plot(imputedMTfromsce).png", height = 5, width = 20)


### -----
# # rm(list = ls()); gc()
# library(HDF5Array)
# library(SingleCellExperiment)
# library(dplyr)
# library(ggplot2)
# 
# # Load the SummarizedExperiment object
# sce <- HDF5Array::loadHDF5SummarizedExperiment("OneDrive - Lowy Medical Research Institute/Database_EDA/analysis_versions/version013/processed_data/sce006.h5")
# 
# # Extract the assays and metadata
# d <- assays(sce)$measure_vars %>% 
#   t() %>% 
#   as.data.frame() %>% 
#   rownames_to_column("subject") %>% 
#   left_join(
#     as.data.frame(colData(sce)),
#     by = "subject"
#   )
# 
# # Select relevant columns
# slope_d <- d %>% 
#   select(subject, month, patID, patID_eye, avg_MT_slope_sans3)
# 
# # Filter distinct rows and group by patID_eye
# slope_q_d <- slope_d %>% 
#   arrange(month) %>% 
#   group_by(patID_eye) %>% 
#   filter(n() == 1) %>% 
#   ungroup() %>% 
#   select(patID_eye, avg_MT_slope_sans3) %>% 
#   distinct()
# 
# # Create slope quantile bins
# slope_q_d <- slope_q_d %>% 
#   mutate(slope_q = cut(avg_MT_slope_sans3, 
#                        breaks = quantile(avg_MT_slope_sans3, seq(0, 1, .2), na.rm = TRUE), 
#                        labels = 1:5, include.lowest = TRUE)) %>% 
#   select(-avg_MT_slope_sans3) %>% 
#   distinct()
# 
# # Join quantile bins with the main dataset
# slope_d <- slope_d %>% 
#   left_join(slope_q_d, by = "patID_eye")
# 
# # Plot the data
# slope_d %>% 
#   ggplot(aes(x = month, y = avg_MT_slope_sans3)) +
#   geom_smooth(method = "lm") +
#   facet_wrap(~slope_q)


# # See if i can use lm to predict route ----
# 
# # rm(list = ls()); gc()
# # load data
# if(!"sce" %in% ls()) sce <-  HDF5Array::loadHDF5SummarizedExperiment("OneDrive - Lowy Medical Research Institute/Database_EDA/analysis_versions/version013/processed_data/sce006.h5")
# 
# # get tsne withi ids
# tsne <- as.data.frame(reducedDim(sce, "TSNE")) %>%
#   rownames_to_column("subject") %>%
#   left_join(
#     as.data.frame(colData(sce)) %>%
#       select(subject, patID, patID_eye, month)
#   ) %>%
#   arrange(patID_eye, month)
# 
# 
# 
# an_eye= "001008_2"
# eye_tsne <- tsne[ tsne$patID_eye == an_eye, ]
# 
# # look at lm line across tsne change
# plot(tsne[, c("TSNE1", "TSNE2")], col = "grey", type = "n", main = an_eye)
# abline(a = 0, b = 0, lwd = 1.5, lty = 2); abline(v = 0, lwd = 1.5, lty = 2)
# points(tsne[, c("TSNE1", "TSNE2")], col = "grey")
# lines(eye_tsne[, c("TSNE1", "TSNE2")], col = "blue", lty = 3, lwd = 6, type = "b")
# abline(lm(TSNE2~TSNE1, eye_tsne), col = "red", lwd = 2)
# 
# 
# # normalize tsnes to center around cluster 2
# tsne1_clust_2_mean <- reducedDim(sce[, sce$km_nn_cluster == 2], "TSNE")[,1] %>% mean()
# tsne2_clust_2_mean <- reducedDim(sce[, sce$km_nn_cluster == 2], "TSNE")[,2] %>% mean()
# norm_tsne <- tsne
# norm_tsne$TSNE1 <- norm_tsne$TSNE1  - tsne1_clust_2_mean
# norm_tsne$TSNE2 <- norm_tsne$TSNE2 - tsne2_clust_2_mean
# 
# # lm for all eyes
# eye_tsne_results <- data.frame();for(i in 1:length(unique(norm_tsne$patID_eye))){
#   eye = unique(norm_tsne$patID_eye)[i]
#   a_result <- broom::tidy(summary(lm(TSNE2~TSNE1, norm_tsne[norm_tsne$patID_eye == eye, ])))
#   a_result$patID_eye <- eye
#   if(i == 1) eye_tsne_results <- a_result else eye_tsne_results <- rbind(eye_tsne_results, a_result)
# }
# 
# # clean results
# eye_tsne_results_wide <- eye_tsne_results %>%
#   mutate(term = case_when(
#     term == "(Intercept)" ~ "intercept",
#     TRUE ~ term
#   )) %>%
#   pivot_wider(names_from = term,
#               id_cols = patID_eye,
#               values_from = c(p.value,estimate ) )
# 
# # filter for sig models; grab those eyes
# eyes <- eye_tsne_results_wide %>%
#   filter(if_all(matches("p.value"), ~ .< .005)) %>%
#   filter((estimate_intercept> 0 & estimate_intercept< 100) &
#            (estimate_TSNE1 > 0 & estimate_TSNE1 < 20)) %>%
#   # filter(`(Intercept)` >0 & TSNE1 > 0) %>%
#   pull(patID_eye) %>%
#   sample(20)
# 
# # plot eyes with line of lm
# par(mfrow = c(2,10));for(an_eye in eyes){
#   eye_tsne <- tsne[ tsne$patID_eye == an_eye, ]
#   plot(tsne[, c("TSNE1", "TSNE2")], col = "grey", type = "n", main = an_eye)
#   abline(a = 0, b = 0, lwd = 1.5, lty = 2); abline(v = 0, lwd = 1.5, lty = 2)
#   points(tsne[, c("TSNE1", "TSNE2")], col = "grey")
#   lines(eye_tsne[, c("TSNE1", "TSNE2")], col = "blue", lty = 3, lwd = 6, type = "b")
#   abline(lm(TSNE2~TSNE1, eye_tsne), col = "red", lwd = 2)
#   # lines(eye_tsne$TSNE1,predict(loess(TSNE2~TSNE1, eye_tsne)), col = "red", lwd = 2)
# };par(mfrow = c(1,1))
# 
# route_data <- as.data.frame(colData(sce)) %>%
#   select(matches("patID|month"), rf_route) %>%
#   left_join(
#     eye_tsne_results_wide
#   ) %>%
#   drop_na()
# 
# truth <- route_data %>%
#   select(-matches("patID|month"))
# train_indices <- createDataPartition(truth$rf_route, p = 0.6, list = FALSE)
# 
# train <- truth[train_indices, ]
# test <- truth[-train_indices , ]
# 
# rf <- randomForest::randomForest(rf_route ~ .,data = train)
# rf2 <- rpart::rpart(rf_route ~ .,
#                     data = train, 
#                     weights = 1+c("A" =2, "B"=2, "C"=1, "D"=1, "X"=0)[match( train$rf_route,c("A", "B", "C", "D", "X"))]
#                     )
# rpart.plot::rpart.plot(rf2, cex = .7)
# 
# calc_acc <- function(sq_m) sum(diag(sq_m))/ sum(sq_m)
# 
# table(train$rf_route, predict(rf, train, type = "class")) %>% calc_acc
# table(train$rf_route, predict(rf2, train, type = "class"))%>% calc_acc
# 
# 
# route_parameter_means <- route_data %>%
#   # filter(if_any(matches("p.val"), ~ . <.01)) %>%
#   select(rf_route, matches("p.val|estim")) %>%
#   group_by(rf_route) %>%
#   summarize_all(~median(.))
# 
# cols = viridis(5)
# for(i in 1:length(route_parameter_means$rf_route)){
#   route = route_parameter_means$rf_route[i]
#   intercept = route_parameter_means$estimate_intercept[route_parameter_means$rf_route==route]
#   estimate = route_parameter_means$estimate_TSNE1[route_parameter_means$rf_route==route]
#   p = route_parameter_means$p.value_TSNE1[route_parameter_means$rf_route==route]
# 
#   if(i == 1) plot(-100:100, -100:100, type = "n")
#   abline(c(intercept, estimate), col = cols[i])
#   text(i,intercept, paste(route, formatC(p, format = "e", digits = 3), sep = ":"))
# };



# ----
# for(i in 1:length(dfs)){
#    df <- dfs[[i]]
#    df_name = names(dfs)[i]
#    cat(paste0(df_name, ": ",paste(dim(df), collapse = ","), "\n"))
# }

# x <- 
#   # matrix(rnorm(1000, mean = 10, sd = 5), ncol = 2) %>% 
#   matrix(sample(runif(1000,min = 0, max = 10), 1000,FALSE), ncol = 2) #%>% round
# colnames(x) <- c("score1", "score2")
# # x <- x[ x[,1] != x[,2]  ,]
# 
# x <- cbind(x, score1minus2 = x[,"score1"] - x[, "score2"])
# x <- cbind(x, normscore1minus2 = (x[,"score1minus2"]) / ((x[,"score1"] + x[, "score2"])))
# 
# plot(data.frame(x))

# PC <- reducedDims(sce)$PCA
# 
# scores_df <- as.data.frame(PC)
# 
# # Extract the first two principal components
# scores_df <- scores_df[, 1:2]
# colnames(scores_df) <- c("PC1", "PC2")
# 
# # Create a biplot using ggplot2 and ggfortify
# p <- ggplot(scores_df, aes(x = PC1, y = PC2)) +
#   geom_point(alpha = 0.5) +  # Plot the principal component scores
#   labs(title = "PCA Biplot", x = "PC1", y = "PC2") +
#   theme_minimal()
# 
# # Add the loadings (rotation) to the plot
# loadings <- attr(PC, "rotation")
# loadings_df <- as.data.frame(loadings[, 1:2])
# colnames(loadings_df) <- c("PC1", "PC2")
# 
# # Scale the loadings for better visualization
# loadings_df <- loadings_df * 20  # Adjust scaling factor as needed
# 
# # Add arrows for the loadings
# p <- p +
#   geom_segment(data = loadings_df, aes(x = 0, y = 0, xend = PC1, yend = PC2),
#                arrow = arrow(length = unit(0.2, "cm")), color = "red") +
#   geom_text(data = loadings_df, aes(x = PC1, y = PC2, label = rownames(loadings_df)),
#             color = "red", vjust = 1, hjust = 1)
# 
# p
# PC <- reducedDims(sce)$PCA
# pcs <- PC %>%
#   as.data.frame() %>%
#   rownames_to_column("subject") %>%
#   left_join(d %>% select(subject, patID_eye, patID))
# 
# 
# # grouping <- "patID_eye"
# grouping <- "patID"
# set.seed(2330)
# pcs <- pcs[pcs[[grouping]] %in% sample(pcs[[grouping]], 100), ]
# 
# cols <- RColorBrewer::brewer.pal(n = 10,name = "Set1")
# cols <- rep(cols, n_distinct(pcs[[grouping]]))
# cols <- cols[1:n_distinct(pcs[[grouping]])]
# 
# names(cols) <- unique(pcs[[grouping]])
# 
# plot(factor(pcs[[grouping]]), 
#      rowMeans(pcs[, c("PC1", "PC2")]),
#      # pcs$PC1,
#      col = cols[match(pcs[[grouping]], names(cols))],
#      pch = 16, 
#      main  = paste("grouping =", grouping)
#      )


# set.seed(3939)
# a_subject <- sample(sce$subject, 1)
# 
# TSNE <- reducedDims(sce)$TSNE 
# PC <- reducedDims(sce)$PCA  
# an_eye_tsnes <- TSNE[ rownames(TSNE) == a_subject , ]
# an_eye_pcs <- PC[ rownames(PC) == a_subject , ]
# rotation <- attributes(reducedDims(sce)$PCA)$rotation
# 
# assay <- assays(sce)$pmm_imp_of_scaled %>% 
#   t() 
# 
# recomposed_subject_values <- an_eye_pcs %*% t( rotation)
# 
# dist(
#   rbind(assay[ rownames(assay) == a_subject,], recomposed_subject_values)
#          )




# df <- assays(sce)$pmm_imp_of_scaled %>% 
#   t() %>% 
#   as.data.frame()
# 
# mod_list <- list(); for(i in 1:ncol(df)){
#   
#   mod <- lm(sce$pseudotime ~ df[,i])
#   
#   mod_list[[names(df)[i]]] <- mod
#   
# }
# 
# summs <- lapply(mod_list, summary)
# 
# results <- data.frame(
#   term = character(), 
#   estimate = numeric(), 
#   pr = numeric()
# )
# 
# # term <- names(summs)[1]
# for(term in names(summs)){
#   
#   coefs <- coef(summs[[term]])
#   
#   a_result <- data.frame(
#     term = c(term), 
#     estimate = coefs[2, "Estimate"], 
#     pr = coefs[, "Pr(>|t|)"], 
#     row.names = NULL
#   )
#   results <- rbind(results, a_result)
#   
# }
# 
# results$q <- p.adjust(results$pr, method =  "BH")
# 
# results %>% 
#   filter(q>0) %>% 
#   ggplot(aes(estimate, -log10(q)))+
#   geom_point()+
#   geom_hline(yintercept = -log10(.05), linetype = 3)+
#   geom_label_repel(aes(label = ifelse(q<.05, term, NA)), size = 5)

# # ratio of late time eyes to total of cluster
# 
# # wd_path = "/Users/christiananderson/Library/CloudStorage/OneDrive-LowyMedicalResearchInstitute/Database_EDA/analysis_versions/version013/"
# # # load experiment
# # sce <- readRDS(paste0(wd_path, "processed_data/sce003.rds"))
# 
# dat <- colData(sce) %>% 
#   as.data.frame() %>% 
#   mutate(month = round(month)) %>% 
#   select(patID_eye, month, km_nn_cluster)
# 
# round_to_year <- function(a_var){
#   cutted <- cut(a_var, 
#     breaks = seq(0, max(a_var)+12, 12), 
#     include.lowest = TRUE, 
#     labels = FALSE,
#     right = FALSE
#     )
#   
#   return(cutted-1)
# }
# # x <- 1:50
# # cbind(x, round_to_year(x))
# 
# dat <- dat %>% 
#   mutate(month = round(month)) %>% 
#   mutate(year = round_to_year(month)) 
# 
# cluster_sizes <-  dat %>% 
#        group_by(km_nn_cluster) %>% 
#        summarize(size = n())
# 
# count_per_month <- dat %>% 
#   group_by(month) %>% 
#   summarize(size = n())
# 
# month_counts <- dat %>% 
#   group_by(km_nn_cluster,month) %>% 
#   summarize(month_count = n()) %>% 
#   ungroup() %>% 
#   left_join(count_per_month) %>% 
#   mutate(ratio = month_count/ size)
#   
# ggplot(dat, aes(x = month, fill = km_nn_cluster, color = km_nn_cluster))+
#   stat_density( kernel = "gaussian", alpha = .01, position = "dodge", linewidth = 1)
#   # stat_density( kernel = "gaussian", alpha  = .4)
# 
# 
# month_counts %>%
#   group_by(km_nn_cluster) %>% 
#   mutate(ratio = ratio - mean(ratio)) %>% 
# ggplot(aes(y = ratio, x = month,  color = km_nn_cluster  , linetype = km_nn_cluster))+
#     geom_smooth(se = F, method = "loess")

## MDS clustering----
# wd_path = "/Users/christiananderson/Library/CloudStorage/OneDrive-LowyMedicalResearchInstitute/Database_EDA/analysis_versions/version013/"
# 
# # load experiment
# sce <- readRDS(paste0(wd_path, "processed_data/sce003.rds"))
# 
# # mds
# x <- assays(sce)$pmm_imp_of_scaled 
# 
# x_dist <- x %>% 
#   dist() 
# 
# x_mds <- cmdscale(x_dist)
# 
# # x_pc <- prcomp(x)
# 
# plot(x_mds, type = "n")
# text(x_mds, label = rownames(x_mds))
# 
# 
# # nmf
# m <- x %>% 
#   t() 
# 
# m_sub <- m[, grep("RW[10]", colnames(m))]
# 
# glimpse(m_sub %>% as.data.frame())
# 
# # Perform NMF
# nmf_result <- nmf(m_sub, rank = 4, method = "brunet", nrun = 50)
# 
# # Extract W and H
# W <- basis(nmf_result)
# H <- coef(nmf_result)
# 
# plot(W)


# ##### test modeling of avg_MT_slope_sans3_rate ----
# sce <- readRDS("/Users/christiananderson/Library/CloudStorage/OneDrive-LowyMedicalResearchInstitute/Database_EDA/analysis_versions/version013/processed_data/sce003.rds")
# 
# data <- metadata(sce)$d %>%
#   as.data.frame()
# 
# plot(density(d$avg_MT_slope_sans3, na.rm = TRUE))
# plot(d$avg_MT_slope_sans3~d$ACU_quintile)
# 
# 
# 
# vars <- sapply(data[, sapply(data, is.numeric)], function (col){
# 
#   (sd(col, na.rm = TRUE) !=0) & (n_distinct(col[!is.na(data$avg_MT_slope_sans3_rate)])>1)
# 
# })
# 
# rates <- data %>%
#   select(patID_eye, month, avg_MT_slope_sans3) %>%
#   arrange(month) %>%
#   group_by(patID_eye) %>%
#   mutate(avg_MT_slope_sans3_rate =
#            (dplyr::last(avg_MT_slope_sans3[!is.na(avg_MT_slope_sans3)]) - dplyr::first(avg_MT_slope_sans3[!is.na(avg_MT_slope_sans3)]))
#             /
#            (dplyr::last(month[!is.na(month)]) - dplyr::first(month[!is.na(month)]))
#   ) %>%
#   ungroup() %>%
#   select(-month, -avg_MT_slope_sans3) %>%
#   na.omit() %>%
#   unique()
# 
# mod_data <- data %>%
#   select(patID_eye, patID, month, all_of(names(vars)[vars] ) ) %>%
#   mutate(across(all_of(names(vars)[vars] ) , ~as.numeric(scale(.)))) %>%
#   select(- avg_MT_slope_sans3_rate) %>%
#   left_join(rates)
# 
# 
# 
# 
# mod_data_first <- mod_data %>%
#   arrange(month) %>%
#   group_by(patID_eye) %>%
#   mutate(across(
#       -c( patID, SEX, month, avg_MT_slope_sans3_rate),
#     ~ dplyr::first(.[!is.na(.)]),
#     .names = "first_{col}"
#   )) %>%
#   select(-month) %>%
#   unique()
# 
# mod_data_first <-   mod_data_first %>%
#   select(patID_eye, patID, SEX, avg_MT_slope_sans3_rate, starts_with("first")) %>%
#   unique()
# 
# mods <- sapply(names(mod_data_first)[4:length(mod_data_first)], function(var){
# 
#   a_formula <- as.formula(paste0(  "avg_MT_slope_sans3_rate ~ SEX + ", var, " + (1|patID)"))
#   a_mod <- lmer( a_formula, data = mod_data_first)
# 
# })
# 
# summs <- lapply(mods, summary)
# 
# results <- data.frame(Model = character(),  Term = character(),Estimate = numeric(),
#                          Std.Error = numeric(),t.value = numeric(),Pr = numeric(),stringsAsFactors = FALSE)
# 
# for (model_name in names(summs)) {
# 
#   model_coef <- coef(summs[[model_name]])
# 
#   a_result_df <- data.frame(Model = model_name,
#                             Term = rownames(model_coef),
#                             Estimate = model_coef[, "Estimate"],
#                             Std.Error = model_coef[, "Std. Error"],
#                             t.value = model_coef[, "t value"],
#                             Pr = model_coef[, "Pr(>|t|)"])
# 
#   results <- rbind(results, a_result_df)
# }
# 
# 
# results <- results %>%
#   filter(!Term %in% c("(Intercept)", "SEX") )
# 
# results_filt <-  results %>%
#   filter(!Term %in% grep("RW|THICK|CENTER|SNELL|LOG|chew|letters|FA_|FA\\d|ACUIT|MT|slope",results$Term, value = TRUE))
# 
# ggplot(results_filt, aes(x = Estimate, y = -log10(Pr) , color = Model))+
#   geom_point( show.legend = FALSE)+
#   geom_hline(yintercept = -log10(.05), linetype = 2, color ="blue")+
#   geom_label_repel(aes(label = ifelse(Pr< .05, Model, NA )  ) ,
#                    max.overlaps = 20,
#                    show.legend = FALSE,
#                    nudge = 3
#                    )+
#   geom_vline(xintercept = 0, linetype = 2, color ="grey")

#### ----

  

# # Plot the scatter plot
# plot(rnorm( length(skw)), skw, type = "n")
# text(rnorm( length(skw)), skw, names(skw), col = "blue", cex = 1.2, srt = 90)
# abline(a = mean(skw), b = 0, lty = 2)


# # Sample data
# x <- 1:10
# y1 <- x^2
# y2 <- x^3
# 
# # Plot the first line
# plot(x, y1, type = "l", col = "blue", ylim = range(c(y1, y2)), xlab = "X-axis", ylab = "Y-axis", main = "Plot with Legend")
# 
# # Add the second line
# lines(x, y2, col = "red")
# 
# # Add a legend
# legend("topright", legend = c("y = x^2", "y = x^3"), col = c("blue", "red"), lty = 1)
# 


# wd_path = "/Users/christiananderson/Library/CloudStorage/OneDrive-LowyMedicalResearchInstitute/Database_EDA/analysis_versions/version012/"
# 
# data_name <- "CF_FFA_DEM_DIAGNOSIS_ENR_FAFGRADE_FAMSTATUS_MH_OCTMAIN_OCTMTCIRRUS_OCTMTSPECTRALIS_OCTMTSTRATUS_VAA_02.rds"
# 
# # joined data full
# d012 <- readRDS(paste0(wd_path,"processed_data/",data_name))
# 
# d012 %>% select(patID, eye_num) %>% unique() %>% nrow()
#  
# d012 %>% 
#   group_by(patID, eye_num, visitID) %>% 
#   filter(n()>1)


# d %>% 
#   group_by(patID, eye_num) %>% 
#   filter(visitID == first(visitID)) %>% 
#   ungroup() %>% 
#   group_by(visitID) %>% 
#   summarize(num_vis = n()) %>% 
#   ungroup()


# grading_diagnosis_info_joined %>%
#   select(patientId.x,eye.x, visitId.x, DATE_period_year.x) %>% 
#   left_join(., new_dfs$VAA , by = join_by(patientId.x == patientId, DATE_period_year.x == DATE_period_year, visitId.x == visitId))


# # Define the heatmap
# m <- d %>% 
#   mutate_all(~as.numeric(!is.na(.))) %>% 
#   as.matrix() 
# 
# 
# unique_ids <- unique(d$visitID) %>% 
#   factor(., levels = sort(.))
# 
# names <- brewer.pal(length(unique_ids), "Paired")
# 
# # Create a named vector
# # named_ids <- setNames(d$visitID, names[match(d$visitID, unique_ids)])
# named_ids <- setNames( names,unique_ids)
# 
# # Create a row annotation with colors
# row_annotation <- rowAnnotation(Identity = d$visitID, col = list(Identity = named_ids))
# 
# # Create and draw the heatmap with row annotations
# heatmap <- Heatmap(m, 
#                    right_annotation = row_annotation, 
#                    show_row_names = FALSE, 
#                    cluster_rows = FALSE, 
#                    cluster_columns  = FALSE, 
#                    col = c("white", "black"))
# 
# draw(heatmap)



# d <- new_dfs$CF_FFA
# 
# d <- d %>% 
#   filter(!is.na(DATE) & !is.na(eye)) %>% 
#   mutate(pat_eye_DATE = paste(patientId, eye, DATE, sep  = "_")) %>% 
#   mutate(pat_eye_DATE_period_year = paste(patientId, eye, DATE_period_year, sep  = "_")) %>% 
#   mutate(pat_eye_visitId = paste(patientId, eye, visitId, sep  = "_")) 
# 
# va <- new_dfs$VAA %>% 
#   mutate(pat_eye_DATE_period_year = paste(patientId, eye, DATE_period_year, sep  = "_")) %>% 
#   mutate(pat_eye_visitId = paste(patientId, eye, visitId, sep  = "_")) %>% 
#   filter(!is.na(DATE) & !is.na(eye))
#   
# 
# d$pat_eye_visitId %>% length()
# va$pat_eye_visitId%>% length()
# 
# setdiff(
# d$pat_eye_visitId %>% unique(),
# va$pat_eye_visitId%>% unique()
# )%>% length()
# setdiff(
# d$pat_eye_DATE_period_year %>% unique(),
# va$pat_eye_DATE_period_year%>% unique()
# ) %>% length()
# setdiff(
# va$pat_eye_visitId%>% unique(),
# d$pat_eye_visitId%>% unique()
# )%>% length()

 # dfs$CF_FFA %>% 
 #  select(matches("PHOTODT")) %>%
 #  filter(if_any(everything(), ~!is.na(.))) %>% 
 #  mutate(diff = difftime(.[[1]], .[[2]], units = "days") %>% as.numeric()) %>% 
 #  filter(abs(diff) > 100)


 
# raw$VAA %>% 
#   select(patientId, eye, visitId, TCVISIT_new, VAA_DATE, VAA_DATE_period_year) %>% 
#   arrange(patientId, VAA_DATE_period_year)

# d %>% 
#   select(matches("period_year")) %>% 
#   filter(spectralis_SCANDT_period_year != SCANDT_period_year) %>% 
#   select(spectralis_SCANDT_period_year, SCANDT_period_year)

  
# ans <- sapply(1:(1*1e6), function(i){
# 
#   return(as.numeric(paste0(i,i,i)) / (i+i+i))
# 
# }, USE.NAMES =  FALSE)
# 
# unique(ans)
# plot(cbind(1:length(ans) , log10(ans) ))
# 
# ans_d <- data.frame(ans = ans, log10_ans= log10(ans)) %>%
#   mutate(n = 1:n()) %>%
#   unique() %>%
#   select(3,1) %>%
#   filter(ans != max(ans))



# d <- readRDS("/Users/christiananderson/Library/CloudStorage/OneDrive-LowyMedicalResearchInstitute/Database_EDA/analysis_versions/version012/processed_data/CF_FFA_DEM_DIAGNOSIS_ENR_FAFGRADE_FAMSTATUS_MH_OCTMAIN_OCTMTCIRRUS_OCTMTSPECTRALIS_OCTMTSTRATUS_VAA_07.rds")
# 
# pdf("/Users/christiananderson/Library/CloudStorage/OneDrive-LowyMedicalResearchInstitute/Database_EDA/analysis_versions/version012/results/VAA_analysis/FA_VA_vars_plot.pdf", height = 30, width =30)
# plot(d %>% select(matches("FA\\d$"),FA_mean, ACUIT, SNELLA, LOGMAR))
# dev.off()
# 
# 
# ggplot(d, aes(x = month, y = FA4))+
#   geom_smooth(method = "lm")+
#   facet_wrap(~FA4_quintile, scales = "free")
# ggplot(d, aes(x = month, y = FA3))+
#   geom_smooth(method = "lm")+
#   facet_wrap(~FA3_quintile, scales = "free")
# ggplot(d, aes(x = month, y = FA2))+
#   geom_smooth(method = "lm")+
#   facet_wrap(~FA2_quintile, scales = "fixed", ncol = 5)

# d %>% 
#   summarize(n_pats = n_distinct(patID), 
#             mean_num_visits = mean(as.numeric(visitID), na.rm =  TRUE),
#             median_num_visits = median(as.numeric(visitID), na.rm =  TRUE)
#             )


# x <- raw$VAA %>% select(matches("vis")) %>% 
#   na.omit %>%
#   mutate(visitID.old_eq_TCVISIT.new = visitId_old == TCVISIT_new, 
#          
#          ) 
# sum(x$visitID.old_eq_TCVISIT.new)
# sum(!x$visitID.old_eq_TCVISIT.new )
#







# l <-  grepl("OCT|BLR", ordering_df$variable); ordering_df$variable[l]

# filtering_d <- d %>%
#   arrange(patID_eye) %>%
#   mutate(across(everything(), ~as.numeric(as.character(.))))
# 
# eye1 <- filtering_d %>%
#   filter(eye_num == 1)
# eye2 <- filtering_d %>%
#   filter(eye_num == 2)
# 
# eye1$eye_num <- NULL
# eye2$eye_num <- NULL
# 
# names(eye1)[3:ncol(eye1)] <- paste(names(eye1)[3:ncol(eye1)], "1", sep = "_")
# names(eye2)[3:ncol(eye2)] <- paste(names(eye2)[3:ncol(eye2)], "2", sep = "_")
# 
# wide <- left_join(eye1, eye2)
# 
# 
# include_vars <- sapply(names(sub)[4:ncol(sub)], function(column){
#   left <- paste0(column, "_2")
#   right <- paste0(column, "_1")
# 
#   var_left <- var(wide[[left]], na.rm = TRUE) %>% as.numeric()
#   var_right <- var(wide[right], na.rm = TRUE)%>% as.numeric()
# 
#   if (var_left == 0 | var_right == 0 | all(wide[[right]] == wide[[left]], na.rm = TRUE)) {
#         return(FALSE)
#     } else {
#         return(TRUE)
#     }
# 
# })
# 
# sub <- d %>% 
#   select(all_of(names(include_vars))) %>%
#   select(where(~is.numeric(.))) %>% 
#   select(where(~sd(., na.rm = TRUE) > 0 && var(., na.rm =TRUE) > 0)) 
# 
# sub_t <- sub %>% 
#   t() %>% 
#   apply(., 2, as.numeric)
# 
# missing_rows <- rowSums(is.na(sub_t)) == ncol(sub_t)
# 
# sub_t_sub <- sub_t[!missing_rows , ]
# sub_t_sub_cor <- cor(sub_t_sub, use = "pairwise.complete.obs" )
# 
# Heatmap( sub_t_sub_cor, 
#         cluster_rows = FALSE, 
#         cluster_columns = FALSE
#         )




# x <- d %>% 
#   arrange(patID_eye) %>% 
#   mutate(across(everything(), ~as.numeric(as.character(.))))
# 
# eye1 <- x %>%
#   filter(eye_num == 1)
# eye2 <- x %>%
#   filter(eye_num == 2)
# 
# eye1$eye_num <- NULL
# eye2$eye_num <- NULL
# 
# names(eye1)[3:ncol(eye1)] <- paste(names(eye1)[3:ncol(eye1)], "1", sep = "_")
# names(eye2)[3:ncol(eye2)] <- paste(names(eye2)[3:ncol(eye2)], "2", sep = "_")
# 
# wide <- left_join(eye1, eye2)
# 
# left <- "stratus_MT04_slope_2"
# right <- "stratus_MT04_slope_1"
# 
# include_vars <- sapply(names(x)[4:ncol(x)], function(column){
#   left <- paste0(column, "_2")
#   right <- paste0(column, "_1")
#   
#   var_left <- var(wide[[left]], na.rm = TRUE) %>% as.numeric()
#   var_right <- var(wide[right], na.rm = TRUE)%>% as.numeric()
#   
#   if (var_left == 0 | var_right == 0 | all(wide[[right]] == wide[[left]], na.rm = TRUE)) {
#         return(FALSE)
#     } else {
#         return(TRUE)
#     }
#   
# })
# 
# 
# include_pattern <- names(include_vars)[include_vars] %>% 
#   paste(., collapse = "|")
# 
# cor_m <- wide %>%
#   select(matches(include_pattern)) %>% 
#   cor(., use = "pairwise.complete.obs")
# 
# pdf("OneDrive - Lowy Medical Research Institute/Database_EDA/analysis_versions/version010/results/structures_ODvOS.pdf", height =100, width = 100)
# Heatmap(cor_m,
#         cluster_rows = FALSE, 
#         cluster_columns = FALSE)
# dev.off()




# values <- 1:100
# 
# n_values <- length(values)
# 
# probs <- seq(from = 10, to = 0, length.out = n_values)
# 
# x <- sample(values, size = 20, replace = TRUE, prob = probs)
# 
# plot(density(sqrt(x)))
# 
# # i = 1
# par(mfrow = c(2, 10));for(i in 1:10) {
#   vec <- rexp(1000, rate = i)
#   plot(density(vec), main = paste("rate = ", i), xlim= c(0, 2), ylim = c(0, 2))
#   qqnorm(vec, ylim = c(0, 3))
# 
# };par(mfrow = c(1,1))

# sum(c(567, 848, 1124, 526,  439, 2659, 1055, 1975, 3243) ) - 3243
# sapply(d, function(x){
#   
#   x <- x %>% as.character()
#   
#   grep("^777$|^7777$|^888$|^8888$|^999$|^9999$", x, value = TRUE)
#   
# })

# setdiff(parafovea_greater_than_one_month, MT5_greater_than_one_month)
# setdiff(parafovea_greater_than_one_month, MT5_rate_greater_than_one_month)
# setdiff(MT5_greater_than_one_month, MT5_rate_greater_than_one_month)
# setdiff(MT5_rate_greater_than_one_month,MT5_greater_than_one_month)
 


# d %>% 
#   filter(!is.na(MT05_rate)) %>% 
#   select(patID_eye, month, NEOVASC) %>% 
#   mutate(NEOVASC = as.numeric(scale(NEOVASC))) %>% 
#   group_by(patID_eye) %>% 
#   mutate(first_NEOVASC = dplyr::first(NEOVASC[!is.na(NEOVASC)]  ))
  
  
  
 
# d %>%
#   filter(!is.na(MT05_rate)) %>% 
#   select(patID_eye, month, all_of(powered_vars)) %>%
#   group_by(patID_eye, month) %>%
#   summarize(across(everything(), ~sum(!is.na(.)))) %>% 
#   ungroup() 
#   
# 
# num_at_each_month <- d %>%
#   filter(!is.na(MT05_rate)) %>% 
#   select(patID_eye, month, all_of(powered_vars)) %>%
#   group_by(month) %>%
#   summarize(across(everything(), ~sum(!is.na(.)))) %>% 
#   ungroup() %>% 
#   column_to_rownames("month")

# # modeling RET_ATROPH by MT5_quintile
# mod_d <- d %>%
#   select(SEX, patID, patID_eye , MT5_quintile , RET_ATROPH) %>%
#   mutate(MT5_quintile = MT5_quintile %>% as.character( ) %>% as.numeric())
#
# mod_d[, c("MT5_quintile", "RET_ATROPH")] <- mod_d[, c("MT5_quintile", "RET_ATROPH")] %>% scale()
#
#
# lmer(RET_ATROPH ~ SEX + (1|patID) + (1|patID_eye) + MT5_quintile, data = mod_d) %>% summary()
# lm( RET_ATROPH~ SEX + MT5_quintile, data = mod_d) %>% summary()
#
#
# count_d <- d %>%
#   select(MT5_quintile, RET_ATROPH) %>%
#   na.omit() %>%
#   dplyr::count(MT5_quintile, RET_ATROPH)
#
# x <- d %>%
#   select(MT5_quintile, RET_ATROPH) %>%
#   na.omit()
# mosaicplot(table(x))
#
# ggplot(count_d, aes(MT5_quintile, RET_ATROPH, fill = n))+
#   geom_tile()



# 
# # look at new dataset and compare missingness to current
# new_CF_FFA <- read.csv("/Users/christiananderson/Downloads/CF-FFA 4.csv")
# 
# names(new_CF_FFA) <- names(new_CF_FFA) %>% 
#   gsub("X\\.\\.", "", .) %>% 
#   gsub("\\.$", "", .) %>% 
#   gsub("\\.", "", .) %>% 
#   gsub("ICLO", "ICL0", .)
# 
# nas <- c( "null", "NULL", "Null", "na" , "Na","", " ", "\t", "\n", "  ",  "777", "888", "999", "7","8","9", "removed")
# new_CF_FFA <- new_CF_FFA %>% 
#   mutate(across(everything(), 
#                 ~gsub('"', "", 
#                       gsub("=", "", 
#                            gsub(" ", "", .x))
#                       )  
#                 
#      )
#   ) %>% 
#   mutate(across(everything(), 
#                 ~ifelse(. %in% nas, NA, .)
#   ))
#   
# 
# # make dataframe binary 
# df_binary <- new_CF_FFA %>% 
#   mutate(month = as.numeric(visitId)) %>% 
#   arrange(month) %>% 
#   mutate_all(~ifelse(is.na(.), 0, 1))
# 
# df_binary_colsums <- colSums(df_binary) %>% 
#   as.data.frame() %>% 
#   rownames_to_column("variable") %>% 
#   rename_with(~all_of(c(variable = "variable", . = "count_present"))) %>% 
#   mutate(index = 1:n()) %>% 
#   arrange(index) 
# 
# df_binary_colsums$variable <- factor(df_binary_colsums$variable, levels = df_binary_colsums$variable)
# 
# colsum_plot <- ggplot(df_binary_colsums, aes(x = variable ,y = count_present, fill = variable))+
#   geom_bar(stat = "identity", show.legend = FALSE, orientation = "x", width = .5, just = 0)+
#   theme_light()+
#   theme(axis.text.x = element_text( hjust = 1, vjust = NULL, size = 10, angle = 90,
#                                     margin = margin(t = -100, r = 0, b = 0, l = 0)
#                                     )
#   )
#   # geom_abline(intercept = nrow(df_binary), col = "blue", linetype = 2, slope = 0)+
#   # annotate("text", x = Inf, y = nrow(df_binary), label = paste("nrow(d):", nrow(df_binary)), hjust = 1, vjust = -0.5, color = "blue")
# 
# colsum_plot
# 
# 
# new_CF_FFA$createdAt %>% substr(., 1, 10) %>% unique() 





# new_CF_FFA$X..NEOVASC. %>% table()
# 
# d$NEOVASC %>% table()
# 
# raw$CF_FFA$NEOVASC %>% table()
# 
# dfs$CF_FFA$NEOVASC %>% table()


# d$patID %>% 
#   grep("\\D", ., value = TRUE, ignore.case =TRUE) %>% 
#   grep("UK", ., invert = TRUE, value = TRUE)
# 
# 
# 
# 
# vis <- as.character(seq(0, 200, 12)) %>%
#   sapply(function(x) {
#     if (nchar(x) < 2) {
#       paste0("00", x)
#     } else if (nchar(x) < 3) {
#       paste0("0", x)
#     } else {
#       x
#     }
#   })
# 
# d$visitID[!d$visitID %in% vis] %>% unique()
# 
# # which datasets have Baseline
# lapply(dfs, function(a_df) grep("Baseline", a_df[["visitId"]])  )
# #^ all
# 
# # # overlaps among mt datasets
# # b <- data.frame(strat = !is.na(joined_clean$stratus_MT05 ),
# #            cirr = !is.na(joined_clean$cirrus_MT05 ),
# #            spec = !is.na(joined_clean$spectralis_MT05 )
# #            )
# # b$sum <- b %>% 
# #   rowSums(.)
# # 
# # b %>% 
# #   filter(sum>2)

 
# pdf(paste0(wd_path, "results/VAA_analysis/imputation_testing/lassoNorm_with_predMatrix.pdf"), height = 15, width = 20)
# par(mfrow = c(6,7))
# 
# dev.off()
# par(mfrow = c(1,1))
# x <- c(d$stratus_MT05, d$cirrus_MT05, d$spectralis_MT05)
# 
# x <- d[, mt_cols] %>%  unlist()
# 
# x <- x[x<900]
# 
# boxplot(x, yaxt = 'n')
# y_range <- range(x, na.rm = TRUE)
# ticks <- seq(from = floor(y_range[1]), to = ceiling(y_range[2]), by = 10) 
# primary_line <- ticks[seq(1, length(ticks), 2)]
# secondary_line <- ticks[seq(2, length(ticks), 2)]
# axis(side = 2, at = ticks, las = 1)
# abline(h = primary_line, col = "grey", lty = 3)
# abline(h = secondary_line, col = "grey", lty = 4)
