# Mon May  5 11:49:35 2025 ------------------------------

rm(list = ls()); gc()

# Packages ----
pacman::p_load(tidyr,tibble, dplyr, reshape2, purrr, ggplot2, gridExtra, ComplexHeatmap, viridis,  grDevices, ggrepel, grid, Metrics, stats,  scater, SingleCellExperiment,  ggfortify, HDF5Array, RColorBrewer, magrittr, patchwork, stringr)

source("/Users/canderson/Documents/LMRI/PRJ2024001/version015/analysis/00_MyFunctions.R")


# Load Data ----
wd_path = "~/Documents/LMRI/PRJ2024001/version015/"

sce <- HDF5Array::loadHDF5SummarizedExperiment(sprintf("%sprocessed_data/sce005.h5", wd_path))

# REDO CLUSTER ORDERING AND PC1 ----
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
  select(my_cluster, pseudoprog) %>% 
  summarize(mean_pp = mean(pseudoprog),median_pp = median(pseudoprog), .by = my_cluster) 
# my_cluster new_cluster
#          3           1
#          6           2
#          4           3
#          9           4
#          1           5
#          2           6
#         12           7
#         10           8
#          11          9
#           8         10
#           5         11
#           7         12
#          13         13
 
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


# plot path
sprintf("%sresults/poster_plots/",wd_path) %>% 
  {if(!dir.exists(.)) dir.create(.)}


p <- plotTSNE(sce[ , !is.na(sce$chew_grade_derived_mine_loose)], colour_by="chew_grade_derived_mine_loose", text_by = "my_cluster", text_colour = "red") +
  ggtitle("TSNE colored by Chew grading score")+
  theme(legend.background = element_rect(fill = "white"),
        legend.margin = ggplot2::margin(t = 4, r = 4, l = 4, b = 4),
        aspect.ratio = 1, 
        title = element_text(color = "white"))

# Extract the data used in the plot
label_data <- ggplot_build(p)$data[[2]]  # Second layer is usually geom_text

# Rebuild plot without the geom_text layer
p$layers[[2]] <- NULL

# Add geom_label manually
p + geom_label(
  data = label_data,
  aes(x = x, y = y, label = label, color = label),
  color = "black",
  fill = "white",
  label.size = 0.2
)

h=6
w=9
ggsave(sprintf("%sresults/poster_plots/chew_TSNE.png",wd_path), height = h, width = w, bg = "black", dpi = 200)


cor_result <- cor.test(as.numeric(sce$chew_grade_derived_mine_loose), sce$pseudoprog)

ggplot(data.frame(colData(sce)), aes(factor(chew_grade_derived_mine_loose), pseudoprog)) +
  geom_boxplot(aes(fill = factor(chew_grade_derived_mine_loose)), show.legend = FALSE) +
  # stat_smooth(aes(group = 1), method = "lm", color = "black", linetype = 1, linewidth = 0.5) +
  scale_fill_viridis_d(begin = .2)+
  theme_bw() +
  xlab("Chew Grade") +
  ylab("Pseudoprogression") +
  labs(fill = "MacTel Stage")+
  ggtitle("MacTel stage (Chew grade) correlation with Pseudoprogression")+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA))+
  theme(aspect.ratio = 1,
        axis.title.x = element_text(),
        axis.title.y = element_text(angle = 90))+
  labs(subtitle =  sprintf("R = %.2f, %s", cor_result$estimate, ifelse(cor_result$p.value<.0001, "p-value < 0.0001",  sprintf("p-value = %.2e", cor_result$p.value))))+
  # ggpubr::stat_compare_means(method = "t.test", comparisons = matrix(c(1:6, 2:7),byrow = T, ncol = 6) %>%split(col(.)), color = rgb(0,0,0,.6))+
  geom_blank()

ggsave(sprintf("%sresults/poster_plots/chew_pp.pdf",wd_path), height = 10, width = 11) 
ggsave(sprintf("%sresults/publication_figures/chew_pp.pdf",wd_path), height = 10, width = 11) 

ggplot(data.frame(colData(sce)), aes(factor(my_cluster), pseudoprog)) +
  geom_boxplot(aes(fill = factor(my_cluster)), show.legend = FALSE) +
  # stat_smooth(aes(group = 1), method = "lm", color = "black", linetype = 1, linewidth = 0.5) +
  scale_fill_manual(values = cluster_cols)+
  theme_bw() +
  xlab("Cluster") +
  ylab("Pseudoprogression") +
  labs(fill = "MacTel Stage")+
  ggtitle("Cluster Pseudoprogression")+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA))+
  theme(aspect.ratio = 1,
        axis.title.x = element_text(),
        axis.title.y = element_text(angle = 90))+
  geom_blank()

ggsave(sprintf("%sresults/poster_plots/cluster_pp.pdf",wd_path), height = 10, width = 11) 
ggsave(sprintf("%sresults/publication_figures/cluster_pp.pdf",wd_path), height = 10, width = 11) 


rdims <- reducedDim(sce, "PCA") %>% 
  data.frame() %>% rownames_to_column("sample") %>% 
  mutate(PC1 = PC1) %>% 
  left_join(reducedDim(sce, 'TSNE') %>% data.frame() %>% rownames_to_column("sample") ) %>% 
  left_join(colData(sce) %>% data.frame() %>% select(-matches("PC",ignore.case = F)), by = "sample") 
  
pc_tree <- metadata(sce)$trajectory_data$pca_line_data %>% 
  mutate(PC1 = PC1) %>% 
  select(edge, PC1,PC2) 

tsne_tree <- metadata(sce)$trajectory_data$tsne_line_data %>% 
  select(edge, TSNE1,TSNE2) 

pc_centers <- pc_tree %>%
  group_by(edge) %>%
  mutate(which = row_number()) %>%
  rowwise() %>%
  mutate(clust = str_extract_all(edge, "\\d+")[[1]][which]) %>%
  ungroup() %>%
  data.frame() %>% 
  mutate(route = case_when(
    clust %in% c(5,8,9,11,10) ~ "a", 
    clust %in% c(12,6,13,4) ~ "b",
    clust %in% c(2,1,3,7) ~ "x", 
    TRUE ~ NA
  ))
# mutate(route = case_when(
#     clust %in% c(1,10,11,5,8) ~ "a", 
#     clust %in% c(7,2,13,9) ~ "b",
#     clust %in% c(6,3,4,12) ~ "x", 
#     TRUE ~ NA
#   ))

pc_tree <- left_join(pc_tree,pc_centers )

tsne_centers <- tsne_tree %>%
  group_by(edge) %>%
  mutate(which = row_number()) %>%
  rowwise() %>%
  mutate(clust = str_extract_all(edge, "\\d+")[[1]][which]) %>%
  ungroup() %>%
  data.frame() %>% 
  mutate(route = case_when(
    clust %in% c(1,10,11,5,8) ~ "a", 
    clust %in% c(7,2,13,9) ~ "b",
    clust %in% c(6,3,4,12) ~ "x", 
    TRUE ~ NA
  ))

tsne_tree <- left_join(tsne_tree,tsne_centers )


  ggplot()+
  # geom_point(data = rdims %>% filter(!is.na(pc_lda_route)), aes(PC1, PC2),color = 'darkgrey', size = 3)+
  geom_point(data = rdims %>% filter(!is.na(pc_lda_route)), aes(PC1, PC2, color = my_cluster), size = 3)+
  geom_line(data = pc_tree, aes(PC1,PC2, group = edge))+
  geom_label(data = pc_centers, aes(PC1,PC2, label = clust))+
  scale_color_manual(values = cluster_cols)+
  xlab(sprintf("PC1 %.0f%%", attributes(reducedDim(sce, 'PCA'))$percentVar[1] ))+
  ylab(sprintf("PC2 %.0f%%", attributes(reducedDim(sce, 'PCA'))$percentVar[2] ))+
  labs(color = "Cluster")+
  ggtitle("Principal components with minimum spanning pc_tree")+
  theme_bw()+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA)
  )+
  # theme_void()+
  theme(aspect.ratio = 1,
          axis.title.x = element_text(),
          axis.title.y = element_text(angle = 90)
        )
  
  ggsave(sprintf("%sresults/poster_plots/cluster_cols_pc_tree_pc.pdf",wd_path), height = 10, width = 11) 
  ggsave(sprintf("%sresults/publication_figures/cluster_cols_pc_tree_pc.pdf",wd_path), height = 10, width = 11) 
  
  
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  route_pal <- c(gg_color_hue(3)[-3], "darkgrey")
  
  
  ggplot()+
  geom_point(data = rdims %>% filter(!is.na(pc_lda_route)) %>% mutate(pc_lda_route = ifelse(as.character(pc_lda_route)=="x", "none", as.character(pc_lda_route))) , aes(PC1, PC2, color = pc_lda_route), size = 3)+
  geom_line(data = pc_tree, aes(PC1,PC2, group = edge))+
  geom_label(data = pc_centers, aes(PC1,PC2, label = clust))+
  xlab(sprintf("PC1 %.0f%%", attributes(reducedDim(sce, 'PCA'))$percentVar[1] ))+
  ylab(sprintf("PC2 %.0f%%", attributes(reducedDim(sce, 'PCA'))$percentVar[2] ))+
  labs(color = "Trajectory")+
  ggtitle("Eye trajectories")+
  theme_bw()+
  scale_color_manual(values = setNames(route_pal, c("a", "b", "none"))) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA)
  )+
  # theme_void()+
  theme(aspect.ratio = 1,
          axis.title.x = element_text(),
          axis.title.y = element_text(angle = 90)
        )
  
ggsave(sprintf("%sresults/poster_plots/pc_lda_route_pc.pdf",wd_path), height = 10, width = 11)
ggsave(sprintf("%sresults/publication_figures/pc_lda_route_pc.pdf",wd_path), height = 10, width = 11)
  
ggplot()+
  geom_point(data = rdims %>% filter(!is.na(pc_lda_route)) %>% mutate(pc_lda_route = ifelse(as.character(pc_lda_route)=="x", "none", as.character(pc_lda_route))) , 
             aes(TSNE1, TSNE2, color = pc_lda_route), size = 3)+
  geom_line(data = tsne_tree, aes(TSNE1,TSNE2, group = edge))+
  geom_label(data = tsne_centers, aes(TSNE1,TSNE2, label = clust))+
  labs(color = "Trajectory")+
  ggtitle("Eye trajectories")+
  theme_bw()+
  scale_color_manual(values = setNames(route_pal, c("a", "b", "none"))) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA),
        legend.box.background = element_rect(fill = "white", color = NA)
        # plot.background = element_rect(fill = "transparent", color = NA),
        # legend.background = element_rect(fill = "transparent", color = NA),
        # legend.box.background = element_rect(fill = "transparent", color = NA)
  )+
  # theme_void()+
  theme(aspect.ratio = 1,
        axis.title.x = element_text(),
        axis.title.y = element_text(angle = 90)
  )

ggsave(sprintf("%sresults/poster_plots/pc_lda_route_tsne.pdf",wd_path), height = 10, width = 11)


ggplot()+
  geom_point(data = rdims %>% filter(!is.na(pc_lda_route)) %>% mutate(pc_lda_route = ifelse(as.character(pc_lda_route)=="x", "none", as.character(pc_lda_route))) , 
             aes(TSNE1, TSNE2, color = pc_lda_route), size = 3)+
  geom_line(data = tsne_tree, aes(TSNE1,TSNE2, group = edge),
            color = "white", linewidth = 2)+
  geom_label(data = tsne_centers, aes(TSNE1,TSNE2, label = clust), size = 5)+
  labs(color = "Trajectory")+
  ggtitle("Eye trajectories")+
  theme_bw()+
  scale_color_manual(values = setNames(route_pal, c("a", "b", "none"))) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        # panel.background = element_rect(fill = "white", color = NA),
        # plot.background = element_rect(fill = "white", color = NA),
        # legend.background = element_rect(fill = "white", color = NA),
        # legend.box.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA)
  )+
  # theme_void()+
  theme(aspect.ratio = 1,
        axis.title.x = element_text(),
        axis.title.y = element_text(angle = 90)
  )
  
ggsave(sprintf("%sresults/poster_plots/pc_lda_route_tsne_no_bg.pdf",wd_path), height = 10, width = 11)
  
  ggplot() +
    geom_point(data = rdims, aes(PC1, PC2), alpha = 0) +
    geom_line(data = pc_tree, aes(PC1, PC2, group = edge, color = route)) +
    geom_point(data = pc_centers , aes(PC1, PC2, color = route),alpha = 1, show.legend = TRUE) +
    geom_label(data = pc_centers , aes(PC1, PC2, label = clust, color = route), show.legend = FALSE) +
    scale_color_manual(values = setNames(route_pal , c("a", "b", "x"))) +
    labs(color = "Trajectory")+
    xlab(sprintf("PC1 %.0f%%", attributes(reducedDim(sce, 'PCA'))$percentVar[1])) +
    ylab(sprintf("PC2 %.0f%%", attributes(reducedDim(sce, 'PCA'))$percentVar[2])) +
    theme_void() +
    theme(aspect.ratio = 1, 
          legend.position.inside =  c(-30,-20))
  ggsave(sprintf("%sresults/poster_plots/pc_lda_route_pc_tree.pdf",wd_path), height = 10, width = 11)
  

  
  # make df of patID, right pc_lda_route and left pc_lda_route; make A and C same, ditto for B,D
  pc_lda_route_eyes <- rdims %>%
    dplyr::select( patID_eye,patID, eye_num, pc_lda_route) %>%
    distinct() %>% 
    mutate(eye_num = ifelse(eye_num==1, "right", ifelse(eye_num==2, "left",eye_num ))) %>% 
    filter(!is.na(eye_num)) %>% 
    mutate(pc_lda_route = as.character(pc_lda_route)) %>% 
    pivot_wider(id_cols = patID, names_from = eye_num, values_from = pc_lda_route) %>% 
    filter(if_all(-patID, ~.!="x"))
  
  # do right route = left routes
  pc_lda_route_table <- table(pc_lda_route_eyes[, c("left", "right")]) 
  pc_lda_chsq <- chisq.test(pc_lda_route_table)
  
  # dimnames(pc_lda_route_table) <- list(paste("Left ", c("a", "b", "x"), sep = ""), paste("Right ", c("a", "b", "x"), sep = ""))
  dimnames(pc_lda_route_table) <- list(paste("Left ", c("a", "b"), sep = ""), paste("Right ", c("a", "b"), sep = ""))
  
  
  pc_lda_route_table_d <- pc_lda_route_table %>% 
    reshape2::melt(varnames = c("Left", "Right")) %>% 
    mutate_all(~gsub("Right|Left| ","",.)) %>% 
    mutate(value = as.numeric(value)) %>% 
    mutate(across(c(Right), ~factor(., levels = c("a", "b")))) %>% 
    mutate(across(c(Left), ~factor(., levels = c("b", "a"))))
  
  pc_lda_route_table_tile <-
    ggplot(pc_lda_route_table_d,aes(Right, Left, fill = value, label = value))+
    geom_tile(show.legend = FALSE)+
    geom_text()+
    theme_minimal()+
    theme(aspect.ratio = 1)+
    scale_fill_continuous(high = "red", low = "white")+
    ggtitle("Patient trajectory eye symetry",sprintf("X^2 = %.0f, p-value = %.02e", pc_lda_chsq$statistic, pc_lda_chsq$p.value) )+
    geom_blank()
  
  ggsave(paste0(wd_path, "results/poster_plots/pc_lda_route_assignment_contingency_table.pdf"), pc_lda_route_table_tile, height = 5, width = 5)
    
  
p1 <-   ggplot()+
  geom_point(data = rdims , aes(PC1, PC2, color = my_cluster), size = 1)+
  geom_line(data = pc_tree, aes(PC1,PC2, group = edge))+
  geom_label(data = pc_centers, aes(PC1,PC2, label = clust))+
  xlab(sprintf("PC1 %.0f%%", attributes(reducedDim(sce, 'PCA'))$percentVar[1] ))+
  ylab(sprintf("PC2 %.0f%%", attributes(reducedDim(sce, 'PCA'))$percentVar[2] ))+
  labs(color = "Cluster")+
  ggtitle("Clusters on PCs with Minimum Spanning Tree")+
  theme_bw()+
  scale_color_manual(values = cluster_cols)+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA))+
  theme(aspect.ratio = 1,
        axis.title.x = element_text(),
        axis.title.y = element_text(angle = 90))
  
p2 <-  ggplot()+
  geom_point(data = rdims , aes(PC1, PC2, color = pseudoprog), size = 1)+
  geom_line(data = pc_tree, aes(PC1,PC2, group = edge))+
  geom_label(data = pc_centers, aes(PC1,PC2, label = clust))+
  xlab(sprintf("PC1 %.0f%%", attributes(reducedDim(sce, 'PCA'))$percentVar[1] ))+
  ylab(sprintf("PC2 %.0f%%", attributes(reducedDim(sce, 'PCA'))$percentVar[2] ))+
  labs(color = "Pseudoprogression")+
  ggtitle("Pseudoprogression on PCs with Minimum Spanning Tree")+
  theme_bw()+
  scale_color_viridis_c()+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA))+
  theme(aspect.ratio = 1,
        axis.title.x = element_text(),
        axis.title.y = element_text(angle = 90))
  
p3 <-   ggplot()+
  geom_point(data = rdims , aes(TSNE1, TSNE2, color = my_cluster), size = 1)+
  geom_line(data = tsne_tree, aes(TSNE1,TSNE2, group = edge))+
  geom_label(data = tsne_centers, aes(TSNE1,TSNE2, label = clust))+
  labs(color = "Cluster")+
  ggtitle("Clusters on TSNEs with Minimum Spanning Tree")+
  theme_bw()+
  scale_color_manual(values = cluster_cols)+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA))+
  theme(aspect.ratio = 1,
        axis.title.x = element_text(),
        axis.title.y = element_text(angle = 90))
  
p4 <-   ggplot()+
  geom_point(data = rdims , aes(TSNE1, TSNE2, color = pseudoprog), size = 1)+
  geom_line(data = tsne_tree, aes(TSNE1,TSNE2, group = edge))+
  geom_label(data = tsne_centers, aes(TSNE1,TSNE2, label = clust))+
  labs(color = "Pseudoprogression")+
  ggtitle("Pseudoprogression on TSNEs with Minimum Spanning Tree")+
  theme_bw()+
  scale_color_viridis_c()+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA))+
  theme(aspect.ratio = 1,
        axis.title.x = element_text(),
        axis.title.y = element_text(angle = 90))

grid <- (p1+p3)/(p2+p4)+plot_layout() & theme(
  plot.background = element_rect(fill = "transparent", color = NA)
)

ggsave(paste0(wd_path, "results/poster_plots/pc_tnse_clust_pp_grid.pdf"),grid, height = 10, width = 12, bg = "transparent")
ggsave(paste0(wd_path, "results/publication_figures/pc_tnse_clust_pp_grid.pdf"),grid, height = 10, width = 12, bg = "transparent")


