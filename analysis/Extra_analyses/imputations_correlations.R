#+ 
#+ this script compiles csv files of imputation R values from all scripts up to and including 09.1_clustering_imp_and_dimred
#+ 
#+ Correlations were found by comparing running pearsons cor on actual versus predicted values for each variable imputed
#+ 

# SETUP ----
pacman::p_load(readr, stringr, magick,tidyr,dplyr,reshape2,purrr,GGally,tidyverse,ggpubr,vcd,ggplot2,tidytext,gridExtra,ComplexHeatmap,mgcv,car,zoo,tools, gt)

rm(list = ls())
gc()


wd_path <- "~/Egnyte - lmri/Shared/LMRI/DS_Group/PRJ2024001/analysis_versions/version015/"


# LOAD DATA ----
dir_path <- sprintf("%sprocessed_data/imputation_Rs", wd_path)
files <- list.files(dir_path)
files_paths <- paste(dir_path, files, sep = "/")

corr_ds <- map(files_paths , ~{ read_csv(.)})

names(corr_ds) <- str_remove_all(files, "\\.csv")


# COMBINE DFS ----

# what vars were imputed more than once??
map(corr_ds, ~{.[['variable']]}) %>% 
  unlist() %>% 
  as.vector() %>% 
  table() %>% 
  sort()

# add imputation var to dfs
corr_ds <- map(names(corr_ds), ~{
  nm <- .
  d <-  corr_ds[[nm]]
  nm <- str_remove_all(nm, "_correlation.*|_R.*")
  d$imputation = nm
  return(d)
})

# add ratio_missing to ds in list
corr_ds <- map(corr_ds, ~{
  x <- .
  if(!"ratio_missing"%in% names(x)){
    x$ratio_missing <- NA
  }
  return(x)
})
# combine all together
corr <- do.call(rbind, corr_ds)

# add var for whether between datatyps or within
corr$`Within Type Imputation` <- ifelse(corr$imputation != "all_vars_impute", T, F)

corr$q_value <- p.adjust(corr$p_value, method = "BH")

corr <- corr %>% 
  mutate(across(matches("_value"), ~case_when(
    # . > 0.05 ~ "not significant", 
    # . < 0.01 ~ "significant",
    # . < 0.05 ~ "moderately significant",
    # . < 0.1 ~ "marginally significant",
    . > 0.05 ~ "not significant",
    . < 0.001 ~ "***",
    . < 0.01 ~ "**",
    . < 0.05 ~ "*",
    . < 0.1 ~ ".",
    TRUE ~ as.character(.)
  ), .names = "{.col}_signif")) %>% 
  mutate(across(matches("signif"), ~factor(., levels = c("***","**", "*", ".", "not significant"))))



# pal <- wesanderson::wes_palette(7)
# ggplot(corr[!corr$`Within Type Imputation`,], aes(R,reorder(variable, R), fill = p_value_signif))+
#   geom_bar(stat = "identity", show.legend = TRUE, position = "dodge", alpha = .5)+
#   geom_text(aes(label = sprintf("%s%%", round(ratio_missing, 3)*100) ))+
#   ggtitle("Imputation performance across variables for between data type imputation", "Percentages indicate percentage missing in full dataset")+
#   scale_fill_viridis_d(begin = 1, end = 0)+
#   theme_bw()
# 
# ggsave(sprintf("%sresults/clustering/09_clustering/imputation_QC/between_imputations_R_bar.png", wd_path), height = 12, width = 10)
### !!! ^ plotting in 9.1_clustering_imp_and_dimred

ggplot(corr[!corr$`Within Type Imputation`,], aes(R,reorder(variable, R), fill = q_value_signif))+
  geom_bar(stat = "identity", show.legend = TRUE, position = "dodge", alpha = .5)+
  geom_text(aes(label = sprintf("%s%%", round(ratio_missing, 3)*100) ))+
  ggtitle("Imputation performance across variables for between data type imputation", "Percentages indicate percentage missing in full dataset")+
  scale_fill_viridis_d(begin = 1, end = 0)+
  theme_bw()

ggplot(corr, aes(R,reorder(variable, R), fill = imputation))+
  geom_bar(stat = "identity", show.legend = TRUE, position = "dodge", alpha = .7)+
  ggtitle("Imputation performance across variables")+
  scale_fill_manual(values = c("all_vars_impute" = pal[1],"MT_rescaling" =  pal[3], "VAA_imputation" = pal[5]))+
  theme_bw()
ggsave(sprintf("%sresults/clustering/09_clustering/imputation_QC/all_imputations_R_bar.png", wd_path), height = 12, width = 10)


ggplot(corr, aes(R, fill = imputation))+
  # stat_density(alpha = .4, bw = .01)+
  geom_density(alpha = .4, bw = .04)+
  ggtitle("Imputation R distributions")+
  theme_bw()
ggsave(sprintf("%sresults/clustering/09_clustering/imputation_QC/all_imputations_R_densities.png", wd_path), height = 8, width = 12)

ggplot(corr[!corr$`Within Type Imputation`,], aes(R))+
  # stat_density(alpha = .4, bw = .01)+
  geom_density(alpha = .4, bw = .04, fill = pal[1])+
  ggtitle("Between data type imputation R distributions")+
  theme_bw()
ggsave(sprintf("%sresults/clustering/09_clustering/imputation_QC/between_imputations_R_density.png", wd_path), height = 8, width = 12)

ggplot(corr, aes(R, fill = imputation))+
  # stat_density(alpha = .4, bw = .01)+
  geom_density(alpha = .4, bw = .04)+
  geom_text_repel(aes(y = -.5, label = variable, jitter(R, 1)), angle = 90, size= 2, hjust = 1,min.segment.length = Inf, max.overlaps = 100, box.padding = .1)+
  ggtitle("Imputation R distributions")+
  theme_bw()
ggsave(sprintf("%sresults/clustering/09_clustering/imputation_QC/all_imputations_R_densities_labeled.png", wd_path), height = 8, width = 12)

# x <- corr  %>% drop_na()
# corr_m <- as.matrix( x%>% select(-c(variable, imputation)))
# rownames(corr_m) <- with(x,paste(variable,imputation, sep = "__"))
# 
# dist <- dist(corr_m)
# 
# # str_dist <- stringdist::stringdistmatrix(matrix(rownames(corr_m)), method = "qgram")
# # combo_dist <- (str_dist/2)+dist
# 
# hc <- hclust(dist, method = "ward.D2")
# labels(hc) <- rownames(corr_m)[hc$order]
# plot(hc)


names(corr) <- str_to_sentence(names(corr))

corr$R <- round(corr$R, 3)
corr <- corr[order(corr$Imputation, -corr$R),]

# make pretty table
gt_table <- corr %>%
  gt() %>% 
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_column_labels()
  ) %>% 
  tab_options(
    table.font.size = px(10),      
    data_row.padding = px(2),      
    column_labels.font.size = px(12),
    column_labels.padding = px(2)
  )

gtsave(gt_table, sprintf("%sresults/clustering/09_clustering/imputation_QC/all_imputations_R_table.pdf", wd_path))

write.xlsx(corr, sprintf("%sprocessed_data/all_imputations_corr.xlsx", wd_path))
