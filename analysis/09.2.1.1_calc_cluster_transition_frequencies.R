
options(max.print= 200)
pacman::p_load(zeallot ,tidyr,dplyr,reshape2,purrr,GGally,tidyverse,ggpubr,vcd,ggplot2,tidytext,gridExtra,ComplexHeatmap,zoo,mgcv,car,zoo,tools,mice,grDevices,ggrepel,grid,gplots,Metrics,stats,lme4,lmerTest,dendextend,magick,umap,scater,SingleCellExperiment,scran,bluster,slingshot,ggfortify,pheatmap,TSCAN,igraph,knitr, qs, patchwork, parallel,txtplot, HDF5Array)

rm(list = ls())
gc()


wd_path = "~/Documents/LMRI/PRJ2024001/version015/"
sce <- HDF5Array::loadHDF5SummarizedExperiment(sprintf("%sprocessed_data/sce002.h5", wd_path))


# \\\\
# \\\\
# Make transition freq table from sce
# \\\\
# \\\\

# Clusters are renamed following this scheme:
# New Cluster scheme:
# old_cluster new_cluster
# 3 1
# 6 2
# 4 3
# 9 4
# 1 5
# 2 6
# 12 7
# 10 8
# 11 9
# 8 10
# 5 11
# 7 12
# 13 13

d = colData(sce) %>% 
  data.frame() %>% 
  select(patID_eye, month, my_cluster)

# make new cluster names
key = tribble(~my_cluster, ~new_cluster,
         3,           1,
         6,           2,
         4,           3,
         9,           4,
         1,           5,
         2,           6,
        12,           7,
        10,           8,
        11,           9,
         8,          10,
         5,          11,
         7,          12,
        13,          13)

# d$new_cluster = d$my_cluster
d$new_cluster = key$new_cluster[match(d$my_cluster,key$my_cluster )]

head(d[,c("my_cluster","new_cluster")])

edges  = d %>% 
  # defactor
  mutate(my_cluster = as.numeric(as.character(my_cluster))) %>%
  group_by(patID_eye) %>% 
  arrange(patID_eye,month) %>% 
  # filter for where > 1 obs
  filter(n()>1) %>% 
  # filter for eyes that exist in > 1 cluster
  filter(n_distinct(new_cluster)>1) %>% 
  # find next cluster eye moves to 
  mutate(next_new_cluster = lead(new_cluster)) %>% 
  # remove rec where no next cluster (at eye's last visit)
  drop_na() %>% 
  ungroup() %>% 
  select(from = new_cluster, to = next_new_cluster) %>% 
  data.matrix()
  

trans = igraph::as_adjacency_matrix(igraph::graph_from_edgelist(edges)) %>% 
  data.matrix()

rownames(trans) = paste("from", 1:nrow(trans), sep ='')
colnames(trans) = paste("to", 1:ncol(trans), sep ='')
  
# \\\\
# \\\\
# Load transition freq table from supp material
# \\\\
# \\\\

# fromxto transition frequencies
t = "337,25,1,18,2,5,2,2,16,6,10,0,0
19,231,0,7,4,0,25,1,25,0,0,2,3
29,3,138,20,1,69,0,4,4,3,1,1,0
47,13,5,302,5,15,1,1,32,2,2,3,0
3,2,0,1,21,1,2,2,1,0,0,1,1
33,6,11,73,1,203,0,1,16,1,4,0,2
2,19,0,0,1,0,34,0,0,0,0,0,2
0,1,0,0,2,0,0,29,0,0,0,0,0
20,48,1,22,2,3,1,0,184,2,1,0,1
1,1,0,0,0,0,0,0,0,41,1,0,0
2,0,0,0,0,0,0,0,0,2,51,0,0
0,2,1,3,0,0,0,0,0,0,0,33,0
0,0,0,0,0,0,2,0,0,0,0,0,22"

M  = data.frame(data.table::fread(t, sep = ',',header = FALSE)) %>% data.matrix()
colnames(M) = paste("to",c("1","2","3","4","5","6","7","8","9","10","11","12","13"), sep = "")
rownames(M)= gsub("to","from", colnames(M))


# rename columns with correct 
for(i in 1:ncol(M)){
  coln = colnames(M)[i]
  rown = rownames(M)[i]
  old_num = gsub("to", "", coln) %>% as.numeric()
  
  if(old_num != as.numeric(gsub("from", "", rown))) stop("Numbers not matching")
  
  new_num = key$new_cluster[key$my_cluster == old_num]
  
  col_new_nm = paste0("to", new_num)
  row_new_nm = paste0("from", new_num)
  colnames(M)[i] = col_new_nm
  rownames(M)[i] = row_new_nm
}

ord <- function(v){
  ord = order(as.numeric(gsub("\\D", "", v)))
  return( ord)
}

M = M[ord(rownames(M)), ord(colnames(M))]


# \\\\
# \\\\
# Compare Transition Matrices
# \\\\
# \\\\
trans
M

sum(trans) - sum(M) # 0

sqrt(sum((trans-M)^2))


# RMSE
sqrt(sum((trans-M)^2)) # 0

#\\\\
#\\\\
# Calcualte transition ratios 
#\\\\
#\\\\


rowsumnorm <- function(X){  X %>% sweep(., 1, rowSums(M, na.rm = TRUE), '/') %>% round(3)}
colsumnorm <- function(X){  X %>% sweep(., 2, colSums(M, na.rm = TRUE), '/') %>% round(3)}


# ratio of leaving
## set diag 0 so not counting self transitions
diag(M) = NA

## remove row/colnames
dimnames(M) = list()

## Of eyes leaving this cluster, what ratio go to each other cluster?
row_sum_norm_M = rowsumnorm(M)
write.table(row_sum_norm_M, '/tmp/rowsumnorm.csv', row.names = FALSE, col.names=FALSE,sep = ',' )  

# ratio of entering
## Of eyes entering this cluster, what ratio come from each other cluster?
col_sum_norm_M = colsumnorm(M)
write.table(col_sum_norm_M, '/tmp/colsumnorm.csv', row.names = FALSE, col.names=FALSE,sep = ',' )  


#\\\\
#\\\\
# Plot transition tile
#\\\\
#\\\\

inds = expand.grid(seq(nrow(M)), seq(nrow(M)))

colnorm_d = data.frame(inds)
rownorm_d = data.frame(inds)

rownorm_d$ratio = sapply(1:nrow(inds), function(i){
  c(i,j) %<-% inds[i,]
  row_sum_norm_M[i,j]
})

colnorm_d$ratio = sapply(1:nrow(inds), function(i){
  c(i,j) %<-% inds[i,]
  col_sum_norm_M[i,j]
})

A = ggplot(rownorm_d, aes(factor(Var1),factor(Var2), fill = ratio))+
  geom_tile()+
  scale_fill_viridis_c(option = 'plasma', na.value = 'grey80')+
  xlab("Destination")+
  ylab("Source")+
  ggtitle("Of eyes leaving each cluster, what ratio go to each other cluster?")+
  theme(aspect.ratio = 1)

B = ggplot(colnorm_d, aes(factor(Var1),factor(Var2), fill = ratio))+
  geom_tile()+
  scale_fill_viridis_c(option = 'plasma', na.value = 'grey80')+
  xlab("Destination")+
  ylab("Source")+
  ggtitle("Of eyes entering each cluster, what ratio come from each other cluster?")+
  theme(aspect.ratio = 1)


ggsave(sprintf("%sresults/publication_figures/cluster_transition_from_to_tile.pdf",wd_path), A+B, height = 7, width = 14)



plot(diag(trans))
