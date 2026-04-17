
pacman::p_load(tidyr, dplyr, reshape2, purrr, GGally, tidyverse, ggplot2, SingleCellExperiment, future, HDF5Array, patchwork, magrittr)


rm(list = ls())
gc()

wd_path = "~/Documents/LMRI/PRJ2024001/version015/"

# load experiment
sce <- loadHDF5SummarizedExperiment(paste0(wd_path, "processed_data/sce003.h5"))

d <- assays(sce)$pmm_imp_of_scaled %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  left_join( colData(sce) %>% data.frame())


v = "chew_grade_derived_mine_loose"
sapply(grep("chew_grade*", names(d), value = TRUE), function(v){
  x = d[, c("patID_eye", "month", "my_cluster", v)]
  contingency = table(x[['my_cluster']], x[[v]])  
  chi = chisq.test(contingency)
  c("X" = chi$statistic, "p" = chi$p.value)
}, simplify = TRUE)

table(d$my_cluster, d$chew_grade_derived_mine_loose)  %>% chisq.test()

# G-test
observed = table(d$my_cluster, d$chew_grade_derived_mine_loose)  
rowsums = rowSums(observed)
colsums = colSums(observed)
expected = outer(rowsums, colsums, "*") / sum(observed)
logratio = log(observed/expected)
crit_val = 2 * sum(ifelse(observed>0, observed * logratio, 0))
df = (length(rowsums)-1)*(length(colsums)-1)
p = pchisq(crit_val, df = df,lower.tail = F)

c(
  G = crit_val,
  p = p, 
  df = df
  )

DescTools::GTest(observed, expected)
