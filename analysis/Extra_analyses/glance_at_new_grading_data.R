
library(ComplexHeatmap)
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(tibble)
library(openxlsx)


rm(list = ls())
gc()




new_cf <- read.csv("/Users/christiananderson/Library/CloudStorage/Egnyte-lmri/Shared/LMRI/DS_Group/PRJ2024001/raw_data/06/CF-FFA.csv")
 

old_cf <- read.csv("/Users/christiananderson/Library/CloudStorage/Egnyte-lmri/Shared/LMRI/DS_Group/PRJ2024001/raw_data/05_imageGradingAndVisualAcuity/consolidated_sheet_15_4_2024/CF-FFA.csv")

clean.it <- function(a_df, rm_bad_image_codes = TRUE){
  
  values_to_NA <- c( "null", "NULL", "Null", "na" , "Na","", " ", "\t", "\n", "  ")
  if(rm_bad_image_codes ) values_to_NA <-  c(values_to_NA,  "777", "7777", "888", "8888", "999", "9999", "7","8","9", "removed")
  
  names(a_df) <- gsub("X..|\\.$", "", names(a_df))
  
  a_df <- a_df %>% 
    mutate_all(~case_when(
      . %in% values_to_NA ~ NA, 
      TRUE ~ .
    )) %>% 
    mutate_all(~gsub('=\"|\"$', "", .))
  
  return(a_df)
}

new_cf <- clean.it(new_cf, rm_bad_image_codes = FALSE)
old_cf <- clean.it(old_cf,rm_bad_image_codes = FALSE)

not.missing.binary <-  function(a_df){
  
  a_df <- a_df %>% 
    mutate_all(~as.numeric(!is.na(.)))
  
  return(a_df)
}

new_cf_obs_binary <- not.missing.binary(new_cf)
old_cf_obs_binary <- not.missing.binary(old_cf)

new_colsums <- colSums(new_cf_obs_binary) %>% 
  as.data.frame() 
old_colsums <- colSums(new_cf_obs_binary) %>% 
  as.data.frame() 

new_colsums$variable <- rownames(new_colsums)
old_colsums$variable <- rownames(old_colsums)

names(new_colsums)[1] <- "new_count_observed"
names(old_colsums)[1] <- "old_count_observed"

setdiff(old_colsums$variable, new_colsums$variable)
setdiff(new_colsums$variable,old_colsums$variable)

colsums <- left_join(new_colsums, old_colsums) %>% 
  pivot_longer(c(old_count_observed, new_count_observed), values_to = "count", names_to = "which_df") %>% 
  mutate(which_df = sapply(str_split(which_df, "_"), `[`, 1))

ggplot(colsums, aes(y = count, x = variable, fill = which_df))+
  geom_bar(stat = "identity", position = "dodge")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

diff <- (colSums(new_cf_obs_binary) -colSums(old_cf_obs_binary)) %>% 
  as.data.frame() 
  
names(diff) <- "count"
diff$variable <- rownames(diff)

ggplot(diff, aes(x =variable, y = count))+ 
  geom_bar(stat = "identity")+ 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,  hjust = 1))+ 
  ggtitle("Difference in non-missing observations between old and new CF table")

