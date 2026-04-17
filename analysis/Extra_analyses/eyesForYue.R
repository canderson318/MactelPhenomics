library(ComplexHeatmap)
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(tibble)

rm(list = ls())
gc()

#+ create a dataframe of  subjects (patient, eye, month) that do not have chew grade
#+  patID, 
#+  eye,
#+  visit 
#+  dates,
#+  scan dates,
#+  type of image, 
#+  type of macular thickness, 
#+  column of whether chew grade matched


# LOAD AND MUNG DATA ------
sce <- readRDS("~/Documents/LMRI/PRJ2024001/version015/processed_data/sce005.rds")

# join data together
colDat.Assay<- left_join(
  colData(sce) %>% as.data.frame(),
  assays(sce)$measure_vars %>% t() %>% as.data.frame() %>% rownames_to_column("subject")
)

# grab date fields
raw <- metadata(sce)$raw

# raw_filt <- lapply(raw, function(df){
#   
#   df <- df %>% 
#     select(patientId, eye, matches("visit"), matches("date"))
#   return(df)
# 
#   })

raw_filt <- raw[ sapply(raw,function(df){grepl("pat_eye_DATE", names(df)) %>%  any()}) ]
raw_filt$VAA <- NULL

raw_joined <- purrr::reduce(raw_filt, left_join,by = "pat_eye_DATE")

# add cirrus
raw_joined <- raw_joined %>% 
  left_join(raw$OCTMTCIRRUS)

raw_joined %>% glimpse()
raw_joined <- raw_joined[, grep("\\.x\\.x|\\.y", names(raw_joined), invert = TRUE)]
names(raw_joined) <- gsub("\\.x", "", names(raw_joined))
names(raw_joined)[duplicated(names(raw_joined))]

old <- c("patientId", "eye", "visitId")
new <- c("patID", "eye_num", "visitID")
for(i in 1:length(old)){
  names(raw_joined)[which(names(raw_joined) == old[i])] <- new[i]
}
names(raw_joined)[duplicated(names(raw_joined))]
names(raw_joined) <- make.unique(names(raw_joined), sep = "_")

raw_joined <- raw_joined %>% 
  unique() 

# add date fields to sub
subj_info <- 
  left_join(
    colDat.Assay %>%  select(subject, patID, patID_eye, eye_num, visitID, month, matches("_date"), chew_grade),
    raw_joined %>% select(patID, visitID, eye_num,matches("date"))
    # by = join_by(patID == patientId, visitID == visitId, eye_num == eye)
    )

# add tsne to sub
subj_info <- subj_info %>% 
  left_join(
    reducedDims(sce)$TSNE %>% as.data.frame() %>%  rownames_to_column("subject")
  ) %>% 
  unique()
subj_info %>% glimpse()

# remove unecessary cols
subj_info <- subj_info %>% 
  select(-matches("HRDEXUDATE"))

# add field of chew_grade_present
subj_info$chew_grade_present <- !is.na(subj_info$chew_grade)

#---#---#---#---
# ADD FIELDS OF HOW MANY VARS FILLED FOR EACH RECORD OF SUBJ_INFO ----
#---#---#---#---
# ## for each df from each grading df df of raw, for each subject, sum number of cols !na 
# raw_grading_names <- c("CF_FFA","FAFGRADE","OCTMAIN","OCTMTCIRRUS","OCTMTSPECTRALIS","OCTMTSTRATUS")
# subj_id_cols <- c("patientId", "eye", "visitId")
# # subj = subj_info$subject[1]
# # a_df_name <- raw_grading_names[1]
# 
# grading_df_subject_counts_observed <- list(); for(a_df_name in raw_grading_names){
#   
#   # grab df
#   a_df <- raw[[a_df_name]]
#   
#   # add count of vars observed for each subject to df of subject ids
#   subject_count_d <- a_df[, subj_id_cols]
#   subject_count_d[, paste0("count_vars_observed_", a_df_name)] <- rowSums(!is.na(a_df[, -which(names(a_df) %in% subj_id_cols )]))
#   subject_count_d[, paste0("ratio_vars_observed_", a_df_name)] <- rowSums(!is.na(a_df[, -which(names(a_df) %in% subj_id_cols )])) / ncol(a_df)
#     
#   
#   # save df of each subjects rowsums(!na) to list of dfs 
#   grading_df_subject_counts_observed[[a_df_name]] <- subject_count_d
# 
# }
# 
# # make df of subject grading data observation counts
# grading_df_subject_counts_observed_d <- purrr::reduce(grading_df_subject_counts_observed, left_join) %>% 
#   unique()
# 
# grading_df_subject_counts_observed_d %>% 
#   group_by(patientId, visitId, eye) %>%
#   filter(n()>1) %>% glimpse()
# # ^ different counts for  same visit, because of VAA wierdness (duplicate subjects)
# 
# # removing duplicated records
# grading_df_subject_counts_observed_d <- grading_df_subject_counts_observed_d %>% 
#   group_by(patientId, visitId, eye) %>%
#   # filter(n()>1) %>% glimpse()
#   # filter(!duplicated(across(everything())) | !duplicated(across(everything()), fromLast = TRUE)) %>%
#   mutate_all(~ifelse(n() > 1, dplyr::first(.) , .)) %>%
#   unique() %>% 
#   ungroup()
# 
# grading_df_subject_counts_observed_d %>% 
#   group_by(patientId, visitId, eye) %>%
#   filter(n()>1) %>% 
#   glimpse()
#   
# # join to subj_info
# subj_info_new <- subj_info %>% 
#   left_join(grading_df_subject_counts_observed_d, by = join_by(patID==patientId, eye_num==eye, visitID==visitId))
# 
# # check for duplicates
# subj_info_new %>%
#   filter(duplicated(across(everything())) | duplicated(across(everything()), fromLast = TRUE)) 
# 
# write.csv(subj_info_new, "/Users/christiananderson/OneDrive - Lowy Medical Research Institute/Database_EDA/Christian's files for Robbie/tsne_eye_feature_info.csv")


#---#---#---#---
# SAVE INFO WITHOUT MODALITY OBSERVATION COUNTS----
#---#---#---#---

names(subj_info)[grep("subject", names(subj_info))] <- "ID"
write.csv(subj_info, "/Users/christiananderson/Library/CloudStorage/Egnyte-lmri/Shared/LMRI/DS_Group/PRJ2024001/Christian's files for Robbie/tsne_eye_feature_info.csv", row.names = FALSE)

# colDat.Assay %>%
#   filter(is.na(eye)) %>%
#   mutate_all(~as.numeric(!is.na(.))) %>%
#   colSums()  %>%
#   barplot(las = 2, horiz = TRUE)

# why is eye missing?; dont know but raw measures and assay measures do not match for some of them>> because I fixed all records that "regress"
# # pats <- colDat.Assay$patID[is.na(colDat.Assay$eye)]
# pat <- sample( colDat.Assay$patID[is.na(colDat.Assay$eye)],1 )
# # pat <- "010036"
# 
# x <- raw_joined %>%
#   # filter(patientId %in% pats & visitId == "000") %>% 
#   filter(patID == pat & visitID == "000") %>% 
#   select(patID, visitID, eye_num, VTV05, LRT05)
# 
# y <- colDat.Assay %>% 
#   # filter(patID %in% x$patientId & visitID == "000") %>% 
#   filter(patID == pat & visitID == "000") %>% 
#   select(subject, eye_num, VTV05, LRT05) 
# 
# x;y
# d <- left_join(as.data.frame(colData(sce)), as.data.frame(t(assays(sce)$pmm_imp_of_scaled)) %>% rownames_to_column("subject"))
# d %>%
#   filter(patID == "010036" & visitID == "000") %>%
#   select(subject, VTV05, LRT05)
# 
# assays(sce)$measure_vars %>%
# # assays(sce)$scaled_measure_vars %>%
#   # assays(sce)$pmm_imp_of_scaled %>%
#   t() %>%
#   as.data.frame() %>%
#   rownames_to_column("subject") %>%
#   filter(subject %in% c("010036_1_0", "010036_2_0")) %>%
#   select(VTV05, LRT05)

