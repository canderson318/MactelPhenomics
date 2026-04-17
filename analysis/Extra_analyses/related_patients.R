# I've uploaded here (Files > Cleaned Data) the cleaned data from Metabolon 2.0.
# Look at the patients file and look for the MACID column.
# This contains the IDs of these patients that should match the ones in the database.
# Divide these patients into families using the data from the database using Proband ID.

library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(purrr)
library(tidyverse)
library(lubridate)
library(knitr)
library(openxlsx)
library(ComplexHeatmap)
library(VennDiagram)

rm(list = ls())
gc()


setwd("~/Documents/LMRI/PRJ2024001/version015/")

# DATA STUFF ---------------------------------------------------------------
# make dir for results
if(!dir.exists("~/Documents/LMRI/PRJ2024001/version015/results/related_metabolon2.0_patients/")) {
  dir.create("~/Documents/LMRI/PRJ2024001/version015/results/related_metabolon2.0_patients/")
}
##### read in dfs 
dfs <- readRDS("~/Documents/LMRI/PRJ2024001/version015/processed_data/dfs_list_01.rds")
famstatus <- dfs$FAMSTATUS
enr <- dfs$ENR

##### read in sample info 
# load in sample info 
sample_d <- read.csv("/Users/christiananderson/Library/CloudStorage/Egnyte-lmri/Shared/LMRI/DS_Group/PRJ2024001/raw_data/metabolon2.0_samples/Samples.csv")

#### clean data
##### join the family info sheets together
fam_d <- full_join(
  famstatus %>% select(patientId, FAMILY_ID, FAMILY_STATUS, matches("pat")) %>% unique,
  enr %>% select(patientId, FAMILYID, PROBANDPW, matches("pat"))%>% unique
)

##### rename family_ID to what it actually is
names(fam_d)[grep("FAMILY_ID", names(fam_d))] <- "FAMILYID_PROBANDPW"

##### fill familyID where missing with the chars before "-" 
###### -> '[' , 1 accesses first element of each element of list outputted by strsplit
fam_d$FAMILYID[is.na(fam_d$FAMILYID)] <- sapply(strsplit(fam_d$FAMILYID_PROBANDPW[is.na(fam_d$FAMILYID)], split = "-"), `[`, 1)

# do same for proband pw from family_id_probandpw
fam_d$PROBANDPW[is.na(fam_d$PROBANDPW)] <- sapply(strsplit(fam_d$FAMILYID_PROBANDPW[is.na(fam_d$PROBANDPW)], split = "-"), `[`, 2)

#### remove leading zeros of patientID
fam_d$patientId_new <- gsub("^0|^00|^000|^0000|^00000|^000000", '', fam_d$patientId) %>% 
  ifelse(nchar(.) == 0, "0", .)



# PATIENTS SHARED BY BOTH DFS  -----------------------------------------
patients_list <- list(sample_d_patients = sample_d$MACID %>% unique, 
                      d_patients_filtered = fam_d$patientId_new[!is.na(fam_d$FAMILYID)] %>% unique,
                      d_patients = fam_d$patientId_new %>% unique
                      )

venn.diagram(
  main = "patients from my df and Samples where FAMILYID is not missing", 
  x = patients_list[c(1,2)],
  category.names = c("sample_d_patients", "d_patients"),
  filename = "results/related_metabolon2.0_patients/patient_venn_FAMILYID_observed.png",
  disable.logging = TRUE,
  output=TRUE, 
  lwd = 2,
  lty = 'blank',
  fill = c("#E69F00", "#56B4E9"),
  cex = .9,
  fontface = "italic",
  cat.dist = c(-0.05, -0.05)
  )

venn.diagram(
  main = "patients from my df and Samples",
  x = patients_list[c(1,3)],
  category.names = c("sample_d_patients", "d_patients"),
  filename = "results/related_metabolon2.0_patients/patient_venn.png",
  disable.logging = TRUE,
  output=TRUE,
  lwd = 2,
  lty = 'blank',
  fill = c("#E69F00", "#56B4E9"),
  cex = .9,
  fontface = "italic",
  cat.dist = c(-0.05, -0.05)
  )
#- - -
# ^ where family ID in famstatus not missing, no overlap between sample_d and my d
# ^ when not filtering out missing FAMILYID, only 105 overlap between these
#- - -



# JOIN FAM_D TO SAMPLE_D --------------------------------------------------

sample_d_new <- left_join(sample_d, 
                          fam_d %>% select(-matches("patId\\d")), 
                          by = join_by(MACID == patientId_new)
                          )
nrow(sample_d) == nrow(sample_d_new)
all(sample_d$MACID == sample_d_new$MACID)

# ANALYZE PATIENT FAMILIES ------------------------------------------------

# these patients not in my data
setdiff(sample_d$MACID , fam_d$patientId_new)
# "T12001"  "DO1-023" "UK223"   "SD2"     "UK232"   "UK01"    "DO1-31"  "MCT011"  "DO1-001" "MCT008"  "UK225"  
# "cy1"     "cy5"     "cy3"     "cy4"     "UK224"   "cy2"     "SD1"     "UK252"   "UK120"   "UK129"   "UK02"   
# "UK03"    "UK04"    "UK05"    "UK06"    "UK07"    "UK10"    "UK68"    "UK69"    "UK70"    "UK180"   "UK181"  
# "UK182"   "UK81"    "UK149"   "UK150"   "UK151"   "UK152"   "UK153"   "UK154"   "UK155"   "UK156"   "UK157"  
# "UK158"   "UK159"   "UK160"   "UK161"   "UK162"   "UK163"   "UK164"   "UK165"   "UK166"   "UK167"   "UK168"  
# "UK169"   "UK170"   "UK171"   "UK172"   "UK173"  

# how many patients in sample_d not in my d
no_match_patients <- sample_d_new$MACID[ is.na(sample_d_new$patientId) ]
no_match_patients %>% unique() %>% length()
# ^ 60 patients don't have match in my dataset

setdiff(sample_d$MACID , fam_d$patientId_new)

# what's up with these patients?
all(no_match_patients == setdiff(sample_d$MACID , fam_d$patientId_new))
no_match_patients 
# ^ mostly UK patients, is it all UK patients in sample_d?

setdiff(
        grep("UK", value = TRUE, sample_d_new$MACID), 
        no_match_patients
)
# ^ not all UK patients are missing from DB

# check other patID variables
pat.table <- sapply(no_match_patients, function(pat) { # for each missing patient

  # for each patient ID column, check if the missing patient is in that column 
  sapply(fam_d[, grep("pat", names(fam_d))], function(pat_col){
    
    return(pat%in%pat_col)
    
  })
})
pat.table %>% any()
# ^ the patients who are not in my data are not found in other patient ID columns

# where no family information assume no relation except wehre notes has "cyp"
cyp_notes <- grep("cyp", sample_d_new$notes, value = TRUE) %>% 
  unique()

# new family data for these relations
cyp_pat_relations <- matrix(ncol = 4, nrow = 5, byrow = TRUE,
                            dimnames =  list(NULL, c("FAMILYID_PROBANDPW", "FAMILY_STATUS", "FAMILYID", "PROBANDPW") ), 
                            data = rep(c("J38-CYP", "related", "J38", "CYP"), 5) )

sample_d_new[ sample_d_new$notes %in% cyp_notes , c("FAMILYID_PROBANDPW", "FAMILY_STATUS", "FAMILYID", "PROBANDPW")] <- cyp_pat_relations

# for where no ```FAMILYID_PROBANDPW FAMILY_STATUS FAMILYID PROBANDPW```, make up some
# ### permute \\D\\d\\d and use one not already in dataset
# faux_family_IDs <- expand.grid(LETTERS, 1:9, 1:9) %>%
#   apply(., 1, function(row) paste0(row, collapse = ""))
# 
# faux_probandPWs <-  expand.grid(LETTERS, LETTERS, LETTERS ) %>%
#   apply(., 1, function(row) paste0(row, collapse = "")) 
# 
# ### filter out existent ones
# faux_probandPWs <- faux_probandPWs %>% 
#   keep(~ !. %in% sample_d_new$PROBANDPW)
# 
# faux_family_IDs <- faux_family_IDs %>% 
#   keep(~ !. %in% sample_d_new$FAMILYID)
# 
# ## fill family data where missing
# ### where is family data missing
# family_cols <- c("FAMILYID_PROBANDPW", "FAMILY_STATUS", "FAMILYID", "PROBANDPW")
# sample_d_new %>% 
#   filter(if_any(all_of(family_cols), ~is.na(.))) %>%
#   select(all_of(family_cols)) %>% 
#   nrow()
# # ^ 169 rows where missing family info, but still family_status for some
# 
# # fill where missing with faux values
# sample_d_new$FAMILYID[is.na(sample_d_new$FAMILYID)] <- faux_family_IDs[1:length(sample_d_new$FAMILYID[is.na(sample_d_new$FAMILYID)])]
# sample_d_new$PROBANDPW[is.na(sample_d_new$PROBANDPW)] <- faux_probandPWs[1:length(sample_d_new$PROBANDPW[is.na(sample_d_new$PROBANDPW)])]
# 
# sample_d_new %>%
#   filter(if_any(all_of(family_cols), ~is.na(.))) %>%
#   select(all_of(family_cols))
# 
# # make familyID_probandpw again
# sample_d_new$FAMILYID_PROBANDPW[is.na(sample_d_new$FAMILYID_PROBANDPW)] <- paste(
#   sample_d_new$FAMILYID[is.na(sample_d_new$FAMILYID_PROBANDPW)],
#   sample_d_new$PROBANDPW[is.na(sample_d_new$FAMILYID_PROBANDPW)],
#   sep = "-"
# )

family_cols <- c("FAMILYID_PROBANDPW", "FAMILY_STATUS", "FAMILYID", "PROBANDPW")
sample_d_new %>%
  filter(if_any(all_of(family_cols), ~is.na(.))) %>%
  select(all_of(family_cols))

# fill where missing with new, unique values
## check if NA\\d\\d ok
lapply(sample_d_new[, c("FAMILYID", "PROBANDPW")], function(x)grep("NA", x, value = TRUE))

new_family_IDs <- make.unique(rep("NA_FID", length(sample_d_new$FAMILYID[is.na(sample_d_new$FAMILYID)]) + 1), sep = "")
new_family_IDs <- new_family_IDs[-1]

new_proband_PWs <- make.unique(rep("NA_PPW", length(sample_d_new$PROBANDPW[is.na(sample_d_new$PROBANDPW)]) + 1), sep = "")
new_proband_PWs <- new_proband_PWs[-1]

# fill missing with new
sample_d_new$FAMILYID[is.na(sample_d_new$FAMILYID)] <- new_family_IDs
sample_d_new$PROBANDPW[is.na(sample_d_new$PROBANDPW)] <- new_proband_PWs

# make familyID_probandPW again
sample_d_new$FAMILYID_PROBANDPW[is.na(sample_d_new$FAMILYID_PROBANDPW)] <- paste(
    sample_d_new$FAMILYID[is.na(sample_d_new$FAMILYID_PROBANDPW)],
    sample_d_new$PROBANDPW[is.na(sample_d_new$FAMILYID_PROBANDPW)],
    sep = "-"
  )

sample_d_new %>%
  filter(if_any(all_of(family_cols), ~is.na(.))) %>%
  select(all_of(family_cols))

# fill family_status with unsure where missing
sample_d_new$FAMILY_STATUS[is.na(sample_d_new$FAMILY_STATUS)] <- "unknown"

#///////////#///////////
#///////////#///////////
# save the original sample_d with family info
write.csv(sample_d_new ,"results/related_metabolon2.0_patients/Samples_02.csv")
#///////////#///////////
#///////////#///////////






