#+ idiot check that missing yrofbrth and sex in version014 13_routeModeling.Rmd is rooted in database issue .
#+ compare new download of dem from db to sce and ensure that they maatch up. if they do, this is an issue with the database
#+ 


pacman::p_load(dplyr, tidyr, ggplot2,reshape2, purrr, stringr, HDF5Array)
rm(list = ls()); gc()

source("/Users/canderson/Documents/LMRI/PRJ2024001/version015/analysis/00_MyFunctions.R")

# LOAD DATA----

wd_path = "~/Documents/LMRI/PRJ2024001/version015/"

# load experiment
sce <- HDF5Array::loadHDF5SummarizedExperiment(paste0(wd_path, "processed_data/sce005.h5"))

# downloaded from 'new db' bookmark
new_DEM <-read.csv(unz("/Users/christiananderson/Library/CloudStorage/Egnyte-lmri/Shared/LMRI/DS_Group/PRJ2024001/raw_data/DEM1.zip", "DEM.csv"))

# load dem used to construct sce
curr_DEM <- read_csv('/Users/christiananderson/Library/CloudStorage/Egnyte-lmri/Shared/LMRI/DS_Group/PRJ2024001/raw_data/06/DEM.csv')

sce_DEM <- select(data.frame(colData(sce)), patID, SEXF, YROFBRTH) %>% distinct()

# put them in list and preprocess
DEMS <- list(curr_DEM = curr_DEM,sce_DEM = sce_DEM, new_DEM = new_DEM) %>% map(clean.it)

sapply(DEMS, dim) %>% 
  {.[,1]-.[,3]}
# -30259     -6

# ((((((((((((((((((()))))))))))))))))))
# ((((((((((((((((((()))))))))))))))))))
# ((((((((((((((((((()))))))))))))))))))
# ((((((((((((((((((()))))))))))))))))))
# COMPARE NEW DOWNLOAD TO SCE ----
# ((((((((((((((((((()))))))))))))))))))
# ((((((((((((((((((()))))))))))))))))))
# ((((((((((((((((((()))))))))))))))))))
# ((((((((((((((((((()))))))))))))))))))

dem <- select(DEMS$new_DEM, patientId, GENDER, YROFBRTH) %>%
  mutate(patID = patientId) %>%
  mutate(SEXF = GENDER=="Female") %>%
  select(-GENDER) %>%
  mutate(YROFBRTH = as.numeric(YROFBRTH)) %>%
  distinct()

inn <- left_join(select(data.frame(colData(sce)), patID, SEXF, YROFBRTH) %>% distinct(),
                dem,
                by = 'patID',
                suffix = c("_sce", "_db"))
inn <- inn[, c("patID", sort(setdiff(names(inn), "patID")))]

# ****** FUNCTION: take any number of vectors and calculate if all missing across rows
all_is.na <- function(...){
  args <- list(...)
  rowAlls(sapply(args,is.na))
}
# inn[all_is.na(inn$YROFBRTH_db, inn$YROFBRTH_sce, inn$SEXF_db), c("YROFBRTH_db", 'YROFBRTH_sce', "SEXF_db")]


# where not all vars -patID are missing and where vars are not identical
inn %>%
  # filter(if_any(-patID, ~!is.na(.))) %>%
  # filter(!identical(SEXF_db,SEXF_sce ) | !identical(YROFBRTH_db, YROFBRTH_sce) & (!all_is.na(SEXF_db,SEXF_sce) |!all_is.na(YROFBRTH_db, YROFBRTH_sce) ))
  ungroup() %>%
  select(-patID) %>%
  mice::md.pattern(rotate.names = T)
# 3133 obs where no sex/yrbrth in sce but in db

# ((((((((((((((((((()))))))))))))))))))
# ((((((((((((((((((()))))))))))))))))))
# ((((((((((((((((((()))))))))))))))))))
# ((((((((((((((((((()))))))))))))))))))
# COMPARE ORIGINAL DOWNLOAD TO SCE ----
# ((((((((((((((((((()))))))))))))))))))
# ((((((((((((((((((()))))))))))))))))))
# ((((((((((((((((((()))))))))))))))))))
# ((((((((((((((((((()))))))))))))))))))
dem <- select(DEMS$curr_DEM, PATID, GENDER, YROFBRTH) %>%
  mutate(patID = PATID) %>%
  mutate(SEXF = GENDER=="Female") %>%
  select(-GENDER, -PATID) %>%
  mutate(YROFBRTH = as.numeric(YROFBRTH)) %>%
  distinct()


inn <- inner_join(select(data.frame(colData(sce)), patID, SEXF, YROFBRTH) %>% distinct(),
                 dem,
                 by = 'patID',
                 suffix = c("_sce", "_db"))

inn <- inn[, c("patID", sort(setdiff(names(inn), "patID")))]

# ****** FUNCTION: take any number of vectors and calculate if all missing across rows
all_is.na <- function(...){
  args <- list(...)
  rowAlls(sapply(args,is.na))
}
# inn[all_is.na(inn$YROFBRTH_db, inn$YROFBRTH_sce, inn$SEXF_db), c("YROFBRTH_db", 'YROFBRTH_sce', "SEXF_db")]


# where not all vars -patID are missing and where vars are not identical
inn %>%
  # filter(if_any(-patID, ~!is.na(.))) %>%
  # filter(!identical(SEXF_db,SEXF_sce ) | !identical(YROFBRTH_db, YROFBRTH_sce) & (!all_is.na(SEXF_db,SEXF_sce) |!all_is.na(YROFBRTH_db, YROFBRTH_sce) ))
  ungroup() %>%
  select(-patID) %>%
  mice::md.pattern(rotate.names = T)
# 3133 obs where no sex/yrbrth in sce but in db
