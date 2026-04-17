
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

########################
#### paths and directories ####
########################
wd_path = "~/Documents/LMRI/PRJ2024001/version015/"

if(!dir.exists(paste0(wd_path, "results/compare_to_master/"))) {
  dir.create(paste0(wd_path, "results/compare_to_master/"), recursive = TRUE)
  
}
########################
#### load data #####
########################

# main d
d <- readRDS(paste0(wd_path, "processed_data/CF_FFA_DEM_DIAGNOSIS_ENRNHOR_ENRNHOS_FAFGRADE_FAMSTATUS_MH_OCTMAIN_OCTMTCIRRUS_OCTMTSPECTRALIS_OCTMTSTRATUS_VAA_lea_01.rds"))

# info 
info <- read.xlsx(paste0(wd_path, "processed_data/database_grading_key_04.xlsx"))


# master excel file
master_xl_path <- "/Users/christiananderson/Library/CloudStorage/Egnyte-lmri/Shared/LMRI/DS_Group/PRJ2024001/raw_data/GRADING_LOCAL_MacTel_GradingData_ready_DEC2021_edited.xlsx"

## getting data from sheets
sheets <- openxlsx::getSheetNames(master_xl_path)
master_dfs <- lapply(sheets, openxlsx::read.xlsx, xlsxFile=master_xl_path)

## assigning names to each data frame
names(master_dfs) <- c("FAF", "OCTMAIN", "OCTMTSPECTRALIS", "OCTMTSTRATUS", "OCTMTCIRRUS", "CF", "FFA") 

########################
#### Master CF and FFA versus CF_FFA ####
########################
# CF and FFA ###

## CF 
CF <- master_dfs$CF

# fix names
CF_names <- names(CF)
CF_first_row <- CF[1,] %>% 
  as.character()

names_d <- cbind(CF_names, CF_first_row)

# where name is "hide" or HRDEXUDATE, make value in first row
new_names <- ifelse(CF_names %in% c("hide", "HRDEXUDATE"), names_d[, 2], names_d[,1])

names(CF) <- new_names

CF <- CF[-1,]

head(CF)

# make duplicate visitID unique; second one with same name
names(CF)[which(names(CF) == "visitID")[2]] <- "old_visitID" 

# remove useless cols
CF <- CF %>% 
  select(!c(X94, COMMCFFFA, `noAbnormality_HISTORICAL;.do.not.import` , `PHOTODT-CF` , `retinal_haem_at`))

# what values to clean for
lapply(CF, unique)

# clean for 7 8 9 and NULL values
values_to_NA <- c( "null", "NULL", "Null", "na" , "Na","", " ", "\t", "\n", "  ",  "777", "7777", "888", "8888", "999", "9999", "7","8","9", "removed")

## what cols have these  to clean for 
few_unique <- names(CF)[sapply(CF, function(x) length(unique(x)) < 30)]

lapply(CF[, few_unique], unique)

# clean for to NA
CF <- CF %>% 
  mutate(across(all_of(few_unique), 
                ~ifelse(. %in% values_to_NA, NA, .)
  )
  )

# FFA
FFA <- master_dfs$FFA

# fix names
FFA_names <- names(FFA)
FFA_first_row <- FFA[1,] %>% 
  as.character()

names_d <- cbind(FFA_names, FFA_first_row)

duplicate_names <- FFA_names[table(FFA_names )>1] %>% 
  unique()

# where name is "hide" or HRDEXUDATE, make value in first row
new_names <- ifelse(FFA_names %in% duplicate_names, FFA_first_row, FFA_names)

names(FFA) <- new_names

FFA <- FFA[-1,]

head(FFA)

# remove useless cols
FFA <- FFA %>% 
  select(!c(COMMCFFFA, `noAbnormality_HISTORICAL;.do.not.import`))

# what values to clean for
lapply(FFA, unique)

## what cols have these  to clean for 
few_unique <- names(FFA)[sapply(FFA, function(x) length(unique(x)) < 30)]

lapply(FFA[, few_unique], unique)

# clean for to NA
FFA <- FFA %>% 
  mutate(across(all_of(few_unique), 
                ~ifelse(. %in% values_to_NA, NA, .)
  )
  )

# join the two together 
CF_FFA <- full_join(CF, FFA)


# filter out controls and family visits
CF_FFA <- CF_FFA %>% 
  filter(
    !VisitID %in% c("control", "Control")
  )

## Rename CF_FFA names to those in d ##
# d subset at the cf-ffa vars
grep_string <- info$variable[info$dataset =="CF_FFA" ] %>% 
  gsub("-", "2345678", .) %>% 
  paste(., collapse = "|")

d_cf.ffa.names <- grep(grep_string, names(d), value = TRUE)

# grab unique variables (excluding zone numbers)
d_cf.ffa_name_search <- d_cf.ffa.names %>% 
  gsub("0\\d$", "0\\\\d", .) %>% 
  unique()


# match vars in d_cf.ffa to vars in CF_FFA
name_match_l <- sapply(d_cf.ffa_name_search, function(name_search){
  grepl(name_search, names(CF_FFA)) %>% any()
}) 

# names to rename in CF_FFA
to_rename <- names(name_match_l)[!name_match_l]
# "RPE_LK0\\d"       "CF_FFA_GRADINGDT" "HRDEXUDATE"       "PHOTODTCF"        "PHOTODTFFA"       "TYP_LEAK_RPE"    
# "NOABNORMCF"       "NOABNORMFL" 

old_new_names_d <- data.frame(old = c("rpe_leakage_", "type_leakage_rpe" ), 
                              new = c("RPE_LK", "TYP_LEAK_RPE"))

# rename CFFFA to d names
for (i in 1:nrow(old_new_names_d)){
  row <- old_new_names_d[i,]
  
  new_names <- names(CF_FFA)[grep(row[1], names(CF_FFA))] %>% 
    gsub(row[1], row[2],. )
  
  names(CF_FFA)[grep(row[1], names(CF_FFA))] <- new_names
  
}

setdiff(d_cf.ffa.names, names(CF_FFA))
# "CF_FFA_GRADINGDT" "HRDEXUDATE"       "PHOTODTCF"        "PHOTODTFFA"       "NOABNORMCF"       "NOABNORMFL" 
# ^ not too important

shared_vars <- intersect(d_cf.ffa.names, names(CF_FFA))[order(intersect(d_cf.ffa.names, names(CF_FFA)))]

# df of var observation counts and source
count_d <- data.frame(d  = colSums(!is.na(d[, shared_vars])), 
                      CF_FFA = colSums(!is.na(CF_FFA[, shared_vars]))) %>% 
  rownames_to_column()


# tidy it
count_d_tidy <- melt(count_d)

names(count_d_tidy) <- c("variable", "dataset", "count_observed")

count_d_tidy$variable <- factor(count_d_tidy$variable, levels = shared_vars)

# barplot of var counts
count_d_tidy %>%
  ggplot(aes(y =variable, x = count_observed, fill = dataset))+
  geom_bar(orientation = "y", stat = "identity", position = "dodge")+
  ggtitle("Comparing master excel file CF_FFA to database CF_FFA")
ggsave(paste0(wd_path, "results/compare_to_master/", "comparing_CF_FFA_vars_counts.png"), height = 15, width = 10)


###########################################################################
#### TAKE REPRESENTATIVE VAR FROM EACH DATAGROUP AND COMPARE TO MAIN ####
###########################################################################

# fix master_dfs names which have "hide"
new_master_dfs <- lapply(master_dfs, function (df){
  
  # if hide, grab from first row
  a_df_names_d <- cbind(names(df), as.character(df[1,]))
  a_df_new_names <- ifelse(names(df) %in% c("HRDEXUDATE", "RPE_LK","hide", paste0("X", as.character(1:100))), a_df_names_d[,2], a_df_names_d[,1])
  
  
  # assign new names
  names(df) <- a_df_new_names
  
  # remove first row
  df <- df[-1, ]
  
  return(df)
})


df_name = "OCTMTSTRATUS"
# add machine names to octmt datasets
# where name is MT0, add machine name to beginning of that name
for(df_name in c("OCTMTSPECTRALIS", "OCTMTSTRATUS","OCTMTCIRRUS") ) {
  
  machine_name <- substr(df_name, 6, nchar(df_name)) %>% 
    tolower()
  
  # stratus names different
  if(df_name == "OCTMTSTRATUS"){
    
    # rename macthickness to MT0...
    new_MT_names <- names(new_master_dfs[[df_name]])[grep("mac_thickness_0", names(new_master_dfs[[df_name]]))] %>% 
      gsub("mac_thickness_0", "MT0", .)
    
    # assign new names
    names(new_master_dfs[[df_name]])[grep("mac_thickness_0", names(new_master_dfs[[df_name]]))] <- paste0(machine_name, "_", new_MT_names)
    
  }else{
    
    names(new_master_dfs[[df_name]])[grep("MT0", names(new_master_dfs[[df_name]]))] <- paste0(machine_name, "_", names(new_master_dfs[[df_name]])[grep("MT0", names(new_master_dfs[[df_name]]))] )
    
  }
}


# FAF, CF_FFA, OCTMAIN, OCTMTSTRATUS, OCTMTCIRRUS, OCTMTSPECTRALIS
sapply(new_master_dfs, names) 

# data frame of var names from both datasets
names_to_compare_d <- data.frame(d = c("RPE_LK05", "ICL05",  "BRV05" ,"cirrus_MT05", "spectralis_MT05", "stratus_MT05", "CTRPOINT_INVOLVD", "AF_INCREASED_1DD"),
                                 master  = c("rpe_leakage_05", "ICL05", "BRV05" ,"cirrus_MT05", "spectralis_MT05", "stratus_MT05", "CTRPOINT_INVOLVD", "AF_INCREASED_1DD")
)

# make df selecting those variables from the master_dfs
master_select_vars <- list()
for(var in names_to_compare_d$master){
  
  # which df has var (returned as list of 1)
  para_df <- new_master_dfs[
    #search dfs for one with var
    sapply(new_master_dfs, function(a_df){ 
      any(grepl(var, names(a_df)))
    }, USE.NAMES = FALSE)
  ]
  
  df <- para_df[[1]]
  
  # grab variable as named vector
  master_select_vars[[var]] <- df[,var]
  
}

# make all items same length
max_length <- sapply(master_select_vars, length, USE.NAMES = FALSE) %>% max()
master_select_vars_lengthened <- lapply(master_select_vars, function(item){
  
  item_length <- length(item)
  if(item_length < max_length){
    item <- c(item, rep(NA, (max_length - item_length) ))
  }
  return(item)
})

# dataframe of select variables to compare to main d
master_select_vars_d <- master_select_vars_lengthened %>% 
  as.data.frame()

sub_master <- master_select_vars_d 

# rename to same as d
names(sub_master) <- names_to_compare_d$d

# remove missing value synonyms
sub_master <- sub_master %>% 
  mutate(across(everything(), 
                ~ifelse(. %in% values_to_NA, NA, .))) %>% 
  mutate(across(everything(), 
                ~as.numeric(.)))


# check for what values now
lapply(sub_master, function(v) head(unique(v), 10))

# select same vars from d
sub_d <- d %>% 
  select(all_of(names_to_compare_d$d)) %>% 
  mutate(across(everything(), ~ as.numeric(.)))

# dataframe of count observed tidy
select_vars_counts_d <- data.frame(
  d = colSums(!is.na(sub_d)),
  master = colSums(!is.na(sub_master))
) %>% 
  rownames_to_column() %>% 
  melt()

names(select_vars_counts_d) <- c("variable", "dataset", "count_observed")

# plot counts bar
select_vars_counts_d %>% 
  ggplot(aes(y =variable, x = count_observed, fill = dataset))+
  geom_bar(orientation = "y", stat = "identity", position = "dodge")+
  ggtitle("Comparing representative variables from master excel database datasets")
ggsave(paste0(wd_path, "results/compare_to_master/", "comparing_representative_vars_counts.png"), height = 15, width = 10)



###########################################################################
#### WHERE PATTERNS WHERE ICL MISSING IN D ####
###########################################################################

## fix visit IDs to match
master_for_var <- CF_FFA %>% 
  filter(!VisitID %in% grep("unsch", CF_FFA$VisitID, value = TRUE)) %>% 
  mutate(visitID =case_when(
             VisitID == "Family Member Visit" ~ "000", 
             VisitID == "Baseline" ~ "000", 
             VisitID == "6 months" ~ "006", 
             VisitID == "1 Year" ~ "012",
             VisitID == "2 Year" ~ "024",
             VisitID == "3 Year" ~ "036",
             VisitID == "4 Year" ~ "048",
             VisitID == "5 Year" ~ "060",
             VisitID == "6 Year" ~ "072",
             VisitID == "8 Year" ~ "072", 
             VisitID == TRUE ~ VisitID
           )
         )
####################
# ICL
####################
# make df of ICL from master and d
master_ICL <- master_for_var %>% 
  select(patientID, eye_num = EYE, VisitID, matches("ICL0")) 

# rename ICL cols
names(master_ICL)[grep("ICL", names(master_ICL))] <- names(master_ICL)[grep("ICL", names(master_ICL))] %>% 
  paste("master", ., sep = "_")
# rename patientID
names(master_ICL)[grep("patientID", names(master_ICL))] <- "patID"


# d subsetted for ICL cols
d_ICL <- d %>% 
  select(patID, eye_num, visitID, matches("ICL0")) 

# rename ICL cols
names(d_ICL)[grep("ICL", names(d_ICL))] <- names(d_ICL)[grep("ICL", names(d_ICL))] %>% 
  paste("d", ., sep = "_")

joined_ICL <- left_join(d_ICL, master_ICL) %>% 
  arrange(visitID, patID)

# bar plot of count observed
joined_ICL %>%
  summarise(across(everything(), ~ sum(!is.na(.)))) %>% 
  melt(., value.name =  "count", variable.name = "variable") %>% 
  ggplot(aes(count, variable, fill = variable))+
    geom_bar(stat = "identity", position = "dodge", orientation = "y", show.legend = FALSE)

# image dataset 
## function to make datasets binary for NA
make_binary_matrix <- function(data = data.frame()){
  new_data <- data %>%
    mutate(across(everything(),
                  ~ifelse(is.na(.), 0, 1)))

  new_data_matrix <- as.matrix(new_data)

  return(new_data_matrix)
}

joined_ICL_binary <- make_binary_matrix(joined_ICL)
rownames(joined_ICL_binary) <- joined_ICL$visitID %>% make.unique(.)

# pdf(paste0(wd_path, "results/compare_to_master/", "ICL_cols_compare.pdf"), height =  100, width = 100)
# Heatmap(
#   joined_ICL_binary,
#   cluster_rows = FALSE,
#   cluster_columns = FALSE,
#   col = c("white", "black"),
#   name = "observed",
#   show_row_names = TRUE,
#   row_names_gp = gpar(fontsize = 1),
# )
# dev.off()


# where missing in d and not in CF_FFA, what subject
d_cols <- grep("d_", names(joined_ICL), value = TRUE)
master_cols <- grep("master_", names(joined_ICL), value = TRUE)

# rows where no d ICL and are master ICL 
idx <- rowSums(!is.na(joined_ICL[ , d_cols])) == 0 & rowSums(!is.na(joined_ICL[ , master_cols])) >0

# what patients who have data in master and not in d for ICL
patients_in_d_missing_ICL <- joined_ICL[ idx , c("patID", "visitID")] %>% 
  unique()

write_delim(patients_in_d_missing_ICL, file = paste0(wd_path, "results/compare_to_master/", "patients_in_d_missing_ICL.txt"), delim = ",")

# rows wehre d ICL and no master ICL
idx <- rowSums(!is.na(joined_ICL[ , d_cols])) > 0 & rowSums(!is.na(joined_ICL[ , master_cols])) == 0

# what patients who have data in master and not in d for ICL
patients_in_master_missing_ICL <- joined_ICL[ idx , c("patID", "visitID")] %>% 
  unique()

write_delim(patients_in_master_missing_ICL, file = paste0(wd_path, "results/compare_to_master/", "patients_in_master_missing_ICL.txt"), delim = ",")

####################
# IR
####################

# # make df of IR from master and d
# master_IR <- master_for_var %>% 
#   select(patientID, eye_num = EYE, VisitID, matches("IR0")) 
# 
# # rename IR cols
# names(master_IR)[grep("IR", names(master_IR))] <- names(master_IR)[grep("IR", names(master_IR))] %>% 
#   paste("master", ., sep = "_")
# 
# # rename patientID
# names(master_IR)[grep("patientID", names(master_IR))] <- "patID"
# 
# 
# # d subsetted for IR cols
# d_IR <- d %>% 
#   select(patID, eye_num, visitID, matches("IR0")) 
# 
# # rename IR cols
# names(d_IR)[grep("IR", names(d_IR))] <- names(d_IR)[grep("IR", names(d_IR))] %>% 
#   paste("d", ., sep = "_")
# 
# joined_IR <- left_join(d_IR, master_IR) %>% 
#   arrange(visitID, patID)
# 
# # bar plot of count observed
# joined_IR %>%
#   summarise(across(everything(), ~ sum(!is.na(.)))) %>% 
#   melt(., value.name =  "count", variable.name = "variable") %>% 
#   ggplot(aes(count, variable, fill = variable))+
#   geom_bar(stat = "identity", position = "dodge", orientation = "y", show.legend = FALSE)
# 
# # image dataset 
# ## function to make datasets binary for NA
# make_binary_matrix <- function(data = data.frame()){
#   new_data <- data %>%
#     mutate(across(everything(),
#                   ~ifelse(is.na(.), 0, 1)))
#   
#   new_data_matrix <- as.matrix(new_data)
#   
#   return(new_data_matrix)
# }
# 
# joined_IR_binary <- make_binary_matrix(joined_IR)
# rownames(joined_IR_binary) <- joined_IR$visitID %>% make.unique(.)
# 
# # pdf(paste0(wd_path, "results/compare_to_master/", "IR_cols_compare.pdf"), height =  100, width = 100)
# Heatmap(
#   joined_IR_binary,
#   cluster_rows = FALSE,
#   cluster_columns = FALSE,
#   col = c("white", "black"),
#   name = "observed",
#   show_row_names = TRUE,
#   row_names_gp = gpar(fontsize = 1),
# )
# # dev.off()



###########################################################################
#### OLD PLOTTING ####
###########################################################################

##### Image datasets ####
# #> function to make datasets binary for NA
# make_binary_matrix <- function(data = data.frame()){
#   new_data <- data %>% 
#     mutate(across(everything(), 
#                   ~ifelse(is.na(.), 0, 1)))
#   
#   new_data_matrix <- as.matrix(new_data)
#   
#   return(new_data_matrix)
# }
# 
# # make dfs binary
# CF_FFA_binary <- CF_FFA %>% 
#   make_binary_matrix(.)
# 
# rownames(CF_FFA_m) <- 1:nrow(CF_FFA_m)
# 
# 
# 
# # make d at CF_FFA vars binary
# d_cf.ffa_binary<- d_cf.ffa %>% 
#   make_binary_matrix(.)
# 
# rownames(d_at_cf.ffa_binary) <- 1:nrow(d_at_cf.ffa_binary)
# 
# 
# 
# 
# #> plot
# h = 100
# w = 100
# yfnt= 1
# xfnt= 5
# # 
# # # master data
# # pdf(paste0(wd_path, "results/compare_to_master/master_CF_FFA_heatmap.pdf"), height = h, width = w )
# # Heatmap(CF_FFA_m ,
# #         show_heatmap_legend = TRUE,
# #         use_raster = F,
# #         column_title  = "CF_FFA columns\n(1 = observed)",
# #         row_title  = "observations",
# #         col = c("white" ,  "navy"),
# #         na_col = "white", 
# #         cluster_rows = FALSE,
# #         cluster_columns = FALSE,
# #         show_row_names = TRUE,
# #         row_names_gp = gpar(fontsize = yfnt),
# #         column_names_gp = gpar( fontsize = xfnt),
# #         show_column_names  = TRUE,
# #         column_names_rot = 45
# # )
# # dev.off()
# # 
# # # production data
# # pdf(paste0(wd_path, "results/compare_to_master/production_heatmap.pdf"), height = h, width = w )
# # Heatmap(d_in_use_CF_FFA_m ,
# #         show_heatmap_legend = TRUE,
# #         use_raster = F,
# #         column_title  = "d in use columns\n(1 = observed)",
# #         row_title  = "observations",
# #         col = c("white" ,  "darkgreen"),
# #         na_col = "white", 
# #         cluster_rows = FALSE,
# #         cluster_columns = FALSE,
# #         row_names_gp = gpar( fontsize = yfnt),
# #         column_names_gp = gpar( fontsize = xfnt),
# #         show_row_names = TRUE,
# #         show_column_names  = TRUE,
# #         column_names_rot = 45
# # )
# # dev.off()
# # 
# # 
# # # new data
# # pdf(paste0(wd_path, "results/compare_to_master/new_CF_FFA_heatmap.pdf"), height = h, width = w )
# # Heatmap(new_CF_FFA_m ,
# #         show_heatmap_legend = TRUE,
# #         use_raster = F,
# #         column_title  = "d in use columns\n(1 = observed)",
# #         row_title  = "observations",
# #         col = c("white" ,  "darkorange"),
# #         na_col = "white", 
# #         cluster_rows = FALSE,
# #         cluster_columns = FALSE,
# #         row_names_gp = gpar( fontsize = yfnt),
# #         column_names_gp = gpar( fontsize = xfnt),
# #         show_row_names = TRUE,
# #         show_column_names  = TRUE,
# #         column_names_rot = 45
# # )
# # dev.off()
# 
# #### COMPARING numbers of observations in each #
# # what variables do the three share, compare their numbers
# shared_vars <- intersect(intersect(colnames(CF_FFA_m), colnames(d_in_use_CF_FFA_m)), colnames(new_CF_FFA_m))
# 
# 
# CF_FFA_counts <- CF_FFA_m[, shared_vars] %>% 
#   colSums() %>% 
#   as.data.frame() %>% 
#   rownames_to_column("feature") %>% 
#   rename_with(~all_of(c("feature", "CF_FFA_count")))
# 
# d_in_use_CF_FFA_counts <- d_in_use_CF_FFA_m[, shared_vars]%>% 
#   colSums() %>% 
#   as.data.frame() %>% 
#   rownames_to_column("feature") %>% 
#   rename_with(~all_of(c("feature", "d_in_use_CF_FFA_count")))
# 
# new_CF_FFA_counts <- new_CF_FFA_m[, shared_vars]%>% 
#   colSums() %>% 
#   as.data.frame() %>% 
#   rownames_to_column("feature") %>% 
#   rename_with(~all_of(c("feature", "new_CF_FFA_count")))
# 
# counts_combined <- left_join(left_join(CF_FFA_counts, d_in_use_CF_FFA_counts), new_CF_FFA_counts) %>% 
#   melt()
# 
# counts_combined_averages <- counts_combined %>% 
#   group_by(variable) %>% 
#   mutate(min_count = min(value)) %>% 
#   mutate(max_count = max(value)) %>% 
#   mutate(average_count = mean(value)) %>% 
#   mutate(median_count = median(value)) %>% 
#   ungroup()
# 
# # how many observed values in each column for each CF_FFA version
# ggplot(counts_combined_averages, aes( value, feature, fill = variable))+
#   geom_bar(stat = "identity", orientation = "y", position = "dodge")+
#   theme_bw()
# # ggsave(paste0(wd_path, "results/compare_to_master/shared_vars_observed_counts_by_df.png"), height = 15, width = 10)
 
