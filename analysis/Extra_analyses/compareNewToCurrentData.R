#+ this script is a contiuation compareDEM.R exploration. that dealt with just DEM, this investigates all tables
#+ compare currrent in use datasets to ones newly donwloaded. 
#+ check that they have same amount of data. if not, then db glitched out with original downloading. 


pacman::p_load(dplyr, tidyr, ggplot2,reshape2, purrr, stringr, HDF5Array, parallel)
rm(list = ls()); gc()


# LOAD DATA----
wd_path = "~/Documents/LMRI/PRJ2024001/version015/"

# my functions
source(sprintf("%sanalysis/00_MyFunctions.R", wd_path))

# get files in directory and extract csvs
orig_path <- "/Users/christiananderson/Library/CloudStorage/Egnyte-lmri/Shared/LMRI/DS_Group/PRJ2024001/raw_data/06" 
new_path <- "/Users/christiananderson/Library/CloudStorage/Egnyte-lmri/Shared/LMRI/DS_Group/PRJ2024001/raw_data/07"

# ****** FUNCTION: grab csv files from directory and clean them
load_and_clean_tables <- function(path_to_directory){
  
  # list of filenames
  files <- dir(path_to_directory)
  csv_files <- files[grep("\\.csv$", files)] 
  
  # load data into list
  raw_dfs <- lapply(paste0(path_to_directory,'/',csv_files), readr::read_csv)
  
  # name items in list 
  names(raw_dfs) <- lapply(csv_files, function(x) substring(x,1, nchar(x) - 4) )
  
  dfs <- raw_dfs
  
  # exclude VAA
  dfs$VAA_emmes <- NULL
  
  # clean chars off
  dfs <- mclapply(dfs, clean.it)
  
  dfs <- map(dfs,janitor::clean_names)
  return(dfs)
}

# for each directory, load tables to list
dfs <- mclapply(list(orig = orig_path, new = new_path), load_and_clean_tables)

## remove VAA, not comparable
dfs$new$VAA <- NULL

# are all dfs accounted for?
sapply(dfs, function(x) sort(names(x))) %>% 
  {.[,1] == .[,2]} %>% 
  all

# COMPARE EACH NEW TO OLD TABLE----

# are all vars there to compare??
Map(function(x,y) {
  # all(x%in%y) & all(y%in%x)
  c(orig_uniques = setdiff(x,y), new_uniques = setdiff(y,x))
  }, map(dfs$orig, names), map(dfs$new, names))
# $DEM
# orig_uniques new_uniques1 new_uniques2 new_uniques3 new_uniques4 new_uniques5 new_uniques6 new_uniques7 
# "site" "patient_id"    "pat_id2"    "pat_id3"    "pat_id4"   "visit_id"        "eye" "created_at"

dfs$orig$DEM[, c("site")] <- NULL
dfs$orig$DEM$patient_id <- dfs$orig$DEM$patid
dfs$new$DEM[, c("patid", "pat_id2","pat_id3","pat_id4","visit_id","eye","created_at")] <- NULL

# dim comparisons
map(dfs,~{sapply(., dim)})
# $orig
# CF-FFA  DEM DIAGNOSIS ENRNHOR ENRNHOS FAFGRADE FAMSTATUS    MH OCTMAIN OCTMTCIRRUS OCTMTSPECTRALIS OCTMTSTRATUS
#  31357 1173     31357   31357   31357    31357     31357 31357   31357       31357           31357        31357
#    98    7        10      28      30       21        11    56      30          22              32           24
# $new
# CF-FFA   DEM DIAGNOSIS ENRNHOR ENRNHOS FAFGRADE FAMSTATUS    MH OCTMAIN OCTMTCIRRUS OCTMTSPECTRALIS OCTMTSTRATUS
#  31432 31432     31432   31432   31432    31432     31432 31432   31432       31432           31432        31432
#     98     6        10      28      30       21        11    56      30          22              32           24
##              ^^ way more dem in new download



# make dfs only id cols
dfs_ids <- lapply(dfs, function(df_list){
  lapply(df_list, function(df){
    df %>% 
      select(any_of(c("visit_id", "eye", "patient_id")))
  })  
})

rm("df_name")

# data availability comparison
for(df_name in names(dfs_ids$orig)){

  if( 
    all(c("visit_id", "eye", "patient_id") %in% names(dfs_ids$orig[[df_name]]))
    ){
    suppressMessages(
      # inner_join(dfs_ids$new[[df_name]],dfs_ids$orig[[df_name]] , by = c("month", "eye", "patient_id"))%>% nrow()
      n <- inner_join(dfs_ids$new[[df_name]],dfs_ids$orig[[df_name]] )%>% nrow()
    )
    }else if(
      all(c("eye", "patient_id") %in% names(dfs_ids$orig[[df_name]]))
      ){
      suppressMessages(
          # inner_join(dfs_ids$new[[df_name]],dfs_ids$orig[[df_name]] , by = c("eye","patient_id"))%>% nrow()
          n <- inner_join(dfs_ids$new[[df_name]],dfs_ids$orig[[df_name]] )%>% nrow()
      )
      } else if(
        all(c( "patient_id") %in% names(dfs_ids$orig[[df_name]]))
        ){
        suppressMessages(
          # inner_join(dfs_ids$new[[df_name]],dfs_ids$orig[[df_name]] , by = c("patient_id"))%>% nrow()
           n <- inner_join(dfs_ids$new[[df_name]],dfs_ids$orig[[df_name]] )%>% nrow()
        )
      }else{
        n <- NA_integer_
      }
  cat(sprintf("%s: %s\n", df_name, n))
}
# CF-FFA: 31327
# DEM: 10677
# DIAGNOSIS: 31327
# ENRNHOR: 31327
# ENRNHOS: 31327
# FAFGRADE: 31327
# FAMSTATUS: 31327
# MH: 31327
# OCTMAIN: 31327
# OCTMTCIRRUS: 31327
# OCTMTSPECTRALIS: 31327
# OCTMTSTRATUS: 31327

