# Wed May  7 11:44:44 2025 ------------------------------

pacman::p_load(tidyr,reshape2,purrr,tidyverse,ggplot2,gridExtra,ComplexHeatmap,viridis,zoo,mgcv,mice,ggrepel,knitr, htmltools,Metrics,psych,car,stats,SingleCellExperiment,dplyr,HDF5Array,parallel,txtplot,  modelsummary, magrittr)

rm(list = ls()); gc()



# script with all functions i have made
source("/Users/canderson/Documents/LMRI/PRJ2024001/version015/analysis/00_MyFunctions.R")

# LOAD DATA ----
wd_path = "~/Documents/LMRI/PRJ2024001/version015/"

sce <-  loadHDF5SummarizedExperiment(paste0(wd_path, "processed_data/sce006.h5"))

pth = "~/Documents/LMRI/PRJ2024001/version015/results/summary.txt"

# COUNT VARIABLES ----
cat("\n", 
    'Number of variables: ',c(names(colData(sce)), rownames(sce)) %>% length(),
    '\n',
    file = pth, append = FALSE, sep = ''
    )


# COUND OBS ----
cat("\n", 
    'Number of observations: ',ncol(sce),
    '\n',
    file = pth, append = TRUE,sep = ''
)

# COUNTS OF PATIENTS/eyes, longitudinal and cross-sectional----

# PATIETNS (((((((((((())))))))))))
pat_visits <- data.frame(colData(sce)) %>% 
  select(patID, month, rs17421627_tmem161b) %>% 
  arrange(patID, month) %>% 
  distinct()

# longitudinal patients
(x <- pat_visits %>% 
  filter(n()>1, .by = patID) %>% 
  distinct(patID) %>% 
  nrow())
# 492
cat("\n", 
    'Longitudinal Patients: ',x,
    '\n',
    file = pth, append = TRUE,sep = ''
)

# single visit patients
(x <- pat_visits %>% 
  filter(n()==1, .by = patID) %>% 
  distinct(patID) %>% 
  nrow())
# 3479
cat("\n", 
    'Single Visit Patients: ',x,
    '\n',
    file = pth, append = TRUE,sep = ''
)

# totla patients
(x = n_distinct(sce$patID))
(y  = 492+3479==x)
cat("\n", 
    'Total Patients: ',x,
    '\n','\n',
    'Single + Longitudinal = total : ',y,
    '\n',
    file = pth, append = TRUE,sep = ''
)


# longitudinal with genetics
(x <- pat_visits %>% 
  filter(!is.na(rs17421627_tmem161b)) %>% 
  filter(n()>1, .by = patID) %>% 
  distinct(patID) %>% 
  nrow())
# 359
cat("\n", 
    'Longitudinal Patients With Genetics: ',x,
    '\n',
    file = pth, append = TRUE,sep = ''
)


# single visit genetics patients
(x <- pat_visits %>% 
  filter(!is.na(rs17421627_tmem161b)) %>% 
  filter(n()==1, .by = patID) %>% 
  distinct(patID) %>% 
  nrow())
# 2463

cat("\n", 
    'Single Visit Patients With Genetics: ',x,
    '\n',
    file = pth, append = TRUE,sep = ''
)

(x <- 2463+359)
# 2822

cat("\n", 
    'Total Patients With Genetics: ',x,
    '\n',
    file = pth, append = TRUE,sep = ''
)



# EYES (((((((((((())))))))))))
eye_visits <- data.frame(colData(sce)) %>% 
  select(patID_eye, month, rs17421627_tmem161b) %>% 
  arrange(patID_eye, month) %>% 
  distinct()

# longitudinal patients
(x <- eye_visits %>% 
  filter(n()>1, .by = patID_eye) %>% 
  distinct(patID_eye) %>% 
  nrow())
# 975

cat("\n", 
    'Longitudinal Eyes: ',x,
    '\n',
    file = pth, append = TRUE,sep = ''
)


# single visit patients
(x <- eye_visits %>% 
  filter(n()==1, .by = patID_eye) %>% 
  distinct(patID_eye) %>% 
  nrow())
# 6913
y = 975+6913
z  = 975+6913==n_distinct(sce$patID_eye)
# 7888

cat("\n", 
    'Single Visit Eyes: ',x,
    '\n','\n',
    'Single plus Longitudinal Eyes: ',y,
    '\n','\n',
    'Total Eyes = : ',z,
    '\n',
    file = pth, append = TRUE,sep = ''
)



# longitudinal with genetics
(x <- eye_visits %>% 
  filter(!is.na(rs17421627_tmem161b)) %>% 
  filter(n()>1, .by = patID_eye) %>% 
  distinct(patID_eye) %>% 
  nrow())
# 712

cat("\n", 
    'Longitudinal Eyes With Genetics: ',y,
    '\n',
    file = pth, append = TRUE,sep = ''
)


# single visit genetics patients
(x <- eye_visits %>% 
  filter(!is.na(rs17421627_tmem161b)) %>% 
  filter(n()==1, .by = patID_eye) %>% 
  distinct(patID_eye) %>% 
  nrow())
# 4899

(y <- 4899+712)
# 5611

cat("\n", 
    'Single Visit Eyes With Genetics: ',x,
    '\n','\n',
    'Total Eyes With Genetics: ',y,
    '\n',
    file = pth, append = TRUE,sep = ''
)



# ((((((((()))))))))
# HOW MANY PEOPLE/EYES WITH/OUT GENETICS ----
(pat_w_genetics <- data.frame(colData(sce)) %>% 
  select(patID, matches("rs")) %>% 
  drop_na() %>% 
  distinct(patID) %>% 
  nrow())
# 2822


(pat_wout_genetics <- data.frame(colData(sce)) %>% 
  select(patID, matches("rs")) %>% 
  filter(if_all(matches("rs"), ~is.na(.))) %>% 
  distinct(patID) %>% 
  nrow())
# 1149

(r <- pat_wout_genetics/(pat_wout_genetics+pat_w_genetics))
# 0.289

cat("\n", 
    'Patients Without Genetics: ',pat_wout_genetics, sprintf("(%.2f%%)", r*100),
    '\n', 
    file = pth, append = TRUE,sep = ''
)


pat_wout_genetics+pat_w_genetics
# 3971
(pat_wout_genetics+pat_w_genetics)==n_distinct(sce$patID)

pat_w_genetics/(pat_wout_genetics+pat_w_genetics)
# 0.71

(eye_w_genetics <- data.frame(colData(sce)) %>% 
  select(patID_eye, matches("rs")) %>% 
  drop_na() %>% 
  distinct(patID_eye) %>% 
  nrow())
# 5611

(eye_wout_genetics <- data.frame(colData(sce)) %>% 
  select(patID_eye, matches("rs")) %>% 
  filter(if_all(matches("rs"), ~is.na(.))) %>% 
  distinct(patID_eye) %>% 
  nrow())
# 2277

(r <- eye_wout_genetics/(eye_wout_genetics+eye_w_genetics))
# 0.289

cat("\n", 
    'Eyes Without Genetics: ',eye_wout_genetics, sprintf("(%.2f%%)", r*100),
    '\n', 
    file = pth, append = TRUE,sep = ''
)


(eye_wout_genetics+eye_w_genetics)
# 7888
(eye_wout_genetics+eye_w_genetics)==n_distinct(sce$patID_eye)
eye_w_genetics/(eye_wout_genetics+eye_w_genetics)
# 0.71

# SEX ----
(x = data.frame(colData(sce)) %>% 
  select(patID, SEXF) %>% 
  group_by(patID) %>% 
  mutate(SEXF =  
           ifelse(is.na(SEXF), dplyr::first(SEXF[!is.na(SEXF)]) , SEXF) ) %>% 
  ungroup() %>% 
  distinct() %>% 
  pull(SEXF) %>% 
  table(useNA = "always") %>% 
  {rbind(count = .,prop = ./sum(.)) %>% round(3) %>% {cbind(., column_total = rowSums(.))}} %>% 
  t() %>% 
   data.frame() %>% 
   rownames_to_column('x'))

cat("\n**SEXF**\n", file = pth, append = TRUE,sep = '')
write_tsv(x, file = pth , append = TRUE)


# TABLE OF DEMOGRAPHICS----
# pat_dat <- colData(sce) %>% 
#   data.frame() %>% 
#   select(patID,matches("SEXF|t2d|SMK|AGE|RACE|ETH|AGE|chew.*loose")) %>% 
#   filter(row_number()==1,.by = patID)
#   
# 
# tabs <- pat_dat %>% 
#   sapply(function(x){
#     if(n_distinct(x)<30){
#       table(x, useNA = 'ifany')
#     }else{
#       summary(x)
#     }
#   })
# 
# tabs <- tabs[order(names(tabs))]
# 
# 
# x <- pat_dat %>% 
#   mutate(t2d = factor(t2d)) %>% 
#   mutate(chew_grade_derived_mine_loose_fac = factor(chew_grade_derived_mine_loose)) %>% 
#   mutate(chew = 
#            case_when(
#     chew_grade_derived_mine_loose <2 ~"alsthn2", 
#     chew_grade_derived_mine_loose<4~"blsthn4", 
#     chew_grade_derived_mine_loose>=4 ~ 'cgrthn4',
#     TRUE~as.character(chew_grade_derived_mine_loose)
#   ) %>% factor()) %>% 
#   mutate(race = case_when(
#    RACE == 1  ~'American Indian or Alaskan Native',
#    RACE == 2  ~'Asian',
#    RACE == 3  ~'Black/African',
#    RACE == 4  ~'Native Hawaiian or Pacific Islander',
#    RACE == 5  ~'White/Caucasian',
#    RACE == 6  ~'Australian Aboriginal or Torres Strait Islander',
#    RACE == 9  ~'Unable to specify',
#    RACE == 99 ~'Other', 
#    TRUE ~ as.character(RACE)
#   ))
# 
# x %>% datasummary_skim()
# 
# x %>%
#   select(SEXF, t2d, AGE, chew_grade_derived_mine_loose_fac, SMKCURRENT) %>%
#   gtsummary::tbl_summary()
# 
# 
# smkrs <- pat_dat$SMKCURRENT %>% table(useNA = "always")
# 
# Psmk <- smkrs[2]/( smkrs[2]+ smkrs[1])
# 
# Psmk_est <- (Psmk*sum(smkrs))/sum(smkrs)
# 
# txtbarchart(x$chew_grade_derived_mine_loose %>% factor())
# 
# dat <- assays(sce)$pmm_imp_of_scaled %>% 
#   t() %>% 
#   data.frame() %>% 
#   rownames_to_column("sample") %>% 
#   left_join(colData(sce) %>% data.frame() %>% select(sample, my_cluster))
# 
# dat_long <- dat %>% 
#   slice_sample(n = 3e3) %>%
#   pivot_longer( -c(sample, my_cluster))
# 
# 


# CHEW FREQUENCIES ----
# paper_chew_freq = tribble(
#   ~stage, ~count,
#   0, 1759,
#   1, 2134,
#   2, 1395,
#   3, 777,
#   4, 3946,
#   5, 1498,
#   6, 767,
# )
# 
# z1 = sum(paper_chew_freq$count[paper_chew_freq$stage %in% 0:2])
# z2 = sum(paper_chew_freq$count[paper_chew_freq$stage %in% 3:4])
# z3 = sum(paper_chew_freq$count[paper_chew_freq$stage %in% 5:6])
# 
# cat('\n', 
#   sprintf("%s from 0 to 2 (%.2f%% of total)", z1, 100*z1/sum(paper_chew_freq$count)),
#   '\n',
#   sprintf("%s from 3 to 4 (%.2f%% of total)", z2, 100*z2/sum(paper_chew_freq$count)),
#   '\n',
#   sprintf("%s from 5 to 6 (%.2f%% of total)", z3, 100*z3/sum(paper_chew_freq$count)),
#   '\n'
# )



chew_freq = data.frame(colData(sce)) %>% 
  select(patID_eye, month, chew_grade_derived_mine_loose) %>% 
  group_by(patID_eye) %>% 
  mutate(across(chew_grade_derived_mine_loose, 
                ~ifelse(is.na(.), dplyr::first(.[!is.na(.)]) , .) )) %>% 
  arrange(month) %>% 
  slice_head(n = 1) %>% 
  pull(chew_grade_derived_mine_loose) %>% 
  table() %>% 
  as.data.frame() %>% 
  rename_with(~c("stage", "count"))

k1 = sum(chew_freq$count[chew_freq$stage %in% 0:2])
k2 = sum(chew_freq$count[chew_freq$stage %in% 3:4])
k3 = sum(chew_freq$count[chew_freq$stage %in% 5:6])

p1 = k1/sum(chew_freq$count)*100
p2 = k2/sum(chew_freq$count)*100
p3 = k3/sum(chew_freq$count)*100



cat(
  '\n', 
  '**Chew Summary**\n',
  sprintf("%s from 0 to 2 (%.2f%% of total)", k1, p1),
  '\n',
  sprintf("%s from 3 to 4 (%.2f%% of total)", k2, p2),
  '\n',
  sprintf("%s from 5 to 6 (%.2f%% of total)", k3, p3),
  '\n', 
  '\nstage\tcount\n', 
  file = pth,   append = TRUE, sep = ''
)
write_tsv(chew_freq, file =pth, append = TRUE)



# >>> 4437 from 0 to 2 (56.25% of total) 
# >>> 2490 from 3 to 4 (31.57% of total) 
# >>> 961 from 5 to 6 (12.18% of total) 

# AGE SUMMARY ----
(x = colData(sce) %>% 
  data.frame() %>% 
  select(patID, month, AGE) %>% 
  group_by(patID) %>% 
  arrange(month) %>% 
  slice_head(n = 1) %>% 
  pull(AGE) %>% 
  summary())

cat("\n", 
    '**Age Summary**',
    '\n',
    'min Q1 med mean Q3 max NAs',
    "\n",
    paste(round(x,2), collapse = '\t'),
    '\n', 
    file = pth, append = TRUE,sep = ''
)



# T2D ----

z <- colData(sce) %>% 
   data.frame() %>% 
   select(patID,patID_eye, t2d) %>% 
   group_by(patID) %>% 
   mutate(t2d = ifelse(is.na(t2d), t2d[!is.na(t2d)] , t2d) )

# # do patients have unique t2d status
# sapply(unique(x$patID), function(pat){
#   n_distinct(x$t2d[x$patID == pat])
# }) %>% 
#   table()

(pat_t2d <- z %>% 
  group_by(patID) %>% 
  slice_head(n = 1) %>% 
  pull(t2d) %>% 
  table(useNA = "always")
  )

(eye_t2d <- z %>% 
  group_by(patID_eye) %>% 
  slice_head(n = 1) %>% 
  pull(t2d) %>% 
  table(useNA = "always")
  )


cat("\n", 
    '**T2D Summary**',
    '\n',
    'Patients without T2D: ', pat_t2d[1] ,
    "\n", 
    'Patients with T2D: ', pat_t2d[2] ,
    "\n", 
    'Patients with missing T2D: ', pat_t2d[3] ,
    "\n",
    "\n",
    'Eyes without T2D: ',      eye_t2d[1] ,
    "\n", 
    'Eyes with T2D: ',         eye_t2d[2] ,
    "\n", 
    'Eyes with missing T2D: ', eye_t2d[3] ,
    "\n", 
    '\n', 
    file = pth, append = TRUE,sep = ''
)
