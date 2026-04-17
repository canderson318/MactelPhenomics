

pacman::p_load(tidyr, dplyr, reshape2, purrr, tidyverse, gridExtra, ComplexHeatmap, viridis, zoo, mgcv, ggrepel, car, caret, stats, ape, SingleCellExperiment, scater, scran, xgboost, mice, parallel, DiagrammeR, patchwork)

rm(list = ls())
gc()

source("/Users/canderson/Documents/LMRI/PRJ2024001/version015/analysis/00_MyFunctions.R")
how.long()


#(((((((((((((((())))))))))))))))
# Load Data ----

wd_path = "~/Documents/LMRI/PRJ2024001/version015/"
# load experiment
sce <-  HDF5Array::loadHDF5SummarizedExperiment(sprintf("%sprocessed_data/sce002.h5", wd_path))


# load field info sheet
fields <- readxl::read_xlsx("/Users/christiananderson/Library/CloudStorage/Egnyte-lmri/Shared/LMRI/DS_Group/PRJ2024001/info/Grading fields and variables_edited.xlsx", 
                            sheet = 1)
# load classification info
chew_class <- readxl::read_xlsx("/Users/christiananderson/Library/CloudStorage/Egnyte-lmri/Shared/LMRI/DS_Group/PRJ2024001/info/Grading fields and variables_edited.xlsx", 
                            sheet = 2)

# load simple classification info
chew_class_simp <- readxl::read_xlsx("/Users/christiananderson/Library/CloudStorage/Egnyte-lmri/Shared/LMRI/DS_Group/PRJ2024001/info/Grading fields and variables_edited.xlsx", 
                            sheet = 3)


a_path <- paste0(wd_path, "results/clustering/09_clustering/9.2.2_deriveChewGrade/")
if(!dir.exists(a_path)) {
  dir.create(a_path, recursive = TRUE)
}

#(((((((((((((((())))))))))))))))
# CLEAN KEYS ----

# ****** FUNCTION: clean chew grade sheets ******
clean_chew <- function(chew_class_df){
  
  # make field field number
  chew_class_df$Field_number <- chew_class_df$Field
  
  # split grade into info and grade
  chew_class_df$Grade_info <- sapply(str_split(chew_class_df$Grade, " = "), `[`, 2)
  chew_class_df$Grade <- as.numeric(sapply(str_split(chew_class_df$Grade, " = "), `[`, 1))
  
  # reorder cols
  v <- c("Grade", "Field_number", "Valid_inputs")
  chew_class_df <- chew_class_df[, c(v, setdiff(names(chew_class_df), v))]
  
  # replicate row however many Field times, e.g. FIeld = '58-66' x9 rows
  new_chew_class_df <- chew_class_df[0,]
  # i =2 
  # j =58
  for(i in 1:nrow(chew_class_df) ){
    a_row <- chew_class_df[i, ]
    a_row$Field_number <- as.character(a_row$Field_number)
    a_field <- a_row["Field_number"]
    
    # if range of fields
    if(grepl("-", a_field)){
      start <- as.numeric(str_split(a_field, "-")[[1]][1])
      end <- as.numeric(str_split(a_field, "-")[[1]][2])
      range <- start:end
      new_a_row <- a_row
      for(j in range){
        new_a_row$Field_number <- as.character(j)
        new_chew_class_df <- new_chew_class_df %>% 
          rbind(new_a_row)
      }
    }else{
      new_chew_class_df <- new_chew_class_df %>% 
        rbind(a_row)
    }
  }
  
  # match field to field_number
  new_chew_class_df$Field <-  fields$Variable_name[match(new_chew_class_df$Field_number, fields$Field_number)]
  new_chew_class_df <- as.data.frame(new_chew_class_df)
  
  return(new_chew_class_df)
}


# shape the grading sheets with function
new_chew_class_orig <- clean_chew(chew_class)
new_chew_class_simp <- clean_chew(chew_class_simp)

# grab what vars used in chew grade
# chew_vars <- new_chew_class_orig$Field %>% unique

chew_vars <- c("ISOS_BREAK", "IR01", "IR02", "IR03", "IR04", "IR05", "IR06", "IR07", "IR08", "IR09", "HR_TEMPORAL", "ISOS_LOC_CEN","ISOS_LOC_JF", "NEOVASC")

# make new one to edit
new_chew_class <- new_chew_class_orig %>% 
  select(-matches("Classification|Grade_info|Field_number"))

# plot what vars for each grade
custom_colors <- c(
  rgb(0, .5, 1),       
  rgb(0, 1, 0),   
  rgb(1, 1, 0),
  rgb(1,.5,0),
  rgb(1, 0, 0)        
)

foo <- function(){
  custom_colors <- c(
    rgb(0, .5, 1),       
    rgb(0, 1, 0),   
    rgb(1, 1, 0),
    rgb(1,.5,0),
    rgb(1, 0, 0)        
  )
      new_chew_class %>%
      # new_chew_class_orig %>%
      # new_chew_class_simp %>%
      mutate(Field = ifelse(is.na(Field), "NA", Field)) %>% 
      select(Grade, Field, Valid_inputs) %>% 
      # separate_rows(Valid_inputs, sep = ",") %>% 
      mutate(Valid_inputs = 
               factor(Valid_inputs, levels = c("0", "1", "2", "1,2", "0,1,2")) 
      ) %>%
      unique() %>% 
      ggplot(aes(factor(Grade), Field, fill = Valid_inputs, label = Field))+
      # geom_line(position = position_jitter(width = .01, height = .05))+
      # geom_tile()+
      scale_fill_manual(values = custom_colors)+
      geom_tile(alpha = .5)+
      geom_text(fontface = "bold")+
      scale_x_discrete(breaks = factor(0:6))+
      theme_bw()+
      theme(axis.text.y = element_blank(), 
            axis.ticks.y =  element_blank(),
            panel.grid = element_blank())+
      ggtitle("Chew Grade Criteria")
    # ggtitle("Chew Grade Simple Criteria")
  }


setdiff(new_chew_class$Field, rowData(sce)$name)
# ^ IR06-IR07 missing because they were removed for being too lossy; ISOS_LOC 1 hot encoded


# # make explicit grading, i.e. where no value for NEOVASC, make 0
# ## Back log 0 for neovasc
# new_neovascs <- new_chew_class[new_chew_class$Field == "NEOVASC", ][rep(1,5),]
# new_neovascs$Grade <- 0:4
# new_chew_class <- new_chew_class %>% rbind(new_neovascs)
# 
# ## fill in grade 3 IR01
# new_IR01 <- new_chew_class[new_chew_class$Grade==3 & new_chew_class$Field== "IR02",]
# new_IR01$Valid_inputs <- 0
# new_IR01$Field <- "IR01"
# new_chew_class <- rbind(new_chew_class, new_IR01)
# 
# ## fill in grade 0 ISOS_LOC
# new_ISOS_LOC <- new_chew_class[new_chew_class$Grade==1 & new_chew_class$Field== "ISOS_LOC",][c(1,1),]
# new_ISOS_LOC$Valid_inputs <- c("0", "0,1,2")
# new_ISOS_LOC$Grade <- c(0, 6)
# new_ISOS_LOC$Field <- "ISOS_LOC"
# new_chew_class <- rbind(new_chew_class, new_ISOS_LOC)
# 
# ## foward fill HR_TEMPORAL
# new_HR_TEMPORAL <- new_chew_class[new_chew_class$Grade==1 & new_chew_class$Field== "HR_TEMPORAL",][rep(1,2),]
# new_HR_TEMPORAL$Grade <- 5:6
# new_HR_TEMPORAL$Valid_inputs <- rep(c("0,1,2"), 2)
# new_chew_class <- rbind(new_chew_class, new_HR_TEMPORAL)
# 
# ## foward fill ISOS_BREAK
# new_ISOS_BREAK <- new_chew_class[new_chew_class$Grade==1 & new_chew_class$Field== "ISOS_BREAK",]
# new_ISOS_BREAK$Grade <- 6
# new_ISOS_BREAK$Valid_inputs <- c("0,1,2")
# new_chew_class <- rbind(new_chew_class, new_ISOS_BREAK)
# 
# Fix IR02-9
# 
## Treat each IR02-9 individually
# ### Forward fill IR02:9
# IR_fields <- grep("IR0[23456789]", new_chew_class$Field, value = TRUE) %>% unique()
# new_IR <- expand.grid(Field = IR_fields, Grade = 5:6)
# new_IR$Valid_inputs <- "0,1,2"
# new_chew_class <- new_chew_class %>% rbind(new_IR)

## Treat IR02-9 as one variable
# ^ the valid input for IR02-9 is the same for each grade
new_chew_class <- new_chew_class %>%
  mutate(Field = Field %>% gsub("IR0[23456789]","IR02-9", .)) %>%
  distinct()

# ### foward fill 'IR02-9' 
# IR_fields <- 'IR02-9'
# new_IR <- expand.grid(Field = IR_fields, Grade = 5:6)
# new_IR$Valid_inputs <- "0,1,2"
# new_chew_class <- new_chew_class %>% rbind(new_IR)

foo()




# (((((((((((((())))))))))))))
# MAKE TRAIN AND TEST DFS----
# (((((((((((((())))))))))))))

# assay of chew vars
GT <-
  assays(sce)$measure_vars %>% 
  data.matrix() %>% 
  t() %>% 
  cbind(chew_grade = sce$chew_grade) %>% 
  as.data.frame() %>% 
  select(any_of( sort(chew_vars , T)), chew_grade) %>% 
  drop_na()  %>% 
  # # fix ISOS_LOC to be 0 not present, 1 juxtafoveal, 2 central
  # mutate(ISOS_LOC = case_when(
  #   ISOS_LOC == 1 ~ 2,
  #   ISOS_LOC == 2 ~ 1,
  #   ISOS_LOC == 0 ~ 0
  # )) %>%
  mutate(chew_grade = as.numeric(chew_grade))

set.seed(12030)
prop_test <- .3
test_inds <- sample(seq(nrow(GT)), prop_test*nrow(GT))
test_inds_samples <- rownames(GT)[test_inds]

# imputation train/test dataset
imp_train_X <- GT
imp_train_X$chew_grade[test_inds] <- NA
test_y <- GT$chew_grade[test_inds]

# XGB train; hold test indices for later
XGB_train_X <- GT[-test_inds, -which(names(GT)== "chew_grade")]
XGB_train_y <- GT$chew_grade[-test_inds]

# XGB test
XGB_test_X <- GT[test_inds, -which(names(GT)== "chew_grade")]
XGB_test_y <- test_y




# (((((((((((((())))))))))))))
# DERIVE CHEW GRADE ----
# (((((((((((((())))))))))))))

#+ 
#+ when no value defined, allow any value?
#+ 

## use ground truth and simplify the vars to ma
derive_assay <-
  assays(sce)$measure_vars %>% 
  data.matrix() %>% 
  t() %>% 
  cbind(chew_grade = sce$chew_grade) %>% 
  as.data.frame() %>% 
  select(any_of( sort(chew_vars , T)), chew_grade) %>% 
  # filter for where all vars not missing allowing  for chew grade to be missing
  drop_na(-chew_grade) %>% 
  data.matrix() %>% 
  # simplify IR02-9 to one var of the max 
  { 
    mutate(as.data.frame(.), 
           `IR02-9` = rowMaxs(as.matrix(select(as.data.frame(.), matches("IR0[23456789]$"))))
    )
  } 

## assume location Juxtafoveal when isos_break observed
derive_assay %>% 
  filter(ISOS_BREAK>0 & ISOS_LOC_JF ==0 & ISOS_LOC_CEN == 0)

# where ISOS_Break and no location, make jf 1; make cen 0 where JF
derive_assay <- derive_assay %>% 
  mutate(ISOS_LOC_JF = ifelse(ISOS_BREAK>0 & ISOS_LOC_JF ==0 & ISOS_LOC_CEN == 0, 1, ISOS_LOC_JF)) %>% 
  mutate(ISOS_LOC_CEN = ifelse(ISOS_LOC_JF >0, 0, ISOS_LOC_CEN )) 

# (((((((((((((())))))))))))))----
# MANUAL APPROACHES----
foo()

# *****
## MINE strict ----
#+ from chew et.al. use their scheme as a progressive linear scale.
#+ Tunde's scheme misinterpreted the '3-- Noncentral pigment/No, non-central, or central EZ/No OCT HR'
#+
#+ 0 - everyone
#+ 1 - ez break =2 + ez location jf =1
#+ 2 - ez break =2 + ez location central =1
#+ 3 - IR02-09
#+ 4 - HR 4
#+ 5 - IR01
#+ 6 - Neov

derive_assay$chew_grade_derived_mine_strict <- NA

# Grade 0
derive_assay$chew_grade_derived_mine_strict[
  derive_assay$ISOS_BREAK == 0 &
    derive_assay$`IR02-9` == 0 &
    derive_assay$IR01 == 0 &
    derive_assay$HR_TEMPORAL == 0
] <- 0


# GRADE 1
derive_assay$chew_grade_derived_mine_strict[
  derive_assay$ISOS_LOC_JF==1 & # 5 1 = Noncentral EZ break/ No pigment/ No OCT HR Location of the IS/OS PR break = valid input 2 (off centre) >> Original ISOS_LOC = 2
    derive_assay$ISOS_BREAK==2
] <- 1

# GRADE 2
derive_assay$chew_grade_derived_mine_strict[
  derive_assay$ISOS_LOC_CEN==1 &  #  2 = Central EZ break/No pigment/ No OCTHR Location of the IS/OS PR break = valid input 1 (central) >> Original ISOS_LOC = 1
    derive_assay$ISOS_BREAK==2
] <- 2


# GRADE 3
derive_assay$chew_grade_derived_mine_strict[
  derive_assay$`IR02-9`==2
] <- 3

# GRADE 4
derive_assay$chew_grade_derived_mine_strict[
  derive_assay$HR_TEMPORAL ==2
] <- 4

# GRADE 5
derive_assay$chew_grade_derived_mine_strict[
  derive_assay$IR01 == 2
] <- 5

# GRADE 6
derive_assay$chew_grade_derived_mine_strict[
  derive_assay$NEOVASC == 2
] <- 6


# grade for unclassifiable obs
derive_assay$chew_grade_derived_mine_strict[with(derive_assay,
  (ISOS_BREAK > 0 |
     `IR02-9` >0 |
     NEOVASC >0 |
     HR_TEMPORAL > 0) &
    is.na(chew_grade_derived_mine_strict)
)]<- -Inf


# *****
## MINE loose ----
#+ from chew et.al. use their scheme as a progressive linear scale.
#+ Tunde's scheme misinterpreted the '3-- Noncentral pigment/No, non-central, or central EZ/No OCT HR'
#+
#+ 0 - everyone
#+ 1 - ez break =2 + ez location jf =1
#+ 2 - ez break =2 + ez location central =1
#+ 3 - IR02-09
#+ 4 - HR 4
#+ 5 - IR01
#+ 6 - Neov

derive_assay$chew_grade_derived_mine_loose <- NA

# Grade 0
derive_assay$chew_grade_derived_mine_loose[
  derive_assay$ISOS_BREAK == 0 &
    derive_assay$`IR02-9` == 0 &
    derive_assay$IR01 == 0 &
    derive_assay$HR_TEMPORAL == 0
] <- 0


# GRADE 1
derive_assay$chew_grade_derived_mine_loose[
  derive_assay$ISOS_LOC_JF==1 &
    derive_assay$ISOS_BREAK%in%c(1,2)
] <- 1

# GRADE 2
derive_assay$chew_grade_derived_mine_loose[
  derive_assay$ISOS_LOC_CEN==1 &
    derive_assay$ISOS_BREAK%in%c(1,2)
] <- 2


# GRADE 3
derive_assay$chew_grade_derived_mine_loose[
  derive_assay$`IR02-9` %in%c(1,2)
] <- 3

# GRADE 4
derive_assay$chew_grade_derived_mine_loose[
  derive_assay$HR_TEMPORAL %in%c(1,2)
] <- 4

# GRADE 5
derive_assay$chew_grade_derived_mine_loose[
  derive_assay$IR01 %in%c(1,2)
] <- 5

# GRADE 6
derive_assay$chew_grade_derived_mine_loose[
  derive_assay$NEOVASC %in%c(1,2)
] <- 6


# grade for unclassifiable obs
derive_assay$chew_grade_derived_mine_loose[with(derive_assay,
  (ISOS_BREAK > 0 |
     `IR02-9` > 0 |
     NEOVASC >0 |
     HR_TEMPORAL > 0) &
    is.na(chew_grade_derived_mine_loose)
)]<- -Inf

derive_assay[,grep("chew", names(derive_assay), value = TRUE)] %>% 
  purrr::map(table, useNA = 'always')

#****
## Yue's ----
#+ function converted to R
#+ ''implementation of tree in paper''
#+

# make new chew var
derive_assay$chew_grade_derived_yue <- NA

derive_assay$chew_grade_derived_yue[with(derive_assay,
  ISOS_BREAK == 0 &
  `IR02-9` ==0 &
  NEOVASC ==0 &
  HR_TEMPORAL == 0
)]<- 0

# derive_assay$chew_grade_derived_yue[with(derive_assay,
#  ISOS_BREAK == 0 &
#  `IR02-9` ==1 &
#  NEOVASC ==0 &
#  HR_TEMPORAL == 0
#   )]<- 0.5

derive_assay$chew_grade_derived_yue[with(derive_assay,
 ISOS_BREAK == 1 &
 `IR02-9`%in% c(0,1) &
 NEOVASC ==0 &
 HR_TEMPORAL == 0
  )]<- 1

derive_assay$chew_grade_derived_yue[with(derive_assay,
 `IR02-9`==0 &
 ISOS_BREAK == 2 &
 NEOVASC ==0 &
 HR_TEMPORAL == 0
  )]<- 2

derive_assay$chew_grade_derived_yue[with(derive_assay,
  `IR02-9` == 1 &
  ISOS_BREAK ==2 &
  HR_TEMPORAL ==0 &
  NEOVASC ==0
  )] <- 3

# derive_assay$chew_grade_derived_yue[with(derive_assay,
#   `IR02-9` == 0 &
#   ISOS_BREAK ==0 &
#   HR_TEMPORAL ==1 &
#   NEOVASC ==0
#   )] <- 3.5


derive_assay$chew_grade_derived_yue[with(derive_assay,
  `IR02-9`==0 &
  ISOS_BREAK %in%c(1,2) &
  NEOVASC ==0 &
  HR_TEMPORAL %in% c(1,2)
  )] <- 4


# derive_assay$chew_grade_derived_yue[with(derive_assay,
#   `IR02-9`== 1 &
#   ISOS_BREAK %in%c(1,2) &
#   NEOVASC ==0 &
#   HR_TEMPORAL %in% c(1,2)
#   )] <- 4.5

derive_assay$chew_grade_derived_yue[with(derive_assay,
  `IR02-9` == 2 &
  NEOVASC == 0
  )]<- 5


derive_assay$chew_grade_derived_yue[with(derive_assay,
  `IR02-9` == 2 &
  NEOVASC == 2
  )]<- 6



# grade for unclassifiable obs
derive_assay$chew_grade_derived_yue[with(derive_assay,
 (ISOS_BREAK > 0 |
    `IR02-9` >0 |
    NEOVASC >0 |
    HR_TEMPORAL > 0) &
  is.na(chew_grade_derived_yue)
  )]<- -Inf

table(derive_assay$chew_grade_derived_yue, useNA = 'ifany')



# derive_assay$chew_grade_derived_yue <- 0
# derive_assay$chew_grade_derived_yue[with(derive_assay,
#   `IR02-9` == 2 & 
#   NEOVASC == 2 
#   )]<- 6
# 
# derive_assay$chew_grade_derived_yue[with(derive_assay,
#   `IR02-9` == 2 & 
#   NEOVASC == 0 
#   )]<- 5
# 
# derive_assay$chew_grade_derived_yue[with(derive_assay,
#   `IR02-9` == 1 & 
#   ISOS_LOC_CEN == 1
#   )] <- 3
# 
# derive_assay$chew_grade_derived_yue[with(derive_assay,
#   `IR02-9`==0 & 
#   ISOS_BREAK %in%c(1,2) & 
#   NEOVASC == 2
#   )] <- 6
# 
# derive_assay$chew_grade_derived_yue[with(derive_assay,
#   `IR02-9`==0 & 
#   ISOS_BREAK %in%c(1,2) & 
#   NEOVASC ==0 & 
#   HR_TEMPORAL %in% c(1,2)
#   )] <- 4
# 
# derive_assay$chew_grade_derived_yue[with(derive_assay,
#   `IR02-9`==0 &
#   ISOS_BREAK == 2 &
#   NEOVASC ==0 &
#   HR_TEMPORAL == 0
#   )]<- 2
# derive_assay$chew_grade_derived_yue[with(derive_assay,
#   `IR02-9`==0 &
#   ISOS_BREAK == 1 &
#   NEOVASC ==0 &
#   HR_TEMPORAL == 0
#   )]<- 1
# 
# # flag unclassified
# derive_assay$chew_grade_derived_yue[with(derive_assay,
#                                          `IR02-9`>0 | 
#                                            ISOS_BREAK >0 |
#                                            NEOVASC >0 | 
#                                            HR_TEMPORAL > 0 
# )]<- -Inf
# 
# table(derive_assay$chew_grade_derived_yue, useNA = "always")


### plot all chew grades ----
tsnes <- sapply(grep("chew", names(derive_assay), value = TRUE), function(nm){
  if(nm == "chew_grade"){
    plotReducedDim(sce[,!is.na(sce$chew_grade)], "TSNE", color_by = "chew_grade")+
      theme(aspect.ratio = 1)+
      scale_color_viridis_c(option = 'plasma')+
      # scale_color_gradientn(colours =  rev(palScripps(10)))+
      ggtitle(nm)
    
  }else{
    samps <- rownames(derive_assay)[!is.na(derive_assay[[nm]]) & !is.infinite(derive_assay[[nm]])]
    plotReducedDim(sce[ ,  colnames(sce)%in% samps], "TSNE", color_by = I(derive_assay[rownames(derive_assay)%in% samps, nm]))+
      theme(aspect.ratio = 1)+
      scale_color_viridis_c(option = 'plasma')+
      # scale_color_gradientn(colours = rev(palScripps(10)))+
      ggtitle(nm)
  }
}, simplify = FALSE)

all_tsnes <- patchwork::wrap_plots(tsnes, ncol = 2)

list(plotTSNE(sce, color_by = "ACUIT")+scale_color_viridis_c(option = "plasma", begin = 1, end = 0)+theme(aspect.ratio = 1) + labs(color = "ACUIT"),
    plotTSNE(sce, color_by = "pseudoprog")+      scale_color_viridis_c(option = "plasma", begin = 0, end = 1)+theme(aspect.ratio = 1) + labs(color = "pseudoprog")) %>% wrap_plots(ncol = 1)


### Disect missing classifications ----
hms <- sapply(grep("chew", names(derive_assay), value = TRUE), function(nm){
  derive_assay %>% 
    filter(if_all(all_of(nm), ~{.<0 | is.na(.)} )) %>%
    # filter(!is.na(chew_grade)) %>%
    select( NEOVASC ,matches("ISOS"), HR_TEMPORAL, `IR02-9`) %>%
    arrange(NEOVASC) %>%
    distinct() %>% 
    data.matrix() %>% 
    t() %>% 
    Heatmap(col = plasma(5), column_title = nm)
}, simplify = FALSE)   


hms_grid <- gridExtra::marrangeGrob(lapply(hms, function(x) grid.grabExpr(draw(x))), ncol = 3, nrow = 2)
ggsave(sprintf("%sresults/clustering/09_clustering/9.2.2_deriveChewGrade/distinct_non_classified_rows.png", wd_path), hms_grid, height = 10, width = 30)


derive_assay %>% 
  select(matches("chew")) %>% 
  mutate_all(~ifelse(is.infinite(.), NA, .)) %>% 
  cor(., use = "pairwise.complete") %>% 
  Heatmap(name = "r", column_title = "Derived Chew Grade Correlations", 
          col = plasma(10, end = .9)
          )


# (((((((((((((())))))))))))))
# Create the flowcharts ----
flowchart_mine_loose <- {grViz("
digraph flowchart {
  graph [layout = dot, rankdir = LR,  label = 'My loose chew grade', labelloc = t, fontsize = 20]  // Horizontal layout (Left-to-Right)

  # Nodes for grades
  node [shape = rectangle, style = filled, color = lightblue]
  grade [label = 'Grade NA']
  grade0 [label = 'Grade 0\nISOS_BREAK == 0\nIR02-9 == 0\nIR01 == 0\nHR_TEMPORAL == 0']
  grade1 [label = 'Grade 1\nISOS_LOC == JF\nISOS_BREAK in {1,2}']
  grade2 [label = 'Grade 2\nISOS_LOC == CEN\nISOS_BREAK in {1,2}']
  grade3 [label = 'Grade 3\nIR02-9 in {1,2}']
  grade4 [label = 'Grade 4\nHR_TEMPORAL in {1,2}']
  grade5 [label = 'Grade 5\nIR01 in {1,2}']
  grade6 [label = 'Grade 6\nNEOVASC in {1,2}']
  gradeInf [label = 'All else -> -Inf\nis.na(Grade)\nISOS_BREAK > 0\nIR02-9 > 0\nIR01 > 0\nHR_TEMPORAL > 0']

  # Edges between nodes
  //start [label = 'Start', shape = ellipse, color = gray]
  grade -> grade0
  grade0 -> grade1
  grade1 -> grade2
  grade2 -> grade3
  grade3 -> grade4
  grade4 -> grade5
  grade5 -> grade6
  grade6 -> gradeInf
}
")}
flowchart_mine_strict <- {grViz("
digraph flowchart {
  graph [layout = dot, rankdir = LR,  label = 'My loose chew grade', labelloc = t, fontsize = 20]  // Horizontal layout (Left-to-Right)

  # Nodes for grades
  node [shape = rectangle, style = filled, color = lightblue]
  grade [label = 'Grade NA']
  grade0 [label = 'Grade 0\nISOS_BREAK == 0\nIR02-9 == 0\nIR01 == 0\nHR_TEMPORAL == 0']
  grade1 [label = 'Grade 1\nISOS_LOC == JF\nISOS_BREAK == 2}']
  grade2 [label = 'Grade 2\nISOS_LOC == CEN\nISOS_BREAK == 2']
  grade3 [label = 'Grade 3\nIR02-9 == 2']
  grade4 [label = 'Grade 4\nHR_TEMPORAL == 2']
  grade5 [label = 'Grade 5\nIR01 == 2']
  grade6 [label = 'Grade 6\nNEOVASC == 2']
  gradeInf [label = 'All else -> -Inf\nis.na(Grade)\nISOS_BREAK > 0\nIR02-9 > 0\nIR01 > 0\nHR_TEMPORAL > 0']

  # Edges between nodes
  //start [label = 'Start', shape = ellipse, color = gray]
  grade -> grade0
  grade0 -> grade1
  grade1 -> grade2
  grade2 -> grade3
  grade3 -> grade4
  grade4 -> grade5
  grade5 -> grade6
  grade6 -> gradeInf
}
")}

flowchart_mine_loose
flowchart_mine_strict


# Evaluate derived chew grade 
## DERIVED chew stats
calc_accu <- function(m){ sum(diag(m))/sum(m)}

# this assay has more data than GT so grab indices that match teh test indices from GT
derive_assay_test_inds <- match( test_inds_samples, rownames(derive_assay))


# performance metrics
(chew_grade_derived_mine_strict_stats <- list(
  CONFUSION = (simple_derived_confusion <- table(predicted = derive_assay$chew_grade_derived_mine_strict[derive_assay_test_inds], 
                                          actual = derive_assay$chew_grade[derive_assay_test_inds]) 
  ),
  N = sum(simple_derived_confusion),
  ACCU = calc_accu(simple_derived_confusion),
  RMSE = RMSE(derive_assay$chew_grade_derived_mine_strict[derive_assay_test_inds ], 
              derive_assay$chew_grade[derive_assay_test_inds], na.rm = TRUE), 
  TEST_COR = cor(derive_assay$chew_grade[derive_assay_test_inds], derive_assay$chew_grade_derived_mine_strict[derive_assay_test_inds], use = "complete.obs"),
  TRAIN_COR = cor(derive_assay$chew_grade[-derive_assay_test_inds], derive_assay$chew_grade_derived_mine_strict[-derive_assay_test_inds], use = "complete.obs"),
  OVERALL_COR = cor(derive_assay$chew_grade, derive_assay$chew_grade_derived_mine_strict, use = "complete.obs")
))

(chew_grade_derived_mine_loose_stats <- list(
  CONFUSION = (strict_derived_confusion <- table(predicted = derive_assay$chew_grade_derived_mine_loose[derive_assay_test_inds], 
                                          actual = derive_assay$chew_grade[derive_assay_test_inds]) 
  ),
  N = sum(strict_derived_confusion),
  ACCU = calc_accu(strict_derived_confusion),
  RMSE = RMSE(derive_assay$chew_grade_derived_mine_loose[derive_assay_test_inds ], 
              derive_assay$chew_grade[derive_assay_test_inds], na.rm = TRUE), 
  TEST_COR = cor(derive_assay$chew_grade[derive_assay_test_inds], derive_assay$chew_grade_derived_mine_loose[derive_assay_test_inds], use = "complete.obs"),
  TRAIN_COR = cor(derive_assay$chew_grade[-derive_assay_test_inds], derive_assay$chew_grade_derived_mine_loose[-derive_assay_test_inds], use = "complete.obs"),
  OVERALL_COR = cor(derive_assay$chew_grade, derive_assay$chew_grade_derived_mine_loose, use = "complete.obs")
))


(chew_grade_derived_yue_stats <- list(
  CONFUSION = (loose_derived_confusion <- table(predicted = derive_assay$chew_grade_derived_yue[derive_assay_test_inds],
                                          actual = derive_assay$chew_grade[derive_assay_test_inds])
  ),
  N = sum(loose_derived_confusion),
  ACCU = calc_accu(loose_derived_confusion),
  RMSE = RMSE(derive_assay$chew_grade_derived_yue[derive_assay_test_inds ],
              derive_assay$chew_grade[derive_assay_test_inds], na.rm = TRUE),
  TEST_COR = cor(derive_assay$chew_grade[derive_assay_test_inds], derive_assay$chew_grade_derived_yue[derive_assay_test_inds], use = "complete.obs"),
  TRAIN_COR = cor(derive_assay$chew_grade[-derive_assay_test_inds], derive_assay$chew_grade_derived_yue[-derive_assay_test_inds], use = "complete.obs"),
  OVERALL_COR = cor(derive_assay$chew_grade, derive_assay$chew_grade_derived_yue, use = "complete.obs")
))

# ^^^ Strict performs better

# fix infinites
derive_assay <- derive_assay %>% 
  mutate_all(~ifelse(is.infinite(.), NA, .))

# add to colData
sce$chew_grade_derived_mine_strict <- as.numeric(derive_assay$chew_grade_derived_mine_strict[match(sce$sample, rownames(derive_assay))])
sce$chew_grade_derived_mine_loose <- as.numeric(derive_assay$chew_grade_derived_mine_loose[match(sce$sample, rownames(derive_assay))])
sce$chew_grade_derived_yue <- as.numeric(derive_assay$chew_grade_derived_yue[match(sce$sample, rownames(derive_assay))])



# # look at how derived differs from actual
# (simple_derived_overall_confusion <- table(predicted = sce$chew_grade_derived_simple, actual = sce$chew_grade))
# sum(simple_derived_overall_confusion)
# calc_accu(simple_derived_overall_confusion)
# 
# (strict_derived_overall_confusion <- table(predicted = sce$chew_grade_derived_strict, actual = sce$chew_grade))
# sum(strict_derived_overall_confusion)
# calc_accu(strict_derived_overall_confusion)
# 
# (loose_derived_overall_confusion <- table(predicted = sce$chew_grade_derived_loose, actual = sce$chew_grade))
# sum(loose_derived_overall_confusion)
# calc_accu(loose_derived_overall_confusion)


# ((((((((((((((((()))))))))))))))))----
# TEST IMPUTING CHEW GRADE WITH MICE ----
# ((((((((((((((((()))))))))))))))))


params <- data.frame(m = 20, method = c("pmm", "rf", "norm.boot",  "lasso.norm", "norm"))
params <- rbind(params, params %>% mutate(m = 30))
results <- data.frame(rmse = numeric(0), cor = numeric(0), accu = numeric(0))

results <- mclapply(seq(nrow(params)), function(i){
  
  method = params[i, "method"]
  m = params[i, "m"]
  type = params[i, "type"]
  
  cat("\n\n", method, "\n\n")
  # # if mice wants a factor type, change to factor
  # if(type == "factor") test[, "chew_grade"] <- factor(test[, "chew_grade"])
  
  # impute
  set.seed(12030)
  imp <- mice(imp_train_X, m = m, method = method)
  
  # average iterations
  imp_sum <- complete(imp);  for(i in 2:imp$m) imp_sum <- imp_sum+ complete(imp, i); imp_mean <- imp_sum / imp$m
  
  
  # grab where matching
  predicted = round(imp_mean[test_inds, "chew_grade"])
  actual = GT[test_inds,"chew_grade"]
  
  # plot(predicted~actual)
  
  rmse <- RMSE(pred = predicted,obs = actual)
  
  cor <- cor(predicted,actual)
  
  conf <- table(predicted , actual)
  
  accu <- sum(diag(conf))/sum(conf)
  
  return(data.frame(rmse = rmse, cor = cor, accu = accu))
  
  }, mc.cores = detectCores()-1) %>% 
  do.call(rbind, .)

# which is best?
all_results <- cbind(params, results)
all_results$score = with(all_results, {(cor*accu*(1/rmse))^1/3})
all_results[order(-all_results$score),]
# RF is better at m = 30 and 20; and still best when correcting ISOS_LOC


# USE best method to impute chew_grade
m = 30; method = "pmm"
imp_mean <- mclapply(seq(m), function(n){
  complete(mice(imp_train_X, m = 1, method = method), 1)
}, mc.cores = detectCores()-1) %>% 
  Reduce('+', .) %>% 
  {./m}


# IMPUTED chew stats
(imputed_stats <- list(
  CONFUSION = (imputed_confusion <- table(predicted = round(imp_mean$chew_grade[test_inds]),
                              actual = test_y)
               ),
  N = sum(imputed_confusion),
  ACCU= calc_accu(imputed_confusion),
  RMSE = RMSE(imp_mean$chew_grade[test_inds], 
              test_y), 
  COR =  cor(imp_mean$chew_grade[test_inds], 
             test_y)
))

# ((((((((((((((()))))))))))))))
# TEST XGBOOST PREDICTION ----
# ((((((((((((((()))))))))))))))


# Create a DMatrix (optimized data structure for XGBoost)
xgb_fit <- xgboost::xgboost(data = data.matrix(XGB_train_X), label = XGB_train_y, nrounds = 15, 
                            params = list(
                              # eta =  .6
                              # max_depth = 3, 
                              # subsample = .5, 
                              # colsample_bytree = 1
                            ))

# Make predictions
XGB_test_pred <- predict(xgb_fit, data.matrix(XGB_test_X))

# XGB chew stats
(xgb_stats <- list(
  CONFUSION = (xgb_confusion <- table(Predicted = round(XGB_test_pred), Actual = test_y)) ,
  N = sum(xgb_confusion),
  ACCU = calc_accu(xgb_confusion),
  RMSE = RMSE(XGB_test_pred, test_y),
  COR = cor(XGB_test_pred, test_y)
))


  

# (((((((((((((((())))))))))))))))
# USING XGB TO impute chew grade ----
# (((((((((((((((())))))))))))))))

# make fit on all available data to deploy on where chew not observed
xgb_fit <- xgboost(
  data = data.matrix(select(GT, -chew_grade)), 
  label = GT$chew_grade,  
  nrounds = 15
)
  
## select vars for chew grade and filter for where they are not missing
xgb_prediction_assay <- assays(sce)$measure_vars %>% 
  data.matrix() %>% 
  t() %>% 
  cbind(chew_grade = sce$chew_grade) %>% 
  as.data.frame() %>% 
  select(any_of( sort(chew_vars , T)), chew_grade) %>% 
  filter(if_all(-chew_grade, ~ !is.na(.))) %>% # ensure all chew_vars filled for prediction
  mutate(chew_grade_xgb = 
           data.matrix(select(., -chew_grade)) %>% 
           predict(xgb_fit, .)
         )

sce$chew_grade_xgb <- xgb_prediction_assay$chew_grade_xgb[match(sce$sample, rownames(xgb_prediction_assay))]

# sce$chew_grade <- factor(sce$chew_grade)
# 
# sce$chew_grade_xgb <- factor(round(sce$chew_grade_xgb))
# 
# sce$chew_grade_derived <- factor(round(sce$chew_grade_derived))



# (((((((((((())))))))))))
# LOOK AT XGB CHEW SHAPS ----
# (((((((((((())))))))))))


params <- list(objective = "multi:softprob", num_class = n_distinct(GT$chew_grade), nthread = 1)
dtrain <- xgb.DMatrix( data.matrix(select(GT,-chew_grade)) , 
                       label = as.integer(GT$chew_grade),
                       nthread = 1)
set.seed(2003)
fit <- xgb.train(params = params, data = dtrain, nrounds = 50)

# Create "mshapviz" object (logit scale)
library(shapviz)
(shap_obj <- shapviz(fit, X_pred = data.matrix(select(GT,-chew_grade)), X = GT))

names(shap_obj) <- as.character(0:6)

# SHAP plots
(beeswarm <- sv_importance(shap_obj, max_display = 20, kind = "beeswarm", alpha = .4))
ggsave(paste0(wd_path, "results/clustering/09_clustering/9.2.2_deriveChewGrade/chewVarShap_beeswarm.pdf"),beeswarm, height = 10, width = 15)

# (((((((((((((((((((())))))))))))))))))))
# MAKE XGB CHEW KEY FROM SHAP VALUES ----
# (((((((((((((((((((())))))))))))))))))))

foo()

#+ for each shap object, grab the matrix of shaps and find which variables have a importance (sd>threshold)
#+ take these vars and average their value where their shap scores are >threshold

## how varibale must the shaps be to be important?
# shap_sd_threshold <- .17
# shap_sd_threshold <- .777773146
shap_sd_threshold <- 0.7708839
# shap_sd_threshold <- 0.73


## what threshold of importance to use when averaging a variables values; values of variable where shap > threshold
# shap_threshold <- .04
# shap_threshold <- -1.7291773
shap_threshold <- -1.9516737
# shap_threshold <- -1.76

# # *******
# objective <- function(par) {
#
#   shap_sd_threshold <- par[1]
#   shap_threshold <- par[2]
# # *******

imp_var_values <- list()
for(obj_name in names(shap_obj)){
  # a shap object
  a_shap <- shap_obj[[obj_name]]

  # which variables are most important to classify into this group?
  most_important <- colnames(a_shap$S)[colSds(a_shap$S)>shap_sd_threshold]

  most_important_values <- vector(length = 0, mode = "numeric")
  # for each variable, average their  values wehre shap >0
  for(imp_var in most_important){
    # what is the average value of this variable where its shap is >0?
    most_important_values[[imp_var]] <- mean(a_shap$X[ a_shap$S[ , imp_var] >shap_threshold , imp_var])
  }

  imp_var_values[[obj_name]] <- most_important_values

}

# grab variables used in classification
overall_imp_vars <- imp_var_values %>% sapply(names ) %>% unlist %>% unique()

# make matrix of var x grade
XGB_grade <- matrix(0, ncol = 7, nrow = length(overall_imp_vars), dimnames = list(sort(overall_imp_vars), as.character(0:6)) )

# fill matrix with averages that classify to each grade
for(coln in colnames(XGB_grade)){
  XGB_grade[ names(imp_var_values[[coln]]) , coln] <- imp_var_values[[coln]]
}

custom_colors <- c(
  rgb(0, .5, 1),
  rgb(0, 1, 0),
  rgb(1, 1, 0)
  # rgb(1,.5,0)
  # rgb(1, 0, 0)
)
custom_colors_gradient <- colorRampPalette(custom_colors)
# custom_colors <- c( "#0080FF","#00FF00","#FFFF00")

XGB_grade %>%
  as.data.frame() %>%
  rownames_to_column("Field") %>%
  pivot_longer(-Field, names_to = "Grade", values_to = "Valid_inputs") %>%
  
  ggplot(aes(factor(Grade), Field, fill = factor(Valid_inputs %>% round()), label = Field))+
  # scale_fill_gradientn(colors = custom_colors_gradient(10))+
  # scale_fill_manual( values= custom_colors)+
  scale_fill_viridis_d(option = "plasma")+
  geom_tile(alpha = .9)+
  geom_text(fontface = "bold", color = "darkgrey")+
  scale_x_discrete(breaks = factor(0:6))+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y =  element_blank(),
        panel.grid = element_blank())+
  labs(fill = "Valid Inputs")+
  ggtitle(paste("XGB Grade Criteria",round(shap_sd_threshold,3) ,round(shap_threshold,3), sep = ','))
ggsave(paste0(wd_path, "results/clustering/09_clustering/9.2.2_deriveChewGrade/XGB_derive_scheme.png"), height = 5, width = 10)

# LOOK AT HOW USING THIS SCHEME CORRELATES TO ACTUAL
## MANUAL APPROACH following XGB shap important scheme

# key for xgb scheme
XGB_grade_key <- round(XGB_grade)
XGB_grade_key[is.nan(XGB_grade_key)] <- 0

# grab vars important to classifying
XGB_derive_assay <-
  assays(sce)$measure_vars %>%
  data.matrix() %>% 
  t() %>%
  cbind(chew_grade = sce$chew_grade) %>%
  as.data.frame() %>%
  select(any_of( sort(chew_vars , T)), chew_grade) %>%
  # filter for where all vars not missing allowing  for chew grade to be missing
  drop_na(-chew_grade)

# reorder columns
XGB_derive_assay <- XGB_derive_assay[ , c(rownames(XGB_grade_key), setdiff(colnames(XGB_derive_assay),rownames(XGB_grade_key)) )  ]

# add new derived chew grade
XGB_derive_assay <- cbind(XGB_derive_assay, chew_grade_XGB_derived = NA)

# suppressWarnings(
#   colnames(XGB_derive_assay)[colnames(XGB_derive_assay) != rownames(XGB_grade_key)]
#   )
#
# head(XGB_derive_assay)

# Iterate through grades
for (grade in 0:6) {
  # Compare each row of the assay with the grade key vector
  key <- XGB_grade_key[, as.character(grade)]
  key <- key[key>0]
  x <- XGB_derive_assay[, names(key), drop = FALSE]

  where <- apply(x, 1, function(row){
    all(row == key)
    }) %>%
    which()

  XGB_derive_assay[where, "chew_grade_XGB_derived"] <- grade
}

# # *******
#   CONFUSION = (xgb_derived_confusion <- table(predicted =  XGB_derive_assay$chew_grade_XGB_derived[derive_assay_test_inds],
#                                               actual = XGB_derive_assay$chew_grade[derive_assay_test_inds]))
#   ACCU = calc_accu(xgb_derived_confusion)
#   RMSE = RMSE(pred = XGB_derive_assay$chew_grade_XGB_derived[derive_assay_test_inds],
#               obs = XGB_derive_assay$chew_grade[derive_assay_test_inds], na.rm = T)
#   TEST_COR = cor(XGB_derive_assay$chew_grade[derive_assay_test_inds],
#                  XGB_derive_assay$chew_grade_XGB_derived[derive_assay_test_inds], use = "complete.obs")
#   # SCORE = (ACCU+TEST_COR)/RMSE
#   SCORE = (TEST_COR)/RMSE
#
#   # return( RMSE )
#   return( -SCORE )
# }
#
# # Define bounds for the parameters
# lower_bounds <- c(shap_sd_threshold = 0, shap_threshold = -2)
# upper_bounds <- c(shap_sd_threshold = 2, shap_threshold = 2)
#
# # Run GenSA
# optimized <- GenSA::GenSA(
#   par = c(shap_sd_threshold = 1, shap_threshold = 0),  # Initial parameter guess
#   lower = lower_bounds,
#   upper = upper_bounds,
#   fn = objective,
#   control = list(max.call = 1000, temperature = 1000)
# )
# optimized
# # shap_sd_threshold    shap_threshold
# # 0.7708839        -1.9516737
# # *******


# make summary list of performance metrics
(xgb_derived_stats <- list(
  CONFUSION = (xgb_derived_confusion <- table(predicted =  XGB_derive_assay$chew_grade_XGB_derived[derive_assay_test_inds],
                                              actual = XGB_derive_assay$chew_grade[derive_assay_test_inds])),
  N = sum(xgb_derived_confusion),
  ACCU = calc_accu(xgb_derived_confusion),
  RMSE = RMSE(pred = XGB_derive_assay$chew_grade_XGB_derived[derive_assay_test_inds],
              obs = XGB_derive_assay$chew_grade[derive_assay_test_inds], na.rm = T),
  TEST_COR = cor(XGB_derive_assay$chew_grade[derive_assay_test_inds],
            XGB_derive_assay$chew_grade_XGB_derived[derive_assay_test_inds], use = "complete.obs"),
  TRAIN_COR = cor(XGB_derive_assay$chew_grade[-derive_assay_test_inds],
            XGB_derive_assay$chew_grade_XGB_derived[-derive_assay_test_inds], use = "complete.obs")
))

# add to colData
sce$chew_grade_XGB_derived <- as.numeric(XGB_derive_assay$chew_grade_XGB_derived[match(sce$sample, rownames(XGB_derive_assay))])

# overall accuracy 
(xgb_overall_confusion <-  table(predicted  = round(sce$chew_grade_XGB_derived), actual = sce$chew_grade))
sum(xgb_overall_confusion)
calc_accu(xgb_overall_confusion)
cor(sce$chew_grade_XGB_derived, sce$chew_grade, use = "complete.obs")
RMSE(sce$chew_grade_XGB_derived, sce$chew_grade, na.rm = TRUE)

# ((((((((()))))))))
# ALL STATS ----
# ((((((((()))))))))
(all_stats_l <- list(
  chew_grade_derived_mine_loose = chew_grade_derived_mine_loose_stats,
  chew_grade_derived_mine_strict = chew_grade_derived_mine_strict_stats,
  chew_grade_derived_yue_= chew_grade_derived_yue_stats,
  imputed = imputed_stats,
  xgb_derived = xgb_derived_stats, 
  xgb= xgb_stats) 
 )


all_stats <-  all_stats_l %>% 
    sapply(function(list){
      l_names <- grep("ACCU|RMSE|^COR$|TEST_COR|^N$" , names(list), value = TRUE)
      m <- as.matrix(list[l_names])
      if("TEST_COR"%in% rownames(m)) {
        m[["COR"]] <- m[["TEST_COR"]]
        m[["TEST_COR"]] <- NULL
      }
      return(data.matrix(m))
    }, USE.NAMES = TRUE, simplify = FALSE) %>% 
  do.call(cbind, .) %>% 
  t() 

rownames(all_stats) <- names(all_stats_l)

all_stats %>% 
  as.data.frame() %>% 
  rownames_to_column("method") %>% 
  mutate_all(unlist) %>% 
  mutate(Reciprocal_RMSE = 1/RMSE) %>% 
  select(-RMSE) %>% 
  melt() %>% 
  ggplot(aes(method, value)) +
  facet_wrap(~variable, ncol = 1, scales = "free") +
  geom_bar(aes(fill = method), stat = "identity") +
  geom_text(aes(label = round(value, 2)), vjust = 1.5) +  # Add value labels
  ggtitle("Different Methods Chew Grade Prediction Results") +
  theme_bw()
ggsave(paste0(wd_path, "results/clustering/09_clustering/9.2.2_deriveChewGrade/Chew_grade_pred_bars.png"), height = 15, width = 10)


# ((((((((()))))))))----
# CALCULATE REMAINDER OF CHEWGRADE_MINE_LOOSE ----
# ((((((((()))))))))

sce$chew_grade_derived_mine_strict %>% is.na() %>% sum
sce$chew_grade_derived_mine_loose%>% is.na() %>% sum

# un-scale vars to orig distrs
means <- metadata(sce)$measure_vars_scale_attr$`scaled:center`[intersect(chew_vars, rownames(sce))]
sds <- metadata(sce)$measure_vars_scale_attr$`scaled:scale`[intersect(chew_vars, rownames(sce))]

# grab vars at missing chew grade
imputed_assay <- assays(sce)$pmm_imp_of_scaled[names(means),is.na(sce$chew_grade_derived_mine_loose)] %>% 
  t()  %>% 
  as.matrix()
  
# un-scale matrix; multiply each column by respective sds; add each column to respective mean; round to nearest int
unscaled_imputed_assay <- imputed_assay %>% 
  sweep(., 2, sds, "*") %>% 
  sweep(., 2, means, "+")  %>% 
  round() %>% 
  as.data.frame()

# simplify IR02-9 to rowmax of zones 2:max zone
unscaled_imputed_assay <- unscaled_imputed_assay %>% 
  mutate(`IR02-9` = rowMaxs(as.matrix(select(as.data.frame(.), matches("IR0[23456789]$"))))) %>%
  select(-c(IR02:IR05))

unscaled_imputed_assay$chew_grade_derived_mine_loose <- NA

# Grade 0
unscaled_imputed_assay$chew_grade_derived_mine_loose[
  unscaled_imputed_assay$ISOS_BREAK == 0 &
    unscaled_imputed_assay$`IR02-9` == 0 &
    unscaled_imputed_assay$IR01 == 0 &
    unscaled_imputed_assay$HR_TEMPORAL == 0
] <- 0


# GRADE 1
unscaled_imputed_assay$chew_grade_derived_mine_loose[
  unscaled_imputed_assay$ISOS_LOC_JF==1 &
    unscaled_imputed_assay$ISOS_BREAK%in%c(1,2)
] <- 1

# GRADE 2
unscaled_imputed_assay$chew_grade_derived_mine_loose[
  unscaled_imputed_assay$ISOS_LOC_CEN==1 &
    unscaled_imputed_assay$ISOS_BREAK%in%c(1,2)
] <- 2


# GRADE 3
unscaled_imputed_assay$chew_grade_derived_mine_loose[
  unscaled_imputed_assay$`IR02-9` %in%c(1,2)
] <- 3

# GRADE 4
unscaled_imputed_assay$chew_grade_derived_mine_loose[
  unscaled_imputed_assay$HR_TEMPORAL %in%c(1,2)
] <- 4

# GRADE 5
unscaled_imputed_assay$chew_grade_derived_mine_loose[
  unscaled_imputed_assay$IR01 %in%c(1,2)
] <- 5

# GRADE 6
unscaled_imputed_assay$chew_grade_derived_mine_loose[
  unscaled_imputed_assay$NEOVASC %in%c(1,2)
] <- 6


# grade for unclassifiable obs
unscaled_imputed_assay$chew_grade_derived_mine_loose[with(unscaled_imputed_assay,
                                                (ISOS_BREAK > 0 |
                                                   `IR02-9` > 0 |
                                                   NEOVASC >0 |
                                                   HR_TEMPORAL > 0) &
                                                  is.na(chew_grade_derived_mine_loose)
)]<- -Inf
  
unscaled_imputed_assay$chew_grade_derived_mine_loose %>% 
  table(useNA = "always") %>% 
  # sum %>% # >> 181
  {}

# add chew_grade to coldata

colData(sce)$chew_grade_derived_mine_loose[match(rownames(unscaled_imputed_assay),sce$sample) ] %>% length()
# 181

# check rownames
all(rownames(colData(sce)[match(rownames(unscaled_imputed_assay),sce$sample), c("chew_grade_derived_mine_loose"), drop = FALSE]) == rownames(unscaled_imputed_assay))


colData(sce)$chew_grade_derived_mine_loose[match(rownames(unscaled_imputed_assay),sce$sample)] <- unscaled_imputed_assay$chew_grade_derived_mine_loose

# (((((((((((())))))))))))----
# SAVE DATA ----
# (((((((((((())))))))))))
# save SCE
saveHDF5SummarizedExperiment(sce, paste0(wd_path, "processed_data/sce002.2.h5"), replace = T)



# MUNG  ----

# library(future)
# plan(multisession)
# how.long(T)
# f1 %<-% {
#   Sys.sleep(4)
#   2
# }
# 
# f2 %<-% {
#   Sys.sleep(4)
#   6
# }
# how.long()
# 
# Sys.sleep(4)
# 
# how.long(T)
# resolved(f1)
# resolved(f2)
# how.long()
# 
# f1+f2
