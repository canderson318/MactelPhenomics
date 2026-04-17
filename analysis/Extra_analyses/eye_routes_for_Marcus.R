

# would you be able to tell me which of the UK patients fall clearly into one of the two progression routes you showed in NYC? I'm setting up my Nightingale ptient selection and want to make sure they are well spread over different disease modalities.  

# Christian, can you extract these?
# Only use the A and the C classifications.
# Ask Lea how to understand if they are form the UK or not.
# Also in the database these people should also have another ID like UK123.
# Please provide Marcus with both the MacTel and the UKid if available in our data.

library(tibble)
library(ggplot2)
library(dplyr)
library(HDF5Array)

rm(list = ls()); gc()

# script with all functions i have made
source("/Users/canderson/Documents/LMRI/PRJ2024001/version015/analysis/00_MyFunctions.R")

# LOAD DATA ----
wd_path = "~/Documents/LMRI/PRJ2024001/version015/"

how.long()
# load experiment
sce <-  loadHDF5SummarizedExperiment(paste0(wd_path, "processed_data/sce006.h5"))
how.long()

# patient IDs sheet
patIDs <- read.csv("/Users/christiananderson/Library/CloudStorage/Egnyte-lmri/Shared/LMRI/DS_Group/PRJ2024001/raw_data/PatientIDs.csv") %>% 
  clean.it() %>% 
  rename_with(~c("patID", "A_number", "Site", "ID4"))
patIDs<- patIDs[-1,]


# GET ROUTES FOR EACH EYE ----
eye_routes <- colData(sce) %>% 
  as.data.frame() %>% 
  select(patID, patID_eye, rf_route_score, rf_route) %>% 
  unique()

# UK Site IDs from raw sheet
raw_UK_pats <- metadata(sce)$raw$DEM %>%
  select(matches("patientId|patId3")) %>%
  rename_with(~c("patID", "Site")) %>%
  drop_na() %>%
  distinct()%>%
  filter(grepl("^0*8", patID)) 

# UK site ids from database
db_UK_pats <- patIDs %>%
  select(patID, Site) %>%
  drop_na() %>%
  distinct() %>%
  filter(grepl("^0*8", patID)) 

sets <- list(
  db_UK_pats = db_UK_pats$Site,
  raw_UK_pats = raw_UK_pats$Site
)

par(mar = rep(0,4))
plot.new()
VennDiagram::venn.diagram(sets, filename = NULL, disable.logging = TRUE) %>%
  grid::grid.draw()
par(mar = rep(4,4))


# # all SItes match where in both
# all(raw_UK_IDs$Site[raw_UK_IDs$patID %in% db_sites$patID] == db_sites$Site[db_sites$patID %in% raw_UK_IDs$patID])

# add site ids and others to eye_routes; filter for uk sites
eye_routes <- eye_routes %>% 
  left_join(db_UK_pats) %>% 
  left_join(patIDs) %>% 
  filter(grepl("^0*8", patID)) %>% 
  arrange(desc(abs(rf_route_score)))
 
eye_routes[eye_routes$rf_route%in% c("A", "B"), c("patID",  "patID_eye", "Site", "rf_route")] %>%
  drop_na() %>%
  unique


# wide table of patients routes and scores
wide_eye_routes <-
  eye_routes %>% 
  drop_na(Site, rf_route) %>% 
  separate(patID_eye, c("patID", "eye")) %>% 
  mutate(eye = case_when(eye == "1" ~ "OD", eye == "2" ~ "OS", TRUE ~ eye)) %>% 
  pivot_wider(names_from = eye, 
              values_from = c(rf_route, rf_route_score), 
              names_prefix = "", 
              values_fill = list(rf_route_score = NA, rf_route = NA)
              ) %>% 
  mutate(rf_route_score_average = rowMeans(across(matches("score")))) %>% 
  arrange(desc(abs(rf_route_score_average)))


# patients who follow A or B
AorB <- wide_eye_routes[wide_eye_routes$rf_route_OD %in% c("A", "B") | wide_eye_routes$rf_route_OS %in% c("A", "B") , ]

# select top ## and add eyes classified as A or B 
patients_for_marcus <- wide_eye_routes %>% 
  rbind(AorB) %>% 
  arrange(desc(abs(rf_route_score_average) ))%>% 
  distinct() %>% 
  mutate(across(where(is.numeric), ~round(.,2))) %>%
  slice_head(n = 50) %>% 
  ungroup()

# grab tsne data and join patient info; find lead tsne x,y for segments; 
tsne <- 
  reducedDim(sce, "TSNE") %>% 
  as.data.frame() %>% 
  rownames_to_column("subject") %>% 
  left_join(
    colData(sce) %>% as.data.frame() %>% select(subject, patID, patID_eye, rf_route, rf_route_score, month)
  ) %>% 
  arrange(patID_eye, month) %>% 
  mutate(across(matches("TSNE"), 
                ~lead(.), 
                .names = "{.col}_end"), 
         .by = patID_eye) 

# select top ## eyes of patients to plot
eyes_of_interest <- eye_routes %>% 
  filter(patID %in% patients_for_marcus$patID) %>% 
  drop_na(rf_route_score,Site) %>% 
  arrange(desc(abs(rf_route_score))) %>% 
  slice_head(n = 30) 
  
p <- tsne %>% 
  filter(patID_eye %in% eyes_of_interest$patID_eye) %>% 
  # filter for where points not too close together
  filter(abs(TSNE1-TSNE1_end) > .5 & abs(TSNE2-TSNE2_end) > .5 & !is.na(TSNE1_end)) %>%
  ggplot(aes(TSNE1, TSNE2, color = rf_route_score, group = patID_eye)) +
  geom_point(data = tsne, aes(TSNE1, TSNE2), color = "grey", alpha = .4) +
  scale_color_gradient2() +
  geom_segment(aes(xend = TSNE1_end, yend = TSNE2_end), 
               # linewidth = 1,
               arrow = grid::arrow(angle = 10, type = "closed", length = unit(.02, "npc"))) +
  geom_point()+
  theme_bw() +
  labs(color = "Route Score") +
  theme(aspect.ratio = 1)+
  ggtitle("Top 10 UK eyes with greatest route certainty")

ggsave("~/Documents/LMRI/PRJ2024001/version015/results/clustering/09_clustering/top_UK_eyes_trajectories.png",
       p,
       height = 8, 
       width = 9, 
       bg = "white", 
       dpi = 150)

names(patients_for_marcus) <- gsub("rf_", "",names(patients_for_marcus))

# save data 
write.csv(patients_for_marcus, "~/Documents/LMRI/PRJ2024001/version015/processed_data/patients_for_marcus.csv",row.names = FALSE)



