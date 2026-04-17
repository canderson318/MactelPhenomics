# Script Analysis Descriptions

## 01_preprocessing.Rmd
In this script tables downloaded from the LMRI db were cleaned and preprocessed to join into on long-ish format dataset where each row is an observation of a patient's eye at a visit. 
- cleaned off artifact characters
- make names cohesive
- remove 'bad image' codes
- format dates
	- collapse multiple dates into one
- lengthen VAA
- join all dfs together on fuzzy date window
- filter out Controls 
- remove observations with only demographic/VAA data

## 02_visits_dates_missingness
This script takes the output of the previous script and cleans for specific missingness and wacky visits. 
- make info data dictionary
- try to infer visitID from date for visits like 084S, 000unsch, ...
- if no visitID could be infered and those visits don't have any salient data remove those observations

## 03_MT_imputation
this script imputes the OCT macular thickness variables and creates a composite variable for each ETDRS zone from the three machine types. 
- map ETDRS zones to right scheme so 5 == TINN, etc.
- remove outlier MT values
- impute within machine using mice pmm
- rescale cirrus/stratus to spectralis distr. 
- take the average of rescaled_cirrus/straatus and spectralis for each zone

## 04_fix_ordinal_vars_and_ETDRS_plot.Rmd
this script fixes some variables incorrect ordinality and makes a mt etdrs plot
- one-hot encode ISOS_loc 
- fix OTHER_DISEASE, FAZ, YELLOWRND, BLRPATTERN, var>2~2
- plot

## 05_longitudinal_MT_EDA.Rmd
this script makes some MT over time plots

## 06_zonal_progression_plot.Rmd
this script makes a zonal progression plot 

## 07_VAA_cleaning_and_FA.Rmd
this script cleans VA data further
- calculate missing logmar/snellen/acuit using their relationships
- impute VAA using mice pmm

## 08_longitudinal_function.Rmd
this script makes some rate variables. not an important script

## 09.1_imp_and_dimred.Rmd
this script takes the clean data set and creates an sce object with the pertinent measure vars as the assay, imputes the whole assay, and calculates reduced dimensions.
- filter for non-redundant structural variables (includes VAA rowscores)
- filter out variables that are too missing and ETDRS>6 (with some exception)
- filter out samples that are too missing
- mean center and scale
- use  mice to impute whole assay with pmm
- fix some vars that don't impute exactly calculate slope from imputed MT
- calculate pca>(umap,tsne)


## 09.1.1_feature_correlation_heatmap.R
this script makes a heatmap of all the feature correlations in sce assay

## 09.2.1_clustering.Rmd
this script takes the imputed data in sce and clusters on PCs. then calculates the MST and infers pseudotime. cluster-cluster transition percentages support MST trajectories. 
- KM-KNN cluster PCA
- combine clusters that are similar: 16 down to 13 clusters
- make cluster phenotypic profile plot as cluster aggregate by mean
- create MST and pseudotime
- calculate from-to and to-from cluster-cluster percentages
- make sequential TSNE/PCA plots

## 9.2.2_deriveChewGrade.R
this script takes the Chew scheme and infers a chew grade from the corresponding features in the data. multiple schemes could make sense, these are compared

## 9.2.3_clusteringRobustnessCheck.Rmd
this script manipulates/subsets the data in several ways to see how much the clustering/trajectories change when using the same pipeline
- permute clustering and overlay MSTs
- run pipeline on longitudinal, singleton, subsampled, noisy,  VAA-less assays
- plot results

## 09.2.3.1.Rmd
building off prev script, this script runs the pipeline stratifying by continent of origin. 

## 09.2.4_networkFlowValidation.Rmd
this script tests how well observed transitions align to the MST, and several hypothetical graphs. 
- make transition matrix from cluster-cluster transitions 
- generate linear graphs based on whatever ordering, e.g. chewgrade_mean. these allow for parralell progression
- make simple, linear only graphs from chew and pp
- compare G-statistic of observed versus expected bifurcated, chew, pp, unordered graphs
- compare these fits to a set semi-random matrices. the rowsums of this matrix match observed rowsums

## 09.3.1_featureClusterShap.R
this script uses a cluster~Features xgboost model to create shap scores that profile clusters phenotypically
- make test and train dataset
- train xgboost model 
- visualize shapscores
- use model to predict cluster probs per sample
- identify highest prob samples per cluster
- plot their pheno profiles

## 09.3._featureMarkerGeneAnalysis.R
this script tries out marker gene analysis to profile clusters. not very useful

## 10_functionAndStructureModeling.Rmd
this notebook models how macular features associate with rate of change in vision/MT
- calculate rates
- model mtrate/varate ~ phenos
- model parafovea/parafovea_rate~phenos
- do the same for slope

## 10.2_functionAndStructureModeling.Rmd
this script runs the same models as before but uses a lagged rate model where inter-visit change versus first-last change models against phenotypes that might be driving progression. *NOT USED

## 11_clusterModeling.Rmd
this script models how phenotypes associate to rate of change in pseudoprogression from one visit to the next using a stepwise rate model (to increases sample size). Also modeling what phenotypes 'push' eyes into their next cluster
- make modeling data to have each visit paired with every subsequent visit. from these pairs model how phenotype presence relates to (vis_i+1 - vis_i) / (t_i+1 - t_i)
- make next_cluster dataseet and model hwat phenotypes are present for each next cluster

## 11.1_cluster_representatives.R
This script is a copy from 09.3.1_featureClusterShap.R that just focuses on identifying samples that are representative of each cluster. 
- make test/train data
- train xgb
- predict cluster and extract probs
- identify highest probs per cluster

## 12_clusterExploration.Rmd
this script classifies each data point to one of three trajectories, A, B, or X for neither A or B. 
- manually label paths
- engineer training data from PCs
- try RF and LDA
- use LDA to predict remaining datapoints traj from PC data

## 13_routeModeling.Rmd
this script regresses trajectory on phenotypes, covariates, and snps; chew, acuit, and pseudoprogression progression rate on phenos; and phenos on snps.
- plot progression bin plots
- correlate/chisq R/L eye trajectories
- model route ~ first_phenos
- model route ~ covariates
- model route ~ snps
- model route ~ snp_covariates
- model acuit_PR ~ first_phenos
- model acuit_PR ~ covs
- model Chew_PR ~ first_phenos
- model Chew_PR ~ covs
- model pseudoprogression_PR ~ first_phenos
- model pseudoprogression_PR ~ covs
- model phenos ~ snps
- make p,q,and pq snp_pheno estimate heatmap

## 13.2_T2D_Contrasts.Rmd
this script tries to associate trajectory to T2D by using contrasts for contrast of Aclusters to Bclusters.
- make different contrasts, endL-endR, allL-allR, ...
- try lmer contrasts
- try glmer contrasts
- try limma contrasts

## 13.3_BMI_Contrasts.Rmd
this script is exactly like last but for BMI differences between trajectories

## 13.3_SNP_Contrasts.Rmd
this script is exactly like last but for BMI differences between trajectories

## 14_GeneticProgressionModeling.Rmd
this script explores the relationship of progression rate to snps and how that relationship changes at different stages of progression. 
- look at demographics
- do families follow trajectories
- model pp_rate ~ Snps
- model pp_rate ~ covariates
- model pp_rate ~ PRS + stage
- model pp_rate ~ pp_rate_PRS + stage

## 14.1_progressionTimeModeling.Rmd
this script tries to estimate how much time passes along pseudoprogression by using a stepwise time change linear model with pp quantiles.
- make chew,next_chew dataframe for stepwise chew rate of change model
- model chew_change ~ year_change*chew_stage
- take the reciprocal of effects for year~chew relatinoship
- model pp_change ~ year_change*bl_pp_quantile_stage
- take reciprocal of effects to estimate year~pp
- calculate cluster pp_averages and estimate time from cluster to cluster

## 14.1.1_pseudoprog_time_shif_brm.Rmd
this script tries to tackle teh same question as the previous using bayesian methods. the model estimates how much time since disease genesis, creating a population disease progression profile.
- estimate individual level random intercepts for temporal shifts using brms
- infer time change along pp using LoBF
- estimate time between clusters from LoBF

## 14.1.1_chew_time_shif_brm.Rmd
same thing as pseudoprog but for chew

## 14.2_geneticContrastsModels
no code, trying to estimate snp relationships to clusters/trajectories

## 15_estimateClusterFeatureDistributions
tries to estimate underlying distributions behind cluster feature means. 
- for each cluster, for each featuere, estimate mean using MLE to find CIs
- plot cluster feature means with CIs











