####################################################
## Read Phenotype data from 2_genetic_data_P_LP.r ##
####################################################
PHENO.ANY_SN <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/2_genetic_data_P_LP.rds")
# PHENO.ANY_SN <- readRDS("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/2_genetic_data_P_LP.rds")


#########################################################
## Read PRS score files and merge them to PHENO.ANY_SN ##
#########################################################
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/prs_out")
# setwd("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/prs_out")

## Meningioma_Claus_prs.profile
meningioma_Claus <- read.table("Meningioma_Claus_prs.profile", header = T)
PHENO.ANY_SN$Meningiom_Claus_PRS <-  meningioma_Claus$SCORE [match(PHENO.ANY_SN$sjlid, meningioma_Claus$IID)]

## Meningioma_Dobbins_prs.profile
Meningioma_Dobbins <- read.table("Meningioma_Dobbins_prs.profile", header = T)
PHENO.ANY_SN$Meningioma_Dobbins_PRS <-  Meningioma_Dobbins$SCORE [match(PHENO.ANY_SN$sjlid, Meningioma_Dobbins$IID)]

## Pleiotropy_Bi_directional_Increasing_prs.profile
Pleiotropy_Bi_directional_Increasing <- read.table("Pleiotropy_Bi_directional_Increasing_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_Bi_directional_Increasing_PRS <-  Pleiotropy_Bi_directional_Increasing$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_Bi_directional_Increasing$IID)]


## Pleiotropy_Bi_directional_Increasing_Significant_prs.profile
Pleiotropy_Bi_directional_Increasing_Sig <- read.table("Pleiotropy_Bi_directional_Increasing_Significant_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_Bi_directional_Increasing_Sig_PRS <-  Pleiotropy_Bi_directional_Increasing_Sig$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_Bi_directional_Increasing_Sig$IID)]


## Pleiotropy_Bi_directional_Decreasing_prs.profile
Pleiotropy_Bi_directional_Decreasing <- read.table("Pleiotropy_Bi_directional_Decreasing_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_Bi_directional_Decreasing_PRS <-  Pleiotropy_Bi_directional_Decreasing$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_Bi_directional_Decreasing$IID)]


## Pleiotropy_Bi_directional_Decreasing_Significant_prs.profile
Pleiotropy_Bi_directional_Decreasing_Sig <- read.table("Pleiotropy_Bi_directional_Decreasing_Significant_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_Bi_directional_Decreasing_Sig_PRS <-  Pleiotropy_Bi_directional_Decreasing_Sig$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_Bi_directional_Decreasing_Sig$IID)]


## Pleiotropy_Meta_analysis_prs.profile
Pleiotropy_Meta_analysis <- read.table("Pleiotropy_Meta_analysis_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_Meta_analysis_PRS <-  Pleiotropy_Meta_analysis$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_Meta_analysis$IID)]

## Pleiotropy_One_cohort_prs.profile
Pleiotropy_Meta_analysis <- read.table("Pleiotropy_One_cohort_prs.profile", header = T)
PHENO.ANY_SN$Pleiotropy_Meta_analysis_PRS <-  Pleiotropy_Meta_analysis$SCORE [match(PHENO.ANY_SN$sjlid, Pleiotropy_Meta_analysis$IID)]


#########################
## Extract Ethnicities ##
#########################
PHENO.ANY_SN.EUR <- PHENO.ANY_SN[PHENO.ANY_SN$PCA.ethnicity == 'EUR', -grep("sjlid|PCA.ethnicity|AGE.ANY_SN", colnames(PHENO.ANY_SN))]



# lapply(list(wgspop, wgsdiag, subneo, radiation, drug, demog, adultbmi, adolhabits, adlthabits), dim)
# save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/demographic.RDATA")