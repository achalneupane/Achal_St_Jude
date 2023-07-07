## Phenotype data for CCSS (to replicate GWAS findings in SJLIFE)
rm(list=ls())
# setwd('/Volumes/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno/')
setwd("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/cardiomyopathy_gwas/yadav_pheno/")
## Phenotype data by Huiqi and Qi
# load('ccss_export_05062022.RData')

## European ancestry
# Original
ccss_org = read.table('CCSS.SJLIFE.ancestry_CCSS.txt', header = TRUE, sep="\t")
ccss_org$SAMPLE <- paste(ccss_org$SAMPLE, ccss_org$SAMPLE, sep = "_")
ccss_org_eur = subset(ccss_org, CEU>=0.8)
# Expansion
ccss_exp = read.table('CCSS_and_1kGP_final_EUR_AFR_EAS.3.Q_ccssexp_samples', header = TRUE)
ccss_exp_eur = subset(ccss_exp, EUR>=0.8)
# CCSS eur
ccss_eur = c(ccss_org_eur$SAMPLE, ccss_exp_eur$ccssid)


## African ancestry
# Original
ccss_org = read.table('CCSS.SJLIFE.ancestry_CCSS.txt', header = TRUE, sep="\t")
ccss_org$SAMPLE <- paste(ccss_org$SAMPLE, ccss_org$SAMPLE, sep = "_")
ccss_org_afr = subset(ccss_org, YRI>=0.8)
# Expansion
ccss_exp = read.table('CCSS_and_1kGP_final_EUR_AFR_EAS.3.Q_ccssexp_samples', header = TRUE)
ccss_exp_afr = subset(ccss_exp, AFR >=0.8)
# CCSS AFR
ccss_afr = c(ccss_org_afr$SAMPLE, ccss_exp_afr$ccssid)


## Asian ancestry
# Original
ccss_org = read.table('CCSS.SJLIFE.ancestry_CCSS.txt', header = TRUE, sep="\t")
ccss_org$SAMPLE <- paste(ccss_org$SAMPLE, ccss_org$SAMPLE, sep = "_")
ccss_org_asa = subset(ccss_org, ASA >=0.8)
# Expansion
ccss_exp = read.table('CCSS_and_1kGP_final_EUR_AFR_EAS.3.Q_ccssexp_samples', header = TRUE)
ccss_exp_asa = subset(ccss_exp, EAS >=0.8)
# CCSS ASA
ccss_asa = c(ccss_org_asa$SAMPLE, ccss_exp_asa$ccssid)

## Read PCA files
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/yutaka_request_09_27_2022/plink_data/RE__Genetic_data_for_follow-up_cardiomyopathy")


ccss_exp_pca.EUR <- read.table("CCSSEXP_EUR_top_20_PCs.eigenvec.ccssid", header = F)
ccss_exp_pca.AFR <- read.table("CCSSEXP_AFR_top_20_PCs.eigenvec.ccssid", header = F)
ccss_exp_pca.EAS <- read.table("CCSSEXP_EAS_top_20_PCs.eigenvec.ccssid", header = F)
ccss_exp_pca.AMR <- read.table("CCSSEXP_AMR_top_20_PCs.eigenvec.ccssid", header = F)

ccss_org_pca.EUR <- read.table("CCSS.AnalysisSet.1-22_EUR_indep_top_10_pcs.eigenvec", header = F)
ccss_org_pca.EUR$V1 <- paste(ccss_org_pca.EUR$V1, ccss_org_pca.EUR$V1, sep = "_")

# library("readxl")

# Read spreadsheet
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/ccss_exp_wgs/preQC_VCF_per_chromosome/yutaka_request_09_27_2022/plink_data/")

## -----1. CCSS_exp_WGS
CCSS_exp_WGS <- read_excel("CCSS_org_and_CCSS_exp_genotype_V2.xlsx", sheet = "CCSS_exp_WGS")

CCSS_exp_WGS$PCA_ethnicity <- NA
CCSS_exp_WGS$PCA_ethnicity [CCSS_exp_WGS$samples %in% ccss_exp_pca.EUR$V1] <- "European"
CCSS_exp_WGS$PCA_ethnicity [CCSS_exp_WGS$samples %in% ccss_exp_pca.AFR$V1] <- "African"
CCSS_exp_WGS$PCA_ethnicity [CCSS_exp_WGS$samples %in% ccss_exp_pca.EAS$V1] <- "Asian"
CCSS_exp_WGS$PCA_ethnicity [CCSS_exp_WGS$samples %in% ccss_exp_pca.AMR$V1] <- "Admixed_American"

CCSS_exp_WGS$Admixture_ethnicity <- NA
CCSS_exp_WGS$Admixture_ethnicity [CCSS_exp_WGS$samples %in% ccss_eur] <- "European"
CCSS_exp_WGS$Admixture_ethnicity [CCSS_exp_WGS$samples %in% ccss_afr] <- "African"
CCSS_exp_WGS$Admixture_ethnicity [CCSS_exp_WGS$samples %in% ccss_asa] <- "Asian"

write.table(CCSS_exp_WGS, "CCSS_EXP_WGS_ethnicity.txt", col.names = T, quote = F)


## -----2. CCSS_exp_WGS_QCed
CCSS_exp_WGS_QCed <- read_excel("CCSS_org_and_CCSS_exp_genotype_V2.xlsx", sheet = "CCSS_exp_WGS_QCed")

CCSS_exp_WGS_QCed$PCA_ethnicity <- NA
CCSS_exp_WGS_QCed$PCA_ethnicity [CCSS_exp_WGS_QCed$samples %in% ccss_exp_pca.EUR$V1] <- "European"
CCSS_exp_WGS_QCed$PCA_ethnicity [CCSS_exp_WGS_QCed$samples %in% ccss_exp_pca.AFR$V1] <- "African"
CCSS_exp_WGS_QCed$PCA_ethnicity [CCSS_exp_WGS_QCed$samples %in% ccss_exp_pca.EAS$V1] <- "Asian"
CCSS_exp_WGS_QCed$PCA_ethnicity [CCSS_exp_WGS_QCed$samples %in% ccss_exp_pca.AMR$V1] <- "Admixed_American"

CCSS_exp_WGS_QCed$Admixture_ethnicity <- NA
CCSS_exp_WGS_QCed$Admixture_ethnicity [CCSS_exp_WGS_QCed$samples %in% ccss_eur] <- "European"
CCSS_exp_WGS_QCed$Admixture_ethnicity [CCSS_exp_WGS_QCed$samples %in% ccss_afr] <- "African"
CCSS_exp_WGS_QCed$Admixture_ethnicity [CCSS_exp_WGS_QCed$samples %in% ccss_asa] <- "Asian"

write.table(CCSS_exp_WGS_QCed, "CCSS_exp_WGS_QCed_ethnicity.txt", col.names = T, quote = F)

## -----3. CCSS_org_hrc
CCSS_org_hrc <- read_excel("CCSS_org_and_CCSS_exp_genotype_V2.xlsx", sheet = "CCSS_org_hrc")

CCSS_org_hrc$PCA_ethnicity <- NA
CCSS_org_hrc$PCA_ethnicity [CCSS_org_hrc$samples %in% ccss_org_pca.EUR$V1] <- "European"

CCSS_org_hrc$Admixture_ethnicity <- NA
CCSS_org_hrc$Admixture_ethnicity [CCSS_org_hrc$samples %in% ccss_eur] <- "European"
CCSS_org_hrc$Admixture_ethnicity [CCSS_org_hrc$samples %in% ccss_afr] <- "African"
CCSS_org_hrc$Admixture_ethnicity [CCSS_org_hrc$samples %in% ccss_asa] <- "Asian"

write.table(CCSS_org_hrc, "CCSS_org_hrc_ethnicity.txt", col.names = T, quote = F)


