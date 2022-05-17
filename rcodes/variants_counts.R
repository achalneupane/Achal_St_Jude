######################
## Achal Neupane    ##
## Date: 05/12/2022 ##
######################

## Read variant lists from Zhaoming et al 
zhaoming.etal.vars <- read.table("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Yadav_Sapkota/additional_papers/Zhaoming_Wang_Genetic_risk_for_subsequent_neoplasms_JCO.2018.77.8589/sjlife-genetics-sn.txt", sep = "\t", header =T)
head(zhaoming.etal.vars)
dim(zhaoming.etal.vars)
var.classes <- table(zhaoming.etal.vars$Class)
var.classes
# FRAMESHIFT      MISSENSE      NONSENSE      promoter    PROTEINDEL    PROTEININS        SILENT        SPLICE SPLICE_REGION 
# 22          2269            25             2            83            24            47           100           557 

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/SJLIFE_WGS/common/sjlife/MERGED_SJLIFE_1_2/annotation/SNPEFF_ANNOTATION/")
## read annotated SJLIFE annotated VCF 
library(data.table)
VCF <- fread("MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr17.PASS.decomposed.vcf-annot-snpeff-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155-FIELDS-simple.txt")
table(VCF$CLNSIG)

WANTED.types <- c("$Pathogenic/Likely_pathogenic$|^Likely_pathogenic$|^Pathogenic/Likely_pathogenic$|^Pathogenic$")
sum(grepl(WANTED.types, VCF$CLNSIG, ignore.case = T))
# 1520

WANTED.vars.clinvar <- VCF[grepl(WANTED.types, VCF$CLNSIG, ignore.case = T),]

