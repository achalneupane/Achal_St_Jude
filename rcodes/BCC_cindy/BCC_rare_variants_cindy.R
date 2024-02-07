## read WES annotation files for POI option 1
rm(list=ls())
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/")

gene_regions <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/all_genes_gene_regions.txt", header = T, sep = "\t")

## 1. clinvar
clinvar <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/all_new_clinvar_P_LP.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(clinvar)
clinvar$SNP <- sub(";.*", "", clinvar$ID)
clinvar <- clinvar[!duplicated(clinvar$SNP),]
clinvar$AF <- as.numeric(clinvar$AF)
clinvar$AF_nfe <- as.numeric(clinvar$AF_nfe)
clinvar$AF_afr <- as.numeric(clinvar$AF_afr)


# make rare; .all is based on allee frequency from gnomAD global population; .eur is gnomAD EUR_nfe; .afr is gnomAD AFR
clinvar.all <- clinvar[which(clinvar$AF <0.01),]
clinvar.eur <- clinvar[which(clinvar$AF <0.01 & clinvar$AF_nfe < 0.01),]
clinvar.afr <- clinvar[which(clinvar$AF <0.01 & clinvar$AF_afr < 0.01),]


loftee <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/loftee/loftee_HC_all_chr_with_gnomad.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(loftee)
loftee$SNP <- sub(";.*", "", loftee$Uploaded_variation)
# LOF variants flagged by LOFTEE as dubious (e.g., affecting poorly conserved exons and splice variants affecting NAGNAG sites or non-canonical splice regions) will be excluded.
loftee <- loftee[!grepl("NAGNAG_SITE|NON_CAN_SPLICE", loftee$LoF_flags),]
table(loftee$LoF_flags)
# - 
# 70494 

loftee <- loftee[!duplicated(loftee$SNP),]


loftee$AF <- as.numeric(loftee$AF.1)
loftee$AF_nfe <- as.numeric(loftee$AF_nfe)
loftee$AF_afr <- as.numeric(loftee$AF_afr)


# make rare
loftee.all <- loftee[which(loftee$AF <0.01),]
loftee.eur <- loftee[which(loftee$AF <0.01 & loftee$AF_nfe < 0.01),]
loftee.afr <- loftee[which(loftee$AF <0.01 & loftee$AF_afr < 0.01),]


snpeff <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/missense_variants_with_overlap_in_more_than_90_percent_of_prediction_tools_all_cols.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(snpeff)
snpeff <- snpeff[!duplicated(snpeff$SNP),]

snpeff$AF <- as.numeric(snpeff$AF)
snpeff$AF_nfe <- as.numeric(snpeff$AF_nfe)
snpeff$AF_afr <- as.numeric(snpeff$AF_afr)


# make rare
snpeff.all <- snpeff[which(snpeff$AF <0.01),]
snpeff.eur <- snpeff[which(snpeff$AF <0.01 & snpeff$AF_nfe < 0.01),]
snpeff.afr <- snpeff[which(snpeff$AF <0.01 & snpeff$AF_afr < 0.01),]




length(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))
# 46076
cc <- as.data.frame(unique(c(clinvar$SNP, loftee$SNP, snpeff$SNP)))
dim(cc)
write.table(cc, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/rare_variants_to_extract.txt", row.names = F, col.names = F, quote = F)
# Dear Achal,
#
# I am writing to ask you for your help. Do you think you might be able to help
# us put the following together over the next 2 weeks to meet the mid-Feb ISLCCC
# abstract deadline?
#
# Basically, we would like to get P/LP variant carrier variables in SJLIFE
# (WES/WGS), CCSS Expansion (WES/WGS), and CCSS Original (WES) considering the
# following 5 gene lists below, using any of the 3 rare variant masks that were
# described in our POI concept proposal (SnpEff; LOFTEE; ClinVar - I think if we
# have to pick one, we should pick the strictest definition, which is probably
# ClinVar). Since this is a quick preliminary pass, maybe we could focus on
# getting at least one cancer susceptibility P/LP variant carrier status
# variable ([1] or [2] below) and at least one BCC-related variable ([3]-[5]
# below). Do you think this would be feasible? Please let me know if you have
# questions or if you would like to meet to discuss.
#
# Thank you so much, Cindy
#
# Gene lists
#
# (1) 60 cancer susceptibility genes - should be Zhaoming's SJCPG60 list. See
# Kim_ST1.txt, use "CSG_60" (genes marked with "x"). From Kim et al evaluating
# cancer susceptibility gene P/LP variants in CCSS, link:
# https://doi.org/10.1093/jncics/pkab007.

kim_ST1 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/rare_variants/Kim_ST1.txt", header = T, sep = "\t")
kim_ST1.csg60 <- kim_ST1[grepl("x", kim_ST1$CSG_60),]

kim_ST1.csg60 <- gene_regions[gene_regions$external_gene_name %in% kim_ST1.csg60$Gene.Name,]


kim_ST1.csg60.clinvar.all <- clinvar.all[na.omit(match(kim_ST1.csg60$Gene.Name, clinvar.all$ANN....GENE)),]
kim_ST1.csg60.clinvar.eur <- clinvar.eur[na.omit(match(kim_ST1.csg60$Gene.Name, clinvar.eur$ANN....GENE)),]
kim_ST1.csg60.clinvar.afr <- clinvar.afr[na.omit(match(kim_ST1.csg60$Gene.Name, clinvar.afr$ANN....GENE)),]


kim_ST1.csg60.loftee.all <- loftee.all[na.omit(match(kim_ST1.csg60$Gene.Name, loftee.all$SYMBOL)),]
kim_ST1.csg60.loftee.eur <- loftee.eur[na.omit(match(kim_ST1.csg60$Gene.Name, loftee.eur$SYMBOL)),]
kim_ST1.csg60.loftee.afr <- loftee.afr[na.omit(match(kim_ST1.csg60$Gene.Name, loftee.afr$SYMBOL)),]


kim_ST1.csg60.snpeff.all <- snpeff.all[na.omit(match(kim_ST1.csg60$Gene.Name, snpeff.all$ANN....GENE)),]
kim_ST1.csg60.snpeff.eur <- snpeff.eur[na.omit(match(kim_ST1.csg60$Gene.Name, snpeff.eur$ANN....GENE)),]
kim_ST1.csg60.snpeff.afr <- snpeff.afr[na.omit(match(kim_ST1.csg60$Gene.Name, snpeff.afr$ANN....GENE)),]


# (2) Expanded list of 172 cancer susceptibility genes. Kim_ST1.txt, use
# "CSG_172" (genes marked with "x"; 172 genes evaluated in CCSS)

kim_ST1.csg172 <- kim_ST1[grepl("x", kim_ST1$CSG_172),]


kim_ST1.csg172.clinvar.all <- clinvar.all[na.omit(match(kim_ST1.csg172$Gene.Name, clinvar.all$ANN....GENE)),]
kim_ST1.csg172.clinvar.eur <- clinvar.eur[na.omit(match(kim_ST1.csg172$Gene.Name, clinvar.eur$ANN....GENE)),]
kim_ST1.csg172.clinvar.afr <- clinvar.afr[na.omit(match(kim_ST1.csg172$Gene.Name, clinvar.afr$ANN....GENE)),]


kim_ST1.csg172.loftee.all <- loftee.all[na.omit(match(kim_ST1.csg172$Gene.Name, loftee.all$SYMBOL)),]
kim_ST1.csg172.loftee.eur <- loftee.eur[na.omit(match(kim_ST1.csg172$Gene.Name, loftee.eur$SYMBOL)),]
kim_ST1.csg172.loftee.afr <- loftee.afr[na.omit(match(kim_ST1.csg172$Gene.Name, loftee.afr$SYMBOL)),]


kim_ST1.csg172.snpeff.all <- snpeff.all[na.omit(match(kim_ST1.csg172$Gene.Name, snpeff.all$ANN....GENE)),]
kim_ST1.csg172.snpeff.eur <- snpeff.eur[na.omit(match(kim_ST1.csg172$Gene.Name, snpeff.eur$ANN....GENE)),]
kim_ST1.csg172.snpeff.afr <- snpeff.afr[na.omit(match(kim_ST1.csg172$Gene.Name, snpeff.afr$ANN....GENE)),]


# (3) Literature-based list of BCC genes. See "NCI table 3" for a list of ~20
# genes (Table 3 in the cancer.gov). These genes are based on what has been
# reported in the literature re: BCC-associated syndromes. Based on Kilgour et
# al, Choquet et al, and NCI cancer.gov (all 3 are very consistent):
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8345475/#B6-cancers-13-03870
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7259534/
# https://www.cancer.gov/types/skin/hp/skin-genetics-pdq

NCI_table3 <- unique(c("PTCH1","PTCH2","SUFU", "CYLD", "XPA","XPB", "ERCC3","XPC","XPD", "ERCC2", "XPE", "DDB2", "XPF", "ERCC4", "XPG", "ERCC5", "POLH"))


NCI_table3.clinvar.all <- clinvar.all[na.omit(match(NCI_table3, clinvar.all$ANN....GENE)),]
NCI_table3.clinvar.eur <- clinvar.eur[na.omit(match(NCI_table3, clinvar.eur$ANN....GENE)),]
NCI_table3.clinvar.afr <- clinvar.afr[na.omit(match(NCI_table3, clinvar.afr$ANN....GENE)),]


NCI_table3.loftee.all <- loftee.all[na.omit(match(NCI_table3, loftee.all$SYMBOL)),]
NCI_table3.loftee.eur <- loftee.eur[na.omit(match(NCI_table3, loftee.eur$SYMBOL)),]
NCI_table3.loftee.afr <- loftee.afr[na.omit(match(NCI_table3, loftee.afr$SYMBOL)),]


NCI_table3.snpeff.all <- snpeff.all[na.omit(match(NCI_table3, snpeff.all$ANN....GENE)),]
NCI_table3.snpeff.eur <- snpeff.eur[na.omit(match(NCI_table3, snpeff.eur$ANN....GENE)),]
NCI_table3.snpeff.afr <- snpeff.afr[na.omit(match(NCI_table3, snpeff.afr$ANN....GENE)),]


# (4) BCC gene panel genes. See "BCC_blueprint_panel.txt" for a list of genes
# tested by Blueprint Genetics panel for BCC:
# https://blueprintgenetics.com/tests/panels/dermatology/hereditary-melanoma-and-skin-cancer-panel/

BCC.panel <- c("BAP1","BRCA1","BRCA2","CDK4","CDKN2A","DDB2","ERCC2","ERCC3","ERCC4","ERCC5","MITF","POT1","PTCH1","PTEN","SUFU","TP53","WRN","XPA","XPC")

BCC.panel.clinvar.all <- clinvar.all[na.omit(match(BCC.panel, clinvar.all$ANN....GENE)),]
BCC.panel.clinvar.eur <- clinvar.eur[na.omit(match(BCC.panel, clinvar.eur$ANN....GENE)),]
BCC.panel.clinvar.afr <- clinvar.afr[na.omit(match(BCC.panel, clinvar.afr$ANN....GENE)),]


BCC.panel.loftee.all <- loftee.all[na.omit(match(BCC.panel, loftee.all$SYMBOL)),]
BCC.panel.loftee.eur <- loftee.eur[na.omit(match(BCC.panel, loftee.eur$SYMBOL)),]
BCC.panel.loftee.afr <- loftee.afr[na.omit(match(BCC.panel, loftee.afr$SYMBOL)),]


BCC.panel.snpeff.all <- snpeff.all[na.omit(match(BCC.panel, snpeff.all$ANN....GENE)),]
BCC.panel.snpeff.eur <- snpeff.eur[na.omit(match(BCC.panel, snpeff.eur$ANN....GENE)),]
BCC.panel.snpeff.afr <- snpeff.afr[na.omit(match(BCC.panel, snpeff.afr$ANN....GENE)),]


# (5) ClinVar BCC-related genes. See "clinvar_result_BCC.txt". These are all
# ClinVar listings for "basal AND cell AND carcinoma" annotated as P or LP and
# with multiple submitters/no conflicts. This list hasn't been reviewed yet in
# terms of conditions, but I think we can use the variant/gene list. For the
# final analysis, I can ask the clinicians to help review the conditions to make
# sure we don't have anything unrelated.

clinvar.BCC.result <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/rare_variants/clinvar_result_BCC.txt", sep = "\t", header = T)
paste0(unique(clinvar.BCC.result$Gene.s.), collapse = "',")
clinvar.BCC.result <- unique(c('ERCC3','XPC','LOC129936244','XPC','BAP1','CCNH','RASA1',
                        'POLH','POLR1C','PTCH1','LOC100507346','PTCH1','LOC130002133',
                        'PTCH1','XPA','RET','LOC130004614','SUFU','SUFU','HRAS','LRRC56',
                        'DDB2','ATM','ATM','C11orf65','BIVM-ERCC5','ERCC5','LOC126861834',
                        'BIVM-ERCC5','ERCC5','BIVM-ERCC5','ERCC5','LOC126861834','NTHL1',
                        'ERCC4','PALB2','TP53','ERCC2','MT-TL1'))


clinvar.BCC.result.clinvar.all <- clinvar.all[na.omit(match(clinvar.BCC.result, clinvar.all$ANN....GENE)),]
clinvar.BCC.result.clinvar.eur <- clinvar.eur[na.omit(match(clinvar.BCC.result, clinvar.eur$ANN....GENE)),]
clinvar.BCC.result.clinvar.afr <- clinvar.afr[na.omit(match(clinvar.BCC.result, clinvar.afr$ANN....GENE)),]


clinvar.BCC.result.loftee.all <- loftee.all[na.omit(match(clinvar.BCC.result, loftee.all$SYMBOL)),]
clinvar.BCC.result.loftee.eur <- loftee.eur[na.omit(match(clinvar.BCC.result, loftee.eur$SYMBOL)),]
clinvar.BCC.result.loftee.afr <- loftee.afr[na.omit(match(clinvar.BCC.result, loftee.afr$SYMBOL)),]


clinvar.BCC.result.snpeff.all <- snpeff.all[na.omit(match(clinvar.BCC.result, snpeff.all$ANN....GENE)),]
clinvar.BCC.result.snpeff.eur <- snpeff.eur[na.omit(match(clinvar.BCC.result, snpeff.eur$ANN....GENE)),]
clinvar.BCC.result.snpeff.afr <- snpeff.afr[na.omit(match(clinvar.BCC.result, snpeff.afr$ANN....GENE)),]


########################
## Get carrier status ##
########################
library("data.table")
raw <- fread("all_BCC_rare_variants_recodeA.raw", header = T)
raw <- as.data.frame(raw)
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
# HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
colnames(raw) = HEADER

raw <- raw[, !(colnames(raw) %in% c("PAT", "MAT", "SEX", "PHENOTYPE"))]

# colnames(raw) # check this with the .bim REF and NON-REFERENCE alleles
dim(raw)
# 46078

sjl <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//Survivor_WES/biallelic/extract_SJLIFE_survivor_iid_fid.txt", header = F)
ccss_exp <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common//Survivor_WES/biallelic/extract_CCSS.samples_iid_fid.txt", header = F)
raw$cohort[(raw$FID %in% sjl$V1)] <- "SJLIFE"
raw$cohort[(raw$FID %in% ccss_exp$V1)] <- "CCSS_exp"

###########
## CSG60 ##
###########
## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% kim_ST1.csg60.clinvar.all$SNP]
raw$kim_ST1.csg60.clinvar.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% kim_ST1.csg60.loftee.all$SNP]
raw$kim_ST1.csg60.loftee.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% kim_ST1.csg60.loftee.all$SNP]
raw$kim_ST1.csg60.loftee.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% kim_ST1.csg60.snpeff.all$SNP]
raw$kim_ST1.csg60.snpeff.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

############
## CSG172 ##
############

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% kim_ST1.csg172.clinvar.all$SNP]
raw$kim_ST1.csg172.clinvar.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% kim_ST1.csg172.loftee.all$SNP]
raw$kim_ST1.csg172.loftee.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% kim_ST1.csg172.loftee.all$SNP]
raw$kim_ST1.csg172.loftee.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% kim_ST1.csg172.snpeff.all$SNP]
raw$kim_ST1.csg172.snpeff.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")


################
## NCI_table3 ##
################

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% NCI_table3.clinvar.all$SNP]
raw$NCI_table3.clinvar.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% NCI_table3.loftee.all$SNP]
raw$NCI_table3.loftee.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% NCI_table3.loftee.all$SNP]
raw$NCI_table3.loftee.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% NCI_table3.snpeff.all$SNP]
raw$NCI_table3.snpeff.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

###############
## BCC.panel ##
###############

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% BCC.panel.clinvar.all$SNP]
raw$BCC.panel.clinvar.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% BCC.panel.loftee.all$SNP]
raw$BCC.panel.loftee.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% BCC.panel.loftee.all$SNP]
raw$BCC.panel.loftee.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% BCC.panel.snpeff.all$SNP]
raw$BCC.panel.snpeff.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

########################
## clinvar.BCC.result ##
########################

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% clinvar.BCC.result.clinvar.all$SNP]
raw$clinvar.BCC.result.clinvar.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% clinvar.BCC.result.loftee.all$SNP]
raw$clinvar.BCC.result.loftee.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% clinvar.BCC.result.loftee.all$SNP]
raw$clinvar.BCC.result.loftee.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

## Extract variant carrier status
col.extract <- colnames(raw)[colnames(raw) %in% clinvar.BCC.result.snpeff.all$SNP]
raw$clinvar.BCC.result.snpeff.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

carrier.status <- raw[,!grepl("chr", colnames(raw))]

carrier.status <- carrier.status[!is.na(carrier.status$cohort),]

##############################
## Add PRS scores to SJLIFE ##
##############################

ST6_prs <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/sjlife/prs_out/ST6_prs.profile", header = T)
PGS000356_prs <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/sjlife/prs_out/PGS000356_prs.profile", header = T)
PGS000454_prs <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/sjlife/prs_out/PGS000454_prs.profile", header = T)
PGS003416_prs <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/sjlife/prs_out/PGS003416_prs.profile", header = T)

SJLIFE.PRS <- cbind.data.frame(IID = ST6_prs$IID, ST6_prs=ST6_prs$SCORE, PGS000356_prs=PGS000356_prs$SCORE, PGS000454_prs=PGS000454_prs$SCORE, PGS003416_prs=PGS003416_prs$SCORE)


################################
## Add PRS scores to ccss_exp ##
################################

ST6_prs <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/ccss_exp/prs_out/ST6_prs.profile", header = T)
PGS000356_prs <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/ccss_exp/prs_out/PGS000356_prs.profile", header = T)
PGS000454_prs <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/ccss_exp/prs_out/PGS000454_prs.profile", header = T)
PGS003416_prs <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/ccss_exp/prs_out/PGS003416_prs.profile", header = T)

CCSS_exp.PRS <- cbind.data.frame(IID = ST6_prs$IID, ST6_prs=ST6_prs$SCORE, PGS000356_prs=PGS000356_prs$SCORE, PGS000454_prs=PGS000454_prs$SCORE, PGS003416_prs=PGS003416_prs$SCORE)


################################
## Add PRS scores to ccss_org ##
################################

ST6_prs <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/ccss_org/prs_out/ST6_prs.profile", header = T)
PGS000356_prs <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/ccss_org/prs_out/PGS000356_prs.profile", header = T)
PGS000454_prs <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/ccss_org/prs_out/PGS000454_prs.profile", header = T)
PGS003416_prs <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/prs/ccss_org/prs_out/PGS003416_prs.profile", header = T)

CCSS_org.PRS <- cbind.data.frame(IID = ST6_prs$IID, ST6_prs=ST6_prs$SCORE, PGS000356_prs=PGS000356_prs$SCORE, PGS000454_prs=PGS000454_prs$SCORE, PGS003416_prs=PGS003416_prs$SCORE)

save(list = c("carrier.status", "SJLIFE.PRS", "CCSS_exp.PRS", "CCSS_org.PRS"), file = "BCC_carrier_and_PRS.RData")





## Check with start position for all genes
all.genes <- unique(c(kim_ST1.csg60$Gene.Name, kim_ST1.csg172$Gene.Name, NCI_table3, BCC.panel, clinvar.BCC.result))

write.table(all.genes, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/all_genes_to_extract.txt", row.names = F, col.names = F, quote = F)                    


###########################################################################################################################

# Hi Achal,
# 
# Attached are the FIDs+IIDs I used for our analysis (recall I arbitrarily removed duplicates; see carrier.status) and the EUR studyids involved in the analysis (subdf).
# 
# Would you have time to help with the following by this Friday noon? 
#   
#   I was hoping I could get an overall frequency count of the genes with the most P/LP variants in our data, using Clinvar, LOFTEE, and Clinvar or LOFTEE as the masks of interest. Please see the figure below from my grant with Zhaoming using CSG60 as an example. We don't need this pretty plot, just tables with the total P/LP counts for maybe the top 20 genes in csg172 and all of the "panel" genes.
# 
# Please let me know if this might be possible. Totally okay if you cannot manage - sorry to ask you on such a tight turnaround.
# 
# Thank you,
# Cindy

# Thank you, Achal! Also, before I forget, please provide the total number of genes in csg172 and panel sets with P/LP variants in our data.
# length(unique(kim_ST1.csg172.clinvar.all$ANN....GENE)) = 84 genes
# length(unique(kim_ST1.csg172.loftee.all$SYMBOL)) = 74

# length(unique(BCC.panel.clinvar.all$ANN....GENE)) = 16
# length(unique(BCC.panel.loftee.all$SYMBOL)) = 12 
###################################
## 1. kim_ST1.csg172.clinvar.all ##
###################################
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/")
load("Achal_carrier.status.RData")

raw.2 <- raw[raw$IID %in% subdf$studyid,]
raw.2 <- raw.2[raw.2$FID %in% carrier.status$FID,]

library(dplyr)
## Make genotype 2 to 1, so we can get the frequency
raw.2 <- raw.2 %>%
  mutate_at(vars(-c(1, 2)), ~ ifelse(. == 2, 1, .))

row.names(raw.2) <- raw.2$IID

kim_ST1.csg172.clinvar.all.raw.2 <- raw.2[kim_ST1.csg172.clinvar.all$SNP]


################################

































































###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

# # install.packages("biomaRt")
# library(biomaRt)
# # Install and load the biomaRt package
# install.packages("biomaRt")
# library(biomaRt)
# 
# # Use the ensembl dataset for human genes (GRCh38)
# ensembl <- useMart("ensembl", host = "www.ensembl.org", dataset = "hsapiens_gene_ensembl")
# # Specify the genes of interest
# genes_of_interest <- all.genes
# 
# # Get gene regions for the specified genes
# gene_regions <- getBM(attributes = c("external_gene_name", "chromosome_name", "start_position", "end_position"),
#                       filters = "external_gene_name",
#                       values = genes_of_interest,
#                       mart = ensembl)
# 
# # Print the results
# print(gene_regions)
# gene_regions <- gene_regions[!grepl("[[:alpha:]]", gene_regions$chromosome_name), ]
# 
# 
# gene_regions
# 
# # cc <- as.data.frame(listAttributes(ensembl))
# 
# all.genes[!all.genes %in% gene_regions$external_gene_name]
# # all.genes[!all.genes %in% gene_regions$external_gene_name]
# # [1] "T"            "XPB"          "XPD"          "XPE"          "XPF"          "XPG"          "LOC129936244"
# # [8] "LOC100507346" "LOC130002133" "LOC130004614" "LOC126861834" "MT-TL1" 
# 
# write.table(gene_regions, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/all_genes_gene_regions.txt", row.names = F, col.names = F, quote = F, sep = "\t")


