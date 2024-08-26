## read WES annotation files for POI option 1
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/")


gene_regions <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/all_genes_gene_regions.txt", header = T, sep = "\t")

## 1. clinvar
clinvar <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/Survivor_WES/annotation/snpEff_round3/clinvar_version_check/clinvar_04_17_2019//all_new_clinvar_P_LP.txt", header = T, sep = "\t", stringsAsFactors = F)
dim(clinvar)
clinvar$SNP <- sub(";.*", "", clinvar$ID)
clinvar <- clinvar[!duplicated(clinvar$SNP),]
clinvar$AF <- as.numeric(clinvar$AF)
clinvar$AF_nfe <- as.numeric(clinvar$AF_nfe)
clinvar$AF_afr <- as.numeric(clinvar$AF_afr)


# make rare; .all is based on allee frequency from gnomAD global population; .eur is gnomAD EUR_nfe; .afr is gnomAD AFR
clinvar.all <- clinvar[which(clinvar$AF <0.01),]
length(unique(clinvar.all$SNP))
# 5223
# clinvar$dbNSFP_ExAC_AF <- as.numeric(clinvar$dbNSFP_ExAC_AF)
# clinvar.all <- clinvar[which(clinvar$dbNSFP_ExAC_AF <0.01),]
cc <- cbind.data.frame(clinvar.all$SNP, clinvar.all$AF, clinvar.all$AF_nfe, clinvar.all$AF_afr, clinvar.all$dbNSFP_ExAC_AF, clinvar.all$CLNSIG)

clinvar.eur <- clinvar[which(clinvar$AF <0.01 & clinvar$AF_nfe < 0.01),]
clinvar.afr <- clinvar[which(clinvar$AF <0.01 & clinvar$AF_afr < 0.01),]

cc <- as.data.frame(unique(c(clinvar.all$SNP)))
dim(cc)
# 4882
write.table(cc, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/clinvar_version_check/04_17_2019/rare_variants_to_extract.txt", row.names = F, col.names = F, quote = F)
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

## Make a function to extract
extract_variants <- function(input_df, reference_df) {
  result <- data.frame()
  
  for (i in seq(nrow(input_df))) {
    # Extract rows from reference_df where chromosome matches and position is within the range
    subset_result <- reference_df[reference_df$CHROM == paste0("chr", input_df$chromosome_name[i]) & 
                                    reference_df$POS >= input_df$start_position[i] & 
                                    reference_df$POS <= input_df$end_position[i], ]
    
    # Append the subset_result to the overall result
    result <- rbind(result, subset_result)
  }
  
  return(result)
}

# Example usage
# result <- extract_variants(kim_ST1.csg60, clinvar.all)




#
# (1) 60 cancer susceptibility genes - should be Zhaoming's SJCPG60 list. See
# Kim_ST1.txt, use "CSG_60" (genes marked with "x"). From Kim et al evaluating
# cancer susceptibility gene P/LP variants in CCSS, link:
# https://doi.org/10.1093/jncics/pkab007.

kim_ST1 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/rare_variants/Kim_ST1.txt", header = T, sep = "\t")
kim_ST1.csg60 <- kim_ST1[grepl("x", kim_ST1$CSG_60),]

library(readxl)
kim_data <- read_excel("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/rare_variants/kim_et_al/Copy of CS20-0146R1 Kim Supp tables1-3_010721.xlsx", sheet = "S2_AN")
# kim_data.clinvar <- kim_data[kim_data$`Method of classification`== "ClinVar" ,]
kim_data.csg60 <- kim_data[kim_data$`Gene Name` %in% kim_ST1.csg60$Gene.Name,]
length(unique(kim_data.csg60$`Gene Name`))

# kim_ST1.csg60 <- gene_regions[gene_regions$external_gene_name %in% kim_ST1.csg60$Gene.Name,]
# kim_ST1.csg60.clinvar.all <- extract_variants(kim_ST1.csg60, clinvar.all)

kim_ST1.csg60.clinvar.all <- clinvar.all[clinvar.all$ANN....GENE %in% kim_ST1.csg60$Gene.Name,]
dim(kim_ST1.csg60.clinvar.all)
# [1] 285 218
length(unique(kim_ST1.csg60.clinvar.all$ANN....GENE))
# 43
kim_ST1.csg60.clinvar.all$KEY <- paste0(kim_ST1.csg60.clinvar.all$CHROM, ":", kim_ST1.csg60.clinvar.all$POS)
# check if they are in CSG60 vars lists
csg.60.vars <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/rare_variants/CSG60_zhaoming/CSG_60_variants_edited.txt", header = T, sep = "\t")

## Keep those in 60 list
csg.60.vars <- csg.60.vars[csg.60.vars$Gene %in% kim_ST1.csg60$Gene.Name,]

csg.60.vars$KEY <- paste0("chr",csg.60.vars$CHROM, ":", csg.60.vars$POS_GRCh38)

sum(csg.60.vars$KEY %in% kim_ST1.csg60.clinvar.all$KEY)
## 79/ 166
length(unique(csg.60.vars$Gene))
## 31
length(unique(kim_ST1.csg60.clinvar.all$SNP))
# 285
length(unique(kim_ST1.csg60.clinvar.all$ANN....GENE))
## 43 genes
sum(unique(csg.60.vars$Gene) %in% unique(kim_ST1.csg60.clinvar.all$ANN....GENE))
## 29

csg.60.vars$duplicated <- ifelse(duplicated(csg.60.vars$KEY), "yes", "no")
csg.60.vars$SNP <- paste0(csg.60.vars$KEY, ":", csg.60.vars$REF, ":", csg.60.vars$ALT)

length(unique(csg.60.vars$SNP))
## 148
csg.60.vars.unique <- csg.60.vars[!duplicated(csg.60.vars$SNP),]
csg.60.vars.unique$matched <- kim_ST1.csg60.clinvar.all$SNP[match(csg.60.vars.unique$KEY, kim_ST1.csg60.clinvar.all$KEY)]
table(is.na(csg.60.vars.unique$matched))
# FALSE  TRUE 
# 67    81 

## 43 0f 60 genes and total of 285 variants found in WES based on clinvar. 29/31 genes match based on clinvar annotation Of these, 67 variants match out of 148 variants in Kim (clinvar set) 
# write.table(csg.60.vars.unique, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/csg.60.vars.unique.txt", row.names = F, col.names = T, quote = F, sep ="\t")

kim_ST1.csg60.loftee.all <- loftee.all[loftee.all$SYMBOL %in% kim_ST1.csg60$Gene.Name,]

kim_ST1.csg60.snpeff.all <- snpeff.all[snpeff.all$ANN....GENE %in% kim_ST1.csg60$Gene.Name,]


# (2) Expanded list of 172 cancer susceptibility genes. Kim_ST1.txt, use
# "CSG_172" (genes marked with "x"; 172 genes evaluated in CCSS)

## get 172 gene list
kim_ST1.csg172 <- kim_ST1[grepl("x", kim_ST1$CSG_172),]


kim_ST1.csg172.clinvar.all <- clinvar.all[clinvar.all$ANN....GENE %in% kim_ST1.csg172$Gene.Name,]
dim(kim_ST1.csg172.clinvar.all) ## These are extracted variants from matched genes
# [1] 473 218

length(unique(kim_ST1.csg172$Gene.Name[!kim_ST1.csg172$Gene.Name %in% kim_ST1.csg172.clinvar.all$ANN....GENE]))
# 88 of 172 genes match

kim_ST1.csg172.clinvar.all$KEY <- paste0(kim_ST1.csg172.clinvar.all$CHROM, ":", kim_ST1.csg172.clinvar.all$POS)
## check if they are in CSG172 vars lists
kim_ST1.csg172.vars <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/rare_variants/kim_et_al/CS20-0146R1 Kim Supp tables1-3_010721.txt", header = T, sep = "\t")
kim_ST1.csg172.vars <- kim_ST1.csg172.vars[grepl("Clinvar", kim_ST1.csg172.vars$Method.of.classification, ignore.case = T),]
grch38 <- read.delim("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/rare_variants/kim_et_al/hglft_genome_d359_661a0.bed", header = F, sep = "\t")

kim_ST1.csg172.vars$KEY <- paste0(kim_ST1.csg172.vars$Chr, ":", kim_ST1.csg172.vars$Position)
grch38$POS37KEY <- sub("-.*", "", grch38$V4)
## GRCh38
kim_ST1.csg172.vars$Position <- grch38$V3[match(kim_ST1.csg172.vars$KEY, grch38$POS37KEY)]
## Keep those in 172 list
kim_ST1.csg172.vars <- kim_ST1.csg172.vars[kim_ST1.csg172.vars$Gene.Name %in% kim_ST1.csg172$Gene.Name,]

kim_ST1.csg172.vars$KEY <- paste0(kim_ST1.csg172.vars$Chr, ":", kim_ST1.csg172.vars$Position)

sum(kim_ST1.csg172.vars$KEY %in% kim_ST1.csg172.clinvar.all$KEY)
## 160/ 383
length(unique(kim_ST1.csg172.vars$Gene.Name))
## 46
length(unique(kim_ST1.csg172.clinvar.all$SNP))
# 473
length(unique(kim_ST1.csg172.clinvar.all$ANN....GENE))
## 84 genes
sum(unique(kim_ST1.csg172.vars$Gene.Name) %in% unique(kim_ST1.csg172.clinvar.all$ANN....GENE))
## 41

kim_ST1.csg172.vars$duplicated <- ifelse(duplicated(kim_ST1.csg172.vars$KEY), "yes", "no")
kim_ST1.csg172.vars$SNP <- paste0(kim_ST1.csg172.vars$KEY, ":", kim_ST1.csg172.vars$Ref, ":", kim_ST1.csg172.vars$Alt)

length(unique(kim_ST1.csg172.vars$SNP))
## 188
kim_ST1.csg172.vars.unique <- kim_ST1.csg172.vars[!duplicated(kim_ST1.csg172.vars$SNP),]
kim_ST1.csg172.vars.unique$matched <- kim_ST1.csg172.clinvar.all$SNP[match(kim_ST1.csg172.vars.unique$KEY, kim_ST1.csg172.clinvar.all$KEY)]
table(is.na(kim_ST1.csg172.vars.unique$matched))
# FALSE  TRUE 
# 64   124 
## 84 0f 172 genes and total of 473 variants found in WES based on clinvar. 41/46 genes match based on clinvar annotation Of these, 64 variants match out of 188 variants in Kim (clinvar set) 

# cc <- kim_ST1.csg172.clinvar.all
# cc$matched <- kim_ST1.csg172.vars.unique$SNP[match(cc$KEY, kim_ST1.csg172.vars.unique$KEY)]
# table(is.na(cc$matched))
# gg <- cbind.data.frame(cc$SNP, cc$matched)

# write.table(gg, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/found_inSJLIFE_but_not_in_kim_ST1.csg172.vars.unique.txt", row.names = F, col.names = T, quote = F, sep ="\t")
# write.table(kim_ST1.csg172.vars.unique, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/kim_ST1.csg172.vars.unique.txt", row.names = F, col.names = T, quote = F, sep ="\t")

kim_ST1.csg172.loftee.all <- loftee.all[loftee.all$SYMBOL %in% kim_ST1.csg172$Gene.Name,]

kim_ST1.csg172.snpeff.all <- snpeff.all[snpeff.all$ANN....GENE %in% kim_ST1.csg172$Gene.Name,]


# (3) Literature-based list of BCC genes. See "NCI table 3" for a list of ~20
# genes (Table 3 in the cancer.gov). These genes are based on what has been
# reported in the literature re: BCC-associated syndromes. Based on Kilgour et
# al, Choquet et al, and NCI cancer.gov (all 3 are very consistent):
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8345475/#B6-cancers-13-03870
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7259534/
# https://www.cancer.gov/types/skin/hp/skin-genetics-pdq

NCI_table3 <- unique(c("PTCH1","PTCH2","SUFU", "CYLD", "XPA","XPB", "ERCC3","XPC","XPD", "ERCC2", "XPE", "DDB2", "XPF", "ERCC4", "XPG", "ERCC5", "POLH"))


NCI_table3.clinvar.all <- clinvar.all[clinvar.all$ANN....GENE %in% NCI_table3,]


NCI_table3.loftee.all <- loftee.all[loftee.all$SYMBOL %in% NCI_table3,]


NCI_table3.snpeff.all <- snpeff.all[snpeff.all$ANN....GENE %in% NCI_table3,]


# (4) BCC gene panel genes. See "BCC_blueprint_panel.txt" for a list of genes
# tested by Blueprint Genetics panel for BCC:
# https://blueprintgenetics.com/tests/panels/dermatology/hereditary-melanoma-and-skin-cancer-panel/

BCC.panel <- c("BAP1","BRCA1","BRCA2","CDK4","CDKN2A","DDB2","ERCC2","ERCC3","ERCC4","ERCC5","MITF","POT1","PTCH1","PTEN","SUFU","TP53","WRN","XPA","XPC")

BCC.panel.clinvar.all <- clinvar.all[clinvar.all$ANN....GENE %in% BCC.panel,]


BCC.panel.loftee.all <- loftee.all[loftee.all$SYMBOL %in% BCC.panel,]


BCC.panel.snpeff.all <- snpeff.all[snpeff.all$ANN....GENE %in% BCC.panel,]


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


clinvar.BCC.result.clinvar.all <- clinvar.all[clinvar.all$ANN....GENE %in% clinvar.BCC.result,]


clinvar.BCC.result.loftee.all <- loftee.all[loftee.all$SYMBOL %in% clinvar.BCC.result,]


clinvar.BCC.result.snpeff.all <- snpeff.all[snpeff.all$ANN....GENE %in% clinvar.BCC.result,]


########################
## Get carrier status ##
########################
library("data.table")
raw <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/clinvar_version_check/04_17_2019/all_BCC_rare_variants_recodeA.raw", header = T)
raw <- as.data.frame(raw)
HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
# HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
colnames(raw) = HEADER

raw <- raw[, !(colnames(raw) %in% c("PAT", "MAT", "SEX", "PHENOTYPE"))]

# colnames(raw) # check this with the .bim REF and NON-REFERENCE alleles
dim(raw)
# [1]  8065 46078

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
col.extract <- colnames(raw)[colnames(raw) %in% clinvar.BCC.result.snpeff.all$SNP]
raw$clinvar.BCC.result.snpeff.all.status <- ifelse(rowSums(raw[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

carrier.status <- raw[,!grepl("chr", colnames(raw))]

carrier.status <- carrier.status[!is.na(carrier.status$cohort),]


raw.sjlife <- raw[which(raw$cohort == "SJLIFE"),1:2]

carrier.status.sjlife <- carrier.status[grepl("SJLIFE", carrier.status$cohort),]
carrier.status.ccss <- carrier.status[grepl("CCSS", carrier.status$cohort),]
sjlife_1 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/sjlife_1/sample_list.txt", header = F)
carrier.status.sjlife2986 <- carrier.status.sjlife[carrier.status.sjlife$IID %in% sjlife_1$V1,]

## SJLIFE
round(table(carrier.status.sjlife$kim_ST1.csg60.clinvar.all.status)[2]/nrow(carrier.status.sjlife)*100,2)
round(table(carrier.status.sjlife$kim_ST1.csg172.clinvar.all.status)[2]/nrow(carrier.status.sjlife)*100,2)
round(table(carrier.status.sjlife$BCC.panel.clinvar.all.status)[2]/nrow(carrier.status.sjlife)*100,2) 
round(table(carrier.status.sjlife$clinvar.BCC.result.clinvar.all.status)[2]/nrow(carrier.status.sjlife)*100,2)
round(table(carrier.status.sjlife$NCI_table3.clinvar.all.status)[2]/nrow(carrier.status.sjlife)*100,2) 

## CCSS
round(table(carrier.status.ccss$kim_ST1.csg60.clinvar.all.status)[2]/nrow(carrier.status.ccss)*100,2) 
round(table(carrier.status.ccss$kim_ST1.csg172.clinvar.all.status)[2]/nrow(carrier.status.ccss)*100,2)
round(table(carrier.status.ccss$BCC.panel.clinvar.all.status)[2]/nrow(carrier.status.ccss)*100,2) 
round(table(carrier.status.ccss$clinvar.BCC.result.clinvar.all.status)[2]/nrow(carrier.status.ccss)*100,2) 
round(table(carrier.status.ccss$NCI_table3.clinvar.all.status)[2]/nrow(carrier.status.ccss)*100,2) 

## SJLIFE1 (2986 samples)
round(table(carrier.status.sjlife2986$kim_ST1.csg60.clinvar.all.status)[2]/nrow(carrier.status.sjlife2986)*100,2) 
round(table(carrier.status.sjlife2986$kim_ST1.csg172.clinvar.all.status)[2]/nrow(carrier.status.sjlife2986)*100,2)
round(table(carrier.status.sjlife2986$BCC.panel.clinvar.all.status)[2]/nrow(carrier.status.sjlife2986)*100,2) 
round(table(carrier.status.sjlife2986$clinvar.BCC.result.clinvar.all.status)[2]/nrow(carrier.status.sjlife2986)*100,2) 
round(table(carrier.status.sjlife2986$NCI_table3.clinvar.all.status)[2]/nrow(carrier.status.sjlife2986)*100,2) 


# IN SJLIFE, 
# 6.15% for CSG60
# 29.62% in CSG172
# 3.12 % BCC panel
# 3.34% clinvar.BCC.result
# 0.96% in NCI_Table3
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

save(list = c("carrier.status", "SJLIFE.PRS", "CCSS_exp.PRS", "CCSS_org.PRS"), file = "BCC_carrier_and_PRS_v2.RData")




## Check with start position for all genes
all.genes <- unique(c(kim_ST1.csg60$Gene.Name, kim_ST1.csg172$Gene.Name, NCI_table3, BCC.panel, clinvar.BCC.result))

write.table(all.genes, "Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/all_genes_to_extract.txt", row.names = F, col.names = F, quote = F)                    

# save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/BCC_rare_variants_cindy_v2.RData")
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
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/BCC_rare_variants_cindy_v2.RData")
# load("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/BCC_rare_variants_cindy_v2.RData")
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/")
# setwd("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy")
load("Achal_carrier.status.RData")

raw.2 <- raw[raw$IID %in% subdf$studyid,]
raw.2 <- raw.2[raw.2$FID %in% carrier.status$FID,]
raw.2 <- raw.2[grepl("FID|IID|chr", colnames(raw.2))]

rownames(raw.2) <- raw.2$IID
raw.2 <- raw.2[-c(1:2)]



## Get the name of all variants and corresponding genes
clinvar.genes <- cbind.data.frame(SNP=clinvar$SNP, GENE=clinvar$ANN....GENE)
snpeff.genes <- cbind.data.frame(SNP=snpeff$SNP, GENE=snpeff$ANN....GENE)
loftee.genes <- cbind.data.frame(SNP=loftee$SNP, GENE=loftee$SYMBOL)

all.tools.genes <- rbind.data.frame(clinvar.genes, snpeff.genes, loftee.genes)
## remove duplicated rows
all.tools.genes <- all.tools.genes %>% distinct()


library(dplyr)
## Make genotype 2 to 1, so we can get the frequency
# raw.2 <- raw.2 %>%
#   mutate_at(vars(-c(1, 2)), ~ ifelse(. == 2, 1, .)) # except first 2 columns
raw.2 <- raw.2 %>% mutate_all(~ ifelse(. == 2, 1, .)) # all columns
# cc <- raw.2[1:50]
# table(cc$`chr1:976215:A:G`)
raw.2 <- t(raw.2)

#############
## CSG 172 ##
#############
## Clinvar 
kim_ST1.csg172.clinvar.all.raw.2 <- as.data.frame(raw.2[kim_ST1.csg172.clinvar.all$SNP,])
kim_ST1.csg172.clinvar.all.raw.2$Carriers_by_variant <- rowSums(kim_ST1.csg172.clinvar.all.raw.2, na.rm = T)
kim_ST1.csg172.clinvar.all.raw.2$SNP <- row.names(kim_ST1.csg172.clinvar.all.raw.2)
kim_ST1.csg172.clinvar.all.raw.2 <- cbind.data.frame(SNP=kim_ST1.csg172.clinvar.all.raw.2$SNP, Carriers_by_variant=kim_ST1.csg172.clinvar.all.raw.2$Carriers_by_variant)
kim_ST1.csg172.clinvar.all.raw.2$GENE <- all.tools.genes$GENE[match(kim_ST1.csg172.clinvar.all.raw.2$SNP, all.tools.genes$SNP)]
kim_ST1.csg172.clinvar.all.raw.2 <- kim_ST1.csg172.clinvar.all.raw.2 %>%
  group_by(GENE) %>%
  mutate(variants_n = n())

# kim_ST1.csg172.clinvar.all.raw.2 <- kim_ST1.csg172.clinvar.all.raw.2[kim_ST1.csg172.clinvar.all.raw.2$Carriers_by_variant != 0,]

write.table(kim_ST1.csg172.clinvar.all.raw.2, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/kim_ST1.csg172.clinvar.all.raw.2.txt", row.names = F, col.names = T, quote = F, sep = "\t")

kim_ST1.csg172.clinvar.all.raw.2_summary <- kim_ST1.csg172.clinvar.all.raw.2 %>%
  group_by(GENE) %>%
  summarize(
    total_carriers = sum(Carriers_by_variant),
    total_variants = first(variants_n)
  )

write.table(kim_ST1.csg172.clinvar.all.raw.2_summary, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/kim_ST1.csg172.clinvar.all.raw.2_summary.txt", row.names = F, col.names = T, quote = F, sep = "\t")


## Check for discrepancy in prevalence between SJLIFE and 
kim_ST1.csg172.vars




## loftee 
kim_ST1.csg172.loftee.all.raw.2 <- as.data.frame(raw.2[kim_ST1.csg172.loftee.all$SNP,])
kim_ST1.csg172.loftee.all.raw.2$Carriers_by_variant <- rowSums(kim_ST1.csg172.loftee.all.raw.2, na.rm = T)
kim_ST1.csg172.loftee.all.raw.2$SNP <- row.names(kim_ST1.csg172.loftee.all.raw.2)
kim_ST1.csg172.loftee.all.raw.2 <- cbind.data.frame(SNP=kim_ST1.csg172.loftee.all.raw.2$SNP, Carriers_by_variant=kim_ST1.csg172.loftee.all.raw.2$Carriers_by_variant)
kim_ST1.csg172.loftee.all.raw.2$GENE <- all.tools.genes$GENE[match(kim_ST1.csg172.loftee.all.raw.2$SNP, all.tools.genes$SNP)]
kim_ST1.csg172.loftee.all.raw.2 <- kim_ST1.csg172.loftee.all.raw.2 %>%
  group_by(GENE) %>%
  mutate(variants_n = n())

# kim_ST1.csg172.loftee.all.raw.2 <- kim_ST1.csg172.loftee.all.raw.2[-1]

write.table(kim_ST1.csg172.loftee.all.raw.2, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/kim_ST1.csg172.loftee.all.raw.2.txt", row.names = F, col.names = T, quote = F, sep = "\t")

kim_ST1.csg172.loftee.all.raw.2_summary <- kim_ST1.csg172.loftee.all.raw.2 %>%
  group_by(GENE) %>%
  summarize(
    total_carriers = sum(Carriers_by_variant),
    total_variants = first(variants_n)
  )

write.table(kim_ST1.csg172.loftee.all.raw.2_summary, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/kim_ST1.csg172.loftee.all.raw.2_summary.txt", row.names = F, col.names = T, quote = F, sep = "\t")


## loftee or clinvar
kim_ST1.csg172.loftee.clinvar.all.raw.2 <- as.data.frame(raw.2[unique(c(kim_ST1.csg172.loftee.all$SNP,kim_ST1.csg172.clinvar.all$SNP)),])
kim_ST1.csg172.loftee.clinvar.all.raw.2$Carriers_by_variant <- rowSums(kim_ST1.csg172.loftee.clinvar.all.raw.2, na.rm = T)
kim_ST1.csg172.loftee.clinvar.all.raw.2$SNP <- row.names(kim_ST1.csg172.loftee.clinvar.all.raw.2)
kim_ST1.csg172.loftee.clinvar.all.raw.2 <- cbind.data.frame(SNP=kim_ST1.csg172.loftee.clinvar.all.raw.2$SNP, Carriers_by_variant=kim_ST1.csg172.loftee.clinvar.all.raw.2$Carriers_by_variant)
kim_ST1.csg172.loftee.clinvar.all.raw.2$GENE <- all.tools.genes$GENE[match(kim_ST1.csg172.loftee.clinvar.all.raw.2$SNP, all.tools.genes$SNP)]
kim_ST1.csg172.loftee.clinvar.all.raw.2 <- kim_ST1.csg172.loftee.clinvar.all.raw.2 %>%
  group_by(GENE) %>%
  mutate(variants_n = n())

# kim_ST1.csg172.loftee.clinvar.all.raw.2 <- kim_ST1.csg172.loftee.clinvar.all.raw.2[-1]

write.table(kim_ST1.csg172.loftee.clinvar.all.raw.2, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/kim_ST1.csg172.loftee.Or.clinvar.all.raw.2.txt", row.names = F, col.names = T, quote = F, sep = "\t")

kim_ST1.csg172.loftee.clinvar.all.raw.2_summary <- kim_ST1.csg172.loftee.clinvar.all.raw.2 %>%
  group_by(GENE) %>%
  summarize(
    total_carriers = sum(Carriers_by_variant),
    total_variants = first(variants_n)
  )

write.table(kim_ST1.csg172.loftee.clinvar.all.raw.2_summary, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/kim_ST1.csg172.loftee.Or.clinvar.all.raw.2_summary.txt", row.names = F, col.names = T, quote = F, sep = "\t")


###############
## BCC panel ##
###############
## Clinvar 
BCC.panel.clinvar.all.raw.2 <- as.data.frame(raw.2[BCC.panel.clinvar.all$SNP,])
BCC.panel.clinvar.all.raw.2$Carriers_by_variant <- rowSums(BCC.panel.clinvar.all.raw.2, na.rm = T)
BCC.panel.clinvar.all.raw.2$SNP <- row.names(BCC.panel.clinvar.all.raw.2)
BCC.panel.clinvar.all.raw.2 <- cbind.data.frame(SNP=BCC.panel.clinvar.all.raw.2$SNP, Carriers_by_variant=BCC.panel.clinvar.all.raw.2$Carriers_by_variant)
BCC.panel.clinvar.all.raw.2$GENE <- all.tools.genes$GENE[match(BCC.panel.clinvar.all.raw.2$SNP, all.tools.genes$SNP)]
BCC.panel.clinvar.all.raw.2 <- BCC.panel.clinvar.all.raw.2 %>%
  group_by(GENE) %>%
  mutate(variants_n = n())

# BCC.panel.clinvar.all.raw.2 <- BCC.panel.clinvar.all.raw.2[BCC.panel.clinvar.all.raw.2$Carriers_by_variant != 0,]

write.table(BCC.panel.clinvar.all.raw.2, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/BCC.panel.clinvar.all.raw.2.txt", row.names = F, col.names = T, quote = F, sep = "\t")

BCC.panel.clinvar.all.raw.2_summary <- BCC.panel.clinvar.all.raw.2 %>%
  group_by(GENE) %>%
  summarize(
    total_carriers = sum(Carriers_by_variant),
    total_variants = first(variants_n)
  )

write.table(BCC.panel.clinvar.all.raw.2_summary, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/BCC.panel.clinvar.all.raw.2_summary.txt", row.names = F, col.names = T, quote = F, sep = "\t")


## loftee 
BCC.panel.loftee.all.raw.2 <- as.data.frame(raw.2[BCC.panel.loftee.all$SNP,])
BCC.panel.loftee.all.raw.2$Carriers_by_variant <- rowSums(BCC.panel.loftee.all.raw.2, na.rm = T)
BCC.panel.loftee.all.raw.2$SNP <- row.names(BCC.panel.loftee.all.raw.2)
BCC.panel.loftee.all.raw.2 <- cbind.data.frame(SNP=BCC.panel.loftee.all.raw.2$SNP, Carriers_by_variant=BCC.panel.loftee.all.raw.2$Carriers_by_variant)
BCC.panel.loftee.all.raw.2$GENE <- all.tools.genes$GENE[match(BCC.panel.loftee.all.raw.2$SNP, all.tools.genes$SNP)]
BCC.panel.loftee.all.raw.2 <- BCC.panel.loftee.all.raw.2 %>%
  group_by(GENE) %>%
  mutate(variants_n = n())

# BCC.panel.loftee.all.raw.2 <- BCC.panel.loftee.all.raw.2[-1]

write.table(BCC.panel.loftee.all.raw.2, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/BCC.panel.loftee.all.raw.2.txt", row.names = F, col.names = T, quote = F, sep = "\t")

BCC.panel.loftee.all.raw.2_summary <- BCC.panel.loftee.all.raw.2 %>%
  group_by(GENE) %>%
  summarize(
    total_carriers = sum(Carriers_by_variant),
    total_variants = first(variants_n)
  )

write.table(BCC.panel.loftee.all.raw.2_summary, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/BCC.panel.loftee.all.raw.2_summary.txt", row.names = F, col.names = T, quote = F, sep = "\t")


## loftee or clinvar
BCC.panel.loftee.clinvar.all.raw.2 <- as.data.frame(raw.2[unique(c(BCC.panel.loftee.all$SNP,BCC.panel.clinvar.all$SNP)),])
BCC.panel.loftee.clinvar.all.raw.2$Carriers_by_variant <- rowSums(BCC.panel.loftee.clinvar.all.raw.2, na.rm = T)
BCC.panel.loftee.clinvar.all.raw.2$SNP <- row.names(BCC.panel.loftee.clinvar.all.raw.2)
BCC.panel.loftee.clinvar.all.raw.2 <- cbind.data.frame(SNP=BCC.panel.loftee.clinvar.all.raw.2$SNP, Carriers_by_variant=BCC.panel.loftee.clinvar.all.raw.2$Carriers_by_variant)
BCC.panel.loftee.clinvar.all.raw.2$GENE <- all.tools.genes$GENE[match(BCC.panel.loftee.clinvar.all.raw.2$SNP, all.tools.genes$SNP)]
BCC.panel.loftee.clinvar.all.raw.2 <- BCC.panel.loftee.clinvar.all.raw.2 %>%
  group_by(GENE) %>%
  mutate(variants_n = n())

# BCC.panel.loftee.clinvar.all.raw.2 <- BCC.panel.loftee.clinvar.all.raw.2[-1]

write.table(BCC.panel.loftee.clinvar.all.raw.2, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/BCC.panel.loftee.Or.clinvar.all.raw.2.txt", row.names = F, col.names = T, quote = F, sep = "\t")

BCC.panel.loftee.clinvar.all.raw.2_summary <- BCC.panel.loftee.clinvar.all.raw.2 %>%
  group_by(GENE) %>%
  summarize(
    total_carriers = sum(Carriers_by_variant),
    total_variants = first(variants_n)
  )

write.table(BCC.panel.loftee.clinvar.all.raw.2_summary, "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis/BCC.panel.loftee.Or.clinvar.all.raw.2_summary.txt", row.names = F, col.names = T, quote = F, sep = "\t")



################################

## Updated 02/08/2024
# Then, Achal, you need to count autosomal dominant and recessive genes
# separately. For dominant genes, those who carry 1 or 2 copies of P/LP variant
# would be carriers. For recessive genes, only those who carry both copies of
# P/LP variants should be designated as carriers. Could you do this and see what
# you get?

## Editing this part based on email on 0.2/08/2024 from Yadav
raw.3 <- raw[raw$IID %in% subdf$studyid,]
raw.3 <- raw.3[raw.3$FID %in% carrier.status$FID,]
raw.3 <- raw.3[grepl("FID|IID|chr", colnames(raw.3))]

rownames(raw.3) <- raw.3$IID
raw.3 <- raw.3[-c(1:2)]

#############
## Clinvar ##
#############
col.extract <- colnames(raw.3)[colnames(raw.3) %in% kim_ST1.csg172.clinvar.all$SNP]


# 1. first check with unknown as AR
kim_ST1.csg172.AR.unknown.clinvar <- kim_ST1.csg172[!grepl("^dominant$", kim_ST1.csg172$Mode.of.Inheritance),]
kim_ST1.csg172.AR.unknown.clinvar.var <- kim_ST1.csg172.clinvar.all$SNP[kim_ST1.csg172.clinvar.all$ANN....GENE %in% kim_ST1.csg172.AR.unknown.clinvar$Gene.Name]


# 2. then unknown as AD
kim_ST1.csg172.AR.clinvar <- kim_ST1.csg172[!grepl("^dominant$|unknown", kim_ST1.csg172$Mode.of.Inheritance),]
kim_ST1.csg172.AR.clinvar.var <- kim_ST1.csg172.clinvar.all$SNP[kim_ST1.csg172.clinvar.all$ANN....GENE %in% kim_ST1.csg172.AR$Gene.Name]


# Unknown as AR
raw.4.AR.unknown.as.recessive.clinvar <- raw.3[, c(col.extract)]
## If recessive is 1, make it zero
sum(raw.4.AR.unknown.as.recessive.clinvar[,kim_ST1.csg172.AR.unknown.clinvar.var] == 1, na.rm = T)
# 183

raw.4.AR.unknown.as.recessive.clinvar[,kim_ST1.csg172.AR.unknown.clinvar.var] <- replace(raw.4.AR.unknown.as.recessive.clinvar[, kim_ST1.csg172.AR.unknown.clinvar.var], raw.4.AR.unknown.as.recessive.clinvar[, kim_ST1.csg172.AR.unknown.clinvar.var] == 1, 0)
raw.4.AR.unknown.as.recessive.clinvar$kim_ST1.csg172.clinvar.all.status <- ifelse(rowSums(raw.4.AR.unknown.as.recessive.clinvar[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

# Unknown as AD
raw.4.AR.only.as.recessive.clinvar <- raw.3[, c(col.extract)]
## If recessive is 1, make it zero
sum(raw.4.AR.only.as.recessive.clinvar[,kim_ST1.csg172.AR.clinvar.var] == 1, na.rm = T)
# 173

raw.4.AR.only.as.recessive.clinvar[,kim_ST1.csg172.AR.clinvar.var] <- replace(raw.4.AR.only.as.recessive.clinvar[, kim_ST1.csg172.AR.clinvar.var], raw.4.AR.only.as.recessive.clinvar[, kim_ST1.csg172.AR.clinvar.var] == 1, 0)
raw.4.AR.only.as.recessive.clinvar$kim_ST1.csg172.clinvar.all.status <- ifelse(rowSums(raw.4.AR.only.as.recessive.clinvar[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

############
## Loftee ##
############
## Extract variant carrier status
col.extract <- colnames(raw.3)[colnames(raw.3) %in% kim_ST1.csg172.loftee.all$SNP]


# 1. first check with unknown as AR
kim_ST1.csg172.AR.unknown.loftee <- kim_ST1.csg172[!grepl("^dominant$", kim_ST1.csg172$Mode.of.Inheritance),]
kim_ST1.csg172.AR.unknown.loftee.var <- kim_ST1.csg172.loftee.all$SNP[kim_ST1.csg172.loftee.all$SYMBOL %in% kim_ST1.csg172.AR.unknown.loftee$Gene.Name]


# 2. then unknown as AD
kim_ST1.csg172.AR.loftee <- kim_ST1.csg172[!grepl("^dominant$|unknown", kim_ST1.csg172$Mode.of.Inheritance),]
kim_ST1.csg172.AR.loftee.var <- kim_ST1.csg172.loftee.all$SNP[kim_ST1.csg172.loftee.all$SYMBOL %in% kim_ST1.csg172.AR.loftee$Gene.Name]

# Unknown as AR
raw.4.AR.unknown.as.recessive.loftee <- raw.3[, c(col.extract)]
## If recessive is 1, make it zero
sum(raw.4.AR.unknown.as.recessive.loftee[,kim_ST1.csg172.AR.unknown.loftee.var] == 1, na.rm = T)
# 258

raw.4.AR.unknown.as.recessive.loftee[,kim_ST1.csg172.AR.unknown.loftee.var] <- replace(raw.4.AR.unknown.as.recessive.loftee[, kim_ST1.csg172.AR.unknown.loftee.var], raw.4.AR.unknown.as.recessive.loftee[, kim_ST1.csg172.AR.unknown.loftee.var] == 1, 0)
raw.4.AR.unknown.as.recessive.loftee$kim_ST1.csg172.loftee.all.status <- ifelse(rowSums(raw.4.AR.unknown.as.recessive.loftee[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

# Unknown as AD
raw.4.AR.only.as.recessive.loftee <- raw.3[, c(col.extract)]
## If recessive is 1, make it zero
sum(raw.4.AR.only.as.recessive.loftee[,kim_ST1.csg172.AR.loftee.var] == 1, na.rm = T)
# 242

raw.4.AR.only.as.recessive.loftee[,kim_ST1.csg172.AR.loftee.var] <- replace(raw.4.AR.only.as.recessive.loftee[, kim_ST1.csg172.AR.loftee.var], raw.4.AR.only.as.recessive.loftee[, kim_ST1.csg172.AR.loftee.var] == 1, 0)
raw.4.AR.only.as.recessive.loftee$kim_ST1.csg172.loftee.all.status <- ifelse(rowSums(raw.4.AR.only.as.recessive.loftee[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

############
## snpeff ##
############

col.extract <- colnames(raw.3)[colnames(raw.3) %in% kim_ST1.csg172.snpeff.all$SNP]

# 1. first check with unknown as AR
kim_ST1.csg172.AR.unknown.snpeff <- kim_ST1.csg172[!grepl("^dominant$", kim_ST1.csg172$Mode.of.Inheritance),]
kim_ST1.csg172.AR.unknown.snpeff.var <- kim_ST1.csg172.snpeff.all$SNP[kim_ST1.csg172.snpeff.all$ANN....GENE %in% kim_ST1.csg172.AR.unknown.snpeff$Gene.Name]


# 2. then unknown as AD
kim_ST1.csg172.AR.snpeff <- kim_ST1.csg172[!grepl("^dominant$|unknown", kim_ST1.csg172$Mode.of.Inheritance),]
kim_ST1.csg172.AR.snpeff.var <- kim_ST1.csg172.snpeff.all$SNP[kim_ST1.csg172.snpeff.all$ANN....GENE %in% kim_ST1.csg172.AR.snpeff$Gene.Name]

# Unknown as AR
raw.4.AR.unknown.as.recessive.snpeff <- raw.3[, c(col.extract)]
## If recessive is 1, make it zero
sum(raw.4.AR.unknown.as.recessive.snpeff[,kim_ST1.csg172.AR.unknown.snpeff.var] == 1, na.rm = T)
# 121

raw.4.AR.unknown.as.recessive.snpeff[,kim_ST1.csg172.AR.unknown.snpeff.var] <- replace(raw.4.AR.unknown.as.recessive.snpeff[, kim_ST1.csg172.AR.unknown.snpeff.var], raw.4.AR.unknown.as.recessive.snpeff[, kim_ST1.csg172.AR.unknown.snpeff.var] == 1, 0)
raw.4.AR.unknown.as.recessive.snpeff$kim_ST1.csg172.snpeff.all.status <- ifelse(rowSums(raw.4.AR.unknown.as.recessive.snpeff[, c(col.extract)], na.rm = T)> 0, "Yes", "No")

# Unknown as AD
raw.4.AR.only.as.recessive.snpeff <- raw.3[, c(col.extract)]
## If recessive is 1, make it zero
sum(raw.4.AR.only.as.recessive.snpeff[,kim_ST1.csg172.AR.snpeff.var] == 1, na.rm = T)
# 109

raw.4.AR.only.as.recessive.snpeff[,kim_ST1.csg172.AR.snpeff.var] <- replace(raw.4.AR.only.as.recessive.snpeff[, kim_ST1.csg172.AR.snpeff.var], raw.4.AR.only.as.recessive.snpeff[, kim_ST1.csg172.AR.snpeff.var] == 1, 0)
raw.4.AR.only.as.recessive.snpeff$kim_ST1.csg172.snpeff.all.status <- ifelse(rowSums(raw.4.AR.only.as.recessive.snpeff[, c(col.extract)], na.rm = T)> 0, "Yes", "No")




carrier_status_unknown.as.AR <- cbind.data.frame(IID = rownames(raw.4.AR.unknown.as.recessive.clinvar), kim_ST1.csg172.clinvar.all.status= raw.4.AR.unknown.as.recessive.clinvar$kim_ST1.csg172.clinvar.all.status, kim_ST1.csg172.loftee.all.status= raw.4.AR.unknown.as.recessive.loftee$kim_ST1.csg172.loftee.all.status, kim_ST1.csg172.snpeff.all.status = raw.4.AR.unknown.as.recessive.snpeff$kim_ST1.csg172.snpeff.all.status)
carrier_status_unknown.as.AD <- cbind.data.frame(IID = rownames(raw.4.AR.only.as.recessive.clinvar), kim_ST1.csg172.clinvar.all.status= raw.4.AR.only.as.recessive.clinvar$kim_ST1.csg172.clinvar.all.status, kim_ST1.csg172.loftee.all.status= raw.4.AR.only.as.recessive.loftee$kim_ST1.csg172.loftee.all.status, kim_ST1.csg172.snpeff.all.status = raw.4.AR.only.as.recessive.snpeff$kim_ST1.csg172.snpeff.all.status)

save(list = c("carrier_status_unknown.as.AR", "carrier_status_unknown.as.AD"), file = "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/BCC/files_shared_by_cindy/analysis//CSG172_carrier_updated_for_AR_genes.RData")


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



################################

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

