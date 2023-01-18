# setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/meta_analysis_with_mr_Mega")
setwd("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/diabetes/meta_analysis_with_mr_Mega/v2_with_preQC_VCF")

library(data.table)

GWAS <- fread("chr1_22_PCA_afr.assoc.logistic.clean.Psorted.withBETA", header = T)
# dim(GWAS)
# [1] 10414749       14
GWAS$KEY1 <- paste0("chr", GWAS$CHR,":", GWAS$BP, ":", GWAS$A2,":", GWAS$A1)
GWAS$KEY2 <- paste0("chr", GWAS$CHR,":", GWAS$BP, ":", GWAS$A1,":", GWAS$A2)
GWAS$KEY3 <- chartr("ACGT", "TGCA", GWAS$KEY1)
GWAS$KEY4 <- chartr("ACGT", "TGCA", GWAS$KEY2)

GWAS$MAF <- NA
GWAS$A1_frq <- NA
GWAS$A2_frq <- NA

for (i in 1:22){
FRQ <- fread(paste0("tmp_AFR_filtered_variants_AFR_freq_chr",i,".frq"), header = T)

# dim(FRQ)
# 1543674       6


# sum(FRQ$SNP %in% GWAS$KEY1)

# 149298 should be



GWAS$MAF[is.na(GWAS$MAF)] <- FRQ$MAF[match(GWAS$KEY1, FRQ$SNP)][is.na(GWAS$MAF)]
GWAS$MAF[is.na(GWAS$MAF)] <- FRQ$MAF[match(GWAS$KEY2, FRQ$SNP)][is.na(GWAS$MAF)]
GWAS$MAF[is.na(GWAS$MAF)] <- FRQ$MAF[match(GWAS$KEY3, FRQ$SNP)][is.na(GWAS$MAF)]
GWAS$MAF[is.na(GWAS$MAF)] <- FRQ$MAF[match(GWAS$KEY4, FRQ$SNP)][is.na(GWAS$MAF)]

GWAS$A1_frq[is.na(GWAS$A1_frq)] <- FRQ$A1[match(GWAS$KEY1, FRQ$SNP)][is.na(GWAS$A1_frq)]
GWAS$A1_frq[is.na(GWAS$A1_frq)] <- FRQ$A1[match(GWAS$KEY2, FRQ$SNP)][is.na(GWAS$A1_frq)]
GWAS$A1_frq[is.na(GWAS$A1_frq)] <- FRQ$A1[match(GWAS$KEY3, FRQ$SNP)][is.na(GWAS$A1_frq)]
GWAS$A1_frq[is.na(GWAS$MAF)] <- FRQ$A1[match(GWAS$KEY4, FRQ$SNP)][is.na(GWAS$A1_frq)]

GWAS$A2_frq[is.na(GWAS$A2_frq)] <- FRQ$A2[match(GWAS$KEY1, FRQ$SNP)][is.na(GWAS$A2_frq)]
GWAS$A2_frq[is.na(GWAS$A2_frq)] <- FRQ$A2[match(GWAS$KEY2, FRQ$SNP)][is.na(GWAS$A2_frq)]
GWAS$A2_frq[is.na(GWAS$A2_frq)] <- FRQ$A2[match(GWAS$KEY3, FRQ$SNP)][is.na(GWAS$A2_frq)]
GWAS$A2_frq[is.na(GWAS$A2_frq)] <- FRQ$A2[match(GWAS$KEY4, FRQ$SNP)][is.na(GWAS$A2_frq)]

print(table(is.na(GWAS$MAF)))
}


write.table(GWAS[,c("CHR", "SNP", "BP", "A1", "TEST", "NMISS", "OR", "SE", "L95", "U95", "STAT", "P", "BETA", "A2", "MAF", "A1_frq", "A2_frq")], "chr1_22_PCA_afr.assoc.logistic.clean.Psorted.withBETA.with.MAF", col.names = T, quote = F, row.names = F)
save.image("AFR.Rdata")

## EUR

GWAS <- fread("chr1_22_PCA_eur.assoc.logistic.clean.Psorted.withBETA", header = T)
# dim(GWAS)
# [1] 10414749       14
GWAS$KEY1 <- paste0("chr", GWAS$CHR,":", GWAS$BP, ":", GWAS$A2,":", GWAS$A1)
GWAS$KEY2 <- paste0("chr", GWAS$CHR,":", GWAS$BP, ":", GWAS$A1,":", GWAS$A2)
GWAS$KEY3 <- chartr("ACGT", "TGCA", GWAS$KEY1)
GWAS$KEY4 <- chartr("ACGT", "TGCA", GWAS$KEY2)

GWAS$MAF <- NA
GWAS$A1_frq <- NA
GWAS$A2_frq <- NA

for (i in 1:22){
  FRQ <- fread(paste0("tmp_EUR_filtered_variants_EUR_freq_chr",i,".frq"), header = T)
  
  # dim(FRQ)
  # 1543674       6
  
  
  # sum(FRQ$SNP %in% GWAS$KEY1)
  
  # 149298 should be
  
  
  
  GWAS$MAF[is.na(GWAS$MAF)] <- FRQ$MAF[match(GWAS$KEY1, FRQ$SNP)][is.na(GWAS$MAF)]
  GWAS$MAF[is.na(GWAS$MAF)] <- FRQ$MAF[match(GWAS$KEY2, FRQ$SNP)][is.na(GWAS$MAF)]
  GWAS$MAF[is.na(GWAS$MAF)] <- FRQ$MAF[match(GWAS$KEY3, FRQ$SNP)][is.na(GWAS$MAF)]
  GWAS$MAF[is.na(GWAS$MAF)] <- FRQ$MAF[match(GWAS$KEY4, FRQ$SNP)][is.na(GWAS$MAF)]
  
  GWAS$A1_frq[is.na(GWAS$A1_frq)] <- FRQ$A1[match(GWAS$KEY1, FRQ$SNP)][is.na(GWAS$A1_frq)]
  GWAS$A1_frq[is.na(GWAS$A1_frq)] <- FRQ$A1[match(GWAS$KEY2, FRQ$SNP)][is.na(GWAS$A1_frq)]
  GWAS$A1_frq[is.na(GWAS$A1_frq)] <- FRQ$A1[match(GWAS$KEY3, FRQ$SNP)][is.na(GWAS$A1_frq)]
  GWAS$A1_frq[is.na(GWAS$MAF)] <- FRQ$A1[match(GWAS$KEY4, FRQ$SNP)][is.na(GWAS$A1_frq)]
  
  GWAS$A2_frq[is.na(GWAS$A2_frq)] <- FRQ$A2[match(GWAS$KEY1, FRQ$SNP)][is.na(GWAS$A2_frq)]
  GWAS$A2_frq[is.na(GWAS$A2_frq)] <- FRQ$A2[match(GWAS$KEY2, FRQ$SNP)][is.na(GWAS$A2_frq)]
  GWAS$A2_frq[is.na(GWAS$A2_frq)] <- FRQ$A2[match(GWAS$KEY3, FRQ$SNP)][is.na(GWAS$A2_frq)]
  GWAS$A2_frq[is.na(GWAS$A2_frq)] <- FRQ$A2[match(GWAS$KEY4, FRQ$SNP)][is.na(GWAS$A2_frq)]
  
  
  print(table(is.na(GWAS$MAF)))
}


write.table(GWAS[,c("CHR", "SNP", "BP", "A1", "TEST", "NMISS", "OR", "SE", "L95", "U95", "STAT", "P", "BETA", "A2", "MAF", "A1_frq", "A2_frq")], "chr1_22_PCA_eur.assoc.logistic.clean.Psorted.withBETA.with.MAF", col.names = T, quote = F, row.names = F)
save.image("EUR.RData")


## Prepare MEGA files
## EUR
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/meta_analysis_with_mr_Mega/v2_with_preQC_VCF")
MEGA.1 <- fread("chr1_22_PCA_eur.assoc.logistic.clean.Psorted.withBETA.with.MAF", header = T)
sum(MEGA.1$A1 == MEGA.1$A1_frq & MEGA.1$A2 == MEGA.1$A2_frq, na.rm = T)
# 10012741
sum(MEGA.1$A1 != MEGA.1$A1_frq | MEGA.1$A2 != MEGA.1$A2_frq, na.rm = T)
# 173540
sum(MEGA.1$A1 == MEGA.1$A2_frq | MEGA.1$A2 == MEGA.1$A1_frq, na.rm = T)
# 173540
## change the maf of those where gwas A1 is not equal to A1 in freq file
MEGA.1$MAF[which(MEGA.1$A1 == MEGA.1$A2_frq | MEGA.1$A2 == MEGA.1$A1_frq)] <- 1-MEGA.1$MAF[which(MEGA.1$A1 == MEGA.1$A2_frq | MEGA.1$A2 == MEGA.1$A1_frq)]

# MEGA.1$KEY <- paste0(MEGA.1$CHR, ":", MEGA.1$BP, ":", MEGA.1$A2_frq, ":", MEGA.1$A1_frq)
# MEGA.1.saved <- MEGA.1

MEGA.FINAL1 <- cbind.data.frame(MARKERNAME= MEGA.1$SNP, EA= MEGA.1$A1, NEA= MEGA.1$A2,OR = MEGA.1$OR, OR_95L= MEGA.1$L95, OR_95U = MEGA.1$U95, EAF = MEGA.1$MAF, N= MEGA.1$NMISS, CHROMOSOME = MEGA.1$CHR, POSITION= MEGA.1$BP)
head(MEGA.FINAL1)
MEGA.FINAL1 <- MEGA.FINAL1[with(MEGA.FINAL1, order(MEGA.FINAL1$CHROMOSOME, MEGA.FINAL1$POSITION)), ]
write.table(MEGA.FINAL1, "MEGA.FINAL.EUR.txt", col.names = T, quote = F, row.names = F)

## AFR
MEGA.1 <- fread("chr1_22_PCA_afr.assoc.logistic.clean.Psorted.withBETA.with.MAF", header = T)
sum(MEGA.1$A1 == MEGA.1$A1_frq & MEGA.1$A2 == MEGA.1$A2_frq, na.rm = T)
# 9496105
sum(MEGA.1$A1 != MEGA.1$A1_frq | MEGA.1$A2 != MEGA.1$A2_frq, na.rm = T)
# 904628
sum(MEGA.1$A1 == MEGA.1$A2_frq | MEGA.1$A2 == MEGA.1$A1_frq, na.rm = T)
# 904628
sum((MEGA.1$A1 == MEGA.1$A2_frq), na.rm = T)
# 904628
## change the maf of those where gwas A1 is not equal to A1 in freq file
MEGA.1$MAF[which(MEGA.1$A1 == MEGA.1$A2_frq | MEGA.1$A2 == MEGA.1$A1_frq)] <- 1-MEGA.1$MAF[which(MEGA.1$A1 == MEGA.1$A2_frq | MEGA.1$A2 == MEGA.1$A1_frq)]

# MEGA.1$KEY <- paste0(MEGA.1$CHR, ":", MEGA.1$BP, ":", MEGA.1$A2_frq, ":", MEGA.1$A1_frq)
# MEGA.2.saved <- MEGA.1
MEGA.FINAL2 <- cbind.data.frame(MARKERNAME= MEGA.1$SNP, EA= MEGA.1$A1, NEA= MEGA.1$A2,OR = MEGA.1$OR, OR_95L= MEGA.1$L95, OR_95U = MEGA.1$U95, EAF = MEGA.1$MAF, N= MEGA.1$NMISS, CHROMOSOME = MEGA.1$CHR, POSITION= MEGA.1$BP)
head(MEGA.FINAL2)
MEGA.FINAL2 <- MEGA.FINAL2[with(MEGA.FINAL2, order(MEGA.FINAL2$CHROMOSOME, MEGA.FINAL2$POSITION)), ]
write.table(MEGA.FINAL2, "MEGA.FINAL.AFR.txt", col.names = T, quote = F, row.names = F)

## KEEP only the common variants between AFR and EUR
sum(MEGA.FINAL1$MARKERNAME %in% MEGA.FINAL2$MARKERNAME)
# 8899222

MEGA.FINAL2.common <- MEGA.FINAL2[match(MEGA.FINAL1$MARKERNAME, MEGA.FINAL2$MARKERNAME),]
MEGA.FINAL2.common<- MEGA.FINAL2.common[!is.na(MEGA.FINAL2.common$MARKERNAME),]

MEGA.FINAL1.common <- MEGA.FINAL1[match(MEGA.FINAL2$MARKERNAME, MEGA.FINAL1$MARKERNAME),]
MEGA.FINAL1.common<- MEGA.FINAL1.common[!is.na(MEGA.FINAL1.common$MARKERNAME),]

sum(MEGA.FINAL1.common$MARKERNAME %in% MEGA.FINAL2.common$MARKERNAME)
# 8899222

## Make same order based on SNP ID in both data
MEGA.FINAL2.common <- MEGA.FINAL2.common[match(MEGA.FINAL1.common$MARKERNAME, MEGA.FINAL2.common$MARKERNAME),]

sum(MEGA.FINAL1.common$MARKERNAME == MEGA.FINAL2.common$MARKERNAME)

write.table(MEGA.FINAL1.common, "MEGA.FINAL.EUR.overlapping.txt", col.names = T, quote = F, row.names = F)
write.table(MEGA.FINAL2.common, "MEGA.FINAL.AFR.overlapping.txt", col.names = T, quote = F, row.names = F)


## Read mega results:
## Overlapping variants
MEGA.res <- fread("MEGA.FINAL.RESULTS_overlapping_vars.txt.result", header = T)
cc <- MEGA.res[which(MEGA.res$`P-value_ancestry_het` < 5e-6),]


# ## Request 2
# (1) Get all variants with P<5e-6 in the AFR-only analysis
# (2) Match up variants in (1) with variants also present in the EUR analysis (i.e., replication): put these in an excel file
# (3) For all variants in (1), please send a rdata file with individual-level data for both AFR and EUR SJLIFE participants (I will find the email I got from Yadav and re-forward it to you so you can match up the variables included)
# Please let me know if this makes sense

GWAS.AFR <- fread("chr1_22_PCA_afr.assoc.logistic.clean.Psorted.withBETA", header = T)
sum(GWAS.AFR$P < 5e-6)
# 130
GWAS.AFR <- GWAS.AFR[GWAS.AFR$P < 5e-6,]
GWAS.AFR$CHR.BP <- paste(GWAS.AFR$CHR, GWAS.AFR$BP, sep = ":")

GWAS.EUR <- fread("chr1_22_PCA_eur.assoc.logistic.clean.Psorted.withBETA", header = T)
GWAS.EUR$CHR.BP <- paste(GWAS.EUR$CHR, GWAS.EUR$BP, sep = ":")

GWAS.EUR <- GWAS.EUR[GWAS.EUR$CHR.BP %in% GWAS.AFR$CHR.BP,]

GWAS.EUR$Exact.Match.in.AFR <- ifelse(GWAS.EUR$SNP %in% GWAS.AFR$SNP, "Yes", "No")

ALL <- cbind.data.frame(GWAS.AFR, GWAS.EUR[match(GWAS.AFR$SNP, GWAS.EUR$SNP),])
dim(ALL)
# 130 

colnames(ALL) <- gsub(".x","_AFR", colnames(ALL))
colnames(ALL) <- gsub(".y","_EUR", colnames(ALL))
write.table(ALL, "TOP_AFR_vars_5e-06_in_EUR_analysis.txt", col.names = T, quote = F, row.names = F)


## Check variants match
plink.extract <- read.table("PLINK_TOP_AFR_vars_5e-06_in_EUR_analysis_VAR_POS.bim")
sum(ALL$SNP %in% plink.extract$V2)
# 130

## Get genotype
raw <- read.table("PLINK_TOP_AFR_vars_5e-06_in_EUR_analysis_VAR_POS_recodeA.raw", header = T)
rownames(raw) <- raw$IID
raw <- raw[!grepl("FID|IID|PAT|MAT|SEX|PHENOTYPE", colnames(raw))]

HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
colnames(raw) = gsub("\\.", ":", HEADER)

sum(colnames(raw) %in% ALL$SNP)
# 130

raw <- raw[colnames(raw) %in% ALL$SNP]

load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/Results_final_to_Cindy/SJLIFE_T2D_GWAS_data.RData")
dat_afr <- dat_afr[c("FID","IID","t2d","gender","agedx","age_last_visit","BMIadj","aa_class_dose_5","maxabdrtdose","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","hemoglobin_a1c","glucose_lvl","insulin","labdt")]
sum(dat_afr$IID %in% row.names(raw))
# 574
dat_afr <- cbind.data.frame(dat_afr, raw[match(dat_afr$IID, row.names(raw)),])


dat_eur <- dat_eur[c("FID","IID","t2d","gender","agedx","age_last_visit","BMIadj","aa_class_dose_5","maxabdrtdose","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","hemoglobin_a1c","glucose_lvl","insulin","labdt")]
sum(dat_eur$IID %in% row.names(raw))
# 3102
dat_eur <- cbind.data.frame(dat_eur, raw[match(dat_eur$IID, row.names(raw)),])

dat_afr_5e_06 <- dat_afr
dat_eur_5e_06 <- dat_eur
# remove all that are rare in EUR 
rare.eur <- c("chr5:98340133:A:ACT", "chr5:98343634:T:C", "chr5:98343890:T:C", "chr5:98343892:C:G", "chr5:98350630:C:G", "chr5:175069792:T:C")
dim(dat_eur_5e_06)
sum(!colnames(dat_eur_5e_06) %in% rare.eur)
# 147
dat_eur_5e_06 <- dat_eur_5e_06[!colnames(dat_eur_5e_06) %in% rare.eur]

# colnames(dat_afr_5e_06)[!colnames(dat_afr_5e_06) %in% colnames(dat_eur_5e_06)]

save.image("TOP_AFR_vars_5e-06_in_EUR_analysis.RDATA")

## Chunk 4
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/METASOFT")
## Prepare Metasoft files
eur.gwas <- fread("chr1_22_PCA_eur.assoc.logistic.clean.Psorted.withBETA", header = T)
afr.gwas <- fread("chr1_22_PCA_afr.assoc.logistic.clean.Psorted.withBETA", header = T)

sum(eur.gwas$SNP %in% afr.gwas$SNP)

eur.gwas$afrBeta <- afr.gwas$BETA[match(eur.gwas$SNP, afr.gwas$SNP)]
eur.gwas$afrSE <- afr.gwas$SE[match(eur.gwas$SNP, afr.gwas$SNP)]


sum(!is.na(eur.gwas$afrBeta))
# 8899222

eur.gwas <- eur.gwas[!is.na(eur.gwas$afrBeta),]
metasoft <- cbind.data.frame(eur.gwas$SNP, eur.gwas$BETA, eur.gwas$SE, eur.gwas$afrBeta, eur.gwas$afrSE)
write.table(metasoft, "METASOFT_input_8899222_overlapping_variants.meta", col.names = F, quote = F, row.names = F)


## Read the METASOFT results
# metasoft.res <- fread("metasoft_res_edited", header = TRUE, sep = " ")
# head(metasoft.res)

metasoft.res <- fread("metasoft_res_edited", sep = "\t", header=F)
header <- c("RSID", "STUDY", "PVALUE_FE", "BETA_FE", "STD_FE", "PVALUE_RE", "BETA_RE", "STD_RE", "PVALUE_RE2", "STAT1_RE2", "STAT2_RE2")
# header <- c("RSID", "STUDY", "PVALUE_FE", "BETA_FE", "STD_FE", "PVALUE_RE", "BETA_RE", "STD_RE", "PVALUE_RE2", "STAT1_RE2", "STAT2_RE2", "PVALUE_BE", "I_SQUARE", "Q", "PVALUE_Q", "TAU_SQUARE", "PVALUES_OF_STUDIES", "MVALUES_OF_STUDIES")
length(header)
metasoft.res <- metasoft.res[,1:11]
colnames(metasoft.res) <- header 
head(metasoft.res)

metasoft.res.FE <- metasoft.res[metasoft.res$PVALUE_FE < 5e-6,c("RSID", "PVALUE_FE", "BETA_FE", "STD_FE")]
metasoft.res.RE <- metasoft.res[metasoft.res$PVALUE_RE < 5e-6,c("RSID", "PVALUE_RE", "BETA_RE", "STD_RE")]
metasoft.res.RE2 <- metasoft.res[metasoft.res$PVALUE_RE2 < 5e-6, c("RSID", "PVALUE_RE2", "STAT1_RE2", "STAT2_RE2")]



load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/Results_final_to_Cindy/TOP_AFR_vars_5e-06_in_EUR_analysis.RDATA")
metal.FE <- fread("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/meta_analysis_with_mr_Mega/v2_with_preQC_VCF/T2D_Meta_analysis_SJLIFE_EUR_and_SJLIFE_AFR_fixed_1.tbl", header = T)
metal.FE <- metal.FE[metal.FE$`P-value` < 5e-6,]


library(purrr)
library(dplyr)
merged_df <- reduce(list(metasoft.res.FE, metasoft.res.RE, metasoft.res.RE2), full_join, by = "RSID")
colnames(merged_df)[1] <- "MarkerName"

merged_df <- reduce(list(merged_df, metal.FE), full_join, by = "MarkerName")

merged_df <- cbind.data.frame(merged_df, eur.gwas[match(merged_df$MarkerName, eur.gwas$SNP),], afr.gwas[match(merged_df$MarkerName, afr.gwas$SNP),])

write.table(merged_df, "METASOFT_FE_RE_RE2_and_metal-FE.meta-analysis.txt", col.names = T, quote = F, row.names = F)

load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/Results_final_to_Cindy/SJLIFE_T2D_GWAS_data.RData")

old.cols <- unique(c(colnames(dat_afr), colnames(dat_eur)))

old.cols <- old.cols[grepl("^chr", old.cols)]

old.cols <- sub("\\.\\.+..*$", "", old.cols)

HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="", old.cols)
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
HEADER = gsub("\\.", ":", HEADER)
  

sum(merged_df$MarkerName %in% unique(c(colnames(dat_afr_5e_06), colnames(dat_eur_5e_06), HEADER)))
# 33

extract <- merged_df$MarkerName[!merged_df$MarkerName %in% unique(c(colnames(dat_afr_5e_06), colnames(dat_eur_5e_06), HEADER))]

write.table(extract, "extract_SNPs.txt", col.names = F, quote = F, row.names = F)



## Get genotype of Meta-analysis variants that were not included in previous Rdata

## Get genotype
raw <- read.table("PLINK_5e-06_meta-analysis_not_in_previous_RDATA_recodeA.raw", header = T)
rownames(raw) <- raw$IID
raw <- raw[!grepl("FID|IID|PAT|MAT|SEX|PHENOTYPE", colnames(raw))]

HEADER = sub(pattern="_[T,A,G,C,*]+",replacement="",colnames(raw))
HEADER = gsub(pattern=";rs\\d+",replacement="",HEADER)
HEADER = gsub("._...DEL", "",HEADER)
HEADER = gsub("....DEL.", ".<*:DEL>",HEADER)
colnames(raw) = gsub("\\.", ":", HEADER)

sum(colnames(raw) %in% merged_df$MarkerName)
# 185

raw$IID <- rownames(raw)
meta_analaysis_vars.AFR <- raw[raw$IID %in% dat_afr$IID,]
meta_analaysis_vars.EUR <- raw[raw$IID %in% dat_eur$IID,]


rm(list=setdiff(ls(), c("meta_analaysis_vars.AFR", "meta_analaysis_vars.EUR")))
save.image("meta-analysis_vars_missed_in_EUR_AFR.RData")



## Cindy also wanted variants first filtered by AFR-only results with P<5e-06 and extract the corresponding variants from EUR and Meta results?
colnames(ALL)[1:14] <- paste0(colnames(ALL)[1:14], ".AFR")
colnames(ALL)[15:30] <- paste0(colnames(ALL)[15:30], ".EUR")
ALL <- ALL[1:29]
colnames(ALL)[2] <- "MarkerName"
colnames(metasoft.res)[1] <- "MarkerName"


AFR.TOP <- reduce(list(ALL, metasoft.res, metal.FE), left_join, by = "MarkerName")

write.table(AFR.TOP, "TOP.AFR.only.with.P.5e-06.and.results.from.meta.txt", col.names = T, quote = F, row.names = F)
