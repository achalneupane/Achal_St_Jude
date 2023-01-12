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

GWAS$A1_frq[is.na(GWAS$MAF)] <- FRQ$A1[match(GWAS$KEY1, FRQ$SNP)][is.na(GWAS$A1_frq)]
GWAS$A1_frq[is.na(GWAS$MAF)] <- FRQ$A1[match(GWAS$KEY2, FRQ$SNP)][is.na(GWAS$A1_frq)]
GWAS$A1_frq[is.na(GWAS$MAF)] <- FRQ$A1[match(GWAS$KEY3, FRQ$SNP)][is.na(GWAS$A1_frq)]
GWAS$A1_frq[is.na(GWAS$MAF)] <- FRQ$A1[match(GWAS$KEY4, FRQ$SNP)][is.na(GWAS$A1_frq)]

GWAS$A2_frq[is.na(GWAS$MAF)] <- FRQ$A2[match(GWAS$KEY1, FRQ$SNP)][is.na(GWAS$A2_frq)]
GWAS$A2_frq[is.na(GWAS$MAF)] <- FRQ$A2[match(GWAS$KEY2, FRQ$SNP)][is.na(GWAS$A2_frq)]
GWAS$A2_frq[is.na(GWAS$MAF)] <- FRQ$A2[match(GWAS$KEY3, FRQ$SNP)][is.na(GWAS$A2_frq)]
GWAS$A2_frq[is.na(GWAS$MAF)] <- FRQ$A2[match(GWAS$KEY4, FRQ$SNP)][is.na(GWAS$A2_frq)]

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
  
  GWAS$A1_frq[is.na(GWAS$MAF)] <- FRQ$A1[match(GWAS$KEY1, FRQ$SNP)][is.na(GWAS$A1_frq)]
  GWAS$A1_frq[is.na(GWAS$MAF)] <- FRQ$A1[match(GWAS$KEY2, FRQ$SNP)][is.na(GWAS$A1_frq)]
  GWAS$A1_frq[is.na(GWAS$MAF)] <- FRQ$A1[match(GWAS$KEY3, FRQ$SNP)][is.na(GWAS$A1_frq)]
  GWAS$A1_frq[is.na(GWAS$MAF)] <- FRQ$A1[match(GWAS$KEY4, FRQ$SNP)][is.na(GWAS$A1_frq)]
  
  GWAS$A2_frq[is.na(GWAS$MAF)] <- FRQ$A2[match(GWAS$KEY1, FRQ$SNP)][is.na(GWAS$A2_frq)]
  GWAS$A2_frq[is.na(GWAS$MAF)] <- FRQ$A2[match(GWAS$KEY2, FRQ$SNP)][is.na(GWAS$A2_frq)]
  GWAS$A2_frq[is.na(GWAS$MAF)] <- FRQ$A2[match(GWAS$KEY3, FRQ$SNP)][is.na(GWAS$A2_frq)]
  GWAS$A2_frq[is.na(GWAS$MAF)] <- FRQ$A2[match(GWAS$KEY4, FRQ$SNP)][is.na(GWAS$A2_frq)]
  
  
  print(table(is.na(GWAS$MAF)))
}


write.table(GWAS[,c("CHR", "SNP", "BP", "A1", "TEST", "NMISS", "OR", "SE", "L95", "U95", "STAT", "P", "BETA", "A2", "MAF", "A1_frq", "A2_frq")], "chr1_22_PCA_eur.assoc.logistic.clean.Psorted.withBETA.with.MAF", col.names = T, quote = F, row.names = F)
save.image("EUR.RData")


## Prepare MEGA files
## EUR
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/diabetes/meta_analysis_with_mr_Mega")
MEGA.1 <- fread("chr1_22_PCA_eur.assoc.logistic.clean.Psorted.withBETA.with.MAF", header = T)
MEGA.FINAL1 <- cbind.data.frame(MARKERNAME= MEGA.1$SNP, EA= MEGA.1$A1, NEA= MEGA.1$A2,OR = MEGA.1$OR, OR_95L= MEGA.1$L95, OR_95U = MEGA.1$U95, EAF = MEGA.1$MAF, N= 3113, CHROMOSOME = MEGA.1$CHR, POSITION= MEGA.1$BP)
head(MEGA.FINAL1)
write.table(MEGA.FINAL1, "MEGA.FINAL.EUR.txt", col.names = T, quote = F, row.names = F)

## AFR
MEGA.1 <- fread("chr1_22_PCA_afr.assoc.logistic.clean.Psorted.withBETA.with.MAF", header = T)
MEGA.FINAL2 <- cbind.data.frame(MARKERNAME= MEGA.1$SNP, EA= MEGA.1$A1, NEA= MEGA.1$A2,OR = MEGA.1$OR, OR_95L= MEGA.1$L95, OR_95U = MEGA.1$U95, EAF = MEGA.1$MAF, N= 575, CHROMOSOME = MEGA.1$CHR, POSITION= MEGA.1$BP)
head(MEGA.FINAL2)
write.table(MEGA.FINAL2, "MEGA.FINAL.AFR.txt", col.names = T, quote = F, row.names = F)


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

ALL <- merge(GWAS.AFR, GWAS.EUR, by = "SNP", all.x = T)
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

save.image("TOP_AFR_vars_5e-06_in_EUR_analysis.RDATA")
