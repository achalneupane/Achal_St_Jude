rm(list=ls())
merged.dat <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ccss_org_ccss_exp_ttn_bag3.pheno", header = T)

sjlife.eur.dat <- merged.dat[merged.dat$cohort==1,]
sjlife.afr.dat <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/sjlife_ttn_bag3_afr.pheno", header = T)

colnames(sjlife.eur.dat)
colnames(sjlife.afr.dat)
sjlife.eur.dat <- sjlife.eur.dat[c("FID", "IID", "CMP", "agedx", "agelstcontact", "gender", "anthra_jco_dose_any", "hrtavg",
                 "ejection_fraction_hrt", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
sjlife.eur.dat$ancestry <- "EUR"
sjlife.afr.dat <- sjlife.afr.dat[c("FID", "IID", "CMP", "agedx", "agelstcontact", "gender", "anthra_jco_dose_any", "hrtavg",
                                   "ejection_fraction_hrt", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")]
sjlife.afr.dat$ancestry <- "AFR"
sjlife <- rbind(sjlife.eur.dat, sjlife.afr.dat)
dim(sjlife)
# 1891
#########################
## Get echo parameters ##
#########################
library(haven)
machine <- read_sas("/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/echo_machine.sas7bdat")
research <- read_sas("/SJShare/SJCOMMON/ECC/SJLife/SJLIFE Data Freeze/2 Final Data SJLIFE/20200430/Clinical Data/echo_research.sas7bdat")
# sum(sjlife$IID %in% research$sjlid)
# # 1640
# research <- research[research$sjlid %in% sjlife$IID ,]
# research.get <- research [c("sjlid", "DateVisitStart", "visittype", "LV_Ejection_Fraction_3D", "LV_End_Diastolic_Volume_3D", "LV_End_Systolic_Volume_3D" ,
#                           "LV_Stroke_Volume_3D", "LVMassMM_Index", "LV_GLPS_AVG", "LV_Relative_Wall_Thickness")]

sum(sjlife$IID %in% machine$sjlid)
# 1891

machine <- machine[machine$sjlid %in% sjlife$IID ,]
machine.get <- machine [c("sjlid", "studydatetime", "visittype", "LV_Ejection_Fraction_3D", "LV_End_Diastolic_Volume_3D", "LV_End_Systolic_Volume_3D" ,
"LV_Stroke_Volume_3D", "LVMassMM_Index", "LV_GLPS_AVG", "LV_Relative_Wall_Thickness")]
table(machine.get$visittype)
# Other Visit SJLIFE Visit 1 SJLIFE Visit 2 SJLIFE Visit 3 SJLIFE Visit 4 SJLIFE Visit 5 SJLIFE Visit 6 SJLIFE Visit 7 
# 221           1891           1557            922            345             93             20              2

## get visit 1 and 2 (complete cases)
machine.get.visit.1.2 <- machine.get[(machine.get$visittype == "SJLIFE Visit 1"| machine.get$visittype == "SJLIFE Visit 2"),]
machine.get.visit.1.2 <- machine.get.visit.1.2[complete.cases(machine.get.visit.1.2),]
table(machine.get.visit.1.2$visittype)

## get visit 1, 2 and 3 (complete cases)
machine.get.visit.1.2.3 <- machine.get[(machine.get$visittype == "SJLIFE Visit 1"| machine.get$visittype == "SJLIFE Visit 2" |
                                        machine.get$visittype == "SJLIFE Visit 3"),]
machine.get.visit.1.2.3 <- machine.get.visit.1.2.3[complete.cases(machine.get.visit.1.2.3),]
table(machine.get.visit.1.2.3$visittype)

## Add TTN and BAG3 top SNPs
TTN_BAG3 <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/sjlife_afr_to_concat_updated_ttn_bag3_top_snps_add.raw", header = T)
sjlife$chr2.178562809.T.C_C <- TTN_BAG3$chr2.178562809.T.C_C[match(sjlife$IID, TTN_BAG3$IID)]
sjlife$chr10.119670121.T.C_C <- TTN_BAG3$chr10.119670121.T.C_C[match(sjlife$IID, TTN_BAG3$IID)]

## Add PLP
sjlife.LMNA.PLP <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes/LMNA_PLP_SJLIFE_recodeA.raw", header = T)
sjlife.LMNA.PLP$carrier <- ifelse(rowSums(sjlife.LMNA.PLP[grepl("chr", colnames(sjlife.LMNA.PLP))]) > 0, 1, 0)

sjlife.TNTT2.PLP <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes/TNTT2_PLP_SJLIFE_recodeA.raw", header = T)
sjlife.TNTT2.PLP$carrier <- ifelse(rowSums(sjlife.TNTT2.PLP[grepl("chr", colnames(sjlife.TNTT2.PLP))]) > 0, 1, 0)

sjlife.TTN.PLP <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes/TTN_PLP_SJLIFE_recodeA.raw", header = T)
sjlife.TTN.PLP$carrier <- ifelse(rowSums(sjlife.TTN.PLP[grepl("chr", colnames(sjlife.TTN.PLP))]) > 0, 1, 0)

sjlife.BAG3.PLP <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes/BAG3_PLP_SJLIFE_recodeA.raw", header = T)
sjlife.BAG3.PLP$carrier <- ifelse(rowSums(sjlife.BAG3.PLP[grepl("chr", colnames(sjlife.BAG3.PLP))]) > 0, 1, 0)

sjlife.ALL.PLP <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3//ALL_P_LP_combinations/rare_variant_analysis_v2/pablo_garcia_et_al_nine_genes/ALL_PLP_SJLIFE_recodeA.raw", header = T)
sjlife.ALL.PLP$carrier <- ifelse(rowSums(sjlife.ALL.PLP[grepl("chr", colnames(sjlife.ALL.PLP))]) > 0, 1, 0)

sjlife$LMNA.PLP_carrier <- sjlife.LMNA.PLP$carrier[match(sjlife$IID, sjlife.LMNA.PLP$IID)]
sjlife$TNTT2.PLP_carrier <- sjlife.TNTT2.PLP$carrier[match(sjlife$IID, sjlife.TNTT2.PLP$IID)]
sjlife$TTN.PLP_carrier <- sjlife.TTN.PLP$carrier[match(sjlife$IID, sjlife.TTN.PLP$IID)]
sjlife$BAG3.PLP_carrier <- sjlife.BAG3.PLP$carrier[match(sjlife$IID, sjlife.BAG3.PLP$IID)]
sjlife$ALL.PLP_carrier <- sjlife.ALL.PLP$carrier[match(sjlife$IID, sjlife.ALL.PLP$IID)]

merged_df <- merge(machine.get, sjlife, by.x = "sjlid", by.y = "IID", all.x = TRUE)
# merged_df <- merged_df[!is.na(merged_df$studydatetime),]

merged_df$studydatetime <- as.Date(merged_df$studydatetime)
merged_df <- merged_df %>%
  arrange(sjlid, desc(studydatetime)) %>%
  group_by(sjlid) %>%
  mutate(RecentVisitNumber = if_else(!is.na(studydatetime), row_number(), NA_integer_)) %>%
  ungroup()

sum(!duplicated(merged_df$sjlid) & is.na(merged_df$studydatetime))
# 4
merged_df$RecentVisitNumber[(!duplicated(merged_df$sjlid) & is.na(merged_df$studydatetime))] <- 1

echo.PLP.eur <- merged_df[merged_df$ancestry=="EUR",]
echo.PLP.afr <- merged_df[merged_df$ancestry=="AFR",]
save(echo.PLP.eur, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/echo.PLP.eur.RData")
save(echo.PLP.afr, file = "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/echo.PLP.afr.RData")

# visit1 <- echo.PLP.eur[which(echo.PLP.eur$RecentVisitNumber == 1),]
# dim(visit1)
# # [1] 1645   37
# table(is.na(visit1$LV_Ejection_Fraction_3D))
# # 1233   412 
# table(is.na(pheno_final$LV_Ejection_Fraction_3D))
# 1233   412 

visitBaseline <- echo.PLP.eur[which(echo.PLP.eur$visittype == "SJLIFE Visit 1"),]
dim(visitBaseline)
table(is.na(visitBaseline$LV_Ejection_Fraction_3D))
# [1] 854   791

# cc <- echo.PLP.eur[which(echo.PLP.eur$RecentVisitNumber == 1),]
# dim(cc)
# 
# table(is.na(cc$LV_Ejection_Fraction_3D))
# table(is.na(pheno_final$LV_Ejection_Fraction_3D))
# pheno_final$IID[!pheno_final$IID %in% cc$sjlid]
# table(pheno_final$IID == cc$sjlid)
# pheno_final$LV_Ejection_Fraction_3D == cc$LV_Ejection_Fraction_3D
# 
# # pheno_final$IID[!pheno_final$IID %in% cc$sjlid]
# # [1] "SJL1083001" "SJL1731909" "SJL4176516" "SJL4829107"