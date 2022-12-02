getwd()

setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno")

## Load ccss and sjlife data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_exp_Genetic_data_P_LP.Rdata")
CCSS_exp.ANY_PHENO <- PHENO.ANY_SN

load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.CCSS_org_Genetic_data_P_LP.Rdata")
CCSS_org.ANY_PHENO <- PHENO.ANY_SN
CCSS_org.ANY_PHENO$ccssid <- sapply(strsplit(CCSS_org.ANY_PHENO$ccssid, "_"), `[`, 1, simplify=T)

PHENO.ANY_SN <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/2_genetic_data_P_LP.rds")
SJLIFE.ANY_PHENO <- PHENO.ANY_SN

rm(list=setdiff(ls(), c("SJLIFE.ANY_PHENO", "CCSS_exp.ANY_PHENO", "CCSS_org.ANY_PHENO")))





# rename columns to make same as sjlife column names
colnames(ccss_org) [match(c("a_dx", "a_end", "SEX", "anth_DED", "HeartAvg", "CMP2plus"), colnames(ccss_org))] <- c("agedx", "agelstcontact", "gender", "anthra_jco_dose_any", "hrtavg", "CMP")



# rename columns to make same as sjlife column names
colnames(ccss_exp) [match(c("a_dx", "a_end", "SEX", "anth_DED", "HeartAvg", "CMP2plus"), colnames(ccss_exp))] <- c("agedx", "agelstcontact", "gender", "anthra_jco_dose_any", "hrtavg", "CMP")



############
## SJLIFE ##
############
## Read pheno file
sjlife <- read.table("sjlife_ttn_bag3.pheno", header = T)
head(sjlife)
dim(sjlife)
# 1645
sum(sjlife$IID %in% SJLIFE.ANY_PHENO$sjlid)
# 1645

## add cancer types
sjlife$diagrp <- SJLIFE.ANY_PHENO$diaggrp[match(sjlife$IID, SJLIFE.ANY_PHENO$sjlid)]



## N
n_sjlife = nrow(sjlife)

## CMP status (sjlife$CMP_EF_HEIRARCHY >= 2 was used for CMP CA/CO variable)
sjlife_CMP2plus_CA <- sum(sjlife$CMP == 2)
sjlife_CMP2plus_CA <- paste0(sjlife_CMP2plus_CA, " (", round((sjlife_CMP2plus_CA/n_sjlife)*100,1), "%)")

sjlife_CMP2plus_CO <- sum(sjlife$CMP == 1)
sjlife_CMP2plus_CO <- paste0(sjlife_CMP2plus_CO, " (", round((sjlife_CMP2plus_CO/n_sjlife)*100,1), "%)")


## Age at last contact: Median (IQR)
sjlife_agelstcontact <- paste0(round(median(sjlife$agelstcontact), 1), " (", round(quantile(sjlife$agelstcontact, prob=c(.25,.5,.75), type=1)[1], 1), "-" , round(quantile(sjlife$agelstcontact, prob=c(.25,.5,.75), type=1)[3],1), ")")

## Age at diagnosis
sjlife_agedx_0_4 <- sum(sjlife$agedx >= 0 & sjlife$agedx < 5)
sjlife_agedx_0_4 <- paste0(sjlife_agedx_0_4, " (", round((sjlife_agedx_0_4/n_sjlife)*100,1), "%)")

sjlife_agedx_5_9 <- sum(sjlife$agedx >= 5 & sjlife$agedx < 10)
sjlife_agedx_5_9 <- paste0(sjlife_agedx_5_9, " (", round((sjlife_agedx_5_9/n_sjlife)*100,1), "%)")

sjlife_agedx_10_14 <- sum(sjlife$agedx >= 10 & sjlife$agedx < 15)
sjlife_agedx_10_14 <- paste0(sjlife_agedx_10_14, " (", round((sjlife_agedx_10_14/n_sjlife)*100,1), "%)")

sjlife_agedx_15_or_plus <- sum(sjlife$agedx >= 15)
sjlife_agedx_15_or_plus <- paste0(sjlife_agedx_15_or_plus, " (", round((sjlife_agedx_15_or_plus/n_sjlife)*100,1), "%)")

## Gender (0 = Male; 1 = Female)
sjlife_male <- sum(sjlife$gender == 0)
sjlife_male <- paste0(sjlife_male, " (", round((sjlife_male/n_sjlife)*100,1), "%)")

sjlife_female <- sum(sjlife$gender == 1)
sjlife_female <- paste0(sjlife_female, " (", round((sjlife_female/n_sjlife)*100,1), "%)")

# ## Ejection fraction
# sjlife_ejectionfraction <- paste0(round(median(sjlife$ejection_fraction_hrt), 3)*100, " (", round(quantile(sjlife$ejection_fraction_hrt, prob=c(.25,.5,.75), type=1)[1], 3)*100, "-" , round(quantile(sjlife$ejection_fraction_hrt, prob=c(.25,.5,.75), type=1)[3],3)*100, ")")


##############
## CCSS_ORG ##
##############

ccss_org <- read.table("ccss_org_eur_cardiotoxic_exposed.pheno", header = T)
head(ccss_org)
dim(ccss_org)
# 3147
sum(ccss_org$IID %in% CCSS_org.ANY_PHENO$ccssid)
# 3147

ccss_org$diagrp <- CCSS_org.ANY_PHENO$diagnose[match(ccss_org$IID, CCSS_org.ANY_PHENO$ccssid)]

## Add age at CHF
ccss_age_CHF <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno/ccss_eur_all.txt", header = T)
ccss_org$IID[!ccss_org$IID %in% ccss_age_CHF$ccssid]
# [1] 15125841 15126401 15166684 15179378 # not present 
sum(ccss_org$IID %in% ccss_age_CHF$ccssid)
# 3143
ccss_org$agevent <- ccss_age_CHF$a_maxCHF15 [match(ccss_org$IID, ccss_age_CHF$ccssid)]
ccss_org$agedx <- ccss_org$a_dx
ccss_org$elapsedAGE <- ccss_org$agevent - ccss_org$agedx

# All elapsed negative age are converted to zero. ageevent probably used floor value so there are some negatives
ccss_org$elapsedAGE[ccss_org$elapsedAGE < 0] <- 0

## N
n_ccss_org = nrow(ccss_org)

## CMP status (CMP2plus was used for CMP CA/CO variable)
ccss_org_CMP2plus_CA <- sum(ccss_org$CMP == 2)
ccss_org_CMP2plus_CA <- paste0(ccss_org_CMP2plus_CA, " (", round((ccss_org_CMP2plus_CA/n_ccss_org)*100,1), "%)")

ccss_org_CMP2plus_CO <- sum(ccss_org$CMP == 1)
ccss_org_CMP2plus_CO <- paste0(ccss_org_CMP2plus_CO, " (", round((ccss_org_CMP2plus_CO/n_ccss_org)*100,1), "%)")


## Age at last contact: Median (IQR)
ccss_org_agelstcontact <- paste0(round(median(ccss_org$agelstcontact), 1), " (", round(quantile(ccss_org$agelstcontact, prob=c(.25,.5,.75), type=1)[1], 1), "-" , round(quantile(ccss_org$agelstcontact, prob=c(.25,.5,.75), type=1)[3],1), ")")

## Age at diagnosis
ccss_org_agedx_0_4 <- sum(ccss_org$agedx >= 0 & ccss_org$agedx < 5)
ccss_org_agedx_0_4 <- paste0(ccss_org_agedx_0_4, " (", round((ccss_org_agedx_0_4/n_ccss_org)*100,1), "%)")

ccss_org_agedx_5_9 <- sum(ccss_org$agedx >= 5 & ccss_org$agedx < 10)
ccss_org_agedx_5_9 <- paste0(ccss_org_agedx_5_9, " (", round((ccss_org_agedx_5_9/n_ccss_org)*100,1), "%)")

ccss_org_agedx_10_14 <- sum(ccss_org$agedx >= 10 & ccss_org$agedx < 15)
ccss_org_agedx_10_14 <- paste0(ccss_org_agedx_10_14, " (", round((ccss_org_agedx_10_14/n_ccss_org)*100,1), "%)")

ccss_org_agedx_15_or_plus <- sum(ccss_org$agedx >= 15)
ccss_org_agedx_15_or_plus <- paste0(ccss_org_agedx_15_or_plus, " (", round((ccss_org_agedx_15_or_plus/n_ccss_org)*100,1), "%)")

## Gender (0 = Male; 1 = Female)
ccss_org_male <- sum(ccss_org$gender == 0)
ccss_org_male <- paste0(ccss_org_male, " (", round((ccss_org_male/n_ccss_org)*100,1), "%)")

ccss_org_female <- sum(ccss_org$gender == 1)
ccss_org_female <- paste0(ccss_org_female, " (", round((ccss_org_female/n_ccss_org)*100,1), "%)")



##############
## CCSS_EXP ##
##############

## CCSS_exp
ccss_exp <- read.table("ccss_exp_eur_cardiotoxic_exposed.pheno", header = T)
head(ccss_exp)
dim(ccss_exp)
# 1484

sum(ccss_exp$IID %in% CCSS_exp.ANY_PHENO$ccssid)
# 1484

ccss_exp$diagrp <- CCSS_exp.ANY_PHENO$diagnose[match(ccss_exp$IID, CCSS_exp.ANY_PHENO$ccssid)]

## N
n_ccss_exp = nrow(ccss_exp)

## CMP status (CMP2plus was used for CMP CA/CO variable)
ccss_exp_CMP2plus_CA <- sum(ccss_exp$CMP == 2)
ccss_exp_CMP2plus_CA <- paste0(ccss_exp_CMP2plus_CA, " (", round((ccss_exp_CMP2plus_CA/n_ccss_exp)*100,1), "%)")

ccss_exp_CMP2plus_CO <- sum(ccss_exp$CMP == 1)
ccss_exp_CMP2plus_CO <- paste0(ccss_exp_CMP2plus_CO, " (", round((ccss_exp_CMP2plus_CO/n_ccss_exp)*100,1), "%)")


## Age at last contact: Median (IQR)
ccss_exp_agelstcontact <- paste0(round(median(ccss_exp$agelstcontact), 1), " (", round(quantile(ccss_exp$agelstcontact, prob=c(.25,.5,.75), type=1)[1], 1), "-" , round(quantile(ccss_exp$agelstcontact, prob=c(.25,.5,.75), type=1)[3],1), ")")

## Age at diagnosis
ccss_exp_agedx_0_4 <- sum(ccss_exp$agedx >= 0 & ccss_exp$agedx < 5)
ccss_exp_agedx_0_4 <- paste0(ccss_exp_agedx_0_4, " (", round((ccss_exp_agedx_0_4/n_ccss_exp)*100,1), "%)")

ccss_exp_agedx_5_9 <- sum(ccss_exp$agedx >= 5 & ccss_exp$agedx < 10)
ccss_exp_agedx_5_9 <- paste0(ccss_exp_agedx_5_9, " (", round((ccss_exp_agedx_5_9/n_ccss_exp)*100,1), "%)")

ccss_exp_agedx_10_14 <- sum(ccss_exp$agedx >= 10 & ccss_exp$agedx < 15)
ccss_exp_agedx_10_14 <- paste0(ccss_exp_agedx_10_14, " (", round((ccss_exp_agedx_10_14/n_ccss_exp)*100,1), "%)")

ccss_exp_agedx_15_or_plus <- sum(ccss_exp$agedx >= 15)
ccss_exp_agedx_15_or_plus <- paste0(ccss_exp_agedx_15_or_plus, " (", round((ccss_exp_agedx_15_or_plus/n_ccss_exp)*100,1), "%)")

## Gender (0 = Male; 1 = Female)
ccss_exp_male <- sum(ccss_exp$gender == 0)
ccss_exp_male <- paste0(ccss_exp_male, " (", round((ccss_exp_male/n_ccss_exp)*100,1), "%)")

ccss_exp_female <- sum(ccss_exp$gender == 1)
ccss_exp_female <- paste0(ccss_exp_female, " (", round((ccss_exp_female/n_ccss_exp)*100,1), "%)")

#######################
## Calculate P-value ##
#######################
library(sjstats)
library(coin)
## P-Age at last contact 

## Note : The Wilcoxon-Mann-Whitney test is for two groups. If you want to
# compare locations of more than two groups with a similar test, then you would
# usually go to the Kruskal-Wallis test. More: https://stats.stackexchange.com/questions/274146/clarification-on-mann-whitney-wilcoxon-test-on-two-to-three-groups

# cc <- cbind.data.frame(age=c(sjlife$agelstcontact, ccss_org$agelstcontact, ccss_exp$agelstcontact),
#                        cohort=c(rep("sjlife", nrow(sjlife)), rep("ccss_org", nrow(ccss_org)), rep("ccss_exp", nrow(ccss_exp))))
# cc <-mannwhitney(cc, age, cohort)
# cc$df$groups <- paste0(cc$df$grp1.label , " Vs ", cc$df$grp2.label)
# 
# # paste(cc$df$groups, cc$df$p, sep = "=")
# cc  = rbind.data.frame(cc$df$groups, cc$df$p)
# colnames(cc) <- cc[1,]
# agelstcontact_p <- cc[-1,]


cc <- cbind.data.frame(age=c(sjlife$agelstcontact, ccss_org$agelstcontact, ccss_exp$agelstcontact),
                       cohort=c(rep("sjlife", nrow(sjlife)), rep("ccss_org", nrow(ccss_org)), rep("ccss_exp", nrow(ccss_exp))))
cc <-kruskal.test(age ~ cohort, data = cc)

## P-Age at diagnosis
cc <- cbind.data.frame(age=c(sjlife$agedx, ccss_org$agedx, ccss_exp$agedx),
                       cohort=c(rep("sjlife", nrow(sjlife)), rep("ccss_org", nrow(ccss_org)), rep("ccss_exp", nrow(ccss_exp))))
cc <-mannwhitney(cc, age, cohort)
cc$df$groups <- paste0(cc$df$grp1.label , " Vs ", cc$df$grp2.label)

# paste(cc$df$groups, cc$df$p, sep = "=")
cc  = rbind.data.frame(cc$df$groups, cc$df$p)
colnames(cc) <- cc[1,]
agedx_p <- cc[-1,]

## P-Sex
cc <- cbind.data.frame(sjlife=table(sjlife$gender), ccss_org=table(ccss_org$gender), ccss_exp=table(ccss_exp$gender))
cc <- as.table(rbind(M=as.numeric(cc[1,-1]), F=as.numeric(cc[2,-1])))
colnames(cc) <- c("SJLIFE", "CCSS_org_hrc", "CCSS_exp_wgs")

sex_p <- chisq.test(cc)
#################################
## Finalized demographic table ##
#################################

SJLIFE_VARS <- c(CA=sjlife_CMP2plus_CA, CO=sjlife_CMP2plus_CO, AGE_AT_LAST_CONTACT = sjlife_agelstcontact, agedx_0_4 = sjlife_agedx_0_4, agedx_5_9 = sjlife_agedx_5_9, agedx_10_14 = sjlife_agedx_10_14, agedx_15orPLUS = sjlife_agedx_15_or_plus, MALE = sjlife_male, FEMALE = sjlife_female)
CCSS_ORG_VARS <- c(CA=ccss_org_CMP2plus_CA, CO=ccss_org_CMP2plus_CO, AGE_AT_LAST_CONTACT = ccss_org_agelstcontact, agedx_0_4 = ccss_org_agedx_0_4, agedx_5_9 = ccss_org_agedx_5_9, agedx_10_14 = ccss_org_agedx_10_14, agedx_15orPLUS = ccss_org_agedx_15_or_plus, MALE = ccss_org_male, FEMALE = ccss_org_female)
CCSS_EXP_VARS <- c(CA=ccss_exp_CMP2plus_CA, CO=ccss_exp_CMP2plus_CO, AGE_AT_LAST_CONTACT = ccss_exp_agelstcontact, agedx_0_4 = ccss_exp_agedx_0_4, agedx_5_9 = ccss_exp_agedx_5_9, agedx_10_14 = ccss_exp_agedx_10_14, agedx_15orPLUS = ccss_exp_agedx_15_or_plus, MALE = ccss_exp_male, FEMALE = ccss_exp_female)

final.df <- cbind.data.frame(SJLIFE_VARS, CCSS_ORG_VARS, CCSS_EXP_VARS)

write.table(final.df, "manuscript_table1.txt", col.names = T, sep = "\t", quote = F)
