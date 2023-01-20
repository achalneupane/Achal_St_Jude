####################################
## Function to create demographic ##
####################################
## NOTE: CCSS CHF data is based on self-report so ejection fraction is not available.
# df <- df
# n = nrow(df)

get_demographic <- function(df, n){
## Ejection fraction
ejection_fraction_hrt <- paste0(round((median(df$ejection_fraction_hrt, na.rm = T)*100), 1), " (", round((quantile(df$ejection_fraction_hrt, prob=c(.25,.5,.75), type=1, na.rm = T)*100)[1], 1), "-" , round((quantile(df$ejection_fraction_hrt, prob=c(.25,.5,.75), type=1, na.rm = T)*100)[3],1), ")")

## Elapsed time
df$elapsedAGE <- df$agevent - df$agedx

# All elapsed negative age are converted to zero. ageevent probably used floor value so there are some negatives
df$elapsedAGE[df$elapsedAGE < 0] <- 0
df_elapsedAGE <- paste0(round(median(df$elapsedAGE,  na.rm = T), 1), " (", round(quantile(df$elapsedAGE, prob=c(.25,.5,.75), type=1, na.rm = T)[1], 1), "-" , round(quantile(df$elapsedAGE, prob=c(.25,.5,.75), type=1,  na.rm = T)[3],1), ")")



## Age at last contact: Median (IQR)
agelstcontact <- paste0(round(median(df$agelstcontact, na.rm = T), 1), " (", round(quantile(df$agelstcontact, prob=c(.25,.5,.75), type=1, na.rm = T)[1], 1), "-" , round(quantile(df$agelstcontact, prob=c(.25,.5,.75), type=1, na.rm = T)[3],1), ")")

## Age at diagnosis
agedx <- paste0(round(median(df$agedx, na.rm = T), 1), " (", round(quantile(df$agedx, prob=c(.25,.5,.75), type=1, na.rm = T)[1], 1), "-" , round(quantile(df$agedx, prob=c(.25,.5,.75), type=1, na.rm = T)[3],1), ")")

## Gender (0 = Male; 1 = Female)
male <- sum(df$gender == 0, na.rm = T)
male <- paste0(male, " (", round((male/n)*100,1), "%)")

female <- sum(df$gender == 1, na.rm = T)
female <- paste0(female, " (", round((female/n)*100,1), "%)")

## Anthracycline dose (%)
anthra_0 <- paste0(sum(df$anthra_jco_dose_any == 0, na.rm = T), 
                   " (", round(sum(df$anthra_jco_dose_any == 0, na.rm = T)/nrow(df)*100, 1), "%)")

anthra_1_250 <- paste0(sum(df$anthra_jco_dose_any >= 1 & df$anthra_jco_dose_any <= 250, na.rm = T),
                       " (", round(sum(df$anthra_jco_dose_any >= 1 & df$anthra_jco_dose_any <= 250, na.rm = T)/nrow(df)*100, 1), "%)")

anthra_gt_250 <- paste0(sum(df$anthra_jco_dose_any > 250, na.rm = T), 
                        " (", round(sum(df$anthra_jco_dose_any > 250, na.rm = T)/nrow(df)*100, 1), "%)")

## Average heart radiation dose
heartrt_0 <- paste0(sum(df$hrtavg == 0, na.rm = T), 
                    " (", round(sum(df$hrtavg == 0, na.rm = T)/nrow(df)*100, 1), "%)")

heartrt_1_25 <- paste0(sum(df$hrtavg >= 1 & df$hrtavg <= 25, na.rm = T), 
                       " (", round(sum(df$hrtavg >= 1 & df$hrtavg <= 25, na.rm = T)/nrow(df)*100, 1), "%)")

heartrt_gt_25 <- paste0(sum(df$hrtavg > 25, na.rm = T), 
                        " (", round(sum(df$hrtavg > 25, na.rm = T)/nrow(df)*100, 1), "%)")

# df$diagrp[grepl("CNS|Nervous", df$diagrp, ignore.case = T)] <- 'Neuroblastoma'
# df$diagrp[grepl("NHL", df$diagrp, ignore.case = T)] <- 'Non-Hodgkin lymphoma'
# df$diagrp[grepl("HD", df$diagrp, ignore.case = T)] <- 'Hodgkin lymphoma'
# df$diagrp[grepl("Bone cancer", df$diagrp, ignore.case = T)] <- 'Osteosarcoma'


out <- setNames(rbind.data.frame(agedx, agelstcontact, df_elapsedAGE, ejection_fraction_hrt, male, female, anthra_0, anthra_1_250, anthra_gt_250, heartrt_0,
                                 heartrt_1_25, heartrt_gt_25), paste0("With CMP (n = ", n, ")"))

rownames(out) <- c("agedx", "agelstcontact", "df_elapsedAGE", "ejection_fraction_hrt", "male", "female", "anthra_0", "anthra_1_250", "anthra_gt_250",
                   "heartrt_0", "heartrt_1_25", "heartrt_gt_25")

out
}

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

rm(list=setdiff(ls(), c("SJLIFE.ANY_PHENO", "CCSS_exp.ANY_PHENO", "CCSS_org.ANY_PHENO", "get_demographic")))



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

load("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/1870_SJLIFE_samples_with_time2echo.RData")

sjlife$time2echo <- d.eur$time2echo[match(sjlife$IID, d.eur$IID)]
sjlife$agevent <- d.eur$agegraded[match(sjlife$IID, d.eur$IID)]
sjlife$agevent[sjlife$CMP == 1] <- NA # nullify controls age at event

# ccss_org_CMP2plus_CA <- sum(ccss_org$CMP == 2)
df <- sjlife[sjlife$CMP == 2,] # With CMP (Cases)
df <- sjlife[sjlife$CMP == 1,] # Without CMP (Controls)


## N
n = nrow(df)


sjlife_with_CMP <- get_demographic(df, n)
sjlife_without_CMP <- get_demographic(df, n)


##############
## CCSS_ORG ##
##############
load("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno/ccss_export_05062022.RData")

ccss_org <- read.table("ccss_org_eur_cardiotoxic_exposed.pheno", header = T)
head(ccss_org)
dim(ccss_org)
# 3147
sum(ccss_org$IID %in% CCSS_org.ANY_PHENO$ccssid)
# 3147

ccss_org$diagrp <- CCSS_org.ANY_PHENO$diagnose[match(ccss_org$IID, CCSS_org.ANY_PHENO$ccssid)]




# rename columns to make same as sjlife column names
colnames(ccss_org) [match(c("a_dx", "a_end", "SEX", "anth_DED", "HeartAvg", "CMP2plus"), colnames(ccss_org))] <- c("agedx", "agelstcontact", "gender", "anthra_jco_dose_any", "hrtavg", "CMP")



## Add age at CHF
ccss_age_CHF <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno/ccss_eur_all.txt", header = T)
ccss_org$IID[!ccss_org$IID %in% ccss_age_CHF$ccssid]
# [1] 15125841 15126401 15166684 15179378 # not present 
sum(ccss_org$IID %in% ccss_age_CHF$ccssid)
# 3143
ccss_org$agevent <- ccss_age_CHF$a_maxCHF15 [match(ccss_org$IID, ccss_age_CHF$ccssid)]


# ccss_org_CMP2plus_CA <- sum(ccss_org$CMP == 2)
df <- ccss_org[ccss_org$CMP == 2,] # With CMP (Cases)
df <- ccss_org[ccss_org$CMP == 1,] # Without CMP (Controls)

## Ejection fraction not available 
df$ejection_fraction_hrt <- NA

## N
n = nrow(df)


CCSS_org_with_CMP <- get_demographic(df, n)
CCSS_org_without_CMP <- get_demographic(df, n)


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

# rename columns to make same as sjlife column names
colnames(ccss_exp) [match(c("a_dx", "a_end", "SEX", "anth_DED", "HeartAvg", "CMP2plus"), colnames(ccss_exp))] <- c("agedx", "agelstcontact", "gender", "anthra_jco_dose_any", "hrtavg", "CMP")


ccss_exp$agevent <- ccss_age_CHF$a_maxCHF15 [match(ccss_exp$IID, ccss_age_CHF$ccssid)]

## Remove younger cases <5 years
ccss_exp.ca.remove <- ccss_exp[(ccss_exp$agevent - ccss_exp$agedx < 5) & ccss_exp$CMP == 2,]


# ccss_org_CMP2plus_CA <- sum(ccss_org$CMP == 2)
# df <- ccss_exp[ccss_exp$CMP == 2,] # With CMP (Cases)
df <- ccss_exp[(ccss_exp$agevent - ccss_exp$agedx > 5) & ccss_exp$CMP == 2,]
df <- ccss_exp[ccss_exp$CMP == 1,] # Without CMP (Controls)

## N
n = nrow(df)


## Ejection fraction not available 
df$ejection_fraction_hrt <- NA


CCSS_exp_with_CMP <- get_demographic(df, n)
CCSS_exp_without_CMP <- get_demographic(df, n)


#######################
## Write final table ##
#######################
final.df <- cbind.data.frame(sjlife_with_CMP , sjlife_without_CMP, CCSS_org_with_CMP, CCSS_org_without_CMP, CCSS_exp_with_CMP , CCSS_exp_without_CMP)

write.table(final.df, "Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/manuscript_table1.txt", col.names = T, sep = "\t", quote = F)








# 
# #######################
# ## Calculate P-value ##
# #######################
# library(sjstats)
# library(coin)
# ## P-Age at last contact 
# 
# ## Note : The Wilcoxon-Mann-Whitney test is for two groups. If you want to
# # compare locations of more than two groups with a similar test, then you would
# # usually go to the Kruskal-Wallis test. More: https://stats.stackexchange.com/questions/274146/clarification-on-mann-whitney-wilcoxon-test-on-two-to-three-groups
# 
# # cc <- cbind.data.frame(age=c(sjlife$agelstcontact, ccss_org$agelstcontact, ccss_exp$agelstcontact),
# #                        cohort=c(rep("sjlife", nrow(sjlife)), rep("ccss_org", nrow(ccss_org)), rep("ccss_exp", nrow(ccss_exp))))
# # cc <-mannwhitney(cc, age, cohort)
# # cc$df$groups <- paste0(cc$df$grp1.label , " Vs ", cc$df$grp2.label)
# # 
# # # paste(cc$df$groups, cc$df$p, sep = "=")
# # cc  = rbind.data.frame(cc$df$groups, cc$df$p)
# # colnames(cc) <- cc[1,]
# # agelstcontact_p <- cc[-1,]
# 
# 
# cc <- cbind.data.frame(age=c(sjlife$agelstcontact, ccss_org$agelstcontact, ccss_exp$agelstcontact),
#                        cohort=c(rep("sjlife", nrow(sjlife)), rep("ccss_org", nrow(ccss_org)), rep("ccss_exp", nrow(ccss_exp))))
# cc <-kruskal.test(age ~ cohort, data = cc)
# 
# ## P-Age at diagnosis
# cc <- cbind.data.frame(age=c(sjlife$agedx, ccss_org$agedx, ccss_exp$agedx),
#                        cohort=c(rep("sjlife", nrow(sjlife)), rep("ccss_org", nrow(ccss_org)), rep("ccss_exp", nrow(ccss_exp))))
# cc <-mannwhitney(cc, age, cohort)
# cc$df$groups <- paste0(cc$df$grp1.label , " Vs ", cc$df$grp2.label)
# 
# # paste(cc$df$groups, cc$df$p, sep = "=")
# cc  = rbind.data.frame(cc$df$groups, cc$df$p)
# colnames(cc) <- cc[1,]
# agedx_p <- cc[-1,]
# 
# ## P-Sex
# cc <- cbind.data.frame(sjlife=table(sjlife$gender), ccss_org=table(ccss_org$gender), ccss_exp=table(ccss_exp$gender))
# cc <- as.table(rbind(M=as.numeric(cc[1,-1]), F=as.numeric(cc[2,-1])))
# colnames(cc) <- c("SJLIFE", "CCSS_org_hrc", "CCSS_exp_wgs")
# 
# sex_p <- chisq.test(cc)
# #################################
# ## Finalized demographic table ##
# #################################
# 
# SJLIFE_VARS <- c(CA=sjlife_CMP2plus_CA, CO=sjlife_CMP2plus_CO, AGE_AT_LAST_CONTACT = sjlife_agelstcontact, agedx_0_4 = sjlife_agedx_0_4, agedx_5_9 = sjlife_agedx_5_9, agedx_10_14 = sjlife_agedx_10_14, agedx_15orPLUS = sjlife_agedx_15_or_plus, MALE = sjlife_male, FEMALE = sjlife_female)
# CCSS_ORG_VARS <- c(CA=ccss_org_CMP2plus_CA, CO=ccss_org_CMP2plus_CO, AGE_AT_LAST_CONTACT = ccss_org_agelstcontact, agedx_0_4 = ccss_org_agedx_0_4, agedx_5_9 = ccss_org_agedx_5_9, agedx_10_14 = ccss_org_agedx_10_14, agedx_15orPLUS = ccss_org_agedx_15_or_plus, MALE = ccss_org_male, FEMALE = ccss_org_female)
# CCSS_EXP_VARS <- c(CA=ccss_exp_CMP2plus_CA, CO=ccss_exp_CMP2plus_CO, AGE_AT_LAST_CONTACT = ccss_exp_agelstcontact, agedx_0_4 = ccss_exp_agedx_0_4, agedx_5_9 = ccss_exp_agedx_5_9, agedx_10_14 = ccss_exp_agedx_10_14, agedx_15orPLUS = ccss_exp_agedx_15_or_plus, MALE = ccss_exp_male, FEMALE = ccss_exp_female)
# 
# final.df <- cbind.data.frame(SJLIFE_VARS, CCSS_ORG_VARS, CCSS_EXP_VARS)


