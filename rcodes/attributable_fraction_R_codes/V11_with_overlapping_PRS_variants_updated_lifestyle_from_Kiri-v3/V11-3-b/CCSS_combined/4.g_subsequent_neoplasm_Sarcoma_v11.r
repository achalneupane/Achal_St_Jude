setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/")
# CCSS.org.ANY_SN <- PHENO.ANY_SN # CCSS org
# CCSS.exp.ANY_SN <- PHENO.ANY_SN ## CCSS exp
# sum(colnames(CCSS.org.ANY_SN) == colnames(CCSS.exp.ANY_SN))
# # 54
# # Since columns are same, we can simply rbind the dataframes
# PHENO.ANY_SN <- rbind.data.frame(CCSS.org.ANY_SN, CCSS.exp.ANY_SN)
# rm(list=setdiff(ls(), c("PHENO.ANY_SN")))
# save.image("00.PHENO.ANY_SARCOMA_CCSS_combined_v11.Rdata")

load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.PHENO.ANY_SARCOMA_CCSS_combined_v11.Rdata")
source("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/attributable_fraction_R_codes/edit_lifestyle_variables.R")
PHENO.ANY_SN <- edit_lifestyle.ccss(PHENO.ANY_SN)

table(PHENO.ANY_SN$CACO)
# 0    1 
# 7839  104





##########################
dat_all = PHENO.ANY_SN
fit_all = glm(formula = SARCOMA ~ Sarcoma_Machiela_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                gender + aa_class_dose_5.category +
                # Current_smoker_yn + PhysicalActivity_yn + RiskyHeavyDrink_yn + Obese_yn +
                EAS + AFR,
              family = binomial,
              data = dat_all)


summary(fit_all)


##########################
## Get predicted values ##
##########################
dat_all$pred_all = predict(fit_all, newdat = dat_all, type = "response")

###############
## Treatment ##
###############

## Move relevant treatment exposures for everyone to no exposure
dat_tx = dat_all

dat_tx$aa_class_dose_5.category [!grepl("Unknown", dat_tx$aa_class_dose_5.category)] = "None"

dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
round(af_by_tx,3)
##################
## P/LP and PRS ##
##################
## P/LP Zhaoming, Qin without Zhaoming and PRS
dat_plp.prs = dat_all
# dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N";
dat_plp.prs$Sarcoma_Machiela_PRS.tertile.category = "1st"

dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs, na.rm = TRUE)
af_by_plp.prs = (N_all - N_no_plp.prs) / N_all
round(af_by_plp.prs,3)
###############
## Lifestyle ##
###############
dat_lifestyle = dat_all

dat_lifestyle$Current_smoker_yn [!grepl("Unknown", dat_lifestyle$Current_smoker_yn)] = "No"
dat_lifestyle$PhysicalActivity_yn [!grepl("Unknown", dat_lifestyle$PhysicalActivity_yn)] = "Yes"
dat_lifestyle$RiskyHeavyDrink_yn [!grepl("Unknown", dat_lifestyle$RiskyHeavyDrink_yn)] = "No"
# dat_lifestyle$HEALTHY_Diet_yn [!grepl("Unknown", dat_lifestyle$HEALTHY_Diet_yn)] = "Yes"
dat_lifestyle$Obese_yn [!grepl("Unknown", dat_lifestyle$Obese_yn)] = "No"


dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category = (N_all - N_no_favorable_lifestyle.category) / N_all
round(af_by_N_no_favorable_lifestyle.category,3)
#################################################
## Treatment, Genetics and Lifestyle, combined ##
#################################################

dat_tx.plp.prs.lifestyle = dat_all

## Nullify Treatment
dat_tx.plp.prs.lifestyle$aa_class_dose_5.category [!grepl("Unknown", dat_tx.plp.prs.lifestyle$aa_class_dose_5.category)] = "None"


## Nullify Genetics
# dat_tx.plp.prs.lifestyle$Zhaoming_carriers = dat_tx.plp.prs.lifestyle$Qin_without_Zhaoming_vars_carriers = "N";
dat_tx.plp.prs.lifestyle$Sarcoma_Machiela_PRS.tertile.category = "1st"

## Nullify Lifestyle
dat_tx.plp.prs.lifestyle$Current_smoker_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$Current_smoker_yn)] = "No"
dat_tx.plp.prs.lifestyle$PhysicalActivity_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$PhysicalActivity_yn)] = "Yes"
dat_tx.plp.prs.lifestyle$RiskyHeavyDrink_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$RiskyHeavyDrink_yn)] = "No"
# dat_tx.plp.prs.lifestyle$HEALTHY_Diet_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$HEALTHY_Diet_yn)] = "Yes"
dat_tx.plp.prs.lifestyle$Obese_yn [!grepl("Unknown", dat_tx.plp.prs.lifestyle$Obese_yn)] = "No"


dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_tx.plp.prs.lifestyle, type = "response")

N_no_favorable_tx.plp.prs.lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_tx.plp.prs.lifestyle.category = (N_all - N_no_favorable_tx.plp.prs.lifestyle.category) / N_all
round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3)
SARCOMA.res <- c(round(af_by_tx,3), round(af_by_plp.prs,3),round(af_by_N_no_favorable_lifestyle.category,3), round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3))
SARCOMA.res

all.res <- cbind.data.frame(SN=SN.res, SMN=SMN.res, NMSC=NMSC.res, BREAST=BREAST.res, THYROID=THYROID.res, MENINGIOMA=MENINGIOMA.res, SARCOMA=SARCOMA.res)
# View(all.res)
View(t(all.res))
