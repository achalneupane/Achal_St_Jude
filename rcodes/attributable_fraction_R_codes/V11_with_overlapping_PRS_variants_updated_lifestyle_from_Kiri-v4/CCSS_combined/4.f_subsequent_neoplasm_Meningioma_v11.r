setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/")
# CCSS.org.ANY_SN <- PHENO.ANY_SN # CCSS org
# CCSS.exp.ANY_SN <- PHENO.ANY_SN ## CCSS exp
# sum(colnames(CCSS.org.ANY_SN) == colnames(CCSS.exp.ANY_SN))
# # 53
# # Since columns are same, we can simply rbind the dataframes
# PHENO.ANY_SN <- rbind.data.frame(CCSS.org.ANY_SN, CCSS.exp.ANY_SN)
# rm(list=setdiff(ls(), c("PHENO.ANY_SN")))
# save.image("00.PHENO.ANY_MENINGIOMA_CCSS_combined_v11.Rdata")

load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/00.PHENO.ANY_MENINGIOMA_CCSS_combined_v11.Rdata")

table(PHENO.ANY_SN$CACO)
# 0    1 
# 7688  255
############################################################################################
############################################################################################
#################################### Work for V11-4 ########################################
############################################################################################
############################################################################################

PHENO.ANY_SN <- PHENO.ANY_SN[!(PHENO.ANY_SN$smoker_former_or_never_yn == "Unknown" &
                                 PHENO.ANY_SN$PhysicalActivity_yn == "Unknown" &
                                 PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn == "Unknown" &
                                 PHENO.ANY_SN$Not_obese_yn == "Unknown" ),]

dim(PHENO.ANY_SN)
# 7821   53

sum((PHENO.ANY_SN$smoker_former_or_never_yn_agesurvey >= 18|
       PHENO.ANY_SN$PhysicalActivity_yn_agesurvey >= 18|
       PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn_agesurvey >= 18|
       PHENO.ANY_SN$Not_obese_yn_agesurvey >= 18), na.rm = T)
# 7821


# ALL.LIFESTYLE.test <- ALL.LIFESTYLE


PHENO.ANY_SN <- PHENO.ANY_SN[which(PHENO.ANY_SN$smoker_former_or_never_yn_agesurvey >= 18 |
                                     PHENO.ANY_SN$PhysicalActivity_yn_agesurvey >= 18 |
                                     PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn_agesurvey >= 18 |
                                     PHENO.ANY_SN$Not_obese_yn_agesurvey >= 18),]


## SN diagnosis after survey age
# sum(((PHENO.ANY_SN$CACO == 1 & PHENO.ANY_SN$smoker_former_or_never_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)|
#   (PHENO.ANY_SN$CACO == 1 & PHENO.ANY_SN$PhysicalActivity_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)|
#   (PHENO.ANY_SN$CACO == 1 & PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)|
#   (PHENO.ANY_SN$CACO == 1 & PHENO.ANY_SN$Not_obese_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), na.rm = T)

sum(((PHENO.ANY_SN$CACO == 1 & PHENO.ANY_SN$smoker_former_or_never_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)&
       (PHENO.ANY_SN$CACO == 1 & PHENO.ANY_SN$PhysicalActivity_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)&
       (PHENO.ANY_SN$CACO == 1 & PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)&
       (PHENO.ANY_SN$CACO == 1 & PHENO.ANY_SN$Not_obese_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)), na.rm = T)
# 0

# PHENO.ANY_SN <- PHENO.ANY_SN[-which((PHENO.ANY_SN$CACO == 1 & PHENO.ANY_SN$smoker_former_or_never_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)|
#     (PHENO.ANY_SN$CACO == 1 & PHENO.ANY_SN$PhysicalActivity_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)|
#     (PHENO.ANY_SN$CACO == 1 & PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)|
#     (PHENO.ANY_SN$CACO == 1 & PHENO.ANY_SN$Not_obese_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)),]

PHENO.ANY_SN <- PHENO.ANY_SN[!((PHENO.ANY_SN$CACO == 1 & PHENO.ANY_SN$smoker_former_or_never_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)&
                                 (PHENO.ANY_SN$CACO == 1 & PHENO.ANY_SN$PhysicalActivity_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)&
                                 (PHENO.ANY_SN$CACO == 1 & PHENO.ANY_SN$NOT_RiskyHeavyDrink_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)&
                                 (PHENO.ANY_SN$CACO == 1 & PHENO.ANY_SN$Not_obese_yn_agesurvey > PHENO.ANY_SN$AGE.ANY_SN)),]

# PHENO.ANY_SN <- PHENO.ANY_SN[PHENO.ANY_SN$sjlid %in% PHENO.ANY_SN$SJLIFEID,]
table(PHENO.ANY_SN$MENINGIOMA)
# 0    1 
# 7606  215

##########################
dat_all = PHENO.ANY_SN
fit_all = glm(formula = MENINGIOMA ~ Meningioma_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + maxsegrtdose.category + epitxn_dose_5.category + 
                smoker_former_or_never_yn + PhysicalActivity_yn + NOT_RiskyHeavyDrink_yn + Not_obese_yn +
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

dat_tx$maxsegrtdose.category =
dat_tx$epitxn_dose_5.category = "None"

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
# dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N"
dat_plp.prs$dat_prs$Meningioma_PRS.tertile.category = "1st"

dat_all$pred_no_plp.prs = predict(fit_all, newdata = dat_plp.prs, type = "response")
N_no_plp.prs = sum(dat_all$pred_no_plp.prs, na.rm = TRUE)
af_by_plp.prs = (N_all - N_no_plp.prs) / N_all
round(af_by_plp.prs,3)
###############
## Lifestyle ##
###############
dat_lifestyle = dat_all

dat_lifestyle$smoker_former_or_never_yn =
dat_lifestyle$PhysicalActivity_yn =
dat_lifestyle$NOT_RiskyHeavyDrink_yn =
dat_lifestyle$Not_obese_yn = "1"

dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category = (N_all - N_no_favorable_lifestyle.category) / N_all
round(af_by_N_no_favorable_lifestyle.category,3)
#################################################
## Treatment, Genetics and Lifestyle, combined ##
#################################################
dat_tx.plp.prs.lifestyle = dat_all

## Nullify Treatment
dat_tx.plp.prs.lifestyle$maxsegrtdose.category =
  dat_tx.plp.prs.lifestyle$epitxn_dose_5.category = "None"


## Nullify Genetics
# dat_tx.plp.prs.lifestyle$Zhaoming_carriers = dat_tx.plp.prs.lifestyle$Qin_without_Zhaoming_vars_carriers = "N";
dat_tx.plp.prs.lifestyle$Meningioma_PRS.tertile.category = "1st"

## Nullify Lifestyle
dat_tx.plp.prs.lifestyle$smoker_former_or_never_yn =
dat_tx.plp.prs.lifestyle$PhysicalActivity_yn =
dat_tx.plp.prs.lifestyle$NOT_RiskyHeavyDrink_yn =
dat_tx.plp.prs.lifestyle$Not_obese_yn = "1"

dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_tx.plp.prs.lifestyle, type = "response")

N_no_favorable_tx.plp.prs.lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_tx.plp.prs.lifestyle.category = (N_all - N_no_favorable_tx.plp.prs.lifestyle.category) / N_all
round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3)
MENINGIOMA.res <- c(round(af_by_tx,3), round(af_by_plp.prs,3),round(af_by_N_no_favorable_lifestyle.category,3), round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3))
MENINGIOMA.res
# 0.870 0.000 0.348 0.937