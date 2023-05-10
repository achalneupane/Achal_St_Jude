# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/6.sjlife_without_diet.SARCOMA.V16.Rdata")

# Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
cc
filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
filtered_cc

# PHENO.ANY_SN$aa_class_dose_5.category <- as.character(PHENO.ANY_SN$aa_class_dose_5.category)
# PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "None"] <- "1st"
# PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "2nd"] <- "2nd-3rd"
# PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "3rd"] <- "2nd-3rd"
# PHENO.ANY_SN$aa_class_dose_5.category <- factor(PHENO.ANY_SN$aa_class_dose_5.category, levels = c("1st", "2nd-3rd"))

######################################
## Attributable fraction of Any SNs ##
######################################
dat_all = PHENO.ANY_SN

fit_all = glm(formula = SARCOMA ~ Sarcoma_Machiela_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                gender + 
                aa_class_dose_5.category +
                Current_smoker_yn + PhysicalActivity_yn + RiskyHeavyDrink_yn + Obese_yn +
                EAS + AFR + 
                any_lifestyle_missing + any_tx_missing,
              family = binomial,
              data = dat_all)

summary(fit_all)


(output <- summary(fit_all)$coefficients)
as.data.frame(apply(output, 2, formatC, format="f", digits=4))
# options(scipen=999)
estimate <- format(round(output[,1],3), nsmall = 3)
std.error <- format(round(output[,2],3), nsmall = 3)
# P.val <- formatC(output[,4], format="G", digits=3)
P.val <- output[,4]
P.val[P.val < 0.001] <- "<0.001"
P.val[!grepl("<", P.val)] <- format(round(as.numeric(P.val[!grepl("<", P.val)]), 3), nsmall = 3)
sn.model <- (setNames(cbind.data.frame(estimate, std.error, P.val
), c("Estimate", "Std.error", "P")))
sn.model <- sn.model[!grepl("AGE_AT_LAST_CONTACT", row.names(sn.model)),]
sn.model


##########################
## Get predicted values ##
##########################
dat_all$pred_all = predict(fit_all, newdat = dat_all, type = "response")

###############
## Treatment ##
###############

## Move relevant treatment exposures for everyone to no exposure
dat_tx = dat_all

dat_tx$any_tx_missing <- "No"

dat_tx$aa_class_dose_5.category = "None"

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

dat_lifestyle$any_lifestyle_missing <- "No"

dat_lifestyle$Current_smoker_yn = "No"
dat_lifestyle$PhysicalActivity_yn = "Yes"
dat_lifestyle$RiskyHeavyDrink_yn = "No"
# dat_lifestyle$HEALTHY_Diet_yn = "Yes"
dat_lifestyle$Obese_yn = "No"


dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category = (N_all - N_no_favorable_lifestyle.category) / N_all
round(af_by_N_no_favorable_lifestyle.category,3)
#################################################
## Treatment, Genetics and Lifestyle, combined ##
#################################################

dat_tx.plp.prs.lifestyle = dat_all

dat_tx.plp.prs.lifestyle$any_tx_missing <- "No"
dat_tx.plp.prs.lifestyle$any_lifestyle_missing <- "No"

## Nullify Treatment
dat_tx.plp.prs.lifestyle$aa_class_dose_5.category = "None"


## Nullify Genetics
# dat_tx.plp.prs.lifestyle$Zhaoming_carriers = dat_tx.plp.prs.lifestyle$Qin_without_Zhaoming_vars_carriers = "N";
dat_tx.plp.prs.lifestyle$Sarcoma_Machiela_PRS.tertile.category = "1st"

## Nullify Lifestyle
dat_tx.plp.prs.lifestyle$Current_smoker_yn = "No"
dat_tx.plp.prs.lifestyle$PhysicalActivity_yn = "Yes"
dat_tx.plp.prs.lifestyle$RiskyHeavyDrink_yn = "No"
# dat_tx.plp.prs.lifestyle$HEALTHY_Diet_yn = "Yes"
dat_tx.plp.prs.lifestyle$Obese_yn = "No"


dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_tx.plp.prs.lifestyle, type = "response")

N_no_favorable_tx.plp.prs.lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_tx.plp.prs.lifestyle.category = (N_all - N_no_favorable_tx.plp.prs.lifestyle.category) / N_all
round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3)
SARCOMA.res <- c(round(af_by_tx,3), round(af_by_plp.prs,3),round(af_by_N_no_favorable_lifestyle.category,3), round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3))
SARCOMA.res
# 0.373 0.023 0.541 0.744

all.res <- cbind.data.frame(SN=SN.res, SMN=SMN.res, NMSC=NMSC.res, BREAST=Breast.res, THYROID=Thyroid.res, MENINGIOMA=Meningioma.res)
# View(all.res)
View(all.res)
