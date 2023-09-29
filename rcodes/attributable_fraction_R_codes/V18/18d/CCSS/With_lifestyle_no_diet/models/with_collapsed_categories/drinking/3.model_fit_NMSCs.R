# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.NMSCs.V18b_without_diet.Rdata")

# Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
cc
filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
filtered_cc

# PHENO.ANY_SN$maxsegrtdose.category <- as.character(PHENO.ANY_SN$maxsegrtdose.category)
# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">0-<18"] <- ">0-<30"
# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">=18-<30"] <- ">0-<30"
# PHENO.ANY_SN$maxsegrtdose.category <- factor(PHENO.ANY_SN$maxsegrtdose.category, levels = c("None", ">0-<30", ">=30"))

# PHENO.ANY_SN$maxpelvisrtdose.category <- as.character(PHENO.ANY_SN$maxpelvisrtdose.category)
# PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$maxpelvisrtdose.category == ">0-<20"] <- "Any"
# PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$maxpelvisrtdose.category == ">=20"] <- "Any"
# PHENO.ANY_SN$maxpelvisrtdose.category <- factor(PHENO.ANY_SN$maxpelvisrtdose.category, levels = c("None", "Any"))

# table(PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$NMSCs == 0])

######################################
## Attributable fraction of Any SNs ##
######################################
dat_all = PHENO.ANY_SN

dat_all=dat_all[dat_all$evt1==1,]
fit_all = glm(formula = event ~ BASALcell_PRS.tertile.category + 
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category +
                Current_smoker_yn + PhysicalActivity_yn + RiskyHeavyDrink_yn + Obese_yn +
                EAS + AFR,
                any_lifestyle_missing + any_rt_missing,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)


##########################
## Get predicted values ##
##########################
dat_all$pred_all = predict(fit_all, newdat = dat_all, type = "response")


## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = T) # Overall
N_all.male = sum(dat_all$pred_all[dat_all$gender == "Male"], na.rm = TRUE) # subset by gender
N_all.female = sum(dat_all$pred_all[dat_all$gender == "Female"], na.rm = TRUE) # subset by gender
## subset by age at diagnosis group
# median(dat_all$AGE_AT_LAST_CONTACT.cs1)
N_all.lt.35 = sum(dat_all$pred_all[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE) # subset by age 35
N_all.gteq.35 = sum(dat_all$pred_all[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE) # subset by age 35


#############
## tx only ##
#############
## There is no chemo in NMSC model
## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
af_by_tx <- "-"

## Male
af_by_tx.male <- "-"

## Female
af_by_tx.female <- "-"

## < 35
af_by_tx.lt.35 <- "-"

## >= 35
af_by_tx.gteq.35 <- "-"

#############
## RT only ##
#############

## Move relevant treatment exposures for everyone to no exposure
dat_rt = dat_all

dat_rt$any_rt_missing <- "No" # **


dat_rt$maxsegrtdose.category =
  dat_rt$maxabdrtdose.category =
  dat_rt$maxpelvisrtdose.category = "None" ## **

dat_all$pred_no_rt = predict(fit_all, newdata = dat_rt, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_rt = sum(dat_all$pred_no_rt, na.rm = TRUE)
af_by_rt = (N_all - N_no_rt) / N_all
af_by_rt <- round(af_by_rt,3)
af_by_rt

## Male
N_no_rt = sum(dat_all$pred_no_rt[dat_all$gender == "Male"], na.rm = TRUE)
af_by_rt.male = (N_all.male - N_no_rt) / N_all.male
af_by_rt.male <- round(af_by_rt.male,3)
af_by_rt.male

## Female
N_no_rt = sum(dat_all$pred_no_rt[dat_all$gender == "Female"], na.rm = TRUE)
af_by_rt.female = (N_all.female - N_no_rt) / N_all.female
af_by_rt.female <- round(af_by_rt.female,3)
af_by_rt.female

## < 35
N_no_rt = sum(dat_all$pred_no_rt[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_rt.lt.35 = (N_all.lt.35 - N_no_rt) / N_all.lt.35
af_by_rt.lt.35 <- round(af_by_rt.lt.35,3)
af_by_rt.lt.35

## >= 35
N_no_rt = sum(dat_all$pred_no_rt[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_rt.gteq.35 = (N_all.gteq.35 - N_no_rt) / N_all.gteq.35
af_by_rt.gteq.35 <- round(af_by_rt.gteq.35,3)
af_by_rt.gteq.35

######################
## Treatment and RT ##
######################

## Move relevant treatment exposures for everyone to no exposure
dat_tx.rt = dat_all

dat_tx.rt$any_rt_missing <- "No" ## **

dat_tx.rt$maxsegrtdose.category =
  dat_tx.rt$maxabdrtdose.category =
  dat_tx.rt$maxpelvisrtdose.category = "None" ## **

dat_all$pred_no_tx.rt = predict(fit_all, newdata = dat_tx.rt, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx.rt = sum(dat_all$pred_no_tx.rt, na.rm = TRUE)
af_by_tx.rt = (N_all - N_no_tx.rt) / N_all
af_by_tx.rt <- round(af_by_tx.rt,3)
af_by_tx.rt

## Male
N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$gender == "Male"], na.rm = TRUE)
af_by_tx.rt.male = (N_all.male - N_no_tx.rt) / N_all.male
af_by_tx.rt.male <- round(af_by_tx.rt.male,3)
af_by_tx.rt.male

## Female
N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$gender == "Female"], na.rm = TRUE)
af_by_tx.rt.female = (N_all.female - N_no_tx.rt) / N_all.female
af_by_tx.rt.female <- round(af_by_tx.rt.female,3)
af_by_tx.rt.female

## < 35
N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_tx.rt.lt.35 = (N_all.lt.35 - N_no_tx.rt) / N_all.lt.35
af_by_tx.rt.lt.35 <- round(af_by_tx.rt.lt.35,3)
af_by_tx.rt.lt.35

## >= 35
N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_tx.rt.gteq.35 = (N_all.gteq.35 - N_no_tx.rt) / N_all.gteq.35
af_by_tx.rt.gteq.35 <- round(af_by_tx.rt.gteq.35,3)
af_by_tx.rt.gteq.35
#########
## PRS ##
#########
dat_prs = dat_all
# dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N";
dat_prs$BASALcell_PRS.tertile.category = "1st"  # **

dat_all$pred_no_prs = predict(fit_all, newdata = dat_prs, type = "response")
N_no_prs = sum(dat_all$pred_no_prs, na.rm = TRUE)
af_by_prs = (N_all - N_no_prs) / N_all
af_by_prs <- round(af_by_prs,3)
af_by_prs

## Male
N_no_prs = sum(dat_all$pred_no_prs[dat_all$gender == "Male"], na.rm = TRUE)
af_by_prs.male = (N_all.male - N_no_prs) / N_all.male
af_by_prs.male <- round(af_by_prs.male,3)
af_by_prs.male

## Female
N_no_prs = sum(dat_all$pred_no_prs[dat_all$gender == "Female"], na.rm = TRUE)
af_by_prs.female = (N_all.female - N_no_prs) / N_all.female
af_by_prs.female <- round(af_by_prs.female,3)
af_by_prs.female

## < 35
N_no_prs = sum(dat_all$pred_no_prs[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_prs.lt.35 = (N_all.lt.35 - N_no_prs) / N_all.lt.35
af_by_prs.lt.35 <- round(af_by_prs.lt.35,3)
af_by_prs.lt.35

## >= 35
N_no_prs = sum(dat_all$pred_no_prs[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_prs.gteq.35 = (N_all.gteq.35 - N_no_prs) / N_all.gteq.35
af_by_prs.gteq.35 <- round(af_by_prs.gteq.35,3)
af_by_prs.gteq.35
###############
## Lifestyle ##
###############
dat_lifestyle = dat_all

dat_lifestyle$any_lifestyle_missing <- "No"

# dat_lifestyle$Current_smoker_yn = "No"
# dat_lifestyle$PhysicalActivity_yn = "Yes"
dat_lifestyle$RiskyHeavyDrink_yn = "No"
# # dat_lifestyle$HEALTHY_Diet_yn = "Yes"
# dat_lifestyle$Obese_yn = "No"

dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_no_favorable_lifestyle.category = (N_all - N_no_favorable_lifestyle.category) / N_all
af_by_no_favorable_lifestyle.category <- round(af_by_no_favorable_lifestyle.category,3)


## Male
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category[dat_all$gender == "Male"], na.rm = TRUE)
af_by_no_favorable_lifestyle.category.male = (N_all.male - N_no_favorable_lifestyle.category) / N_all.male
af_by_no_favorable_lifestyle.category.male <- round(af_by_no_favorable_lifestyle.category.male,3)
af_by_no_favorable_lifestyle.category.male

## Female
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category[dat_all$gender == "Female"], na.rm = TRUE)
af_by_no_favorable_lifestyle.category.female = (N_all.female - N_no_favorable_lifestyle.category) / N_all.female
af_by_no_favorable_lifestyle.category.female <- round(af_by_no_favorable_lifestyle.category.female,3)
af_by_no_favorable_lifestyle.category.female

## < 35
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_no_favorable_lifestyle.category.lt.35 = (N_all.lt.35 - N_no_favorable_lifestyle.category) / N_all.lt.35
af_by_no_favorable_lifestyle.category.lt.35 <- round(af_by_no_favorable_lifestyle.category.lt.35,3)
af_by_no_favorable_lifestyle.category.lt.35

## >= 35
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_no_favorable_lifestyle.category.gteq.35 = (N_all.gteq.35 - N_no_favorable_lifestyle.category) / N_all.gteq.35
af_by_no_favorable_lifestyle.category.gteq.35 <- round(af_by_no_favorable_lifestyle.category.gteq.35,3)
af_by_no_favorable_lifestyle.category.gteq.35

#################################################
## Treatment, Genetics and Lifestyle, combined ##
#################################################

dat_tx.prs.lifestyle = dat_all

dat_tx.prs.lifestyle$any_rt_missing <- "No" ## **

dat_tx.prs.lifestyle$any_lifestyle_missing <- "No"

## Nullify Treatment
dat_tx.prs.lifestyle$maxsegrtdose.category =
  dat_tx.prs.lifestyle$maxabdrtdose.category =
  dat_tx.prs.lifestyle$maxpelvisrtdose.category = "None" ## **

## Nullify Genetics
# dat_tx.plp.prs.lifestyle$Zhaoming_carriers = dat_tx.plp.prs.lifestyle$Qin_without_Zhaoming_vars_carriers = "N";
dat_tx.prs.lifestyle$BASALcell_PRS.tertile.category = "1st" ## **

## Nullify Lifestyle
# dat_tx.prs.lifestyle$Current_smoker_yn = "No"
# dat_tx.prs.lifestyle$PhysicalActivity_yn = "Yes"
dat_tx.prs.lifestyle$RiskyHeavyDrink_yn = "No"
# # dat_tx.prs.lifestyle$HEALTHY_Diet_yn = "Yes"
# dat_tx.prs.lifestyle$Obese_yn = "No"


dat_all$pred_no_combined = predict(fit_all, newdata = dat_tx.prs.lifestyle, type = "response")

N_no_combined = sum(dat_all$pred_no_combined, na.rm = TRUE)
af_by_combined = (N_all - N_no_combined) / N_all
af_by_combined <- round(af_by_combined,3)
af_by_combined


## Male
N_no_combined = sum(dat_all$pred_no_combined[dat_all$gender == "Male"], na.rm = TRUE)
af_by_combined.male = (N_all.male - N_no_combined) / N_all.male
af_by_combined.male <- round(af_by_combined.male,3)
af_by_combined.male

## Female
N_no_combined = sum(dat_all$pred_no_combined[dat_all$gender == "Female"], na.rm = TRUE)
af_by_combined.female = (N_all.female - N_no_combined) / N_all.female
af_by_combined.female <- round(af_by_combined.female,3)
af_by_combined.female

## < 35
N_no_combined = sum(dat_all$pred_no_combined[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_combined.lt.35 = (N_all.lt.35 - N_no_combined) / N_all.lt.35
af_by_combined.lt.35 <- round(af_by_combined.lt.35,3)
af_by_combined.lt.35

## >= 35
N_no_combined = sum(dat_all$pred_no_combined[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_combined.gteq.35 = (N_all.gteq.35 - N_no_combined) / N_all.gteq.35
af_by_combined.gteq.35 <- round(af_by_combined.gteq.35,3)
af_by_combined.gteq.35

##
NMSC.res <- data.frame(
  Variable = c("Radiation", "Chemo", "All treatments", "PRS", "Lifestyle", "Combined"),
  Overall = c(af_by_rt, af_by_tx, af_by_tx.rt, af_by_prs, af_by_no_favorable_lifestyle.category, af_by_combined),
  Female = c(af_by_rt.female, af_by_tx.female, af_by_tx.rt.female, af_by_prs.female, af_by_no_favorable_lifestyle.category.female, af_by_combined.female),
  Male = c(af_by_rt.male, af_by_tx.male, af_by_tx.rt.male, af_by_prs.male, af_by_no_favorable_lifestyle.category.male, af_by_combined.male),
  age.lt35 = c(af_by_rt.lt.35, af_by_tx.lt.35, af_by_tx.rt.lt.35, af_by_prs.lt.35, af_by_no_favorable_lifestyle.category.lt.35, af_by_combined.lt.35),
  age.gteq = c(af_by_rt.gteq.35, af_by_tx.gteq.35, af_by_tx.rt.gteq.35, af_by_prs.gteq.35, af_by_no_favorable_lifestyle.category.gteq.35, af_by_combined.gteq.35)
)

# View(NMSC.res)
# 
# 
# #########################################
# ## Check PRS and treatment interaction ##
# #########################################
# dat_all=PHENO.ANY_SN[PHENO.ANY_SN$evt1==1,]
# fit_all = glm(formula = event ~ BASALcell_PRS.tertile.category + 
#                 AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
#                 gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category +
#                 Current_smoker_yn + PhysicalActivity_yn + RiskyHeavyDrink_yn + Obese_yn +
#                 EAS + AFR,
#                 any_lifestyle_missing + any_rt_missing +
#                 maxsegrtdose.category*BASALcell_PRS.tertile.category + 
#                 maxabdrtdose.category*BASALcell_PRS.tertile.category + 
#                 maxpelvisrtdose.category*BASALcell_PRS.tertile.category, 
#               family = "poisson", offset = log(dat_all$PY), data = dat_all)
# 
# summary(fit_all)
# 
# (output <- summary(fit_all)$coefficients)
# as.data.frame(apply(output, 2, formatC, format="f", digits=4))
# # options(scipen=999)
# estimate <- format(round(output[,1],3), nsmall = 3)
# std.error <- format(round(output[,2],3), nsmall = 3)
# # P.val <- formatC(output[,4], format="G", digits=3)
# P.val <- output[,4]
# P.val[P.val < 0.001] <- "<0.001"
# P.val[!grepl("<", P.val)] <- format(round(as.numeric(P.val[!grepl("<", P.val)]), 3), nsmall = 3)
# sn.model <- (setNames(cbind.data.frame(estimate, std.error, P.val
# ), c("Estimate", "Std.error", "P")))
# sn.model <- sn.model[!grepl("AGE_AT_LAST_CONTACT", row.names(sn.model)),]
# sn.model
# View(sn.model)