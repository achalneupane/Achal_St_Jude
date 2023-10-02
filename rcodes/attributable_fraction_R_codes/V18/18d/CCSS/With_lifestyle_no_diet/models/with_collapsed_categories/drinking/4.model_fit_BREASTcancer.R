# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.BREASTcancer.V18b_without_diet.Rdata")

# Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
cc
filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
filtered_cc

# PHENO.ANY_SN$maxchestrtdose.category <- as.character(PHENO.ANY_SN$maxchestrtdose.category)
# PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$maxchestrtdose.category == ">0-<20"] <- "Any"
# PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$maxchestrtdose.category == ">=20"] <- "Any"
# PHENO.ANY_SN$maxchestrtdose.category <- factor(PHENO.ANY_SN$maxchestrtdose.category, levels = c("None", "Any"))

# table(PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$BREASTcancer == 1])

# PHENO.ANY_SN$anthra_jco_dose_5.category <- as.character(PHENO.ANY_SN$anthra_jco_dose_5.category)
# PHENO.ANY_SN$anthra_jco_dose_5.category[PHENO.ANY_SN$anthra_jco_dose_5.category == "1st"] <- "1st-2nd"
# PHENO.ANY_SN$anthra_jco_dose_5.category[PHENO.ANY_SN$anthra_jco_dose_5.category == "2nd"] <- "1st-2nd"
# PHENO.ANY_SN$anthra_jco_dose_5.category <- factor(PHENO.ANY_SN$anthra_jco_dose_5.category, levels = c("None", "1st-2nd", "3rd"))

# table(PHENO.ANY_SN$anthra_jco_dose_5.category[PHENO.ANY_SN$BREASTcancer == 1])

######################################
######################################
dat_all = PHENO.ANY_SN
dat_all=dat_all[dat_all$evt1==1,]

fit_all = glm(formula = event ~ 
                Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + 
                AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
                maxchestrtdose.category + anthra_jco_dose_5.category +
                Current_smoker_yn + PhysicalActivity_yn + RiskyHeavyDrink_yn + Obese_yn +
                EAS + AFR,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

##########################
##########################
dat_all$pred_all = predict(fit_all, newdat = dat_all, type = "response")


N_all = sum(dat_all$pred_all, na.rm = T) # Overall
N_all.male = sum(dat_all$pred_all[dat_all$gender == "Male"], na.rm = TRUE) # subset by gender
N_all.female = sum(dat_all$pred_all[dat_all$gender == "Female"], na.rm = TRUE) # subset by gender
# median(dat_all$AGE_AT_LAST_CONTACT.cs1)
N_all.lt.35 = sum(dat_all$pred_all[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE) # subset by age 35
N_all.gteq.35 = sum(dat_all$pred_all[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE) # subset by age 35


#############
#############
dat_tx = dat_all

dat_tx$any_chemo_missing <- "No" # **
dat_tx$anthra_jco_dose_5.category = "None" ## **

dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
af_by_tx <- round(af_by_tx,3)
af_by_tx

N_no_tx = sum(dat_all$pred_no_tx[dat_all$gender == "Male"], na.rm = TRUE)
af_by_tx.male = (N_all.male - N_no_tx) / N_all.male
af_by_tx.male <- round(af_by_tx.male,3)
af_by_tx.male

N_no_tx = sum(dat_all$pred_no_tx[dat_all$gender == "Female"], na.rm = TRUE)
af_by_tx.female = (N_all.female - N_no_tx) / N_all.female
af_by_tx.female <- round(af_by_tx.female,3)
af_by_tx.female

N_no_tx = sum(dat_all$pred_no_tx[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_tx.lt.35 = (N_all.lt.35 - N_no_tx) / N_all.lt.35
af_by_tx.lt.35 <- round(af_by_tx.lt.35,3)
af_by_tx.lt.35

N_no_tx = sum(dat_all$pred_no_tx[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_tx.gteq.35 = (N_all.gteq.35 - N_no_tx) / N_all.gteq.35
af_by_tx.gteq.35 <- round(af_by_tx.gteq.35,3)
af_by_tx.gteq.35

#############
#############

dat_rt = dat_all

dat_rt$any_rt_missing <- "No" # **

dat_rt$maxchestrtdose.category = "None" ## **

dat_all$pred_no_rt = predict(fit_all, newdata = dat_rt, type = "response")

N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_rt = sum(dat_all$pred_no_rt, na.rm = TRUE)
af_by_rt = (N_all - N_no_rt) / N_all
af_by_rt <- round(af_by_rt,3)
af_by_rt

N_no_rt = sum(dat_all$pred_no_rt[dat_all$gender == "Male"], na.rm = TRUE)
af_by_rt.male = (N_all.male - N_no_rt) / N_all.male
af_by_rt.male <- round(af_by_rt.male,3)
af_by_rt.male

N_no_rt = sum(dat_all$pred_no_rt[dat_all$gender == "Female"], na.rm = TRUE)
af_by_rt.female = (N_all.female - N_no_rt) / N_all.female
af_by_rt.female <- round(af_by_rt.female,3)
af_by_rt.female

N_no_rt = sum(dat_all$pred_no_rt[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_rt.lt.35 = (N_all.lt.35 - N_no_rt) / N_all.lt.35
af_by_rt.lt.35 <- round(af_by_rt.lt.35,3)
af_by_rt.lt.35

N_no_rt = sum(dat_all$pred_no_rt[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_rt.gteq.35 = (N_all.gteq.35 - N_no_rt) / N_all.gteq.35
af_by_rt.gteq.35 <- round(af_by_rt.gteq.35,3)
af_by_rt.gteq.35

######################
######################

dat_tx.rt = dat_all

dat_tx.rt$any_chemo_missing <- "No" ## **
dat_tx.rt$any_rt_missing <- "No" ## **

dat_tx.rt$maxchestrtdose.category =
  dat_tx.rt$anthra_jco_dose_5.category = "None" ## **

dat_all$pred_no_tx.rt = predict(fit_all, newdata = dat_tx.rt, type = "response")

N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx.rt = sum(dat_all$pred_no_tx.rt, na.rm = TRUE)
af_by_tx.rt = (N_all - N_no_tx.rt) / N_all
af_by_tx.rt <- round(af_by_tx.rt,3)
af_by_tx.rt

N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$gender == "Male"], na.rm = TRUE)
af_by_tx.rt.male = (N_all.male - N_no_tx.rt) / N_all.male
af_by_tx.rt.male <- round(af_by_tx.rt.male,3)
af_by_tx.rt.male

N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$gender == "Female"], na.rm = TRUE)
af_by_tx.rt.female = (N_all.female - N_no_tx.rt) / N_all.female
af_by_tx.rt.female <- round(af_by_tx.rt.female,3)
af_by_tx.rt.female

N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_tx.rt.lt.35 = (N_all.lt.35 - N_no_tx.rt) / N_all.lt.35
af_by_tx.rt.lt.35 <- round(af_by_tx.rt.lt.35,3)
af_by_tx.rt.lt.35

N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_tx.rt.gteq.35 = (N_all.gteq.35 - N_no_tx.rt) / N_all.gteq.35
af_by_tx.rt.gteq.35 <- round(af_by_tx.rt.gteq.35,3)
af_by_tx.rt.gteq.35
#########
#########
dat_prs = dat_all
# dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N";
dat_prs$Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category = "1st"  # **

dat_all$pred_no_prs = predict(fit_all, newdata = dat_prs, type = "response")
N_no_prs = sum(dat_all$pred_no_prs, na.rm = TRUE)
af_by_prs = (N_all - N_no_prs) / N_all
af_by_prs <- round(af_by_prs,3)
af_by_prs

N_no_prs = sum(dat_all$pred_no_prs[dat_all$gender == "Male"], na.rm = TRUE)
af_by_prs.male = (N_all.male - N_no_prs) / N_all.male
af_by_prs.male <- round(af_by_prs.male,3)
af_by_prs.male

N_no_prs = sum(dat_all$pred_no_prs[dat_all$gender == "Female"], na.rm = TRUE)
af_by_prs.female = (N_all.female - N_no_prs) / N_all.female
af_by_prs.female <- round(af_by_prs.female,3)
af_by_prs.female

N_no_prs = sum(dat_all$pred_no_prs[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_prs.lt.35 = (N_all.lt.35 - N_no_prs) / N_all.lt.35
af_by_prs.lt.35 <- round(af_by_prs.lt.35,3)
af_by_prs.lt.35

N_no_prs = sum(dat_all$pred_no_prs[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_prs.gteq.35 = (N_all.gteq.35 - N_no_prs) / N_all.gteq.35
af_by_prs.gteq.35 <- round(af_by_prs.gteq.35,3)
af_by_prs.gteq.35
###############
###############
dat_lifestyle = dat_all

dat_lifestyle$any_lifestyle_missing <- "No"

# dat_lifestyle$Current_smoker_yn = "No"
# dat_lifestyle$PhysicalActivity_yn = "Yes"
dat_lifestyle$RiskyHeavyDrink_yn = "No"
# dat_lifestyle$Obese_yn = "No"

dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_no_favorable_lifestyle.category = (N_all - N_no_favorable_lifestyle.category) / N_all
af_by_no_favorable_lifestyle.category <- round(af_by_no_favorable_lifestyle.category,3)
af_by_no_favorable_lifestyle.category

N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category[dat_all$gender == "Male"], na.rm = TRUE)
af_by_no_favorable_lifestyle.category.male = (N_all.male - N_no_favorable_lifestyle.category) / N_all.male
af_by_no_favorable_lifestyle.category.male <- round(af_by_no_favorable_lifestyle.category.male,3)
af_by_no_favorable_lifestyle.category.male

N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category[dat_all$gender == "Female"], na.rm = TRUE)
af_by_no_favorable_lifestyle.category.female = (N_all.female - N_no_favorable_lifestyle.category) / N_all.female
af_by_no_favorable_lifestyle.category.female <- round(af_by_no_favorable_lifestyle.category.female,3)
af_by_no_favorable_lifestyle.category.female

N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_no_favorable_lifestyle.category.lt.35 = (N_all.lt.35 - N_no_favorable_lifestyle.category) / N_all.lt.35
af_by_no_favorable_lifestyle.category.lt.35 <- round(af_by_no_favorable_lifestyle.category.lt.35,3)
af_by_no_favorable_lifestyle.category.lt.35

N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_no_favorable_lifestyle.category.gteq.35 = (N_all.gteq.35 - N_no_favorable_lifestyle.category) / N_all.gteq.35
af_by_no_favorable_lifestyle.category.gteq.35 <- round(af_by_no_favorable_lifestyle.category.gteq.35,3)
af_by_no_favorable_lifestyle.category.gteq.35

#################################################
#################################################
dat_tx.prs.lifestyle = dat_all

dat_tx.prs.lifestyle$any_chemo_missing <- "No" ## **
dat_tx.prs.lifestyle$any_rt_missing <- "No" ## **

dat_tx.prs.lifestyle$any_lifestyle_missing <- "No"

dat_tx.prs.lifestyle$maxchestrtdose.category =
  dat_tx.prs.lifestyle$anthra_jco_dose_5.category = "None" ## **

# dat_tx.plp.prs.lifestyle$Zhaoming_carriers = dat_tx.plp.prs.lifestyle$Qin_without_Zhaoming_vars_carriers = "N";
dat_tx.prs.lifestyle$Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category = "1st" ## **

# dat_tx.prs.lifestyle$Current_smoker_yn = "No"
# dat_tx.prs.lifestyle$PhysicalActivity_yn = "Yes"
dat_tx.prs.lifestyle$RiskyHeavyDrink_yn = "No"
# dat_tx.prs.lifestyle$Obese_yn = "No"


dat_all$pred_no_combined = predict(fit_all, newdata = dat_tx.prs.lifestyle, type = "response")

N_no_combined = sum(dat_all$pred_no_combined, na.rm = TRUE)
af_by_combined = (N_all - N_no_combined) / N_all
af_by_combined <- round(af_by_combined,3)
af_by_combined


N_no_combined = sum(dat_all$pred_no_combined[dat_all$gender == "Male"], na.rm = TRUE)
af_by_combined.male = (N_all.male - N_no_combined) / N_all.male
af_by_combined.male <- round(af_by_combined.male,3)
af_by_combined.male

N_no_combined = sum(dat_all$pred_no_combined[dat_all$gender == "Female"], na.rm = TRUE)
af_by_combined.female = (N_all.female - N_no_combined) / N_all.female
af_by_combined.female <- round(af_by_combined.female,3)
af_by_combined.female

N_no_combined = sum(dat_all$pred_no_combined[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_combined.lt.35 = (N_all.lt.35 - N_no_combined) / N_all.lt.35
af_by_combined.lt.35 <- round(af_by_combined.lt.35,3)
af_by_combined.lt.35

N_no_combined = sum(dat_all$pred_no_combined[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_combined.gteq.35 = (N_all.gteq.35 - N_no_combined) / N_all.gteq.35
af_by_combined.gteq.35 <- round(af_by_combined.gteq.35,3)
af_by_combined.gteq.35

##
Breast.res <- data.frame(
  Variable = c("Radiation", "Chemo", "All treatments", "PRS", "Lifestyle", "Combined"),
  Overall = c(af_by_rt, af_by_tx, af_by_tx.rt, af_by_prs, af_by_no_favorable_lifestyle.category, af_by_combined),
  Female = c(af_by_rt.female, af_by_tx.female, af_by_tx.rt.female, af_by_prs.female, af_by_no_favorable_lifestyle.category.female, af_by_combined.female),
  Male = c(af_by_rt.male, af_by_tx.male, af_by_tx.rt.male, af_by_prs.male, af_by_no_favorable_lifestyle.category.male, af_by_combined.male),
  age.lt35 = c(af_by_rt.lt.35, af_by_tx.lt.35, af_by_tx.rt.lt.35, af_by_prs.lt.35, af_by_no_favorable_lifestyle.category.lt.35, af_by_combined.lt.35),
  age.gteq = c(af_by_rt.gteq.35, af_by_tx.gteq.35, af_by_tx.rt.gteq.35, af_by_prs.gteq.35, af_by_no_favorable_lifestyle.category.gteq.35, af_by_combined.gteq.35)
)

# View(Breast.res)
# 
# 
# 
# 
# #########################################
# ## Check PRS and treatment interaction ##
# #########################################
# dat_all = PHENO.ANY_SN
# dat_all=dat_all[dat_all$evt1==1,]
# 
# fit_all = glm(formula = event ~ 
#                 Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category +
#                 AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + 
#                 AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
#                 maxchestrtdose.category + anthra_jco_dose_5.category +
#                 Current_smoker_yn + PhysicalActivity_yn + RiskyHeavyDrink_yn + Obese_yn +
#                 EAS + AFR,
#                 maxchestrtdose.category * Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category +
#                 anthra_jco_dose_5.category * Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category,
#               family = "poisson", offset = log(dat_all$PY), data = dat_all)
# 
# summary(fit_all)
# 
# (output <- summary(fit_all)$coefficients)
# as.data.frame(apply(output, 2, formatC, format="f", digits=4))
# estimate <- format(round(output[,1],3), nsmall = 3)
# std.error <- format(round(output[,2],3), nsmall = 3)
# P.val <- output[,4]
# P.val[P.val < 0.001] <- "<0.001"
# P.val[!grepl("<", P.val)] <- format(round(as.numeric(P.val[!grepl("<", P.val)]), 3), nsmall = 3)
# sn.model <- (setNames(cbind.data.frame(estimate, std.error, P.val
# ), c("Estimate", "Std.error", "P")))
# sn.model <- sn.model[!grepl("AGE_AT_LAST_CONTACT", row.names(sn.model)),]
# sn.model