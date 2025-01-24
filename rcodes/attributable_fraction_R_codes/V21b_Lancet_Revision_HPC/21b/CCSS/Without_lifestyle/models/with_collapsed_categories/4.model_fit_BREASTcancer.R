# load ANY SN data
# load("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.BREASTcancer_without_lifestyle.V20b.Rdata")
load("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.BREASTcancer_without_lifestyle.V20b.Rdata")

# Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
cc
filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
filtered_cc

## Admixture classification
admixture <- read.table("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
admixture$INDIVIDUAL <- sapply(strsplit(admixture$INDIVIDUAL,"_"), `[`, 1)
EUR.admix <- admixture$INDIVIDUAL[admixture$EUR > 0.8]
AFR.admix <- admixture$INDIVIDUAL[admixture$AFR > 0.6]

PHENO.ANY_SN$admixture <- NA
PHENO.ANY_SN$admixture [PHENO.ANY_SN$ccssid %in% EUR.admix] <- "EUR"
PHENO.ANY_SN$admixture [PHENO.ANY_SN$ccssid %in% AFR.admix] <- "AFR"

PHENO.ANY_SN$maxchestrtdose.category <- as.character(PHENO.ANY_SN$maxchestrtdose.category)
PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$maxchestrtdose.category == ">0-<20"] <- "Any"
PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$maxchestrtdose.category == ">=20"] <- "Any"
PHENO.ANY_SN$maxchestrtdose.category <- factor(PHENO.ANY_SN$maxchestrtdose.category, levels = c("None", "Any"))

# table(PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$BREASTcancer == 1])

PHENO.ANY_SN$anthra_jco_dose_5.category <- as.character(PHENO.ANY_SN$anthra_jco_dose_5.category)
PHENO.ANY_SN$anthra_jco_dose_5.category[PHENO.ANY_SN$anthra_jco_dose_5.category == "1st"] <- "1st-2nd"
PHENO.ANY_SN$anthra_jco_dose_5.category[PHENO.ANY_SN$anthra_jco_dose_5.category == "2nd"] <- "1st-2nd"
PHENO.ANY_SN$anthra_jco_dose_5.category <- factor(PHENO.ANY_SN$anthra_jco_dose_5.category, levels = c("None", "1st-2nd", "3rd"))

# table(PHENO.ANY_SN$anthra_jco_dose_5.category[PHENO.ANY_SN$BREASTcancer == 1])
######################################
## Attributable fraction of Any SNs ##
######################################
dat_all = PHENO.ANY_SN
dat_all = dat_all[dat_all$gender == "Female",]
dat_all=dat_all[dat_all$evt1==1,]

fit_all = glm(formula = event ~ 
                Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + 
                AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
                maxchestrtdose.category + anthra_jco_dose_5.category +
                EAS + AFR +
                any_chemo_missing + any_rt_missing,
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
## Subset by ancestry
N_all.EUR = sum(dat_all$pred_all[dat_all$admixture == "EUR"], na.rm = TRUE) # subset by ancestry
N_all.AFR = sum(dat_all$pred_all[dat_all$admixture == "AFR"], na.rm = TRUE) # subset by ancestry


#############
## tx only ##
#############
## Move relevant treatment exposures for everyone to no exposure
dat_tx = dat_all

dat_tx$any_chemo_missing <- "No" # **
dat_tx$anthra_jco_dose_5.category = "None" ## **

dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
# af_by_tx <- round(af_by_tx,3)
af_by_tx

## Male
N_no_tx = sum(dat_all$pred_no_tx[dat_all$gender == "Male"], na.rm = TRUE)
af_by_tx.male = (N_all.male - N_no_tx) / N_all.male
# af_by_tx.male <- round(af_by_tx.male,3)
af_by_tx.male

## Female
N_no_tx = sum(dat_all$pred_no_tx[dat_all$gender == "Female"], na.rm = TRUE)
af_by_tx.female = (N_all.female - N_no_tx) / N_all.female
# af_by_tx.female <- round(af_by_tx.female,3)
af_by_tx.female

## < 35
N_no_tx = sum(dat_all$pred_no_tx[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_tx.lt.35 = (N_all.lt.35 - N_no_tx) / N_all.lt.35
# af_by_tx.lt.35 <- round(af_by_tx.lt.35,3)
af_by_tx.lt.35

## >= 35
N_no_tx = sum(dat_all$pred_no_tx[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_tx.gteq.35 = (N_all.gteq.35 - N_no_tx) / N_all.gteq.35
# af_by_tx.gteq.35 <- round(af_by_tx.gteq.35,3)
af_by_tx.gteq.35

## EUR
N_no_tx = sum(dat_all$pred_no_tx[dat_all$admixture == "EUR"], na.rm = TRUE)
af_by_tx.EUR = (N_all.EUR - N_no_tx) / N_all.EUR
# af_by_tx.EUR <- round(af_by_tx.EUR,3)
af_by_tx.EUR

## AFR
N_no_tx = sum(dat_all$pred_no_tx[dat_all$admixture == "AFR"], na.rm = TRUE)
af_by_tx.AFR = (N_all.AFR - N_no_tx) / N_all.AFR
# af_by_tx.AFR <- round(af_by_tx.AFR,3)
af_by_tx.AFR

#############
## RT only ##
#############

## Move relevant treatment exposures for everyone to no exposure
dat_rt = dat_all

dat_rt$any_rt_missing <- "No" # **

dat_rt$maxchestrtdose.category = "None" ## **

dat_all$pred_no_rt = predict(fit_all, newdata = dat_rt, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_rt = sum(dat_all$pred_no_rt, na.rm = TRUE)
af_by_rt = (N_all - N_no_rt) / N_all
# af_by_rt <- round(af_by_rt,3)
af_by_rt

## Male
N_no_rt = sum(dat_all$pred_no_rt[dat_all$gender == "Male"], na.rm = TRUE)
af_by_rt.male = (N_all.male - N_no_rt) / N_all.male
# af_by_rt.male <- round(af_by_rt.male,3)
af_by_rt.male

## Female
N_no_rt = sum(dat_all$pred_no_rt[dat_all$gender == "Female"], na.rm = TRUE)
af_by_rt.female = (N_all.female - N_no_rt) / N_all.female
# af_by_rt.female <- round(af_by_rt.female,3)
af_by_rt.female

## < 35
N_no_rt = sum(dat_all$pred_no_rt[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_rt.lt.35 = (N_all.lt.35 - N_no_rt) / N_all.lt.35
# af_by_rt.lt.35 <- round(af_by_rt.lt.35,3)
af_by_rt.lt.35

## >= 35
N_no_rt = sum(dat_all$pred_no_rt[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_rt.gteq.35 = (N_all.gteq.35 - N_no_rt) / N_all.gteq.35
# af_by_rt.gteq.35 <- round(af_by_rt.gteq.35,3)
af_by_rt.gteq.35

## EUR
N_no_rt = sum(dat_all$pred_no_rt[dat_all$admixture == "EUR"], na.rm = TRUE)
af_by_rt.EUR = (N_all.EUR - N_no_rt) / N_all.EUR
# af_by_rt.EUR <- round(af_by_rt.EUR,3)
af_by_rt.EUR

## AFR
N_no_rt = sum(dat_all$pred_no_rt[dat_all$admixture == "AFR"], na.rm = TRUE)
af_by_rt.AFR = (N_all.AFR - N_no_rt) / N_all.AFR
# af_by_rt.AFR <- round(af_by_rt.AFR,3)
af_by_rt.AFR

######################
## Treatment and RT ##
######################

## Move relevant treatment exposures for everyone to no exposure
dat_tx.rt = dat_all

dat_tx.rt$any_chemo_missing <- "No" ## **
dat_tx.rt$any_rt_missing <- "No" ## **

dat_tx.rt$maxchestrtdose.category =
  dat_tx.rt$anthra_jco_dose_5.category = "None" ## **

dat_all$pred_no_tx.rt = predict(fit_all, newdata = dat_tx.rt, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx.rt = sum(dat_all$pred_no_tx.rt, na.rm = TRUE)
af_by_tx.rt = (N_all - N_no_tx.rt) / N_all
# af_by_tx.rt <- round(af_by_tx.rt,3)
af_by_tx.rt

## Male
N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$gender == "Male"], na.rm = TRUE)
af_by_tx.rt.male = (N_all.male - N_no_tx.rt) / N_all.male
# af_by_tx.rt.male <- round(af_by_tx.rt.male,3)
af_by_tx.rt.male

## Female
N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$gender == "Female"], na.rm = TRUE)
af_by_tx.rt.female = (N_all.female - N_no_tx.rt) / N_all.female
# af_by_tx.rt.female <- round(af_by_tx.rt.female,3)
af_by_tx.rt.female

## < 35
N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_tx.rt.lt.35 = (N_all.lt.35 - N_no_tx.rt) / N_all.lt.35
# af_by_tx.rt.lt.35 <- round(af_by_tx.rt.lt.35,3)
af_by_tx.rt.lt.35

## >= 35
N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_tx.rt.gteq.35 = (N_all.gteq.35 - N_no_tx.rt) / N_all.gteq.35
# af_by_tx.rt.gteq.35 <- round(af_by_tx.rt.gteq.35,3)
af_by_tx.rt.gteq.35

## EUR
N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$admixture == "EUR"], na.rm = TRUE)
af_by_tx.rt.EUR = (N_all.EUR - N_no_tx.rt) / N_all.EUR
# af_by_tx.rt.EUR <- round(af_by_tx.rt.EUR,3)
af_by_tx.rt.EUR

## AFR
N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$admixture == "AFR"], na.rm = TRUE)
af_by_tx.rt.AFR = (N_all.AFR - N_no_tx.rt) / N_all.AFR
# af_by_tx.rt.AFR <- round(af_by_tx.rt.AFR,3)
af_by_tx.rt.AFR
#########
## PRS ##
#########
## P/LP Zhaoming, Qin without Zhaoming and PRS
dat_prs = dat_all
# dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N";
dat_prs$Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category = "1st"  # **

dat_all$pred_no_prs = predict(fit_all, newdata = dat_prs, type = "response")
N_no_prs = sum(dat_all$pred_no_prs, na.rm = TRUE)
af_by_prs = (N_all - N_no_prs) / N_all
# af_by_prs <- round(af_by_prs,3)
af_by_prs

## Male
N_no_prs = sum(dat_all$pred_no_prs[dat_all$gender == "Male"], na.rm = TRUE)
af_by_prs.male = (N_all.male - N_no_prs) / N_all.male
# af_by_prs.male <- round(af_by_prs.male,3)
af_by_prs.male

## Female
N_no_prs = sum(dat_all$pred_no_prs[dat_all$gender == "Female"], na.rm = TRUE)
af_by_prs.female = (N_all.female - N_no_prs) / N_all.female
# af_by_prs.female <- round(af_by_prs.female,3)
af_by_prs.female

## < 35
N_no_prs = sum(dat_all$pred_no_prs[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_prs.lt.35 = (N_all.lt.35 - N_no_prs) / N_all.lt.35
# af_by_prs.lt.35 <- round(af_by_prs.lt.35,3)
af_by_prs.lt.35

## >= 35
N_no_prs = sum(dat_all$pred_no_prs[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_prs.gteq.35 = (N_all.gteq.35 - N_no_prs) / N_all.gteq.35
# af_by_prs.gteq.35 <- round(af_by_prs.gteq.35,3)
af_by_prs.gteq.35

## EUR
N_no_prs = sum(dat_all$pred_no_prs[dat_all$admixture == "EUR"], na.rm = TRUE)
af_by_prs.EUR = (N_all.EUR - N_no_prs) / N_all.EUR
# af_by_prs.EUR <- round(af_by_prs.EUR,3)
af_by_prs.EUR

## AFR
N_no_prs = sum(dat_all$pred_no_prs[dat_all$admixture == "AFR"], na.rm = TRUE)
af_by_prs.AFR = (N_all.AFR - N_no_prs) / N_all.AFR
# af_by_prs.AFR <- round(af_by_prs.AFR,3)
af_by_prs.AFR
###############
## Lifestyle ##
###############
af_by_no_favorable_lifestyle.category <- "-"

## Male
af_by_no_favorable_lifestyle.category.male <- "-"

## Female
af_by_no_favorable_lifestyle.category.female <- "-"

## < 35
af_by_no_favorable_lifestyle.category.lt.35 <- "-"

## >= 35
af_by_no_favorable_lifestyle.category.gteq.35 <- "-"

## EUR
af_by_no_favorable_lifestyle.category.EUR <- "-"

## AFR
af_by_no_favorable_lifestyle.category.AFR <- "-"

#################################################
## Treatment, Genetics and Lifestyle, combined ##
#################################################
dat_tx.prs.lifestyle = dat_all

dat_tx.prs.lifestyle$any_chemo_missing <- "No" ## **
dat_tx.prs.lifestyle$any_rt_missing <- "No" ## **


## Nullify Treatment
dat_tx.prs.lifestyle$maxchestrtdose.category =
  dat_tx.prs.lifestyle$anthra_jco_dose_5.category = "None" ## **

## Nullify Genetics
# dat_tx.plp.prs.lifestyle$Zhaoming_carriers = dat_tx.plp.prs.lifestyle$Qin_without_Zhaoming_vars_carriers = "N";
dat_tx.prs.lifestyle$Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category = "1st" ## **


dat_all$pred_no_combined = predict(fit_all, newdata = dat_tx.prs.lifestyle, type = "response")

N_no_combined = sum(dat_all$pred_no_combined, na.rm = TRUE)
af_by_combined = (N_all - N_no_combined) / N_all
# af_by_combined <- round(af_by_combined,3)
af_by_combined


## Male
N_no_combined = sum(dat_all$pred_no_combined[dat_all$gender == "Male"], na.rm = TRUE)
af_by_combined.male = (N_all.male - N_no_combined) / N_all.male
# af_by_combined.male <- round(af_by_combined.male,3)
af_by_combined.male

## Female
N_no_combined = sum(dat_all$pred_no_combined[dat_all$gender == "Female"], na.rm = TRUE)
af_by_combined.female = (N_all.female - N_no_combined) / N_all.female
# af_by_combined.female <- round(af_by_combined.female,3)
af_by_combined.female

## < 35
N_no_combined = sum(dat_all$pred_no_combined[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
af_by_combined.lt.35 = (N_all.lt.35 - N_no_combined) / N_all.lt.35
# af_by_combined.lt.35 <- round(af_by_combined.lt.35,3)
af_by_combined.lt.35

## >= 35
N_no_combined = sum(dat_all$pred_no_combined[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
af_by_combined.gteq.35 = (N_all.gteq.35 - N_no_combined) / N_all.gteq.35
# af_by_combined.gteq.35 <- round(af_by_combined.gteq.35,3)
af_by_combined.gteq.35

## EUR
N_no_combined = sum(dat_all$pred_no_combined[dat_all$admixture == "EUR"], na.rm = TRUE)
af_by_combined.EUR = (N_all.EUR - N_no_combined) / N_all.EUR
# af_by_combined.EUR <- round(af_by_combined.EUR,3)
af_by_combined.EUR

## AFR
N_no_combined = sum(dat_all$pred_no_combined[dat_all$admixture == "AFR"], na.rm = TRUE)
af_by_combined.AFR = (N_all.AFR - N_no_combined) / N_all.AFR
# af_by_combined.AFR <- round(af_by_combined.AFR,3)
af_by_combined.AFR

##
Breast.res <- data.frame(
  Variable = c("Radiation", "Chemo", "All treatments", "PRS", "Lifestyle", "Combined"),
  Overall = c(af_by_rt, af_by_tx, af_by_tx.rt, af_by_prs, af_by_no_favorable_lifestyle.category, af_by_combined),
  Female = c(af_by_rt.female, af_by_tx.female, af_by_tx.rt.female, af_by_prs.female, af_by_no_favorable_lifestyle.category.female, af_by_combined.female),
  Male = c(af_by_rt.male, af_by_tx.male, af_by_tx.rt.male, af_by_prs.male, af_by_no_favorable_lifestyle.category.male, af_by_combined.male),
  age.lt35 = c(af_by_rt.lt.35, af_by_tx.lt.35, af_by_tx.rt.lt.35, af_by_prs.lt.35, af_by_no_favorable_lifestyle.category.lt.35, af_by_combined.lt.35),
  age.gteq = c(af_by_rt.gteq.35, af_by_tx.gteq.35, af_by_tx.rt.gteq.35, af_by_prs.gteq.35, af_by_no_favorable_lifestyle.category.gteq.35, af_by_combined.gteq.35),
  EUR = c(af_by_rt.EUR, af_by_tx.EUR, af_by_tx.rt.EUR, af_by_prs.EUR, af_by_no_favorable_lifestyle.category.EUR, af_by_combined.EUR),
  AFR = c(af_by_rt.AFR, af_by_tx.AFR, af_by_tx.rt.AFR, af_by_prs.AFR, af_by_no_favorable_lifestyle.category.AFR, af_by_combined.AFR)
)

# # View(Breast.res)
# 
# 
# 
# 
######################################
#########################################
options(warn = -1)
library(dplyr)
lastrow=dat_all %>% 
  group_by(ccssid) %>% 
  slice(which.max(AGE_AT_LAST_CONTACT.cs1))
### Qi note: I think AGE_AT_LAST_CONTACT.cs1 is your age spline 1 and hence it is the originally age row. The above code take the last row in each person, so I can see how many unique ids.
all_ids=lastrow[,"ccssid"]


af_by_tx.save = NULL
af_by_tx.save = rbind(af_by_tx.save, af_by_tx)

af_by_rt.save = NULL
af_by_rt.save = rbind(af_by_rt.save, af_by_rt)

af_by_tx.rt.save = NULL
af_by_tx.rt.save = rbind(af_by_tx.rt.save, af_by_tx.rt)

af_by_prs.save = NULL
af_by_prs.save = rbind(af_by_prs.save, af_by_prs)

af_by_no_favorable_lifestyle.category.save = NULL
af_by_no_favorable_lifestyle.category.save = rbind(af_by_no_favorable_lifestyle.category.save, af_by_no_favorable_lifestyle.category)

af_by_combined.save = NULL
af_by_combined.save = rbind(af_by_combined.save, af_by_combined)

af_by_rt.female.save = NULL
af_by_rt.female.save = rbind(af_by_rt.female.save, af_by_rt.female)

af_by_tx.female.save = NULL
af_by_tx.female.save = rbind(af_by_tx.female.save, af_by_tx.female)

af_by_tx.rt.female.save = NULL
af_by_tx.rt.female.save = rbind(af_by_tx.rt.female.save, af_by_tx.rt.female)

af_by_prs.female.save = NULL
af_by_prs.female.save = rbind(af_by_prs.female.save, af_by_prs.female)

af_by_no_favorable_lifestyle.category.female.save = NULL
af_by_no_favorable_lifestyle.category.female.save = rbind(af_by_no_favorable_lifestyle.category.female.save, af_by_no_favorable_lifestyle.category.female)

af_by_combined.female.save = NULL
af_by_combined.female.save = rbind(af_by_combined.female.save, af_by_combined.female)

af_by_rt.male.save = NULL
af_by_rt.male.save = rbind(af_by_rt.male.save, af_by_rt.male)

af_by_tx.male.save = NULL
af_by_tx.male.save = rbind(af_by_tx.male.save, af_by_tx.male)

af_by_tx.rt.male.save = NULL
af_by_tx.rt.male.save = rbind(af_by_tx.rt.male.save, af_by_tx.rt.male)

af_by_prs.male.save = NULL
af_by_prs.male.save = rbind(af_by_prs.male.save, af_by_prs.male)

af_by_no_favorable_lifestyle.category.male.save = NULL
af_by_no_favorable_lifestyle.category.male.save = rbind(af_by_no_favorable_lifestyle.category.male.save, af_by_no_favorable_lifestyle.category.male)

af_by_combined.male.save = NULL
af_by_combined.male.save = rbind(af_by_combined.male.save, af_by_combined.male)

af_by_rt.lt.35.save = NULL
af_by_rt.lt.35.save = rbind(af_by_rt.lt.35.save, af_by_rt.lt.35)

af_by_tx.lt.35.save = NULL
af_by_tx.lt.35.save = rbind(af_by_tx.lt.35.save, af_by_tx.lt.35)

af_by_tx.rt.lt.35.save = NULL
af_by_tx.rt.lt.35.save = rbind(af_by_tx.rt.lt.35.save, af_by_tx.rt.lt.35)

af_by_prs.lt.35.save = NULL
af_by_prs.lt.35.save = rbind(af_by_prs.lt.35.save, af_by_prs.lt.35)

af_by_no_favorable_lifestyle.category.lt.35.save = NULL
af_by_no_favorable_lifestyle.category.lt.35.save = rbind(af_by_no_favorable_lifestyle.category.lt.35.save, af_by_no_favorable_lifestyle.category.lt.35)

af_by_combined.lt.35.save = NULL
af_by_combined.lt.35.save = rbind(af_by_combined.lt.35.save, af_by_combined.lt.35)

af_by_rt.gteq.35.save = NULL
af_by_rt.gteq.35.save = rbind(af_by_rt.gteq.35.save, af_by_rt.gteq.35)

af_by_tx.gteq.35.save = NULL
af_by_tx.gteq.35.save = rbind(af_by_tx.gteq.35.save, af_by_tx.gteq.35)

af_by_tx.rt.gteq.35.save = NULL
af_by_tx.rt.gteq.35.save = rbind(af_by_tx.rt.gteq.35.save, af_by_tx.rt.gteq.35)

af_by_prs.gteq.35.save = NULL
af_by_prs.gteq.35.save = rbind(af_by_prs.gteq.35.save, af_by_prs.gteq.35)

af_by_no_favorable_lifestyle.category.gteq.35.save = NULL
af_by_no_favorable_lifestyle.category.gteq.35.save = rbind(af_by_no_favorable_lifestyle.category.gteq.35.save, af_by_no_favorable_lifestyle.category.gteq.35)

af_by_combined.gteq.35.save = NULL
af_by_combined.gteq.35.save = rbind(af_by_combined.gteq.35.save, af_by_combined.gteq.35)

af_by_rt.EUR.save = NULL
af_by_rt.EUR.save = rbind(af_by_rt.EUR.save, af_by_rt.EUR)

af_by_tx.EUR.save = NULL
af_by_tx.EUR.save = rbind(af_by_tx.EUR.save, af_by_tx.EUR)

af_by_tx.rt.EUR.save = NULL
af_by_tx.rt.EUR.save = rbind(af_by_tx.rt.EUR.save, af_by_tx.rt.EUR)

af_by_prs.EUR.save = NULL
af_by_prs.EUR.save = rbind(af_by_prs.EUR.save, af_by_prs.EUR)

af_by_no_favorable_lifestyle.category.EUR.save = NULL
af_by_no_favorable_lifestyle.category.EUR.save = rbind(af_by_no_favorable_lifestyle.category.EUR.save, af_by_no_favorable_lifestyle.category.EUR)

af_by_combined.EUR.save = NULL
af_by_combined.EUR.save = rbind(af_by_combined.EUR.save, af_by_combined.EUR)

af_by_rt.AFR.save = NULL
af_by_rt.AFR.save = rbind(af_by_rt.AFR.save, af_by_rt.AFR)

af_by_tx.AFR.save = NULL
af_by_tx.AFR.save = rbind(af_by_tx.AFR.save, af_by_tx.AFR)

af_by_tx.rt.AFR.save = NULL
af_by_tx.rt.AFR.save = rbind(af_by_tx.rt.AFR.save, af_by_tx.rt.AFR)

af_by_prs.AFR.save = NULL
af_by_prs.AFR.save = rbind(af_by_prs.AFR.save, af_by_prs.AFR)

af_by_no_favorable_lifestyle.category.AFR.save = NULL
af_by_no_favorable_lifestyle.category.AFR.save = rbind(af_by_no_favorable_lifestyle.category.AFR.save, af_by_no_favorable_lifestyle.category.AFR)

af_by_combined.AFR.save = NULL
af_by_combined.AFR.save = rbind(af_by_combined.AFR.save, af_by_combined.AFR)

options(warn = -1)
dat_all2 <- dat_all
tryCatch({
  for(boot in 1:1000){
    #   boot=1;
    print(boot)
    set.seed(boot);
    rows_sample=sample(1:nrow(all_ids), nrow(all_ids), replace = TRUE, prob = NULL)           
    length(rows_sample)
    length(unique(rows_sample))
    sample_ids=all_ids[rows_sample,]
    ##make fake ids -- consider each sampled id a new person
    sample_ids$newid=seq(1:length(rows_sample))
    
    #        aa=table(sample_ids)
    # table(aa)
    ### with boot=1 as an example, 1651 ccssids were seleted once, there was 4 ids were selected 6 times. For example SJL1264701
    # aa[aa==6]
    
    ##### many to many merge to make new data
    library(dplyr)
    dat_all <- dat_all2 %>% 
      left_join(sample_ids, by = "ccssid")
    
    fit_all = glm(formula = event ~ 
                    Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category +
                    AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + 
                    AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
                    maxchestrtdose.category + anthra_jco_dose_5.category +
                    EAS + AFR +
                    any_chemo_missing + any_rt_missing,
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
    ## Subset by ancestry
    N_all.EUR = sum(dat_all$pred_all[dat_all$admixture == "EUR"], na.rm = TRUE) # subset by ancestry
    N_all.AFR = sum(dat_all$pred_all[dat_all$admixture == "AFR"], na.rm = TRUE) # subset by ancestry
    
    
    #############
    ## tx only ##
    #############
    ## Move relevant treatment exposures for everyone to no exposure
    dat_tx = dat_all
    
    dat_tx$any_chemo_missing <- "No" # **
    dat_tx$anthra_jco_dose_5.category = "None" ## **
    
    dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")
    
    ## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
    N_all = sum(dat_all$pred_all, na.rm = TRUE)
    N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
    af_by_tx = (N_all - N_no_tx) / N_all
    # af_by_tx <- round(af_by_tx,3)
    af_by_tx
    
    ## Male
    N_no_tx = sum(dat_all$pred_no_tx[dat_all$gender == "Male"], na.rm = TRUE)
    af_by_tx.male = (N_all.male - N_no_tx) / N_all.male
    # af_by_tx.male <- round(af_by_tx.male,3)
    af_by_tx.male
    
    ## Female
    N_no_tx = sum(dat_all$pred_no_tx[dat_all$gender == "Female"], na.rm = TRUE)
    af_by_tx.female = (N_all.female - N_no_tx) / N_all.female
    # af_by_tx.female <- round(af_by_tx.female,3)
    af_by_tx.female
    
    ## < 35
    N_no_tx = sum(dat_all$pred_no_tx[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
    af_by_tx.lt.35 = (N_all.lt.35 - N_no_tx) / N_all.lt.35
    # af_by_tx.lt.35 <- round(af_by_tx.lt.35,3)
    af_by_tx.lt.35
    
    ## >= 35
    N_no_tx = sum(dat_all$pred_no_tx[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
    af_by_tx.gteq.35 = (N_all.gteq.35 - N_no_tx) / N_all.gteq.35
    # af_by_tx.gteq.35 <- round(af_by_tx.gteq.35,3)
    af_by_tx.gteq.35
    
    ## EUR
    N_no_tx = sum(dat_all$pred_no_tx[dat_all$admixture == "EUR"], na.rm = TRUE)
    af_by_tx.EUR = (N_all.EUR - N_no_tx) / N_all.EUR
    # af_by_tx.EUR <- round(af_by_tx.EUR,3)
    af_by_tx.EUR
    
    ## AFR
    N_no_tx = sum(dat_all$pred_no_tx[dat_all$admixture == "AFR"], na.rm = TRUE)
    af_by_tx.AFR = (N_all.AFR - N_no_tx) / N_all.AFR
    # af_by_tx.AFR <- round(af_by_tx.AFR,3)
    af_by_tx.AFR
    
    #############
    ## RT only ##
    #############
    
    ## Move relevant treatment exposures for everyone to no exposure
    dat_rt = dat_all
    
    dat_rt$any_rt_missing <- "No" # **
    
    dat_rt$maxchestrtdose.category = "None" ## **
    
    dat_all$pred_no_rt = predict(fit_all, newdata = dat_rt, type = "response")
    
    ## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
    N_all = sum(dat_all$pred_all, na.rm = TRUE)
    N_no_rt = sum(dat_all$pred_no_rt, na.rm = TRUE)
    af_by_rt = (N_all - N_no_rt) / N_all
    # af_by_rt <- round(af_by_rt,3)
    af_by_rt
    
    ## Male
    N_no_rt = sum(dat_all$pred_no_rt[dat_all$gender == "Male"], na.rm = TRUE)
    af_by_rt.male = (N_all.male - N_no_rt) / N_all.male
    # af_by_rt.male <- round(af_by_rt.male,3)
    af_by_rt.male
    
    ## Female
    N_no_rt = sum(dat_all$pred_no_rt[dat_all$gender == "Female"], na.rm = TRUE)
    af_by_rt.female = (N_all.female - N_no_rt) / N_all.female
    # af_by_rt.female <- round(af_by_rt.female,3)
    af_by_rt.female
    
    ## < 35
    N_no_rt = sum(dat_all$pred_no_rt[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
    af_by_rt.lt.35 = (N_all.lt.35 - N_no_rt) / N_all.lt.35
    # af_by_rt.lt.35 <- round(af_by_rt.lt.35,3)
    af_by_rt.lt.35
    
    ## >= 35
    N_no_rt = sum(dat_all$pred_no_rt[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
    af_by_rt.gteq.35 = (N_all.gteq.35 - N_no_rt) / N_all.gteq.35
    # af_by_rt.gteq.35 <- round(af_by_rt.gteq.35,3)
    af_by_rt.gteq.35
    
    ## EUR
    N_no_rt = sum(dat_all$pred_no_rt[dat_all$admixture == "EUR"], na.rm = TRUE)
    af_by_rt.EUR = (N_all.EUR - N_no_rt) / N_all.EUR
    # af_by_rt.EUR <- round(af_by_rt.EUR,3)
    af_by_rt.EUR
    
    ## AFR
    N_no_rt = sum(dat_all$pred_no_rt[dat_all$admixture == "AFR"], na.rm = TRUE)
    af_by_rt.AFR = (N_all.AFR - N_no_rt) / N_all.AFR
    # af_by_rt.AFR <- round(af_by_rt.AFR,3)
    af_by_rt.AFR
    
    ######################
    ## Treatment and RT ##
    ######################
    
    ## Move relevant treatment exposures for everyone to no exposure
    dat_tx.rt = dat_all
    
    dat_tx.rt$any_chemo_missing <- "No" ## **
    dat_tx.rt$any_rt_missing <- "No" ## **
    
    dat_tx.rt$maxchestrtdose.category =
      dat_tx.rt$anthra_jco_dose_5.category = "None" ## **
    
    dat_all$pred_no_tx.rt = predict(fit_all, newdata = dat_tx.rt, type = "response")
    
    ## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
    N_all = sum(dat_all$pred_all, na.rm = TRUE)
    N_no_tx.rt = sum(dat_all$pred_no_tx.rt, na.rm = TRUE)
    af_by_tx.rt = (N_all - N_no_tx.rt) / N_all
    # af_by_tx.rt <- round(af_by_tx.rt,3)
    af_by_tx.rt
    
    ## Male
    N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$gender == "Male"], na.rm = TRUE)
    af_by_tx.rt.male = (N_all.male - N_no_tx.rt) / N_all.male
    # af_by_tx.rt.male <- round(af_by_tx.rt.male,3)
    af_by_tx.rt.male
    
    ## Female
    N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$gender == "Female"], na.rm = TRUE)
    af_by_tx.rt.female = (N_all.female - N_no_tx.rt) / N_all.female
    # af_by_tx.rt.female <- round(af_by_tx.rt.female,3)
    af_by_tx.rt.female
    
    ## < 35
    N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
    af_by_tx.rt.lt.35 = (N_all.lt.35 - N_no_tx.rt) / N_all.lt.35
    # af_by_tx.rt.lt.35 <- round(af_by_tx.rt.lt.35,3)
    af_by_tx.rt.lt.35
    
    ## >= 35
    N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
    af_by_tx.rt.gteq.35 = (N_all.gteq.35 - N_no_tx.rt) / N_all.gteq.35
    # af_by_tx.rt.gteq.35 <- round(af_by_tx.rt.gteq.35,3)
    af_by_tx.rt.gteq.35
    
    ## EUR
    N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$admixture == "EUR"], na.rm = TRUE)
    af_by_tx.rt.EUR = (N_all.EUR - N_no_tx.rt) / N_all.EUR
    # af_by_tx.rt.EUR <- round(af_by_tx.rt.EUR,3)
    af_by_tx.rt.EUR
    
    ## AFR
    N_no_tx.rt = sum(dat_all$pred_no_tx.rt[dat_all$admixture == "AFR"], na.rm = TRUE)
    af_by_tx.rt.AFR = (N_all.AFR - N_no_tx.rt) / N_all.AFR
    # af_by_tx.rt.AFR <- round(af_by_tx.rt.AFR,3)
    af_by_tx.rt.AFR
    #########
    ## PRS ##
    #########
    ## P/LP Zhaoming, Qin without Zhaoming and PRS
    dat_prs = dat_all
    # dat_plp.prs$Zhaoming_carriers = dat_plp.prs$Qin_without_Zhaoming_vars_carriers = "N";
    dat_prs$Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category = "1st"  # **
    
    dat_all$pred_no_prs = predict(fit_all, newdata = dat_prs, type = "response")
    N_no_prs = sum(dat_all$pred_no_prs, na.rm = TRUE)
    af_by_prs = (N_all - N_no_prs) / N_all
    # af_by_prs <- round(af_by_prs,3)
    af_by_prs
    
    ## Male
    N_no_prs = sum(dat_all$pred_no_prs[dat_all$gender == "Male"], na.rm = TRUE)
    af_by_prs.male = (N_all.male - N_no_prs) / N_all.male
    # af_by_prs.male <- round(af_by_prs.male,3)
    af_by_prs.male
    
    ## Female
    N_no_prs = sum(dat_all$pred_no_prs[dat_all$gender == "Female"], na.rm = TRUE)
    af_by_prs.female = (N_all.female - N_no_prs) / N_all.female
    # af_by_prs.female <- round(af_by_prs.female,3)
    af_by_prs.female
    
    ## < 35
    N_no_prs = sum(dat_all$pred_no_prs[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
    af_by_prs.lt.35 = (N_all.lt.35 - N_no_prs) / N_all.lt.35
    # af_by_prs.lt.35 <- round(af_by_prs.lt.35,3)
    af_by_prs.lt.35
    
    ## >= 35
    N_no_prs = sum(dat_all$pred_no_prs[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
    af_by_prs.gteq.35 = (N_all.gteq.35 - N_no_prs) / N_all.gteq.35
    # af_by_prs.gteq.35 <- round(af_by_prs.gteq.35,3)
    af_by_prs.gteq.35
    
    ## EUR
    N_no_prs = sum(dat_all$pred_no_prs[dat_all$admixture == "EUR"], na.rm = TRUE)
    af_by_prs.EUR = (N_all.EUR - N_no_prs) / N_all.EUR
    # af_by_prs.EUR <- round(af_by_prs.EUR,3)
    af_by_prs.EUR
    
    ## AFR
    N_no_prs = sum(dat_all$pred_no_prs[dat_all$admixture == "AFR"], na.rm = TRUE)
    af_by_prs.AFR = (N_all.AFR - N_no_prs) / N_all.AFR
    # af_by_prs.AFR <- round(af_by_prs.AFR,3)
    af_by_prs.AFR
    ###############
    ## Lifestyle ##
    ###############
    af_by_no_favorable_lifestyle.category <- "-"
    
    ## Male
    af_by_no_favorable_lifestyle.category.male <- "-"
    
    ## Female
    af_by_no_favorable_lifestyle.category.female <- "-"
    
    ## < 35
    af_by_no_favorable_lifestyle.category.lt.35 <- "-"
    
    ## >= 35
    af_by_no_favorable_lifestyle.category.gteq.35 <- "-"
    
    ## EUR
    af_by_no_favorable_lifestyle.category.EUR <- "-"
    
    ## AFR
    af_by_no_favorable_lifestyle.category.AFR <- "-"
    
    #################################################
    ## Treatment, Genetics and Lifestyle, combined ##
    #################################################
    dat_tx.prs.lifestyle = dat_all
    
    dat_tx.prs.lifestyle$any_chemo_missing <- "No" ## **
    dat_tx.prs.lifestyle$any_rt_missing <- "No" ## **
    
    
    ## Nullify Treatment
    dat_tx.prs.lifestyle$maxchestrtdose.category =
      dat_tx.prs.lifestyle$anthra_jco_dose_5.category = "None" ## **
    
    ## Nullify Genetics
    # dat_tx.plp.prs.lifestyle$Zhaoming_carriers = dat_tx.plp.prs.lifestyle$Qin_without_Zhaoming_vars_carriers = "N";
    dat_tx.prs.lifestyle$Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category = "1st" ## **
    
    
    dat_all$pred_no_combined = predict(fit_all, newdata = dat_tx.prs.lifestyle, type = "response")
    
    N_no_combined = sum(dat_all$pred_no_combined, na.rm = TRUE)
    af_by_combined = (N_all - N_no_combined) / N_all
    # af_by_combined <- round(af_by_combined,3)
    af_by_combined
    
    
    ## Male
    N_no_combined = sum(dat_all$pred_no_combined[dat_all$gender == "Male"], na.rm = TRUE)
    af_by_combined.male = (N_all.male - N_no_combined) / N_all.male
    # af_by_combined.male <- round(af_by_combined.male,3)
    af_by_combined.male
    
    ## Female
    N_no_combined = sum(dat_all$pred_no_combined[dat_all$gender == "Female"], na.rm = TRUE)
    af_by_combined.female = (N_all.female - N_no_combined) / N_all.female
    # af_by_combined.female <- round(af_by_combined.female,3)
    af_by_combined.female
    
    ## < 35
    N_no_combined = sum(dat_all$pred_no_combined[dat_all$AGE_AT_LAST_CONTACT.cs1 < 35], na.rm = TRUE)
    af_by_combined.lt.35 = (N_all.lt.35 - N_no_combined) / N_all.lt.35
    # af_by_combined.lt.35 <- round(af_by_combined.lt.35,3)
    af_by_combined.lt.35
    
    ## >= 35
    N_no_combined = sum(dat_all$pred_no_combined[dat_all$AGE_AT_LAST_CONTACT.cs1 >= 35], na.rm = TRUE)
    af_by_combined.gteq.35 = (N_all.gteq.35 - N_no_combined) / N_all.gteq.35
    # af_by_combined.gteq.35 <- round(af_by_combined.gteq.35,3)
    af_by_combined.gteq.35
    
    ## EUR
    N_no_combined = sum(dat_all$pred_no_combined[dat_all$admixture == "EUR"], na.rm = TRUE)
    af_by_combined.EUR = (N_all.EUR - N_no_combined) / N_all.EUR
    # af_by_combined.EUR <- round(af_by_combined.EUR,3)
    af_by_combined.EUR
    
    ## AFR
    N_no_combined = sum(dat_all$pred_no_combined[dat_all$admixture == "AFR"], na.rm = TRUE)
    af_by_combined.AFR = (N_all.AFR - N_no_combined) / N_all.AFR
    # af_by_combined.AFR <- round(af_by_combined.AFR,3)
    af_by_combined.AFR
    
    #########################################
    af_by_tx.save = rbind(af_by_tx.save, af_by_tx)
    
    af_by_rt.save = rbind(af_by_rt.save, af_by_rt)
    af_by_tx.save = rbind(af_by_tx.save, af_by_tx)
    af_by_tx.rt.save = rbind(af_by_tx.rt.save, af_by_tx.rt)
    af_by_prs.save = rbind(af_by_prs.save, af_by_prs)
    af_by_no_favorable_lifestyle.category.save = rbind(af_by_no_favorable_lifestyle.category.save, af_by_no_favorable_lifestyle.category)
    af_by_combined.save = rbind(af_by_combined.save, af_by_combined)
    
    af_by_rt.female.save = rbind(af_by_rt.female.save, af_by_rt.female)
    af_by_tx.female.save = rbind(af_by_tx.female.save, af_by_tx.female)
    af_by_tx.rt.female.save = rbind(af_by_tx.rt.female.save, af_by_tx.rt.female)
    af_by_prs.female.save = rbind(af_by_prs.female.save, af_by_prs.female)
    af_by_no_favorable_lifestyle.category.female.save = rbind(af_by_no_favorable_lifestyle.category.female.save, af_by_no_favorable_lifestyle.category.female)
    af_by_combined.female.save = rbind(af_by_combined.female.save, af_by_combined.female)
    
    af_by_rt.male.save = rbind(af_by_rt.male.save, af_by_rt.male)
    af_by_tx.male.save = rbind(af_by_tx.male.save, af_by_tx.male)
    af_by_tx.rt.male.save = rbind(af_by_tx.rt.male.save, af_by_tx.rt.male)
    af_by_prs.male.save = rbind(af_by_prs.male.save, af_by_prs.male)
    af_by_no_favorable_lifestyle.category.male.save = rbind(af_by_no_favorable_lifestyle.category.male.save, af_by_no_favorable_lifestyle.category.male)
    af_by_combined.male.save = rbind(af_by_combined.male.save, af_by_combined.male)
    
    af_by_rt.lt.35.save = rbind(af_by_rt.lt.35.save, af_by_rt.lt.35)
    af_by_tx.lt.35.save = rbind(af_by_tx.lt.35.save, af_by_tx.lt.35)
    af_by_tx.rt.lt.35.save = rbind(af_by_tx.rt.lt.35.save, af_by_tx.rt.lt.35)
    af_by_prs.lt.35.save = rbind(af_by_prs.lt.35.save, af_by_prs.lt.35)
    af_by_no_favorable_lifestyle.category.lt.35.save = rbind(af_by_no_favorable_lifestyle.category.lt.35.save, af_by_no_favorable_lifestyle.category.lt.35)
    af_by_combined.lt.35.save = rbind(af_by_combined.lt.35.save, af_by_combined.lt.35)
    
    af_by_rt.gteq.35.save = rbind(af_by_rt.gteq.35.save, af_by_rt.gteq.35)
    af_by_tx.gteq.35.save = rbind(af_by_tx.gteq.35.save, af_by_tx.gteq.35)
    af_by_tx.rt.gteq.35.save = rbind(af_by_tx.rt.gteq.35.save, af_by_tx.rt.gteq.35)
    af_by_prs.gteq.35.save = rbind(af_by_prs.gteq.35.save, af_by_prs.gteq.35)
    af_by_no_favorable_lifestyle.category.gteq.35.save = rbind(af_by_no_favorable_lifestyle.category.gteq.35.save, af_by_no_favorable_lifestyle.category.gteq.35)
    af_by_combined.gteq.35.save = rbind(af_by_combined.gteq.35.save, af_by_combined.gteq.35)
    
    af_by_rt.EUR.save = rbind(af_by_rt.EUR.save, af_by_rt.EUR)
    af_by_tx.EUR.save = rbind(af_by_tx.EUR.save, af_by_tx.EUR)
    af_by_tx.rt.EUR.save = rbind(af_by_tx.rt.EUR.save, af_by_tx.rt.EUR)
    af_by_prs.EUR.save = rbind(af_by_prs.EUR.save, af_by_prs.EUR)
    af_by_no_favorable_lifestyle.category.EUR.save = rbind(af_by_no_favorable_lifestyle.category.EUR.save, af_by_no_favorable_lifestyle.category.EUR)
    af_by_combined.EUR.save = rbind(af_by_combined.EUR.save, af_by_combined.EUR)
    
    af_by_rt.AFR.save = rbind(af_by_rt.AFR.save, af_by_rt.AFR)
    af_by_tx.AFR.save = rbind(af_by_tx.AFR.save, af_by_tx.AFR)
    af_by_tx.rt.AFR.save = rbind(af_by_tx.rt.AFR.save, af_by_tx.rt.AFR)
    af_by_prs.AFR.save = rbind(af_by_prs.AFR.save, af_by_prs.AFR)
    af_by_no_favorable_lifestyle.category.AFR.save = rbind(af_by_no_favorable_lifestyle.category.AFR.save, af_by_no_favorable_lifestyle.category.AFR)
    af_by_combined.AFR.save = rbind(af_by_combined.AFR.save, af_by_combined.AFR)
    
  }
  
  
  #############################################
  
  
  #############################################
  
  
  af_by_rt = paste0(round(af_by_rt.save[1,], 3), "_", 
                    round(quantile(af_by_rt.save[-1,], probs = c(0.025,0.975))[1], 3), "-",
                    round(quantile(af_by_rt.save[-1,], probs = c(0.025,0.975))[2],3))
  
  af_by_rt.female = paste0(round(af_by_rt.female.save[1,], 3), "_", 
                           round(quantile(af_by_rt.female.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                           round(quantile(af_by_rt.female.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  # af_by_rt.male = paste0(round(af_by_rt.male.save[1,], 3), "_", 
  #                        round(quantile(af_by_rt.male.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
  #                        round(quantile(af_by_rt.male.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  
  af_by_rt.lt.35 = paste0(round(af_by_rt.lt.35.save[1,], 3), "_", 
                          round(quantile(af_by_rt.lt.35.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                          round(quantile(af_by_rt.lt.35.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  af_by_rt.gteq.35 = paste0(round(af_by_rt.gteq.35.save[1,], 3), "_", 
                            round(quantile(af_by_rt.gteq.35.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                            round(quantile(af_by_rt.gteq.35.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  
  af_by_rt.EUR = paste0(round(af_by_rt.EUR.save[1,], 3), "_", 
                        round(quantile(af_by_rt.EUR.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                        round(quantile(af_by_rt.EUR.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  af_by_rt.AFR = paste0(round(af_by_rt.AFR.save[1,], 3), "_", 
                        round(quantile(af_by_rt.AFR.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                        round(quantile(af_by_rt.AFR.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  
  ##
  
  af_by_tx = paste0(round(af_by_tx.save[1,], 3), "_", 
                    round(quantile(af_by_tx.save[-1,], probs = c(0.025,0.975))[1], 3), "-",
                    round(quantile(af_by_tx.save[-1,], probs = c(0.025,0.975))[2],3))
  
  af_by_tx.female = paste0(round(af_by_tx.female.save[1,], 3), "_", 
                           round(quantile(af_by_tx.female.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                           round(quantile(af_by_tx.female.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  # af_by_tx.male = paste0(round(af_by_tx.male.save[1,], 3), "_", 
  #                        round(quantile(af_by_tx.male.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
  #                        round(quantile(af_by_tx.male.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  
  af_by_tx.lt.35 = paste0(round(af_by_tx.lt.35.save[1,], 3), "_", 
                          round(quantile(af_by_tx.lt.35.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                          round(quantile(af_by_tx.lt.35.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  af_by_tx.gteq.35 = paste0(round(af_by_tx.gteq.35.save[1,], 3), "_", 
                            round(quantile(af_by_tx.gteq.35.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                            round(quantile(af_by_tx.gteq.35.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  
  af_by_tx.EUR = paste0(round(af_by_tx.EUR.save[1,], 3), "_", 
                        round(quantile(af_by_tx.EUR.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                        round(quantile(af_by_tx.EUR.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  af_by_tx.AFR = paste0(round(af_by_tx.AFR.save[1,], 3), "_", 
                        round(quantile(af_by_tx.AFR.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                        round(quantile(af_by_tx.AFR.save[-1,], probs = c(0.025, 0.975))[2], 3))
  ##
  
  
  af_by_tx.rt = paste0(round(af_by_tx.rt.save[1,], 3), "_", 
                       round(quantile(af_by_tx.rt.save[-1,], probs = c(0.025,0.975))[1], 3), "-",
                       round(quantile(af_by_tx.rt.save[-1,], probs = c(0.025,0.975))[2],3))
  
  af_by_tx.rt.female = paste0(round(af_by_tx.rt.female.save[1,], 3), "_", 
                              round(quantile(af_by_tx.rt.female.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                              round(quantile(af_by_tx.rt.female.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  # af_by_tx.rt.male = paste0(round(af_by_tx.rt.male.save[1,], 3), "_", 
  #                           round(quantile(af_by_tx.rt.male.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
  #                           round(quantile(af_by_tx.rt.male.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  
  af_by_tx.rt.lt.35 = paste0(round(af_by_tx.rt.lt.35.save[1,], 3), "_", 
                             round(quantile(af_by_tx.rt.lt.35.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                             round(quantile(af_by_tx.rt.lt.35.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  af_by_tx.rt.gteq.35 = paste0(round(af_by_tx.rt.gteq.35.save[1,], 3), "_", 
                               round(quantile(af_by_tx.rt.gteq.35.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                               round(quantile(af_by_tx.rt.gteq.35.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  
  af_by_tx.rt.EUR = paste0(round(af_by_tx.rt.EUR.save[1,], 3), "_", 
                           round(quantile(af_by_tx.rt.EUR.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                           round(quantile(af_by_tx.rt.EUR.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  af_by_tx.rt.AFR = paste0(round(af_by_tx.rt.AFR.save[1,], 3), "_", 
                           round(quantile(af_by_tx.rt.AFR.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                           round(quantile(af_by_tx.rt.AFR.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  
  ##
  
  af_by_prs = paste0(round(af_by_prs.save[1,], 3), "_", 
                     round(quantile(af_by_prs.save[-1,], probs = c(0.025,0.975))[1], 3), "-",
                     round(quantile(af_by_prs.save[-1,], probs = c(0.025,0.975))[2],3))
  
  af_by_prs.female = paste0(round(af_by_prs.female.save[1,], 3), "_", 
                            round(quantile(af_by_prs.female.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                            round(quantile(af_by_prs.female.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  # af_by_prs.male = paste0(round(af_by_prs.male.save[1,], 3), "_", 
  #                         round(quantile(af_by_prs.male.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
  #                         round(quantile(af_by_prs.male.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  
  af_by_prs.lt.35 = paste0(round(af_by_prs.lt.35.save[1,], 3), "_", 
                           round(quantile(af_by_prs.lt.35.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                           round(quantile(af_by_prs.lt.35.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  af_by_prs.gteq.35 = paste0(round(af_by_prs.gteq.35.save[1,], 3), "_", 
                             round(quantile(af_by_prs.gteq.35.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                             round(quantile(af_by_prs.gteq.35.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  
  af_by_prs.EUR = paste0(round(af_by_prs.EUR.save[1,], 3), "_", 
                         round(quantile(af_by_prs.EUR.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                         round(quantile(af_by_prs.EUR.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  af_by_prs.AFR = paste0(round(af_by_prs.AFR.save[1,], 3), "_", 
                         round(quantile(af_by_prs.AFR.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                         round(quantile(af_by_prs.AFR.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  # ##
  # 
  # af_by_no_favorable_lifestyle.category = paste0(round(af_by_no_favorable_lifestyle.category.save[1,], 3), "_", 
  #                                                round(quantile(af_by_no_favorable_lifestyle.category.save[-1,], probs = c(0.025,0.975))[1], 3), "-",
  #                                                round(quantile(af_by_no_favorable_lifestyle.category.save[-1,], probs = c(0.025,0.975))[2],3))
  # 
  # af_by_no_favorable_lifestyle.category.female = paste0(round(af_by_no_favorable_lifestyle.category.female.save[1,], 3), "_", 
  #                                                       round(quantile(af_by_no_favorable_lifestyle.category.female.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
  #                                                       round(quantile(af_by_no_favorable_lifestyle.category.female.save[-1,], probs = c(0.025, 0.975))[2], 3))
  # 
  # af_by_no_favorable_lifestyle.category.male = paste0(round(af_by_no_favorable_lifestyle.category.male.save[1,], 3), "_", 
  #                                                     round(quantile(af_by_no_favorable_lifestyle.category.male.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
  #                                                     round(quantile(af_by_no_favorable_lifestyle.category.male.save[-1,], probs = c(0.025, 0.975))[2], 3))
  # 
  # 
  # af_by_no_favorable_lifestyle.category.lt.35 = paste0(round(af_by_no_favorable_lifestyle.category.lt.35.save[1,], 3), "_", 
  #                                                      round(quantile(af_by_no_favorable_lifestyle.category.lt.35.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
  #                                                      round(quantile(af_by_no_favorable_lifestyle.category.lt.35.save[-1,], probs = c(0.025, 0.975))[2], 3))
  # 
  # af_by_no_favorable_lifestyle.category.gteq.35 = paste0(round(af_by_no_favorable_lifestyle.category.gteq.35.save[1,], 3), "_", 
  #                                                        round(quantile(af_by_no_favorable_lifestyle.category.gteq.35.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
  #                                                        round(quantile(af_by_no_favorable_lifestyle.category.gteq.35.save[-1,], probs = c(0.025, 0.975))[2], 3))
  # 
  # 
  # af_by_no_favorable_lifestyle.category.EUR = paste0(round(af_by_no_favorable_lifestyle.category.EUR.save[1,], 3), "_", 
  #                                                    round(quantile(af_by_no_favorable_lifestyle.category.EUR.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
  #                                                    round(quantile(af_by_no_favorable_lifestyle.category.EUR.save[-1,], probs = c(0.025, 0.975))[2], 3))
  # 
  # af_by_no_favorable_lifestyle.category.AFR = paste0(round(af_by_no_favorable_lifestyle.category.AFR.save[1,], 3), "_", 
  #                                                    round(quantile(af_by_no_favorable_lifestyle.category.AFR.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
  #                                                    round(quantile(af_by_no_favorable_lifestyle.category.AFR.save[-1,], probs = c(0.025, 0.975))[2], 3))
  # 
  # ##
  
  
  af_by_combined = paste0(round(af_by_combined.save[1,], 3), "_", 
                          round(quantile(af_by_combined.save[-1,], probs = c(0.025,0.975))[1], 3), "-",
                          round(quantile(af_by_combined.save[-1,], probs = c(0.025,0.975))[2],3))
  
  af_by_combined.female = paste0(round(af_by_combined.female.save[1,], 3), "_", 
                                 round(quantile(af_by_combined.female.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                                 round(quantile(af_by_combined.female.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  # af_by_combined.male = paste0(round(af_by_combined.male.save[1,], 3), "_", 
  #                              round(quantile(af_by_combined.male.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
  #                              round(quantile(af_by_combined.male.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  
  af_by_combined.lt.35 = paste0(round(af_by_combined.lt.35.save[1,], 3), "_", 
                                round(quantile(af_by_combined.lt.35.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                                round(quantile(af_by_combined.lt.35.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  af_by_combined.gteq.35 = paste0(round(af_by_combined.gteq.35.save[1,], 3), "_", 
                                  round(quantile(af_by_combined.gteq.35.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                                  round(quantile(af_by_combined.gteq.35.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  
  af_by_combined.EUR = paste0(round(af_by_combined.EUR.save[1,], 3), "_", 
                              round(quantile(af_by_combined.EUR.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                              round(quantile(af_by_combined.EUR.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  af_by_combined.AFR = paste0(round(af_by_combined.AFR.save[1,], 3), "_", 
                              round(quantile(af_by_combined.AFR.save[-1,], probs = c(0.025, 0.975))[1], 3), "-", 
                              round(quantile(af_by_combined.AFR.save[-1,], probs = c(0.025, 0.975))[2], 3))
  
  
  
  
  
  SN.res <- data.frame(
    Variable = c("Radiation", "Chemo", "All treatments", "PRS", "Lifestyle", "Combined"),
    Overall = c(af_by_rt, af_by_tx, af_by_tx.rt, af_by_prs, af_by_no_favorable_lifestyle.category, af_by_combined),
    Female = c(af_by_rt.female, af_by_tx.female, af_by_tx.rt.female, af_by_prs.female, af_by_no_favorable_lifestyle.category.female, af_by_combined.female),
    Male = c(af_by_rt.male, af_by_tx.male, af_by_tx.rt.male, af_by_prs.male, af_by_no_favorable_lifestyle.category.male, af_by_combined.male),
    age.lt35 = c(af_by_rt.lt.35, af_by_tx.lt.35, af_by_tx.rt.lt.35, af_by_prs.lt.35, af_by_no_favorable_lifestyle.category.lt.35, af_by_combined.lt.35),
    age.gteq = c(af_by_rt.gteq.35, af_by_tx.gteq.35, af_by_tx.rt.gteq.35, af_by_prs.gteq.35, af_by_no_favorable_lifestyle.category.gteq.35, af_by_combined.gteq.35),
    EUR = c(af_by_rt.EUR, af_by_tx.EUR, af_by_tx.rt.EUR, af_by_prs.EUR, af_by_no_favorable_lifestyle.category.EUR, af_by_combined.EUR),
    AFR = c(af_by_rt.AFR, af_by_tx.AFR, af_by_tx.rt.AFR, af_by_prs.AFR, af_by_no_favorable_lifestyle.category.AFR, af_by_combined.AFR)
  )
  
  
  write.table(SN.res, file = "/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/ANALYSIS/confidence_interval/CCSS_BREASTcancer_without_lifestyle_res_output.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  # write.table(SN.res, file = "Z:/ResearchHome/Groups/sapkogrp//projects/Genomics/common/attr_fraction/ANALYSIS/confidence_interval/any_SN_without_res_output.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  
}, error = function(e) {
  # Error handling: log the error and continue
  cat("Error:", e$message, "\n")
})


#########################################
## Check PRS and treatment interaction ##
#########################################
dat_all = PHENO.ANY_SN
dat_all = dat_all[dat_all$gender == "Female",]
dat_all=dat_all[dat_all$evt1==1,]

fit_all = glm(formula = event ~
                Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 +
                AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
                maxchestrtdose.category + anthra_jco_dose_5.category +
                EAS + AFR +
                any_chemo_missing + any_rt_missing +
                maxchestrtdose.category * Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category +
                anthra_jco_dose_5.category * Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

(output <- summary(fit_all)$coefficients)
as.data.frame(apply(output, 2, formatC, format="f", digits=4))
# options(scipen=999)
# estimate <- format(round(output[,1],3), nsmall = 3)
# std.error <- format(round(output[,2],3), nsmall = 3)
# P.val <- formatC(output[,4], format="G", digits=3)
P.val <- output[,4]
P.val[P.val < 0.001] <- "<0.001"
# P.val[!grepl("<", P.val)] <- format(round(as.numeric(P.val[!grepl("<", P.val)]), 3), nsmall = 3)
sn.model <- (setNames(cbind.data.frame(estimate, std.error, P.val
), c("Estimate", "Std.error", "P")))
sn.model <- sn.model[!grepl("AGE_AT_LAST_CONTACT", row.names(sn.model)),]
sn.model




##############################################
## Check ancestry and treatment interaction ##
##############################################
PHENO.ANY_SN$ancestry <- "Other"
PHENO.ANY_SN$ancestry [PHENO.ANY_SN$EUR > 0.8] <- "EUR"
PHENO.ANY_SN$ancestry [PHENO.ANY_SN$AFR > 0.6] <- "AFR"
PHENO.ANY_SN$ancestry <- factor(PHENO.ANY_SN$ancestry, levels = c("EUR", "AFR", "Other"))

dat_all = PHENO.ANY_SN
dat_all = dat_all[dat_all$gender == "Female",]
dat_all=dat_all[dat_all$evt1==1,]

fit_all = glm(formula = event ~
                Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 +
                AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
                maxchestrtdose.category + anthra_jco_dose_5.category +
                ancestry +
                any_chemo_missing + any_rt_missing +
                maxchestrtdose.category * ancestry +
                anthra_jco_dose_5.category * ancestry,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

(output <- summary(fit_all)$coefficients)
as.data.frame(apply(output, 2, formatC, format="f", digits=4))
# options(scipen=999)
# estimate <- format(round(output[,1],3), nsmall = 3)
# std.error <- format(round(output[,2],3), nsmall = 3)
# P.val <- formatC(output[,4], format="G", digits=3)
P.val <- output[,4]
P.val[P.val < 0.001] <- "<0.001"
# P.val[!grepl("<", P.val)] <- format(round(as.numeric(P.val[!grepl("<", P.val)]), 3), nsmall = 3)
sn.model <- (setNames(cbind.data.frame(estimate, std.error, P.val
), c("Estimate", "Std.error", "P")))
sn.model <- sn.model[!grepl("AGE_AT_LAST_CONTACT", row.names(sn.model)),]
sn.model