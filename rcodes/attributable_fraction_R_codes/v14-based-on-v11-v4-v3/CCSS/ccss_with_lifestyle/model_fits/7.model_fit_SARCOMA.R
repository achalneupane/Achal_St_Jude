# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.SARCOMA.V14.Rdata")

# Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
cc
filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
filtered_cc

PHENO.ANY_SN$aa_class_dose_5.category <- as.character(PHENO.ANY_SN$aa_class_dose_5.category)
PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "1st"] <- "1st-2nd"
PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "2nd"] <- "1st-2nd"
PHENO.ANY_SN$aa_class_dose_5.category <- factor(PHENO.ANY_SN$aa_class_dose_5.category, levels = c("None", "1st-2nd", "3rd"))

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

dat_tx$aa_class_dose_5.category = "None"

dat_all$pred_no_tx = predict(fit_all, newdata = dat_tx, type = "response")

## Attributable fraction calculation. First get the "predicted" number of SNs based on the model including all variables
N_all = sum(dat_all$pred_all, na.rm = TRUE)
N_no_tx = sum(dat_all$pred_no_tx, na.rm = TRUE)
af_by_tx = (N_all - N_no_tx) / N_all
round(af_by_tx,3)
# 0.369
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
# 0.067
###############
## Lifestyle ##
###############
dat_lifestyle = dat_all

dat_lifestyle$Current_smoker_yn = "No"
dat_lifestyle$PhysicalActivity_yn = "Yes"
dat_lifestyle$RiskyHeavyDrink_yn = "No"
# dat_lifestyle$HEALTHY_Diet_yn = "Yes"
dat_lifestyle$Obese_yn = "No"


dat_all$pred_no_favorable_lifestyle.category = predict(fit_all, newdata = dat_lifestyle, type = "response")
N_no_favorable_lifestyle.category = sum(dat_all$pred_no_favorable_lifestyle.category, na.rm = TRUE)
af_by_N_no_favorable_lifestyle.category = (N_all - N_no_favorable_lifestyle.category) / N_all
round(af_by_N_no_favorable_lifestyle.category,3)
# 0.774
#################################################
## Treatment, Genetics and Lifestyle, combined ##
#################################################

dat_tx.plp.prs.lifestyle = dat_all

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
# 0.877
SARCOMA.res <- c(round(af_by_tx,3), round(af_by_plp.prs,3),round(af_by_N_no_favorable_lifestyle.category,3), round(af_by_N_no_favorable_tx.plp.prs.lifestyle.category,3))
SARCOMA.res
# 0.204 -0.006  0.147  0.324
# 0.272 0.064 0.109 0.402

all.res <- cbind.data.frame(SN=SN.res, SMN=SMN.res, NMSC=NMSC.res, BREAST=BREAST.res, THYROID=THYROID.res, MENINGIOMA=MENINGIOMA.res, SARCOMA=SARCOMA.res)
# View(all.res)
View(t(all.res))


ChangeNames <- function(x) {
  row.names(x)[grepl("PRS.tertile.category2nd", row.names(x))] <- "PRS.tertile.category2nd"
  row.names(x)[grepl("PRS.tertile.category3rd", row.names(x))] <- "PRS.tertile.category3rd"
  return(x)
}

df_list <- list(sn.model=sn.model, smn.model=smn.model, nmsc.model=nmsc.model, breast.model=breast.model
                ,thyroid.model=thyroid.model, meningioma.model = meningioma.model, sarcoma.model = sarcoma.model)
df_list <- lapply(df_list, ChangeNames)
list2env(df_list ,.GlobalEnv)
merge.all <- function(x, ..., by = "row.names") {
  L <- list(...)
  for (i in seq_along(L)) {
    x <- merge(x, L[[i]], by = by, all = TRUE)
    rownames(x) <- x$Row.names
    x$Row.names <- NULL
  }
  return(x)
}

df <- merge.all(sn.model, smn.model,nmsc.model, breast.model, thyroid.model, meningioma.model, sarcoma.model)
df <- df[!grepl("AGE_AT_LAST_CONTACT", row.names(df)),]
dim(df)

target <- c("(Intercept)", "genderFemale", "AGE_AT_DIAGNOSIS5-9", "AGE_AT_DIAGNOSIS10-14", "AGE_AT_DIAGNOSIS>=15", "AFR", "EAS",
            "PRS.tertile.category2nd","PRS.tertile.category3rd",
            "epitxn_dose_5.category1st","epitxn_dose_5.category2nd","epitxn_dose_5.category3rd","epitxn_dose_5.categoryUnknown",
            "aa_class_dose_5.category1st","aa_class_dose_5.category2nd","aa_class_dose_5.category3rd","aa_class_dose_5.categoryUnknown",
            "anthra_jco_dose_5.category1st","anthra_jco_dose_5.category2nd","anthra_jco_dose_5.category3rd","anthra_jco_dose_5.categoryUnknown",
            "maxabdrtdose.category>0-<30", "maxabdrtdose.category>=30","maxabdrtdose.categoryUnknown",
            "maxchestrtdose.category>0-<20","maxchestrtdose.category>=20", "maxchestrtdose.categoryUnknown",
            "maxneckrtdose.category>0-<11", "maxneckrtdose.category>=11-<20","maxneckrtdose.category>=20-<30","maxneckrtdose.category>=30", "maxneckrtdose.categoryUnknown",
            "maxpelvisrtdose.category>0-<20", "maxpelvisrtdose.category>=20","maxpelvisrtdose.categoryUnknown",
            "maxsegrtdose.category>0-<18", "maxsegrtdose.category>=18-<30", "maxsegrtdose.category>=30","maxsegrtdose.categoryUnknown",
            "Current_smoker_ynYes", "Current_smoker_ynUnknown", 
            "RiskyHeavyDrink_ynYes", "RiskyHeavyDrink_ynUnknown",
            "PhysicalActivity_ynNo", "PhysicalActivity_ynUnknown",
            "Obese_ynYes", "Obese_ynUnknown")

df <- df[match(target, rownames(df)),]