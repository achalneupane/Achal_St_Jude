
## 1.------------------------------------- SJLIFE
rm(list=ls())
######################################
## Attributable fraction of Any SNs ##
######################################
# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/6.sjlife_without_lifestyle.Any_SNs.V18.Rdata")

dat_all = PHENO.ANY_SN
dat_all=dat_all[dat_all$evt1==1,]


# Fit the Poisson regression model
fit_all <- glm(formula = event ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                 AGE_AT_DIAGNOSIS + gender + 
                 maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category +
                 EAS + AFR,
               family = "poisson", offset = log(dat_all$PY), data = dat_all)

# Get the summary of the Poisson regression model
summary_fit_all <- summary(fit_all)

results <- data.frame(
  Predictor = rownames(coef(summary_fit_all)),
  "RR (95% CI)" = sprintf("%.2f (%.2f-%.2f)",
                          exp(coef(summary_fit_all)[, "Estimate"]),
                          exp(coef(summary_fit_all)[, "Estimate"] - 1.96 * coef(summary_fit_all)[, "Std. Error"]),
                          exp(coef(summary_fit_all)[, "Estimate"] + 1.96 * coef(summary_fit_all)[, "Std. Error"])
  ),
  "P-Value" = format(coef(summary_fit_all)[, "Pr(>|z|)"], nsmall = 4)
)



ANY_SN.vars <- results[grepl("diagnosis|gender|maxsegrt|maxabdrtdose|pelvis|chest|neck|aa_class_dose_5|anthra|epitxn|prs", results$Predictor, ignore.case = T),]
ANY_SN.vars
###################################
## Attributable fraction of SMNs ##
###################################
# load SMN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/6.sjlife_without_lifestyle.SMNs.V18.Rdata")

dat_all = PHENO.ANY_SN
dat_all=dat_all[dat_all$evt1==1,]
fit_all = glm(formula = event ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender +
                maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category +
                EAS + AFR,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)


# Get the summary of the Poisson regression model
summary_fit_all <- summary(fit_all)

results <- data.frame(
  Predictor = rownames(coef(summary_fit_all)),
  "RR (95% CI)" = sprintf("%.2f (%.2f-%.2f)",
                          exp(coef(summary_fit_all)[, "Estimate"]),
                          exp(coef(summary_fit_all)[, "Estimate"] - 1.96 * coef(summary_fit_all)[, "Std. Error"]),
                          exp(coef(summary_fit_all)[, "Estimate"] + 1.96 * coef(summary_fit_all)[, "Std. Error"])
  ),
  "P-Value" = format(coef(summary_fit_all)[, "Pr(>|z|)"], nsmall = 4)
)



SMN.vars <- results[grepl("diagnosis|gender|maxsegrt|maxabdrtdose|pelvis|chest|neck|aa_class_dose_5|anthra|epitxn|prs", results$Predictor, ignore.case = T),]
SMN.vars
########################################
## Attributable fraction of Any NMSCs ##
########################################

# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/6.sjlife_without_lifestyle.NMSCs.V18.Rdata")

dat_all = PHENO.ANY_SN

# dat_all=dat_all[dat_all$evt1==1,]
fit_all = glm(formula = event ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category +
                EAS + AFR,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

# Get the summary of the Poisson regression model
summary_fit_all <- summary(fit_all)

results <- data.frame(
  Predictor = rownames(coef(summary_fit_all)),
  "RR (95% CI)" = sprintf("%.2f (%.2f-%.2f)",
                          exp(coef(summary_fit_all)[, "Estimate"]),
                          exp(coef(summary_fit_all)[, "Estimate"] - 1.96 * coef(summary_fit_all)[, "Std. Error"]),
                          exp(coef(summary_fit_all)[, "Estimate"] + 1.96 * coef(summary_fit_all)[, "Std. Error"])
  ),
  "P-Value" = format(coef(summary_fit_all)[, "Pr(>|z|)"], nsmall = 4)
)



NMSC.vars <- results[grepl("diagnosis|gender|maxsegrt|maxabdrtdose|pelvis|chest|neck|aa_class_dose_5|anthra|epitxn|prs", results$Predictor, ignore.case = T),]
NMSC.vars
################################################
## Attributable fraction of Any Breast cancer ##
################################################
# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/6.sjlife_without_lifestyle.BREASTcancer.V18.Rdata")


PHENO.ANY_SN <- PHENO.ANY_SN[PHENO.ANY_SN$gender == "Female",]
dat_all = PHENO.ANY_SN
dat_all=dat_all[dat_all$evt1==1,]


# dat_all <- dat_all %>%
#   mutate(id = as.numeric(factor(sjlid)))


fit_all = glm(formula = event ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + 
                AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
                maxchestrtdose.category + anthra_jco_dose_5.category +
                EAS + AFR,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

# Get the summary of the Poisson regression model
summary_fit_all <- summary(fit_all)

results <- data.frame(
  Predictor = rownames(coef(summary_fit_all)),
  "RR (95% CI)" = sprintf("%.2f (%.2f-%.2f)",
                          exp(coef(summary_fit_all)[, "Estimate"]),
                          exp(coef(summary_fit_all)[, "Estimate"] - 1.96 * coef(summary_fit_all)[, "Std. Error"]),
                          exp(coef(summary_fit_all)[, "Estimate"] + 1.96 * coef(summary_fit_all)[, "Std. Error"])
  ),
  "P-Value" = format(coef(summary_fit_all)[, "Pr(>|z|)"], nsmall = 4)
)



BREASTCANCER.vars <- results[grepl("diagnosis|gender|maxsegrt|maxabdrtdose|pelvis|chest|neck|aa_class_dose_5|anthra|epitxn|prs", results$Predictor, ignore.case = T),]
BREASTCANCER.vars
##########################################
## Attributable fraction of Any THYROID ##
##########################################

# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/6.sjlife_without_lifestyle.THYROIDcancer.V18.Rdata")


dat_all = PHENO.ANY_SN
dat_all=dat_all[dat_all$evt1==1,]

fit_all = glm(formula = event ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + 
                maxneckrtdose.category + epitxn_dose_5.category,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

# Get the summary of the Poisson regression model
summary_fit_all <- summary(fit_all)

results <- data.frame(
  Predictor = rownames(coef(summary_fit_all)),
  "RR (95% CI)" = sprintf("%.2f (%.2f-%.2f)",
                          exp(coef(summary_fit_all)[, "Estimate"]),
                          exp(coef(summary_fit_all)[, "Estimate"] - 1.96 * coef(summary_fit_all)[, "Std. Error"]),
                          exp(coef(summary_fit_all)[, "Estimate"] + 1.96 * coef(summary_fit_all)[, "Std. Error"])
  ),
  "P-Value" = format(coef(summary_fit_all)[, "Pr(>|z|)"], nsmall = 4)
)



THYROID.vars <- results[grepl("diagnosis|gender|maxsegrt|maxabdrtdose|pelvis|chest|neck|aa_class_dose_5|anthra|epitxn|prs", results$Predictor, ignore.case = T),]
THYROID.vars

#############################################
## Attributable fraction of Any Meningioma ##
#############################################
# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/6.sjlife_without_lifestyle.MENINGIOMA.V18.Rdata")

dat_all = PHENO.ANY_SN
dat_all=dat_all[dat_all$evt1==1,]

fit_all = glm(formula = event ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + 
                maxsegrtdose.category + epitxn_dose_5.category,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

# Get the summary of the Poisson regression model
summary_fit_all <- summary(fit_all)

results <- data.frame(
  Predictor = rownames(coef(summary_fit_all)),
  "RR (95% CI)" = sprintf("%.2f (%.2f-%.2f)",
                          exp(coef(summary_fit_all)[, "Estimate"]),
                          exp(coef(summary_fit_all)[, "Estimate"] - 1.96 * coef(summary_fit_all)[, "Std. Error"]),
                          exp(coef(summary_fit_all)[, "Estimate"] + 1.96 * coef(summary_fit_all)[, "Std. Error"])
  ),
  "P-Value" = format(coef(summary_fit_all)[, "Pr(>|z|)"], nsmall = 4)
)



MENINGIOMA.vars <- results[grepl("diagnosis|gender|maxsegrt|maxabdrtdose|pelvis|chest|neck|aa_class_dose_5|anthra|epitxn|prs", results$Predictor, ignore.case = T),]
MENINGIOMA.vars

##########################################
## Attributable fraction of Any SARCOMA ##
##########################################

# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/6.sjlife_without_lifestyle.SARCOMA.V17.Rdata")

# # Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
# cc
# filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
# filtered_cc
# 
# 
# 
# PHENO.ANY_SN$aa_class_dose_5.category <- as.character(PHENO.ANY_SN$aa_class_dose_5.category)
# PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "None"] <- "1st"
# PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "2nd"] <- "2nd-3rd"
# PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "3rd"] <- "2nd-3rd"
# PHENO.ANY_SN$aa_class_dose_5.category <- factor(PHENO.ANY_SN$aa_class_dose_5.category, levels = c("1st", "2nd-3rd"))


dat_all = PHENO.ANY_SN
dat_all=dat_all[dat_all$evt1==1,]

fit_all = glm(formula = event ~ AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                gender + 
                aa_class_dose_5.category +
                EAS + AFR,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

summary_fit_all <- summary(fit_all)

results <- data.frame(
  Predictor = rownames(coef(summary_fit_all)),
  "RR (95% CI)" = sprintf("%.2f (%.2f-%.2f)",
                          exp(coef(summary_fit_all)[, "Estimate"]),
                          exp(coef(summary_fit_all)[, "Estimate"] - 1.96 * coef(summary_fit_all)[, "Std. Error"]),
                          exp(coef(summary_fit_all)[, "Estimate"] + 1.96 * coef(summary_fit_all)[, "Std. Error"])
  ),
  "P-Value" = format(coef(summary_fit_all)[, "Pr(>|z|)"], nsmall = 4)
)



SARCOMA.vars <- results[grepl("diagnosis|gender|maxsegrt|maxabdrtdose|pelvis|chest|neck|aa_class_dose|anthra|epitxn|prs", results$Predictor, ignore.case = T),]
SARCOMA.vars

## Combine all


dataframes <- list(ANY_SN=ANY_SN.vars, Meningioma=MENINGIOMA.vars, Sarcoma=SARCOMA.vars, Breast=BREASTCANCER.vars, Thyroid=THYROID.vars, NMSC=NMSC.vars, SMN=SMN.vars)
for (i in 1:length(dataframes)) {
  print(i)
  df <- dataframes[[i]]
  rownames(df) <- NULL
  # Extract the part of the column names after "PRS.tertile.category"
  df$Predictor <- sub(".*PRS.tertile.category", "PRS_tertile_", df$Predictor, ignore.case = T)
  colnames(df)[-1] <- paste0(c("RR (95% CI)", "P"), "_", names(dataframes)[i])
  dataframes[[i]] <- df
}
# dataframes[[2]]

empty_dataframe <- Reduce(function(x, y) merge(x, y, by = "Predictor", all = TRUE), dataframes)
       


# Print the resulting dataframe
print(empty_dataframe)

# Find the columns that start with "P_"
p_columns <- empty_dataframe[, grepl("^P_", colnames(empty_dataframe))]

# Convert the selected columns to numeric
for (col in colnames(p_columns)) {
  empty_dataframe[, col] <- as.numeric(empty_dataframe[, col])
  empty_dataframe[, col][empty_dataframe[, col] < 0.0001] <- "< .0001"
}


## Sort
empty_dataframe.ordered <- empty_dataframe[match(c(NA, "AGE_AT_DIAGNOSIS5-9", "AGE_AT_DIAGNOSIS10-14", "AGE_AT_DIAGNOSIS>=15", "genderFemale", 
                            NA, "maxsegrtdose.category>0-<18", "maxsegrtdose.category>=18-<30", "maxsegrtdose.category>=30", NA,
                            NA, "maxabdrtdose.category>0-<30", "maxabdrtdose.category>=30", NA,
                            NA, "maxpelvisrtdose.category>0-<20", "maxpelvisrtdose.category>=20", NA,
                            NA, "maxchestrtdose.category>0-<20", "maxchestrtdose.category>=20", NA, 
                            NA, "maxneckrtdose.category>0-<11", "maxneckrtdose.category>=11-<20", "maxneckrtdose.category>=20-<30", "maxneckrtdose.category>=30", NA,
                            NA, "aa_class_dose_5.category1st", "aa_class_dose_5.category2nd", "aa_class_dose_5.category3rd",
                            NA, "anthra_jco_dose_5.category1st",  "anthra_jco_dose_5.category2nd", "anthra_jco_dose_5.category3rd",
                            NA, "epitxn_dose_5.category1st", "epitxn_dose_5.category2nd", "epitxn_dose_5.category3rd",
                            NA, "PRS_tertile_2nd", "PRS_tertile_3rd"), empty_dataframe$Predictor),]


# Fill missing values with NA
empty_dataframe.ordered[is.na(empty_dataframe.ordered)] <- ""                  









## 2.------------------------------------- CCSS
rm(list=ls())
######################################
## Attributable fraction of Any SNs ##
######################################
# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.Any_SNs_without_lifestyle.V17b.Rdata")

# Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
# cc
# filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
# filtered_cc
# 
# PHENO.ANY_SN$maxsegrtdose.category <- as.character(PHENO.ANY_SN$maxsegrtdose.category)
# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">0-<18"] <- ">0-<30"
# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">=18-<30"] <- ">0-<30"
# PHENO.ANY_SN$maxsegrtdose.category <- factor(PHENO.ANY_SN$maxsegrtdose.category, levels = c("None", ">0-<30", ">=30"))
# 
# table(PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$event == 1]) ## **


dat_all = PHENO.ANY_SN
dat_all=dat_all[dat_all$evt1==1,]

# Fit the Poisson regression model
fit_all <- glm(formula = event ~ Pleiotropy_PRSWEB_PRS.tertile.category +
                 AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                 AGE_AT_DIAGNOSIS + gender + 
                 maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category +
                 EAS + AFR + 
                 any_chemo_missing + any_rt_missing,
               family = "poisson", offset = log(dat_all$PY), data = dat_all)

# Get the summary of the Poisson regression model
summary_fit_all <- summary(fit_all)

results <- data.frame(
  Predictor = rownames(coef(summary_fit_all)),
  "RR (95% CI)" = sprintf("%.2f (%.2f-%.2f)",
                          exp(coef(summary_fit_all)[, "Estimate"]),
                          exp(coef(summary_fit_all)[, "Estimate"] - 1.96 * coef(summary_fit_all)[, "Std. Error"]),
                          exp(coef(summary_fit_all)[, "Estimate"] + 1.96 * coef(summary_fit_all)[, "Std. Error"])
  ),
  "P-Value" = format(coef(summary_fit_all)[, "Pr(>|z|)"], nsmall = 4)
)



ANY_SN.vars <- results[grepl("diagnosis|gender|maxsegrt|maxabdrtdose|pelvis|chest|neck|aa_class_dose_5|anthra|epitxn|prs", results$Predictor, ignore.case = T),]

###################################
## Attributable fraction of SMNs ##
###################################
# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.SMNs_without_lifestyle.V17b.Rdata")


# # Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
# cc
# filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
# filtered_cc
# 
# PHENO.ANY_SN$maxsegrtdose.category <- as.character(PHENO.ANY_SN$maxsegrtdose.category)
# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">0-<18"] <- ">0-<30"
# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">=18-<30"] <- ">0-<30"
# PHENO.ANY_SN$maxsegrtdose.category <- factor(PHENO.ANY_SN$maxsegrtdose.category, levels = c("None", ">0-<30", ">=30"))


dat_all = PHENO.ANY_SN
dat_all=dat_all[dat_all$evt1==1,]
fit_all = glm(formula = event ~ Pleiotropy_PRSWEB_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender +
                maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category +
                EAS + AFR + 
                any_chemo_missing + any_rt_missing,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)


# Get the summary of the Poisson regression model
summary_fit_all <- summary(fit_all)

results <- data.frame(
  Predictor = rownames(coef(summary_fit_all)),
  "RR (95% CI)" = sprintf("%.2f (%.2f-%.2f)",
                          exp(coef(summary_fit_all)[, "Estimate"]),
                          exp(coef(summary_fit_all)[, "Estimate"] - 1.96 * coef(summary_fit_all)[, "Std. Error"]),
                          exp(coef(summary_fit_all)[, "Estimate"] + 1.96 * coef(summary_fit_all)[, "Std. Error"])
  ),
  "P-Value" = format(coef(summary_fit_all)[, "Pr(>|z|)"], nsmall = 4)
)



SMN.vars <- results[grepl("diagnosis|gender|maxsegrt|maxabdrtdose|pelvis|chest|neck|aa_class_dose_5|anthra|epitxn|prs", results$Predictor, ignore.case = T),]

########################################
## Attributable fraction of Any NMSCs ##
########################################

# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.NMSCs_without_lifestyle.V17b.Rdata")

# # Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
# cc
# filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
# filtered_cc
# 
# PHENO.ANY_SN$maxsegrtdose.category <- as.character(PHENO.ANY_SN$maxsegrtdose.category)
# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">0-<18"] <- ">0-<30"
# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">=18-<30"] <- ">0-<30"
# PHENO.ANY_SN$maxsegrtdose.category <- factor(PHENO.ANY_SN$maxsegrtdose.category, levels = c("None", ">0-<30", ">=30"))
# 
# PHENO.ANY_SN$maxpelvisrtdose.category <- as.character(PHENO.ANY_SN$maxpelvisrtdose.category)
# PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$maxpelvisrtdose.category == ">0-<20"] <- "Any"
# PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$maxpelvisrtdose.category == ">=20"] <- "Any"
# PHENO.ANY_SN$maxpelvisrtdose.category <- factor(PHENO.ANY_SN$maxpelvisrtdose.category, levels = c("None", "Any"))
# 
# # table(PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$NMSCs == 0])


dat_all = PHENO.ANY_SN

dat_all=dat_all[dat_all$evt1==1,]
fit_all = glm(formula = event ~ BASALcell_PRS.tertile.category + 
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category +
                EAS + AFR +
                any_rt_missing,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

# Get the summary of the Poisson regression model
summary_fit_all <- summary(fit_all)

results <- data.frame(
  Predictor = rownames(coef(summary_fit_all)),
  "RR (95% CI)" = sprintf("%.2f (%.2f-%.2f)",
                          exp(coef(summary_fit_all)[, "Estimate"]),
                          exp(coef(summary_fit_all)[, "Estimate"] - 1.96 * coef(summary_fit_all)[, "Std. Error"]),
                          exp(coef(summary_fit_all)[, "Estimate"] + 1.96 * coef(summary_fit_all)[, "Std. Error"])
  ),
  "P-Value" = format(coef(summary_fit_all)[, "Pr(>|z|)"], nsmall = 4)
)



NMSC.vars <- results[grepl("diagnosis|gender|maxsegrt|maxabdrtdose|pelvis|chest|neck|aa_class_dose_5|anthra|epitxn|prs", results$Predictor, ignore.case = T),]

################################################
## Attributable fraction of Any Breast cancer ##
################################################
# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.BREASTcancer_without_lifestyle.V17b.Rdata")

# # Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
# cc
# filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
# filtered_cc
# 
# PHENO.ANY_SN$maxchestrtdose.category <- as.character(PHENO.ANY_SN$maxchestrtdose.category)
# PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$maxchestrtdose.category == ">0-<20"] <- "Any"
# PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$maxchestrtdose.category == ">=20"] <- "Any"
# PHENO.ANY_SN$maxchestrtdose.category <- factor(PHENO.ANY_SN$maxchestrtdose.category, levels = c("None", "Any"))
# 
# # table(PHENO.ANY_SN$maxchestrtdose.category[PHENO.ANY_SN$BREASTcancer == 1])
# 
# PHENO.ANY_SN$anthra_jco_dose_5.category <- as.character(PHENO.ANY_SN$anthra_jco_dose_5.category)
# PHENO.ANY_SN$anthra_jco_dose_5.category[PHENO.ANY_SN$anthra_jco_dose_5.category == "1st"] <- "1st-2nd"
# PHENO.ANY_SN$anthra_jco_dose_5.category[PHENO.ANY_SN$anthra_jco_dose_5.category == "2nd"] <- "1st-2nd"
# PHENO.ANY_SN$anthra_jco_dose_5.category <- factor(PHENO.ANY_SN$anthra_jco_dose_5.category, levels = c("None", "1st-2nd", "3rd"))
# 
# # table(PHENO.ANY_SN$anthra_jco_dose_5.category[PHENO.ANY_SN$BREASTcancer == 1])

dat_all = PHENO.ANY_SN
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

# Get the summary of the Poisson regression model
summary_fit_all <- summary(fit_all)

results <- data.frame(
  Predictor = rownames(coef(summary_fit_all)),
  "RR (95% CI)" = sprintf("%.2f (%.2f-%.2f)",
                          exp(coef(summary_fit_all)[, "Estimate"]),
                          exp(coef(summary_fit_all)[, "Estimate"] - 1.96 * coef(summary_fit_all)[, "Std. Error"]),
                          exp(coef(summary_fit_all)[, "Estimate"] + 1.96 * coef(summary_fit_all)[, "Std. Error"])
  ),
  "P-Value" = format(coef(summary_fit_all)[, "Pr(>|z|)"], nsmall = 4)
)



BREASTCANCER.vars <- results[grepl("diagnosis|gender|maxsegrt|maxabdrtdose|pelvis|chest|neck|aa_class_dose_5|anthra|epitxn|prs", results$Predictor, ignore.case = T),]

##########################################
## Attributable fraction of Any THYROID ##
##########################################
# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.THYROIDcancer_without_lifestyle.V17b.Rdata")

# # Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
# cc
# filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
# filtered_cc
# 
# PHENO.ANY_SN$maxneckrtdose.category <- as.character(PHENO.ANY_SN$maxneckrtdose.category)
# PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == ">0-<11"] <- ">0-<30"
# PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == ">=11-<20"] <- ">0-<30"
# PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == ">=20-<30"] <- ">0-<30"
# PHENO.ANY_SN$maxneckrtdose.category <- factor(PHENO.ANY_SN$maxneckrtdose.category, levels = c("None", ">0-<30", ">=30"))
# 
# table(PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$THYROIDcancer == 1])



dat_all = PHENO.ANY_SN
dat_all=dat_all[dat_all$evt1==1,]

fit_all = glm(formula = event ~ Thyroid_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + 
                maxneckrtdose.category + epitxn_dose_5.category + 
                EAS + AFR +
                any_chemo_missing + any_rt_missing,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

# Get the summary of the Poisson regression model
summary_fit_all <- summary(fit_all)

results <- data.frame(
  Predictor = rownames(coef(summary_fit_all)),
  "RR (95% CI)" = sprintf("%.2f (%.2f-%.2f)",
                          exp(coef(summary_fit_all)[, "Estimate"]),
                          exp(coef(summary_fit_all)[, "Estimate"] - 1.96 * coef(summary_fit_all)[, "Std. Error"]),
                          exp(coef(summary_fit_all)[, "Estimate"] + 1.96 * coef(summary_fit_all)[, "Std. Error"])
  ),
  "P-Value" = format(coef(summary_fit_all)[, "Pr(>|z|)"], nsmall = 4)
)



THYROID.vars <- results[grepl("diagnosis|gender|maxsegrt|maxabdrtdose|pelvis|chest|neck|aa_class_dose_5|anthra|epitxn|prs", results$Predictor, ignore.case = T),]


#############################################
## Attributable fraction of Any Meningioma ##
#############################################
# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.MENINGIOMA_without_lifestyle.V17b.Rdata")

# # Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
# cc
# filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
# filtered_cc
# 
# PHENO.ANY_SN$maxsegrtdose.category <- as.character(PHENO.ANY_SN$maxsegrtdose.category)
# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == "None"] <- "<30"
# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">0-<18"] <- "<30"
# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">=18-<30"] <- "<30"
# PHENO.ANY_SN$maxsegrtdose.category <- factor(PHENO.ANY_SN$maxsegrtdose.category, levels = c("<30", ">=30"))
# 
# table(PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$MENINGIOMA == 1])
# 
# PHENO.ANY_SN$epitxn_dose_5.category <- as.character(PHENO.ANY_SN$epitxn_dose_5.category)
# PHENO.ANY_SN$epitxn_dose_5.category[PHENO.ANY_SN$epitxn_dose_5.category == "2nd"] <- "2nd-3rd"
# PHENO.ANY_SN$epitxn_dose_5.category[PHENO.ANY_SN$epitxn_dose_5.category == "3rd"] <- "2nd-3rd"
# PHENO.ANY_SN$epitxn_dose_5.category <- factor(PHENO.ANY_SN$epitxn_dose_5.category, levels = c("None", "1st", "2nd-3rd"))


dat_all = PHENO.ANY_SN
dat_all=dat_all[dat_all$evt1==1,]

fit_all = glm(formula = event ~ Meningioma_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + 
                maxsegrtdose.category + epitxn_dose_5.category + 
                EAS + AFR +
                any_chemo_missing + any_rt_missing,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

# Get the summary of the Poisson regression model
summary_fit_all <- summary(fit_all)

results <- data.frame(
  Predictor = rownames(coef(summary_fit_all)),
  "RR (95% CI)" = sprintf("%.2f (%.2f-%.2f)",
                          exp(coef(summary_fit_all)[, "Estimate"]),
                          exp(coef(summary_fit_all)[, "Estimate"] - 1.96 * coef(summary_fit_all)[, "Std. Error"]),
                          exp(coef(summary_fit_all)[, "Estimate"] + 1.96 * coef(summary_fit_all)[, "Std. Error"])
  ),
  "P-Value" = format(coef(summary_fit_all)[, "Pr(>|z|)"], nsmall = 4)
)



MENINGIOMA.vars <- results[grepl("diagnosis|gender|maxsegrt|maxabdrtdose|pelvis|chest|neck|aa_class_dose_5|anthra|epitxn|prs", results$Predictor, ignore.case = T),]


##########################################
## Attributable fraction of Any SARCOMA ##
##########################################

# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.SARCOMA_without_lifestyle.V17b.Rdata")

# # Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
# cc
# filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
# filtered_cc
# 
# 
# 
# PHENO.ANY_SN$aa_class_dose_5.category <- as.character(PHENO.ANY_SN$aa_class_dose_5.category)
# PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "None"] <- "1st"
# PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "2nd"] <- "2nd-3rd"
# PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "3rd"] <- "2nd-3rd"
# PHENO.ANY_SN$aa_class_dose_5.category <- factor(PHENO.ANY_SN$aa_class_dose_5.category, levels = c("1st", "2nd-3rd"))


dat_all = PHENO.ANY_SN
dat_all=dat_all[dat_all$evt1==1,]

fit_all = glm(formula = event ~ Sarcoma_Machiela_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                gender + 
                aa_class_dose_5.category +
                EAS + AFR + 
                any_chemo_missing,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

summary_fit_all <- summary(fit_all)

results <- data.frame(
  Predictor = rownames(coef(summary_fit_all)),
  "RR (95% CI)" = sprintf("%.2f (%.2f-%.2f)",
                          exp(coef(summary_fit_all)[, "Estimate"]),
                          exp(coef(summary_fit_all)[, "Estimate"] - 1.96 * coef(summary_fit_all)[, "Std. Error"]),
                          exp(coef(summary_fit_all)[, "Estimate"] + 1.96 * coef(summary_fit_all)[, "Std. Error"])
  ),
  "P-Value" = format(coef(summary_fit_all)[, "Pr(>|z|)"], nsmall = 4)
)



SARCOMA.vars <- results[grepl("diagnosis|gender|maxsegrt|maxabdrtdose|pelvis|chest|neck|aa_class_dose|anthra|epitxn|prs", results$Predictor, ignore.case = T),]


## Combine all


dataframes <- list(ANY_SN=ANY_SN.vars, Meningioma=MENINGIOMA.vars, Sarcoma=SARCOMA.vars, Breast=BREASTCANCER.vars, Thyroid=THYROID.vars, NMSC=NMSC.vars, SMN=SMN.vars)
for (i in 1:length(dataframes)) {
  print(i)
  df <- dataframes[[i]]
  rownames(df) <- NULL
  # Extract the part of the column names after "PRS.tertile.category"
  df$Predictor <- sub(".*PRS.tertile.category", "PRS_tertile_", df$Predictor, ignore.case = T)
  colnames(df)[-1] <- paste0(c("RR (95% CI)", "P"), "_", names(dataframes)[i])
  dataframes[[i]] <- df
}
# dataframes[[2]]

empty_dataframe <- Reduce(function(x, y) merge(x, y, by = "Predictor", all = TRUE), dataframes)



# Print the resulting dataframe
print(empty_dataframe)

# Find the columns that start with "P_"
p_columns <- empty_dataframe[, grepl("^P_", colnames(empty_dataframe))]

# Convert the selected columns to numeric
for (col in colnames(p_columns)) {
  empty_dataframe[, col] <- as.numeric(empty_dataframe[, col])
  empty_dataframe[, col][empty_dataframe[, col] < 0.0001] <- "< .0001"
}


## Sort
empty_dataframe.ordered <- empty_dataframe[match(c(NA, "AGE_AT_DIAGNOSIS5-9", "AGE_AT_DIAGNOSIS10-14", "AGE_AT_DIAGNOSIS>=15", "genderFemale", 
                                                   NA, "maxsegrtdose.category>0-<18", "maxsegrtdose.category>=18-<30", "maxsegrtdose.category>=30", NA,
                                                   NA, "maxabdrtdose.category>0-<30", "maxabdrtdose.category>=30", NA,
                                                   NA, "maxpelvisrtdose.category>0-<20", "maxpelvisrtdose.category>=20", NA,
                                                   NA, "maxchestrtdose.category>0-<20", "maxchestrtdose.category>=20", NA, 
                                                   NA, "maxneckrtdose.category>0-<11", "maxneckrtdose.category>=11-<20", "maxneckrtdose.category>=20-<30", "maxneckrtdose.category>=30", NA,
                                                   NA, "aa_class_dose_5.category1st", "aa_class_dose_5.category2nd", "aa_class_dose_5.category3rd",
                                                   NA, "anthra_jco_dose_5.category1st",  "anthra_jco_dose_5.category2nd", "anthra_jco_dose_5.category3rd",
                                                   NA, "epitxn_dose_5.category1st", "epitxn_dose_5.category2nd", "epitxn_dose_5.category3rd",
                                                   NA, "PRS_tertile_2nd", "PRS_tertile_3rd"), empty_dataframe$Predictor),]


# Fill missing values with NA
empty_dataframe.ordered[is.na(empty_dataframe.ordered)] <- ""      



#####################################
## Compare Qi categories with ours ##
#####################################
