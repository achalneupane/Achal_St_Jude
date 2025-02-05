# load ANY SN data
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/6.sjlife_without_lifestyle.THYROIDcancer.V20b.Rdata")

# Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
cc
filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
filtered_cc

## Admixture classification
admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
EUR.admix <- admixture$INDIVIDUAL[admixture$EUR > 0.8]
AFR.admix <- admixture$INDIVIDUAL[admixture$AFR > 0.6]

PHENO.ANY_SN$admixture <- NA
PHENO.ANY_SN$admixture [PHENO.ANY_SN$sjlid %in% EUR.admix] <- "EUR"
PHENO.ANY_SN$admixture [PHENO.ANY_SN$sjlid %in% AFR.admix] <- "AFR"

PHENO.ANY_SN$maxneckrtdose.category <- as.character(PHENO.ANY_SN$maxneckrtdose.category)
PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == ">0-<11"] <- ">0-<30"
PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == ">=11-<20"] <- ">0-<30"
PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$maxneckrtdose.category == ">=20-<30"] <- ">0-<30"
PHENO.ANY_SN$maxneckrtdose.category <- factor(PHENO.ANY_SN$maxneckrtdose.category, levels = c("None", ">0-<30", ">=30"))

table(PHENO.ANY_SN$maxneckrtdose.category[PHENO.ANY_SN$THYROIDcancer == 1])


PHENO.ANY_SN$epitxn_dose_5.category <- as.character(PHENO.ANY_SN$epitxn_dose_5.category)
PHENO.ANY_SN$epitxn_dose_5.category[PHENO.ANY_SN$epitxn_dose_5.category == "1st"] <- "Any"
PHENO.ANY_SN$epitxn_dose_5.category[PHENO.ANY_SN$epitxn_dose_5.category == "2nd"] <- "Any"
PHENO.ANY_SN$epitxn_dose_5.category[PHENO.ANY_SN$epitxn_dose_5.category == "3rd"] <- "Any"
PHENO.ANY_SN$epitxn_dose_5.category <- factor(PHENO.ANY_SN$epitxn_dose_5.category, levels = c("None", "Any"))

######################################
## Attributable fraction of Any SNs ##
######################################
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

THYROIDcancer.sjlife <- {}
# Model with interaction between PRS and maxneckrtdose.category
fit_neckrt <- glm(formula = event ~ Thyroid_PRS.tertile.category +
                    AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                    AGE_AT_DIAGNOSIS + gender + 
                    maxneckrtdose.category + epitxn_dose_5.category + 
                    EAS + AFR +
                    any_chemo_missing + any_rt_missing + 
                    Thyroid_PRS.tertile.category * maxneckrtdose.category,
                  family = "poisson", offset = log(dat_all$PY), data = dat_all)

coefficients <- coef(fit_neckrt)

# Extract interaction term coefficients
interaction_coefficients <- coefficients[grep(":", names(coefficients))]
# Print the interaction coefficients
cc <- as.data.frame(interaction_coefficients)
THYROIDcancer.sjlife <- rbind.data.frame(THYROIDcancer.sjlife, cc)


# Model with interaction between PRS and epitxn_dose_5.category
fit_epitxnrt <- glm(formula = event ~ Thyroid_PRS.tertile.category +
                      AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                      AGE_AT_DIAGNOSIS + gender + 
                      maxneckrtdose.category + epitxn_dose_5.category + 
                      EAS + AFR +
                      any_chemo_missing + any_rt_missing + 
                      Thyroid_PRS.tertile.category * epitxn_dose_5.category,
                    family = "poisson", offset = log(dat_all$PY), data = dat_all)

coefficients <- coef(fit_epitxnrt)

# Extract interaction term coefficients
interaction_coefficients <- coefficients[grep(":", names(coefficients))]
# Print the interaction coefficients
cc <- as.data.frame(interaction_coefficients)
THYROIDcancer.sjlife <- rbind.data.frame(THYROIDcancer.sjlife, cc)

# ANOVA for likelihood ratio tests
anova_neckrt <- anova(fit_all, fit_neckrt, test = "Chisq")
anova_epitxnrt <- anova(fit_all, fit_epitxnrt, test = "Chisq")

# Create a results table
results <- data.frame(
  Model = c("Interaction: PRS * maxneckrtdose.category", "Interaction: PRS * epitxn_dose_5.category"),
  Df = c(anova_neckrt$Df[2], anova_epitxnrt$Df[2]),
  Deviance = c(anova_neckrt$Deviance[2], anova_epitxnrt$Deviance[2]),
  Chisq = c(anova_neckrt$`Resid. Dev`[1] - anova_neckrt$`Resid. Dev`[2], 
            anova_epitxnrt$`Resid. Dev`[1] - anova_epitxnrt$`Resid. Dev`[2]),
  P_value = c(anova_neckrt$`Pr(>Chi)`[2], anova_epitxnrt$`Pr(>Chi)`[2])
)

print(results)


##############################################
## Check ancestry and treatment interaction ##
##############################################
PHENO.ANY_SN$ancestry <- "Other"
PHENO.ANY_SN$ancestry [PHENO.ANY_SN$EUR > 0.8] <- "EUR"
PHENO.ANY_SN$ancestry [PHENO.ANY_SN$AFR > 0.6] <- "AFR"
PHENO.ANY_SN$ancestry <- factor(PHENO.ANY_SN$ancestry, levels = c("EUR", "AFR", "Other"))

dat_all = PHENO.ANY_SN
dat_all=dat_all[dat_all$evt1==1,]

fit_all = glm(formula = event ~ Thyroid_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender +
                maxneckrtdose.category + epitxn_dose_5.category +
                ancestry +
                any_chemo_missing + any_rt_missing,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

# Define the models with interactions
fit_maxneck_interaction <- glm(formula = event ~ 
                                 Thyroid_PRS.tertile.category +
                                 AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                                 AGE_AT_DIAGNOSIS + gender +
                                 maxneckrtdose.category + epitxn_dose_5.category +
                                 ancestry +
                                 any_chemo_missing + any_rt_missing +
                                 maxneckrtdose.category * ancestry,
                               family = "poisson", offset = log(dat_all$PY), data = dat_all)

coefficients <- coef(fit_maxneck_interaction)

# Extract interaction term coefficients
interaction_coefficients <- coefficients[grep(":", names(coefficients))]
# Print the interaction coefficients
cc <- as.data.frame(interaction_coefficients)
THYROIDcancer.sjlife <- rbind.data.frame(THYROIDcancer.sjlife, cc)


fit_epitxn_interaction <- glm(formula = event ~ 
                                Thyroid_PRS.tertile.category +
                                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                                AGE_AT_DIAGNOSIS + gender +
                                maxneckrtdose.category + epitxn_dose_5.category +
                                ancestry +
                                any_chemo_missing + any_rt_missing +
                                epitxn_dose_5.category * ancestry,
                              family = "poisson", offset = log(dat_all$PY), data = dat_all)

coefficients <- coef(fit_epitxn_interaction)

# Extract interaction term coefficients
interaction_coefficients <- coefficients[grep(":", names(coefficients))]
# Print the interaction coefficients
cc <- as.data.frame(interaction_coefficients)
THYROIDcancer.sjlife <- rbind.data.frame(THYROIDcancer.sjlife, cc)


# Perform likelihood ratio tests
lrt_maxneck <- anova(fit_all, fit_maxneck_interaction, test = "Chisq")
lrt_epitxn <- anova(fit_all, fit_epitxn_interaction, test = "Chisq")

# Create a results table
results <- data.frame(
  Interaction = c("ancestry * maxneckrtdose.category", 
                  "ancestry * epitxn_dose_5.category"),
  Df = c(lrt_maxneck$Df[2], 
         lrt_epitxn$Df[2]),
  Deviance = c(lrt_maxneck$Deviance[2], 
               lrt_epitxn$Deviance[2]),
  Chisq = c(lrt_maxneck$`Resid. Dev`[1] - lrt_maxneck$`Resid. Dev`[2],
            lrt_epitxn$`Resid. Dev`[1] - lrt_epitxn$`Resid. Dev`[2]),
  P_value = c(lrt_maxneck$`Pr(>Chi)`[2], 
              lrt_epitxn$`Pr(>Chi)`[2])
)

# Print results
print(results)

