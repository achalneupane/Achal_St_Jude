# load ANY SN data
load("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/6.sjlife_without_lifestyle.Any_SNs.V20b.Rdata")

# Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
cc
filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
filtered_cc

## Admixture classification
admixture <- read.table("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
EUR.admix <- admixture$INDIVIDUAL[admixture$EUR > 0.8]
AFR.admix <- admixture$INDIVIDUAL[admixture$AFR > 0.6]

PHENO.ANY_SN$admixture <- NA
PHENO.ANY_SN$admixture [PHENO.ANY_SN$sjlid %in% EUR.admix] <- "EUR"
PHENO.ANY_SN$admixture [PHENO.ANY_SN$sjlid %in% AFR.admix] <- "AFR"

PHENO.ANY_SN$maxsegrtdose.category <- as.character(PHENO.ANY_SN$maxsegrtdose.category)
PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">0-<18"] <- ">0-<30"
PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">=18-<30"] <- ">0-<30"
PHENO.ANY_SN$maxsegrtdose.category <- factor(PHENO.ANY_SN$maxsegrtdose.category, levels = c("None", ">0-<30", ">=30"))


######################################
## Attributable fraction of Any SNs ##
######################################
dat_all = PHENO.ANY_SN
dat_all=dat_all[dat_all$evt1==1,]

# Null model
fit_all = glm(formula = event ~ Pleiotropy_PRSWEB_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender + 
                maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category +
                EAS + AFR + 
                any_chemo_missing + any_rt_missing,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

# Model with interaction between PRS and maxsegrtdose.category
fit_segrt <- glm(formula = event ~ Pleiotropy_PRSWEB_PRS.tertile.category +
                   AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                   AGE_AT_DIAGNOSIS + gender + 
                   maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category +
                   EAS + AFR + 
                   any_chemo_missing + any_rt_missing + 
                   Pleiotropy_PRSWEB_PRS.tertile.category * maxsegrtdose.category,
                 family = "poisson", offset = log(dat_all$PY), data = dat_all)

# Model with interaction between PRS and maxabdrtdose.category
fit_abdrt <- glm(formula = event ~ Pleiotropy_PRSWEB_PRS.tertile.category +
                   AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                   AGE_AT_DIAGNOSIS + gender + 
                   maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category +
                   EAS + AFR + 
                   any_chemo_missing + any_rt_missing + 
                   Pleiotropy_PRSWEB_PRS.tertile.category * maxabdrtdose.category,
                 family = "poisson", offset = log(dat_all$PY), data = dat_all)

# Model with interaction between PRS and maxchestrtdose.category
fit_chest <- glm(formula = event ~ Pleiotropy_PRSWEB_PRS.tertile.category +
                   AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                   AGE_AT_DIAGNOSIS + gender + 
                   maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category +
                   EAS + AFR + 
                   any_chemo_missing + any_rt_missing + 
                   Pleiotropy_PRSWEB_PRS.tertile.category * maxchestrtdose.category,
                 family = "poisson", offset = log(dat_all$PY), data = dat_all)

# Model with interaction between PRS and epitxn_dose_5.category
fit_epitxn <- glm(formula = event ~ Pleiotropy_PRSWEB_PRS.tertile.category +
                    AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                    AGE_AT_DIAGNOSIS + gender + 
                    maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category +
                    EAS + AFR + 
                    any_chemo_missing + any_rt_missing + 
                    Pleiotropy_PRSWEB_PRS.tertile.category * epitxn_dose_5.category,
                  family = "poisson", offset = log(dat_all$PY), data = dat_all)

# ANOVA for likelihood ratio tests
anova_segrt <- anova(fit_all, fit_segrt, test = "Chisq")
anova_abdrt <- anova(fit_all, fit_abdrt, test = "Chisq")
anova_chest <- anova(fit_all, fit_chest, test = "Chisq")
anova_epitxn <- anova(fit_all, fit_epitxn, test = "Chisq")

# Create a results table
results <- data.frame(
  Model = c("Interaction: PRS * maxsegrtdose.category", "Interaction: PRS * maxabdrtdose.category", 
            "Interaction: PRS * maxchestrtdose.category", "Interaction: PRS * epitxn_dose_5.category"),
  Df = c(anova_segrt$Df[2], anova_abdrt$Df[2], anova_chest$Df[2], anova_epitxn$Df[2]),
  Deviance = c(anova_segrt$Deviance[2], anova_abdrt$Deviance[2], anova_chest$Deviance[2], anova_epitxn$Deviance[2]),
  Chisq = c(anova_segrt$`Resid. Dev`[1] - anova_segrt$`Resid. Dev`[2], 
            anova_abdrt$`Resid. Dev`[1] - anova_abdrt$`Resid. Dev`[2], 
            anova_chest$`Resid. Dev`[1] - anova_chest$`Resid. Dev`[2], 
            anova_epitxn$`Resid. Dev`[1] - anova_epitxn$`Resid. Dev`[2]),
  P_value = c(anova_segrt$`Pr(>Chi)`[2], anova_abdrt$`Pr(>Chi)`[2], anova_chest$`Pr(>Chi)`[2], anova_epitxn$`Pr(>Chi)`[2])
)

print(results)



##############################################
## Check ancestry and treatment interaction ##
##############################################
PHENO.ANY_SN$ancestry <- "Other"
PHENO.ANY_SN$ancestry [PHENO.ANY_SN$EUR > 0.8] <- "EUR"
PHENO.ANY_SN$ancestry [PHENO.ANY_SN$AFR > 0.6] <- "AFR"
PHENO.ANY_SN$ancestry <- factor(PHENO.ANY_SN$ancestry, levels = c("EUR", "AFR", "Other"))

dat_all=PHENO.ANY_SN[PHENO.ANY_SN$evt1==1,]
fit_all = glm(formula = event ~ Pleiotropy_PRSWEB_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender +
                maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category +
                ancestry +
                any_chemo_missing + any_rt_missing,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

# Interaction models
fit_segrt_interaction <- glm(formula = event ~ Pleiotropy_PRSWEB_PRS.tertile.category +
                               AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                               AGE_AT_DIAGNOSIS + gender +
                               maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category +
                               ancestry +
                               any_chemo_missing + any_rt_missing + 
                               ancestry * maxsegrtdose.category,
                             family = "poisson", offset = log(dat_all$PY), data = dat_all)

fit_abdrt_interaction <- glm(formula = event ~ Pleiotropy_PRSWEB_PRS.tertile.category +
                               AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                               AGE_AT_DIAGNOSIS + gender +
                               maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category +
                               ancestry +
                               any_chemo_missing + any_rt_missing + 
                               ancestry * maxabdrtdose.category,
                             family = "poisson", offset = log(dat_all$PY), data = dat_all)

fit_chestrt_interaction <- glm(formula = event ~ Pleiotropy_PRSWEB_PRS.tertile.category +
                                 AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                                 AGE_AT_DIAGNOSIS + gender +
                                 maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category +
                                 ancestry +
                                 any_chemo_missing + any_rt_missing + 
                                 ancestry * maxchestrtdose.category,
                               family = "poisson", offset = log(dat_all$PY), data = dat_all)

fit_epitxn_interaction <- glm(formula = event ~ Pleiotropy_PRSWEB_PRS.tertile.category +
                                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                                AGE_AT_DIAGNOSIS + gender +
                                maxsegrtdose.category + maxabdrtdose.category + maxchestrtdose.category + epitxn_dose_5.category +
                                ancestry +
                                any_chemo_missing + any_rt_missing + 
                                ancestry * epitxn_dose_5.category,
                              family = "poisson", offset = log(dat_all$PY), data = dat_all)

# Perform likelihood ratio tests
lrt_segrt <- anova(fit_all, fit_segrt_interaction, test = "Chisq")
lrt_abdrt <- anova(fit_all, fit_abdrt_interaction, test = "Chisq")
lrt_chestrt <- anova(fit_all, fit_chestrt_interaction, test = "Chisq")
lrt_epitxn <- anova(fit_all, fit_epitxn_interaction, test = "Chisq")

results <- data.frame(
  Interaction = c("ancestry * maxsegrtdose.category", 
                  "ancestry * maxabdrtdose.category", 
                  "ancestry * maxchestrtdose.category", 
                  "ancestry * epitxn_dose_5.category"),
  Df = c(lrt_segrt$Df[2], 
         lrt_abdrt$Df[2], 
         lrt_chestrt$Df[2], 
         lrt_epitxn$Df[2]),
  Deviance = c(lrt_segrt$Deviance[2], 
               lrt_abdrt$Deviance[2], 
               lrt_chestrt$Deviance[2], 
               lrt_epitxn$Deviance[2]),
  Chisq = c(lrt_segrt$`Resid. Dev`[1] - lrt_segrt$`Resid. Dev`[2],
            lrt_abdrt$`Resid. Dev`[1] - lrt_abdrt$`Resid. Dev`[2],
            lrt_chestrt$`Resid. Dev`[1] - lrt_chestrt$`Resid. Dev`[2],
            lrt_epitxn$`Resid. Dev`[1] - lrt_epitxn$`Resid. Dev`[2]),
  P_value = c(lrt_segrt$`Pr(>Chi)`[2], 
              lrt_abdrt$`Pr(>Chi)`[2], 
              lrt_chestrt$`Pr(>Chi)`[2], 
              lrt_epitxn$`Pr(>Chi)`[2])
)

print(results)
