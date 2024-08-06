# load ANY SN data
# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.NMSCs_without_lifestyle.V20b.Rdata")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.NMSCs_without_lifestyle.V20b.Rdata")

# Yutaka's email on 03/16/2023:  It seems maxsegrtdose 0-18 Gy is a very small group and perhaps needs to be combined with 18-30 Gy
cc
filtered_cc <- cc[cc[, 2] < 10 | cc[, 3] < 10, 1]
filtered_cc

## Admixture classification
admixture <- read.table("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/admixture/merged.ancestry.file.txt", header = T)
admixture$INDIVIDUAL <- sapply(strsplit(admixture$INDIVIDUAL,"_"), `[`, 1)
EUR.admix <- admixture$INDIVIDUAL[admixture$EUR > 0.8]
AFR.admix <- admixture$INDIVIDUAL[admixture$AFR > 0.6]

PHENO.ANY_SN$admixture <- NA
PHENO.ANY_SN$admixture [PHENO.ANY_SN$ccssid %in% EUR.admix] <- "EUR"
PHENO.ANY_SN$admixture [PHENO.ANY_SN$ccssid %in% AFR.admix] <- "AFR"

PHENO.ANY_SN$maxsegrtdose.category <- as.character(PHENO.ANY_SN$maxsegrtdose.category)
PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">0-<18"] <- ">0-<30"
PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">=18-<30"] <- ">0-<30"
PHENO.ANY_SN$maxsegrtdose.category <- factor(PHENO.ANY_SN$maxsegrtdose.category, levels = c("None", ">0-<30", ">=30"))

PHENO.ANY_SN$maxpelvisrtdose.category <- as.character(PHENO.ANY_SN$maxpelvisrtdose.category)
PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$maxpelvisrtdose.category == ">0-<20"] <- "Any"
PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$maxpelvisrtdose.category == ">=20"] <- "Any"
PHENO.ANY_SN$maxpelvisrtdose.category <- factor(PHENO.ANY_SN$maxpelvisrtdose.category, levels = c("None", "Any"))

# table(PHENO.ANY_SN$maxpelvisrtdose.category[PHENO.ANY_SN$NMSCs == 0])

######################################
## Attributable fraction of Any SNs ##
######################################
dat_all = PHENO.ANY_SN

dat_all=dat_all[dat_all$evt1==1,]
fit_all = glm(formula = event ~ BASALcell_PRS.tertile.category + 
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category +
                EAS + AFR +
                any_rt_missing,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

# Load necessary libraries
library(stats)

# Model with interaction between PRS and maxsegrtdose.category
fit_segrt <- glm(formula = event ~ BASALcell_PRS.tertile.category + 
                   AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                   gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category +
                   EAS + AFR +
                   any_rt_missing + 
                   BASALcell_PRS.tertile.category * maxsegrtdose.category,
                 family = "poisson", offset = log(dat_all$PY), data = dat_all)

# Model with interaction between PRS and maxabdrtdose.category
fit_abdrt <- glm(formula = event ~ BASALcell_PRS.tertile.category + 
                   AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                   gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category +
                   EAS + AFR +
                   any_rt_missing + 
                   BASALcell_PRS.tertile.category * maxabdrtdose.category,
                 family = "poisson", offset = log(dat_all$PY), data = dat_all)

# Model with interaction between PRS and maxpelvisrtdose.category
fit_pelvis <- glm(formula = event ~ BASALcell_PRS.tertile.category + 
                    AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                    gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category +
                    EAS + AFR +
                    any_rt_missing + 
                    BASALcell_PRS.tertile.category * maxpelvisrtdose.category,
                  family = "poisson", offset = log(dat_all$PY), data = dat_all)

# ANOVA for likelihood ratio tests
anova_segrt <- anova(fit_all, fit_segrt, test = "Chisq")
anova_abdrt <- anova(fit_all, fit_abdrt, test = "Chisq")
anova_pelvis <- anova(fit_all, fit_pelvis, test = "Chisq")

# Create a results table
results <- data.frame(
  Model = c("Interaction: PRS * maxsegrtdose.category", "Interaction: PRS * maxabdrtdose.category", 
            "Interaction: PRS * maxpelvisrtdose.category"),
  Df = c(anova_segrt$Df[2], anova_abdrt$Df[2], anova_pelvis$Df[2]),
  Deviance = c(anova_segrt$Deviance[2], anova_abdrt$Deviance[2], anova_pelvis$Deviance[2]),
  Chisq = c(anova_segrt$`Resid. Dev`[1] - anova_segrt$`Resid. Dev`[2], 
            anova_abdrt$`Resid. Dev`[1] - anova_abdrt$`Resid. Dev`[2], 
            anova_pelvis$`Resid. Dev`[1] - anova_pelvis$`Resid. Dev`[2]),
  P_value = c(anova_segrt$`Pr(>Chi)`[2], anova_abdrt$`Pr(>Chi)`[2], anova_pelvis$`Pr(>Chi)`[2])
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
fit_all = glm(formula = event ~ BASALcell_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category +
                ancestry +
                any_rt_missing,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)


# Define the models with interactions
fit_segrt_interaction <- glm(formula = event ~ BASALcell_PRS.tertile.category +
                               AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                               gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category +
                               ancestry +
                               any_rt_missing +
                               maxsegrtdose.category * ancestry,
                             family = "poisson", offset = log(dat_all$PY), data = dat_all)

fit_abdrt_interaction <- glm(formula = event ~ BASALcell_PRS.tertile.category +
                               AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                               gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category +
                               ancestry +
                               any_rt_missing +
                               maxabdrtdose.category * ancestry,
                             family = "poisson", offset = log(dat_all$PY), data = dat_all)

fit_pelvisrt_interaction <- glm(formula = event ~ BASALcell_PRS.tertile.category +
                                  AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                                  gender + maxsegrtdose.category + maxabdrtdose.category + maxpelvisrtdose.category +
                                  ancestry +
                                  any_rt_missing +
                                  maxpelvisrtdose.category * ancestry,
                                family = "poisson", offset = log(dat_all$PY), data = dat_all)

# Perform likelihood ratio tests
lrt_segrt <- anova(fit_all, fit_segrt_interaction, test = "Chisq")
lrt_abdrt <- anova(fit_all, fit_abdrt_interaction, test = "Chisq")
lrt_pelvisrt <- anova(fit_all, fit_pelvisrt_interaction, test = "Chisq")

# Create a results table
results <- data.frame(
  Interaction = c("ancestry * maxsegrtdose.category", 
                  "ancestry * maxabdrtdose.category", 
                  "ancestry * maxpelvisrtdose.category"),
  Df = c(lrt_segrt$Df[2], 
         lrt_abdrt$Df[2], 
         lrt_pelvisrt$Df[2]),
  Deviance = c(lrt_segrt$Deviance[2], 
               lrt_abdrt$Deviance[2], 
               lrt_pelvisrt$Deviance[2]),
  Chisq = c(lrt_segrt$`Resid. Dev`[1] - lrt_segrt$`Resid. Dev`[2],
            lrt_abdrt$`Resid. Dev`[1] - lrt_abdrt$`Resid. Dev`[2],
            lrt_pelvisrt$`Resid. Dev`[1] - lrt_pelvisrt$`Resid. Dev`[2]),
  P_value = c(lrt_segrt$`Pr(>Chi)`[2], 
              lrt_abdrt$`Pr(>Chi)`[2], 
              lrt_pelvisrt$`Pr(>Chi)`[2])
)

# Print results
print(results)
