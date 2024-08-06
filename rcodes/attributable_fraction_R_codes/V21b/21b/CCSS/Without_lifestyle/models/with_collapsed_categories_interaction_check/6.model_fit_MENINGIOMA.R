# load ANY SN data
# load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.MENINGIOMA_without_lifestyle.V20b.Rdata")
load("Z:/ResearchHome/Groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.MENINGIOMA_without_lifestyle.V20b.Rdata")

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

# PHENO.ANY_SN$maxsegrtdose.category <- as.character(PHENO.ANY_SN$maxsegrtdose.category)
# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == "None"] <- "<30"
# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">0-<18"] <- "<30"
# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">=18-<30"] <- "<30"
# PHENO.ANY_SN$maxsegrtdose.category <- factor(PHENO.ANY_SN$maxsegrtdose.category, levels = c("<30", ">=30"))

PHENO.ANY_SN$maxsegrtdose.category <- as.character(PHENO.ANY_SN$maxsegrtdose.category)
PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == "None"] <- "<18"
PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">0-<18"] <- "<18"
# PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$maxsegrtdose.category == ">=18-<30"] <- "<18"
PHENO.ANY_SN$maxsegrtdose.category <- factor(PHENO.ANY_SN$maxsegrtdose.category, levels = c("<18", ">=18-<30", ">=30"))

table(PHENO.ANY_SN$maxsegrtdose.category[PHENO.ANY_SN$MENINGIOMA == 1])


PHENO.ANY_SN$epitxn_dose_5.category <- as.character(PHENO.ANY_SN$epitxn_dose_5.category)
PHENO.ANY_SN$epitxn_dose_5.category[PHENO.ANY_SN$epitxn_dose_5.category == "2nd"] <- "2nd-3rd"
PHENO.ANY_SN$epitxn_dose_5.category[PHENO.ANY_SN$epitxn_dose_5.category == "3rd"] <- "2nd-3rd"
PHENO.ANY_SN$epitxn_dose_5.category <- factor(PHENO.ANY_SN$epitxn_dose_5.category, levels = c("None", "1st", "2nd-3rd"))
######################################
## Attributable fraction of Any SNs ##
######################################
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


# Model with interaction between PRS and maxsegrtdose.category
fit_segrt <- glm(formula = event ~ Meningioma_PRS.tertile.category +
                   AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                   AGE_AT_DIAGNOSIS + gender + 
                   maxsegrtdose.category + epitxn_dose_5.category + 
                   EAS + AFR +
                   any_chemo_missing + any_rt_missing + 
                   Meningioma_PRS.tertile.category * maxsegrtdose.category,
                 family = "poisson", offset = log(dat_all$PY), data = dat_all)

# Model with interaction between PRS and epitxn_dose_5.category
fit_epitxnrt <- glm(formula = event ~ Meningioma_PRS.tertile.category +
                      AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                      AGE_AT_DIAGNOSIS + gender + 
                      maxsegrtdose.category + epitxn_dose_5.category + 
                      EAS + AFR +
                      any_chemo_missing + any_rt_missing + 
                      Meningioma_PRS.tertile.category * epitxn_dose_5.category,
                    family = "poisson", offset = log(dat_all$PY), data = dat_all)

# ANOVA for likelihood ratio tests
anova_segrt <- anova(fit_all, fit_segrt, test = "Chisq")
anova_epitxnrt <- anova(fit_all, fit_epitxnrt, test = "Chisq")

# Create a results table
results <- data.frame(
  Model = c("Interaction: PRS * maxsegrtdose.category", "Interaction: PRS * epitxn_dose_5.category"),
  Df = c(anova_segrt$Df[2], anova_epitxnrt$Df[2]),
  Deviance = c(anova_segrt$Deviance[2], anova_epitxnrt$Deviance[2]),
  Chisq = c(anova_segrt$`Resid. Dev`[1] - anova_segrt$`Resid. Dev`[2], 
            anova_epitxnrt$`Resid. Dev`[1] - anova_epitxnrt$`Resid. Dev`[2]),
  P_value = c(anova_segrt$`Pr(>Chi)`[2], anova_epitxnrt$`Pr(>Chi)`[2])
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

fit_all = glm(formula = event ~ Meningioma_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                AGE_AT_DIAGNOSIS + gender +
                maxsegrtdose.category + epitxn_dose_5.category +
                ancestry +
                any_chemo_missing + any_rt_missing,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

# Define the models with interactions
fit_maxsegrt_interaction <- glm(formula = event ~ 
                                  Meningioma_PRS.tertile.category +
                                  AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                                  AGE_AT_DIAGNOSIS + gender +
                                  maxsegrtdose.category + epitxn_dose_5.category +
                                  ancestry +
                                  any_chemo_missing + any_rt_missing +
                                  maxsegrtdose.category * ancestry,
                                family = "poisson", offset = log(dat_all$PY), data = dat_all)

fit_epitxn_interaction <- glm(formula = event ~ 
                                Meningioma_PRS.tertile.category +
                                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                                AGE_AT_DIAGNOSIS + gender +
                                maxsegrtdose.category + epitxn_dose_5.category +
                                ancestry +
                                any_chemo_missing + any_rt_missing +
                                epitxn_dose_5.category * ancestry,
                              family = "poisson", offset = log(dat_all$PY), data = dat_all)

# Perform likelihood ratio tests
lrt_maxsegrt <- anova(fit_all, fit_maxsegrt_interaction, test = "Chisq")
lrt_epitxn <- anova(fit_all, fit_epitxn_interaction, test = "Chisq")

# Create a results table
results <- data.frame(
  Interaction = c("ancestry * maxsegrtdose.category", 
                  "ancestry * epitxn_dose_5.category"),
  Df = c(lrt_maxsegrt$Df[2], 
         lrt_epitxn$Df[2]),
  Deviance = c(lrt_maxsegrt$Deviance[2], 
               lrt_epitxn$Deviance[2]),
  Chisq = c(lrt_maxsegrt$`Resid. Dev`[1] - lrt_maxsegrt$`Resid. Dev`[2],
            lrt_epitxn$`Resid. Dev`[1] - lrt_epitxn$`Resid. Dev`[2]),
  P_value = c(lrt_maxsegrt$`Pr(>Chi)`[2], 
              lrt_epitxn$`Pr(>Chi)`[2])
)

# Print results
print(results)

