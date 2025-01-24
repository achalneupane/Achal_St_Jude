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


# Model with interaction between PRS and maxchestrtdose.category
fit_chestrt <- glm(formula = event ~ 
                     Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category +
                     AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + 
                     AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
                     maxchestrtdose.category + anthra_jco_dose_5.category +
                     EAS + AFR +
                     any_chemo_missing + any_rt_missing + 
                     Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category * maxchestrtdose.category,
                   family = "poisson", offset = log(dat_all$PY), data = dat_all)

# Model with interaction between PRS and anthra_jco_dose_5.category
fit_anthrart <- glm(formula = event ~ 
                      Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category +
                      AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + 
                      AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
                      maxchestrtdose.category + anthra_jco_dose_5.category +
                      EAS + AFR +
                      any_chemo_missing + any_rt_missing + 
                      Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category * anthra_jco_dose_5.category,
                    family = "poisson", offset = log(dat_all$PY), data = dat_all)

# ANOVA for likelihood ratio tests
anova_chestrt <- anova(fit_all, fit_chestrt, test = "Chisq")
anova_anthrart <- anova(fit_all, fit_anthrart, test = "Chisq")

# Create a results table
results <- data.frame(
  Model = c("Interaction: PRS * maxchestrtdose.category", "Interaction: PRS * anthra_jco_dose_5.category"),
  Df = c(anova_chestrt$Df[2], anova_anthrart$Df[2]),
  Deviance = c(anova_chestrt$Deviance[2], anova_anthrart$Deviance[2]),
  Chisq = c(anova_chestrt$`Resid. Dev`[1] - anova_chestrt$`Resid. Dev`[2], 
            anova_anthrart$`Resid. Dev`[1] - anova_anthrart$`Resid. Dev`[2]),
  P_value = c(anova_chestrt$`Pr(>Chi)`[2], anova_anthrart$`Pr(>Chi)`[2])
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
dat_all = dat_all[dat_all$gender == "Female",]
dat_all=dat_all[dat_all$evt1==1,]

fit_all = glm(formula = event ~
                Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 +
                AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
                maxchestrtdose.category + anthra_jco_dose_5.category +
                ancestry +
                any_chemo_missing + any_rt_missing,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)

# Define the models with interactions
fit_chestrt_interaction <- glm(formula = event ~ 
                                 Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category +
                                 AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 +
                                 AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
                                 maxchestrtdose.category + anthra_jco_dose_5.category +
                                 ancestry +
                                 any_chemo_missing + any_rt_missing +
                                 maxchestrtdose.category * ancestry,
                               family = "poisson", offset = log(dat_all$PY), data = dat_all)

fit_anthra_interaction <- glm(formula = event ~ 
                                Mavaddat_2019_ER_OVERALL_Breast_PRS.tertile.category +
                                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 +
                                AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 + AGE_AT_DIAGNOSIS +
                                maxchestrtdose.category + anthra_jco_dose_5.category +
                                ancestry +
                                any_chemo_missing + any_rt_missing +
                                anthra_jco_dose_5.category * ancestry,
                              family = "poisson", offset = log(dat_all$PY), data = dat_all)

# Perform likelihood ratio tests
lrt_chestrt <- anova(fit_all, fit_chestrt_interaction, test = "Chisq")
lrt_anthra <- anova(fit_all, fit_anthra_interaction, test = "Chisq")

# Create a results table
results <- data.frame(
  Interaction = c("ancestry * maxchestrtdose.category", 
                  "ancestry * anthra_jco_dose_5.category"),
  Df = c(lrt_chestrt$Df[2], 
         lrt_anthra$Df[2]),
  Deviance = c(lrt_chestrt$Deviance[2], 
               lrt_anthra$Deviance[2]),
  Chisq = c(lrt_chestrt$`Resid. Dev`[1] - lrt_chestrt$`Resid. Dev`[2],
            lrt_anthra$`Resid. Dev`[1] - lrt_anthra$`Resid. Dev`[2]),
  P_value = c(lrt_chestrt$`Pr(>Chi)`[2], 
              lrt_anthra$`Pr(>Chi)`[2])
)

# Print results
print(results)
