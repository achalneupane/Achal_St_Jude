# load ANY SN data
# load("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.SARCOMA_without_lifestyle.V20b.Rdata")
load("/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/PHENOTYPE/ccss.SARCOMA_without_lifestyle.V20b.Rdata")

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

PHENO.ANY_SN$aa_class_dose_5.category <- as.character(PHENO.ANY_SN$aa_class_dose_5.category)
PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "None"] <- "1st"
PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "2nd"] <- "2nd-3rd"
PHENO.ANY_SN$aa_class_dose_5.category[PHENO.ANY_SN$aa_class_dose_5.category == "3rd"] <- "2nd-3rd"
PHENO.ANY_SN$aa_class_dose_5.category <- factor(PHENO.ANY_SN$aa_class_dose_5.category, levels = c("1st", "2nd-3rd"))

######################################
## Attributable fraction of Any SNs ##
######################################
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

# Model with interaction between PRS and aa_class_dose_5.category
fit_interaction <- glm(formula = event ~ Sarcoma_Machiela_PRS.tertile.category +
                         AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                         gender + 
                         aa_class_dose_5.category +
                         EAS + AFR +
                         any_chemo_missing + 
                         Sarcoma_Machiela_PRS.tertile.category * aa_class_dose_5.category,
                       family = "poisson", offset = log(dat_all$PY), data = dat_all)

# ANOVA for likelihood ratio test
anova_result <- anova(fit_all, fit_interaction, test = "Chisq")

# Create a results table
results <- data.frame(
  Model = "Interaction: PRS * aa_class_dose_5.category",
  Df = anova_result$Df[2],
  Deviance = anova_result$Deviance[2],
  Chisq = anova_result$`Resid. Dev`[1] - anova_result$`Resid. Dev`[2],
  P_value = anova_result$`Pr(>Chi)`[2]
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

fit_all = glm(formula = event ~ Sarcoma_Machiela_PRS.tertile.category +
                AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                gender +
                aa_class_dose_5.category +
                ancestry +
                any_chemo_missing,
              family = "poisson", offset = log(dat_all$PY), data = dat_all)

summary(fit_all)


# Define the model with interaction
fit_aa_class_interaction <- glm(formula = event ~ 
                                  Sarcoma_Machiela_PRS.tertile.category +
                                  AGE_AT_LAST_CONTACT.cs1 + AGE_AT_LAST_CONTACT.cs2 + AGE_AT_LAST_CONTACT.cs3 + AGE_AT_LAST_CONTACT.cs4 +
                                  gender +
                                  aa_class_dose_5.category +
                                  ancestry +
                                  any_chemo_missing +
                                  aa_class_dose_5.category * ancestry,
                                family = "poisson", offset = log(dat_all$PY), data = dat_all)

# Perform likelihood ratio test
lrt_aa_class <- anova(fit_all, fit_aa_class_interaction, test = "Chisq")

# Create a results table
results <- data.frame(
  Interaction = "ancestry * aa_class_dose_5.category",
  Df = lrt_aa_class$Df[2],
  Deviance = lrt_aa_class$Deviance[2],
  Chisq = lrt_aa_class$`Resid. Dev`[1] - lrt_aa_class$`Resid. Dev`[2],
  P_value = lrt_aa_class$`Pr(>Chi)`[2]
)

# Print results
print(results)
