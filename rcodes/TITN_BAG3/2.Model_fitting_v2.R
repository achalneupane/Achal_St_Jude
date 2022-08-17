load("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/analysis_model_fitting/phenotype_carrier_status_v2.RDATA")

colnames(pheno.ccss_exp_eur)
pheno.ccss_exp_eur$CMP2plus <- factor(ifelse(pheno.ccss_exp_eur$CMP2plus == 1, 0,1))
colnames(pheno.sjlife_ttn_bag3)
pheno.sjlife_ttn_bag3$CMP <- factor(ifelse(pheno.sjlife_ttn_bag3$CMP == 1, 0,1))
#################################################
## Contingency table for fisher or Chi-sq test ##
## > 5 Chi; < 5 Fisher                         ##
#################################################

## 1. With Clinvar and Lof only in CCSS_exp
## -------------------------------
# a. TITN + BAG3 
TITN.BAG3.C.LoF <- table(pheno.ccss_exp_eur$CMP2plus, pheno.ccss_exp_eur$TITN_BAG3.clin.LoF.carrierSTATUS)
TITN.BAG3.C.LoF.test <- fisher.test(TITN.BAG3.C.LoF)
TITN.BAG3.C.LoF.test.adj <- glm(formula = CMP2plus ~  TITN_BAG3.clin.LoF.carrierSTATUS + a_dx + a_end + SEX + anth_DED + HeartAvg + PC1 + PC2 + PC3 +
                                PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
                                data = pheno.ccss_exp_eur)
# TITN.BAG3.C.LoF.test.adj <- coef(summary(TITN.BAG3.C.LoF.test.adj))[2,4]



# b. TITN only
TITN.C.LoF <- table(pheno.ccss_exp_eur$CMP2plus, pheno.ccss_exp_eur$TITN.Clinv.LoF.carrierSTATUS)
TITN.C.LoF.test <- fisher.test(TITN.C.LoF) # Got warning for chi-square
TITN.C.LoF.test.adj <- glm(formula = CMP2plus ~  TITN.Clinv.LoF.carrierSTATUS + a_dx + a_end + SEX + anth_DED + HeartAvg + PC1 + PC2 + PC3 +
                                  PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
                                data = pheno.ccss_exp_eur)

# c. BAG3 only
BAG3.C.LoF <- table(pheno.ccss_exp_eur$CMP2plus, pheno.ccss_exp_eur$BAG3.Clinv.LoF.carrierSTATUS)
BAG3.C.LoF.test <- fisher.test(BAG3.C.LoF) # Got error; no carriers
BAG3.C.LoF.test.adj <- glm(formula = CMP2plus ~  BAG3.Clinv.LoF.carrierSTATUS + a_dx + a_end + SEX + anth_DED + HeartAvg + PC1 + PC2 + PC3 +
                             PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
                           data = pheno.ccss_exp_eur)

# ## 2. With Clinvar, Lof, REVEL in CCSS_exp
# ## -------------------------------
# # a. TITN + BAG3 
# TITN.BAG3.C.LoF.REVEL <- table(pheno.ccss_exp_eur$CMP2plus, pheno.ccss_exp_eur$TITN_BAG3.clin.LoF.REVEL.carrierSTATUS)
# TITN.BAG3.C.LoF.REVEL.test <- fisher.test(TITN.BAG3.C.LoF.REVEL)
# TITN.BAG3.C.LoF.REVEL.test.adj <- glm(formula = CMP2plus ~  TITN_BAG3.clin.LoF.REVEL.carrierSTATUS + a_dx + a_end + SEX + anth_DED + HeartAvg + PC1 + PC2 + PC3 +
#                              PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
#                            data = pheno.ccss_exp_eur)
# 
# # b. TITN only
# TITN.C.LoF.REVEL <- table(pheno.ccss_exp_eur$CMP2plus, pheno.ccss_exp_eur$TITN.Clinv.LoF.REVEL.carrierSTATUS)
# TITN.C.LoF.REVEL.test <- fisher.test(TITN.C.LoF.REVEL) 
# TITN.C.LoF.REVEL.test.adj <- glm(formula = CMP2plus ~  TITN.Clinv.LoF.REVEL.carrierSTATUS + a_dx + a_end + SEX + anth_DED + HeartAvg + PC1 + PC2 + PC3 +
#                                         PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
#                                       data = pheno.ccss_exp_eur)

# # c. BAG3 only
# BAG3.C.LoF.REVEL <- table(pheno.ccss_exp_eur$CMP2plus, pheno.ccss_exp_eur$BAG3.Clinv.LoF.REVEL.carrierSTATUS)
# BAG3.C.LoF.REVEL.test <- fisher.test(BAG3.C.LoF.REVEL) # Got warning for chi-square
# BAG3.C.LoF.REVEL.test.adj <- glm(formula = CMP2plus ~  BAG3.Clinv.LoF.REVEL.carrierSTATUS + a_dx + a_end + SEX + anth_DED + HeartAvg + PC1 + PC2 + PC3 +
#                                    PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
#                                  data = pheno.ccss_exp_eur)

## 3. With Clinvar and Lof only in SJLIFE
## -------------------------------
# a. TITN + BAG3 
TITN.BAG3.C.LoF.sjlife <- table(pheno.sjlife_ttn_bag3$CMP, pheno.sjlife_ttn_bag3$TITN_BAG3.clin.LoF.carrierSTATUS)
TITN.BAG3.C.LoF.test.sjlife <- fisher.test(TITN.BAG3.C.LoF.sjlife)
TITN.BAG3.C.LoF.test.sjlife.adj <- glm(formula = CMP ~  TITN_BAG3.clin.LoF.carrierSTATUS + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                                   PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
                                 data = pheno.sjlife_ttn_bag3)

# b. TITN only
TITN.C.LoF.sjlife <- table(pheno.sjlife_ttn_bag3$CMP, pheno.sjlife_ttn_bag3$TITN.Clinv.LoF.carrierSTATUS)
TITN.C.LoF.test.sjlife <- fisher.test(TITN.C.LoF.sjlife) # Got warning for chi-square
TITN.C.LoF.test.sjlife.adj <- glm(formula = CMP ~  TITN.Clinv.LoF.carrierSTATUS + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                                         PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
                                       data = pheno.sjlife_ttn_bag3)

# c. BAG3 only
BAG3.C.LoF.sjlife <- table(pheno.sjlife_ttn_bag3$CMP, pheno.sjlife_ttn_bag3$BAG3.Clinv.LoF.carrierSTATUS)
BAG3.C.LoF.test.sjlife <- fisher.test(BAG3.C.LoF.sjlife) # Got error; no carriers
BAG3.C.LoF.test.sjlife.adj <- glm(formula = CMP ~  BAG3.Clinv.LoF.carrierSTATUS + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                                    PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
                                  data = pheno.sjlife_ttn_bag3)

# ## 4. With Clinvar, Lof, REVEL in SJLIFE
# ## -------------------------------
# # a. TITN + BAG3 
# TITN.BAG3.C.LoF.REVEL.sjlife <- table(pheno.sjlife_ttn_bag3$CMP, pheno.sjlife_ttn_bag3$TITN_BAG3.clin.LoF.REVEL.carrierSTATUS)
# TITN.BAG3.C.LoF.REVEL.test.sjlife <- fisher.test(TITN.BAG3.C.LoF.REVEL.sjlife)
# TITN.BAG3.C.LoF.REVEL.test.sjlife.adj <- glm(formula = CMP ~  TITN_BAG3.clin.LoF.REVEL.carrierSTATUS + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
#                                     PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
#                                   data = pheno.sjlife_ttn_bag3)
# 
# # b. TITN only
# TITN.C.LoF.REVEL.sjlife <- table(pheno.sjlife_ttn_bag3$CMP, pheno.sjlife_ttn_bag3$TITN.Clinv.LoF.REVEL.carrierSTATUS)
# TITN.C.LoF.REVEL.test.sjlife <- fisher.test(TITN.C.LoF.REVEL.sjlife) 
# TITN.C.LoF.REVEL.test.sjlife.adj <- glm(formula = CMP ~  TITN.Clinv.LoF.REVEL.carrierSTATUS + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
#                                           PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
#                                         data = pheno.sjlife_ttn_bag3)
# 
# # c. BAG3 only
# BAG3.C.LoF.REVEL.sjlife <- table(pheno.sjlife_ttn_bag3$CMP, pheno.sjlife_ttn_bag3$BAG3.Clinv.LoF.REVEL.carrierSTATUS)
# BAG3.C.LoF.REVEL.test.sjlife <- fisher.test(BAG3.C.LoF.REVEL.sjlife) # Got warning for chi-square
# BAG3.C.LoF.REVEL.test.sjlife.adj <- glm(formula = CMP ~  BAG3.Clinv.LoF.REVEL.carrierSTATUS + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
#                                           PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
#                                         data = pheno.sjlife_ttn_bag3)

#############
## Results ##
#############
# Clinvar and LoF only
cbind.data.frame(setNames(rbind.data.frame(c(TITN.BAG3.C.LoF.test$p.value, coef(summary(TITN.BAG3.C.LoF.test.adj))[2,4]),
c(TITN.C.LoF.test$p.value, coef(summary(TITN.C.LoF.test.adj))[2,4]),
c(NA, NA)
), c("Unadj.CCSS", "adj.CCSS")), 
setNames(rbind.data.frame(c(TITN.BAG3.C.LoF.test.sjlife$p.value, coef(summary(TITN.BAG3.C.LoF.test.sjlife.adj))[2,4]),
                          c(TITN.C.LoF.test.sjlife$p.value, coef(summary(TITN.C.LoF.test.sjlife.adj))[2,4]),
                          c(NA, NA)
), c("Unadj.SJLIFE", "adj.SJLIFE")))


####################################################################################
## We have decided to fit logistic P with P/LP from clinvar and LoF for the final ##
## table. We may use logistic with P/LP from clinvar, LoF and REVEL as            ##
## supplementary table.                                                           ##
####################################################################################
CI1 <- paste0(" (",round(exp(c(coef(summary(TITN.BAG3.C.LoF.test.adj))[2,1]))-1.96*as.data.frame(summary(TITN.BAG3.C.LoF.test.adj)$coefficients)[c(1,2)][2,2],2), " to ", round(exp(c(coef(summary(TITN.BAG3.C.LoF.test.adj))[2,1]))+1.96*as.data.frame(summary(TITN.BAG3.C.LoF.test.adj)$coefficients)[c(1,2)][2,2], 2),")")
CI2 <- paste0(" (",round(exp(c(coef(summary(TITN.C.LoF.test.adj))[2,1]))-1.96*as.data.frame(summary(TITN.C.LoF.test.adj)$coefficients)[c(1,2)][2,2],2), " to ", round(exp(c(coef(summary(TITN.C.LoF.test.adj))[2,1]))+1.96*as.data.frame(summary(TITN.C.LoF.test.adj)$coefficients)[c(1,2)][2,2], 2),")")
# CI3 <- paste0(" (",round(exp(c(coef(summary(BAG3.C.LoF.test.adj))[2,1]))-1.96*as.data.frame(summary(BAG3.C.LoF.test.adj)$coefficients)[c(1,2)][2,2],2), " to ", round(exp(c(coef(summary(BAG3.C.LoF.test.adj))[2,1]))+1.96*as.data.frame(summary(BAG3.C.LoF.test.adj)$coefficients)[c(1,2)][2,2], 2),")")

CCSS.Clin.LoF.df <- cbind.data.frame(setNames(rbind.data.frame(c("TTN_BAG3", TITN.BAG3.C.LoF[1,2]/sum(TITN.BAG3.C.LoF[1,]), TITN.BAG3.C.LoF[2,2]/sum(TITN.BAG3.C.LoF[2,]), paste0(round(exp(c(coef(summary(TITN.BAG3.C.LoF.test.adj))[2,1])),2), CI1), coef(summary(TITN.BAG3.C.LoF.test.adj))[2,4]),
                                           c("TTN", TITN.C.LoF[1,2]/sum(TITN.C.LoF[1,]), TITN.C.LoF[2,2]/sum(TITN.C.LoF[2,]), paste0(round(exp(c(coef(summary(TITN.C.LoF.test.adj))[2,1])),2), CI2), coef(summary(TITN.C.LoF.test.adj))[2,4]),
                                           c("BAG3", 0, 0, NA, NA)
), c("Gene", "prevalence.in.CCSS.exp.CO", "prevalence.in.CCSS.exp.CA","OR.CCSS (CI)", "adj.P.CCSS")))

write.table(CCSS.Clin.LoF.df, "CCSS.Clin.LoF.df.txt", sep = "\t", col.names = T, row.names = F, quote = F)


CI1 <- paste0(" (",round(exp(c(coef(summary(TITN.BAG3.C.LoF.test.sjlife.adj))[2,1]))-1.96*as.data.frame(summary(TITN.BAG3.C.LoF.test.sjlife.adj)$coefficients)[c(1,2)][2,2],2), " to ", round(exp(c(coef(summary(TITN.BAG3.C.LoF.test.sjlife.adj))[2,1]))+1.96*as.data.frame(summary(TITN.BAG3.C.LoF.test.sjlife.adj)$coefficients)[c(1,2)][2,2], 2),")")
CI2 <- paste0(" (",round(exp(c(coef(summary(TITN.C.LoF.test.sjlife.adj))[2,1]))-1.96*as.data.frame(summary(TITN.C.LoF.test.sjlife.adj)$coefficients)[c(1,2)][2,2],2), " to ", round(exp(c(coef(summary(TITN.C.LoF.test.sjlife.adj))[2,1]))+1.96*as.data.frame(summary(TITN.C.LoF.test.sjlife.adj)$coefficients)[c(1,2)][2,2], 2),")")
# CI3 <- paste0(" (",round(exp(c(coef(summary(BAG3.C.LoF.test.sjlife.adj))[2,1]))-1.96*as.data.frame(summary(BAG3.C.LoF.test.sjlife.adj)$coefficients)[c(1,2)][2,2],2), " to ", round(exp(c(coef(summary(BAG3.C.LoF.test.sjlife.adj))[2,1]))+1.96*as.data.frame(summary(BAG3.C.LoF.test.sjlife.adj)$coefficients)[c(1,2)][2,2], 2),")")

# prop.test(TITN.BAG3.C.LoF.sjlife[1,2],sum(TITN.BAG3.C.LoF.sjlife[1,]))
SJLIFE.Clin.LoF.df <- cbind.data.frame(setNames(rbind.data.frame(c("TTN_BAG3", TITN.BAG3.C.LoF.sjlife[1,2]/sum(TITN.BAG3.C.LoF.sjlife[1,]), TITN.BAG3.C.LoF.sjlife[2,2]/sum(TITN.BAG3.C.LoF.sjlife[2,]), paste0(round(exp(c(coef(summary(TITN.BAG3.C.LoF.test.sjlife.adj))[2,1])),2), CI1), coef(summary(TITN.BAG3.C.LoF.test.sjlife.adj))[2,4]),
                                                               c("TTN", TITN.C.LoF.sjlife[1,2]/sum(TITN.C.LoF.sjlife[1,]), TITN.C.LoF.sjlife[2,2]/sum(TITN.C.LoF.sjlife[2,]), paste0(round(exp(c(coef(summary(TITN.C.LoF.test.sjlife.adj))[2,1])),2), CI2), coef(summary(TITN.C.LoF.test.sjlife.adj))[2,4]),
                                                               c("BAG3", 0, 0, NA, NA)
), c("Gene", "prevalence.in.SJLIFE.exp.CO", "prevalence.in.SJLIFE.exp.CA","OR.SJLIFE (CI)", "adj.P.SJLIFE")))

write.table(SJLIFE.Clin.LoF.df, "SJLIFE.Clin.LoF.df.txt", sep = "\t", col.names = T, row.names = F, quote = F)


save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/analysis_model_fitting/2.Model_fitting_v2.RDATA")
