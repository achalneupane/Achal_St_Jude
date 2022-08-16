load("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/analysis_model_fitting/phenotype_carrier_status.RDATA")

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

## 2. With Clinvar, Lof, REVEL in CCSS_exp
## -------------------------------
# a. TITN + BAG3 
TITN.BAG3.C.LoF.REVEL <- table(pheno.ccss_exp_eur$CMP2plus, pheno.ccss_exp_eur$TITN_BAG3.clin.LoF.REVEL.carrierSTATUS)
TITN.BAG3.C.LoF.REVEL.test <- fisher.test(TITN.BAG3.C.LoF.REVEL)
TITN.BAG3.C.LoF.REVEL.test.adj <- glm(formula = CMP2plus ~  TITN_BAG3.clin.LoF.REVEL.carrierSTATUS + a_dx + a_end + SEX + anth_DED + HeartAvg + PC1 + PC2 + PC3 +
                             PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
                           data = pheno.ccss_exp_eur)

# b. TITN only
TITN.C.LoF.REVEL <- table(pheno.ccss_exp_eur$CMP2plus, pheno.ccss_exp_eur$TITN.Clinv.LoF.REVEL.carrierSTATUS)
TITN.C.LoF.REVEL.test <- fisher.test(TITN.C.LoF.REVEL) 
TITN.C.LoF.REVEL.test.adj <- glm(formula = CMP2plus ~  TITN.Clinv.LoF.REVEL.carrierSTATUS + a_dx + a_end + SEX + anth_DED + HeartAvg + PC1 + PC2 + PC3 +
                                        PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
                                      data = pheno.ccss_exp_eur)

# c. BAG3 only
BAG3.C.LoF.REVEL <- table(pheno.ccss_exp_eur$CMP2plus, pheno.ccss_exp_eur$BAG3.Clinv.LoF.REVEL.carrierSTATUS)
BAG3.C.LoF.REVEL.test <- fisher.test(BAG3.C.LoF.REVEL) # Got warning for chi-square
BAG3.C.LoF.REVEL.test.adj <- glm(formula = CMP2plus ~  BAG3.Clinv.LoF.REVEL.carrierSTATUS + a_dx + a_end + SEX + anth_DED + HeartAvg + PC1 + PC2 + PC3 +
                                   PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
                                 data = pheno.ccss_exp_eur)

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

## 4. With Clinvar, Lof, REVEL in SJLIFE
## -------------------------------
# a. TITN + BAG3 
TITN.BAG3.C.LoF.REVEL.sjlife <- table(pheno.sjlife_ttn_bag3$CMP, pheno.sjlife_ttn_bag3$TITN_BAG3.clin.LoF.REVEL.carrierSTATUS)
TITN.BAG3.C.LoF.REVEL.test.sjlife <- fisher.test(TITN.BAG3.C.LoF.REVEL.sjlife)
TITN.BAG3.C.LoF.REVEL.test.sjlife.adj <- glm(formula = CMP ~  TITN_BAG3.clin.LoF.REVEL.carrierSTATUS + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                                    PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
                                  data = pheno.sjlife_ttn_bag3)

# b. TITN only
TITN.C.LoF.REVEL.sjlife <- table(pheno.sjlife_ttn_bag3$CMP, pheno.sjlife_ttn_bag3$TITN.Clinv.LoF.REVEL.carrierSTATUS)
TITN.C.LoF.REVEL.test.sjlife <- fisher.test(TITN.C.LoF.REVEL.sjlife) 
TITN.C.LoF.REVEL.test.sjlife.adj <- glm(formula = CMP ~  TITN.Clinv.LoF.REVEL.carrierSTATUS + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                                          PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
                                        data = pheno.sjlife_ttn_bag3)

# c. BAG3 only
BAG3.C.LoF.REVEL.sjlife <- table(pheno.sjlife_ttn_bag3$CMP, pheno.sjlife_ttn_bag3$BAG3.Clinv.LoF.REVEL.carrierSTATUS)
BAG3.C.LoF.REVEL.test.sjlife <- fisher.test(BAG3.C.LoF.REVEL.sjlife) # Got warning for chi-square
BAG3.C.LoF.REVEL.test.sjlife.adj <- glm(formula = CMP ~  BAG3.Clinv.LoF.REVEL.carrierSTATUS + agedx + agelstcontact + gender + anthra_jco_dose_any + hrtavg + PC1 + PC2 + PC3 +
                                          PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = binomial,
                                        data = pheno.sjlife_ttn_bag3)

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
#           Unadj.CCSS  adj.CCSS Unadj.SJLIFE adj.SJLIFE
# TTN+BAG3  0.4826113 0.3087160    0.9860529  0.6778491
# TTN       0.3353264 0.2975230    0.9860529  0.6778491
# BAG3      1.0000000 0.9844187    1.0000000         NA


# Clinvar, LoF and REVEL
cbind.data.frame(setNames(rbind.data.frame(c(TITN.BAG3.C.LoF.REVEL.test$p.value, coef(summary(TITN.BAG3.C.LoF.REVEL.test.adj))[2,4]),
                                           c(TITN.C.LoF.REVEL.test$p.value, coef(summary(TITN.C.LoF.REVEL.test.adj))[2,4]),
                                           c(BAG3.C.LoF.REVEL.test$p.value, coef(summary(BAG3.C.LoF.REVEL.test.adj))[2,4])
), c("Unadj.CCSS", "adj.CCSS")), 
setNames(rbind.data.frame(c(TITN.BAG3.C.LoF.REVEL.test.sjlife$p.value, coef(summary(TITN.BAG3.C.LoF.REVEL.test.sjlife.adj))[2,4]),
                          c(TITN.C.LoF.REVEL.test.sjlife$p.value, coef(summary(TITN.C.LoF.REVEL.test.sjlife.adj))[2,4]),
                          c(BAG3.C.LoF.REVEL.test.sjlife$p.value,  coef(summary(BAG3.C.LoF.REVEL.test.sjlife.adj))[2,4])
), c("Unadj.SJLIFE", "adj.SJLIFE")))
#          Unadj.CCSS  adj.CCSS Unadj.SJLIFE adj.SJLIFE
# TTN+BAG3  0.6919650 0.6507032    0.9088692  0.8266780
# TTN       0.7672838 0.7591472    1.0000000  0.9184862
# BAG3      1.0000000 0.9868182    0.6050729  0.9804864

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

SJLIFE.Clin.LoF.df <- cbind.data.frame(setNames(rbind.data.frame(c("TTN_BAG3", TITN.BAG3.C.LoF[1,2]/sum(TITN.BAG3.C.LoF[1,]), TITN.BAG3.C.LoF[2,2]/sum(TITN.BAG3.C.LoF[2,]), paste0(round(exp(c(coef(summary(TITN.BAG3.C.LoF.test.sjlife.adj))[2,1])),2), CI1), coef(summary(TITN.BAG3.C.LoF.test.sjlife.adj))[2,4]),
                                                               c("TTN", TITN.C.LoF[1,2]/sum(TITN.C.LoF[1,]), TITN.C.LoF[2,2]/sum(TITN.C.LoF[2,]), paste0(round(exp(c(coef(summary(TITN.C.LoF.test.sjlife.adj))[2,1])),2), CI2), coef(summary(TITN.C.LoF.test.sjlife.adj))[2,4]),
                                                               c("BAG3", 0, 0, NA, NA)
), c("Gene", "prevalence.in.SJLIFE.exp.CO", "prevalence.in.SJLIFE.exp.CA","OR.SJLIFE (CI)", "adj.P.SJLIFE")))

write.table(SJLIFE.Clin.LoF.df, "SJLIFE.Clin.LoF.df.txt", sep = "\t", col.names = T, row.names = F, quote = F)

## With Clinvar, LoF, and REVEL
CI1 <- paste0(" (",round(exp(c(coef(summary(TITN.BAG3.C.LoF.REVEL.test.adj))[2,1]))-1.96*as.data.frame(summary(TITN.BAG3.C.LoF.REVEL.test.adj)$coefficients)[c(1,2)][2,2],2), " to ", round(exp(c(coef(summary(TITN.BAG3.C.LoF.REVEL.test.adj))[2,1]))+1.96*as.data.frame(summary(TITN.BAG3.C.LoF.REVEL.test.adj)$coefficients)[c(1,2)][2,2], 2),")")
CI2 <- paste0(" (",round(exp(c(coef(summary(TITN.C.LoF.REVEL.test.adj))[2,1]))-1.96*as.data.frame(summary(TITN.C.LoF.REVEL.test.adj)$coefficients)[c(1,2)][2,2],2), " to ", round(exp(c(coef(summary(TITN.C.LoF.REVEL.test.adj))[2,1]))+1.96*as.data.frame(summary(TITN.C.LoF.REVEL.test.adj)$coefficients)[c(1,2)][2,2], 2),")")
CI3 <- paste0(" (",round(exp(c(coef(summary(BAG3.C.LoF.REVEL.test.adj))[2,1]))-1.96*as.data.frame(summary(BAG3.C.LoF.REVEL.test.adj)$coefficients)[c(1,2)][2,2],2), " to ", round(exp(c(coef(summary(BAG3.C.LoF.REVEL.test.adj))[2,1]))+1.96*as.data.frame(summary(BAG3.C.LoF.REVEL.test.adj)$coefficients)[c(1,2)][2,2], 2),")")

CCSS.Clin.LoF.REVEL.df <- cbind.data.frame(setNames(rbind.data.frame(c("TTN_BAG3", TITN.BAG3.C.LoF.REVEL[1,2]/sum(TITN.BAG3.C.LoF.REVEL[1,]), TITN.BAG3.C.LoF.REVEL[2,2]/sum(TITN.BAG3.C.LoF.REVEL[2,]), paste0(round(exp(c(coef(summary(TITN.BAG3.C.LoF.REVEL.test.adj))[2,1])),2), CI1), coef(summary(TITN.BAG3.C.LoF.REVEL.test.adj))[2,4]),
                                                               c("TTN", TITN.C.LoF.REVEL[1,2]/sum(TITN.C.LoF.REVEL[1,]), TITN.C.LoF.REVEL[2,2]/sum(TITN.C.LoF.REVEL[2,]), paste0(round(exp(c(coef(summary(TITN.C.LoF.REVEL.test.adj))[2,1])),2), CI2), coef(summary(TITN.C.LoF.REVEL.test.adj))[2,4]),
                                                               c("BAG3", 0, 0, paste0(round(exp(c(coef(summary(BAG3.C.LoF.REVEL.test.adj))[2,1])),2), CI3), coef(summary(BAG3.C.LoF.REVEL.test.adj))[2,4])
), c("Gene", "prevalence.in.CCSS.exp.CO", "prevalence.in.CCSS.exp.CA","OR.CCSS (CI)", "adj.P.CCSS")))

write.table(CCSS.Clin.LoF.REVEL.df, "CCSS.Clin.LoF.REVEL.df.txt", sep = "\t", col.names = T, row.names = F, quote = F)

CI1 <- paste0(" (",round(exp(c(coef(summary(TITN.BAG3.C.LoF.REVEL.test.sjlife.adj))[2,1]))-1.96*as.data.frame(summary(TITN.BAG3.C.LoF.REVEL.test.sjlife.adj)$coefficients)[c(1,2)][2,2],2), " to ", round(exp(c(coef(summary(TITN.BAG3.C.LoF.REVEL.test.sjlife.adj))[2,1]))+1.96*as.data.frame(summary(TITN.BAG3.C.LoF.REVEL.test.sjlife.adj)$coefficients)[c(1,2)][2,2], 2),")")
CI2 <- paste0(" (",round(exp(c(coef(summary(TITN.C.LoF.REVEL.test.sjlife.adj))[2,1]))-1.96*as.data.frame(summary(TITN.C.LoF.REVEL.test.sjlife.adj)$coefficients)[c(1,2)][2,2],2), " to ", round(exp(c(coef(summary(TITN.C.LoF.REVEL.test.sjlife.adj))[2,1]))+1.96*as.data.frame(summary(TITN.C.LoF.REVEL.test.sjlife.adj)$coefficients)[c(1,2)][2,2], 2),")")
CI3 <- paste0(" (",round(exp(c(coef(summary(BAG3.C.LoF.REVEL.test.sjlife.adj))[2,1]))-1.96*as.data.frame(summary(BAG3.C.LoF.REVEL.test.sjlife.adj)$coefficients)[c(1,2)][2,2],2), " to ", round(exp(c(coef(summary(BAG3.C.LoF.REVEL.test.sjlife.adj))[2,1]))+1.96*as.data.frame(summary(BAG3.C.LoF.REVEL.test.sjlife.adj)$coefficients)[c(1,2)][2,2], 2),")")

SJLIFE.Clin.LoF.REVEL.df <- cbind.data.frame(setNames(rbind.data.frame(c("TTN_BAG3", TITN.BAG3.C.LoF.REVEL.sjlife[1,2]/sum(TITN.BAG3.C.LoF.REVEL.sjlife[1,]), TITN.BAG3.C.LoF.REVEL.sjlife[2,2]/sum(TITN.BAG3.C.LoF.REVEL.sjlife[2,]), paste0(round(exp(c(coef(summary(TITN.BAG3.C.LoF.REVEL.test.sjlife.adj))[2,1])),2), CI1), coef(summary(TITN.BAG3.C.LoF.REVEL.test.sjlife.adj))[2,4]),
                                                               c("TTN", TITN.C.LoF.REVEL.sjlife[1,2]/sum(TITN.C.LoF.REVEL.sjlife[1,]), TITN.C.LoF.REVEL.sjlife[2,2]/sum(TITN.C.LoF.REVEL.sjlife[2,]), paste0(round(exp(c(coef(summary(TITN.C.LoF.REVEL.test.sjlife.adj))[2,1])),2), CI2), coef(summary(TITN.C.LoF.REVEL.test.sjlife.adj))[2,4]),
                                                               c("BAG3", BAG3.C.LoF.REVEL.sjlife[1,2]/sum(BAG3.C.LoF.REVEL.sjlife[1,]), BAG3.C.LoF.REVEL.sjlife[2,2]/sum(BAG3.C.LoF.REVEL.sjlife[2,]), paste0(round(exp(c(coef(summary(BAG3.C.LoF.REVEL.test.sjlife.adj))[2,1])),2), CI3), coef(summary(BAG3.C.LoF.REVEL.test.sjlife.adj))[2,4])
), c("Gene", "prevalence.in.SJLIFE.exp.CO", "prevalence.in.SJLIFE.exp.CA","OR.SJLIFE (CI)", "adj.P.SJLIFE")))

write.table(SJLIFE.Clin.LoF.REVEL.df, "SJLIFE.Clin.LoF.REVEL.df.txt", sep = "\t", col.names = T, row.names = F, quote = F)


save.image("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/analysis_model_fitting/2.Model_fitting_v2.RDATA")
