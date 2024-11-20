setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/")
## Load echo data
load('Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/DNUT/echo.PLP.eur.RData')
load('Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/Rcodes/DNUT/echo.PLP.afr.RData')
echo.PLP.eur$chr2.178562809.T.C_C_yn = ifelse(echo.PLP.eur$chr2.178562809.T.C_C>0, "Yes", "No")
echo.PLP.eur$chr10.119670121.T.C_C_yn = ifelse(echo.PLP.eur$chr10.119670121.T.C_C>0, "Yes", "No")
echo.PLP.afr$chr2.178562809.T.C_C_yn = ifelse(echo.PLP.afr$chr2.178562809.T.C_C>0, "Yes", "No")
echo.PLP.afr$chr10.119670121.T.C_C_yn = ifelse(echo.PLP.afr$chr10.119670121.T.C_C>0, "Yes", "No")
# Clean the outliers
# Some EF values are in fractions; change them to percentages
# EUR
echo.PLP.eur$LV_Ejection_Fraction_3D = ifelse(echo.PLP.eur$LV_Ejection_Fraction_3D<1, echo.PLP.eur$LV_Ejection_Fraction_3D*100,
                                              echo.PLP.eur$LV_Ejection_Fraction_3D)
echo.PLP.eur$LV_End_Diastolic_Volume_3D[echo.PLP.eur$LV_End_Diastolic_Volume_3D<10] = NA
echo.PLP.eur$LV_End_Systolic_Volume_3D[echo.PLP.eur$LV_End_Systolic_Volume_3D<10 | echo.PLP.eur$LV_End_Systolic_Volume_3D>150] = NA
echo.PLP.eur$LV_Stroke_Volume_3D[echo.PLP.eur$LV_Stroke_Volume_3D<10] = NA
echo.PLP.eur$LVMassMM_Index[echo.PLP.eur$LVMassMM_Index<10] = NA
# AFR
echo.PLP.afr$LV_Ejection_Fraction_3D = ifelse(echo.PLP.afr$LV_Ejection_Fraction_3D<1, echo.PLP.afr$LV_Ejection_Fraction_3D*100,
                                              echo.PLP.afr$LV_Ejection_Fraction_3D)
echo.PLP.afr$LV_End_Diastolic_Volume_3D[echo.PLP.afr$LV_End_Diastolic_Volume_3D<10] = NA
echo.PLP.afr$LV_End_Systolic_Volume_3D[echo.PLP.afr$LV_End_Systolic_Volume_3D<10 | echo.PLP.afr$LV_End_Systolic_Volume_3D>150] = NA
echo.PLP.afr$LV_Stroke_Volume_3D[echo.PLP.afr$LV_Stroke_Volume_3D<10] = NA
echo.PLP.afr$LVMassMM_Index[echo.PLP.afr$LVMassMM_Index<10] = NA

## Visualize the data
library(ggplot2)
library(rstatix)
echo.PLP.eur = subset(echo.PLP.eur, !is.na(chr10.119670121.T.C_C_yn) & !is.na(chr2.178562809.T.C_C_yn))
# Only keep SJLIFE visits 1, 2 and 3 due to numbers
echo.PLP.eur = subset(echo.PLP.eur, visittype=="SJLIFE Visit 1" | visittype=="SJLIFE Visit 2")
# pdf('echo_parameters_BAG3_genotype_SJLIFE_eur.pdf')
ggplot(echo.PLP.eur, aes(chr10.119670121.T.C_C_yn, LV_Ejection_Fraction_3D, fill = visittype)) + geom_boxplot()
ggplot(echo.PLP.eur, aes(chr10.119670121.T.C_C_yn, LV_End_Diastolic_Volume_3D, fill = visittype)) + geom_boxplot()
ggplot(echo.PLP.eur, aes(chr10.119670121.T.C_C_yn, LV_End_Systolic_Volume_3D, fill = visittype)) + geom_boxplot()
ggplot(echo.PLP.eur, aes(chr10.119670121.T.C_C_yn, LV_GLPS_AVG, fill = visittype)) + geom_boxplot()
# dev.off()
# Summary statistics
group_by(echo.PLP.eur, chr10.119670121.T.C_C_yn, visittype) %>% get_summary_stats(LV_Ejection_Fraction_3D, show = c("mean", "median", "sd"))
group_by(echo.PLP.eur, chr10.119670121.T.C_C_yn, visittype) %>% get_summary_stats(LV_End_Diastolic_Volume_3D, show = c("mean", "median", "sd"))
group_by(echo.PLP.eur, chr10.119670121.T.C_C_yn, visittype) %>% get_summary_stats(LV_End_Systolic_Volume_3D, show = c("mean", "median", "sd"))
group_by(echo.PLP.eur, chr10.119670121.T.C_C_yn, visittype) %>% get_summary_stats(LV_GLPS_AVG, show = c("mean", "median", "sd"))

## add diagnosis era
diagDT <- readRDS("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/pheno/diagDT.rds")

## Baseline
eur_baseline = subset(echo.PLP.eur, visittype=="SJLIFE Visit 1")
dim(eur_baseline)
eur_baseline$era_numeric <- factor(diagDT$era_numeric[match(eur_baseline$sjlid, diagDT$IID)])

afr_baseline = subset(echo.PLP.afr, visittype=="SJLIFE Visit 1")
dim(afr_baseline)
afr_baseline$era_numeric <- factor(diagDT$era_numeric[match(afr_baseline$sjlid, diagDT$IID)])

## Fit regression models for each visit (EUR)
param = colnames(eur_baseline)[grep("^LV", colnames(eur_baseline))]
snps = c('chr2.178562809.T.C_C', 'chr10.119670121.T.C_C')
dat = eur_baseline
res = NULL
for (p in 1:length(param)){
  echo_p  = param[p]
  for (s in (1:2)){
    snp_s = snps[s]
    fit = lm(scale(dat[[echo_p]])~dat[[snp_s]]+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=dat)
    beta = summary(fit)$coef[2,1]
    se = summary(fit)$coef[2,2]
    pval = summary(fit)$coef[2,4]
    res = rbind(res, data.frame(parameter = echo_p, snp = snp_s, beta, se, pval))
  }
}
print(res)

res.ttn <- res[grepl("chr2", res$snp),]
res.bag3 <- res[grepl("chr10", res$snp),]



## Fit regression models for each visit (AFR)
param = colnames(afr_baseline)[grep("^LV", colnames(afr_baseline))]
snps = c('chr2.178562809.T.C_C', 'chr10.119670121.T.C_C')
dat = afr_baseline
res = NULL
for (p in 1:length(param)){
  echo_p  = param[p]
  for (s in (1:2)){
    snp_s = snps[s]
    fit = lm(scale(dat[[echo_p]])~dat[[snp_s]]+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 + era_numeric, data=dat)
    beta = summary(fit)$coef[2,1]
    se = summary(fit)$coef[2,2]
    pval = summary(fit)$coef[2,4]
    res = rbind(res, data.frame(parameter = echo_p, snp = snp_s, beta, se, pval))
  }
}
print(res)

res.ttn <- res[grepl("chr2", res$snp),]
res.bag3 <- res[grepl("chr10", res$snp),]

##########
## Male ##
##########
## Baseline
eur_baseline = subset(echo.PLP.eur, visittype=="SJLIFE Visit 1")
## 0 is  male, 1 is female
eur_baseline <- eur_baseline[eur_baseline$gender == 0,]
dim(eur_baseline)
afr_baseline = subset(echo.PLP.afr, visittype=="SJLIFE Visit 1")
## 0 is  male, 1 is female
afr_baseline <- afr_baseline[afr_baseline$gender == 0,]
dim(afr_baseline)

## Fit regression models for each visit (EUR)
param = colnames(eur_baseline)[grep("^LV", colnames(eur_baseline))]
snps = c('chr2.178562809.T.C_C', 'chr10.119670121.T.C_C')
dat = eur_baseline
res = NULL
for (p in 1:length(param)){
  echo_p  = param[p]
  for (s in (1:2)){
    snp_s = snps[s]
    fit = lm(scale(dat[[echo_p]])~dat[[snp_s]]+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 + era_numeric, data=dat)
    beta = summary(fit)$coef[2,1]
    se = summary(fit)$coef[2,2]
    pval = summary(fit)$coef[2,4]
    res = rbind(res, data.frame(parameter = echo_p, snp = snp_s, beta, se, pval))
  }
}
print(res)

res.ttn <- res[grepl("chr2", res$snp),]
res.bag3 <- res[grepl("chr10", res$snp),]



## Fit regression models for each visit (AFR)
param = colnames(afr_baseline)[grep("^LV", colnames(afr_baseline))]
snps = c('chr2.178562809.T.C_C', 'chr10.119670121.T.C_C')
dat = afr_baseline
res = NULL
for (p in 1:length(param)){
  echo_p  = param[p]
  for (s in (1:2)){
    snp_s = snps[s]
    fit = lm(scale(dat[[echo_p]])~dat[[snp_s]]+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 + era_numeric, data=dat)
    beta = summary(fit)$coef[2,1]
    se = summary(fit)$coef[2,2]
    pval = summary(fit)$coef[2,4]
    res = rbind(res, data.frame(parameter = echo_p, snp = snp_s, beta, se, pval))
  }
}
print(res)

res.ttn <- res[grepl("chr2", res$snp),]
res.bag3 <- res[grepl("chr10", res$snp),]



##########
## Female ##
##########
## Baseline
eur_baseline = subset(echo.PLP.eur, visittype=="SJLIFE Visit 1")
## 0 is  male, 1 is female
eur_baseline <- eur_baseline[eur_baseline$gender == 1,]
dim(eur_baseline)
afr_baseline = subset(echo.PLP.afr, visittype=="SJLIFE Visit 1")
## 0 is  male, 1 is female
afr_baseline <- afr_baseline[afr_baseline$gender == 1,]
dim(afr_baseline)

## Fit regression models for each visit (EUR)
param = colnames(eur_baseline)[grep("^LV", colnames(eur_baseline))]
snps = c('chr2.178562809.T.C_C', 'chr10.119670121.T.C_C')
dat = eur_baseline
res = NULL
for (p in 1:length(param)){
  echo_p  = param[p]
  for (s in (1:2)){
    snp_s = snps[s]
    fit = lm(scale(dat[[echo_p]])~dat[[snp_s]]+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 + era_numeric, data=dat)
    beta = summary(fit)$coef[2,1]
    se = summary(fit)$coef[2,2]
    pval = summary(fit)$coef[2,4]
    res = rbind(res, data.frame(parameter = echo_p, snp = snp_s, beta, se, pval))
  }
}
print(res)

res.ttn <- res[grepl("chr2", res$snp),]
res.bag3 <- res[grepl("chr10", res$snp),]



## Fit regression models for each visit (AFR)
param = colnames(afr_baseline)[grep("^LV", colnames(afr_baseline))]
snps = c('chr2.178562809.T.C_C', 'chr10.119670121.T.C_C')
dat = afr_baseline
res = NULL
for (p in 1:length(param)){
  echo_p  = param[p]
  for (s in (1:2)){
    snp_s = snps[s]
    fit = lm(scale(dat[[echo_p]])~dat[[snp_s]]+agedx+agelstcontact+gender+anthra_jco_dose_any+hrtavg+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10 + era_numeric, data=dat)
    beta = summary(fit)$coef[2,1]
    se = summary(fit)$coef[2,2]
    pval = summary(fit)$coef[2,4]
    res = rbind(res, data.frame(parameter = echo_p, snp = snp_s, beta, se, pval))
  }
}
print(res)

res.ttn <- res[grepl("chr2", res$snp),]
res.bag3 <- res[grepl("chr10", res$snp),]


# # Regression for change in echo parameters
# library(dplyr)
# p = 'LV_Ejection_Fraction_3D'
# newdat = echo.PLP.eur %>% group_by(sjlid) %>% mutate(diff=LV_Ejection_Fraction_3D-lag(LV_Ejection_Fraction_3D, default=first(LV_Ejection_Fraction_3D)))
