setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/")
## Load echo data
load('Rcodes/EUR_common_variants.RData')
load('Rcodes/AFR_common_variants.RData')

## Need to merge echo data for SJLIFE EUR and SJLIFE AFR cohorts and we only want to analyze for:
# c("LV_Ejection_Fraction_3D", "LV_End_Diastolic_Volume_3D", "LV_End_Systolic_Volume_3D",
# "LV_Stroke_Volume_3D", "LVMassMM_Index", "LV_GLPS_AVG", "LV_Relative_Wall_Thickness")


EUR_common_variants$rs3829746_yn = ifelse(EUR_common_variants$rs3829746>0, "Yes", "No")
EUR_common_variants$rs2234962_yn = ifelse(EUR_common_variants$rs2234962>0, "Yes", "No")
AFR_common_variants$rs3829746_yn = ifelse(AFR_common_variants$rs3829746>0, "Yes", "No")
AFR_common_variants$rs2234962_yn = ifelse(AFR_common_variants$rs2234962>0, "Yes", "No")
# Clean the outliers
# Some EF values are in fractions; change them to percentages
# EUR
EUR_common_variants$LV_Ejection_Fraction_3D = ifelse(EUR_common_variants$LV_Ejection_Fraction_3D<1, EUR_common_variants$LV_Ejection_Fraction_3D*100,
                                              EUR_common_variants$LV_Ejection_Fraction_3D)
EUR_common_variants$LV_End_Diastolic_Volume_3D[EUR_common_variants$LV_End_Diastolic_Volume_3D<10] = NA
EUR_common_variants$LV_End_Systolic_Volume_3D[EUR_common_variants$LV_End_Systolic_Volume_3D<10 | EUR_common_variants$LV_End_Systolic_Volume_3D>150] = NA
EUR_common_variants$LV_Stroke_Volume_3D[EUR_common_variants$LV_Stroke_Volume_3D<10] = NA
EUR_common_variants$LVMassMM_Index[EUR_common_variants$LVMassMM_Index<10] = NA
# AFR
AFR_common_variants$LV_Ejection_Fraction_3D = ifelse(AFR_common_variants$LV_Ejection_Fraction_3D<1, AFR_common_variants$LV_Ejection_Fraction_3D*100,
                                              AFR_common_variants$LV_Ejection_Fraction_3D)
AFR_common_variants$LV_End_Diastolic_Volume_3D[AFR_common_variants$LV_End_Diastolic_Volume_3D<10] = NA
AFR_common_variants$LV_End_Systolic_Volume_3D[AFR_common_variants$LV_End_Systolic_Volume_3D<10 | AFR_common_variants$LV_End_Systolic_Volume_3D>150] = NA
AFR_common_variants$LV_Stroke_Volume_3D[AFR_common_variants$LV_Stroke_Volume_3D<10] = NA
AFR_common_variants$LVMassMM_Index[AFR_common_variants$LVMassMM_Index<10] = NA

## Visualize the data
library(ggplot2)
library(rstatix)
EUR_common_variants = subset(EUR_common_variants, !is.na(rs2234962_yn) & !is.na(rs3829746_yn))
# Only keep SJLIFE visits 1, 2 and 3 due to numbers
EUR_common_variants = subset(EUR_common_variants, visittype=="SJLIFE Visit 1" | visittype=="SJLIFE Visit 2")
# pdf('echo_parameters_BAG3_genotype_SJLIFE_eur.pdf')
ggplot(EUR_common_variants, aes(rs2234962_yn, LV_Ejection_Fraction_3D, fill = visittype)) + geom_boxplot()
ggplot(EUR_common_variants, aes(rs2234962_yn, LV_End_Diastolic_Volume_3D, fill = visittype)) + geom_boxplot()
ggplot(EUR_common_variants, aes(rs2234962_yn, LV_End_Systolic_Volume_3D, fill = visittype)) + geom_boxplot()
ggplot(EUR_common_variants, aes(rs2234962_yn, LV_GLPS_AVG, fill = visittype)) + geom_boxplot()
# dev.off()
# Summary statistics
group_by(EUR_common_variants, rs2234962_yn, visittype) %>% get_summary_stats(LV_Ejection_Fraction_3D, show = c("mean", "median", "sd"))
group_by(EUR_common_variants, rs2234962_yn, visittype) %>% get_summary_stats(LV_End_Diastolic_Volume_3D, show = c("mean", "median", "sd"))
group_by(EUR_common_variants, rs2234962_yn, visittype) %>% get_summary_stats(LV_End_Systolic_Volume_3D, show = c("mean", "median", "sd"))
group_by(EUR_common_variants, rs2234962_yn, visittype) %>% get_summary_stats(LV_GLPS_AVG, show = c("mean", "median", "sd"))

## Baseline
eur_baseline = subset(EUR_common_variants, visittype=="SJLIFE Visit 1")
dim(eur_baseline)
afr_baseline = subset(AFR_common_variants, visittype=="SJLIFE Visit 1")
dim(afr_baseline)

## Fit regression models for each visit (EUR)
param = colnames(eur_baseline)[grep("^LV", colnames(eur_baseline))]
snps = c('rs3829746', 'rs2234962')
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
snps = c('rs3829746', 'rs2234962')
dat = afr_baseline
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

##########
## Male ##
##########
## Baseline
eur_baseline = subset(EUR_common_variants, visittype=="SJLIFE Visit 1")
## 0 is  male, 1 is female
eur_baseline <- eur_baseline[eur_baseline$gender == 0,]
dim(eur_baseline)
afr_baseline = subset(AFR_common_variants, visittype=="SJLIFE Visit 1")
## 0 is  male, 1 is female
afr_baseline <- afr_baseline[afr_baseline$gender == 0,]
dim(afr_baseline)

## Fit regression models for each visit (EUR)
param = colnames(eur_baseline)[grep("^LV", colnames(eur_baseline))]
snps = c('rs3829746', 'rs2234962')
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
snps = c('rs3829746', 'rs2234962')
dat = afr_baseline
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



############
## Female ##
############
## Baseline
eur_baseline = subset(EUR_common_variants, visittype=="SJLIFE Visit 1")
## 0 is  male, 1 is female
eur_baseline <- eur_baseline[eur_baseline$gender == 1,]
dim(eur_baseline)
afr_baseline = subset(AFR_common_variants, visittype=="SJLIFE Visit 1")
## 0 is  male, 1 is female
afr_baseline <- afr_baseline[afr_baseline$gender == 1,]
dim(afr_baseline)

## Fit regression models for each visit (EUR)
param = colnames(eur_baseline)[grep("^LV", colnames(eur_baseline))]
snps = c('rs3829746', 'rs2234962')
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
snps = c('rs3829746', 'rs2234962')
dat = afr_baseline
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


# # Regression for change in echo parameters
# library(dplyr)
# p = 'LV_Ejection_Fraction_3D'
# newdat = EUR_common_variants %>% group_by(sjlid) %>% mutate(diff=LV_Ejection_Fraction_3D-lag(LV_Ejection_Fraction_3D, default=first(LV_Ejection_Fraction_3D)))
