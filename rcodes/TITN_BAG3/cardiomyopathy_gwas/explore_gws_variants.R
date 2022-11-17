## Phenotype data for SJLIFE
rm(list=ls())
setwd('/Volumes/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas')
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas")
source('mafCalc.R')

## Phenotype data used in GWAS
pheno_gwas = read.table('pheno/sjlife_eur_dox_only_pcs.pheno', header = TRUE)

## Genotypes of genome-wide significant SNPs
geno_chr3 = read.table('sjlife_chr3.raw', header = TRUE)
geno_chr3[c(2:6)] = NULL
geno_chr16 = read.table('sjlife_chr16.raw', header = TRUE)
geno_chr16[c(2:6)] = NULL
geno = merge(geno_chr3, geno_chr16, by="FID")

## Merge phenotype and genotype data
dat = merge(pheno_gwas, geno, by="FID")
dat$CMP2plus = dat$CMP2plus - 1

########## Check if variants on chr16 are indepedent
results = NULL
vars = colnames(dat)[grep("^chr", colnames(dat))]
for (i in 1:length(vars)){
  varname = vars[i]
  varname_split = unlist(strsplit(varname, "\\.|_"))
  varname_split = varname_split[varname_split != ""]
  chr = varname_split[1]
  bp = varname_split[2]
  ea = varname_split[5]
  nea = varname_split[6]
  maf_cases = mafcalc(dat[[varname]][dat$CMP2plus==1])
  maf_cont = mafcalc(dat[[varname]][dat$CMP2plus==0])
  fit = glm(CMP2plus~dat[[varname]]+agedx+as.factor(gender)+agelstcontact+anthra_jco_dose_any+
                      PC1+PC2+PC3+PC4+PC5, data = dat, family = binomial)
  logodds = summary(fit)$coef[2,1]
  se = summary(fit)$coef[2,2]
  or = exp(logodds)
  l95 = exp(logodds - (1.96*se))
  u95 = exp(logodds + (1.96*se))
  p = summary(fit)$coef[2,4]
  results = rbind(results, data.frame(chr, bp, ea, nea, maf_cases, maf_cont, or, l95, u95, p))
}
print(results)

# Firth's bias-reduced logistic regression
library(logistf)
summary(logistf(CMP2plus~dat$chr16.25611595.C.T_T..C.+agedx+as.factor(gender)+agelstcontact+anthra_jco_dose_any+PC1+PC2+PC3+PC4+PC5, data = dat, family = binomial))
