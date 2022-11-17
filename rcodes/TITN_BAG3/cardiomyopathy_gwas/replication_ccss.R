rm(list=ls())
# setwd('/Volumes/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas')
setwd("Z:/ResearchHome/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas")
source('mafCalc.R')

## Phenotype data in CCSS for replication
pheno = read.table('pheno/ccss_eur_dox_only.txt', header = TRUE)

## Add top 5 PCs
pcs = read.table('pheno/pca/ccss/org_plus_exp_eur_dox_only_top_20_pc.eigenvec', header = TRUE)
pcs_5 = pcs[c(1,3:7)]

## Merge phenotype with top 5 pcs
pheno_pcs = merge(pheno, pcs_5, by = 'FID')

## Genotype data
# Org
org_3 = read.table('ccss_org_chr3.raw', header = TRUE)
org_3[c(2:6)] = NULL
org_16 = read.table('ccss_org_chr16.raw', header = TRUE)
org_16[c(2:6)] = NULL
org_3_16 = merge(org_3, org_16, by = "FID")
# Exp
exp_3 = read.table('ccss_exp_chr3.raw', header = TRUE)
exp_3[c(2:6)] = NULL
exp_16 = read.table('ccss_exp_chr16.raw', header = TRUE)
exp_16[c(2:6)] = NULL
exp_3_16 = merge(exp_3, exp_16, by = "FID")
# There are 5 variants in both org and exp data and they are in the same order; so replace column names of org data with those of exp, to allow combine
colnames(org_3_16) = colnames(exp_3_16)
geno = rbind(org_3_16, exp_3_16)

## Merge with phenotype data
dat = merge(pheno_pcs, geno, by="FID")

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
  fit = glm(CMP3plus~dat[[varname]]+a_dx+as.factor(SEX)+a_end+anth_DED+
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
summary(logistf(CMP2plus~dat$chr16.25611595.C.T_T..C.+a_dx+as.factor(SEX)+a_end+anth_DED+PC1+PC2+PC3+PC4+PC5, data = dat, family = binomial))
