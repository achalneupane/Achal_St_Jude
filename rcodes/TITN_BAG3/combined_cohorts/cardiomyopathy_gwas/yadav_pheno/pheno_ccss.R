## Phenotype data for CCSS (to replicate GWAS findings in SJLIFE)
rm(list=ls())
# setwd('/Volumes/Groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno/')
setwd("Z:/ResearchHome/ClusterHome/aneupane/St_Jude/Achal_St_Jude/rcodes/TITN_BAG3/cardiomyopathy_gwas/yadav_pheno/")
## Phenotype data by Huiqi and Qi
load('ccss_export_05062022.RData')
pheno = export
# format the phenotype
pheno$SEX = pheno$SEX - 1
pheno$CMP2plus = ifelse(pheno$maxCHF15>1, 1, ifelse(pheno$maxCHF15==0, 0, NA))
pheno$CMP3plus = ifelse(pheno$maxCHF15>2, 1, ifelse(pheno$maxCHF15==0, 0, NA))

## Exclude overlapping SJLIFE samples in GWAS
sjlife_gwas = read.table('sjlife_eur_dox_only.txt', header = TRUE)
sjlife_overlap_ccss = read.table('SJLIFE_WGS_samples_overlap_with_CCSS_org_SNP_samples.txt', header = TRUE)
sjlife_overlap_ccss_gwas = subset(sjlife_overlap_ccss, V1 %in% sjlife_gwas$FID)
pheno_uniq = subset(pheno, !(ccssid %in% sjlife_overlap_ccss_gwas$ccssid))

## European ancestry
# Original
ccss_org = read.table('CCSS.SJLIFE.ancestry_CCSS.txt', header = TRUE, sep="\t")
ccss_org_eur = subset(ccss_org, CEU>=0.8)
# Expansion
ccss_exp = read.table('CCSS_and_1kGP_final_EUR_AFR_EAS.3.Q_ccssexp_samples', header = TRUE)
ccss_exp_eur = subset(ccss_exp, EUR>=0.8)
# CCSS eur
ccss_eur = c(ccss_org_eur$SAMPLE, ccss_exp_eur$ccssid)

## Phenotype data in Europeans only
final = subset(pheno_uniq, ccssid %in% ccss_eur)
# Save everything for future re-analyses
write.table(final, 'ccss_eur_all.txt', row.names=F, quote=F, sep="\t")

## Based on previous JNCI analyses, genetic effects were stronger among survivors exposed to doxorubicin only (no chest RT and daunorubicin)
final_dox = subset(final, doxo_totdosepersqm>0 & daun_totdosepersqm==0 & chestrt_yn==2)

## Subset including required variables only and format to PLINK to allow addition of top PCs only
plinkfile = final_dox[c('ccssid', 'ccssid', 'CMP2plus', 'a_dx','SEX','a_end', 'anth_DED')]
plinkfile = plinkfile[complete.cases(plinkfile),]
# Add CMP3plus variable
cmp3plus = pheno[c('ccssid', 'CMP3plus')]
plinkfile = merge(plinkfile, cmp3plus, by='ccssid')
colnames(plinkfile)[c(1,2)] = c('FID', 'IID')
write.table(plinkfile, 'ccss_eur_dox_only.txt', row.names = FALSE, quote = FALSE)
