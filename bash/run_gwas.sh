#!/bin/bash

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3

## Run GWAS for each cohort
module load plink/1.90b
# SJLIFE
 plink \
 --allow-no-sex \
 --bfile sjlife \
 --maf 0.01 \
 --hwe 1e-06 \
 --logistic hide-covar \
 --ci 0.95 \
 --pheno pheno/sjlife_ttn_bag3.pheno \
 --pheno-name CMP \
 --covar pheno/sjlife_ttn_bag3.pheno \
 --covar-name agedx,agelstcontact,gender,anthra_jco_dose_any,hrtavg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 --out sjlife_results
# CCSS org
 plink \
 --allow-no-sex \
 --bfile ccss_org \
 --maf 0.01 \
 --hwe 1e-06 \
 --logistic hide-covar \
 --ci 0.95 \
 --pheno pheno/ccss_org_eur_cardiotoxic_exposed.pheno \
 --pheno-name CMP2plus \
 --covar pheno/ccss_org_eur_cardiotoxic_exposed.pheno \
 --covar-name a_dx,a_end,SEX,anth_DED,HeartAvg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 --out ccss_org_results
# CCSS exp
 plink \
 --allow-no-sex \
 --bfile ccss_exp \
 --maf 0.01 \
 --hwe 1e-06 \
 --logistic hide-covar \
 --ci 0.95 \
 --pheno pheno/ccss_exp_eur_cardiotoxic_exposed.pheno \
 --pheno-name CMP2plus \
 --covar pheno/ccss_exp_eur_cardiotoxic_exposed.pheno \
 --covar-name a_dx,a_end,SEX,anth_DED,HeartAvg,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
 --out ccss_exp_results

## Process the results
cat sjlife_results.assoc.logistic | grep -vE 'NA|Inf' | awk '{$1=$1; print $0}' | awk '!a[$0]++' | sort -k12,12g > sjlife_results.assoc.logistic.clean.Psorted
cat ccss_org_results.assoc.logistic | grep -vE 'NA|Inf' | awk '{$1=$1; print $0}' | awk '!a[$0]++' | sort -k12,12g > ccss_org_results.assoc.logistic.clean.Psorted
cat ccss_exp_results.assoc.logistic | grep -vE 'NA|Inf' | awk '{$1=$1; print $0}' | awk '!a[$0]++' | sort -k12,12g > ccss_exp_results.assoc.logistic.clean.Psorted

