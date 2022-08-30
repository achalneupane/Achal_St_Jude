#!/usr/bin/bash

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemapc
ln -s ../../../../Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr16.PASS.decomposed.vcf.gz .
PHENO=../pheno/sjlife_eur_dox_only_pcs.pheno

# extract CMP samples from $PHENO and variants from ../CMP_chr1_22.assoc.logistic.clean.Psorted_chr16_25611595_250kb
awk '{print $2}' $PHENO > samples.list
awk '{print $2}' ../CMP_chr1_22.assoc.logistic.clean.Psorted_chr16_25611595_250kb > CMP_SNP.list
samplesnp_TITN_gt_MAF_1_perc_vars.list
