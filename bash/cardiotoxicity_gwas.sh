#!/usr/bin/bash

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap
ln -s ../../../../Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr16.PASS.decomposed.vcf.gz .
PHENO=../pheno/sjlife_eur_dox_only_pcs.pheno

# extract CMP samples from $PHENO and variants from ../CMP_chr1_22.assoc.logistic.clean.Psorted_chr16_25611595_250kb

awk '{print $2}' $PHENO > samples.list

module load plink/1.90b
plink --chr 16 --from-bp 25361595 --make-bed --out sjlife_CMP --to-bp 25861595 --vcf MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr16.PASS.decomposed.vcf.gz

awk '{print $2}' ../CMP_chr1_22.assoc.logistic.clean.Psorted_chr16_25611595_250kb > CMP_SNP.list
# keep samples and variants
plink --bfile sjlife --extract samplesnp.list --keep samples_for_maf.txt --keep-allele-order --make-bed --out samplesnp.dat

