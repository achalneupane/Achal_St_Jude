#!/usr/bin/bash

cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/finemap
ln -s ../../../../Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr16.PASS.decomposed.vcf.gz .
PHENO=../pheno/sjlife_eur_dox_only_pcs.pheno

# extract CMP samples from $PHENO and variants from ../CMP_chr1_22.assoc.logistic.clean.Psorted_chr16_25611595_250kb

awk '{print $1"\t"$2}' $PHENO > samples.list

module load plink/1.90b
plink --chr 16 --from-bp 25361595 --make-bed --out sjlife_CMP --to-bp 25861595 --vcf MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr16.PASS.decomposed.vcf.gz

plink --bfile sjlife_CMP --freq --keep-allele-order --out sjlife.freq.out
awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' sjlife.freq.out.frq > sjlife.freq.out.frq_edited1

awk '{print $2}' ../CMP_chr1_22.assoc.logistic.clean.Psorted_chr16_25611595_250kb > samplesnp.list
# keep samples and variants
plink --bfile sjlife_CMP --extract samplesnp.list --keep samples.list --keep-allele-order --make-bed --out samplesnp.dat


## Create LD matrix
plink --r2 --bfile samplesnp.dat --matrix --out samplesnp_ld_matrix


# # chmod 777 *
# # dos2unix samplesnp.z

## Without any assumption for causal SNPs
/research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/ttn_bag3/FINEMAP/finemap_v1.4.1_x86_64/finemap_v1.4.1_x86_64 --sss --log --in-files samplesnp --dataset 1