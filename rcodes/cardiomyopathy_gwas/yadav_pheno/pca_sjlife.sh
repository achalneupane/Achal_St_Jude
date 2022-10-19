#!/bin/bash

# Calculate top 10 PCs for analytical study samples
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno
mkdir -p pca
mkdir -p pca/sjlife

module load plink/1.90b

# Extract LD-pruned SNPs
for chr in {1..22};do
 bsub -P gwas -J pca -eo pca/sjlife/${chr}.err -oo pca/sjlife/${chr}.out -R "rusage[mem=10000]" \
 "plink --vcf ../../../../Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${chr}.PASS.decomposed.vcf.gz --maf 0.1 --geno 0.05 --hwe 1e-06 --indep-pairwise 100 25 0.2 --out pca/sjlife/chr${chr};
 plink --vcf ../../../../Genomics/common/sjlife/MERGED_SJLIFE_1_2/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${chr}.PASS.decomposed.vcf.gz --extract pca/sjlife/chr${chr}.prune.in --make-bed --out pca/sjlife/chr${chr}_pruned"
done

# Merge
for chr in {1..22};do
 echo pca/sjlife/chr${chr}_pruned.bed pca/sjlife/chr${chr}_pruned.bim pca/sjlife/chr${chr}_pruned.fam >> pca/sjlife/merge_list.txt
done
plink --merge-list pca/sjlife/merge_list.txt --make-bed --out pca/sjlife/chr1_22_pruned

# Calculate top 20 PCs for the GWAS analytical sample (europeans doxorubicin-exposed only; no exposure to chest RT and daunorubicin; motivated by results in the JNCI paper)
plink --bfile pca/sjlife/chr1_22_pruned --keep sjlife_eur_dox_only.txt --maf 0.1 --geno 0.05 --hwe 1e-06 --pca 20 'header' --out pca/sjlife/chr1_22_pruned_eur_dox_only_top_20_pc

