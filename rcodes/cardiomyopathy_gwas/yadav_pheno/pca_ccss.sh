#!/bin/bash

# Calculate top 10 PCs for analytical study samples
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Cardiotoxicity/common/gwas/pheno
mkdir -p pca
mkdir -p pca/ccss

module load plink/1.90b

### CCSS original cohort
# Extract LD-pruned SNPs; use genotyped SNPs for CCSS original cohort
plink --bfile ../../../../../../yasuigrp/projects/CCSSGWAS/common/CCSS_plink_Final_Analysis_Set/CCSS.AnalysisSet.1-22 --maf 0.1 --geno 0.05 --hwe 1e-06 --indep-pairwise 100 25 0.2 --out pca/ccss/org_chr_1_22
plink --bfile ../../../../../../yasuigrp/projects/CCSSGWAS/common/CCSS_plink_Final_Analysis_Set/CCSS.AnalysisSet.1-22 --extract pca/ccss/org_chr_1_22.prune.in --make-bed --out pca/ccss/org_chr_1_22_pruned
# Lift over hg19 positions to hg38
awk '{print "chr"$1"\t"$4-1"\t"$4"\t"$1"_"$4}' pca/ccss/org_chr_1_22_pruned.bim > pca/ccss/org_chr_1_22_pruned_hg19.bed
module load liftover
liftOver pca/ccss/org_chr_1_22_pruned_hg19.bed ../../../../Genomics/common/ccss_org_hrc/hg19ToHg38.over.chain pca/ccss/org_chr_1_22_pruned_hg38.bed pca/ccss/org_chr_1_22_pruned_hg19.unmapped
# Drop unmapped varints
grep -v ^# pca/ccss/org_chr_1_22_pruned_hg19.unmapped | sed 's/^chr//g' | awk '{print $1"_"$3}' > pca/ccss/drop_variants_unmapped.txt
plink --bfile pca/ccss/org_chr_1_22_pruned --exclude pca/ccss/drop_variants_unmapped.txt --make-bed --out pca/ccss/org_chr_1_22_pruned_mapped
# Update coordinates and name
sed 's/^chr//g' pca/ccss/org_chr_1_22_pruned_hg38.bed | awk '{print $4, $3}' > pca/ccss/org_chr_1_22_pruned_update_pos.txt
sed 's/^chr//g' pca/ccss/org_chr_1_22_pruned_hg38.bed | awk '{print $4, $1":"$3}' > pca/ccss/org_chr_1_22_pruned_update_name.txt
plink --bfile pca/ccss/org_chr_1_22_pruned_mapped --update-map pca/ccss/org_chr_1_22_pruned_update_pos.txt --make-bed --out pca/ccss/org_chr_1_22_pruned_mapped_hg38
plink --bfile pca/ccss/org_chr_1_22_pruned_mapped_hg38 --update-name pca/ccss/org_chr_1_22_pruned_update_name.txt --make-bed --out pca/ccss/org_chr_1_22_pruned_mapped_hg38_varname_updated

### CCSS expansion cohort
# Extract LD-pruned SNPs
for chr in {1..22};do
 bsub -P gwas -J pca -eo pca/ccss/${chr}.err -oo pca/ccss/${chr}.out -R "rusage[mem=10000]" \
 "plink --vcf ../../../../Genomics/common/ccss_exp_wgs/CCSS.GATKv3.4.VQSR_chr${chr}.PASS.decomposed.ccssid.vcf.gz --maf 0.1 --geno 0.05 --hwe 1e-06 --indep-pairwise 100 25 0.2 --out pca/ccss/exp_chr${chr};
 plink --vcf ../../../../Genomics/common/ccss_exp_wgs/CCSS.GATKv3.4.VQSR_chr${chr}.PASS.decomposed.ccssid.vcf.gz --extract pca/ccss/exp_chr${chr}.prune.in --make-bed --out pca/ccss/exp_chr${chr}_pruned"
done
# Merge
for chr in {1..22};do
 echo pca/ccss/exp_chr${chr}_pruned.bed pca/ccss/exp_chr${chr}_pruned.bim pca/ccss/exp_chr${chr}_pruned.fam >> pca/ccss/exp_merge_list.txt
done
plink --merge-list pca/ccss/exp_merge_list.txt --make-bed --out pca/ccss/exp_chr_1_22_pruned
# Update variant names to match with those in the CCSS original cohort while excluding duplicates based on chr:pos
awk '{print $2, $1":"$4}' pca/ccss/exp_chr_1_22_pruned.bim > pca/ccss/exp_chr_1_22_pruned_update_name.txt
awk '!a[$1":"$4]++{print $2}' pca/ccss/exp_chr_1_22_pruned.bim > pca/ccss/exp_chr_1_22_pruned_keep.txt
plink --bfile pca/ccss/exp_chr_1_22_pruned --extract pca/ccss/exp_chr_1_22_pruned_keep.txt --make-bed --out pca/ccss/exp_chr_1_22_pruned_uniq
plink --bfile pca/ccss/exp_chr_1_22_pruned_uniq --update-name pca/ccss/exp_chr_1_22_pruned_update_name.txt --make-bed --out pca/ccss/exp_chr_1_22_pruned_varname_updated

## Combine data from CCSS original and expansion cohorts
plink --bfile pca/ccss/org_chr_1_22_pruned_mapped_hg38_varname_updated --bmerge pca/ccss/exp_chr_1_22_pruned_varname_updated --make-bed --out pca/ccss/org_plus_exp
# Some variants need to be flipped
plink --bfile pca/ccss/org_chr_1_22_pruned_mapped_hg38_varname_updated --flip pca/ccss/org_plus_exp-merge.missnp --make-bed --out pca/ccss/org_chr_1_22_pruned_mapped_hg38_varname_updated_flipped
# Merge again
plink --bfile pca/ccss/org_chr_1_22_pruned_mapped_hg38_varname_updated_flipped --bmerge pca/ccss/exp_chr_1_22_pruned_varname_updated --make-bed --out pca/ccss/org_plus_exp
# There are still variants with mismatching alleles between the 2 datasets; drop them
plink --bfile pca/ccss/org_chr_1_22_pruned_mapped_hg38_varname_updated_flipped --exclude pca/ccss/org_plus_exp-merge.missnp --make-bed --out pca/ccss/org_chr_1_22_pruned_mapped_hg38_varname_updated_flipped_dropped
# Merge
plink --bfile pca/ccss/org_chr_1_22_pruned_mapped_hg38_varname_updated_flipped_dropped --bmerge pca/ccss/exp_chr_1_22_pruned_varname_updated --geno 0.05 --hwe 1e-06 --maf 0.1 --make-bed --out pca/ccss/org_plus_exp

### Now compute top PCs for the replication analytical sample (europeans doxorubicin-exposed only; no exposure to chest RT and daunorubicin; motivated by results in the JNCI paper)
plink --bfile pca/ccss/org_plus_exp --keep ccss_eur_dox_only.txt --maf 0.1 --geno 0.05 --hwe 1e-06 --pca 20 'header' --out pca/ccss/org_plus_exp_eur_dox_only_top_20_pc

