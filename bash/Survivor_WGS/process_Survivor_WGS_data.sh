#!/bin/bash

module load vcftools/0.1.13
module load bcftools/1.17
#module load vt/2016.11.07
module load plink/1.90b
# module load R/3.4.0
module load R/4.2.2-rhel8
module load bedtools/2.25.0
module load htslib/1.3.1

# cd /home/ysapkota/Work/WGS_SJLIFE
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS/chr*.Survivor_WGS.GATK4180.hg38.vcf.gz .
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS/chr*.Survivor_WGS.GATK4180.hg38.vcf.gz.tbi .
cp /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS/wgs_sampleIDs_SJLIDmatch_23March2023.txt ./Phenotypes

########################################
## 1. Merge VCF for missingness check ##
########################################
module load gatk/4.1.8.0 
#!/bin/bash
# Create a list of VCF files to merge (same samples but split by chromosomes)
vcf_files=""
for chr in {1..22} X Y; do
    vcf_files+="-I chr${chr}.Survivor_WGS.GATK4180.hg38.vcf.gz "
done

# Merge VCF files using GATK GatherVCF to create a giant VCF with all chromosomes
gatk GatherVcfs $vcf_files -O chrALL.Survivor_WGS.GATK4180.hg38.vcf.gz
# create index for merged VCF
gatk IndexFeatureFile -I chrALL.Survivor_WGS.GATK4180.hg38.vcf.gz



## Get VCF samples names from one of the VCFs
 # zcat chrY.Survivor_WGS.GATK4180.hg38.vcf.gz | head -10000 | grep "^#CHROM" | awk '{for (i=10; i<=NF; i++) print $i}' > VCFsample_names.txt
 #---------------- Sample QC --------------------------

## Rename 2 samples in the master VCF file, these are incorrectly labelled for wrong samples
bcftools reheader -s Phenotypes/master_vcf_rename.txt VCF_original/chrALL.Survivor_WGS.GATK4180.hg38.vcf.gz > VCF_original/chrALL.Survivor_WGS.GATK4180.hg38_renamed2.vcf.gz
tabix -pvcf VCF_original/chrALL.Survivor_WGS.GATK4180.hg38_renamed.vcf.gz

WORKDIR="/research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed" 
## Per-sample missingness, using the GIANT VCF file
## This command is used to calculate the per-sample and per-variant missingness in a large VCF file (chrALL.Survivor_WGS.GATK4180.hg38.vcf.gz) using PLINK. It runs on a high-performance computing (HPC) cluster, utilizing 100 threads to expedite the process. The job is submitted to the priority queue with a high memory allocation (200000 MB) to handle the large dataset. The output, which will include detailed statistics on missing genotypes, is directed to a specified file (QC/chrALL.Survivor_WGS.GATK4180.hg38_missingness). The purpose of this analysis is to identify and quantify the extent of missing data, which is crucial for ensuring data quality and informing further steps in the analysis pipeline, such as filtering out low-quality samples or variants.
bsub -q priority -P SurvivorGWAS -J missing -eo log/missing.err -oo log/missing.out -R "rusage[mem=200000]" "plink --nonfounders \
 --vcf VCF_original/chrALL.Survivor_WGS.GATK4180.hg38_renamed.vcf.gz \
 --threads 100 \
 --missing \
 --out QC/chrALL.Survivor_WGS.GATK4180.hg38_missingness2"
awk 'FNR==1 || $6 > 0.05 {print $1"_"$2}' QC/chrALL.Survivor_WGS.GATK4180.hg38_missingness.imiss > QC/chrALL.Survivor_WGS.GATK4180.hg38_missingness.imiss.0.05plus



## Download 
Please note the exclusion of high-LD regions before performing admixture analyses (this is also true for performing PCA and IBD analyses). 
See https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD) for high LD regions based on hg38 coordinates to exclude.

cp /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS/hg38.fa* .
gatk CreateSequenceDictionary -R hg38.fa -O hg38.dict

## For chr1-22, perform sequence-level QC, fix multiallelic variants, run gross HWE check, and run LD pruning
for chr in {1..22} X Y;do
  export chr=${chr}
  bsub -q priority -P SurvivorGWAS -J QC1.chr$chr -eo log/QC1.$chr.err -oo log/QC1.$chr.out -R "rusage[mem=20000]" bash process_Survivor_WGS_data_perchr.sh $chr
done

## Process these files to include biallelic markers only to perform additional sample QCs
cd QC/
for chr in {1..22};do
 # Get biallelic variants
 awk 'a[$4]++' Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep.bim \
 | awk 'NR==FNR{a[$4];next}!($4 in a){print $2}' \
 - Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep.bim \
 > Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic
 # Harmonize variant names for chr:pos
 sed 's/^chr//g' Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic \
 | awk -F':' '{print $1":"$2}' \
 | paste  Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic \
 - > Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic_update_var_names
 # Make a list to extract
 awk '{print $2}' Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic_update_var_names \
  > Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic_update_var_names.list
 plink --nonfounders \
  --bfile Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep \
  --extract Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic \
  --make-bed \
  --out Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep_biallelic
 plink --nonfounders \
  --bfile Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep_biallelic \
  --update-name Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic_update_var_names \
  --make-bed \
  --out Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated
done
# Make one file including chr1-22
for chr in {2..22}; do
 echo Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated.bed \
 Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated.bim \
 Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated.fam >> merge_list_SJLIFE.txt
done
plink --nonfounders \
 --bfile Survivor_WGS.GATK4180.hg38_renamed_chr1.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated \
 --merge-list merge_list_SJLIFE.txt \
 --make-bed \
 --out Survivor_WGS.GATK4180.hg38_renamed_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated

## Sex discordance check is done by Zhaoming [as he said], and we do not have genotype data for X-linked markers, this step is skipped and assumed there are no sex discordances.

## Per-sample heterozygosity check
plink --nonfounders \
 --bfile Survivor_WGS.GATK4180.hg38_renamed_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated \
 --het \
 --out Survivor_WGS.GATK4180.hg38_renamed_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated_het
R CMD BATCH plot-imiss-vs-het_Survivor_WGS.R
# This will produce a list of samples with outlying heterozygosity
# Samples with excess heterozygosity can be due to sample contamination and hence make a list with >0.4 hetrate for exclusion from the main dataset
# Samples with excess homozygoisity can be due to inbreeding or other stuffs but may not necessarily due to sample contamination, and hence make a list with <=0.4 hetrate as a flag list
mkdir -p dropsamples
awk 'FNR>1 && $NF>0.4{print $1, $2}' Survivor_WGS_Per_sample_heterozygosity_outlier_check.txt > dropsamples/excessHet.drop.samples
awk 'FNR>1 && $NF<=0.4{print $1, $2}' Survivor_WGS_Per_sample_heterozygosity_outlier_check.txt > dropsamples/excessHom.flag.samples

## Cryptic relatedness check
plink --nonfounders \
 --bfile Survivor_WGS.GATK4180.hg38_renamed_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated \
 --remove dropsamples/excessHet.drop.samples \
 --genome \
 --out Survivor_WGS.GATK4180.hg38_renamed_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated_ibd


# Potentially related pairs of samples - keep them as a flag list only
awk 'FNR==1 || $10 > 0.2' Survivor_WGS.GATK4180.hg38_renamed_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated_ibd.genome \
 > dropsamples/relatedness.flag.samples

#---------------- PCA --------------------------

mkdir -p pca
cd pca

# Extract genotype data for these SNPs from 1000G data - use hg38 because Survivor_WGS data is in hg38 coordinates
for chr in {1..22};do
 plink --nonfounders \
  --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/1kGP/plink/ALL.chr${chr}_GRCh38.genotypes.20170504_biallelic_uniq_chrpos \
  --extract ../Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic_update_var_names.list \
  --make-bed \
  --out ALL.chr${chr}_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SurvivorWGS
done

# Combine genotype data from chr1-22 into one file
rm merge_list_1kGP.txt
for chr in {2..22}; do
 echo ALL.chr${chr}_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SurvivorWGS.bed \
 ALL.chr${chr}_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SurvivorWGS.bim \
 ALL.chr${chr}_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SurvivorWGS.fam >> merge_list_1kGP.txt
done
plink --nonfounders \
 --bfile ALL.chr1_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SurvivorWGS \
 --merge-list merge_list_1kGP.txt \
 --make-bed \
 --out ALL.chr1-22_GRCh38.genotypes.20240923_biallelic_uniq_chrpos_Survivor_WGS
# Merge SJLFIE and 1kGP datasets
plink --nonfounders \
 --bfile ../Survivor_WGS.GATK4180.hg38_renamed_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated \
 --bmerge ALL.chr1-22_GRCh38.genotypes.20240923_biallelic_uniq_chrpos_Survivor_WGS.bed ALL.chr1-22_GRCh38.genotypes.20240923_biallelic_uniq_chrpos_Survivor_WGS.bim ALL.chr1-22_GRCh38.genotypes.20240923_biallelic_uniq_chrpos_Survivor_WGS.fam \
 --make-bed \
 --out Survivor_WGS_and_1kGP_tmp1
# Exclude variants with multiple alleles
for file in ../Survivor_WGS.GATK4180.hg38_renamed_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated ALL.chr1-22_GRCh38.genotypes.20240923_biallelic_uniq_chrpos_Survivor_WGS; do
 plink --nonfounders \
  --bfile $file \
  --exclude Survivor_WGS_and_1kGP_tmp1-merge.missnp \
  --make-bed \
  --out $file.clean
done
# Now merge again
plink --nonfounders \
 --bfile ../Survivor_WGS.GATK4180.hg38_renamed_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated.clean \
 --bmerge \
 ALL.chr1-22_GRCh38.genotypes.20240923_biallelic_uniq_chrpos_Survivor_WGS.clean.bed ALL.chr1-22_GRCh38.genotypes.20240923_biallelic_uniq_chrpos_Survivor_WGS.clean.bim ALL.chr1-22_GRCh38.genotypes.20240923_biallelic_uniq_chrpos_Survivor_WGS.clean.fam \
 --geno 0.1 \
 --make-bed \
 --out Survivor_WGS_and_1kGP_final

# Perform PCA and obtain top 20 PCs
plink --nonfounders \
 --bfile Survivor_WGS_and_1kGP_final \
 --remove ../dropsamples/excessHet.drop.samples \
 --exclude range /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/high-LD-regions-GRCh38.txt \
 --range \
 --pca 20 \
 --out Survivor_WGS_and_1kGP_final_top_20_PCs
# Plot and identfiication of genetic ancestry is done using R script [~/bin/PCA_plot_for_SurvivorWGS_samples.R] - this script generates list of individuals with EUR and AFR acnestries (based on genotype data).
Rscript ~/bin/PCA_plot_for_SurvivorWGS_samples.R

# Compute top 20 PCs for Europeans and Africans alone
plink --nonfounders \
 --bfile Survivor_WGS_and_1kGP_final \
 --keep Survivor_WGS_EUR_based_on_1kGP_Phase_3_data.txt \
 --exclude range /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/high-LD-regions-GRCh38_no_chr.txt \
 --pca 20 \
 --out Survivor_WGS_and_1kGP_final_EUR_top_20_PCs
 
 plink --nonfounders \
 --bfile Survivor_WGS_and_1kGP_final \
 --keep Survivor_WGS_AFR_based_on_1kGP_Phase_3_data.txt \
 --exclude range /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/high-LD-regions-GRCh38_no_chr.txt \
 --pca 20 \
 --out Survivor_WGS_and_1kGP_final_AFR_top_20_PCs

# Also compute top 20 PCs for all the 2986 Survivor_WGS samples, just in case; however analysis of all samples is not recommended
plink --nonfounders \
 --bfile Survivor_WGS_and_1kGP_final \
 --keep ../Survivor_WGS.GATK4180.hg38_renamed_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated.clean.fam \
 --remove ../dropsamples/excessHet.drop.samples \
 --pca 20 \
 --out Survivor_WGS_and_1kGP_final_all_Survivor_WGS_samples_20_PCs

#---------------- Rename PCGP ids in the VCF files to Survivor_WGS ids --------------------------

# Rename sample ids in the VCF file to Survivor_WGS ids, while excluding the 20 samples with excess heetrozygosity
cd ../../VCF_original
awk 'FNR>1{if($4!~/SJ/) vcfid=$NF; else vcfid=$3"-"$2;print vcfid, $1}' ../Phenotypes/VCF_ids_to_SJLIDs.map_edited > ../Phenotypes/SurvivorWGS_3006_vcfid_to_sjlid.linkfile
awk '{print $1"_"$2}' ../QC/dropsamples/excessHet.drop.samples \
 | grep -F -w -f - ../Phenotypes/SurvivorWGS_3006_vcfid_to_sjlid.linkfile \
 | awk '{print $2}' \
 > ../Phenotypes/SurvivorWGS_20_excessHet_to_drop_sjlid.sample
# Also prepare PC files including SJLIDS
awk '{if($1~/SJ/) $1=$1"_"$2; else $1=$1; print}' ../QC/pca/Survivor_WGS_and_1kGP_final_EUR_top_20_PCs.eigenvec \
 | awk 'NR==FNR{a[$1]=$2;next}($1 in a){$1=$2=a[$1]; print}' ../Phenotypes/SurvivorWGS_3006_vcfid_to_sjlid.linkfile - \
 > ../Phenotypes/Survivor_WGS_and_1kGP_final_EUR_top_20_PCs.eigenvec.sjlid
awk '{if($1~/SJ/) $1=$1"_"$2; else $1=$1; print}' ../QC/pca/Survivor_WGS_and_1kGP_final_AFR_top_20_PCs.eigenvec \
 | awk 'NR==FNR{a[$1]=$2;next}($1 in a){$1=$2=a[$1]; print}' ../Phenotypes/SurvivorWGS_3006_vcfid_to_sjlid.linkfile - \
 > ../Phenotypes/Survivor_WGS_and_1kGP_final_AFR_top_20_PCs.eigenvec.sjlid 
 
for chr in {1..22}; do
 bsub -q priority -P SurvivorGWAS -J SurvivorWGS_rename$chr -eo ../log/SurvivorWGS_rename$chr.err -oo ../log/SurvivorWGS_rename$chr.out -R "rusage[mem=3000]" \
 "bcftools reheader -s ../Phenotypes/SurvivorWGS_3006_vcfid_to_sjlid.linkfile Survivor_WGS.GATK4180.hg38_renamed_chr${chr}.PASS.decomposed.vcf.gz \
  | bcftools view -S ^../Phenotypes/SurvivorWGS_20_excessHet_to_drop.sample.sjlid \
  > Survivor_WGS.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.vcf.gz"
done
# Exclude variants with more than 10% missingness
for chr in {1..22}; do 
 bsub -q priority -P SurvivorGWAS -J SurvivorWGS_all_missing$chr -eo ../log/SurvivorWGS_all_missing$chr.err -oo ../log/SurvivorWGS_all_missing$chr.out -R "rusage[mem=3000]" \
 "vcftools \
 --gzvcf Survivor_WGS.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.vcf.gz \
 --max-missing 0.90 \
 --recode \
 --stdout \
 | bgzip \
 > Survivor_WGS.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.qced.vcf.gz"
done
 
# Write stats
for chr in {1..22}; do 
 bsub -P SurvivorGWAS -J SurvivorWGS_all_qced_stats$chr -eo ../log/SurvivorWGS_all_qced_stats$chr.err -oo ../log/SurvivorWGS_all_qced_stats$chr.out -R "rusage[mem=3000]" \
  "bcftools stats Survivor_WGS.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.qced.vcf.gz \
  > ../QC/stats/Survivor_WGS.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.qced.bcftools_stats"
done

# Index the VCF files
for chr in {1..22}; do
 bsub -P SurvivorGWAS -J SurvivorWGS_idx_all$chr -eo ../log/SurvivorWGS_idx_all$chr.err -oo ../log/SurvivorWGS_idx_all$chr.out -R "rusage[mem=3000]" \
  "tabix -pvcf Survivor_WGS.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.qced.vcf.gz"
done

#---------------- Prepare VCF files for EUR and AFR separately, excluding the 20 samples with excess heterozygosity --------------------------
mkdir -p EUR
mkdir -p AFR

# Prepare sample files for both EUR and AFR samples
awk 'FNR>1{if($1!~/SJ/) $1=$1; else $1=$1"_"$2; print $1}' ../QC/pca/SurvivorWGS_EUR_based_on_1kGP_Phase_3_data.txt \
 | awk 'NR==FNR{a[$1];next}($1 in a){print $2}' - ../Phenotypes/SurvivorWGS_3006_vcfid_to_sjlid.linkfile \
 > ../Phenotypes/SurvivorWGS_EUR_PCA_1kGP.sjlid
awk 'FNR>1{if($1!~/SJ/) $1=$1; else $1=$1"_"$2; print $1}' ../QC/pca/SurvivorWGS_AFR_based_on_1kGP_Phase_3_data.txt \
 | awk 'NR==FNR{a[$1];next}($1 in a){print $2}' - ../Phenotypes/SurvivorWGS_3006_vcfid_to_sjlid.linkfile \
 > ../Phenotypes/SurvivorWGS_AFR_PCA_1kGP.sjlid 

# Now subset VCF files for EUR and AFR using the above sample files, while excluding variants with >10% missingness and not in HWE (P < 1e-10)
for chr in {1..22}; do
 bsub -P SurvivorGWAS -J SurvivorWGS_eur$chr -eo ../log/SurvivorWGS_eur$chr.err -oo ../log/SurvivorWGS_eur$chr.out -R "rusage[mem=3000]" \
  "vcftools \
   --gzvcf Survivor_WGS.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.vcf.gz \
   --keep ../Phenotypes/SurvivorWGS_EUR_PCA_1kGP.sjlid \
   --max-missing 0.90 \
   --hwe 1e-10 \
   --recode \
   --stdout \
   | bgzip \
   > EUR/Survivor_WGS.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_eur.qced.vcf.gz"
 bsub -P SurvivorGWAS -J SurvivorWGS_afr$chr -eo ../log/SurvivorWGS_afr$chr.err -oo ../log/SurvivorWGS_afr$chr.out -R "rusage[mem=3000]" \
  "vcftools \
   --gzvcf Survivor_WGS.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.vcf.gz \
   --keep ../Phenotypes/SurvivorWGS_AFR_PCA_1kGP.sjlid \
   --max-missing 0.90 \
   --hwe 1e-10 \
   --recode \
   --stdout \
   | bgzip \
   > AFR/Survivor_WGS.GERMLINE.412.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_afr.qced.vcf.gz"
done

# Write stats
for chr in {1..22}; do 
 bsub -P SurvivorGWAS -J SurvivorWGS_eur_qced_stats$chr -eo ../log/SurvivorWGS_eur_qced_stats$chr.err -oo ../log/SurvivorWGS_eur_qced_stats$chr.out -R "rusage[mem=3000]" \
  "bcftools stats EUR/Survivor_WGS.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_eur.qced.vcf.gz \
  > ../QC/stats/Survivor_WGS.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_eur.qced.bcftools_stats"
  bsub -P SurvivorGWAS -J SurvivorWGS_afr_qced_stats$chr -eo ../log/SurvivorWGS_afr_qced_stats$chr.err -oo ../log/SurvivorWGS_afr_qced_stats$chr.out -R "rusage[mem=3000]" \
  "bcftools stats AFR/Survivor_WGS.GERMLINE.412.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_afr.qced.vcf.gz \
  > ../QC/stats/Survivor_WGS.GERMLINE.412.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_afr.qced.bcftools_stats"
done

# Index the VCF files
for chr in {1..22}; do
 bsub -P SurvivorGWAS -J SurvivorWGS_idx_eur$chr -eo ../log/SurvivorWGS_idx_eur$chr.err -oo ../log/SurvivorWGS_idx_eur$chr.out -R "rusage[mem=3000]" \
  "tabix -pvcf EUR/Survivor_WGS.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_eur.qced.vcf.gz"
 bsub -P SurvivorGWAS -J SurvivorWGS_idx_afr -eo ../log/SurvivorWGS_idx_afr.err -oo ../log/SurvivorWGS_idx_afr.out -R "rusage[mem=3000]" \
  "tabix -pvcf AFR/Survivor_WGS.GERMLINE.412.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_afr.qced.vcf.gz" 
done

#---------------- Prepare corresponding PLINK files ----------------
mkdir -p PLINK/EUR
mkdir -p PLINK/AFR
mkdir -p PLINK/ALL

for chr in {1..22}; do
 bsub -q priority -P SurvivorGWAS -J SurvivorWGS_plink$chr -eo ../log/SurvivorWGS_plink$chr.err -oo ../log/SurvivorWGS_plink$chr.out -R "rusage[mem=3000]" \
  "plink --vcf Survivor_WGS.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.qced.vcf.gz \
   --keep-allele-order \
   --allow-extra-chr \
   --make-bed \
   --out PLINK/ALL/Survivor_WGS.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid;
  plink --vcf EUR/Survivor_WGS.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_eur.qced.vcf.gz \
   --keep-allele-order \
   --allow-extra-chr \
   --make-bed \
   --out PLINK/EUR/Survivor_WGS.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_eur;
  plink --vcf AFR/Survivor_WGS.GERMLINE.412.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_afr.qced.vcf.gz \
   --keep-allele-order \
   --allow-extra-chr \
   --make-bed \
   --out PLINK/AFR/Survivor_WGS.GERMLINE.412.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_afr" 
done

#---------------- Prepare group files ----------------
cd ../Analyses
bash scripts/preprocess.sh




# # Also remove
# SJNORM041129_G1 REMOVE-not SJLIFE       REMOVE  NA      SJNORM041129_G1-TB-15-5985      REMOVE
# SJNORM041129_G2 REMOVE-not SJLIFE       REMOVE  NA      SJNORM041129_G2-TB-15-6901      REMOVE
# SJNORM041130_G1 REMOVE-not SJLIFE       REMOVE  NA      SJNORM041130_G1-TB-15-6225      REMOVE
# SJNORM041130_G2 REMOVE-not SJLIFE       REMOVE  NA      SJNORM041130_G2-TB-15-6902      REMOVE
# SJNORM041131_G1 REMOVE-not SJLIFE       REMOVE  NA      SJNORM041131_G1-TB-15-6369      REMOVE
# SJNORM041131_G2 REMOVE-not SJLIFE       REMOVE  NA      SJNORM041131_G2-TB-15-6903      REMOVE
# SJNORM041132_G1 REMOVE-not SJLIFE       REMOVE  NA      SJNORM041132_G1-TB-15-6370      REMOVE
# SJNORM041132_G2 REMOVE-not SJLIFE       REMOVE  NA      SJNORM041132_G2-TB-15-6904      REMOVE
# SJNORM041133_G1 REMOVE-not SJLIFE       REMOVE  NA      SJNORM041133_G1-TB-15-6399      REMOVE
# SJNORM041133_G2 REMOVE-not SJLIFE       REMOVE  NA      SJNORM041133_G2-TB-15-6905      REMOVE






# ###################
# ## Run admixture ##
# ###################

# module load plink/1.90b

# #SJLIFE data
# bsub -R "rusage[mem=320000]" -P plinkmergeSJLIFE plink --memory 320000 --keep-allele-order --merge-list allfiles.txt --maf 0.05 --hwe 1e-06 --geno 0.05 --make-bed --out SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05
# bsub -R "rusage[mem=120000]" -P plinkmergeSJLIFE plink --bfile chr.ALL.Survivor_WES.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_ID_updated --keep-allele-order --allow-extra-chr --exclude range /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS_QCed/high-LD-regions-GRCh38.txt --geno 0.01 --hwe 0.0001 --maf 0.05 --make-bed --indep-pairwise 100 25 0.2 --out All.survivors_indep_pairwise
# #1000Genome Data
# # bsub -R "rusage[mem=320000]" -P plinkmergeSJLIFE plink --memory 320000 --keep-allele-order --allow-extra-chr --geno 0.05 --hwe 1e-06 --maf 0.05 --make-bed --out 1000genomes_merged_maf0.05_hwe1e-06_geno0.05 --vcf merged1000genomes.vcf --vcf-half-call h
# # bsub -R "rusage[mem=320000]" -P plinkmergeSJLIFE plink --bfile 1000genomes_merged_maf0.05_hwe1e-06_geno0.05 --memory 320000 --keep-allele-order --allow-extra-chr --make-bed --indep-pairwise 1500 30 0.3 --out 1000genomes_merged_maf0.05_hwe1e-06_geno0.05_indep_pairwise

# #Merge 1000 genomes data with SJLIFE data
# bsub -R "rusage[mem=320000]" -P plinkmergeSJLIFE plink --bfile All.survivors_indep_pairwise --bmerge /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/1kGP/1000genomes_merged_maf0.05_hwe1e-06_geno0.05 --keep-allele-order --allow-extra-chr --make-bed --memory 320000 --out 1000genomes_All.survivors_indep_pairwise_maf0.05_geno0.05_indep_pairwise
# '''
# 1031903 MB RAM detected; reserving 320000 MB for main workspace.
# 4481 people loaded from
# SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05_indep_pairwise.fam.
# 2504 people to be merged from
# ../../../1kGP/1000genomes_merged_maf0.05_hwe1e-06_geno0.05.fam.
# Of these, 2504 are new, while 0 are present in the base dataset.
# 5185765 markers loaded from
# SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05_indep_pairwise.bim.
# 5120204 markers to be merged from
# ../../../1kGP/1000genomes_merged_maf0.05_hwe1e-06_geno0.05.bim.
# Of these, 5119558 are new, while 646 are present in the base dataset.
# Error: 560 variants with 3+ alleles present.
# * If you believe this is due to strand inconsistency, try --flip with
#   1000genomes_SJLIFE_chromosomes_combine_maf0.05_hwe1e-06_geno0.05_indep_pairwise-merge.missnp.
#   (Warning: if this seems to work, strand errors involving SNPs with A/T or C/G
#   alleles probably remain in your data.  If LD between nearby SNPs is high,
#   --flip-scan should detect them.)
# * If you are dealing with genuine multiallelic variants, we recommend exporting
#   that subset of the data to VCF (via e.g. '--recode vcf'), merging with
#   another tool/script, and then importing the result; PLINK is not yet suited
#   to handling them.
# '''