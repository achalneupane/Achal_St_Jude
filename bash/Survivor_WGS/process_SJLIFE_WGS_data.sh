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
 --vcf chrALL.Survivor_WGS.GATK4180.hg38_renamed.vcf.gz \
 --threads 100 \
 --missing \
 --out QC/chrALL.Survivor_WGS.GATK4180.hg38_missingness"
awk 'FNR==1 || $6 > 0.05 {print $1"_"$2}' QC/chrALL.Survivor_WGS.GATK4180.hg38_missingness.imiss > QC/chrALL.Survivor_WGS.GATK4180.hg38_missingness.imiss.0.05plus



## Download 
Please note the exclusion of high-LD regions before performing admixture analyses (this is also true for performing PCA and IBD analyses). 
See https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD) for high LD regions based on hg38 coordinates to exclude.

cp /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WGS/hg38.fa* .
gatk CreateSequenceDictionary -R hg38.fa -O hg38.dict

## For chr1-22, perform sequence-level QC, fix multiallelic variants, run gross HWE check, and run LD pruning
for chr in {1..22} X Y;do
  export chr=${chr}
  bsub -q priority -P SurvivorGWAS -J QC1.chr$chr -eo log/QC1.$chr.err -oo log/QC1.$chr.out -R "rusage[mem=20000]" bash process_SJLIFE_WGS_data_perchr.sh $chr
done

## Process these files to include biallelic markers only to perform additional sample QCs
cd QC/
for chr in {1..22};do
 # Get biallelic variants
 awk 'a[$4]++' SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep.bim \
 | awk 'NR==FNR{a[$4];next}!($4 in a){print $2}' \
 - SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep.bim \
 > SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic
 # Harmonize variant names for chr:pos
 sed 's/^chr//g' SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic \
 | awk -F':' '{print $1":"$2}' \
 | paste  SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic \
 - > SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic_update_var_names
 # Make a list to extract
 awk '{print $2}' SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic_update_var_names \
  > SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic_update_var_names.list
 plink --nonfounders \
  --bfile SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep \
  --extract SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic \
  --make-bed \
  --out SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep_biallelic
 plink --nonfounders \
  --bfile SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep_biallelic \
  --update-name SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic_update_var_names \
  --make-bed \
  --out SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated
done
# Make one file including chr1-22
for chr in {2..22}; do
 echo SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated.bed \
 SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated.bim \
 SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated.fam >> merge_list_SJLIFE.txt
done
plink --nonfounders \
 --bfile SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr1.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated \
 --merge-list merge_list_SJLIFE.txt \
 --make-bed \
 --out SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated

## Sex discordance check is done by Zhaoming [as he said], and we do not have genotype data for X-linked markers, this step is skipped and assumed there are no sex discordances.

## Per-sample heterozygosity check
plink --nonfounders \
 --bfile SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated \
 --het \
 --out SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated_het
R CMD BATCH plot-imiss-vs-het_SJLIFE_WGS.R
# This will produce a list of samples with outlying heterozygosity
# Samples with excess heterozygosity can be due to sample contamination and hence make a list with >0.4 hetrate for exclusion from the main dataset
# Samples with excess homozygoisity can be due to inbreeding or other stuffs but may not necessarily due to sample contamination, and hence make a list with <=0.4 hetrate as a flag list
mkdir -p dropsamples
awk 'FNR>1 && $NF>0.4{print $1, $2}' SJLIFE_WGS_Per_sample_heterozygosity_outlier_check.txt > dropsamples/excessHet.drop.samples
awk 'FNR>1 && $NF<=0.4{print $1, $2}' SJLIFE_WGS_Per_sample_heterozygosity_outlier_check.txt > dropsamples/excessHom.flag.samples

## Cryptic relatedness check
plink --nonfounders \
 --bfile SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated \
 --remove dropsamples/excessHet.drop.samples \
 --genome \
 --out SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated_ibd
# Potentially related pairs of samples - keep them as a flag list only
awk 'FNR==1 || $10 > 0.2' SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated_ibd.genome \
 > dropsamples/relatedness.flag.samples

#---------------- PCA --------------------------

mkdir -p pca
cd pca

# Extract genotype data for these SNPs from 1000G data - use hg38 because SJLIFE data is in hg38 coordinates
for chr in {1..22};do
 plink --nonfounders \
  --bfile ~/Work/1000G_phase_3/20130502/hg38/PLINKformat/ALL.chr${chr}_GRCh38.genotypes.20170504_biallelic_uniq_chrpos \
  --extract ../SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed_common_pruned_indep.bim_biallelic_update_var_names.list \
  --make-bed \
  --out ALL.chr${chr}_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SJLIFEWGS
done
# Combine genotype data from chr1-22 into one file
for chr in {2..22}; do
 echo ALL.chr${chr}_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SJLIFEWGS.bed \
 ALL.chr${chr}_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SJLIFEWGS.bim \
 ALL.chr${chr}_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SJLIFEWGS.fam >> merge_list_1kGP.txt
done
plink --nonfounders \
 --bfile ALL.chr1_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SJLIFEWGS \
 --merge-list merge_list_1kGP.txt \
 --make-bed \
 --out ALL.chr1-22_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SJLIFEWGS
# Merge SJLFIE and 1kGP datasets
plink --nonfounders \
 --bfile ../SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated \
 --bmerge ALL.chr1-22_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SJLIFEWGS.bed ALL.chr1-22_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SJLIFEWGS.bim ALL.chr1-22_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SJLIFEWGS.fam \
 --make-bed \
 --out SJLIFE_and_1kGP_tmp1
# Exclude variants with multiple alleles
for file in ../SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated ALL.chr1-22_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SJLIFEWGS; do
 plink --nonfounders \
  --bfile $file \
  --exclude SJLIFE_and_1kGP_tmp1-merge.missnp \
  --make-bed \
  --out $file.clean
done
# Now merge again
plink --nonfounders \
 --bfile ../SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated.clean \
 --bmerge \
 ALL.chr1-22_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SJLIFEWGS.clean.bed ALL.chr1-22_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SJLIFEWGS.clean.bim ALL.chr1-22_GRCh38.genotypes.20170504_biallelic_uniq_chrpos_SJLIFEWGS.clean.fam \
 --geno 0.1 \
 --make-bed \
 --out SJLIFE_and_1kGP_final

# Perform PCA and obtain top 20 PCs
plink --nonfounders \
 --bfile SJLIFE_and_1kGP_final \
 --remove ../dropsamples/excessHet.drop.samples \
 --exclude range ~/bin/high-ld_regions_to_exclude_hg38.txt \
 --range \
 --pca 20 \
 --out SJLIFE_and_1kGP_final_top_20_PCs
# Plot and identfiication of genetic ancestry is done using R script [~/bin/PCA_plot_for_SJLIFEWGS_samples.R] - this script generates list of individuals with EUR and AFR acnestries (based on genotype data).
Rscript ~/bin/PCA_plot_for_SJLIFEWGS_samples.R

# Compute top 20 PCs for Europeans and Africans alone
plink --nonfounders \
 --bfile SJLIFE_and_1kGP_final \
 --keep SJLIFEWGS_EUR_based_on_1kGP_Phase_3_data.txt \
 --exclude range ~/bin/high-ld_regions_to_exclude_hg38.txt \
 --pca 20 \
 --out SJLIFE_and_1kGP_final_EUR_top_20_PCs
 
 plink --nonfounders \
 --bfile SJLIFE_and_1kGP_final \
 --keep SJLIFEWGS_AFR_based_on_1kGP_Phase_3_data.txt \
 --exclude range ~/bin/high-ld_regions_to_exclude_hg38.txt \
 --pca 20 \
 --out SJLIFE_and_1kGP_final_AFR_top_20_PCs

# Also compute top 20 PCs for all the 2986 SJLIFE samples, just in case; however analysis of all samples is not recommended
plink --nonfounders \
 --bfile SJLIFE_and_1kGP_final \
 --keep ../SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr1-22.PASS.decomposed_common_pruned_indep_biallelic_varnames_updated.clean.fam \
 --remove ../dropsamples/excessHet.drop.samples \
 --pca 20 \
 --out SJLIFE_and_1kGP_final_all_SJLIFE_samples_20_PCs

#---------------- Rename PCGP ids in the VCF files to SJLIFE ids --------------------------

# Rename sample ids in the VCF file to SJLIFE ids, while excluding the 20 samples with excess heetrozygosity
cd ../../VCF_original
awk 'FNR>1{if($4!~/SJ/) vcfid=$NF; else vcfid=$3"-"$2;print vcfid, $1}' ../Phenotypes/VCF_ids_to_SJLIDs.map_edited > ../Phenotypes/SJLIFEWGS_3006_vcfid_to_sjlid.linkfile
awk '{print $1"_"$2}' ../QC/dropsamples/excessHet.drop.samples \
 | grep -F -w -f - ../Phenotypes/SJLIFEWGS_3006_vcfid_to_sjlid.linkfile \
 | awk '{print $2}' \
 > ../Phenotypes/SJLIFEWGS_20_excessHet_to_drop_sjlid.sample
# Also prepare PC files including SJLIDS
awk '{if($1~/SJ/) $1=$1"_"$2; else $1=$1; print}' ../QC/pca/SJLIFE_and_1kGP_final_EUR_top_20_PCs.eigenvec \
 | awk 'NR==FNR{a[$1]=$2;next}($1 in a){$1=$2=a[$1]; print}' ../Phenotypes/SJLIFEWGS_3006_vcfid_to_sjlid.linkfile - \
 > ../Phenotypes/SJLIFE_and_1kGP_final_EUR_top_20_PCs.eigenvec.sjlid
awk '{if($1~/SJ/) $1=$1"_"$2; else $1=$1; print}' ../QC/pca/SJLIFE_and_1kGP_final_AFR_top_20_PCs.eigenvec \
 | awk 'NR==FNR{a[$1]=$2;next}($1 in a){$1=$2=a[$1]; print}' ../Phenotypes/SJLIFEWGS_3006_vcfid_to_sjlid.linkfile - \
 > ../Phenotypes/SJLIFE_and_1kGP_final_AFR_top_20_PCs.eigenvec.sjlid 
 
for chr in {1..22}; do
 bsub -q priority -P SurvivorGWAS -J sjlifewgs_rename$chr -eo ../log/sjlifewgs_rename$chr.err -oo ../log/sjlifewgs_rename$chr.out -R "rusage[mem=3000]" \
 "bcftools reheader -s ../Phenotypes/SJLIFEWGS_3006_vcfid_to_sjlid.linkfile SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.vcf.gz \
  | bcftools view -S ^../Phenotypes/SJLIFEWGS_20_excessHet_to_drop.sample.sjlid \
  > SJLIFE.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.vcf.gz"
done
# Exclude variants with more than 10% missingness
for chr in {1..22}; do 
 bsub -q priority -P SurvivorGWAS -J sjlifewgs_all_missing$chr -eo ../log/sjlifewgs_all_missing$chr.err -oo ../log/sjlifewgs_all_missing$chr.out -R "rusage[mem=3000]" \
 "vcftools \
 --gzvcf SJLIFE.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.vcf.gz \
 --max-missing 0.90 \
 --recode \
 --stdout \
 | bgzip \
 > SJLIFE.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.qced.vcf.gz"
done
 
# Write stats
for chr in {1..22}; do 
 bsub -P SurvivorGWAS -J sjlifewgs_all_qced_stats$chr -eo ../log/sjlifewgs_all_qced_stats$chr.err -oo ../log/sjlifewgs_all_qced_stats$chr.out -R "rusage[mem=3000]" \
  "bcftools stats SJLIFE.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.qced.vcf.gz \
  > ../QC/stats/SJLIFE.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.qced.bcftools_stats"
done

# Index the VCF files
for chr in {1..22}; do
 bsub -P SurvivorGWAS -J sjlifewgs_idx_all$chr -eo ../log/sjlifewgs_idx_all$chr.err -oo ../log/sjlifewgs_idx_all$chr.out -R "rusage[mem=3000]" \
  "tabix -pvcf SJLIFE.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.qced.vcf.gz"
done

#---------------- Prepare VCF files for EUR and AFR separately, excluding the 20 samples with excess heterozygosity --------------------------
mkdir -p EUR
mkdir -p AFR

# Prepare sample files for both EUR and AFR samples
awk 'FNR>1{if($1!~/SJ/) $1=$1; else $1=$1"_"$2; print $1}' ../QC/pca/SJLIFEWGS_EUR_based_on_1kGP_Phase_3_data.txt \
 | awk 'NR==FNR{a[$1];next}($1 in a){print $2}' - ../Phenotypes/SJLIFEWGS_3006_vcfid_to_sjlid.linkfile \
 > ../Phenotypes/SJLIFEWGS_EUR_PCA_1kGP.sjlid
awk 'FNR>1{if($1!~/SJ/) $1=$1; else $1=$1"_"$2; print $1}' ../QC/pca/SJLIFEWGS_AFR_based_on_1kGP_Phase_3_data.txt \
 | awk 'NR==FNR{a[$1];next}($1 in a){print $2}' - ../Phenotypes/SJLIFEWGS_3006_vcfid_to_sjlid.linkfile \
 > ../Phenotypes/SJLIFEWGS_AFR_PCA_1kGP.sjlid 

# Now subset VCF files for EUR and AFR using the above sample files, while excluding variants with >10% missingness and not in HWE (P < 1e-10)
for chr in {1..22}; do
 bsub -P SurvivorGWAS -J sjlifewgs_eur$chr -eo ../log/sjlifewgs_eur$chr.err -oo ../log/sjlifewgs_eur$chr.out -R "rusage[mem=3000]" \
  "vcftools \
   --gzvcf SJLIFE.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.vcf.gz \
   --keep ../Phenotypes/SJLIFEWGS_EUR_PCA_1kGP.sjlid \
   --max-missing 0.90 \
   --hwe 1e-10 \
   --recode \
   --stdout \
   | bgzip \
   > EUR/SJLIFE.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_eur.qced.vcf.gz"
 bsub -P SurvivorGWAS -J sjlifewgs_afr$chr -eo ../log/sjlifewgs_afr$chr.err -oo ../log/sjlifewgs_afr$chr.out -R "rusage[mem=3000]" \
  "vcftools \
   --gzvcf SJLIFE.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.vcf.gz \
   --keep ../Phenotypes/SJLIFEWGS_AFR_PCA_1kGP.sjlid \
   --max-missing 0.90 \
   --hwe 1e-10 \
   --recode \
   --stdout \
   | bgzip \
   > AFR/SJLIFE.GERMLINE.412.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_afr.qced.vcf.gz"
done

# Write stats
for chr in {1..22}; do 
 bsub -P SurvivorGWAS -J sjlifewgs_eur_qced_stats$chr -eo ../log/sjlifewgs_eur_qced_stats$chr.err -oo ../log/sjlifewgs_eur_qced_stats$chr.out -R "rusage[mem=3000]" \
  "bcftools stats EUR/SJLIFE.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_eur.qced.vcf.gz \
  > ../QC/stats/SJLIFE.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_eur.qced.bcftools_stats"
  bsub -P SurvivorGWAS -J sjlifewgs_afr_qced_stats$chr -eo ../log/sjlifewgs_afr_qced_stats$chr.err -oo ../log/sjlifewgs_afr_qced_stats$chr.out -R "rusage[mem=3000]" \
  "bcftools stats AFR/SJLIFE.GERMLINE.412.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_afr.qced.vcf.gz \
  > ../QC/stats/SJLIFE.GERMLINE.412.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_afr.qced.bcftools_stats"
done

# Index the VCF files
for chr in {1..22}; do
 bsub -P SurvivorGWAS -J sjlifewgs_idx_eur$chr -eo ../log/sjlifewgs_idx_eur$chr.err -oo ../log/sjlifewgs_idx_eur$chr.out -R "rusage[mem=3000]" \
  "tabix -pvcf EUR/SJLIFE.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_eur.qced.vcf.gz"
 bsub -P SurvivorGWAS -J sjlifewgs_idx_afr -eo ../log/sjlifewgs_idx_afr.err -oo ../log/sjlifewgs_idx_afr.out -R "rusage[mem=3000]" \
  "tabix -pvcf AFR/SJLIFE.GERMLINE.412.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_afr.qced.vcf.gz" 
done

#---------------- Prepare corresponding PLINK files ----------------
mkdir -p PLINK/EUR
mkdir -p PLINK/AFR
mkdir -p PLINK/ALL

for chr in {1..22}; do
 bsub -q priority -P SurvivorGWAS -J sjlifewgs_plink$chr -eo ../log/sjlifewgs_plink$chr.err -oo ../log/sjlifewgs_plink$chr.out -R "rusage[mem=3000]" \
  "plink --vcf SJLIFE.GERMLINE.2986.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid.qced.vcf.gz \
   --keep-allele-order \
   --allow-extra-chr \
   --make-bed \
   --out PLINK/ALL/SJLIFE.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid;
  plink --vcf EUR/SJLIFE.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_eur.qced.vcf.gz \
   --keep-allele-order \
   --allow-extra-chr \
   --make-bed \
   --out PLINK/EUR/SJLIFE.GERMLINE.2364.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_eur;
  plink --vcf AFR/SJLIFE.GERMLINE.412.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_afr.qced.vcf.gz \
   --keep-allele-order \
   --allow-extra-chr \
   --make-bed \
   --out PLINK/AFR/SJLIFE.GERMLINE.412.GATKv3.4.vqsr.release.0714_chr${chr}.PASS.decomposed.sjlid_afr" 
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
