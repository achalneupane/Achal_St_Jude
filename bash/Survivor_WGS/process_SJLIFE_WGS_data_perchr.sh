#!/bin/bash

# Processes performed to QC the SJLIFE WGS data provided by Comp. Bio. ["SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714.vcf.gz"]

# Load software packages
module load vcftools/0.1.13
module load bcftools/1.4.1
#module load vt/2016.11.07
module load plink/1.90b
module load R/3.4.0
module load bedtools/2.25.0
module load htslib/1.3.1

# Define/create files/directories
INDIR=/home/ysapkota/Work/WGS_SJLIFE/VCF_original

OUTDIR=/home/ysapkota/Work/WGS_SJLIFE/QC

VCFROOT="SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714"

cd $OUTDIR

# Process data per chromosome
chr=$1

# Split data for each chromosome
vcftools --gzvcf ${INDIR}/${VCFROOT}.vcf.gz --remove SJLIFE.GERMLINE.3006.GATKv3.4.vqsr.release.0714_missingness.imiss.0.05plus --chr chr$chr --recode  --stdout | bgzip  > ${INDIR}/${VCFROOT}_chr${chr}.vcf.gz

# Index the VCF file
tabix -pvcf ${INDIR}/${VCFROOT}_chr${chr}.vcf.gz

# Write stats
bcftools stats ${INDIR}/${VCFROOT}_chr${chr}.vcf.gz > stats/${VCFROOT}_chr${chr}.bcftools_stats

# Run basic QC - sequence level
vcftools --gzvcf ${INDIR}/${VCFROOT}_chr${chr}.vcf.gz \
--chr chr$chr \
--keep-filtered PASS \
--minGQ 20 \
--min-meanDP 10 \
--recode \
--stdout \
| bgzip > ${INDIR}/${VCFROOT}_chr${chr}.PASS.vcf.gz

# Write stats
bcftools stats ${INDIR}/${VCFROOT}_chr${chr}.PASS.vcf.gz > stats/${VCFROOT}_chr${chr}.PASS.bcftools_stats

# Fix multi-allelic markers and normalize
bcftools norm -Ou -m -any ${INDIR}/${VCFROOT}_chr${chr}.PASS.vcf.gz \
 | bcftools norm -Ou -f hg38.fa \
 | bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' \
 | bcftools view -Oz \
 > ${INDIR}/${VCFROOT}_chr${chr}.PASS.decomposed.vcf.gz

# Write stats
bcftools stats ${INDIR}/${VCFROOT}_chr${chr}.PASS.decomposed.vcf.gz > stats/${VCFROOT}_chr${chr}.PASS.decomposed.bcftools_stats

tabix -pvcf ${INDIR}/${VCFROOT}_chr${chr}.PASS.decomposed.vcf.gz

# Also write PLINK files for further processing
plink --vcf ${INDIR}/${VCFROOT}_chr${chr}.PASS.decomposed.vcf.gz \
 --keep-allele-order \
 --allow-extra-chr 0 \
 --hwe 1e-10 \
 --make-bed --out ${VCFROOT}_chr${chr}.PASS.decomposed.qced

# Perform LD pruning to obtain independent set of markers
plink --nonfounders \
 --bfile ${VCFROOT}_chr${chr}.PASS.decomposed.qced \
 --exclude range ~/bin/high-ld_regions_to_exclude_hg38.txt \
 --geno 0.01 \
 --hwe 0.0001 \
 --maf 0.1 \
 --indep-pairwise 100 25 0.2 \
 --out ${VCFROOT}_chr${chr}.PASS.decomposed.qced_common_pruned
plink --nonfounders \
 --bfile ${VCFROOT}_chr${chr}.PASS.decomposed.qced \
 --extract ${VCFROOT}_chr${chr}.PASS.decomposed.qced_common_pruned.prune.in \
 --make-bed \
 --out ${VCFROOT}_chr${chr}.PASS.decomposed_common_pruned_indep
