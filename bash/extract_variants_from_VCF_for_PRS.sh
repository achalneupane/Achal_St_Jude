#!/usr/bin/bash
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/prs/plink_data

module load bcftools/1.9
module load plink/1.90b

grep -w chr${CHR} PRS_vars.bed | sort -V > PRS_vars_chr${CHR}.bed
bcftools view -Oz MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic_renamed_ID_edited.vcf.gz --threads ${THREADS} -R PRS_vars_chr${CHR}.bed > PRS_chr${CHR}.vcf.gz
plink --vcf PRS_chr${CHR}.vcf.gz --double-id --vcf-half-call m --keep-allele-order --threads ${THREADS} --make-bed --out PRS_chr${CHR}
