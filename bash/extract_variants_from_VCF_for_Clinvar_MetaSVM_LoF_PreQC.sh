#!/usr/bin/bash
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/attr_fraction/annotated_variants/annotated_vars_from_PreQC_VCF_Clinvar_MetaSVM_LoF

module load bcftools/1.9
module load plink/1.90b

grep -w chr${CHR} Variants_from_annotation_Clinvar_MetaSVM_LoF_PreQC.bed | sort -V > Annotated_Pathogenic_vars_chr${CHR}.bed
sed -i "s/\r//g" Annotated_Pathogenic_vars_chr${CHR}.bed
bcftools view -Oz MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr${CHR}.preQC_biallelic_renamed_ID_edited.vcf.gz --threads ${THREADS} -R Annotated_Pathogenic_vars_chr${CHR}.bed > Annotated_Pathogenic_PreQC_chr${CHR}.vcf.gz
plink --vcf Annotated_Pathogenic_PreQC_chr${CHR}.vcf.gz --double-id --vcf-half-call m --keep-allele-order --threads ${THREADS} --make-bed --out Annotated_Pathogenic_PreQC_chr${CHR}
