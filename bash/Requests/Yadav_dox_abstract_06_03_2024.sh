# Can you please do the following when you get a chance?

# Extract all the non-synonymous variants within the four genes (Gnptab, Slc25a3, Uhrf1bp1l, and Chpt1) among all SJLIFE survivors. They are all on chromosome 12. Please use the WGS data.
# Perform association analysis with cardiomyopathy using Fisher’s Exact test (without any covariates). This will be just 2 by 3 table. Please do this for Grade 2 or higher vs. Grade 0, as well as Grade 3 or higher vs. Grade 0. Survivors should be exposed to anthracyclines alone without chest RT.
 
# Can you then repeat above 1) and 2) for pathogenic/likely pathogenic variants only?


cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff
head -1 MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr12.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-latest-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.gnomAD-FIELDS-simple.txt > yadav_dox_abstract_06_03_2024/header
grep -wEi 'Gnptab|Slc25a3|Uhrf1bp1l|Chpt1' MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr12.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-latest-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.gnomAD-FIELDS-simple.txt > yadav_dox_abstract_06_03_2024/dox_abstract_tmp.txt

cd yadav_dox_abstract_06_03_2024
cat header dox_abstract_tmp.txt  > dox_abstract.txt 

awk '{split($3, a, ";"); print a[1]}' dox_abstract.txt > extract_vars_dox.txt

## Run the R code dox_analysis.R
module load plink/1.90b
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr12.preQC_biallelic_renamed_ID_edited.vcf.gz --keep-allele-order --extract extract_vars_dox.txt --make-bed --out extract_vars_dox

## recode
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr12.preQC_biallelic_renamed_ID_edited.vcf.gz --keep-allele-order --extract extract_vars_dox.txt --make-bed --out extract_vars_dox
plink --bfile extract_vars_dox --keep-allele-order --recodeA --out extract_vars_dox_recodeA
plink --bfile extract_vars_dox -freq --out extract_vars_dox_freq
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6}'  extract_vars_dox_freq.frq >  extract_vars_dox_fre_tabsep.frq