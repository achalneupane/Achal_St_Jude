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

## Run the R code dox_analysis_missense_only_v2.R
module load plink/1.90b
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr12.preQC_biallelic_renamed_ID_edited.vcf.gz --keep-allele-order --extract extract_vars_dox.txt --make-bed --out extract_vars_dox

## recode
plink --bfile /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr12.preQC_biallelic_renamed_ID_edited.vcf.gz --keep-allele-order --extract extract_vars_dox.txt --make-bed --out extract_vars_dox
plink --bfile extract_vars_dox --keep-allele-order --recodeA --out extract_vars_dox_recodeA
plink --bfile extract_vars_dox -freq --out extract_vars_dox_freq
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6}'  extract_vars_dox_freq.frq >  extract_vars_dox_fre_tabsep.frq


## round 2: Email 09/03/2024
# Hi Achal,
# Similar to what you did before, can you please look at these additional genes? Please look at the genes together if they are on the same chromosome but separately for those in different chromosomes. Also, look at all of them together.
# This is related to a grant submission next month so can you please do as soon as you finish the analyses/table etc. for TTN paper?
# Thanks,
# Yadav
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/annotation/snpEff
mkdir yadav_dox_abstract_09_04_2024

head -1 MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr12.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-latest-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.gnomAD-FIELDS-simple.txt > yadav_dox_abstract_09_04_2024/header
grep -wEi 'Hspa12a|Rbm20|Adrb1|Pdzd8|Pdcd4|Shoc2|Afap1l2|Ablim1|Slc18a2|Add3|Dusp5|Bbip1|Adra2a|Acsl5|Vti1a|Tcf7l2|Nrap|Casp7|Nhlrc2|Gfra1|Pnliprp2|Rab11fip2|Mycbp2|Abcc4|Slain1|Ednrb|Gpc6|Tgds|Gpr180|Uggt2|Mbnl2' MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr*.preQC_biallelic_renamed_ID_edited.vcf-annot-snpeff-latest-dbnsfp-ExAC.0.3-clinvar.GRCh38.vcf.dbSNP155.gnomAD-FIELDS-simple.txt > yadav_dox_abstract_09_04_2024/dox_abstract_tmp.txt
cd yadav_dox_abstract_09_04_2024
cat header dox_abstract_tmp.txt  > dox_abstract.txt

## Run the R code dox_analysis_missense_only_v2_round2.R
