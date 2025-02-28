# 2/24/2025: Hi Achal,
# Can you please look into this variant with respect to cardiovascular outcome in AA survivors in SJLIFE? Please try looking at any CV late effect and then by individual late effect as outcomes.

# variant: rs76992529 (chr18:31598655:G:A)
#1. https://pmc.ncbi.nlm.nih.gov/articles/PMC4874558/
#2. https://www.nature.com/articles/s41598-021-91113-6
 
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/MERGED_sjlife1_2_PreQC/cleaned/nihms

ln -s ../MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr18.preQC_biallelic_renamed_ID_edited.vcf.gz.*

module load plink/1.90b
plink --bfile MERGED.SJLIFE.1.2.GATKv3.4.VQSR.chr18.preQC_biallelic_renamed_ID_edited.vcf.gz --extract extract_variants.txt --keep-allele-order --make-bed --out subbed
plink --bfile subbed --recodeA --keep-allele-order --out subbed_recodeA 