## Extracting variants on completely QCed (including GQ, DP and VQSR) data
cd /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/WES_rare_variant/sjlife
## This is the QC'ed data for Survivors
ln -s /research_jude/rgs01_jude/groups/sapkogrp/projects/Genomics/common/Survivor_WES_QC/biallelic2/plink_all/chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated.* .

# Get rare_variants_to_extract.txt from extract_PLP.R
plink --bfile chrALL.SJLIFE_CCSS_WES_101724.GATK4180.hg38_biallelic.geno.0.1.hwe.1e-15.LCR.removed.MAC.ge.1_updated --extract rare_variants_to_extract.txt --keep sjlife_samples.txt --keep-allele-order --make-bed -out all_rare_variants_all_sjlife
plink --bfile all_rare_variants_all_sjlife --extract rare_variants_to_extract.txt --max-maf 0.01 --keep-allele-order --make-bed -out all_rare_variants_maf0.01_all_sjlife

## recode to A
plink --bfile all_rare_variants_all_sjlife --recodeA --out all_rare_variants_all_sjlife_recodeA
plink --bfile all_rare_variants_maf0.01_all_sjlife --recodeA --out all_rare_variants_maf0.01_all_sjlife_recodeA
## recode to A maf 0.0001
plink --bfile all_rare_variants_all_sjlife --extract rare_variants_to_extract.txt --max-maf 0.0001 --keep-allele-order --make-bed -out all_rare_variants_maf0.0001_all_sjlife
plink --bfile all_rare_variants_maf0.0001_all_sjlife --recodeA --out all_rare_variants_maf0.0001_all_sjlife_recodeA